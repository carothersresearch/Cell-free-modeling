"""
Monte Carlo kinetic-parameter jitter engine for Figure 6 (malate design campaign +
citrate-synthase/dilution synergy analysis), reusing scripts/sim2.py's validated
Design class and model-building logic. Same jitter methodology as
scripts/jitter_lib.py (Figure 5): N=30 draws around trainparams[60], sigma_log10 =
1/1.645 (calibrated from "90% of predictions within one order of magnitude"),
applied only to the 427 kinetic-prefixed parameters (Kcat_V_*, Keq_*, Km_*, Ki_*,
Gnc_*, Gc_*, kdeg_*) -- NOT the v-prefixed exogenous targets or rel0-3.

Does NOT touch the heterologous-enzyme-dilution convention (still holds heterologous
enzyme concentration constant regardless of biosynthesis D, i.e. the pre-fix
behavior) -- Figure 6 explicitly keeps that convention unchanged per the user, unlike
Figure 5F's jitter_lib.py which has the corrected undiluted-TXTL-expression formula.

A separate, real importable module (not inline notebook/exec code) so Windows
multiprocessing (spawn) can pickle worker functions -- each worker imports this
module directly and builds its own RoadRunner model per simulate call via
_fresh_runner_jittered(), so no RoadRunner/SWIG object crosses a process boundary.
"""
import sys
import os
import pickle
import numpy as np
import antimony
import roadrunner
from roadrunner import Config, Logger

Logger.disableLogging()
Config.setValue(Config.ROADRUNNER_DISABLE_PYTHON_DYNAMIC_PROPERTIES, True)
Config.setValue(Config.LOADSBMLOPTIONS_RECOMPILE, False)
Config.setValue(Config.LLJIT_OPTIMIZATION_LEVEL, 4)
Config.setValue(Config.LLVM_SYMBOL_CACHE, True)
Config.setValue(Config.LOADSBMLOPTIONS_OPTIMIZE_GVN, True)
Config.setValue(Config.LOADSBMLOPTIONS_OPTIMIZE_CFG_SIMPLIFICATION, True)
Config.setValue(Config.LOADSBMLOPTIONS_OPTIMIZE_INSTRUCTION_COMBINING, True)
Config.setValue(Config.LOADSBMLOPTIONS_OPTIMIZE_DEAD_INST_ELIMINATION, True)
Config.setValue(Config.LOADSBMLOPTIONS_OPTIMIZE_DEAD_CODE_ELIMINATION, True)
Config.setValue(Config.LOADSBMLOPTIONS_OPTIMIZE_INSTRUCTION_SIMPLIFIER, True)
Config.setValue(Config.SIMULATEOPTIONS_COPY_RESULT, True)

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import sim2
from sim2 import Design

SIGMA_LOG10 = 1.0 / 1.645
KINETIC_PREFIXES = ('Kcat_V_', 'Keq_', 'Km_', 'Ki_', 'Gnc_', 'Gc_', 'kdeg_')
BASE_INDEX = 60
N_JITTER_DEFAULT = 30

_BASE_PARAMS = None
_KINETIC_KEYS = None


def _load_base():
    global _BASE_PARAMS, _KINETIC_KEYS
    if _BASE_PARAMS is None:
        with open(sim2.PARAMS_PATH, 'rb') as f:
            trainparams = pickle.load(f)
        _BASE_PARAMS = dict(trainparams[BASE_INDEX][2])
        _KINETIC_KEYS = [k for k in _BASE_PARAMS if any(k.startswith(p) for p in KINETIC_PREFIXES)]
    return _BASE_PARAMS, _KINETIC_KEYS


def get_base_params():
    base, _ = _load_base()
    return base


def jittered_params(seed):
    """Draw one jittered parameter dict from the base (index-60) parameter set.
    `seed` fully determines the draw for reproducibility across workers."""
    base, kinetic_keys = _load_base()
    rng = np.random.default_rng(seed)
    params = dict(base)
    for k in kinetic_keys:
        z = rng.standard_normal()
        params[k] = params[k] * (10 ** (z * SIGMA_LOG10))
    return params


_RFIT = None


def _get_cached_rfit():
    """Build the RoadRunner model ONCE per worker process and reuse it across every
    design/draw via reset() -- rebuilding from antimony source (as sim2._fresh_runner
    and an earlier version of this function did) costs ~8s per call (JIT compilation),
    which would dominate wall-clock time across thousands of jittered-draw x design
    combinations. Matches the caching pattern in scripts/jitter_lib.py::_get_rfit."""
    global _RFIT
    if _RFIT is None:
        antimony.clearPreviousLoads()
        code = antimony.loadAntimonyString(sim2._FULL_TEXT)
        if code < 0:
            raise RuntimeError(antimony.getLastError())
        mod = antimony.getMainModuleName()
        sbml = antimony.getSBMLString(mod)
        _RFIT = roadrunner.RoadRunner(sbml)
    return _RFIT


def run_design_jittered(design: Design, params, n_points=400, extra_selections=()):
    """Same construction logic as sim2.run_design, but takes an explicit
    (possibly jittered) params dict instead of the module-level FITTED_PARAMS, and
    accepts extra reaction/species selections (e.g. R46, R46_f, R46_r, AcetylCoA,
    Oxaloacetate) needed by the citrate-synthase synergy panels. Reapplies the full
    jittered kinetic-parameter dict on every call (cheap -- no recompilation) against
    a process-cached model, since every design-relevant knob below is also
    unconditionally overwritten every call, so no state can leak between designs
    without an explicit resetToOrigin()."""
    r = _get_cached_rfit()
    for name, val in params.items():
        r.model[name] = val
    dil = 1.0 / design.dilution_factor_D
    r.model["dilution_factor"] = dil

    for sp, info in sim2.SPECIES_MAP.items():
        if info["driver"].startswith("p_eEC") and sp != "eEC11137":
            mult = design.endog_multipliers.get(sp, 1.0)
            r.model[info["driver"]] = mult

    mdh_info = sim2.SPECIES_MAP["eEC11137"]
    mdh_native_mult = design.endog_multipliers.get("eEC11137", 1.0)
    p_mdh = mdh_native_mult + (design.eEC11137_overexpr_uM / 1000.0) / (mdh_info["ref"] * dil * r.model["rel0"])
    r.model["p_eEC11137"] = p_mdh

    for hid in sim2.HETEROLOGOUS_ENZYMES:
        target_mM = design.het_enzymes_uM.get(hid, 0.0) / 1000.0
        p_val = target_mM / (sim2.HET_REF * dil * r.model[sim2.HET_REL[hid]])
        r.model[f"p_{hid}"] = p_val

    for sp, info in sim2.SPECIES_MAP.items():
        if not info["driver"].startswith("v"):
            continue
        v_default = params.get(info["driver"], None)
        if v_default is None:
            continue
        add_mM = design.exogenous_mM.get(sp, None)
        if add_mM is not None:
            r.model[info["driver"]] = v_default + add_mM / (info["ref"] * dil)
        else:
            r.model[info["driver"]] = v_default

    r.model["v73"] = 10
    r.reset()

    if design.extra_species_overrides_mM:
        for sid, val in design.extra_species_overrides_mM.items():
            r.model[sid] = val

    base_selections = ["time", "SMalate", "Cum_rTCA_Malate", "Cum_TCA_Malate"]
    selections = base_selections + [s for s in extra_selections if s not in base_selections]
    r.selections = selections
    n_pts = max(200, n_points)
    try:
        result = r.simulate(0, design.duration_sec, n_pts)
    except RuntimeError:
        # aggressive jitter draws occasionally produce a genuine CVODE stiffness
        # failure (same failure mode seen throughout Figure 5's jitter work) --
        # treat as a failed draw (NaN) rather than crashing the whole batch.
        out = {
            "final_malate_mM": np.nan, "cum_rTCA_malate_mM": np.nan,
            "cum_TCA_malate_mM": np.nan, "rTCA_to_TCA_ratio": np.nan,
        }
        for s in extra_selections:
            out[s] = np.full(n_pts, np.nan)
        if extra_selections:
            out["time"] = np.linspace(0, design.duration_sec, n_pts)
        return out

    res = np.array(result)
    cols = result.colnames
    idx = {c: i for i, c in enumerate(cols)}

    smalate_final = res[-1, idx["SMalate"]]
    cum_rtca = res[-1, idx["Cum_rTCA_Malate"]]
    cum_tca = res[-1, idx["Cum_TCA_Malate"]]

    # physical-plausibility filter: aggressive jitter draws occasionally "succeed"
    # numerically but integrate to a nonphysical value (seen up to ~1e14 in Figure
    # 5's jitter work) without raising a solver error -- flag as a failed draw.
    if not np.isfinite(smalate_final) or smalate_final < 0 or smalate_final > 100:
        smalate_final = np.nan

    out = {
        "final_malate_mM": float(smalate_final),
        "cum_rTCA_malate_mM": float(cum_rtca),
        "cum_TCA_malate_mM": float(cum_tca),
        "rTCA_to_TCA_ratio": float(cum_rtca / cum_tca) if cum_tca > 1e-12 else float("inf"),
    }
    for s in extra_selections:
        out[s] = res[:, idx[s]]
    if extra_selections:
        out["time"] = res[:, idx["time"]]
    return out


def run_parallel(worker_fn, n_draws, n_workers=6, base_seed=1000, initializer=None, initargs=()):
    import multiprocessing as mp
    ctx = mp.get_context('spawn')
    seeds = [base_seed + i for i in range(n_draws)]
    with ctx.Pool(processes=n_workers, initializer=initializer, initargs=initargs) as pool:
        results = pool.map(worker_fn, seeds)
    return results
