"""
Shared Monte Carlo parameter-jitter engine for Figure 5 (malate/rTCA analysis) and, as
of the Figure 4 jitter migration, Figure 4 (dilution/flux-by-module analysis) too.

Jitters the 427 fitted kinetic parameters (Kcat_V_*, Keq_*, Km_*, Ki_*, Gnc_*, Gc_*,
kdeg_*) of trainparams[60] independently, log-normal, sigma_log10 calibrated from the
CatPred confidence figure ("90% of predictions in the highest-confidence subset are
within an order of magnitude error" -> 1.645*sigma = 1 -> sigma_log10 = 1/1.645).
Does NOT jitter the 71 v-prefixed exogenous-target parameters (experimenter-controlled,
not measurement-uncertain) or rel0-3 (already part of trainparams, treated as fixed
scaling factors here, same as the cluster4-ensemble analysis treats them).

Built as a real importable module (not inline notebook/exec code) specifically so that
multiprocessing on Windows (spawn start method) can pickle the worker functions -- each
worker process imports this module directly and gets its own RoadRunner model instance
via _get_rfit(), so no RoadRunner/SWIG object ever crosses a process boundary.

Includes a physical-plausibility filter baked into every panel worker: any single-species
concentration readout (SMalate, CO2) outside [0, PLAUSIBLE_CONC_MAX] mM, or any flux-
integral readout (R20/R25 amounts) with |value| > PLAUSIBLE_FLUX_MAX mM, is treated as a
failed draw (NaN) alongside genuine CVODE RuntimeErrors -- these arise from a small
fraction (~2-3%) of aggressive jitter draws producing numerically "successful" but
biologically nonsensical results (observed magnitudes up to ~1e14) that don't raise
exceptions on their own.
"""
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

MODEL_PATH = 'model/iMC057.txt'
PARAMS_PATH = 'parameters/final_params.pkl'
BASE_INDEX = 60

SIGMA_LOG10 = 1.0 / 1.645  # ~0.6079, from "90% within 1 order of magnitude" -> 1.645*sigma=1
KINETIC_PREFIXES = ('Kcat_V_', 'Keq_', 'Km_', 'Ki_', 'Gnc_', 'Gc_', 'kdeg_')
N_JITTER_DEFAULT = 30

PLAUSIBLE_CONC_MAX = 100.0   # mM, generous vs. all substrate pools being low single-digit mM
PLAUSIBLE_FLUX_MAX = 100.0   # mM, magnitude bound on R20/R25 time-integrated amounts

D_BASE = 1 / 200

# ---------------------------------------------------------------------------
# Model loading (once per process)
# ---------------------------------------------------------------------------
def loada(path):
    antimony.clearPreviousLoads()
    code = antimony.loadAntimonyFile(path)
    if code < 0:
        raise RuntimeError(antimony.getLastError())
    module_name = antimony.getMainModuleName()
    sbml = antimony.getSBMLString(module_name)
    return roadrunner.RoadRunner(sbml)


_RFIT = None
_BASE_PARAMS = None
_KINETIC_KEYS = None
_REACTION_IDS = None


def _get_rfit():
    """Lazily create this process's own RoadRunner instance. Safe across
    multiprocessing workers since each worker has a separate Python interpreter."""
    global _RFIT, _BASE_PARAMS, _KINETIC_KEYS, _REACTION_IDS
    if _RFIT is None:
        _RFIT = loada(MODEL_PATH)
        with open(PARAMS_PATH, 'rb') as f:
            trainparams = pickle.load(f)
        _BASE_PARAMS = dict(trainparams[BASE_INDEX][2])
        _KINETIC_KEYS = [k for k in _BASE_PARAMS if any(k.startswith(p) for p in KINETIC_PREFIXES)]
        _REACTION_IDS = _RFIT.model.getReactionIds()
    return _RFIT


def get_reaction_ids():
    _get_rfit()
    return _REACTION_IDS


def get_param_dicts_from_params(rfit, params):
    """Figure 4's notebook-native equivalent of get_param_dicts(rfit, entry), but takes
    a raw {param_name: value} dict (e.g. from jittered_params(seed)) instead of a
    (id, loss, params_dict) tuple from trainparams -- same index-based
    setGlobalParameterValues([*labels.values(), ...], [*values.values(), ...]) pattern
    used throughout Figure_4.ipynb, unchanged."""
    labels, values = {}, {}
    for k, pid in enumerate(rfit.model.getGlobalParameterIds()):
        if pid in params:
            values[pid] = params[pid]
            labels[pid] = k
    return labels, values


def get_base_params():
    _get_rfit()
    return _BASE_PARAMS


def get_kinetic_keys():
    _get_rfit()
    return _KINETIC_KEYS


def jittered_params(seed):
    """Draw one jittered parameter dict from the base (index-60) parameter set.
    `seed` fully determines the draw so results are reproducible and each
    multiprocessing worker can independently generate draw i without shared state."""
    base = get_base_params()
    kinetic_keys = get_kinetic_keys()
    rng = np.random.default_rng(seed)
    params = dict(base)
    for k in kinetic_keys:
        z = rng.standard_normal()
        params[k] = params[k] * (10 ** (z * SIGMA_LOG10))
    return params


def _apply_params(rfit, params):
    for name, val in params.items():
        rfit.model[name] = val


def _check_conc(x):
    if x is None or not np.isfinite(x) or x < 0 or x > PLAUSIBLE_CONC_MAX:
        return np.nan
    return x


def _check_flux(x):
    if x is None or not np.isfinite(x) or abs(x) > PLAUSIBLE_FLUX_MAX:
        return np.nan
    return x


# ---------------------------------------------------------------------------
# Shared background-setup pattern (native, zero-heterologous-enzyme baseline)
# ---------------------------------------------------------------------------
def _native_background(rfit, d, params):
    _apply_params(rfit, params)
    rfit.model['dilution_factor'] = d
    rfit.model['p_hEC6411'] = 0
    rfit.model['p_hEC11719'] = 0
    rfit.model['p_hEC43117'] = 0
    rfit.reset()
    rel = (rfit.model['rel0'], rfit.model['rel1'], rfit.model['rel2'], rfit.model['rel3'])
    pre = {
        'NADH': rfit.getValue('NADH'), 'NADPH': rfit.getValue('NADPH'),
        'Formate': rfit.getValue('Formate'), 'ATP': rfit.getValue('ATP'),
        'AcetylCoA': rfit.getValue('AcetylCoA'), 'HCO3': rfit.getValue('HCO3'),
        'eEC11137': rfit.getValue('eEC11137'), 'Pyruvate': rfit.getValue('Pyruvate'),
        'Oxaloacetate': rfit.getValue('Oxaloacetate'), 'L_Serine': rfit.getValue('L_Serine'),
    }
    return rel, pre


def _set_common_cofactors(rfit, d, pre, substrate=None, substrate_conc=1.0):
    rfit.model['v73'] = 0
    rfit.model['v44'] = (1 / 0.083 / d) + (pre['NADH'] / 0.083 / d)
    rfit.model['v46'] = (1 / 0.12 / d) + (pre['NADPH'] / 0.12 / d)
    rfit.model['v74'] = (10 / 0.001 / d) + (pre['Formate'] / 0.001 / d)
    rfit.model['v6'] = (1 / 1.5 / d) + (pre['ATP'] / 1.5 / d)
    rfit.model['v10'] = (1 / 0.61 / d) + (pre['AcetylCoA'] / 0.61 / d)
    rfit.model['v35'] = (10 / 0.015 / d) + (pre['HCO3'] / 0.015 / d)
    rfit.model['v53'] = (substrate_conc / 0.001 / d) + (pre['Pyruvate'] / 0.001 / d) if substrate == 'Pyruvate' else 0
    rfit.model['v50'] = (substrate_conc / 0.001 / d) + (pre['Oxaloacetate'] / 0.001 / d) if substrate == 'Oxaloacetate' else 0
    rfit.model['v38'] = (substrate_conc / 0.001 / d) + (pre['L_Serine'] / 0.001 / d) if substrate == 'Serine' else 0


def _set_enzymes(rfit, d, rel, pre, hEC6411, hEC11719, hEC43117, eEC11137_incr):
    rel0, rel1, rel2, rel3 = rel
    rfit.model['dilution_factor'] = d
    rfit.model['p_hEC6411'] = hEC6411 / (0.01 * d * rel1) if hEC6411 else 0
    rfit.model['p_hEC11719'] = hEC11719 / (0.01 * d * rel2) if hEC11719 else 0
    rfit.model['p_hEC43117'] = hEC43117 / (0.01 * d * rel3) if hEC43117 else 0
    if eEC11137_incr:
        rfit.model['p_eEC11137'] = (eEC11137_incr / (0.012601796 * d * rel0)) + (pre['eEC11137'] / (0.012601796 * d * rel0))


# ---------------------------------------------------------------------------
# Figure 5B: 8 named conditions
# ---------------------------------------------------------------------------
CONDITIONS_5B = [
    ('- Pathway + Oxaloacetate', (0, 0, 0, 0), 'Oxaloacetate'),
    ('mdh', (0, 0, 0, 0.025), 'Oxaloacetate'),
    ('fdh', (0, 0.025, 0, 0), 'Oxaloacetate'),
    ('mdh+fdh', (0, 0.025, 0, 0.025), 'Oxaloacetate'),
    ('- Pathway + Pyruvate', (0, 0, 0, 0), 'Pyruvate'),
    ('mdh+fdh+pyc', (0.050, 0.025, 0, 0.025), 'Pyruvate'),
    ('- Pathway + Serine', (0, 0, 0, 0), 'Serine'),
    ('mdh+fdh+pyc+sds', (0.050, 0.025, 0.050, 0.025), 'Serine'),
]


def run_5b_draw(params, d=D_BASE):
    rfit = _get_rfit()
    out = {}
    for label, enzymes, substrate in CONDITIONS_5B:
        rel, pre = _native_background(rfit, d, params)
        _set_common_cofactors(rfit, d, pre, substrate=substrate, substrate_conc=1.0)
        _set_enzymes(rfit, d, rel, pre, *enzymes)
        rfit.reset()
        try:
            results = rfit.simulate(0, 4 * 60 * 60, selections=['time', 'SMalate'])
            val = _check_conc(results[-1][1])
        except RuntimeError:
            val = np.nan
        rfit.resetToOrigin()
        out[label] = val
    return out


# ---------------------------------------------------------------------------
# Figure 5C: enzyme titrations (4 enzymes x 3 concs x 3 substrates = 36 conditions)
# ---------------------------------------------------------------------------
ENZYMES_5C = ['fdh', 'mdh', 'pyc', 'sds']
CONCS_5C = [0, 0.001, 0.005]
SUBSTRATES_5C = ['Oxaloacetate', 'Pyruvate', 'Serine']
BASE_DOSE_5C = {'fdh': 0.025, 'mdh': 0.025, 'pyc': 0.050, 'sds': 0.050}


def _run_5c_condition(rfit, d, rel, pre, enzyme, conc, substrate):
    rel0, rel1, rel2, rel3 = rel
    dose = {k: (conc if k == enzyme else BASE_DOSE_5C[k]) for k in ENZYMES_5C}
    _set_common_cofactors(rfit, d, pre, substrate=substrate, substrate_conc=1.0)
    _set_enzymes(rfit, d, rel, pre, dose['pyc'], dose['fdh'], dose['sds'], dose['mdh'])
    rfit.reset()
    try:
        results = rfit.simulate(0, 4 * 60 * 60, selections=['time', 'SMalate'])
        val = _check_conc(results[-1][1])
    except RuntimeError:
        val = np.nan
    rfit.resetToOrigin()
    return val


def run_5c_draw(params, d=D_BASE):
    rfit = _get_rfit()
    out = {}
    for substrate in SUBSTRATES_5C:
        for enzyme in ENZYMES_5C:
            for conc in CONCS_5C:
                rel, pre = _native_background(rfit, d, params)
                val = _run_5c_condition(rfit, d, rel, pre, enzyme, conc, substrate)
                out[f"{substrate}|{enzyme}|{conc}"] = val
    return out


# ---------------------------------------------------------------------------
# Figure 5D: substrate titration (3 substrates x 10 concs = 30 conditions), full
# base-case enzyme dose held constant throughout.
# ---------------------------------------------------------------------------
TITRATION_CONCS_5D = [0, 0.01, 0.5, 1, 5, 10, 20, 50, 100, 1000]
SUBSTRATES_5D = ['Serine', 'Pyruvate', 'Oxaloacetate']


def _run_5d_condition(rfit, d, rel, pre, substrate, conc):
    _set_common_cofactors(rfit, d, pre, substrate=substrate, substrate_conc=conc)
    _set_enzymes(rfit, d, rel, pre, 0.050, 0.025, 0.050, 0.025)
    rfit.reset()
    try:
        results = rfit.simulate(0, 4 * 60 * 60, selections=['time', 'SMalate'])
        val = _check_conc(results[-1][1])
    except RuntimeError:
        val = np.nan
    rfit.resetToOrigin()
    return val


def run_5d_draw(params, d=D_BASE):
    rfit = _get_rfit()
    out = {}
    for substrate in SUBSTRATES_5D:
        for conc in TITRATION_CONCS_5D:
            rel, pre = _native_background(rfit, d, params)
            val = _run_5d_condition(rfit, d, rel, pre, substrate, conc)
            out[f"{substrate}|{conc}"] = val
    return out


# ---------------------------------------------------------------------------
# Figure 5E: dilution titration (Serine-fed full-pathway recipe), reading final
# [Malate], final [CO2], initial [CO2].
# ---------------------------------------------------------------------------
DILFACS_5E = [200, 100, 50, 20, 10]


def run_5e_draw(params):
    rfit = _get_rfit()
    out = {}
    for dilfac in DILFACS_5E:
        d_local = 1 / dilfac
        rel, pre = _native_background(rfit, d_local, params)
        _set_common_cofactors(rfit, d_local, pre, substrate='Serine', substrate_conc=1.0)
        _set_enzymes(rfit, d_local, rel, pre, 0.050, 0.025, 0.050, 0.025)
        rfit.reset()
        try:
            results = rfit.simulate(0, 4 * 60 * 60, selections=['time', 'SMalate', 'CO2'])
            fm = _check_conc(results[-1][1])
            fc = _check_conc(results[-1][2])
            ic = _check_conc(results[0][2])
        except RuntimeError:
            fm, fc, ic = np.nan, np.nan, np.nan
        rfit.resetToOrigin()
        out[dilfac] = (fm, fc, ic)
    return out


# ---------------------------------------------------------------------------
# Figure 5F: R20 (TCA)/R25 (rTCA) time-integrated malate amounts across dilution,
# FAST dense-simulate method (single simulate() call per draw) -- used ONLY for the
# jitter sub-analysis. The main (non-jittered, cluster4) Figure 5F keeps the original
# chunked-loop method untouched, per earlier explicit user request.
# ---------------------------------------------------------------------------
DILFACS_5F = [10, 20, 50, 100, 200]
N_TIME_POINTS_5F = 961  # matches original np.arange(0, 4*60*60+1, 15) resolution
TIME_POINTS_5F = np.arange(0, 4 * 60 * 60 + 1, 15)
REACTIONS_OF_INTEREST_5F = [19, 24]  # R20, R25

# (hEC6411/pyc, hEC11719/fdh, hEC43117/sds, eEC11137 increment/mdh) in mM
ENZYMES_5F_SERINE = (0.050, 0.025, 0.050, 0.025)     # full pathway, matches Figure 5B's "mdh+fdh+pyc+sds"
ENZYMES_5F_PYRUVATE = (0.050, 0.025, 0, 0.025)       # SDS not needed w/o Serine, matches "mdh+fdh+pyc"

# Heterologous enzymes (and the MDH overexpression increment) are expressed once from
# plasmids during the UNDILUTED cell-free gene-expression (TXTL) phase -- the 150 uM
# budget, and each enzyme's target dose, is defined AT D=10, not at whatever dilution
# factor the biosynthesis phase later uses. When the lysate (now containing that newly
# expressed protein) is subsequently diluted with buffer/substrates/cofactors for
# biosynthesis, that protein dilutes right along with everything else in the lysate --
# it should NOT be held at a constant concentration regardless of biosynthesis D.
D_UNDILUTED = 10
D_UNDILUTED_FACTOR = 1 / D_UNDILUTED  # 0.1


def _set_enzymes_expressed_undiluted(rfit, d, rel, pre, hEC6411, hEC11719, hEC43117, eEC11137_incr):
    """Correct, dilution-aware version of _set_enzymes for Figure 5F: solves for
    p_hEC*/p_eEC11137-increment ONCE at the fixed undiluted reference (D=10) rather than
    at the currently-simulated `d`, so the resulting protein concentration naturally
    scales down as (D_UNDILUTED/D) for any more-dilute biosynthesis condition -- matching
    the real two-phase protocol (undiluted TXTL expression, then dilution for
    biosynthesis) instead of artificially holding concentration constant across D.
    The native MDH background component (`pre['eEC11137']`, already measured at the
    current simulated `d`) is unaffected -- only the plasmid-driven increment changes."""
    rel0, rel1, rel2, rel3 = rel
    rfit.model['dilution_factor'] = d
    rfit.model['p_hEC6411'] = hEC6411 / (0.01 * D_UNDILUTED_FACTOR * rel1) if hEC6411 else 0
    rfit.model['p_hEC11719'] = hEC11719 / (0.01 * D_UNDILUTED_FACTOR * rel2) if hEC11719 else 0
    rfit.model['p_hEC43117'] = hEC43117 / (0.01 * D_UNDILUTED_FACTOR * rel3) if hEC43117 else 0
    if eEC11137_incr:
        rfit.model['p_eEC11137'] = (eEC11137_incr / (0.012601796 * D_UNDILUTED_FACTOR * rel0)) + (pre['eEC11137'] / (0.012601796 * d * rel0))


def run_5f_draw_fast_variant(params, substrate, enzymes):
    """Generalized fast dense-simulate R20/R25 method (one simulate() call per
    dilution factor), parameterized by substrate and enzyme doses."""
    rfit = _get_rfit()
    out = {}
    for dilfac in DILFACS_5F:
        d_local = 1 / dilfac
        rel, pre = _native_background(rfit, d_local, params)
        _set_common_cofactors(rfit, d_local, pre, substrate=substrate, substrate_conc=1.0)
        _set_enzymes_expressed_undiluted(rfit, d_local, rel, pre, *enzymes)
        rfit.reset()
        try:
            results = rfit.simulate(0, 4 * 60 * 60, N_TIME_POINTS_5F,
                                     selections=['time', 'R20', 'R25'])
            t = np.asarray(results['time'])
            r20_flux = -np.asarray(results['R20'])
            r25_flux = -np.asarray(results['R25'])
            r20 = _check_flux(np.trapz(r20_flux, t))
            r25 = _check_flux(np.trapz(r25_flux, t))
        except RuntimeError:
            r20, r25 = np.nan, np.nan
        rfit.resetToOrigin()
        out[dilfac] = (r20, r25)
    return out


def run_5f_draw_fast(params):
    return run_5f_draw_fast_variant(params, substrate='Serine', enzymes=ENZYMES_5F_SERINE)


def run_5f_draw_fast_pyruvate(params):
    return run_5f_draw_fast_variant(params, substrate='Pyruvate', enzymes=ENZYMES_5F_PYRUVATE)


def run_5f_reference_chunked(params, dilfac, substrate, enzymes):
    """Original chunked-simulate-loop + raw-trapz R20/R25 method (matches the notebook's
    run_5f_member exactly), single point (no jitter), parameterized by substrate/enzymes.
    Used only for the cheap single-point reference, not the (expensive) jittered ensemble."""
    rfit = _get_rfit()
    d_local = 1 / dilfac
    rel, pre = _native_background(rfit, d_local, params)
    _set_common_cofactors(rfit, d_local, pre, substrate=substrate, substrate_conc=1.0)
    _set_enzymes_expressed_undiluted(rfit, d_local, rel, pre, *enzymes)
    rfit.reset()

    fluxes = {rid: [] for rid in REACTIONS_OF_INTEREST_5F}
    for j in range(len(TIME_POINTS_5F) - 1):
        try:
            rfit.simulate(int(TIME_POINTS_5F[j]), int(TIME_POINTS_5F[j + 1]), 2)
            reaction_fluxes = rfit.getReactionRates()
            for rid in REACTIONS_OF_INTEREST_5F:
                fluxes[rid].append(-(reaction_fluxes[rid]))
        except RuntimeError:
            rfit.reset()
            rfit.simulate(0, int(TIME_POINTS_5F[j]))
            rfit.simulate(start=int(TIME_POINTS_5F[j]), end=int(TIME_POINTS_5F[j + 1]), steps=2)
            reaction_fluxes = rfit.getReactionRates()
            for rid in REACTIONS_OF_INTEREST_5F:
                fluxes[rid].append(-(reaction_fluxes[rid]))

    product_amounts = {}
    for rid in REACTIONS_OF_INTEREST_5F:
        v = np.array(fluxes[rid])
        product_amounts[rfit.model.getReactionIds()[rid]] = np.trapz(v, TIME_POINTS_5F[1:])
    rfit.resetToOrigin()
    return product_amounts['R20'], product_amounts['R25']


# ---------------------------------------------------------------------------
# Figure 5G: MCA-inspired control coefficients, restricted to a caller-supplied
# species list (the ~18-20 species that already survive the |C|>1e-3 filter in the
# existing cluster4 point-estimate run -- passed in explicitly, not all 71).
# ---------------------------------------------------------------------------
def _load_species_lookup():
    import re
    with open(MODEL_PATH, 'r') as f:
        text = f.read()
    pattern = re.compile(r"(?P<species>\S+)\s*=\s*(?P<v>v\d+)\*(?P<val>[0-9.eE+-]+)\*dilution_factor;")
    lookup = {}
    for match in pattern.finditer(text):
        try:
            val = float(match.group("val"))
        except ValueError:
            val = None
        lookup[match.group("species")] = {"v": match.group("v"), "val": val}
    return lookup


_LOOKUP = None


def get_species_lookup():
    global _LOOKUP
    if _LOOKUP is None:
        _LOOKUP = _load_species_lookup()
    return _LOOKUP


def _run_5g_perturbation(rfit, d, params, s, direction, lookup):
    rel, pre = _native_background(rfit, d, params)
    pre_s = rfit.getValue(s)
    _set_common_cofactors(rfit, d, pre, substrate='Serine', substrate_conc=1.0)

    v_name = lookup[s]["v"]
    ref_val = lookup[s]["val"]
    multiplier = 100 if direction > 0 else (1 / 100)
    rfit.model[v_name] = pre_s * multiplier / ref_val / d

    _set_enzymes(rfit, d, rel, pre, 0.050, 0.025, 0.050, 0.025)
    rfit.reset()
    try:
        results = rfit.simulate(0, 4 * 60 * 60, selections=['time', 'SMalate'])
        val = _check_conc(results[-1][1])
    except RuntimeError:
        val = np.nan
    rfit.resetToOrigin()
    return val


def run_5g_draw(params, species_list, d=D_BASE):
    rfit = _get_rfit()
    lookup = get_species_lookup()
    out = {}
    for s in species_list:
        malate_up = _run_5g_perturbation(rfit, d, params, s, +1, lookup)
        malate_down = _run_5g_perturbation(rfit, d, params, s, -1, lookup)
        if (np.isnan(malate_up) or np.isnan(malate_down) or malate_up <= 0 or malate_down <= 0):
            out[s] = np.nan
        else:
            out[s] = (np.log(malate_up) - np.log(malate_down)) / np.log(1e4)
    return out


# ---------------------------------------------------------------------------
# Figure 4: dilution/flux-by-module analysis. Ports Figure_4.ipynb's per-cluster4-member
# simulation logic to per-jittered-draw, verbatim except for swapping the ensemble
# source (26 independently-fit cluster4 parameter sets -> N jittered draws around
# trainparams[60]) -- same index-based setGlobalParameterValues([2,49,51,53]) native
# "no heterologous enzyme" background pattern used throughout the original notebook.
# ---------------------------------------------------------------------------
PATHWAY_DICT = {
    'TCA Cycle': [
        'R1', 'R7', 'R8', 'R20', 'R25', 'R27', 'R28', 'R41', 'R42', 'R46', 'R47'
    ],
    'Glycolysis / Gluconeogenesis': [
        'R10', 'R13', 'R18', 'R19', 'R21', 'R22', 'R23', 'R30', 'R32', 'R33',
        'R34', 'R35', 'R37', 'R38', 'R45'
    ],
    'Pentose Phosphate Pathway': [
        'R39', 'R40', 'R43', 'R44', 'R54',
    ],
    'Fatty Acid Biosynthesis': [
        'R2', 'R3', 'R4', 'R6', 'R11', 'R12', 'R14', 'R24', 'R26'
    ],
    'Fatty Acid Degradation (β-Oxidation)': [
        'R5', 'R15', 'R16', 'R17'
    ],
    'Pyruvate Metabolism': [
        'R31', 'R36', 'R48', 'R50', 'R53'
    ],
    'Glycine & Serine Metabolism': [
        'R49', 'R52'
    ],
    'Energy & Nucleotide Metabolism + CO2/Bicarb Conversion': [
        'R9', 'R29', 'R51'
    ]
}

# ---- Figure 4A: module flux fraction + total flux over time, across 4 dilutions ----
D_FACTORS_4A = [1 / 200, 1 / 100, 1 / 50, 1 / 10]
TIME_POINTS_4A = np.arange(0, 4 * 60 * 60 + 1, 60 * 60)  # 4 one-hour intervals


def run_4a_draw(params):
    """Note: unlike the original cluster4 loop (26 pre-validated good fits, never
    needed a try/except here), jittered draws span a much wider parameter space and a
    fraction fail to integrate (CVODE errors) -- caught per-dilution-factor and filled
    with NaN, matching the failure-handling convention used throughout Figure 5's
    draw functions (run_5b_draw et al.)."""
    rfit = _get_rfit()
    reaction_ids = get_reaction_ids()
    labels, values = get_param_dicts_from_params(rfit, params)
    n_steps = len(TIME_POINTS_4A) - 1
    out = {}
    for d in D_FACTORS_4A:
        try:
            rfit.model.setGlobalParameterValues(
                [*labels.values(), 2, 49, 51, 53],
                [*values.values(), d, 0, 0, 0]
            )
            rfit.reset()

            module_fluxes = {m: [] for m in PATHWAY_DICT}
            total_flux_list = []
            for j in range(len(TIME_POINTS_4A) - 1):
                rfit.simulate(int(TIME_POINTS_4A[j]), int(TIME_POINTS_4A[j + 1]))
                reaction_fluxes = np.abs(rfit.getReactionRates())
                total_flux = np.sum(reaction_fluxes)
                total_flux_list.append(float(total_flux))

                step = {}
                for module, labs in PATHWAY_DICT.items():
                    size = len(labs) if len(labs) > 0 else 1
                    flux = sum(reaction_fluxes[reaction_ids.index(l)] for l in labs if l in reaction_ids)
                    step[module] = flux / size
                total_at_t = sum(step.values()) if total_flux > 0 else 1
                for module in step:
                    module_fluxes[module].append(float(step[module] / total_at_t))

            out[str(d)] = {
                "frac": {m: module_fluxes[m] for m in PATHWAY_DICT},
                "total_flux": total_flux_list,
            }
        except RuntimeError:
            out[str(d)] = {
                "frac": {m: [float('nan')] * n_steps for m in PATHWAY_DICT},
                "total_flux": [float('nan')] * n_steps,
            }
        rfit.resetToOrigin()
    return out


# ---- Figure 4B: total integrated flux per reaction, across 7 dilutions ----
D_FACTORS_4B = np.array([1 / 200, 1 / 100, 1 / 50, 1 / 30, 1 / 20, 1 / 15, 1 / 10])


def _simulate_and_record_rates_4b(rfit, d, labels, values):
    totaltime = 60
    addedtime = 60 * 5

    rfit.model.setGlobalParameterValues(
        [*labels.values(), 2, 49, 51, 53],
        [*values.values(), d, 0, 0, 0]
    )
    rfit.reset()
    rfit.simulate(0, totaltime)
    ratesarray = np.atleast_2d(rfit.getReactionRates()).T

    while totaltime <= 4 * 60 * 60:
        rfit.model.setGlobalParameterValues(
            [*labels.values(), 2, 49, 51, 53],
            [*values.values(), d, 0, 0, 0]
        )
        rfit.reset()
        rfit.simulate(0, totaltime + addedtime)
        rates = np.atleast_2d(rfit.getReactionRates()).T
        ratesarray = np.concatenate((ratesarray, rates), axis=1)
        totaltime += addedtime
        rfit.resetToOrigin()

    rfit.resetToOrigin()
    return ratesarray


def _calc_total_flux_4b(ratesarray):
    addedtime = 300
    return np.sum(np.abs(ratesarray), axis=1) * addedtime


def run_4b_draw(params):
    rfit = _get_rfit()
    n_reactions = len(get_reaction_ids())
    labels, values = get_param_dicts_from_params(rfit, params)
    member_flux = []
    for d in D_FACTORS_4B:
        try:
            rates_array = _simulate_and_record_rates_4b(rfit, d, labels, values)
            flux = _calc_total_flux_4b(rates_array).tolist()
        except RuntimeError:
            flux = [float('nan')] * n_reactions
            rfit.resetToOrigin()
        member_flux.append(flux)
    return member_flux  # (7, n_reactions)


# ---- Figure 4C: reaction x dilution heatmap, integrated flux over 4 dilutions ----
D_FACTORS_4C = [1 / 10, 1 / 50, 1 / 100, 1 / 200]
TIME_POINTS_4C = np.arange(0, 4 * 60 * 60 + 1, 30 * 60)  # 30-min steps over 4 h
DT_4C = 30 * 60


def run_4c_draw(params):
    rfit = _get_rfit()
    labels, values = get_param_dicts_from_params(rfit, params)
    member_rows = []
    for d in D_FACTORS_4C:
        try:
            vcopy = values.copy()
            vcopy['v73'] = 0  # hydroxycitrate placeholder, see CLAUDE.md note on v73/REF=1e-7
            rfit.model.setGlobalParameterValues(
                [*labels.values(), 2, 49, 51, 53],
                [*vcopy.values(), d, 0, 0, 0]
            )
            rfit.reset()
            rate_records = []
            for j in range(len(TIME_POINTS_4C) - 1):
                rfit.simulate(int(TIME_POINTS_4C[j]), int(TIME_POINTS_4C[j + 1]))
                rate_records.append(rfit.getReactionRates()[:57])
            rate_records = np.array(rate_records)
            integrated_flux = np.trapz(rate_records, dx=DT_4C, axis=0)
            member_rows.append(integrated_flux.tolist())
        except RuntimeError:
            member_rows.append([float('nan')] * 57)
        rfit.resetToOrigin()
    return member_rows  # (4, 57)


# ---- Figure 4D: metabolite-spike sink analysis, across 3 spikes x 2 dilutions ----
SPIKES_4D = {"Malate": "v56", "Pyruvate": "v53", "Serine": "v38"}
DILUTIONS_4D = {"10X": 1 / 10, "200X": 1 / 200}

_SUBSET_NAMES_4D = None
_SUBSET_IDX_4D = None


def _get_subset_info_4d(rfit):
    global _SUBSET_NAMES_4D, _SUBSET_IDX_4D
    if _SUBSET_NAMES_4D is None:
        probe = rfit.simulate(0, 1)
        ids = [c.strip('[]') for c in probe.colnames]
        idx = list(range(51, 81)) + list(range(82, 95)) + list(range(97, 122)) + [124, 125]
        _SUBSET_NAMES_4D = [ids[i] for i in idx]
        _SUBSET_IDX_4D = idx
        rfit.resetToOrigin()
    return _SUBSET_NAMES_4D, _SUBSET_IDX_4D


def get_subset_names_4d():
    rfit = _get_rfit()
    names, _ = _get_subset_info_4d(rfit)
    return names


N_TIME_DEFAULT_4D = 51  # RoadRunner's default simulate(start, end) output resolution


def run_4d_draw(params):
    rfit = _get_rfit()
    labels, values_base = get_param_dicts_from_params(rfit, params)
    subset_names, subset_idx = _get_subset_info_4d(rfit)
    n_subset = len(subset_idx)

    out = {}
    for spike_name, spike_param in SPIKES_4D.items():
        out[spike_name] = {}
        for d_name, d in DILUTIONS_4D.items():
            try:
                param_values = dict(values_base)

                # negative control (no spike), all three heterologous enzymes zeroed
                param_values[spike_param] = 0
                rfit.model.setGlobalParameterValues(
                    [*labels.values(), 2, 49, 51, 53],
                    [*param_values.values(), d, 0, 0, 0]
                )
                rfit.reset()
                sim_neg = rfit.simulate(0, 4 * 60 * 60)

                if spike_name == 'Malate':
                    param_values[spike_param] = 5 / 1.7 / d
                else:
                    param_values[spike_param] = 5 / 0.001 / d
                rfit.model.setGlobalParameterValues(
                    [*labels.values(), 2, 49, 51, 53],
                    [*param_values.values(), d, 0, 0, 0]
                )
                rfit.reset()
                sim_spike = rfit.simulate(0, 4 * 60 * 60)

                change = np.asarray(sim_spike) - np.asarray(sim_neg)
                subset_changes = change[:, subset_idx]

                out[spike_name][d_name] = {
                    "tlist": [float(t) for t in sim_spike['time']],
                    "changes": subset_changes.tolist(),
                }
            except RuntimeError:
                out[spike_name][d_name] = {
                    "tlist": [float('nan')] * N_TIME_DEFAULT_4D,
                    "changes": [[float('nan')] * n_subset for _ in range(N_TIME_DEFAULT_4D)],
                }
            rfit.resetToOrigin()
    return out


# ---------------------------------------------------------------------------
# Parallel driver
# ---------------------------------------------------------------------------
def _worker_5b(seed):
    return run_5b_draw(jittered_params(seed))


def _worker_5c(seed):
    return run_5c_draw(jittered_params(seed))


def _worker_5d(seed):
    return run_5d_draw(jittered_params(seed))


def _worker_5e(seed):
    return run_5e_draw(jittered_params(seed))


def _worker_5f(seed):
    return run_5f_draw_fast(jittered_params(seed))


def _worker_5f_pyruvate(seed):
    return run_5f_draw_fast_pyruvate(jittered_params(seed))


def _worker_4a(seed):
    return run_4a_draw(jittered_params(seed))


def _worker_4b(seed):
    return run_4b_draw(jittered_params(seed))


def _worker_4c(seed):
    return run_4c_draw(jittered_params(seed))


def _worker_4d(seed):
    return run_4d_draw(jittered_params(seed))


_WORKER_5G_SPECIES = None


def _init_5g_worker(species_list):
    global _WORKER_5G_SPECIES
    _WORKER_5G_SPECIES = species_list
    _get_rfit()  # pre-load model in this worker


def _worker_5g(seed):
    return run_5g_draw(jittered_params(seed), _WORKER_5G_SPECIES)


def run_parallel(worker_fn, n_draws, n_workers=6, base_seed=1000, initializer=None, initargs=()):
    import multiprocessing as mp
    ctx = mp.get_context('spawn')
    seeds = [base_seed + i for i in range(n_draws)]
    with ctx.Pool(processes=n_workers, initializer=initializer, initargs=initargs) as pool:
        results = pool.map(worker_fn, seeds)
    return results
