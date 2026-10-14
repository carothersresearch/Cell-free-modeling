"""
Simulation harness for the UPDATED cell-free rTCA malate production model
(model/iMC057.txt + parameters/final_params.pkl), matching the user's own
parameterization pattern given in CLAUDE.md's "Loading parameter set" /
"An example of changing parameters and simulating" sections.

KEY MECHANICS (reverse-engineered directly from the model file + cross
-validated numerically against roadrunner's own computed defaults -- see
scripts/species_param_map.json, built by parsing every line of the form
`species = p_eEC_X * rel0 * REF * dilution_factor;` /
`species = p_hEC_X * relN * REF * dilution_factor;` /
`species = vN * REF * dilution_factor;` in the model file):

  eEC_X  = p_eEC_X * rel0        * REF_X * dilution_factor
  hEC_X  = p_hEC_X * rel{1,2,3}  * 0.01  * dilution_factor   (rel1->hEC6411, rel2->hEC11719, rel3->hEC43117)
  metab  = v_i                   * REF_i * dilution_factor

All of Kcat_V_*, Keq_*, Km_*, Ki_*, Gnc_*, Gc_*, kdeg_*, v1-v74, and
rel0-rel3 come from parameters/final_params.pkl (trainparams[60][2], the
exact index used in the user's own example code -- loss=456, rank 27/100
by fit quality; a better fit exists at index 95 (loss=23) if the user
prefers). p_eEC*, p_hEC*, and dilution_factor are NOT in trainparams and
are set separately, exactly as in the user's example.

dilution_factor convention (CONFIRMED WITH USER, superseding an earlier
wrong guess of 10/D): the raw antimony parameter multiplies concentration
directly (dilution_factor UP => MORE concentrated), and
    dilution_factor_antimony = 1 / D_userfacing
where D_userfacing is CLAUDE.md's user-facing dilution factor (bounds
10-400, undiluted=10, larger=more dilute). So D_userfacing=10 (undiluted)
=> antimony param = 0.1, and D_userfacing=200 => antimony param = 0.005
(exactly matching the user's own "d=1/200 for a 200X dilution" example).
Note antimony param = 1.0 (the file's literal default) is NOT a physically
achievable/meaningful state -- some dilution is always needed to add
metabolites to the lysate, so D_userfacing=1 never occurs.

TIME UNITS: the user's own example simulates 4 hours as
`rfit.simulate(0, 4*60*60, ...)` i.e. SECONDS, not the minutes used by the
old model file's simulate(0, 4*60) convention. This harness uses SECONDS
throughout (duration_sec).

Heterologous/overexpression targets: p_hEC_X = target_mM / (0.01 * dilution_factor * rel_i)
eEC11137 (MDH) overexpression on top of native: p_eEC11137 = 1 + increment_mM / (REF_eEC11137 * dilution_factor)
CONFIRMED WITH USER: the 150 uM heterologous budget covers FOUR items, not
three -- hEC6411 (PYC), hEC11719 (FDH), hEC43117 (SDS), and the eEC11137
(MDH) overexpression INCREMENT (its native background does not count
against the budget, only the plasmid-driven "extra" does). Base-case split
is PYC=50, FDH=25, SDS=50, MDH increment=25 uM (50+25+50+25=150).
Exogenous metabolite spikes (added on top of diluted native lysate background):
  v_i_new = exogenous_addition_mM / (REF_i * dilution_factor) + v_i_trainparams_default
Endogenous enzyme knockdown/overexpression multiplier (incl. eEC2331/CS as the
hydroxycitrate proxy per CLAUDE.md): p_eEC_X = multiplier directly (native = 1.0).

Cumulative malate production via each route tracked exactly via extra ODE
species (same technique as before): Cum_rTCA_Malate' = piecewise(-R25, R25<0, 0),
Cum_TCA_Malate' = piecewise(-R20, R20<0, 0).
"""
import antimony
import roadrunner
import pickle
import json
import os
import numpy as np

_DIR = os.path.dirname(os.path.abspath(__file__))
MODEL_PATH = os.path.join(_DIR, "..", "model", "iMC057.txt")
PARAMS_PATH = os.path.join(_DIR, "..", "parameters", "final_params.pkl")
MAP_PATH = os.path.join(_DIR, "species_param_map.json")

TRAINPARAMS_INDEX = 60  # matches the user's own example code exactly

EXTRA_TRACKING = """
species Cum_rTCA_Malate = 0, Cum_TCA_Malate = 0;
Cum_rTCA_Malate' = piecewise(-R25, R25 < 0, 0);
Cum_TCA_Malate' = piecewise(-R20, R20 < 0, 0);
"""

EXOGENOUS_BASE_MM = {
    "ATP": 1.0,
    "AcetylCoA": 1.0,
    "HCO3": 10.0,
    "NADH": 1.0,
    "NADPH": 1.0,
    "Formate": 10.0,
    "L_Serine": 1.0,
}

# Hydroxycitrate (v73) is intentionally NOT in EXOGENOUS_BASE_MM / handled by the
# generic reciprocal-target formula. BUG FOUND AND FIXED this session: applying
# that formula to Hydroxycitrate (REF=1e-7) computed a v73 large enough to hit a
# REAL ~10 mM Hydroxycitrate concentration -- comparable to its own inhibition Ki
# (~10.4 mM), so NOT actually inert, contradicting both (a) CLAUDE.md's own note
# that hydroxycitrate has no apparent effect in the model, and (b) the validated
# reference notebook, which sets v73=10 literally (giving a negligible ~5e-9 mM
# Hydroxycitrate concentration) and relies on eEC2331 (citrate synthase) as the
# actual modeled lever instead. v73 is fixed at this literal placeholder value
# below, matching the notebook exactly, regardless of dilution or design.

HETEROLOGOUS_ENZYMES = ["hEC6411", "hEC11719", "hEC43117"]
HET_REL = {"hEC6411": "rel1", "hEC11719": "rel2", "hEC43117": "rel3"}
HET_REF = 0.01


# dilution_factor_antimony = 1 / D_userfacing (CONFIRMED with user 2024: "d=1/10 is
# undiluted", i.e. D_userfacing=10 -> raw param 0.1; raw param=1 is not physically
# achievable since some dilution is always needed to add metabolites).

with open(MODEL_PATH) as f:
    BASE_TEXT = f.read()
_FULL_TEXT = BASE_TEXT + EXTRA_TRACKING

with open(PARAMS_PATH, "rb") as f:
    _TRAINPARAMS = pickle.load(f)
FITTED_PARAMS = _TRAINPARAMS[TRAINPARAMS_INDEX][2]

with open(MAP_PATH) as f:
    SPECIES_MAP = json.load(f)  # species -> {driver, rel, ref}


def _fresh_runner():
    antimony.clearPreviousLoads()
    code = antimony.loadAntimonyString(_FULL_TEXT)
    if code < 0:
        raise RuntimeError(antimony.getLastError())
    mod = antimony.getMainModuleName()
    sbml = antimony.getSBMLString(mod)
    r = roadrunner.RoadRunner(sbml)
    for name, val in FITTED_PARAMS.items():
        r.model[name] = val
    return r


class Design:
    def __init__(
        self,
        dilution_factor_D=200.0,         # user-facing D (10-400, least-dilute bound=10; BASE BIOPRODUCTION SCHEME USES 200)
        duration_sec=14400.0,            # 4 h default
        het_enzymes_uM=None,             # dict override for hEC6411/hEC11719/hEC43117 (uM)
        eEC11137_overexpr_uM=25.0,       # additive MDH overexpression target (uM), drawn from the 150 uM budget
        endog_multipliers=None,          # dict p_eEC_X multiplier (native = 1.0)
        exogenous_overrides_mM=None,     # dict overriding EXOGENOUS_BASE_MM additions
        extra_species_overrides_mM=None, # dict absolute initial-concentration overrides (applied last)
    ):
        self.dilution_factor_D = dilution_factor_D
        self.duration_sec = duration_sec
        # 150 uM total heterologous budget: PYC(hEC6411)=50, FDH(hEC11719)=25,
        # SDS(hEC43117)=50, MDH(eEC11137) overexpression INCREMENT=25 -- confirmed
        # with user: MDH overexpression is plasmid-driven and draws from the SAME
        # 150 uM budget as the three hEC enzymes (50+25+50+25=150).
        self.het_enzymes_uM = {"hEC6411": 50.0, "hEC11719": 25.0, "hEC43117": 50.0}
        if het_enzymes_uM:
            self.het_enzymes_uM.update(het_enzymes_uM)
        self.eEC11137_overexpr_uM = eEC11137_overexpr_uM
        self.endog_multipliers = dict(endog_multipliers or {})
        self.exogenous_mM = dict(EXOGENOUS_BASE_MM)
        if exogenous_overrides_mM:
            self.exogenous_mM.update(exogenous_overrides_mM)
        self.extra_species_overrides_mM = dict(extra_species_overrides_mM or {})


def run_design(design: Design, n_points=400):
    r = _fresh_runner()

    dil = 1.0 / design.dilution_factor_D  # antimony dilution_factor param

    r.model["dilution_factor"] = dil

    # 1) endogenous enzymes: p_eEC_X = multiplier (native=1.0), except eEC11137 handled specially below
    for sp, info in SPECIES_MAP.items():
        if info["driver"].startswith("p_eEC") and sp != "eEC11137":
            mult = design.endog_multipliers.get(sp, 1.0)
            r.model[info["driver"]] = mult

    # 2) eEC11137 (MDH): native (p=1) contribution + additive overexpression target
    # BUG FOUND AND FIXED this session: this formula was missing the `* rel0`
    # divisor present in every other species formula (species = p * rel0 * REF *
    # dilution_factor), inflating p_eEC11137 by 1/rel0 (~5x) and cascading into
    # every downstream result. Confirmed against the validated manual/notebook
    # computation, which divides by (REF * dil * rel0), not just (REF * dil).
    mdh_info = SPECIES_MAP["eEC11137"]
    mdh_native_mult = design.endog_multipliers.get("eEC11137", 1.0)
    p_mdh = mdh_native_mult + (design.eEC11137_overexpr_uM / 1000.0) / (mdh_info["ref"] * dil * r.model["rel0"])
    r.model["p_eEC11137"] = p_mdh

    # 3) heterologous enzymes: absolute target concentration (mM), independent of native background
    for hid in HETEROLOGOUS_ENZYMES:
        target_mM = design.het_enzymes_uM.get(hid, 0.0) / 1000.0
        p_val = target_mM / (HET_REF * dil * r.model[HET_REL[hid]])
        r.model[f"p_{hid}"] = p_val

    # 4) metabolites: v_i = v_i_default (trainparams) [+ exogenous addition / (REF*dil)]
    for sp, info in SPECIES_MAP.items():
        if not info["driver"].startswith("v"):
            continue
        v_default = FITTED_PARAMS.get(info["driver"], None)
        if v_default is None:
            continue  # not part of fitted set; leave file default in place
        add_mM = design.exogenous_mM.get(sp, None)
        if add_mM is not None:
            r.model[info["driver"]] = v_default + add_mM / (info["ref"] * dil)
        else:
            r.model[info["driver"]] = v_default

    # Hydroxycitrate: literal placeholder, matching the validated reference
    # notebook exactly -- NOT solved for a real mM target (see EXOGENOUS_BASE_MM
    # comment above for why).
    r.model["v73"] = 10

    r.reset()

    # 5) final explicit absolute species overrides (applied after reset, then reset again)
    if design.extra_species_overrides_mM:
        for sid, val in design.extra_species_overrides_mM.items():
            r.model[sid] = val

    species_ids = list(r.model.getFloatingSpeciesIds())
    r.selections = ["time"] + species_ids + ["Cum_rTCA_Malate", "Cum_TCA_Malate", "R25", "R20"]
    result = r.simulate(0, design.duration_sec, max(200, n_points))

    res = np.array(result)
    cols = result.colnames
    idx = {c: i for i, c in enumerate(cols)}

    smalate_init = res[0, idx["SMalate"]]
    smalate_final = res[-1, idx["SMalate"]]
    cum_rtca = res[-1, idx["Cum_rTCA_Malate"]]
    cum_tca = res[-1, idx["Cum_TCA_Malate"]]

    return {
        "final_malate_mM": float(smalate_final),
        "initial_malate_mM": float(smalate_init),
        "delta_malate_mM": float(smalate_final - smalate_init),
        "cum_rTCA_malate_mM": float(cum_rtca),
        "cum_TCA_malate_mM": float(cum_tca),
        "rTCA_to_TCA_ratio": float(cum_rtca / cum_tca) if cum_tca > 1e-12 else float("inf"),
        "rTCA_fraction": float(cum_rtca / (cum_rtca + cum_tca)) if (cum_rtca + cum_tca) > 1e-12 else float("nan"),
        "raw_result": res,
        "colnames": cols,
        "final_state": {sid: float(res[-1, idx[sid]]) for sid in species_ids},
    }
