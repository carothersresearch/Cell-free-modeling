"""
Parameter-ensemble simulation engine for Figure_4.ipynb and Figure_5.ipynb.

Two ensemble sources are supported, both summarized in the notebooks as median + 25th-75th
percentile (IQR):

* 'jitter'  -- Monte Carlo parameter jitter around the representative parameter set
               (trainparams[60]). Each draw perturbs all fitted kinetic parameters
               (Kcat_V_*, Keq_*, Km_*, Ki_*, Gnc_*, Gc_*, kdeg_*) independently and
               log-normally: value *= 10**(z * sigma), z ~ N(0, 1), sigma = 1/1.645 (from
               CatPred's "90% of highest-confidence predictions are within one order of
               magnitude" -> 1.645*sigma = 1). Exogenous-target parameters (v*) and
               rel0-rel3 are not jittered. Draw i uses seed base_seed + i, so every
               ensemble is exactly reproducible.
* 'cluster' -- the independently fitted parameter sets of a trainparams cluster
               (e.g. cluster 4 = trainparams[37:63], which contains index 60).

Used by: Figure 4B (jitter), Figures 4A/4C/4D (cluster 4), Figures 5B-5G (jitter).

This is a real importable module (not notebook code) so that multiprocessing on Windows
(spawn start method) can pickle the worker function -- each worker process imports this
module and builds its own RoadRunner instance, so no RoadRunner object crosses a process
boundary. Paths are relative to the working directory, which must contain iMC057.txt and
filtering_parameter_sets/ (i.e. the PaperFigures directory).

Draws that fail to integrate (CVODE RuntimeError) are recorded as NaN. For Figure 5
readouts, numerically "successful" but nonphysical results (a concentration outside
[0, 100] mM, or |R20/R25 integrated amount| > 100 mM) are also recorded as NaN; these come
from a small fraction of aggressive jitter draws (observed magnitudes up to ~1e14).
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

MODEL_PATH = 'iMC057.txt'
PARAMS_PATH = 'filtering_parameter_sets/final_params.pkl'
BASE_INDEX = 60
CLUSTER4 = range(37, 63)  # trainparams[37:63], contains BASE_INDEX

SIGMA_LOG10 = 1.0 / 1.645  # ~0.6079
KINETIC_PREFIXES = ('Kcat_V_', 'Keq_', 'Km_', 'Ki_', 'Gnc_', 'Gc_', 'kdeg_')

PLAUSIBLE_CONC_MAX = 100.0   # mM
PLAUSIBLE_FLUX_MAX = 100.0   # mM, bound on |R20/R25 time-integrated amount|

D_BASE = 1 / 200


# ---------------------------------------------------------------------------
# Model / parameter loading (once per process)
# ---------------------------------------------------------------------------
def loada(path):
    antimony.clearPreviousLoads()
    code = antimony.loadAntimonyFile(path)
    if code < 0:
        raise RuntimeError(antimony.getLastError())
    return roadrunner.RoadRunner(antimony.getSBMLString(antimony.getMainModuleName()))


_RFIT = None
_TRAINPARAMS = None
_KINETIC_KEYS = None
_REACTION_IDS = None


def _get_rfit():
    global _RFIT, _TRAINPARAMS, _KINETIC_KEYS, _REACTION_IDS
    if _RFIT is None:
        _RFIT = loada(MODEL_PATH)
        with open(PARAMS_PATH, 'rb') as f:
            _TRAINPARAMS = pickle.load(f)
        base = _TRAINPARAMS[BASE_INDEX][2]
        _KINETIC_KEYS = [k for k in base if any(k.startswith(p) for p in KINETIC_PREFIXES)]
        _REACTION_IDS = _RFIT.model.getReactionIds()
    return _RFIT


def get_reaction_ids():
    _get_rfit()
    return _REACTION_IDS


def get_params(index):
    """Fitted parameter dict of trainparams[index]."""
    _get_rfit()
    return dict(_TRAINPARAMS[index][2])


def get_base_params():
    return get_params(BASE_INDEX)


def jittered_params(seed):
    """One jittered parameter dict around trainparams[BASE_INDEX]; fully determined by seed."""
    params = get_base_params()
    rng = np.random.default_rng(seed)
    for k in _KINETIC_KEYS:
        z = rng.standard_normal()
        params[k] = params[k] * (10 ** (z * SIGMA_LOG10))
    return params


def get_param_dicts_from_params(rfit, params):
    """Index-based {param: position} / {param: value} dicts for
    rfit.model.setGlobalParameterValues([*labels.values(), ...], [*values.values(), ...])."""
    labels, values = {}, {}
    for k, pid in enumerate(rfit.model.getGlobalParameterIds()):
        if pid in params:
            values[pid] = params[pid]
            labels[pid] = k
    return labels, values


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
# Figure 4: endogenous lysate (no heterologous enzymes), across dilution factors.
# Index-based setGlobalParameterValues([..., 2, 49, 51, 53], [..., d, 0, 0, 0]) zeroes
# dilution_factor -> d and p_hEC6411 / p_hEC11719 / p_hEC43117 -> 0.
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

# ---- Figure 4A: module flux fraction + total flux over time ----
D_FACTORS_4A = [1 / 200, 1 / 100, 1 / 50, 1 / 10]
TIME_POINTS_4A = np.arange(0, 4 * 60 * 60 + 1, 60 * 60)


def run_4a_draw(params):
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
            for j in range(n_steps):
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

            out[d] = {"frac": module_fluxes, "total_flux": total_flux_list}
        except RuntimeError:
            out[d] = {"frac": {m: [np.nan] * n_steps for m in PATHWAY_DICT},
                      "total_flux": [np.nan] * n_steps}
        rfit.resetToOrigin()
    return out


# ---- Figure 4B: total integrated |flux| per reaction, across 7 dilutions ----
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


def run_4b_draw(params):
    rfit = _get_rfit()
    n_reactions = len(get_reaction_ids())
    labels, values = get_param_dicts_from_params(rfit, params)
    member_flux = []
    for d in D_FACTORS_4B:
        try:
            rates_array = _simulate_and_record_rates_4b(rfit, d, labels, values)
            flux = (np.sum(np.abs(rates_array), axis=1) * 300).tolist()
        except RuntimeError:
            flux = [np.nan] * n_reactions
            rfit.resetToOrigin()
        member_flux.append(flux)
    return member_flux  # (7, n_reactions)


# ---- Figure 4C: integrated net flux per reaction, across 4 dilutions ----
D_FACTORS_4C = [1 / 10, 1 / 50, 1 / 100, 1 / 200]
TIME_POINTS_4C = np.arange(0, 4 * 60 * 60 + 1, 30 * 60)
DT_4C = 30 * 60


def run_4c_draw(params):
    rfit = _get_rfit()
    labels, values = get_param_dicts_from_params(rfit, params)
    member_rows = []
    for d in D_FACTORS_4C:
        try:
            vcopy = values.copy()
            vcopy['v73'] = 0
            rfit.model.setGlobalParameterValues(
                [*labels.values(), 2, 49, 51, 53],
                [*vcopy.values(), d, 0, 0, 0]
            )
            rfit.reset()
            rate_records = []
            for j in range(len(TIME_POINTS_4C) - 1):
                rfit.simulate(int(TIME_POINTS_4C[j]), int(TIME_POINTS_4C[j + 1]))
                rate_records.append(rfit.getReactionRates()[:57])
            integrated_flux = np.trapz(np.array(rate_records), dx=DT_4C, axis=0)
            member_rows.append(integrated_flux.tolist())
        except RuntimeError:
            member_rows.append([np.nan] * 57)
        rfit.resetToOrigin()
    return member_rows  # (4, 57)


# ---- Figure 4D: metabolite-spike sink analysis, 3 spikes x 2 dilutions ----
SPIKES_4D = {"Malate": "v56", "Pyruvate": "v53", "Serine": "v38"}
DILUTIONS_4D = {"10X": 1 / 10, "200X": 1 / 200}
SUBSET_IDX_4D = list(range(51, 81)) + list(range(82, 95)) + list(range(97, 122)) + [124, 125]
N_TIME_DEFAULT_4D = 51  # RoadRunner's default simulate(start, end) output resolution


def get_subset_names_4d():
    rfit = _get_rfit()
    probe = rfit.simulate(0, 1)
    ids = [c.strip('[]') for c in probe.colnames]
    rfit.resetToOrigin()
    return [ids[i] for i in SUBSET_IDX_4D]


def run_4d_draw(params):
    rfit = _get_rfit()
    labels, values_base = get_param_dicts_from_params(rfit, params)
    n_subset = len(SUBSET_IDX_4D)

    out = {}
    for spike_name, spike_param in SPIKES_4D.items():
        out[spike_name] = {}
        for d_name, d in DILUTIONS_4D.items():
            try:
                param_values = dict(values_base)

                # negative control (no spike)
                param_values[spike_param] = 0
                rfit.model.setGlobalParameterValues(
                    [*labels.values(), 2, 49, 51, 53],
                    [*param_values.values(), d, 0, 0, 0]
                )
                rfit.reset()
                sim_neg = rfit.simulate(0, 4 * 60 * 60)

                # 5 mM spike
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
                out[spike_name][d_name] = {
                    "tlist": [float(t) for t in sim_spike['time']],
                    "changes": change[:, SUBSET_IDX_4D].tolist(),
                }
            except RuntimeError:
                out[spike_name][d_name] = {
                    "tlist": [np.nan] * N_TIME_DEFAULT_4D,
                    "changes": [[np.nan] * n_subset for _ in range(N_TIME_DEFAULT_4D)],
                }
            rfit.resetToOrigin()
    return out


# ---------------------------------------------------------------------------
# Figure 5: heterologous malate pathway. Base recipe (mM): +1 NADH, +1 NADPH, +10 formate,
# +1 ATP, +1 acetyl-CoA, +10 HCO3, +1 substrate; enzymes PYC=0.050, FDH=0.025, SDS=0.050,
# MDH overexpression increment=0.025 (150 uM total heterologous expression budget).
# ---------------------------------------------------------------------------
def _native_background(rfit, d, params):
    _apply_params(rfit, params)
    rfit.model['dilution_factor'] = d
    rfit.model['p_hEC6411'] = 0
    rfit.model['p_hEC11719'] = 0
    rfit.model['p_hEC43117'] = 0
    rfit.reset()
    rel = (rfit.model['rel0'], rfit.model['rel1'], rfit.model['rel2'], rfit.model['rel3'])
    pre = {s: rfit.getValue(s) for s in
           ['NADH', 'NADPH', 'Formate', 'ATP', 'AcetylCoA', 'HCO3', 'eEC11137',
            'Pyruvate', 'Oxaloacetate', 'L_Serine']}
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


def _final_malate(rfit):
    rfit.reset()
    try:
        results = rfit.simulate(0, 4 * 60 * 60, selections=['time', 'SMalate'])
        val = _check_conc(results[-1][1])
    except RuntimeError:
        val = np.nan
    rfit.resetToOrigin()
    return val


# ---- Figure 5B: 8 named conditions; enzymes = (pyc, fdh, sds, mdh increment) in mM ----
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
        out[label] = _final_malate(rfit)
    return out


# ---- Figure 5C: enzyme titrations (4 enzymes x 3 concs x 3 substrates) ----
ENZYMES_5C = ['fdh', 'mdh', 'pyc', 'sds']
CONCS_5C = [0, 0.001, 0.005]
SUBSTRATES_5C = ['Oxaloacetate', 'Pyruvate', 'Serine']
BASE_DOSE_5C = {'fdh': 0.025, 'mdh': 0.025, 'pyc': 0.050, 'sds': 0.050}


def run_5c_draw(params, d=D_BASE):
    rfit = _get_rfit()
    out = {}
    for substrate in SUBSTRATES_5C:
        for enzyme in ENZYMES_5C:
            for conc in CONCS_5C:
                rel, pre = _native_background(rfit, d, params)
                dose = {k: (conc if k == enzyme else BASE_DOSE_5C[k]) for k in ENZYMES_5C}
                _set_common_cofactors(rfit, d, pre, substrate=substrate, substrate_conc=1.0)
                _set_enzymes(rfit, d, rel, pre, dose['pyc'], dose['fdh'], dose['sds'], dose['mdh'])
                out[f"{substrate}|{enzyme}|{conc}"] = _final_malate(rfit)
    return out


# ---- Figure 5D: substrate titrations, full base-case enzyme set ----
TITRATION_CONCS_5D = [0, 0.01, 0.5, 1, 5, 10, 20, 50, 100, 1000]
SUBSTRATES_5D = ['Serine', 'Pyruvate', 'Oxaloacetate']


def run_5d_draw(params, d=D_BASE):
    rfit = _get_rfit()
    out = {}
    for substrate in SUBSTRATES_5D:
        for conc in TITRATION_CONCS_5D:
            rel, pre = _native_background(rfit, d, params)
            _set_common_cofactors(rfit, d, pre, substrate=substrate, substrate_conc=conc)
            _set_enzymes(rfit, d, rel, pre, 0.050, 0.025, 0.050, 0.025)
            out[f"{substrate}|{conc}"] = _final_malate(rfit)
    return out


# ---- Figure 5E: dilution titration (Serine-fed full pathway): final malate, final/initial CO2 ----
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


# ---- Figure 5F: malate produced via R20 (TCA, fumarase) vs R25 (rTCA, MDH) across dilution ----
DILFACS_5F = [10, 20, 50, 100, 200]
N_TIME_POINTS_5F = 961  # 15 s resolution over 4 h
ENZYMES_5F_SERINE = (0.050, 0.025, 0.050, 0.025)
ENZYMES_5F_PYRUVATE = (0.050, 0.025, 0, 0.025)  # SDS not needed without serine

# Heterologous enzymes (and the MDH increment) are expressed during the undiluted TXTL
# phase, so their target doses are defined at D=10 and then dilute along with the rest of
# the lysate for biosynthesis at the simulated D (concentration scales as 10/D).
D_UNDILUTED_FACTOR = 1 / 10


def _set_enzymes_expressed_undiluted(rfit, d, rel, pre, hEC6411, hEC11719, hEC43117, eEC11137_incr):
    rel0, rel1, rel2, rel3 = rel
    rfit.model['dilution_factor'] = d
    rfit.model['p_hEC6411'] = hEC6411 / (0.01 * D_UNDILUTED_FACTOR * rel1) if hEC6411 else 0
    rfit.model['p_hEC11719'] = hEC11719 / (0.01 * D_UNDILUTED_FACTOR * rel2) if hEC11719 else 0
    rfit.model['p_hEC43117'] = hEC43117 / (0.01 * D_UNDILUTED_FACTOR * rel3) if hEC43117 else 0
    if eEC11137_incr:
        rfit.model['p_eEC11137'] = (eEC11137_incr / (0.012601796 * D_UNDILUTED_FACTOR * rel0)) + (pre['eEC11137'] / (0.012601796 * d * rel0))


def _run_5f(params, substrate, enzymes):
    """Returns {dilfac: (R20 amount, R25 amount)}, each the time integral of the reaction's
    malate-producing (negative) net rate over 4 h."""
    rfit = _get_rfit()
    out = {}
    for dilfac in DILFACS_5F:
        d_local = 1 / dilfac
        rel, pre = _native_background(rfit, d_local, params)
        _set_common_cofactors(rfit, d_local, pre, substrate=substrate, substrate_conc=1.0)
        _set_enzymes_expressed_undiluted(rfit, d_local, rel, pre, *enzymes)
        rfit.reset()
        try:
            results = rfit.simulate(0, 4 * 60 * 60, N_TIME_POINTS_5F, selections=['time', 'R20', 'R25'])
            t = np.asarray(results['time'])
            r20 = _check_flux(np.trapz(-np.asarray(results['R20']), t))
            r25 = _check_flux(np.trapz(-np.asarray(results['R25']), t))
        except RuntimeError:
            r20, r25 = np.nan, np.nan
        rfit.resetToOrigin()
        out[dilfac] = (r20, r25)
    return out


def run_5f_draw(params):
    return _run_5f(params, 'Serine', ENZYMES_5F_SERINE)


def run_5f_pyruvate_draw(params):
    return _run_5f(params, 'Pyruvate', ENZYMES_5F_PYRUVATE)


# ---- Figure 5G: MCA-inspired control coefficients (+/-100-fold initial concentration) ----
_LOOKUP = None


def get_species_lookup():
    """species -> {'v': its v_i parameter, 'val': its REF constant}, parsed from the model."""
    global _LOOKUP
    if _LOOKUP is None:
        import re
        with open(MODEL_PATH, 'r') as f:
            text = f.read()
        pattern = re.compile(r"(?P<species>\S+)\s*=\s*(?P<v>v\d+)\*(?P<val>[0-9.eE+-]+)\*dilution_factor;")
        _LOOKUP = {}
        for match in pattern.finditer(text):
            try:
                val = float(match.group("val"))
            except ValueError:
                val = None
            _LOOKUP[match.group("species")] = {"v": match.group("v"), "val": val}
    return _LOOKUP


def _run_5g_perturbation(rfit, d, params, s, direction, lookup):
    rel, pre = _native_background(rfit, d, params)
    pre_s = rfit.getValue(s)
    _set_common_cofactors(rfit, d, pre, substrate='Serine', substrate_conc=1.0)
    multiplier = 100 if direction > 0 else (1 / 100)
    rfit.model[lookup[s]["v"]] = pre_s * multiplier / lookup[s]["val"] / d
    _set_enzymes(rfit, d, rel, pre, 0.050, 0.025, 0.050, 0.025)
    return _final_malate(rfit)


def run_5g_draw(params, species_list, d=D_BASE):
    rfit = _get_rfit()
    lookup = get_species_lookup()
    out = {}
    for s in species_list:
        malate_up = _run_5g_perturbation(rfit, d, params, s, +1, lookup)
        malate_down = _run_5g_perturbation(rfit, d, params, s, -1, lookup)
        if np.isnan(malate_up) or np.isnan(malate_down) or malate_up <= 0 or malate_down <= 0:
            out[s] = np.nan
        else:
            out[s] = (np.log(malate_up) - np.log(malate_down)) / np.log(1e4)
    return out


def get_species_5g():
    """The 71 tracked species perturbed in Figure 5G (same id slice as the original notebook)."""
    rfit = _get_rfit()
    lookup = get_species_lookup()
    return [s for s in rfit.getIds()[50:121] if s in lookup]


# ---------------------------------------------------------------------------
# Parallel ensemble driver
# ---------------------------------------------------------------------------
PANELS = {
    '4a': run_4a_draw, '4b': run_4b_draw, '4c': run_4c_draw, '4d': run_4d_draw,
    '5b': run_5b_draw, '5c': run_5c_draw, '5d': run_5d_draw, '5e': run_5e_draw,
    '5f': run_5f_draw, '5f_pyruvate': run_5f_pyruvate_draw, '5g': run_5g_draw,
}


def _worker(task):
    panel, source, key, kwargs = task
    params = jittered_params(key) if source == 'jitter' else get_params(key)
    return PANELS[panel](params, **kwargs)


def run_ensemble(panel, source='jitter', n_draws=30, base_seed=1000, indices=CLUSTER4,
                 n_workers=6, **kwargs):
    """Run one panel's simulation across an ensemble, in parallel.

    source='jitter':  n_draws jittered draws around trainparams[BASE_INDEX], seeds
                      base_seed, base_seed+1, ...
    source='cluster': one run per trainparams index in `indices` (default cluster 4).
    Extra kwargs are passed to the panel function (e.g. species_list for '5g').
    Returns a list with one result per ensemble member, in order.
    """
    import multiprocessing as mp
    keys = [base_seed + i for i in range(n_draws)] if source == 'jitter' else list(indices)
    tasks = [(panel, source, k, kwargs) for k in keys]
    with mp.get_context('spawn').Pool(processes=n_workers) as pool:
        return pool.map(_worker, tasks)
