# Metabolic Model Simulation Project

## Overview
This project uses Tellurium to simulate a metabolic model stored as Antimony/SBML files in `model/`. This Antimony file contains a kinetic metabolic model of E. coli derived crude lysate. The model is meant to be representative of the endogenous metabolism, as well as inform engineering strategies of a heterologous bioproduction pathway to produce malate using rTCA. To that point, the model contains three heterologous enzymes, hEC6411 (pyruvate carboxylase, PYC), hEC11719 (formate dehydrogenase, FDH), and hEC43117 (serine deaminase/ammonia lyase, SDS). The endogenous enzyme malate dehydrogenase (eEC11137, MDH) is also overexpressed in the engineered pathway, from a plasmid. In the model file, the heterologous pathway enzyme concentrations are set to zero, to reflect the native conditions, but in reality this cell free system is capable of producing approximately 150 uM of protein total, which would need to be split amongst the heterologous proteins. **The 150 uM budget covers FOUR items, not three: PYC (hEC6411), FDH (hEC11719), SDS (hEC43117), and the MDH (eEC11137) overexpression INCREMENT — MDH is plasmid-driven, so the "extra" MDH on top of its small native background draws from the same 150 uM heterologous expression budget as the other three.** In the base case, MDH and FDH are expressed half as much as PYC and SDS: PYC=50 uM, SDS=50 uM, FDH=25 uM, MDH overexpression increment=25 uM (50+50+25+25=150 uM exactly). This pathway is intended to hook up to a separate bioproduction pathway which generates Serine from formate (CO2 equivalent). All concentration values are in mM. The concentrations of endogenous enzymes are tunable only by engineering the lysate. Readily tunable parameters in the system include the dilution factor of the crude lysate, all heterologously expressed enzymes, and the small molecule inhibitor, hydroxycitrate. Additionally, the starting concentrations of nearly all metabolites could be adjusted in vitro. In our base bioproduction schema, we only add in exogenous 1 mM ATP, 1 mM acetyl-CoA, 10mM sodium bicarbonate, 10 mM hydroxycitrate, 1 mM NADH, 1 mM NADPH, 10 mM formate, and 1 mM serine. Similarly to the heterologous enzymes, the values for concentrations of these exogenously added metabolites in the model does not reflect the experimentally added concentrations, only the concentration of that metabolite within the plain crude lysate. It should be noted that the model is parameterized with our best guesses, and we have been unable to recapitulate all experimentally observed results using the model as is.

**Current model file: `model/iMC057.txt`** (an earlier file, `model/cell_free_model.txt`, was used in a prior iteration of this project and is now retired — do not use it. The reaction network and the R20/R25 semantics described below are unchanged between the two files, but `iMC057.txt` is the one wired up to `parameters/final_params.pkl` and is the one described by every section below.)

## Variables you can change
- Concentrations of exogenous cofactors and metabolites.
- Heterologous enzyme expression levels (remember you have a budget of total heterologous expression capacity of 150 uM of protein, covering PYC + FDH + SDS + the MDH overexpression increment — see Overview).
- Knockdown or overexpression of endogenous enzymes.
- Dilution of the lysate. **The lysate can be diluted no more than 10-fold less than the base recipe — dilution_factor=10 is the LEAST DILUTE the lysate can physically be (a bound, not a recipe you'd normally run), since some dilution is always necessary to add metabolites.** Larger values of dilution_factor indicate more dilute lysate (ie. dilution_factor=200 represents a 200X dilution of the lysate). The bounds of this variable are 10-400. **The base bioproduction scheme itself uses dilution_factor=200** — do not treat dilution_factor=10 as the base case; it is the most-concentrated end of the allowed exploration range, and is itself a strong lever (see "Key implementation notes" below — reducing dilution from the 200x base case toward the 10x floor is, by a wide margin, the single most impactful lever found in this system).
- Addition of other metabolites found in the model but not currently added exogenousely.
- Inhibition of citrate synthase. While the hydroxycitrate small molecule inhibitor does not seem to have an apparent effect in the model (its concentration barely affects anything through the wiring described below — do not spend design effort trying to set Hydroxycitrate's own concentration to a specific target mM value, it isn't a meaningful lever via that route), literature and prior experiments prove that addition of hydroxycitrate strongly decreases the activity of citrate synthase in all directions. Change the concentration of citrate synthase (eEC2331) as a way to simulate the addition of hydroxycitrate.
- Duration of experiment. Our baseline conditions simulate biosynthesis for four hours. Your bounds are 30 minutes to 24 hours. **Note the model's time unit is SECONDS (see below) — 4 hours = 14400, 30 min = 1800, 24 h = 86400.**

Remember when choosing designs that more interventions can lead to more work for the engineer implementing these in vitro, and the cost of substrate can significantly increase the cost of the production platform. Do your best to keep designs feasible, both economically and experimentally.

## Objectives
We want to optimize the production of malate through a carbon-conserving route. This means that we want to produce as much malate as possible, but we also want to maximize the ratio of malate being produced through reaction R25 (representing rTCA, via eEC11137/MDH running in reverse: Oxaloacetate + NADH -> SMalate + NAD) compared to through reaction R20 (representing TCA, via eEC4212/fumarase running in reverse: Fumarate + H2O -> SMalate). The malate species is `SMalate` in the model. R25 and R20 are the only two reactions that produce or consume SMalate. Each is written as a single net reaction (`R25 = R25_f - R25_r`, positive = malate-CONSUMING direction), so cumulative malate PRODUCED via each route over a simulation must be tracked by integrating the negative-going part of each net rate — see the working pattern in "Key implementation notes" below.

## Outputs
Return a structured JSON file detailing the conditions you changed from the base case. This needs to be machine readable so that I can load it with Python and re-simulate the model myself to validate your simulations.

Return your simulated values for final malate concentration, ratio of cumulative malate production from rTCA vs TCA, and the difference in malate concentration from start to end.

For each design, I want a structured output of Design ID, the intervention(s), the variables changed and their proposed value(s), the rationale, and expected benefit and tradeoff.

## Model format
Models are in Antimony format. The documented load pattern is:
```python
import tellurium as te
r = te.loada("model/iMC057.txt")
```
**However, `import tellurium` currently fails in this project's default (anaconda) Python environment** — `libroadrunner` is compiled against numpy<2 and this environment has numpy 2.x, causing a DLL init failure; separately, `tellurium`'s optional `rrplugins` submodule fails to import even when numpy is downgraded. Workaround, validated this session:
1. Create an isolated venv with numpy<2: `python -m venv .venv_sim` then `.venv_sim/Scripts/python.exe -m pip install "numpy<2" tellurium` (do this once; do not touch the system anaconda numpy, many other packages depend on numpy 2.x).
2. Bypass the broken `tellurium` wrapper and use `antimony` + `roadrunner` directly (this is all `te.loada` does internally anyway):
```python
import antimony, roadrunner

def loada(path):
    antimony.clearPreviousLoads()
    code = antimony.loadAntimonyFile(path)
    if code < 0:
        raise RuntimeError(antimony.getLastError())
    module_name = antimony.getMainModuleName()
    sbml = antimony.getSBMLString(module_name)
    return roadrunner.RoadRunner(sbml)

r = loada("model/iMC057.txt")
```
A full working harness built on this pattern (with the fitted-parameter loading, dilution/enzyme/metabolite formulas, and cumulative rTCA/TCA tracking all implemented and validated) exists at `scripts/sim2.py` — prefer reusing/extending that file over re-deriving this from scratch.

## Parameterization mechanism (how `model/iMC057.txt` actually works)

Unlike a simpler Antimony file where species initial values might just be literal numbers, every species in `model/iMC057.txt` is computed from a formula wired into the model itself:
```
eEC_X  = p_eEC_X * rel0        * REF_X * dilution_factor      # endogenous enzymes (rel0 shared by all of them)
hEC_X  = p_hEC_X * rel{1,2,3}  * 0.01  * dilution_factor      # heterologous enzymes (rel1->hEC6411, rel2->hEC11719, rel3->hEC43117)
metab  = v_i                   * REF_i * dilution_factor      # metabolites
```
`REF_X`/`REF_i` are per-species constants baked into the model file (parsed once this session into `scripts/species_param_map.json` and cross-validated against roadrunner's own computed defaults with zero mismatches — reuse that file rather than re-parsing). `p_eEC*`, `p_hEC*`, `dilution_factor`, `rel0-rel3`, and `v1-v74` are all ordinary global parameters you set directly (e.g. `r.model['p_eEC4212'] = 0.1`, or via `r.model.setGlobalParameterValues([...], [...])` by index as in the example below — both are equivalent).

**Fitted kinetic parameters** (`Kcat_V_*`, `Keq_*`, `Km_*`, `Ki_*`, `Gnc_*`, `Gc_*`, `kdeg_*`, plus `v1-v74` and `rel0-rel3`) are NOT the file's own literal defaults — they must be loaded from `parameters/final_params.pkl`:
```python
import pickle
with open('parameters/final_params.pkl', 'rb') as f:
    trainparams = pickle.load(f)
fitted = trainparams[60][2]   # trainparams is a list of 100 (id, loss, params_dict) tuples
for name, val in fitted.items():
    r.model[name] = val
```
Index 60 is the one used in this project's own analysis notebooks and should be used by default for consistency. Note it is NOT the best-fitting set in the pickle (loss=456, rank 27 of 100 by fit quality — index 95 has the lowest loss, 23, if a better fit is ever wanted instead). `p_eEC*`, `p_hEC*`, and `dilution_factor` are NOT part of `trainparams` and must always be set separately (as below).

**dilution_factor convention** — this tripped up analysis multiple times this session, get it right: the raw antimony parameter `dilution_factor` multiplies concentration directly (bigger value = MORE concentrated), and
```
dilution_factor_antimony = 1 / D
```
where `D` is the "Variables you can change" dilution factor described above (bounds 10-400). So `D=10` (least dilute) → `dilution_factor_antimony=0.1`; the actual base case `D=200` → `dilution_factor_antimony=0.005`. A raw value of 1.0 (the file's own literal default) does not correspond to any physically achievable `D` and should never be used directly.

**Time units are SECONDS**, not minutes — confirmed directly. 4 hours = `r.simulate(0, 14400)`.

**Hydroxycitrate / v73**: setting `v73` doesn't meaningfully change `Hydroxycitrate`'s actual modeled concentration (its `REF` constant is `1e-7`, so any reasonable `v73` value produces a near-zero concentration) — this is consistent with CLAUDE.md's existing note that hydroxycitrate's own concentration has no apparent effect in the model. Don't try to hit a specific mM target for `Hydroxycitrate` via `v73`; the working notebook pattern just sets `v73=10` as an inert placeholder and relies on `p_eEC2331` (citrate synthase) as the actual, effective lever instead.

### Validated end-to-end example (base case, D=200)

This reproduces the project's own reference notebook to 8-9 significant figures (validated twice this session against two different versions of it) and reflects the corrected base-case enzyme/dilution values:

```python
import antimony, roadrunner, pickle

antimony.clearPreviousLoads()
antimony.loadAntimonyFile('model/iMC057.txt')
sbml = antimony.getSBMLString(antimony.getMainModuleName())
r = roadrunner.RoadRunner(sbml)

with open('parameters/final_params.pkl', 'rb') as f:
    trainparams = pickle.load(f)
for name, val in trainparams[60][2].items():
    r.model[name] = val

D = 200          # base bioproduction scheme dilution factor (NOT 10 -- see above)
d = 1 / D

# first pass: apply dilution, zero out heterologous enzymes, to read off the
# native ("pre") background at this dilution before adding anything exogenous
r.model['dilution_factor'] = d
r.model['p_hEC6411'] = 0
r.model['p_hEC11719'] = 0
r.model['p_hEC43117'] = 0
r.reset()

rel0, rel1, rel2, rel3 = r.rel0, r.rel1, r.rel2, r.rel3
preNADH, preNADPH   = r.getValue('NADH'), r.getValue('NADPH')
preFormate, preATP  = r.getValue('Formate'), r.getValue('ATP')
preAcCoA, preHCO3    = r.getValue('AcetylCoA'), r.getValue('HCO3')
preMdh, preSerine    = r.getValue('eEC11137'), r.getValue('L_Serine')

# base-case heterologous/MDH targets (mM) -- 50/25/50 uM PYC/FDH/SDS + 25 uM MDH increment
hEC6411, hEC11719, hEC43117, eEC11137_increment = 0.050, 0.025, 0.050, 0.025

p_hEC6411  = hEC6411  / (0.01 * d * rel1)
p_hEC11719 = hEC11719 / (0.01 * d * rel2)
p_hEC43117 = hEC43117 / (0.01 * d * rel3)
p_eEC11137 = (eEC11137_increment / (0.012601796 * d * rel0)) + (preMdh / (0.012601796 * d * rel0))

r.model['v73'] = 10  # inert placeholder, see note above
r.model['v44'] = (1 / 0.083 / d) + (preNADH / 0.083 / d)      # NADH,   +1 mM exogenous
r.model['v46'] = (1 / 0.12 / d)  + (preNADPH / 0.12 / d)      # NADPH,  +1 mM exogenous
r.model['v74'] = (10 / 0.001 / d) + (preFormate / 0.001 / d)  # Formate,+10 mM exogenous
r.model['v6']  = (1 / 1.5 / d)   + (preATP / 1.5 / d)         # ATP,    +1 mM exogenous
r.model['v10'] = (1 / 0.61 / d)  + (preAcCoA / 0.61 / d)      # AcetylCoA, +1 mM exogenous
r.model['v35'] = (10 / 0.015 / d) + (preHCO3 / 0.015 / d)     # HCO3,   +10 mM exogenous
r.model['v38'] = (1 / 0.001 / d) + (preSerine / 0.001 / d)    # Serine, +1 mM exogenous

r.model['dilution_factor'] = d
r.model['p_eEC11137']  = p_eEC11137
r.model['p_hEC6411']   = p_hEC6411
r.model['p_hEC11719']  = p_hEC11719
r.model['p_hEC43117']  = p_hEC43117

r.reset()
results = r.simulate(0, 4 * 60 * 60)   # 4 hours, in seconds
print(results['[SMalate]'][0], results['[SMalate]'][-1])
# -> 0.001056..., 0.106795...  (base case: final malate ~0.107 mM)
```

To change any endogenous enzyme's expression (knockdown or overexpression, e.g. eEC4212/fumarase or eEC2331/citrate synthase), set `r.model['p_eEC4212'] = <multiplier>` where 1.0 = native/unmodified — no dilution-aware math needed, the formula already handles it. To change a metabolite's exogenous addition amount, repeat the `preX`/`v_i` pattern above with a different mM target. To change dilution alone, change `D` (and therefore `d`) and recompute every `v_i`/`p_hEC*`/`p_eEC11137` line above, since they all depend on `d`.

### Tracking cumulative malate production via R25 (rTCA) vs R20 (TCA)

Append two extra species + rate rules to the Antimony source text before loading (referencing the reaction fluxes `R25`/`R20` directly by name, which Antimony permits):
```
species Cum_rTCA_Malate = 0, Cum_TCA_Malate = 0;
Cum_rTCA_Malate' = piecewise(-R25, R25 < 0, 0);
Cum_TCA_Malate' = piecewise(-R20, R20 < 0, 0);
```
This gives an exact ODE-integrated (not post-hoc trapezoidal) measure of gross malate produced via each route; `Cum_rTCA_Malate[-1] / Cum_TCA_Malate[-1]` is the ratio the Objectives section asks for. A full working implementation of this, plus the whole parameterization pattern above wrapped in a reusable `Design`/`run_design` API, is in `scripts/sim2.py`.

## Key implementation notes / corrections learned (read before starting a new analysis)

These were each discovered and validated the hard way this session — get them right from the start next time:

1. **Base case dilution is D=200, not D=10.** D=10 is only the least-dilute bound of the allowed 10-400 range. This was the single most consequential mistake made this session: it changes base-case yield by roughly 8x and is by far the most powerful lever in the whole design space (reducing dilution from 200 toward 10 alone increases malate yield ~7-8x, at only a modest selectivity cost).
2. **The 150 uM heterologous budget covers PYC+FDH+SDS+MDH-increment (four items), not just the three hEC* enzymes.** Base split: PYC=50, SDS=50, FDH=25, MDH increment=25 uM.
3. **Time unit is seconds, not minutes.** `r.simulate(0, 4*60*60)` for 4 hours.
4. **`dilution_factor_antimony = 1/D`**, not `D/10` or `10/D` — confirmed directly; do not re-derive this from the model file's own unhelpful literal default (1.0), which corresponds to no physically achievable `D`.
5. Selectivity (rTCA:TCA ratio) is EXTREMELY sensitive to dilution: at the base D=200 the TCA/fumarase route is almost completely suppressed (ratio in the millions); at D=10 it's still excellent but much more finite (~37:1 in the unmodified base recipe). This creates a genuine yield-vs-selectivity tradeoff axis via dilution alone that's worth exploring explicitly in any new design work, not just treating dilution=10 as a free win.
6. Citrate synthase (eEC2331) and fumarase (eEC4212) knockdowns are cheap (no reagent cost), low-risk, and reliably improve the rTCA:TCA ratio; their effect on raw yield is regime-dependent (sometimes small, sometimes a meaningful +20-30% depending on the dilution/enzyme setting) — always test empirically for the specific condition being evaluated rather than assuming a fixed effect size.
7. Serine dose-response is strongly positive and often the single largest yield lever available, especially at lower dilution (more concentrated lysate) — but is a real, scaling reagent cost, unlike the enzyme knockdowns.
8. Extra NADH/NADPH and extra HCO3/Formate supplementation (beyond the base 1x/10x recipe) have historically shown inconsistent, regime-dependent effects (sometimes mildly positive, sometimes negative) across different (mistaken) base-case parameterizations tried this session — do not assume either direction; re-test whenever the base case changes.
9. When in doubt about what "base case" numerically means, ask rather than infer — this session got the dilution value wrong twice in a row (once assuming D=10 was base, once by not distinguishing "d=1/10 is undiluted" from "the base recipe uses d=1/10") before the user gave the explicit, unambiguous answer (D=200 for dilution, 50/50/25/25 uM for the four expression items). Cross-check any new "base case" assumption against an actual reference run before building an extensive analysis on top of it.
10. **Hydroxycitrate (`v73`) must be a literal placeholder (`v73=10`), never solved for a real mM target.** An early version of `scripts/sim2.py` applied the same generic exogenous-target formula used for ATP/serine/etc. to Hydroxycitrate too — since its `REF` constant is `1e-7`, hitting a literal 10 mM target actually produced a real 10 mM Hydroxycitrate concentration, comparable to its own inhibition Ki (~10.4 mM), so it was NOT inert as intended and quietly suppressed citrate synthase/aconitase in every design (base-case yield was off by ~2.3% from this alone). Always special-case `Hydroxycitrate`/`v73` — set it to the literal placeholder value used by the reference notebook, and use `p_eEC2331` (citrate synthase) as the actual lever, per CLAUDE.md's own guidance.
11. **Every species formula includes a `rel` term in the denominator when solving for a target concentration — do not drop it.** `species = p * rel_N * REF * dilution_factor` means solving for `p` given a target concentration requires dividing by `rel_N` too (`p = target / (REF * dilution_factor * rel_N)`), not just by `REF * dilution_factor`. An early version of `scripts/sim2.py` dropped `* rel0` from the eEC11137 (MDH) overexpression formula specifically, silently delivering only ~20% (`rel0`) of the intended overexpression increment in every design. This class of bug is easy to introduce when writing one formula per enzyme by hand — when adding a new lever, cross-check its formula against the validated worked example above (or against `scripts/species_param_map.json` + the general formula) rather than deriving it fresh.
12. **After writing up a session's example code for documentation purposes, actually run it and check the printed value against an independent computation before trusting it.** Both bugs above were only caught because the worked example in this file was executed and its output compared against a hand-verified reference computation, rather than assumed correct because `scripts/sim2.py`'s own earlier design outputs looked internally self-consistent. Internal consistency (a script agreeing with itself across many runs) is not validation — a script can be internally consistent and wrong.

## Working files from this session
- `scripts/sim2.py` — validated simulation harness (`Design` class + `run_design` function) implementing everything above, including cumulative rTCA/TCA tracking. Both bugs in note 10/11 above are fixed as of the current version.
- `scripts/species_param_map.json` — parsed, validated species-to-parameter-formula mapping for every species in `model/iMC057.txt`.
- `scripts/run_designs5.py` — reproduces the latest 5-design analysis end to end (supersedes `run_designs4.py`, `run_designs3.py`, `run_designs2.py`, `run_designs.py`, all of which predate one or both bug fixes above).
- `results/malate_rtca_designs_v5.json` — latest design analysis output (structured per the Outputs section above), including a written log of exactly what was wrong in each prior version and why. Earlier versions (`malate_rtca_designs.json`, `_v2`, `_v3`, `_v4`) in `results/` are all superseded and kept only for history — do not use their numbers.

## Heterologous pathway diversification (10 substrates x 3 pathways, FRENDA-parameterized)

A separate sub-project (not the malate/rTCA work above): 30 new heterologous pathway
designs (10 starting substrates x 3 pathways each) that plug into the same base
crude-lysate model, screened for product yield within the same 150 uM heterologous
expression budget. Unlike the malate pathway, these use **real kinetics fetched via
FRENDA** (`CellFreeModeling/`, see that directory's own CLAUDE.md) rather than
hand-fitted trainparams -- BRENDA-experimental or CatPred-predicted Kcat/Km per EC
number, eQuilibrator Keq, converted into the same Common Modular Rate Law Antimony
style as `model/iMC057.txt` (see the R25/R55 pattern in that file).

### To resimulate or change parameters in a new session

```bash
cd LLM_CFS_2
.venv_sim/Scripts/python.exe scripts/run_pathway_sims.py       # regenerates all 30 model/pathways/*.txt AND resimulates them
.venv_sim/Scripts/python.exe scripts/compile_results_table.py  # rebuilds results/pathway_designs_table.csv + .json from the sim output
```
(`run_pathway_sims.py <N>` regenerates+resimulates only the first N pathways -- useful for quick iteration, but running it with no argument is required to refresh the full `results/pathway_sim_results.json` before recompiling the table, or the table will silently reflect a partial batch.)

**Where to change what:**
- Per-pathway topology (enzymes, substrate, dose, cofactor additions, module, notes, MSP): `scripts/pathway_designs.py` (the `PATHWAYS` list, built via the `add(...)` helper). The `HCO3_2_3HP` design is special-cased in `scripts/pathway_generator.py::get_enzyme_steps` (split into two real EC-backed steps, MCR_a/MCR_b) rather than read from its `PATHWAYS` entry directly.
- enzyme_name -> EC number mapping: `scripts/ec_map.py` (`ENZYME_TO_EC`). Also holds `COLLIDING_ECS`/`RESERVED_LABELS` (6.4.1.1 and 1.17.1.9 collide with the base model's own hEC6411/hEC11719 -- the generator auto-appends an "X" suffix to keep new pathways' enzyme pools isolated from the original malate pathway's R55/R56 reactions).
- Real kinetics source: `parameters/pathway_frenda_reactions.csv` (copied from `CellFreeModeling/my_pathways_output/Reaction.csv` for durability -- the FRENDA tool's own output directory is not guaranteed stable across unrelated future FRENDA runs). Parsed by `scripts/frenda_kinetics.py::load_kinetics()`. To re-fetch from FRENDA with different ECs/organisms, rerun FRENDA on `parameters/pathway_frenda_proteome_input.csv` (or an edited copy) and overwrite `pathway_frenda_reactions.csv` with the new `Reaction.csv`.
- Manual kinetics overrides (for ECs FRENDA/KEGG can't resolve, e.g. 1.17.1.9 formate dehydrogenase reuses the base model's own fitted hEC11719 values, 2.3.1.264 PHB synthase and the spontaneous formaldehyde+THF condensation use literature-typical estimates): `scripts/frenda_kinetics.py::MANUAL_KINETICS_OVERRIDE`. Also in that file: `REASONABLE_KM_DEFAULT`/`REASONABLE_KCAT_DEFAULT`, substituted whenever FRENDA's own `Km_Imputed`/`Kcat_Imputed` columns flag a value as a crude generic fallback (confirmed ~40% of fetched values were exactly this) rather than a real per-enzyme number.
- Global constants (150 uM budget, D=10 dilution, 4 h duration, enzyme-budget-split-per-pathway): top of `scripts/pathway_generator.py` (`TOTAL_BUDGET_MM`, `DILUTION_D`, `DURATION_S`). Changing `DILUTION_D`/`DURATION_S` only requires rerunning `run_pathway_sims.py`; changing enzyme topology/EC mapping/kinetics requires the same (it regenerates `model/pathways/*.txt` from scratch every run, so files there are always reproducible, not hand-edited).
- Yield metric: reported as both a 4h-endpoint value AND a peak-over-the-run value (`scripts/run_pathway_sims.py::simulate_one`) -- several designs are limited by a one-time NADH/NADPH bolus that depletes well before 4h, after which the reversible reaction runs backward and erodes the endpoint reading, so peak is the more realistic "harvested at the right time" number.

### Outputs
- `model/pathways/<Design_ID>.txt` — 30 standalone Antimony models (base model + spliced-in pathway), regenerated fresh by every `pathway_generator.py`/`run_pathway_sims.py` run.
- `results/pathway_designs_table.csv` / `results/pathway_designs.json` — the deliverable table (enzymes + expression levels, enzyme count, new substrate/cofactor additions, module classification, peak/final yield, product MSP, notes).
- `results/pathway_sim_results.json` — raw per-design simulation output (product initial/final/peak concentration, peak time, yield).

