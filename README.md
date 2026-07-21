# PaperFigures

Code to reproduce the figures in the manuscript, plus two supporting pieces of work: a parameter-set filtering pipeline that produces the kinetic parameter ensemble used throughout the figures, and a synthetic-data experiment used to validate the fitting/inference approach.

All notebooks in this folder are self-contained — they load the model and support files (`iMC057.txt`, `Reactions.csv`, `SpeciesBaseMechanisms.csv`, `labels.csv`, `filtering_parameter_sets/final_params.pkl`) via **relative paths**, so they must be run with this folder (`PaperFigures/`) as the working directory.

## Contents

```
PaperFigures/
├── Figure_2.ipynb                  # KEGG pathway map layout of the model
├── Figure_3.ipynb                  # empty — not yet populated
├── Figure_4.ipynb                  # flux/module analysis, heatmaps, sinks
├── Figure_5.ipynb                  # malate production, titrations, dilution, sensitivity
├── Figure_6.ipynb                  # candidate engineering strategy comparison
├── iMC057.txt                      # base kinetic model (Antimony)
├── iMC057_plus_ATP_regen.txt       # iMC057 variant with an added ATP regeneration system
├── Reactions.csv                   # reaction/kinetic parameter table backing the model
├── SpeciesBaseMechanisms.csv       # species/enzyme metadata backing the model
├── labels.csv                      # KEGG compound ID -> display label mapping (used by Figure_2)
├── ko01100.xml                     # KEGG global metabolic pathway map (layout source for Figure_2)
├── filtering_parameter_sets/       # kinetic parameter-set filtering pipeline
└── synthetic_data_experiment/      # synthetic-data generation/validation experiment
```

## Figures

- **Figure_2** — Maps the model's reactions onto KEGG pathway map (`ko01100.xml`) coordinates to lay out the reaction network diagram.
- **Figure_3** — Placeholder; notebook is currently empty.
- **Figure_4** — 4A: flux through biosynthetic modules over time/dilution factor. 4B: per-module flux zoom-ins. 4C: heatmap of reaction/metabolite changes across conditions. Also includes a parameter-set clustering sidebar. 4D: sink/output analysis.
- **Figure_5** — 5B: malate production across experimental conditions. 5C: enzyme titrations. 5D: substrate titrations. 5E: effect of dilution on malate/CO2 production. 5F: relative malate contribution from TCA vs. rTCA across dilution factors. 5G: MCA-inspired sensitivity coefficient analysis.
- **Figure_6** — Simulates the base pathway plus ten candidate engineering strategies (e.g. increasing pyruvate/serine/bicarbonate, omitting Fdh or HCT, reallocating protein budget, overexpressing MDH, adding ATP regeneration) and plots the comparison.

## `filtering_parameter_sets/`

`filter.ipynb` runs a series of consistency checks ("dummy tests") against sampled kinetic parameter sets — e.g. HCT increasing citrate, NADH consumption, glycine production from serine, NADH regeneration, and overall pathway efficacy — and filters down to the parameter sets that pass. Outputs:

- `all_labels_parameters_error.dat` / `top3_labels_parameters_error.dat` — per-parameter-set error/labels from the filtering runs.
- `final_params.pkl` — the filtered parameter set ensemble consumed by Figures 4–6.
- `final_params.tar` / `final_params.7z` — compressed archives of the same parameter set (for storage/transfer).

## `synthetic_data_experiment/`

A validation experiment, separate from the main figures, that generates synthetic datasets from the kinetic model:

- **make_perfect_dataset.ipynb** — simulates the base model (`iMC057.txt`) from a grid of starting conditions, including specific enzyme knockout/overexpression variants (e.g. `200X_mdh_neg`, `200X_pyc`), to produce noiseless ("perfect") ground-truth datasets.
- **make_diff_datasets.ipynb** — perturbs kinetic parameters/enzyme/metabolite priors and regenerates an Antimony model from those priors, then simulates high/medium/low-quantity experiment sets with high/medium/low observation noise added, to test how dataset size and quality affect downstream fitting.
- `perfect_datasets.zip` / `synthetic_datasets.zip` — generated output datasets from the two notebooks above.

## Setup

From the repository root:

```sh
conda env create -f environment.yml
conda activate <env name from environment.yml>
```

Then open notebooks in this folder directly (not from `notebooks/` or another working directory), since all file references are relative to `PaperFigures/`.
