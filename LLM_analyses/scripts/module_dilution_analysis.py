"""
Why don't biosynthetic modules dilute linearly?

Classifies all 54 native reactions (R1-R54; R55-R57 are the heterologous
engineered pathway, excluded here) into biosynthetic modules, sweeps
dilution factor D from 10 to 400 for a PURE native-lysate condition (no
heterologous enzymes, no exogenous cofactor/substrate spikes -- i.e. every
species scales with dilution exactly as encoded in the model), and measures
each module's instantaneous flux (via roadrunner's getReactionRates(), no
time integration needed -- this isolates the pure enzyme-kinetics dilution
effect from any downstream depletion/accumulation dynamics) as a function
of D.

For each reaction, also parses its Antimony rate-law block directly to find
every Km_<Substrate>_<Enzyme> parameter it references, and computes that
substrate's concentration/Km ratio at the native (D=10, undiluted) reference
condition. This tests the mechanistic hypothesis: reactions running well
below Km for their substrate (S/Km << 1) should show closer-to-QUADRATIC
flux suppression under dilution (both enzyme E and substrate S dilute
together: v ~ kcat*E*S/Km => v ~ 1/D^2), while reactions running near/above
Km (S/Km >~ 1, e.g. because the substrate pool is naturally very large, like
water/bicarbonate/ammonia) should show closer-to-LINEAR suppression
(v ~ kcat*E => v ~ 1/D), since they're not simultaneously substrate-starved.

Outputs:
  results/analysis/module_dilution_scan.csv       (raw per-reaction, per-D rates)
  results/analysis/module_km_analysis.csv         (per-reaction S/Km + fitted alpha)
  results/analysis/10_module_flux_share_vs_dilution.png
  results/analysis/11_module_dilution_exponent.png
  results/analysis/12_alpha_vs_substrate_saturation.png
  results/analysis/module_dilution_summary.txt

Usage:
    ./.venv_sim/Scripts/python.exe scripts/module_dilution_analysis.py
"""
import sys
import os
import re
import json
from pathlib import Path
from collections import defaultdict

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from sim2 import Design, _fresh_runner, SPECIES_MAP, FITTED_PARAMS, MODEL_PATH

ROOT = Path(__file__).resolve().parent.parent
OUT_DIR = ROOT / "results" / "analysis"
OUT_DIR.mkdir(exist_ok=True)

sns.set_theme(style="white", font_scale=1.0)
OKABE_ITO = ["#E69F00", "#56B4E9", "#009E73", "#D55E00", "#CC79A7", "#0072B2", "#F0E442", "#999999",
             "#8B4513", "#4B0082", "#2E8B57", "#B22222"]

# ---------------------------------------------------------------------------
# Reaction -> biosynthetic module classification (from enzyme EC identities
# and reaction stoichiometry read directly off the model file)
# ---------------------------------------------------------------------------
MODULE_OF_REACTION = {
    # TCA cycle
    "R1": "TCA cycle", "R7": "TCA cycle", "R8": "TCA cycle", "R20": "TCA cycle",
    "R25": "TCA cycle", "R27": "TCA cycle", "R28": "TCA cycle", "R41": "TCA cycle",
    "R42": "TCA cycle", "R46": "TCA cycle",
    # Glycolysis / gluconeogenesis
    "R10": "Glycolysis/Gluconeogenesis", "R13": "Glycolysis/Gluconeogenesis",
    "R18": "Glycolysis/Gluconeogenesis", "R19": "Glycolysis/Gluconeogenesis",
    "R21": "Glycolysis/Gluconeogenesis", "R22": "Glycolysis/Gluconeogenesis",
    "R33": "Glycolysis/Gluconeogenesis", "R34": "Glycolysis/Gluconeogenesis",
    "R35": "Glycolysis/Gluconeogenesis", "R45": "Glycolysis/Gluconeogenesis",
    "R32": "Glycolysis/Gluconeogenesis",
    # Pyruvate metabolism / anaplerosis
    "R36": "Pyruvate metabolism", "R53": "Pyruvate metabolism", "R52": "Pyruvate metabolism",
    "R37": "Pyruvate metabolism", "R38": "Pyruvate metabolism",
    "R30": "Pyruvate metabolism", "R31": "Pyruvate metabolism",
    # Pentose phosphate pathway
    "R39": "Pentose phosphate pathway", "R40": "Pentose phosphate pathway",
    "R43": "Pentose phosphate pathway", "R44": "Pentose phosphate pathway",
    "R54": "Pentose phosphate pathway",
    # Fatty acid synthesis
    "R2": "Fatty acid synthesis", "R3": "Fatty acid synthesis", "R4": "Fatty acid synthesis",
    "R5": "Fatty acid synthesis", "R6": "Fatty acid synthesis", "R11": "Fatty acid synthesis",
    "R12": "Fatty acid synthesis", "R14": "Fatty acid synthesis", "R26": "Fatty acid synthesis",
    # Fatty acid beta-oxidation / activation
    "R15": "Fatty acid oxidation", "R16": "Fatty acid oxidation", "R17": "Fatty acid oxidation",
    "R24": "Fatty acid oxidation",
    # Glyoxylate shunt
    "R47": "Glyoxylate shunt",
    # Acetate metabolism
    "R48": "Acetate metabolism", "R50": "Acetate metabolism",
    # One-carbon / serine metabolism
    "R49": "One-carbon/serine metabolism",
    # Nucleotide / energy metabolism
    "R9": "Nucleotide/energy metabolism", "R29": "Nucleotide/energy metabolism",
    # CO2/bicarbonate interconversion
    "R51": "CO2/bicarbonate interconversion",
    # Glycogen/starch metabolism
    "R23": "Glycogen/starch metabolism",
}
NATIVE_REACTIONS = list(MODULE_OF_REACTION.keys())

# ---------------------------------------------------------------------------
# Parse each reaction's rate-law block from the raw Antimony text to find
# every Km_<Substrate>_<Enzyme> parameter it references
# ---------------------------------------------------------------------------
with open(MODEL_PATH) as f:
    model_text = f.read()

# split into per-reaction blocks: text from "R{n} :" up to (not including) the next "R{n+1} :"
block_pattern = re.compile(r"^(R\d+) :.*?(?=^R\d+ :|\Z)", re.MULTILINE | re.DOTALL)
reaction_blocks = {m.group(1): m.group(0) for m in block_pattern.finditer(model_text)}

km_pattern = re.compile(r"Km_(\w+?)_(eEC\d+|hEC\d+)\b")

def parse_km_terms(rxn_id):
    block = reaction_blocks.get(rxn_id, "")
    matches = km_pattern.findall(block)
    # dedupe while preserving order
    seen = []
    for metab, enz in matches:
        key = (metab, f"Km_{metab}_{enz}")
        if key not in seen:
            seen.append(key)
    return seen  # list of (metabolite_name, Km_param_name)

# ---------------------------------------------------------------------------
# Native (D=10, no exogenous spikes, no heterologous) reference concentrations
# ---------------------------------------------------------------------------
def build_native_runner(D):
    design = Design(
        dilution_factor_D=float(D),
        het_enzymes_uM={"hEC6411": 0.0, "hEC11719": 0.0, "hEC43117": 0.0},
        eEC11137_overexpr_uM=0.0,
        exogenous_overrides_mM={k: 0.0 for k in
                                 ["ATP", "AcetylCoA", "HCO3", "NADH", "NADPH", "Formate", "L_Serine"]},
    )
    r = _fresh_runner()
    dil = 1.0 / design.dilution_factor_D
    r.model["dilution_factor"] = dil
    for sp, info in SPECIES_MAP.items():
        if info["driver"].startswith("p_eEC") and sp != "eEC11137":
            r.model[info["driver"]] = 1.0
    r.model["p_eEC11137"] = 1.0
    for hid in ["hEC6411", "hEC11719", "hEC43117"]:
        r.model[f"p_{hid}"] = 0.0
    for sp, info in SPECIES_MAP.items():
        if info["driver"].startswith("v"):
            v_default = FITTED_PARAMS.get(info["driver"])
            if v_default is not None:
                r.model[info["driver"]] = v_default
    r.model["v73"] = 10  # inert Hydroxycitrate placeholder (see CLAUDE.md)
    r.reset()
    return r


r_native = build_native_runner(10.0)
native_conc = {sid: r_native.model[sid] for sid in r_native.model.getFloatingSpeciesIds()}

km_analysis_rows = []
for rxn_id in NATIVE_REACTIONS:
    terms = parse_km_terms(rxn_id)
    for metab, km_name in terms:
        km_val = FITTED_PARAMS.get(km_name)
        conc = native_conc.get(metab)
        if km_val is None or conc is None:
            continue
        km_analysis_rows.append(dict(
            reaction=rxn_id, module=MODULE_OF_REACTION[rxn_id],
            substrate=metab, km_param=km_name, km_value=km_val,
            native_conc_mM=conc, S_over_Km=conc / km_val,
        ))
km_df = pd.DataFrame(km_analysis_rows)

# per-reaction "most-limiting" S/Km = minimum across its substrate terms
limiting = km_df.groupby("reaction")["S_over_Km"].min().rename("min_S_over_Km")

# ---------------------------------------------------------------------------
# Dilution sweep: instantaneous reaction rates at t=0 (no integration)
# ---------------------------------------------------------------------------
D_VALUES = np.geomspace(10, 400, 15)
rate_rows = []
for D in D_VALUES:
    r = build_native_runner(D)
    rate_ids = list(r.model.getReactionIds())
    rates = r.model.getReactionRates()
    rate_map = dict(zip(rate_ids, rates))
    for rxn_id in NATIVE_REACTIONS:
        rate_rows.append(dict(dilution_factor_D=D, reaction=rxn_id,
                               module=MODULE_OF_REACTION[rxn_id],
                               rate=rate_map.get(rxn_id, np.nan)))
    print(f"D={D:7.2f}  done")

rate_df = pd.DataFrame(rate_rows)
rate_df["abs_rate"] = rate_df["rate"].abs()
rate_df.to_csv(OUT_DIR / "module_dilution_scan.csv", index=False)

# ---------------------------------------------------------------------------
# Aggregate by module: total |flux| per module per D, and module share of
# total native flux
# ---------------------------------------------------------------------------
module_flux = rate_df.groupby(["dilution_factor_D", "module"])["abs_rate"].sum().reset_index()
total_flux = module_flux.groupby("dilution_factor_D")["abs_rate"].sum().rename("total")
module_flux = module_flux.merge(total_flux, on="dilution_factor_D")
module_flux["share"] = module_flux["abs_rate"] / module_flux["total"]

# fit log-log slope (dilution exponent alpha, flux ~ D^-alpha) per module
alphas = {}
for mod, sub in module_flux.groupby("module"):
    sub = sub[sub["abs_rate"] > 1e-12]
    if len(sub) < 5:
        continue
    slope, intercept = np.polyfit(np.log(sub["dilution_factor_D"]), np.log(sub["abs_rate"]), 1)
    alphas[mod] = -slope
alpha_s = pd.Series(alphas).sort_values()

# merge with reaction-level S/Km (module-level: min over member reactions' min_S_over_Km)
module_min_skm = {}
for rxn_id, mod in MODULE_OF_REACTION.items():
    if rxn_id in limiting.index:
        module_min_skm.setdefault(mod, []).append(limiting[rxn_id])
module_min_skm = {m: min(v) for m, v in module_min_skm.items()}

summary_lines = []
def log(msg=""):
    print(msg)
    summary_lines.append(msg)

log("=" * 78)
log("Module dilution-scaling analysis (native lysate only, D=10 to D=400)")
log("=" * 78)
log("\nDilution exponent alpha per module (flux ~ D^-alpha; alpha=1 expected if")
log("only enzyme dilutes; alpha=2 if BOTH enzyme and a limiting substrate dilute")
log("together; alpha<1 = 'protected', gains relative share as dilution increases):\n")
for mod, a in alpha_s.items():
    skm = module_min_skm.get(mod, np.nan)
    log(f"  {mod:32s} alpha={a:.2f}   most-limiting native S/Km = {skm:.3g}" if not np.isnan(skm)
        else f"  {mod:32s} alpha={a:.2f}   (no Km data)")

log("\nModule share of total native flux, low dilution (D=10) vs high dilution (D=400):")
share_lo = module_flux[module_flux["dilution_factor_D"] == module_flux["dilution_factor_D"].min()].set_index("module")["share"]
share_hi = module_flux[module_flux["dilution_factor_D"] == module_flux["dilution_factor_D"].max()].set_index("module")["share"]
share_cmp = pd.DataFrame({"share_D10": share_lo, "share_D400": share_hi})
share_cmp["fold_change"] = share_cmp["share_D400"] / share_cmp["share_D10"]
share_cmp = share_cmp.sort_values("fold_change", ascending=False)
for mod, row in share_cmp.iterrows():
    log(f"  {mod:32s} {row['share_D10']*100:6.2f}% -> {row['share_D400']*100:6.2f}%  "
        f"({row['fold_change']:.2f}x relative share)")

# ---------------------------------------------------------------------------
# Plot 10: module flux share vs dilution (stream/line plot)
# ---------------------------------------------------------------------------
modules_sorted = share_cmp.index.tolist()
fig, ax = plt.subplots(figsize=(10, 7))
for i, mod in enumerate(modules_sorted):
    sub = module_flux[module_flux["module"] == mod].sort_values("dilution_factor_D")
    ax.plot(sub["dilution_factor_D"], sub["share"] * 100, "-o", color=OKABE_ITO[i % len(OKABE_ITO)],
            markersize=4, linewidth=1.8, label=mod)
ax.set_xscale("log")
ax.set_yscale("log")
ax.set_xlabel("Dilution factor D (log scale)")
ax.set_ylabel("Share of total native flux (%, log scale)")
ax.set_title("Which modules gain relative share of total flux as dilution increases?\n"
             "(log y-axis: Glycolysis dominates in absolute terms, but relative TRENDS for\n"
             "every other module are only visible on a log scale)", fontsize=12)
ax.legend(frameon=False, fontsize=8, loc="center left", bbox_to_anchor=(1.01, 0.5))
sns.despine(ax=ax)
plt.tight_layout()
plt.savefig(OUT_DIR / "10_module_flux_share_vs_dilution.png", dpi=180, bbox_inches="tight")
plt.close()

# ---------------------------------------------------------------------------
# Plot 11: dilution exponent alpha per module (bar chart)
# ---------------------------------------------------------------------------
fig, ax = plt.subplots(figsize=(9, 6.5))
colors = [OKABE_ITO[2] if a < 0.9 else (OKABE_ITO[5] if a < 1.5 else OKABE_ITO[3]) for a in alpha_s.values]
ax.barh(range(len(alpha_s)), alpha_s.values, color=colors)
ax.set_yticks(range(len(alpha_s)))
ax.set_yticklabels(alpha_s.index)
ax.set_ylim(-0.7, len(alpha_s) + 0.9)
ax.axvline(1.0, color="black", linewidth=1, linestyle="--", alpha=0.6)
ax.axvline(2.0, color="black", linewidth=1, linestyle=":", alpha=0.4)
ax.text(1.02, len(alpha_s) + 0.55, "linear\n(enzyme only)", fontsize=8, color="#555555", va="top")
ax.text(2.02, len(alpha_s) + 0.55, "quadratic\n(enzyme + substrate)", fontsize=8, color="#555555", va="top")
for i, v in enumerate(alpha_s.values):
    ax.text(v + 0.05, i, f"{v:.2f}", va="center", fontsize=8)
ax.set_xlim(0, 3.6)
ax.set_xlabel("Fitted dilution exponent alpha  (flux ~ D^-alpha)")
ax.set_title("Dilution sensitivity by biosynthetic module\n(lower = more 'protected' / relatively more active when dilute)", fontsize=12, pad=14)
sns.despine(ax=ax)
plt.tight_layout()
plt.savefig(OUT_DIR / "11_module_dilution_exponent.png", dpi=180, bbox_inches="tight")
plt.close()

# ---------------------------------------------------------------------------
# Plot 12: alpha vs native substrate saturation (mechanistic validation)
# ---------------------------------------------------------------------------
plot_data = pd.DataFrame({"alpha": alpha_s, "min_S_over_Km": pd.Series(module_min_skm)}).dropna()
fig, ax = plt.subplots(figsize=(8.5, 7))
for i, (mod, row) in enumerate(plot_data.iterrows()):
    ax.scatter(row["min_S_over_Km"], row["alpha"], s=110, color=OKABE_ITO[i % len(OKABE_ITO)], zorder=3,
               edgecolor="white", linewidth=0.8)
    ax.annotate(mod, (row["min_S_over_Km"], row["alpha"]), fontsize=8, xytext=(6, 4), textcoords="offset points")
ax.set_xscale("log")
ax.axhline(1.0, color="gray", linewidth=1, linestyle="--", alpha=0.5)
ax.axhline(2.0, color="gray", linewidth=1, linestyle=":", alpha=0.4)
ax.axvline(1.0, color="gray", linewidth=1, linestyle="-", alpha=0.3)
ax.set_xlabel("Most-limiting substrate's native S/Km (log scale)\n(<1 = substrate-limited, >1 = substrate-saturated)")
ax.set_ylabel("Fitted dilution exponent alpha")
ax.set_title("Mechanistic check: substrate saturation predicts dilution sensitivity", fontsize=12)
sns.despine(ax=ax)
plt.tight_layout()
plt.savefig(OUT_DIR / "12_alpha_vs_substrate_saturation.png", dpi=180, bbox_inches="tight")
plt.close()

km_df.to_csv(OUT_DIR / "module_km_analysis.csv", index=False)
with open(OUT_DIR / "module_dilution_summary.txt", "w") as f:
    f.write("\n".join(summary_lines))

print(f"\nWrote outputs to {OUT_DIR}")
