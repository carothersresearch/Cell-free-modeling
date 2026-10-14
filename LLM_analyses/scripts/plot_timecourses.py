"""
Time-course comparison for all 30 pathways: product concentration vs. time,
un-optimized (D=10, original dose) overlaid with optimized (best-found D,
dose) on the same axes. Substrates ordered by ascending feedstock cost
(matching 06_pathway_map.png); within each substrate's row, the 3 pathways
ordered left-to-right by descending optimized peak yield.

Split into two 5 (substrate) x 3 (pathway) grids (15 panels each, covering
the cheapest/priciest 5 substrates respectively) so each half fits on one
A4 page with 1" margins on all sides (6.27" x 9.69" usable area) at a
legible 12pt Arial throughout -- the original single 10x3 grid was too
small to read once printed/placed in a document.
"""
import os
import re
import sys
import json

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(__file__))
from loader import loada
from pathway_generator import DURATION_S, NATIVE_OVERLAP_SUBSTRATES, ORTHOGONAL_SUBSTRATES, BASE_HET_PARAMS

ROOT = r"C:\Users\mkcoo\Documents\UW\LLM_CFS_2"
OUT_DIR = os.path.join(ROOT, "results", "pathway_analysis")
N_POINTS = 200

# A4 (8.27" x 11.69") minus 1" margins on all sides
A4_USABLE_W = 8.27 - 2.0
A4_USABLE_H = 11.69 - 2.0

plt.rcParams["font.family"] = "Arial"
plt.rcParams["font.size"] = 12
FIG_FONT = 12

# shorter panel titles for products whose raw names are too long/underscored
# to read as a panel title in the ~1.7"-wide column this layout allows
PRODUCT_RENAME = {
    "DGlucose6phosphate": "G6P", "PHB_monomer": "PHB monomer",
    "L_Alanine": "L-Alanine", "L_Serine": "L-Serine",
}

with open(os.path.join(ROOT, "results", "pathway_sim_results_unoptimized.json"), encoding="utf-8") as f:
    baseline_results = {r["pid"]: r for r in json.load(f) if "meta" in r}
with open(os.path.join(ROOT, "results", "pathway_optimization.json"), encoding="utf-8") as f:
    opt_results = {r["pid"]: r for r in json.load(f)}

subs_cost = pd.read_csv(os.path.join(ROOT, "results", "substrate_candidates.csv")).dropna(subset=["substrate"])


def parse_cost_midpoint(s):
    nums = re.findall(r"\d+\.?\d*", str(s).split("(")[0].split(";")[0])
    nums = [float(n) for n in nums[:2]]
    return sum(nums) / len(nums) if nums else np.nan


subs_cost = subs_cost.set_index("substrate")
subs_cost["cost_mid"] = subs_cost["est_MSP_usd_per_kg"].map(parse_cost_midpoint)
SUBSTRATE_RENAME = {"HCO3": "Bicarbonate/CO2", "S_Lactate": "S-Lactate", "D_Xylose": "D-Xylose"}


def substrate_cost(key):
    return subs_cost.loc[SUBSTRATE_RENAME.get(key, key), "cost_mid"]


def run_timecourse(meta, D, dose):
    r = loada(meta["file"])
    for pgene in BASE_HET_PARAMS:
        setattr(r, pgene, 0.0)
    r.dilution_factor = 1.0 / D
    r.reset()

    substrate = meta["substrate"]
    if substrate in NATIVE_OVERLAP_SUBSTRATES:
        r[substrate] = float(r[substrate]) + dose
    elif substrate in ORTHOGONAL_SUBSTRATES:
        r[substrate] = dose

    for cofactor, cdose in meta["cofactor_additions"].items():
        r[cofactor] = float(r[cofactor]) + cdose

    product = meta["product"]
    result = r.simulate(0, DURATION_S, N_POINTS)
    colnames = [c.strip("[]") for c in result.colnames]
    idx = colnames.index(product)
    t_h = [result[j][0] / 3600.0 for j in range(N_POINTS)]
    y = [float(result[j][idx]) for j in range(N_POINTS)]
    return t_h, y


# ---------------------------------------------------------------- gather all timecourses
rows = []
for pid, base in baseline_results.items():
    meta = base["meta"]
    opt = opt_results[pid]
    t_base, y_base = run_timecourse(meta, 10.0, meta["substrate_dose_mM"])
    t_opt, y_opt = run_timecourse(meta, opt["optimal_D"], opt["optimal_dose_mM"])
    rows.append(dict(
        pid=pid, substrate_key=meta["substrate"], product=meta["product"],
        t_base=t_base, y_base=y_base, t_opt=t_opt, y_opt=y_opt,
        opt_yield=opt["optimal_yield_pct_peak"] if opt["optimal_yield_pct_peak"] is not None else -1,
        cost=substrate_cost(meta["substrate"]),
    ))
    print(f"{pid}: done")

df = pd.DataFrame(rows)

# ---------------------------------------------------------------- layout ordering
substrate_order = (df[["substrate_key", "cost"]].drop_duplicates()
                   .sort_values("cost")["substrate_key"].tolist())
halves = [substrate_order[:5], substrate_order[5:]]


def make_grid(substrates, part_label, out_path):
    n_rows = len(substrates)
    fig, axes = plt.subplots(n_rows, 3, figsize=(A4_USABLE_W, A4_USABLE_H), sharex=True)

    for i, sub in enumerate(substrates):
        sub_rows = df[df["substrate_key"] == sub].sort_values("opt_yield", ascending=False).reset_index(drop=True)
        for j in range(3):
            ax = axes[i, j]
            row = sub_rows.iloc[j]
            ax.plot(row["t_base"], row["y_base"], color="#a0a0a0", linewidth=1.6, linestyle="--",
                    label="Un-optimized" if (i == 0 and j == 0) else None)
            ax.plot(row["t_opt"], row["y_opt"], color="#D6604D", linewidth=1.9,
                    label="Optimized" if (i == 0 and j == 0) else None)
            title = PRODUCT_RENAME.get(row["product"], row["product"]).replace("_", " ")
            ax.set_title(title, fontsize=FIG_FONT, pad=3)
            ax.tick_params(labelsize=FIG_FONT)
            ax.xaxis.set_major_locator(plt.MaxNLocator(3))
            ax.yaxis.set_major_locator(plt.MaxNLocator(3))
            # near-zero panels (e.g. BDO23 ~1e-5 to 1e-8 mM) would otherwise get a
            # separate scientific-notation "offset text" floating just above the
            # axes, colliding with the title directly above it -- an explicit
            # per-tick formatter avoids that extra element entirely
            maxabs = max([abs(v) for v in row["y_base"]] + [abs(v) for v in row["y_opt"]], default=0)
            if 0 < maxabs < 1e-3:
                ax.yaxis.set_major_formatter(plt.FuncFormatter(lambda v, pos: f"{v:.1e}"))
            ax.grid(alpha=0.25, linewidth=0.5)
            for spine in ("top", "right"):
                ax.spines[spine].set_visible(False)
            ax.margins(y=0.15)
            if j == 0:
                ax.set_ylabel(sub.replace("_", " "), fontsize=FIG_FONT, fontweight="bold")
            if i == n_rows - 1:
                ax.set_xlabel("Time (h)", fontsize=FIG_FONT)

    handles, labels = axes[0, 0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="upper center", bbox_to_anchor=(0.5, 1.0), ncol=2,
               fontsize=FIG_FONT, frameon=False)
    fig.text(0.005, 0.5, "Product concentration (mM)", va="center", rotation="vertical", fontsize=FIG_FONT)

    plt.tight_layout(rect=[0.02, 0, 1, 0.965])
    plt.savefig(out_path)
    plt.close()
    print(f"Wrote {out_path}")


make_grid(halves[0], "Part 1 of 2 (lowest-cost 5 substrates)",
          os.path.join(OUT_DIR, "08_timecourses_before_after_part1.svg"))
make_grid(halves[1], "Part 2 of 2 (highest-cost 5 substrates)",
          os.path.join(OUT_DIR, "08_timecourses_before_after_part2.svg"))
