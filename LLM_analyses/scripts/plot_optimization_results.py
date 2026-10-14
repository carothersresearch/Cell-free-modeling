"""
Dumbbell chart: baseline (D=10, original dose) vs. optimized (D, dose swept
over 10-400 / 1-100 mM) peak product titer for all 30 pathways, sorted by
fold-improvement. Annotates the optimal (D, dose) found for each.
"""
import os
import json

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = r"C:\Users\mkcoo\Documents\UW\LLM_CFS_2"
OUT_PATH = os.path.join(ROOT, "results", "pathway_analysis", "07_optimization_improvement.png")

with open(os.path.join(ROOT, "results", "pathway_optimization.json"), encoding="utf-8") as f:
    data = json.load(f)

FLOOR = 1e-5  # mM, log-scale display floor for near-zero values
for d in data:
    d["baseline_disp"] = max(d["baseline_peak_mM"], FLOOR)
    d["optimal_disp"] = max(d["optimal_peak_mM"], FLOOR)
    d["fold"] = d["optimal_disp"] / d["baseline_disp"]

data.sort(key=lambda d: d["fold"])

n = len(data)
fig, ax = plt.subplots(figsize=(13, 0.42 * n + 1.5))
ax.set_xscale("log")

y = np.arange(n)
baseline_x = [d["baseline_disp"] for d in data]
optimal_x = [d["optimal_disp"] for d in data]

for i, d in enumerate(data):
    ax.plot([d["baseline_disp"], d["optimal_disp"]], [i, i], color="#999999", linewidth=2, zorder=1,
             solid_capstyle="round")

ax.scatter(baseline_x, y, s=70, color="#b0b0b0", edgecolor="#555555", linewidth=0.8, zorder=3, label="Baseline (D=10, original dose)")
ax.scatter(optimal_x, y, s=70, color="#2E86AB", edgecolor="#164a63", linewidth=0.8, zorder=3, label="Optimized (best D, dose)")

for i, d in enumerate(data):
    fold_label = f"{d['fold']:.1f}x" if d["fold"] < 1000 else f"{d['fold']:.0f}x"
    param_label = f"(D={d['optimal_D']:.0f}, dose={d['optimal_dose_mM']:.0f}mM)"
    ax.annotate(fold_label, (300, i), xytext=(0, 0), textcoords="offset points", va="center",
                ha="left", fontsize=8.5,
                fontweight="bold" if d["fold"] >= 2 else "normal",
                color="#164a63" if d["fold"] >= 2 else "#777777")
    ax.annotate(param_label, (300, i), xytext=(38, 0), textcoords="offset points",
                va="center", ha="left", fontsize=7.5, color="#888888", style="italic")

ax.set_yticks(y)
ax.set_yticklabels([d["pid"].replace("_", " ") for d in data], fontsize=9)
ax.set_xlim(FLOOR * 0.5, 3000)
ax.set_ylim(-1, n)
ax.set_xlabel("Peak product titer (mM, log scale)")
ax.set_title("Dilution/loading optimization: baseline vs. best-found (D, substrate dose) per pathway\n"
              "(cofactor doses unchanged; grid search over D=10-400, dose=1-100 mM; sorted by fold-improvement)",
              fontsize=12, pad=12)
ax.legend(loc="lower right", fontsize=9, frameon=True)
ax.grid(axis="x", which="both", alpha=0.25)
ax.grid(axis="y", alpha=0)

plt.tight_layout()
plt.savefig(OUT_PATH, dpi=180, bbox_inches="tight", facecolor="white")
plt.close()
print(f"Wrote {OUT_PATH}")

# summary stats
folds = [d["fold"] for d in data]
print(f"\nMedian fold-improvement: {np.median(folds):.2f}x")
print(f"Pathways with >=2x improvement: {sum(f >= 2 for f in folds)}/{n}")
print(f"Pathways with ~no improvement (<1.1x): {sum(f < 1.1 for f in folds)}/{n}")
