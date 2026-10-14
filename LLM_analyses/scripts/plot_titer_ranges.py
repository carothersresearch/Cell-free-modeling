"""
Consumes results/pathway_titer_grid.json (full D x dose grid of peak product
titer, per pathway, from titer_range_study.py) and produces:

  results/pathway_titer_range_summary.csv  -- per-pathway min/median/max/dynamic-range
  results/pathway_analysis/09_titer_range_dumbbell.svg -- range plot, one row per pathway

Answers: "what's the range of achievable titers per pathway, regardless of
D/dose tuning" -- i.e. some pathways are simply capped low no matter what,
others are highly tunable across orders of magnitude.
"""
import os
import sys
import json

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(__file__))
from paper_style import MODULE_COLORS

ROOT = r"C:\Users\mkcoo\Documents\UW\LLM_CFS_2"
GRID_PATH = os.path.join(ROOT, "results", "pathway_titer_grid.json")
META_PATH = os.path.join(ROOT, "results", "pathway_analysis", "merged_pathway_analysis_data.csv")
OUT_CSV = os.path.join(ROOT, "results", "pathway_titer_range_summary.csv")
OUT_SVG = os.path.join(ROOT, "results", "pathway_analysis", "09_titer_range_dumbbell.svg")

# Regime classification thresholds (see CLAUDE.md "Titer-range study" writeup):
#   - a "viable ceiling" is a grid-max titer that is actually a meaningful amount to have made
#   - "sensitive" means the D/dose choice swings peak titer by >=10x across the grid
CEILING_THRESHOLD_MM = 0.01
FOLD_THRESHOLD = 10.0

REGIME_HIGH = "1. High & stable (near ceiling everywhere)"
REGIME_TUNABLE = "2. Wide dynamic range (D/dose choice matters)"
REGIME_LOW = "3. Low & capped (blocked regardless of D/dose)"
REGIME_ORDER = [REGIME_HIGH, REGIME_TUNABLE, REGIME_LOW]
REGIME_BAND_COLOR = {REGIME_HIGH: "#ffffff", REGIME_TUNABLE: "#f2f2f2", REGIME_LOW: "#e4e4e4"}


def classify_regime(max_mM, dynamic_range_fold):
    if max_mM < CEILING_THRESHOLD_MM:
        # never reaches a practically useful titer, however large the fold looks --
        # a 1000x jump between two near-zero numbers is not "tunable" in any useful sense
        return REGIME_LOW
    return REGIME_TUNABLE if dynamic_range_fold >= FOLD_THRESHOLD else REGIME_HIGH


def main():
    with open(GRID_PATH, encoding="utf-8") as f:
        grid_data = json.load(f)
    meta = pd.read_csv(META_PATH)
    meta_by_id = meta.set_index("design_id").to_dict("index")

    rows = []
    for entry in grid_data["results"]:
        pid = entry["pid"]
        vals = np.array([g["peak_mM"] for g in entry["grid"] if not g.get("failed") and g["peak_mM"] is not None])
        m = meta_by_id.get(pid, {})
        floor = max(vals.min(), 1e-6)
        rows.append(dict(
            pid=pid,
            module_category=m.get("module_category", "?"),
            substrate=entry["substrate"],
            product=entry["product"],
            n_enzymes=m.get("n_enzymes"),
            n_grid_points=len(vals),
            min_mM=vals.min(),
            p25_mM=np.percentile(vals, 25),
            median_mM=np.median(vals),
            p75_mM=np.percentile(vals, 75),
            max_mM=vals.max(),
            dynamic_range_fold=vals.max() / floor,
            baseline_mM=entry["best"].get("peak_mM") if False else None,  # placeholder, filled below
        ))

    df = pd.DataFrame(rows)

    # attach baseline (D=10, original dose) and optimal (best-in-grid) values from the
    # existing optimization output, for reference markers on the plot
    with open(os.path.join(ROOT, "results", "pathway_optimization.json"), encoding="utf-8") as f:
        opt = {r["pid"]: r for r in json.load(f)}
    df["baseline_mM"] = df["pid"].map(lambda p: opt[p]["baseline_peak_mM"])
    df["optimal_mM"] = df["pid"].map(lambda p: opt[p]["optimal_peak_mM"])
    df["regime"] = [classify_regime(m, f) for m, f in zip(df["max_mM"], df["dynamic_range_fold"])]

    # group by regime (in REGIME_ORDER), sorted by max_mM descending within each group
    df["_regime_rank"] = df["regime"].map({r: i for i, r in enumerate(REGIME_ORDER)})
    df = df.sort_values(["_regime_rank", "max_mM"], ascending=[True, False]).reset_index(drop=True)
    df.drop(columns="_regime_rank").to_csv(OUT_CSV, index=False)
    print(f"Wrote {OUT_CSV}")

    # --- plot: dumbbell/range chart, grouped into 3 regime bands, colored by module ---
    modules = sorted(df["module_category"].unique())
    color_map = {mod: MODULE_COLORS[mod] for mod in modules}

    GROUP_GAP = 1.5  # extra vertical spacing (in row units) inserted between regime groups

    y_positions = []
    band_ranges = {}   # regime -> (y_min, y_max) for axhspan
    header_y = {}       # regime -> y position for the section header text
    cur_y = 0.0
    prev_regime = None
    for _, row in df.iterrows():
        if prev_regime is not None and row["regime"] != prev_regime:
            cur_y += GROUP_GAP
        if row["regime"] not in band_ranges:
            band_ranges[row["regime"]] = [cur_y - 0.6, cur_y]
            # inverted y-axis: smaller y renders higher on the page, so the header sits
            # just ABOVE the group's first row (nearer the band's top edge)
            header_y[row["regime"]] = cur_y - 0.45
        band_ranges[row["regime"]][1] = cur_y + 0.6
        y_positions.append(cur_y)
        cur_y += 1.0
        prev_regime = row["regime"]
    df["_y"] = y_positions

    fig_h = 0.34 * (cur_y + 2) + 1.9
    fig, ax = plt.subplots(figsize=(11.5, fig_h))

    floor_plot = 1e-4  # mM, plotting floor for log scale / near-zero pathways

    for regime, (y0, y1) in band_ranges.items():
        ax.axhspan(y0, y1, color=REGIME_BAND_COLOR[regime], zorder=0)

    for _, row in df.iterrows():
        i = row["_y"]
        lo = max(row["min_mM"], floor_plot)
        hi = max(row["max_mM"], floor_plot)
        c = color_map[row["module_category"]]
        ax.plot([lo, hi], [i, i], color=c, lw=2.2, alpha=0.55, solid_capstyle="round", zorder=2)
        ax.scatter(max(row["min_mM"], floor_plot), i, color=c, s=22, zorder=3, edgecolor="white", linewidth=0.4)
        ax.scatter(max(row["max_mM"], floor_plot), i, color=c, s=42, zorder=3, marker="D", edgecolor="white", linewidth=0.4)
        ax.scatter(max(row["baseline_mM"], floor_plot), i, color="black", s=26, zorder=4, marker="|")

    n_by_regime = df["regime"].value_counts()
    for regime in REGIME_ORDER:
        if regime not in header_y:
            continue
        ax.text(floor_plot * 0.55, header_y[regime], f"{regime}  (n={n_by_regime[regime]})",
                fontsize=9.5, fontweight="bold", color="#333333", va="bottom", ha="left", zorder=5)

    ax.set_yticks(y_positions)
    ax.set_yticklabels(df["pid"], fontsize=8)
    ax.set_xscale("log")
    ax.set_xlabel("Peak product titer across full D x dose grid (mM, log scale)")
    ax.set_xlim(floor_plot * 0.5, 30)
    ax.set_ylim(-0.8, cur_y - 1.0 + 0.8)
    ax.invert_yaxis()
    ax.axvline(floor_plot, color="#cccccc", lw=0.8, zorder=1)
    ax.axvline(CEILING_THRESHOLD_MM, color="#999999", lw=0.9, ls="--", zorder=1)
    # label the threshold line in the empty gap between the 1st and 2nd regime bands
    # (guaranteed clear of both row content and the left-anchored section headers)
    gap0_y = (band_ranges[REGIME_HIGH][1] + band_ranges[REGIME_TUNABLE][0]) / 2
    ax.text(CEILING_THRESHOLD_MM, gap0_y, "viable-ceiling threshold (0.01 mM)",
            fontsize=7, color="#777777", va="center", ha="center",
            bbox=dict(facecolor="white", edgecolor="none", pad=1.5))
    ax.set_title("Range of achievable titers per pathway (49-point D x dose grid, D: 10-400, dose: 1-100 mM)\n"
                  "circle = grid min, diamond = grid max, black tick = original baseline (D=10, design dose)\n"
                  "grouped into 3 regimes: ceiling >=0.01 mM & fold <10x = stable-high; ceiling >=0.01 mM & fold >=10x = tunable; ceiling <0.01 mM = capped",
                  fontsize=9.5)

    handles = [plt.Line2D([0], [0], color=color_map[m], lw=3, label=m) for m in modules]
    ax.legend(handles=handles, loc="lower right", fontsize=7.5, title="module", title_fontsize=8, framealpha=0.95)

    ax.grid(axis="x", which="major", color="#dddddd", zorder=1)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig.tight_layout()
    fig.savefig(OUT_SVG)
    print(f"Wrote {OUT_SVG}")

    # --- quick console summary ---
    print("\nTop 5 highest achievable ceiling (max_mM):")
    print(df.sort_values("max_mM", ascending=False)[["pid", "min_mM", "max_mM", "dynamic_range_fold"]].head(5).to_string(index=False))
    print("\nBottom 5 lowest ceiling (max_mM) -- capped regardless of D/dose:")
    print(df.sort_values("max_mM", ascending=True)[["pid", "min_mM", "max_mM", "dynamic_range_fold"]].head(5).to_string(index=False))
    print("\nMost D/dose-sensitive (largest dynamic_range_fold):")
    print(df.sort_values("dynamic_range_fold", ascending=False)[["pid", "min_mM", "max_mM", "dynamic_range_fold"]].head(8).to_string(index=False))
    print("\nLeast sensitive (smallest dynamic_range_fold, i.e. flat regardless of D/dose):")
    print(df.sort_values("dynamic_range_fold", ascending=True)[["pid", "min_mM", "max_mM", "dynamic_range_fold"]].head(8).to_string(index=False))


if __name__ == "__main__":
    main()
