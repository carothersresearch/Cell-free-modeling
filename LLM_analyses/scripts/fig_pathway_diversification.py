"""
Combined publication figure for the heterologous pathway diversification study,
styled to match CellFree_Fig4_v3.png (see paper_style.py for the sampled
palette/rcParams -- same MODULE_COLORS mapping used in both panels here so line
color means the same thing across the whole figure):

  Panel A: substrate -> product map for all 30 pathway designs (same geometry as
           scripts/plot_pathway_map.py).
  Panel B: range of achievable titer per pathway across the full D x dose grid,
           grouped into 3 regimes (same geometry as scripts/plot_titer_ranges.py).

Usage:
    ./.venv_sim/Scripts/python.exe scripts/fig_pathway_diversification.py
"""
import os
import re
import sys
import json

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

sys.path.insert(0, os.path.dirname(__file__))
from paper_style import MODULE_ORDER, MODULE_COLORS, FIG4_GRAY, add_panel_label

ROOT = r"C:\Users\mkcoo\Documents\UW\LLM_CFS_2"
OUT_PATH = os.path.join(ROOT, "results", "pathway_analysis", "Fig_pathway_diversification.png")


# =============================================================================
# Panel A: substrate -> product map (scripts/plot_pathway_map.py, recolored)
# =============================================================================
def panel_a(fig, ax):
    df = pd.read_csv(os.path.join(ROOT, "results", "pathway_designs_table.csv"))
    subs_cost = pd.read_csv(os.path.join(ROOT, "results", "substrate_candidates.csv"))
    subs_cost = subs_cost.dropna(subset=["substrate"])

    SUBSTRATE_RENAME = {"HCO3": "Bicarbonate/CO2", "S_Lactate": "S-Lactate", "D_Xylose": "D-Xylose"}
    df["substrate_key"] = df["substrate"].map(lambda s: SUBSTRATE_RENAME.get(s, s))

    def parse_cost_midpoint(s):
        nums = re.findall(r"\d+\.?\d*", str(s).split("(")[0].split(";")[0])
        nums = [float(n) for n in nums[:2]]
        return sum(nums) / len(nums) if nums else np.nan

    subs_cost = subs_cost.set_index("substrate")
    subs_cost["cost_mid"] = subs_cost["est_MSP_usd_per_kg"].map(parse_cost_midpoint)
    df["substrate_cost"] = df["substrate_key"].map(subs_cost["cost_mid"])

    df["product_clean"] = df["product"].str.replace(r"\s*\((native|novel)\)", "", regex=True)
    df["product_is_native"] = df["product"].str.contains(r"\(native\)")

    def categorize_module(m):
        m = m.lower()
        if "one-carbon" in m or "serine" in m:
            return "C1/serine metabolism"
        if "tca" in m or "anaplerosis" in m:
            return "TCA/anaplerosis"
        if "fatty acid" in m or "ketone" in m:
            return "Fatty acid/ketone metabolism"
        if "amino acid" in m:
            return "Amino acid biosynthesis"
        if "pentose phosphate" in m:
            return "Pentose phosphate pathway"
        if "glycolysis" in m:
            return "Glycolysis"
        if "pyruvate" in m:
            return "Pyruvate metabolism"
        return "Other"

    df["module_category"] = df["module"].map(categorize_module)

    substrate_order = (df[["substrate_key", "substrate_cost"]].drop_duplicates()
                       .sort_values("substrate_cost").reset_index(drop=True))
    n_sub = len(substrate_order)
    substrate_y = {row.substrate_key: i / (n_sub - 1) for i, row in substrate_order.iterrows()}
    substrate_order["y"] = substrate_order["substrate_key"].map(substrate_y)

    prod_weight = df.groupby("product_clean").apply(
        lambda g: np.average([substrate_y[s] for s in g["substrate_key"]], weights=None))
    product_order = prod_weight.sort_values()
    n_prod = len(product_order)
    product_y = {p: i / (n_prod - 1) for i, p in enumerate(product_order.index)}

    df["y0"] = df["substrate_key"].map(substrate_y)
    df["y1"] = df["product_clean"].map(product_y)

    def smoothstep(t):
        return 3 * t**2 - 2 * t**3

    def curve_xy(y0, y1, n=150):
        t = np.linspace(0, 1, n)
        return t, y0 + (y1 - y0) * smoothstep(t)

    MAX_YIELD = df["yield_peak_pct"].max()
    WIDTH_MIN, WIDTH_MAX = 0.6, 6.0

    def yield_to_width(y):
        y = max(y, 0.0)
        return WIDTH_MIN + (WIDTH_MAX - WIDTH_MIN) * np.sqrt(y) / np.sqrt(MAX_YIELD)

    for _, row in df.sort_values("yield_peak_pct").iterrows():
        x, y = curve_xy(row["y0"], row["y1"])
        w = yield_to_width(row["yield_peak_pct"])
        color = MODULE_COLORS[row["module_category"]]
        ax.plot(x, y, color=color, linewidth=w, alpha=0.75, solid_capstyle="round", zorder=2)

        n_enz = int(row["n_enzymes"])
        t_dots = np.linspace(0.22, 0.78, n_enz)
        yd = row["y0"] + (row["y1"] - row["y0"]) * smoothstep(t_dots)
        ax.scatter(t_dots, yd, s=20, facecolor="white", edgecolor="#2b2b2b", linewidth=0.8, zorder=4)

    for _, row in substrate_order.iterrows():
        ax.scatter([0], [row["y"]], s=110, facecolor="#2b2b2b", edgecolor="white", linewidth=1.0, zorder=5)
        ax.annotate(f"{row['substrate_key']}  ${row['substrate_cost']:.2f}/kg", (0, row["y"]),
                    xytext=(-10, 0), textcoords="offset points", ha="right", va="center",
                    fontsize=8, fontweight="bold")

    for p, y in product_y.items():
        is_native = df.loc[df["product_clean"] == p, "product_is_native"].iloc[0]
        marker = "D" if is_native else "o"
        ax.scatter([1], [y], s=85, facecolor="#2b2b2b", edgecolor="white", linewidth=1.0,
                   marker=marker, zorder=5)
        ax.annotate(p, (1, y), xytext=(10, 0), textcoords="offset points",
                    ha="left", va="center", fontsize=7.5)

    ax.set_xlim(-0.50, 1.30)
    ax.set_ylim(-0.06, 1.06)
    ax.invert_yaxis()
    ax.axis("off")
    ax.set_title("30 heterologous pathway designs: substrate " r"$\rightarrow$" " product",
                 fontsize=12.5, fontweight="bold", pad=10)
    ax.text(0, -0.055, "SUBSTRATE (cheapest at top)", fontsize=8, style="italic", color="#555555", ha="left")
    ax.text(1, -0.055, "PRODUCT", fontsize=8, style="italic", color="#555555", ha="left")

    # anchor legends in FIGURE-fraction coords (not axes-fraction) so they land inside
    # the reserved left margin regardless of axA's own width/position in the gridspec
    module_handles = [Line2D([0], [0], color=MODULE_COLORS[m], linewidth=3.5) for m in MODULE_ORDER]
    leg1 = fig.legend(module_handles, MODULE_ORDER, title="Biosynthetic module", loc="upper left",
                       bbox_to_anchor=(0.008, 0.93), bbox_transform=fig.transFigure,
                       fontsize=7.5, title_fontsize=8.5, frameon=False)
    fig.add_artist(leg1)

    yield_examples = [0.1, 3, 15]
    width_handles = [Line2D([0], [0], color="#555555", linewidth=yield_to_width(v)) for v in yield_examples]
    leg2 = fig.legend(width_handles, [f"{v:g}% peak yield" for v in yield_examples], title="Line thickness",
                       loc="upper left", bbox_to_anchor=(0.008, 0.52), bbox_transform=fig.transFigure,
                       fontsize=7.5, title_fontsize=8.5, frameon=False)
    fig.add_artist(leg2)

    dot_handle = [Line2D([0], [0], marker="o", color="none", markerfacecolor="white",
                          markeredgecolor="#2b2b2b", markersize=6.5, linewidth=0)]
    diamond_handle = [Line2D([0], [0], marker="D", color="none", markerfacecolor="#2b2b2b",
                              markeredgecolor="white", markersize=6.5, linewidth=0)]
    circle_handle = [Line2D([0], [0], marker="o", color="none", markerfacecolor="#2b2b2b",
                             markeredgecolor="white", markersize=6.5, linewidth=0)]
    fig.legend(dot_handle + diamond_handle + circle_handle,
               ["1 heterologous enzyme", "product already native", "product is novel"],
               loc="upper left", bbox_to_anchor=(0.008, 0.34), bbox_transform=fig.transFigure,
               fontsize=7.5, frameon=False)


# =============================================================================
# Panel B: titer-range dumbbell, grouped into 3 regimes (scripts/plot_titer_ranges.py, recolored)
# =============================================================================
CEILING_THRESHOLD_MM = 0.01
FOLD_THRESHOLD = 10.0
REGIME_HIGH = "1. High & stable"
REGIME_TUNABLE = "2. Wide dynamic range"
REGIME_LOW = "3. Low & capped"
REGIME_ORDER = [REGIME_HIGH, REGIME_TUNABLE, REGIME_LOW]
REGIME_BAND_COLOR = {REGIME_HIGH: "#ffffff", REGIME_TUNABLE: "#f2f2f2", REGIME_LOW: "#e4e4e4"}


def classify_regime(max_mM, fold):
    if max_mM < CEILING_THRESHOLD_MM:
        return REGIME_LOW
    return REGIME_TUNABLE if fold >= FOLD_THRESHOLD else REGIME_HIGH


def panel_b(fig, ax):
    with open(os.path.join(ROOT, "results", "pathway_titer_grid.json"), encoding="utf-8") as f:
        grid_data = json.load(f)
    meta = pd.read_csv(os.path.join(ROOT, "results", "pathway_analysis", "merged_pathway_analysis_data.csv"))
    meta_by_id = meta.set_index("design_id").to_dict("index")
    with open(os.path.join(ROOT, "results", "pathway_optimization.json"), encoding="utf-8") as f:
        opt = {r["pid"]: r for r in json.load(f)}

    rows = []
    for entry in grid_data["results"]:
        pid = entry["pid"]
        vals = np.array([g["peak_mM"] for g in entry["grid"] if not g.get("failed") and g["peak_mM"] is not None])
        m = meta_by_id.get(pid, {})
        floor = max(vals.min(), 1e-6)
        rows.append(dict(pid=pid, module_category=m.get("module_category", "Other"),
                          min_mM=vals.min(), max_mM=vals.max(),
                          dynamic_range_fold=vals.max() / floor,
                          baseline_mM=opt[pid]["baseline_peak_mM"]))
    df = pd.DataFrame(rows)
    df["regime"] = [classify_regime(m, f) for m, f in zip(df["max_mM"], df["dynamic_range_fold"])]
    df["_regime_rank"] = df["regime"].map({r: i for i, r in enumerate(REGIME_ORDER)})
    df = df.sort_values(["_regime_rank", "max_mM"], ascending=[True, False]).reset_index(drop=True)

    GROUP_GAP = 1.5
    y_positions, band_ranges = [], {}
    cur_y = 0.0
    prev_regime = None
    for _, row in df.iterrows():
        if prev_regime is not None and row["regime"] != prev_regime:
            cur_y += GROUP_GAP
        if row["regime"] not in band_ranges:
            band_ranges[row["regime"]] = [cur_y - 0.6, cur_y]
        band_ranges[row["regime"]][1] = cur_y + 0.6
        y_positions.append(cur_y)
        cur_y += 1.0
        prev_regime = row["regime"]
    df["_y"] = y_positions

    floor_plot = 1e-4
    for regime, (y0, y1) in band_ranges.items():
        ax.axhspan(y0, y1, color=REGIME_BAND_COLOR[regime], zorder=0)

    for _, row in df.iterrows():
        i = row["_y"]
        lo = max(row["min_mM"], floor_plot)
        hi = max(row["max_mM"], floor_plot)
        c = MODULE_COLORS[row["module_category"]]
        ax.plot([lo, hi], [i, i], color=c, lw=2.0, alpha=0.6, solid_capstyle="round", zorder=2)
        ax.scatter(lo, i, color=c, s=16, zorder=3, edgecolor="white", linewidth=0.4)
        ax.scatter(hi, i, color=c, s=34, zorder=3, marker="D", edgecolor="white", linewidth=0.4)
        ax.scatter(max(row["baseline_mM"], floor_plot), i, color="black", s=20, zorder=4, marker="|")

    n_by_regime = df["regime"].value_counts()
    for regime in REGIME_ORDER:
        if regime not in band_ranges:
            continue
        header_y = band_ranges[regime][0] + 0.3
        ax.text(floor_plot * 0.55, header_y, f"{regime}  (n={n_by_regime[regime]})",
                fontsize=8, fontweight="bold", color="#333333", va="center", ha="left", zorder=5)

    ax.set_yticks(y_positions)
    ax.set_yticklabels(df["pid"], fontsize=6.6)
    ax.set_xscale("log")
    ax.set_xlabel("Peak product titer across full D x dose grid (mM)")
    ax.set_xlim(floor_plot * 0.5, 30)
    ax.set_ylim(-0.8, cur_y - 1.0 + 0.8)
    ax.invert_yaxis()
    ax.axvline(CEILING_THRESHOLD_MM, color=FIG4_GRAY, lw=0.9, ls="--", zorder=1)
    ax.set_title("Range of achievable titer per pathway\n(49-pt D x dose grid; circle=min, diamond=max, |=baseline)",
                 fontsize=11, fontweight="bold")


def main():
    fig = plt.figure(figsize=(18, 13))
    gs = fig.add_gridspec(1, 2, width_ratios=[1.25, 1.0], wspace=0.32,
                           left=0.155, right=0.98, top=0.95, bottom=0.05)
    axA = fig.add_subplot(gs[0, 0])
    axB = fig.add_subplot(gs[0, 1])
    panel_a(fig, axA)
    panel_b(fig, axB)
    add_panel_label(fig, axA, "A.")
    add_panel_label(fig, axB, "B.")
    fig.savefig(OUT_PATH, dpi=220, facecolor="white")
    print(f"Wrote {OUT_PATH}")


if __name__ == "__main__":
    main()
