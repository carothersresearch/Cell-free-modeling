"""
Combined publication figure for the malate/rTCA design campaign, styled to match
CellFree_Fig4_v3.svg (see paper_style.py for the sampled palette/rcParams):

  Panel A: which parameters drive each objective (RF feature importance) --
           recomputed from results/design_campaign_final.csv, same model as
           scripts/analyze_designs.py Plot 3.
  Panel B: is design T5 additive? (base + 4 single-lever deltas + synergy = T5) --
           recomputed from results/analysis/factorial_design_results.csv, same
           decomposition as scripts/factorial_analysis.py.
  Panel C: citrate synthase (R46) net-synthesis flux vs. lysate dilution D --
           recomputed from results/analysis/citrate_synthase_synergy_dilution_scan.csv.
  Panel D: final malate vs. dilution D, base (native CS) vs. CS knockdown -- same
           CSV, single representative (index-60) parameter set, no jitter -- with
           the % improvement at D=10 and D=200 called out as a vertical arrow
           between the two curves at each point.
  Panel E: final malate, CS native vs. knockdown, faceted at D=200 and D=10 --
           same CSV, the concrete numbers behind panels C/D.
  Panel F: AcetylCoA time course (first 1.5 h), CS native vs. knockdown, faceted
           at D=200 and D=10 -- recomputed from
           results/analysis/citrate_synthase_synergy_timecourses.csv, the
           mechanistic evidence that CS competes with the engineered pathway for
           the same AcetylCoA/oxaloacetate pool once dilution is reduced.

C-F together visualize the "masked secondary bottleneck" argument: citrate
synthase carries negligible flux at the D=200 base case (so knocking it down
barely helps), but reducing dilution toward D=10 raises enzyme concentrations
enough to activate it as a real competing sink -- at which point knockdown
becomes strongly synergistic. See scripts/citrate_synthase_synergy_analysis.py
for the simulation that generated the panel C-F source data.

Usage:
    ./.venv_sim/Scripts/python.exe scripts/fig_malate_campaign.py
"""
import os
import sys
import json

import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestRegressor
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker

sys.path.insert(0, os.path.dirname(__file__))
from paper_style import FIG4_RED, FIG4_BLUE, FIG4_GRAY, FIG4_MODULE_COLORS, add_panel_label

ROOT = r"C:\Users\mkcoo\Documents\UW\LLM_CFS_2"
CAMPAIGN_CSV = os.path.join(ROOT, "results", "design_campaign_final.csv")
FACTORIAL_CSV = os.path.join(ROOT, "results", "analysis", "factorial_design_results.csv")
CS_SCAN_CSV = os.path.join(ROOT, "results", "analysis", "citrate_synthase_synergy_dilution_scan.csv")
CS_TIMECOURSE_CSV = os.path.join(ROOT, "results", "analysis", "citrate_synthase_synergy_timecourses.csv")
OUT_PATH = os.path.join(ROOT, "results", "analysis", "Fig_malate_campaign.svg")

# N=30 Monte Carlo kinetic-parameter jitter ensembles (scripts/campaign_jitter_lib.py,
# same sigma_log10=1/1.645 methodology as Figure 5's scripts/jitter_lib.py), precomputed
# by scripts/run_jitter_fig6a.py / _fig6b.py / _fig6cdef.py. Panels A, B, C, E, F plot
# the ensemble median with a 25th-75th percentile (IQR) band/error bar instead of the
# single index-60 point estimate, consistent with Figure 5. Panel D is the deliberate
# exception, by explicit request: it uses only the single representative (index-60)
# parameter set, no jitter, reading directly from CS_SCAN_CSV rather than JITTER_CDEF_PATH.
JITTER_A_PATH = os.path.join(ROOT, "results", "jitter_fig6a.json")
JITTER_B_PATH = os.path.join(ROOT, "results", "jitter_fig6b.json")
JITTER_CDEF_PATH = os.path.join(ROOT, "results", "jitter_fig6cdef.json")

# native = blue, CS-knockdown = red -- reused from the existing Fig4 palette
# (FIG4_BLUE/FIG4_RED) rather than introducing a third accent color, so panels
# C-F read as part of the same figure family as A/B.
NATIVE_COLOR = FIG4_BLUE
KD_COLOR = FIG4_RED

PARAM_COLS = [
    "dilution_factor_D", "hEC6411_PYC_uM", "hEC11719_FDH_uM", "hEC43117_SDS_uM",
    "eEC11137_MDH_overexpr_uM", "eEC2331_CS_mult", "eEC4212_Fumarase_mult",
    "eEC1241_PDH_mult", "eEC11142_IDH_mult", "eEC41131_PEPsynthase_mult",
    "eEC41149_PEPCK_mult", "eEC27140_PyruvateKinase_mult", "eEC11127_LDH_mult",
    "ATP_mM", "AcetylCoA_mM", "HCO3_mM", "NADH_mM", "NADPH_mM", "Formate_mM",
    "L_Serine_mM", "Glycine_added_mM", "CoA_added_mM", "duration_hr",
]
SHORT_NAME = {
    "dilution_factor_D": "Dilution (D)", "hEC6411_PYC_uM": "PYC (hEC6411)",
    "hEC11719_FDH_uM": "FDH (hEC11719)", "hEC43117_SDS_uM": "SDS (hEC43117)",
    "eEC11137_MDH_overexpr_uM": "MDH overexpr.", "eEC2331_CS_mult": "Citrate synthase x",
    "eEC4212_Fumarase_mult": "Fumarase x", "eEC1241_PDH_mult": "PDH x",
    "eEC11142_IDH_mult": "IDH x", "eEC41131_PEPsynthase_mult": "PEP synthase x",
    "eEC41149_PEPCK_mult": "PEPCK x", "eEC27140_PyruvateKinase_mult": "Pyruvate kinase x",
    "eEC11127_LDH_mult": "LDH x", "ATP_mM": "ATP (mM)", "AcetylCoA_mM": "AcetylCoA (mM)",
    "HCO3_mM": "HCO3 (mM)", "NADH_mM": "NADH (mM)", "NADPH_mM": "NADPH (mM)",
    "Formate_mM": "Formate (mM)", "L_Serine_mM": "Serine (mM)",
    "Glycine_added_mM": "Glycine added (mM)", "CoA_added_mM": "CoA added (mM)",
    "duration_hr": "Duration (h)",
}


def panel_a(ax):
    with open(JITTER_A_PATH) as f:
        j = json.load(f)
    imp_malate_ens = j["imp_malate_ensemble"]
    imp_ratio_ens = j["imp_ratio_ensemble"]

    imp_malate_med = pd.Series({c: np.nanmedian(imp_malate_ens[c]) for c in PARAM_COLS})
    imp_ratio_med = pd.Series({c: np.nanmedian(imp_ratio_ens[c]) for c in PARAM_COLS})
    imp_malate_lo = pd.Series({c: np.nanpercentile(imp_malate_ens[c], 25) for c in PARAM_COLS})
    imp_malate_hi = pd.Series({c: np.nanpercentile(imp_malate_ens[c], 75) for c in PARAM_COLS})
    imp_ratio_lo = pd.Series({c: np.nanpercentile(imp_ratio_ens[c], 25) for c in PARAM_COLS})
    imp_ratio_hi = pd.Series({c: np.nanpercentile(imp_ratio_ens[c], 75) for c in PARAM_COLS})

    imp_df = pd.DataFrame({"Final malate": imp_malate_med, "log10(rTCA:TCA)": imp_ratio_med})
    imp_df["combined"] = imp_df.mean(axis=1)
    order = imp_df.sort_values("combined", ascending=True).index
    imp_df_plot = imp_df.loc[order].drop(columns="combined").rename(index=SHORT_NAME)

    y = np.arange(len(imp_df_plot))
    bar_h = 0.38
    malate_err = [(imp_malate_med - imp_malate_lo).loc[order].values, (imp_malate_hi - imp_malate_med).loc[order].values]
    ratio_err = [(imp_ratio_med - imp_ratio_lo).loc[order].values, (imp_ratio_hi - imp_ratio_med).loc[order].values]
    ax.barh(y + bar_h / 2, imp_df_plot["Final malate"], height=bar_h, color=FIG4_RED, label="Final malate (mM)",
            xerr=malate_err, capsize=1.5, error_kw={"linewidth": 0.7})
    ax.barh(y - bar_h / 2, imp_df_plot["log10(rTCA:TCA)"], height=bar_h, color=FIG4_BLUE, label="log10(rTCA:TCA ratio)",
            xerr=ratio_err, capsize=1.5, error_kw={"linewidth": 0.7})
    ax.set_yticks(y)
    ax.set_yticklabels(imp_df_plot.index, fontsize=6.5)
    ax.tick_params(axis="x", labelsize=7.5)
    ax.tick_params(axis="y", length=2, pad=2)
    ax.set_xlabel("Random Forest feature importance (median, IQR; n=30 jittered draws)", fontsize=8)
    ax.set_title("Which parameters drive each objective?", fontsize=10, fontweight="bold")
    ax.legend(frameon=False, loc="lower right", fontsize=7)
    ax.set_ylim(-0.8, len(imp_df_plot) - 0.2)


def panel_b(ax):
    with open(JITTER_B_PATH) as f:
        j = json.load(f)
    factor_names = ["Dilution", "CS_KD", "Fum_KD", "Serine"]

    base_arr = np.array(j["base_malate"])
    t5_arr = np.array(j["t5_malate"])
    sum_singles_arr = np.array(j["sum_singles"])
    synergy_arr = np.array(j["synergy"])
    single_arrs = {f: np.array(j["single_deltas"][f]) for f in factor_names}

    def med_iqr(arr):
        return np.nanmedian(arr), np.nanpercentile(arr, 25), np.nanpercentile(arr, 75)

    base_malate, base_lo, base_hi = med_iqr(base_arr)
    t5_malate, t5_lo, t5_hi = med_iqr(t5_arr)
    sum_singles, _, _ = med_iqr(sum_singles_arr)
    synergy, syn_lo, syn_hi = med_iqr(synergy_arr)
    synergy_pct = synergy / sum_singles * 100

    single_med = {f: med_iqr(single_arrs[f])[0] for f in factor_names}
    single_lo = {f: med_iqr(single_arrs[f])[1] for f in factor_names}
    single_hi = {f: med_iqr(single_arrs[f])[2] for f in factor_names}

    # narrative order (matches design story / CLAUDE.md emphasis): the two
    # dominant single levers first (Dilution, Serine), then the two cheap
    # secondary knockdown levers (CS KD, Fumarase KD). Each step's value is the
    # ensemble median; err = (lo, hi) IQR bounds for that step's own quantity.
    steps = [
        ("Base\n(D=200)", base_malate, base_lo, base_hi, "base"),
        ("+Dilution\nD=10", single_med["Dilution"], single_lo["Dilution"], single_hi["Dilution"], "lever"),
        ("+Serine\n50mM", single_med["Serine"], single_lo["Serine"], single_hi["Serine"], "lever"),
        ("+CS KD\n1%", single_med["CS_KD"], single_lo["CS_KD"], single_hi["CS_KD"], "lever"),
        ("+Fum KD\n1%", single_med["Fum_KD"], single_lo["Fum_KD"], single_hi["Fum_KD"], "lever"),
        ("Synergy\n(interact.)", synergy, syn_lo, syn_hi, "synergy"),
        ("T5\n(all four)", t5_malate, t5_lo, t5_hi, "total"),
    ]

    color_by_kind = {"base": FIG4_GRAY, "lever": FIG4_MODULE_COLORS["Glycolysis/Gluconeogenesis"],
                      "synergy": FIG4_RED, "total": FIG4_MODULE_COLORS["TCA Cycle"]}

    cum = 0.0
    labels, lows, heights, colors, err_lo, err_hi = [], [], [], [], [], []
    for label, val, lo, hi, kind in steps:
        labels.append(label)
        colors.append(color_by_kind[kind])
        if kind in ("base", "total"):
            lows.append(0.0)
            heights.append(val)
            err_lo.append(val - lo)
            err_hi.append(hi - val)
            cum = val
        else:
            lows.append(cum)
            heights.append(val)
            err_lo.append(val - lo)
            err_hi.append(hi - val)
            cum += val

    x = np.arange(len(steps))
    ax.bar(x, heights, bottom=lows, color=colors, width=0.62, edgecolor="none", zorder=3)
    ax.errorbar(x, [lows[i] + heights[i] for i in range(len(steps))], yerr=[err_lo, err_hi],
                fmt="none", ecolor="#333333", elinewidth=0.8, capsize=2, zorder=4)

    for i in range(len(steps) - 1):
        top_i = lows[i] + max(heights[i], 0) if heights[i] >= 0 else lows[i]
        ax.plot([x[i] + 0.31, x[i + 1] - 0.31], [top_i, top_i], color="#888888", lw=1.0, ls=":", zorder=2)

    for i, (label, val, lo, hi, kind) in enumerate(steps):
        top = lows[i] + heights[i] + err_hi[i]
        va = "bottom" if heights[i] >= 0 else "top"
        txt = f"{val:.3f}" if kind == "base" else f"{val:+.3f}"
        ax.text(x[i], top + (0.05 if heights[i] >= 0 else -0.05), txt, ha="center", va=va,
                fontsize=6.5, fontweight="bold" if kind in ("base", "total") else "normal")

    ax.set_xticks(x)
    ax.set_xticklabels(labels, fontsize=6.5, rotation=0, ha="center")
    ax.tick_params(axis="x", pad=3)
    ax.tick_params(axis="y", labelsize=7)
    ax.set_ylabel("Final malate (mM)", fontsize=8)
    ax.set_title(f"Is T5 additive? sum of parts = {sum_singles:.3f} mM vs. actual T5 = "
                 f"{t5_malate:.3f} mM ({synergy_pct:+.1f}% synergy)\n(median, IQR; n=30 jittered draws)",
                 fontsize=8.5, fontweight="bold")
    ax.set_ylim(0, (t5_malate + err_hi[-1]) * 1.2)


def _annotate_focal_points(ax, x_vals, y_vals, focal_x, fmt):
    """Label just the D=10 and D=200 points on a C/D line plot. Both curves are
    decaying from a max near D=10, so the y-axis headroom (set by the caller to
    1.25x the max) always has room for an upward label -- only the horizontal
    anchor needs to flip for the D=10 point, which sits near the left edge."""
    x_vals = list(x_vals)
    for fx in focal_x:
        idx = x_vals.index(fx)
        xv, yv = x_vals[idx], y_vals[idx]
        color = ax.lines[-1].get_color()
        ax.scatter([xv], [yv], s=30, color=color, zorder=4, edgecolor="white", linewidth=0.8)
        near_left_edge = idx == 0
        ha = "left" if near_left_edge else "center"
        xytext = (4, 8) if near_left_edge else (0, 8)
        ax.annotate(fmt(yv), (xv, yv), textcoords="offset points", xytext=xytext,
                    ha=ha, fontsize=7.5, fontweight="bold")


def _load_cdef():
    with open(JITTER_CDEF_PATH) as f:
        return json.load(f)


def panel_c(ax):
    """Citrate synthase net-synthesis flux vs. dilution D -- median + IQR band."""
    j = _load_cdef()
    D_SCAN = sorted(j["D_SCAN"])
    med = np.array([np.nanmedian(j["scan_ensemble"][str(D)]["mean_CS_sink_flux_wt"]) for D in D_SCAN])
    lo = np.array([np.nanpercentile(j["scan_ensemble"][str(D)]["mean_CS_sink_flux_wt"], 25) for D in D_SCAN])
    hi = np.array([np.nanpercentile(j["scan_ensemble"][str(D)]["mean_CS_sink_flux_wt"], 75) for D in D_SCAN])
    x = np.array(D_SCAN)
    ax.plot(x, med, color=NATIVE_COLOR, lw=1.8, marker="o", ms=3.5, zorder=3)
    ax.fill_between(x, lo, hi, color=NATIVE_COLOR, alpha=0.2, linewidth=0, zorder=2)
    _annotate_focal_points(ax, x, med, [10, 200], lambda v: f"{v:.2e}")
    ax.set_xlabel("Dilution factor D", fontsize=9)
    ax.set_ylabel("Mean CS synthesis flux (mM/s)", fontsize=8.5)
    ax.set_title("CS flux rises sharply as\ndilution decreases\n(median, IQR; n=30 jittered draws)", fontsize=8.5, fontweight="bold")
    # explicit tick values + a formatter (rather than ax.ticklabel_format's
    # auto-placed "1e-3" offset text) -- the offset text anchors above the
    # y-axis in the same spot as the panel label and collided with it.
    ymax = max(np.nanmax(hi), np.nanmax(med)) * 1.35
    ax.set_yticks(np.linspace(0, ymax, 4))
    ax.yaxis.set_major_formatter(mticker.FuncFormatter(lambda v, pos: "0" if v == 0 else f"{v * 1000:.1f}e-3"))
    ax.tick_params(axis="both", labelsize=7.5)
    ax.set_xlim(0, 410)
    ax.set_ylim(0, ymax)


D_POINTS_PANEL_D = [10, 50, 100, 150, 200]


def panel_d(ax):
    """Final malate, base (native CS) vs. CS knockdown, at five dilution factors
    (10X, 50X, 100X, 150X, 200X) -- paired bar chart, single representative
    (index-60) parameter set, no jitter. % improvement is annotated above each
    pair of bars."""
    df = pd.read_csv(CS_SCAN_CSV).set_index("dilution_D")
    wt = np.array([df.loc[d, "final_malate_mM_wt"] for d in D_POINTS_PANEL_D])
    kd = np.array([df.loc[d, "final_malate_mM_kd"] for d in D_POINTS_PANEL_D])
    pct = np.array([df.loc[d, "cs_kd_pct_improvement"] for d in D_POINTS_PANEL_D])

    x = np.arange(len(D_POINTS_PANEL_D))
    width = 0.35
    ax.bar(x - width / 2, wt, width, color=NATIVE_COLOR, zorder=3, label="Base")
    ax.bar(x + width / 2, kd, width, color=KD_COLOR, zorder=3, label="CS KD")

    ymax = max(wt.max(), kd.max()) * 1.3
    for i, (w, k, p) in enumerate(zip(wt, kd, pct)):
        ax.text(x[i], max(w, k) + ymax * 0.03, f"+{p:.1f}%", ha="center", va="bottom",
                fontsize=7.5, fontweight="bold", color="#333333")

    ax.set_xticks(x)
    ax.set_xticklabels([f"{d}X" for d in D_POINTS_PANEL_D], fontsize=7.5)
    ax.set_xlabel("Dilution factor", fontsize=9)
    ax.set_ylabel("Final malate (mM)", fontsize=8.5)
    ax.set_title("...and that flux translates into a\nreal knockdown benefit at low dilution\n(single representative parameter set)", fontsize=8.5, fontweight="bold")
    ax.tick_params(axis="both", labelsize=7.5)
    ax.set_ylim(0, ymax)
    ax.legend(fontsize=7, frameon=False, loc="upper right", handlelength=1.5)


def panel_e(fig, spec):
    """Final malate, native vs. CS-knockdown, faceted at D=200 and D=10 -- median + IQR error bar."""
    j = _load_cdef()
    sub = spec.subgridspec(1, 2, wspace=0.55)
    axes = [fig.add_subplot(sub[0, 0]), fig.add_subplot(sub[0, 1])]
    facets = [(200, "D = 200\n(base case)"), (10, "D = 10\n(most concentrated)")]

    for ax, (dval, title) in zip(axes, facets):
        row = j["scan_ensemble"][str(dval)]
        wt = np.array(row["final_malate_wt"])
        kd = np.array(row["final_malate_kd"])
        med = [np.nanmedian(wt), np.nanmedian(kd)]
        lo = [np.nanpercentile(wt, 25), np.nanpercentile(kd, 25)]
        hi = [np.nanpercentile(wt, 75), np.nanpercentile(kd, 75)]
        err = [[m - l for m, l in zip(med, lo)], [h - m for m, h in zip(med, hi)]]
        colors = [NATIVE_COLOR, KD_COLOR]
        bars = ax.bar([0, 1], med, yerr=err, capsize=3, width=0.6, color=colors, zorder=3)
        for b, v in zip(bars, med):
            ax.text(b.get_x() + b.get_width() / 2, v, f"{v:.3f}", ha="center", va="bottom", fontsize=7.5)
        pct_med = np.nanmedian(row["cs_kd_pct_improvement"])
        ax.text(0.5, max(hi) * 1.22, f"+{pct_med:.1f}%",
                ha="center", fontsize=8, fontweight="bold", color="#333333")
        ax.set_xticks([0, 1])
        ax.set_xticklabels(["Native", "Knockdown"], fontsize=7.5)
        ax.set_title(title, fontsize=8.5, fontweight="bold")
        ax.set_ylim(0, max(hi) * 1.4)
        ax.tick_params(axis="y", labelsize=7.5)
    axes[0].set_ylabel("Final malate (mM)\n(median, IQR; n=30 jittered draws)", fontsize=7.5)
    return axes


def panel_f(fig, spec):
    """AcetylCoA time course, native vs. CS-knockdown, faceted at D=200/D=10 -- median + IQR band."""
    j = _load_cdef()
    sub = spec.subgridspec(1, 2, wspace=0.55)
    axes = [fig.add_subplot(sub[0, 0]), fig.add_subplot(sub[0, 1])]
    facets = [(200, "D = 200"), (10, "D = 10")]

    for ax, (dval, title) in zip(axes, facets):
        native_draws = j["timecourse_ensemble"][f"D{dval}_native"]
        kd_draws = j["timecourse_ensemble"][f"D{dval}_kd"]
        t = np.array(native_draws[0]["time_h"])
        mask = t <= 1.5

        for draws, color, label in [(native_draws, NATIVE_COLOR, "Native"), (kd_draws, KD_COLOR, "Knockdown")]:
            vals = np.array([draw["AcetylCoA_mM"] for draw in draws])  # (n_draws, n_time)
            med = np.nanmedian(vals, axis=0)
            lo = np.nanpercentile(vals, 25, axis=0)
            hi = np.nanpercentile(vals, 75, axis=0)
            ax.plot(t[mask], med[mask], color=color, lw=1.6, label=label)
            ax.fill_between(t[mask], lo[mask], hi[mask], color=color, alpha=0.2, linewidth=0)

        ax.set_title(title, fontsize=8.5, fontweight="bold")
        ax.set_xlabel("time (h)", fontsize=8)
        ax.set_ylim(0, 1.05)
        ax.tick_params(axis="both", labelsize=7.5)
    axes[0].set_ylabel("AcetylCoA (mM)\n(median, IQR; n=30 jittered draws)", fontsize=7.5)
    # lower-left is empty in every facet: both curves start high (near t=0) and
    # only fall over time, never approaching y=0 -- upper-right/upper-left both
    # collide with the curves' early, still-high values.
    axes[1].legend(frameon=False, fontsize=7, loc="lower left", handlelength=1.4)
    return axes


def main():
    # Portrait layout sized to fit within 1" margins of an A4 page:
    # 8.27 - 2*1 = 6.27" wide, 11.69 - 2*1 = 9.69" tall. Six panels stacked in four
    # rows -- A (RF importance, 22 param labels) and B (waterfall) each need the
    # full width to stay readable; C/D and E/F are two-up rows of the new
    # dilution-dependent citrate-synthase panels.
    #
    # fig.tight_layout() is deliberately NOT used: it warned "Axes that are not
    # compatible with tight_layout" for panel_e/f's nested subgridspecs and, in
    # practice, produced a badly broken layout (panel A's 22 labels crushed into
    # a sliver, panel B's tick labels colliding with panel C/D's titles). A
    # GridSpec with explicit left/right/top/bottom/hspace/wspace instead gives
    # every axes a deterministic position up front, which nested subgridspecs
    # (panel_e/f) inherit correctly with no post-hoc reflow needed.
    fig = plt.figure(figsize=(6.27, 9.69))
    gs = fig.add_gridspec(
        4, 2, height_ratios=[2.05, 0.85, 0.80, 0.80],
        left=0.155, right=0.985, top=0.975, bottom=0.045,
        hspace=0.70, wspace=0.45,
    )

    axA = fig.add_subplot(gs[0, :])
    axB = fig.add_subplot(gs[1, :])
    axC = fig.add_subplot(gs[2, 0])
    axD = fig.add_subplot(gs[2, 1])

    panel_a(axA)
    panel_b(axB)
    panel_c(axC)
    panel_d(axD)
    axesE = panel_e(fig, gs[3, 0])
    axesF = panel_f(fig, gs[3, 1])

    def add_small_panel_label(ax, letter):
        # paper_style.add_panel_label's 18pt label is sized for the full-width
        # A/B panels; at that size it visually intrudes on the smaller C-F
        # quadrant panels, so those get a scaled-down label instead. The
        # topmost y-tick in C/D/E/F sits right at the axes' top edge (their
        # ylim max doubles as a tick), so the label needs real clearance
        # above bbox.y1, not just a token offset, to avoid sitting on top of it.
        bbox = ax.get_position()
        fig.text(bbox.x0 - 0.028, bbox.y1 + 0.026, letter, fontsize=12, fontweight="bold", ha="left", va="top")

    add_panel_label(fig, axA, "A.")
    add_panel_label(fig, axB, "B.")
    add_small_panel_label(axC, "C.")
    add_small_panel_label(axD, "D.")
    add_small_panel_label(axesE[0], "E.")
    add_small_panel_label(axesF[0], "F.")

    fig.savefig(OUT_PATH, dpi=220, facecolor="white")
    print(f"Wrote {OUT_PATH}")


if __name__ == "__main__":
    main()
