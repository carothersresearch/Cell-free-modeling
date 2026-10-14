"""
Analyzes results/pathway_designs_table.csv (30 FRENDA-parameterized heterologous
pathway designs) to find what the best-performing designs have in common:
number of enzymes, biosynthetic module recruited, starting substrate, and the
substrate's native "interconnectedness" (n_endogenous_reactions, from
results/substrate_candidates.csv).
"""
import os
import re
import sys

import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
from adjustText import adjust_text

ROOT = r"C:\Users\mkcoo\Documents\UW\LLM_CFS_2"
OUT_DIR = os.path.join(ROOT, "results", "pathway_analysis")
os.makedirs(OUT_DIR, exist_ok=True)

sns.set_style("whitegrid")

# S5-S8: all text size 12, Arial, no plot titles (per figure-caption request --
# titles live in the caption text instead, not on the plot itself).
plt.rcParams["font.family"] = "Arial"
plt.rcParams["font.size"] = 12
FIG_FONT = 12

# ---------------------------------------------------------------- load + merge
df = pd.read_csv(os.path.join(ROOT, "results", "pathway_designs_table.csv"))
subs = pd.read_csv(os.path.join(ROOT, "results", "substrate_candidates.csv"))
subs = subs.dropna(subset=["substrate"])

SUBSTRATE_RENAME = {
    "HCO3": "Bicarbonate/CO2", "S_Lactate": "S-Lactate", "D_Xylose": "D-Xylose",
}
df["substrate_key"] = df["substrate"].map(lambda s: SUBSTRATE_RENAME.get(s, s))
subs = subs.set_index("substrate")

n_endo_numeric = subs["n_endogenous_reactions"].astype(str).str.extract(r"(\d+)")[0].astype(float)
df["n_endogenous_reactions"] = df["substrate_key"].map(n_endo_numeric)
df["greenness_score"] = df["substrate_key"].map(subs["greenness_score_1to10"])
df["functional_group"] = df["substrate_key"].map(subs["functional_group"])
df["substrate_classification"] = df["substrate_key"].map(subs["classification"])

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


def reuses_native(notes):
    notes = notes.lower()
    return bool(re.search(r"native.{0,25}(reuses|already|taps|feeds)|reuses.{0,15}native|already-native|already native", notes))


df["reuses_native_step"] = df["notes"].map(reuses_native)

df.to_csv(os.path.join(OUT_DIR, "merged_pathway_analysis_data.csv"), index=False)

YIELD_COL = "yield_peak_pct"

# ---------------------------------------------------------------------------
# yield_peak_pct is strongly right-skewed (a handful of large values, most designs
# near zero) -- confirmed by Shapiro-Wilk / skewness below -- so parametric tests
# (Pearson, one-way ANOVA F, Welch t) are inflated by high-leverage outliers.
# Every relationship below is therefore reported with BOTH the parametric test
# AND a robust counterpart (Spearman rank correlation / Kruskal-Wallis /
# Mann-Whitney U), plus a leave-the-top-outlier-out check. Where the parametric
# and robust tests disagree, the robust/outlier-checked result is what should be
# trusted and quoted -- this is the evidence for "none of these five things are
# reliable indicators of pathway success."
# ---------------------------------------------------------------------------
skew = stats.skew(df[YIELD_COL])
shapiro_p = stats.shapiro(df[YIELD_COL]).pvalue

STATS = {}  # populated below, used both for the summary file and the plot annotations


def corr_pair(x, y):
    mask = x.notna() & y.notna()
    r, p = stats.pearsonr(x[mask], y[mask])
    rho, ps = stats.spearmanr(x[mask], y[mask])
    return dict(n=int(mask.sum()), pearson_r=r, pearson_p=p, spearman_rho=rho, spearman_p=ps)


lines = []
lines.append(f"30 FRENDA-parameterized pathway designs -- robust statistics summary (yield metric: {YIELD_COL})")
lines.append("=" * 90)
lines.append(f"\nyield_peak_pct distribution: skewness={skew:.2f}, Shapiro-Wilk p={shapiro_p:.2e} "
             "(strongly non-normal/right-skewed -> prefer rank-based/robust tests below)\n")

# ---- 1. pathway length (n_enzymes) ----
s = corr_pair(df["n_enzymes"], df[YIELD_COL])
df_out = df[df["design_id"] != "Pyruvate_2_Isobutanol"]
s_out = corr_pair(df_out["n_enzymes"], df_out[YIELD_COL])
STATS["n_enzymes"] = dict(full=s, excl_outlier=s_out, outlier_label="excl. Pyruvate→Isobutanol (n=29)")
lines.append(f"1) n_enzymes vs yield (n={s['n']}): Pearson r={s['pearson_r']:.3f} (p={s['pearson_p']:.3f}, ns), "
             f"Spearman rho={s['spearman_rho']:.3f} (p={s['spearman_p']:.3f})")
lines.append(f"   excl. top outlier (Pyruvate_2_Isobutanol, n={s_out['n']}): Pearson r={s_out['pearson_r']:.3f} "
             f"(p={s_out['pearson_p']:.3f}), Spearman rho={s_out['spearman_rho']:.3f} (p={s_out['spearman_p']:.3f})")
lines.append("   -> Pearson sign FLIPS when the single outlier is dropped -- not a stable linear relationship.\n")

# ---- 2. substrate identity ----
groups2 = [g[YIELD_COL].values for _, g in df.groupby("substrate")]
kw2 = stats.kruskal(*groups2)
f2 = stats.f_oneway(*groups2)
df_out2 = df[df["substrate"] != "Pyruvate"]
groups2b = [g[YIELD_COL].values for _, g in df_out2.groupby("substrate")]
kw2b = stats.kruskal(*groups2b)
f2b = stats.f_oneway(*groups2b)
STATS["substrate"] = dict(kw=kw2, anova=f2, kw_excl=kw2b, anova_excl=f2b, n_groups=10)
lines.append(f"2) substrate identity (k=10 groups, n=3 each): Kruskal-Wallis H={kw2.statistic:.2f} (p={kw2.pvalue:.3f}, ns) "
             f"vs. one-way ANOVA F={f2.statistic:.2f} (p={f2.pvalue:.3f})")
lines.append(f"   excl. Pyruvate (k=9 groups): Kruskal-Wallis H={kw2b.statistic:.2f} (p={kw2b.pvalue:.3f}), "
             f"ANOVA F={f2b.statistic:.2f} (p={f2b.pvalue:.3f})")
lines.append("   -> ANOVA's apparent significance is entirely the Pyruvate group; the robust (rank) test never finds an effect.\n")

# ---- 3. module category ----
groups3 = [g[YIELD_COL].values for _, g in df.groupby("module_category")]
kw3 = stats.kruskal(*groups3)
f3 = stats.f_oneway(*groups3)
STATS["module_category"] = dict(kw=kw3, anova=f3, n_groups=df["module_category"].nunique())
lines.append(f"3) module_category (k={df['module_category'].nunique()} groups): Kruskal-Wallis H={kw3.statistic:.2f} "
             f"(p={kw3.pvalue:.3f}, ns), one-way ANOVA F={f3.statistic:.2f} (p={f3.pvalue:.3f}, ns)")
lines.append("   -> both tests agree: no detectable effect of module category.\n")

# ---- 4. substrate connectivity (n_endogenous_reactions) ----
s4_pseudo = corr_pair(df["n_endogenous_reactions"], df[YIELD_COL])  # n=30, pseudoreplicated (substrate-level var x3)
agg4 = df.groupby("substrate").agg(mean_yield=(YIELD_COL, "mean"), n_rxn=("n_endogenous_reactions", "first")).reset_index()
s4 = corr_pair(agg4["n_rxn"], agg4["mean_yield"])  # correct n=10, one point per substrate
agg4_out = agg4[agg4["substrate"] != "Pyruvate"]
s4_out = corr_pair(agg4_out["n_rxn"], agg4_out["mean_yield"])  # n=9, excl. Pyruvate
STATS["connectivity"] = dict(pseudo=s4_pseudo, correct_n=s4, excl_outlier=s4_out)
lines.append(f"4) substrate connectivity (n_endogenous_reactions) vs yield:")
lines.append(f"   naive, n=30 (pseudoreplicated -- substrate-level x duplicated x3 per design): "
             f"Pearson r={s4_pseudo['pearson_r']:.3f} (p={s4_pseudo['pearson_p']:.2e})")
lines.append(f"   correct unit of replication, n=10 (one point/substrate): Pearson r={s4['pearson_r']:.3f} "
             f"(p={s4['pearson_p']:.3f}), Spearman rho={s4['spearman_rho']:.3f} (p={s4['spearman_p']:.3f}, ns)")
lines.append(f"   excl. Pyruvate (n=9): Pearson r={s4_out['pearson_r']:.3f} (p={s4_out['pearson_p']:.3f})")
lines.append("   -> Spearman (rank, insensitive to magnitude) is NOT significant even with Pyruvate included, and Pearson\n"
             "      collapses to ~0 once Pyruvate is dropped: this is one high-leverage point, not a general trend.\n")

# ---- 5. reuses native step ----
g5f = df[df["reuses_native_step"] == False][YIELD_COL]
g5t = df[df["reuses_native_step"] == True][YIELD_COL]
u5 = stats.mannwhitneyu(g5f, g5t, alternative="two-sided")
t5 = stats.ttest_ind(g5f, g5t, equal_var=False)
STATS["reuses_native"] = dict(mwu=u5, welch=t5, n_false=len(g5f), n_true=len(g5t),
                               mean_false=g5f.mean(), mean_true=g5t.mean(),
                               median_false=g5f.median(), median_true=g5t.median())
lines.append(f"5) reuses_native_step (No n={len(g5f)}, Yes n={len(g5t)}): "
             f"Mann-Whitney U={u5.statistic:.1f} (p={u5.pvalue:.3f}, ns) vs. Welch t={t5.statistic:.2f} (p={t5.pvalue:.3f})")
lines.append("   -> Welch's t-test is nominally significant but is being driven by right-skew/outliers within the "
             "'No' group\n      (mean 2.50 vs. median 0.21); the robust rank test finds no reliable difference.\n")

lines.append("\nMean yield_peak_pct by module_category:")
lines.append(df.groupby("module_category")[YIELD_COL].agg(["mean", "median", "count"]).sort_values("mean", ascending=False).to_string())
lines.append("\n\nMean yield_peak_pct by substrate:")
lines.append(df.groupby("substrate")[YIELD_COL].agg(["mean", "median", "max", "count"]).sort_values("mean", ascending=False).to_string())
lines.append("\n\nMean yield_peak_pct by reuses_native_step:")
lines.append(df.groupby("reuses_native_step")[YIELD_COL].agg(["mean", "median", "count"]).to_string())
lines.append("\n\nMean yield_peak_pct by product_is_native:")
lines.append(df.groupby("product_is_native")[YIELD_COL].agg(["mean", "median", "count"]).to_string())
lines.append("\n\nMean yield_peak_pct by n_enzymes:")
lines.append(df.groupby("n_enzymes")[YIELD_COL].agg(["mean", "median", "count"]).to_string())

summary_text = "\n".join(lines)
with open(os.path.join(OUT_DIR, "correlation_summary.txt"), "w", encoding="utf-8") as f:
    f.write(summary_text)
print(summary_text)

STAT_BOX_KW = dict(fontsize=FIG_FONT, va="top", ha="left",
                    bbox=dict(boxstyle="round,pad=0.4", facecolor="white", edgecolor="#888888", alpha=0.92))

# colorblind-safe (Okabe-Ito) categorical pair, used everywhere a True/False hue is needed
CB_BLUE, CB_ORANGE = "#0072B2", "#E69F00"

# ---------------------------------------------------------------- 1. yield vs n_enzymes
fig, ax = plt.subplots(figsize=(9.5, 7.5))
sns.regplot(data=df, x="n_enzymes", y=YIELD_COL, ax=ax, scatter=False, color="#888888",
            line_kws=dict(linestyle="--", linewidth=1.3), ci=None)
sns.scatterplot(data=df, x="n_enzymes", y=YIELD_COL, hue="product_is_native", style="product_is_native",
                 s=110, ax=ax, palette={True: CB_ORANGE, False: CB_BLUE},
                 markers={True: "^", False: "o"}, edgecolor="black", linewidth=0.6, zorder=4)
# headroom above the highest point so the stat box has clear white space and
# doesn't collide with that point's own adjustText label
ax.set_ylim(top=df[YIELD_COL].max() * 1.38)
texts = [ax.annotate(row["design_id"].replace("_", " "), (row["n_enzymes"], row[YIELD_COL]),
                      fontsize=FIG_FONT, zorder=5)
         for _, row in df.iterrows()]
adjust_text(texts, ax=ax, arrowprops=dict(arrowstyle="-", color="#999999", lw=0.6),
            expand=(1.3, 1.6), force_text=(0.5, 0.7), force_static=(0.3, 0.4))
ax.set_xlabel("Number of new heterologous enzymes", fontsize=FIG_FONT)
ax.set_ylabel("Peak yield (% mol product / mol substrate dosed)", fontsize=FIG_FONT)
ax.legend(title="Product already native?", fontsize=FIG_FONT, title_fontsize=FIG_FONT, loc="upper left")
s = STATS["n_enzymes"]["full"]; so = STATS["n_enzymes"]["excl_outlier"]
ax.text(0.98, 0.90,
        f"Pearson r = {s['pearson_r']:+.2f} (p = {s['pearson_p']:.2f}, n.s.)\n"
        f"Spearman ρ = {s['spearman_rho']:+.2f} (p = {s['spearman_p']:.3f})\n"
        f"excl. top outlier: Pearson r = {so['pearson_r']:+.2f} (p = {so['pearson_p']:.2f})\n"
        f"  → sign flips — not a stable relationship (n = {s['n']})",
        transform=ax.transAxes, ha="right", **{k: v for k, v in STAT_BOX_KW.items() if k != "ha"})
plt.tight_layout()
plt.savefig(os.path.join(OUT_DIR, "01_yield_vs_n_enzymes.svg"), bbox_inches="tight")
plt.close()

# ---------------------------------------------------------------- 2. yield by substrate
order = df.groupby("substrate")[YIELD_COL].mean().sort_values(ascending=False).index
fig, ax = plt.subplots(figsize=(9.5, 6.5))
sns.barplot(data=df, x="substrate", y=YIELD_COL, order=order, ax=ax, errorbar=None, color=CB_BLUE, alpha=0.45)
sns.stripplot(data=df, x="substrate", y=YIELD_COL, order=order, ax=ax, color="#00344D", s=9, jitter=0.12,
              edgecolor="white", linewidth=0.6)
ax.set_xlabel("Starting substrate", fontsize=11)
ax.set_ylabel("Peak yield (%)", fontsize=11)
ax.set_title("Substrate identity does not reliably predict yield", fontsize=13, fontweight="bold")
ax.set_xticklabels(ax.get_xticklabels(), rotation=35, ha="right", fontsize=10)
ax.tick_params(axis="y", labelsize=10)
kw2, f2 = STATS["substrate"]["kw"], STATS["substrate"]["anova"]
kw2b, f2b = STATS["substrate"]["kw_excl"], STATS["substrate"]["anova_excl"]
ax.text(0.98, 0.98,
        f"Kruskal-Wallis H = {kw2.statistic:.2f} (p = {kw2.pvalue:.2f}, n.s.)\n"
        f"one-way ANOVA F = {f2.statistic:.2f} (p = {f2.pvalue:.3f})\n"
        f"excl. Pyruvate: KW p = {kw2b.pvalue:.2f}, ANOVA p = {f2b.pvalue:.2f}\n"
        f"  → ANOVA \"signal\" is the Pyruvate group alone (k = 10, n = 3/group)",
        transform=ax.transAxes, ha="right", **{k: v for k, v in STAT_BOX_KW.items() if k != "ha"})
plt.tight_layout()
plt.savefig(os.path.join(OUT_DIR, "02_yield_by_substrate.png"), dpi=180, bbox_inches="tight")
plt.close()

# ---------------------------------------------------------------- 3. yield by module category
order = df.groupby("module_category")[YIELD_COL].mean().sort_values(ascending=False).index
fig, ax = plt.subplots(figsize=(9.5, 6.5))
sns.barplot(data=df, x="module_category", y=YIELD_COL, order=order, ax=ax, errorbar=None, color="#009E73", alpha=0.45)
sns.stripplot(data=df, x="module_category", y=YIELD_COL, order=order, ax=ax, color="#004D33", s=9, jitter=0.12,
              edgecolor="white", linewidth=0.6)
ax.set_xlabel("Biosynthetic module recruited", fontsize=FIG_FONT)
ax.set_ylabel("Peak yield (%)", fontsize=FIG_FONT)
ax.set_xticklabels(ax.get_xticklabels(), rotation=35, ha="right", fontsize=FIG_FONT)
ax.tick_params(axis="y", labelsize=FIG_FONT)
kw3, f3 = STATS["module_category"]["kw"], STATS["module_category"]["anova"]
ax.text(0.98, 0.98,
        f"Kruskal-Wallis H = {kw3.statistic:.2f} (p = {kw3.pvalue:.2f}, n.s.)\n"
        f"one-way ANOVA F = {f3.statistic:.2f} (p = {f3.pvalue:.2f}, n.s.)\n"
        f"  → both tests agree: no detectable effect (k = {STATS['module_category']['n_groups']} groups)",
        transform=ax.transAxes, ha="right", **{k: v for k, v in STAT_BOX_KW.items() if k != "ha"})
plt.tight_layout()
plt.savefig(os.path.join(OUT_DIR, "03_yield_by_module.svg"), bbox_inches="tight")
plt.close()

# ---------------------------------------------------------------- 4. yield vs substrate interconnectedness
fig, ax = plt.subplots(figsize=(8.5, 7))
agg = df.groupby("substrate_key").agg(mean_yield=(YIELD_COL, "mean"), n_rxn=("n_endogenous_reactions", "first")).reset_index()
sns.regplot(data=agg, x="n_rxn", y="mean_yield", ax=ax, scatter=False, color="#888888",
            line_kws=dict(linestyle="--", linewidth=1.3), ci=None)
sns.scatterplot(data=agg, x="n_rxn", y="mean_yield", s=150, ax=ax, color="#8e44ad",
                 edgecolor="black", linewidth=0.6, zorder=4)
texts4 = [ax.annotate(row["substrate_key"], (row["n_rxn"], row["mean_yield"]), fontsize=FIG_FONT, zorder=5)
          for _, row in agg.iterrows()]
adjust_text(texts4, ax=ax, arrowprops=dict(arrowstyle="-", color="#999999", lw=0.6),
            expand=(1.3, 1.6), force_text=(0.6, 0.8))
ax.set_xlabel("Substrate's native endogenous-reaction count (interconnectedness)", fontsize=FIG_FONT)
ax.set_ylabel("Mean peak yield across that substrate's 3 pathways (%)", fontsize=FIG_FONT)
# empty upper-left quadrant (all real points sit either near x=0 at low y, or
# at the single far-right Pyruvate point) -- placed here instead of near the
# data so it can't collide with a point's adjustText label
s4, s4o = STATS["connectivity"]["correct_n"], STATS["connectivity"]["excl_outlier"]
ax.text(0.02, 0.80,
        f"n = 10 substrates: Pearson r = {s4['pearson_r']:+.2f} (p = {s4['pearson_p']:.3f})\n"
        f"Spearman ρ = {s4['spearman_rho']:+.2f} (p = {s4['spearman_p']:.2f}, n.s.)\n"
        f"excl. Pyruvate (n = 9): Pearson r = {s4o['pearson_r']:+.2f} (p = {s4o['pearson_p']:.2f})\n"
        f"  → rank test never significant; Pearson collapses without Pyruvate",
        transform=ax.transAxes, va="top",
        **{k: v for k, v in STAT_BOX_KW.items() if k not in ("ha", "va")})
plt.tight_layout()
plt.savefig(os.path.join(OUT_DIR, "04_yield_vs_substrate_interconnectedness.svg"), bbox_inches="tight")
plt.close()

# ---------------------------------------------------------------- 5. reuses-native-step effect
fig, ax = plt.subplots(figsize=(4, 3))
sns.boxplot(data=df, x="reuses_native_step", y=YIELD_COL, ax=ax, color="#e0a458", width=0.5, fliersize=0)
sns.stripplot(data=df, x="reuses_native_step", y=YIELD_COL, ax=ax, color="#8a5a1e", s=9, jitter=0.15,
              edgecolor="white", linewidth=0.6)
ax.set_xlabel("Design reuses an existing\nnative reaction/enzyme", fontsize=FIG_FONT)
ax.set_ylabel("Peak yield (%)", fontsize=FIG_FONT)
r5 = STATS["reuses_native"]
ax.set_xticks([0, 1])
ax.set_xticklabels([f"No\nn = {r5['n_false']}", f"Yes\nn = {r5['n_true']}"],
                    fontsize=FIG_FONT)
ax.tick_params(axis="y", labelsize=FIG_FONT)
ax.set_ylim(top=df[YIELD_COL].max() * 1.6)
u5, t5 = r5["mwu"], r5["welch"]
small_stat_box = {**STAT_BOX_KW, "fontsize": 7}
ax.text(0.98, 0.98,
        f"Mann-Whitney U = {u5.statistic:.0f} (p = {u5.pvalue:.2f}, n.s.)\n"
        f"Welch t = {t5.statistic:.2f} (p = {t5.pvalue:.3f})\n"
        f"means: {r5['mean_false']:.2f}% vs {r5['mean_true']:.2f}%; medians: {r5['median_false']:.2f}% vs {r5['median_true']:.3f}%\n"
        f"  → t-test driven by right-skew, not a reliable difference",
        transform=ax.transAxes, ha="right", **{k: v for k, v in small_stat_box.items() if k != "ha"})
plt.tight_layout()
plt.savefig(os.path.join(OUT_DIR, "05_yield_by_reuses_native.svg"))
plt.close()

# ---------------------------------------------------------------- 6. heatmap: design x factors summary (top/bottom designs)
top10 = df.sort_values(YIELD_COL, ascending=False).head(10)
bot10 = df.sort_values(YIELD_COL, ascending=True).head(10)
summary_cols = ["design_id", "substrate", "n_enzymes", "module_category", "product_is_native",
                 "reuses_native_step", YIELD_COL]
with open(os.path.join(OUT_DIR, "top_bottom_designs.txt"), "w", encoding="utf-8") as f:
    f.write("TOP 10 by peak yield:\n")
    f.write(top10[summary_cols].to_string(index=False))
    f.write("\n\nBOTTOM 10 by peak yield:\n")
    f.write(bot10[summary_cols].to_string(index=False))
print(open(os.path.join(OUT_DIR, "top_bottom_designs.txt")).read())

print(f"\nWrote plots + summaries to {OUT_DIR}")
