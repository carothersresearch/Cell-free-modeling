"""
Clustering / parameter-variability / correlation analysis of the malate-rTCA
engineering-strategy design campaign (results/design_campaign_final.csv).

Answers:
  1) Do the parameters separate designs into natural clusters, and how do
     those data-driven clusters compare to the mechanistic tags assigned
     during design generation?
  2) Which parameters are "narrow range for success" (only a tight band of
     values appears among high performers) vs "loosely deterministic" (high
     performers exist across nearly the whole explored range)?
  3) Which parameters matter most for each objective (variability-driven
     importance, via Random Forest)?
  4) Which parameter pairs trade off against each other (correlated, esp.
     negatively, across the design set)?

Produces PNGs in results/analysis/ and a text summary
(results/analysis/summary.txt) with the numeric findings behind each plot.

Usage:
    ./.venv_sim/Scripts/python.exe scripts/analyze_designs.py
"""
import warnings
warnings.filterwarnings("ignore")

from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import seaborn as sns

from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score, adjusted_rand_score
from sklearn.ensemble import RandomForestRegressor
import umap

# ---------------------------------------------------------------------------
# Setup
# ---------------------------------------------------------------------------
ROOT = Path(__file__).resolve().parent.parent
DATA_PATH = ROOT / "results" / "design_campaign_final.csv"
OUT_DIR = ROOT / "results" / "analysis"
OUT_DIR.mkdir(exist_ok=True)

# Okabe-Ito colorblind-safe categorical palette (fixed order, never cycled/reassigned)
OKABE_ITO = ["#E69F00", "#56B4E9", "#009E73", "#D55E00", "#CC79A7", "#0072B2", "#F0E442", "#999999"]
DIVERGING_CMAP = "vlag"  # blue-white-red, neutral midpoint

sns.set_theme(style="white", font_scale=1.0)

PARAM_COLS = [
    "dilution_factor_D", "hEC6411_PYC_uM", "hEC11719_FDH_uM", "hEC43117_SDS_uM",
    "eEC11137_MDH_overexpr_uM", "eEC2331_CS_mult", "eEC4212_Fumarase_mult",
    "eEC1241_PDH_mult", "eEC11142_IDH_mult", "eEC41131_PEPsynthase_mult",
    "eEC41149_PEPCK_mult", "eEC27140_PyruvateKinase_mult", "eEC11127_LDH_mult",
    "ATP_mM", "AcetylCoA_mM", "HCO3_mM", "NADH_mM", "NADPH_mM", "Formate_mM",
    "L_Serine_mM", "Glycine_added_mM", "CoA_added_mM", "duration_hr",
]
SHORT_NAME = {
    "dilution_factor_D": "Dilution (D)",
    "hEC6411_PYC_uM": "PYC (hEC6411)",
    "hEC11719_FDH_uM": "FDH (hEC11719)",
    "hEC43117_SDS_uM": "SDS (hEC43117)",
    "eEC11137_MDH_overexpr_uM": "MDH overexpr.",
    "eEC2331_CS_mult": "Citrate synthase x",
    "eEC4212_Fumarase_mult": "Fumarase x",
    "eEC1241_PDH_mult": "PDH x",
    "eEC11142_IDH_mult": "IDH x",
    "eEC41131_PEPsynthase_mult": "PEP synthase x",
    "eEC41149_PEPCK_mult": "PEPCK x",
    "eEC27140_PyruvateKinase_mult": "Pyruvate kinase x",
    "eEC11127_LDH_mult": "LDH x",
    "ATP_mM": "ATP (mM)",
    "AcetylCoA_mM": "AcetylCoA (mM)",
    "HCO3_mM": "HCO3 (mM)",
    "NADH_mM": "NADH (mM)",
    "NADPH_mM": "NADPH (mM)",
    "Formate_mM": "Formate (mM)",
    "L_Serine_mM": "Serine (mM)",
    "Glycine_added_mM": "Glycine added (mM)",
    "CoA_added_mM": "CoA added (mM)",
    "duration_hr": "Duration (h)",
}

# ---------------------------------------------------------------------------
# Load + prepare
# ---------------------------------------------------------------------------
df = pd.read_csv(DATA_PATH)
base_mask = df["design_id"].str.contains("base")
base_row = df[base_mask].iloc[0]
d = df[~base_mask].reset_index(drop=True).copy()

# cap infinite ratios for numeric handling (use 10x the max finite value observed)
finite_ratios = d.loc[np.isfinite(d["rTCA_to_TCA_ratio"]), "rTCA_to_TCA_ratio"]
cap = finite_ratios.max() * 10
d["ratio_capped"] = d["rTCA_to_TCA_ratio"].replace([np.inf], cap)
d["log10_ratio"] = np.log10(d["ratio_capped"] + 1)
d["log10_malate"] = np.log10(d["final_malate_mM"])

d["n_tags"] = d["tags"].str.split("|").apply(len)
d["primary_tag"] = d["tags"].str.split("|").str[0]
lever_group = pd.cut(d["n_tags"], bins=[0, 1, 2, 3, 4, 100], labels=["1 lever", "2 levers", "3 levers", "4 levers", "5+ levers"])
d["lever_group"] = lever_group

# log-transform every parameter (all are naturally multiplicative / span orders
# of magnitude: dilution factor, fold-change multipliers, and concentrations
# that were dosed across wide ranges); log1p for anything that can be exactly 0.
X_log = pd.DataFrame(index=d.index)
for c in PARAM_COLS:
    vals = d[c].values.astype(float)
    X_log[c] = np.log1p(vals) if (vals <= 0).any() else np.log(vals)

scaler = StandardScaler()
X_scaled = scaler.fit_transform(X_log.values)
X_scaled_df = pd.DataFrame(X_scaled, columns=PARAM_COLS, index=d.index)

summary_lines = []
def log(msg=""):
    print(msg)
    summary_lines.append(msg)

log("=" * 70)
log("Malate/rTCA design campaign -- clustering & variability analysis")
log(f"n designs analyzed: {len(d)} (base case excluded from clustering)")
log("=" * 70)

# ---------------------------------------------------------------------------
# PCA + UMAP
# ---------------------------------------------------------------------------
pca = PCA(n_components=2, random_state=42)
pca_xy = pca.fit_transform(X_scaled)
log(f"\nPCA explained variance ratio: PC1={pca.explained_variance_ratio_[0]:.3f}, "
    f"PC2={pca.explained_variance_ratio_[1]:.3f} (total={pca.explained_variance_ratio_[:2].sum():.3f})")

loadings = pd.DataFrame(pca.components_.T, index=PARAM_COLS, columns=["PC1", "PC2"])
log("\nTop PC1 loadings (what PC1 mostly represents):")
for name, val in loadings["PC1"].reindex(loadings["PC1"].abs().sort_values(ascending=False).index).head(6).items():
    log(f"  {SHORT_NAME[name]:22s} {val:+.3f}")
log("Top PC2 loadings:")
for name, val in loadings["PC2"].reindex(loadings["PC2"].abs().sort_values(ascending=False).index).head(6).items():
    log(f"  {SHORT_NAME[name]:22s} {val:+.3f}")

reducer = umap.UMAP(n_components=2, random_state=42, n_neighbors=15, min_dist=0.3)
umap_xy = reducer.fit_transform(X_scaled)

# ---------------------------------------------------------------------------
# KMeans clustering, k chosen by silhouette score
# ---------------------------------------------------------------------------
sil_scores = {}
labels_by_k = {}
for k in range(2, 9):  # capped at 8 to match the fixed 8-color categorical palette
    km = KMeans(n_clusters=k, random_state=42, n_init=10)
    labels = km.fit_predict(X_scaled)
    sil_scores[k] = silhouette_score(X_scaled, labels)
    labels_by_k[k] = labels
best_k = max(sil_scores, key=sil_scores.get)
d["cluster"] = labels_by_k[best_k]

log(f"\nSilhouette scores by k: " + ", ".join(f"k={k}:{s:.3f}" for k, s in sil_scores.items()))
log(f"Selected k={best_k} (highest silhouette={sil_scores[best_k]:.3f})")
log("\nCluster sizes and mean objectives:")
for c in sorted(d["cluster"].unique()):
    sub = d[d["cluster"] == c]
    log(f"  cluster {c}: n={len(sub):3d}  mean malate={sub['final_malate_mM'].mean():.3f} mM  "
        f"median log10(ratio)={sub['log10_ratio'].median():.2f}  "
        f"top tags={sub['tags'].str.split('|').explode().value_counts().head(3).to_dict()}")

# ---------------------------------------------------------------------------
# Plot 1: PCA + UMAP scatter (color=cluster, size=malate, shape=lever count)
# ---------------------------------------------------------------------------
marker_list = ["o", "s", "^", "D", "P", "X"]
lever_categories = list(d["lever_group"].cat.categories)
marker_map = {cat: marker_list[i % len(marker_list)] for i, cat in enumerate(lever_categories)}
cluster_colors = {c: OKABE_ITO[i % len(OKABE_ITO)] for i, c in enumerate(sorted(d["cluster"].unique()))}

size_min, size_max = 40, 500
malate_norm = (d["final_malate_mM"] - d["final_malate_mM"].min()) / (d["final_malate_mM"].max() - d["final_malate_mM"].min())
sizes = size_min + malate_norm * (size_max - size_min)

fig, (legend_ax, ax0, ax1) = plt.subplots(1, 3, figsize=(17, 6.5), gridspec_kw={"width_ratios": [0.42, 1, 1]})
legend_ax.axis("off")
axes = [ax0, ax1]
for ax, (xy, title) in zip(axes, [(pca_xy, f"PCA (PC1 {pca.explained_variance_ratio_[0]*100:.0f}%, PC2 {pca.explained_variance_ratio_[1]*100:.0f}%)"),
                                   (umap_xy, "UMAP")]):
    for lever_cat in lever_categories:
        mask = d["lever_group"] == lever_cat
        if mask.sum() == 0:
            continue
        colors = [cluster_colors[c] for c in d.loc[mask, "cluster"]]
        ax.scatter(xy[mask.values, 0], xy[mask.values, 1], s=sizes[mask], c=colors,
                   marker=marker_map[lever_cat], alpha=0.8, edgecolor="white", linewidth=0.6,
                   label=lever_cat)
    ax.set_title(title, fontsize=12)
    ax.set_xlabel("Dim 1")
    ax.set_ylabel("Dim 2")
    sns.despine(ax=ax)

cluster_legend = [Line2D([0], [0], marker="o", color="w", markerfacecolor=cluster_colors[c],
                          markersize=10, label=f"Cluster {c}") for c in sorted(cluster_colors)]
shape_legend = [Line2D([0], [0], marker=marker_map[cat], color="gray", linestyle="",
                        markersize=9, label=str(cat)) for cat in lever_categories]
size_legend = [Line2D([0], [0], marker="o", color="gray", linestyle="", markersize=np.sqrt(s / 3.5),
                       label=lbl) for s, lbl in [(size_min, f"{d['final_malate_mM'].min():.2f} mM"),
                                                   (size_max, f"{d['final_malate_mM'].max():.2f} mM")]]
leg1 = legend_ax.legend(handles=cluster_legend, loc="upper left", bbox_to_anchor=(0.0, 1.0),
                         title="Color = cluster", frameon=False, fontsize=9, title_fontsize=10)
legend_ax.add_artist(leg1)
leg2 = legend_ax.legend(handles=shape_legend, loc="upper left", bbox_to_anchor=(0.0, 0.55),
                         title="Shape = # levers changed", frameon=False, fontsize=9, title_fontsize=10)
legend_ax.add_artist(leg2)
legend_ax.legend(handles=size_legend, loc="upper left", bbox_to_anchor=(0.0, 0.20),
                  title="Size = final malate", frameon=False, fontsize=9, title_fontsize=10)
fig.suptitle("Design-space embedding: cluster identity, lever count, and malate yield", fontsize=13)
plt.tight_layout()
plt.savefig(OUT_DIR / "01_pca_umap_scatter.png", dpi=180, bbox_inches="tight")
plt.close()

# ---------------------------------------------------------------------------
# Plot 2: clustered heatmap of parameters
# ---------------------------------------------------------------------------
row_colors = d["cluster"].map(cluster_colors)
tag_group_colors_map = {cat: OKABE_ITO[(i + 3) % len(OKABE_ITO)] for i, cat in enumerate(lever_categories)}
row_colors2 = d["lever_group"].map(tag_group_colors_map)
row_color_df = pd.DataFrame({"Cluster": row_colors.values, "# levers": row_colors2.values}, index=d.index)

heat_df = X_scaled_df.rename(columns=SHORT_NAME)
g = sns.clustermap(heat_df, cmap=DIVERGING_CMAP, center=0, vmin=-2.5, vmax=2.5,
                    row_colors=row_color_df, figsize=(13, 16),
                    dendrogram_ratio=(0.12, 0.12), cbar_pos=(0.02, 0.83, 0.02, 0.15),
                    xticklabels=True, yticklabels=False)
g.ax_heatmap.set_xlabel("Parameter (z-scored, log-transformed)")
g.ax_heatmap.set_ylabel(f"Designs (n={len(d)}, hierarchically clustered)")
plt.setp(g.ax_heatmap.get_xticklabels(), rotation=45, ha="right")
g.fig.suptitle("Clustered heatmap of design parameters", y=1.02, fontsize=13)

cluster_patches = [plt.matplotlib.patches.Patch(facecolor=cluster_colors[c], label=f"Cluster {c}") for c in sorted(cluster_colors)]
lever_patches = [plt.matplotlib.patches.Patch(facecolor=tag_group_colors_map[cat], label=str(cat)) for cat in lever_categories]
leg_a = g.fig.legend(handles=cluster_patches, title="Cluster", loc="upper left", bbox_to_anchor=(1.005, 0.98), frameon=False, fontsize=9, title_fontsize=10)
g.fig.add_artist(leg_a)
g.fig.legend(handles=lever_patches, title="# levers", loc="upper left", bbox_to_anchor=(1.005, 0.66), frameon=False, fontsize=9, title_fontsize=10)

g.savefig(OUT_DIR / "02_parameter_clustermap.png", dpi=180, bbox_inches="tight")
plt.close()

# ---------------------------------------------------------------------------
# Plot 3: parameter importance (Random Forest, variability-driven) for each objective
# ---------------------------------------------------------------------------
rf_malate = RandomForestRegressor(n_estimators=500, random_state=42, max_depth=6)
rf_malate.fit(X_scaled, d["final_malate_mM"])
imp_malate = pd.Series(rf_malate.feature_importances_, index=PARAM_COLS)

rf_ratio = RandomForestRegressor(n_estimators=500, random_state=42, max_depth=6)
rf_ratio.fit(X_scaled, d["log10_ratio"])
imp_ratio = pd.Series(rf_ratio.feature_importances_, index=PARAM_COLS)

imp_df = pd.DataFrame({"Final malate": imp_malate, "log10(rTCA:TCA)": imp_ratio})
imp_df["combined"] = imp_df.mean(axis=1)
imp_df = imp_df.sort_values("combined", ascending=True)
imp_df_plot = imp_df.drop(columns="combined").rename(index=SHORT_NAME)

fig, ax = plt.subplots(figsize=(9, 9))
y = np.arange(len(imp_df_plot))
bar_h = 0.38
ax.barh(y + bar_h / 2, imp_df_plot["Final malate"], height=bar_h, color=OKABE_ITO[0], label="Final malate (mM)")
ax.barh(y - bar_h / 2, imp_df_plot["log10(rTCA:TCA)"], height=bar_h, color=OKABE_ITO[1], label="log10(rTCA:TCA ratio)")
ax.set_yticks(y)
ax.set_yticklabels(imp_df_plot.index)
ax.set_xlabel("Random Forest feature importance (variability-driven)")
ax.set_title("Which parameters drive each objective?", fontsize=13)
ax.legend(frameon=False, loc="lower right")
sns.despine(ax=ax)
plt.tight_layout()
plt.savefig(OUT_DIR / "03_parameter_importance.png", dpi=180, bbox_inches="tight")
plt.close()

log("\nTop 5 parameters by importance for final malate:")
for name, val in imp_malate.sort_values(ascending=False).head(5).items():
    log(f"  {SHORT_NAME[name]:22s} {val:.3f}")
log("Top 5 parameters by importance for rTCA:TCA ratio:")
for name, val in imp_ratio.sort_values(ascending=False).head(5).items():
    log(f"  {SHORT_NAME[name]:22s} {val:.3f}")

# ---------------------------------------------------------------------------
# Plot 4: "narrow range for success" -- top-performer range vs full range
# ---------------------------------------------------------------------------
combined_score = d["final_malate_mM"].rank(pct=True) + d["log10_ratio"].rank(pct=True)
top_mask = combined_score >= combined_score.quantile(0.75)  # top quartile by combined rank
log(f"\n'Top performer' set for narrow-range analysis: top quartile by combined rank "
    f"(n={top_mask.sum()} of {len(d)})")

narrowness = {}
for c in PARAM_COLS:
    full_range = X_log[c].max() - X_log[c].min()
    top_range = X_log.loc[top_mask, c].max() - X_log.loc[top_mask, c].min()
    narrowness[c] = (top_range / full_range) if full_range > 0 else np.nan
narrow_s = pd.Series(narrowness).sort_values()

fig, ax = plt.subplots(figsize=(9, 9))
colors = [OKABE_ITO[3] if v < 0.35 else (OKABE_ITO[5] if v < 0.7 else OKABE_ITO[2]) for v in narrow_s.values]
ax.barh(range(len(narrow_s)), narrow_s.values, color=colors)
ax.set_yticks(range(len(narrow_s)))
ax.set_yticklabels([SHORT_NAME[c] for c in narrow_s.index])
for i, v in enumerate(narrow_s.values):
    ax.text(max(v, 0.012) + 0.012, i, f"{v:.2f}", va="center", fontsize=8, color="#333333")
ax.invert_yaxis()  # most-constrained (smallest ratio) at the top
ax.axvline(1.0, color="black", linewidth=1, linestyle="--", alpha=0.5)
ax.set_xlim(0, 1.12)
ax.set_xlabel("Top-quartile range / full explored range  (log scale of parameter)")
ax.set_title("Narrow range for success (top) vs loosely deterministic (bottom, near 1)\n"
             "Near-zero: top performers converge on one value for this parameter", fontsize=12)
legend_elems = [
    Line2D([0], [0], marker="s", color="w", markerfacecolor=OKABE_ITO[3], markersize=12, label="Narrow range required (<0.35)"),
    Line2D([0], [0], marker="s", color="w", markerfacecolor=OKABE_ITO[5], markersize=12, label="Moderately constrained (0.35-0.7)"),
    Line2D([0], [0], marker="s", color="w", markerfacecolor=OKABE_ITO[2], markersize=12, label="Loosely deterministic (>=0.7)"),
]
ax.legend(handles=legend_elems, loc="upper right", frameon=False, fontsize=9)
sns.despine(ax=ax)
plt.tight_layout()
plt.savefig(OUT_DIR / "04_narrow_range_for_success.png", dpi=180, bbox_inches="tight")
plt.close()

log("\nMost constrained parameters (narrow range required among top performers):")
for name, val in narrow_s.head(5).items():
    log(f"  {SHORT_NAME[name]:22s} ratio={val:.2f}")
log("Most loosely-deterministic parameters (top performers span nearly the full range):")
for name, val in narrow_s.tail(5).items():
    log(f"  {SHORT_NAME[name]:22s} ratio={val:.2f}")

# ---------------------------------------------------------------------------
# Plot 5: parameter correlation heatmap (+ objectives), Spearman
# ---------------------------------------------------------------------------
corr_input = X_log.copy()
corr_input["Final malate (mM)"] = d["final_malate_mM"].values
corr_input["log10(rTCA:TCA)"] = d["log10_ratio"].values
corr = corr_input.rename(columns=SHORT_NAME).corr(method="spearman")

fig, ax = plt.subplots(figsize=(13, 11))
mask = np.triu(np.ones_like(corr, dtype=bool), k=1)
sns.heatmap(corr, mask=mask, cmap=DIVERGING_CMAP, center=0, vmin=-1, vmax=1,
            square=True, linewidths=0.4, linecolor="white",
            cbar_kws={"label": "Spearman correlation", "shrink": 0.6}, ax=ax)
ax.set_title("Parameter-parameter and parameter-objective correlations\n(negative = tradeoff, positive = co-varying)", fontsize=12)
plt.setp(ax.get_xticklabels(), rotation=45, ha="right")
plt.tight_layout()
plt.savefig(OUT_DIR / "05_correlation_heatmap.png", dpi=180, bbox_inches="tight")
plt.close()

# report strongest off-diagonal parameter-parameter correlations (tradeoffs)
param_only_corr = corr.loc[[SHORT_NAME[c] for c in PARAM_COLS], [SHORT_NAME[c] for c in PARAM_COLS]]
pairs = []
cols = param_only_corr.columns
for i in range(len(cols)):
    for j in range(i + 1, len(cols)):
        pairs.append((cols[i], cols[j], param_only_corr.iloc[i, j]))
pairs_df = pd.DataFrame(pairs, columns=["param_a", "param_b", "spearman_r"]).sort_values("spearman_r")

log("\nStrongest NEGATIVE parameter-parameter correlations (candidate tradeoffs):")
for _, row in pairs_df.head(6).iterrows():
    log(f"  {row['param_a']:22s} vs {row['param_b']:22s} r={row['spearman_r']:+.2f}")
log("Strongest POSITIVE parameter-parameter correlations (co-varying by design):")
for _, row in pairs_df.tail(6).iloc[::-1].iterrows():
    log(f"  {row['param_a']:22s} vs {row['param_b']:22s} r={row['spearman_r']:+.2f}")

log("\nParameter correlation with objectives:")
for obj in ["Final malate (mM)", "log10(rTCA:TCA)"]:
    log(f"  vs {obj}:")
    s = corr.loc[[SHORT_NAME[c] for c in PARAM_COLS], obj].sort_values(key=abs, ascending=False)
    for name, val in s.head(5).items():
        log(f"    {name:22s} r={val:+.2f}")

# ---------------------------------------------------------------------------
# Plot 6: tags vs data-driven clusters
# ---------------------------------------------------------------------------
exploded = d[["cluster", "tags"]].copy()
exploded["tag"] = exploded["tags"].str.split("|")
exploded = exploded.explode("tag").reset_index(drop=True)
crosstab = pd.crosstab(exploded["tag"].to_numpy(), exploded["cluster"].to_numpy())
crosstab_norm = crosstab.div(crosstab.sum(axis=1), axis=0)  # row-normalized: % of that tag's designs in each cluster
crosstab_norm = crosstab_norm.loc[crosstab.sum(axis=1).sort_values(ascending=False).index]

fig, ax = plt.subplots(figsize=(9, 11))
sns.heatmap(crosstab_norm, cmap="Blues", vmin=0, vmax=1, annot=crosstab.loc[crosstab_norm.index],
            fmt="d", cbar_kws={"label": "Fraction of tag's designs in cluster"}, linewidths=0.4, linecolor="white", ax=ax)
ax.set_xlabel("Data-driven cluster")
ax.set_ylabel("Mechanistic tag (assigned during design generation)")
ax.set_title("Do data-driven clusters recover the manual mechanistic tags?\n(cell = count, shade = row-normalized fraction)", fontsize=12)
plt.tight_layout()
plt.savefig(OUT_DIR / "06_tag_vs_cluster.png", dpi=180, bbox_inches="tight")
plt.close()

# Adjusted Rand Index using primary_tag as a single-label proxy for the manual scheme
ari = adjusted_rand_score(d["primary_tag"], d["cluster"])
log(f"\nAdjusted Rand Index (primary tag vs data-driven cluster): {ari:.3f}")
log("(0 = no better than random agreement, 1 = perfect agreement)")

lever_ari = adjusted_rand_score(d["lever_group"].astype(str), d["cluster"])
log(f"Adjusted Rand Index (# levers changed vs data-driven cluster): {lever_ari:.3f}")

# ---------------------------------------------------------------------------
# save summary + augmented CSV
# ---------------------------------------------------------------------------
with open(OUT_DIR / "summary.txt", "w") as f:
    f.write("\n".join(summary_lines))

d_out = d.copy()
d_out.to_csv(OUT_DIR / "designs_with_clusters.csv", index=False)

print(f"\nAll outputs written to {OUT_DIR}")
