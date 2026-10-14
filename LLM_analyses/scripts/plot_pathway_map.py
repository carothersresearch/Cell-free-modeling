"""
One figure summarizing all 30 heterologous pathway designs: substrates on the
left (ordered cheapest-to-priciest, top to bottom), products on the right,
one curved line per pathway connecting them. Dots along each line mark the
enzymatic steps in that pathway -- open circles for the added heterologous
enzymes, closed circles for native/endogenous enzymes the design explicitly
reuses (per NATIVE_STEP_INSERTS below, hand-curated from pathway_designs.py's
own notes text, e.g. "reuses native thiolase (R5)"); line thickness encodes
peak yield; line color encodes the biosynthetic module recruited.
"""
import os
import sys

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

sys.path.insert(0, os.path.dirname(__file__))
from paper_style import MODULE_ORDER, MODULE_COLORS

ROOT = r"C:\Users\mkcoo\Documents\UW\LLM_CFS_2"
OUT_PATH = os.path.join(ROOT, "results", "pathway_analysis", "06_pathway_map.svg")

df = pd.read_csv(os.path.join(ROOT, "results", "pathway_designs_table.csv"))
subs_cost = pd.read_excel(os.path.join(ROOT, "results", "pathway_analysis", "FeedstockCosts_Manual.xlsx"))

SUBSTRATE_RENAME = {"HCO3": "Bicarbonate/CO2", "S_Lactate": "S-Lactate", "D_Xylose": "D-Xylose"}
df["substrate_key"] = df["substrate"].map(lambda s: SUBSTRATE_RENAME.get(s, s))

# FeedstockCosts_Manual.xlsx uses plain substrate names (Bicarbonate, Xylose, Lactate)
# rather than this project's raw substrate codes (HCO3, D_Xylose, S_Lactate) -- map
# xlsx rows onto the raw substrate codes used in pathway_designs_table.csv.
XLSX_SUBSTRATE_TO_RAW = {"Bicarbonate": "HCO3", "Xylose": "D_Xylose", "Lactate": "S_Lactate"}
subs_cost["substrate_raw"] = subs_cost["Substrate"].map(lambda s: XLSX_SUBSTRATE_TO_RAW.get(s, s))
subs_cost = subs_cost.set_index("substrate_raw")
df["substrate_cost"] = df["substrate"].map(subs_cost["$/kg"])

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

# Which designs explicitly reuse a specific named native/endogenous enzyme step
# (parsed by hand from the "reuses native X (Rnn)" mentions in each design's notes
# in pathway_designs.py -- generic "feeds into native glycolysis/PPP" mentions that
# don't name one specific reused enzyme are intentionally excluded). Each value is
# a list of insertion points: the integer is how many of the design's own
# heterologous enzymes come before that native step (0 = native step is first,
# n_enzymes = native step is last); repeat the same integer for consecutive
# native steps at the same point in the chain.
NATIVE_STEP_INSERTS = {
    "Formate_1_Serine": [2],       # native SHMT (R49), after FTL+MTHFD
    "Acetate_1_PHB": [0, 0],       # native acetate activation (R48) + thiolase (R5), before PhaB/PhaC
    "Acetate_2_Acetone": [0],      # native thiolase (R5), before AtoDA
    "Acetate_3_Ethanol": [0],      # native acetate activation (R48), before ALDH/ADH
    "HCO3_1_Malate": [1],          # native MDH, after PYC_new
    "HCO3_2_3HP": [0],             # native acetyl-CoA carboxylase (R6), before MCR
    "Lactate_1_BDO": [0],          # native LDH (R53), before ALS3/ALDC2/BDH2
    "Lactate_2_Alanine": [0],      # native LDH (R53), before AlaDH3
    "Methanol_3_Serine": [2],      # native SHMT (R49), after MDH_sol3+FTHF_cond
    "Glucose_2_Mannitol": [1],     # native PGI (R21), between HK2 and MtDH
    "Ethanol_1_PHB": [2, 2],       # native acetate activation (R48) + thiolase (R5), between ALDH2 and PhaB2
    "Ethanol_3_Acetone": [2],      # native thiolase (R5), between ALDH3 and AtoDA2
}

# ---------------------------------------------------------------- node ordering
substrate_order = (df[["substrate_key", "substrate_cost"]].drop_duplicates()
                   .sort_values("substrate_cost").reset_index(drop=True))
n_sub = len(substrate_order)
substrate_y = {row.substrate_key: i / (n_sub - 1) for i, row in substrate_order.iterrows()}
substrate_order["y"] = substrate_order["substrate_key"].map(substrate_y)

# order products by the (count-weighted) average y of substrates that feed them, to reduce crossings
prod_weight = df.groupby("product_clean").apply(
    lambda g: np.average([substrate_y[s] for s in g["substrate_key"]], weights=None)
)
product_order = prod_weight.sort_values()
n_prod = len(product_order)
product_y = {p: i / (n_prod - 1) for i, p in enumerate(product_order.index)}

df["y0"] = df["substrate_key"].map(substrate_y)
df["y1"] = df["product_clean"].map(product_y)

# ---------------------------------------------------------------- geometry helpers
def smoothstep(t):
    return 3 * t**2 - 2 * t**3


def curve_xy(y0, y1, n=150):
    t = np.linspace(0, 1, n)
    x = t
    y = y0 + (y1 - y0) * smoothstep(t)
    return x, y


MAX_YIELD = df["yield_peak_pct"].max()
WIDTH_MIN, WIDTH_MAX = 0.8, 8.0


def yield_to_width(y):
    y = max(y, 0.0)
    return WIDTH_MIN + (WIDTH_MAX - WIDTH_MIN) * np.sqrt(y) / np.sqrt(MAX_YIELD)


# ---------------------------------------------------------------- plot
fig_h = max(n_prod * 0.275, 4.5)
fig, ax = plt.subplots(figsize=(13, fig_h))
fig.patch.set_facecolor("white")
ax.set_facecolor("white")

for _, row in df.sort_values("yield_peak_pct").iterrows():
    x, y = curve_xy(row["y0"], row["y1"])
    w = yield_to_width(row["yield_peak_pct"])
    color = MODULE_COLORS[row["module_category"]]
    ax.plot(x, y, color=color, linewidth=w, alpha=0.72, solid_capstyle="round", zorder=2)

    n_enz = int(row["n_enzymes"])
    native_inserts = NATIVE_STEP_INSERTS.get(row["design_id"], [])
    insert_counts = {}
    for k in native_inserts:
        insert_counts[k] = insert_counts.get(k, 0) + 1
    is_heterologous = []
    for k in range(n_enz + 1):
        is_heterologous.extend([False] * insert_counts.get(k, 0))
        if k < n_enz:
            is_heterologous.append(True)
    is_heterologous = np.array(is_heterologous, dtype=bool)

    t_dots = np.linspace(0.22, 0.78, len(is_heterologous))
    xd = t_dots
    yd = row["y0"] + (row["y1"] - row["y0"]) * smoothstep(t_dots)
    ax.scatter(xd[is_heterologous], yd[is_heterologous], s=34, facecolor="white",
               edgecolor="#2b2b2b", linewidth=1.0, zorder=4)
    ax.scatter(xd[~is_heterologous], yd[~is_heterologous], s=34, facecolor="#2b2b2b",
               edgecolor="white", linewidth=1.0, zorder=4)

# substrate nodes + labels
for _, row in substrate_order.iterrows():
    ax.scatter([0], [row["y"]], s=170, facecolor="#2b2b2b", edgecolor="white", linewidth=1.2, zorder=5)
    ax.annotate(f"{row['substrate_key']}\n${row['substrate_cost']:.2f}/kg", (0, row["y"]),
                xytext=(-14, 0), textcoords="offset points", ha="right", va="center",
                fontsize=10, fontweight="bold")

# product nodes + labels
for p, y in product_y.items():
    is_native = df.loc[df["product_clean"] == p, "product_is_native"].iloc[0]
    marker = "D" if is_native else "o"
    ax.scatter([1], [y], s=130, facecolor="#2b2b2b", edgecolor="white", linewidth=1.2,
               marker=marker, zorder=5)
    ax.annotate(p, (1, y), xytext=(14, 0), textcoords="offset points",
                ha="left", va="center", fontsize=9.5)

ax.set_xlim(-0.34, 1.30)
ax.set_ylim(-0.06, 1.06)
ax.invert_yaxis()
ax.axis("off")
ax.set_title("30 heterologous pathway designs: substrate -> product\n"
             "(substrates ordered cheapest-to-priciest top-to-bottom; open dots = heterologous enzymes, "
             "closed dots = reused native enzymes;\nline thickness = peak yield; color = biosynthetic module)",
             fontsize=12.5, pad=14)
ax.text(0, -0.055, "SUBSTRATE (cheapest at top)", fontsize=9, style="italic", color="#555555", ha="left")
ax.text(1, -0.055, "PRODUCT", fontsize=9, style="italic", color="#555555", ha="left")

# --- legends
module_handles = [Line2D([0], [0], color=MODULE_COLORS[m], linewidth=4) for m in MODULE_ORDER]
leg1 = ax.legend(module_handles, MODULE_ORDER, title="Biosynthetic module", loc="upper left",
                  bbox_to_anchor=(-0.34, 1.0), fontsize=8.5, title_fontsize=9, frameon=False)
ax.add_artist(leg1)

yield_examples = [0.1, 3, 15]
width_handles = [Line2D([0], [0], color="#555555", linewidth=yield_to_width(v)) for v in yield_examples]
width_labels = [f"{v:g}% peak yield" for v in yield_examples]
leg2 = ax.legend(width_handles, width_labels, title="Line thickness", loc="upper left",
                  bbox_to_anchor=(-0.34, 0.62), fontsize=8.5, title_fontsize=9, frameon=False)
ax.add_artist(leg2)

het_dot_handle = [Line2D([0], [0], marker="o", color="none", markerfacecolor="white",
                          markeredgecolor="#2b2b2b", markersize=7, linewidth=0)]
native_dot_handle = [Line2D([0], [0], marker="o", color="none", markerfacecolor="#2b2b2b",
                             markeredgecolor="white", markersize=7, linewidth=0)]
diamond_handle = [Line2D([0], [0], marker="D", color="none", markerfacecolor="#2b2b2b",
                          markeredgecolor="white", markersize=8, linewidth=0)]
circle_handle = [Line2D([0], [0], marker="o", color="none", markerfacecolor="#2b2b2b",
                         markeredgecolor="white", markersize=10, linewidth=0)]
leg3 = ax.legend(het_dot_handle + native_dot_handle + diamond_handle + circle_handle,
                  ["heterologous enzyme (open)", "native enzyme, reused (closed)",
                   "product already native", "product is novel"],
                  loc="upper left", bbox_to_anchor=(-0.34, 0.40), fontsize=8.5, frameon=False)

plt.tight_layout()
plt.savefig(OUT_PATH, dpi=180, bbox_inches="tight", facecolor="white")
png_path = OUT_PATH.replace(".svg", ".png")
plt.savefig(png_path, dpi=180, bbox_inches="tight", facecolor="white")
plt.close()
print(f"Wrote {OUT_PATH}")
print(f"Wrote {png_path}")
