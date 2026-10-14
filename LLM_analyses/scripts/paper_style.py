"""
Shared style constants for the combined publication figures (fig_malate_campaign.py,
fig_pathway_diversification.py), reverse-engineered from CellFree_Fig4_v3.png by
sampling its legend-swatch pixels directly (scripts/paper_style.py itself does not
regenerate that sampling -- see chat history / git blame if the source values below
ever need re-deriving).

Everything downstream (module colors, panel-label style, spine style) is centralized
here so the two figure scripts, and any future panel added to either, stay visually
consistent with each other and with CellFree_Fig4_v3.png.
"""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# Fig4 uses a full-box (4-spine) axis style throughout, not the despined look used
# elsewhere in this project's exploratory analysis plots -- match it here.
plt.rcParams.update({
    "font.family": "Arial",
    "font.size": 10,
    "axes.linewidth": 1.1,
    "axes.edgecolor": "#000000",
    "axes.spines.top": True,
    "axes.spines.right": True,
    "xtick.major.width": 1.0,
    "ytick.major.width": 1.0,
    "axes.grid": False,
    "figure.facecolor": "white",
    "savefig.facecolor": "white",
})

PANEL_LABEL_KW = dict(fontsize=18, fontweight="bold", ha="left", va="top")

# --- colors sampled directly from CellFree_Fig4_v3.png ---------------------------
FIG4_MODULE_COLORS = {
    "TCA Cycle": "#2A7F62",
    "Glycolysis/Gluconeogenesis": "#8390FA",
    "Pentose Phosphate Pathway": "#FAC748",
    "Fatty Acid Biosynthesis": "#F06C9B",
    "Fatty Acid Degradation": "#B6244F",
    "Pyruvate Metabolism": "#FF9913",
    "Glycine & Serine Metabolism": "#C8A2C8",
    "Regeneration Reactions": "#7F00FF",
}
FIG4_RED = "#CC4C44"    # panel D "increase" bars
FIG4_BLUE = "#65A9CF"   # panel D "decrease" bars
FIG4_GRAY = "#999999"   # neutral / "Other" / half-max reference line

# Maps this project's pathway_analysis module_category vocabulary (see
# plot_pathway_map.py / plot_titer_ranges.py) onto the Fig4 palette above, so the
# two new pathway-diversification panels read as part of the same figure family as
# the original malate/rTCA native-module analysis (CellFree_Fig4_v3.png panels A-C).
MODULE_ORDER = ["C1/serine metabolism", "TCA/anaplerosis", "Fatty acid/ketone metabolism",
                 "Amino acid biosynthesis", "Pentose phosphate pathway", "Glycolysis",
                 "Pyruvate metabolism", "Other"]
MODULE_COLORS = {
    "C1/serine metabolism": FIG4_MODULE_COLORS["Glycine & Serine Metabolism"],
    "TCA/anaplerosis": FIG4_MODULE_COLORS["TCA Cycle"],
    "Fatty acid/ketone metabolism": FIG4_MODULE_COLORS["Fatty Acid Biosynthesis"],
    "Amino acid biosynthesis": FIG4_MODULE_COLORS["Regeneration Reactions"],
    "Pentose phosphate pathway": FIG4_MODULE_COLORS["Pentose Phosphate Pathway"],
    "Glycolysis": FIG4_MODULE_COLORS["Glycolysis/Gluconeogenesis"],
    "Pyruvate metabolism": FIG4_MODULE_COLORS["Pyruvate Metabolism"],
    "Other": FIG4_GRAY,
}


def add_panel_label(fig, ax, letter):
    """Place a bold 'A.'-style panel label just outside the axes' top-left corner,
    matching CellFree_Fig4_v3.png's panel-label placement/weight."""
    bbox = ax.get_position()
    fig.text(bbox.x0 - 0.018, bbox.y1 + 0.012, letter, **PANEL_LABEL_KW)
