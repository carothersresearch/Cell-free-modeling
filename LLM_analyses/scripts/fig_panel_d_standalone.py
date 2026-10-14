"""
Standalone export of Figure 6 Panel D (final malate vs. dilution D, base vs. CS
knockdown, single representative parameter set, no jitter). Reuses panel_d()
from fig_malate_campaign.py so the two stay in sync -- edit panel_d() there, not
here.

Usage:
    ./.venv_sim/Scripts/python.exe scripts/fig_panel_d_standalone.py
"""
import os
import sys

import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(__file__))
from fig_malate_campaign import panel_d

ROOT = r"C:\Users\mkcoo\Documents\UW\LLM_CFS_2"
OUT_PATH = os.path.join(ROOT, "results", "analysis", "Fig_panel_D_standalone.svg")

if __name__ == "__main__":
    fig, ax = plt.subplots(figsize=(3.5, 3.2))
    panel_d(ax)
    plt.tight_layout()
    plt.savefig(OUT_PATH)
    print("saved", OUT_PATH)
