"""
enzyme_name (as used in pathway_designs.py enzyme tuples) -> EC number
(as used in CellFreeModeling/my_pathways_proteome.csv, fetched via FRENDA).

Many enzyme_names across different pathways are the same catalytic activity
(e.g. ALS, ALS2, ALS3... are all acetolactate synthase) -- they share one EC
and therefore one FRENDA-sourced kinetics record. Where a pathway_designs.py
enzyme tuple used a single generic "dehydrogenase_NADP" template step to stand
in for what is actually a two-step reaction (MCR: malonyl-CoA -> malonate
semialdehyde -> 3-HP), the EC map/generator instead splits it into two real
enzyme steps (MCR_a, MCR_b) to match the two ECs actually fetched from FRENDA.
"""

ENZYME_TO_EC = {
    "FTL": "6.3.4.3",
    "MTHFD": "1.5.1.5",
    "PFL_rev": "2.3.1.54",
    "FDH2": "1.17.1.9",
    "FDH3": "1.17.1.9",
    "AlaDH": "1.4.1.1",
    "AlaDH2": "1.4.1.1",
    "AlaDH3": "1.4.1.1",
    "AlaDH4": "1.4.1.1",
    "PhaB": "1.1.1.36",
    "PhaB2": "1.1.1.36",
    "PhaC": "2.3.1.264",
    "PhaC2": "2.3.1.264",
    "AtoDA": "2.8.3.8",
    "AtoDA2": "2.8.3.8",
    "AdvB": "4.1.1.4",
    "AdvB2": "4.1.1.4",
    "ALDH": "1.2.1.10",   # acylating: AcetylCoA <-> Acetaldehyde + CoA
    "ALDH2": "1.2.1.3",   # non-acylating: Acetaldehyde -> Acetate
    "ALDH3": "1.2.1.3",
    "ADH": "1.1.1.1",
    "ADH2": "1.1.1.1",
    "ADH3": "1.1.1.1",
    "ADH4": "1.1.1.1",
    "ADH5": "1.1.1.1",
    "ALS": "2.2.1.6",
    "ALS2": "2.2.1.6",
    "ALS3": "2.2.1.6",
    "ALS4": "2.2.1.6",
    "ALS5": "2.2.1.6",
    "ALS6": "2.2.1.6",
    "ALDC": "4.1.1.5",
    "ALDC2": "4.1.1.5",
    "ALDC3": "4.1.1.5",
    "ALDC4": "4.1.1.5",
    "ALDC5": "4.1.1.5",
    "BDH": "1.1.1.4",
    "BDH2": "1.1.1.4",
    "BDH3": "1.1.1.4",
    "BDH4": "1.1.1.4",
    "BDH5": "1.1.1.4",
    "KARI": "1.1.1.86",
    "DHAD": "4.2.1.9",
    "KDC": "4.1.1.72",
    "PYC_new": "6.4.1.1",     # collides with base model's hEC6411 -- generator renames
    "PPC": "4.1.1.31",
    "AspDH": "1.4.1.21",
    "PCT": "2.8.3.1",
    "LCD": "4.2.1.54",
    "TE": "3.1.2.20",
    "GlpK": "2.7.1.30",
    "GlpK2": "2.7.1.30",
    "GPD1": "1.1.1.8",
    "GPD1b": "1.1.1.8",
    "MgsA": "4.2.3.3",
    "FucO": "1.1.1.77",
    "GlyDH": "1.1.1.6",
    "XR": "1.1.1.307",
    "XDH": "1.1.1.175",
    "XylA": "5.3.1.5",
    "XylB": "2.7.1.17",
    "MDH_sol": "1.1.1.244",
    "MDH_sol2": "1.1.1.244",
    "MDH_sol3": "1.1.1.244",
    # formaldehyde + THF -> methylene-THF is largely non-enzymatic/spontaneous
    # in vivo (no widely-recognized standalone EC number distinct from the
    # oxidative MTHFD direction, EC 1.5.1.5) -- modeled with a dedicated
    # pseudo-EC key carrying fast literature-typical spontaneous-condensation
    # kinetics (see frenda_kinetics.MANUAL_KINETICS_OVERRIDE).
    "FTHF_cond": "spontaneous_FTHF",
    "Hps": "4.1.2.43",
    "Phi": "5.3.1.27",
    "HK": "2.7.1.1",
    "HK2": "2.7.1.1",
    "MtDH": "1.1.1.67",
    "GDH": "1.1.1.47",
    # MCR is modeled as two real steps in the generator (not a pathway_designs.py
    # enzyme_name -- handled specially for the HCO3_2_3HP design):
    "MCR_a": "1.2.1.75",   # malonyl-CoA reductase (aldehyde-forming)
    "MCR_b": "1.1.1.298",  # malonate-semialdehyde reductase (-> 3-HP)
}

# EC numbers that collide with symbols already declared in model/iMC057.txt
# (base model's engineered malate pathway). The generator appends a suffix to
# any new hEC label built from these ECs to keep the new reaction's enzyme
# pool fully isolated from the original R55/R56 reactions.
COLLIDING_ECS = {"6.4.1.1", "1.17.1.9"}
RESERVED_LABELS = {"eEC11137", "hEC6411", "hEC11719", "hEC43117"}
