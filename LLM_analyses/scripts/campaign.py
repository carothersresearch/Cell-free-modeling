"""
100-design exploration campaign for malate/rTCA production.

Implements the propose -> simulate -> accept/reject loop as a single script
(the only practical way to run ~100+ full ODE simulations rigorously rather
than fabricating numbers). For each candidate design, in order:
  - simulate it via scripts/sim2.py's validated Design/run_design
  - accept if it either (a) beats at least one objective of any already-
    accepted design, or (b) introduces at least one mechanistic tag not yet
    represented among accepted designs AND beats the base case on at least
    one objective
  - otherwise reject and move to the next candidate

All designs report values for every parameter touched anywhere in the
campaign (base-case value if that particular design didn't change it), per
the requested CSV format.
"""
import sys
import os
import csv

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from sim2 import Design, run_design

# ---------------------------------------------------------------------------
# Base-case parameter values (every column that appears anywhere below)
# ---------------------------------------------------------------------------
BASE = dict(
    dilution_factor_D=200.0,
    hEC6411_PYC_uM=50.0,
    hEC11719_FDH_uM=25.0,
    hEC43117_SDS_uM=50.0,
    eEC11137_MDH_overexpr_uM=25.0,
    eEC2331_CS_mult=1.0,
    eEC4212_Fumarase_mult=1.0,
    eEC1241_PDH_mult=1.0,
    eEC11142_IDH_mult=1.0,
    eEC41131_PEPsynthase_mult=1.0,
    eEC41149_PEPCK_mult=1.0,
    eEC27140_PyruvateKinase_mult=1.0,
    eEC11127_LDH_mult=1.0,
    ATP_mM=1.0,
    AcetylCoA_mM=1.0,
    HCO3_mM=10.0,
    NADH_mM=1.0,
    NADPH_mM=1.0,
    Formate_mM=10.0,
    L_Serine_mM=1.0,
    Glycine_added_mM=0.0,
    CoA_added_mM=0.0,
    duration_hr=4.0,
)

ENDOG_KEY_TO_SPECIES = {
    "eEC2331_CS_mult": "eEC2331",
    "eEC4212_Fumarase_mult": "eEC4212",
    "eEC1241_PDH_mult": "eEC1241",
    "eEC11142_IDH_mult": "eEC11142",
    "eEC41131_PEPsynthase_mult": "eEC41131",
    "eEC41149_PEPCK_mult": "eEC41149",
    "eEC27140_PyruvateKinase_mult": "eEC27140",
    "eEC11127_LDH_mult": "eEC11127",
}


def make_design(params):
    endog = {}
    for key, sp in ENDOG_KEY_TO_SPECIES.items():
        if params[key] != 1.0:
            endog[sp] = params[key]
    exo = {}
    if params["ATP_mM"] != BASE["ATP_mM"]:
        exo["ATP"] = params["ATP_mM"]
    if params["AcetylCoA_mM"] != BASE["AcetylCoA_mM"]:
        exo["AcetylCoA"] = params["AcetylCoA_mM"]
    if params["HCO3_mM"] != BASE["HCO3_mM"]:
        exo["HCO3"] = params["HCO3_mM"]
    if params["NADH_mM"] != BASE["NADH_mM"]:
        exo["NADH"] = params["NADH_mM"]
    if params["NADPH_mM"] != BASE["NADPH_mM"]:
        exo["NADPH"] = params["NADPH_mM"]
    if params["Formate_mM"] != BASE["Formate_mM"]:
        exo["Formate"] = params["Formate_mM"]
    if params["L_Serine_mM"] != BASE["L_Serine_mM"]:
        exo["L_Serine"] = params["L_Serine_mM"]
    if params["Glycine_added_mM"] > 0:
        exo["Glycine"] = params["Glycine_added_mM"]
    if params["CoA_added_mM"] > 0:
        exo["CoA"] = params["CoA_added_mM"]

    return Design(
        dilution_factor_D=params["dilution_factor_D"],
        duration_sec=params["duration_hr"] * 3600.0,
        het_enzymes_uM={
            "hEC6411": params["hEC6411_PYC_uM"],
            "hEC11719": params["hEC11719_FDH_uM"],
            "hEC43117": params["hEC43117_SDS_uM"],
        },
        eEC11137_overexpr_uM=params["eEC11137_MDH_overexpr_uM"],
        endog_multipliers=endog,
        exogenous_overrides_mM=exo,
    )


def p(**overrides):
    d = dict(BASE)
    d.update(overrides)
    return d


# ---------------------------------------------------------------------------
# Candidate pool: (tags, params)
# ---------------------------------------------------------------------------
candidates = []

def add(tags, **overrides):
    candidates.append((tags if isinstance(tags, list) else [tags], p(**overrides)))

# 1) Pure dilution reduction (yield-focused axis)
for D in [175, 150, 125, 100, 80, 65, 50, 40, 32, 25, 20, 16, 13, 10]:
    add(["dilution_reduction"], dilution_factor_D=float(D))

# 2) Pure dilution increase (selectivity-focused axis)
for D in [250, 300, 350, 400]:
    add(["dilution_increase"], dilution_factor_D=float(D))

# 3) Citrate synthase knockdown alone
for m in [0.5, 0.3, 0.1, 0.05, 0.02, 0.005]:
    add(["citrate_synthase_knockdown"], eEC2331_CS_mult=m)

# 4) Fumarase knockdown alone
for m in [0.5, 0.3, 0.1, 0.05, 0.02, 0.005]:
    add(["fumarase_knockdown"], eEC4212_Fumarase_mult=m)

# 5) Combined CS + fumarase knockdown
for m in [0.3, 0.1, 0.05, 0.02, 0.005]:
    add(["citrate_synthase_knockdown", "fumarase_knockdown"], eEC2331_CS_mult=m, eEC4212_Fumarase_mult=m)

# 6) Serine supplementation alone
for s in [2, 5, 10, 15, 20, 25, 35, 50, 75, 100]:
    add(["serine_supplementation"], L_Serine_mM=float(s))

# 7) HCO3/Formate co-supplementation
for f in [15, 20, 30, 50, 75, 100]:
    add(["hco3_formate_supplementation"], HCO3_mM=float(f), Formate_mM=float(f))

# 8) ATP/AcetylCoA supplementation
for a in [2, 3, 5, 8, 12]:
    add(["atp_acetylcoa_supplementation"], ATP_mM=float(a), AcetylCoA_mM=float(a))

# 9) NADH/NADPH supplementation
for n in [2, 3, 5, 8, 12]:
    add(["nadh_nadph_supplementation"], NADH_mM=float(n), NADPH_mM=float(n))

# 10) Heterologous budget reallocation (150 uM total, held fixed dilution=200)
realloc = [
    (75, 25, 50), (50, 50, 50), (100, 25, 25), (25, 25, 100),
    (25, 100, 25), (110, 15, 25), (25, 15, 110), (60, 30, 60),
    (90, 10, 50), (50, 10, 90),
]
for pyc, fdh, sds in realloc:
    add(["heterologous_budget_reallocation"], hEC6411_PYC_uM=float(pyc), hEC11719_FDH_uM=float(fdh), hEC43117_SDS_uM=float(sds))

# 11) MDH overexpression dose (rebalance PYC/SDS to keep 150 uM total: PYC=SDS=(150-FDH-MDH)/2, FDH=25 fixed)
for mdh in [5, 10, 40, 60, 90, 125]:
    remaining = (150 - 25 - mdh) / 2.0
    add(["mdh_overexpression_dose"], eEC11137_MDH_overexpr_uM=float(mdh), hEC6411_PYC_uM=remaining, hEC43117_SDS_uM=remaining)

# 12) Other competing-enzyme knockdowns (individually)
add(["pdh_knockdown"], eEC1241_PDH_mult=0.1)
add(["pdh_knockdown"], eEC1241_PDH_mult=0.01)
add(["idh_knockdown"], eEC11142_IDH_mult=0.1)
add(["idh_knockdown"], eEC11142_IDH_mult=0.01)
add(["pep_synthase_knockdown"], eEC41131_PEPsynthase_mult=0.1)
add(["pepck_knockdown"], eEC41149_PEPCK_mult=0.1)
add(["pyruvate_kinase_knockdown"], eEC27140_PyruvateKinase_mult=0.1)
add(["ldh_knockdown"], eEC11127_LDH_mult=0.1)
# overexpression variants (opposite direction, exploratory)
add(["pep_synthase_overexpression"], eEC41131_PEPsynthase_mult=5.0)
add(["pepck_overexpression"], eEC41149_PEPCK_mult=5.0)

# 13) Duration changes alone
for h in [0.5, 1, 2, 8, 16, 24]:
    add(["duration_change"], duration_hr=float(h))

# 14) New metabolite additions
for g in [2, 5, 10, 20]:
    add(["glycine_supplementation"], Glycine_added_mM=float(g))
for c in [1, 3, 6]:
    add(["coa_supplementation"], CoA_added_mM=float(c))

# 15) Multi-lever combinations (dilution + genetic KD + serine, various strengths)
combos = [
    (10, 0.5, 0.5, 1),
    (10, 0.1, 0.1, 5),
    (10, 0.1, 0.1, 10),
    (10, 0.1, 0.1, 25),
    (10, 0.05, 0.05, 25),
    (10, 0.01, 0.01, 50),
    (10, 0.01, 0.01, 100),
    (20, 0.1, 0.1, 10),
    (25, 0.1, 0.1, 10),
    (25, 0.01, 0.01, 25),
    (40, 0.1, 0.1, 25),
    (50, 0.05, 0.05, 25),
    (15, 0.1, 0.1, 15),
    (13, 0.05, 0.05, 40),
    (10, 0.1, 0.1, 50),
]
for D, cs, fum, ser in combos:
    add(["dilution_reduction", "citrate_synthase_knockdown", "fumarase_knockdown", "serine_supplementation"],
        dilution_factor_D=float(D), eEC2331_CS_mult=cs, eEC4212_Fumarase_mult=fum, L_Serine_mM=float(ser))

# 16) Dilution + heterologous budget reallocation
add(["dilution_reduction", "heterologous_budget_reallocation"], dilution_factor_D=10.0, hEC6411_PYC_uM=90.0, hEC11719_FDH_uM=10.0, hEC43117_SDS_uM=50.0)
add(["dilution_reduction", "heterologous_budget_reallocation"], dilution_factor_D=10.0, hEC6411_PYC_uM=25.0, hEC11719_FDH_uM=25.0, hEC43117_SDS_uM=100.0)
add(["dilution_reduction", "mdh_overexpression_dose"], dilution_factor_D=10.0, eEC11137_MDH_overexpr_uM=90.0, hEC6411_PYC_uM=17.5, hEC43117_SDS_uM=17.5)
add(["dilution_reduction", "atp_acetylcoa_supplementation"], dilution_factor_D=10.0, ATP_mM=5.0, AcetylCoA_mM=5.0)
add(["dilution_reduction", "hco3_formate_supplementation"], dilution_factor_D=10.0, HCO3_mM=30.0, Formate_mM=30.0)

# 17) Extreme aggressive stacks
add(["extreme_aggressive", "dilution_reduction", "citrate_synthase_knockdown", "fumarase_knockdown",
     "serine_supplementation", "idh_knockdown"],
    dilution_factor_D=10.0, eEC2331_CS_mult=0.005, eEC4212_Fumarase_mult=0.005,
    eEC11142_IDH_mult=0.1, L_Serine_mM=100.0)
add(["extreme_aggressive", "dilution_reduction", "citrate_synthase_knockdown", "fumarase_knockdown",
     "serine_supplementation", "pdh_knockdown"],
    dilution_factor_D=10.0, eEC2331_CS_mult=0.01, eEC4212_Fumarase_mult=0.01,
    eEC1241_PDH_mult=0.1, L_Serine_mM=75.0)
add(["extreme_aggressive", "dilution_reduction", "citrate_synthase_knockdown", "fumarase_knockdown",
     "serine_supplementation", "mdh_overexpression_dose"],
    dilution_factor_D=10.0, eEC2331_CS_mult=0.01, eEC4212_Fumarase_mult=0.01,
    L_Serine_mM=50.0, eEC11137_MDH_overexpr_uM=60.0, hEC6411_PYC_uM=32.5, hEC43117_SDS_uM=32.5)
add(["extreme_aggressive", "dilution_increase", "citrate_synthase_knockdown", "fumarase_knockdown",
     "serine_supplementation"],
    dilution_factor_D=25.0, eEC2331_CS_mult=0.01, eEC4212_Fumarase_mult=0.01, L_Serine_mM=25.0)
add(["extreme_aggressive", "dilution_reduction", "citrate_synthase_knockdown", "fumarase_knockdown",
     "serine_supplementation", "glycine_supplementation"],
    dilution_factor_D=10.0, eEC2331_CS_mult=0.05, eEC4212_Fumarase_mult=0.05, L_Serine_mM=25.0, Glycine_added_mM=10.0)

print(f"Total candidates: {len(candidates)}")

# ---------------------------------------------------------------------------
# Run base case
# ---------------------------------------------------------------------------
base_design = make_design(BASE)
base_res = run_design(base_design)
base_final = base_res["final_malate_mM"]
base_ratio = base_res["rTCA_to_TCA_ratio"]
print(f"BASE: final_malate={base_final:.5f} ratio={base_ratio}")

# ---------------------------------------------------------------------------
# Propose -> simulate -> accept/reject loop
# ---------------------------------------------------------------------------
accepted = []
used_tags = set()

for i, (tags, params) in enumerate(candidates):
    design = make_design(params)
    res = run_design(design)
    final = res["final_malate_mM"]
    ratio = res["rTCA_to_TCA_ratio"]

    novel_mechanism = any(t not in used_tags for t in tags)
    improves_base = (final > base_final) or (ratio > base_ratio)
    better_than_existing = len(accepted) == 0 or any(
        (final > a["final_malate_mM"]) or (ratio > a["rTCA_to_TCA_ratio"]) for a in accepted
    )

    accept = better_than_existing or (novel_mechanism and improves_base)

    print(f"[{i:3d}] tags={tags} final={final:.5f} ratio={ratio:.2f} -> {'ACCEPT' if accept else 'reject'}")

    if accept:
        row = dict(params)
        row["design_id"] = f"D{len(accepted)+1:03d}"
        row["final_malate_mM"] = final
        row["rTCA_to_TCA_ratio"] = ratio
        row["tags"] = "|".join(tags)
        accepted.append(row)
        used_tags.update(tags)

print(f"\nTotal accepted: {len(accepted)}")

# ---------------------------------------------------------------------------
# Write CSV
# ---------------------------------------------------------------------------
columns = ["design_id", "final_malate_mM", "rTCA_to_TCA_ratio"] + list(BASE.keys()) + ["tags"]

out_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "results", "design_campaign.csv")
with open(out_path, "w", newline="") as f:
    writer = csv.DictWriter(f, fieldnames=columns)
    writer.writeheader()
    # also write base case as row 0 for reference
    base_row = dict(BASE)
    base_row["design_id"] = "D000_base"
    base_row["final_malate_mM"] = base_final
    base_row["rTCA_to_TCA_ratio"] = base_ratio
    base_row["tags"] = "base_case"
    writer.writerow(base_row)
    for row in accepted:
        writer.writerow(row)

print(f"Wrote {out_path}")
