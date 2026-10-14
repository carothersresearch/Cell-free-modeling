"""
Generates 30 new Antimony model files (one per pathway design in
pathway_designs.PATHWAYS), using REAL FRENDA/BRENDA/CatPred-sourced kinetics
(scripts/frenda_kinetics.py) applied to the hand-designed reaction topology
in pathway_designs.py, spliced into a copy of the base crude-lysate model
(model/iMC057.txt). Then simulates each and extracts yield metrics.

Rate law: identical Common Modular Rate Law (MRL) macro-expansion style used
throughout model/iMC057.txt (see R25/R55 for the reference pattern this
mirrors exactly).
"""
import os
import re
import sys
import json

sys.path.insert(0, os.path.dirname(__file__))
sys.path.insert(0, r"C:\Users\mkcoo\Documents\UW\LLM_CFS_2")

from pathway_designs import PATHWAYS
from ec_map import ENZYME_TO_EC, RESERVED_LABELS
from frenda_kinetics import load_kinetics, get_km_for_role
from loader import loada

BASE_MODEL_PATH = r"C:\Users\mkcoo\Documents\UW\LLM_CFS_2\model\iMC057.txt"
OUT_DIR = r"C:\Users\mkcoo\Documents\UW\LLM_CFS_2\model\pathways"
RESULTS_PATH = r"C:\Users\mkcoo\Documents\UW\LLM_CFS_2\results\pathway_sim_results.json"

TOTAL_BUDGET_MM = 0.150   # 150 uM total heterologous expression budget per design
DILUTION_D = 10.0         # least-dilute bound; matches original malate pathway convention
DURATION_S = 4 * 60 * 60  # 4 h

NATIVE_OVERLAP_SUBSTRATES = {"Formate", "Acetate", "Pyruvate", "HCO3", "S_Lactate"}
ORTHOGONAL_SUBSTRATES = {"Glycerol", "D_Xylose", "Methanol", "Glucose", "Ethanol"}

BASE_HET_PARAMS = ["p_hEC6411", "p_hEC11719", "p_hEC43117"]  # zeroed for isolation


def get_base_species(base_text):
    return set(re.findall(r'^species\s+([A-Za-z_][A-Za-z0-9_]*)\s*;', base_text, re.MULTILINE))


def ec_label(ec):
    return "hEC" + ec.replace(".", "")


def unique_label(ec, used_labels):
    lbl = ec_label(ec)
    if lbl in RESERVED_LABELS or lbl in used_labels:
        lbl = lbl + "X"
    used_labels.add(lbl)
    return lbl


def mrl_block(rid, enzyme_label, subs, prods, kcat_v, keq):
    """subs/prods: list of (species_name, km_value, stoich_coeff). Returns (equation_text,
    param_lines). Stoichiometric coefficient generalizes the standard Common Modular Rate
    Law construction (Liebermeister & Klipp): a co-substrate/co-product present at
    stoichiometry n contributes (1+S/Km)^n to the saturation term and S^n/Km^n to the
    mass-action term, rather than the implicit n=1 used when every reaction was assumed
    to consume/produce exactly one copy of each listed species."""
    sub_names = [s[0] for s in subs]
    prod_names = [p[0] for p in prods]
    sub_coeffs = [s[2] for s in subs]
    prod_coeffs = [p[2] for p in prods]
    sub_kms = [f"Km_{n}_{enzyme_label}" for n in sub_names]
    prod_kms = [f"Km_{n}_{enzyme_label}" for n in prod_names]

    lhs = " + ".join(f"{c} {n}" for n, c in zip(sub_names, sub_coeffs))
    rhs = " + ".join(f"{c} {n}" for n, c in zip(prod_names, prod_coeffs))
    eq = f"{rid} : {lhs} + {enzyme_label} -> {enzyme_label} + {rhs};"

    D_sub = "*".join(f"((1+{n}/{k})^{c})" for n, k, c in zip(sub_names, sub_kms, sub_coeffs))
    D_prod = "*".join(f"((1+{n}/{k})^{c})" for n, k, c in zip(prod_names, prod_kms, prod_coeffs))
    D_line = f"{rid}_D := ({D_sub} + {D_prod} -1);"

    sub_km_prod = "*".join(f"({k}^{c})" for k, c in zip(sub_kms, sub_coeffs))
    prod_km_prod = "*".join(f"({k}^{c})" for k, c in zip(prod_kms, prod_coeffs))
    hf = f"hKcat_F_{rid} := Kcat_V_{rid}*(Keq_{rid}*({sub_km_prod})/({prod_km_prod}))^(0.5);"
    hr = f"hKcat_R_{rid} := Kcat_V_{rid}*(Keq_{rid}*({sub_km_prod})/({prod_km_prod}))^(-0.5);"

    sub_conc_prod = "*".join(f"({n}^{c})" for n, c in zip(sub_names, sub_coeffs))
    prod_conc_prod = "*".join(f"({n}^{c})" for n, c in zip(prod_names, prod_coeffs))
    f_line = f"{rid}_f := {enzyme_label} * (hKcat_F_{rid}*({sub_conc_prod})/({sub_km_prod}))/{rid}_D;"
    r_line = f"{rid}_r := {enzyme_label} * (hKcat_R_{rid}*({prod_conc_prod})/({prod_km_prod}))/{rid}_D;"
    rate_line = f"{rid} = {rid}_f - {rid}_r;"

    eq_block = "\n".join([eq, D_line, hf, hr, f_line, r_line, rate_line])

    params = []
    for k, (n, v, c) in zip(sub_kms, subs):
        params.append(f"{k} = {v};")
    for k, (n, v, c) in zip(prod_kms, prods):
        params.append(f"{k} = {v};")
    params.append(f"Kcat_V_{rid} = {kcat_v};")
    params.append(f"Keq_{rid} = {keq};")
    return eq_block, params


def get_enzyme_steps(pathway, kin_lookup, used_labels):
    """Returns list of dicts: enzyme_label, ec, rid, subs:[(name,km,coeff)], prods:[(name,km,coeff)].

    Substrate-list entries are normally 2-tuples (species_name, km_role_key), implicitly
    stoichiometry 1. A 3rd element gives an explicit stoichiometric coefficient (e.g.
    ("Pyruvate", "default", 2) for a co-substrate present at stoichiometry 2) -- used for
    acetolactate synthase (ALS), the only enzyme in this library whose real mechanism
    consumes two copies of the same substrate (2 Pyruvate -> Acetolactate + CO2) rather
    than one copy each of two distinct substrates; see frenda_kinetics.FORCED_OVERRIDE's
    own comment for the real reaction this stoichiometry corrects. Product lists are plain
    name strings and are always stoichiometry 1 (no product in this library needs
    otherwise)."""
    pid = pathway["id"]
    steps_raw = pathway["enzymes"]

    if pid == "HCO3_2_3HP":
        steps_raw = [
            ("MCR_a", None, "N1", [("Malonyl_CoA", None), ("NADPH", None)], ["MalonateSemialdehyde", "CoA", "NADP"]),
            ("MCR_b", None, "N2", [("MalonateSemialdehyde", None), ("NADPH", None)], ["ThreeHP", "NADP"]),
        ]

    steps = []
    for enzyme_name, _template, rid, sub_list, prod_names in steps_raw:
        ec = ENZYME_TO_EC[enzyme_name]
        kin = kin_lookup.get(ec)
        if kin is None:
            raise KeyError(f"No FRENDA kinetics found for EC {ec} ({enzyme_name}) in pathway {pid}")
        label = unique_label(ec, used_labels)
        sub_names = [s[0] for s in sub_list]
        sub_coeffs = [s[2] if len(s) > 2 else 1 for s in sub_list]
        subs = [(n, get_km_for_role(kin, "km_subs", i), c) for i, (n, c) in enumerate(zip(sub_names, sub_coeffs))]
        prods = [(n, get_km_for_role(kin, "km_prods", i), 1) for i, n in enumerate(prod_names)]
        steps.append(dict(enzyme_name=enzyme_name, enzyme_label=label, ec=ec, rid=f"{pid}_{rid}",
                           subs=subs, prods=prods, kcat_v=kin["kcat_v"], keq=kin["keq"]))
    return steps


def build_pathway_model(pathway, base_text, base_species, kin_lookup):
    pid = pathway["id"]
    used_labels = set()
    steps = get_enzyme_steps(pathway, kin_lookup, used_labels)

    n_enz = len(steps)
    budget_share = TOTAL_BUDGET_MM / n_enz

    declared_new = set()
    species_decls = []
    param_lines = []
    reaction_blocks = []
    degradation_blocks = []
    enzyme_labels_conc = {}

    # enzyme species
    for st in steps:
        lbl = st["enzyme_label"]
        if lbl not in declared_new:
            species_decls.append(f"species {lbl};\n{lbl} = {budget_share};")
            param_lines.append(f"kdeg_{lbl} = 0.0001;")
            degradation_blocks.append(f"D_{lbl} : {lbl} -> ;\nD_{lbl} = kdeg_{lbl}*{lbl};")
            declared_new.add(lbl)
            enzyme_labels_conc[lbl] = budget_share

    # substrate/product species not already native
    orthogonal_sub = pathway["substrate"] if pathway["substrate"] in ORTHOGONAL_SUBSTRATES else None
    for st in steps:
        for n, _km, _coeff in st["subs"] + st["prods"]:
            if n in base_species or n in declared_new:
                continue
            if n == orthogonal_sub:
                conc = pathway["substrate_dose_mM"]
            else:
                conc = 0
            species_decls.append(f"species {n};\n{n} = {conc};")
            declared_new.add(n)

    # reaction blocks
    for st in steps:
        eq_block, plines = mrl_block(st["rid"], st["enzyme_label"], st["subs"], st["prods"],
                                      st["kcat_v"], st["keq"])
        reaction_blocks.append(eq_block)
        param_lines.extend(plines)

    new_text = "\n\n".join(
        ["\n".join(species_decls), "\n".join(param_lines),
         "\n\n".join(reaction_blocks), "\n\n".join(degradation_blocks)]
    )

    full_text = base_text + "\n\n// ===== Heterologous pathway: " + pid + " =====\n\n" + new_text + "\n"

    meta = dict(pid=pid, product=pathway["product"], substrate=pathway["substrate"],
                substrate_dose_mM=pathway["substrate_dose_mM"],
                orthogonal_new_substrate=orthogonal_sub is not None,
                enzymes=[dict(name=s["enzyme_name"], label=s["enzyme_label"], ec=s["ec"],
                               conc_mM=enzyme_labels_conc[s["enzyme_label"]]) for s in steps],
                n_enzymes=n_enz, cofactor_additions=pathway["cofactor_additions"],
                module=pathway["module"], notes=pathway["notes"], msp=pathway["msp"],
                product_native=pathway["product_native"])
    return full_text, meta


def run_all(limit=None):
    os.makedirs(OUT_DIR, exist_ok=True)
    with open(BASE_MODEL_PATH, encoding="utf-8") as f:
        base_text = f.read()
    base_species = get_base_species(base_text)
    kin_lookup = load_kinetics()

    pathways = PATHWAYS if limit is None else PATHWAYS[:limit]
    all_meta = []
    for pathway in pathways:
        full_text, meta = build_pathway_model(pathway, base_text, base_species, kin_lookup)
        out_path = os.path.join(OUT_DIR, f"{pathway['id']}.txt")
        with open(out_path, "w", encoding="utf-8") as f:
            f.write(full_text)
        meta["file"] = out_path
        all_meta.append(meta)
        print(f"Wrote {out_path} ({meta['n_enzymes']} enzymes)")
    return all_meta


if __name__ == "__main__":
    limit = int(sys.argv[1]) if len(sys.argv) > 1 else None
    meta = run_all(limit)
    print(f"\nGenerated {len(meta)} pathway models.")
