"""
Merges results/pathway_sim_results.json (simulation output) with
pathway_designs.py (design intent/rationale) into the final deliverable:
results/pathway_designs_table.csv and results/pathway_designs.json.
"""
import os
import sys
import json
import csv

sys.path.insert(0, os.path.dirname(__file__))
from pathway_designs import PATHWAYS

RESULTS_PATH = r"C:\Users\mkcoo\Documents\UW\LLM_CFS_2\results\pathway_sim_results.json"
OUT_CSV = r"C:\Users\mkcoo\Documents\UW\LLM_CFS_2\results\pathway_designs_table.csv"
OUT_JSON = r"C:\Users\mkcoo\Documents\UW\LLM_CFS_2\results\pathway_designs.json"

DESIGN_BY_ID = {p["id"]: p for p in PATHWAYS}


def fmt_enzymes(meta):
    return "; ".join(f"{e['name']} (EC {e['ec']}) = {e['conc_mM']*1000:.1f} uM"
                      for e in meta["enzymes"])


def fmt_new_additions(meta):
    parts = []
    if meta["orthogonal_new_substrate"]:
        parts.append(f"{meta['substrate']} {meta['substrate_dose_mM']:g} mM (new feedstock)")
    else:
        parts.append(f"{meta['substrate']} +{meta['substrate_dose_mM']:g} mM (native pool top-up)")
    for cof, dose in meta["cofactor_additions"].items():
        parts.append(f"{cof} +{dose:g} mM")
    return "; ".join(parts)


def main():
    with open(RESULTS_PATH, encoding="utf-8") as f:
        results = json.load(f)

    rows = []
    for res in results:
        pid = res["pid"]
        design = DESIGN_BY_ID[pid]
        if "error" in res:
            rows.append(dict(
                design_id=pid, substrate=design["substrate"], product=design["product"],
                heterologous_enzymes="ERROR", n_enzymes=len(design["enzymes"]),
                new_substrates_cofactors="", module=design["module"],
                yield_peak_pct=None, yield_final_pct=None, product_peak_mM=None,
                product_final_mM=None, peak_time_h=None,
                msp_usd_per_kg=f"{design['msp'][0]}-{design['msp'][1]}",
                notes=f"SIMULATION FAILED: {res['error']}",
            ))
            continue
        meta = res["meta"]
        rows.append(dict(
            design_id=pid,
            substrate=meta["substrate"],
            product=meta["product"] + (" (native)" if meta["product_native"] else " (novel)"),
            heterologous_enzymes=fmt_enzymes(meta),
            n_enzymes=meta["n_enzymes"],
            new_substrates_cofactors=fmt_new_additions(meta),
            module=meta["module"],
            yield_peak_pct=round(res["yield_pct_peak_mol_mol"], 3) if res["yield_pct_peak_mol_mol"] is not None else None,
            yield_final_pct=round(res["yield_pct_mol_mol"], 3) if res["yield_pct_mol_mol"] is not None else None,
            product_peak_mM=round(res["product_peak_mM"], 5),
            product_final_mM=round(res["product_final_mM"], 5),
            peak_time_h=round(res["peak_time_h"], 2),
            msp_usd_per_kg=f"{meta['msp'][0]}-{meta['msp'][1]}",
            notes=meta["notes"],
        ))

    fieldnames = ["design_id", "substrate", "product", "heterologous_enzymes", "n_enzymes",
                  "new_substrates_cofactors", "module", "yield_peak_pct", "yield_final_pct",
                  "product_peak_mM", "product_final_mM", "peak_time_h", "msp_usd_per_kg", "notes"]
    os.makedirs(os.path.dirname(OUT_CSV), exist_ok=True)
    with open(OUT_CSV, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)

    with open(OUT_JSON, "w", encoding="utf-8") as f:
        json.dump(rows, f, indent=2, default=str)

    print(f"Wrote {OUT_CSV} ({len(rows)} rows)")
    print(f"Wrote {OUT_JSON}")


if __name__ == "__main__":
    main()
