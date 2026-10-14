"""
Loads each generated pathway Antimony model (model/pathways/<id>.txt),
simulates it in isolation (original malate-pathway heterologous enzymes
zeroed, D=10, 4 h), and extracts yield metrics for that pathway's product.
"""
import os
import sys
import json

sys.path.insert(0, os.path.dirname(__file__))
from loader import loada
from pathway_generator import (run_all, DILUTION_D, DURATION_S,
                                NATIVE_OVERLAP_SUBSTRATES, BASE_HET_PARAMS,
                                RESULTS_PATH)

N_POINTS = 400


def simulate_one(meta):
    r = loada(meta["file"])
    for p in BASE_HET_PARAMS:
        try:
            setattr(r, p, 0.0)
        except Exception:
            pass
    r.dilution_factor = 1.0 / DILUTION_D
    r.reset()

    product = meta["product"]
    product_initial = float(r[product])

    substrate = meta["substrate"]
    if substrate in NATIVE_OVERLAP_SUBSTRATES:
        r[substrate] = float(r[substrate]) + meta["substrate_dose_mM"]

    for cofactor, dose in meta["cofactor_additions"].items():
        r[cofactor] = float(r[cofactor]) + dose

    result = r.simulate(0, DURATION_S, N_POINTS)
    colnames = [c.strip("[]") for c in result.colnames]
    idx = colnames.index(product)
    times = [result[j][0] for j in range(N_POINTS)]
    values = [float(result[j][idx]) for j in range(N_POINTS)]
    product_final = values[-1]

    # Several designs are cofactor-bolus-limited (a one-time NADH/NADPH
    # addition gets depleted well before 4 h, after which the reversible MRL
    # reaction can run backward and erode the endpoint reading). Peak
    # concentration reached during the run is the more realistic "achievable
    # yield if harvested at the right time" metric -- reported alongside the
    # raw 4 h endpoint for transparency.
    peak_idx = max(range(N_POINTS), key=lambda j: values[j])
    product_peak = values[peak_idx]
    peak_time_h = times[peak_idx] / 3600.0

    delta_final = product_final - product_initial
    delta_peak = product_peak - product_initial
    dose = meta["substrate_dose_mM"]
    yield_pct_final = 100.0 * delta_final / dose if dose else None
    yield_pct_peak = 100.0 * delta_peak / dose if dose else None

    return dict(pid=meta["pid"], product=product, product_initial_mM=product_initial,
                product_final_mM=product_final, delta_mM=delta_final,
                product_peak_mM=product_peak, delta_peak_mM=delta_peak, peak_time_h=peak_time_h,
                substrate_dose_mM=dose, yield_pct_mol_mol=yield_pct_final,
                yield_pct_peak_mol_mol=yield_pct_peak)


def run_all_sims(limit=None):
    all_meta = run_all(limit)
    results = []
    for meta in all_meta:
        try:
            res = simulate_one(meta)
            res["meta"] = meta
            results.append(res)
            print(f"{meta['pid']}: final {res['product']} = {res['product_final_mM']:.4f} mM "
                  f"(yield {res['yield_pct_mol_mol']:.2f}%) | peak = {res['product_peak_mM']:.4f} mM "
                  f"@ {res['peak_time_h']:.2f}h (yield {res['yield_pct_peak_mol_mol']:.2f}%)")
        except Exception as e:
            print(f"{meta['pid']}: FAILED -- {type(e).__name__}: {e}")
            results.append(dict(pid=meta["pid"], error=str(e), meta=meta))

    os.makedirs(os.path.dirname(RESULTS_PATH), exist_ok=True)
    with open(RESULTS_PATH, "w", encoding="utf-8") as f:
        json.dump(results, f, indent=2, default=str)
    print(f"\nWrote {RESULTS_PATH}")
    return results


if __name__ == "__main__":
    limit = int(sys.argv[1]) if len(sys.argv) > 1 else None
    run_all_sims(limit)
