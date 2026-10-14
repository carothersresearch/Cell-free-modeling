"""
For each of the 30 pathway designs (model/pathways/<id>.txt, kinetics-fixed
versions), grid-searches dilution factor D (10-400, per CLAUDE.md bounds) and
substrate loading dose (1-100 mM) to maximize PEAK product titer (mM) -- not
yield_pct, since yield_pct = 100*delta/dose has dose in the denominator and
would spuriously reward tiny, impractical doses. Cofactor additions are left
at their originally-designed values (only D and substrate dose are searched,
per the user's request).

Reuses each pathway's already-generated Antimony file (model/pathways/*.txt)
without regeneration -- D is applied via r.dilution_factor after reset(), and
substrate dose is applied post-reset exactly as in run_pathway_sims.py
(native-overlap substrates: add on top of the freshly-computed native
baseline at that D; orthogonal substrates: overwrite directly), so no new
model files are needed for the sweep itself.
"""
import os
import sys
import json
import time

import numpy as np

sys.path.insert(0, os.path.dirname(__file__))
from loader import loada
from pathway_generator import DURATION_S, NATIVE_OVERLAP_SUBSTRATES, ORTHOGONAL_SUBSTRATES, BASE_HET_PARAMS

ROOT = r"C:\Users\mkcoo\Documents\UW\LLM_CFS_2"
UNOPT_RESULTS = os.path.join(ROOT, "results", "pathway_sim_results_unoptimized.json")
OUT_PATH = os.path.join(ROOT, "results", "pathway_optimization.json")

D_GRID = np.geomspace(10, 400, 7)
DOSE_GRID = np.geomspace(1, 100, 7)
N_POINTS_SWEEP = 120     # coarse resolution during grid search, for speed
N_POINTS_FINAL = 400     # high resolution for the confirmed optimum


def evaluate(r, meta, D, dose, n_points):
    for p in BASE_HET_PARAMS:
        setattr(r, p, 0.0)
    r.dilution_factor = 1.0 / D
    r.reset()

    substrate = meta["substrate"]
    if substrate in NATIVE_OVERLAP_SUBSTRATES:
        r[substrate] = float(r[substrate]) + dose
    elif substrate in ORTHOGONAL_SUBSTRATES:
        r[substrate] = dose
    else:
        raise ValueError(f"unknown substrate class for {substrate}")

    for cofactor, cdose in meta["cofactor_additions"].items():
        r[cofactor] = float(r[cofactor]) + cdose

    product = meta["product"]
    product_initial = float(r[product])

    result = r.simulate(0, DURATION_S, n_points)
    colnames = [c.strip("[]") for c in result.colnames]
    idx = colnames.index(product)
    values = [float(result[j][idx]) for j in range(n_points)]
    peak_idx = max(range(n_points), key=lambda j: values[j])
    peak = values[peak_idx]
    peak_time_h = float(result[peak_idx][0]) / 3600.0
    final = values[-1]
    delta_peak = peak - product_initial
    yield_pct_peak = 100.0 * delta_peak / dose if dose else None
    return dict(peak_mM=peak, final_mM=final, initial_mM=product_initial,
                delta_peak_mM=delta_peak, peak_time_h=peak_time_h, yield_pct_peak=yield_pct_peak)


def optimize_one(meta):
    r = loada(meta["file"])
    best = None
    grid_log = []
    n_failed = 0
    for D in D_GRID:
        for dose in DOSE_GRID:
            try:
                res = evaluate(r, meta, D, dose, N_POINTS_SWEEP)
            except RuntimeError:
                # numerically stiff at this (D, dose) extreme -- skip, don't crash the sweep
                n_failed += 1
                grid_log.append(dict(D=float(D), dose=float(dose), peak_mM=None, failed=True))
                r = loada(meta["file"])  # reload -- a failed integration can leave state inconsistent
                continue
            grid_log.append(dict(D=float(D), dose=float(dose), peak_mM=res["peak_mM"]))
            if best is None or res["peak_mM"] > best["peak_mM"]:
                best = dict(D=float(D), dose=float(dose), **res)

    if n_failed:
        print(f"    ({n_failed}/{len(D_GRID)*len(DOSE_GRID)} grid points numerically unstable, skipped)")

    # confirm at high resolution -- if THIS specific point is also unstable at higher
    # resolution, fall back to the coarse-resolution result already in hand.
    try:
        final = evaluate(r, meta, best["D"], best["dose"], N_POINTS_FINAL)
        best.update(final)
    except RuntimeError:
        print(f"    (high-res confirmation at optimum unstable; keeping coarse-resolution values)")
    return best, grid_log


def main():
    with open(UNOPT_RESULTS, encoding="utf-8") as f:
        baseline_results = {r["pid"]: r for r in json.load(f) if "meta" in r}

    out = []
    t0 = time.time()
    for i, (pid, base) in enumerate(baseline_results.items()):
        meta = base["meta"]
        t1 = time.time()
        best, grid_log = optimize_one(meta)
        elapsed = time.time() - t1
        out.append(dict(
            pid=pid, product=meta["product"], substrate=meta["substrate"],
            baseline_D=10.0, baseline_dose_mM=meta["substrate_dose_mM"],
            baseline_peak_mM=base["product_peak_mM"], baseline_yield_pct_peak=base["yield_pct_peak_mol_mol"],
            optimal_D=best["D"], optimal_dose_mM=best["dose"],
            optimal_peak_mM=best["peak_mM"], optimal_yield_pct_peak=best["yield_pct_peak"],
            optimal_peak_time_h=best["peak_time_h"],
            fold_improvement=(best["peak_mM"] / base["product_peak_mM"]) if base["product_peak_mM"] > 1e-9 else None,
        ))
        print(f"[{i+1}/30] {pid}: baseline peak={base['product_peak_mM']:.5f} mM (D=10, dose={meta['substrate_dose_mM']}) "
              f"-> optimal peak={best['peak_mM']:.5f} mM (D={best['D']:.1f}, dose={best['dose']:.1f})  [{elapsed:.1f}s]")

    total = time.time() - t0
    print(f"\nTotal optimization time: {total:.1f}s")

    os.makedirs(os.path.dirname(OUT_PATH), exist_ok=True)
    with open(OUT_PATH, "w", encoding="utf-8") as f:
        json.dump(out, f, indent=2, default=str)
    print(f"Wrote {OUT_PATH}")


if __name__ == "__main__":
    main()
