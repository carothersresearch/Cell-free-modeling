"""
Titer-range study: for each of the 30 FRENDA-parameterized pathway designs,
sweep the FULL dilution (D) x substrate-dose grid (same D_GRID/DOSE_GRID as
optimize_pathways.py: 7x7=49 points, D in [10,400] log-spaced, dose in
[1,100] mM log-spaced) and record peak product titer at EVERY grid point,
not just the optimum. This answers "what's the range of achievable titers
for this pathway" rather than "what's the best achievable titer."

Reuses optimize_one() from optimize_pathways.py unchanged -- it already
computes a full grid_log internally during its search, it just wasn't being
persisted. This script captures that grid_log for every pathway.
"""
import os
import sys
import json
import time

sys.path.insert(0, os.path.dirname(__file__))
from optimize_pathways import optimize_one, D_GRID, DOSE_GRID

ROOT = r"C:\Users\mkcoo\Documents\UW\LLM_CFS_2"
UNOPT_RESULTS = os.path.join(ROOT, "results", "pathway_sim_results_unoptimized.json")
OUT_PATH = os.path.join(ROOT, "results", "pathway_titer_grid.json")


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
            best=best, grid=grid_log,
        ))
        n_valid = sum(1 for g in grid_log if not g.get("failed"))
        print(f"[{i+1}/30] {pid}: {n_valid}/{len(grid_log)} valid grid points, "
              f"best peak={best['peak_mM']:.5g} mM at D={best['D']:.1f}, dose={best['dose']:.1f}  [{elapsed:.1f}s]")

    total = time.time() - t0
    print(f"\nTotal sweep time: {total:.1f}s")

    os.makedirs(os.path.dirname(OUT_PATH), exist_ok=True)
    with open(OUT_PATH, "w", encoding="utf-8") as f:
        json.dump(dict(D_GRID=[float(d) for d in D_GRID], DOSE_GRID=[float(d) for d in DOSE_GRID], results=out), f, indent=2, default=str)
    print(f"Wrote {OUT_PATH}")


if __name__ == "__main__":
    main()
