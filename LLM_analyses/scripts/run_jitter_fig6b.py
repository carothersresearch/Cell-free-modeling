"""
N=30 Monte Carlo kinetic-parameter jitter for Figure 6 Panel B (T5 additivity
waterfall). Re-simulates the FIXED 16 conditions of the 2^4 factorial design
(Dilution x CS_KD x Fum_KD x Serine, see scripts/factorial_analysis.py::build_design)
under each jittered draw, and recomputes the same additivity decomposition
(4 single-lever deltas, sum of parts, actual T5 delta, synergy) per draw.
"""
import sys, os, time, json, itertools
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import campaign_jitter_lib as cjl
from sim2 import Design

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
OUT_PATH = os.path.join(ROOT, "results", "jitter_fig6b.json")

N_DRAWS = 30
N_WORKERS = 6
FACTOR_NAMES = ["Dilution", "CS_KD", "Fum_KD", "Serine"]


def build_design(bits):
    dil, cs, fum, ser = bits
    endog = {}
    if cs:
        endog["eEC2331"] = 0.01
    if fum:
        endog["eEC4212"] = 0.01
    exo = {}
    if ser:
        exo["L_Serine"] = 50.0
    return Design(
        dilution_factor_D=10.0 if dil else 200.0,
        endog_multipliers=endog,
        exogenous_overrides_mM=exo,
    )


ALL_BITS = list(itertools.product([0, 1], repeat=4))


def _worker(seed):
    params = cjl.jittered_params(seed)
    malate_by_bits = {}
    for bits in ALL_BITS:
        design = build_design(bits)
        res = cjl.run_design_jittered(design, params, n_points=100)
        malate_by_bits[bits] = res["final_malate_mM"]

    base_bits = (0, 0, 0, 0)
    t5_bits = (1, 1, 1, 1)
    base_malate = malate_by_bits[base_bits]
    t5_malate = malate_by_bits[t5_bits]

    single_deltas = {}
    for i, f in enumerate(FACTOR_NAMES):
        bits = [0, 0, 0, 0]
        bits[i] = 1
        single_deltas[f] = malate_by_bits[tuple(bits)] - base_malate

    sum_singles = sum(single_deltas.values())
    actual_delta = t5_malate - base_malate
    synergy = actual_delta - sum_singles

    return {
        "base_malate": base_malate, "t5_malate": t5_malate,
        "single_deltas": single_deltas,
        "sum_singles": sum_singles, "actual_delta": actual_delta, "synergy": synergy,
    }


if __name__ == "__main__":
    t0 = time.time()
    draws = cjl.run_parallel(_worker, N_DRAWS, n_workers=N_WORKERS, base_seed=9000)
    elapsed = time.time() - t0
    print(f"ran {N_DRAWS} jittered draws x 16 factorial conditions in {elapsed:.1f}s ({N_WORKERS} workers)")

    out = {
        "n_draws": N_DRAWS, "elapsed_s": elapsed,
        "base_malate": [dr["base_malate"] for dr in draws],
        "t5_malate": [dr["t5_malate"] for dr in draws],
        "sum_singles": [dr["sum_singles"] for dr in draws],
        "actual_delta": [dr["actual_delta"] for dr in draws],
        "synergy": [dr["synergy"] for dr in draws],
        "single_deltas": {f: [dr["single_deltas"][f] for dr in draws] for f in FACTOR_NAMES},
    }
    with open(OUT_PATH, "w") as f:
        json.dump(out, f, indent=1)
    print("saved", OUT_PATH)

    print(f"base malate: median={np.median(out['base_malate']):.4f}  IQR=[{np.percentile(out['base_malate'],25):.4f}, {np.percentile(out['base_malate'],75):.4f}]")
    print(f"T5 malate:   median={np.median(out['t5_malate']):.4f}  IQR=[{np.percentile(out['t5_malate'],25):.4f}, {np.percentile(out['t5_malate'],75):.4f}]")
    for f in FACTOR_NAMES:
        vals = np.array(out["single_deltas"][f])
        print(f"  {f:10s} delta: median={np.median(vals):.4f}  IQR=[{np.percentile(vals,25):.4f}, {np.percentile(vals,75):.4f}]")
    syn = np.array(out["synergy"])
    print(f"synergy: median={np.median(syn):.4f}  IQR=[{np.percentile(syn,25):.4f}, {np.percentile(syn,75):.4f}]")
