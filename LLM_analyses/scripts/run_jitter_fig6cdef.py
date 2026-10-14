"""
N=30 Monte Carlo kinetic-parameter jitter for Figure 6 Panels C/D/E/F (citrate
synthase / dilution synergy). Re-simulates, under each jittered draw:
  - the 12-point dilution scan (D in {10,15,20,30,50,75,100,125,150,200,300,400})
    x {CS native, CS knockdown} used by panels C/D/E (CS synthesis flux vs D,
    % malate gain from CS knockdown vs D, final malate native-vs-KD at D=200/D=10)
  - the 4 full AcetylCoA time courses (D in {200,10} x {native,KD}) used by panel F
matching scripts/citrate_synthase_synergy_analysis.py's exact conditions.
"""
import sys, os, time, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import campaign_jitter_lib as cjl
from sim2 import Design

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
OUT_PATH = os.path.join(ROOT, "results", "jitter_fig6cdef.json")

N_DRAWS = 30
N_WORKERS = 6

D_SCAN = [10, 15, 20, 30, 50, 75, 100, 125, 150, 200, 300, 400]
CS_NATIVE = {}
CS_KD = {"eEC2331": 0.01}
BASE_KW = dict(
    het_enzymes_uM={"hEC6411": 50.0, "hEC11719": 25.0, "hEC43117": 50.0},
    eEC11137_overexpr_uM=25.0,
    duration_sec=14400.0,
)
TC_DILFACS = [200.0, 10.0]
N_POINTS_SCAN = 200
N_POINTS_TC = 300


def _worker(seed):
    params = cjl.jittered_params(seed)

    # --- C/D/E: dilution scan, mean CS synthesis flux + final malate, native vs KD ---
    scan = {}
    for D in D_SCAN:
        row = {}
        for tag, endog in [("wt", CS_NATIVE), ("kd", CS_KD)]:
            design = Design(dilution_factor_D=D, endog_multipliers=dict(endog), **BASE_KW)
            res = cjl.run_design_jittered(design, params, n_points=N_POINTS_SCAN,
                                           extra_selections=["R46_f", "R46_r"])
            t = res["time"]
            sink = res["R46_r"] - res["R46_f"]
            mean_sink = float(np.trapz(sink, t) / (t[-1] - t[0]))
            row[f"mean_CS_sink_flux_{tag}"] = mean_sink
            row[f"final_malate_{tag}"] = res["final_malate_mM"]
        row["cs_kd_delta_malate"] = row["final_malate_kd"] - row["final_malate_wt"]
        denom = row["final_malate_wt"]
        row["cs_kd_pct_improvement"] = 100.0 * row["cs_kd_delta_malate"] / denom if denom else float("nan")
        scan[D] = row

    # --- F: full AcetylCoA time courses, 4 conditions ---
    timecourses = {}
    for D in TC_DILFACS:
        for tag, endog in [("native", CS_NATIVE), ("kd", CS_KD)]:
            design = Design(dilution_factor_D=D, endog_multipliers=dict(endog), **BASE_KW)
            res = cjl.run_design_jittered(design, params, n_points=N_POINTS_TC,
                                           extra_selections=["AcetylCoA"])
            key = f"D{int(D)}_{tag}"
            timecourses[key] = {"time_h": (res["time"] / 3600.0).tolist(), "AcetylCoA_mM": res["AcetylCoA"].tolist()}

    return {"scan": scan, "timecourses": timecourses}


if __name__ == "__main__":
    t0 = time.time()
    draws = cjl.run_parallel(_worker, N_DRAWS, n_workers=N_WORKERS, base_seed=10000)
    elapsed = time.time() - t0
    print(f"ran {N_DRAWS} jittered draws x (12 D-scan pts x 2 CS states + 4 timecourses) in {elapsed:.1f}s ({N_WORKERS} workers)")

    scan_ensemble = {str(D): {k: [dr["scan"][D][k] for dr in draws] for k in draws[0]["scan"][D]} for D in D_SCAN}
    tc_ensemble = {key: [dr["timecourses"][key] for dr in draws] for key in draws[0]["timecourses"]}

    out = {"n_draws": N_DRAWS, "elapsed_s": elapsed, "D_SCAN": D_SCAN,
           "scan_ensemble": scan_ensemble, "timecourse_ensemble": tc_ensemble}
    with open(OUT_PATH, "w") as f:
        json.dump(out, f, indent=1)
    print("saved", OUT_PATH)

    for D in D_SCAN:
        vals = np.array(scan_ensemble[str(D)]["cs_kd_pct_improvement"])
        print(f"  D={D:4d}  cs_kd_pct_improvement median={np.nanmedian(vals):.2f}%  IQR=[{np.nanpercentile(vals,25):.2f}, {np.nanpercentile(vals,75):.2f}]")
