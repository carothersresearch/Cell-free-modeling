"""
N=30 Monte Carlo kinetic-parameter jitter for Figure 6 Panel A (RF feature
importance). Re-simulates the FIXED 100 already-accepted design-campaign
conditions (results/design_campaign_final.csv, base case excluded) under each
jittered draw, refits the same two RandomForestRegressor models (final malate,
log10 rTCA:TCA ratio) per draw, and saves the per-draw feature importances so
Panel A can plot median+IQR instead of a single point estimate.

The design SET itself (which 100 conditions were accepted) is not re-decided
per draw -- that selection was already made from a much larger candidate pool
using the point-estimate simulation and is treated as fixed; only the
downstream simulated outcome (and therefore the RF fit) varies with jitter.
"""
import sys, os, time, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestRegressor
import campaign_jitter_lib as cjl
from sim2 import Design

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
CAMPAIGN_CSV = os.path.join(ROOT, "results", "design_campaign_final.csv")
OUT_PATH = os.path.join(ROOT, "results", "jitter_fig6a.json")

N_DRAWS = 30
N_WORKERS = 6

PARAM_COLS = [
    "dilution_factor_D", "hEC6411_PYC_uM", "hEC11719_FDH_uM", "hEC43117_SDS_uM",
    "eEC11137_MDH_overexpr_uM", "eEC2331_CS_mult", "eEC4212_Fumarase_mult",
    "eEC1241_PDH_mult", "eEC11142_IDH_mult", "eEC41131_PEPsynthase_mult",
    "eEC41149_PEPCK_mult", "eEC27140_PyruvateKinase_mult", "eEC11127_LDH_mult",
    "ATP_mM", "AcetylCoA_mM", "HCO3_mM", "NADH_mM", "NADPH_mM", "Formate_mM",
    "L_Serine_mM", "Glycine_added_mM", "CoA_added_mM", "duration_hr",
]

ENDOG_KEY_TO_SPECIES = {
    "eEC2331_CS_mult": "eEC2331", "eEC4212_Fumarase_mult": "eEC4212",
    "eEC1241_PDH_mult": "eEC1241", "eEC11142_IDH_mult": "eEC11142",
    "eEC41131_PEPsynthase_mult": "eEC41131", "eEC41149_PEPCK_mult": "eEC41149",
    "eEC27140_PyruvateKinase_mult": "eEC27140", "eEC11127_LDH_mult": "eEC11127",
}


def row_to_design(row):
    endog = {sp: float(row[key]) for key, sp in ENDOG_KEY_TO_SPECIES.items()}
    exo = {
        "ATP": float(row["ATP_mM"]), "AcetylCoA": float(row["AcetylCoA_mM"]),
        "HCO3": float(row["HCO3_mM"]), "NADH": float(row["NADH_mM"]),
        "NADPH": float(row["NADPH_mM"]), "Formate": float(row["Formate_mM"]),
        "L_Serine": float(row["L_Serine_mM"]),
    }
    if float(row["Glycine_added_mM"]) > 0:
        exo["Glycine"] = float(row["Glycine_added_mM"])
    if float(row["CoA_added_mM"]) > 0:
        exo["CoA"] = float(row["CoA_added_mM"])
    return Design(
        dilution_factor_D=float(row["dilution_factor_D"]),
        duration_sec=float(row["duration_hr"]) * 3600.0,
        het_enzymes_uM={"hEC6411": float(row["hEC6411_PYC_uM"]), "hEC11719": float(row["hEC11719_FDH_uM"]),
                         "hEC43117": float(row["hEC43117_SDS_uM"])},
        eEC11137_overexpr_uM=float(row["eEC11137_MDH_overexpr_uM"]),
        endog_multipliers=endog,
        exogenous_overrides_mM=exo,
    )


def _worker(seed):
    df = pd.read_csv(CAMPAIGN_CSV)
    base_mask = df["design_id"].str.contains("base")
    d = df[~base_mask].reset_index(drop=True).copy()
    designs = [row_to_design(row) for _, row in d.iterrows()]

    params = cjl.jittered_params(seed)
    final_malate = []
    ratio = []
    for design in designs:
        res = cjl.run_design_jittered(design, params)
        final_malate.append(res["final_malate_mM"])
        ratio.append(res["rTCA_to_TCA_ratio"])
    d["final_malate_mM"] = final_malate
    d["rTCA_to_TCA_ratio"] = ratio
    # drop designs where this jittered draw failed/was implausible (NaN final malate) --
    # RandomForestRegressor can't fit against a NaN target; the remaining ~85-100 of
    # 100 designs are still a fine sample for one draw's importance estimate.
    d = d[np.isfinite(d["final_malate_mM"])].reset_index(drop=True)

    finite_ratios = d.loc[np.isfinite(d["rTCA_to_TCA_ratio"]), "rTCA_to_TCA_ratio"]
    cap = finite_ratios.max() * 10 if len(finite_ratios) else 1e6
    d["log10_ratio"] = np.log10(d["rTCA_to_TCA_ratio"].replace([np.inf], cap) + 1)

    X_log = pd.DataFrame(index=d.index)
    for c in PARAM_COLS:
        vals = d[c].values.astype(float)
        X_log[c] = np.log1p(vals) if (vals <= 0).any() else np.log(vals)
    X_scaled = (X_log - X_log.mean()) / X_log.std()

    rf_malate = RandomForestRegressor(n_estimators=500, random_state=42, max_depth=6)
    rf_malate.fit(X_scaled, d["final_malate_mM"])
    imp_malate = dict(zip(PARAM_COLS, rf_malate.feature_importances_.tolist()))

    rf_ratio = RandomForestRegressor(n_estimators=500, random_state=42, max_depth=6)
    rf_ratio.fit(X_scaled, d["log10_ratio"])
    imp_ratio = dict(zip(PARAM_COLS, rf_ratio.feature_importances_.tolist()))

    return {"imp_malate": imp_malate, "imp_ratio": imp_ratio}


if __name__ == "__main__":
    t0 = time.time()
    draws = cjl.run_parallel(_worker, N_DRAWS, n_workers=N_WORKERS, base_seed=8000)
    elapsed = time.time() - t0
    print(f"ran {N_DRAWS} jittered draws x 100 designs x 2 RF fits in {elapsed:.1f}s ({N_WORKERS} workers)")

    imp_malate_ens = {c: [dr["imp_malate"][c] for dr in draws] for c in PARAM_COLS}
    imp_ratio_ens = {c: [dr["imp_ratio"][c] for dr in draws] for c in PARAM_COLS}

    out = {"n_draws": N_DRAWS, "elapsed_s": elapsed, "param_cols": PARAM_COLS,
           "imp_malate_ensemble": imp_malate_ens, "imp_ratio_ensemble": imp_ratio_ens}
    with open(OUT_PATH, "w") as f:
        json.dump(out, f, indent=1)
    print("saved", OUT_PATH)

    for c in PARAM_COLS:
        vals = np.array(imp_malate_ens[c])
        print(f"  {c:30s} malate_imp median={np.median(vals):.4f}  IQR=[{np.percentile(vals,25):.4f}, {np.percentile(vals,75):.4f}]")
