"""
Data generation for the "masked secondary bottleneck" claim about citrate
synthase (eEC2331 / reaction R46):

  At the baseline dilution (D=200), enzyme concentrations are low enough that
  citrate synthase carries little flux, so its knockdown barely moves malate
  yield. Reducing dilution to D=10 raises both heterologous and endogenous
  enzyme concentrations, activating citrate synthase as a real competing sink
  for oxaloacetate (OAA) and acetyl-CoA (AcCoA) -- the same two substrates the
  engineered rTCA route needs. Knocking down citrate synthase at D=10 then
  redirects carbon into the engineered pathway, producing a strong synergy
  that is invisible at D=200.

R46 (citrate synthase) is written in the model as:
  R46 : Citrate + CoA + eEC2331 -> eEC2331 + AcetylCoA + H2O + Oxaloacetate
  R46 = R46_f - R46_r
so R46_f is citrate LYSIS (citrate -> AcCoA + OAA) and R46_r is citrate
SYNTHESIS (AcCoA + OAA -> citrate), the classical CS reaction and the
competing sink described in the paragraph. Keq_R46 = 1.6e-7 (see
model/iMC057.txt line 750), so R46_r dominates essentially everywhere R46_f
is negligible -- net synthesis flux is well approximated by R46_r alone, but
we track both plus the net R46 for transparency.

R25 (eEC11137/MDH reverse, the rTCA route that makes malate) is:
  R25 : SMalate + NAD + eEC11137 -> eEC11137 + Oxaloacetate + NADH
  R25 = R25_f - R25_r
R25_r is OAA + NADH -> SMalate + NAD, i.e. the malate-PRODUCING direction and
the other consumer of the same OAA pool that R46_r draws down -- included so
the OAA-competition mechanism is visible directly, not just asserted.

Produces:
  results/analysis/citrate_synthase_synergy_timecourses.csv   (per-timepoint, 4 conditions)
  results/analysis/citrate_synthase_synergy_dilution_scan.csv (per-D, WT vs KD, 12 D values)
  results/analysis/citrate_synthase_synergy_summary.json      (headline numbers for the writeup)
"""
import sys
import os
import json
from pathlib import Path
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from sim2 import (
    Design, _fresh_runner, SPECIES_MAP, HETEROLOGOUS_ENZYMES, HET_REL,
    HET_REF, FITTED_PARAMS,
)

ROOT = Path(__file__).resolve().parent.parent
OUT_DIR = ROOT / "results" / "analysis"
OUT_DIR.mkdir(exist_ok=True)

EXTRA_SELECTIONS = ["R46", "R46_f", "R46_r", "R25", "R25_f", "R25_r",
                    "Oxaloacetate", "AcetylCoA", "Citrate", "SMalate"]


def run_with_flux_tracking(design: Design, n_points=400):
    r = _fresh_runner()
    dil = 1.0 / design.dilution_factor_D
    r.model["dilution_factor"] = dil

    for sp, info in SPECIES_MAP.items():
        if info["driver"].startswith("p_eEC") and sp != "eEC11137":
            mult = design.endog_multipliers.get(sp, 1.0)
            r.model[info["driver"]] = mult

    mdh_info = SPECIES_MAP["eEC11137"]
    mdh_native_mult = design.endog_multipliers.get("eEC11137", 1.0)
    p_mdh = mdh_native_mult + (design.eEC11137_overexpr_uM / 1000.0) / (mdh_info["ref"] * dil * r.model["rel0"])
    r.model["p_eEC11137"] = p_mdh

    for hid in HETEROLOGOUS_ENZYMES:
        target_mM = design.het_enzymes_uM.get(hid, 0.0) / 1000.0
        p_val = target_mM / (HET_REF * dil * r.model[HET_REL[hid]])
        r.model[f"p_{hid}"] = p_val

    for sp, info in SPECIES_MAP.items():
        if not info["driver"].startswith("v"):
            continue
        v_default = FITTED_PARAMS.get(info["driver"], None)
        if v_default is None:
            continue
        add_mM = design.exogenous_mM.get(sp, None)
        if add_mM is not None:
            r.model[info["driver"]] = v_default + add_mM / (info["ref"] * dil)
        else:
            r.model[info["driver"]] = v_default

    r.model["v73"] = 10
    r.reset()

    if design.extra_species_overrides_mM:
        for sid, val in design.extra_species_overrides_mM.items():
            r.model[sid] = val

    r.selections = ["time"] + EXTRA_SELECTIONS
    result = r.simulate(0, design.duration_sec, max(200, n_points))
    res = np.array(result)
    cols = result.colnames
    idx = {c: i for i, c in enumerate(cols)}

    out = {c: res[:, idx[c]] for c in cols}
    out["eEC2331_uM"] = None  # filled by caller if needed
    return out, r


CS_NATIVE = {}
CS_KD = {"eEC2331": 0.01}  # matches factorial_analysis.py's CS knockdown convention

BASE_KW = dict(
    het_enzymes_uM={"hEC6411": 50.0, "hEC11719": 25.0, "hEC43117": 50.0},
    eEC11137_overexpr_uM=25.0,
    duration_sec=14400.0,  # 4 h base case
)

# ---------------------------------------------------------------------------
# 1) Full time courses, 4 conditions: D in {200, 10} x CS in {native, KD}
# ---------------------------------------------------------------------------
conditions = [
    ("D=200, CS native", 200.0, CS_NATIVE),
    ("D=200, CS knockdown", 200.0, CS_KD),
    ("D=10, CS native", 10.0, CS_NATIVE),
    ("D=10, CS knockdown", 10.0, CS_KD),
]

tc_rows = []
for label, D, endog in conditions:
    design = Design(dilution_factor_D=D, endog_multipliers=dict(endog), **BASE_KW)
    out, r = run_with_flux_tracking(design)
    n = len(out["time"])
    for i in range(n):
        tc_rows.append({
            "condition": label,
            "dilution_D": D,
            "cs_knockdown": bool(endog),
            "time_s": out["time"][i],
            "time_h": out["time"][i] / 3600.0,
            "R46_net": out["R46"][i],
            "R46_f_citrate_lysis": out["R46_f"][i],
            "R46_r_citrate_synthesis": out["R46_r"][i],
            "CS_sink_flux": out["R46_r"][i] - out["R46_f"][i],  # net carbon draw into citrate
            "R25_net": out["R25"][i],
            "R25_r_rTCA_malate_production": out["R25_r"][i],
            "Oxaloacetate_mM": out["Oxaloacetate"][i],
            "AcetylCoA_mM": out["AcetylCoA"][i],
            "Citrate_mM": out["Citrate"][i],
            "SMalate_mM": out["SMalate"][i],
        })
    print(f"{label}: final CS sink flux={tc_rows[-1]['CS_sink_flux']:.6e} mM/s, "
          f"final malate={tc_rows[-1]['SMalate_mM']:.5f} mM")

tc_df = pd.DataFrame(tc_rows)
tc_df.to_csv(OUT_DIR / "citrate_synthase_synergy_timecourses.csv", index=False)
print(f"Wrote {OUT_DIR / 'citrate_synthase_synergy_timecourses.csv'} ({len(tc_df)} rows)")

# ---------------------------------------------------------------------------
# 2) Dilution scan: for each D, run CS-native and CS-KD, record peak/mean CS
#    sink flux and final malate, to show WHERE the flux (and the synergy)
#    turns on as D is reduced from 200 toward 10.
# ---------------------------------------------------------------------------
D_SCAN = [10, 15, 20, 30, 50, 75, 100, 125, 150, 200, 300, 400]

scan_rows = []
for D in D_SCAN:
    row = {"dilution_D": D}
    for tag, endog in [("wt", CS_NATIVE), ("kd", CS_KD)]:
        design = Design(dilution_factor_D=D, endog_multipliers=dict(endog), **BASE_KW)
        out, r = run_with_flux_tracking(design)
        sink = out["R46_r"] - out["R46_f"]
        t = out["time"]
        mean_sink = float(np.trapz(sink, t) / (t[-1] - t[0]))
        row[f"peak_CS_sink_flux_mM_per_s_{tag}"] = float(np.max(sink))
        row[f"mean_CS_sink_flux_mM_per_s_{tag}"] = mean_sink
        row[f"final_malate_mM_{tag}"] = float(out["SMalate"][-1])
        row[f"final_OAA_mM_{tag}"] = float(out["Oxaloacetate"][-1])
        row[f"final_AcCoA_mM_{tag}"] = float(out["AcetylCoA"][-1])
    row["cs_kd_delta_malate_mM"] = row["final_malate_mM_kd"] - row["final_malate_mM_wt"]
    row["cs_kd_pct_improvement"] = 100.0 * row["cs_kd_delta_malate_mM"] / row["final_malate_mM_wt"]
    scan_rows.append(row)
    print(f"D={D:>4.0f}  mean_sink_wt={row['mean_CS_sink_flux_mM_per_s_wt']:.4e}  "
          f"mean_sink_kd={row['mean_CS_sink_flux_mM_per_s_kd']:.4e}  "
          f"malate_wt={row['final_malate_mM_wt']:.4f}  malate_kd={row['final_malate_mM_kd']:.4f}  "
          f"KD improvement={row['cs_kd_pct_improvement']:+.1f}%")

scan_df = pd.DataFrame(scan_rows)
scan_df.to_csv(OUT_DIR / "citrate_synthase_synergy_dilution_scan.csv", index=False)
print(f"Wrote {OUT_DIR / 'citrate_synthase_synergy_dilution_scan.csv'} ({len(scan_df)} rows)")

# ---------------------------------------------------------------------------
# 3) Headline summary numbers
# ---------------------------------------------------------------------------
d200 = scan_df[scan_df.dilution_D == 200].iloc[0]
d10 = scan_df[scan_df.dilution_D == 10].iloc[0]

summary = {
    "D200": {
        "mean_CS_sink_flux_wt_mM_per_s": d200["mean_CS_sink_flux_mM_per_s_wt"],
        "mean_CS_sink_flux_kd_mM_per_s": d200["mean_CS_sink_flux_mM_per_s_kd"],
        "final_malate_wt_mM": d200["final_malate_mM_wt"],
        "final_malate_kd_mM": d200["final_malate_mM_kd"],
        "cs_kd_pct_improvement": d200["cs_kd_pct_improvement"],
    },
    "D10": {
        "mean_CS_sink_flux_wt_mM_per_s": d10["mean_CS_sink_flux_mM_per_s_wt"],
        "mean_CS_sink_flux_kd_mM_per_s": d10["mean_CS_sink_flux_mM_per_s_kd"],
        "final_malate_wt_mM": d10["final_malate_mM_wt"],
        "final_malate_kd_mM": d10["final_malate_mM_kd"],
        "cs_kd_pct_improvement": d10["cs_kd_pct_improvement"],
    },
    "fold_increase_in_mean_CS_sink_flux_D200_to_D10_native": float(
        d10["mean_CS_sink_flux_mM_per_s_wt"] / d200["mean_CS_sink_flux_mM_per_s_wt"]
    ),
    "fold_increase_in_CS_KD_benefit_D200_to_D10": float(
        d10["cs_kd_pct_improvement"] / d200["cs_kd_pct_improvement"]
    ) if d200["cs_kd_pct_improvement"] != 0 else None,
}
with open(OUT_DIR / "citrate_synthase_synergy_summary.json", "w") as f:
    json.dump(summary, f, indent=2)
print(f"\nWrote {OUT_DIR / 'citrate_synthase_synergy_summary.json'}")
print(json.dumps(summary, indent=2))
