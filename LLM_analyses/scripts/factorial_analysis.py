"""
Full 2^4 factorial decomposition of design T5 (dilution D=10, CS KD 1%,
fumarase KD 1%, serine 50 mM) relative to the true base case (D=200, no KD,
serine 1 mM), to determine whether the four levers combine additively or
show synergy/antagonism, and to quantify exactly which lever PAIRS interact.

Runs all 16 on/off combinations of the 4 factors (holding heterologous
enzyme levels, duration, and all other exogenous additions at the base
recipe throughout), then:
  1) Checks simple additivity: does the sum of the 4 individual (single-
     lever) deltas equal the full 4-lever (T5) delta?
  2) Fits the exact (saturated, no residual) factorial effects model on
     coded +-1 factors, giving every main effect and 2-/3-/4-way
     interaction effect on final malate.
  3) Pulls intermediate-pool diagnostics (Pyruvate, Oxaloacetate, NADH,
     L_Serine) at t=end for a few key conditions to explain WHY the
     interactions have the sign/magnitude they do.

Usage:
    ./.venv_sim/Scripts/python.exe scripts/factorial_analysis.py
"""
import sys
import os
import itertools
from pathlib import Path
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from sim2 import Design, run_design

ROOT = Path(__file__).resolve().parent.parent
OUT_DIR = ROOT / "results" / "analysis"
OUT_DIR.mkdir(exist_ok=True)

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

rows = []
diagnostics = {}
species_of_interest = ["Pyruvate", "Oxaloacetate", "L_Serine", "NADH", "SMalate"]

for bits in itertools.product([0, 1], repeat=4):
    design = build_design(bits)
    res = run_design(design)
    cols = res["colnames"]
    idx = {c: i for i, c in enumerate(cols)}
    raw = res["raw_result"]
    row = dict(zip(FACTOR_NAMES, bits))
    row["final_malate_mM"] = res["final_malate_mM"]
    row["delta_malate_mM"] = res["delta_malate_mM"]
    row["rTCA_to_TCA_ratio"] = res["rTCA_to_TCA_ratio"]
    rows.append(row)
    diagnostics[bits] = {sp: (raw[0, idx[sp]], raw[-1, idx[sp]]) for sp in species_of_interest}
    print(f"D={bits[0]} CS={bits[1]} Fum={bits[2]} Ser={bits[3]}  "
          f"final_malate={res['final_malate_mM']:.5f}  ratio={res['rTCA_to_TCA_ratio']:.2f}")

df = pd.DataFrame(rows)
df.to_csv(OUT_DIR / "factorial_design_results.csv", index=False)

base_row = df[(df[FACTOR_NAMES] == 0).all(axis=1)].iloc[0]
t5_row = df[(df[FACTOR_NAMES] == 1).all(axis=1)].iloc[0]
base_malate = base_row["final_malate_mM"]
t5_malate = t5_row["final_malate_mM"]

print("\n" + "=" * 70)
print("SIMPLE ADDITIVITY CHECK")
print("=" * 70)
print(f"Base case (all off):        final malate = {base_malate:.5f} mM")
print(f"T5 (all four on):           final malate = {t5_malate:.5f} mM")
print(f"Actual combined improvement:              {t5_malate - base_malate:+.5f} mM")

single_deltas = {}
for i, f in enumerate(FACTOR_NAMES):
    bits = [0, 0, 0, 0]
    bits[i] = 1
    single_row = df[(df[FACTOR_NAMES] == bits).all(axis=1)].iloc[0]
    delta = single_row["final_malate_mM"] - base_malate
    single_deltas[f] = delta
    print(f"  {f:10s} alone:            final malate = {single_row['final_malate_mM']:.5f} mM  (delta = {delta:+.5f})")

sum_singles = sum(single_deltas.values())
actual_delta = t5_malate - base_malate
gap = actual_delta - sum_singles
print(f"\nSum of 4 individual deltas:                {sum_singles:+.5f} mM")
print(f"Actual T5 delta:                            {actual_delta:+.5f} mM")
print(f"Interaction gap (actual - sum of singles):  {gap:+.5f} mM  "
      f"({'SUPER-additive / synergy' if gap > 0 else 'SUB-additive / antagonism'}, "
      f"{gap/sum_singles*100:+.1f}% relative to the naive sum)")

# ---------------------------------------------------------------------------
# Exact saturated factorial-effects model (coded +-1), via Yates-style
# contrasts: effect_S = (1/8) * sum over all 16 runs of (product of signs
# for factors in S) * y_run   [2^4=16 runs, 2^(4-1)=8 comparisons per effect]
# ---------------------------------------------------------------------------
coded = df[FACTOR_NAMES].apply(lambda c: 2 * c - 1)  # 0/1 -> -1/+1
y = df["final_malate_mM"].values

def effect(subset):
    sign = np.ones(len(df))
    for f in subset:
        sign = sign * coded[f].values
    return float(np.mean(sign * y) * 2)  # standard factorial effect estimate

print("\n" + "=" * 70)
print("FULL FACTORIAL EFFECTS DECOMPOSITION (exact, on final malate mM)")
print("=" * 70)
print("(main effect = average change in malate from turning that factor on,")
print(" averaged over all 8 settings of the other three factors)")

main_effects = {f: effect([f]) for f in FACTOR_NAMES}
for f, e in sorted(main_effects.items(), key=lambda x: -abs(x[1])):
    print(f"  Main effect  {f:20s} {e:+.5f} mM")

print("\n2-way interaction effects (>0 = synergy between that pair, <0 = antagonism):")
two_way = {}
for a, b in itertools.combinations(FACTOR_NAMES, 2):
    e = effect([a, b])
    two_way[(a, b)] = e
for (a, b), e in sorted(two_way.items(), key=lambda x: -abs(x[1])):
    print(f"  {a:10s} x {b:10s}          {e:+.5f} mM")

print("\n3-way interaction effects:")
three_way = {}
for combo in itertools.combinations(FACTOR_NAMES, 3):
    e = effect(list(combo))
    three_way[combo] = e
for combo, e in sorted(three_way.items(), key=lambda x: -abs(x[1])):
    print(f"  {' x '.join(combo):32s} {e:+.5f} mM")

four_way = effect(FACTOR_NAMES)
print(f"\n4-way interaction effect (all four): {four_way:+.5f} mM")

grand_mean = float(np.mean(y))
print(f"\nGrand mean across all 16 runs: {grand_mean:.5f} mM")
print(f"Sum of all main effects:       {sum(main_effects.values()):+.5f} mM")
print(f"Sum of all 2-way effects:      {sum(two_way.values()):+.5f} mM")
print(f"Sum of all 3-way effects:      {sum(three_way.values()):+.5f} mM")
print(f"4-way effect:                  {four_way:+.5f} mM")

# ---------------------------------------------------------------------------
# Mechanistic diagnostics: intermediate pool behavior for key conditions
# ---------------------------------------------------------------------------
print("\n" + "=" * 70)
print("MECHANISTIC DIAGNOSTICS (init -> final, mM)")
print("=" * 70)
key_conditions = {
    "Base (0,0,0,0)": (0, 0, 0, 0),
    "Dilution only (1,0,0,0)": (1, 0, 0, 0),
    "CS+Fum KD only (0,1,1,0)": (0, 1, 1, 0),
    "Serine only (0,0,0,1)": (0, 0, 0, 1),
    "Dilution+Serine (1,0,0,1)": (1, 0, 0, 1),
    "T5 -- all four (1,1,1,1)": (1, 1, 1, 1),
}
for label, bits in key_conditions.items():
    print(f"\n{label}:")
    for sp in species_of_interest:
        init, final = diagnostics[bits][sp]
        print(f"  {sp:14s} {init:.5f} -> {final:.5f}")

print(f"\nWrote {OUT_DIR / 'factorial_design_results.csv'}")
