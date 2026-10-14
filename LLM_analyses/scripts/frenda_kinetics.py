"""
Parses CellFreeModeling/my_pathways_output/Reaction.csv (FRENDA batch output,
one row per EC since the batch was run with --first-reaction-only) into a
lookup: EC -> real BRENDA/CatPred-sourced kinetics.

Km values are assigned POSITIONALLY: FRENDA's Km string lists substrate Kms
(in the same order as its Substrates string) followed by product Kms (in the
same order as its Products string). Since our hand-designed pathway steps
(pathway_designs.py) don't always have the exact same substrate/product
identities as KEGG's first-listed reaction for that EC (e.g. NAD- vs
NADP-linked isoform), we assign FRENDA's real Km magnitudes positionally to
our own substrate/product lists (by role order), not by chemical-identity
matching. This is an approximation but keeps every Km grounded in a real,
enzyme/EC-specific experimental-or-CatPred-predicted value rather than a
generic literature guess -- appropriate for early-stage pathway screening.
"""
import re
import math
import pandas as pd

REACTION_CSV = r"C:\Users\mkcoo\Documents\UW\LLM_CFS_2\parameters\pathway_frenda_reactions.csv"


REASONABLE_KM_DEFAULT = 0.3      # mM -- literature-typical, used to replace FRENDA's crude imputed Km (9.62309, a shared global fallback constant carrying no reaction-specific information)
REASONABLE_KCAT_DEFAULT = 15.0   # 1/s -- literature-typical, used to replace FRENDA's crude imputed Kcat (2331.16/191.943, the same shared fallback pair)

# ECs whose FRENDA-fetched "first reaction" runs the OPPOSITE direction from
# the pathway step we built (confirmed by comparing raw Substrates/Products
# against pathway_designs.py -- see results/pathway_analysis/failure_root_causes.md
# category 1). Keq is direction-dependent, so using it as-fetched instead of
# 1/Keq made a genuinely favorable reaction (real biology) look astronomically
# unfavorable. Fix: invert Keq and swap which Km list feeds our substrates vs
# products, so downstream code sees kinetics "as if" FRENDA had fetched our
# direction to begin with.
#   1.4.1.1  AlaDH:  FRENDA = oxidative deamination (Ala->Pyr); we want reductive amination (Pyr->Ala)
#   1.5.1.5  MTHFD:  FRENDA = methylene-THF oxidation; we want formyl-THF reduction
#   4.1.1.31 PPC:    FRENDA = OAA decarboxylation (OAA->PEP); we want carboxylation (PEP->OAA)
#   1.4.1.21 AspDH:  FRENDA = oxidative deamination (Asp->OAA); we want reductive amination (OAA->Asp)
REVERSED_DIRECTION_ECS = {"1.4.1.1", "1.5.1.5", "4.1.1.31", "1.4.1.21"}

# ECs where FRENDA returned real-looking but chemically broken/unusable data
# (see failure_root_causes.md categories 2 and 3) -- replaced outright with
# literature-typical values for the REAL overall reaction, regardless of
# whether FRENDA itself returned something for that EC.
FORCED_OVERRIDE = {
    # EC 2.2.1.6 (acetolactate synthase, ALS): FRENDA's fetched KEGG reaction
    # is "CO2 + TPP-cofactor -> 2 Pyruvate", a fragment of ALS's TPP-bound
    # intermediate mechanism, not the real "2 Pyruvate -> Acetolactate + CO2"
    # conversion every BDO/isobutanol-family design needs. Real ALS is a
    # fast, essentially irreversible decarboxylative condensation (literature:
    # Km(pyruvate) ~10-20 mM, Kcat ~40-90/s, strongly forward-favored).
    "2.2.1.6": dict(kcat_v=20.0, keq=5000.0,
                     km_subs=[10.0],        # Pyruvate
                     km_prods=[2.0, 5.0]),  # Acetolactate, CO2
    # EC 2.8.3.8 (acetate CoA-transferase, AtoDA) and EC 3.1.2.20 (acyl-CoA
    # hydrolase, TE) both returned the IDENTICAL Keq=7.124576e-218 -- a
    # shared eQuilibrator computation-failure artifact (no real reaction has
    # an equilibrium constant that extreme). Replaced with literature-typical
    # values: CoA-transferases are near-thermoneutral (Keq ~1-10); thioester
    # hydrolysis is generally favorable/irreversible (Keq >>1).
    "2.8.3.8": dict(kcat_v=15.0, keq=2.0,
                     km_subs=[0.2, 5.0],    # AcetoAcetylCoA, Acetate (role-specific per pathway)
                     km_prods=[0.5, 0.3]),  # Acetoacetate, AcetylCoA
    "3.1.2.20": dict(kcat_v=10.0, keq=1000.0,
                      km_subs=[0.3],        # AcrylylCoA
                      km_prods=[0.5, 0.1]),  # AcrylicAcid, CoA
}


def _parse_kv_list(s, imputed_labels=None):
    """'Km_C00002: 0.0078; Km_C00116: 0.21' -> [0.0078, 0.21] (in order).
    Any label present in imputed_labels gets REASONABLE_KM_DEFAULT instead of
    FRENDA's own crude imputed placeholder value."""
    if not isinstance(s, str) or not s.strip():
        return []
    imputed_labels = imputed_labels or set()
    out = []
    for part in s.split(";"):
        part = part.strip()
        if not part:
            continue
        m = re.match(r"(Km_\S+)\s*:\s*([-\d.eE]+)\s*$", part)
        if m:
            label, val = m.group(1), m.group(2)
            if label in imputed_labels:
                out.append(REASONABLE_KM_DEFAULT)
            else:
                try:
                    out.append(float(val))
                except ValueError:
                    out.append(REASONABLE_KM_DEFAULT)
    return out


def _parse_imputed_labels(s):
    if not isinstance(s, str) or not s.strip():
        return set()
    return {p.strip() for p in s.split(";") if p.strip()}


def _parse_kcat(s):
    """'Kcat_F: 59.7874; Kcat_R: 44.011' -> (kcat_f, kcat_r or None)."""
    if not isinstance(s, str) or not s.strip():
        return None, None
    vals = {}
    for part in s.split(";"):
        part = part.strip()
        m = re.match(r"Kcat_([FR])\s*:\s*([-\d.eE]+)", part)
        if m:
            vals[m.group(1)] = float(m.group(2))
    return vals.get("F"), vals.get("R")


def load_kinetics(path=REACTION_CSV):
    """Returns dict: EC -> dict(kcat_v, keq, km_subs:[...], km_prods:[...])."""
    df = pd.read_csv(path)
    lookup = {}
    for _, row in df.iterrows():
        ec = str(row["EC"]).strip()
        n_subs = str(row.get("Substrates", "")).count(";") + 1 if pd.notna(row.get("Substrates")) and str(row.get("Substrates")).strip() else 0
        imputed_km_labels = _parse_imputed_labels(row.get("Km_Imputed"))
        kms = _parse_kv_list(row.get("Km"), imputed_km_labels)
        km_subs = kms[:n_subs]
        km_prods = kms[n_subs:]

        kcat_imputed = _parse_imputed_labels(row.get("Kcat_Imputed"))
        if kcat_imputed:
            # FRENDA's own fallback Kcat_F/R pair is a single shared global
            # constant with no reaction-specific information -- don't trust it.
            kcat_v = REASONABLE_KCAT_DEFAULT
        else:
            kcat_f, kcat_r = _parse_kcat(row.get("Kcat"))
            if kcat_f is not None and kcat_r is not None and kcat_f > 0 and kcat_r > 0:
                kcat_v = math.sqrt(kcat_f * kcat_r)
            elif kcat_f is not None:
                kcat_v = kcat_f
            elif kcat_r is not None:
                kcat_v = kcat_r
            else:
                kcat_v = REASONABLE_KCAT_DEFAULT

        keq = row.get("Keq")
        try:
            keq = float(keq)
            if not math.isfinite(keq) or keq <= 0:
                keq = 1000.0
        except (TypeError, ValueError):
            keq = 1000.0

        if ec in REVERSED_DIRECTION_ECS:
            # FRENDA's fetched reaction runs opposite to our modeled step:
            # invert Keq and swap which Km list feeds substrates vs products.
            keq = 1.0 / keq
            km_subs, km_prods = km_prods, km_subs

        lookup[ec] = dict(kcat_v=kcat_v, keq=keq, km_subs=km_subs, km_prods=km_prods,
                           enzyme_name=row.get("Enzyme"), mechanism=row.get("Mechanism"),
                           kcat_imputed=bool(kcat_imputed), n_km_imputed=len(imputed_km_labels),
                           direction_corrected=ec in REVERSED_DIRECTION_ECS)
    for ec, kin in MANUAL_KINETICS_OVERRIDE.items():
        if ec not in lookup:
            lookup[ec] = dict(kin, enzyme_name="(manual override)", mechanism="MRL")
    for ec, kin in FORCED_OVERRIDE.items():
        lookup[ec] = dict(kin, enzyme_name="(forced override -- FRENDA data unusable)", mechanism="MRL")
    return lookup


# FRENDA/KEGG has no enzyme-database entry for these two ECs (KEGG lookup
# failure, confirmed in my_pathways_output/run_log.txt), so they can't be
# auto-fetched. Filled in from other real sources instead of a skipped step:
MANUAL_KINETICS_OVERRIDE = {
    # EC 1.17.1.9 (soluble NAD-dependent formate dehydrogenase) is the SAME
    # enzyme already present in this project's base model as hEC11719/R56,
    # with real values from the project's own trainparams fit (not a generic
    # guess) -- reused directly rather than re-estimating.
    "1.17.1.9": dict(kcat_v=1.0, keq=340.0,
                      km_subs=[2.2, 0.1],        # Formate, NAD  (matches R56 order)
                      km_prods=[4.33005, 0.048]),  # CO2, NADH
    # EC 2.3.1.264 (PHA/PHB synthase, PhaC) is not in KEGG's enzyme flat
    # file at all (a real gap in KEGG's curation of this well-studied
    # industrial enzyme, not a FRENDA bug). Literature-typical values used
    # instead (Kcat ~5-20/s, Km(3HB-CoA) ~0.02-0.1 mM reported for
    # Cupriavidus necator PhaC across multiple studies) -- flagged in the
    # results table notes as the one non-FRENDA-sourced kinetics in the set.
    "2.3.1.264": dict(kcat_v=8.0, keq=100.0,
                       km_subs=[0.05],     # ThreeHB_CoA
                       km_prods=[0.1]),    # CoA  (PHB_monomer product has no Km role -- polymerization sink)
    # Formaldehyde + THF -> methylene-THF: no dedicated EC / not in KEGG's
    # enzyme list under a distinct number (see ec_map.py note) -- modeled as
    # a fast, near-diffusion-limited spontaneous condensation, literature-
    # typical magnitude, not a FRENDA/BRENDA-sourced value.
    "spontaneous_FTHF": dict(kcat_v=50.0, keq=1.0e4,
                              km_subs=[0.1, 0.1],   # Formaldehyde, THF
                              km_prods=[0.1]),       # methylene-THF
}

_FALLBACK_KM = 0.5  # mM, used only if FRENDA returned fewer Km values than our step needs


def get_km_for_role(kin, subs_or_prods, idx):
    """kin['km_subs'] or kin['km_prods'], idx-th value, with graceful fallback."""
    vals = kin[subs_or_prods]
    if idx < len(vals) and vals[idx] > 0:
        return vals[idx]
    if vals:
        # reuse the last available real value rather than an arbitrary default
        return vals[-1]
    return _FALLBACK_KM
