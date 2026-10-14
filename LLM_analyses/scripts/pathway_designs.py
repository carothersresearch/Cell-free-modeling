"""
30 heterologous pathway designs (10 starting substrates x 3 pathways each) for
the crude-lysate model (model/iMC057.txt).

IMPORTANT METHODOLOGICAL NOTE: none of the new enzymes below have real fitted
kinetic parameters (unlike every other enzyme in this project, which used
parameters/final_params.pkl). Kcat/Km values are literature-typical
ORDER-OF-MAGNITUDE estimates assigned by enzyme class (see KINETIC_TEMPLATES),
not experimentally characterized values. All new reactions are modeled as
irreversible saturating (product-of-Michaelis-Menten-terms) kinetics -- a
standard simplification for early-stage pathway screening, chosen because
most of these steps are genuinely near-irreversible in vivo (decarboxylations,
strongly exergonic reductions) or are pulled forward by downstream removal.
Treat all yield numbers as illustrative/comparative, not validated predictions.

All enzymes are cytoplasmic/soluble (no membrane-associated steps: no
PQQ-dependent periplasmic methanol dehydrogenase, no membrane-bound
oxidases/dehydrogenases, no ETC components, no transporters).

Enzyme species get a flat, dilution-independent initial concentration (mM),
consistent with how this project has always modeled newly-added heterologous
enzymes (added post-dilution to hit a target titer). New non-native product/
intermediate species start at 0 mM. All designs run at D=10 (least dilute
allowed) with a 4 h (14400 s) duration, matching how the original malate
pathway was engineered.
"""

# Generic kinetic parameter templates by enzyme class (Kcat in 1/s, Km in mM).
# All new reactions use: v = Kcat * E * PRODUCT_i(S_i / (Km_i + S_i))
KINETIC_TEMPLATES = {
    "kinase_ATP":        dict(kcat=20.0, km_default=0.10, km_ATP=0.15),
    "dehydrogenase_NAD":  dict(kcat=15.0, km_default=0.30, km_NAD=0.15, km_NADH=0.15),
    "dehydrogenase_NADP": dict(kcat=12.0, km_default=0.30, km_NADP=0.05, km_NADPH=0.05),
    "decarboxylase":     dict(kcat=10.0, km_default=0.50),
    "condensing_lyase":  dict(kcat=8.0,  km_default=0.50),
    "isomerase":         dict(kcat=25.0, km_default=0.20),
    "CoA_transferase":   dict(kcat=10.0, km_default=0.30, km_acceptor=1.00),
    "carboxylase":       dict(kcat=8.0,  km_default=0.30, km_HCO3=2.00),
    "reductive_aminase":  dict(kcat=12.0, km_default=0.30, km_NH3=5.00, km_NADH=0.15),
    "oxidase_direct":    dict(kcat=15.0, km_default=0.40, km_NAD=0.15, km_NADH=0.15),
    "condensing_C1":     dict(kcat=8.0,  km_default=0.40),
}

# Each pathway dict:
#   id, substrate, product, product_is_native (bool), enzymes: list of
#   (enzyme_name, template, reaction_id, substrates:[(name,km_key)], products:[names])
#   cofactor_additions: dict species->mM to add exogenously (beyond starting substrate)
#   substrate_dose_mM: how much of the starting substrate to add
#   module: classification string
#   notes: caveats
#   est_product_msp: (low, high) $/kg

PATHWAYS = []

def add(pid, substrate, substrate_dose, product, product_native, enzymes, cofactors, module, notes, msp):
    PATHWAYS.append(dict(id=pid, substrate=substrate, substrate_dose_mM=substrate_dose,
                          product=product, product_native=product_native, enzymes=enzymes,
                          cofactor_additions=cofactors, module=module, notes=notes, msp=msp))

# ===================== 1. FORMATE =====================
add("Formate_1_Serine", "Formate", 10.0, "L_Serine", True,
    enzymes=[
        ("FTL", "condensing_C1", "N1", [("Formate", "default"), ("Tetrahydrofolate", "default"), ("ATP", "km_ATP")], ["FormylTHF", "ADP", "Orthophosphate"]),
        ("MTHFD", "dehydrogenase_NADP", "N2", [("FormylTHF", "default"), ("NADPH", "km_NADPH")], ["_5_10_Methylenetetrahydrofolate", "NADP"]),
    ],
    cofactors={"ATP": 1.0, "NADPH": 1.0, "Tetrahydrofolate": 1.0, "Glycine": 1.0},
    module="One-carbon/serine metabolism",
    notes="Feeds native SHMT (R49) via the newly-generated methylene-THF pool. Simulation revealed native free THF AND native Glycine are both essentially zero (~0.0001 mM each, trace placeholders, not real standing pools) -- SHMT cannot run forward to make serine without both, so this design requires small exogenous top-ups of both, not just THF. This is a real, simulation-derived correction to the original 'reuses native Glycine' assumption, not a design choice.",
    msp=(5, 10))

add("Formate_2_Pyruvate", "Formate", 10.0, "Pyruvate", True,
    enzymes=[
        ("PFL_rev", "condensing_lyase", "N1", [("Formate", "default"), ("AcetylCoA", "default")], ["Pyruvate", "CoA"]),
    ],
    cofactors={"AcetylCoA": 1.0},
    module="Pyruvate metabolism",
    notes="Reverse pyruvate formate-lyase. Real PFL is an O2-sensitive glycyl-radical enzyme in vivo -- functional expression/activity in an aerobic lysate prep is a genuine open risk, not just a modeling detail.",
    msp=(3, 8))

add("Formate_3_Alanine", "Formate", 10.0, "L_Alanine", False,
    enzymes=[
        ("FDH2", "dehydrogenase_NAD", "N1", [("Formate", "default"), ("NAD", "km_NAD")], ["CO2", "NADH"]),
        ("AlaDH", "reductive_aminase", "N2", [("Pyruvate", "default"), ("Ammonia", "km_NH3"), ("NADH", "km_NADH")], ["L_Alanine", "NAD"]),
    ],
    cofactors={},
    module="Pyruvate metabolism / amino acid biosynthesis",
    notes="Uses the lysate's huge native ammonia pool as the amino donor; FDH regenerates NADH in-line so no exogenous NADH is needed, sidestepping the native NADH scarcity issue.",
    msp=(2, 4))

# ===================== 2. ACETATE =====================
add("Acetate_1_PHB", "Acetate", 10.0, "PHB_monomer", False,
    enzymes=[
        ("PhaB", "dehydrogenase_NADP", "N1", [("AcetoAcetylCoA", "default"), ("NADPH", "km_NADPH")], ["ThreeHB_CoA", "NADP"]),
        ("PhaC", "condensing_lyase", "N2", [("ThreeHB_CoA", "default")], ["PHB_monomer", "CoA"]),
    ],
    cofactors={"NADPH": 1.0},
    module="Fatty acid synthesis (novel polymer branch)",
    notes="Reuses native acetate activation (R48) and native thiolase (R5) already in the model. PHB is an insoluble polymer -- PHB_monomer tracks flux into the product, not actual polymer accumulation/precipitation, which this ODE framework cannot represent.",
    msp=(2, 4))

add("Acetate_2_Acetone", "Acetate", 10.0, "Acetone", False,
    enzymes=[
        ("AtoDA", "CoA_transferase", "N1", [("AcetoAcetylCoA", "default"), ("Acetate", "km_acceptor")], ["Acetoacetate", "AcetylCoA"]),
        ("AdvB", "decarboxylase", "N2", [("Acetoacetate", "default")], ["Acetone", "CO2"]),
    ],
    cofactors={},
    module="Fatty acid/ketone-body-adjacent metabolism",
    notes="Classic Clostridium acetone pathway, reusing native thiolase (R5). A second acetate molecule is consumed as the CoA acceptor -- net 2 acetate in per acetone out.",
    msp=(0.8, 1.5))

add("Acetate_3_Ethanol", "Acetate", 10.0, "Ethanol", False,
    enzymes=[
        ("ALDH", "dehydrogenase_NAD", "N1", [("AcetylCoA", "default"), ("NADH", "km_NADH")], ["Acetaldehyde", "CoA", "NAD"]),
        ("ADH", "dehydrogenase_NAD", "N2", [("Acetaldehyde", "default"), ("NADH", "km_NADH")], ["Ethanol", "NAD"]),
    ],
    cofactors={"NADH": 1.0},
    module="Pyruvate/fermentation metabolism",
    notes="Reuses native acetate activation (R48) to reach AcetylCoA. Cheapest-value product of the ten Acetate/Formate/Pyruvate/Lactate options -- included for completeness and as a simple, well-validated benchmark route.",
    msp=(0.5, 0.8))

# ===================== 3. PYRUVATE =====================
add("Pyruvate_1_BDO", "Pyruvate", 10.0, "BDO23", False,
    enzymes=[
        ("ALS", "condensing_lyase", "N1", [("Pyruvate", "default", 2)], ["Acetolactate", "CO2"]),
        ("ALDC", "decarboxylase", "N2", [("Acetolactate", "default")], ["Acetoin", "CO2"]),
        ("BDH", "dehydrogenase_NAD", "N3", [("Acetoin", "default"), ("NADH", "km_NADH")], ["BDO23", "NAD"]),
    ],
    cofactors={"NADH": 1.0},
    module="Pyruvate metabolism",
    notes="Classic, well-validated 3-enzyme industrial route (2,3-butanediol). ALS formally needs 2 pyruvate per acetolactate.",
    msp=(1, 2))

add("Pyruvate_2_Isobutanol", "Pyruvate", 10.0, "Isobutanol", False,
    enzymes=[
        ("ALS2", "condensing_lyase", "N1", [("Pyruvate", "default", 2)], ["Acetolactate", "CO2"]),
        ("KARI", "dehydrogenase_NADP", "N2", [("Acetolactate", "default"), ("NADPH", "km_NADPH")], ["DihydroxyIsovalerate", "NADP"]),
        ("DHAD", "isomerase", "N3", [("DihydroxyIsovalerate", "default")], ["KetoIsovalerate"]),
        ("KDC", "decarboxylase", "N4", [("KetoIsovalerate", "default")], ["Isobutyraldehyde", "CO2"]),
        ("ADH2", "dehydrogenase_NAD", "N5", [("Isobutyraldehyde", "default"), ("NADH", "km_NADH")], ["Isobutanol", "NAD"]),
    ],
    cofactors={"NADPH": 1.0, "NADH": 1.0},
    module="Pyruvate / branched-chain amino acid-adjacent metabolism",
    notes="The Atsumi/Liao valine-pathway-derived isobutanol route -- 5 enzymes, the most complex design in this set. DHAD uses an Fe-S cluster (O2-sensitive in vivo, though not membrane-associated) -- a real activity risk to flag, not a blocker.",
    msp=(1, 1.5))

add("Pyruvate_3_Alanine", "Pyruvate", 10.0, "L_Alanine", False,
    enzymes=[
        ("FDH3", "dehydrogenase_NAD", "N1", [("Formate", "default"), ("NAD", "km_NAD")], ["CO2", "NADH"]),
        ("AlaDH2", "reductive_aminase", "N2", [("Pyruvate", "default"), ("Ammonia", "km_NH3"), ("NADH", "km_NADH")], ["L_Alanine", "NAD"]),
    ],
    cofactors={"Formate": 10.0},
    module="Amino acid biosynthesis / pyruvate metabolism",
    notes="Simplest possible pyruvate valorization (1 net conversion), but still needs an NADH source -- a small formate co-feed is included specifically to regenerate NADH in-line rather than assuming native NADH (which is scarce, ~0.003 mM) is sufficient alone.",
    msp=(2, 4))

# ===================== 4. BICARBONATE/CO2 =====================
add("HCO3_1_Malate", "HCO3", 10.0, "SMalate", True,
    enzymes=[
        ("PYC_new", "carboxylase", "N1", [("Pyruvate", "default"), ("ATP", "default"), ("HCO3", "km_HCO3")], ["Oxaloacetate", "ADP"]),
    ],
    cofactors={"ATP": 1.0, "Pyruvate": 0.0},
    module="TCA cycle / anaplerosis (reductive route)",
    notes="Single new enzyme (pyruvate carboxylase); reuses native MDH (eEC11137, already in the model) to finish OAA->malate. This recreates a simplified, standalone version of the project's original engineered route, now correctly scoped to only what's native plus this one addition -- no assumption of the other heterologous enzymes.",
    msp=(2, 3))

add("HCO3_2_3HP", "HCO3", 10.0, "ThreeHP", False,
    enzymes=[
        ("MCR", "dehydrogenase_NADP", "N1", [("Malonyl_CoA", "default"), ("NADPH", "km_NADPH")], ["ThreeHP", "CoA", "NADP"]),
    ],
    cofactors={"NADPH": 2.0, "AcetylCoA": 1.0, "ATP": 1.0},
    module="Fatty acid synthesis-adjacent (malonyl-CoA branch)",
    notes="Reuses native acetyl-CoA carboxylase (R6, already in the model) for the actual HCO3-fixation step -- only the downstream bifunctional malonyl-CoA reductase is new. 3-HP is a top-tier DOE platform chemical (precursor to acrylic acid).",
    msp=(2, 4))

add("HCO3_3_Aspartate", "HCO3", 10.0, "L_Aspartate", False,
    enzymes=[
        ("PPC", "carboxylase", "N1", [("Phosphoenolpyruvate", "default"), ("HCO3", "km_HCO3")], ["Oxaloacetate", "Orthophosphate"]),
        ("AspDH", "reductive_aminase", "N2", [("Oxaloacetate", "default"), ("Ammonia", "km_NH3"), ("NADH", "km_NADH")], ["L_Aspartate", "NAD"]),
    ],
    cofactors={"NADH": 1.0},
    module="Amino acid biosynthesis / anaplerosis",
    notes="PEP carboxylase is a cleaner, higher-affinity carboxylase than PYC and is not the same enzyme used in the malate route. Draws on native PEP and the lysate's abundant native ammonia pool for the amination step.",
    msp=(2, 4))

# ===================== 5. S-LACTATE =====================
add("Lactate_1_BDO", "S_Lactate", 10.0, "BDO23", False,
    enzymes=[
        ("ALS3", "condensing_lyase", "N1", [("Pyruvate", "default", 2)], ["Acetolactate", "CO2"]),
        ("ALDC2", "decarboxylase", "N2", [("Acetolactate", "default")], ["Acetoin", "CO2"]),
        ("BDH2", "dehydrogenase_NAD", "N3", [("Acetoin", "default"), ("NADH", "km_NADH")], ["BDO23", "NAD"]),
    ],
    cofactors={"NADH": 1.0},
    module="Pyruvate metabolism (via native LDH, R53)",
    notes="Lactate feeds pyruvate via the ALREADY-NATIVE LDH (R53) -- no new enzyme needed for the entry step itself, only the same 3-enzyme BDO exit used for pyruvate. Demonstrates how tightly lactate and pyruvate are linked in this network.",
    msp=(1, 2))

add("Lactate_2_Alanine", "S_Lactate", 10.0, "L_Alanine", False,
    enzymes=[
        ("AlaDH3", "reductive_aminase", "N1", [("Pyruvate", "default"), ("Ammonia", "km_NH3"), ("NADH", "km_NADH")], ["L_Alanine", "NAD"]),
    ],
    cofactors={"NADH": 1.0},
    module="Amino acid biosynthesis / pyruvate metabolism",
    notes="Simplest design in the whole set (1 new enzyme) -- native LDH (R53) already converts lactate to pyruvate; only the reductive amination is new.",
    msp=(2, 4))

add("Lactate_3_AcrylicAcid", "S_Lactate", 10.0, "AcrylicAcid", False,
    enzymes=[
        ("PCT", "CoA_transferase", "N1", [("S_Lactate", "default"), ("AcetylCoA", "km_acceptor")], ["LactylCoA", "Acetate"]),
        ("LCD", "isomerase", "N2", [("LactylCoA", "default")], ["AcrylylCoA"]),
        ("TE", "condensing_lyase", "N3", [("AcrylylCoA", "default")], ["AcrylicAcid", "CoA"]),
    ],
    cofactors={"AcetylCoA": 1.0},
    module="Fatty acid/CoA-thioester metabolism (novel branch)",
    notes="The Clostridium propionicum-derived route -- the one lactate-specific design that does NOT just funnel through pyruvate. Lactoyl-CoA dehydratase is a genuinely difficult, oxygen-sensitive radical enzyme in real biology -- flag this as the highest-risk expression/activity item in the whole 30-design set, not a minor caveat.",
    msp=(1.5, 2))

# ===================== 6. GLYCEROL =====================
add("Glycerol_1_PDO12", "Glycerol", 10.0, "PDO12", False,
    enzymes=[
        ("GlpK", "kinase_ATP", "N1", [("Glycerol", "default"), ("ATP", "km_ATP")], ["Glycerol3P", "ADP"]),
        ("GPD1", "dehydrogenase_NAD", "N2", [("Glycerol3P", "default"), ("NAD", "km_NAD")], ["Glycerone_phosphate", "NADH"]),
        ("MgsA", "isomerase", "N3", [("Glycerone_phosphate", "default")], ["Methylglyoxal", "Orthophosphate"]),
        ("FucO", "dehydrogenase_NAD", "N4", [("Methylglyoxal", "default"), ("NADH", "km_NADH")], ["PDO12", "NAD"]),
    ],
    cofactors={"ATP": 1.0, "NADH": 1.0},
    module="Glycolysis-adjacent (most dilution-protected native module)",
    notes="Enters at DHAP (Glycerone_phosphate), the node found to be the most dilution-robust in the whole network (alpha=1.42). Methylglyoxal is a reactive electrophile -- real risk of off-target lysate protein damage if it accumulates faster than FucO consumes it.",
    msp=(1, 2))

add("Glycerol_2_BDO", "Glycerol", 10.0, "BDO23", False,
    enzymes=[
        ("GlpK2", "kinase_ATP", "N1", [("Glycerol", "default"), ("ATP", "km_ATP")], ["Glycerol3P", "ADP"]),
        ("GPD1b", "dehydrogenase_NAD", "N2", [("Glycerol3P", "default"), ("NAD", "km_NAD")], ["Glycerone_phosphate", "NADH"]),
        ("ALS4", "condensing_lyase", "N3", [("Pyruvate", "default", 2)], ["Acetolactate", "CO2"]),
        ("ALDC3", "decarboxylase", "N4", [("Acetolactate", "default")], ["Acetoin", "CO2"]),
        ("BDH3", "dehydrogenase_NAD", "N5", [("Acetoin", "default"), ("NADH", "km_NADH")], ["BDO23", "NAD"]),
    ],
    cofactors={"ATP": 1.0, "NADH": 1.0},
    module="Glycolysis + pyruvate metabolism",
    notes="Longest route in the set (5 enzymes): glycerol enters at DHAP, flows through the rest of native glycolysis to pyruvate, then exits via the standard BDO route. Tests whether a fully-integrated entry+exit design pays off relative to the shorter, single-node options.",
    msp=(1, 2))

add("Glycerol_3_DHA", "Glycerol", 10.0, "DHA", False,
    enzymes=[
        ("GlyDH", "dehydrogenase_NAD", "N1", [("Glycerol", "default"), ("NAD", "km_NAD")], ["DHA", "NADH"]),
    ],
    cofactors={},
    module="Glycolysis-adjacent (novel terminal product)",
    notes="Single-enzyme design that deliberately AVOIDS ATP entirely (NAD+-dependent oxidation, not a kinase) -- sidesteps the native ATP-scarcity constraint identified for the other two glycerol routes. DHA has no native consuming reaction, so it simply accumulates as a clean terminal product.",
    msp=(5, 8))

# ===================== 7. D-XYLOSE =====================
add("Xylose_1_Xylitol", "D_Xylose", 10.0, "Xylitol", False,
    enzymes=[
        ("XR", "dehydrogenase_NADP", "N1", [("D_Xylose", "default"), ("NADPH", "km_NADPH")], ["Xylitol", "NADP"]),
    ],
    cofactors={"NADPH": 1.0},
    module="Pentose phosphate pathway-adjacent (novel terminal product)",
    notes="Single-enzyme design (xylose reductase). Xylitol has no native consuming reaction and just accumulates. NADPH is scarce natively (~0.003 mM) -- the exogenous NADPH addition is not optional here, it's load-bearing.",
    msp=(3, 5))

add("Xylose_2_XylonicAcid", "D_Xylose", 10.0, "XylonicAcid", False,
    enzymes=[
        ("XDH", "dehydrogenase_NAD", "N1", [("D_Xylose", "default"), ("NAD", "km_NAD")], ["XylonicAcid", "NADH"]),
    ],
    cofactors={},
    module="Pentose phosphate pathway-adjacent (novel terminal product)",
    notes="Single-enzyme oxidative route (xylose dehydrogenase) -- the xylose analogue of glucose's gluconic acid design. Net NADH producer rather than consumer, unlike most of this set.",
    msp=(2, 4))

add("Xylose_3_BDO", "D_Xylose", 10.0, "BDO23", False,
    enzymes=[
        ("XylA", "isomerase", "N1", [("D_Xylose", "default")], ["D_Xylulose_5_phosphate"]),
        ("XylB", "kinase_ATP", "N2", [("D_Xylulose_5_phosphate", "default"), ("ATP", "km_ATP")], ["D_Ribulose_5_phosphate", "ADP"]),
        ("ALS5", "condensing_lyase", "N3", [("Pyruvate", "default", 2)], ["Acetolactate", "CO2"]),
        ("ALDC4", "decarboxylase", "N4", [("Acetolactate", "default")], ["Acetoin", "CO2"]),
        ("BDH4", "dehydrogenase_NAD", "N5", [("Acetoin", "default"), ("NADH", "km_NADH")], ["BDO23", "NAD"]),
    ],
    cofactors={"ATP": 1.0, "NADH": 1.0},
    module="Pentose phosphate pathway + pyruvate metabolism",
    notes="XylA/XylB feed the native PPP (moderately dilution-sensitive, alpha=2.37), which connects through to native glycolysis and pyruvate; exits via the standard BDO route. Note: XylA/XylB as written here are simplified stand-ins for the real isomerase+kinase pair and don't exactly reproduce native PPP stoichiometry.",
    msp=(1, 2))

# ===================== 8. METHANOL =====================
add("Methanol_1_G6P", "Methanol", 10.0, "DGlucose6phosphate", True,
    enzymes=[
        ("MDH_sol", "dehydrogenase_NAD", "N1", [("Methanol", "default"), ("NAD", "km_NAD")], ["Formaldehyde", "NADH"]),
        ("Hps", "condensing_C1", "N2", [("Formaldehyde", "default"), ("D_Ribulose_5_phosphate", "default")], ["Hexulose6P"]),
        ("Phi", "isomerase", "N3", [("Hexulose6P", "default")], ["D_Fructose_6_phosphate"]),
    ],
    cofactors={},
    module="Glycolysis/PPP (synthetic methylotrophy, RuMP pathway)",
    notes="Uses the SOLUBLE NAD-dependent methanol dehydrogenase variant (e.g. Bacillus methanolicus-type), deliberately NOT the classic PQQ-dependent periplasmic/membrane-associated MDH, per the no-membrane-enzyme constraint. Classic RuMP formaldehyde-fixation pathway; feeds native Fructose-6-P, tracked here via the connected native G6P pool.",
    msp=(10, 30))

add("Methanol_2_Alanine", "Methanol", 10.0, "L_Alanine", False,
    enzymes=[
        ("MDH_sol2", "dehydrogenase_NAD", "N1", [("Methanol", "default"), ("NAD", "km_NAD")], ["Formaldehyde", "NADH"]),
        ("AlaDH4", "reductive_aminase", "N2", [("Pyruvate", "default"), ("Ammonia", "km_NH3"), ("NADH", "km_NADH")], ["L_Alanine", "NAD"]),
    ],
    cofactors={},
    module="Pyruvate metabolism / C1 metabolism",
    notes="Methanol oxidation directly supplies the NADH the reductive amination needs -- no exogenous NADH required, similar in spirit to the Formate-to-Alanine design. Formaldehyde is a reactive, toxic intermediate -- same caveat as the glycerol/methylglyoxal route.",
    msp=(2, 4))

add("Methanol_3_Serine", "Methanol", 10.0, "L_Serine", True,
    enzymes=[
        ("MDH_sol3", "dehydrogenase_NAD", "N1", [("Methanol", "default"), ("NAD", "km_NAD")], ["Formaldehyde", "NADH"]),
        ("FTHF_cond", "condensing_C1", "N2", [("Formaldehyde", "default"), ("Tetrahydrofolate", "default")], ["_5_10_Methylenetetrahydrofolate"]),
    ],
    cofactors={"Tetrahydrofolate": 1.0, "Glycine": 1.0},
    module="One-carbon/serine metabolism",
    notes="Parallels the Formate-to-Serine design one step upstream. Feeds native SHMT (R49) via the newly-generated methylene-THF pool. Same native THF-and-Glycine scarcity finding as Formate_1_Serine applies here (both are ~0.0001 mM trace placeholders, not real pools) -- both are topped up exogenously for the same reason. The two C1-feedstock routes to the same product (Formate and Methanol) are a natural side-by-side comparison in the results table.",
    msp=(5, 10))

# ===================== 9. GLUCOSE =====================
add("Glucose_1_BDO", "Glucose", 10.0, "BDO23", False,
    enzymes=[
        ("HK", "kinase_ATP", "N1", [("Glucose", "default"), ("ATP", "km_ATP")], ["DGlucose6phosphate", "ADP"]),
        ("ALS6", "condensing_lyase", "N2", [("Pyruvate", "default", 2)], ["Acetolactate", "CO2"]),
        ("ALDC5", "decarboxylase", "N3", [("Acetolactate", "default")], ["Acetoin", "CO2"]),
        ("BDH5", "dehydrogenase_NAD", "N4", [("Acetoin", "default"), ("NADH", "km_NADH")], ["BDO23", "NAD"]),
    ],
    cofactors={"ATP": 1.0, "NADH": 1.0},
    module="Glycolysis + pyruvate metabolism",
    notes="The classic, industrially-proven glucose fermentation route. Hexokinase feeds directly into the most dilution-protected native module (glycolysis, alpha=1.42) found in this system.",
    msp=(1, 2))

add("Glucose_2_Mannitol", "Glucose", 10.0, "Mannitol", False,
    enzymes=[
        ("HK2", "kinase_ATP", "N1", [("Glucose", "default"), ("ATP", "km_ATP")], ["DGlucose6phosphate", "ADP"]),
        ("MtDH", "dehydrogenase_NAD", "N2", [("D_Fructose_6_phosphate", "default"), ("NADH", "km_NADH")], ["Mannitol", "NAD"]),
    ],
    cofactors={"ATP": 1.0, "NADH": 1.0},
    module="Glycolysis-adjacent (novel terminal branch)",
    notes="Only 2 enzymes; taps the native glycolytic Fructose-6-P pool (reached via native PGI, R21, from the hexokinase product) for a direct reduction to the sugar-alcohol product, sidestepping the multi-step BDO/pyruvate exit route entirely.",
    msp=(3, 5))

add("Glucose_3_GluconicAcid", "Glucose", 10.0, "GluconicAcid", False,
    enzymes=[
        ("GDH", "dehydrogenase_NAD", "N1", [("Glucose", "default"), ("NAD", "km_NAD")], ["GluconicAcid", "NADH"]),
    ],
    cofactors={},
    module="Novel direct-oxidation branch (bypasses native glycolysis)",
    notes="Single enzyme, and the ONLY glucose design that needs no ATP at all -- direct oxidation instead of phosphorylation, deliberately avoiding the native ATP-scarcity bottleneck identified for the kinase-dependent routes.",
    msp=(1, 2))

# ===================== 10. ETHANOL =====================
add("Ethanol_1_PHB", "Ethanol", 10.0, "PHB_monomer", False,
    enzymes=[
        ("ADH3", "dehydrogenase_NAD", "N1", [("Ethanol", "default"), ("NAD", "km_NAD")], ["Acetaldehyde", "NADH"]),
        ("ALDH2", "dehydrogenase_NAD", "N2", [("Acetaldehyde", "default"), ("NAD", "km_NAD")], ["Acetate", "NADH"]),
        ("PhaB2", "dehydrogenase_NADP", "N3", [("AcetoAcetylCoA", "default"), ("NADPH", "km_NADPH")], ["ThreeHB_CoA", "NADP"]),
        ("PhaC2", "condensing_lyase", "N4", [("ThreeHB_CoA", "default")], ["PHB_monomer", "CoA"]),
    ],
    cofactors={"NADPH": 1.0},
    module="Fatty acid synthesis-adjacent (novel polymer branch)",
    notes="Longest ethanol route (4 enzymes): full oxidation to acetate, then reuses native acetate activation (R48) and native thiolase (R5) exactly as the Acetate_1_PHB design does. Net NADH producer (2x) upstream, NADPH consumer downstream.",
    msp=(2, 4))

add("Ethanol_2_Acetaldehyde", "Ethanol", 10.0, "Acetaldehyde", False,
    enzymes=[
        ("ADH4", "dehydrogenase_NAD", "N1", [("Ethanol", "default"), ("NAD", "km_NAD")], ["Acetaldehyde", "NADH"]),
    ],
    cofactors={},
    module="Simple oxidation / pyruvate-adjacent",
    notes="The simplest design in the entire 30-pathway set (1 enzyme, direct oxidation). Acetaldehyde is volatile and reactive -- real product-recovery and off-target-reactivity concerns despite the trivial pathway design.",
    msp=(1, 1.5))

add("Ethanol_3_Acetone", "Ethanol", 10.0, "Acetone", False,
    enzymes=[
        ("ADH5", "dehydrogenase_NAD", "N1", [("Ethanol", "default"), ("NAD", "km_NAD")], ["Acetaldehyde", "NADH"]),
        ("ALDH3", "dehydrogenase_NAD", "N2", [("Acetaldehyde", "default"), ("NAD", "km_NAD")], ["Acetate", "NADH"]),
        ("AtoDA2", "CoA_transferase", "N3", [("AcetoAcetylCoA", "default"), ("Acetate", "km_acceptor")], ["Acetoacetate", "AcetylCoA"]),
        ("AdvB2", "decarboxylase", "N4", [("Acetoacetate", "default")], ["Acetone", "CO2"]),
    ],
    cofactors={},
    module="Fatty acid/ketone-adjacent metabolism",
    notes="Ethanol oxidized to acetate first, then follows the same acetone route as Acetate_2_Acetone -- again reusing native thiolase (R5). Direct comparison point against the Acetate-fed version of the same product.",
    msp=(0.8, 1.5))

print(f"Total pathways defined: {len(PATHWAYS)}")
