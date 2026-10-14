"""
KEGG compound ID -> existing native species name in model/iMC057.txt.

Used to splice FRENDA/ODBM-generated reactions (which name metabolites by
raw KEGG ID, e.g. "C00002") into the base crude-lysate model so that shared
metabolites (ATP, NADH, Pyruvate, etc.) correctly connect to the SAME pools
already tracked natively, rather than creating disconnected duplicate species.

Verified against FRENDA's actual output where noted (glycerol kinase test:
C00002=ATP, C00008=ADP, C00116=Glycerol, C00093=Glycerol-3-phosphate;
formate-THF ligase test: C00143=5,10-Methylenetetrahydrofolate).
"""

KEGG_TO_NATIVE = {
    "C00001": "H2O",
    "C00002": "ATP",
    "C00003": "NAD",
    "C00004": "NADH",
    "C00005": "NADPH",
    "C00006": "NADP",
    "C00008": "ADP",
    "C00009": "Orthophosphate",
    "C00010": "CoA",
    "C00011": "CO2",
    "C00014": "Ammonia",
    "C00020": "AMP",
    "C00022": "Pyruvate",
    "C00024": "AcetylCoA",
    "C00033": "Acetate",
    "C00036": "Oxaloacetate",
    "C00037": "Glycine",
    "C00058": "Formate",
    "C00065": "L_Serine",
    "C00074": "Phosphoenolpyruvate",
    "C00091": "SuccinylCoA",
    "C00122": "Fumarate",
    "C00143": "_5_10_Methylenetetrahydrofolate",
    "C00149": "SMalate",
    "C00186": "S_Lactate",
    "C00092": "DGlucose6phosphate",
    "C00288": "HCO3",
    "C00332": "AcetoAcetylCoA",
    "C00083": "Malonyl_CoA",
    "C00111": "Glycerone_phosphate",  # DHAP
    "C00101": "Tetrahydrofolate",
    "C00354": "D_Fructose_1_6_bisphosphate",
    "C00085": "D_Fructose_6_phosphate",
    "C00117": "D_ribose_5_phosphate",
    "C00199": "D_Xylulose_5_phosphate",
    "C00231": "D_Xylulose_5_phosphate",
    "C00309": "D_Ribulose_5_phosphate",
    # deliberately NOT mapped (genuinely new/orthogonal starting substrates --
    # these get their own new native-like species, no existing counterpart):
    #   C00116 Glycerol, C00181 D-Xylose, C00132 Methanol, C00031 Glucose,
    #   C00469 Ethanol
}
