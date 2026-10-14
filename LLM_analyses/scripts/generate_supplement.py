"""
Builds the supplementary-methods Word document summarizing the crude-lysate
cell-free malate/rTCA optimization work and the 30-pathway FRENDA-parameterized
substrate-diversification study, for use as a paper supplement.

Usage: ./.venv_sim/Scripts/python.exe scripts/generate_supplement.py
Output: results/Supplementary_Methods_and_Results.docx
"""
import os
from docx import Document
from docx.shared import Pt, Inches, RGBColor
from docx.enum.text import WD_ALIGN_PARAGRAPH
from docx.enum.table import WD_TABLE_ALIGNMENT
from docx.oxml.ns import qn
from docx.oxml import OxmlElement

ROOT = r"C:\Users\mkcoo\Documents\UW\LLM_CFS_2"
OUT_PATH = os.path.join(ROOT, "results", "Supplementary_Methods_and_Results.docx")

doc = Document()

# ---------------------------------------------------------------- styles
normal = doc.styles["Normal"]
normal.font.name = "Calibri"
normal.font.size = Pt(10.5)

for i in range(1, 4):
    doc.styles[f"Heading {i}"].font.color.rgb = RGBColor(0x1F, 0x3B, 0x57)

CODE_FONT = "Consolas"


def add_title(text):
    p = doc.add_paragraph()
    p.alignment = WD_ALIGN_PARAGRAPH.CENTER
    r = p.add_run(text)
    r.bold = True
    r.font.size = Pt(18)
    r.font.color.rgb = RGBColor(0x1F, 0x3B, 0x57)


def add_subtitle(text):
    p = doc.add_paragraph()
    p.alignment = WD_ALIGN_PARAGRAPH.CENTER
    r = p.add_run(text)
    r.italic = True
    r.font.size = Pt(11)
    r.font.color.rgb = RGBColor(0x55, 0x55, 0x55)


def h1(text):
    doc.add_heading(text, level=1)


def h2(text):
    doc.add_heading(text, level=2)


def p(text, bold=False, italic=False, size=10.5):
    para = doc.add_paragraph()
    r = para.add_run(text)
    r.bold = bold
    r.italic = italic
    r.font.size = Pt(size)
    return para


def bullet(text):
    para = doc.add_paragraph(style="List Bullet")
    para.add_run(text)
    return para


def code_block(text):
    para = doc.add_paragraph()
    para.paragraph_format.left_indent = Inches(0.3)
    shade = OxmlElement("w:shd")
    shade.set(qn("w:val"), "clear")
    shade.set(qn("w:fill"), "F2F2F2")
    para.paragraph_format.element.get_or_add_pPr().append(shade)
    for i, line in enumerate(text.strip("\n").split("\n")):
        if i > 0:
            para.add_run().add_break()
        r = para.add_run(line)
        r.font.name = CODE_FONT
        r.font.size = Pt(9)
    return para


def script_ref(text):
    para = doc.add_paragraph()
    r = para.add_run("Script(s): ")
    r.bold = True
    r.font.size = Pt(9.5)
    r2 = para.add_run(text)
    r2.font.name = CODE_FONT
    r2.font.size = Pt(9.5)
    r2.font.color.rgb = RGBColor(0x2E, 0x5C, 0x8A)
    return para


def add_figure(path, caption, width=6.0):
    doc.add_picture(path, width=Inches(width))
    last_p = doc.paragraphs[-1]
    last_p.alignment = WD_ALIGN_PARAGRAPH.CENTER
    cap = doc.add_paragraph()
    cap.alignment = WD_ALIGN_PARAGRAPH.CENTER
    r = cap.add_run(caption)
    r.italic = True
    r.font.size = Pt(9)
    r.font.color.rgb = RGBColor(0x55, 0x55, 0x55)
    doc.add_paragraph()


def make_table(headers, rows, col_widths=None, font_size=9):
    table = doc.add_table(rows=1, cols=len(headers))
    table.style = "Light Grid Accent 1"
    table.alignment = WD_TABLE_ALIGNMENT.CENTER
    hdr_cells = table.rows[0].cells
    for i, htext in enumerate(headers):
        hdr_cells[i].text = ""
        r = hdr_cells[i].paragraphs[0].add_run(htext)
        r.bold = True
        r.font.size = Pt(font_size)
    for row in rows:
        cells = table.add_row().cells
        for i, val in enumerate(row):
            cells[i].text = ""
            r = cells[i].paragraphs[0].add_run(str(val))
            r.font.size = Pt(font_size)
    if col_widths:
        for i, w in enumerate(col_widths):
            for row in table.rows:
                row.cells[i].width = Inches(w)
    doc.add_paragraph()
    return table


# ================================================================== TITLE PAGE
add_title("Supplementary Methods and Results")
add_subtitle("Kinetic-Model-Guided Design of rTCA Malate Production and\n"
             "FRENDA-Parameterized Heterologous Pathway Diversification\n"
             "in an E. coli Crude-Lysate Cell-Free System")
doc.add_paragraph()
p("This supplement documents the computational workflow, analyses, and design decisions "
  "underlying (1) optimization of malate production via a reductive TCA (rTCA) heterologous "
  "pathway in an E. coli crude-lysate cell-free system, and (2) a subsequent exploration of "
  "30 additional heterologous pathway designs spanning 10 candidate feedstocks, parameterized "
  "with literature/database-derived enzyme kinetics via the FRENDA pipeline. All results were "
  "produced using a mechanistic ODE kinetic model (model/iMC057.txt, Antimony/SBML format) "
  "simulated with libRoadRunner. Every analysis referenced here is backed by a named, "
  "version-controlled Python script in the project repository so that all reported numbers are "
  "independently reproducible; script names are given at the end of each subsection.")
doc.add_page_break()

# ================================================================== S1
h1("S1. Overview and Model System")
p("The base kinetic model (model/iMC057.txt) is a Common-Modular-Rate-Law (Liebermeister & "
  "Klipp, 2006) ODE representation of central metabolism in an E. coli-derived crude lysate, "
  "built from proteomic abundance data (EC numbers + relative abundances; see "
  "CellFreeModeling/examples/iMC057.csv) and fit to experimental time-course data. Fitted "
  "kinetic parameters (Kcat, Km, Keq, and the relative-density scaling terms rel0-rel3) are "
  "loaded from parameters/final_params.pkl, a set of 100 independent parameter fits; fit index "
  "60 is used throughout for consistency with the project's prior analysis notebooks (loss=456, "
  "rank 27/100 by fit quality -- not the single best fit, index 95 at loss=23, but the "
  "established reference point for this body of work).")
p("The base model natively encodes a heterologous reductive-TCA (rTCA) pathway for malate "
  "production, comprising three heterologous enzymes -- pyruvate carboxylase (hEC6411, "
  "“PYC”), formate dehydrogenase (hEC11719, “FDH”), and serine "
  "deaminase/ammonia-lyase (hEC43117, “SDS”) -- plus plasmid-driven overexpression of "
  "the endogenous malate dehydrogenase (eEC11137, “MDH”) run in the malate-producing "
  "direction. The competing, carbon-inefficient route to the same product is native fumarase "
  "(eEC4212) run in reverse (TCA direction). Total heterologous expression capacity is fixed at "
  "150 μM of protein, split in the base recipe as PYC=50, SDS=50, FDH=25, "
  "MDH-overexpression-increment=25 μM.")
p("Two independent decision variables recur throughout this work and are defined precisely "
  "because they were each a source of early confusion in this project (see §S2): the "
  "lysate dilution factor D (bounds 10-400; D=10 is the least-dilute, most-concentrated setting "
  "physically achievable, D=200 is the actual base bioproduction recipe -- D=10 is NOT the base "
  "case), and the simulation duration (bounds 30 min-24 h; base case is 4 h = 14400 s, the "
  "model's native time unit being seconds).")

# ================================================================== S2
h1("S2. Model Access, Parameterization Convention, and Reproducibility")
p("tellurium's Python import fails in the project's default (Anaconda) environment: "
  "libroadrunner is compiled against numpy<2 while the environment ships numpy 2.x, causing a "
  "DLL initialization failure; separately, tellurium's optional rrplugins submodule fails "
  "independently even after a numpy downgrade. The workaround, used for every simulation in "
  "this study, is to bypass tellurium entirely and call the antimony and roadrunner packages "
  "it wraps directly, inside an isolated virtual environment pinned to numpy<2 (.venv_sim).")
code_block('''import antimony, roadrunner

def loada(path):
    antimony.clearPreviousLoads()
    code = antimony.loadAntimonyFile(path)
    if code < 0:
        raise RuntimeError(antimony.getLastError())
    module_name = antimony.getMainModuleName()
    sbml = antimony.getSBMLString(module_name)
    return roadrunner.RoadRunner(sbml)

r = loada("model/iMC057.txt")''')
p("Every species in model/iMC057.txt is computed from a formula wired into the model itself "
  "rather than a literal default value: endogenous enzymes as "
  "eEC_X = p_eEC_X × rel0 × REF_X × dilution_factor, heterologous enzymes as "
  "hEC_X = p_hEC_X × rel{1,2,3} × 0.01 × dilution_factor, and metabolites as "
  "metab = v_i × REF_i × dilution_factor. Critically, the raw Antimony parameter "
  "dilution_factor multiplies concentration directly (larger = more concentrated), so the "
  "user-facing dilution factor D relates to it by dilution_factor = 1/D. This inverse "
  "relationship, and the fact that the base bioproduction recipe uses D=200 (not the D=10 "
  "floor), were each independently mis-inferred earlier in this project before being confirmed "
  "directly against a reference notebook; both are called out explicitly here because they are "
  "the single most consequential parameterization details in the whole study (§S3).")
p("The validated, end-to-end reproduction protocol for the base case (reproduces the project's "
  "reference notebook to 8-9 significant figures) is:")
code_block('''import antimony, roadrunner, pickle

antimony.clearPreviousLoads()
antimony.loadAntimonyFile('model/iMC057.txt')
sbml = antimony.getSBMLString(antimony.getMainModuleName())
r = roadrunner.RoadRunner(sbml)

with open('parameters/final_params.pkl', 'rb') as f:
    trainparams = pickle.load(f)
for name, val in trainparams[60][2].items():
    r.model[name] = val

D = 200                       # base bioproduction scheme dilution factor
d = 1 / D
r.model['dilution_factor'] = d
r.model['p_hEC6411']  = 0.050 / (0.01 * d * r.rel1)   # PYC  = 50 uM
r.model['p_hEC11719'] = 0.025 / (0.01 * d * r.rel2)   # FDH  = 25 uM
r.model['p_hEC43117'] = 0.050 / (0.01 * d * r.rel3)   # SDS  = 50 uM
# eEC11137 (MDH) overexpression increment (25 uM) added on top of its
# native background; exogenous cofactor/substrate additions (ATP, AcetylCoA,
# HCO3, NADH, NADPH, Formate, L-Serine) solved the same way -- see
# scripts/sim2.py for the complete, validated formula set.

r.reset()
results = r.simulate(0, 4 * 60 * 60)   # 4 h, in seconds
print(results['[SMalate]'][0], results['[SMalate]'][-1])
# -> 0.001056 mM, 0.106795 mM  (base case: final malate ~0.107 mM)''')
script_ref("scripts/loader.py (loada helper); scripts/sim2.py (full Design/run_design API "
           "implementing every formula above, with cumulative rTCA-vs-TCA flux tracking).")

# ================================================================== S3
h1("S3. Malate/rTCA Pathway: Identifying the Strongest Production Levers")
p("Design space: five candidate interventions were simulated on top of the corrected base "
  "case (T0), building up from the single highest-leverage change (lysate dilution) toward a "
  "multi-intervention optimum. All values below are from a single, internally consistent "
  "simulation batch (results/malate_rtca_designs_v5.json).")

make_table(
    ["ID", "Intervention(s)", "Final malate (mM)", "vs. T0", "rTCA:TCA ratio", "rTCA fraction"],
    [
        ["T0", "Base case (D=200, standard recipe)", "0.107", "—", "2.57 × 10⁶", "99.9999%"],
        ["T1", "D: 200 → 10", "0.854", "+700%", "37.5:1", "97.4%"],
        ["T2", "T1 + CS/fumarase knockdown to 10%", "1.043", "+877%", "381:1", "99.7%"],
        ["T3*", "T2 + serine 1→10 mM", "1.693", "+1485%", "717:1", "99.86%"],
        ["T4", "As T3 but D=25 (not 10)", "1.119", "+947%", "23,886:1", "99.996%"],
        ["T5", "D=10, CS/fumarase to 1%, serine→50 mM", "2.835", "+2555%", "16,741:1", "99.994%"],
    ],
    col_widths=[0.4, 2.0, 1.0, 0.7, 1.0, 1.0],
)
p("*T3 (“balanced”) was the top overall recommendation from this batch: best combination "
  "of absolute yield, selectivity, and cost/complexity. T4 trades 34% of T3's yield for a 33x "
  "better selectivity ratio via a more conservative dilution choice, illustrating that "
  "selectivity is extremely sensitive to dilution specifically (not to the other three "
  "interventions, which were held fixed between T3 and T4). T5 is the ceiling of the tested "
  "design space but is dominated by T4 on selectivity despite far greater cost.")

h2("S3.1 Lysate Dilution: The Single Highest-Leverage Variable")
p("Reducing dilution from the base D=200 to the D=10 floor alone (T1) increases yield "
  "~8-fold, the largest single-intervention effect found anywhere in this study. The mechanism "
  "was resolved by classifying all 54 native reactions into biosynthetic modules and fitting "
  "each module's instantaneous flux to a power law, flux ∝ D−α, across a pure "
  "native-lysate (no exogenous additions) dilution sweep (D=10 to D=400):")
make_table(
    ["Module", "α (dilution exponent)", "Interpretation"],
    [
        ["Glycolysis/gluconeogenesis", "1.42", "near-linear; substrate near/above Km, protected"],
        ["CO2/bicarbonate interconversion", "1.91", "intermediate"],
        ["TCA cycle", "2.01", "quadratic; enzyme AND substrate co-dilute"],
        ["Pentose phosphate pathway", "2.37", "quadratic-plus"],
        ["Acetate metabolism", "2.83", "steep"],
        ["One-carbon/serine metabolism", "3.15", "steepest; most dilution-sensitive module found"],
    ],
    col_widths=[2.2, 1.5, 3.0],
)
p("α=1 is expected if only enzyme concentration dilutes; α=2 if a limiting substrate "
  "co-dilutes with the enzyme (Michaelis-Menten term v ≈ kcat·E·S/Km when S≪Km, "
  "so v ∝ D−2 since both E and S individually scale as 1/D). The native competing route "
  "to malate (fumarase/TCA, α≈2) is suppressed far faster than the engineered rTCA route "
  "as dilution increases, which is the mechanistic reason selectivity improves dramatically at "
  "high dilution (base case D=200: rTCA fraction 99.9999%) even though absolute engineered "
  "flux is also suppressed — concentrating the lysate (lowering D) recovers absolute yield "
  "at a real but bounded selectivity cost (T1: still 97.4% rTCA at D=10).")
add_figure(os.path.join(ROOT, "results", "analysis", "11_module_dilution_exponent.png"),
           "Figure S1. Fitted dilution exponent α per biosynthetic module (native lysate, "
           "no exogenous additions, D=10-400). Modules with α≈2 (TCA cycle, pyruvate "
           "metabolism) lose flux far faster under dilution than glycolysis (α=1.42), "
           "explaining why increasing dilution improves rTCA:TCA selectivity.")
script_ref("scripts/module_dilution_analysis.py -> results/analysis/module_dilution_summary.txt, "
           "10_module_flux_share_vs_dilution.png, 11_module_dilution_exponent.png, "
           "12_alpha_vs_substrate_saturation.png.")

h2("S3.2 Serine Addition: The Largest Reagent-Cost Lever")
p("L-Serine is the pathway's carbon-entry substrate (consumed by hEC43117/SDS). A full dose-"
  "response sweep (1-100 mM) at the D=10 floor (results/analysis/serine_dilution_sweep.csv) "
  "shows a classic saturating response: final malate rises from 0.854 mM (1 mM serine, the base "
  "recipe's addition) to 1.479 mM at 10 mM, 2.011 mM at 20 mM, and plateaus near 2.46-2.50 mM "
  "beyond ~50-75 mM, with essentially no further gain from 75→100 mM. Marginal yield-on-"
  "serine (mol malate produced per mol serine actually consumed) falls monotonically from 9.66 "
  "at the 1 mM dose to 0.46 at 100 mM, quantifying the diminishing-returns point directly: doses "
  "above ~25-35 mM buy little additional titer per unit of (real, scaling) reagent cost. This "
  "dose-response was found to hold in the same qualitative shape across every dilution tested "
  "(D=10, 25, 50, 100, 200), with the saturating plateau value itself falling as the lysate is "
  "diluted (e.g. ~1.99 mM ceiling at D=100 vs. ~2.50 mM at D=10), confirming dilution and "
  "serine dose act as largely independent, additive levers rather than substituting for one "
  "another.")
add_figure(os.path.join(ROOT, "results", "analysis", "07_pareto_malate_vs_ratio_by_dilution.png"),
           "Figure S2. Final malate titer vs. rTCA:TCA selectivity ratio across the "
           "serine dose sweep (1-100 mM), overlaid for 5 dilution factors. Each curve traces "
           "the yield/selectivity tradeoff frontier at a fixed dilution as serine dose "
           "increases.")
script_ref("scripts/pareto_serine_dilution.py -> results/analysis/serine_dilution_sweep.csv, "
           "07_pareto_malate_vs_ratio_by_dilution.png, "
           "08_pareto_serine_yield_vs_ratio_by_dilution.png.")

h2("S3.3 Endogenous Knockdowns (Citrate Synthase, Fumarase)")
p("Knocking down endogenous citrate synthase (eEC2331) and fumarase (eEC4212) to 10% of their "
  "native expression level (simulating hydroxycitrate addition and direct genetic "
  "downregulation, respectively -- hydroxycitrate's own modeled concentration has no "
  "appreciable direct effect, see §S6 caveat below) is a zero-reagent-cost, low-risk lever "
  "that reliably improves selectivity (T1→T2: ratio 37.5→381, a 10× improvement) with a "
  "modest yield co-benefit (+877% vs. T0 at T2, vs. +700% at T1 alone).")

script_ref("scripts/run_designs5.py (reproduces the T0-T5 batch end to end) -> "
           "results/malate_rtca_designs_v5.json.")

doc.add_page_break()

# ================================================================== S4
h1("S4. Substrate Screening for Pathway Diversification")
p("To extend the platform beyond malate, 10 candidate starting substrates were evaluated as "
  "entry points for new heterologous pathways -- a mix of substrates already present in the "
  "model (Formate, Acetate, Pyruvate, Bicarbonate/CO2, S-Lactate) and substrates entirely "
  "absent/orthogonal to it (Glycerol, D-Xylose, Methanol, Glucose, Ethanol). Each was scored on "
  "estimated market feedstock price (MSP), classification/functional group, a qualitative "
  "sustainability (“greenness”) score (1-10), and native metabolic connectivity -- the "
  "number of endogenous model reactions already touching that substrate, and their "
  "instantaneous activity at the undiluted (D=10), no-exogenous-addition baseline.")
make_table(
    ["Substrate", "Est. MSP ($/kg)", "Greenness", "Native rxns"],
    [
        ["Formate", "0.4-0.8", "8", "0"],
        ["Acetate", "0.4-0.6", "7", "2"],
        ["Pyruvate", "3-8", "5", "6"],
        ["Bicarbonate/CO2", "0.1-0.3", "10", "2"],
        ["S-Lactate", "1.5-3", "7", "1"],
        ["Glycerol", "0.2-0.4", "9", "0"],
        ["D-Xylose", "1.0-3.0", "8", "0"],
        ["Methanol", "0.3-0.5", "6", "0"],
        ["Glucose", "0.3-0.6", "4", "0"],
        ["Ethanol", "0.5-0.8", "6", "0"],
    ],
    col_widths=[1.6, 1.6, 1.2, 1.2],
)
p("These 10 substrates were carried forward unchanged into pathway design (§S5); their "
  "native-connectivity values re-emerge as an important, initially counter-intuitive predictor "
  "of pathway yield in §S7 once downstream kinetic-modeling artifacts were corrected.")
script_ref("Compiled via direct simulation queries against model/iMC057.txt "
           "-> results/substrate_candidates.csv (no standalone generation script; ad hoc "
           "analysis session, cross-checked against a dedicated undiluted/no-heterologous-"
           "enzyme control run).")

# ================================================================== S5
h1("S5. Design and FRENDA-Based Kinetic Parameterization of 30 Heterologous Pathways")
p("For each of the 10 substrates, three distinct heterologous pathways to different products "
  "were designed (30 pathways total), subject to: (i) a shared 150 μM total heterologous "
  "expression budget per pathway, split evenly across that pathway's enzyme count; (ii) "
  "exclusion of membrane-associated enzyme isoforms (e.g. periplasmic PQQ-dependent methanol "
  "dehydrogenase, membrane-bound glycerol-3-phosphate dehydrogenase), since these are known to "
  "be difficult to functionally express in a soluble lysate system -- soluble/cytoplasmic "
  "alternatives were selected instead (e.g. Bacillus methanolicus-type NAD-dependent methanol "
  "dehydrogenase, Saccharomyces cerevisiae GPD1); and (iii) explicit listing of any newly "
  "required exogenous substrate/cofactor additions.")
p("Kinetic parameters for each new enzyme were obtained via FRENDA (Fast Retrieval of ENzyme "
  "DAta; CellFreeModeling/), a pipeline that maps EC numbers to KEGG reactions and retrieves "
  "experimentally measured Km/Kcat from BRENDA where available, falling back to CatPred "
  "(ML-based prediction from protein sequence + substrate SMILES) otherwise, with reaction "
  "equilibrium constants (Keq) computed via eQuilibrator. 41 unique EC numbers were submitted "
  "(covering all enzymes across all 30 designs); FRENDA successfully resolved 39 (2 EC numbers "
  "-- formate dehydrogenase and PHA synthase -- have no entry in KEGG's enzyme database under "
  "that number and were parameterized manually instead, see §S6).")
p("Because FRENDA's auto-retrieved KEGG reaction for a given EC does not always match the "
  "exact substrate/product role structure of the hand-designed pathway step it was applied to "
  "(e.g. a different cofactor specificity, or the reverse reaction direction), a hybrid "
  "parameterization strategy was used: the real, FRENDA/BRENDA/CatPred-sourced Kcat, Km, and "
  "Keq values were applied to the project's own hand-designed reaction topology (same "
  "Common-Modular-Rate-Law Antimony macro-expansion style as the base model, see §S1), with "
  "Km values assigned positionally (substrate-role-by-substrate-role) rather than by strict "
  "chemical-identity matching. This preserved full control over pathway connectivity to native "
  "metabolite pools while still grounding every kinetic parameter in a real, enzyme/EC-specific "
  "source rather than a generic literature guess. The consequences and correction of this "
  "approximation are the subject of §S6.")
script_ref("scripts/pathway_designs.py (all 30 pathway topologies/substrates/cofactors/"
           "products/module classifications/rationale); scripts/ec_map.py (enzyme-name -> EC "
           "mapping, 41 ECs); scripts/frenda_kinetics.py (parses FRENDA output into a "
           "Kcat/Km/Keq lookup); scripts/pathway_generator.py (splices FRENDA-parameterized "
           "reactions into a copy of model/iMC057.txt per pathway -> model/pathways/*.txt); "
           "scripts/run_pathway_sims.py (simulates all 30, D=10, 4 h); "
           "scripts/compile_results_table.py (-> results/pathway_designs_table.csv). Raw FRENDA "
           "output preserved at parameters/pathway_frenda_reactions.csv, "
           "pathway_frenda_species.csv, pathway_frenda_proteome_input.csv.")

doc.add_page_break()

# ================================================================== S6
h1("S6. Diagnosis and Correction of FRENDA-Derived Kinetic Modeling Artifacts")
p("An initial simulation pass of all 30 pathways (D=10, 4 h, each design's own substrate dose "
  "and cofactor additions) produced uniformly low yields (peak yield, defined as 100 × "
  "Δproduct / substrate dose, mol product per mol substrate dosed): median well under 1%, "
  "18 of 30 designs at or near 0%. Rather than accept this at face value, each near-zero design "
  "was diagnosed programmatically -- comparing each new reaction's initial forward/reverse flux, "
  "Keq, and cofactor/substrate concentration trajectories against FRENDA's raw fetched "
  "Substrates/Products strings -- revealing five distinct, largely non-overlapping root causes:")
make_table(
    ["#", "Root cause", "Nature", "Designs affected"],
    [
        ["1", "Keq applied backwards (FRENDA's fetched reaction ran opposite to the modeled "
              "step; e.g. alanine dehydrogenase, methylenetetrahydrofolate dehydrogenase, PEP "
              "carboxylase, aspartate dehydrogenase)", "Modeling artifact",
         "4 (Formate_3_Alanine, Lactate_2_Alanine, HCO3_3_Aspartate, +partial Formate_1_Serine)"],
        ["2", "Acetolactate synthase (EC 2.2.1.6) fetched a non-representative KEGG fragment "
              "reaction (a TPP-intermediate half-reaction) instead of the real overall "
              "2 Pyruvate -> Acetolactate + CO2 conversion", "Data-selection artifact",
         "5 (every 2,3-BDO/isobutanol-family design)"],
        ["3", "Two unrelated ECs (acetate CoA-transferase, acyl-CoA hydrolase) returned the "
              "identical, chemically implausible Keq = 7.12 x 10−218", "Computation-failure "
              "artifact", "3 (Acetate_2_Acetone, Ethanol_3_Acetone, Lactate_3_AcrylicAcid)"],
        ["4", "Methanol oxidation (NAD-dependent methanol dehydrogenase) is genuinely "
              "thermodynamically uphill; FRENDA's fetched direction matched the modeled step "
              "and its Keq (4.5 x 10−6) is consistent with published values", "Real biology, "
              "not a bug", "3 (all Methanol-fed designs)"],
        ["5", "One-time cofactor bolus (ATP, NADH) depleted by the broader native reaction "
              "network well before the entry step could sustain flux; directly confirmed via "
              "concentration trajectories", "Real engineering constraint",
         "Compounds several designs, most cleanly Glycerol_2_BDO/Xylose_3_BDO/Glucose_1_BDO"],
    ],
    col_widths=[0.3, 2.6, 1.3, 1.9],
)
p("Causes 1-3 were corrected in scripts/frenda_kinetics.py: for the four reversed-direction "
  "ECs, Keq was inverted (1/Keq) and the substrate/product Km-role lists swapped to match; for "
  "acetolactate synthase and the two corrupted-Keq ECs, FRENDA's returned values were replaced "
  "outright with literature-typical values for the real overall reaction. Causes 4 and 5 were "
  "left uncorrected, as they reflect genuine chemistry/engineering constraints rather than "
  "software or data-retrieval defects; cause 5 was independently confirmed by re-simulating "
  "one affected design (Glucose_1_BDO) after the cause-2 fix and observing ATP still crash to "
  "~0 mM with native pyruvate never accumulating -- i.e. the fixed enzyme simply had nothing to "
  "act on.")
p("Re-simulating all 30 designs after the causes-1-3 fixes produced substantial, selective "
  "improvement exactly where predicted:")
make_table(
    ["Design", "Peak yield before fix", "Peak yield after fix", "Cause fixed"],
    [
        ["Pyruvate_2_Isobutanol", "0.29%", "17.85% (best design overall)", "#2 (ALS)"],
        ["Pyruvate_3_Alanine", "0.56%", "15.00%", "#1 (AlaDH direction)"],
        ["Lactate_3_AcrylicAcid", "0.00%", "5.86%", "#3 (corrupted Keq)"],
        ["Pyruvate_1_BDO", "0.02%", "2.30%", "#2 (ALS)"],
        ["HCO3_3_Aspartate", "0.00%", "0.28%", "#1 (PPC + AspDH direction)"],
        ["Formate_3_Alanine", "0.00%", "0.08%", "#1 (AlaDH direction)"],
    ],
    col_widths=[1.8, 1.4, 2.2, 1.3],
)
p("Designs whose bottleneck was cause 4 or 5 (not addressed by these fixes) remained at or "
  "near their pre-fix values, as expected -- a useful internal check that the diagnosis "
  "correctly separated fixable software/data artifacts from genuine biological/engineering "
  "limits.")
script_ref("scripts/diagnose_failures.py (-> results/pathway_analysis/failure_diagnosis.json, "
           "failure_root_causes.md); fixes implemented in scripts/frenda_kinetics.py "
           "(REVERSED_DIRECTION_ECS, FORCED_OVERRIDE); results regenerated via "
           "scripts/run_pathway_sims.py and scripts/compile_results_table.py.")

# ================================================================== S7
h1("S7. Design-Space Analysis of the 30 Pathways")
p("Following the kinetics correction, the 30-design results table was merged with the "
  "substrate-screening data (§S4) and analyzed for predictors of yield: number of new "
  "heterologous enzymes (pathway length), biosynthetic module recruited, starting substrate "
  "identity, and the substrate's native metabolic connectivity (endogenous reaction count).")
p("Pathway length was NOT a significant predictor of yield after the kinetics fix (Pearson "
  "r=+0.17, p=0.38), a substantial reversal from the pre-fix analysis (r=−0.45, p=0.014), "
  "indicating that most of the original apparent “shorter pathways perform better” "
  "signal was actually “pathways that happened to avoid the broken acetolactate-synthase or "
  "alanine-dehydrogenase steps perform better.” Substrate native connectivity became the "
  "strongest single correlate of yield post-fix (Pearson r=0.71 for mean yield vs. endogenous "
  "reaction count), driven substantially by pyruvate (6 native reactions, the most-connected "
  "substrate screened) whose two pyruvate-fed designs went from among the weakest performers "
  "pre-fix to the two best performers in the entire 30-design set post-fix -- though this "
  "correlation is concentrated in that one substrate as a leverage point (Spearman "
  "ρ=0.29, not significant) rather than a smooth trend across all 10 substrates.")
p("Methods: module category was assigned by keyword classification of each design's free-text "
  "module description into 7 buckets (C1/serine, TCA/anaplerosis, fatty-acid/ketone, amino-acid "
  "biosynthesis, pentose phosphate pathway, glycolysis, pyruvate metabolism); correlation "
  "statistics (Pearson and Spearman) were computed with scipy.stats. Five summary figures were "
  "produced: yield vs. enzyme count (scatter, regression line), yield by substrate (bar + "
  "strip), yield by module (bar + strip), yield vs. substrate connectivity (scatter, "
  "regression line), and yield by native-step-reuse (box + strip).")
add_figure(os.path.join(ROOT, "results", "pathway_analysis", "01_yield_vs_n_enzymes.png"),
           "Figure S3. Peak yield vs. number of new heterologous enzymes, post-kinetics-fix "
           "(Pearson r=+0.17, p=0.38 -- not significant).", width=5.5)
add_figure(os.path.join(ROOT, "results", "pathway_analysis", "04_yield_vs_substrate_interconnectedness.png"),
           "Figure S4. Mean peak yield (per substrate, across its 3 pathways) vs. the "
           "substrate's native endogenous-reaction count, post-kinetics-fix (Pearson r=0.71, "
           "p<0.001; driven substantially by the pyruvate outlier).", width=5.5)
script_ref("scripts/analyze_pathway_space.py -> results/pathway_analysis/"
           "{01_yield_vs_n_enzymes, 02_yield_by_substrate, 03_yield_by_module, "
           "04_yield_vs_substrate_interconnectedness, 05_yield_by_reuses_native}.png, "
           "correlation_summary.txt, top_bottom_designs.txt, merged_pathway_analysis_data.csv.")

h2("S7.1 Summary Figure: The 30-Pathway Design Map")
p("A single summary figure maps all 30 designs simultaneously: candidate substrates positioned "
  "vertically by ascending market feedstock cost, product nodes on the right (marker shape "
  "distinguishing products already native to the base model from fully novel ones), one curved "
  "line per pathway connecting substrate to product, dot markers along each line indicating the "
  "number of heterologous enzymes in that step, line thickness (square-root scaled, so "
  "near-zero-yield pathways remain visible rather than disappearing) encoding peak yield, and "
  "line color encoding the recruited biosynthetic module.")
script_ref("scripts/plot_pathway_map.py -> results/pathway_analysis/06_pathway_map.png.")

add_figure(os.path.join(ROOT, "results", "pathway_analysis", "06_pathway_map.png"),
           "Figure S5. All 30 pathway designs: substrate (left, ordered by ascending "
           "feedstock cost) to product (right, diamond=already native to the base model, "
           "circle=novel), one line per pathway, dots=enzyme count, line thickness=peak "
           "yield, color=biosynthetic module.", width=6.3)

doc.add_page_break()

# ================================================================== S8
h1("S8. Dilution- and Substrate-Loading Optimization of the 30 Pathways")
p("For each of the 30 kinetics-corrected pathway designs, a two-dimensional grid search over "
  "lysate dilution factor D (7 log-spaced points, 10-400, the project's own physical bounds) "
  "and substrate loading dose (7 log-spaced points, 1-100 mM) was performed to maximize peak "
  "product titer (mM) -- deliberately not yield percent, since yield% = 100 × Δproduct / "
  "dose has dose in the denominator and would spuriously reward vanishingly small, impractical "
  "doses. Cofactor additions (ATP, NADH, NADPH, etc.) were held fixed at each design's "
  "originally specified value; only D and substrate dose were searched, consistent with the "
  "scope of this analysis. A coarse-resolution (120-point) simulation was used across the "
  "49-point grid per pathway for speed, with the identified optimum re-confirmed at full "
  "(400-point) resolution. A small number of grid corners (extreme D combined with high dose) "
  "were numerically stiff (CVODE integration failure); these were caught and skipped rather "
  "than allowed to abort the search.")
p("Across all 30 pathways, the median fold-improvement in peak titer was only 1.01× -- most "
  "designs were already close to their D=10/original-dose optimum, consistent with many being "
  "limited by total enzyme turnover capacity (bounded by the shared 150 μM budget divided "
  "across each pathway's enzyme count) rather than by substrate availability or dilution "
  "per se. Six of 30 pathways improved ≥ 2-fold:")
make_table(
    ["Design", "Baseline (D=10, orig. dose)", "Optimum found", "Fold"],
    [
        ["Formate_2_Pyruvate", "8.2 × 10−3 mM", "D=400, dose=100 mM → 0.163 mM", "19.8×"],
        ["Glycerol_2_BDO", "4 × 10−5 mM", "D=117, dose=100 mM → 1.8 × 10−4 mM", "4.6×"],
        ["Xylose_3_BDO", "4 × 10−5 mM", "D=117, dose=100 mM → 1.9 × 10−4 mM", "4.4×"],
        ["Glucose_1_BDO", "5 × 10−5 mM", "D=117, dose=100 mM → 1.9 × 10−4 mM", "3.6×"],
        ["Glycerol_3_DHA", "0.027 mM", "D=10, dose=100 mM → 0.075 mM", "2.8×"],
        ["Pyruvate_1_BDO", "0.230 mM", "D=10, dose=46 mM → 0.538 mM", "2.3×"],
    ],
    col_widths=[1.8, 1.9, 2.1, 0.8],
)
p("The Formate_2_Pyruvate result is notable in direction: its optimum lies at maximum dilution "
  "(D=400), the opposite of the malate pathway's own dominant lever (§S3.1) -- plausibly "
  "because dilution simultaneously suppresses the native reactions competing for the same "
  "pyruvate/formate pools. The three cofactor-depletion-limited designs (§S6, cause 5) "
  "improved a real but practically insignificant 3.6-4.6×, remaining far below any of the "
  "kinetics-fixed pathways in absolute titer -- confirming that D/dose tuning alone cannot "
  "rescue a design whose true bottleneck is a non-regenerating cofactor bolus; that would "
  "require a dosing-schedule or cofactor-regeneration-system redesign, outside this analysis's "
  "scope.")
p("The pre-optimization (“un-optimized”) model files and results were preserved for "
  "comparison rather than overwritten: model/pathways_unoptimized/, "
  "results/pathway_sim_results_unoptimized.json, results/pathway_designs_table_unoptimized.csv, "
  "results/pathway_designs_unoptimized.json, results/pathway_analysis_unoptimized/.")
script_ref("scripts/optimize_pathways.py (-> results/pathway_optimization.json, "
           "pathway_optimization_table.csv); scripts/plot_optimization_results.py -> "
           "results/pathway_analysis/07_optimization_improvement.png (log-scale dumbbell chart, "
           "one row per pathway sorted by fold-improvement, baseline vs. optimum, annotated "
           "with the (D, dose) that achieved it).")

add_figure(os.path.join(ROOT, "results", "pathway_analysis", "07_optimization_improvement.png"),
           "Figure S6. Baseline (D=10, original dose) vs. optimized (best-found D, dose) "
           "peak product titer for all 30 pathways, sorted by fold-improvement.", width=6.3)

p("Full time-course product-formation curves (un-optimized overlaid with optimized, all 30 "
  "pathways, laid out by substrate x descending optimized yield) are provided as Figure S7.")
add_figure(os.path.join(ROOT, "results", "pathway_analysis", "08_timecourses_before_after.png"),
           "Figure S7. Product concentration vs. time (0-4 h), un-optimized (gray dashed) vs. "
           "optimized (red solid), for all 30 pathways. Rows: substrate, ascending feedstock "
           "cost top-to-bottom. Columns: pathway, descending optimized peak yield left-to-"
           "right. Note independent y-axis scaling per panel (titers span ~4 orders of "
           "magnitude across the design set) and that ordering is by yield percent (molar "
           "efficiency), not absolute titer -- e.g. the HCO3 row's leftmost panel (ThreeHP) "
           "has a far smaller absolute titer than its neighbor (SMalate) but a marginally "
           "higher molar yield, since it required a 100x smaller substrate dose.", width=6.3)
script_ref("scripts/plot_timecourses.py -> results/pathway_analysis/"
           "08_timecourses_before_after.png.")

# ================================================================== S9
h1("S9. Limitations and Caveats")
bullet("Fitted kinetic parameters for the base model (parameters/final_params.pkl, index 60) "
       "are one of 100 independent fits to experimental data, not the single best fit by loss "
       "(index 95); index 60 was used throughout for consistency with prior project analyses.")
bullet("New-pathway Km values were assigned positionally (substrate-role order) rather than by "
       "strict chemical-identity matching against FRENDA's fetched reaction, because the "
       "hand-designed reaction topology does not always exactly match KEGG's first-listed "
       "reaction for a given EC. This is a defensible approximation for early-stage screening "
       "but not a substitute for reaction-specific experimental characterization.")
bullet("All 30 new pathway reactions use simplified 1:1 substrate:product stoichiometry, even "
       "where the real reaction is not 1:1 (e.g. acetolactate synthase formally consumes 2 "
       "pyruvate per acetolactate).")
bullet("Cofactor and substrate additions are modeled as a single, non-regenerating bolus at "
       "t=0, not a continuous or fed-batch supply; several designs (§S6, S8) are "
       "demonstrably limited by this modeling choice rather than by enzyme kinetics per se.")
bullet("PHB-producing designs track flux into a soluble “PHB_monomer” proxy species; "
       "the model cannot represent real polymer accumulation/precipitation.")
bullet("Enzymes with known oxygen-sensitivity in vivo (pyruvate formate-lyase, dihydroxy-acid "
       "dehydratase, lactoyl-CoA dehydratase) are modeled with the same kinetics regardless of "
       "aerobic/anaerobic lysate preparation conditions; real expression/activity risk for "
       "these steps is flagged in the per-design notes (results/pathway_designs_table.csv) but "
       "not mechanistically modeled.")
bullet("Two EC numbers (formate dehydrogenase, EC 1.17.1.9; PHA synthase, EC 2.3.1.264) have no "
       "entry in KEGG's enzyme database and could not be auto-parameterized by FRENDA; the "
       "former reused this project's own previously-fitted hEC11719 parameters (same enzyme, "
       "already present in the base model), the latter used literature-typical values.")
bullet("Hydroxycitrate's own modeled concentration has no appreciable direct effect in this "
       "model (its REF scaling constant is extremely small); its physiological inhibition of "
       "citrate synthase is instead simulated via direct modulation of citrate synthase "
       "(eEC2331) expression level, per project convention.")

# ================================================================== S10
h1("S10. Data and Code Availability")
p("All scripts referenced in this document reside under scripts/ in the project repository. "
  "Key input/output data files:")
make_table(
    ["File", "Contents"],
    [
        ["model/iMC057.txt", "Base kinetic model (Antimony format)"],
        ["parameters/final_params.pkl", "100 independent fitted parameter sets"],
        ["model/pathways/*.txt", "30 kinetics-corrected pathway models (post-§S6 fix)"],
        ["model/pathways_unoptimized/*.txt", "30 pathway models prior to D/dose optimization (§S8)"],
        ["parameters/pathway_frenda_*.csv", "Raw FRENDA output and proteomic input for the 30-pathway EC set"],
        ["results/malate_rtca_designs_v5.json", "T0-T5 malate/rTCA design batch (§S3)"],
        ["results/analysis/serine_dilution_sweep.csv", "Serine x dilution Pareto sweep (§S3.2)"],
        ["results/analysis/module_dilution_summary.txt", "Per-module dilution exponents (§S3.1)"],
        ["results/substrate_candidates.csv", "10-substrate screening table (§S4)"],
        ["results/pathway_designs_table.csv", "Final 30-pathway results table (post-fix, post-optimization inputs)"],
        ["results/pathway_analysis/failure_root_causes.md", "Full failure-mode diagnosis writeup (§S6)"],
        ["results/pathway_optimization_table.csv", "D/dose optimization results (§S8)"],
        ["results/pathway_analysis/*.png", "All summary and diagnostic figures"],
    ],
    col_widths=[2.4, 4.0],
    font_size=8.5,
)

doc.save(OUT_PATH)
print(f"Wrote {OUT_PATH}")
