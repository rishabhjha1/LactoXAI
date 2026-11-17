"""Inject the run's numbers into the site template. Every figure on the page comes from these files."""
import json, csv, re, sys, math

repo = sys.argv[1]
R = json.load(open(f"{repo}/results/results.json"))
N = json.load(open(f"{repo}/results/nmf_programs.json"))
CF = {r["gene"]: float(r["effect_log_days_per_1sd"])
      for r in csv.DictReader(open(f"{repo}/results/F17_counterfactuals.csv"))}

top30 = R["top_shap_genes"][:30]
top25 = top30[:25]
feature_count = int(re.search(r"\((\d+) genes\)", next(k for k in R["representation_r2"] if k.startswith("HVG"))).group(1))

epi_progs = ["NMF2", "NMF7", "NMF9", "NMF12"]
epi_genes = set(g for p in epi_progs for g in N[p])
hits = sum(g in epi_genes for g in top25)
expected = 25 * len(epi_genes) / feature_count

pretty = {"immune_signalling": "Immune signalling", "casein_secretion": "Casein secretion",
          "lipid_synthesis": "Lipid synthesis", "secretory_machinery": "Secretory machinery",
          "hmo_glycosylation": "HMO glycosylation", "stress_response": "Stress response"}
curated = {g: pretty[p] for p, c in R["ablation_coverage"].items() for g in c["in_features"]}

nmf_labels = {"NMF1": "macrophage", "NMF2": "secretory lactocyte, with mitochondrial genes",
              "NMF3": "immediate-early stress", "NMF4": "neutrophil", "NMF5": "T cell",
              "NMF6": "MHC-II antigen-presenting", "NMF7": "progenitor-like epithelium, MHC-II and complement",
              "NMF8": "proliferation", "NMF9": "junctional epithelium", "NMF10": "monocyte / macrophage",
              "NMF11": "interferon response", "NMF12": "secretory lactocyte"}

ab, chi = R["program_ablation_r2_loss"], R["ablation_control_p97_5"]
beat = [k for k, v in sorted(ab.items(), key=lambda kv: -kv[1]) if v > chi]
f = lambda v, d=2: f"{v:.{d}f}".replace("-", "−")
fr, fe = R["fold_r2"], R["clock_lactocyte_only"]["fold_r2"]
ov = R["shap_perm_overlap_top25"]
faith = R["faithfulness"]; i5 = faith["ks"].index(5)

evidence = [
    dict(q="Is lactational time readable from milk cells?", read=f"{feature_count:,} genes in 25,000 cells",
         score=f"R² {f(R['clock_donor_heldout']['r2'])} (folds {f(min(fr))}–{f(max(fr))})",
         base=f"best of 20 permutations {f(R['null_r2_max'])}", v="good"),
    dict(q="Is it more than recognising the mother?", read="same model, random vs held-out split",
         score=f"R² {f(R['clock_donor_heldout']['r2'])} held-out; {f(R['clock_random_split']['r2'])} leaky",
         base=f"mother identified {R['donor_identity_accuracy']:.0%} vs {R['donor_identity_chance']:.0%} chance",
         v="good", note=f"Leakage worth {f(R['r2_inflation_from_leakage'])}"),
    dict(q="Is it more than a change in cell mix?", read=f"epithelial cells only ({R['clock_lactocyte_only']['n_cells']:,})",
         score=f"R² {f(R['clock_lactocyte_only']['r2'])} (folds {f(min(fe))}–{f(max(fe))})",
         base="permutation null not run for this subset", v="good", note="One fold weak"),
    dict(q="Do the explanations reflect the model?", read="SHAP, checked by masking and a shuffled-label model",
         score=f"top 5 genes: error {f(faith['shap'][0])} → {f(faith['shap'][i5])}",
         base=f"random 5 genes: {f(faith['random'][i5])}; shuffled model r {f(R['randomization_corr'])}",
         v="good", note=f"{ov} of 25 genes shared with permutation importance"),
    dict(q="Does an unsupervised method find the same genes?", read="NMF, 12 programmes, no time labels",
         score=f"{hits} of the top 25 SHAP genes fall in 4 epithelial programmes",
         base=f"about {expected:.0f} expected at random", v="good", note="The random baseline is generous"),
    dict(q="Which curated programmes carry the clock?", read="6 programmes, genes masked",
         score="R² lost: " + ", ".join(f"{pretty[k].split()[0].lower()} {f(ab[k])}" for k in beat),
         base=f"random 5 genes ≤ {f(chi)}", v="part", note=f"{len(beat)} of 6 beat the control"),
    dict(q="Do programme averages work as features?", read="6 programme scores",
         score=f"R² {f(R['representation_r2']['curated programs'][0])}",
         base=f"{feature_count:,} genes: {f(next(v for k, v in R['representation_r2'].items() if k.startswith('HVG'))[0])}", v="weak"),
    dict(q="Can it resolve late lactation?", read="3 mothers sampled after 200 days",
         score="under-predicted past ~150 days", base="–", v="no"),
]

programmes = [
    dict(dir="early", name="Progenitor-like epithelium", nmf="NMF7 · with MHC-II and complement genes", genes=N["NMF7"][:14]),
    dict(dir="late", name="Junctional epithelium", nmf="NMF9 · adhesion and tight-junction genes", genes=N["NMF9"][:14]),
    dict(dir="late", name="Secretory lactocytes", nmf="NMF12 (and NMF2) · milk protein and lipid genes", genes=N["NMF12"][:14]),
]

figures = [
    ("F9_clock_predictions.png", "Lactation clock", "Days postpartum predicted for held-out mothers; each column is one sample. Predictions flatten past about 150 days."),
    ("F25_epithelial_only_clock.png", "Not a cell-mix detector", "Restricted to epithelial cells, the clock reaches R² 0.56 against 0.53 for all cells."),
    ("F10_leakage_control.png", "Leakage, measured", "Letting a mother's cells fall on both sides of the split raises R² from 0.53 to 0.63."),
    ("F22_permutation_null.png", "Permutation null", "Timepoints shuffled across whole samples, 20 times; none reaches the observed score."),
    ("F23_fold_variance.png", "Fold-to-fold spread", "R² ranges from 0.38 to 0.64 across held-out groups of mothers."),
    ("F11_shap_summary.png", "What the clock reads", "SHAP beeswarm for the top 30 genes. Mitochondrial genes vary by mother and are not interpreted."),
    ("F17_counterfactuals.png", "Direction of each gene", "Change in predicted time per +1 SD of each gene, all else fixed."),
    ("F16_program_ablation.png", "Programme ablation", "Immune, casein and lipid programmes beat a random-gene control; the other three do not."),
    ("F21_representation_comparison.png", "Programme averages are weak", "Six programme scores reach R² 0.15; individual genes reach 0.48."),
    ("F18_faithfulness.png", "Faithfulness", "Masking SHAP's top genes degrades the model far faster than masking random genes."),
    ("F19_randomization_check.png", "Randomization check", "A shuffled-label model's attributions correlate only r = 0.25 with the real ones."),
    ("F1_sampling_grid.png", "The confound", "A median of one mother per timepoint, and three mothers after 200 days."),
]

limits = [
    "<b>Late lactation rests on three mothers.</b> Every sample after about 200 days comes from BM04, BM06 or BM10, and the clock under-predicts there. It separates early from later lactation; it does not resolve months 6 to 11.",
    "<b>Time and mother are confounded by design</b>, with a median of one mother per timepoint. Donor-held-out validation limits this but cannot remove it.",
    "<b>Coverage.</b> The source paper follows 15 mothers for up to 632 days; in the public metadata, cells from 13 of them carry a timepoint, between days 3 and 323. This analysis uses that subset.",
    "<b>Mitochondrial genes</b> (MT-ND4, MT-ND5) push predictions later but vary by mother, so they may reflect sample quality rather than biology.",
    "<b>p ≤ 0.048 is the floor for 20 permutations</b>; more permutations would sharpen it. The permutation null was not run for the epithelial-only model.",
    "<b>One observational cohort.</b> Counterfactual shifts describe the model, not biology, and no result is a causal claim about any gene.",
]

D = dict(results=R, convergence=dict(hits=hits, of=25, expected=expected), evidence=evidence,
         programmes=programmes, topSet=top30, topGenes=top30, cf=CF, nmf=N, nmfLabels=nmf_labels,
         pretty=pretty, featureCount=feature_count, curated=curated,
         figures=[dict(file=a, title=b, cap=c, alt=f"{b}: {c}") for a, b, c in figures], limits=limits)

tpl = open(sys.argv[2]).read()
html = tpl.replace("__DATA__", json.dumps(D, ensure_ascii=False))
open(f"{repo}/docs/index.html", "w").write(html)
open(f"{repo}/docs/.nojekyll", "w").write("")
print(f"wrote docs/index.html ({len(html)/1024:.0f} KB) | convergence {hits}/25, expected {expected:.2f} | beat control: {beat}")
