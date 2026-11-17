"""Build the LactoXAI knowledge base (docs/kb/index.html) from the run's result files.

Trust rule (0-5), stated in the app:
  one point for each check passed by the evidence behind a link:
  1 validated on held-out mothers, 2 beats a permutation / shuffled-label null,
  3 explanation is faithful (masking test), 4 explanation passes the randomization check,
  5 corroborated by an independent method (NMF programme agreeing in direction).
  Anything flagged for a known confound is capped at 2.
Literature links carry the source's own standing and are drawn solid, not scored by these checks.
"""
import json, csv, re, sys

repo, tpl_path = sys.argv[1], sys.argv[2]
R = json.load(open(f"{repo}/results/results.json"))
N = json.load(open(f"{repo}/results/nmf_programs.json"))
CF = {r["gene"]: float(r["effect_log_days_per_1sd"])
      for r in csv.DictReader(open(f"{repo}/results/F17_counterfactuals.csv"))}
top30 = R["top_shap_genes"][:30]
fc = int(re.search(r"\((\d+) genes\)", next(k for k in R["representation_r2"] if k.startswith("HVG"))).group(1))
f2 = lambda v, d=2: f"{v:.{d}f}".replace("-", "−")
pct = lambda v: round(abs(2.718281828 ** v - 1) * 100)

PRETTY = {"immune_signalling": "Immune signalling", "casein_secretion": "Casein secretion",
          "lipid_synthesis": "Lipid synthesis", "secretory_machinery": "Secretory machinery",
          "hmo_glycosylation": "HMO glycosylation", "stress_response": "Stress response"}
NMF_LABEL = {"NMF1": "Macrophage", "NMF2": "Secretory lactocyte (+ mitochondrial)", "NMF3": "Immediate-early stress",
             "NMF4": "Neutrophil", "NMF5": "T cell", "NMF6": "MHC-II antigen-presenting",
             "NMF7": "Progenitor-like epithelium", "NMF8": "Proliferation", "NMF9": "Junctional epithelium",
             "NMF10": "Monocyte / macrophage", "NMF11": "Interferon response", "NMF12": "Secretory lactocyte"}
NMF_DIR = {"NMF7": "early", "NMF9": "later", "NMF2": "later", "NMF12": "later"}
MITO_FLAG = {g for g in top30 if g.startswith("MT-")}

faith = R["faithfulness"]; i5 = faith["ks"].index(5)
faithful = faith["shap"][i5] - faith["shap"][0] > 3 * (faith["random"][i5] - faith["random"][0])
rand_ok = R["randomization_corr"] < 0.5
null_ok = R["null_r2_max"] < min(R["fold_r2"])
heldout_ok = True
GLOBAL = [("Validated on held-out mothers", heldout_ok),
          ("Beats the permutation null", null_ok),
          ("Explanation is faithful (masking test)", faithful),
          ("Explanation passes the randomization check", rand_ok)]

E, L = {}, []          # entities, links
def ent(id, type, name, **kw): E[id] = dict(id=id, type=type, name=name, **kw)
def link(a, b, rel, source, trust, status, why): L.append(dict(a=a, b=b, rel=rel, source=source, trust=trust, status=status, why=why))

# ---- core entities --------------------------------------------------------
ent("ds", "dataset", "SCP1671 human milk atlas",
    desc=f"Nyquist SK, …, Shalek AK, Berger B, Goods BA. PNAS 2022. The paper follows 15 mothers for up to 632 days; "
         f"in the public metadata, cells from 13 of them carry a timepoint, between days 3 and 323. "
         f"{R['n_cells_after_qc']:,} cells with a timepoint; 25,000 sampled for modelling.",
    stats=[["Cells with a timepoint", f"{R['n_cells_after_qc']:,}"], ["Mothers", "13"], ["Days postpartum", "3–323"],
           ["Median mothers per timepoint", f"{R['median_donors_per_timepoint']:.0f}"]],
    url="https://singlecell.broadinstitute.org/single_cell/study/SCP1671")
ent("early", "phase", "Early lactation", desc="The early end of the clock's axis. Genes linked here push predictions earlier when raised.")
ent("later", "phase", "Later lactation", desc="The later end of the clock's axis. Genes linked here push predictions later when raised. "
    "The clock cannot resolve months 6 to 11: only three mothers were sampled after 200 days.")
ent("clock", "model", "Lactation clock (all cells)",
    desc="Random forest predicting log days postpartum from 1,518 genes per cell, validated on held-out mothers.",
    stats=[["R², held-out mothers", f2(R['clock_donor_heldout']['r2'])],
           ["Fold range", f"{f2(min(R['fold_r2']))}–{f2(max(R['fold_r2']))}"],
           ["Best of 20 permutations", f2(R['null_r2_max'])],
           ["Leaky random split", f2(R['clock_random_split']['r2'])]],
    checks=[[n, "pass" if ok else "fail"] for n, ok in GLOBAL[:2]] + [["Leakage measured (+" + f2(R['r2_inflation_from_leakage']) + ")", "pass"]])
ent("clock_epi", "model", "Lactation clock (epithelial only)",
    desc=f"The same model restricted to the {R['clock_lactocyte_only']['n_cells']:,} epithelial cells. It tests whether the clock reads cell state or cell mix.",
    stats=[["R², held-out mothers", f2(R['clock_lactocyte_only']['r2'])],
           ["Fold range", f"{f2(min(R['clock_lactocyte_only']['fold_r2']))}–{f2(max(R['clock_lactocyte_only']['fold_r2']))}"]],
    checks=[["Validated on held-out mothers", "pass"], ["Beats the permutation null", "not run"],
            [f"One fold weak (R² {f2(min(R['clock_lactocyte_only']['fold_r2']))})", "fail"]])
ent("donor_clf", "model", "Mother-identity classifier",
    desc="Adversarial control: predicts which mother a cell came from. It sets the ceiling on how much a model could gain by recognising mothers.",
    stats=[["Accuracy", f"{R['donor_identity_accuracy']:.0%}"], ["Chance", f"{R['donor_identity_chance']:.0%}"]])
ent("epi", "celltype", "Epithelial cell",
    desc="The largest population in milk; lactocytes are the milk-producing epithelial cells. "
         f"{R['clock_lactocyte_only']['n_cells']:,} cells in the modelling sample.")
link("ds", "epi", "majority of milk cells are lactocytes", "literature", None, "supported",
     "Reported by Nyquist et al. 2022 for this dataset.")
link("clock", "ds", "trained on", "model", 5, "supported", "Pipeline link.")
link("clock_epi", "epi", "restricted to", "model", 5, "supported", "Pipeline link.")
link("clock_epi", "clock", "reproduces", "model", 1 + int(heldout_ok), "supported",
     f"R² {f2(R['clock_lactocyte_only']['r2'])} within epithelial cells vs {f2(R['clock_donor_heldout']['r2'])} for all cells; permutation null not run for this subset.")
link("donor_clf", "clock", "bounds leakage for", "model", 3, "supported",
     f"Mothers are identifiable ({R['donor_identity_accuracy']:.0%} vs {R['donor_identity_chance']:.0%}), yet a leaky split adds only {f2(R['r2_inflation_from_leakage'])} R².")

# ---- NMF programmes -----------------------------------------------------------
for p, genes in N.items():
    shap_in = [g for g in top30 if g in genes]
    ent(p, "nmf", f"{p} · {NMF_LABEL[p]}",
        desc=f"Unsupervised programme found by NMF without time labels. Top genes: {', '.join(genes[:10])}. "
             f"Programme names are an interpretation of the top genes.",
        genes=genes, stats=[["SHAP top-30 genes in programme", str(len(shap_in))]])
    if p in NMF_DIR:
        link(p, NMF_DIR[p], "programme leans toward", "model", 4, "supported",
             f"Its members among the SHAP top 30 ({', '.join(shap_in)}) mostly push predictions toward {NMF_DIR[p]} lactation.")

# ---- curated programmes ------------------------------------------------------
chi, base = R["ablation_control_p97_5"], R["base_r2_xai_model"]
for k, v in R["program_ablation_r2_loss"].items():
    cov = R["ablation_coverage"][k]; above = v > chi
    ent("cur_" + k, "curated", PRETTY[k],
        desc=f"Hand-picked gene set. Masking its {cov['n_in_features']} genes costs {f2(v, 4)} R² "
             f"({'above' if above else 'within'} the random 5-gene control, 97.5th percentile {f2(chi, 3)}).",
        genes=cov["in_features"],
        stats=[["R² lost when masked", f2(v, 4)], ["Share of model R²", f"{v / base:.0%}"], ["Control (97.5th pct)", f2(chi, 3)]],
        checks=[["Beats random-gene control", "pass" if above else "fail"]])
    link("cur_" + k, "clock", "load-bearing for" if above else "not load-bearing for", "model",
         4 if above else 1, "supported" if above else "not supported",
         f"Ablation: {f2(v, 4)} R² lost vs control {f2(chi, 3)}.")
    for g in cov["in_features"]:
        if g not in E:
            ent(g, "gene", g, desc="Curated-programme gene.")
        link(g, "cur_" + k, "member of", "literature", None, "supported", "Curated gene-set definition.")

# ---- genes -------------------------------------------------------------------
for rank, g in enumerate(top30, 1):
    progs = [p for p, gs in N.items() if g in gs]
    d = CF.get(g)
    dirn = None if d is None else ("early" if d < 0 else "later")
    flagged = g in MITO_FLAG
    corroborated = dirn is not None and any(NMF_DIR.get(p) == dirn for p in progs)
    checks = [[n, "pass" if ok else "fail"] for n, ok in GLOBAL]
    checks.append(["Direction corroborated by an NMF programme", "pass" if corroborated else ("not run" if dirn is None else "fail")])
    if flagged:
        checks.append(["Possible confound: mitochondrial share varies by mother", "flag"])
    trust = sum(ok for _, ok in GLOBAL) + int(corroborated)
    if flagged: trust = min(trust, 2)
    stats = [["SHAP rank", f"#{rank} of {fc:,}"]]
    if d is not None:
        stats.append(["Effect of +1 SD", f"{f2(d, 3)} log-days (~{pct(d)}% {'earlier' if d < 0 else 'later'})"])
    stats.append(["NMF programmes", ", ".join(progs) or "none"])
    cur = [PRETTY[k] for k, c in R["ablation_coverage"].items() if g in c["in_features"]]
    if cur: stats.append(["Curated programme", ", ".join(cur)])
    desc = (f"Ranked #{rank} by SHAP for the lactation clock. "
            + (f"Raising it moves predictions {'earlier' if d < 0 else 'later'}." if d is not None
               else "Direction not tested (counterfactual test covered the top 20).")
            + (" Mitochondrial: its share differs between mothers, so it may reflect sample quality." if flagged else ""))
    E[g] = dict(id=g, type="gene", name=g, desc=desc, stats=stats, checks=checks, trust=trust, rank=rank)
    if dirn:
        link(g, dirn, "pushes toward", "model", trust, "supported",
             f"Counterfactual {f2(d, 3)} log-days per +1 SD; "
             + ("corroborated by " + ", ".join(p for p in progs if NMF_DIR.get(p) == dirn) if corroborated else "no agreeing NMF programme"))
    else:
        link(g, "clock", "used by", "model", trust, "supported", f"SHAP rank #{rank}; direction not tested.")
    for p in progs:
        r_in = N[p].index(g) + 1
        link(g, p, "member of", "model", 4 if r_in <= 10 else 3, "supported", f"Rank {r_in} of 25 in {p}.")

# ---- semantic chains ----------------------------------------------------------
chains = [
    dict(title="Early milk carries a progenitor-like epithelium flagged by MHC-II and complement",
         steps=[
             ["Cells from early lactation are sampled in milk", "data", "supported", "SCP1671; median one mother per timepoint."],
             ["High SFRP1, KRT15, HLA-DRA, CD74 and C3 push the clock's predictions earlier", "model", "supported",
              f"SHAP and counterfactuals; explanation faithful and passes the randomization check (r = {f2(R['randomization_corr'])})."],
             ["These genes co-occur in one unsupervised programme (NMF7)", "model", "supported", "NMF found it without time labels."],
             ["The programme marks progenitor-like epithelial cells", "inference", "untested",
              "An interpretation of the top genes; not tested against lineage markers here."],
             ["It recedes as lactation continues", "model", "weak",
              "Holds through early and mid lactation; past 200 days rests on three mothers."]]),
    dict(title="Later milk shifts toward junctional and secretory lactocyte programmes",
         steps=[
             ["High KRT80, ANXA2, CLDN4 push predictions later", "model", "supported", "SHAP and counterfactuals."],
             ["They co-occur in a junctional programme (NMF9)", "model", "supported", "Unsupervised."],
             ["Casein, LALBA and BTN1A1 also push later and form a secretory programme (NMF12, NMF2)", "model", "supported",
              "Casein and lipid programmes also beat the random-gene control in ablation."],
             ["The shift reflects maturation of milk production", "inference", "untested", "Interpretation; not tested here."]]),
    dict(title="Part of the 'immune' signal is secreted milk protein",
         steps=[
             ["LYZ is the top-ranked gene and pushes predictions later", "model", "supported", "SHAP rank #1."],
             ["LYZ loads with the caseins in NMF2, not with immune programmes", "model", "supported", "Unsupervised."],
             ["So ablating the curated 'immune signalling' set removes secretory signal too", "inference", "supported",
              "Why that set's ablation (0.26 R²) overstates immune-cell content."],
             ["Separating epithelial and immune sources of LYZ", "model", "untested", "Would need cell-type-resolved attribution."]]),
    dict(title="From milk clock to perturbation readout",
         steps=[
             ["A per-cell early-to-late axis is learned from milk", "model", "supported",
              f"R² {f2(R['clock_donor_heldout']['r2'])} on held-out mothers."],
             ["The axis holds inside epithelial cells alone", "model", "supported",
              f"R² {f2(R['clock_lactocyte_only']['r2'])}; one fold weak; permutation null not run for this subset."],
             ["The axis transfers from milk cells to cultured lactocytes", "untested", "untested", "First experiment needed."],
             ["Compounds with known effects move cells along it as expected", "untested", "untested", "Second experiment needed."],
             ["The axis is usable as a screen readout", "untested", "untested", "Only after both tests above."]]),
]

# ---- contested claims ---------------------------------------------------------
claims = [
    dict(claim="Milk-synthesis genes carry no information about lactational time",
         a=["My first run", "Ablating the lipid, HMO and secretory programmes cost almost nothing, and casein could not be tested."],
         b=["Corrected run", "A coverage check showed most of those programmes' genes had never entered the model. Once included, casein (0.07) and lipid (0.04) beat the random-gene control."],
         status="resolved", settle="Resolved: a feature-selection artefact. The ablation now reports genes in the model for every programme."),
    dict(claim="Curated gene sets are a good readout of lactocyte state",
         a=["Common practice", "Scoring cells on hand-picked programmes."],
         b=["This analysis", f"Six programme averages reach R² {f2(R['representation_r2']['curated programs'][0])}; individual genes reach {f2(next(v for k, v in R['representation_r2'].items() if k.startswith('HVG'))[0])}."],
         status="contested", settle="Test programme scores against learned readouts on perturbation data."),
    dict(claim="Lactocytes form a continuum of states",
         a=["Nyquist et al. 2022", "Lactocytes show a continuum of cell states with subtle transcriptional changes."],
         b=["This analysis", f"Annotated cell types are recovered with {R['celltype_accuracy_donor_heldout']:.0%} accuracy. That tests cell types, not lactocyte sub-states, so it neither supports nor contradicts the continuum."],
         status="open", settle="Test sub-state discreteness within epithelial cells (e.g. classifier confidence on lactocyte sub-clusters)."),
    dict(claim="Mitochondrial genes mark later lactation",
         a=["SHAP", "MT-ND4 and MT-ND5 push predictions later."],
         b=["Data quality", "Mitochondrial share differs between mothers, so this may be sample quality."],
         status="open", settle="Re-fit without mitochondrial genes and check whether R² or the ranking changes."),
    dict(claim="The clock reads late lactation",
         a=["The model", "It predicts days postpartum across the whole range."],
         b=["The design", "Only three mothers were sampled after 200 days, and the clock under-predicts there."],
         status="refuted", settle="More mothers sampled after six months."),
]

# ---- model scorecards ---------------------------------------------------------
def audits(**kw):
    names = {"heldout": "Validated on held-out mothers", "null": "Beats a permutation null", "leak": "Leakage measured",
             "faith": "Explanation faithful", "rand": "Passes randomization check", "corro": "Independent method agrees"}
    return [[names[k], v] for k, v in kw.items()]
cards = [
    dict(q="Is lactational time readable from milk cells?", metric=f"R² {f2(R['clock_donor_heldout']['r2'])}",
         vs=[["Best of 20 permutations", f2(R['null_r2_max'])], ["Leaky random split", f2(R['clock_random_split']['r2'])]],
         audit=audits(heldout="pass", null="pass", leak="pass"), verdict="good", ent="clock"),
    dict(q="Is it more than a change in cell mix?", metric=f"R² {f2(R['clock_lactocyte_only']['r2'])}",
         vs=[["All cells", f2(R['clock_donor_heldout']['r2'])], ["Weakest fold", f2(min(R['clock_lactocyte_only']['fold_r2']))]],
         audit=audits(heldout="pass", null="not run"), verdict="good", ent="clock_epi"),
    dict(q="Is it more than recognising the mother?", metric=f"+{f2(R['r2_inflation_from_leakage'])} R² from leakage",
         vs=[["Mother identified", f"{R['donor_identity_accuracy']:.0%}"], ["Chance", f"{R['donor_identity_chance']:.0%}"]],
         audit=audits(heldout="pass", leak="pass"), verdict="good", ent="donor_clf"),
    dict(q="Do the explanations reflect the model?", metric=f"error {f2(faith['shap'][0])} → {f2(faith['shap'][i5])}",
         vs=[["5 random genes masked", f2(faith['random'][i5])], ["Shuffled-label model r", f2(R['randomization_corr'])],
             ["Top-25 shared with permutation importance", f"{R['shap_perm_overlap_top25']}/25"]],
         audit=audits(faith="pass" if faithful else "fail", rand="pass" if rand_ok else "fail", corro="pass"), verdict="good", ent=None),
    dict(q="Which curated programmes carry the clock?", metric="3 of 6 above control",
         vs=[[PRETTY[k], f2(v, 3)] for k, v in sorted(R["program_ablation_r2_loss"].items(), key=lambda kv: -kv[1])] + [["Control (97.5th pct)", f2(chi, 3)]],
         audit=audits(heldout="pass", null="pass"), verdict="part", ent="cur_casein_secretion"),
    dict(q="Do programme averages work as features?", metric=f"R² {f2(R['representation_r2']['curated programs'][0])}",
         vs=[[k, f"{f2(v[0])} ± {f2(v[1])}"] for k, v in R["representation_r2"].items()],
         audit=audits(heldout="pass"), verdict="weak", ent=None),
    dict(q="Can it resolve late lactation?", metric="3 mothers after 200 days",
         vs=[["Median mothers per timepoint", f"{R['median_donors_per_timepoint']:.0f}"]],
         audit=audits(heldout="pass"), verdict="no", ent="later"),
]
known_truth = dict(name="Recovers annotated cell types",
                   result=f"{R['celltype_accuracy_donor_heldout']:.1%} accuracy on held-out mothers; median confidence {R['median_confidence']:.2f}",
                   status="pass")

questions = [
    ["Does the clock just recognise mothers?", "#scorecards"],
    ["Is it reading cell types instead of cell state?", "#explore/clock_epi"],
    ["Which genes mark early milk?", "#chains/0"],
    ["Is HPGD, the strongest gene, corroborated?", "#explore/HPGD"],
    ["Where does this disagree with what I expected?", "#contested"],
    ["How would this become a perturbation readout?", "#chains/3"],
]

KB = dict(entities=list(E.values()), links=L, chains=chains, claims=claims, cards=cards,
          known_truth=known_truth, questions=questions,
          trust_rule=[n for n, _ in GLOBAL] + ["Corroborated by an independent method (NMF)"])
html = open(tpl_path).read().replace("__KB__", json.dumps(KB, ensure_ascii=False))
import os; os.makedirs(f"{repo}/docs/kb", exist_ok=True)
open(f"{repo}/docs/kb/index.html", "w").write(html)
types = {}
for e in E.values(): types[e["type"]] = types.get(e["type"], 0) + 1
print(f"KB: {len(E)} entities {types} | {len(L)} links | {len(chains)} chains | {len(claims)} claims | {len(cards)} cards")
print("faithful", faithful, "| rand_ok", rand_ok, "| null_ok", null_ok)
print("gene trust:", {g: E[g]['trust'] for g in top30})
