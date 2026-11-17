"""Redraw F16 from the run's results.json numbers (the original's auto-title was wrong)."""
import json, sys
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

R = json.load(open(sys.argv[1]))
out = sys.argv[2]
abl, cov = R["program_ablation_r2_loss"], R["ablation_coverage"]
ctrl, base = R["ablation_control_p97_5"], R["base_r2_xai_model"]
pretty = {"immune_signalling": "Immune signalling", "casein_secretion": "Casein secretion",
          "lipid_synthesis": "Lipid synthesis", "secretory_machinery": "Secretory machinery",
          "hmo_glycosylation": "HMO glycosylation", "stress_response": "Stress response"}

items = sorted(abl.items(), key=lambda kv: kv[1])
names = [f"{pretty[k]}  ({cov[k]['n_in_features']} genes)" for k, _ in items]
vals = [v for _, v in items]
above = [v > ctrl for v in vals]

plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10.5,
                     "axes.spines.top": False, "axes.spines.right": False})
fig, ax = plt.subplots(figsize=(8.6, 4.4), dpi=200)
ax.axvspan(0, ctrl, color="#e5e7eb", zorder=0)
bars = ax.barh(names, vals, height=0.6, zorder=3,
               color=["#2563eb" if a else "#9ca3af" for a in above])
ax.axvline(0, color="#111", lw=0.9, zorder=4)
for b, v, a in zip(bars, vals, above):
    txt = f"{v:.3f}  ({v/base:.0%})" if a else f"{v:.4f}  within control"
    ax.text(max(v, 0) + 0.004, b.get_y() + b.get_height() / 2, txt, va="center",
            fontsize=9.5, color="#111" if a else "#6b7280", zorder=5)
ax.text(ctrl + 0.003, -0.75, f"random 5-gene sets\n(97.5th percentile {ctrl:.3f})",
        fontsize=8.5, color="#6b7280", va="bottom")
ax.set_xlim(-0.005, max(vals) * 1.32)
ax.set_ylim(-0.9, len(vals) - 0.45)
ax.set_xlabel(f"R² lost when the programme's genes are masked (model R² {base:.2f})")
if "notitle" not in sys.argv:
    ax.set_title("Immune, casein and lipid programmes beat a random-gene control;\n"
                 "HMO glycosylation, secretory machinery and stress response do not",
                 loc="left", fontsize=11.5)
fig.tight_layout()
fig.savefig(out, dpi=200)
print("wrote", out)
