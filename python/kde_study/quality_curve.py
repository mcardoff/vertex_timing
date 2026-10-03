"""Figure: how often the provided time is right vs how often a time is provided, as the Q threshold moves.
One panel per sample; Athena's own working point marked. Output: figs/quality/quality_curve.png (+ csv table)."""
import os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from qual import *

BLUE, ORANGE, INK, MUTED = "#2a78d6", "#eb6834", "#1f1f1e", "#6b6a63"
TITLES = {'zjets_novbs': 'Z+jets', 'dijet_novbs': 'dijet', 'vbf_novbs': 'VBF H→inv', 'ttbar_novbs': 'ttbar'}
out = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "figs", "quality")
os.makedirs(out, exist_ok=True)
fig, axs = plt.subplots(1, 4, figsize=(15, 4.2), sharey=False)
rows = ["sample,threshold,time_provided_pct,correct_of_provided_pct"]
for ax, (n, m) in zip(axs, SAMPLES):
    d = get(n, m); q = qvars(d)
    Q = q['S1mS2'] / np.sqrt(2 * q['S1'] - q['S1mS2'])
    o = np.argsort(-Q); ok = q['ok'][o]
    acc = np.arange(1, len(ok) + 1) / len(ok); pur = np.cumsum(ok) / np.arange(1, len(ok) + 1)
    s = acc > 0.05
    ax.plot(100 * acc[s], 100 * pur[s], color=BLUE, lw=2)
    for thr in (1.0, 2.0, 3.0):
        k = Q >= thr; a, p = 100 * k.mean(), 100 * q['ok'][k].mean()
        ax.plot(a, p, 'o', ms=8, color=BLUE, mec='white', mew=2)
        ax.annotate(f"Q ≥ {thr:g}", (a, p), textcoords="offset points", xytext=(-6, 8), ha='right', fontsize=9, color=INK)
        rows.append(f"{TITLES[n]},Q>={thr:g},{a:.1f},{p:.1f}")
    v, aok = athena(d); a, p = 100 * v.mean(), 100 * aok[v].mean()
    ax.plot(a, p, 's', ms=9, color=ORANGE, mec='white', mew=2)
    ax.annotate("Athena", (a, p), textcoords="offset points", xytext=(-8, -14), ha="right", fontsize=9, color=INK)
    rows.append(f"{TITLES[n]},Athena,{a:.1f},{p:.1f}")
    rows.append(f"{TITLES[n]},no cut,100.0,{100*q['ok'].mean():.1f}")
    ax.set_title(TITLES[n], fontsize=11, color=INK, loc='left')
    ax.set_xlabel("events with a time provided [%]", color=MUTED, fontsize=9)
    ax.set_xlim(0, 102); ax.grid(True, color="#e6e5df", lw=0.6); ax.set_axisbelow(True)
    for sp in ('top', 'right'): ax.spines[sp].set_visible(False)
    for sp in ('left', 'bottom'): ax.spines[sp].set_color("#c9c8c0")
    ax.tick_params(colors=MUTED, labelsize=8)
axs[0].set_ylabel("provided times within 60 ps of truth [%]", color=MUTED, fontsize=9)
fig.suptitle("Time-quality flag: purity of the provided time against how often one is provided (blue: new method as the Q threshold moves; orange: Athena)",
             fontsize=11, color=INK, x=0.01, ha='left')
fig.tight_layout(rect=(0, 0, 1, 0.94))
fig.savefig(os.path.join(out, "quality_curve.png"), dpi=130, facecolor="#fcfcfb")
open(os.path.join(out, "quality_curve.csv"), "w").write("\n".join(rows) + "\n")
print("\n".join(rows))
