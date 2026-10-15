#!/usr/bin/env python3
"""How many tracks must a jet's time be built from? The timing-veto summary of
vbs_time_veto_plot.py at --min-trk k = 1..N (a jet with fewer tracks in its time counts
as untimed, so neither test uses it), per case, drawn as
  top     S/sqrt(B) change vs no veto, per method
  bottom  background efficiency (filled) and signal efficiency (open), per method

Input: <dir>/<case>_k<k>_summary.md, written by
    vbs_time_veto_plot.py ... --acceptance {both,one,all} --min-trk k --out <dir>/<case>_k<k>

    PYTHONNOUSERSITE=1 ~/.venv-hgtd/bin/python python/vbs_mintrk_scan.py \\
        --dir figs/time_veto/mintrk_scan --out figs/time_veto/mintrk_scan/mintrk_scan
"""
import argparse, os, re
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
ap.add_argument("--dir", required=True)
ap.add_argument("--kmax", type=int, default=5)
ap.add_argument("--out", required=True)
args = ap.parse_args()
plt.rcParams.update({"font.size": 12})

CASES = [("both", "Forward-Forward"), ("one", "Central-Forward"), ("all", "All selected events")]
METHODS = [("jj", "Jet time compatibility", "black", "o"), ("t0_trkptz", "TRKPTZ $t_0$", "#C0392B", "s"),
           ("t0_waves", "WAVeS $t_0$", "#E69F00", "^"), ("t0_hgtd", "HGTD $t_0$ (Athena)", "#3B8ED0", "v"),
           ("t0_truth", "truth $t_0$", "#8C5A3C", "P")]
ROW = re.compile(r"^\|\s*(\S+)\s*\|\s*([\d.]+)\s*\+-\s*[\d.]+%\s*\|\s*([\d.]+)\s*\+-\s*[\d.]+%\s*\|.*\|\s*([\d.]+)\s*\|\s*$")


def read(path):
    """{method: (eps_S, eps_B, S/sqrt(B) ratio)} from the FIRST method table of a summary"""
    out, started = {}, False
    for ln in open(path):
        if ln.startswith("| method |"):
            if started: break
            started = True; continue
        if started:
            m = ROW.match(ln.strip())
            if m: out[m.group(1)] = (float(m.group(2)), float(m.group(3)), float(m.group(4)))
            elif not ln.startswith("|"): break
    return out


ks = np.arange(1, args.kmax + 1)
fig, axes = plt.subplots(2, 3, figsize=(15, 7.6), sharex=True, gridspec_kw=dict(hspace=0.06, wspace=0.22))
for c, (case, title) in enumerate(CASES):
    data = [read(os.path.join(args.dir, f"{case}_k{k}_summary.md")) for k in ks]
    a0, a1 = axes[0, c], axes[1, c]
    for key, lab, col, mk in METHODS:
        if key not in data[0] or (case == "one" and key == "jj"):
            continue
        z = [100 * (d[key][2] - 1) for d in data]
        a0.plot(ks, z, marker=mk, color=col, lw=2, ms=7, label=lab)
        a1.plot(ks, [d[key][1] for d in data], marker=mk, color=col, lw=2, ms=7)
        a1.plot(ks, [d[key][0] for d in data], marker=mk, color=col, lw=1.2, ms=7, ls="--", mfc="white")
    a0.axhline(0, color="#999999", lw=0.8, ls=":")
    a0.set_title(title, fontsize=13, fontweight="bold", loc="left")
    a1.set_xlabel("tracks a jet time must be built from ($\\geq k$)")
    a1.set_xticks(ks); a1.set_xticklabels([f"$\\geq${k}" for k in ks])
    a0.grid(alpha=0.3); a1.grid(alpha=0.3)
    a1.set_ylim(40, 101)
    if c == 0:
        a0.set_ylabel("$S/\\sqrt{B}$ change vs no veto [%]")
        a1.set_ylabel("efficiency [%]\nfilled: $Z$+jets, open: VBF")
axes[0, 0].legend(fontsize=10, frameon=False, loc="upper right")
fig.suptitle("Timing veto vs the jet-time track requirement: VBF $H\\to$inv. against $Z(\\to\\ell\\ell)$+jets, "
             "JVT + fJVT loose, recoil > 100 GeV   ATLAS Simulation Internal", fontsize=11.5, x=0.08, ha="left", y=0.97)
os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)
for ext in ("pdf", "png"):
    fig.savefig(f"{args.out}.{ext}", dpi=150, bbox_inches="tight")
print("wrote", args.out + ".pdf / .png")
