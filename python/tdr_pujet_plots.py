#!/usr/bin/env python3
"""TDR Sec. 3.3.3 figures from util/tdr_pujet's per-jet tree.

  <out>_roc      TDR Fig 3.20: pileup-jet rejection (1 / pileup efficiency) against hard-scatter
                 efficiency (0.80-1.0), 2.4 < |eta| < 4.0, for 30-50 GeV and > 50 GeV jets; every
                 tagger against ITk only, with the ratio to ITk underneath
  <out>_eta      TDR Fig 3.21: pileup-jet rate relative to ITk only, against |eta|, at 85% hard-scatter
                 efficiency (the working point set per |eta| bin), 30-50 GeV
  <out>_wp85.json  every tagger's R_pT threshold at 85% hard-scatter efficiency per pT bin (30-50,
                 > 50 GeV) and region (central |eta| < 2.4, forward 2.4-4.0): the TDR preselection's
                 "every jet passes the ITk-only R_pT tagger at 85%", and the timing taggers that
                 replace it, read by tdr_vbf_plots.py

Jets are counted unweighted, as the TDR's per-jet performance. A jet passes a tagger when its R_pT is
>= the threshold; the curve is traced over every threshold, and rejection is read at matched
hard-scatter efficiency, never by snapping to bins. Between two thresholds the pileup efficiency is
interpolated linearly: that is the tagger accepting the jets between them at random. Above the last
threshold before R_pT = 0 (every tagger has a spike of hard-scatter jets at R_pT = 0, larger with a
wrong t0) that interpolation is the only way to reach the efficiency, so the curve there is a
mixture, not a cut. The working points in the json are real cuts: the highest threshold keeping >= 85%.

    cd python && PYTHONNOUSERSITE=1 ~/.venv-hgtd/bin/python tdr_pujet_plots.py \\
        ../figs/hists/tdr_pujet.root --label "VBF H#to inv. (local)" --out ../figs/tdr/local_tdr_pujet
"""
import argparse, json, os
import numpy as np
import uproot
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
ap.add_argument("tree", help="tdr_pujet output")
ap.add_argument("--label", default="VBF $H\\to$inv.")
ap.add_argument("--out", required=True)
args = ap.parse_args()
plt.rcParams.update({"font.size": 12})
os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)

# taggers: (column, legend, colour, marker, line style); ITk first (the reference)
TAGGERS = [("rpt_itk", "ITk only", "#000000", "", (0, (6, 3))),
           ("rpt_self", "ITk + HGTD, self-tagging", "#2E9E4F", "^", "-"),
           ("rpt_t0_hgtd", "ITk + HGTD, $t_0$: HGTD (Athena)", "#3B8ED0", "v", "-"),
           ("rpt_comb", "ITk + HGTD, combined (Athena $t_0$, else self)", "#7B3294", "D", "-"),
           ("rpt_t0_trkptz", "ITk + HGTD, $t_0$: TRKPTZ", "#C0392B", "s", "-"),
           ("rpt_t0_waves", "ITk + HGTD, $t_0$: WAVeS", "#E69F00", "o", "-"),
           ("rpt_t0_truth", "ITk + HGTD, truth $t_0$ (reference)", "#A0522D", "", (0, (2, 2)))]
EFF = np.round(np.arange(0.800, 0.9751, 0.025), 3)    # the sampled working points

t = uproot.open(args.tree)
J = t["jets"].arrays(library="np")
meta = t["meta"].arrays(["gate_sigma", "n_read"], library="np")
gate = float(np.unique(meta["gate_sigma"])[0])
ae = np.abs(J["eta"])


def curve(hs, pu):
    """the ROC: for each distinct hard-scatter efficiency (ascending), the lowest pileup efficiency
    reaching it and the threshold (the highest one) that does; pass = R_pT >= threshold"""
    th = np.unique(np.concatenate([hs, pu]))
    e_hs = 1.0 - np.searchsorted(np.sort(hs), th, side="left") / len(hs)
    e_pu = 1.0 - np.searchsorted(np.sort(pu), th, side="left") / len(pu)
    u, inv = np.unique(e_hs, return_inverse=True)
    pmin = np.full(len(u), np.inf); np.minimum.at(pmin, inv, e_pu)
    tmax = np.full(len(u), -np.inf); np.maximum.at(tmax, inv, th)
    return u, pmin, tmax


def at_eff(hs, pu, effs):
    """pileup efficiency and its binomial error at each hard-scatter efficiency (interpolated)"""
    e_hs, e_pu, _ = curve(hs, pu)
    pe = np.interp(effs, e_hs, e_pu)
    return pe, np.sqrt(np.clip(pe * (1 - pe), 0, None) / len(pu))


def working_point(hs, pu, eff):
    """a real cut: the highest threshold keeping >= eff of the hard-scatter jets"""
    e_hs, e_pu, th = curve(hs, pu)
    i = np.flatnonzero(e_hs >= eff)[0]
    return float(th[i]), float(e_hs[i]), float(e_pu[i])


def sel(pt_lo, pt_hi, eta_lo, eta_hi):
    return (J["pt"] > pt_lo) & (J["pt"] < pt_hi) & (ae >= eta_lo) & (ae < eta_hi)


SLICES = [("30 < $p_T$ < 50 GeV", 30.0, 50.0), ("$p_T$ > 50 GeV", 50.0, 1e9)]

# ── Fig 3.20 analogue ─────────────────────────────────────────────────────────
fig, axes = plt.subplots(2, 2, figsize=(14, 7.4), sharex=True, gridspec_kw=dict(height_ratios=[2.2, 1], hspace=0.05, wspace=0.18))
print(f"{args.tree}: gate {gate:g} sigma\n\nrejection at matched hard-scatter efficiency, 2.4 < |eta| < 4.0 (ratio to ITk in brackets)")
for c, (lab, lo, hi) in enumerate(SLICES):
    m = sel(lo, hi, 2.4, 4.0)
    nhs, npu = int(np.sum(m & (J["hs"] == 1))), int(np.sum(m & (J["pu"] == 1)))
    ref = None
    print(f"\n  {lab.replace('$', '')}: {nhs} hard-scatter / {npu} pileup jets")
    print("  " + f"{'tagger':48s}" + "".join(f"{e:>14.3f}" for e in EFF))
    for col, name, colour, mk, ls in TAGGERS:
        hs, pu = J[col][m & (J["hs"] == 1)], J[col][m & (J["pu"] == 1)]
        dense = np.linspace(0.80, 1.0, 401)
        pd_ = at_eff(hs, pu, dense)[0]
        axes[0, c].plot(dense, np.where(pd_ > 0, 1.0 / pd_, np.nan), color=colour, lw=1.8, ls=ls, label=name)
        pe, perr = at_eff(hs, pu, EFF)
        rej = np.where(pe > 0, 1.0 / pe, np.nan)
        rerr = np.where(pe > 0, perr / pe ** 2, np.nan)
        if mk: axes[0, c].errorbar(EFF, rej, yerr=rerr, fmt=mk, color=colour, ms=6, capsize=2, lw=1.2)
        if ref is None: ref = rej
        ratio = rej / ref
        axes[1, c].plot(EFF, ratio, color=colour, lw=1.6, ls=ls, marker=mk or None, ms=5)
        print("  " + f"{name.replace('$', ''):48s}" + "".join(f"{r:7.1f} ({q:4.2f})" for r, q in zip(rej, ratio)))
    axes[0, c].set_yscale("log"); axes[0, c].set_xlim(0.80, 1.0)
    axes[0, c].set_title(f"{lab}, 2.4 < $|\\eta|$ < 4.0", fontsize=12)
    axes[0, c].grid(alpha=0.25, which="both"); axes[1, c].grid(alpha=0.25)
    axes[1, c].axhline(1.0, color="#888888", lw=0.8)
    axes[1, c].set_xlabel("hard-scatter efficiency")
    axes[0, c].text(0.97, 0.95, f"{nhs:,} hard-scatter / {npu:,} pileup jets", transform=axes[0, c].transAxes, fontsize=9.5, color="#555555", ha="right", va="top")
axes[0, 0].set_ylabel("pileup-jet rejection"); axes[1, 0].set_ylabel("ratio to ITk")
axes[0, 1].legend(fontsize=9, frameon=False, loc="lower left")
fig.text(0.5, 0.005, f"{args.label}, $\\langle\\mu\\rangle$ = 200; $R_{{p_T}}$ tagger, $t_0$ gate {gate:g}$\\sigma$; markers at the sampled working points   "
         "ATLAS Simulation Internal", ha="center", fontsize=9.5, color="#555555")
for ext in ("png", "pdf"): fig.savefig(f"{args.out}_roc.{ext}", dpi=140, bbox_inches="tight")
print("\nwrote", args.out + "_roc.png")

# ── Fig 3.21 analogue ─────────────────────────────────────────────────────────
EB = np.array([2.4, 2.8, 3.2, 3.6, 4.0])
fig, ax = plt.subplots(figsize=(8.4, 5.2))
rel, nbin = {}, []
for col, name, colour, mk, ls in TAGGERS:
    v = []
    for b in range(len(EB) - 1):
        m = sel(30.0, 50.0, EB[b], EB[b + 1])
        hs, pu = J[col][m & (J["hs"] == 1)], J[col][m & (J["pu"] == 1)]
        hi, pi = J["rpt_itk"][m & (J["hs"] == 1)], J["rpt_itk"][m & (J["pu"] == 1)]
        pe = at_eff(hs, pu, np.array([0.85]))[0][0]
        pr = at_eff(hi, pi, np.array([0.85]))[0][0]
        v.append(pe / pr if pr > 0 else np.nan)
        if col == "rpt_itk": nbin.append((len(hs), len(pu)))
    rel[col] = v
    ax.stairs(v, EB, color=colour, lw=1.8, ls=ls, label=name, baseline=None)
ax.set_xlim(2.4, 4.0); ax.set_ylim(0, max(1.6, 1.1 * np.nanmax([x for v in rel.values() for x in v])))
ax.set_xlabel("jet $|\\eta|$"); ax.set_ylabel("relative pileup-jet rate (to ITk)")
ax.legend(fontsize=8.6, frameon=False, loc="upper left", ncol=1)
ax.grid(alpha=0.25)
ax.text(0.02, 0.04, "85% hard-scatter efficiency (per $|\\eta|$ bin), 30 < $p_T$ < 50 GeV", transform=ax.transAxes, fontsize=9.5, color="#555555")
fig.text(0.5, -0.01, f"{args.label}, $t_0$ gate {gate:g}$\\sigma$   ATLAS Simulation Internal", ha="center", fontsize=9.5, color="#555555")
for ext in ("png", "pdf"): fig.savefig(f"{args.out}_eta.{ext}", dpi=140, bbox_inches="tight")
print("\nrelative pileup-jet rate at 85%, 30-50 GeV, |eta| bins " + " ".join(f"{EB[b]:.1f}-{EB[b + 1]:.1f}" for b in range(len(EB) - 1)))
for col, name, *_ in TAGGERS:
    print(f"  {name.replace('$', ''):48s}" + "".join(f"{x:8.2f}" for x in rel[col]))
print(f"  {'(hard-scatter / pileup jets per bin)':48s}" + "".join(f"  {h}/{p}" for h, p in nbin))
print("wrote", args.out + "_eta.png")

# ── the ITk-only 85% working points, for the TDR preselection ────────────────
wp = {}
print("\nR_pT working points at 85% hard-scatter efficiency (threshold / achieved HS eff / pileup eff):")
for col, name, *_ in TAGGERS:
    wp[col] = {}
    line = f"  {name.replace('$', ''):48s}"
    for region, (elo, ehi) in {"central": (0.0, 2.4), "forward": (2.4, 4.0)}.items():
        for lab, lo, hi in SLICES:
            m = sel(lo, hi, elo, ehi)
            hs, pu = J[col][m & (J["hs"] == 1)], J[col][m & (J["pu"] == 1)]
            thr, ehs, epu = working_point(hs, pu, 0.85)
            key = f"{region}_{int(lo)}_{'inf' if hi > 1e8 else int(hi)}"
            wp[col][key] = dict(threshold=thr, hs_eff=ehs, pu_eff=epu, n_hs=int(len(hs)), n_pu=int(len(pu)))
            line += f"  {key} {thr:.3f}/{ehs:.3f}/{epu:.4f}"
    print(line)
json.dump(dict(tree=os.path.abspath(args.tree), hs_eff=0.85, gate_sigma=gate, working_points=wp), open(args.out + "_wp85.json", "w"), indent=1)
print("wrote", args.out + "_wp85.json")
