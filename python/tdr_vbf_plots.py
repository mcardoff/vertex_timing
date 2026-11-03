#!/usr/bin/env python3
"""TDR Sec. 3.4.1 figures from util/tdr_vbf's event tree.

The TDR preselection: >= 2 jets, leading pT > 75 GeV, second > 50 GeV, |Delta eta(j1, j2)| > 3, after
"all jets are required to pass an ITk-only R_pT pileup jet tagger operating at 85% hard-scatter
efficiency". The tagger is applied here as a jet cleaning: a jet failing it is dropped, and the
selection runs on the jets that remain (jets beyond |eta| = 4.0 have no tracks and are kept). A timing
tagger replaces the ITk one at its own 85% working point, from tdr_pujet_plots.py's json.

  <out>_cat_mjj   TDR Fig 3.24: the fraction of preselected events with m_jj above a threshold (dashed)
                  and its split into central-central / central-forward / forward-forward (forward:
                  |eta| > 2.4), ITk-only tagger; signal solid, background (--bkg) dotted
  <out>_pupair    per category and tagger: the share of preselected events whose VBF pair (the two
                  leading surviving jets) contains a pileup jet (TDR labels) -- hatched where that
                  jet is at |eta| >= 4.0, beyond tracking and HGTD, where no R_pT tagger can act --
                  and the preselected events relative to the ITk-only tagger

Events are weighted by `weight`; errors are binomial on the effective entries.

    cd python && PYTHONNOUSERSITE=1 ~/.venv-hgtd/bin/python tdr_vbf_plots.py --sig ../figs/hists/tdr_vbf.root \\
        --wp ../figs/tdr/local_tdr_pujet_wp85.json --out ../figs/tdr/local_tdr_vbf
"""
import argparse, json, os
import numpy as np
import uproot
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches

ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
ap.add_argument("--sig", required=True, help="tdr_vbf output, VBF H->inv")
ap.add_argument("--bkg", help="tdr_vbf output, Z+jets (optional)")
ap.add_argument("--wp", required=True, help="tdr_pujet_plots.py's <out>_wp85.json")
ap.add_argument("--sig-label", default="VBF $H\\to$inv.")
ap.add_argument("--bkg-label", default="Z+jets")
ap.add_argument("--note", default="", help="appended to the footer, e.g. '(local sample)'")
ap.add_argument("--out", required=True)
args = ap.parse_args()
plt.rcParams.update({"font.size": 12})
os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)

J1, J2, DETA = 75.0, 50.0, 3.0
# (key, legend, colour) -- the colours of tdr_pujet_plots.py; None = no pileup tagging at all
TAGGERS = [(None, "no pileup-jet tagging", "#B8B8B8"),
           ("rpt_itk", "ITk only (TDR)", "#000000"),
           ("rpt_t0_hgtd", "ITk + HGTD, $t_0$: HGTD (Athena)", "#3B8ED0"),
           ("rpt_t0_trkptz", "ITk + HGTD, $t_0$: TRKPTZ", "#C0392B"),
           ("rpt_t0_waves", "ITk + HGTD, $t_0$: WAVeS", "#E69F00"),
           ("rpt_t0_truth", "ITk + HGTD, truth $t_0$ (reference)", "#A0522D")]
CATS = ["CC", "CF", "FF"]
WP = json.load(open(args.wp))
gate_wp = WP["gate_sigma"]


def load(path):
    t = uproot.open(path)
    E = t["events"].arrays(library="np")
    g = float(np.unique(t["meta"].arrays(["gate_sigma"], library="np")["gate_sigma"])[0])
    if abs(g - gate_wp) > 1e-9:
        raise SystemExit(f"{path}: gate {g} sigma, but the working points were made at {gate_wp}")
    return E


def thresholds(col):
    w = WP["working_points"][col]
    return {(r, s): w[f"{r}_{s}"]["threshold"] for r in ("central", "forward") for s in ("30_50", "50_inf")}


def select(E, col):
    """per event: passes, m_jj, category, pileup jet in the pair, one beyond |eta| = 4.0 --
    after dropping the jets failing `col`"""
    n = len(E["weight"])
    ok = np.zeros(n, bool); mjj = np.zeros(n); cat = np.full(n, -1); pupair = np.zeros(n, bool); pu40 = np.zeros(n, bool)
    th = thresholds(col) if col else None
    for i in range(n):
        pt, eta, phi = E["jet_pt"][i], E["jet_eta"][i], E["jet_phi"][i]
        keep = np.ones(len(pt), bool)
        if col:
            r = E["jet_" + col][i]
            ae = np.abs(eta)
            thr = np.where(ae < 2.4, np.where(pt < 50, th["central", "30_50"], th["central", "50_inf"]),
                           np.where(pt < 50, th["forward", "30_50"], th["forward", "50_inf"]))
            keep = (r < 0) | (r >= thr)                       # r < 0: beyond |eta| = 4.0, no tracks
        s = np.flatnonzero(keep)                              # the tree's jets are pT-ordered
        if len(s) < 2: continue
        a, b = s[0], s[1]
        if pt[a] <= J1 or pt[b] <= J2 or abs(eta[a] - eta[b]) <= DETA: continue
        ok[i] = True
        v = [pt[a] * np.cosh(eta[a]), pt[b] * np.cosh(eta[b])]
        px = pt[a] * np.cos(phi[a]) + pt[b] * np.cos(phi[b]); py = pt[a] * np.sin(phi[a]) + pt[b] * np.sin(phi[b])
        pz = pt[a] * np.sinh(eta[a]) + pt[b] * np.sinh(eta[b])
        mjj[i] = np.sqrt(max((v[0] + v[1]) ** 2 - px * px - py * py - pz * pz, 0.0))
        cat[i] = int(abs(eta[a]) > 2.4) + int(abs(eta[b]) > 2.4)
        pupair[i] = bool(E["jet_pu"][i][a] or E["jet_pu"][i][b])
        pu40[i] = any(E["jet_pu"][i][j] and abs(eta[j]) >= 4.0 for j in (a, b))
    return ok, mjj, cat, pupair, pu40


samples = [("sig", args.sig_label, load(args.sig), "-")]
if args.bkg: samples.append(("bkg", args.bkg_label, load(args.bkg), ":"))
R = {}
for key, lab, E, _ in samples:
    for col, *_ in TAGGERS:
        R[key, col] = select(E, col)
    # the all-jets selection must reproduce tdr_vbf's own flag
    assert np.array_equal(R[key, None][0], E["presel_all"]), "offline preselection disagrees with presel_all"

foot = f"$t_0$ gate {gate_wp:g}$\\sigma$, R$_{{p_T}}$ taggers at 85% hard-scatter efficiency   ATLAS Simulation Internal"
foot = (args.note + "; " if args.note else "") + foot

# ── Fig 3.24 analogue ─────────────────────────────────────────────────────────
CC = {"CC": "#222222", "CF": "#C0392B", "FF": "#3B8ED0"}       # the TDR's black / red / blue
X = np.linspace(0, 3000, 301)
fig, ax = plt.subplots(figsize=(8.6, 5.6))
print("Fig 3.24 analogue (ITk-only tagger): fraction of preselected events with m_jj above the threshold")
for key, lab, E, ls in samples:
    ok, mjj, cat, *_ = R[key, "rpt_itk"]
    w = E["weight"].astype(float)
    tot = w[ok].sum()
    frac = lambda m: np.array([w[m & (mjj > x)].sum() for x in X]) / tot
    ax.plot(X, frac(ok), color="#777777", lw=2.0, ls=(0, (6, 3)) if key == "sig" else (0, (1, 2)), label=f"{lab}, all")
    for c, name in enumerate(CATS):
        ax.plot(X, frac(ok & (cat == c)), color=CC[name], lw=2.0, ls=ls, label=f"{lab}, {name}")
    print(f"  {lab}: {ok.sum()} events;  " + "  ".join(f"{name} {100 * w[ok & (cat == c)].sum() / tot:.1f}%" for c, name in enumerate(CATS))
          + "  | m_jj > 1 TeV: " + "  ".join(f"{name} {100 * w[ok & (cat == c) & (mjj > 1000)].sum() / tot:.1f}%" for c, name in enumerate(CATS)))
ax.set_xlim(0, 3000); ax.set_ylim(0, 1.05)
ax.set_xlabel("$m_{jj}$ threshold [GeV]"); ax.set_ylabel("fraction of preselected events")
ax.grid(alpha=0.25); ax.legend(fontsize=9.5, frameon=False, loc="upper right")
ax.text(0.98, 0.62, "VBF preselection, ITk-only R$_{p_T}$ tagger at 85%\nforward: $|\\eta|$ > 2.4",
        transform=ax.transAxes, fontsize=9.5, color="#555555", ha="right", va="top")
fig.text(0.5, -0.01, foot, ha="center", fontsize=9, color="#555555")
for ext in ("png", "pdf"): fig.savefig(f"{args.out}_cat_mjj.{ext}", dpi=140, bbox_inches="tight")
print("wrote", args.out + "_cat_mjj.png")

# ── pileup jets in the VBF pair, per category and tagger ──────────────────────
for key, lab, E, _ in samples:
    w = E["weight"].astype(float)
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(14, 5.4), gridspec_kw=dict(wspace=0.22))
    nT = len(TAGGERS); bw = 0.8 / nT
    ref = R[key, "rpt_itk"]
    print(f"\n{lab}: per category -- events (weighted share with a pileup jet in the pair; of which beyond |eta| 4.0) [relative to ITk only]")
    print("  " + f"{'tagger':40s}" + "".join(f"{c:>30s}" for c in CATS))
    for k, (col, name, colour) in enumerate(TAGGERS):
        ok, _, cat, pp, p40 = R[key, col]
        frac, f40, err, rel, rerr, line = [], [], [], [], [], ""
        for c in range(3):
            m = ok & (cat == c)
            sw, sw2 = w[m].sum(), (w[m] ** 2).sum()
            f = w[m & pp].sum() / sw if sw > 0 else np.nan
            neff = sw * sw / sw2 if sw2 > 0 else 0
            frac.append(100 * f); err.append(100 * np.sqrt(f * (1 - f) / neff) if neff > 0 else np.nan)
            f40.append(100 * w[m & p40].sum() / sw if sw > 0 else np.nan)
            sref = w[ref[0] & (ref[2] == c)].sum()
            rel.append(sw / sref if sref > 0 else np.nan)
            rerr.append(np.sqrt(sw2) / sref if sref > 0 else np.nan)
            line += f"{int(m.sum()):>7d} ({100 * f:5.2f}%; {f40[-1]:5.2f}%) [{rel[-1]:5.3f}]"
        x = np.arange(3) + (k - (nT - 1) / 2) * bw
        a1.bar(x, frac, bw * 0.9, yerr=err, color=colour, edgecolor="white", lw=1, label=name, capsize=2, error_kw=dict(lw=1, ecolor="#444444"))
        a1.bar(x, f40, bw * 0.9, color="none", edgecolor="white", lw=0, hatch="////")
        if col is None:
            noTag = rel
            continue
        a2.bar(x, rel, bw * 0.9, yerr=rerr, color=colour, edgecolor="white", lw=1, capsize=2, error_kw=dict(lw=1, ecolor="#444444"))
        print("  " + f"{name.replace('$', ''):40s}" + line)
    for a in (a1, a2):
        a.set_xticks(range(3)); a.set_xticklabels(["central-central", "central-forward", "forward-forward"])
        a.grid(alpha=0.25, axis="y")
    a1.set_ylabel("events with a pileup jet in the VBF pair [%]")
    a2.set_ylabel("preselected events, relative to ITk only")
    a2.axhline(1.0, color="#888888", lw=0.8)
    a2.set_ylim(0.95, 1.05)
    a2.text(0.02, 0.97, "without any tagging: " + " / ".join(f"{v:.2f}" for v in noTag) + " (CC / CF / FF)",
            transform=a2.transAxes, fontsize=9.5, color="#555555", va="top")
    a1.set_ylim(0, 1.12 * np.nanmax([b.get_height() for b in a1.patches]))
    h, l = a1.get_legend_handles_labels()
    h.append(matplotlib.patches.Patch(facecolor="#B8B8B8", edgecolor="white", hatch="////")); l.append("hatched: that pileup jet is at $|\\eta|$ $\\geq$ 4.0 (no tracks)")
    fig.legend(h, l, fontsize=9.5, frameon=False, loc="lower center", ncol=4, bbox_to_anchor=(0.5, -0.1))
    a1.set_title(f"{lab}: pileup jets in the selected pair", fontsize=12)
    a2.set_title(f"{lab}: preselected events", fontsize=12)
    fig.text(0.5, -0.13, foot, ha="center", fontsize=9, color="#555555")
    for ext in ("png", "pdf"): fig.savefig(f"{args.out}_pupair_{key}.{ext}", dpi=140, bbox_inches="tight")
    print("wrote", f"{args.out}_pupair_{key}.png")
