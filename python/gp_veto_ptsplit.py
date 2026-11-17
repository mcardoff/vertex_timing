#!/usr/bin/env python3
"""Timing only where JVT acts (pT < 60 GeV) against the veto on every leg (results/veto_ptsplit.md): the S/sqrt(B)
change of the jet-vs-t0 veto (WAVeS, truth t0) and of the perfect HGTD veto when only legs in a pT range may veto,
before the BDT and at its 80% cut. Same inputs as tdr_fwdpair_veto_ptsplit.py.

  cd python && PYTHONNOUSERSITE=1 ~/.venv-hgtd/bin/python gp_veto_ptsplit.py <vbf tree> <zjets tree> <kinbdt npz> \\
      --out ../figs/tdr/room/veto_ptsplit
"""
import argparse

import numpy as np
import uproot

import tdr_fwdpair_common as FC
import tdr_root_style as S

ap = argparse.ArgumentParser()
ap.add_argument("sig"); ap.add_argument("bkg"); ap.add_argument("bdt")
ap.add_argument("--out", required=True)
args = ap.parse_args()
R = S.setup()
Es = uproot.open(args.sig)["events"].arrays(library="np")
Eb = uproot.open(args.bkg)["events"].arrays(library="np")
D = np.load(args.bdt)
f = FC.calibration(Es, verbose=False)
wS, wB = Es["weight"].astype(float), Eb["weight"].astype(float)
o = np.argsort(-D["sig_score"]); cs = np.cumsum(wS[o]) / wS.sum()
CUTS = [("Before the BDT", None), ("BDT at #varepsilon_{S} = 80%", D["sig_score"][o][np.searchsorted(cs, 0.8)])]
RANGES = [("every leg", 0, 1e9), ("< 60 GeV", 0, 60), ("< 50 GeV", 0, 50), ("#geq 60 GeV", 60, 1e9)]
SRCS = [("waves", "WAVeS t_{0}", S.COL["waves"], 20), ("truth", "Truth t_{0}", S.COL["truth"], 25),
        ("perfect", "Perfect HGTD veto", "#000000", 34)]
fwd = lambda E, k: np.abs(E[f"j{k}_eta"]) > 2.4


def vetoed(E, src, lo, hi):
    v = np.zeros(len(E["weight"]), bool)
    P = None if src == "perfect" else FC.pulls(E, f)
    for k in (1, 2):
        inr = (E[f"j{k}_pt"] >= lo) & (E[f"j{k}_pt"] < hi)
        v |= inr & (((E[f"j{k}_pu"] == 1) & fwd(E, k)) if src == "perfect" else (np.nan_to_num(P[f"{src}_{k}"], nan=0.0) >= 3.0))
    return v


cv = S.keep(R.TCanvas("ptsplit", "", 1600, 750))
for c, (clab, thr) in enumerate(CUTS):
    bS = np.ones(len(wS), bool) if thr is None else D["sig_score"] >= thr
    bB = np.ones(len(wB), bool) if thr is None else D["bkg_score"] >= thr
    cv.cd()
    S.pad(f"p{c}", 0.5 * c, 0.0, 0.5 * c + 0.5, 1.0, top=0.06, bottom=0.15, right=0.04)
    vals = {}
    for s, *_ in SRCS:
        vals[s] = []
        for _, lo, hi in RANGES:
            vS, vB = vetoed(Es, s, lo, hi), vetoed(Eb, s, lo, hi)
            vals[s].append(100 * ((wS[bS & ~vS].sum() / wS[bS].sum()) / np.sqrt(wB[bB & ~vB].sum() / wB[bB].sum()) - 1))
    ymax = 1.75 * max(max(v) for v in vals.values())
    fr = S.keep(R.gPad.DrawFrame(0, 0, len(RANGES), ymax))
    S.axes(fr, "Legs allowed to veto", "S/#sqrt{B} change [%]", size=0.05, ytoff=1.3, xtoff=1.2)
    fr.GetXaxis().SetNdivisions(0)
    for i, (lab, *_) in enumerate(RANGES):
        S.text(0.16 + 0.80 * (i + 0.5) / len(RANGES), 0.10, lab, 0.042, align=21)
    leg = S.legend(0.40, 0.47, 0.78, 0.65, size=0.042, margin=0.3)
    for j, (s, lab, col, mk) in enumerate(SRCS):
        g = S.style(S.graph(np.arange(len(RANGES)) + 0.35 + 0.15 * j, vals[s]), col, mk, 1, 2, 1.8)
        g.Draw("P SAME"); leg.AddEntry(g, lab, "P")
    leg.Draw()
    pts = np.concatenate([Eb[f"j{k}_pt"][(Eb[f"j{k}_pu"] == 1) & fwd(Eb, k) & bB] for k in (1, 2)])
    wl = np.concatenate([wB[(Eb[f"j{k}_pu"] == 1) & fwd(Eb, k) & bB] for k in (1, 2)])
    S.header(0.20, 0.87, ["#sqrt{s} = 14 TeV, #LT#mu#GT = 200, VBF H#rightarrowinv. vs Z+jets", f"JVT to |#eta| < 4.0 (p_{{T}} < 60), {clab}",
                          f"Pileup legs in HGTD #geq 60 GeV: {100 * wl[pts >= 60].sum() / wl.sum():.0f}%"], size=0.042)
    S.end_pad()
    print(clab, {s: [f"{v:+.0f}%" for v in vals[s]] for s in vals})
S.save(cv, args.out)
