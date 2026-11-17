#!/usr/bin/env python3
"""The R_pT guard on the timing veto (results/jet_core_variants.md): the S/sqrt(B) change against the threshold below
which a leg's ITk R_pT must be for its time to veto, with the standard core jet time and without the 30 GeV track cap
(nocap), WAVeS and truth t0, before the BDT and at its 80% cut. Inputs: tdr_fwdpair trees with the core-time variants
(commit 014c810 on) and the kinematic BDT.

  cd python && PYTHONNOUSERSITE=1 ~/.venv-hgtd/bin/python gp_veto_rptguard.py <vbf tree> <zjets tree> <kinbdt npz> \\
      --out ../figs/tdr/room/veto_rptguard
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
S0 = uproot.open(args.sig)["events"].arrays(library="np")
B0 = uproot.open(args.bkg)["events"].arrays(library="np")
D = np.load(args.bdt)
GUARDS = [0.05, 0.1, 0.15, 0.2, 0.3, 0.5, np.inf]
LAB = ["0.05", "0.1", "0.15", "0.2", "0.3", "0.5", "none"]
ROWS = [("", "waves", "WAVeS t_{0}", S.COL["waves"], 20, 1), ("nocap", "waves", "WAVeS t_{0}, no 30 GeV cap", S.COL["waves"], 24, 2),
        ("", "truth", "Truth t_{0}", S.COL["truth"], 21, 1), ("nocap", "truth", "Truth t_{0}, no 30 GeV cap", S.COL["truth"], 25, 2)]
wS, wB = S0["weight"].astype(float), B0["weight"].astype(float)
o = np.argsort(-D["sig_score"]); cs = np.cumsum(wS[o]) / wS.sum()
CUTS = [("Before the BDT", None), ("BDT at #varepsilon_{S} = 80%", D["sig_score"][o][np.searchsorted(cs, 0.8)])]


def view(E, v):
    if not v: return E
    X = dict(E)
    for k in (1, 2):
        for q in ("t", "sig", "n"): X[f"j{k}_core_{q}"] = E[f"j{k}_core_{q}_{v}"]
    return X


cache = {}
for v in ("", "nocap"):
    Sv, Bv = view(S0, v), view(B0, v)
    f = FC.calibration(Sv, verbose=False)
    cache[v] = (Sv, Bv, FC.pulls(Sv, f), FC.pulls(Bv, f))
veto = lambda E, P, src, X: np.any([(np.nan_to_num(P[f"{src}_{k}"], nan=0) >= 3) & (E[f"j{k}_rpt_itk"] < X) for k in (1, 2)], axis=0)

cv = S.keep(R.TCanvas("rptguard", "", 1600, 750))
for c, (clab, thr) in enumerate(CUTS):
    bS = np.ones(len(wS), bool) if thr is None else D["sig_score"] >= thr
    bB = np.ones(len(wB), bool) if thr is None else D["bkg_score"] >= thr
    vals = []
    for v, src, *_ in ROWS:
        Sv, Bv, PS, PB = cache[v]
        vals.append([100 * ((wS[bS & ~veto(Sv, PS, src, X)].sum() / wS[bS].sum())
                            / np.sqrt(wB[bB & ~veto(Bv, PB, src, X)].sum() / wB[bB].sum()) - 1) for X in GUARDS])
    cv.cd()
    S.pad(f"p{c}", 0.5 * c, 0.0, 0.5 * c + 0.5, 1.0, top=0.06, bottom=0.20, right=0.04)
    lo, hi = min(min(x) for x in vals), max(max(x) for x in vals)
    fr = S.keep(R.gPad.DrawFrame(-0.5, lo - 2, len(GUARDS) - 0.5, hi + 1.0 * (hi - lo) + 5))
    S.axes(fr, "A leg may veto only if its ITk R_{p_{T}} <", "S/#sqrt{B} change [%]", size=0.05, ytoff=1.3, xtoff=1.5)
    fr.GetXaxis().SetNdivisions(0)
    for i, l in enumerate(LAB):
        S.text(0.16 + 0.80 * (i + 0.5) / len(GUARDS), 0.155, l, 0.042, align=21)
    S.text(0.96, 0.06, "A leg may veto only if its ITk R_{p_{T}} <", 0.045, align=31)   # the 0-division axis drops its title
    leg = S.legend(0.20, 0.55, 0.62, 0.72, size=0.036, margin=0.3)
    for (v, src, lab, col, mk, ls), y in zip(ROWS, vals):
        g = S.style(S.graph(np.arange(len(GUARDS)), y), col, mk, ls, 2, 1.5); g.Draw("LP SAME"); leg.AddEntry(g, lab, "LP")
    leg.Draw()
    S.header(0.20, 0.87, ["#sqrt{s} = 14 TeV, #LT#mu#GT = 200, VBF H#rightarrowinv. vs Z+jets", f"JVT to |#eta| < 4.0 (p_{{T}} < 60), {clab}"], size=0.042)
    S.end_pad()
    print(clab, {r[2]: [f"{x:+.1f}" for x in y] for r, y in zip(ROWS, vals)})
S.save(cv, args.out)
