#!/usr/bin/env python3
"""Section 3 with each JVT variant as the baseline (2026-11-09 mini deck, slides/slides_jvt_variants.tex): S/sqrt(B)
relative to the analysis baseline (JVT to |eta| < 4.0, pT < 60 GeV, no fJVT: jvt4pt60) without the veto, against the
kinematic BDT's signal efficiency, for each variant alone and with the 3-sigma WAVeS t0 veto on top.

Each variant has its own tdr_fwdpair --pair=one --jvt=<wp> trees (JVT removes jets BEFORE the pair is formed, so the
pairs differ) and its own kinematic BDT (tdr_fwdpair_kinbdt.py), cut at the signal efficiency on the x axis relative to
that variant's events. S and B are absolute weighted yields, so a variant that removes signal at the JVT step pays for
it here (the method of tdr_fwdpair_baselines.py). Errors: weighted Poisson on B only.

  --strat LABEL:WP:COLOUR ...   (trees ../condor/<s>/<s>_onefwd_<WP>_tdr_fwdpair.root, BDT
                                 ../condor/zjets/vbf_vs_zjets_onefwd_<WP>_kinbdt.npz; the first is the reference)
"""
import argparse

import numpy as np
import uproot

import tdr_fwdpair_common as FC
import tdr_root_style as S

ap = argparse.ArgumentParser()
ap.add_argument("--strat", nargs="+", required=True, help="LABEL:WP:COLOUR, the reference (jvt4pt60) first")
ap.add_argument("--title", default="")
ap.add_argument("--out", required=True)
args = ap.parse_args()
R = S.setup()
EFFS = np.round(np.arange(1.0, 0.249, -0.05), 2)


def load(wp):
    D = np.load(f"../condor/zjets/vbf_vs_zjets_onefwd_{wp}_kinbdt.npz")
    Es, Eb = (uproot.open(f"../condor/{s}/{s}_onefwd_{wp}_tdr_fwdpair.root")["events"].arrays(library="np") for s in ("vbf", "zjets"))
    cal = FC.calibration(Es, verbose=False)
    pu = lambda E: ((E["j1_pu"] == 1) & (np.abs(E["j1_eta"]) > 2.4)) | ((E["j2_pu"] == 1) & (np.abs(E["j2_eta"]) > 2.4))
    return dict(ss=D["sig_score"], sb=D["bkg_score"], ws=Es["weight"].astype(float), wb=Eb["weight"].astype(float),
                Ks=FC.vetoes(Es, cal)["waves"], Kb=FC.vetoes(Eb, cal)["waves"], Ps=~pu(Es), Pb=~pu(Eb))


def yields(T, veto, e):
    ss, ws = T["ss"], T["ws"]
    o = np.argsort(-ss)
    thr = -np.inf if e >= 1.0 else ss[o][np.searchsorted(np.cumsum(ws[o]) / ws.sum(), e)]
    ms, mb = ss >= thr, T["sb"] >= thr
    if veto == "waves": ms, mb = ms & T["Ks"], mb & T["Kb"]
    if veto == "perfect": ms, mb = ms & T["Ps"], mb & T["Pb"]
    return T["ws"][ms].sum(), T["wb"][mb].sum(), np.sum(T["wb"][mb] ** 2), int(mb.sum())


strats = [s.split(":") for s in args.strat]
T = {wp: load(wp) for _, wp, _ in strats}
ref = [yields(T[strats[0][1]], None, e) for e in EFFS]
rows = []
print("S/sqrt(B) relative to the reference without the veto, at BDT signal efficiency " + " ".join(f"{e:.2f}" for e in EFFS))
for lab, wp, col in strats:
    for veto in (None, "waves", "perfect"):
        y, ey, nb = [], [], []
        for (s0, b0, _, _), e in zip(ref, EFFS):
            s, b, vb, n = yields(T[wp], veto, e)
            f = (s / s0) / np.sqrt(b / b0)
            y.append(f); ey.append(f * 0.5 * np.sqrt(vb) / b); nb.append(n)
        rows.append((lab, wp, col, veto, np.array(y), np.array(ey), nb))
        print(f"  {lab.replace('#', ''):34s} {'+ ' + veto if veto else '':10s} " + " ".join(f"{a:.3f}" for a in y)
              + f"   (Z+jets MC {nb[0]} / {nb[4]} / {nb[10]})")
s0, b0 = ref[0][0], ref[0][1]
for lab, wp, col in strats:
    s, b, _, _ = yields(T[wp], None, 1.0)
    print(f"  {lab.replace('#', ''):34s} before the BDT: S {100 * s / s0:.1f}%, B {100 * b / b0:.1f}% of the reference's")

MAIN, RAT = 0.68, 0.32
TS, TSR = 0.05 * 600 / (800 * MAIN), 0.05 * 600 / (800 * RAT)
cv = S.keep(R.TCanvas("jvtvar", "", 800, 800))
cv.cd(); S.pad("jv_m", 0.0, RAT, 1.0, 1.0, top=0.05 / MAIN * 0.6, bottom=0.02)
show = [r for r in rows if r[3] != "perfect"]
top = max(np.nanmax(r[4] + r[5]) for r in show)
bot = min(0.8, 0.05 * np.floor(20 * (min(np.nanmin(r[4] - r[5]) for r in show) - 0.03)))
fr = S.keep(R.gPad.DrawFrame(0.225, bot, 1.025, bot + 2.2 * (top - bot)))   # headroom: header + legend
S.axes(fr, "", "S/#sqrt{B} relative to the baseline", size=TS, ytoff=1.35)
fr.GetXaxis().SetLabelSize(0)
one = S.keep(R.TLine(0.225, 1.0, 1.025, 1.0)); one.SetLineStyle(2); one.SetLineColor(R.kGray + 2); one.Draw()
leg = S.legend(0.19, 0.50, 0.92, 0.66, size=0.55 * TS, ncol=2, margin=0.25)
for lab, wp, col, veto, y, ey, _ in show:
    g = S.style(S.graph(EFFS, y, ey), col, 20 if veto is None else 24, 1 if veto is None else 2, 2, 1.1); g.Draw("LP SAME")
    leg.AddEntry(g, lab + (" + WAVeS veto" if veto else ""), "LP")
leg.Draw()
S.header(0.19, 0.88, ["#sqrt{s} = 14 TeV, HL-LHC, #LT#mu#GT = 200", "S: VBF H#rightarrowinv., B: Z+jets, #geq 1 VBS jet in HGTD",
                      args.title or "each JVT variant's own pairs and kinematic BDT"], size=0.8 * TS)
R.gPad.RedrawAxis()
cv.cd(); S.pad("jv_r", 0.0, 0.0, 1.0, RAT, top=0.03, bottom=0.36)
fr2 = S.keep(R.gPad.DrawFrame(0.225, 0.95, 1.025, 1.4))
S.axes(fr2, "Signal efficiency of the kinematic BDT cut (1 = no BDT)", "Veto: with / without", size=TSR, ytoff=1.35 * TS / TSR, xtoff=1.2, ndiv_y=503)
one2 = S.keep(R.TLine(0.225, 1.0, 1.025, 1.0)); one2.SetLineStyle(2); one2.SetLineColor(R.kGray + 2); one2.Draw()
by = {(r[1], r[3]): r for r in rows}
for lab, wp, col in strats:
    g2 = S.style(S.graph(EFFS, by[(wp, "waves")][4] / by[(wp, None)][4]), col, 24, 2, 2, 1.0); g2.Draw("LP SAME")
S.save(cv, args.out)
