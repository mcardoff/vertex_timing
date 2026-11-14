#!/usr/bin/env python3
"""Where the analysis veto's room is (results/veto_room.md; slides/slides_timing_followups.tex), on tdr_fwdpair trees
(--pair=one --jvt=jvt4pt60):
  left   Z+jets pileup VBS legs in HGTD with a core jet time: t_core - t_HS, against a Gaussian of 248 ps (two
         independent collisions, each with a 175 ps time spread) of the same area. The legs the veto misses are the ones
         inside its ~90 ps window -- collisions that happened in time with the hard scatter.
  right  VBF hard-scatter legs in HGTD: the veto pull against the truth t0, by core track count. One- and two-track
         cores carry the tails that veto signal.

  cd python && PYTHONNOUSERSITE=1 ~/.venv-hgtd/bin/python gp_veto_room.py \\
      ../condor/vbf/vbf_onefwd_jvt4pt60_wfix_tdr_fwdpair.root ../condor/zjets/zjets_onefwd_jvt4pt60_wfix_tdr_fwdpair.root \\
      --out ../figs/tdr/room/veto_room
"""
import argparse

import numpy as np
import uproot

import tdr_fwdpair_common as FC
import tdr_root_style as S

ap = argparse.ArgumentParser()
ap.add_argument("sig")
ap.add_argument("bkg")
ap.add_argument("--out", required=True)
args = ap.parse_args()
R = S.setup()
Es = uproot.open(args.sig)["events"].arrays(library="np")
Eb = uproot.open(args.bkg)["events"].arrays(library="np")
f = FC.calibration(Es, verbose=False)
fwd = lambda E, k: np.abs(E[f"j{k}_eta"]) > 2.4
DT_PU = 247.8   # ps: t_PU - t_HS for independent collisions (memory pileup-time-model)

cv = S.keep(R.TCanvas("veto_room", "", 1600, 750))
TS = 0.042
# ── left: pileup legs ──
cv.cd()
S.pad("p1", 0.0, 0.0, 0.5, 1.0, top=0.06, bottom=0.15, right=0.04)
dt, w = [], []
for k in (1, 2):
    m = (Eb[f"j{k}_pu"] == 1) & fwd(Eb, k) & (Eb[f"j{k}_core_n"] >= 1)
    dt.append((Eb[f"j{k}_core_t"] - Eb["t0_truth"])[m]); w.append(Eb["weight"][m].astype(float))
dt, w = np.concatenate(dt), np.concatenate(w)
h = S.keep(R.TH1D("h_pu", "", 40, -800, 800))
for x, wi in zip(dt, w): h.Fill(min(max(x, -799.0), 799.0), wi)
h.Scale(1.0 / h.Integral())
S.style(h, "kP10Red", 0, 1, 3)
g = S.keep(R.TF1("g", f"[0]*exp(-0.5*(x/{DT_PU})^2)", -800, 800))
g.SetParameter(0, h.GetBinWidth(1) / (np.sqrt(2 * np.pi) * DT_PU))
g.SetLineColor(R.kGray + 2); g.SetLineStyle(2); g.SetLineWidth(3); g.SetNpx(400)
ymax = 1.7 * max(h.GetMaximum(), g.GetMaximum())
fr = S.keep(R.gPad.DrawFrame(-800, 0, 800, ymax))
S.axes(fr, "t_{core} #minus t_{HS} [ps]", "Fraction of legs", size=0.05, ytoff=1.5, xtoff=1.3)
box = S.keep(R.TBox(-100, 0, 100, ymax / 1.7)); box.SetFillColorAlpha(R.kGray, 0.35); box.Draw()
h.Draw("HIST SAME"); g.Draw("SAME")
inw = 100 * np.average(np.abs(dt) < 100, weights=w)
from math import erf, sqrt
exp_in = 100 * erf(100 / (sqrt(2) * DT_PU))
leg = S.legend(0.55, 0.55, 0.88, 0.70, size=TS, margin=0.3)
leg.AddEntry(h, "Pileup legs", "L"); leg.AddEntry(g, "Gaussian, 248 ps", "L"); leg.Draw()
S.header(0.20, 0.87, ["#sqrt{s} = 14 TeV, #LT#mu#GT = 200, Z+jets", "Timed pileup VBS legs, 2.4 < |#eta| < 3.8",
                      f"Within #pm100 ps (grey): {inw:.0f}% (independent: {exp_in:.0f}%)"], size=TS)
S.end_pad()
# ── right: hard-scatter legs' pull ──
cv.cd()
S.pad("p2", 0.5, 0.0, 1.0, 1.0, top=0.06, bottom=0.15, right=0.04)
R.gPad.SetLogy()
P = FC.pulls(Es, f)
pl, nn = [], []
for k in (1, 2):
    m = (Es[f"j{k}_hs"] == 1) & fwd(Es, k) & (Es[f"j{k}_core_n"] >= 1)
    pl.append(np.nan_to_num(P[f"truth_{k}"][m], nan=0)); nn.append(Es[f"j{k}_core_n"][m])
pl, nn = np.concatenate(pl), np.concatenate(nn)
GR = [("1 core track", nn == 1, "kP10Red", 1), ("2 core tracks", nn == 2, "kP10Orange", 2), ("#geq 3 core tracks", nn >= 3, "kP10Blue", 7)]
hs_, ymax = [], 0
for name, m, col, ls in GR:
    hh = S.keep(R.TH1D(f"h_{ls}", "", 40, 0, 10))
    for x in pl[m]: hh.Fill(min(x, 9.99))
    hh.Scale(1.0 / max(m.sum(), 1)); S.style(hh, col, 0, ls, 3)
    hs_.append((hh, name, 100 * np.mean(pl[m] >= 3), 100 * np.mean(m)))
    ymax = max(ymax, hh.GetMaximum())
fr2 = S.keep(R.gPad.DrawFrame(0, 1e-4, 10, 100 * ymax))
S.axes(fr2, "|t_{core} #minus t_{0}^{truth}| / #sigma", "Fraction of legs", size=0.05, ytoff=1.5, xtoff=1.3)
cut = S.keep(R.TLine(3, 1e-4, 3, ymax)); cut.SetLineStyle(2); cut.SetLineColor(R.kGray + 2); cut.Draw()
leg2 = S.legend(0.56, 0.50, 0.93, 0.67, size=TS, margin=0.3)
for hh, name, out, share in hs_:
    hh.Draw("HIST SAME"); leg2.AddEntry(hh, name, "L")
leg2.Draw()
S.header(0.20, 0.87, ["#sqrt{s} = 14 TeV, #LT#mu#GT = 200, VBF H#rightarrowinv.",
                      "Hard-scatter VBS legs in HGTD, truth t_{0}",
                      "Pull #geq 3: " + " / ".join(f"{o:.1f}%" for _, _, o, _ in hs_) + " (legs: " + " / ".join(f"{sh:.0f}%" for _, _, _, sh in hs_) + ")"], size=TS)
S.end_pad()
S.save(cv, args.out)
print(f"PU legs within 100 ps: {inw:.1f}% (independent collisions {exp_in:.1f}%)")
for hh, name, out, share in hs_: print(f"HS legs, {name}: {share:.0f}% of legs, {out:.1f}% with pull >= 3")
