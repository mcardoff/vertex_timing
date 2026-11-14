#!/usr/bin/env python3
"""Why neither selection lever moves the > 50 GeV capture (results/wsel_levers.md, todo vt-5t0;
slides/slides_timing_followups.tex).

In the selection failures -- the WAVeS t0 wrong (>= 60 ps) while the hard scatter's cluster (the perfect-selection fix's
choice) carries a right WAVeS in-jet time -- the ratio of WAVeS's score of its pick to its score of the hard-scatter
cluster, factor by factor (left), and the share of each cluster's pT the vertex fit gave PV0 (right):
  sum pT_trk / max(dR, 0.05)   the bracket without the jet pT (what wins)
  pT_jet (bracket-weighted)    what lever 1 (pT_jet^alpha) acts on: 1 when both clusters' tracks share one nearest jet
  exp(-1.5 |dz|)               the cluster z term
Input: the csv of util/scratch/wsel_levers_diag.cxx (git-ignored scratch; its build line is in results/wsel_levers.md).

  cd python && PYTHONNOUSERSITE=1 ~/.venv-hgtd/bin/python gp_wsel_levers_why.py <wsel_selfail.csv> \\
      --out ../figs/tdr/wsel/local_wsel_why
"""
import argparse

import numpy as np

import tdr_root_style as S

ap = argparse.ArgumentParser()
ap.add_argument("csv")
ap.add_argument("--pucat", type=int, default=2, help="the event's hardest forward PU jet: 0 none, 1 30-50 GeV, 2 > 50 GeV")
ap.add_argument("--label", default="local VBF")
ap.add_argument("--out", required=True)
args = ap.parse_args()
R = S.setup()
D = np.genfromtxt(args.csv, delimiter=",", names=True)
D = D[D["pucat"] == args.pucat]
CAT = {0: "no forward pileup jet", 1: "forward pileup jet 30-50 GeV", 2: "forward pileup jet > 50 GeV"}[args.pucat]
print(f"{len(D)} selection failures, {CAT}")

cv = S.keep(R.TCanvas("why", "", 1600, 750))
TS = 0.05
# ── left: the score's factors, WAVeS pick over the hard-scatter cluster (log x) ──
cv.cd()
S.pad("p1", 0.0, 0.0, 0.5, 1.0, top=0.06, bottom=0.15, right=0.04)
R.gPad.SetLogx()
edges = np.logspace(-1, 1, 25)
FAC = [("B0", "#Sigma p_{T}^{trk} / #DeltaR", "kP10Red", 1), ("J", "p_{T}^{jet}", "kP10Violet", 2),
       ("Z", "e^{-1.5|#Deltaz|}", "kP10Blue", 7)]
hists, ymax = [], 0.0
for k, name, col, ls in FAC:
    r = np.clip(D[k + "W"] / D[k + "H"], edges[0] * 1.0001, edges[-1] * 0.9999)   # under/overflow into the edge bins
    h = S.keep(R.TH1D(f"h_{k}", "", len(edges) - 1, edges))
    for v in r: h.Fill(v)
    h.Scale(1.0 / len(D))
    S.style(h, col, 0, ls, 3)
    hists.append((h, name, r))
    ymax = max(ymax, h.GetMaximum())
fr = S.keep(R.gPad.DrawFrame(0.1, 0.0, 10.0, 1.9 * ymax))
S.axes(fr, "WAVeS pick / hard-scatter cluster", "Fraction of events", size=TS, ytoff=1.5, xtoff=1.3)
fr.GetXaxis().SetMoreLogLabels()
one = S.keep(R.TLine(1.0, 0.0, 1.0, 1.1 * ymax)); one.SetLineStyle(2); one.SetLineColor(R.kGray + 2); one.Draw()
leg = S.legend(0.62, 0.50, 0.95, 0.70, size=0.042, margin=0.35)
for h, name, r in hists:
    h.Draw("HIST SAME")
    leg.AddEntry(h, name, "L")
leg.Draw()
ahead = lambda k: 100 * np.mean(D[k + "W"] > D[k + "H"] * (1 + 1e-6))   # 1e-6: one shared jet makes p_T^jet tie up to rounding
S.header(0.20, 0.87, [f"#sqrt{{s}} = 14 TeV, #LT#mu#GT = 200, {args.label}",
                      f"WAVeS selection failures, {CAT}",
                      f"Pick ahead: #Sigma p_{{T}}/#DeltaR {ahead('B0'):.0f}%, p_{{T}}^{{jet}} {ahead('J'):.0f}%, z {ahead('Z'):.0f}%"],
         size=0.042)
S.end_pad()
# ── right: the vertex fit's PV0 share of each cluster ──
cv.cd()
S.pad("p2", 0.5, 0.0, 1.0, 1.0, top=0.06, bottom=0.15, right=0.04)
SH = [("fpvW", "WAVeS pick", "kP10Orange", 1), ("fpvH", "Hard-scatter cluster", "kP10Green", 2)]
hs, ymax = [], 0.0
for k, name, col, ls in SH:
    h = S.keep(R.TH1D(f"h_{k}", "", 10, 0.0, 1.0001))
    for v in D[k]: h.Fill(v)
    h.Scale(1.0 / len(D))
    S.style(h, col, 0, ls, 3)
    hs.append((h, name, D[k]))
    ymax = max(ymax, h.GetMaximum())
fr2 = S.keep(R.gPad.DrawFrame(0.0, 0.0, 1.0, 1.9 * ymax))
S.axes(fr2, "p_{T} share the vertex fit gave PV0", "Fraction of events", size=TS, ytoff=1.5, xtoff=1.3)
leg2 = S.legend(0.45, 0.52, 0.88, 0.68, size=0.042, margin=0.35)
for h, name, v in hs:
    h.Draw("HIST SAME")
    leg2.AddEntry(h, name, "L")
leg2.Draw()
S.header(0.20, 0.87, [f"#sqrt{{s}} = 14 TeV, #LT#mu#GT = 200, {args.label}",
                      f"Median: pick {np.median(D['fpvW']):.2f}, hard-scatter cluster {np.median(D['fpvH']):.2f}",
                      f"Hard-scatter cluster higher in {100 * np.mean(D['fpvH'] > D['fpvW']):.0f}% of events"], size=0.042)
S.end_pad()
S.save(cv, args.out)
