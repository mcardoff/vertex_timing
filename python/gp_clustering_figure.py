#!/usr/bin/env python3
"""The clustering and the WAVeS choice in one event, for the Global Performance talk's algorithm slide
(slides/slides_gp_timing.tex): the C++ clusters of python/clustering_dump.cxx (talk inputs: TDR-cleaned times, 2.5
sigma z association), drawn in the ATLAS style.

  upper   the event's forward tracks with a time at (t, pT), coloured by cluster; filled = hard scatter, open =
          pileup (truth); horizontal bars = sigma_t; each cluster's time +- sigma as a band; truth t_HS dashed,
          the WAVeS t0 solid (and Athena's t0 dotted when it has one)
  lower   each cluster's WAVeS score at its time; the chosen cluster outlined, the one with the largest pT sum marked

  PYTHONNOUSERSITE=1 ~/.venv-hgtd/bin/python gp_clustering_figure.py --entry 6087 --out ../figs/tdr/gp/gp_clustering_event
"""
import argparse
import json
import os
import subprocess

import numpy as np

import tdr_root_style as S

ap = argparse.ArgumentParser()
ap.add_argument("--file", default="/Users/mcard/project/ntuple-hgtd/user.mcardiff.45809429.Output._000001.SuperNtuple.root")
ap.add_argument("--entry", type=int, required=True)
ap.add_argument("--out", required=True)
ap.add_argument("--label", default="VBF H#rightarrowinv.")
args = ap.parse_args()

js = f"{args.out}.json"
if not os.path.exists(js):
    subprocess.run(["root", "-l", "-b", "-q", f'clustering_dump.cxx("{args.file}", {args.entry}, "{js}")'],
                   check=True, capture_output=True, text=True)
E = json.load(open(js))
C = E["clusters"]
order = np.argsort([c["t"] for c in C])                      # colours left to right in time
PAL = ["#3F90DA", "#FFA90E", "#BD1F01", "#94A4A2", "#832DB6", "#A96B59", "#E76300", "#B9AC70", "#717581", "#92DADD"]
iW, iP, tHS = E["i_waves"], E["i_sumpt"], E["t_hs"]
ts = [tr["t"] for c in C for tr in c["tracks"]]
x0, x1 = min(ts) - 60, max(ts) + 60
pts = [tr["pt"] for c in C for tr in c["tracks"]]
print(f"entry {E['entry']}: {len(C)} clusters, t_HS {tHS:.0f} ps, WAVeS t0 {E['t0_waves']:.0f} ps, "
      f"Athena {E['athena_t0']:.0f} ps (valid {E['athena_valid']}); chosen {iW} (HS cluster {E['i_hs']}), largest pT sum {iP}")
for k in order:
    c = C[k]
    print(f"  cluster {k}: t {c['t']:7.1f} +- {c['sigma']:4.1f} ps, {len(c['tracks'])} tracks, sum pT {c['sumpt']:5.1f}, "
          f"purity {c['purity']:.2f}, WAVeS {c['waves']:.3g}")

R = S.setup()
col = {int(k): S.colour(PAL[i % len(PAL)]) for i, k in enumerate(order)}
def pale(hexc, a):   # the colour mixed with white: ROOT's PNG output drops fill transparency
    h = hexc.lstrip("#"); rgb = [int(h[j:j + 2], 16) for j in (0, 2, 4)]
    return R.TColor.GetColor(*[int(round(255 - a * (255 - v))) for v in rgb])
band = {int(k): pale(PAL[i % len(PAL)], 0.22) for i, k in enumerate(order)}
band_sel = {int(k): pale(PAL[i % len(PAL)], 0.45) for i, k in enumerate(order)}
cv = S.keep(R.TCanvas("clev", "", 1000, 860))
# ---- upper: the tracks
cv.cd(); S.pad("u", 0, 0.34, 1, 1, top=0.17, bottom=0.02, left=0.12, right=0.04)
R.gPad.SetLogy()
y0, y1 = 0.7, max(pts) * 3.0
fr = S.keep(R.gPad.DrawFrame(x0, y0, x1, y1))
S.axes(fr, "", "Track p_{T} [GeV]", size=0.055, ytoff=0.95)
fr.GetXaxis().SetLabelSize(0)
for k in order:                                               # cluster bands under the points
    c = C[k]
    b = S.keep(R.TBox(c["t"] - c["sigma"], y0, c["t"] + c["sigma"], y1))
    b.SetFillColor(band_sel[k] if k == iW else band[k]); b.SetLineWidth(0); b.Draw()
for t, sty, w, c_ in ((tHS, 2, 2, R.kBlack), (E["t0_waves"], 1, 3, S.colour(S.COL["waves"]))):
    l = S.keep(R.TLine(t, y0, t, y1)); l.SetLineStyle(sty); l.SetLineWidth(w); l.SetLineColor(c_); l.Draw()
if E["athena_valid"] and x0 < E["athena_t0"] < x1:
    l = S.keep(R.TLine(E["athena_t0"], y0, E["athena_t0"], y1)); l.SetLineStyle(3); l.SetLineWidth(3)
    l.SetLineColor(S.colour(S.COL["hgtd"])); l.Draw()
for k in order:
    for hs in (1, 0):
        tr = [t for t in C[k]["tracks"] if t["hs"] == hs]
        if not tr: continue
        g = S.keep(R.TGraphErrors(len(tr)))
        for i, t in enumerate(tr): g.SetPoint(i, t["t"], t["pt"]); g.SetPointError(i, t["res"], 0)
        g.SetMarkerStyle(20 if hs else 24); g.SetMarkerSize(1.3); g.SetMarkerColor(col[k]); g.SetLineColor(col[k])
        g.SetLineWidth(2); g.Draw("P SAME")
leg = S.legend(0.56, 0.835, 0.97, 0.995, size=0.033, ncol=2, margin=0.22)
mh = S.keep(R.TGraph(1)); mh.SetMarkerStyle(20); mh.SetMarkerSize(1.3)
mp = S.keep(R.TGraph(1)); mp.SetMarkerStyle(24); mp.SetMarkerSize(1.3)
bx = S.keep(R.TH1D("bx", "", 1, 0, 1)); bx.SetFillColor(pale(PAL[0], 0.45)); bx.SetLineWidth(0)
lt = S.keep(R.TLine()); lt.SetLineStyle(2); lt.SetLineWidth(2)
lw = S.keep(R.TLine()); lw.SetLineWidth(3); lw.SetLineColor(S.colour(S.COL["waves"]))
leg.AddEntry(mh, "hard-scatter track", "P"); leg.AddEntry(lt, "truth t_{HS}", "L")
leg.AddEntry(mp, "pileup track", "P"); leg.AddEntry(lw, "WAVeS t_{0}", "L")
leg.AddEntry(bx, "cluster t #pm #sigma", "F")
if E["athena_valid"] and x0 < E["athena_t0"] < x1:
    la = S.keep(R.TLine()); la.SetLineStyle(3); la.SetLineWidth(3); la.SetLineColor(S.colour(S.COL["hgtd"]))
    leg.AddEntry(la, "Athena t_{0}", "L")
leg.Draw()
S.header(0.12, 0.945, [f"#sqrt{{s}} = 14 TeV, #LT#mu#GT = 200, {args.label}", "forward tracks with a time, coloured by cluster"],
         size=0.038)
R.gPad.RedrawAxis()
# ---- lower: the WAVeS score per cluster
cv.cd(); S.pad("d", 0, 0, 1, 0.34, top=0.03, bottom=0.30, left=0.12, right=0.04)
R.gPad.SetLogy()
sc = [max(c["waves"], 1e-3) for c in C]
lo, hi = min(sc) / 5, max(sc) * 40
fd = S.keep(R.gPad.DrawFrame(x0, lo, x1, hi))
S.axes(fd, "Track time [ps]", "WAVeS score", size=0.11, ytoff=0.45, xtoff=1.05, ndiv_y=503)
for k in order:
    c = C[k]; w = max(c["sigma"], 6.0)
    b = S.keep(R.TBox(c["t"] - w, lo, c["t"] + w, sc[k])); b.SetFillColor(col[k]); b.SetLineColor(R.kBlack)
    b.SetLineWidth(3 if k == iW else 1); b.Draw("L")
    if k == iW:
        t = S.keep(R.TLatex(c["t"], sc[k] * 2.0, "chosen")); t.SetTextAlign(21); t.SetTextSize(0.10); t.Draw()
    if k == iP and iP != iW:
        t = S.keep(R.TLatex(c["t"], sc[k] * 2.0, "largest #Sigmap_{T}")); t.SetTextAlign(21); t.SetTextSize(0.10)
        t.SetTextColor(R.kGray + 2); t.Draw()
R.gPad.RedrawAxis()
S.save(cv, args.out)
