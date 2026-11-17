#!/usr/bin/env python3
"""Correct and misassigned track times against |eta|, as in the ntuples and after the TDR's cleaning -- TDR
Figs. 3.15 / 3.16 on our ntuples, the figure of slides_tdr_timing's "Track times" slide (the user's request,
2026-11-05).

Tracks: pT > 1 GeV, 2.4 < |eta| < 4.0, Track_quality (the analysis's track selection; --all-tracks drops it).
A track's time is
  as in the ntuples   Track_hasValidTime == 1 (>= 1 HGTD hit, hits consistent in time: the TDR's outlier removal,
                      step (c), looks applied upstream -- notes/tdr_presentation.md, "What Track_hasValidTime is")
  cleaned             and the TDR's cuts (a) Track_hitInLastLayer, (b) >= 2 HGTD hits at 3.5 < |eta| < 3.9:
                      util/tdr_common.h's timed()
and is correct when more than 50% of its HGTD hits come from the track's own particle (Track_nHGTDPrimaryHits /
Track_nHGTDHits > 0.5), misassigned otherwise -- the TDR's definition (Sec. 3.2.1, Fig. 3.16). Rates are fractions
of ALL selected tracks in the |eta| bin, as the TDR's "time association rate". Lower panel: the misassigned share
of the timed tracks.

Needs Track_hitInLastLayer, which the skims drop: read the local VBF ntuples or the raw grid ones.

  PYTHONNOUSERSITE=1 ~/.venv-hgtd/bin/python tdr_track_times.py --files 4 --out ../figs/tdr/local_tdr_track_times
"""
import argparse
import glob

import awkward as ak
import numpy as np
import uproot

import tdr_root_style as S

ap = argparse.ArgumentParser()
ap.add_argument("--dir", default="/Users/mcard/project/ntuple-hgtd/")
ap.add_argument("--files", type=int, default=4, help="the first N files of the sorted list")
ap.add_argument("--all-tracks", action="store_true", help="do not require Track_quality")
ap.add_argument("--out", required=True)
args = ap.parse_args()

FILES = sorted(glob.glob(args.dir.rstrip("/") + "/*.root"))[: args.files]
EDGES = np.round(np.arange(2.4, 4.0001, 0.1), 2)
BR = ["Track_pt", "Track_eta", "Track_quality", "Track_hasValidTime", "Track_hitInLastLayer",
      "Track_nHGTDHits", "Track_nHGTDPrimaryHits"]
H = {k: np.zeros(len(EDGES) - 1) for k in ("all", "raw_ok", "raw_bad", "cl_ok", "cl_bad")}
nev = 0
for A in uproot.iterate([f + ":ntuple" for f in FILES], BR, step_size="200 MB"):
    nev += len(A)
    F = {b: ak.to_numpy(ak.flatten(A[b])) for b in BR}
    ae = np.abs(F["Track_eta"])
    sel = (F["Track_pt"] > 1.0) & (ae > 2.4) & (ae < 4.0)             # GeV
    if not args.all_tracks: sel &= F["Track_quality"].astype(bool)
    ae, nh, npr = ae[sel], F["Track_nHGTDHits"][sel], F["Track_nHGTDPrimaryHits"][sel]
    raw = F["Track_hasValidTime"][sel] == 1
    cl = raw & F["Track_hitInLastLayer"][sel].astype(bool) & ~((ae > 3.5) & (ae < 3.9) & (nh < 2))
    ok = npr > 0.5 * nh                                                 # > 50% primary hits
    h = lambda m: np.histogram(ae[m], EDGES)[0]
    H["all"] += h(np.ones(len(ae), bool))
    H["raw_ok"] += h(raw & ok); H["raw_bad"] += h(raw & ~ok)
    H["cl_ok"] += h(cl & ok); H["cl_bad"] += h(cl & ~ok)

tot = {k: v.sum() for k, v in H.items()}
print(f"{len(FILES)} files, {nev:,} events, {tot['all']:,.0f} tracks (pT > 1 GeV, 2.4 < |eta| < 4.0"
      f"{'' if args.all_tracks else ', quality'})")
for lab, a, b in (("as in the ntuples", "raw_ok", "raw_bad"), ("cleaned (a) + (b)", "cl_ok", "cl_bad")):
    print(f"  {lab:18s} correct {tot[a] / tot['all']:6.1%}  misassigned {tot[b] / tot['all']:6.1%}"
          f"  misassigned share of timed {tot[b] / (tot[a] + tot[b]):6.1%}")
print(f"  cleaning keeps {tot['cl_ok'] / tot['raw_ok']:.1%} of the correct times and "
      f"{tot['cl_bad'] / tot['raw_bad']:.1%} of the misassigned ones")
for lo, hi in ((2.4, 3.0), (3.0, 3.5), (3.5, 3.9), (3.9, 4.0)):
    m = (EDGES[:-1] >= lo - 1e-9) & (EDGES[1:] <= hi + 1e-9)
    s = lambda k: H[k][m].sum()
    print(f"  |eta| {lo}-{hi}: misassigned share {s('raw_bad') / (s('raw_ok') + s('raw_bad')):.1%} -> "
          f"{s('cl_bad') / (s('cl_ok') + s('cl_bad')):.1%};  correct rate {s('raw_ok') / s('all'):.1%} -> {s('cl_ok') / s('all'):.1%}")

R = S.setup()
x = 0.5 * (EDGES[:-1] + EDGES[1:])
rate = lambda k: H[k] / H["all"]
share = lambda a, b: H[b] / (H[a] + H[b])
GREEN, MAG = "#2E8B57", "kP10Violet"   # the TDR's green (correct) and magenta (misassigned)
MAIN, RAT = 0.68, 0.32
TS, TSR = 0.05 * 600 / (800 * MAIN), 0.05 * 600 / (800 * RAT)
cv = S.keep(R.TCanvas("trktimes", "", 800, 800))
cv.cd(); S.pad("tt_m", 0.0, RAT, 1.0, 1.0, top=0.05 / MAIN * 0.6, bottom=0.02)
fr = S.keep(R.gPad.DrawFrame(2.4, 0.0, 4.0, 1.25))
S.axes(fr, "", "Time association rate", size=TS, ytoff=1.35)
fr.GetXaxis().SetLabelSize(0)
leg = S.legend(0.19, 0.555, 0.92, 0.72, size=0.62 * TS, ncol=2, margin=0.22)
for lab, k, colour, mk, ls in (("Correct, as in the ntuples", "raw_ok", GREEN, 24, 2),
                               ("Correct, cleaned", "cl_ok", GREEN, 20, 1),
                               ("Misassigned, as in the ntuples", "raw_bad", MAG, 25, 2),
                               ("Misassigned, cleaned", "cl_bad", MAG, 21, 1)):
    g = S.style(S.graph(x, rate(k)), colour, mk, ls, 2, 1.1); g.Draw("LP SAME")
    leg.AddEntry(g, lab, "LP")
leg.Draw()
S.header(0.19, 0.88, ["VBF H#rightarrowinv., #LT#mu#GT = 200, tracks p_{T} > 1 GeV",
                      "Correct: > 50% of the HGTD hits from the track's particle",
                      "Cleaned: last ITk hit on the last layer; #geq 2 hits at 3.5 < |#eta| < 3.9"], size=0.8 * TS)
R.gPad.RedrawAxis()
cv.cd(); S.pad("tt_r", 0.0, 0.0, 1.0, RAT, top=0.03, bottom=0.36)
fr2 = S.keep(R.gPad.DrawFrame(2.4, 0.05, 4.0, 0.58))   # from 0.05: a "0" label would sit on the x axis's "2.4"
S.axes(fr2, "Track |#eta|", "Misassigned / timed", size=TSR, ytoff=1.35 * TS / TSR, xtoff=1.2, ndiv_y=503)
fr2.GetYaxis().SetTitleSize(0.8 * TSR); fr2.GetYaxis().SetTitleOffset(1.35 * TS / (0.8 * TSR))
for k, mk, ls in ((("raw_ok", "raw_bad"), 25, 2), (("cl_ok", "cl_bad"), 21, 1)):
    S.style(S.graph(x, share(*k)), MAG, mk, ls, 2, 1.0).Draw("LP SAME")
R.gPad.RedrawAxis()
S.save(cv, args.out)
