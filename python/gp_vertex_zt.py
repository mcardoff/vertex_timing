#!/usr/bin/env python3
"""Why time: one bunch crossing's interactions in the z-t plane, for the Global Performance talk's slide 3
(slides/slides_gp_timing.tex). A new version of the HGTD "Vertex t vs Vertex z" plot (ATLAS Simulation Preliminary,
HGTD talks around 2019) from our own local VBF sample, which keeps every truth vertex (~200 per event; the grid
samples keep only the hard scatter's).

  ellipses   the truth interactions at (z, t), +-1 mm x +-30 ps: a 1 GeV forward track's z0 resolution (|eta| ~ 2.5,
             HGTD TDR Fig. 2.6) and HGTD's per-track time resolution -- in z such a track is compatible with several
             interactions, in time with one
  dashed     the reconstructed vertices' z (ITk only)
  red        the hard scatter
The event is the one, among the first --scan events, with the most pileup interactions within 1 mm of the hard
scatter in z (a crowded one, as the point of the figure).

  PYTHONNOUSERSITE=1 ~/.venv-hgtd/bin/python gp_vertex_zt.py --out ../figs/tdr/gp/gp_vertex_zt
"""
import argparse
import glob

import numpy as np
import uproot

import tdr_root_style as S

ap = argparse.ArgumentParser()
ap.add_argument("--ntuples", default="/Users/mcard/project/ntuple-hgtd/")
ap.add_argument("--scan", type=int, default=400, help="events scanned for the most crowded hard scatter")
ap.add_argument("--sigz", type=float, default=1.0, help="ellipse half-width in z [mm]")
ap.add_argument("--sigt", type=float, default=30.0, help="ellipse half-height in t [ps]")
ap.add_argument("--out", required=True)
args = ap.parse_args()

f = sorted(glob.glob(args.ntuples.rstrip("/") + "/*.root"))[0]
A = uproot.open(f)["ntuple"].arrays(["TruthVtx_z", "TruthVtx_time", "TruthVtx_isHS", "RecoVtx_z"],
                                     entry_stop=args.scan, library="np")
best, bestn = 0, -1
for i in range(len(A["TruthVtx_z"])):
    z, t = A["TruthVtx_z"][i], A["TruthVtx_time"][i]
    if len(z) < 150 or abs(z[0]) > 40 or abs(t[0]) > 250: continue
    n = int(np.sum(np.abs(z[1:] - z[0]) < 1.0))
    if n > bestn: best, bestn = i, n
z, t, rz = A["TruthVtx_z"][best], A["TruthVtx_time"][best], A["RecoVtx_z"][best]
hs = 0
near = np.abs(z - z[hs]) < 1.0; near[hs] = False
print(f"{f.split('/')[-1]} entry {best}: {len(z)} truth interactions, {len(rz)} reco vertices; hard scatter at "
      f"z = {z[hs]:.2f} mm, t = {t[hs]:.0f} ps; {near.sum()} pileup interactions within 1 mm in z, their |dt| to it: "
      + ", ".join(f"{abs(v):.0f}" for v in (t[near] - t[hs])) + " ps")

R = S.setup()
# the original's layout (user, 2026-11-11): no legend or header, only the ATLAS label, top right inside the frame; the
# axis runs to 500 ps with the lines stopping at 400, which leaves the label clear of the interactions
cv = S.keep(R.TCanvas("zt", "", 1000, 580))   # wide: stacked under TDR Fig. 2.6 on slide 3, height is the limit
S.pad("p", 0, 0, 1, 1, top=0.05, bottom=0.16, left=0.10, right=0.03)
fr = S.keep(R.gPad.DrawFrame(-100, -400, 100, 500))
S.axes(fr, "Vertex z [mm]", "Vertex t [ps]", size=0.058, ytoff=0.8, xtoff=1.1)
fr.GetYaxis().SetNdivisions(509)
for zz in rz:                                       # reconstructed vertices: z only
    if -100 < zz < 100:
        l = S.keep(R.TLine(zz, -400, zz, 400)); l.SetLineStyle(2); l.SetLineColor(R.kBlack); l.SetLineWidth(1); l.Draw()
PU_FILL, PU_LINE, HS_COL = R.TColor.GetColor("#7C93D9"), R.TColor.GetColor("#1F3A93"), R.TColor.GetColor("#D62728")
lh = S.keep(R.TLine(z[hs], -400, z[hs], 400)); lh.SetLineColor(HS_COL); lh.SetLineWidth(2); lh.Draw()   # under the ellipses
for i in range(len(z)):
    if i == hs or not (-100 < z[i] < 100 and -400 < t[i] < 400): continue
    e = S.keep(R.TEllipse(z[i], t[i], args.sigz, args.sigt)); e.SetFillColor(PU_FILL); e.SetLineColor(PU_LINE); e.SetLineWidth(1); e.Draw()
eh = S.keep(R.TEllipse(z[hs], t[hs], args.sigz, args.sigt)); eh.SetFillColor(HS_COL); eh.SetLineColor(R.kBlack); eh.SetLineWidth(2); eh.Draw()
S.header(0.585, 0.875, [], size=0.058)
R.gPad.RedrawAxis()
S.save(cv, args.out)
