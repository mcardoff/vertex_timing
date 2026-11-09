#!/usr/bin/env python3
"""Why WAVeS beats TRKPTZ on VBF and trails it on Z+jets, for the Global Performance talk (slides/slides_gp_timing.tex):
the core-fraction difference WAVeS - TRKPTZ split exactly into its two parts, event by event.

  selection   WAVeS's chosen cluster (its full-cluster time, dt_waves_full_clean) against TRKPTZ's -- the jet-weighted
              score picking a different cluster
  re-timing   WAVeS's time (recomputed from the cluster's tracks near forward jets) against its own cluster's full time
total = selection + re-timing, per bin. Upper panels: the three against the number of forward hard-scatter tracks with a
time; the text gives the inclusive gains and losses of each part. util/tdr_t0 trees (commit c413fbc), cleaned times.

  PYTHONNOUSERSITE=1 ~/.venv-hgtd/bin/python gp_waves_vs_trkptz.py --trees ../condor/vbf/vbf_tdr_t0_v3.root \\
      ../condor/zjets/zjets_tdr_t0_v3.root --out ../figs/tdr/gp/gp_waves_vs_trkptz
"""
import argparse

import numpy as np
import uproot

import tdr_root_style as S

ap = argparse.ArgumentParser()
ap.add_argument("--trees", nargs=2, required=True, metavar=("VBF", "ZJETS"))
ap.add_argument("--out", required=True)
args = ap.parse_args()
FOLD, PASS = 15, 60.0
R = S.setup()
cv = S.keep(R.TCanvas("wvt", "", 1400, 700))
for c, (path, lab) in enumerate(zip(args.trees, ("VBF H#rightarrowinv.", "Z+jets"))):
    E = uproot.open(path)["events"].arrays(library="np")
    w = E["weight"].astype(float); W = w.sum()
    T, F, V = (np.abs(E[k]) < PASS for k in ("dt_trkptz_clean", "dt_waves_full_clean", "dt_waves_clean"))
    n = np.minimum(E["n_hs_trk"], FOLD)
    pc = lambda m: 100 * w[m].sum() / W
    sel_g, sel_l = pc(F & ~T), pc(T & ~F)
    ret_g, ret_l = pc(V & ~F), pc(F & ~V)
    print(f"{lab.replace('#rightarrow', '->')}: total {pc(V) - pc(T):+.2f} = selection {sel_g - sel_l:+.2f} (+{sel_g:.2f} / -{sel_l:.2f}) "
          f"+ re-timing {ret_g - ret_l:+.2f} (+{ret_g:.2f} / -{ret_l:.2f})")
    xs = np.arange(FOLD + 1)
    d = lambda a, b: np.array([100 * (w[(n == k) & a].sum() - w[(n == k) & b].sum()) / max(w[n == k].sum(), 1e-12) for k in xs])
    curves = [("total: WAVeS #minus TRKPTZ", d(V, T), "#000000", 20, 1),
              ("selection: the jet-weighted choice", d(F, T), "kP10Blue", 21, 1),
              ("re-timing: the in-jet time", d(V, F), "kP10Yellow", 22, 1)]
    cv.cd(); S.pad(f"wv_{c}", 0.5 * c, 0.0, 0.5 * c + 0.5, 1.0, top=0.06, bottom=0.13, left=0.14, right=0.04)
    lo = min(-8.0, min(np.min(y) for _, y, *_ in curves) - 1); hi = max(8.0, max(np.max(y) for _, y, *_ in curves) + 1)
    fr = S.keep(R.gPad.DrawFrame(-0.5, lo, FOLD + 0.5, hi + 0.9 * (hi - lo)))
    S.axes(fr, f"Forward hard-scatter tracks with a time (last bin: #geq {FOLD})", "Core fraction difference [points]", size=0.042, ytoff=1.4, xtoff=1.2)
    zero = S.keep(R.TLine(-0.5, 0, FOLD + 0.5, 0)); zero.SetLineStyle(2); zero.SetLineColor(R.kGray + 2); zero.Draw()
    leg = S.legend(0.18, 0.58, 0.62, 0.72, size=0.033)
    for name, y, col, mk, ls in curves:
        g = S.style(S.graph(xs, y), col, mk, ls, 2, 1.1); g.Draw("LP SAME")
        leg.AddEntry(g, name, "LP")
    leg.Draw()
    S.header(0.18, 0.88, ["#sqrt{s} = 14 TeV, #LT#mu#GT = 200",
                          f"{lab}: all events {pc(V) - pc(T):+.1f} = selection {sel_g - sel_l:+.1f} + re-timing {ret_g - ret_l:+.1f}",
                          f"selection: +{sel_g:.1f} / #minus{sel_l:.1f};  re-timing: +{ret_g:.1f} / #minus{ret_l:.1f}  (fixed / broken)"], size=0.033)
    R.gPad.RedrawAxis()
S.save(cv, args.out)
