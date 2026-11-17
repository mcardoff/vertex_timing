#!/usr/bin/env python3
"""The t0 pull behind the R_pT gate's calibration factor f, for the Global Performance talk's backup
(slides/slides_gp_timing.tex, linked from "The R_pT tagger, with a t0 gate").

(t0 - t_HS) / (f sigma_t0) for the WAVeS and Athena t0s of the jet section's events (util/tdr_fwdpair --pair=one trees:
the talk's t0s, cleaned track times, 2.5 sigma association), with f the per-sample inflation the gate uses
(util/rpt_v5_common.h's table, read from the tree), against a unit Gaussian. If f is right the core is a unit Gaussian.
Weighted (generator weights). Each legend gives the core width fitted to |pull| < 2 (a Gaussian), the share beyond
3 and the median f sigma_t0. The style follows python/vbs_time_veto_plot.py's t0_figure (the time-veto study's version,
which also drew TRKPTZ and predates the cleaning).

  PYTHONNOUSERSITE=1 ~/.venv-hgtd/bin/python gp_t0_pull.py --vbf ../condor/vbf/vbf_onefwd_tdr_fwdpair.root \\
      --zjets ../condor/zjets/zjets_onefwd_tdr_fwdpair.root --out ../figs/tdr/gp/gp_t0_pull
"""
import argparse

import numpy as np
import uproot

import tdr_root_style as S

ap = argparse.ArgumentParser()
ap.add_argument("--vbf", required=True)
ap.add_argument("--zjets", required=True)
ap.add_argument("--out", required=True)
args = ap.parse_args()
R = S.setup()
SRC = [("waves", "WAVeS t_{0}", S.COL["waves"], 1), ("hgtd", "Athena t_{0}", S.COL["hgtd"], 7)]
cv = S.keep(R.TCanvas("pull", "", 1600, 650))
for p, (lab, path) in enumerate((("VBF H#rightarrowinv.", args.vbf), ("Z+jets", args.zjets))):
    E = uproot.open(path)["events"].arrays(library="np")
    w = E["weight"].astype(float)
    cv.cd(); S.pad(f"p{p}", 0.5 * p, 0.0, 0.5 * p + 0.5, 1.0, top=0.05, bottom=0.14, left=0.14, right=0.04)
    R.gPad.SetLogy()
    fr = S.keep(R.gPad.DrawFrame(-8, 2e-5, 8, 300))
    S.axes(fr, "(t_{0} - t_{HS}) / (f #sigma_{t_{0}})   (overflow in the edge bins)", "Fraction of events with that t_{0}",
           size=0.045, ytoff=1.4, xtoff=1.3)
    leg = S.legend(0.17, 0.60, 0.93, 0.76, size=0.031)
    for src, name, col, ls in SRC:
        ok = E[f"ok_{src}"].astype(bool) & (E[f"sig0_{src}"] > 0)
        sg = E[f"infl_{src}"][ok] * E[f"sig0_{src}"][ok]
        pull = (E[f"t0_{src}"][ok] - E["t0_truth"][ok]) / sg
        ww = w[ok]
        h = S.keep(R.TH1D(f"h{p}{src}", "", 64, -8, 8))
        for x, wi in zip(np.clip(pull, -7.99, 7.99), ww): h.Fill(x, wi)
        h.Scale(1.0 / h.Integral())
        S.style(h, col, 0, ls, 3); h.Draw("HIST SAME")
        g = R.TF1(f"g{p}{src}", "gaus", -2, 2); h.Fit(g, "RQ0N")
        core = abs(g.GetParameter(2))
        beyond = 100 * np.sum(ww[np.abs(pull) >= 3]) / np.sum(ww)
        f = float(np.unique(E[f"infl_{src}"][ok])[0])
        leg.AddEntry(h, f"{name}: f = {f:.2f};  core #sigma = {core:.2f};  {beyond:.0f}% beyond 3", "L")
        print(f"{lab.replace('#rightarrow', '->'):14s} {src:5s}: f {f:.2f}, fitted core sigma {core:.3f}, {beyond:.1f}% beyond 3, "
              f"median f*sigma {np.median(sg):.1f} ps, valid in {100 * ok.mean():.0f}% of {len(w):,} events")
    u = S.keep(R.TF1(f"u{p}", "[0]*exp(-0.5*x*x)", -8, 8)); u.SetParameter(0, 0.25 / np.sqrt(2 * np.pi))
    u.SetLineColor(R.kBlack); u.SetLineStyle(2); u.SetLineWidth(2); u.Draw("SAME")
    leg.AddEntry(u, "unit Gaussian", "L")
    leg.Draw()
    S.header(0.17, 0.88, [f"#sqrt{{s}} = 14 TeV, #LT#mu#GT = 200, {lab}: the jet section's events"], size=0.042)
    R.gPad.RedrawAxis()
S.save(cv, args.out)
