#!/usr/bin/env python3
"""The captured t0, and what leaving the jet out does to it (results/lojo_t0.md; slides/slides_timing_followups.tex).

For the forward pileup jets the time gate has to remove -- those passing the ITk-only R_pT cut at 85% hard-scatter
efficiency, per pT bin (2.4 < |eta| < 4.0, 30-50 and > 50 GeV) -- the fraction of their R_pT the gate keeps (gated
R_pT / ITk-only R_pT), with the event WAVeS t0, the two leave-one-jet-out t0s and the truth t0. A jet the gate keeps
whole (the last bin) sits at the time of the t0 that judged it: with the event t0, its own collision's time (the
capture). Unweighted, as tdr_pujet_plots.py; needs a tdr_pujet --lojo tree.

  cd python && PYTHONNOUSERSITE=1 ~/.venv-hgtd/bin/python gp_lojo_capture.py \\
      ../condor/vbf/vbf_lojo_rptptmaxnone_tdr_pujet.root --out ../figs/tdr/lojo/vbf_lojo_capture
"""
import argparse

import numpy as np
import uproot

import tdr_root_style as S

ap = argparse.ArgumentParser()
ap.add_argument("tree")
ap.add_argument("--label", default="VBF H#rightarrowinv.")
ap.add_argument("--out", required=True)
args = ap.parse_args()
R = S.setup()
J = uproot.open(args.tree)["jets"].arrays(["pt", "eta", "hs", "pu", "rpt_itk", "rpt_t0_waves", "rpt_t0_waves_lojet",
                                           "rpt_t0_waves_lojo", "ninjet_lojo", "rpt_t0_truth"], library="np")
J["guarded"] = np.where(J["ninjet_lojo"] >= 2, J["rpt_t0_waves_lojo"], J["rpt_t0_waves"])
SRC = [("rpt_t0_waves", "WAVeS t_{0}", S.COL["waves"], 1), ("rpt_t0_waves_lojet", "Jet left out", "kP10Green", 1),
       ("guarded", "Guarded", "kP10Violet", 1), ("rpt_t0_truth", "Truth t_{0}", S.COL["truth"], 3)]


def itk_cut(hs, pu, eff=0.85):
    """the highest ITk-only R_pT threshold keeping >= eff of the hard-scatter jets (tdr_pujet_plots.working_point)"""
    th = np.unique(np.concatenate([hs, pu]))
    e_hs = 1.0 - np.searchsorted(np.sort(hs), th, side="left") / len(hs)
    return float(th[np.flatnonzero(e_hs >= eff)[-1]])


ae = np.abs(J["eta"])
cv = S.keep(R.TCanvas("cap", "", 1600, 700))
for p, (lab, lo, hi) in enumerate((("30 < p_{T} < 50 GeV", 30, 50), ("p_{T} > 50 GeV", 50, 1e9))):
    m = (J["pt"] > lo) & (J["pt"] < hi) & (ae >= 2.4) & (ae < 4.0)
    hs, pu = m & (J["hs"] == 1), m & (J["pu"] == 1)
    cut = itk_cut(J["rpt_itk"][hs], J["rpt_itk"][pu])
    sel = pu & (J["rpt_itk"] >= cut)
    cv.cd(); S.pad(f"p{p}", 0.5 * p, 0.0, 0.5 * p + 0.5, 1.0, top=0.05, bottom=0.14, left=0.14, right=0.04)
    fr = S.keep(R.gPad.DrawFrame(0.0, 0.0, 1.0, 1.0))
    S.axes(fr, "Fraction of the jet's R_{pT} the t_{0} gate keeps", "Fraction of these pileup jets", size=0.045,
           ytoff=1.4, xtoff=1.3)
    leg = S.legend(0.17, 0.52, 0.80, 0.72, size=0.036, margin=0.12)
    print(f"{lab.replace('_{T}', 'T')}: {sel.sum()} pileup jets pass ITk-only R_pT >= {cut:.4f} (85% HS)")
    for col, name, colour, ls in SRC:
        f = np.clip(J[col][sel] / J["rpt_itk"][sel], 0.0, 0.99999)
        h = S.keep(R.TH1D(f"h{p}{col}", "", 10, 0.0, 1.0))
        for x in f: h.Fill(x)
        h.Scale(1.0 / max(h.Integral(), 1))
        S.style(h, colour, 0, ls, 3); h.Draw("HIST SAME")
        whole = 100 * np.mean(f >= 0.9)
        leg.AddEntry(h, f"{name}: keeps #geq 90% of {whole:.0f}%", "L")
        print(f"  {name.replace('_{0}', '0'):14s} keeps >= 90% of the R_pT of {whole:.1f}%, < 10% of {100 * np.mean(f < 0.1):.1f}%")
    leg.Draw()
    S.header(0.17, 0.88, [f"#sqrt{{s}} = 14 TeV, #LT#mu#GT = 200, {args.label}",
                          f"{lab}, 2.4 < |#eta| < 4.0: {sel.sum():,} pileup jets", "passing the ITk-only cut at 85% HS efficiency"],
             size=0.042)
    R.gPad.RedrawAxis()
    S.end_pad()
S.save(cv, args.out)
