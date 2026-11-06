#!/usr/bin/env python3
"""Why jet vs jet trails the jet-vs-t0 veto, on the forward-forward pairs where both can act (the user, 2026-11-06:
"for the forward-forward cases ... a plot of the pull distribution between the two jets with a time, just to demonstrate
this point").

Events: tdr_fwdpair --pair=one --jvt=loose trees, forward-forward pairs (cat 2) with both VBS jets timed (core time from
>= 1 track), optionally after the kinematic BDT at a signal efficiency (--bdt-eff, needs --bdt). Pulls as the veto's
(tdr_fwdpair_common, jet-time calibration from the VBF hard-scatter legs):
  left   jet vs jet    |t_1 - t_2| / sqrt(s_1^2 + s_2^2)
  right  jet vs t0     the larger of the two jets' |t_jet - t0| / sqrt(s_jet^2 + (f sigma_t0)^2), WAVeS t0
         (0 where the event has no t0: the veto keeps it)
By the pair's truth labels: two hard-scatter jets from VBF H->inv. (Z+jets has 42 such events), one or two pileup jets
from Z+jets. Unit area per category, weighted; the dashed line is the 3 sigma veto; the legend gives the fraction vetoed.

  PYTHONNOUSERSITE=1 ~/.venv-hgtd/bin/python tdr_fwdpair_pulls.py --sig ../condor/vbf/vbf_onefwd_jvtLoose_tdr_fwdpair.root \\
      --bkg ../condor/zjets/zjets_onefwd_jvtLoose_tdr_fwdpair.root --out ../figs/tdr/grid/onefwd_jvt/grid_onefwd_jvt_ff_pulls
"""
import argparse

import numpy as np
import uproot

import tdr_fwdpair_common as FC
import tdr_root_style as S

ap = argparse.ArgumentParser()
ap.add_argument("--sig", required=True); ap.add_argument("--bkg", required=True); ap.add_argument("--out", required=True)
ap.add_argument("--bdt", help="the kinematic BDT's scores (tdr_fwdpair_kinbdt.py npz)")
ap.add_argument("--bdt-eff", type=float, help="cut the BDT at this signal efficiency")
ap.add_argument("--nsigma", type=float, default=3.0)
args = ap.parse_args()
R = S.setup()
E = {"vbf": uproot.open(args.sig)["events"].arrays(library="np"), "zjets": uproot.open(args.bkg)["events"].arrays(library="np")}
f = FC.calibration(E["vbf"], verbose=False)
sel = {s: (e["cat"] == 2) & (e["j1_core_n"] >= 1) & (e["j2_core_n"] >= 1) for s, e in E.items()}
cutline = "No BDT"
if args.bdt_eff:
    D = np.load(args.bdt); sc = {"vbf": D["sig_score"], "zjets": D["bkg_score"]}
    w0 = E["vbf"]["weight"].astype(float); o = np.argsort(-sc["vbf"]); cs = np.cumsum(w0[o]) / w0.sum()
    thr = sc["vbf"][o][np.searchsorted(cs, args.bdt_eff)]
    sel = {s: m & (sc[s] >= thr) for s, m in sel.items()}
    cutline = f"BDT at #varepsilon_{{S}} = {100 * args.bdt_eff:.0f}%"
CATS = [("HS + HS, VBF H#rightarrowinv.", "vbf", 2, 0, "#2E8B57", 1),
        ("HS + PU, Z+jets", "zjets", 1, 1, S.COL["hgtd"], 1),
        ("PU + PU, Z+jets", "zjets", 0, 2, S.COL["trkptz"], 1)]
NB, XMAX = 40, 20.0
cv = S.keep(R.TCanvas("ffpulls", "", 1600, 600))
for side, (key, xtitle) in enumerate((("jj", "Jet vs jet pull  |t_{1} #minus t_{2}| / #sqrt{#sigma_{1}^{2} + #sigma_{2}^{2}}"),
                                       ("t0", "Jet vs WAVeS t_{0} pull (the larger jet's)"))):
    cv.cd(); S.pad(f"pp{side}", 0.5 * side, 0.0, 0.5 * side + 0.5, 1.0, top=0.05, bottom=0.16, left=0.13, right=0.04)
    R.gPad.SetLogy()
    fr = S.keep(R.gPad.DrawFrame(0.0, 1.5e-4, XMAX, 300.0))   # headroom on the log scale: header, then the legend, then curves
    S.axes(fr, xtitle, "Fraction of events / 0.5", size=0.05, ytoff=1.15, xtoff=1.3)
    ln = S.keep(R.TLine(args.nsigma, 1.5e-4, args.nsigma, 0.6)); ln.SetLineStyle(2); ln.SetLineColor(R.kGray + 2); ln.SetLineWidth(2); ln.Draw()
    leg = S.legend(0.40, 0.55, 0.94, 0.70, size=0.040, margin=0.13)
    print(f"\n{key} pull, forward-forward pairs, both jets timed, {cutline.replace('#varepsilon_{S}', 'eps_S')}:")
    for name, s, nhs, npu, colour, ls in CATS:
        e = E[s]; w = e["weight"].astype(float)
        m = sel[s] & (e["j1_hs"] + e["j2_hs"] == nhs) & (e["j1_pu"] + e["j2_pu"] == npu)
        P = FC.pulls(e, f)
        x = P["jj"] if key == "jj" else np.fmax(np.nan_to_num(P["waves_1"], nan=0.0), np.nan_to_num(P["waves_2"], nan=0.0))
        x = np.nan_to_num(x, nan=0.0)[m]; ww = w[m]
        h = S.keep(R.TH1D(f"h_{key}_{nhs}{npu}", "", NB, 0.0, XMAX))
        for v, wv in zip(np.minimum(x, XMAX - 1e-6), ww): h.Fill(v, wv)
        h.Scale(1.0 / h.Integral())
        h.SetLineColor(S.colour(colour)); h.SetLineWidth(3); h.SetLineStyle(ls); h.Draw("HIST SAME")
        frac = ww[x >= args.nsigma].sum() / ww.sum()
        leg.AddEntry(h, f"{name}: {100 * frac:.0f}% vetoed", "L")
        print(f"  {name.replace('#rightarrow', '->'):30s} {m.sum():7d} MC events, vetoed at {args.nsigma:g} sigma: {frac:.1%}")
    leg.Draw()
    lines = ["#sqrt{s} = 14 TeV, HL-LHC, #LT#mu#GT = 200", "Forward-forward pairs, both jets timed", f"JVT + fJVT loose, {cutline}"] if side == 0 else \
            [f"{args.nsigma:g}#sigma veto: right of the dashed line", "Last bin: overflow"]
    S.header(0.17, 0.88, lines, size=0.045)
    R.gPad.RedrawAxis()
S.save(cv, args.out)
