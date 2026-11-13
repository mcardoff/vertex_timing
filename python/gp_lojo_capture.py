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
ap.add_argument("--rpt", action="store_true", help="the R_pT distributions of these pileup jets, ITk only and gated by each t0")
ap.add_argument("--bars", action="store_true", help="one stacked bar per t0 instead of the distribution: the share of these jets the "
                "gate keeps whole (>= 90%% of R_pT; the capture), trims, removes (< 10%%), or cannot act on (no timed R_pT track)")
args = ap.parse_args()
R = S.setup()
J = uproot.open(args.tree)["jets"].arrays(["pt", "eta", "hs", "pu", "ntimed", "rpt_itk", "rpt_t0_waves", "rpt_t0_waves_lojet",
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
if args.rpt:
    ROWS = [("rpt_itk", "ITk only", "#000000", 2), ("rpt_t0_waves", "WAVeS t_{0}", S.COL["waves"], 1),
            ("rpt_t0_waves_lojet", "Jet left out", "kP10Green", 1), ("guarded", "Guarded", "kP10Violet", 1),
            ("rpt_t0_truth", "Truth t_{0}", S.COL["truth"], 3)]
    cv = S.keep(R.TCanvas("capr", "", 1600, 700))
    for p, (lab, lo, hi) in enumerate((("30 < p_{T} < 50 GeV", 30, 50), ("p_{T} > 50 GeV", 50, 1e9))):
        m = (J["pt"] > lo) & (J["pt"] < hi) & (ae >= 2.4) & (ae < 4.0)
        hs, pu = m & (J["hs"] == 1), m & (J["pu"] == 1)
        cut = itk_cut(J["rpt_itk"][hs], J["rpt_itk"][pu])
        sel = pu & (J["rpt_itk"] >= cut)
        cv.cd(); S.pad(f"r{p}", 0.5 * p, 0.0, 0.5 * p + 0.5, 1.0, top=0.05, bottom=0.14, left=0.14, right=0.04)
        hists = []
        for col, name, colour, ls in ROWS:
            h = S.keep(R.TH1D(f"r{p}{col}", "", 20, 0.0, 1.0))
            for x in np.clip(J[col][sel], 0.0, 0.9999): h.Fill(x)
            h.Scale(1.0 / max(h.Integral(), 1)); hists.append(h)
        top = max(h.GetMaximum() for h in hists)
        fr = S.keep(R.gPad.DrawFrame(0.0, 0.0, 1.0, 2.4 * top))
        S.axes(fr, "R_{pT}   (overflow in the last bin)", "Fraction of these pileup jets", size=0.045, ytoff=1.45, xtoff=1.3)
        ln = S.keep(R.TLine(cut, 0.0, cut, 1.45 * top)); ln.SetLineStyle(2); ln.SetLineColor(R.kGray + 2); ln.Draw()
        t = S.keep(R.TLatex(cut - 0.01, 1.38 * top, f"ITk cut {cut:.2f}")); t.SetTextFont(42); t.SetTextSize(0.034); t.SetTextAlign(31); t.SetTextColor(R.kGray + 2); t.Draw()
        leg = S.legend(0.47, 0.36, 0.95, 0.64, size=0.036, margin=0.16)
        leg.SetHeader("still above the ITk cut after the gate:")
        print(f"{lab.replace('_{T}', 'T')}: {sel.sum()} pileup jets with ITk-only R_pT >= {cut:.3f}")
        for (col, name, colour, ls), h in zip(ROWS, hists):
            S.style(h, colour, 0, ls, 3); h.Draw("HIST SAME")
            leg.AddEntry(h, f"{name}: {100 * np.mean(J[col][sel] >= cut):.0f}%" if col != "rpt_itk" else "ITk only: 100% (by selection)", "L")
            print(f"  {name.replace('_{0}', '0'):12s} still >= the ITk cut: {100 * np.mean(J[col][sel] >= cut):5.1f}%, median R_pT {np.median(J[col][sel]):.3f}")
        leg.Draw()
        S.header(0.18, 0.88, [f"#sqrt{{s}} = 14 TeV, #LT#mu#GT = 200, {args.label}",
                              f"{lab}, 2.4 < |#eta| < 4.0: the {sel.sum():,} pileup jets",
                              "passing the ITk-only cut at 85% HS efficiency"], size=0.040)
        R.gPad.RedrawAxis()
        S.end_pad()
    S.save(cv, args.out)
    raise SystemExit
if args.bars:
    # categories, left to right: the capture against the axis, so it compares directly between the t0s
    CATS = [("Kept whole (#geq 90%)", "kP10Red", R.kWhite), ("Trimmed", "#9CC3E6", R.kBlack),
            ("Removed (< 10%)", "kP10Blue", R.kWhite), ("No timed track", S.COL["none"], R.kBlack)]
    ROWS = [("rpt_t0_waves", "WAVeS t_{0}"), ("rpt_t0_waves_lojet", "Jet left out"), ("guarded", "Guarded"),
            ("rpt_t0_truth", "Truth t_{0}")]
    def dtext(x, y, txt, size, align, colour=R.kBlack):   # in the frame's coordinates (S.text is NDC)
        t = R.TLatex(); t.SetTextFont(42); t.SetTextSize(size); t.SetTextAlign(align); t.SetTextColor(colour)
        return S.keep(t.DrawLatex(x, y, txt))
    cv = S.keep(R.TCanvas("capb", "", 1600, 680))
    for p, (lab, lo, hi) in enumerate((("30 < p_{T} < 50 GeV", 30, 50), ("p_{T} > 50 GeV", 50, 1e9))):
        m = (J["pt"] > lo) & (J["pt"] < hi) & (ae >= 2.4) & (ae < 4.0)
        hs, pu = m & (J["hs"] == 1), m & (J["pu"] == 1)
        sel = pu & (J["rpt_itk"] >= itk_cut(J["rpt_itk"][hs], J["rpt_itk"][pu]))
        cv.cd(); S.pad(f"b{p}", 0.5 * p, 0.0, 0.5 * p + 0.5, 1.0, top=0.40, bottom=0.14, left=0.21, right=0.04)
        fr = S.keep(R.gPad.DrawFrame(0.0, 0.0, 1.0, len(ROWS)))
        S.axes(fr, "Fraction of these pileup jets", "", size=0.05, xtoff=1.2)
        fr.GetYaxis().SetLabelSize(0); fr.GetYaxis().SetTickLength(0); fr.GetXaxis().SetNdivisions(505)
        untimed = J["ntimed"][sel] == 0
        print(f"{lab.replace('_{T}', 'T')}: {sel.sum()} pileup jets pass the ITk-only cut; no timed R_pT track {100 * untimed.mean():.1f}%")
        for i, (col, name) in enumerate(ROWS):
            f = J[col][sel] / J["rpt_itk"][sel]
            share = [np.mean(~untimed & (f >= 0.9)), np.mean(~untimed & (f >= 0.1) & (f < 0.9)),
                     np.mean(~untimed & (f < 0.1)), np.mean(untimed)]
            y = len(ROWS) - 1 - i
            x = 0.0
            for (cname, ccol, tcol), w in zip(CATS, share):
                b = S.keep(R.TBox(x, y + 0.18, x + w, y + 0.82)); b.SetFillColor(S.colour(ccol)); b.SetLineWidth(0); b.Draw()
                gap = S.keep(R.TLine(x + w, y + 0.18, x + w, y + 0.82)); gap.SetLineColor(R.kWhite); gap.SetLineWidth(3); gap.Draw()
                if w >= 0.07:
                    dtext(x + w / 2, y + 0.5, f"{100 * w:.0f}%", 0.042, 22, tcol)
                x += w
            dtext(-0.02, y + 0.5, name, 0.045, 32)
            print(f"  {name.replace('_{0}', '0'):12s} " + " | ".join(f"{c[0].split(' (')[0]} {100 * w:.1f}%" for c, w in zip(CATS, share)))
        leg = S.legend(0.21, 0.615, 0.96, 0.715, size=0.040, ncol=2, margin=0.12)
        for cname, ccol, tcol in CATS:
            bx = S.keep(R.TBox()); bx.SetFillColor(S.colour(ccol)); bx.SetLineWidth(0)
            leg.AddEntry(bx, cname, "F")
        leg.Draw()
        S.header(0.21, 0.94, [f"#sqrt{{s}} = 14 TeV, #LT#mu#GT = 200, {args.label}",
                              f"{lab}, 2.4 < |#eta| < 4.0: the {sel.sum():,} pileup jets",
                              "passing the ITk-only cut at 85% HS efficiency"], size=0.042)
        R.gPad.RedrawAxis()
        S.end_pad()
    S.save(cv, args.out)
    raise SystemExit
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
