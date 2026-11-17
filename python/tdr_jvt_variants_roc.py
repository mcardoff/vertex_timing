#!/usr/bin/env python3
"""Per-jet performance of JVT's R_pT on forward jets, plain and time-gated (2026-11-09 mini deck,
slides/slides_jvt_variants.tex): hard-scatter vs pileup efficiency as the threshold is scanned.

Jets: util/tdr_pujet trees with the JVT columns (commit c182466): jvt_rpt (the jet's ghost tracks the vertex fit gave
PV0, over its pT) and jvt_rpt_t0_<t0> (the same with the section-2 time gate on those tracks: 2 sigma, cleaned times,
untimed tracks count). Forward 2.5 <= |eta| < 4.0 and 30 <= pT < 60 GeV: where the baseline (jvt4pt60) applies JVT
beyond its Run 2 range. TDR truth labels; generator weights. Markers: the loose threshold 0.0167 on every curve, and the
plain R_pT at the tighter thresholds the analysis reruns use (0.03, 0.05, 0.10).

  PYTHONNOUSERSITE=1 ~/.venv-hgtd/bin/python tdr_jvt_variants_roc.py \\
      --jets ../condor/vbf/vbf_tdr_pujet_jvt.root ../condor/zjets/zjets_tdr_pujet_jvt.root --out ../figs/tdr/jvt_variants/jvt_roc
"""
import argparse

import numpy as np
import uproot

import tdr_root_style as S

ap = argparse.ArgumentParser()
ap.add_argument("--jets", nargs=2, required=True, metavar=("VBF", "ZJETS"))
ap.add_argument("--out", required=True)
ap.add_argument("--pt", nargs=2, type=float, default=(30.0, 60.0))
args = ap.parse_args()
LO, HI = args.pt
THR = np.concatenate([[0.0], np.geomspace(1e-3, 0.6, 400)])
VARS = [("jvt_rpt", "JVT R_{pT} (vertex fit)", "#000000", 1),
        ("jvt_rpt_t0_waves", "time-gated, WAVeS t_{0}", S.COL["waves"], 1),
        ("jvt_rpt_t0_truth", "time-gated, truth t_{0}", S.COL["truth"], 2)]
MARK = [0.03, 0.05, 0.10]


def eff(x, w, t):
    o = np.argsort(x); xs, cw = x[o], np.concatenate([[0.0], np.cumsum(w[o])])
    return 1.0 - cw[np.searchsorted(xs, t, side="left")] / w.sum()   # weight with x >= t


R = S.setup()
cv = S.keep(R.TCanvas("jvtroc", "", 1400, 700))
for c, (path, lab) in enumerate(zip(args.jets, ("VBF H#rightarrowinv.", "Z+jets"))):
    J = uproot.open(path)["jets"].arrays(["weight", "pt", "eta", "hs", "pu"] + [v for v, *_ in VARS], library="np")
    ae, w = np.abs(J["eta"]), J["weight"].astype(float)
    sel = (ae >= 2.5) & (ae < 4.0) & (J["pt"] >= LO) & (J["pt"] < HI)
    hs, pu = sel & (J["hs"] == 1), sel & (J["pu"] == 1)
    print(f"\n{lab.replace('#rightarrow', '->')}: forward jets {LO:.0f}-{HI:.0f} GeV: {hs.sum():,} HS, {pu.sum():,} PU")
    cv.cd(); S.pad(f"roc_{c}", 0.5 * c, 0.34, 0.5 * c + 0.5, 1.0, top=0.06, bottom=0.02, left=0.15, right=0.04)
    R.gPad.SetLogy()
    fr = S.keep(R.gPad.DrawFrame(0.55, 3e-3, 1.0, 1.0))
    S.axes(fr, "", "Pileup jet efficiency", size=0.055, ytoff=1.25)
    fr.GetXaxis().SetLabelSize(0)
    leg = S.legend(0.50, 0.06, 0.95, 0.36, size=0.042)   # lower right: the curves are above 0.03 there
    curves = {}
    for v, name, col, ls in VARS:
        eh, ep = eff(J[v][hs], w[hs], THR), eff(J[v][pu], w[pu], THR)
        curves[v] = (eh, ep)
        g = S.style(S.graph(eh, ep), col, 0, ls, 3); g.Draw("L SAME")
        leg.AddEntry(g, name, "L")
        m = S.style(S.graph([eff(J[v][hs], w[hs], 0.0167)], [eff(J[v][pu], w[pu], 0.0167)]), col, 20, ls, 2, 1.6); m.Draw("P SAME")
        print(f"  {name.replace('#', '').replace('_{pT}', 'pT').replace('_{0}', '0'):28s} at 0.0167: HS {eff(J[v][hs], w[hs], 0.0167):.3f}  PU {eff(J[v][pu], w[pu], 0.0167):.4f}")
    for t in MARK:
        eh, ep = eff(J["jvt_rpt"][hs], w[hs], t), eff(J["jvt_rpt"][pu], w[pu], t)
        m = S.style(S.graph([eh], [ep]), "#000000", 25, 1, 2, 1.5); m.Draw("P SAME")
        S.keep(R.TLatex(eh - 0.008, ep * 1.25, f"{t:g}")).SetTextSize(0.03)
        print(f"  plain JVT R_pT >= {t:<5g}            HS {eh:.3f}  PU {ep:.4f}")
    # the plain threshold with the same HS efficiency as the WAVeS-gated loose cut, and its PU efficiency
    target = eff(J["jvt_rpt_t0_waves"][hs], w[hs], 0.0167)
    eh, ep = curves["jvt_rpt"]
    i = int(np.argmin(np.abs(eh - target)))
    print(f"  at the WAVeS-gated loose cut's HS efficiency ({target:.3f}): plain threshold {THR[i]:.4f}, PU {ep[i]:.4f} "
          f"(gated: {eff(J['jvt_rpt_t0_waves'][pu], w[pu], 0.0167):.4f})")
    dummy = S.keep(R.TGraph()); S.style(dummy, "#000000", 20, 1, 2, 1.4)
    leg.AddEntry(dummy, "R_{pT} #geq 0.0167 (loose)", "P")
    d2 = S.keep(R.TGraph()); S.style(d2, "#000000", 25, 1, 2, 1.4)
    leg.AddEntry(d2, "plain R_{pT} #geq 0.03, 0.05, 0.10", "P")
    leg.Draw()
    S.header(0.19, 0.88, ["#sqrt{s} = 14 TeV, HL-LHC, #LT#mu#GT = 200",
                          f"{lab}, jets 2.5 #leq |#eta| < 4.0, {LO:.0f} #leq p_{{T}} < {HI:.0f} GeV"], size=0.044)
    R.gPad.RedrawAxis()
    # lower panel: each variant's pileup efficiency over plain JVT's at the same hard-scatter efficiency
    cv.cd(); S.pad(f"rat_{c}", 0.5 * c, 0.0, 0.5 * c + 0.5, 0.34, top=0.03, bottom=0.32, left=0.15, right=0.04)
    fr2 = S.keep(R.gPad.DrawFrame(0.55, 0.3, 1.0, 1.5))
    S.axes(fr2, "Hard-scatter jet efficiency", "PU eff. / plain JVT's", size=0.105, ytoff=0.62, xtoff=1.05, ndiv_y=504)
    one = S.keep(R.TLine(0.55, 1.0, 1.0, 1.0)); one.SetLineStyle(2); one.SetLineColor(R.kGray + 2); one.Draw()
    eh0, ep0 = curves["jvt_rpt"]
    o0 = np.argsort(eh0)
    xs = np.linspace(0.56, eff(J["jvt_rpt"][hs], w[hs], 0.0167), 120)   # up to the loose point: below 0.0167 the
                                                                        # plain curve has almost no working points
    base = np.interp(xs, eh0[o0], ep0[o0])
    for v, name, col, ls in VARS[1:]:
        eh, ep = curves[v]; o = np.argsort(eh)
        g = S.style(S.graph(xs, np.interp(xs, eh[o], ep[o]) / base), col, 0, ls, 3); g.Draw("L SAME")
    R.gPad.RedrawAxis()
S.save(cv, args.out)
