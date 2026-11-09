#!/usr/bin/env python3
"""The t0 core fraction by VBS pair type, Athena / TRKPTZ / WAVeS, for the Global Performance talk
(slides/slides_gp_timing.tex): WAVeS gains when the forward VBS jets are hard scatter and is pulled away when the only
forward jet is pileup.

Events: section 3's selection before pileup tagging (tdr_fwdpair --pair=one, --jvt=none: two VBS jets, >= 1 in HGTD),
truth labels (TDR), generator weights. Core fraction = |t0 - t_HS| < 60 ps over all events (no t0 = failure).

  PYTHONNOUSERSITE=1 ~/.venv-hgtd/bin/python gp_t0_by_pair.py --trees ../condor/vbf/vbf_onefwd_tdr_fwdpair.root \\
      ../condor/zjets/zjets_onefwd_tdr_fwdpair.root --out ../figs/tdr/gp/gp_t0_by_pair
"""
import argparse

import numpy as np
import uproot

import tdr_root_style as S

ap = argparse.ArgumentParser()
ap.add_argument("--trees", nargs=2, required=True, metavar=("VBF", "ZJETS"))
ap.add_argument("--out", required=True)
args = ap.parse_args()
T0S = [("hgtd", "HGTD (Athena) t_{0}", S.COL["hgtd"]), ("trkptz", "TRKPTZ t_{0}", S.COL["trkptz"]), ("waves", "WAVeS t_{0}", S.COL["waves"])]
LABELS = ["both HS", "fwd HS, cen PU", "fwd PU, cen HS", "both PU"]

R = S.setup()
cv = S.keep(R.TCanvas("bypair", "", 1400, 700))
for c, (path, lab) in enumerate(zip(args.trees, ("VBF H#rightarrowinv.", "Z+jets"))):
    E = uproot.open(path)["events"].arrays(library="np")
    w = E["weight"].astype(float)
    fwd = lambda k: np.abs(E[f"j{k}_eta"]) > 2.4
    hs = lambda k: E[f"j{k}_hs"] == 1
    pu = lambda k: E[f"j{k}_pu"] == 1
    masks = [hs(1) & hs(2),
             (E["cat"] == 1) & ((fwd(1) & hs(1) & pu(2)) | (fwd(2) & hs(2) & pu(1))),
             (E["cat"] == 1) & ((fwd(1) & pu(1) & hs(2)) | (fwd(2) & pu(2) & hs(1))),
             pu(1) & pu(2)]
    cv.cd(); S.pad(f"bp_{c}", 0.5 * c, 0.0, 0.5 * c + 0.5, 1.0, top=0.06, bottom=0.14, left=0.13, right=0.04)
    fr = S.keep(R.TH1D(f"fr_{c}", "", len(masks), 0, len(masks)))
    fr.SetMinimum(0.0); fr.SetMaximum(1.75)
    for i, l in enumerate(LABELS): fr.GetXaxis().SetBinLabel(i + 1, l)
    fr.GetXaxis().SetLabelSize(0.044); fr.GetXaxis().LabelsOption("h"); fr.GetYaxis().SetTitle("Core fraction (|t_{0} #minus t_{HS}| < 60 ps)")
    fr.GetYaxis().SetTitleSize(0.045); fr.GetYaxis().SetLabelSize(0.04); fr.GetYaxis().SetTitleOffset(1.3)
    fr.Draw("AXIS")
    leg = S.legend(0.17, 0.665, 0.93, 0.715, size=0.034, ncol=3, margin=0.25)
    print(f"{lab.replace('#rightarrow', '->')}:")
    for k, (t0, name, col) in enumerate(T0S):
        h = S.keep(R.TH1D(f"h_{c}_{t0}", "", len(masks), 0, len(masks)))
        ok = (E[f"ok_{t0}"] == 1) & (np.abs(E[f"t0_{t0}"] - E["t0_truth"]) < 60)
        for i, m in enumerate(masks):
            f = w[m & ok].sum() / w[m].sum()
            h.SetBinContent(i + 1, f)
        h.SetFillColor(S.colour(col)); h.SetLineColor(S.colour(col)); h.SetMarkerSize(0); h.SetBarWidth(0.26); h.SetBarOffset(0.11 + 0.26 * k)
        h.Draw("BAR SAME")
        leg.AddEntry(h, name, "F")
        print(f"  {t0:7s} " + "  ".join(f"{LABELS[i]} {100 * h.GetBinContent(i + 1):.1f}" for i in range(len(masks))))
    for i, m in enumerate(masks):
        sh = 100 * w[m].sum() / w.sum()
        t = S.keep(R.TLatex(i + 0.5, 1.02, f"{sh:.0f}% of events" if sh >= 1 else f"{sh:.1f}% of events"))
        t.SetTextAlign(21); t.SetTextSize(0.032); t.SetTextColor(R.kGray + 2); t.Draw()
    leg.Draw()
    S.header(0.17, 0.88, ["#sqrt{s} = 14 TeV, #LT#mu#GT = 200",
                          f"{lab}: two VBS jets, #geq 1 in HGTD; HS / PU: truth labels"], size=0.034)
    R.gPad.RedrawAxis()
S.save(cv, args.out)
