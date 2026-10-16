#!/usr/bin/env python3
"""The jet-time cone: t_jet - t_HS^truth of the signal's truth-matched hard-scatter tagging jets
for two util/vbs_time_veto outputs that differ only in --jet-time-dr (default 0.2 against e.g.
0.4), jets timed with >= --min-trk tracks. The legend gives, per cone, the share of HS legs that
get a usable time, the core sigma (narrow width of a double-Gaussian fit, FIT_DBLGAUS's recipe)
and the share beyond 60 ps.

    PYTHONNOUSERSITE=1 PYTHONPATH=/opt/homebrew/Cellar/root/6.40.04/lib/root \\
    ~/.venv-hgtd/bin/python python/vbs_jet_time_cone.py \\
        condor/vbf/vbf_jvtLoose_vbs_time_veto.root condor/vbf/vbf_jvtLoose_jtdr0p4_vbs_time_veto.root \\
        --out figs/time_veto_k2/jet_time_cone_dt
"""
import argparse
import numpy as np
import uproot
import vbs_time_veto_common as C

ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
ap.add_argument("trees", nargs="+", help="vbs_time_veto outputs (signal), one per cone")
ap.add_argument("--min-trk", type=int, default=2)
ap.add_argument("--out", required=True)
args = ap.parse_args()

ROOT = C.setup_root(args.out)
COLS = ["#9E9E9E", "kP10Blue", "#C0392B", "#E69F00"]
c = ROOT.TCanvas("c", "", 800, 650); C.KEEP.append(c)
c.SetLeftMargin(0.14); c.SetRightMargin(0.05); c.SetLogy(True)
leg = ROOT.TLegend(0.16, 0.735 - 0.06 * len(args.trees), 0.60, 0.735); ROOT.StyleLegend(leg, 0.03); C.KEEP.append(leg)
for i, path in enumerate(args.trees):
    t = uproot.open(path)
    dr = float(np.unique(t["meta"].arrays(["jet_time_dr"], library="np")["jet_time_dr"])[0])
    e = t["events"].arrays(["t0_truth"] + [f"{L}_{k}" for L in "ab" for k in ["hs", "core_n", "core_t"]], library="np")
    dts, nhs, ntm = [], 0, 0
    for L in "ab":
        hs = e[f"{L}_hs"]; tm = hs & (e[f"{L}_core_n"] >= args.min_trk)
        dts.append((e[f"{L}_core_t"] - e["t0_truth"])[tm]); nhs += hs.sum(); ntm += tm.sum()
    dt = np.concatenate(dts)
    h = ROOT.TH1D(f"h{i}", "", 80, -400, 400); C.KEEP.append(h)
    for x in np.clip(dt, -399.9, 399.9): h.Fill(x)
    h.Scale(1.0 / h.Integral())
    h.SetLineColor(C.colour(COLS[i % len(COLS)])); h.SetLineWidth(2)
    if i == 0: h.SetFillColor(C.colour(COLS[0])); h.SetFillStyle(3004)
    if i == 0:
        h.SetMinimum(2e-5); h.SetMaximum(2000.0)
        h.GetXaxis().SetTitle("t_{jet} - t_{HS}^{truth} [ps]   (overflow in the edge bins)")
        h.GetYaxis().SetTitle("Fraction of truth-matched HS tagging jets")
        h.Draw("HIST")
    else:
        h.Draw("HIST SAME")
    leg.AddEntry(h, f"#DeltaR < {dr:.1f}: {100 * ntm / nhs:.0f}% of HS jets timed; core #sigma {C.dgaus_core(dt)[0]:.1f} ps, "
                 f"{100 * np.mean(np.abs(dt) > 60):.1f}% beyond 60 ps", "f" if i == 0 else "l")
leg.Draw()
ROOT.ATLASLabel(0.17, 0.88, "")
tx = ROOT.TLatex(); tx.SetNDC(); tx.SetTextFont(42); tx.SetTextSize(0.04); C.KEEP.append(tx)
tx.DrawLatex(0.31, 0.88, "Simulation Internal")
t2 = ROOT.TLatex(); t2.SetNDC(); t2.SetTextFont(42); t2.SetTextSize(0.028); t2.SetTextColor(ROOT.kGray + 2); C.KEEP.append(t2)
t2.DrawLatex(0.17, 0.835, "VBF H#rightarrowinv., truth-matched hard-scatter tagging jets, JVT + fJVT loose")
t2.DrawLatex(0.17, 0.80, f"jet time from #geq {args.min_trk} timed ghost tracks within #DeltaR of the axis, 2#sigma time clustering")
t2.DrawLatex(0.17, 0.765, "core #sigma: narrow width of a double-Gaussian fit")
c.Print(args.out + ".pdf"); c.Print(args.out + ".png")
print("wrote", args.out + ".pdf / .png")
