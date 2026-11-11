#!/usr/bin/env python3
"""What fixing each t0 failure mode is worth, for the Global Performance talk (slides/slides_gp_timing.tex).

Reads clustering_plot's failure_fixes.json (one per sample): the WAVeS t0 with its failure modes repaired one after
another with truth, in the failure ladder's order (src/failure_fixes.h, MyUtl::Fix):
  re-timing       the full-cluster time when the in-jet subset is the less pure
  HS track times  the hard-scatter tracks >= 3 sigma from t_HS (misassigned HGTD hits) dropped, re-clustered
  selection       the hard scatter's own cluster chosen (most HS pT)
  misclustering   that cluster timed from its HS tracks alone
  HS track eff.   every untimed HS track in HGTD given a time, Gaus(t_HS, 30 ps)
  HS time res.    every HS track at t_HS (1 ps): what is left has no forward HS track at all
Upper panel: core fraction (|t0 - t_HS| < 60 ps, all events), a waterfall of the cumulative chain; each step is labelled
with its gain in the chain and, in brackets, the gain of that fix alone on WAVeS. Lower panel: the resolution
(FIT_PEAKTAIL's peak width) along the chain; dashed: TRKPTZ and Athena.

  PYTHONNOUSERSITE=1 ~/.venv-hgtd/bin/python gp_failure_fixes.py --json ../condor/vbf/fixes/vbf_cleant0_assoc2p5_failure_fixes.json \\
      ../condor/zjets/fixes/zjets_cleant0_assoc2p5_failure_fixes.json --out ../figs/tdr/gp/gp_failure_fixes
"""
import argparse
import json

import tdr_root_style as S

ap = argparse.ArgumentParser()
ap.add_argument("--json", nargs="+", required=True)
ap.add_argument("--labels", nargs="+", default=None)
ap.add_argument("--out", required=True)
args = ap.parse_args()
LABELS = args.labels or [{"vbf": "VBF H#rightarrowinv.", "zjets": "Z+jets", "local": "VBF H#rightarrowinv. (local)"}
                         .get(json.load(open(p))["sample"], json.load(open(p))["sample"]) for p in args.json]
# the chain: (short name of the cumulative row, label, colour (the ladder's), short name of the fix alone)
STEPS = [("WFIX_RETIME", "#splitline{re-}{timing}", "#E69F00", "WFIX_RETIME"),
         ("WFIX_C2", "#splitline{HS track}{times}", "#5A3E2B", "WFIX_MISTIMED"),
         ("WFIX_C3", "selection", "#3F90DA", "WFIX_SELECT"),
         ("WFIX_C4", "#splitline{mis-}{clustering}", "#B3231C", "WFIX_MISCL"),
         ("WFIX_C5", "#splitline{HS track}{efficiency}", "#8E5DB5", "WFIX_TEFF"),
         ("WFIX_C6", "#splitline{HS time}{resolution}", "#C98B5E", None)]
plain = lambda n: n.replace("#splitline{", "").replace("}{", " ").replace("}", "").replace("_{0", "0").replace("#", "")
FLOOR_PS = 3.0   # FIT_PEAKTAIL's lower limit on the peak width

R = S.setup()
ncol = len(args.json)
cv = S.keep(R.TCanvas("fixes", "", 850 * ncol, 860))
for c, (path, lab) in enumerate(zip(args.json, LABELS)):
    rows = {r["short"]: r for r in json.load(open(path))["rows"]}
    core = lambda k: 100 * rows[k]["core"]
    xs = list(range(len(STEPS) + 1))
    names = ["#splitline{WAVeS}{t_{0}}"] + [s[1] for s in STEPS]
    levels = [core("WAVES")] + [core(s[0]) for s in STEPS]
    lo = 5 * int(min(core("HGTD"), levels[0]) // 5) - 5
    print(f"{lab.replace('#rightarrow', '->')}: WAVeS {levels[0]:.2f}  TRKPTZ {core('TRKPTZ'):.2f}  Athena {core('HGTD'):.2f}")
    # --- upper: core-fraction waterfall
    cv.cd(); S.pad(f"up_{c}", c / ncol, 0.38, (c + 1) / ncol, 1.0, top=0.25, bottom=0.02, left=0.13, right=0.04)
    fr = S.keep(R.TH1D(f"fu_{c}", "", len(xs), -0.5, len(xs) - 0.5))
    fr.SetMinimum(lo); fr.SetMaximum(100.5)
    fr.GetYaxis().SetTitle("Core fraction [%]")
    fr.GetYaxis().SetTitleSize(0.055); fr.GetYaxis().SetLabelSize(0.050); fr.GetYaxis().SetTitleOffset(1.05)
    fr.GetXaxis().SetLabelSize(0); fr.Draw("AXIS")
    base = S.keep(R.TBox(-0.35, lo, 0.35, levels[0])); base.SetFillColor(S.colour(S.COL["waves"])); base.SetLineColor(R.kBlack); base.Draw("L")
    base2 = S.keep(R.TBox(-0.35, lo, 0.35, levels[0])); base2.SetFillStyle(0); base2.SetLineColor(R.kBlack); base2.Draw()
    for i, (key, name, col, alone) in enumerate(STEPS, start=1):
        y0, y1 = levels[i - 1], levels[i]
        b = S.keep(R.TBox(i - 0.35, min(y0, y1), i + 0.35, max(y0, y1, min(y0, y1) + 0.003 * (100 - lo))))
        b.SetFillColor(S.colour(col)); b.SetLineColor(S.colour(col)); b.Draw()
        ln = S.keep(R.TLine(i - 1 + 0.35, y0, i - 0.35, y0)); ln.SetLineStyle(3); ln.SetLineColor(R.kGray + 2); ln.Draw()
        # under the box: the chain only rises, so the space below each step is always empty
        for k, txt in enumerate([f"+{y1 - y0:.1f}"] + ([f"({core(alone) - levels[0]:+.1f})"] if alone else [])):
            t = S.keep(R.TLatex(i, min(y0, y1) - (0.025 + 0.065 * k) * (100 - lo), txt))
            t.SetTextAlign(23); t.SetTextSize(0.042); t.SetTextColor(R.kBlack if k == 0 else R.kGray + 2); t.Draw()
        name = plain(name)
        print(f"   + {name:15s} {y1:6.2f} ({y1 - y0:+.2f} in the chain" + (f", {core(alone) - levels[0]:+.2f} alone)" if alone else ")"))
    tw = S.keep(R.TLatex(0, levels[0] + 0.02 * (100 - lo), f"{levels[0]:.1f}")); tw.SetTextAlign(21); tw.SetTextSize(0.040); tw.Draw()
    for k, nm, style in (("TRKPTZ", "TRKPTZ", 2), ("HGTD", "Athena", 3)):
        l = S.keep(R.TLine(-0.5, core(k), len(xs) - 0.5, core(k))); l.SetLineStyle(style); l.SetLineColor(R.kGray + 2); l.Draw()
        tt = S.keep(R.TLatex(len(xs) - 0.55, core(k) - 0.045 * (100 - lo), f"{nm} {core(k):.1f}"))
        tt.SetTextAlign(31); tt.SetTextSize(0.036); tt.SetTextColor(R.kGray + 2); tt.Draw()
    hundred = S.keep(R.TLine(-0.5, 100, len(xs) - 0.5, 100)); hundred.SetLineColor(R.kGray); hundred.Draw()
    S.header(0.13, 0.955, ["#sqrt{s} = 14 TeV, #LT#mu#GT = 200;  |t_{0} #minus t_{HS}| < 60 ps, all events",
                           f"{lab}: WAVeS t_{{0}}, its failure modes fixed with truth in turn",
                           "under each step: its gain in the chain;  (grey): the fix alone on WAVeS"], size=0.042)
    R.gPad.RedrawAxis()
    # --- lower: resolution along the chain
    cv.cd(); S.pad(f"dn_{c}", c / ncol, 0.0, (c + 1) / ncol, 0.38, top=0.03, bottom=0.30, left=0.13, right=0.04)
    sig = [rows["WAVES"]["sigma_peak"]] + [rows[s[0]]["sigma_peak"] for s in STEPS]
    err = [rows["WAVES"]["sigma_peak_err"]] + [rows[s[0]]["sigma_peak_err"] for s in STEPS]
    ref = [rows[k]["sigma_peak"] for k in ("TRKPTZ", "HGTD")]
    fd = S.keep(R.TH1D(f"fd_{c}", "", len(xs), -0.5, len(xs) - 0.5))
    fd.SetMinimum(0.0); fd.SetMaximum(1.25 * max(sig + ref))
    fd.GetXaxis().SetLabelSize(0); fd.GetXaxis().SetNdivisions(0)
    fd.GetYaxis().SetTitle("#sigma_{peak} [ps]"); fd.GetYaxis().SetTitleSize(0.085); fd.GetYaxis().SetLabelSize(0.075)
    fd.GetYaxis().SetTitleOffset(0.72); fd.GetYaxis().SetNdivisions(504); fd.Draw("AXIS")
    for i, n in enumerate(names):
        tl = S.keep(R.TLatex(i, -0.09 * fd.GetMaximum(), n)); tl.SetTextAlign(23); tl.SetTextSize(0.068); tl.Draw()
    g = S.keep(R.TGraphErrors(len(xs)))
    for i, (s_, e_) in enumerate(zip(sig, err)): g.SetPoint(i, i, s_); g.SetPointError(i, 0, max(e_, 0))
    g.SetMarkerStyle(20); g.SetMarkerSize(1.3); g.SetLineWidth(2); g.Draw("LP SAME")
    for k, style in (("TRKPTZ", 2), ("HGTD", 3)):
        l = S.keep(R.TLine(-0.5, rows[k]["sigma_peak"], len(xs) - 0.5, rows[k]["sigma_peak"])); l.SetLineStyle(style); l.SetLineColor(R.kGray + 2); l.Draw()
        nm = {"TRKPTZ": "TRKPTZ", "HGTD": "Athena"}[k]
        tt = S.keep(R.TLatex(-0.4, rows[k]["sigma_peak"] + 0.03 * fd.GetMaximum(), nm)); tt.SetTextAlign(11); tt.SetTextSize(0.062); tt.SetTextColor(R.kGray + 2); tt.Draw()
    if sig[-1] <= FLOOR_PS + 0.01:
        t = S.keep(R.TLatex(len(xs) - 1.15, sig[-1], "fit floor")); t.SetTextAlign(32); t.SetTextSize(0.062); t.Draw()
    print("   sigma_peak: " + "  ".join(f"{plain(n)} {s_:.1f}" for n, s_ in zip(names, sig))
          + f" | TRKPTZ {ref[0]:.1f} Athena {ref[1]:.1f}")
    R.gPad.RedrawAxis()
S.save(cv, args.out)
