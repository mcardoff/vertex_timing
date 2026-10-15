#!/usr/bin/env python3
"""HGTD acceptance of the VBS tagging pair, defined by the legs' TRACKS, not their
axes: a tagging jet is in acceptance when >= --min-trk of its ghost tracks pass
the track selection inside HGTD (vbs_time_veto's hgtd_n: 1 < pT < 30 GeV,
quality, 2.38 < |eta| < 4.0, timed or not). A jet whose axis is past |eta| = 4.0
(or just inside 2.38) can still put tracks there.

Reads util/vbs_time_veto outputs (the region-plot selection, not --sr) and writes
  <out>_vs_eta.(pdf|png)  share of tagging jets with >= --min-trk (1) and >= --good-trk
                          (4) TIMED tracks in HGTD, vs the jet axis |eta|
  <out>_summary.md        weighted share of signal and background, after the
                          recoil cuts, with 2 / 1 / 0 legs in acceptance, by
                          tracks and by axis (2.38 < |eta_axis| < 4.0)

    PYTHONNOUSERSITE=1 PYTHONPATH=/opt/homebrew/Cellar/root/6.40.04/lib/root \\
    ~/.venv-hgtd/bin/python python/vbs_hgtd_acceptance.py \\
        --sig condor/vbf/vbf_jvtLoose_vbs_time_veto.root \\
        --bkg condor/zjets/zjets_jvtLoose_vbs_time_veto.root \\
        --met-min 100 --zpt-min 100 --out figs/time_veto/jvtLoose_recoil100_hgtdacc
"""
import argparse, os
import numpy as np
import uproot
import vbs_time_veto_common as C

ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
ap.add_argument("--sig", required=True)
ap.add_argument("--bkg")
ap.add_argument("--met-min", type=float, default=100.0)
ap.add_argument("--zpt-min", type=float, default=100.0)
ap.add_argument("--min-trk", type=int, default=1, help="the loose case: >= this many TIMED HGTD tracks")
ap.add_argument("--good-trk", type=int, default=4, help="the 'good time' case: >= this many TIMED HGTD tracks")
ap.add_argument("--time-trk", type=int, default=0, help="if > 0, the second curve is instead 'a jet time built from "
                ">= this many tracks' (core_n: the dR < 0.2 time cluster) -- the follow-up's jet-time requirement")
ap.add_argument("--tag", default="JVT + fJVT loose")
ap.add_argument("--out", required=True)
args = ap.parse_args()

COLS = ["weight", "met_truth", "z_pt"] + [f"{L}_{k}" for L in "ab" for k in ["eta", "hgtd_n", "hgtd_ntimed", "hs", "pu", "core_n"]]
ETA_AXIS = (2.38, 4.0)
EDGES = np.round(np.arange(1.6, 4.61, 0.1), 2)


def load(path, col, cut):
    e = uproot.open(path)["events"].arrays(COLS, library="np")
    if cut > 0:
        if np.all(e[col] < 0): raise SystemExit(f"{path}: {col} unavailable")
        k = e[col] > cut
        e = {c: v[k] for c, v in e.items()}
    return e


S = {"sig": load(args.sig, "met_truth", args.met_min)}
if args.bkg: S["bkg"] = load(args.bkg, "z_pt", args.zpt_min)
NAME = {"sig": "VBF H#rightarrowinv.", "bkg": "Z(#rightarrowll)+jets"}
MDNAME = {"sig": "VBF H->inv", "bkg": "Z(->ll)+jets"}

# ── category table ───────────────────────────────────────────────────────────
md = [f"# HGTD acceptance of the tagging pair: {os.path.basename(args.out)}", "",
      f"A leg is in acceptance by TRACKS when enough of its ghost tracks pass the track selection in "
      f"2.38 < |eta| < 4.0 (rows: >= {args.min_trk} / >= {args.good_trk} with a valid HGTD time, or >= {args.min_trk} timed or "
      f"not); by AXIS when {ETA_AXIS[0]} < |eta_jet| < {ETA_AXIS[1]}. "
      f"Recoil: VBF truth E_T^miss > {args.met_min:g} GeV" + (f", Z+jets pT(ll) > {args.zpt_min:g} GeV" if args.bkg else "") + ". "
      "Shares are generator-weighted; MC event counts in brackets.", "",
      "| sample | definition | both legs | one leg | neither |", "|---|---|---|---|---|"]
for r, e in S.items():
    w = e["weight"].astype(float)
    for lab, inA, inB in [(f">= {args.min_trk} timed track(s)", e["a_hgtd_ntimed"] >= args.min_trk, e["b_hgtd_ntimed"] >= args.min_trk),
                          (f">= {args.good_trk} timed tracks", e["a_hgtd_ntimed"] >= args.good_trk, e["b_hgtd_ntimed"] >= args.good_trk),
                          (f">= {args.min_trk} track(s), timed or not", e["a_hgtd_n"] >= args.min_trk, e["b_hgtd_n"] >= args.min_trk),
                          ("axis", (np.abs(e["a_eta"]) > ETA_AXIS[0]) & (np.abs(e["a_eta"]) < ETA_AXIS[1]),
                                   (np.abs(e["b_eta"]) > ETA_AXIS[0]) & (np.abs(e["b_eta"]) < ETA_AXIS[1]))]:
        n = inA.astype(int) + inB.astype(int)
        cells = [f"{100 * w[n == k].sum() / w.sum():.1f}% ({int((n == k).sum()):,})" for k in (2, 1, 0)]
        md.append(f"| {MDNAME[r]} | {lab} | " + " | ".join(cells) + " |")
md.append("")
# the pair's truth composition inside the both-by-tracks category
md += [f"**Inside 'both legs with >= {args.min_trk} timed track(s)': the pair's truth labels** (weighted share)", "",
       "| sample | HS+HS | HS+PU | PU+PU | other |", "|---|---|---|---|---|"]
for r, e in S.items():
    w = e["weight"].astype(float)
    m = (e["a_hgtd_ntimed"] >= args.min_trk) & (e["b_hgtd_ntimed"] >= args.min_trk)
    hA, hB, pA, pB = e["a_hs"], e["b_hs"], e["a_pu"], e["b_pu"]
    comp = [hA & hB, (hA & pB) | (pA & hB), pA & pB]
    sh = [w[m & c].sum() / w[m].sum() for c in comp]
    md.append(f"| {MDNAME[r]} | " + " | ".join(f"{100 * x:.1f}%" for x in sh) + f" | {100 * (1 - sum(sh)):.1f}% |")
md.append("")
os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)
open(args.out + "_summary.md", "w").write("\n".join(md) + "\n")
print("\n".join(md))

# ── share of tagging jets in acceptance vs axis |eta| ─────────────────────────
ROOT = C.setup_root(args.out)
c = ROOT.TCanvas("c_acc", "", 800, 650); C.KEEP.append(c)
c.SetLeftMargin(0.13); c.SetRightMargin(0.04); c.SetTopMargin(0.05)
fr = ROOT.TH1D("fr_acc", "", 1, EDGES[0], EDGES[-1]); C.KEEP.append(fr)
fr.SetMinimum(0.0); fr.SetMaximum(2.0)
fr.GetXaxis().SetTitle("tagging-jet axis |#eta|")
fr.GetYaxis().SetTitle("Fraction of tagging jets in HGTD acceptance")
fr.Draw("AXIS")
box = ROOT.TBox(ETA_AXIS[0], 0.0, ETA_AXIS[1], 1.0); box.SetFillColorAlpha(ROOT.kGray, 0.25); box.Draw()
one = ROOT.TLine(EDGES[0], 1.0, EDGES[-1], 1.0); one.SetLineStyle(2); one.SetLineColor(ROOT.kGray + 1); one.Draw()
C.KEEP.extend([box, one])
leg = ROOT.TLegend(0.17, 0.62, 0.93, 0.745); ROOT.StyleLegend(leg, 0.030); leg.SetNColumns(2); C.KEEP.append(leg)
STY = {("sig", "n"): ("#000000", 20, 1), ("sig", "t"): ("#000000", 24, 2),
       ("bkg", "n"): ("kP10Blue", 21, 1), ("bkg", "t"): ("kP10Blue", 25, 2)}
for r, e in S.items():
    eta = np.concatenate([np.abs(e["a_eta"]), np.abs(e["b_eta"])])
    w = np.concatenate([e["weight"], e["weight"]]).astype(float)
    for key, k, lab in [("n", args.min_trk, "#geq 1 timed track" if args.min_trk == 1 else f"#geq {args.min_trk} timed tracks"),
                        ("t", args.time_trk or args.good_trk,
                         f"jet time from #geq {args.time_trk} tracks" if args.time_trk else f"#geq {args.good_trk} timed tracks")]:
        col = "core_n" if (key == "t" and args.time_trk) else "hgtd_ntimed"
        has = np.concatenate([e[f"a_{col}"], e[f"b_{col}"]]) >= k
        g = ROOT.TGraphAsymmErrors(); C.KEEP.append(g)
        idx = np.digitize(eta, EDGES) - 1
        for b in range(len(EDGES) - 1):
            m = idx == b
            if m.sum() < 10: continue
            v, _, lo, hi = C.eff(w[m], has[m])
            i = g.GetN(); g.SetPoint(i, 0.5 * (EDGES[b] + EDGES[b + 1]), v); g.SetPointError(i, 0, 0, lo, hi)
        cc, mk, ls = STY[(r, key)]
        g.SetMarkerColor(C.colour(cc)); g.SetLineColor(C.colour(cc)); g.SetMarkerStyle(mk); g.SetLineStyle(ls)
        g.SetMarkerSize(1.0); g.SetLineWidth(2)
        g.Draw("PL SAME")
        leg.AddEntry(g, f"{NAME[r]}: {lab}", "pl")
leg.Draw()
ROOT.ATLASLabel(0.17, 0.89, "")
t = ROOT.TLatex(); t.SetNDC(); t.SetTextFont(42); t.SetTextSize(0.04); C.KEEP.append(t)
t.DrawLatex(0.30, 0.89, "Simulation Internal")
t2 = ROOT.TLatex(); t2.SetNDC(); t2.SetTextFont(42); t2.SetTextSize(0.028); t2.SetTextColor(ROOT.kGray + 2); C.KEEP.append(t2)
t2.DrawLatex(0.17, 0.845, f"#sqrt{{s}} = 14 TeV, HL-LHC; VBS pair, m_{{jj}} > 200 GeV, {args.tag}, recoil > {args.met_min:g} GeV")
t2.DrawLatex(0.17, 0.81, "in acceptance: the jet's timed ghost tracks (1 < p_{T} < 30 GeV) in 2.38 < |#eta_{trk}| < 4.0")
t2.DrawLatex(0.17, 0.775, "shaded: the jet axis inside HGTD")
c.RedrawAxis()
c.Print(args.out + "_vs_eta.pdf"); c.Print(args.out + "_vs_eta.png")
print("wrote", args.out + "_vs_eta.pdf / .png and _summary.md")
