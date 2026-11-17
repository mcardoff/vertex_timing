#!/usr/bin/env python3
"""Seven-band VBS pair composition vs m_jj, ROOT/ATLAS style, from a
vbs_region_diag output -- untagged or JVT/fJVT-tagged alike.

Every selected event's chosen pair (wide_ block: forward = |eta| >= 2.4, no
upper edge) is put in exactly one band by the truth identity and region of
its two legs, binned in m_jj (overflow folded into the top bin), and every
column is normalised to 1.

    PYTHONNOUSERSITE=1 ~/.venv-hgtd/bin/python python/vbs_region_stack.py \
        condor/zjets/zjets_jvtLoose_vbs_region_diag.root \
        --label "#sqrt{s} = 14 TeV, HL-LHC, Z+jets" \
        --tag "JVT + fJVT loose applied before pairing" \
        --out condor/zjets/zjets_jvtLoose_vbs_region_stack
"""
import argparse, os
import numpy as np
import uproot
import ROOT

SPLIT     = 2.4
DETA_MIN  = 2.5
MJJ_EDGES = [500, 750, 1000, 1250, 1500, 1750, 2000, 2500, 3000, 4000]

# band index, label, fill colour (hex), in stack order bottom -> top
BANDS = [
    ("r1",       "R1: fwd HS + fwd PU",       "#3F8FFF"),
    ("r2",       "R2: fwd PU + central HS",   "#B3231C"),
    ("hspu_oth", "HS + PU, other #eta config", "#F5A02A"),
    ("hshs",     "HS + HS (genuine pair)",    "#7BCB6B"),
    ("pupu_ff",  "PU + PU, both forward",     "#636363"),
    ("pupu_oth", "PU + PU, other #eta config", "#CFCFCF"),
    ("neither",  "A leg carries neither label", "#C79BF5"),
]

ap = argparse.ArgumentParser()
ap.add_argument("file")
ap.add_argument("--label", default="#sqrt{s} = 14 TeV, HL-LHC, VBF H#rightarrowinv.")
ap.add_argument("--tag", default="", help="third label line suffix, e.g. the working point")
ap.add_argument("--deta", type=float, default=DETA_MIN)
ap.add_argument("--zpt-min", type=float, default=0.0,
                help="require the dilepton pT above this (GeV): the Z->ll stand-in for the SR's MET cut; "
                     "needs the z_pt column (Z+jets outputs from 2026-09-18 on)")
ap.add_argument("--opphemi", action="store_true",
                help="require the two legs in opposite hemispheres (the Run-2 H->inv cut; "
                     "always true for the m_jj picker, a real cut for --pair=lead outputs)")
ap.add_argument("--out", required=True, help="output stem (.pdf and .png)")
args = ap.parse_args()

t = uproot.open(args.file)["events"]
cols = ["wide_pair_mjj", "wide_pair_deta",
        "wide_legA_hs", "wide_legA_pu", "wide_legA_abseta",
        "wide_legB_hs", "wide_legB_pu", "wide_legB_abseta"]
has_hemi = "wide_pair_same_hemi" in t.keys()
has_zpt  = "z_pt" in t.keys()
a = t.arrays(cols + (["wide_pair_same_hemi"] if has_hemi else []) + (["z_pt"] if has_zpt else []), library="np")
if args.zpt_min > 0 and not has_zpt:
    raise SystemExit("--zpt-min needs the z_pt column; this file predates it")
sel = (a["wide_pair_mjj"] >= MJJ_EDGES[0]) & (a["wide_pair_deta"] > args.deta)
n_same = int((sel & (a["wide_pair_same_hemi"] > 0.5)).sum()) if has_hemi else 0
if args.opphemi and has_hemi:
    sel &= a["wide_pair_same_hemi"] < 0.5
n_before_zpt = int(sel.sum())
if args.zpt_min > 0:
    sel &= a["z_pt"] >= args.zpt_min
fA, fB = a["wide_legA_abseta"] >= SPLIT, a["wide_legB_abseta"] >= SPLIT
hA, hB = a["wide_legA_hs"] > .5, a["wide_legB_hs"] > .5
pA, pB = a["wide_legA_pu"] > .5, a["wide_legB_pu"] > .5
xA, xB = ~hA & ~pA, ~hB & ~pB
r1   = fA & fB & ((hA & pB) | (hB & pA))
r2   = (fA & pA & ~fB & hB) | (fB & pB & ~fA & hA)
hspu = ((hA & pB) | (hB & pA)) & ~r1 & ~r2
hshs = hA & hB
ppff = pA & pB & fA & fB
ppo  = pA & pB & ~ppff
nei  = xA | xB
band = np.full(len(sel), -1)
for i, m in enumerate([r1, r2, hspu, hshs, ppff, ppo, nei]):
    band[(band < 0) & m] = i
assert (band[sel] >= 0).all(), "band ladder not exhaustive"
mjj = np.minimum(a["wide_pair_mjj"], 0.5 * (MJJ_EDGES[-2] + MJJ_EDGES[-1]))  # fold overflow

ROOT.gROOT.SetBatch(True)
ROOT.gInterpreter.Declare('#include "src/AtlasStyle.h"\n#include "src/AtlasLabels.h"')
ROOT.SetAtlasStyle()
ROOT.gStyle.SetOptStat(0)
edges = np.array(MJJ_EDGES, dtype=float)
tot = ROOT.TH1D("tot", "", len(edges) - 1, edges)
for v in mjj[sel]: tot.Fill(v)
hists, stk = [], ROOT.THStack("stk", "")
for i, (key, lab, hexc) in enumerate(BANDS):
    h = ROOT.TH1D(f"h_{key}", "", len(edges) - 1, edges)
    for v in mjj[sel & (band == i)]: h.Fill(v)
    for b in range(1, h.GetNbinsX() + 1):
        n, T = h.GetBinContent(b), tot.GetBinContent(b)
        h.SetBinContent(b, n / T if T > 0 else 0.0); h.SetBinError(b, 0.0)
    h.SetFillColor(ROOT.TColor.GetColor(hexc)); h.SetLineColor(ROOT.kBlack); h.SetLineWidth(1)
    hists.append(h); stk.Add(h)

# Two pads: the composition on top, the raw event population per column
# underneath on a log scale -- without it the eye reads a clean fraction in a
# column holding three events as though it meant something.
c = ROOT.TCanvas("c", "c", 800, 780)
pTop = ROOT.TPad("pTop", "", 0.0, 0.30, 1.0, 1.0)
pBot = ROOT.TPad("pBot", "", 0.0, 0.00, 1.0, 0.30)
pTop.SetBottomMargin(0.025); pTop.SetTopMargin(0.06)
pBot.SetTopMargin(0.04);     pBot.SetBottomMargin(0.40)
pTop.Draw(); pBot.Draw()

pTop.cd()
stk.Draw("HIST")
stk.GetXaxis().SetLabelSize(0)
stk.GetXaxis().SetTitle("")
stk.GetYaxis().SetTitle("Fraction of events")
stk.GetYaxis().SetTitleOffset(1.25)
stk.SetMaximum(2.2); stk.SetMinimum(0.0)
pTop.Modified()
n = int(sel.sum())
# ATLASLabel's text offset scales with the pad aspect, so on the tall pad the
# gap after "ATLAS" opens up; draw the word with the helper and the rest by hand.
ROOT.ATLASLabel(0.18, 0.895, "")
al = ROOT.TLatex(); al.SetNDC(); al.SetTextFont(42); al.SetTextSize(0.045)
al.DrawLatex(0.305, 0.895, "Simulation Internal")
ROOT.ATLASEnergyLabel(0.18, 0.85, args.label)
p = ROOT.TLatex(); p.SetNDC(); p.SetTextFont(42); p.SetTextSize(0.028); p.SetTextColor(ROOT.kGray + 2)
line3 = f"|#Delta#eta| > {args.deta:g}" + (",  #eta_{1}#eta_{2} < 0" if args.opphemi else "") \
        + (f",  p_{{T}}^{{ll}} > {args.zpt_min:g} GeV" if args.zpt_min > 0 else "") \
        + f",   forward = |#eta| #geq {SPLIT} (no upper edge, split at {SPLIT})"
p.DrawLatex(0.18, 0.81, line3)
p.DrawLatex(0.18, 0.78, (args.tag + "   " if args.tag else "") + f"#font[62]{{{n:,}}} events in the plot")
leg = ROOT.TLegend(0.17, 0.545, 0.93, 0.765)
ROOT.StyleLegend(leg); leg.SetNColumns(2)
for h, (_, lab, _) in zip(hists, BANDS): leg.AddEntry(h, lab, "f")
leg.Draw()

pBot.cd(); pBot.SetLogy(True)
tot.SetFillColor(ROOT.TColor.GetColor("#DCE6F2")); tot.SetLineColor(ROOT.kBlack); tot.SetLineWidth(1)
tot.GetXaxis().SetTitle("m_{jj} [GeV]   (last bin includes overflow)")
tot.GetYaxis().SetTitle("Events")
tot.GetXaxis().SetTitleSize(0.13); tot.GetXaxis().SetLabelSize(0.11); tot.GetXaxis().SetTitleOffset(1.25)
tot.GetYaxis().SetTitleSize(0.11); tot.GetYaxis().SetLabelSize(0.10); tot.GetYaxis().SetTitleOffset(0.5)
tot.GetYaxis().SetNdivisions(503)
lo = max(1.0, tot.GetMinimum(0.0) / 3.0); hi = tot.GetMaximum() * 60.0
tot.SetMinimum(lo); tot.SetMaximum(hi)
tot.Draw("HIST")
# the per-column count, written above each bar
# staggered on alternate bins so the narrow 250 GeV columns do not collide
q = ROOT.TLatex(); q.SetTextFont(42); q.SetTextSize(0.075); q.SetTextAlign(21)
for b in range(1, tot.GetNbinsX() + 1):
    v = tot.GetBinContent(b)
    if v <= 0: continue
    # first column: left-anchored just inside the frame, so a five-digit
    # count is not clipped by the axis
    q.SetTextAlign(11 if b == 1 else 21)
    x = tot.GetXaxis().GetBinLowEdge(b) + 15 if b == 1 else tot.GetXaxis().GetBinCenter(b)
    q.DrawLatex(x, v * (5.5 if b % 2 else 2.0), f"{int(v):,}")
pBot.RedrawAxis()

c.Print(args.out + ".pdf"); c.Print(args.out + ".png")
print(f"{args.file}: {n:,} pairs with m_jj >= {MJJ_EDGES[0]} and |dEta| > {args.deta:g}"
      + (f" (pT(ll) > {args.zpt_min:g} keeps {n:,} of {n_before_zpt:,})" if args.zpt_min > 0 else "")
      + (f" ({n_same:,} same-hemisphere pairs {'removed' if args.opphemi else 'KEPT'})" if has_hemi else "") + "; "
      + ", ".join(f"{k} {100.0*(band[sel]==i).mean():.1f}%" for i, (k, _, _) in enumerate(BANDS)))
print("wrote", args.out + ".pdf / .png")
