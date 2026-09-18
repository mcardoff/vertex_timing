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
ap.add_argument("--out", required=True, help="output stem (.pdf and .png)")
args = ap.parse_args()

t = uproot.open(args.file)["events"]
a = t.arrays(["wide_pair_mjj", "wide_pair_deta",
              "wide_legA_hs", "wide_legA_pu", "wide_legA_abseta",
              "wide_legB_hs", "wide_legB_pu", "wide_legB_abseta"], library="np")
sel = (a["wide_pair_mjj"] >= MJJ_EDGES[0]) & (a["wide_pair_deta"] > args.deta)
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

c = ROOT.TCanvas("c", "c", 800, 600)
stk.Draw("HIST")
stk.GetXaxis().SetTitle("m_{jj} [GeV]   (last bin includes overflow)")
stk.GetYaxis().SetTitle("Fraction of events")
stk.SetMaximum(2.2); stk.SetMinimum(0.0)
c.Modified()
ROOT.ATLASLabel(0.18, 0.90, "Simulation Internal")
ROOT.ATLASEnergyLabel(0.18, 0.855, args.label)
p = ROOT.TLatex(); p.SetNDC(); p.SetTextFont(42); p.SetTextSize(0.028); p.SetTextColor(ROOT.kGray + 2)
line3 = f"|#Delta#eta| > {args.deta:g},   forward = |#eta| #geq {SPLIT} (no upper edge, split at {SPLIT})"
p.DrawLatex(0.18, 0.815, line3)
if args.tag:
    p.DrawLatex(0.18, 0.785, args.tag)
leg = ROOT.TLegend(0.17, 0.545, 0.93, 0.775)
ROOT.StyleLegend(leg); leg.SetNColumns(2)
for h, (_, lab, _) in zip(hists, BANDS): leg.AddEntry(h, lab, "f")
leg.Draw()
c.Print(args.out + ".pdf"); c.Print(args.out + ".png")
n = int(sel.sum())
print(f"{args.file}: {n:,} pairs with m_jj >= {MJJ_EDGES[0]} and |dEta| > {args.deta:g}; "
      + ", ".join(f"{k} {100.0*(band[sel]==i).mean():.1f}%" for i, (k, _, _) in enumerate(BANDS)))
print("wrote", args.out + ".pdf / .png")
