#!/usr/bin/env python3
"""The VBS-pair region composition (R1, R2, HS+HS, ...) per m_jj column, in exactly the
events and binning of the timing-veto efficiency figures: util/vbs_time_veto trees, the
recoil cuts (VBF truth E_T^miss, Z+jets pT(ll)), and the HGTD-acceptance case (--acceptance,
timed tracks as vbs_time_veto_plot.py). Generator-weighted fractions; the MC event count of
each column underneath. The region ladder and colours are vbs_region_stack.py's (forward =
|eta_axis| >= 2.4), so the plots read as the region talk's.

    PYTHONNOUSERSITE=1 PYTHONPATH=/opt/homebrew/Cellar/root/6.40.04/lib/root \\
    ~/.venv-hgtd/bin/python python/vbs_case_regions.py \\
        --sig condor/vbf/vbf_jvtLoose_vbs_time_veto.root --bkg condor/zjets/zjets_jvtLoose_vbs_time_veto.root \\
        --acceptance both --out figs/time_veto/jvtLoose_recoil100_accboth_t_regions
writes <out>_{sig,bkg}.(pdf|png).
"""
import argparse, os
import numpy as np
import uproot
import vbs_time_veto_common as C

BANDS = [   # vbs_region_stack.py's, bottom -> top
    ("R1: fwd HS + fwd PU",        "#3F8FFF"),
    ("R2: fwd PU + central HS",    "#B3231C"),
    ("HS + PU, other #eta config", "#F5A02A"),
    ("HS + HS (genuine pair)",     "#7BCB6B"),
    ("PU + PU, both forward",      "#636363"),
    ("PU + PU, other #eta config", "#CFCFCF"),
    ("A leg carries neither label", "#C79BF5"),
]
SPLIT = 2.4
MJJ_EDGES = [200, 500, 750, 1000, 1500, 2000, 2500, 3000, 4000]   # vbs_time_veto_plot.py's

ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
ap.add_argument("--sig", required=True)
ap.add_argument("--bkg", required=True)
ap.add_argument("--met-min", type=float, default=100.0)
ap.add_argument("--zpt-min", type=float, default=100.0)
ap.add_argument("--out", required=True)
C.add_acceptance_args(ap)
args = ap.parse_args()


def bands(e):
    fA, fB = np.abs(e["a_eta"]) >= SPLIT, np.abs(e["b_eta"]) >= SPLIT
    hA, hB, pA, pB = e["a_hs"], e["b_hs"], e["a_pu"], e["b_pu"]
    r1 = fA & fB & ((hA & pB) | (hB & pA))
    r2 = (fA & pA & ~fB & hB) | (fB & pB & ~fA & hA)
    hspu = ((hA & pB) | (hB & pA)) & ~r1 & ~r2
    hshs = hA & hB
    ppff = pA & pB & fA & fB
    ppo = pA & pB & ~ppff
    nei = (~hA & ~pA) | (~hB & ~pB)
    b = np.full(len(fA), -1)
    for i, m in enumerate([r1, r2, hspu, hshs, ppff, ppo, nei]):
        b[(b < 0) & m] = i
    assert (b >= 0).all()
    return b


ROOT = C.setup_root(args.out)
edges = np.array(MJJ_EDGES, float)
col = C.acc_column(args.acc_count)
for role, path, rc, cut, label in [
        ("sig", args.sig, "met_truth", args.met_min, "#sqrt{s} = 14 TeV, HL-LHC, VBF H#rightarrowinv."),
        ("bkg", args.bkg, "z_pt", args.zpt_min, "#sqrt{s} = 14 TeV, HL-LHC, Z(#rightarrowll)+jets")]:
    cols = ["weight", "mjj", rc, "a_eta", "b_eta", "a_hs", "b_hs", "a_pu", "b_pu", "a_" + col, "b_" + col]
    e = uproot.open(path)["events"].arrays(cols, library="np")
    keep = (e[rc] > cut) & C.acceptance_mask(e["a_" + col], e["b_" + col], args.acceptance, args.acc_min_trk)
    e = {k: v[keep] for k, v in e.items()}
    b, w = bands(e), e["weight"].astype(float)
    idx = C.mjj_index(e["mjj"], edges)
    c = ROOT.TCanvas("c_" + role, "", 800, 780); C.KEEP.append(c)
    pTop = ROOT.TPad("pTop" + role, "", 0.0, 0.30, 1.0, 1.0); pBot = ROOT.TPad("pBot" + role, "", 0.0, 0.0, 1.0, 0.30)
    pTop.SetBottomMargin(0.025); pTop.SetTopMargin(0.06); pBot.SetTopMargin(0.04); pBot.SetBottomMargin(0.40)
    pTop.Draw(); pBot.Draw(); C.KEEP.extend([pTop, pBot])
    pTop.cd()
    stk = ROOT.THStack("stk" + role, ""); C.KEEP.append(stk)
    hs = []
    for i, (lab, hexc) in enumerate(BANDS):
        h = ROOT.TH1D(f"h{role}{i}", "", len(edges) - 1, edges); C.KEEP.append(h)
        for k in range(len(edges) - 1):
            m = idx == k
            W = w[m].sum()
            h.SetBinContent(k + 1, w[m & (b == i)].sum() / W if W > 0 else 0.0)
        h.SetFillColor(ROOT.TColor.GetColor(hexc)); h.SetLineColor(ROOT.kBlack); h.SetLineWidth(1)
        stk.Add(h); hs.append(h)
    stk.Draw("HIST")
    stk.GetXaxis().SetLabelSize(0); stk.GetYaxis().SetTitle("Fraction of events (weighted)")
    stk.GetYaxis().SetTitleOffset(1.25); stk.SetMaximum(1.85); stk.SetMinimum(0.0)
    ytxt = C.atlas_block(label, ["**" + C.CASE_TEXT[args.acceptance]])
    leg = ROOT.TLegend(0.17, ytxt - 0.22, 0.93, ytxt - 0.01); ROOT.StyleLegend(leg); leg.SetNColumns(2); C.KEEP.append(leg)
    for h, (lab, _) in zip(hs, BANDS): leg.AddEntry(h, lab, "f")
    leg.Draw(); pTop.RedrawAxis()
    pBot.cd(); pBot.SetLogy(True)
    tot = ROOT.TH1D("tot" + role, "", len(edges) - 1, edges); C.KEEP.append(tot)
    counts = [int((idx == k).sum()) for k in range(len(edges) - 1)]
    for k, n in enumerate(counts): tot.SetBinContent(k + 1, n)
    tot.SetFillColor(ROOT.TColor.GetColor("#DCE6F2")); tot.SetLineColor(ROOT.kBlack)
    tot.GetXaxis().SetTitle("m_{jj} [GeV]   (last bin includes overflow)"); tot.GetYaxis().SetTitle("MC events")
    tot.GetXaxis().SetTitleSize(0.13); tot.GetXaxis().SetLabelSize(0.11); tot.GetXaxis().SetTitleOffset(1.25)
    tot.GetYaxis().SetTitleSize(0.11); tot.GetYaxis().SetLabelSize(0.10); tot.GetYaxis().SetTitleOffset(0.5)
    tot.GetYaxis().SetNdivisions(503)
    tot.SetMinimum(max(0.5, min([n for n in counts if n > 0] or [1]) / 3.0)); tot.SetMaximum(max(counts + [1]) * 60.0)
    tot.Draw("HIST")
    q = ROOT.TLatex(); q.SetTextFont(42); q.SetTextSize(0.075); C.KEEP.append(q)
    top = max(counts + [1])
    for k, n in enumerate(counts):
        if n <= 0: continue
        q.SetTextAlign(11 if k == 0 else 21)
        x = edges[0] + 0.015 * (edges[-1] - edges[0]) if k == 0 else 0.5 * (edges[k] + edges[k + 1])
        q.DrawLatex(x, top * (12.0 if k % 2 == 0 else 3.0), f"{n:,}")
    pBot.RedrawAxis()
    c.Print(f"{args.out}_{role}.pdf"); c.Print(f"{args.out}_{role}.png")
    shares = ", ".join(f"{lab.split(':')[0]} {100 * w[b == i].sum() / w.sum():.1f}%" for i, (lab, _) in enumerate(BANDS))
    print(f"wrote {args.out}_{role}.pdf / .png  ({int(keep.sum()):,} MC events; {shares})")
