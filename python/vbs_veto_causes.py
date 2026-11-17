#!/usr/bin/env python3
"""What the jet-vs-t0 timing veto removes, and what it leaves, per m_jj column.

--what removed : the fraction of each m_jj column's events the veto removes, stacked by
                 the cause (exclusive, first match wins):
                   1. a pileup tagging jet fails the test
                   2. only hard-scatter jets fail, and the t0 is wrong (> 3 sigma_t0 from truth)
                   3. only hard-scatter jets fail, t0 right, the failing jet's time is set by
                      pileup tracks (< 50% of its time cluster's pT from the hard scatter)
                   4. only hard-scatter jets fail, t0 right, its time from its own tracks
                      (the resolution tail)
                   5. only jets without a truth label fail
--what kept    : the events the veto keeps that still have a pileup tagging jet, and why
                 each such jet survived (per pileup leg, weighted):
                   no timed track in HGTD, axis central (|eta| < 2.38)
                   no timed track in HGTD, axis past |eta| = 4
                   axis outside HGTD; timed tracks only at its edge (none within dR < 0.2: no jet time)
                   jet time from one track (below --min-trk, not used)
                   timed, within 3 sigma of t0
The veto is vbs_time_veto_plot.py's (vbs_time_veto_common.t0_veto), one jet-time
calibration from the signal's HS legs for both samples; generator-weighted.

    PYTHONNOUSERSITE=1 PYTHONPATH=/opt/homebrew/Cellar/root/6.40.04/lib/root \\
    ~/.venv-hgtd/bin/python python/vbs_veto_causes.py \\
        --sig condor/vbf/vbf_jvtLoose_vbs_time_veto.root --bkg condor/zjets/zjets_jvtLoose_vbs_time_veto.root \\
        --t0 waves --min-trk 2 --what removed --out figs/time_veto_k2/jvtLoose_recoil100_veto_causes_waves
writes <out>_{sig,bkg}.(pdf|png) and <out>.md (the numbers).
"""
import argparse, os
import numpy as np
import uproot
import vbs_time_veto_common as C

ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
ap.add_argument("--sig", required=True)
ap.add_argument("--bkg", required=True)
ap.add_argument("--t0", default="waves", choices=["trkptz", "waves", "hgtd", "truth"])
ap.add_argument("--min-trk", type=int, default=2)
ap.add_argument("--nsigma", type=float, default=3.0)
ap.add_argument("--met-min", type=float, default=100.0)
ap.add_argument("--zpt-min", type=float, default=100.0)
ap.add_argument("--mjj-edges", default="200,1000,2000,4000")
ap.add_argument("--what", choices=["removed", "kept"], default="removed")
ap.add_argument("--tag", default="JVT + fJVT loose")
ap.add_argument("--out", required=True)
args = ap.parse_args()

T0NAME = {"trkptz": "TRKPTZ", "waves": "WAVeS", "hgtd": "HGTD (Athena)", "truth": "truth"}[args.t0]
REMOVED = [("Rejected a PU jet",                                  "#555555"),
           ("Rejected a HS jet due to wrong t_{0}",               "#D55E00"),
           ("Rejected a HS jet due to PU tracks",                 "#CC79A7"),
           ("Rejected a HS jet due to track time misassignment",  "#2E9E4F"),
           ("only jets without a truth label fail",               "#C79BF5")]
KEPT = [("no timed track in HGTD: axis central (|#eta| < 2.38)",         "#CFCFCF"),
        ("no timed track in HGTD: axis past |#eta| = 4",                 "#8C8C8C"),
        ("axis outside HGTD, timed tracks only at its edge: no jet time", "#56B4E9"),
        (f"jet time from 1 track (below #geq {args.min_trk}, not used)",  "#E69F00"),
        (f"timed, within {args.nsigma:g}#sigma of t_{{0}}",               "#B3231C")]
edges = np.array([float(x) for x in args.mjj_edges.split(",")])
nb = len(edges) - 1

cols_needed = None
def load(path, rc, cut):
    t = uproot.open(path)["events"]
    e = t.arrays([k for k in t.keys() if not k.startswith("jet_")], library="np")
    k = e[rc] > cut
    return {c: v[k] for c, v in e.items()}, e


S, S_all = load(args.sig, "met_truth", args.met_min)
B, _ = load(args.bkg, "z_pt", args.zpt_min)
cal = C.measure(S_all, True, "signal HS legs")


def removed_classes(e):
    vet, fa, fb = C.t0_veto(e, cal, args.t0, args.min_trk, args.nsigma)
    pu = (fa & e["a_pu"]) | (fb & e["b_pu"])
    hsf = {L: (f & e[f"{L}_hs"]) for L, f in (("a", fa), ("b", fb))}
    hs_any = hsf["a"] | hsf["b"]
    if args.t0 == "truth":
        t0bad = np.zeros(len(vet), bool)
    else:
        s0 = e[f"infl_{args.t0}"] * e[f"sig0_{args.t0}"]
        t0bad = np.abs(e[f"t0_{args.t0}"] - e["t0_truth"]) > args.nsigma * np.where(s0 > 0, s0, np.inf)
    by_pu = hsf["a"] & (e["a_core_hsfrac"] >= 0) & (e["a_core_hsfrac"] < 0.5) | \
            hsf["b"] & (e["b_core_hsfrac"] >= 0) & (e["b_core_hsfrac"] < 0.5)
    c = np.full(len(vet), -1)
    for i, m in enumerate([pu, ~pu & hs_any & t0bad, ~pu & hs_any & ~t0bad & by_pu,
                           ~pu & hs_any & ~t0bad & ~by_pu, ~pu & ~hs_any]):
        c[(c < 0) & vet & m] = i
    assert ((c >= 0) == vet).all()
    return c


def kept_legs(e):
    """(event index, class) per pileup leg of a kept event"""
    vet, _, _ = C.t0_veto(e, cal, args.t0, args.min_trk, args.nsigma)
    out_i, out_c = [], []
    for L in "ab":
        m = ~vet & e[f"{L}_pu"]
        eta, nt, nc, n = np.abs(e[f"{L}_eta"]), e[f"{L}_hgtd_ntimed"], e[f"{L}_core_ncand"], e[f"{L}_core_n"]
        c = np.full(len(m), -1)
        for i, mm in enumerate([(nt == 0) & (eta < 2.38), (nt == 0) & (eta >= 2.38),
                                (nt > 0) & (nc == 0), (nc > 0) & (n < args.min_trk), n >= args.min_trk]):
            c[(c < 0) & mm] = i
        sel = np.flatnonzero(m)
        assert (c[sel] >= 0).all()
        out_i.append(sel); out_c.append(c[sel])
    return np.concatenate(out_i), np.concatenate(out_c)


ROOT = C.setup_root(args.out)
md = [f"# {os.path.basename(args.out)}: {'what the veto removes' if args.what == 'removed' else 'pileup legs the veto keeps'}"
      f" ({T0NAME} t0, jet time from >= {args.min_trk} tracks, {args.nsigma:g} sigma)", ""]
for role, e, label in [("sig", S, "#sqrt{s} = 14 TeV, HL-LHC, VBF H#rightarrowinv."),
                       ("bkg", B, "#sqrt{s} = 14 TeV, HL-LHC, Z(#rightarrowll)+jets")]:
    w = e["weight"].astype(float)
    col = C.mjj_index(e["mjj"], edges)
    if args.what == "removed":
        cls, bands = removed_classes(e), REMOVED
        # fraction of the COLUMN's events vetoed for each cause
        frac = np.array([[w[(col == k) & (cls == i)].sum() / max(w[col == k].sum(), 1e-12) for k in range(nb)]
                         for i in range(len(bands))])
        nmc = [int(((col == k) & (cls >= 0)).sum()) for k in range(nb)]
        tot_mc = [int((col == k).sum()) for k in range(nb)]
        ytitle = "Fraction of the events vetoed"
        allfrac = [w[cls == i].sum() / w.sum() for i in range(len(bands))]
        md += [f"## {role}: vetoes {100 * w[cls >= 0].sum() / w.sum():.1f}% of the events ({int((cls >= 0).sum())} MC)",
               "| cause | of the vetoed | " + " | ".join(f"{edges[k]:.0f}-{edges[k + 1]:.0f}" for k in range(nb)) + " |",
               "|---|---|" + "---|" * nb]
        V = w[cls >= 0].sum()
        for i, (lab, _) in enumerate(bands):
            md.append(f"| {lab} | {100 * w[cls == i].sum() / V:.0f}% | "
                      + " | ".join(f"{100 * frac[i, k]:.1f}%" for k in range(nb)) + " |")
        md.append("| (column total vetoed) | | " + " | ".join(f"{100 * frac[:, k].sum():.1f}% ({nmc[k]}/{tot_mc[k]} MC)" for k in range(nb)) + " |")
        md.append("")
    else:
        ii, cc = kept_legs(e)
        bands = KEPT
        wl, cl = w[ii], col[ii]
        frac = np.array([[wl[(cl == k) & (cc == i)].sum() / max(wl[cl == k].sum(), 1e-12) for k in range(nb)]
                         for i in range(len(bands))])
        nmc = [int((cl == k).sum()) for k in range(nb)]
        tot_mc = nmc
        vet, _, _ = C.t0_veto(e, cal, args.t0, args.min_trk, args.nsigma)
        anypu = e["a_pu"] | e["b_pu"]
        share = [w[(col == k) & ~vet & anypu].sum() / max(w[(col == k) & ~vet].sum(), 1e-12) for k in range(nb)]
        ytitle = "Fraction of the kept pileup tagging jets"
        md += [f"## {role}: kept events with a pileup tagging jet: {100 * w[~vet & anypu].sum() / w[~vet].sum():.0f}% of the kept weight",
               "| why the pileup leg survived | all | " + " | ".join(f"{edges[k]:.0f}-{edges[k + 1]:.0f}" for k in range(nb)) + " |",
               "|---|---|" + "---|" * nb]
        for i, (lab, _) in enumerate(bands):
            md.append(f"| {lab} | {100 * wl[cc == i].sum() / wl.sum():.0f}% | " + " | ".join(f"{100 * frac[i, k]:.0f}%" for k in range(nb)) + " |")
        md.append("| (kept events with a pileup leg) | | " + " | ".join(f"{100 * share[k]:.0f}% ({nmc[k]} legs)" for k in range(nb)) + " |")
        md.append("")

    c = ROOT.TCanvas("c_" + role, "", 800, 780); C.KEEP.append(c)
    pTop = ROOT.TPad("pTop" + role, "", 0.0, 0.30, 1.0, 1.0); pBot = ROOT.TPad("pBot" + role, "", 0.0, 0.0, 1.0, 0.30)
    pTop.SetBottomMargin(0.025); pTop.SetTopMargin(0.06); pBot.SetTopMargin(0.04); pBot.SetBottomMargin(0.40)
    pTop.Draw(); pBot.Draw(); C.KEEP.extend([pTop, pBot])
    pTop.cd()
    stk = ROOT.THStack("stk" + role, ""); C.KEEP.append(stk)
    hs = []
    for i, (lab, hexc) in enumerate(bands):
        h = ROOT.TH1D(f"h{role}{i}", "", nb, edges); C.KEEP.append(h)
        for k in range(nb): h.SetBinContent(k + 1, frac[i, k])
        h.SetFillColor(ROOT.TColor.GetColor(hexc)); h.SetLineColor(ROOT.kBlack); h.SetLineWidth(1)
        stk.Add(h); hs.append(h)
    top = frac.sum(axis=0).max()
    stk.Draw("HIST")
    stk.GetXaxis().SetLabelSize(0); stk.GetYaxis().SetTitle(ytitle); stk.GetYaxis().SetTitleOffset(1.25)
    ymax = top * 3.0 if args.what == "removed" else 2.25   # the legend sits above the stacks
    stk.SetMinimum(0.0); stk.SetMaximum(ymax)
    head = (f"**What the {T0NAME} t_{{0}} veto removes" if args.what == "removed"
            else f"**Pileup tagging jets the {T0NAME} t_{{0}} veto keeps")
    C.atlas_block(label, [head, f"{args.tag}, jet time from #geq {args.min_trk} tracks, all selected events"])
    ly = (0.41, 0.67) if args.what == "removed" else (0.475, 0.735)
    leg = ROOT.TLegend(0.17, ly[0], 0.93, ly[1]); ROOT.StyleLegend(leg); leg.SetTextSize(0.030); C.KEEP.append(leg)
    for h, (lab, _) in zip(hs, bands): leg.AddEntry(h, lab, "f")
    leg.Draw()
    if args.what == "removed":   # the column's total vetoed fraction, on top of each stack
        q = ROOT.TLatex(); q.SetTextFont(42); q.SetTextSize(0.034); q.SetTextAlign(21); C.KEEP.append(q)
        for k in range(nb):
            q.DrawLatex(0.5 * (edges[k] + edges[k + 1]), frac[:, k].sum() + 0.02 * ymax, f"{100 * frac[:, k].sum():.1f}%")
    pTop.RedrawAxis()
    pBot.cd()
    tot = ROOT.TH1D("tot" + role, "", nb, edges); C.KEEP.append(tot)
    sc = 1e-3 if max(nmc + [1]) >= 10000 else 1.0                     # in 10^3 when large
    for k, n in enumerate(nmc): tot.SetBinContent(k + 1, n * sc)
    tot.SetFillColor(ROOT.TColor.GetColor("#DCE6F2")); tot.SetLineColor(ROOT.kBlack)
    tot.GetXaxis().SetTitle("m_{jj} [GeV]   (last bin includes overflow)")
    tot.GetYaxis().SetTitle(("MC vetoed" if args.what == "removed" else "MC jets") + (" [10^{3}]" if sc != 1.0 else ""))
    tot.GetXaxis().SetTitleSize(0.13); tot.GetXaxis().SetLabelSize(0.11); tot.GetXaxis().SetTitleOffset(1.25)
    tot.GetYaxis().SetTitleSize(0.10); tot.GetYaxis().SetLabelSize(0.10); tot.GetYaxis().SetTitleOffset(0.55)
    tot.SetMinimum(0.0); tot.SetMaximum(max(nmc + [1]) * sc * 2.2); tot.GetYaxis().SetNdivisions(203)
    tot.Draw("HIST")
    q2 = ROOT.TLatex(); q2.SetTextFont(42); q2.SetTextSize(0.085); q2.SetTextAlign(21); C.KEEP.append(q2)
    fmt = (lambda n: f"{n / 1e3:.1f}k") if sc != 1.0 else (lambda n: f"{n:,}")
    for k, n in enumerate(nmc):
        q2.DrawLatex(0.5 * (edges[k] + edges[k + 1]), max(nmc + [1]) * sc * 1.3,
                     f"{fmt(n)} of {fmt(tot_mc[k])}" if args.what == "removed" else fmt(n))
    pBot.RedrawAxis()
    c.Print(f"{args.out}_{role}.pdf"); c.Print(f"{args.out}_{role}.png")
open(args.out + ".md", "w").write("\n".join(md) + "\n")
print("\n".join(md))
