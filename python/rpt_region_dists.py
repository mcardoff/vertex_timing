#!/usr/bin/env python3
"""R_pT in the wide-window VBS regions: the z-only (ITk-only) distribution
against the timed one, one set of pages per timing case, from the per-event
tree rpt_v5_hist writes beside its histograms (<prefix>rpt_v5_regions.root).

The regions are rpt_v5_hist's WIDE window (see RegionRow there): the max-m_jj
VBS pair, forward = |eta| > 2.4 with NO upper edge, central = |eta| < 2.4,
paper HS/PU labels. Legs past |eta| 3.8 are included on purpose. Past 4.0 a
jet has no tracks and sits at R_pT = 0 under every case, and the |eta|-band
page separates those legs out.

    PYTHONNOUSERSITE=1 ~/.venv-hgtd/bin/python python/rpt_region_dists.py \\
        condor/vbf/vbf_rpt_v5_regions.root --sample vbf \\
        --hist-file condor/vbf/vbf_rpt_v5_hist.root

--hist-file also runs the standing consistency check. The tree's `core` rows
(the narrow 2.4-3.8 region the _r1/_r2 histograms use) must reproduce
HS_*_r1 / PU_*_r1 / PU_*_r2 bin for bin, for every scenario, or nothing is
drawn.

Run from the repo root (the ATLAS style headers are included by relative path).
"""
import argparse, os
import numpy as np
import uproot
import ROOT

SCEN = ["zonly", "hgtd", "trkptz", "waves", "waves_ideal", "truth", "tzp"]

# (scenario row, page title, recipe line, colour) in presentation order. The
# last two reuse rpt_v5's existing rows, each of which idealises TWO things,
# and every page they appear on says so.
CASES = [
    ("trkptz",      "TRKPTZ t_{0}",                "",                                                        ROOT.kP10Red),
    ("waves",       "WAVeS t_{0}",                 "",                                                        ROOT.kP10Yellow),
    ("hgtd",        "HGTD t_{0} (Athena)",         "",                                                        ROOT.kP10Blue),
    ("truth",       "Truth vertex t_{0}",          "t_{0} = truth #oplus 10 ps,  track times = truth #oplus 30 ps",       ROOT.kP10Brown),
    ("waves_ideal", "Ideal track time assignment", "track times = truth #oplus 30 ps,  t_{0} = cluster closest to truth", ROOT.kP10Violet),
]

# (key, title, region, R_pT column prefix, |eta| column of that leg)
POPS = [
    ("r1_hs", "forward HS leg", 1, "rpt_hs", "hs_eta"),
    ("r1_pu", "forward PU leg", 1, "rpt_pu", "pu_eta"),
    ("r2_pu", "forward PU leg", 2, "rpt_pu", "pu_eta"),
]
REGION_TEXT = {1: "VBS R1: HS + PU leg, both |#eta| > 2.4",
               2: "VBS R2: PU leg |#eta| > 2.4 + HS leg |#eta| < 2.4"}

# |eta| bands of the plotted leg: the band page's two columns, and the finer
# split the summary table carries (ITk tracks, and so R_pT, end at 4.0).
PLOT_BANDS  = [("2.4 < |#eta| < 3.8", 2.4, 3.8), ("|#eta| > 3.8", 3.8, np.inf)]
TABLE_BANDS = [("all", 0.0, np.inf), ("2.4-3.8", 2.4, 3.8),
               ("3.8-4.0", 3.8, 4.0), (">4.0", 4.0, np.inf)]

ZONLY_FILL, ZONLY_LINE = "#D9D9D9", "#8C8C8C"

ap = argparse.ArgumentParser(description=__doc__,
                             formatter_class=argparse.RawDescriptionHelpFormatter)
ap.add_argument("file", help="merged <prefix>rpt_v5_regions.root")
ap.add_argument("--sample", required=True, help="output-name prefix, e.g. vbf / zjets / local")
ap.add_argument("--hist-file", help="matching rpt_v5_hist.root: energy label + core-subset check")
ap.add_argument("--label", help="energy label (default: from --hist-file, else the VBF one)")
ap.add_argument("--out-dir", default="figs/rpt_regions")
ap.add_argument("--bin-width", type=float,
                help="display bin width (default 0.02 for vbf/local, 0.05 otherwise)")
ap.add_argument("--xmax", type=float, default=1.5)
ap.add_argument("--require-fwd-acc-jet", action="store_true",
                help="keep only events with >= 1 jet > 30 GeV in 2.38 < |eta| < 4.0: "
                     "vbs_region_diag's preselection, which reproduces the composition-plot "
                     "event sets exactly (off by default: the regions are defined by the pair)")
args = ap.parse_args()

# vbf / local, or a tagged variant such as vbf_jvtLoose
width = args.bin_width or (0.02 if args.sample.split("_")[0] in ("vbf", "local") else 0.05)
os.makedirs(os.path.join(args.out_dir, "png"), exist_ok=True)

tree = uproot.open(args.file)["regions"]
# jvt_wp exists from 2026-09-25; an older tree was necessarily untagged.
has_jvt = "jvt_wp" in tree.keys()
cols = (["region", "core", "hs_eta", "pu_eta", "n_jets_fwd_acc",
         "vbs_mjj_cut", "vbs_deta_cut", "gate_sigma"] + (["jvt_wp"] if has_jvt else [])
        + [f"rpt_{leg}_{s}" for leg in ("hs", "pu") for s in SCEN])
a = tree.arrays(cols, library="np")
if len(a["region"]) == 0:
    raise SystemExit(f"{args.file}: the regions tree is empty")

# ── Standing check: the tree's narrow subset IS the _r1/_r2 histograms ──────
label = args.label
if args.hist_file:
    hf = uproot.open(args.hist_file)
    if not label and "meta_energy_label" in hf:
        label = str(hf["meta_energy_label"])
    bad = []
    for s in SCEN:
        for key, reg, col in [(f"HS_{s}_r1", 1, f"rpt_hs_{s}"),
                              (f"PU_{s}_r1", 1, f"rpt_pu_{s}"),
                              (f"PU_{s}_r2", 2, f"rpt_pu_{s}")]:
            H = hf[key]
            edges, ref = H.axis().edges(), H.values(flow=True)
            v = a[col][(a["region"] == reg) & a["core"]]
            inside = v[(v >= edges[0]) & (v < edges[-1])]
            cnt, _ = np.histogram(inside, bins=edges)
            mine = np.concatenate([[(v < edges[0]).sum()], cnt, [(v >= edges[-1]).sum()]])
            if not np.array_equal(mine, ref):
                bad.append(key)
    if bad:
        raise SystemExit(f"core-subset check FAILED against {args.hist_file}: {bad}")
    print(f"core-subset check: tree reproduces all {3 * len(SCEN)} region histograms "
          f"of {args.hist_file} bin for bin")
label = label or "#sqrt{s} = 14 TeV, HL-LHC, VBF H#rightarrowinv."

sel = np.ones(len(a["region"]), dtype=bool)
if args.require_fwd_acc_jet:
    sel &= a["n_jets_fwd_acc"] >= 1
mjj_cut, deta_cut = float(a["vbs_mjj_cut"][0]), float(a["vbs_deta_cut"][0])
sel_text = f"max-m_{{jj}} pair, m_{{jj}} > {mjj_cut:g} GeV"
if deta_cut > 0:
    sel_text += f", |#Delta#eta| > {deta_cut:g}"
sel_text += ", p_{T}^{jet} > 30 GeV"
if args.require_fwd_acc_jet:
    sel_text += ", #geq1 jet in 2.38 < |#eta| < 4.0"
jvt_name = ["none", "loose", "tight"][int(a["jvt_wp"][0])] if has_jvt else "none"
if jvt_name != "none":
    # every page must say the pair was formed from tagger-passing jets only
    sel_text += f", JVT + fJVT {jvt_name} before pairing"


def mask(region, eta_col=None, lo=0.0, hi=np.inf):
    m = sel & (a["region"] == region)
    if eta_col is not None:
        e = np.abs(a[eta_col])
        m &= (e >= lo) & (e < hi)
    return m


def stats(z, t):
    """Per-population numbers; z and t are the SAME jets, z-only and timed."""
    n = len(z)
    if n == 0:
        return dict(n=0, mz=np.nan, mt=np.nan, f0z=np.nan, f0t=np.nan, low=np.nan, zer=np.nan)
    return dict(n=n, mz=z.mean(), mt=t.mean(),
                f0z=(z == 0).mean(), f0t=(t == 0).mean(),
                low=(t < z - 1e-12).mean(), zer=((t == 0) & (z > 0)).mean())


keep = []  # PyROOT objects must outlive the canvas Print


def th1(vals, name):
    nb = int(round(args.xmax / width))
    h = ROOT.TH1D(name, "", nb, 0.0, args.xmax)
    cnt, _ = np.histogram(vals[vals < args.xmax], bins=nb, range=(0.0, args.xmax))
    for b, n in enumerate(cnt, start=1):
        h.SetBinContent(b, float(n))
        h.SetBinError(b, float(np.sqrt(n)))
    keep.append(h)
    return h


def style_pair(hz, ht, colour):
    hz.SetFillColor(ROOT.TColor.GetColor(ZONLY_FILL))
    hz.SetLineColor(ROOT.TColor.GetColor(ZONLY_LINE))
    hz.SetLineWidth(1)
    ht.SetLineColor(colour)
    ht.SetLineWidth(3)
    ht.SetFillStyle(0)


ROOT.gROOT.SetBatch(True)
ROOT.gErrorIgnoreLevel = ROOT.kWarning  # one "Info in <TCanvas::Print>" per page otherwise
ROOT.gInterpreter.Declare('#include "src/AtlasStyle.h"\n#include "src/AtlasLabels.h"')
ROOT.SetAtlasStyle()
ROOT.gStyle.SetOptStat(0)

pdf = os.path.join(args.out_dir, f"{args.sample}_rpt_region_dists.pdf")
c = ROOT.TCanvas("c", "c", 800, 820)
c.Print(pdf + "[")
table = []  # (case title, pop, band, stats)


def header(x, y, size=0.045):
    # ATLASLabel's gap after "ATLAS" scales with the pad aspect, so draw the
    # word with the helper and the rest by hand (as vbs_region_stack.py does).
    ROOT.ATLASLabel(x, y, "")
    t = ROOT.TLatex(); t.SetNDC(); t.SetTextFont(42); t.SetTextSize(size)
    t.DrawLatex(x + 0.125, y, "Simulation Internal")
    ROOT.ATLASEnergyLabel(x, y - 0.045, label)


def main_page(case, pop):
    s, title, recipe, colour = case
    key, ptitle, region, prefix, eta_col = pop
    m = mask(region)
    z, t = a[f"{prefix}_zonly"][m], a[f"{prefix}_{s}"][m]
    st = stats(z, t)
    hz, ht = th1(z, f"hz_{s}_{key}"), th1(t, f"ht_{s}_{key}")
    style_pair(hz, ht, colour)

    c.Clear()
    top = ROOT.TPad(f"top_{s}_{key}", "", 0.0, 0.30, 1.0, 1.0)
    bot = ROOT.TPad(f"bot_{s}_{key}", "", 0.0, 0.00, 1.0, 0.30)
    top.SetBottomMargin(0.02); top.SetTopMargin(0.05)
    bot.SetTopMargin(0.03);    bot.SetBottomMargin(0.38)
    for p in (top, bot):
        p.SetLeftMargin(0.14); p.SetRightMargin(0.05)
    top.Draw(); bot.Draw(); keep.extend([top, bot])

    top.cd(); top.SetLogy(True)
    # The text block fills the upper half of the pad from x ~ R_pT 0.06
    # rightwards. Bins left of that (the R_pT = 0 spike) may reach the top;
    # everything to the right has to stay under the text on the log scale.
    left = right = 1.0
    for b in range(1, hz.GetNbinsX() + 1):
        v = max(hz.GetBinContent(b), ht.GetBinContent(b))
        if hz.GetXaxis().GetBinLowEdge(b) < 0.06: left = max(left, v)
        else:                                     right = max(right, v)
    hz.SetMaximum(max(3.0 * left, 5e3 * right)); hz.SetMinimum(0.5)
    hz.GetYaxis().SetTitle(f"Jets / {width:g}")
    hz.GetYaxis().SetTitleOffset(1.35)
    hz.GetXaxis().SetLabelSize(0)
    hz.Draw("HIST")
    ht.Draw("HIST SAME")
    top.RedrawAxis()

    header(0.18, 0.885)
    tx = ROOT.TLatex(); tx.SetNDC(); tx.SetTextFont(42); tx.SetTextSize(0.034)
    tx.DrawLatex(0.18, 0.785, f"{REGION_TEXT[region]}:  {ptitle}")
    tx.SetTextSize(0.029); tx.SetTextColor(ROOT.kGray + 2)
    tx.DrawLatex(0.18, 0.745, sel_text)
    tx.SetTextColor(ROOT.kBlack); tx.SetTextSize(0.036); tx.SetTextFont(62)
    tx.DrawLatex(0.18, 0.695, title)
    if recipe:
        tx.SetTextFont(42); tx.SetTextSize(0.029); tx.SetTextColor(ROOT.kGray + 2)
        tx.DrawLatex(0.18, 0.657, recipe)
    over = (z >= args.xmax).mean() if len(z) else 0.0
    tx.SetTextFont(42); tx.SetTextSize(0.027); tx.SetTextColor(ROOT.kGray + 2)
    tx.DrawLatex(0.18, 0.617 if recipe else 0.657,
                 f"{st['n']:,} jets" + (f";  {100 * over:.1f}% above R_{{pT}} = {args.xmax:g} (not drawn)"
                                        if over > 0 else ""))
    leg = ROOT.TLegend(0.17, 0.44, 0.93, 0.58 if recipe else 0.62)
    ROOT.StyleLegend(leg, 0.031)
    leg.AddEntry(hz, f"ITk-only (z only):  #LTR_{{pT}}#GT = {st['mz']:.3f},  "
                     f"R_{{pT}} = 0 for {100 * st['f0z']:.1f}%", "f")
    leg.AddEntry(ht, f"{title}:  #LTR_{{pT}}#GT = {st['mt']:.3f},  "
                     f"R_{{pT}} = 0 for {100 * st['f0t']:.1f}%", "l")
    leg.Draw(); keep.append(leg)

    bot.cd()
    r = ht.Clone(f"ratio_{s}_{key}"); keep.append(r)
    r.Divide(hz)
    for b in range(1, r.GetNbinsX() + 1):
        r.SetBinError(b, 0.0)
    r.SetMinimum(0.0); r.SetMaximum(1.99)
    r.SetLineColor(colour); r.SetLineWidth(2)
    r.GetYaxis().SetTitle("timed / z-only")
    r.GetXaxis().SetTitle("R_{pT}")
    for ax, (ls, ts) in ((r.GetXaxis(), (0.11, 0.13)), (r.GetYaxis(), (0.10, 0.10))):
        ax.SetLabelSize(ls); ax.SetTitleSize(ts)
    r.GetXaxis().SetTitleOffset(1.2)
    r.GetYaxis().SetTitleOffset(0.6); r.GetYaxis().SetNdivisions(504)
    r.Draw("HIST")
    one = ROOT.TLine(0.0, 1.0, args.xmax, 1.0); one.SetLineStyle(2); one.SetLineColor(ROOT.kGray + 1)
    one.Draw(); keep.append(one)

    c.Print(pdf)
    c.Print(os.path.join(args.out_dir, "png", f"{args.sample}_{s}_{key}.png"))


def band_page(case):
    """3 x 2 grid: rows R1-HS / R1-PU / R2-PU, columns the leg's |eta| band."""
    s, title, recipe, colour = case
    c.Clear()
    head = ROOT.TPad(f"head_{s}", "", 0.0, 0.92, 1.0, 1.0); head.Draw(); keep.append(head)
    head.cd()
    t = ROOT.TLatex(); t.SetNDC(); t.SetTextSize(0.30)
    t.SetTextFont(72); t.DrawLatex(0.03, 0.55, "ATLAS")
    t.SetTextFont(42); t.DrawLatex(0.135, 0.55, "Simulation Internal")
    t.SetTextFont(62); t.DrawLatex(0.50, 0.55, title)
    t.SetTextFont(42); t.SetTextSize(0.22); t.SetTextColor(ROOT.kGray + 2)
    t.DrawLatex(0.03, 0.15, label.replace("#sqrt{s} = 14 TeV, ", "") + ";   " + sel_text
                + (";   " + recipe if recipe else ""))
    c.cd()
    grid = ROOT.TPad(f"grid_{s}", "", 0.0, 0.0, 1.0, 0.92); grid.Draw(); keep.append(grid)
    grid.Divide(2, 3, 0.002, 0.002)
    for row, pop in enumerate(POPS):
        key, ptitle, region, prefix, eta_col = pop
        for col, (btxt, lo, hi) in enumerate(PLOT_BANDS):
            pad = grid.cd(1 + 2 * row + col)
            pad.SetLogy(True); pad.SetTopMargin(0.12); pad.SetBottomMargin(0.16)
            pad.SetLeftMargin(0.15); pad.SetRightMargin(0.04)
            m = mask(region, eta_col, lo, hi)
            z, x = a[f"{prefix}_zonly"][m], a[f"{prefix}_{s}"][m]
            st = stats(z, x)
            hz, ht = th1(z, f"bz_{s}_{key}_{col}"), th1(x, f"bt_{s}_{key}_{col}")
            style_pair(hz, ht, colour)
            ymax = max(hz.GetMaximum(), ht.GetMaximum(), 1.0)
            hz.SetMaximum(ymax * 30); hz.SetMinimum(0.5)
            hz.GetXaxis().SetTitle("R_{pT}"); hz.GetYaxis().SetTitle(f"Jets / {width:g}")
            for ax in (hz.GetXaxis(), hz.GetYaxis()):
                ax.SetLabelSize(0.06); ax.SetTitleSize(0.065)
            hz.GetYaxis().SetTitleOffset(1.1); hz.GetXaxis().SetTitleOffset(1.0)
            hz.Draw("HIST"); ht.Draw("HIST SAME"); pad.RedrawAxis()
            tt = ROOT.TLatex(); tt.SetNDC(); tt.SetTextFont(42); tt.SetTextSize(0.062)
            tt.DrawLatex(0.17, 0.90, f"R{region} {ptitle}, {btxt}")
            tt.SetTextSize(0.055); tt.SetTextColor(ROOT.kGray + 2)
            if st["n"]:
                tt.DrawLatex(0.17, 0.81, f"{st['n']:,} jets;  R_{{pT}}=0: {100 * st['f0z']:.0f}% #rightarrow "
                                         f"{100 * st['f0t']:.0f}%;  lowered {100 * st['low']:.1f}%")
            else:
                tt.DrawLatex(0.17, 0.81, "no jets")
            if row == 0 and col == 0:
                leg = ROOT.TLegend(0.55, 0.52, 0.95, 0.74)
                ROOT.StyleLegend(leg, 0.055)
                leg.AddEntry(hz, "ITk-only", "f"); leg.AddEntry(ht, "timed", "l")
                leg.Draw(); keep.append(leg)
    c.Print(pdf)
    c.Print(os.path.join(args.out_dir, "png", f"{args.sample}_{s}_etabands.png"))


for case in CASES:
    for pop in POPS:
        main_page(case, pop)
    band_page(case)
    s = case[0]
    for pop in POPS:
        key, ptitle, region, prefix, eta_col = pop
        for bname, lo, hi in TABLE_BANDS:
            m = mask(region, eta_col, lo, hi)
            table.append((case[1], f"R{region} {ptitle}", bname,
                          stats(a[f"{prefix}_zonly"][m], a[f"{prefix}_{s}"][m])))
c.Print(pdf + "]")

# ── Summary table ────────────────────────────────────────────────────────────
def latex_to_md(txt):
    return (txt.replace("t_{0}", "t0").replace("#oplus", "(+)")
               .replace("#LT", "<").replace("#GT", ">"))


n_r1, n_r2 = int((sel & (a["region"] == 1)).sum()), int((sel & (a["region"] == 2)).sum())
lines = [f"# R_pT in the wide VBS regions: {args.sample}", "",
         f"Source: `{args.file}`  ",
         f"Selection: {latex_to_md(sel_text).replace('_{jj}', '_jj').replace('_{T}^{jet}', 'T(jet)')}; "
         f"forward = |eta| > 2.4 (no upper edge); gate {float(a['gate_sigma'][0]):g} sigma  ",
         f"Events: R1 {n_r1:,}, R2 {n_r2:,}"
         + (f" (without the forward-acceptance-jet requirement: R1 {int((a['region'] == 1).sum()):,}, "
            f"R2 {int((a['region'] == 2).sum()):,})" if args.require_fwd_acc_jet else
            f" (with vbs_region_diag's >=1 jet in 2.38<|eta|<4.0: "
            f"R1 {int(((a['region'] == 1) & (a['n_jets_fwd_acc'] >= 1)).sum()):,}, "
            f"R2 {int(((a['region'] == 2) & (a['n_jets_fwd_acc'] >= 1)).sum()):,})"), "",
         "`lowered` = timed R_pT below z-only for the same jet; `zeroed` = timed R_pT = 0 where "
         "z-only was not. Bands are the plotted leg's |eta|.", ""]
cur = None
for title, pop, band, st in table:
    if title != cur:
        cur = title
        lines += ["", f"## {latex_to_md(title)}", "",
                  "| population | |eta| band | jets | <R_pT> z-only | <R_pT> timed | "
                  "R_pT=0 z-only | R_pT=0 timed | lowered | zeroed |",
                  "|---|---|---:|---:|---:|---:|---:|---:|---:|"]
    if st["n"] == 0:
        lines.append(f"| {pop} | {band} | 0 | | | | | | |")
        continue
    lines.append(f"| {pop} | {band} | {st['n']:,} | {st['mz']:.3f} | {st['mt']:.3f} | "
                 f"{100 * st['f0z']:.1f}% | {100 * st['f0t']:.1f}% | "
                 f"{100 * st['low']:.1f}% | {100 * st['zer']:.1f}% |")
md = os.path.join(args.out_dir, f"{args.sample}_rpt_region_summary.md")
with open(md, "w") as f:
    f.write("\n".join(lines) + "\n")
print("\n".join(lines))
print(f"\nwrote {pdf}\n      {md}\n      {os.path.join(args.out_dir, 'png')}/{args.sample}_*.png")
