"""Shared by the util/vbs_time_veto plotting scripts: python/vbs_time_veto_plot.py
(the region-plot selection) and python/vbs_time_veto_sr.py (the VBF H->inv
signal region). Sample normalisation, the jet-time calibration, weighted
efficiencies, and the ROOT figure layouts. ROOT is imported by setup_root(),
not at import time."""
import os
import numpy as np
import awkward as ak

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

# Normalisation to --lumi (3000 fb^-1 by default): yield = w * L * sigma * eff / sum(w_raw).
# sigma x filter eff from AMI's 14 TeV evgen metadata (e8481), sum(w) over EVERY event of
# the RAW (pre-skim) sample: util/scratch/region_weights/sumw.txt, and the table in
# results/rpt_region_distributions.md. The raw sums are what the normalisation needs --
# a skimmed run's own meta sum is of the skim, not the sample.
NORM = {
    "vbf":   dict(xs_pb=4.4316 * 0.5148835, sumw=6497726.944462,
                  note="600026 VBF H->ZZ->4nu MET75, sigma_VBF only: B(H->inv) = 100%"),
    "zjets": dict(xs_pb=2067.2 + 2062.1, sumw=2067163337.673096 + 2061509864.110840,
                  note="601189 + 601190 Z->ee + Z->mumu"),
}

# What a COMPLETE run reads, so a missing shard is caught before it shows up as a
# yield: VBF must be read from the RAW ntuples (the truth record is needed for
# E_T^miss), so its meta sum of weights must equal the raw sum above; Z+jets reads
# the skims, whose event count after the vertex + Z->ll cuts is fixed.
EXPECTED_COMPLETE = {"vbf": dict(n_read=1464400, sumw=6497726.944462), "zjets": dict(n_read=530962),
                     "dijet": dict(n_read=699500)}   # dijet: the RAW JZ0-JZ4 ntuples (--sr-met=jets)

T0_SRC = ["trkptz", "waves", "hgtd", "tzp", "truth"]
# (legend, colour, marker): the region talk's case colours; jet-jet in black.
# Event vetoes filled, jet removal open in the same colour.
STYLE = {
    "jj":        ("Jet vs jet (no t_{0})",                     "#000000",    20),
    "jj_all":    ("Jet vs jet, all timed ghost tracks",        "#8C8C8C",    24),
    "t0_trkptz": ("Jet vs t_{0}: TRKPTZ",                      "kP10Red",    21),
    "t0_waves":  ("Jet vs t_{0}: WAVeS",                       "kP10Yellow", 22),
    "t0_hgtd":   ("Jet vs t_{0}: HGTD (Athena)",               "kP10Blue",   23),
    "t0_tzp":    ("Jet vs t_{0}: TZP",                         "kP10Green",  33),
    "t0_truth":  ("Jet vs t_{0}: truth (perfect t_{0})",       "kP10Brown",  34),
    "rp_jj":     ("Re-pair: time-compatible jets",             "#000000",    24),
    "rp_trkptz": ("Jet removal: TRKPTZ t_{0}",                 "kP10Red",    25),
    "rp_waves":  ("Jet removal: WAVeS t_{0}",                  "kP10Yellow", 26),
    "rp_hgtd":   ("Jet removal: HGTD t_{0} (Athena)",          "kP10Blue",   32),
    "rp_tzp":    ("Jet removal: TZP t_{0}",                    "kP10Green",  27),
    "rp_truth":  ("Jet removal: truth t_{0}",                  "kP10Brown",  28),
    "none":      ("No timing veto",                            "#7F7F7F",    25),
}

# bins of n = tracks in the core jet-time estimate: the pull width falls from
# 2.3 (one track) to 1.3 (nine or more) on VBF, which a single factor averages over
NBINS = [1, 2, 3, 4, 6, 9, 10**6]


def mad(x):
    return 1.4826 * float(np.median(np.abs(x)))


def hs_pulls(e_all):
    """(t_jet - t_HS) / sigma_jet and n for every paper-HS leg with a core time"""
    p, n = [], []
    for L in "ab":
        m = e_all[f"{L}_hs"] & (e_all[f"{L}_core_n"] > 0)
        p.append((e_all[f"{L}_core_t"][m] - e_all["t0_truth"][m]) / e_all[f"{L}_core_sig"][m])
        n.append(e_all[f"{L}_core_n"][m])
    return np.concatenate(p), np.concatenate(n)


class Calib:
    """sigma_jet -> f(n) x sigma_jet; `fac` None means one global factor."""
    def __init__(self, glob, fac=None, source="fixed", nlegs=0):
        self.glob, self.fac, self.source, self.nlegs = glob, fac, source, nlegs

    def __call__(self, n):
        n = np.asarray(n)
        if self.fac is None: return np.full(len(n), self.glob)
        out = np.full(len(n), self.glob)
        for (lo, hi), f in zip(zip(NBINS[:-1], NBINS[1:]), self.fac):
            out[(n >= lo) & (n < hi)] = f
        return out

    def jagged(self, n):
        """the same, for a jagged (per-jet) array of n"""
        return ak.unflatten(self(ak.to_numpy(ak.flatten(n))), ak.num(n))

    def short(self):
        if self.fac is None: return f"{self.glob:.2f}"
        return f"f(n_{{trk}}): {self.fac[0]:.2f} (1 trk) - {self.fac[-1]:.2f} (#geq {NBINS[-2]})"

    def md(self):
        if self.fac is None: return f"x {self.glob:.3f} ({self.source})"
        bins = ", ".join(f"n={lo}" + ("" if hi == lo + 1 else (f"-{hi - 1}" if hi < 10**6 else "+")) + f": {f:.2f}"
                         for (lo, hi), f in zip(zip(NBINS[:-1], NBINS[1:]), self.fac))
        return f"x f(n) from {self.source} ({self.nlegs:,} legs; {bins}; global {self.glob:.3f})"


def measure(e_all, binned, source):
    p, n = hs_pulls(e_all)
    g = mad(p)
    if not binned: return Calib(g, None, source, len(p))
    fac = []
    for lo, hi in zip(NBINS[:-1], NBINS[1:]):
        m = (n >= lo) & (n < hi)
        fac.append(mad(p[m]) if m.sum() >= 20 else g)
    return Calib(g, fac, source, len(p))


def eff(w, kept):
    """Weighted efficiency, its weighted-binomial error, and a Wilson interval.

    The error is sqrt(sum_kept w^2 (1-e)^2 + sum_rejected w^2 e^2) / sum w. It
    vanishes at e = 0 or 1, which in a column of a few Z+jets events is a lie,
    so the drawn interval is Wilson's (1 sigma) at the effective event count
    that reproduces that variance -- or, at the boundary, at (sum w)^2 / sum w^2.
    Returns (e, err, lo, hi) with lo/hi the distances below/above e."""
    W = w.sum()
    if W <= 0: return np.nan, np.nan, np.nan, np.nan
    e = w[kept].sum() / W
    var = ((w[kept] ** 2).sum() * (1 - e) ** 2 + (w[~kept] ** 2).sum() * e ** 2) / W ** 2
    n = e * (1 - e) / var if (0 < e < 1 and var > 0) else W ** 2 / (w ** 2).sum()
    ec = min(max(e, 0.0), 1.0)
    cen = (ec + 0.5 / n) / (1 + 1 / n)
    half = np.sqrt(ec * (1 - ec) / n + 0.25 / n ** 2) / (1 + 1 / n)
    return e, np.sqrt(var), max(e - (cen - half), 0.0), max(cen + half - e, 0.0)


def mjj_index(mjj, edges):
    """column of each event; overflow into the top bin, NaN (not kept) -> -1"""
    idx = np.digitize(np.minimum(mjj, 0.5 * (edges[-2] + edges[-1])), edges) - 1
    return np.where(np.isfinite(mjj), idx, -1)


# ── ROOT drawing ─────────────────────────────────────────────────────────────
ROOT = None
KEEP = []   # PyROOT objects must outlive the canvas print


def setup_root(out_stem):
    global ROOT
    import ROOT as R
    ROOT = R
    ROOT.gROOT.SetBatch(True)
    ROOT.gInterpreter.Declare(f'#include "{REPO}/src/AtlasStyle.h"\n#include "{REPO}/src/AtlasLabels.h"')
    ROOT.SetAtlasStyle()
    ROOT.gStyle.SetOptStat(0)
    ROOT.gErrorIgnoreLevel = ROOT.kWarning
    os.makedirs(os.path.dirname(os.path.abspath(out_stem)), exist_ok=True)
    return ROOT


def colour(c):
    return getattr(ROOT, c) if c.startswith("k") else ROOT.TColor.GetColor(c)


def graph(series, k, n_meth, key, edges, style):
    """markers dodged inside each m_jj column so the methods do not overlap"""
    g = ROOT.TGraphAsymmErrors()
    for b, (v, lo, hi, _) in enumerate(series):
        if not np.isfinite(v): continue
        x0, x1 = edges[b], edges[b + 1]
        x = x0 + (x1 - x0) * (0.5 + 0.55 * ((k + 0.5) / n_meth - 0.5))
        i = g.GetN()
        g.SetPoint(i, x, v)
        g.SetPointError(i, 0, 0, lo, hi)
    lab, c, mk = style[key]
    g.SetMarkerColor(colour(c)); g.SetLineColor(colour(c)); g.SetMarkerStyle(mk)
    g.SetMarkerSize(1.15 if key == "jj" else 1.0); g.SetLineWidth(2)
    KEEP.append(g)
    return g


def atlas_block(label, lines):
    ROOT.ATLASLabel(0.18, 0.895, "")
    al = ROOT.TLatex(); al.SetNDC(); al.SetTextFont(42); al.SetTextSize(0.045)
    al.DrawLatex(0.305, 0.895, "Simulation Internal")
    ROOT.ATLASEnergyLabel(0.18, 0.85, label)
    p = ROOT.TLatex(); p.SetNDC(); p.SetTextFont(42); p.SetTextSize(0.028); p.SetTextColor(ROOT.kGray + 2)
    for i, ln in enumerate(lines):
        p.DrawLatex(0.18, 0.81 - 0.036 * i, ln)
    KEEP.extend([al, p])
    return 0.81 - 0.036 * len(lines)


def y_range(seriess, methods, ref=1.0, pct=False):
    """(ymin, ytop) of column_figure's frame for one or more series dicts: pass the
    series of every figure that shares a slide, and they get the same axis"""
    lows = [v - lo for series in seriess for m in methods for v, lo, _, _ in series[m] if np.isfinite(v)]
    highs = [v + hi for series in seriess for m in methods for v, _, hi, _ in series[m] if np.isfinite(v)]
    ytop = max([max(highs) if highs else 1.0] + ([ref] if ref is not None else []))
    if pct:               # signed change in %: round down to the next 10%, the ref line in range
        ymin = np.floor((min(lows + [ref]) - 3.0) / 10) * 10 if lows else -10.0
    elif ref is not None: # efficiency-like: round down to the next 0.1
        ymin = max(0.0, np.floor((min(lows) - 0.03) * 10) / 10) if lows else 0.0
    else:                 # absolute: a margin below the lowest point
        ymin = max(0.0, min(lows) - 0.06 * (ytop - min(lows))) if lows else 0.0
    return ymin, ytop


def column_figure(stem, ylab, series, legend_extra, label, lines, counts, count_title, edges, style, methods,
                  ref=1.0, xnote="last bin includes overflow", pct=False, yrange=None):
    """The region plots' layout: per-m_jj-column values on top (one dodged marker
    per method), the MC event count of each column underneath on a log scale.
    `ref` draws a dashed reference line (None: none); `pct`: the values are a
    signed change in %, so the axis may go below zero; `yrange`: (ymin, ytop)
    from y_range() over several figures, to share their axis."""
    c = ROOT.TCanvas("c_" + os.path.basename(stem), "", 800, 780)
    pTop = ROOT.TPad("pTop", "", 0.0, 0.30, 1.0, 1.0)
    pBot = ROOT.TPad("pBot", "", 0.0, 0.00, 1.0, 0.30)
    pTop.SetBottomMargin(0.025); pTop.SetTopMargin(0.06)
    pBot.SetTopMargin(0.04); pBot.SetBottomMargin(0.40)
    pTop.Draw(); pBot.Draw(); KEEP.extend([c, pTop, pBot])
    pTop.cd()
    frame = ROOT.TH1D("frame_" + os.path.basename(stem), "", len(edges) - 1, edges); KEEP.append(frame)
    ymin, ytop = yrange if yrange is not None else y_range([series], methods, ref, pct)
    frame.SetMinimum(ymin); frame.SetMaximum(ytop + (ytop - ymin) * (0.95 + 0.13 * max(0, len(lines) - 3)))
    frame.GetXaxis().SetLabelSize(0); frame.GetYaxis().SetTitle(ylab); frame.GetYaxis().SetTitleOffset(1.25)
    frame.Draw("AXIS")
    if ref is not None:
        one = ROOT.TLine(edges[0], ref, edges[-1], ref); one.SetLineStyle(2); one.SetLineColor(ROOT.kGray + 1); one.Draw()
        KEEP.append(one)
    ytxt = atlas_block(label, lines)
    nrow = (len(methods) + 1) // 2
    leg = ROOT.TLegend(0.17, ytxt - 0.055 * nrow, 0.93, ytxt - 0.005); ROOT.StyleLegend(leg); leg.SetNColumns(2); KEEP.append(leg)
    leg.SetMargin(0.15)   # TLegend's 0.25 marker column runs the long right-hand entries into the frame
    for k, m in enumerate(methods):
        g = graph(series[m], k, len(methods), m, edges, style)
        g.Draw("P SAME")
        leg.AddEntry(g, style[m][0] + legend_extra.get(m, ""), "p")
    leg.Draw()
    pTop.RedrawAxis()
    pBot.cd(); pBot.SetLogy(True)
    tot = ROOT.TH1D("tot_" + os.path.basename(stem), "", len(edges) - 1, edges); KEEP.append(tot)
    for b, n in enumerate(counts): tot.SetBinContent(b + 1, n)
    tot.SetFillColor(ROOT.TColor.GetColor("#DCE6F2")); tot.SetLineColor(ROOT.kBlack); tot.SetLineWidth(1)
    tot.GetXaxis().SetTitle(f"m_{{jj}} [GeV]   ({xnote})")
    tot.GetYaxis().SetTitle(count_title)
    tot.GetXaxis().SetTitleSize(0.13); tot.GetXaxis().SetLabelSize(0.11); tot.GetXaxis().SetTitleOffset(1.25)
    tot.GetYaxis().SetTitleSize(0.11); tot.GetYaxis().SetLabelSize(0.10); tot.GetYaxis().SetTitleOffset(0.5)
    tot.GetYaxis().SetNdivisions(503)
    top = max(max(counts), 1)
    tot.SetMinimum(max(0.5, min([n for n in counts if n > 0] or [1]) / 3.0)); tot.SetMaximum(top * 60.0)
    tot.Draw("HIST")
    q = ROOT.TLatex(); q.SetTextFont(42); q.SetTextSize(0.075); KEEP.append(q)
    for b, n in enumerate(counts):
        if n <= 0: continue
        q.SetTextAlign(11 if b == 0 else 21)
        x = edges[0] + 0.015 * (edges[-1] - edges[0]) if b == 0 else 0.5 * (edges[b] + edges[b + 1])
        q.DrawLatex(x, top * (12.0 if b % 2 == 0 else 3.0), f"{int(n):,}")
    pBot.RedrawAxis()
    c.Print(stem + ".pdf"); c.Print(stem + ".png")
    print("wrote", stem + ".pdf / .png")


def roc_figure(stem, xy, methods, style, mark, lines, xlab, ylab, dashed=(), xmax=1.005):
    """One curve per method through its (x, y) points, the point at index `mark`
    enlarged; the curves fall from the top left to (1, 1), so the legend goes on
    the right under the label block. `lines`: up to three grey text lines."""
    c = ROOT.TCanvas("c_" + os.path.basename(stem), "", 800, 700); KEEP.append(c)
    c.SetLeftMargin(0.15); c.SetRightMargin(0.05)
    xs = np.concatenate([xy[m][0] for m in methods]); ys = np.concatenate([xy[m][1] for m in methods])
    fin = np.isfinite(xs) & np.isfinite(ys)   # a point with no background left has no rejection
    xs, ys = xs[fin], ys[fin]
    ylo = ys.min() - 0.05
    hi = xmax if xs.max() <= xmax else xs.max() + 0.01   # a yield ratio can pass 1
    fr = ROOT.TH1D("fr_" + os.path.basename(stem), "", 10, max(0.0, xs.min() - 0.03), hi)
    KEEP.append(fr)
    fr.SetMinimum(ylo); fr.SetMaximum(ys.max() + (ys.max() - ylo) * (0.9 + 0.08 * max(0, len(methods) - 5)))
    fr.GetXaxis().SetTitle(xlab)
    fr.GetYaxis().SetTitle(ylab)
    fr.Draw("AXIS")
    leg = ROOT.TLegend(0.50, 0.77 - 0.04 * len(methods), 0.93, 0.77); ROOT.StyleLegend(leg, 0.026); KEEP.append(leg)
    leg.SetMargin(0.15)
    for m in methods:
        lab, col, mk = style[m]
        x, y = np.asarray(xy[m][0], float), np.asarray(xy[m][1], float)
        ok = np.isfinite(x) & np.isfinite(y)
        g = ROOT.TGraph(int(ok.sum()), x[ok].copy(), y[ok].copy()); KEEP.append(g)
        g.SetLineColor(colour(col)); g.SetLineWidth(2); g.SetMarkerColor(colour(col)); g.SetMarkerStyle(mk); g.SetMarkerSize(0.8)
        if m in dashed: g.SetLineStyle(2)
        g.Draw("LP SAME")
        if ok[mark]:
            g3 = ROOT.TGraph(1, x[mark:mark + 1].copy(), y[mark:mark + 1].copy()); KEEP.append(g3)
            g3.SetMarkerColor(colour(col)); g3.SetMarkerStyle(mk); g3.SetMarkerSize(2.0); g3.Draw("P SAME")
        leg.AddEntry(g, lab, "lp")
    leg.Draw()
    ROOT.ATLASLabel(0.18, 0.895, "")
    t = ROOT.TLatex(); t.SetNDC(); t.SetTextFont(42); t.SetTextSize(0.04); KEEP.append(t)
    t.DrawLatex(0.315, 0.895, "Simulation Internal")
    t2 = ROOT.TLatex(); t2.SetNDC(); t2.SetTextFont(42); t2.SetTextSize(0.026); t2.SetTextColor(ROOT.kGray + 2); KEEP.append(t2)
    for i, ln in enumerate(lines):
        t2.DrawLatex(0.18, 0.855 - 0.03 * i, ln)
    c.Print(stem + ".pdf"); c.Print(stem + ".png")
    print("wrote", stem + ".pdf / .png")
