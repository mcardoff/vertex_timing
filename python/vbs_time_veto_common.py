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
    "jj":        ("Jet time compatibility (no t_{0})",         "#000000",    20),
    "jj_all":    ("Jet time compatibility, all ghost tracks",  "#8C8C8C",    24),
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


# HGTD acceptance of the tagging pair, from the legs' TRACKS, not the jet axes: an
# axis past |eta| = 4.0 can still put tracks in HGTD. vbs_time_veto counts each
# leg's ghost tracks passing the track selection in 2.38 < |eta| < 4.0 -- with a
# valid HGTD time (hgtd_ntimed, --acc-count timed, the default: the jet can be
# tested) or without (hgtd_n, --acc-count any). Reco-only. The category counts
# the legs with >= min_trk such tracks: "both" (2), "one" (1), "neither" (0);
# "all" is no restriction.
ACCEPTANCE = {"all": None, "both": 2, "one": 1, "neither": 0}
ACC_TEXT = {"both": "both tagging jets", "one": "exactly one tagging jet", "neither": "neither tagging jet"}
# the talk's names for the categories (vbs_time_veto_plot.py --compact headline)
CASE_TEXT = {"both": "Forward-Forward: both tagging jets in HGTD",
             "one": "Central-Forward: one tagging jet in HGTD",
             "neither": "neither tagging jet in HGTD", "all": "All selected events"}


def add_acceptance_args(ap):
    ap.add_argument("--acceptance", default="all", choices=list(ACCEPTANCE),
                    help="restrict to events whose (no-timing) tagging pair has this many legs in HGTD "
                         "acceptance, by tracks: >= --acc-min-trk ghost tracks in 2.38 < |eta| < 4.0, timed or not")
    ap.add_argument("--acc-min-trk", type=int, default=1, help="HGTD-acceptance tracks a leg needs (--acceptance)")
    ap.add_argument("--acc-count", default="timed", choices=["timed", "any"],
                    help="--acceptance counts the leg's HGTD tracks with a valid time (default) or all of them")


def acceptance_mask(nA, nB, mode, min_trk):
    """events whose pair has ACCEPTANCE[mode] legs with >= min_trk HGTD tracks"""
    want = ACCEPTANCE[mode]
    if want is None: return np.ones(len(nA), bool)
    return ((np.asarray(nA) >= min_trk).astype(int) + (np.asarray(nB) >= min_trk).astype(int)) == want


def acc_column(count):
    """the per-leg column --acc-count reads (prefix a_ / b_)"""
    return "hgtd_ntimed" if count == "timed" else "hgtd_n"


def acceptance_line(mode, min_trk, count="timed"):
    """the figures' label line; '' for no restriction"""
    if ACCEPTANCE[mode] is None: return ""
    what = "timed track" if count == "timed" else "track"
    trk = f"a {what}" if min_trk == 1 else f"#geq {min_trk} {what}s"
    return f"{ACC_TEXT[mode]} with {trk} in HGTD (2.38 < |#eta_{{trk}}| < 4.0)"


def mad(x):
    return 1.4826 * float(np.median(np.abs(x)))


_DG_N = [0]


def dgaus_core(values, lo=-400.0, hi=400.0, nbins=400):
    """the project's core resolution: the narrower width of a double-Gaussian fit, the exact
    recipe of FIT_DBLGAUS in src/plotting_utilities.h (shared mean within +-50 ps, sigma1
    seeded at 13 ps, both widths >= 5 ps, fit option RQI over the histogram range; the smaller
    sigma is the core, as fillFromModel's swapBySize). Entries outside [lo, hi] are not
    filled. Returns (sigma_core, its error)."""
    R = ROOT
    _DG_N[0] += 1
    h = R.TH1D(f"h_dgfit_{_DG_N[0]}", "", nbins, lo, hi)
    for x in np.asarray(values, float): h.Fill(x)
    f = R.TF1(f"dgaus_fit_{_DG_N[0]}", "[1]*TMath::Exp(-0.5*((x-[0])/[3])^2) + [2]*TMath::Exp(-0.5*((x-[0])/[4])^2)", lo, hi)
    f.SetParameters(h.GetMean(), 0.8 * h.GetMaximum(), 1e-2 * h.GetMaximum(), 13.0, h.GetStdDev())
    f.SetParLimits(0, -50.0, 50.0); f.SetParLimits(1, 0, 1e9); f.SetParLimits(2, 0, 1e9)
    f.SetParLimits(3, 5.0, 1e6); f.SetParLimits(4, 5.0, 1e6)
    h.Fit(f, "RQIN0")
    s1, s2 = abs(f.GetParameter(3)), abs(f.GetParameter(4))
    k = 3 if s1 <= s2 else 4
    KEEP.extend([h, f])
    return abs(f.GetParameter(k)), f.GetParError(k)


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


def legs(path):
    """the leg columns hs_pulls reads, from a vbs_time_veto output"""
    import uproot
    cols = ["t0_truth"] + [f"{L}_{k}" for L in "ab" for k in ["hs", "core_n", "core_t", "core_sig"]]
    return uproot.open(path)["events"].arrays(cols, library="np")


def measure(e_all, binned, source):
    p, n = hs_pulls(e_all)
    g = mad(p)
    if not binned: return Calib(g, None, source, len(p))
    fac = []
    for lo, hi in zip(NBINS[:-1], NBINS[1:]):
        m = (n >= lo) & (n < hi)
        fac.append(mad(p[m]) if m.sum() >= 20 else g)
    return Calib(g, fac, source, len(p))


def t0_veto(e, cal, src, min_trk, nsigma=3.0):
    """the jet-vs-t0 event veto exactly as vbs_time_veto_plot.decisions(): a leg is tested when its
    jet time is built from >= min_trk tracks and the t0 exists; it fails when
    |t_jet - t0| / sqrt((f(n) sigma_jet)^2 + (infl sigma_t0)^2) >= nsigma.
    Returns (vetoed, failA, failB) per event."""
    ok, t0, s0 = e[f"ok_{src}"].astype(bool), e[f"t0_{src}"], e[f"infl_{src}"] * e[f"sig0_{src}"]
    fails = []
    for L in "ab":
        sg = cal(e[f"{L}_core_n"]) * e[f"{L}_core_sig"]
        with np.errstate(divide="ignore", invalid="ignore"):
            pull = np.abs(e[f"{L}_core_t"] - t0) / np.sqrt(sg ** 2 + s0 ** 2)
        fails.append(ok & (e[f"{L}_core_n"] >= min_trk) & (pull >= nsigma))
    return fails[0] | fails[1], fails[0], fails[1]


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
    # stars (29/30) and diamonds (33/27) draw small at the common size
    g.SetMarkerSize(1.15 if key == "jj" else 1.7 if mk in (29, 30, 33, 27) else 1.0); g.SetLineWidth(2)
    KEEP.append(g)
    return g


def atlas_block(label, lines):
    ROOT.ATLASLabel(0.18, 0.895, "")
    al = ROOT.TLatex(); al.SetNDC(); al.SetTextFont(42); al.SetTextSize(0.045)
    al.DrawLatex(0.305, 0.895, "Simulation Internal")
    ROOT.ATLASEnergyLabel(0.18, 0.85, label)
    p = ROOT.TLatex(); p.SetNDC(); p.SetTextFont(42); p.SetTextSize(0.028); p.SetTextColor(ROOT.kGray + 2)
    h = ROOT.TLatex(); h.SetNDC(); h.SetTextFont(62); h.SetTextSize(0.036)   # a line starting "**": a headline
    y = 0.81
    for ln in lines:
        if ln.startswith("**"):
            h.DrawLatex(0.18, y - 0.004, ln[2:]); y -= 0.046
        else:
            p.DrawLatex(0.18, y, ln); y -= 0.036
    KEEP.extend([al, p, h])
    return y


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


def _count_panel(pad, stem, counts, count_title, count_parts, edges, xnote, xaxis, ts, legend, labels=None):
    """the MC count per m_jj column in `pad`: on a log scale, or -- with count_parts,
    [(label, colour, per-column counts)] -- a linear stack, in 10^3 when large, the
    total's outline on top. `ts` scales the text to the pad's height; `xaxis` draws
    the m_jj axis; `legend` the parts' legend"""
    pad.cd(); pad.SetLogy(True)
    name = os.path.basename(stem)
    tot = ROOT.TH1D("tot_" + name, "", len(edges) - 1, edges); KEEP.append(tot)
    for b, n in enumerate(counts): tot.SetBinContent(b + 1, n)
    tot.SetFillColor(ROOT.TColor.GetColor("#DCE6F2")); tot.SetLineColor(ROOT.kBlack); tot.SetLineWidth(1)
    tot.GetXaxis().SetTitle(f"m_{{jj}} [GeV]   ({xnote})" if xaxis else "")
    tot.GetYaxis().SetTitle(count_title)
    tot.GetXaxis().SetTitleSize(0.13 * ts); tot.GetXaxis().SetLabelSize(0.11 * ts if xaxis else 0)
    tot.GetXaxis().SetTitleOffset(1.25); tot.GetXaxis().SetTickLength(0.03 * ts)
    tot.GetYaxis().SetTitleSize(0.11 * ts); tot.GetYaxis().SetLabelSize(0.10 * ts); tot.GetYaxis().SetTitleOffset(0.5 if ts == 1.0 else 0.42)
    tot.GetYaxis().SetNdivisions(503)
    top = max(max(counts), 1)
    sc = 1.0
    tot.SetMinimum(max(0.5, min([n for n in counts if n > 0] or [1]) / 3.0)); tot.SetMaximum(top * 60.0)
    if count_parts:   # the count split by case: a stack (linear, so the shares read), the total's outline on top
        pad.SetLogy(False)
        sc = 1e-3 if top >= 10000 else 1.0   # thousands: the stack in units of 10^3 (the printed counts stay exact)
        if sc != 1.0:
            for b in range(len(counts)): tot.SetBinContent(b + 1, counts[b] * sc)
            tot.GetYaxis().SetTitle(count_title + " [10^{3}]")
        tot.SetMinimum(0.0); tot.SetMaximum(top * sc * 2.45)
        ROOT.TGaxis.SetMaxDigits(6)   # full numbers on the axis: the x10^3 exponent sits on the pad edge
        tot.GetYaxis().SetNdivisions(203 if ts == 1.0 else 202)
        tot.SetFillStyle(0)
        tot.Draw("HIST")
        stk = ROOT.THStack("stk_" + name, ""); KEEP.append(stk)
        lg = ROOT.TLegend(0.16, 1 - pad.GetTopMargin() - 0.14 * ts, 0.97, 1 - pad.GetTopMargin() - 0.005)
        lg.SetNColumns(len(count_parts) + (1 if labels else 0)); KEEP.append(lg)
        lg.SetBorderSize(0); lg.SetFillStyle(0); lg.SetTextSize(0.075 * ts); lg.SetTextFont(42)
        for k, (lab, col, vals) in enumerate(count_parts):
            h = ROOT.TH1D(f"cp_{k}_" + name, "", len(edges) - 1, edges); KEEP.append(h)
            for b, n in enumerate(vals): h.SetBinContent(b + 1, n * sc)
            h.SetFillColor(ROOT.TColor.GetColor(col)); h.SetLineColor(ROOT.kBlack); h.SetLineWidth(1)   # as the region stacks
            stk.Add(h)
            if sum(vals) > 0: lg.AddEntry(h, lab, "f")
        if labels: lg.AddEntry(ROOT.nullptr, "numbers: S / B", "")
        stk.Draw("HIST SAME")
        tot.Draw("HIST SAME")
        if legend: lg.Draw()
    else:
        tot.Draw("HIST")
    q = ROOT.TLatex(); q.SetTextFont(42); q.SetTextSize((0.064 if labels else 0.075) * ts); KEEP.append(q)
    for b, n in enumerate(counts):
        if n <= 0: continue
        q.SetTextAlign(11 if b == 0 else 21)
        x = edges[0] + 0.015 * (edges[-1] - edges[0]) if b == 0 else 0.5 * (edges[b] + edges[b + 1])
        q.DrawLatex(x, top * (sc * (1.34 if b % 2 == 0 else 1.06) if count_parts else (12.0 if b % 2 == 0 else 3.0)),
                    labels[b] if labels else f"{int(n):,}")
    pad.RedrawAxis()


def column_figure(stem, ylab, series, legend_extra, label, lines, counts, count_title, edges, style, methods,
                  ref=1.0, xnote="last bin includes overflow", pct=False, yrange=None, compact=False, count_parts=None,
                  sig_counts=None):
    """The region plots' layout: per-m_jj-column values on top (one dodged marker
    per method), the MC event count of each column underneath on a log scale.
    `ref` draws a dashed reference line (None: none); `pct`: the values are a
    signed change in %, so the axis may go below zero; `yrange`: (ymin, ytop)
    from y_range() over several figures, to share their axis."""
    c = ROOT.TCanvas("c_" + os.path.basename(stem), "", 800, 780)
    pTop = ROOT.TPad("pTop", "", 0.0, 0.30, 1.0, 1.0)
    pBot = ROOT.TPad("pBot", "", 0.0, 0.00, 1.0, 0.30)
    pBot.SetTopMargin(0.04); pBot.SetBottomMargin(0.40)
    pTop.SetBottomMargin(0.025); pTop.SetTopMargin(0.06)
    pTop.Draw(); pBot.Draw(); KEEP.extend([c, pTop, pBot])
    pTop.cd()
    frame = ROOT.TH1D("frame_" + os.path.basename(stem), "", len(edges) - 1, edges); KEEP.append(frame)
    ymin, ytop = yrange if yrange is not None else y_range([series], methods, ref, pct)
    # headroom for the label block and legend; at least enough that the legend's
    # bottom edge (pad NDC; frame from 0.025 to 0.94) stays above the highest point
    legbot = 0.81 - 0.036 * len(lines) - 0.055 * ((len(methods) + 1) // 2)
    top = ymin + 0.915 * (ytop - ymin) / (legbot - 0.04)   # the legend's bottom edge clears the highest point
    if not compact:   # the long label block needs room of its own
        top = max(top, ytop + (ytop - ymin) * (0.95 + 0.13 * max(0, len(lines) - 3)))
    frame.SetMinimum(ymin); frame.SetMaximum(top)
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
    if sig_counts is None:
        _count_panel(pBot, stem, counts, count_title, count_parts, edges, xnote, True, 1.0, True)
    else:   # S/B-type figures: signal + background MC in ONE panel, case-stacked; labels "S / B"
        scounts, sparts = sig_counts
        tot = [a + b for a, b in zip(scounts, counts)]
        parts = None
        if count_parts and sparts:
            parts = [(lab, col, [x + y for x, y in zip(sv, bv)]) for (lab, col, sv), (_, _, bv) in zip(sparts, count_parts)]
        _count_panel(pBot, stem, tot, "S + B MC", parts, edges, xnote, True, 1.0, True,
                     labels=[(f"{a / 1000:.1f}k" if a >= 10000 else f"{a:,}") + f" / {b:,}" for a, b in zip(scounts, counts)])
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
