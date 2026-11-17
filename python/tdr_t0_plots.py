#!/usr/bin/env python3
"""The vertex t0 figures of the TDR-structured presentation (util/tdr_t0's tree and fits), ATLAS style.

  <out>_resid     t0 - t_HS for TRKPTZ (left) and WAVeS (right), from the ntuple's track times (dashed)
                  and with the TDR's track-time cleaning (solid), Athena's HGTD t0 (grey, cannot be
                  re-cleaned) as the reference; normalised to the selected events, so a missing t0
                  lowers the area. Lower panels: cleaned / as is in 20 ps bins.
  <out>_tdr       the TDR Fig 3.19 analogue: t0 - t_HS for the presentation's t0s -- HGTD (Athena), and
                  TRKPTZ and WAVeS from cleaned track times (the default everywhere) -- one panel, the
                  legend giving each one's core fraction and FIT_PEAKTAIL core width
  <out>_core_nhs  core fraction (|t0 - t_HS| < 60 ps over all events, a missing t0 a failure) against
                  the number of forward hard-scatter tracks with a time (folded at 20, as the main
                  analysis); lower panel: cleaned / as is (errors from the events that move, paired).

  <out>_waves_resid     clustering_plot's inclusive-resolution style, WAVeS only (the user's choice for the
                        slides, 2026-11-04): t0 - t_HS with cleaned track times (filled, with its peak + tail
                        fit) against the ntuple's times (dashed), log scale, core fraction and width in the legend
  <out>_waves_core_nhs  clustering_plot's core-fraction style, WAVeS only: points with binomial errors and no
                        connecting lines, cleaned (filled) vs as is (open), over the n-forward-HS-track
                        distribution (grey, arbitrary scale); no ratio panel

    cd python && PYTHONNOUSERSITE=1 ~/.venv-hgtd/bin/python tdr_t0_plots.py ../figs/hists/tdr_t0.root \\
        --label "VBF H#rightarrowinv. (local sample)" --out ../figs/tdr/local_tdr_t0
"""
import argparse
import numpy as np
import uproot
import tdr_root_style as S

ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
ap.add_argument("tree", help="tdr_t0 output (its fits are read from <stem>_fits.root)")
ap.add_argument("--label", default="VBF H#rightarrowinv.", help="sample, ROOT TLatex")
ap.add_argument("--out", required=True)
ap.add_argument("--t0s", default=None, help="comma list of the t0s to draw (hgtd, trkptz, waves); default all. "
                "The Global Performance talk shows WAVeS against Athena only (user, 2026-11-11): --t0s hgtd,waves,truth")
args = ap.parse_args()
SHOW = lambda name: not args.t0s or name.split("_")[0] in args.t0s.split(",")
R = S.setup()
PASS = 60.0
FOLD = 20
ENERGY = "#sqrt{s} = 14 TeV, #LT#mu#GT = 200"

E = uproot.open(args.tree)["events"].arrays(library="np")
N = len(E["n_hs_trk"])
fits = R.TFile.Open(args.tree.replace(".root", "_fits.root"))
sig = lambda n: fits.Get(f"sigcore_{n}").GetVal()
core = {n: np.abs(E[f"dt_{n}"]) < PASS for n in ("hgtd", "trkptz_raw", "trkptz_clean", "waves_raw", "waves_clean")}   # NaN -> False
frac = {n: core[n].mean() for n in core}

print(f"{args.tree}: {N} events (main-analysis gate)\n  {'t0':14s} {'has a t0':>9s} {'core fraction':>14s} {'sigma_core':>11s}")
for n in core:
    print(f"  {n:14s} {100 * np.isfinite(E['dt_' + n]).mean():8.1f}% {100 * frac[n]:13.2f}% {sig(n):8.1f} ps")
for s in ("trkptz", "waves"):
    a, b = core[f"{s}_raw"], core[f"{s}_clean"]
    gain, loss = np.sum(~a & b), np.sum(a & ~b)
    print(f"  {s}: cleaning moves {gain} events into the core and {loss} out: {100 * (gain - loss) / N:+.2f} +- {100 * np.sqrt(gain + loss) / N:.2f} points (paired)")

T0S = [t for t in [("trkptz", "TRKPTZ t_{0}", S.COL["trkptz"]), ("waves", "WAVeS t_{0}", S.COL["waves"])] if SHOW(t[0])]
MAIN, RAT = 0.68, 0.32
TS_M, TS_R = 0.05 * 600 / (800 * MAIN), 0.05 * 600 / (800 * RAT)

# ── residuals ─────────────────────────────────────────────────────────────────
cv = S.keep(R.TCanvas("resid", "", 1600, 800))
for c, (key, title, colour) in enumerate(T0S):
    cv.cd(); S.pad(f"m{c}", 0.5 * c, RAT, 0.5 * c + 0.5, 1.0, top=0.05 / MAIN * 0.6, bottom=0.02)
    R.gPad.SetLogy()
    hs = {}
    for name, col, ls, w in [("hgtd", S.COL["hgtd"], 1, 2), (f"{key}_raw", colour, 2, 3), (f"{key}_clean", colour, 1, 3)]:
        h = S.keep(fits.Get(f"h_{name}").Clone(f"c_{name}_{c}"))
        h.GetListOfFunctions().Clear()
        h.Scale(1.0 / N)
        S.style(h, col, 0, ls, w)
        hs[name] = h
    fr = hs["hgtd"]
    fr.SetMinimum(2e-5); fr.SetMaximum(300.0)
    S.axes(fr, "", "Fraction of events / 4 ps", size=TS_M, ytoff=1.35)
    fr.GetXaxis().SetLabelSize(0); fr.GetXaxis().SetNdivisions(505)
    fr.Draw("HIST")
    hs[f"{key}_raw"].Draw("HIST SAME"); hs[f"{key}_clean"].Draw("HIST SAME")
    leg = S.legend(0.20, 0.60, 0.93, 0.78, size=0.75 * TS_M, margin=0.12)
    labels = {"hgtd": "HGTD (Athena) t_{0}", f"{key}_raw": f"{title}, as is", f"{key}_clean": f"{title}, cleaned (TDR)"}
    for name, h in hs.items():
        leg.AddEntry(h, f"{labels[name]}: {100 * frac[name]:.2f}%, #sigma_{{core}} {sig(name):.1f} ps", "L")
    leg.Draw()
    S.header(0.20, 0.87, [f"{ENERGY}, {args.label}"], size=TS_M)
    S.text(0.20, 0.06, "Core fraction: |t_{0} #minus t_{HS}| < 60 ps, all events", 0.75 * TS_M)
    R.gPad.RedrawAxis()
    cv.cd(); S.pad(f"r{c}", 0.5 * c, 0.0, 0.5 * c + 0.5, RAT, top=0.03, bottom=0.36)
    num = S.keep(hs[f"{key}_clean"].Clone(f"num{c}")); den = S.keep(hs[f"{key}_raw"].Clone(f"den{c}"))
    num.Rebin(5); den.Rebin(5); num.Divide(den)
    num.SetLineStyle(1); num.SetLineWidth(2); num.SetMarkerStyle(20); num.SetMarkerSize(0.8)
    num.SetMarkerColor(S.colour(colour)); num.SetLineColor(S.colour(colour))
    num.SetMinimum(0.6); num.SetMaximum(1.4)
    S.axes(num, "t_{0} #minus t_{HS} [ps]", "Cleaned / as is", size=TS_R, ytoff=1.35 * TS_M / TS_R, xtoff=1.2, ndiv_y=504)
    num.GetXaxis().SetNdivisions(505); num.GetYaxis().ChangeLabel(1, -1, 0)
    num.Draw("P")
    one = S.keep(R.TLine(-400, 1, 400, 1)); one.SetLineStyle(2); one.SetLineColor(R.kGray + 2); one.Draw()
    S.text(0.20, 0.82, "20 ps bins", 0.75 * TS_R)
S.save(cv, f"{args.out}_resid")

# ── the TDR Fig 3.19 analogue: the cleaned t0s with Athena's ──────────────────
cv = S.keep(R.TCanvas("tdr319", "", 800, 600))
R.gPad.SetLogy()
leg = S.legend(0.20, 0.57, 0.80, 0.78, size=0.036, margin=0.12)
leg.SetHeader("Core fraction (|t_{0} #minus t_{HS}| < 60 ps), core width:")
first = True
for name, lab, col, w in [("hgtd", "HGTD (Athena) t_{0}", S.COL["hgtd"], 3), ("trkptz_clean", "TRKPTZ t_{0}", S.COL["trkptz"], 3),
                          ("waves_clean", "WAVeS t_{0}", S.COL["waves"], 3)]:
    if not SHOW(name): continue
    h = S.keep(fits.Get(f"h_{name}").Clone(f"t319_{name}"))
    h.GetListOfFunctions().Clear()
    h.Scale(1.0 / N)
    S.style(h, col, 0, 1, w)
    if first:
        h.SetMinimum(2e-5); h.SetMaximum(3000.0)   # headroom: the legend sits above the peak
        S.axes(h, "t_{0} #minus t_{HS} [ps]", "Fraction of events / 4 ps")
        h.GetXaxis().SetNdivisions(505)
        h.Draw("HIST"); first = False
    else:
        h.Draw("HIST SAME")
    leg.AddEntry(h, f"{lab}: {100 * frac[name]:.1f}%, #sigma_{{core}} {sig(name):.1f} ps", "L")
leg.Draw()
S.header(0.20, 0.87, [f"{ENERGY}, {args.label}"])
R.gPad.RedrawAxis()
S.save(cv, f"{args.out}_tdr")

# ── core fraction against forward hard-scatter tracks ─────────────────────────
x = np.minimum(E["n_hs_trk"], FOLD)
bins = np.arange(0, FOLD + 1)
def per_bin(m):
    f, e = [], []
    for b in bins:
        sel = x == b; n = sel.sum()
        p = m[sel].mean() if n else np.nan
        f.append(p); e.append(np.sqrt(p * (1 - p) / n) if n else np.nan)
    return np.array(f), np.array(e)

cv = S.keep(R.TCanvas("corenhs", "", 800, 800))
cv.cd(); S.pad("m", 0.0, RAT, 1.0, 1.0, top=0.05 / MAIN * 0.6, bottom=0.02)
fr = S.keep(R.gPad.DrawFrame(-0.5, 0.0, FOLD + 0.5, 1.45))
S.axes(fr, "", "Core fraction", size=TS_M, ytoff=1.35)
fr.GetXaxis().SetLabelSize(0); fr.GetYaxis().ChangeLabel(1, -1, 0)
leg = S.legend(0.45, 0.30, 0.93, 0.60, size=0.75 * TS_M, margin=0.18)
for name, lab, col, mk, ls in [("hgtd", "HGTD (Athena) t_{0}", S.COL["hgtd"], 23, 1),
                               ("trkptz_raw", "TRKPTZ t_{0}, as is", S.COL["trkptz"], 25, 2), ("trkptz_clean", "TRKPTZ t_{0}, cleaned", S.COL["trkptz"], 21, 1),
                               ("waves_raw", "WAVeS t_{0}, as is", S.COL["waves"], 24, 2), ("waves_clean", "WAVeS t_{0}, cleaned", S.COL["waves"], 20, 1)]:
    if not SHOW(name): continue
    f, e = per_bin(core[name])
    g = S.style(S.graph(bins, f, e), col, mk, ls, 2, 1.1); g.Draw("LP SAME")
    leg.AddEntry(g, lab, "LP")
leg.Draw()
S.header(0.20, 0.87, [f"{ENERGY}, {args.label}", "Core fraction: |t_{0} #minus t_{HS}| < 60 ps, all events"], size=TS_M)
cv.cd(); S.pad("r", 0.0, 0.0, 1.0, RAT, top=0.03, bottom=0.36)
fr2 = S.keep(R.gPad.DrawFrame(-0.5, 0.92, FOLD + 0.5, 1.12))
S.axes(fr2, f"Forward hard-scatter tracks with a time (last bin: #geq {FOLD})", "Cleaned / as is", size=TS_R, ytoff=1.35 * TS_M / TS_R, xtoff=1.2, ndiv_y=504)
fr2.GetYaxis().ChangeLabel(1, -1, 0)
one = S.keep(R.TLine(-0.5, 1, FOLD + 0.5, 1)); one.SetLineStyle(2); one.SetLineColor(R.kGray + 2); one.Draw()
print(f"\ncore fraction, cleaned / as is, by forward HS tracks (0..{FOLD - 1}, >= {FOLD}); events per bin:")
print("  " + " ".join(f"{np.sum(x == b):6d}" for b in bins))
for k, (key, title, colour) in enumerate(T0S):
    a, b_ = core[f"{key}_raw"], core[f"{key}_clean"]
    d, de = [], []
    for b in bins:
        sel = x == b; n = sel.sum()
        gain, loss, base = np.sum(sel & ~a & b_), np.sum(sel & a & ~b_), np.sum(sel & a)
        # paired: cleaned / as is = 1 + (gain - loss) / base, the error from the events that moved
        d.append(1 + (gain - loss) / base if base else np.nan); de.append(np.sqrt(gain + loss) / base if base else np.nan)
    g = S.style(S.graph(bins + (k - 0.5) * 0.25, d, de), colour, 21 if key == "trkptz" else 20, 1, 2, 1.0); g.Draw("P SAME")
    print(f"  {key:7s} " + " ".join(f"{v:6.3f}" for v in d))
S.save(cv, f"{args.out}_core_nhs")


# ── clustering_plot style, WAVeS only (the slides' section 1) ─────────────────
# Each with a ratio pad, cleaned / as in the ntuple (the user's request, 2026-11-04).
LIGHT = R.TColor.GetColor("#FFD891")                     # a light WAVeS fill
GREY = R.TColor.GetColor("#4D4D4D")
MW, RW = 0.70, 0.30                                      # main / ratio pad heights
TSW, TSR_ = 0.052, 0.052 * MW / RW


def two_pads(name):
    cv_ = S.keep(R.TCanvas(name, "", 800, 860))
    cv_.cd(); m = S.pad(name + "_m", 0.0, RW, 1.0, 1.0, top=0.04, bottom=0.02, left=0.15, right=0.04)
    cv_.cd(); r = S.pad(name + "_r", 0.0, 0.0, 1.0, RW, top=0.03, bottom=0.36, left=0.15, right=0.04)
    return cv_, m, r


def ratio_frame(x0, x1, lo, hi, xtitle, ndiv_y=504):
    """a fixed range, the same on every sample, chosen so the labels sit clear of the frame's edges"""
    f = S.keep(R.gPad.DrawFrame(x0, lo, x1, hi))
    S.axes(f, xtitle, "Cleaned / as is", size=TSR_, ytoff=1.4 * TSW / TSR_, xtoff=1.15, ndiv_y=ndiv_y)
    one = S.keep(R.TLine(x0, 1.0, x1, 1.0)); one.SetLineStyle(2); one.SetLineColor(R.kGray + 2); one.Draw()
    return f


# t0 - t_HS, WAVeS: cleaned (filled, its fit) vs ntuple times (dashed); ratio in 20 ps bins
cvw, pm, pr = two_pads("waves_resid")
pm.cd(); pm.SetLogy()
hc = S.keep(fits.Get("h_waves_clean").Clone("w_clean")); hr = S.keep(fits.Get("h_waves_raw").Clone("w_raw"))
fc = hc.GetListOfFunctions().FindObject("peaktail_fit")
fc = S.keep(fc.Clone("w_fit")) if fc else None
hc.GetListOfFunctions().Clear(); hr.GetListOfFunctions().Clear()
peak = max(hc.GetMaximum(), hr.GetMaximum())
frw = S.keep(pm.DrawFrame(-400.0, 0.5, 400.0, 100.0 * peak))
S.axes(frw, "", f"Events / {hc.GetBinWidth(1):.0f} ps", size=TSW, ytoff=1.4)
frw.GetXaxis().SetLabelSize(0)
hc.SetFillColor(LIGHT); hc.SetLineColor(S.colour(S.COL["waves"])); hc.SetLineWidth(2); hc.Draw("HIST SAME")
hr.SetFillStyle(0); hr.SetLineColor(GREY); hr.SetLineWidth(2); hr.SetLineStyle(2); hr.Draw("HIST SAME")
if fc:
    fc.SetLineColor(S.colour(S.COL["waves"])); fc.SetLineWidth(2); fc.SetNpx(800); fc.Draw("SAME")
legw = S.legend(0.66, 0.73, 0.95, 0.91, size=0.040, margin=0.25)   # top right: clear of the narrow peak
legw.AddEntry(hc, "Cleaned", "F")
legw.AddEntry(hr, "As in the ntuple", "L")
if fc: legw.AddEntry(fc, "Fit (cleaned)", "L")
legw.Draw()
S.header(0.19, 0.92, [f"{ENERGY}, {args.label}", "WAVeS t_{0}, track times:", "core fraction (< 60 ps), core width",
                      f"cleaned: {100 * frac['waves_clean']:.1f}%, {sig('waves_clean'):.1f} ps",
                      f"as in the ntuple: {100 * frac['waves_raw']:.1f}%, {sig('waves_raw'):.1f} ps"], size=0.040)
pm.RedrawAxis()
pr.cd()
REB = 5                                                  # 4 ps -> 20 ps bins for the ratio
rc, rr = S.keep(hc.Clone("w_clean_r")), S.keep(hr.Clone("w_raw_r"))
rc.Rebin(REB); rr.Rebin(REB)
xs, ys, es = [], [], []
for b in range(1, rc.GetNbinsX() + 1):
    a, c_ = rr.GetBinContent(b), rc.GetBinContent(b)
    if a > 0 and c_ > 0:
        xs.append(rc.GetBinCenter(b)); ys.append(c_ / a)
        es.append(c_ / a * np.sqrt(1.0 / c_ + 1.0 / a))  # independent Poisson: conservative (the same events)
ratio_frame(-400.0, 400.0, 0.7, 1.3, "t_{0} #minus t_{HS} [ps]", ndiv_y=503)
S.style(S.graph(np.array(xs), np.array(ys), np.array(es)), S.COL["waves"], 20, 1, 2, 0.9).Draw("P SAME")
S.save(cvw, f"{args.out}_waves_resid")

# core fraction vs forward HS tracks, WAVeS: cleaned vs ntuple times; paired ratio
cvn, pm, pr = two_pads("waves_core_nhs")
pm.cd()
frn = S.keep(pm.DrawFrame(-0.5, 0.0, FOLD + 0.5, 1.75))
S.axes(frn, "", "Core fraction (|t_{0} #minus t_{HS}| < 60 ps)", size=TSW, ytoff=1.4)
frn.GetXaxis().SetLabelSize(0)
counts = np.array([np.sum(x == k) for k in range(FOLD + 1)], float)
hn = S.keep(R.TH1D("nhs_dist", "", FOLD + 1, -0.5, FOLD + 0.5))
for k, c_ in enumerate(counts): hn.SetBinContent(k + 1, 0.5 * c_ / counts.max())   # arbitrary scale: tallest bin at 0.5
hn.SetFillColor(R.TColor.GetColor("#BDBDC4")); hn.SetLineColor(R.TColor.GetColor("#BDBDC4")); hn.Draw("HIST SAME")
legn = S.legend(0.42, 0.62, 0.93, 0.78, size=0.040, margin=0.2)   # above the points (<= 1), below the header
for name, txt, mk in (("waves_clean", "WAVeS t_{0}, cleaned", 20), ("waves_raw", "WAVeS t_{0}, as in the ntuple", 24)):
    f, e = per_bin(core[name])
    g = S.style(S.graph(np.arange(FOLD + 1, dtype=float), np.array(f), np.array(e)), S.COL["waves"] if mk == 20 else "#4D4D4D", mk, 1, 2, 1.3)
    g.Draw("P SAME")                                       # points and error bars only: no connecting lines
    legn.AddEntry(g, txt, "PE")
legn.AddEntry(hn, "Events (arb. scale)", "F")
legn.Draw()
S.header(0.19, 0.92, [f"{ENERGY}, {args.label}", "Core fraction over all events (no t_{0}: a failure)"], size=0.040)
pm.RedrawAxis()
pr.cd()
a, b_ = core["waves_raw"], core["waves_clean"]
d, de = [], []
for b in bins:
    sel = x == b
    gain, loss, base = np.sum(sel & ~a & b_), np.sum(sel & a & ~b_), np.sum(sel & a)
    # paired, as <out>_core_nhs: cleaned / as is = 1 + (gain - loss) / base, the error from the events that moved
    d.append(1 + (gain - loss) / base if base else np.nan); de.append(np.sqrt(gain + loss) / base if base else np.nan)
ratio_frame(-0.5, FOLD + 0.5, 0.93, 1.10, f"Forward hard-scatter tracks with a time (last bin: #geq {FOLD})")
S.style(S.graph(np.arange(FOLD + 1, dtype=float), np.array(d), np.array(de)), S.COL["waves"], 20, 1, 2, 1.1).Draw("P SAME")
S.save(cvn, f"{args.out}_waves_core_nhs")
