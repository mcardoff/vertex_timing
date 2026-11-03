#!/usr/bin/env python3
"""The vertex t0 figures of the TDR-structured presentation (util/tdr_t0's tree and fits), ATLAS style.

  <out>_resid     t0 - t_HS for TRKPTZ (left) and WAVeS (right), from the ntuple's track times (dashed)
                  and with the TDR's track-time cleaning (solid), Athena's HGTD t0 (grey, cannot be
                  re-cleaned) as the reference; normalised to the selected events, so a missing t0
                  lowers the area. Lower panels: cleaned / as is in 20 ps bins.
  <out>_core_nhs  core fraction (|t0 - t_HS| < 60 ps over all events, a missing t0 a failure) against
                  the number of forward hard-scatter tracks with a time (folded at 20, as the main
                  analysis); lower panel: cleaned minus as is.

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
args = ap.parse_args()
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

T0S = [("trkptz", "TRKPTZ t_{0}", "#C0392B"), ("waves", "WAVeS t_{0}", "#E69F00")]
MAIN, RAT = 0.68, 0.32
TS_M, TS_R = 0.05 * 600 / (800 * MAIN), 0.05 * 600 / (800 * RAT)

# ── residuals ─────────────────────────────────────────────────────────────────
cv = S.keep(R.TCanvas("resid", "", 1600, 800))
for c, (key, title, colour) in enumerate(T0S):
    cv.cd(); S.pad(f"m{c}", 0.5 * c, RAT, 0.5 * c + 0.5, 1.0, top=0.05 / MAIN * 0.6, bottom=0.02)
    R.gPad.SetLogy()
    hs = {}
    for name, col, ls, w in [("hgtd", "#8C8C8C", 1, 2), (f"{key}_raw", colour, 2, 3), (f"{key}_clean", colour, 1, 3)]:
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
for name, lab, col, mk, ls in [("hgtd", "HGTD (Athena) t_{0}", "#8C8C8C", 22, 1),
                               ("trkptz_raw", "TRKPTZ t_{0}, as is", "#C0392B", 25, 2), ("trkptz_clean", "TRKPTZ t_{0}, cleaned", "#C0392B", 21, 1),
                               ("waves_raw", "WAVeS t_{0}, as is", "#E69F00", 24, 2), ("waves_clean", "WAVeS t_{0}, cleaned", "#E69F00", 20, 1)]:
    f, e = per_bin(core[name])
    g = S.style(S.graph(bins, f, e), col, mk, ls, 2, 1.1); g.Draw("LP SAME")
    leg.AddEntry(g, lab, "LP")
leg.Draw()
S.header(0.20, 0.87, [f"{ENERGY}, {args.label}", "Core fraction: |t_{0} #minus t_{HS}| < 60 ps, all events"], size=TS_M)
cv.cd(); S.pad("r", 0.0, 0.0, 1.0, RAT, top=0.03, bottom=0.36)
fr2 = S.keep(R.gPad.DrawFrame(-0.5, -0.04, FOLD + 0.5, 0.06))
S.axes(fr2, f"Forward hard-scatter tracks with a time (last bin: #geq {FOLD})", "Cleaned #minus as is", size=TS_R, ytoff=1.35 * TS_M / TS_R, xtoff=1.2, ndiv_y=505)
zero = S.keep(R.TLine(-0.5, 0, FOLD + 0.5, 0)); zero.SetLineStyle(2); zero.SetLineColor(R.kGray + 2); zero.Draw()
print(f"\ncore fraction, cleaned minus as is, by forward HS tracks (0..{FOLD - 1}, >= {FOLD}); events per bin:")
print("  " + " ".join(f"{np.sum(x == b):6d}" for b in bins))
for k, (key, title, colour) in enumerate(T0S):
    a, b_ = core[f"{key}_raw"], core[f"{key}_clean"]
    d, de = [], []
    for b in bins:
        sel = x == b; n = sel.sum()
        gain, loss = np.sum(sel & ~a & b_), np.sum(sel & a & ~b_)
        d.append((gain - loss) / n if n else np.nan); de.append(np.sqrt(gain + loss) / n if n else np.nan)
    g = S.style(S.graph(bins + (k - 0.5) * 0.25, d, de), colour, 21 if key == "trkptz" else 20, 1, 2, 1.0); g.Draw("P SAME")
    print(f"  {key:7s} " + " ".join(f"{100 * v:+6.2f}" for v in d))
S.save(cv, f"{args.out}_core_nhs")
