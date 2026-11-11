#!/usr/bin/env python3
"""The failure-mode fixes against the number of forward hard-scatter tracks with a time, for the Global Performance talk
(slides/slides_gp_timing.tex): core fraction (upper) and resolution (lower) per bin, WAVeS and its fixed rows.

Reads clustering_plot's failure_fixes.json ("nhs" per row: the main program's own n_HS bins, its underflow -- n <= 2,
below HS_TRACK_MIN, where most failures are -- as the first bin, the last bin the fold >= 20). sigma_peak is FIT_PEAKTAIL
on the bin's t0 - t_HS, as the main program's slices; bins under 100 events have none.
  --mode alone   each fix alone on WAVeS (src/failure_fixes.h)
  --mode chain   the fixes added one after another, in the failure ladder's order

  PYTHONNOUSERSITE=1 ~/.venv-hgtd/bin/python gp_failure_fixes_nhs.py --mode alone \\
      --json ../condor/vbf/fixes/vbf_cleant0_assoc2p5_failure_fixes.json \\
      ../condor/zjets/fixes/zjets_cleant0_assoc2p5_failure_fixes.json --out ../figs/tdr/gp/gp_fix_alone_nhs
"""
import argparse
import json

import numpy as np

import tdr_root_style as S

ap = argparse.ArgumentParser()
ap.add_argument("--json", nargs="+", required=True)
ap.add_argument("--labels", nargs="+", default=None)
ap.add_argument("--mode", choices=["alone", "chain"], required=True)
ap.add_argument("--out", required=True)
args = ap.parse_args()
NAMES = {"vbf": "VBF H#rightarrowinv.", "zjets": "Z+jets", "local": "VBF H#rightarrowinv. (local sample)"}
LABELS = args.labels or [NAMES.get(json.load(open(p))["sample"], json.load(open(p))["sample"]) for p in args.json]
# (row, legend, colour, marker); WAVeS in black here: its deck colour is too close to re-timing's
BASE = ("WAVES", "WAVeS t_{0}", "#000000", 20)
ROWS = {"alone": [BASE,
                  ("WFIX_RETIME", "+ re-timing fixed", "#E69F00", 21),
                  ("WFIX_MISTIMED", "+ misassigned HS times fixed", "#5A3E2B", 22),
                  ("WFIX_SELECT", "+ selection fixed", "#3F90DA", 23),
                  ("WFIX_MISCL", "+ misclustering fixed", "#B3231C", 33),
                  ("WFIX_TEFF", "+ HS track efficiency fixed", "#8E5DB5", 34)],
        "chain": [BASE,
                  ("WFIX_RETIME", "re-timing", "#E69F00", 21),
                  ("WFIX_C2", "+ misassigned HS times", "#5A3E2B", 22),
                  ("WFIX_C3", "+ selection", "#3F90DA", 23),
                  ("WFIX_C4", "+ misclustering", "#B3231C", 33),
                  ("WFIX_C5", "+ HS track efficiency", "#8E5DB5", 34),
                  ("WFIX_C6", "+ HS time resolution", "#C98B5E", 29)]}[args.mode]
TITLE = {"alone": "each failure mode fixed alone, with truth", "chain": "failure modes fixed in turn, with truth"}[args.mode]

R = S.setup()
ncol = len(args.json)
cv = S.keep(R.TCanvas("fixnhs", "", 800 * ncol, 900))
for c, (path, lab) in enumerate(zip(args.json, LABELS)):
    rows = {r["short"]: r for r in json.load(open(path))["rows"]}
    nb = rows["WAVES"]["nhs"]
    x = np.array(nb["x"], float)
    # the first entry is the underflow: n <= 2 in files written while HS_TRACK_MIN was 2.5 (drawn at 2, "<= 2"),
    # empty since the axis starts at 0 (2026-11-11; dropped)
    keep = np.ones(len(x), bool)
    if nb["total"][0] > 0: x[0] = 2.0
    else: keep[0] = False
    LUMPED = nb["total"][0] > 0
    x = x[keep]
    sel = lambda a: np.asarray(a, float)[keep]
    tot0 = sel(nb["total"])
    w_core = 100 * rows["WAVES"]["core"]
    off = np.linspace(-0.24, 0.24, len(ROWS))
    # --- upper: core fraction
    cv.cd(); S.pad(f"u_{c}", c / ncol, 0.40, (c + 1) / ncol, 1.0, top=0.05, bottom=0.02, left=0.14, right=0.04)
    allc = [sel(rows[k]["nhs"]["pass"]) / np.maximum(sel(rows[k]["nhs"]["total"]), 1) for k, *_ in ROWS]
    lo = max(0.0, np.floor(10 * min(np.min(a[tot0 > 0]) for a in allc)) / 10 - 0.1)
    hi = 1.0 + 0.42 * (1.0 - lo)
    fr = S.keep(R.gPad.DrawFrame(x[0] - 0.6, lo, x[-1] + 0.6, hi))
    S.axes(fr, "", "Core fraction (|t_{0} #minus t_{HS}| < 60 ps)", size=0.048, ytoff=1.35)
    fr.GetXaxis().SetLabelSize(0)
    one = S.keep(R.TLine(x[0] - 0.6, 1.0, x[-1] + 0.6, 1.0)); one.SetLineColor(R.kGray); one.Draw()
    dist = S.keep(R.TGraph(len(x)))
    for i, (xi, t) in enumerate(zip(x, tot0)): dist.SetPoint(i, xi, lo + 0.16 * (1.0 - lo) * t / tot0.max())
    dist.SetLineStyle(2); dist.SetLineColor(R.kGray + 1); dist.SetLineWidth(2); dist.Draw("L")
    # under the curves: from ~5 tracks on they sit near 1, the event distribution low
    leg = S.legend(0.36, 0.545 - 0.052 * (len(ROWS) + 1), 0.94, 0.545, size=0.036, margin=0.12)
    leg.SetHeader("in brackets: gain over all events [points]")
    print(f"{lab.replace('#rightarrow', '->')} ({args.mode}): events per bin " + " ".join(f"{int(t)}" for t in tot0))
    for j, (k, name, col, mk) in enumerate(ROWS):
        r = rows[k]["nhs"]; p, t = sel(r["pass"]), sel(r["total"])
        m = t > 0
        eff = np.where(m, p / np.maximum(t, 1), np.nan); err = np.where(m, np.sqrt(eff * (1 - eff) / np.maximum(t, 1)), np.nan)
        g = S.style(S.graph(x + off[j], eff, err), col, mk, 1, 2, 1.1); g.Draw("PL")
        gain = 100 * rows[k]["core"] - w_core
        leg.AddEntry(g, f"{name}  ({100 * rows[k]['core']:.1f}%)" if k == "WAVES" else f"{name}  ({gain:+.1f})", "LP")
        print(f"   {k:14s} " + " ".join(f"{100 * e:5.1f}" for e in eff[:8]))
    leg.Draw()
    S.header(0.18, 0.90, ["#sqrt{s} = 14 TeV, #LT#mu#GT = 200;  grey dashed: events (a.u.)",
                          f"{lab}: WAVeS t_{{0}}, {TITLE}"], size=0.038)
    R.gPad.RedrawAxis()
    # --- lower: resolution
    cv.cd(); S.pad(f"d_{c}", c / ncol, 0.0, (c + 1) / ncol, 0.40, top=0.03, bottom=0.25, left=0.14, right=0.04)
    # a failed or unconstrained slice fit (error > half its value) is not drawn; the HS-time-resolution row sits at
    # the fit's 3 ps floor in every bin and is left out
    def good(k):
        s_, e_ = sel(rows[k]["nhs"]["sigma_peak"]), sel(rows[k]["nhs"]["sigma_peak_err"])
        ok = (s_ > 0) & (e_ > 0) & (e_ < 0.5 * s_)
        return np.where(ok, s_, np.nan), np.where(ok, e_, 0.0)
    SROWS = [r for r in ROWS if r[0] != "WFIX_C6"]
    sig = [good(k)[0] for k, *_ in SROWS]
    smax = max(np.nanmax(s) for s in sig if np.any(np.isfinite(s)))
    fd = S.keep(R.gPad.DrawFrame(x[0] - 0.6, 0.0, x[-1] + 0.6, 1.15 * smax))
    S.axes(fd, "Forward hard-scatter tracks with a time", "#sigma_{peak} [ps]", size=0.072, ytoff=0.85, xtoff=1.05, ndiv_y=505)
    fd.GetXaxis().SetLabelSize(0)
    for xi in range(int(x[0]), int(x[-1]) + 1, 2):
        t = S.keep(R.TLatex(xi, -0.05 * 1.15 * smax, "#leq 2" if (xi == 2 and LUMPED) else ("#geq 20" if xi == int(x[-1]) else str(xi))))
        t.SetTextAlign(23); t.SetTextSize(0.068); t.Draw()
    for j, (k, name, col, mk) in enumerate(ROWS):
        if k == "WFIX_C6": continue
        s_, e_ = good(k)
        g = S.style(S.graph(x + off[j], s_, e_), col, mk, 1, 2, 1.1); g.Draw("PL")
    if any(r[0] == "WFIX_C6" for r in ROWS):
        t = S.keep(R.TLatex(x[-1] + 0.4, 1.08 * smax, "HS time resolution: at the fit's 3 ps floor, not drawn"))
        t.SetTextAlign(33); t.SetTextSize(0.060); t.SetTextColor(R.kGray + 2); t.Draw()
    R.gPad.RedrawAxis()
S.save(cv, args.out)
