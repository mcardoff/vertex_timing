#!/usr/bin/env python3
"""Compare tdr_pujet runs: each t0-gated tagger's pileup jet rejection relative to ITk only, against
hard-scatter efficiency, 2.4 < |eta| < 4.0, 30-50 GeV (left) and > 50 GeV (right), ATLAS style.

Colour = the t0 (HGTD, TRKPTZ, WAVeS, truth); line style = the run (first run dashed, second solid,
any further ones dotted and dash-dotted). Each run's ratio is to its own ITk-only curve. The console
prints the rejection ratios at the sampled working points and at 85% per |eta| bin, for every run.
A quick-test tool (e.g. track-time cleaning on/off, the gate width), not a TDR figure.

    cd python && PYTHONNOUSERSITE=1 ~/.venv-hgtd/bin/python tdr_pujet_compare.py \\
        ../figs/hists/tdr_pujet.root:"As in the ntuple" ../figs/hists/cleantimes_tdr_pujet.root:"TDR cleaning" \\
        --label "VBF H#rightarrowinv. (local sample)" --out ../figs/tdr/local_tdr_pujet_cleaning
"""
import argparse
import numpy as np
import uproot
import tdr_root_style as S

ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
ap.add_argument("runs", nargs="+", help="tdr_pujet output:legend label")
ap.add_argument("--label", default="VBF H#rightarrowinv.", help="sample, ROOT TLatex")
ap.add_argument("--plot-runs", type=int, default=2, help="draw only the first N runs (all are printed)")
ap.add_argument("--out", required=True)
args = ap.parse_args()
R = S.setup()

T0S = [("rpt_t0_hgtd", "HGTD (Athena) t_{0}", "#3B8ED0", 23), ("rpt_t0_trkptz", "TRKPTZ t_{0}", "#C0392B", 21),
       ("rpt_t0_waves", "WAVeS t_{0}", "#E69F00", 20), ("rpt_t0_truth", "Truth t_{0}", "#A0522D", 25)]
STYLES = [2, 1, 3, 4]
EFF = np.round(np.arange(0.800, 0.9751, 0.025), 3)
SLICES = [("30 < p_{T} < 50 GeV", 30.0, 50.0), ("p_{T} > 50 GeV", 50.0, 1e9)]
plain = lambda s: s.replace("_{0}", "0").replace("#", "")


def curve(hs, pu):
    th = np.unique(np.concatenate([hs, pu]))
    e_hs = 1.0 - np.searchsorted(np.sort(hs), th, side="left") / len(hs)
    e_pu = 1.0 - np.searchsorted(np.sort(pu), th, side="left") / len(pu)
    u, inv = np.unique(e_hs, return_inverse=True)
    p = np.full(len(u), np.inf); np.minimum.at(p, inv, e_pu)
    return u, p


def rejection(hs, pu, effs):
    pe = np.interp(effs, *curve(hs, pu))
    return np.where(pe > 0, 1.0 / pe, np.nan)


runs = []
for spec in args.runs:
    path, _, lab = spec.partition(":")
    t = uproot.open(path)
    J = t["jets"].arrays(library="np")
    meta = t["meta"].arrays(library="np")
    gate = float(np.unique(meta["gate_sigma"])[0])
    clean = bool(np.unique(meta["clean_times"])[0]) if "clean_times" in meta else False
    runs.append((lab or path, J, gate, clean))
    print(f"{lab or path}: {path}  (gate {gate:g} sigma, track-time cleaning {'on' if clean else 'off'})")


def ratio(J, col, lo, hi, elo, ehi, effs):
    ae = np.abs(J["eta"])
    m = (J["pt"] > lo) & (J["pt"] < hi) & (ae >= elo) & (ae < ehi)
    H, P = m & (J["hs"] == 1), m & (J["pu"] == 1)
    return rejection(J[col][H], J[col][P], effs) / rejection(J["rpt_itk"][H], J["rpt_itk"][P], effs)


cv = S.keep(R.TCanvas("cmp", "", 1600, 600))
for c, (lab, lo, hi) in enumerate(SLICES):
    print(f"\n{plain(lab)}, 2.4 < |eta| < 4.0: rejection relative to ITk at " + " / ".join(f"{e:.3f}" for e in EFF))
    vals = {}
    for k, (lab_r, J, *_ ) in enumerate(runs):
        for col, name, *_ in T0S:
            vals[k, col] = ratio(J, col, lo, hi, 2.4, 4.0, EFF)
            print(f"  {lab_r[:24]:24s} {plain(name):20s} " + " ".join(f"{x:5.2f}" for x in vals[k, col]))
    cv.cd(); S.pad(f"p{c}", 0.5 * c, 0.0, 0.5 * c + 0.5, 1.0)
    shown = [v for (k, _), v in vals.items() if k < args.plot_runs]
    ymax = max(2.0, 1.45 * np.nanmax(shown))
    fr = S.keep(R.gPad.DrawFrame(0.80, 0.0, 1.0, ymax))
    S.axes(fr, "Hard-scatter jet efficiency", "Rejection relative to ITk only")
    fr.GetXaxis().SetNdivisions(505); fr.GetYaxis().ChangeLabel(1, -1, 0)
    one = S.keep(R.TLine(0.80, 1.0, 1.0, 1.0)); one.SetLineColor(R.kGray + 2); one.SetLineStyle(2); one.Draw()
    leg = S.legend(0.62, 0.50, 0.93, 0.68, size=0.036)
    for k in range(min(args.plot_runs, len(runs))):
        for col, name, colour, mk in T0S:
            g = S.style(S.graph(EFF, vals[k, col]), colour, mk, STYLES[k], 2, 1.1 if k else 0.9)
            g.Draw("LP SAME")
            if k == 0: leg.AddEntry(g, name, "LP")
    leg.Draw()
    leg2 = S.legend(0.20, 0.17, 0.60, 0.17 + 0.055 * min(args.plot_runs, len(runs)), size=0.036)
    for k in range(min(args.plot_runs, len(runs))):
        d = S.keep(R.TGraph(1)); d.SetLineStyle(STYLES[k]); d.SetLineWidth(2); d.SetLineColor(R.kGray + 3)
        leg2.AddEntry(d, runs[k][0], "L")
    leg2.Draw()
    S.header(0.20, 0.87, ["#sqrt{s} = 14 TeV, HL-LHC, #LT#mu#GT = 200", args.label, f"{lab}, 2.4 < |#eta| < 4.0"])
S.save(cv, args.out)

EB = [2.4, 2.8, 3.2, 3.6, 4.0]
print("\n30-50 GeV, rejection relative to ITk at 85%, |eta| bins " + "  ".join(f"{EB[b]}-{EB[b + 1]}" for b in range(4)))
for k, (lab_r, J, *_) in enumerate(runs):
    for col, name, *_ in T0S:
        print(f"  {lab_r[:24]:24s} {plain(name):20s} " + "  ".join(f"{ratio(J, col, 30, 50, EB[b], EB[b + 1], np.array([0.85]))[0]:7.2f}" for b in range(4)))
