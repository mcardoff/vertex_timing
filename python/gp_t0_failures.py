#!/usr/bin/env python3
"""Where the clustering t0 fails, for the Global Performance talk (slides/slides_gp_timing.tex): every selected event
put in one class, against the number of forward hard-scatter tracks with a time.

  core              |t0 - t_HS| < 60 ps
  selection         a time cluster within 60 ps of t_HS existed (its full-cluster time), the score chose another
  re-timing         (WAVeS) the chosen cluster is within 60 ps, its in-jet recomputed time is not
  no good cluster   no cluster within 60 ps: the hard scatter's tracks are mistimed, too few, or merged with pileup
  no t0             no cluster at all: no timed forward track compatible with the vertex
Inputs: util/tdr_t0 trees with the T0Extra columns (commit 57b5438): dt_<t0>_clean, n_clusters_clean, dt_best_clean,
dt_waves_full_clean. Cleaned track times, generator weights, the t0 selection of the talk.

  PYTHONNOUSERSITE=1 ~/.venv-hgtd/bin/python gp_t0_failures.py --trees ../condor/vbf/vbf_tdr_t0_v2.root \\
      ../condor/zjets/zjets_tdr_t0_v2.root --out ../figs/tdr/gp/gp_t0_failures
"""
import argparse

import numpy as np
import uproot

import tdr_root_style as S

ap = argparse.ArgumentParser()
ap.add_argument("--trees", nargs=2, required=True, metavar=("VBF", "ZJETS"))
ap.add_argument("--out", required=True)
args = ap.parse_args()
PASS, FOLD = 60.0, 20
CLASSES = [("core", "within 60 ps of t_{HS}", "#7BCB6B"),
           ("selection", "selection: a good cluster, another chosen", "#3F90DA"),
           ("retime", "re-timing: right cluster, in-jet time off", "#E69F00"),
           ("nogood", "no cluster within 60 ps of t_{HS}", "#B3231C"),
           ("not0", "no t_{0}: no timed forward track", "#9C9C9C")]


def classes(E, t0):
    dt = E[f"dt_{t0}_clean"]
    has = E["n_clusters_clean"] > 0
    core = has & (np.abs(dt) < PASS)
    good = has & (np.abs(E["dt_best_clean"]) < PASS)
    chosen_ok = np.abs(E["dt_waves_full_clean"]) < PASS if t0 == "waves" else np.abs(dt) < PASS
    retime = has & ~core & chosen_ok if t0 == "waves" else np.zeros(len(dt), bool)
    sel = has & ~core & ~retime & good
    nogood = has & ~core & ~good
    return {"core": core, "selection": sel, "retime": retime, "nogood": nogood, "not0": ~has}


R = S.setup()
cv = S.keep(R.TCanvas("fail", "", 1400, 700))
for c, (path, lab) in enumerate(zip(args.trees, ("VBF H#rightarrowinv.", "Z+jets"))):
    E = uproot.open(path)["events"].arrays(library="np")
    w = E["weight"].astype(float)
    n = np.minimum(E["n_hs_trk"], FOLD)
    print(f"\n{lab.replace('#rightarrow', '->')}: {len(w):,} events")
    for t0 in ("trkptz", "waves"):
        C = classes(E, t0)
        tot = sum(m.sum() for m in C.values())
        assert tot == len(w), (t0, tot, len(w))
        fails = w.sum() - w[C["core"]].sum()
        print(f"  {t0:7s} " + "  ".join(f"{k} {100 * w[m].sum() / w.sum():5.1f}%" for k, m in C.items())
              + "   | of the failures: " + "  ".join(f"{k} {100 * w[m].sum() / fails:4.0f}%" for k, m in C.items() if k != "core"))
        for lo, hi in ((0, 2), (3, 6), (7, FOLD)):
            mm = (n >= lo) & (n <= hi)
            print(f"      n_hs {lo}-{hi}{'+' if hi == FOLD else ''} ({100 * w[mm].sum() / w.sum():4.1f}% of events): "
                  + "  ".join(f"{k} {100 * w[m & mm].sum() / w[mm].sum():5.1f}%" for k, m in C.items()))
    C = classes(E, "waves")
    cv.cd(); S.pad(f"f_{c}", 0.5 * c, 0.0, 0.5 * c + 0.5, 1.0, top=0.06, bottom=0.13, left=0.13, right=0.04)
    st = S.keep(R.THStack(f"st_{c}", ""))
    leg = S.legend(0.17, 0.585, 0.80, 0.785, size=0.032, margin=0.12)
    hs = []
    den = np.histogram(n, np.arange(-0.5, FOLD + 1.5), weights=w)[0]
    for key, name, col in CLASSES:
        h = S.keep(R.TH1D(f"h_{c}_{key}", "", FOLD + 1, -0.5, FOLD + 0.5))
        num = np.histogram(n[C[key]], np.arange(-0.5, FOLD + 1.5), weights=w[C[key]])[0]
        for i, v in enumerate(num / np.maximum(den, 1e-12)): h.SetBinContent(i + 1, v)
        h.SetFillColor(S.colour(col)); h.SetLineColor(R.kBlack); h.SetLineWidth(1)
        st.Add(h); hs.append((h, name))
    fr = S.keep(R.gPad.DrawFrame(-0.5, 0.0, FOLD + 0.5, 1.85))   # headroom: header + legend above the stack
    S.axes(fr, "Forward hard-scatter tracks with a time (last bin: #geq 20)", "Fraction of events (WAVeS t_{0})", size=0.042, ytoff=1.3, xtoff=1.2)
    st.Draw("HIST SAME")
    for h, name in hs: leg.AddEntry(h, name, "F")
    leg.Draw()
    dist = S.keep(R.TH1D(f"d_{c}", "", FOLD + 1, -0.5, FOLD + 0.5))
    for i, v in enumerate(den / den.max() * 0.25): dist.SetBinContent(i + 1, v)
    dist.SetLineColor(R.kBlack); dist.SetLineStyle(2); dist.SetLineWidth(2); dist.Draw("HIST SAME")
    S.header(0.17, 0.88, ["#sqrt{s} = 14 TeV, #LT#mu#GT = 200", f"{lab}; dashed: the events' distribution (a.u.)"], size=0.034)
    R.gPad.RedrawAxis()
S.save(cv, args.out)
