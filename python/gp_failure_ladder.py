#!/usr/bin/env python3
"""util/failure_decomposition's principal-cause ladder, applied offline (weighted) to its per-event trees, for the
Global Performance talk (slides/slides_gp_timing.tex).

The program's ladder (first match wins), its own definitions:
  1. timing misassignment  (a) the ideal-resolution world (every track at its truth time, 1 ps) passes, or (b) > 50%
                           of the hard scatter's timed pT is >= 3 sigma from the truth vertex time
  2. selection             another cluster passes
  3. misclustering         a cluster's hard-scatter core alone (its pileup tracks removed) passes
  4. track efficiency      the ideal world passes once the untimed tracks get a time too
  5. irreducible
The talk's version adds two things measured since: (0) re-timing first -- WAVeS's chosen cluster passes on its full
time and fails on its in-jet one -- and selection by the hard scatter's OWN cluster (most HS pT, > 50% HS) instead of
"another cluster passes", which a pileup cluster satisfies by chance (~21% per cluster). Timing is split into its two
routes: (b) "the hard scatter's track times are wrong" and (a) only, "resolution".
Run with --clean-times --assoc-nsigma=2.5 --waves-time=injet the program reproduces the talk's t0s event by event.
Those options tag its output, <sample>_cleant0_assoc2p5_wavesinjet_failure_decomposition.root, and are recorded per
row (clean_times, assoc_nsigma, waves_time_injet). Trees made before 2026-11-17 carry neither: the talk's were written
untagged, as <sample>_failure_decomposition.root.

  PYTHONNOUSERSITE=1 ~/.venv-hgtd/bin/python gp_failure_ladder.py \\
      --trees ../condor/vbf/vbf_cleant0_assoc2p5_wavesinjet_failure_decomposition.root \\
      ../condor/zjets/zjets_cleant0_assoc2p5_wavesinjet_failure_decomposition.root --out ../figs/tdr/gp/gp_failure_ladder
"""
import argparse

import numpy as np
import uproot

import tdr_root_style as S

ap = argparse.ArgumentParser()
ap.add_argument("--trees", nargs="+", required=True, help="one panel per tree (the talk: VBF, or VBF and Z+jets)")
ap.add_argument("--out", required=True)
ap.add_argument("--no-ideal", action="store_true", help="drop the ideal-time oracles (timing route (a), track efficiency): on the "
                "grid samples pileup tracks have no truth link (vertex or particle), so the ideal world draws each one an "
                "independent random time, dissolving the pileup clusters -- the oracle passes ~16 points more often than on "
                "the local sample, where it is valid. What remains is 'no usable HS cluster'.")
ap.add_argument("--labels", nargs="+", default=("VBF H#rightarrowinv.", "Z+jets"))
ap.add_argument("--failures-only", action="store_true", help="the WAVeS failures alone: each bin's causes stack to 1 "
                "(fraction of the failures), the dashed line the failures (a.u.) instead of the events")
ap.add_argument("--fold", type=int, default=20, help="last bin: >= this many tracks (15 for --failures-only: few failures above it)")
args = ap.parse_args()
NP = len(args.trees)
if len(args.labels) < NP: raise SystemExit("one --labels entry per tree")
FOLD = args.fold
TALK = [("core", "within 60 ps of t_{HS}", "#7BCB6B"),
        ("retime", "re-timing: right cluster, in-jet time off", "#E69F00"),
        ("badtime", "misassigned HS track times (removing them fixes it)", "#5A3E2B"),
        ("resol", "timing resolution (ideal times recover it)", "#C98B5E"),
        ("selection", "selection: the HS cluster was right, not chosen", "#3F90DA"),
        ("miscl", "misclustering: pileup pulls the HS cluster off", "#B3231C"),
        ("teff", "track efficiency: untimed HS tracks", "#8E5DB5"),
        ("irred", "irreducible", "#9C9C9C")]
TALK_NOIDEAL = [x for x in TALK if x[0] not in ("resol", "teff", "irred")] + [("nousable", "no usable HS cluster (too few / off HS tracks)", "#8E5DB5")]


def ladder(E, s, talk):
    """principal cause per event for score s ('waves' / 'trkptz'); talk = the talk's version"""
    fail = ~E[f"pass_{s}"]
    out, left = {"core": ~fail}, fail.copy()
    def take(key, m):
        out[key] = left & m
        left[:] = left & ~m
    if talk:
        take("retime", E["retime_waves"] if s == "waves" else np.zeros(len(fail), bool))
        # misassigned HS times caused it: most of the HS pT mistimed, OR removing the mistimed HS tracks fixes it
        # (fix_mistimed_*, trees from 2026-11-10 on; older trees have only the majority rung)
        take("badtime", E["badtime"] | (E[f"fix_mistimed_{s}"] if f"fix_mistimed_{s}" in E else False))
        if not args.no_ideal: take("resol", E[f"ideal_{s}"])
        take("selection", E[f"sel_hs_{s}"])
        if args.no_ideal:
            take("miscl", E["miscl"])
            out["nousable"] = left.copy()
            return out
    else:
        take("timing", E["badtime"] | E[f"ideal_{s}"])
        take("selection", E[f"sel_any_{s}"])
    take("miscl", E["miscl"])
    take("teff", E[f"teff_{s}"])
    out["irred"] = left.copy()
    return out


R = S.setup()
cv = S.keep(R.TCanvas("ladder", "", 700 * NP, 760))
for c, (path, lab) in enumerate(zip(args.trees, args.labels)):
    E = uproot.open(path)["events"].arrays(library="np")
    w = E["weight"].astype(float)
    n = np.minimum(E["n_hs_trk"], FOLD)
    print(f"\n{lab.replace('#rightarrow', '->')}: {len(w):,} events")
    for s in ("waves", "trkptz"):
        for talk in (False, True):
            C = ladder(E, s, talk)
            fails = w[~E[f"pass_{s}"]].sum()
            print(f"  {s:7s} {'talk' if talk else 'original':8s} core {100 * w[C['core']].sum() / w.sum():5.1f}% | of the failures: "
                  + "  ".join(f"{k} {100 * w[m].sum() / fails:4.1f}%" for k, m in C.items() if k != "core"))
    fw = ~E["pass_waves"]
    print(f"  WAVeS failures with > 50% of HS pT mistimed (any rung): {100 * w[fw & E['badtime']].sum() / w[fw].sum():.1f}%; "
          f"'another cluster passes' {100 * w[fw & E['sel_any_waves']].sum() / w[fw].sum():.1f}% vs the HS cluster passes "
          f"{100 * w[fw & E['sel_hs_waves']].sum() / w[fw].sum():.1f}%")
    C = ladder(E, "waves", True)
    cv.cd(); S.pad(f"ld_{c}", c / NP, 0.0, (c + 1) / NP, 1.0, top=0.06, bottom=0.12, left=0.13, right=0.04)
    st = S.keep(R.THStack(f"lst_{c}", ""))
    fsel = ~E["pass_waves"] if args.failures_only else np.ones(len(w), bool)
    den = np.histogram(n[fsel], np.arange(-0.5, FOLD + 1.5), weights=w[fsel])[0]
    if args.failures_only:
        nf = np.histogram(n[fsel], np.arange(-0.5, FOLD + 1.5))[0]
        print("  failures per bin (unweighted): " + " ".join(f"{i}:{v}" for i, v in enumerate(nf)))
    leg = S.legend(0.17, 0.52, 0.93 if args.failures_only else 0.80, 0.755, size=0.027, margin=0.10 if args.failures_only else 0.12)
    wf = w[~E["pass_waves"]].sum()
    for key, name, col in (TALK_NOIDEAL if args.no_ideal else TALK):
        if args.failures_only and key == "core": continue
        h = S.keep(R.TH1D(f"lh_{c}_{key}", "", FOLD + 1, -0.5, FOLD + 0.5))
        num = np.histogram(n[C[key]], np.arange(-0.5, FOLD + 1.5), weights=w[C[key]])[0]
        for i, v in enumerate(num / np.maximum(den, 1e-12)): h.SetBinContent(i + 1, v)
        h.SetFillColor(S.colour(col)); h.SetLineColor(R.kBlack); h.SetLineWidth(1)
        # failures only: each cause's share of the WAVeS failures (weighted) and its event count, in the legend
        lab_ = f"{100 * w[C[key]].sum() / wf:.1f}% ({C[key].sum() / 1e3:.1f}k)  {name}" if args.failures_only else name
        st.Add(h); leg.AddEntry(h, lab_, "F")
    fr = S.keep(R.gPad.DrawFrame(-0.5, 0.0, FOLD + 0.5, 2.25))   # headroom: header + legend above the stack
    S.axes(fr, f"Forward hard-scatter tracks with a time (last bin: #geq {FOLD})",
           "Fraction of WAVeS t_{0} failures" if args.failures_only else "Fraction of events (WAVeS t_{0})",
           size=0.042, ytoff=1.3, xtoff=1.15)
    st.Draw("HIST SAME")
    leg.Draw()
    dist = S.keep(R.TH1D(f"ld_{c}", "", FOLD + 1, -0.5, FOLD + 0.5))
    for i, v in enumerate(den / den.max() * 0.25): dist.SetBinContent(i + 1, v)
    dist.SetLineColor(R.kBlack); dist.SetLineStyle(2); dist.SetLineWidth(2); dist.Draw("HIST SAME")
    S.header(0.17, 0.88, ["#sqrt{s} = 14 TeV, #LT#mu#GT = 200", lab,
                          "failure_decomposition's ladder, first match; dashed: "
                          + ("failures (a.u.)" if args.failures_only else "events (a.u.)")], size=0.030)
    R.gPad.RedrawAxis()
S.save(cv, args.out)
