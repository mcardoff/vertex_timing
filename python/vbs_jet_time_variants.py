#!/usr/bin/env python3
"""Jet-time estimators side by side: util/vbs_time_veto runs that differ only in how the tagging jet's
time is built (--jet-time-dr cone, --jet-time-score sumpt|waves), all at the deck's settings (JVT +
fJVT loose, jet time from >= 2 tracks, 3 sigma, recoil > 100 GeV, WAVeS t0). Each variant's jet-time
sigma calibration f(n) is measured on its own signal HS legs.

Per variant:
  jet times   (truth-matched tagging jets with >= 1 timed track in HGTD, all selected events)
              HS jets timed, pileup jets timed, timed HS jets with a mostly-pileup cluster / > 60 ps,
              HS core sigma (narrow width of a double-Gaussian fit)
  the veto    S/sqrt(B) change vs no veto, all events and per case (Forward-Forward / Central-Forward),
              its difference to the reference variant (paired bootstrap: the same events in every
              variant), signal and Z+jets kept

    cd python && PYTHONNOUSERSITE=1 PYTHONPATH=/opt/homebrew/Cellar/root/6.40.04/lib/root \\
    ~/.venv-hgtd/bin/python vbs_jet_time_variants.py \\
        "dR<0.2, sum pT=jvtLoose" "no cone, WAVeS=jvtLoose_jtdr1p0_jtwaves" ...
(each argument: label=<tag>, the trees being condor/{vbf,zjets}/<sample>_<tag>_vbs_time_veto.root; the
first is the reference)
"""
import argparse, os
import numpy as np
import uproot
import ROOT
import vbs_time_veto_common as C

ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
ap.add_argument("variants", nargs="+", help="label=tag")
ap.add_argument("--dir", default=os.path.join(os.path.dirname(os.path.abspath(__file__)), "../condor"))
ap.add_argument("--nsigma", type=float, default=3.0)
ap.add_argument("--min-trk", type=int, default=2)
ap.add_argument("--nboot", type=int, default=300)
args = ap.parse_args()
C.ROOT = ROOT
ROOT.gErrorIgnoreLevel = ROOT.kWarning


def load(p):
    ev = uproot.open(p)["events"]
    e = ev.arrays([k for k in ev.keys() if not k.startswith("jet_")], library="np")
    o = np.lexsort((e["entry"], e["file_path"]))
    return {k: v[o] for k, v in e.items()}


V, KEY = [], None
for spec in args.variants:
    lab, tag = spec.rsplit("=", 1)
    S_all = load(f"{args.dir}/vbf/vbf_{tag}_vbs_time_veto.root")
    B_all = load(f"{args.dir}/zjets/zjets_{tag}_vbs_time_veto.root")
    key = (S_all["file_path"], S_all["entry"], B_all["file_path"], B_all["entry"])
    if KEY is None: KEY = key
    assert all(np.array_equal(a, b) for a, b in zip(KEY, key)), f"{lab}: event lists differ"
    cal = C.measure(S_all, True, f"HS legs, {lab}")
    # jet times, all selected events
    def legs(e, lab_):
        g = lambda k: np.concatenate([e[f"{L}_{k}"][e[f"{L}_{lab_}"]] for L in "ab"])
        return dict(w=np.concatenate([e["weight"][e[f"{L}_{lab_}"]] for L in "ab"]).astype(float), inh=g("hgtd_ntimed") >= 1,
                    n=g("core_n"), hsf=g("core_hsfrac"), dt=g("core_t") - np.concatenate([e["t0_truth"][e[f"{L}_{lab_}"]] for L in "ab"]))
    hs, pu = legs(S_all, "hs"), legs(B_all, "pu")
    frac = lambda w, num, den: 100 * w[num & den].sum() / w[den].sum()
    tm = hs["n"] >= args.min_trk
    J = dict(hs_timed=frac(hs["w"], tm, hs["inh"]), pu_timed=frac(pu["w"], pu["n"] >= args.min_trk, pu["inh"]),
             mostly_pu=frac(hs["w"], hs["hsf"] < 0.5, hs["inh"] & tm), tail=frac(hs["w"], np.abs(hs["dt"]) > 60, hs["inh"] & tm),
             core=C.dgaus_core(hs["dt"][hs["inh"] & tm])[0])
    ks, kb = S_all["met_truth"] > 100, B_all["z_pt"] > 100
    S = {k: v[ks] for k, v in S_all.items()}; B = {k: v[kb] for k, v in B_all.items()}
    vs = C.t0_veto(S, cal, "waves", args.min_trk, args.nsigma)[0]
    vb = C.t0_veto(B, cal, "waves", args.min_trk, args.nsigma)[0]
    V.append((lab, J, S, B, vs, vb))

S, B = V[0][2], V[0][3]
ws, wb = S["weight"].astype(float), B["weight"].astype(float)
nin = lambda e: (e["a_hgtd_ntimed"] >= 1).astype(int) + (e["b_hgtd_ntimed"] >= 1).astype(int)
CASES = {"all events": (np.ones(len(ws), bool), np.ones(len(wb), bool)),
         "Forward-Forward": (nin(S) == 2, nin(B) == 2), "Central-Forward": (nin(S) == 1, nin(B) == 1)}
rng = np.random.default_rng(5)
BS = [np.bincount(rng.integers(0, len(ws), len(ws)), minlength=len(ws)).astype(float) for _ in range(args.nboot)]
BB = [np.bincount(rng.integers(0, len(wb), len(wb)), minlength=len(wb)).astype(float) for _ in range(args.nboot)]


def z(vs, vb, ms, mb, bs=1.0, bb=1.0):
    a, b = ws * bs, wb * bb
    return 100 * ((a[ms & ~vs].sum() / a[ms].sum()) / np.sqrt(b[mb & ~vb].sum() / b[mb].sum()) - 1)


print(f"\n{'jet times':44s} {'HS timed':>9s} {'PU timed':>9s} {'HS mostly-PU':>13s} {'HS > 60 ps':>11s} {'HS core sigma':>14s}")
for lab, J, *_ in V:
    print(f"{lab:44s} {J['hs_timed']:8.1f}% {J['pu_timed']:8.1f}% {J['mostly_pu']:12.1f}% {J['tail']:10.1f}% {J['core']:11.1f} ps")
for cname, (ms, mb) in CASES.items():
    print(f"\nthe veto, {cname}: S/sqrt(B) change vs no veto (difference to '{V[0][0]}', paired), VBF kept, Z+jets kept")
    r0 = z(V[0][4], V[0][5], ms, mb); r0b = np.array([z(V[0][4], V[0][5], ms, mb, bs, bb) for bs, bb in zip(BS, BB)])
    for lab, J, _, _, vs, vb in V:
        r = z(vs, vb, ms, mb); rb = np.array([z(vs, vb, ms, mb, bs, bb) for bs, bb in zip(BS, BB)])
        print(f"  {lab:42s} {r:+6.1f}%  ({r - r0:+5.1f} +- {np.std(rb - r0b):.1f})   VBF {100 * ws[ms & ~vs].sum() / ws[ms].sum():5.1f}%  "
              f"Z+jets {100 * wb[mb & ~vb].sum() / wb[mb].sum():5.1f}%")
