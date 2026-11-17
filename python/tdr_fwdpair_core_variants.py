"""The veto's jet time with more tracks, or a z window (results/jet_core_variants.md, todo vt-r2k.21 / vt-r2k.20).

For each core jet time on tdr_fwdpair trees -- the standard one (dR < 0.2, 1 < pT < 30 GeV) and the variants
j*_core_{t,sig,n}_<v> (r04, nocap, r04nocap, zwin, zpref; util/tdr_fwdpair.cxx CV_NAMES) -- the timing veto rebuilt on it,
with its own calibration f(n) from the VBF hard-scatter legs (tdr_fwdpair_common.calibration), and:
  - the hard-scatter legs in HGTD: timed share, median core tracks, out of time (pull >= 3 vs the truth t0), and the
    out-of-time rate of single-track cores;
  - the pileup legs in HGTD (Z+jets): timed share, median core tracks, the median 3 sigma window, and the share inside;
  - S/sqrt(B) with the WAVeS and the truth t0, before the BDT and at its cuts (--bdt, matched by row order).
Without a Z+jets tree (--bkg) only the hard-scatter side is shown.

    cd python && PYTHONNOUSERSITE=1 ~/.venv-hgtd/bin/python tdr_fwdpair_core_variants.py <vbf tree> \\
        [--bkg <zjets tree>] [--bdt <kinbdt npz>] [--effs 0.8]
"""
import argparse

import numpy as np
import uproot

import tdr_fwdpair_common as FC

VARIANTS = [("", "dR < 0.2, 1-30 GeV (standard)"), ("r04", "dR < 0.4"), ("nocap", "no 30 GeV cap"),
            ("r04nocap", "dR < 0.4, no cap"), ("zwin", "dR < 0.2, PV0 z window"),
            ("zpref", "dR < 0.2, z-compatible cluster preferred")]
ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
ap.add_argument("sig")
ap.add_argument("--bkg")
ap.add_argument("--bdt")
ap.add_argument("--effs", default="0.8")
args = ap.parse_args()
S0 = uproot.open(args.sig)["events"].arrays(library="np")
B0 = uproot.open(args.bkg)["events"].arrays(library="np") if args.bkg else None
D = np.load(args.bdt) if args.bdt else None
if D is not None:
    assert len(D["sig_score"]) == len(S0["weight"]) and len(D["bkg_score"]) == len(B0["weight"]), "BDT does not match"
fwd = lambda E, k: np.abs(E[f"j{k}_eta"]) > 2.4


def view(E, v):
    """E with the core-time columns replaced by variant v's"""
    if not v: return E
    X = dict(E)
    for k in (1, 2):
        for q in ("t", "sig", "n"): X[f"j{k}_core_{q}"] = E[f"j{k}_core_{q}_{v}"]
    return X


def eff(E, kept, base):
    w = E["weight"].astype(float)
    return w[base & kept].sum() / w[base].sum()


cuts = [("no BDT", None)]
if D is not None:
    wS = S0["weight"].astype(float); o = np.argsort(-D["sig_score"]); cs = np.cumsum(wS[o]) / wS.sum()
    cuts += [(f"BDT {e:.0%}", D["sig_score"][o][np.searchsorted(cs, e)]) for e in (float(x) for x in args.effs.split(","))]

for v, lab in VARIANTS:
    S = view(S0, v)
    f = FC.calibration(S, verbose=False)
    P = FC.pulls(S, f)
    m = [(S[f"j{k}_hs"] == 1) & fwd(S, k) for k in (1, 2)]
    n = np.concatenate([S[f"j{k}_core_n"][m[k - 1]] for k in (1, 2)])
    p = np.concatenate([np.nan_to_num(P[f"truth_{k}"][m[k - 1]], nan=0.0) for k in (1, 2)])
    t = n >= 1
    print(f"\n== {lab} ==")
    print(f"  HS legs in HGTD: timed {100 * t.mean():.1f}%, median core tracks {np.median(n[t]):.0f}, single-track "
          f"{100 * np.mean(n[t] == 1):.0f}%; out of time (truth t0) {100 * np.mean(p[t] >= 3):.2f}% "
          f"(single-track {100 * np.mean(p[t & (n == 1)] >= 3):.1f}%, >= 3 tracks {100 * np.mean(p[t & (n >= 3)] >= 3):.2f}%)")
    if B0 is None:
        continue
    B = view(B0, v)
    PB = FC.pulls(B, f)
    mb = [(B[f"j{k}_pu"] == 1) & fwd(B, k) for k in (1, 2)]
    nb = np.concatenate([B[f"j{k}_core_n"][mb[k - 1]] for k in (1, 2)])
    wb = np.concatenate([B["weight"][mb[k - 1]].astype(float) for k in (1, 2)])
    with np.errstate(invalid="ignore", divide="ignore"):
        den = np.concatenate([np.abs(B[f"j{k}_core_t"] - B["t0_waves"])[mb[k - 1]] / PB[f"waves_{k}"][mb[k - 1]] for k in (1, 2)])
    pw = np.concatenate([np.nan_to_num(PB[f"waves_{k}"][mb[k - 1]], nan=0.0) for k in (1, 2)])
    tb = nb >= 1
    ok = tb & np.isfinite(den)
    print(f"  PU legs in HGTD: timed {100 * np.average(tb, weights=wb):.1f}%, median core tracks {np.median(nb[tb]):.0f}, "
          f"median 3-sigma window (WAVeS) {3 * np.median(den[ok]):.0f} ps; timed and inside it {100 * np.average(pw[tb] < 3, weights=wb[tb]):.1f}%")
    for cl, thr in cuts:
        bS = np.ones(len(S["weight"]), bool) if thr is None else D["sig_score"] >= thr
        bB = np.ones(len(B["weight"]), bool) if thr is None else D["bkg_score"] >= thr
        KS, KB = FC.vetoes(S, f), FC.vetoes(B, f)
        out = []
        for src in ("waves", "truth"):
            eS, eB = eff(S, KS[src], bS), eff(B, KB[src], bB)
            out.append(f"{src} {100 * (eS / np.sqrt(eB) - 1):+.1f}% (eps_S {eS:.3f}, eps_B {eB:.3f})")
        print(f"  S/sqrt(B), {cl}: " + "; ".join(out))
