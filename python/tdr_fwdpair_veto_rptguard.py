"""An R_pT guard on the timing veto (results/jet_core_variants.md, todo vt-r2k.21).

The hard-scatter legs the veto wrongly removes are mostly jets whose core time is a PILEUP track (71% of out-of-time
hard-scatter cores on local VBF; util/scratch/core_tails_diag.cxx): their own tracks sit at PV0 in z but carry no time.
Their ITk R_pT stays high, a pileup jet's is usually low -- and above 60 GeV, where most of the veto's gain is
(results/veto_ptsplit.md), JVT never looks at R_pT. So: a leg may veto on its time only if its ITk R_pT (j*_rpt_itk) is
below X. Scanned on tdr_fwdpair trees with the analysis BDT, as tdr_fwdpair_kinplots.py --veto.

    cd python && PYTHONNOUSERSITE=1 ~/.venv-hgtd/bin/python tdr_fwdpair_veto_rptguard.py <vbf tree> <zjets tree> <kinbdt npz>
"""
import argparse

import numpy as np
import uproot

import tdr_fwdpair_common as FC

ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
ap.add_argument("sig"); ap.add_argument("bkg"); ap.add_argument("bdt")
ap.add_argument("--guards", default="inf,0.5,0.3,0.2,0.15,0.1,0.05", help="R_pT thresholds below which a leg may veto")
ap.add_argument("--core", default="", help="a core-time variant suffix (tdr_fwdpair CV_NAMES), default the standard core")
args = ap.parse_args()
S = uproot.open(args.sig)["events"].arrays(library="np")
B = uproot.open(args.bkg)["events"].arrays(library="np")
for E in (S, B):
    if args.core:
        for k in (1, 2):
            for q in ("t", "sig", "n"): E[f"j{k}_core_{q}"] = E[f"j{k}_core_{q}_{args.core}"]
D = np.load(args.bdt)
assert len(D["sig_score"]) == len(S["weight"]) and len(D["bkg_score"]) == len(B["weight"])
f = FC.calibration(S, verbose=False)
wS, wB = S["weight"].astype(float), B["weight"].astype(float)
o = np.argsort(-D["sig_score"]); cs = np.cumsum(wS[o]) / wS.sum()
cuts = [("no BDT", None), ("BDT 80%", D["sig_score"][o][np.searchsorted(cs, 0.8)])]
fwd = lambda E, k: np.abs(E[f"j{k}_eta"]) > 2.4
PS, PB = FC.pulls(S, f), FC.pulls(B, f)


def vetoed(E, P, src, X):
    v = np.zeros(len(E["weight"]), bool)
    for k in (1, 2):
        v |= (np.nan_to_num(P[f"{src}_{k}"], nan=0.0) >= 3.0) & (E[f"j{k}_rpt_itk"] < X)
    return v


print(f"core time: {args.core or 'standard'}")
for X in (float(x) for x in args.guards.split(",")):
    line = []
    for cl, thr in cuts:
        bS = np.ones(len(wS), bool) if thr is None else D["sig_score"] >= thr
        bB = np.ones(len(wB), bool) if thr is None else D["bkg_score"] >= thr
        for src in ("waves", "truth"):
            eS = wS[bS & ~vetoed(S, PS, src, X)].sum() / wS[bS].sum()
            eB = wB[bB & ~vetoed(B, PB, src, X)].sum() / wB[bB].sum()
            line.append(f"{src} {cl}: {100 * (eS / np.sqrt(eB) - 1):+5.1f}% ({eS:.3f}/{eB:.3f})")
    # what the guard does per leg (WAVeS, no BDT): HS legs in HGTD vetoed, PU legs in HGTD vetoed
    hs = [(S[f"j{k}_hs"] == 1) & fwd(S, k) for k in (1, 2)]
    pu = [(B[f"j{k}_pu"] == 1) & fwd(B, k) for k in (1, 2)]
    vh = np.concatenate([((np.nan_to_num(PS[f"waves_{k}"], nan=0) >= 3) & (S[f"j{k}_rpt_itk"] < X))[hs[k - 1]] for k in (1, 2)])
    vp = np.concatenate([((np.nan_to_num(PB[f"waves_{k}"], nan=0) >= 3) & (B[f"j{k}_rpt_itk"] < X))[pu[k - 1]] for k in (1, 2)])
    wp = np.concatenate([wB[pu[k - 1]] for k in (1, 2)])
    print(f"R_pT < {X:>4}: HS legs vetoed {100 * vh.mean():.2f}%, PU legs vetoed {100 * np.average(vp, weights=wp):.1f}% | " + "; ".join(line))
