"""How much of the timing veto's gain lives in jets JVT covers (pT < 60 GeV)? (results/veto_ptsplit.md, 2026-11-14)

The user: fold timing into the JVT option -- it acts on jets below 60 GeV only, "so we are not missing too much of our
gain". Tested here on tdr_fwdpair trees (--pair=one --jvt=jvt4pt60) with the analysis' kinematic BDT, exactly as
tdr_fwdpair_kinplots.py --veto (its S/sqrt(B) factors reproduce): the jet-vs-t0 veto with only the legs in a pT range
allowed to veto, against the veto on every leg; the perfect HGTD veto likewise restricted. Also where the pileup legs
the veto has to find sit in pT, after JVT.

    cd python && PYTHONNOUSERSITE=1 ~/.venv-hgtd/bin/python tdr_fwdpair_veto_ptsplit.py \\
        ../condor/vbf/vbf_onefwd_jvt4pt60_wfix_tdr_fwdpair.root ../condor/zjets/zjets_onefwd_jvt4pt60_wfix_tdr_fwdpair.root \\
        ../condor/zjets/vbf_vs_zjets_onefwd_jvt4pt60_kinbdt.npz
"""
import argparse

import numpy as np
import uproot

import tdr_fwdpair_common as FC

ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
ap.add_argument("sig")
ap.add_argument("bkg")
ap.add_argument("bdt", help="the kinematic BDT's scores, matched to the trees by row order")
ap.add_argument("--effs", default="0.8", help="BDT signal efficiencies besides no BDT")
args = ap.parse_args()
S = uproot.open(args.sig)["events"].arrays(library="np")
B = uproot.open(args.bkg)["events"].arrays(library="np")
D = np.load(args.bdt)
assert len(D["sig_score"]) == len(S["weight"]) and len(D["bkg_score"]) == len(B["weight"]), "BDT scores do not match the trees"
f = FC.calibration(S, verbose=False)
wS, wB = S["weight"].astype(float), B["weight"].astype(float)
o = np.argsort(-D["sig_score"]); cs = np.cumsum(wS[o]) / wS.sum()   # tdr_fwdpair_kinplots' BDT cut
cuts = [(None, None)] + [(e, D["sig_score"][o][np.searchsorted(cs, e)]) for e in (float(x) for x in args.effs.split(","))]
fwd = lambda E, k: np.abs(E[f"j{k}_eta"]) > 2.4
RANGES = [("every leg", 0, 1e9), ("legs < 60 GeV (JVT's range)", 0, 60), ("legs < 50 GeV", 0, 50), ("legs >= 60 GeV", 60, 1e9)]


def veto(E, src, lo, hi):
    """vetoed per event: a leg with lo <= pT < hi, timed, >= 3 sigma from the t0 (tdr_fwdpair_common.vetoes otherwise)"""
    P = FC.pulls(E, f)
    v = np.zeros(len(E["weight"]), bool)
    for k in (1, 2):
        v |= (E[f"j{k}_pt"] >= lo) & (E[f"j{k}_pt"] < hi) & (np.nan_to_num(P[f"{src}_{k}"], nan=0.0) >= 3.0)
    return v


def perfect(E, lo, hi):
    """the perfect HGTD veto restricted to pileup legs in the pT range (both samples, as kinplots)"""
    v = np.zeros(len(E["weight"]), bool)
    for k in (1, 2):
        v |= (E[f"j{k}_pu"] == 1) & fwd(E, k) & (E[f"j{k}_pt"] >= lo) & (E[f"j{k}_pt"] < hi)
    return v


for eff, thr in cuts:
    bS = np.ones(len(wS), bool) if thr is None else D["sig_score"] >= thr
    bB = np.ones(len(wB), bool) if thr is None else D["bkg_score"] >= thr
    print(f"\n== {'before the BDT' if thr is None else f'BDT at eps_S = {eff:.0%}'}: {bS.sum()} VBF / {bB.sum()} Z+jets MC events ==")
    # where the pileup legs in HGTD sit in pT, after JVT (and the BDT)
    ph = [(B[f"j{k}_pu"] == 1) & fwd(B, k) & bB for k in (1, 2)]
    pts = np.concatenate([B[f"j{k}_pt"][ph[k - 1]] for k in (1, 2)])
    wl = np.concatenate([wB[ph[k - 1]] for k in (1, 2)])
    print(f"  Z+jets pileup legs in HGTD: {len(pts)} MC legs; pT < 50 / 50-60 / >= 60 GeV: "
          + " / ".join(f"{100 * wl[(pts >= a) & (pts < b)].sum() / wl.sum():.0f}%" for a, b in ((0, 50), (50, 60), (60, 1e9))))
    print(f"  {'veto':38s} {'eps_S':>6s} {'eps_B':>6s} {'S/sqrt(B)':>10s} {'share of the gain':>18s}")
    for src in ("waves", "truth", "perfect"):
        full = None
        for name, lo, hi in RANGES:
            vS = perfect(S, lo, hi) if src == "perfect" else veto(S, src, lo, hi)
            vB = perfect(B, lo, hi) if src == "perfect" else veto(B, src, lo, hi)
            eS = wS[bS & ~vS].sum() / wS[bS].sum()
            eB = wB[bB & ~vB].sum() / wB[bB].sum()
            fac = eS / np.sqrt(eB)
            if full is None: full = fac
            lab = {"waves": "jet vs WAVeS t0", "truth": "jet vs truth t0", "perfect": "perfect HGTD veto"}[src]
            print(f"  {lab + ', ' + name:38s} {eS:6.3f} {eB:6.3f} {fac:9.3f}  {100 * (fac - 1) / (full - 1):16.0f}%")
