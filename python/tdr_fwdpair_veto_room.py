"""Where the analysis veto's room is (results/veto_room.md, 2026-11-14): the jet-vs-t0 timing veto against the perfect
HGTD veto (truth labels), on tdr_fwdpair trees (--pair=one --jvt=jvt4pt60), before the BDT, weighted as
tdr_fwdpair_kinplots.py --veto (whose S/sqrt(B) factors it reproduces: WAVeS +26%, truth +32%, perfect +48%).
  - the gap split into signal loss (hard-scatter legs out of time) and background missed (pileup legs untimed / in time);
  - counterfactuals: hard-scatter legs never veto; better jet-time resolution (pileup-leg pulls / r); a leg vetoes only
    with >= m core tracks;
  - the hard-scatter legs' out-of-time rate by core track count.

    cd python && PYTHONNOUSERSITE=1 ~/.venv-hgtd/bin/python tdr_fwdpair_veto_room.py \\
        ../condor/vbf/vbf_onefwd_jvt4pt60_wfix_tdr_fwdpair.root ../condor/zjets/zjets_onefwd_jvt4pt60_wfix_tdr_fwdpair.root
"""
import argparse

import numpy as np
import uproot

import tdr_fwdpair_common as FC

ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
ap.add_argument("sig", help="tdr_fwdpair tree, VBF H->inv (also the jet-time calibration)")
ap.add_argument("bkg", help="tdr_fwdpair tree, Z+jets")
args = ap.parse_args()
S = uproot.open(args.sig)["events"].arrays(library="np")
B = uproot.open(args.bkg)["events"].arrays(library="np")
f = FC.calibration(S, verbose=False)
fwd = lambda E, k: np.abs(E[f"j{k}_eta"]) > 2.4
pu_hgtd = lambda E: ((E["j1_pu"] == 1) & fwd(E, 1)) | ((E["j2_pu"] == 1) & fwd(E, 2))


def eff(E, kept):
    w = E["weight"].astype(float)
    return w[kept].sum() / w.sum()


print("S/sqrt(B) factor, before the BDT (eps_S, eps_B)")
for src in ("waves", "truth"):
    kS, kB = FC.vetoes(S, f)[src], FC.vetoes(B, f)[src]
    print(f"  jet vs {src:6s} t0: eps_S {eff(S, kS):.3f}  eps_B {eff(B, kB):.3f}  -> x{eff(S, kS) / np.sqrt(eff(B, kB)):.3f}")
kB, kS = ~pu_hgtd(B), ~pu_hgtd(S)   # as tdr_fwdpair_kinplots: the perfect veto removes PU-in-HGTD events in BOTH samples
print(f"  perfect veto      : eps_S {eff(S, kS):.3f}  eps_B {eff(B, kB):.3f}  -> x{eff(S, kS) / np.sqrt(eff(B, kB)):.3f}")
# counterfactual: hard-scatter legs never veto (perfect HS jet times), per t0 source
for src in ("waves", "truth"):
    out = {}
    for name, E in (("S", S), ("B", B)):
        P = FC.pulls(E, f)
        v = np.zeros(len(E["weight"]), bool)
        for k in (1, 2):
            v |= (E[f"j{k}_hs"] != 1) & (np.nan_to_num(P[f"{src}_{k}"], nan=0) >= 3)
        out[name] = eff(E, ~v)
    print(f"  jet vs {src:6s} t0, HS legs never veto: eps_S {out['S']:.3f}  eps_B {out['B']:.3f}  -> x{out['S'] / np.sqrt(out['B']):.3f}")


def legs(E, src):
    """per leg k: timed (core_n >= 1), pull against the t0 source (NaN untimed)"""
    P = FC.pulls(E, f)
    return {k: (E[f"j{k}_core_n"] >= 1, P[f"{src}_{k}"]) for k in (1, 2)}


for src in ("truth", "waves"):
    print(f"\n== jet vs {src} t0 ==")
    # background: the events the perfect veto removes (a PU leg in HGTD), and what the timing veto does with them
    L = legs(B, src)
    w = B["weight"].astype(float)
    ph = pu_hgtd(B)
    veto = ~FC.vetoes(B, f)[src]
    # classify each event by its PU legs in HGTD: caught (some PU-HGTD leg pull >= 3), else untimed (every PU-HGTD leg
    # untimed), else in time (a timed PU-HGTD leg with pull < 3)
    caught = np.zeros(len(w), bool); untimed_all = ph.copy(); intime = np.zeros(len(w), bool)
    for k in (1, 2):
        isph = (B[f"j{k}_pu"] == 1) & fwd(B, k)
        t, p = L[k]
        caught |= isph & t & (np.nan_to_num(p, nan=0) >= 3)
        untimed_all &= ~(isph & t)
        intime |= isph & t & (np.nan_to_num(p, nan=0) < 3)
    W = w[ph].sum()
    print(f"Z+jets events with a PU leg in HGTD: {ph.sum()} ({100 * w[ph].sum() / w.sum():.1f}% of the weight)")
    print(f"  vetoed by timing: {100 * w[ph & veto].sum() / W:.1f}%   (a PU-HGTD leg out of time: {100 * w[ph & caught].sum() / W:.1f}%)")
    print(f"  missed: every PU-HGTD leg untimed {100 * w[ph & ~veto & untimed_all].sum() / W:.1f}%, "
          f"a timed PU-HGTD leg in time {100 * w[ph & ~veto & ~untimed_all].sum() / W:.1f}%")
    other = ~ph
    print(f"Z+jets events with NO PU leg in HGTD: {other.sum()} ({100 * w[other].sum() / w.sum():.1f}%), vetoed anyway "
          f"{100 * w[other & veto].sum() / w[other].sum():.1f}% (the perfect veto keeps all)")
    # the PU legs in HGTD that are timed: their pull distribution
    pulls = np.concatenate([L[k][1][((B[f"j{k}_pu"] == 1) & fwd(B, k)) & L[k][0]] for k in (1, 2)])
    ws = np.concatenate([w[((B[f"j{k}_pu"] == 1) & fwd(B, k)) & L[k][0]] for k in (1, 2)])
    dt = np.concatenate([(B[f"j{k}_core_t"] - B["t0_truth"])[((B[f"j{k}_pu"] == 1) & fwd(B, k)) & L[k][0]] for k in (1, 2)])
    n = np.concatenate([B[f"j{k}_core_n"][((B[f"j{k}_pu"] == 1) & fwd(B, k)) & L[k][0]] for k in (1, 2)])
    print(f"  timed PU legs in HGTD: {len(pulls)}; pull < 3 for {100 * np.average(pulls < 3, weights=ws):.1f}%; "
          f"|t_core - t_HS| < 100 ps for {100 * np.average(np.abs(dt) < 100, weights=ws):.1f}%, median |dt| {np.median(np.abs(dt)):.0f} ps; "
          f"core tracks median {np.median(n):.0f}, 1 track {100 * np.mean(n == 1):.0f}%")
    # signal: what the veto costs
    LS = legs(S, src)
    wS = S["weight"].astype(float)
    vS = ~FC.vetoes(S, f)[src]
    print(f"VBF: vetoed {100 * wS[vS].sum() / wS.sum():.1f}%; of those, by an HS leg {100 * wS[vS & ((S['j1_hs'] == 1) & (np.nan_to_num(LS[1][1], nan=0) >= 3) | (S['j2_hs'] == 1) & (np.nan_to_num(LS[2][1], nan=0) >= 3))].sum() / max(wS[vS].sum(), 1e-12):.0f}%")
    sel = np.ones(len(wS), bool)
    nleg = sum(((S[f"j{k}_hs"] == 1) & fwd(S, k)).sum() for k in (1, 2))
    nt = sum(((S[f"j{k}_hs"] == 1) & fwd(S, k) & LS[k][0]).sum() for k in (1, 2))
    print(f"  HS legs in HGTD timed: {100 * nt / nleg:.1f}% of {nleg}")

print("\n== jet-time resolution counterfactual: PU-leg pulls / r (their dt is the collision offset), HS legs unchanged ==")
for src in ("waves", "truth"):
    for r in (1.0, 0.75, 0.5):
        out = {}
        for name, E in (("S", S), ("B", B)):
            P = FC.pulls(E, f)
            v = np.zeros(len(E["weight"]), bool)
            for k in (1, 2):
                p = np.nan_to_num(P[f"{src}_{k}"], nan=0)
                v |= np.where(E[f"j{k}_pu"] == 1, p / r, p) >= 3
            out[name] = eff(E, ~v)
        print(f"  jet vs {src:6s} t0, resolution x{r:.2f}: eps_S {out['S']:.3f}  eps_B {out['B']:.3f}  -> x{out['S'] / np.sqrt(out['B']):.3f}")
# the veto window for timed PU legs in HGTD: 3 x the pull's denominator
for src in ("waves", "truth"):
    P = FC.pulls(B, f)
    with np.errstate(invalid="ignore", divide="ignore"):   # a zero pull (|dt| = 0) gives 0/0, dropped below
        den = np.concatenate([np.abs(B[f"j{k}_core_t"] - B[f"t0_{src}"])[(B[f"j{k}_pu"] == 1) & fwd(B, k) & (B[f"j{k}_core_n"] >= 1)]
                              / P[f"{src}_{k}"][(B[f"j{k}_pu"] == 1) & fwd(B, k) & (B[f"j{k}_core_n"] >= 1)] for k in (1, 2)])
    den = den[np.isfinite(den)]
    print(f"  {src}: 3-sigma window for timed PU legs: median {3 * np.median(den):.0f} ps (16-84%: {3 * np.percentile(den, 16):.0f}-{3 * np.percentile(den, 84):.0f})")

print("\n== HS legs in HGTD out of time (pull >= 3) by core track count, VBF ==")
for src in ("waves", "truth"):
    P = FC.pulls(S, f)
    n_all, p_all = [], []
    for k in (1, 2):
        m = (S[f"j{k}_hs"] == 1) & fwd(S, k) & (S[f"j{k}_core_n"] >= 1)
        n_all.append(S[f"j{k}_core_n"][m]); p_all.append(np.nan_to_num(P[f"{src}_{k}"][m], nan=0))
    n_all, p_all = np.concatenate(n_all), np.concatenate(p_all)
    out = []
    for lo, hi in ((1, 1), (2, 2), (3, 4), (5, 10**6)):
        m = (n_all >= lo) & (n_all <= hi)
        out.append(f"n={lo}{'' if hi == lo else ('+' if hi > 100 else f'-{hi}')}: {100 * np.mean(m):.0f}% of legs, {100 * np.mean(p_all[m] >= 3):.1f}% out")
    print(f"  {src}: all {100 * np.mean(p_all >= 3):.1f}% out; " + "; ".join(out))

print("\n== reco-only variant: a leg vetoes only with >= m core tracks ==")
for src in ("waves", "truth"):
    for mmin in (1, 2, 3):
        out = {}
        for name, E in (("S", S), ("B", B)):
            P = FC.pulls(E, f)
            v = np.zeros(len(E["weight"]), bool)
            for k in (1, 2):
                v |= (E[f"j{k}_core_n"] >= mmin) & (np.nan_to_num(P[f"{src}_{k}"], nan=0) >= 3)
            out[name] = eff(E, ~v)
        print(f"  jet vs {src:6s} t0, core n >= {mmin}: eps_S {out['S']:.3f}  eps_B {out['B']:.3f}  -> x{out['S'] / np.sqrt(out['B']):.3f}")
