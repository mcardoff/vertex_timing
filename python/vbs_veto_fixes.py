#!/usr/bin/env python3
"""Reco-only fixes to the WAVeS jet-vs-t0 veto, one per step of the money plot's chain, measured
against the veto as it is (JVT + fJVT loose, recoil > 100 GeV, jet time from >= 2 timed ghost tracks
within dR < 0.2, 3 sigma). Every fix uses only reconstructed quantities; truth only scores it.

  onetrack  out-of-time pileup (coverage): a jet time from ONE track is used too, at a tighter
            threshold k sigma (its sigma already carries f(1) = 2.3)
  fallback  out-of-time pileup (coverage): the dR < 0.2 time when it has >= 2 tracks, otherwise the
            dR < 0.4 time (>= 2 tracks; f(n) calibrated on the 0.4 cone's own HS legs)
  etacone   the same, but the 0.4 cone only for a jet axis well inside HGTD (2.5 < |eta| < 3.8)
  dominant  false vetoes: a jet time vetoes only when its time cluster holds >= a fraction f of the
            cone's timed tracks (one collision dominates the cone)
  tight     acceptance: on top of the timing veto, a tagging jet failing JVT tight (R_pT < 0.0767,
            |eta| < 2.5, pT < 60) or fJVT tight (fJVT > 0.4, 2.5 <= |eta| < 4.5, pT < 120) vetoes

For each: S/sqrt(B) change vs no veto (all events) and vs the WAVeS veto (paired bootstrap over the
same events); the signal and Z+jets kept; and, scored with truth, the share of Z+jets events with an
out-of-time pileup tagging jet that are vetoed, and the share of VBF hard-scatter pairs vetoed (false
vetoes). Needs the --pure-times trees (for the truth scoring) and the --jet-time-dr=0.4 trees.

    cd python && PYTHONNOUSERSITE=1 PYTHONPATH=/opt/homebrew/Cellar/root/6.40.04/lib/root \\
    ~/.venv-hgtd/bin/python vbs_veto_fixes.py
"""
import argparse, os
import numpy as np
import uproot
import ROOT
import vbs_time_veto_common as C

ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
ap.add_argument("--dir", default=os.path.join(os.path.dirname(os.path.abspath(__file__)), "../condor"))
ap.add_argument("--nsigma", type=float, default=3.0)
ap.add_argument("--nboot", type=int, default=300)
ap.add_argument("--out", help="figure stem: the fixes against each chain step's available gain")
args = ap.parse_args()
C.ROOT = ROOT
ROOT.gErrorIgnoreLevel = ROOT.kWarning


def load(p):
    ev = uproot.open(p)["events"]
    e = ev.arrays([k for k in ev.keys() if not k.startswith("jet_")], library="np")
    o = np.lexsort((e["entry"], e["file_path"]))
    return {k: v[o] for k, v in e.items()}


D = {}
for smp, rc in (("vbf", "met_truth"), ("zjets", "z_pt")):
    e2 = load(f"{args.dir}/{smp}/{smp}_jvtLoose_pure_vbs_time_veto.root")
    e4 = load(f"{args.dir}/{smp}/{smp}_jvtLoose_jtdr0p4_vbs_time_veto.root")
    assert np.array_equal(e2["entry"], e4["entry"]) and np.array_equal(e2["file_path"], e4["file_path"])
    D[smp] = (e2, e4, e2[rc] > 100)
cal2 = C.measure(D["vbf"][0], True, "HS legs, dR < 0.2")
cal4 = C.measure(D["vbf"][1], True, "HS legs, dR < 0.4")
_p = {"t0_truth": D["vbf"][0]["t0_truth"]}
for L in "ab":
    _p.update({f"{L}_hs": D["vbf"][0][f"{L}_hs"], f"{L}_core_n": D["vbf"][0][f"{L}_pure_hs_n"],
               f"{L}_core_t": D["vbf"][0][f"{L}_pure_hs_t"], f"{L}_core_sig": D["vbf"][0][f"{L}_pure_hs_sig"]})
pcal = C.measure(_p, True, "pure HS legs")
S2, S4, ks = D["vbf"]; B2, B4, kb = D["zjets"]
S2 = {k: v[ks] for k, v in S2.items()}; S4 = {k: v[ks] for k, v in S4.items()}
B2 = {k: v[kb] for k, v in B2.items()}; B4 = {k: v[kb] for k, v in B4.items()}
ws, wb = S2["weight"].astype(float), B2["weight"].astype(float)


def leg_fail(e, L, cal, t, sig, n, nmin, nsig):
    ok = e["ok_waves"].astype(bool)
    s0 = e["infl_waves"] * e["sig0_waves"]
    sg = cal(n) * sig
    with np.errstate(divide="ignore", invalid="ignore"):
        pull = np.abs(t - e["t0_waves"]) / np.sqrt(sg ** 2 + s0 ** 2)
    return ok & (n >= nmin) & (pull >= nsig)


def base_leg(e, L, nmin=2, nsig=None):
    return leg_fail(e, L, cal2, e[f"{L}_core_t"], e[f"{L}_core_sig"], e[f"{L}_core_n"], nmin, nsig or args.nsigma)


def waves(e, e4):
    return base_leg(e, "a") | base_leg(e, "b")


def onetrack(k):
    def f(e, e4):
        v = np.zeros(len(e["weight"]), bool)
        for L in "ab":
            n = e[f"{L}_core_n"]
            v |= base_leg(e, L) | ((n == 1) & base_leg(e, L, nmin=1, nsig=k))
        return v
    return f


def fallback(inside_only):
    def f(e, e4):
        v = np.zeros(len(e["weight"]), bool)
        for L in "ab":
            n2 = e[f"{L}_core_n"]
            use4 = n2 < 2
            if inside_only:
                ae = np.abs(e[f"{L}_eta"]); use4 &= (ae > 2.5) & (ae < 3.8)
            f4 = leg_fail(e4, L, cal4, e4[f"{L}_core_t"], e4[f"{L}_core_sig"], e4[f"{L}_core_n"], 2, args.nsigma)
            v |= base_leg(e, L) | (use4 & f4)
        return v
    return f


def dominant(frac):
    def f(e, e4):
        v = np.zeros(len(e["weight"]), bool)
        for L in "ab":
            n, nc = e[f"{L}_core_n"], e[f"{L}_core_ncand"]
            dom = n >= frac * np.maximum(nc, 1)
            v |= base_leg(e, L) & dom
        return v
    return f


def agree(e, e4):
    """false vetoes: a dR < 0.2 failure counts only if the dR < 0.4 time (when it has >= 2 tracks) fails too"""
    v = np.zeros(len(e["weight"]), bool)
    for L in "ab":
        n4 = e4[f"{L}_core_n"]
        f4 = leg_fail(e4, L, cal4, e4[f"{L}_core_t"], e4[f"{L}_core_sig"], n4, 2, args.nsigma)
        v |= base_leg(e, L) & ((n4 < 2) | f4)
    return v


def sumpt(frac):
    def f(e, e4):
        v = np.zeros(len(e["weight"]), bool)
        for L in "ab":
            v |= base_leg(e, L) & (e[f"{L}_core_sumpt"] >= frac * e[f"{L}_pt"])
        return v
    return f


def sumpt_ok(e, e4, frac):
    """the event is vetoed by a leg whose cluster carries >= frac of its pT: approximated per event as
    at least one tagging jet passing the sum-pT requirement (coverage-fallback legs included)"""
    return np.logical_or.reduce([e[f"{L}_core_sumpt"] >= frac * e[f"{L}_pt"] for L in "ab"] +
                                [e4[f"{L}_core_sumpt"] >= frac * e4[f"{L}_pt"] for L in "ab"])


def combo(*fns):
    return lambda e, e4: np.logical_or.reduce([fn(e, e4) for fn in fns])


def tight(e):
    v = np.zeros(len(e["weight"]), bool)
    for L in "ab":
        ae, pt = np.abs(e[f"{L}_eta"]), e[f"{L}_pt"]
        jvt_win = (ae < 2.5) & (pt < 60); fjvt_win = (ae >= 2.5) & (ae < 4.5) & (pt < 120)
        v |= (jvt_win & (e[f"{L}_jvt_rpt"] < 0.0767)) | (fjvt_win & (e[f"{L}_fjvt"] > 0.4))
    return v


def with_tight(fn):
    return lambda e, e4: fn(e, e4) | tight(e)


def classes(e):
    oot = np.zeros(len(e["weight"]), bool)
    for L in "ab":
        n = e[f"{L}_pure_pu_n"]; sg = pcal(n) * e[f"{L}_pure_pu_sig"]
        with np.errstate(divide="ignore", invalid="ignore"):
            pull = np.abs(e[f"{L}_pure_pu_t"] - e["t0_truth"]) / sg
        oot |= e[f"{L}_pu"] & (n >= 1) & (pull >= args.nsigma)
    reach_none = np.ones(len(e["weight"]), bool)
    for L in "ab":
        reach_none &= ~(e[f"{L}_pu"] & (e[f"{L}_hgtd_n"] >= 1))
    pu_unreach = (e["a_pu"] | e["b_pu"]) & reach_none
    itm = np.zeros(len(e["weight"]), bool)
    for L in "ab":
        n = e[f"{L}_pure_pu_n"]; sg = pcal(n) * e[f"{L}_pure_pu_sig"]
        with np.errstate(divide="ignore", invalid="ignore"):
            pull = np.abs(e[f"{L}_pure_pu_t"] - e["t0_truth"]) / sg
        itm |= e[f"{L}_pu"] & (n >= 1) & (pull < args.nsigma)
    return oot, e["a_hs"] & e["b_hs"], pu_unreach, itm & ~oot


OOT_B, HSHS_S, UNR_B, ITM_B = classes(B2)[0], classes(S2)[1], classes(B2)[2], classes(B2)[3]
rng = np.random.default_rng(11)
BS = [np.bincount(rng.integers(0, len(ws), len(ws)), minlength=len(ws)).astype(float) for _ in range(args.nboot)]
BB = [np.bincount(rng.integers(0, len(wb), len(wb)), minlength=len(wb)).astype(float) for _ in range(args.nboot)]


def z(vs, vb, bs=1.0, bb=1.0):
    a, b = ws * bs, wb * bb
    return 100 * ((a[~vs].sum() / a.sum()) / np.sqrt(b[~vb].sum() / b.sum()) - 1)


P = lambda w, m: 100 * w[m].sum() / w.sum()
v0s, v0b = waves(S2, S4), waves(B2, B4)
z0 = z(v0s, v0b); z0b = np.array([z(v0s, v0b, bs, bb) for bs, bb in zip(BS, BB)])
FIXES = [("WAVeS veto (as it is)", waves)] + \
        [(f"onetrack: 1-track times at >= {k:g} sigma", onetrack(k)) for k in (3, 4, 5, 6)] + \
        [("fallback: dR < 0.4 time when dR < 0.2 has < 2 tracks", fallback(False)),
         ("etacone: the same, axis 2.5 < |eta| < 3.8 only", fallback(True))] + \
        [(f"dominant: veto only if the cluster holds >= {int(100 * f)}% of the cone", dominant(f)) for f in (0.5, 0.67, 0.8)] + \
        [("agree: a dR < 0.2 failure needs the dR < 0.4 time to fail too", agree)] + \
        [(f"sumpt: veto only if the cluster carries >= {int(100 * f)}% of the jet pT", sumpt(f)) for f in (0.05, 0.1)] + \
        [("COMBINED: etacone + onetrack 4 sigma", combo(fallback(True), onetrack(4))),
         ("COMBINED + sumpt 5%", lambda e, e4: combo(fallback(True), onetrack(4))(e, e4) & sumpt_ok(e, e4, 0.05)),
         ("COMBINED + JVT/fJVT tight", with_tight(combo(fallback(True), onetrack(4)))),
         ("tight: + JVT/fJVT tight on the tagging jets", with_tight(waves)),
         ("JVT/fJVT tight alone (no timing)", lambda e, e4: tight(e))]
print(f"{'fix':58s} {'S/sqrtB':>8s} {'vs WAVeS':>14s} {'VBF kept':>9s} {'Z kept':>7s} {'OOT-PU Z vetoed':>16s} {'HS+HS VBF vetoed':>17s} {'unreach. PU Z vetoed':>21s} {'in-time PU Z vetoed':>20s}")
R = {}
for lab, fn in FIXES:
    vs, vb = fn(S2, S4), fn(B2, B4)
    zz = z(vs, vb); zb = np.array([z(vs, vb, bs, bb) for bs, bb in zip(BS, BB)])
    R[lab] = (zz, zz - z0, np.std(zb - z0b))
    print(f"{lab:58s} {zz:+7.1f}% {zz - z0:+6.1f} +- {np.std(zb - z0b):.1f} {100 - P(ws, vs):8.1f}% {100 - P(wb, vb):6.1f}% "
          f"{100 * wb[OOT_B & vb].sum() / wb[OOT_B].sum():15.0f}% {100 * ws[HSHS_S & vs].sum() / ws[HSHS_S].sum():16.1f}% "
          f"{100 * wb[UNR_B & vb].sum() / wb[UNR_B].sum():20.0f}% {100 * wb[ITM_B & vb].sum() / wb[ITM_B].sum():19.0f}%")

# ── the figure: each fix's gain over the WAVeS veto, against the gain its chain step makes available ──
if args.out:
    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    plt.rcParams.update({"font.size": 12})
    STEPS = [("+ perfect $t_0$", 0.9, "#A96B59", []),
             ("+ out-of-time\npileup caught", 8.4, "#1B7837",
              [("$\\Delta R < 0.4$ cone where the axis is inside HGTD", "etacone: the same, axis 2.5 < |eta| < 3.8 only"),
               ("one-track jet times, at $\\geq 4\\sigma$", "onetrack: 1-track times at >= 4 sigma"),
               ("both", "COMBINED: etacone + onetrack 4 sigma"),
               ("$\\Delta R < 0.4$ cone everywhere", "fallback: dR < 0.4 time when dR < 0.2 has < 2 tracks")]),
             ("+ in-time\npileup caught", 5.4, "#3F90DA",
              [("JVT + fJVT tight on the tagging jets (costs 13 points of signal)", "tight: + JVT/fJVT tight on the tagging jets")]),
             ("+ untimed\npileup caught", 0.2, "#7B3294", []),
             ("+ false vetoes\nundone", 1.7, "#717581",
              [("veto only if the time cluster has $\\geq 5\\%$ of the jet $p_T$", "sumpt: veto only if the cluster carries >= 5% of the jet pT"),
               ("veto only if 1 cluster holds $\\geq 50\\%$ of the cone", "dominant: veto only if the cluster holds >= 50% of the cone"),
               ("veto only if the $\\Delta R < 0.4$ time agrees", "agree: a dR < 0.2 failure needs the dR < 0.4 time to fail too")]),
             ("ceiling: pileup\nHGTD cannot reach", 7.0, "#000000",
              [("JVT + fJVT tight (as above: it reaches both)", "tight: + JVT/fJVT tight on the tagging jets")])]
    rows = []
    for name, avail, col, fixes in STEPS:
        rows.append(("step", name, avail, col))
        for lab, key in fixes: rows.append(("fix", lab, R[key][1], R[key][2]))
        if not fixes: rows.append(("none", "the vertex-time scores: not repeated here" if avail > 0.5 else "not worth one", 0, 0))
    rows.append(("step", "all of it\ntogether", None, "#444444"))
    rows.append(("fix", "out-of-time fixes + the sum-$p_T$ requirement", R["COMBINED + sumpt 5%"][1], R["COMBINED + sumpt 5%"][2]))
    rows.append(("fix", "out-of-time fixes + JVT/fJVT tight", R["COMBINED + JVT/fJVT tight"][1], R["COMBINED + JVT/fJVT tight"][2]))
    fig, ax = plt.subplots(figsize=(17.5, 9.0))
    y = 0; ys = []; labels = []
    cur_col = "#444444"
    for kind, lab, v, x in rows:
        if kind == "step":
            y -= 0.5
            cur_col = x
            if v is not None:
                ax.barh(y, v, height=0.62, color=cur_col, alpha=0.22, edgecolor=cur_col, lw=1.2)
                ax.text(v + 0.12, y, f"available: +{v:.1f}", va="center", fontsize=13, color=cur_col, fontweight="bold")
            ax.text(-0.25, y, lab, va="center", ha="right", fontsize=14, color=cur_col, fontweight="bold", transform=ax.transData)
        elif kind == "fix":
            ax.errorbar(v, y, xerr=x, fmt="o", color=cur_col, ms=7, capsize=3, lw=1.6)
            ax.text(max(v + x, 0) + 0.15, y, f"{lab}:  {v:+.1f} $\\pm$ {x:.1f}", va="center", fontsize=13)
        else:
            ax.text(0.15, y, lab, va="center", fontsize=13, color="#888888", style="italic")
        ys.append(y); y -= 1
    ax.axvline(0, color="#888888", lw=0.9)
    ax.set_yticks([]); ax.set_ylim(y + 0.3, 0.3)
    ax.set_xlim(-4.6, 11.5)
    ax.set_xlabel("$S/\\sqrt{B}$ change from the WAVeS veto [% points]", fontsize=14)
    ax.tick_params(labelsize=13)
    ax.grid(alpha=0.25, axis="x")
    for sp in ("left", "right", "top"): ax.spines[sp].set_visible(False)
    fig.text(0.5, -0.02, "all selected events, recoil > 100 GeV, JVT + fJVT loose, WAVeS $t_0$, 3$\\sigma$; every fix uses reconstructed quantities only   "
             "ATLAS Simulation Internal", ha="center", fontsize=9.5, color="#555555")
    for ext in ("png", "pdf"): fig.savefig(f"{args.out}.{ext}", dpi=150, bbox_inches="tight")
    print("wrote", args.out + ".png")
