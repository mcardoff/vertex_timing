#!/usr/bin/env python3
"""
qoverp_term.py -- can cluster_qOverP_sigma improve the selection?

cluster_qOverP_sigma = 1/sqrt(SUM_t 1/sigma_qOverP_t^2): the inverse-variance-
combined momentum uncertainty of the cluster (export_training_data.cxx:1131).
Per-track sigma_qOverP falls with momentum, so the combined value falls with both
track COUNT and track MOMENTUM -- i.e. it plausibly restates n_tracks and sumpt,
which the score already carries. Three questions, whole population, both samples:

  1. Redundancy: correlation of log(sigma_qP) with log(n_tracks), log(sumpt) and
     the score itself, plus its tie-break accuracy CONDITIONAL on the pair having
     equal n_tracks (does it know anything track count does not?).
  2. Switch accounting (the established protocol): re-rank ambiguous top-2 pairs
     by it, one global direction, count better vs worse.
  3. Multiplicative integration: score' = tz * sigma_qP^(-beta), one global beta,
     core fraction with the guarded timing on both samples.
"""
import argparse, sys, os
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from tiebreak_scan import load, auc, toppairs, EVT

PASS = 60.0
V = "cluster_qOverP_sigma"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--vbf", default="figs/hists/vbf_mjj500p0_training.root")
    ap.add_argument("--zjets",
                    default="/Users/mcard/project/vertex_timing/training_new/zjets_novbs_training.root")
    a = ap.parse_args()
    S = {}
    for tag, path in (("vbf", a.vbf), ("zjets", a.zjets)):
        print(f"  loading {tag} ...", flush=True)
        S[tag] = load(path)

    # ---- 1. redundancy ----------------------------------------------------
    print(f"\n  1. REDUNDANCY -- what does {V} know that the score does not?")
    for tag in ("vbf", "zjets"):
        c, _ = S[tag]
        multi = c.groupby(EVT, sort=False)["cluster_idx"].transform("size") > 1
        cm = c[multi]
        m = np.isfinite(cm[V]) & (cm[V] > 0) & (cm["tz"] > 0)
        lx = np.log(cm.loc[m, V].to_numpy(float))
        corr = {q: np.corrcoef(lx, np.log(np.maximum(cm.loc[m, q].to_numpy(float), 1e-9)))[0, 1]
                for q in ("n_tracks", "sumpt", "tz")}
        print(f"     {tag:6s} corr(log sigma_qP, log n_trk) = {corr['n_tracks']:+.3f}   "
              f"log sumpt = {corr['sumpt']:+.3f}   log score = {corr['tz']:+.3f}")
        aa, bb = toppairs(c)
        dom0 = aa["is_dom"].to_numpy()
        va = np.where(dom0, aa[V].to_numpy(float), bb[V].to_numpy(float))
        vf = np.where(dom0, bb[V].to_numpy(float), aa[V].to_numpy(float))
        eqn = aa["n_tracks"].to_numpy() == bb["n_tracks"].to_numpy()
        ok = np.isfinite(va) & np.isfinite(vf)
        # direction: smaller sigma marks HS, so accuracy = P(v_true < v_false)
        acc_all = np.mean(va[ok] < vf[ok]) + 0.5 * np.mean(va[ok] == vf[ok])
        k = ok & eqn
        acc_eq = (np.mean(va[k] < vf[k]) + 0.5 * np.mean(va[k] == vf[k])) if k.sum() else np.nan
        print(f"            tie-break acc: all pairs {acc_all:.3f}   "
              f"pairs with EQUAL n_tracks ({int(k.sum()):,}) {acc_eq:.3f}")

    # ---- 2. switch accounting --------------------------------------------
    print(f"\n  2. SWITCH ACCOUNTING -- re-rank ambiguous pair by smaller {V}")
    print(f"     {'f':>4s} " + " ".join(f"{t + ' ' + k:>13s}" for t in ("vbf", "zjets")
                                        for k in ("net", "bet:wor")))
    for f in (0.5, 0.8):
        line = f"     {f:4.1f} "
        for tag in ("vbf", "zjets"):
            c, truth = S[tag]
            t2 = c.sort_values(EVT + ["tz"], ascending=[True, True, True, False]).copy()
            t2["rk"] = t2.groupby(EVT, sort=False).cumcount()
            t2 = t2[t2["rk"] < 2]
            n2 = t2.groupby(EVT, sort=False)["rk"].transform("size") == 2
            a1 = t2[(t2["rk"] == 0) & n2].set_index(
                pd.MultiIndex.from_frame(t2[(t2["rk"] == 0) & n2][EVT]))
            b1 = t2[(t2["rk"] == 1) & n2].set_index(
                pd.MultiIndex.from_frame(t2[(t2["rk"] == 1) & n2][EVT]))
            b1 = b1.reindex(a1.index)
            y = truth.reindex(a1.index).to_numpy()
            amb = b1["tz"].to_numpy() > f * a1["tz"].to_numpy()
            va, vb = a1[V].to_numpy(float), b1[V].to_numpy(float)
            swap = amb & np.isfinite(va) & np.isfinite(vb) & (vb < va)
            tA, tB = a1["t_out"].to_numpy(), b1["t_out"].to_numpy()
            ok0 = np.abs(tA - y) < PASS
            ok1 = np.abs(np.where(swap, tB, tA) - y) < PASS
            bet = int((swap & ok1 & ~ok0).sum()); wor = int((swap & ~ok1 & ok0).sum())
            line += f" {100.0*(bet-wor)/len(truth):+12.2f} {bet:6d}:{wor:<6d}"
        print(line)

    # ---- 3. multiplicative term ------------------------------------------
    print(f"\n  3. MULTIPLICATIVE -- score' = tz * {V}^(-beta), one global beta")
    print(f"     {'beta':>6s} " + " ".join(f"{t:>9s}" for t in ("vbf", "zjets")) + "   worst d")
    base = {}
    for beta in (0.0, 0.05, 0.1, 0.2, 0.4, 0.8, 1.6):
        vals = {}
        for tag in ("vbf", "zjets"):
            c, truth = S[tag]
            sig = np.maximum(c[V].to_numpy(float), 1e-9)
            c["_sc"] = c["tz"].to_numpy() * np.power(sig, -beta)
            ip = c.groupby(EVT, sort=False)["_sc"].idxmax()
            p = c.loc[ip].set_index(pd.MultiIndex.from_frame(c.loc[ip][EVT]))
            p.index.names = EVT
            vals[tag] = 100.0 * (np.abs(p["t_out"].reindex(truth.index).to_numpy()
                                        - truth.to_numpy()) < PASS).mean()
        if beta == 0.0: base = dict(vals)
        d = [vals[t] - base[t] for t in ("vbf", "zjets")]
        print(f"     {beta:6.2f} " + " ".join(f"{vals[t]:8.2f}%" for t in ("vbf", "zjets"))
              + f"   {min(d):+6.2f}" + ("   safe on both" if min(d) >= -0.05 else ""))


if __name__ == "__main__":
    main()
