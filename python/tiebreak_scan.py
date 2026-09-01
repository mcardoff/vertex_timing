#!/usr/bin/env python3
"""
tiebreak_scan.py -- what ELSE in the reco record could separate the true cluster
from the one we wrongly pick, tested on the WHOLE population from the start?

The jet-fraction lesson: a variable ranked on the loss subset overstates itself
(0.90 there, 0.612 population-wide on Z+jets). So this scan never conditions on
losing. Three stages, each population-wide and run on BOTH samples:

  1. Population AUC: over every cluster in every multi-cluster event, how well
     does each reco column identify the HS-dominant cluster?
  2. Tie-break accuracy: restrict to the reco-defined AMBIGUOUS pairs -- the
     top-2 clusters by our score, in events where the HS-dominant cluster is one
     of them -- and ask how often each variable ranks the true one first. The
     conditioning is on the score's own top-2, which is reco-only, so this is
     not the loss-subset trap. The bar to beat is the score's own accuracy on
     the same pairs.
  3. Switch accounting for the survivors: when the score is ambiguous
     (score_2nd > f * score_1st), re-rank the pair by the candidate variable
     (ONE global direction, same threshold everywhere) and count better vs
     worse on core fraction. A candidate is worth developing only if the ratio
     is decisively > 1 on BOTH samples.

Also computes one column the export lacks: |t_cluster - t0_classical|, the
agreement with the model-free aggregate -- the classical analogue of the ML
selector's top feature.
"""
import argparse, gc, sys, os
import numpy as np
import pandas as pd
import uproot

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from classical_t0 import wmedian, winsor

EVT = ["sample_id", "file_idx", "event_num"]
KEY = EVT + ["cluster_idx"]
A_DZ, GUARD, PASS = 0.7, 3, 60.0
SKIP = set(KEY + ["weight", "delta_t", "t_truth", "dzcl", "s", "tz", "t_out",
                  "t_inj", "n_inj", "is_dom", "hgtd_valid"])


def auc(x, y):
    m = np.isfinite(x); x, y = x[m], y[m].astype(bool)
    if y.sum() == 0 or (~y).sum() == 0: return np.nan
    r = pd.Series(x).rank().to_numpy()
    return (r[y].sum() - y.sum() * (y.sum() + 1) / 2.0) / (y.sum() * (~y).sum())


def load(path):
    c = uproot.open(path)["clusters"].arrays(library="pd")
    for q in c.columns:
        if c[q].dtype == np.float64: c[q] = c[q].astype(np.float32)
    t = uproot.open(path)["tracks"].arrays(
        KEY + ["pt", "time", "timeRes", "z0_pull_pv", "sigma_z0",
               "dr_nearest_fwdjet"], library="pd")
    t = t.dropna(subset=["time", "timeRes", "z0_pull_pv", "sigma_z0"])
    t = t[t["timeRes"] > 0]
    t["dz"] = t["z0_pull_pv"].abs() * t["sigma_z0"]
    t["w"] = t["pt"] * np.exp(-A_DZ * t["dz"])
    t["iv"] = 1.0 / t["timeRes"] ** 2
    t["wiv"] = t["w"] * t["iv"]
    t["inj"] = (t["dr_nearest_fwdjet"] < 0.4).fillna(False)

    c["t_truth"] = c["cluster_time"] - c["delta_t"]
    truth = c.groupby(EVT, sort=False)["t_truth"].first(); truth.index.names = EVT
    mi = pd.MultiIndex.from_frame(c[KEY])
    g = t[t["inj"]].groupby(KEY, sort=False)
    tij = g.apply(lambda x: (x["iv"] * x["time"]).sum() / x["iv"].sum(),
                  include_groups=False).reindex(mi).to_numpy()
    c["t_inj"] = np.where(np.isnan(tij), c["cluster_time"], tij)
    c["n_inj"] = g.size().reindex(mi).fillna(0).to_numpy()
    s = t.groupby(KEY, sort=False)["w"].sum().rename("s").reset_index()
    c = c.merge(s, on=KEY, how="left")
    c["tz"] = c["s"].fillna(0.0) * np.exp(-1.5 * c["delta_z"].abs()).astype(np.float32)
    c["t_out"] = np.where(c["n_inj"] >= GUARD, c["t_inj"], c["cluster_time"])
    # the one computed column the export lacks: agreement with the classical
    # aggregate t0 (weighted median + double truncation, model-free)
    t0 = winsor(t, winsor(t, wmedian(t, "w", "time"), 100.0, "wiv"), 60.0, "wiv")
    c["dt_to_agg_cl"] = np.abs(
        c["cluster_time"].to_numpy()
        - t0.reindex(pd.MultiIndex.from_frame(c[EVT])).to_numpy())
    del t; gc.collect()

    hb = c.sort_values(["truth_n_hs_tracks", "n_tracks"]).groupby(EVT, sort=False).tail(1)
    c["is_dom"] = c.index.isin(hb.index)
    return c, truth


def toppairs(c):
    """Top-2 clusters by tz per multi-cluster event, as an aligned wide frame,
    restricted to events where the HS-dominant cluster is one of the two."""
    t2 = c.sort_values(EVT + ["tz"], ascending=[True, True, True, False]).copy()
    t2["rk"] = t2.groupby(EVT, sort=False).cumcount()
    t2 = t2[t2["rk"] < 2]
    n = t2.groupby(EVT, sort=False)["rk"].transform("size")
    nd = t2.groupby(EVT, sort=False)["is_dom"].transform("sum")
    t2 = t2[(n == 2) & (nd == 1)]
    a = t2[t2["rk"] == 0].set_index(pd.MultiIndex.from_frame(t2[t2["rk"] == 0][EVT]))
    b = t2[t2["rk"] == 1].set_index(pd.MultiIndex.from_frame(t2[t2["rk"] == 1][EVT]))
    b = b.reindex(a.index)
    return a, b


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

    cand = [q for q in S["vbf"][0].columns
            if q not in SKIP and not q.startswith("truth_")
            and S["vbf"][0][q].dtype != object and q in S["zjets"][0].columns]

    # ---- stages 1+2 -------------------------------------------------------
    rows = []
    P = {}
    for tag in ("vbf", "zjets"):
        c, _ = S[tag]
        multi = c.groupby(EVT, sort=False)["cluster_idx"].transform("size") > 1
        cm = c[multi]
        lab = cm["is_dom"].to_numpy()
        aa, bb = toppairs(c)
        P[tag] = (aa, bb)
        dom0 = aa["is_dom"].to_numpy()
        base_acc = dom0.mean()          # score's own tie-break accuracy
        print(f"\n  {tag}: {int(multi.sum()):,} clusters in multi-cluster events; "
              f"{len(aa):,} ambiguous top-2 pairs; score's own accuracy {100*base_acc:.2f}%")
        for q in cand:
            pa = auc(cm[q].to_numpy(float), lab)
            va = np.where(dom0, aa[q].to_numpy(float), bb[q].to_numpy(float))
            vf = np.where(dom0, bb[q].to_numpy(float), aa[q].to_numpy(float))
            ok = np.isfinite(va) & np.isfinite(vf)
            acc = (np.mean(va[ok] > vf[ok]) + 0.5 * np.mean(va[ok] == vf[ok])) if ok.sum() else np.nan
            rows.append((q, tag, pa, acc))
    R = pd.DataFrame(rows, columns=["var", "smp", "pop", "tie"]).pivot(
        index="var", columns="smp", values=["pop", "tie"])
    R.columns = [f"{x}_{y}" for x, y in R.columns]
    # one global direction from the pooled population AUC; a variable whose
    # direction flips between samples is unusable by construction
    R["dir"] = np.sign(R[["pop_vbf", "pop_zjets"]].mean(axis=1) - 0.5)
    for col in ("pop_vbf", "pop_zjets", "tie_vbf", "tie_zjets"):
        R["f" + col] = 0.5 + R["dir"] * (R[col] - 0.5)
    R["worst_tie"] = R[["ftie_vbf", "ftie_zjets"]].min(axis=1)
    R = R.sort_values("worst_tie", ascending=False)

    print(f"\n  ===== ranked by WORST-sample tie-break accuracy (one global direction)")
    print(f"  {'variable':<28s} {'pop vbf':>8s} {'pop zj':>7s} {'tie vbf':>8s} "
          f"{'tie zj':>7s} {'worst':>6s}")
    for v, r in R.head(18).iterrows():
        print(f"  {v:<28s} {r['fpop_vbf']:8.3f} {r['fpop_zjets']:7.3f} "
              f"{r['ftie_vbf']:8.3f} {r['ftie_zjets']:7.3f} {r['worst_tie']:6.3f}")

    # ---- stage 3: switch accounting for the survivors ---------------------
    top = [v for v in R.index[:6] if v not in ("trkptz_score", "sumpt", "sumpt2",
                                               "trkptz_rank", "trkptz_ratio_to_max",
                                               "is_max_trkptz", "is_max_sumpt",
                                               "sumpt_rank", "sumpt_ratio_to_max",
                                               "sumpt_frac_of_event")][:3]
    print(f"\n  ===== switch accounting: when tz_2nd > f*tz_1st, re-rank the pair by v")
    print(f"  {'variable':<24s} {'f':>4s} " +
          " ".join(f"{t + ' ' + k:>12s}" for t in ("vbf", "zjets")
                   for k in ("net", "bet:wor")))
    for v in top:
        d = float(R.loc[v, "dir"])
        for f in (0.5, 0.8):
            line = f"  {v:<24s} {f:4.1f} "
            for tag in ("vbf", "zjets"):
                c, truth = S[tag]
                aa, bb = None, None
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
                amb = (b1["tz"].to_numpy() > f * a1["tz"].to_numpy())
                va, vb = d * a1[v].to_numpy(float), d * b1[v].to_numpy(float)
                swap = amb & np.isfinite(va) & np.isfinite(vb) & (vb > va)
                tA = a1["t_out"].to_numpy(); tB = b1["t_out"].to_numpy()
                ok0 = np.abs(tA - y) < PASS
                ok1 = np.abs(np.where(swap, tB, tA) - y) < PASS
                bet = int((swap & ok1 & ~ok0).sum()); wor = int((swap & ~ok1 & ok0).sum())
                # events whose baseline pick was NOT the top-2 leader are
                # untouched; net is over all events for scale
                nev = len(truth)
                net = 100.0 * (bet - wor) / nev
                line += f" {net:+11.2f} {bet:6d}:{wor:<5d}"
            print(line)


if __name__ == "__main__":
    main()
