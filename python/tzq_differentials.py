#!/usr/bin/env python3
"""
tzq_differentials.py -- the five-method ladder differentially on zjets, dijet
and ttbar (novbs exports; the local C++ pipeline only reaches VBF ntuples).

Methods, matching the C++ registry rows exactly:
  TRKPTZ    argmax trkptz_score, full-cluster time
  TZ        argmax SUM pT e^(-0.7|dz_t|) x e^(-1.5|dz_cl|), full-cluster time
  TZ_GIJ    TZ selection, in-jet time when the cluster holds >= 3 in-jet tracks
  TZQ       TZ x (1/sigma_qP)^0.8 selection, guarded in-jet time
  WAVES     argmax waves_score, in-jet time ALWAYS (deployed behaviour)

Differential axes: event-level n truth forward HS tracks (sum of the
per-cluster column over the event -- the corrected aggregation), and an event
PU-fraction PROXY, 1 - nHS/n_fwd_tracks_reco (the exports do not carry the
C++ puRatio directly; labelled as a proxy for that reason).
"""
import argparse, gc
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import uproot

EVT = ["sample_id", "file_idx", "event_num"]
KEY = EVT + ["cluster_idx"]
A_DZ, GUARD, PASS, BETA = 0.7, 3, 60.0, 0.8


def load(path):
    c = uproot.open(path)["clusters"].arrays(
        KEY + ["cluster_time", "delta_z", "delta_t", "trkptz_score", "waves_score",
               "cluster_qOverP_sigma", "truth_n_hs_tracks", "n_fwd_tracks_reco"],
        library="pd")
    t = uproot.open(path)["tracks"].arrays(
        KEY + ["pt", "time", "timeRes", "z0_pull_pv", "sigma_z0",
               "dr_nearest_fwdjet"], library="pd")
    t = t.dropna(subset=["time", "timeRes", "z0_pull_pv", "sigma_z0"])
    t = t[t["timeRes"] > 0]
    t["dz"] = t["z0_pull_pv"].abs() * t["sigma_z0"]
    t["w"] = t["pt"] * np.exp(-A_DZ * t["dz"])
    t["iv"] = 1.0 / t["timeRes"] ** 2
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
    c["tz"] = c["s"].fillna(0.0) * np.exp(-1.5 * c["delta_z"].abs())
    c["tzq"] = c["tz"] * np.power(
        np.maximum(c["cluster_qOverP_sigma"].to_numpy(float), 1e-9), -BETA)
    c["t_out"] = np.where(c["n_inj"] >= GUARD, c["t_inj"], c["cluster_time"])
    del t; gc.collect()
    return c, truth


def okvec(c, truth, scol, tcol):
    ip = c.groupby(EVT, sort=False)[scol].idxmax()
    p = c.loc[ip].set_index(pd.MultiIndex.from_frame(c.loc[ip][EVT]))
    p.index.names = EVT
    return (np.abs(p[tcol].reindex(truth.index).to_numpy()
                   - truth.to_numpy()) < PASS)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dir", default="/Users/mcard/project/vertex_timing/training_new")
    ap.add_argument("--samples", default="zjets,dijet,ttbar")
    ap.add_argument("--out", default="figs/tzq_differentials.pdf")
    a = ap.parse_args()
    smps = a.samples.split(",")
    METH = [("TRKPTZ", "trkptz_score", "cluster_time"),
            ("TZ", "tz", "cluster_time"),
            ("TZ_GIJ", "tz", "t_out"),
            ("TZQ", "tzq", "t_out"),
            ("WAVES", "waves_score", "t_inj")]
    NHSB = [(0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (5, 6), (6, 8), (8, 10),
            (10, 14), (14, 20), (20, 10**9)]
    PUB = [(0.0, .5), (.5, .6), (.6, .7), (.7, .8), (.8, .85), (.85, .9),
           (.9, .95), (.95, 1.001)]

    fig, axes = plt.subplots(2, len(smps), figsize=(5.2 * len(smps), 8.5),
                             sharey="row")
    col = {"TRKPTZ": "#d62728", "TZ": "#2ca02c", "TZ_GIJ": "#1f77b4",
           "TZQ": "#9467bd", "WAVES": "#bcbd22"}
    for j, smp in enumerate(smps):
        print(f"\n===== {smp} (novbs)", flush=True)
        c, truth = load(f"{a.dir}/{smp}_novbs_training.root")
        OK = {m: okvec(c, truth, scol, tcol) for m, scol, tcol in METH}
        gev = c.groupby(EVT, sort=False)
        nhs = gev["truth_n_hs_tracks"].sum().reindex(truth.index).to_numpy()
        nrec = gev["n_fwd_tracks_reco"].first().reindex(truth.index).to_numpy()
        puf = np.clip(1.0 - nhs / np.maximum(nrec, 1), 0, 1)
        print("  inclusive: " + "  ".join(f"{m} {100*OK[m].mean():.2f}%" for m, _, _ in METH))

        for row, (v, bins, lab) in enumerate(
                ((nhs, NHSB, "n truth fwd HS tracks (event)"),
                 (puf, PUB, "event PU fraction proxy [1 - nHS/nReco]"))):
            print(f"\n  vs {lab}")
            print(f"    {'bin':>10s} {'events':>8s} " +
                  " ".join(f"{m:>8s}" for m, _, _ in METH) + f" {'TZQ-WAV':>8s}")
            xs, ys = [], {m: [] for m, _, _ in METH}
            for lo, hi in bins:
                msk = (v >= lo) & (v < hi)
                if msk.sum() < 200: continue
                blab = (f"{lo:g}" if hi == lo + 1 and row == 0 else
                        (f">={lo:g}" if hi > 10**8 else f"{lo:g}-{hi:g}"))
                r = {m: 100 * OK[m][msk].mean() for m, _, _ in METH}
                print(f"    {blab:>10s} {int(msk.sum()):8d} " +
                      " ".join(f"{r[m]:7.2f}%" for m, _, _ in METH) +
                      f" {r['TZQ']-r['WAVES']:+8.2f}")
                xs.append(0.5 * (lo + min(hi, (v.max() if hi > 10**8 else hi))))
                for m, _, _ in METH: ys[m].append(r[m])
            ax = axes[row, j]
            for m, _, _ in METH:
                ax.plot(xs, ys[m], "o-", ms=3, lw=1.4, color=col[m], label=m)
            ax.set_xlabel(lab, fontsize=9)
            if row == 0: ax.set_title(f"{smp}  ({len(truth):,} events)", fontsize=11)
            if j == 0: ax.set_ylabel("core fraction [%]", fontsize=10)
            ax.grid(alpha=0.3)
        del c; gc.collect()
    axes[0, 0].legend(fontsize=8)
    fig.suptitle("Five-method ladder, differential (novbs exports)", fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    fig.savefig(a.out)
    print(f"\n  wrote {a.out}")


if __name__ == "__main__":
    main()
