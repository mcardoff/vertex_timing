#!/usr/bin/env python3
"""
jet_term_consistency.py -- is the bounded jet term broadly useful, or does it only
help in the rare subpopulation it was derived from?

The alpha term was motivated by an AUC ranking computed on the 543 events
TRKPTZ_TZ loses to WAVeS. That is a SELECTED subset, and a variable that
separates there need not separate anywhere else. If the term flips picks rarely,
or flips them roughly as often for the worse as for the better, it is a
coincidence of that subset and is not worth carrying.

Three population-wide checks, none restricted to the losses:

  1. Population AUC. Over EVERY cluster in EVERY event, how well does
     frac_pt_in_fwdjet identify the HS-dominant cluster, against the variables
     already in the score? If it only separates on the 543, it will read ~0.5
     here.
  2. Switch accounting. Between alpha = 0 and alpha > 0, how many events change
     their pick at all, and of those, how many improve vs degrade? A term worth
     having must win decisively more often than it loses -- not merely net
     positive by a hair.
  3. Differential. Is the improvement spread across multiplicity, or concentrated
     in one corner?
"""
import argparse, gc
import numpy as np
import pandas as pd
import uproot

EVT = ["sample_id", "file_idx", "event_num"]
KEY = EVT + ["cluster_idx"]
A_DZ, GUARD, PASS = 0.7, 3, 60.0


def auc(x, y):
    m = np.isfinite(x); x, y = x[m], y[m].astype(bool)
    if y.sum() == 0 or (~y).sum() == 0: return np.nan
    r = pd.Series(x).rank().to_numpy()
    return (r[y].sum() - y.sum() * (y.sum() + 1) / 2.0) / (y.sum() * (~y).sum())


def load(path):
    c = uproot.open(path)["clusters"].arrays(
        KEY + ["cluster_time", "delta_z", "delta_t", "frac_pt_in_fwdjet",
               "n_tracks", "sumpt", "truth_n_hs_tracks", "truth_purity"], library="pd")
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
    c["dzcl"] = np.exp(-1.5 * c["delta_z"].abs())
    truth = c.groupby(EVT, sort=False)["t_truth"].first(); truth.index.names = EVT
    mi = pd.MultiIndex.from_frame(c[KEY])
    g = t[t["inj"]].groupby(KEY, sort=False)
    tij = g.apply(lambda x: (x["iv"] * x["time"]).sum() / x["iv"].sum(),
                  include_groups=False).reindex(mi).to_numpy()
    c["t_inj"] = np.where(np.isnan(tij), c["cluster_time"], tij)
    c["n_inj"] = g.size().reindex(mi).fillna(0).to_numpy()
    t["_s"] = t["w"]
    s = t.groupby(KEY, sort=False)["_s"].sum().rename("s").reset_index()
    c = c.merge(s, on=KEY, how="left")
    c["base"] = c["s"].fillna(0.0) * c["dzcl"]
    c["t_out"] = np.where(c["n_inj"] >= GUARD, c["t_inj"], c["cluster_time"])
    c["fjet"] = c["frac_pt_in_fwdjet"].fillna(0.0).clip(0, 1)
    del t; gc.collect()
    return c, truth


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dir", default="/Users/mcard/project/vertex_timing/training_new")
    ap.add_argument("--samples", default="vbf,zjets")
    ap.add_argument("--sel", default="novbs")
    a = ap.parse_args()

    for smp in a.samples.split(","):
        print(f"\n{'='*72}\n===== {smp}", flush=True)
        c, truth = load(f"{a.dir}/{smp}_{a.sel}_training.root")
        y = truth.to_numpy()

        # --- 1. population AUC over EVERY cluster ------------------------
        # label: is this the HS-dominant cluster of its event?
        idx = c.groupby(EVT, sort=False)["truth_n_hs_tracks"].idxmax()
        lab = np.zeros(len(c), bool); lab[c.index.get_indexer(idx)] = True
        multi = c.groupby(EVT, sort=False)["cluster_idx"].transform("size") > 1
        print(f"\n  1. POPULATION AUC -- identify the HS-dominant cluster among all "
              f"{int(multi.sum()):,} clusters in multi-cluster events")
        for q in ("fjet", "base", "sumpt", "n_tracks"):
            v = c.loc[multi, q].to_numpy(float)
            print(f"     {q:<12s} {auc(v, lab[multi.to_numpy()]):.3f}")
        v = np.abs(c.loc[multi, "delta_z"].to_numpy(float))
        print(f"     {'|delta_z|':<12s} {1-auc(v, lab[multi.to_numpy()]):.3f}   (inverted: "
              f"smaller |dz| marks the HS cluster)")

        # --- 2. switch accounting ---------------------------------------
        def pickvec(alpha):
            c["_sc"] = c["base"] * (1.0 + alpha * c["fjet"])
            ip = c.groupby(EVT, sort=False)["_sc"].idxmax()
            p = c.loc[ip].set_index(pd.MultiIndex.from_frame(c.loc[ip][EVT]))
            p.index.names = EVT; p = p.reindex(truth.index)
            return p["cluster_idx"].to_numpy(), np.abs(p["t_out"].to_numpy() - y) < PASS
        ci0, ok0 = pickvec(0.0)
        print(f"\n  2. SWITCH ACCOUNTING (baseline alpha=0: {100*ok0.mean():.2f}%)")
        print(f"     {'alpha':>6s} {'switched':>9s} {'better':>8s} {'worse':>7s} "
              f"{'neutral':>8s} {'ratio':>7s} {'net':>7s}")
        for al in (0.25, 0.5, 1.0, 2.0, 4.0):
            ci, ok = pickvec(al)
            sw = ci != ci0
            better = int((sw & ok & ~ok0).sum()); worse = int((sw & ~ok & ok0).sum())
            neut = int(sw.sum()) - better - worse
            r = better / worse if worse else np.inf
            print(f"     {al:6.2f} {100*sw.mean():8.2f}% {better:8d} {worse:7d} "
                  f"{neut:8d} {r:7.2f} {100*(ok.mean()-ok0.mean()):+6.2f}")

        # --- 3. differential ---------------------------------------------
        ci1, ok1 = pickvec(1.0)
        nhs = c.groupby(EVT, sort=False)["truth_n_hs_tracks"].sum().reindex(truth.index).to_numpy()
        print(f"\n  3. DIFFERENTIAL of alpha=1 vs alpha=0, by n forward HS tracks")
        print(f"     {'bin':>7s} {'events':>8s} {'a=0':>8s} {'a=1':>8s} {'delta':>7s} "
              f"{'switched':>9s}")
        for lo, hi in [(0,2),(2,4),(4,6),(6,9),(9,13),(13,20),(20,10**9)]:
            m = (nhs >= lo) & (nhs < hi)
            if m.sum() < 300: continue
            lab_ = f"{lo}-{hi-1}" if hi < 10**8 else f">={lo}"
            print(f"     {lab_:>7s} {int(m.sum()):8d} {100*ok0[m].mean():7.2f}% "
                  f"{100*ok1[m].mean():7.2f}% {100*(ok1[m].mean()-ok0[m].mean()):+6.2f} "
                  f"{100*(ci1[m]!=ci0[m]).mean():8.2f}%")
        del c; gc.collect()


if __name__ == "__main__":
    main()
