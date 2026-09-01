#!/usr/bin/env python3
"""
missing_information.py -- what reco information separates the cluster we WRONGLY
pick from the true hard-scatter cluster, on the events we lose?

Restricted to the clean losses (a right answer existed, the two methods picked
different clusters, ours is genuinely pileup-like). Within those events there are
exactly two clusters of interest: the one we picked (label 0) and the
HS-dominant one (label 1). Every reco cluster column is ranked by how well it
separates them, which is a direct answer to "what are we not using".

Columns already in our score (sumpt, delta_z) are reported too, as the control:
they should sit at or below 0.5, since by construction our score preferred the
wrong one.

Truth columns are excluded from the ranking but `truth_purity` is kept visible as
a sanity ceiling -- it is the label in disguise and must come out on top.
"""
import argparse
import numpy as np
import pandas as pd
import uproot

EVT = ["sample_id", "file_idx", "event_num"]
KEY = EVT + ["cluster_idx"]
A_DZ, GUARD, PASS = 0.7, 3, 60.0
OURS = {"sumpt", "delta_z", "trkptz_score", "tz"}


def auc(x, y):
    m = np.isfinite(x)
    x, y = x[m], y[m].astype(bool)
    if y.sum() == 0 or (~y).sum() == 0: return np.nan
    r = pd.Series(x).rank().to_numpy()
    return (r[y].sum() - y.sum() * (y.sum() + 1) / 2.0) / (y.sum() * (~y).sum())


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--file", default="figs/hists/vbf_mjj500p0_training.root")
    a = ap.parse_args()

    c = uproot.open(a.file)["clusters"].arrays(library="pd")
    t = uproot.open(a.file)["tracks"].arrays(
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
    c["tz"] = c["s"].fillna(0.0) * c["dzcl"]
    c["t_ours"] = np.where(c["n_inj"] >= GUARD, c["t_inj"], c["cluster_time"])

    def pick(col):
        ip = c.groupby(EVT, sort=False)[col].idxmax()
        p = c.loc[ip].set_index(pd.MultiIndex.from_frame(c.loc[ip][EVT]))
        p.index.names = EVT
        return p.reindex(truth.index)

    po, pw = pick("tz"), pick("waves_score")
    hb = c.sort_values(["truth_n_hs_tracks", "n_tracks"]).groupby(EVT, sort=False).tail(1)
    hb = hb.set_index(pd.MultiIndex.from_frame(hb[EVT])); hb.index.names = EVT
    hb = hb.reindex(truth.index)
    y = truth.to_numpy()
    ok_o = np.abs(po["t_ours"].to_numpy() - y) < PASS
    ok_w = np.abs(pw["t_inj"].to_numpy() - y) < PASS
    lose = (ok_w & ~ok_o
            & (np.abs(hb["t_ours"].to_numpy() - y) < PASS)
            & (po["cluster_idx"].to_numpy() != pw["cluster_idx"].to_numpy())
            & (po["truth_n_hs_tracks"].to_numpy() <= 2)
            & (po["truth_purity"].to_numpy() <= 0.50)
            & (po["n_tracks"].to_numpy() >= 3))
    print(f"  {int(lose.sum())} clean losses\n")

    bad, good = po[lose], hb[lose]
    num = [q for q in c.columns
           if c[q].dtype != object and q not in KEY + ["t_truth", "dzcl", "s"]
           and not q.startswith("truth_")]
    rows = []
    for q in num:
        x = np.r_[bad[q].to_numpy(float), good[q].to_numpy(float)]
        lab = np.r_[np.zeros(len(bad)), np.ones(len(good))]
        v = auc(x, lab)
        if not np.isfinite(v): continue
        rows.append((q, v, abs(v - 0.5)))
    R = pd.DataFrame(rows, columns=["var", "auc", "sep"]).sort_values("sep", ascending=False)
    R["used"] = R["var"].isin(OURS)
    print("  separating the cluster WE PICKED (label 0) from the true HS cluster (label 1)")
    print("  AUC > 0.5 => higher value marks the TRUE cluster; < 0.5 => marks OUR wrong one\n")
    print(f"  {'variable':<30s} {'AUC':>6s}  {'note':<28s}")
    for _, r in R.head(22).iterrows():
        note = "ALREADY IN OUR SCORE" if r["used"] else ""
        print(f"  {r['var']:<30s} {r['auc']:6.3f}  {note}")
    print(f"\n  our score's own inputs, for contrast:")
    for q in ("sumpt", "delta_z", "trkptz_score", "tz", "n_tracks"):
        if q in R["var"].values:
            print(f"    {q:<28s} {float(R[R['var']==q].auc.iloc[0]):6.3f}")
    print(f"\n  truth_purity (the label in disguise, as a ceiling): "
          f"{auc(np.r_[bad['truth_purity'].to_numpy(float), good['truth_purity'].to_numpy(float)], np.r_[np.zeros(len(bad)), np.ones(len(good))]):.3f}")


if __name__ == "__main__":
    main()
