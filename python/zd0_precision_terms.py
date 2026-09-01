#!/usr/bin/env python3
"""
zd0_precision_terms.py -- the qP-precision construction applied to z0 and d0.

score' = base * (SUM_t 1/var_x)^{beta/2}   for x in {z0, d0, qOverP}

cluster_z_sigma / cluster_d0_sigma / cluster_qOverP_sigma are all
1/sqrt(SUM 1/var) over the same constituents, so the three are siblings and
likely highly correlated -- the scan therefore reports:

  0. mutual correlations of the three log-precision sums;
  1. each term ALONE on the TZ + guarded-timing base, beta scanned with ONE
     global value on all four samples;
  2. each of z0/d0 STACKED on the shipped qP term (does it add anything the
     qP factor does not already carry?).

Verdict rule unchanged: a term must be >= -0.05 on every sample to be safe,
and is only worth carrying if its worst-sample gain beats what is already
there.
"""
import argparse, gc
import numpy as np
import pandas as pd
import uproot

EVT = ["sample_id", "file_idx", "event_num"]
KEY = EVT + ["cluster_idx"]
A_DZ, GUARD, PASS, BQP = 0.7, 3, 60.0, 0.8
SIG = {"z0": "cluster_z_sigma", "d0": "cluster_d0_sigma", "qP": "cluster_qOverP_sigma"}


def load(path):
    c = uproot.open(path)["clusters"].arrays(
        KEY + ["cluster_time", "delta_z", "delta_t"] + list(SIG.values()), library="pd")
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
    c["t_out"] = np.where(c["n_inj"] >= GUARD, c["t_inj"], c["cluster_time"])
    for k, col in SIG.items():
        c[f"p_{k}"] = 1.0 / np.maximum(c[col].to_numpy(float), 1e-12)  # 1/sigma
    del t; gc.collect()
    return c, truth


def cf(c, truth, sc):
    c["_sc"] = sc
    ip = c.groupby(EVT, sort=False)["_sc"].idxmax()
    p = c.loc[ip].set_index(pd.MultiIndex.from_frame(c.loc[ip][EVT]))
    p.index.names = EVT
    return 100.0 * (np.abs(p["t_out"].reindex(truth.index).to_numpy()
                           - truth.to_numpy()) < PASS).mean()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dir", default="/Users/mcard/project/vertex_timing/training_new")
    ap.add_argument("--samples", default="vbf,zjets,dijet,ttbar")
    a = ap.parse_args()
    smps = a.samples.split(",")
    BETAS = (0.0, 0.2, 0.4, 0.8, 1.6)
    R = {}     # (mode, key, beta) -> {sample: cf}
    corr = {}
    for smp in smps:
        print(f"  loading {smp} ...", flush=True)
        c, truth = load(f"{a.dir}/{smp}_novbs_training.root")
        # 0. correlations on multi-cluster events
        multi = c.groupby(EVT, sort=False)["cluster_idx"].transform("size") > 1
        L = {k: np.log(c.loc[multi, f"p_{k}"].to_numpy(float)) for k in SIG}
        corr[smp] = {(x, y): float(np.corrcoef(L[x], L[y])[0, 1])
                     for x, y in (("z0", "d0"), ("z0", "qP"), ("d0", "qP"))}
        base = c["tz"].to_numpy(float)
        for k in SIG:
            for b in BETAS:
                R.setdefault(("alone", k, b), {})[smp] = cf(
                    c, truth, base * np.power(c[f"p_{k}"].to_numpy(float), b))
        qp = base * np.power(c["p_qP"].to_numpy(float), BQP)
        for k in ("z0", "d0"):
            for b in BETAS:
                R.setdefault(("stack", k, b), {})[smp] = cf(
                    c, truth, qp * np.power(c[f"p_{k}"].to_numpy(float), b))
        del c; gc.collect()

    print("\n  0. mutual correlations of log(1/sigma) (multi-cluster events)")
    for smp in smps:
        print(f"     {smp:6s} " + "  ".join(f"{x}-{y} {corr[smp][(x,y)]:+.3f}"
                                            for x, y in (("z0","d0"),("z0","qP"),("d0","qP"))))

    for mode, title, keys in (("alone", "1. ALONE on TZ + guarded timing "
                               "(beta=0 is TZ_GIJ)", list(SIG)),
                              ("stack", "2. STACKED on the shipped qP term "
                               "(beta=0 is TZQ)", ["z0", "d0"])):
        print(f"\n  {title}")
        for k in keys:
            print(f"     (1/sigma_{k})^beta")
            print(f"     {'beta':>6s} " + " ".join(f"{s:>9s}" for s in smps) + "   worst d")
            for b in BETAS:
                v = R[(mode, k, b)]
                d = [v[s] - R[(mode, k, 0.0)][s] for s in smps]
                print(f"     {b:6.2f} " + " ".join(f"{v[s]:8.2f}%" for s in smps)
                      + f"   {min(d):+6.2f}"
                      + ("   safe on all" if min(d) >= -0.05 and b > 0 else ""))


if __name__ == "__main__":
    main()
