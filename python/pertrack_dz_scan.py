#!/usr/bin/env python3
"""
pertrack_dz_scan.py -- does weighting TRKPTZ's pT sum by each track's OWN
distance to the primary vertex beat the incumbent cluster-level exp(-1.5|dz|)?

Incumbent:   score = (SUM_t pT_t) * exp(-1.5 |z_cluster - z_PV|)
Per-track:   score = SUM_t pT_t * f(x_t) [* exp(-1.5 |dz_cluster|) if `keep_cl`]

with x_t either the raw distance |z0_t - z_PV| in mm or the significance
|z0_t - z_PV| / sigma_z0 (the exported `z0_pull_pv`).

Prior evidence says this should LOSE: the WAVES comment in
src/clustering_structs.h states the cluster-level z-term beats per-track z-pull
weighting "because it averages over track-by-track z noise". This scans it
rather than taking that on trust, and on Z+jets rather than VBF.

Metric is core fraction: fraction of events whose argmax cluster has
|delta_t| < 60 ps. Every variant is evaluated on the SAME events and the SAME
clusters, so the comparison is paired and only the ranking changes.
"""
import argparse
import numpy as np
import pandas as pd
import uproot

EVT = ["sample_id", "file_idx", "event_num"]
KEY = EVT + ["cluster_idx"]


def load(path):
    c = uproot.open(path)["clusters"].arrays(
        KEY + ["sumpt", "delta_z", "delta_t", "trkptz_score"], library="pd")
    t = uproot.open(path)["tracks"].arrays(
        KEY + ["pt", "z0_pull_pv", "sigma_z0"], library="pd")
    t = t.dropna(subset=["z0_pull_pv", "sigma_z0"])
    t["sig"] = t["z0_pull_pv"].abs()
    t["dz"] = t["sig"] * t["sigma_z0"]          # back out the raw |z0 - z_PV|
    c["ok"] = c["delta_t"].abs() < 60
    c["dzcl"] = np.exp(-1.5 * c["delta_z"].abs())
    return c, t


def score_cf(c, t, wcol, keep_cl):
    """Aggregate a per-track weight into a cluster score, take the argmax per
    event, and return the core fraction."""
    t["_w"] = t["pt"] * wcol
    s = t.groupby(KEY, sort=False)["_w"].sum().rename("s").reset_index()
    m = c.merge(s, on=KEY, how="left")
    m["s"] = m["s"].fillna(0.0)
    if keep_cl:
        m["s"] *= m["dzcl"]
    pick = m.loc[m.groupby(EVT, sort=False)["s"].idxmax()]
    return 100.0 * pick["ok"].mean()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--file", required=True)
    ap.add_argument("--tag", default="")
    a = ap.parse_args()
    c, t = load(a.file)
    nev = c[EVT].drop_duplicates().shape[0]
    print(f"\n===== {a.tag or a.file}  |  {nev:,} events, {len(c):,} clusters, "
          f"{len(t):,} timed tracks")

    # incumbent, recomputed the same way every variant is, so the comparison is
    # like-for-like rather than against the stored column
    base = score_cf(c, t, 1.0, keep_cl=True)
    print(f"\n  TRKPTZ incumbent  (cluster-level exp(-1.5|dz_cl|))   {base:6.2f}%")
    flat = score_cf(c, t, 1.0, keep_cl=False)
    print(f"  no z term at all  (plain sum pT)                      {flat:6.2f}%"
          f"   {flat - base:+.2f}")

    # (label, per-track weight) -- built once, evaluated under both keep_cl modes
    variants = (
        [(f"exp(-{a_}*dz)",     np.exp(-a_ * t["dz"]))          for a_ in (0.5, 1.5, 3.0, 6.0)]
      + [(f"1/max(dz,{f_})",    1.0 / np.maximum(t["dz"], f_))  for f_ in (0.01, 0.05, 0.10)]
      + [(f"exp(-{b_}*sig)",    np.exp(-b_ * t["sig"]))         for b_ in (0.2, 0.5, 1.0)]
      + [(f"1/max(sig,{f_})",   1.0 / np.maximum(t["sig"], f_)) for f_ in (0.5, 1.0, 2.0)]
    )
    out = [(name, k, score_cf(c, t, w, k))
           for k in (False, True) for name, w in variants]

    for keep_cl in (False, True):
        lbl = ("replaces the cluster term" if not keep_cl
               else "multiplies the cluster term")
        print(f"\n  per-track weight {lbl}:")
        for name, k, cf in out:
            if k != keep_cl: continue
            print(f"    {name:<20s} {cf:6.2f}%   {cf - base:+.2f}")


if __name__ == "__main__":
    main()
