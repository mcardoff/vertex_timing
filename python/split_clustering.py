#!/usr/bin/env python3
"""
split_clustering.py -- cluster the in-jet and out-of-jet tracks SEPARATELY,
rather than clustering everything together and subsetting afterwards.

Everything tested so far builds one cluster collection from all tracks and then
picks a time from a subset. This instead runs the clustering twice on disjoint
track sets, giving two independent collections and therefore two independent
candidate times whose AGREEMENT is new information.

The clustering is a faithful reimplementation of MyUtl::doIterativeClustering
(src/clustering_functions.h), 1-D in time:
  - every track starts as its own cluster, value = time, sigma = timeRes
  - seed on the highest-pT unconsumed ORIGINAL single-track cluster
  - grow by repeatedly absorbing the single nearest unconsumed original cluster
    with d < distCut, recomputing the centroid after each absorption
  - d(a,b) = |t_a - t_b| / sqrt(sigma_a^2 + sigma_b^2)
  - merge: inverse-variance mean, sigma_new = s1*s2/sqrt(s1^2+s2^2)

VALIDATED before use: run on the full track set it must reproduce the partition
the C++ actually produced, which the export records as `cluster_idx`. The script
refuses to report subset results unless that check passes.
"""
import argparse
import numpy as np
import pandas as pd
import uproot

EVT = ["sample_id", "file_idx", "event_num"]
KEY = EVT + ["cluster_idx"]
DIST_CUT = 3.0
A_DZ = 0.7


def cluster_1d(t, s, pt, distcut=DIST_CUT):
    """Returns a label per track. Mirrors doIterativeClustering exactly."""
    n = len(t)
    lab = np.full(n, -1, np.int32)
    consumed = np.zeros(n, bool)
    k = 0
    while True:
        free = ~consumed
        if not free.any(): break
        seed = np.flatnonzero(free)[np.argmax(pt[free])]   # highest-pT unconsumed
        consumed[seed] = True
        lab[seed] = k
        cv, cs = t[seed], s[seed]
        while True:
            free = ~consumed
            if not free.any(): break
            idx = np.flatnonzero(free)
            d = np.abs(cv - t[idx]) / np.sqrt(cs * cs + s[idx] * s[idx])
            j = np.argmin(d)
            if d[j] >= distcut: break                       # strict: absorb if d < cut
            b = idx[j]
            consumed[b] = True
            lab[b] = k
            w1, w2 = 1.0 / (cs * cs), 1.0 / (s[b] * s[b])
            cv = (cv * w1 + t[b] * w2) / (w1 + w2)
            cs = (cs * s[b]) / np.sqrt(cs * cs + s[b] * s[b])
        k += 1
    return lab


def same_partition(a, b):
    """Partition equality, independent of label numbering."""
    return (pd.crosstab(a, b) > 0).sum(axis=1).max() == 1 and \
           (pd.crosstab(b, a) > 0).sum(axis=1).max() == 1


def run(df, mask=None):
    """Cluster each event (optionally a track subset) and return per-track labels."""
    out = np.full(len(df), -1, np.int32)
    sub = np.ones(len(df), bool) if mask is None else mask.to_numpy()
    tt = df["time"].to_numpy(); ss = df["timeRes"].to_numpy(); pp = df["pt"].to_numpy()
    for _, idx in df.groupby(EVT, sort=False).indices.items():
        i = idx[sub[idx]]
        if len(i) == 0: continue
        out[i] = cluster_1d(tt[i], ss[i], pp[i])
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--file", required=True)
    ap.add_argument("--tag", default="")
    ap.add_argument("--validate-events", type=int, default=4000)
    a = ap.parse_args()

    c = uproot.open(a.file)["clusters"].arrays(
        KEY + ["cluster_time", "delta_t"], library="pd")
    t = uproot.open(a.file)["tracks"].arrays(
        KEY + ["track_idx", "pt", "time", "timeRes", "z0_pull_pv", "sigma_z0",
               "dr_nearest_fwdjet", "truth_is_hs"], library="pd")
    t = t.dropna(subset=["time", "timeRes", "z0_pull_pv", "sigma_z0"])
    t = t[t["timeRes"] > 0].reset_index(drop=True)
    t["sdz"] = t["z0_pull_pv"] * t["sigma_z0"]          # SIGNED z0 - z_PV
    t["ivz"] = 1.0 / t["sigma_z0"] ** 2
    t["iv"] = 1.0 / t["timeRes"] ** 2
    t["inj"] = (t["dr_nearest_fwdjet"] < 0.4).fillna(False)
    c["t_truth"] = c["cluster_time"] - c["delta_t"]
    truth = c.groupby(EVT, sort=False)["t_truth"].first(); truth.index.names = EVT

    # ---- validation against the C++ partition -----------------------------
    keys = t[EVT].drop_duplicates().head(a.validate_events)
    v = t.merge(keys, on=EVT, how="inner")
    lab = run(v)
    ok = bad = 0
    for _, g in v.assign(mine=lab).groupby(EVT, sort=False):
        if same_partition(g["mine"], g["cluster_idx"]): ok += 1
        else: bad += 1
    print(f"\n  clustering validation: {ok}/{ok+bad} events reproduce the C++ "
          f"partition exactly ({100*ok/(ok+bad):.2f}%)")
    if ok / (ok + bad) < 0.98:
        print("  REFUSING to report subset results -- reimplementation is not faithful.")
        return

    # ---- separate clusterings ---------------------------------------------
    print(f"\n===== {a.tag or a.file}  |  clustering in-jet and out-of-jet separately")
    res = {}
    for lab_, mask in (("all", None), ("in-jet", t["inj"]), ("out-jet", ~t["inj"])):
        t[f"L_{lab_}"] = run(t, mask)
        d = t[t[f"L_{lab_}"] >= 0].copy()
        g = d.groupby(EVT + [f"L_{lab_}"], sort=False)
        cl = pd.DataFrame({
            "tt": g.apply(lambda x: (x["iv"] * x["time"]).sum() / x["iv"].sum(),
                          include_groups=False),
            "dz": g.apply(lambda x: (x["ivz"] * x["sdz"]).sum() / x["ivz"].sum(),
                          include_groups=False),
            "sumw": g.apply(lambda x: (x["pt"] * np.exp(-A_DZ * x["sdz"].abs())).sum(),
                            include_groups=False),
            "n": g.size(),
            "pur": g.apply(lambda x: (x["pt"] * (x["truth_is_hs"] > 0.5)).sum() / x["pt"].sum(),
                           include_groups=False)}).reset_index()
        cl["score"] = cl["sumw"] * np.exp(-1.5 * cl["dz"].abs())
        pick = cl.loc[cl.groupby(EVT, sort=False)["score"].idxmax()]
        pick = pick.set_index(pd.MultiIndex.from_frame(pick[EVT])); pick.index.names = EVT
        res[lab_] = pick.reindex(truth.index)
        dt = res[lab_]["tt"] - truth
        print(f"    {lab_:<8s} clusters/event {cl.groupby(EVT).size().mean():5.2f}   "
              f"core frac {100*(dt.abs()<60).mean():6.2f}%   "
              f"sel. purity {100*res[lab_]['pur'].mean():5.2f}%   "
              f"n trk {res[lab_]['n'].mean():5.1f}")

    dA = (res["all"]["tt"] - truth).to_numpy()
    dI = (res["in-jet"]["tt"] - truth).to_numpy()
    dO = (res["out-jet"]["tt"] - truth).to_numpy()
    sep = np.abs(res["in-jet"]["tt"].to_numpy() - res["out-jet"]["tt"].to_numpy())
    print(f"\n  agreement between the two independent collections:")
    for cut in (30, 60, 100, 200):
        m = sep < cut
        if m.sum() < 100: continue
        print(f"    |t_in - t_out| < {cut:3d} ps : {100*m.mean():5.1f}% of events, "
              f"and there all-track core frac is {100*(np.abs(dA[m])<60).mean():6.2f}% "
              f"(vs {100*(np.abs(dA[~m])<60).mean():6.2f}% when they disagree)")
    rng = np.random.default_rng(0)
    s1 = np.nanstd(dI - dA); s2 = np.nanstd(dO - dA)
    N = np.vstack([dA, dA + rng.normal(0, s1, len(dA)), dA + rng.normal(0, s2, len(dA))])
    nullb = N[np.argmin(np.abs(N), axis=0), np.arange(N.shape[1])]
    O = np.vstack([dA, dI, dO])
    orc = O[np.argmin(np.abs(O), axis=0), np.arange(O.shape[1])]
    print(f"\n  oracle over the three collections : {100*(np.abs(orc)<60).mean():6.2f}%")
    print(f"  matched-noise null                : {100*(np.abs(nullb)<60).mean():6.2f}%")
    print(f"  REAL complementarity              : "
          f"{100*(np.abs(orc)<60).mean()-100*(np.abs(nullb)<60).mean():+6.2f}")


if __name__ == "__main__":
    main()
