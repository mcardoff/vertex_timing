#!/usr/bin/env python3
"""
injet_decomposition.py -- separate what WAVeS's SELECTION is worth from what its
in-jet RE-TIMING is worth, and test both against the aggregate t0.

The deployed `Score::WAVES` row does two independent things:
  1. ranks clusters by the WAVeS score, and
  2. recomputes the chosen cluster's time from only its tracks within
     dR < 0.4 of a qualifying forward jet (Cluster::calculateTime).
Every published WAVeS number folds the two together. This runs the full
selection x timing grid so they can be read apart, and adds TRKPTZ and
TRKPTZ_TZ under the same in-jet timing -- i.e. does the re-timing transfer to a
selector that never looked at a jet?

The in-jet rule is reproduced exactly: `dr_nearest_fwdjet` in the track tree
comes from `nearestForwardJet`, which applies the SAME jet qualification as
calculateTime's WAVES branch (isJetRemoved, pT > MIN_JET_PT, forward eta band),
so `dr_nearest_fwdjet < 0.4` is the same predicate. The fallback matches too:
when no track survives the filter, calculateTime returns values[0], the
full-cluster time.
"""
import argparse
import numpy as np
import pandas as pd
import uproot

EVT = ["sample_id", "file_idx", "event_num"]
KEY = EVT + ["cluster_idx"]
A_DZ = 0.7
DR_JET = 0.4


def wmedian(df, wcol, vcol):
    d = df.sort_values(EVT + [vcol])
    codes, uniq = pd.factorize(pd.MultiIndex.from_frame(d[EVT]), sort=False)
    if not isinstance(uniq, pd.MultiIndex):
        uniq = pd.MultiIndex.from_tuples([tuple(x) for x in uniq], names=EVT)
    w = d[wcol].to_numpy(float)
    tot = np.bincount(codes, weights=w)
    start = np.r_[0, np.flatnonzero(np.diff(codes)) + 1]
    cw = np.cumsum(w)
    off = np.zeros(len(d)); off[start[1:]] = cw[start[1:] - 1]
    cw = cw - np.maximum.accumulate(off)
    hit = cw >= 0.5 * tot[codes]
    idx = np.arange(len(d))
    firsts = np.minimum.reduceat(np.where(hit, idx, len(d) + 1), start)
    s = pd.Series(d[vcol].to_numpy()[np.clip(firsts, 0, len(d) - 1)], index=uniq)
    s.index.names = EVT
    return s


def winsor(t, centre, win, wcol):
    j = t.merge(centre.rename("c0").rename_axis(EVT).reset_index(), on=EVT, how="left")
    j = j[(j["time"] - j["c0"]).abs() < win]
    g = j.groupby(EVT, sort=False)
    r = g.apply(lambda x: (x[wcol] * x["time"]).sum() / x[wcol].sum(), include_groups=False)
    r.index.names = EVT
    return r.reindex(centre.index).fillna(centre)


def agg_t0(t, sub=None):
    """Weighted median + double truncation over `t`, optionally restricted to a
    subset of tracks (in-jet). Events the subset empties fall back to NaN and are
    filled by the caller."""
    x = t if sub is None else t[sub]
    if len(x) == 0: return None
    return winsor(x, winsor(x, wmedian(x, "w", "time"), 100.0, "wiv"), 60.0, "wiv")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--file", required=True)
    ap.add_argument("--tag", default="")
    a = ap.parse_args()

    c = uproot.open(a.file)["clusters"].arrays(
        KEY + ["cluster_time", "delta_z", "delta_t", "trkptz_score", "waves_score",
               "n_forward_jets", "truth_n_hs_tracks"], library="pd")
    t = uproot.open(a.file)["tracks"].arrays(
        KEY + ["pt", "time", "timeRes", "z0_pull_pv", "sigma_z0",
               "dr_nearest_fwdjet"], library="pd")
    t = t.dropna(subset=["z0_pull_pv", "sigma_z0", "time", "timeRes"])
    t = t[t["timeRes"] > 0]
    t["dz"] = t["z0_pull_pv"].abs() * t["sigma_z0"]
    t["w"] = t["pt"] * np.exp(-A_DZ * t["dz"])
    t["wiv"] = t["w"] / t["timeRes"] ** 2
    t["iv"] = 1.0 / t["timeRes"] ** 2
    t["injet"] = (t["dr_nearest_fwdjet"] < DR_JET).fillna(False)

    c["t_truth"] = c["cluster_time"] - c["delta_t"]
    c["dzcl"] = np.exp(-1.5 * c["delta_z"].abs())
    truth = c.groupby(EVT, sort=False)["t_truth"].first(); truth.index.names = EVT
    nev = len(truth)
    print(f"\n===== {a.tag or a.file}  |  {nev:,} events, {len(c):,} clusters, "
          f"{len(t):,} timed tracks  ({100*t['injet'].mean():.1f}% in-jet)")

    # per-cluster in-jet time (inverse-variance over in-jet tracks), with the
    # same fallback to the full-cluster time the C++ uses
    ij = t[t["injet"]]
    g = ij.groupby(KEY, sort=False)
    tij = (g.apply(lambda x: (x["iv"] * x["time"]).sum() / x["iv"].sum(),
                   include_groups=False).rename("t_injet"))
    c = c.merge(tij.reset_index(), on=KEY, how="left")
    c["t_injet"] = c["t_injet"].fillna(c["cluster_time"])

    # TRKPTZ_TZ score
    t["_s"] = t["w"]
    s = t.groupby(KEY, sort=False)["_s"].sum().rename("tzsum").reset_index()
    c = c.merge(s, on=KEY, how="left")
    c["tz_score"] = c["tzsum"].fillna(0.0) * c["dzcl"]

    def pick(col):
        ip = c.groupby(EVT, sort=False)[col].idxmax()
        p = c.loc[ip].set_index(pd.MultiIndex.from_frame(c.loc[ip][EVT]))
        p.index.names = EVT
        return p.reindex(truth.index)

    SEL = [("TRKPTZ", "trkptz_score"), ("TRKPTZ_TZ", "tz_score"), ("WAVeS", "waves_score")]
    print("\n  selection x timing grid  (rows = which cluster is chosen, "
          "cols = what time it reports)")
    print(f"    {'':<12s} {'standard':>10s} {'in-jet':>10s} {'re-timing':>11s}")
    picks = {}
    for lab, col in SEL:
        p = pick(col); picks[lab] = p
        a_ = 100 * ((p["cluster_time"] - truth).abs() < 60).mean()
        b_ = 100 * ((p["t_injet"] - truth).abs() < 60).mean()
        print(f"    {lab:<12s} {a_:9.2f}% {b_:9.2f}% {b_ - a_:+10.2f}")

    # ---- aggregate t0, all tracks and in-jet only -------------------------
    t0_all = agg_t0(t)
    t0_ij = agg_t0(t, t["injet"])
    t0_ij = t0_ij.reindex(truth.index).fillna(t0_all.reindex(truth.index))
    okA = ((t0_all.reindex(truth.index) - truth).abs() < 60).fillna(False)
    okI = ((t0_ij - truth).abs() < 60).fillna(False)
    print(f"\n  aggregate t0, ALL timed tracks           {100*okA.mean():6.2f}%")
    print(f"  aggregate t0, IN-JET tracks only         {100*okI.mean():6.2f}%"
          f"   {100*okI.mean()-100*okA.mean():+.2f}")

    # ---- differential performance of the aggregate ------------------------
    E = pd.DataFrame(index=truth.index)
    E["agg"] = okA.to_numpy()
    E["agg_ij"] = okI.to_numpy()
    for lab in ("TRKPTZ", "TRKPTZ_TZ", "WAVeS"):
        E[lab] = ((picks[lab]["cluster_time"] - truth).abs() < 60).to_numpy()
    E["nhs"] = c.groupby(EVT, sort=False)["truth_n_hs_tracks"].first().reindex(E.index).to_numpy()
    E["njet"] = c.groupby(EVT, sort=False)["n_forward_jets"].first().reindex(E.index).to_numpy()

    for var, title, edges in [
            ("nhs", "n truth forward HS tracks", [0, 1, 2, 3, 4, 5, 6, 8, 10, 14, 20, 999]),
            ("njet", "n forward jets", [0, 1, 2, 3, 4, 999])]:
        print(f"\n  aggregate t0 differential vs {title}")
        print(f"    {'bin':>10s} {'events':>8s} {'TRKPTZ':>8s} {'TZ':>8s} "
              f"{'WAVeS':>8s} {'agg t0':>8s} {'agg-TZ':>8s}")
        v = E[var].to_numpy()
        for lo, hi in zip(edges[:-1], edges[1:]):
            m = (v >= lo) & (v < hi)
            if m.sum() < 150: continue
            r = {k: 100 * E.loc[m, k].mean() for k in
                 ("TRKPTZ", "TRKPTZ_TZ", "WAVeS", "agg")}
            lab = f"{lo}" if hi == lo + 1 else (f">={lo}" if hi == 999 else f"{lo}-{hi-1}")
            print(f"    {lab:>10s} {m.sum():8d} {r['TRKPTZ']:7.2f}% {r['TRKPTZ_TZ']:7.2f}% "
                  f"{r['WAVeS']:7.2f}% {r['agg']:7.2f}% {r['agg']-r['TRKPTZ_TZ']:+8.2f}")


if __name__ == "__main__":
    main()
