#!/usr/bin/env python3
"""
classical_t0.py -- a fully classical (no model, no truth) event-level t0, and
the cluster score that measures agreement with it.

Motivation: in the learned selector the two "distance to the event aggregate"
features were 42.9% of total gain, the top two by a wide margin. Nothing about
that construction needs a model -- the learned per-track P(HS) can be replaced
by the weight already validated for TRKPTZ_TZ,

    w_t = pT_t * exp(-TRACK_DZ_WEIGHT * |z0_t - z_PV|)

which uses only reconstructed quantities. This measures:

  A) the classical aggregate t0 as a DIRECT answer  (|t0 - t_truth| < 60 ps),
     i.e. how much of the ML pipeline's aggregate-t0 result is the median-over-
     all-tracks structure rather than the learned weights; and
  B) TRKPTZ_TZ multiplied by an agreement term exp(-|t_clus - t0|/tau), i.e.
     using the aggregate only to RANK clusters, which keeps the reported time a
     real cluster time and changes selection alone.

Estimator is a weighted MEDIAN, then Winsorised twice (drop tracks >100 ps from
the running estimate and recompute, then >60 ps and recompute) -- the same
closed-form recipe as the ML pipeline's trunc-tagw, with w_t in place of p_t.
The median rather than the inverse-variance mean because ~13.5% of genuine HS
tracks sit >200 ps from truth against a ~25 ps quoted resolution, and the
inverse-variance mean is the estimator most exposed to that tail.
"""
import argparse
import numpy as np
import pandas as pd
import uproot

EVT = ["sample_id", "file_idx", "event_num"]
KEY = EVT + ["cluster_idx"]
A_DZ = 0.7          # TRACK_DZ_WEIGHT, as validated for TRKPTZ_TZ


def wmedian(df, wcol, vcol):
    """Weighted median of `vcol` per event. Sorted-cumulative-weight crossing of
    half the per-event total; no interpolation, so the result is always an
    observed track time."""
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
    """One truncation pass: drop tracks more than `win` from the current centre,
    recompute the weighted mean of what survives, fall back to the centre for
    events left with nothing."""
    j = t.merge(centre.rename("c0").rename_axis(EVT).reset_index(), on=EVT, how="left")
    j = j[(j["time"] - j["c0"]).abs() < win]
    g = j.groupby(EVT, sort=False)
    num = g.apply(lambda x: (x[wcol] * x["time"]).sum(), include_groups=False)
    den = g[wcol].sum()
    r = (num / den)
    r.index.names = EVT
    return r.reindex(centre.index).fillna(centre)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--file", required=True)
    ap.add_argument("--tag", default="")
    a = ap.parse_args()

    c = uproot.open(a.file)["clusters"].arrays(
        KEY + ["cluster_time", "delta_z", "delta_t"], library="pd")
    t = uproot.open(a.file)["tracks"].arrays(
        KEY + ["pt", "time", "timeRes", "z0_pull_pv", "sigma_z0"], library="pd")
    t = t.dropna(subset=["z0_pull_pv", "sigma_z0", "time", "timeRes"])
    t = t[t["timeRes"] > 0]
    t["dz"] = t["z0_pull_pv"].abs() * t["sigma_z0"]
    t["w"] = t["pt"] * np.exp(-A_DZ * t["dz"])          # the TRKPTZ_TZ weight
    t["wiv"] = t["w"] / t["timeRes"] ** 2                # + inverse-variance

    c["ok"] = c["delta_t"].abs() < 60
    c["dzcl"] = np.exp(-1.5 * c["delta_z"].abs())
    c["t_truth"] = c["cluster_time"] - c["delta_t"]
    truth = c.groupby(EVT, sort=False)["t_truth"].first()
    truth.index.names = EVT
    nev = len(truth)
    print(f"\n===== {a.tag or a.file}  |  {nev:,} events, {len(c):,} clusters, "
          f"{len(t):,} timed tracks")

    # ---- baselines --------------------------------------------------------
    def cf_sel(extra=None):
        t["_s"] = t["w"]
        s = t.groupby(KEY, sort=False)["_s"].sum().rename("s").reset_index()
        m = c.merge(s, on=KEY, how="left")
        m["s"] = m["s"].fillna(0.0) * m["dzcl"]
        if extra is not None:
            m["s"] = m["s"] * extra(m)
        pick = m.loc[m.groupby(EVT, sort=False)["s"].idxmax()]
        return 100.0 * pick["ok"].mean()

    t["_p"] = t["pt"]
    sp = t.groupby(KEY, sort=False)["_p"].sum().rename("s").reset_index()
    mb = c.merge(sp, on=KEY, how="left")
    mb["s"] = mb["s"].fillna(0.0) * mb["dzcl"]
    base = 100.0 * mb.loc[mb.groupby(EVT, sort=False)["s"].idxmax()]["ok"].mean()
    tz = cf_sel()
    print(f"\n  TRKPTZ                                   {base:6.2f}%")
    print(f"  TRKPTZ_TZ  (per-track dz weight)         {tz:6.2f}%   {tz - base:+.2f}")

    # ---- A) classical aggregate t0 as a DIRECT answer ---------------------
    med = wmedian(t, "w", "time")
    t0 = winsor(t, winsor(t, med, 100.0, "wiv"), 60.0, "wiv")
    okm = ((med.reindex(truth.index) - truth).abs() < 60).fillna(False)
    ok0 = ((t0.reindex(truth.index) - truth).abs() < 60).fillna(False)
    # inverse-variance mean over all tracks, as the estimator control
    ivm = (t.assign(x=t.wiv * t.time).groupby(EVT, sort=False)["x"].sum()
           / t.groupby(EVT, sort=False)["wiv"].sum())
    ivm.index.names = EVT
    oki = ((ivm.reindex(truth.index) - truth).abs() < 60).fillna(False)
    print(f"\n  --- A) the aggregate as the ANSWER (no cluster chosen at all)")
    print(f"  weighted median                          {100*okm.mean():6.2f}%")
    print(f"  + double Winsorisation                   {100*ok0.mean():6.2f}%")
    print(f"  inv-variance mean over all tracks        {100*oki.mean():6.2f}%"
          f"   <- estimator control")

    # ---- B) the aggregate as a RANKING term -------------------------------
    print(f"\n  --- B) the aggregate used only to RANK clusters")
    cc = c.merge(t0.rename("t0").rename_axis(EVT).reset_index(), on=EVT, how="left")
    dt = (cc["cluster_time"] - cc["t0"]).abs()
    for tau in (20.0, 40.0, 60.0, 100.0, 150.0, 250.0):
        w = np.exp(-dt / tau).fillna(1.0).to_numpy()
        m = cc.copy()
        t["_s"] = t["w"]
        s = t.groupby(KEY, sort=False)["_s"].sum().rename("s").reset_index()
        m = m.merge(s, on=KEY, how="left")
        m["s"] = m["s"].fillna(0.0) * m["dzcl"] * w
        cf = 100.0 * m.loc[m.groupby(EVT, sort=False)["s"].idxmax()]["ok"].mean()
        print(f"  x exp(-|t_clus - t0| / {tau:5.0f} ps)          {cf:6.2f}%"
              f"   {cf - tz:+.2f} vs TZ")

    # union: how often is at least one of (selected cluster, aggregate) right
    t["_s"] = t["w"]
    s = t.groupby(KEY, sort=False)["_s"].sum().rename("s").reset_index()
    m = c.merge(s, on=KEY, how="left"); m["s"] = m["s"].fillna(0.0) * m["dzcl"]
    pick = m.loc[m.groupby(EVT, sort=False)["s"].idxmax()].set_index(
        pd.MultiIndex.from_frame(m.loc[m.groupby(EVT, sort=False)["s"].idxmax()][EVT]))
    pick.index.names = EVT
    u = (pick["ok"].reindex(truth.index).fillna(False) | ok0)
    print(f"\n  union(TRKPTZ_TZ cluster, classical t0)   {100*u.mean():6.2f}%"
          f"   <- ceiling for any chooser between the two")


if __name__ == "__main__":
    main()
