#!/usr/bin/env python3
"""
timesource_oracle.py -- DIAGNOSTIC ONLY. For the cluster a real selector chose,
compute all three candidate times (raw / in-jet / out-of-jet) and report the one
closest to truth. This is a truth-using oracle and is NOT a method; it measures
how much the choice of TIME SOURCE could be worth if something could make it
per-event.

Selection is held fixed at TRKPTZ_TZ throughout, so the cluster never changes --
only which subset of its tracks the reported time is built from. That also means
the SELECTED-CLUSTER purity is identical across all three rows by construction;
what differs is the purity of the track subset actually used for the time, which
is reported separately.

Includes a matched-noise NULL. Best-of-N over N candidates is upward-biased even
when the candidates carry no complementary information, and this codebase has
been burned by exactly that before (a naive best-cluster oracle read 92.6%
against a 53.4% null). The null here replaces the two alternative times with
Gaussian perturbations of the raw time whose spreads match the observed
raw-to-in-jet and raw-to-out-jet spreads, then takes the same best-of-3. Any
oracle gain above the null is real complementarity; the rest is the statistic.
"""
import argparse
import numpy as np
import pandas as pd
import uproot

EVT = ["sample_id", "file_idx", "event_num"]
KEY = EVT + ["cluster_idx"]
A_DZ = 0.7
PASS = 60.0


def stats(dt):
    """dt = t - t_truth, one entry per event."""
    d = dt[np.isfinite(dt)]
    core = d[np.abs(d) < PASS]
    return dict(
        n=len(d),
        core_frac=100.0 * len(core) / len(d),
        core_rms=float(np.sqrt(np.mean(core ** 2))) if len(core) else float("nan"),
        rms=float(np.sqrt(np.mean(d ** 2))),
        # robust width: half the central 68% interval, insensitive to the tail
        iqr68=float(0.5 * (np.percentile(d, 84.135) - np.percentile(d, 15.865))),
        median_abs=float(np.median(np.abs(d))))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--file", required=True)
    ap.add_argument("--tag", default="")
    a = ap.parse_args()

    c = uproot.open(a.file)["clusters"].arrays(
        KEY + ["cluster_time", "delta_z", "delta_t", "truth_n_hs_tracks"], library="pd")
    t = uproot.open(a.file)["tracks"].arrays(
        KEY + ["pt", "time", "timeRes", "z0_pull_pv", "sigma_z0",
               "dr_nearest_fwdjet", "truth_is_hs"], library="pd")
    t = t.dropna(subset=["z0_pull_pv", "sigma_z0", "time", "timeRes"])
    t = t[t["timeRes"] > 0]
    t["dz"] = t["z0_pull_pv"].abs() * t["sigma_z0"]
    t["w"] = t["pt"] * np.exp(-A_DZ * t["dz"])
    t["iv"] = 1.0 / t["timeRes"] ** 2
    t["inj"] = (t["dr_nearest_fwdjet"] < 0.4).fillna(False)
    t["pths"] = t["pt"] * (t["truth_is_hs"] > 0.5)

    c["t_truth"] = c["cluster_time"] - c["delta_t"]
    c["dzcl"] = np.exp(-1.5 * c["delta_z"].abs())
    truth = c.groupby(EVT, sort=False)["t_truth"].first(); truth.index.names = EVT

    # per-cluster subset times + the purity of the tracks each one uses
    mi = pd.MultiIndex.from_frame(c[KEY])
    for lab, mask in (("raw", None), ("inj", t["inj"]), ("outj", ~t["inj"])):
        x = t if mask is None else t[mask]
        g = x.groupby(KEY, sort=False)
        tt = (g.apply(lambda d: (d["iv"] * d["time"]).sum() / d["iv"].sum(),
                      include_groups=False).reindex(mi).to_numpy())
        pu = (g["pths"].sum() / g["pt"].sum()).reindex(mi).to_numpy()
        c[f"t_{lab}"] = np.where(np.isnan(tt), c["cluster_time"], tt)
        c[f"pur_{lab}"] = pu
    # cluster-level purity (same for all three by construction)
    gall = t.groupby(KEY, sort=False)
    c["pur_cluster"] = (gall["pths"].sum() / gall["pt"].sum()).reindex(mi).to_numpy()

    # TRKPTZ_TZ selection
    t["_s"] = t["w"]
    s = t.groupby(KEY, sort=False)["_s"].sum().rename("s").reset_index()
    c = c.merge(s, on=KEY, how="left")
    c["sc"] = c["s"].fillna(0.0) * c["dzcl"]
    ip = c.groupby(EVT, sort=False)["sc"].idxmax()
    pk = c.loc[ip].set_index(pd.MultiIndex.from_frame(c.loc[ip][EVT]))
    pk.index.names = EVT
    pk = pk.reindex(truth.index)
    y = truth.to_numpy()

    D = {k: pk[f"t_{k}"].to_numpy() - y for k in ("raw", "inj", "outj")}
    A3 = np.vstack([D["raw"], D["inj"], D["outj"]])
    win = np.argmin(np.abs(A3), axis=0)
    D["ORACLE3"] = A3[win, np.arange(A3.shape[1])]
    A2 = np.vstack([D["raw"], D["inj"]])
    D["ORACLE2 (raw+inj)"] = A2[np.argmin(np.abs(A2), axis=0), np.arange(A2.shape[1])]

    # matched-noise null
    rng = np.random.default_rng(0)
    s1 = np.nanstd(D["inj"] - D["raw"]); s2 = np.nanstd(D["outj"] - D["raw"])
    N3 = np.vstack([D["raw"],
                    D["raw"] + rng.normal(0, s1, len(y)),
                    D["raw"] + rng.normal(0, s2, len(y))])
    D["NULL best-of-3"] = N3[np.argmin(np.abs(N3), axis=0), np.arange(N3.shape[1])]

    print(f"\n===== {a.tag or a.file}   {len(y):,} events, selection fixed at TRKPTZ_TZ")
    print(f"  matched-noise null spreads: raw->in-jet {s1:.1f} ps, raw->out-jet {s2:.1f} ps\n")
    print(f"  {'time source':<20s} {'core frac':>10s} {'core RMS':>9s} {'68% half':>9s} "
          f"{'|dt| med':>9s} {'full RMS':>10s}")
    order = ["raw", "inj", "outj", "ORACLE2 (raw+inj)", "ORACLE3", "NULL best-of-3"]
    S = {}
    for k in order:
        st = stats(D[k]); S[k] = st
        print(f"  {k:<20s} {st['core_frac']:9.2f}% {st['core_rms']:8.2f} "
              f"{st['iqr68']:8.2f} {st['median_abs']:8.2f} {st['rms']:9.1f}")
    print(f"\n  oracle3 gain over raw      : {S['ORACLE3']['core_frac']-S['raw']['core_frac']:+.2f}")
    print(f"  null  gain over raw        : {S['NULL best-of-3']['core_frac']-S['raw']['core_frac']:+.2f}"
          f"   <- pure best-of-N inflation")
    print(f"  REAL complementarity       : "
          f"{S['ORACLE3']['core_frac']-S['NULL best-of-3']['core_frac']:+.2f}")

    print(f"\n  purity of the tracks used for the time (pT-weighted HS fraction):")
    for lab, nm in (("cluster", "whole selected cluster"), ("raw", "raw = whole cluster"),
                    ("inj", "in-jet subset"), ("outj", "out-of-jet subset")):
        v = pk[f"pur_{lab}"].to_numpy()
        print(f"    {nm:<24s} {100*np.nanmean(v):6.2f}%")
    print(f"    (selected-cluster purity is identical across the three rows by"
          f" construction -- only the SUBSET used differs)")

    print(f"\n  which source wins the oracle: "
          + "  ".join(f"{n} {100*(win==i).mean():.1f}%"
                      for i, n in enumerate(("raw", "in-jet", "out-jet"))))

    # ---- differential vs event-level n HS tracks --------------------------
    nhs = c.groupby(EVT, sort=False)["truth_n_hs_tracks"].sum().reindex(truth.index).to_numpy()
    print(f"\n  differential vs event-level n forward HS tracks")
    print(f"    {'bin':>7s} {'events':>7s} | {'raw':>7s} {'in-jet':>7s} {'out-jet':>7s} "
          f"{'ORACLE3':>8s} {'NULL':>7s} | {'real':>6s} | {'win: raw/inj/outj':>20s}")
    for lo, hi in [(0,1),(1,2),(2,3),(3,4),(4,5),(5,6),(6,8),(8,10),(10,14),(14,20),(20,10**9)]:
        m = (nhs >= lo) & (nhs < hi)
        if m.sum() < 200: continue
        lab = f"{lo}" if hi == lo+1 else (f">={lo}" if hi > 10**8 else f"{lo}-{hi-1}")
        cf = {k: stats(D[k][m])["core_frac"] for k in
              ("raw", "inj", "outj", "ORACLE3", "NULL best-of-3")}
        w = win[m]
        print(f"    {lab:>7s} {m.sum():7d} | {cf['raw']:6.2f}% {cf['inj']:6.2f}% "
              f"{cf['outj']:6.2f}% {cf['ORACLE3']:7.2f}% {cf['NULL best-of-3']:6.2f}% | "
              f"{cf['ORACLE3']-cf['NULL best-of-3']:+5.2f} | "
              f"{100*(w==0).mean():5.1f}/{100*(w==1).mean():4.1f}/{100*(w==2).mean():4.1f}")


if __name__ == "__main__":
    main()
