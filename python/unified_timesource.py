#!/usr/bin/env python3
"""
unified_timesource.py -- can ONE per-event rule decide raw-vs-in-jet timing, with
the SAME threshold on every sample?

Z+jets is the background to VBF H->inv, so a per-SAMPLE configuration is not
implementable: signal and background would get different algorithms and any
resulting separation would be manufactured. A per-EVENT function of reconstructed
quantities is legitimate even though its distribution differs by sample -- the
algorithm is fixed, the events differ.

Baselines to beat, TRKPTZ_TZ selection throughout:
                        VBF      Z+jets
  raw (full cluster)   91.49%    63.89%     <- must not lose on either
  in-jet always        92.31%    62.46%     <- +0.82 / -1.43, the sample-dependent
                                               result we are trying to obtain
                                               legitimately

A rule is only interesting if a SINGLE threshold is >= raw on BOTH samples. Every
candidate below is computed from reco alone (jet dR, track times, timeRes).
"""
import argparse
import numpy as np
import pandas as pd
import uproot

EVT = ["sample_id", "file_idx", "event_num"]
KEY = EVT + ["cluster_idx"]
A_DZ = 0.7


def load(path):
    c = uproot.open(path)["clusters"].arrays(
        KEY + ["cluster_time", "delta_z", "delta_t"], library="pd")
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
    g_all = t.groupby(KEY, sort=False)
    ij = t[t["inj"]]
    g_ij = ij.groupby(KEY, sort=False)
    tij = (g_ij.apply(lambda x: (x["iv"] * x["time"]).sum() / x["iv"].sum(),
                      include_groups=False).reindex(mi).to_numpy())
    c["t_inj"] = np.where(np.isnan(tij), c["cluster_time"], tij)
    # per-event-reco quantities describing how much timing information is in-jet
    c["n_inj"] = g_ij.size().reindex(mi).fillna(0).to_numpy()
    c["n_all"] = g_all.size().reindex(mi).fillna(0).to_numpy()
    c["f_inj"] = c["n_inj"] / c["n_all"].clip(lower=1)
    iv_ij = g_ij["iv"].sum().reindex(mi).fillna(0).to_numpy()
    iv_all = g_all["iv"].sum().reindex(mi).replace(0, np.nan).to_numpy()
    c["fw_inj"] = iv_ij / iv_all
    # spread of the in-jet times, a proxy for whether they agree with each other
    c["sd_inj"] = g_ij["time"].std().reindex(mi).to_numpy()

    t["_s"] = t["w"]
    s = t.groupby(KEY, sort=False)["_s"].sum().rename("s").reset_index()
    c = c.merge(s, on=KEY, how="left")
    c["sc"] = c["s"].fillna(0.0) * c["dzcl"]
    ip = c.groupby(EVT, sort=False)["sc"].idxmax()
    pk = c.loc[ip].set_index(pd.MultiIndex.from_frame(c.loc[ip][EVT]))
    pk.index.names = EVT
    return pk.reindex(truth.index), truth


def cf(t, truth):
    return 100.0 * (np.abs(np.asarray(t) - truth.to_numpy()) < 60).mean()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--vbf", required=True)
    ap.add_argument("--zjets", required=True)
    a = ap.parse_args()
    S = {"VBF": load(a.vbf), "Z+jets": load(a.zjets)}

    print(f"\n  {'':<34s} {'VBF':>9s} {'Z+jets':>9s}")
    base = {}
    for k, (pk, tr) in S.items(): base[k] = cf(pk["t_inj"] * 0 + pk["cluster_time"], tr)
    print(f"  {'raw (full cluster)':<34s} "
          + " ".join(f"{base[k]:8.2f}%" for k in S))
    alw = {k: cf(S[k][0]["t_inj"], S[k][1]) for k in S}
    print(f"  {'in-jet always':<34s} "
          + " ".join(f"{alw[k]:8.2f}%" for k in S)
          + "   <- sample-dependent target")
    for k in S: print(f"     {k} in-jet-always delta: {alw[k]-base[k]:+.2f}")

    RULES = [
        ("n_inj",  "n in-jet tracks in cluster >=", [1, 2, 3, 4, 5, 6, 8, 10]),
        ("f_inj",  "in-jet track FRACTION >=",      [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7]),
        ("fw_inj", "in-jet timing-WEIGHT frac >=",  [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7]),
    ]
    for col, label, cuts in RULES:
        print(f"\n  use in-jet time when {label} X, else raw:")
        print(f"    {'X':>6s} {'VBF':>9s} {'d':>7s} {'used%':>6s} | "
              f"{'Z+jets':>9s} {'d':>7s} {'used%':>6s} | verdict")
        for x in cuts:
            row = f"    {x:>6} "
            okboth = True
            for k in S:
                pk, tr = S[k]
                use = (pk[col] >= x).fillna(False).to_numpy()
                v = np.where(use, pk["t_inj"], pk["cluster_time"])
                s_ = cf(v, tr); d = s_ - base[k]
                if d < -0.05: okboth = False
                row += f"{s_:8.2f}% {d:+6.2f} {100*use.mean():5.1f}% | "
            print(row + ("BOTH >= raw" if okboth else ""))


if __name__ == "__main__":
    main()
