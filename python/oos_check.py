#!/usr/bin/env python3
"""
oos_check.py -- out-of-sample check of the final TZP configuration
(alpha, beta, gamma) = (1.0, 0.6, 0.45), guard >= 3.

Every constant was tuned on the four novbs mu=200 samples. Two held-out tests:

  mjj500p0 (mu=200): the VBS-enriched selection. Honest framing: these EVENTS
      overlap the novbs tuning sets (mjj500 is a stricter selection, so a
      subset), but the population MIX is very different -- this tests that the
      constants survive reweighting toward VBS topology, not statistical
      independence.
  mu0 (novbs): genuinely disjoint events at zero pileup. TRKPTZ is ~99.9%
      there; the only question is NO HARM -- a per-track z-weight, a softened
      cluster term and a d0 power must not break events with a single vertex.

Pass criteria: TZP - TRKPTZ >= 0 on every mjj500 sample; |TZP - TRKPTZ| within
noise on every mu0 sample.
"""
import argparse, gc, os
import numpy as np
import pandas as pd
import uproot

EVT = ["sample_id", "file_idx", "event_num"]
KEY = EVT + ["cluster_idx"]
ALPHA, BETA, GAMMA = 1.0, 0.6, 0.45
GUARD, PASS = 3, 60.0


def run(path):
    c = uproot.open(path)["clusters"].arrays(
        KEY + ["cluster_time", "delta_z", "delta_t", "cluster_d0_sigma",
               "trkptz_score", "waves_score"], library="pd")
    t = uproot.open(path)["tracks"].arrays(
        KEY + ["pt", "time", "timeRes", "z0_pull_pv", "sigma_z0",
               "dr_nearest_fwdjet"], library="pd")
    t = t.dropna(subset=["time", "timeRes", "z0_pull_pv", "sigma_z0"])
    t = t[t["timeRes"] > 0]
    t["dz"] = t["z0_pull_pv"].abs() * t["sigma_z0"]
    t["wa"] = t["pt"] * np.exp(-ALPHA * t["dz"])
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
    s = t.groupby(KEY, sort=False)["wa"].sum().rename("sa").reset_index()
    c = c.merge(s, on=KEY, how="left")
    c["tzp"] = (c["sa"].fillna(0.0) * np.exp(-BETA * c["delta_z"].abs())
                * np.power(1.0 / np.maximum(c["cluster_d0_sigma"].to_numpy(float), 1e-12), GAMMA))
    c["t_out"] = np.where(c["n_inj"] >= GUARD, c["t_inj"], c["cluster_time"])
    y = truth.to_numpy()

    out = {"n": len(truth)}
    for lab, scol, tcol in (("TRKPTZ", "trkptz_score", "cluster_time"),
                            ("TZP", "tzp", "t_out"),
                            ("WAVES", "waves_score", "t_inj")):
        ip = c.groupby(EVT, sort=False)[scol].idxmax()
        p = c.loc[ip].set_index(pd.MultiIndex.from_frame(c.loc[ip][EVT]))
        p.index.names = EVT; p = p.reindex(truth.index)
        out[lab] = 100.0 * (np.abs(p[tcol].to_numpy() - y) < PASS).mean()
        # binomial error on the delta is ~ sqrt(2 p q / n); report per-method err
        out[f"e_{lab}"] = 100.0 * np.sqrt(out[lab] / 100 * (1 - out[lab] / 100) / len(truth))
    del c, t; gc.collect()
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dir", default="/Users/mcard/project/vertex_timing/training_new")
    a = ap.parse_args()
    GROUPS = [
        ("mjj500 (mu=200): VBS-enriched reweighting of the tuning events",
         [("vbf", "vbf_mjj500p0"), ("zjets", "zjets_mjj500p0"),
          ("dijet", "dijet_mjj500p0"), ("ttbar", "ttbar_mjj500p0")]),
        ("mu0 (novbs): disjoint zero-pileup controls -- no-harm test",
         [("vbf_mu0", "vbf_mu0_novbs"), ("zeejets_mu0", "zeejets_mu0_novbs"),
          ("ttbar_mu0", "ttbar_mu0_novbs")]),
    ]
    for title, rows in GROUPS:
        print(f"\n===== {title}")
        print(f"  {'sample':<12s} {'events':>9s} {'TRKPTZ':>9s} {'TZP':>9s} "
              f"{'WAVES':>9s} {'TZP-TRKPTZ':>11s} {'+-':>5s}")
        for lab, stem in rows:
            path = f"{a.dir}/{stem}_training.root"
            if not os.path.exists(path):
                print(f"  {lab:<12s}  (file missing: {stem})"); continue
            r = run(path)
            err = np.sqrt(r["e_TRKPTZ"] ** 2 + r["e_TZP"] ** 2)
            print(f"  {lab:<12s} {r['n']:9,d} {r['TRKPTZ']:8.2f}% {r['TZP']:8.2f}% "
                  f"{r['WAVES']:8.2f}% {r['TZP']-r['TRKPTZ']:+11.2f} {err:5.2f}")


if __name__ == "__main__":
    main()
