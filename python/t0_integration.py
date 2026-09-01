#!/usr/bin/env python3
"""
t0_integration.py -- integrate the classical aggregate t0 into the TZP pipeline,
scanned on all four samples with one global parameter set.

The aggregate (per-track w = pT e^{-|dz|}-weighted median of HGTD times, double
truncation at 100/60 ps with w/sigma_t^2 weights) beats TRKPTZ standalone and
its ML analogue was the learned selector's top feature. Two sample-independent
entry points:

  A. SCORE term: S' = S_TZP * f(|t_C - t0|), f in {Cauchy 1/(1+(dt/tau)^2),
     exp(-dt/tau)}, tau scanned. Circularity caveat: t0 is built from the same
     tracks, and the largest cluster dominates the median, so for the leading
     cluster the term is partly self-agreement; for the others it is real
     corroboration.
  B. TIMING fallback: keep the TZP pick, but REPORT t0 when the picked
     cluster's (guarded) time sits further than X from it. Changes no
     selection; the reported time is no longer any cluster's time in the
     fallback branch.

Both are evaluated against the TZP joint-tune baseline (1.0, 0.6, 0.45) with
guarded in-jet timing, and the union(TZP pick, t0) ceiling is printed for
scale.
"""
import argparse, gc, sys, os
import numpy as np
import pandas as pd
import uproot

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from classical_t0 import wmedian, winsor

EVT = ["sample_id", "file_idx", "event_num"]
KEY = EVT + ["cluster_idx"]
ALPHA, BETA, GAMMA = 1.0, 0.6, 0.45
GUARD, PASS = 3, 60.0


def load(path):
    c = uproot.open(path)["clusters"].arrays(
        KEY + ["cluster_time", "delta_z", "delta_t", "cluster_d0_sigma"],
        library="pd")
    t = uproot.open(path)["tracks"].arrays(
        KEY + ["pt", "time", "timeRes", "z0_pull_pv", "sigma_z0",
               "dr_nearest_fwdjet"], library="pd")
    t = t.dropna(subset=["time", "timeRes", "z0_pull_pv", "sigma_z0"])
    t = t[t["timeRes"] > 0]
    t["dz"] = t["z0_pull_pv"].abs() * t["sigma_z0"]
    t["w"] = t["pt"] * np.exp(-0.7 * t["dz"])        # aggregate keeps ITS OWN tuned weight
    t["wa"] = t["pt"] * np.exp(-ALPHA * t["dz"])     # score's per-track sum
    t["iv"] = 1.0 / t["timeRes"] ** 2
    t["wiv"] = t["w"] * t["iv"]
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
    c["score"] = (c["sa"].fillna(0.0) * np.exp(-BETA * c["delta_z"].abs())
                  * np.power(1.0 / np.maximum(c["cluster_d0_sigma"].to_numpy(float), 1e-12),
                             GAMMA))
    c["t_out"] = np.where(c["n_inj"] >= GUARD, c["t_inj"], c["cluster_time"])

    t0 = winsor(t, winsor(t, wmedian(t, "w", "time"), 100.0, "wiv"), 60.0, "wiv")
    c["t0"] = t0.reindex(pd.MultiIndex.from_frame(c[EVT])).to_numpy()
    del t; gc.collect()
    return c, truth


def pick(c, truth, sc):
    c["_s"] = sc
    ip = c.groupby(EVT, sort=False)["_s"].idxmax()
    p = c.loc[ip].set_index(pd.MultiIndex.from_frame(c.loc[ip][EVT]))
    p.index.names = EVT
    return p.reindex(truth.index)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dir", default="/Users/mcard/project/vertex_timing/training_new")
    ap.add_argument("--samples", default="vbf,zjets,dijet,ttbar")
    a = ap.parse_args()
    smps = a.samples.split(",")
    S = {}
    for smp in smps:
        print(f"  loading {smp} ...", flush=True)
        S[smp] = load(f"{a.dir}/{smp}_novbs_training.root")

    base, unio, okt0 = {}, {}, {}
    P = {}
    for smp in smps:
        c, truth = S[smp]
        p = pick(c, truth, c["score"].to_numpy()); P[smp] = p
        y = truth.to_numpy()
        okb = np.abs(p["t_out"].to_numpy() - y) < PASS
        ok0 = np.abs(p["t0"].to_numpy() - y) < PASS
        base[smp] = 100 * okb.mean(); okt0[smp] = 100 * ok0.mean()
        unio[smp] = 100 * (okb | ok0).mean()
    print("\n  baselines")
    print("    " + " ".join(f"{s:>9s}" for s in smps))
    for lab, v in (("TZP pick", base), ("t0 alone", okt0), ("union", unio)):
        print(f"    {lab:<10s}" + " ".join(f"{v[s]:8.2f}%" for s in smps))

    print("\n  A. SCORE agreement term  S' = S * f(|t_C - t0|)")
    for fname, f in (("cauchy", lambda dt, tau: 1.0 / (1.0 + (dt / tau) ** 2)),
                     ("exp", lambda dt, tau: np.exp(-dt / tau))):
        print(f"     f = {fname}(tau)")
        print(f"     {'tau':>6s} " + " ".join(f"{s:>9s}" for s in smps) + "   worst d")
        for tau in (60, 120, 250, 500, 1000):
            vals = {}
            for smp in smps:
                c, truth = S[smp]
                dt = np.abs(c["cluster_time"].to_numpy() - c["t0"].to_numpy())
                p = pick(c, truth, c["score"].to_numpy() * f(dt, tau))
                vals[smp] = 100 * (np.abs(p["t_out"].to_numpy() - truth.to_numpy()) < PASS).mean()
            d = [vals[s] - base[s] for s in smps]
            print(f"     {tau:6d} " + " ".join(f"{vals[s]:8.2f}%" for s in smps)
                  + f"   {min(d):+6.2f}" + ("   safe on all" if min(d) >= -0.05 else ""))

    print("\n  B. TIMING fallback: report t0 when |t_pick - t0| > X")
    print(f"     {'X':>6s} " + " ".join(f"{s:>9s}" for s in smps)
          + "   worst d   (fires)")
    for X in (60, 100, 150, 250, 400):
        vals, fr = {}, {}
        for smp in smps:
            p = P[smp]; y = S[smp][1].to_numpy()
            far = np.abs(p["t_out"].to_numpy() - p["t0"].to_numpy()) > X
            tv = np.where(far, p["t0"].to_numpy(), p["t_out"].to_numpy())
            vals[smp] = 100 * (np.abs(tv - y) < PASS).mean()
            fr[smp] = 100 * far.mean()
        d = [vals[s] - base[s] for s in smps]
        print(f"     {X:6d} " + " ".join(f"{vals[s]:8.2f}%" for s in smps)
              + f"   {min(d):+6.2f}   " + "/".join(f"{fr[s]:.0f}%" for s in smps)
              + ("   safe on all" if min(d) >= -0.05 else ""))

    print("\n  A+B combined at each one's best (if both have a safe point)")


if __name__ == "__main__":
    main()
