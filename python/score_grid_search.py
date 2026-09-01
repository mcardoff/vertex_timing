#!/usr/bin/env python3
"""
score_grid_search.py -- joint grid search over the three score parameters,

    S(C) = [SUM_t pT_t e^{-ALPHA |z0_t - z_PV|}] * e^{-BETA |z_C - z_PV|}
           * (1/sigma_d0)^{GAMMA}

with the guarded in-jet timing held fixed (it does not depend on the score).
ALPHA = 0.7 and GAMMA = 1.2 were each tuned 1-D with the others at their
defaults; BETA = 1.5 is inherited from original TRKPTZ and has never been
scanned at all. A joint optimum can differ from the product of 1-D optima if
the terms interact.

Objective: the WORST-SAMPLE gain over {vbf, zjets, dijet, ttbar} (novbs)
relative to the current point (0.7, 1.5, 1.2) -- one global parameter set,
per the sample-independence requirement. Per-sample optima are reported too,
as the measure of how much a per-sample tune would buy (and therefore of how
much the constraint costs).

Overfitting guard: the 1-D slices through the chosen point are printed so a
sharp (non-plateau) optimum is visible rather than silently trusted.
"""
import argparse, gc
import numpy as np
import pandas as pd
import uproot

EVT = ["sample_id", "file_idx", "event_num"]
KEY = EVT + ["cluster_idx"]
GUARD, PASS = 3, 60.0
ALPHAS = (0.35, 0.5, 0.7, 1.0, 1.4)
BETAS = (0.75, 1.0, 1.5, 2.0, 3.0)
GAMMAS = (0.6, 0.9, 1.2, 1.6, 2.0)
CUR = (0.7, 1.5, 1.2)


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
    t["iv"] = 1.0 / t["timeRes"] ** 2
    t["inj"] = (t["dr_nearest_fwdjet"] < 0.4).fillna(False)

    c["t_truth"] = c["cluster_time"] - c["delta_t"]
    truth_ev = c.groupby(EVT, sort=False)["t_truth"].first()
    mi = pd.MultiIndex.from_frame(c[KEY])
    g = t[t["inj"]].groupby(KEY, sort=False)
    tij = g.apply(lambda x: (x["iv"] * x["time"]).sum() / x["iv"].sum(),
                  include_groups=False).reindex(mi).to_numpy()
    t_inj = np.where(np.isnan(tij), c["cluster_time"], tij)
    n_inj = g.size().reindex(mi).fillna(0).to_numpy()
    t_out = np.where(n_inj >= GUARD, t_inj, c["cluster_time"].to_numpy())
    tmi = c[EVT].merge(truth_ev.rename("y").reset_index(), on=EVT, how="left")["y"].to_numpy()
    ok_out = np.abs(t_out - tmi) < PASS

    # per-cluster ingredients: track sums per alpha via bincount on factorized keys
    ccode, cuq = pd.factorize(pd.MultiIndex.from_frame(t[KEY]), sort=False)
    ncl_t = len(cuq)
    sums = {}
    pt = t["pt"].to_numpy(float); dz = t["dz"].to_numpy(float)
    for a_ in ALPHAS:
        s = np.bincount(ccode, weights=pt * np.exp(-a_ * dz), minlength=ncl_t)
        sums[a_] = pd.Series(s, index=cuq).reindex(mi).fillna(0.0).to_numpy()
    del t; gc.collect()

    ecode, _ = pd.factorize(pd.MultiIndex.from_frame(c[EVT]), sort=False)
    dzc = np.abs(c["delta_z"].to_numpy(float))
    lp = np.log(1.0 / np.maximum(c["cluster_d0_sigma"].to_numpy(float), 1e-12))
    nev = int(ecode.max()) + 1
    # segment ends for last-of-segment picking (events are contiguous in file order)
    ends = np.r_[np.flatnonzero(np.diff(ecode)) , len(ecode) - 1]
    return dict(sums=sums, dzc=dzc, lp=lp, ecode=ecode, ends=ends,
                ok=ok_out, nev=nev)


def core(D, a_, b_, g_):
    with np.errstate(divide="ignore"):
        L = np.log(np.maximum(D["sums"][a_], 1e-300)) - b_ * D["dzc"] + g_ * D["lp"]
    # argmax per event: stable lexsort by (event, score); last row of each
    # segment is that event's maximum
    perm = np.lexsort((L, D["ecode"]))
    picked = perm[D["ends"]]
    return 100.0 * D["ok"][picked].mean()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dir", default="/Users/mcard/project/vertex_timing/training_new")
    ap.add_argument("--samples", default="vbf,zjets,dijet,ttbar")
    a = ap.parse_args()
    smps = a.samples.split(",")
    D = {}
    for smp in smps:
        print(f"  loading {smp} ...", flush=True)
        D[smp] = load(f"{a.dir}/{smp}_novbs_training.root")

    R = {}
    for a_ in ALPHAS:
        for b_ in BETAS:
            for g_ in GAMMAS:
                R[(a_, b_, g_)] = {s: core(D[s], a_, b_, g_) for s in smps}
        print(f"  alpha={a_} done", flush=True)

    cur = R[CUR]
    print(f"\n  current point (alpha,beta,gamma) = {CUR}: "
          + "  ".join(f"{s} {cur[s]:.2f}%" for s in smps))

    rows = [(k, min(v[s] - cur[s] for s in smps), np.mean([v[s] - cur[s] for s in smps]))
            for k, v in R.items()]
    rows.sort(key=lambda r: -r[1])
    print(f"\n  TOP 12 by WORST-sample gain vs current  (a, b, g)")
    print(f"  {'a':>5s} {'b':>5s} {'g':>5s} " + " ".join(f"{s:>8s}" for s in smps)
          + f" {'worst':>7s} {'mean':>7s}")
    for k, worst, mean in rows[:12]:
        v = R[k]
        print(f"  {k[0]:5.2f} {k[1]:5.2f} {k[2]:5.2f} "
              + " ".join(f"{v[s]-cur[s]:+7.2f}" for s in smps)
              + f" {worst:+7.2f} {mean:+7.2f}")

    print(f"\n  PER-SAMPLE optima (what a per-sample tune would buy)")
    for s in smps:
        k = max(R, key=lambda k: R[k][s])
        print(f"    {s:6s} best {k} at {R[k][s]:.2f}%  (+{R[k][s]-cur[s]:.2f} vs current)")

    best = rows[0][0]
    print(f"\n  1-D slices through the joint optimum {best} (worst-sample gain)")
    for i, (name, vals) in enumerate((("alpha", ALPHAS), ("beta", BETAS), ("gamma", GAMMAS))):
        line = f"    {name:>5s}: "
        for x in vals:
            k = list(best); k[i] = x; k = tuple(k)
            w = min(R[k][s] - cur[s] for s in smps)
            line += f"{x:g}:{w:+.2f}  "
        print(line)


if __name__ == "__main__":
    main()
