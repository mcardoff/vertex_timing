#!/usr/bin/env python3
"""
jet_fraction_term.py -- add a BOUNDED jet-association term to the score.

missing_information.py showed that on the events we lose, our own inputs carry
nothing (delta_z AUC 0.502, sumpt 0.466) while every jet-association variable
separates at 0.86-0.90, against a 0.978 truth ceiling. So the missing information
is jet association -- which WAVeS already uses, and which is exactly what makes
WAVeS fail on Z+jets.

The hypothesis is that the problem is WAVeS's FUNCTIONAL FORM rather than the
information. WAVeS multiplies by jet pT (canonical: pT_trk^2 pT_jet^2 / dR), so a
single high-pT pileup jet can dominate the whole score -- the documented "captured
by the fake" mode. A BOUNDED fraction cannot do that: frac_pt_in_fwdjet lives in
[0, 1], so a jet can at most multiply a cluster's score by (1 + alpha), never
scale it without limit.

    score = SUM_t pT_t exp(-0.7|dz_t|) * exp(-1.5|dz_cl|) * (1 + alpha * f_jet)

alpha = 0 is TRKPTZ_TZ exactly. Scanned on all four mu=200 samples with ONE alpha,
since a per-sample alpha is not implementable (Z+jets is the VBF H->inv
background).
"""
import argparse, gc
import numpy as np
import pandas as pd
import uproot

EVT = ["sample_id", "file_idx", "event_num"]
KEY = EVT + ["cluster_idx"]
A_DZ, GUARD, PASS = 0.7, 3, 60.0


def load(path):
    c = uproot.open(path)["clusters"].arrays(
        KEY + ["cluster_time", "delta_z", "delta_t", "waves_score",
               "frac_pt_in_fwdjet", "n_tracks_in_fwdjet", "min_dr_to_fwdjet"],
        library="pd")
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
    g = t[t["inj"]].groupby(KEY, sort=False)
    tij = g.apply(lambda x: (x["iv"] * x["time"]).sum() / x["iv"].sum(),
                  include_groups=False).reindex(mi).to_numpy()
    c["t_inj"] = np.where(np.isnan(tij), c["cluster_time"], tij)
    c["n_inj"] = g.size().reindex(mi).fillna(0).to_numpy()
    t["_s"] = t["w"]
    s = t.groupby(KEY, sort=False)["_s"].sum().rename("s").reset_index()
    c = c.merge(s, on=KEY, how="left")
    c["base"] = c["s"].fillna(0.0) * c["dzcl"]
    c["t_out"] = np.where(c["n_inj"] >= GUARD, c["t_inj"], c["cluster_time"])
    c["fjet"] = c["frac_pt_in_fwdjet"].fillna(0.0).clip(0, 1)
    del t; gc.collect()
    return c, truth


def cf(c, truth, alpha, col="fjet"):
    c["_sc"] = c["base"] * (1.0 + alpha * c[col])
    ip = c.groupby(EVT, sort=False)["_sc"].idxmax()
    p = c.loc[ip].set_index(pd.MultiIndex.from_frame(c.loc[ip][EVT]))
    p.index.names = EVT
    return 100.0 * (np.abs(p["t_out"].reindex(truth.index).to_numpy()
                           - truth.to_numpy()) < PASS).mean()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dir", default="/Users/mcard/project/vertex_timing/training_new")
    ap.add_argument("--samples", default="vbf,zjets,dijet,ttbar")
    ap.add_argument("--sel", default="novbs")
    a = ap.parse_args()
    smps = a.samples.split(",")
    ALPHAS = [0.0, 0.25, 0.5, 1.0, 2.0, 4.0, 8.0, 16.0]
    res, wav = {}, {}
    for smp in smps:
        print(f"  loading {smp} ...", flush=True)
        c, truth = load(f"{a.dir}/{smp}_{a.sel}_training.root")
        res[smp] = {al: cf(c, truth, al) for al in ALPHAS}
        c["_sc"] = c["waves_score"]
        ip = c.groupby(EVT, sort=False)["_sc"].idxmax()
        p = c.loc[ip].set_index(pd.MultiIndex.from_frame(c.loc[ip][EVT])); p.index.names = EVT
        wav[smp] = 100.0 * (np.abs(p["t_inj"].reindex(truth.index).to_numpy()
                                   - truth.to_numpy()) < PASS).mean()
        del c; gc.collect()

    print(f"\n===== score = base * (1 + alpha * frac_pt_in_fwdjet), guarded in-jet time")
    print(f"  alpha = 0 is TRKPTZ_TZ + guard exactly\n")
    print(f"  {'alpha':>6s} " + " ".join(f"{s:>9s}" for s in smps) + "   worst delta")
    for al in ALPHAS:
        d = [res[s][al] - res[s][0.0] for s in smps]
        print(f"  {al:6.2f} " + " ".join(f"{res[s][al]:8.2f}%" for s in smps)
              + f"   {min(d):+6.2f}" + ("   safe on all" if min(d) >= -0.05 else ""))
    print(f"\n  {'deployed WAVeS':>6s} " + " ".join(f"{wav[s]:8.2f}%" for s in smps))
    best = max(ALPHAS, key=lambda al: min(res[s][al] - res[s][0.0] for s in smps)
               if min(res[s][al] - res[s][0.0] for s in smps) >= -0.05 else -99)
    print(f"\n  best alpha safe on every sample: {best}")
    print(f"  vs deployed WAVeS: " + " ".join(f"{s} {res[s][best]-wav[s]:+.2f}" for s in smps))


if __name__ == "__main__":
    main()
