#!/usr/bin/env python3
"""
guarded_allsamples.py -- the full selection x timing grid on every mu=200 sample,
to test whether the sample-independent guard generalises past VBF and Z+jets.

Two questions:
  1. Does the >= 3 in-jet-track guard stay >= the raw baseline on dijet and ttbar
     too? A rule safe on two samples is not demonstrated safe on four.
  2. Does the guard help WAVeS? Deployed WAVeS applies its in-jet re-timing
     UNCONDITIONALLY, and that re-timing was measured at +0.82 on VBF and -1.43
     on Z+jets -- so WAVeS should inherit the same Z+jets deficit, and guarding
     it should recover that without costing VBF much.

Selection and timing are varied independently, so every cell is comparable:
selections are {TRKPTZ, TRKPTZ_TZ, WAVeS}, timings are {raw, in-jet always,
in-jet guarded at >= N}. Deployed WAVeS is the (WAVeS, in-jet always) cell.
"""
import argparse, gc
import numpy as np
import pandas as pd
import uproot

EVT = ["sample_id", "file_idx", "event_num"]
KEY = EVT + ["cluster_idx"]
A_DZ = 0.7
GUARD = 3


def load(path):
    c = uproot.open(path)["clusters"].arrays(
        KEY + ["cluster_time", "delta_z", "delta_t", "trkptz_score", "waves_score"],
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
    ij = t[t["inj"]]
    g = ij.groupby(KEY, sort=False)
    tij = (g.apply(lambda x: (x["iv"] * x["time"]).sum() / x["iv"].sum(),
                   include_groups=False).reindex(mi).to_numpy())
    c["t_inj"] = np.where(np.isnan(tij), c["cluster_time"], tij)
    c["n_inj"] = g.size().reindex(mi).fillna(0).to_numpy()
    t["_s"] = t["w"]
    s = t.groupby(KEY, sort=False)["_s"].sum().rename("s").reset_index()
    c = c.merge(s, on=KEY, how="left")
    c["tz_score"] = c["s"].fillna(0.0) * c["dzcl"]
    frac_inj = float(t["inj"].mean())
    del t; gc.collect()
    return c, truth, frac_inj


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dir", required=True)
    ap.add_argument("--samples", default="vbf,zjets,dijet,ttbar")
    ap.add_argument("--sel", default="novbs")
    a = ap.parse_args()

    SEL = [("TRKPTZ", "trkptz_score"), ("TRKPTZ_TZ", "tz_score"), ("WAVeS", "waves_score")]
    rows, fracs, nev = {}, {}, {}
    for smp in a.samples.split(","):
        path = f"{a.dir}/{smp}_{a.sel}_training.root"
        print(f"  loading {smp} ...", flush=True)
        c, truth, fj = load(path)
        fracs[smp] = fj; nev[smp] = len(truth)
        y = truth.to_numpy()
        for sname, scol in SEL:
            ip = c.groupby(EVT, sort=False)[scol].idxmax()
            pk = c.loc[ip].set_index(pd.MultiIndex.from_frame(c.loc[ip][EVT]))
            pk.index.names = EVT; pk = pk.reindex(truth.index)
            raw = pk["cluster_time"].to_numpy()
            inj = pk["t_inj"].to_numpy()
            nij = pk["n_inj"].to_numpy()
            for tname, tv in (("raw", raw),
                              ("in-jet always", inj),
                              (f"in-jet if >={GUARD}", np.where(nij >= GUARD, inj, raw))):
                rows.setdefault((sname, tname), {})[smp] = 100.0 * (np.abs(tv - y) < 60).mean()
        del c; gc.collect()

    smps = a.samples.split(",")
    print(f"\n===== selection x timing, all mu=200 samples ({a.sel} selection)")
    print("  in-jet TRACK fraction: " + "  ".join(f"{s} {100*fracs[s]:.1f}%" for s in smps))
    print("  events:                " + "  ".join(f"{s} {nev[s]:,}" for s in smps))
    print(f"\n  {'selection':<11s} {'timing':<16s} " + " ".join(f"{s:>9s}" for s in smps))
    for sname, _ in SEL:
        base = rows[(sname, "raw")]
        for tname in ("raw", "in-jet always", f"in-jet if >={GUARD}"):
            v = rows[(sname, tname)]
            line = f"  {sname:<11s} {tname:<16s} " + " ".join(f"{v[s]:8.2f}%" for s in smps)
            if tname != "raw":
                line += "   d " + " ".join(f"{v[s]-base[s]:+6.2f}" for s in smps)
            print(line)
        print()

    print("  verdicts")
    for sname, _ in SEL:
        base = rows[(sname, "raw")]
        for tname in ("in-jet always", f"in-jet if >={GUARD}"):
            d = {s: rows[(sname, tname)][s] - base[s] for s in smps}
            worst = min(d.values()); bad = [s for s in smps if d[s] < -0.05]
            print(f"    {sname:<11s} {tname:<16s} worst {worst:+6.2f}"
                  + (f"   LOSES on {','.join(bad)}" if bad else "   safe on all"))


if __name__ == "__main__":
    main()
