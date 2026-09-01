#!/usr/bin/env python3
"""
pick_display_events.py -- select event-display candidates that illustrate a
genuine SELECTION difference, in both directions.

Screening on the time criterion alone (|t - t_truth| < 60) is not enough: it
admits events where the loser picked a lone genuine hard-scatter track carrying a
mis-measured HGTD time. Those are the irreducible ~13.5% mis-measurement tail,
not a selection failure, and they flatter whichever method happens to avoid them.

A genuine selection difference requires all of:
  - the event HAS a right answer: the HS-dominant cluster (most truth HS tracks)
    is itself within 60 ps of truth, so the winner had something real to find;
  - the WINNER picks a passing cluster;
  - the LOSER picks a DIFFERENT cluster that fails AND is not the HS-dominant one
    AND is genuinely pileup-like (few HS tracks, low purity) -- which is what
    rules out the mis-measured-HS-track artifact.

Both directions are produced. WAVeS wins slightly more often than the guarded
method on VBF (1077 vs 1001), so showing only one direction misrepresents it.
"""
import argparse, glob
import numpy as np
import pandas as pd
import uproot

EVT = ["sample_id", "file_idx", "event_num"]
KEY = EVT + ["cluster_idx"]
A_DZ, GUARD, PASS = 0.7, 3, 60.0
MAX_HS_LOSER = 2      # loser's cluster must hold at most this many HS tracks
MAX_PUR_LOSER = 0.50  # ...and be at most this pure


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--file", default="figs/hists/vbf_mjj500p0_training.root")
    ap.add_argument("--ntuple-dir", default="/Users/mcard/project/ntuple-hgtd")
    ap.add_argument("--n", type=int, default=3)
    a = ap.parse_args()

    c = uproot.open(a.file)["clusters"].arrays(
        KEY + ["cluster_time", "delta_z", "delta_t", "waves_score", "trkptz_score",
               "n_tracks", "truth_purity", "truth_n_hs_tracks"], library="pd")
    t = uproot.open(a.file)["tracks"].arrays(
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
    c["tz"] = c["s"].fillna(0.0) * c["dzcl"]
    # each method's reported time for every cluster
    c["t_ours"] = np.where(c["n_inj"] >= GUARD, c["t_inj"], c["cluster_time"])
    c["t_wav"] = c["t_inj"]

    def pick(col):
        ip = c.groupby(EVT, sort=False)[col].idxmax()
        p = c.loc[ip].set_index(pd.MultiIndex.from_frame(c.loc[ip][EVT]))
        p.index.names = EVT
        return p.reindex(truth.index)

    po, pw = pick("tz"), pick("waves_score")
    # the HS-dominant cluster: most truth HS tracks (ties -> most tracks)
    hb = c.sort_values(["truth_n_hs_tracks", "n_tracks"]).groupby(EVT, sort=False).tail(1)
    hb = hb.set_index(pd.MultiIndex.from_frame(hb[EVT])); hb.index.names = EVT
    hb = hb.reindex(truth.index)

    y = truth.to_numpy()
    E = pd.DataFrame(index=truth.index)
    E["fi"] = truth.index.get_level_values(1); E["ev"] = truth.index.get_level_values(2)
    E["truth"] = y
    E["t_o"] = po["t_ours"].to_numpy(); E["t_w"] = pw["t_wav"].to_numpy()
    E["ok_o"] = np.abs(E["t_o"] - y) < PASS
    E["ok_w"] = np.abs(E["t_w"] - y) < PASS
    E["ci_o"] = po["cluster_idx"].to_numpy(); E["ci_w"] = pw["cluster_idx"].to_numpy()
    E["hs_ci"] = hb["cluster_idx"].to_numpy()
    # was a right answer even available?
    E["hs_ok"] = np.abs(hb["t_ours"].to_numpy() - y) < PASS
    E["hs_n"] = hb["truth_n_hs_tracks"].to_numpy()
    for tag, p in (("o", po), ("w", pw)):
        E[f"n_{tag}"] = p["n_tracks"].to_numpy()
        E[f"pur_{tag}"] = p["truth_purity"].to_numpy()
        E[f"hs_{tag}"] = p["truth_n_hs_tracks"].to_numpy()
        E[f"inj_{tag}"] = p["n_inj"].to_numpy()

    print(f"  {len(E):,} events | ours {100*E.ok_o.mean():.2f}%  WAVeS {100*E.ok_w.mean():.2f}%"
          f"  |  ours-only {int((E.ok_o&~E.ok_w).sum())}  WAVeS-only {int((~E.ok_o&E.ok_w).sum())}")

    files = sorted(glob.glob(f"{a.ntuple_dir}/*.root"))
    for wtag, ltag, label in (("o", "w", "OURS wins"), ("w", "o", "WAVeS wins")):
        m = (E[f"ok_{wtag}"] & ~E[f"ok_{ltag}"]
             & E["hs_ok"] & (E["hs_n"] >= 4)                 # a real answer existed
             & (E[f"ci_{wtag}"] != E[f"ci_{ltag}"])          # they disagree on the cluster
             & (E[f"hs_{ltag}"] <= MAX_HS_LOSER)             # loser's cluster is pileup-like
             & (E[f"pur_{ltag}"] <= MAX_PUR_LOSER)
             & (E[f"n_{ltag}"] >= 3)                         # NOT a lone mis-timed track
             & (E[f"n_{wtag}"] >= 5) & (E[f"n_{wtag}"] <= 25))
        D = E[m].copy()
        D["miss"] = np.abs(D[f"t_{ltag}"] - D["truth"])
        D = D.sort_values("miss", ascending=False)
        print(f"\n  === {label}: {len(D)} clean candidates "
              f"(loser picked a genuinely pileup-like cluster)")
        for _, r in D.head(a.n).iterrows():
            f = files[int(r.fi)]
            print(f"    {f.split('/')[-1]}  ev {int(r.ev):5d}")
            print(f"      truth {r.truth:8.1f} | ours {r.t_o:8.1f} (cl {int(r.ci_o)}, "
                  f"{int(r.n_o)} trk, pur {r.pur_o:.2f}, {int(r.hs_o)} HS, {int(r.inj_o)} in-jet)")
            print(f"      {'':13s}| WAVeS {r.t_w:7.1f} (cl {int(r.ci_w)}, "
                  f"{int(r.n_w)} trk, pur {r.pur_w:.2f}, {int(r.hs_w)} HS)")
            print(f"      --extra_time {r['t_o' if wtag=='o' else 't_w']:.1f}")


if __name__ == "__main__":
    main()
