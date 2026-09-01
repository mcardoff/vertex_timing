#!/usr/bin/env python3
"""
failure_modes.py -- why does TRKPTZ_TZ + guarded in-jet lose the events it loses?

Our score is  SUM_t pT_t exp(-0.7|dz_t|) * exp(-1.5|dz_cluster|)  -- purely
pT and z. WAVeS adds jet proximity, which is information our score does not have
at all. The question is whether the events we lose are the ones where z and pT
genuinely cannot separate, i.e. whether the loss is structural rather than a
tuning problem.

Compares, event by event, the cluster WE picked against the cluster WAVeS picked,
restricted to the clean disagreement sets from pick_display_events.py (a right
answer existed; the loser's cluster is genuinely pileup-like). Both directions
are characterised so the comparison is symmetric.
"""
import argparse, glob
import numpy as np
import pandas as pd
import uproot

EVT = ["sample_id", "file_idx", "event_num"]
KEY = EVT + ["cluster_idx"]
A_DZ, GUARD, PASS = 0.7, 3, 60.0


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--file", default="figs/hists/vbf_mjj500p0_training.root")
    ap.add_argument("--ntuple-dir", default="/Users/mcard/project/ntuple-hgtd")
    a = ap.parse_args()

    EX = ["cluster_time", "delta_z", "delta_t", "waves_score", "trkptz_score",
          "n_tracks", "sumpt", "truth_purity", "truth_n_hs_tracks",
          "n_tracks_in_fwdjet", "frac_pt_in_fwdjet", "min_dr_to_fwdjet",
          "pt_of_nearest_fwdjet"]
    c = uproot.open(a.file)["clusters"].arrays(KEY + EX, library="pd")
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
    c["t_ours"] = np.where(c["n_inj"] >= GUARD, c["t_inj"], c["cluster_time"])

    def pick(col):
        ip = c.groupby(EVT, sort=False)[col].idxmax()
        p = c.loc[ip].set_index(pd.MultiIndex.from_frame(c.loc[ip][EVT]))
        p.index.names = EVT
        return p.reindex(truth.index)

    po, pw = pick("tz"), pick("waves_score")
    hb = c.sort_values(["truth_n_hs_tracks", "n_tracks"]).groupby(EVT, sort=False).tail(1)
    hb = hb.set_index(pd.MultiIndex.from_frame(hb[EVT])); hb.index.names = EVT
    hb = hb.reindex(truth.index)
    y = truth.to_numpy()

    ok_o = np.abs(po["t_ours"].to_numpy() - y) < PASS
    ok_w = np.abs(pw["t_inj"].to_numpy() - y) < PASS
    hs_ok = np.abs(hb["t_ours"].to_numpy() - y) < PASS
    diff = po["cluster_idx"].to_numpy() != pw["cluster_idx"].to_numpy()

    print(f"  {len(y):,} events   ours {100*ok_o.mean():.2f}%   WAVeS {100*ok_w.mean():.2f}%")
    print(f"  raw disagreements: ours-only {int((ok_o&~ok_w).sum())}, "
          f"WAVeS-only {int((~ok_o&ok_w).sum())}")

    for wtag, W, L, label in (("w", pw, po, "WAVeS WINS  (we lose)"),
                              ("o", po, pw, "OURS WINS   (WAVeS loses)")):
        win = (ok_w & ~ok_o) if wtag == "w" else (ok_o & ~ok_w)
        m = (win & hs_ok & diff
             & (L["truth_n_hs_tracks"].to_numpy() <= 2)
             & (L["truth_purity"].to_numpy() <= 0.50)
             & (L["n_tracks"].to_numpy() >= 3))
        n = int(m.sum())
        print(f"\n  ===== {label}: {n} clean events "
              f"({100*n/len(y):.2f}% of all, {100*n/max(int(win.sum()),1):.0f}% of its raw disagreements)")
        w_, l_ = W[m], L[m]
        print(f"    {'quantity':<26s} {'loser cluster':>14s} {'winner cluster':>15s}")
        for q, lab in (("n_tracks", "n tracks"), ("sumpt", "sum pT [GeV]"),
                       ("delta_z", "|delta z| [mm]"), ("n_inj", "n in-jet tracks"),
                       ("frac_pt_in_fwdjet", "frac pT in fwd jet"),
                       ("min_dr_to_fwdjet", "min dR to fwd jet"),
                       ("truth_purity", "truth purity"),
                       ("truth_n_hs_tracks", "n truth HS tracks")):
            lv = l_[q].abs() if q == "delta_z" else l_[q]
            wv = w_[q].abs() if q == "delta_z" else w_[q]
            print(f"    {lab:<26s} {np.nanmedian(lv):14.3f} {np.nanmedian(wv):15.3f}")
        # failure-mode taxonomy on the LOSER's pick vs the WINNER's
        bigger = (l_["sumpt"].to_numpy() > w_["sumpt"].to_numpy())
        closer = (l_["delta_z"].abs().to_numpy() < w_["delta_z"].abs().to_numpy())
        jetblind = (w_["n_inj"].to_numpy() > l_["n_inj"].to_numpy())
        nojet = (l_["n_inj"].to_numpy() == 0)
        print(f"\n    failure modes of the loser's pick (not exclusive):")
        print(f"      higher sum pT than the winner's cluster : {100*bigger.mean():5.1f}%")
        print(f"      CLOSER to the PV in z than the winner's : {100*closer.mean():5.1f}%")
        print(f"      winner's cluster has more in-jet tracks : {100*jetblind.mean():5.1f}%")
        print(f"      loser's cluster has NO in-jet tracks    : {100*nojet.mean():5.1f}%")
        print(f"      bigger AND closer in z (z+pT degenerate): {100*(bigger&closer).mean():5.1f}%")
        if wtag == "w":
            files = sorted(glob.glob(f"{a.ntuple_dir}/*.root"))
            D = pd.DataFrame({
                "fi": truth.index.get_level_values(1)[m],
                "ev": truth.index.get_level_values(2)[m],
                "big": bigger, "close": closer, "nojet": nojet,
                "t_ours": l_["t_ours"].to_numpy(), "truth": y[m],
                "n_l": l_["n_tracks"].to_numpy(), "pur_l": l_["truth_purity"].to_numpy(),
                "n_w": w_["n_tracks"].to_numpy(), "hs_w": w_["truth_n_hs_tracks"].to_numpy()})
            for cond, name in ((D.big & D.close & D.nojet, "bigger + closer in z + no jet content"),
                               (D.big & ~D.close, "bigger in pT but FURTHER in z")):
                sub = D[cond].head(2)
                print(f"\n    display candidates -- {name} ({int(cond.sum())} events):")
                for _, r in sub.iterrows():
                    print(f"      {files[int(r.fi)].split('/')[-1]}  ev {int(r.ev):5d}  "
                          f"--extra_time {r.t_ours:.1f}   (ours {int(r.n_l)} trk pur {r.pur_l:.2f}"
                          f" | WAVeS {int(r.n_w)} trk {int(r.hs_w)} HS | truth {r.truth:.1f})")


if __name__ == "__main__":
    main()
