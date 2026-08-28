#!/usr/bin/env python3
"""
aggregate_t0.py — estimate the vertex t0 directly, instead of selecting a cluster.

Every model in this study so far picks ONE cluster and reports its time. Nothing
in the +-60 ps criterion requires that; it only requires a number. Measured
ceilings on canonical Z+jets say the distinction is worth a lot:

    TRKPTZ, selects a cluster                     61.9%
    Deep Sets, selects a cluster                  67.3%
    realistic selection ceiling                  ~74.3%
    median over PERFECTLY tagged HS tracks        81.7%
    truncated mean, same tracks                   82.3%

The gap exists for a specific reason. 13.7% of true hard-scatter tracks sit more
than 200 ps from truth while their quoted timeRes is ~25 ps, and
Cluster::calculateTime uses an inverse-variance mean -- the estimator MOST
sensitive to that tail, because it trusts sigma. On identical tracks the median
scores 81.7% against the inverse-variance mean's 65.3%.

Swapping the estimator inside a cluster does NOT recover this (61.9 -> 61.6):
clusterTracksInTime already groups by time agreement, so within one cluster the
times concur and median == mean. The mis-timed tracks are split into their OWN
clusters, which is exactly why no single-cluster selection can reach them and
why aggregating across the event can.

So: tag tracks, then aggregate robustly.

  stage 1   per-track P(hard scatter). ~1.2M labelled tracks per sample, a far
            better posed problem than ranking ~5.7 candidates per event.
  stage 2   robust weighted estimate of t0 over the tagged tracks.

A blind median over the raw associated list is NOT a baseline worth beating on
its own terms -- that list is only 20.3% hard scatter on Z+jets and scores 44.9%
-- but it is reported so the tagger's contribution is visible.

TRUTH USE. truth_is_hs is the stage-1 training label and nothing else. It is
never a model input, and the reported core fraction is always computed against
the true HS time on held-out events.
"""
import argparse, glob, os, sys

import numpy as np
import pandas as pd
import uproot

EVT = ["sample_id", "file_idx", "event_num"]
KEY = EVT + ["cluster_idx"]
PASS_PS = 60.0
SAMPLE_NAME = {-1.0: "local", 0.0: "vbf", 1.0: "zjets", 2.0: "dijet", 3.0: "zeejets",
               4.0: "vbf_mu0", 5.0: "zeejets_mu0", 6.0: "ttbar_mu0", 7.0: "ttbar"}

# Per-track inputs. Deliberately excludes every truth_* column; t_pull_cluster
# and loo_* are functions of the cluster's own reconstructed time, not of truth.
FEATS = ["pt", "eta", "z0", "d0", "sigma_z0", "sigma_d0", "time", "timeRes",
         "time_valid", "quality", "nhgtd_hits", "z0_pull_pv", "t_pull_cluster",
         "loo_shift", "loo_pull", "dr_nearest_fwdjet", "pt_nearest_fwdjet",
         "is_ghost_of_nearest", "is_lepton", "cluster_time", "cluster_delta_z",
         "dr_nearest_anyjet", "pt_nearest_anyjet", "is_in_any_jet",
         "dz_to_nearest_pu_vtx_trk", "closer_to_pu_than_pv",
         "cluster_dz_to_nearest_pu_vtx", "cluster_pv_pu_dz_ratio",
         "cluster_nearest_vtx_waves_rank", "cluster_nearest_vtx_waves_frac",
         "cluster_closest_vtx_is_pv"]


def log(m): print(m, flush=True)


def load(path, max_events):
    """Tracks + the event's true HS time, for one sample."""
    c = uproot.open(path)["clusters"].arrays(
        KEY + ["delta_t", "cluster_time", "trkptz_score", "waves_score"], library="pd")
    c["t_truth"] = c["cluster_time"] - c["delta_t"]
    truth = c.groupby(EVT, sort=False)["t_truth"].first().rename("t_truth")
    if max_events and len(truth) > max_events:
        keep = truth.index[:max_events]
        truth = truth.loc[keep]
        c = c.set_index(EVT).loc[keep].reset_index()
    have = set(uproot.open(path)["tracks"].keys())
    cols = [f for f in FEATS if f in have]
    miss = [f for f in FEATS if f not in have]
    t = uproot.open(path)["tracks"].arrays(
        KEY + ["track_idx", "truth_is_hs"] + cols, library="pd")
    t = t.merge(truth.reset_index(), on=EVT, how="inner")
    return c, t, cols, miss


def agg_rules(tk, truth, n_ev, pcol=None, thresholds=(0.3, 0.5, 0.7)):
    """t0 per event under several selection + estimator combinations."""
    out = {}
    tk = tk[(tk.time_valid > 0) & (tk.timeRes > 0)]

    def score(series, lab):
        s = series.dropna()
        ok = (s - truth.reindex(s.index)).abs() < PASS_PS
        out[lab] = 100.0 * ok.sum() / n_ev

    def trunc(d):
        m = d.time.median()
        k = d[(d.time - m).abs() < 100]
        return k.time.mean() if len(k) else m

    score(tk.groupby(EVT, sort=False)["time"].median(), "median, ALL tracks (no tagger)")
    if pcol is not None:
        for th in thresholds:
            s = tk[tk[pcol] > th]
            if not len(s):
                continue
            g = s.groupby(EVT, sort=False)
            score(g["time"].median(), f"median, p>{th}")
            score(g.apply(trunc, include_groups=False), f"truncated mean, p>{th}")
            score(g.apply(lambda d: (d.time / d.timeRes**2).sum()
                                    / (1 / d.timeRes**2).sum(), include_groups=False),
                  f"inv-var mean, p>{th}")
        # weighted median using the tagger probability as the weight
        def wmed(d):
            o = d.sort_values("time")
            cw = o[pcol].cumsum()
            tot = o[pcol].sum()
            if tot <= 0:
                return np.nan
            return o.time.to_numpy()[np.searchsorted(cw.to_numpy(), 0.5 * tot)]
        score(tk.groupby(EVT, sort=False).apply(wmed, include_groups=False),
              "prob-weighted median (all tracks)")
    hs = tk[tk.truth_is_hs > 0.5]
    g = hs.groupby(EVT, sort=False)
    score(g["time"].median(), "CEILING: median, perfect tagger")
    score(g.apply(trunc, include_groups=False), "CEILING: truncated mean, perfect")
    return out


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--input-dir", required=True)
    p.add_argument("--samples", default="vbf,zjets,dijet,ttbar")
    p.add_argument("--max-events", type=int, default=60000,
                   help="per sample, before the train/test split")
    p.add_argument("--seed", type=int, default=0)
    p.add_argument("--test-frac", type=float, default=0.4)
    args = p.parse_args()

    C, T, cols, miss = {}, {}, None, None
    for s in args.samples.split(","):
        hit = (glob.glob(os.path.join(args.input_dir, f"{s}_*training.root"))
               or glob.glob(os.path.join(args.input_dir, s, f"{s}_*training.root")))
        if not hit:
            log(f"  {s:8s} NOT FOUND"); continue
        c, t, cols, miss = load(hit[0], args.max_events)
        C[s], T[s] = c, t
        log(f"  {s:8s} {c[EVT].drop_duplicates().shape[0]:>7,} ev  "
            f"{len(t):>9,} track rows  HS {100*t.truth_is_hs.mean():.1f}%")
    if not C:
        sys.exit("no samples loaded")
    if miss:
        log(f"\n  !! {len(miss)} feature(s) absent from this export, dropped: {miss}")

    # Event-level split, shared across samples via the same RNG stream.
    rng = np.random.default_rng(args.seed)
    for s in T:
        ev = T[s][EVT].drop_duplicates().reset_index(drop=True)
        ev["fold"] = rng.random(len(ev))
        T[s] = T[s].merge(ev, on=EVT, how="left")
        C[s] = C[s].merge(ev, on=EVT, how="left")

    tr = pd.concat([T[s][T[s].fold >= args.test_frac] for s in T], ignore_index=True)
    log(f"\nstage 1: per-track HS tagger on {len(tr):,} training tracks "
        f"({100*tr.truth_is_hs.mean():.1f}% positive)")
    try:
        import xgboost as xgb
        m = xgb.XGBClassifier(n_estimators=400, learning_rate=0.06, max_depth=7,
                              subsample=0.8, colsample_bytree=0.8, min_child_weight=10,
                              eval_metric="logloss", n_jobs=-1, tree_method="hist")
    except ImportError:
        from sklearn.ensemble import HistGradientBoostingClassifier
        m = HistGradientBoostingClassifier(max_iter=400, learning_rate=0.06, random_state=0)
    m.fit(tr[cols].to_numpy(np.float32), tr.truth_is_hs.to_numpy(int))

    log("\n" + "=" * 78)
    log("AGGREGATE t0 -- core fraction on HELD-OUT events")
    log("=" * 78)
    rows, names = {}, []
    for s in sorted(T):
        te = T[s][T[s].fold < args.test_frac].copy()
        ce = C[s][C[s].fold < args.test_frac]
        truth = ce.groupby(EVT, sort=False)["t_truth"].first()
        n_ev = len(truth)
        if n_ev < 200:
            continue
        te["p"] = m.predict_proba(te[cols].to_numpy(np.float32))[:, 1]
        te = te.drop_duplicates(EVT + ["track_idx"])
        r = agg_rules(te, truth, n_ev, pcol="p")
        ce = ce.copy()
        ce["ok"] = ce["delta_t"].abs() < PASS_PS
        for lab, col in (("TRKPTZ (selects a cluster)", "trkptz_score"),
                         ("WAVeS (selects a cluster)", "waves_score")):
            pick = ce.loc[ce.groupby(EVT, sort=False)[col].idxmax()]
            r[lab] = 100.0 * pick["ok"].sum() / n_ev
        r["oracle over clusters"] = 100.0 * ce.groupby(EVT, sort=False)["ok"].max().sum() / n_ev
        rows[s] = r
        names.append(s)
        # tagger quality, for the record
        auc = ""
        try:
            from sklearn.metrics import roc_auc_score
            auc = f"   tagger AUC {roc_auc_score(te.truth_is_hs, te.p):.3f}"
        except Exception:
            pass
        log(f"\n  {s}  ({n_ev:,} test events){auc}")
    order = ["TRKPTZ (selects a cluster)", "WAVeS (selects a cluster)",
             "median, ALL tracks (no tagger)", "prob-weighted median (all tracks)"]
    order += [k for k in rows[names[0]] if k.startswith(("median, p", "truncated", "inv-var"))]
    order += ["CEILING: median, perfect tagger", "CEILING: truncated mean, perfect",
              "oracle over clusters"]
    log("\n" + f"  {'estimator':38s}" + "".join(f"{n:>10s}" for n in names))
    log("  " + "-" * (38 + 10 * len(names)))
    for k in order:
        if k not in rows[names[0]]:
            continue
        log(f"  {k:38s}" + "".join(f"{rows[n].get(k, float('nan')):>9.1f}%" for n in names))


if __name__ == "__main__":
    main()
