#!/usr/bin/env python3
"""
pairwise_probe.py — study D. Is the information there at all?

THE QUESTION. Take the Z+jets events where a good candidate exists and Deep Sets
picked a different one. Hand a classifier just those two clusters and ask which
has the smaller |dt|. If a strong classifier cannot beat a coin flip, then the
reconstructed inputs do not contain what is needed to tell them apart, and NO
selector -- a better Deep Sets, a top-3 re-ranker, a pairwise ranker -- can close
that part of the gap. If it can, there is exploitable signal and the ranking
architecture is worth revisiting.

This is the experiment that turns "we tried things and they did not help" into a
measured statement about the information available.

THE ONE THING THAT WOULD INVALIDATE IT. The pairs are constructed as
  A = the cluster Deep Sets scored highest (and got wrong)
  B = the cluster with the smallest |dt| (the right answer)
so by construction ds_score(A) > ds_score(B), ALWAYS. Feeding the Deep Sets score
in would let the probe score ~100% with the rule "pick whichever has the LOWER
Deep Sets score" -- perfectly predictive of the label and completely useless,
since applying it to the full population would mean always overturning the model.
That is leakage of the construction, not information about the physics, so the
Deep Sets score is excluded. TRKPTZ and WAVeS stay: they take no part in
selecting the pair, so they carry no such tell.

CONTROLS. A negative result is only meaningful if the probe can learn anything at
all, so two controls run on the identical pipeline:
  n_tracks   predict which cluster has more tracks -- trivially in the features,
             so this must come out near 100% or the probe itself is broken.
  any-pair   the same |dt| question on RANDOM cluster pairs from the same events
             rather than the hard ones. This is the "how much signal exists in
             general" reference; the gap between it and the main task is the part
             that is specifically hard.
"""
import argparse, importlib.util, json, os, sys

import numpy as np
import pandas as pd
import torch

_spec = importlib.util.spec_from_file_location(
    "tds", os.path.join(os.path.dirname(os.path.abspath(__file__)), "train_deepsets.py"))
tds = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(tds)
log = tds.log

# Cluster-level features. Reco-only; no truth, and deliberately no Deep Sets
# score (see the docstring). trkptz/waves stay -- they do not select the pair.
FEATS = [
    "n_tracks", "sumpt", "sumpt2", "maxpt", "pt_2nd", "lead_pt_frac", "meanpt", "medianpt",
    "delta_z", "delta_z_resunits", "cluster_z_sigma", "z_chi2_ndf", "max_abs_zpull",
    "z_spread", "median_z0_pull",
    "cluster_time", "cluster_time_sigma", "time_chi2_ndf", "max_abs_tpull", "time_spread",
    "n_valid_time", "frac_valid_time", "sumpt_valid_time", "n_tpull_gt2",
    "min_timeRes", "median_timeRes", "max_timeRes",
    "mean_nhgtd", "frac_nhgtd_ge2", "sumpt_w_nhgtd",
    "mean_quality", "mean_abs_eta", "eta_spread", "mean_d0_signif",
    "n_tracks_in_fwdjet", "frac_pt_in_fwdjet", "min_dr_to_fwdjet", "pt_of_nearest_fwdjet",
    "n_ghost_tracks", "sumpt_ghost", "trkptz_score", "waves_score",
    "dt_to_nearest_cluster", "dz_to_nearest_cluster", "n_clusters",
    "hgtd_time", "hgtd_valid", "hgtd_time_res", "dt_cluster_to_hgtd",
]
_leak = [f for f in FEATS if f.startswith("truth_") or f in ("delta_t", "abs_dt",
                                                             "within60", "deepsets")]
assert not _leak, f"LEAKAGE IN PROBE FEATURES: {_leak}"


def fit_probe(Xtr, ytr, Xte, yte):
    """Return test accuracy. xgboost when available, sklearn otherwise."""
    try:
        import xgboost as xgb
        m = xgb.XGBClassifier(n_estimators=400, learning_rate=0.06, max_depth=6,
                              subsample=0.8, colsample_bytree=0.8, min_child_weight=5,
                              eval_metric="logloss", n_jobs=-1, tree_method="hist")
    except ImportError:
        from sklearn.ensemble import HistGradientBoostingClassifier
        m = HistGradientBoostingClassifier(max_iter=300, learning_rate=0.06, random_state=0)
    m.fit(Xtr, ytr)
    return float((m.predict(Xte) == yte).mean())


def build_pairs(pairs, rng, label_col):
    """Randomise which of the two clusters is presented first.

    Without this the label would be perfectly predicted by position rather than
    by physics, and the probe would report 100% while learning nothing."""
    first_is_better = rng.random(len(pairs)) < 0.5
    P = np.where(first_is_better[:, None], pairs["A"].to_numpy(), pairs["B"].to_numpy())
    Q = np.where(first_is_better[:, None], pairs["B"].to_numpy(), pairs["A"].to_numpy())
    X = np.hstack([P, Q, P - Q]).astype(np.float32)
    y = first_is_better.astype(int) if label_col == "B_better" else None
    return X, y


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--input-dir", required=True)
    p.add_argument("--model", required=True)
    p.add_argument("--selection", choices=["canonical", "loose"], default="canonical")
    p.add_argument("--seed", type=int, default=0)
    p.add_argument("--test-events", type=int, default=100_000)
    p.add_argument("--step", default="300 MB")
    p.add_argument("--out", default="")
    args = p.parse_args()

    dev = torch.device("cpu")
    files = tds.find_files(args.input_dir, args.selection)
    ck = torch.load(os.path.join(args.model, "best_model.pt"),
                    map_location=dev, weights_only=False)

    cols = sorted(set(tds.CLUSTER_COLS + FEATS))
    df = pd.concat([tds.read_tree(pth, "clusters", cols) for pth in files.values()],
                   ignore_index=True)
    df["abs_dt"] = df["delta_t"].abs()
    df["within60"] = (df["abs_dt"] < tds.PASS_PS).astype(np.float32)

    rng = np.random.default_rng(args.seed)
    ev = df[tds.EVT].drop_duplicates().reset_index(drop=True)
    ev["fold"] = rng.random(len(ev))
    fold = ev.set_index(tds.EVT)["fold"]
    df = df.merge(ev, on=tds.EVT, how="left")

    _, TE_T = tds.stream_tracks(files, fold, {s: 0 for s in files},
                                args.test_events // max(1, len(files)), args.step)
    labels = df[tds.KEY + ["within60"]].drop_duplicates(tds.KEY)
    TST = tds.make_tensors(TE_T, pd.Series(ck["mu"]), pd.Series(ck["sd"]),
                           ck["na_cols"], labels, dev)
    pool = ck["key"][0] if isinstance(ck.get("key"), (list, tuple)) else "sum"
    net = tds.DeepSets(len(ck["features"]) + len(ck["na_cols"]),
                       ck["hidden"], ck["embed"], pool=pool)
    net.load_state_dict(ck["state_dict"]); net.eval()
    with torch.no_grad():
        s = net(TST["X"], TST["c"], TST["ncl"]).cpu().numpy()
    ds = TST["m"][tds.KEY].copy(); ds["deepsets"] = s
    t = df.merge(ds, on=tds.KEY, how="inner")

    results = {}
    prng = np.random.default_rng(12345)
    for sample in ("zjets", "vbf", "dijet", "ttbar"):
        sid = [k for k, v in tds.SAMPLE_NAME.items() if v == sample]
        sub = t[t["sample_id"] == sid[0]] if sid else t.iloc[:0]
        if not len(sub):
            continue

        g = sub.groupby(tds.EVT, sort=False)
        pick = sub.loc[g["deepsets"].idxmax().to_numpy()]          # A: what DS chose
        best = sub.loc[g["abs_dt"].idxmin().to_numpy()]            # B: the right answer
        m = pick[tds.EVT + ["within60"]].merge(
            best[tds.EVT + ["within60"]], on=tds.EVT, suffixes=("_pick", "_best"))
        # DS failures on RECOVERABLE events: a good candidate existed, DS missed it.
        fail = m[(m["within60_best"] > 0) & (m["within60_pick"] < 1)][tds.EVT]
        if len(fail) < 400:
            log(f"  {sample}: only {len(fail)} failure pairs, skipping")
            continue

        A = pick.merge(fail, on=tds.EVT)[FEATS].to_numpy(np.float32)
        B = best.merge(fail, on=tds.EVT)[FEATS].to_numpy(np.float32)
        A = np.nan_to_num(A); B = np.nan_to_num(B)

        n = len(A); cut = int(0.6 * n)
        flip = prng.random(n) < 0.5
        P = np.where(flip[:, None], A, B); Q = np.where(flip[:, None], B, A)
        X = np.hstack([P, Q, P - Q]); y = (~flip).astype(int)   # 1 => P is the better one
        acc = fit_probe(X[:cut], y[:cut], X[cut:], y[cut:])

        # control 1: a question that IS in the features
        ntr_i = FEATS.index("n_tracks")
        y2 = (P[:, ntr_i] > Q[:, ntr_i]).astype(int)
        acc_ctrl = fit_probe(X[:cut], y2[:cut], X[cut:], y2[cut:])

        # control 2: the same |dt| question on RANDOM pairs, not the hard ones
        two = sub.groupby(tds.EVT).filter(lambda d: len(d) >= 2)
        samp = two.groupby(tds.EVT, sort=False).sample(n=2, random_state=7)
        gg = samp.groupby(tds.EVT, sort=False)
        R1 = gg.nth(0); R2 = gg.nth(1)
        k = min(len(R1), len(R2), 40000)
        F1 = np.nan_to_num(R1[FEATS].to_numpy(np.float32)[:k])
        F2 = np.nan_to_num(R2[FEATS].to_numpy(np.float32)[:k])
        better = (R1["abs_dt"].to_numpy()[:k] < R2["abs_dt"].to_numpy()[:k]).astype(int)
        Xr = np.hstack([F1, F2, F1 - F2]); cr = int(0.6 * k)
        acc_any = fit_probe(Xr[:cr], better[:cr], Xr[cr:], better[cr:])

        results[sample] = {"n_failure_pairs": int(n), "acc_failure": 100 * acc,
                           "acc_control_ntracks": 100 * acc_ctrl,
                           "acc_random_pairs": 100 * acc_any}
        log(f"  {sample:8s} pairs {n:>6,}   HARD {100*acc:5.1f}%   "
            f"random-pair {100*acc_any:5.1f}%   control(n_tracks) {100*acc_ctrl:5.1f}%")

    log("\n" + "=" * 76)
    log(f"{'sample':10s}{'pairs':>9s}{'HARD pairs':>13s}{'random pairs':>14s}{'control':>11s}")
    log("-" * 76)
    for k, v in results.items():
        log(f"{k:10s}{v['n_failure_pairs']:>9,}{v['acc_failure']:>12.1f}%"
            f"{v['acc_random_pairs']:>13.1f}%{v['acc_control_ntracks']:>10.1f}%")
    log("=" * 76)
    log("HARD pairs at ~50% with the control near 100% is the information-limit")
    log("verdict: the probe works, and the discriminating information is absent.")
    log("Well above 50% means signal exists and the ranking stage is worth revisiting.")

    if args.out:
        os.makedirs(os.path.dirname(os.path.abspath(args.out)) or ".", exist_ok=True)
        json.dump({"args": vars(args), "results": results}, open(args.out, "w"), indent=2)
        log(f"\nwrote {args.out}")


if __name__ == "__main__":
    main()
