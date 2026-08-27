#!/usr/bin/env python3
"""
compare_models.py — do Deep Sets and a cluster-level GBDT make the SAME mistakes?

Item 7 of the agreed plan. The two models see genuinely different representations
of the same event: Deep Sets pools over a cluster's constituent TRACKS, the GBDT
reads per-cluster SUMMARY features. If their errors are partly decorrelated there
is something to gain by combining them; if they fail on the same events there is
not, and the question is closed rather than left open.

That distinction is worth measuring rather than assuming, because the two
plausible stories point opposite ways. Different representations suggest
independent errors. But the failure autopsy says misses are events where two
clusters "look the same in every reco summary, one of which happens to be HS" --
and a track-pooling model reading the same underlying tracks may well find them
equally indistinguishable. Only the measurement separates those.

WHAT IT REPORTS

  per sample, on the same test events:
    core fraction   Deep Sets / GBDT / TRKPTZ / WAVeS / oracle
    agreement       both right | DS only | GBDT only | both wrong
    phi             correlation of the two error indicators. 0 = independent
                    errors (stacking has room), 1 = identical failures (closed).
    union ceiling   fraction where AT LEAST ONE model is right. This is the
                    HARD CEILING on any combination of the two, and it is the
                    number that decides whether item 7 goes anywhere.
    stack           a concrete combiner (rank-average of the two scores), so the
                    ceiling is reported next to what is actually realised.

The union ceiling is the point of the script. A stack can only ever pick between
what the two models offer, so if the union is barely above the better model
alone, no combiner -- however clever -- can help, and no further work on this is
warranted.

USAGE
  # after a train_deepsets.py run has written best_model.pt
  python compare_models.py --input-dir <exports> --model runs/canonical \\
      --selection canonical --out runs/canonical/compare.json

The Deep Sets side is the SHIPPED model, loaded from its checkpoint -- not a
reimplementation -- and the event split is reproduced from the same --seed, so
the two models are compared on identical test events.
"""
import argparse, importlib.util, json, os, sys

import numpy as np
import pandas as pd
import torch

# Reuse the trainer's data path wholesale: same loader, same key, same split,
# same core-fraction definition. Duplicating any of it would let the two drift,
# and a comparison run on a subtly different split measures nothing.
_spec = importlib.util.spec_from_file_location(
    "tds", os.path.join(os.path.dirname(os.path.abspath(__file__)), "train_deepsets.py"))
tds = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(tds)
log = tds.log


# Cluster-level features for the GBDT. Deliberately the per-cluster SUMMARY
# view -- this is the representation whose ceiling the failure autopsy measured
# at ~67% on Z+jets, and the point here is to compare representations, not to
# hand the GBDT the tracks that make Deep Sets what it is.
#
# The two hand-designed scores are INCLUDED. They are legitimate inference-time
# inputs, and excluding them would handicap the GBDT into a strawman.
GBDT_FEATURES = [
    "n_tracks", "sumpt", "sumpt2", "maxpt", "pt_2nd", "lead_pt_frac", "meanpt", "medianpt",
    "delta_z", "delta_z_resunits", "cluster_z_sigma", "z_chi2_ndf", "max_abs_zpull",
    "z_spread", "median_z0_pull",
    "cluster_time_sigma", "time_chi2_ndf", "max_abs_tpull", "time_spread",
    "n_valid_time", "frac_valid_time", "sumpt_valid_time",
    "mean_nhgtd", "frac_nhgtd_ge2", "sumpt_w_nhgtd",
    "mean_quality", "mean_abs_eta", "eta_spread", "mean_d0_signif",
    "n_tracks_in_fwdjet", "frac_pt_in_fwdjet", "min_dr_to_fwdjet", "pt_of_nearest_fwdjet",
    "n_ghost_tracks", "sumpt_ghost",
    "trkptz_score", "waves_score",
    "n_clusters", "sumpt_rank", "sumpt_frac_of_event", "sumpt_ratio_to_max",
    "trkptz_rank", "trkptz_ratio_to_max", "is_max_sumpt", "is_max_trkptz",
    "dt_to_nearest_cluster", "dz_to_nearest_cluster", "n_clusters_within_60ps",
    "n_forward_jets", "lead_jet_pt", "sublead_jet_pt", "n_fwd_tracks_reco", "event_sumpt",
    "hgtd_time", "hgtd_valid", "hgtd_time_res", "dt_cluster_to_hgtd",
]
_leak = [f for f in GBDT_FEATURES if f.startswith("truth_") or f in ("delta_t", "abs_dt", "within60")]
assert not _leak, f"TRUTH IN GBDT FEATURES: {_leak}"


def fit_gbdt(fit, val, features):
    """Pointwise P(within60) on clusters; selection is argmax within event.

    Pointwise rather than listwise because the objective is not the variable
    under test here -- the representation is. Both models are then read the same
    way (argmax over the event), so a difference in their mistakes is a
    difference in what they can see, not in how they were trained."""
    X, y = fit[features].to_numpy(np.float32), fit["within60"].to_numpy(np.float32)
    Xv, yv = val[features].to_numpy(np.float32), val["within60"].to_numpy(np.float32)
    try:
        import xgboost as xgb
        m = xgb.XGBClassifier(
            n_estimators=2000, learning_rate=0.05, max_depth=6,
            subsample=0.8, colsample_bytree=0.8, min_child_weight=5,
            eval_metric="logloss", early_stopping_rounds=50, n_jobs=-1, tree_method="hist")
        m.fit(X, y, eval_set=[(Xv, yv)], verbose=False)
        return m, f"xgboost {xgb.__version__} ({m.best_iteration + 1} trees)"
    except ImportError:
        # Documented fallback: the AF's Jupyter kernel has had no xgboost. This
        # is a different implementation, not a different question.
        from sklearn.ensemble import HistGradientBoostingClassifier
        import sklearn
        m = HistGradientBoostingClassifier(
            max_iter=600, learning_rate=0.06, max_leaf_nodes=31,
            early_stopping=True, validation_fraction=0.15, random_state=0)
        m.fit(np.vstack([X, Xv]), np.concatenate([y, yv]))
        return m, f"sklearn HistGB {sklearn.__version__} ({m.n_iter_} iters)"


def phi_coefficient(a, b):
    """Correlation of two binary vectors. Here: the two models' ERROR indicators.

    Reported instead of raw agreement because agreement is dominated by the easy
    events both get right -- on VBF that is ~90% of them, so two models could
    agree 90% of the time while failing on completely disjoint events."""
    a, b = np.asarray(a, bool), np.asarray(b, bool)
    n11 = np.sum(a & b); n10 = np.sum(a & ~b)
    n01 = np.sum(~a & b); n00 = np.sum(~a & ~b)
    den = np.sqrt(float(n11 + n10) * (n01 + n00) * (n11 + n01) * (n10 + n00))
    return float((n11 * n00 - n10 * n01) / den) if den > 0 else float("nan")


def picks(frame, score_col):
    """Index of each event's argmax-scoring cluster, and whether it is within 60 ps."""
    idx = frame.groupby(tds.EVT, sort=False)[score_col].idxmax()
    chosen = frame.loc[idx.to_numpy()]
    return chosen.set_index(pd.MultiIndex.from_frame(chosen[tds.EVT]))["within60"].astype(bool)


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--input-dir", required=True)
    p.add_argument("--model", required=True, help="dir holding best_model.pt from train_deepsets")
    p.add_argument("--selection", choices=["canonical", "loose"], default="canonical")
    p.add_argument("--seed", type=int, default=0, help="MUST match the train_deepsets run")
    p.add_argument("--val-lo", type=float, default=0.60)
    p.add_argument("--step", default="300 MB")
    p.add_argument("--test-events", type=int, default=100_000)
    p.add_argument("--train-events", type=int, default=150_000)
    p.add_argument("--out", default="")
    args = p.parse_args()

    dev = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    files = tds.find_files(args.input_dir, args.selection)
    if not files:
        sys.exit(f"FATAL: no exports under {args.input_dir}")

    ck_path = os.path.join(args.model, "best_model.pt")
    if not os.path.exists(ck_path):
        sys.exit(f"FATAL: no best_model.pt in {args.model} -- run train_deepsets.py first")
    ck = torch.load(ck_path, map_location=dev, weights_only=False)

    # Cluster table: GBDT inputs plus everything the trainer's own loader needs.
    want = sorted(set(tds.CLUSTER_COLS + GBDT_FEATURES))
    parts = []
    for name, path in files.items():
        d = tds.read_tree(path, "clusters", want)
        log(f"  {name:12s} {len(d):>9,} clusters")
        parts.append(d)
    df = pd.concat(parts, ignore_index=True)
    df["abs_dt"] = df["delta_t"].abs()
    df["within60"] = (df["abs_dt"] < tds.PASS_PS).astype(np.float32)

    # Reproduce the trainer's split EXACTLY: same seed, same event key, same 0.7
    # boundary. Anything else compares the two models on different events.
    rng = np.random.default_rng(args.seed)
    ev = df[tds.EVT].drop_duplicates().reset_index(drop=True)
    ev["fold"] = rng.random(len(ev))
    df = df.merge(ev, on=tds.EVT, how="left")
    fit_c = df[df.fold < args.val_lo]
    val_c = df[(df.fold >= args.val_lo) & (df.fold < 0.7)]
    test_c = df[df.fold >= 0.7].copy()
    log(f"\nclusters: fit {len(fit_c):,}  val {len(val_c):,}  test {len(test_c):,}"
        f"   ({test_c.groupby(tds.EVT).ngroups:,} test events)")

    # ---- GBDT ---------------------------------------------------------------
    feats = [f for f in GBDT_FEATURES if f in df.columns]
    if len(feats) < len(GBDT_FEATURES):
        log(f"  note: {len(GBDT_FEATURES) - len(feats)} feature(s) absent from these exports")
    gb, gb_desc = fit_gbdt(fit_c, val_c, feats)
    log(f"GBDT: {gb_desc}, {len(feats)} features")
    test_c["gbdt"] = gb.predict_proba(test_c[feats].to_numpy(np.float32))[:, 1]

    # ---- Deep Sets: the shipped checkpoint, not a reimplementation -----------
    fold = ev.set_index(tds.EVT)["fold"]
    per_train = {s: 0 for s in files}          # eval only; no fitting here
    _, TE_T = tds.stream_tracks(files, fold, per_train,
                                args.test_events // max(1, len(files)), args.step)
    # The checkpoint's "key" is the (pool, lambda) pair the run selected on
    # validation, so the pool is recovered from it rather than assumed -- loading
    # a gated model's weights into a sum-pooled net would score silently wrong.
    pool = ck["key"][0] if isinstance(ck.get("key"), (list, tuple)) else "sum"
    net = tds.DeepSets(len(ck["features"]) + len(ck["na_cols"]),
                       ck["hidden"], ck["embed"], pool=pool)
    net.load_state_dict(ck["state_dict"])
    net = net.to(dev).eval()
    log(f"Deep Sets: pool={pool}, {ck['hidden']}/{ck['embed']}, "
        f"{len(ck['features'])} features + {len(ck['na_cols'])} indicators")

    labels = df[tds.KEY + ["within60"]].drop_duplicates(tds.KEY)
    TST = tds.make_tensors(TE_T, pd.Series(ck["mu"]), pd.Series(ck["sd"]),
                           ck["na_cols"], labels, dev)
    with torch.no_grad():
        scores = net(TST["X"], TST["c"], TST["ncl"]).cpu().numpy()
    # make_tensors returns `m`: one row per cluster, in the same order as the
    # scores, carrying KEY. That is the join back onto the cluster table.
    ds = TST["m"][tds.KEY].copy()
    ds["deepsets"] = scores
    test_c = test_c.merge(ds, on=tds.KEY, how="inner")
    log(f"scored {test_c.groupby(tds.EVT).ngroups:,} events with both models")

    # ---- compare ------------------------------------------------------------
    rows, out = [], {}
    for sid, sub in list(test_c.groupby("sample_id")) + [(None, test_c)]:
        nm = "ALL" if sid is None else tds.SAMPLE_NAME.get(sid, str(sid))
        d_ok, g_ok = picks(sub, "deepsets"), picks(sub, "gbdt")
        d_ok, g_ok = d_ok.align(g_ok, join="inner")
        n = len(d_ok)
        if n == 0:
            continue
        both = float((d_ok & g_ok).mean()); dsonly = float((d_ok & ~g_ok).mean())
        gbonly = float((~d_ok & g_ok).mean()); neither = float((~d_ok & ~g_ok).mean())
        union = both + dsonly + gbonly
        # A concrete combiner, so the ceiling is quoted beside what is realised.
        sub = sub.copy()
        for c in ("deepsets", "gbdt"):
            sub[f"r_{c}"] = sub.groupby(tds.EVT, sort=False)[c].rank(pct=True)
        sub["stack"] = sub["r_deepsets"] + sub["r_gbdt"]
        st = float(picks(sub, "stack").mean())
        rec = {"n_events": int(n),
               "deepsets": 100 * float(d_ok.mean()), "gbdt": 100 * float(g_ok.mean()),
               "trkptz": tds.core_fraction(sub, "trkptz_score"),
               "waves": tds.core_fraction(sub, "waves_score"),
               "oracle": tds.core_fraction(sub, "abs_dt", False),
               "both_right": 100 * both, "ds_only": 100 * dsonly,
               "gbdt_only": 100 * gbonly, "both_wrong": 100 * neither,
               "union_ceiling": 100 * union, "stack": 100 * st,
               "phi_errors": phi_coefficient(~d_ok, ~g_ok)}
        out[nm] = rec
        rows.append((nm, rec))

    log("\n" + "=" * 108)
    log(f"{'sample':12s}{'n':>8s}{'TRKPTZ':>8s}{'WAVeS':>8s}{'DeepSets':>10s}"
        f"{'GBDT':>8s}{'stack':>8s}{'union':>8s}{'oracle':>8s}{'phi(err)':>10s}")
    log("-" * 108)
    for nm, r in rows:
        log(f"{nm:12s}{r['n_events']:>8,}{r['trkptz']:>7.1f}%{r['waves']:>7.1f}%"
            f"{r['deepsets']:>9.1f}%{r['gbdt']:>7.1f}%{r['stack']:>7.1f}%"
            f"{r['union_ceiling']:>7.1f}%{r['oracle']:>7.1f}%{r['phi_errors']:>10.3f}")
    log("=" * 108)
    log(f"\n{'sample':12s}{'both right':>12s}{'DS only':>10s}{'GBDT only':>11s}{'both wrong':>12s}")
    for nm, r in rows:
        log(f"{nm:12s}{r['both_right']:>11.1f}%{r['ds_only']:>9.1f}%"
            f"{r['gbdt_only']:>10.1f}%{r['both_wrong']:>11.1f}%")

    log("\nreading this:")
    log("  union - max(DeepSets, GBDT) is the MOST any combiner could add. If that")
    log("  gap is small, item 7 is closed regardless of how the stack performs.")
    log("  phi near 1 means the two fail on the same events, i.e. the per-cluster")
    log("  and per-track views are hitting the same wall rather than different ones.")

    if args.out:
        os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)
        with open(args.out, "w") as fh:
            json.dump({"args": vars(args), "gbdt": gb_desc,
                       "features": feats, "per_sample": out}, fh, indent=2)
        log(f"\nwrote {args.out}")


if __name__ == "__main__":
    main()
