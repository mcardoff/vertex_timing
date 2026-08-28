#!/usr/bin/env python3
"""
pairwise_ranker.py — a deployed ranker over candidate pairs.

Study D established that the information is there: given the cluster Deep Sets
chose and the one it should have, a GBDT on cluster summaries separates them
84.2% of the time on Z+jets against a 50% coin flip, with a 100% control. So the
ordering signal exists and the current model is not using all of it.

This trains the model that tries to use it:

    P(C_A better than C_B)

over candidate pairs, aggregated back to a per-event ordering.

WHY THIS IS NOT STUDY D AGAIN, AND WHY IT MIGHT STILL FAIL

D is a DIAGNOSTIC on curated pairs -- exactly one right answer, handed to the
classifier. Deployment is different in a way that has already burned this study
once: a previous stage-2 re-ranker scored -0.6 against apparently similar
headroom. The difference is the prior. In deployment the ranker sees EVERY pair,
the vast majority easy, and it has to not break those while fixing the hard ones.

Two consequences for the design:

  Trained on ALL pairs within an event, not only the confusable ones. Training
  only on hard pairs would produce a model that has never seen an easy pair and
  cannot be trusted to leave one alone.

  Evaluated END-TO-END ON CORE FRACTION, never on pair accuracy. Pair accuracy is
  the quantity D already measured and it does not translate; the only number that
  counts is whether the event-level pick improves.

AGGREGATION. Each cluster's score is the sum over its pairings of P(it wins) --
a Borda count. Summing win-probabilities rather than counting wins keeps the
margin information: beating a rival 0.9 should count for more than 0.51, and with
a handful of candidates per event the difference decides ties.

ANTI-SYMMETRY. The model is fed each pair in BOTH orders and the two predictions
averaged as P(A>B) and 1-P(B>A). A pairwise model has no built-in guarantee that
swapping the inputs flips the answer, and an asymmetric model produces an
ordering that depends on how pairs were enumerated -- an artefact that would look
like a result.

The Deep Sets score is available as a feature: unlike in study D, it is not
leakage here, because pairs are enumerated exhaustively rather than constructed
around what Deep Sets got wrong. Both variants are reported, since a ranker that
only works by copying Deep Sets is not adding anything.
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

BASE_FEATS = [
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
# Present only in exports made after study E; picked up when available so this
# runs against both generations.
LOO_FEATS = ["max_abs_loo_shift", "mean_abs_loo_shift", "loo_shift_lead_pt"]


def make_pairs(sub, feats, rng, max_pairs=6):
    """All within-event pairs (capped), as (X, y, event_id, i, j).

    Capped because an event with n candidates has n(n-1)/2 pairs and the tail of
    high-multiplicity events would otherwise dominate the training set purely by
    combinatorics."""
    Xs, ys, ev, ii, jj = [], [], [], [], []
    F = sub[feats].to_numpy(np.float32)
    F = np.nan_to_num(F)
    dt = sub["abs_dt"].to_numpy()
    eid = sub["_eid"].to_numpy()
    start = 0
    for e, n in zip(*np.unique(eid, return_counts=True)):
        idx = np.arange(start, start + n); start += n
        if n < 2:
            continue
        combos = [(a, b) for k, a in enumerate(idx) for b in idx[k+1:]]
        if len(combos) > max_pairs:
            combos = [combos[k] for k in rng.choice(len(combos), max_pairs, replace=False)]
        for a, b in combos:
            Xs.append(np.concatenate([F[a], F[b], F[a] - F[b]]))
            ys.append(1 if dt[a] < dt[b] else 0)     # 1 => A is the better one
            ev.append(e); ii.append(a); jj.append(b)
    if not Xs:
        return None
    return (np.vstack(Xs).astype(np.float32), np.array(ys),
            np.array(ev), np.array(ii), np.array(jj))


def fit(Xtr, ytr):
    try:
        import xgboost as xgb
        m = xgb.XGBClassifier(n_estimators=500, learning_rate=0.06, max_depth=6,
                              subsample=0.8, colsample_bytree=0.8, min_child_weight=10,
                              eval_metric="logloss", n_jobs=-1, tree_method="hist")
    except ImportError:
        from sklearn.ensemble import HistGradientBoostingClassifier
        m = HistGradientBoostingClassifier(max_iter=350, learning_rate=0.06, random_state=0)
    m.fit(Xtr, ytr)
    return m


def borda(model, sub, feats):
    """Per-cluster Borda score: SUM over pairings of P(this one wins).

    Every pair scored in both orders and averaged, so the ordering cannot depend
    on enumeration order (see the anti-symmetry note in the docstring)."""
    F = np.nan_to_num(sub[feats].to_numpy(np.float32))
    eid = sub["_eid"].to_numpy()
    score = np.zeros(len(sub))
    start = 0
    for e, n in zip(*np.unique(eid, return_counts=True)):
        idx = np.arange(start, start + n); start += n
        if n < 2:
            continue
        A, B = [], []
        for k, a in enumerate(idx):
            for b in idx[k+1:]:
                A.append((a, b)); B.append((b, a))
        rows = np.vstack([np.concatenate([F[a], F[b], F[a] - F[b]]) for a, b in A + B])
        pr = model.predict_proba(rows.astype(np.float32))[:, 1]
        half = len(A)
        for k, (a, b) in enumerate(A):
            p = 0.5 * (pr[k] + (1.0 - pr[half + k]))   # anti-symmetrised
            score[a] += p
            score[b] += (1.0 - p)
    return score


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--input-dir", required=True)
    p.add_argument("--model", required=True, help="Deep Sets run dir, for its scores")
    p.add_argument("--selection", choices=["canonical", "loose"], default="canonical")
    p.add_argument("--seed", type=int, default=0)
    p.add_argument("--test-events", type=int, default=100_000)
    p.add_argument("--step", default="300 MB")
    p.add_argument("--max-pairs", type=int, default=6)
    p.add_argument("--out", default="")
    args = p.parse_args()

    dev = torch.device("cpu")
    files = tds.find_files(args.input_dir, args.selection)
    ck = torch.load(os.path.join(args.model, "best_model.pt"),
                    map_location=dev, weights_only=False)

    import uproot
    avail = set(uproot.open(list(files.values())[0])["clusters"].keys())
    have_loo = all(c in avail for c in LOO_FEATS)
    feats = BASE_FEATS + (LOO_FEATS if have_loo else [])
    note = ("PRESENT" if have_loo else
            "ABSENT -- this export predates study E; rerun the exporter to include them")
    log(f"features: {len(feats)}   leave-one-out columns: {note}")

    cols = sorted(set(tds.CLUSTER_COLS + feats))
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
    dsc = TST["m"][tds.KEY].copy(); dsc["deepsets"] = s
    t = df.merge(dsc, on=tds.KEY, how="inner").sort_values(tds.KEY).reset_index(drop=True)
    t["_eid"] = t.groupby(tds.EVT, sort=False).ngroup()

    prng = np.random.default_rng(7)
    out = {}
    for variant, fs in (("without DS score", feats), ("with DS score", feats + ["deepsets"])):
        log(f"\n### variant: {variant} ({len(fs)} features) ###")
        res = {}
        for sample in ("zjets", "vbf", "dijet", "ttbar"):
            sid = [k for k, v in tds.SAMPLE_NAME.items() if v == sample]
            sub = t[t["sample_id"] == sid[0]] if sid else t.iloc[:0]
            if len(sub) < 2000:
                continue
            sub = sub.sort_values("_eid").reset_index(drop=True)
            evs = sub["_eid"].unique()
            cut = evs[int(0.6 * len(evs))]
            tr = sub[sub["_eid"] < cut].reset_index(drop=True)
            te = sub[sub["_eid"] >= cut].reset_index(drop=True)
            pk = make_pairs(tr, fs, prng, args.max_pairs)
            if pk is None:
                continue
            m = fit(pk[0], pk[1])
            te = te.copy(); te["rank_score"] = borda(m, te, fs)

            def cf(col, hi=True):
                g = te.groupby(tds.EVT, sort=False)[col]
                idx = (g.idxmax() if hi else g.idxmin()).to_numpy()
                return float(100 * te.loc[idx, "within60"].mean())

            res[sample] = {"ranker": cf("rank_score"), "deepsets": cf("deepsets"),
                           "trkptz": cf("trkptz_score"), "oracle": cf("abs_dt", False),
                           "n_test_events": int(te.groupby(tds.EVT).ngroups)}
            r = res[sample]
            log(f"  {sample:8s} ranker {r['ranker']:5.1f}%   DeepSets {r['deepsets']:5.1f}%"
                f"   delta {r['ranker']-r['deepsets']:+5.1f}   (TRKPTZ {r['trkptz']:5.1f}"
                f"  oracle {r['oracle']:5.1f})")
        out[variant] = res

    log("\n" + "=" * 74)
    log("END-TO-END CORE FRACTION -- the only number that counts")
    log("=" * 74)
    for variant, res in out.items():
        log(f"{variant}:")
        for k, r in res.items():
            log(f"   {k:8s} {r['ranker']:6.1f}%  vs DeepSets {r['deepsets']:6.1f}%"
                f"   {r['ranker']-r['deepsets']:+.1f}")
    log("=" * 74)
    log("A positive delta means the pairwise ordering beats the Deep Sets ordering")
    log("on the SAME events. Anything else means study D's 84% did not survive")
    log("contact with the full candidate population -- which is exactly what the")
    log("previous stage-2 re-ranker's -0.6 was.")

    if args.out:
        os.makedirs(os.path.dirname(os.path.abspath(args.out)) or ".", exist_ok=True)
        json.dump({"args": vars(args), "have_loo": bool(have_loo), "results": out},
                  open(args.out, "w"), indent=2)
        log(f"\nwrote {args.out}")


if __name__ == "__main__":
    main()
