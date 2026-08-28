#!/usr/bin/env python3
"""
guarded_transfer.py -- the guarded selection metric, trained on VBF only,
applied frozen to every other sample.

The question this answers is different from guarded_selection.py's. There the
metric was fit and evaluated on the same sample (cross-fitted); here NOTHING
from the test samples touches training. Both stages -- P(within 60 ps) and
P(switching helps) -- are fit on grid VBF canonical alone, then applied frozen
to zjets / dijet / ttbar and the three mu=0 samples. That makes the mu=0 rows a
pure robustness control: with ~no pileup TRKPTZ is nearly perfect, there is
almost nothing to gain, and a healthy guard should simply stay quiet (lost ~ 0)
rather than fire on distributions it has never seen.

Feature discipline: the feature list is the INTERSECTION of reco cluster
columns across all listed samples, so a column missing anywhere (e.g. the
vertex-fit block absent from the grid exports) is used nowhere. truth_*,
delta_t and weight are excluded as always.

Stage-1 p for the pair-building on VBF itself is 2-fold out-of-fold, so stage 2
never sees an optimistic probability; test samples are scored with the
full-VBF stage-1 fit.
"""
import argparse, glob, json, os, sys

import numpy as np
import pandas as pd
import uproot

EVT = ["sample_id", "file_idx", "event_num"]
PASS_PS = 60.0


def log(m): print(m, flush=True)


def find(input_dir, sample, tag):
    hit = glob.glob(os.path.join(input_dir, f"{sample}_{tag}_training.root"))
    return hit[0] if hit else None


def cluster_cols(path):
    return set(uproot.open(path)["clusters"].keys())


def load(path, cols):
    c = uproot.open(path)["clusters"].arrays(list(cols), library="pd")
    c["ok"] = c["delta_t"].abs() < PASS_PS
    return c


def reco_features(cols):
    drop = {"ok", "k", "p", "delta_t", "weight", "cluster_idx"}
    return sorted(x for x in cols
                  if x not in EVT and x not in drop and not x.startswith("truth_"))


def build_pairs(df, FE, K):
    g = df.groupby(EVT, sort=False)
    ib = g["trkptz_score"].idxmax()
    o = df.loc[ib].copy(); o["_eid"] = np.arange(len(o))
    emap = o.set_index(pd.MultiIndex.from_frame(o[EVT]))[["_eid"]]
    d2 = df.loc[df.index.difference(ib)].copy()
    d2["_r"] = d2.groupby(EVT, sort=False)["p"].rank(ascending=False, method="first")
    ch = d2[d2._r <= K].copy()
    ch["_eid"] = emap.reindex(pd.MultiIndex.from_frame(ch[EVT]))["_eid"].to_numpy()
    ch = ch[np.isfinite(ch._eid)]
    ch["_eid"] = ch._eid.astype(int)
    O = o.set_index("_eid")
    Xo = O.loc[ch._eid, FE].to_numpy(np.float32)
    Xn = ch[FE].to_numpy(np.float32)
    X = np.hstack([Xo, Xn, Xn - Xo,
                   O.loc[ch._eid, ["p"]].to_numpy(np.float32),
                   ch[["p"]].to_numpy(np.float32)])
    return (X, O.loc[ch._eid, "ok"].to_numpy(bool), ch["ok"].to_numpy(bool),
            ch["_eid"].to_numpy(), O)


def evaluate(df, FE, K, m2, thresholds):
    X, okO, okN, eid, O = build_pairs(df, FE, K)
    pickok = O["ok"].to_numpy(bool)
    n = len(O); trk = int(pickok.sum())
    ps = m2.predict_proba(X)[:, 1]
    order = np.lexsort((-ps, eid))
    first = np.unique(eid[order], return_index=True)[1]
    be, bp, bok = eid[order][first], ps[order][first], okN[order][first]
    out = {}
    for t in thresholds:
        ok = pickok.copy()
        sw = bp > t
        ok[be[sw]] = bok[sw]
        out[t] = dict(lost=int((pickok & ~ok).sum()), gained=int((~pickok & ok).sum()))
    return n, trk, out


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--input-dir", required=True)
    p.add_argument("--train-sample", default="vbf")
    p.add_argument("--test-samples",
                   default="zjets,dijet,ttbar,vbf_mu0,zeejets_mu0,ttbar_mu0")
    p.add_argument("--tag", default="mjj500p0")
    p.add_argument("--challengers", type=int, default=2)
    p.add_argument("--thresholds", default="0.6,0.9,0.95,0.98")
    p.add_argument("--seed", type=int, default=0)
    p.add_argument("--out", default="")
    args = p.parse_args()
    THR = [float(t) for t in args.thresholds.split(",")]

    import xgboost as xgb
    names = [args.train_sample] + args.test_samples.split(",")
    paths = {}
    for s in names:
        pth = find(args.input_dir, s, args.tag)
        if pth is None:
            log(f"  !! {s}: no {args.tag} export found, skipped")
            continue
        paths[s] = pth
    common = None
    for s, pth in paths.items():
        cols = cluster_cols(pth)
        common = cols if common is None else (common & cols)
    FE = reco_features(common)
    need = set(EVT + ["delta_t", "trkptz_score", "cluster_idx"]) | set(FE)
    log(f"{len(paths)} samples, {len(FE)} common reco features\n")

    # ---- train on VBF -------------------------------------------------------
    tr = load(paths[args.train_sample], need & common | {"delta_t"})
    rng = np.random.default_rng(args.seed)
    ev = tr[EVT].drop_duplicates().reset_index(drop=True)
    ev["k"] = rng.integers(0, 2, len(ev))
    tr = tr.merge(ev, on=EVT, how="left")
    log(f"train: {args.train_sample} {ev.shape[0]:,} events, {len(tr):,} clusters")

    def fit_m1(d):
        m = xgb.XGBClassifier(n_estimators=500, learning_rate=0.05, max_depth=6,
                              subsample=0.8, colsample_bytree=0.8, min_child_weight=20,
                              eval_metric="logloss", n_jobs=-1, tree_method="hist",
                              random_state=args.seed)
        m.fit(d[FE].to_numpy(np.float32), d.ok.to_numpy(int))
        return m

    tr["p"] = np.nan
    for kt in (0, 1):                        # out-of-fold p for pair building
        m1k = fit_m1(tr[tr.k != kt])
        msk = tr.k == kt
        tr.loc[msk, "p"] = m1k.predict_proba(tr.loc[msk, FE].to_numpy(np.float32))[:, 1]
    log("stage 1: out-of-fold p on train done")
    m1_full = fit_m1(tr)                     # scorer for the test samples

    X, okO, okN, _, _ = build_pairs(tr, FE, args.challengers)
    dec = okO != okN
    m2 = xgb.XGBClassifier(n_estimators=500, learning_rate=0.05, max_depth=5,
                           subsample=0.8, colsample_bytree=0.8, min_child_weight=10,
                           eval_metric="logloss", n_jobs=-1, tree_method="hist",
                           random_state=args.seed)
    m2.fit(X[dec], (okN & ~okO)[dec].astype(int))
    log(f"stage 2: trained on {int(dec.sum()):,} decided pairs "
        f"({int((okN & ~okO)[dec].sum()):,} switch-helps)\n")
    del tr, X, okO, okN

    # ---- apply frozen -------------------------------------------------------
    results = {}
    for s in args.test_samples.split(","):
        if s not in paths:
            continue
        te = load(paths[s], need & cluster_cols(paths[s]) | {"delta_t"})
        te["p"] = m1_full.predict_proba(te[FE].to_numpy(np.float32))[:, 1]
        n, trk, out = evaluate(te, FE, args.challengers, m2, THR)
        results[s] = dict(n=n, trkptz=trk, thr=out)
        line = f"  {s:12s} {n:>8,} ev  TRKPTZ {100*trk/n:6.2f}%"
        for t in THR:
            l, g = out[t]["lost"], out[t]["gained"]
            line += f"   [t={t}] {100*(trk-l+g)/n:6.2f}% (-{l}/+{g})"
        log(line)
        del te

    if args.out:
        json.dump({"args": vars(args), "results": results}, open(args.out, "w"),
                  indent=2, default=int)
        log(f"\nwrote {args.out}")


if __name__ == "__main__":
    main()
