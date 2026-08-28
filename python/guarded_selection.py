#!/usr/bin/env python3
"""
guarded_selection.py -- a cluster-selection metric that strictly improves on
TRKPTZ, by construction.

THE ASYMMETRIC OBJECTIVE. "Strictly improves" means two things at once: events
TRKPTZ already gets right must stay right (minimise passed-TRKPTZ & fails-new),
and events it gets wrong should be recovered (maximise failed-TRKPTZ &
passes-new). A replacement score cannot promise the first -- any free
reordering risks the 90.8% TRKPTZ already wins. So the metric is GUARDED:
TRKPTZ's cluster stands unless a challenger convinces a dedicated switch model.

TWO STAGES, BOTH CROSS-FITTED (every verdict below is out-of-fold):

  stage 1  P(cluster within 60 ps): XGBoost over the reco cluster columns.
           Its top features are trkptz_rank and trkptz_ratio_to_max (~87%),
           i.e. the model IS TRKPTZ-plus-corrections -- the corrections being
           jet association, timing self-consistency and event context.
  stage 2  P(switching helps): XGBoost on (pick, challenger) PAIRS -- both
           clusters' features, their differences, both stage-1 probabilities.
           Challengers are the top-K non-pick clusters by stage-1 p (K=2).
           Trained only on decided pairs (exactly one of the two is correct),
           which is the event-level asymmetric objective verbatim.

The stage-2 threshold is the loss/gain knob. Because deviation needs evidence,
losses fall monotonically as it rises. Cross-fitted over ALL 39,681 local VBF
events (TRKPTZ 90.83%, cluster oracle 98.35%):

    thr    pass%    lost  gained   gain:loss
    0.60   93.87%    299    1506       5.0     <- max net
    0.90   93.12%     55     965      17.5
    0.95   92.66%     18     746      41.4
    0.98   92.01%      7     476      68.0     <- near-strict

At thr=0.98 the metric loses 7 of 36,042 previously-passing events (0.019%)
while recovering 476 of 3,639 failures. The structural ceiling of a top-2
challenger set is 97.7%, so decision quality, not candidate coverage, is what
separates 93.9 from the oracle.

History: a single-stage additive guard on stage-1 p alone (delta=0.4) reached
92.84% at ~20:1 -- in git history; this two-stage version dominates it at every
operating point. A switch model trained on only ONE fold's pairs is starved
(~590 positives) and no better than the additive guard; feeding it three folds
of out-of-fold pairs is what moved the frontier. The best no-ML formula
(trkptz x (1+8*frac_pt_in_fwdjet) x exp(-0.15*(maxpull-2)+), margin 2.0)
manages 47 lost / 230 gained per ~20k events and remains the fallback if a
closed formula is required.

Truth discipline: truth_* columns, delta_t and the event weight are excluded
from every feature list; truth enters only as training labels and evaluation.
"""
import argparse, glob, json, os, sys

import numpy as np
import pandas as pd
import uproot

EVT = ["sample_id", "file_idx", "event_num"]
PASS_PS = 60.0


def log(m): print(m, flush=True)


def load(input_dir, sample):
    hit = (glob.glob(os.path.join(input_dir, f"{sample}_*training.root"))
           or glob.glob(os.path.join(input_dir, sample, f"{sample}_*training.root"))
           or glob.glob(os.path.join(input_dir, "*training.root")))
    if not hit:
        sys.exit(f"no export found for {sample} under {input_dir}")
    c = uproot.open(hit[0])["clusters"].arrays(library="pd")
    c["ok"] = c["delta_t"].abs() < PASS_PS
    return c, hit[0]


def features(c):
    drop = {"ok", "k", "p", "delta_t", "weight", "cluster_idx"}
    return [x for x in c.columns
            if x not in EVT and x not in drop
            and not x.startswith("truth_") and c[x].dtype != object]


def build_pairs(df, FE, K):
    """(pick, challenger) rows: top-K non-pick clusters by stage-1 p."""
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


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--input-dir", required=True)
    p.add_argument("--sample", default="vbf")
    p.add_argument("--threshold", type=float, default=0.95,
                   help="stage-2 switch threshold; the loss/gain knob")
    p.add_argument("--challengers", type=int, default=2)
    p.add_argument("--folds", type=int, default=4)
    p.add_argument("--seed", type=int, default=0)
    p.add_argument("--out", default="")
    args = p.parse_args()

    import xgboost as xgb
    c, path = load(args.input_dir, args.sample)
    rng = np.random.default_rng(args.seed)
    ev = c[EVT].drop_duplicates().reset_index(drop=True)
    ev["k"] = rng.integers(0, args.folds, len(ev))
    c = c.merge(ev, on=EVT, how="left")
    FE = features(c)
    log(f"{path}")
    log(f"{c[EVT].drop_duplicates().shape[0]:,} events, {len(FE)} reco features")

    # ---- stage 1: out-of-fold P(within 60 ps) for every cluster -------------
    c["p"] = np.nan
    for kt in range(args.folds):
        tr = c[c.k != kt]
        m1 = xgb.XGBClassifier(n_estimators=500, learning_rate=0.05, max_depth=6,
                               subsample=0.8, colsample_bytree=0.8, min_child_weight=20,
                               eval_metric="logloss", n_jobs=-1, tree_method="hist",
                               random_state=args.seed)
        m1.fit(tr[FE].to_numpy(np.float32), tr.ok.to_numpy(int))
        msk = c.k == kt
        c.loc[msk, "p"] = m1.predict_proba(c.loc[msk, FE].to_numpy(np.float32))[:, 1]
    log("stage 1: out-of-fold p done")

    # ---- stage 2: cross-fitted switch decision ------------------------------
    X, okO, okN, eid, O = build_pairs(c, FE, args.challengers)
    kfold = O["k"].to_numpy(int)[eid]
    pickok = O["ok"].to_numpy(bool)
    evk_all = O["k"].to_numpy(int)
    n_all = len(O)
    tot_trk = int(pickok.sum())

    lost = gained = 0
    for kt in range(args.folds):
        trm = kfold != kt
        dec = trm & (okO != okN)
        m2 = xgb.XGBClassifier(n_estimators=500, learning_rate=0.05, max_depth=5,
                               subsample=0.8, colsample_bytree=0.8, min_child_weight=10,
                               eval_metric="logloss", n_jobs=-1, tree_method="hist",
                               random_state=args.seed)
        m2.fit(X[dec], (okN & ~okO)[dec].astype(int))
        tem = ~trm
        ps = m2.predict_proba(X[tem])[:, 1]
        e_t, okN_t = eid[tem], okN[tem]
        # best challenger per event = highest switch probability
        order = np.lexsort((-ps, e_t))
        first = np.unique(e_t[order], return_index=True)[1]
        be, bp, bok = e_t[order][first], ps[order][first], okN_t[order][first]
        ok = pickok.copy()
        sw = bp > args.threshold
        ok[be[sw]] = bok[sw]
        evk = evk_all == kt
        lost += int((pickok[evk] & ~ok[evk]).sum())
        gained += int((~pickok[evk] & ok[evk]).sum())

    log(f"\nALL {n_all:,} events (every verdict out-of-fold), threshold {args.threshold}")
    log(f"  TRKPTZ  : {tot_trk:,} pass ({100*tot_trk/n_all:.2f}%)")
    log(f"  guarded : {tot_trk-lost+gained:,} pass ({100*(tot_trk-lost+gained)/n_all:.2f}%)")
    log(f"  lost {lost} ({100*lost/max(tot_trk,1):.3f}% of passing)   gained {gained}   "
        f"ratio {gained/max(lost,1):.1f}")
    if args.out:
        json.dump({"args": vars(args), "n": n_all, "trkptz": tot_trk,
                   "lost": lost, "gained": gained}, open(args.out, "w"), indent=2)
        log(f"wrote {args.out}")


if __name__ == "__main__":
    main()
