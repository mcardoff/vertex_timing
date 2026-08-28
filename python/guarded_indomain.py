#!/usr/bin/env python3
"""
guarded_indomain.py -- retrain the SWITCH stage in-domain, keeping stage 1
VBF-trained; plus the fully in-domain variant as the reference point.

Three variants of the guarded metric, forming an ablation of where the domain
knowledge lives:

  A  both stages VBF-frozen        (measured by guarded_transfer.py)
  B  stage-1 VBF, stage-2 in-domain   <- this run's request
  C  both stages in-domain            (guarded_selection.py recipe)

B answers: is the transfer gap in the probability model or in the switch
decision? Stage 2 is cheap to retrain per sample (it sees only ~1 pair per
event), so if B closes most of the A->C gap, a per-sample switch on a shared
probability model is the deployable shape.

All in-domain training is 4-fold cross-fitted; every quoted verdict is
out-of-fold. Stage-2 training pairs are built from honest probabilities:
variant B's p comes from the VBF model (out-of-domain, hence honest), variant
C's from in-domain out-of-fold fits.
"""
import argparse, glob, json, os, sys

import numpy as np
import pandas as pd
import uproot

EVT = ["sample_id", "file_idx", "event_num"]
PASS_PS = 60.0


def log(m): print(m, flush=True)


def find(d, s, tag):
    h = glob.glob(os.path.join(d, f"{s}_{tag}_training.root"))
    return h[0] if h else sys.exit(f"missing {s}")


def load(path):
    c = uproot.open(path)["clusters"].arrays(library="pd")
    c["ok"] = c["delta_t"].abs() < PASS_PS
    return c


def reco_features(cols):
    drop = {"ok", "k", "p", "delta_t", "weight", "cluster_idx"}
    return sorted(x for x in cols
                  if x not in EVT and x not in drop and not x.startswith("truth_"))


def xgb_clf(depth, mcw, seed):
    import xgboost as xgb
    return xgb.XGBClassifier(n_estimators=500, learning_rate=0.05, max_depth=depth,
                             subsample=0.8, colsample_bytree=0.8, min_child_weight=mcw,
                             eval_metric="logloss", n_jobs=-1, tree_method="hist",
                             random_state=seed)


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


def crossfit_switch(df, FE, K, THR, seed):
    """m2 4-fold on this sample's pairs; df must carry honest p and fold k."""
    X, okO, okN, eid, O = build_pairs(df, FE, K)
    kf = O["k"].to_numpy(int)[eid]
    pickok = O["ok"].to_numpy(bool)
    evk_all = O["k"].to_numpy(int)
    n, trk = len(O), int(pickok.sum())
    tot = {t: [0, 0] for t in THR}
    for kt in range(4):
        trm = kf != kt
        dec = trm & (okO != okN)
        m2 = xgb_clf(5, 10, seed)
        m2.fit(X[dec], (okN & ~okO)[dec].astype(int))
        tem = ~trm
        ps = m2.predict_proba(X[tem])[:, 1]
        e_t, okN_t = eid[tem], okN[tem]
        order = np.lexsort((-ps, e_t))
        first = np.unique(e_t[order], return_index=True)[1]
        be, bp, bok = e_t[order][first], ps[order][first], okN_t[order][first]
        evk = evk_all == kt
        for t in THR:
            ok = pickok.copy()
            sw = bp > t
            ok[be[sw]] = bok[sw]
            tot[t][0] += int((pickok[evk] & ~ok[evk]).sum())
            tot[t][1] += int((~pickok[evk] & ok[evk]).sum())
    return n, trk, tot


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--input-dir", required=True)
    p.add_argument("--samples", default="zjets,ttbar")
    p.add_argument("--tag", default="mjj500p0")
    p.add_argument("--challengers", type=int, default=2)
    p.add_argument("--seed", type=int, default=0)
    p.add_argument("--out", default="")
    args = p.parse_args()
    THR = (0.6, 0.9, 0.95, 0.98)

    # ---- stage-1 VBF model (shared across variant B) -------------------------
    vbf = load(find(args.input_dir, "vbf", args.tag))
    FE = reco_features(set(vbf.columns))
    log(f"{len(FE)} reco features (grid schema)")
    m1_vbf = xgb_clf(6, 20, args.seed)
    m1_vbf.fit(vbf[FE].to_numpy(np.float32), vbf.ok.to_numpy(int))
    log(f"stage-1 VBF model trained on {len(vbf):,} clusters")
    del vbf

    rng = np.random.default_rng(args.seed)
    results = {}
    for s in args.samples.split(","):
        c = load(find(args.input_dir, s, args.tag))
        ev = c[EVT].drop_duplicates().reset_index(drop=True)
        ev["k"] = rng.integers(0, 4, len(ev))
        c = c.merge(ev, on=EVT, how="left")
        log(f"\n=== {s}: {len(ev):,} events ===")

        # variant B: p from the VBF model
        c["p"] = m1_vbf.predict_proba(c[FE].to_numpy(np.float32))[:, 1]
        nB, trkB, totB = crossfit_switch(c, FE, args.challengers, THR, args.seed)
        log(f"  B (stage-1 VBF, switch in-domain):")
        for t in THR:
            l, g = totB[t]
            log(f"    t={t:4}  {100*(trkB-l+g)/nB:6.2f}%  (-{l}/+{g})")

        # variant C: fully in-domain (out-of-fold stage-1 too)
        c["p"] = np.nan
        for kt in range(4):
            m1k = xgb_clf(6, 20, args.seed)
            tr = c[c.k != kt]
            m1k.fit(tr[FE].to_numpy(np.float32), tr.ok.to_numpy(int))
            msk = c.k == kt
            c.loc[msk, "p"] = m1k.predict_proba(c.loc[msk, FE].to_numpy(np.float32))[:, 1]
        nC, trkC, totC = crossfit_switch(c, FE, args.challengers, THR, args.seed)
        log(f"  C (fully in-domain):")
        for t in THR:
            l, g = totC[t]
            log(f"    t={t:4}  {100*(trkC-l+g)/nC:6.2f}%  (-{l}/+{g})")
        results[s] = dict(n=nB, trkptz=trkB,
                          B={t: totB[t] for t in THR}, C={t: totC[t] for t in THR})
        del c

    if args.out:
        json.dump({"args": vars(args), "results": results},
                  open(args.out, "w"), indent=2, default=int)
        log(f"\nwrote {args.out}")


if __name__ == "__main__":
    main()
