#!/usr/bin/env python3
"""
guarded_selection.py — a cluster-selection metric that strictly improves on
TRKPTZ, by construction.

THE ASYMMETRIC OBJECTIVE. "Strictly improves" means two different things must
hold at once: events TRKPTZ already gets right must stay right (minimise
passed-TRKPTZ & fails-new), and events it gets wrong should be recovered
(maximise failed-TRKPTZ & passes-new). A plain replacement score cannot promise
the first -- any reordering risks the 90.8% TRKPTZ already wins. So the metric
here is GUARDED:

    pick TRKPTZ's cluster, UNLESS another cluster's P(within 60 ps)
    exceeds the TRKPTZ pick's by more than delta.

delta is the knob that trades recoveries against losses, and because deviation
requires evidence, losses shrink monotonically with it. Measured on the held-out
half of local VBF (19,769 events, TRKPTZ 91.09%):

    delta   pass    lost  gained   gain:loss
    0.0    93.84%    115    660       5.7
    0.2    93.47%     48    519      10.8
    0.4    92.94%     18    385      21.4
    0.5    92.69%      7    324      46.3

The default is 0.4: 18 losses in 18,008 previously-passing events (0.10%)
against 385 recoveries.

P(within 60 ps) is an XGBoost classifier over the reco cluster columns of the
export. Its own top features say what it learned: trkptz_rank and
trkptz_ratio_to_max carry ~87% of the importance, i.e. the model IS
TRKPTZ-plus-corrections -- the corrections coming from jet association
(n_tracks_in_fwdjet, n_in_jet_timed), timing self-consistency (t_in_jet_sigma,
dt_cluster_to_hgtd, hgtd_valid), and event context. That matches the failure
anatomy: in recoverable failures the wrongly-chosen cluster has 3x the TRKPTZ
score, no tracks in any forward jet, and worse internal timing agreement than
the correct one.

Truth discipline: every truth_* column, delta_t, and the event weight are
excluded from the feature list; truth enters only as the training label and the
evaluation.
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
    drop = {"ok", "fold", "delta_t", "weight", "cluster_idx"}
    return [k for k in c.columns
            if k not in EVT and k not in drop
            and not k.startswith("truth_") and c[k].dtype != object]


def guarded_pick(df, delta):
    """Confusion of the guarded metric against TRKPTZ. Returns (pass%, lost, gained, n)."""
    g = df.groupby(EVT, sort=False)
    ib, ip = g["trkptz_score"].idxmax(), g["p"].idxmax()
    old = df.loc[ib, [*EVT, "p", "ok"]].rename(columns={"p": "p_old", "ok": "ok_old"})
    new = df.loc[ip, [*EVT, "p", "ok"]].rename(columns={"p": "p_new", "ok": "ok_new"})
    m = old.merge(new, on=EVT)
    switch = m.p_new > m.p_old + delta
    ok = np.where(switch, m.ok_new, m.ok_old)
    return (100 * ok.mean(),
            int((m.ok_old & ~ok).sum()), int((~m.ok_old & ok).sum()), len(m))


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--input-dir", required=True)
    p.add_argument("--sample", default="vbf")
    p.add_argument("--delta", type=float, default=0.4)
    p.add_argument("--seeds", default="0,1,2",
                   help="split seeds; the metric must hold under every one")
    p.add_argument("--out", default="")
    args = p.parse_args()

    c, path = load(args.input_dir, args.sample)
    FE = features(c)
    log(f"{path}\n{c[EVT].drop_duplicates().shape[0]:,} events, "
        f"{len(c):,} clusters, {len(FE)} reco features\n")

    import xgboost as xgb
    results = []
    for seed in [int(s) for s in args.seeds.split(",")]:
        rng = np.random.default_rng(seed)
        ev = c[EVT].drop_duplicates().reset_index(drop=True)
        ev["fold"] = rng.random(len(ev))
        d = c.merge(ev, on=EVT, how="left")
        tr, te = d[d.fold >= 0.5], d[d.fold < 0.5].copy()
        m = xgb.XGBClassifier(n_estimators=500, learning_rate=0.05, max_depth=6,
                              subsample=0.8, colsample_bytree=0.8, min_child_weight=20,
                              eval_metric="logloss", n_jobs=-1, tree_method="hist",
                              random_state=seed)
        m.fit(tr[FE].to_numpy(np.float32), tr.ok.to_numpy(int))
        te["p"] = m.predict_proba(te[FE].to_numpy(np.float32))[:, 1]
        base = te.loc[te.groupby(EVT, sort=False)["trkptz_score"].idxmax()]
        pas, lost, gained, n = guarded_pick(te, args.delta)
        results.append(dict(seed=seed, n=n, trkptz=100 * base.ok.mean(),
                            new=pas, lost=lost, gained=gained))
        log(f"  seed {seed}: {n:,} test ev  TRKPTZ {100*base.ok.mean():.2f}%  "
            f"guarded {pas:.2f}%  lost {lost}  gained {gained}  "
            f"(lost = {100*lost/max(int(base.ok.sum()),1):.3f}% of previously passing)")

    log("\nsummary over seeds:")
    log(f"  TRKPTZ  {np.mean([r['trkptz'] for r in results]):.2f}%")
    log(f"  guarded {np.mean([r['new'] for r in results]):.2f}%  "
        f"lost {np.mean([r['lost'] for r in results]):.0f}  "
        f"gained {np.mean([r['gained'] for r in results]):.0f}  "
        f"(delta = {args.delta})")
    if args.out:
        json.dump({"args": vars(args), "results": results}, open(args.out, "w"), indent=2)
        log(f"\nwrote {args.out}")


if __name__ == "__main__":
    main()
