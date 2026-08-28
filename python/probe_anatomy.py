#!/usr/bin/env python3
"""
probe_anatomy.py — what carries the 84.2%, and where does the other 15.8% sit?

Study D established that a GBDT separates the chosen-but-wrong candidate from the
correct one 84.2% of the time on Z+jets (50% chance, 100% control). Two follow-ups
on that same measurement, neither needing a re-export.

1. FEATURE ABLATION. Retrain the probe on feature groups. This is targeted
   discovery on a measured effect, not generic feature engineering: the question
   is WHICH information carries the discrimination, and the answer says what the
   Deep Sets representation is failing to encode. If timing alone carries it while
   everything else sits near chance, the deficiency is in how timing reliability
   is represented. If timing is near chance and pT/z structure carries it, the
   direction is completely different.

   The two hand-designed scores are held out of every physics group and given
   their own row. They are functions of pT and z, so leaving them in would let
   "z only" smuggle in a pT-weighted score and blur exactly the split being
   measured.

2. ERROR ANATOMY. Where does the 15.8% live? Two very different worlds:

     concentrated at the 60 ps boundary  -> intrinsic ambiguity. A correct
       candidate at 55 ps against a wrong one at 65 ps is a distinction the
       metric cares about and the physics barely does; little to win.
     spread across obvious cases         -> a model deficiency. If the probe
       misses pairs where the wrong candidate is at 200 ps, there is real signal
       being left behind.

   Reported as accuracy binned by the |dt| of the WRONG candidate (how obviously
   wrong it is) and by the separation between the two.
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

GROUPS = {
    "timing": ["cluster_time", "cluster_time_sigma", "time_chi2_ndf", "max_abs_tpull",
               "time_spread", "n_valid_time", "frac_valid_time", "sumpt_valid_time",
               "n_tpull_gt2", "min_timeRes", "median_timeRes", "max_timeRes",
               "dt_to_nearest_cluster", "hgtd_time", "hgtd_valid", "hgtd_time_res",
               "dt_cluster_to_hgtd"],
    "z / vertex": ["delta_z", "delta_z_resunits", "cluster_z_sigma", "z_chi2_ndf",
                   "max_abs_zpull", "z_spread", "median_z0_pull", "dz_to_nearest_cluster"],
    "pT / multiplicity": ["n_tracks", "sumpt", "sumpt2", "maxpt", "pt_2nd", "lead_pt_frac",
                          "meanpt", "medianpt", "n_clusters"],
    "jets": ["n_tracks_in_fwdjet", "frac_pt_in_fwdjet", "min_dr_to_fwdjet",
             "pt_of_nearest_fwdjet", "n_ghost_tracks", "sumpt_ghost"],
    "track quality / HGTD hits": ["mean_nhgtd", "frac_nhgtd_ge2", "sumpt_w_nhgtd",
                                  "mean_quality", "mean_abs_eta", "eta_spread",
                                  "mean_d0_signif"],
    # Held separate: TRKPTZ and WAVeS are functions of pT and z, so leaving them
    # inside a physics group would let that group smuggle in a hand-tuned
    # combination and blur the very split being measured.
    "hand-designed scores": ["trkptz_score", "waves_score"],
}
ALL = sorted({f for v in GROUPS.values() for f in v})


def fit_acc(Xtr, ytr, Xte, yte):
    try:
        import xgboost as xgb
        m = xgb.XGBClassifier(n_estimators=400, learning_rate=0.06, max_depth=6,
                              subsample=0.8, colsample_bytree=0.8, min_child_weight=5,
                              eval_metric="logloss", n_jobs=-1, tree_method="hist")
    except ImportError:
        from sklearn.ensemble import HistGradientBoostingClassifier
        m = HistGradientBoostingClassifier(max_iter=300, learning_rate=0.06, random_state=0)
    m.fit(Xtr, ytr)
    return m, float((m.predict(Xte) == yte).mean())


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

    cols = sorted(set(tds.CLUSTER_COLS + ALL))
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
    t = df.merge(dsc, on=tds.KEY, how="inner")

    prng = np.random.default_rng(12345)
    out = {}
    for sample in ("zjets", "vbf", "dijet", "ttbar"):
        sid = [k for k, v in tds.SAMPLE_NAME.items() if v == sample]
        sub = t[t["sample_id"] == sid[0]] if sid else t.iloc[:0]
        if not len(sub):
            continue
        g = sub.groupby(tds.EVT, sort=False)
        pick = sub.loc[g["deepsets"].idxmax().to_numpy()]
        best = sub.loc[g["abs_dt"].idxmin().to_numpy()]
        m = pick[tds.EVT + ["within60"]].merge(
            best[tds.EVT + ["within60"]], on=tds.EVT, suffixes=("_pick", "_best"))
        fail = m[(m["within60_best"] > 0) & (m["within60_pick"] < 1)][tds.EVT]
        if len(fail) < 400:
            continue

        A = pick.merge(fail, on=tds.EVT)          # chosen, wrong
        B = best.merge(fail, on=tds.EVT)          # correct
        n = len(A); cut = int(0.6 * n)
        flip = prng.random(n) < 0.5
        y = (~flip).astype(int)                    # 1 => the first slot is the better one

        def build(feats):
            fa = np.nan_to_num(A[feats].to_numpy(np.float32))
            fb = np.nan_to_num(B[feats].to_numpy(np.float32))
            P = np.where(flip[:, None], fa, fb); Q = np.where(flip[:, None], fb, fa)
            return np.hstack([P, Q, P - Q]).astype(np.float32)

        log(f"\n{'='*72}\n{sample}  ({n:,} failure pairs)\n{'='*72}")
        res = {}
        Xall = build(ALL)
        model_all, acc_all = fit_acc(Xall[:cut], y[:cut], Xall[cut:], y[cut:])
        res["ALL"] = 100 * acc_all
        log(f"  {'ALL features':30s}{100*acc_all:7.1f}%   ({len(ALL)} features)")
        for gname, feats in GROUPS.items():
            X = build(feats)
            _, a = fit_acc(X[:cut], y[:cut], X[cut:], y[cut:])
            res[gname] = 100 * a
            log(f"  {gname:30s}{100*a:7.1f}%   ({len(feats)})")
        notime = [f for f in ALL if f not in GROUPS["timing"]]
        X = build(notime); _, a = fit_acc(X[:cut], y[:cut], X[cut:], y[cut:])
        res["ALL except timing"] = 100 * a
        log(f"  {'ALL except timing':30s}{100*a:7.1f}%   ({len(notime)})")

        # ---- error anatomy, on the held-out half ----
        pr = model_all.predict_proba(Xall[cut:])[:, 1]
        correct = (pr > 0.5).astype(int) == y[cut:]
        dt_wrong = A["abs_dt"].to_numpy()[cut:]     # |dt| of the chosen (wrong) one
        dt_right = B["abs_dt"].to_numpy()[cut:]     # |dt| of the correct one
        sep = dt_wrong - dt_right

        log(f"\n  error anatomy -- accuracy by how obviously wrong the chosen one is")
        edges = [60, 80, 120, 200, 400, np.inf]
        lab = ["60-80", "80-120", "120-200", "200-400", ">400"]
        anat = {}
        for lo, hi, l in zip(edges[:-1], edges[1:], lab):
            msk = (dt_wrong >= lo) & (dt_wrong < hi)
            if msk.sum() < 20:
                continue
            anat[l] = {"acc": float(100 * correct[msk].mean()), "n": int(msk.sum())}
            log(f"    |dt| of wrong candidate {l:>9s} ps : {100*correct[msk].mean():5.1f}%"
                f"   (n={msk.sum():,})")
        res["_anatomy_dt_wrong"] = anat

        log(f"\n  accuracy by separation |dt_wrong| - |dt_right|")
        sedges = [0, 25, 50, 100, 200, np.inf]; slab = ["<25", "25-50", "50-100", "100-200", ">200"]
        sanat = {}
        for lo, hi, l in zip(sedges[:-1], sedges[1:], slab):
            msk = (sep >= lo) & (sep < hi)
            if msk.sum() < 20:
                continue
            sanat[l] = {"acc": float(100 * correct[msk].mean()), "n": int(msk.sum())}
            log(f"    separation {l:>8s} ps : {100*correct[msk].mean():5.1f}%   (n={msk.sum():,})")
        res["_anatomy_separation"] = sanat
        out[sample] = res

    log("\n" + "=" * 72)
    log("Reading the ablation: a group near 50% carries nothing on its own; the")
    log("gap between ALL and ALL-except-timing is what timing adds that nothing")
    log("else supplies. Reading the anatomy: accuracy RISING with |dt| of the")
    log("wrong candidate means the residual is boundary ambiguity (little to")
    log("win); FLAT means the probe misses obvious cases too (signal left behind).")

    if args.out:
        os.makedirs(os.path.dirname(os.path.abspath(args.out)) or ".", exist_ok=True)
        json.dump({"args": vars(args), "results": out}, open(args.out, "w"), indent=2)
        log(f"\nwrote {args.out}")


if __name__ == "__main__":
    main()
