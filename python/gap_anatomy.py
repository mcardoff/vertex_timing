#!/usr/bin/env python3
"""
gap_anatomy.py — what is the 25.5-point Z+jets gap actually made of?

Studies B and C of the forward plan, in one pass because both need the same
expensive load (cluster tree + streamed tracks + the trained Deep Sets model).

B. THE SELECTOR LADDER, re-derived at the canonical selection.

   The decomposition everything currently rests on -- "perfect HS track identity
   with real times still stops at 80.2%" -- was measured at the OLD fiducial
   region, before the m_jj >= 500 selection and before the group-key fix. It is
   the number that splits the gap into TRACK IDENTITY versus TIME QUALITY, and
   quoting it against canonical results is exactly the kind of cross-selection
   mixing this project has a banner about. So it is recomputed here.

   The ladder, each rung a selector applied to the same events:

     TRKPTZ / WAVeS      the hand-designed scores
     Deep Sets           the trained model
     SUM pT . truth      PERFECT per-track HS identity, pooled by pT, REAL times.
                         Needs no model: it is what a flawless tagger would give
                         you under the current pooling rule.
     SUM pT . truth, timed
                         the same, but restricted to tracks with a valid time --
                         separates "the tagger is perfect" from "the tracks it
                         picks have no time to contribute".
     oracle              argmin |dt| over candidates: does a good cluster EXIST.

   The rung-to-rung differences are the answer. If Deep Sets -> SUM pT.truth is
   large, the remaining problem is identifying HS tracks. If that step is small
   and SUM pT.truth -> oracle is large, then perfect track identity does not
   help and the residual is the QUALITY of the times themselves -- which no
   selector can fix.

C. RANK AND CANDIDATE SEPARATION.

   Among events where a good candidate exists, where does Deep Sets rank it? The
   previous round found the right cluster inside the top 3 for 87% of events
   while picking it 67% of the time, which frames the problem as fine-grained
   ranking rather than gross misidentification -- but that was the GBDT, at the
   old selection. Repeated here for Deep Sets at canonical.

   Paired with efficiency against the z-distance to the nearest other candidate.
   If failures concentrate where two candidates sit close in z, that is a
   concrete mechanism rather than a correlation, and it is one a selector may
   genuinely be unable to resolve if the times themselves overlap.

USAGE
  python gap_anatomy.py --input-dir /data/mcardiff/training \\
      --model condor/runs/canonical --selection canonical --out condor/runs/anatomy.json
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


def pick_frac(frame, col, higher=True):
    """Argmax (or argmin) within event, then fraction inside 60 ps -- per sample."""
    g = frame.groupby(tds.EVT, sort=False)[col]
    idx = (g.idxmax() if higher else g.idxmin()).to_numpy()
    chosen = frame.loc[idx]
    return {tds.SAMPLE_NAME.get(s, str(s)): float(100 * sub["within60"].mean())
            for s, sub in chosen.groupby("sample_id")}


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--input-dir", required=True)
    p.add_argument("--model", required=True)
    p.add_argument("--selection", choices=["canonical", "loose"], default="canonical")
    p.add_argument("--seed", type=int, default=0, help="must match the model's --seed")
    p.add_argument("--test-events", type=int, default=100_000)
    p.add_argument("--step", default="300 MB")
    p.add_argument("--out", default="")
    args = p.parse_args()

    dev = torch.device("cpu")
    files = tds.find_files(args.input_dir, args.selection)
    if not files:
        sys.exit(f"FATAL: no exports under {args.input_dir}")

    ck = torch.load(os.path.join(args.model, "best_model.pt"),
                    map_location=dev, weights_only=False)

    # Cluster table, plus the context columns study C needs.
    extra = ["dz_to_nearest_cluster", "dt_to_nearest_cluster", "n_clusters",
             "n_valid_time", "cluster_time_sigma"]
    df = pd.concat([tds.read_tree(pth, "clusters",
                                  sorted(set(tds.CLUSTER_COLS + extra)))
                    for pth in files.values()], ignore_index=True)
    df["abs_dt"] = df["delta_t"].abs()
    df["within60"] = (df["abs_dt"] < tds.PASS_PS).astype(np.float32)

    rng = np.random.default_rng(args.seed)
    ev = df[tds.EVT].drop_duplicates().reset_index(drop=True)
    ev["fold"] = rng.random(len(ev))
    fold = ev.set_index(tds.EVT)["fold"]
    df = df.merge(ev, on=tds.EVT, how="left")

    # Test tracks, streamed exactly as the trainer does, so the event set matches.
    _, TE_T = tds.stream_tracks(files, fold, {s: 0 for s in files},
                                args.test_events // max(1, len(files)), args.step)

    # ---- Deep Sets scores on the test set --------------------------------
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

    # ---- the PERFECT-TAGGER rungs, computed from the track tree ----------
    # Sum pT over a cluster's tracks that truly come from the HS vertex. This is
    # what a flawless per-track tagger would hand the current pooling rule -- no
    # model involved, so it is a ceiling on tagging, not an estimate of one.
    TE_T = TE_T.copy()
    TE_T["_hs_pt"] = TE_T["pt"] * TE_T[tds.TRACK_LABEL]
    TE_T["_hs_pt_timed"] = TE_T["_hs_pt"] * (TE_T["time_valid"] > 0)
    agg = (TE_T.groupby(tds.KEY, sort=False)[["_hs_pt", "_hs_pt_timed"]]
                .sum().reset_index())

    t = df.merge(ds, on=tds.KEY, how="inner").merge(agg, on=tds.KEY, how="left")
    t[["_hs_pt", "_hs_pt_timed"]] = t[["_hs_pt", "_hs_pt_timed"]].fillna(0.0)
    log(f"\ntest set: {len(t):,} clusters / {t.groupby(tds.EVT).ngroups:,} events")

    names = ["vbf", "zjets", "dijet", "ttbar"]
    rungs = [
        ("TRKPTZ",                 "trkptz_score",  True),
        ("WAVeS",                  "waves_score",   True),
        ("Deep Sets",              "deepsets",      True),
        ("SUM pT . truth",         "_hs_pt",        True),
        ("SUM pT . truth, timed",  "_hs_pt_timed",  True),
        ("oracle  (|dt| argmin)",  "abs_dt",        False),
    ]
    log("\n" + "=" * 78)
    log("B -- SELECTOR LADDER at the canonical selection")
    log("=" * 78)
    log(f"{'rung':26s}" + "".join(f"{n:>12s}" for n in names))
    log("-" * 78)
    ladder = {}
    for lab, col, hi in rungs:
        r = pick_frac(t, col, hi)
        ladder[lab] = r
        log(f"{lab:26s}" + "".join(f"{r.get(n, float('nan')):>11.1f}%" for n in names))
    log("=" * 78)
    z = ladder["Deep Sets"].get("zjets"), ladder["SUM pT . truth"].get("zjets"), \
        ladder["oracle  (|dt| argmin)"].get("zjets")
    if all(v is not None for v in z):
        log(f"\nZ+jets split of the gap:")
        log(f"  Deep Sets -> perfect tagger : {z[1]-z[0]:+5.1f}  (TRACK IDENTITY)")
        log(f"  perfect tagger -> oracle    : {z[2]-z[1]:+5.1f}  (TIME QUALITY)")
        log("  The larger of the two is where the remaining gap actually lives.")

    # ---- C: rank of the best valid candidate ------------------------------
    log("\n" + "=" * 78)
    log("C -- RANK of the best valid candidate under Deep Sets")
    log("   (events that HAVE a valid candidate; rank 1 = we already pick it)")
    log("=" * 78)
    t["_rank"] = t.groupby(tds.EVT, sort=False)["deepsets"].rank(ascending=False, method="first")
    good = t[t["within60"] > 0]
    # EVT ALREADY BEGINS WITH sample_id -- appending it again duplicates the
    # column and reset_index() refuses to insert it.
    best_rank = good.groupby(tds.EVT, sort=False)["_rank"].min().reset_index()
    log(f"{'sample':10s}{'n events':>10s}" + "".join(f"{'top-'+str(k):>9s}" for k in (1,2,3,5,10)))
    ranks = {}
    for sid, sub in best_rank.groupby("sample_id"):
        nm = tds.SAMPLE_NAME.get(sid, str(sid))
        cum = {k: float(100 * (sub["_rank"] <= k).mean()) for k in (1,2,3,5,10)}
        ranks[nm] = cum
        log(f"{nm:10s}{len(sub):>10,}" + "".join(f"{cum[k]:>8.1f}%" for k in (1,2,3,5,10)))

    # ---- C: efficiency vs distance to the nearest other candidate ---------
    log("\n" + "=" * 78)
    log("C -- Deep Sets efficiency vs z-distance to the NEAREST other candidate")
    log("=" * 78)
    sel = t.loc[t.groupby(tds.EVT, sort=False)["deepsets"].idxmax().to_numpy()].copy()
    edges = [0, 0.2, 0.5, 1.0, 2.0, 5.0, np.inf]
    sel["_bin"] = pd.cut(sel["dz_to_nearest_cluster"], edges)
    lab = ["<0.2", "0.2-0.5", "0.5-1", "1-2", "2-5", ">5"]
    log(f"{'sample':10s}" + "".join(f"{l:>10s}" for l in lab))
    zsep = {}
    for sid, sub in sel.groupby("sample_id"):
        nm = tds.SAMPLE_NAME.get(sid, str(sid))
        row, counts = [], []
        for b in sub["_bin"].cat.categories:
            m = sub[sub["_bin"] == b]
            row.append(float(100 * m["within60"].mean()) if len(m) else float("nan"))
            counts.append(int(len(m)))
        zsep[nm] = {"eff": row, "n": counts}
        log(f"{nm:10s}" + "".join(f"{v:>9.1f}%" for v in row))
        log(f"{'  (n)':10s}" + "".join(f"{c:>10,}" for c in counts))
    log("=" * 78)
    log("A collapse in the left-hand bins is the concrete mechanism: two candidates")
    log("close in z, whose times may genuinely overlap within resolution.")

    if args.out:
        os.makedirs(os.path.dirname(os.path.abspath(args.out)) or ".", exist_ok=True)
        json.dump({"args": vars(args), "ladder": ladder, "rank": ranks,
                   "z_separation": zsep}, open(args.out, "w"), indent=2)
        log(f"\nwrote {args.out}")


if __name__ == "__main__":
    main()
