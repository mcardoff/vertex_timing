#!/usr/bin/env python3
"""
dump_failures.py — export the Z+jets events Deep Sets gets wrong, with every
candidate, for visual inspection.

Writes one JSON holding, per selected event, the full candidate table plus the
identifiers needed to reopen it: file_idx maps back to a real ntuple filename
(the exporter's file_idx is the position in the SORTED file list, which is how
setupChain builds the chain), and event_num is the entry within that file. That
pair is exactly what python/event_display.py takes, so anything picked here can
be reproduced as a full display rather than only inspected as numbers.

CATEGORIES. Failures are not all the same thing, and lumping them hides the
question worth asking:

  recoverable   a candidate within 60 ps EXISTS and Deep Sets picked a different
                one. A selection failure -- the interesting case.
  no-candidate  no cluster is within 60 ps at all. No selector can fix these;
                they are a clustering or timing limitation, and if they carry a
                common signature that is a statement about the input data rather
                than about the model.

Both are dumped and tagged so the two can be compared. Events are SAMPLED from
each category rather than taken in order: file order correlates with nothing
physical, but it does correlate with which files happened to be processed first.
"""
import argparse, glob, importlib.util, json, os, sys

import numpy as np
import pandas as pd
import torch

_spec = importlib.util.spec_from_file_location(
    "tds", os.path.join(os.path.dirname(os.path.abspath(__file__)), "train_deepsets.py"))
tds = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(tds)
log = tds.log

# Per-candidate columns worth having in front of you when diagnosing by eye.
CAND = ["cluster_idx", "cluster_time", "delta_t", "cluster_time_sigma",
        "delta_z", "delta_z_resunits", "cluster_z_sigma",
        "n_tracks", "sumpt", "maxpt", "n_valid_time", "frac_valid_time",
        "time_chi2_ndf", "max_abs_tpull", "time_spread",
        "min_dr_to_fwdjet", "n_tracks_in_fwdjet",
        "trkptz_score", "waves_score",
        "truth_purity", "truth_n_hs_tracks", "truth_hs_frac_tracks",
        "dz_to_nearest_cluster", "dt_to_nearest_cluster"]
EVENT = ["n_clusters", "n_forward_jets", "lead_jet_pt", "sublead_jet_pt",
         "n_fwd_tracks_reco", "event_sumpt", "n_reco_vertices", "reco_vtx_z",
         "hgtd_time", "hgtd_valid", "hgtd_time_res"]


def file_map(ntuple_dir):
    """Position in the sorted file list -> path, mirroring setupChain.

    setupChain descends ONE level into subdirectories (some samples store their
    ROOT files a directory deeper) and sorts the combined list; file_idx indexes
    that sorted list, so this has to reproduce the same walk or the mapping back
    to a real file is silently wrong."""
    paths = []
    if not os.path.isdir(ntuple_dir):
        return {}
    for e in sorted(os.listdir(ntuple_dir)):
        p = os.path.join(ntuple_dir, e)
        if os.path.isdir(p):
            for n in sorted(os.listdir(p)):
                q = os.path.join(p, n)
                if not os.path.isdir(q):
                    paths.append(q)
        else:
            paths.append(p)
    paths.sort()
    return {i: p for i, p in enumerate(paths)}


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--input-dir", required=True)
    p.add_argument("--model", required=True)
    p.add_argument("--sample", default="zjets")
    p.add_argument("--ntuple-dir", default="/data/mcardiff/exotic_superntuples/zjets/")
    p.add_argument("--selection", choices=["canonical", "loose"], default="canonical")
    p.add_argument("--seed", type=int, default=0)
    p.add_argument("--n-recoverable", type=int, default=100)
    p.add_argument("--n-nocandidate", type=int, default=30)
    p.add_argument("--test-events", type=int, default=100_000)
    p.add_argument("--step", default="300 MB")
    p.add_argument("--out", required=True)
    args = p.parse_args()

    dev = torch.device("cpu")
    files = tds.find_files(args.input_dir, args.selection)
    ck = torch.load(os.path.join(args.model, "best_model.pt"),
                    map_location=dev, weights_only=False)

    cols = sorted(set(tds.CLUSTER_COLS + CAND + EVENT))
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

    sid = [k for k, v in tds.SAMPLE_NAME.items() if v == args.sample][0]
    t = df.merge(dsc, on=tds.KEY, how="inner")
    t = t[t["sample_id"] == sid].copy()
    if not len(t):
        sys.exit(f"FATAL: no {args.sample} clusters in the test set")

    g = t.groupby(tds.EVT, sort=False)
    pick = t.loc[g["deepsets"].idxmax().to_numpy()]
    best = t.loc[g["abs_dt"].idxmin().to_numpy()]
    st = pick[tds.EVT + ["within60", "cluster_idx"]].merge(
        best[tds.EVT + ["within60", "cluster_idx"]], on=tds.EVT,
        suffixes=("_pick", "_best"))

    n_ev = len(st)
    ok = st["within60_pick"] > 0
    recoverable = st[(st["within60_best"] > 0) & (st["within60_pick"] < 1)]
    nocand = st[st["within60_best"] < 1]
    log(f"\n{args.sample} test events: {n_ev:,}")
    log(f"  Deep Sets correct          : {int(ok.sum()):,}  ({100*ok.mean():.1f}%)")
    log(f"  FAILURE, recoverable       : {len(recoverable):,}  ({100*len(recoverable)/n_ev:.1f}%)")
    log(f"  FAILURE, no valid candidate: {len(nocand):,}  ({100*len(nocand)/n_ev:.1f}%)")

    fmap = file_map(args.ntuple_dir)
    log(f"  ntuple files discovered    : {len(fmap):,}"
        + ("" if fmap else "  (ntuple dir unreadable -- file names will be null)"))

    srng = np.random.default_rng(99)
    out = {"sample": args.sample, "selection": args.selection,
           "counts": {"test_events": int(n_ev), "correct": int(ok.sum()),
                      "recoverable_failures": int(len(recoverable)),
                      "no_candidate_failures": int(len(nocand))},
           "events": []}

    for tag, tbl, want in (("recoverable", recoverable, args.n_recoverable),
                           ("no_candidate", nocand, args.n_nocandidate)):
        if not len(tbl):
            continue
        take = tbl.iloc[srng.choice(len(tbl), min(want, len(tbl)), replace=False)]
        for _, row in take.iterrows():
            key = {c: float(row[c]) for c in tds.EVT}
            cl = t[(t["sample_id"] == row["sample_id"]) &
                   (t["file_idx"] == row["file_idx"]) &
                   (t["event_num"] == row["event_num"])].sort_values("cluster_idx")
            fidx = int(row["file_idx"])
            path = fmap.get(fidx)
            rec = {
                "category": tag,
                "key": key,
                "file_idx": fidx,
                "entry": int(row["event_num"]),
                "ntuple_file": os.path.basename(path) if path else None,
                "picked_cluster": int(row["cluster_idx_pick"]),
                "correct_cluster": int(row["cluster_idx_best"]),
                "event": {c: (None if c not in cl.columns or not np.isfinite(cl.iloc[0][c])
                              else float(cl.iloc[0][c])) for c in EVENT},
                "candidates": [],
            }
            for _, c in cl.iterrows():
                rec["candidates"].append(
                    {k: (None if (k not in cl.columns or not np.isfinite(c[k])) else float(c[k]))
                     for k in CAND}
                    | {"deepsets": float(c["deepsets"]),
                       "is_picked": bool(int(c["cluster_idx"]) == int(row["cluster_idx_pick"])),
                       "is_correct": bool(int(c["cluster_idx"]) == int(row["cluster_idx_best"])),
                       "within60": bool(c["within60"] > 0)})
            out["events"].append(rec)

    os.makedirs(os.path.dirname(os.path.abspath(args.out)) or ".", exist_ok=True)
    json.dump(out, open(args.out, "w"), indent=1)
    log(f"\nwrote {args.out}  ({len(out['events'])} events, "
        f"{sum(len(e['candidates']) for e in out['events']):,} candidates)")


if __name__ == "__main__":
    main()
