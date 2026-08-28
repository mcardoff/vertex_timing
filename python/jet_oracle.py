#!/usr/bin/env python3
"""
jet_oracle.py — if you knew perfectly which forward jets were hard-scatter,
how much of the Z+jets gap would close?

This is the cheap gate on a whole line of work. WAVeS weights tracks by proximity
to forward jets on the assumption that those jets are hard-scatter objects. In
Z+jets they are mostly NOT -- the canonical measurement puts the forward dijet at
~94% pileup-supplied -- which is why WAVeS gains on VBF and loses 1.9 points on
Z+jets. So the question is whether a jet HS/PU discriminator would fix it.

Answering that with a MODEL first would be backwards. This computes the ceiling
directly, using truth jet identity, with no model at all:

  SUM pT, HS jets     sum pT over the cluster's tracks whose nearest forward jet
                      is truth-matched to a hard-scatter jet
  SUM pT, PU jets     the same for pileup jets -- the control. If jet identity
                      carries the signal, these two must separate; if BOTH score
                      well, the discrimination is coming from jet proximity
                      generally rather than from HS/PU identity, and a tagger
                      would buy nothing.
  WAVeS-form, HS only the actual WAVeS weighting (pT_trk * pT_jet / max(dR, floor),
                      damped by exp(-1.5|dz|)) restricted to HS-jet tracks. This
                      is what WAVeS would score if its central assumption held.

Read against the existing rungs. If "WAVeS-form, HS only" lands near the
Sum-pT-truth rung (79.4% on Z+jets), jet identity recovers most of the
track-identity half and a reco jet tagger is worth building. If it lands near
Deep Sets, knowing the jets does not help and neither attention supervision nor a
tagger is justified.

TRUTH IS USED ONLY TO COMPUTE THE CEILING. Nothing here is a model input; the
point is to measure what a perfect tagger would be worth before anyone builds an
imperfect one.
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

WAVES_DR_FLOOR = 0.05      # matches clustering_constants.h
EXTRA_TRK = ["truth_nearest_fwdjet_is_hs"]


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

    df = pd.concat([tds.read_tree(pth, "clusters", tds.CLUSTER_COLS + ["delta_z"])
                    for pth in files.values()], ignore_index=True)
    df["abs_dt"] = df["delta_t"].abs()
    df["within60"] = (df["abs_dt"] < tds.PASS_PS).astype(np.float32)

    rng = np.random.default_rng(args.seed)
    ev = df[tds.EVT].drop_duplicates().reset_index(drop=True)
    ev["fold"] = rng.random(len(ev))
    fold = ev.set_index(tds.EVT)["fold"]
    df = df.merge(ev, on=tds.EVT, how="left")

    # stream_tracks reads TRACK_FEATURES, which deliberately excludes truth. The
    # oracle needs one truth column, so extend the list for the READ and restore
    # it before make_tensors -- otherwise the extra column would enter the model's
    # input tensor, which is exactly the leak this project guards against.
    orig = list(tds.TRACK_FEATURES)
    tds.TRACK_FEATURES = orig + EXTRA_TRK
    try:
        _, TE_T = tds.stream_tracks(files, fold, {s: 0 for s in files},
                                    args.test_events // max(1, len(files)), args.step)
    finally:
        tds.TRACK_FEATURES = orig
    if EXTRA_TRK[0] not in TE_T.columns:
        sys.exit(f"FATAL: {EXTRA_TRK[0]} absent -- this export predates it; re-export needed")

    labels = df[tds.KEY + ["within60"]].drop_duplicates(tds.KEY)
    TST = tds.make_tensors(TE_T, pd.Series(ck["mu"]), pd.Series(ck["sd"]),
                           ck["na_cols"], labels, dev)
    assert TST["X"].shape[1] == len(orig) + len(ck["na_cols"]), \
        "truth column leaked into the model input tensor"
    pool = ck["key"][0] if isinstance(ck.get("key"), (list, tuple)) else "sum"
    net = tds.DeepSets(len(ck["features"]) + len(ck["na_cols"]),
                       ck["hidden"], ck["embed"], pool=pool)
    net.load_state_dict(ck["state_dict"]); net.eval()
    with torch.no_grad():
        s = net(TST["X"], TST["c"], TST["ncl"]).cpu().numpy()
    dsc = TST["m"][tds.KEY].copy(); dsc["deepsets"] = s

    # ---- per-cluster jet-oracle sums, straight from the track table ----
    T = TE_T.copy()
    isHS = (T[EXTRA_TRK[0]] > 0.5).astype(float)
    inJet = T["dr_nearest_fwdjet"].notna().astype(float)
    T["_hs_pt"] = T["pt"] * isHS * inJet
    T["_pu_pt"] = T["pt"] * (1.0 - isHS) * inJet
    T["_hs_tot"] = T["pt"] * T[tds.TRACK_LABEL]          # per-TRACK truth, for reference
    dr = T["dr_nearest_fwdjet"].clip(lower=WAVES_DR_FLOOR)
    T["_hs_waves"] = (T["pt"] * T["pt_nearest_fwdjet"] / dr).fillna(0.0) * isHS * inJet
    agg = (T.groupby(tds.KEY, sort=False)[["_hs_pt", "_pu_pt", "_hs_tot", "_hs_waves"]]
             .sum().reset_index())

    t = df.merge(dsc, on=tds.KEY, how="inner").merge(agg, on=tds.KEY, how="left")
    for c in ("_hs_pt", "_pu_pt", "_hs_tot", "_hs_waves"):
        t[c] = t[c].fillna(0.0)
    # WAVeS damps by the cluster's z compatibility; apply the same factor so the
    # comparison is against WAVeS's actual functional form, not a piece of it.
    t["_hs_waves_dz"] = t["_hs_waves"] * np.exp(-1.5 * t["delta_z"].abs())

    def cf(col, hi=True):
        g = t.groupby(tds.EVT, sort=False)[col]
        idx = (g.idxmax() if hi else g.idxmin()).to_numpy()
        ch = t.loc[idx]
        return {tds.SAMPLE_NAME.get(k, str(k)): float(100 * v["within60"].mean())
                for k, v in ch.groupby("sample_id")}

    rungs = [("TRKPTZ", "trkptz_score", True),
             ("WAVeS (as shipped)", "waves_score", True),
             ("Deep Sets", "deepsets", True),
             ("SUM pT, tracks in PU jets", "_pu_pt", True),
             ("SUM pT, tracks in HS jets", "_hs_pt", True),
             ("WAVeS-form, HS jets only", "_hs_waves_dz", True),
             ("SUM pT . truth (per track)", "_hs_tot", True),
             ("oracle (|dt| argmin)", "abs_dt", False)]
    names = ["vbf", "zjets", "dijet", "ttbar"]
    log("\n" + "=" * 80)
    log("JET-IDENTITY ORACLE -- what perfect HS/PU jet knowledge would be worth")
    log("=" * 80)
    log(f"{'rung':30s}" + "".join(f"{n:>12s}" for n in names))
    log("-" * 80)
    out = {}
    for lab, col, hi in rungs:
        r = cf(col, hi); out[lab] = r
        log(f"{lab:30s}" + "".join(f"{r.get(n, float('nan')):>11.1f}%" for n in names))
    log("=" * 80)
    z = {k: v.get("zjets") for k, v in out.items()}
    if z.get("Deep Sets") and z.get("WAVeS-form, HS jets only"):
        log(f"\nZ+jets: Deep Sets {z['Deep Sets']:.1f}  ->  "
            f"WAVeS-form on HS jets {z['WAVeS-form, HS jets only']:.1f}  "
            f"({z['WAVeS-form, HS jets only']-z['Deep Sets']:+.1f})")
        log(f"        per-track perfect tagger {z.get('SUM pT . truth (per track)', float('nan')):.1f}"
            f"   oracle {z.get('oracle (|dt| argmin)', float('nan')):.1f}")
        log("\nIf the HS-jet rung is near the per-track tagger, jet identity carries")
        log("most of the recoverable half and a reco jet tagger is worth building.")
        log("If it is near Deep Sets, it is not. And if the PU-jet rung scores")
        log("comparably, the signal is jet PROXIMITY, not jet IDENTITY.")

    if args.out:
        os.makedirs(os.path.dirname(os.path.abspath(args.out)) or ".", exist_ok=True)
        json.dump({"args": vars(args), "ladder": out}, open(args.out, "w"), indent=2)
        log(f"\nwrote {args.out}")


if __name__ == "__main__":
    main()
