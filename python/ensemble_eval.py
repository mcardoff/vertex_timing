#!/usr/bin/env python3
"""
ensemble_eval.py — does averaging several seeds' scores select better clusters?

WHY THIS IS WORTH A SHOT, and why it is not a violation of the frozen
architecture. Four models differing ONLY in seed spread 2.4 points of Z+jets
core fraction (65.1-67.5), against an infrastructure floor of exactly 0.000 --
four runs at one seed are bit-identical. So that 2.4 points is not measurement
noise, it is genuine variance between equally-valid models: each lands in a
different minimum and gets a different set of events wrong. That is the textbook
condition under which averaging their scores recovers accuracy, and it is a
change to the SELECTION PROCEDURE, not to the model. The architecture is
untouched, and item 9's bar -- "reopen only with a concrete reason from a study"
-- is met by the measurement above.

LEAKAGE IS THE WHOLE DESIGN PROBLEM. --seed in train_deepsets.py fixes the event
SPLIT, so two runs at different --seed have different test folds and member B
will have trained on events sitting in member A's test fold. Averaging those and
quoting the result is self-deception. Members must therefore share --seed and
differ only in --init-seed, which fixes weights and batch order. This script
REFUSES to combine models whose recorded --seed differs, rather than trusting the
caller to have got it right.

COMBINING RULE. Raw scores are not comparable across members: the listwise
softmax is invariant to a per-event additive shift, so each member's scores live
on its own arbitrary offset and scale. Averaging them raw lets whichever member
happens to have the widest spread dominate. Two rules are reported:

  rank  -- average of each member's within-event percentile rank. Scale-free,
           so no member can dominate; the natural choice for a pure ranking task.
  zsum  -- average of each member's within-event standardised score. Keeps some
           magnitude information, so a member that is confident about one cluster
           can outvote members that are indifferent.

Both are reported against every single member and against the best member, so
the gain (if any) is quoted against the right baseline -- ensembling has to beat
the BEST member to be worth its cost, not the average one.

USAGE
  python ensemble_eval.py --input-dir /data/mcardiff/training \\
      --models runs/ens_i0,runs/ens_i1,runs/ens_i2,runs/ens_i3 \\
      --selection canonical --out runs/ensemble.json
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


def core_fraction_by_sample(frame, col):
    """Argmax within event, then the fraction inside 60 ps -- per sample.

    Identical rule to the trainer's own reporting, so ensemble and member numbers
    are directly comparable rather than nearly so."""
    out = {}
    for sid, sub in frame.groupby("sample_id"):
        pick = sub.loc[sub.groupby(tds.EVT, sort=False)[col].idxmax()]
        out[tds.SAMPLE_NAME.get(sid, str(sid))] = float(100 * pick["within60"].mean())
    return out


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--input-dir", required=True)
    p.add_argument("--models", required=True,
                   help="comma-separated run directories, each with best_model.pt")
    p.add_argument("--selection", choices=["canonical", "loose"], default="canonical")
    p.add_argument("--test-events", type=int, default=100_000)
    p.add_argument("--step", default="300 MB")
    p.add_argument("--out", default="")
    args = p.parse_args()

    dirs = [d.strip() for d in args.models.split(",") if d.strip()]
    if len(dirs) < 2:
        sys.exit("FATAL: --models needs at least two run directories")

    dev = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    cks, split_seeds = [], set()
    for d in dirs:
        f = os.path.join(d, "best_model.pt")
        rj = os.path.join(d, "results.json")
        if not os.path.exists(f):
            sys.exit(f"FATAL: no best_model.pt in {d}")
        ck = torch.load(f, map_location=dev, weights_only=False)
        a = json.load(open(rj))["args"] if os.path.exists(rj) else {}
        split_seeds.add(a.get("seed"))
        cks.append((d, ck, a))
        log(f"  {d:24s} split-seed {a.get('seed')}  init-seed {a.get('init_seed')}"
            f"  selection {a.get('selection')}")

    # The refusal that keeps this honest -- see the leakage note in the docstring.
    if len(split_seeds) > 1:
        sys.exit(f"FATAL: members disagree on --seed {sorted(split_seeds)}. They have "
                 "different test folds, so combining them leaks training events into "
                 "the evaluation. Retrain with one --seed and differing --init-seed.")
    sels = {a.get("selection") for _, _, a in cks}
    if len(sels) > 1:
        sys.exit(f"FATAL: members trained on different selections {sorted(sels)}")

    seed = split_seeds.pop()
    log(f"\nall {len(cks)} members share split-seed {seed} -- test folds are identical")

    # Rebuild the shared split exactly as the trainer does.
    files = tds.find_files(args.input_dir, args.selection)
    df = tds.load_clusters(files)
    rng = np.random.default_rng(seed)
    ev = df[tds.EVT].drop_duplicates().reset_index(drop=True)
    ev["fold"] = rng.random(len(ev))
    fold = ev.set_index(tds.EVT)["fold"]
    df = df.merge(ev, on=tds.EVT, how="left")
    labels = df[tds.KEY + ["within60"]].drop_duplicates(tds.KEY)

    _, TE_T = tds.stream_tracks(files, fold, {s: 0 for s in files},
                                args.test_events // max(1, len(files)), args.step)

    base = None
    for d, ck, a in cks:
        TST = tds.make_tensors(TE_T, pd.Series(ck["mu"]), pd.Series(ck["sd"]),
                               ck["na_cols"], labels, dev)
        pool = ck["key"][0] if isinstance(ck.get("key"), (list, tuple)) else "sum"
        net = tds.DeepSets(len(ck["features"]) + len(ck["na_cols"]),
                           ck["hidden"], ck["embed"], pool=pool)
        net.load_state_dict(ck["state_dict"]); net = net.to(dev).eval()
        with torch.no_grad():
            s = net(TST["X"], TST["c"], TST["ncl"]).cpu().numpy()
        col = TST["m"][tds.KEY].copy()
        col[f"s_{os.path.basename(d)}"] = s
        base = col if base is None else base.merge(col, on=tds.KEY, how="inner")

    scored = base.merge(labels, on=tds.KEY, how="left")
    members = [c for c in scored.columns if c.startswith("s_")]

    # Per-event normalisation, then average. See the combining-rule note above.
    g = scored.groupby(tds.EVT, sort=False)
    for c in members:
        scored[f"r_{c}"] = g[c].rank(pct=True)
        # A single-cluster event has std = NaN (ddof=1 over one value), not 0, so
        # replace(0, ...) alone leaves it NaN, the whole group becomes NaN, and
        # idxmax returns NaN rather than a row label. Such an event is trivially
        # its own argmax, so any finite constant is correct; 1.0 makes z = 0.
        mu_ = g[c].transform("mean")
        sd_ = g[c].transform("std").replace(0, 1.0).fillna(1.0)
        scored[f"z_{c}"] = (scored[c] - mu_) / sd_
    scored["ens_rank"] = scored[[f"r_{c}" for c in members]].mean(axis=1)
    scored["ens_zsum"] = scored[[f"z_{c}" for c in members]].mean(axis=1)

    per_member = {c[2:]: core_fraction_by_sample(scored, c) for c in members}
    ens_rank = core_fraction_by_sample(scored, "ens_rank")
    ens_zsum = core_fraction_by_sample(scored, "ens_zsum")
    names = sorted({s for v in per_member.values() for s in v})

    log("\n" + "=" * 84)
    log(f"{'model':30s}" + "".join(f"{n:>12s}" for n in names))
    log("-" * 84)
    for m, v in per_member.items():
        log(f"{m:30s}" + "".join(f"{v.get(n, float('nan')):>11.1f}%" for n in names))
    best = {n: max(v.get(n, float("-inf")) for v in per_member.values()) for n in names}
    log(f"{'BEST single member':30s}" + "".join(f"{best[n]:>11.1f}%" for n in names))
    log("-" * 84)
    log(f"{'ensemble (rank avg)':30s}" + "".join(f"{ens_rank.get(n, float('nan')):>11.1f}%" for n in names))
    log(f"{'ensemble (z avg)':30s}" + "".join(f"{ens_zsum.get(n, float('nan')):>11.1f}%" for n in names))
    log("-" * 84)
    log(f"{'gain of best ens vs best':30s}"
        + "".join(f"{max(ens_rank.get(n, -1e9), ens_zsum.get(n, -1e9)) - best[n]:>+11.1f} " for n in names))
    log("=" * 84)
    log("\nThe gain line is against the BEST single member, not the average one:")
    log("an ensemble that only beats the mean member is not worth N times the")
    log("compute, since you could have trained one model and kept it.")

    if args.out:
        os.makedirs(os.path.dirname(os.path.abspath(args.out)) or ".", exist_ok=True)
        json.dump({"args": vars(args), "split_seed": seed, "members": per_member,
                   "best_single": best, "ensemble_rank": ens_rank,
                   "ensemble_zsum": ens_zsum}, open(args.out, "w"), indent=2)
        log(f"\nwrote {args.out}")


if __name__ == "__main__":
    main()
