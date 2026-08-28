#!/usr/bin/env python3
"""
e2e_t0.py — learn per-track weights and the vertex t0 TOGETHER, against the
60 ps criterion itself.

WHY THIS EXISTS. Three things have now failed in this study for the same reason,
and it is not the architecture:

  |dt|-ordered pairwise ranker   -4.1   trained to order by |dt| when the metric
                                        rewards ANY candidate inside 60 ps
  purity-filtered labels         -2.8   removed positives the metric counts
  two-stage tagger + median       ---   stage 1 maximises AUC, nothing in the
                                        pipeline ever optimises the metric

The two-stage evidence is explicit: on Z+jets, tightening the tagger threshold
makes t0 WORSE (60.0 -> 53.4 -> 40.6) even though tighter means purer, because
purity starves the estimator of tracks. AUC and the objective actively disagree.

So this model optimises the objective directly. It emits a weight per track,
combines the track times into one t0 with a differentiable robust estimator, and
is trained on a smooth surrogate of "is |t0 - t_truth| < 60 ps".

THE ESTIMATOR IS THE POINT. An inverse-variance mean scores 65.3% on perfectly
tagged HS tracks while a plain median scores 81.7%, because 13.7% of true HS
tracks sit >200 ps from truth against a ~25 ps quoted resolution. So the
combiner here is a WEIGHTED L1 MEDIAN, reached by unrolled IRLS:

    t <- SUM (w_i / |t_i - t|) t_i  /  SUM (w_i / |t_i - t|)

which is differentiable in both w and t, and whose fixed point is the weighted
median. Five unrolled steps from a weighted-mean start.

Ceilings for reference (canonical Z+jets, held out):
    TRKPTZ 61.9 | Deep Sets 67.3 | two-stage aggregate 63.0 | perfect tagger 81.7

TRUTH USE. truth_is_hs is NOT used at all -- there is no tagging label here, only
the event's t0. The reported core fraction is |t0 - t_truth| < 60 ps on held-out
events.
"""
import argparse, glob, json, os, sys

import numpy as np
import pandas as pd
import torch
import torch.nn as nn
import uproot

EVT = ["sample_id", "file_idx", "event_num"]
KEY = EVT + ["cluster_idx"]
PASS_PS = 60.0
SAMPLE_NAME = {-1.0: "local", 0.0: "vbf", 1.0: "zjets", 2.0: "dijet", 3.0: "zeejets",
               4.0: "vbf_mu0", 5.0: "zeejets_mu0", 6.0: "ttbar_mu0", 7.0: "ttbar"}

FEATS = ["pt", "eta", "z0", "d0", "sigma_z0", "sigma_d0", "time", "timeRes",
         "time_valid", "quality", "nhgtd_hits", "z0_pull_pv", "t_pull_cluster",
         "loo_shift", "loo_pull", "dr_nearest_fwdjet", "pt_nearest_fwdjet",
         "is_ghost_of_nearest", "is_lepton", "cluster_time", "cluster_delta_z",
         "dr_nearest_anyjet", "pt_nearest_anyjet", "is_in_any_jet",
         "dz_to_nearest_pu_vtx_trk", "closer_to_pu_than_pv",
         "cluster_dz_to_nearest_pu_vtx", "cluster_pv_pu_dz_ratio",
         "cluster_nearest_vtx_waves_rank", "cluster_nearest_vtx_waves_frac",
         "cluster_closest_vtx_is_pv", "on_pv", "on_pu_vtx", "vtx_unassigned",
         "vtx_weight", "vtx_sumpt2_frac", "cluster_frac_trk_on_pv",
         "cluster_dominant_vtx_is_pv", "cluster_dominant_vtx_sumpt2_frac"]


def log(m): print(m, flush=True)


class WeightNet(nn.Module):
    """Per-track weight, then a weighted L1 median over each event's tracks.

    The weight is softplus so it is strictly positive: a track can be made
    irrelevant but never negative, which would let the model place t0 outside
    the range of its own measurements.
    """
    def __init__(self, n_in, hidden=96, iters=5, floor=10.0):
        super().__init__()
        self.iters = iters
        # Residual floor for the IRLS reweighting, in ps. NOT a numerical
        # epsilon: at 1e-3 a track sitting on the current estimate gets weight
        # w/1e-3 and the gradient explodes -- that made validation oscillate
        # between 37% and 74% and settle 20 points below its own best. A floor
        # at the timing-resolution scale Winsorises the estimator instead,
        # putting it between the L1 median (81.7% on perfect tracks) and the
        # truncated mean (82.3%), which is where the good estimators live.
        self.floor = floor
        self.f = nn.Sequential(
            nn.Linear(n_in, hidden), nn.ReLU(),
            nn.Linear(hidden, hidden), nn.ReLU(),
            nn.Linear(hidden, 1))

    def forward(self, x, t, eidx, n_ev):
        w = torch.nn.functional.softplus(self.f(x)).squeeze(-1) + 1e-6
        num = torch.zeros(n_ev, device=x.device).index_add_(0, eidx, w * t)
        den = torch.zeros(n_ev, device=x.device).index_add_(0, eidx, w)
        t0 = num / den.clamp_min(1e-9)                      # weighted mean start
        for _ in range(self.iters):
            # IRLS step toward the weighted L1 median. The residual is DETACHED:
            # differentiating through the reweighting as well as through the
            # estimate makes an unrolled fixed-point iteration ill-conditioned.
            # Detaching gives the standard fixed-point gradient and is what
            # turns this from divergent into stable.
            r = (t - t0[eidx]).abs().detach().clamp_min(self.floor)
            u = w / r
            num = torch.zeros(n_ev, device=x.device).index_add_(0, eidx, u * t)
            den = torch.zeros(n_ev, device=x.device).index_add_(0, eidx, u)
            t0 = num / den.clamp_min(1e-9)
        return t0, w


def window_loss(t0, truth, tau):
    """Smooth surrogate for |t0 - truth| < PASS_PS.

    Anneal tau from wide to narrow: early on the gradient must reach events that
    are hundreds of ps out, later it should concentrate on the boundary. A pure
    step function has zero gradient everywhere and cannot train at all.
    """
    return (1.0 - torch.sigmoid((PASS_PS - (t0 - truth).abs()) / tau)).mean()


def load(path, max_events, cols):
    c = uproot.open(path)["clusters"].arrays(
        KEY + ["delta_t", "cluster_time", "trkptz_score", "waves_score"], library="pd")
    c["t_truth"] = c["cluster_time"] - c["delta_t"]
    truth = c.groupby(EVT, sort=False)["t_truth"].first()
    if max_events and len(truth) > max_events:
        truth = truth.iloc[:max_events]
        c = c.set_index(EVT).loc[truth.index].reset_index()
    t = uproot.open(path)["tracks"].arrays(KEY + ["track_idx"] + cols, library="pd")
    t = t.merge(truth.rename("t_truth").reset_index(), on=EVT, how="inner")
    t = t.drop_duplicates(EVT + ["track_idx"])
    t = t[(t.time_valid > 0) & (t.timeRes > 0)]
    return c, t


def tensors(df, cols, mu, sd, dev):
    """Flat track tensor + an event index, so no padding is needed."""
    d = df.sort_values(EVT).reset_index(drop=True)
    codes, uniq = pd.factorize(pd.MultiIndex.from_frame(d[EVT]), sort=False)
    x = ((d[cols] - mu) / sd).fillna(0.0).to_numpy(np.float32)
    return (torch.from_numpy(x).to(dev),
            torch.from_numpy(d["time"].to_numpy(np.float32)).to(dev),
            torch.from_numpy(codes.astype(np.int64)).to(dev),
            torch.from_numpy(d.groupby(EVT, sort=False)["t_truth"].first()
                              .to_numpy(np.float32)).to(dev),
            len(uniq), pd.MultiIndex.from_tuples(uniq, names=EVT))


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--input-dir", required=True)
    p.add_argument("--samples", default="vbf,zjets,dijet,ttbar")
    p.add_argument("--max-events", type=int, default=40000)
    p.add_argument("--epochs", type=int, default=40)
    p.add_argument("--lr", type=float, default=5e-4)
    p.add_argument("--hidden", type=int, default=96)
    p.add_argument("--iters", type=int, default=5)
    p.add_argument("--floor", type=float, default=10.0,
                   help="IRLS residual floor in ps; see WeightNet")
    p.add_argument("--clip", type=float, default=1.0)
    p.add_argument("--batch-events", type=int, default=2048)
    p.add_argument("--seed", type=int, default=0)
    p.add_argument("--init-seed", type=int, default=None)
    p.add_argument("--torch-threads", type=int, default=8)
    p.add_argument("--out", default="")
    args = p.parse_args()
    if args.init_seed is None:
        args.init_seed = args.seed
    torch.set_num_threads(args.torch_threads)
    dev = torch.device("cpu")

    first = None
    for s in args.samples.split(","):
        hit = (glob.glob(os.path.join(args.input_dir, f"{s}_*training.root"))
               or glob.glob(os.path.join(args.input_dir, s, f"{s}_*training.root")))
        if hit:
            first = hit[0]; break
    if first is None:
        sys.exit("no samples found")
    have = set(uproot.open(first)["tracks"].keys())
    cols = [f for f in FEATS if f in have]
    miss = [f for f in FEATS if f not in have]
    if miss:
        log(f"  !! absent from this export, dropped: {miss}")

    C, T = {}, {}
    for s in args.samples.split(","):
        hit = (glob.glob(os.path.join(args.input_dir, f"{s}_*training.root"))
               or glob.glob(os.path.join(args.input_dir, s, f"{s}_*training.root")))
        if not hit:
            continue
        C[s], T[s] = load(hit[0], args.max_events, cols)
        log(f"  {s:8s} {C[s][EVT].drop_duplicates().shape[0]:>7,} ev  {len(T[s]):>9,} timed tracks")

    rng = np.random.default_rng(args.seed)
    for s in T:
        ev = T[s][EVT].drop_duplicates().reset_index(drop=True)
        ev["fold"] = rng.random(len(ev))
        T[s] = T[s].merge(ev, on=EVT, how="left")
        C[s] = C[s].merge(ev, on=EVT, how="left")
    TR = pd.concat([T[s][T[s].fold >= 0.4] for s in T], ignore_index=True)
    VA = pd.concat([T[s][(T[s].fold >= 0.3) & (T[s].fold < 0.4)] for s in T], ignore_index=True)

    mu = TR[cols].mean(); sd = TR[cols].std().replace(0, 1.0)
    Xtr, ttr, etr, ytr, ntr, _ = tensors(TR, cols, mu, sd, dev)
    Xva, tva, eva, yva, nva, _ = tensors(VA, cols, mu, sd, dev)
    log(f"\n{len(cols)} features | train {ntr:,} ev / {len(TR):,} tracks | val {nva:,} ev")

    torch.manual_seed(args.init_seed)
    net = WeightNet(len(cols), args.hidden, args.iters, args.floor).to(dev)
    opt = torch.optim.Adam(net.parameters(), lr=args.lr)
    log(f"parameters: {sum(p.numel() for p in net.parameters()):,}\n")

    best, best_state = -1.0, None
    for ep in range(1, args.epochs + 1):
        # Wide early so distant events still produce gradient, narrow later so
        # the model works the actual boundary.
        tau = float(np.interp(ep, [1, args.epochs], [80.0, 30.0]))
        net.train()
        perm = torch.randperm(ntr)
        for i in range(0, ntr, args.batch_events):
            sel = perm[i:i + args.batch_events]
            mask = torch.isin(etr, sel)
            if not mask.any():
                continue
            remap = torch.full((ntr,), -1, dtype=torch.long)
            remap[sel] = torch.arange(len(sel))
            t0, _ = net(Xtr[mask], ttr[mask], remap[etr[mask]], len(sel))
            loss = window_loss(t0, ytr[sel], tau)
            opt.zero_grad(); loss.backward()
            torch.nn.utils.clip_grad_norm_(net.parameters(), args.clip)
            opt.step()
        net.eval()
        with torch.no_grad():
            t0v, _ = net(Xva, tva, eva, nva)
            acc = float(((t0v - yva).abs() < PASS_PS).float().mean() * 100)
        star = ""
        if acc > best:
            best, best_state = acc, {k: v.clone() for k, v in net.state_dict().items()}
            star = " *"
        if ep % 5 == 0 or ep == 1 or star:
            log(f"  epoch {ep:3d}/{args.epochs}  tau {tau:5.1f}  VAL {acc:5.2f}%{star}")
    net.load_state_dict(best_state)
    log(f"\nrestored best val {best:.2f}%\n")

    log("=" * 74)
    log("END-TO-END t0 -- core fraction on HELD-OUT events")
    log("=" * 74)
    res = {}
    names = sorted(T)
    for s in names:
        te = T[s][T[s].fold < 0.3]
        ce = C[s][C[s].fold < 0.3].copy()
        if te.empty:
            continue
        X, tt, ee, yy, nn_, idx = tensors(te, cols, mu, sd, dev)
        net.eval()
        with torch.no_grad():
            t0, w = net(X, tt, ee, nn_)
        ce["ok"] = ce["delta_t"].abs() < PASS_PS
        n_ev = ce[EVT].drop_duplicates().shape[0]
        r = {"end-to-end t0": float(((t0 - yy).abs() < PASS_PS).float().sum()) / n_ev * 100}
        for lab, col in (("TRKPTZ", "trkptz_score"), ("WAVeS", "waves_score")):
            pick = ce.loc[ce.groupby(EVT, sort=False)[col].idxmax()]
            r[lab] = 100.0 * pick["ok"].sum() / n_ev
        r["oracle over clusters"] = 100.0 * ce.groupby(EVT, sort=False)["ok"].max().sum() / n_ev
        res[s] = r
    log(f"  {'method':26s}" + "".join(f"{n:>10s}" for n in names))
    log("  " + "-" * (26 + 10 * len(names)))
    for k in ("TRKPTZ", "WAVeS", "end-to-end t0", "oracle over clusters"):
        log(f"  {k:26s}" + "".join(f"{res[n][k]:>9.1f}%" for n in names if n in res))
    if args.out:
        os.makedirs(os.path.dirname(os.path.abspath(args.out)) or ".", exist_ok=True)
        json.dump({"args": vars(args), "results": res, "val": best}, open(args.out, "w"), indent=2)
        log(f"\nwrote {args.out}")


if __name__ == "__main__":
    main()
