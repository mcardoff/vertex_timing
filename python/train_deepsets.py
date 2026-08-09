#!/usr/bin/env python3
"""
train_deepsets.py — cluster selection via Deep Sets over a cluster's tracks.

Standalone port of training.ipynb §7-§8 so the training can run as a condor job
(or over ssh in a screen) instead of a live notebook. Same model, same losses,
same event-disjoint validation; the notebook stays the place to look at numbers.

  cluster_score = rho( POOL_tracks phi(track) )

trained directly on WHICH CLUSTER TO PICK -- a multi-positive listwise softmax
over the event's clusters, where every cluster inside PASS_PS counts as correct,
so the loss is the physics metric rather than a proxy for it.

Optionally adds an auxiliary per-track head:

  loss = listwise_CE(selection) + lam * BCE(phi_aux(track), truth_is_hs)

lam is SCANNED (default includes 0) rather than set. Adding information has
repeatedly hurt in this study, and truth_is_hs is a demonstrably imperfect
target: a *perfect* per-track tagger under pT pooling caps at 80.2% against a
92.2% oracle. lam=0 reproduces the plain model exactly, so it is the internal
control and the worst case is "0 wins".

Writes <out>/results.json (all metrics, per lambda, per sample) and
<out>/best_model.pt. Nothing is selected on the test set: the epoch and the
lambda are both chosen on a validation slice carved out of TRAIN.

Examples
--------
  # everything under a directory, standard selection
  python train_deepsets.py --input-dir /data/mcardiff/training --out runs/tight

  # the loosened VBS export, quick smoke run
  python train_deepsets.py --input-dir /data/mcardiff/training --selection loose \\
      --epochs 2 --train-events 3000 --test-events 1500 --lambdas 0 --out /tmp/smoke
"""
import argparse, glob, json, os, sys, time, copy

import numpy as np
import pandas as pd
import uproot
import torch
import torch.nn as nn
import torch.nn.functional as F

PASS_PS = 60.0                                   # PASS_SIGMA in clustering_constants.h
EVT     = ["sample_id", "event_num"]
KEY     = ["sample_id", "event_num", "cluster_idx"]
SAMPLES = ["vbf", "zjets", "dijet"]
SAMPLE_NAME = {0.0: "vbf", 1.0: "zjets", 2.0: "dijet", 3.0: "local"}

TRACK_LABEL    = "truth_is_hs"
TRACK_FEATURES = ["pt", "eta", "theta", "z0", "d0", "qOverP",
                  "sigma_z0", "sigma_d0", "sigma_qOverP",
                  "time", "timeRes", "time_valid", "quality", "nhgtd_hits",
                  "z0_pull_pv", "t_pull_cluster",
                  "dr_nearest_fwdjet", "pt_nearest_fwdjet", "is_ghost_of_nearest",
                  "is_lepton", "cluster_time", "cluster_delta_z"]
# Truth can never be an input. TRACK_FEATURES is an allowlist, but assert anyway --
# nHGTDPrimaryHits is truth (Athena's per-hit m_isprime) and was worth ~14 points on
# Z+jets the one time it leaked in.
_bad = [f for f in TRACK_FEATURES if f.startswith("truth_")
        or f in ("mean_nhgtd_primary", "nhgtd_primary")]
assert not _bad, f"TRUTH IN TRACK FEATURES: {_bad}"

CLUSTER_COLS = KEY + ["delta_t", "trkptz_score", "waves_score"]


def log(msg):
    print(msg, flush=True)                       # condor logs are read while running


# ------------------------------------------------------------------ data ----
def find_files(input_dir, selection):
    """<s>_deta0p0_training.root when loose (falling back per sample), else <s>_training.root."""
    tag, out = "deta0p0", {}
    for s in SAMPLES:
        cand = []
        if selection == "loose":
            cand += [f"{s}_{tag}_training.root"]
        cand += [f"{s}_training.root"]
        for base in cand:
            hit = (glob.glob(os.path.join(input_dir, base))
                   or glob.glob(os.path.join(input_dir, s, base)))
            if hit:
                out[s] = os.path.abspath(hit[0])
                log(f"  {s:6s} [{'loose ' if tag in base else 'tight '}] -> {out[s]}")
                break
    if not out:
        hit = glob.glob(os.path.join(input_dir, "training.root"))
        if hit:
            out["local"] = os.path.abspath(hit[0])
            log(f"  {'local':6s} [tight ] -> {out['local']}")
    return out


def load_clusters(files):
    parts = []
    for name, path in files.items():
        d = uproot.open(path)["clusters"].arrays(CLUSTER_COLS, library="pd")
        log(f"  {name:6s} {len(d):>9,} clusters  {d.groupby(EVT).ngroups:>8,} events")
        parts.append(d)
    df = pd.concat(parts, ignore_index=True)
    # (sample_id, event_num, cluster_idx) uniqueness is load-bearing: it is the event
    # key for the split and the listwise groups. It breaks if a tight and a loose
    # export of the SAME sample are loaded together -- loose is a strict superset with
    # identical event_num, so nothing else would complain.
    dup = int(df.duplicated(KEY).sum())
    if dup:
        sys.exit(f"FATAL: {dup:,} duplicate {KEY} rows -- tight and loose exports of "
                 f"the same sample loaded together?")
    df["abs_dt"]   = df["delta_t"].abs()
    df["within60"] = (df["abs_dt"] < PASS_PS).astype(np.float32)
    return df


def core_fraction(frame, col, higher_is_better=True):
    s = frame[col].fillna(-np.inf if higher_is_better else np.inf)
    g = s.groupby([frame[c] for c in EVT], sort=False)
    idx = g.idxmax() if higher_is_better else g.idxmin()
    return 100.0 * (frame.loc[idx.to_numpy(), "abs_dt"] < PASS_PS).mean()


# ------------------------------------------------------------------ model ---
class DeepSets(nn.Module):
    """pool='sum'  : pooled = SUM_i phi(x_i)
       pool='gate' : pooled = SUM_i sigmoid(a(x_i)) phi(x_i) -- an UNNORMALISED
                     attention. Softmax attention would normalise the weights to
                     sum to 1 per cluster, making a 20-track and a 2-track cluster
                     identical; Sum pT is one of the strongest signals here, so
                     destroying magnitude is the wrong inductive bias."""

    def __init__(self, d_in, d_h=96, d_e=64, pool="sum", aux=False):
        super().__init__()
        self.pool, self.use_aux = pool, aux
        self.phi = nn.Sequential(nn.Linear(d_in, d_h), nn.ReLU(),
                                 nn.Linear(d_h, d_h), nn.ReLU(), nn.Linear(d_h, d_e))
        if pool == "gate":
            self.att = nn.Sequential(nn.Linear(d_in, d_h), nn.ReLU(), nn.Linear(d_h, 1))
        self.rho = nn.Sequential(nn.Linear(d_e, d_h), nn.ReLU(), nn.Linear(d_h, 1))
        if aux:
            # reads phi's OWN embedding, so its gradient lands on the shared
            # per-track representation -- that is the entire point of the head
            self.aux = nn.Linear(d_e, 1)

    def forward(self, x, cidx, ncl, want_aux=False):
        h = self.phi(x)
        a = self.aux(h).squeeze(-1) if (want_aux and self.use_aux) else None
        if self.pool == "gate":
            h = h * torch.sigmoid(self.att(x))
        pooled = torch.zeros(ncl, h.shape[1], dtype=h.dtype,
                             device=h.device).index_add_(0, cidx, h)
        sc = self.rho(pooled).squeeze(-1)
        return (sc, a) if want_aux else sc


def listwise_ce(scores, eidx, pos, nev):
    """-log( softmax mass on ACCEPTABLE clusters ), averaged over events.
    Multi-positive: any cluster inside the window is a correct pick, which IS the
    physics metric rather than a proxy for it."""
    m = torch.full((nev,), -1e30, dtype=scores.dtype,
                   device=scores.device).scatter_reduce(0, eidx, scores,
                                                        reduce="amax", include_self=True)
    e = torch.exp(scores - m[eidx])
    Z = torch.zeros(nev, dtype=scores.dtype, device=scores.device).index_add_(0, eidx, e)
    P = torch.zeros(nev, dtype=scores.dtype, device=scores.device).index_add_(0, eidx, e * pos)
    ok = P > 0
    return -(torch.log(P[ok] + 1e-12) - torch.log(Z[ok] + 1e-12)).mean()


# ------------------------------------------------------------------ prep ----
def stream_tracks(paths, fold, per_sample_train, per_sample_test, step):
    """Per-sample quotas, NOT a global cap. A global cap filled sequentially reads the
    first file only and silently trains on one topology -- that bug appeared three
    separate times in this study before it was caught."""
    keep_tr, keep_te = [], []
    for name, path in paths.items():
        ntr = nte = 0
        for b in uproot.iterate(path, sorted(set(TRACK_FEATURES + KEY + [TRACK_LABEL])),
                                library="pd", step_size=step):
            b = b.assign(fold=fold.reindex(pd.MultiIndex.from_arrays(
                [b["sample_id"], b["event_num"]])).to_numpy())
            for is_tr in (True, False):
                cap = per_sample_train if is_tr else per_sample_test
                got = ntr if is_tr else nte
                if got >= cap:
                    continue
                sel = b[(b.fold < 0.7) if is_tr else (b.fold >= 0.7)]
                if not len(sel):
                    continue
                ev = sel[EVT].drop_duplicates().iloc[:max(0, cap - got)]
                sel = sel.merge(ev, on=EVT, how="inner")
                (keep_tr if is_tr else keep_te).append(sel)
                if is_tr:
                    ntr += len(ev)
                else:
                    nte += len(ev)
            if ntr >= per_sample_train and nte >= per_sample_test:
                break
        log(f"  {name:6s} train {ntr:>7,} ev   test {nte:>7,} ev"
            + ("" if ntr >= per_sample_train else "   (all available)"))
    return (pd.concat(keep_tr, ignore_index=True),
            pd.concat(keep_te, ignore_index=True))


def make_tensors(frame, mu, sd, na_cols, labels, dev):
    """Sorted by KEY so tracks -> clusters -> events are contiguous runs."""
    f = frame.sort_values(KEY, kind="mergesort").reset_index(drop=True)
    X = ((f[TRACK_FEATURES] - mu) / sd)
    M = f[na_cols].isna().astype(np.float32).to_numpy() if na_cols else None
    # NaN means "no valid time" and 0.0 is a legitimate time, so standardise, fill with
    # 0 (= the mean), and hand the net an explicit missing-indicator -- otherwise it
    # cannot distinguish absent from average.
    X = np.nan_to_num(X.to_numpy(dtype=np.float32), nan=0.0, posinf=0.0, neginf=0.0)
    if M is not None:
        X = np.hstack([X, M])
    cl_id, cl_key = pd.factorize(pd.MultiIndex.from_frame(f[KEY]), sort=False)
    cl = pd.DataFrame(list(cl_key), columns=KEY)
    ev_id, _ = pd.factorize(pd.MultiIndex.from_frame(cl[EVT]), sort=False)
    m = cl.merge(labels, on=KEY, how="left")
    return dict(X=torch.from_numpy(X).to(dev),
                c=torch.from_numpy(cl_id.astype(np.int64)).to(dev),
                e=torch.from_numpy(ev_id.astype(np.int64)).to(dev),
                y=torch.from_numpy(m["within60"].fillna(0).to_numpy(np.float32)).to(dev),
                a=torch.from_numpy(f[TRACK_LABEL].to_numpy(np.float32)).to(dev),
                m=m, ncl=len(m), nev=int(ev_id.max()) + 1)


def spans(d):
    c, e = d["c"].cpu().numpy(), d["e"].cpu().numpy()
    return (np.searchsorted(c, np.arange(d["ncl"])),
            np.searchsorted(c, np.arange(d["ncl"]), side="right"),
            np.searchsorted(e, np.arange(d["nev"])),
            np.searchsorted(e, np.arange(d["nev"]), side="right"))


def evaluate(net, d):
    net.eval()
    with torch.no_grad():
        s = net(d["X"], d["c"], d["ncl"]).cpu().numpy()
    net.train()
    out = d["m"].copy()
    out["ds"] = s
    per = {}
    for sid, sub in out.groupby("sample_id"):
        pick = sub.loc[sub.groupby(EVT, sort=False)["ds"].idxmax()]
        per[SAMPLE_NAME.get(sid, str(sid))] = float(100 * pick["within60"].mean())
    pick = out.loc[out.groupby(EVT, sort=False)["ds"].idxmax()]
    return per, float(100 * pick["within60"].mean()), out


def train_one(pool, lam, FIT, VAL, sp, args, dev):
    torch.manual_seed(args.seed)
    net = DeepSets(FIT["X"].shape[1], args.hidden, args.embed, pool, aux=lam > 0).to(dev)
    opt = torch.optim.Adam(net.parameters(), lr=args.lr)
    sched = torch.optim.lr_scheduler.CosineAnnealingLR(opt, T_max=args.epochs,
                                                       eta_min=args.lr / 20)
    CS, CE, ES, EE = sp
    tag = f"{pool} lam={lam:g}"
    log(f"\n[{tag}] parameters: {sum(p.numel() for p in net.parameters()):,}")
    rng = np.random.default_rng(args.seed)
    order = np.arange(FIT["nev"])
    best, best_state, best_ep, hist = -1.0, None, -1, []
    for ep in range(args.epochs):
        rng.shuffle(order)
        tot, nb, t0 = 0.0, 0, time.time()
        for i in range(0, len(order), args.batch_events):
            evs = order[i:i + args.batch_events]
            cl_idx = np.concatenate([np.arange(ES[e], EE[e]) for e in evs])
            tk_idx = np.concatenate([np.arange(CS[c], CE[c]) for c in cl_idx])
            loc_c = np.repeat(np.arange(len(cl_idx)), CE[cl_idx] - CS[cl_idx])
            loc_e = np.repeat(np.arange(len(evs)), EE[evs] - ES[evs])
            tk_t = torch.from_numpy(tk_idx).to(dev)
            cl_t = torch.from_numpy(cl_idx).to(dev)
            # want_aux=True ALWAYS: forward() returns a bare tensor when it is False,
            # so `sc, aux = <tensor>` would raise and the lam=0 control -- the baseline
            # the whole scan is measured against -- would never run.
            sc, aux = net(FIT["X"][tk_t], torch.from_numpy(loc_c).to(dev),
                          len(cl_idx), want_aux=True)
            loss = listwise_ce(sc, torch.from_numpy(loc_e).to(dev),
                               FIT["y"][cl_t], len(evs))
            if lam > 0:
                loss = loss + lam * F.binary_cross_entropy_with_logits(aux, FIT["a"][tk_t])
            opt.zero_grad()
            loss.backward()
            opt.step()
            tot += loss.detach().item()
            nb += 1
        sched.step()
        vper, _, _ = evaluate(net, VAL)
        macro = float(np.mean(list(vper.values())))   # macro so vbf cannot dominate
        star = ""
        if macro > best:
            best, best_state, best_ep, star = macro, copy.deepcopy(net.state_dict()), ep + 1, "  *"
        hist.append({"epoch": ep + 1, "loss": tot / max(nb, 1), "val": vper, "macro": macro})
        log(f"  [{tag}] epoch {ep+1:>2}/{args.epochs}  loss {tot/max(nb,1):.4f}  "
            f"[{time.time()-t0:.0f}s]  VAL "
            + "  ".join(f"{k} {v:.1f}%" for k, v in vper.items())
            + f"  macro {macro:.2f}{star}")
    net.load_state_dict(best_state)
    log(f"  [{tag}] restored epoch {best_ep} (VAL macro {best:.2f})")
    return net, best, best_ep, hist


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--input-dir", required=True, help="directory holding *_training.root")
    p.add_argument("--out", required=True, help="output directory for results.json / best_model.pt")
    p.add_argument("--selection", choices=["tight", "loose"], default="loose")
    p.add_argument("--pools", default="sum", help="comma-separated: sum,gate")
    p.add_argument("--lambdas", default="0,0.1,0.3,1.0",
                   help="auxiliary-head weights to scan; KEEP 0 as the control")
    p.add_argument("--train-events", type=int, default=150_000, help="total, split per sample")
    p.add_argument("--test-events", type=int, default=100_000)
    p.add_argument("--epochs", type=int, default=15)
    p.add_argument("--batch-events", type=int, default=512)
    p.add_argument("--lr", type=float, default=2e-3)
    p.add_argument("--hidden", type=int, default=96)
    p.add_argument("--embed", type=int, default=64)
    p.add_argument("--val-lo", type=float, default=0.60,
                   help="train-fold events >= this are the validation slice")
    p.add_argument("--step", default="300 MB", help="uproot.iterate batch size")
    p.add_argument("--seed", type=int, default=0)
    args = p.parse_args()

    os.makedirs(args.out, exist_ok=True)
    dev = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    log(f"torch {torch.__version__}   device {dev}"
        + (f" ({torch.cuda.get_device_name(0)})" if dev.type == "cuda" else ""))
    log(f"selection: {args.selection}")

    files = find_files(args.input_dir, args.selection)
    if not files:
        sys.exit(f"FATAL: no *_training.root under {args.input_dir}")
    df = load_clusters(files)
    log(f"total {len(df):,} clusters / {df.groupby(EVT).ngroups:,} events")

    # Split BY EVENT -- two clusters of one event must never straddle train/test.
    rng = np.random.default_rng(args.seed)
    ev = df[EVT].drop_duplicates().reset_index(drop=True)
    ev["fold"] = rng.random(len(ev))
    fold = ev.set_index(EVT)["fold"]
    df = df.merge(ev, on=EVT, how="left")
    test_df = df[df.fold >= 0.7]

    log("\nreference selectors (test events):")
    ref = {}
    for sid, sub in test_df.groupby("sample_id"):
        nm = SAMPLE_NAME.get(sid, str(sid))
        ref[nm] = {"TRKPTZ": core_fraction(sub, "trkptz_score"),
                   "WAVeS":  core_fraction(sub, "waves_score"),
                   "oracle": core_fraction(sub, "abs_dt", False)}
        log(f"  {nm:6s} TRKPTZ {ref[nm]['TRKPTZ']:5.1f}%   WAVeS {ref[nm]['WAVeS']:5.1f}%"
            f"   oracle {ref[nm]['oracle']:5.1f}%")

    log("\nstreaming tracks:")
    n_s = max(1, len(files))
    TR_T, TE_T = stream_tracks({k: f"{v}:tracks" for k, v in files.items()}, fold,
                               args.train_events // n_s, args.test_events // n_s, args.step)
    if TR_T["sample_id"].nunique() != len(files):
        sys.exit("FATAL: not every sample reached the track training set -- quota logic broke")

    # Model selection needs its own split; picking the epoch on TEST would be peeking.
    FIT_T = TR_T[TR_T.fold < args.val_lo]
    VAL_T = TR_T[TR_T.fold >= args.val_lo]
    log(f"tracks: fit {len(FIT_T):,} rows / {FIT_T.groupby(EVT).ngroups:,} ev   "
        f"val {len(VAL_T):,} / {VAL_T.groupby(EVT).ngroups:,} ev   "
        f"test {len(TE_T):,} / {TE_T.groupby(EVT).ngroups:,} ev")

    mu = FIT_T[TRACK_FEATURES].mean()                 # standardise on FIT only
    sd = FIT_T[TRACK_FEATURES].std().replace(0, 1.0)
    na_cols = [c for c in TRACK_FEATURES if FIT_T[c].isna().any()]
    labels = df[KEY + ["within60"]].drop_duplicates(KEY)
    log(f"{len(TRACK_FEATURES)} features + {len(na_cols)} missing-indicators")

    FIT = make_tensors(FIT_T, mu, sd, na_cols, labels, dev)
    VAL = make_tensors(VAL_T, mu, sd, na_cols, labels, dev)
    TST = make_tensors(TE_T, mu, sd, na_cols, labels, dev)
    sp = spans(FIT)

    results, best_key, best_macro, best_net = {}, None, -1.0, None
    for pool in args.pools.split(","):
        for lam in [float(x) for x in args.lambdas.split(",")]:
            net, vmacro, vep, hist = train_one(pool.strip(), lam, FIT, VAL, sp, args, dev)
            per, pooled, _ = evaluate(net, TST)
            key = f"{pool.strip()} lam={lam:g}"
            results[key] = {"val_macro": vmacro, "best_epoch": vep,
                            "test_per_sample": per, "test_pooled": pooled,
                            "history": hist}
            if vmacro > best_macro:
                best_key, best_macro, best_net = key, vmacro, net

    log("\n" + "=" * 66)
    names = sorted({n for r in results.values() for n in r["test_per_sample"]})
    log(f"{'selector':26s}" + "".join(f"{n:>10s}" for n in names))
    for k, r in results.items():
        log(f"{'Deep Sets [' + k + ']':26s}"
            + "".join(f"{r['test_per_sample'].get(n, float('nan')):9.1f}%" for n in names))
    log(f"{'TRKPTZ (reference)':26s}" + "".join(f"{ref[n]['TRKPTZ']:9.1f}%" for n in names))
    log(f"{'oracle (ceiling)':26s}" + "".join(f"{ref[n]['oracle']:9.1f}%" for n in names))
    # Chosen on VALIDATION, never on test.
    log(f"\nselected on validation: {best_key}  (macro {best_macro:.2f})")
    log("lambda scan (val macro):  "
        + "   ".join(f"{k.split('lam=')[1]}: {r['val_macro']:.2f}" for k, r in results.items()))
    if best_key.endswith("lam=0"):
        log("  -> lam=0 won: the auxiliary target does not improve the representation, "
            "consistent with its measured 80.2% ceiling")
    else:
        log("  -> the auxiliary head helped: truth_is_hs shapes phi usefully")

    with open(os.path.join(args.out, "results.json"), "w") as fh:
        json.dump({"args": vars(args), "reference": ref, "results": results,
                   "selected": best_key, "selected_val_macro": best_macro}, fh, indent=2)
    torch.save({"state_dict": best_net.state_dict(), "key": best_key,
                "mu": mu.to_dict(), "sd": sd.to_dict(), "na_cols": na_cols,
                "features": TRACK_FEATURES, "hidden": args.hidden, "embed": args.embed},
               os.path.join(args.out, "best_model.pt"))
    log(f"\nwrote {args.out}/results.json and {args.out}/best_model.pt")


if __name__ == "__main__":
    main()
