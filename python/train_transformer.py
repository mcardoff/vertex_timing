#!/usr/bin/env python3
"""
train_transformer.py — vertex time by DIRECT REGRESSION, no clustering.

Every earlier approach picks a time-cluster and uses its time, so it is bounded
by what the clustering happened to produce: the "oracle" 92.2% on Z+jets means
7.8% of events contain NO acceptable cluster and are unwinnable by construction,
and S14's 80.2% is an artefact of pT pooling. This model is bounded by neither.
It reads the event's forward tracks directly and predicts t_HS.

  tokens  = the event's tracks (~29, max ~101) + a CLS token carrying
            event-level context (vertex z, jet summaries, multiplicities)
  encoder = TransformerEncoder -- tracks attend to EACH OTHER, which a Deep Sets
            sum cannot express
  head    = MIXTURE DENSITY, K components (weight, mu, sigma)

WHY A MIXTURE AND NOT A PLAIN REGRESSION. The posterior over t_HS is genuinely
multimodal: one mode per candidate vertex. A point estimate trained on MSE learns
the MEAN of those modes, which is wrong for every one of them -- two candidates at
-80 and +80 ps give a prediction of 0, missing both. It would look like the
architecture failing when it is the output parameterisation. The mixture lets the
model keep the candidates separate and commit to one; taking the argmax-weight
component is "learned clustering + learned selection", end to end and
differentiable. It also yields a per-event uncertainty, which no earlier model
produced. The report prints the mean-prediction score alongside, which makes the
size of this effect visible rather than assumed.

METRIC is |t_pred - t_truth| < PASS_PS, i.e. exactly the core fraction every other
number in this study is quoted in -- but with no clustering ceiling above it.

Optional AUXILIARY supervision, lambda-scanned (never fixed), with the same
noise-aware verdict as train_deepsets.py:
  * truth_is_hs                 -- per-track, always available
  * truth_nearest_fwdjet_is_hs  -- per-track jet identity, used when the exporter
                                   provides it (see the TODO in
                                   util/export_training_data.cxx); silently
                                   skipped when absent.

Examples
--------
  python train_transformer.py --input-dir /data/mcardiff/training --out runs/tf
  python train_transformer.py --input-dir figs/hists --epochs 3 \\
      --train-events 3000 --test-events 1500 --lambdas 0 --out /tmp/tf_smoke
"""
import argparse, glob, json, os, sys, time, copy, statistics

import numpy as np
import pandas as pd
import uproot
import torch
import torch.nn as nn
import torch.nn.functional as F

PASS_PS = 60.0
EVT = ["sample_id", "event_num"]
KEY = ["sample_id", "event_num", "cluster_idx"]
SAMPLES = ["vbf", "zjets", "dijet"]
SAMPLE_NAME = {0.0: "vbf", 1.0: "zjets", 2.0: "dijet", 3.0: "local"}

TRACK_LABEL = "truth_is_hs"
JET_LABEL   = "truth_nearest_fwdjet_is_hs"       # optional; see exporter TODO
TRACK_FEATURES = ["pt", "eta", "theta", "z0", "d0", "qOverP",
                  "sigma_z0", "sigma_d0", "sigma_qOverP",
                  "time", "timeRes", "time_valid", "quality", "nhgtd_hits",
                  "z0_pull_pv", "t_pull_cluster",
                  "dr_nearest_fwdjet", "pt_nearest_fwdjet", "is_ghost_of_nearest",
                  "is_lepton"]
# NB cluster_time / cluster_delta_z are deliberately EXCLUDED: they are outputs of
# the clustering this model exists to do without, so feeding them back would
# reintroduce the very dependence being removed.
EVENT_FEATURES = ["event_sumpt", "n_clusters", "reco_vtx_z", "n_reco_vertices",
                  "n_forward_jets", "lead_jet_pt", "sublead_jet_pt",
                  "lead_jet_abseta", "sublead_jet_abseta", "n_fwd_tracks_reco"]
CLUSTER_COLS = KEY + ["delta_t", "cluster_time", "trkptz_score", "waves_score"]

_bad = [f for f in TRACK_FEATURES + EVENT_FEATURES if f.startswith("truth_")]
assert not _bad, f"TRUTH IN FEATURES: {_bad}"


def log(m):
    print(m, flush=True)


# ------------------------------------------------------------------- data ---
def find_files(input_dir, selection):
    tags, out = ("novbs", "deta0p0"), {}
    for s in SAMPLES:
        cand = ([f"{s}_{t}_training.root" for t in tags] if selection == "loose" else [])
        cand += [f"{s}_training.root"]
        for base in cand:
            hit = (glob.glob(os.path.join(input_dir, base))
                   or glob.glob(os.path.join(input_dir, s, base)))
            if hit:
                out[s] = os.path.abspath(hit[0])
                log(f"  {s:6s} [{next((t for t in tags if t in base), 'tight'):7s}] -> {out[s]}")
                break
    if not out:
        hit = glob.glob(os.path.join(input_dir, "training.root"))
        if hit:
            out["local"] = os.path.abspath(hit[0])
            log(f"  {'local':6s} [tight  ] -> {out['local']}")
    return out


def load_events(files):
    """One row per EVENT: the regression target plus event-level context.

    t_truth is recovered as cluster_time - delta_t, which is identical for every
    cluster of an event (delta_t is defined against the same truth HS time), so
    any cluster reconstructs it. Verified below rather than assumed."""
    parts = []
    for name, path in files.items():
        cols = CLUSTER_COLS + [c for c in EVENT_FEATURES
                               if c in uproot.open(path)["clusters"].keys()]
        d = uproot.open(path)["clusters"].arrays(cols, library="pd")
        log(f"  {name:6s} {len(d):>9,} clusters  {d.groupby(EVT).ngroups:>8,} events")
        parts.append(d)
    df = pd.concat(parts, ignore_index=True)
    dup = int(df.duplicated(KEY).sum())
    if dup:
        sys.exit(f"FATAL: {dup:,} duplicate {KEY} rows -- tight and loose exports "
                 f"of one sample loaded together?")
    df["abs_dt"]  = df["delta_t"].abs()
    df["t_truth"] = df["cluster_time"] - df["delta_t"]

    spread = df.groupby(EVT)["t_truth"].agg(lambda x: x.max() - x.min())
    bad = int((spread > 1e-3).sum())
    if bad:
        sys.exit(f"FATAL: t_truth is not constant within {bad:,} events -- "
                 f"cluster_time - delta_t is not the event's truth HS time")
    log(f"  t_truth consistent across clusters in all {len(spread):,} events")
    return df


def core_fraction(sub, col, hib=True):
    s = sub[col].fillna(-np.inf if hib else np.inf)
    g = s.groupby([sub[c] for c in EVT], sort=False)
    idx = g.idxmax() if hib else g.idxmin()
    return 100.0 * (sub.loc[idx.to_numpy(), "abs_dt"] < PASS_PS).mean()


def stream_tracks(paths, fold, per_train, per_test, step, read_cols):
    keep_tr, keep_te = [], []
    for name, path in paths.items():
        ntr = nte = 0
        for b in uproot.iterate(path, read_cols, library="pd", step_size=step):
            b = b.assign(fold=fold.reindex(pd.MultiIndex.from_arrays(
                [b["sample_id"], b["event_num"]])).to_numpy())
            for is_tr in (True, False):
                cap, got = (per_train, ntr) if is_tr else (per_test, nte)
                if got >= cap:
                    continue
                sel = b[(b.fold < 0.7) if is_tr else (b.fold >= 0.7)]
                if not len(sel):
                    continue
                ev = sel[EVT].drop_duplicates().iloc[:max(0, cap - got)]
                sel = sel.merge(ev, on=EVT, how="inner")
                (keep_tr if is_tr else keep_te).append(sel)
                if is_tr: ntr += len(ev)
                else:     nte += len(ev)
            if ntr >= per_train and nte >= per_test:
                break
        log(f"  {name:6s} train {ntr:>7,} ev   test {nte:>7,} ev"
            + ("" if ntr >= per_train else "   (all available)"))
    return (pd.concat(keep_tr, ignore_index=True),
            pd.concat(keep_te, ignore_index=True))


def build(frame, ev_df, mu, sd, emu, esd, na_cols, dev, aux_cols, ty=(0.0, 1.0)):
    """Pad each event's tracks into [n_events, max_len, n_feat] with a key mask."""
    f = frame.sort_values(EVT + ["track_idx"], kind="mergesort").reset_index(drop=True)
    X = np.nan_to_num(((f[TRACK_FEATURES] - mu) / sd).to_numpy(np.float32),
                      nan=0.0, posinf=0.0, neginf=0.0)
    if na_cols:
        X = np.hstack([X, f[na_cols].isna().astype(np.float32).to_numpy()])
    ev_id, ev_key = pd.factorize(pd.MultiIndex.from_frame(f[EVT]), sort=False)
    nev = int(ev_id.max()) + 1
    counts = np.bincount(ev_id, minlength=nev)
    L = int(counts.max())

    Xp = np.zeros((nev, L, X.shape[1]), np.float32)
    Mp = np.ones((nev, L), bool)                       # True == PAD (torch convention)
    Ap = {c: np.zeros((nev, L), np.float32) for c in aux_cols}
    pos = np.zeros(nev, np.int64)
    for i, e in enumerate(ev_id):
        j = pos[e]; Xp[e, j] = X[i]; Mp[e, j] = False
        for c in aux_cols:
            Ap[c][e, j] = f[c].iat[i]
        pos[e] = j + 1

    ek = pd.DataFrame(list(ev_key), columns=EVT).merge(ev_df, on=EVT, how="left")
    E = np.nan_to_num(((ek[EVENT_FEATURES] - emu) / esd).to_numpy(np.float32),
                      nan=0.0, posinf=0.0, neginf=0.0)
    return dict(X=torch.from_numpy(Xp).to(dev),
                M=torch.from_numpy(Mp).to(dev),
                E=torch.from_numpy(E).to(dev),
                y=torch.from_numpy(ek["t_truth"].to_numpy(np.float32)).to(dev),
                # STANDARDISED target for the loss. Training on raw ps is a trap:
                # t_truth spans hundreds of ps while a fresh head emits mu~0 with
                # sigma=exp(0)=1 ps, so the initial NLL is enormous and the model
                # escapes by inflating sigma and predicting the mean -- it then
                # stalls there and looks like a failed architecture. Predictions
                # are converted back to ps for the metric.
                ys=torch.from_numpy(((ek["t_truth"] - ty[0]) / ty[1])
                                    .to_numpy(np.float32)).to(dev),
                ty=ty,
                A={c: torch.from_numpy(v).to(dev) for c, v in Ap.items()},
                meta=ek, nev=nev, L=L)


# ------------------------------------------------------------------ model ---
class TimeFormer(nn.Module):
    def __init__(self, d_trk, d_evt, d_model=96, nhead=4, layers=3, K=4,
                 aux_names=()):
        super().__init__()
        self.K = K
        self.trk = nn.Linear(d_trk, d_model)
        self.cls = nn.Linear(d_evt, d_model)           # event context AS the CLS token
        enc = nn.TransformerEncoderLayer(d_model, nhead, 4 * d_model, batch_first=True,
                                         dropout=0.0, norm_first=True)
        self.enc = nn.TransformerEncoder(enc, layers)
        self.head = nn.Linear(d_model, 3 * K)          # (logit_w, mu, log_sigma) x K
        self.aux = nn.ModuleDict({c: nn.Linear(d_model, 1) for c in aux_names})

    def forward(self, X, M, E, want_aux=False):
        h = torch.cat([self.cls(E).unsqueeze(1), self.trk(X)], 1)
        m = torch.cat([torch.zeros(M.shape[0], 1, dtype=torch.bool, device=M.device), M], 1)
        h = self.enc(h, src_key_padding_mask=m)
        o = self.head(h[:, 0])                          # read the mixture off CLS
        logit_w, mu, log_sigma = o.split(self.K, dim=-1)
        a = {c: self.aux[c](h[:, 1:]).squeeze(-1) for c in self.aux} if want_aux else {}
        return logit_w, mu, log_sigma.clamp(-6.0, 4.0), a


def mdn_nll(logit_w, mu, log_sigma, y):
    """-log sum_k w_k N(y; mu_k, sigma_k). Keeps candidate vertices as separate
    modes instead of averaging them into a time that fits none."""
    logw = F.log_softmax(logit_w, -1)
    z = (y.unsqueeze(-1) - mu) / log_sigma.exp()
    logN = -0.5 * z ** 2 - log_sigma - 0.5 * float(np.log(2 * np.pi))
    return -(torch.logsumexp(logw + logN, -1)).mean()


def predict(logit_w, mu):
    """argmax-weight component = commit to one candidate. Also return the mixture
    MEAN, which is what a plain regression would learn -- reported so the cost of
    averaging over modes is measured rather than argued about."""
    w = F.softmax(logit_w, -1)
    return mu.gather(1, w.argmax(1, keepdim=True)).squeeze(1), (w * mu).sum(1)


def evaluate(net, d, bs=512):
    net.eval()
    best, mean = [], []
    with torch.no_grad():
        for i in range(0, d["nev"], bs):
            sl = slice(i, i + bs)
            lw, mu, _, _ = net(d["X"][sl], d["M"][sl], d["E"][sl])
            b, m = predict(lw, mu)
            t0_, ts_ = d["ty"]
            b, m = b * ts_ + t0_, m * ts_ + t0_
            best.append(b.cpu().numpy()); mean.append(m.cpu().numpy())
    net.train()
    out = d["meta"].copy()
    out["pred"] = np.concatenate(best)
    out["pred_mean"] = np.concatenate(mean)
    out["res"] = (out["pred"] - out["t_truth"]).abs()
    out["res_mean"] = (out["pred_mean"] - out["t_truth"]).abs()
    per, per_mean = {}, {}
    for sid, sub in out.groupby("sample_id"):
        n = SAMPLE_NAME.get(sid, str(sid))
        per[n] = float(100 * (sub["res"] < PASS_PS).mean())
        per_mean[n] = float(100 * (sub["res_mean"] < PASS_PS).mean())
    return per, per_mean, float(100 * (out["res"] < PASS_PS).mean()), out


def train_one(lam, aux_cols, FIT, VAL, args, dev):
    torch.manual_seed(args.seed)
    net = TimeFormer(FIT["X"].shape[2], FIT["E"].shape[1], args.d_model, args.nhead,
                     args.layers, args.mixtures,
                     aux_names=aux_cols if lam > 0 else ()).to(dev)
    opt = torch.optim.AdamW(net.parameters(), lr=args.lr, weight_decay=1e-2)
    sch = torch.optim.lr_scheduler.CosineAnnealingLR(opt, T_max=args.epochs,
                                                     eta_min=args.lr / 20)
    tag = f"lam={lam:g}"
    log(f"\n[{tag}] parameters: {sum(p.numel() for p in net.parameters()):,}"
        + (f"   aux: {','.join(aux_cols)}" if lam > 0 and aux_cols else ""))
    rng = np.random.default_rng(args.seed)
    order = np.arange(FIT["nev"])
    best, best_state, best_ep, hist = -1.0, None, -1, []
    for ep in range(args.epochs):
        rng.shuffle(order)
        tot, nb, t0 = 0.0, 0, time.time()
        for i in range(0, len(order), args.batch_events):
            b = torch.from_numpy(order[i:i + args.batch_events]).to(dev)
            lw, mu, ls, a = net(FIT["X"][b], FIT["M"][b], FIT["E"][b],
                                want_aux=lam > 0)
            loss = mdn_nll(lw, mu, ls, FIT["ys"][b])
            if lam > 0 and a:
                keep = ~FIT["M"][b]
                for c, logits in a.items():
                    loss = loss + (lam / len(a)) * F.binary_cross_entropy_with_logits(
                        logits[keep], FIT["A"][c][b][keep])
            opt.zero_grad(); loss.backward()
            nn.utils.clip_grad_norm_(net.parameters(), 5.0)
            opt.step()
            tot += loss.detach().item(); nb += 1
        sch.step()
        vper, _, _, _ = evaluate(net, VAL)
        macro = float(np.mean(list(vper.values())))
        star = ""
        if macro > best:
            best, best_state, best_ep, star = macro, copy.deepcopy(net.state_dict()), ep + 1, "  *"
        hist.append(macro)
        log(f"  [{tag}] epoch {ep+1:>2}/{args.epochs}  nll {tot/max(nb,1):.4f}  "
            f"[{time.time()-t0:.0f}s]  VAL "
            + "  ".join(f"{k} {v:.1f}%" for k, v in vper.items())
            + f"  macro {macro:.2f}{star}")
    net.load_state_dict(best_state)
    log(f"  [{tag}] restored epoch {best_ep} (VAL macro {best:.2f})")
    return net, best, hist


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--input-dir", required=True)
    p.add_argument("--out", required=True)
    p.add_argument("--selection", choices=["tight", "loose"], default="tight")
    p.add_argument("--lambdas", default="0,0.1,0.3",
                   help="auxiliary weights to scan; KEEP 0 as the control")
    p.add_argument("--train-events", type=int, default=150_000)
    p.add_argument("--test-events", type=int, default=100_000)
    p.add_argument("--epochs", type=int, default=20)
    p.add_argument("--batch-events", type=int, default=256)
    p.add_argument("--lr", type=float, default=3e-4)
    p.add_argument("--d-model", type=int, default=96)
    p.add_argument("--nhead", type=int, default=4)
    p.add_argument("--layers", type=int, default=3)
    p.add_argument("--mixtures", type=int, default=4)
    p.add_argument("--val-lo", type=float, default=0.60)
    p.add_argument("--step", default="300 MB")
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
    df = load_events(files)

    rng = np.random.default_rng(args.seed)
    ev = df[EVT].drop_duplicates().reset_index(drop=True)
    ev["fold"] = rng.random(len(ev))
    fold = ev.set_index(EVT)["fold"]
    df = df.merge(ev, on=EVT, how="left")

    log("\nreference selectors (test events, cluster-based):")
    ref = {}
    for sid, sub in df[df.fold >= 0.7].groupby("sample_id"):
        n = SAMPLE_NAME.get(sid, str(sid))
        ref[n] = {"TRKPTZ": core_fraction(sub, "trkptz_score"),
                  "WAVeS": core_fraction(sub, "waves_score"),
                  "cluster_oracle": core_fraction(sub, "abs_dt", False)}
        log(f"  {n:6s} TRKPTZ {ref[n]['TRKPTZ']:5.1f}%   WAVeS {ref[n]['WAVeS']:5.1f}%"
            f"   cluster-oracle {ref[n]['cluster_oracle']:5.1f}%  <- a CLUSTERING ceiling,"
            f" not a ceiling on this model")

    # which auxiliary targets actually exist in this export?
    avail = set(uproot.open(list(files.values())[0])["tracks"].keys())
    aux_cols = [c for c in (TRACK_LABEL, JET_LABEL) if c in avail]
    missing = [c for c in (TRACK_LABEL, JET_LABEL) if c not in avail]
    if missing:
        log(f"\nauxiliary targets absent from this export (skipped): {missing}"
            f"   -- see the TODO in util/export_training_data.cxx")
    log(f"auxiliary targets in use: {aux_cols or 'none'}")

    read = sorted(set(TRACK_FEATURES + EVT + ["track_idx"] + aux_cols))
    log("\nstreaming tracks:")
    n_s = max(1, len(files))
    TR, TE = stream_tracks({k: f"{v}:tracks" for k, v in files.items()}, fold,
                           args.train_events // n_s, args.test_events // n_s,
                           args.step, read)
    if TR["sample_id"].nunique() != len(files):
        sys.exit("FATAL: not every sample reached the training set -- quota logic broke")

    FIT_T, VAL_T = TR[TR.fold < args.val_lo], TR[TR.fold >= args.val_lo]
    log(f"tracks: fit {len(FIT_T):,} / {FIT_T.groupby(EVT).ngroups:,} ev   "
        f"val {len(VAL_T):,} / {VAL_T.groupby(EVT).ngroups:,} ev   "
        f"test {len(TE):,} / {TE.groupby(EVT).ngroups:,} ev")

    mu = FIT_T[TRACK_FEATURES].mean(); sd = FIT_T[TRACK_FEATURES].std().replace(0, 1.0)
    na_cols = [c for c in TRACK_FEATURES if FIT_T[c].isna().any()]
    ev_df = df.drop_duplicates(EVT)[EVT + ["t_truth"]
                                    + [c for c in EVENT_FEATURES if c in df.columns]]
    fit_ev = ev_df.merge(FIT_T[EVT].drop_duplicates(), on=EVT)
    emu = fit_ev[EVENT_FEATURES].mean(); esd = fit_ev[EVENT_FEATURES].std().replace(0, 1.0)
    log(f"{len(TRACK_FEATURES)} track features + {len(na_cols)} missing-indicators, "
        f"{len(EVENT_FEATURES)} event features")

    fit_t = ev_df.merge(FIT_T[EVT].drop_duplicates(), on=EVT)["t_truth"]
    TY = (float(fit_t.mean()), float(fit_t.std()) or 1.0)
    log(f"target t_truth: mean {TY[0]:.1f} ps  sd {TY[1]:.1f} ps  "
        f"(standardised for the loss, converted back for the metric)")
    FIT = build(FIT_T, ev_df, mu, sd, emu, esd, na_cols, dev, aux_cols, TY)
    VAL = build(VAL_T, ev_df, mu, sd, emu, esd, na_cols, dev, aux_cols, TY)
    TST = build(TE,    ev_df, mu, sd, emu, esd, na_cols, dev, aux_cols, TY)
    log(f"padded sequences: fit {tuple(FIT['X'].shape)}  (max tracks/event {FIT['L']})")

    results, nets = {}, {}
    for lam in [float(x) for x in args.lambdas.split(",")]:
        net, vmacro, hist = train_one(lam, aux_cols, FIT, VAL, args, dev)
        per, per_mean, pooled, out = evaluate(net, TST)
        k = f"lam={lam:g}"
        nets[k] = net
        results[k] = {"val_macro": vmacro, "history": hist, "test_per_sample": per,
                      "test_per_sample_meanpred": per_mean, "test_pooled": pooled,
                      "median_residual_ps": float(out["res"].median())}

    names = sorted({n for r in results.values() for n in r["test_per_sample"]})
    log("\n" + "=" * 72)
    log(f"{'selector':34s}" + "".join(f"{n:>10s}" for n in names))
    log(f"{'TRKPTZ (cluster-based)':34s}" + "".join(f"{ref[n]['TRKPTZ']:9.1f}%" for n in names))
    for k, r in results.items():
        log(f"{'Transformer+MDN [' + k + ']':34s}"
            + "".join(f"{r['test_per_sample'].get(n, float('nan')):9.1f}%" for n in names))
    log(f"\n{'(mixture MEAN, i.e. what a plain':34s}")
    for k, r in results.items():
        log(f"{'  regression would learn) [' + k + ']':34s}"
            + "".join(f"{r['test_per_sample_meanpred'].get(n, float('nan')):9.1f}%" for n in names))
    log(f"\n{'cluster oracle (NOT our ceiling)':34s}"
        + "".join(f"{ref[n]['cluster_oracle']:9.1f}%" for n in names))

    # lambda chosen on validation, and only when it clears the run's own noise
    log("\nlambda scan (val macro):  "
        + "   ".join(f"{k.split('lam=')[1]}: {r['val_macro']:.3f}" for k, r in results.items()))
    noise = max(statistics.pstdev(r["history"][-8:]) for r in results.values()
                if len(r["history"]) >= 3)
    zero = next((k for k in results if k.endswith("lam=0")), None)
    best_key = max(results, key=lambda k: results[k]["val_macro"])
    if zero is not None:
        nz = {k: v for k, v in results.items() if k != zero}
        if nz:
            top = max(nz, key=lambda k: nz[k]["val_macro"])
            margin = nz[top]["val_macro"] - results[zero]["val_macro"]
            log(f"  epoch-to-epoch noise (sd, last 8): {noise:.3f}   "
                f"best non-zero minus lam=0: {margin:+.3f}")
            if margin < 2 * noise:
                best_key = zero
                log("  -> WITHIN NOISE (margin < 2 sd). Selecting lam=0 on parsimony.")
            else:
                log(f"  -> {top} clears the noise floor: the auxiliary labels help")
    log(f"\nselected: {best_key}  (val macro {results[best_key]['val_macro']:.3f})")

    with open(os.path.join(args.out, "results.json"), "w") as fh:
        json.dump({"args": vars(args), "reference": ref, "results": results,
                   "selected": best_key, "aux_used": aux_cols}, fh, indent=2)
    torch.save({"state_dict": nets[best_key].state_dict(), "key": best_key,
                "track_features": TRACK_FEATURES, "event_features": EVENT_FEATURES,
                "na_cols": na_cols, "mu": mu.to_dict(), "sd": sd.to_dict(),
                "emu": emu.to_dict(), "esd": esd.to_dict(),
                "t_mu": TY[0], "t_sd": TY[1],
                "d_model": args.d_model, "nhead": args.nhead, "layers": args.layers,
                "mixtures": args.mixtures},
               os.path.join(args.out, "best_model.pt"))
    log(f"\nwrote {args.out}/results.json and {args.out}/best_model.pt")


if __name__ == "__main__":
    main()
