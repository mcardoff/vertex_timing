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

THE ARCHITECTURE IS FROZEN. Capacity, pooling form, the auxiliary head and the
feature sets were all scanned and all came back null or negative; the only thing
that ever moved the number was the training objective above. The defaults in
parse_args() ARE that frozen configuration -- see the block there for the
evidence, and reopen a lever only with a concrete reason from a study.

An auxiliary per-track head is still available:

  loss = listwise_CE(selection) + lam * BCE(phi_aux(track), truth_is_hs)

but --lambdas now defaults to "0" (head OFF) rather than to a scan. The scan is
what a previous run mistook for a win, on a +0.006 margin against a 0.05
within-run sd. truth_is_hs is also a demonstrably imperfect target: a *perfect*
per-track tagger under pT pooling caps at 80.2% against a 92.2% oracle.

Two study knobs that do NOT change the model:
  --train-samples   fit on a subset, evaluate on all -> topology transfer matrix
  --sample-frac     scale one sample's train quota   -> learning curve

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
# The event key gained file_idx when the exporter learned --file-shard. Each shard
# builds its own TChain and restarts event_num at 0, so (sample_id, event_num)
# stops identifying an event the moment shards are hadd'd together -- it silently
# merges unrelated events into one listwise group, and every row still looks valid.
# file_idx is the file's index in the sample's FULL sorted list, so the key is
# invariant under shard count. See util/export_training_data.cxx.
EVT     = ["sample_id", "file_idx", "event_num"]
KEY     = EVT + ["cluster_idx"]

# TRAINING samples: mu = 200 only. The mu=0 samples (vbf_mu0, zeejets_mu0,
# ttbar_mu0) are deliberately NOT here. results/baseline_deta0p0_mjj500p0.md
# measures every method at ~99.9% core fraction at mu=0 on every process: one
# vertex means one cluster and a trivial label, so they contribute essentially no
# gradient while diluting the mix with free wins. They are floor/control
# MEASUREMENTS, not training data -- which is where this problem differs from the
# HGTD track-vertex-association BDT, which does train on mu=0 usefully.
SAMPLES = ["vbf", "zjets", "dijet", "ttbar", "zeejets"]

# Must match MyUtl::sampleId() in src/sample_config.h, which is append-only.
# Renumbering there silently relabels every dataset already on disk.
SAMPLE_NAME = {-1.0: "local", 0.0: "vbf", 1.0: "zjets", 2.0: "dijet",
               3.0: "zeejets", 4.0: "vbf_mu0", 5.0: "zeejets_mu0",
               6.0: "ttbar_mu0", 7.0: "ttbar"}

TRACK_LABEL    = "truth_is_hs"
TRACK_FEATURES = ["pt", "eta", "theta", "z0", "d0", "qOverP",
                  "sigma_z0", "sigma_d0", "sigma_qOverP",
                  "time", "timeRes", "time_valid", "quality", "nhgtd_hits",
                  "z0_pull_pv", "t_pull_cluster",
                  "dr_nearest_fwdjet", "pt_nearest_fwdjet", "is_ghost_of_nearest",
                  "is_lepton", "cluster_time", "cluster_delta_z",
                  # Added after the vertex-relational study. RecoVtx_z is a full
                  # array and only index 0 was ever read; these recover the
                  # pileup vertices ITk already reconstructed. Conditioned on
                  # |delta_z| the PU-vertex distance still lifts discrimination
                  # 1.78-2.34x in the displaced bins -- which is exactly where
                  # the model picks pileup clusters. closest_vtx_is_pv is NOT
                  # here: conditioned on delta_z its lift collapses to ~1.2 and
                  # inverts at large |delta_z|, so it is a restatement of
                  # delta_z rather than new information.
                  "dz_to_nearest_pu_vtx_trk", "closer_to_pu_than_pv",
                  "cluster_dz_to_nearest_pu_vtx", "cluster_pv_pu_dz_ratio",
                  # Jet association independent of the forward band, and
                  # leave-one-out time stability (studies E and in/out-of-jet).
                  "is_in_any_jet", "dr_nearest_anyjet", "loo_pull"]
# Truth can never be an input. TRACK_FEATURES is an allowlist, but assert anyway --
# nHGTDPrimaryHits is truth (Athena's per-hit m_isprime) and was worth ~14 points on
# Z+jets the one time it leaked in.
_bad = [f for f in TRACK_FEATURES if f.startswith("truth_")
        or f in ("mean_nhgtd_primary", "nhgtd_primary")]
assert not _bad, f"TRUTH IN TRACK FEATURES: {_bad}"

CLUSTER_COLS = KEY + ["delta_t", "trkptz_score", "waves_score"]


def log(msg):
    print(msg, flush=True)                       # condor logs are read while running


def read_tree(path, tree, cols):
    """Read `cols`, tolerating exports that predate file_idx joining the event key.

    A pre-fix file is always a single unsharded export, where event_num was the
    chain-global entry number and so already unique on its own. Substituting a
    constant 0 therefore makes the three-column key exactly equivalent to the old
    two-column one FOR THAT FILE, and keeps the duplicate guard meaningful. Any
    other missing column is fatal -- silently dropping a feature would change what
    the model trains on without changing anything visible."""
    t = uproot.open(path)[tree]
    have = set(t.keys())
    missing = [c for c in cols if c not in have]
    if [c for c in missing if c != "file_idx"]:
        sys.exit(f"FATAL: {path}:{tree} missing {[c for c in missing if c != 'file_idx']}")
    d = t.arrays([c for c in cols if c in have], library="pd")
    if "file_idx" in missing:
        log(f"    note: {os.path.basename(path)} predates file_idx; assuming unsharded")
        d["file_idx"] = np.float32(0)
    return d


# ------------------------------------------------------------------ data ----
# Which export a --selection maps to, most-preferred first. The tag is part of
# the exported FILENAME (see MyUtl::resolveSelection), so the selection a run
# trained on is recoverable from the file it read rather than from memory.
#
#   canonical  --vbs-mjj=500     the analysis definition every baseline in
#                                results/baseline_deta0p0_mjj500p0.md is quoted at
#   loose      --vbs-deta=-1     VBS candidate-PAIR requirement dropped entirely;
#                                the wider population, for the data-limitation study
#
# "deta0p0" and the untagged name are LEGACY fallbacks: deta0p0 was an earlier
# loosening attempt that only dropped the |Deta| magnitude and left the pair
# requirement standing (worth +4% on Z+jets, not the several-fold intended), and
# an untagged file predates the selection being tagged at all. Both are kept so
# old exports still load, and both are labelled as such in the log -- mixing one
# into a canonical comparison is the mistake this labelling exists to prevent.
SELECTION_TAGS = {
    "canonical": ("mjj500p0",),
    "loose":     ("novbs",),
}
LEGACY_TAGS = ("deta0p0",)


def find_files(input_dir, selection):
    """Locate one export per sample for the requested selection."""
    want = SELECTION_TAGS[selection]
    out, stale = {}, []
    for s in SAMPLES:
        cand  = [(f"{s}_{t}_training.root", t, False) for t in want]
        cand += [(f"{s}_{t}_training.root", t, True) for t in LEGACY_TAGS]
        cand += [(f"{s}_training.root", "untagged", True)]
        for base, which, is_legacy in cand:
            hit = (glob.glob(os.path.join(input_dir, base))
                   or glob.glob(os.path.join(input_dir, s, base)))
            if hit:
                out[s] = os.path.abspath(hit[0])
                if is_legacy:
                    stale.append(f"{s} ({which})")
                log(f"  {s:12s} [{which:9s}{'!' if is_legacy else ' '}] -> {out[s]}")
                break
    if not out:
        hit = glob.glob(os.path.join(input_dir, "training.root"))
        if hit:
            out["local"] = os.path.abspath(hit[0])
            log(f"  {'local':12s} [{'untagged':9s}!] -> {out['local']}")
            stale.append("local (untagged)")
    if stale:
        # Loud, but not fatal: a legacy file is usable, it just is not the
        # selection asked for, so any number from it is not comparable with the
        # baseline. Silently mixing fiducial regions is the failure this catches.
        log(f"\n  !! {len(stale)} sample(s) fell back to a NON-'{selection}' export: "
            + ", ".join(stale))
        log("     Their events are from a different fiducial region. Re-export "
            "before quoting these against results/baseline_deta0p0_mjj500p0.md.")
    return out


def load_clusters(files):
    parts = []
    for name, path in files.items():
        d = read_tree(path, "clusters", CLUSTER_COLS)
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

def iter_tracks(path, cols, step):
    """uproot.iterate over the tracks tree, with the same pre-file_idx tolerance
    as read_tree -- otherwise an old export would load its clusters fine and then
    fail here, which is a worse failure than not supporting it at all."""
    have = set(uproot.open(path)["tracks"].keys())
    want = [c for c in cols if c in have]
    # The ":tracks" object path is appended HERE, not by the caller: a bare
    # filename is now ambiguous to uproot.iterate, since the export carries three
    # TTrees (clusters, tracks, jets) rather than two.
    for b in uproot.iterate({path: "tracks"}, want, library="pd", step_size=step):
        if "file_idx" not in have:
            b = b.assign(file_idx=np.float32(0))
        yield b


def stream_tracks(paths, fold, per_sample_train, per_sample_test, step):
    """Per-sample quotas, NOT a global cap. A global cap filled sequentially reads the
    first file only and silently trains on one topology -- that bug appeared three
    separate times in this study before it was caught.

    per_sample_train is a dict {sample: cap}; a cap of 0 streams that sample's TEST
    events but fits on none of it, which is what --train-samples uses to build a
    transfer matrix in a single pass. per_sample_test may be an int (same for all)."""
    keep_tr, keep_te = [], []
    for name, path in paths.items():
        cap_tr = per_sample_train[name] if isinstance(per_sample_train, dict) else per_sample_train
        cap_te = per_sample_test[name] if isinstance(per_sample_test, dict) else per_sample_test
        ntr = nte = 0
        for b in iter_tracks(path, sorted(set(TRACK_FEATURES + KEY + [TRACK_LABEL])),
                             step):
            b = b.assign(fold=fold.reindex(pd.MultiIndex.from_arrays(
                [b[c] for c in EVT])).to_numpy())
            for is_tr in (True, False):
                cap = cap_tr if is_tr else cap_te
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
            if ntr >= cap_tr and nte >= cap_te:
                break
        note = ""
        if cap_tr == 0:
            note = "   (eval only)"
        elif ntr < cap_tr:
            note = "   (all available)"
        log(f"  {name:12s} train {ntr:>7,} ev   test {nte:>7,} ev" + note)
    return (pd.concat(keep_tr, ignore_index=True) if keep_tr else None,
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


def val_loss(net, d):
    """The TRAINING objective evaluated on a held-out slice.

    This is what the epoch is selected on, in place of validation core fraction.
    Core fraction on the validation slice was measured to be unreliable: on VBF
    it collapsed to 17.6% and DEGRADED with training while the same weights
    scored 91.7% on test, and TRKPTZ pushed through the identical tensor scored
    90.0% -- so the slice and its labels are sound and the fault is in the
    argmax-of-model-scores step on that slice specifically. Selecting on it
    restored epoch 1 and cost ~1.3 points of Z+jets test core fraction.

    Loss avoids that failure entirely: it reads the score distribution rather
    than only its argmax, it is the quantity actually being optimised, and it
    moves smoothly. Core fraction is still REPORTED every epoch, because the two
    disagreeing is the signature of the open bug and it should stay visible."""
    net.eval()
    with torch.no_grad():
        sc = net(d["X"], d["c"], d["ncl"])
        l = listwise_ce(sc, d["e"], d["y"], d["nev"]).item()
    net.train()
    return l


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
    # init_seed, NOT seed: --seed fixes the event SPLIT, --init-seed fixes the
    # weights and the batch order. Keeping them separate is what makes an
    # ensemble honest -- members must share one split, or member B has trained on
    # events that are in member A's test fold and the combined number is leaked.
    torch.manual_seed(args.init_seed)
    net = DeepSets(FIT["X"].shape[1], args.hidden, args.embed, pool, aux=lam > 0).to(dev)
    opt = torch.optim.Adam(net.parameters(), lr=args.lr)
    sched = torch.optim.lr_scheduler.CosineAnnealingLR(opt, T_max=args.epochs,
                                                       eta_min=args.lr / 20)
    CS, CE, ES, EE = sp
    tag = f"{pool} lam={lam:g}"
    log(f"\n[{tag}] parameters: {sum(p.numel() for p in net.parameters()):,}")
    rng = np.random.default_rng(args.init_seed)
    order = np.arange(FIT["nev"])
    best_vl, best_state, best_ep, hist = float("inf"), None, -1, []
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
        vl = val_loss(net, VAL)
        # SELECT ON VAL LOSS, not on val macro core fraction -- see val_loss().
        # Lower is better, hence the sign flip against the old `macro > best`.
        star = ""
        if vl < best_vl:
            best_vl, best_state, best_ep, star = vl, copy.deepcopy(net.state_dict()), ep + 1, "  *"
        hist.append({"epoch": ep + 1, "loss": tot / max(nb, 1), "val": vper,
                     "macro": macro, "val_loss": vl})
        log(f"  [{tag}] epoch {ep+1:>2}/{args.epochs}  loss {tot/max(nb,1):.4f}  "
            f"vloss {vl:.4f}  [{time.time()-t0:.0f}s]  VAL "
            + "  ".join(f"{k} {v:.1f}%" for k, v in vper.items())
            + f"  macro {macro:.2f}{star}")
    net.load_state_dict(best_state)
    sel_macro = next(h["macro"] for h in hist if h["epoch"] == best_ep)
    log(f"  [{tag}] restored epoch {best_ep} (val loss {best_vl:.4f}, "
        f"its val macro {sel_macro:.2f})")
    return net, sel_macro, best_ep, hist


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--input-dir", required=True, help="directory holding *_training.root")
    p.add_argument("--out", required=True, help="output directory for results.json / best_model.pt")
    p.add_argument("--selection", choices=["canonical", "loose"], default="canonical")
    # ---- THE ARCHITECTURE IS FROZEN ------------------------------------------
    # These four defaults ARE the frozen configuration -- run 1's 24k-parameter
    # model. The search that closed is documented in results/README.md; briefly,
    # four levers were scanned and every one came back null or negative:
    #   capacity   5.7x parameters (hidden/embed 256/128) scored 0.2 LOWER on all
    #              three samples, and cost 44 s/epoch against 25.
    #   pooling    gated vs sum: +0.1 on zjets, i.e. inside run-to-run noise.
    #   aux head   a per-track truth_is_hs head, lambda scanned: within noise every
    #              time. Hence --lambdas defaults to "0" (head OFF), not a scan --
    #              the scan is what a previous run mistook for a +0.006 win against
    #              a 0.05 within-run sd.
    #   features   cluster features and per-track tagger features, both null.
    # The ONLY change that ever moved the number was the training objective
    # (multi-positive listwise), which is not a flag -- it is the loss.
    #
    # Reopen one of these only with a concrete reason from a study, not on spec.
    # Overriding a default here is how you run that test; leaving them alone is
    # how you reproduce the frozen model.
    p.add_argument("--pools", default="sum", help="comma-separated: sum,gate")
    p.add_argument("--lambdas", default="0",
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
    p.add_argument("--seed", type=int, default=0,
                   help="fixes the EVENT SPLIT. Ensemble members must share this.")
    # Separated from --seed so several models can be trained on ONE split and
    # then combined. The measured seed-to-seed spread is 2.4 points on Z+jets
    # (65.1-67.5 over four seeds) against an infrastructure floor of exactly 0 --
    # that much variance between equivalent models is precisely the condition
    # under which averaging their scores buys real accuracy.
    p.add_argument("--init-seed", type=int, default=None,
                   help="fixes WEIGHTS and batch order; defaults to --seed")
    # torch's intra-op thread count. Left at 0 (torch's own default) unless set.
    #
    # This is a REPRODUCIBILITY knob, not a speed one. torch parallelises its
    # reductions across threads, and the order a reduction is summed in changes
    # its last bits -- so the same seed on the same data gives bitwise different
    # results at different thread counts. A condor run at request_cpus=4 and a
    # rerun on a busy login node diverged from epoch 1 (loss 0.4345 vs 0.4242)
    # for exactly this reason, and the two disagreed by 18 points on one sample's
    # validation number.
    #
    # Set it equal to request_cpus in the submit files. Two runs then differ only
    # by what you meant to change, which is what makes a learning curve or a
    # transfer matrix readable at all.
    p.add_argument("--torch-threads", type=int, default=0,
                   help="torch.set_num_threads(); 0 = leave torch's default")

    # ---- study knobs (do not change the model) -------------------------------
    # --train-samples: fit on a SUBSET, evaluate on everything. This is the
    # topology-transfer test: `--train-samples vbf` then reading the zjets column
    # of the per-sample table answers "does a selector learned on VBF transfer to
    # a sample where the forward jets are mostly pileup?" -- the question WAVeS
    # fails (+0.7 vbf, -1.9 zjets). Excluded samples still stream their TEST
    # events, so every column of the transfer matrix is filled from one run.
    p.add_argument("--train-samples", default="",
                   help="comma-separated subset to FIT on (default: all found). "
                        "Test/eval always covers every sample loaded.")
    # --sample-frac: scale one sample's training quota, for a learning curve.
    # `--sample-frac zjets=0.25` fits on a quarter of the Z+jets training events
    # with everything else unchanged. Repeat at 0.25/0.5/1.0 and the slope says
    # whether Z+jets is data-limited -- it contributed every one of its available
    # train-fold events in the last run, and the canonical selection cut it to
    # ~36.5k, so this is the question that decides whether more MC would help.
    p.add_argument("--sample-frac", default="",
                   help="comma-separated <sample>=<frac>, e.g. zjets=0.25")
    args = p.parse_args()

    if args.init_seed is None:
        args.init_seed = args.seed
    os.makedirs(args.out, exist_ok=True)
    if args.torch_threads > 0:
        torch.set_num_threads(args.torch_threads)
    log(f"torch threads: {torch.get_num_threads()}"
        + ("" if args.torch_threads > 0 else "  (DEFAULT -- runs are not comparable "
                                             "across machines; see --torch-threads)"))
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

    # ---- per-sample training quotas ------------------------------------------
    # --train-samples zeroes the excluded samples' TRAIN quota while leaving their
    # TEST quota intact, so one pass fills every column of a transfer matrix.
    # --sample-frac scales one sample's quota for the learning curve.
    train_on = ([s.strip() for s in args.train_samples.split(",") if s.strip()]
                or list(files))
    unknown = [s for s in train_on if s not in files]
    if unknown:
        sys.exit(f"FATAL: --train-samples names {unknown}; loaded samples are {list(files)}")
    fracs = {}
    for tok in filter(None, (t.strip() for t in args.sample_frac.split(","))):
        k, _, v = tok.partition("=")
        k = k.strip()
        if k not in files:
            sys.exit(f"FATAL: --sample-frac names '{k}'; loaded samples are {list(files)}")
        if k not in train_on:
            sys.exit(f"FATAL: --sample-frac scales '{k}', which --train-samples excludes")
        fracs[k] = float(v)

    # The budget is divided over the samples actually being FIT on, not over every
    # sample loaded -- otherwise excluding a sample would silently shrink the
    # training set instead of just changing its composition, and a transfer run
    # would be confounded by having less data as well as different data.
    base = args.train_events // max(1, len(train_on))
    per_train = {s: (int(round(base * fracs.get(s, 1.0))) if s in train_on else 0)
                 for s in files}
    per_test = args.test_events // max(1, len(files))

    log(f"\ntrain on: {', '.join(train_on)}"
        + (f"   (eval-only: {', '.join(s for s in files if s not in train_on)})"
           if len(train_on) < len(files) else ""))
    if fracs:
        log("sample-frac: " + ", ".join(f"{k}={v:g}" for k, v in fracs.items()))

    log("\nstreaming tracks:")
    TR_T, TE_T = stream_tracks(files, fold, per_train, per_test, args.step)
    if TR_T is None:
        sys.exit("FATAL: --train-samples left nothing to fit on")
    got = TR_T["sample_id"].map(SAMPLE_NAME).nunique()
    if got != len(train_on):
        sys.exit(f"FATAL: {got} of {len(train_on)} requested train samples reached the "
                 "track training set -- quota logic broke")

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

    # Composition of each tensor, per sample. The epoch is selected on VAL macro,
    # so a VAL slice that is not comparable to TEST silently corrupts that choice
    # -- and a per-sample val/test difference is invisible in the macro. The
    # oracle column is the tell: it is a property of the EVENTS alone, so if VAL
    # and TEST disagree there, the slices differ in what is even achievable and
    # no model number between them is comparable.
    # `trk` is a CONTROL: TRKPTZ scored through this tensor's own m/argmax
    # machinery. Its value is known independently -- it is in the reference table
    # printed above -- so if it disagrees there, the tensor is wrong rather than
    # the model. Without a control, a collapsing model number is ambiguous
    # between "the model is bad on these events" and "these rows are misaligned",
    # and those need completely different fixes.
    trk = df[KEY + ["trkptz_score"]].drop_duplicates(KEY)
    log("\ntensor composition (orc = ceiling; trk = TRKPTZ through THIS tensor's "
        "own argmax -- must match the reference table above):")
    for sid in sorted(set(FIT["m"]["sample_id"]) | set(TST["m"]["sample_id"])):
        cells = []
        for d in (FIT, VAL, TST):
            sub = d["m"][d["m"]["sample_id"] == sid]
            nev = sub.groupby(EVT, sort=False).ngroups
            if not nev:
                cells.append(f"{0:>8,} ev  orc   nan  trk   nan"); continue
            orc = 100.0 * sub.groupby(EVT, sort=False)["within60"].max().mean()
            s2 = sub.merge(trk, on=KEY, how="left")
            pick = s2.loc[s2.groupby(EVT, sort=False)["trkptz_score"].idxmax()]
            cells.append(f"{nev:>8,} ev  orc {orc:5.1f}%  trk {100*pick['within60'].mean():5.1f}%")
        log(f"  {SAMPLE_NAME.get(sid, str(sid)):12s}" + "  ".join(cells))

    results, nets, best_key, best_macro, best_net = {}, {}, None, -1.0, None
    for pool in args.pools.split(","):
        for lam in [float(x) for x in args.lambdas.split(",")]:
            net, vmacro, vep, hist = train_one(pool.strip(), lam, FIT, VAL, sp, args, dev)
            per, pooled, _ = evaluate(net, TST)
            key = f"{pool.strip()} lam={lam:g}"
            nets[key] = net
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
    # Chosen on VALIDATION, never on test -- and only when the margin CLEARS THE
    # RUN'S OWN NOISE. An earlier version declared victory for whichever lambda
    # scored highest, and duly reported "the auxiliary head helped" on a margin of
    # +0.006 macro points against an epoch-to-epoch sd of ~0.05 within a single
    # run: a tenth of the noise. The noise floor is measured from the last 8
    # epochs of each run rather than hardcoded, so it self-calibrates.
    log("\nlambda scan (val macro):  "
        + "   ".join(f"{k.split('lam=')[1]}: {r['val_macro']:.3f}" for k, r in results.items()))
    import statistics as _st
    zero = next((k for k in results if k.endswith("lam=0")), None)
    if zero is not None:
        nz = {k: v for k, v in results.items() if k != zero}
        if nz:
            # Computed HERE, not above, for two reasons. It is only meaningful
            # when there is a non-zero lambda to compare against -- and with the
            # frozen default (--lambdas 0) there never is, so computing it
            # unconditionally crashed every single-lambda run on an empty max().
            #
            # The default is +inf, not 0: with fewer than 3 epochs there is no
            # noise estimate, and the honest response to "I cannot measure the
            # noise" is to refuse to declare a winner, which +inf does by making
            # every margin fall inside it. Defaulting to 0 would do the opposite
            # and wave through any margin at all -- reintroducing exactly the
            # +0.006-against-0.05-sd mistake this block exists to prevent.
            noise = max((_st.pstdev([h["macro"] for h in r["history"]][-8:])
                         for r in results.values() if len(r["history"]) >= 3),
                        default=float("inf"))
            top = max(nz, key=lambda k: nz[k]["val_macro"])
            margin = nz[top]["val_macro"] - results[zero]["val_macro"]
            log(f"  epoch-to-epoch noise (sd of last 8): {noise:.3f}   "
                f"best non-zero minus lam=0: {margin:+.3f}")
            if margin < 2 * noise:
                best_key, best_macro, best_net = zero, results[zero]["val_macro"], nets[zero]
                log(f"  -> WITHIN NOISE (margin < 2 sd). Selecting {zero} on parsimony: the "
                    f"auxiliary target does not improve the representation, consistent with "
                    f"its measured 80.2% ceiling.")
            else:
                log(f"  -> {top} clears the noise floor: truth_is_hs shapes phi usefully")
    log(f"\nselected: {best_key}  (val macro {best_macro:.3f})")

    # The learning curve's x-axis must be the events ACTUALLY fitted on, not the
    # requested fraction. Z+jets is the sample the curve is about and it is the
    # one that runs out: it contributed every available train-fold event in the
    # last run, so --sample-frac=1.0 can silently mean "all there was". Plotting
    # against the request would then show a flat top that is a quota artefact
    # rather than saturation -- the exact thing the study is trying to measure.
    fitted = (TR_T.assign(_n=TR_T["sample_id"].map(SAMPLE_NAME))
                  .groupby("_n")[EVT].apply(lambda d: len(d.drop_duplicates()))
                  .to_dict())
    log("\ntrain events actually fitted: "
        + ", ".join(f"{k} {v:,}" for k, v in sorted(fitted.items())))

    with open(os.path.join(args.out, "results.json"), "w") as fh:
        json.dump({"args": vars(args), "reference": ref, "results": results,
                   "train_events_fitted": fitted,
                   "selected": best_key, "selected_val_macro": best_macro}, fh, indent=2)
    torch.save({"state_dict": best_net.state_dict(), "key": best_key,
                "mu": mu.to_dict(), "sd": sd.to_dict(), "na_cols": na_cols,
                "features": TRACK_FEATURES, "hidden": args.hidden, "embed": args.embed},
               os.path.join(args.out, "best_model.pt"))
    log(f"\nwrote {args.out}/results.json and {args.out}/best_model.pt")


if __name__ == "__main__":
    main()
