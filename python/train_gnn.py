#!/usr/bin/env python3
"""
train_gnn.py -- cluster selection AND t0 as a graph over the event's tracks.

Nodes are TRACKS IN THE EVENT, not tracks in a cluster. Edges are all-pairs
within the event, carrying two types:

    same_cluster      i and j landed in the same cluster of the existing partition
    cross_cluster     they did not

so the clustering enters as a PRIOR the network can override, never as the graph
topology. Building edges only inside clusters would disconnect the graph, and the
dominant failure channel is selection -- comparing candidates against each other
(26.5% of principal causes, and results/ml_overview.md finding 1: an acceptable
cluster exists in the collection for the large majority of failures). A model that
cannot see across the partition is blind to exactly that.

WHY A GRAPH AT ALL. This is not a capacity play -- capacity is a measured null
(5.7x parameters scored 0.2 LOWER, see train_deepsets.py). It is that ONE hand-built
round of message passing is already the second most productive thing in the study.
final_method_eval.py:148-165 does an all-pairs join within the event, gates the
neighbour's HS probability by hard boxes (|dt|<60, |dz|<2.0), and sums:

    nb_phzt = SUM_{j != i} ph_j * (|dz_ij| < 2.0) * (|dt_ij| < 60)

That is one GNN layer with a fixed, non-learned edge kernel. Its six columns take
37% of the round-2 tagger's gain (nb_phzt alone 23.76, against 25.98 for the top
feature), and nb_ph60 is the #2 feature in the entire e2e t0 model. What the hand
version discards, and what this exists to test:

  1. the kernels are hard boxes on RAW dt, so per-track resolution is thrown away.
     The clustering next to it uses |dt|/sqrt(s_i^2+s_j^2) < 3. Edges here carry
     the pull, not the picosecond.
  2. one hop, never iterated. HS tracks form a mutually-reinforcing clique.
  3. sum-of-counts is degree-dominated -- which is why the unweighted nb_n60/nb_n30
     rank LAST of the six and the ph-weighted ones rank first.
  4. frozen bootstrap: round 1 -> freeze ph -> round 2 cannot correct round 1.

THE TRAP THIS AVOIDS. A transformer over event tracks is a fully-connected GNN with
learned attention, and that already collapsed on VBF (77%, plateaued by epoch 4;
finding 9). The diagnosis was that it had to LEARN the averaging that cluster
methods get structurally. So the averaging stays structural here: the network emits
per-track weights and the time is

    t_hat_c = SUM_{i in c} w_i t_i / SUM_{i in c} w_i

It never re-partitions. Re-partitioning is separately measured and it loses:
WAVES_RECLUST is -13.5 points on zjets (results/baseline_deta0p0_mjj500p0.md),
split clustering is -3.14 / 25.34% with -2.61 real complementarity
(results/split_clustering.md), misclustering is 3.1% of principal causes, and the
RpT ladder puts clustering at ~5-6% against ~22% for times+selection. The mechanism
is statistics: with ~29 timed tracks and a 60 ps window, a soft weight can push a
contaminant to 0.05 and keep the other 28 contributing, while a partition boundary
throws it out and takes the statistics with it. Every gain this project has banked
works INSIDE a fixed partition -- in-jet re-timing (+0.8), tagger-probability
weights over 1/sigma^2 (+1.72), double-Winsorisation (+0.19).

RESIDUAL FROM THE INCUMBENT. Both heads are parameterised so that at initialisation
the model IS the current algorithm, exactly:

    score_c = tau * log(trkptz_score_c) + f(pool_c)      f's last layer zero-init
    w_i     = (1/sigma_t,i^2) * exp(g(h_i))              g's last layer zero-init

so f = g = 0 at epoch 0 gives argmax = TRKPTZ and t_hat = the inverse-variance mean.
Training is then strictly a refinement of the incumbent, and "did the GNN help?"
becomes readable off the weights against their own initialisation instead of off a
delta that has to clear a +-1.9 split band. It also makes the TRKPTZ control
structural rather than bolted on: --epochs 0 MUST reproduce the reference table.
Ablate with --no-residual.

Written against train_deepsets.py and importing its data path wholesale, so the
event key, the fold split, the per-sample quotas, the multi-positive listwise loss
and the val-LOSS epoch selection are the same code, not a re-implementation.

Examples
--------
  # local VBF draft run
  PYTHONNOUSERSITE=1 ~/.venv-hgtd/bin/python python/train_gnn.py \
      --input-dir /Users/mcard/vertex_timing_backup_20260828/exports \
      --samples vbf --train-events 12000 --test-events 6000 --epochs 8 \
      --out runs/gnn_vbf

  # the control: does the learned edge kernel beat the hand-built hard boxes?
  ... --edge-kernel box --out runs/gnn_vbf_boxkernel
"""
import argparse, copy, json, os, sys, time

import numpy as np
import pandas as pd
import torch
import torch.nn as nn
import torch.nn.functional as F

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import train_deepsets as ds
from train_deepsets import (EVT, KEY, PASS_PS, SAMPLE_NAME, TRACK_LABEL,
                            core_fraction, find_files, listwise_ce, log)

# Rebound in main() to whatever the loaded export actually carries -- see there.
TRACK_FEATURES = list(ds.TRACK_FEATURES)

# cluster_time joins train_deepsets' list so the truth HS time is recoverable as
# cluster_time - delta_t. It is constant within an event by construction (both are
# per-cluster measurements of the same event), which main() asserts.
CLUSTER_COLS = KEY + ["delta_t", "trkptz_score", "waves_score", "truth_purity",
                      "cluster_time"]

# Per-track columns the GRAPH needs, beyond TRACK_FEATURES. These build edges, so
# they are read raw and unstandardised -- an edge pull is only meaningful in
# physical units.
GEOM = ["time", "timeRes", "z0", "sigma_z0", "eta", "phi"]

EDGE_FEATURES = ["dt_pull", "dt_raw", "dz_pull", "dz_raw", "dr", "same_cluster"]


# ------------------------------------------------------------------ graph ---

def event_edges(counts):
    """All-pairs (i != j) node indices for a batch of events with `counts` tracks each.

    Returned indices are LOCAL to the concatenated batch. At 29.4 tracks/event
    (p99 61, max 106) a fully-connected event is ~860 edges, so no kNN
    sparsification, no neighbour sampling, and no edge budget: the complete graph
    is cheaper here than the pandas all-pairs join that produced nb_phzt, which
    needed 60k-event chunking to avoid OOM on ttbar."""
    src, dst, base = [], [], 0
    for n in counts:
        idx = np.arange(base, base + n, dtype=np.int64)
        s = np.repeat(idx, n)
        d = np.tile(idx, n)
        m = s != d
        src.append(s[m]); dst.append(d[m])
        base += n
    if not src:
        return np.zeros(0, np.int64), np.zeros(0, np.int64)
    return np.concatenate(src), np.concatenate(dst)


def edge_attr(G, src, dst):
    """Edge features from the RAW geometry tensor G = [time, timeRes, z0, sigma_z0,
    eta, phi, cluster_id].

    dt_pull is the quantity doIterativeClustering itself cuts on at 3 sigma, and it
    is the one the hand-built nb_* columns throw away by boxing raw |dt| at 60 ps
    while per-track timeRes runs 17.5-35.0 ps. dz_pull likewise. The raw versions
    ride along because the z-association study found raw |dz| in mm beats z0
    significance as a per-track weight (+1.76 vs +1.21 on zjets), so which scale
    wins is an open question and the model gets both."""
    t, tr, z, sz, eta, phi, cid = (G[:, i] for i in range(7))
    dt = t[src] - t[dst]
    dz = z[src] - z[dst]
    dphi = torch.remainder(phi[src] - phi[dst] + np.pi, 2 * np.pi) - np.pi
    deta = eta[src] - eta[dst]
    return torch.stack([
        dt / torch.sqrt(tr[src] ** 2 + tr[dst] ** 2 + 1e-6),
        dt / 50.0,                                        # ps, O(1)-scaled
        dz / torch.sqrt(sz[src] ** 2 + sz[dst] ** 2 + 1e-6),
        dz,                                               # mm
        torch.sqrt(deta ** 2 + dphi ** 2),
        (cid[src] == cid[dst]).to(t.dtype),
    ], dim=1)


class EdgeGate(nn.Module):
    """kernel='learned' : sigmoid(MLP(edge_attr)) -- a smooth, resolution-aware
                          version of the hand-built box.
       kernel='box'     : the hand-built kernel itself, (|dt|<60) & (|dz|<2.0),
                          NOT learned. This is the control that decides whether the
                          graph bought anything beyond what round 2 already had. If
                          learned ~ box, stop.

    The gate is UNNORMALISED (sigmoid, not softmax over neighbours) for the same
    reason DeepSets pools with sigmoid gates rather than attention: normalising
    would make a 20-neighbour and a 2-neighbour node identical, and magnitude --
    how much well-measured pT agrees with you -- is one of the strongest signals
    here. A degree column is handed to the node update so the model can normalise
    if it decides it wants to; it is not forced to."""

    def __init__(self, d_e, d_h, kernel):
        super().__init__()
        self.kernel = kernel
        if kernel == "learned":
            self.net = nn.Sequential(nn.Linear(d_e, d_h), nn.ReLU(), nn.Linear(d_h, 1))

    def forward(self, ea):
        if self.kernel == "box":
            return ((ea[:, 1].abs() * 50.0 < 60.0) &
                    (ea[:, 3].abs() < 2.0)).to(ea.dtype).unsqueeze(1)
        return torch.sigmoid(self.net(ea))


class RelLayer(nn.Module):
    """One message-passing round with SEPARATE weights per edge type.

        h_i <- h_i + U([h_i, m_intra_i, m_inter_i, deg_i])
        m_*_i = SUM_{j in *} gate(e_ij) * M_*([h_j, e_ij])

    Two edge types is what carries "the existing partition is a prior": intra-cluster
    messages get their own transform from cross-cluster ones, so the model can trust
    the partition, discount it, or invert it -- and ablating an edge type gives the
    failure-mode decorrelation test with ONE model instead of two."""

    def __init__(self, d_h, d_e, kernel):
        super().__init__()
        self.m_intra = nn.Sequential(nn.Linear(d_h + d_e, d_h), nn.ReLU(),
                                     nn.Linear(d_h, d_h))
        self.m_inter = nn.Sequential(nn.Linear(d_h + d_e, d_h), nn.ReLU(),
                                     nn.Linear(d_h, d_h))
        self.g_intra = EdgeGate(d_e, d_h, kernel)
        self.g_inter = EdgeGate(d_e, d_h, kernel)
        self.upd = nn.Sequential(nn.Linear(3 * d_h + 2, d_h), nn.ReLU(),
                                 nn.Linear(d_h, d_h))
        self.norm = nn.LayerNorm(d_h)

    def forward(self, h, src, dst, ea, same):
        n = h.shape[0]
        msg = torch.cat([h[dst], ea], dim=1)
        acc, deg = [], []
        for sel, mlp, gate in ((same, self.m_intra, self.g_intra),
                               (~same, self.m_inter, self.g_inter)):
            s, m = src[sel], msg[sel]
            w = gate(ea[sel])
            a = torch.zeros(n, h.shape[1], dtype=h.dtype, device=h.device)
            a.index_add_(0, s, mlp(m) * w)
            d = torch.zeros(n, 1, dtype=h.dtype, device=h.device)
            d.index_add_(0, s, w)
            acc.append(a); deg.append(torch.log1p(d))
        return self.norm(h + self.upd(torch.cat([h] + acc + deg, dim=1)))


class EventGNN(nn.Module):
    def __init__(self, d_in, d_e, d_h=96, layers=2, kernel="learned", residual=True):
        super().__init__()
        self.residual = residual
        self.enc = nn.Sequential(nn.Linear(d_in, d_h), nn.ReLU(), nn.Linear(d_h, d_h))
        self.layers = nn.ModuleList([RelLayer(d_h, d_e, kernel) for _ in range(layers)])
        # Cluster score: pooled with a sigmoid gate (unnormalised, as in DeepSets --
        # destroying magnitude is the wrong bias when Sum pT is this strong).
        self.att = nn.Sequential(nn.Linear(d_h, d_h), nn.ReLU(), nn.Linear(d_h, 1))
        self.rho = nn.Sequential(nn.Linear(d_h, d_h), nn.ReLU(), nn.Linear(d_h, 1))
        # Per-track log-weight. Zero-init last layer => w = 1/sigma_t^2 exactly.
        self.wt = nn.Sequential(nn.Linear(d_h, d_h), nn.ReLU(), nn.Linear(d_h, 1))
        self.tau = nn.Parameter(torch.ones(1))
        if residual:
            nn.init.zeros_(self.rho[-1].weight); nn.init.zeros_(self.rho[-1].bias)
            nn.init.zeros_(self.wt[-1].weight);  nn.init.zeros_(self.wt[-1].bias)

    def forward(self, X, G, cidx, ncl, src, dst, base_score):
        ea = edge_attr(G, src, dst)
        same = ea[:, 5] > 0.5
        h = self.enc(X)
        for L in self.layers:
            h = L(h, src, dst, ea, same)

        pooled = torch.zeros(ncl, h.shape[1], dtype=h.dtype, device=h.device)
        pooled.index_add_(0, cidx, h * torch.sigmoid(self.att(h)))
        sc = self.rho(pooled).squeeze(-1)
        if self.residual:
            sc = self.tau * base_score + sc

        # t_hat: inverse-variance mean MODULATED by a learned per-track factor.
        # exp() keeps w > 0, so this can never flip a track's sign or produce a
        # negative denominator; at init the exponent is 0 and w is exactly 1/s^2.
        t, tres = G[:, 0], G[:, 1]
        lw = self.wt(h).squeeze(-1).clamp(-6.0, 6.0)
        w = torch.exp(lw) / (tres ** 2 + 1e-6)
        num = torch.zeros(ncl, dtype=h.dtype, device=h.device).index_add_(0, cidx, w * t)
        den = torch.zeros(ncl, dtype=h.dtype, device=h.device).index_add_(0, cidx, w)
        return sc, num / (den + 1e-12), lw


# ------------------------------------------------------------------ data ----

def make_tensors(frame, mu, sd, na_cols, labels, dev):
    """As train_deepsets.make_tensors, plus the raw geometry block and the event
    span table the graph builder needs. Sorted by KEY, so tracks -> clusters ->
    events are contiguous runs and an event's tracks are one slice."""
    f = frame.sort_values(KEY, kind="mergesort").reset_index(drop=True)
    X = (f[TRACK_FEATURES] - mu) / sd
    M = f[na_cols].isna().astype(np.float32).to_numpy() if na_cols else None
    X = np.nan_to_num(X.to_numpy(dtype=np.float32), nan=0.0, posinf=0.0, neginf=0.0)
    if M is not None:
        X = np.hstack([X, M])

    cl_id, cl_key = pd.factorize(pd.MultiIndex.from_frame(f[KEY]), sort=False)
    cl = pd.DataFrame(list(cl_key), columns=KEY)
    ev_id, _ = pd.factorize(pd.MultiIndex.from_frame(cl[EVT]), sort=False)
    m = cl.merge(labels, on=KEY, how="left")

    G = np.stack([f["time"].to_numpy(np.float32), f["timeRes"].to_numpy(np.float32),
                  f["z0"].to_numpy(np.float32), f["sigma_z0"].to_numpy(np.float32),
                  f["eta"].to_numpy(np.float32), f["phi"].to_numpy(np.float32),
                  cl_id.astype(np.float32)], axis=1)
    G = np.nan_to_num(G, nan=0.0, posinf=0.0, neginf=0.0)
    # timeRes = 0 would make the inverse-variance weight infinite. It is 17.5-35.0 ps
    # in every export seen, but a zero here would be silent and catastrophic.
    G[:, 1] = np.where(G[:, 1] > 1e-3, G[:, 1], 25.0)

    base = np.log(np.maximum(m["trkptz_score"].fillna(0).to_numpy(np.float32), 1e-6))
    return dict(
        X=torch.from_numpy(X).to(dev), G=torch.from_numpy(G).to(dev),
        c=torch.from_numpy(cl_id.astype(np.int64).copy()).to(dev),
        e=torch.from_numpy(ev_id.astype(np.int64).copy()).to(dev),
        y=torch.from_numpy(m["within60"].fillna(0).to_numpy(np.float32).copy()).to(dev),
        tt=torch.from_numpy(m["truth_t"].fillna(0).to_numpy(np.float32).copy()).to(dev),
        tfit=torch.from_numpy(m["time_fit"].fillna(0).to_numpy(np.float32).copy()).to(dev),
        base=torch.from_numpy(base).to(dev),
        a=torch.from_numpy(f[TRACK_LABEL].to_numpy(np.float32).copy()).to(dev),
        m=m, ncl=len(m), nev=int(ev_id.max()) + 1)


def spans(d):
    """Track spans per cluster, cluster spans per event, AND track spans per event.
    The last is new: the graph is built over an EVENT's tracks, so the batcher needs
    to slice tracks by event directly, not by walking its clusters."""
    c, e = d["c"].cpu().numpy(), d["e"].cpu().numpy()
    ev_of_tk = e[c]
    return dict(
        CS=np.searchsorted(c, np.arange(d["ncl"])),
        CE=np.searchsorted(c, np.arange(d["ncl"]), side="right"),
        ES=np.searchsorted(e, np.arange(d["nev"])),
        EE=np.searchsorted(e, np.arange(d["nev"]), side="right"),
        TS=np.searchsorted(ev_of_tk, np.arange(d["nev"])),
        TE=np.searchsorted(ev_of_tk, np.arange(d["nev"]), side="right"))


def batch_of(d, sp, evs, dev):
    """Gather one batch of whole events. Tracks are taken per EVENT (contiguous),
    and cluster ids are relabelled local-to-batch so index_add_ targets are dense."""
    cl_idx = np.concatenate([np.arange(sp["ES"][e], sp["EE"][e]) for e in evs])
    tk_idx = np.concatenate([np.arange(sp["TS"][e], sp["TE"][e]) for e in evs])
    counts = sp["TE"][evs] - sp["TS"][evs]
    src, dst = event_edges(counts)

    tk = torch.from_numpy(tk_idx).to(dev)
    cl = torch.from_numpy(cl_idx).to(dev)
    # global cluster id -> position within this batch
    remap = torch.full((d["ncl"],), -1, dtype=torch.long, device=dev)
    remap[cl] = torch.arange(len(cl_idx), device=dev)
    loc_c = remap[d["c"][tk]]
    loc_e = torch.from_numpy(
        np.repeat(np.arange(len(evs)), sp["EE"][evs] - sp["ES"][evs])).to(dev)

    G = d["G"][tk].clone()
    G[:, 6] = loc_c.to(G.dtype)          # same_cluster must compare LOCAL ids
    return dict(X=d["X"][tk], G=G, c=loc_c, ncl=len(cl_idx), e=loc_e, nev=len(evs),
                src=torch.from_numpy(src).to(dev), dst=torch.from_numpy(dst).to(dev),
                y=d["y"][cl], tt=d["tt"][cl], base=d["base"][cl])


# ------------------------------------------------------------- train/eval ---

def run_full(net, d, sp, dev, chunk=2048):
    """Forward the whole tensor in event chunks; returns scores and learned times."""
    net.eval()
    SC, TH = [], []
    with torch.no_grad():
        for i in range(0, d["nev"], chunk):
            b = batch_of(d, sp, np.arange(i, min(i + chunk, d["nev"])), dev)
            sc, th, _ = net(b["X"], b["G"], b["c"], b["ncl"], b["src"], b["dst"],
                            b["base"])
            SC.append(sc.cpu()); TH.append(th.cpu())
    net.train()
    return torch.cat(SC).numpy(), torch.cat(TH).numpy()


def tensor_ref(d):
    """TRKPTZ and the oracle through THIS tensor's own events and argmax machinery.

    The reference table printed off `df` covers the WHOLE test fold; a tensor holds a
    --test-events subsample of it, and at 1,000 events the binomial spread at ~90% is
    ~1 point. Comparing epoch 0 against the full-fold number therefore fires the
    incumbent check on sampling noise -- it did, on the first run of this script. The
    control has to be computed on the events the model is actually scored on, which is
    the same reason train_deepsets.py prints a per-tensor `trk` column."""
    out = d["m"]
    per = {}
    for sid, sub in out.groupby("sample_id"):
        g = sub.groupby(EVT, sort=False)
        pick = sub.loc[g["trkptz_score"].idxmax()]
        per[SAMPLE_NAME.get(sid, str(sid))] = {
            "nev": int(g.ngroups),
            "trkptz": float(100 * pick["within60"].mean()),
            "oracle": float(100 * g["within60"].max().mean())}
    return per


def evaluate(net, d, sp, dev):
    """Two core fractions per sample, and they answer different questions.

      cf_sel   learned SELECTION, incumbent cluster time -- comparable to
               train_deepsets.py's number, which is selection-only.
      cf_full  learned selection AND learned t0 -- the deliverable.

    Reporting only cf_full would confound the two heads; reporting only cf_sel
    would hide the head this script exists to test."""
    s, th = run_full(net, d, sp, dev)
    out = d["m"].copy()
    out["sc"], out["t_hat"] = s, th
    out["adt_fit"] = (out["time_fit"] - out["truth_t"]).abs()
    out["adt_gnn"] = (out["t_hat"] - out["truth_t"]).abs()
    per = {}
    for sid, sub in out.groupby("sample_id"):
        pick = sub.loc[sub.groupby(EVT, sort=False)["sc"].idxmax()]
        per[SAMPLE_NAME.get(sid, str(sid))] = {
            "cf_sel":  float(100 * (pick["adt_fit"] < PASS_PS).mean()),
            "cf_full": float(100 * (pick["adt_gnn"] < PASS_PS).mean())}
    return per, out


def val_loss(net, d, sp, dev, lam_t, band):
    net.eval()
    tot, n = 0.0, 0
    with torch.no_grad():
        for i in range(0, d["nev"], 2048):
            b = batch_of(d, sp, np.arange(i, min(i + 2048, d["nev"])), dev)
            sc, th, _ = net(b["X"], b["G"], b["c"], b["ncl"], b["src"], b["dst"],
                            b["base"])
            l = listwise_ce(sc, b["e"], b["y"], b["nev"])
            l = l + lam_t * time_loss(th, b["tt"], band)
            tot += l.item(); n += 1
    net.train()
    return tot / max(n, 1)


def time_loss(th, tt, band):
    """Huber on (t_hat - truth), restricted to clusters already within `band` of the
    truth HS time.

    The restriction matters. Asking a pure-pileup cluster to output the hard-scatter
    time is asking for something unlearnable -- it has no HS tracks to weight up. The
    band keeps the target to clusters that plausibly ARE the HS one, and it is wider
    than the 60 ps window on purpose: the events re-weighting can actually flip are
    the near-misses just outside it, not the ones already passing."""
    d = (th - tt) / PASS_PS          # residual in units of the pass window
    m = d.abs() * PASS_PS < band
    if not bool(m.any()):
        return th.sum() * 0.0
    # delta=0.5 window-units == 30 ps: the same physical quadratic/linear transition
    # the raw-ps version used, but the loss is now O(1) instead of O(10^3). In raw
    # picoseconds the Huber term outweighed the listwise CE ~1000:1 and the selection
    # head was effectively not being trained at lam_time=1.
    return F.huber_loss(d[m], torch.zeros_like(d[m]), delta=0.5)


def dump(args, ctx, hist, best_ep, final=None):
    """Rewrite results.json. Called after EVERY epoch, not only at the end.

    A run that dies at epoch 9 of 12 -- OOM, a killed shell, a laptop lid -- used to
    lose its entire history along with the epochs that did finish, since the single
    write happened after the final test evaluation. `in_progress` says whether the
    test block is there yet, so a partial file is never mistaken for a finished run."""
    d = {"args": vars(args), **ctx, "history": hist, "best_epoch": best_ep,
         "in_progress": final is None}
    if final is not None:
        d.update(final)
    tmp = os.path.join(args.out, "results.json.tmp")
    with open(tmp, "w") as fh:
        json.dump(d, fh, indent=2)
    os.replace(tmp, os.path.join(args.out, "results.json"))   # atomic; never truncated


def train(net, FIT, VAL, spF, spV, args, dev, ctx):
    opt = torch.optim.Adam(net.parameters(), lr=args.lr)
    sched = torch.optim.lr_scheduler.CosineAnnealingLR(opt, T_max=max(args.epochs, 1),
                                                       eta_min=args.lr / 20)
    rng = np.random.default_rng(args.init_seed)
    order = np.arange(FIT["nev"])
    best_vl, best_state, best_ep, hist = float("inf"), copy.deepcopy(net.state_dict()), 0, []
    for ep in range(args.epochs):
        rng.shuffle(order)
        tot = nb = 0
        t0 = time.time()
        for i in range(0, len(order), args.batch_events):
            b = batch_of(FIT, spF, order[i:i + args.batch_events], dev)
            sc, th, lw = net(b["X"], b["G"], b["c"], b["ncl"], b["src"], b["dst"],
                             b["base"])
            loss = (listwise_ce(sc, b["e"], b["y"], b["nev"])
                    + args.lam_time * time_loss(th, b["tt"], args.time_band))
            opt.zero_grad(); loss.backward()
            nn.utils.clip_grad_norm_(net.parameters(), 5.0)
            opt.step()
            tot += loss.detach().item(); nb += 1
        sched.step()
        vl = val_loss(net, VAL, spV, dev, args.lam_time, args.time_band)
        vper, _ = evaluate(net, VAL, spV, dev)
        star = ""
        # SELECT ON VAL LOSS, never on val core fraction -- train_deepsets.val_loss()
        # documents the collapse (17.6% on weights that scored 91.7% on test).
        if vl < best_vl:
            best_vl, best_state, best_ep, star = vl, copy.deepcopy(net.state_dict()), ep + 1, "  *"
        hist.append({"epoch": ep + 1, "loss": tot / max(nb, 1), "val_loss": vl,
                     "val": vper})
        log(f"  epoch {ep+1:>2}/{args.epochs}  loss {tot/max(nb,1):.4f}  "
            f"vloss {vl:.4f}  [{time.time()-t0:.0f}s]  VAL "
            + "  ".join(f"{k} sel {v['cf_sel']:.1f}% full {v['cf_full']:.1f}%"
                        for k, v in vper.items()) + star)
        dump(args, ctx, hist, best_ep)
    net.load_state_dict(best_state)
    log(f"  restored epoch {best_ep} (val loss {best_vl:.4f})")
    return best_ep, hist


def main():
    global TRACK_FEATURES
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--input-dir", required=True)
    p.add_argument("--out", required=True)
    p.add_argument("--selection", choices=["canonical", "loose"], default="canonical")
    p.add_argument("--samples", default="",
                   help="comma-separated samples to LOAD (default: all). Loading one "
                        "sample is how the local box runs this at all -- each export "
                        "is 1.4-2.2 GB.")
    p.add_argument("--layers", type=int, default=2)
    p.add_argument("--hidden", type=int, default=96)
    p.add_argument("--edge-kernel", choices=["learned", "box"], default="learned",
                   help="'box' freezes the gate to the hand-built (|dt|<60)&(|dz|<2) "
                        "kernel. THE control: if learned ~ box, the graph bought "
                        "nothing that nb_phzt did not already have.")
    p.add_argument("--no-residual", action="store_true",
                   help="drop the incumbent initialisation (score from scratch, "
                        "w from scratch). The ablation for 'start at TRKPTZ'.")
    p.add_argument("--lam-time", type=float, default=1.0,
                   help="weight on the t0 Huber term; 0 = selection only")
    p.add_argument("--time-band", type=float, default=150.0,
                   help="ps; clusters further than this from truth are not time targets")
    p.add_argument("--train-events", type=int, default=40_000)
    p.add_argument("--test-events", type=int, default=20_000)
    p.add_argument("--epochs", type=int, default=10)
    p.add_argument("--batch-events", type=int, default=256)
    p.add_argument("--lr", type=float, default=2e-3)
    p.add_argument("--val-lo", type=float, default=0.60)
    p.add_argument("--step", default="300 MB")
    p.add_argument("--seed", type=int, default=0, help="fixes the EVENT SPLIT")
    p.add_argument("--init-seed", type=int, default=None, help="fixes WEIGHTS/batch order")
    p.add_argument("--torch-threads", type=int, default=0)
    p.add_argument("--drop-features", default="",
                   help="comma-separated TRACK_FEATURES to remove. The A/B knob for a "
                        "feature ablation: run once with the list and once without at "
                        "the same --seed, so the split and the events are IDENTICAL and "
                        "only the inputs differ. Use this rather than comparing two "
                        "exports -- the local and grid VBF files are different "
                        "productions, so a cross-file comparison confounds the feature "
                        "set with the sample.")
    p.add_argument("--device", default="auto", choices=["auto", "cpu", "cuda", "mps"])
    args = p.parse_args()
    if args.init_seed is None:
        args.init_seed = args.seed
    os.makedirs(args.out, exist_ok=True)
    if args.torch_threads > 0:
        torch.set_num_threads(args.torch_threads)
    if args.device == "auto":
        dev = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    else:
        dev = torch.device(args.device)
    log(f"torch {torch.__version__}  device {dev}  threads {torch.get_num_threads()}")

    # find_files reads train_deepsets.SAMPLES; narrowing it is how --samples works.
    if args.samples:
        want = [s.strip() for s in args.samples.split(",") if s.strip()]
        bad = [s for s in want if s not in ds.SAMPLES]
        if bad:
            sys.exit(f"FATAL: --samples names {bad}; known: {ds.SAMPLES}")
        ds.SAMPLES = want
        log(f"loading only: {', '.join(want)}")

    files = find_files(args.input_dir, args.selection)
    if not files:
        sys.exit(f"FATAL: no *_training.root under {args.input_dir}")

    parts = []
    for name, path in files.items():
        d = ds.read_tree(path, "clusters", CLUSTER_COLS)
        log(f"  {name:6s} {len(d):>9,} clusters  {d.groupby(EVT).ngroups:>8,} events")
        parts.append(d)
    df = pd.concat(parts, ignore_index=True)
    dup = int(df.duplicated(KEY).sum())
    if dup:
        sys.exit(f"FATAL: {dup:,} duplicate {KEY} rows")
    df["abs_dt"] = df["delta_t"].abs()
    df["within60"] = (df["abs_dt"] < PASS_PS).astype(np.float32)
    df["time_fit"] = df["cluster_time"]
    df["truth_t"] = df["cluster_time"] - df["delta_t"]

    # truth_t must be one number per event: it is the truth HS vertex time, and both
    # inputs are per-cluster measurements OF THAT SAME EVENT. A spread here means the
    # key is wrong (the shard-invariance bug this key exists to prevent) and every
    # listwise group downstream is a mixture of two events.
    spread = df.groupby(EVT)["truth_t"].agg(lambda s: s.max() - s.min())
    if float(spread.max()) > 1e-3:
        sys.exit(f"FATAL: truth_t varies within {int((spread > 1e-3).sum()):,} events "
                 f"(max spread {float(spread.max()):.4g} ps) -- event key is broken")
    log(f"total {len(df):,} clusters / {df.groupby(EVT).ngroups:,} events "
        f"(truth_t event-constant: OK)")

    rng = np.random.default_rng(args.seed)
    ev = df[EVT].drop_duplicates().reset_index(drop=True)
    ev["fold"] = rng.random(len(ev))
    fold = ev.set_index(EVT)["fold"]
    df = df.merge(ev, on=EVT, how="left")

    log("\nreference selectors (test events):")
    ref = {}
    for sid, sub in df[df.fold >= 0.7].groupby("sample_id"):
        nm = SAMPLE_NAME.get(sid, str(sid))
        ref[nm] = {"TRKPTZ": core_fraction(sub, "trkptz_score"),
                   "WAVeS": core_fraction(sub, "waves_score"),
                   "oracle": core_fraction(sub, "abs_dt", False)}
        log(f"  {nm:6s} TRKPTZ {ref[nm]['TRKPTZ']:5.1f}%  WAVeS {ref[nm]['WAVeS']:5.1f}%"
            f"  oracle {ref[nm]['oracle']:5.1f}%")

    need = sorted(set(TRACK_FEATURES + GEOM + KEY + [TRACK_LABEL]))
    per_train = {s: args.train_events // max(1, len(files)) for s in files}
    per_test = args.test_events // max(1, len(files))
    log("\nstreaming tracks:")
    TR_T, TE_T = _stream(files, fold, per_train, per_test, args.step, need)
    if TR_T is None:
        sys.exit("FATAL: nothing to fit on")

    # An export older than the extended AF re-export simply does not have the
    # Group-1/Group-3 columns (on_pv, chi2_ndf, btag_*, the vertex-relational block).
    # iter_tracks drops absent columns SILENTLY, so without this the model would
    # train on a quietly different feature set than the one named in the source --
    # and two runs against two exports would not be comparable while looking it.
    # Narrow explicitly, name every casualty, and record the survivors in the
    # checkpoint so what a model actually saw is recoverable from the file.
    if args.drop_features:
        drop = [f.strip() for f in args.drop_features.split(",") if f.strip()]
        unknown = [f for f in drop if f not in TRACK_FEATURES]
        if unknown:
            sys.exit(f"FATAL: --drop-features names features not in TRACK_FEATURES: "
                     f"{unknown}")
        TRACK_FEATURES = [f for f in TRACK_FEATURES if f not in drop]
        log(f"\ndropped {len(drop)} feature(s) by request; {len(TRACK_FEATURES)} remain")

    absent = [f for f in TRACK_FEATURES if f not in TR_T.columns]
    if absent:
        TRACK_FEATURES = [f for f in TRACK_FEATURES if f in TR_T.columns]
        log(f"\n  !! this export lacks {len(absent)} of the declared track features; "
            f"{len(TRACK_FEATURES)} remain. NOT comparable with a run on an export "
            f"that has them.")
        log("     missing: " + ", ".join(absent))

    FIT_T = TR_T[TR_T.fold < args.val_lo]
    VAL_T = TR_T[TR_T.fold >= args.val_lo]
    log(f"tracks: fit {len(FIT_T):,} / {FIT_T.groupby(EVT).ngroups:,} ev   "
        f"val {len(VAL_T):,} / {VAL_T.groupby(EVT).ngroups:,} ev   "
        f"test {len(TE_T):,} / {TE_T.groupby(EVT).ngroups:,} ev")

    mu = FIT_T[TRACK_FEATURES].mean()
    sd = FIT_T[TRACK_FEATURES].std().replace(0, 1.0)
    na_cols = [c for c in TRACK_FEATURES if FIT_T[c].isna().any()]
    labels = df[KEY + ["within60", "truth_t", "time_fit", "trkptz_score"]].drop_duplicates(KEY)

    FIT = make_tensors(FIT_T, mu, sd, na_cols, labels, dev)
    VAL = make_tensors(VAL_T, mu, sd, na_cols, labels, dev)
    TST = make_tensors(TE_T, mu, sd, na_cols, labels, dev)
    spF, spV, spT = spans(FIT), spans(VAL), spans(TST)

    torch.manual_seed(args.init_seed)
    net = EventGNN(FIT["X"].shape[1], len(EDGE_FEATURES), args.hidden, args.layers,
                   args.edge_kernel, residual=not args.no_residual).to(dev)
    log(f"\nparameters: {sum(p.numel() for p in net.parameters()):,}   "
        f"layers {args.layers}  kernel {args.edge_kernel}  "
        f"residual {not args.no_residual}")

    # THE STRUCTURAL CONTROL. With the residual parameterisation, the untrained
    # network IS the incumbent: score = tau*log(trkptz_score), w = 1/sigma_t^2. So
    # cf_sel here must equal the TRKPTZ reference above, and cf_full must equal
    # TRKPTZ-selection-with-inverse-variance-timing. A mismatch is a tensor bug, not
    # a model result -- exactly the ambiguity train_deepsets' trk control column
    # exists to remove, made structural instead of bolted on.
    # Composition of each tensor. orc is a property of the EVENTS alone, so a VAL/TEST
    # disagreement there means the slices differ in what is even achievable and no
    # model number between them is comparable. trk is TRKPTZ through that tensor's own
    # argmax -- the number epoch 0 must reproduce.
    log("\ntensor composition (trk = TRKPTZ on THESE events; orc = ceiling):")
    trefs = {}
    for nm, d in (("fit", FIT), ("val", VAL), ("test", TST)):
        trefs[nm] = tensor_ref(d)
        log(f"  {nm:5s} " + "   ".join(
            f"{k} {v['nev']:>6,} ev  trk {v['trkptz']:5.1f}%  orc {v['oracle']:5.1f}%"
            for k, v in trefs[nm].items()))
    tref = trefs["test"]

    per0, _ = evaluate(net, TST, spT, dev)
    log("\nepoch 0 (untrained == incumbent by construction):")
    for k, v in per0.items():
        log(f"  {k:6s} cf_sel {v['cf_sel']:5.1f}%  (TRKPTZ on these events "
            f"{tref[k]['trkptz']:5.1f}%)   cf_full {v['cf_full']:5.1f}%")
        if not args.no_residual and abs(v["cf_sel"] - tref[k]["trkptz"]) > 0.05:
            log(f"  !! {k}: epoch-0 selection does not reproduce TRKPTZ on its own "
                f"tensor -- that is a tensor/wiring bug, not a model result")

    ctx = {"reference": ref, "tensor_ref": trefs, "features": TRACK_FEATURES,
           "epoch0": per0}
    log("")
    best_ep, hist = train(net, FIT, VAL, spF, spV, args, dev, ctx)
    per, _ = evaluate(net, TST, spT, dev)

    log("\n" + "=" * 72)
    names = sorted(per)
    log(f"{'selector':30s}" + "".join(f"{n:>11s}" for n in names))
    log(f"{'TRKPTZ (this tensor)':30s}" + "".join(f"{tref[n]['trkptz']:10.1f}%" for n in names))
    log(f"{'GNN [sel only]':30s}" + "".join(f"{per[n]['cf_sel']:10.1f}%" for n in names))
    log(f"{'GNN [sel + learned t0]':30s}" + "".join(f"{per[n]['cf_full']:10.1f}%" for n in names))
    log(f"{'oracle (this tensor)':30s}" + "".join(f"{tref[n]['oracle']:10.1f}%" for n in names))
    log(f"\n(full test fold, for scale: "
        + "; ".join(f"{n} TRKPTZ {ref[n]['TRKPTZ']:.1f}% WAVeS {ref[n]['WAVeS']:.1f}%"
                    for n in names) + ")")

    dump(args, ctx, hist, best_ep, final={"test": per})
    torch.save({"state_dict": net.state_dict(), "mu": mu.to_dict(), "sd": sd.to_dict(),
                "na_cols": na_cols, "features": TRACK_FEATURES, "hidden": args.hidden,
                "layers": args.layers, "kernel": args.edge_kernel,
                "residual": not args.no_residual},
               os.path.join(args.out, "best_model.pt"))
    log(f"\nwrote {args.out}/results.json and {args.out}/best_model.pt")


def _stream(files, fold, per_train, per_test, step, need):
    """train_deepsets.stream_tracks with the column list widened to include GEOM.

    Its own version closes over TRACK_FEATURES at call time, which would drop the raw
    time/z0/eta/phi the graph builds edges from. Same per-sample-quota logic, which is
    the part that matters -- a global cap filled sequentially reads the first file only
    and silently trains on one topology, a bug that appeared three times in this study."""
    keep_tr, keep_te = [], []
    for name, path in files.items():
        cap_tr, cap_te = per_train[name], per_test
        ntr = nte = 0
        for b in ds.iter_tracks(path, need, step):
            b = b.assign(fold=fold.reindex(pd.MultiIndex.from_arrays(
                [b[c] for c in EVT])).to_numpy())
            for is_tr in (True, False):
                cap, got = (cap_tr, ntr) if is_tr else (cap_te, nte)
                if got >= cap:
                    continue
                sel = b[(b.fold < 0.7) if is_tr else (b.fold >= 0.7)]
                if not len(sel):
                    continue
                evs = sel[EVT].drop_duplicates().iloc[:max(0, cap - got)]
                sel = sel.merge(evs, on=EVT, how="inner")
                (keep_tr if is_tr else keep_te).append(sel)
                if is_tr:
                    ntr += len(evs)
                else:
                    nte += len(evs)
            if ntr >= cap_tr and nte >= cap_te:
                break
        log(f"  {name:12s} train {ntr:>7,} ev   test {nte:>7,} ev"
            + ("   (all available)" if ntr < cap_tr else ""))
    return (pd.concat(keep_tr, ignore_index=True) if keep_tr else None,
            pd.concat(keep_te, ignore_index=True))


if __name__ == "__main__":
    main()
