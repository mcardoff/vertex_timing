#!/usr/bin/env python3
"""
model_inputs_report.py -- dump the exact input list and per-feature importance
for each of the five models in the zjets_seventy / final_method_eval pipeline.

This is a DOCUMENTATION run, not a performance measurement: each model is fit
ONCE on the full sample with the production hyperparameters rather than
4-fold x n-seed cross-fitted, because importances do not need out-of-fold
predictions and the full protocol costs 4-20x more. Every number here describes
which inputs a model leans on, never how well it scores -- the scores in
results/final_method_eval.json are the cross-fitted ones and stand unchanged.

Fits, in pipeline order:
  1  tagger round 1   xgb(7,10,n=600)  on track rows            -> gain by feature
  2  tagger round 2   same + 6 neighbor-context columns
  3  stage-1 score    xgb(6,20,n=500)  on cluster rows, label ok_tagw
  4  e2e WeightNet    torch, permutation importance on core fraction
  5  meta chooser     xgb(4,10,n=500)  on the disagreement subset

Writes JSON to --out. Run against training_new/ (the live export).
"""
import argparse, gc, json, os, subprocess, sys

import numpy as np
import pandas as pd
import uproot
import xgboost as xgb

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from final_method_eval import EVT, KEY, shrink, keyfold, log

TOPN = 15


def xgbc(depth, mcw, n, lr, jobs):
    return xgb.XGBClassifier(n_estimators=n, learning_rate=lr, max_depth=depth,
                             subsample=0.8, colsample_bytree=0.8,
                             min_child_weight=mcw, eval_metric="logloss",
                             n_jobs=jobs, tree_method="hist", random_state=0)


def imps(model, feats):
    """Normalised gain, descending. Gain is the metric that answers 'how much
    did splitting on this reduce the loss', which is the question here; the
    default 'weight' just counts splits and over-rewards high-cardinality
    continuous columns."""
    g = model.get_booster().get_score(importance_type="gain")
    v = np.array([g.get(f"f{i}", 0.0) for i in range(len(feats))])
    v = v / v.sum() if v.sum() > 0 else v
    o = np.argsort(-v)
    return [{"f": feats[i], "imp": round(100 * float(v[i]), 2)} for i in o]


def auc(x, y):
    """Rank-based (Mann-Whitney) AUC of a single variable against a binary
    label. NaNs are dropped rather than imputed -- a sentinel value would land
    at one end of the ranking and manufacture separation."""
    m = np.isfinite(x)
    x, y = x[m], y[m].astype(bool)
    ns, nb = int(y.sum()), int((~y).sum())
    if ns == 0 or nb == 0: return float("nan")
    r = pd.Series(x).rank().to_numpy()
    return float((r[y].sum() - ns * (ns + 1) / 2.0) / (ns * nb))


def panel(x, y, name, step, split, sig_lab, bkg_lab, nbins=44):
    """One signal-vs-background comparison. Few distinct values -> categorical
    bars; otherwise a density histogram over the combined [0.5, 99.5] range so
    one long tail cannot flatten the whole picture. Both classes are normalised
    to unit area, since the classes differ in size by up to 30x and raw counts
    would only show which is bigger."""
    x = np.asarray(x, float); y = np.asarray(y).astype(bool)
    m = np.isfinite(x); x, y = x[m], y[m]
    a = auc(x, y)
    if a == a and a < 0.5:        # orient so signal sits high; note the flip
        a, flipped = 1.0 - a, True
    else:
        flipped = False
    uq = np.unique(x)
    out = {"step": step, "var": name, "split": split, "sig": sig_lab,
           "bkg": bkg_lab, "auc": None if a != a else round(a, 3),
           "flipped": flipped, "n_sig": int(y.sum()), "n_bkg": int((~y).sum())}
    if len(uq) <= 6:
        out["kind"] = "cat"
        out["labels"] = [f"{v:g}" for v in uq]
        out["s"] = [round(float((x[y] == v).mean()), 5) for v in uq]
        out["b"] = [round(float((x[~y] == v).mean()), 5) for v in uq]
        return out
    lo, hi = np.percentile(x, [0.5, 99.5])
    if hi <= lo: hi = lo + 1e-6
    e = np.linspace(lo, hi, nbins + 1)
    xc = np.clip(x, lo, hi)
    hs, _ = np.histogram(xc[y], bins=e, density=True)
    hb, _ = np.histogram(xc[~y], bins=e, density=True)
    out.update(kind="hist", edges=[round(float(v), 5) for v in e],
               s=[round(float(v), 8) for v in hs],
               b=[round(float(v), 8) for v in hb],
               clip_hi=round(100 * float((x > hi).mean()), 2),
               median_s=round(float(np.median(x[y])), 4),
               median_b=round(float(np.median(x[~y])), 4))
    return out


E2E_IMP = r'''
import json,sys,numpy as np,pandas as pd,torch,torch.nn as nn
wd=sys.argv[1]
cfg=json.load(open(f"{wd}/e2e_cfg.json")); EVT=cfg["evt"]; TF2=cfg["feats"]
tz2=pd.read_pickle(f"{wd}/e2e_in.pkl"); tz2=tz2.loc[:,~tz2.columns.duplicated()]
torch.set_num_threads(int(sys.argv[2]))
class WNet(nn.Module):
    def __init__(self,n,h=96,it=5,floor=10.0):
        super().__init__(); self.it=it; self.floor=floor
        self.f=nn.Sequential(nn.Linear(n,h),nn.ReLU(),nn.Linear(h,h),nn.ReLU(),nn.Linear(h,1))
    def forward(self,x,t,e,n_ev):
        w=torch.nn.functional.softplus(self.f(x)).squeeze(-1)+1e-6
        num=torch.zeros(n_ev).index_add_(0,e,w*t); den=torch.zeros(n_ev).index_add_(0,e,w)
        t0=num/den.clamp_min(1e-9)
        for _ in range(self.it):
            r=(t-t0[e]).abs().detach().clamp_min(self.floor); u=w/r
            num=torch.zeros(n_ev).index_add_(0,e,u*t); den=torch.zeros(n_ev).index_add_(0,e,u)
            t0=num/den.clamp_min(1e-9)
        return t0
mu=tz2[TF2].mean(); sd=tz2[TF2].std().replace(0,1.0)
def pack(df):
    cod,uq=pd.factorize(pd.MultiIndex.from_frame(df[EVT]),sort=False)
    X=torch.from_numpy(((df[TF2]-mu)/sd).fillna(0).to_numpy(np.float32))
    return (X,torch.from_numpy(df["time"].to_numpy(np.float32)),
            torch.from_numpy(cod.astype(np.int64)),
            torch.from_numpy(df.groupby(EVT,sort=False)["yt"].first().to_numpy(np.float32)),
            len(uq))
tr_=tz2[tz2.k!=0].sort_values(EVT); te_=tz2[tz2.k==0].sort_values(EVT)
Xa,ta,ea,ya,na=pack(tr_); Xe,te2,ee,ye,ne_=pack(te_)
torch.manual_seed(0)
net=WNet(len(TF2)); opt=torch.optim.Adam(net.parameters(),lr=2e-3)
best=None;bacc=-1
for ep in range(1,201):
    tau=float(np.interp(ep,[1,200],[80.0,25.0])); net.train()
    t0b=net(Xa,ta,ea,na)
    loss=(1.0-torch.sigmoid((60.0-(t0b-ya).abs())/tau)).mean()
    opt.zero_grad(); loss.backward()
    torch.nn.utils.clip_grad_norm_(net.parameters(),1.0); opt.step()
    if ep%5==0 or ep==200:
        net.eval()
        with torch.no_grad(): acc=float(((net(Xe,te2,ee,ne_)-ye).abs()<60).float().mean())
        if acc>bacc: bacc=acc; best={k:v.clone() for k,v in net.state_dict().items()}
net.load_state_dict(best); net.eval()
with torch.no_grad(): base=float(((net(Xe,te2,ee,ne_)-ye).abs()<60).float().mean())
# permutation importance: shuffle one column of the held-out fold, re-run the
# whole weight+IRLS forward pass, record the core-fraction drop in points.
rng=np.random.default_rng(0); out=[]
for i,f in enumerate(TF2):
    Xp=Xe.clone(); Xp[:,i]=Xe[rng.permutation(len(Xe)),i]
    with torch.no_grad(): a=float(((net(Xp,te2,ee,ne_)-ye).abs()<60).float().mean())
    out.append({"f":f,"drop":round(100*(base-a),3)})
out.sort(key=lambda d:-d["drop"])
json.dump({"base":round(100*base,2),"perm":out},open(f"{wd}/e2e_imp.json","w"),indent=1)
print(f"e2e fold0 val {100*base:.2f}%",flush=True)
'''


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--input-dir", required=True)
    ap.add_argument("--workdir", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--python", default=sys.executable)
    ap.add_argument("--jobs", type=int, default=6)
    a = ap.parse_args()
    D, W, J = a.input_dir, a.workdir, a.jobs
    os.makedirs(W, exist_ok=True)
    R = {}

    files = [f"{D}/zjets_novbs_training.root", f"{D}/zjets_mjj500p0_training.root"]
    log("loading zjets (novbs + canonical)")
    cs = [shrink(uproot.open(f)["clusters"].arrays(library="pd")) for f in files]
    k0 = set(map(tuple, cs[0][EVT].drop_duplicates().to_numpy()))
    cs[1] = cs[1][~pd.MultiIndex.from_frame(cs[1][EVT]).isin(k0)]
    c = pd.concat(cs, ignore_index=True); del cs
    c["ok"] = c["delta_t"].abs() < 60
    c["t_truth"] = c["cluster_time"] - c["delta_t"]
    c["k"] = keyfold(c)
    ts = []
    for f in files:
        x = shrink(uproot.open(f)["tracks"].arrays(library="pd"))
        ts.append(x[(x.time_valid > 0) & (x.timeRes > 0)])
    ts[1] = ts[1][~pd.MultiIndex.from_frame(ts[1][EVT]).isin(
        set(map(tuple, ts[0][EVT].drop_duplicates().to_numpy())))]
    t = pd.concat(ts, ignore_index=True).drop_duplicates(EVT + ["track_idx"]); del ts
    t["k"] = keyfold(t)
    log(f"  {len(c):,} clusters  {len(t):,} timed tracks  "
        f"{c[EVT].drop_duplicates().shape[0]:,} events")

    # ---- 1. tagger round 1 -------------------------------------------------
    TFE = [x for x in t.columns if x not in EVT and x not in
           {"track_idx", "cluster_idx", "k", "weight"} and not x.startswith("truth_")
           and t[x].dtype != object]
    log(f"[1] tagger round 1: {len(TFE)} inputs, {len(t):,} rows")
    m = xgbc(7, 10, 600, 0.05, J)
    m.fit(t[TFE].to_numpy(np.float32), t.truth_is_hs.to_numpy(int))
    R["tagger_r1"] = {"n_features": len(TFE), "features": TFE,
                      "n_rows": int(len(t)), "top": imps(m, TFE)[:TOPN],
                      "hs_frac": round(100 * float(t.truth_is_hs.mean()), 2)}
    t["ph"] = m.predict_proba(t[TFE].to_numpy(np.float32))[:, 1]
    del m; gc.collect()

    # ---- 2. tagger round 2 (neighbor context) ------------------------------
    nb = t[EVT + ["track_idx", "time", "z0", "ph"]].copy()
    ekeys = nb[EVT].drop_duplicates().reset_index(drop=True)
    NBs = []
    for i0 in range(0, len(ekeys), 60000):
        sub = nb.merge(ekeys.iloc[i0:i0 + 60000], on=EVT, how="inner")
        j = sub.merge(sub, on=EVT, suffixes=("", "_o"))
        j = j[j.track_idx != j.track_idx_o]
        dt_ = np.abs(j["time"] - j["time_o"]); dz_ = np.abs(j["z0"] - j["z0_o"])
        gI = [j[k] for k in EVT] + [j["track_idx"]]
        NBs.append(pd.DataFrame({
            "nb_n60": (dt_ < 60).groupby(gI).sum(),
            "nb_ph60": (j["ph_o"] * (dt_ < 60)).groupby(gI).sum(),
            "nb_n30": (dt_ < 30).groupby(gI).sum(),
            "nb_ph30": (j["ph_o"] * (dt_ < 30)).groupby(gI).sum(),
            "nb_phz2": (j["ph_o"] * (dz_ < 2.0)).groupby(gI).sum(),
            "nb_phzt": (j["ph_o"] * ((dz_ < 2.0) & (dt_ < 60))).groupby(gI).sum()}))
        del sub, j; gc.collect()
    NB = pd.concat(NBs); NB.index.names = EVT + ["track_idx"]; del NBs
    t = t.merge(NB.reset_index(), on=EVT + ["track_idx"], how="left")
    for f in NB.columns: t[f] = t[f].fillna(0)
    TFE2 = TFE + list(NB.columns)
    log(f"[2] tagger round 2: {len(TFE2)} inputs (+{len(NB.columns)} neighbor)")
    m = xgbc(7, 10, 600, 0.05, J)
    m.fit(t[TFE2].to_numpy(np.float32), t.truth_is_hs.to_numpy(int))
    ii = imps(m, TFE2)
    R["tagger_r2"] = {"n_features": len(TFE2), "new_features": list(NB.columns),
                      "top": ii[:TOPN],
                      "neighbor_share": round(sum(d["imp"] for d in ii
                                                  if d["f"] in set(NB.columns)), 2)}
    t["ph"] = m.predict_proba(t[TFE2].to_numpy(np.float32))[:, 1]
    del m, NB; gc.collect()

    # ---- t0 aggregate + per-cluster tag-weighted time ----------------------
    d = t.sort_values(EVT + ["time"])
    codes, uniq = pd.factorize(pd.MultiIndex.from_frame(d[EVT]), sort=False)
    if not isinstance(uniq, pd.MultiIndex):
        uniq = pd.MultiIndex.from_tuples([tuple(x) for x in uniq], names=EVT)
    w = d["ph"].to_numpy(); tm = d["time"].to_numpy()
    tot = np.bincount(codes, weights=w)
    start = np.r_[0, np.flatnonzero(np.diff(codes)) + 1]
    cw = np.cumsum(w); off = np.zeros(len(d)); off[start[1:]] = cw[start[1:] - 1]
    cw = cw - np.maximum.accumulate(off)
    hit = cw >= 0.5 * tot[codes]; idx = np.arange(len(d))
    firsts = np.minimum.reduceat(np.where(hit, idx, len(d) + 1), start)
    med = pd.Series(tm[np.clip(firsts, 0, len(d) - 1)], index=uniq); med.index.names = EVT

    def trunc(center, win):
        d2 = t.merge(center.rename("c0").rename_axis(EVT).reset_index(), on=EVT, how="left")
        d2 = d2[(d2["time"] - d2["c0"]).abs() < win]
        d2["wp"] = d2["ph"] / d2["timeRes"] ** 2
        r = d2.groupby(EVT, sort=False).apply(
            lambda x: (x.wp * x.time).sum() / x.wp.sum(), include_groups=False)
        r.index.names = EVT
        return r.reindex(med.index).fillna(med)
    t0f = trunc(trunc(med, 100), 60)

    mi = pd.MultiIndex.from_frame(c[EVT])
    c["t0_agg"] = med.reindex(mi).to_numpy()
    c["sum_ph"] = pd.Series(tot, index=med.index).reindex(mi).to_numpy()
    c["dt_to_agg"] = (c["cluster_time"] - c["t0_agg"]).abs()
    t["ph_pt"] = t["ph"] * t["pt"]
    t["w_t"] = t["ph"] / t["timeRes"] ** 2
    t["wt_t"] = t["w_t"] * t["time"]
    gcl = t.groupby(KEY, sort=False)
    cand = pd.DataFrame({"ph_sum": gcl["ph"].sum(), "ph_mean": gcl["ph"].mean(),
                         "ph_max": gcl["ph"].max(), "ph_pt_sum": gcl["ph_pt"].sum(),
                         "n_ph05": gcl["ph"].apply(lambda x: (x > 0.5).sum()),
                         "t_tagw": gcl["wt_t"].sum() / gcl["w_t"].sum()})
    tw0 = cand["t_tagw"].copy()
    tw = cand["t_tagw"]
    for win in (100.0, 60.0):
        tj = t.merge(tw.rename("tw").reset_index(), on=KEY, how="left")
        tj = tj[(tj["time"] - tj["tw"]).abs() < win]
        gg2 = tj.groupby(KEY, sort=False)
        tw = (gg2["wt_t"].sum() / gg2["w_t"].sum()).reindex(cand.index).fillna(tw)
    cand["t_tagw"] = tw
    cand["ph_pt_frac"] = cand["ph_pt_sum"] / gcl["pt"].sum()
    c = c.merge(cand.reset_index(), on=KEY, how="left")
    c["dt_tagw_agg"] = (c["t_tagw"] - c["t0_agg"]).abs()
    c["t_tagw_raw"] = tw0.reindex(pd.MultiIndex.from_frame(c[KEY])).to_numpy()
    tw_eff = np.where(np.isnan(c["t_tagw"]), c["cluster_time"], c["t_tagw"])
    c["ok_tagw"] = np.abs(tw_eff - c["t_truth"]) < 60

    # ---- 3. stage-1 cluster score -----------------------------------------
    drop = {"ok", "ok_tagw", "k", "p", "p_old", "delta_t", "weight", "cluster_idx",
            "t_truth", "t_tagw_raw"}
    FEA = [x for x in c.columns if x not in EVT and x not in drop
           and not x.startswith("truth_") and c[x].dtype != object]
    NEWC = ["t0_agg", "sum_ph", "dt_to_agg", "ph_sum", "ph_mean", "ph_max",
            "ph_pt_sum", "n_ph05", "t_tagw", "ph_pt_frac", "dt_tagw_agg"]
    log(f"[3] stage-1 cluster score: {len(FEA)} inputs, {len(c):,} rows")
    m = xgbc(6, 20, 500, 0.05, J)
    m.fit(c[FEA].to_numpy(np.float32), c["ok_tagw"].to_numpy(int))
    ii = imps(m, FEA)
    R["stage1"] = {"n_features": len(FEA), "features": FEA, "n_rows": int(len(c)),
                   "top": ii[:TOPN],
                   "tagger_derived": [x for x in NEWC if x in FEA],
                   "tagger_share": round(sum(d["imp"] for d in ii
                                             if d["f"] in set(NEWC)), 2),
                   "trkptz_share": round(sum(d["imp"] for d in ii
                                             if d["f"] == "trkptz_score"), 2)}
    c["p"] = m.predict_proba(c[FEA].to_numpy(np.float32))[:, 1]
    del m; gc.collect()

    # step-3 arithmetic ablation: what each Winsorisation pass is worth, on the
    # TRKPTZ-picked cluster, so it is a pure re-timing delta with no selection change
    g = c.groupby(EVT, sort=False)
    tp = c.loc[g["trkptz_score"].idxmax()]
    ok = lambda col: round(100 * float((np.abs(
        np.where(np.isnan(tp[col]), tp["cluster_time"], tp[col]) - tp["t_truth"])
        < 60).mean()), 2)
    R["step3_ablation"] = {"cluster_time_invvar": ok("cluster_time"),
                           "tagw_no_winsor": ok("t_tagw_raw"),
                           "tagw_double_winsor": ok("t_tagw"),
                           "note": "on the TRKPTZ-chosen cluster; selection identical, "
                                   "only the reported time differs"}

    # ---- 4. e2e WeightNet permutation importance ---------------------------
    TF2 = [x for x in t.columns if x not in EVT and x not in
           {"track_idx", "cluster_idx", "k", "weight", "w_t", "wt_t", "ph_pt"}
           and not x.startswith("truth_") and t[x].dtype != object]
    truth_ev = c.groupby(EVT, sort=False)["t_truth"].first(); truth_ev.index.names = EVT
    tz2 = t[EVT + ["k", "time"] + [f for f in TF2 if f != "time"]].merge(
        truth_ev.rename("yt").rename_axis(EVT).reset_index(), on=EVT, how="inner")
    tz2.to_pickle(f"{W}/e2e_in.pkl")
    json.dump({"feats": TF2, "evt": EVT}, open(f"{W}/e2e_cfg.json", "w"))
    open(f"{W}/e2e_imp.py", "w").write(E2E_IMP)
    log(f"[4] e2e WeightNet: {len(TF2)} inputs, permutation importance (subprocess)")
    r = subprocess.run([a.python, f"{W}/e2e_imp.py", W, str(J)],
                       capture_output=True, text=True, timeout=10800)
    print(r.stdout, flush=True)
    if r.returncode == 0 and os.path.exists(f"{W}/e2e_imp.json"):
        e = json.load(open(f"{W}/e2e_imp.json"))
        R["e2e"] = {"n_features": len(TF2), "features": TF2,
                    "fold0_core_frac": e["base"], "top": e["perm"][:TOPN]}
    else:
        log(r.stderr[-1500:]); R["e2e"] = {"error": r.stderr[-400:]}

    # ---- 5. meta chooser ---------------------------------------------------
    from final_method_eval import etab
    c["p_old"] = c["p"]  # etab wants the column; not used for importances
    E, MFX = etab(c, t, t0f)
    dec = E.ok_base != E.ok_t0
    log(f"[5] meta chooser: {len(MFX)} inputs, {int(dec.sum()):,} disagreeing events "
        f"of {len(E):,}")
    m = xgbc(4, 10, 500, 0.05, J)
    m.fit(E.loc[dec, MFX].to_numpy(np.float32),
          (E.loc[dec, "ok_t0"] & ~E.loc[dec, "ok_base"]).astype(int))
    R["meta"] = {"n_features": len(MFX), "features": MFX,
                 "n_events": int(len(E)), "n_decisive": int(dec.sum()),
                 "decisive_frac": round(100 * float(dec.mean()), 2),
                 "t0_wins_frac": round(100 * float(
                     (E.loc[dec, "ok_t0"] & ~E.loc[dec, "ok_base"]).mean()), 2),
                 "top": imps(m, MFX)[:TOPN]}

    json.dump(R, open(a.out, "w"), indent=1)
    log(f"wrote {a.out}")

    # ---- signal-vs-background for the top 3 of each step -------------------
    # Each step gets the split its OWN model was trained against, which is not
    # the same question in each case -- a track-level HS/PU split for the
    # tagger, an in-window/out-of-window cluster split for the selector, and a
    # which-answer-wins split for the chooser. Comparing across steps therefore
    # compares separations of different quantities, deliberately.
    log("\ncomputing signal-vs-background distributions")
    NEED = ["time", "timeRes", "ph", "truth_is_hs", "nb_phzt", "nb_ph60",
            "closer_to_pu_than_pv", "cluster_closest_vtx_is_pv"]
    tv = t[EVT + NEED].merge(
        truth_ev.rename("t_truth_ev").rename_axis(EVT).reset_index(),
        on=EVT, how="left")
    tv["abs_dt_truth"] = (tv["time"] - tv["t_truth_ev"]).abs()
    yt = tv["truth_is_hs"].to_numpy().astype(bool)
    TL, BL = "hard scatter", "pileup"
    P = []
    for v in ("closer_to_pu_than_pv", "nb_phzt", "cluster_closest_vtx_is_pv"):
        P.append(panel(tv[v], yt, v, "1", "truth_is_hs", TL, BL))
    yc = c["ok_tagw"].to_numpy().astype(bool)
    for v in ("dt_tagw_agg", "dt_to_agg", "n_ph05"):
        P.append(panel(c[v], yc, v, "2", "ok_tagw",
                       "cluster lands within 60 ps", "cluster misses"))
    # step 3 has no fitted model, so "top 3" is the three quantities its
    # arithmetic consumes: the new weight, the old weight, and the residual
    # both are trying to estimate.
    for v in ("ph", "timeRes", "abs_dt_truth"):
        P.append(panel(tv[v], yt, v, "3", "truth_is_hs", TL, BL))
    for v in ("ph", "nb_ph60", "timeRes"):
        P.append(panel(tv[v], yt, v, "4", "truth_is_hs", TL, BL))
    dec = (E.ok_base != E.ok_t0).to_numpy()
    ye = (E["ok_t0"] & ~E["ok_base"]).to_numpy().astype(bool)[dec]
    for v in ("dt_tagw_agg", "abs_dt_base_t0", "p_gap"):
        P.append(panel(E[v].to_numpy()[dec], ye, v, "5", "switch helps",
                       "aggregate wins", "cluster wins"))
    # track-time outlier rate, quoted on the page as the reason for a median
    hs = tv.loc[yt, "abs_dt_truth"]
    D3 = {"panels": P,
          "hs_tail_gt200ps": round(100 * float((hs > 200).mean()), 2),
          "hs_median_abs_dt": round(float(hs.median()), 2),
          "hs_median_timeRes": round(float(tv.loc[yt, "timeRes"].median()), 2)}
    dout = a.out.replace(".json", "_dists.json")
    json.dump(D3, open(dout, "w"), indent=1)
    log(f"wrote {dout}  ({len(P)} panels)")


if __name__ == "__main__":
    main()
