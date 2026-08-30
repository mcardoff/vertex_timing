#!/usr/bin/env python3
"""
final_method_eval.py -- the zjets_seventy method, run in-domain on every mu=200
sample plus the three mu=0 controls, for the cross-sample report.

Same recipe per sample as python/zjets_seventy.py (two-round tagger, trunc-tagw
base answer, e2e t0 answer, pooled two-way meta), with two hardening changes:
the round-2 neighbor self-merge is CHUNKED by event block (the unchunked form
OOM'd on ttbar at ~12 GB), and all subprocess/e2e state lives under --workdir
rather than /tmp (a /tmp purge killed the first evaluation run mid-flight and
took its log with it).

Sample-specific choices, deliberate:
  zjets   novbs-augmented training (94k events, key-hashed folds), 5-seed
          stage-1 -- exactly the shipped 70.00% configuration
  dijet   3-seed stage-1 (small sample)
  vbf, ttbar  1-seed stage-1 (large, stable; seed spread measured ~0.05 there)
  mu0     controls: no e2e (nothing to learn at ~99.9% TRKPTZ), formula t0,
          meta trained on the mu=200 pool applied FROZEN -- the question is
          only whether the method breaks where it has nothing to gain.
"""
import argparse, gc, json, os, subprocess, sys

import numpy as np
import pandas as pd
import uproot
import xgboost as xgb

EVT = ["sample_id", "file_idx", "event_num"]
KEY = EVT + ["cluster_idx"]
THR = (0.35, 0.4, 0.45, 0.5, 0.55, 0.6)


def log(m): print(m, flush=True)


def xgbc(depth, mcw, n=500, lr=0.05, seed=0):
    return xgb.XGBClassifier(n_estimators=n, learning_rate=lr, max_depth=depth,
                             subsample=0.8, colsample_bytree=0.8, min_child_weight=mcw,
                             eval_metric="logloss", n_jobs=-1, tree_method="hist",
                             random_state=seed)


def shrink(df):
    """uproot hands back float64/int64; halve the footprint. Keys stay int64
    (keyfold multiplies into uint64 range)."""
    for col in df.columns:
        if col in EVT: continue
        if df[col].dtype == np.float64: df[col] = df[col].astype(np.float32)
        elif df[col].dtype == np.int64 and col not in ("track_idx", "cluster_idx"):
            df[col] = df[col].astype(np.int32)
    return df


def keyfold(df):
    h = (df["file_idx"].to_numpy(np.int64) * 1000003 + df["event_num"].to_numpy(np.int64))
    h = (h.astype(np.uint64) * np.uint64(2654435761))
    return ((h >> np.uint64(13)) % np.uint64(4)).astype(int)


E2E_SUB = r'''
import json,sys,numpy as np,pandas as pd,torch,torch.nn as nn
wd=sys.argv[1]
cfg=json.load(open(f"{wd}/e2e_cfg.json")); EVT=cfg["evt"]; TF2=cfg["feats"]
tz2=pd.read_pickle(f"{wd}/e2e_in.pkl")
tz2=tz2.loc[:,~tz2.columns.duplicated()]
torch.set_num_threads(8)
class WNet(nn.Module):
    def __init__(self,n,h=96,it=5,floor=10.0):
        super().__init__(); self.it=it; self.floor=floor
        self.f=nn.Sequential(nn.Linear(n,h),nn.ReLU(),nn.Linear(h,h),nn.ReLU(),nn.Linear(h,1))
    def forward(self,x,t,e,n_ev):
        w=torch.nn.functional.softplus(self.f(x)).squeeze(-1)+1e-6
        num=torch.zeros(n_ev).index_add_(0,e,w*t); den=torch.zeros(n_ev).index_add_(0,e,w)
        t0=num/den.clamp_min(1e-9)
        for _ in range(self.it):
            r=(t-t0[e]).abs().detach().clamp_min(self.floor)
            u=w/r
            num=torch.zeros(n_ev).index_add_(0,e,u*t); den=torch.zeros(n_ev).index_add_(0,e,u)
            t0=num/den.clamp_min(1e-9)
        return t0
mu=tz2[TF2].mean(); sd=tz2[TF2].std().replace(0,1.0)
out=[]
def pack(df):
    cod,uq=pd.factorize(pd.MultiIndex.from_frame(df[EVT]),sort=False)
    if not isinstance(uq,pd.MultiIndex):
        uq=pd.MultiIndex.from_tuples([tuple(x) for x in uq],names=EVT)
    X=torch.from_numpy(((df[TF2]-mu)/sd).fillna(0).to_numpy(np.float32))
    return (X,torch.from_numpy(df["time"].to_numpy(np.float32)),
            torch.from_numpy(cod.astype(np.int64)),
            torch.from_numpy(df.groupby(EVT,sort=False)["yt"].first().to_numpy(np.float32)),
            len(uq),uq)
for kt in range(4):
    tr_=tz2[tz2.k!=kt].sort_values(EVT); te_=tz2[tz2.k==kt].sort_values(EVT)
    Xa,ta,ea,ya,na,_=pack(tr_); Xe,te2,ee,ye,ne_,uqe=pack(te_)
    torch.manual_seed(0)
    net=WNet(len(TF2)); opt=torch.optim.Adam(net.parameters(),lr=2e-3)
    best=None;bacc=-1
    for ep in range(1,201):
        tau=float(np.interp(ep,[1,200],[80.0,25.0]))
        net.train()
        t0b=net(Xa,ta,ea,na)
        loss=(1.0-torch.sigmoid((60.0-(t0b-ya).abs())/tau)).mean()
        opt.zero_grad(); loss.backward()
        torch.nn.utils.clip_grad_norm_(net.parameters(),1.0); opt.step()
        if ep%5==0 or ep==200:
            net.eval()
            with torch.no_grad():
                tv=net(Xe,te2,ee,ne_)
                acc=float(((tv-ye).abs()<60).float().mean())
            if acc>bacc: bacc=acc; best={k:v.clone() for k,v in net.state_dict().items()}
    net.load_state_dict(best)
    with torch.no_grad(): tv=net(Xe,te2,ee,ne_).numpy()
    out.append(pd.Series(tv,index=uqe))
    print(f"  e2e fold {kt}: val-best {100*bacc:.2f}%",flush=True)
pd.concat(out).to_pickle(f"{wd}/e2e_out.pkl")
'''


def build(files, n_seeds):
    cs = [shrink(uproot.open(f)["clusters"].arrays(library="pd")) for f in files]
    if len(cs) > 1:
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
    if len(ts) > 1:
        ts[1] = ts[1][~pd.MultiIndex.from_frame(ts[1][EVT]).isin(
            set(map(tuple, ts[0][EVT].drop_duplicates().to_numpy())))]
    t = pd.concat(ts, ignore_index=True).drop_duplicates(EVT + ["track_idx"]); del ts
    t["k"] = keyfold(t)
    TFE = [x for x in t.columns if x not in EVT and x not in
           {"track_idx", "cluster_idx", "k", "weight"} and not x.startswith("truth_")
           and t[x].dtype != object]
    t["ph"] = np.nan
    for kt in range(4):
        m = xgbc(7, 10, n=600, lr=0.05); tr = t[t.k != kt]
        m.fit(tr[TFE].to_numpy(np.float32), tr.truth_is_hs.to_numpy(int))
        msk = t.k == kt
        t.loc[msk, "ph"] = m.predict_proba(t.loc[msk, TFE].to_numpy(np.float32))[:, 1]
        del m, tr; gc.collect()
    # round 2: neighbor context, chunked by event block (ttbar OOM'd unchunked)
    nb = t[EVT + ["track_idx", "time", "z0", "ph"]].copy()
    ekeys = nb[EVT].drop_duplicates().reset_index(drop=True)
    NBs = []
    for i0 in range(0, len(ekeys), 60000):
        blk = ekeys.iloc[i0:i0 + 60000]
        sub = nb.merge(blk, on=EVT, how="inner")
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
    TFE2 = TFE + list(NB.columns); del NB; gc.collect()
    for kt in range(4):
        m = xgbc(7, 10, n=600, lr=0.05); tr = t[t.k != kt]
        m.fit(tr[TFE2].to_numpy(np.float32), tr.truth_is_hs.to_numpy(int))
        msk = t.k == kt
        t.loc[msk, "ph"] = m.predict_proba(t.loc[msk, TFE2].to_numpy(np.float32))[:, 1]
        del m, tr; gc.collect()
    # event t0: prob-weighted median, double-truncated
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
        gg = d2.groupby(EVT, sort=False)
        r = gg.apply(lambda x: (x.wp * x.time).sum() / x.wp.sum(), include_groups=False)
        r.index.names = EVT
        return r.reindex(med.index).fillna(med)
    t0f = trunc(trunc(med, 100), 60)
    mi = pd.MultiIndex.from_frame(c[EVT])
    c["t0_agg"] = med.reindex(mi).to_numpy()
    c["sum_ph"] = pd.Series(tot, index=med.index).reindex(mi).to_numpy()
    c["dt_to_agg"] = (c["cluster_time"] - c["t0_agg"]).abs()
    t["ph_pt"] = t["ph"] * t["pt"]; t["w_t"] = t["ph"] / t["timeRes"] ** 2; t["wt_t"] = t["w_t"] * t["time"]
    gcl = t.groupby(KEY, sort=False)
    cand = pd.DataFrame({"ph_sum": gcl["ph"].sum(), "ph_mean": gcl["ph"].mean(),
                         "ph_max": gcl["ph"].max(), "ph_pt_sum": gcl["ph_pt"].sum(),
                         "n_ph05": gcl["ph"].apply(lambda x: (x > 0.5).sum()),
                         "t_tagw": gcl["wt_t"].sum() / gcl["w_t"].sum()})
    tw = cand["t_tagw"]
    for win in (100.0, 60.0):
        tj = t.merge(tw.rename("tw").reset_index(), on=KEY, how="left")
        tj = tj[(tj["time"] - tj["tw"]).abs() < win]
        gg2 = tj.groupby(KEY, sort=False)
        upd = gg2["wt_t"].sum() / gg2["w_t"].sum()
        tw = upd.reindex(cand.index).fillna(tw)
    cand["t_tagw"] = tw
    cand["ph_pt_frac"] = cand["ph_pt_sum"] / gcl["pt"].sum()
    c = c.merge(cand.reset_index(), on=KEY, how="left")
    c["dt_tagw_agg"] = (c["t_tagw"] - c["t0_agg"]).abs()
    tw_eff = np.where(np.isnan(c["t_tagw"]), c["cluster_time"], c["t_tagw"])
    c["ok_tagw"] = np.abs(tw_eff - c["t_truth"]) < 60
    drop = {"ok", "ok_tagw", "k", "p", "p_old", "delta_t", "weight", "cluster_idx", "t_truth"}
    FEA = [x for x in c.columns if x not in EVT and x not in drop
           and not x.startswith("truth_") and c[x].dtype != object]
    for lab, col in (("ok_tagw", "p"), ("ok", "p_old")):
        ps = np.zeros(len(c))
        for sd in range(n_seeds):
            pc = np.full(len(c), np.nan)
            for kt in range(4):
                m1 = xgbc(6, 20, seed=sd); tr = c[c.k != kt]
                m1.fit(tr[FEA].to_numpy(np.float32), tr[lab].to_numpy(int))
                msk = (c.k == kt).to_numpy()
                pc[msk] = m1.predict_proba(c.loc[msk, FEA].to_numpy(np.float32))[:, 1]
                del m1, tr; gc.collect()
            ps += pc
        c[col] = ps / n_seeds
    return c, t, t0f


def run_e2e(t, c, wd, py):
    TF2 = [x for x in t.columns if x not in EVT and x not in
           {"track_idx", "cluster_idx", "k", "weight", "w_t", "wt_t", "ph_pt"}
           and not x.startswith("truth_") and t[x].dtype != object]
    truth_ev = c.groupby(EVT, sort=False)["t_truth"].first(); truth_ev.index.names = EVT
    tz2 = t[EVT + ["k", "time"] + [f for f in TF2 if f != "time"]].merge(
        truth_ev.rename("yt").rename_axis(EVT).reset_index(), on=EVT, how="inner")
    tz2.to_pickle(f"{wd}/e2e_in.pkl")
    json.dump({"feats": TF2, "evt": EVT}, open(f"{wd}/e2e_cfg.json", "w"))
    if os.path.exists(f"{wd}/e2e_out.pkl"): os.remove(f"{wd}/e2e_out.pkl")
    open(f"{wd}/e2e_sub.py", "w").write(E2E_SUB)
    r = subprocess.run([py, f"{wd}/e2e_sub.py", wd], capture_output=True, text=True, timeout=10800)
    print(r.stdout, flush=True)
    if r.returncode != 0:
        print(r.stderr[-1500:], flush=True); return None
    s = pd.read_pickle(f"{wd}/e2e_out.pkl"); s.index.names = EVT
    return s.reindex(truth_ev.index)


def etab(c, t, t0):
    g = c.groupby(EVT, sort=False)
    truth = g["t_truth"].first(); truth.index.names = EVT
    ip = g["p"].idxmax()
    pick = c.loc[ip].set_index(pd.MultiIndex.from_frame(c.loc[ip, EVT]))
    trkpick = c.loc[g["trkptz_score"].idxmax()]
    c2 = c.loc[c.index.difference(ip)]
    p2 = c2.groupby(EVT, sort=False)["p"].max(); p2.index.names = EVT
    cc = c.merge(t0.rename("t0ev").rename_axis(EVT).reset_index(), on=EVT, how="left")
    cc["dt_c_t0"] = (cc["cluster_time"] - cc["t0ev"]).abs()
    nearc = cc[cc["dt_c_t0"] < 60].groupby(EVT, sort=False)
    corro = pd.DataFrame({"maxp_near_t0": nearc["p"].max(), "n_cl_near_t0": nearc.size(),
                          "sump_near_t0": nearc["p"].sum()})
    corro.index.names = EVT
    E = pd.DataFrame(index=pick.index)
    E["k"] = pick["k"].to_numpy()
    E["truth"] = truth.reindex(E.index).to_numpy()
    tw = pick["t_tagw"].to_numpy(); ctm = pick["cluster_time"].to_numpy()
    E["base_t"] = np.where(np.isnan(tw), ctm, tw)
    E["ok_base"] = (np.abs(E["base_t"] - E["truth"]) < 60)
    E["ok_ctime"] = (np.abs(ctm - E["truth"]) < 60)
    E["ok_trkptz"] = trkpick["ok"].to_numpy()
    E["t0"] = t0.reindex(E.index).to_numpy()
    E["ok_t0"] = ((E["t0"] - E["truth"]).abs() < 60).fillna(False)
    MF = ["p", "trkptz_score", "waves_score", "n_tracks", "sumpt", "delta_z",
          "cluster_time_sigma", "frac_valid_time", "time_chi2_ndf", "max_abs_tpull",
          "n_clusters", "t0_agg", "sum_ph", "dt_to_agg", "ph_sum", "ph_mean",
          "ph_pt_frac", "n_ph05", "t_tagw", "dt_tagw_agg"]
    for f in MF: E[f] = pick[f].to_numpy()
    E["p_old"] = pick["p_old"].to_numpy()
    E["p2"] = p2.reindex(E.index).fillna(0).to_numpy()
    E["p_gap"] = E["p"] - E["p2"]
    E["abs_dt_base_t0"] = np.abs(E["base_t"] - E["t0"])
    E["abs_dt_base_ctime"] = np.abs(E["base_t"] - ctm)
    for f in corro.columns: E[f] = corro[f].reindex(E.index).fillna(0).to_numpy()
    if "nb_ph60" in t.columns:
        tpk = t.merge(pick[EVT + ["cluster_idx"]].reset_index(drop=True),
                      on=EVT + ["cluster_idx"], how="inner")
        gp = tpk.groupby(EVT, sort=False)
        E["pick_nb_ph60_mean"] = gp["nb_ph60"].mean().reindex(E.index).to_numpy()
        E["pick_nb_ph60_max"] = gp["nb_ph60"].max().reindex(E.index).to_numpy()
    else:
        E["pick_nb_ph60_mean"] = np.nan; E["pick_nb_ph60_max"] = np.nan
    tt = t[EVT + ["time", "ph"]].merge(E[["t0", "base_t"]].rename_axis(EVT).reset_index(),
                                       on=EVT, how="inner")
    for lab, col in (("t0", "t0"), ("base", "base_t")):
        nr = (tt["time"] - tt[col]).abs() < 60
        E[f"n_near_{lab}"] = nr.groupby([tt[k] for k in EVT]).sum().reindex(E.index).fillna(0).to_numpy()
        E[f"ph_near_{lab}"] = (tt["ph"] * nr).groupby([tt[k] for k in EVT]).sum().reindex(E.index).fillna(0).to_numpy()
    MFX = MF + ["p_old", "p2", "p_gap", "abs_dt_base_t0", "abs_dt_base_ctime",
                "maxp_near_t0", "n_cl_near_t0", "sump_near_t0", "n_near_t0",
                "ph_near_t0", "n_near_base", "ph_near_base",
                "pick_nb_ph60_mean", "pick_nb_ph60_max"]
    return E, MFX


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--input-dir", required=True)
    ap.add_argument("--workdir", required=True)
    ap.add_argument("--python", default=sys.executable)
    ap.add_argument("--out", default="")
    args = ap.parse_args()
    D = args.input_dir; os.makedirs(args.workdir, exist_ok=True)
    CFG = {
        "vbf":   dict(files=[f"{D}/vbf_mjj500p0_training.root"], seeds=1, e2e=True),
        "zjets": dict(files=[f"{D}/zjets_novbs_training.root",
                             f"{D}/zjets_mjj500p0_training.root"], seeds=5, e2e=True),
        "dijet": dict(files=[f"{D}/dijet_mjj500p0_training.root"], seeds=3, e2e=True),
        "ttbar": dict(files=[f"{D}/ttbar_mjj500p0_training.root"], seeds=1, e2e=True),
        "vbf_mu0":     dict(files=[f"{D}/vbf_mu0_mjj500p0_training.root"], seeds=1, e2e=False),
        "zeejets_mu0": dict(files=[f"{D}/zeejets_mu0_mjj500p0_training.root"], seeds=1, e2e=False),
        "ttbar_mu0":   dict(files=[f"{D}/ttbar_mu0_mjj500p0_training.root"], seeds=1, e2e=False),
    }
    ES, comp = {}, {}
    MFX = None
    for s, cfg in CFG.items():
        log(f"\n===== {s} =====")
        c, t, t0f = build(cfg["files"], cfg["seeds"])
        truth = c.groupby(EVT, sort=False)["t_truth"].first(); truth.index.names = EVT
        if s == "zjets":
            can = uproot.open(f"{D}/zjets_mjj500p0_training.root")["clusters"].arrays(EVT, library="pd")
            cm = truth.index.isin(pd.MultiIndex.from_frame(can.drop_duplicates()))
        else:
            cm = np.ones(len(truth), bool)
        okf = ((t0f.reindex(truth.index) - truth).abs() < 60).fillna(False)
        t0 = t0f; lab_f = 100 * okf[cm].mean(); lab_e = float("nan")
        if cfg["e2e"]:
            e2 = run_e2e(t, c, args.workdir, args.python)
            if e2 is not None:
                oke = ((e2.reindex(truth.index) - truth).abs() < 60).fillna(False)
                lab_e = 100 * oke[cm].mean()
                if lab_e > lab_f: t0 = e2
        E, MFX = etab(c, t, t0)
        E["is_can"] = cm
        ES[s] = E
        comp[s] = dict(n=int(cm.sum()),
                       trkptz=100 * E.loc[cm, "ok_trkptz"].mean(),
                       ctime=100 * E.loc[cm, "ok_ctime"].mean(),
                       base=100 * E.loc[cm, "ok_base"].mean(),
                       t0_formula=lab_f, t0_e2e=lab_e,
                       t0_used=100 * E.loc[cm, "ok_t0"].mean(),
                       union=100 * (E.loc[cm, "ok_base"] | E.loc[cm, "ok_t0"]).mean())
        log(f"{s}: TRKPTZ {comp[s]['trkptz']:.2f}  ctime {comp[s]['ctime']:.2f}  "
            f"base {comp[s]['base']:.2f}  t0f {lab_f:.2f}  t0e2e {lab_e:.2f}  "
            f"union {comp[s]['union']:.2f}")
        del c, t; gc.collect()

    MU200 = ["vbf", "zjets", "dijet", "ttbar"]
    ALL = pd.concat([ES[s] for s in MU200])
    log("\n===== pooled meta (mu200), per-sample evaluation =====")
    res = {s: {x: 0 for x in THR} for s in ES}
    for kt in range(4):
        trm = ALL["k"] != kt; dec = trm & (ALL.ok_base != ALL.ok_t0)
        mods = []
        for sd in range(3):
            mm = xgbc(4, 10, n=500, seed=sd)
            mm.fit(ALL.loc[dec, MFX].to_numpy(np.float32),
                   (ALL.loc[dec, "ok_t0"] & ~ALL.loc[dec, "ok_base"]).astype(int))
            mods.append(mm)
        for s, E in ES.items():
            tem = (E["k"] == kt) & E["is_can"]
            if not tem.any(): continue
            ps = np.mean([mm.predict_proba(E.loc[tem, MFX].to_numpy(np.float32))[:, 1]
                          for mm in mods], axis=0)
            for x in THR:
                sw = ps > x
                ok = np.where(sw, E.loc[tem, "ok_t0"], E.loc[tem, "ok_base"])
                res[s][x] += int(ok.sum())
    table = {}
    for s in ES:
        n = comp[s]["n"]
        table[s] = {str(x): 100 * res[s][x] / n for x in THR}
        line = f"  {s:12s} TRKPTZ {comp[s]['trkptz']:6.2f}"
        for x in THR: line += f"  [{x}] {table[s][str(x)]:6.2f}"
        log(line)
    if args.out:
        json.dump({"components": comp, "final": table}, open(args.out, "w"), indent=1)
        log(f"wrote {args.out}")


if __name__ == "__main__":
    main()
