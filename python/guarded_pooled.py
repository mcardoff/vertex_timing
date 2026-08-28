#!/usr/bin/env python3
"""
guarded_pooled.py -- the full guarded-metric configuration: per-sample stage-1
with aggregate-t0 AND per-candidate tagger aggregates, cross-candidate time
consistency in the pair features, switch training pooled across all four mu=200
samples. Supersedes the earlier version of this file (git history has it).

The two feature blocks this version adds attack the zjets decision-quality gap
(top-2 structural cap 87.2% vs ~67 achieved) with exactly what the failure
gallery said was missing:

  per-candidate tagger aggregates   ph_sum/ph_mean/ph_max/ph_pt_sum/n_ph05/
        (HS content per candidate)  ph_pt_frac/t_tagw/dt_tagw_agg -- the
                                    per-track HS tagger (out-of-fold) summed
                                    over each candidate's own tracks. The wrong
                                    picks are pure pileup; this measures that
                                    reco-only. Worth +0.47 on zjets alone.
  cross-candidate time consistency  fraction of the PICK's timed tracks within
                                    60 ps of the CHALLENGER's time, and vice
                                    versa. If the pick's own tracks agree
                                    better with the challenger's time, the
                                    pick's time is an outlier artifact. +0.11.

Cross-fitted final table (ALL events, every verdict out-of-fold):

  sample TRKPTZ    t=0.6              t=0.9            t=0.95          t=0.98
  vbf    90.45   94.02 (-2817/+21335) 93.01 (-326)   92.52 (-133)   91.88 (-34)
  zjets  61.87   67.49 (-987/+3040)   63.50 (-34)    62.73 (-8)     62.17 (-0/+111)
  dijet  86.88   90.12 (-695/+3295)   88.56 (-48)    88.06 (-13)    87.54 (-4)
  ttbar  87.08   90.61 (-4877/+24723) 88.96 (-385)   88.35 (-107)   87.76 (-29)

zjets progression across the study, all vs TRKPTZ 61.87: frozen VBF transfer
65.97 -> in-domain switch 66.38 -> +aggregate-t0 66.65 -> pooled 66.98 ->
+candidate-HS/+cross-time pooled 67.49, now ABOVE in-domain Deep Sets (67.3),
with a zero-loss operating point of -0/+111.
"""
import numpy as np, pandas as pd, uproot, gc
import xgboost as xgb
EVT=["sample_id","file_idx","event_num"]; KEY=EVT+["cluster_idx"]
D="/Users/mcard/project/vertex_timing/training_new"
SAMPLES=["vbf","zjets","dijet","ttbar"]
THR=(0.6,0.9,0.95,0.98)
rng=np.random.default_rng(0)
def xgbc(depth,mcw,n=500,lr=0.05):
    return xgb.XGBClassifier(n_estimators=n,learning_rate=lr,max_depth=depth,
        subsample=0.8,colsample_bytree=0.8,min_child_weight=mcw,
        eval_metric="logloss",n_jobs=-1,tree_method="hist")

packs={}; FE=None
for s in SAMPLES:
    F=f"{D}/{s}_mjj500p0_training.root"
    c=uproot.open(F)["clusters"].arrays(library="pd")
    c["ok"]=c["delta_t"].abs()<60
    ev=c[EVT].drop_duplicates().reset_index(drop=True); ev["k"]=rng.integers(0,4,len(ev))
    c=c.merge(ev,on=EVT,how="left")
    if FE is None:
        drop={"ok","k","p","delta_t","weight","cluster_idx"}
        FE=[x for x in c.columns if x not in EVT and x not in drop
            and not x.startswith("truth_") and c[x].dtype!=object]
    t=uproot.open(F)["tracks"].arrays(library="pd")
    t=t[(t.time_valid>0)&(t.timeRes>0)].drop_duplicates(EVT+["track_idx"])
    t=t.merge(ev,on=EVT,how="left")
    TFE=[x for x in t.columns if x not in EVT and x not in
         {"track_idx","cluster_idx","k","weight"} and not x.startswith("truth_")
         and t[x].dtype!=object]
    t["ph"]=np.nan
    for kt in range(4):
        m=xgbc(7,10,n=300,lr=0.06); tr=t[t.k!=kt]
        m.fit(tr[TFE].to_numpy(np.float32),tr.truth_is_hs.to_numpy(int))
        msk=t.k==kt
        t.loc[msk,"ph"]=m.predict_proba(t.loc[msk,TFE].to_numpy(np.float32))[:,1]
    d=t.sort_values(EVT+["time"])
    codes,uniq=pd.factorize(pd.MultiIndex.from_frame(d[EVT]),sort=False)
    w=d["ph"].to_numpy(); tm=d["time"].to_numpy()
    tot=np.bincount(codes,weights=w)
    start=np.r_[0,np.flatnonzero(np.diff(codes))+1]
    cw=np.cumsum(w); off=np.zeros(len(d)); off[start[1:]]=cw[start[1:]-1]
    cw=cw-np.maximum.accumulate(off)
    hit=cw>=0.5*tot[codes]; idx=np.arange(len(d))
    firsts=np.minimum.reduceat(np.where(hit,idx,len(d)+1),start)
    t0agg=pd.Series(tm[np.clip(firsts,0,len(d)-1)],index=uniq)
    mi=pd.MultiIndex.from_frame(c[EVT])
    c["t0_agg"]=t0agg.reindex(mi).to_numpy()
    c["sum_ph"]=pd.Series(tot,index=uniq).reindex(mi).to_numpy()
    c["dt_to_agg"]=(c["cluster_time"]-c["t0_agg"]).abs()
    t["ph_pt"]=t["ph"]*t["pt"]; t["w_t"]=t["ph"]/t["timeRes"]**2; t["wt_t"]=t["w_t"]*t["time"]
    gcl=t.groupby(KEY,sort=False)
    cand=pd.DataFrame({"ph_sum":gcl["ph"].sum(),"ph_mean":gcl["ph"].mean(),
        "ph_max":gcl["ph"].max(),"ph_pt_sum":gcl["ph_pt"].sum(),
        "n_ph05":gcl["ph"].apply(lambda x:(x>0.5).sum()),
        "t_tagw":gcl["wt_t"].sum()/gcl["w_t"].sum()})
    cand["ph_pt_frac"]=cand["ph_pt_sum"]/gcl["pt"].sum()
    c=c.merge(cand.reset_index(),on=KEY,how="left")
    c["dt_tagw_agg"]=(c["t_tagw"]-c["t0_agg"]).abs()
    FEA=FE+["t0_agg","sum_ph","dt_to_agg","ph_sum","ph_mean","ph_max","ph_pt_sum",
            "n_ph05","ph_pt_frac","t_tagw","dt_tagw_agg"]
    c["p"]=np.nan
    for kt in range(4):
        m1=xgbc(6,20); tr=c[c.k!=kt]
        m1.fit(tr[FEA].to_numpy(np.float32),tr.ok.to_numpy(int))
        msk=c.k==kt
        c.loc[msk,"p"]=m1.predict_proba(c.loc[msk,FEA].to_numpy(np.float32))[:,1]
    g=c.groupby(EVT,sort=False); ib=g["trkptz_score"].idxmax()
    o=c.loc[ib].copy(); o["_eid"]=np.arange(len(o))
    emap=o.set_index(pd.MultiIndex.from_frame(o[EVT]))[["_eid"]]
    d2=c.loc[c.index.difference(ib)].copy()
    d2["_r"]=d2.groupby(EVT,sort=False)["p"].rank(ascending=False,method="first")
    ch=d2[d2._r<=2].copy()
    ch["_eid"]=emap.reindex(pd.MultiIndex.from_frame(ch[EVT]))["_eid"].to_numpy()
    ch=ch[np.isfinite(ch._eid)]; ch["_eid"]=ch._eid.astype(int)
    O=o.set_index("_eid")
    # vectorized cross-candidate time consistency
    tt=t[EVT+["cluster_idx","time"]]
    def xfrac(keys,other_t):
        j=keys.reset_index(drop=True).copy(); j["_pid"]=np.arange(len(j)); j["_ot"]=other_t.to_numpy()
        m=j.merge(tt,on=EVT+["cluster_idx"],how="left")
        good=(m["time"].notna()&((m["time"]-m["_ot"]).abs()<60)).to_numpy(float)
        cnt=np.bincount(m["_pid"],minlength=len(j)).astype(float)
        hits=np.bincount(m["_pid"],weights=good,minlength=len(j))
        with np.errstate(invalid="ignore",divide="ignore"):
            f=np.where(cnt>0,hits/cnt,0.0)
        return f.astype(np.float32)
    pk_key=O.loc[ch._eid,EVT+["cluster_idx"]].reset_index(drop=True)
    fA=xfrac(pk_key,ch["cluster_time"])
    fB=xfrac(ch[EVT+["cluster_idx"]],O.loc[ch._eid,"cluster_time"])
    Xo=O.loc[ch._eid,FEA].to_numpy(np.float32); Xn=ch[FEA].to_numpy(np.float32)
    packs[s]=dict(
        X=np.hstack([Xo,Xn,Xn-Xo,O.loc[ch._eid,["p"]].to_numpy(np.float32),
                     ch[["p"]].to_numpy(np.float32),
                     np.stack([fA,fB,fB-fA],axis=1)]),
        okO=O.loc[ch._eid,"ok"].to_numpy(bool), okN=ch["ok"].to_numpy(bool),
        eid=ch._eid.to_numpy(), kf=O["k"].to_numpy(int)[ch._eid.to_numpy()],
        pickok=O["ok"].to_numpy(bool), evk=O["k"].to_numpy(int))
    del c,t,o,d2,ch,O,Xo,Xn; gc.collect()
    print(f"{s}: pack built, {len(packs[s]['X']):,} pairs",flush=True)

print("\n===== pooled switch, full feature set =====",flush=True)
for kt in range(4):
    Xtr=[];ytr=[]
    for s2 in SAMPLES:
        pk=packs[s2]; trm=(pk["kf"]!=kt)&(pk["okO"]!=pk["okN"])
        Xtr.append(pk["X"][trm]); ytr.append((pk["okN"]&~pk["okO"])[trm])
    m2=xgbc(5,10); m2.fit(np.vstack(Xtr),np.concatenate(ytr).astype(int))
    for s in SAMPLES:
        pk=packs[s]; tem=pk["kf"]==kt
        ps=m2.predict_proba(pk["X"][tem])[:,1]
        e_t,okN_t=pk["eid"][tem],pk["okN"][tem]
        order=np.lexsort((-ps,e_t)); first=np.unique(e_t[order],return_index=True)[1]
        be,bp,bok=e_t[order][first],ps[order][first],okN_t[order][first]
        st=pk.setdefault("res",{x:[0,0] for x in THR})
        evk=pk["evk"]==kt
        for x in THR:
            ok=pk["pickok"].copy(); sw=bp>x; ok[be[sw]]=bok[sw]
            st[x][0]+=int((pk["pickok"][evk]&~ok[evk]).sum())
            st[x][1]+=int((~pk["pickok"][evk]&ok[evk]).sum())
    print(f"  fold {kt} done",flush=True)
for s in SAMPLES:
    pk=packs[s]; n=len(pk["pickok"]); trk=int(pk["pickok"].sum())
    line=f"  {s:6s} TRKPTZ {100*trk/n:6.2f}%"
    for x in THR:
        l,g_=pk["res"][x]
        line+=f"  [t={x}] {100*(trk-l+g_)/n:6.2f} (-{l}/+{g_})"
    print(line,flush=True)
