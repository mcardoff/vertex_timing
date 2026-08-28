#!/usr/bin/env python3
"""
guarded_pooled.py -- two upgrades to the guarded metric, measured together and
separately: aggregate-t0 features in both stages, and POOLING the switch
training across all four mu=200 samples.

WHY POOL. The zjets switch model is starved -- ~36k events yield only ~2k
switch-helps positives -- and the ceiling measurement says decision quality is
what binds there: the top-2 challenger structural cap on zjets is 87.2% against
~66.7 achieved (on VBF the same gap is ~3 points). ttbar alone contributes ~24k
positives, so pooling multiplies zjets's effective switch training by ~10x.
Stage 1 and the aggregate-t0 tagger stay per-sample (in-domain, out-of-fold);
only the switch stage pools. Evaluation is per-sample and fold-honest.

AGGREGATE-t0 FEATURES (t0_agg, sum_ph, dt_to_agg): the prob-weighted-median
event time from the per-track HS tagger, i.e. event-level timing context that
no single cluster row carries. Measured worth ~+0.3-0.4 on zjets at moderate
operating points, null-not-negative on VBF.

RESULT (see results/guarded_pooled.log for full tables): pooling helps exactly
the small samples and costs nothing on the large ones. zjets max-net rises to
66.98 with FEWER losses (-898 vs -1057), its zero-loss row rises +45 -> +74,
and at matched ~25 lost pooled gains +538 vs own ~+380. dijet: same gains at
-656 vs -809 lost. vbf/ttbar: unchanged within noise (they dominate the pool).

zjets progression across this study, all vs TRKPTZ 61.87:
  frozen VBF transfer 65.97 -> in-domain switch 66.38 -> +agg 66.65 ->
  pooled+agg 66.98, with a strict zero-loss operating point (-0/+74) available.
"""
import numpy as np, pandas as pd, uproot, gc
import xgboost as xgb
EVT=["sample_id","file_idx","event_num"]
D="/Users/mcard/project/vertex_timing/training_new"
SAMPLES=["vbf","zjets","dijet","ttbar"]
rng=np.random.default_rng(0)
THR=(0.6,0.9,0.95,0.98)

def xgbc(depth,mcw):
    return xgb.XGBClassifier(n_estimators=500,learning_rate=0.05,max_depth=depth,
        subsample=0.8,colsample_bytree=0.8,min_child_weight=mcw,
        eval_metric="logloss",n_jobs=-1,tree_method="hist")

packs={}
FE=None
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
    # ---- aggregate t0 (in-domain tagger, out-of-fold) ----------------------
    t=uproot.open(F)["tracks"].arrays(library="pd")
    t=t[(t.time_valid>0)&(t.timeRes>0)].drop_duplicates(EVT+["track_idx"])
    t=t.merge(ev,on=EVT,how="left")
    TFE=[x for x in t.columns if x not in EVT and x not in
         {"track_idx","cluster_idx","k","weight"} and not x.startswith("truth_")
         and t[x].dtype!=object]
    t["ph"]=np.nan
    for kt in range(4):
        m=xgbc(7,10); m.set_params(n_estimators=300,learning_rate=0.06)
        tr=t[t.k!=kt]; m.fit(tr[TFE].to_numpy(np.float32),tr.truth_is_hs.to_numpy(int))
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
    sumw=pd.Series(tot,index=uniq)
    mi=pd.MultiIndex.from_frame(c[EVT])
    c["t0_agg"]=t0agg.reindex(mi).to_numpy()
    c["sum_ph"]=sumw.reindex(mi).to_numpy()
    c["dt_to_agg"]=(c["cluster_time"]-c["t0_agg"]).abs()
    del t,d; gc.collect()
    # ---- in-domain out-of-fold stage-1 p ------------------------------------
    FEA=FE+["t0_agg","sum_ph","dt_to_agg"]
    c["p"]=np.nan
    for kt in range(4):
        m1=xgbc(6,20); tr=c[c.k!=kt]
        m1.fit(tr[FEA].to_numpy(np.float32),tr.ok.to_numpy(int))
        msk=c.k==kt
        c.loc[msk,"p"]=m1.predict_proba(c.loc[msk,FEA].to_numpy(np.float32))[:,1]
    # ---- pairs ---------------------------------------------------------------
    g=c.groupby(EVT,sort=False); ib=g["trkptz_score"].idxmax()
    o=c.loc[ib].copy(); o["_eid"]=np.arange(len(o))
    emap=o.set_index(pd.MultiIndex.from_frame(o[EVT]))[["_eid"]]
    d2=c.loc[c.index.difference(ib)].copy()
    d2["_r"]=d2.groupby(EVT,sort=False)["p"].rank(ascending=False,method="first")
    ch=d2[d2._r<=2].copy()
    ch["_eid"]=emap.reindex(pd.MultiIndex.from_frame(ch[EVT]))["_eid"].to_numpy()
    ch=ch[np.isfinite(ch._eid)]; ch["_eid"]=ch._eid.astype(int)
    O=o.set_index("_eid")
    Xo=O.loc[ch._eid,FEA].to_numpy(np.float32); Xn=ch[FEA].to_numpy(np.float32)
    packs[s]=dict(
        X=np.hstack([Xo,Xn,Xn-Xo,O.loc[ch._eid,["p"]].to_numpy(np.float32),
                     ch[["p"]].to_numpy(np.float32)]),
        okO=O.loc[ch._eid,"ok"].to_numpy(bool), okN=ch["ok"].to_numpy(bool),
        eid=ch._eid.to_numpy(), kf=O["k"].to_numpy(int)[ch._eid.to_numpy()],
        pickok=O["ok"].to_numpy(bool), evk=O["k"].to_numpy(int))
    del c,o,d2,ch,O,Xo,Xn; gc.collect()
    print(f"{s}: pairs {len(packs[s]['X']):,}",flush=True)

for mode in ("own","pooled"):
    print(f"\n===== switch training: {mode} =====",flush=True)
    for kt in range(4):
        Xtr=[];ytr=[]
        srcs=SAMPLES if mode=="pooled" else None
        for s in SAMPLES:
            pk=packs[s]
            use=[s] if mode=="own" else SAMPLES
            # training pairs come from `use` samples, folds != kt
        # build once per fold
        for s2 in (SAMPLES if mode=="pooled" else []):
            pk=packs[s2]; trm=(pk["kf"]!=kt)&(pk["okO"]!=pk["okN"])
            Xtr.append(pk["X"][trm]); ytr.append((pk["okN"]&~pk["okO"])[trm])
        if mode=="own":
            pass
        else:
            m2=xgbc(5,10); m2.fit(np.vstack(Xtr),np.concatenate(ytr).astype(int))
        for s in SAMPLES:
            pk=packs[s]
            if mode=="own":
                trm=(pk["kf"]!=kt)&(pk["okO"]!=pk["okN"])
                m2s=xgbc(5,10); m2s.fit(pk["X"][trm],(pk["okN"]&~pk["okO"])[trm].astype(int))
                mm=m2s
            else:
                mm=m2
            tem=pk["kf"]==kt
            ps=mm.predict_proba(pk["X"][tem])[:,1]
            e_t,okN_t=pk["eid"][tem],pk["okN"][tem]
            order=np.lexsort((-ps,e_t)); first=np.unique(e_t[order],return_index=True)[1]
            be,bp,bok=e_t[order][first],ps[order][first],okN_t[order][first]
            st=pk.setdefault("res_"+mode,{x:[0,0] for x in THR})
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
            l,g_=pk["res_"+mode][x]
            line+=f"  [t={x}] {100*(trk-l+g_)/n:6.2f} (-{l}/+{g_})"
        print(line,flush=True)
