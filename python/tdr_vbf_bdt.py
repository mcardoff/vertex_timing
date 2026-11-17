#!/usr/bin/env python3
"""A VBF H->inv vs Z+jets BDT at the TDR VBF preselection (Ariel's suggestion, vt-r2k.7), to see what
background remains after an analysis-like cut: tdr_vbf_plots.py --bdt draws the m_jj stacks after it.

Events: tdr_vbf's trees, preselected offline exactly as the stacks are (tdr_vbf_common.select with the
ITk-only R_pT tagger at 85%: the TDR's preselection). Inputs (tdr_vbf_common.pair_features): the VBF pair's
m_jj, |Delta eta|, |Delta phi|, both jets' pT and |eta|, eta1 * eta2, the surviving jets, and the recoil --
the truth MET for VBF, the Z pT for Z+jets (its leptons stand in for Z->nunu). The user's choice
(2026-11-04): kinematics + recoil. NB the VBF sample carries a generator MET > 75 GeV filter, so the recoil
separates strongly. No pileup-tagging or timing input: the BDT is the analysis, timing is studied after it.
Two-fold: events split by the parity of (file_idx * 1000003 + entry); each half is scored by the BDT trained
on the other. Training weights |w|, equal class totals, mean 1. Writes, per sample, one score per tree row
(NaN for events failing the preselection) to an .npz.

Run with ~/.venv-vtxml (xgboost; no ROOT):
  PYTHONNOUSERSITE=1 ~/.venv-vtxml/bin/python tdr_vbf_bdt.py --sig <vbf>_tdr_vbf.root --bkg <zjets>_tdr_vbf.root \\
      --wp <grid>_tdr_pujet_wp85.json --out <x>.npz
"""
import argparse

import numpy as np
import xgboost as xgb
from sklearn.metrics import roc_auc_score

import tdr_vbf_common as C

ap = argparse.ArgumentParser()
ap.add_argument("--sig", required=True)
ap.add_argument("--bkg", required=True)
ap.add_argument("--wp", required=True)
ap.add_argument("--out", required=True)
ap.add_argument("--seed", type=int, default=20261104)
args = ap.parse_args()

WP = C.load_wp(args.wp)
S, B = C.load(args.sig, WP), C.load(args.bkg, WP)
FS, FB = C.pair_features(S, WP, "rpt_itk"), C.pair_features(B, WP, "rpt_itk")
okS, okB = np.isfinite(FS["mjj"]), np.isfinite(FB["mjj"])
for name, F, ok in (("signal", FS, okS), ("background", FB, okB)):
    bad = [k for k in C.FEATURES if not np.all(np.isfinite(F[k][ok]))]
    if bad: raise SystemExit(f"{name}: non-finite inputs {bad} in preselected events")
XS = np.column_stack([FS[k] for k in C.FEATURES]).astype(np.float32)
XB = np.column_stack([FB[k] for k in C.FEATURES]).astype(np.float32)
par = lambda E: ((E["file_idx"].astype(np.int64) * 1000003 + E["entry"]) % 2).astype(bool)
pS, pB = par(S), par(B)
ws, wb = np.abs(S["weight"].astype(float)), np.abs(B["weight"].astype(float))
print(f"preselected (ITk-only R_pT at 85%): VBF {okS.sum()} (halves {(okS & ~pS).sum()} / {(okS & pS).sum()}), "
      f"Z+jets {okB.sum()} (halves {(okB & ~pB).sum()} / {(okB & pB).sum()})")

sS, sB = np.full(len(okS), np.nan), np.full(len(okB), np.nan)
imp = np.zeros(len(C.FEATURES))
aucs = []
for fold in (0, 1):
    trS, trB = okS & (pS == bool(fold)), okB & (pB == bool(fold))
    teS, teB = okS & (pS != bool(fold)), okB & (pB != bool(fold))
    n_tr = trS.sum() + trB.sum()
    w = np.concatenate([ws[trS] / ws[trS].sum(), wb[trB] / wb[trB].sum()]) * 0.5 * n_tr
    X = np.vstack([XS[trS], XB[trB]])
    y = np.concatenate([np.ones(trS.sum()), np.zeros(trB.sum())])
    # ~300 background events per fold, each weighing ~175 after the class balancing: a leaf must hold several
    # of them (min_child_weight is a hessian sum, w p (1 - p)), shallow trees, L2. Chosen by held-out AUC
    # (2026-11-04 scan): depth 3 / min_child_weight 5 gave train 0.966 vs held out 0.852; this 0.911 vs 0.866
    clf = xgb.XGBClassifier(n_estimators=150, max_depth=2, learning_rate=0.05, subsample=0.8, colsample_bytree=0.8,
                            min_child_weight=1000.0, reg_lambda=20.0, tree_method="hist", random_state=args.seed + fold, n_jobs=4)
    clf.fit(X, y, sample_weight=w)
    sS[teS] = clf.predict_proba(XS[teS])[:, 1]
    sB[teB] = clf.predict_proba(XB[teB])[:, 1]
    aucs.append(roc_auc_score(y, clf.predict_proba(X)[:, 1], sample_weight=w))
    imp += clf.feature_importances_
yt = np.concatenate([np.ones(okS.sum()), np.zeros(okB.sum())])
wt = np.concatenate([ws[okS] / ws[okS].sum(), wb[okB] / wb[okB].sum()])
auc = roc_auc_score(yt, np.concatenate([sS[okS], sB[okB]]), sample_weight=wt)
print(f"AUC: training {np.mean(aucs):.4f}, held out (both folds) {auc:.4f}")
print("importance (gain share, mean of the folds): " + ", ".join(f"{k} {v / 2:.2f}" for k, v in sorted(zip(C.FEATURES, imp), key=lambda x: -x[1])))
# background left at fixed signal efficiencies (held-out scores, weighted)
o = np.argsort(-sS[okS]); cs = np.cumsum(S["weight"][okS][o]) / S["weight"][okS].sum()
for eff in (0.8, 0.5, 0.2):
    thr = sS[okS][o][np.searchsorted(cs, eff)]
    m = okB & (sB >= thr)
    print(f"  signal efficiency {eff:.0%}: score > {thr:.4f}, Z+jets kept {B['weight'][m].sum() / B['weight'][okB].sum():.4f} "
          f"({int(m.sum())} MC events)")
np.savez_compressed(args.out, sig_score=sS, bkg_score=sB, features=np.array(C.FEATURES))
print(f"wrote {args.out}")
