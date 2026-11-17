#!/usr/bin/env python3
"""A VBF H->inv vs Z+jets BDT on the 2022 PUB note's loose selection (util/tdr_fwdpair.cxx's events: the two
leading jets in 2.4 < |eta| < 3.8, opposite eta, pT > 50 / 30 GeV) -- the user's alternative (2026-11-04) to the
TDR preselection + recoil BDT (tdr_vbf_bdt.py), whose recoil left too little Z+jets. No MET: the inputs are the
two leading jets only -- m_jj, |Delta eta|, |Delta phi|, both pT, both |eta|. No pileup-tagging or timing input.
Every event of each sample is used (signal: all VBF events, background: all Z+jets events). Two-fold: events
split by the parity of (file_idx * 1000003 + entry), each half scored by the BDT trained on the other; weights
|w| class-balanced, mean 1. Settings from a held-out AUC scan (2026-11-04): train 0.896, held out 0.884.
Writes one score per tree row to an .npz for tdr_fwdpair_kinplots.py.

  PYTHONNOUSERSITE=1 ~/.venv-vtxml/bin/python tdr_fwdpair_kinbdt.py --sig <vbf>_tdr_fwdpair.root \\
      --bkg <zjets>_tdr_fwdpair.root --out <x>.npz
"""
import argparse

import numpy as np
import uproot
import xgboost as xgb
from sklearn.metrics import roc_auc_score

ap = argparse.ArgumentParser()
ap.add_argument("--sig", required=True)
ap.add_argument("--bkg", required=True)
ap.add_argument("--out", required=True)
ap.add_argument("--seed", type=int, default=20261104)
args = ap.parse_args()

FEATURES = ["mjj", "deta", "dphi", "pt1", "pt2", "aeta1", "aeta2"]


def load(path):
    E = uproot.open(path)["events"].arrays(library="np")
    dphi = np.abs((E["j1_phi"] - E["j2_phi"] + np.pi) % (2 * np.pi) - np.pi)
    X = np.column_stack([E["mjj"], E["deta"], dphi, E["j1_pt"], E["j2_pt"], np.abs(E["j1_eta"]), np.abs(E["j2_eta"])])
    return E, X.astype(np.float32)


S, XS = load(args.sig)
B, XB = load(args.bkg)
par = lambda E: ((E["file_idx"].astype(np.int64) * 1000003 + E["entry"]) % 2).astype(bool)
pS, pB = par(S), par(B)
ws, wb = np.abs(S["weight"].astype(float)), np.abs(B["weight"].astype(float))
print(f"VBF {len(ws)} events (halves {(~pS).sum()} / {pS.sum()}), Z+jets {len(wb)} (halves {(~pB).sum()} / {pB.sum()})")
sS, sB = np.full(len(ws), np.nan), np.full(len(wb), np.nan)
imp, aucs = np.zeros(len(FEATURES)), []
for fold in (0, 1):
    a, b = pS == bool(fold), pB == bool(fold)
    n = a.sum() + b.sum()
    w = np.concatenate([ws[a] / ws[a].sum(), wb[b] / wb[b].sum()]) * 0.5 * n
    X = np.vstack([XS[a], XB[b]])
    y = np.concatenate([np.ones(a.sum()), np.zeros(b.sum())])
    clf = xgb.XGBClassifier(n_estimators=200, max_depth=3, learning_rate=0.05, subsample=0.8, colsample_bytree=0.8,
                            min_child_weight=200.0, reg_lambda=10.0, tree_method="hist", random_state=args.seed + fold, n_jobs=4)
    clf.fit(X, y, sample_weight=w)
    sS[~a] = clf.predict_proba(XS[~a])[:, 1]
    sB[~b] = clf.predict_proba(XB[~b])[:, 1]
    aucs.append(roc_auc_score(y, clf.predict_proba(X)[:, 1], sample_weight=w))
    imp += clf.feature_importances_
auc = roc_auc_score(np.r_[np.ones(len(ws)), np.zeros(len(wb))], np.r_[sS, sB], sample_weight=np.r_[ws / ws.sum(), wb / wb.sum()])
print(f"AUC: training {np.mean(aucs):.4f}, held out (both folds) {auc:.4f}")
print("importance (mean of the folds): " + ", ".join(f"{k} {v / 2:.2f}" for k, v in sorted(zip(FEATURES, imp), key=lambda x: -x[1])))
np.savez_compressed(args.out, sig_score=sS, bkg_score=sB, features=np.array(FEATURES))
print(f"wrote {args.out}")
