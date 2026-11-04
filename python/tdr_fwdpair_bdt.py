#!/usr/bin/env python3
"""ATL-HGTD-PUB-2022-001 Sec. 6, the BDTs (its Table 1): train on half of the events, score the other half.

Reads util/tdr_fwdpair.cxx's trees. Signal: the VBF sample's events whose two jets are both hard-scatter;
background: the background sample's events with at least one pileup jet (the note's Sec. 3). One BDT per
row, each on its own inputs:
  itk         JVF(j1), JVF(j2)                                    the note's "BDT ITk"
  t0_<s>      JVF(j1, t0), JVF(j2, t0) for s = hgtd (Athena: the note's reconstructed t0), trkptz, waves,
              truth (the note's t0^true)
  dtime       JVF(j1), JVF(j2), |t_j1 - t_j2| / sqrt(s_j1^2 + s_j2^2), |t_lead1 - t_lead2| / sqrt(...)
              -- the note's "BDT ITk+HGTD(Delta time)", no t0: jet time = weighted mean over the jet's
              timed tracks, s_j^2 their weighted variance; missing (NaN) when a jet has no timed track,
              so the BDT falls back on the JVFs as the note does
  dtime_hs    beyond the note: as dtime, but the jet time over the jet's PV0-associated timed tracks
              (z only, still no t0) with the error of the mean -- at mu = 200 the all-track spread is
              pileup-dominated and the note's pull separates nothing (notes/tdr_presentation.md)
Split: events with an even (file_idx * 1000003 + entry + jz) train, odd ones are scored. Training weights:
|w| normalised to equal totals per class (and, for dijet, per JZ slice: the slices have no cross-sections,
notes/tdr_presentation.md). Scores, signed weights and labels go to an .npz for tdr_fwdpair_plots.py.

Run with ~/.venv-vtxml (xgboost; no ROOT needed):
  PYTHONNOUSERSITE=1 ~/.venv-vtxml/bin/python tdr_fwdpair_bdt.py --sig <vbf>.root --bkg <zjets>.root --out <x>.npz
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

T0S = ("hgtd", "trkptz", "waves", "truth")


def load(path):
    return uproot.open(path)["events"].arrays(library="np")


def pull(a, va, b, vb, ok):
    with np.errstate(invalid="ignore", divide="ignore"):
        p = np.abs(a - b) / np.sqrt(va + vb)
    return np.where(ok & np.isfinite(p), p, np.nan)


def features(E):
    f = {"itk": [E["j1_jvf_itk"], E["j2_jvf_itk"]]}
    for s in T0S:
        f["t0_" + s] = [E[f"j1_jvf_t0_{s}"], E[f"j2_jvf_t0_{s}"]]
    timed = (E["j1_ntimed"] > 0) & (E["j2_ntimed"] > 0)
    p_jet = pull(E["j1_t_jet"], E["j1_var_spread"], E["j2_t_jet"], E["j2_var_spread"], timed)
    p_lead = pull(E["j1_t_lead"], E["j1_sig_lead"] ** 2, E["j2_t_lead"], E["j2_sig_lead"] ** 2, timed)
    timed_hs = (E["j1_ntimed_hs"] > 0) & (E["j2_ntimed_hs"] > 0)
    p_hs = pull(E["j1_t_hs"], E["j1_var_hs_mean"], E["j2_t_hs"], E["j2_var_hs_mean"], timed_hs)
    f["dtime"] = [E["j1_jvf_itk"], E["j2_jvf_itk"], p_jet, p_lead]
    f["dtime_hs"] = [E["j1_jvf_itk"], E["j2_jvf_itk"], p_hs, p_lead]
    return {k: np.column_stack(v).astype(np.float32) for k, v in f.items()}


S, B = load(args.sig), load(args.bkg)
sig = (S["j1_hs"] == 1) & (S["j2_hs"] == 1)
bkg = (B["j1_pu"] == 1) | (B["j2_pu"] == 1)
FS, FB = features(S), features(B)
odd = lambda E: ((E["file_idx"].astype(np.int64) * 1000003 + E["entry"] + E["jz"]) % 2).astype(bool)
oS, oB = odd(S), odd(B)
is_dijet = bool(np.any(B["jz"] >= 0))

# training weights: |w|, equal class totals; dijet slices equal within the background
ws = np.abs(S["weight"].astype(float))
wb = np.abs(B["weight"].astype(float))
if is_dijet:
    for z in np.unique(B["jz"][bkg]):
        m = bkg & (B["jz"] == z)
        wb[m] /= wb[m & ~oB].sum() or 1.0
trS, trB = sig & ~oS, bkg & ~oB
teS, teB = sig & oS, bkg & oB
# equal class totals, mean weight 1 (xgboost's min_child_weight is in these units: weights summing to 1
# per class forbid every split and leave a constant classifier, AUC 0.5)
n_tr = trS.sum() + trB.sum()
ws_tr = ws[trS] / ws[trS].sum() * 0.5 * n_tr
wb_tr = wb[trB] / wb[trB].sum() * 0.5 * n_tr
print(f"signal (both HS): {sig.sum()} of {len(sig)} events; train {trS.sum()}, score {teS.sum()}")
print(f"background (>= 1 PU): {bkg.sum()} of {len(bkg)} events; train {trB.sum()}, score {teB.sum()}"
      + ("  per JZ slice: " + ", ".join(f"JZ{z} {int(np.sum(bkg & (B['jz'] == z)))}" for z in np.unique(B["jz"])) if is_dijet else ""))

out = {"sig_weight": S["weight"][teS].astype(float), "bkg_weight": B["weight"][teB].astype(float),
       "bkg_jz": B["jz"][teB].astype(int)}
for name in FS:
    X = np.vstack([FS[name][trS], FB[name][trB]])
    y = np.concatenate([np.ones(trS.sum()), np.zeros(trB.sum())])
    w = np.concatenate([ws_tr, wb_tr])
    clf = xgb.XGBClassifier(n_estimators=300, max_depth=3, learning_rate=0.05, subsample=0.8,
                            tree_method="hist", random_state=args.seed, missing=np.nan, n_jobs=4)
    clf.fit(X, y, sample_weight=w)
    s_te, b_te = clf.predict_proba(FS[name][teS])[:, 1], clf.predict_proba(FB[name][teB])[:, 1]
    yt = np.concatenate([np.ones(len(s_te)), np.zeros(len(b_te))])
    wt = np.concatenate([np.abs(out["sig_weight"]) / np.abs(out["sig_weight"]).sum(),
                         wb[teB] / wb[teB].sum()])
    auc_te = roc_auc_score(yt, np.concatenate([s_te, b_te]), sample_weight=wt)
    auc_tr = roc_auc_score(y, clf.predict_proba(X)[:, 1], sample_weight=w)
    print(f"  BDT {name:9s} AUC train {auc_tr:.4f}  scored half {auc_te:.4f}")
    out[f"sig_{name}"] = s_te
    out[f"bkg_{name}"] = b_te
np.savez_compressed(args.out, **out)
print(f"wrote {args.out}")
