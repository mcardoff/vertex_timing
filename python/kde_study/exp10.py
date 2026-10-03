"""Round 8: last selection extras on top of (p=.5, hgtd .5, time pT^.5)."""
import numpy as np, uproot
from exp7 import run, kde, DS, SAMPLES
from wb import D

B = dict(p=0.5, hg=0.5, pt_t=0.5)
print("% of fails removed:                                                      zjets   dijet   vbf    ttbar")
run("reference", lambda d: kde(d, **B))
for pd in (0.15, 0.3):
    run(f"x n_eff^{pd}", lambda d, pd=pd: kde(d, pd=pd, **B))
for p in (0.35, 0.65):
    run(f"p {p}", lambda d, p=p: kde(d, p=p, hg=0.5, pt_t=0.5))
run("p .5 + cap 10", lambda d: kde(d, cap=10, **B))
for ptt in (0.0, 0.25):
    run(f"time pT^{ptt}", lambda d, ptt=ptt: kde(d, p=0.5, hg=0.5, pt_t=ptt))
for hg in (0.3, 0.8):
    run(f"hgtd {hg}", lambda d, hg=hg: kde(d, p=0.5, hg=hg, pt_t=0.5))
# leptons on zjets
a = uproot.open(D + "zjets_novbs_training.root")['tracks'].arrays(['file_idx', 'event_num', 'is_lepton'], library='np')
key = a['file_idx'].astype(np.int64) * 10_000_000 + a['event_num'].astype(np.int64)
isl = a['is_lepton'][np.argsort(key, kind='stable')]
d = DS['zjets_novbs']
print("zjets: lepton timed-track fraction", isl.mean(), " events with one", np.mean(np.bincount(d['tracks']['ev'], isl, d['nev']) > 0))
