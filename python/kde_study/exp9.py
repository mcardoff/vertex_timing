"""Anatomy of the new method vs TZP by n timed HS tracks; other windows; candidate oracle."""
import numpy as np
from wb import evmax_pick
from exp1 import core
from exp7 import kde, DS, SAMPLES

for n, _ in SAMPLES:
    d = DS[n]; c = d['clusters']; nev = d['nev']
    S, m = kde(d, p=0.5, hg=0.5, pt_t=0.5)
    p = evmax_pick(c['ev'], S, nev)
    dt_new = m[p] - d['ttruth']; dt_old = d['T0'][d['pick']] - d['ttruth']
    print(f"\n{n}: core60 old {np.mean(np.abs(dt_old)<60)*100:.2f} new {np.mean(np.abs(dt_new)<60)*100:.2f}")
    for w in (20, 30, 45, 60, 90, 150):
        fo = np.mean(np.abs(dt_old) >= w); fn = np.mean(np.abs(dt_new) >= w)
        print(f"   window {w:3d}: fail old {fo*100:5.2f} new {fn*100:5.2f}  removed {100*(fo-fn)/fo:+.1f}%")
    for lab, x in (('old', dt_old), ('new', dt_new)):
        core_ = x[np.abs(x) < 60]
        print(f"   {lab}: RMS inside 60ps {core_.std():.2f} ps, median |dt| {np.median(np.abs(x)):.2f}")
    nhs = np.bincount(c['ev'], c['truth_n_hs_tracks'], nev)
    cand_ok = np.bincount(c['ev'], np.abs(m - d['ttruth'][c['ev']]) < 60, nev) > 0
    print(f"   any-mode-passes {cand_ok.mean()*100:.2f}; picks changed {np.mean(p!=d['pick'])*100:.1f}%; "
          f"gained {np.mean((np.abs(dt_new)<60)&(np.abs(dt_old)>=60))*100:.2f} lost {np.mean((np.abs(dt_new)>=60)&(np.abs(dt_old)<60))*100:.2f}")
    for lo, hi in ((0, 0), (1, 1), (2, 3), (4, 7), (8, 999)):
        s = (nhs >= lo) & (nhs <= hi)
        print(f"   nHS {lo}-{hi}: {s.mean()*100:5.1f}% ev  old {np.mean(np.abs(dt_old[s])<60)*100:5.1f}  new {np.mean(np.abs(dt_new[s])<60)*100:5.1f}  any-mode {cand_ok[s].mean()*100:5.1f}")
