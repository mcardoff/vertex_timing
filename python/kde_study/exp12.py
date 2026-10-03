"""mu=0 harm diagnosis: which half (selection or time) loses events, and what the lost events look like."""
import numpy as np
from wb import evmax_pick
from exp1 import prep, tzp_score, tzp_time, core
from exp6 import ms_all

for name, mx in (('zeejets_mu0_novbs', None), ('vbf_mu0_novbs', 20)):
    d = prep(name, mx); c = d['clusters']; t = d['tracks']; nev = d['nev']
    S0 = tzp_score(d); T0 = tzp_time(c); d['T0'] = T0
    ts = np.searchsorted(t['ev'], np.arange(nev)); te = np.searchsorted(t['ev'], np.arange(nev), side='right')
    cnt = (te - ts)[c['ev']]
    d['pc'] = np.repeat(np.arange(len(c['ev'])), cnt)
    d['pt'] = np.repeat(ts[c['ev']], cnt) + np.arange(cnt.sum()) - np.repeat(np.cumsum(cnt) - cnt, cnt)
    d['info'] = np.bincount(t['crow'], 1 / t['sigma_d0'] ** 2, len(c['ev']))
    from exp7 import kde
    p0 = evmax_pick(c['ev'], S0, nev)
    base = core(d, p0, T0)
    S, m = kde(d, p=0.5, hg=0.0, pt_t=0.5)
    p1 = evmax_pick(c['ev'], S, nev)
    print(name, nev, "clusters/ev", len(c['ev']) / nev)
    print("  TZP", base.sum(), " new sel+new time", core(d, p1, m).sum(), " TZP sel+new time", core(d, p0, m).sum(), " new sel + T0", core(d, p1, T0).sum(),
          " TZP sel + plain cluster time", core(d, p0, c['cluster_time']).sum())
    lost = base & ~core(d, p0, m)
    print("  lost by re-timing alone:", lost.sum(), " gained:", (~base & core(d, p0, m)).sum())
    il = np.where(lost)[0][:6]
    for e in il:
        s = t['ev'] == e
        print(f"   ev {e}: truth {d['ttruth'][e]:.1f}  T0 {T0[p0[e]]:.1f}  ms {m[p0[e]]:.1f}  ncl {np.sum(c['ev']==e)}")
        o = np.argsort(t['time'][s])
        print("      t:", np.round(t['time'][s][o], 0), "\n      pt:", np.round(t['pt'][s][o], 1), "\n      cl:", t['cluster_idx'][s][o].astype(int), " picked", int(c['cluster_idx'][p0[e]]),
              "\n      hs:", t['truth_is_hs'][s][o].astype(int), "\n      sig:", np.round(t['timeRes'][s][o], 0))
    # alternatives for the time
    wt = t['pt'] ** 0.5 * np.exp(-t['dz']) / t['timeRes'] ** 2
    seed = np.bincount(t['crow'], wt * t['time'], len(c['ev'])) / np.bincount(t['crow'], wt, len(c['ev']))
    print("  weighted in-cluster time (no mean-shift):", core(d, p0, seed).sum())
    for W in (30., 45., 60.):
        for ss in (0., 1.0, 2.0):
            mm = ms_all(d, wt, seed, W, 3, ss)[0]
            print(f"  ms W {W} sigscale {ss}: {core(d, p0, np.where(np.isfinite(mm), mm, T0)).sum()}")
