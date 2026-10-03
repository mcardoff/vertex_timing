"""Frozen configuration with a single-cluster guard: keep the plain cluster time when the event has one cluster."""
import numpy as np
from wb import evmax_pick
from exp1 import prep, tzp_score, tzp_time, core

SETS = [('zjets_novbs', None, 0), ('dijet_novbs', None, 0), ('vbf_novbs', 120, 0), ('ttbar_novbs', 60, 0),
        ('vbf_mu0_novbs', 20, 0), ('ttbar_mu0_novbs', 20, 0), ('zeejets_mu0_novbs', None, 0)]
for name, mx, mn in SETS:
    d = prep(name, mx, mn); c = d['clusters']; t = d['tracks']; nev = d['nev']
    S0 = tzp_score(d); T0 = tzp_time(c); d['T0'] = T0
    base = core(d, evmax_pick(c['ev'], S0, nev), T0)
    ts = np.searchsorted(t['ev'], np.arange(nev)); te = np.searchsorted(t['ev'], np.arange(nev), side='right')
    cnt = (te - ts)[c['ev']]
    d['pc'] = np.repeat(np.arange(len(c['ev'])), cnt)
    d['pt'] = np.repeat(ts[c['ev']], cnt) + np.arange(cnt.sum()) - np.repeat(np.cumsum(cnt) - cnt, cnt)
    d['info'] = np.bincount(t['crow'], 1 / t['sigma_d0'] ** 2, len(c['ev']))
    from exp7 import kde
    S, m = kde(d, p=0.5, hg=0.5, pt_t=0.5)
    m = np.where(np.isfinite(m), m, T0)
    p = evmax_pick(c['ev'], S, nev)
    ncl = np.bincount(c['ev'], minlength=nev)
    a = core(d, p, m)
    for k in (1, 2):
        mg = np.where(ncl[c['ev']] <= k, T0, m)
        g = core(d, p, mg)
        print(f"{name:20s} nev {nev:7d}  TZP {base.sum():7d} ({base.mean()*100:.2f})  new {a.sum():7d} ({a.mean()*100:.2f})  "
              f"guard ncl<={k}: {g.sum():7d} ({g.mean()*100:.2f})  [{np.mean(ncl<=k)*100:.1f}% of ev]")
