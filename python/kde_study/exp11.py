"""Out-of-sample checks of the frozen configuration: disjoint file ranges, mjj500 selections, mu=0 controls."""
import numpy as np
from wb import evmax_pick
from exp1 import prep, tzp_score, tzp_time, core
import exp6

SETS = [('vbf_novbs', 300, 120), ('ttbar_novbs', 180, 60), ('zjets_mjj500p0', None, 0), ('dijet_mjj500p0', None, 0),
        ('vbf_mjj500p0', 500, 380), ('ttbar_mjj500p0', 400, 300),
        ('vbf_mu0_novbs', 20, 0), ('ttbar_mu0_novbs', 20, 0), ('zeejets_mu0_novbs', None, 0)]
print(f"{'set':34s} {'events':>8s} {'TZP':>7s} {'new':>7s} {'no-HGTD':>8s}  fails removed (new / no-HGTD)")
for name, mx, mn in SETS:
    d = prep(name, mx, mn); c = d['clusters']; t = d['tracks']
    S0 = tzp_score(d); T0 = tzp_time(c); d['T0'] = T0
    base = core(d, evmax_pick(c['ev'], S0, d['nev']), T0)
    ts = np.searchsorted(t['ev'], np.arange(d['nev'])); te = np.searchsorted(t['ev'], np.arange(d['nev']), side='right')
    cnt = (te - ts)[c['ev']]
    d['pc'] = np.repeat(np.arange(len(c['ev'])), cnt)
    d['pt'] = np.repeat(ts[c['ev']], cnt) + np.arange(cnt.sum()) - np.repeat(np.cumsum(cnt) - cnt, cnt)
    d['info'] = np.bincount(t['crow'], 1 / t['sigma_d0'] ** 2, len(c['ev']))
    from exp7 import kde
    r = []
    for hg in (0.5, 0.0):
        S, m = kde(d, p=0.5, hg=hg, pt_t=0.5)
        m = np.where(np.isfinite(m), m, T0)
        r.append(core(d, evmax_pick(c['ev'], S, d['nev']), m).mean())
    b = base.mean()
    print(f"{name+f'[{mn}:{mx}]':34s} {d['nev']:8d} {b*100:7.2f} {r[0]*100:7.2f} {r[1]*100:8.2f}  {(r[0]-b)/(1-b)*100:+6.1f}% / {(r[1]-b)/(1-b)*100:+6.1f}%"
          f"   (+-{np.sqrt(b*(1-b)/d['nev'])*100:.2f})")
    del d
