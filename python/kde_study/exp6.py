"""Round 5: mean-shift re-timing details + kernel-density (soft-cluster) SELECTION."""
import sys, numpy as np
from wb import evmax_pick
from exp1 import prep, tzp_score, tzp_time, core

SAMPLES = [('zjets_novbs', None), ('dijet_novbs', None), ('vbf_novbs', 120), ('ttbar_novbs', 60)]
DS = {n: prep(n, m) for n, m in SAMPLES}
for n, _ in SAMPLES:
    d = DS[n]; c = d['clusters']; t = d['tracks']
    d['S0'] = tzp_score(d); d['T0'] = tzp_time(c)
    d['pick'] = evmax_pick(c['ev'], d['S0'], d['nev'])
    d['base_ok'] = core(d, d['pick'], d['T0'])
    assert (np.diff(t['ev']) >= 0).all()
    # (cluster,track) pairs within an event
    ts = np.searchsorted(t['ev'], np.arange(d['nev'])); te = np.searchsorted(t['ev'], np.arange(d['nev']), side='right')
    cnt = (te - ts)[c['ev']]
    pc = np.repeat(np.arange(len(c['ev'])), cnt)
    off = np.arange(cnt.sum()) - np.repeat(np.cumsum(cnt) - cnt, cnt)
    pt_ = np.repeat(ts[c['ev']], cnt) + off
    d['pc'] = pc; d['pt'] = pt_
    d['info'] = np.bincount(t['crow'], 1 / t['sigma_d0'] ** 2, len(c['ev']))


def ms_all(d, wt, seed, W, nit=3, sigscale=0.0, kern='gauss'):
    """mean-shift from every cluster's seed time over all event tracks. returns (mode time, density) per cluster"""
    t = d['tracks']; pc = d['pc']; pt_ = d['pt']; nc = len(seed)
    m = seed.copy(); tt = t['time'][pt_]; w = wt[pt_]
    wid2 = W ** 2 + (sigscale * t['timeRes'][pt_]) ** 2
    for i in range(nit + 1):
        d2 = (tt - m[pc]) ** 2 / wid2
        k = np.exp(-0.5 * d2) if kern == 'gauss' else (d2 < 1).astype(float)
        den = np.bincount(pc, w * k, nc)
        if i == nit: break
        num = np.bincount(pc, w * k * tt, nc)
        with np.errstate(invalid='ignore', divide='ignore'):
            m2 = num / den
        m = np.where(np.isfinite(m2), m2, m)
    return m, den


def run(label, fn):
    """fn(d) -> (score per cluster, time per cluster)"""
    out = []
    for n, _ in SAMPLES:
        d = DS[n]; c = d['clusters']
        S, T = fn(d)
        T = np.where(np.isfinite(T), T, d['T0'])
        ok = core(d, evmax_pick(c['ev'], S, d['nev']), T)
        base = d['base_ok']
        out.append((ok.mean() - base.mean()) * 100 / (1 - base.mean()))
    print(f"{label:66s} " + "  ".join(f"{x:+6.2f}" for x in out) + f"   worst {min(out):+.2f}  mean {np.mean(out):+.2f}")
    return out


def wtime(t):
    return t['pt'] * np.exp(-t['dz']) / t['timeRes'] ** 2


def seed_time(d):
    t = d['tracks']; nc = len(d['clusters']['ev']); w = wtime(t)
    return np.bincount(t['crow'], w * t['time'], nc) / np.bincount(t['crow'], w, nc)


if __name__ == "__main__":
    print("% of fails removed:                                                      zjets   dijet   vbf    ttbar")
    print("-- timing only (TZP selection) --")
    for W in (45., 55., 60., 65., 75.):
        run(f"mean-shift gauss W {W} nit 3", lambda d, W=W: (d['S0'], ms_all(d, wtime(d['tracks']), seed_time(d), W)[0]))
    for nit in (1, 2, 6):
        run(f"mean-shift gauss W 60 nit {nit}", lambda d, nit=nit: (d['S0'], ms_all(d, wtime(d['tracks']), seed_time(d), 60., nit)[0]))
    run("mean-shift W 60, seed = plain cluster time", lambda d: (d['S0'], ms_all(d, wtime(d['tracks']), d['clusters']['cluster_time'].astype(float), 60.)[0]))
    for W, ss in ((50., 1.0), (40., 1.5), (0., 2.0), (0., 2.5), (0., 3.0)):
        run(f"mean-shift width^2 = {W}^2 + ({ss} sig_t)^2", lambda d, W=W, ss=ss: (d['S0'], ms_all(d, wtime(d['tracks']), seed_time(d), W, 3, ss)[0]))
    for W in (60., 90., 120.):
        run(f"mean-shift BOX half-width {W}", lambda d, W=W: (d['S0'], ms_all(d, wtime(d['tracks']), seed_time(d), W, 3, 0, 'box')[0]))
    run("mean-shift W 60, weight w/o 1/sig_t^2", lambda d: (d['S0'], ms_all(d, d['tracks']['pt'] * np.exp(-d['tracks']['dz']), seed_time(d), 60.)[0]))
    print("-- KDE selection: score = density at the mode (x cluster envelope), time = mode --")
    for W in (40., 60., 80.):
        def f(d, W=W, env=True):
            t = d['tracks']; c = d['clusters']
            m, _ = ms_all(d, wtime(t), seed_time(d), 60.)
            _, den = ms_all(d, t['pt'] * np.exp(-t['dz']), m, W, 0)
            e = np.exp(-0.6 * np.abs(c['delta_z'])) * d['info'] ** 0.225 if env else 1.0
            return den * e, m
        run(f"KDE-select W {W} + env", f)
        run(f"KDE-select W {W} no env", lambda d, W=W: f(d, W, False))
