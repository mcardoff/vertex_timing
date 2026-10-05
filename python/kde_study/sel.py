"""Selection study on the (t,z) + kernel-density baseline: per-CANDIDATE kernel variables, top-2 structure,
runner-up / aggregate fallbacks as a function of the quality flag."""
from reclus import *


def cands(d, W=60., sfac=None, tfac=None, seedfac=None):
    """per-cluster (candidate) arrays on the TZP_KDE_TZ configuration, plus pair index arrays"""
    t = d['tracks']; nev = d['nev']
    lab = cluster(d, cut=3.0, useZ=True, seedw=t['w_tzp'] * (1.0 if seedfac is None else seedfac))
    c, crow = features(d, lab); nc = len(c['ev'])
    cnt = (d['off'][1:] - d['off'][:-1])[c['ev']]
    pc = np.repeat(np.arange(nc), cnt)
    pt_ = np.repeat(d['off'][:-1][c['ev']], cnt) + np.arange(cnt.sum()) - np.repeat(np.cumsum(cnt) - cnt, cnt)
    base = t['pt_d'] ** 0.5 * np.exp(-t['dz_d']); wt = base / t['timeRes_d'] ** 2
    if tfac is not None: wt = wt * tfac
    sbase = base if sfac is None else base * sfac
    m = np.bincount(crow, wt * t['time_d'], nc) / np.bincount(crow, wt, nc)
    tt = t['time_d'][pt_]; w = wt[pt_]
    for _ in range(3):
        k = w * np.exp(-0.5 * ((tt - m[pc]) / W) ** 2)
        den = np.bincount(pc, k, nc); num = np.bincount(pc, k * tt, nc)
        m = np.where(den > 0, num / np.where(den > 0, den, 1), m)
    ncl = np.bincount(c['ev'], minlength=nev)
    c['T'] = np.where(ncl[c['ev']] <= 1, c['t_tzp'], m)
    c['m'] = m
    dt = tt - m[pc]
    K = np.exp(-0.5 * (dt / W) ** 2)
    c['dens'] = np.bincount(pc, sbase[pt_] * K, nc)
    c['S'] = c['dens'] * c['env']
    c['ok'] = np.abs(c['T'] - d['ttruth'][c['ev']]) < 60
    d['_pairs'] = (pc, pt_, dt, K, base, wt)
    c['crow'] = crow
    return c


def top2(d, c):
    """winner and best DISTINCT-mode runner-up per event (row indices; -1 if none)"""
    nev = d['nev']
    p1 = evmax_pick(c['ev'], c['S'], nev)
    far = np.abs(c['T'] - c['T'][p1][c['ev']]) > 60
    s2 = np.where(far, c['S'], -1.0)
    p2 = evmax_pick(c['ev'], s2, nev)
    has2 = s2[p2] > 0
    p2 = np.where(has2, p2, -1)
    S1 = c['S'][p1]; S2 = np.where(has2, c['S'][np.maximum(p2, 0)], 0.0)
    Q = (S1 - S2) / np.sqrt(S1 + S2)
    return p1, p2, Q


def aggregate_t0(d, p=1.0):
    """classical aggregate: weighted median of all event tracks, then two Winsorised weighted means"""
    t = d['tracks']; nev = d['nev']
    w = t['pt_d'] ** p * np.exp(-t['dz_d']) / t['timeRes_d'] ** 2
    o = np.lexsort((t['time_d'], t['ev'])); ev = t['ev'][o]; tm = t['time_d'][o]; ww = w[o]
    cw = np.cumsum(ww); tot = np.bincount(ev, ww, nev)
    start = np.r_[0, np.cumsum(tot)[:-1]]
    frac = (cw - start[ev]) / tot[ev]
    first = np.r_[True, (frac[1:] >= 0.5) & ((frac[:-1] < 0.5) | (ev[1:] != ev[:-1]))] & (frac >= 0.5)
    med = np.zeros(nev); idx = np.where(first)[0]
    med[ev[idx][::-1]] = tm[idx][::-1]          # first crossing per event wins
    m = med
    for win in (90., 60.):
        cl = np.clip(t['time_d'], m[t['ev']] - win, m[t['ev']] + win)
        m = np.bincount(t['ev'], w * cl, nev) / np.bincount(t['ev'], w, nev)
    return m


if __name__ == "__main__":
    edges = [-1, 0.5, 1.0, 1.5, 2.0, 3.0, 4.0, 1e9]
    for n, m in SAMPLES:
        d = get(n, m); c = cands(d); p1, p2, Q = top2(d, c)
        ok1 = c['ok'][p1]; ok2 = np.where(p2 >= 0, c['ok'][np.maximum(p2, 0)], False)
        agg = aggregate_t0(d); oka = np.abs(agg - d['ttruth']) < 60
        anyc = np.bincount(c['ev'], c['ok'], d['nev']) > 0
        print(f"\n{n}: pick right {ok1.mean()*100:.2f} | runner-up right {ok2.mean()*100:.2f} | aggregate right {oka.mean()*100:.2f} | "
              f"pick or runner-up {np.mean(ok1|ok2)*100:.2f} | pick or aggregate {np.mean(ok1|oka)*100:.2f} | any candidate {anyc.mean()*100:.2f}")
        print(f"  {'Q bin':11s} {'share':>6s} {'pick':>6s} {'2nd':>6s} {'agg':>6s} {'pick|2nd':>9s} {'pick|agg':>9s} {'any cand':>9s}   of pick-wrong: 2nd right / agg right / another cand right")
        for i in range(len(edges) - 1):
            s = (Q >= edges[i]) & (Q < edges[i + 1]); w_ = s & ~ok1
            print(f"  {edges[i]:4.1f}-{min(edges[i+1],9.9):4.1f} {s.mean()*100:6.1f} {ok1[s].mean()*100:6.1f} {ok2[s].mean()*100:6.1f} {oka[s].mean()*100:6.1f} "
                  f"{np.mean(ok1[s]|ok2[s])*100:9.1f} {np.mean(ok1[s]|oka[s])*100:9.1f} {anyc[s].mean()*100:9.1f}   "
                  f"{ok2[w_].mean()*100:5.1f} / {oka[w_].mean()*100:5.1f} / {anyc[w_].mean()*100:5.1f}")
