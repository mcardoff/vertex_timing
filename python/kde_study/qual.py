"""Quality variables for the TZP_KDE_TZ pick: per-event reco-only candidates + truth labels."""
from reclus import *


def qvars(d):
    """returns dict of per-event arrays for the TZP_KDE_TZ pick"""
    t = d['tracks']; nev = d['nev']
    lab = cluster(d, cut=3.0, useZ=True, seedw=t['w_tzp'])
    c, crow = features(d, lab); nc = len(c['ev'])
    S, T = kde_eval(d, c, crow)
    p = evmax_pick(c['ev'], S, nev)
    m = T[p]                                        # reported time per event
    q = {'ok': np.abs(m - d['ttruth']) < 60, 'absdt': np.abs(m - d['ttruth'])}
    S1 = S[p]
    # runner-up among candidates whose time is a DIFFERENT mode (> 60 ps away)
    far = np.abs(T - m[c['ev']]) > 60
    S2 = np.zeros(nev); np.maximum.at(S2, c['ev'][far], S[far])
    q['S1'] = S1; q['margin'] = 1 - S2 / S1
    q['logS1'] = np.log(S1); q['logS2'] = np.log(S2 + 1e-3); q['S1mS2'] = S1 - S2
    # kernel sums at the reported time, over all event tracks
    k = np.exp(-0.5 * ((t['time_d'] - m[t['ev']]) / 60.) ** 2)
    base = t['pt_d'] ** 0.5 * np.exp(-t['dz_d'])
    wt = base / t['timeRes_d'] ** 2
    dens = np.bincount(t['ev'], base * k, nev); tot = np.bincount(t['ev'], base, nev)
    q['dens'] = dens; q['frac'] = dens / tot
    q['neff'] = np.bincount(t['ev'], k, nev)
    q['ntrk'] = np.bincount(t['ev'], minlength=nev).astype(float)
    q['nfrac'] = q['neff'] / q['ntrk']
    sw = np.bincount(t['ev'], wt * k, nev)
    q['sigma_m'] = 1 / np.sqrt(np.bincount(t['ev'], k / t['timeRes_d'] ** 2, nev))
    q['krms'] = np.sqrt(np.bincount(t['ev'], wt * k * (t['time_d'] - m[t['ev']]) ** 2, nev) / sw)
    q['wz'] = np.bincount(t['ev'], np.exp(-t['dz_d']) * k, nev)                 # z-weighted track count
    q['sumpt_k'] = np.bincount(t['ev'], t['pt_d'] * np.exp(-t['dz_d']) * k, nev)
    q['n_injet_k'] = np.bincount(t['ev'], t['injet'] * k, nev)
    q['n_close'] = np.bincount(t['ev'], (t['dz_d'] < 0.5) * k, nev)            # tracks within 0.5 mm of the PV
    q['prec_k'] = np.bincount(t['ev'], k / t['sigma_z0_d'] ** 2 * (t['dz_d'] < 3 * t['sigma_z0_d']), nev)
    # hard cluster of the pick
    q['ncl_trk'] = c['n'][p].astype(float); q['absdz'] = np.abs(c['dz'][p]); q['env'] = c['env'][p]; q['info'] = c['info'][p]
    q['shift'] = np.abs(m - c['time'][p])
    tch = np.bincount(crow, (t['time_d'] - c['time'][crow]) ** 2 / t['timeRes_d'] ** 2, nc) / np.maximum(c['n'] - 1, 1)
    q['tchi2'] = tch[p]
    q['nclus'] = np.bincount(c['ev'], minlength=nev).astype(float)
    # in-jet vs out-of-jet kernel times around the mode
    for nm, sel in (('in', t['injet']), ('out', ~t['injet'])):
        den = np.bincount(t['ev'], wt * k * sel, nev)
        with np.errstate(invalid='ignore', divide='ignore'):
            q['t_' + nm] = np.bincount(t['ev'], wt * k * sel * t['time_d'], nev) / den
    q['inout'] = np.abs(q['t_in'] - q['t_out'])
    # truth purities
    hs = t['truth_is_hs'] > 0.5
    picked = p[t['ev']] == crow
    q['purity_cluster'] = np.bincount(t['ev'], t['pt_d'] * hs * picked, nev) / np.bincount(t['ev'], t['pt_d'] * picked, nev)
    q['purity_kernel'] = np.bincount(t['ev'], wt * k * hs, nev) / sw
    q['m'] = m
    return q


def baseline_tzp(d):
    t = d['tracks']; c, crow = features(d, cluster(d))
    p = evmax_pick(c['ev'], c['S_tzp'], d['nev'])
    hs = t['truth_is_hs'] > 0.5; picked = p[t['ev']] == crow
    pur = np.bincount(t['ev'], t['pt_d'] * hs * picked, d['nev']) / np.bincount(t['ev'], t['pt_d'] * picked, d['nev'])
    return np.abs(c['t_tzp'][p] - d['ttruth']) < 60, pur


def athena(d):
    c = d['clusters']; nev = d['nev']
    valid = np.zeros(nev, bool); tm = np.zeros(nev)
    valid[c['ev']] = c['hgtd_valid'] > 0.5; tm[c['ev']] = c['hgtd_time']
    return valid, np.abs(tm - d['ttruth']) < 60


def auc(x, y):
    """P(x_good > x_bad)"""
    o = np.argsort(x, kind='stable'); r = np.empty(len(x)); r[o] = np.arange(len(x))
    n1 = y.sum(); n0 = len(y) - n1
    return (r[y].sum() - n1 * (n1 - 1) / 2) / (n1 * n0)


def purity_at(x, ok, keep):
    """core fraction of the `keep` fraction of events with the largest x"""
    thr = np.quantile(x, 1 - keep)
    s = x >= thr
    return ok[s].mean()


if __name__ == "__main__":
    VARS = ['S1', 'margin', 'S1mS2', 'dens', 'frac', 'neff', 'nfrac', 'ntrk', 'sigma_m', 'krms', 'wz', 'sumpt_k', 'n_injet_k', 'n_close',
            'prec_k', 'ncl_trk', 'absdz', 'env', 'info', 'shift', 'tchi2', 'nclus']
    out = {}
    for n, m in SAMPLES:
        d = get(n, m); q = qvars(d); ok = q['ok']
        v, aok = athena(d)
        print(f"\n{n}: core {ok.mean()*100:.2f}  cluster purity {np.nanmean(q['purity_cluster'])*100:.1f}  kernel purity {np.nanmean(q['purity_kernel'])*100:.1f}"
              f"   | Athena: time provided {v.mean()*100:.1f}%, core among provided {aok[v].mean()*100:.1f}%, ours on the same events {ok[v].mean()*100:.1f}%")
        print(f"{'variable':10s} {'AUC':>6s}  core among best 90% / 80% / 70% / 50% of events")
        for k in VARS:
            x = np.nan_to_num(q[k], nan=-1e9)
            a = auc(x, ok); sgn = 1 if a >= 0.5 else -1
            print(f"{k:10s} {max(a,1-a):6.3f}  " + " / ".join(f"{purity_at(sgn * x, ok, f)*100:5.1f}" for f in (.9, .8, .7, .5)) + ("" if sgn > 0 else "   (low is good)"))
        has = np.isfinite(q['inout'])
        print(f"inout defined {has.mean()*100:.1f}%: core | <30ps {ok[has & (q['inout']<30)].mean()*100:.1f} ({np.mean(has & (q['inout']<30))*100:.1f}% of ev), "
              f"30-60 {ok[has & (q['inout']>=30) & (q['inout']<60)].mean()*100:.1f}, >60 {ok[has & (q['inout']>=60)].mean()*100:.1f}, undefined {ok[~has].mean()*100:.1f}")
