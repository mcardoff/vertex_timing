"""With a PERFECT partition the score still picks wrong. What separates the pure HS cluster from the cluster picked instead?"""
from reclus import *
sw = lambda d: d['tracks']['w_tzp']
for n, m in (SAMPLES[0], SAMPLES[2]):
    d = get(n, m); t = d['tracks']; hs = t['truth_is_hs'] > 0.5
    intime = hs & (np.abs(t['time_d'] - d['ttruth'][t['ev']]) < 3 * t['timeRes_d'])
    lab = np.where(intime, 4000, cluster(d, cut=3.0, useZ=True, seedw=sw(d))).astype(np.int32)
    c, crow = features(d, lab)
    S, T = kde_eval(d, c, crow)
    nc = len(c['ev'])
    key = t['ev'].astype(np.int64) * 4096 + lab; uk = np.unique(key); is_hs = (uk % 4096) == 4000
    izv = 1 / t['sigma_z0_d'] ** 2
    zc = c['dz'] + d['vz'][c['ev']]
    c['sumw'] = np.bincount(crow, t['w_tzp'], nc)
    c['sumpt'] = np.bincount(crow, t['pt_d'], nc)
    c['sig_zc'] = 1 / np.sqrt(np.bincount(crow, izv, nc))
    c['dz_signif'] = np.abs(c['dz']) / c['sig_zc']
    c['zchi2_centroid'] = np.bincount(crow, (t['z0_d'] - zc[crow]) ** 2 * izv, nc) / np.maximum(c['n'] - 1, 1)
    c['zchi2_pv'] = np.bincount(crow, t['dz_d'] ** 2 * izv, nc) / c['n']
    c['mean_dz'] = np.bincount(crow, t['dz_d'], nc) / c['n']
    c['min_sigz'] = np.full(nc, 99.); np.minimum.at(c['min_sigz'], crow, t['sigma_z0_d'])
    c['tchi2'] = np.bincount(crow, (t['time_d'] - c['time'][crow]) ** 2 / t['timeRes_d'] ** 2, nc) / np.maximum(c['n'] - 1, 1)
    c['absdz'] = np.abs(c['dz']); c['S'] = S
    p = evmax_pick(c['ev'], S, d['nev'])
    hsrow = np.full(d['nev'], -1); hsrow[c['ev'][is_hs]] = np.where(is_hs)[0]
    good = (hsrow >= 0)
    good[good] &= np.abs(c['time'][hsrow[good]] - d['ttruth'][good]) < 60
    lose = good & (p != hsrow) & (np.abs(T[p] - d['ttruth']) >= 60)
    win = good & (p == hsrow)
    print(f"\n{n}: HS cluster exists & in window {good.mean()*100:.1f}% of events; picked {win.mean()*100:.1f}%; lost to another cluster {lose.mean()*100:.1f}%")
    print(f"{'median':16s} {'HS (lost)':>10s} {'picked PU':>10s} | {'HS (won)':>10s}   AUC(HS-lost vs picked)")
    for k in ('n', 'sumpt', 'sumw', 'absdz', 'dz_signif', 'sig_zc', 'zchi2_centroid', 'zchi2_pv', 'mean_dz', 'min_sigz', 'info', 'tchi2', 'S'):
        a = c[k][hsrow[lose]]; b = c[k][p[lose]]
        auc = (np.mean(a > b) + 0.5 * np.mean(a == b))
        print(f"{k:16s} {np.median(a):10.3f} {np.median(b):10.3f} | {np.median(c[k][hsrow[win]]):10.3f}   P(HS > picked) = {auc:.3f}")
