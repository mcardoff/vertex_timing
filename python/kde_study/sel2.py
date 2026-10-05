"""Top-2 tie-break scan. For events where exactly one of (winner, distinct runner-up) is right:
how often does 'take the candidate with the larger X' choose correctly, against the score's own accuracy?
Population-wide: all such events, and separately the low-Q half where a switch could matter."""
from sel import *


def cvars(d, c):
    t = d['tracks']; pc, pt_, dt, K, base, wt = d['_pairs']; nc = len(c['ev']); crow = c['crow']
    v = {}
    g = lambda x: x[pt_]
    v['S'] = c['S']; v['dens'] = c['dens']; v['env'] = c['env']
    v['neff'] = np.bincount(pc, K, nc)
    v['wz'] = np.bincount(pc, g(np.exp(-t['dz_d'])) * K, nc)
    v['sumptz'] = np.bincount(pc, g(t['pt_d'] * np.exp(-t['dz_d'])) * K, nc)
    v['meanw'] = v['dens'] / v['neff']
    v['zcompact'] = v['wz'] / v['neff']
    sw = np.bincount(pc, g(wt) * K, nc)
    v['-krms'] = -np.sqrt(np.bincount(pc, g(wt) * K * dt ** 2, nc) / sw)
    for Wn in (20., 30.):
        v[f'peak{int(Wn)}'] = np.bincount(pc, base[pt_] * np.exp(-0.5 * (dt / Wn) ** 2), nc) / c['dens']
    side = np.bincount(pc, base[pt_] * (np.exp(-0.5 * ((dt - 150) / 60.) ** 2) + np.exp(-0.5 * ((dt + 150) / 60.) ** 2)), nc)
    v['contrast'] = c['dens'] / (c['dens'] + 0.5 * side)
    v['dens-0.5side'] = c['dens'] - 0.5 * side
    sig = g(t['timeRes_d'])
    Ks = np.exp(-0.5 * (dt / (2.0 * sig)) ** 2)
    v['dens_sig2'] = np.bincount(pc, base[pt_] * Ks, nc)
    v['-tchi2k'] = -np.bincount(pc, base[pt_] * K * (dt / sig) ** 2, nc) / c['dens']
    v['n_injet'] = np.bincount(pc, g(t['injet'].astype(float)) * K, nc)
    v['injet_frac'] = np.bincount(pc, base[pt_] * g(t['injet'].astype(float)) * K, nc) / c['dens']
    v['prec'] = np.bincount(pc, K / g(t['sigma_z0_d']) ** 2, nc)
    v['n_close'] = np.bincount(pc, g((t['dz_d'] < 0.5).astype(float)) * K, nc)
    v['maxw'] = np.zeros(nc); np.maximum.at(v['maxw'], pc, base[pt_] * K)
    v['-topfrac'] = -v['maxw'] / c['dens']
    v['nhgtd'] = np.bincount(pc, base[pt_] * K * g(t['nhgtd_hits'].astype(float)), nc) / c['dens']
    v['-timeres'] = -np.bincount(pc, base[pt_] * K * sig, nc) / c['dens']
    v['abseta'] = np.bincount(pc, base[pt_] * K * np.abs(g(t['eta'].astype(float))), nc) / c['dens']
    zk = np.bincount(pc, base[pt_] * K * g(t['z0_d']), nc) / c['dens'] - d['vz'][c['ev']]
    v['-|zk|'] = -np.abs(zk)
    # z-coherence of the kernel tracks about their own kernel centroid
    izv = 1 / g(t['sigma_z0_d']) ** 2
    zc = np.bincount(pc, K * izv * g(t['z0_d']), nc) / np.bincount(pc, K * izv, nc)
    v['-zchi2k'] = -np.bincount(pc, K * izv * (g(t['z0_d']) - zc[pc]) ** 2, nc) / np.maximum(v['neff'], 1e-9)
    v['-zsig_k'] = -np.abs(zc - d['vz'][c['ev']]) * np.sqrt(np.bincount(pc, K * izv, nc))
    # hard cluster
    v['ncl'] = c['n'].astype(float); v['-absdz'] = -np.abs(c['dz']); v['info'] = c['info']
    v['-tchi2'] = -np.bincount(crow, (t['time_d'] - c['time'][crow]) ** 2 / t['timeRes_d'] ** 2, nc) / np.maximum(c['n'] - 1, 1)
    v['-|t|'] = -np.abs(c['T'])
    return v


if __name__ == "__main__":
    res = {}
    for n, m in SAMPLES:
        d = get(n, m); c = cands(d); p1, p2, Q = top2(d, c); v = cvars(d, c)
        has = p2 >= 0; a = p1[has]; b = p2[has]; q = Q[has]
        dec = c['ok'][a] ^ c['ok'][b]                      # exactly one right
        first = c['ok'][a]                                 # the score's pick is the right one
        lowq = q < 1.0
        print(f"\n{n}: events with a runner-up {has.mean()*100:.1f}%; decisive (exactly one right) {dec.mean()*100:.1f}% of those; "
              f"score accuracy on decisive: all {first[dec].mean()*100:.1f}, Q<1 {first[dec & lowq].mean()*100:.1f}")
        for k, x in v.items():
            d1 = x[a] > x[b]; tie = x[a] == x[b]
            acc = lambda s: (np.mean(np.where(tie[s], 0.5, d1[s] == first[s])) * 100)
            res.setdefault(k, []).append((acc(dec), acc(dec & lowq)))
    print("\naccuracy of 'take the larger X' on decisive pairs:  all / Q<1      (zjets | dijet | vbf | ttbar)")
    order = sorted(res, key=lambda k: -min(r[1] for r in res[k]))
    for k in order:
        print(f"{k:14s} " + " | ".join(f"{a:5.1f} / {b:5.1f}" for a, b in res[k]))
