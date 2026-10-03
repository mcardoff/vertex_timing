"""Round 1: cluster-level score variants against TZP, all four samples."""
import sys, numpy as np
from wb import evmax_pick
from exp1 import prep, tzp_score, tzp_time, core

SAMPLES = [('zjets_novbs', None), ('dijet_novbs', None), ('vbf_novbs', 120), ('ttbar_novbs', 60)]
DS = {n: prep(n, m) for n, m in SAMPLES}


def run(label, fn):
    """fn(d) -> (score per cluster row, time per cluster row)"""
    out = []
    for n, _ in SAMPLES:
        d = DS[n]; c = d['clusters']
        S, T = fn(d)
        p = evmax_pick(c['ev'], S, d['nev'])
        ok = core(d, p, T)
        base = d['base_ok']
        out.append((ok.mean() - base.mean()) * 100 / (1 - base.mean()) / 1.0)
    print(f"{label:46s} " + "  ".join(f"{x:+6.2f}" for x in out) + f"   worst {min(out):+.2f}")
    return out


for n, _ in SAMPLES:
    d = DS[n]; c = d['clusters']; t = d['tracks']
    d['S0'] = tzp_score(d); d['T0'] = tzp_time(c)
    d['base_ok'] = core(d, evmax_pick(c['ev'], d['S0'], d['nev']), d['T0'])
    nc = len(c['ev'])
    d['info'] = np.bincount(t['crow'], 1 / t['sigma_d0'] ** 2, nc)
    d['env'] = np.exp(-0.6 * np.abs(c['delta_z'])) * d['info'] ** 0.225
    print(n, d['nev'], f"base {d['base_ok'].mean()*100:.2f}", "track time mean/std", t['time'].mean(), t['time'].std(),
          "truth t std", d['ttruth'].std())

print("\n% of failing events removed (negative = more fails):       zjets   dijet   vbf    ttbar")
# V1 pT power
for p in (0.5, 0.75, 1.25, 1.5):
    def f(d, p=p):
        t = d['tracks']; nc = len(d['clusters']['ev'])
        return np.bincount(t['crow'], t['pt'] ** p * np.exp(-1.0 * t['dz']), nc) * d['env'], d['T0']
    run(f"pT^{p}", f)
# V1b pt cap
for cap in (5, 10, 15):
    def f(d, cap=cap):
        t = d['tracks']; nc = len(d['clusters']['ev'])
        return np.bincount(t['crow'], np.minimum(t['pt'], cap) * np.exp(-1.0 * t['dz']), nc) * d['env'], d['T0']
    run(f"pT cap {cap}", f)
# V2 background subtraction: expected PU weight under the cluster's time
for lam in (0.5, 1.0, 2.0, 4.0):
    for sig in (150., 200.):
        def f(d, lam=lam, sig=sig):
            t = d['tracks']; c = d['clusters']; nc = len(c['ev'])
            w = t['pt'] * np.exp(-1.0 * t['dz'])
            Wev = np.bincount(t['ev'], w, d['nev'])
            # window ~ 6 sigma_cluster-ish: use fixed 90 ps effective width
            g = np.exp(-0.5 * (c['cluster_time'] / sig) ** 2) / (sig * np.sqrt(2 * np.pi)) * 90.0
            s = np.bincount(t['crow'], w, nc) - lam * Wev[c['ev']] * g
            return np.maximum(s, 1e-9 * (1 + np.bincount(t['crow'], w, nc))) * d['env'], d['T0']
        run(f"bkg-sub lam {lam} sig {sig}", f)
# V2b multiplicative tail boost: S * exp(+k t^2/2sig^2)
for k in (0.1, 0.25, 0.5, 1.0):
    run(f"S / G(t)^{k} (sig 180)", lambda d, k=k: (d['S0'] * np.exp(k * 0.5 * (d['clusters']['cluster_time'] / 180.) ** 2), d['T0']))
# V3 HGTD agreement
for kap in (0.5, 1.0, 2.0, 5.0):
    for tau in (30., 60.):
        def f(d, kap=kap, tau=tau):
            c = d['clusters']
            dt = np.abs(c['cluster_time'] - c['hgtd_time'])
            a = np.where(c['hgtd_valid'] > 0.5, np.exp(-0.5 * (dt / tau) ** 2), 0.0)
            return d['S0'] * (1 + kap * a), d['T0']
        run(f"HGTD agree kappa {kap} tau {tau}", f)
# V4 compactness
for q in (0.1, 0.25, 0.5):
    run(f"S / (1+chi2ndf)^{q}", lambda d, q=q: (d['S0'] / (1 + np.nan_to_num(d['clusters']['time_chi2_ndf'])) ** q, d['T0']))
# V5 multiplicity
for q in (-0.25, 0.25, 0.5):
    run(f"S * n^{q}", lambda d, q=q: (d['S0'] * d['clusters']['n_tracks'] ** q, d['T0']))
# V6 time sigma
for q in (-0.5, 0.25, 0.5):
    run(f"S * sigma_t^{q}", lambda d, q=q: (d['S0'] * d['clusters']['cluster_time_sigma'] ** q, d['T0']))
