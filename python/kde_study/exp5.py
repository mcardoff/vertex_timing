"""Round 4: richer per-track timing weights + windowed (mean-shift) re-timing. Selection fixed = TZP."""
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
    t['injet'] = (t['dr_nearest_fwdjet'] < 0.4) & np.isfinite(t['dr_nearest_fwdjet'])
    t['picked'] = d['pick'][t['ev']] == t['crow']


def run(label, tfn_ev):
    """tfn_ev(d) -> per-EVENT time (nan = fall back to baseline)"""
    out = []
    for n, _ in SAMPLES:
        d = DS[n]
        Te = tfn_ev(d)
        base_t = d['T0'][d['pick']]
        Te = np.where(np.isfinite(Te), Te, base_t)
        ok = np.abs(Te - d['ttruth']) < 60
        base = d['base_ok']
        out.append((ok.mean() - base.mean()) * 100 / (1 - base.mean()))
    print(f"{label:62s} " + "  ".join(f"{x:+6.2f}" for x in out) + f"   worst {min(out):+.2f}  mean {np.mean(out):+.2f}")
    return out


def clus_time(d, w):
    t = d['tracks']; s = t['picked']
    num = np.bincount(t['ev'][s], (w * t['time'])[s], d['nev']); den = np.bincount(t['ev'][s], w[s], d['nev'])
    with np.errstate(invalid='ignore', divide='ignore'):
        return num / den


def base_w(t, p=1.0, a=1.0):
    return t['pt'] ** p * np.exp(-a * t['dz']) / t['timeRes'] ** 2


if __name__ == "__main__":
    print("% of fails removed:                                                  zjets   dijet   vbf    ttbar")
    run("TZP-weighted (p=1,a=1)", lambda d: clus_time(d, base_w(d['tracks'])))
    for kj in (2, 4, 8, 16):
        run(f"x (1 + {kj}*injet)", lambda d, kj=kj: clus_time(d, base_w(d['tracks']) * (1 + kj * d['tracks']['injet'])))
    for kj in (1, 3):
        def f(d, kj=kj):
            t = d['tracks']
            jw = np.where(t['injet'], kj * t['pt_nearest_fwdjet'] / 30.0 / np.maximum(t['dr_nearest_fwdjet'], 0.05) * 0.1, 0)
            return clus_time(d, base_w(t) * (1 + jw))
        run(f"x (1 + {kj}*WAVeS-like pTjet/dR)", f)
    for b in (0.5, 1.0, 2.0):
        run(f"x sigma_z0^-{b}", lambda d, b=b: clus_time(d, base_w(d['tracks']) * d['tracks']['sigma_z0'] ** -b))
        run(f"x sigma_d0^-{b}", lambda d, b=b: clus_time(d, base_w(d['tracks']) * d['tracks']['sigma_d0'] ** -b))
    for cap in (5, 10):
        run(f"pT cap {cap}", lambda d, cap=cap: clus_time(d, np.minimum(d['tracks']['pt'], cap) * np.exp(-d['tracks']['dz']) / d['tracks']['timeRes'] ** 2))
    run("no 1/sig_t^2", lambda d: clus_time(d, d['tracks']['pt'] * np.exp(-d['tracks']['dz'])))
    for k in (1.0,):
        for cc in (0.3, 1.0):
            def f(d, k=k, cc=cc):
                t = d['tracks']; s = k * t['sigma_z0']
                L = np.exp(-0.5 * (t['dz'] / s) ** 2) / (np.sqrt(2 * np.pi) * s)
                return clus_time(d, t['pt'] * L / (L + cc) / t['timeRes'] ** 2)
            run(f"LR prob k {k} c {cc}", f)
    # windowed mean-shift over ALL event tracks, seeded at the weighted cluster time
    for W in (40., 60., 90.):
        for kern in ('gauss',):
            def f(d, W=W):
                t = d['tracks']
                w = base_w(t)
                m = clus_time(d, w)
                for _ in range(3):
                    k = np.exp(-0.5 * ((t['time'] - m[t['ev']]) / W) ** 2)
                    k = np.where(np.isfinite(k), k, 0)
                    num = np.bincount(t['ev'], w * k * t['time'], d['nev']); den = np.bincount(t['ev'], w * k, d['nev'])
                    with np.errstate(invalid='ignore', divide='ignore'):
                        m2 = num / den
                    m = np.where(np.isfinite(m2), m2, m)
                return m
            run(f"mean-shift all tracks, gauss W {W}", f)
