"""Round 6: stack terms on the KDE selection + mean-shift time."""
import sys, numpy as np
from wb import evmax_pick
from exp1 import core
from exp6 import DS, SAMPLES, ms_all, wtime, seed_time


def run(label, fn, show=True):
    out = []
    for n, _ in SAMPLES:
        d = DS[n]; c = d['clusters']
        S, T = fn(d)
        T = np.where(np.isfinite(T), T, d['T0'])
        ok = core(d, evmax_pick(c['ev'], S, d['nev']), T)
        base = d['base_ok']
        out.append((ok.mean() - base.mean()) * 100 / (1 - base.mean()))
    if show:
        print(f"{label:66s} " + "  ".join(f"{x:+6.2f}" for x in out) + f"   worst {min(out):+.2f}  mean {np.mean(out):+.2f}")
    return out


def kde(d, W=60., Wt=60., cap=None, p=1.0, a=1.0, at=1.0, pt_t=1.0, inj=0.0, b=0.6, g=0.45, soft_env=False, hg=0.0, tau=30., nit=3, pd=0.0):
    t = d['tracks']; c = d['clusters']; nc = len(c['ev'])
    pt = t['pt'] if cap is None else np.minimum(t['pt'], cap)
    injet = (t['dr_nearest_fwdjet'] < 0.4)
    wt = t['pt'] ** pt_t * np.exp(-at * t['dz']) / t['timeRes'] ** 2 * (1 + inj * injet)
    seed = np.bincount(t['crow'], wt * t['time'], nc) / np.bincount(t['crow'], wt, nc)
    m, _ = ms_all(d, wt, seed, Wt, nit)
    ws = pt ** p * np.exp(-a * t['dz'])
    _, den = ms_all(d, ws, m, W, 0)
    if soft_env:
        pc = d['pc']; pt_ = d['pt']
        k = np.exp(-0.5 * ((t['time'][pt_] - m[pc]) / W) ** 2)
        zc = np.bincount(pc, k * t['pt'][pt_] * t['z0'][pt_], nc) / np.bincount(pc, k * t['pt'][pt_], nc)
        info = np.bincount(pc, k / t['sigma_d0'][pt_] ** 2, nc)
        e = np.exp(-b * np.abs(zc - d['vz'][c['ev']])) * info ** (0.5 * g)
    else:
        e = np.exp(-b * np.abs(c['delta_z'])) * d['info'] ** (0.5 * g)
    S = den * e
    if pd:
        _, n_ = ms_all(d, np.ones(len(t['pt'])), m, W, 0)
        S = S * n_ ** pd
    if hg:
        dt = np.abs(m - c['hgtd_time'])
        S = S * (1 + hg * np.where(c['hgtd_valid'] > 0.5, np.exp(-0.5 * (dt / tau) ** 2), 0.0))
    return S, m


if __name__ == "__main__":
    print("% of fails removed:                                                      zjets   dijet   vbf    ttbar")
    run("KDE W60 + env (reference)", lambda d: kde(d))
    run("soft env", lambda d: kde(d, soft_env=True))
    for cap in (4, 6, 10):
        run(f"sel pT cap {cap}", lambda d, cap=cap: kde(d, cap=cap))
    for p in (0.5, 0.75):
        run(f"sel pT^{p}", lambda d, p=p: kde(d, p=p))
    for a in (0.6, 1.5):
        run(f"sel alpha {a}", lambda d, a=a: kde(d, a=a))
    for b in (0.3, 1.0):
        run(f"env beta {b}", lambda d, b=b: kde(d, b=b))
    for g in (0.2, 0.8):
        run(f"env gamma {g}", lambda d, g=g: kde(d, g=g))
    for hg in (0.5, 1.0, 2.0):
        run(f"HGTD agree {hg}", lambda d, hg=hg: kde(d, hg=hg))
    for inj in (1, 3):
        run(f"time w x(1+{inj} injet)", lambda d, inj=inj: kde(d, inj=inj))
    for at in (0.5, 2.0):
        run(f"time alpha {at}", lambda d, at=at: kde(d, at=at))
    run("time pT^0.5", lambda d: kde(d, pt_t=0.5))
    for Wt in (45., 75.):
        run(f"time W {Wt}", lambda d, Wt=Wt: kde(d, Wt=Wt))
    run("nit 6", lambda d: kde(d, nit=6))
