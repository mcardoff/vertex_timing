"""Re-cluster the exported tracks with an arbitrary algorithm, rebuild every cluster-level quantity from the
tracks alone, and evaluate TZP and the kernel-density score (no Athena term) on the new partition."""
import ctypes, os, subprocess, numpy as np
from wb import evmax_pick
from exp1 import prep

HERE = os.path.dirname(os.path.abspath(__file__))
LIB = os.path.join(HERE, "libclus.dylib")
if not os.path.exists(LIB) or os.path.getmtime(LIB) < os.path.getmtime(os.path.join(HERE, "clus.cpp")):
    subprocess.check_call(["clang++", "-O3", "-std=c++17", "-shared", "-fPIC", "-o", LIB, os.path.join(HERE, "clus.cpp")])
_lib = ctypes.CDLL(LIB)
_dp = np.ctypeslib.ndpointer(np.float64, flags="C")
_lib.cluster.argtypes = [ctypes.c_int, np.ctypeslib.ndpointer(np.int64, flags="C"), _dp, _dp, _dp, _dp, _dp, _dp,
                         ctypes.c_int, ctypes.c_double, ctypes.c_int, ctypes.c_double, ctypes.c_double,
                         np.ctypeslib.ndpointer(np.int32, flags="C")]

SAMPLES = [('zjets_novbs', None), ('dijet_novbs', None), ('vbf_novbs', 120), ('ttbar_novbs', 60)]
ITER, SIMUL, CONE, MSHIFT, WINDOW = 0, 1, 2, 3, 4


def setup(name, mx=None, mn=0):
    d = prep(name, mx, mn); t = d['tracks']; nev = d['nev']
    assert (np.diff(t['ev']) >= 0).all()
    d['off'] = np.searchsorted(t['ev'], np.arange(nev + 1)).astype(np.int64)
    for k in ('time', 'timeRes', 'z0', 'sigma_z0', 'pt', 'sigma_d0', 'dz'):
        t[k + '_d'] = np.ascontiguousarray(t[k], dtype=np.float64)
    t['injet'] = np.nan_to_num(t['dr_nearest_fwdjet'], nan=9.) < 0.4
    t['w_tzp'] = t['pt_d'] * np.exp(-t['dz_d'])
    return d


def cluster(d, method=ITER, cut=3.0, seedw=None, kw=None, useZ=False, floor=0.0, par=0.0, szscale=1.0, szfloor=0.0):
    t = d['tracks']
    sz = np.sqrt((szscale * t['sigma_z0_d']) ** 2 + szfloor ** 2)
    seedw = t['pt_d'] if seedw is None else np.ascontiguousarray(seedw, dtype=np.float64)
    kw = t['w_tzp'] if kw is None else np.ascontiguousarray(kw, dtype=np.float64)
    lab = np.zeros(len(t['ev']), np.int32)
    _lib.cluster(d['nev'], d['off'], t['time_d'], t['timeRes_d'], t['z0_d'], sz, seedw, kw,
                 method, cut, int(useZ), floor, par, lab)
    return lab


def features(d, lab, b=0.6, g=0.45, zw='prec'):
    """cluster table from a labelling: dict of per-cluster arrays + per-track cluster row"""
    t = d['tracks']
    key = t['ev'].astype(np.int64) * 4096 + lab
    uk, crow = np.unique(key, return_inverse=True)
    nc = len(uk); ev = (uk // 4096).astype(np.int64)
    iv = 1 / t['timeRes_d'] ** 2
    c = {'ev': ev, 'n': np.bincount(crow, minlength=nc)}
    c['time'] = np.bincount(crow, iv * t['time_d'], nc) / np.bincount(crow, iv, nc)
    izv = 1 / t['sigma_z0_d'] ** 2
    if zw == 'pt': izv = t['pt_d']
    elif zw == 'ptprec': izv = t['pt_d'] * izv
    c['dz'] = np.bincount(crow, izv * t['z0_d'], nc) / np.bincount(crow, izv, nc) - d['vz'][ev]
    c['info'] = np.bincount(crow, 1 / t['sigma_d0_d'] ** 2, nc)
    nin = np.bincount(crow, t['injet'], nc)
    with np.errstate(invalid='ignore', divide='ignore'):
        tin = np.bincount(crow, iv * t['time_d'] * t['injet'], nc) / np.bincount(crow, iv * t['injet'], nc)
    c['t_tzp'] = np.where(nin >= 3, tin, c['time'])
    c['env'] = np.exp(-b * np.abs(c['dz'])) * c['info'] ** (0.5 * g)
    c['S_tzp'] = np.bincount(crow, t['w_tzp'], nc) * c['env']
    return c, crow


def kde_eval(d, c, crow, W=60., p=0.5, nit=3, zk=0.0, zk_time=False, zfloor=0.0):
    """kernel-density score + mean-shift time for every cluster (pairs built per call)"""
    t = d['tracks']; nc = len(c['ev']); nev = d['nev']
    cnt = (d['off'][1:] - d['off'][:-1])[c['ev']]
    pc = np.repeat(np.arange(nc), cnt)
    pt_ = np.repeat(d['off'][:-1][c['ev']], cnt) + np.arange(cnt.sum()) - np.repeat(np.cumsum(cnt) - cnt, cnt)
    base = t['pt_d'] ** p * np.exp(-t['dz_d'])
    wt = base / t['timeRes_d'] ** 2
    m = np.bincount(crow, wt * t['time_d'], nc) / np.bincount(crow, wt, nc)
    tt = t['time_d'][pt_]; w = wt[pt_]
    kz = 1.0
    if zk > 0:  # z-compatibility of each track with the cluster's z centroid
        zc = c['dz'] + d['vz'][c['ev']]
        kz = np.exp(-0.5 * (t['z0_d'][pt_] - zc[pc]) ** 2 / ((zk * t['sigma_z0_d'][pt_]) ** 2 + zfloor ** 2))
        if zk_time: w = w * kz
    for _ in range(nit):
        k = w * np.exp(-0.5 * ((tt - m[pc]) / W) ** 2)
        den = np.bincount(pc, k, nc); num = np.bincount(pc, k * tt, nc)
        m = np.where(den > 0, num / np.where(den > 0, den, 1), m)
    dens = np.bincount(pc, base[pt_] * kz * np.exp(-0.5 * ((tt - m[pc]) / W) ** 2), nc)
    ncl = np.bincount(c['ev'], minlength=nev)
    tk = np.where(ncl[c['ev']] <= 1, c['t_tzp'], m)
    return dens * c['env'], tk


def ok_of(d, c, S, T):
    p = evmax_pick(c['ev'], S, d['nev'])
    return np.abs(T[p] - d['ttruth']) < 60


def evaluate(d, lab, kde=True, feat=None, **kk):
    c, crow = features(d, lab, **(feat or {}))
    out = {'tzp': ok_of(d, c, c['S_tzp'], c['t_tzp']), 'ncl': len(c['ev']) / d['nev']}
    if kde:
        S, T = kde_eval(d, c, crow, **kk)
        out['kde'] = ok_of(d, c, S, T)
    return out


_DS = {}


def get(n, m):
    if n not in _DS:
        d = setup(n, m)
        d['ref'] = evaluate(d, cluster(d))
        _DS[n] = d
    return _DS[n]


def run(label, kde_kw=None, **kw):
    """print % of the production-TZP fails removed by (TZP on new clusters) and (KDE on new clusters)"""
    a, b, ncl = [], [], []; kde_kw = kde_kw or {}
    for n, m in SAMPLES:
        d = get(n, m)
        r = evaluate(d, cluster(d, **{k: (v(d) if callable(v) else v) for k, v in kw.items()}), **kde_kw)
        f0 = 1 - d['ref']['tzp'].mean()
        a.append((r['tzp'].mean() - d['ref']['tzp'].mean()) / f0 * 100)
        b.append((r['kde'].mean() - d['ref']['tzp'].mean()) / f0 * 100)
        ncl.append(r['ncl'])
    print(f"{label:44s} TZP: " + " ".join(f"{x:+6.2f}" for x in a) + "  | KDE: " + " ".join(f"{x:+6.2f}" for x in b)
          + f"  | worst {min(b):+.2f} mean {np.mean(b):+.2f} | ncl {ncl[0]:.1f}/{ncl[2]:.1f}")
    return a, b


if __name__ == "__main__":
    for n, m in SAMPLES:
        d = get(n, m); c0 = d['clusters']
        lab = cluster(d)
        same = np.mean(np.bincount(d['tracks']['ev'][lab != d['tracks']['cluster_idx'].astype(int)], minlength=d['nev']) == 0)
        print(f"{n}: production partition reproduced exactly in {same*100:.2f}% of events (label-for-label); "
              f"TZP on re-clustered {d['ref']['tzp'].mean()*100:.2f}, KDE {d['ref']['kde'].mean()*100:.2f}")
