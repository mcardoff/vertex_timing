"""Track-list study on a dump_tracks file: where do timed hard-scatter tracks go, and what does a different
selection do to the full TZP_KDE_TZ chain (reliability weights on)?
   wide.py <trackdump.root> [max_file_idx]"""
import sys, numpy as np, uproot
from wb import evmax_pick
import reclus
from sel import cands

NSIG0, PTLO0, PTHI0 = 3.0, 1.0, 30.0


def load(path, maxfile=None):
    f = uproot.open(path)
    cut = None if maxfile is None else f"file_idx<{maxfile}"
    # the dump is very loose (345 tracks/event); keep what the scans below can use, or memory explodes
    tcut = cut
    e = f['events'].arrays(library='np', cut=cut); t = f['tracks'].arrays(library='np', cut=tcut)
    ek = e['file_idx'].astype(np.int64) * 10_000_000 + e['event_num'].astype(np.int64)
    keep = e['n_nominal'] > 0                       # the nominal event population
    o = np.argsort(ek[keep]); uk = ek[keep][o]
    ev = {k: v[keep][o] for k, v in e.items()}
    tk = t['file_idx'].astype(np.int64) * 10_000_000 + t['event_num'].astype(np.int64)
    pos = np.clip(np.searchsorted(uk, tk), 0, len(uk) - 1); ok = uk[pos] == tk
    o2 = np.argsort(pos[ok], kind='stable')
    tr = {k: v[ok][o2] for k, v in t.items()}; tr['ev'] = pos[ok][o2]
    tr['dz'] = np.abs(tr['z0'] - ev['vz'][tr['ev']])
    return ev, tr


def make_d(ev, tr, mask):
    t = {k: v[mask] for k, v in tr.items()}
    nev = len(ev['vz'])
    d = {'tracks': t, 'nev': nev, 'vz': ev['vz'].astype(float), 'ttruth': ev['ttruth'].astype(float)}
    d['off'] = np.searchsorted(t['ev'], np.arange(nev + 1)).astype(np.int64)
    for k in ('time', 'timeRes', 'z0', 'sigma_z0', 'pt', 'sigma_d0', 'dz'):
        t[k + '_d'] = np.ascontiguousarray(t[k], dtype=np.float64)
    t['injet'] = t['dr_nearest_fwdjet'] < 0.4
    t['w_tzp'] = t['pt_d'] * np.exp(-t['dz_d'])
    return d


def nominal(tr, nsig=NSIG0, ptlo=PTLO0, pthi=PTHI0, quality=True):
    m = (tr['nsig'] < nsig) & (tr['pt'] > ptlo) & (tr['pt'] < pthi)
    return m & (tr['quality'] > 0.5) if quality else m


def evaluate(ev, tr, mask, zscale=2.0, hit=0.5, extra=None):
    d = make_d(ev, tr, mask); t = d['tracks']
    zf = 1 / (1 + t['sigma_z0_d'] / zscale)
    rel = np.where(t['nhgtd_hits'] >= 2, 1.0, hit) * zf
    if extra is not None:
        e = extra(t); rel = rel * e; zf = zf * e
    c = cands(d, sfac=rel, tfac=zf)
    p = evmax_pick(c['ev'], c['S'], d['nev'])
    ok = np.zeros(d['nev'], bool); has = p >= 0
    ok[has] = c['ok'][p[has]]
    return ok


