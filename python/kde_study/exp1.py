"""Baseline TZP reproduction + failure anatomy on the offline exports."""
import sys, numpy as np
from wb import load, evmax_pick


def prep(name, maxfiles=None, minfiles=0):
    d = load(name, maxfiles, minfiles)
    c = d['clusters']; t = d['tracks']; nev = len(d['ekeys'])
    # cluster row index for each track: key (ev, cluster_idx)
    ckey = c['ev'].astype(np.int64) * 1000 + c['cluster_idx'].astype(np.int64)
    o = np.argsort(ckey); cs = ckey[o]
    tkey = t['ev'].astype(np.int64) * 1000 + t['cluster_idx'].astype(np.int64)
    pos = np.clip(np.searchsorted(cs, tkey), 0, len(cs) - 1)
    assert (cs[pos] == tkey).all()
    t['crow'] = o[pos]
    # per-event PV z and truth time
    vz = np.zeros(nev); vz[c['ev']] = c['reco_vtx_z']
    tt = np.zeros(nev); tt[c['ev']] = c['cluster_time'] - c['delta_t']
    d['vz'] = vz; d['ttruth'] = tt; d['nev'] = nev
    t['dz'] = np.abs(t['z0'] - vz[t['ev']])
    return d


def tzp_score(d, a=1.0, b=0.6, g=0.45):
    c = d['clusters']; t = d['tracks']; nc = len(c['ev'])
    s = np.bincount(t['crow'], t['pt'] * np.exp(-a * t['dz']), nc)
    info = np.bincount(t['crow'], 1.0 / t['sigma_d0'] ** 2, nc)
    return s * np.exp(-b * np.abs(c['delta_z'])) * info ** (0.5 * g)


def tzp_time(c, guard=3):
    use = (c['n_in_jet_timed'] >= guard) & np.isfinite(c['t_in_jet'])
    return np.where(use, c['t_in_jet'], c['cluster_time'])


def core(d, pick, tsel):
    """pick: cluster row per event; tsel: per-cluster-row time"""
    c = d['clusters']
    ok = pick >= 0
    dt = np.full(d['nev'], 1e9); dt[ok] = tsel[pick[ok]] - d['ttruth'][ok]
    return np.abs(dt) < 60


if __name__ == "__main__":
    for s in sys.argv[1:]:
        name, _, mf = s.partition(':')
        d = prep(name, int(mf) if mf else None)
        c = d['clusters']; nev = d['nev']
        S = tzp_score(d); T = tzp_time(c)
        p = evmax_pick(c['ev'], S, nev)
        ok = core(d, p, T)
        okraw = core(d, p, c['cluster_time'])
        print(f"{s}: nev {nev}  TZP {ok.mean()*100:.2f}  (full-cluster time {okraw.mean()*100:.2f})")
        # anatomy
        cl_ok = np.abs(T - d['ttruth'][c['ev']]) < 60
        any_ok = np.bincount(c['ev'], cl_ok, nev) > 0
        fail = ~ok
        # HS-dominant cluster per event = max truth_n_hs_tracks
        hsbest = evmax_pick(c['ev'], c['truth_n_hs_tracks'] + 1e-3 * c['truth_purity'], nev)
        nhs_ev = np.bincount(c['ev'], c['truth_n_hs_tracks'], nev)
        picked_hs = (p == hsbest)
        print(f"  oracle any-cluster-passes {any_ok.mean()*100:.2f}")
        print(f"  fails: {fail.mean()*100:.2f}% of events; of them: another cluster passes {np.mean(any_ok[fail])*100:.1f}%, "
              f"picked IS the max-HS cluster {np.mean(picked_hs[fail])*100:.1f}%, no HS timed track in event {np.mean(nhs_ev[fail]==0)*100:.1f}%")
        for lo, hi in ((0, 0), (1, 1), (2, 3), (4, 7), (8, 999)):
            m = (nhs_ev >= lo) & (nhs_ev <= hi)
            print(f"   nHS timed {lo}-{hi}: {m.mean()*100:5.1f}% of ev, core {ok[m].mean()*100:5.1f}, share of fails {fail[m].sum()/fail.sum()*100:5.1f}%, "
                  f"any-ok {any_ok[m].mean()*100:5.1f}, maxHS-cluster core {core(d, hsbest, T)[m].mean()*100:5.1f}")
