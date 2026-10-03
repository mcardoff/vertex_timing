"""Offline workbench: load a training export into flat numpy arrays grouped by event."""
import uproot, numpy as np, os, pickle, tempfile
D = os.environ.get("TRAINING_DIR", "/Users/mcard/project/vertex_timing/training_data/")
CACHE = os.environ.get("WB_CACHE", os.path.join(tempfile.gettempdir(), "vt_kde_cache"))
os.makedirs(CACHE, exist_ok=True)
CC = ['event_num', 'file_idx', 'cluster_idx', 'cluster_time', 'delta_z', 'cluster_d0_sigma', 'cluster_z_sigma', 'sumpt',
      'cluster_time_sigma', 'n_tracks', 't_in_jet', 'n_in_jet_timed', 'trkptz_score', 'waves_score', 'hgtd_time',
      'hgtd_valid', 'hgtd_time_res', 'delta_t', 'truth_purity', 'truth_n_hs_tracks', 'n_reco_vertices', 'reco_vtx_z',
      'n_forward_jets', 'frac_pt_in_fwdjet', 'time_chi2_ndf', 'local_vtx_density']
TC = ['event_num', 'file_idx', 'cluster_idx', 'pt', 'eta', 'phi', 'z0', 'sigma_z0', 'sigma_d0', 'd0', 'time', 'timeRes',
      'nhgtd_hits', 'z0_pull_pv', 'dr_nearest_fwdjet', 'pt_nearest_fwdjet', 'dr_nearest_anyjet', 'is_in_any_jet',
      'truth_is_hs', 'closer_to_pu_than_pv', 'quality']


def load(name, maxfiles=None, minfiles=0):
    cache = f"{CACHE}/{name}_{maxfiles}.pkl" if not minfiles else f"{CACHE}/{name}_{minfiles}_{maxfiles}.pkl"
    if os.path.exists(cache):
        return pickle.load(open(cache, 'rb'))
    f = uproot.open(D + name + "_training.root")
    out = {}
    for tn, cols in (('clusters', CC), ('tracks', TC)):
        t = f[tn]
        cut = None if maxfiles is None else f"(file_idx<{maxfiles})&(file_idx>={minfiles})"
        a = t.arrays(cols, cut=cut, library='np')
        key = a['file_idx'].astype(np.int64) * 10_000_000 + a['event_num'].astype(np.int64)
        o = np.argsort(key, kind='stable')
        a = {k: v[o] for k, v in a.items()}
        a['key'] = key[o]
        out[tn] = a
    ck = out['clusters']['key']
    uk = np.unique(ck)
    out['ekeys'] = uk
    out['clusters']['ev'] = np.searchsorted(uk, ck)
    tk = out['tracks']['key']
    ev = np.clip(np.searchsorted(uk, tk), 0, len(uk) - 1)
    ok = uk[ev] == tk
    assert ok.all(), ok.mean()
    out['tracks']['ev'] = ev
    pickle.dump(out, open(cache, 'wb'), protocol=4)
    return out


def evmax_pick(ev, score, nev):
    """index (into the row arrays) of the argmax-score row per event; -1 if none"""
    o = np.lexsort((score, ev))
    last = np.r_[ev[o][1:] != ev[o][:-1], True]
    pick = np.full(nev, -1, dtype=np.int64)
    pick[ev[o][last]] = o[last]
    return pick


if __name__ == "__main__":
    import sys
    for s in sys.argv[1:]:
        d = load(s)
        c = d['clusters']; t = d['tracks']; nev = len(d['ekeys'])
        print(s, nev, len(c['ev']), len(t['ev']))
        p = evmax_pick(c['ev'], c['trkptz_score'], nev)
        print(' TRKPTZ core', np.mean(np.abs(c['delta_t'][p]) < 60) * 100)
        print(' |z0pull*sig| median', np.nanmedian(np.abs(t['z0_pull_pv'] * t['sigma_z0'])))
        print(' tracks w/o cluster', np.mean(t['cluster_idx'] < 0), 'timeRes q', np.quantile(t['timeRes'], [.1, .5, .9]))
