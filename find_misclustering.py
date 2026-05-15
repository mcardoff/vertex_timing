#!/usr/bin/env python3
"""
find_misclustering.py
Scans events to find extreme misclustering cases where HS tracks are split
across multiple clusters and no cluster has high purity.

Uses uproot.iterate + awkward arrays to vectorize event pre-filtering;
the clustering itself is per-event (inherently sequential).
"""

import numpy as np
import uproot
import awkward as ak
from pathlib import Path

# ── Mirror C++ constants ──────────────────────────────────────────────────────
DIST_CUT           = 3.0
MIN_PT             = 1.0
MAX_PT             = 30.0
MIN_ABS_ETA        = 2.38
MAX_ABS_ETA        = 4.0
MAX_NSIGMA_Z       = 3.0
MIN_CLUSTER_TRACKS = 3
MIN_FWD_JET_PT     = 30.0

FILES = sorted(Path('../ntuple-hgtd').glob('*.root'))

BRANCHES = [
    'RecoVtx_z',
    'TruthVtx_isHS', 'TruthVtx_time',
    'AntiKt4EMTopoJets_pt', 'AntiKt4EMTopoJets_eta',
    'Track_z0', 'Track_var_z0', 'Track_pt', 'Track_eta',
    'Track_time', 'Track_timeRes',
    'Track_hasValidTime', 'Track_quality', 'Track_truthVtx_idx',
]

# Max allowed |cluster_time - truth_HS_time| for HS-containing clusters.
# Clusters further than this are flagged as misassignment, not misclustering.
MAX_MISASSIGN_DT = 150.0  # ps

# ── Clustering (per-event, inherently sequential) ─────────────────────────────
def cluster_event(tracks):
    """tracks: list of dicts with t, tres, pt, is_hs"""
    N = len(tracks)
    if N == 0:
        return []

    C = [dict(t=tr['t'], tres=tr['tres'], pt=tr['pt'],
              hs_pt=tr['pt'] if tr['is_hs'] else 0.0,
              n=1)
         for tr in tracks]

    def tdist(a, b):
        return abs(a['t'] - b['t']) / np.sqrt(a['tres']**2 + b['tres']**2)

    def merge(a, b):
        wa, wb = 1/a['tres']**2, 1/b['tres']**2
        return dict(
            t    = (a['t']*wa + b['t']*wb) / (wa + wb),
            tres = np.sqrt(1/(wa + wb)),
            pt   = a['pt'] + b['pt'],
            hs_pt= a['hs_pt'] + b['hs_pt'],
            n    = a['n'] + b['n'],
        )

    done   = [False] * N
    finals = []

    while True:
        available = [i for i in range(N) if not done[i]]
        if not available:
            break
        seed = max(available, key=lambda i: C[i]['pt'])
        while True:
            neighbors = [j for j in range(N) if j != seed and not done[j]]
            if not neighbors:
                break
            nn = min(neighbors, key=lambda j: tdist(C[seed], C[j]))
            if tdist(C[seed], C[nn]) < DIST_CUT:
                C[seed] = merge(C[seed], C[nn])
                done[nn] = True
            else:
                break
        finals.append(C[seed])
        done[seed] = True

    return [c for c in finals if c['n'] >= MIN_CLUSTER_TRACKS]


def purity(c):
    return c['hs_pt'] / c['pt'] if c['pt'] > 0 else 0.0


candidates = []

for fpath in FILES:
    fnum = fpath.stem.split('._')[1].split('.')[0]
    print(f'Scanning {fnum} ...', flush=True)

    n_scanned = 0
    for batch, report in uproot.iterate(
            f'{fpath}:ntuple', BRANCHES,
            step_size=500, report=True):

        offset = report.start   # global entry offset for this batch

        # ── Vectorized event pre-filter ───────────────────────────────────────
        # 1. Need at least one reco vertex
        has_vtx = ak.num(batch['RecoVtx_z']) > 0

        # 2. Need >= 1 forward jet with pT > 30 GeV
        jet_pt  = batch['AntiKt4EMTopoJets_pt']
        jet_eta = batch['AntiKt4EMTopoJets_eta']
        fwd_jet = (np.abs(jet_eta) > MIN_ABS_ETA) & (np.abs(jet_eta) < MAX_ABS_ETA) & (jet_pt > MIN_FWD_JET_PT)
        has_fwd_jet = ak.any(fwd_jet, axis=1)

        # 3. Need >= 2 HGTD tracks (rough: valid time + forward eta)
        trk_eta   = batch['Track_eta']
        trk_valid = batch['Track_hasValidTime'] == 1
        trk_fwd   = (np.abs(trk_eta) > MIN_ABS_ETA) & (np.abs(trk_eta) < MAX_ABS_ETA)
        n_hgtd    = ak.sum(trk_fwd & trk_valid, axis=1)
        has_tracks = n_hgtd >= 2

        keep = has_vtx & has_fwd_jet & has_tracks
        indices = np.where(ak.to_numpy(keep))[0]

        n_scanned += len(batch['RecoVtx_z'])

        for local_i in indices:
            evt = offset + local_i

            reco_z = float(batch['RecoVtx_z'][local_i][0])
            is_hs_vtx   = ak.to_numpy(batch['TruthVtx_isHS'][local_i]).astype(bool)
            truth_hs_t  = float(ak.to_numpy(batch['TruthVtx_time'][local_i])[is_hs_vtx][0])

            # Per-track numpy arrays
            eta   = ak.to_numpy(batch['Track_eta'][local_i])
            pt    = ak.to_numpy(batch['Track_pt'][local_i])
            z0    = ak.to_numpy(batch['Track_z0'][local_i])
            vz0   = ak.to_numpy(batch['Track_var_z0'][local_i])
            t     = ak.to_numpy(batch['Track_time'][local_i])
            tres  = ak.to_numpy(batch['Track_timeRes'][local_i])
            valid = ak.to_numpy(batch['Track_hasValidTime'][local_i]) == 1
            qual  = ak.to_numpy(batch['Track_quality'][local_i]) == 1
            tvtx  = ak.to_numpy(batch['Track_truthVtx_idx'][local_i]).astype(int)
            is_hs = is_hs_vtx[tvtx]

            # Vectorized track selection
            in_hgtd  = (np.abs(eta) > MIN_ABS_ETA) & (np.abs(eta) < MAX_ABS_ETA)
            nsigma_z = np.abs(z0 - reco_z) / np.sqrt(np.maximum(vz0, 1e-12))
            sel = in_hgtd & (pt >= MIN_PT) & (pt < MAX_PT) & (nsigma_z < MAX_NSIGMA_Z) & qual & valid

            idx = np.where(sel)[0]
            if len(idx) == 0:
                continue

            n_hs = int(np.sum(is_hs[idx]))
            if n_hs < 5:
                continue

            tracks = [dict(t=t[i], tres=tres[i], pt=pt[i], is_hs=bool(is_hs[i]))
                      for i in idx]

            finals = cluster_event(tracks)
            if not finals:
                continue

            hs_clusters = [c for c in finals if c['hs_pt'] > 0]
            if len(hs_clusters) < 2:
                continue

            # Reject misassignment: if any HS-containing cluster is far from
            # the truth HS time, the split is due to bad timing, not geometry.
            if any(abs(c['t'] - truth_hs_t) > MAX_MISASSIGN_DT for c in hs_clusters):
                continue

            purities  = [purity(c) for c in hs_clusters]
            max_pur   = max(purities)

            candidates.append(dict(
                file=fnum, evt=evt,
                max_purity=max_pur,
                n_hs_clusters=len(hs_clusters),
                n_hs_tracks=n_hs,
                purities=sorted(purities, reverse=True),
                total_clusters=len(finals),
            ))

    print(f'  scanned {n_scanned} events, {len(candidates)} candidates so far')

# Sort: lowest max purity first, then most HS tracks, then most HS clusters
candidates.sort(key=lambda x: (x['max_purity'], -x['n_hs_tracks'], -x['n_hs_clusters']))

print(f'\n{"="*72}')
print(f'Top 50 extreme misclustering candidates')
print(f'{"="*72}')
print(f'{"File":>8}  {"Evt":>6}  {"MaxPur":>7}  {"#HSClust":>8}  '
      f'{"#HSTrk":>7}  {"#Clust":>7}  Purities')
print(f'{"-"*72}')
for c in candidates[:50]:
    pur_str = '  '.join(f'{p:.2f}' for p in c['purities'])
    print(f'{c["file"]:>8}  {c["evt"]:>6}  {c["max_purity"]:>7.3f}  '
          f'{c["n_hs_clusters"]:>8}  {c["n_hs_tracks"]:>7}  '
          f'{c["total_clusters"]:>7}  [{pur_str}]')

if candidates:
    best = candidates[0]
    print(f'\nBest candidate commands:')
    print(f'  python3 clustering_animation.py --file_num {best["file"]} --event_num {best["evt"]}')
    print(f'  python3 event_display.py --file_num {best["file"]} --event_num {best["evt"]} --extra_time 0.0')
