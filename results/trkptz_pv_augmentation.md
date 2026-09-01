# Augmenting TRKPTZ with `closer_to_pu_than_pv` — local VBF, 2026-09-01

Proposal: restrict TRKPTZ's pT sum with the per-track vertex-proximity flag,

    score = exp(-1.5 |dz_cluster|) * SUM_tracks { pT if <flag matches>, else 0 }

with `RecoVtx_z[0]` as the primary vertex. Implemented in the main program as
three registry scores (`src/clustering_constants.h`, ids 24/25/26), computed in
`Cluster::updateScores` beside TRKPTZ itself, sharing the identical `exp(-1.5|dz|)`
factor so ONLY the pT sum differs. Timing is untouched: none of the three appears
in `Cluster::calculateTime`'s score list, so all report the standard
inverse-variance mean, and the comparison is pure selection.

Flag helper: `BranchPointerWrapper::closerToPuThanPv` (`src/clustering_structs.h`),
mirroring `vertexRelation()` in `util/export_training_data.cxx`. Reco-only —
`RecoVtx_z[0]` is the primary, every other entry is pileup, no truth consulted.

## Result — `clustering_hist` + `clustering_plot`, local VBF (45,335 events)

Core fraction at PASS_SIGMA = 60 ps, summed over all n-forward-HS-track bins.

| score | core frac | vs TRKPTZ |
|---|---:|---:|
| WAVES | 94.11% | +1.08 |
| HGTD | 93.32% | +0.29 |
| **TRKPTZ_PUW** — PU-side down-weighted, w = 0.4 | **93.31%** | **+0.28** |
| TRKPTZ — baseline | 93.04% | — |
| **TRKPTZ_PV** — hard veto, as proposed | **91.76%** | **-1.27** |
| TRKPTZ_PU — orientation control | 84.00% | -9.03 |

**The hard veto as specified costs 1.27 points.** The orientation control settles
that this is not a sign error: keeping the PU-side tracks instead is 9 points
worse, so the flag carries real information and `_PV` is the correct reading of
it. The veto is simply too blunt to use.

## Why it fails (local VBF, 15.25M timed tracks)

| quantity | value |
|---|---:|
| P(flag = 1 \| genuine HS track) | **46.8%** |
| P(flag = 1 \| PU track) | 88.0% |
| clusters vetoed to exactly zero | **56.0%** |
| median cluster pT surviving the veto | 0.0% |
| events where every cluster scores zero (argmax arbitrary) | 1.02% |

The flag removes 88% of pileup but **also 47% of the hard scatter**. TRKPTZ ranks
clusters by hard-scatter pT content, so destroying half of that signal costs more
than the pileup suppression buys — especially since `exp(-1.5|dz|)` already
penalises off-vertex clusters, which is the same information in smoother form.

The degeneracy is a second, smaller problem: a cluster with every track flagged
scores exactly 0, and in 1.02% of events that is true of *every* cluster, leaving
`chooseCluster` to keep `collection[0]` — an arbitrary pick.

## Softening it does help

Replacing the veto with a down-weight, `SUM_PV pT + w * SUM_PU pT`, so w = 1 is
TRKPTZ and w = 0 is the veto. Scanned in python against the exported local VBF
(m_jj >= 500; TRKPTZ reproduced from `sumpt` and `delta_z` to max relative error
9.6e-05 over 2.93M clusters, and `sumpt` from the track table to 1.2e-07, before
trusting any of it):

| w | 0.0 | 0.1 | 0.2 | 0.3 | **0.4** | 0.5 | 0.6 | 0.8 | 1.0 |
|---|---|---|---|---|---|---|---|---|---|
| core frac | 89.22 | 90.34 | 90.65 | 90.74 | **90.76** | 90.73 | 90.70 | 90.58 | 90.45 |

Smooth and single-peaked, +0.31 over the w = 1 baseline, flat to within 0.03
across w in [0.2, 0.6]. `PU_SIDE_PT_WEIGHT = 0.4` is now that peak.

Two independent confirmations that this is not a scan artefact: the python
w = 0 point (-1.23) reproduces the C++ `TRKPTZ_PV` row (-1.27) across a
different selection and a separate implementation, and the C++ `TRKPTZ_PUW` row
lands at +0.28 against the python scan's predicted +0.31.

## Standing caveats

- **Local VBF only.** The weight has NOT been re-scanned on zjets/dijet/ttbar,
  and zjets is where the headroom is. Do that before treating 0.4 as general.
- **+0.28 is well below WAVeS's +1.08** on the same sample, so this is a cheap
  refinement of the baseline rather than a competitive selector.
- The degeneracy at w = 0 disappears for any w > 0 (the PU-side sum is never
  identically zero for a non-empty cluster), so the down-weighted form has no
  arbitrary-argmax population at all.
