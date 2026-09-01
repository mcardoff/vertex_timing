# Joint (alpha, beta, gamma) grid search for the TZP score

`python/score_grid_search.py` + edge extension. The score

    S = [SUM_t pT e^{-alpha |z0_t - z_PV|}] * e^{-beta |z_C - z_PV|} * (1/sigma_d0)^gamma

had alpha = 0.7 and gamma = 1.2 from 1-D tunes and **beta = 1.5 inherited from
original TRKPTZ, never scanned by anyone**. 5x5x5 grid + extensions past both
edges, all four novbs samples, objective = worst-sample gain (one global point).

## The joint optimum: (1.0, 0.6, 0.45)

| | vbf | zjets | dijet | ttbar | worst |
|---|---:|---:|---:|---:|---:|
| vs (0.7, 1.5, 1.2) | +0.16 | **+0.38** | +0.15 | +0.15 | **+0.14** |
| absolute | 91.66% | 64.71% | 88.97% | 89.08% | |

Interior on every axis, on a plateau:

    alpha: 0.7:+0.08   1.0:+0.14   1.4:+0.08          (clean peak)
    beta:  0.15:-0.24  0.3:+0.00  0.45:+0.11  0.6:+0.14  0.75:+0.14  1.0:+0.06
    gamma: 0.15:+0.11  0.3:+0.14  0.45:+0.14  0.6:+0.12  0.9:+0.09   (flat)

## The physics of the rebalance

The three terms overlap, so joint tuning trades them off: a stronger per-track
z-term (0.7 -> 1.0) makes the steep cluster-level z-term partly redundant
(1.5 -> 0.6, a 2.5x softening of a constant that predates this whole effort),
and less precision weight is then needed (1.2 -> 0.45). The 1-D optima were
each conditioned on the others' old values.

## The sample-independence constraint costs almost nothing

Per-sample optima buy +0.14 (ttbar) to +0.41 (zjets) over the current point --
i.e. a per-sample tune would gain at most ~0.3 more than the global one. The
constraint is cheap here.

## Implementation

TZP-only constants `TZP_TRACK_DZ_WEIGHT = 1.0`, `TZP_CLUSTER_DZ_WEIGHT = 0.6`,
`TZP_D0_PRECISION = 0.45`; `TRACK_DZ_WEIGHT` and the 1.5 are untouched so
TRKPTZ and the historical ladder rows (TZ, TZ_GIJ) are bit-identical. C++
validation on local VBF (mjj>=200, all 48,314 events):

| | core frac | vs TRKPTZ |
|---|---:|---:|
| **TZP, joint tune** | **91.99%** | **+1.64** |
| TZP, previous (0.7, 1.5, 1.2) | 91.83% | +1.47 |
| WAVES (deployed) | 91.71% | +1.36 |
| TZ_GIJ | 91.58% | +1.22 |
| TRKPTZ | 90.36% | -- |

The +0.16 step matches the python prediction exactly; WAVES reproduces at 91.71.

## Standing caveat, restated

This is the fourth round of in-sample optimisation on the same four novbs
samples. Every choice is plateau-supported and uniform across samples, which is
the best available internal evidence -- but a genuinely out-of-sample
confirmation (mjj500 selections, mu0 no-harm controls) is the cheap way to close
the accumulated-tuning concern before this is quoted anywhere.
