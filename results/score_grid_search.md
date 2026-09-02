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

---

# Out-of-sample check (`python/oos_check.py`)

## mjj500 (mu=200) -- the VBS-enriched selection

Honest framing: these events OVERLAP the novbs tuning sets (mjj500 is stricter,
so a subset), but the population mix is very different -- this tests survival
under reweighting toward VBS topology, not statistical independence.

| sample | events | TRKPTZ | TZP | WAVES | TZP-TRKPTZ | +- |
|---|---:|---:|---:|---:|---:|---:|
| vbf | 519,115 | 90.45 | **92.08** | 91.91 | **+1.63** | 0.06 |
| zjets | 36,531 | 61.87 | **64.31** | 57.94 | **+2.44** | 0.36 |
| dijet | 80,383 | 86.88 | **88.58** | 86.08 | **+1.70** | 0.16 |
| ttbar | 561,745 | 87.08 | **88.84** | 86.08 | **+1.77** | 0.06 |

Gains of the same size as on the tuning selection (+1.6 to +2.4), and TZP beats
deployed WAVeS on ALL FOUR (vbf +0.17, zjets +6.37, dijet +2.50, ttbar +2.76).
On the canonical zjets selection specifically -- the headline population of the
whole effort -- TZP reaches **64.31% against the 61.87% TRKPTZ baseline**, a
purely classical +2.44.

## mu0 (novbs) -- disjoint zero-pileup no-harm controls

| sample | events | TRKPTZ | TZP | delta |
|---|---:|---:|---:|---:|
| vbf_mu0 | 2,741,680 | 99.87 | 99.86 | -0.01 |
| zeejets_mu0 | 34,310 | 99.76 | 99.76 | -0.01 |
| ttbar_mu0 | 2,128,740 | 99.87 | 99.86 | -0.01 |

4.9M genuinely disjoint events: deltas of -0.01 everywhere, i.e. ~300 events out
of 4.9M. The full stack of changes (per-track z-weight, softened cluster term,
d0 power, guarded in-jet timing) is inert where there is nothing to gain.
Deployed WAVeS by contrast costs 0.3-0.8 even at mu0.

## Verdict

The accumulated-tuning concern is closed to the extent local data allows: the
constants survive a strong reweighting with undiminished gains, and do no harm
on disjoint zero-pileup data. What this does NOT test: a statistically
independent mu=200 dataset. That would need new MC.
