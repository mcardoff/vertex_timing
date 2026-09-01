# Clustering in-jet and out-of-jet tracks separately — Z+jets

`python/split_clustering.py`. Everything before this built ONE collection from
all tracks and subset the time afterwards. This runs the clustering twice on
disjoint track sets, giving two independent collections.

**The clustering is a validated reimplementation**, not an approximation:
`doIterativeClustering` is 1-D in time here (seed on highest-pT unconsumed track,
absorb the nearest unconsumed original cluster while
`|dt|/sqrt(s_a^2+s_b^2) < 3.0`, recompute the centroid after each absorption).
Run on the full track set it reproduces the C++ partition **3000/3000 events
exactly**, and the script refuses to report subset results below 98%.

## As a replacement estimator: no

| collection | clusters/event | core frac | sel. purity | n trk in pick |
|---|---:|---:|---:|---:|
| all tracks (baseline) | 5.60 | **63.89%** | 42.24% | 10.3 |
| out-of-jet only | 5.52 | 60.75% (**-3.14**) | 39.31% | 9.6 |
| in-jet only | 1.47 | **25.34%** | 39.37% | **1.8** |

The in-jet collection collapses because Z+jets has only ~1.3 in-jet tracks per
event (5.0% of ~26). There is nothing to cluster: the selected "cluster"
averages 1.8 tracks, so its time is essentially one track's HGTD measurement.

Oracle over all three collections **69.39%** against a matched-noise null of
**72.00%** -- real complementarity **-2.61**. Combining the collections, even
with truth, buys nothing beyond the best-of-N statistic.

## As a confidence flag: yes, and it is the strongest reco-only quality
## variable found so far

Only **53.8%** of Z+jets events have any in-jet track, so `t_in` is undefined for
the rest. Restricting to events where the comparison is actually defined:

| | fraction of defined events | core frac |
|---|---:|---:|
| `\|t_in - t_out\| < 30 ps` | 26.0% | **81.60%** |
| `\|t_in - t_out\| < 60 ps` | 40.1% | **79.42%** |
| disagree | 59.9% | 58.36% |
| *(events with no in-jet track at all)* | *46.2% of all* | *60.51%* |

**Gap +21.05 points**, and it is not multiplicity in disguise -- the separation
holds inside every track-multiplicity bin:

| n timed trk | events | agree | core \| agree | core \| disagree | gap |
|---|---:|---:|---:|---:|---:|
| 0-9 | 8163 | 8.5% | 88.33% | 61.32% | **+27.01** |
| 10-19 | 22548 | 17.1% | 84.05% | 61.31% | +22.73 |
| 20-29 | 29240 | 21.3% | 79.67% | 59.16% | +20.52 |
| 30-44 | 26936 | 26.9% | 77.35% | 58.81% | +18.54 |
| >=45 | 7090 | 31.5% | 74.63% | 56.02% | +18.61 |

Two independent collections agreeing on a time is strong evidence the time is
right, even though neither collection is individually better than clustering
everything together. That is a **confidence** statement, not an estimator
improvement -- consistent with the oracle finding no complementarity.

## Why this is worth having

Nothing else measured in this work identifies a 40%-of-events subset at 79% core
fraction against a 62% average, from reco alone. The obvious consumer is the RpT
side: the documented efficiency ceiling there comes from applying the time gate
to jets whose t0 is untrustworthy (a gate that empties a jet sends it to
R_pT = 0 and it becomes indistinguishable from pileup). A per-event confidence
flag lets the gate be applied only where it is earned.

Caveats: defined for 53.8% of Z+jets events only; VBF not yet run; and the flag
costs a second clustering pass, though on the in-jet subset that is ~1.3 tracks
and therefore nearly free.
