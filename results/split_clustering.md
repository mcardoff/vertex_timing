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

---

# VBF (local export, 39,681 events) — the same study

Clustering validation again **3000/3000 exact**.

## As a replacement estimator

| collection | clusters/event | core frac | sel. purity | n trk in pick |
|---|---:|---:|---:|---:|
| all tracks | 5.64 | **91.49%** | 72.66% | 13.8 |
| in-jet only | 1.77 | 86.38% | **86.03%** | 5.2 |
| out-of-jet only | 5.39 | 66.00% | 44.70% | 9.4 |

Still no replacement -- all-track wins -- but the in-jet collection is *viable*
here (86.38% at 86% purity on 5.2 tracks) where in Z+jets it collapsed to 25.34%
on 1.8 tracks. And the ordering INVERTS between samples: in VBF the in-jet
collection is the good one and out-of-jet is the poor one (66.00%), which is the
exact reverse of Z+jets (25.34% vs 60.75%).

## And the oracle flips sign

| | Z+jets | VBF |
|---|---:|---:|
| oracle over the three collections | 69.39% | 95.62% |
| matched-noise null | 72.00% | 93.67% |
| **real complementarity** | **-2.61** | **+1.95** |

**On VBF the separate collections carry genuinely complementary information**
(~2 points of real, non-artifact headroom); on Z+jets they carry none. Same
mechanism as everywhere else in this document: the jet partition is informative
where the jets hold the hard scatter (VBF, 47% of HS timing weight) and
uninformative where they do not (Z+jets, 10%).

## The confidence flag works on both, with different reach

| | Z+jets | VBF |
|---|---:|---:|
| events with `t_in` defined | 53.8% | **93.3%** |
| agree (`< 60 ps`), of defined | 40.1% | 61.2% |
| core \| agree | 79.42% | **97.90%** |
| core \| disagree | 58.36% | 85.49% |
| **gap** | **+21.05** | **+12.40** |
| core \| no in-jet track at all | 60.51% | **69.17%** |

Larger absolute gap on Z+jets, but VBF's baseline is 28 points higher so there
is less room. Not multiplicity in disguise on either sample -- on VBF the gap
runs +13.1 to +21.7 across track-multiplicity bins.

Note the last row: on VBF, having **no in-jet track at all** is itself a strong
negative indicator (69.17% against a 91.49% average) and covers only 6.7% of
events. On Z+jets that condition covers 46.2% of events and barely
discriminates (60.51% vs 63.89%), so it is only usable in VBF-like topologies.
