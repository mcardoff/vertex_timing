# WAVeS decomposed: selection vs in-jet re-timing, and the aggregate t0 differential

`python/injet_decomposition.py`. The deployed `Score::WAVES` row does two
independent things -- ranks by the WAVeS score, AND recomputes the chosen
cluster's time from only its tracks within dR < 0.4 of a qualifying forward jet.
Every published WAVeS number folds the two together. The in-jet predicate is
reproduced exactly: `dr_nearest_fwdjet` comes from `nearestForwardJet`, which
applies the same jet qualification as `calculateTime`'s WAVES branch, and the
no-surviving-track fallback to the full-cluster time matches too.

## The re-timing is orthogonal to the selection, and worth more than it

**Local VBF** (39,681 events, 19.6% of tracks in-jet):

| selection | standard time | in-jet time | re-timing |
|---|---:|---:|---:|
| TRKPTZ | 90.83% | 91.66% | **+0.83** |
| TRKPTZ_TZ | 91.49% | **92.31%** | **+0.82** |
| WAVeS | 91.35% | 92.18% | **+0.82** |

**The re-timing gain is identical (+0.82-0.83) for all three selectors**,
including two that never look at a jet. It is a pure additive term, completely
separable from ranking.

So WAVeS decomposes as **selection +0.52, re-timing +0.82** -- the re-timing is
worth ~1.6x the selection, and every WAVeS number quoted to date has been
crediting the ranking for it.

**The best combination is TRKPTZ_TZ + in-jet timing at 92.31%, which beats
deployed WAVeS (92.18%)** -- i.e. a z-only selector wearing WAVeS's timing trick
beats WAVeS itself.

## But the re-timing flips sign on Z+jets

**Z+jets** (93,977 events, **5.0%** of tracks in-jet):

| selection | standard time | in-jet time | re-timing |
|---|---:|---:|---:|
| TRKPTZ | 62.12% | 60.75% | **-1.37** |
| TRKPTZ_TZ | **63.89%** | 62.46% | **-1.43** |
| WAVeS | 60.21% | 58.81% | **-1.40** |

Again identical across selectors (-1.37/-1.43/-1.40), and again a pure additive
term -- but negative. Two documented reasons, both pointing the same way: only
5.0% of Z+jets tracks are in-jet against 19.6% in VBF, so the filter discards
95% of the timing information; and Z+jets forward jets are usually PILEUP, so
what survives the filter is preferentially pileup-proximate.

**This is a sample-dependent sign flip, not a magnitude difference.** Any
deployment of the in-jet re-timing needs a gate, and the in-jet track fraction
is the obvious reco-only candidate for it.

Restricting the AGGREGATE to in-jet tracks is worse still: -0.98 on VBF,
**-9.79** on Z+jets (62.89 -> 53.10). The aggregate's whole advantage is
statistics, and the jet filter removes 80-95% of it.

## Aggregate t0, differential

> **Corrected 2026-09-01.** The first version of this section binned on
> `truth_n_hs_tracks` via `.first()` per event. That column is **per cluster**
> ("HS tracks in this cluster"), so `.first()` took whichever cluster happened
> to be listed first, not the event total. It inflated the zero-HS-track
> population from 9.41% to 35.45% on Z+jets (0.86% -> 5.71% on VBF) and
> exaggerated every entry. The event-level count is the SUM over the event's
> clusters; tables below use that.

**Local VBF**, vs n truth forward HS tracks (event-level):

| n HS trk | events | TRKPTZ | TZ | WAVeS | agg t0 | agg - TZ |
|---|---:|---:|---:|---:|---:|---:|
| 0 | 343 | 22.16 | 21.28 | 21.87 | 25.66 | **+4.37** |
| 1 | 695 | 44.17 | 46.91 | 48.06 | 49.78 | +2.88 |
| 2 | 1223 | 60.34 | 62.47 | 64.27 | 64.02 | +1.55 |
| 3 | 1775 | 73.41 | 74.42 | 75.77 | 76.51 | +2.08 |
| 4 | 2378 | 80.82 | 83.10 | 83.43 | 83.47 | +0.38 |
| 5 | 2920 | 86.44 | 87.60 | 88.25 | 87.84 | +0.24 |
| 6-7 | 6469 | 90.99 | 92.18 | 92.21 | 92.55 | +0.37 |
| 8-9 | 6475 | 95.18 | 95.61 | 95.41 | 95.29 | -0.32 |
| 10-13 | 9882 | 97.54 | 97.65 | 96.86 | 98.04 | +0.38 |
| >=14 | 7521 | 99.5 | 99.5 | 98.8 | 99.3 | ~-0.1 |

**Z+jets**:

| n HS trk | events | TRKPTZ | TZ | WAVeS | agg t0 | agg - TZ |
|---|---:|---:|---:|---:|---:|---:|
| 0 | 8844 | 19.69 | 18.69 | 19.41 | 21.30 | **+2.61** |
| 1 | 9501 | 33.38 | 34.91 | 32.86 | 35.89 | +0.98 |
| 2 | 9568 | 44.00 | 45.55 | 43.29 | 46.06 | +0.51 |
| 3 | 9430 | 53.63 | 55.30 | 52.29 | 55.28 | -0.02 |
| 4 | 9008 | 61.49 | 64.29 | 59.45 | 62.71 | -1.58 |
| 5 | 8140 | 67.99 | 70.90 | 66.25 | 68.92 | -1.98 |
| 6-7 | 13857 | 75.77 | 78.31 | 72.89 | 75.23 | **-3.08** |
| 8-9 | 10167 | 82.99 | 85.46 | 79.64 | 82.88 | -2.59 |
| 10-13 | 10706 | 90.06 | 91.55 | 86.78 | 89.34 | -2.20 |
| 14-19 | 4160 | 95.22 | 96.27 | 92.69 | 94.30 | -1.97 |

**The aggregate wins at LOW multiplicity and loses at high, on both samples**,
with the crossover at ~3 HS tracks on Z+jets and ~8 on VBF. It wins the 0-track
bin on both (+4.37 / +2.61) -- the opposite of what the buggy version showed.
Magnitudes are smaller than first reported (peak +4.4 rather than +5.7).

**TRKPTZ_TZ is the best selector in almost every Z+jets bin**, and its margin
over the aggregate grows with multiplicity, which is what makes the crossover
usable as a gate.

## Z+jets is not merely sparse -- it is worse at matched multiplicity

| | events | median n HS trk | zero-HS | TRKPTZ |
|---|---:|---:|---:|---:|
| local VBF | 39,681 | 9 | 0.86% | 90.83% |
| Z+jets | 93,977 | 5 | 9.41% | 62.12% |

Reweighting Z+jets' per-multiplicity rates onto VBF's multiplicity spectrum:

| | core fraction |
|---|---:|
| Z+jets, own spectrum | 62.12% |
| **Z+jets rates, VBF spectrum** | **79.86%** |
| VBF, actual | 90.83% |

So of the 28.7-point gap, **17.7 points (62%) is sparsity** and **11.0 points
(38%) is Z+jets being genuinely harder at the same forward-HS-track count**.
At 6-7 HS tracks VBF reaches 92.18% where Z+jets manages 78.31%. Sparsity is the
larger term but it is not the whole story, and the residual 11 points is not
addressed by anything studied this session.

## Implications

1. **The in-jet re-timing is the single most transferable thing WAVeS has**, and
   it belongs in `calculateTime` behind a sample-appropriate gate rather than
   welded to one score. It is +0.82 on VBF for ANY selector, and -1.4 on Z+jets
   for any selector.
2. **TRKPTZ_TZ is the most robust selector tested**: best or near-best in almost
   every Z+jets bin, and +0.66 on VBF.
3. **A multiplicity-gated chooser between TRKPTZ_TZ and the aggregate** is
   motivated by a real crossover, but the win is small (a few points in the
   sparsest bins, which are a minority of events). The binning variable is truth
   (`truth_n_hs_tracks`); `n_fwd_tracks_reco` is the reco proxy and must be
   validated before this is actionable.
4. **Z+jets' residual 11-point deficit at matched multiplicity is unexplained**
   and is a better target than further selector work.
