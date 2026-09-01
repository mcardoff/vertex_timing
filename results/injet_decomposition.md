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

**Local VBF**, vs n truth forward HS tracks:

| n HS trk | events | TRKPTZ | TZ | WAVeS | agg t0 | agg - TZ |
|---|---:|---:|---:|---:|---:|---:|
| 0 | 2265 | 47.64 | 51.26 | 52.49 | 44.50 | **-6.75** |
| 1 | 2256 | 59.09 | 58.78 | 57.18 | 64.01 | **+5.23** |
| 2 | 1725 | 72.81 | 74.26 | 72.70 | 79.94 | **+5.68** |
| 3 | 2157 | 82.71 | 84.05 | 84.47 | 87.81 | **+3.76** |
| 4 | 2789 | 91.07 | 92.40 | 92.51 | 94.05 | +1.65 |
| 5 | 3146 | 94.53 | 95.33 | 95.84 | 95.39 | +0.06 |
| 6-7 | 6803 | 97.31 | 98.02 | 97.88 | 97.50 | -0.51 |
| 8+ | 18540 | ~99 | ~99 | ~99 | ~99 | ~0 |

**Z+jets**:

| n HS trk | events | TRKPTZ | TZ | WAVeS | agg t0 | agg - TZ |
|---|---:|---:|---:|---:|---:|---:|
| 0 | **33315** | 36.31 | 38.32 | 35.40 | 35.37 | -2.95 |
| 1 | 13981 | 47.71 | 49.49 | 45.86 | 50.55 | +1.06 |
| 2 | 7549 | 60.88 | 62.96 | 58.93 | 66.83 | **+3.87** |
| 3 | 6316 | 74.48 | 76.63 | 71.53 | 77.91 | +1.28 |
| 4 | 5836 | 83.16 | 85.20 | 80.35 | 84.44 | -0.75 |
| 5 | 5464 | 88.21 | 90.32 | 84.99 | 88.49 | -1.83 |
| 6-7 | 8553 | 92.82 | 94.42 | 89.98 | 92.46 | -1.96 |
| 10-13 | 5461 | 98.54 | 98.96 | 96.32 | 97.73 | -1.23 |

**The aggregate wins in a low-but-NONZERO multiplicity band and loses outside
it, on both samples.** VBF: wins 1-4 (up to +5.7), loses at 0 (-6.75). Z+jets:
wins 1-3 (up to +3.87), loses at 0 and everywhere above 4.

This CORRECTS an earlier reading in `classical_scoring_brainstorm.md` that the
aggregate is "relatively stronger where tracks are plentiful". Within a sample
the opposite is true -- it is strongest in sparse events, provided there is at
least one hard-scatter track to find. Z+jets looks worse inclusively only
because **35% of its events sit in the 0-HS-track bin**, where nothing works
(every method 35-38%) and the aggregate is worst.

Against n forward jets the differential is nearly flat (VBF +0.10 to +0.64,
Z+jets -0.93 to -2.27), so jet count is not the axis that separates them --
consistent with every other jet-based handle tested this session.

## Implications

1. **The in-jet re-timing is the single most transferable thing WAVeS has**, and
   it belongs in `calculateTime` behind a sample-appropriate gate rather than
   welded to one score. It is +0.82 on VBF for ANY selector.
2. **TRKPTZ_TZ is the most robust selector tested**: best or near-best in almost
   every Z+jets bin, and +0.66 on VBF.
3. **A multiplicity-gated chooser between TRKPTZ_TZ and the aggregate is
   well-motivated by the differential** -- the crossover is sharp and sits at the
   same place on both samples. The binning variable here is truth
   (`truth_n_hs_tracks`); `n_fwd_tracks_reco` is the reco proxy and must be
   checked before this is actionable.
4. **The 0-HS-track bin is the real Z+jets problem**: 35% of events, every
   method at 35-38%. No selector or estimator studied this session moves it.
