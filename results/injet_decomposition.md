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

## Event-level n forward HS tracks — the distribution

Sum of the per-cluster `truth_n_hs_tracks` over each event.

| n HS | VBF n | VBF frac | VBF cum | zjets n | zjets frac | zjets cum |
|---:|---:|---:|---:|---:|---:|---:|
| 0 | 343 | 0.86% | 0.9% | 8844 | **9.41%** | 9.4% |
| 1 | 695 | 1.75% | 2.6% | 9501 | 10.11% | 19.5% |
| 2 | 1223 | 3.08% | 5.7% | 9568 | 10.18% | 29.7% |
| 3 | 1775 | 4.47% | 10.2% | 9430 | 10.03% | 39.7% |
| 4 | 2378 | 5.99% | 16.2% | 9008 | 9.59% | 49.3% |
| 5 | 2920 | 7.36% | 23.5% | 8140 | 8.66% | 58.0% |
| 6 | 3159 | 7.96% | 31.5% | 7507 | 7.99% | 66.0% |
| 7 | 3310 | 8.34% | 39.8% | 6350 | 6.76% | 72.7% |
| 8 | 3351 | 8.44% | 48.3% | 5559 | 5.92% | 78.6% |
| 9 | 3124 | 7.87% | 56.1% | 4608 | 4.90% | 83.5% |
| 10 | 2951 | 7.44% | 63.6% | 3686 | 3.92% | 87.5% |
| 11-15 | 9882 | 24.9% | 88.5% | 9303 | 9.9% | 97.4% |
| 16-19 | 3102 | 7.82% | 96.3% | 1877 | 2.00% | 99.4% |
| >=20 | 1468 | 3.70% | 100% | 596 | 0.63% | 100% |

| | mean | median | p25 | p75 | <=2 | <=5 |
|---|---:|---:|---:|---:|---:|---:|
| VBF | 9.35 | 9 | 6 | 12 | 5.70% | 23.52% |
| **zjets** | **5.40** | **5** | **2** | **8** | **29.70%** | **57.98%** |

**The shapes differ in kind, not just in mean.** VBF is peaked at 7-8. Z+jets is
almost FLAT from 0 to 4 — each of those bins holds 9.5-10.2% of the sample — then
falls. So Z+jets is not a shifted VBF; it has a large, structurally sparse
population that VBF simply does not have. **58% of Z+jets events have <= 5
forward HS tracks**, against 23.5% in VBF, and in that regime nothing exceeds
~70% core fraction.

## The reco proxy is too weak to gate on

A multiplicity-gated chooser needs a reconstructable stand-in for the truth
count. `n_fwd_tracks_reco` is not it:

| | mean n reco | median | **corr(n_reco, n_HS)** | zero-reco |
|---|---:|---:|---:|---:|
| VBF | 29.2 | 28 | **0.421** | 0.00% |
| zjets | 25.8 | 25 | **0.345** | 0.00% |

At correlation 0.35 on Z+jets — and with 25.8 reco forward tracks per event
against 5.4 genuine HS ones — the gate would be badly smeared, and no event ever
has zero reco tracks to key on. **This substantially weakens implication 3
above**: the truth-gated crossover is real but there is currently no reco
variable that locates it. Finding one is a prerequisite, not a detail.

---

# TimeSource implemented in the main program

`Score::timeSource` (`enum class TimeSource { RAW, IN_JET, OUT_JET }`) now
selects which subset of the CHOSEN cluster's tracks its reported time is built
from. `calculateTime` dispatches on the flag instead of the old hard-coded
`score == Score::WAVES || ...` list, so the re-timing is a knob any score can
set rather than a property welded to one. OUT_JET is the exact complement of the
same predicate. All three fall back to the full-cluster time on an empty subset,
matching the previous WAVES behaviour.

WAVES / WAVES_MISCL / WAVES_MISAS / VBF_R1 / VBF_R2 carry `TimeSource::IN_JET`
so the refactor is behaviour-preserving. **Regression check: WAVES reproduces at
94.11% exactly, unchanged from the pre-refactor run.**

`TRKPTZ_TZ_IJ` (28) and `TRKPTZ_TZ_OJ` (29) copy TRKPTZ_TZ's score value
verbatim, so selection is bit-identical across the three and any separation is
the TIMING alone. `makeComparisonPlots` emits `timesource_<key>.pdf`.

## Local VBF (clustering_hist, m_jj >= 200, 45,335 events)

| | core frac | vs TRKPTZ |
|---|---:|---:|
| **TRKPTZ_TZ, in-jet t** | **94.35%** | **+1.32** |
| WAVES (deployed: WAVeS selection + in-jet t) | 94.11% | +1.08 |
| TRKPTZ_TZ, raw t | 93.68% | +0.64 |
| TRKPTZ | 93.04% | — |
| TRKPTZ_TZ, out-jet t | 91.60% | -1.44 |

**TRKPTZ_TZ + in-jet timing beats deployed WAVeS**, using a selector that never
looks at a jet.

## Where the timing actually lives

| | HS tracks that are in-jet | HS purity in-jet | HS purity out-jet | enrichment |
|---|---:|---:|---:|---:|
| local VBF | **47.0%** | 76.7% | 21.1% | 3.6x |
| Z+jets | **10.2%** | 42.4% | 19.8% | 2.1x |

TRKPTZ_TZ selection, three timings (export, m_jj >= 500):

| | VBF | Z+jets |
|---|---:|---:|
| raw (full cluster) | 91.49% | **63.89%** |
| in-jet | **92.31%** (+0.82) | 62.46% (-1.43) |
| out-jet | 89.31% (-2.18) | 63.61% (-0.27) |

**The mechanism is statistics, not purity.** Z+jets forward jets are still
HS-*enriched* (42.4% vs 19.8%, 2.1x) -- the earlier shorthand that they are
"usually pileup" was too strong. The problem is that they contain only **10.2%
of the hard-scatter tracks**, so the filter trades away 90% of the signal for a
2x purity gain. In VBF the same filter keeps 47% of the signal at 3.6x
enrichment, which is why the identical operation is +0.82 there and -1.43 here.

Corollary: on Z+jets the out-of-jet time is nearly identical to the raw time
(-0.27), because the out-of-jet subset IS ~95% of the tracks and ~90% of the
hard scatter. There is no useful third option there -- only in VBF is the
partition informative.

# Multiplicity-gated chooser: measured, and it does not work

Z+jets, TRKPTZ_TZ cluster answer vs classical aggregate t0 (union 72.31%,
so **+8.43 is available**):

| gate | best | vs cluster-only |
|---|---:|---:|
| **truth** n_HS <= 2 *(unreachable upper bound)* | 64.28% | **+0.40** |
| reco n_fwd_tracks <= 22 | 64.28% | **+0.40** |
| coarse 2-bin reco gate, out-of-fold | 64.26% | +0.37 |

**The truth-gated version is no better than the reco-gated one.** So the weak
`corr(n_reco, n_HS) = 0.345` is not the binding constraint -- multiplicity
simply is not the variable that says which answer to trust. It captures 4.7% of
the available headroom, and a perfect multiplicity oracle would capture no more.

This closes the idea. It also confirms the ratio it was premised on: genuine HS
tracks are **20.9%** of reconstructed forward tracks on Z+jets (VBF ~32%), i.e.
about 1 in 5.

---

# Time-source oracle (DIAGNOSTIC ONLY) — Z+jets

For the cluster TRKPTZ_TZ chose, compute all three times and report whichever is
closest to truth. Selection never changes, so this isolates the value of picking
the TIME SOURCE per event. Truth-using; not a method.

**With a matched-noise null**, because best-of-N is upward-biased even on
uninformative candidates: the null replaces the two alternatives with Gaussian
perturbations of the raw time at the observed raw->in-jet (18.6 ps) and
raw->out-jet (4.7 ps) spreads, then takes the same best-of-3.

| time source | core frac | core RMS | 68% half | median \|dt\| | full RMS |
|---|---:|---:|---:|---:|---:|
| raw | 63.89% | 23.40 | 84.69 | 29.05 | 182.4 |
| in-jet | 62.46% | 24.37 | 88.81 | 32.25 | 183.2 |
| out-of-jet | 63.61% | 23.93 | 85.35 | 30.33 | 182.6 |
| ORACLE2 (raw+in-jet) | 64.88% | 22.01 | 80.52 | 25.69 | 180.9 |
| **ORACLE3** | **65.01%** | 21.75 | 79.82 | 25.20 | 180.6 |
| **NULL best-of-3** | **65.55%** | 20.55 | 76.78 | 22.58 | 178.1 |

| | |
|---|---:|
| oracle3 gain over raw | +1.12 |
| null gain over raw | **+1.67** |
| **real complementarity** | **-0.54** |

**The null BEATS the real oracle.** The three time sources carry less
complementary information than three independent Gaussians of the same width --
they are highly correlated (raw and out-of-jet differ by only 4.7 ps RMS,
because out-of-jet IS ~90% of the tracks). The apparent +1.12 is entirely a
best-of-N artifact and then some.

Purity of the tracks used (pT-weighted HS fraction). The SELECTED-CLUSTER purity
is identical across rows by construction; only the subset differs:

| subset used for the time | purity |
|---|---:|
| whole cluster (= raw) | 42.24% |
| in-jet subset | **53.23%** |
| out-of-jet subset | 40.00% |

Note the in-jet subset IS purer (53.2 vs 42.2) and still times worse -- the
statistics argument again, not a purity one.

Oracle wins: raw 57.5%, in-jet 28.0%, out-of-jet 14.4%.

Differential vs event-level n forward HS tracks -- "real" is the oracle minus
the null, i.e. the part that is not best-of-N inflation:

| n HS | events | raw | in-jet | out-jet | ORACLE3 | NULL | **real** |
|---|---:|---:|---:|---:|---:|---:|---:|
| 0 | 8844 | 18.69 | 18.52 | 18.70 | 19.87 | 21.29 | **-1.42** |
| 1 | 9501 | 34.91 | 34.82 | 34.79 | 36.41 | 37.80 | -1.39 |
| 2 | 9568 | 45.55 | 45.22 | 45.21 | 47.18 | 48.00 | -0.83 |
| 3 | 9430 | 55.30 | 54.06 | 54.83 | 56.72 | 57.34 | -0.62 |
| 4 | 9008 | 64.29 | 62.70 | 63.77 | 65.48 | 66.03 | -0.56 |
| 5 | 8140 | 70.90 | 69.47 | 70.47 | 72.11 | 72.49 | -0.38 |
| 6-7 | 13857 | 78.31 | 76.03 | 78.01 | 79.38 | 79.58 | -0.20 |
| 8-9 | 10167 | 85.46 | 83.09 | 85.26 | 86.32 | 86.38 | -0.06 |
| 10-13 | 10706 | 91.55 | 89.30 | 91.36 | 92.18 | 92.19 | -0.01 |
| 14-19 | 4160 | 96.27 | 93.89 | 96.15 | 96.54 | 96.49 | +0.05 |

Negative in every bin with statistics. **There is nothing to gain from
per-event time-source selection on Z+jets** -- do not build a chooser for it.

## VBF control — where the same diagnostic IS positive

| | raw | in-jet | out-jet | ORACLE3 | NULL | **real** |
|---|---:|---:|---:|---:|---:|---:|
| core frac | 91.49% | 92.31% | 89.31% | 92.92% | 92.30% | **+0.62** |
| core RMS | 17.72 | 16.42 | 22.10 | 12.48 | 14.36 | |

Oracle wins: raw 28.4%, **in-jet 49.4%**, out-of-jet 22.2%; in-jet subset purity
88.25% against 72.66% for the whole cluster.

So the diagnostic is not vacuous -- it finds real complementarity on VBF
(+0.62) and none on Z+jets (-0.54). Consistent with everything else in this
document: the jet partition is informative where the jets hold the hard scatter
(VBF, 47% of HS timing weight) and uninformative where they do not (Z+jets,
10%).
