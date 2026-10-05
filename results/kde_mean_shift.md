# TZP_KDE: kernel-density selection, mean-shift time, (t,z) clustering (2026-10-03)

> **Update 2026-10-05.** Two time-reliability factors were added to the track
> weights (`selection_reliability.md`): fails removed are now 19.8 / 11.1 /
> 16.5 / 18.8% (vbf / zjets / dijet / ttbar, grid C++). The numbers below are
> the 2026-10-03 state without them.

`Score::TZP_KDE` (id 34, main collection) and `Score::TZP_KDE_TZ` (id 35, its
own (t,z) collection). `applyKernelDensityScore` and `ITERATIVE_ZSEED` in
`src/clustering_functions.h`. Classical, closed-form, one configuration for
every sample, and **no input from Athena's `RecoVtx_time`**.

Study scripts, `python/kde_study/`:
- `wb.py`, `exp1.py`: load an export, reproduce TZP.
- `exp2..13.py`: the scoring and timing rounds (some carry the Athena
  agreement term, `hg`, that was later removed).
- `clus.cpp`, `reclus.py`: re-cluster the exported tracks with any algorithm
  and rebuild every cluster quantity from the tracks. Run with the production
  settings it reproduces the exported partition label for label in 100.00% of
  events on all four samples, and TZP to 0.00-0.02.
- `clus1..7.py`: the clustering rounds.
- `hist_core.py`: inclusive core fractions from a `clustering_hist` file.

Goal: remove >= 10% of the events TZP fails (|t0 - t_truth| >= 60 ps) without
losing what TZP already gets right.

## The algorithm

    clustering : iterative, jointly in (t, z0): absorb the nearest track while
                 sqrt(dt^2/(s_tA^2+s_tB^2) + dz^2/(s_zA^2+s_zB^2)) < 3
                 seeds taken in order of pT e^{-|z0 - z_PV|}          [TZP_KDE_TZ]
                 (TZP_KDE keeps the production time-only clustering)
    w_t        = sqrt(pT) * exp(-1.0 |z0_t - z_PV|)                   per track
    t_C        : seed = SUM_C (w/s_t^2) t / SUM_C (w/s_t^2), then 3 passes of
                 t <- SUM_all (w/s_t^2) K60(t_i - t) t_i / SUM_all (w/s_t^2) K60(t_i - t)
    S(C)       = [SUM_all w K60(t_i - t_C)]
                 * exp(-0.6 |z_C - z_PV|) * (SUM_C 1/var_d0)^0.225    (TZP's envelope, hard cluster)
    t0         = t_C of the argmax; if the event has ONE cluster, TZP's time

K60 is a Gaussian of width 60 ps over the times of ALL clustered tracks.

## Result

Fraction of TZP's failing events removed, training exports (novbs):

| sample | events | TZP | TZP_KDE | TZP_KDE_TZ | fails removed (KDE / KDE_TZ) |
|---|---:|---:|---:|---:|---:|
| zjets (all) | 93,977 | 64.71 | 67.54 | **67.96** | 8.0% / **9.2%** |
| dijet (all) | 132,138 | 88.97 | 90.35 | **90.54** | 12.6% / **14.3%** |
| vbf (files 0-119) | 58,162 | 91.67 | 92.84 | **92.93** | 14.1% / **15.1%** |
| ttbar (files 0-59) | 61,129 | 88.93 | 90.44 | **90.70** | 13.6% / **16.0%** |

C++ (`clustering_hist`, local VBF, m_jj >= 200, 48,323 events):

| | pass | fail | core | fails removed |
|---|---:|---:|---:|---:|
| HGTD (Athena) | 41,007 | 7,316 | 84.86 | |
| TRKPTZ | 43,656 | 4,667 | 90.34 | |
| WAVeS | 44,311 | 4,012 | 91.70 | |
| TZP | 44,446 | 3,877 | 91.98 | |
| TZP_KDE | 44,992 | 3,331 | 93.11 | 14.1% |
| **TZP_KDE_TZ** | **45,059** | **3,264** | **93.25** | **15.8%** |

The C++ matches the offline prediction (14.1% and 15.1-17.2%), and the
pre-existing rows are unchanged.

**Z+jets misses the 10% bar (9.2%; 9.8% on the canonical mjj500 selection).**
21.5% of its TZP fails are events with no timed hard-scatter track, which no
selector can recover; of the rest, 11.7% are removed.

## Out of sample

Constants were chosen on the first four rows. Nothing re-tuned below.

| set | events | TZP | TZP_KDE_TZ | fails removed |
|---|---:|---:|---:|---:|
| vbf novbs, files 120-299 | 87,187 | 91.65 | 93.09 | 17.2% |
| ttbar novbs, files 60-179 | 130,330 | 89.11 | 90.85 | 16.0% |
| zjets mjj500 | 36,531 | 64.31 | 67.81 | 9.8% |
| dijet mjj500 | 80,383 | 88.58 | 90.19 | 14.0% |
| vbf mjj500, files 380-499 | 42,548 | 91.96 | 93.28 | 16.3% |
| ttbar mjj500, files 300-399 | 69,305 | 88.60 | 90.40 | 15.8% |

The zjets and dijet mjj500 rows overlap their novbs tuning sets (a subset under
a stricter selection); the vbf and ttbar rows are disjoint files.

mu = 0 controls, passing events:

| set | events | TZP | TZP_KDE | TZP_KDE_TZ |
|---|---:|---:|---:|---:|
| vbf_mu0, files 0-19 | 57,670 | 57,592 | 57,611 | 57,606 |
| ttbar_mu0, files 0-19 | 42,756 | 42,692 | 42,703 | 42,697 |
| zeejets_mu0 | 34,310 | 34,227 | 34,237 | 34,228 |

No harm. The single-cluster guard matters here: without it the weighted time
loses 10-23 events per sample, since ~90% of mu = 0 events have one cluster and
nothing to purify. At mu = 200 the guard touches under 1% of events.

## Performance is maintained

Core fraction by number of timed HS tracks, TZP -> TZP_KDE_TZ:

| n HS timed | zjets | dijet | vbf | ttbar |
|---|---|---|---|---|
| 0 | 19.5 -> 18.9 | 20.4 -> 20.7 | 23.9 -> 22.2 | 20.1 -> 19.0 |
| 1 | 34.7 -> 37.4 | 37.9 -> 40.8 | 43.6 -> 47.3 | 38.7 -> 43.1 |
| 2-3 | 51.5 -> 55.6 | 60.8 -> 64.6 | 71.3 -> 74.6 | 62.5 -> 66.7 |
| 4-7 | 73.3 -> 77.5 | 85.4 -> 88.0 | 90.3 -> 92.1 | 85.3 -> 87.9 |
| >= 8 | 90.8 -> 93.8 | 97.4 -> 98.3 | 98.1 -> 98.8 | 97.1 -> 98.2 |

The 0-track row is chance (there is no right answer in the event), and it
drifts by about a point either way.

Fails removed at other windows (30 / 60 / 90 / 150 ps): zjets 9.7 / 9.2 / 7.4 /
5.4%, dijet 15.8 / 14.3 / 11.2 / 7.1%, vbf 14.8 / 15.1 / 10.8 / 4.0%, ttbar
18.2 / 16.0 / 12.3 / 6.9%. RMS inside the 60 ps core: zjets 23.2 -> 21.8 ps,
dijet 16.7 -> 15.6, vbf 16.3 -> 15.6, ttbar 17.5 -> 16.0.

## Where the gain comes from (percent of TZP fails removed)

| step | zjets | dijet | vbf | ttbar |
|---|---:|---:|---:|---:|
| TZP selection, time re-weighted by pT e^-dz / s_t^2 inside the cluster | 3.5 | 5.4 | 8.6 | 7.2 |
| + mean-shift over all event tracks (60 ps) | 5.6 | 9.1 | 13.1 | 11.3 |
| + select on the kernel density at the mode | 6.5 | 10.6 | 13.6 | 12.3 |
| + sqrt(pT) in both weights (= TZP_KDE) | 8.0 | 12.6 | 14.1 | 13.6 |
| + (t,z) clustering, z-weighted seeds (= TZP_KDE_TZ) | 9.2 | 14.3 | 15.1 | 16.0 |

**The time estimator is the largest single piece.** With the selection
untouched, the inverse-variance cluster mean leaves 6-13% of the fails on the
table; an oracle timing the selected cluster from its truth-HS tracks removes
13-20%.

Plateaus: kernel width 55-75 ps, 3 vs 6 passes identical, pT power 0.35-0.65 or
a 5 GeV cap. Re-scanning the envelope's beta and gamma on the (t,z) clusters
returns (0.6, 0.45) again.

## The clustering scan

Each row re-clusters the same tracks and scores them two ways: TZP as is (hard
cluster sums, guarded in-jet time), and the kernel-density score. Numbers are
percent of the production-TZP fails removed, zjets / dijet / vbf / ttbar.

| clustering | TZP score | kernel-density score |
|---|---|---|
| iterative 3.0 sigma (production) | 0 / 0 / 0 / 0 | 8.0 / 12.6 / 14.1 / 13.6 |
| iterative 2.0 | +1.2 / +0.2 / -1.4 / +1.6 | 8.3 / 12.9 / 15.1 / 15.0 |
| iterative 2.5 | +1.1 / +1.0 / +0.4 / +2.3 | 8.4 / 13.0 / 14.7 / 14.4 |
| iterative 4.0 | -5.2 / -8.4 / -8.3 / -10.6 | 7.4 / 11.3 / 13.1 / 13.0 |
| simultaneous (agglomerative) 3.0 | -0.8 / -9.6 / -14.2 / -7.0 | 8.5 / 13.3 / 15.4 / 15.0 |
| cone 2.0 | +0.3 / +1.2 / +4.6 / +2.7 | 8.1 / 12.5 / 15.0 / 14.7 |
| cone 3.0 | -4.4 / -3.7 / +2.8 / -3.3 | 7.5 / 11.6 / 13.8 / 13.2 |
| iterative 3.0, seeds by pT e^-dz | +1.2 / +0.9 / +0.6 / +1.0 | 8.1 / 12.7 / 13.8 / 13.9 |
| iterative 3.0 + 30 ps floor in the distance | -8.0 / -14.4 / -15.2 / -18.9 | 6.7 / 10.4 / 11.4 / 11.9 |
| fixed 60 ps window, weighted centroid | +3.1 / +4.8 / +6.1 / +5.9 | 8.6 / 13.0 / 15.5 / 14.6 |
| fixed 90 ps window, weighted centroid | +3.4 / +7.6 / +10.8 / +8.0 | 7.9 / 12.5 / 14.6 / 14.5 |
| mean-shift modes, W = 20 ps | -1.1 / -6.9 / -11.6 / -6.0 | 8.3 / 13.0 / 14.4 / 14.0 |
| mean-shift modes, W = 45 ps | -11.0 / -25.4 / -32.3 / -33.4 | 6.8 / 10.2 / 11.1 / 11.9 |
| iterative (t,z) 2.5 | +2.5 / +2.8 / +3.1 / +4.0 | 9.1 / 14.1 / 15.4 / 15.7 |
| iterative (t,z) 3.0 | +2.6 / +4.2 / +4.2 / +4.9 | 9.2 / 14.3 / 15.1 / 15.5 |
| **iterative (t,z) 3.0, seeds by pT e^-dz** | **+3.7 / +5.6 / +5.0 / +6.3** | **9.2 / 14.3 / 15.1 / 16.0** |
| iterative (t,z) 4.0 | -1.0 / -1.1 / -2.0 / -2.3 | 8.1 / 12.2 / 13.7 / 13.7 |
| (t,z) 3.0, sigma_z x 0.5 | +0.1 / -2.0 / -1.9 / -1.5 | 8.0 / 12.0 / 14.1 / 13.8 |
| (t,z) 3.0, sigma_z x 1.5 | +1.5 / +3.2 / +3.9 / +3.6 | 8.6 / 13.4 / 14.4 / 14.5 |
| simultaneous (t,z) 3.0 | -0.3 / -9.5 / -11.3 / -9.1 | 8.9 / 13.2 / 14.8 / 14.7 |

Three things this shows.

- **The kernel-density score barely cares which partition supplies its
  candidates** (zjets 6.7-9.2 across every algorithm tried). It sums over all
  tracks, so the partition only fixes the seeds and the cluster envelope.
- **For the hard-cluster TZP score the clustering matters a lot**, mostly
  downward. The two changes that help it are a weighted centroid with a fixed
  time window (up to +10.8 on VBF) and clustering jointly in z (+4 to +6).
- **Joint (t,z) clustering is the one change that helps both.** Clusters become
  z-coherent, so the envelope's z centroid and d0 sum describe one vertex
  instead of a time coincidence. A z-aware kernel on top (down-weighting tracks
  z-incompatible with the cluster) costs 0.2-3 points at every width.

## The clustering is now within 2 points of its ceiling

Truth-assisted partitions, scored and timed exactly as above:

| partition | zjets | dijet | vbf | ttbar |
|---|---:|---:|---:|---:|
| reco (t,z) clusters (TZP_KDE_TZ) | 9.2 | 14.3 | 15.1 | 16.0 |
| in-time truth-HS tracks forced into one pure cluster | 10.7 | 16.1 | 17.7 | 18.3 |
| same partition, TZP score instead | -1.8 | 0.0 | +2.3 | +0.6 |
| same partition + ORACLE choice of that cluster | 63.3 | 77.0 | 80.6 | 81.9 |

A perfect partition is worth 1.5-2.6 points more than the reco one. Choosing
the pure HS cluster once it exists is worth 63-82. **What is left is selection,
not clustering.** Forcing ALL truth-HS tracks into one cluster (mis-timed ones
included) is worse than the reco partition (7.7 / 9.3 / 11.4 / 9.8).

Why the selection loses, Z+jets, events where the pure HS cluster exists and
passes but another cluster is picked (19.7% of events; VBF 5.2%):

| median | HS cluster (lost) | picked cluster | HS cluster when it wins |
|---|---:|---:|---:|
| n tracks | 2 | 8 | 5 |
| sum pT [GeV] | 3.3 | 12.4 | 10.6 |
| \|z_C - z_PV\| [mm] | 0.65 | 0.45 | 0.19 |
| sigma of z centroid [mm] | 0.66 | 0.33 | 0.28 |
| z chi2/ndf about its centroid | 0.42 | 1.54 | 1.09 |
| t chi2/ndf | 0.37 | 1.46 | 1.04 |

The lost hard-scatter cluster is two or three soft tracks against an
eight-track pileup coincidence that is larger, harder and closer in z. The only
variables pointing the right way are the internal chi2s, and as score terms
those were worth under 1% (they are loss-subset rankings; see the methodology
note in CLAUDE.MD).

## Closed negatives from this study

- Expected-pileup subtraction under the cluster and a tail boost `S / G(t)^k`:
  strongly negative at every strength.
- Per-track z0 precision in the selection weight (`sigma_z0^-b`, a z0
  likelihood ratio): negative or neutral, although HS purity rises from 24% at
  sigma_z0 > 1.6 mm to 98% below 0.1 mm.
- A kernel-weighted cluster envelope instead of the hard cluster's: costs 4-7
  points of the gain.
- A z-aware kernel, in the density or in the mean-shift.
- In-jet boosts of the time weight, a sigma_t-scaled kernel, a box kernel,
  1.5-2 sigma trimming with inverse-variance weights.
- chi2/ndf, multiplicity and sigma_t powers on the score: within +-1%.
- Distance floors, wider cuts, agglomerative and mean-shift-mode partitions
  (table above).

## Removed on request: the Athena agreement term

An earlier version multiplied the score by
`1 + 0.5 exp(-((t_C - t_HGTD)/30 ps)^2/2)` when `RecoVtx_time` was valid. It was
worth 0.5 (zjets) to 4 (vbf) more points of the fails. It is gone: the score is
meant to replace Athena's time, not depend on it.

## Grid run of the C++ (clusters 3248506-09, default selection m_jj >= 200)

46 shards, none held; merged with `hist_merge`. Inclusive core fraction:

| | vbf | zjets | dijet | ttbar |
|---|---:|---:|---:|---:|
| events | 630,720 | 64,950 | 105,137 | 863,552 |
| HGTD (Athena) | 84.53 | 41.60 | 77.88 | 77.24 |
| TRKPTZ | 90.02 | 61.94 | 86.66 | 86.71 |
| WAVeS | 91.49 | 58.50 | 85.98 | 85.98 |
| TZP | 91.69 | 64.48 | 88.30 | 88.50 |
| TZP_KDE | 92.94 | 67.35 | 89.74 | 90.14 |
| **TZP_KDE_TZ** | **93.03** | **67.78** | **89.94** | **90.29** |
| TZP fails removed by TZP_KDE_TZ | 16.2% | 9.3% | 14.0% | 15.6% |

Differential plots: `clustering_plot --sample=<s> --hist-file=<merged>` writes
`<s>/comparisons/<s>_newmethod_<key>.pdf` (HGTD / TRKPTZ / TZP / TZP_KDE_TZ)
and `<s>_kde_<key>.pdf` (TRKPTZ / WAVeS / TZP / TZP_KDE / TZP_KDE_TZ) for the
five keys. `python/kde_study/grid_core.py` prints the table above.

## Not done

- `rpt_v5_hist`, `export_training_data` and the event display still use TZP.
