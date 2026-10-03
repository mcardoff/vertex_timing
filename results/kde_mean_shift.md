# TZP_KDE: kernel-density selection + mean-shift time (2026-10-03)

`Score::TZP_KDE` (id 34), `applyKernelDensityScore` in `src/clustering_functions.h`.
Study scripts: `python/kde_study/` (`wb.py` loads an export, `exp1..13.py` are the
rounds below, `hist_core.py` reads inclusive core fractions from a
`clustering_hist` file). Classical, closed-form, one configuration for every
sample.

Goal: remove >= 10% of the events TZP fails (|t0 - t_truth| >= 60 ps) without
losing what TZP already gets right.

## The algorithm

The iterative 3 sigma clustering is unchanged and supplies the candidates. What
changes is how a candidate is timed and scored: both now use a Gaussian kernel
over ALL clustered tracks of the event rather than the cluster's own members.

    w_t   = sqrt(pT) * exp(-1.0 |z0_t - z_PV|)                 per track
    t_C   : seed = SUM_C (w/sigma_t^2) t / SUM_C (w/sigma_t^2), then 3 passes of
            t <- SUM_all (w/sigma_t^2) K60(t_i - t) t_i / SUM_all (w/sigma_t^2) K60(t_i - t)
    S(C)  = [SUM_all w K60(t_i - t_C)]
            * exp(-0.6 |z_C - z_PV|) * (SUM_C 1/var_d0)^0.225     (TZP's envelope, hard cluster)
            * [1 + 0.5 exp(-((t_C - t_HGTD)/30 ps)^2 / 2)]         (only if RecoVtx_time is valid)
    t0    = t_C of the argmax; if the event has ONE cluster, TZP's time instead

K60 is a Gaussian of width 60 ps.

## Result

Fraction of TZP's failing events removed. Offline = the `training_data` exports
(novbs), python reimplementation; TZP reproduces to 0.03.

| sample | events | TZP | TZP_KDE | fails removed |
|---|---:|---:|---:|---:|
| zjets (all) | 93,977 | 64.69 | **67.69** | **8.5%** |
| dijet (all) | 132,138 | 88.98 | **90.50** | **13.8%** |
| vbf (files 0-119) | 58,162 | 91.69 | **93.19** | **18.1%** |
| ttbar (files 0-59) | 61,129 | 88.95 | **90.62** | **15.1%** |

C++ (`clustering_hist`, local VBF, m_jj >= 200, 48,323 events):

| | pass | fail | core |
|---|---:|---:|---:|
| TRKPTZ | 43,656 | 4,667 | 90.34 |
| WAVeS | 44,311 | 4,012 | 91.70 |
| TZP | 44,446 | 3,877 | 91.98 |
| **TZP_KDE** | **45,164** | **3,159** | **93.46** |

That is 18.5% of TZP's fails removed, against 18.1% predicted offline. Every
other row is unchanged by the code change.

**Z+jets misses the 10% bar (8.5%).** 21.5% of its TZP fails are events with no
timed hard-scatter track at all, which no selector can recover; of the rest,
10.8% are removed. The other three samples clear the bar.

## Out of sample

Constants were chosen on the four rows above. Disjoint file ranges and the
m_jj >= 500 exports, nothing re-tuned:

| set | events | TZP | TZP_KDE | fails removed | without HGTD term |
|---|---:|---:|---:|---:|---:|
| vbf novbs, files 120-299 | 87,187 | 91.68 | 93.31 | 19.6% | 15.8% |
| ttbar novbs, files 60-179 | 130,330 | 89.12 | 90.89 | 16.2% | 14.4% |
| zjets mjj500 | 36,531 | 64.29 | 67.51 | 9.0% | 8.5% |
| dijet mjj500 | 80,383 | 88.59 | 90.13 | 13.5% | 12.3% |
| vbf mjj500, files 380-499 | 42,548 | 91.99 | 93.57 | 19.7% | 15.8% |
| ttbar mjj500, files 300-399 | 69,305 | 88.62 | 90.45 | 16.1% | 14.4% |

The zjets and dijet mjj500 rows overlap their novbs tuning sets (a subset under
a stricter selection); the vbf and ttbar rows are disjoint files.

mu = 0 controls (passing events, single-cluster guard on):

| set | events | TZP | TZP_KDE |
|---|---:|---:|---:|
| vbf_mu0, files 0-19 | 57,670 | 57,593 (99.87) | 57,611 (99.90) |
| ttbar_mu0, files 0-19 | 42,756 | 42,691 (99.85) | 42,703 (99.88) |
| zeejets_mu0 | 34,310 | 34,227 (99.76) | 34,237 (99.79) |

Without the guard the same three lose 12, 10 and 23 events: 90% of mu = 0 events
have one cluster, and there the weighted time is a slightly noisier estimator
than the inverse-variance mean with nothing to purify. At mu = 200 the guard
touches 0.6-1.0% of events and moves nothing.

## Performance is maintained everywhere it was measured

Z+jets, by number of timed HS tracks (others have the same shape):

| n HS timed | share of events | TZP | TZP_KDE |
|---|---:|---:|---:|
| 0 | 9.4% | 19.5 | 19.6 |
| 1 | 10.1% | 34.7 | 37.2 |
| 2-3 | 20.2% | 51.5 | 54.8 |
| 4-7 | 33.0% | 73.3 | 77.2 |
| >= 8 | 27.3% | 90.7 | 93.6 |

Not specific to the 60 ps window. Fails removed at other windows:

| window | zjets | dijet | vbf | ttbar |
|---|---:|---:|---:|---:|
| 20 ps | 7.1% | 11.5% | 9.5% | 13.2% |
| 30 ps | 9.2% | 15.4% | 15.5% | 17.1% |
| 45 ps | 9.2% | 15.2% | 18.4% | 17.6% |
| 60 ps | 8.5% | 13.8% | 18.1% | 15.1% |
| 90 ps | 6.9% | 10.9% | 14.2% | 11.8% |
| 150 ps | 5.3% | 7.6% | 9.6% | 6.8% |

RMS inside the 60 ps core also narrows: zjets 23.2 -> 21.8 ps, dijet 16.7 -> 15.6,
vbf 16.2 -> 15.6, ttbar 17.4 -> 16.1.

The pick changes in 2.8-3.4% of events (8.6% on zjets); gained:lost is about
2:1 on every sample (zjets 6.08% : 3.08%, vbf 2.54% : 1.04%).

## Where the gain comes from (percent of TZP fails removed)

| step | zjets | dijet | vbf | ttbar |
|---|---:|---:|---:|---:|
| TZP selection, time re-weighted by pT e^-dz / sigma_t^2 inside the cluster | 3.5 | 5.4 | 8.6 | 7.2 |
| + mean-shift over all event tracks (60 ps) | 5.6 | 9.1 | 13.1 | 11.3 |
| + select on the kernel density at the mode | 6.5 | 10.6 | 13.6 | 12.3 |
| + sqrt(pT) in both weights | 8.1 | 12.5 | 13.9 | 13.5 |
| + HGTD agreement term | 8.5 | 13.8 | 18.1 | 15.1 |

**The time estimator is the larger half.** With the selection untouched, the
plain inverse-variance cluster mean leaves 6-13% of the fails on the table. An
oracle that times the selected cluster from its truth-HS tracks only removes
13-20%, so the mean-shift captures roughly two thirds of what re-timing can
give. This is item 1 of `classical_scoring_brainstorm.md`'s list, measured.

Plateaus (worst-sample / mean): kernel width 55-75 ps, 3 vs 6 passes identical,
pT power 0.35-0.65 or a 5 GeV cap, HGTD bonus 0.3-0.8, tau 30 vs 60. Re-scanning
TZP's alpha/beta/gamma on top moves nothing.

## Closed negatives from this study

- Expected-pileup subtraction under the cluster (`S - lambda W_ev G(t)`) and a
  tail boost `S / G(t)^k`: strongly negative at every strength.
- Per-track precision in the selection weight (`sigma_z0^-b`, a z0 likelihood
  ratio): negative or neutral. HS purity does rise steeply with z0 precision
  (98% at sigma_z0 < 0.1 mm against 24% above 1.6 mm), but the selection is
  already on its plateau in those weights.
- A kernel-weighted ("soft") cluster envelope instead of the hard cluster's:
  costs 4-7 points of the gain. The hard cluster's z centroid and d0 sum are
  better.
- In-jet boosts of the time weight, a sigma_t-scaled kernel, a box kernel,
  trimming at 1.5-2 sigma with inverse-variance weights: none beat the Gaussian.
- chi2/ndf, multiplicity and sigma_t powers on the score: within +-1%.

## Not done

- No grid run of the C++ (`clustering_hist --sample=`). The offline numbers
  come from the grid exports, and the C++ agrees with them on local VBF.
- `rpt_v5_hist`, `export_training_data` and the event display still use TZP.
- The HGTD term reads Athena's `RecoVtx_time`; the row without it is quoted
  above for a t0 that must not depend on Athena's.
