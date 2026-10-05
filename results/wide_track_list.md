# The wide track list for TZP_KDE_TZ (2026-10-05)

Follows `selection_reliability.md`. Question: about a tenth of Z+jets events
have no right-timed hard-scatter track in the clustering list. Are those tracks
absent, or cut?

Tools: `util/dump_tracks.cxx` (+ `condor/dump_tracks.sub`) writes every timed
track in HGTD acceptance inside a loose window (8 sigma, 4.5 below 1 GeV,
pT > 0.5 GeV, no ceiling, no quality cut) with the cut variables as columns.
`python/kde_study/wide.py` rebuilds any track list from it and reruns the full
chain (clustering, reliability weights, kernel score, mean-shift time);
`wide_tab.py`, `wide2..5.py` are the scans. With the nominal list it reproduces
the export result exactly (Z+jets 64,489 / 93,977). The evaluation is chunked
by event: unchunked it took 12-16 GB on 55k events.

## The rule

`getWideTracks` (`TrackFilterType::WIDE`), used by `TZP_KDE_TZ` and by the
`kde*` scenarios of `rpt_v5_hist`. Every other score keeps the nominal list.

| | nominal list | wide list |
|---|---|---|
| pT ceiling | 30 GeV | none |
| z0 significance, pT >= 1 GeV | < 3 | < 4 |
| pT floor | 1 GeV | 0.5 GeV, and below 1 GeV: >= 2 HGTD hits, significance < 3 |
| tracks per event | 26-30 | 73-77 |

## Result

Grid C++ (clusters 3262799-806, 92 shards, none held; m_jj >= 200):

| | vbf | zjets | dijet | ttbar |
|---|---:|---:|---:|---:|
| events | 630,720 | 64,950 | 105,137 | 863,552 |
| Athena HGTD | 84.53 | 41.60 | 77.88 | 77.24 |
| TZP | 91.69 | 64.48 | 88.30 | 88.50 |
| TZP_KDE_TZ, nominal list (previous) | 93.33 | 68.43 | 90.23 | 90.66 |
| **TZP_KDE_TZ, wide list** | **94.10** | **70.81** | **90.95** | **91.67** |
| TZP fails removed, previous | 19.8% | 11.1% | 16.5% | 18.8% |
| **TZP fails removed, now** | **29.0%** | **17.8%** | **22.6%** | **27.6%** |

Local VBF C++: 3,877 -> 2,796 fails (27.9%).

## Where the hard-scatter tracks were

Per event, timed tracks in acceptance (novbs; Z+jets / local VBF, one number
where they agree):

| category | tracks/event | HS purity | HS time right | of events with NO right-timed HS track in the nominal list, has one here |
|---|---:|---:|---:|---:|
| nominal list | 26 / 29 | 21% / 32% | 81% / 86% | 0 |
| pT > 30 GeV | 0.04 / 0.24 | 79% / 93% | 98% / 96% | 2% / 7% |
| z0 significance 3-4 | 6.8 | 1.9% / 3.7% | 77% / 86% | 3% / 6% |
| z0 significance 4-8 | 26 | 0.5% / 1.1% | 73% / 82% | 3% / 6% |
| pT 0.5-1 GeV, significance < 3 | 59 | 8% | 69% | **60% / 64%** |
| pT 0.5-1 GeV, significance 3-4.5 | 26 | 0.4% | 56% | 2% |

- **The 30 GeV ceiling removed the best tracks in the sample.** A quarter of a
  track per VBF event, 93% hard-scatter, time right 96%. Lifting it alone
  removes 9.0% of VBF's remaining fails (ttbar 4.5%, zjets and dijet 1.7-1.9%).
  The ceiling dates from linear-pT sums; with sqrt(pT) weights one hard track
  can no longer outvote a cluster.
- **Sub-GeV tracks are where the "no hard-scatter track" events were hiding.**
  They double the list for 8% purity, which a hard-cluster pT sum could not
  afford. Events with no right-timed HS track fall from 12.2% to 5.7% on
  Z+jets (dijet 2.3 -> 0.9, vbf 1.4 -> 0.6, ttbar 1.7 -> 0.5).
- Looser z association above 1 GeV is a small positive up to 4 sigma and flat
  beyond; with the ceiling in place it was worth nothing on Z+jets.

Steps, percent of the nominal-list configuration's fails removed (novbs dumps):

| | zjets | dijet | vbf | ttbar |
|---|---:|---:|---:|---:|
| no pT ceiling | 1.7 | 1.9 | 9.0 | 4.5 |
| + significance < 4 above 1 GeV | 1.8 | 2.1 | 10.7 | 6.7 |
| **+ 0.5-1 GeV tracks, >= 2 hits, significance < 3** | **7.7** | **7.1** | **12.2** | **12.0** |
| same, any hit count | 8.1 | 6.9 | 10.8 | 11.8 |

The two-hit requirement is the sample-independent choice: single-hit sub-GeV
tracks help Z+jets by 0.4 and cost VBF 1.4. Down-weighting sub-GeV tracks
(x0.5, x0.7) helps VBF by 0.3-0.6 and costs Z+jets 0.2-0.7; not adopted.
dijet and ttbar were not used to choose the rule.

By nominal-list HS track count, core fraction before -> after:

| n HS | zjets | dijet | vbf | ttbar |
|---|---|---|---|---|
| 0 | 19.0 -> 25.9 | 20.7 -> 26.7 | 22.1 -> 32.7 | 20.6 -> 31.9 |
| 1 | 38.3 -> 43.0 | 41.6 -> 45.9 | 48.1 -> 56.0 | 45.8 -> 53.1 |
| 2-3 | 57.1 -> 60.5 | 66.2 -> 69.2 | 76.1 -> 79.6 | 68.5 -> 72.8 |
| 4-7 | 78.1 -> 79.9 | 88.5 -> 89.3 | 92.6 -> 93.5 | 88.4 -> 89.8 |
| >= 8 | 94.0 -> 94.2 | 98.3 -> 98.3 | 98.9 -> 98.9 | 98.3 -> 98.4 |

mu = 0 (dumps of vbf_mu0, 138,245 events, and zeejets_mu0, 34,310): passing
events 138,029 -> 138,098 and 34,226 -> 34,257. A gain, not a cost.

The change is not free event by event: on Z+jets 6.0% of events are gained and
3.6% lost (vbf 1.9 / 1.1).

## Quality flag and RpT after the change

Grid C++, Q >= 1.5 (threshold unchanged):

| | vbf | zjets | dijet | ttbar |
|---|---:|---:|---:|---:|
| time provided | 91.7% | 62.3% | 87.6% | 87.2% |
| right among provided | 97.42 | 83.37 | 95.79 | 96.45 |
| core over all events (withheld = fail) | 89.31 | 51.95 | 83.87 | 84.15 |
| same, previous | 87.36 | 47.49 | 82.20 | 81.80 |
| Athena: provided / right / over all | 91.6 / 92.2 / 84.5 | 56.0 / 74.3 / 41.6 | 86.9 / 89.6 / 77.9 | 85.7 / 90.1 / 77.2 |
| mean cluster purity, no cut -> cut | 71.7 -> 75.1 | 40.3 -> 49.8 | 66.9 -> 71.9 | 66.5 -> 71.2 |

The flag now passes more events at the same purity. Mean cluster purity is 1-2
points lower than before on vbf / dijet / ttbar: the selected hard cluster
absorbs soft pileup tracks, which the time weights discount.

RpT, forward jets > 40 GeV, rejection at matched 0.80 efficiency (tzp / kde):
vbf 156 / 159, zjets 77.7 / 78.8, dijet 84.6 / 84.1, ttbar 87.8 / 91.4. VBF R1
162 / 172, R2 125 / 134. Level to +8%; the tagger still sees little of the
core-fraction gain. Measured kde inflation on this run: vbf 1.36, zjets 1.58,
ttbar 1.41 (in use 1.33 / 1.65 / 1.46).

## Not done

- Tracks below 0.5 GeV (7.7 per event, 2.6% HS, time right 60%) were not tried.
- The wide list feeds TZP_KDE_TZ only. Whether TZP itself wants any of it is
  untested; the ceiling almost certainly, the sub-GeV tracks probably not.
- `export_training_data` and the event display still use the nominal list.
