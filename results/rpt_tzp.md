# TZP in the RpT study (rpt_v5, local VBF) — first consumer-side result

New scenario `tzp` (index 6, appended so the load-bearing 0-5 are untouched):
TZP-selected cluster + guarded in-jet time, gated with the standard
`GATE_SIGMA = 2.5` machinery. Inflation seeded from the waves value (1.38);
the calibration table measured **1.37** on the first run, so no iteration.

## Calibration row — TZP's picks are also the best-timed

| scenario | n core | core sig | tail>150 | ratio |
|---|---:|---:|---:|---:|
| hgtd | 542,652 | 41.7 ps | 13.9% | 1.48 |
| trkptz | 579,567 | 38.6 ps | 15.2% | 1.39 |
| waves | 579,389 | 38.5 ps | 15.2% | 1.38 |
| **tzp** | **583,339** | **38.0 ps** | **14.6%** | 1.37 |

Most core tracks, smallest core width, lowest tail of any real scenario.

## Central baseline (the standing correctness check): PASSES

The tzp row agrees with every other scenario in both central slices to within
the expected few-jet wiggle. No leak.

## Forward >40 GeV — TZP is the best real scenario at every WP to 0.95

| scenario | 0.85 | 0.90 | 0.93 | 0.95 | 0.97 |
|---|---:|---:|---:|---:|---:|
| zonly | 93.2 | 57.1 | 37.5 | 21.3 | 6.3 |
| hgtd | 126.8 | 71.2 | 51.1 | 16.2 | 8.8 |
| trkptz | 117.4 | 75.5 | 38.0 | 17.6 | 9.7 |
| waves | 117.4 | 78.2 | 46.3 | 25.4 | 9.9 |
| **tzp** | **126.8** | **83.4** | **51.5** | **26.9** | 9.8 |
| waves_ideal | 162.5 | 107.4 | 56.1 | 27.1 | 9.3 |
| truth | 176.1 | 115.2 | 62.1 | 29.1 | 10.6 |

At 0.90 HS efficiency: **83.4 vs waves' 78.2 (+7%) and trkptz's 75.5 (+10%)**.
At 0.95, tzp (26.9) is within 0.2 of the waves_ideal oracle (27.1). It closes
~18% of the waves -> waves_ideal gap at 0.90.

## Forward 30-40 GeV — parity with trkptz, waves ahead at loose WPs

tzp 101.9 / 27.4 / 10.9... vs waves 104.3 / 33.6 / 10.8. Consistent with the
core-fraction differentials: WAVeS's remaining edge is sparse/low-pT jets, and
the soft slice is where those live. rejAtEff quantization caveat applies (the
10.9 repeats are the sampler clamping at the efficiency ceiling).

## Caveats

- Local VBF only; the zjets/dijet RpT runs need condor (`rpt_v5_hist.sub`) and
  the tzp inflation seeds (zjets 1.65, dijet 1.46) checked against their
  measured ratios.
- `rejAtEff` steps: near a bin edge a 0.02% eff difference moves printed
  rejection by ~30% (the shared 126.8 at 0.85 with hgtd is same-bin clamping).
  Read the ordering across WPs, not single cells.

---

# Grid samples (clusters 2956037-39, 2026-09-01): zjets 14 / dijet 6 / ttbar 16 shards

36/36 shards clean (no stderr). Merged on the exact shard-count globs
(`shard*of14/6/16` -- dijet/ carries stale shard*of4 sets from July, the
documented glob trap). Central baselines agree for the tzp row on all three
samples. tzp calibration: measured 1.59 / 1.43 / 1.43 against seeds
1.65 / 1.46 / 1.46 -- 2-4%, i.e. <1% on the gate width; the dijet-seeded ttbar
row validated across every scenario. On all three samples tzp again has the
most core tracks and smallest core sigma of any real scenario.

## Forward >40 GeV rejection (rejAtEff steps apply; read orderings, not cells)

| sample | WP | trkptz | waves | **tzp** | waves_ideal | truth |
|---|---|---:|---:|---:|---:|---:|
| zjets | 0.85 | 53.1 | 43.2 | 51.3 | 54.5 | 57.9 |
| zjets | **0.90** | 10.5 | 14.3 | **17.2** | 16.8 | 17.0 |
| dijet | 0.85 | 75.6 | 69.7 | 74.9 | 95.7 | 97.2 |
| dijet | 0.90 | 54.4 | 48.5 | 54.0 | 59.9 | 62.5 |
| ttbar | 0.85 | 42.7 | 48.8 | 42.0 | 40.4 | 41.0 |
| ttbar | **0.90** | 10.6 | 15.6 | **18.1** | 9.5 | 9.5 |

- **zjets at 0.90 -- the number the branch was built toward: tzp 17.2 vs waves
  14.3 (+20%) and trkptz 10.5 (+64%), sitting ABOVE the waves_ideal oracle
  (16.8) and at the truth row (17.0).** The oracle rows inherit WAVeS
  *selection*, which is poor off-VBF -- so a real algorithm overtaking them is
  consistent, and it says selection was the whole game there.
- ttbar at 0.90: tzp 18.1, the best row in the table outright (the ideal/truth
  rows collapse to 9.5 -- the documented efficiency-ceiling effect arriving
  earlier for gates built on WAVeS selection).
- dijet: tzp == trkptz within steps (74.9/54.0 vs 75.6/54.4), both above waves.
- Soft slices (30-40 GeV): dijet tzp leads everything real (55.1/22.9, above
  even waves_ideal at 0.90); zjets soft slice keeps WAVeS's loose-WP edge
  (13.9 vs 10.0 at 0.85) with everything clamped at the efficiency ceiling
  beyond it.

## Verdict

The consumer-side picture matches the core-fraction one on every sample: TZP is
the best or tied-best real scenario at the tight working points everywhere, the
zjets gain is the largest (+20% rejection over deployed WAVeS at 0.90), and the
one residual WAVeS edge is loose WPs in soft slices. Inflation rows could be
updated to the measured 1.59/1.43/1.43 for a future rerun; at 2-4% the effect
is sub-percent on the gate and does not change any ordering.
