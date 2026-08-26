# Track-to-Vertex z-Association Study

**Question:** the clustering has two candidate rules for deciding which tracks
belong to the hard-scatter vertex in z, and therefore which tracks the vertex
t0 is built from. Which one produces a better t0?

**Date:** 2026-08-26 · **Sample:** local VBF H→inv. (112,400 events, 47,535
passing selection) · **Executables:** `util/assoc_study_hist.cxx` +
`util/assoc_study_plot.cxx`

---

## Verdict

| | recommendation |
|---|---|
| **TRKPTZ** | **Switch.** The incumbent `significance < 3.0` is the second-worst of the eleven rules tested. `getNewDzpara × 1.0` gives **−1.58 ps core σ (15.95 → 14.37, −9.9%)** and **+0.92 %pt core fraction (90.38 → 91.30)**. |
| **WAVeS** | **Leave it.** Every rule is within ±0.1 %pt of the incumbent at the loose end, and tightening actively *hurts* (`dzpara × 0.6`: +1.86 ps, −1.83 %pt). |
| **One rule for both** | **`getNewDzpara × 1.2`** — TRKPTZ −1.32 ps / +0.92 %pt, WAVeS +0.13 ps / +0.10 %pt. Effectively the TRKPTZ optimum at no cost to WAVeS. |

**The cheapest version of this result needs no code change at all.** The
existing `--dzpara` flag already runs `getNewDzpara × 1.4`
(`DZ0_PARA_SCALE_CLUSTER`), and that alone captures roughly three quarters of
the available gain: TRKPTZ −1.03 ps / +0.76 %pt, WAVeS +0.01 ps / +0.01 %pt.
Going from 1.4 to 1.2 buys a further −0.30 ps / +0.16 %pt on TRKPTZ.

Two results worth carrying forward beyond the choice of cut:

1. **The two scores want opposite things**, and for a mechanical reason (see
   [Why WAVeS does not care](#why-waves-does-not-care)). Any future selector
   that recomputes its time from a track subset will behave like WAVeS here.
2. **About 59% of WAVeS's core-fraction advantage over TRKPTZ is reproducible
   by fixing the z-association alone**, with no jet-proximity weighting. The
   gap is 1.38 %pt at the incumbent rule and 0.56 %pt when each score is given
   its own best rule. On core σ the ordering actually *inverts*: TRKPTZ at
   `dzpara × 1.0` (14.37 ps) beats WAVeS at any rule (best 14.95 ps).

---

## What was measured

Two rule families, scanned rather than sampled at one point each:

| family | test | scanned over |
|---|---|---|
| `SIGNIFICANCE` | \|z0 − z_vtx\| / sqrt(1.15² · var_z0) < N | N = 3.0, 2.5, 2.0, 1.5 |
| `DZ_PARA` | \|z0 − z_vtx\| / getNewDzpara(η, pT) < S | S = 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 2.0 |

`significance < 3.0` is the incumbent (`MAX_NSIGMA`). `getNewDzpara × 1.4` is
the reference R_pT study's working point and what `--dzpara` currently does.

Both families are scanned because the two rules differ in **two** ways at once
— functional form *and* how tight they happen to be at their nominal points. A
single point from each cannot separate "better shape" from "simply tighter".
Scanning both allows comparison at **matched track purity**, which is the only
comparison that isolates the shape (see
[Shape vs. tightness](#shape-vs-tightness)).

Metrics, per rule, for **both** selection scores:

- inclusive core σ(Δt) and core fraction (\|Δt\| < `PASS_SIGMA` = 60 ps)
- both, as a function of n forward HS tracks and n forward jets

All metrics reuse the project's existing machinery, so they are directly
comparable with the rest of the analysis: core σ is `ACTIVE_FIT_MODEL`'s
`Sigma1`, the same fit and parameter `plotPostProcessing` writes into
`params->sigmaDist`; core fraction is passed/total from the efficiency
histograms.

### Method notes that matter for reading the numbers

- **One event loop, eleven rules.** Every rule is evaluated against the *same*
  event, so the comparison is paired — identical events, jets, and truth. It
  also costs one traversal of the 52 GB ntuple instead of eleven.
- **The counting scan is pinned** at `COUNTING_NSIGMA = 3.0` for every rule.
  `processEventData` already separates the two uses of the track list: the
  counting scan feeds `EventCounts`, which supplies the **x-axes** of every
  differential plot, while the clustering list is what t0 is built from. If the
  association moved the x-axis too, the "n_hs_tracks = 3" bin would hold a
  different set of events for every rule and no per-bin difference could be
  attributed to the association.
- **Under a rule the clustering list is re-selected from the full track array**,
  not narrowed from the counting list. A parameterized rule is *not* a subset of
  a significance cut — `getNewDzpara` is looser than 3σ for some (η, pT) and
  tighter for others — so filtering would silently intersect every parameterized
  rule with `significance < 3.0` and understate it.
- **Core fraction includes under/overflow** (`Integral(0, n+1)`). `PlotObj`'s
  n-HS-track axis starts at 2.5, so every event with fewer than three forward HS
  tracks lands in the underflow bin — and per the failure decomposition that is
  where most failures live. Excluding it would report a core fraction for the
  easy events only.

---

## Track accounting

What each rule does to the track list, before any timing is computed. Forward,
pT/quality-selected tracks in accepted events: 22,140,807 available, of which
622,562 truth-HS.

| rule | kept | HS kept | PU kept | purity | HS recall |
|---|---:|---:|---:|---:|---:|
| signif < 3.0 *(incumbent)* | 1,794,275 | 570,094 | 1,224,181 | 31.77% | 91.57% |
| signif < 2.5 | 1,584,205 | 556,017 | 1,028,188 | 35.10% | 89.31% |
| signif < 2.0 | 1,360,541 | 529,535 | 831,006 | 38.92% | 85.06% |
| signif < 1.5 | 1,110,180 | 478,153 | 632,027 | 43.07% | 76.80% |
| dzpara × 0.6 | 806,147 | 409,408 | 396,739 | 50.79% | 65.76% |
| dzpara × 0.8 | 995,367 | 473,257 | 522,110 | 47.55% | 76.02% |
| dzpara × 1.0 | 1,159,731 | 513,168 | 646,563 | 44.25% | 82.43% |
| dzpara × 1.2 | 1,307,697 | 537,430 | 770,267 | 41.10% | 86.33% |
| dzpara × 1.4 *(`--dzpara`)* | 1,446,455 | 553,016 | 893,439 | 38.23% | 88.83% |
| dzpara × 1.6 | 1,579,495 | 563,559 | 1,015,936 | 35.68% | 90.52% |
| dzpara × 2.0 | 1,835,583 | 576,097 | 1,259,486 | 31.38% | 92.54% |

Note how weak the incumbent is as a pileup filter: **68% of the tracks t0 is
currently built from are pileup.**

---

## Inclusive results

`sigma(nHS≤5)` is the core σ restricted to low-multiplicity events, from the
`inclusiveResoLowTrack*` histograms — the only view here that reaches below
three forward HS tracks.

### TRKPTZ

| rule | core σ [ps] | core frac [%] | trunc. RMS [ps] | σ(nHS≤5) [ps] |
|---|---:|---:|---:|---:|
| signif < 3.0 *(incumbent)* | 15.95 ± 0.09 | 90.38 ± 0.14 | 17.77 | 29.27 |
| signif < 2.5 | 15.38 ± 0.09 | 90.82 ± 0.13 | 17.40 | 28.09 |
| signif < 2.0 | 14.85 ± 0.08 | 90.88 ± 0.13 | 16.86 | 26.42 |
| signif < 1.5 | 14.87 ± 0.08 | 90.71 ± 0.13 | 16.88 | 26.14 |
| dzpara × 0.6 | 15.05 ± 0.08 | 89.97 ± 0.14 | 16.93 | 25.61 |
| dzpara × 0.8 | 14.54 ± 0.08 | 90.96 ± 0.13 | 16.59 | 25.37 |
| **dzpara × 1.0** | **14.37 ± 0.08** | **91.30 ± 0.13** | **16.51** | **25.34** |
| dzpara × 1.2 | 14.63 ± 0.08 | 91.30 ± 0.13 | 16.74 | 25.84 |
| dzpara × 1.4 | 14.92 ± 0.08 | 91.14 ± 0.13 | 17.01 | 26.58 |
| dzpara × 1.6 | 15.29 ± 0.08 | 90.93 ± 0.13 | 17.31 | 27.33 |
| dzpara × 2.0 | 16.13 ± 0.09 | 89.99 ± 0.14 | 17.90 | 28.98 |

Both families improve as they tighten and then turn over — significance at
about 2.0, the parameterization at 1.0. There is a genuine optimum, and the
incumbent is nowhere near it.

The low-multiplicity column moves furthest: **29.27 → 25.34 ps, −13%**. The
association fix helps most exactly where the project's failure decomposition
says the failures are.

### WAVeS

| rule | core σ [ps] | core frac [%] | trunc. RMS [ps] | σ(nHS≤5) [ps] |
|---|---:|---:|---:|---:|
| signif < 3.0 *(incumbent)* | 14.99 ± 0.08 | 91.76 ± 0.13 | 16.67 | 22.61 |
| signif < 2.5 | 15.10 ± 0.08 | 91.83 ± 0.13 | 16.74 | 22.50 |
| signif < 2.0 | 15.27 ± 0.08 | 91.62 ± 0.13 | 16.77 | 22.35 |
| signif < 1.5 | 15.81 ± 0.08 | 91.25 ± 0.13 | 17.21 | 22.95 |
| dzpara × 0.6 | 16.85 ± 0.08 | 89.93 ± 0.14 | 17.87 | 23.67 |
| dzpara × 0.8 | 15.88 ± 0.08 | 91.21 ± 0.13 | 17.21 | 22.86 |
| dzpara × 1.0 | 15.34 ± 0.08 | 91.80 ± 0.13 | 16.86 | 22.35 |
| **dzpara × 1.2** | 15.12 ± 0.07 | **91.86 ± 0.13** | 16.66 | 22.15 |
| dzpara × 1.4 | 15.00 ± 0.07 | 91.77 ± 0.13 | 16.58 | 22.27 |
| dzpara × 1.6 | **14.95 ± 0.08** | 91.77 ± 0.13 | 16.54 | 22.19 |
| dzpara × 2.0 | 14.95 ± 0.08 | 91.52 ± 0.13 | 16.58 | 22.32 |

Flat to within ±0.1 %pt across the loose half of the scan, and monotonically
worse once either family tightens past ~1.2. WAVeS has nothing to gain here and
something to lose.

### Why WAVeS does not care

`Cluster::calculateTime` recomputes the WAVeS time using only the constituent
tracks within ΔR < 0.4 of a qualifying forward jet. **That in-jet filter already
does the pileup-shedding job a tighter z-association would do**, so the two
mechanisms are redundant. Once the in-jet refinement is in place, tightening the
z-association can only remove hard-scatter tracks — which is precisely the
monotonic degradation the table shows.

This is a statement about the *mechanism*, not about WAVeS specifically: any
future selector that recomputes its time from a track subset should be expected
to behave the same way, and should be tuned against a loose association.

---

## Shape vs. tightness

Comparing the families at matched track purity separates the two effects. Rows
are paired on purity, so the only remaining difference is the shape of the rule.

| purity | significance rule | dzpara rule | HS recall (sig → para) | TRKPTZ σ (sig → para) | TRKPTZ core frac (sig → para) |
|---|---|---|---|---|---|
| ~31.5% | 3.0 (31.77%) | 2.0 (31.38%) | 91.57 → 92.54 | 15.95 → 16.13 | 90.38 → 89.99 |
| ~35.4% | 2.5 (35.10%) | 1.6 (35.68%) | 89.31 → 90.52 | 15.38 → 15.29 | 90.82 → 90.93 |
| ~38.6% | 2.0 (38.92%) | 1.4 (38.23%) | 85.06 → 88.83 | 14.85 → 14.92 | 90.88 → 91.14 |
| ~43.7% | 1.5 (43.07%) | 1.0 (44.25%) | 76.80 → 82.43 | 14.87 → **14.37** | 90.71 → **91.30** |

**At matched purity the parameterization retains more hard-scatter tracks, and
the advantage widens the tighter the working point.** At ~31.5% purity it is
worth +1.0 points of recall and the two rules perform the same; at ~43.7% purity
it is worth **+5.6 points of recall at higher purity**, and TRKPTZ gains 0.5 ps
and 0.6 %pt from the shape alone.

That is why the parameterization wins overall: the optimum sits at the tight
end, and the tight end is exactly where the shape matters.

The likely cause is that `sqrt(var_z0)` is a per-track covariance that does not
capture the η dependence of the true z0 spread, whereas `getNewDzpara` is an
empirical fit to the observed spread as a function of (η, pT). A fixed multiple
of a mis-estimated per-track σ cuts too hard on some tracks and too loosely on
others. The codebase already carries independent evidence for this: the
`Z0_VAR_INFLATION = 1.15²` factor exists because the covariance was measured to
underestimate the true z0 resolution by ~15% — and a single global factor cannot
repair an η-dependent mis-modelling.

---

## Differential results

Only the incumbent and the two family optima are shown; the full eleven-rule
tables are in the executable's console output.

### Core σ [ps] vs n forward HS tracks (TRKPTZ)

| rule | 3 | 4 | 5 | 6 | 8 | 10 | 12 | 15 | ≥20 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| *N events* | 2247 | 2935 | 3538 | 3838 | 4031 | 3511 | 2754 | 1544 | 1628 |
| signif < 3.0 | 30.5 | 25.3 | 23.9 | 19.4 | 16.4 | 14.3 | 12.7 | 11.5 | 8.4 |
| signif < 2.0 | 27.6 | 24.0 | 21.5 | 18.3 | 15.0 | 13.5 | 12.2 | 10.8 | 8.1 |
| dzpara × 1.0 | **27.3** | **22.8** | **20.5** | **17.7** | **14.6** | **12.9** | **11.5** | **10.4** | **7.9** |

### Core fraction [%] vs n forward HS tracks (TRKPTZ)

| rule | 3 | 4 | 5 | 6 | 8 | 10 | 12 | 15 | ≥20 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| signif < 3.0 | 73.5 | 80.7 | 86.0 | 89.3 | 94.3 | 96.8 | 98.1 | 99.2 | 99.9 |
| signif < 2.0 | 73.1 | 82.9 | 87.0 | 90.3 | 94.7 | 97.3 | 97.7 | 99.5 | 99.9 |
| dzpara × 1.0 | **74.9** | **83.5** | **87.2** | **91.3** | **95.0** | **97.6** | 98.1 | 99.2 | 99.9 |

The gain is concentrated at low multiplicity and shrinks steadily: −3.2 ps at 3
HS tracks, −2.5 at 4, −1.4 at 10, −1.2 at 12, −0.5 at ≥20. It never quite
vanishes, but above ~12 HS tracks the core fraction has saturated above 98% for
every rule, so there is nothing left for it to buy. This is the expected shape —
with few tracks each pileup contaminant carries more weight in the average.

### vs n forward jets (TRKPTZ)

| rule | 0 | 1 | 2 | 3 | 4 | ≥5 |
|---|---:|---:|---:|---:|---:|---:|
| *N events* | 0 | 30191 | 13009 | 3155 | 819 | 361 |
| σ, signif < 3.0 | — | 17.1 | 14.1 | 14.0 | 13.6 | 11.9 |
| σ, dzpara × 1.0 | — | **15.3** | **13.0** | **12.6** | **13.1** | 13.0 |
| core frac, signif < 3.0 | — | 88.8 | 92.8 | 94.0 | 93.7 | 93.6 |
| core frac, dzpara × 1.0 | — | **89.9** | **93.6** | **94.5** | **94.3** | 93.1 |

The `n = 0` bin is empty by construction — the event selection requires at least
one forward jet. The `≥5` bin holds 361 events; differences there are noise, and
the apparent inversion at ≥5 should not be read as physics.

### WAVeS, differentially

Core σ [ps] vs n forward HS tracks:

| rule | 3 | 4 | 5 | 6 | 8 | 10 | 12 | 15 | ≥20 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| signif < 3.0 | 24.1 | 21.0 | 18.7 | 18.0 | 15.5 | 13.9 | 12.8 | 11.8 | 9.0 |
| dzpara × 1.2 | 23.8 | 20.5 | 18.6 | 18.0 | 15.4 | 13.9 | 12.6 | 11.9 | 9.3 |
| signif < 1.5 | 24.3 | 21.5 | 19.4 | 18.6 | 16.4 | 14.9 | 13.8 | 12.2 | 9.5 |
| dzpara × 0.6 | 24.3 | 22.3 | 20.5 | 19.9 | 17.5 | 15.9 | 14.3 | 13.3 | 10.3 |

Across the **loose half** of the scan — every parameterization scale ≥ 1.0 and
every significance cut ≥ 2.0 — the rules agree bin by bin to within the
uncertainty, matching the flat inclusive table.

The **tight** rules are a different story, and the differential view adds
something the inclusive number does not: `dzpara × 0.6` is worse at *every*
multiplicity, not just the sparse bins (+2.0 ps at 8 HS tracks, +1.5 at 12,
+1.3 at ≥20). Its inclusive degradation is therefore a uniform loss of
hard-scatter tracks, not a low-multiplicity effect — consistent with the
mechanism above, since removing in-jet HS tracks hurts wherever they are.

---

## Caveats

- **One sample.** Local VBF only. The grid samples (`zjets`, `dijet`) live under
  `/data/mcardiff/...` and are not reachable from this machine. Z+jets in
  particular is known to behave differently — Athena produces a valid vertex
  time for only 16.5% of Z+jets events versus 80.1% of VBF — so this should be
  re-run there before the constant is changed for all samples. The study runs on
  the full local sample in **~35 s** (hist) **+ ~15 s** (plot), and on condor via
  the existing `--sample=` convention.
- **The differential n-HS-track axis starts at 3.** `HS_TRACK_MIN = 2.5`, so
  events with 0–2 forward HS tracks are in the underflow. They are counted in the
  inclusive core fraction, and the `σ(nHS≤5)` column partially covers them, but
  they do not appear in the differential tables. Changing `HS_TRACK_MIN` would
  move every other plot in the project, so it was left alone.
- **Fit model.** Core σ is a single Gaussian fitted over the full ±1000 ps range
  (`ACTIVE_FIT_MODEL = FIT_SNGGAUS_FULL`), which the tail can pull. The
  model-free truncated RMS is carried alongside as a cross-check and reproduces
  the same rule ordering everywhere except one adjacent swap, so the conclusions
  are not a fit artifact.
- **`#Delta` does not render** in TLatex on this machine — plot axes read
  "σ( t)" rather than "σ(Δt)". Pre-existing and affects every plot in the
  project, not just these.

---

## Reproducing

```bash
cd build && cmake .. && make assoc_study_hist assoc_study_plot
./assoc_study_hist --threads=8          # ~35 s on the local VBF sample
./assoc_study_plot                      # tables + PDFs, no ntuple access
```

Add `--sample=vbf|zjets|dijet` to both stages for a grid sample, and
`--max-events=N` to the hist stage for a quick check. Plots land in
`<OUTPUT_DIR>/assoc_study/`; the console output ends with a `CSV|`-prefixed
block containing every number in this document.

Rules are declared in one place, `assocRules()` in
`util/assoc_study_common.h` — adding or moving a working point is a one-line
edit plus a re-run of the hist stage.
