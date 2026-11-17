# CLAUDE.MD - Vertex Timing Analysis Project

## Project Overview

This is a High-Granularity Timing Detector (HGTD) analysis project for ATLAS particle physics research. The project focuses on vertex timing reconstruction and clustering algorithms for analyzing collision events in high-pileup environments.

## Core Functionality

The project performs timing-based vertex reconstruction using:
- HGTD detector timing information
- Track-based timing measurements
- Jet clustering and filtering algorithms
- Multi-variable analysis for vertex identification

## Project Structure

Run `ls` for the layout; each executable's header describes it. Non-obvious points:
- `src/clustering_dt.cxx` and `util/rpt_v5.cxx` are legacy monolithic executables kept only as
  validation ground truth; the `*_hist` / `*_plot` pairs are the primary workflow.
- `python/` scripts are run from within `python/` (the event display invokes its macro by relative name).
- `slides/`: hand-written Beamer talks reading figures in place from `../figs` and `../condor`
  (`vvhbeamer.sty` is a copy of `~/project/VBS_HIGGS/slides`' style). Build with
  `cd slides && latexmk -pdf -outdir=build <deck>.tex`; `grep -c 'missing plot' build/<deck>.log`
  must be 0 — a missing PNG is a grey box, not an error.
- `results/`: one Markdown/JSON/HTML record per study, positive AND negative results.
- `share/models/`, `share/scripts/`: historical ML artifacts with no active consumer.
- `util/scratch/`: one-off macros (e.g. the region-tree friend-file macros in `util/scratch/region_weights/`).
- `condor/vbs_region_diag.sub` needs its `jvt`/`jvttag` macros to agree; `condor/vbs_time_veto.sub`
  must read VBF from the RAW ntuples (`inputargs=--ntuple-dir=...`), Z+jets from its skim.

## Main Executables

Since the histogramming/plotting split (see "Histogramming vs. Plotting
Split" below), the two primary analyses each ship as a **pair** of
executables: a `*_hist` stage that runs the `TTreeProcessorMT` event loop
against the ntuples and writes every histogram to a ROOT file, and a
`*_plot` stage that reads that file back and produces the PDFs, touching no
ntuples at all. The original single-executable versions (`clustering_dt`,
`rpt_v5`) are kept building as a reference/fallback but are no longer the
primary workflow.

### Primary Analysis
- **clustering_hist.cxx** / **clustering_plot.cxx** (`src/`): Split of the
  main clustering analysis. `clustering_hist` runs the event loop (only the
  HGTD timing scenario is currently built/processed — `Scenario::IDEAL_RES`/
  `IDEAL_EFF` and `buildAnalysisMap`'s branching for them still exist for
  future reuse, but `main()` no longer instantiates or iterates those maps)
  and writes its histogram file via `MyUtl::histFilePath()` (see
  "Sample-prefixed output paths" below). `clustering_plot`
  reads that file and produces the `theplot_*`, `waves_oracle_*`, and
  `hgtd_trkptz_waves_*` comparison PDF groups (see Score Registry section),
  plus inclusive resolution/pull plots and shape-comparison plots — no
  ntuple access, runs in seconds. Active scores: `TRKPTZ`, `WAVES`,
  `JET_T_REFINED`, `TEST_MISAS`, `WAVES_MISCL`, `WAVES_MISAS`, plus `HGTD`.
  Iterative clustering (`ClusteringMethod::ITERATIVE`, `DIST_CUT_CONE =
  3.0σ`) is the default main collection. Shared booking/plotting logic
  (`Scenario` enum, `buildAnalysisMap`, `makeComparisonPlots`) lives in
  `src/clustering_dt_common.h`.
- **clustering_dt.cxx** (`src/`): Legacy monolithic executable (event loop +
  plotting in one run) that the split above was validated against
  bin-for-bin/pixel-for-pixel; still builds, but `clustering_hist` +
  `clustering_plot` is the primary workflow going forward.

### Diagnostics & Failure Analysis
- **failure_decomposition.cxx** (`util/`): Assigns every TRKPTZ-failing event
  to a single principal cause via a fixed precedence ladder (see Failure
  Decomposition section below).

### RpT / Report Studies
- **rpt_v5_hist.cxx** / **rpt_v5_plot.cxx** (`util/`): Split of the current
  RpT (jet-level pileup-discrimination) study. `rpt_v5_hist` runs the event
  loop (same `TTreeProcessorMT`/`ThreadState` structure as the legacy
  version) and writes its histogram file via `MyUtl::histFilePath()` (see
  "Sample-prefixed output paths" below). `rpt_v5_plot` reads that file and
  does `generate_roc_grid`, `rejAtEff`, the harmonized `[0.80, 1.0]` ROC
  window (see below), the console summary tables, and all PDF drawing — no
  ntuple access. Shared `Scenario`/binning/save-load logic lives in
  `util/rpt_v5_common.h`.

  **Seven scenarios**, in `makeScenarios()` order — the index is load-bearing,
  since `fillJets` fills `sv[0..6]` positionally:

  | idx | name | what it is |
  |---|------|------------|
  | 0 | `zonly`       | ITk-only, no time gate |
  | 1 | `hgtd`        | Athena `RecoVtx_time` |
  | 2 | `trkptz`      | TRKPTZ-selected cluster time |
  | 3 | `waves`       | WAVeS-selected cluster time (in-jet refined) |
  | 4 | `waves_ideal` | ORACLE: re-clustered with 30 ps track times AND the cluster closest in time to truth |
  | 5 | `truth`       | vertex = Gaus(truth HS vtx, 10 ps), tracks = Gaus(own truth time, 30 ps) |
  | 6 | `tzp`         | TZP-selected cluster + guarded in-jet time (appended last so 0–5 stay stable) |

  "Own truth time" is the truth PARTICLE's production time where the track has
  one, else its truth vertex time, else a Gaus(t_HS, 175 ps) pileup draw.
  `waves_ideal` never uses the WAVeS score: it picks the smeared-world cluster
  closest to truth, and its t0 is that cluster's plain mean (no in-jet
  refinement). Both idealised worlds are generated in `src/idealised_timing.h`
  (seed, draw order, re-clustering, truth-closest pick), shared with the
  event-display macro so a display replays exactly the times the histograms
  used.

  **The `truth`/`waves_ideal` rows are NOT upper bounds for the real rows.**
  Both live in a re-smeared world — every timed track redrawn at a flat 30 ps,
  which is WORSE than HGTD's real per-hit resolutions (35/25/21 ps for 1/2/3
  hits) — and each idealises two things at once (`truth`: t0 AND track times;
  `waves_ideal`: track times AND selection). On ttbar the real `tzp` row
  genuinely sits above both at matched efficiency.

  **Never compare `rejAtEff` cells across scenarios near the efficiency
  ceiling.** The cells snap to bins whose ACTUAL efficiencies can sit up to
  ~1.6% apart, which near the ceiling moves rejection by factors — an early
  "+20% over waves on zjets" headline from those cells inverted to "~10–15%
  BEHIND waves" once traced at matched efficiency (results/rpt_tzp.md carries
  the correction). Trace the (eff, rej) curves and compare at matched
  efficiency instead.

  **Core fraction and RpT-gate quality are NOT the same objective.** On zjets
  WAVeS has the worst core fraction of any selector (57.9%) and the best real
  matched-efficiency RpT rejection: the 2.5σ gate passes untimed tracks
  unconditionally and empties the timed component on a catastrophically wrong
  t0, so it forgives WAVeS's failure mode and punishes slightly-wrong t0s.
  What the RpT gate actually rewards in a t0 is an open question.

  **Four jet populations**, each filled for all seven scenarios: two forward pT
  slices (`_lo` 30–40 GeV, `_hi` >40 GeV, both 2.4 < |eta| < 3.8), the same
  two at central |eta| < 2.4 (`_lo_cen`/`_hi_cen`), and the two VBS regions
  (`_r1`/`_r2`, both forward legs 2.4 < |eta| < 3.8 -- the NARROW window).
  The central pair is an ITk-only baseline AND a standing
  correctness check: |eta| < 2.4 is outside HGTD acceptance, so all rows
  must agree there. That check has caught two real bugs (see Physics
  Findings).

  The wide-window region tree (`<prefix>rpt_v5_regions.root`, TTree-only: merge with `hadd`),
  its friend files, `--jvt=none|loose|tight`, and the region R_pT distribution / display scripts are
  documented in the `vbs-timing-studies` skill.

  Both idealised rows only smear tracks with `trackTimeValid == 1` — idealising
  the time of a track HGTD never measured invents coverage it does not have,
  which the central baseline exposes immediately. Their smeared times are
  generated ONCE per track and shared by the re-clustering and both gates, so
  a row never clusters on one draw and gates on another. The per-thread RNG is
  re-seeded per event from event-stable quantities, so results do not depend on
  how `TTreeProcessorMT` distributes entries. `clusterTracksInTime`'s own
  `useSmearedTimes` path is NOT used: it calls `getSmearedTrackTime`, which
  draws from the global `gRandom` and is a data race under MT.
- **generate_rpt.cxx**, **rpt_v2.cxx**..**rpt_v4.cxx** (`util/`): Earlier
  iterations of the same study (z-only, real-HGTD, TRKPTZ-selected timing
  gates), monolithic (event loop + plotting), not split.
- **rpt_v5.cxx** (`util/`): Legacy monolithic executable that
  `rpt_v5_hist`/`rpt_v5_plot` was validated against; still builds, but the
  split pair is the primary workflow going forward.

### Event-level timing veto (2026-09-28)
- **vbs_time_veto.cxx** (`util/`) + `python/vbs_time_veto_{plot,sr,common}.py`, `vbs_hspu_study.py`,
  `vbs_sr_minimal.py`, `vbf_hinv_sr.py`: method, validation and results are in the
  `vbs-timing-studies` skill and `results/vbs_time_veto*.md`. Two rules that must not be missed:
  VBF must be read from the RAW ntuples (its recoil cut is truth MET and the skims drop the truth
  record), and S/sqrt(B) only ever compares S and B measured in the same selection.

### ML Data Export
- **export_training_data.cxx** (`util/`): Exports the cluster-SELECTION
  training data as a ROOT file holding **three** TTrees — `clusters` (one row
  per candidate cluster, ~96 columns), `tracks` (one row per track, ~31), and
  `jets` (one row per forward jet, 13). Event selection is reco-only, so the
  training population is reproducible at inference time. Every truth-derived
  column is `truth_`-prefixed and must be held out of the model inputs;
  `delta_t` (cluster time − truth HS time) is the target.

  **Group key is `(sample_id, file_idx, event_num)`.** Each model here ranks
  the candidates of ONE event, so a colliding key silently merges two events'
  candidates into one ranking group — a training bug with no symptom, since
  every row is individually valid. `event_num` is the entry within its file
  and `file_idx` is that file's index in the sample's **full sorted list**, so
  the key is shard-count invariant. Two prior bugs this fixes: the sample tag
  was a hardcoded `vbf?0:zjets?1:dijet?2:3` ladder that collapsed all five
  later-registered samples onto id 3 (it now comes from `MyUtl::sampleId()`
  beside the registry, append-only), and `event_num` was the chain-global
  entry number, which restarts at 0 in every shard.

  Supports `--file-shard=<i>/<N>`; **merge the shards with `hadd`, not
  `hist_merge`** — the inverse of the rule for the `*_hist` outputs, because
  this file holds only TTrees (which `hadd` concatenates correctly) and none
  of the `TParameter` scalars that make `hadd` wrong there. Validated on the
  local 33-file VBF sample: a 4-way sharded run, `hadd`'d, is bit-identical to
  the unsharded run in all three trees with row counts summing exactly.
  Submit template: `condor/export_training_data_sharded.sub`.

### z-Association Study
- **assoc_study_hist.cxx** / **assoc_study_plot.cxx** (`util/`): Compares the
  two track-to-vertex z-association rules (z0-significance cut vs. the
  `getNewDzpara` parameterization) by measuring the vertex t0 each produces,
  for both `TRKPTZ` and `WAVES`. Reports inclusive core σ and core fraction,
  plus both as a function of n forward HS tracks and n forward jets, plus
  track-level accounting (kept / HS purity / HS recall). Same hist/plot split
  as the primary analyses; the rule list lives in `assocRules()`
  (`util/assoc_study_common.h`). Full write-up and results:
  `share/docs/Z_ASSOCIATION_STUDY.md`; headline findings in Physics Findings
  below.

  All rules are evaluated against the SAME event in one pass, which is what
  makes the comparison paired. Two structural points make that work:
  - The counting scan stays pinned at `COUNTING_NSIGMA` (3.0) for every rule,
    so `EventCounts` — and therefore every plot x-axis — is common. Otherwise
    the "n_hs_tracks = 3" bin would hold a different set of events per rule.
  - Under a rule the clustering list is **re-selected from the full track
    array**, not narrowed from the counting list: a parameterized rule is not
    a subset of a significance cut (`getNewDzpara` is looser than 3σ for some
    (η, pT) and tighter for others), so narrowing would silently intersect
    every parameterized rule with `significance < 3.0`.

### Utilities
- **maxpt_jet_test.cxx** (`util/`): Diagnostic comparing WAVeS/TRKPTZ pass
  rates under different `maxTrkPt` caps (current `MAX_TRACK_PT`, leading-jet
  pT, subleading-jet pT). Pure diagnostic, changes nothing.
- **track_dt.cxx** (`util/`): Inclusive per-track HGTD-time-vs-truth-particle
  residual, no event selection.
- **runHGTD_Clustering.cxx** (`python/`): Interactive per-event clustering.
  Output is parsed by `event_display.py` to produce a 2-page PDF.

## Histogramming vs. Plotting Split

Both primary analyses follow the FastFrames pattern: run the event loop once
against the ntuples, save every histogram needed for plotting into a ROOT
file, then do all plot styling/derivation/output as a separate, fast, local,
iterable step that never touches the ntuples again. This exists because plot
tweaks (axis range, binning, a new working-point table, a style fix) used to
force a full high-stats rerun (zjets alone is 8M events, multi-minute local /
condor-scale for real work) even though the underlying data hadn't changed.

- **`src/histogram_io.h`**: shared `MyUtl::HistWriter`/`MyUtl::HistReader` —
  a thin RAII wrapper around a ROOT `TFile` for writing/reading a named set
  of histograms plus scalar run metadata. Naming convention: every object
  keeps its own existing ROOT name (`GetName()`) as the TFile key — no
  `TDirectory` nesting — since every histogram booked in this codebase
  already embeds enough context (score/variable/scenario/slice) to be
  globally unique within a run. Scalars use `TParameter<Long64_t>`/
  `TParameter<Double_t>` under a `meta_`-prefixed key (event counts,
  accumulator doubles); the energy label is a `TObjString` under
  `meta_energy_label`.
- **Sample-prefixed output paths** (`src/sample_config.h`'s `histFilePath()`/
  `plotFilePath()`): when a run is given `--sample=<name>`, every output file
  name is prefixed with `<name>_` so files stay distinguishable if copied out
  of their per-sample directory into a shared one. For histograms this also
  flattens away the old `hists/` subdirectory — a single sample-prefixed file
  directly under `OUTPUT_DIR` no longer needs it to stay unambiguous — e.g.
  `vbf/vbf_clustering_hist.root`, `vbf/vbf_rpt_v5_hist.root`. Plots keep their
  existing subdirectory layout (`comparisons/`, `inclusive/`, `fullplots/`,
  `rpt_plots/`, ...); only the file name gets prefixed, e.g.
  `vbf/comparisons/vbf_theplot_pu_frac.pdf`,
  `vbf/rpt_plots/vbf_rpt_v5.pdf`. For the local/default run (no `--sample`),
  both helpers fall back to the original unprefixed conventions
  (`<OUTPUT_DIR>/hists/<name>.root`, `<OUTPUT_DIR>/<subdir>/<name>.pdf`) —
  `MyUtl::SAMPLE_NAME` is only set to a non-empty value by the four primary
  split executables (`clustering_hist`/`clustering_plot`/`rpt_v5_hist`/
  `rpt_v5_plot`) when `resolveSample()` matched a real `--sample=` flag; the
  legacy monolithic executables (`clustering_dt.cxx`, `rpt_v5.cxx`) never set
  it and so keep their original unprefixed, un-flattened output paths
  unconditionally, even when run with `--sample=`.
- **Reconstruction pattern**: the plotting stage always "construct-empty-
  then-load" — it re-runs the *same* booking logic the histogramming stage
  used (`buildAnalysisMap()` / `makeScenarios()`) to get correctly-named,
  correctly-binned empty histograms, then calls `loadFrom()`/`loadScenarios()`
  to fill them from the file. This reuses the existing binning/title/name
  construction code instead of duplicating it in a second place, and doubles
  as a binning-mismatch guard (`HistReader::LoadInto` throws if a saved
  histogram's binning doesn't match the freshly-booked destination's).
- **`AnalysisObj`/`PlotObj::saveTo`/`loadFrom`** (`src/plotting_utilities.h`):
  mirror `mergeFrom`'s exact member list (`WriteHist`/`LoadInto` in place of
  `Add`) — deliberately exclude `params`/`effEstimate`/`slicesHists`/
  `slicesFits`/`slicesModels`/`efficiency`, since those are purely derived by
  `plotPostProcessing()` from the saved raw histograms and must be
  recomputed after `loadFrom()`, not persisted.
- **rpt_v5's ordering constraint, preserved structurally**: `generate_roc()`/
  `rejAtEff()` use explicit `Integral(1, bin+1)` bin ranges specifically so
  they're immune to any prior `SetRangeUser()` restriction; the console
  summary's no-arg `Integral()`/`GetMean()` calls are not immune, so they
  must run before `styleScen()`'s cosmetic `SetRangeUser(0.0, 1.5)`. In
  `rpt_v5_plot.cxx` this is enforced structurally, not just by comment:
  `styleScen` is a local lambda whose *definition* is placed lexically after
  the entire derived-numbers block, so calling it earlier would be a compile
  error (undeclared identifier), not a silent runtime bug.
- **Harmonized ROC window**: the ROC x-axis ("Hard Scatter Efficiency") is a
  fixed `[0.80, 1.0]` window for every sample and pT-slice, so plots are
  directly comparable across vbf/zjets/dijet and against the reference study,
  whose ROCs use the same range. 0.80 is also where `rocEffGrid()`'s first
  sampled working point sits — an earlier `[0.85, 1.0]` zoom cropped that
  point off the drawn window. Tunable by editing
  `roc_xmin_default`/`roc_xmax_default` in `rpt_v5_plot.cxx` and rerunning
  `rpt_v5_plot --sample=<sample>` locally against an already-saved histogram
  file — no condor resubmission.
- **Fixed-efficiency ROC sampling**: `generate_roc_grid()` samples each curve
  at the working points in `rocEffGrid()` (0.800 → 0.975 in 0.025 steps)
  rather than once per histogram bin, so every scenario shares identical
  x-points and the ratio panels are exact. Beware when reading the console
  tables: `rejAtEff` returns the first bin with eff <= target, so near a step
  a 0.02% efficiency difference can move the printed rejection by 30%.
  Repeated identical values across working points mean the curve's ceiling
  has been hit and the sampler is clamping, not that the physics is flat.
- **Per-panel y-maxima**: the inclusive slices, the central baseline, and the
  R1/R2 regions each compute their own. They are far apart in scale — with
  the oracle rows, R1 peaks near 880 and R2 near 1180 against an inclusive
  ceiling near 640 — so a shared maximum clips panels.
- **Selection labels wrap** on an explicit `@@` delimiter supplied by the
  caller. The legend is an opaque box whose left edge is at NDC 0.58 and which
  spans the same y range as the label, so a long single line is painted over.
  Splitting on a width estimate does not work (TLatex markup means glyph count
  != string length) and `|` cannot be the delimiter (it appears in every
  `|#eta|` label).
- **Validated** (bin-for-bin histogram identity + pixel-for-pixel PDF
  identity) against the local default VBF sample for both split pairs before
  ever touching a grid sample — see git history for the validation method
  (temporary `saveTo()` instrumentation in the legacy executables, `uproot`
  bin-content diff, `pdftoppm`+ImageMagick `compare` page-by-page PNG diff).

## WAVeS Clustering / Scoring

`Score::WAVES` (jet-proximity weighted score) is the current alternative to
the baseline `Score::TRKPTZ`, implemented in `Cluster::updateScores`
(`src/clustering_structs.h`) and `Cluster::calculateTime`/`calculatePurity`
(`src/event_processing.h`):

- **Score**: `Σ_i pT_i · pT_jet(i) / max(ΔR_i, WAVES_DR_FLOOR)`, summed over
  the cluster's constituent tracks (nearest qualifying forward jet per
  track), then multiplied by `exp(−1.5·|Δz_cluster|)`. Falls back to the
  TRKPTZ score if no qualifying forward jets exist in the event.
- **Time/purity**: for the `WAVES` row specifically (not its oracle
  variants), both the reported time and purity are recomputed using only the
  subset of the selected cluster's tracks that fall within `ΔR < 0.4` of a
  qualifying forward jet — this purifies the time by dropping absorbed PU
  tracks that the jet-agnostic clustering step picked up.
- **`Score::JET_T_REFINED`** (short name `WAVES_RECLUST`, id 19): a
  *dedicated* collection — tracks are first filtered to those within
  `ΔR < 0.4` of a forward jet (`TrackFilterType::JET`), then re-clustered
  with `ClusteringMethod::ITERATIVE` at the tighter `DIST_CUT_T_REFINED =
  2.0σ`. The best cluster in that collection is selected by the TRKPTZ score
  (via the aux-collection path in `selectClusters`).
- **`Score::WAVES_MISCL`** / **`Score::WAVES_MISAS`**: oracle variants that
  reuse the WAVeS selection and time, gated at fill time on cluster purity
  (`> 0.75`) / event-level HS timing purity (`≥ 0.95`) respectively — see
  `requiresPurity` in the Score Registry section. Unlike the `WAVES` row,
  these use the default full-cluster purity (not the in-jet-only
  recomputation), so their purity gate stays comparable to `TEST_MISAS`.

`Score::CONE` (legacy cone clustering, `ClusteringMethod::CONE` at
`DIST_CUT_CONE = 3.0σ`) still exists in the registry as a historical
baseline but is not currently emplaced in `buildAnalysisMap`.

## Failure Decomposition

`util/failure_decomposition.cxx` assigns every event where the selection
score fails the efficiency window to a **single principal cause**, via a
fixed precedence ladder (first match wins) — this superseded an earlier
non-exclusive four-category scheme:

| Precedence | Name | Condition |
|---|------|-----------|
| 1 | Timing misassignment | Ideal Resolution scenario passes, **or** a majority of HS pT is mis-timed vs. truth (`|pull| ≥ 3σ`) — the latter isn't fixed by ideal resolution since the times are offset, not just smeared |
| 2 | Selection failure | Another qualifying cluster in the same real-time collection passes `passEfficiency()` — the right answer existed but wasn't chosen |
| 3 | Misclustering | Grouping *all* truth-HS tracks into one cluster (real times, no PU contamination) passes the window, but the algorithm's chosen cluster does not |
| 4 | Track efficiency loss | Ideal Res. + Ideal Eff. passes but Ideal Res. alone does not — giving currently time-less HS tracks an ideal time is what recovers the event |
| 5 | Irreducible | None of the above — not an algorithm failure |

## Key Analysis Parameters

Values, with their rationale comments, live in `src/clustering_constants.h`; the notes below are
what those comments do not carry.

**Track-to-vertex association.** Three ways to express it, all in
`clustering_constants.h`:
- `MAX_NSIGMA` — the incumbent z0-significance cut. **Since 2026-11-16 (f897b1e) every z0
  significance uses the RAW z0 covariance** (`Z0_VAR_INFLATION = 1.15²` removed, the user's
  decision): an "N sigma" in any earlier output or note is 1.15 N sigma raw (the talk's 2.5 ->
  2.875, the HS-track count's 3 -> 3.45). See notes/gp_talk.md.
- `USE_DZ_PARA` (runtime, from `--dzpara`) — swaps in `getNewDzpara ×
  DZ0_PARA_SCALE_CLUSTER` instead; the caller's significance argument is then
  ignored, since the two rules are not on a common scale.
- `struct AssocRule` — the same information as an explicit *value* (kind +
  cut), so several rules can be evaluated against one event in a single pass.
  `passTrackVertexAssociation` / `getAssociatedTracks` have AssocRule
  overloads that deliberately do NOT consult the globals, and
  `processEventData` takes an optional `assoc` argument. `nullptr` (the
  default) reproduces the global-driven behaviour exactly, so every executable
  that does not pass one is bit-for-bit unchanged. Only
  `util/assoc_study_hist.cxx` uses it today.
- `AssocRule::orGhost` — unions the rule with ghost association (accept a
  track if it passes the z-test OR is ghost-associated to a qualifying forward
  jet, via `forwardJetGhostTracks`). **This is a SET-level property, not a
  per-track one**: it needs the event's jets, so it is applied inside
  `getAssociatedTracks` and is *not* implemented by
  `passTrackVertexAssociation`. Anything reproducing a rule's track list must
  therefore call `getAssociatedTracks` rather than looping the per-track
  function, or it will silently get the z-only list. (Studied and rejected on
  physics grounds — see Physics Findings.)

The kinematic half of the track selection (η window, pT window, quality flag)
lives in `passTrackKinematics` so the two `getAssociatedTracks` overloads
cannot drift apart — only the association test differs between them.

The z-association study (see below) found the incumbent `MAX_NSIGMA = 3.0` is
**not** the best choice for TRKPTZ. It has not been changed, because the study
has so far only been run on the local VBF sample.

RpT-specific constants (`JET_ETA_MIN/MAX`, `RPT_TRACK_JET_DR`, `DZ0_PARA_SCALE`, `CENTRAL_Z_SIGNIF`,
`GATE_SIGMA`) are at the top of `util/rpt_v5_hist.cxx`; the idealised smears in `src/idealised_timing.h`.

**The two z-association rules are not interchangeable.** `getNewDzpara` comes
from a FORWARD-region reference study and does not extrapolate inward: it
returns sigma_z0 = 31 um at eta = 0 rising to 4.4 mm at eta = 3.8, and the
central end is several times tighter than ITk z0 resolution at 1 GeV. Using it
centrally roughly halves ITk-only rejection (163 -> 77 at 0.800 HS efficiency).
5.0 is not a guess: it is the value the reference macro itself carries as a
commented alternative, and an independent scan of ours reproduces it as the
best working point.

**Vertex-time inflation is per sample AND per scenario**, resolved by
`inflationFor(SAMPLE_NAME)` in `main()` before any worker thread starts. The
table lives in `util/rpt_v5_common.h`, shared with `vbs_time_veto`:

| scenario | vbf | zjets | dijet | ttbar |
|---|---|---|---|---|
| hgtd | 1.48 | 1.61 | 1.53 | (1.53) |
| trkptz | 1.39 | 1.63 | 1.45 | (1.45) |
| waves | 1.38 | 1.65 | 1.46 | (1.46) |
| tzp | 1.38 | 1.65 | 1.46 | (1.46) |

The ttbar column (parenthesised) is SEEDED from dijet, not measured — the
2026-09-01 grid run measured tzp at 1.59/1.43/1.43 vs the seeds (2–4% off,
sub-percent on the gate width), and the dijet-seeded ttbar row validated
across every scenario.

Measured off `PRINT_PULL_DIAG`'s ratio column; the run echoes `[calibration]`
at startup and prints the freshly measured ratio beside the value in use, so a
stale entry is visible in the same log. The measurement is NOT circular (raw
`var_vtx`, ungated track list, vertex times independent of the gate), so one
pass suffices. Re-measure for any NEW sample. See Physics Findings for why
this matters far less than it looks.
`CONE_ITER_K` also still exists as a constant but its only historical
consumer (a `cone_iter_k_sweep.cxx` sweep executable) is not present in this
tree, so it's currently unused.

## Build System

TMVA is intentionally excluded (see `CMakeLists.txt:34`) — `HGTD_SORT`/
`CONE_BDT` are registered `Score`s but their selector
(`chooseHGTDSortCluster`) is currently a stub returning an empty `Cluster`,
and neither is emplaced in `buildAnalysisMap`.

Flags are `-O3 -march=native -g` (no `-DNDEBUG`, no LTO) — see the SIGILL caveat under Grid Submission.

## Data Format

ROOT ntuples, tree `"ntuple"`; the branches are bound in `BranchPointerWrapper` (`src/clustering_structs.h`).

**The grid productions keep ONLY the hard-scatter truth vertex.**
- `TruthVtx_*` has exactly one entry per event on grid VBF and Z+jets, in the
  raw and skimmed copies alike.
- So 98% of tracks there have `Track_truthVtx_idx = -1`: pileup tracks are not
  linked to anything.
- The local VBF sample keeps ~200 truth vertices per event (4% unlinked).
- Test HS as `Track_truthVtx_idx == 0`. Never index `TruthVtx_*` with a raw
  track index without a `>= 0` guard: in Python `[-1]` wraps round to the HS
  vertex. That made `event_display.py` draw every grid pileup track as HS
  until 2026-09-25.

### Skimmed + slimmed copies (`/data/mcardiff/skimmed_ntuples/`)

`util/skim_ntuples.cxx` writes a reduced copy of each sample, once, and every
analysis pass reads that copy rather than the originals.

**This is now the registry default, not an opt-in.** `resolveSample` resolves
`vbf`, `zjets`, `dijet` and `ttbar` to `/data/mcardiff/skimmed_ntuples/<name>/`,
so a bare `--sample=zjets` already reads the skim. `--ntuple-dir=` is the escape
hatch BACK to the originals -- it overrides ONLY the input directory, leaving
the sample's label, output paths, sample name and overlap-removal flag alone:

```bash
./build/rpt_v5_hist --sample=zjets      # skimmed copy (default)
./build/rpt_v5_hist --sample=zjets \
    --ntuple-dir=/data/mcardiff/exotic_superntuples/zjets/   # originals
```

Every run's `[chain] N files added from <dir>` line names the directory it
actually read, so which tree produced a given output is answerable from its log
alone. The raw path for each swapped sample is in a comment beside its registry
entry.

**Four samples are NOT skimmed** and still resolve to the originals: the three
mu=0 entries, plus `zeejets` -- which cannot be skimmed the same way, because it
selects a channel by pointing at ONE DSID subdirectory of the raw zjets tree,
and the skim flattens a sample's whole input tree into one output directory. It
needs its own `--sample=zeejets` skim run, not a path swap.

**`skim_ntuples` refuses to skim a skim.** With the registry pointing at the
copies, a bare `skim_ntuples --sample=zjets` would otherwise read
`skimmed_ntuples/zjets` and write back into it, applying the selection twice.
It now hard-errors on two independent checks -- input and output resolving to
the same directory, and the first input file already carrying `skim_n_input` --
so `condor/skim_ntuples.sub` passes an explicit `--ntuple-dir=$(rawdir)`.
Override `rawdir` for `vbf` (`highstats_vbf`, not `vbf`).

**Consequence for `export_training_data`.** Its group key is
`(sample_id, file_idx, event_num)`, where `file_idx` indexes the sample's full
sorted file list and `event_num` is the entry within that file. Both are defined
against whichever tree was read, and the skim has 48 files where the raw sample
had 1001. An export made entirely from skimmed input is self-consistent; an
export that MIXES pre- and post-switch shards is not, and the collision is
silent. Re-export whole, do not append.

| sample | files | size | vs. original |
|---|---|---|---|
| vbf   | 1001 -> **48** | 448.5 -> **173.6 GB** | 20.9x fewer files, 2.58x smaller |
| zjets | 2000 -> **80** | 612.0 -> **76.7 GB**  | 25.0x fewer files, 7.98x smaller |
| dijet | 1400 -> **60** | 214.0 -> **47.7 GB**  | 23.3x fewer files, 4.49x smaller |
| ttbar | 1000 -> **48** | 855.4 -> **389.3 GB** | 20.8x fewer files, 2.20x smaller |

195 branches -> 58. Z+jets gains most because its Z->ll requirement keeps only
26.5% of events; ttbar gains least because it has no lepton selection, so only
the slim applies.

**The file-count reduction matters more than the size reduction.** The cost of
a pass is dominated by per-file read round trips, not bytes: `/data/mcardiff`
is a symlink to `/cold/user/mcardiff` (CephFS, cold tier), a worker issues
~155 read syscalls per file at ~45 KB each, and that syscall count at ceph
latency IS the 0.3-1.2 s/file in "Job cost" below. Page cache is worth 26x on
the same files (8.1 s/file cold vs 0.31 warm).

**What the skim cuts is the LOOSEST COMMON DENOMINATOR**, and this is the part
to be careful about, since a skimmed-away event is invisible to every
executable forever:

    a reco and a truth vertex exist, AND
    |TruthVtx_z[0] - RecoVtx_z[0]| < MAX_VTX_DZ, AND
    passLeptonSelection()

`passBasicCuts()` ALSO requires `>= MIN_JETS` jets, and putting that in the
skim was measurably wrong: `util/rpt_v5_hist.cxx` deliberately does not apply
it, so the jet cut silently removed 20 of its events and shifted 28 of its
histograms on a 3-file test. It is counted and reported, not cut on. Nothing
topological is applied either -- freezing one study's VBS cuts into every other
study's input is the trap the whole design avoids. **Any new cut added here
must be one that every consumer already applies.**

**Event counts change.** `meta_n_total` after the skim is the skimmed count, so
an "N% of events pass" figure computed against it is a fraction of the skim,
not of the sample. Each output file carries `skim_n_input`,
`skim_n_pass_jets`, `skim_n_pass_vtx`, `skim_n_pass_lepton` and
`skim_n_output` as `TParameter<Long64_t>` so the original denominators travel
with the data.

**The keep list is a contract with `BranchPointerWrapper`.** A branch added to
that constructor must be added to `KEEP_ALWAYS` / `KEEP_LEPTON` /
`KEEP_EXTENDED` in `util/skim_ntuples.cxx`, or the skimmed sample fails at
first access. Branches absent from a given production (`Track_btagIp_*` is
absent from every grid sample) are reported and skipped, matching
`BranchPointerWrapper`'s own `hasBranch()` tolerance.

Validated original-vs-skimmed on the same input files: `vbs_region_diag` all
89 columns bit-identical (row-by-row after sorting), `rpt_v5_hist` all 84
histograms bit-identical with 17 of 19 scalars equal -- the 2 that differ being
`meta_n_total`/`meta_n_events_total`, as above. Regenerate with
`condor/skim_ntuples.sub` (66 jobs, ~1 h; writes straight to `/data` and
transfers nothing back).

### dijet is JZ-sliced, and PV0 is not the hard scatter in the soft slices

`dijet` is five JZ slices in one directory (JZ0 500 files, JZ1 400, JZ2 200,
JZ3 200, JZ4 100). The analysis assumes `RecoVtx_z[0]` is the hard-scatter
vertex, and in the soft slices it simply is not:

| slice | `RecoVtx_z[0]` within 2 mm of truth HS | nearest reco vtx IS index 0 |
|---|---|---|
| JZ0 | 2.9% | 0.6% |
| JZ1 | 38.9% | 35.2% |
| JZ2 | 90.1% | 87.5% |
| JZ3 | 98.4% | 97.9% |
| JZ4 | 99.8% | 99.3% |

File-count weighted that is 46.2%, matching the 46.5% `meta_n_pass_basic` of
the saved rpt_v5 dijet run. So **the dijet result is effectively a JZ2+
measurement** -- `passBasicCuts` discards ~97% of JZ0 and ~61% of JZ1. It is
not that the hard scatter is unreconstructed: SOME reco vertex is within 2 mm
in 94.7% of dijet events, it is just not at index 0.

Two practical consequences. A dijet diagnostic run over the first few files
sees JZ0 only (`setupChain` sorts by path and JZ0 sorts first) and will look
catastrophically broken -- that is how this was found. And `--file-shard`
partitions dijet unevenly in physics content, since a shard's files are spread
across slices with very different survival rates.

## Data Locations (UChicago AF)

**Rule: bulk data lives on `/data`, never in `$HOME`.** `/home` is a shared NFS
export (23 T total, chronically ~80% full) and is quota'd per user; `/data` is
CephFS with petabytes free. Merged products belong on `/data`; per-shard
intermediates are disposable and should be deleted once their merge is
verified. See "Cleaning up after a sharded campaign" in the `condor-grid-submission` skill for
the procedure and the 2026-09-08 precedent (50 G recovered from one directory,
which was half the `$HOME` quota).

| what | where |
|---|---|
| Ntuples the registry reads | `/data/mcardiff/skimmed_ntuples/<sample>/` for vbf/zjets/dijet/ttbar; `/data/mcardiff/exotic_superntuples/<sample>/` for the three mu=0 samples and `zeejets` — see "Skimmed + slimmed copies" above |
| ML exports, current | `/data/mcardiff/training_new/` — 14 files, 19 G; post-vertex-re-export column set |
| ML exports, previous | `/data/mcardiff/training/` — 14 files, 15 G; the pre-re-export column set |

Both ML export generations carry the same 14
`<sample>_{mjj500p0,novbs}_training.root` names, so always say which directory
you mean. `training_new` files are ~1.33x larger than their `training`
counterparts — that ratio is the added vertex-fit/track-quality columns, not
more events.

## Event Display

```bash
cd python && python event_display.py --file_path <FILE.root> --event_num <ENTRY> --extra_time <dt>
```

Internally calls `runHGTD_Clustering.cxx` and parses its stdout to produce a
3-page PDF. Script lives in **`python/`** (not `share/scripts/`), and must be
run from there (the macro is invoked by relative name). `--file_num` no longer
exists; `--file_path` takes the full path of any sample's file.

**Region displays reproduce one rpt_v5 row exactly** (2026-09-25). The
defaults are unchanged, but a display can be given that row's exact gate:
- `--t0_sigma`: σ_t0 before inflation. Without it the display rebuilds σ from
  its own cluster nearest `--extra_time`, which is only right for WAVeS.
- `--infl`, `--gate_sigma`.
- `--assoc_nsigma 2.5`: rpt_v5_hist's clustering input. The display default
  is 3.0.
- `--ideal_times --ideal_t0 truth|cluster`: the macro replays the idealised
  world through `src/idealised_timing.h` and the gate uses those times.
- `--legs HS,PU --leg_labels HS,PU`: the analysis pair, coloured by its paper
  labels. The Python pair search ignores the Z+jets overlap removal, so pass
  the pair explicitly.
- `--expect_rpt`: the display recomputes both legs' R_pT and prints
  `RPT_CHECK OK|MISMATCH`.

`python/rpt_region_displays.py render` passes all of these from the region
tree and fails on any mismatch, so do not hand-assemble region commands.

The VBS candidate pair is marked with a muted steel-blue cross-hatch drawn
UNDER the jet wedges (fill colour is already spoken for: green = truth HS,
grey = PU), and each tagging jet's label carries `R_pT = <z-only> -> <timed>`.
The timed value applies rpt_v5's gate, so it needs `--infl` set to the
sample's WAVeS inflation (vbf 1.38, zjets 1.65, dijet 1.46). t0 comes from
`--extra_time` directly — NOT from the matched cluster's own weighted-mean
time, which differs from WAVeS's in-jet refined time and made the displayed
R_pT disagree with the analysis.

`--jet_idx`/`--jet_label` paint a target jet orange, which OVERRIDES the
truth-identity colouring. Do not pass them for region displays: which leg is
HS and which is PU is the whole point of those plots.

**The local box cannot run this without help.** Its ROOT needs python3.14
while `~/.local/share/python/lib/python/site-packages` (a version-agnostic
user site dir) holds python3.9's numpy, whose compiled extensions fail to load
under 3.14. Set `PYTHONNOUSERSITE=1` to skip it and use `~/.venv-hgtd`, which
has a matching matplotlib. On the AF use ROOT 6.38.04, not an LCG_106 view:
that ships ROOT 6.32, where `TF1::GetParError` has no `const char*` overload
and `src/plotting_utilities.h` fails to compile.

## Machine Learning

No ML is active in the C++ analysis. The history — cluster-selection GBDT / Deep Sets / transformer
(2026-08, concluded), the event-graph GNN (2026-09, exploratory), and the legacy DNN/ONNX artifacts in
`share/models/` — is in the `physics-findings` skill, `results/ml_overview.md` and
`results/gnn_event_graph.md`. Removed study scripts: recover from git history, do not rewrite from scratch.

## Grid Submission (condor)

The full procedure — submit templates, `transfer_output_remaps`, sharding, selection tags, job cost,
cleanup — is in the `condor-grid-submission` skill; load it before any condor work. Rules that must
never be missed:
- **Merge `*_hist` shards with `build/hist_merge`, never `hadd`**: `hadd` keeps one shard's
  `TParameter` scalars beside histograms summed over all shards, silently. TTree-only files (ML
  exports, `*_rpt_v5_regions*.root`) are the opposite: `hadd -f`.
- **Glob on the shard count** (`*.shard*of12.root`), never `shard*`: a directory can hold several
  shardings of the same data, and the loose glob double-counts every event.
- Outputs land in `condor/<sample>/`, not `<repo-root>/<sample>/`; give the plot stage `--hist-file=`.
- The build is `-march=native`: an execute host with a different CPU dies with SIGILL.
- Verify a merge by entry counts, never by file sizes, before deleting shards.

## Score Registry

All scores are defined as `static const Score` members of the `Score` struct in `src/clustering_constants.h` and listed in `SCORE_REGISTRY`. Each carries its own collection spec: `distCut`, `method` (`CONE`/`ITERATIVE`/`SIMULTANEOUS`/`ITERATIVE_SPLIT`), `useZ0`, `filter` (`TrackFilterType`: `ALL`, `JET`, `HS_ONLY`), plus `buildsCollection()` to indicate it manages its own cluster collection (only true when `usesOwnCollection && distCut >= 0`).

The id / short-name / long-name table is `SCORE_REGISTRY` itself. Score 33 (TZP) is still the C++
symbol `Score::TRKPTZ_TZQ`, renamed after the d0-supersedes-qP swap.

Two fields added to `Score` (2026-09-01): **`timeSource`**
(`TimeSource::{FULL, IN_JET, OUT_JET}`) — which subset of the SELECTED
cluster's tracks `calculateTime` builds the reported time from (`FULL`
deliberately not `RAW`: macOS's `<sys/ioctl_compat.h>` defines a `RAW` macro
that only bites when ROOT/cling compiles `runHGTD_Clustering.cxx`, silently
breaking `event_display.py`) — and **`minSubsetTracks`** (fall back to the
full-cluster time below it; 0 = historical behaviour, which the WAVES family
keeps). `calculateTime` dispatches on the flag instead of the old hard-coded
WAVES score list. **Score long names must never contain `/`** — the histogram
key embeds them, and `TFile::Get`/uproot treat `/` as a directory separator:
the write succeeds and every readback fails.

Oracle scores (`requiresPurity=true`) gate their denominator: events outside the gate are not counted. **Important:** `calcHSTimingPurity` (event-level, used by `TEST_MISAS`/`WAVES_MISAS`) returns `0.0f` — not `1.0f` — when no qualifying HS tracks exist, ensuring that events with zero forward HS tracks are excluded from oracle denominators.

`clustering_hist`/`clustering_plot`'s currently active scores (built by `buildAnalysisMap`, `src/clustering_dt_common.h`) are `TRKPTZ`, `WAVES`, `JET_T_REFINED`, `TEST_MISAS`, `WAVES_MISCL`, `WAVES_MISAS`, `VBF_R1/R2`, the full id 24–33 ladder above, plus `HGTD`.

Comparison PDF groups produced by `clustering_plot` from histograms written by `clustering_hist` (`makeComparisonPlots`, one set per `KEYS = {pu_frac, fjet, ftrack, hs_track, truthjets}`):
- `theplot_*.pdf` — HGTD vs TRKPTZ vs TEST_MISAS
- `waves_oracle_*.pdf` — HGTD vs TRKPTZ vs WAVES vs WAVES_MISAS
- `hgtd_trkptz_waves_*.pdf` — HGTD vs TRKPTZ vs WAVES (no oracles)
- `vbf_regions_*.pdf` — WAVES vs VBF_R1 vs VBF_R2 (topology, algorithm fixed)
- `trkptz_tz_*.pdf` — TRKPTZ vs TRKPTZ_TZ vs WAVES (the per-track Δz term isolated)
- `timesource_*.pdf` — TRKPTZ_TZ raw vs in-jet vs out-of-jet time (selection bit-identical)
- `guarded_*.pdf` — the deployable set: TRKPTZ / TZ / TZ_GIJ / **TZP** / WAVES

Plus inclusive resolution/pull plots (`inclusive/inclusivereso_*`,
`inclusive/inclusivepull_*`, log/lin) and shape-comparison plots
(`inclusive/shape_comparison*`) over `{HGTD, TRKPTZ, WAVES, WAVES_MISAS}`
(resolution set) / `{TRKPTZ, WAVES, WAVES_MISAS}` (shape set).

## plotting_utilities.h Notes

- **`moneyPlot` / `plotWithLegend`**: the `plotWithLegend` lambda accepts optional legend coordinates `(legX1, legY1, legX2, legY2)` so each call can position the legend independently (bottom-right for efficiency plots, top-right for resolution plots).
- **TEfficiency rendering**: `TEfficiency` objects are converted to `TGraphAsymmErrors` with x-errors explicitly set to 0 at construction — this avoids ROOT's repaint-time x-error restoration in `PaintGraph`. All efficiency plots are drawn this way.
- **Fit quality guard**: `plotPostProcessing` skips slices with too few entries and only calls `fillEach` when the relative uncertainty on the core σ is `< 50%`. Bins failing these guards are left empty rather than showing garbage fit results.
- **y-axis labels**: efficiency-style plots use **"Core Fraction"**; the pass window is `PASS_SIGMA = 60 ps`.
- **`AtlasStyle.h`**: calls `atlasStyle->SetErrorX(0)` globally (effective for TH1/TGraphErrors; TEfficiency uses the conversion approach above).
- **Histogram directory registration**: `clustering_hist.cxx`, `rpt_v5_hist.cxx` (and the legacy `clustering_dt.cxx`/`rpt_v5.cxx`) call `TH1::AddDirectory(kFALSE)` early in `main()`, before any histogram construction — required since all four parallelize their event loop (see Multicore below): every worker thread's lazily-built histogram set is identically named, and without this call they'd race on `gDirectory`'s registry. `PlotObj`/`AnalysisObj`'s own constructors still don't manage `gDirectory` themselves; they rely on the caller having set this first. The `*_plot` executables don't need this — they never build histograms concurrently, just load into ones already booked single-threaded.
- **`AnalysisObj`/`PlotObj::saveTo`/`loadFrom`** (adjacent to `mergeFrom`/`plotPostProcessing`): round-trip the raw, event-loop-filled histograms through `src/histogram_io.h`'s `HistWriter`/`HistReader`. Mirror `mergeFrom`'s exact member list by construction — if a new histogram member is added to `AnalysisObj`, add it to `mergeFrom`, `saveTo`, and `loadFrom` together. See "Histogramming vs. Plotting Split" above for the full design.

## Physics Findings

> ### ⚠ Selection change, 2026-08-26
> The default VBS selection changed from **|Deta| >= 3, m_jj >= 0** to
> **no |Deta|, m_jj >= 200 GeV** (`VBS_JET_D_ETA_DEFAULT` /
> `VBS_JET_MJJ_DEFAULT`). |Deta| was measured purity-*negative* on Z+jets
> (2.69% -> 2.51%) and dijet (40.1% -> 35.0%) and neutral on VBF, while m_jj
> is at-least-as-good everywhere — see the constants' comments for the
> `vbs_mjj_diag` numbers behind the swap.
>
> **Every result below this banner was produced under the OLD default**, except
> where a section states its own selection. They are not directly comparable
> with anything produced after the change, and the z-association study in
> particular (`share/docs/Z_ASSOCIATION_STUDY.md`) is entirely pre-change.
> Re-run before quoting a number against new output.
>
> Reproducing a pre-change result now needs `--vbs-deta=3 --vbs-mjj=0`, which
> tags its output `deta3p0_mjj0p0_`. Histogram files written from 2026-08-26
> carry `meta_vbs_deta`/`meta_vbs_mjj`, and the `*_plot` stages abort on a
> mismatch (`HistReader::CheckSelection`); files older than that cannot be
> checked and produce a warning instead — those are exactly the ones whose
> untagged name no longer implies their selection.

### TZP: the classical selector (2026-09-01) — supersedes WAVeS as the best non-ML score

`Score::TRKPTZ_TZQ` (short name `TRKPTZ_TZP`, id 33). One closed-form,
sample-independent algorithm:

```
S(C)  =  [ Σ_t pT_t · exp(−1.0·|z0_t − z_PV|) ] · exp(−0.6·|z_C − z_PV|) · [ Σ_t 1/σ²_d0,t ]^0.225
t0    =  in-jet inverse-variance mean if the cluster holds ≥3 in-jet tracks, else full-cluster mean
```

Core fraction (novbs; mjj500 in parens; μ=0 controls −0.01 everywhere over 4.9M events):

| | vbf | zjets | dijet | ttbar |
|---|---|---|---|---|
| TZP | 91.66 (92.08) | 64.71 (**64.31**) | 88.97 (88.58) | 89.08 (88.84) |
| vs TRKPTZ | +1.6 | **+2.4** | +1.7 | +1.8 |
| vs deployed WAVeS | +0.2 | **+6.4** | +2.5 | +2.8 |

Canonical zjets 61.87 → **64.31** with no ML and no per-sample configuration.
Consumer side (rpt_v5, matched-efficiency reading — see the caveats in the
rpt_v5 section): VBF >40 GeV rejection +3–4% over deployed WAVeS at 0.90;
ttbar a wash among real scenarios but above the truth row; zjets RpT stays
~10–15% behind WAVeS at matched eff despite the +6.4 core-fraction margin —
the core-fraction/RpT dissociation above.
Mechanisms and the full results trail are in the `physics-findings` skill.

- **Closed negatives** (do not retry): canonical vertex-WAVeS weighting (its
  score is ~10⁴:1 concentrated on the PV, so it cancels in the within-event
  argmax); jet-fraction score term (population coin-flip on zjets, 0.612 AUC vs
  0.90 on the loss subset); aggregate-t0 as score term or timing fallback
  (aggregate right only ~27% of disagreements); multiplicity-gated chooser
  (truth gate no better than reco gate — multiplicity is not the variable);
  per-event adaptive α (optimum is 0.7–1.0 in every n-vertices bin).
- **What remains uncaptured**: union(TZP, aggregate t0) is 3–7 points above TZP
  on every sample; only a chooser reaches it and no classical one works. The
  `|t_in − t_out|` split-clustering agreement flag is the strongest reco-only
  confidence variable found (79% vs 58% core at +19–27 in every multiplicity
  bin) — unused by TZP, natural consumer is the RpT gate.
- **Methodology, hard-won**: (1) rank variables on the events a method LOSES
  and they overstate themselves — generate hypotheses there, but believe only
  population-wide switch accounting; (2) any best-of-N oracle needs a
  matched-noise null (the 3-timing oracle read +1.12 real, null +1.67);
  (3) event-level quantities cancel exactly in a within-event argmax;
  (4) `truth_n_hs_tracks` is PER CLUSTER — event totals are the SUM, and
  `.first()` silently reads one cluster's count.

### Detailed findings (lazy-loaded)

The VBS pair composition / JVT study is in the `vbs-timing-studies` skill; the loss decomposition,
topology gain, R2 capture, gate severity, efficiency ceiling, calibration, z-association study and
reference replication are in the `physics-findings` skill.

**Anything setting `EXTENDED_BRANCHES` must call `recordAvailableBranches`
(src/event_processing.h) after `setupChain` and before constructing the
wrapper.** The grid samples lack `Track_btagIp_*`; a TTreeReader with a
missing branch iterates ZERO entries and exits 0. The first Z+jets submission
of this diagnostic produced three "successful" empty files that way.

### Methodological warnings

- **Never rank a selection oracle by cluster purity.** `calcPurity` is a
  pT-weighted fraction of tracks from truth vertex 0 — it measures where a
  track CAME FROM, not whether its time is right. A lone HS track carrying a
  mis-assigned HGTD hit forms a 100%-pure cluster sitting at a wrong time, so
  max-purity selects exactly the dominant failure mode (HS jets driven to
  R_pT = 0 rise 3.0% -> 11.8%). Rank on `|t_cluster - t_truth|` instead.
- **When selecting display candidates for a failure mode, require the
  failure.** Gating on low HS timing purity alone yields events where the
  tracks are badly mis-timed and the algorithm recovers the right time anyway.
- **A flat 30 ps track smear is WORSE than HGTD's real timing** (Initial
  scenario: 35/25/21 ps for 1/2/3 hits, ~26 ps effective), so the reference's
  30 ps scenario is not an upper bound on what the simulated detector
  delivers.

## Development Notes

- **Namespace**: all utilities live in `namespace MyUtl`
- **Debug mode**: `MyUtl::DEBUG` (`src/clustering_constants.h`) for verbose per-event logging
- **Multicore**: `clustering_hist.cxx`/`rpt_v5_hist.cxx` (and the legacy `clustering_dt.cxx`/`rpt_v5.cxx`) parallelize their per-event loop via `ROOT::TTreeProcessorMT` (`--threads=N`, `MyUtl::resolveThreads()` in `src/sample_config.h`, default `min(hardware_concurrency(), 8)`) — `ROOT::EnableImplicitMT()` alone only parallelises ROOT-internal operations, not the event loop itself. Per-thread accumulation state (`AnalysisObj` map / RpT `ThreadState`) is built lazily once per worker thread and merged into one result after the loop (`AnalysisObj::mergeFrom()` in `src/plotting_utilities.h`). See `util/ttreemt_prototype.cxx` for the standalone prototype this was validated against, and the "TTreeProcessorMT prototype" note below on a real thread-safety bug found via stress-testing. The `*_plot` executables are single-threaded — no event loop, no ntuple access, just histogram load + fit/derive + draw.
- **Constants**: modify analysis cuts in `src/clustering_constants.h`; always rebuild after header changes
- **Plot styling**: ROOT `TColor`/`TCanvas`; ATLAS style via `AtlasStyle.h`/`AtlasLabels.h`
- **Event loop**: `TTreeProcessorMT::Process()` hands each task its own `TTreeReader`; a fresh `BranchPointerWrapper` is bound to it per invocation (cheap — pointer setup only)
- **File identification**: inside the `TTreeProcessorMT` lambda, `reader.GetTree()->GetCurrentFile()->GetName()` gives the current file path, and `reader.GetTree()->GetReadEntry()` gives the local (0-indexed) entry number directly — both needed to reconstruct an event display command. This differs from the old sequential loop's `chain.GetReadEntry() - chain.GetChainOffset()` trick, which needed an outer `TChain` not available inside the lambda; a `TChain`-bound *sequential* `TTreeReader`'s `GetTree()->GetReadEntry()` returns the chain-global entry number, not local, which is exactly why that trick existed.
- **TTreeProcessorMT thread-safety gotcha**: `AnalysisObj`'s constructor calls `SetFillColorAlpha()`, which lazily registers a new `TColor` into ROOT's *global* color table the first time a given (color, alpha) pair is used — not thread-safe. Worker threads racing to build their first `AnalysisObj`/`ThreadState` concurrently corrupted the heap (intermittent `SIGABRT`, found via stress-testing, not a one-off). Fixed by serializing each thread's entire lazy-init construction (not just registry bookkeeping) behind its mutex — see the per-thread registry pattern in both executables' `main()`.
