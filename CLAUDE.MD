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

```
vertex_timing/
├── src/                           # Core library headers + main analysis executables
│   ├── clustering_constants.h    # Analysis constants, Score struct + SCORE_REGISTRY, ClusteringMethod enum
│   ├── clustering_functions.h    # Clustering algorithms (clusterTracksInTime, chooseCluster, ...)
│   ├── clustering_includes.h     # Common ROOT/C++ includes
│   ├── clustering_structs.h      # BranchPointerWrapper, Cluster (incl. updateScores/calcFeatures), EventCounts
│   ├── sample_config.h           # --sample=vbf|zjets|dijet + --hist-file= resolution (default: local ../../ntuple-hgtd/)
│   ├── histogram_io.h            # HistWriter/HistReader: TFile-based histogram+scalar I/O shared by the
│   │                             #   *_hist/*_plot split executables (see "Main Executables" below)
│   ├── event_processing.h        # setupChain, getAssociatedTracks, selectClusters, processEventData
│   ├── plotting_utilities.h      # PlotObj / AnalysisObj + all plotting helpers (moneyPlot, inclusivePlot, ...)
│   ├── AtlasLabels.h / AtlasStyle.h  # ATLAS plot style utilities
│   ├── clustering_dt_common.h    # Shared by clustering_hist.cxx/clustering_plot.cxx: Scenario enum,
│   │                             #   buildAnalysisMap, makeComparisonPlots
│   ├── clustering_hist.cxx       # Histogramming stage: event loop only, writes hists/clustering_hist.root
│   ├── clustering_plot.cxx       # Plotting stage: reads that file, produces all PDFs, no ntuple access
│   └── clustering_dt.cxx         # Legacy monolithic executable (event loop + plotting in one run); retained
│                                 #   only as the validation ground truth for the split above, see Main Executables
├── util/                          # Diagnostic / report / ML-data-export executables
│   ├── failure_decomposition.cxx # Per-event principal-cause failure classification (see section below)
│   ├── export_training_data.cxx  # Exports per-cluster CSV rows (8 features + purity/labels) for future DNN retraining
│   ├── generate_rpt.cxx          # RpT-style timing/discrimination report (older version)
│   ├── rpt_v2.cxx .. rpt_v4.cxx  # RpT discrimination study iterations
│   ├── rpt_v5_common.h           # Shared by rpt_v5_hist.cxx/rpt_v5_plot.cxx: Scenario struct, RpT histogram
│   │                             #   binning/booking, saveScenarios/loadScenarios
│   ├── rpt_v5_hist.cxx           # Histogramming stage: event loop only, writes hists/rpt_v5_hist.root
│   ├── rpt_v5_plot.cxx           # Plotting stage: ROC/console-summary/PDF from that file, no ntuple access
│   ├── rpt_v5.cxx                # Legacy monolithic executable; retained only as validation ground truth
│   ├── maxpt_jet_test.cxx        # A/B/C test of the MAX_TRACK_PT cap vs per-event leading-jet pT
│   ├── track_dt.cxx              # Inclusive per-track HGTD-time-vs-truth-particle-time residual
│   └── sweep_utilities.h         # Shared helpers for parameter-sweep executables (no sweep .cxx files currently present)
├── python/                        # Interactive/visualisation scripts (run from within this dir)
│   ├── event_display.py          # Event visualisation — calls runHGTD_Clustering.cxx, produces 2-page PDF
│   ├── runHGTD_Clustering.cxx    # Interactive per-event clustering
│   ├── clustering_animation.py   # Animates the iterative clustering algorithm for one event
│   └── event_tinkering.py        # Quick event inspection and jet printing via uproot
├── share/
│   ├── models/                   # ML artifacts from a prior (currently inactive) DNN/BDT effort — see
│   │                             #   Machine Learning section below
│   ├── scripts/                  # ONNX export/inspection Python utilities (historical; no active consumer in src/)
│   └── docs/                     # Design docs, incl. older ML integration guides
├── condor/                        # HTCondor grid submission templates (see Grid Submission)
│   ├── run_analysis.sh           # Generic job executable (any --sample/--threads-aware target)
│   ├── clustering_dt.sub         # Submit file: one job per sample, for the legacy clustering_dt
│   ├── clustering_hist.sub       # Submit file: one job per sample, for clustering_hist (histogramming stage)
│   ├── rpt_v5.sub                # Submit file: one job per sample, for the legacy rpt_v5
│   ├── rpt_v5_hist.sub           # Submit file: one job per sample, for rpt_v5_hist (histogramming stage)
│   └── logs/                     # condor stdout/stderr/log output
├── CMakeLists.txt
└── build/                        # Build output (created by `cmake ..`)
```

Note: several directories/files described in earlier versions of this doc
(`vertex_time.cxx`, `check_jet_sorting.cxx`, `test_vbf_*.cxx`, `DBSCAN.ipynb`,
`figs/`, `old/`, `src/ml_model.h`, and most of the `util/` diagnostics and
parameter-sweep executables) do not exist in the current tree — this section
reflects what's actually present.

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
  and writes `<OUTPUT_DIR>/hists/clustering_hist.root`. `clustering_plot`
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
  version) and writes `<OUTPUT_DIR>/hists/rpt_v5_hist.root` (5 scenarios × 2
  pT slices × {HS, PU} histograms, plus event-count/floor-counter scalars).
  `rpt_v5_plot` reads that file and does `generate_roc`, `rejAtEff`,
  `pickXRange` (including the `MIN_HS_EFF = 0.5` ROC-window floor — see
  below), the console summary tables, and all PDF drawing — no ntuple
  access. Shared `Scenario`/binning/save-load logic lives in
  `util/rpt_v5_common.h`.
- **generate_rpt.cxx**, **rpt_v2.cxx**..**rpt_v4.cxx** (`util/`): Earlier
  iterations of the same study (z-only, real-HGTD, TRKPTZ-selected timing
  gates), monolithic (event loop + plotting), not split.
- **rpt_v5.cxx** (`util/`): Legacy monolithic executable that
  `rpt_v5_hist`/`rpt_v5_plot` was validated against; still builds, but the
  split pair is the primary workflow going forward.

### ML Data Export
- **export_training_data.cxx** (`util/`): Exports one CSV row per
  cluster/passing event (8 features + purity/labels) mirroring
  `clustering_dt`'s event/track selection exactly. Produces training data for
  a possible future DNN, but no model currently consumes it in `src/` — see
  Machine Learning Model section.

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
- **ROC-window floor fix**: `pickXRange()`'s populated-but-shifted-low
  fallback previously had no lower bound — for Z+jets (many forward "HS"
  jets with ~no tracks) it could widen down past 50% HS efficiency, a
  physically meaningless "high-efficiency working point". `rpt_v5_plot.cxx`
  now clamps that fallback's lower edge to `MIN_HS_EFF = 0.5`. Tunable purely
  by editing the constant and rerunning `rpt_v5_plot --sample=zjets` locally
  against an already-saved `rpt_v5_hist.root` — no condor resubmission.
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

From `src/clustering_constants.h`:
```cpp
MIN_TRACK_PT          = 1.0   // GeV — track selection floor
MAX_TRACK_PT          = 30.0  // GeV — track selection ceiling
DIST_CUT_CONE         = 3.0   // σ — main iterative clustering distance (also CONE baseline)
DIST_CUT_T_REFINED    = 2.0   // σ — JET_T_REFINED (WAVES_RECLUST) stage-2 iterative distance
MIN_CLUSTER_TRACKS    = 0     // minimum tracks required to select a cluster
IDEAL_TRACK_RES       = 1.0   // ps — flat per-track resolution for ideal scenarios
PASS_SIGMA            = 60.0  // ps — ±window for "Core Fraction" efficiency oracle
TRUTH_PULL_CUT        = 3.0   // σ — |pull| threshold for oracle timing purity gates
WAVES_DR_FLOOR        = 0.05  // minimum ΔR for the WAVeS 1/ΔR weight (prevents divergence)
```
`CONE_ITER_K` also still exists as a constant but its only historical
consumer (a `cone_iter_k_sweep.cxx` sweep executable) is not present in this
tree, so it's currently unused.

## Build System

### Requirements
- CMake ≥ 3.10
- C++17 compiler
- ROOT (RIO, Core, Gpad, Graf, Hist, Tree)
- Boost (filesystem)

TMVA is intentionally excluded (see `CMakeLists.txt:34`) — `HGTD_SORT`/
`CONE_BDT` are registered `Score`s but their selector
(`chooseHGTDSortCluster`) is currently a stub returning an empty `Cluster`,
and neither is emplaced in `buildAnalysisMap`.

### Build Commands
```bash
cd build
cmake ..
make [<target>]   # e.g. make failure_decomposition
```

### Compiler Flags
`-O3 -march=native -g` (no `-DNDEBUG`, no LTO — `CMakeLists.txt:37`).

## Data Format

ROOT ntuples (tree name `"ntuple"`) with the branches bound in
`BranchPointerWrapper` (`src/clustering_structs.h`):

| Group | Branches |
|-------|----------|
| Tracks | `Track_z0/d0/pt/eta/phi/theta/qOverP`, `Track_var_z0/d0/qOverP/theta/phi0`, `Track_time/timeRes/hasValidTime`, `Track_quality`, `Track_truthVtx_idx/truthPart_idx`, `Track_nHGTDHits/nHGTDPrimaryHits` |
| Vertices | `TruthVtx_z/time/isHS`, `RecoVtx_z/time/timeRes/hasValidTime` |
| Reco jets | `AntiKt4EMTopoJets_pt/eta/phi`, `AntiKt4EMTopoJets_truthHSJet_idx`, `AntiKt4EMTopoJets_ghostTrack_idx` |
| Truth jets | `TruthHSJet_pt/eta/phi`, `TruthITPUJet_pt/eta/phi`, `TruthOOTPUJet_pt/eta/phi` |
| Particles | `TruthPart_prodVtx_time` (bound as `particleT`) |
| Event | `weight` |

## Event Display

```bash
cd python && python event_display.py --file_num <NUM> --event_num <ENTRY> --extra_time <dt>
```

Internally calls `runHGTD_Clustering.cxx` and parses its stdout to produce a 2-page PDF (R-Z geometry view + time histogram). Script lives in **`python/`** (not `share/scripts/`).

## Machine Learning Model

**Currently inactive.** There is no `src/ml_model.h`, no ONNX runtime
integration, and no `Score::TEST_ML` in the current `Score` registry — the
WAVeS jet-proximity score (see section above) is the current alternative to
the TRKPTZ baseline. What remains from the prior ML effort:
- `share/models/` still holds artifacts (`model.onnx`, `model_weights.json`,
  `normalization_params.json`, a TMVA `.weights.xml`) from that earlier
  work, but nothing in `src/` currently loads or runs inference against
  them.
- `share/scripts/export_onnx_model.py` (SavedModel → ONNX → JSON weights
  pipeline) and `share/docs/ML_MODEL_INTEGRATION.md` /
  `ONNX_INTEGRATION_GUIDE.md` document that pipeline for reference, but it
  has no active consumer.
- `util/export_training_data.cxx` still exports per-cluster training CSVs
  (8 features, purity/labels) mirroring the current `clustering_dt`
  selection, in case retraining is picked back up later.
- `src/clustering_structs.h`'s `Cluster::calcFeatures()` still computes and
  normalises the same 8-feature vector (used historically as DNN input); its
  `rawDeltaZ`/`rawDeltaZResunits` outputs are reused by `updateScores` for
  TRKPTZ, but the normalised `features` vector itself is otherwise unused
  today.

## Common Workflows

```bash
# Full analysis (two-step: histogram once, then plot — see "Histogramming
# vs. Plotting Split" above)
cd build && ./clustering_hist && ./clustering_plot

# Iterate on plotting only (styling, binning, a new comparison group) without
# re-running the event loop, as long as ../figs/hists/clustering_hist.root
# is already there:
cd build && ./clustering_plot

# Failure decomposition
cd build && ./failure_decomposition

# Event display for a specific event
cd python && python event_display.py --file_num 000009 --event_num 1856 --extra_time 0.0

# RpT study (WAVeS-based), two-step
cd build && ./rpt_v5_hist && ./rpt_v5_plot

# Iterate on rpt_v5's plotting only (e.g. tuning pickXRange's MIN_HS_EFF)
cd build && ./rpt_v5_plot

# Rebuild after header change
cd build && make
```

## Grid Submission (condor)

Templates in `condor/` run `clustering_hist` / `rpt_v5_hist` on the UChicago AF HTCondor pool, one job per `--sample` value (`vbf`, `zjets`, `dijet`). Condor only runs the **histogramming** stage — the plotting stage is fast/local and runs afterward on whatever machine the `.root` file lands on, with no condor submission needed:

- **`condor/run_analysis.sh`** — generic executable, invoked as `run_analysis.sh <executable> <sample> [<threads>]`. Sources `atlasLocalSetup.sh` + `lsetup "root <version>"` (relaxing `set -e`/`-u` around just that call — those scripts reference unset internal variables like `ALRB_frontlineSite` and aren't strict-mode safe), then runs the transferred `<executable>` **directly from the job's scratch directory** with `--sample=<sample>` and, if given, `--threads=<N>`. No longer moves the binary into a fake `build/` subdir first (see `src/sample_config.h`'s output-path note below) — that indirection existed only to make the old `../<sample>` output convention resolve to the scratch root; the grid-sample output paths are cwd-relative now, so it's unnecessary and was a latent risk (if that step were ever skipped, `../<sample>` would climb out of the scratch directory into its parent, which may not be writable). Unchanged by the histogramming/plotting split otherwise — generic enough for both the legacy and `*_hist` executables.
- **`condor/clustering_hist.sub`** / **`condor/rpt_v5_hist.sub`** — one `.sub` file per histogramming executable, each queuing 3 jobs via `queue sample in (vbf, zjets, dijet)`. Submit each separately (own `condor_submit` call) to get independent job clusters. `arguments` passes `4` as the `<threads>` arg, matching `request_cpus = 4` — both `*_hist` executables parallelize their event loop via `TTreeProcessorMT` (see Multicore in Development Notes), and condor's cgroup throttles any threads spawned beyond `request_cpus`, so the two must track each other.
- **`condor/clustering_dt.sub`** / **`condor/rpt_v5.sub`** — legacy templates for the monolithic executables; still work, kept for reference/fallback, not the primary workflow.

Execute hosts do **not** share a filesystem with the submit host on this pool — confirmed via a `TARGET.FileSystemDomain == MY.FileSystemDomain` match failure when `should_transfer_files = NO` was tried first (0 of 157 slots matched). The `*_hist.sub` files use `transfer_input_files` (stages the prebuilt binary in) and `transfer_output_files = $(sample)/hists/<name>.root` (stages back a single small histogram file — much smaller than the full PDF tree the legacy `.sub` files' `transfer_output_files = $(sample)` transfers). Ntuple inputs are read via the absolute `/data/mcardiff/exotic_superntuples/...` path baked into `--sample=` (see `src/sample_config.h`), so they don't need transferring.

**`preserve_relative_paths = True` is required** on both `*_hist.sub` files — condor's default output-transfer behavior flattens any nested `transfer_output_files` path to its basename when staging back, it does not recreate `$(sample)/hists/` in the submit directory on its own. Without it, all 3 queued jobs (one per `--sample`) each try to write back a file literally named `clustering_hist.root`/`rpt_v5_hist.root` directly into the submit directory, silently overwriting each other — confirmed by running the real jobs once without this flag: only the last-to-finish sample's histogram survived, and no `vbf/`, `zjets/`, `dijet/` subdirectories were created at all. (The legacy `.sub` files' `transfer_output_files = $(sample)` isn't affected — that argument has no nested subdirectory, so its basename already equals the whole thing they want back.)

```bash
cd build && cmake .. && make clustering_hist rpt_v5_hist   # build before submitting; condor doesn't rebuild
cd condor
condor_submit -dump dryrun.ad clustering_hist.sub          # dry-run: validate without queuing
condor_submit clustering_hist.sub
condor_submit rpt_v5_hist.sub                               # separate job cluster
condor_q
# After jobs complete, the .root files land in ../<sample>/hists/. Plot
# locally — fast, no further condor submission:
cd ../build
./clustering_plot --sample=<sample>
./rpt_v5_plot --sample=<sample>
```

Logs land in `condor/logs/clustering_hist_<sample>.<ClusterId>.<ProcId>.{out,err,log}` (prefixed `rpt_v5_hist_` for the rpt_v5_hist cluster).

**Caveat:** the build uses `-march=native` (see Compiler Flags above), tuned to whichever machine ran `make`. If execute hosts have a different CPU microarchitecture than the build host, jobs can crash with an "Illegal instruction" (SIGILL) rather than a normal error.

## Score Registry

All scores are defined as `static const Score` members of the `Score` struct in `src/clustering_constants.h` and listed in `SCORE_REGISTRY`. Each carries its own collection spec: `distCut`, `method` (`CONE`/`ITERATIVE`/`SIMULTANEOUS`/`ITERATIVE_SPLIT`), `useZ0`, `filter` (`TrackFilterType`: `ALL`, `JET`, `HS_ONLY`), plus `buildsCollection()` to indicate it manages its own cluster collection (only true when `usesOwnCollection && distCut >= 0`).

| id | Short | Long name | Notes |
|----|-------|-----------|-------|
| 0  | HGTD  | HGTD Algorithm | Own collection; selected by proximity to `RecoVtx_time` |
| 1  | TRKPT | Σ pT | Plain pT sum |
| 2  | TRKPTZ | Σ pTe^{−1.5\|Δz\|} | Primary/baseline score |
| 3  | PASS  | Pass Cluster | First cluster passing quality cuts (not currently emplaced in `buildAnalysisMap`) |
| 7  | CONE  | Cone | Legacy cone baseline; own collection at `DIST_CUT_CONE`, `ClusteringMethod::CONE` (not currently emplaced) |
| 8  | CONE_BDT | Cone (BDT) | Threshold 0.3; selector is a stub (TMVA excluded) — currently non-functional |
| 9  | FILTJET | Filter Tracks in Jets | JET filter; own collection at `DIST_CUT_CONE`, `ClusteringMethod::CONE` (not currently emplaced) |
| 10 | HGTD_SORT | HGTD BDT (pT-sorted) | Threshold 0.3; selector is a stub (TMVA excluded) — currently non-functional |
| 13 | TEST_MISAS | ΣpTe^{−1.5\|Δz\|} [Events with Perfect Timing] | Oracle: `requiresPurity`, gated on event-level HS timing purity ≥ 0.95 |
| 14 | TEST_HS | ΣpTe^{−1.5\|Δz\|} (HS only) | Own collection; `HS_ONLY` filter, `ClusteringMethod::CONE` at `DIST_CUT_CONE` (not currently emplaced) |
| 18 | WAVES | WAVeS Score | Jet-proximity weighted score; main collection (see WAVeS section above) |
| 19 | JET_T_REFINED | WAVeS 2σ t Re-clustering | Own collection (`WAVES_RECLUST` short name); `JET` filter, `ClusteringMethod::ITERATIVE` at `DIST_CUT_T_REFINED`; selected by TRKPTZ |
| 20 | WAVES_MISCL | WAVeS [Events with Pure Clusters] | Oracle: `requiresPurity`, gated on cluster purity > 0.75 |
| 21 | WAVES_MISAS | WAVeS [Events with Perfect Timing] | Oracle: `requiresPurity`, gated on event-level HS timing purity ≥ 0.95 |

Oracle scores (`requiresPurity=true`) gate their denominator: events outside the gate are not counted. **Important:** `calcHSTimingPurity` (event-level, used by `TEST_MISAS`/`WAVES_MISAS`) returns `0.0f` — not `1.0f` — when no qualifying HS tracks exist, ensuring that events with zero forward HS tracks are excluded from oracle denominators.

`clustering_hist`/`clustering_plot`'s currently active scores (built by `buildAnalysisMap`, `src/clustering_dt_common.h`) are `TRKPTZ`, `WAVES`, `JET_T_REFINED`, `TEST_MISAS`, `WAVES_MISCL`, `WAVES_MISAS`, plus `HGTD`.

Comparison PDF groups produced by `clustering_plot` from histograms written by `clustering_hist` (`makeComparisonPlots`, one set per `KEYS = {pu_frac, fjet, ftrack, hs_track, truthjets}`):
- `theplot_*.pdf` — HGTD vs TRKPTZ vs TEST_MISAS
- `waves_oracle_*.pdf` — HGTD vs TRKPTZ vs WAVES vs WAVES_MISAS
- `hgtd_trkptz_waves_*.pdf` — HGTD vs TRKPTZ vs WAVES (no oracles)

Plus inclusive resolution/pull plots (`inclusive/inclusivereso_*`,
`inclusive/inclusivepull_*`, log/lin) and shape-comparison plots
(`inclusive/shape_comparison*`) over `{HGTD, TRKPTZ, WAVES, WAVES_MISAS}`
(resolution set) / `{TRKPTZ, WAVES, WAVES_MISAS}` (shape set).

## event_processing.h Internals

Key structures and entry points:

- **`EventResult`**:
  ```cpp
  struct EventResult {
    int    code          = -1;   // -1: rejected by selection; 0: normal; 3: MISAS fail
    double time          = -1.0; // TRKPTZ-selected cluster time
    int    nFwdHS        =  0;   // n forward HS tracks (3σ counting step)
    bool   trkptzPass    = false;
    bool   misasInDenom  = false; // event entered the TEST_MISAS denominator
    bool   misasPass     = false; // in denominator AND passed
    double clusQuality   =  0.0;  // combined quality of the TRKPTZ-selected cluster
  };
  ```
  `processEventData` returns `EventResult`; early-exit returns `{}` (code=-1).

- **Section E** (auxiliary collections): an `applyFilter` lambda dispatches on `TrackFilterType`, then a `SCORE_REGISTRY` loop calls `clusterTracksInTime` for every score where `s.buildsCollection()` is true, populating `unordered_map<int, vector<Cluster>> auxCollections`. No named per-score collection variables.

- **`selectClusters`**: takes the main clusters, tracks, branch, the active `analyses` map, and `auxCollections`; returns `unordered_map<int,Cluster>` keyed on `Score::id`. Handles the main-collection scores via the all-scores `chooseCluster` overload, aux-collection scores via TRKPTZ selection, and `HGTD`/`HGTD_SORT`/`CONE_BDT` via their own dedicated logic.

## plotting_utilities.h Notes

- **`moneyPlot` / `plotWithLegend`**: the `plotWithLegend` lambda accepts optional legend coordinates `(legX1, legY1, legX2, legY2)` so each call can position the legend independently (bottom-right for efficiency plots, top-right for resolution plots).
- **TEfficiency rendering**: `TEfficiency` objects are converted to `TGraphAsymmErrors` with x-errors explicitly set to 0 at construction — this avoids ROOT's repaint-time x-error restoration in `PaintGraph`. All efficiency plots are drawn this way.
- **Fit quality guard**: `plotPostProcessing` skips slices with too few entries and only calls `fillEach` when the relative uncertainty on the core σ is `< 50%`. Bins failing these guards are left empty rather than showing garbage fit results.
- **y-axis labels**: efficiency-style plots use **"Core Fraction"**; the pass window is `PASS_SIGMA = 60 ps`.
- **`AtlasStyle.h`**: calls `atlasStyle->SetErrorX(0)` globally (effective for TH1/TGraphErrors; TEfficiency uses the conversion approach above).
- **Histogram directory registration**: `clustering_hist.cxx`, `rpt_v5_hist.cxx` (and the legacy `clustering_dt.cxx`/`rpt_v5.cxx`) call `TH1::AddDirectory(kFALSE)` early in `main()`, before any histogram construction — required since all four parallelize their event loop (see Multicore below): every worker thread's lazily-built histogram set is identically named, and without this call they'd race on `gDirectory`'s registry. `PlotObj`/`AnalysisObj`'s own constructors still don't manage `gDirectory` themselves; they rely on the caller having set this first. The `*_plot` executables don't need this — they never build histograms concurrently, just load into ones already booked single-threaded.
- **`AnalysisObj`/`PlotObj::saveTo`/`loadFrom`** (adjacent to `mergeFrom`/`plotPostProcessing`): round-trip the raw, event-loop-filled histograms through `src/histogram_io.h`'s `HistWriter`/`HistReader`. Mirror `mergeFrom`'s exact member list by construction — if a new histogram member is added to `AnalysisObj`, add it to `mergeFrom`, `saveTo`, and `loadFrom` together. See "Histogramming vs. Plotting Split" above for the full design.

## Development Notes

- **Namespace**: all utilities live in `namespace MyUtl`
- **Debug mode**: `MyUtl::DEBUG` (`src/clustering_constants.h`) for verbose per-event logging
- **Multicore**: `clustering_hist.cxx`/`rpt_v5_hist.cxx` (and the legacy `clustering_dt.cxx`/`rpt_v5.cxx`) parallelize their per-event loop via `ROOT::TTreeProcessorMT` (`--threads=N`, `MyUtl::resolveThreads()` in `src/sample_config.h`, default `min(hardware_concurrency(), 8)`) — `ROOT::EnableImplicitMT()` alone only parallelises ROOT-internal operations, not the event loop itself. Per-thread accumulation state (`AnalysisObj` map / RpT `ThreadState`) is built lazily once per worker thread and merged into one result after the loop (`AnalysisObj::mergeFrom()` in `src/plotting_utilities.h`). See `util/ttreemt_prototype.cxx` for the standalone prototype this was validated against, and the "TTreeProcessorMT prototype" note below on a real thread-safety bug found via stress-testing. The `*_plot` executables are single-threaded — no event loop, no ntuple access, just histogram load + fit/derive + draw.
- **Constants**: modify analysis cuts in `src/clustering_constants.h`; always rebuild after header changes
- **Plot styling**: ROOT `TColor`/`TCanvas`; ATLAS style via `AtlasStyle.h`/`AtlasLabels.h`
- **Event loop**: `TTreeProcessorMT::Process()` hands each task its own `TTreeReader`; a fresh `BranchPointerWrapper` is bound to it per invocation (cheap — pointer setup only)
- **File identification**: inside the `TTreeProcessorMT` lambda, `reader.GetTree()->GetCurrentFile()->GetName()` gives the current file path, and `reader.GetTree()->GetReadEntry()` gives the local (0-indexed) entry number directly — both needed to reconstruct an event display command. This differs from the old sequential loop's `chain.GetReadEntry() - chain.GetChainOffset()` trick, which needed an outer `TChain` not available inside the lambda; a `TChain`-bound *sequential* `TTreeReader`'s `GetTree()->GetReadEntry()` returns the chain-global entry number, not local, which is exactly why that trick existed.
- **TTreeProcessorMT thread-safety gotcha**: `AnalysisObj`'s constructor calls `SetFillColorAlpha()`, which lazily registers a new `TColor` into ROOT's *global* color table the first time a given (color, alpha) pair is used — not thread-safe. Worker threads racing to build their first `AnalysisObj`/`ThreadState` concurrently corrupted the heap (intermittent `SIGABRT`, found via stress-testing, not a one-off). Fixed by serializing each thread's entire lazy-init construction (not just registry bookkeeping) behind its mutex — see the per-thread registry pattern in both executables' `main()`.
