---
name: physics-findings
description: Recorded physics results for vertex timing — TZP mechanisms, failure decomposition, where HGTD's gain lives (R1/R2), gate severity, efficiency ceiling, inflation calibration, the z-association study, reference-macro replication, and the ML / GNN history. Load before quoting a number, proposing a new score, or re-running a closed study.
---

## TZP: mechanisms and the trail behind it

Full trail: `results/{trkptz_pv_augmentation, injet_decomposition,
sample_independent_timing, guarded_allsamples, qoverp_term, score_grid_search,
t0_integration, failure_modes, tiebreak_scan, jet_term_consistency,
classical_scoring_brainstorm, split_clustering}.md`. Key mechanisms:

- **WAVeS decomposes into selection (+0.5) and in-jet re-timing (+0.8), and the
  re-timing is a pure additive term** — identical to two decimals across
  TRKPTZ/TZ/WAVeS selectors on every sample tested. Unconditionally applied
  (deployed WAVeS) it is NEGATIVE on 3 of 4 samples (zjets −1.4); the ≥3-track
  guard makes it safe everywhere because the sample dependence then emerges
  from the data (fires on 70% of VBF events, 9% of zjets) instead of being
  configured in. The sign is set by the share of HS *timing weight* inside
  jets: 48% VBF vs 10% zjets — a statistics effect, not purity (zjets in-jet
  purity is still 2.1× enriched and the subset still times worse).
- **Gentle-beats-steep, uniformly.** Every steep/binary z- or agreement-based
  penalty lost (PV-proximity veto −1.27, vertex-WAVeS multiplier −5.7, t0
  fallback −2.8); the gentle continuous versions of the same information won
  (+1.76). At ~106 vertices/event the HS cluster is often NOT nearest the PV.
- **TRKPTZ's cluster-level 1.5 was stale.** Never tuned by anyone; jointly with
  the per-track term the optimum softens it 2.5× to 0.6 (α 0.7→1.0, γ 1.2→0.45,
  each 1-D optimum having been conditioned on the others' old values).
- **The three precision sums (z0/d0/qP) are one signal** (mutual corr
  0.86–0.95); d0 is the best packaging and stacking a second adds nothing.
- **Raw |Δz| in mm beats z0 significance** as the per-track weight (+1.76 vs
  +1.21 on zjets) — same conclusion as the z-association study reached about
  `sqrt(var_z0)`.

## Machine Learning

Two distinct ML eras, neither active in the C++ analysis today.

### Cluster-selection ML (2026-08, concluded — see `results/ml_overview.md`)

A GBDT / Deep Sets / transformer effort against the same cluster-selection
problem TZP now solves classically. **`results/ml_overview.md` is the
consolidated record** — the endpoint ladder (TRKPTZ 61.87 → TZP 64.31 → Deep
Sets selector 66.94 → full guarded pipeline 69.79 on canonical zjets, against
an honest ceiling of ~74.3 after the best-of-N oracle correction), the
measured mechanism (the per-track HS tagger is load-bearing; the learned
cluster score is 43% agreement-with-aggregate-t0 and 0.19% raw trkptz_score),
the transfer matrix (learned selection transfers across topologies, ttbar
alone ≈ pooled; WAVeS doesn't transfer), the null levers, and the methodology
rules (AUC anti-correlated with core fraction; null-test best-of-N oracles;
split-dominated variance). Individual run records live beside it in
`results/` (`PLAN.md`, `README.md`, `oracle_is_inflated.md`,
`feature_ablation_4sample.md`, `final_method_eval_results.md`, the `.json`
histories, the `.html` reports).

Retained infrastructure: `util/export_training_data.cxx` (the exporter — ML
Data Export section above), `python/train_deepsets.py` /
`train_transformer.py` / `ensemble_eval.py` / `compare_models.py` /
`final_method_eval.py`, and `condor/run_training.sh` + `condor/train_*.sub` +
`condor/probe_env.*` (GPU-less AF training jobs). The ~35 one-off study
scripts that produced the `results/` documents were **removed from `python/`
once their findings were captured** — recover from git history if a
measurement needs rerunning; do not rewrite them from scratch.

### Event-graph GNN (2026-09, exploratory — see `results/gnn_event_graph.md`)

`python/train_gnn.py` represents an event as a graph over its tracks (all-pairs
edges, typed same-cluster / cross-cluster, so the existing partition is a prior
the model can override rather than the topology). It never re-partitions; the
time head emits per-track weights and `t_hat = Σ w_i t_i / Σ w_i` keeps the
averaging structural. Both heads are residual from the incumbent with zero-init
last layers, so epoch 0 IS TRKPTZ selection with inverse-variance timing.

**Finding (VBF, local export, 5-seed paired sweep): message passing is worth
+0.80 ± 0.14 on the timing head and +0.26 ± 0.10 on selection.** Its value is
trust estimation — a track's own features cannot say whether its time is
reliable (timeRes separates HS from PU at AUC 0.508), but its neighbourhood
can. Two thirds of the total gain over TRKPTZ is the listwise objective and
pooled selector, not the graph.

**Not established, and why this is exploratory rather than concluded:** VBF
only, on the small local export (23k fit events). VBF carries 47% of HS timing
weight inside jets against zjets' 10%, which is exactly the quantity governing
whether jet-neighbourhood information means anything — expect the +0.80 to
SHRINK on zjets. That is the falsifying test and it has not been run. Nor has
`--no-residual`. The residual-init control fired on first use against the wrong
denominator (whole test fold vs the 1,000-event tensor) — compute controls on
the tensor's own events.

### Legacy DNN/TMVA artifacts (historical)

There is no `src/ml_model.h`, no ONNX runtime integration, and no
`Score::TEST_ML` in the current `Score` registry. What remains from that
earlier effort:
- `share/models/` still holds artifacts (`model.onnx`, `model_weights.json`,
  `normalization_params.json`, a TMVA `.weights.xml`) from that earlier
  work, but nothing in `src/` currently loads or runs inference against
  them.
- `share/scripts/export_onnx_model.py` (SavedModel → ONNX → JSON weights
  pipeline) and `share/docs/ML_MODEL_INTEGRATION.md` /
  `ONNX_INTEGRATION_GUIDE.md` document that pipeline for reference, but it
  has no active consumer.
- `src/clustering_structs.h`'s `Cluster::calcFeatures()` still computes and
  normalises the same 8-feature vector (used historically as DNN input); its
  `rawDeltaZ`/`rawDeltaZResunits` outputs are reused by `updateScores` for
  TRKPTZ, but the normalised `features` vector itself is otherwise unused
  today.


## rpt_v5 / failure-decomposition findings

Results as of 2026-08-20, from `rpt_v5` (grid samples, |Deta| >= 0 and
m_jj >= 500 GeV) and `failure_decomposition` (local VBF). Numbers are pileup
rejection (1/mistag) at fixed hard-scatter efficiency unless stated.

### Where the loss actually lives

⚠️ Updated 2026-11-09 (notes/gp_talk.md, "failure_decomposition rerun on the talk's t0s"): with the TDR's track-time cleaning, the timing rung falls from 67% to 49% of WAVeS failures (local VBF), and the wrong-time route proper ("> 50% of the HS pT >= 3 sigma off") from 18-20% to 10%. Counting also the failures that removing the mistimed HS tracks fixes (oracle fix_mistimed_*, 2026-11-10), misassigned HS times cause 35% of local-VBF WAVeS failures with the ntuple times and 17.6% with cleaning: the largest single identifiable cause before cleaning. Most of the timing rung is the ideal-time oracle, which overlaps heavily with selection and misclustering, and which is BIASED on the grid samples (pileup tracks have no truth link). "Another cluster passes" is inflated by pileup accidents: the hard scatter's own cluster passes in 42%, not 83%. The table below is the August, uncleaned picture.

⚠️ Measured as repairs, 2026-11-11 (notes/gp_talk.md, "Fixing each failure mode in the main program"): clustering_hist's WFIX_* rows fix each mode with HS truth on the same events. Local VBF, talk inputs, WAVeS 92.71 -> 99.49 in the ladder's order: re-timing +0.24, misassigned HS times +0.87, selection +2.31, misclustering +2.25, HS track efficiency +1.01 (+3.05 alone), HS time resolution +0.10. The resolution moves only when pileup leaves the cluster (17.4 -> 13.1 ps). The gated oracle rows (WAVES_MISAS, WAVES_MISCL) drop 26-35% of events and overstate their mode; do not quote them as a fix's worth.

`failure_decomposition` assigns each failing event one principal cause:

| cause | share |
|---|---|
| Timing misassignment | **67.4%** |
| Selection failure | 26.5% |
| Misclustering | **3.1%** |
| Track efficiency loss | 0.7% |
| Irreducible | 2.3% |

**The clustering is not the bottleneck.** Grouping all truth-HS tracks into
one perfect cluster almost never rescues a lost event, so the good cluster is
essentially always present and simply not chosen. Selection dominates in
sparse events (52.9% at 0–1 forward HS tracks) and that is where nearly all
failures live (78% fail rate at 0–1 tracks, ~1% above 14); timing
misassignment dominates in track-rich ones.

The scenario ladder says the same thing independently. VBF forward >40 GeV at
0.90 efficiency: `zonly 65.0 -> waves 87.1 -> waves_ideal 106.2 -> truth
111.3`. So **idealising times and selection together is worth ~22%, and the
clustering leaves only ~5% more on the table**. dijet agrees (63.9 vs 68.0,
~6%). Two independent measurements plus the 3.1% above.

Caveat: part of the residual `waves_ideal -> truth` gap is not a clustering
error at all — truth is handed a 10 ps vertex time, while even a perfectly
chosen cluster must estimate t0 from a handful of 30 ps tracks. Separating
those is open.

### The gain is concentrated in the topology

ITk -> WAVeS at 0.80 efficiency:

| selection | VBF | dijet | Z+jets |
|---|---|---|---|
| Inclusive forward | x1.3 | x1.2 | x1.2 |
| R1 (both tags forward) | **x1.87** | **x2.53** | x1.19 |
| R2 (fwd PU + central HS) | x1.17 | x1.18 | x1.04 |

The inclusive number understates HGTD roughly two-fold, because it averages
over jets where timing is irrelevant: pileup jets with no PV0 tracks to begin
with, and HS jets with ample track support.

### R2 fails because WAVeS is captured by the fake

The score ordering **inverts** between regions, consistently across all three
samples — WAVeS wins R1, the jet-agnostic TRKPTZ wins R2 (VBF x1.32 vs x1.17;
dijet x1.44 vs x1.18; Z+jets x1.09 vs x1.04).

WAVeS weights tracks by proximity to forward jets. In R2 *the only forward jet
is the pileup one*, so WAVeS up-weights exactly the tracks of the fake vertex,
drags t0 toward it, and the fake then passes its own gate. The fix is either
regional (use TRKPTZ in R2) or structural (compute t0 from tracks outside
jets, so no single jet can capture the score).

### Timing discriminates by severity, not frequency

Per forward jet, VBF, WAVeS gate vs ITk-only:

| jet | n | raised | lowered | -> zero | mean loss |
|---|---|---|---|---|---|
| HS | 56,342 | **0** | 26.5% | 0.6% | 18.4% |
| PU | 15,307 | **0** | 7.3% | 6.0% | **90.4%** |

The gate fires on HS jets more often than PU ones, but a touched PU jet loses
90% of its R_pT against 18% for an HS jet. Most PU jets are untouched because
they had nothing to remove (they pass the z-association at 0.8%). "Raised" is
0 by construction: `applyTimeGate` returns a SUBSET, so R_pT can only fall.

### The efficiency ceiling is the binding constraint

A jet the gate empties goes to R_pT = 0 and becomes indistinguishable from
pileup, so every timed scenario has a LOWER maximum reachable HS efficiency
than ITk-only (VBF 30–40: 93.8% -> 91.0%; Z+jets: 87.6% -> 82.4%). This, not
resolution, is why curves collapse above ~0.93 and why ratio panels fall back
to 1. Proposed fix, NOT implemented: fall back to the jet's ITk-only R_pT when
the gate empties it — the jet-level analogue of TDR self-tagging.

### Calibration is correct but near-irrelevant

The pull denominator is sqrt(f^2 sigma_vtx^2 + sigma_trk^2) with sigma_vtx
~9 ps against sigma_trk ~26 ps, so a **17% change in inflation moves the
effective gate width by 3%**. Worth having right; not worth another iteration.
Offsets are 0.6–1.5 ps (negligible). Non-Gaussian tails outside a 150 ps core
are 14–15% (VBF), 16–18% (dijet), 21–25% (Z+jets) — the tail, not the core, is
what a Gaussian pull cannot handle.

Z+jets is limited by t0 AVAILABILITY, not by the algorithm: Athena produces a
valid vertex time for 80.1% of VBF events, 73.6% of dijet, and **16.5%** of
Z+jets.

### The z-association is mistuned for TRKPTZ, and irrelevant for WAVeS

From `assoc_study_hist`/`assoc_study_plot`, local VBF, 47,535 selected events,
11 association rules evaluated against the same events (2026-08-26). Full
write-up: `share/docs/Z_ASSOCIATION_STUDY.md`.

The incumbent `significance < 3.0` keeps a track list that is **68% pileup**,
and it is the second-worst of the eleven rules for TRKPTZ:

| score | best rule | vs. incumbent |
|---|---|---|
| TRKPTZ | `getNewDzpara × 1.0` | **−1.58 ps core σ** (15.95 → 14.37), **+0.92 %pt core fraction** (90.38 → 91.30) |
| WAVeS  | `getNewDzpara × 1.2` | +0.13 ps, +0.10 %pt — i.e. nothing |

Both families improve as they tighten and then turn over (significance at
~2.0, the parameterization at ~1.0), so there is a real optimum and the
incumbent is not near it. Simply passing the existing `--dzpara` flag
(`getNewDzpara × 1.4`) already captures ~3/4 of the gain with no code change.
The TRKPTZ gain is concentrated at low multiplicity — σ(nHS≤5) improves
29.3 → 25.3 ps (−13%), which is exactly where the failure decomposition says
the failures live.

**The two scores want opposite things, for a mechanical reason.**
`Cluster::calculateTime` recomputes the WAVeS time from only the tracks within
ΔR < 0.4 of a forward jet, and that in-jet filter already does the
pileup-shedding a tighter z-association would do. The two mechanisms are
redundant, so tightening can only cost WAVeS hard-scatter tracks — and it does,
monotonically (`dzpara × 0.6`: +1.86 ps, −1.83 %pt, and worse at *every*
multiplicity, not just sparse events). Expect any future selector that
recomputes its time from a track subset to behave the same way, and to want a
LOOSE association.

**~59% of WAVeS's advantage over TRKPTZ is just the z-association.** The
core-fraction gap is 1.38 %pt at the incumbent rule but 0.56 %pt when each
score gets its own best rule. On core σ the ordering inverts outright: TRKPTZ
at `dzpara × 1.0` (14.37 ps) beats WAVeS at any rule (best 14.95 ps). WAVeS's
residual edge is in core FRACTION, not core width — i.e. it suppresses gross
misassignment, consistent with the failure decomposition, not resolution.

**The parameterization is a better shape, not merely a tighter cut.** Compared
at matched track purity the two families separate, and the gap widens the
tighter the working point: at ~31.5% purity the parameterization is worth +1.0
points of HS recall and the two perform identically; at ~43.7% purity it is
worth **+5.6 points of recall at higher purity**, and TRKPTZ gains 0.5 ps and
0.6 %pt from the shape alone. Likely cause: `sqrt(var_z0)` is a per-track
covariance that does not capture the η dependence of the true z0 spread, and
`Z0_VAR_INFLATION = 1.15²` is independent evidence it is mis-estimated — a
single global factor cannot repair an η-dependent mis-modelling. (The factor was
removed 2026-11-16, f897b1e; this study's significance cuts N are 1.15 N of the raw covariance.)

**Union with ghost association: tested and rejected.** The idea was that ghost
association, being derived independently of z, might recover HS tracks the
z-cut misses. OR-ing it in costs **20–24 %pt of core fraction** for BOTH
scores (TRKPTZ 91.30 → 68.47, WAVeS 91.80 → 72.98, on a `dzpara × 1.0` base).

The failure mode is worth remembering because it is not the obvious one.
**Core σ barely moves (~1 ps) while core fraction collapses** — so the added
tracks are not dragging the HS cluster's time, they are causing the wrong
cluster to be *selected*. Confirmed differentially: every z-only rule
saturates at ~99.9% core fraction once there are ≥20 forward HS tracks, while
the union plateaus at **84%** and never saturates. It picks wrong one time in
six even in the most track-rich events.

Why: the added pileup is **time-random** (mean |Δt| to the truth HS vertex of
224 ps, against 229 ps for a track drawn at random from the event — measured
in `util/scratch/ghost_assoc_diag.cxx`), but sits inside forward jets and so
carries real pT. It builds *decoy clusters* rather than poisoning the HS one,
feeding the selection channel that the failure decomposition already names as
dominant. WAVeS is hit as hard as TRKPTZ here despite being immune to every
z-only change, because jet-proximate tracks are exactly what its score
up-weights.

Also established while testing this: **ghost association is numerically the
same as `filterTracksInJets` at ΔR < 0.4** (1,142,857 tracks / 24.44% HS vs
1,139,170 / 24.54%), so it is not new information relative to what WAVeS
already uses.

**Not yet acted on.** `MAX_NSIGMA` is unchanged: this has only been run on the
local VBF sample. Z+jets in particular should be checked first, since Athena
produces a valid vertex time for only 16.5% of its events. Rerunning is cheap
(~35 s hist + ~15 s plot on the full local sample; condor via `--sample=`).

### Reference-study replication

Central region at the reference configuration (`util/myJet_ana_fr.C`, present
in the main working tree but not committed — it covers both regions via
commented-out eta lines and carries the z0-significance cut as a commented
alternative): the
**timing improvement factor reproduces** (ours 1.63 vs their ~1.6 at 0.800)
but both absolute rejections sit ~2.3x low. A uniform offset cannot be caused
by the time gate, so the open question is the pileup population or its
normalisation, not the timing model. Jet collection (EMTopo) and the
isHS/isPU definitions are confirmed to match.

