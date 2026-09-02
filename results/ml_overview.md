# Cluster-selection ML: the overall result (2026-08-27 → 2026-08-30)

The consolidated record of the machine-learning cluster-selection effort. The
effort is **concluded**: its endpoint numbers below stand as the measured
ceiling for what reconstruction information contains, and its mechanism
analysis directly seeded the classical TZP selector that was deployed instead
(see `results/score_grid_search.md` and the TZP section of CLAUDE.md). Every
number here is at the canonical selection (`--vbs-mjj=500`) with the
shard-invariant `(sample_id, file_idx, event_num)` group key unless stated;
metric is core fraction (|t − t_truth| < 60 ps) on held-out test folds.

## The ladder, end to end

| method | vbf | **zjets** | dijet | ttbar |
|---|---:|---:|---:|---:|
| TRKPTZ (baseline) | 90.45 | 61.87 | 86.88 | 87.08 |
| WAVeS | 91.1 | 59.4 | 86.2 | 86.3 |
| **TZP** (classical, deployed) | 92.08 | 64.31 | 88.58 | 88.84 |
| Deep Sets selector alone (22 feat) | 93.30 | 66.94 | 90.00 | 90.30 |
| full guarded pipeline (pooled meta, thr 0.5) | **95.26** | **69.79** | **91.00** | **91.75** |
| naive best-|Δt| oracle | 98.3 | 92.6 | 97.7 | 98.0 |
| **honest ceiling** (see "oracle is inflated") | — | **~74.3** | — | — |

Reading it: the classical TZP captures **roughly half of the Deep Sets
selector's gain on every sample** (zjets 2.44 of 5.07 points, vbf/dijet/ttbar
54–57%) with zero deployment cost — no model, no inference, no training data.
The full pipeline holds a further ~+5.5 on zjets over TZP, at the cost of a
four-model stack (per-track tagger → learned cluster score → learned t0 →
meta chooser). μ=0 controls sit at 99.9% for everything — those samples are
excluded from training (one vertex ⇒ trivial label, no gradient).

Raw run records: `final_method_eval_results.md` (the condor evaluation, cluster
2951035), `zjets_seventy_results.md` (the run that first crossed 70% on
zjets), `feature_ablation_4sample.md`, `guarded_pooled_results.md`,
`deepsets_tight.json` / `transformer_tight.json` (full epoch histories),
`final_method_eval.json`, `guarded_indomain.json`, `guarded_transfer_vbf.json`.
Narrative pages: `final_method_report.html`, `method_explainer.html`,
`status_page.html` + `zjets_failure_gallery.html` (both predate the oracle
correction below — their "25-point gap" framing is measured against the
inflated target), `README.md` (previous round, pre-canonical selection —
banner explains why its numbers are not comparable), `PLAN.md` (the nine-item
runbook and its outcomes).

## What the pipeline actually leans on (measured importances)

From `model_inputs_report.py` → `model_inputs.json`, fitting each model on
zjets (full analysis in `method_explainer.html`):

- **The per-track HS tagger is the load-bearing component.** Permuting its
  output costs the e2e t0 net 21.8 points (next feature: 6.9); the learned
  cluster score is 66% tagger-derived. The round-1 tagger is structurally
  organised around `closer_to_pu_than_pv` (65.8% of gain — largely restating
  the z-association); round 2's neighbour columns take 37% collectively and
  demote it to 26%, which is displacement, i.e. real added information.
- **The learned cluster score is not a better ΣpT — it is a consistency check
  against the aggregate t0.** Raw `trkptz_score` carries 0.19% of its gain;
  its top two features are distance-to-aggregate terms (42.9% combined). The
  aggregate t0 is therefore a *prerequisite* for the learned score, not an
  optional add-on — partial adoption that takes the score without the
  aggregate is incoherent.
- **`timeRes` separates HS from PU at AUC 0.508 — a coin flip** — yet it is
  the only weight the inverse-variance mean uses. That is the mechanism behind
  the re-timing gain, decomposed exactly: inverse-variance 62.12 →
  tagger-probability-weighted 63.84 (+1.72) → double-Winsorised 64.03 (+0.19).
- **The meta chooser is the weak link**: it sees the 7.35% of events where
  cluster and aggregate disagree, and there the aggregate is right only 26.7%
  of the time. Its single-variable AUCs are the weakest of the five models —
  the quantitative form of "two-thirds of the union spread goes uncaptured",
  which is the same union(cluster, aggregate) frontier TZP left open.

## Established findings

1. **The problem is selection, not clustering.** An acceptable cluster exists
   in the collection for the large majority of failures; agrees with the
   clustering-side failure decomposition (67.4% timing misassignment, 26.5%
   selection, 3.1% misclustering).
2. **The naive oracle overstated the zjets headroom ~2.5×**
   (`oracle_is_inflated.md`). Best-|Δt| is a best-of-N statistic at N≈5.7: a
   null with event structure destroyed still scores 53.4% (real: 92.6%), and
   13.6% of oracle picks have truth purity exactly 0. The genuinely learnable
   headroom over TRKPTZ is ~12.5 points (realistic ceiling ~74.3%), of which
   Deep Sets delivers 5.4 — 44%. Rule: **null-test any best-of-N oracle
   before quoting it as headroom.**
3. **The learned selector transfers across topologies; WAVeS does not.** 15 of
   16 train→test cells beat TRKPTZ; ttbar alone is statistically
   indistinguishable from training on all four samples; the sole failing cell
   is zjets→vbf. WAVeS's +0.7 vbf / −1.9 zjets topology-boundness is confirmed
   on three independent legs.
4. **zjets is information-limited, not data-limited.** 5.5× the training
   events moves core fraction 0.2 points inside a ±0.8 band. More Z+jets MC
   buys nothing; the residual needs new information (the time-quality term).
5. **The training objective was the only lever that ever moved the number**
   (multi-positive listwise selection loss, +3.5 zjets over P(HS)+pooling).
   Null levers, all measured: capacity (5.7×, slightly worse), pooling form,
   auxiliary per-track head (every λ), cluster-feature engineering, tagger
   features into the selector, the jet-association block (dropping it *helps*
   zjets), stage-2 re-ranking on per-cluster features (−0.6 with 21 points
   available), seed ensembling (structurally confined to init variance, which
   is the small term; −0.0..+0.1), and a reco-feature cluster-trust classifier
   (cannot flag confidently-wrong selections).
6. **Feature groups are sample-antagonistic** (`feature_ablation_4sample.md`):
   vertex-relation features help only zjets (+0.34, 5/5 seeds) and are mildly
   negative elsewhere; jet/LOO features do the reverse (+0.12 vbf, −0.16
   zjets). A single feature set is a compromise, and a local-VBF-only ablation
   is actively misleading.
7. **Run-to-run variance is split-dominated** (zjets ±1.9 split+init vs ±0.5
   init-only; condor itself is bit-deterministic). Any claimed effect must
   clear the split band — multi-seed or it did not happen.
8. **Epoch selection by validation loss, never by val-argmax core fraction** —
   the val-macro path collapsed to 17.6% on weights that scored 91.7% on test
   (root cause still open; a TRKPTZ control column through every tensor guards
   the recurrence). Val-loss selection is better on every sample.
9. **Transformer/MDN**: the mixture argmax beats the mixture mean everywhere
   (+6.8 zjets) — the MDN is doing its job — but the no-clustering transformer
   collapses on vbf (77%), plateaued by epoch 4; leading hypothesis is that it
   must *learn* the averaging that cluster methods get structurally. Untested
   fix: emit per-track weights, output Σwᵢtᵢ.

## Methodology (the transferable part)

- AUC is not the figure of merit — it saturated at 1.000 with 6.4 points of
  core fraction left, and across five tagger variants the two were
  *anti-correlated*.
- Quote per sample, never pooled: VBF's small headroom hides everything.
- The label is `|Δt| < 60 ps`, never cluster purity — purity measures
  provenance, not time correctness, and max-purity selection amplifies the
  dominant failure mode.
- Gain importance ≠ standalone separation (65.8% gain at AUC 0.664 vs a
  neighbour column at 0.754): pull both before concluding either.
- The group key must be shard-invariant; `hadd` for exports, `hist_merge` for
  histogram files; μ=0 stays out of the training mix (unlike the HGTD
  track-association BDT — do not copy that recipe).

## Where this leaves the ML option

Reopen only with a concrete reason (PLAN item 9). The two open technical
threads are the val-macro selection bug and the transformer/vbf collapse.
The strategic thread is the same one the classical side ends on: the
union(cluster answer, aggregate t0) spread that no chooser — learned meta or
classical — has captured, and whose strongest known reco-only correlate is
the `|t_in − t_out|` split-clustering agreement flag. Group-1/Group-3 columns
from the extended AF re-export (`on_pv`, `chi2_ndf`, `btag_*`) were never fed
to any model and remain untested.

Training/eval infrastructure retained in-tree: `python/train_deepsets.py`,
`python/train_transformer.py`, `python/ensemble_eval.py`,
`python/compare_models.py`, `python/final_method_eval.py`, driven by
`condor/run_training.sh` + `condor/train_*.sub`; exporter
`util/export_training_data.cxx` (see CLAUDE.md, ML Data Export). The ~35
one-off study scripts behind the results in this directory were removed from
`python/` after their findings were captured here and in the per-study
documents — recover any of them from git history if a measurement needs
rerunning.
