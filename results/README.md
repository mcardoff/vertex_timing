# Cluster-selection / vertex-timing ML results

State of the work as of the last condor runs. Raw metrics are in the `.json`
files here (per lambda, per sample, plus the full epoch history); the code is
`python/train_deepsets.py` and `python/train_transformer.py`, driven by
`condor/train_*.sub`. The working notebook is `training.ipynb`; superseded
exploration with its findings table is `training_archive.ipynb`.

Metric throughout is **core fraction**: the fraction of events whose assigned
vertex time lands within `PASS_SIGMA = 60 ps` of the truth hard-scatter time.

## Where things stand

| selector | dijet | vbf | **zjets** |
|---|---|---|---|
| TRKPTZ (the baseline being replaced) | 87.0% | 90.1% | 62.1% |
| WAVeS | 86.3% | 90.7% | 60.2% |
| **Deep Sets** (cluster selection, end-to-end) | **90.2%** | **93.0%** | **68.1%** |
| Transformer + MDN (no clustering) | 89.8% | 77.0% | 65.3% |
| — its mixture MEAN (= a plain regression) | 87.5% | 74.6% | 58.5% |
| constant-prediction floor | 27.2% | 26.8% | 26.9% |
| cluster oracle (a CLUSTERING ceiling) | 97.8% | 98.2% | 92.2% |

**Deep Sets is the best selector**, +6.0 on Z+jets over TRKPTZ and best on all
three topologies. Z+jets is the only sample with real headroom (62.1% against a
92.2% ceiling); vbf and dijet start at ~90%.

## What is established

- **The problem is SELECTION, not clustering.** An acceptable cluster exists in
  92.2% of Z+jets events; TRKPTZ picks it 62% of the time.
- **WAVeS is topology-bound** (+0.7 vbf, −1.9 zjets). A learned selector is not:
  every off-diagonal of the train/test topology matrix beats TRKPTZ.
- **Changing the OBJECTIVE is the only thing that ever worked.** Training
  end-to-end on the selection objective gained +3.5 on zjets over learning
  `P(HS)` and pooling it with a fixed `pT` rule.
- **Six levers were null**: cluster-feature engineering, per-track tagger
  features, model capacity (5.7×, slightly worse), pooling form (gated vs sum),
  the jet-association block (dropping it *helps* zjets), and an auxiliary
  per-track head on `truth_is_hs` (λ scanned, within noise every time).
- **The residual is not reachable from per-cluster features.** A stage-2
  re-ranker with 21 points available scored −0.6.

## Decomposition of the remaining Z+jets gap

| | |
|---|---|
| 62.1 → 68.1 | what selection recovers from reco features — **achieved** |
| 68.1 → 80.2 | what a *perfect* per-track HS tagger would add — blocked; §15 showed these features cannot tag better |
| 80.2 → 92.2 | **time quality**, not track identity — `Σ pT·truth_is_hs` uses perfect HS identity with *real* times and still stops at 80.2% |

The last row is why the hypothetical ~90 ps calorimeter timing matters: it is
exogenous information the reco features do not contain, and the argument for it
is now a measurement rather than an inference.

## Open threads

- **Transformer/VBF collapse.** 77.0% where every other method reaches 90%+, and
  it plateaued by epoch 4, so not undertraining. Best hypothesis: the cluster
  methods output a cluster's *mean track time* (precision ~σ/√N, and VBF's HS
  clusters have the most tracks), while the transformer regresses a free scalar
  and must learn that averaging. Proposed fix, untested: emit per-track weights
  and output `Σ wᵢ·tᵢ` instead of an unconstrained μ.
- **The MDN is working** — argmax-weight beats the mixture mean everywhere, by
  +6.8 on zjets. A plain regression would have lost that to mode-averaging.
- **`novbs` loose export ready, not yet trained on**: 94,187 events vs 67,074
  tight (+40%). The VBS jet-pair requirement was worth that much; the *forward
  jet* requirement turns out to be the real gatekeeper.
- **μ0 / ttbar samples requested.** μ0 VBF gives the intrinsic timing floor;
  Z+jets and ttbar at μ0 measure how much of the forward-jet population is
  pileup (their selection efficiency should collapse if it is).
- **`Track_leptonID` provenance still unverified.** It defines the Z→ℓℓ fiducial
  region via `passLeptonSelection`. The leakage sentinel shows its derived
  features are worth ~0, so it can be dropped from the inputs at no cost — but
  if it is truth-matched, the *selection* is partly truth-based.

## Methodology notes worth not re-learning

- **AUC is not the figure of merit.** It saturates at 1.000 while core fraction
  still has 6.4 points to climb, and across five tagger variants AUC and core
  fraction were *anti*-correlated.
- **Ablate per sample, never pooled.** VBF has ~8 points of headroom and hides
  everything; a truth feature worth ~14 points on Z+jets ablated to exactly 0.0
  on VBF.
- **A scan must clear its own noise.** Picking the best λ by bare argmax once
  reported "the auxiliary head helped" on a +0.006 margin against a 0.05
  within-run sd. Both trainers now require >2 sd or select the control.
- **Correlate a suspect feature against the LABEL**, not against whatever truth
  quantity is nearest to hand. `mean_nhgtd_primary` had ρ=−0.065 with
  `truth_purity` but −0.44 with `abs_dt`.
