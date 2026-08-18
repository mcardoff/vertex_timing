# Baseline at the canonical selection (`--vbs-deta=0.0 --vbs-mjj=500`)

Tag: `deta0p0_mjj500p0`. Supersedes every number in `results/README.md`,
which was taken on the old fiducial region (|Deta| >= 3.0, no mjj cut) and
should be read as historical only.

Sources: `clustering_hist` clusters 2883059 (mu=200) and 2883159 (mu=0).
Core fraction = fraction of accepted events whose assigned vertex time is
within +/- `PASS_SIGMA` = 60 ps of the truth HS time, integrated over the
sig/mix/bkg composition categories of `reso_*_hgtdtimes_<score>`.

## Core fraction

| sample | mu | N | HGTD | TRKPTZ | WAVES | WAVES_RECLUST | TEST_MISAS | WAVES_MISCL | WAVES_MISAS |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|
| vbf         | 200 |  519115 | 92.4% | 90.5% | 91.9% | 92.1% | 92.9% | 97.8% | 95.0% |
| vbf_mu0     |   0 | 2097581 | 99.8% | 99.9% | 99.6% | 99.3% | 99.9% | 99.6% | 99.7% |
| zjets       | 200 |   36531 | 73.8% | 61.9% | 57.9% | 44.4% | 64.2% | 93.4% | 61.2% |
| zeejets_mu0 |   0 |    4524 | 99.8% | 99.9% | 99.3% | 98.9% |100.0% | 99.3% | 99.5% |
| ttbar_mu0   |   0 |  941956 | 99.8% | 99.9% | 99.5% | 99.1% |100.0% | 99.5% | 99.6% |
| dijet       | 200 |   80383 | 89.9% | 86.9% | 86.1% | 85.0% | 89.4% | 97.7% | 89.6% |

### Reading

1. **At mu=0 every method is ~99.9%, on every process.** The 60 ps window is
   loose compared with the per-track HGTD resolution once a cluster is pure,
   so detector resolution is NOT what limits us. The entire mu=200 deficit
   (VBF 90.5%, Z+jets 61.9%) is pileup-induced contamination and
   misassociation. This turns "the residual needs better timing" from an
   elimination argument into a measurement: it does not.

2. **mu=0 samples carry no signal for the cluster-selection problem.** With
   one vertex there is one cluster and the label is trivial, so adding them
   to selection training contributes ~no gradient while diluting the mix with
   trivially-easy events. This is the point where our problem differs from
   the HGTD track-vertex-association BDT, which trains on mu=0 usefully
   (clean positive track-vertex pairs). Do not copy that sample mix over.
   They remain valuable exactly as used here: as a floor/control measurement.

3. WAVeS is still below TRKPTZ on Z+jets (57.9% vs 61.9%) and its
   re-clustering variant is far below (44.4%) -- unchanged in direction from
   the old selection, so the canonical cut did not create that gap.

## Topology acceptance: is the forward dijet real or pileup?

Rates below are conditional on passing lepton selection and the vertex/MIN_JETS
requirement, so they isolate the jet-pT/VBS-topology cut itself.

| sample | mu | lepton-pass | + vtx/jets | accepted | topology pass |
|---|---|---:|---:|---:|---:|
| vbf         | 200 | 1464400 | 1184243 |  519188 | 43.84% |
| vbf_mu0     |   0 | 4754000 | 4357023 | 2103384 | 48.28% |
| zjets       | 200 |  601799 |  519105 |   36592 |  7.05% |
| zeejets_mu0 |   0 | 1053408 | 1028598 |    4538 |  0.44% |
| ttbar_mu0   |   0 | 4934000 | 4932079 |  943820 | 19.14% |
| dijet       | 200 |  699500 |  323858 |   80402 | 24.83% |

- **VBF: mu200/mu0 = 0.91x.** Pileup does not help a genuine VBF topology
  pass -- it very slightly hurts, presumably by hijacking the candidate pair.
- **Z+jets: mu200/mu0 = 16.0x.** Only ~6.3% of the Z+jets events in our
  selection would be there without pileup. The forward dijet is supplied by
  pileup in ~94% of them.

This is the direct measurement of what was previously only inferred from the
loosened-cut scan: the Z+jets "signal region" is largely a pileup-jet region,
which is why a jet-proximity score (WAVeS) built on the assumption that the
forward jets are hard-scatter objects underperforms there while transferring
poorly from VBF.

### Caveat

`zjets` is DSID 601189 (Zee) + 601190 (Zmumu) combined; `zeejets_mu0` is
601189 only. The comparison is made on rates *conditional* on lepton
selection, where the hadronic activity of the two channels is equivalent, so
the ratio is fair. The like-for-like partner (`zeejets`, 601189 at mu=200) is
still to be run and would remove the caveat entirely.
