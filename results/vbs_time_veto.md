# Event-level timing veto on the VBS pair: jet vs jet, or jet vs t0

`util/vbs_time_veto.cxx` + `python/vbs_time_veto_plot.py`, 2026-09-28.
Local VBF only so far. The grid run (VBF signal, Z+jets background) is
pending; commands are at the end.

## The question

This comes from the feedback of Friday 2026-09-25. In R1 (both tagging jets
forward, one HS and one PU), no vertex t0 is needed to reject the event. The two
jets can be tested against each other. The comparison asked for is signal (VBF
H→inv) and background (Z+jets) yields, weighted, after the analysis recoil
cut, under two vetoes:

| method | veto the event when | can act on |
|---|---|---|
| **jet vs jet** | both legs are timed and \|t_A − t_B\| / √(σ_A² + σ_B²) ≥ 3 | pairs with both legs timed (R1 geometry) |
| **jet vs t0** | the event has a t0 and **any** timed leg has \|t_leg − t0\| / √(σ_leg² + (f σ_t0)²) ≥ 3 | every pair with a timed leg: R1, R2, others |

Each jet vs t0 row uses a different t0: TRKPTZ, WAVeS, HGTD (Athena
`RecoVtx_time`, only where valid), TZP, and the truth HS time (perfect t0,
σ = 0). f is rpt_v5's per-sample t0 inflation (`inflationFor`, now in
`util/rpt_v5_common.h`).

**Selection:** the region plots' own, so an efficiency here describes a
region-plot column (`vbs_region_diag` + `vbs_region_stack.py`):
- vertex quality; Z→ℓℓ and lepton-jet overlap removal (Z+jets only);
- ≥ 2 jets > 30 GeV, ≥ 1 of them in 2.38 < |η| < 4.0;
- JVT + fJVT loose before pairing;
- the max-m_jj opposite-hemisphere pair, m_jj ≥ 200 GeV, no |Δη| cut.

**Recoil cut** (applied offline):
- VBF: truth E_T^miss = |Σ p_T^ν| > 200 GeV. This needs the truth record,
  which the skims drop, so VBF is read from the raw ntuples.
- Z+jets: pT(ll) > 200 GeV, as in `vbs_region_diag`.

**Weights:** every efficiency is Σw(kept) / Σw with the generator weight.
Yields use the normalisation in `rpt_region_distributions.md`.

## The jet time is the whole problem

"The ghost-associated tracks with valid time inside each jet", averaged as they
are, do not measure the jet's time at μ = 200. Local VBF, forward jets
above 30 GeV:

- A forward jet's ghost area holds **~31 tracks, ~23 of them timed**.
- **Only 21% of an HS jet's timed ghost tracks come from the hard scatter.**
  The rest, ~24 tracks per jet, is pileup that falls inside the jet's
  catchment area.

So the inverse-variance mean of all of them (`all_`, kept in the tree for the
record) is a pileup average:
- χ²/ndf ≈ 65;
- the core spread about truth is 23× its quoted σ;
- 11% of HS jets are within 3σ of t_HS;
- **25% of HS+HS pairs pass the jet vs jet test, against 21% of R1 pairs.** It
  cannot tell them apart.

The default (`core_`) estimator:
- keeps timed ghost tracks with the standard kinematics (1 < pT < 30 GeV,
  quality flag) within ΔR < 0.2 of the axis (the R_pT cone);
- clusters them in time (iterative, 2σ);
- takes the inverse-variance mean of the highest-ΣpT cluster.

It raises the chosen tracks' HS share to ~73%, and puts 87.7% of HS legs
within 3σ of t_HS (raw σ, full local sample). It was chosen on local-VBF truth (19,547 events, both legs in
2.4–4.0, raw σ). HS+HS pass rate at 3σ, varying one ingredient at a time from
the default:

| variant | HS+HS kept | R1 kept |
|---|---:|---:|
| **ΔR 0.2, 2σ, highest-ΣpT cluster (default)** | **83.3%** | 20.2% |
| ΔR 0.4 | 77.1% | 19.0% |
| 3σ clustering | 81.6% | 19.4% |
| cluster with the most tracks | 82.3% | 18.1% |
| no clustering (ΔR 0.2, pT > 1 GeV, plain mean) | 39.4% | 23.0% |
| no clustering, all timed ghost tracks | 24.7% | 21.3% |

**Its quoted σ is still too small.** It is a statistical error that knows
nothing of the pileup left in the cluster or of mis-assigned hits. Over
40,555 paper-HS legs (full local sample):
- the per-jet pull |t_jet − t_HS| / σ has a Gaussian-core width of **1.50**
  (1.4826 × median);
- rpt_v5's RMS/quoted t0 definition gives 1.80;
- 5.9% of HS legs are more than 150 ps off.

The plotting script applies σ_jet × (that width, measured per file) by
default, as rpt_v5 applies its t0 inflation. `--jet-infl 1` gives the raw
prescription.

## Validation (all pass)

- **Event set.**
  - Local VBF with `--jvt=loose` selects 38,524 events. That is exactly the
    m_jj ≥ 200 set of `vbs_region_diag --jvt=loose`: 0 events differ.
  - m_jj, |Δη|, both legs' pT and |η|, and the paper labels agree event by
    event (m_jj to 2×10⁻⁴ GeV, the diag's float storage).
  - Without the recoil cut, R1 / R2 are 4,183 / 2,579, the known local
    composition counts.
- **t0.**
  - Against rpt_v5_hist's `--jvt=loose` region tree, every one of the 6,762
    region events passing the composition preselection is present.
  - t0, σ_t0, validity and inflation are identical: Δt0 = 0 for hgtd, trkptz,
    waves and tzp.
  - The other 1,537 region rows are exactly those with `n_jets_fwd_acc` = 0.
- **Z+jets path** (lepton selection, overlap removal, pT(ll), Z+jets
  inflation), on the 109 grid Z+jets events picked for the region displays:
  - 107 pass the selection, and all 107 are in the grid `jvtLoose` region
    tree;
  - their t0s are identical;
  - pT(ll) equals `zjets_jvtLoose_zpt.root` to 5×10⁻⁵ GeV.
- **Truth MET:** 99.7% of local VBF events are above 75 GeV, matching the
  generator MET75 filter.

## Local VBF: signal efficiency

JVT + fJVT loose, m_jj > 200 GeV, E_T^miss(truth) > 200 GeV: 6,290 MC
events. The local sample stands in for the grid production. Figure:
`figs/time_veto/local_jvtLoose_sig_eff_mjj.(pdf|png)`.

| veto | ε_S, σ_jet × 1.50 (default) | raw σ | ≥ 2 tracks per jet time |
|---|---:|---:|---:|
| **jet vs jet** | **94.2 ± 0.3%** | 92.8% | 97.1% |
| jet vs t0: TRKPTZ | 84.7 ± 0.5% | 81.1% | 91.0% |
| jet vs t0: WAVeS | 86.0 ± 0.4% | 83.6% | 92.1% |
| jet vs t0: HGTD (Athena) | 86.3 ± 0.4% | 82.6% | 91.5% |
| jet vs t0: TZP | 85.4 ± 0.5% | 82.6% | 91.5% |
| jet vs t0: truth (perfect t0) | 85.4 ± 0.5% | 80.5% | 91.9% |
| jet vs jet, all timed ghost tracks | 76.6 ± 0.5% | 71.6% | 78.2% |

- **The jet vs t0 veto costs ~15% of the signal. The perfect t0 does no
  better than the real ones.**
  - That loss is set by the jet time, not by the t0. It tests every timed
    leg, and an HS leg's core time has a non-Gaussian tail: 7.9% of HS legs
    fail 3σ against the true t_HS even after calibration.
  - It pays that tail on up to two legs in the 93% of events with at least
    one timed leg.
- **Jet vs jet pays it only where both legs are timed**, in 24% of the
  (weighted) events.

By pair composition (default column; share = weighted fraction):

| composition | share | jet vs jet | TRKPTZ | WAVeS | HGTD | truth t0 |
|---|---:|---:|---:|---:|---:|---:|
| HS + HS | 78.8% | 96.7 | 90.5 | 91.9 | 90.5 | 91.4 |
| R1: fwd HS + fwd PU | 8.6% | 71.1 | 60.8 | 62.0 | 62.9 | 61.6 |
| R2: fwd PU + central HS | 6.1% | 97.1 | 56.1 | 57.4 | 71.0 | 54.8 |
| PU + PU, both forward | 1.5% | 92.5 | 58.1 | 57.0 | 73.1 | 57.0 |

- **R2 is where the two differ by construction.** Jet vs jet cannot act (its
  central HS leg has no time) and keeps 97%. The t0 rows reject it like a
  fake, ~45%. HGTD rejects less there because its t0 is valid less often.
- **R1 is kept at 71% by jet vs jet.** That is mostly because 59% (weighted) of its
  PU legs are past |η| 4.0, carry no tracks, and leave the pair untestable. The
  t0 rows lose the same R1 events and also test the HS leg.
- The raw-σ and ≥ 2-track columns are the same trade-off at other working
  points. They move signal efficiency and pileup rejection together. The
  `--min-trk 2` column keeps 83% of R1 against 71%, because single-track jet
  times carry much of the pileup rejection.

## What is still missing: the background

None of this says which veto is better. That needs Z+jets, and two warnings
apply to reading it:

1. **At a fixed 3σ the vetoes sit at different signal efficiencies.** Compare
   them at matched ε_S. `<out>_tradeoff.(pdf|png)` scans the threshold from
   1σ to 12σ and draws ε_B against ε_S for each method.
2. **Z+jets at pT(ll) > 200 GeV is thin.** It is ~200 MC events with JVT loose
   (204 pairs in the region stacks), a few percent statistical error inclusive
   and far more per column. `--zpt-min 0` or `100` gives higher-statistics
   references with a different pair composition.

The background is the Z→ℓℓ stand-in for Z→νν, and the signal assumes
B(H→inv) = 100%. The relative S/√B change (`vs no veto` in the summary
table, and `<out>_s_over_sqrtb`) depends on neither.

## Reproducing

```bash
# local (laptop): the 33-file local VBF sample, 112,400 events
cd build && make vbs_time_veto && ./vbs_time_veto --jvt=loose   # -> ../figs/hists/jvtLoose_vbs_time_veto.root
cd .. && PYTHONNOUSERSITE=1 PYTHONPATH=/opt/homebrew/Cellar/root/6.40.04/lib/root \
  ~/.venv-hgtd/bin/python python/vbs_time_veto_plot.py --sig figs/hists/jvtLoose_vbs_time_veto.root \
  --sig-label "#sqrt{s} = 14 TeV, HL-LHC, VBF H#rightarrowinv. (local)" --tag "JVT + fJVT loose" \
  --out figs/time_veto/local_jvtLoose

# AF: VBF from the RAW ntuples (truth MET), Z+jets from its skim; from condor/
condor_submit -a sample=vbf -a nshards=20 -a 'selargs=--jvt=loose' -a seltag=jvtLoose_ \
  -a 'inputargs=--ntuple-dir=/data/mcardiff/exotic_superntuples/highstats_vbf/' vbs_time_veto.sub
condor_submit -a sample=zjets -a nshards=8 -a 'selargs=--jvt=loose' -a seltag=jvtLoose_ vbs_time_veto.sub
# then, from the repo root (TTrees only: hadd, pinned shard counts)
hadd -f condor/vbf/vbf_jvtLoose_vbs_time_veto.root condor/vbf/vbf_jvtLoose_vbs_time_veto.shard*of20.root
hadd -f condor/zjets/zjets_jvtLoose_vbs_time_veto.root condor/zjets/zjets_jvtLoose_vbs_time_veto.shard*of8.root
python python/vbs_time_veto_plot.py --sig condor/vbf/vbf_jvtLoose_vbs_time_veto.root \
  --bkg condor/zjets/zjets_jvtLoose_vbs_time_veto.root --tag "JVT + fJVT loose" --out figs/time_veto/jvtLoose
```

The plotting script checks completeness. VBF must sum to 1,464,400 events and
the raw Σw of 6,497,726.94; Z+jets to 530,962 skim events. It refuses a
E_T^miss cut on a file that has no truth record.
