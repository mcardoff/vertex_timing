# Event-level timing veto on the VBS pair: jet vs jet, or jet vs t0

`util/vbs_time_veto.cxx` + `python/vbs_time_veto_plot.py`, 2026-09-28.

The grid samples ran as two condor clusters:
- VBF: raw ntuples (for the truth MET), cluster 3189030, 20 shards,
  ~11 min each;
- Z+jets: skim, cluster 3189031, 8 shards, ~3 min each.

Figures are in `figs/time_veto/jvtLoose_*`. The full tables are in
`figs/time_veto/jvtLoose_summary.md`.

## The question

This comes from the feedback of Friday 2026-09-25. In R1 (both tagging jets
forward, one HS and one PU), no vertex t0 is needed to reject the event: the
two jets can be tested against each other. The comparison asked for is signal
(VBF H→inv) and background (Z+jets) yields, weighted, after the analysis
recoil cut:

| method | veto the event when | can act on |
|---|---|---|
| **jet vs jet** | both legs are timed and \|t_A − t_B\| / √(σ_A² + σ_B²) ≥ 3 | pairs with both legs timed (R1 geometry) |
| **jet vs t0** | the event has a t0 and **any** timed leg has \|t_leg − t0\| / √(σ_leg² + (f σ_t0)²) ≥ 3 | every pair with a timed leg: R1, R2, others |

- **t0 rows.** Each jet vs t0 row uses a different t0: TRKPTZ, WAVeS, HGTD
  (Athena `RecoVtx_time`, where valid), TZP, and the truth HS time (a perfect
  t0, σ = 0).
- **f** is rpt_v5's per-sample t0 inflation (`inflationFor`, now in
  `util/rpt_v5_common.h`).
- **The jet time and its σ** are described in "The jet time" below.

**Selection:** the region plots' own, so every efficiency describes a
region-plot column (`vbs_region_diag` + `vbs_region_stack.py`):
- vertex quality; Z→ℓℓ and lepton-jet overlap removal (Z+jets only);
- ≥ 2 jets > 30 GeV, ≥ 1 of them in 2.38 < |η| < 4.0;
- **JVT + fJVT loose before pairing**;
- the max-m_jj opposite-hemisphere pair, m_jj ≥ 200 GeV, no |Δη| cut.

**Recoil cut:**
- VBF: truth E_T^miss = |Σ p_T^ν| > 200 GeV. The skims drop the truth record,
  hence the raw ntuples.
- Z+jets: pT(ll) > 200 GeV.

**Weights:** efficiencies are Σw(kept) / Σw with the generator weight.
Yields are at 3000 fb⁻¹ with the normalisation of
`rpt_region_distributions.md`.

## Answer

**Neither veto improves S/√B meaningfully after the recoil cut.**
- Jet vs jet lowers it.
- Jet vs t0 is within a few percent of no veto, which is at the edge of the
  Z+jets MC statistics.

| veto | ε_S (VBF) | ε_B (Z+jets) | S | B | S/√B vs no veto |
|---|---:|---:|---:|---:|---:|
| none | — | — | 383,009 | 1,251,702 ± 89,388 | 1 |
| **jet vs jet** | **95.1 ± 0.1%** | **96.0 ± 1.4%** | 364,413 | 1,201,649 | **0.971** |
| jet vs t0: TRKPTZ | 85.8 ± 0.1% | 70.5 ± 3.2% | 328,650 | 882,493 | 1.022 |
| jet vs t0: WAVeS | 87.3 ± 0.1% | 71.5 ± 3.2% | 334,345 | 895,002 | 1.032 |
| jet vs t0: HGTD (Athena) | 86.9 ± 0.1% | 78.5 ± 2.9% | 332,790 | 982,598 | 0.981 |
| jet vs t0: TZP | 86.6 ± 0.1% | 69.5 ± 3.3% | 331,682 | 869,984 | 1.039 |
| jet vs t0: truth (perfect t0) | 86.6 ± 0.1% | 69.0 ± 3.3% | 331,799 | 863,714 | 1.043 |
| jet vs jet, all timed ghost tracks | 76.1 ± 0.2% | 86.5 ± 2.4% | 291,373 | 1,082,734 | 0.818 |

- **The background rests on 204 Z+jets MC events.** That is ±3% (absolute)
  on ε_B, so each jet vs t0 S/√B ratio carries about ±2.4%. The rows share
  events, so their differences are better known than that.
- **Jet vs jet can only lose.**
  - At 3σ it costs 4.9% of the signal and removes 4.0% of the background.
  - Its trade-off curve saturates at ε_B ≈ 0.95 even at a 1σ threshold.
  - At no threshold does it raise S/√B above 1.
- **At matched signal efficiency, jet vs t0 beats jet vs jet.**
  - Tuned to jet vs jet's 3σ point (ε_S = 95.1%), the t0 vetoes keep 84.5%
    (truth t0), 86.0% (TZP), 86.5% (WAVeS), 88.0% (TRKPTZ) and 89.5% (HGTD)
    of the background, against 96.0% (`_tradeoff` figure).
  - Scanning the threshold, the best jet vs t0 gains are +2 to +5% in S/√B.
    Those optima are tuned on the same 204 events.
- **HGTD's t0 is the weakest t0**, because it is valid less often (the same
  pattern as in the R_pT study).

**S/√B vs no veto, by m_jj threshold** (the analysis cuts on m_jj; Z+jets MC
events in brackets):

| m_jj > | jet vs jet | TRKPTZ | WAVeS | HGTD | TZP | truth t0 |
|---|---:|---:|---:|---:|---:|---:|
| 200 GeV (204) | 0.971 | 1.022 | 1.032 | 0.981 | 1.039 | 1.043 |
| 500 GeV (163) | 0.974 | 1.055 | 1.063 | 0.998 | 1.070 | 1.086 |
| 1000 GeV (89) | 0.970 | 1.070 | 1.068 | 0.997 | 1.080 | 1.090 |
| 1500 GeV (42) | 0.937 | 0.986 | 0.984 | 0.937 | 0.994 | 1.027 |
| 2000 GeV (21) | 0.925 | 0.913 | 0.923 | 0.901 | 0.920 | 0.976 |

Above 1.5 TeV the background has too few events to say anything.

## Why jet vs jet cannot reach the background

After pT(ll) > 200 the Z recoils against a real jet, so the background is
mostly pairs with a central leg:

| weighted share | R1 | R2 | HS + HS | PU + PU, both fwd | neither label | both legs timed |
|---|---:|---:|---:|---:|---:|---:|
| VBF | 8.6% | 6.5% | 77.6% | 1.7% | 3.3% | 23.8% |
| Z+jets | 7.5% | **34.5%** | 35.0% | 5.5% | 12.0% | **6.5%** |

- **R2 is the largest fake class.** Its HS leg is central and has no time,
  so jet vs jet keeps 98.6% of it.
- **Only 6.5% of background events have two timed legs** to compare, against
  23.8% of the signal. In the signal, those events include every genuine
  HS+HS forward pair whose jet time has a tail.

Efficiency (%) by composition (MC events in brackets for Z+jets):

| composition | VBF: jet vs jet | VBF: TRKPTZ / WAVeS / HGTD / truth | Z+jets: jet vs jet | Z+jets: TRKPTZ / WAVeS / HGTD / truth |
|---|---:|---:|---:|---:|
| HS + HS | 97.2 | 91.0 / 92.7 / 90.6 / 92.3 | 100.0 (72) | 90.0 / 92.9 / 90.0 / 92.9 |
| R1 | 76.6 | 66.9 / 68.2 / 69.2 / 67.2 | 73.3 (15) | 53.3 / 60.0 / 60.0 / 60.0 |
| R2 | 97.4 | 60.1 / 59.9 / 73.3 / 57.8 | 98.6 (69) | 53.6 / 55.1 / 72.5 / 52.2 |
| PU + PU, both fwd | 93.5 | 63.4 / 62.6 / 73.7 / 60.5 | 81.8 (11) | 54.5 / 54.5 / 63.6 / 27.3 |
| neither label | 90.8 | 72.1 / 73.1 / 76.0 / 71.5 | 95.8 (26) | 70.8 / 62.5 / 75.0 / 58.3 |

## Why jet vs t0 gains so little

- **It removes the fakes**: a third to a half of R1 and R2 (HGTD less).
  But **it also removes 7–10% of genuine HS+HS pairs, in signal and background
  alike**.
  - In VBF, HS+HS is 78% of the events. In Z+jets after the recoil cut it is
    35%, and the fakes are not much more than half.
  - So the veto keeps ~87% of signal against ~70% of background.
- **The perfect t0 is no better than TZP or WAVeS.** The HS+HS loss is set by
  the jet time, not the t0.
  - Even after calibration, 7.1% of grid-VBF HS legs fail 3σ against the
    TRUE t_HS: pileup left in the time cluster, or a pileup cluster that
    outweighs the HS one.
  - Single-track jet times are 14% of HS legs, fail 24% of the time, and
    cause half of those failures. Jets with four or more tracks fail ≤ 2%.
  - The veto pays that tail on every timed leg, in 93% of the signal
    events.
- **A better jet time is the lever**, not a better t0.

## Robustness

Changing one convention at a time; S/√B relative to no veto (ε_S / ε_B for
the jet vs t0 rows in brackets):

| convention | jet vs jet | TRKPTZ | WAVeS | HGTD | TZP | truth t0 |
|---|---:|---:|---:|---:|---:|---:|
| default: σ_jet × f(n) from VBF, both samples | 0.971 | 1.022 (85.8 / 70.5) | 1.032 | 0.981 | 1.039 | 1.043 (86.6 / 69.0) |
| one global factor from VBF (1.507), both | 0.966 | 1.039 | 1.049 | 0.985 | 1.052 | 1.061 |
| each sample its own global factor (1.51 / 1.83) | 0.961 | 1.009 | 1.023 | 0.979 | 1.029 | 1.037 |
| raw σ (no calibration) | 0.953 | 1.001 (80.7 / 65.0) | 1.036 | 0.954 | 1.027 | 1.030 |
| ≥ 2 tracks for a jet time | 0.985 | 1.024 (90.8 / 78.5) | 1.032 | 0.988 | 1.036 | 1.041 |

Jet vs jet is below 1 in every row. Jet vs t0 stays within −5% to +6%.

**A looser Z+jets recoil cut gives far more statistics, but a different
background.** It swaps in the pileup pairs that the recoil cut removes, and
the signal side stays at MET > 200:

| Z+jets cut | MC events | ε_B: jet vs jet | ε_B: TRKPTZ / WAVeS / HGTD / TZP / truth |
|---|---:|---:|---:|
| pT(ll) > 200 (default) | 204 | 96.0% | 70.5 / 71.5 / 78.5 / 69.5 / 69.0% |
| pT(ll) > 100 | 1,372 | 93.8% | 64.3 / 65.7 / 75.6 / 64.7 / 64.1% |
| none | 26,883 | 92.0% | 50.9 / 51.9 / 71.8 / 50.9 / 47.2% |

The t0 veto removes half of the uncut Z+jets pairs, which are mostly pileup.
It removes 30% after the analysis recoil cut. **The recoil cut takes away
most of what timing is good at.**

## Second method: remove the failing jets and re-pair (grid re-run pending)

**Why.** With both recoil cuts at 100 GeV (`figs/time_veto/jvtLoose_recoil100_*`,
1,372 Z+jets MC events), jet vs jet still lowers S/√B (×0.97), while the jet
vs t0 vetoes gain 4–6%. The event veto has a structural problem:
- R1 is as common in the signal (10.7%) as in Z+jets (9.9%). In VBF, the
  max-m_jj pair picks up a forward pileup jet just as often.
- Rejecting the EVENT therefore throws away signal whose genuine VBF jets are
  still there.

**The method.** Remove the jets that fail the timing test, then re-form the
pair from the remaining jets, as JVT/fJVT do.
- `util/vbs_time_veto` now also writes every jet the pair was chosen from:
  the pT-passing, post-JVT jets, with kinematics, labels and core jet time
  (`jet_*` arrays).
- The plotting script re-pairs offline:
  - `rp_<src>` removes every timed jet incompatible with that t0 at 3σ. It
    then requires ≥ 2 jets, the max-m_jj opposite-hemisphere pair, and
    m_jj ≥ 200.
  - `rp_jj` removes nothing; it takes the max-m_jj pair whose two jets are
    time-compatible.
- **The forward-jet preselection is not re-applied after removal.** The
  region plots' "≥ 1 jet in 2.38 < |η| < 4.0" defines the starting sample.
  It is not an analysis cut, and re-applying it loses 4.1% of the signal even
  though a valid pair survives (`--repair-reapply-fwd` re-applies it).
- **Plot binning.** Re-paired events move in m_jj, so efficiencies are binned
  in each event's no-veto m_jj. Yields and S/√B are binned in the final
  pair's m_jj.

**Validation (local VBF).**
- The pair re-formed from the arrays equals the stored pair in all events.
- At an infinite threshold every re-pair variant reproduces the no-veto
  yields exactly.
- The event-veto tables are unchanged, identical on the grid files.

**Local VBF signal, truth MET > 100** (29,070 MC events):

| t0 | event veto | jets removed, pair re-formed |
|---|---:|---:|
| TRKPTZ | 83.7% | 90.2% |
| WAVeS | 85.6% | 91.9% |
| HGTD (Athena) | 85.6% | 91.2% |
| truth | 84.7% | 91.3% |

Jet removal recovers about 6.5 points of signal, mostly in R1: 80–83% kept,
against 63–67%. What it still loses is structural:
- In ~6% of events, a genuine VBF jet fails the test and no other jet is left
  to pair with.
- ~1% have no opposite-hemisphere pair.

The background needs the grid re-run, since the current grid files predate
the jet arrays.

## The jet time

"The ghost-associated tracks with valid time inside each jet", averaged as they
are, do not measure the jet's time at μ = 200. Local VBF, forward jets above
30 GeV:
- A forward jet's ghost area holds **~31 tracks, ~23 of them timed**.
- **Only 21% of an HS jet's timed ghost tracks come from the hard scatter.**
  The rest, ~24 tracks per jet, is pileup that falls inside the jet's
  catchment area.

So their inverse-variance mean (`all_`, kept in the tree for the record) is a
pileup average:
- χ²/ndf ≈ 65;
- the core spread about truth is 23× its quoted σ;
- 11% of HS jets are within 3σ of t_HS;
- 25% of HS+HS pairs pass the jet vs jet test, against 21% of R1 pairs.

It cannot tell them apart. On the grid it is the last row of the answer table:
it keeps 76% of the signal.

The default (`core_`) estimator:
- keeps timed ghost tracks with the standard kinematics (1 < pT < 30 GeV,
  quality flag) within ΔR < 0.2 of the axis (the R_pT cone);
- clusters them in time (iterative, 2σ);
- takes the inverse-variance mean of the highest-ΣpT cluster.

It raises the HS share of the tracks used to ~73%, and puts 87.7% of HS legs
within 3σ of t_HS (raw σ, local VBF). It was chosen on local-VBF truth (19,547
events, both legs in 2.4–4.0, raw σ). HS+HS and R1 pass rates at 3σ, one
ingredient changed at a time:

| variant | HS+HS kept | R1 kept |
|---|---:|---:|
| **ΔR 0.2, 2σ, highest-ΣpT cluster (default)** | **83.3%** | 20.2% |
| ΔR 0.4 | 77.1% | 19.0% |
| 3σ clustering | 81.6% | 19.4% |
| cluster with the most tracks | 82.3% | 18.1% |
| no clustering (ΔR 0.2, pT > 1 GeV, plain mean) | 39.4% | 23.0% |
| no clustering, all timed ghost tracks | 24.7% | 21.3% |

## The jet-time calibration

The core estimator's quoted σ is a statistical error. It knows nothing of the
pileup left in the cluster or of mis-assigned hits, so it is too small, as the
vertex-t0 σ is.
- The Gaussian-core width of the HS-leg pull |t_jet − t_HS| / σ (1.4826 ×
  median) is **1.507 on grid VBF (526,456 legs) and 1.833 on Z+jets (3,830)**.
- **It is mostly a function of the track count n of the estimate.** On VBF it
  is 2.34 for n = 1, then 1.64, 1.46, 1.37 (n = 4–5), 1.34 (n = 6–8) and 1.31
  (n ≥ 9).
- Z+jets HS jets are softer: 22% single-track against 14%. That accounts for
  most of its larger global width.
- With VBF's n-binned factors, Z+jets' width drops to 1.15. A residual is
  left: at fixed n, its HS jets still time worse.

**Default: σ_jet × f(n), measured once on VBF, applied to both samples.**
- The veto is then one procedure for signal and background, as an analysis
  would apply it.
- Per-sample calibration would give the background a looser veto than the
  signal, from truth the analysis does not have.
- The residual width is a real property of the veto: softer HS jets fail it
  more often. Z+jets' own calibration is the "each sample its own" robustness
  row, and moves nothing.
- The t0 inflation stays rpt_v5's per-sample table. Using VBF's values for
  Z+jets too changes ε_B by ≤ 0.5 points.

## Validation (all pass)

- **Completeness.**
  - VBF read all 1,464,400 raw events, and its summed weight equals the raw
    sample's Σw exactly (6,497,726.944), the normalisation denominator.
  - Z+jets read all 530,962 skim events.
  - All 28 jobs exited 0 with empty stderr.
  - The hadd'ed rows equal the shard sums (504,362 / 26,883) with no
    duplicate events.
- **Event set.**
  - Z+jets selects 26,883 events, and 204 pass pT(ll) > 200. Both equal
    the JVT-loose Z+jets `vbs_region_diag` run of 2026-09-18, the region
    stacks' input.
  - On local VBF the selected set equals `vbs_region_diag --jvt=loose`'s
    event by event (38,524; legs, m_jj and labels identical).
- **t0 against rpt_v5_hist's JVT-loose region trees.**
  - Every region event passing the composition preselection is found: 91,624
    VBF and 5,528 Z+jets.
  - t0 and validity are identical for hgtd, trkptz, waves and tzp (Δt0 = 0).
  - None of the region rows that fail the preselection is selected.
  - The VBF match is by m_jj and leg pT, since this run read the raw files
    and the tree the skims. So raw and skim give identical t0s.
- **Recoil variables against last week's friend files.**
  - VBF truth MET equals `vbf_jvtLoose_truth_met.root` (key lookup into the
    raw ntuples) to 6×10⁻⁵ GeV.
  - pT(ll) equals `zjets_jvtLoose_zpt.root`.

## Caveats

- **The background is Z→ℓℓ** with pT(ll) as the stand-in for Z→νν E_T^miss.
  The signal is at B(H→inv) = 100%.
- The relative S/√B columns depend on neither.
- **B is only 204 MC events** after the recoil cut.

## Figures

`figs/time_veto/`:
- `jvtLoose_sig_eff_mjj`, `jvtLoose_bkg_eff_mjj`: efficiency per m_jj column,
  region-plot layout, MC events underneath.
- `jvtLoose_s_over_sqrtb`: ε_S / √ε_B per column. That is S/√B divided by
  its no-veto value, so it cancels the two normalisation stand-ins
  (B(H→inv) = 100%, Z→ℓℓ for Z→νν).
- `jvtLoose_s_over_sqrtb_yields`: S/√B from the yields at 3000 fb⁻¹ per
  column, with the no-veto value as its own marker. `--bf-hinv` scales the
  signal to a realistic branching ratio.
- `jvtLoose_tradeoff`: ε_B against ε_S, threshold scanned from 1σ to 12σ.
- `local_jvtLoose_*`: the local-VBF development run.

## Reproducing

```bash
# AF: VBF from the RAW ntuples (truth MET), Z+jets from its skim; from condor/
condor_submit -a sample=vbf -a nshards=20 -a 'selargs=--jvt=loose' -a seltag=jvtLoose_ \
  -a 'inputargs=--ntuple-dir=/data/mcardiff/exotic_superntuples/highstats_vbf/' vbs_time_veto.sub
condor_submit -a sample=zjets -a nshards=8 -a 'selargs=--jvt=loose' -a seltag=jvtLoose_ vbs_time_veto.sub
# from the repo root (TTrees only: hadd, pinned shard counts)
hadd -f condor/vbf/vbf_jvtLoose_vbs_time_veto.root condor/vbf/vbf_jvtLoose_vbs_time_veto.shard*of20.root
hadd -f condor/zjets/zjets_jvtLoose_vbs_time_veto.root condor/zjets/zjets_jvtLoose_vbs_time_veto.shard*of8.root

# laptop, after scp of the two merged files into condor/{vbf,zjets}/
PYTHONNOUSERSITE=1 PYTHONPATH=/opt/homebrew/Cellar/root/6.40.04/lib/root \
~/.venv-hgtd/bin/python python/vbs_time_veto_plot.py \
  --sig condor/vbf/vbf_jvtLoose_vbs_time_veto.root --bkg condor/zjets/zjets_jvtLoose_vbs_time_veto.root \
  --tag "JVT + fJVT loose" --out figs/time_veto/jvtLoose
# robustness: --jet-infl sig | per-file | 1, --min-trk 2, --zpt-min 100 / 0
# jet removal + re-pairing against the event veto, recoil cuts at 100 GeV:
#   ... --met-min 100 --zpt-min 100 --methods rp_trkptz,rp_waves,rp_hgtd,rp_truth \
#       --compare trkptz,waves,hgtd,truth --out figs/time_veto/jvtLoose_recoil100_repair

# local VBF development run: 33 files, 112,400 events
cd build && ./vbs_time_veto --jvt=loose        # -> ../figs/hists/jvtLoose_vbs_time_veto.root
```

The plotting script checks completeness:
- VBF must sum to 1,464,400 events and the raw Σw;
- Z+jets to 530,962 skim events.

It refuses an E_T^miss cut on a file with no truth record.
