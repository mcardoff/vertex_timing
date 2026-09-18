# The VBS pair with pileup-jet tagging applied first

`util/vbs_region_diag.cxx --jvt=none|loose|tight`, 2026-09-18: local VBF
(33 files, 112,400 events) and the skimmed grid Z+jets sample (80 files,
530,962 events, condor clusters 3050292-94). Follows
`vbs_pair_composition.md`, on Ariel's point that an analysis applies JVT and
fJVT to the jets BEFORE it forms a tagging pair, so a composition measured on
untagged jets describes a population no analysis sees.

**Headline, Z+jets: the standard taggers do not turn Z+jets into an R2
sample. They remove the CENTRAL pileup and leave the forward pileup, so the
surviving pairs are forward-forward pileup.** Under the loose point 59% of
selected events are dropped outright (77% tight); the F-PU + C-PU cell
collapses from 37.8% to 14.5% to 5.7% of pairs; but F-PU + F-PU GROWS as a
share, 27.8 → 38.8 → 42.0%, because fJVT at μ = 200 keeps two thirds of
forward pileup jets (PU efficiency 67% loose, 47% tight) while the central
JVT keeps 15% / 1.8%. R2 does rise, 7.8 → 15.1 → 20.3% of pairs, but by
removing its competition rather than by creation: 78% of post-loose R2
events were R2 before the tagger, and the absolute R2 count falls 5,064 →
4,055 → 3,143. Even tight leaves 60% of Z+jets pairs with no hard-scatter
leg at all. The previous write-up's conclusion — the central HS leg is
mostly not there — stands after tagging, with a second reason on top: the
forward fake cannot be tagged away at this pileup.

**VBF, the control: the taggers remove pileup and leave R1/R2 alone.** What
they remove is 98.5% paper-PU on the JVT side and 86% on the fJVT side; the
HS+HS share of the chosen pair rises from 56% to 71% (loose) and 78%
(tight); pileup-only pairs fall from 5.2% to 2.8% and 1.9%. R1 and R2 shrink
only modestly (11.9 → 10.9 → 9.4%, 6.9 → 6.7 → 6.1%), and event by event 70%
of R1 events and 67% of R2 events are still R1/R2 after the loose point.

dijet has not been run (same three submissions with `sample=dijet`).

## What is computed, and what it is not

The SuperNtuples carry no Jvt/fJvt decoration (241 branches, none match, on
either the EMTopo or EMPFlow collection), so both are computed in the
diagnostic from the vertex fit's own track assignment (`Track_recoVtx_idx`,
which is -1 for 47% of tracks and 0 for the primary vertex) and each jet's
ghost-associated tracks:

| quantity | definition | source |
|---|---|---|
| R_pT | Σ pT(ghost tracks fitted to vertex 0) / pT_jet | arXiv:1510.03823 |
| corrJVF | pT_PV / (pT_PV + pT_PU / (k n_PU^trk)), k = 0.01 | stored, not cut on |
| fJVT | max over pileup vertices i of (p_T^miss,i · ĵ_T) / pT_jet | arXiv:1705.02211 |
| p_T^miss,i | −½ (Σ pT of central tracks fitted to i + Σ pT of central jets whose dominant ghost vertex is i) | |

**JVT proper is a k-NN likelihood over (corrJVF, R_pT) and cannot be rebuilt
from the ntuple**, so the "JVT" cut is R_pT alone. Its thresholds were set by
`python/vbs_jvt_calibrate.py` on the `jets` tree of a `--jvt=none` run so
that the paper-HS efficiency of central jets in the JVT window (|η| < 2.5,
30 < pT < 60 GeV) matches the published working-point efficiencies. fJVT
thresholds are the published ones. Every threshold and window, in one place:

| | JVT (R_pT proxy) | fJVT |
|---|---|---|
| window | \|η\| < 2.5, pT < 60 GeV | 2.5 ≤ \|η\| < 4.5, pT < 120 GeV |
| loose | R_pT ≥ 0.0167  (↔ Run-2 "Medium", 92% HS eff) | fJVT ≤ 0.5  (published Loose) |
| tight | R_pT ≥ 0.0767  (↔ Run-2 "Tight", 85% HS eff) | fJVT ≤ 0.4  (published Tight) |

A jet outside both windows is never removed. A jet with no ghost tracks has
R_pT = 0 and fails JVT inside its window (none occur above 30 GeV centrally).
Overlap-removed jets are skipped only in the p_T^miss jet assignment.

**Two things to be honest about.** The loose R_pT threshold, 0.0167, sits
just above the 8% of paper-HS central jets that have NO primary-vertex ghost
track at all, so "loose JVT" here is close to "has any PV track". Those 8% are
mostly label noise: a paper-HS label only needs a truth HS jet above 10 GeV
within ΔR 0.3, so a 30 GeV pileup jet sitting on a 10 GeV HS jet is labelled
HS. And fJVT is weaker at μ = 200 than in its paper: the max over ~100 pileup
vertices of a noisy projection inflates the HS jets' values (median 0.20,
90th percentile 0.44, where the paper's HS distribution peaks near zero), so
the published thresholds give HS efficiencies of 93.5% / 86.4% against
pileup efficiencies of 68% / 49% in the forward window. It is the standard
tool applied as published, not re-tuned; the `jets` tree carries every jet's
discriminants so the thresholds can be scanned offline.

Both discriminants separate in the right direction: over random (HS, PU)
pairs in their windows, R_pT is higher for the HS jet 91.5% of the time and
fJVT is higher for the PU jet 78.4% of the time.

## What the working points do (VBF, wide_ block, m_jj ≥ 200 GeV)

| | no tagger | loose | tight |
|---|---:|---:|---:|
| selected events | 48,640 | 38,673 | 33,623 |
| pairs in window | 48,323 | 38,524 | 33,525 |
| R1 (F-HS + F-PU) | 11.94% | 10.86% | 9.39% |
| R2 (F-PU + C-HS) | 6.93% | 6.69% | 6.07% |
| no HS leg | 6.45% | 3.56% | 2.39% |
| one HS leg | 37.50% | 25.66% | 19.98% |
| two HS legs | 56.05% | 70.78% | 77.64% |
| PU+PU, any region | 5.17% | 2.77% | 1.88% |
| both legs forward | 33.80% | 34.80% | 33.42% |

Jets removed, over the selected events (HS/PU by paper label; the remainder
are "neither"):

| | loose | tight |
|---|---:|---:|
| by JVT | 44,442 (HS 635, PU 40,616) | 43,799 (HS 972, PU 38,619) |
| by fJVT | 6,734 (HS 851, PU 5,128) | 9,361 (HS 1,219, PU 7,107) |
| events dropped below the jet requirements | 8,788 | 14,236 |

At the VBS-topology selection (m_jj ≥ 500 GeV, |Δη| ≥ 2.5) the picture is
the same: R1 14.5 → 12.8 → 10.9%, R2 6.6 → 6.1 → 5.4%, two-HS 59.8 → 71.7
→ 77.5%, no-HS 6.3 → 3.8 → 2.6%.

The 21-cell figure: `figs/local_vbs_jvt.png` (m_jj ≥ 200) and
`figs/local_vbs_jvt_mjj500.png`, from `python/vbs_jvt_plot.py`.

## Where the R1/R2 events go (event by event, m_jj ≥ 200)

`python/vbs_jvt_migration.py` joins the none and tagged outputs on
(file, entry). Rows are the region under no tagger; columns, where the same
event lands after it.

none → loose:

| | R1 | R2 | both HS | other | dropped | n |
|---|---:|---:|---:|---:|---:|---:|
| R1 | **69.7%** | 0.4% | 13.2% | 3.5% | 13.2% | 5,768 |
| R2 | 0.0% | **66.5%** | 6.8% | 0.5% | 26.1% | 3,350 |
| both HS | 0.1% | 0.1% | 92.6% | 0.5% | 6.7% | 27,087 |
| other | 1.2% | 2.5% | 9.9% | 34.1% | 52.0% | 12,118 |

none → tight: R1 stays R1 51.1% (20.3% → both HS, 25.4% dropped); R2 stays
R2 48.1% (10.3% → both HS, 41.0% dropped); "other" is dropped 67.1%.

Two readings:

- **R1 → both-HS is the fJVT doing its job on the pileup leg.** The forward
  pileup jet is tagged, the event still has a second HS jet, and the picker
  re-pairs to it. The R1 events that DROP are the ones where the tagged leg
  was the second jet; they leave the selection rather than becoming a
  different region.
- **R2 drops more than R1 (26% vs 13%) and the difference is the central
  jet count.** An R2 event is a forward fake plus one central HS jet; tag the
  fake and there is usually nothing left to pair. That is the same statement
  as the previous study's "R2's scarcity is the central hard-scatter leg
  mostly not being there", seen from the other side.
- **The post-tagger R1/R2 are almost entirely the pre-tagger R1/R2** (96% and
  86% under loose). The taggers do not manufacture new R1/R2 events out of
  "Other"; they thin the existing ones.

## Z+jets (skimmed grid sample, wide_ block, m_jj ≥ 200 GeV)

The no-tagger column reproduces the previous write-up's 64,950 pairs exactly,
so the skim and the earlier raw-sample run agree.

| | no tagger | loose | tight |
|---|---:|---:|---:|
| selected events | 69,701 | 28,293 | 16,114 |
| pairs in window | 64,950 | 26,883 | 15,517 |
| F-PU + C-PU | 37.8% | 14.5% | 5.7% |
| F-PU + F-PU | 27.8% | **38.8%** | **42.0%** |
| R2 (F-PU + C-HS) | 7.8% | 15.1% | 20.3% |
| R1 (F-HS + F-PU) | 3.3% | 5.5% | 6.8% |
| F-HS + C-HS | 1.7% | 4.2% | 6.8% |
| no HS leg | 79.9% | 68.6% | 60.3% |
| two HS legs | 2.0% | 4.8% | 7.6% |
| both legs forward | 36.3% | **51.8%** | **57.2%** |

Jets removed (HS/PU by paper label): loose, JVT 37,294 (HS 479, PU 33,549),
fJVT 14,080 (HS 291, PU 12,679); tight, JVT 24,484 (HS 542, PU 21,064), fJVT
13,725 (HS 353, PU 12,305). The tight JVT count is LOWER than loose because
the tight run drops more events before the count is taken.

**Why the composition tilts forward.** The two taggers are not equally
effective, and the gap is larger on Z+jets than on VBF:

| Z+jets, fixed thresholds | HS eff | PU eff |
|---|---:|---:|
| JVT loose (R_pT ≥ 0.0167, central) | 90.8% | 15.3% |
| JVT tight (R_pT ≥ 0.0767) | 80.9% | 1.8% |
| fJVT loose (≤ 0.5, forward) | 84.1% | 66.8% |
| fJVT tight (≤ 0.4) | 70.3% | 47.4% |

So a central pileup jet is removed 85–98% of the time and a forward one
33–53% of the time. Every cell with a central PU leg drains (F-PU + C-PU
37.8 → 5.7%, F-HS + C-PU 4.8 → 1.1%, F-X + C-PU 3.2 → 0.6%) and every cell
whose pileup is forward-only holds or grows as a share (F-PU + F-PU, F-PU +
F-X, R1). "Both legs forward" goes from 36% to 57% of pairs. The thresholds
are the VBF-calibrated ones held fixed; on Z+jets they sit ~1 point (loose)
and ~4 points (tight) below their VBF hard-scatter efficiencies, and the
fJVT forward HS efficiency is 9–16 points below VBF's, because Z+jets
forward HS jets are softer and their fJVT median is 0.30 against VBF's 0.20.

**Migration, none → loose** (rows: region before; columns: after):

| | R1 | R2 | both HS | other | dropped | n |
|---|---:|---:|---:|---:|---:|---:|
| R1 | 63.3% | 0.8% | 2.7% | 6.5% | 25.9% | 2,165 |
| R2 | 0.0% | 62.8% | 0.4% | 0.6% | 35.9% | 5,064 |
| both HS | 0.0% | 0.3% | 83.6% | 0.4% | 15.6% | 1,279 |
| other | 0.2% | 1.5% | 0.3% | 35.2% | 61.9% | 56,442 |

Post-loose R2 is 78% pre-tagger R2 and 21% pre-tagger "other" (a pileup pair
whose central PU leg was tagged and re-paired to a central HS jet); post-loose
R1 is 93% pre-tagger R1. Under tight, R1 stays R1 43% and R2 stays R2 45%,
with 47% / 54% dropped. The taggers thin R1/R2 on Z+jets faster than on VBF
(26% / 36% dropped against 13% / 26%) because a Z+jets R1/R2 event rarely has
a spare jet to re-pair to.

The original two-sample figure (ranked bars + 6×6 matrices, m_jj ≥ 500,
|Δη| ≥ 2.5) rebuilt on the tagged data:
`figs/vbs_pair_composition_jvtLoose.png` / `_jvtTight.png`, from
`python/vbs_pair_composition_plot.py`. Under loose, F-PU + F-PU is 57.4% of
Z+jets pairs on its own, against VBF's 47.2% in F-HS + C-HS.

The seven-band stack vs m_jj (ROOT/ATLAS style, |Δη| > 2.5, columns
normalised to 1) for each working point:
`condor/zjets/zjets[_jvtLoose|_jvtTight]_vbs_region_stack.png` and
`figs/local[...]_vbs_region_stack.png`, from `python/vbs_region_stack.py`.
On Z+jets the "PU + PU, both forward" band is the only one that grows with
tagging and with m_jj; above 2 TeV it is ~85% of every column.

Figures: `condor/zjets/zjets_vbs_jvt.png` (m_jj ≥ 200) and
`condor/zjets/zjets_vbs_jvt_mjj500.png` (m_jj ≥ 500, |Δη| ≥ 2.5, where R1 /
R2 go 5.8 → 8.0 → 9.2% / 5.5 → 8.5 → 10.4% and no-HS-leg 83.4 → 77.4 →
73.3%).

## What this does NOT settle

- **dijet.** Not run. Same three submissions with `sample=dijet`, then
  `python/vbs_jvt_plot.py --dir condor/dijet --sample dijet`.
- **Whether fJVT's weakness here is the proxy or the pileup.** The fJVT
  computed is the published definition, but its inputs are the vertex fit's
  track assignment at μ = 200 with ~100 vertices, and the max over that many
  vertices inflates every jet's value. A real Athena fJVT on one file would
  say whether 67% PU efficiency at the loose cut is what the tool does at
  this pileup or an artefact of the reconstruction here. Until then the
  Z+jets forward-forward residual should be read as "not removable by fJVT
  as computed here", not as a property of the sample.
- **The JVT proxy is not JVT.** If the numbers are going to be quoted as
  "with JVT", the ntuple production should be asked for the Jvt/fJvt
  decorations (or NNJvt for Run 3), and the proxies validated against them on
  one file. The `jets` tree has everything needed for that comparison.
- **Thresholds are calibrated on local VBF.** A different sample's HS jets
  have a different R_pT distribution; the thresholds are held FIXED across
  samples on purpose (a working point is a cut, not an efficiency), but the
  per-sample HS efficiency should be read off the jets tree, not assumed.

## Reproduction

```bash
cd build && make vbs_region_diag
./vbs_region_diag                 # --jvt=none, writes ../figs/local_vbs_region_diag.root
./vbs_region_diag --jvt=loose     # ../figs/local_jvtLoose_vbs_region_diag.root
./vbs_region_diag --jvt=tight     # ../figs/local_jvtTight_vbs_region_diag.root
cd .. && PYTHONNOUSERSITE=1 ~/.venv-hgtd/bin/python python/vbs_jvt_calibrate.py figs/local_vbs_region_diag.root
PYTHONNOUSERSITE=1 ~/.venv-hgtd/bin/python python/vbs_jvt_plot.py --dir figs --sample local --label "VBF (local)"
PYTHONNOUSERSITE=1 ~/.venv-hgtd/bin/python python/vbs_jvt_migration.py --dir figs --sample local
```

Each local run is ~25 s; a Z+jets condor job is ~19 min single-threaded.
The `--jvt=none` output was verified bit-identical to the pre-JVT diagnostic
on all 95 pre-existing columns (48,640 rows), so nothing in the earlier
write-up moved.

**A trap that cost one submission.** The first three Z+jets jobs wrote 10 KB
files with exit code 0 and zero events read: with `EXTENDED_BRANCHES` on, the
wrapper bound `Track_btagIp_*`, which local VBF has and every grid sample
lacks, and a TTreeReader with a missing branch iterates nothing rather than
failing. `recordAvailableBranches` (shared with the exporter) now runs before
binding, and the diagnostic exits 2 on zero events.
