# The VBS pair with pileup-jet tagging applied first

`util/vbs_region_diag.cxx --jvt=none|loose|tight`, local VBF (33 files,
112,400 events), 2026-09-18. Follows `vbs_pair_composition.md`, on Ariel's
point that an analysis applies JVT and fJVT to the jets BEFORE it forms a
tagging pair, so a composition measured on untagged jets describes a
population no analysis sees.

**Headline (VBF): the taggers remove pileup, not R1/R2.** What they remove is
98.5% paper-PU on the JVT side and 86% on the fJVT side; the HS+HS share of
the chosen pair rises from 56% to 71% (loose) and 78% (tight); and the
pileup-only pairs fall from 5.2% to 2.8% and 1.9%. R1 and R2 shrink only
modestly (11.9 → 10.9 → 9.4%, 6.9 → 6.7 → 6.1%), and event by event 70% of the
R1 events and 67% of the R2 events are still R1/R2 after the loose point.
The R1/R2 population Ariel asked about survives the taggers; it is the
pileup-only "Other" that goes.

**This is VBF only.** The whole reason the previous study exists is that
Z+jets is a PU+PU plot (70.7% of its pairs), and that is where the taggers
should act hardest. The Z+jets and dijet runs need condor
(`condor/vbs_region_diag.sub`) and have not been submitted.

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

## What this does NOT settle

- **Z+jets.** Everything above is VBF, where 56% of pairs were already HS+HS
  and the taggers had little to remove. On Z+jets, 79.7% of pairs had no HS
  leg; the loose point should remove most of the 65.6% in the two pileup-only
  cells, and the question that matters is what fraction of the SURVIVING pairs
  is R2. Submit all three working points per sample:
  `condor_submit -a sample=zjets [-a jvt=loose -a jvttag=_jvtLoose] vbs_region_diag.sub`,
  then `python/vbs_jvt_plot.py --dir condor/zjets --sample zjets`.
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

Each local run is ~25 s. The `--jvt=none` output was verified bit-identical
to the pre-JVT diagnostic on all 95 pre-existing columns (48,640 rows), so
nothing in the earlier write-up moved.
