# R_pT in the VBS regions: z-only against timed, per timing case

`util/rpt_v5_hist.cxx` region tree + `python/rpt_region_dists.py`,
2026-09-25. Grid VBF (skim: 1,187,999 events passing the vertex cut, of
1,464,400 generated; condor cluster 3143674, 8 shards) and grid Z+jets
(skim: 530,962 events passing vertex + Z->ll, of 2,000,000; cluster
3143675, 14 shards), default selection (max-m_jj pair, m_jj >= 200 GeV, no
|Deta|, no JVT/fJVT), GATE_SIGMA = 3.0.

The regions are shown as **distributions** — the same jets' R_pT without
timing (ITk-only, z association only) and with it — instead of ROCs, for five
timing cases. Legs past |eta| 3.8 are included. See "Definitions" below.

**The presentation version applies JVT + fJVT (loose) before the VBS pair is
formed**, as an analysis would. See "With JVT + fJVT loose" below, clusters
3144162/3. It removes 27% of VBF R1 events and 16% of R2, but leaves every
conclusion below unchanged to within a point. The untagged numbers that follow
remain the tagger-free reference.

**Headline, VBF.**
- **Timeable core (2.4 < |eta| < 3.8):**
  - With z only, 83.7% of forward PU legs already have R_pT = 0.
  - Every timing case raises that to 88–91% (both regions). Of the PU legs
    that had any R_pT, timing zeroes 34–43% in R1.
  - The price is on the HS leg: 27–30% of HS legs lose some R_pT, and 0.5–1.4%
    lose all of it. Mean HS-leg R_pT goes from 0.374 to 0.353–0.360.
- **Legs past |eta| 4.0 carry no tracks:**
  - 86% of such HS legs and 99% of such PU legs have R_pT = 0 with z only, and
    no case changes that.
  - These legs are **16% of R1 HS legs, 58% of R1 PU legs and 57% of R2 PU
    legs**.
  - Including them therefore mostly adds an untouchable spike at R_pT = 0.
    The |eta|-band page separates it out.

**R2 shows the known WAVeS failure, now in distribution form.** In R2 the only
forward jet is the pileup one, and WAVeS up-weights exactly its tracks. Of the
forward PU legs with z-only R_pT > 0, share zeroed:

| | TRKPTZ | WAVeS | HGTD t0 |
|---|---:|---:|---:|
| VBF | 40% | 26% | 28% |
| Z+jets | 39% | 20% | 25% |

In R1, where both legs are forward, WAVeS is level with TRKPTZ and HGTD on
suppression (41% vs 43% / 40%), and it is the gentlest on the HS leg
(26.6% lowered and 0.5% zeroed, against 29.6% and 1.4% for TRKPTZ).

**The two idealised rows are not upper bounds here either.** "Truth vertex t0"
and "Ideal track time assignment" zero fewer core PU legs in R1 (34%) than any
real t0 (40–43%). They lower fewer HS legs than TRKPTZ or HGTD (27.0–27.4%
vs 29.4–29.6%), about as many as WAVeS (26.6%).
- Both redraw every timed track at a flat 30 ps: coarser than HGTD's 21–25 ps
  for 2–3-hit tracks.
- The truth row gates at ±94.9 ps.
- So these rows move along the same efficiency/rejection trade-off rather than
  beyond it. CLAUDE.MD's caveat on these rows applies unchanged.

## Definitions

**Regions.** Both come from the max-m_jj opposite-hemisphere pair of jets
above 30 GeV, with paper truth labels: HS means ΔR < 0.3 to a truth HS jet
above 10 GeV; PU means ΔR > 0.6 from every truth HS jet above 4 GeV.
- **R1:** both legs forward, one HS and one PU.
- **R2:** forward PU leg plus central HS leg.
- Forward is |eta| > 2.4 with **no upper edge**; central is |eta| < 2.4.

These are the composition plots' event sets
(`results/vbs_pair_composition.md`).

The region is defined by the pair alone. The composition plots additionally
require ≥1 jet above 30 GeV in 2.38 < |eta| < 4.0 (vbs_region_diag's
preselection); applying that reproduces their counts exactly:

| | R1 | R2 | R1, composition preselection | R2, composition preselection | R1 ∩ narrow | R2 ∩ narrow |
|---|---:|---:|---:|---:|---:|---:|
| VBF | 82,519 | 62,081 | 79,054 | 43,613 | 23,742 | 23,243 |
| Z+jets | 2,324 | 8,161 | 2,165 | 5,064 | 749 | 3,653 |

"Narrow" is the old R_pT region: both forward legs in 2.4–3.8. It still drives
the `_r1/_r2` histograms and the region ROCs. The R2 difference is almost
entirely forward PU legs past |eta| 4.0 in events with no other forward jet
in HGTD acceptance.

**Cases** (existing `rpt_v5` rows; the last two were chosen to be reused
as they are):

| case | row | t0 | track times |
|---|---|---|---|
| TRKPTZ t0 | `trkptz` | TRKPTZ-selected cluster, σ_t0 × 1.39 (VBF) / 1.63 (Z+jets) | real |
| WAVeS t0 | `waves` | WAVeS cluster, in-jet refined, × 1.38 / 1.65 | real |
| HGTD t0 (Athena) | `hgtd` | `RecoVtx_time`, × 1.48 / 1.61; no gate where invalid | real |
| Truth vertex t0 | `truth` | Gaus(t_HS, 10 ps), × 1.0 | every timed track at Gaus(truth time, 30 ps) |
| Ideal track time assignment | `waves_ideal` | smeared-world cluster closest to truth, × 1.38 / 1.65 | every timed track at Gaus(truth time, 30 ps) |

A timed track is kept iff |t_trk − t0| / sqrt(infl² σ_t0² + σ_trk²) < 3.
Untimed tracks are always kept, so timed R_pT ≤ z-only R_pT for every jet.

## Numbers

Legs in 2.4 < |eta| < 3.8; all legs in parentheses. The per-band tables,
including 3.8–4.0 and > 4.0 separately, are in
`figs/rpt_regions/<sample>_rpt_region_summary.md`.

**VBF**

| case | R1 HS: ⟨R_pT⟩ z → t | R1 HS lowered / zeroed | R1 PU: R_pT=0 z → t | R2 PU: R_pT=0 z → t |
|---|---|---|---|---|
| TRKPTZ t0 | 0.374 → 0.353 | 29.6% / 1.4% (25.0% / 1.6%) | 83.7% → 90.6% (92.5% → 95.7%) | 83.7% → 90.2% (92.3% → 95.4%) |
| WAVeS t0 | 0.374 → 0.358 | 26.6% / 0.5% (22.4% / 0.8%) | 83.7% → 90.3% (92.5% → 95.6%) | 83.7% → 87.8% (92.3% → 94.4%) |
| HGTD t0 (Athena) | 0.374 → 0.353 | 29.4% / 1.0% (24.5% / 1.0%) | 83.7% → 90.3% (92.5% → 95.5%) | 83.7% → 88.3% (92.3% → 94.5%) |
| Truth vertex t0 | 0.374 → 0.360 | 27.0% / 0.5% (23.0% / 0.6%) | 83.7% → 89.2% (92.5% → 95.1%) | 83.7% → 89.4% (92.3% → 95.1%) |
| Ideal track time assignment | 0.374 → 0.360 | 27.4% / 0.5% (23.2% / 0.6%) | 83.7% → 89.2% (92.5% → 95.1%) | 83.7% → 89.2% (92.3% → 95.0%) |

**Z+jets**

| case | R1 HS: ⟨R_pT⟩ z → t | R1 HS lowered / zeroed | R1 PU: R_pT=0 z → t | R2 PU: R_pT=0 z → t |
|---|---|---|---|---|
| TRKPTZ t0 | 0.298 → 0.281 | 21.2% / 3.1% (16.6% / 2.7%) | 82.5% → 89.7% (91.2% → 94.7%) | 83.2% → 89.7% (90.9% → 94.4%) |
| WAVeS t0 | 0.298 → 0.287 | 16.7% / 1.2% (13.3% / 1.2%) | 82.5% → 88.9% (91.2% → 94.4%) | 83.2% → 86.4% (90.9% → 92.7%) |
| HGTD t0 (Athena) | 0.298 → 0.283 | 18.0% / 1.8% (14.1% / 1.6%) | 82.5% → 88.5% (91.2% → 94.3%) | 83.2% → 87.4% (90.9% → 93.1%) |
| Truth vertex t0 | 0.298 → 0.287 | 16.8% / 1.4% (13.4% / 1.3%) | 82.5% → 89.0% (91.2% → 94.4%) | 83.2% → 88.7% (90.9% → 93.9%) |
| Ideal track time assignment | 0.298 → 0.287 | 16.4% / 1.4% (13.0% / 1.3%) | 82.5% → 88.5% (91.2% → 94.2%) | 83.2% → 88.3% (90.9% → 93.8%) |

Z+jets core R1 has only 171 PU legs with any z-only R_pT (VBF: 4,896), so its
per-case differences are at the few-leg level. Its HS legs start lower
(⟨R_pT⟩ 0.298 vs 0.374, 8.6% at zero vs 2.4%).

## With JVT + fJVT loose

`rpt_v5_hist --jvt=loose`, clusters 3144162 (VBF, 8 shards) and 3144163
(Z+jets, 14 shards); outputs tagged `jvtLoose_`. The in-code taggers of
`vbs_region_diag --jvt` (`src/jet_tagging.h`) remove failing jets from the
pT-passing list **before** the max-m_jj pair is formed:
- JVT: R_pT proxy (vertex-fit assignment) ≥ 0.0167 for |eta| < 2.5 and
  pT < 60 GeV, the 92% HS-efficiency point.
- fJVT: ≤ 0.5 for 2.5 ≤ |eta| < 4.5 and pT < 120 GeV.

Only the regions change; the inclusive forward slices are tagger-free by
design. A tagger removes at least one pT-passing jet in 69% of VBF region
events (JVT) and 22% (fJVT).

| events | R1 | R2 | R1, composition preselection | R2, composition preselection |
|---|---:|---:|---:|---:|
| VBF, no tagger | 82,519 | 62,081 | 79,054 | 43,613 |
| VBF, JVT+fJVT loose | **60,050** | **52,098** | 57,152 | 34,472 |
| Z+jets, no tagger | 2,324 | 8,161 | 2,165 | 5,064 |
| Z+jets, JVT+fJVT loose | **1,606** | **7,274** | 1,473 | 4,055 |

**The per-leg picture does not change.** fJVT does not use the leg's own
tracks, and JVT acts only below |eta| 2.5. So the taggers change which events
and pairs are in the regions, not what timing does to a leg that is there.
Legs in 2.4 < |eta| < 3.8; all legs in parentheses:

**VBF, JVT + fJVT loose**

| case | R1 HS: ⟨R_pT⟩ z → t | R1 HS lowered / zeroed | R1 PU: R_pT=0 z → t | R2 PU: R_pT=0 z → t |
|---|---|---|---|---|
| TRKPTZ t0 | 0.370 → 0.350 | 29.6% / 1.4% (25.1% / 1.6%) | 82.8% → 90.2% (92.5% → 95.8%) | 82.8% → 89.8% (92.4% → 95.5%) |
| WAVeS t0 | 0.370 → 0.354 | 26.7% / 0.5% (22.5% / 0.7%) | 82.8% → 89.9% (92.5% → 95.7%) | 82.8% → 87.4% (92.4% → 94.5%) |
| HGTD t0 (Athena) | 0.370 → 0.349 | 29.6% / 1.0% (24.8% / 1.0%) | 82.8% → 89.9% (92.5% → 95.6%) | 82.8% → 87.7% (92.4% → 94.6%) |
| Truth vertex t0 | 0.370 → 0.357 | 27.0% / 0.5% (23.1% / 0.6%) | 82.8% → 88.6% (92.5% → 95.1%) | 82.8% → 88.9% (92.4% → 95.1%) |
| Ideal track time assignment | 0.370 → 0.356 | 27.4% / 0.5% (23.4% / 0.6%) | 82.8% → 88.6% (92.5% → 95.1%) | 82.8% → 88.7% (92.4% → 95.0%) |

**Z+jets, JVT + fJVT loose**

| case | R1 HS: ⟨R_pT⟩ z → t | R1 HS lowered / zeroed | R1 PU: R_pT=0 z → t | R2 PU: R_pT=0 z → t |
|---|---|---|---|---|
| TRKPTZ t0 | 0.296 → 0.279 | 21.3% / 3.3% (16.9% / 2.8%) | 82.8% → 89.7% (92.0% → 95.3%) | 82.4% → 89.5% (91.2% → 94.7%) |
| WAVeS t0 | 0.296 → 0.286 | 16.4% / 1.0% (13.1% / 0.9%) | 82.8% → 89.2% (92.0% → 95.2%) | 82.4% → 86.3% (91.2% → 93.2%) |
| HGTD t0 (Athena) | 0.296 → 0.282 | 17.8% / 1.4% (14.2% / 1.4%) | 82.8% → 88.2% (92.0% → 94.8%) | 82.4% → 87.5% (91.2% → 93.6%) |
| Truth vertex t0 | 0.296 → 0.285 | 16.7% / 1.3% (13.3% / 1.2%) | 82.8% → 89.5% (92.0% → 95.1%) | 82.4% → 88.2% (91.2% → 94.1%) |
| Ideal track time assignment | 0.296 → 0.285 | 16.1% / 1.1% (12.9% / 1.1%) | 82.8% → 88.9% (92.0% → 94.8%) | 82.4% → 87.9% (91.2% → 93.9%) |

- **R2 WAVeS weakness survives tagging.** Of the core PU legs with z-only
  R_pT > 0, share zeroed:

  | | TRKPTZ | WAVeS | HGTD t0 |
  |---|---:|---:|---:|
  | VBF | 40% | 27% | 28% |
  | Z+jets | 40% | 22% | 29% |

- **The untimeable share grows slightly.** Legs past |eta| 4.0 are now 60% of
  VBF R1 PU legs and 59% of R2 PU legs (58% / 57% untagged). fJVT reaches
  4.5, but only below 120 GeV, and nothing tags beyond 4.5.

**Checks (all pass):**
- **Region-only changes.** In each tagged histogram file only the 21 region
  histograms differ from the untagged one. The other 63 and every scalar are
  identical; the only new key is `meta_vbs_jvt_wp`.
- **Tree against histograms.** The tagged trees reproduce their 21 region
  histograms bin for bin.
- **Independent implementation.** With the composition preselection, the
  Z+jets tree equals `vbs_region_diag --jvt=loose`'s grid composition exactly
  (1,473 / 4,055), as it does on local VBF (4,183 / 2,579).
- **Merges.** hadd rows equal the shard sums (VBF 112,148; Z+jets 8,880).
- **Job logs.** All 22 jobs: working point logged, zero subset violations,
  empty stderr.

The plots are `figs/rpt_regions/{vbf,zjets}_jvtLoose_rpt_region_dists.pdf`,
PNGs in `figs/rpt_regions/png/`.

## Why the gate gets tracks wrong

`python/rpt_region_pies.py`, from the region tree's attribution columns
(rpt_v5_hist c57a46f; clusters 3144559/60, JVT + fJVT loose). Figures:
`figs/rpt_regions/{vbf,zjets}_jvtLoose_{hs_pt_removed,pu_pt_survives}_pie.(pdf|png)`.

Every call is judged against the truth HS vertex time only, with one
counterfactual: the same gate, with t0 moved to t_HS.

**HS pT removed from HS jets** (R1 forward HS legs; share of the removed pT):

| | TRKPTZ | WAVeS | HGTD t0 |
|---|---:|---:|---:|
| VBF: assigned track time incorrect | **68%** | **79%** | **61%** |
| VBF: incorrect t0, wrong cluster / small offset | 24% / 8% | 17% / 5% | 28% / 11% |
| Z+jets: assigned track time incorrect | **57%** | **64%** | **54%** |
| Z+jets: incorrect t0, wrong cluster / small offset | 34% / 9% | 29% / 7% | 41% / 5% |

The removed pT is 3–5% of the HS legs' HS-track pT. Most of it is lost to the
track's own time, not to the vertex t0: the case for track-time-assignment
studies. Z+jets shifts toward the t0, as expected for a sample whose
hard-scatter tracks are fewer.

**PU pT surviving in PU jets** (R1 + R2 forward PU legs; share of the kept
PU-track pT):

| | TRKPTZ | WAVeS | HGTD t0 |
|---|---:|---:|---:|
| VBF: no HGTD time / no vertex t0 | 36% / – | 32% / – | 33% / 14% |
| VBF: time compatible with t_HS | 32% | 30% | 28% |
| VBF: incorrect t0, wrong cluster / small offset | 29% / 3% | 35% / 2% | 23% / 2% |
| Z+jets: no HGTD time / no vertex t0 | 34% / – | 27% / – | 29% / 23% |
| Z+jets: time compatible with t_HS | 29% | 25% | 21% |
| Z+jets: incorrect t0, wrong cluster / small offset | 35% / 3% | 46% / 2% | 27% / 0.5% |

About a third of the surviving pileup pT carries no HGTD time at all. Another
quarter to a third is in time with the hard scatter, where no t0 could remove
it. WAVeS' larger wrong-cluster share on Z+jets (46%) is its R2 weakness
again.

## Event displays

There are 20 R1 and 20 R2 displays per case and sample: 400 in total, in
`figs/rpt_regions/event_displays/<sample>/<case>/<r1|r2>/`. Each is the full
3-page PDF plus a PNG of page 2, the page used on slides.
- **R1** is ranked by |change in the HS-minus-PU R_pT margin| and tagged
  `helped` / `hurt`.
- **R2** is ranked by the drop in the forward PU leg's R_pT (`suppressed`).
- **The file** records the rank, the tag and the source skim file and entry.
- **The title** gives the case, t0 and the change.

**Every display is checked against the analysis.** It recomputes both legs'
z-only and timed R_pT from that case's exact t0, σ_t0, inflation and gate,
and prints `RPT_CHECK`. The two idealised cases also replay the smeared times
through `src/idealised_timing.h` and print `REPLAY_CHECK`: the smeared vertex
time, or the truth-closest smeared cluster, must equal the tree's t0.
**400/400 pass both.** Z+jets pairs are passed explicitly (`--legs`), because
the display's own pair search does not apply the lepton overlap removal.

**Read the R1 lists as mostly failures.**

| R1 displays tagged `hurt` (of 20) | TRKPTZ | WAVeS | HGTD | Truth t0 | Ideal |
|---|---:|---:|---:|---:|---:|
| VBF | 19 | 16 | 19 | 12 | 15 |
| Z+jets | 18 | 18 | 18 | 14 | 14 |

A ranking by |Δmargin| favours large changes. The largest come from a wrong
t0 emptying the HS leg (an HS leg can lose 0.8 of R_pT). Suppressing the PU
leg gains at most that leg's own small R_pT. A balanced helped/hurt split
needs a separate ranking; the tree has everything required.

The candidate lists are
`condor/<sample>/<sample>_region_display_candidates.csv`. They carry every
render argument, so a single display can be regenerated from its row.

**JVT + fJVT loose displays.** 400 more are in
`figs/rpt_regions/event_displays/<sample>_jvtLoose/<case>/<r1|r2>/`, selected
from the tagged trees in the same way. **400/400 pass RPT_CHECK and
REPLAY_CHECK.**
- Jets the tagger removed before pairing are drawn faded with a dotted edge
  and labelled `[fails JVT]` / `[fails fJVT]`. This explains, on the display
  itself, why a visibly larger jet is not a leg.
- 134 of 200 VBF displays and 120 of 200 Z+jets displays show at least one
  such jet (33 each show an fJVT removal).
- The second title line names the working point.
- The R1 lists stay mostly failures:

  | R1 displays tagged `hurt` (of 20) | TRKPTZ | WAVeS | HGTD | Truth t0 | Ideal |
  |---|---:|---:|---:|---:|---:|
  | VBF | 19 | 19 | 19 | 13 | 18 |
  | Z+jets | 18 | 17 | 18 | 15 | 14 |

## Checks

- **Refactor.** Moving the idealised timing into `src/idealised_timing.h` left
  `rpt_v5_hist` bit-identical on local VBF: 84/84 histograms, 18/18 scalars.
- **Grid regression.**
  - VBF is identical to the unmerged 2026-09-18 run (same code path, same
    skim) in every histogram and scalar.
  - Z+jets is identical to the 2026-09-01 run in all 84 histograms. The only
    differences are the 6 scalars counted before the skim's vertex+lepton cuts
    (that run read the original ntuples).
- **Tree against histograms.** The tree's narrow (`core`) rows reproduce all
  21 region histograms bin for bin on both samples.
  `rpt_region_dists.py --hist-file` reruns this check every time.
- **Tree against an independent implementation.** With the composition
  preselection, the tree's counts equal `vbs_region_diag`'s exactly on both
  samples (table above).
- **Merges.** hadd'd tree rows equal the sum over shards: VBF 144,600;
  Z+jets 10,485.
- **Job logs.** Zero narrow-not-in-wide events in all 22 jobs; empty stderr.

## Reproducing

```bash
# AF: histograms + tree (both merged; see CLAUDE.MD "Sharding a sample")
condor_submit -a sample=vbf -a nshards=8 condor/rpt_v5_hist.sub     # from condor/
./build/hist_merge condor/vbf/vbf_rpt_v5_hist.root condor/vbf/vbf_rpt_v5_hist.shard*of8.root
hadd -f condor/vbf/vbf_rpt_v5_regions.root condor/vbf/vbf_rpt_v5_regions.shard*of8.root

# laptop: distributions (+ the core-subset check)
export PYTHONNOUSERSITE=1 PYTHONPATH=/opt/homebrew/Cellar/root/6.40.04/lib/root
~/.venv-hgtd/bin/python python/rpt_region_dists.py condor/vbf/vbf_rpt_v5_regions.root \
    --sample vbf --hist-file condor/vbf/vbf_rpt_v5_hist.root

# JVT + fJVT loose: the same, with the tag (outputs vbf_jvtLoose_*; condor from condor/):
#   condor_submit -a sample=vbf -a nshards=8 -a 'selargs=--jvt=loose' -a seltag=jvtLoose_ rpt_v5_hist.sub
#   ... then the same merge/plot/select/pick/render with vbf_jvtLoose in every name and --sample vbf_jvtLoose

# displays: select (laptop) -> pick (AF, needs the skims) -> render (laptop)
~/.venv-hgtd/bin/python python/rpt_region_displays.py select condor/vbf/vbf_rpt_v5_regions.root --sample vbf
#   AF, after lsetup root: use the LCG python by absolute path and check `import ROOT` first --
#   `lsetup` intermittently leaves ROOTSYS unset, and bash can keep /bin/python hashed.
#   $PY python/rpt_region_displays.py pick condor/vbf/vbf_region_display_candidates.csv   (~6 min, cold ceph)
~/.venv-hgtd/bin/python python/rpt_region_displays.py render condor/vbf/vbf_region_display_candidates.csv \
    --picked condor/vbf/vbf_region_display_events.root \
    --pick-map condor/vbf/vbf_region_display_events_map.csv --jobs 6          # ~2.5 min / 200
```
