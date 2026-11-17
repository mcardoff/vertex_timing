---
name: vbs-timing-studies
description: VBS-region timing studies — rpt_v5's wide-window region tree and region R_pT plots/displays, the vbs_time_veto event veto / jet removal (incl. the VBF H->inv SR, minimal/relaxed sets, HS+PU study), and the VBS pair composition with in-code JVT/fJVT. Load before running, plotting or quoting any of them.
---

## rpt_v5_hist's wide-window region tree, and the region R_pT distributions/displays

**Wide-window region tree** (2026-09-25): beside the histogram file,
`rpt_v5_hist` writes `<prefix>rpt_v5_regions.root`, tree `regions`, one row
per event whose max-m_jj pair is R1/R2 with forward = |eta| > 2.4 and **no
upper edge**. These are the same event sets as the composition plots, legs
past 3.8 included. Per row it stores:
- every scenario's per-leg R_pT (`rpt_hs_<s>`, `rpt_pu_<s>`; R2's
  `rpt_hs_*` = −1, since its HS leg is central);
- each scenario's exact gate inputs (`t0_<s>`, `sig0_<s>` before inflation,
  `infl_<s>`, `ok_<s>`);
- leg kinematics, and the z-associated / timed track counts in each leg's
  R_pT cone;
- `core`, flagging the rows that are also in the narrow `_r1/_r2` region
  (a strict subset);
- `n_jets_fwd_acc`, the jets >30 GeV in 2.38 < |eta| < 4.0.
`vbs_region_diag`'s preselection requires ≥1 such jet, and requiring it
here reproduces the composition counts exactly. Without it R2 is ~40%
larger (local VBF: 4,744 vs 3,350), almost all of that from PU legs past
|eta| 4.0 with R_pT = 0 in every row.
**`--jvt=none|loose|tight`** (2026-09-25) applies the in-code JVT/fJVT
taggers (`src/jet_tagging.h`, identical to `vbs_region_diag --jvt`) to the
pT-passing jets BEFORE the VBS pair is formed. It gates the regions only
(narrow histograms + wide tree), like the `--vbs-*` knobs; the inclusive
forward slices stay tagger-free. A tagged run:
- turns on EXTENDED_BRANCHES (+ `recordAvailableBranches`), and refuses a
  sample without `Track_recoVtx_idx` rather than removing every jet;
- tags its outputs `jvtLoose_`/`jvtTight_` (SELECTION_TAG; condor
  `selargs=--jvt=loose seltag=jvtLoose_`);
- writes `meta_vbs_jvt_wp`;
- fills `jvt_wp`, each leg's JVT R_pT / fJVT, and the removed jets
  (`rm_jvt`, `rm_fjvt`, drawn by the displays) in the tree.
Validated on local VBF: `--jvt=none` is bit-identical to before, and
`--jvt=loose` reproduces `vbs_region_diag --jvt=loose`'s composition counts
exactly (R1 4,183, R2 2,579 with its preselection).
Merge shards with `hadd` (TTrees only); `hist_merge` refuses the key type,
which is why it is a separate file.
**The tree has no event weight, pT(ll) or MET column.** For the grid
trees, three per-row friend files (same row order) sit beside them in
`condor/<sample>/`:
- `<stem>_region_weights.root`: the generator weight;
- `<stem>_zpt.root`: Z+jets pT(ll), `vbs_region_diag`'s `z_pt` definition;
- `<stem>_truth_met.root`: VBF truth |Σ p_T^ν|, recovered from the raw
  ntuples by a (weight, RecoVtx_z[0], TruthVtx_z[0]) key, since the skims
  drop the truth record.
The macros are in `util/scratch/region_weights/`. The method and the
3000 fb⁻¹ yields are in `results/rpt_region_distributions.md`. Writing these
columns in `rpt_v5_hist` itself would retire the lookups. Standing checks, all passing on local
VBF: the core rows reproduce all 21 region histograms bin for bin
(`python/rpt_region_dists.py --hist-file` runs this every time), and zero
narrow events fall outside the wide region.
The narrow `_r1/_r2` histograms are kept for the region ROCs. A wide ROC is
not meaningful: HS legs past |eta| 4.0 have no tracks, which puts the
HS-efficiency ceiling below the 0.80 window.

**Region R_pT distributions and displays** (from that tree):
- `python/rpt_region_dists.py` draws, per timing case (TRKPTZ, WAVeS, HGTD,
  and the `truth` / `waves_ideal` rows presented as "Truth vertex t0" /
  "Ideal track time assignment"):
  - R1-HS, R1-PU and R2-PU pages, z-only (filled) vs timed, with a ratio pad;
  - an |eta|-band page (2.4–3.8 vs >3.8);
  - a markdown summary table.
- `python/rpt_region_displays.py select|pick|render` produces 20 R1 + 20 R2
  event displays per case. `pick` copies the events out of the skims (AF
  only) so `render` can run on the laptop. Every display must pass
  RPT_CHECK, and the idealised rows REPLAY_CHECK (see Event Display).
- Results, grid VBF + Z+jets, untagged and JVT+fJVT loose:
  `results/rpt_region_distributions.md`.
  - Legs past |eta| 4.0 carry no tracks and are 57–60% of VBF PU legs.
  - In R2, WAVeS zeroes far fewer fakes than TRKPTZ (26–27% vs 40% on VBF).
  - The loose tagger removes 27% of VBF R1 events and 16% of R2 but leaves
    the per-leg picture unchanged.
  - A |Δmargin| ranking of R1 displays is mostly "timing hurt" events.


## Event-level timing veto (2026-09-28)
- **vbs_time_veto.cxx** (`util/`) + **python/vbs_time_veto_plot.py**: Friday
  2026-09-25's question. R1 needs no t0: veto the event when the two tagging
  jets' times are incompatible, |t_A - t_B| / sqrt(sA^2 + sB^2) >= 3. That is
  compared with testing EACH timed leg against the event t0, which also reaches
  R2. The metric is weighted signal (VBF, truth MET > 200) and background
  (Z+jets, pT(ll) > 200) efficiency. Results: `results/vbs_time_veto.md`.
  - **Answer (grid, JVT + fJVT loose): neither veto improves S/sqrt(B) after
    the recoil cut.**
    - Jet vs jet: eps_S 95.1%, eps_B 96.0%, S/sqrt(B) x 0.97. Only 6.5% of the
      background has two timed legs, and R2 (a central leg) is 35% of it.
    - Jet vs t0: eps_S ~86%, eps_B ~70%, x 1.02-1.04 (HGTD 0.98). That is
      within the +-2.4% that 204 Z+jets MC events allow, but it beats jet vs
      jet at matched eps_S.
    - A PERFECT t0 does no better. The signal loss is the HS jet time's tail:
      7% of HS legs fail 3 sigma against the true t_HS, half of them
      single-track.
  - **The executable writes ingredients, never decisions.** For every event in
    the region plots' selection it stores both legs' jet times, all five t0s
    with sigma and inflation, the weight, and both recoil proxies. Threshold,
    jet-time calibration and cuts are applied offline, with weights.
  - **The jet time is the crux.** Averaging ALL timed ghost tracks, the literal
    prescription, gives a pileup average: only 21% of an HS jet's are
    hard-scatter, chi2/ndf ~ 65. It passes 25% of HS+HS pairs, against 21% of
    R1. The default `core_` time keeps timed ghost tracks passing
    passTrackKinematics within dR < 0.2, clusters them in time (iterative,
    2 sigma), and takes the highest-sum-pT cluster.
  - **The jet-time calibration is ONE procedure for S and B.** sigma_jet is
    scaled by f(n tracks), measured on VBF HS legs: 2.34 at n = 1 falling to
    1.31 at n >= 9. The same factors apply to both samples. A per-sample
    calibration (global 1.51 VBF / 1.83 Z+jets) would give the background a
    looser veto than the signal; `--jet-infl per-file` shows it moves nothing.
  - **VBF must be read from the RAW ntuples.** Its recoil cut is truth
    |sum pT(nu)| > 200 GeV, and the skims drop the truth record.
    `condor/vbs_time_veto.sub` passes it via `inputargs=--ntuple-dir=...`. On a
    skim, met_truth = -1 and the plot script refuses the cut. It also checks
    completeness: VBF 1,464,400 events and the raw sum of weights (matched
    exactly); Z+jets 530,962 skim events.
  - **Validated.**
    - Every reco t0 equals rpt_v5's JVT-loose region trees event by event
      (91,624 VBF, 5,528 Z+jets). The VBF run read the raw files and the trees
      the skims.
    - Truth MET and pT(ll) equal last week's friend files.
    - Z+jets' 26,883 selected / 204 after pT(ll) > 200 equal the region
      stacks' `vbs_region_diag` counts.
    - Local VBF's event set equals `vbs_region_diag --jvt=loose`'s exactly.
  - **Job cost:** VBF raw, 20 shards, ~11 min each; Z+jets skim, 8 shards,
    ~3 min each.
  - **Second method, jet removal + re-pairing (`rp_*`).** The event veto
    rejects signal R1 events, and R1 is as common in VBF as in Z+jets. So the
    tree now also carries every candidate jet (`jet_*`), and the plot script
    removes the jets failing the t0 test and re-forms the pair. The forward-jet
    preselection is not re-applied after removal.
  - **Result (grid, recoil cuts at 100 GeV): jet removal does not beat the
    event veto.** It keeps 90-92% of the signal (event veto 84-86%) but
    78-85% of Z+jets (event veto 64-76%). S/sqrt(B) is x1.03-1.04 against
    x1.04-1.06. Re-pairing recovers background as readily as signal.
  - The event veto with a TRKPTZ / WAVeS / TZP t0 stays the best use of the
    timing, about +5%.
  - **In the VBF H->inv SIGNAL REGION (2026-09-29, `--sr` +
    `python/vbs_time_veto_sr.py`, `results/vbs_time_veto_sr.md`).** The
    selection is ANA-EXOT-2020-11-INT1's SR, emulated offline; the internal
    note is `~/Downloads/ATL-COM-PHYS-2020-310.pdf`. It includes JVT, fJVT on
    j1/j2 with the E_T^miss-dependent working point, 2-4 jets, the j3/j4
    centrality and m_rel cuts, E_T^jet,no-jvt, the dphi/deta/m_jj cuts and the
    16-bin union.
    - The SR flips the jet-removal conclusion on the signal side. The event
      veto costs 7-9% of the signal: SR pairs are 99.4% HS+HS. Removing
      failing jets first GAINS 2-5%: removed out-of-time pileup jets no longer
      fail the extra-jet cuts.
    - Pileup jets cost the SR 42% of its signal (truth ceiling x1.71). Only
      x1.14 of that lies in HGTD acceptance.
    - **The background cannot be measured: 4 Z+jets MC events survive the
      SR.** Inclusive Powheg DY has too little VBF topology; the whole sample
      gives 223 candidates.
    - `--sr` exists because 14% of SR candidates fail the region-plot
      preselection. It stores every jet > 20 GeV with the JVT R_pT, fJVT and
      windows.
    - EMTopo, not the note's PFlow: this production's PFlow is miscalibrated
      forward (HS response 1.14 at 2.4-4.0, 1.31 beyond).
    - `python/vbs_time_veto_common.py` is shared by both plot scripts, which
      is validated byte- and pixel-identical.
    - Job cost: VBF raw ~22 min per shard (it reads the truth record for
      every event), Z+jets ~5 min.
    - **HS+PU study (`python/vbs_hspu_study.py`).** It asks how often a
      pileup jet is a tagging jet, and whether the t0 test catches it.
      - Multijet = dijet raw with `--sr --sr-met=jets`: `met_jets` is
        |sum pT| of the JVT-passing jets > 20 GeV, since truth MET is ~0.
      - The JZ slices carry pT-biased weights (5e-9 to 1), so quote
        effective counts; they are not normalised to each other.
      - The multijet MC reaches the tagging-pair stages only: nothing
        survives |dphi_jj| < 2 plus the jet veto.
      - At the pair stage (MET > 100) the HS+PU share is Z+jets 44%,
        JZ2 53%, JZ3 10%, JZ4 0%, and VBF 5%.
      - Of the PU tagging jets, 1/3-1/2 lie beyond |eta| 4. The t0 test
        catches ~2/3 of those in HGTD.
      - **The event veto pays when the pileup-tagged background share f
        exceeds f* ~ 30%.** That is with eps_S, c and h all in the minimal
        set: WAVeS 31%, truth 30%, TRKPTZ 35%, HGTD 47%. An earlier 20%
        mixed selections (eps_S from the full SR, c and h at the pair
        stage).
      - Our Z+jets at 34-44% overstates f (Powheg DY). Run 2 measured
        6-13%.
    - **WITHDRAWN 2026-10-01: the "SR background, estimated" numbers.** They
      took S from the full SR, and B from Z+jets at +dphi times a transfer
      factor, so S and B sat in different selections. Its "jet removal
      +10 to +12%" came from that factor. **S/sqrt(B) only ever compares S
      and B measured in the same selection.**
    - **Minimal cut set** (`python/vbs_sr_minimal.py`, 2026-10-01). Only
      B(timing)/B(none) enters S/sqrt(B), so a cut uncorrelated with the
      timing decision can be dropped. Each cut is tested pass vs fail on f
      and on the veto's kept fraction.
      - Keep: the pair (j1/j2 pT, deta, m_jj), fJVT, and E_T^miss > 100.
        Also keep |dphi_jj| < 2: a back-to-back multijet pair passes it only
        with a PU leg (+5 sigma).
      - Drop: E_T^miss > 160/200 with the bins, E_T^jet,no-jvt, and the jet
        veto (<= 1.6 sigma everywhere). That leaves Z+jets 52 MC events
        (SR: 4).
      - The pair cuts and the floor cannot be loosened. f depends on m_jj,
        deta and j2 pT within the SR range. Below pT(ll) = 100, f is 81-85%
        (`--sr-recoil-min=0`, tagged `sr_rec0p0`).
      - S/sqrt(B), with S and B both in the minimal set
        (`vbs_time_veto_sr.py --region minimal`, figs/time_veto/srmin_*):
        - event veto: +3.6% WAVeS, +2.9% TZP, +2.0% truth, +0.6% TRKPTZ,
          -3.8% HGTD (+-4%);
        - jet removal: 0.8 points more. Without the jet veto it is just the
          veto plus re-pairing.
      - At Run 2's f = 6-13% the event veto costs 3-4%.
    - **Relaxed set** (2026-10-01, REGIONS["relaxed"]): the minimal set with
      j1 > 40, j2 > 30 GeV and |m_ll - 91.2| < 40 GeV, on the `vbs_time_veto
      --sr --sr-j1-pt=40 --sr-j2-pt=30` files (`sr_j1pt40p0_j2pt30p0`).
      - Z+jets: 110 MC events.
      - S/sqrt(B), S and B both in it: +2.0% WAVeS, +1.9% TZP, +2.2% truth,
        +0.2% TRKPTZ, -2.8% HGTD (+-2.8%). Jet removal is the same.
      - The added j2 < 50 GeV Z+jets pairs are 82% pileup-tagged, so f rises
        31% -> 55%, and f* = 40%.
      - Validated: under the 80/50 file's calibration (`--calib-file`), those
        files reproduce the SR and minimal-set results exactly.
      - eps_B equals f x the signal's per-type response within 1 sigma. So
        the background MC is needed for f (and where its PU legs sit), not
        for the response.


## VBS pair composition, and what JVT/fJVT do to it (2026-09-02 / 2026-09-18)

`util/vbs_region_diag` classifies every selected event's chosen tagging pair
by the truth identity and detector region of each leg (21 cells; R1 and R2
are two of them). `results/vbs_pair_composition.md`: **Z+jets has no
hard-scatter leg to pair with** — at m_jj >= 200, 79.7% of its pairs carry no
paper-HS jet and the two pileup-only cells hold 65.6%, while VBF puts an HS
jet on both legs 55% of the time on the same code. R2's scarcity on Z+jets is
the central HS leg mostly not being there, not the m_jj picker outranking it.

The SuperNtuples carry **no Jvt/fJvt branch**, so `--jvt=loose|tight`
computes both in-code (`src/jet_tagging.h`, since 2026-09-25 shared with
`rpt_v5_hist --jvt`) from `Track_recoVtx_idx` + ghost tracks (R_pT as the
JVT proxy, calibrated to the published 92% / 85% HS efficiencies on local
VBF; fJVT per arXiv:1705.02211 with the published 0.5 / 0.4 cuts) and
applies them BEFORE pairing. `results/vbs_pair_composition_jvt.md`:

- **Z+jets: the taggers do not make it an R2 sample.** They remove central
  pileup (JVT PU eff 15% / 1.8%) far better than forward pileup (fJVT PU eff
  67% / 47% at μ = 200), so F-PU + C-PU drains 37.8 → 5.7% while F-PU + F-PU
  GROWS 27.8 → 42.0% of pairs. R2 rises 7.8 → 20.3% by losing its
  competition, not by creation (78% of post-loose R2 was R2 before); even
  tight leaves 60% of pairs with no HS leg. 59% / 77% of selected events are
  dropped. The no-tagger run reproduces the earlier 64,950 pairs exactly.
- **VBF, the control:** removal is 98.5% pileup-pure on the JVT side, HS+HS
  share of the pair rises 56 → 71 → 78%, R1/R2 shrink only modestly and 70% /
  67% of R1 / R2 events are still R1 / R2 after the loose point.
- **The forward-pileup dominance is a selection artefact plus a sample
  limitation.** The diagnostic has no recoil requirement; 64% of its selected
  Z+jets events have no HS jet above 30 GeV, and two forward pileup jets pass
  the VBS cuts alone in 3.7% of all μ = 200 events. With pT(ll) > 200 (the
  analysis' MET cut in the Z→ll CR) the survivors are R2 36% + HS+HS 34% and
  F-PU+F-PU 7% — but only 0.9% of events survive, since inclusive Powheg DY
  has almost no Z above 200 GeV, and the sample has ≤1 ME parton and no EW
  Zjj. `z_pt`/`z_mass` are stored per event; `vbs_region_stack.py --zpt-min`.
- dijet not run. `--jvt=none` is bit-identical to the pre-JVT diagnostic.

