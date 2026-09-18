// -----------------------------------------------------------------------------
// vbs_region_diag -- why does the Z+jets VBS pair almost never land in R2?
//
// The four-band composition (results/vbs_region_mjj_fourband.md) puts 92.5% of
// Z+jets pairs above m_jj 500 in "Other" and only 3.46% in R2, against the
// physical expectation that Z+jets should be DOMINATED by R2 -- a forward
// pileup fake paired with the central hard-scatter jet the Z recoils against.
// "Other" was never decomposed, so the disagreement had three live
// explanations that no existing diagnostic could separate:
//
//   (a) the topology is absent -- the event has no central hard-scatter jet
//       above MIN_JET_PT to be the second leg, so no pair could be R2;
//   (b) the topology is present but the PAIR PICKER misses it --
//       calcBestVbsPair takes the largest-m_jj opposite-hemisphere pair, and
//       m_jj grows like cosh(dEta), so a wider pair outranks a genuine
//       forward-pileup + central-hard-scatter pair that is actually there;
//   (c) the winning pair is not an ANALYSIS pair at all -- calcBestVbsPair
//       ranks over every jet above MIN_JET_PT at any |eta|, including jets
//       beyond MAX_ABS_ETA_JET; and "forward" is currently HGTD acceptance
//       (MIN_ABS_ETA_JET < |eta| < MAX_ABS_ETA_JET), so a jet past 4.0 is
//       neither forward nor central and can only ever land the pair in Other.
//
// So every quantity here is computed THREE times over the same event, from the
// same jets, differing only in what counts as a taggable jet and where
// "forward" ends:
//
//   all_   every pT-passing jet; forward = MIN_ABS_ETA_JET..MAX_ABS_ETA_JET
//          -- exactly what calcBestVbsPair and classifyVbsRegion do today.
//   acc_   only pT-passing jets with |eta| < MAX_ABS_ETA_JET; same forward
//          window -- analysis-quality tagging jets only, so a jet the analysis
//          would never tag on cannot win the m_jj ranking.
//   wide_  every pT-passing jet; forward = anything above MIN_ABS_ETA_JET with
//          NO upper edge -- "forward" as a topology statement rather than a
//          detector-acceptance one. A tagging jet at |eta| 4.3 is still a
//          forward tagging jet; HGTD simply cannot time it.
//
// (c) is then the difference between the blocks, and (a) vs (b) is a two-way
// table within any one of them.
//
// Each block also carries, per pair SHAPE, the best-m_jj pair of that shape
// that EXISTS in the event -- with its legs -- whether or not the picker chose
// it. Comparing that pair's kinematics against the chosen one is what shows
// which jets the m_jj ranking is actually reaching for.
//
// No clustering and no time gates: the question is purely about jet content,
// and dropping them keeps a full Z+jets pass cheap.
//
// JVT / fJVT (2026-09-18, on Ariel's feedback). An analysis applies the
// standard pileup-jet taggers BEFORE it forms a tagging pair, so the
// composition above describes a jet population no analysis actually sees. The
// SuperNtuples carry no Jvt/fJvt decoration (241 branches, none match, on
// either jet collection), so both are computed here from the vertex fit's own
// track assignment (Track_recoVtx_idx) and each jet's ghost-associated tracks:
//
//   R_pT       = sum pT(ghost tracks fitted to vertex 0) / pT_jet   [1510.03823]
//   corrJVF    = pT_PV / (pT_PV + pT_PU / (k n_PU^trk)),  k = 0.01  [stored only]
//   fJVT       = max_{i>0} (p_T^miss,i . j_T) / pT_jet              [1705.02211]
//   p_T^miss,i = -1/2 ( sum_{trk fitted to i, |eta|<2.5} p_T
//                     + sum_{central jets whose dominant ghost vertex is i} p_T )
//
// JVT proper is a k-NN likelihood over (corrJVF, R_pT) that cannot be rebuilt
// from the ntuple, so the "JVT" cut here is R_pT alone, with thresholds set so
// the paper-HS efficiency on local-VBF central jets inside the JVT window
// reproduces the published working-point efficiencies (JVT_WPS below). fJVT
// thresholds are the published ones. Windows: JVT |eta| < 2.5 and pT < 60 GeV;
// fJVT 2.5 <= |eta| < 4.5 and pT < 120 GeV. A jet outside both windows is never
// removed; a jet with no ghost tracks has R_pT = 0 and fails JVT in its window.
//
// --jvt=none|loose|tight picks the working point. The default, none, is
// bit-identical to the pre-JVT diagnostic apart from the added columns. The
// `jets` tree holds every pT-passing jet's discriminants BEFORE the cut, so a
// threshold can be re-derived (python/vbs_jvt_calibrate.py) or scanned offline
// without a rerun.
// -----------------------------------------------------------------------------
#include <TChain.h>
#include <TFile.h>
#include <TLorentzVector.h>
#include <TTree.h>
#include <TTreeReader.h>
#include <TVector2.h>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <string>
#include <vector>

#include "clustering_constants.h"
#include "sample_config.h"
#include "clustering_includes.h"
#include "clustering_structs.h"
#include "clustering_functions.h"
#include "event_processing.h"

using namespace MyUtl;

namespace {

// Zone boundaries: exactly classifyVbsRegion's isFwd/isCentral as
// vbs_region_mjj calls it (central = |eta| < MIN_ABS_ETA_JET), so a category
// here cannot mean something different from the four-band plot's.
// fwdMax is the upper edge of "forward"; pass a huge value for the wide_
// convention, where forward has no upper edge and zone 2 is unreachable.
int zoneOf(double absEta, double fwdMax) {
  if (absEta < MIN_ABS_ETA_JET) return 0;   // central
  if (absEta < fwdMax)          return 1;   // forward
  return 2;                                 // beyond "forward"
}
const double FWD_MAX_UNBOUNDED = 1e9;

// ── JVT / fJVT proxies (see the file header) ────────────────────────────────
const double JVT_ETA_MAX         = 2.5;    // JVT window: |eta| < 2.5 ...
const double JVT_PT_MAX          = 60.0;   // ... and pT < 60 GeV
const double FJVT_ETA_MIN        = 2.5;    // fJVT window: 2.5 <= |eta| < 4.5 ...
const double FJVT_ETA_MAX        = 4.5;
const double FJVT_PT_MAX         = 120.0;  // ... and pT < 120 GeV
const double FJVT_CEN_JET_PT_MIN = 20.0;   // central jets entering p_T^miss,i
const double FJVT_TRK_ETA_MAX    = 2.5;    // tracks entering p_T^miss,i
const double CORRJVF_K           = 0.01;

struct JvtWP {
  const char* name;     // --jvt= value
  const char* tag;      // output-file tag
  double      rptMin;   // JVT proxy: keep a JVT-window jet iff R_pT >= this
  double      fjvtMax;  // keep an fJVT-window jet iff fJVT <= this
};
// R_pT thresholds: calibrated on local-VBF paper-HS jets in the JVT window
// (|eta| < 2.5, 30 < pT < 60 GeV) to the published JVT working-point HS
// efficiencies -- loose <-> the default Run-2 "Medium" point (92%), tight <->
// "Tight" (85%) -- by python/vbs_jvt_calibrate.py on the `jets` tree of a
// --jvt=none run. fJVT: published Loose 0.5 / Tight 0.4.
const JvtWP JVT_WPS[] = {
  {"none",  "",          -1.0, 1e9},
  {"loose", "_jvtLoose", 0.0167, 0.5},   // 92.0% HS eff on local VBF, PU eff 13.2%
  {"tight", "_jvtTight", 0.0767, 0.4},   // 85.0% HS eff on local VBF, PU eff  1.9%
};

// One jet's discriminants, computed for EVERY jet in the array (not just the
// pT-passing ones) because the fJVT of a forward jet depends on the central
// jets' vertex assignment.
struct JetDisc {
  float rpt = 0.f, corrjvf = -1.f, fjvt = 0.f;
  int   nGhost = 0, nGhostPV = 0;
  bool  jvtWin = false, fjvtWin = false;
};

std::vector<JetDisc> computeJetDiscriminants(const BranchPointerWrapper& b) {
  const int nJ = (int)b.topoJetPt.GetSize();
  const int nV = (int)b.recoVtxZ.GetSize();
  std::vector<JetDisc> D(nJ);

  // Per-vertex transverse momentum of the tracks the fit assigned to it,
  // central tracks only (fJVT's p_T^miss is a central quantity). Index -1 is
  // "fitted to no vertex" (47% of tracks) and contributes nowhere.
  std::vector<double> vpx(nV, 0.0), vpy(nV, 0.0);
  int nPUtrk = 0;
  for (int t = 0; t < (int)b.trackPt.GetSize(); ++t) {
    const int v = b.recoVtxOf(t);
    if (v < 0 || v >= nV) continue;
    if (v > 0) ++nPUtrk;
    if (std::abs((double)b.trackEta[t]) >= FJVT_TRK_ETA_MAX) continue;
    vpx[v] += b.trackPt[t] * std::cos(b.trackPhi[t]);
    vpy[v] += b.trackPt[t] * std::sin(b.trackPhi[t]);
  }

  // Per-jet ghost-track sums split by vertex; a central jet is then assigned
  // to whichever vertex dominates its ghost pT, and if that is a PILEUP vertex
  // the jet's pT enters that vertex's p_T^miss. PV-dominated (i.e. hard
  // scatter) jets enter no PU vertex's balance. Overlap-removed jets are
  // leptons and are skipped for the assignment only.
  std::vector<double> perV(nV);
  for (int j = 0; j < nJ; ++j) {
    std::fill(perV.begin(), perV.end(), 0.0);
    double sumPV = 0.0, sumPU = 0.0;
    JetDisc& d = D[j];
    d.nGhost = (int)b.topoJetGhostTrackIdx[j].size();
    for (int idx : b.topoJetGhostTrackIdx[j]) {
      const int v = b.recoVtxOf(idx);
      if (v < 0 || v >= nV) continue;
      perV[v] += b.trackPt[idx];
      if (v == 0) { sumPV += b.trackPt[idx]; ++d.nGhostPV; }
      else          sumPU += b.trackPt[idx];
    }
    const double pt = b.topoJetPt[j], aeta = std::abs((double)b.topoJetEta[j]);
    d.rpt     = (float)(sumPV / pt);
    d.corrjvf = (sumPV + sumPU > 0.0)
              ? (float)(sumPV / (sumPV + sumPU / (CORRJVF_K * std::max(nPUtrk, 1))))
              : -1.f;
    d.jvtWin  = aeta < JVT_ETA_MAX && pt < JVT_PT_MAX;
    d.fjvtWin = aeta >= FJVT_ETA_MIN && aeta < FJVT_ETA_MAX && pt < FJVT_PT_MAX;
    if (aeta < FJVT_TRK_ETA_MAX && pt > FJVT_CEN_JET_PT_MIN && !b.isJetRemoved(j)) {
      int vBest = -1; double best = 0.0;
      for (int v = 0; v < nV; ++v) if (perV[v] > best) { best = perV[v]; vBest = v; }
      if (vBest > 0) {
        vpx[vBest] += pt * std::cos(b.topoJetPhi[j]);
        vpy[vBest] += pt * std::sin(b.topoJetPhi[j]);
      }
    }
  }

  // fJVT: the pileup vertex whose missing transverse momentum points most
  // along the jet, normalised to the jet pT. Computed for every jet; the
  // window is applied by the caller.
  for (int j = 0; j < nJ; ++j) {
    const double pt = b.topoJetPt[j];
    const double ux = std::cos(b.topoJetPhi[j]), uy = std::sin(b.topoJetPhi[j]);
    double best = -std::numeric_limits<double>::infinity();
    for (int v = 1; v < nV; ++v) {
      const double proj = -0.5 * (vpx[v] * ux + vpy[v] * uy) / pt;
      if (proj > best) best = proj;
    }
    D[j].fjvt = (nV > 1) ? (float)best : 0.f;
  }
  return D;
}

// Everything the diagnostic asks of one jet collection. Filled twice per event
// -- see the file header -- so the two must stay structurally identical.
struct PairBlock {
  float n_jets;
  float pair_mjj, pair_deta;
  float legA_zone, legA_hs, legA_pu, legA_pt, legA_abseta;
  float legB_zone, legB_hs, legB_pu, legB_pt, legB_abseta;
  // Best pair of each shape present in the event, chosen or not.
  float alt_r1_mjj, alt_r2_mjj, alt_pupu_mjj, alt_bothhs_mjj;
  // ...and its legs, named by the ROLE each plays in the shape rather than by
  // index: for R2 "which leg is the fake" is the whole question, and an (i,j)
  // ordering would not say.
  float alt_r2_fwd_abseta, alt_r2_fwd_pt, alt_r2_cen_abseta, alt_r2_cen_pt;
  float alt_r1_hs_abseta,  alt_r1_hs_pt,  alt_r1_pu_abseta,  alt_r1_pu_pt;
  // The chosen legs' JVT/fJVT discriminants, filled under every working point
  // (including none) so the offline side can see what a cut WOULD do.
  float legA_rpt, legA_fjvt, legB_rpt, legB_fjvt;
};

struct Row {
  float file_idx, entry;
  // Event-wide jet content under the paper labels (independent of any pairing).
  float n_all_jets, n_fwd_jets, n_beyond_jets;
  float n_fwd_hs, n_fwd_pu, n_cen_hs, n_cen_pu, n_any_hs, n_neither;
  // Same labels for jets past MAX_ABS_ETA_JET, counted separately so the wide_
  // convention's forward counts are just fwd + beyond rather than a rerun.
  float n_beyond_hs, n_beyond_pu;
  float n_fwd_hs_trk;          // forward truth-HS tracks (the t0 question)
  float lead_pt, lead_abseta;
  // Truth-side jet content, independent of any reco match: answers whether
  // "no paper-HS reco jet anywhere" (n_any_hs == 0) means no truth jet
  // existed, or one existed and was too soft / unmatched. NOT gated on
  // passPtIdx -- this is every TruthHSJet in the event, reco selection or not.
  float n_truth_hs10, lead_truth_hs_pt;
  // For the leading (highest-pT) truth HS jet specifically: does a RAW reco
  // jet -- from the full AntiKt4EMTopoJets array, before the pT>30 cut and
  // BEFORE lepton-overlap removal -- sit within dR<0.3 of it? And if so, was
  // that raw jet the one OR stripped? Answers whether an "unmatched hard
  // truth jet" (Z+jets only) is jet-finding failing to place a jet there at
  // all, versus OR removing the jet that would have matched.
  float lead_truth_raw_match, lead_truth_raw_or_removed;
  // What the active --jvt working point removed from the pT-passing list, by
  // truth identity. All zero under --jvt=none.
  float n_rm_jvt, n_rm_jvt_hs, n_rm_jvt_pu, n_rm_fjvt, n_rm_fjvt_hs, n_rm_fjvt_pu;
  PairBlock all, acc, wide;
};

// One row per pT-passing jet BEFORE the working-point cut, so thresholds can
// be derived or scanned offline against the same population the cut sees.
struct JetRow {
  float file_idx, entry, pt, abseta, hs, pu;
  float rpt, corrjvf, fjvt, n_ghost, n_ghost_pv, jvt_win, fjvt_win, removed;
};

}  // namespace

int main(int argc, char** argv) {
  SampleConfig cfg = resolveSample(argc, argv);
  ENERGY_LABEL    = cfg.energyLabel;
  OUTPUT_DIR      = cfg.outputDir;
  SAMPLE_NAME     = cfg.sampleName;
  OVERLAP_REMOVAL = cfg.overlapRemoval;
  // --no-or: diagnostic-only override, so region composition can be measured
  // with lepton-jet overlap removal switched off. Does NOT touch OUTPUT_DIR
  // or the output filename -- rerun with and without and compare the two
  // trees directly rather than relying on naming to keep them apart.
  bool noOR = false;
  for (int i = 1; i < argc; ++i) if (std::string(argv[i]) == "--no-or") noOR = true;
  if (noOR) OVERLAP_REMOVAL = false;
  // --jvt=none|loose|tight: pileup-jet tagging applied to the pT-passing list
  // before any pairing (see the file header). Tags the output file.
  const JvtWP* wpSel = &JVT_WPS[0];
  for (int i = 1; i < argc; ++i) {
    const std::string a = argv[i];
    if (a.rfind("--jvt=", 0) != 0) continue;
    const std::string w = a.substr(6);
    const JvtWP* hit = nullptr;
    for (const auto& wp : JVT_WPS) if (w == wp.name) hit = &wp;
    if (!hit) { std::cerr << "[diag] unknown --jvt=" << w << " (none|loose|tight)\n"; return 1; }
    wpSel = hit;
  }
  const JvtWP& WP = *wpSel;
  // The discriminants need the vertex fit's track assignment, which lives in
  // the extended branch set. Must be set before the BranchPointerWrapper binds.
  EXTENDED_BRANCHES = true;
  const Long64_t maxEvents = resolveMaxEvents(argc, argv);
  // MUST be set explicitly -- setupChain reads MyUtl::FILE_SHARD, and nothing
  // populates it as a side effect of resolveSample. Omitting this does not
  // error: --file-shard is silently ignored and every shard reads the WHOLE
  // sample, so a sharded run quietly does N times the work and, if the outputs
  // are merged, produces N copies of every row. Found exactly that way.
  MyUtl::FILE_SHARD = MyUtl::resolveShard(argc, argv);

  // As vbs_region_mjj: strip the topology selection so the diagnostic sees the
  // whole paired population rather than redrawing its own cut.
  VBS_JET_MJJ   = 0.0;
  VBS_JET_D_ETA = 0.0;

  std::string outPath = OUTPUT_DIR + "/" +
                        (SAMPLE_NAME.empty() ? std::string("local") : SAMPLE_NAME) +
                        (noOR ? "_noOR" : "") + WP.tag + "_vbs_region_diag.root";
  for (int i = 1; i < argc; ++i) {
    std::string a = argv[i];
    if (a.rfind("--out=", 0) == 0) outPath = a.substr(6);
  }
  boost::filesystem::create_directories(OUTPUT_DIR);
  std::cout << "[diag] sample=" << (SAMPLE_NAME.empty() ? "local" : SAMPLE_NAME)
            << (noOR ? " (OR disabled)" : "")
            << " jvt=" << WP.name;
  if (WP.rptMin >= 0.0)
    std::cout << " (R_pT >= " << WP.rptMin << " for |eta| < " << JVT_ETA_MAX
              << ", pT < " << JVT_PT_MAX << "; fJVT <= " << WP.fjvtMax << " for "
              << FJVT_ETA_MIN << " <= |eta| < " << FJVT_ETA_MAX << ", pT < " << FJVT_PT_MAX << ")";
  std::cout << " out=" << outPath << "\n";

  TChain chain("ntuple");
  setupChain(chain, cfg.ntupleDir.c_str(), MyUtl::FILE_SHARD);
  TTreeReader reader(&chain);
  BranchPointerWrapper branch(reader);
  if (!branch.trackRecoVtxIdx) {
    std::cerr << "[diag] Track_recoVtx_idx is not in this sample -- the JVT/fJVT "
                 "discriminants need the vertex fit's assignment; a silent fallback "
                 "would give every jet R_pT = 0 and remove them all\n";
    return 1;
  }

  TFile out(outPath.c_str(), "RECREATE");
  TTree tree("events", "per-event VBS pair composition diagnostic");
  Row R{};
#define BR(n) tree.Branch(#n, &R.n)
  BR(file_idx); BR(entry);
  BR(n_all_jets); BR(n_fwd_jets); BR(n_beyond_jets);
  BR(n_fwd_hs); BR(n_fwd_pu); BR(n_cen_hs); BR(n_cen_pu); BR(n_any_hs); BR(n_neither);
  BR(n_beyond_hs); BR(n_beyond_pu);
  BR(n_fwd_hs_trk); BR(lead_pt); BR(lead_abseta);
  BR(n_truth_hs10); BR(lead_truth_hs_pt);
  BR(lead_truth_raw_match); BR(lead_truth_raw_or_removed);
  BR(n_rm_jvt); BR(n_rm_jvt_hs); BR(n_rm_jvt_pu);
  BR(n_rm_fjvt); BR(n_rm_fjvt_hs); BR(n_rm_fjvt_pu);
#undef BR
  TTree jtree("jets", "per-jet JVT/fJVT discriminants, every pT-passing jet before the cut");
  JetRow J{};
#define JB(n) jtree.Branch(#n, &J.n)
  JB(file_idx); JB(entry); JB(pt); JB(abseta); JB(hs); JB(pu);
  JB(rpt); JB(corrjvf); JB(fjvt); JB(n_ghost); JB(n_ghost_pv);
  JB(jvt_win); JB(fjvt_win); JB(removed);
#undef JB
  // Two identical blocks, distinguished only by their branch prefix.
  auto branchBlock = [&](const char* pfx, PairBlock& B) {
    auto nm = [&](const char* f) { return std::string(pfx) + f; };
#define BB(f) tree.Branch(nm(#f).c_str(), &B.f)
    BB(n_jets); BB(pair_mjj); BB(pair_deta);
    BB(legA_zone); BB(legA_hs); BB(legA_pu); BB(legA_pt); BB(legA_abseta);
    BB(legB_zone); BB(legB_hs); BB(legB_pu); BB(legB_pt); BB(legB_abseta);
    BB(alt_r1_mjj); BB(alt_r2_mjj); BB(alt_pupu_mjj); BB(alt_bothhs_mjj);
    BB(alt_r2_fwd_abseta); BB(alt_r2_fwd_pt); BB(alt_r2_cen_abseta); BB(alt_r2_cen_pt);
    BB(alt_r1_hs_abseta);  BB(alt_r1_hs_pt);  BB(alt_r1_pu_abseta);  BB(alt_r1_pu_pt);
    BB(legA_rpt); BB(legA_fjvt); BB(legB_rpt); BB(legB_fjvt);
#undef BB
  };
  branchBlock("all_",  R.all);
  branchBlock("acc_",  R.acc);
  branchBlock("wide_", R.wide);

  long nSeen = 0, nSel = 0, nNoPairAll = 0, nNoPairAcc = 0, nDisagree = 0;
  long nJetsPre = 0, nRmJvt = 0, nRmFjvt = 0, nEvtLostToWP = 0;

  while (reader.Next()) {
    if (maxEvents > 0 && nSeen >= maxEvents) break;
    ++nSeen;

    // Same preamble as vbs_region_mjj, so the two diagnostics describe the same
    // event population.
    branch.computeOverlapRemoval();
    if (!branch.passLeptonSelection()) continue;
    if (branch.vetoLeptonOverlap())    continue;
    if (!branch.passBasicCuts())       continue;

    std::vector<int> passPtIdx;
    int nPt = 0, nPtEta = 0;
    branch.collectPtPassingJets(passPtIdx, nPt, nPtEta);
    if (nPt    < MIN_PASSPT_JETS)  continue;
    if (nPtEta < MIN_PASSETA_JETS) continue;

    R = Row{};
    const std::string fp = reader.GetTree()->GetCurrentFile()->GetName();
    { size_t p = fp.rfind('_');
      size_t d = (p == std::string::npos) ? std::string::npos : fp.find('.', p);
      if (p != std::string::npos && d != std::string::npos)
        R.file_idx = (float)std::atoi(fp.substr(p + 1, d - p - 1).c_str()); }
    R.entry = (float)reader.GetTree()->GetReadEntry();

    // ── JVT / fJVT: jets tree on the PRE-cut list, then the cut itself ──────
    // The event-level jet requirements are re-applied to the post-cut list: an
    // analysis counts only tagger-passing jets, so an event whose second jet
    // the tagger removes has no pair to classify. Under --jvt=none both checks
    // are the same check and nothing changes.
    const std::vector<JetDisc> D = computeJetDiscriminants(branch);
    std::vector<int> keptIdx; keptIdx.reserve(passPtIdx.size());
    int nPtWP = 0, nPtEtaWP = 0;
    for (int j : passPtIdx) {
      const JetDisc& d = D[j];
      const double eta = branch.topoJetEta[j], phi = branch.topoJetPhi[j];
      const double aeta = std::abs(eta);
      const bool hs = branch.isJetPaperHS(eta, phi);
      const bool pu = branch.isJetPaperPU(eta, phi);
      const bool rmJvt  = d.jvtWin  && d.rpt  <  WP.rptMin;
      const bool rmFjvt = d.fjvtWin && d.fjvt >  WP.fjvtMax;   // windows are disjoint in |eta|
      J = JetRow{};
      J.file_idx = R.file_idx; J.entry = R.entry;
      J.pt = branch.topoJetPt[j]; J.abseta = (float)aeta;
      J.hs = hs ? 1.f : 0.f; J.pu = pu ? 1.f : 0.f;
      J.rpt = d.rpt; J.corrjvf = d.corrjvf; J.fjvt = d.fjvt;
      J.n_ghost = (float)d.nGhost; J.n_ghost_pv = (float)d.nGhostPV;
      J.jvt_win = d.jvtWin ? 1.f : 0.f; J.fjvt_win = d.fjvtWin ? 1.f : 0.f;
      J.removed = (rmJvt || rmFjvt) ? 1.f : 0.f;
      jtree.Fill();
      ++nJetsPre;
      if (rmJvt)  { ++nRmJvt;  ++R.n_rm_jvt;  if (hs) ++R.n_rm_jvt_hs;  if (pu) ++R.n_rm_jvt_pu;  continue; }
      if (rmFjvt) { ++nRmFjvt; ++R.n_rm_fjvt; if (hs) ++R.n_rm_fjvt_hs; if (pu) ++R.n_rm_fjvt_pu; continue; }
      keptIdx.push_back(j);
      ++nPtWP;
      if (aeta > MIN_ABS_ETA_JET && aeta < MAX_ABS_ETA_JET) ++nPtEtaWP;
    }
    if (nPtWP < MIN_PASSPT_JETS || nPtEtaWP < MIN_PASSETA_JETS) { ++nEvtLostToWP; continue; }
    passPtIdx.swap(keptIdx);

    // ── Event-wide paper-label content ──────────────────────────────────────
    // The same loop classifyEventRegion runs, kept here rather than called so
    // the per-zone pileup count and the "neither" count are available;
    // classifyEventRegion returns only the four its ladder needs.
    for (int j : passPtIdx) {
      const double eta = branch.topoJetEta[j], phi = branch.topoJetPhi[j];
      const bool hs = branch.isJetPaperHS(eta, phi);
      const bool pu = branch.isJetPaperPU(eta, phi);
      const int  z  = zoneOf(std::abs(eta), MAX_ABS_ETA_JET);
      if (branch.topoJetPt[j] > R.lead_pt) {
        R.lead_pt = branch.topoJetPt[j]; R.lead_abseta = (float)std::abs(eta);
      }
      if      (hs && z == 1) ++R.n_fwd_hs;
      else if (hs && z == 0) ++R.n_cen_hs;
      if      (pu && z == 1) ++R.n_fwd_pu;
      else if (pu && z == 0) ++R.n_cen_pu;
      if (hs)        ++R.n_any_hs;
      if (!hs && !pu) ++R.n_neither;
      if (z == 1)    ++R.n_fwd_jets;
      if (z == 2) {
        ++R.n_beyond_jets;
        if (hs) ++R.n_beyond_hs;
        if (pu) ++R.n_beyond_pu;
      }
    }
    R.n_all_jets = (float)passPtIdx.size();

    int leadTruthIdx = -1;
    for (int t = 0; t < (int)branch.truthHSJetPt.GetSize(); ++t) {
      const float tpt = branch.truthHSJetPt[t];
      if (tpt <= 10.0f) continue;
      ++R.n_truth_hs10;
      if (tpt > R.lead_truth_hs_pt) { R.lead_truth_hs_pt = tpt; leadTruthIdx = t; }
    }
    // Raw-match test for the leading truth HS jet: scan EVERY reco jet in the
    // full array (not passPtIdx), so this is independent of both the pT cut
    // and lepton-overlap removal -- isJetRemoved is checked separately below,
    // on whichever raw jet actually matches.
    if (leadTruthIdx >= 0) {
      const double tEta = branch.truthHSJetEta[leadTruthIdx];
      const double tPhi = branch.truthHSJetPhi[leadTruthIdx];
      for (int j = 0; j < (int)branch.topoJetPt.GetSize(); ++j) {
        if (branch.topoJetPt[j] <= MIN_JET_PT) continue;
        const double deta = branch.topoJetEta[j] - tEta;
        const double dphi = TVector2::Phi_mpi_pi(branch.topoJetPhi[j] - tPhi);
        if (std::hypot(deta, dphi) < 0.3) {
          R.lead_truth_raw_match = 1.0f;
          if (branch.isJetRemoved(j)) R.lead_truth_raw_or_removed = 1.0f;
          break;
        }
      }
    }

    for (int t = 0; t < (int)branch.trackPt.GetSize(); ++t) {
      if (branch.trackToTruthvtx[t] != 0)    continue;
      if (branch.trackPt[t] <= MIN_TRACK_PT) continue;
      const double ae = std::abs((double)branch.trackEta[t]);
      if (ae > MIN_ABS_ETA_JET && ae < MAX_ABS_ETA_JET) ++R.n_fwd_hs_trk;
    }

    // ── One jet collection in, one filled PairBlock out ─────────────────────
    auto fillBlock = [&](const std::vector<int>& idx, double fwdMax,
                         PairBlock& B) -> bool {
      B = PairBlock{};
      B.n_jets = (float)idx.size();
      B.pair_mjj = B.pair_deta = -1.f;
      B.alt_r1_mjj = B.alt_r2_mjj = B.alt_pupu_mjj = B.alt_bothhs_mjj = -1.f;
      B.alt_r2_fwd_abseta = B.alt_r2_cen_abseta = -1.f;
      B.alt_r1_hs_abseta  = B.alt_r1_pu_abseta  = -1.f;

      std::vector<char> hs(idx.size()), pu(idx.size());
      std::vector<int>  zn(idx.size());
      for (size_t k = 0; k < idx.size(); ++k) {
        const double eta = branch.topoJetEta[idx[k]], phi = branch.topoJetPhi[idx[k]];
        hs[k] = branch.isJetPaperHS(eta, phi) ? 1 : 0;
        pu[k] = branch.isJetPaperPU(eta, phi) ? 1 : 0;
        zn[k] = zoneOf(std::abs(eta), fwdMax);
      }

      int bestA = -1, bestB = -1;
      double bestM = -1.0;
      for (size_t a = 0; a < idx.size(); ++a) {
        for (size_t b = a + 1; b < idx.size(); ++b) {
          const int i = idx[a], j = idx[b];
          const float ei = branch.topoJetEta[i], ej = branch.topoJetEta[j];
          if (ei * ej >= 0) continue;             // opposite hemispheres
          TLorentzVector vi, vj;
          vi.SetPtEtaPhiM(branch.topoJetPt[i], ei, branch.topoJetPhi[i], 0.0);
          vj.SetPtEtaPhiM(branch.topoJetPt[j], ej, branch.topoJetPhi[j], 0.0);
          const double m = (vi + vj).M();
          if (m > bestM) {
            bestM = m; bestA = (int)a; bestB = (int)b;
          }
          const bool r1 = (zn[a] == 1 && zn[b] == 1) &&
                          ((hs[a] && pu[b]) || (hs[b] && pu[a]));
          const bool r2 = (zn[a] == 1 && pu[a] && zn[b] == 0 && hs[b]) ||
                          (zn[b] == 1 && pu[b] && zn[a] == 0 && hs[a]);
          if (r1 && m > B.alt_r1_mjj) {
            B.alt_r1_mjj = (float)m;
            const bool iIsHS = hs[a] && pu[b];
            const int jhs = iIsHS ? i : j, jpu = iIsHS ? j : i;
            B.alt_r1_hs_abseta = (float)std::abs((double)branch.topoJetEta[jhs]);
            B.alt_r1_hs_pt     = branch.topoJetPt[jhs];
            B.alt_r1_pu_abseta = (float)std::abs((double)branch.topoJetEta[jpu]);
            B.alt_r1_pu_pt     = branch.topoJetPt[jpu];
          }
          if (r2 && m > B.alt_r2_mjj) {
            B.alt_r2_mjj = (float)m;
            const bool iIsFwd = (zn[a] == 1 && pu[a]);
            const int jfwd = iIsFwd ? i : j, jcen = iIsFwd ? j : i;
            B.alt_r2_fwd_abseta = (float)std::abs((double)branch.topoJetEta[jfwd]);
            B.alt_r2_fwd_pt     = branch.topoJetPt[jfwd];
            B.alt_r2_cen_abseta = (float)std::abs((double)branch.topoJetEta[jcen]);
            B.alt_r2_cen_pt     = branch.topoJetPt[jcen];
          }
          if (pu[a] && pu[b] && m > B.alt_pupu_mjj)   B.alt_pupu_mjj   = (float)m;
          if (hs[a] && hs[b] && m > B.alt_bothhs_mjj) B.alt_bothhs_mjj = (float)m;
        }
      }
      if (bestA < 0) return false;

      auto leg = [&](size_t k, float& z, float& h, float& p, float& pt, float& ae,
                     float& rpt, float& fjvt) {
        z = (float)zn[k]; h = hs[k] ? 1.f : 0.f; p = pu[k] ? 1.f : 0.f;
        pt = branch.topoJetPt[idx[k]];
        ae = (float)std::abs((double)branch.topoJetEta[idx[k]]);
        rpt = D[idx[k]].rpt; fjvt = D[idx[k]].fjvt;
      };
      leg(bestA, B.legA_zone, B.legA_hs, B.legA_pu, B.legA_pt, B.legA_abseta,
          B.legA_rpt, B.legA_fjvt);
      leg(bestB, B.legB_zone, B.legB_hs, B.legB_pu, B.legB_pt, B.legB_abseta,
          B.legB_rpt, B.legB_fjvt);
      B.pair_mjj  = (float)bestM;
      B.pair_deta = (float)std::abs((double)branch.topoJetEta[idx[bestA]] -
                                    (double)branch.topoJetEta[idx[bestB]]);
      return true;
    };

    // Block 1: every pT-passing jet -- reproduces calcBestVbsPair exactly.
    const bool okAll = fillBlock(passPtIdx, MAX_ABS_ETA_JET, R.all);
    if (!okAll) ++nNoPairAll;

    // Block 2: analysis-quality tagging jets only -- pT-passing AND inside the
    // jet acceptance. A jet beyond MAX_ABS_ETA_JET is not something the
    // analysis tags on and is not something HGTD can time, so letting one win
    // the m_jj ranking describes a pair the analysis would never form.
    std::vector<int> accIdx;
    for (int j : passPtIdx)
      if (std::abs((double)branch.topoJetEta[j]) < MAX_ABS_ETA_JET) accIdx.push_back(j);
    const bool okAcc = fillBlock(accIdx, MAX_ABS_ETA_JET, R.acc);
    if (!okAcc) ++nNoPairAcc;

    // Block 3: every pT-passing jet again, but with no upper edge on
    // "forward". Same pair the all_ block picks -- the m_jj ranking does not
    // consult the zones -- so any difference between all_ and wide_ is purely
    // the region LABELLING, not a different pair.
    const bool okWide = fillBlock(passPtIdx, FWD_MAX_UNBOUNDED, R.wide);

    if (!okAll && !okAcc) continue;
    ++nSel;

    // Cross-check against the shared classifier: what the clustering side
    // actually fills must match what this
    // reconstructs from the legs, or the diagnostic describes another region.
    // Cross-check against the SHARED classifier. That classifier now runs at
    // VBS_FWD_ETA_MAX, so the block it must agree with is wide_, not all_ --
    // pointing it at all_ would report a disagreement on every event whose
    // pair holds a leg past 4.00, which is exactly the population this
    // diagnostic exists to study.
    const bool locR1 = (R.wide.legA_zone == 1 && R.wide.legB_zone == 1) &&
                       ((R.wide.legA_hs && R.wide.legB_pu) ||
                        (R.wide.legB_hs && R.wide.legA_pu));
    const bool locR2 = (R.wide.legA_zone == 1 && R.wide.legA_pu &&
                        R.wide.legB_zone == 0 && R.wide.legB_hs) ||
                       (R.wide.legB_zone == 1 && R.wide.legB_pu &&
                        R.wide.legA_zone == 0 && R.wide.legA_hs);
    // Only meaningful under --jvt=none: the shared classifier pairs from the
    // UNTAGGED jet list, so under a working point it disagrees on exactly the
    // events whose pair the tagger changed -- which is the measurement, not a
    // bug. Gated rather than reinterpreted so the counter keeps meaning "bug".
    if (WP.rptMin < 0.0) {
      const VbsRegion shared =
          branch.classifyVbsRegion(MIN_ABS_ETA_JET, VBS_FWD_ETA_MAX, MIN_ABS_ETA_JET);
      if (okWide && (locR1 != (shared == VbsRegion::R1) ||
                     locR2 != (shared == VbsRegion::R2))) ++nDisagree;
    }

    tree.Fill();
  }

  std::cout << "\n[diag] seen " << nSeen << ", selected " << nSel
            << ", no pair (all jets) " << nNoPairAll
            << ", no pair (acceptance jets) " << nNoPairAcc
            << ", R1/R2 cross-check disagreements "
            << (WP.rptMin < 0.0 ? std::to_string(nDisagree) : std::string("n/a under a working point")) << "\n";
  std::cout << "[diag] jvt=" << WP.name << ": " << nJetsPre << " pT-passing jets, removed "
            << nRmJvt << " by JVT and " << nRmFjvt << " by fJVT; "
            << nEvtLostToWP << " events dropped below the jet requirements\n";
  tree.Write();
  jtree.Write();
  out.Close();
  std::cout << "[diag] wrote " << outPath << "\n";
  return 0;
}
