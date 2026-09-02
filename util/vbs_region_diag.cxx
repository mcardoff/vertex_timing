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
  PairBlock all, acc, wide;
};

}  // namespace

int main(int argc, char** argv) {
  SampleConfig cfg = resolveSample(argc, argv);
  ENERGY_LABEL    = cfg.energyLabel;
  OUTPUT_DIR      = cfg.outputDir;
  SAMPLE_NAME     = cfg.sampleName;
  OVERLAP_REMOVAL = cfg.overlapRemoval;
  const Long64_t maxEvents = resolveMaxEvents(argc, argv);

  // As vbs_region_mjj: strip the topology selection so the diagnostic sees the
  // whole paired population rather than redrawing its own cut.
  VBS_JET_MJJ   = 0.0;
  VBS_JET_D_ETA = 0.0;

  std::string outPath = OUTPUT_DIR + "/" +
                        (SAMPLE_NAME.empty() ? std::string("local") : SAMPLE_NAME) +
                        "_vbs_region_diag.root";
  for (int i = 1; i < argc; ++i) {
    std::string a = argv[i];
    if (a.rfind("--out=", 0) == 0) outPath = a.substr(6);
  }
  boost::filesystem::create_directories(OUTPUT_DIR);
  std::cout << "[diag] sample=" << (SAMPLE_NAME.empty() ? "local" : SAMPLE_NAME)
            << " out=" << outPath << "\n";

  TChain chain("ntuple");
  setupChain(chain, cfg.ntupleDir.c_str(), MyUtl::FILE_SHARD);
  TTreeReader reader(&chain);
  BranchPointerWrapper branch(reader);

  TFile out(outPath.c_str(), "RECREATE");
  TTree tree("events", "per-event VBS pair composition diagnostic");
  Row R{};
#define BR(n) tree.Branch(#n, &R.n)
  BR(file_idx); BR(entry);
  BR(n_all_jets); BR(n_fwd_jets); BR(n_beyond_jets);
  BR(n_fwd_hs); BR(n_fwd_pu); BR(n_cen_hs); BR(n_cen_pu); BR(n_any_hs); BR(n_neither);
  BR(n_beyond_hs); BR(n_beyond_pu);
  BR(n_fwd_hs_trk); BR(lead_pt); BR(lead_abseta);
#undef BR
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
#undef BB
  };
  branchBlock("all_",  R.all);
  branchBlock("acc_",  R.acc);
  branchBlock("wide_", R.wide);

  long nSeen = 0, nSel = 0, nNoPairAll = 0, nNoPairAcc = 0, nDisagree = 0;

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

      auto leg = [&](size_t k, float& z, float& h, float& p, float& pt, float& ae) {
        z = (float)zn[k]; h = hs[k] ? 1.f : 0.f; p = pu[k] ? 1.f : 0.f;
        pt = branch.topoJetPt[idx[k]];
        ae = (float)std::abs((double)branch.topoJetEta[idx[k]]);
      };
      leg(bestA, B.legA_zone, B.legA_hs, B.legA_pu, B.legA_pt, B.legA_abseta);
      leg(bestB, B.legB_zone, B.legB_hs, B.legB_pu, B.legB_pt, B.legB_abseta);
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
    fillBlock(passPtIdx, FWD_MAX_UNBOUNDED, R.wide);

    if (!okAll && !okAcc) continue;
    ++nSel;

    // Cross-check the unrestricted block against the shared classifier: what
    // rpt_v5 and the clustering side actually fill must match what this
    // reconstructs from the legs, or the diagnostic describes another region.
    const bool locR1 = (R.all.legA_zone == 1 && R.all.legB_zone == 1) &&
                       ((R.all.legA_hs && R.all.legB_pu) ||
                        (R.all.legB_hs && R.all.legA_pu));
    const bool locR2 = (R.all.legA_zone == 1 && R.all.legA_pu &&
                        R.all.legB_zone == 0 && R.all.legB_hs) ||
                       (R.all.legB_zone == 1 && R.all.legB_pu &&
                        R.all.legA_zone == 0 && R.all.legA_hs);
    const VbsRegion shared =
        branch.classifyVbsRegion(MIN_ABS_ETA_JET, MAX_ABS_ETA_JET, MIN_ABS_ETA_JET);
    if (okAll && (locR1 != (shared == VbsRegion::R1) ||
                  locR2 != (shared == VbsRegion::R2))) ++nDisagree;

    tree.Fill();
  }

  std::cout << "\n[diag] seen " << nSeen << ", selected " << nSel
            << ", no pair (all jets) " << nNoPairAll
            << ", no pair (acceptance jets) " << nNoPairAcc
            << ", R1/R2 cross-check disagreements " << nDisagree << "\n";
  tree.Write();
  out.Close();
  std::cout << "[diag] wrote " << outPath << "\n";
  return 0;
}
