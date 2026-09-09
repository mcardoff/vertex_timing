// -----------------------------------------------------------------------------
// skim_ntuples -- write a skimmed + slimmed copy of a sample, once, so every
// later pass reads a fraction of the bytes and a fraction of the files.
//
// Why this exists. The AF ntuples live on CephFS -- /data/mcardiff is a symlink
// to /cold/user/mcardiff, the cold tier -- and the cost of a full pass is
// dominated by per-file read round trips, not by CPU or bandwidth. Measured on
// zjets: 2000 files x 306 MB = 612 GB, 306 KB per event, 195 branches, of which
// this code binds 62; a worker issues ~155 read syscalls per file at ~45 KB
// each, which at ceph latency IS the 0.3-1.2 s/file recorded in CLAUDE.md. Two
// independent factors are therefore worth attacking, and this attacks both:
//
//   SKIM   drop events no executable can use. Every consumer applies
//          passBasicCuts() and passLeptonSelection() before doing anything, so
//          an event failing either is read and discarded on every single pass.
//          On zjets that alone is 73.5% of the sample (530,962 of 2,000,000
//          survive the Z->ll requirement).
//
//   SLIM   drop branches nothing binds. BranchPointerWrapper binds 37 branches
//          always, Track_leptonID under OVERLAP_REMOVAL, and 24 more under
//          EXTENDED_BRANCHES (export_training_data only) -- 47.7% of the file's
//          bytes, consistently across all four samples. The other 52% is read
//          past on every pass and never used.
//
// ...and a third, implicitly: output is chunked so one output file replaces
// many input files, which is a direct multiplier on the per-file round trips
// that dominate the cost in the first place.
//
// SELECTION IS THE LOOSEST COMMON DENOMINATOR, not "the" event selection. Only
// a cut EVERY consumer applies may go in here, because an event dropped by the
// skim is invisible to all of them forever. That is:
//
//     a reco and a truth vertex exist, AND
//     |TruthVtx_z[0] - RecoVtx_z[0]| < MAX_VTX_DZ, AND
//     passLeptonSelection()            (a no-op unless OVERLAP_REMOVAL)
//
// Note what is NOT here. passBasicCuts() ALSO requires >= MIN_JETS jets, and
// skimming on that is wrong: util/rpt_v5_hist.cxx deliberately does not apply
// it ("every forward jet in the acceptance contributes an independent R_pT
// measurement"), so the jet cut silently removed 20 of its events on a 3-file
// local test and shifted 28 of its histograms. It is only 0.26% of events, and
// consumers that want it apply it themselves. Nothing topological is applied
// either: the VBS knobs (--vbs-mjj/--vbs-deta) are a per-study choice and
// freezing one study's topology into every other study's input is exactly the
// trap this comment exists to prevent. vetoLeptonOverlap() is likewise left
// out -- currently universal, cheap downstream, and not worth the coupling.
//
// The lepton and vertex cuts ARE evaluated through BranchPointerWrapper rather
// than reimplemented, so those cannot drift from what the analysis does.
//
// WHAT THIS COSTS YOU. Event counts change: n_total after the skim is the
// skimmed count, not the generated count, so any "N% of events pass" figure
// computed against it is a fraction of the skim, not of the sample. Each output
// file therefore carries the pre-skim counts as TParameters (skim_n_input,
// skim_n_pass_jets, skim_n_pass_vtx, skim_n_pass_lepton, skim_n_output) so the
// original denominators stay recoverable. And a study that later needs a
// dropped branch has to go back to the full sample -- the keep list below is
// the contract.
// -----------------------------------------------------------------------------
#include <TChain.h>
#include <TEntryList.h>
#include <TFile.h>
#include <TParameter.h>
#include <TTree.h>
#include <TTreeReader.h>

#include <algorithm>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

#include "clustering_constants.h"
#include "sample_config.h"
#include "clustering_includes.h"
#include "clustering_structs.h"
#include "event_processing.h"

using namespace MyUtl;

namespace {

// The keep list IS the contract with BranchPointerWrapper -- these are exactly
// the branches its constructor binds. Anything added there must be added here
// or the skimmed sample will fail at first access with a missing branch.
// Grouped the same way the constructor groups them so the two can be diffed by
// eye.
const std::vector<std::string> KEEP_ALWAYS = {
  "weight",
  "Track_d0", "Track_z0", "Track_pt", "Track_eta", "Track_qOverP",
  "Track_theta", "Track_phi",
  "Track_var_d0", "Track_var_z0", "Track_var_qOverP", "Track_var_theta",
  "Track_var_phi0",
  "Track_time", "Track_timeRes", "Track_hasValidTime", "Track_quality",
  "Track_truthVtx_idx", "Track_truthPart_idx",
  "Track_nHGTDHits", "Track_nHGTDPrimaryHits",
  "TruthVtx_z", "TruthVtx_time", "TruthVtx_isHS",
  "RecoVtx_z", "RecoVtx_time", "RecoVtx_timeRes", "RecoVtx_hasValidTime",
  "AntiKt4EMTopoJets_pt", "AntiKt4EMTopoJets_eta", "AntiKt4EMTopoJets_phi",
  "TruthHSJet_pt", "TruthHSJet_eta", "TruthHSJet_phi",
  "AntiKt4EMTopoJets_truthHSJet_idx", "AntiKt4EMTopoJets_ghostTrack_idx",
  "TruthPart_prodVtx_time",
};

// Bound only when OVERLAP_REMOVAL is set, but kept unconditionally: the skimmed
// sample must not depend on which flag the SKIM happened to run with.
const std::vector<std::string> KEEP_LEPTON = { "Track_leptonID" };

// EXTENDED_BRANCHES (export_training_data). Kept by default -- a skim that
// silently breaks the ML exporter is a trap -- at a cost of 47.7% of the file
// against 31.1% without them. --no-extended drops them if space is tight.
const std::vector<std::string> KEEP_EXTENDED = {
  "Track_recoVtx_idx", "Track_recoVtx_weight", "RecoVtx_sumPt2",
  "RecoVtx_isHS", "RecoVtx_isPU",
  "Track_chi2", "Track_ndof",
  "Track_numberOfPixelHits", "Track_numberOfPixelHoles",
  "Track_numberOfPixelSharedHits", "Track_numberOfPixelSplitHits",
  "Track_numberOfSCTHits", "Track_numberOfSCTHoles", "Track_numberOfSCTSharedHits",
  "Track_numberOfInnermostPixelLayerHits",
  "Track_numberOfNextToInnermostPixelLayerHits",
  "Track_btagIp_d0", "Track_btagIp_d0Uncertainty",
  "Track_btagIp_z0SinTheta", "Track_btagIp_z0SinThetaUncertainty",
  "Track_cov_d0z0", "Track_cov_z0theta", "Track_cov_z0phi0", "Track_cov_z0qOverP",
};

struct Counts {
  Long64_t input = 0, passJets = 0, passVtx = 0, passLepton = 0, output = 0;
};

}  // namespace

int main(int argc, char** argv) {
  SampleConfig cfg = resolveSample(argc, argv);
  SAMPLE_NAME     = cfg.sampleName;
  OVERLAP_REMOVAL = cfg.overlapRemoval;
  const Long64_t maxEvents = resolveMaxEvents(argc, argv);
  // MUST be set explicitly -- setupChain reads MyUtl::FILE_SHARD, and nothing
  // populates it as a side effect of resolveSample. Omitting this does not
  // error: --file-shard is silently ignored and every shard reads the WHOLE
  // sample, so a sharded run quietly does N times the work and, if the outputs
  // are merged, produces N copies of every row. Found exactly that way.
  MyUtl::FILE_SHARD = MyUtl::resolveShard(argc, argv);

  std::string outDir = "/data/mcardiff/skimmed_ntuples/" +
                       (SAMPLE_NAME.empty() ? std::string("local") : SAMPLE_NAME);
  int  filesPerOutput = 40;
  bool withExtended   = true;
  for (int i = 1; i < argc; ++i) {
    std::string a = argv[i];
    if      (a.rfind("--out-dir=", 0) == 0)         outDir = a.substr(10);
    else if (a.rfind("--files-per-output=", 0) == 0) filesPerOutput = std::stoi(a.substr(19));
    else if (a == "--no-extended")                   withExtended = false;
  }

  std::vector<std::string> keep = KEEP_ALWAYS;
  keep.insert(keep.end(), KEEP_LEPTON.begin(), KEEP_LEPTON.end());
  if (withExtended) keep.insert(keep.end(), KEEP_EXTENDED.begin(), KEEP_EXTENDED.end());

  // Build the input file list exactly as setupChain does -- same sort, same
  // nested-directory walk, same --file-shard semantics -- so a sharded skim
  // partitions the sample the same way a sharded analysis run does.
  TChain probe("ntuple");
  setupChain(probe, cfg.ntupleDir.c_str(), MyUtl::FILE_SHARD);
  std::vector<std::string> files;
  if (auto* fl = probe.GetListOfFiles())
    for (int i = 0; i < fl->GetEntries(); ++i)
      files.push_back(fl->At(i)->GetTitle());

  // ---------------------------------------------------------------------
  // Refuse to skim a skim.
  //
  // The registry now resolves vbf/zjets/dijet/ttbar to their SKIMMED copies
  // (see sample_config.h), so a bare `skim_ntuples --sample=zjets` would read
  // skimmed_ntuples/zjets and write straight back into it -- reading and
  // rewriting the same directory, with the second selection pass silently
  // applied on top of the first. Re-skimming needs an explicit
  // `--ntuple-dir=<original>`.
  //
  // Checked two ways, because either alone has a hole: the path test misses an
  // out-dir that merely aliases the input, and the marker test misses an empty
  // or unreadable first file.
  {
    namespace fs = boost::filesystem;
    boost::system::error_code ec1, ec2;
    const fs::path inCanon  = fs::weakly_canonical(fs::path(cfg.ntupleDir), ec1);
    const fs::path outCanon = fs::weakly_canonical(fs::path(outDir), ec2);
    if (!ec1 && !ec2 && inCanon == outCanon) {
      std::cerr << "[skim] REFUSING: input and output are the same directory ("
                << inCanon.string() << ").\n"
                   "       Pass --ntuple-dir=<original sample dir> to re-skim.\n";
      return 1;
    }
    if (!files.empty()) {
      std::unique_ptr<TFile> f(TFile::Open(files.front().c_str(), "READ"));
      if (f && !f->IsZombie() && f->Get("skim_n_input")) {
        std::cerr << "[skim] REFUSING: " << files.front()
                  << " already carries skim_n_input, so this input is itself a "
                     "skim.\n"
                     "       Pass --ntuple-dir=<original sample dir> to "
                     "re-skim from the originals.\n";
        return 1;
      }
    }
  }

  boost::filesystem::create_directories(outDir);
  const std::string tag = MyUtl::FILE_SHARD.active()
      ? "_shard" + std::to_string(MyUtl::FILE_SHARD.index) + "of" +
        std::to_string(MyUtl::FILE_SHARD.count)
      : "";
  std::cout << "[skim] sample=" << (SAMPLE_NAME.empty() ? "local" : SAMPLE_NAME)
            << "  files=" << files.size()
            << "  keep=" << keep.size() << " branches"
            << (withExtended ? " (with extended)" : " (no extended)")
            << "\n[skim] out=" << outDir << "  " << filesPerOutput
            << " input files per output file\n";

  Counts tot;
  int chunk = 0;

  // Work in chunks of filesPerOutput input files, one output file each: the
  // per-file read round trip is the cost this whole executable exists to
  // reduce, so collapsing many inputs into one output is not just tidiness.
  //
  // Each chunk is copied with TChain::SetEntryList + CopyTree(""), ROOT's own
  // selective-copy path. The hand-rolled CloneTree(0) + CopyAddresses loop this
  // replaced segfaulted inside CopyAddresses' undo when the output tree had
  // been cloned from the very tree being undone -- and getting that ownership
  // dance right is exactly what CopyTree already does.
  for (size_t base = 0; base < files.size(); base += (size_t)filesPerOutput) {
    if (maxEvents > 0 && tot.input >= maxEvents) break;
    const size_t hi = std::min(files.size(), base + (size_t)filesPerOutput);
    Counts c;

    TChain chain("ntuple");
    for (size_t i = base; i < hi; ++i) chain.Add(files[i].c_str());

    // ── Pass 1: selection, through BranchPointerWrapper ───────────────────
    // TTreeReader is lazy, so only the branches passBasicCuts /
    // passLeptonSelection actually touch are read. Scoped so the reader's
    // branch addresses are gone before the copy sets its own.
    auto elist = std::make_unique<TEntryList>("skim_list", "skim_list");
    {
      TTreeReader reader(&chain);
      BranchPointerWrapper branch(reader);
      Long64_t e = -1;
      while (reader.Next()) {
        ++e;
        if (maxEvents > 0 && tot.input >= maxEvents) break;
        ++tot.input; ++c.input;
        // Counted, NOT cut on -- see the header. Reported so the cost of
        // adding it later is visible without another pass.
        if ((int)branch.topoJetPt.GetSize() >= MIN_JETS) { ++tot.passJets; ++c.passJets; }
        if (branch.recoVtxZ.GetSize() == 0 || branch.truthVtxZ.GetSize() == 0) continue;
        if (std::abs(branch.truthVtxZ[0] - branch.recoVtxZ[0]) > MAX_VTX_DZ) continue;
        ++tot.passVtx; ++c.passVtx;
        branch.computeOverlapRemoval();
        if (!branch.passLeptonSelection()) continue;
        ++tot.passLepton; ++c.passLepton;
        elist->Enter(e, &chain);
      }
    }
    chain.ResetBranchAddresses();

    // ── Pass 2: copy the kept entries, kept branches only ─────────────────
    // The chunk is in page cache from pass 1, which touched only a handful of
    // small branches, so this second read is the cheap one.
    chain.SetBranchStatus("*", 0);
    for (const auto& b : keep) {
      if (!chain.GetBranch(b.c_str())) {
        // Productions differ -- Track_btagIp_* is absent from every grid
        // sample, and BranchPointerWrapper already tolerates that via
        // hasBranch(). Report once rather than aborting.
        if (base == 0) std::cout << "[skim]   (branch absent in this sample: " << b << ")\n";
        continue;
      }
      chain.SetBranchStatus(b.c_str(), 1);
    }
    chain.SetEntryList(elist.get());

    const std::string path = outDir + "/" +
        (SAMPLE_NAME.empty() ? std::string("local") : SAMPLE_NAME) +
        "_skim" + tag + "_" + std::to_string(chunk++) + ".root";
    std::unique_ptr<TFile> out(TFile::Open(path.c_str(), "RECREATE"));
    out->cd();
    TTree* outTree = chain.CopyTree("");
    c.output = outTree ? outTree->GetEntries() : 0;
    tot.output += c.output;

    // Counters travel with the data, so the pre-skim denominators survive even
    // when the file is copied away from this log.
    auto put = [&](const char* n, Long64_t v) { TParameter<Long64_t> p(n, v); p.Write(); };
    put("skim_n_input",       c.input);
    put("skim_n_pass_jets",   c.passJets);
    put("skim_n_pass_vtx",    c.passVtx);
    put("skim_n_pass_lepton", c.passLepton);
    put("skim_n_output",      c.output);
    if (outTree) outTree->Write("", TObject::kOverwrite);
    const double mb = out->GetSize() / 1e6;
    out.reset();
    chain.SetEntryList(nullptr);

    std::cout << "[skim] files " << base << "-" << (hi - 1) << " -> "
              << path.substr(path.rfind('/') + 1) << "  "
              << c.output << "/" << c.input << " events, " << mb << " MB\n" << std::flush;
  }

  auto pct = [&](Long64_t n) { return tot.input ? 100.0 * n / tot.input : 0.0; };
  std::cout << "\n[skim] input " << tot.input
            << "\n[skim]   (>= " << MIN_JETS << " jets, NOT cut on) " << tot.passJets
            << " (" << pct(tot.passJets) << "%)"
            << "\n[skim]   vertex dz       " << tot.passVtx
            << " (" << pct(tot.passVtx) << "%)"
            << "\n[skim]   lepton sel      " << tot.passLepton
            << " (" << pct(tot.passLepton) << "%)"
            << "\n[skim] written " << tot.output << " events in " << chunk
            << " files (" << pct(tot.output) << "% of input)\n";
  return 0;
}
