// ---------------------------------------------------------------------------
// vbs_region_mjj.cxx — VBS-topology composition of the selection vs m_jj.
//
// The analysis argues that HGTD's gain is concentrated in the VBS topologies,
// but nothing so far said what FRACTION of the selection it can act on at all.
// This answers that: a stacked composition plot, one column per m_jj bin, every
// column normalised to 1, with the three reachable topologies at the bottom.
//
//   R1  both tags forward, one HS + one PU  -- timing says WHICH is real
//   R2  forward PU + central HS             -- timing rejects the fake
//   R3  forward HS + central PU             -- the mirror of R2: the GENUINE
//       tag is the timeable one, so timing confirms rather than rejects
//
// Reading it: the y axis is "given an event at this m_jj, what topology is it?"
// Because each column is normalised independently, the falling m_jj spectrum
// cannot wash out the tail. The raw spectrum is on page 2 so low-statistics
// columns -- whose fractions are noisy however clean they look -- stay
// identifiable.
//
// Denominator: every event the rest of the selection admits -- lepton selection
// / overlap veto / passBasicCuts / >=MIN_PASSPT_JETS above MIN_JET_PT /
// >=MIN_PASSETA_JETS forward -- that HAS an opposite-hemisphere candidate pair,
// with NO |Deta| requirement. The m_jj axis starts at the 500 GeV working point
// and events below it are dropped and counted, so the plot describes the sample
// the analysis actually runs on rather than one it never sees.
//
// Pure diagnostic: reads only, changes no analysis behaviour.
//
// Output: <OUTPUT_DIR>/diagnostics/vbs_region_mjj.pdf   (2 pages)
//         <OUTPUT_DIR>/vbs_region_mjj.root              (raw + normalised)
//         + a console composition table
// ---------------------------------------------------------------------------

#include <TCanvas.h>
#include <TChain.h>
#include <TFile.h>
#include <TH1.h>
#include <TLegend.h>
#include <TLine.h>
#include <TLatex.h>
#include <THStack.h>
#include <TStyle.h>
#include <TTreeReader.h>
#include <ROOT/TTreeProcessorMT.hxx>

#include <boost/filesystem.hpp>
#include <array>
#include <cstdio>
#include <iostream>
#include <string>
#include <atomic>
#include <memory>
#include <mutex>
#include <vector>

#include "clustering_constants.h"
#include "event_processing.h"
#include "plotting_utilities.h"
#include "sample_config.h"
#include "AtlasStyle.h"
#include "AtlasLabels.h"

using namespace MyUtl;

namespace {

  // ── Categories ───────────────────────────────────────────────────────────
  // Ordered to answer one question: in what fraction of events can HGTD
  // contribute at all? The first three are the topologies where it can; the
  // rest are the ways an event puts itself out of reach.
  //
  //   R1  both tags forward, one HS + one PU  -- timing must say WHICH is real
  //   R2  forward PU + central HS             -- timing must reject the fake
  //   R3  forward HS + central PU             -- the mirror of R2: the GENUINE
  //       tag is the timeable one, so timing confirms rather than rejects. A
  //       different use of the same information, and on local VBF roughly 3x
  //       the size of R2, which is why it is broken out rather than buried in a
  //       residual.
  //
  // EXHAUSTIVE by construction: per-column normalisation is meaningless
  // otherwise, so the ladder ends in a catch-all and main() checks the totals.
  //
  // NO_ACC ("neither tag in HGTD acceptance") deliberately sits ABOVE the
  // truth-content categories in the ladder: for the question being asked, an
  // event with both tags out of acceptance is unreachable regardless of whether
  // its tags are truth-HS, and saying so is more informative than sorting it by
  // a truth label HGTD will never get to use. R1/R2/R3 all require a forward
  // leg, so they cannot collide with it.
  //
  // OTHER is not padding. isJetPaperHS and isJetPaperPU are NOT complements --
  // HS wants a match within dR 0.3 of a truth HS jet above 10 GeV, PU wants NO
  // match within dR 0.6 above 4 GeV -- so a jet between the two cones is
  // neither, and those events have to land somewhere. It also holds the mixed
  // pairs whose second leg is beyond |eta| = MAX_ABS_ETA_JET.
  enum Cat { R1 = 0, R2, CAN_HELP, NO_T0, MAY_NOT, NCAT };

  const std::array<const char*, NCAT> CAT_KEY = {
    "r1", "r2", "can_help", "no_t0", "may_not"
  };

  const std::array<const char*, NCAT> CAT_LABEL = {
    "R1: fwd HS + fwd PU",
    "R2: fwd PU + central HS",
    "Other: HGTD can help",
    "fwd PU only: t_{0} set by the fake",
    "Other: HGTD may not help"
  };

  // The project palette, in registry order.
  const std::array<Color_t, NCAT> CAT_COLOR = { C01, C02, C03, C05, C04 };

  // ── m_jj binning ─────────────────────────────────────────────────────────
  // Variable width: narrow where the spectrum is steep and every sample has
  // events, widening into the tail so Z+jets columns are not built from a
  // handful of events. 500 is deliberately an edge -- it is the analysis working
  // point, marked with a line on the plot, so the columns to its left are
  // exactly the events the cut throws away.
  const std::vector<double> MJJ_EDGES = {
    0, 250, 500, 750, 1000, 1250, 1500, 1750, 2000, 2500, 3000, 4000
  };

  // ── Region eta window ────────────────────────────────────────────────────
  // Uses clustering_constants.h's MIN_ABS_ETA_JET / MAX_ABS_ETA_JET (2.38-4.00),
  // the same window event_processing.h classifies the VBF_R1 / VBF_R2 scores
  // with -- see the classifyVbsRegion call in main().
  //
  // Worth knowing when comparing against the R_pT study: util/rpt_v5_hist.cxx
  // classifies the SAME regions with its own 2.4 < |eta| < 3.8 instead, and the
  // difference is not cosmetic. On local VBF above m_jj 500 this window calls
  // 2348 events R1 where rpt_v5's calls 1830 -- a factor 1.3 -- because it also
  // admits pairs with a leg between 3.8 and 4.0. So the R1/R2 fractions here are
  // the clustering analysis's regions, and will read HIGHER than the region
  // yields rpt_v5_plot prints for the same sample. That is a real difference in
  // definition between the two halves of the analysis, not a bug in either.

  // Fold overflow into the top bin so the multi-TeV VBF tail is counted rather
  // than silently dropped from a column that would otherwise never sum to 1.
  double folded(const TH1D* h, double v) {
    return std::min(v, h->GetXaxis()->GetBinCenter(h->GetNbinsX()));
  }

// -----------------------------------------------------------------------------
// ThreadState -- one per TTreeProcessorMT worker, merged after the loop.
//
// Same hand-rolled mutex-guarded-registry + thread_local-pointer-cache pattern
// as util/rpt_v5_hist.cxx and src/clustering_hist.cxx, and for the same reason:
// this struct is trivially copy-constructible via its raw TH1D*, so
// ROOT::TThreadedObject would silently SHALLOW-copy the pointers into every
// "cloned" slot and have all workers Fill() the same histograms unsynchronised.
//
// Histograms are constructed here rather than in main() so each worker owns its
// own set; TH1::AddDirectory(kFALSE) in main() keeps the identically-named sets
// from racing on ROOT's global directory registry.
// -----------------------------------------------------------------------------
struct ThreadState {
  std::array<TH1D*, NCAT> raw{};
  TH1D* hTot = nullptr;

  long nSeen = 0, nSel = 0, nNoPair = 0, nBelowMjj = 0;
  std::array<long, NCAT> nCat{};
  long nHelpFwdHS = 0, nHelpNoHS = 0, nNoT0 = 0;
  long nMayNotHsBeyond = 0, nMayNotNoFwd = 0;
  long nAnyFwd = 0;

  ThreadState() {
    const int    nbin = (int)MJJ_EDGES.size() - 1;
    const double* ed  = MJJ_EDGES.data();
    for (int c = 0; c < NCAT; ++c)
      raw[c] = new TH1D((std::string("mjj_raw_") + CAT_KEY[c]).c_str(), "", nbin, ed);
    hTot = new TH1D("mjj_raw_total", "", nbin, ed);
  }
};

}  // namespace

int main(int argc, char** argv) {
  SetAtlasStyle();
  gStyle->SetOptStat(0);
  gErrorIgnoreLevel = kFatal;
  // REQUIRED before any histogram is constructed: every worker thread builds
  // an identically-named ThreadState histogram set, and without this they race
  // on ROOT's global gDirectory registry. Same reason clustering_hist.cxx and
  // rpt_v5_hist.cxx call it.
  TH1::AddDirectory(kFALSE);
  // The PDF page comes out portrait with the landscape canvas in one corner --
  // ROOT's default paper, shared by every plot in this project. SetPaperSize
  // does NOT override it here (tried); pdf_to_jpg.sh trims the border on the
  // way to a slide, so it is cosmetic only. Don't "fix" it again.

  SampleConfig cfg = resolveSample(argc, argv);
  ENERGY_LABEL    = cfg.energyLabel;
  OUTPUT_DIR      = cfg.outputDir;
  SAMPLE_NAME     = cfg.sampleName;
  OVERLAP_REMOVAL = cfg.overlapRemoval;
  const Long64_t maxEvents = resolveMaxEvents(argc, argv);
  const unsigned nThreads  = resolveThreads(argc, argv);

  // Strip the VBF topology selection. classifyVbsRegion applies VBS_JET_MJJ and
  // VBS_JET_D_ETA itself and returns NONE when either fails, so leaving them at
  // their defaults would classify every low-m_jj event as NONE and produce a
  // plot that merely redraws its own cut.
  //
  // Set here unconditionally rather than expecting --vbs-mjj=0 --vbs-deta=0 on
  // the command line: a forgotten flag on a grid submission would not error, it
  // would return a plausible-looking all-NONE plot. resolveSelection() is
  // deliberately NOT called, which also leaves SELECTION_TAG empty so the
  // output names stay unadorned.
  VBS_JET_MJJ   = 0.0;
  VBS_JET_D_ETA = 0.0;

  boost::filesystem::create_directories(OUTPUT_DIR + "/diagnostics");
  if (SAMPLE_NAME.empty())
    boost::filesystem::create_directories(OUTPUT_DIR + "/hists");

  TChain chain("ntuple");
  setupChain(chain, cfg.ntupleDir.c_str());
  if (chain.GetListOfFiles()->GetEntries() == 0) {
    std::cerr << "No ROOT files found.  Aborting.\n";
    return 1;
  }
  ROOT::EnableImplicitMT(nThreads);

  // Per-thread state registry, merged into one after the loop. The mutex
  // serializes the whole ThreadState construction, not just the push_back --
  // see the note on the struct, and the TColor race this pattern was
  // originally hardened against in clustering_dt.cxx.
  std::mutex stateRegistryMutex;
  std::vector<std::unique_ptr<ThreadState>> stateRegistry;
  std::atomic<Long64_t> progressCounter{0};

  // One chain scan for the progress denominator -- chain.GetEntries() opens
  // every file, so it is called exactly once here and never again (see
  // setupChain's own note on why its guard uses GetListOfFiles instead).
  const Long64_t nEntries = chain.GetEntries();
  const Long64_t denom    = (maxEvents > 0 && maxEvents < nEntries) ? maxEvents : nEntries;
  std::cout << "Starting Event Loop (" << nThreads << " threads)\n";

  ROOT::TTreeProcessorMT proc(chain, nThreads);
  proc.Process([&](TTreeReader& reader) {
    // Fresh per invocation: a worker can be handed a different TTreeReader
    // across task ranges, so this cannot be thread_local.
    BranchPointerWrapper branch(reader);

    thread_local ThreadState* tlState = nullptr;
    if (!tlState) {
      std::lock_guard<std::mutex> lock(stateRegistryMutex);
      stateRegistry.push_back(std::make_unique<ThreadState>());
      tlState = stateRegistry.back().get();
    }
    ThreadState& state = *tlState;

  while (reader.Next()) {
    // Shared counter, same pattern as clustering_hist.cxx: stop this worker's
    // range once the global cap is crossed rather than capping per-thread.
    const Long64_t n = ++progressCounter;
    if (maxEvents > 0 && n > maxEvents) break;
    ++state.nSeen;
    if (n % 20000 == 0)
      std::cout << "Progress: " << n << "/" << denom << "\r" << std::flush;

    // Same preamble as processEventData / vbs_mjj_diag, minus the topology cut.
    branch.computeOverlapRemoval();
    if (!branch.passLeptonSelection()) continue;
    if (branch.vetoLeptonOverlap())    continue;
    if (!branch.passBasicCuts())       continue;

    std::vector<int> passPtIdx;
    int nPt = 0, nPtEta = 0;
    branch.collectPtPassingJets(passPtIdx, nPt, nPtEta);
    if (nPt    < MIN_PASSPT_JETS)  continue;
    if (nPtEta < MIN_PASSETA_JETS) continue;

    // No opposite-hemisphere pair means no m_jj, so the event cannot be placed
    // on this axis at all. Excluded from the denominator and counted, so the
    // exclusion shows up in the summary instead of quietly shrinking the sample.
    auto pair = branch.calcBestVbsPair(passPtIdx);
    if (!pair.valid()) { ++state.nNoPair; continue; }
    ++state.nSel;

    // Event-level region from the SHARED classifier (clustering_structs.h),
    // so vbs_region_mjj and rpt_v5_hist cannot drift on what the regions mean.
    // Membership is decided over EVERY jet in the event via isJetTruthHS, not
    // over the VBS candidate pair -- m_jj below is only the x-axis variable.
    int cnt[5] = {0, 0, 0, 0, 0};  // nFwdHS, nFwdPU, nCenHS, nAnyHS, nFwdHSTrk
    const EventRegion region =
        branch.classifyEventRegion(MIN_ABS_ETA_JET, MAX_ABS_ETA_JET,
                                   MIN_ABS_ETA_JET, cnt);

    Cat cat;
    switch (region) {
      case EventRegion::R1:       cat = R1;       break;
      case EventRegion::R2:       cat = R2;       break;
      case EventRegion::CAN_HELP: cat = CAN_HELP; break;
      case EventRegion::NO_T0:    cat = NO_T0;    break;
      default:                    cat = MAY_NOT;  break;
    }

    // Sub-case bookkeeping for the summary, mirroring the ladder in
    // classifyEventRegion so the printed breakdown cannot disagree with it.
    if (region == EventRegion::NO_T0) {
      ++state.nNoT0;
    } else if (region == EventRegion::CAN_HELP) {
      ++state.nHelpNoHS;                          // fake, no HS jet, but HS tracks exist
    } else if (region == EventRegion::MAY_NOT) {
      if      (cnt[0] >= 1) ++state.nHelpFwdHS;   // fwd HS, no competing fwd PU
      else if (cnt[1] >= 1) ++state.nMayNotHsBeyond;  // fwd PU, HS beyond acceptance
      else                  ++state.nMayNotNoFwd;     // nothing forward at all
    }
    if (cnt[0] >= 1 || cnt[1] >= 1) ++state.nAnyFwd;

    // Below the axis floor the event cannot be placed on this plot. Dropped
    // rather than piled into the first column, which would misreport it, and
    // counted so the drop is visible in the summary.    // Below the axis floor the event cannot be placed on this plot. Dropped
    // rather than piled into the first column, which would misreport it, and
    // counted so the drop is visible in the summary.
    if (pair.mjj < MJJ_EDGES.front()) { ++state.nBelowMjj; continue; }

    const double m = folded(state.hTot, pair.mjj);
    state.raw[cat]->Fill(m);
    state.hTot->Fill(m);
    ++state.nCat[cat];
  }
  });
  std::cout << "\nFINISHED PROCESSING\n";

  // --- Merge per-thread state into one ---
  if (stateRegistry.empty()) {
    std::cerr << "No events processed.  Aborting.\n";
    return 1;
  }
  ThreadState& merged = *stateRegistry.front();
  for (size_t i = 1; i < stateRegistry.size(); ++i) {
    ThreadState& o = *stateRegistry[i];
    for (int c = 0; c < NCAT; ++c) {
      merged.raw[c]->Add(o.raw[c]);
      merged.nCat[c] += o.nCat[c];
    }
    merged.hTot->Add(o.hTot);
    merged.nSeen     += o.nSeen;      merged.nSel        += o.nSel;
    merged.nNoPair   += o.nNoPair;    merged.nBelowMjj   += o.nBelowMjj;
    merged.nHelpFwdHS += o.nHelpFwdHS; merged.nHelpNoHS  += o.nHelpNoHS;
    merged.nNoT0 += o.nNoT0;
    merged.nMayNotHsBeyond += o.nMayNotHsBeyond;
    merged.nMayNotNoFwd    += o.nMayNotNoFwd;
    merged.nAnyFwd   += o.nAnyFwd;
  }

  // Local aliases so the plotting/summary code below reads unchanged.
  auto& raw  = merged.raw;
  auto* hTot = merged.hTot;
  auto& nCat = merged.nCat;
  const long nSeen = merged.nSeen, nSel = merged.nSel;
  const long nNoPair = merged.nNoPair, nBelowMjj = merged.nBelowMjj;
  const long nHelpFwdHS = merged.nHelpFwdHS, nHelpNoHS = merged.nHelpNoHS;
  const long nNoT0 = merged.nNoT0;
  const long nMayNotHsBeyond = merged.nMayNotHsBeyond;
  const long nMayNotNoFwd = merged.nMayNotNoFwd;
  const long nAnyFwd = merged.nAnyFwd;
  (void)nSeen;
  const int nbin = (int)MJJ_EDGES.size() - 1;

  // ── Per-column normalisation ─────────────────────────────────────────────
  // Empty columns are left empty rather than filled with NaN: a bin with no
  // events has no composition, and 0/0 would poison the stack's axis range.
  std::array<TH1D*, NCAT> frac{};
  for (int c = 0; c < NCAT; ++c) {
    frac[c] = (TH1D*)raw[c]->Clone((std::string("mjj_frac_") + CAT_KEY[c]).c_str());
    frac[c]->Reset();
  }
  for (int b = 1; b <= nbin; ++b) {
    const double tot = hTot->GetBinContent(b);
    if (tot <= 0) continue;
    for (int c = 0; c < NCAT; ++c) {
      const double n = raw[c]->GetBinContent(b);
      frac[c]->SetBinContent(b, n / tot);
      // Binomial error on the fraction -- the honest uncertainty here, since
      // the six categories are a partition of the same column, not independent.
      frac[c]->SetBinError(b, std::sqrt(std::max(0.0, n * (1.0 - n / tot))) / tot);
    }
  }

  // ── Plots ────────────────────────────────────────────────────────────────
  const std::string pdf = plotFilePath("diagnostics", "vbs_region_mjj.pdf");
  TCanvas* c1 = new TCanvas("c1", "c1", 800, 600);
  c1->Print((pdf + "[").c_str());

  // Page 1: the composition stack.
  c1->Clear();
  c1->SetLogy(false);
  THStack* stk = new THStack("stk_region_mjj", "");
  for (int c = 0; c < NCAT; ++c) {
    frac[c]->SetFillColor(CAT_COLOR[c]);
    frac[c]->SetLineColor(kBlack);
    frac[c]->SetLineWidth(1);
    stk->Add(frac[c]);
  }
  stk->Draw("HIST");
  stk->GetXaxis()->SetTitle("m_{jj} [GeV]   (last bin includes overflow)");
  stk->GetYaxis()->SetTitle("Fraction of events");
  // Every column reaches exactly 1, so the headroom for the labels has to come
  // from the axis maximum -- there is no empty corner inside a full stack.
  stk->SetMaximum(1.9);
  stk->SetMinimum(0.0);
  c1->Modified();

  // Mark the analysis working point. 500 is a bin edge (see MJJ_EDGES), so the
  // line falls between columns rather than cutting one in half.
  TLine* cut = new TLine(500.0, 0.0, 500.0, 1.0);
  cut->SetLineStyle(2);
  cut->SetLineWidth(2);
  cut->SetLineColor(kBlack);
  cut->Draw();

  // ATLAS labels on top, legend beneath them.
  ATLASLabel(0.18, 0.90, "Simulation Internal");
  ATLASEnergyLabel(0.18, 0.855, MyUtl::ENERGY_LABEL.c_str());

  TLegend* leg = new TLegend(0.17, 0.61, 0.93, 0.82);
  StyleLegend(leg);
  leg->SetNColumns(2);
  for (int c = 0; c < NCAT; ++c) leg->AddEntry(frac[c], CAT_LABEL[c], "f");
  leg->AddEntry(cut, "m_{jj} = 500 GeV (analysis cut)", "l");
  leg->Draw();

  c1->Print(pdf.c_str());

  // Page 2: the raw spectrum behind it. Without this the eye reads a clean
  // fraction in a column holding three events as though it meant something.
  c1->Clear();
  c1->SetLogy(true);
  hTot->SetLineColor(kBlack);
  hTot->SetLineWidth(2);
  hTot->GetXaxis()->SetTitle("m_{jj} [GeV]   (last bin includes overflow)");
  hTot->GetYaxis()->SetTitle("Events per bin");
  hTot->SetMinimum(0.5);
  hTot->SetMaximum(400.0 * hTot->GetMaximum());  // headroom for the taller legend
  hTot->Draw("HIST");
  for (int c = 0; c < NCAT; ++c) {
    raw[c]->SetLineColor(CAT_COLOR[c]);
    raw[c]->SetLineWidth(2);
    raw[c]->Draw("HIST SAME");
  }
  ATLASLabel(0.18, 0.90, "Simulation Internal");
  ATLASEnergyLabel(0.18, 0.855, MyUtl::ENERGY_LABEL.c_str());
  TLegend* leg2 = new TLegend(0.17, 0.59, 0.93, 0.82);
  StyleLegend(leg2);
  leg2->SetNColumns(2);
  leg2->AddEntry(hTot, "All (denominator)", "l");
  for (int c = 0; c < NCAT; ++c) leg2->AddEntry(raw[c], CAT_LABEL[c], "l");
  leg2->Draw();
  c1->Print(pdf.c_str());

  c1->Print((pdf + "]").c_str());

  // ── ROOT output ──────────────────────────────────────────────────────────
  // Raw as well as normalised: the fractions alone cannot be re-binned, summed
  // across samples, or given an error bar after the fact.
  const std::string rootPath = histFilePath("vbs_region_mjj.root");
  {
    TFile f(rootPath.c_str(), "RECREATE");
    for (int c = 0; c < NCAT; ++c) { raw[c]->Write(); frac[c]->Write(); }
    hTot->Write();
    f.Close();
  }

  // ── Console summary ──────────────────────────────────────────────────────
  printf("\n================================================================\n");
  printf("  VBS topology vs m_jj: where can HGTD contribute?\n");
  printf("================================================================\n");
  printf("  Events read                      : %8ld\n", nSeen);
  printf("  Passing selection, with a pair   : %8ld  (%.2f%% of read)\n",
         nSel, nSeen ? 100.0 * nSel / nSeen : 0.0);
  printf("    excluded: no opp.-hemisphere pair : %8ld\n", nNoPair);
  // Only meaningful when the axis has a floor above 0; silent otherwise rather
  // than printing a guaranteed zero.
  if (MJJ_EDGES.front() > 0.0)
    printf("    excluded: m_jj below %.0f GeV      : %8ld\n",
           MJJ_EDGES.front(), nBelowMjj);
  printf("  Plotted (denominator of every column): %8ld\n", nSel - nBelowMjj);
  const long nPlot = nSel - nBelowMjj;
  printf("\n  %-34s %10s %9s\n", "category", "events", "of plotted");
  printf("  %s\n", std::string(56, '-').c_str());
  long sum = 0;
  for (int c = 0; c < NCAT; ++c) {
    sum += nCat[c];
    printf("  %-34s %10ld %8.2f%%\n", CAT_LABEL[c], nCat[c],
           nPlot ? 100.0 * nCat[c] / nPlot : 0.0);
  }
  printf("  %s\n", std::string(56, '-').c_str());
  printf("  %-34s %10ld %8.2f%%\n", "TOTAL", sum,
         nPlot ? 100.0 * sum / nPlot : 0.0);
  if (sum != nPlot)
    printf("\n  *** WARNING: categories sum to %ld but %ld events were plotted.\n"
           "      The category ladder is not exhaustive; every column's stack\n"
           "      is therefore normalised against the wrong denominator.\n",
           sum, nPlot);

  // Both "Other" bands broken out, so neither is a black box on the plot.
  printf("\n  \"Other: HGTD can help\" is made of:\n");
  printf("      fwd PU, no HS jet, and    %8ld (%.2f%%)  reject the fake, and a real\n"
         "        >=1 forward HS TRACK                             t0 exists to do it with\n",
         nHelpNoHS, nPlot ? 100.0*nHelpNoHS/nPlot : 0.0);

  // Split out because "can help" was claiming a job the event cannot support:
  // with no forward hard-scatter track the only forward tracks belong to the
  // pileup jet, so the forward t0 IS that pileup vertex's time.
  printf("\n  \"fwd PU only: t0 set by the fake\":\n");
  printf("      fwd PU, no HS jet, and    %8ld (%.2f%%)  no forward hard-scatter\n"
         "        ZERO forward HS tracks                           track to build a t0 from\n",
         nNoT0, nPlot ? 100.0*nNoT0/nPlot : 0.0);

  printf("\n  \"Other: HGTD may not help\" is made of:\n");
  printf("      fwd HS, no fwd PU         %8ld (%.2f%%)  genuine tag, but no competing\n"
         "                                                   fake -- nothing to disambiguate\n",
         nHelpFwdHS, nPlot ? 100.0*nHelpFwdHS/nPlot : 0.0);
  printf("      no forward jet at all     %8ld (%.2f%%)  nothing timeable\n",
         nMayNotNoFwd, nPlot ? 100.0*nMayNotNoFwd/nPlot : 0.0);
  printf("      fwd PU, HS only beyond |eta| %.1f %5ld (%.2f%%)  edge case, kept visible\n",
         MAX_ABS_ETA_JET, nMayNotHsBeyond,
         nPlot ? 100.0*nMayNotHsBeyond/nPlot : 0.0);

  // The headline this plot exists to give.
  const long nReach = nCat[R1] + nCat[R2] + nCat[CAN_HELP];
  printf("\n  HGTD can contribute: %ld  (%.2f%% of plotted)\n",
         nReach, nPlot ? 100.0 * nReach / nPlot : 0.0);

  // A different denominator: not "of all selected events" but "of events with
  // any forward jet at all" -- i.e. conditioned on HGTD having something in
  // acceptance to measure in the first place. R1, R2 and CAN_HELP are each by
  // construction a subset of this (every matching condition requires a forward
  // jet), so their shares say what fraction of the reachable-in-principle
  // population each region actually accounts for.
  printf("\n  Of events with >=1 forward jet (N=%ld, %.2f%% of plotted):\n",
         nAnyFwd, nPlot ? 100.0 * nAnyFwd / nPlot : 0.0);
  printf("    1. %-38s %8ld (%.2f%%)\n", CAT_LABEL[R1],
         nCat[R1], nAnyFwd ? 100.0 * nCat[R1] / nAnyFwd : 0.0);
  printf("    2. %-38s %8ld (%.2f%%)\n", CAT_LABEL[R2],
         nCat[R2], nAnyFwd ? 100.0 * nCat[R2] / nAnyFwd : 0.0);
  printf("    3. %-38s %8ld (%.2f%%)\n", CAT_LABEL[CAN_HELP],
         nCat[CAN_HELP], nAnyFwd ? 100.0 * nCat[CAN_HELP] / nAnyFwd : 0.0);
  printf("       %-38s %8ld (%.2f%%)\n", "(forward jet but unreachable)",
         nAnyFwd - nCat[R1] - nCat[R2] - nCat[CAN_HELP],
         nAnyFwd ? 100.0 * (nAnyFwd - nCat[R1] - nCat[R2] - nCat[CAN_HELP]) / nAnyFwd : 0.0);

  printf("\n  Wrote %s\n", pdf.c_str());
  printf("  Wrote %s\n", rootPath.c_str());
  return 0;
}
