// ---------------------------------------------------------------------------
// vbs_region_mjj.cxx — VBS-topology composition of the selection vs m_jj.
//
// The analysis argues that HGTD's gain is concentrated in the VBS topologies,
// but nothing so far said what FRACTION of the selection it can act on at all.
// This answers that: a stacked composition plot, one column per m_jj bin, every
// column normalised to 1, with the two reachable topologies at the bottom.
//
//   R1        both tags forward, one HS + one PU  -- timing says WHICH is real
//   R2        forward PU + central HS             -- timing rejects the fake
//   BOTH_HS   both VBS jets truth hard-scatter    -- nothing to reject
//   OTHER     everything else
//
// Reading it: the y axis is "given an event at this m_jj, what topology is it?"
// Because each column is normalised independently, the falling m_jj spectrum
// cannot wash out the tail. The raw spectrum is on page 2 so low-statistics
// columns -- whose fractions are noisy however clean they look -- stay
// identifiable.
//
// Denominator: every event the rest of the selection admits -- lepton selection
// / overlap veto / passBasicCuts / >=MIN_PASSPT_JETS above MIN_JET_PT /
// >=MIN_PASSETA_JETS forward -- that HAS an opposite-hemisphere candidate pair
// with |Delta eta| >= 2.5. The m_jj axis starts at the 500 GeV working point
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

#include <boost/filesystem.hpp>
#include <array>
#include <cstdio>
#include <iostream>
#include <string>
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
  // Four bands, per the simplified scheme:
  //
  //   R1       both tags forward, one HS + one PU -- timing says WHICH is real
  //   R2       forward PU + central HS            -- timing rejects the fake
  //   BOTH_HS  both VBS jets truth hard-scatter (any eta) -- a genuine VBS
  //            pair; nothing for timing to reject
  //   OTHER    everything else: mixed pairs in other eta configurations
  //            (incl. the old R3), both-PU pairs, out-of-acceptance pairs, and
  //            pairs with a neither-label jet (isJetPaperHS / isJetPaperPU are
  //            NOT complements, so a jet between the two cones is neither)
  //
  // EXHAUSTIVE by construction: per-column normalisation is meaningless
  // otherwise, so the ladder ends in a catch-all and main() checks the totals.
  enum Cat { R1 = 0, R2, BOTH_HS, OTHER, NCAT };

  const std::array<const char*, NCAT> CAT_KEY = {
    "r1", "r2", "both_hs", "other"
  };

  const std::array<const char*, NCAT> CAT_LABEL = {
    "R1: both tags fwd, HS + PU",
    "R2: fwd PU + central HS",
    "Both VBS jets truth HS",
    "Other"
  };

  // The two regions timing can act on are saturated; the rest muted.
  const std::array<Color_t, NCAT> CAT_COLOR = {
    C01,        // R1       blue  -+ HGTD can contribute
    C02,        // R2       red   -+
    kGray,      // both HS  (dominant, lightest neutral)
    kGray + 2   // other    dark grey
  };

  // ── m_jj binning ─────────────────────────────────────────────────────────
  // Starts at the 500 GeV analysis working point: below it the sample is not
  // the one the analysis runs on, and including it compressed everything above
  // into half the axis. Events below 500 are dropped from the plot and counted
  // separately in the summary, so the exclusion is stated rather than implied.
  //
  // Variable width: narrow where the spectrum is steep and every sample has
  // events, widening into the tail so Z+jets columns are not built from a
  // handful of events.
  const std::vector<double> MJJ_EDGES = {
    500, 750, 1000, 1250, 1500, 1750, 2000, 2500, 3000, 4000
  };

  // VBS candidate-pair |Delta eta| requirement for THIS plot. Applied as an
  // explicit counted exclusion in the loop (like the m_jj axis floor), NOT via
  // the VBS_JET_D_ETA global -- classifyVbsRegion consults that global
  // internally and returns NONE on failure, which would silently break the
  // local-vs-shared R1/R2 cross-check. The pair itself is still
  // calcBestVbsPair's highest-m_jj opposite-hemisphere pair; the cut is on the
  // chosen pair, it does not re-choose one that passes.
  const double PLOT_MIN_DETA = 2.5;

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

}  // namespace

int main(int argc, char** argv) {
  SetAtlasStyle();
  gStyle->SetOptStat(0);
  gErrorIgnoreLevel = kFatal;
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
  TTreeReader reader(&chain);
  BranchPointerWrapper branch(reader);

  const int    nbin  = (int)MJJ_EDGES.size() - 1;
  const double* eddy = MJJ_EDGES.data();

  std::array<TH1D*, NCAT> raw{};
  for (int c = 0; c < NCAT; ++c)
    raw[c] = new TH1D((std::string("mjj_raw_") + CAT_KEY[c]).c_str(), "", nbin, eddy);

  TH1D* hTot = new TH1D("mjj_raw_total", "", nbin, eddy);

  long nSeen = 0, nSel = 0, nNoPair = 0, nBelowMjj = 0, nBelowDeta = 0;
  long nDisagreeR1 = 0, nDisagreeR2 = 0;
  std::array<long, NCAT> nCat{};

  // Where OTHER's mixed-truth events actually sit, so the residual does not
  // become the black box that "other eta config" was before R3 was split out of
  // it. Indexed [HS leg zone][PU leg zone]: 0 = central |eta| < MIN_ABS_ETA_JET,
  // 1 = HGTD acceptance, 2 = beyond MAX_ABS_ETA_JET.

  const Long64_t nEntries = chain.GetEntries();
  const Long64_t denom    = (maxEvents > 0 && maxEvents < nEntries) ? maxEvents : nEntries;

  while (reader.Next()) {
    if (maxEvents > 0 && nSeen >= maxEvents) break;
    ++nSeen;
    if (nSeen % 20000 == 0)
      std::cout << "Progress: " << nSeen << "/" << denom << "\r" << std::flush;

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
    if (!pair.valid()) { ++nNoPair; continue; }
    ++nSel;

    // R1/R2 from the SHARED classifier, at the same unbounded-forward window
    // the clustering analysis now uses (VBS_FWD_ETA_MAX). rpt_v5 stays on its
    // own tighter window on purpose -- see the constant.
    const VbsRegion region =
        branch.classifyVbsRegion(MIN_ABS_ETA_JET, VBS_FWD_ETA_MAX, MIN_ABS_ETA_JET);

    const int a = pair.idxI, b = pair.idxJ;

    // Zones match classifyVbsRegion's own isFwd/isCentral exactly, including
    // the strict inequalities: a jet sitting precisely on MIN_ABS_ETA_JET is
    // neither central nor in acceptance there, so it must be neither here.
    auto zone = [&](int j) {
      const double e = std::abs((double)branch.topoJetEta[j]);
      if (e < MIN_ABS_ETA_JET)                          return 0;  // central
      if (e > MIN_ABS_ETA_JET && e < VBS_FWD_ETA_MAX)   return 1;  // forward
      return 2;                                                    // beyond
    };

    // Same truth helpers classifyVbsRegion uses, on the same pair.
    const bool hsA = branch.isJetPaperHS(branch.topoJetEta[a], branch.topoJetPhi[a]);
    const bool puA = branch.isJetPaperPU(branch.topoJetEta[a], branch.topoJetPhi[a]);
    const bool hsB = branch.isJetPaperHS(branch.topoJetEta[b], branch.topoJetPhi[b]);
    const bool puB = branch.isJetPaperPU(branch.topoJetEta[b], branch.topoJetPhi[b]);
    const int  zA  = zone(a), zB = zone(b);

    // ── The two regions, defined together ───────────────────────────────────
    //
    // isR1/isR2 must still agree with classifyVbsRegion, since that is what
    // rpt_v5 and the clustering analysis fill -- cross-checked per event below,
    // so a future edit to either definition cannot drift silently.
    const bool isR1 = (zA == 1 && zB == 1) && ((hsA && puB) || (hsB && puA));
    const bool isR2 = (zA == 1 && puA && zB == 0 && hsB) ||
                      (zB == 1 && puB && zA == 0 && hsA);
    if (isR1 != (region == VbsRegion::R1)) ++nDisagreeR1;
    if (isR2 != (region == VbsRegion::R2)) ++nDisagreeR2;

    Cat cat;
    if      (isR1)       cat = R1;
    else if (isR2)       cat = R2;
    // Truth-content band regardless of eta: a genuine VBS pair is a genuine
    // VBS pair wherever its legs sit, and the question the plot answers is
    // "what is this event", not "is it in acceptance".
    else if (hsA && hsB) cat = BOTH_HS;
    else                 cat = OTHER;

    // Below the axis floor the event cannot be placed on this plot. Dropped
    // rather than piled into the first column, which would misreport it, and
    // counted so the drop is visible in the summary.
    if (pair.dEta < PLOT_MIN_DETA)    { ++nBelowDeta; continue; }
    if (pair.mjj < MJJ_EDGES.front()) { ++nBelowMjj; continue; }

    const double m = folded(hTot, pair.mjj);
    raw[cat]->Fill(m);
    hTot->Fill(m);
    ++nCat[cat];
  }
  std::cout << "\nFINISHED PROCESSING\n";

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

  // ATLAS labels on top, legend beneath them. The 500 GeV cut line that used to
  // sit here is gone: the axis now starts at the cut, so it would just redraw
  // the frame edge.
  ATLASLabel(0.18, 0.90, "Simulation Internal");
  ATLASEnergyLabel(0.18, 0.855, MyUtl::ENERGY_LABEL.c_str());

  TLegend* leg = new TLegend(0.17, 0.61, 0.93, 0.82);
  StyleLegend(leg);
  leg->SetNColumns(2);
  for (int c = 0; c < NCAT; ++c) leg->AddEntry(frac[c], CAT_LABEL[c], "f");
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
  printf("    excluded: |dEta| below %.1f          : %8ld\n",
         PLOT_MIN_DETA, nBelowDeta);
  printf("    excluded: m_jj below %.0f GeV        : %8ld\n",
         MJJ_EDGES.front(), nBelowMjj);
  printf("  Plotted (denominator of every column): %8ld\n",
         nSel - nBelowDeta - nBelowMjj);
  const long nPlot = nSel - nBelowDeta - nBelowMjj;
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

  // The headline this plot exists to give.
  const long nReach = nCat[R1] + nCat[R2];
  printf("\n  HGTD can contribute (R1 + R2): %ld  (%.2f%% of plotted)\n",
         nReach, nPlot ? 100.0 * nReach / nPlot : 0.0);

  if (nDisagreeR1 || nDisagreeR2)
    printf("\n  *** WARNING: local R1/R2 disagree with classifyVbsRegion on\n"
           "      %ld / %ld events. The regions here no longer mean what\n"
           "      rpt_v5 and the clustering analysis fill.\n",
           nDisagreeR1, nDisagreeR2);
  else
    printf("  (local R1/R2 agree with classifyVbsRegion on every event)\n");

  printf("\n  Wrote %s\n", pdf.c_str());
  printf("  Wrote %s\n", rootPath.c_str());
  return 0;
}
