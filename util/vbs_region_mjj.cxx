// ---------------------------------------------------------------------------
// vbs_region_mjj.cxx — VBS-topology composition of the selection vs m_jj.
//
// The analysis argues that HGTD's gain is concentrated in two VBS topologies
// (see classifyVbsRegion in clustering_structs.h):
//
//   R1  both candidate legs forward, one paper-HS + one paper-PU
//   R2  forward paper-PU leg + central paper-HS leg
//
// but nothing so far says what FRACTION of the selected sample actually lands
// in them, or how that fraction depends on m_jj -- which is the knob the
// selection uses to buy topology purity. This answers both: a stacked
// composition plot, one column per m_jj bin, every column normalised to 1.
//
// Reading it: the y axis is "given an event at this m_jj, what topology is it?"
// Because each column is normalised independently, the falling m_jj spectrum
// cannot wash out the tail, where R1 matters most. The raw (un-normalised)
// spectrum is drawn on a second page so low-statistics columns -- whose
// fractions are noisy however clean they look -- stay identifiable.
//
// Denominator: every event the rest of the selection admits -- lepton selection
// / overlap veto / passBasicCuts / >=MIN_PASSPT_JETS above MIN_JET_PT /
// >=MIN_PASSETA_JETS forward -- that HAS an opposite-hemisphere candidate pair,
// but WITHOUT any m_jj or |Deta| requirement. That is the whole point: the plot
// shows what a cut on m_jj would be selecting FROM, so it must not itself apply
// one.
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
  // EXHAUSTIVE by construction: per-column normalisation is only meaningful if
  // every plotted event lands in exactly one of these, so the ladder below ends
  // in a catch-all rather than a final `if`.
  //
  // AMBIG is not a padding category. isJetPaperHS and isJetPaperPU
  // (clustering_structs.h) are NOT complements -- HS wants a match within
  // dR 0.3 of a truth HS jet above 10 GeV, PU wants NO match within dR 0.6
  // above 4 GeV -- so a jet landing between those two cones is neither. Folding
  // those events into "both PU" would misreport them; dropping them would make
  // the stack quietly sum to less than 1.
  enum Cat { R1 = 0, R2, MIXED_ETA, BOTH_HS, BOTH_PU, AMBIG, NCAT };

  const std::array<const char*, NCAT> CAT_KEY = {
    "r1", "r2", "mixed_eta", "both_hs", "both_pu", "ambig"
  };

  const std::array<const char*, NCAT> CAT_LABEL = {
    "R1: both tags fwd, HS + PU",
    "R2: fwd PU + central HS",
    "HS + PU, other #eta config",
    "Both tags truth-HS",
    "Both tags pileup",
    "Ambiguous truth match"
  };

  const std::array<Color_t, NCAT> CAT_COLOR = {
    C01,  // R1     blue
    C07,  // R2     orange
    C03,  // mixed  yellow
    C08,  // bothHS green
    C09,  // bothPU ash
    C06   // ambig  brown
  };

  // ── m_jj binning ─────────────────────────────────────────────────────────
  // Variable width: narrow where the spectrum is steep and every sample has
  // events, widening into the tail so Z+jets columns are not single-event
  // fractions. 500 is deliberately an edge -- it is the analysis working point,
  // so the reader can see the cut position without a drawn line.
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

  long nSeen = 0, nSel = 0, nNoPair = 0;
  std::array<long, NCAT> nCat{};

  // Breakdown of MIXED_ETA (one paper-HS leg + one paper-PU leg, but not R1/R2)
  // by where each leg actually sits: 0 = central |eta| < MIN_ABS_ETA_JET,
  // 1 = forward MIN_ABS_ETA_JET..MAX_ABS_ETA_JET, 2 = outside both (|eta| beyond
  // MAX_ABS_ETA_JET, or exactly on a boundary). Indexed [hs leg][pu leg], so the
  // asymmetry between "forward HS + central PU" and R2's "forward PU + central
  // HS" is visible rather than averaged away.
  long nMix[3][3] = {};

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

    // R1/R2 from the SHARED classifier, so the regions here cannot mean
    // something different from the ones rpt_v5 and the clustering analysis fill.
    const VbsRegion region =
        branch.classifyVbsRegion(MIN_ABS_ETA_JET, MAX_ABS_ETA_JET, MIN_ABS_ETA_JET);

    Cat cat;
    if (region == VbsRegion::R1) {
      cat = R1;
    } else if (region == VbsRegion::R2) {
      cat = R2;
    } else {
      // Same helpers classifyVbsRegion uses, on the same pair, so the remainder
      // is split on exactly the truth definitions that defined R1/R2.
      const int a = pair.idxI, b = pair.idxJ;
      const bool hsA = branch.isJetPaperHS(branch.topoJetEta[a], branch.topoJetPhi[a]);
      const bool puA = branch.isJetPaperPU(branch.topoJetEta[a], branch.topoJetPhi[a]);
      const bool hsB = branch.isJetPaperHS(branch.topoJetEta[b], branch.topoJetPhi[b]);
      const bool puB = branch.isJetPaperPU(branch.topoJetEta[b], branch.topoJetPhi[b]);

      auto zone = [&](int j) {
        const double e = std::abs((double)branch.topoJetEta[j]);
        if (e < MIN_ABS_ETA_JET)                              return 0;  // central
        if (e > MIN_ABS_ETA_JET && e < MAX_ABS_ETA_JET)       return 1;  // forward
        return 2;                                                        // neither
      };

      if ((hsA && puB) || (hsB && puA)) {
        cat = MIXED_ETA;
        const int hsIdx = (hsA && puB) ? a : b;
        const int puIdx = (hsA && puB) ? b : a;
        ++nMix[zone(hsIdx)][zone(puIdx)];
      }
      else if (hsA && hsB) cat = BOTH_HS;
      else if (puA && puB) cat = BOTH_PU;
      else                 cat = AMBIG;
    }

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
  // Every column reaches exactly 1, so the headroom for the legend and the
  // ATLAS labels has to come from the axis maximum -- there is no empty corner
  // inside a stack that is full by construction.
  stk->SetMaximum(1.8);
  stk->SetMinimum(0.0);
  c1->Modified();

  // Mark the analysis working point. 500 is a bin edge (see MJJ_EDGES), so the
  // line falls between columns rather than cutting one in half.
  TLine* cut = new TLine(500.0, 0.0, 500.0, 1.0);
  cut->SetLineStyle(2);
  cut->SetLineWidth(2);
  cut->SetLineColor(kBlack);
  cut->Draw();
  // The line goes in the legend rather than getting a floating TLatex label:
  // every column reaches 1, so there is no gap over the stack to write into
  // without colliding with the ATLAS labels.
  TLegend* leg = new TLegend(0.17, 0.72, 0.93, 0.92);
  StyleLegend(leg);
  leg->SetNColumns(2);
  for (int c = 0; c < NCAT; ++c) leg->AddEntry(frac[c], CAT_LABEL[c], "f");
  leg->AddEntry(cut, "m_{jj} = 500 GeV (analysis cut)", "l");
  leg->Draw();

  ATLASLabel(0.18, 0.66, "Simulation Internal");
  ATLASEnergyLabel(0.18, 0.61, MyUtl::ENERGY_LABEL.c_str());
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
  hTot->SetMaximum(20.0 * hTot->GetMaximum());
  hTot->Draw("HIST");
  for (int c = 0; c < NCAT; ++c) {
    raw[c]->SetLineColor(CAT_COLOR[c]);
    raw[c]->SetLineWidth(2);
    raw[c]->Draw("HIST SAME");
  }
  TLegend* leg2 = new TLegend(0.17, 0.74, 0.93, 0.91);
  StyleLegend(leg2);
  leg2->SetNColumns(2);
  leg2->AddEntry(hTot, "All (denominator)", "l");
  for (int c = 0; c < NCAT; ++c) leg2->AddEntry(raw[c], CAT_LABEL[c], "l");
  leg2->Draw();
  ATLASLabel(0.18, 0.68, "Simulation Internal");
  ATLASEnergyLabel(0.18, 0.63, MyUtl::ENERGY_LABEL.c_str());
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
  printf("  VBS-topology composition vs m_jj  (no m_jj / |Deta| cut applied)\n");
  printf("================================================================\n");
  printf("  Events read                      : %8ld\n", nSeen);
  printf("  Passing selection, with a pair   : %8ld  (%.2f%% of read)\n",
         nSel, nSeen ? 100.0 * nSel / nSeen : 0.0);
  printf("    excluded: no opp.-hemisphere pair : %8ld\n", nNoPair);
  printf("\n  %-34s %10s %9s\n", "category", "events", "of sel.");
  printf("  %s\n", std::string(56, '-').c_str());
  long sum = 0;
  for (int c = 0; c < NCAT; ++c) {
    sum += nCat[c];
    printf("  %-34s %10ld %8.2f%%\n", CAT_LABEL[c], nCat[c],
           nSel ? 100.0 * nCat[c] / nSel : 0.0);
  }
  printf("  %s\n", std::string(56, '-').c_str());
  printf("  %-34s %10ld %8.2f%%\n", "TOTAL", sum,
         nSel ? 100.0 * sum / nSel : 0.0);
  if (sum != nSel)
    printf("\n  *** WARNING: categories sum to %ld but %ld events were selected.\n"
           "      The category ladder is not exhaustive; every column's stack\n"
           "      is therefore normalised against the wrong denominator.\n",
           sum, nSel);

  // Where the "HS + PU, other eta config" events actually sit. This is the
  // largest non-both-HS category, so it is worth saying what it is rather than
  // leaving it as a residual.
  const char* ZONE[3] = { "central", "forward", "outside" };
  printf("\n  \"HS + PU, other eta config\" breakdown (rows = HS leg, cols = PU leg):\n");
  printf("    %-10s %10s %10s %10s %10s\n", "HS \\ PU", ZONE[0], ZONE[1], ZONE[2], "total");
  printf("    %s\n", std::string(54, '-').c_str());
  for (int i = 0; i < 3; ++i) {
    long rowsum = nMix[i][0] + nMix[i][1] + nMix[i][2];
    printf("    %-10s %10ld %10ld %10ld %10ld\n",
           ZONE[i], nMix[i][0], nMix[i][1], nMix[i][2], rowsum);
  }
  printf("    (forward/forward is R1 and forward-PU/central-HS is R2, so both are\n"
         "     absent here by construction: central = |eta| < %.2f, forward =\n"
         "     %.2f-%.2f, outside = beyond that.)\n",
         MIN_ABS_ETA_JET, MIN_ABS_ETA_JET, MAX_ABS_ETA_JET);

  printf("\n  Wrote %s\n", pdf.c_str());
  printf("  Wrote %s\n", rootPath.c_str());
  return 0;
}
