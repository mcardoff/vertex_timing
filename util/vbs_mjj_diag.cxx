// ---------------------------------------------------------------------------
// vbs_mjj_diag.cxx — m_jj distribution of the VBS candidate pair.
//
// The selection currently defines its "VBS topology" with |Deta| >= 3 on the
// candidate pair (VBS_JET_D_ETA, runtime-settable via --vbs-deta=). m_jj is the
// standard VBS discriminant, and calcBestVbsPair already computes it while
// finding the pair. Before swapping the cut, this dumps what the distribution
// actually looks like.
//
// Denominator: every event the rest of the selection admits -- lepton
// selection / overlap veto / passBasicCuts / >=MIN_PASSPT_JETS above
// MIN_JET_PT / >=MIN_PASSETA_JETS forward -- with an opposite-hemisphere pair,
// but WITHOUT any |Deta| or m_jj requirement. That is the population a cut on
// either variable gets chosen from.
//
// Every distribution is split by whether both legs of the candidate pair are
// truth-HS-matched (isJetTruthHS). A pair with a pileup leg means the event
// entered the topology because of pileup rather than because of its physics --
// the pathology documented on VBS_JET_D_ETA in clustering_constants.h -- so the
// split is what says whether m_jj is a cleaner topology definition than |Deta|
// or merely a different one.
//
// Pure diagnostic: reads only, changes no analysis behaviour.
//
// Output: <OUTPUT_DIR>/diagnostics/vbs_mjj_diag.pdf  (3 pages) + console tables
// ---------------------------------------------------------------------------

#include <TCanvas.h>
#include <TChain.h>
#include <TH1.h>
#include <TH2.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TTreeReader.h>

#include <boost/filesystem.hpp>
#include <cstdio>
#include <iostream>
#include <vector>

#include "clustering_constants.h"
#include "event_processing.h"
#include "plotting_utilities.h"
#include "sample_config.h"
#include "AtlasStyle.h"
#include "AtlasLabels.h"

using namespace MyUtl;

namespace {

  // m_jj axis: 0-4000 GeV in 40 GeV bins, overflow folded into the last bin so
  // the multi-TeV VBF tail stays visible without stretching the axis.
  const int    MJJ_NBIN = 100;
  const double MJJ_MIN  = 0.0, MJJ_MAX = 4000.0;
  const int    DETA_NBIN = 90;
  const double DETA_MIN  = 0.0, DETA_MAX = 9.0;

  void fillFolded(TH1D* h, double v) {
    const double last = h->GetXaxis()->GetBinCenter(h->GetNbinsX());
    h->Fill(std::min(v, last));
  }

  // One (threshold -> retained) scan line for either variable.
  struct ScanRow {
    double thr;
    long   nAll, nHS;   // events kept; those whose candidate pair is both-legs-HS
  };

  std::vector<ScanRow> scan(const std::vector<double>& thresholds,
                            const std::vector<double>& vals,
                            const std::vector<char>&   bothHS) {
    std::vector<ScanRow> rows;
    rows.reserve(thresholds.size());
    for (double t : thresholds) {
      ScanRow r{t, 0, 0};
      for (size_t i = 0; i < vals.size(); ++i) {
        if (vals[i] < t) continue;
        ++r.nAll;
        if (bothHS[i]) ++r.nHS;
      }
      rows.push_back(r);
    }
    return rows;
  }

  void printScan(const char* title, const char* unit,
                 const std::vector<ScanRow>& rows, long nTot, long nTotHS) {
    printf("\n  %s\n", title);
    printf("    %-10s %10s %8s %10s %8s %8s\n",
           "cut", "kept", "of all", "genuine", "purity", "HS eff");
    printf("    %s\n", std::string(58, '-').c_str());
    for (const auto& r : rows) {
      printf("    >= %-7.4g %10ld %7.2f%% %10ld %7.2f%% %7.2f%%\n",
             r.thr, r.nAll,
             nTot   ? 100.0 * r.nAll / nTot   : 0.0,
             r.nHS,
             r.nAll ? 100.0 * r.nHS  / r.nAll : 0.0,
             nTotHS ? 100.0 * r.nHS  / nTotHS : 0.0);
    }
    (void)unit;
  }

  void drawTriple(TCanvas* c, TH1D* hAll, TH1D* hHS, TH1D* hPU,
                  const char* xtitle, bool logy) {
    c->Clear();
    c->SetLogy(logy);
    hAll->SetLineColor(kBlack);  hAll->SetLineWidth(2);
    hHS ->SetLineColor(C01);     hHS ->SetLineWidth(2);
    hPU ->SetLineColor(C02);     hPU ->SetLineWidth(2);
    hPU ->SetLineStyle(2);

    hAll->GetXaxis()->SetTitle(xtitle);
    hAll->GetYaxis()->SetTitle("Events");
    hAll->SetMaximum((logy ? 20.0 : 1.5) * hAll->GetMaximum());
    hAll->SetMinimum(logy ? 0.5 : 0.0);
    hAll->Draw("HIST");
    hHS ->Draw("HIST SAME");
    hPU ->Draw("HIST SAME");

    TLegend* leg = new TLegend(0.55, 0.72, 0.92, 0.88);
    StyleLegend(leg);
    leg->AddEntry(hAll, "All candidate pairs",   "l");
    leg->AddEntry(hHS,  "Both legs truth-HS",    "l");
    leg->AddEntry(hPU,  "#geq1 pileup leg",      "l");
    leg->Draw();

    ATLASLabel(0.18, 0.88, "Simulation Internal");
    ATLASEnergyLabel(0.18, 0.82, MyUtl::ENERGY_LABEL.c_str());
    c->Print(plotFilePath("diagnostics", "vbs_mjj_diag.pdf").c_str());
  }

}  // namespace

// Axis label for |Deta_jj|. Written out as a difference of the two etas rather
// than with #Delta: uppercase #Delta renders as a BLANK glyph under ROOT 6.40 on
// macOS (lowercase #delta, #Sigma and #eta are all fine), which silently eats
// the symbol. Same reason m_jj's fold is spelled out below instead of annotated.
const char* DETA_TITLE = "|#eta_{j1} - #eta_{j2}|";

int main(int argc, char** argv) {
  SetAtlasStyle();
  gStyle->SetOptStat(0);
  gErrorIgnoreLevel = kFatal;

  SampleConfig cfg = resolveSample(argc, argv);
  ENERGY_LABEL = cfg.energyLabel;
  OUTPUT_DIR   = cfg.outputDir;
  SAMPLE_NAME  = cfg.sampleName;
  OVERLAP_REMOVAL = cfg.overlapRemoval;
  const Long64_t maxEvents = resolveMaxEvents(argc, argv);

  boost::filesystem::create_directories(OUTPUT_DIR + "/diagnostics");

  TChain chain("ntuple");
  setupChain(chain, cfg.ntupleDir.c_str());
  TTreeReader reader(&chain);
  BranchPointerWrapper branch(reader);

  TH1D* hMjjAll  = new TH1D("mjj_all",  "", MJJ_NBIN,  MJJ_MIN,  MJJ_MAX);
  TH1D* hMjjHS   = new TH1D("mjj_hs",   "", MJJ_NBIN,  MJJ_MIN,  MJJ_MAX);
  TH1D* hMjjPU   = new TH1D("mjj_pu",   "", MJJ_NBIN,  MJJ_MIN,  MJJ_MAX);
  TH1D* hDetaAll = new TH1D("deta_all", "", DETA_NBIN, DETA_MIN, DETA_MAX);
  TH1D* hDetaHS  = new TH1D("deta_hs",  "", DETA_NBIN, DETA_MIN, DETA_MAX);
  TH1D* hDetaPU  = new TH1D("deta_pu",  "", DETA_NBIN, DETA_MIN, DETA_MAX);
  TH2D* h2       = new TH2D("mjj_vs_deta",
                            (std::string(";") + DETA_TITLE + ";m_{jj} [GeV]").c_str(),
                            45, DETA_MIN, DETA_MAX, 50, MJJ_MIN, MJJ_MAX);

  // Per-event values retained for the threshold scans (cheap: one double + one
  // char per selected event, and the scans need the full list to sweep cuts).
  std::vector<double> mjjVals, detaVals;
  std::vector<char>   bothHSVals;

  long nSeen = 0, nSel = 0, nNoPair = 0, nBothHS = 0, nFwdLeg = 0, nFwdLegHS = 0;

  const Long64_t nEntries = chain.GetEntries();
  const Long64_t denom    = (maxEvents > 0 && maxEvents < nEntries) ? maxEvents : nEntries;

  while (reader.Next()) {
    if (maxEvents > 0 && nSeen >= maxEvents) break;
    ++nSeen;
    if (nSeen % 20000 == 0)
      std::cout << "Progress: " << nSeen << "/" << denom << "\r" << std::flush;

    // Same preamble as processEventData, minus the |Deta| requirement.
    branch.computeOverlapRemoval();
    if (!branch.passLeptonSelection()) continue;
    if (branch.vetoLeptonOverlap())    continue;
    if (!branch.passBasicCuts())       continue;

    std::vector<int> passPtIdx;
    int nPt = 0, nPtEta = 0;
    branch.collectPtPassingJets(passPtIdx, nPt, nPtEta);
    if (nPt    < MIN_PASSPT_JETS)  continue;
    if (nPtEta < MIN_PASSETA_JETS) continue;

    auto pair = branch.calcBestVbsPair(passPtIdx);
    if (!pair.valid()) { ++nNoPair; continue; }  // no opposite-hemisphere pair
    ++nSel;

    const bool bothHS = branch.isJetTruthHS(pair.idxI) &&
                        branch.isJetTruthHS(pair.idxJ);
    if (bothHS) ++nBothHS;

    // Is either leg of the pair inside HGTD acceptance? That is the leg whose
    // timing the detector can actually contribute to.
    auto isFwd = [&](int j) {
      double e = std::abs(branch.topoJetEta[j]);
      return e > MIN_ABS_ETA_JET && e < MAX_ABS_ETA_JET;
    };
    const bool fwdLeg = isFwd(pair.idxI) || isFwd(pair.idxJ);
    if (fwdLeg) { ++nFwdLeg; if (bothHS) ++nFwdLegHS; }

    fillFolded(hMjjAll, pair.mjj);
    fillFolded(bothHS ? hMjjHS : hMjjPU, pair.mjj);
    hDetaAll->Fill(pair.dEta);
    (bothHS ? hDetaHS : hDetaPU)->Fill(pair.dEta);
    h2->Fill(pair.dEta, std::min(pair.mjj, MJJ_MAX - 1.0));

    mjjVals.push_back(pair.mjj);
    detaVals.push_back(pair.dEta);
    bothHSVals.push_back(bothHS ? 1 : 0);
  }
  std::cout << "\nFINISHED PROCESSING\n";

  // ── Plots ─────────────────────────────────────────────────────────────────
  const std::string pdf = plotFilePath("diagnostics", "vbs_mjj_diag.pdf");
  TCanvas* c = new TCanvas("c", "c", 800, 600);
  c->Print((pdf + "[").c_str());
  drawTriple(c, hMjjAll,  hMjjHS,  hMjjPU,
             "m_{jj} [GeV]   (last bin = overflow)", true);
  drawTriple(c, hDetaAll, hDetaHS, hDetaPU, DETA_TITLE, false);
  c->Clear();
  c->SetLogy(false);
  c->SetLogz(true);
  c->SetRightMargin(0.16);
  h2->Draw("COLZ");
  ATLASLabel(0.18, 0.88, "Simulation Internal");
  ATLASEnergyLabel(0.18, 0.82);
  c->Print(pdf.c_str());
  c->Print((pdf + "]").c_str());

  // ── Console summary ───────────────────────────────────────────────────────
  printf("\n================================================================\n");
  printf("  VBS candidate pair: m_jj vs |Deta| as the topology cut\n");
  printf("================================================================\n");
  printf("  Events read                    : %8ld\n", nSeen);
  printf("  Passing selection minus |Deta| : %8ld  (%.2f%% of read)\n",
         nSel, nSeen ? 100.0 * nSel / nSeen : 0.0);
  printf("    ... rejected: no opp.-hemisphere pair : %8ld\n", nNoPair);
  printf("    ... both pair legs truth-HS   : %8ld  (%.2f%%)\n",
         nBothHS, nSel ? 100.0 * nBothHS / nSel : 0.0);
  printf("    ... >=1 pair leg in HGTD accept.: %6ld  (%.2f%%), of which both-HS %ld (%.2f%%)\n",
         nFwdLeg, nSel ? 100.0 * nFwdLeg / nSel : 0.0,
         nFwdLegHS, nFwdLeg ? 100.0 * nFwdLegHS / nFwdLeg : 0.0);
  printf("\n  \"genuine\" = both legs of the candidate pair truth-HS-matched.\n");
  printf("  purity = genuine/kept.  HS eff = genuine kept / all genuine.\n");

  const std::vector<double> mjjThr  = {0, 100, 200, 300, 400, 500, 600, 800,
                                       1000, 1250, 1500, 2000};
  const std::vector<double> detaThr = {0, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 5.0};
  auto mjjRows  = scan(mjjThr,  mjjVals,  bothHSVals);
  auto detaRows = scan(detaThr, detaVals, bothHSVals);
  // |Deta| working point this diagnostic compares against. Deliberately a
  // LOCAL constant and NOT VBS_JET_D_ETA_DEFAULT: this tool exists to compare
  // candidate working points against each other, so its reference has to stay
  // put when the analysis default moves. It moved on 2026-08-26 (3.0 -> 0.0,
  // partly because of THIS diagnostic's output); tracking it would have made
  // the comparison below "|Deta| >= 0", i.e. no cut at all, silently turning
  // the equal-efficiency test into a comparison of the baseline with itself.
  constexpr double DETA_REFERENCE = 3.0;

  printScan("m_jj cut scan [GeV]", "GeV", mjjRows,  nSel, nBothHS);
  printScan("|Deta| cut scan",     "",    detaRows, nSel, nBothHS);

  // Equal-efficiency comparison: find the m_jj threshold that keeps the same
  // fraction of genuine pairs as the |Deta| >= DETA_REFERENCE
  // working point, then compare purities. This is the apples-to-apples question
  // -- at matched signal efficiency, which cut admits less pileup topology?
  long detaKeptHS = 0, detaKept = 0;
  for (size_t i = 0; i < detaVals.size(); ++i) {
    if (detaVals[i] < DETA_REFERENCE) continue;
    ++detaKept;
    if (bothHSVals[i]) ++detaKeptHS;
  }
  if (detaKeptHS > 0 && nBothHS > 0) {
    const double targetEff = (double)detaKeptHS / nBothHS;
    double bestThr = 0.0, bestKeptHS = nBothHS, bestKept = nSel;
    for (double t = 0.0; t <= 3000.0; t += 10.0) {
      long kHS = 0, k = 0;
      for (size_t i = 0; i < mjjVals.size(); ++i) {
        if (mjjVals[i] < t) continue;
        ++k;
        if (bothHSVals[i]) ++kHS;
      }
      if ((double)kHS / nBothHS < targetEff) break;
      bestThr = t; bestKeptHS = kHS; bestKept = k;
    }
    printf("\n  Equal-genuine-efficiency comparison\n");
    printf("    |Deta| >= %.1f  : eff %.2f%%  purity %.2f%%  (kept %ld)\n",
           DETA_REFERENCE, 100.0 * targetEff,
           detaKept ? 100.0 * detaKeptHS / detaKept : 0.0, detaKept);
    printf("    m_jj  >= %-5.0f : eff %.2f%%  purity %.2f%%  (kept %ld)\n",
           bestThr, nBothHS ? 100.0 * bestKeptHS / nBothHS : 0.0,
           bestKept ? 100.0 * bestKeptHS / bestKept : 0.0, (long)bestKept);
  }
  printf("\n  Wrote %s\n", pdf.c_str());
  printf("================================================================\n");
  return 0;
}
