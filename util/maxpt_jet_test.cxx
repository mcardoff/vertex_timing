// ---------------------------------------------------------------------------
// maxpt_jet_test.cxx
//
// Quick A/B/C test: today every track entering clustering must satisfy a HARD
// upper cut pT < MAX_TRACK_PT (30 GeV).  This drops genuine high-pT HS tracks.
// Here we compare the WAVeS / TRKPTZ efficiency under three maxTrkPt choices:
//
//   Baseline    : MAX_TRACK_PT          (30 GeV, current)
//   Leading jet : leading reco-jet pT   (per event)
//   Subleading  : 2nd-highest reco-jet pT (per event)
//
// Reports pass rates, net change, and rescued (fail->pass) / regressed
// (pass->fail) counts vs baseline.  Pure diagnostic — changes nothing.
// ---------------------------------------------------------------------------

#include <TROOT.h>
#include "clustering_constants.h"
#include "event_processing.h"

using namespace MyUtl;

// Cluster `tracks` (real times, iterative) and report whether the best cluster
// by each score passes the efficiency window.
static void evalCap(const std::vector<int>& tracks, BranchPointerWrapper* branch,
                    bool& wavesPass, bool& trkptzPass) {
  wavesPass = trkptzPass = false;
  if (tracks.empty()) return;
  auto cl = clusterTracksInTime(
    tracks, branch, DIST_CUT_CONE,
    /*useSmearedTimes=*/false, /*checkTimeValid=*/true, /*smearRes=*/-1.0,
    ClusteringMethod::ITERATIVE, /*useZ0=*/false,
    /*sortTracks=*/false, /*calcPurityFlag=*/false);
  if (cl.empty()) return;

  const Cluster *bw = &cl[0], *bt = &cl[0];
  double mw = bw->scores.at(Score::WAVES.id), mt = bt->scores.at(Score::TRKPTZ.id);
  for (const auto& c : cl) {
    double sw = c.scores.at(Score::WAVES.id);
    double st = c.scores.at(Score::TRKPTZ.id);
    if (sw > mw) { mw = sw; bw = &c; }
    if (st > mt) { mt = st; bt = &c; }
  }
  wavesPass  = bw->passEfficiency(branch);
  trkptzPass = bt->passEfficiency(branch);
}

// Pass counts + rescued/regressed vs baseline, for one (cap, score) combo.
struct Tally { long pass = 0, resc = 0, regr = 0; };
static void acc(Tally& t, bool base, bool var) {
  t.pass += var;
  if (!base && var) t.resc++;
  if (base && !var) t.regr++;
}

auto main() -> int {
  TChain chain("ntuple");
  setupChain(chain, "../../ntuple-hgtd/");
  TTreeReader reader(&chain);
  BranchPointerWrapper branch(reader);
  ROOT::EnableImplicitMT();
  gErrorIgnoreLevel = kFatal;

  long total = 0;
  long wBaseP = 0, tBaseP = 0;            // baseline pass counts
  Tally wLead, wSub, tLead, tSub;         // variant tallies
  long subTighter = 0, subLooser = 0;
  double sumLead = 0, sumSub = 0;

  const Long64_t N = chain.GetEntries();
  while (reader.Next()) {
    const Long64_t r = chain.GetReadEntry() + 1;
    if (r % 100 == 0) std::cout << "Progress: " << r << "/" << N << "\r" << std::flush;

    if (!branch.passBasicCuts()) continue;
    if (!branch.passJetPtCut())  continue;

    // Top-two reco-jet pT.
    double j1 = 0., j2 = 0.;
    for (int j = 0; j < (int)branch.topoJetPt.GetSize(); ++j) {
      double p = branch.topoJetPt[j];
      if (p > j1)      { j2 = j1; j1 = p; }
      else if (p > j2) { j2 = p; }
    }
    const double leadCap = j1 > 0. ? j1 : MAX_TRACK_PT;
    const double subCap  = j2 > 0. ? j2 : MAX_TRACK_PT;  // <2 jets → no change
    sumLead += leadCap; sumSub += subCap;
    if (subCap < MAX_TRACK_PT) subTighter++;
    else if (subCap > MAX_TRACK_PT) subLooser++;

    auto trkBase = getAssociatedTracks(&branch, MIN_TRACK_PT, MAX_TRACK_PT, 3.0);
    auto trkLead = getAssociatedTracks(&branch, MIN_TRACK_PT, leadCap,      3.0);
    auto trkSub  = getAssociatedTracks(&branch, MIN_TRACK_PT, subCap,       3.0);

    bool wB, tB, wL, tL, wS, tS;
    evalCap(trkBase, &branch, wB, tB);
    evalCap(trkLead, &branch, wL, tL);
    evalCap(trkSub,  &branch, wS, tS);

    total++;
    wBaseP += wB; tBaseP += tB;
    acc(wLead, wB, wL); acc(wSub, wB, wS);
    acc(tLead, tB, tL); acc(tSub, tB, tS);
  }
  std::cout << "\nFINISHED PROCESSING\n";

  auto pct = [&](long n) { return total > 0 ? 100.*n/total : 0.; };
  printf("\n================================================================\n");
  printf("  MAX_TRACK_PT  ->  reco-jet pT   (efficiency A/B/C test)\n");
  printf("================================================================\n");
  printf("  Total events:             %8ld   (baseline cap %.0f GeV)\n", total, MAX_TRACK_PT);
  printf("  Mean leading-jet cap:     %8.1f GeV\n", total > 0 ? sumLead/total : 0.);
  printf("  Mean subleading-jet cap:  %8.1f GeV\n", total > 0 ? sumSub/total  : 0.);
  printf("  Subleading cap < %.0f GeV:  %8ld   (tighter than baseline)\n", MAX_TRACK_PT, subTighter);
  printf("  Subleading cap > %.0f GeV:  %8ld   (looser  than baseline)\n", MAX_TRACK_PT, subLooser);

  auto report = [&](const char* name, long baseP, const Tally& lead, const Tally& sub) {
    printf("\n  --- %s  (baseline %ld pass, %.3f%%) ---\n", name, baseP, pct(baseP));
    auto line = [&](const char* tag, const Tally& t) {
      printf("    %-14s %8ld pass  (%.3f%%)  net %+5ld (%+.3f%%)   resc %4ld  regr %4ld\n",
             tag, t.pass, pct(t.pass), t.pass - baseP, pct(t.pass) - pct(baseP), t.resc, t.regr);
    };
    line("leading jet",    lead);
    line("subleading jet", sub);
  };
  report("WAVeS",  wBaseP, wLead, wSub);
  report("TRKPTZ", tBaseP, tLead, tSub);
  printf("================================================================\n");
  return 0;
}
