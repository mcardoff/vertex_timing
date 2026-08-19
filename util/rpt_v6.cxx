// ---------------------------------------------------------------------------
// rpt_v6.cxx — DIRECT PORT of util/myJet_ana_fr.C (a colleague's RpT study).
//
// This is NOT a new version of rpt_v5 and deliberately shares no code with it.
// The point is to reproduce the reference analysis as literally as possible on
// our ntuples, so that when its output disagrees with rpt_v5's the difference
// can be attributed to the analysis rather than to a porting shortcut. Every
// cut value, formula, histogram binning, ROC construction and canvas below is
// the original's; only the branch reads are remapped.
//
// ── Branch remapping (their EventTree -> our "ntuple") ──────────────────────
//   jet_pt/eta/phi        -> AntiKt4EMTopoJets_pt/eta/phi
//   jet_tracks_idx        -> AntiKt4EMTopoJets_ghostTrack_idx
//   track_z0              -> Track_z0
//   track_var_z0          -> Track_var_z0
//   track_theta/phi       -> Track_theta / Track_phi
//   track_qOverP          -> Track_qOverP
//   track_t30             -> Track_time          (see NOTE 2)
//   recovertex_z          -> RecoVtx_z
//   truthvertex_z/t       -> TruthVtx_z / TruthVtx_time
//   jet_isHS / jet_isPU   -> NOT AVAILABLE       (see NOTE 1)
//
// ── Substitutions the port could not avoid ─────────────────────────────────
// NOTE 1 (jet labels). Their ntuple carries per-jet jet_isHS/jet_isPU flags;
//   ours does not. Substituted with the ATL-HGTD-PUB-2022-001 Sec. 3 cone
//   definitions used throughout this codebase: HS = dR<0.3 from a truth HS jet
//   with pT>10 GeV, PU = dR>0.6 from every truth HS jet with pT>4 GeV. Their
//   "(isHS==1 && isPU==0)" / "(isHS==0 && isPU==1)" tests map onto
//   isPaperHS && !isPaperPU / !isPaperHS && isPaperPU respectively, preserving
//   the original's exclusivity. Jets satisfying neither are dropped, as there.
// NOTE 2 (track time). track_t30 is NOT a reconstructed time. Slide 2 of
//   jet_pileup_studies_5th_October_2023.pdf states it outright:
//     "Track time is retrieved from truth level information / Assign time to
//      tracks from truth vertex info with a smearing: track -> particle_truth
//      -> vertex_truth / Different smearings (gaussian) are considered:
//      30, 60, 90ps."
//   So t30 = Gaus(truth VERTEX time of the track's truth vertex, 30 ps), and
//   the literal /30.0 in the pull is simply dividing by that same smearing
//   width. Mapped here onto Track_truthVtx_idx -> TruthVtx_time (96.6% of
//   tracks carry the link; the rest cannot be assigned a time and are skipped
//   for the timed sums only).
//   Our Track_time (the real HGTD reconstructed time) is NOT used and must not
//   be: an earlier revision of this port mapped t30 onto it, which silently fed
//   the ~59% of tracks with Track_hasValidTime == 0 -- all carrying
//   Track_time == 0.0 exactly -- into the cut, and made the timed ROC come out
//   WORSE than the untimed one. The whole reference study is a truth-smearing
//   exercise; there is no reconstructed track time anywhere in it.
// NOTE 3 (track pT). They recompute pT from qOverP/theta; we read Track_pt
//   (already GeV). Same quantity, read instead of derived. Their second form
//   pt2 = |1e-3/qOverP|*sin(theta) is likewise Track_pt, and is what the RpT
//   numerator accumulates.
//
// ── Things that differ from rpt_v5 and are LIKELY the source of disagreement ─
//   * NEITHER time in this study is reconstructed. Track times are truth vertex
//     times smeared by SMEAR_PS (30/60/90), and the "reco" vertex time is the
//     truth vertex time smeared by 10 ps. rpt_v5 uses real HGTD track times
//     (with their real ~41% validity and non-Gaussian tails) and real
//     HGTD/WAVeS/TRKPTZ vertex times. The reference therefore measures an
//     IDEALISED CEILING, not a competing reconstruction -- which is the single
//     biggest reason its numbers sit above rpt_v5's.
//   * Vertex-quality cut is |z_reco - z_truth| < 0.5 mm, not MAX_VTX_DZ = 2.0.
//   * Jet window is 30 < pT < 50 GeV (a single slice), not 30-40 / >40.
//   * Track pT window is 1 < pT < 45 GeV, not MIN/MAX_TRACK_PT = 1/30.
//   * Track-to-vertex association is the getNewDzpara parameterization at
//     |z|/dz0 < 1.4, and the jet cone is dR < 0.2.
//   Single-threaded on purpose: the original is sequential, and the per-track
//   TRandom draw makes the result order-dependent.
//
// Output: <OUTPUT_DIR>/rpt_v6/*.pdf  (the same five the macro writes)
// ---------------------------------------------------------------------------

#include <TCanvas.h>
#include <TChain.h>
#include <TF1.h>
#include <TGraph.h>
#include <TH1.h>
#include <TH2.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TLine.h>
#include <TPad.h>
#include <TRandom.h>
#include <TStyle.h>
#include <TTreeReader.h>
#include <TTreeReaderArray.h>
#include <TVector2.h>

#include <boost/filesystem.hpp>
#include <cmath>
#include <cstdio>
#include <iostream>
#include <algorithm>
#include <sstream>
#include <string>
#include <vector>

#include "clustering_constants.h"
#include "sample_config.h"
#include "event_processing.h"   // setupChain only

// ---------------------------------------------------------------------------
// getNewDzpara — verbatim from myJet_ana_fr.C, including the non-else-if
// chain (a value landing exactly on a bin boundary takes the LAST matching
// branch; preserved deliberately).
// ---------------------------------------------------------------------------
// Track cutflow diagnostic. rpt_v6 runs single-threaded (no TTreeProcessorMT),
// so plain counters are safe. Counts tracks surviving each stage of the per-jet
// selection so a track deficit can be pinned on a specific cut rather than
// inferred from the RpT distribution.
struct TrkCutflow {
  long long jets = 0, ghost = 0, idx_ok = 0, pt_eta = 0, zcut = 0, dr = 0;
};
// Split by jet truth-type: the aggregate is dominated by pile-up jets, whose
// tracks legitimately fail a PV0 association, so a low overall pass rate says
// nothing on its own. The HS row is the one that carries information.
static TrkCutflow CF_HS, CF_PU;

static double getNewDzpara(double ETA, double PT) {
  ETA = std::fabs(ETA);
  double p_v[7] = {0, 0, 0, 0, 0, 0, 0};
  auto put = [&p_v](const double (&s)[7]) { std::copy(s, s + 7, p_v); };

  if (PT <= 1.5) {
    static const double t[7] = { 0.0314036, 0.790955, -2.65987, 3.62073, -2.18228, 0.614866, -0.0634521 };
    put(t);
  }
  if ((PT > 1.5) && (PT <= 2.5)) {
    static const double t[7] = { 0.0229273, 0.540101, -1.80727, 2.45187, -1.47382, 0.414345, -0.0426769 };
    put(t);
  }
  if ((PT > 2.5) && (PT <= 5.0)) {
    static const double t[7] = { 0.0163773, 0.345112, -1.14474, 1.54382, -0.923523, 0.258617, -0.0265446 };
    put(t);
  }
  if ((PT > 5) && (PT <= 10)) {
    static const double t[7] = { 0.010919, 0.179329, -0.581971, 0.773186, -0.45679, 0.126608, -0.012875 };
    put(t);
  }
  if (PT > 10) {
    static const double t[7] = { 0.00835945, 0.0957783, -0.299255, 0.38722, -0.22351, 0.0607521, -0.00606524 };
    put(t);
  }

  double Dzpara = p_v[0] + p_v[1]*ETA + p_v[2]*pow(ETA,2) + p_v[3]*pow(ETA,3) + p_v[4]*pow(ETA,4);
  Dzpara += p_v[5]*pow(ETA,5) + p_v[6]*pow(ETA,6);
  return std::fabs(Dzpara);
}

static void myText(Double_t x, Double_t y, Color_t color, const char* text) {
  TLatex l;
  l.SetNDC();
  l.SetTextColor(color);
  l.SetTextSize(0.04);
  l.DrawLatex(x, y, text);
}

static void ATLAS_LABEL(Double_t x, Double_t y, Color_t /*color*/) {
  TLatex* tex = new TLatex(x, y, "ATLAS Simulation Internal");
  tex->SetNDC();
  tex->SetTextFont(72);
  tex->SetLineWidth(2);
  tex->Draw();
}

static void SetCanvasAttr(TCanvas* tc) {
  tc->Range(-28, -70.10724, 149.1636, 56.43432);
  tc->SetFillColor(0);
  tc->SetBorderMode(0);
  tc->SetBorderSize(2);
  tc->SetTickx(1);
  tc->SetTicky(1);
  tc->SetLeftMargin(0.158046);
  tc->SetRightMargin(0.05172414);
  tc->SetTopMargin(0.05084746);
  tc->SetBottomMargin(0.1588983);
  tc->SetFrameBorderMode(0);
}

static inline double dR_(double e1, double p1, double e2, double p2) {
  double de = e1 - e2;
  double dp = TVector2::Phi_mpi_pi(p1 - p2);
  return std::sqrt(de*de + dp*dp);
}

int main(int argc, char** argv) {
  gStyle->SetOptStat(0);

  auto sample = MyUtl::resolveSample(argc, argv);
  MyUtl::ENERGY_LABEL = sample.energyLabel;
  MyUtl::OUTPUT_DIR   = sample.outputDir;
  MyUtl::SAMPLE_NAME  = sample.sampleName;
  const Long64_t maxEvents = MyUtl::resolveMaxEvents(argc, argv);

  // Track-time smearing width in ps. The deck studies 30/60/90 (slide 2); the
  // macro's active path is t30, with its t60/t90 variants sitting commented out
  // beside it. The same value is used both to smear and to divide, exactly as
  // the macro does (track_t30 with /30.0, track_t60 with /60.0, ...).
  // Time-pull threshold at which the self-tagging clustering splits a jet.
  const double SELFTAG_NSIGMA = 3.0;
  // Jet window / acceptance. Defaults are the macro's; the TDR figure needs
  // --ptmin=50 --ptmax=1e9 for its second panel and --etamax=4.0 for both.
  double JPT_MIN = 30.0, JPT_MAX = 50.0, JETA_MIN = 2.4, JETA_MAX = 3.8;
  // Global-t0 window. ATLAS-TDR-031 Sec. 3.3.3: "RpT can be recomputed after
  // removing tracks outside a 2*sigma_t window around the reconstructed time of
  // the hard-scatter vertex." Two sigma, not the reference macro's three.
  // Overridable via --t0sigma=<n>: our reconstructed vertex time cannot
  // actually support 2 sigma (the HS-track pull against it has a core width of
  // 1.76, so a nominal 2 sigma keeps only ~58% of genuine HS tracks), and being
  // able to scan this separates "the window is wrong for our t0" from "the
  // method does not work".
  double T0_NSIGMA = 2.0;
  // Vertex-time error inflation, the timing analogue of Z0_VAR_INFLATION.
  // Measured on this sample: sigma_trk ~ 27.0 ps and sigma_vtx ~ 9.1 ps give a
  // quoted uncertainty of 28.6 ps, but the CORE spread of (t_trk - t_vtx) is
  // 50.8 ps -- the quoted error understates the truth by 1.78x, so a nominal
  // 2 sigma behaves like ~1.1 sigma. Since sigma_trk matches the TDR's
  // "Initial" scenario exactly, the excess sits in the vertex term: making the
  // core pull unit-width needs sigma_vtx_eff = sqrt(50.8^2 - 27.0^2) = 43 ps,
  // i.e. a factor ~4.7 on the quoted 9.1 ps.
  double VTXT_INFL = 1.0;
  // Dominant levers on the ITk-only baseline: a tight z0 window or a tight
  // track-jet cone strips PU tracks and inflates rejection. The macro's values
  // are 1.4 and 0.2; the jet radius is 0.4.
  double DZ0_SCALE = 1.4;
  double TRK_DR    = 0.2;
  // Track-to-vertex association mode.
  //   para   (default) : |z0 - PVz| / getNewDzpara(eta,pT) < DZ0_SCALE, the
  //                      reference macro's active path.
  //   signif           : |z0 - PVz| / sqrt(var_z0) < ZSIG, the standard
  //                      significance cut. The macro computes this quantity
  //                      (signi_cut) and the slide deck defines it, but its
  //                      use is commented out there; --zcut=signif switches to
  //                      it, which is the honest way to test whether the tight
  //                      parameterized window is what inflates the ITk-only
  //                      baseline.
  //   atlas            : |(z0 - PVz) * sin(theta)| < ZSINTH mm, the ACTUAL
  //                      ATLAS association. Published in arXiv:1510.03823
  //                      (the JVT/RpT paper, footnote 2): "Tracks are assigned
  //                      to vertices by requiring |dz x sin(theta)| < 1 mm. In
  //                      cases where more than one vertex satisfies this
  //                      criterion, ambiguity is resolved by choosing the
  //                      vertex with the largest summed pT^2 of tracks."
  //                      The sin(theta) factor is the crux: it makes the
  //                      window in z very LOOSE forward (10 mm at |eta| = 3,
  //                      27 mm at 4) precisely where the other two forms are
  //                      tightest (~0.4 mm), which is why they suppress
  //                      pileup contamination the TDR expects to be present.
  //                      The reference macro computes this quantity
  //                      (delz_cut) but leaves its use commented out.
  // Jet collection. The reference macro reads generic jet_pt/eta/phi, so which
  // collection it used is not recoverable from the source -- and it matters for
  // R_pT, whose DENOMINATOR is the reco jet pT. PFlow improves central jet pT
  // resolution substantially while being near-equivalent forward, which is the
  // shape of our discrepancy (forward reproduces the TDR, central does not).
  std::string JETCOLL = "topo";   // topo | pflow
  std::string ZCUT = "para";
  double ZSIG = 3.0;
  double ZSINTH = 1.0;   // mm
  // Event-level vertex-quality cut |z_reco - z_truth|. The macro uses 0.5 mm;
  // MAX_VTX_DZ elsewhere in this codebase (and the usual ATLAS value) is 2.0.
  // A tight value keeps only events whose PV0 position is essentially exact,
  // which flatters any z-based association.
  double MAX_DZ = 0.5;   // mm
  int    NBINS  = 100;   // R_pT histogram bins (macro value); --nbins to refine
  double SMEAR_PS = 30.0;
  for (int i = 1; i < argc; ++i) {
    std::string a(argv[i]);
    if (a.rfind("--smear=", 0) == 0) SMEAR_PS = std::stod(a.substr(8));
    if (a.rfind("--ptmin=", 0) == 0) JPT_MIN  = std::stod(a.substr(8));
    if (a.rfind("--ptmax=", 0) == 0) JPT_MAX  = std::stod(a.substr(8));
    if (a.rfind("--etamax=",0) == 0) JETA_MAX = std::stod(a.substr(9));
    if (a.rfind("--etamin=",0) == 0) JETA_MIN = std::stod(a.substr(9));
    if (a.rfind("--jets=",0)   == 0) JETCOLL  = a.substr(7);
    if (a.rfind("--t0sigma=",0)== 0) T0_NSIGMA= std::stod(a.substr(10));
    if (a.rfind("--vtxtinfl=",0)==0) VTXT_INFL= std::stod(a.substr(11));
    if (a.rfind("--dz0scale=",0)==0) DZ0_SCALE= std::stod(a.substr(11));
    if (a.rfind("--trkdr=",0)  == 0) TRK_DR   = std::stod(a.substr(8));
    if (a.rfind("--zcut=",0)   == 0) ZCUT     = a.substr(7);
    if (a.rfind("--zsig=",0)   == 0) ZSIG     = std::stod(a.substr(7));
    if (a.rfind("--zsinth=",0) == 0) ZSINTH   = std::stod(a.substr(9));
    if (a.rfind("--maxdz=",0)  == 0) MAX_DZ   = std::stod(a.substr(8));
    if (a.rfind("--nbins=",0)  == 0) NBINS    = std::stoi(a.substr(8));
  }
  std::cout << "[rpt_v6] track-time smearing: " << SMEAR_PS << " ps"
            << "  (t30 = Gaus(truth vertex time, " << SMEAR_PS << "))\n";

  const std::string outDir = MyUtl::OUTPUT_DIR + "/rpt_v6";
  boost::filesystem::create_directories(outDir);
  std::string pre = outDir + "/" +
      (MyUtl::SAMPLE_NAME.empty() ? std::string("") : MyUtl::SAMPLE_NAME + "_");
  {   // keep each configuration's output distinct
    std::ostringstream os;
    if (SMEAR_PS != 30.0) os << "t" << (int)SMEAR_PS << "_";
    if (JPT_MIN != 30.0 || JPT_MAX != 50.0)
      os << "pt" << (int)JPT_MIN << (JPT_MAX > 1e8 ? "plus" : "") << "_";
    if (JETA_MAX != 3.8) os << "eta" << (int)(JETA_MAX*10) << "_";
    if (JETA_MIN != 2.4) os << "etamin" << (int)(JETA_MIN*10) << "_";
    if (JETCOLL != "topo") os << JETCOLL << "_";
    if (ZCUT == "signif") os << "zsig" << (int)ZSIG << "_";
    if (ZCUT == "atlas") {
      os << "zatlas";
      if (ZSINTH != 1.0) { std::ostringstream z; z << ZSINTH; std::string t=z.str();
                           std::replace(t.begin(), t.end(), '.', 'p'); os << t; }
      os << "_";
    }
    if (NBINS != 100) os << "nb" << NBINS << "_";
    if (VTXT_INFL != 1.0) { std::ostringstream v; v << "vi" << VTXT_INFL;
      std::string t=v.str(); std::replace(t.begin(),t.end(),'.','p'); os << t << "_"; }
    if (MAX_DZ != 0.5) { std::ostringstream d; d << "dz" << MAX_DZ; std::string t=d.str();
                         std::replace(t.begin(), t.end(), '.', 'p'); os << t << "_"; }
    pre += os.str();
  }

  TChain chain("ntuple");
  MyUtl::setupChain(chain, sample.ntupleDir.c_str());
  TTreeReader reader(&chain);

  const std::string JC = (JETCOLL == "pflow") ? "AntiKt4EMPFlowJets" : "AntiKt4EMTopoJets";
  TTreeReaderArray<float> jet_pt   (reader, (JC + "_pt").c_str());
  TTreeReaderArray<float> jet_eta  (reader, (JC + "_eta").c_str());
  TTreeReaderArray<float> jet_phi  (reader, (JC + "_phi").c_str());
  TTreeReaderArray<std::vector<int>> jet_tracks_idx(reader, (JC + "_ghostTrack_idx").c_str());
  TTreeReaderArray<float> track_pt (reader, "Track_pt");
  TTreeReaderArray<float> track_eta(reader, "Track_eta");
  TTreeReaderArray<float> track_phi(reader, "Track_phi");
  TTreeReaderArray<float> track_z0 (reader, "Track_z0");
  TTreeReaderArray<float> track_var_z0(reader, "Track_var_z0");
  TTreeReaderArray<int>   track_truthVtx_idx(reader, "Track_truthVtx_idx");
  // Real HGTD quantities, for the "realistic" scenario (see REALISTIC below).
  TTreeReaderArray<float> track_time    (reader, "Track_time");
  TTreeReaderArray<float> track_timeRes (reader, "Track_timeRes");
  TTreeReaderArray<int>   track_timeValid(reader, "Track_hasValidTime");
  TTreeReaderArray<float> recoVtx_time   (reader, "RecoVtx_time");
  TTreeReaderArray<float> recoVtx_timeRes(reader, "RecoVtx_timeRes");
  TTreeReaderArray<int>   recoVtx_valid  (reader, "RecoVtx_hasValidTime");
  TTreeReaderArray<float> recovertex_z (reader, "RecoVtx_z");
  TTreeReaderArray<float> truthvertex_z(reader, "TruthVtx_z");
  TTreeReaderArray<float> truthvertex_t(reader, "TruthVtx_time");
  TTreeReaderArray<float> truthHSJet_pt (reader, "TruthHSJet_pt");
  TTreeReaderArray<float> truthHSJet_eta(reader, "TruthHSJet_eta");
  TTreeReaderArray<float> truthHSJet_phi(reader, "TruthHSJet_phi");

  // R_pT histogram binning. The macro uses 100 bins over [-0.02, 2.0], which
  // samples the ROC too coarsely to resolve the high-efficiency end -- the
  // 0.90-0.95 region comes out as a few long straight segments rather than the
  // plateau-then-fall the TDR shows. The ROC construction reads GetNbinsX(), so
  // raising this refines the curve without touching any logic.
  const int nbins_histo = NBINS;
  TH1F* h_Rpt          = new TH1F("h_Rpt",         "", nbins_histo, -0.02, 2.0);
  TH1F* h_Rpt_pu       = new TH1F("h_Rpt_pu",      "", nbins_histo, -0.02, 2.0);
  TH1F* h_Rpt_t        = new TH1F("h_Rpt_t",       "", nbins_histo, -0.02, 2.0);
  TH1F* h_Rpt_pu_t     = new TH1F("h_Rpt_pu_t",    "", nbins_histo, -0.02, 2.0);
  TH1F* h_Rpt_t_truth  = new TH1F("h_Rpt_t_truth", "", nbins_histo, -0.02, 2.0);
  TH1F* h_Rpt_pu_t_truth = new TH1F("h_Rpt_pu_t_truth", "", nbins_histo, -0.02, 2.0);
  TH1F* h_Rpt_reco     = new TH1F("h_Rpt_reco",    "", nbins_histo, -0.02, 2.0);
  TH1F* h_Rpt_pu_reco  = new TH1F("h_Rpt_pu_reco", "", nbins_histo, -0.02, 2.0);
  TH1F* h_Rpt_self     = new TH1F("h_Rpt_self",    "", nbins_histo, -0.02, 2.0);
  TH1F* h_Rpt_pu_self  = new TH1F("h_Rpt_pu_self", "", nbins_histo, -0.02, 2.0);
  TH1F* h_Rpt_full     = new TH1F("h_Rpt_full",    "", nbins_histo, -0.02, 2.0);
  TH1F* h_Rpt_pu_full  = new TH1F("h_Rpt_pu_full", "", nbins_histo, -0.02, 2.0);
  TH1F* h_Rpt_perft0    = new TH1F("h_Rpt_perft0",    "", nbins_histo, -0.02, 2.0);
  TH1F* h_Rpt_pu_perft0 = new TH1F("h_Rpt_pu_perft0", "", nbins_histo, -0.02, 2.0);

  int num_DZ = 0;
  long n_jets_filled_hs = 0, n_jets_filled_pu = 0;
  // Realistic-scenario bookkeeping: how often the HGTD actually supplies a
  // usable time, i.e. how often the timing cut does any work at all.
  long n_trk_timed = 0, n_trk_zonly_trk = 0, n_trk_zonly_vtx = 0;
  long n_evt_vtx_timed = 0;
  // Self-tagging applicability, the TDR's stated limitation made measurable.
  long n_selftag_applied = 0, n_selftag_toofew = 0, n_selftag_split = 0;
  // How often the COMBINED method actually differs from t0-only: it does so
  // exactly on the jets whose event has no vertex time, since elsewhere the two
  // are identical by construction.
  long n_jets_total = 0, n_jets_fallback = 0, n_jets_fallback_active = 0;
  // Stage-by-stage RpT decomposition: where do the tracks actually go?
  //   0 = pT/eta only, 1 = + z-association, 2 = + track-jet cone
  double sHS0=0,sHS1=0,sHS2=0, sPU0=0,sPU1=0,sPU2=0;
  long   cHS=0, cPU=0, nghost_hs=0, nghost_pu=0, njet_hs=0, njet_pu=0;

  TRandom* r = new TRandom();

  Long64_t seen = 0;
  const Long64_t nEntries = chain.GetEntries();
  std::cout << "Number of events: " << nEntries << std::endl;

  // ######|| Event selections ||######
  while (reader.Next()) {
    if (maxEvents > 0 && seen >= maxEvents) break;
    ++seen;
    if (seen % 20000 == 0) std::cout << "Event #" << seen << "\r" << std::flush;

    if (recovertex_z.GetSize() < 1 || truthvertex_z.GetSize() < 1) continue;
    float DZ = std::fabs(truthvertex_z[0] - recovertex_z[0]);
    if (DZ > MAX_DZ) continue;

    num_DZ = num_DZ + 1;

    // REALISTIC scenario inputs: the actual HGTD vertex time. When the vertex
    // has no valid time the timing cut cannot run at all, and every track in
    // the event falls back to z-information only (kept, never discarded).
    const bool vtx_time_ok = (recoVtx_valid.GetSize() > 0 && recoVtx_valid[0] == 1);
    if (vtx_time_ok) ++n_evt_vtx_timed;
    const double vtx_t   = vtx_time_ok ? recoVtx_time[0]    : 0.0;
    const double vtx_res = vtx_time_ok ? recoVtx_timeRes[0] * VTXT_INFL : 0.0;

    // Paper HS/PU labels (NOTE 1), evaluated per jet below.
    auto isPaperHS = [&](double je, double jp) {
      for (int t = 0; t < (int)truthHSJet_pt.GetSize(); ++t) {
        if (truthHSJet_pt[t] < 10.0) continue;
        if (dR_(je, jp, truthHSJet_eta[t], truthHSJet_phi[t]) < 0.3) return true;
      }
      return false;
    };
    auto isPaperPU = [&](double je, double jp) {
      for (int t = 0; t < (int)truthHSJet_pt.GetSize(); ++t) {
        if (truthHSJet_pt[t] < 4.0) continue;
        if (dR_(je, jp, truthHSJet_eta[t], truthHSJet_phi[t]) < 0.6) return false;
      }
      return true;
    };

    // ######|| Jets selections ||######
    for (int i = 0; i < (int)jet_tracks_idx.GetSize(); ++i) {
      if (jet_pt[i] < JPT_MIN || jet_pt[i] > JPT_MAX) continue;
      if ((std::fabs(jet_eta[i]) < JETA_MIN) || (std::fabs(jet_eta[i]) > JETA_MAX)) continue;

      const bool cf_isHS = isPaperHS(jet_eta[i], jet_phi[i]);
      const bool cf_isPU = isPaperPU(jet_eta[i], jet_phi[i]);
      TrkCutflow* CF = cf_isHS ? &CF_HS : (cf_isPU ? &CF_PU : nullptr);
      if (CF) { CF->jets++; CF->ghost += (long long)jet_tracks_idx[i].size(); }

      float trackPT_0 = 0, trackPT_1 = 0, trackPT_2 = 0, trackPT_3 = 0, trackPT_4 = 0;
      float trackPT_5 = 0;   // REALISTIC / t0-only: real HGTD vertex + track times
      float trackPT_6 = 0;   // SELF-TAGGING only
      float trackPT_7 = 0;   // t0 + self-tagging (full ITk+HGTD)
      // DIAGNOSTIC: real HGTD track times gated against a PERFECT vertex time
      // (the truth vertex t0, zero error). Isolates how much of the gap to the
      // TDR is our reconstructed t0 rather than the track timing -- our
      // per-track resolutions already match the TDR's "Initial" scenario
      // (35/25/21 ps for 1/2/3 hits), so if this recovers the TDR's improvement
      // factor, the vertex time is the entire remaining deficit.
      float trackPT_8 = 0;
      // Per-jet record of the tracks surviving z0 + dR, for the self-tagging
      // pass below. Untimed tracks are recorded with timed=false: they can
      // never be assigned to a time sub-jet, and are kept unconditionally.
      struct JTrk { double t, res, pt; bool timed, passT0; };
      std::vector<JTrk> jtrks;

      // ######|| Tracks selections ||######
      for (int j = 0; j < (int)jet_tracks_idx[i].size(); ++j) {
        int idex = jet_tracks_idx[i][j];
        if (idex < 0 || idex >= (int)track_pt.GetSize()) continue;
        if (CF) CF->idx_ok++;

        float eta = track_eta[idex];
        float pt  = track_pt[idex];      // GeV (NOTE 3)
        float pt2 = track_pt[idex];      // GeV (NOTE 3)

        float z = track_z0[idex] - recovertex_z[0];

        float vertex_time_truth = truthvertex_t[0];
        float vertex_time_reco  = r->Gaus(vertex_time_truth, 10.0);

        // t30: the track's OWN truth vertex time, smeared by SMEAR_PS. This is
        // what carries the HS/PU separation -- a pileup track inherits its own
        // vertex's time, ~175 ps away from the hard scatter's on average.
        const int tvtx = track_truthVtx_idx[idex];
        const bool has_truth_time = (tvtx >= 0 && tvtx < (int)truthvertex_t.GetSize());
        float track_t30_val = 0.f;
        if (has_truth_time)
          track_t30_val = r->Gaus(truthvertex_t[tvtx], SMEAR_PS);

        float time_cut_reco  = (track_t30_val - vertex_time_reco) / SMEAR_PS;
        float time_cut_truth = (track_t30_val - truthvertex_t[0])  / SMEAR_PS;

        float Dz0para = getNewDzpara(eta, pt);

        float deta = jet_eta[i] - eta;
        float dphi = jet_phi[i] - track_phi[idex];
        if (dphi > TMath::Pi()) dphi -= 2*TMath::Pi();
        float Dr = std::sqrt(deta*deta + dphi*dphi);

        if (pt < 1.0 || pt > 45.) continue;
        if (std::fabs(eta) > 4.0) continue;

        if (CF) CF->pt_eta++;
        trackPT_0 += pt2;

        if (ZCUT == "atlas") {
          // sin(theta) = 1/cosh(eta) exactly, so no extra branch is needed.
          const double sinth = 1.0 / std::cosh(eta);
          if (std::fabs(z * sinth) > ZSINTH) continue;
        } else if (ZCUT == "signif") {
          const double vz = track_var_z0[idex];
          if (!(vz > 0)) continue;
          if (std::fabs(z) / std::sqrt(vz) > ZSIG) continue;
        } else {
          if ((std::fabs(z)/Dz0para) > DZ0_SCALE) continue;   // 0.8 is better (theirs)
        }

        if (CF) CF->zcut++;
        trackPT_1 += pt2;

        if (Dr > TRK_DR) continue;

        if (CF) CF->dr++;
        trackPT_2 += pt2;

        // Record for the self-tagging pass (runs after this track loop).
        {
          const bool trk_ok = (track_timeValid[idex] == 1);
          bool passT0 = true;   // untimed / no vertex time -> kept on z-info
          if (vtx_time_ok && trk_ok) {
            const double dtv = track_time[idex] - vtx_t;
            const double sgv = std::sqrt(track_timeRes[idex]*track_timeRes[idex]
                                         + vtx_res*vtx_res);
            passT0 = (sgv > 0 && std::fabs(dtv / sgv) < T0_NSIGMA);
          }
          jtrks.push_back({trk_ok ? (double)track_time[idex] : 0.0,
                           trk_ok ? (double)track_timeRes[idex] : 0.0,
                           (double)pt2, trk_ok, passT0});
        }

        // Tracks with no truth vertex link (~3.4%) cannot be assigned a t30
        // at all, so they contribute to the untimed sums only. The original has
        // no such case because its ntuple ships a t30 for every track.
        // ── REALISTIC: real HGTD vertex time + real HGTD track time ───────
        // Not part of the reference macro. Same |pull| < 3 structure as the
        // idealised cut, but the pull is formed from the ACTUAL per-track and
        // per-vertex resolutions rather than a flat smearing width.
        //
        // Fallback rule (as specified): if EITHER the vertex time or this
        // track's time is unavailable, the track is judged on z-information
        // alone and KEPT -- never discarded. That is the physically honest
        // statement (the detector cannot reject what it did not measure) and
        // it matches rpt_v5's applyTimeGate convention, which is what makes
        // this scenario the like-for-like comparison against rpt_v5.
        const bool trk_time_ok = (track_timeValid[idex] == 1);
        // Perfect-t0 variant: same track times and same fallback rule, but the
        // reference is the truth vertex time with no error of its own.
        if (!trk_time_ok) {
          trackPT_8 += pt2;
        } else {
          const double dt8 = track_time[idex] - truthvertex_t[0];
          const double sg8 = track_timeRes[idex];
          if (sg8 > 0 && std::fabs(dt8 / sg8) < T0_NSIGMA) trackPT_8 += pt2;
        }
        if (!vtx_time_ok) {
          ++n_trk_zonly_vtx;
          trackPT_5 += pt2;                       // z-info only
        } else if (!trk_time_ok) {
          ++n_trk_zonly_trk;
          trackPT_5 += pt2;                       // z-info only
        } else {
          ++n_trk_timed;
          const double dt   = track_time[idex] - vtx_t;
          const double sig  = std::sqrt(track_timeRes[idex]*track_timeRes[idex]
                                        + vtx_res*vtx_res);
          if (sig > 0 && std::fabs(dt / sig) < T0_NSIGMA) trackPT_5 += pt2;
        }

        if (!has_truth_time) continue;   // no truth vertex -> no t30 to cut on
        if (std::fabs(time_cut_reco)  < 3.0) trackPT_3 += pt2;   // cut 3 is better (theirs)
        if (std::fabs(time_cut_truth) < 3.0) trackPT_4 += pt2;
      }  // track_loop ends

      // ── SELF-TAGGING (ATLAS-TDR-031) ──────────────────────────────────────
      // "does not require the knowledge of the hard-scatter time. The key idea
      //  is to check the consistency of the measured production time for all
      //  tracks associated to the same physics object among themselves ...
      //  finding clusters of tracks within a jet that have compatible times,
      //  and splitting the jet into smaller sub-jets with consistent times."
      //
      // Implemented as single-linkage clustering of the jet's TIMED tracks
      // along the time axis: sort by time, start a new sub-jet whenever the
      // pull to the previous track exceeds SELFTAG_NSIGMA. The surviving
      // sub-jet is the one carrying the most pT -- the jet's own time core --
      // and the out-of-time stragglers are dropped.
      //
      // The TDR's stated limitation is respected: self-tagging "requires
      // physics objects to have at least two tracks with time assigned". With
      // fewer than two, no consistency statement is possible and every track is
      // kept (z-information only), never discarded. The TDR predicts this
      // blunts the method against pileup jets specifically, since most of them
      // have only one track in HGTD acceptance -- n_selftag_* below measures
      // exactly that.
      {
        std::vector<const JTrk*> timed;
        for (const auto& jt : jtrks) if (jt.timed) timed.push_back(&jt);

        double untimed_pt = 0.0, untimed_pt_t0 = 0.0;
        for (const auto& jt : jtrks)
          if (!jt.timed) { untimed_pt += jt.pt; untimed_pt_t0 += jt.pt; }

        if (timed.size() < 2) {
          // No self-tagging possible: keep everything.
          ++n_selftag_toofew;
          for (const auto& jt : jtrks) trackPT_6 += jt.pt;
        } else {
          ++n_selftag_applied;
          std::sort(timed.begin(), timed.end(),
                    [](const JTrk* a, const JTrk* b) { return a->t < b->t; });

          // Single-linkage split on the time pull between adjacent tracks.
          std::vector<std::vector<const JTrk*>> subjets;
          subjets.emplace_back();
          subjets.back().push_back(timed[0]);
          for (size_t q = 1; q < timed.size(); ++q) {
            const double sg = std::sqrt(timed[q]->res*timed[q]->res
                                        + timed[q-1]->res*timed[q-1]->res);
            const double pull = (sg > 0) ? (timed[q]->t - timed[q-1]->t) / sg : 0.0;
            if (std::fabs(pull) > SELFTAG_NSIGMA) subjets.emplace_back();
            subjets.back().push_back(timed[q]);
          }

          // Keep the highest-sum-pT sub-jet.
          size_t best = 0; double bestpt = -1.0;
          for (size_t q = 0; q < subjets.size(); ++q) {
            double sp = 0.0;
            for (auto* jt : subjets[q]) sp += jt->pt;
            if (sp > bestpt) { bestpt = sp; best = q; }
          }
          if (subjets.size() > 1) ++n_selftag_split;

          // TDR: "The RpT variable is then recomputed for each in-time cluster
          // plus the additional tracks that have no time assigned. Then, the jet
          // RpT is defined as the maximum RpT of all sub-jets." Every sub-jet
          // shares the same untimed term and the same denominator, so the
          // maximum-RpT sub-jet is exactly the maximum-sum-pT one selected above.
          trackPT_6 += untimed_pt;      // untimed tracks always survive
          for (auto* jt : subjets[best]) trackPT_6 += jt->pt;
        }
      }

      // COMBINED (ATLAS-TDR-031): "It is also possible to combine both
      // approaches such that when no t0 is found, the self-tagging method is
      // used." A fallback, NOT an intersection -- an earlier revision took the
      // AND of the two cuts, which over-rejected and pushed the combined curve
      // BELOW t0-only, the opposite of the TDR's figure.
      trackPT_7 = vtx_time_ok ? trackPT_5 : trackPT_6;
      ++n_jets_total;
      if (!vtx_time_ok) {
        ++n_jets_fallback;
        // ...and the fallback only CHANGES anything when self-tagging actually
        // removed something, i.e. when the jet had >=2 timed tracks to split.
        if (trackPT_6 < trackPT_2 - 1e-6) ++n_jets_fallback_active;
      }

      float Rpt2 = trackPT_2 / jet_pt[i];
      float Rpt3 = trackPT_3 / jet_pt[i];
      float Rpt4 = trackPT_4 / jet_pt[i];
      float Rpt5 = trackPT_5 / jet_pt[i];
      float Rpt6 = trackPT_6 / jet_pt[i];
      float Rpt7 = trackPT_7 / jet_pt[i];
      float Rpt8 = trackPT_8 / jet_pt[i];

      bool hs = isPaperHS(jet_eta[i], jet_phi[i]);
      bool pu = isPaperPU(jet_eta[i], jet_phi[i]);
      if (hs && !pu) {
        sHS0 += trackPT_0/jet_pt[i]; sHS1 += trackPT_1/jet_pt[i];
        sHS2 += trackPT_2/jet_pt[i]; ++cHS;
        nghost_hs += (long)jet_tracks_idx[i].size(); ++njet_hs;
      }
      if (!hs && pu) {
        sPU0 += trackPT_0/jet_pt[i]; sPU1 += trackPT_1/jet_pt[i];
        sPU2 += trackPT_2/jet_pt[i]; ++cPU;
        nghost_pu += (long)jet_tracks_idx[i].size(); ++njet_pu;
      }

      if (hs && !pu) {
        h_Rpt->Fill(Rpt2);
        h_Rpt_t->Fill(Rpt3);
        h_Rpt_t_truth->Fill(Rpt4);
        h_Rpt_reco->Fill(Rpt5);
        h_Rpt_self->Fill(Rpt6);
        h_Rpt_full->Fill(Rpt7);
        h_Rpt_perft0->Fill(Rpt8);
        ++n_jets_filled_hs;
      }
      if (!hs && pu) {
        h_Rpt_pu->Fill(Rpt2);
        h_Rpt_pu_t->Fill(Rpt3);
        h_Rpt_pu_t_truth->Fill(Rpt4);
        h_Rpt_pu_reco->Fill(Rpt5);
        h_Rpt_pu_self->Fill(Rpt6);
        h_Rpt_pu_full->Fill(Rpt7);
        h_Rpt_pu_perft0->Fill(Rpt8);
        ++n_jets_filled_pu;
      }
    }  // Jet_loop ends
  }  // event loop
  std::cout << "\nnum_DZ events# : " << num_DZ << std::endl;

  // ######|| Save RpT distributions ||######
  gStyle->SetPalette(1);
  TCanvas* c0 = new TCanvas("c0", "c0", 15, 34, 800, 600);
  SetCanvasAttr(c0);
  c0->SetGridy();

  h_Rpt->SetLineColor(2);  h_Rpt->SetMarkerColor(2);  h_Rpt->SetLineWidth(2);
  h_Rpt->SetYTitle("Arbitrary units"); h_Rpt->SetXTitle("R_{p_{T}}");
  h_Rpt_pu->SetLineColor(1); h_Rpt_pu->SetMarkerColor(1); h_Rpt_pu->SetLineWidth(2);
  h_Rpt_pu->SetFillColorAlpha(1, 0.5);
  h_Rpt_pu->SetYTitle("Arbitrary units"); h_Rpt_pu->SetXTitle("R_{p_{T}}");

  h_Rpt_pu->DrawNormalized();
  h_Rpt->DrawNormalized("same");

  TLegend* leg2 = new TLegend(0.35, 0.73, 0.60, 0.85, NULL, "brNDC");
  leg2->SetFillColor(0); leg2->SetBorderSize(0);
  leg2->SetTextFont(42); leg2->SetTextSize(0.04);
  leg2->AddEntry(h_Rpt, "HS jet");
  leg2->AddEntry(h_Rpt_pu, "PU jet");
  leg2->Draw();
  ATLAS_LABEL(0.18, 0.88, 1);
  c0->Print((pre + "RptDist_without_time.pdf").c_str());

  // ######|| Save RpT distributions with time ||######
  TCanvas* c0_1 = new TCanvas("c0_1", "c0_1", 15, 34, 800, 600);
  c0_1->SetGridy();
  h_Rpt_t->SetLineColor(2); h_Rpt_t->SetMarkerColor(2); h_Rpt_t->SetLineWidth(2);
  h_Rpt_t->SetYTitle("Arbitrary units"); h_Rpt_t->SetXTitle("R_{p_{T}}");
  h_Rpt_pu_t->SetLineColor(1); h_Rpt_pu_t->SetMarkerColor(1); h_Rpt_pu_t->SetLineWidth(2);
  h_Rpt_pu_t->SetFillColorAlpha(1, 0.5);
  h_Rpt_pu_t->SetYTitle("Arbitrary units"); h_Rpt_pu_t->SetXTitle("R_{p_{T}}");
  h_Rpt_pu_t->DrawNormalized();
  h_Rpt_t->DrawNormalized("same");

  TLegend* leg2_1 = new TLegend(0.35, 0.73, 0.60, 0.85, NULL, "brNDC");
  leg2_1->SetFillColor(0); leg2_1->SetBorderSize(0);
  leg2_1->SetTextFont(42); leg2_1->SetTextSize(0.04);
  leg2_1->AddEntry(h_Rpt_t, "HS jet");
  leg2_1->AddEntry(h_Rpt_pu_t, "PU jet");
  leg2_1->Draw();
  ATLAS_LABEL(0.18, 0.88, 1);
  c0_1->Print((pre + "RptDist_with_time.pdf").c_str());

  // ── ROC curves (their exact construction: Integral(i, nbins), i from 0) ────
  std::cout << "----------------" << std::endl << "Without time" << std::endl;
  TGraph* g = new TGraph();
  float Ns = h_Rpt->Integral(), Nb = h_Rpt_pu->Integral();
  int nbins = h_Rpt->GetNbinsX();
  int k = 0;
  for (int i = 0; i <= nbins; i++) {
    float S_eff = h_Rpt->Integral(i, nbins) / Ns;
    float B_eff = h_Rpt_pu->Integral(i, nbins) / Nb;
    if (B_eff == 0) continue;
    g->SetPoint(k, S_eff, 1./B_eff);
    k++;
  }
  g->SetLineWidth(4); g->SetLineColor(1); g->SetFillColor(0);

  std::cout << "----------------" << std::endl << "With truth time" << std::endl;
  TGraph* g2 = new TGraph();
  float Ns2 = h_Rpt_t_truth->Integral(), Nb2 = h_Rpt_pu_t_truth->Integral();
  int nbins2 = h_Rpt_t->GetNbinsX();
  int k2 = 0;
  for (int i = 0; i <= nbins2; i++) {
    float S_eff_tr = h_Rpt_t_truth->Integral(i, nbins2) / Ns2;
    float B_eff_tr = h_Rpt_pu_t_truth->Integral(i, nbins2) / Nb2;
    if (B_eff_tr == 0) continue;
    g2->SetPoint(k2, S_eff_tr, 1./B_eff_tr);
    k2++;
  }
  g2->SetLineWidth(4); g2->SetLineColor(2); g2->SetLineStyle(5); g2->SetFillColor(0);

  std::cout << "----------------" << std::endl << "With reco time" << std::endl;
  TGraph* g3 = new TGraph();
  float Ns3 = h_Rpt_t->Integral(), Nb3 = h_Rpt_pu_t->Integral();
  int nbins3 = h_Rpt_t->GetNbinsX();
  int k3 = 0;
  for (int i = 0; i <= nbins3; i++) {
    float S_eff_t = h_Rpt_t->Integral(i, nbins3) / Ns3;
    float B_eff_t = h_Rpt_pu_t->Integral(i, nbins3) / Nb3;
    if (B_eff_t == 0) continue;
    g3->SetPoint(k3, S_eff_t, 1./B_eff_t);
    k3++;
  }
  g3->SetLineWidth(4); g3->SetLineColor(1); g3->SetLineStyle(2); g3->SetFillColor(0);

  std::cout << "----------------" << std::endl << "With REAL HGTD times" << std::endl;
  TGraph* g4 = new TGraph();
  float Ns4 = h_Rpt_reco->Integral(), Nb4 = h_Rpt_pu_reco->Integral();
  int nbins4 = h_Rpt_reco->GetNbinsX();
  int k4 = 0;
  for (int i = 0; i <= nbins4; i++) {
    float S_eff_r = h_Rpt_reco->Integral(i, nbins4) / Ns4;
    float B_eff_r = h_Rpt_pu_reco->Integral(i, nbins4) / Nb4;
    if (B_eff_r == 0) continue;
    g4->SetPoint(k4, S_eff_r, 1./B_eff_r);
    k4++;
  }
  g4->SetLineWidth(4); g4->SetLineColor(kBlue+1); g4->SetFillColor(0);

  // ── ROC sampled at fixed HS-efficiency working points ─────────────────────
  // The TDR figure carries ~7 markers evenly spaced from 0.80 to ~0.98, not one
  // per histogram bin. Scanning bins (makeRoc below) instead samples wherever
  // the R_pT distribution happens to put thresholds, which on our data clusters
  // points unevenly and leaves a hard break across the R_pT = 0 atom (85% of
  // pileup jets sit there, so no threshold exists inside it).
  //
  // Sampling on efficiency instead fixes both: walk the cumulative
  // distributions and linearly interpolate 1/B_eff at each target S_eff. The
  // interpolation across the zero-atom is the physically right thing -- it is
  // what breaking ties at random inside that bin would deliver -- and it makes
  // every curve share identical x-points, so the ratio panel becomes exact
  // instead of relying on TGraph::Eval.
  auto makeRocGrid = [](TH1F* hs, TH1F* hp, const std::vector<double>& targets) {
    TGraph* gr = new TGraph();
    const int nb = hs->GetNbinsX();
    const double Ns = hs->Integral(), Nb = hp->Integral();
    if (Ns <= 0 || Nb <= 0) return gr;

    std::vector<double> S(nb + 2), B(nb + 2);
    for (int i = 0; i <= nb; ++i) {
      S[i] = hs->Integral(i, nb) / Ns;
      B[i] = hp->Integral(i, nb) / Nb;
    }
    int k = 0;
    for (double e : targets) {
      if (e > S[0] || e <= 0) continue;              // efficiency unreachable
      int i = 0;
      while (i <= nb && S[i] > e) ++i;               // first bin at/below target
      if (i == 0 || i > nb) continue;
      const double s_hi = S[i-1], s_lo = S[i];
      const double b_hi = B[i-1], b_lo = B[i];
      const double f = (s_hi > s_lo) ? (s_hi - e) / (s_hi - s_lo) : 0.0;
      const double b = b_hi + f * (b_lo - b_hi);
      if (b <= 0) continue;                          // no background left: stop
      gr->SetPoint(k++, e, 1.0 / b);
    }
    return gr;
  };

  // Working points: 0.80 to 0.975 in steps of 0.025, matching the TDR's marker
  // spacing. Curves terminate naturally wherever the background empties out.
  std::vector<double> effGrid;
  for (double e = 0.800; e <= 0.9751; e += 0.025) effGrid.push_back(e);

  auto makeRoc = [](TH1F* hs, TH1F* hp) {
    TGraph* gr = new TGraph();
    float Ns_ = hs->Integral(), Nb_ = hp->Integral();
    int nb_ = hs->GetNbinsX(), kk = 0;
    for (int i = 0; i <= nb_; i++) {
      float se = hs->Integral(i, nb_) / Ns_;
      float be = hp->Integral(i, nb_) / Nb_;
      if (be == 0) continue;
      gr->SetPoint(kk++, se, 1./be);
    }
    return gr;
  };
  TGraph* g5 = makeRoc(h_Rpt_self, h_Rpt_pu_self);   // self-tagging only
  TGraph* g6 = makeRoc(h_Rpt_full, h_Rpt_pu_full);   // t0 + self-tagging
  TGraph* g7 = makeRoc(h_Rpt_perft0, h_Rpt_pu_perft0);  // real trk t, perfect t0
  g5->SetLineWidth(4); g5->SetLineColor(kGreen+2); g5->SetLineStyle(9); g5->SetFillColor(0);
  g6->SetLineWidth(4); g6->SetLineColor(kRed);     g6->SetFillColor(0);

  // Ratio graphs (reco/no-time and truth/no-time), as in the original.
  TGraph* ratioGraph  = new TGraph();
  TGraph* ratioGraph2 = new TGraph();
  for (int i = 0; i < g3->GetN(); ++i) {
    double x1, y1, x2, y2, x3, y3;
    g ->GetPoint(i, x1, y1);
    g3->GetPoint(i, x3, y3);
    g2->GetPoint(i, x2, y2);
    if (y1 == 0) continue;
    ratioGraph ->SetPoint(i, x3, y3/y1);
    ratioGraph2->SetPoint(i, x2, y2/y1);
  }

  TH2F* h2 = new TH2F("h2", "", 125, 0.8, 1, 100, 0.0, 170);
  h2->GetXaxis()->SetTitle("Efficiency for hard-scatter jets");
  h2->GetYaxis()->SetTitle("Pile-up jet rejection");
  h2->GetXaxis()->SetTitleOffset(1.25);
  h2->SetStats(0);

  // Labels track the ACTUAL configuration -- these were hardcoded to the
  // macro's 2.4-3.8 / 30-50, which silently mislabels any --etamax/--ptmin run.
  char lbl_eta[64], lbl_pt[64];
  std::snprintf(lbl_eta, sizeof(lbl_eta), "%.1f <|#eta|< %.1f", JETA_MIN, JETA_MAX);
  if (JPT_MAX > 1e8) std::snprintf(lbl_pt, sizeof(lbl_pt), "p_{T}> %.0f GeV", JPT_MIN);
  else               std::snprintf(lbl_pt, sizeof(lbl_pt), "%.0f <p_{T}< %.0f GeV", JPT_MIN, JPT_MAX);

  auto drawLabels = [&](double x = 0.60) {
    ATLAS_LABEL(0.18, 0.85, 1);
    myText(x, 0.80, 1, "#sqrt{s}=14 TeV, <#mu>=200");
    myText(x, 0.75, 1, "VBF H #rightarrow invisible");
    myText(x, 0.70, 1, lbl_eta);
    myText(x, 0.65, 1, lbl_pt);
  };

  // ── ROC: no-time vs reco-time ─────────────────────────────────────────────
  TCanvas* c1 = new TCanvas("c1", "c1", 15, 34, 800, 600);
  SetCanvasAttr(c1);
  h2->Draw();
  g->Draw("same");
  g3->Draw("same");
  {
    TLegend* L = new TLegend(0.22, 0.70, 0.59, 0.85, NULL, "brNDC");
    L->SetFillColor(kWhite); L->SetBorderSize(0);
    L->SetTextFont(42); L->SetTextSize(0.04); L->SetFillColor(0);
    L->AddEntry(g,  "R_pT");
    L->AddEntry(g3, "R_pT 30 ps");
    L->Draw();
    drawLabels();
    c1->Print((pre + "t30_reco.pdf").c_str());
  }

  // ── ROC: all three curves ─────────────────────────────────────────────────
  TCanvas* c2 = new TCanvas("c2", "c2", 15, 34, 800, 600);
  SetCanvasAttr(c2);
  h2->Draw();
  g->Draw("same");
  g2->Draw("same");
  g3->Draw("same");
  {
    TLegend* L = new TLegend(0.22, 0.70, 0.59, 0.85, NULL, "brNDC");
    L->SetFillColor(kWhite); L->SetBorderSize(0);
    L->SetTextFont(42); L->SetTextSize(0.04); L->SetFillColor(0);
    L->AddEntry(g,  "R_pT");
    L->AddEntry(g2, "R_pT 30 ps (truth)");
    L->AddEntry(g3, "R_pT 30 ps (reco) ");
    L->Draw();
    drawLabels();
    c2->Print((pre + "truth_reco_t30.pdf").c_str());
  }

  // ── ROC: IDEALISED vs REALISTIC ───────────────────────────────────────────
  // The comparison that sets the baseline: idealised (truth times smeared at
  // SMEAR_PS, vertex smeared at 10 ps, ~97% of tracks timed) against realistic
  // (real HGTD track + vertex times, ~41% of tracks timed, the rest kept on
  // z-information alone). The gap between them is the cost of real
  // reconstruction, and is the number rpt_v5 should be judged against.
  TCanvas* c4 = new TCanvas("c4", "c4", 15, 34, 800, 600);
  SetCanvasAttr(c4);
  h2->Draw();
  g->Draw("same");
  g3->Draw("same");
  g4->Draw("same");
  {
    TLegend* L = new TLegend(0.22, 0.68, 0.62, 0.85, NULL, "brNDC");
    L->SetFillColor(kWhite); L->SetBorderSize(0);
    L->SetTextFont(42); L->SetTextSize(0.035); L->SetFillColor(0);
    L->AddEntry(g,  "ITk only (no time)");
    L->AddEntry(g3, "Idealised (smeared truth t)");
    L->AddEntry(g4, "Realistic (real HGTD t)");
    L->Draw();
    drawLabels();
    c4->Print((pre + "idealised_vs_realistic.pdf").c_str());
  }

  // ── TDR-style figure: the four ATLAS-TDR-031 p.45 curves + ratio panel ────
  {
    TCanvas* ct = new TCanvas("ct", "ct", 15, 34, 800, 700);
    TPad* tp1 = new TPad("tp1", "tp1", 0., 0.3, 1., 1.);
    tp1->SetBottomMargin(0.02); tp1->Draw(); tp1->cd();
    TGraph* gi = makeRocGrid(h_Rpt,      h_Rpt_pu,      effGrid);
    TGraph* gs = makeRocGrid(h_Rpt_self, h_Rpt_pu_self, effGrid);
    TGraph* gt = makeRocGrid(h_Rpt_reco, h_Rpt_pu_reco, effGrid);
    TGraph* gf = makeRocGrid(h_Rpt_full, h_Rpt_pu_full, effGrid);
    gi->SetLineColor(kBlack);   gi->SetLineStyle(2); gi->SetLineWidth(4);
    gs->SetLineColor(kGreen+2); gs->SetLineStyle(9); gs->SetLineWidth(4);
    gt->SetLineColor(kBlue);    gt->SetLineStyle(3); gt->SetLineWidth(4);
    gf->SetLineColor(kRed);     gf->SetLineStyle(1); gf->SetLineWidth(4);
    for (auto* gg : {gi, gs, gt, gf}) { gg->SetMarkerStyle(20); gg->SetMarkerSize(0.7); }
    h2->Draw(); gi->Draw("same"); gs->Draw("same"); gt->Draw("same"); gf->Draw("same");
    h2->GetXaxis()->SetLabelSize(0);   // x labels live on the ratio pad
    h2->GetXaxis()->SetTitleSize(0);
    TLegend* L = new TLegend(0.56, 0.58, 0.95, 0.86, NULL, "brNDC");
    L->SetFillStyle(0);          // transparent: an opaque box here hid the
    L->SetBorderSize(0);         // curves between eff 0.85 and 0.90
    L->SetTextFont(42); L->SetTextSize(0.036);
    L->AddEntry(gi, "ITk", "l");
    L->AddEntry(gs, "ITk+HGTD (self-tagging only)", "l");
    L->AddEntry(gt, "ITk+HGTD (t_{0} only)", "l");
    L->AddEntry(gf, "ITk+HGTD", "l");
    L->Draw();
    // Labels upper-LEFT, as in the TDR figure (legend now occupies upper-right).
    // Kept left of x=0.55 so the upper-right legend never overlaps it.
    { TLatex t; t.SetNDC(); t.SetTextFont(72); t.SetTextSize(0.045);
      t.DrawLatex(0.20, 0.86, "ATLAS"); t.SetTextFont(42);
      t.DrawLatex(0.325, 0.86, "Simulation Internal"); }
    myText(0.20, 0.795, 1, "VBF H #rightarrow invisible, <#mu>=200");
    myText(0.20, 0.735, 1, lbl_pt);
    myText(0.20, 0.675, 1, lbl_eta);

    ct->cd();
    TPad* tp2 = new TPad("tp2", "tp2", 0, 0.05, 1, 0.3);
    tp2->SetTopMargin(0); tp2->SetBottomMargin(0.25); tp2->Draw(); tp2->cd();
    TH2F* hr = new TH2F("hr", "", 125, 0.8, 1, 100, 0.9, 1.6);
    hr->GetXaxis()->SetTitle("HS efficiency");
    hr->GetYaxis()->SetTitle("ratio");
    hr->GetYaxis()->SetTitleSize(0.11); hr->GetXaxis()->SetTitleSize(0.11);
    hr->GetXaxis()->SetLabelSize(0.09);  hr->GetYaxis()->SetLabelSize(0.09);
    hr->GetYaxis()->SetTitleOffset(0.33); hr->SetStats(0);
    hr->Draw();
    // Every curve shares gi's x-points by construction, so the ratio is exact
    // (no TGraph::Eval interpolation) -- and it stops where either curve stops.
    auto ratioTo = [&](TGraph* num, Color_t col, Style_t sty) {
      TGraph* rg = new TGraph(); int n = 0;
      for (int i = 0; i < gi->GetN() && i < num->GetN(); ++i) {
        double x, y, xn, yn;
        gi ->GetPoint(i, x,  y);
        num->GetPoint(i, xn, yn);
        if (y <= 0 || std::fabs(x - xn) > 1e-9) continue;
        rg->SetPoint(n++, x, yn / y);
      }
      rg->SetLineColor(col); rg->SetLineStyle(sty); rg->SetLineWidth(3);
      rg->SetMarkerStyle(20); rg->SetMarkerSize(0.6); rg->SetMarkerColor(col);
      rg->Draw("same");
      return rg;
    };
    ratioTo(gs, kGreen+2, 9);
    ratioTo(gt, kBlue,    3);
    ratioTo(gf, kRed,     1);
    TLine* one = new TLine(0.8, 1.0, 1.0, 1.0);
    one->SetLineStyle(2); one->Draw("same");
    ct->Print((pre + "tdr_style.pdf").c_str());
    h2->GetXaxis()->SetLabelSize(0.04);   // restore for the canvases below
    h2->GetXaxis()->SetTitleSize(0.04);
  }

  // ── ROC + ratio panel ─────────────────────────────────────────────────────
  TCanvas* c3 = new TCanvas("c3", "c3", 15, 34, 800, 600);
  TPad* p1 = new TPad("p1", "p1", 0., 0.3, 1., 1.);
  p1->Draw();
  p1->SetBottomMargin(0.0);
  p1->cd();
  h2->Draw();
  g->Draw("same");
  g2->Draw("same");
  g3->Draw("same");
  {
    TLegend* L = new TLegend(0.22, 0.70, 0.59, 0.85, NULL, "brNDC");
    L->SetFillColor(kWhite); L->SetBorderSize(0);
    L->SetTextFont(42); L->SetTextSize(0.04); L->SetFillColor(0);
    L->AddEntry(g,  "R_pT");
    L->AddEntry(g2, "R_pT 30 ps (truth)");
    L->AddEntry(g3, "R_pT 30 ps (reco)");
    L->Draw();
    drawLabels();
  }
  c3->cd();
  TPad* p2 = new TPad("pad2", "pad2", 0, 0.05, 1, 0.3);
  p2->SetTopMargin(0);
  p2->SetBottomMargin(0.2);
  p2->Draw();
  p2->cd();
  TH2F* h3 = new TH2F("h3", "", 125, 0.8, 1, 125, 0.0, 6.0);
  h3->GetXaxis()->SetTitle("Efficiency for hard-scatter jets");
  h3->GetYaxis()->SetTitle("Ratio to ITK");
  h3->GetXaxis()->SetRangeUser(0.8, 1.0);
  h3->GetYaxis()->SetTitleSize(0.10);
  h3->GetXaxis()->SetTitleSize(0.10);
  h3->GetXaxis()->SetLabelSize(0.075);
  h3->GetYaxis()->SetLabelSize(0.075);
  h3->GetYaxis()->SetTitleOffset(0.35);
  h3->SetStats(0);
  ratioGraph->SetLineStyle(2);
  ratioGraph2->SetLineStyle(5);
  ratioGraph2->SetLineColor(2);
  h3->Draw();
  ratioGraph->Draw("same");
  ratioGraph2->Draw("same");
  c3->Print((pre + "ratio_t30.pdf").c_str());

  // ── Console summary ───────────────────────────────────────────────────────
  auto rejAt = [](TGraph* gr, double eff) {
    double best = -1, bx = 0, by = 0;
    for (int i = 0; i < gr->GetN(); ++i) {
      double x, y; gr->GetPoint(i, x, y);
      if (x >= eff && (best < 0 || x < bx)) { best = 1; bx = x; by = y; }
    }
    return best < 0 ? 0.0 : by;
  };
  std::printf("\n================================================================\n");
  std::printf("  rpt_v6 — direct port of myJet_ana_fr.C\n");
  std::printf("================================================================\n");
  std::printf("  Events read              : %8lld\n", (long long)seen);
  std::printf("  Events passing |dz|<%.1fmm: %8d\n", MAX_DZ, num_DZ);
  std::printf("  HS jets filled           : %8ld\n", n_jets_filled_hs);
  std::printf("  PU jets filled           : %8ld\n", n_jets_filled_pu);
  std::printf("\n  Stage decomposition of <R_pT> (mean over jets):\n");
  std::printf("    %-10s %10s %10s %10s %12s\n", "", "pT/eta", "+z assoc", "+dR cone", "<ghost trk>");
  if (cHS) std::printf("    %-10s %10.3f %10.3f %10.3f %12.1f\n", "HS jets",
                       sHS0/cHS, sHS1/cHS, sHS2/cHS, njet_hs ? (double)nghost_hs/njet_hs : 0.0);
  if (cPU) std::printf("    %-10s %10.3f %10.3f %10.3f %12.1f\n", "PU jets",
                       sPU0/cPU, sPU1/cPU, sPU2/cPU, njet_pu ? (double)nghost_pu/njet_pu : 0.0);

  std::printf("\n  R_pT (no timing) shape -- what sets the ITk-only baseline:\n");
  std::printf("    HS jets: <R_pT> = %.3f   frac in R_pT~0 bin = %.3f\n",
              h_Rpt->GetMean(),
              h_Rpt->GetEntries() ? h_Rpt->GetBinContent(h_Rpt->FindBin(0.0))/h_Rpt->GetEntries() : 0.0);
  std::printf("    PU jets: <R_pT> = %.3f   frac in R_pT~0 bin = %.3f\n",
              h_Rpt_pu->GetMean(),
              h_Rpt_pu->GetEntries() ? h_Rpt_pu->GetBinContent(h_Rpt_pu->FindBin(0.0))/h_Rpt_pu->GetEntries() : 0.0);
  std::printf("    (dz0scale=%.2f, trk-jet dR<%.2f)\n", DZ0_SCALE, TRK_DR);
  std::printf("\n  PU rejection (1/mistag) at fixed HS efficiency:\n");
  std::printf("    %-22s %8s %8s %8s %8s\n", "scenario", "0.85", "0.90", "0.93", "0.95");
  std::printf("    %-22s %8.1f %8.1f %8.1f %8.1f\n", "R_pT (no time)",
              rejAt(g,0.85), rejAt(g,0.90), rejAt(g,0.93), rejAt(g,0.95));
  std::printf("    %-22s %8.1f %8.1f %8.1f %8.1f\n", "R_pT 30ps (reco)",
              rejAt(g3,0.85), rejAt(g3,0.90), rejAt(g3,0.93), rejAt(g3,0.95));
  std::printf("    %-22s %8.1f %8.1f %8.1f %8.1f\n", "R_pT 30ps (truth)",
              rejAt(g2,0.85), rejAt(g2,0.90), rejAt(g2,0.93), rejAt(g2,0.95));
  std::printf("    %-22s %8.1f %8.1f %8.1f %8.1f\n", "ITk+HGTD (t0 only)",
              rejAt(g4,0.85), rejAt(g4,0.90), rejAt(g4,0.93), rejAt(g4,0.95));
  std::printf("    %-22s %8.1f %8.1f %8.1f %8.1f\n", "ITk+HGTD (self-tag)",
              rejAt(g5,0.85), rejAt(g5,0.90), rejAt(g5,0.93), rejAt(g5,0.95));
  std::printf("    %-22s %8.1f %8.1f %8.1f %8.1f\n", "ITk+HGTD (full)",
              rejAt(g6,0.85), rejAt(g6,0.90), rejAt(g6,0.93), rejAt(g6,0.95));
  std::printf("    %-22s %8.1f %8.1f %8.1f %8.1f\n", "  [diag] perfect t0",
              rejAt(g7,0.85), rejAt(g7,0.90), rejAt(g7,0.93), rejAt(g7,0.95));
  {
    const long st = n_selftag_applied + n_selftag_toofew;
    std::printf("\n  Self-tagging applicability (TDR: needs >=2 timed tracks):\n");
    std::printf("    jets with >=2 timed tracks : %8ld / %-8ld (%.1f%%)\n",
                n_selftag_applied, st, st ? 100.0*n_selftag_applied/st : 0.0);
    std::printf("    jets with <2  (no self-tag): %8ld (%.1f%%)\n",
                n_selftag_toofew, st ? 100.0*n_selftag_toofew/st : 0.0);
    std::printf("    jets actually split in time: %8ld (%.1f%% of applicable)\n",
                n_selftag_split,
                n_selftag_applied ? 100.0*n_selftag_split/n_selftag_applied : 0.0);
    std::printf("\n  Combined-vs-t0only reach (fallback fires only without a vertex t0):\n");
    std::printf("    jets in events with NO vertex t0 : %8ld / %-8ld (%.1f%%)\n",
                n_jets_fallback, n_jets_total,
                n_jets_total ? 100.0*n_jets_fallback/n_jets_total : 0.0);
    std::printf("    ...where self-tagging changed RpT: %8ld (%.1f%% of all jets)\n",
                n_jets_fallback_active,
                n_jets_total ? 100.0*n_jets_fallback_active/n_jets_total : 0.0);
  }
  {
    const long tot = n_trk_timed + n_trk_zonly_trk + n_trk_zonly_vtx;
    std::printf("\n  Realistic-scenario timing availability:\n");
    std::printf("    events with a valid HGTD vertex time : %8ld / %-8d (%.1f%%)\n",
                n_evt_vtx_timed, num_DZ, num_DZ ? 100.0*n_evt_vtx_timed/num_DZ : 0.0);
    std::printf("    tracks actually time-gated           : %8ld / %-8ld (%.1f%%)\n",
                n_trk_timed, tot, tot ? 100.0*n_trk_timed/tot : 0.0);
    std::printf("    tracks on z-info only (no track time): %8ld (%.1f%%)\n",
                n_trk_zonly_trk, tot ? 100.0*n_trk_zonly_trk/tot : 0.0);
    std::printf("    tracks on z-info only (no vtx time)  : %8ld (%.1f%%)\n",
                n_trk_zonly_vtx, tot ? 100.0*n_trk_zonly_vtx/tot : 0.0);
  }
  {
    // Dump the working-point grid actually realised, so a missing point is
    // visible as a missing row rather than silently as a gap in the curve.
    auto dump = [&](const char* name, TH1F* hs, TH1F* hp) {
      TGraph* gg = makeRocGrid(hs, hp, effGrid);
      std::printf("    %-18s", name);
      for (int i = 0; i < gg->GetN(); ++i) { double x,y; gg->GetPoint(i,x,y);
        std::printf(" %.3f:%.0f", x, y); }
      std::printf("   (%d pts)\n", gg->GetN());
    };
    {
    auto pct = [](long long a, long long b) { return b > 0 ? 100.0 * a / b : 0.0; };
    auto cfRow = [&](const char* nm, const TrkCutflow& C) {
      std::printf("    %-3s jets %8lld | ghost %6.2f/jet | pT,eta %5.1f%% | z-assoc %5.1f%%"
                  " | dR %5.1f%%  ->  %.2f trk/jet\n",
                  nm, C.jets, C.jets ? double(C.ghost) / C.jets : 0.0,
                  pct(C.pt_eta, C.idx_ok), pct(C.zcut, C.pt_eta), pct(C.dr, C.zcut),
                  C.jets ? double(C.dr) / C.jets : 0.0);
    };
    std::printf("\n  Track cutflow by jet truth-type (zcut=%s):\n", ZCUT.c_str());
    cfRow("HS", CF_HS);
    cfRow("PU", CF_PU);
  }

  std::printf("\n  ROC working points (eff:rejection):\n");
    dump("ITk",        h_Rpt,      h_Rpt_pu);
    // The reference macro's own timed curve: t30 = Gaus(track truth vtx time,
    // 30 ps) against a vertex time of Gaus(truth, 10 ps), cut at a FLAT
    // |dt|/30 < 3 -- not a quadrature-sum pull and not inflated. This is the
    // dotted "R_pT 30 ps" on the slides, so it is the only one of these curves
    // directly comparable to them.
    dump("R_pT 30ps",  h_Rpt_t,    h_Rpt_pu_t);
    dump("self-tag",   h_Rpt_self, h_Rpt_pu_self);
    dump("t0 only",    h_Rpt_reco, h_Rpt_pu_reco);
    dump("full",       h_Rpt_full, h_Rpt_pu_full);
  }
  std::printf("\n  Wrote %s*.pdf\n", pre.c_str());
  std::printf("================================================================\n");
  return 0;
}
