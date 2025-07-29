#include <Rtypes.h>
#include <TRandom1.h>
#include <TCanvas.h>
#include <TChain.h>
#include <TColor.h>
#include <TDirectory.h>
#include <TEfficiency.h>
#include <TFile.h>
#include <TF1.h>
#include <TGaxis.h>
#include <TGraph.h>
#include <TGraphAsymmErrors.h>
#include <TGraphErrors.h>
#include <TH1.h>
#include <TH2.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TLine.h> 
#include <TProfile.h>
#include <TMath.h>
#include <TTreeReader.h>
#include <TTreeReaderArray.h>
#include <TString.h>
#include <TStyle.h>
#include <TVector2.h>

// Boost Headers
#include <boost/filesystem.hpp>

// Standard C++ Library Headers
#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <ostream>
#include <utility>
#include <vector>
#include <map>    
#include <memory> 

namespace myutl {
  static bool debug = false;
  static int c01          = kP10Blue  ;  
  static int c02          = kP10Yellow;  
  static int c03          = kP10Red   ;  
  static int c04          = kP10Gray  ; 
  static int c05          = kP10Violet;  
  static int c06          = kP10Brown ;  
  static int c07          = kP10Orange;  
  static int c08          = kP10Green ; 
  static int c09          = kP10Ash   ;  
  static int c10          = kP10Cyan  ; 
  static std::vector<int> colors = {c01, c02, c03, c04, c05, c06, c07, c08, c09, c10};

  // cut variables
  static int    min_jets          =  1;    // min number of jets
  static double min_jetpt         =  30.0; // self explanatory
  static double min_abs_eta_jet   =  2.38;  // min eta for a jet to be considered "forward"
  static double min_abs_eta_track =  2.38;  // min eta for a track to be considered "forward"
  static double min_track_pt      =  1.0;  // track_pt > 1.0 GeV
  static double max_vtx_dz        =  2.0;  // max error for reco HS vertex z
  static double max_nsigma        =  3.0;  // how close a track can be to PV

  const double min_hgtd_eta = 2.38, max_hgtd_eta = 4.0;

  const double diff_min = -1000.0, diff_max = 1000.0;
  const double diff_width = 2.0;

  const double purity_min = 0, purity_max = 1;
  const double purity_width = 0.05;

  const double fjet_min = 1, fjet_max = 31, fold_fjet = 5;
  const double fjet_width = 1.0;

  const double track_min = 0, track_max = 100, fold_track = 20;
  const double track_width = 1.0;

  const double pu_track_min = track_min, pu_track_max = track_max, fold_hs_track = fold_track;
  const double pu_track_width = track_width;

  const double hs_track_min = track_min, hs_track_max = track_max, fold_pu_track = fold_track;
  const double hs_track_width = track_width;

  const double pu_frac_min = 0, pu_frac_max = 1 , fold_pu_frac = pu_frac_max;
  const double pu_frac_width = 0.05;

  const double z_min = -200, z_max = 200, fold_z = 100;
  const double z_width = 10.0;

  
  static TF1 *dgaus_fitfunc = new TF1("dgaus", "[1]*TMath::Exp(-0.5 * ((x - [0]) / [3])^2) + [2]*TMath::Exp(-0.5 * ((x - [0]) / [4])^2)");

  static TF1* create_fit(TH1D *hist, bool fixbkg) {
    dgaus_fitfunc->SetParNames("Mean", "Norm1", "Norm2", "Sigma1", "Sigma2");
    TF1* fit = new TF1("dgaus_fit", "dgaus", hist->GetXaxis()->GetXmin(), hist->GetXaxis()->GetXmax());
    fit->SetParameters(0, 0.8*hist->GetMaximum(), 1e-2*hist->GetMaximum(), 26.0, 175.0);
    fit->FixParameter(0, 0);
    fit->SetParLimits(1, 0, 1.E6); // can only be positive
    fit->SetParLimits(2, 0, 1.E6); // can only be positive
    fit->SetParLimits(3, 0, 1.E6); // can only be positive
    if (fixbkg)
      fit->FixParameter(4, 175);
    else
      fit->SetParLimits(4, 0, 1.E6); // realistic values
    fit->SetLineColor(kRed);
    fit->SetLineWidth(2);
    fit->SetNpx(1000);
    hist->Fit(fit,"RQ");
    return fit;
  }
  
  /// enums and utilities for them
  enum ScoreType {
    HGTD = 0, SIG = 1, SUMPT2 = 3, TRKPT = 2, SIGTRKPT = 4,
    MAXPT = 5, IDEAL = 6, MAXHS = 7, JETPTDR = 8, TRKPTDR = 9,
    INVALID = -99 };

  static std::vector<ScoreType> enum_vec = {
    HGTD, TRKPT, MAXHS,
  }; // valid values

  static const char* toString(myutl::ScoreType score) {
    switch (score) {
    case ScoreType::HGTD:     return "HGTD Algorithm";
    case ScoreType::SIG:      return "exp(-|s|)";
    case ScoreType::SUMPT2:   return "Sum(p_{T})^{2}";
    case ScoreType::TRKPT:    return "Track p_{T}";
    case ScoreType::SIGTRKPT: return "Track p_{T}*exp(-|s|)";
    case ScoreType::MAXPT:    return "Cluster with max p_{T}";
    case ScoreType::IDEAL:    return "Ideal cluster choice";
    case ScoreType::MAXHS:    return "Max HS Track Clsuter";
    case ScoreType::JETPTDR:  return "exp(-min dR)*Jet p_{T}";
    case ScoreType::TRKPTDR:  return "exp(-min dR)*Track p_{T}";
    default:                  return "INVALID";
    }
  }

  enum FitParamFields { MEAN = 0, SIGMA = 1, CORE = 2, BKG = 3, RATIO = 4 };
  static std::vector<FitParamFields> fitparam_vec = {SIGMA, CORE, BKG, RATIO, MEAN};
  
  static const char* toString(FitParamFields key) {
    switch (key) {
    case FitParamFields::SIGMA: return "SIGMA";
    case FitParamFields::CORE:  return "CORE";
    case FitParamFields::BKG:   return "BKG";
    case FitParamFields::RATIO: return "RATIO";
    case FitParamFields::MEAN:  return "MEAN";
    default:                    return "INVALID";
    }
  }

  struct Cluster {
    std::vector<double> values;
    std::vector<double> sigmas;
    std::vector<double> all_times;
    std::vector<int> track_indices;
    std::map<ScoreType,double> scores;
    double purity = 0.0;
    bool max_pt_cluster = false;
    bool wasMerged = false;
    int nConstituents=1;

    bool operator==(const Cluster& other) {
      bool same_values = values.at(0) == other.values.at(0);
      bool same_sigmas = sigmas.at(0) == other.sigmas.at(0);
      bool same_consts = nConstituents == other.nConstituents;
      // this SHOULD be sufficienct
      return same_consts and same_sigmas and same_values;
    }

    bool operator!=(const Cluster& other) {
      return !(*this == other);
    }      
  };

  struct BranchPointerWrapper {
    TTreeReader& reader;

    TTreeReaderValue<float> weight;

    TTreeReaderArray<float> track_z0;
    TTreeReaderArray<float> track_pt;
    TTreeReaderArray<float> track_eta;
    TTreeReaderArray<float> track_phi;
    TTreeReaderArray<float> track_var_z0;
    TTreeReaderArray<float> track_time;
    TTreeReaderArray<float> track_time_res;
    TTreeReaderArray<int>   track_time_valid;
    TTreeReaderArray<int>   track_to_truthvtx;
    TTreeReaderArray<int>   track_to_particle;
    TTreeReaderArray<bool>  track_quality;
    TTreeReaderArray<int>   track_hgtd_hits;
    TTreeReaderArray<int>   track_prim_hits;
    TTreeReaderArray<float> track_near_idx;
    TTreeReaderArray<float> track_near_sig;

    TTreeReaderArray<float> truth_vtx_z;
    TTreeReaderArray<float> truth_vtx_time;
    TTreeReaderArray<bool>  truth_vtx_ishs;

    TTreeReaderArray<float> reco_vtx_z;
    TTreeReaderArray<float> reco_vtx_time;
    TTreeReaderArray<float> reco_vtx_timeRes;
    TTreeReaderArray<int>   reco_vtx_valid;

    TTreeReaderArray<float> topojet_pt;
    TTreeReaderArray<float> topojet_eta;
    TTreeReaderArray<float> topojet_phi;

    TTreeReaderArray<float> truthhsjet_pt;
    TTreeReaderArray<float> truthhsjet_eta;

    TTreeReaderArray<float> particle_t;

  BranchPointerWrapper(TTreeReader& r)
    : reader           (r                          ),
      weight           (r, "weight"                ),
      track_z0         (r, "Track_z0"              ),
      track_pt         (r, "Track_pt"              ),
      track_eta        (r, "Track_eta"             ),
      track_phi        (r, "Track_phi"             ),
      track_var_z0     (r, "Track_var_z0"          ),
      track_time       (r, "Track_time"            ),
      track_time_res   (r, "Track_timeRes"         ),
      track_time_valid (r, "Track_hasValidTime"    ),
      track_quality    (r, "Track_quality"         ),
      track_to_truthvtx(r, "Track_truthVtx_idx"    ),
      track_to_particle(r, "Track_truthPart_idx"   ),
      track_hgtd_hits  (r, "Track_nHGTDHits"       ),
      track_prim_hits  (r, "Track_nHGTDPrimaryHits"),
      track_near_idx   (r, "Track_nearestVtx_idx"  ),
      track_near_sig   (r, "Track_nearestVtx_sig"  ),
      truth_vtx_z      (r, "TruthVtx_z"            ),
      truth_vtx_time   (r, "TruthVtx_time"         ),
      truth_vtx_ishs   (r, "TruthVtx_isHS"         ),
      reco_vtx_z       (r, "RecoVtx_z"             ),
      reco_vtx_time    (r, "RecoVtx_time"          ),
      reco_vtx_timeRes (r, "RecoVtx_timeRes"       ),
      reco_vtx_valid   (r, "RecoVtx_hasValidTime"  ),
      topojet_pt       (r, "AntiKt4EMTopoJets_pt"  ),
      topojet_eta      (r, "AntiKt4EMTopoJets_eta" ),
      topojet_phi      (r, "AntiKt4EMTopoJets_phi" ),
      truthhsjet_pt    (r, "TruthHSJet_pt"         ),
      truthhsjet_eta   (r, "TruthHSJet_eta"        ),
      particle_t       (r, "TruthPart_prodVtx_time")
    {}

    bool pass_basic_cuts() {
      if (topojet_pt.GetSize() < min_jets) {
	if (debug) std::cout << "Skipping low jet event" << std::endl;
        return false;
      }
    
      // check reco HS vertex is with 2mm of truth HS vertex
      if(std::abs(truth_vtx_z[0] - reco_vtx_z[0]) > max_vtx_dz) {
	if(debug)	std::cout << "Skipping event due to incorrect HS vertex" << std::endl;
	return false;
      }
      
      return true;
    }

    bool pass_jet_pt_cut() { 
      int passptcount = 0, passptetacount = 0;
      for(int jet_idx = 0; jet_idx < truthhsjet_eta.GetSize(); ++jet_idx) {
	float eta = truthhsjet_eta[jet_idx], pt = truthhsjet_pt[jet_idx];
	passptcount += (pt > min_jetpt) ? 1 : 0;
	passptetacount += (pt > min_jetpt) and (eta > min_abs_eta_jet) ? 1 : 0;
        if (passptcount >= 2 and passptetacount >= 1)
          return true;
      }
      if (debug) std::cout << "Skipping pt cut event" << std::endl;
      return false;
    }

    void count_forward_jets(int& nForwardJet) {
    nForwardJet = 0;
    for(int jet_idx = 0; jet_idx < this->topojet_eta.GetSize(); ++jet_idx) {
      float jet_eta = this->topojet_eta[jet_idx];
      if (std::abs(jet_eta) > min_abs_eta_jet)
	nForwardJet++;
    }
  }
    
    void count_forward_tracks(int& nFTrack, int& nFTrack_HS, int& nFTrack_PU,
			      std::vector<int> tracks) {
      nFTrack = 0; nFTrack_HS = 0; nFTrack_PU = 0;
      for(auto trk_idx: tracks) {
	auto eta = this->track_eta[trk_idx];
	auto valid = this->track_time_valid[trk_idx] == 1;
	auto quality = this->track_quality[trk_idx] == true;
	auto pt = this->track_pt[trk_idx];
	auto truthvtx = this->track_to_truthvtx[trk_idx];
	// already know these pass association
	if (std::abs(eta) > min_abs_eta_track and
	    pt > min_track_pt and valid and quality) {
	  nFTrack++;
	  if (truthvtx != -1 and this->truth_vtx_ishs[truthvtx])
	    nFTrack_HS++;
	  else
	    nFTrack_PU++;
	}
      }
    }

    double calc_jetpt_dr_score(int trk_idx) {
      double
	trk_eta = this->track_eta[trk_idx],
	trk_phi = this->track_phi[trk_idx];


      double min_dR = 1e6;
      int min_idx = -1;
      for (int jet_idx=0; jet_idx < this->topojet_eta.GetSize(); ++jet_idx) {
	double
	  jet_eta = this->topojet_eta[jet_idx],
	  jet_phi = this->topojet_phi[jet_idx];

	double
	  deta = jet_eta-trk_eta,
	  dphi = TVector2::Phi_mpi_pi(jet_phi - trk_phi);

	double this_dR = std::sqrt(deta*deta + dphi*dphi);
	if (this_dR < min_dR) {
	  min_dR = this_dR;
	  min_idx = jet_idx;
	}
      }
      double returnScore = this->topojet_pt[min_idx]*std::exp(-min_dR);
      return returnScore;
    }

    double calc_trkpt_dr_score(int trk_idx) {
      double
	trk_eta = this->track_eta[trk_idx],
	trk_phi = this->track_phi[trk_idx];


      double min_dR = 1e6;
      for (int jet_idx=0; jet_idx < this->topojet_eta.GetSize(); ++jet_idx) {
	double
	  jet_eta = this->topojet_eta[jet_idx],
	  jet_phi = this->topojet_phi[jet_idx];

	double
	  deta = jet_eta-trk_eta,
	  dphi = TVector2::Phi_mpi_pi(jet_phi - trk_phi);

	double this_dR = std::sqrt(deta*deta + dphi*dphi);
	if (this_dR < min_dR) {
	  min_dR = this_dR;
	}
      }
      double returnScore = this->track_pt[trk_idx]*std::exp(-min_dR);
      return returnScore;
    }
  };

  struct FitParams {
    TH1D* mean_dist;
    TH1D* sigma_dist;
    TH1D* core_amp_dist;
    TH1D* back_amp_dist;
    TH1D* amp_ratio_dist;

    FitParams(const char* title, const char* times, const char* name, 
	      double x_min, double fold_val, double fold_max,
	      double x_width, int color) {
      auto full_title = Form("%%s vs %s (%s);%s;%%s", title, times, title);
      mean_dist      = new TH1D(Form("mean_dist_%s", name),
				Form(full_title, "Common #mu", "#mu"),
				(int)((fold_max-x_min)/x_width), x_min, fold_max);
    
      sigma_dist     = new TH1D(Form("sigma_dist_%s", name),
				Form(full_title, "Core #sigma", "Core #sigma"),
				(int)((fold_max-x_min)/x_width), x_min, fold_max);
    
      core_amp_dist  = new TH1D(Form("amp1_dist_%s", name),
				Form(full_title, "Core Amplitude", "Amplitude"),
				(int)((fold_max-x_min)/x_width), x_min, fold_max);
    
      back_amp_dist  = new TH1D(Form("amp2_dist_%s", name),
				Form(full_title, "Bkg Amplitude", "Amplitude"),
				(int)((fold_max-x_min)/x_width), x_min, fold_max);
    
      amp_ratio_dist = new TH1D(Form("ampratio_dist_%s", name),
				Form(full_title, "Core Amp/Bkg Amp", "A.U."),
				(int)((fold_max-x_min)/x_width), x_min, fold_max);

      sigma_dist->SetMaximum(175);
      sigma_dist->SetMinimum(5);

      // mean_dist->SetMaximum( 1.0);
      // mean_dist->SetMinimum(-1.0);
    
      for (auto hist: {mean_dist, sigma_dist, core_amp_dist, back_amp_dist, amp_ratio_dist}) {
	hist->SetLineWidth(2);
	hist->SetLineColor(color);
      }
    }

    void fill_each(int idx, TF1* fit) {
      double sigma1 = std::min(fit->GetParameter("Sigma1"), fit->GetParameter("Sigma2"));
      // if (fit->GetChisquare()/fit->GetNDF() < 0.4)
      // 	return;
      if (sigma1 < 1) sigma1 = 175.0;
      
      sigma_dist->SetBinContent(idx+1,sigma1);
      sigma_dist->SetBinError(idx+1,fit->GetParError("Sigma1"));
    
      double amp1     = fit->GetParameter("Norm1"), amp2     = fit->GetParameter("Norm2");
      double amp1_err = fit->GetParError ("Norm1"), amp2_err = fit->GetParError ("Norm2");
      double ampratio     = amp1/amp2;
      double ampratio_err = std::sqrt(pow((1.0/amp2)*amp1_err,2)+pow(-amp2_err*ampratio/amp2,2));

      core_amp_dist->SetBinContent(idx+1,amp1);   
      core_amp_dist->SetBinError(idx+1,amp1_err); 
    
      back_amp_dist->SetBinContent(idx+1,amp2);   
      back_amp_dist->SetBinError(idx+1,amp2_err); 

      amp_ratio_dist->SetBinContent(idx+1,ampratio);   
      amp_ratio_dist->SetBinError(idx+1,ampratio_err); 
    }

    TH1D* from_enum(FitParamFields fit) {
      switch (fit) {
      case FitParamFields::MEAN:  return mean_dist;
      case FitParamFields::SIGMA: return sigma_dist;
      case FitParamFields::CORE:  return core_amp_dist;
      case FitParamFields::BKG:   return back_amp_dist;
      case FitParamFields::RATIO: return amp_ratio_dist;
      default:                    return nullptr;
      }
    }
  };

  struct PlotObj {
    double bin_width;
    double fold_value, fold_max;
    TString fname;
    const char* xtitle;
    const char* ytitle;
    const char*  times;
    double x_min, x_max, x_wid;
    double y_min, y_max, y_wid;
    std::map<ScoreType, FitParams> params; 
    std::unique_ptr<TH1D> eff_total;
    std::map<ScoreType, std::unique_ptr<TH1D>> eff_pass;
    std::map<ScoreType, std::unique_ptr<TEfficiency>> efficiency; 
    std::map<ScoreType, std::unique_ptr<TH2D>> hist;
    std::map<ScoreType, std::unique_ptr<TH2D>> purity;
    std::map<ScoreType, std::vector<std::unique_ptr<TH1D>>> slices_hists;
    std::map<ScoreType, std::vector<std::unique_ptr<TF1>>> slices_fits;
    std::unique_ptr<TLegend> algolegend;

    PlotObj(const char* title, const char* times,
	    const char* fname,
	    double x_min, double x_max, double x_wid,
	    double y_min, double y_max, double y_wid,
	    double p_min, double p_max, double p_wid,
	    double fold_val, double fold_max)
      : bin_width(x_wid), fname(fname),
	xtitle(title), ytitle("#Delta t[ps]"), times(times),
	x_min(x_min), x_max(x_max), x_wid(x_wid),
	y_min(y_min), y_max(y_max), y_wid(y_wid),
	fold_value(fold_val), fold_max(fold_max)
    {// Allocate with make_unique for owned objects
      algolegend = std::make_unique<TLegend>(0.6, 0.7, 0.9, 0.9);
      algolegend->SetTextSize(0.03);
      TString name(title);
      name.ReplaceAll(" ", "");
      for (ScoreType score : enum_vec) {
	hist[score]    = std::make_unique<TH2D>(Form("%s_%s", name.Data(), toString(score)),
						Form("%s t_{0} - TruthVtx t_{0} vs %s (%s);%s;#Delta t[ps]",
						     toString(score), title, times, title),
						(int)((x_max - x_min) / x_wid), x_min, x_max,
						(int)((y_max - y_min) / y_wid), y_min, y_max);

	purity[score]  = std::make_unique<TH2D>(Form("cluster_purity_%s_%s", name.Data(), toString(score)),
						Form("%s Cluster Purity vs %s (%s);%s;Cluster Purity",
						     toString(score), title, times, title),
						(int)((x_max - x_min) / x_wid), x_min, x_max,
						(int)((p_max - p_min) / p_wid), p_min, p_max);

	eff_pass[score] = std::make_unique<TH1D>(Form("good_%s_%s", name.Data(), toString(score)),
						 Form("%s Eff. Events Count vs %s (%s);Entries;%s",
						      toString(score), title, times, title),
						 (int)((fold_max - x_min) / x_wid), x_min, fold_max);

	purity[score]->SetLineColor(colors[score]);

	this->params.emplace(score, FitParams(title, times, Form("%s_%s", name.Data(), toString(score)),
					      x_min, fold_val, fold_max, x_wid,colors[score]));

	algolegend->AddEntry(params.at(score).sigma_dist, toString(score));
      }

      eff_total = std::make_unique<TH1D>(Form("flat_%s", name.Data()),
					 Form("n Events vs %s (%s);Entries;%s",
					      title, times, title),
					 (int)((fold_max - x_min) / x_wid), x_min, fold_max);
    }

    // Explicitly delete copy constructor and copy assignment operator
    // This is crucial because unique_ptr is non-copyable.
    PlotObj(const PlotObj&) = delete;
    PlotObj& operator=(const PlotObj&) = delete;

    // Explicitly define move constructor and move assignment operator
    // These will be implicitly generated by the compiler if all members are movable
    // and no copy operations are user-defined/deleted.
    // However, explicitly defining them makes it clear and robust.
    PlotObj(PlotObj&& other) noexcept = default;
    PlotObj& operator=(PlotObj&& other) noexcept = default;


    // Destructor (optional with unique_ptr, but good practice for clarity or if
    // you have non-unique_ptr cleanup)
    ~PlotObj() {
        // unique_ptr members will be automatically deleted.
        // If there were any raw pointers that were NOT owned by unique_ptr,
        // you would delete them here.
    }


    inline void set_param_maxes() {
      double max_val = -10, min_val = 1e60;
      for (auto key: fitparam_vec) {
	for (auto& [score, fitparam]: params) {
	  TH1D *hist = fitparam.from_enum(key);
	  hist->GetXaxis()->SetNdivisions(510);
	  double this_max = hist->GetMaximum();
	  if (this_max > max_val)
	    max_val = this_max;

	  double this_min = hist->GetMinimum();
	  if (this_min < min_val)
	    min_val = this_min;
	}

	if (debug) std::cout << "Max val " << toString(key) << ": " << max_val << std::endl;
	if (debug) std::cout << "Min val " << toString(key) << ": " << min_val << std::endl;
      
	for (auto& [score, fitparam]: params) {
	  double pad = 0.5*(1-0.75) * (max_val-min_val);
	  fitparam.from_enum(key)->GetYaxis()->SetRangeUser(min_val-pad, max_val+pad);
	}
	// reset values
	max_val = -10;
	min_val = 1e60;
      }
    }

    inline void plot_postprocessing() {
      for (ScoreType score: enum_vec) {
	TH2D* big_hist = hist.at(score).get();
	for(int j = 0; j < big_hist->GetNbinsX(); ++j) {
	  auto color = colors[j % colors.size()]; // Assuming colors is std::vector or similar, and size() is valid
	  auto name = Form("%s_slice%d",big_hist->GetName(),j);

	  int right = big_hist->GetXaxis()->GetBinLowEdge(j+1) >= fold_value
	    ? big_hist->GetNbinsX()+1
	    : j+1;

	  std::unique_ptr<TH1D> hSlice = std::unique_ptr<TH1D>((TH1D*)big_hist->ProjectionY(name, j+1, right)->Clone());

	  if (hSlice->Integral() == 0 || hSlice->GetEntries() < 6) {
	    continue;
	  }

	  hSlice->SetLineColor(color);
	  hSlice->SetLineWidth(2);

	  double leftEdge  = big_hist->GetXaxis()->GetBinLowEdge(j+1);
	  double rightEdge = big_hist->GetXaxis()->GetBinLowEdge(right+1);

	  if (leftEdge < 1.0 && leftEdge > 0.0)
	    leftEdge *= 100;
	  if (rightEdge < 1.0 && rightEdge > 0.0)
	    rightEdge *= 100;

	  hSlice->SetTitle(Form("%s #in [%d,%d);%s;A.U.",
				big_hist->GetTitle(),
				(int)leftEdge,(int)rightEdge,
				this->ytitle));
      
	  std::unique_ptr<TF1> dgaus_fit = std::unique_ptr<TF1>(create_fit(hSlice.get(),true));
	  params.at(score).fill_each(j, dgaus_fit.get());
	  slices_hists[score].push_back(std::move(hSlice));
	  slices_fits[score].push_back(std::move(dgaus_fit));
	
	  if (right == big_hist->GetNbinsX()+1)
	    break;
	}

	std::unique_ptr<TEfficiency> eff = std::make_unique<TEfficiency>(*this->eff_pass.at(score), *this->eff_total);
	eff->SetStatisticOption(TEfficiency::kFNormal);  // use normal approximation
	eff->SetTitle(Form("Efficiency vs %s (%s);%s;Efficiency", this->xtitle, this->times, this->xtitle));
        eff->SetLineColor(myutl::colors[score]);
        eff->SetLineWidth(2);
	efficiency[score] = std::move(eff);
      }
      this->set_param_maxes();
    }
    
    inline void plot_logic(TCanvas* canvas) {
      // Draw 2Ds
      canvas->Print(Form("%s[",fname.Data()));
      for (const auto& [score,entry_ptr] : this->hist) {
	entry_ptr->GetXaxis()->SetRangeUser(entry_ptr->GetXaxis()->GetXmin(), this->fold_max);
	entry_ptr->Draw("COLZ");
	canvas->Print(fname);
      }

      // Draw Purities
      double purity_ymin = 0, purity_ymax = 1.4;
      bool first = true;
      // canvas->SetLogy(true);
      for (const auto& [score,entry] : purity) {
	TProfile* prof = entry->ProfileX();
	prof->SetMarkerColor(colors[score]);
	prof->SetLineColor(colors[score]);
	prof->SetLineWidth(2);
	prof->GetXaxis()->SetRangeUser(x_min, fold_max);
	prof->GetXaxis()->SetNdivisions(510);  // 5 major, 1 minor subdivision
	prof->GetYaxis()->SetRangeUser(purity_ymin,purity_ymax);
	if (first) {
	  prof->SetTitle(Form("Avg. Cluster Purity vs. %s (%s);%s;%s",
			      entry->GetXaxis()->GetTitle(),this->times,
			      entry->GetXaxis()->GetTitle(),"Purity"));
	  prof->Draw("E1");  // 'E' to draw error bars
	  first = false;
	} else
	  prof->Draw("E1 SAME");  // 'E' to draw error bars
      }
      algolegend->Draw("SAME");
      canvas->Print(fname);
      // canvas->SetLogy(false);
      
      // Draw Efficiencies
      first = true;
      double eff_ymin = 0, eff_ymax = 1.4;
      for (auto score : enum_vec) {
	auto eff = efficiency.at(score).get();
	if (first) {
	  eff->Draw("AP");
	  gPad->Update();  // Ensure painting
	  eff->GetPaintedGraph()->GetYaxis()->SetRangeUser(eff_ymin, eff_ymax);
	  eff->GetPaintedGraph()->GetXaxis()->SetLimits(x_min, fold_max);
	  eff->GetPaintedGraph()->GetXaxis()->SetRangeUser(x_min, fold_max);
	  eff->GetPaintedGraph()->GetXaxis()->SetNdivisions(510);
	  gPad->Update();
	  first = false;
	} else
	  eff->Draw("P SAME");
      }

      TLegend* efflegend = (TLegend*)algolegend->Clone();
      TLine *max_eff_line = new TLine(x_min, 0.99, fold_max, 0.99);
      max_eff_line->SetLineColor(kRed);
      max_eff_line->SetLineWidth(2);
      max_eff_line->SetLineStyle(4);
      efflegend->AddEntry(max_eff_line,"99% Efficiency");
      max_eff_line->Draw("SAME");
      efflegend->Draw("SAME");
      canvas->Print(fname);

      // Draw Resolution
      for (auto key: fitparam_vec) {
	first = true;
	for (auto& [score, param]: params) {
	  TH1D *hist = param.from_enum(key);
	  if (first) {
	    hist->Draw();
	    first = false;
	  } else
	    hist->Draw("SAME");
	}

	TLegend* resolegend = (TLegend*)algolegend->Clone();
	if (key == SIGMA) {
	  TLine *min_exp_res = new TLine(x_min, 15, fold_max, 15);
	  min_exp_res->SetLineColor(kRed);
	  min_exp_res->SetLineWidth(2);   
	  min_exp_res->SetLineStyle(4);

	  resolegend->AddEntry(min_exp_res,"15ps");
	  min_exp_res->Draw("SAME");
	}
	resolegend->Draw("SAME");
	canvas->Print(fname);
      }

      // draw slices o7
      TLatex latex;
      latex.SetTextSize(0.04);
      latex.SetTextAlign(13);  // Align left-top
      for (auto score: enum_vec) {
	const auto& fits = slices_fits[score];
	const auto& slices = slices_hists[score];
	if (slices.empty())
	  continue;

	for (int i=0; i < fits.size(); i++) {
	  auto fit = fits[i].get();
	  auto hSlice = slices[i].get();
	  std::unique_ptr<TLegend> thislegend = std::make_unique<TLegend>(0.65, 0.75, 0.9, 0.9);
	  auto res1text = Form("#sigma_{1}^{dgaus}=%.2f",fit->GetParameter(3));
	  auto res2text = Form("#sigma_{2}^{dgaus}=%.2f",fit->GetParameter(4));
	  auto chi2text = Form("#chi^{2}=%.2f",fit->GetChisquare()/fit->GetNDF());
	  thislegend->AddEntry(hSlice,"Histogram");
	  thislegend->AddEntry(fit,"Double Gaussian Fit");

	  hSlice->GetYaxis()->SetTitleOffset(1.4);
	  hSlice->Draw("HIST");
	  fit->Draw("SAME");
	  hSlice->GetXaxis()->SetRangeUser(-100, 100);
	  latex.DrawLatexNDC(0.18, 0.90, res1text);
	  latex.DrawLatexNDC(0.18, 0.85, res2text);
	  latex.DrawLatexNDC(0.18, 0.80, chi2text);
	  thislegend->Draw("SAME");
	  canvas->Print(fname); // slices
	
	  // with log scale
	  canvas->SetLogy(true);
	  hSlice->Draw("HIST");
	  fit->Draw("SAME");
	  hSlice->GetXaxis()->SetRangeUser(-400, 400);
	  latex.DrawLatexNDC(0.18, 0.90, res1text);
	  latex.DrawLatexNDC(0.18, 0.85, res2text);
	  latex.DrawLatexNDC(0.18, 0.80, chi2text);
	  thislegend->Draw("SAME");
	  canvas->Print(fname); // slices
	  canvas->SetLogy(false);
	}
      }
      canvas->Print(Form("%s]",fname.Data()));
    }
  };
  
  static double getSmearedTrackTime(int idx, double res, BranchPointerWrapper *branch) {
    // know that it has a valid time in HGTD and has a valid truth link

    int particle_idx = branch->track_to_particle[idx]; // this is anything
    int truthvtx_idx = branch->track_to_truthvtx[idx]; // this is anything
    if (debug) std::cout << "Particle_idx: " << particle_idx << std::endl;
    if (debug) std::cout << "TruthVtx_idx: " << truthvtx_idx << std::endl;
    double t_part;
    if (particle_idx == -1) { // no particle link
      if (truthvtx_idx != -1) { // some valid truth link, smear that vertex time
        if (debug) std::cout << "PATH A" << std::endl;
        t_part = branch->truth_vtx_time[truthvtx_idx];
      } else { // some random pileup, assume
        if (debug) std::cout << "PATH B" << std::endl;
        t_part = gRandom->Gaus(branch->truth_vtx_time[0],175);
      }
    } else {
      if (debug) std::cout << "PATH C" << std::endl;
      t_part = branch->particle_t[particle_idx];
    }
    if (debug) std::cout << "DONE" << std::endl;
    return gRandom->Gaus(t_part,res);
  }

  static bool passTrackVertexAssociation(int track_idx, int vertex_idx, BranchPointerWrapper *branch, double significance_cut) {
    double
      trk_z0    = branch->track_z0[track_idx],
      trk_z_var = branch->track_var_z0[track_idx],
      trk_pt    = branch->track_pt[track_idx],
      trk_eta   = branch->track_eta[track_idx],
      vx_z      = branch->reco_vtx_z[vertex_idx],
      vx_z_var  = 0.0;
    int near = (int)branch->track_near_idx[track_idx];
    if (near != vertex_idx and branch->track_near_sig[track_idx] < 1.0)
      return false;
    
    double nsigma = std::abs(trk_z0 - vx_z) / std::sqrt(trk_z_var + vx_z_var);
    return nsigma < significance_cut;
  }

  static double getDistanceBetweenClusters(Cluster a, Cluster b) {
    std::vector<double> a_v = a.values;
    std::vector<double> a_s = a.sigmas;
    std::vector<double> b_v = b.values;
    std::vector<double> b_s = b.sigmas;

    if (a_v.size() != b_v.size())
      std::cout << "Uh ohhh" << std::endl;

    std::vector<double> distances(a_v.size());
    for (int i=0; i < a_v.size(); i++) {
    
      double diff_i = (a_v.at(i) - b_v.at(i));
      double denom_i = sqrt(a_s.at(i)*a_s.at(i) + b_s.at(i)*b_s.at(i));
      distances.at(i) = diff_i/denom_i;
    }

    double dsqr = 0.0;
    for (double d : distances){
      dsqr += d*d;
    }
    return sqrt(dsqr);
  }

  static Cluster mergeClusters(Cluster a, Cluster b) {
    Cluster merged_cluster;
    merged_cluster.wasMerged = true;
  
    for (size_t i = 0; i < a.values.size(); i++) {
      double v1 = a.values.at(i);
      double v2 = b.values.at(i);
      double s1 = a.sigmas.at(i);
      double s2 = b.sigmas.at(i);

      double new_cluster_value =
        (v1 / (s1 * s1) + v2 / (s2 * s2)) / (1. / (s1 * s1) + 1. / (s2 * s2));
      double new_cluster_sigma =
	sqrt((s2 * s1) * (s2 * s1) / (s2 * s2 + s1 * s1));
      merged_cluster.values.push_back(new_cluster_value);
      merged_cluster.sigmas.push_back(new_cluster_sigma);
    }

    for (auto t: a.all_times)
      merged_cluster.all_times.push_back(t);
    for (auto t: b.all_times)
      merged_cluster.all_times.push_back(t);

    for (auto i: a.track_indices)
      merged_cluster.track_indices.push_back(i);
    for (auto i: b.track_indices)
      merged_cluster.track_indices.push_back(i);

    merged_cluster.nConstituents = a.nConstituents+b.nConstituents;
    merged_cluster.max_pt_cluster = a.max_pt_cluster || b.max_pt_cluster;

    for (ScoreType score: enum_vec)
      merged_cluster.scores[score] = a.scores[score]+b.scores[score];

    merged_cluster.scores[ScoreType::IDEAL] = merged_cluster.values.at(0);

    if (merged_cluster.max_pt_cluster)
      merged_cluster.scores[ScoreType::MAXPT] = 1e6;

    return merged_cluster;
  }

  static double calc_pt_purity(Cluster cluster, BranchPointerWrapper *branch) {
    double num = 0.0, denom = 0.0;
    for (auto trk: cluster.track_indices) {
      if (branch->track_to_truthvtx[trk] == 0) num += branch->track_pt[trk];
      denom += branch->track_pt[trk];
    }
    return num/denom;
  }

  static std::vector<int> getAssociatedTracks(BranchPointerWrapper *branch, double min_trk_pt) {
    std::vector<int> good_tracks;

    for (int trk = 0; trk < branch->track_z0.GetSize(); ++trk) {
      double
	trk_eta = branch->track_eta[trk],
	trk_pt  = branch->track_pt[trk],
	trk_quality = branch->track_quality[trk];

      if (min_hgtd_eta > std::abs(trk_eta) or std::abs(trk_eta) > max_hgtd_eta)
        continue;

      if (trk_pt < min_trk_pt)
	continue;

      if (not trk_quality)
	continue;

      if (passTrackVertexAssociation(trk, 0, branch, 3.0))
	good_tracks.push_back(trk);
    }

    return good_tracks;
  }

  static std::vector<Cluster> makeSimpleClusters(
    std::vector<int> track_indices,
    BranchPointerWrapper *branch,
    bool use_smeared_times,                         // if true, use smeared_times map
    const std::map<int, double>& smeared_times_map, // map of smeared times
    double fixed_res,                               // fixed resolution for smeared times
    bool check_time_valid                           // if true, apply branch->track_time_valid check
  ) {
    std::vector<Cluster> simpleClusters;
    for (auto idx: track_indices) {
      if (!check_time_valid || branch->track_time_valid[idx] == 1) {
	double trk_time;
	double trk_res; 
	if (use_smeared_times) {
	  trk_time = smeared_times_map.at(idx);
	  trk_res  = fixed_res;
	} else {
	  trk_time = branch->track_time[idx];   
	  trk_res  = branch->track_time_res[idx];
	}

	double trk_pt = branch->track_pt[idx];
	double z_significance = std::abs(branch->track_z0[idx] - branch->reco_vtx_z[0]) / std::sqrt(branch->track_var_z0[idx]);

	std::map<ScoreType,double> scores;
	scores[ScoreType::HGTD] = 0; // do not use this
	scores[ScoreType::SIG] = exp(-z_significance);
	scores[ScoreType::MAXPT] = 0;
	scores[ScoreType::TRKPT] = trk_pt;
	scores[ScoreType::SUMPT2] = trk_pt*trk_pt;
	scores[ScoreType::SIGTRKPT] = trk_pt*exp(-z_significance);
	scores[ScoreType::MAXHS] = branch->track_to_truthvtx[idx] == 0 ? 1 : 0;
	scores[ScoreType::JETPTDR] = branch->calc_jetpt_dr_score(idx);
	scores[ScoreType::TRKPTDR] = branch->calc_trkpt_dr_score(idx);
	scores[ScoreType::IDEAL] = trk_time;

	Cluster cluster = {{trk_time}, {trk_res}, {trk_time}, {idx}, scores};
	simpleClusters.push_back(cluster);
      }
    }
    return simpleClusters;
  }

  static void doSimultaneousClustering(std::vector<Cluster> *collection, double dist_cut) {
    double distance = 1.e30;
    while (collection->size() > 1) {
      // std::cout << "entering while loop" << std::endl;

      int i0 = 0;
      int j0 = 0;

      distance = getDistanceBetweenClusters(collection->at(0), collection->at(1));

      for (size_t i = 0; i < collection->size(); i++) {
	for (size_t j = i + 1; j < collection->size(); j++) {
	  Cluster a = collection->at(i);
	  Cluster b = collection->at(j);

	  if (a.wasMerged or b.wasMerged)
	    continue;
	
	  double current_distance =
	    getDistanceBetweenClusters(a, b);
	  if (current_distance <= distance) {
	    distance = current_distance;
	    i0 = i;
	    j0 = j;
	  }
	} // j loop
      } // i loop

      // fuse closest two vertices
      if (distance < dist_cut and i0 != j0) {
	Cluster new_cluster = mergeClusters(collection->at(i0),collection->at(j0));
	collection->erase(collection->begin()+j0);
	if (i0 < j0)
	  collection->erase(collection->begin()+i0);
	else
	  collection->erase(collection->begin()+(i0-1));

	collection->push_back(new_cluster);
      } else {
	if (std::find_if(collection->begin(), collection->end(),
			 [](Cluster a) {return a.wasMerged;}) != collection->end()) {
	  for (int idx=0; idx < collection->size(); ++idx) {
	    collection->at(idx).wasMerged = false;
	  }
	
	} else {
	  break;
	}
      }
      // break;
    } // While
  }

  static std::vector<Cluster> doConeClustering(
    std::vector<int> input_indices,             
    BranchPointerWrapper *branch,
    float nsigma_cut_val,                           // Size of window in t to use
    bool use_smeared_times,                         // if true, use smeared_times map
    const std::map<int, double>& smeared_times_map, // map of smeared times
    double fixed_res,                               // fixed resolution for smeared times
    bool check_time_valid                           // if true, apply branch->track_time_valid check
  ) {
    if (debug) std::cout << "Starting n-sigma cone clustering with " << input_indices.size() << " tracks" << std::endl;
	
    std::vector<Cluster> clusters;
    if (input_indices.empty()) {
      if (debug) std::cout << "No tracks for cone clustering" << std::endl;
      return clusters;
    }
    
    std::vector<bool> processed_track_indices(branch->track_pt.GetSize(), false);

    // Step 0: Filter input indices based on time validity if requested
    std::vector<int> timefiltered_indices;
    for (auto idx : input_indices) {
      if (!check_time_valid || branch->track_time_valid[idx] == 1) { // Apply check only if requested
	timefiltered_indices.push_back(idx);
      }
    }
    
    while (true) {
      // Step 1: Find the highest-pt unprocessed track from the time-filtered list
      int max_pt_seed_idx = -1;
      double max_pt_val = -1.0;

      for (int idx : timefiltered_indices) {
	if (processed_track_indices[idx]) continue; // Skip if already processed

	double this_pt = branch->track_pt[idx];
	if (this_pt > max_pt_val) {
	  max_pt_val = this_pt;
	  max_pt_seed_idx = idx;
	}
      }

      if (max_pt_seed_idx == -1) break; // all tracks have been checked

      // Step 2: Build cluster around max_pt_seed_idx
      processed_track_indices[max_pt_seed_idx] = true; // Mark seed as processed

      // Determine seed time and resolution based on `use_smeared_times`
      double seed_time;
      double seed_time_res;

      if (use_smeared_times) {
	seed_time = smeared_times_map.at(max_pt_seed_idx);
	seed_time_res = fixed_res; // Use fixed_res for smeared times
      } else {
	seed_time = branch->track_time[max_pt_seed_idx];
	seed_time_res = branch->track_time_res[max_pt_seed_idx];
      }

      std::vector<int> tracks_to_merge_into_cluster = {max_pt_seed_idx};

      for (int candidate_idx : timefiltered_indices) {
	if (processed_track_indices[candidate_idx]) continue; // Skip if already processed

	// Determine candidate time and resolution based on `use_smeared_times`
	double candidate_time;
	double candidate_time_res;

	if (use_smeared_times) {
	  candidate_time = smeared_times_map.at(candidate_idx);
	  candidate_time_res = fixed_res; // Use fixed_res for smeared times
	} else {
	  candidate_time = branch->track_time[candidate_idx];
	  candidate_time_res = branch->track_time_res[candidate_idx];
	}

	double dt_diff_abs = std::abs(seed_time - candidate_time);
	double dt_res_combined = std::sqrt(seed_time_res*seed_time_res + candidate_time_res*candidate_time_res);
	double nsigma = dt_diff_abs / dt_res_combined;

	if (nsigma < nsigma_cut_val) { // Use the passed nsigma_cut_val
	  tracks_to_merge_into_cluster.push_back(candidate_idx);
	  processed_track_indices[candidate_idx] = true; // Mark as processed
	}
      }
      
      std::vector<Cluster> initial_simple_clusters_for_seed =
	makeSimpleClusters(tracks_to_merge_into_cluster, branch,
			   use_smeared_times, smeared_times_map,
			   fixed_res, check_time_valid);

      Cluster final_seed_cluster;
      final_seed_cluster = initial_simple_clusters_for_seed[0];
      for (size_t k = 1; k < initial_simple_clusters_for_seed.size(); ++k) {
	final_seed_cluster = mergeClusters(final_seed_cluster, initial_simple_clusters_for_seed[k]);
      }
      clusters.push_back(final_seed_cluster);
      
    }
      
    for (Cluster& cluster: clusters)
      cluster.purity = calc_pt_purity(cluster, branch);
    
    return clusters;
  }

  static std::vector<Cluster> clusterTracksInTime(
     std::vector<int> track_indices,
     BranchPointerWrapper *branch,
     double distance_cut,    // Cone size/distance cut
     double smear_res,       // fixed resolution for smeared times
     bool use_smeared_times, // for ideal cases
     bool check_time_valid,  // for ideal efficiency
     bool isCone
  ) {
    if (track_indices.empty() && debug)
      std::cout << "EMPTY!!!!" << std::endl;

    std::map<int,double> smeared_times_map;

    // if use smeared times, generated
    if (use_smeared_times) {
      for (auto idx: track_indices) {
	if (!check_time_valid || branch->track_time_valid[idx] == 1) {
	  double smeared_time = getSmearedTrackTime(idx, smear_res, branch);
	  smeared_times_map.emplace(idx, smeared_time);
	}
      }
    }

    std::vector<Cluster> collection;
    
    if (not isCone) {
      collection = makeSimpleClusters(track_indices, branch, use_smeared_times, smeared_times_map, smear_res, check_time_valid);
      doSimultaneousClustering(&collection, distance_cut);
    } else {
      collection = doConeClustering(track_indices, branch, distance_cut, use_smeared_times, smeared_times_map, smear_res, check_time_valid);
    }

    for (Cluster& cluster: collection)
      cluster.purity = calc_pt_purity(cluster, branch);

    return collection;
  }
  
  static Cluster chooseHGTDCluster(std::vector<Cluster> collection, BranchPointerWrapper *branch) {
    double min_diff = 1e50, reco_time = branch->reco_vtx_time[0];
    Cluster min_cluster;
    for (Cluster cluster: collection) {
      double this_diff = std::abs(cluster.values.at(0) - reco_time);
      if (this_diff < min_diff) {
	min_diff = this_diff;
	min_cluster = cluster;
      }
    }
    return min_cluster;
  }

  static Cluster chooseBestHGTDCluster(std::vector<Cluster> collection, BranchPointerWrapper *branch) {
    double min_diff = 1e50, truth_time = branch->truth_vtx_time[0];
    Cluster min_cluster;
    for (Cluster cluster: collection) {
      double this_diff = std::abs(cluster.values.at(0) - truth_time);
      if (this_diff < min_diff) {
	min_diff = this_diff;
	min_cluster = cluster;
      }
    }
    return min_cluster;
  }
  
  static std::map<ScoreType, Cluster> chooseCluster(std::vector<Cluster> collection, BranchPointerWrapper *branch) {
    std::map<ScoreType,Cluster> output;
    // std::cout << "GAY LOSER DETECTED " << collection.size()<< "\n"; // 71354
  
    for (ScoreType score: enum_vec) {
      if (score == ScoreType::HGTD) {
	continue;
      }
      
      output[score] = collection[0]; // final time we are giving to the user
      double max_score = score != ScoreType::IDEAL
	? output[score].scores.at(score)
	: std::abs(output[score].scores.at(score) - branch->truth_vtx_time[0]);
    
      for (Cluster& cluster: collection) {
	double comp_score = score != ScoreType::IDEAL
	  ? cluster.scores[score]
	  : std::abs(cluster.scores[score] - branch->truth_vtx_time[0]);

	bool query = score != ScoreType::IDEAL
	  ? comp_score > max_score
	  : comp_score < max_score;

	if (query) {
	  max_score = comp_score;
	  output[score] = cluster;
	}
      }
    }
  
    return output;
  }

  static bool passEfficiency(Cluster cluster, BranchPointerWrapper *branch) {
    double diff = std::abs(cluster.values.at(0)-branch->truth_vtx_time[0]);
    if (diff > 60)
      return false;

    int nHSTrack = 0;
    for (auto idx: cluster.track_indices) {
      if (branch->track_to_truthvtx[idx] == 0)
	nHSTrack++;
    }

    return nHSTrack > 0;
  }

  static void plot_inclusive(const char* fname, bool logscale, bool fixbkg, double x_min, double x_max,
			     TCanvas *canvas, std::map<ScoreType,TH1D*> inclusive_resos) {
    TLatex latex;
    latex.SetTextSize(0.04);
    latex.SetTextAlign(13); 
    canvas->Print(Form("%s[",fname));
    canvas->SetLogy(logscale);
    for (auto pair: inclusive_resos) {
      TH1D *hist = pair.second;
      TF1* fit = create_fit(hist,fixbkg);
      TLegend* inclusive_legend = new TLegend(0.65, 0.75, 0.9, 0.9);

      inclusive_legend->AddEntry(hist,"Histogram");
      inclusive_legend->AddEntry(fit,"Double Gaussian Fit");
    
      hist->GetXaxis()->SetRangeUser(x_min, x_max);
      hist->Draw("HIST");
      fit->Draw("SAME");
      inclusive_legend->Draw("SAME");
    
      double dg_sigma = fit->GetParameter(3);
      latex.DrawLatexNDC(0.18, 0.90,Form("#sigma_{1}^{dgaus}=%.2f",dg_sigma));
      canvas->Print(fname);
    }
    canvas->Print(Form("%s]",fname));
  }

  static void plot_purity(const char* fname, bool logscale, TCanvas* canvas, TLegend* legend, std::map<ScoreType,TH1D*> purity_map) {
    canvas->Print(Form("%s[",fname));
    double maxval = 0.7;
    bool first = true;
    canvas->SetLogy(logscale);
    for (auto pair: purity_map) {
      TH1D *hist = pair.second;
      hist->SetLineColor(colors[pair.first % colors.size()]);
      hist->SetLineWidth(2);
      hist->SetBinContent(hist->GetNbinsX(), 
			  hist->GetBinContent(hist->GetNbinsX()) + hist->GetBinContent(hist->GetNbinsX()+1));
      hist->SetBinContent(hist->GetNbinsX()+1, 0); // clear overflow
      hist->Scale(1/hist->Integral());
      hist->SetMaximum(maxval);
      hist->SetMinimum(0.00);
      if (first) {
	hist->Draw("HIST");
	first = false;
      } else
	hist->Draw("HIST SAME");
    }
    canvas->SetLogy(false);
    legend->Draw("SAME");
    canvas->Print(Form("%s",fname));
    canvas->Print(Form("%s]",fname));
  }
  
  static void setup_chain(TChain &chain, const char* ntuple_dir) {
    // VBF H->Invisible sample
    for (const auto& entry : boost::filesystem::directory_iterator(ntuple_dir)) {
      if (entry.is_directory()) continue;
      if(debug) {std::cout << "Adding file: " << entry.path() << std::endl;}
      chain.Add(entry.path().c_str());
      // break;
    }
  
    if (chain.GetEntries() == 0) {
      std::cerr << "No ROOT files found in directory: " << std::endl;
      return;
    }
  }

  static void setup_chain(TChain &chain, std::string Number) {
    chain.Add(Form("../ntuple/user.scheong.42871997.Output._%s.SuperNtuple.root", Number.c_str()));
}
  
  template <typename T>
  static T folded(T raw, T fold) {
    return (raw >= fold) ? fold : raw;
  }

  static void process_event_data(
    BranchPointerWrapper *branch,
    bool use_smeared_times,
    bool check_valid_times,
    std::map<ScoreType,TH1D*>& inclusive_resos,
    std::map<ScoreType,TH1D*>& inclusive_purity,
    PlotObj& fjet,
    PlotObj& ftrack,
    PlotObj& pu_frac,
    PlotObj& hs_track,
    PlotObj& pu_track,
    PlotObj& recovtx_z,
    TH2D* hs_pu_inclusive
  ) {
    // check if vertex selection is correct & number of jets
    if (not branch->pass_basic_cuts()) return;

    // check if there is one forward jet with pt > 30 GeV
    if (not branch->pass_jet_pt_cut()) return;

    int nForwardJet=0;
    branch->count_forward_jets(nForwardJet);
    
    std::vector<int> tracks = getAssociatedTracks(branch, min_track_pt);

    int nForwardTrack=0, nForwardTrack_HS=0, nForwardTrack_PU=0;
    branch->count_forward_tracks(nForwardTrack,nForwardTrack_HS,nForwardTrack_PU,tracks);

    hs_pu_inclusive->Fill(nForwardTrack_HS,nForwardTrack_PU);
    
    double pu_ratio = (double)nForwardTrack_PU / (double)nForwardTrack;
    double reco_z = branch->reco_vtx_z[0];

    auto eff_fill_val_fjet     = folded(nForwardJet     , (int)fold_fjet) ;
    auto eff_fill_val_track    = folded(nForwardTrack   , (int)fold_track);
    auto eff_fill_val_hs_track = folded(nForwardTrack_HS, (int)fold_hs_track);
    auto eff_fill_val_pu_track = folded(nForwardTrack_PU, (int)fold_pu_track);
    auto eff_fill_val_z        = reco_z;
    auto eff_fill_val_pu_ratio = pu_ratio;

    if (debug) {
      std::cout << "nForwardJet = " << nForwardJet << std::endl;
      std::cout << "nForwardTrack = " << nForwardTrack << std::endl;
      std::cout << "nForwardTrack_HS = " << nForwardTrack_HS << std::endl;
      std::cout << "nForwardTrack_PU = " << nForwardTrack_PU << std::endl;
      std::cout << "Vertex_z = " << reco_z << std::endl;
    }
    
    std::vector<Cluster> cluster =
      clusterTracksInTime(tracks, branch, 3.0, 10.0, use_smeared_times, check_valid_times, true);

    fjet.eff_total->     Fill(eff_fill_val_fjet    );
    ftrack.eff_total->   Fill(eff_fill_val_track   );
    pu_frac.eff_total->  Fill(eff_fill_val_pu_ratio);
    hs_track.eff_total-> Fill(eff_fill_val_hs_track);
    pu_track.eff_total-> Fill(eff_fill_val_pu_track);
    recovtx_z.eff_total->Fill(eff_fill_val_z       );
    
    if (cluster.size() == 0) return;

    std::map<ScoreType,Cluster> chosen =
      chooseCluster(cluster, branch);

    // run HGTD Clustering (simultaneous)
    std::vector<Cluster> hgtd_clusters =
      clusterTracksInTime(tracks, branch, 3.0, -1, false, true, false);

    if (branch->reco_vtx_valid[0] == 1)
      chosen[ScoreType::HGTD] = chooseHGTDCluster(hgtd_clusters, branch);
    
    // bool eventquery =
    //   passEfficiency(chosen[ScoreType::IDEAL], branch) &&
    //   not passEfficiency(chosen[ScoreType::TRKPT], branch);

    // if (eventquery) {
    //   Long64_t event_num = branch->reader.GetTree()->GetReadEvent()-branch->reader.GetTree()->GetChainOffset();
    //   TString fileName   = (branch->reader.GetTree()->GetCurrentFile()->GetName());
    //   TString file =  fileName(40,6);
    //   std::cout <<
    // 	Form("python event_display_VBF_R25.py --file_num %s --event_num %lld --vtxID 0",
    // 	     file.Data(), event_num) << std::endl;
    // }
    
    for (ScoreType score : enum_vec) {
      
      if (branch->reco_vtx_valid[0] == 0 and score == HGTD)
	continue;
      
      Cluster scored = chosen[score];
      double score_based_time = scored.values.at(0);
      double cluster_purity = scored.purity;
      double diff = score_based_time - branch->truth_vtx_time[0];

      inclusive_resos.at(score)->Fill(diff);
      inclusive_purity.at(score)->Fill(cluster_purity);
      
      if (passEfficiency(scored, branch)) {
	fjet.eff_pass.at(score)->     Fill(eff_fill_val_fjet    );
	ftrack.eff_pass.at(score)->   Fill(eff_fill_val_track   );
	pu_frac.eff_pass.at(score)->  Fill(eff_fill_val_pu_ratio);
	hs_track.eff_pass.at(score)-> Fill(eff_fill_val_hs_track);
	pu_track.eff_pass.at(score)-> Fill(eff_fill_val_pu_track);
	recovtx_z.eff_pass.at(score)->Fill(eff_fill_val_z       );
      }
      
      // fill diff hists
      fjet.hist.at(score)->     Fill(nForwardJet     , diff);
      ftrack.hist.at(score)->   Fill(nForwardTrack   , diff);
      pu_frac.hist.at(score)->  Fill(pu_ratio        , diff);
      hs_track.hist.at(score)-> Fill(nForwardTrack_HS, diff);
      pu_track.hist.at(score)-> Fill(nForwardTrack_PU, diff);
      recovtx_z.hist.at(score)->Fill(reco_z          , diff);

      // fill purities
      fjet.purity.at(score)->     Fill(nForwardJet     , cluster_purity);
      ftrack.purity.at(score)->   Fill(nForwardTrack   , cluster_purity);
      pu_frac.purity.at(score)->  Fill(pu_ratio        , cluster_purity);
      hs_track.purity.at(score)-> Fill(nForwardTrack_HS, cluster_purity);
      pu_track.purity.at(score)-> Fill(nForwardTrack_PU, cluster_purity);
      recovtx_z.purity.at(score)->Fill(reco_z          , cluster_purity);
    }
  }
}
