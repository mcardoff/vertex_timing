#include <TCanvas.h>
#include <TColor.h>
#include <TChain.h>
#include <TEfficiency.h>
#include <TF1.h>
#include <TGraphAsymmErrors.h>
#include <TH1.h>
#include <TH2.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TLine.h>
#include <TProfile.h>
#include <TTreeReaderArray.h>
#include <memory>

#include <boost/filesystem.hpp>

namespace myutl {
  static bool debug = false;
  static int c01 = TColor::GetColor("#e41a1c"); 
  static int c02 = TColor::GetColor("#377eb8"); 
  static int c03 = TColor::GetColor("#4daf4a"); 
  static int c04 = TColor::GetColor("#984ea3");
  static int c05 = TColor::GetColor("#ff7f00"); 
  static int c06 = TColor::GetColor("#ffff33"); 
  static int c07 = TColor::GetColor("#a65628"); 
  static int c08 = TColor::GetColor("#f781bf");
  static int c09 = TColor::GetColor("#999999"); 
  static int c10 = TColor::GetColor("#92dadd");
  static std::vector<int> colors = {c01, c02, c03, c04, c05, c06, c07, c08, c09, c10};

  // cut variables
  static int    min_jets          =  1;    // min number of jets
  static double min_jetpt         =  30.0; // self explanatory
  static double min_abs_eta_jet   =  2.0;  // min eta for a jet to be considered "forward"
  static double min_abs_eta_track =  2.4;  // min eta for a track to be considered "forward"
  static double min_track_pt      =  1.0;  // track_pt > 1.0 GeV
  static double max_vtx_dz        =  2.0;  // max error for reco HS vertex z
  static double max_nsigma        =  3.0;  // how close a track can be to PV

  
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
  enum ScoreType { HGTD = 0, SIG = 1, SUMPT2 = 3, TRKPT = 2, SIGTRKPT = 4, MAXPT = 5, IDEAL = 6, INVALID = -99 };
  // std::vector<ScoreType> enum_vec = {HGTD,  SIG, SUMPT2, TRKPT, SIGTRKPT, IDEAL}; // valid values
  static std::vector<ScoreType> enum_vec = {HGTD, TRKPT, IDEAL}; // valid values

  static const char* toString(myutl::ScoreType score) {
    switch (score) {
    case ScoreType::HGTD:     return "HGTD Algorithm";
    case ScoreType::SIG:      return "exp(-|s|)";
    case ScoreType::SUMPT2:   return "Sum(p_{T})^{2}";
    case ScoreType::TRKPT:    return "Track p_{T}";
    case ScoreType::SIGTRKPT: return "Track p_{T}*exp(-|s|)";
    case ScoreType::MAXPT:    return "Cluster with max p_{T}";
    case ScoreType::IDEAL:    return "Ideal cluster choice";
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
    TTreeReaderArray<float> track_var_z0;
    TTreeReaderArray<float> track_time;
    TTreeReaderArray<float> track_time_res;
    TTreeReaderArray<int>   track_time_valid;
    TTreeReaderArray<int>   track_to_truthvtx;
    TTreeReaderArray<int>   track_to_particle;
    // TTreeReaderArray<int>   track_hgtd_hits;
    // TTreeReaderArray<int>   track_prim_hits;

    TTreeReaderArray<float> truth_vtx_z;
    TTreeReaderArray<float> truth_vtx_time;
    TTreeReaderArray<bool>  truth_vtx_ishs;

    TTreeReaderArray<float> reco_vtx_z;
    TTreeReaderArray<float> reco_vtx_time;
    TTreeReaderArray<float> reco_vtx_timeRes;
    TTreeReaderArray<int>   reco_vtx_valid;

    TTreeReaderArray<float> topojet_pt;
    TTreeReaderArray<float> topojet_eta;

    TTreeReaderArray<float> truthhsjet_pt;
    TTreeReaderArray<float> truthhsjet_eta;

    TTreeReaderArray<float> particle_t;

  BranchPointerWrapper(TTreeReader& r)
    : reader           (r                          ),
      weight           (r, "weight"                ),
      track_z0         (r, "Track_z0"              ),
      track_pt         (r, "Track_pt"              ),
      track_eta        (r, "Track_eta"             ),
      track_var_z0     (r, "Track_var_z0"          ),
      track_time       (r, "Track_time"            ),
      track_time_res   (r, "Track_timeRes"         ),
      track_time_valid (r, "Track_hasValidTime"    ),
      track_to_truthvtx(r, "Track_truthVtx_idx"    ),
      track_to_particle(r, "Track_truthPart_idx"   ),
      // track_hgtd_hits  (r, "Track_nHGTDHits"       ),
      // track_prim_hits  (r, "Track_nHGTDPrimaryHits"),
      truth_vtx_z      (r, "TruthVtx_z"            ),
      truth_vtx_time   (r, "TruthVtx_time"         ),
      truth_vtx_ishs   (r, "TruthVtx_isHS"         ),
      reco_vtx_z       (r, "RecoVtx_z"             ),
      reco_vtx_time    (r, "RecoVtx_time"          ),
      reco_vtx_timeRes (r, "RecoVtx_timeRes"       ),
      reco_vtx_valid   (r, "RecoVtx_hasValidTime"  ),
      topojet_pt       (r, "AntiKt4EMTopoJets_pt"  ),
      topojet_eta      (r, "AntiKt4EMTopoJets_eta" ),
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
	return false;;
      }
      
      return true;
    }

    bool pass_jet_pt_cut() { 
      int passptcount = 0;
      for(int jet_idx = 0; jet_idx < truthhsjet_eta.GetSize(); ++jet_idx) {
	float eta = truthhsjet_eta[jet_idx], pt = truthhsjet_pt[jet_idx];
	passptcount += pt > min_jetpt ? 1 : 0;
	if (passptcount >= 2)
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
	(nForwardJet)++;
    }
  }
    
    void count_forward_tracks(int& nFTrack, int& nFTrack_HS, int& nFTrack_PU,
				   std::vector<int> tracks) {
      nFTrack = 0; nFTrack_HS = 0; nFTrack_PU = 0;
      for(auto trk_idx: tracks) {
	auto eta = this->track_eta[trk_idx];
	auto valid = this->track_time_valid[trk_idx];
	auto pt = this->track_pt[trk_idx];
	auto truthvtx = this->track_to_truthvtx[trk_idx];
	// already know these pass association
	if (std::abs(eta) > min_abs_eta_track) {
	  nFTrack++;
	  if (truthvtx != -1 and this->truth_vtx_ishs[truthvtx])
	    nFTrack_HS++;
	  else
	    nFTrack_PU++;
	}
      }
    }
  };

  struct FitParams {
    TH1D* mean_dist;
    TH1D* sigma_dist;
    TH1D* core_amp_dist;
    TH1D* back_amp_dist;
    TH1D* amp_ratio_dist;

    FitParams(const char* title, const char* name,
	      double x_min, double fold_val, double fold_max,
	      double x_width, int color) {
      auto full_title = Form("%%s vs %s;%s;%%s", title, title);
      mean_dist = new TH1D(Form("mean_dist_%s", name),
			   Form(full_title, "Common #mu", "#mu"),
			   (int)((fold_max-x_min)/x_width), x_min, fold_max);
    
      sigma_dist = new TH1D(Form("sigma_dist_%s", name),
			    Form(full_title, "Core #sigma", "Core #sigma"),
			    (int)((fold_max-x_min)/x_width), x_min, fold_max);
    
      core_amp_dist = new TH1D(Form("amp1_dist_%s", name),
			       Form(full_title, "Core Amplitude", "Amplitude"),
			       (int)((fold_max-x_min)/x_width), x_min, fold_max);
    
      back_amp_dist = new TH1D(Form("amp2_dist_%s", name),
			       Form(full_title, "Bkg Amplitude", "Amplitude"),
			       (int)((fold_max-x_min)/x_width), x_min, fold_max);
    
      amp_ratio_dist = new TH1D(Form("ampratio_dist_%s", name),
				Form(full_title, "Core Amp/Bkg Amp", "A.U."),
				(int)((fold_max-x_min)/x_width), x_min, fold_max);

      sigma_dist->SetMaximum(50);
      sigma_dist->SetMinimum(5);

      // mean_dist->SetMaximum( 1.0);
      // mean_dist->SetMinimum(-1.0);
    
      for (auto hist: {mean_dist, sigma_dist, core_amp_dist, back_amp_dist, amp_ratio_dist}) {
	hist->SetLineWidth(2);
	hist->SetLineColor(color);
      }
    }

    void fill_each(int idx, TF1* fit) {
      sigma_dist->SetBinContent(idx+1,fit->GetParameter("Sigma1"));   // maybe should be i *shrug*
      sigma_dist->SetBinError(idx+1,fit->GetParError("Sigma1")); // maybe should be i *shrug*
    
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
    const char* fname;
    const char* xtitle;
    const char* ytitle;
    double x_min, x_max, x_wid;
    double y_min, y_max, y_wid;
    std::unique_ptr<TH1D> eff_total;
    std::map<ScoreType, std::unique_ptr<TH1D>> eff_pass;
    std::map<ScoreType, std::unique_ptr<TH2D>> hist;
    std::map<ScoreType, std::unique_ptr<TH2D>> purity;
    std::map<ScoreType, FitParams> params; // FitParams might handle its own pointers, or be values
    std::map<ScoreType, std::vector<std::unique_ptr<TH1D>>> slices_hists;
    std::map<ScoreType, std::vector<std::unique_ptr<TF1>>> slices_fits;
    std::map<ScoreType, std::unique_ptr<TEfficiency>> efficiency; // This needs to be carefully managed if created from eff_pass and eff_total
    std::unique_ptr<TLegend> algolegend;

  PlotObj(const char* title, const char* name, const char* fname,
	  double x_min, double x_max, double x_wid,
	  double y_min, double y_max, double y_wid,
	  double p_min, double p_max, double p_wid,
	  double fold_val, double fold_max)
    : bin_width(x_wid), fname(fname),
      xtitle(title), ytitle("#Delta t[ps]"),
      fold_value(fold_val), fold_max(fold_max),
      x_min(x_min), x_max(x_max), x_wid(x_wid),
      y_min(y_min), y_max(y_max), y_wid(y_wid)
    {// Allocate with make_unique for owned objects
      algolegend = std::make_unique<TLegend>(0.6, 0.7, 0.9, 0.9);
      algolegend->SetTextSize(0.03);

      for (ScoreType score : enum_vec) {
	hist[score] = std::make_unique<TH2D>(Form("%s_%s", name, toString(score)),
					     Form("%s t_{0} - TruthVtx t_{0} vs %s ;%s;#Delta t[ps]", toString(score), title, title),
					     (int)((x_max - x_min) / x_wid), x_min, x_max,
					     (int)((y_max - y_min) / y_wid), y_min, y_max);

	purity[score] = std::make_unique<TH2D>(Form("cluster_purity_%s_%s", name, toString(score)),
					       Form("%s Cluster Purity vs %s ;%s;Cluster Purity", toString(score), title, title),
					       (int)((x_max - x_min) / x_wid), x_min, x_max,
					       (int)((p_max - p_min) / p_wid), p_min, p_max);

	eff_pass[score] = std::make_unique<TH1D>(Form("good_%s_%s", name, toString(score)),
						 Form("%s Eff Events ;Entries;%s", title, title),
						 (int)((fold_max - x_min) / x_wid), x_min, fold_max);

	purity[score]->SetLineColor(colors[score]);

	this->params.emplace(score,
			     FitParams(title, Form("%s_%s", name, toString(score)),
				       x_min, fold_val, fold_max, x_wid,
				       colors[score]));
        // slices_hists[score] = {};
        // slices_fits[score] = {}; 

	algolegend->AddEntry(params.at(score).sigma_dist, toString(score));
      }

      eff_total = std::make_unique<TH1D>(Form("flat_%s", name),
					 Form("%s All Events ;Entries;%s", title, title),
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
	eff->SetTitle(Form("Efficiency vs %s;%s;Efficiency", this->xtitle, this->xtitle));
        eff->SetLineColor(myutl::colors[score]);
        eff->SetLineWidth(2);
	efficiency[score] = std::move(eff);
      }
      this->set_param_maxes();
    }
    
    inline void plot_logic(TCanvas* canvas) {
      // Draw 2Ds
      canvas->Print(Form("%s[",fname));
      for (const auto& [score,entry_ptr] : this->hist) {
	entry_ptr->GetXaxis()->SetRangeUser(entry_ptr->GetXaxis()->GetXmin(), this->fold_max);
	entry_ptr->Draw("COLZ");
	canvas->Print(fname);
      }

      // Draw Purities
      double purity_ymin = 0, purity_ymax = 1.4;
      bool first = true;
      for (const auto& [score,entry] : purity) {
	TProfile* prof = entry->ProfileX();
	prof->SetMarkerColor(colors[score]);
	prof->SetLineColor(colors[score]);
	prof->SetLineWidth(2);
	prof->GetXaxis()->SetRangeUser(x_min, fold_max);
	prof->GetXaxis()->SetNdivisions(510);  // 5 major, 1 minor subdivision
	prof->GetYaxis()->SetRangeUser(purity_ymin,purity_ymax);
	if (first) {
	  prof->SetTitle(Form("Avg. Cluster Purity vs. %s",entry->GetXaxis()->GetTitle()));
	  prof->Draw("E1");  // 'E' to draw error bars
	  first = false;
	} else
	  prof->Draw("E1 SAME");  // 'E' to draw error bars
      }
      algolegend->Draw("SAME");
      canvas->Print(fname);

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
	  thislegend->AddEntry(hSlice,"Histogram");
	  thislegend->AddEntry(fit,"Double Gaussian Fit");

	  hSlice->GetYaxis()->SetTitleOffset(1.4);
	  hSlice->Draw("HIST");
	  fit->Draw("SAME");
	  hSlice->GetXaxis()->SetRangeUser(-100, 100);
	  latex.DrawLatexNDC(0.18, 0.90, res1text);
	  latex.DrawLatexNDC(0.18, 0.85, res2text);
	  thislegend->Draw("SAME");
	  canvas->Print(fname); // slices
	
	  // with log scale
	  canvas->SetLogy(true);
	  hSlice->Draw("HIST");
	  fit->Draw("SAME");
	  hSlice->GetXaxis()->SetRangeUser(-400, 400);
	  latex.DrawLatexNDC(0.18, 0.90, res1text);
	  latex.DrawLatexNDC(0.18, 0.85, res2text);
	  thislegend->Draw("SAME");
	  canvas->Print(fname); // slices
	  canvas->SetLogy(false);
	}
      }
      canvas->Print(Form("%s]",fname));
    }
  };

  static bool passTrackVertexAssociation(int track_idx, int vertex_idx, BranchPointerWrapper *branch, double min_trk_pt, double significance_cut) {
    double
      trk_z0    = branch->track_z0[track_idx],
      trk_z_var = branch->track_var_z0[track_idx],
      trk_pt    = branch->track_pt[track_idx],
      trk_eta   = branch->track_eta[track_idx],
      vx_z      = branch->reco_vtx_z[vertex_idx],
      vx_z_var  = 0.0;
  
    if (trk_pt < min_trk_pt)
      return false;

    // if (std::abs(trk_eta) < 2.4 || std::abs(trk_eta) > 4.0)
    //   return false;

    if (branch->track_time_valid[track_idx] == 0)
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

    for (ScoreType score: enum_vec) {
      merged_cluster.scores[score] = a.scores[score]+b.scores[score];
    }

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
	trk_pt  = branch->track_pt[trk];
      if (std::abs(trk_eta) < 2.4 || std::abs(trk_eta) > 4.0)
	continue;

      if (trk_pt < 1.0)
	continue;

      if (passTrackVertexAssociation(trk, 0, branch, min_trk_pt, 3.0))
	good_tracks.push_back(trk);
    }

    return good_tracks;
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

  static bool passEfficiency(double diff, Cluster cluster, BranchPointerWrapper *branch) {
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

  template <typename T>
  static T folded(T raw, T fold) {
    return (raw >= fold) ? fold : raw;
  }
}
