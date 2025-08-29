#ifndef PLOTTING_UTILITIES_H
#define PLOTTING_UTILITIES_H

#include "clustering_includes.h"
#include "clustering_constants.h"
#include "clustering_structs.h"
#include <algorithm>
#include <memory>

namespace myutl {

  static TF1 *dgaus_fitfunc = new TF1("dgaus", "[1]*TMath::Exp(-0.5 * ((x - [0]) / [3])^2) + [2]*TMath::Exp(-0.5 * ((x - [0]) / [4])^2)");
  static TF1 *sgaus_fitfunc = new TF1("sgaus", "[0]*TMath::Exp(-0.5 * ((x - [1]) / [2])^2)");

  static TF1* create_dgaus_fit(
    TH1D *hist, bool fixbkg
  ) {
    dgaus_fitfunc->SetParNames("Mean", "Norm1", "Norm2", "Sigma1", "Sigma2");
    TF1* fit = new TF1("dgaus_fit", "dgaus", hist->GetXaxis()->GetXmin(), hist->GetXaxis()->GetXmax());
    fit->SetParameters(0, 0.8*hist->GetMaximum(), 1e-2*hist->GetMaximum(), 26.0, 175.0);
    fit->FixParameter(0, 0);
    fit->SetParLimits(1, 0, 1.E6); // can only be positive
    fit->SetParLimits(2, 0, 1.E6); // can only be positive
    fit->SetParLimits(3, 0.1, 1.E6); // can only be positive
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

  static TF1* create_sgaus_fit(
    TH1D *hist
  ) {
    TF1* fit = new TF1("gaus_fit", "sgaus", hist->GetXaxis()->GetXmin(), hist->GetXaxis()->GetXmax());
    fit->SetParameters(hist->GetMaximum(), 0, 30);
    fit->SetParLimits(0, 0, hist->GetMaximum()); // can only be positive
    fit->FixParameter(1, 0);
    fit->SetParLimits(2, 0, 1.E4); // can only be positive
    // fit->FixParameter(2, 175);
    fit->SetLineColor(kBlue);
    fit->SetLineWidth(2);
    fit->SetNpx(1000);
    hist->Fit(fit,"RQ");
    return fit;
  }
  
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
      int nbins = (int)((fold_max - x_min) / x_width);
      mean_dist      = new TH1D(Form("mean_dist_%s", name),
				Form(full_title, "Common #mu", "#mu"),
				nbins, x_min, fold_max);

      sigma_dist     = new TH1D(Form("sigma_dist_%s", name),
				Form(full_title, "Core #sigma", "Core #sigma"),
				nbins, x_min, fold_max);

      core_amp_dist  = new TH1D(Form("amp1_dist_%s", name),
				Form(full_title, "Core Amplitude", "Amplitude"),
				nbins, x_min, fold_max);

      back_amp_dist  = new TH1D(Form("amp2_dist_%s", name),
				Form(full_title, "Bkg Amplitude", "Amplitude"),
				nbins, x_min, fold_max);

      amp_ratio_dist = new TH1D(Form("ampratio_dist_%s", name),
				Form(full_title, "Core Amp/Bkg Amp", "A.U."),
				nbins, x_min, fold_max);

      for (const auto& hist: {mean_dist, sigma_dist, core_amp_dist, back_amp_dist, amp_ratio_dist}) {
	hist->SetLineWidth(2);
	hist->SetLineColor(color);
      }
    }

    void fillEach(int idx, TF1* fit) {
      double sigma1 = fit->GetParameter("Sigma1");
      
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

    void fillGaus(int idx, TF1* fit) {
      double sigma = fit->GetParameter(2);

      sigma_dist->SetBinContent(idx+1,sigma);
      sigma_dist->SetBinError(idx+1,fit->GetParError(2));

      double amp2     = fit->GetParameter(0);
      double amp2_err = fit->GetParError (0);
    
      back_amp_dist->SetBinContent(idx+1,amp2);   
      back_amp_dist->SetBinError(idx+1,amp2_err); 
    }

    TH1D* fromEnum(FitParamFields fit) {
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
    ScoreType score_to_use;
    double bin_width;
    double fold_value, fold_max;
    TString fname;
    const char* xtitle;
    const char* ytitle;
    const char* times;
    double x_min, x_max, x_wid;
    double y_min, y_max, y_wid;
    std::unique_ptr<FitParams> params;
    std::unique_ptr<TH1D> eff_total;
    std::unique_ptr<TH1D> eff_pass;
    std::unique_ptr<TEfficiency> efficiency;
    std::unique_ptr<TH2D> hist;
    std::unique_ptr<TH2D> purity;
    std::vector<std::unique_ptr<TH1D>> slices_hists;
    std::vector<std::unique_ptr<TF1>> slices_fits;

  PlotObj(const char* title, const char* times,
	  const char* fname, ScoreType score,
	  double x_min, double x_max, double x_wid,
	  double y_min, double y_max, double y_wid,
	  double p_min, double p_max, double p_wid,
	  double fold_val, double fold_max)
    : bin_width(x_wid), fname(fname), score_to_use(score),
      xtitle(title), ytitle("#Delta t[ps]"), times(times),
      x_min(x_min), x_max(x_max), x_wid(x_wid),
      y_min(y_min), y_max(y_max), y_wid(y_wid),
      fold_value(fold_val), fold_max(fold_max)
    {
      TString name(title);
      name.ReplaceAll(" ", "");
      int xbins = (int)((x_max - x_min) / x_wid);
      int ybins = (int)((y_max - y_min) / y_wid);
      int pbins = (int)((p_max - p_min) / p_wid);
      int fbins = (int)((fold_max - x_min) / x_wid);
      hist = std::make_unique<TH2D>(
        Form("%s_%s_%s", name.Data(), toString(score),times),
	Form("%s t_{0} - TruthVtx t_{0} vs %s (%s);%s;#Delta t[ps]",
	toString(score), title, times, title),
	xbins, x_min, x_max, ybins, y_min, y_max);

      eff_pass = std::make_unique<TH1D>(
        Form("good_%s_%s_%s", name.Data(), toString(score), times),
	Form("%s Eff. Events Count vs %s (%s);Entries;%s", 
	     toString(score), title, times, title),
	fbins, x_min, fold_max);
      
      eff_total = std::make_unique<TH1D>(
        Form("flat_%s_%s_%s", name.Data(), toString(score), times),
	Form("%s Eff. Events Count vs %s (%s);Entries;%s", 
	     toString(score), title, times, title),
	fbins, x_min, fold_max);
      
      params = std::make_unique<FitParams>(
        title, times, Form("%s_%s_%s", name.Data(), toString(score),times),
	x_min, fold_val, fold_max, x_wid,colors[score]);
      
      purity = std::make_unique<TH2D>(
        Form("cluster_purity_%s_%s_%s", name.Data(), toString(score),times),
	Form("%s Cluster Purity vs %s (%s);%s;Cluster Purity",
	     toString(score), title, times, title),
	xbins, x_min, x_max, pbins, p_min, p_max);
      purity->SetLineColor(colors[score]);
    }

    PlotObj(const PlotObj&) = delete;
    PlotObj& operator=(const PlotObj&) = delete;
    PlotObj(PlotObj&& other) noexcept = default;
    PlotObj& operator=(PlotObj&& other) noexcept = default;
    ~PlotObj() {}

    void FillPurity(double& x, double&y) { this->purity->Fill(x,y); }
    void FillPurity(int& x, double&y) { this->purity->Fill(x,y); }

    void FillDiff(double& x, double&y) { this->hist->Fill(x,y); }
    void FillDiff(int& x, double&y) { this->hist->Fill(x,y); }

    void FillPass(double& val) { this->eff_pass->Fill(val); }
    void FillPass(int& val) { this->eff_pass->Fill(val); }

    void FillTotal(double& val) { this->eff_total->Fill(val); }
    void FillTotal(int& val) { this->eff_total->Fill(val); }
    
    inline void setParamMaxes() {
      double max_val = -10, min_val = 1e60;
      for (auto key: fitparam_vec) {
	TH1D *hist = params->fromEnum(key);
	hist->GetXaxis()->SetNdivisions(510);
	double this_max = hist->GetMaximum();
	if (this_max > max_val)
	  max_val = this_max;

	double this_min = hist->GetMinimum();
	if (this_min < min_val)
	  min_val = this_min;

	if (key == FitParamFields::SIGMA)
	  max_val = 40.0, min_val = 10;

	double pad = 0.5*(1-0.75) * (max_val-min_val);
	params->fromEnum(key)->GetYaxis()->SetRangeUser(min_val-pad, max_val+pad);
      }
    }

    inline void plotPostProcessing() {
      for(int j = 0; j < hist->GetNbinsX(); ++j) {
	auto color = colors[j % colors.size()];
	auto name = Form("%s_slice%d",hist->GetName(),j);

	int right = hist->GetXaxis()->GetBinLowEdge(j+1) >= fold_value
	  ? hist->GetNbinsX()+1 : j+1;

	std::unique_ptr<TH1D> hSlice = std::unique_ptr<TH1D>((TH1D*)hist->ProjectionY(name, j+1, right)->Clone());

	if (hSlice->Integral() == 0 || hSlice->GetEntries() < 6) {
	  continue;
	}

	hSlice->SetLineColor(color);
	hSlice->SetLineWidth(2);

	double leftEdge  = hist->GetXaxis()->GetBinLowEdge(j+1);
	double rightEdge = hist->GetXaxis()->GetBinLowEdge(right+1);

	if (leftEdge < 1.0 && leftEdge > 0.0)
	  leftEdge *= 100;
	if (rightEdge < 1.0 && rightEdge > 0.0)
	  rightEdge *= 100;

	hSlice->SetTitle(Form("%s #in [%d,%d);%s;A.U.",
			      hist->GetTitle(),
			      (int)leftEdge,(int)rightEdge,
			      this->ytitle));

	std::unique_ptr<TF1> dgaus_fit = std::unique_ptr<TF1>(create_dgaus_fit(hSlice.get(),true));
	std::unique_ptr<TF1> sgaus_fit = std::unique_ptr<TF1>(create_sgaus_fit(hSlice.get()));

	double perc_err = dgaus_fit->GetParError("Sigma1")/dgaus_fit->GetParameter("Sigma1");

	if (perc_err < 0.10) {
	  params->fillEach(j, dgaus_fit.get());
	  slices_fits.push_back(std::move(dgaus_fit));
	} else {
	  params->fillGaus(j, sgaus_fit.get());
	  slices_fits.push_back(std::move(sgaus_fit));
	}
	slices_hists.push_back(std::move(hSlice));

	if (right == hist->GetNbinsX()+1)
	  break;
      }

      std::unique_ptr<TEfficiency> eff = std::make_unique<TEfficiency>(*this->eff_pass, *this->eff_total);
      eff->SetStatisticOption(TEfficiency::kFNormal);
      eff->SetTitle(Form("Efficiency vs %s (%s);%s;Efficiency", this->xtitle, this->times, this->xtitle));
      eff->SetLineColor(myutl::colors[score_to_use]);
      eff->SetLineWidth(2);
      efficiency = std::move(eff);

      this->setParamMaxes();
    }

    inline void printEfficiencyStats() {
      std::cout << this->times << " Efficiency for dt vs. " << this->xtitle << " ScoreType: " << toString(this->score_to_use) << std::endl;
      for (int i=1; i<eff_pass->GetNbinsX(); i++) {
	double num_pass = eff_pass->GetBinContent(i), num_total = eff_total->GetBinContent(i);
	std::cout << toString(this->score_to_use) << " "
		  << xtitle << " in [" << eff_pass->GetBinLowEdge(i)
		  << ", " << eff_pass->GetBinLowEdge(i+1) << ")" << std::endl;
	std::cout << "Num. Passing: " << num_pass << std::endl;
	std::cout << "Num. Total : " << num_total << std::endl;
	std::cout << "Num. Failing: " << num_total - num_pass << std::endl;
      }
      std::cout << "--------------------------------------------------" << std::endl;
      std::cout << "Total Num. Passing: " << eff_pass->Integral() << std::endl;
      std::cout << "Total Num. Failing: " << eff_total->Integral() - eff_pass->Integral() << std::endl;
      std::cout << "Overall Efficiency: " << eff_pass->Integral()/eff_total->Integral() << std::endl;
      std::cout << "--------------------------------------------------" << std::endl;
    }
  };

  struct AnalysisObj {
    // PlotObjs which store histograms vs event-level variables
    std::map<std::string, std::unique_ptr<PlotObj>> data_objects;
    // inclusive resolution for the scoretype provided
    std::unique_ptr<TH1D> inclusive_reso;
    std::unique_ptr<TH1D> inclusive_purity;
    
    // Non-const version for modification
    std::unique_ptr<PlotObj>& operator[](const std::string& key) {
      return data_objects[key];
    }

    const std::unique_ptr<PlotObj>& operator[](const std::string& key) const {
      auto it = data_objects.find(key);
      if (it == data_objects.end())
        throw std::out_of_range("Key not found in AnalysisObj");
      return it->second;
    }

    const std::unique_ptr<PlotObj>& get(const std::string& key) const {
        return data_objects.at(key);
    }
    
    AnalysisObj(
      const char* filename_ider, // hgtdtimes v idealres, 
      const char* timetype_ider, // HGTD Times v Ideal Res. HGTD
      ScoreType score
    ) {
      data_objects["fjet"] = std::make_unique<PlotObj>(
        "n Forward Jets", timetype_ider,
	Form("figs/%s_nfjet.pdf",filename_ider),
	score,
	fjet_min  , fjet_max  , fjet_width  ,
	diff_min  , diff_max  , diff_width  ,
	purity_min, purity_max, purity_width,
	fold_fjet, fold_fjet+fjet_width     );
  
      data_objects["ftrack"] = std::make_unique<PlotObj>(
        "n Forward Tracks", timetype_ider,
	Form("figs/%s_ntrack.pdf",filename_ider),
	score,
	track_min , track_max , track_width ,
	diff_min  , diff_max  , diff_width  ,
	purity_min, purity_max, purity_width,
	fold_track, fold_track+track_width  );
  
      data_objects["pu_frac"] = std::make_unique<PlotObj>(
        "Pile Up Fraction", timetype_ider, 
        Form("figs/%s_pufrac.pdf",filename_ider),
	score,
	pu_frac_min, pu_frac_max, pu_frac_width,
	diff_min   , diff_max   , diff_width   ,
	purity_min , purity_max , purity_width ,
	fold_pu_frac, pu_frac_max              );
  
      data_objects["hs_track"] = std::make_unique<PlotObj>(
        "n Forward HS Tracks", timetype_ider, 
	Form("figs/%s_nhstrack.pdf",filename_ider),
	score,
	hs_track_min, hs_track_max, hs_track_width ,
	diff_min    , diff_max    , diff_width     ,
	purity_min  , purity_max  , purity_width   ,
	fold_hs_track, fold_hs_track+hs_track_width);
  
      data_objects["pu_track"] = std::make_unique<PlotObj>(
        "n Forward PU Tracks", timetype_ider, 
	Form("figs/%s_nputrack.pdf",filename_ider),
	score,
	pu_track_min, pu_track_max, pu_track_width ,
	diff_min    , diff_max    , diff_width     ,
	purity_min  , purity_max  , purity_width   ,
	fold_pu_track, fold_pu_track+pu_track_width);

      inclusive_reso = std::make_unique<TH1D>(
        Form("%s_reso_%s",toString(score), filename_ider),
	Form("Inclusive %s t_{0} - TruthVtx t_{0} (%s);#Delta t[ps];Entries",
	     toString(score), timetype_ider),
	(int)((diff_max-diff_min)/diff_width), diff_min, diff_max);

      inclusive_purity = std::make_unique<TH1D>(
        Form("%s_purity_%s", toString(score), filename_ider),
	Form("%s Purity (%s);Purity;Entries", toString(score), timetype_ider),
	(int)((purity_max-purity_min)/purity_width), purity_min, purity_max);
    }

    inline void postProcessing() {
      for (auto &[str, plt]: this->data_objects)
	plt->plotPostProcessing();
    }
  };
  
  static void moneyPlot(
    const char* fname,
    const std::string key,
    TCanvas* canvas,
    const std::vector<AnalysisObj*>& plts // these are by score now
  ) {
    if (plts.size() == 0)
      return;

    // Plot Key Metrics for each of these PlotObjs on the same plot
    canvas->Print(Form("%s[",fname));
    double
      x_min = plts[0]->get(key)->x_min,
      fold_max = plts[0]->get(key)->fold_max;
    const char * xtitle = plts[0]->get(key)->xtitle;
    int counter = 0;
    
    bool first = true;
    double res_ymin = 0.0, res_ymax = 40;
    TLegend* reslegend = new TLegend(0.65, 0.65, 0.9, 0.9);
    for (const auto &ana: plts) {
      const auto& plt = ana->get(key);
      const auto& params = plt->params;
      auto res = params->sigma_dist;
      res->SetLineColor(colors[counter % colors.size()]);
      reslegend->AddEntry(res, Form("%s: %s", toString(plt->score_to_use), plt->times));
      counter++;
      if (first) {
	res->SetTitle(Form("Core #sigma vs %s",xtitle));
	res->GetYaxis()->SetRangeUser(res_ymin,res_ymax);
	res->Draw("E1");
	first = false;
      } else
	res->Draw("E1 SAME");
    }
    reslegend->Draw("SAME");
    canvas->Print(fname);

    // Plot Efficiency
    first = true;
    double eff_ymin = 0.0, eff_ymax = 1.5;
    counter = 0;
    TLegend* efflegend = new TLegend(0.65, 0.65, 0.9, 0.9);
    for (const auto& ana: plts) {
      const auto& plt = ana->get(key);
      const auto& efficiency = plt->efficiency;	
	efficiency->SetLineColor(colors[counter % colors.size()]);
	efflegend->AddEntry(efficiency.get(), Form("%s: %s", toString(plt->score_to_use), plt->times));
	counter++;
	if (first) {
	  efficiency->SetTitle(Form("Efficiency vs %s",xtitle));
	  efficiency->Draw("AP");
	  gPad->Update();  // Ensure painting
	  efficiency->GetPaintedGraph()->GetYaxis()->SetRangeUser(eff_ymin, eff_ymax);
	  efficiency->GetPaintedGraph()->GetXaxis()->SetLimits(x_min, fold_max);
	  efficiency->GetPaintedGraph()->GetXaxis()->SetRangeUser(x_min, fold_max);
	  efficiency->GetPaintedGraph()->GetXaxis()->SetNdivisions(510);
	  auto* total = (TH1D*)plt->eff_total->Clone();
	  total->SetLineColor(colors[8]);
	  total->SetLineWidth(2);
	  total->Scale(1.0 / total->Integral());
	  total->Scale((0.2 - eff_ymin) / total->GetMaximum());
	  total->Draw("HIST SAME"); 
	  gPad->Update();
	  first = false;
	} else
	  efficiency->Draw("P SAME");
    }

    TLine *max_eff_line = new TLine(x_min, 0.99, fold_max, 0.99);
    max_eff_line->SetLineColor(kRed);
    max_eff_line->SetLineWidth(2);
    max_eff_line->SetLineStyle(4);
    efflegend->AddEntry(max_eff_line,"99% Efficiency");
    max_eff_line->Draw("SAME");
    efflegend->Draw("SAME");
    canvas->Print(fname);

    // Plot Purity
    first = true;
    counter = 0;
    double pur_ymin = 0.0, pur_ymax = 1.5;
    TLegend* purlegend = new TLegend(0.65, 0.65, 0.9, 0.9);
    for (const auto& ana: plts) {
      const auto& plt = ana->get(key);
      const auto& purity = plt->purity;
      auto pur = purity->ProfileX();
      pur->SetLineWidth(2);
      pur->GetXaxis()->SetRangeUser(x_min, fold_max);
      pur->GetYaxis()->SetRangeUser(pur_ymin, pur_ymax);
      pur->SetLineColor(colors[counter % colors.size()]);
      purlegend->AddEntry(pur, Form("%s: %s", toString(plt->score_to_use), plt->times));
      counter++;
      if (first) {
	pur->SetTitle(Form("Avg. Cluster Purity vs. %s;%s;%s",
			   purity->GetXaxis()->GetTitle(),
			   purity->GetXaxis()->GetTitle(),"Purity"));
	pur->Draw("E1");  // 'E' to draw error bars
	first = false;
      } else
	pur->Draw("E1 SAME");
    }
    TLine *max_pur_line = new TLine(x_min, 1.0, fold_max, 1.0);
    max_pur_line->SetLineColor(kRed);
    max_pur_line->SetLineWidth(2);
    max_pur_line->SetLineStyle(4);
    purlegend->AddEntry(max_pur_line,"100% Purity");
    max_pur_line->Draw("SAME");
    purlegend->Draw("SAME");
    canvas->Print(fname);
    canvas->Print(Form("%s]",fname));
  }

  static void inclusivePlot(
    const char* fname,
    bool logscale, bool fixbkg,
    double x_min, double x_max,
    TCanvas *canvas,
    std::map<ScoreType,TH1D*> inclusive_resos
  ) {
    TLatex latex;
    latex.SetTextSize(0.04);
    latex.SetTextAlign(13); 
    canvas->Print(Form("%s[",fname));
    canvas->SetLogy(logscale);
    for (auto pair: inclusive_resos) {
      TH1D *hist = pair.second;
      TF1* fit1 = create_dgaus_fit(hist,fixbkg);
      TLegend* inclusive_legend = new TLegend(0.65, 0.75, 0.9, 0.9);

      inclusive_legend->AddEntry(hist,"Histogram");
      inclusive_legend->AddEntry(fit1,"Double Gaussian Fit");
    
      hist->GetXaxis()->SetRangeUser(x_min, x_max);
      hist->Draw("HIST");
      fit1->Draw("SAME");
      inclusive_legend->Draw("SAME");
    
      double dg_sigma = fit1->GetParameter(3);
      latex.DrawLatexNDC(0.18, 0.90,Form("#sigma_{1}^{dgaus}=%.2f",dg_sigma));
      canvas->Print(fname);
    }
    canvas->Print(Form("%s]",fname));
  }

  static void purityPlot(
    const char* fname,
    bool logscale,
    TCanvas* canvas,
    TLegend* legend,
    std::map<ScoreType,TH1D*> purity_map
  ) {
    canvas->Print(Form("%s[",fname));
    double maxval = 0.4;
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
}
#endif // PLOTTING_UTILITIES_H
