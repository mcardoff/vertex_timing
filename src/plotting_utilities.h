#ifndef PLOTTING_UTILITIES_H
#define PLOTTING_UTILITIES_H

#include "clustering_includes.h"
#include "clustering_constants.h"
#include <TEfficiency.h>
#include <iostream>
#include <memory>

namespace MyUtl {

  TF1 *dgausFitFunc = new TF1("dgaus", "[1]*TMath::Exp(-0.5 * ((x - [0]) / [3])^2) + [2]*TMath::Exp(-0.5 * ((x - [0]) / [4])^2)");
  TF1 *sgausFitFunc = new TF1("sgaus", "[0]*TMath::Exp(-0.5 * ((x - [1]) / [2])^2)");

  auto createDblFit(
    TH1D* hist, bool fixbkg
  ) -> TF1* {
    dgausFitFunc->SetParNames("Mean", "Norm1", "Norm2", "Sigma1", "Sigma2");
    TF1* fit = new TF1("dgaus_fit", "dgaus", hist->GetXaxis()->GetXmin(), hist->GetXaxis()->GetXmax());
    fit->SetParameters(0, 0.8*hist->GetMaximum(), 1e-2*hist->GetMaximum(), 26.0, 175.0);
    fit->FixParameter(0, 0);
    fit->SetParLimits(1, 0, 1.E6); // can only be positive
    fit->SetParLimits(2, 0, 1.E6); // can only be positive
    fit->SetParLimits(3, 0.1, 1.E6); // can only be positive
    if (fixbkg) { fit->FixParameter(4, PILEUP_SMEAR); }
    else { fit->SetParLimits(4, 0, 1.E6); } // realistic values
    fit->SetLineColor(kRed);
    fit->SetLineWidth(2);
    fit->SetNpx(1000);
    hist->Fit(fit,"RQ");
    return fit;
  }

  auto createSngFit(
    TH1D* hist
  ) -> TF1* {
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
    TH1D* meanDist;
    TH1D* sigmaDist;
    TH1D* coreAmpDist;
    TH1D* backAmpDist;
    TH1D* ampRatioDist;

    FitParams(const char* title, const char* times, const char* name, 
	      double xMin, double foldVal, double foldMax,
	      double xWidth, int color) {
      auto fullTitle = Form("%%s vs %s (%s);%s;%%s", title, times, title);
      int nbins = (int)((foldMax - xMin) / xWidth);
      meanDist      = new TH1D(Form("mean_dist_%s", name),
			       Form(fullTitle, "Common #mu", "#mu"),
			       nbins, xMin, foldMax);

      sigmaDist     = new TH1D(Form("sigma_dist_%s", name),
			       Form(fullTitle, "Core #sigma", "Core #sigma"),
			       nbins, xMin, foldMax);

      coreAmpDist  = new TH1D(Form("amp1_dist_%s", name),
			      Form(fullTitle, "Core Amplitude", "Amplitude"),
			      nbins, xMin, foldMax);

      backAmpDist  = new TH1D(Form("amp2_dist_%s", name),
			      Form(fullTitle, "Bkg Amplitude", "Amplitude"),
			      nbins, xMin, foldMax);

      ampRatioDist = new TH1D(Form("ampratio_dist_%s", name),
			      Form(fullTitle, "Core Amp/Bkg Amp", "A.U."),
			      nbins, xMin, foldMax);

      for (const auto& hist: {meanDist, sigmaDist, coreAmpDist, backAmpDist, ampRatioDist}) {
	hist->SetLineWidth(2);
	hist->SetLineColor(color);
      }
    }

    void fillEach(int idx, TF1* fit) const {
      double sigma1 = fit->GetParameter("Sigma1");
      
      sigmaDist->SetBinContent(idx+1,sigma1);
      sigmaDist->SetBinError(idx+1,fit->GetParError("Sigma1"));
    
      double amp1     = fit->GetParameter("Norm1");
      double amp2     = fit->GetParameter("Norm2");
      double amp1Err  = fit->GetParError ("Norm1");
      double amp2Err  = fit->GetParError ("Norm2");
      double ampRatio = amp1/amp2;
      double ampRatioErr = std::sqrt(pow((1.0/amp2)*amp1Err,2)+pow(-amp2Err*ampRatio/amp2,2));

      coreAmpDist->SetBinContent(idx+1,amp1);   
      coreAmpDist->SetBinError(idx+1,amp1Err); 
    
      backAmpDist->SetBinContent(idx+1,amp2);   
      backAmpDist->SetBinError(idx+1,amp2Err); 

      ampRatioDist->SetBinContent(idx+1,ampRatio);
      ampRatioDist->SetBinError(idx+1,ampRatioErr); 
    }

    void fillGaus(int idx, TF1* fit) const {
      double sigma = fit->GetParameter(2);

      sigmaDist->SetBinContent(idx+1,sigma);
      sigmaDist->SetBinError(idx+1,fit->GetParError(2));

      double amp2    = fit->GetParameter(0);
      double amp2Err = fit->GetParError (0);
    
      backAmpDist->SetBinContent(idx+1,amp2);   
      backAmpDist->SetBinError(idx+1,amp2Err); 
    }

    auto fromEnum(FitParamFields fit) const -> TH1D* {
      switch (fit) {
      case FitParamFields::MEAN:  return meanDist;
      case FitParamFields::SIGMA: return sigmaDist;
      case FitParamFields::CORE:  return coreAmpDist;
      case FitParamFields::BKG:   return backAmpDist;
      case FitParamFields::RATIO: return ampRatioDist;
      default:                    return nullptr;
      }
    }
  };

  struct PlotObj {
    Score scoreToUse;
    double binWidth;
    double foldValue, foldMax;
    TString fname;
    const char* xtitle;
    const char* ytitle;
    const char* times;
    double xMin, xMax, xWid;
    double yMin, yMax, yWid;
    std::unique_ptr<FitParams> params;
    std::unique_ptr<TH1D> effTotal;
    std::unique_ptr<TH1D> effPass;
    std::unique_ptr<TEfficiency> efficiency;
    std::unique_ptr<TH2D> hist;
    std::unique_ptr<TH2D> purity;
    std::vector<std::unique_ptr<TH1D>> slicesHists;
    std::vector<std::unique_ptr<TF1>> slicesFits;

  PlotObj(const char* title, const char* times,
	  const char* fname, Score score,
	  double xMin, double xMax, double xWid,
	  double yMin, double yMax, double yWid,
	  double pMin, double pMax, double pWid,
	  double foldVal, double foldMax)
    : binWidth(xWid), fname(fname), scoreToUse(score),
      xtitle(title), ytitle("#Delta t[ps]"), times(times),
      xMin(xMin), xMax(xMax), xWid(xWid),
      yMin(yMin), yMax(yMax), yWid(yWid),
      foldValue(foldVal), foldMax(foldMax)
    {
      TString name(title);
      name.ReplaceAll(" ", "_");
      name.ReplaceAll("{", "");
      name.ReplaceAll("}", "");
      int xbins = (int)((xMax - xMin) / xWid);
      int ybins = (int)((yMax - yMin) / yWid);
      int pbins = (int)((pMax - pMin) / pWid);
      int fbins = (int)((foldMax - xMin) / xWid);
      hist = std::make_unique<TH2D>(
        Form("%s_%s_%s", name.Data(), toString(score),times),
	Form("%s t_{0} - TruthVtx t_{0} vs %s (%s);%s;#Delta t[ps]",
	toString(score), title, times, title),
	xbins, xMin, xMax, ybins, yMin, yMax);

      effPass = std::make_unique<TH1D>(
        Form("good_%s_%s_%s", name.Data(), toString(score), times),
	Form("%s Eff. Events Count vs %s (%s);Entries;%s", 
	     toString(score), title, times, title),
	fbins, xMin, foldMax);
      
      effTotal = std::make_unique<TH1D>(
        Form("flat_%s_%s_%s", name.Data(), toString(score), times),
	Form("%s Eff. Events Count vs %s (%s);Entries;%s", 
	     toString(score), title, times, title),
	fbins, xMin, foldMax);
      
      params = std::make_unique<FitParams>(
        title, times, Form("%s_%s_%s", name.Data(), toString(score),times),
	xMin, foldVal, foldMax, xWid,COLORS[score]);
      
      purity = std::make_unique<TH2D>(
        Form("cluster_purity_%s_%s_%s", name.Data(), toString(score),times),
	Form("%s Cluster Purity vs %s (%s);%s;Cluster Purity",
	     toString(score), title, times, title),
	xbins, xMin, xMax, pbins, pMin, pMax);
      purity->SetLineColor(COLORS[score]);
    }

    PlotObj(const PlotObj&) = delete;
    PlotObj& operator=(const PlotObj&) = delete;
    PlotObj(PlotObj&& other) noexcept = default;
    PlotObj& operator=(PlotObj&& other) noexcept = default;
    ~PlotObj() {}

    void fillPurity(double& x, double&y) { this->purity->Fill(x,y); }
    void fillPurity(int& x, double&y) { this->purity->Fill(x,y); }

    void fillDiff(double& x, double&y) { this->hist->Fill(x,y); }
    void fillDiff(int& x, double&y) { this->hist->Fill(x,y); }

    void fillPass(double& val) { this->effPass->Fill(val); }
    void fillPass(int& val) { this->effPass->Fill(val); }

    void fillTotal(double& val) { this->effTotal->Fill(val); }
    void fillTotal(int& val) { this->effTotal->Fill(val); }
    
    inline void setParamMaxes() {
      double maxVal = -10, minVal = 1e60;
      for (auto key: FITPARAM_VEC) {
	TH1D *hist = params->fromEnum(key);
	hist->GetXaxis()->SetNdivisions(510);
	double thisMax = hist->GetMaximum();
	if (thisMax > maxVal)
	  maxVal = thisMax;

	double thisMin = hist->GetMinimum();
	if (thisMin < minVal)
	  minVal = thisMin;

	if (key == FitParamFields::SIGMA)
	  maxVal = 40.0, minVal = 10;

	double pad = 0.5*(1-0.75) * (maxVal-minVal);
	params->fromEnum(key)->GetYaxis()->SetRangeUser(minVal-pad, maxVal+pad);
      }
    }

    inline void plotPostProcessing() {
      for(int j = 0; j < hist->GetNbinsX(); ++j) {
	auto color = COLORS[j % COLORS.size()];
	auto name = Form("%s_slice%d",hist->GetName(),j);

	int right = hist->GetXaxis()->GetBinLowEdge(j+1) >= foldValue
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

	std::unique_ptr<TF1> dgausFit = std::unique_ptr<TF1>(createDblFit(hSlice.get(),true));
	std::unique_ptr<TF1> sgausFit = std::unique_ptr<TF1>(createSngFit(hSlice.get()));

	double percErr = dgausFit->GetParError("Sigma1")/dgausFit->GetParameter("Sigma1");

	if (percErr < 0.10) {
	  params->fillEach(j, dgausFit.get());
	  slicesFits.push_back(std::move(dgausFit));
	} else {
	  params->fillGaus(j, sgausFit.get());
	  slicesFits.push_back(std::move(sgausFit));
	}
	slicesHists.push_back(std::move(hSlice));

	if (right == hist->GetNbinsX()+1)
	  break;
      }

      std::unique_ptr<TEfficiency> eff = std::make_unique<TEfficiency>(*this->effPass, *this->effTotal);
      eff->SetStatisticOption(TEfficiency::kFNormal);
      eff->SetTitle(Form("Efficiency vs %s (%s);%s;Efficiency", this->xtitle, this->times, this->xtitle));
      eff->SetLineColor(MyUtl::COLORS[this->scoreToUse]);
      eff->SetLineWidth(2);
      efficiency = std::move(eff);

      this->setParamMaxes();
    }

    inline void printEfficiencyStats() {
      std::cout << this->times << " Efficiency for dt vs. " << this->xtitle << " ScoreType: " << toString(this->scoreToUse) << '\n';
      for (int i=1; i<effPass->GetNbinsX(); i++) {
	double numPass = effPass->GetBinContent(i), numTotal = effTotal->GetBinContent(i);
	std::cout << toString(this->scoreToUse) << " "
		  << xtitle << " in [" << effPass->GetBinLowEdge(i)
		  << ", " << effPass->GetBinLowEdge(i+1) << ")" << '\n';
	std::cout << "Num. Passing: " << numPass << '\n';
	std::cout << "Num. Total : " << numTotal << '\n';
	std::cout << "Num. Failing: " << numTotal - numPass << '\n';
      }
      std::cout << "--------------------------------------------------\n";
      std::cout << "Total Num. Passing: " << effPass->Integral() << '\n';
      std::cout << "Total Num. Failing: " << effTotal->Integral() - effPass->Integral() << '\n';
      std::cout << "Overall Efficiency: " << effPass->Integral()/effTotal->Integral() << '\n';
      std::cout << "--------------------------------------------------\n";
    }
  };

  struct AnalysisObj {
    // PlotObjs which store histograms vs event-level variables
    std::map<std::string, std::unique_ptr<PlotObj>> dataObjects;
    // inclusive resolution for the scoretype provided
    std::unique_ptr<TH1D> inclusiveReso;
    std::unique_ptr<TH1D> inclusivePurity;
    
    // Non-const version for modification
    std::unique_ptr<PlotObj>& operator[](const std::string& key) {
      return dataObjects[key];
    }

    const std::unique_ptr<PlotObj>& operator[](const std::string& key) const {
      auto it = dataObjects.find(key);
      if (it == dataObjects.end())
        throw std::out_of_range("Key not found in AnalysisObj");
      return it->second;
    }

    const std::unique_ptr<PlotObj>& get(const std::string& key) const {
        return dataObjects.at(key);
    }
    
    AnalysisObj(
      const char* timetypeIDer, // HGTD Times v Ideal Res. HGTD
      Score score
    ) {
      TString filenameIDer(timetypeIDer);
      filenameIDer.ReplaceAll(" ", "");
      filenameIDer.ReplaceAll(".", "");
      filenameIDer.ToLower();
      dataObjects["fjet"] = std::make_unique<PlotObj>(
        "n Forward Jets", timetypeIDer,
	Form("figs/%s_nfjet.pdf",filenameIDer.Data()),
	score,
	FJET_MIN  , FJET_MAX  , FJET_WIDTH  ,
	DIFF_MIN  , DIFF_MAX  , DIFF_WIDTH  ,
	PURITY_MIN, PURITY_MAX, PURITY_WIDTH,
	FOLD_FJET, FOLD_FJET+FJET_WIDTH     );
  
      dataObjects["ftrack"] = std::make_unique<PlotObj>(
        "n Forward Tracks", timetypeIDer,
	Form("figs/%s_ntrack.pdf",filenameIDer.Data()),
	score,
	TRACK_MIN , TRACK_MAX , TRACK_WIDTH ,
	DIFF_MIN  , DIFF_MAX  , DIFF_WIDTH  ,
	PURITY_MIN, PURITY_MAX, PURITY_WIDTH,
	FOLD_TRACK, FOLD_TRACK+TRACK_WIDTH  );
  
      dataObjects["pu_frac"] = std::make_unique<PlotObj>(
        "Pile Up Fraction", timetypeIDer, 
        Form("figs/%s_pufrac.pdf",filenameIDer.Data()),
	score,
	PU_FRAC_MIN, PU_FRAC_MAX, PU_FRAC_WIDTH,
	DIFF_MIN   , DIFF_MAX   , DIFF_WIDTH   ,
	PURITY_MIN , PURITY_MAX , PURITY_WIDTH ,
	FOLD_PU_FRAC, PU_FRAC_MAX              );
  
      dataObjects["hs_track"] = std::make_unique<PlotObj>(
        "n Forward HS Tracks", timetypeIDer, 
	Form("figs/%s_nhstrack.pdf",filenameIDer.Data()),
	score,
	HS_TRACK_MIN, HS_TRACK_MAX, HS_TRACK_WIDTH ,
	DIFF_MIN    , DIFF_MAX    , DIFF_WIDTH     ,
	PURITY_MIN  , PURITY_MAX  , PURITY_WIDTH   ,
	FOLD_HS_TRACK, FOLD_HS_TRACK+HS_TRACK_WIDTH);
  
      dataObjects["pu_track"] = std::make_unique<PlotObj>(
        "n Forward PU Tracks", timetypeIDer, 
	Form("figs/%s_nputrack.pdf",filenameIDer.Data()),
	score,
	PU_TRACK_MIN, PU_TRACK_MAX, PU_TRACK_WIDTH ,
	DIFF_MIN    , DIFF_MAX    , DIFF_WIDTH     ,
	PURITY_MIN  , PURITY_MAX  , PURITY_WIDTH   ,
	FOLD_PU_TRACK, FOLD_PU_TRACK+PU_TRACK_WIDTH);

      inclusiveReso = std::make_unique<TH1D>(
        Form("%s_reso_%s",toString(score), filenameIDer.Data()),
	Form("Inclusive %s t_{0} - TruthVtx t_{0} (%s);#Delta t[ps];Entries",
	     toString(score), timetypeIDer),
	(int)((DIFF_MAX-DIFF_MIN)/DIFF_WIDTH), DIFF_MIN, DIFF_MAX);

      inclusivePurity = std::make_unique<TH1D>(
        Form("%s_purity_%s", toString(score), filenameIDer.Data()),
	Form("%s Purity (%s);Purity;Entries", toString(score), timetypeIDer),
	(int)((PURITY_MAX-PURITY_MIN)/PURITY_WIDTH), PURITY_MIN, PURITY_MAX);
    }

    inline void postProcessing() {
      for (auto &[str, plt]: this->dataObjects)
	plt->plotPostProcessing();
    }
  };

  void moneyPlot(
    const char* fname,
    const std::string& key,
    TCanvas* canvas,
    const std::vector<AnalysisObj*>& plts
  ) {
    if (plts.empty()) return;

    double xMin = plts[0]->get(key)->xMin;
    double xMax = plts[0]->get(key)->foldMax;

    auto plotWithLegend = [&](
        auto getter, const char* title,
	double yMin, double yMax, double
	refVal = -1, const char* refLabel = nullptr
    ) {
      TLegend* legend = new TLegend(0.55, 0.65, 0.9, 0.9);
      bool first = true;
      int counter = 0;

      for (const auto& ana : plts) {
	const auto& plt = ana->get(key);
	auto obj = getter(plt);
	obj->SetLineColor(COLORS[counter++ % COLORS.size()]);
	legend->AddEntry(obj, Form("%s: %s", toString(plt->scoreToUse), plt->times));

	if (first) {
	  obj->SetTitle(Form("%s vs %s", title, plt->xtitle));
	  obj->Draw("E1");
	  if constexpr (std::is_same_v<decltype(obj), TEfficiency*> || std::is_same_v<decltype(obj), std::shared_ptr<TEfficiency>>) {
	    std::cout << "DOING EFFICIENCY STUFF" << std::endl;
	    gPad->Update();
	    obj->GetPaintedGraph()->GetYaxis()->SetRangeUser(yMin, yMax);
	    obj->GetPaintedGraph()->GetXaxis()->SetRangeUser(xMin, xMax);
	  }
	  first = false;
	} else
	  obj->Draw("E1 SAME");
      }

      if (refVal != -1) {
	TLine* refLine = new TLine(xMin, refVal, xMax, refVal);
	refLine->SetLineColor(kRed);
	refLine->SetLineWidth(2);
	refLine->SetLineStyle(4);
	legend->AddEntry(refLine, refLabel);
	refLine->Draw("SAME");
      }

      legend->Draw("SAME");
      gPad->Update();
    };

    canvas->Print(Form("%s[", fname));

    // Plot Resolution
    plotWithLegend(
        [](auto& plt) { return plt->params->sigmaDist; },
	"Core #sigma", RES_YMIN, RES_YMAX,
	15.0, "15ps");
    canvas->Print(fname);

    // Plot Efficiency
    plotWithLegend(
        [](auto& plt) { return plt->efficiency.get(); },
	"Efficiency", EFF_YMIN, EFF_YMAX,
	0.99, "99% Efficiency");

    canvas->Print(fname);
    canvas->Print(Form("%s]", fname));
  }


  void inclusivePlot(
    const char* fname,
    bool logScale, bool fixBkg,
    double xMin, double xMax,
    TCanvas *canvas,
    const std::vector<AnalysisObj*>& plts
  ) {
    TLatex latex;
    latex.SetTextSize(0.04);
    latex.SetTextAlign(13); 
    canvas->Print(Form("%s[",fname));
    canvas->SetLogy(logScale);
    for (const auto& plt: plts) {
      TH1D* hist = (TH1D*)plt->inclusiveReso->Clone();
      TF1* fit1 = createDblFit(hist,fixBkg);
      TLegend* inclusiveLegend = new TLegend(0.65, 0.75, 0.9, 0.9);

      inclusiveLegend->AddEntry(hist,"Histogram");
      inclusiveLegend->AddEntry(fit1,"Double Gaussian Fit");
    
      hist->GetXaxis()->SetRangeUser(xMin, xMax);
      hist->Draw("HIST");
      fit1->Draw("SAME");
      inclusiveLegend->Draw("SAME");
    
      double dgSigma = fit1->GetParameter(3);
      latex.DrawLatexNDC(0.18, 0.90,Form("#sigma_{1}^{dgaus}=%.2f",dgSigma));
      canvas->Print(fname);
    }
    canvas->Print(Form("%s]",fname));
  }

  void purityPlot(
    const char* fname,
    bool logscale,
    TCanvas* canvas,
    TLegend* legend,
    const std::map<Score,TH1D*>& purityMap
  ) {
    canvas->Print(Form("%s[",fname));
    double maxval = 0.4;
    bool first = true;
    canvas->SetLogy(logscale);
    for (auto pair: purityMap) {
      TH1D *hist = pair.second;
      hist->SetLineColor(COLORS[pair.first % COLORS.size()]);
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
