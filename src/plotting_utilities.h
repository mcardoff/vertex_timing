#ifndef PLOTTING_UTILITIES_H
#define PLOTTING_UTILITIES_H

// ---------------------------------------------------------------------------
// plotting_utilities.h
//   All histogram booking, filling, fitting, and plotting logic.  Three
//   structs and three free functions are defined here, all inside MyUtl.
//
//   Global fit functions (module-level TF1 objects):
//     tgausFitFunc — triple Gaussian (3 cores, shared mean)
//     dgausFitFunc — double Gaussian (core + background, shared mean)
//     sgausFitFunc — single Gaussian
//
//   Fit factory functions:
//     createTrpFit — configure and fit tgausFitFunc to a histogram
//     createDblFit — configure and fit dgausFitFunc to a histogram
//     createSngFit — configure and fit sgausFitFunc to a histogram
//
//   Structs:
//     FitParams   — collection of per-bin resolution and amplitude histograms
//                   extracted from double-Gaussian fits to timing residuals
//     PlotObj     — one (variable, score) analysis cell: owns the 2D timing
//                   residual histogram, efficiency TH1Ds, purity 2D histogram,
//                   per-bin slice histograms and their fits, and the
//                   TEfficiency object; drives all per-cell plotting
//     AnalysisObj — collection of PlotObjs keyed by variable name, plus the
//                   inclusive timing residual and purity histograms; top-level
//                   container for one (scenario, score) combination
//
//   Free plot functions:
//     moneyPlot     — multi-algorithm comparison: core σ, background σ,
//                     efficiency vs one kinematic variable
//     inclusivePlot — inclusive timing residual plot with double-Gaussian fit
//     purityPlot    — normalised cluster purity distributions
// ---------------------------------------------------------------------------

#include "clustering_includes.h"
#include "clustering_constants.h"
#include <TEfficiency.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <iomanip>
#include <iostream>
#include <memory>
#include <sstream>

namespace MyUtl {

  // ---------------------------------------------------------------------------
  // Global fit function objects
  //   Allocated once at program start and reused across all fits.
  //   tgausFitFunc — triple Gaussian: shared mean [0], three amplitudes [1-3],
  //                  three widths [4-6].  Width [6] is fixed to PILEUP_SMEAR.
  //   dgausFitFunc — double Gaussian: shared mean [0], two amplitudes [1-2],
  //                  two widths [3-4].  Width [4] optionally fixed to PILEUP_SMEAR.
  //   sgausFitFunc — single Gaussian: amplitude [0], mean [1] (fixed to 0),
  //                  width [2].
  // ---------------------------------------------------------------------------
  TF1 *tgausFitFunc = new TF1("tgaus", "[1]*TMath::Exp(-0.5 * ((x - [0]) / [4])^2) + "
			               "[2]*TMath::Exp(-0.5 * ((x - [0]) / [5])^2) + "
			               "[3]*TMath::Exp(-0.5 * ((x - [0]) / [6])^2)"  );
  TF1 *dgausFitFunc = new TF1("dgaus", "[1]*TMath::Exp(-0.5 * ((x - [0]) / [3])^2) + "
                                       "[2]*TMath::Exp(-0.5 * ((x - [0]) / [4])^2)");
  TF1 *sgausFitFunc = new TF1("sgaus", "[0]*TMath::Exp(-0.5 * ((x - [1]) / [2])^2)");

  // ---------------------------------------------------------------------------
  // createTrpFit
  //   Configures tgausFitFunc with sensible starting parameters relative to
  //   the histogram maximum, fixes the mean to 0 and the third width to
  //   PILEUP_SMEAR, constrains all amplitudes and widths to be positive, then
  //   fits the function to hist over its full range.  Returns the fitted TF1*.
  //   The caller owns the returned pointer (it is managed by ROOT's memory).
  // ---------------------------------------------------------------------------
  auto createTrpFit(
    TH1D* hist
  ) -> TF1* {
    tgausFitFunc->SetParNames("Mean", "Norm1", "Norm2", "Norm3", "Sigma1", "Sigma2", "Sigma3");
    TF1* fit = new TF1("tgaus_fit", "tgaus", hist->GetXaxis()->GetXmin(), hist->GetXaxis()->GetXmax());
    fit->SetParameters(0,
		       0.80*hist->GetMaximum(),
		       0.40*hist->GetMaximum(),
		       0.01*hist->GetMaximum(),
		       13.0, 26.0, 175.0);
    fit->FixParameter(0, 0);
    fit->FixParameter(6, PILEUP_SMEAR);
    fit->SetParLimits(1, 0, 1.E6); // can only be positive
    fit->SetParLimits(2, 0, 1.E6); // can only be positive
    fit->SetParLimits(3, 0, 1.E6); // can only be positive
    fit->SetParLimits(4, 0.1, PILEUP_SMEAR); // can only be positive
    fit->SetParLimits(5, 0.1, PILEUP_SMEAR); // can only be positive
    fit->SetLineColor(kRed);
    fit->SetLineWidth(2);
    fit->SetNpx(1000);
    hist->Fit(fit,"RQ");
    return fit;
  }

  // ---------------------------------------------------------------------------
  // createDblFit
  //   Configures dgausFitFunc for a double-Gaussian fit.  When fixbkg is true
  //   the background width is fixed to PILEUP_SMEAR (models a flat pileup
  //   component); when false it is free within [1, PILEUP_SMEAR].  Returns
  //   the fitted TF1*.  This is the primary fit used in plotPostProcessing
  //   for per-bin timing residual slices and in inclusivePlot.
  // ---------------------------------------------------------------------------
  auto createDblFit(
    TH1D* hist, bool fixbkg
  ) -> TF1* {
    dgausFitFunc->SetParNames("Mean", "Norm1", "Norm2", "Sigma1", "Sigma2");
    // TF1* fit = new TF1("dgaus_fit", "dgaus", hist->GetXaxis()->GetXmin(), hist->GetXaxis()->GetXmax());
    TF1* fit = new TF1("dgaus_fit", "dgaus", hist->GetXaxis()->GetXmin(), hist->GetXaxis()->GetXmax());
    fit->SetParameters(0, 0.8*hist->GetMaximum(), 1e-2*hist->GetMaximum(), 26.0, 175.0);
    fit->FixParameter(0, 0);
    fit->SetParLimits(1, 0, 1.E6); // can only be positive
    fit->SetParLimits(2, 0, 1.E6); // can only be positive
    fit->SetParLimits(3, 0.1, 1.E6); // can only be positive
    if (fixbkg) { fit->FixParameter(4, PILEUP_SMEAR); }
    else { fit->SetParLimits(4, 1.0, PILEUP_SMEAR); } // realistic values
    fit->SetLineColor(kRed);
    fit->SetLineWidth(2);
    fit->SetNpx(1000);
    hist->Fit(fit,"RQ");
    return fit;
  }

  // ---------------------------------------------------------------------------
  // createSngFit
  //   Configures sgausFitFunc for a single-Gaussian fit with the mean fixed
  //   to 0.  Used as a fallback when the double-Gaussian fit is poorly
  //   constrained (e.g. very few entries in a slice).  Returns the fitted TF1*.
  // ---------------------------------------------------------------------------
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
  
  // ---------------------------------------------------------------------------
  // FitParams
  //   Holds six TH1D* distributions (one per FitParamFields enum value)
  //   that collect fitted double-Gaussian parameters as a function of the
  //   plot's x-axis variable (e.g. n forward jets).  Each bin of these
  //   histograms corresponds to one timing-residual slice that was fitted
  //   with createDblFit.
  //
  //   Member methods:
  //     FitParams() constructor — books all six histograms with consistent
  //                               binning and colour
  //     fillEach               — extract core/background σ and amplitudes
  //                               from a double-Gaussian fit result (5-param)
  //     fillGaus               — extract σ and amplitude from a single-
  //                               Gaussian fit result (3-param fallback)
  //     fromEnum               — return the TH1D* for a given FitParamFields
  //                               key (used when drawing resolution plots)
  // ---------------------------------------------------------------------------
  struct FitParams {
    TH1D* meanDist;
    TH1D* sigmaDist;
    TH1D* bkgSigmaDist;
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

      bkgSigmaDist  = new TH1D(Form("bkg_sigma_dist_%s", name),
			       Form(fullTitle, "Background #sigma", "Background #sigma"),
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

      for (const auto& hist: {meanDist, bkgSigmaDist, sigmaDist, coreAmpDist, backAmpDist, ampRatioDist}) {
	hist->SetLineWidth(2);
	hist->SetLineColor(color);
      }
    }

    // -----------------------------------------------------------------------
    // fillEach
    //   Fills the core σ, background σ, core amplitude, background
    //   amplitude, and amplitude-ratio histograms for bin idx from a
    //   double-Gaussian (5-parameter) fit.  The smaller of Sigma1/Sigma2
    //   is treated as the core and the larger as the background.
    // -----------------------------------------------------------------------
    void fillEach(int idx, TF1* fit) const {
      double sigma1 = fit->GetParameter("Sigma1");
      double sigma2 = fit->GetParameter("Sigma2");

      if (sigma1 > sigma2) {
	sigmaDist->SetBinContent(idx+1,sigma2);
	sigmaDist->SetBinError(idx+1,fit->GetParError("Sigma2"));
	
	bkgSigmaDist->SetBinContent(idx+1,sigma1);
	bkgSigmaDist->SetBinError(idx+1,fit->GetParError("Sigma1"));
      } else {
	sigmaDist->SetBinContent(idx+1,sigma1);
	sigmaDist->SetBinError(idx+1,fit->GetParError("Sigma1"));
	
	bkgSigmaDist->SetBinContent(idx+1,sigma2);
	bkgSigmaDist->SetBinError(idx+1,fit->GetParError("Sigma2"));
      }
      
    
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

    // -----------------------------------------------------------------------
    // fillGaus
    //   Fallback for single-Gaussian fits: fills only the sigma and
    //   background-amplitude histograms (the core σ and amplitude are not
    //   meaningful for a single Gaussian).
    // -----------------------------------------------------------------------
    void fillGaus(int idx, TF1* fit) const {
      double sigma = fit->GetParameter("Sigma1");

      sigmaDist->SetBinContent(idx+1,sigma);
      sigmaDist->SetBinError(idx+1,fit->GetParError("Sigma1"));

      double amp2    = fit->GetParameter("Norm1");
      double amp2Err = fit->GetParError ("Norm1");
    
      backAmpDist->SetBinContent(idx+1,amp2);   
      backAmpDist->SetBinError(idx+1,amp2Err); 
    }

    // -----------------------------------------------------------------------
    // fromEnum
    //   Returns the TH1D* corresponding to a FitParamFields enumerator.
    //   Used in plotLogic to iterate over fit parameter distributions in
    //   the canonical FITPARAM_VEC order.
    // -----------------------------------------------------------------------
    auto fromEnum(FitParamFields fit) const -> TH1D* {
      switch (fit) {
      case FitParamFields::MEAN:   return meanDist;
      case FitParamFields::SIGMA:  return sigmaDist;
      case FitParamFields::BSIGMA: return bkgSigmaDist;
      case FitParamFields::CORE:   return coreAmpDist;
      case FitParamFields::BKG:    return backAmpDist;
      case FitParamFields::RATIO:  return ampRatioDist;
      default:                     return nullptr;
      }
    }
  };

  // ---------------------------------------------------------------------------
  // PlotObj
  //   One (kinematic variable, Score) analysis cell.  Owns:
  //     hist        — 2D timing residual histogram (x = kinematic var, y = Δt)
  //     effPass     — 1D histogram counting events passing the efficiency test
  //     effTotal    — 1D histogram counting all events (denominator)
  //     efficiency  — TEfficiency built from effPass / effTotal after the loop
  //     purity      — 2D cluster-purity histogram (x = kinematic var, y = purity)
  //     params      — FitParams collecting per-bin fit results
  //     slicesHists — per-bin 1D timing residual projections
  //     slicesFits  — double-Gaussian fits to each slice
  //
  //   The foldValue / foldMax mechanism collapses high-multiplicity overflow
  //   into the last visible bin so that sparse tails do not dominate plots.
  //
  //   Member methods:
  //     PlotObj constructor  — books all histograms with consistent naming
  //     fillPurity           — fill the purity 2D histogram
  //     fillDiff             — fill the timing residual 2D histogram
  //     fillPass / fillTotal — fill efficiency numerator / denominator
  //     plotPostProcessing   — slice the 2D hist, fit each slice, build
  //                            TEfficiency; called after the event loop
  //     plotLogic            — draw and print all pages to a multi-page PDF
  //     printEfficiencyStats — print a formatted per-bin efficiency table
  // ---------------------------------------------------------------------------
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
	xMin, foldVal, foldMax, xWid,COLORS[score % COLORS.size()]);
      
      purity = std::make_unique<TH2D>(
        Form("cluster_purity_%s_%s_%s", name.Data(), toString(score),times),
	Form("%s Cluster Purity vs %s (%s);%s;Cluster Purity",
	     toString(score), title, times, title),
	xbins, xMin, xMax, pbins, pMin, pMax);
      purity->SetLineColor(COLORS[score % COLORS.size()]);
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
    
    // -----------------------------------------------------------------------
    // plotPostProcessing
    //   Post-event-loop processing: projects the 2D residual histogram into
    //   per-bin slices (collapsing bins beyond foldValue), fits each slice
    //   with createDblFit, stores the results in FitParams, and constructs
    //   the TEfficiency object.  Must be called before plotLogic.
    // -----------------------------------------------------------------------
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
	
	std::unique_ptr<TF1> dgausFit = std::unique_ptr<TF1>(createDblFit(hSlice.get(),true)); // varied bg
	// std::unique_ptr<TF1> dgausFit = std::unique_ptr<TF1>(createTrpFit(hSlice.get())); // TRIPLE GAUS :O
	std::unique_ptr<TF1> sgausFit = std::unique_ptr<TF1>(createSngFit(hSlice.get()));

	// double percErr = dgausFit->GetParError("Sigma1")/dgausFit->GetParameter("Sigma1");

	// if (percErr < 0.10) {
	params->fillEach(j, dgausFit.get());
	slicesFits.push_back(std::move(dgausFit));
	// } else {
	  // params->fillGaus(j, sgausFit.get());
	  // slicesFits.push_back(std::move(sgausFit));
	// }
	slicesHists.push_back(std::move(hSlice));

	if (right == hist->GetNbinsX()+1)
	  break;
      }

      std::unique_ptr<TEfficiency> eff = std::make_unique<TEfficiency>(*this->effPass, *this->effTotal);
      eff->SetStatisticOption(TEfficiency::kFNormal);
      eff->SetTitle(Form("Efficiency vs %s (%s);%s;Efficiency", this->xtitle, this->times, this->xtitle));
      eff->SetLineColor(MyUtl::COLORS[this->scoreToUse % COLORS.size()]);
      eff->SetLineWidth(2);
      efficiency = std::move(eff);

    }

    // -----------------------------------------------------------------------
    // plotLogic
    //   Renders and prints all pages for this PlotObj to a multi-page PDF
    //   (fname).  Page order:
    //     1. 2D timing residual (COLZ)
    //     2. TEfficiency vs x with 99% reference line
    //     3–N. Per-FitParamField distributions (σ_core, σ_bkg, amplitudes…)
    //     N+. Per-bin timing residual slices with fit overlay, linear and
    //         log scale
    //   Requires plotPostProcessing() to have been called first.
    // -----------------------------------------------------------------------
    inline void plotLogic(TCanvas *canvas) {
      canvas->Print(Form("%s[",fname.Data()));
      TH2D* hist = (TH2D*)this->hist->Clone();
      hist->GetXaxis()->SetRangeUser(xMin, foldMax);
      hist->Draw("COLZ");
      canvas->Print(fname);

      // Draw Efficiencies
      efficiency->Draw("AP");
      gPad->Update();  // Ensure painting
      efficiency->GetPaintedGraph()->GetYaxis()->SetRangeUser(EFF_YMIN, EFF_YMAX);
      efficiency->GetPaintedGraph()->GetXaxis()->SetLimits(xMin, foldMax);
      efficiency->GetPaintedGraph()->GetXaxis()->SetRangeUser(xMin, foldMax);
      efficiency->GetPaintedGraph()->GetXaxis()->SetNdivisions(510);
      gPad->Update();

      TLegend* efflegend = new TLegend(0.2, 0.7, 0.4, 0.9);
      efflegend->AddEntry(efficiency.get(), "Algorithm Efficiency");
      TLine *maxEffLine = new TLine(xMin, 0.99, foldMax, 0.99);
      maxEffLine->SetLineColor(kRed);
      maxEffLine->SetLineWidth(2);
      maxEffLine->SetLineStyle(4);
      efflegend->AddEntry(maxEffLine,"99% Efficiency");
      maxEffLine->Draw("SAME");
      efflegend->Draw("SAME");
      canvas->Print(fname);

      // Draw Resolution
      for (auto key: FITPARAM_VEC) {
	TH1D *hist = params->fromEnum(key);
	hist->Draw();
	hist->GetXaxis()->SetNdivisions(510);
	if (key == SIGMA) {
	  hist->GetYaxis()->SetRangeUser(RES_YMIN, RES_YMAX);
	  TLegend* resoLegend = new TLegend(0.2, 0.7, 0.4, 0.9);
	  TLine *minExpRes = new TLine(xMin, 15, foldMax, 15);
	  minExpRes->SetLineColor(kRed);
	  minExpRes->SetLineWidth(2);   
	  minExpRes->SetLineStyle(4);
	  resoLegend->AddEntry(minExpRes,"15 ps");
	  minExpRes->Draw("SAME");
	  resoLegend->Draw("SAME");
	} else if (key == BSIGMA) {
	  hist->GetYaxis()->SetRangeUser(RES_YMIN, 2*RES_YMAX);
	  TLegend* resoLegend = new TLegend(0.6, 0.8, 0.8, 0.9);
	  TLine *minExpRes = new TLine(xMin, 30, foldMax, 30);
	  minExpRes->SetLineColor(kRed);
	  minExpRes->SetLineWidth(2);   
	  minExpRes->SetLineStyle(4);
	  resoLegend->AddEntry(minExpRes,"30 ps");
	  minExpRes->Draw("SAME");
	  resoLegend->Draw("SAME");
	}
	canvas->Print(fname);
      }

      // draw slices o7
      TLatex latex;
      latex.SetTextSize(0.04);
      latex.SetTextAlign(13); 
      for (int iSlice=0; iSlice < slicesHists.size(); ++iSlice) {
	const auto& hSlice = slicesHists[iSlice];
	const auto& fit = slicesFits[iSlice];
	TLegend* thislegend = new TLegend(0.65, 0.75, 0.9, 0.9);
	TString restext;
	if (fit->GetNpar() == 5) { // Double Gaussian
	  restext = Form("#sigma_{1}^{dgaus}=%.2f(%.2f%%), #sigma_{2}^{dgaus}=%.2f(%.2f%%)",
			 fit->GetParameter("Sigma1"),100*fit->GetParError("Sigma1")/fit->GetParameter("Sigma1"),
			 fit->GetParameter("Sigma2"),100*fit->GetParError("Sigma2")/fit->GetParameter("Sigma2"));
	  thislegend->AddEntry(fit.get(),"Double Gaussian Fit");
	} else {
	  restext = Form("#sigma_{1}^{gaus}=%.2f(%.2f%%)",
			 fit->GetParameter(2),100*fit->GetParError(2)/fit->GetParameter(2));
	  thislegend->AddEntry(fit.get(),"Gaussian Fit");
	}
	
	thislegend->AddEntry(hSlice.get(),"Histogram");

	auto chi2text = Form("#chi^{2}=%.2f",fit->GetChisquare()/fit->GetNDF());
	  
	hSlice->GetYaxis()->SetTitleOffset(1.4);
	hSlice->Draw("HIST");
	fit->Draw("SAME");
	hSlice->GetXaxis()->SetRangeUser(-100, 100);
	latex.DrawLatexNDC(0.18, 0.90, restext);
	latex.DrawLatexNDC(0.18, 0.85, chi2text);
	thislegend->Draw("SAME");
	canvas->Print(fname); // slices
	
	// with log scale
	canvas->SetLogy(true);
	hSlice->Draw("HIST");
	fit->Draw("SAME");
	hSlice->GetXaxis()->SetRangeUser(-400, 400);
	latex.DrawLatexNDC(0.18, 0.90, restext);
	latex.DrawLatexNDC(0.18, 0.85, chi2text);
	thislegend->Draw("SAME");
	canvas->Print(fname); // slices
	canvas->SetLogy(false);
      }
      canvas->Print(Form("%s]",fname.Data()));
    }

    // -----------------------------------------------------------------------
    // printEfficiencyStats
    //   Prints a formatted per-bin efficiency table to stdout showing the
    //   scenario and score label, then one row per x-axis bin with the
    //   pass count, total count, fail count, and efficiency percentage.
    //   A summary footer shows the inclusive (integral) efficiency.
    // -----------------------------------------------------------------------
    inline void printEfficiencyStats() {
      const int W_BIN  = 22;
      const int W_COL  = 10;
      const std::string SEP(W_BIN + W_COL * 4 + 2, '-');

      // Header
      std::cout << '\n' << SEP << '\n';
      std::string headerLabel = std::string(this->times) + "  " + toStringShort(this->scoreToUse);
      // The data columns ("Pass", "Total", "Fail", "Eff (%)") are each w_col wide
      // and must end at the same position as the separator (w_bin + w_col*4).
      // Compute the width of the first column header field so it starts right
      // after the label and the whole row still aligns with the separator.
      int firstColW = (W_BIN + W_COL) - (int)headerLabel.size();
      if (firstColW < 1) firstColW = 1;
      std::cout << headerLabel
                << std::right
                << std::setw(firstColW) << "Pass"
                << std::setw(W_COL) << "Total"
                << std::setw(W_COL) << "Fail"
                << std::setw(W_COL) << "Eff (%)"
                << '\n';
      std::cout << std::left << std::setw(W_BIN) << (std::string("  dt vs. ") + this->xtitle)
                << '\n';
      std::cout << SEP << '\n';

      // Per-bin rows
      for (int i = 1; i <= effPass->GetNbinsX(); i++) {
        double numPass  = effPass->GetBinContent(i);
        double numTotal = effTotal->GetBinContent(i);
        double numFail  = numTotal - numPass;
        double eff      = (numTotal > 0) ? 100.0 * numPass / numTotal : 0.0;

        std::ostringstream binLabel;
        binLabel << std::fixed << std::setprecision(2)
                 << "[" << effPass->GetBinLowEdge(i)
                 << ", " << effPass->GetBinLowEdge(i + 1) << ")";

        std::cout << std::left  << std::setw(W_BIN) << binLabel.str()
                  << std::right << std::fixed << std::setprecision(1)
                  << std::setw(W_COL) << numPass
                  << std::setw(W_COL) << numTotal
                  << std::setw(W_COL) << numFail
                  << std::setw(W_COL) << eff
                  << '\n';
      }

      // Summary footer
      double totalPass = effPass->Integral();
      double totalAll  = effTotal->Integral();
      double overallEff = (totalAll > 0) ? 100.0 * totalPass / totalAll : 0.0;
      std::cout << SEP << '\n';
      std::cout << std::left  << std::setw(W_BIN) << "Total"
                << std::right << std::fixed << std::setprecision(1)
                << std::setw(W_COL) << totalPass
                << std::setw(W_COL) << totalAll
                << std::setw(W_COL) << (totalAll - totalPass)
                << std::setw(W_COL) << overallEff
                << '\n';
      std::cout << SEP << '\n';
    }
  };

  // ---------------------------------------------------------------------------
  // AnalysisObj
  //   Top-level container for one (timing scenario, Score) combination.
  //   Holds a map of PlotObjs keyed by variable name ("fjet", "ftrack",
  //   "pu_frac", "hs_track", "pu_track", "vtx_dz") plus an inclusive
  //   timing-residual histogram and an inclusive purity histogram.
  //
  //   The constructor books all PlotObjs with the histogram binning defined
  //   in clustering_constants.h.  File and histogram names are derived from
  //   the scenario label and Score enum so they are unique across scenarios.
  //
  //   Member methods:
  //     operator[]        — access a PlotObj by variable name (mutable / const)
  //     get               — const access with bounds checking
  //     postProcessing    — calls plotPostProcessing() on all PlotObjs
  //     fullPlotting      — calls plotLogic() on all PlotObjs
  //     printEfficiencyStats — delegates to the named PlotObj's method
  // ---------------------------------------------------------------------------
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
      TString timeIDer(timetypeIDer);
      timeIDer.ReplaceAll(" ", "");
      timeIDer.ReplaceAll(".", "");
      timeIDer.ReplaceAll("+", "");
      timeIDer.ToLower();

      TString scoreIDer(toString(score));
      scoreIDer.ReplaceAll("_", "");
      scoreIDer.ReplaceAll(" ", "_");
      scoreIDer.ReplaceAll(".", "");
      scoreIDer.ReplaceAll("{", "");
      scoreIDer.ReplaceAll("}", "");
      scoreIDer.ReplaceAll("(", "");
      scoreIDer.ReplaceAll(")", "");
      scoreIDer.ReplaceAll("#", "");
      scoreIDer.ReplaceAll("-", "");
      scoreIDer.ReplaceAll("|", "");
      scoreIDer.ToLower();

      TString filenameIDer = timeIDer + "_" + scoreIDer;
  
      dataObjects["fjet"] = std::make_unique<PlotObj>(
        "n Forward Jets", timetypeIDer,
	Form("../figs/fullplots/%s_nfjet.pdf",filenameIDer.Data()),
	score,
	FJET_MIN  , FJET_MAX  , FJET_WIDTH  ,
	DIFF_MIN  , DIFF_MAX  , DIFF_WIDTH  ,
	PURITY_MIN, PURITY_MAX, PURITY_WIDTH,
	FOLD_FJET, FOLD_FJET+FJET_WIDTH     );

      dataObjects["vtx_dz"] = std::make_unique<PlotObj>(
        "|Reco HS z - Truth HS z| (mm)", timetypeIDer,
	Form("../figs/fullplots/%s_vtx_dz.pdf",filenameIDer.Data()),
	score,
	VTX_DZ_MIN  , VTX_DZ_MAX  , VTX_DZ_WIDTH,
	DIFF_MIN    , DIFF_MAX    , DIFF_WIDTH  ,
	PURITY_MIN  , PURITY_MAX  , PURITY_WIDTH,
	FOLD_VTX_DZ , FOLD_VTX_DZ+VTX_DZ_WIDTH  );
  
      dataObjects["ftrack"] = std::make_unique<PlotObj>(
        "n Forward Tracks", timetypeIDer,
	Form("../figs/fullplots/%s_ntrack.pdf",filenameIDer.Data()),
	score,
	TRACK_MIN , TRACK_MAX , TRACK_WIDTH ,
	DIFF_MIN  , DIFF_MAX  , DIFF_WIDTH  ,
	PURITY_MIN, PURITY_MAX, PURITY_WIDTH,
	FOLD_TRACK, FOLD_TRACK+TRACK_WIDTH  );
  
      dataObjects["pu_frac"] = std::make_unique<PlotObj>(
        "Pile Up Fraction", timetypeIDer, 
        Form("../figs/fullplots/%s_pufrac.pdf",filenameIDer.Data()),
	score,
	PU_FRAC_MIN , PU_FRAC_MAX, PU_FRAC_WIDTH,
	DIFF_MIN    , DIFF_MAX   , DIFF_WIDTH   ,
	PURITY_MIN  , PURITY_MAX , PURITY_WIDTH ,
	FOLD_PU_FRAC, PU_FRAC_MAX              );
  
      dataObjects["hs_track"] = std::make_unique<PlotObj>(
        "n Forward HS Tracks", timetypeIDer, 
	Form("../figs/fullplots/%s_nhstrack.pdf",filenameIDer.Data()),
	score,
	HS_TRACK_MIN, HS_TRACK_MAX, HS_TRACK_WIDTH ,
	DIFF_MIN    , DIFF_MAX    , DIFF_WIDTH     ,
	PURITY_MIN  , PURITY_MAX  , PURITY_WIDTH   ,
	FOLD_HS_TRACK, FOLD_HS_TRACK+HS_TRACK_WIDTH);
  
      dataObjects["pu_track"] = std::make_unique<PlotObj>(
        "n Forward PU Tracks", timetypeIDer, 
	Form("../figs/fullplots/%s_nputrack.pdf",filenameIDer.Data()),
	score,
	PU_TRACK_MIN, PU_TRACK_MAX, PU_TRACK_WIDTH ,
	DIFF_MIN    , DIFF_MAX    , DIFF_WIDTH     ,
	PURITY_MIN  , PURITY_MAX  , PURITY_WIDTH   ,
	FOLD_PU_TRACK, FOLD_PU_TRACK+PU_TRACK_WIDTH);

      inclusiveReso = std::make_unique<TH1D>(
        Form("reso_%s", filenameIDer.Data()),
	Form("Inclusive %s t_{0} - TruthVtx t_{0} (%s);#Delta t[ps];Entries",
	     toString(score), timetypeIDer),
	(int)((DIFF_MAX-DIFF_MIN)/DIFF_WIDTH), DIFF_MIN, DIFF_MAX);

      inclusivePurity = std::make_unique<TH1D>(
	Form("purity_%s", filenameIDer.Data()),
	Form("%s Purity (%s);Purity;Entries", toString(score), timetypeIDer),
	(int)((PURITY_MAX-PURITY_MIN)/PURITY_WIDTH), PURITY_MIN, PURITY_MAX);
    }

    // -----------------------------------------------------------------------
    // postProcessing
    //   Calls plotPostProcessing() on every PlotObj in dataObjects.
    //   Must be called after the event loop and before fullPlotting.
    // -----------------------------------------------------------------------
    inline void postProcessing() {
      for (auto &[str, plt]: this->dataObjects)
	plt->plotPostProcessing();
    }

    // -----------------------------------------------------------------------
    // fullPlotting
    //   Calls plotLogic() on every PlotObj in dataObjects, writing one
    //   multi-page PDF per (variable, Score) combination to figs/.
    // -----------------------------------------------------------------------
    inline void fullPlotting(TCanvas *canvas) {
      for (const auto& [str, plt]: this->dataObjects)
	plt->plotLogic(canvas);
    }

    // -----------------------------------------------------------------------
    // printEfficiencyStats
    //   Delegates to PlotObj::printEfficiencyStats() for the PlotObj
    //   identified by key (e.g. "pu_frac").
    // -----------------------------------------------------------------------
    inline void printEfficiencyStats(std::string key) {
      dataObjects[key]->printEfficiencyStats();
    }
  };

  // ---------------------------------------------------------------------------
  // moneyPlot
  //   Generates a three-page comparison PDF for a given kinematic variable
  //   (key) across the set of AnalysisObj* in plts.  Pages:
  //     1. Core timing resolution σ vs x (with 15 ps reference line)
  //     2. Background timing resolution σ_bkg vs x (with 15 ps reference line)
  //     3. Efficiency vs x (with 99% reference line)
  //   Each page overlays all provided AnalysisObj results with distinct colours
  //   and adds a normalised event-count histogram (from plts[0]) in the
  //   background to show where statistics are concentrated.
  // ---------------------------------------------------------------------------
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
	double yMin, double yMax, double refVal = -1,
	const char* refLabel = nullptr
    ) {
      TLegend* legend = new TLegend(0.55, 0.65, 0.9, 0.9);
      bool first = true;
      int counter = 0;
      bool isEfficiency = false;
      for (const auto& ana : plts) {
	const auto& plt = ana->get(key);
	auto obj = getter(plt);
	obj->SetLineColor(COLORS[counter++ % COLORS.size()]);
	legend->AddEntry(obj, Form("%s: %s", toString(plt->scoreToUse), plt->times));

	if (first) {
	  obj->SetTitle(Form("%s vs %s", title, plt->xtitle));
	  obj->Draw("E1");
	  if constexpr (std::is_same_v<decltype(obj), TEfficiency*> || std::is_same_v<decltype(obj), std::shared_ptr<TEfficiency>>) {
	    gPad->Update();
	    obj->GetPaintedGraph()->GetXaxis()->SetRangeUser(xMin, xMax);
	    obj->GetPaintedGraph()->GetYaxis()->SetRangeUser(yMin, yMax);
	  } else {
	    obj->GetXaxis()->SetRangeUser(xMin, xMax);
	    obj->GetYaxis()->SetRangeUser(yMin, yMax);
	  }
	  first = false;
	} else
	  obj->Draw("E1 SAME");
      }

      auto* total = (TH1D*)plts[0]->get(key)->effTotal->Clone();
      total->SetLineColorAlpha(COLORS[8],0.5);
      total->SetFillColorAlpha(COLORS[8],0.5);
      total->SetLineWidth(2);
      total->Scale(1.0 / total->Integral());
      total->Scale((0.2*yMax - yMin) / total->GetMaximum());
      total->Draw("HIST SAME");

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

    // Plot Background Resolution
    plotWithLegend(
        [](auto& plt) { return plt->params->bkgSigmaDist; },
	"Background #sigma", RES_YMIN, RES_YMAX,
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


  // ---------------------------------------------------------------------------
  // inclusivePlot
  //   Plots the inclusive timing residual histogram for each AnalysisObj in
  //   plts, overlaying a double-Gaussian fit (createDblFit) and annotating
  //   with σ_1, σ_2 from the fit.  Supports linear and logarithmic y-axis
  //   (controlled by logScale) and custom x-axis range (xMin, xMax).  Each
  //   AnalysisObj gets its own page in the output PDF.
  // ---------------------------------------------------------------------------
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
      TF1* fit1 = createDblFit(hist,true);
      TLegend* inclusiveLegend = new TLegend(0.65, 0.75, 0.9, 0.9);

      inclusiveLegend->AddEntry(hist,"Histogram");
      inclusiveLegend->AddEntry(fit1,"Double Gaussian Fit");
    
      hist->GetXaxis()->SetRangeUser(xMin, xMax);
      hist->Draw("HIST");
      fit1->Draw("SAME");
      inclusiveLegend->Draw("SAME");
    
      double dgSigma1 = fit1->GetParameter("Sigma1");
      double dgSigma2 = fit1->GetParameter("Sigma2");
      latex.DrawLatexNDC(0.18, 0.90,Form("#sigma_{1}^{dgaus}=%.2f",dgSigma1));
      latex.DrawLatexNDC(0.18, 0.85,Form("#sigma_{2}^{dgaus}=%.2f",dgSigma2));
      canvas->Print(fname);
    }
    canvas->Print(Form("%s]",fname));
  }

  // ---------------------------------------------------------------------------
  // purityPlot
  //   Draws normalised cluster-purity distributions for each Score in
  //   purityMap on a single canvas page, using the COLORS palette.  Overflow
  //   is folded into the last bin before normalisation.  logscale toggles
  //   the y-axis.  The provided legend is drawn on the same pad.
  // ---------------------------------------------------------------------------
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
