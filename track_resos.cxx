#include <Rtypes.h>
#include <TBranch.h>
#include <TChain.h>
#include <TFile.h>
#include <TDirectory.h>
#include <TCanvas.h>
#include <TString.h>
#include <TStyle.h>
#include <TColor.h>
#include <TMath.h>
#include <TTreeReader.h>
#include <TTreeReaderArray.h>
#include <TLatex.h>
#include <TF1.h>
#include <TH1.h>
#include <TH2.h>
#include <TLegend.h>
#include <TProfile.h>
#include <TGraphErrors.h>
#include <TGraphAsymmErrors.h>
#include <TEfficiency.h>
#include <TGraph.h>
#include <TGaxis.h>

#include <bits/types/struct_sched_param.h>
#include <boost/filesystem.hpp>
#include <cassert>
#include <cmath>
#include <iostream>
#include <ostream>
#include <utility>
#include <vector>

#include "./clustering_utilities.h"

#define debug false

using namespace myutl;

struct TrackFitParams {
  TH1D* mean_dist;
  TH1D* sigma_dist;
  TH1D* core_amp_dist;
  TH1D* back_amp_dist;
  TH1D* amp_ratio_dist;

  TrackFitParams(const char* title, const char* name,
	    double x_min, double fold_val, double fold_max,
	    double x_width, int color) {
    auto full_title = Form("%%s vs %s;%s;%%s", title, title);
    std::cout << "FOLD VAL: " << fold_val << "\n";
    std::cout << "X MIN: " << x_min << "\n";
    std::cout << "X WID: " << x_width << "\n";
    mean_dist = new TH1D(Form("mean_dist_%s", name),
			 Form(full_title, "Common #mu", "#mu"),
			 (int)((fold_val-x_min)/x_width), x_min, fold_max);
    
    sigma_dist = new TH1D(Form("sigma_dist_%s", name),
			  Form(full_title, "Core #sigma", "Core #sigma"),
			  (int)((fold_val-x_min)/x_width), x_min, fold_max);
    
    core_amp_dist = new TH1D(Form("amp1_dist_%s", name),
			     Form(full_title, "Core Amplitude", "Amplitude"),
			     (int)((fold_val-x_min)/x_width), x_min, fold_max);
    
    back_amp_dist = new TH1D(Form("amp2_dist_%s", name),
			     Form(full_title, "Bkg Amplitude", "Amplitude"),
			     (int)((fold_val-x_min)/x_width), x_min, fold_max);
    
    amp_ratio_dist = new TH1D(Form("ampratio_dist_%s", name),
			      Form(full_title, "Core Amp/Bkg Amp", "A.U."),
			      (int)((fold_val-x_min)/x_width), x_min, fold_max);

    std::cout << "N BINS SIGMA: " << sigma_dist->GetNbinsX() << "\n";
    // mean_dist->SetMaximum(60.0);
    // mean_dist->SetMinimum(00.0);
    
    mean_dist->SetMaximum( 1.0);
    mean_dist->SetMinimum(-1.0);
    
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

  void set_param_maxes() {
    for (auto hist: {sigma_dist, mean_dist, core_amp_dist, back_amp_dist, amp_ratio_dist}) {
      std::cout << hist->GetMaximum() << "\n";
      hist->GetYaxis()->SetRangeUser(0.5*hist->GetMinimum() ,1.5*hist->GetMaximum());
    }
  }

  TH1D* from_enum(FitParamFields fit) {
    switch (fit) {
    case MEAN:  return mean_dist;
    case SIGMA: return sigma_dist;
    case CORE:  return core_amp_dist;
    case BKG:   return back_amp_dist;
    case RATIO: return amp_ratio_dist;
    default:    return nullptr;
    }
  }
};

struct TrackPlotObj {
  double bin_width;
  double fold_value, fold_max;
  const char* fname;
  double x_min, x_max, x_wid;
  double y_min, y_max, y_wid;
  TH1D* eff_total; // doesnt need to be a map, same for all
  TH1D* eff_pass;
  TH2D* hist;
  TH2D* purity;
  TrackFitParams params;
  std::vector<std::pair<TF1*,TH1D*>> slices;
  TEfficiency* efficiency;

  TrackPlotObj(const char* title, const char* name, const char* fname,
	  double x_min, double x_max, double x_wid,
	  double y_min, double y_max, double y_wid,
	  double fold_val, double fold_max, int color)
    : bin_width(x_wid), fname(fname),
      fold_value(fold_val), fold_max(fold_max),
      x_min(x_min), x_max(x_max), x_wid(x_wid),
      y_min(y_min), y_max(y_max), y_wid(y_wid),
      params(title, name, x_min, fold_val, fold_max, x_wid, color)
  {
    hist = new TH2D(name,
		    Form("Track t - Truth Particle t vs %s;%s;#Delta t[ps]", title,title),
		    (int)((x_max-x_min)/x_wid), x_min, x_max,(int)((y_max-y_min)/y_wid), y_min, y_max);
    
    eff_pass = new TH1D(Form("good_%s",name),
			Form("%s Eff Events;Entries;%s",title,title),
			(int)((fold_val-x_min)/x_wid), x_min, fold_max);

    std::cout << "N BINS SIGMA: " << eff_pass->GetNbinsX() << "\n";

    eff_total = new TH1D(Form("flat_%s",name),
			 Form("%s All Events;Entries;%s",title,title),
			 (int)((fold_val-x_min)/x_wid), x_min, fold_max);
    
    slices = std::vector<std::pair<TF1*,TH1D*>>{};
  }
};

void track_resos() {
  gStyle->SetOptStat(0);
  // gPad->SetLeftMargin(0.1);

  TChain chain ("ntuple");
  setup_chain(chain, "../ntuple-hgtd");
  TTreeReader reader(&chain);
  BranchPointerWrapper branch(reader);

  TCanvas *canvas = new TCanvas("canvas", "Histograms", 800, 600);
  canvas->SetLeftMargin(0.15);
  // gPad->SetLeftMargin(0.1);

  double diff_min = -7000.0, diff_max = 7000.0;
  double diff_width = 2.0;

  double trk_eta_min = 2.4, trk_eta_max = 4.0 , fold_eta = trk_eta_max;
  double trk_eta_width = 0.1;

  double trk_pt_min = 0.0, trk_pt_max = 500.0 , fold_pt = 10.0;
  double trk_pt_width = 1.0;

  double trk_hits_min = 1, trk_hits_max = 4 , fold_hits = 4;
  double trk_hits_width = 1;

  double z_min = -200, z_max = 200, fold_z = 100;
  double z_width = 10.0;

  // All Needed Histogramming information

  TH1D* inclusive_resos = new TH1D("inclusivereso","Inclusive Track t - Truth Particle t;#Delta t[ps];Entries",
				   (int)((diff_max-diff_min)/diff_width), diff_min, diff_max);
  TH1D* lasttrackinglayer_resos = new TH1D("lastlayerhit","Tracks with hit in last layer Track t - Truth Particle t;#Delta t[ps];Entries",
					   (int)((diff_max-diff_min)/diff_width), diff_min, diff_max);
    
  TrackPlotObj track_eta("Track #eta" , "trk_eta"   , "figs/track_eta.pdf",
			 trk_eta_min  , trk_eta_max , trk_eta_width       ,
			 diff_min     , diff_max    , diff_width          ,
			 fold_eta     , trk_eta_max , colors[0]           );
  
  TrackPlotObj track_pt ("Track p_{T}", "trk_pt"    , "figs/track_pt.pdf",
			 trk_pt_min   , trk_pt_max  , trk_pt_width       ,
			 diff_min     , diff_max    , diff_width         ,
			 fold_pt      , fold_pt+trk_pt_width  , colors[1]);
  
  TrackPlotObj nhits    ("n HGTD Hits", "hits_hgtd" , "figs/track_hgtdhits.pdf",
			 trk_hits_min , trk_hits_max, trk_hits_width           ,
			 diff_min     , diff_max    , diff_width               ,
			 fold_hits    , fold_hits   , colors[2]    );
  
  TrackPlotObj recovtx_z("Reco Vtx z (mm)", "z", "figs/track_vtx_z.pdf",
			 z_min        , z_max       , z_width          ,
			 diff_min     , diff_max    , diff_width       ,
			 z_max        , z_max       , colors[3]        );
  
  while (reader.Next()) {
    std::string filename = chain.GetFile()->GetName(); // file we're in
    Long64_t this_evnt = chain.GetReadEntry() - chain.GetChainOffset();
    
    if (branch.topojet_pt.GetSize() < 1) {
      if(debug)
	std::cout << "Skipping low jet event" << std::endl;
      continue;  // Skip if no jets are present
    }
    
    // check reco HS vertex is with 2mm of truth HS vertex
    if(std::abs(branch.truth_vtx_z[0] - branch.reco_vtx_z[0]) > max_vtx_dz) {
      if(debug)
	std::cout << "Skipping event due to incorrect HS vertex" << std::endl;
      continue;
    }

    // check if there is one forward jet with pt > 30 GeV
    bool pass_jet_pt_cut = false;
    int nForwardJet = 0;
    for(int idx = 0; idx < branch.topojet_eta.GetSize(); ++idx) {
      float eta = branch.topojet_eta[idx], pt = branch.topojet_pt[idx];
      if (std::abs(eta) > min_abs_eta_jet) {
	nForwardJet++;
	if (!pass_jet_pt_cut && pt > min_jetpt) {
	  // check if the pt > 30
	  pass_jet_pt_cut = true;
	}
      }
    }

    for (int i=0; i < branch.track_eta.GetSize(); i++) {
      
      if (branch.track_pt[i] < 1.0)
	continue;

      double
	this_track_pt  = branch.track_pt[i],
	this_track_eta = std::abs(branch.track_eta[i]),
	this_track_hgtd_hits = branch.track_hgtd_hits[i],
	this_track_prim_hits = branch.track_prim_hits[i];

      int
	part = branch.track_to_particle[i],
	t_vtx = branch.track_to_truthvtx[i];

      if (part == -1)
	continue;

      // NOW We go no further
      if (!pass_jet_pt_cut)
	continue;

      double reco_z = t_vtx != -1 ? branch.truth_vtx_z[t_vtx] : -1e100;
      
      auto eff_fill_val_hgtd_hits  = this_track_hgtd_hits;
      auto eff_fill_val_track_pt   = (this_track_pt  >= fold_pt ) ? fold_pt : this_track_pt ;
      auto eff_fill_val_track_eta  = this_track_eta;
      auto eff_fill_val_z          = reco_z;
    
      nhits.eff_total->    Fill(eff_fill_val_hgtd_hits);
      track_pt.eff_total-> Fill(eff_fill_val_track_pt );
      track_eta.eff_total->Fill(eff_fill_val_track_eta);
      recovtx_z.eff_total->Fill(eff_fill_val_z        );

      if (branch.track_time_valid[i] == 0)
	continue;

      float diff = branch.track_time[i] - branch.particle_t[part];
      // if (diff_min > diff || diff > diff_max)
	// std::cout << "DIFF IS MISSED :C " << diff << "\n";
      inclusive_resos->Fill(diff);
      
      if (std::abs(diff) < 90) {
	nhits.eff_pass->    Fill(eff_fill_val_hgtd_hits);
	track_pt.eff_pass-> Fill(eff_fill_val_track_pt );
	track_eta.eff_pass->Fill(eff_fill_val_track_eta);
	recovtx_z.eff_pass->Fill(eff_fill_val_z        );
      }

      // fill diff hists
      nhits.hist->    Fill(this_track_hgtd_hits, diff);
      track_eta.hist->Fill(this_track_eta      , diff);
      track_pt.hist-> Fill(this_track_pt       , diff);
      recovtx_z.hist->Fill(reco_z              , diff);
    }    
  }

  // TF1 *dgaus_fitfunc = new TF1("dgaus", "[1] * TMath::Exp(-0.5 * ((x - [0]) / [3])^2) / ([3] * TMath::Sqrt(2*TMath::Pi())) + [2] * TMath::Exp(-0.5 * ((x - [0]) / [4])^2) / ([4] * TMath::Sqrt(2*TMath::Pi()))", diff_min, diff_max);
  TF1 *dgaus_fitfunc = new TF1("dgaus", "[1]*TMath::Exp(-0.5 * ((x - [0]) / [3])^2) + [2]*TMath::Exp(-0.5 * ((x - [0]) / [4])^2)", diff_min, diff_max);
  dgaus_fitfunc->SetParNames("Mean", "Norm1", "Norm2", "Sigma1", "Sigma2");

  
  TLegend *algolegend = new TLegend(0.2, 0.7, 0.4, 0.9);

  // add all objects to plotobjs before drawing everything
  std::vector<TrackPlotObj> plots = {track_eta, track_pt, recovtx_z, nhits};
  for (int i = 0; i < plots.size(); i++) {
    auto fname = plots[i].fname;
    double fold_val = plots[i].fold_value;
    // Draw Comparisons of all the resos
    TH2D* hist = plots[i].hist;
    TH1D *pass = plots[i].eff_pass;
    TH1D *total = plots[i].eff_total;

    // these are all same through diff scores
    auto xtitle = hist->GetXaxis()->GetTitle(), ytitle = hist->GetYaxis()->GetTitle();
    int xmin = total->GetXaxis()->GetXmin(), xmax = total->GetXaxis()->GetXmax();
    auto xbin = total->GetNbinsX();

    for(int j = 0; j < hist->GetNbinsX(); ++j) {
      auto color = colors[j % colors.size()];
      auto name = Form("%s_slice%d",hist->GetName(),j);
      
      int right = hist->GetXaxis()->GetBinLowEdge(j+1) >= fold_val
	? hist->GetNbinsX()+1
	: j+1;
      
      TH1D* hSlice = hist->ProjectionY(name, j+1, right);

      if (hSlice->Integral() == 0 || hSlice->GetEntries() < 6)
	continue;

      hSlice->SetLineColor(color);
      hSlice->SetLineWidth(2);

      double leftEdge  = hist->GetXaxis()->GetBinLowEdge(j+1);
      double rightEdge = hist->GetXaxis()->GetBinLowEdge(right+1);
      double bin_center = hist->GetXaxis()->GetBinCenter(j+1);

      if (leftEdge < 1.0 && leftEdge > 0.0)
	leftEdge *= 100;
      if (rightEdge < 1.0 && rightEdge > 0.0)
	rightEdge *= 100;

      // set correct title
      hSlice->SetTitle(Form("%s #in [%d,%d);%s;A.U.",
			    hist->GetTitle(),
			    (int)leftEdge,(int)rightEdge,
			    ytitle));
      
      // Perform Double Gaussian Fit
      TF1* dgaus_fit = create_dgaus_fit(hSlice,false);
      
      plots[i].params.fill_each(j, dgaus_fit);
      plots[i].slices.push_back(std::make_pair(dgaus_fit,hSlice));
	
      if (right == hist->GetNbinsX()+1)
	break;
    }

    TEfficiency *efficiency_plot = new TEfficiency(*pass, *total);
    efficiency_plot->SetStatisticOption(TEfficiency::kFNormal);  // use normal approximation
    efficiency_plot->SetTitle(Form("Efficiency vs %s;%s;Efficiency", xtitle, xtitle));
    efficiency_plot->SetLineColor(colors[i]);
    efficiency_plot->SetLineWidth(2);
    plots[i].efficiency = efficiency_plot;
    plots[i].params.set_param_maxes();
  }
  
  std::cout << "FINISHED CREATING " << std::endl;
  
  for (int i = 0; i < plots.size(); i++) {
    // Draw 2Ds
    TrackPlotObj plot = plots[i];
    TH2D *hist = plot.hist;
    canvas->Print(Form("%s[",plot.fname));
    hist->GetXaxis()->SetRangeUser(hist->GetXaxis()->GetXmin(), plot.fold_max);
    hist->Draw("COLZ");
    canvas->Print(Form("%s",plot.fname));

    // Draw Efficiencies
    auto  eff = plot.efficiency;
    eff->SetMarkerColor(colors[i]);
    eff->SetLineColor(colors[i]);
    eff->SetLineWidth(2);
    eff->Draw("AP");
    gPad->Update();  // Ensure painting
    eff->GetPaintedGraph()->GetYaxis()->SetRangeUser(0.0, 1.1);
    // std::cout << "FOLD MAX" << plot.fold_max << std::endl;
    eff->GetPaintedGraph()->GetXaxis()->SetLimits(plot.x_min, plot.fold_max);
    eff->GetPaintedGraph()->GetXaxis()->SetRangeUser(plot.x_min, plot.fold_max);
    eff->GetPaintedGraph()->GetXaxis()->SetNdivisions(510);
    eff->Print("ALL");
    gPad->Update();
    
    TLegend* efflegend = (TLegend*)algolegend->Clone();
    efflegend->AddEntry(eff,"Algorithm Efficiency");
    TLine *max_eff_line = new TLine(plot.x_min, 0.99, plot.fold_max, 0.99);
    max_eff_line->SetLineColor(kRed);
    max_eff_line->SetLineWidth(2);
    max_eff_line->SetLineStyle(4);
    efflegend->AddEntry(max_eff_line,"99% Efficiency");
    max_eff_line->Draw("SAME");
    efflegend->Draw("SAME");
    canvas->Print(Form("%s",plot.fname));

    // Draw Resolution
    for (auto key: fitparam_vec) {
      TH1D *hist = plot.params.from_enum(key);
      hist->Draw();
      hist->GetXaxis()->SetNdivisions(510);
      
      
      if (key == SIGMA) {
	TLegend* resolegend = (TLegend*)algolegend->Clone();
	TLine *min_exp_res = new TLine(plot.x_min, 30, plot.fold_max, 30);
	min_exp_res->SetLineColor(kRed);
	min_exp_res->SetLineWidth(2);   
	min_exp_res->SetLineStyle(4);

	resolegend->AddEntry(min_exp_res,"30ps");
	min_exp_res->Draw("SAME");
	resolegend->Draw("SAME");
      }
      canvas->Print(Form("%s",plot.fname));
    }

    // draw slices o7
    TLatex latex;
    latex.SetTextSize(0.04);
    latex.SetTextAlign(13);  // Align left-top
    for (const auto& [fit, hSlice]: plot.slices) {
      TLegend* thislegend = new TLegend(0.65, 0.75, 0.9, 0.9);
      auto res1text = Form("#sigma_{core}^{dgaus}=%.2f",fit->GetParameter(3));
      auto res2text = Form("#sigma_{bkg}^{dgaus}=%.2f ",fit->GetParameter(4));
      thislegend->AddEntry(hSlice,"Histogram");
      thislegend->AddEntry(fit,"Double Gaussian Fit");

      hSlice->GetYaxis()->SetTitleOffset(1.4);
      hSlice->Draw("HIST");
      fit->Draw("SAME");
      hSlice->GetXaxis()->SetRangeUser(-100, 100);
      latex.DrawLatexNDC(0.18, 0.90, res1text);
      latex.DrawLatexNDC(0.18, 0.85, res2text);
      thislegend->Draw("SAME");
      canvas->Print(Form("%s",plot.fname)); // slices
	
      // with log scale
      canvas->SetLogy(true);
      hSlice->Draw("HIST");
      fit->Draw("SAME");
      hSlice->GetXaxis()->SetRangeUser(-1000, 1000);
      latex.DrawLatexNDC(0.18, 0.90, res1text);
      latex.DrawLatexNDC(0.18, 0.85, res2text);
      thislegend->Draw("SAME");
      canvas->Print(Form("%s",plot.fname)); // slices
      canvas->SetLogy(false);
    }
    canvas->Print(Form("%s]",plot.fname));
  }  

  auto inclusive_fname = "figs/inclusive_track_resos_logscale.pdf";
  TLatex latex;
  latex.SetTextSize(0.04);
  latex.SetTextAlign(13);  // Align left-top
  // plot inclusive fits
  canvas->Print(Form("%s[",inclusive_fname));
  canvas->SetLogy(true);
  std::vector<std::tuple<TH1D*,TF1*,TLegend*>> hist_fit_vec;
  TH1D *hist = inclusive_resos;
  TF1* fit = create_dgaus_fit(hist,false);
  TLegend* inclusive_legend = new TLegend(0.65, 0.75, 0.9, 0.9);

  inclusive_legend->AddEntry(hist,"Histogram");
  inclusive_legend->AddEntry(fit,"Double Gaussian Fit");
    
  hist->GetXaxis()->SetRangeUser(-400, 400);
  hist->Draw("HIST");
  fit->Draw("SAME");
  inclusive_legend->Draw("SAME");
    
  double dg_sigma = fit->GetParameter(3);
  latex.DrawLatexNDC(0.18, 0.90,Form("#sigma_{1}^{dgaus}=%.2f",dg_sigma));
  canvas->Print(Form("%s",inclusive_fname));
    
  canvas->Print(Form("%s]",inclusive_fname));

  // linear scale
  inclusive_fname = "figs/inclusive_track_resos_linscale.pdf";
  canvas->SetLogy(false);
  canvas->Print(Form("%s[",inclusive_fname));
  hist = inclusive_resos;
  hist->GetXaxis()->SetRangeUser(-200, 200);
  hist->Draw("HIST");
  fit->Draw("SAME");
  inclusive_legend->Draw("SAME");
  
  latex.DrawLatexNDC(0.18, 0.90, Form("#sigma_{1}^{dgaus}=%.2f", dg_sigma));
  canvas->Print(Form("%s",inclusive_fname));
  canvas->Print(Form("%s]",inclusive_fname));
}
