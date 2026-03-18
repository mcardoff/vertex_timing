#pragma once
//
// ATLAS Style, based on a style file from BaBar
// M.Sutton, ATLAS Collaboration 2010
//
// Converted to a single header-only file (inline functions) for use in
// compiled C++ projects without CINT / gROOT::LoadMacro.
//

#include "TROOT.h"
#include "TStyle.h"

/// Constructs and returns a new TStyle configured with ATLAS plot conventions.
inline TStyle* AtlasStyle() {
  TStyle* atlasStyle = new TStyle("ATLAS", "Atlas style");

  // — Background colours (all white) ——————————————————————————————————————
  const Int_t icol = 0;
  atlasStyle->SetFrameBorderMode(icol);
  atlasStyle->SetFrameFillColor(icol);
  atlasStyle->SetCanvasBorderMode(icol);
  atlasStyle->SetCanvasColor(icol);
  atlasStyle->SetPadBorderMode(icol);
  atlasStyle->SetPadColor(icol);
  atlasStyle->SetStatColor(icol);

  // — Paper size (A4) ———————————————————————————————————————————————————————
  atlasStyle->SetPaperSize(20, 26);

  // — Pad margins ——————————————————————————————————————————————————————————
  atlasStyle->SetPadTopMargin(0.05);
  atlasStyle->SetPadRightMargin(0.05);
  atlasStyle->SetPadBottomMargin(0.16);
  atlasStyle->SetPadLeftMargin(0.16);

  // — Axis title offsets ————————————————————————————————————————————————————
  atlasStyle->SetTitleXOffset(1.4);
  atlasStyle->SetTitleYOffset(1.4);

  // — Fonts (42 = Helvetica, precision 2) ——————————————————————————————————
  const Int_t    font  = 42;
  const Double_t tsize = 0.05;
  atlasStyle->SetTextFont(font);
  atlasStyle->SetTextSize(tsize);
  for (const char* ax : {"x", "y", "z"}) {
    atlasStyle->SetLabelFont(font, ax);
    atlasStyle->SetTitleFont(font, ax);
    atlasStyle->SetLabelSize(tsize, ax);
    atlasStyle->SetTitleSize(tsize, ax);
  }

  // — Markers & lines ———————————————————————————————————————————————————————
  atlasStyle->SetMarkerSize(0.8);
  atlasStyle->SetHistLineWidth(2.);
  atlasStyle->SetLineStyleString(2, "[12 12]");  // dashed line style 2

  // — Error bars ————————————————————————————————————————————————————————————
  atlasStyle->SetErrorX(0.5);         // half-bin-width horizontal error bars
  atlasStyle->SetEndErrorSize(0.);

  // — Suppress statistics/title/fit boxes ———————————————————————————————————
  atlasStyle->SetOptTitle(0);
  atlasStyle->SetOptStat(0);
  atlasStyle->SetOptFit(0);

  // — Tick marks on all four sides of each pad ——————————————————————————————
  atlasStyle->SetPadTickX(1);
  atlasStyle->SetPadTickY(1);

  return atlasStyle;
}

/// Apply the ATLAS style globally.  Safe to call multiple times.
inline void SetAtlasStyle() {
  static TStyle* atlasStyle = nullptr;
  if (!atlasStyle) atlasStyle = AtlasStyle();
  gROOT->SetStyle("ATLAS");
  gROOT->ForceStyle();
}
