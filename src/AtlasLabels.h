#pragma once
//
// ATLAS plot labels
// Based on AtlasLabels.C from the ATLAS Collaboration
//
// Converted to a single header-only file (inline functions) for use in
// compiled C++ projects without CINT / gROOT::LoadMacro.
//

#include "TLatex.h"
#include "TLegend.h"

/// Apply the ATLAS-style legend appearance (border, white fill, Helvetica font).
/// Call immediately after `new TLegend(...)` before adding entries.
/// @param textSize  Label text size (default 0.04 — slightly smaller than main text)
inline void StyleLegend(TLegend* leg, float textSize = 0.03) {
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetFillStyle(0);
  leg->SetTextFont(42);
  leg->SetTextSize(textSize);
}

/// Draw "ATLAS <text>" on the current pad using NDC coordinates.
///
/// @param x      NDC x-position of the left edge of "ATLAS"
/// @param y      NDC y-position (baseline)
/// @param text   Text to draw after "ATLAS" (e.g. "Simulation Internal").
///               Pass nullptr to draw only "ATLAS".
/// @param color  Text colour (default: black = 1)
///
/// Typical usage (after drawing histogram/graph, before canvas->Print):
///   ATLASLabel(0.20, 0.88, "Simulation Internal");
inline void ATLASLabel(Double_t x, Double_t y,
                       const char* text = nullptr, Color_t color = 1) {
  TLatex l;
  l.SetNDC();
  l.SetTextFont(72);   // bold italic — the standard ATLAS font
  l.SetTextColor(color);

  // Horizontal gap between "ATLAS" and the trailing text, scaled to the
  // current pad's aspect ratio so it looks right on any canvas size.
  const double delx = 0.115 * 696 * gPad->GetWh() / (472 * gPad->GetWw());
  l.DrawLatex(x, y, "ATLAS");

  if (text) {
    TLatex p;
    p.SetNDC();
    p.SetTextFont(42);  // regular Helvetica for the trailing label
    p.SetTextColor(color);
    p.DrawLatex(x + delx, y, text);
  }
}

/// Draw a secondary energy / luminosity label (e.g. "#sqrt{s} = 14 TeV")
/// below the main ATLASLabel.  Typically called with y = ATLASLabel_y - 0.06.
///
/// @param x     NDC x-position
/// @param y     NDC y-position
/// @param text  LaTeX string (default: "#sqrt{s} = 14 TeV")
/// @param color Text colour (default: black = 1)
inline void ATLASEnergyLabel(Double_t x, Double_t y,
                             const char* text = "#sqrt{s} = 14 TeV HL-LHC",
                             Color_t color = 1) {
  TLatex p;
  p.SetNDC();
  p.SetTextFont(42);
  p.SetTextSize(0.038);  // slightly smaller than the main ATLAS label
  p.SetTextColor(color);
  p.DrawLatex(x, y, text);
}
