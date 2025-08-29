#ifndef COMMON_INCLUDES_H
#define COMMON_INCLUDES_H

// ROOT Headers (ordered alphabetically or by common usage)
#include <Rtypes.h>
#include <TRandom1.h> // Only if used by multiple files, otherwise keep in specific .cxx
#include <TBranch.h>
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
#include <TLine.h> // Only if used by multiple files, otherwise keep in specific .cxx
#include <TProfile.h>
#include <TMath.h>
#include <TTreeReader.h>
#include <TTreeReaderArray.h>
#include <TString.h>
#include <TStyle.h>


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

// Project-Specific Headers (always include after common ones)
#include "./clustering_utilities.h" // Assuming my_utilities.h is in the same directory or adjust path

// You might also want to centralize #define debug false if it's consistently used
// or move it to a global configuration. For now, let's keep it here if it's always false.

#endif // COMMON_INCLUDES_H
