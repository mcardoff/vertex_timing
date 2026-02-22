#ifndef CLUSTERING_INCLUDES_H
#define CLUSTERING_INCLUDES_H

// ---------------------------------------------------------------------------
// clustering_includes.h
//   Single aggregation point for all external headers required by the
//   analysis.  Every other source file includes only this header (plus
//   the project-specific headers it depends on) so that third-party
//   include paths are managed in one place.  Sections:
//     1. ROOT headers  — framework types, I/O, histogramming, fitting,
//                        graphics, and TTree reading
//     2. Boost headers — filesystem iteration for ntuple directory scan
//     3. Standard C++ library headers
// ---------------------------------------------------------------------------

// ---------------------------------------------------------------------------
// 1. ROOT headers
// ---------------------------------------------------------------------------
#include <Rtypes.h>
#include <RtypesCore.h>
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
#include "TMVA/Tools.h"
#include "TMVA/Reader.h"
#include <TTreeReader.h>
#include <TTreeReaderArray.h>
#include <TString.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TVector2.h>
#include <TLorentzVector.h>

// ---------------------------------------------------------------------------
// 2. Boost headers
// ---------------------------------------------------------------------------
// Boost Headers
#include <boost/filesystem.hpp>

// ---------------------------------------------------------------------------
// 3. Standard C++ library headers
// ---------------------------------------------------------------------------
// Standard C++ Library Headers
#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <numeric>
#include <ostream>
#include <utility>
#include <vector>
#include <map>
#include <memory>

#endif // CLUSTERING_INCLUDES_H
