#ifndef CLUSTERING_CONSTANTS_H
#define CLUSTERING_CONSTANTS_H

// ---------------------------------------------------------------------------
// clustering_constants.h
//   Central repository for all compile-time constants, enums, and small
//   utility functions shared across the analysis.  Everything lives inside
//   the MyUtl namespace so that name collisions with ROOT globals are
//   avoided.  Sections:
//     1. Debug flag
//     2. Plot colour palette
//     3. Event / track selection cuts
//     4. Histogram axis ranges and fold values
//     5. Score enum + string converters
//     6. FitParamFields enum + string converter
//     7. folded() helper
// ---------------------------------------------------------------------------

#include "clustering_includes.h"

namespace MyUtl {

// ---------------------------------------------------------------------------
// 1. Debug flag
//   Set to true to enable verbose per-event stdout logging throughout the
//   analysis.  Compiled out entirely when false (no runtime cost).
// ---------------------------------------------------------------------------
  const bool DEBUG = false  ;
// ---------------------------------------------------------------------------
// 2. Plot colour palette
//   Eleven colours drawn from ROOT's P10 and P6 palettes.  COLORS is the
//   ordered vector used when cycling through scores or histogram slices.
// ---------------------------------------------------------------------------
  const Color_t C01 = kP10Blue  ;
  const Color_t C02 = kP10Red   ;
  const Color_t C03 = kP10Yellow;
  const Color_t C04 = kP10Gray  ;
  const Color_t C05 = kP10Violet;
  const Color_t C06 = kP10Brown ;
  const Color_t C07 = kP10Orange;
  const Color_t C08 = kP10Green ;
  const Color_t C09 = kP10Ash   ;
  const Color_t C10 = kP10Cyan  ;
  const Color_t C11 = kP6Blue   ;
  const std::vector<Color_t> COLORS = {C01, C02, C03, C04, C05, C06,
                                       C07, C08, C09, C10, C11};  

// ---------------------------------------------------------------------------
// 3. Event / track selection cuts
//   All kinematic and association thresholds used by the event selection
//   and track clustering pipeline.  Changing a cut here propagates
//   automatically to every function that references these constants.
// ---------------------------------------------------------------------------
  // cut variables
  const int    MIN_JETS           = 2;     // min n jets
  const int    MIN_PASSPT_JETS    = 2;     // min n jets >30 GeV
  const int    MIN_PASSETA_JETS   = 1;     // min n forward jets >30GeV
  const int    MIN_NHS_TRACK      = 2;     // testing only
  const int    MAX_NHS_TRACK      = 6;     // testing only
  const double VBS_JET_D_ETA      = 3.0;   // min eta separation for VBS Jets
  const double MIN_JETPT          = 30.0;  // self explanatory
  const double MIN_ABS_ETA_JET    = 2.38;  // min eta for a "forward" jet
  const double MAX_ABS_ETA_JET    = 4.00;  // max eta for a "forward" jet
  const double MIN_ABS_ETA_TRACK  = 2.38;  // min eta for a "forward" track 
  const double MAX_ABS_ETA_TRACK  = 4.00;  // max eta for a "forward" track
  const double MIN_TRACK_PT       = 1.0;   // clustered track_pt > 0.5 GeV
  const double MAX_TRACK_PT       = 30.0;  // clustered track_pt < 30.0 GeV
  const double MIN_TRACK_PT_COUNT = 1.0;   // track_pt > 1.0 GeV for counting purposes
  const double MAX_VTX_DZ         = 2.0;   // max error for reco HS vertex z
  const double MAX_NSIGMA         = 3.0;   // how close a track can be to PV
  const double MAX_TRK_VTX_SIG    = 3.0;   // Pileup removal sigma
  const double PASS_SIGMA         = 20.0;  // Pass threshold for efficiency
  const double PILEUP_SMEAR       = 175.0; // Pileup track resolution
  const double MIN_HGTD_ETA       = 2.38;  // Min eta of HGTD
  const double MAX_HGTD_ETA       = 4.0;   // Max eta of HGTD

// ---------------------------------------------------------------------------
// 4. Histogram axis ranges and fold values
//   xMin/xMax define the full histogram axis.  FOLD_* values mark the
//   point at which overflow is collapsed into the last visible bin so that
//   sparse high-multiplicity tails don't dominate the plots.
// ---------------------------------------------------------------------------
  // plotting ranges
  const double DIFF_MIN = -1000.0, DIFF_MAX = 1000.0;
  const double DIFF_WIDTH = 2.0;

  const double PURITY_MIN = 0, PURITY_MAX = 1;
  const double PURITY_WIDTH = 0.05;

  const double FJET_MIN = 0, FJET_MAX = 31, FOLD_FJET = 5;
  const double FJET_WIDTH = 1.0;

  const double VTX_DZ_MIN = 0, VTX_DZ_MAX = 5.0, FOLD_VTX_DZ = 2.0;
  const double VTX_DZ_WIDTH = 0.1;

  const double TRACK_MIN = 0, TRACK_MAX = 100, FOLD_TRACK = 25;
  const double TRACK_WIDTH = 1.0;

  const double PU_TRACK_MIN = TRACK_MIN, PU_TRACK_MAX = TRACK_MAX, FOLD_HS_TRACK = FOLD_TRACK;  
  const double PU_TRACK_WIDTH = TRACK_WIDTH;

  const double HS_TRACK_MIN = TRACK_MIN, HS_TRACK_MAX = TRACK_MAX, FOLD_PU_TRACK = FOLD_TRACK;
  const double HS_TRACK_WIDTH = TRACK_WIDTH;

  const double PU_FRAC_WIDTH = 0.05;
  const double PU_FRAC_MIN = 0, PU_FRAC_MAX = 1.0 + PU_FRAC_WIDTH, FOLD_PU_FRAC = 1.0;

  const double Z_MIN = -200, Z_MAX = 200, FOLD_Z = 100;
  const double Z_WIDTH = 10.0;

  const double EFF_YMIN = 0.0, EFF_YMAX = 1.5;
  const double PUR_YMIN = 0.0, PUR_YMAX = 1.5;
  const double RES_YMIN = 0.0, RES_YMAX = 40;

// ---------------------------------------------------------------------------
// 5. Score enum + string converters
//   Each enumerator identifies one cluster-selection algorithm.  ENUM_VEC
//   is the canonical iteration order used throughout the analysis.
//   toString()      returns a ROOT-LaTeX label suitable for plot titles.
//   toStringShort() returns a compact identifier used in histogram names
//                   and on-screen tables.
// ---------------------------------------------------------------------------
  /// enums and utilities for them
  enum Score {
    HGTD = 0,
    TRKPT = 1,
    TRKPTZ = 2,
    PASS = 3,
    CALO90 = 4,
    CALO60 = 5,
    JUST90 = 6,
    JUST60 = 7,
    FILT90 = 8,
    FILT60 = 9,
    FILTJET = 10,
    TESTML = 11,
    TEST_MISCL = 12,  // TESTML + purity cut (misclustering study)
  };

  const std::vector<Score> ENUM_VEC = {
    Score::HGTD, Score::PASS, Score::TRKPT, Score::TRKPTZ, Score::CALO60,
    Score::CALO90, Score::JUST60, Score::JUST90, Score::FILT90, Score::FILT60,
    Score::FILTJET, Score::TESTML, Score::TEST_MISCL,
  }; // valid values

  auto toString(
    MyUtl::Score score
  ) -> const char* {
    switch (score) {
    case Score::HGTD:    return "HGTD Algorithm";
    case Score::TRKPTZ:  return "Track p_{T} exp(-|#Delta z|)";
    case Score::TRKPT:   return "Track p_{T}";
    case Score::CALO90:  return "Calo Time Exclusion (90 ps)";
    case Score::CALO60:  return "Calo Time Exclusion (60 ps)";
    case Score::JUST90:  return "Calo Time (90 ps)";
    case Score::JUST60:  return "Calo Time (60 ps)";
    case Score::PASS:    return "Pass Cluster";
    case Score::FILT90:  return "90 ps Track Filter";
    case Score::FILT60:  return "60 ps Track Filter";
    case Score::FILTJET: return "Only Tracks in Jets";
    case Score::TESTML:     return "DNN";
    case Score::TEST_MISCL: return "TRKPTZ (pure clusters)";
    default:                return "INVALID";
    }
  }

  auto toStringShort(
    MyUtl::Score score
  ) -> const char* {
    switch (score) {
    case Score::HGTD:       return "HGTD";
    case Score::TRKPTZ:     return "TRKPTZ";
    case Score::TRKPT:      return "TRKPT";
    case Score::CALO90:     return "CALO90";
    case Score::CALO60:     return "CALO60";
    case Score::PASS:       return "PASS";
    case Score::JUST90:     return "JUST90";
    case Score::JUST60:     return "JUST60";
    case Score::FILT90:     return "FILT90";
    case Score::FILT60:     return "FILT60";
    case Score::FILTJET:    return "FILTJET";
    case Score::TESTML:     return "TESTML";
    case Score::TEST_MISCL: return "MISCL";
    default:                return "INVALID";
    }
  }

// ---------------------------------------------------------------------------
// 6. FitParamFields enum + string converter
//   Indexes the individual parameters extracted from a double-Gaussian fit:
//   MEAN, core SIGMA, background BSIGMA, amplitudes CORE / BKG, and their
//   RATIO.  Used by FitParams::fromEnum() to select which distribution to
//   draw or fill.
// ---------------------------------------------------------------------------
  enum FitParamFields {
    MEAN = 0, SIGMA = 1, CORE = 2, BKG = 3, RATIO = 4, BSIGMA = 5,
  };
  const std::vector<FitParamFields> FITPARAM_VEC =
    { FitParamFields::SIGMA, FitParamFields::BSIGMA, FitParamFields::CORE, FitParamFields::BKG, FitParamFields::RATIO, FitParamFields::MEAN };

  auto toString(
    FitParamFields key
  ) -> const char* {
    switch (key) {
    case FitParamFields::SIGMA:  return "SIGMA";
    case FitParamFields::BSIGMA: return "BSIGMA";
    case FitParamFields::CORE:   return "CORE";
    case FitParamFields::BKG:    return "BKG";
    case FitParamFields::RATIO:  return "RATIO";
    case FitParamFields::MEAN:   return "MEAN";
    default:                     return "INVALID";
    }
  }

// ---------------------------------------------------------------------------
// 7. folded() helper
//   Clamps raw to fold when raw >= fold, collapsing overflow into the last
//   visible histogram bin.  Used when filling efficiency/resolution plots
//   to avoid sparse tails at high multiplicity.
// ---------------------------------------------------------------------------
  template <typename T>
  auto folded(T raw, T fold) -> T {
    return (raw >= fold) ? fold : raw;
  }
  
}
#endif // CLUSTERING_CONSTANTS_H
