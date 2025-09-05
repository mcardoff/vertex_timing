#ifndef CLUSTERING_CONSTANTS_H
#define CLUSTERING_CONSTANTS_H

#include "clustering_includes.h"

namespace MyUtl {
  const bool DEBUG = false  ;
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
  const std::vector<Color_t> COLORS = {C01, C02, C03, C04, C05, C06, C07, C08, C09, C10, C11};

  // cut variables
  const int    MIN_JETS          = 1;    // min number of jets
  const int    MIN_PASSPT_JETS   = 0;    // min number of jets required to be >30 GeV
  const int    MIN_PASSETA_JETS  = 1;    // min number of forward jets required to be >30 GeV
  const int    MIN_NHS_TRACK     = 2;    // testing only
  const int    MAX_NHS_TRACK     = 6;    // testing only
  const double MIN_JETPT         = 30.0; // self explanatory
  const double MIN_ABS_ETA_JET   = 2.38; // min eta for a jet to be considered "forward"
  const double MIN_ABS_ETA_TRACK = 2.38; // min eta for a track to be considered "forward"
  const double MIN_TRACK_PT      = 1.0;  // track_pt > 1.0 GeV
  const double MAX_TRACK_PT      = 30.0; // track_pt < 30.0 GeV
  const double MAX_VTX_DZ        = 2.0;  // max error for reco HS vertex z
  const double MAX_NSIGMA        = 3.0;  // how close a track can be to PV
  const double MAX_TRK_VTX_SIG   = 3.0;  // How close a track can be to another vertex
  const double PILEUP_SMEAR = 175.0;

  const double MIN_HGTD_ETA = 2.38;
  const double MAX_HGTD_ETA = 4.0;

  const double DIFF_MIN = -1000.0, DIFF_MAX = 1000.0;
  const double DIFF_WIDTH = 2.0;

  const double PURITY_MIN = 0, PURITY_MAX = 1;
  const double PURITY_WIDTH = 0.05;

  const double FJET_MIN = 0, FJET_MAX = 31, FOLD_FJET = 5;
  const double FJET_WIDTH = 1.0;

  const double TRACK_MIN = 0, TRACK_MAX = 100, FOLD_TRACK = 20;
  const double TRACK_WIDTH = 1.0;

  const double PU_TRACK_MIN = TRACK_MIN, PU_TRACK_MAX = TRACK_MAX, FOLD_HS_TRACK = FOLD_TRACK;
  const double PU_TRACK_WIDTH = TRACK_WIDTH;

  const double HS_TRACK_MIN = TRACK_MIN, HS_TRACK_MAX = TRACK_MAX, FOLD_PU_TRACK = FOLD_TRACK;
  const double HS_TRACK_WIDTH = TRACK_WIDTH;

  const double PU_FRAC_MIN = 0, PU_FRAC_MAX = 1 , FOLD_PU_FRAC = PU_FRAC_MAX;
  const double PU_FRAC_WIDTH = 0.05;

  const double Z_MIN = -200, Z_MAX = 200, FOLD_Z = 100;
  const double Z_WIDTH = 10.0;

  const double EFF_YMIN = 0.0, EFF_YMAX = 1.5;
  const double PUR_YMIN = 0.0, PUR_YMAX = 1.5;
  const double RES_YMIN = 0.0, RES_YMAX = 40;

  /// enums and utilities for them
  enum Score {
    HGTD   = 0, TRKPT  = 1, TRKPTZ = 2, PASS   = 3,
    CALO90 = 4, CALO60 = 5, JUST90 = 6, JUST60 = 7,
    FILT90 = 8, FILT60 = 9, FILTJET = 10,
  };

  const std::vector<Score> ENUM_VEC = {
    Score::HGTD  ,
    Score::TRKPT , Score::TRKPTZ,
    // Score::CALO90, Score::CALO60,
    // Score::JUST60, Score::JUST90, Score::PASS,
    // Score::FILT90, Score::FILT60, Score::FILTJET,
    
  }; // valid values

  auto toString(
    MyUtl::Score score
  ) -> const char* {
    switch (score) {
    case Score::HGTD:    return "HGTD Algorithm";
    case Score::TRKPTZ:  return "Track p_{T}exp(-|#Deltaz|)";
    case Score::TRKPT:   return "Track p_{T}";
    case Score::CALO90:  return "Calo Time Exclusion (90 ps)";
    case Score::CALO60:  return "Calo Time Exclusion (60 ps)";
    case Score::JUST90:  return "Calo Time (90 ps)";
    case Score::JUST60:  return "Calo Time (60 ps)";
    case Score::PASS:    return "Pass Cluster";
    case Score::FILT90:  return "90 ps Track Filter";
    case Score::FILT60:  return "60 ps Track Filter";
    case Score::FILTJET: return "Only Tracks in Jets";
    default:             return "INVALID";
    }
  }

  enum FitParamFields {
    MEAN = 0, SIGMA = 1, CORE = 2, BKG = 3, RATIO = 4
  };
  const std::vector<FitParamFields> FITPARAM_VEC =
    { FitParamFields::SIGMA, FitParamFields::CORE, FitParamFields::BKG, FitParamFields::RATIO, FitParamFields::MEAN };

  auto toString(
    FitParamFields key
  ) -> const char* {
    switch (key) {
    case FitParamFields::SIGMA: return "SIGMA";
    case FitParamFields::CORE:  return "CORE";
    case FitParamFields::BKG:   return "BKG";
    case FitParamFields::RATIO: return "RATIO";
    case FitParamFields::MEAN:  return "MEAN";
    default:                    return "INVALID";
    }
  }

  template <typename T>
  auto folded(T raw, T fold) -> T {
    return (raw >= fold) ? fold : raw;
  }
  
}
#endif // CLUSTERING_CONSTANTS_H
