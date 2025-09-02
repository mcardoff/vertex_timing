#ifndef CLUSTERING_CONSTANTS_H
#define CLUSTERING_CONSTANTS_H

#include "clustering_includes.h"

namespace myutl {
  static bool debug = false  ;
  static int c01 = kP10Blue  ;
  static int c02 = kP10Red   ;
  static int c03 = kP10Yellow;
  static int c04 = kP10Gray  ;
  static int c05 = kP10Violet;
  static int c06 = kP10Brown ;
  static int c07 = kP10Orange;
  static int c08 = kP10Green ;
  static int c09 = kP10Ash   ;
  static int c10 = kP10Cyan  ;
  static std::vector<int> colors = {c01, c02, c03, c04, c05, c06, c07, c08, c09, c10};

  // cut variables
  static int    min_jets          = 1;    // min number of jets
  static int    min_passpt_jets   = 0;    // min number of jets required to be >30 GeV
  static int    min_passeta_jets  = 1;    // min number of forward jets required to be >30 GeV
  static int    min_nhs_track     = 2;    // testing only
  static int    max_nhs_track     = 6;    // testing only
  static double min_jetpt         = 30.0; // self explanatory
  static double min_abs_eta_jet   = 2.38; // min eta for a jet to be considered "forward"
  static double min_abs_eta_track = 2.38; // min eta for a track to be considered "forward"
  static double min_track_pt      = 1.0;  // track_pt > 1.0 GeV
  static double max_track_pt      = 30.0; // track_pt < 30.0 GeV
  static double max_vtx_dz        = 2.0;  // max error for reco HS vertex z
  static double max_nsigma        = 3.0;  // how close a track can be to PV
  static double max_trk_vtx_sig   = 3.0;  // How close a track can be to another vertex

  const double min_hgtd_eta = 2.38, max_hgtd_eta = 4.0;

  const double diff_min = -1000.0, diff_max = 1000.0;
  const double diff_width = 2.0;

  const double purity_min = 0, purity_max = 1;
  const double purity_width = 0.05;

  const double fjet_min = 0, fjet_max = 31, fold_fjet = 5;
  const double fjet_width = 1.0;

  const double track_min = 0, track_max = 100, fold_track = 20;
  const double track_width = 1.0;

  const double pu_track_min = track_min, pu_track_max = track_max, fold_hs_track = fold_track;
  const double pu_track_width = track_width;

  const double hs_track_min = track_min, hs_track_max = track_max, fold_pu_track = fold_track;
  const double hs_track_width = track_width;

  const double pu_frac_min = 0, pu_frac_max = 1 , fold_pu_frac = pu_frac_max;
  const double pu_frac_width = 0.05;

  const double z_min = -200, z_max = 200, fold_z = 100;
  const double z_width = 10.0;

  /// enums and utilities for them
  enum ScoreType {
    HGTD = 0, TRKPT = 1, TRKPTZ = 2, SUMPT2 = 3, DZ = 4,
    PASS = 5, IDEAL = 6, MAXHS = 7, JETPTDR = 8, TRKPTDR = 9,
    INVALID = -99 };

  static std::vector<ScoreType> enum_vec = {
    HGTD, TRKPT, TRKPTZ, PASS, //MAXHS, 
  }; // valid values

  static const char* toString(
    myutl::ScoreType score
  ) {
    switch (score) {
    case ScoreType::HGTD:     return "HGTD Algorithm";
    case ScoreType::TRKPTZ:   return "Track p_{T}exp(-|#Deltaz|)";
    case ScoreType::DZ:       return "exp(-|#Deltaz|)";
    case ScoreType::SUMPT2:   return "Sum(p_{T})^{2}";
    case ScoreType::TRKPT:    return "Track p_{T}";
    // case ScoreType::MAXPT:    return "Cluster with max p_{T}";
    case ScoreType::IDEAL:    return "Ideal cluster choice";
    case ScoreType::PASS:     return "Pass Cluster";
    case ScoreType::MAXHS:    return "Max HS Track Cluster";
    case ScoreType::JETPTDR:  return "exp(-min dR)*Jet p_{T}";
    case ScoreType::TRKPTDR:  return "exp(-min dR)*Track p_{T}";
    default:                  return "INVALID";
    }
  }

  enum FitParamFields { MEAN = 0, SIGMA = 1, CORE = 2, BKG = 3, RATIO = 4 };
  static std::vector<FitParamFields> fitparam_vec = { SIGMA, CORE, BKG, RATIO, MEAN };

  static const char* toString(
    FitParamFields key
  ) {
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
  static T folded(T raw, T fold) {
    return (raw >= fold) ? fold : raw;
  }
  
}
#endif // CLUSTERING_CONSTANTS_H
