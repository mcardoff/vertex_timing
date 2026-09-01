#ifndef CLUSTERING_DT_COMMON_H
#define CLUSTERING_DT_COMMON_H

// ---------------------------------------------------------------------------
// clustering_dt_common.h
//   Shared between clustering_hist.cxx (event loop, writes histograms to a
//   ROOT file) and clustering_plot.cxx (reads that file, produces PDFs) --
//   split out of the former single clustering_dt.cxx so the plotting stage
//   can build a std::map<Score, AnalysisObj> (via buildAnalysisMap) and
//   reproduce the comparison plots (via makeComparisonPlots) without
//   depending on anything ntuple-related.
// ---------------------------------------------------------------------------

#include "clustering_constants.h"
#include "plotting_utilities.h"

using namespace MyUtl;

// ---------------------------------------------------------------------------
// Scenario
//   HGTD is the only scenario currently built/processed; IDEAL_RES/IDEAL_EFF
//   remain for future reuse (see buildAnalysisMap's branching below).
// ---------------------------------------------------------------------------
enum class Scenario { HGTD, IDEAL_RES, IDEAL_EFF };

// ---------------------------------------------------------------------------
// Helper: build one analysis map for a given scenario label.
//   The active scores are listed explicitly so adding/removing a score
//   only requires changing this one function.
// ---------------------------------------------------------------------------
inline auto buildAnalysisMap(
   Scenario scenario
) -> std::map<Score, AnalysisObj> {
  const char* label = [&]() -> const char* {
    switch (scenario) {
      case Scenario::HGTD:      return "HGTD Times";
      case Scenario::IDEAL_RES: return "Ideal Res.";
      case Scenario::IDEAL_EFF: return "Ideal Res.+Eff.";
    }
    return "";
  }();

  std::map<Score, AnalysisObj> m;

  // Scores active in all scenarios
  m.emplace(Score::TRKPTZ,     AnalysisObj(label, Score::TRKPTZ    ));
  m.emplace(Score::WAVES,   AnalysisObj(label, Score::WAVES  ));
  m.emplace(Score::JET_T_REFINED, AnalysisObj(label, Score::JET_T_REFINED));
  // m.emplace(Score::PASS,       AnalysisObj(label, Score::PASS      ));
  m.emplace(Score::TEST_MISAS, AnalysisObj(label, Score::TEST_MISAS));
  m.emplace(Score::WAVES_MISCL, AnalysisObj(label, Score::WAVES_MISCL));
  m.emplace(Score::WAVES_MISAS, AnalysisObj(label, Score::WAVES_MISAS));
  // VBS topology regions: WAVeS selection restricted to R1 / R2 events. Read
  // against the plain WAVES row above to see what the topology costs or buys,
  // with the algorithm held fixed.
  m.emplace(Score::VBF_R1,      AnalysisObj(label, Score::VBF_R1));
  m.emplace(Score::VBF_R2,      AnalysisObj(label, Score::VBF_R2));
  // TRKPTZ with the pT sum split on the per-track PV/PU proximity flag. Both
  // orientations, so the flag's sense is settled by the run rather than assumed.
  m.emplace(Score::TRKPTZ_PV,   AnalysisObj(label, Score::TRKPTZ_PV));
  m.emplace(Score::TRKPTZ_PU,   AnalysisObj(label, Score::TRKPTZ_PU));
  m.emplace(Score::TRKPTZ_PUW,  AnalysisObj(label, Score::TRKPTZ_PUW));
  m.emplace(Score::TRKPTZ_TZ,   AnalysisObj(label, Score::TRKPTZ_TZ));
  m.emplace(Score::TRKPTZ_TZ_IJ, AnalysisObj(label, Score::TRKPTZ_TZ_IJ));
  m.emplace(Score::TRKPTZ_TZ_OJ, AnalysisObj(label, Score::TRKPTZ_TZ_OJ));
  m.emplace(Score::TRKPTZ_TZ_GIJ, AnalysisObj(label, Score::TRKPTZ_TZ_GIJ));
  m.emplace(Score::WAVES_GIJ,   AnalysisObj(label, Score::WAVES_GIJ));

  // Scores active only in the real-HGTD scenario
  if (scenario == Scenario::HGTD) {
    m.emplace(Score::HGTD,       AnalysisObj(label, Score::HGTD      ));
    // m.emplace(Score::TEST_HS,    AnalysisObj(label, Score::TEST_HS    ));
    // m.emplace(Score::HGTD_SORT, AnalysisObj(label, Score::HGTD_SORT));
    // m.emplace(Score::CONE,       AnalysisObj(label, Score::CONE      ));
    // m.emplace(Score::CONE_BDT,  AnalysisObj(label, Score::CONE_BDT ));
    // m.emplace(Score::FILTJET,   AnalysisObj(label, Score::FILTJET  ));
  }

  return m;
}

// ---------------------------------------------------------------------------
// Helper: generate all per-KEY comparison plots for one variable KEY.
// ---------------------------------------------------------------------------
inline void makeComparisonPlots(
  const char* key,
  TCanvas* canvas,
  std::map<Score, AnalysisObj>& mapHGTD
) {
  moneyPlot(MyUtl::plotFilePath("comparisons", TString::Format("theplot_%s.pdf", key).Data()).c_str(), key, canvas,
	    {
	      &mapHGTD.at(Score::HGTD),
	      &mapHGTD.at(Score::TRKPTZ),
	      &mapHGTD.at(Score::TEST_MISAS),
	    });

  // WAVeS oracle comparison: selection by WAVeS score, gated like the TRKPTZ oracles.
  // Focus on the timing-misassignment oracle: HGTD, TRKPTZ, WAVeS, and WAVeS
  // with ideal timing (WAVES_MISAS).  Misclustering case removed.
  // Color override: WAVES yellow, ideal timing violet.
  moneyPlot(MyUtl::plotFilePath("comparisons", TString::Format("waves_oracle_%s.pdf", key).Data()).c_str(), key, canvas,
	    {
	      &mapHGTD.at(Score::HGTD),
	      &mapHGTD.at(Score::TRKPTZ),
	      &mapHGTD.at(Score::WAVES),
	      &mapHGTD.at(Score::WAVES_MISAS),
	    },
	    {C01, C02, C03, C04});

  // Simple three-way comparison: HGTD vs TRKPTZ vs WAVeS (no oracles).
  moneyPlot(MyUtl::plotFilePath("comparisons", TString::Format("hgtd_trkptz_waves_%s.pdf", key).Data()).c_str(), key, canvas,
	    {
	      &mapHGTD.at(Score::HGTD),
	      &mapHGTD.at(Score::TRKPTZ),
	      &mapHGTD.at(Score::WAVES),
	    },
	    {C01, C02, C03});

  // VBS topology regions: the same WAVeS algorithm measured inclusively, in R1
  // (both candidate legs forward), and in R2 (forward PU leg + central HS leg).
  // WAVES is included as the reference the two regions are read against -- the
  // algorithm is identical across all three curves, so any separation is the
  // topology talking, not the selector.
  moneyPlot(MyUtl::plotFilePath("comparisons", TString::Format("vbf_regions_%s.pdf", key).Data()).c_str(), key, canvas,
	    {
	      &mapHGTD.at(Score::WAVES),
	      &mapHGTD.at(Score::VBF_R1),
	      &mapHGTD.at(Score::VBF_R2),
	    },
	    {C03, C08, C02});

  // Per-track dz weighting against the baseline it modifies, with WAVeS as the
  // reference for what a jet-aware selector buys. TRKPTZ and TRKPTZ_TZ share
  // the same clusters, the same cluster-level exp(-1.5|dz|) envelope and the
  // same time estimator, and differ ONLY in the per-track pT weight -- so the
  // separation between those two curves is exactly what the per-track term is
  // worth, with nothing else moving.
  moneyPlot(MyUtl::plotFilePath("comparisons", TString::Format("trkptz_tz_%s.pdf", key).Data()).c_str(), key, canvas,
	    {
	      &mapHGTD.at(Score::TRKPTZ),
	      &mapHGTD.at(Score::TRKPTZ_TZ),
	      &mapHGTD.at(Score::WAVES),
	    },
	    {C02, C08, C03});

  // Where does the usable timing live? The same selector (TRKPTZ_TZ) reporting
  // the full-cluster time, the in-jet subset, and its complement. Selection is
  // bit-identical across all three -- updateScores copies one value into all
  // three score ids -- so any separation is the TIMING alone.
  moneyPlot(MyUtl::plotFilePath("comparisons", TString::Format("timesource_%s.pdf", key).Data()).c_str(), key, canvas,
	    {
	      &mapHGTD.at(Score::TRKPTZ_TZ),
	      &mapHGTD.at(Score::TRKPTZ_TZ_IJ),
	      &mapHGTD.at(Score::TRKPTZ_TZ_OJ),
	    },
	    {C08, C01, C02});

  // The deployable comparison: baseline, the per-track dz weighting, the same
  // plus the SAMPLE-INDEPENDENT guarded in-jet re-timing, and WAVeS as the
  // incumbent to beat. Every curve here is a single algorithm with no per-sample
  // configuration, so this is the set that could actually ship.
  moneyPlot(MyUtl::plotFilePath("comparisons", TString::Format("guarded_%s.pdf", key).Data()).c_str(), key, canvas,
	    {
	      &mapHGTD.at(Score::TRKPTZ),
	      &mapHGTD.at(Score::TRKPTZ_TZ),
	      &mapHGTD.at(Score::TRKPTZ_TZ_GIJ),
	      &mapHGTD.at(Score::WAVES),
	      &mapHGTD.at(Score::WAVES_GIJ),
	    },
	    {C02, C08, C01, C03, C04});

}

#endif  // CLUSTERING_DT_COMMON_H
