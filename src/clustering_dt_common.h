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

}

#endif  // CLUSTERING_DT_COMMON_H
