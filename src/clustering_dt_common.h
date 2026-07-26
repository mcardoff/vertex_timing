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
  // WAVeS fix variants — each isolates one change vs Score::WAVES (see the Score
  // definitions in clustering_constants.h). Compared in the waves_fixes_* group.
  m.emplace(Score::WAVES_LOJO,   AnalysisObj(label, Score::WAVES_LOJO  ));
  m.emplace(Score::WAVES_JETCAP, AnalysisObj(label, Score::WAVES_JETCAP));
  m.emplace(Score::WAVES_KERNEL, AnalysisObj(label, Score::WAVES_KERNEL));

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

  // --- WAVeS fix comparison -------------------------------------------------
  // The three self-selection fixes against the baseline they modify. Each variant
  // differs from Score::WAVES by exactly one change (LOJO: a jet cannot vote via
  // its own ghost tracks; JETCAP: nearJetPt capped; KERNEL: bounded ΔR kernel in
  // place of 1/ΔR), and all share WAVES's in-jet time/purity, so any separation
  // here is attributable to the selection change alone.
  // TRKPTZ is kept as the common reference line.
  moneyPlot(MyUtl::plotFilePath("comparisons", TString::Format("waves_fixes_%s.pdf", key).Data()).c_str(), key, canvas,
	    {
	      &mapHGTD.at(Score::TRKPTZ),
	      &mapHGTD.at(Score::WAVES),
	      &mapHGTD.at(Score::WAVES_LOJO),
	      &mapHGTD.at(Score::WAVES_JETCAP),
	      &mapHGTD.at(Score::WAVES_KERNEL),
	    },
	    {C02, C03, C05, C07, C08});

  // Same fixes measured against the achievable ceiling (WAVES_MISAS oracle), so the
  // plot answers "how much of the gap to perfect timing does each fix actually close?"
  moneyPlot(MyUtl::plotFilePath("comparisons", TString::Format("waves_fixes_oracle_%s.pdf", key).Data()).c_str(), key, canvas,
	    {
	      &mapHGTD.at(Score::WAVES),
	      &mapHGTD.at(Score::WAVES_LOJO),
	      &mapHGTD.at(Score::WAVES_JETCAP),
	      &mapHGTD.at(Score::WAVES_KERNEL),
	      &mapHGTD.at(Score::WAVES_MISAS),
	    },
	    {C03, C05, C07, C08, C04});
}

#endif  // CLUSTERING_DT_COMMON_H
