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
  m.emplace(Score::TRKPTZ_TZJ,  AnalysisObj(label, Score::TRKPTZ_TZJ));
  m.emplace(Score::TRKPTZ_TZQ,  AnalysisObj(label, Score::TRKPTZ_TZQ));
  // Failure-mode fixes on WAVeS (src/failure_fixes.h): each alone, then the
  // cumulative chain in the failure ladder's order.
  for (const Score* s : {&Score::WFIX_RETIME, &Score::WFIX_MISTIMED, &Score::WFIX_SELECT,
                         &Score::WFIX_MISCL, &Score::WFIX_TEFF, &Score::WFIX_C3,
                         &Score::WFIX_C4, &Score::WFIX_C5, &Score::WFIX_C6,
                         &Score::WSEL_OK, &Score::WSEL_OK_MISCL,
                         &Score::WLOO_ALL, &Score::WLOO_MISTIMED, &Score::WLOO_SELECT, &Score::WLOO_MISCL, &Score::WLOO_TEFF})
    m.emplace(*s, AnalysisObj(label, *s));

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
  // Failure-mode fixes (src/failure_fixes.h). One file per failure mode: TRKPTZ
  // and WAVeS, then WAVeS with that one mode fixed -- the gap to WAVeS is what
  // the mode costs. Then the cumulative chain, in the failure ladder's order.
  const std::vector<std::pair<const char*, const Score*>> FIXES = {
    {"retime", &Score::WFIX_RETIME}, {"mistimed", &Score::WFIX_MISTIMED}, {"selection", &Score::WFIX_SELECT},
    {"miscl", &Score::WFIX_MISCL}, {"teff", &Score::WFIX_TEFF}};
  if (mapHGTD.count(Score::WFIX_RETIME)) {
    for (const auto& [name, sc] : FIXES)
      moneyPlot(MyUtl::plotFilePath("fixes", TString::Format("fix_%s_%s.pdf", name, key).Data()).c_str(), key, canvas,
                {&mapHGTD.at(Score::TRKPTZ), &mapHGTD.at(Score::WAVES), &mapHGTD.at(*sc)}, {C02, C03, C05});
    moneyPlot(MyUtl::plotFilePath("fixes", TString::Format("fix_chain_%s.pdf", key).Data()).c_str(), key, canvas,
              {&mapHGTD.at(Score::WAVES), &mapHGTD.at(Score::WFIX_MISTIMED), &mapHGTD.at(Score::WFIX_C3),
               &mapHGTD.at(Score::WFIX_C4), &mapHGTD.at(Score::WFIX_C5), &mapHGTD.at(Score::WFIX_C6)},
              {C03, C06, C01, C02, C04, C09});
    // The talk's version (slide 13, user 2026-11-11): the differential gain of
    // each fix, WAVeS then the chain through HS track efficiency, core
    // fraction from 0. The gap between two curves is that mode's contribution.
    // Slide 14 (user, 2026-11-12): P(good time | everything fixed but x), all
    // four fixed as the reference; n_HS >= 3, both pages (core width, fraction).
    // Named from "perfect" (user, 2026-11-12): each row is perfect with ONE failure mode left real; WAVeS as is.
    moneyPlot(MyUtl::plotFilePath("talk", TString::Format("talk_fix_allbut_%s.pdf", key).Data()).c_str(), key, canvas,
              {&mapHGTD.at(Score::WAVES), &mapHGTD.at(Score::WLOO_ALL), &mapHGTD.at(Score::WLOO_MISTIMED),
               &mapHGTD.at(Score::WLOO_SELECT), &mapHGTD.at(Score::WLOO_MISCL), &mapHGTD.at(Score::WLOO_TEFF)},
              {C03, kBlack, C06, C01, C02, C04}, -1.0, -1.0, 2.5,
              {"WAVeS", "Perfect (all four fixed with truth)", "Perfect + real HS track times",
               "Perfect + real cluster choice", "Perfect + real clustering", "Perfect + real HS track-time efficiency"});
    // Perfect timing / selection / both (a reviewer, 2026-11-13: "the plot where you were considering perfect timing
    // for tracks, and perfect selection"; the user: "a core fraction plot, with just 3 cases"): WAVeS as is, the
    // misassigned HS track times removed, the HS cluster chosen, and both. Page 2 is the core fraction.
    moneyPlot(MyUtl::plotFilePath("talk", TString::Format("talk_perfect_%s.pdf", key).Data()).c_str(), key, canvas,
              {&mapHGTD.at(Score::WAVES), &mapHGTD.at(Score::WFIX_MISTIMED), &mapHGTD.at(Score::WFIX_SELECT),
               &mapHGTD.at(Score::WFIX_C3)},
              {C03, C06, C01, kBlack}, -1.0, -1.0, 2.5,
              {"WAVeS", "Perfect timing", "Perfect selection", "Perfect timing + selection"});   // defined on the slide
    // ... and the same as the gap itself (user, 2026-11-12): each failure's cost,
    // read straight off the plot; WAVeS's whole gap in grey.
    gapPlot(MyUtl::plotFilePath("talk", TString::Format("talk_fix_cost_%s.pdf", key).Data()).c_str(), key, canvas,
            &mapHGTD.at(Score::WLOO_ALL),
            {&mapHGTD.at(Score::WLOO_MISTIMED), &mapHGTD.at(Score::WLOO_SELECT), &mapHGTD.at(Score::WLOO_MISCL),
             &mapHGTD.at(Score::WLOO_TEFF)},
            {C06, C01, C02, C04},
            {"misassigned HS track times", "the wrong cluster chosen", "pileup in the cluster", "untimed HS tracks"},
            &mapHGTD.at(Score::WAVES), "WAVeS: all of them, as now", 2.5, /*zoomCoreFraction=*/true);
    // Slide 14 before (user, 2026-11-11): each fix ALONE on WAVeS, to see which are
    // worth fixing at reco level first; n_HS >= 3 as the June talk.
    moneyPlot(MyUtl::plotFilePath("talk", TString::Format("talk_fix_alone_%s.pdf", key).Data()).c_str(), key, canvas,
              {&mapHGTD.at(Score::WAVES), &mapHGTD.at(Score::WFIX_MISTIMED), &mapHGTD.at(Score::WFIX_SELECT),
               &mapHGTD.at(Score::WFIX_MISCL), &mapHGTD.at(Score::WFIX_TEFF)},
              {C03, C06, C01, C02, C04}, -1.0, -1.0, 2.5);
    moneyPlot(MyUtl::plotFilePath("talk", TString::Format("talk_fix_chain_%s.pdf", key).Data()).c_str(), key, canvas,
              {&mapHGTD.at(Score::WAVES), &mapHGTD.at(Score::WFIX_MISTIMED), &mapHGTD.at(Score::WFIX_C3),
               &mapHGTD.at(Score::WFIX_C4), &mapHGTD.at(Score::WFIX_C5)},
              {C03, C06, C01, C02, C04}, 0.0, 1.3);
  }

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

  // Athena's t0 alone (the talk's "vertex time today" slide, user 2026-11-12): core
  // width and core fraction against the forward HS tracks, n_HS >= 3.
  moneyPlot(MyUtl::plotFilePath("talk", TString::Format("talk_athena_%s.pdf", key).Data()).c_str(), key, canvas,
	    {&mapHGTD.at(Score::HGTD)}, {C01}, -1.0, -1.0, 2.5);

  // The talk's t0s (the Global Performance talk, slides 8/9): WAVeS against
  // HGTD (Athena) -- WAVeS is the one algorithm presented (user, 2026-11-11) --
  // with the core-fraction page down to 0 -- Z+jets sits at
  // 0.2-0.7 below 5 HS tracks, under the default floor.
  moneyPlot(MyUtl::plotFilePath("talk", TString::Format("talk_t0_%s.pdf", key).Data()).c_str(), key, canvas,
	    {&mapHGTD.at(Score::HGTD), &mapHGTD.at(Score::WAVES)},
	    {C01, C03}, -1.0, -1.0, 2.5);   // n_HS >= 3 as in the June GPER talk (user, 2026-11-11); y range from the points

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
	      &mapHGTD.at(Score::TRKPTZ_TZQ),
	    },
	    {C02, C08, C01, C03, C04});

}

#endif  // CLUSTERING_DT_COMMON_H
