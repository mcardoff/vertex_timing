#ifndef ASSOC_STUDY_COMMON_H
#define ASSOC_STUDY_COMMON_H

// ---------------------------------------------------------------------------
// assoc_study_common.h
//   Shared between assoc_study_hist.cxx (event loop, writes histograms to a
//   ROOT file) and assoc_study_plot.cxx (reads that file, produces the report
//   and PDFs) -- the same hist/plot split the two primary analyses use.
//
//   WHAT THIS STUDY ASKS
//   The clustering has two candidate rules for deciding which tracks belong to
//   the hard-scatter vertex in z, and therefore which tracks t0 is built from:
//
//     SIGNIFICANCE  |z0 - z_vtx| / sqrt(1.15^2 * var_z0) < N
//                   the rule this analysis has always used, at N = 3.0.
//     DZ_PARA       |z0 - z_vtx| / getNewDzpara(eta, pT)  < S
//                   the reference R_pT study's empirical parameterization of
//                   the observed z0 spread, at S = 1.4 there.
//
//   They are not on a common scale, so "which is tighter" is not even
//   well-posed globally -- getNewDzpara is looser than 3σ for some (eta, pT)
//   and tighter for others. The only way to answer "which one produces a
//   better vertex t0" is to measure the t0 both ways.
//
//   METRICS (per rule, for BOTH selection scores)
//     inclusive core σ(Δt)     — ACTIVE_FIT_MODEL fit to Δt = t_sel - t_truth
//     inclusive core fraction  — fraction of events with |Δt| < PASS_SIGMA
//     core σ vs n forward HS tracks / n forward jets
//     core fraction vs the same two
//
//   Both scores are reported because they consume the track list differently:
//   TRKPTZ weights by pT and |Δz| only, so it sees the association directly;
//   WAVeS additionally re-weights by jet proximity and recomputes its time
//   from in-jet tracks alone, which can mask or amplify what the association
//   did. A rule that helps one need not help the other.
//
//   WHY ONE EVENT LOOP RATHER THAN ONE RUN PER RULE
//   Every rule is evaluated against the SAME event, so the comparison is
//   paired: the rules see identical events, identical jets, identical truth,
//   and the per-bin populations of the differential plots are common by
//   construction (see COUNTING_NSIGMA in clustering_constants.h). It also
//   costs one pass over a 52 GB ntuple instead of seven.
// ---------------------------------------------------------------------------

#include "clustering_constants.h"
#include "plotting_utilities.h"

#include <string>
#include <vector>

using namespace MyUtl;

// ---------------------------------------------------------------------------
// The rules under study.
//
//   The significance family is scanned as well as the parameterized one, and
//   that is not padding. The two rules differ in TWO ways at once -- their
//   functional form AND how tight they happen to be at their nominal working
//   points -- so a single point from each cannot separate "the parameterization
//   is a better shape" from "the parameterization is simply tighter here".
//   Scanning both families lets rules be compared at MATCHED track purity,
//   which is the only comparison that isolates the shape.
//
//   The parameterization is scanned around the reference study's 1.4 working
//   point, extending well below it because 1.4 was tuned for a different
//   objective -- jet-level R_pT pileup rejection, where shedding pileup tracks
//   is worth losing hard-scatter ones -- and there is no reason that optimum
//   should transfer to estimating a vertex time, where lost HS tracks directly
//   degrade the average t0 is built from.
//
//   `tag` doubles as the AnalysisObj label, so it must be filename-safe and
//   unique: AnalysisObj lowercases it and strips spaces/dots/plus to build
//   every histogram name, and the whole study writes into one ROOT file whose
//   keys are those names.
// ---------------------------------------------------------------------------
inline const std::vector<AssocRule>& assocRules() {
  static const std::vector<AssocRule> rules = {
    { AssocRule::Kind::SIGNIFICANCE, 3.0, "signif30", "|#Deltaz_{0}|/#sigma_{z_{0}} < 3.0"  },
    { AssocRule::Kind::SIGNIFICANCE, 2.5, "signif25", "|#Deltaz_{0}|/#sigma_{z_{0}} < 2.5"  },
    { AssocRule::Kind::SIGNIFICANCE, 2.0, "signif20", "|#Deltaz_{0}|/#sigma_{z_{0}} < 2.0"  },
    { AssocRule::Kind::SIGNIFICANCE, 1.5, "signif15", "|#Deltaz_{0}|/#sigma_{z_{0}} < 1.5"  },
    { AssocRule::Kind::DZ_PARA,      0.6, "dzpara06", "|#Deltaz_{0}|/#deltaz(#eta,p_{T}) < 0.6" },
    { AssocRule::Kind::DZ_PARA,      0.8, "dzpara08", "|#Deltaz_{0}|/#deltaz(#eta,p_{T}) < 0.8" },
    { AssocRule::Kind::DZ_PARA,      1.0, "dzpara10", "|#Deltaz_{0}|/#deltaz(#eta,p_{T}) < 1.0" },
    { AssocRule::Kind::DZ_PARA,      1.2, "dzpara12", "|#Deltaz_{0}|/#deltaz(#eta,p_{T}) < 1.2" },
    { AssocRule::Kind::DZ_PARA,      1.4, "dzpara14", "|#Deltaz_{0}|/#deltaz(#eta,p_{T}) < 1.4" },
    { AssocRule::Kind::DZ_PARA,      1.6, "dzpara16", "|#Deltaz_{0}|/#deltaz(#eta,p_{T}) < 1.6" },
    { AssocRule::Kind::DZ_PARA,      2.0, "dzpara20", "|#Deltaz_{0}|/#deltaz(#eta,p_{T}) < 2.0" },

    // Union with ghost association: accept a track if it passes the z-test OR
    // is ghost-associated to a qualifying forward jet. Tests the premise that
    // ghost association picks up tracks carrying hard-scatter information the
    // z-cut misses -- a WIDER list, not a tighter one, and so the opposite
    // direction from everything above.
    //
    // Two base rules: the incumbent, and the best single rule for TRKPTZ, so
    // the union can be read against both a loose and a tight starting point.
    //
    // What the track-level numbers say to expect (util/scratch/
    // ghost_assoc_diag.cxx): the added tracks are essentially time-random --
    // mean |Delta t| to the truth HS vertex of 224 ps against 229 ps for a
    // track picked at random from the event -- so most of them cannot be
    // absorbed into the HS cluster and instead form their own. But 24% of
    // them do land inside the clustering window, which takes the pileup
    // fraction of the tracks actually inside the HS cluster from 29% to 44%
    // (dzpara 1.0 base). Whether that costs core fraction is what these rows
    // measure; a track-level projection cannot answer it.
    { AssocRule::Kind::SIGNIFICANCE, 3.0, "signif30ghost",
      "|#Deltaz_{0}|/#sigma_{z_{0}} < 3.0 #cup ghost", true },
    { AssocRule::Kind::DZ_PARA,      1.0, "dzpara10ghost",
      "|#Deltaz_{0}|/#deltaz(#eta,p_{T}) < 1.0 #cup ghost", true },
  };
  return rules;
}

// Plain-ASCII rule name for console tables and the markdown report, where the
// TLatex `label` above would be unreadable.
inline std::string ruleAscii(const AssocRule& r) {
  char buf[64];
  if (r.kind == AssocRule::Kind::SIGNIFICANCE)
    std::snprintf(buf, sizeof buf, "signif < %.1f", r.cut);
  else
    std::snprintf(buf, sizeof buf, "dzpara x %.1f", r.cut);
  return std::string(buf) + (r.orGhost ? " OR ghost" : "");
}

// The two selection scores reported for every rule. TRKPTZ is the baseline
// algorithm; WAVES is the jet-proximity-weighted alternative currently under
// study. Neither builds its own cluster collection (usesOwnCollection is false
// for both), so a rule's analysis map needs no auxiliary collections and the
// event loop stays cheap enough to run seven of them per event.
inline const std::vector<Score>& assocScores() {
  static const std::vector<Score> scores = { Score::TRKPTZ, Score::WAVES };
  return scores;
}

// ---------------------------------------------------------------------------
// buildAssocMap
//   One analysis map per rule. Called by BOTH stages -- the plotting stage
//   re-runs this same booking to get correctly-named, correctly-binned empty
//   histograms and then loads the saved contents into them, which is the
//   "construct-empty-then-load" pattern the other split analyses use (and
//   doubles as a binning-mismatch guard via HistReader::LoadInto).
// ---------------------------------------------------------------------------
inline auto buildAssocMap(const AssocRule& rule) -> std::map<Score, AnalysisObj> {
  std::map<Score, AnalysisObj> m;
  for (const auto& s : assocScores())
    m.emplace(s, AnalysisObj(rule.tag.c_str(), s));
  return m;
}

// Base name of the histogram file the two stages exchange. Goes through
// MyUtl::histFilePath() so it follows the same sample-prefixed convention as
// clustering_hist.root / rpt_v5_hist.root.
inline const char* ASSOC_HIST_BASENAME = "assoc_study_hist.root";

// Scalar keys written by the hist stage and read back by the plot stage: the
// per-rule track accounting that gives the resolution numbers their context
// (how many tracks a rule keeps, and how clean they are).
inline std::string assocScalarKey(const AssocRule& r, const char* field) {
  return "meta_" + r.tag + "_" + field;
}

#endif  // ASSOC_STUDY_COMMON_H
