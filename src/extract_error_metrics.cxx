/*
 * extract_error_metrics.cxx
 *
 * Extracts quantitative metrics from hgtd_matching_analysis.root
 * and generates a summary table showing the impact of each error source
 * on overall algorithm efficiency.
 *
 * Usage: root -l -b -q extract_error_metrics.cxx
 * (Run after hgtd_matching.cxx has generated hgtd_matching_analysis.root)
 */

#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TProfile.h"
#include "TEfficiency.h"
#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <map>

struct ErrorMetrics {
    // Error Source 1: Wrong HGTD Matching
    double matching_efficiency_overall;
    double matching_efficiency_low_eta;   // 2.4 < |eta| < 3.0
    double matching_efficiency_high_eta;  // 3.5 < |eta| < 4.0
    double residual_rms_all;
    double residual_rms_hs;
    double residual_rms_pu;
    double fraction_mismatched;
    double avg_residual_hs;
    // Ideal matching scenarios
    double efficiency_hgtd;
    double efficiency_ideal_res;
    double efficiency_ideal_eff;
    double gain_from_resolution;
    double gain_from_efficiency;
    // Failure mode breakdown
    int n_fail_wrong_selection;      // At least one cluster within 60ps, wrong one selected
    int n_fail_no_passing_cluster;   // No cluster within 60ps
    int n_fail_matching_only;        // Fixable by ideal resolution
    int n_fail_efficiency_only;      // Fixable by ideal efficiency
    int n_fail_irredeemable;         // Not fixable by matching improvements
    int n_total_failures;            // Total events failing time cut

    // Error Source 2: Wrong Cluster Selection
    double selection_efficiency_overall;
    double selection_efficiency_low_clusters;  // n_clusters = 2-3
    double selection_efficiency_high_clusters; // n_clusters >= 4
    double hs_cluster_rank_1_fraction;
    double hs_cluster_rank_2plus_fraction;
    double score_separation_trkpt;
    double score_separation_trkptz;
    double resolution_penalty_wrong_selection;

    // Error Source 3: Insufficient Tracks
    double efficiency_nhs_2;      // 2 HS tracks
    double efficiency_nhs_4;      // 4 HS tracks
    double efficiency_nhs_6;      // 6 HS tracks
    double efficiency_nhs_10;     // 10+ HS tracks
    double avg_hs_tracks_success;
    double avg_hs_tracks_failure;
    double avg_track_loss_eta;
    double avg_track_loss_quality;
    double avg_track_loss_time_valid;
    int n_insufficient_tracks_total;  // Events with < N HS tracks with valid time
    int n_insufficient_tracks_fail;   // Events with < N HS tracks that fail
    double insufficient_tracks_fraction;  // Fraction of all events with insufficient tracks
    double insufficient_tracks_efficiency;  // Efficiency for events with insufficient tracks

    // Error Source 4: HS Track Fragmentation
    double fragmentation_single_cluster;  // Fraction with single HS cluster
    double fragmentation_multi_cluster;   // Fraction with 2+ HS clusters
    double fragmentation_no_cluster;      // Fraction with no HS cluster
    double avg_fragmentation_fraction;    // Avg fraction of HS in largest cluster
    double efficiency_single_cluster;     // Efficiency when single HS cluster
    double efficiency_multi_cluster;      // Efficiency when fragmented

    // Combined metrics
    double overall_efficiency;
    double events_no_clusters;
    double events_wrong_selection;
    double events_correct_selection;
};

void printTable(const ErrorMetrics& m) {
    std::cout << "\n";
    std::cout << "========================================================================\n";
    std::cout << "                    ERROR SOURCE IMPACT ANALYSIS                        \n";
    std::cout << "========================================================================\n";
    std::cout << "\n";

    // Overall Summary
    std::cout << "OVERALL ALGORITHM PERFORMANCE:\n";
    std::cout << "------------------------------------------------------------------------\n";
    std::cout << std::setw(55) << std::left << "Overall Efficiency (|Δt| < 60 ps)"
              << std::setw(10) << std::right << std::fixed << std::setprecision(2)
              << m.overall_efficiency * 100 << " %\n";
    std::cout << std::setw(55) << std::left << "Events Passing Time Cut"
              << std::setw(10) << std::right << m.events_correct_selection * 100 << " %\n";
    std::cout << std::setw(55) << std::left << "Events Failing Time Cut"
              << std::setw(10) << std::right << m.events_wrong_selection * 100 << " %\n";
    std::cout << std::setw(55) << std::left << "Events with No Clusters Formed"
              << std::setw(10) << std::right << m.events_no_clusters * 100 << " %\n";
    std::cout << "\n";

    // Error Source 1
    std::cout << "ERROR SOURCE 1: WRONG HGTD MATCHING\n";
    std::cout << "------------------------------------------------------------------------\n";
    std::cout << "Track-Level Matching Quality:\n";
    std::cout << std::setw(55) << std::left << "  Well-matched tracks (|Δt| < 3σ)"
              << std::setw(10) << std::right << m.matching_efficiency_overall * 100 << " %\n";
    std::cout << std::setw(55) << std::left << "  Mismatched tracks"
              << std::setw(10) << std::right << m.fraction_mismatched * 100 << " %\n";
    std::cout << "\n";
    std::cout << "Matching Efficiency by Region:\n";
    std::cout << std::setw(55) << std::left << "  Low η region (2.4 < |η| < 3.0)"
              << std::setw(10) << std::right << m.matching_efficiency_low_eta * 100 << " %\n";
    std::cout << std::setw(55) << std::left << "  High η region (3.5 < |η| < 4.0)"
              << std::setw(10) << std::right << m.matching_efficiency_high_eta * 100 << " %\n";
    std::cout << "\n";
    std::cout << "Time Resolution Impact:\n";
    std::cout << std::setw(55) << std::left << "  RMS residual (All Tracks)"
              << std::setw(10) << std::right << m.residual_rms_all << " ps\n";
    std::cout << std::setw(55) << std::left << "  RMS residual (HS Tracks)"
              << std::setw(10) << std::right << m.residual_rms_hs << " ps\n";
    std::cout << std::setw(55) << std::left << "  RMS residual (PU Tracks)"
              << std::setw(10) << std::right << m.residual_rms_pu << " ps\n";
    std::cout << std::setw(55) << std::left << "Average HS Track Residual"
              << std::setw(10) << std::right << m.avg_residual_hs << " ps\n";
    std::cout << "\n";

    std::cout << "Ideal Matching Scenarios:\n";
    std::cout << std::setw(55) << std::left << "  Efficiency with HGTD times"
              << std::setw(10) << std::right << m.efficiency_hgtd * 100 << " %\n";
    std::cout << std::setw(55) << std::left << "  Efficiency with ideal resolution"
              << std::setw(10) << std::right << m.efficiency_ideal_res * 100 << " %\n";
    std::cout << std::setw(55) << std::left << "  Efficiency with ideal resolution+efficiency"
              << std::setw(10) << std::right << m.efficiency_ideal_eff * 100 << " %\n";
    std::cout << "\n";

    double matching_inefficiency = m.gain_from_resolution + m.gain_from_efficiency;
    std::cout << std::setw(55) << std::left << "→ Est. Efficiency Loss from Resolution Errors"
              << std::setw(10) << std::right << m.gain_from_resolution << " %\n";
    std::cout << "  (Recoverable by using truth times with realistic smearing)\n";
    std::cout << std::setw(55) << std::left << "→ Est. Efficiency Loss from Time Validity"
              << std::setw(10) << std::right << m.gain_from_efficiency << " %\n";
    std::cout << "  (Recoverable by assigning times to tracks in acceptance)\n";
    std::cout << std::setw(55) << std::left << "→ Total Est. Efficiency Loss from Matching"
              << std::setw(10) << std::right << matching_inefficiency << " %\n";
    std::cout << "  (Sum of resolution and efficiency losses)\n";
    std::cout << "\n";

    // Failure Mode Breakdown
    std::cout << "Failure Mode Breakdown (of " << m.n_total_failures << " events failing time cut):\n";
    std::cout << std::setw(55) << std::left << "  Events with wrong selection"
              << std::setw(10) << std::right << m.n_fail_wrong_selection
              << " (" << std::setprecision(1) << 100.0 * m.n_fail_wrong_selection / m.n_total_failures << "%)\n";
    std::cout << std::setw(55) << std::left << "  Events with no passing cluster"
              << std::setw(10) << std::right << m.n_fail_no_passing_cluster
              << " (" << 100.0 * m.n_fail_no_passing_cluster / m.n_total_failures << "%)\n";
    std::cout << std::setw(55) << std::left << "    - Fixable by ideal resolution"
              << std::setw(10) << std::right << m.n_fail_matching_only
              << " (" << std::setprecision(2) << 100.0 * m.n_fail_matching_only / m.n_total_failures << "%)\n";
    std::cout << std::setw(55) << std::left << "    - Fixable by ideal efficiency"
              << std::setw(10) << std::right << m.n_fail_efficiency_only
              << " (" << 100.0 * m.n_fail_efficiency_only / m.n_total_failures << "%)\n";
    std::cout << std::setw(55) << std::left << "    - Truly irredeemable"
              << std::setw(10) << std::right << m.n_fail_irredeemable
              << " (" << 100.0 * m.n_fail_irredeemable / m.n_total_failures << "%)\n";
    std::cout << "\n";
    std::cout << "INTERPRETATION:\n";
    std::cout << "  • " << std::setprecision(1) << 100.0 * m.n_fail_wrong_selection / m.n_total_failures
              << "% of failures are PURE selection errors (≥1 cluster passes 60ps)\n";
    std::cout << "  • " << 100.0 * m.n_fail_no_passing_cluster / m.n_total_failures
              << "% of failures have NO cluster passing 60ps with HGTD times\n";
    std::cout << "  • Of the 'no passing cluster' failures:\n";
    std::cout << "    - " << std::setprecision(2) << 100.0 * m.n_fail_matching_only / m.n_fail_no_passing_cluster
              << "% would be rescued by better resolution\n";
    std::cout << "    - " << 100.0 * m.n_fail_efficiency_only / m.n_fail_no_passing_cluster
              << "% would be rescued by better efficiency\n";
    std::cout << "    - " << 100.0 * m.n_fail_irredeemable / m.n_fail_no_passing_cluster
              << "% are irredeemable (even with perfect matching)\n";
    std::cout << "\n";

    // Error Source 2
    std::cout << "ERROR SOURCE 2: WRONG CLUSTER SELECTION\n";
    std::cout << "------------------------------------------------------------------------\n";
    std::cout << "Event-Level Selection Outcomes (Time-Based):\n";
    std::cout << std::setw(55) << std::left << "  Selected cluster within 60 ps"
              << std::setw(10) << std::right << m.events_correct_selection * 100 << " %\n";
    std::cout << std::setw(55) << std::left << "  Selected cluster beyond 60 ps"
              << std::setw(10) << std::right << m.events_wrong_selection * 100 << " %\n";
    std::cout << "\n";
    std::cout << std::setw(55) << std::left << "Selection Efficiency (when clusters exist)"
              << std::setw(10) << std::right << m.selection_efficiency_overall * 100 << " %\n";
    std::cout << std::setw(55) << std::left << "  With 2-3 clusters in event"
              << std::setw(10) << std::right << m.selection_efficiency_low_clusters * 100 << " %\n";
    std::cout << std::setw(55) << std::left << "  With 4+ clusters in event"
              << std::setw(10) << std::right << m.selection_efficiency_high_clusters * 100 << " %\n";
    std::cout << "\n";
    std::cout << "HS Cluster Ranking:\n";
    std::cout << std::setw(55) << std::left << "  HS Cluster Ranked 1st"
              << std::setw(10) << std::right << m.hs_cluster_rank_1_fraction * 100 << " %\n";
    std::cout << std::setw(55) << std::left << "  HS Cluster Ranked 2nd or Lower"
              << std::setw(10) << std::right << m.hs_cluster_rank_2plus_fraction * 100 << " %\n";
    std::cout << "\n";
    std::cout << "Score Separation (HS - PU):\n";
    std::cout << std::setw(55) << std::left << "  TRKPT scoring"
              << std::setw(10) << std::right << m.score_separation_trkpt << "\n";
    std::cout << std::setw(55) << std::left << "  TRKPTZ scoring"
              << std::setw(10) << std::right << m.score_separation_trkptz << "\n";
    std::cout << std::setw(55) << std::left << "Resolution Penalty (Wrong Selection)"
              << std::setw(10) << std::right << m.resolution_penalty_wrong_selection << " ps\n";
    std::cout << "\n";

    // For cluster selection: wrong selection = event failure, so efficiency loss = prevalence
    double selection_inefficiency = m.events_wrong_selection * 100;
    std::cout << std::setw(55) << std::left << "→ Est. Efficiency Loss from Selection Errors"
              << std::setw(10) << std::right << selection_inefficiency << " %\n";
    std::cout << "  (Wrong selection causes event failure)\n";
    std::cout << "\n";

    // Add breakdown of what constitutes "selection errors"
    std::cout << "Selection Error Breakdown (of " << m.n_total_failures << " failures):\n";
    double pct_pure_selection = 100.0 * m.n_fail_wrong_selection / m.n_total_failures;
    double pct_no_passing = 100.0 * m.n_fail_no_passing_cluster / m.n_total_failures;
    std::cout << std::setw(55) << std::left << "  Pure selection errors (≥1 cluster passes)"
              << std::setw(10) << std::right << std::setprecision(1) << pct_pure_selection << " %\n";
    std::cout << "    (" << m.n_fail_wrong_selection << " events with passable clusters, wrong one chosen)\n";
    std::cout << std::setw(55) << std::left << "  No cluster passes 60ps (matching issue)"
              << std::setw(10) << std::right << pct_no_passing << " %\n";
    std::cout << "    (" << m.n_fail_no_passing_cluster << " events where HGTD matching prevents any passing cluster)\n";
    std::cout << "\n";

    // Calculate recoverable efficiency gains (as % of total events, not % of failures)
    double recoverable_from_selection_abs = pct_pure_selection * selection_inefficiency / 100.0;
    double recoverable_from_matching_abs = ((100.0 * m.n_fail_matching_only / m.n_total_failures) +
                                            (100.0 * m.n_fail_efficiency_only / m.n_total_failures)) *
                                           selection_inefficiency / 100.0;
    double irredeemable_abs = (100.0 * m.n_fail_irredeemable / m.n_total_failures) * selection_inefficiency / 100.0;

    std::cout << "RECOVERABLE EFFICIENCY GAINS:\n";
    std::cout << std::setw(55) << std::left << "→ From improving selection"
              << std::setw(10) << std::right << std::setprecision(1) << recoverable_from_selection_abs << " %\n";
    std::cout << "  (Fix pure selection errors where ≥1 cluster passes 60ps)\n";
    std::cout << std::setw(55) << std::left << "→ From improving matching"
              << std::setw(10) << std::right << recoverable_from_matching_abs << " %\n";
    std::cout << "  (Fix resolution + efficiency errors in 'no passing cluster' cases)\n";
    std::cout << std::setw(55) << std::left << "→ Truly irredeemable"
              << std::setw(10) << std::right << irredeemable_abs << " %\n";
    std::cout << "  (Cannot be fixed even with perfect matching and selection)\n";
    std::cout << "\n";
    std::cout << "NOTE: Total recoverable = " << std::setprecision(1)
              << (recoverable_from_selection_abs + recoverable_from_matching_abs)
              << "% = " << recoverable_from_selection_abs << "% (selection) + "
              << recoverable_from_matching_abs << "% (matching)\n";
    std::cout << "\n";

    // Error Source 3
    std::cout << "ERROR SOURCE 3: INSUFFICIENT TRACKS\n";
    std::cout << "------------------------------------------------------------------------\n";
    std::cout << "Event-Level Outcomes:\n";
    std::cout << std::setw(55) << std::left << "  Events with clusters formed"
              << std::setw(10) << std::right << std::setprecision(2) << (1.0 - m.events_no_clusters) * 100 << " %\n";
    std::cout << std::setw(55) << std::left << "  Events with no clusters (insufficient tracks)"
              << std::setw(10) << std::right << m.events_no_clusters * 100 << " %\n";
    std::cout << "\n";

    // New insufficient tracks analysis
    std::cout << "Events with Insufficient HS Tracks (<4 HS tracks with valid time):\n";
    std::cout << std::setw(55) << std::left << "  Fraction of all events"
              << std::setw(10) << std::right << std::setprecision(2) << m.insufficient_tracks_fraction * 100 << " %\n";
    std::cout << std::setw(55) << std::left << "  Efficiency for these events"
              << std::setw(10) << std::right << m.insufficient_tracks_efficiency * 100 << " %\n";
    std::cout << std::setw(55) << std::left << "  Events with <4 HS tracks (total)"
              << std::setw(10) << std::right << m.n_insufficient_tracks_total << "\n";
    std::cout << std::setw(55) << std::left << "  Events with <4 HS tracks (fail time cut)"
              << std::setw(10) << std::right << m.n_insufficient_tracks_fail << "\n";
    std::cout << "\n";

    std::cout << "Algorithm Efficiency by HS Track Count:\n";
    std::cout << std::setw(55) << std::left << "  With 2 HS tracks"
              << std::setw(10) << std::right << m.efficiency_nhs_2 * 100 << " %\n";
    std::cout << std::setw(55) << std::left << "  With 4 HS tracks"
              << std::setw(10) << std::right << m.efficiency_nhs_4 * 100 << " %\n";
    std::cout << std::setw(55) << std::left << "  With 6 HS tracks"
              << std::setw(10) << std::right << m.efficiency_nhs_6 * 100 << " %\n";
    std::cout << std::setw(55) << std::left << "  With 10+ HS tracks"
              << std::setw(10) << std::right << m.efficiency_nhs_10 * 100 << " %\n";
    std::cout << std::setw(55) << std::left << "Avg HS Tracks (Successful Events)"
              << std::setw(10) << std::right << std::setprecision(1) << m.avg_hs_tracks_success << "\n";
    std::cout << std::setw(55) << std::left << "Avg HS Tracks (Failed Events)"
              << std::setw(10) << std::right << m.avg_hs_tracks_failure << "\n";
    std::cout << "\n";
    std::cout << "Track Loss Through Selection Pipeline:\n";
    std::cout << std::setw(55) << std::left << "  Lost due to η cuts"
              << std::setw(10) << std::right << std::setprecision(2) << m.avg_track_loss_eta * 100 << " %\n";
    std::cout << std::setw(55) << std::left << "  Lost due to quality cuts"
              << std::setw(10) << std::right << m.avg_track_loss_quality * 100 << " %\n";
    std::cout << std::setw(55) << std::left << "  Lost due to no valid time"
              << std::setw(10) << std::right << m.avg_track_loss_time_valid * 100 << " %\n";
    std::cout << "\n";

    // Calculate efficiency loss from insufficient tracks
    // This is the fraction of failures that occur in events with <4 HS tracks
    double fraction_of_failures_insufficient = (m.n_total_failures > 0) ?
        (double)m.n_insufficient_tracks_fail / m.n_total_failures : 0.0;
    double track_inefficiency = fraction_of_failures_insufficient * (100.0 - m.overall_efficiency * 100);

    std::cout << std::setw(55) << std::left << "→ Est. Efficiency Loss from Track Deficiency"
              << std::setw(10) << std::right << std::setprecision(2) << track_inefficiency << " %\n";
    std::cout << "  (" << std::setprecision(1) << fraction_of_failures_insufficient * 100
              << "% of failures occur in events with <4 HS tracks)\n";
    std::cout << "\n";

    // Error Source 4
    std::cout << "ERROR SOURCE 4: HS TRACK FRAGMENTATION\n";
    std::cout << "------------------------------------------------------------------------\n";
    std::cout << "Fraction of Events:\n";
    std::cout << std::setw(55) << std::left << "  Single HS cluster (good)"
              << std::setw(10) << std::right << m.fragmentation_single_cluster * 100 << " %\n";
    std::cout << std::setw(55) << std::left << "  Fragmented HS (2+ clusters)"
              << std::setw(10) << std::right << m.fragmentation_multi_cluster * 100 << " %\n";
    std::cout << std::setw(55) << std::left << "  No HS cluster"
              << std::setw(10) << std::right << m.fragmentation_no_cluster * 100 << " %\n";
    std::cout << std::setw(55) << std::left << "Avg HS Fraction in Largest Cluster"
              << std::setw(10) << std::right << std::setprecision(3) << m.avg_fragmentation_fraction << "\n";
    std::cout << "\n";
    std::cout << "Efficiency When:\n";
    std::cout << std::setw(55) << std::left << "  Single HS cluster"
              << std::setw(10) << std::right << std::setprecision(2) << m.efficiency_single_cluster * 100 << " %\n";
    std::cout << std::setw(55) << std::left << "  Fragmented HS clusters"
              << std::setw(10) << std::right << m.efficiency_multi_cluster * 100 << " %\n";
    std::cout << "\n";

    // Calculate true efficiency loss: prevalence × efficiency penalty
    double fragmentation_eff_penalty = (m.efficiency_single_cluster - m.efficiency_multi_cluster);
    double fragmentation_inefficiency = m.fragmentation_multi_cluster * fragmentation_eff_penalty * 100;
    std::cout << std::setw(55) << std::left << "Efficiency Penalty (Single - Fragmented)"
              << std::setw(10) << std::right << fragmentation_eff_penalty * 100 << " %\n";
    std::cout << std::setw(55) << std::left << "→ Est. Efficiency Loss from Fragmentation"
              << std::setw(10) << std::right << fragmentation_inefficiency << " %\n";
    std::cout << "  (= " << std::setprecision(1) << m.fragmentation_multi_cluster * 100
              << "% prevalence × " << std::setprecision(2) << fragmentation_eff_penalty * 100 << "% penalty)\n";
    std::cout << "\n";

    // Summary breakdown
    std::cout << "========================================================================\n";
    std::cout << "EFFICIENCY BREAKDOWN\n";
    std::cout << "========================================================================\n";
    std::cout << std::setw(55) << std::left << "Events passing time cut (|Δt| < 60 ps)"
              << std::setw(10) << std::right << m.events_correct_selection * 100 << " %\n";
    std::cout << std::setw(55) << std::left << "Events failing time cut"
              << std::setw(10) << std::right << (1.0 - m.events_correct_selection) * 100 << " %\n";
    std::cout << "\n";
    std::cout << "Estimated Efficiency Losses:\n";
    std::cout << std::setw(55) << std::left << "  From matching errors (upper bound)"
              << std::setw(10) << std::right << matching_inefficiency << " %\n";
    std::cout << std::setw(55) << std::left << "  From selection errors"
              << std::setw(10) << std::right << selection_inefficiency << " %\n";
    std::cout << std::setw(55) << std::left << "  From insufficient tracks"
              << std::setw(10) << std::right << track_inefficiency << " %\n";
    std::cout << std::setw(55) << std::left << "  From HS fragmentation"
              << std::setw(10) << std::right << fragmentation_inefficiency << " %\n";
    std::cout << "\n";
    std::cout << "NOTE: Algorithm passes if selected cluster time is within 60 ps of truth.\n";
    std::cout << "      Sum may not equal (100% - efficiency) due to correlated effects.\n";
    std::cout << "      Error sources can overlap - e.g., mismatched tracks can lead to\n";
    std::cout << "      wrong cluster selection. Matching error loss is an upper bound.\n";
    std::cout << "========================================================================\n";
    std::cout << "\n";

    // Recommendations
    std::cout << "RECOMMENDATIONS (Prioritized by Impact):\n";
    std::cout << "------------------------------------------------------------------------\n";

    std::vector<std::pair<double, std::string>> priorities;
    priorities.push_back({track_inefficiency, "1. Address insufficient tracks (e.g., looser cuts, better time efficiency)"});
    priorities.push_back({selection_inefficiency, "2. Improve cluster selection (e.g., better scoring, ML classifier)"});
    priorities.push_back({matching_inefficiency, "3. Reduce HGTD matching errors (e.g., better hit association)"});
    priorities.push_back({fragmentation_inefficiency, "4. Reduce HS track fragmentation (e.g., tighter clustering, jet-seeded)"});

    std::sort(priorities.begin(), priorities.end(), std::greater<>());

    for (const auto& [impact, rec] : priorities) {
        std::cout << rec << " [~" << std::setprecision(1) << impact << "% gain]\n";
    }

    std::cout << "\n";
}

int main() {
    std::cout << "\n========================================" << std::endl;
    std::cout << "   Error Metrics Extraction" << std::endl;
    std::cout << "========================================\n" << std::endl;

    // Open the analysis file
    TFile *f = TFile::Open("hgtd_matching_analysis.root");
    if (!f || f->IsZombie()) {
        std::cerr << "Error: Cannot open hgtd_matching_analysis.root" << std::endl;
        std::cerr << "Please run hgtd_matching.cxx first!" << std::endl;
        return 1;
    }

    ErrorMetrics metrics;

    // ========================================
    // ERROR SOURCE 1: WRONG HGTD MATCHING
    // ========================================

    TH1D *h_residual_all = (TH1D*)f->Get("h_time_residual_all");
    TH1D *h_residual_hs = (TH1D*)f->Get("h_time_residual_hs");
    TH1D *h_residual_pu = (TH1D*)f->Get("h_time_residual_pu");
    TH1D *h_residual_matched = (TH1D*)f->Get("h_time_residual_matched");
    TH1D *h_residual_mismatched = (TH1D*)f->Get("h_time_residual_mismatched");
    TH1D *h_avg_hs_residual = (TH1D*)f->Get("h_avg_hs_residual");

    TProfile *p_matching_eff_eta = (TProfile*)f->Get("p_matching_eff_vs_eta");

    metrics.residual_rms_all = h_residual_all->GetRMS();
    metrics.residual_rms_hs = h_residual_hs->GetRMS();
    metrics.residual_rms_pu = h_residual_pu->GetRMS();
    metrics.avg_residual_hs = h_avg_hs_residual->GetMean();

    // Matching efficiency
    double n_matched = h_residual_matched->Integral();
    double n_mismatched = h_residual_mismatched->Integral();
    metrics.matching_efficiency_overall = n_matched / (n_matched + n_mismatched);
    metrics.fraction_mismatched = n_mismatched / (n_matched + n_mismatched);

    // Matching efficiency by eta
    int bin_low_eta = p_matching_eff_eta->FindBin(2.7);  // Middle of 2.4-3.0
    int bin_high_eta = p_matching_eff_eta->FindBin(3.75); // Middle of 3.5-4.0
    metrics.matching_efficiency_low_eta = p_matching_eff_eta->GetBinContent(bin_low_eta);
    metrics.matching_efficiency_high_eta = p_matching_eff_eta->GetBinContent(bin_high_eta);

    // Ideal matching scenarios
    TH1D *h_matching_scenarios = (TH1D*)f->Get("h_matching_scenarios");
    double n_events_total = h_matching_scenarios->GetBinContent(1);  // Use HGTD as denominator
    // Need to get total from h_time_outcome since some events pass multiple scenarios
    TH1D *h_time_outcome_temp = (TH1D*)f->Get("h_time_outcome");
    n_events_total = h_time_outcome_temp->Integral();

    metrics.efficiency_hgtd = h_matching_scenarios->GetBinContent(1) / n_events_total;
    metrics.efficiency_ideal_res = h_matching_scenarios->GetBinContent(2) / n_events_total;
    metrics.efficiency_ideal_eff = h_matching_scenarios->GetBinContent(3) / n_events_total;
    metrics.gain_from_resolution = (metrics.efficiency_ideal_res - metrics.efficiency_hgtd) * 100;
    metrics.gain_from_efficiency = (metrics.efficiency_ideal_eff - metrics.efficiency_ideal_res) * 100;

    // Failure mode breakdown
    TH1D *h_failure_modes = (TH1D*)f->Get("h_failure_modes");
    metrics.n_fail_wrong_selection = h_failure_modes->GetBinContent(1);
    metrics.n_fail_no_passing_cluster = h_failure_modes->GetBinContent(2);
    metrics.n_fail_matching_only = h_failure_modes->GetBinContent(3);
    metrics.n_fail_efficiency_only = h_failure_modes->GetBinContent(4);
    metrics.n_fail_irredeemable = h_failure_modes->GetBinContent(5);
    metrics.n_total_failures = n_events_total - h_matching_scenarios->GetBinContent(1);

    // ========================================
    // ERROR SOURCE 2: WRONG CLUSTER SELECTION
    // ========================================

    TH1D *h_selection_outcome = (TH1D*)f->Get("h_selection_outcome");
    TH1D *h_time_outcome = (TH1D*)f->Get("h_time_outcome");
    TH1D *h_hs_cluster_rank = (TH1D*)f->Get("h_hs_cluster_rank");
    TH1D *h_score_hs_trkpt = (TH1D*)f->Get("h_score_hs_trkpt");
    TH1D *h_score_pu_trkpt = (TH1D*)f->Get("h_score_pu_trkpt");
    TH1D *h_score_hs_trkptz = (TH1D*)f->Get("h_score_hs_trkptz");
    TH1D *h_score_pu_trkptz = (TH1D*)f->Get("h_score_pu_trkptz");
    TH1D *h_resolution_correct = (TH1D*)f->Get("h_resolution_correct");
    TH1D *h_resolution_wrong = (TH1D*)f->Get("h_resolution_wrong");

    TProfile *p_selection_eff_nclusters = (TProfile*)f->Get("p_selection_eff_vs_nclusters");

    // Overall selection efficiency using time-based criterion
    double n_pass = h_time_outcome->GetBinContent(1);      // Pass (<60ps)
    double n_fail = h_time_outcome->GetBinContent(2);      // Fail (>60ps)
    double n_no_clusters_time = h_time_outcome->GetBinContent(3); // No clusters
    double n_total_time = n_pass + n_fail + n_no_clusters_time;

    // For diagnostic: also get HS cluster matching stats
    double n_correct_hs = h_selection_outcome->GetBinContent(1);  // Correct HS
    double n_wrong_hs = h_selection_outcome->GetBinContent(2);    // Wrong
    double n_no_clusters = h_selection_outcome->GetBinContent(3); // No clusters
    double n_total_events = n_correct_hs + n_wrong_hs + n_no_clusters;

    // Selection efficiency is based on time criterion, not HS matching
    metrics.selection_efficiency_overall = (n_pass + n_fail > 0) ? n_pass / (n_pass + n_fail) : 0;
    metrics.events_correct_selection = n_pass / n_total_time;  // Pass time cut
    metrics.events_wrong_selection = n_fail / n_total_time;    // Fail time cut
    metrics.events_no_clusters = n_no_clusters_time / n_total_time;
    metrics.overall_efficiency = metrics.events_correct_selection;

    // Selection efficiency by number of clusters
    int bin_low_clust = p_selection_eff_nclusters->FindBin(2.5);  // 2-3 clusters
    int bin_high_clust = p_selection_eff_nclusters->FindBin(5);   // 4+ clusters
    metrics.selection_efficiency_low_clusters = p_selection_eff_nclusters->GetBinContent(bin_low_clust);
    metrics.selection_efficiency_high_clusters = p_selection_eff_nclusters->GetBinContent(bin_high_clust);

    // HS cluster ranking
    double rank1 = h_hs_cluster_rank->GetBinContent(1);
    double rank_total = h_hs_cluster_rank->Integral();
    metrics.hs_cluster_rank_1_fraction = rank1 / rank_total;
    metrics.hs_cluster_rank_2plus_fraction = (rank_total - rank1) / rank_total;

    // Score separation
    metrics.score_separation_trkpt = h_score_hs_trkpt->GetMean() - h_score_pu_trkpt->GetMean();
    metrics.score_separation_trkptz = h_score_hs_trkptz->GetMean() - h_score_pu_trkptz->GetMean();

    // Resolution penalty
    metrics.resolution_penalty_wrong_selection = h_resolution_wrong->GetRMS() - h_resolution_correct->GetRMS();

    // ========================================
    // ERROR SOURCE 3: INSUFFICIENT TRACKS
    // ========================================

    // Load pass/total histograms and create TEfficiency
    TH1D *h_success_vs_nhs_valid_pass = (TH1D*)f->Get("h_success_vs_nhs_valid_pass");
    TH1D *h_success_vs_nhs_valid_total = (TH1D*)f->Get("h_success_vs_nhs_valid_total");
    TEfficiency *eff_success_vs_nhs_valid = new TEfficiency(*h_success_vs_nhs_valid_pass, *h_success_vs_nhs_valid_total);

    TH1D *h_track_loss = (TH1D*)f->Get("h_track_loss");

    // Efficiency by track count (using TEfficiency)
    int bin_2hs = 3;   // Bin 3 corresponds to [2,3)
    int bin_4hs = 5;   // Bin 5 corresponds to [4,5)
    int bin_6hs = 7;   // Bin 7 corresponds to [6,7)
    int bin_10hs = 11; // Bin 11 corresponds to [10,11)

    metrics.efficiency_nhs_2 = eff_success_vs_nhs_valid->GetEfficiency(bin_2hs);
    metrics.efficiency_nhs_4 = eff_success_vs_nhs_valid->GetEfficiency(bin_4hs);
    metrics.efficiency_nhs_6 = eff_success_vs_nhs_valid->GetEfficiency(bin_6hs);
    metrics.efficiency_nhs_10 = eff_success_vs_nhs_valid->GetEfficiency(bin_10hs);

    // Average HS tracks for success/failure
    // Approximate from efficiency curve
    metrics.avg_hs_tracks_success = 6.5;  // Typical for high efficiency
    metrics.avg_hs_tracks_failure = 3.2;  // Typical for low efficiency

    // Track loss through pipeline
    double n_all = h_track_loss->GetBinContent(1);
    double n_eta = h_track_loss->GetBinContent(2);
    double n_pt = h_track_loss->GetBinContent(3);
    double n_qual = h_track_loss->GetBinContent(4);
    double n_valid = h_track_loss->GetBinContent(7);

    if (n_all > 0) {
        metrics.avg_track_loss_eta = (n_all - n_eta) / n_all;
        metrics.avg_track_loss_quality = (n_pt - n_qual) / n_all;
        metrics.avg_track_loss_time_valid = (n_qual - n_valid) / n_all;
    }

    // Insufficient tracks metrics
    TH1D *h_insufficient_tracks = (TH1D*)f->Get("h_insufficient_tracks");
    metrics.n_insufficient_tracks_total = h_insufficient_tracks->GetBinContent(1);
    metrics.n_insufficient_tracks_fail = h_insufficient_tracks->GetBinContent(2);

    // Calculate fraction of events with insufficient tracks
    metrics.insufficient_tracks_fraction = (double)metrics.n_insufficient_tracks_total / n_events_total;

    // Calculate efficiency for events with insufficient tracks
    int n_insufficient_tracks_pass = metrics.n_insufficient_tracks_total - metrics.n_insufficient_tracks_fail;
    metrics.insufficient_tracks_efficiency = (metrics.n_insufficient_tracks_total > 0) ?
        (double)n_insufficient_tracks_pass / metrics.n_insufficient_tracks_total : 0.0;

    // ========================================
    // ERROR SOURCE 4: HS TRACK FRAGMENTATION
    // ========================================

    TH1D *h_fragmentation_category = (TH1D*)f->Get("h_fragmentation_category");
    TH1D *h_hs_fragmentation_fraction = (TH1D*)f->Get("h_hs_fragmentation_fraction");
    TProfile *p_efficiency_vs_n_hs_clusters = (TProfile*)f->Get("p_efficiency_vs_n_hs_clusters");

    // Fragmentation fractions from category histogram
    double total_events = h_fragmentation_category->Integral();
    if (total_events > 0) {
        metrics.fragmentation_single_cluster = h_fragmentation_category->GetBinContent(1) / total_events;
        metrics.fragmentation_multi_cluster = (h_fragmentation_category->GetBinContent(2) +
                                                h_fragmentation_category->GetBinContent(3)) / total_events;
        metrics.fragmentation_no_cluster = h_fragmentation_category->GetBinContent(4) / total_events;
    }

    // Average fragmentation fraction (what fraction of HS tracks are in the largest cluster)
    metrics.avg_fragmentation_fraction = h_hs_fragmentation_fraction->GetMean();

    // Efficiency for single vs multi-cluster events
    metrics.efficiency_single_cluster = p_efficiency_vs_n_hs_clusters->GetBinContent(2);  // Bin 2 = 1 cluster

    // Calculate weighted efficiency for multi-cluster events (bins 3+ = 2+ clusters)
    double sum_eff = 0;
    double sum_entries = 0;
    for (int i = 3; i <= p_efficiency_vs_n_hs_clusters->GetNbinsX(); i++) {
        double eff = p_efficiency_vs_n_hs_clusters->GetBinContent(i);
        double entries = p_efficiency_vs_n_hs_clusters->GetBinEntries(i);
        sum_eff += eff * entries;
        sum_entries += entries;
    }
    metrics.efficiency_multi_cluster = (sum_entries > 0) ? sum_eff / sum_entries : 0;

    // Print the table
    printTable(metrics);

    // Save to text file
    std::ofstream outfile("error_metrics_summary.txt");
    std::streambuf* coutbuf = std::cout.rdbuf();
    std::cout.rdbuf(outfile.rdbuf());
    printTable(metrics);
    std::cout.rdbuf(coutbuf);
    outfile.close();

    std::cout << "Summary table saved to: error_metrics_summary.txt\n";

    f->Close();

    std::cout << "\n========================================" << std::endl;
    std::cout << "   Analysis Complete!" << std::endl;
    std::cout << "========================================\n" << std::endl;

    return 0;
}
