/*
 * hgtd_matching.cxx
 *
 * Comprehensive HGTD matching error analysis for vertex timing reconstruction
 *
 * Quantifies three major error sources:
 * 1. Wrong HGTD matching (incorrect track time assignments)
 * 2. Wrong cluster selection (incorrect hard scatter identification)
 * 3. Insufficient tracks (too few tracks for reliable clustering)
 *
 * Usage: root -l -b -q hgtd_matching.cxx
 */

#include "clustering_includes.h"
#include "clustering_constants.h"
#include "clustering_structs.h"
#include "clustering_functions.h"
#include "event_processing.h"

using namespace MyUtl;

int main() {
  std::cout << "\n========================================" << std::endl;
  std::cout << "   HGTD Matching Error Analysis" << std::endl;
  std::cout << "========================================\n" << std::endl;

  // Setup chain
  TChain chain("ntuple");
  setupChain(chain, "../../ntuple-hgtd/");
  TTreeReader reader(&chain);
  BranchPointerWrapper branch(reader);

  // Output file
  TFile *fout = new TFile("hgtd_matching_analysis.root", "RECREATE");

  // ========================================
  // ERROR SOURCE 1: WRONG HGTD MATCHING
  // ========================================

  // 1.1 Track-Level Time Residuals
  TH1D *h_time_residual_all = new TH1D("h_time_residual_all",
				       "All Track Time Residuals;#Deltat = t_{track} - t_{truth} [ps];Fraction of tracks",
				       200, -500, 500);

  TH1D *h_time_residual_hs = new TH1D("h_time_residual_hs",
				      "HS Track Time Residuals;#Deltat [ps];Fraction of tracks",
				      200, -500, 500);

  TH1D *h_time_residual_pu = new TH1D("h_time_residual_pu",
				      "PU Track Time Residuals;#Deltat [ps];Fraction of tracks",
				      200, -500, 500);

  TH1D *h_time_residual_matched = new TH1D("h_time_residual_matched",
					   "Well-Matched Tracks (|#Deltat| < 3#sigma);#Deltat [ps];Fraction of tracks",
					   200, -500, 500);

  TH1D *h_time_residual_mismatched = new TH1D("h_time_residual_mismatched",
					      "Mis-Matched Tracks (|#Deltat| > 3#sigma);#Deltat [ps];Fraction of tracks",
					      200, -500, 500);

  // 1.2 Matching Quality vs Kinematics
  TProfile *p_matching_eff_vs_eta = new TProfile("p_matching_eff_vs_eta",
						 "Matching Efficiency vs |#eta|;|#eta|;Efficiency",
						 30, 2.4, 4.0);

  TProfile *p_matching_eff_vs_pt = new TProfile("p_matching_eff_vs_pt",
						"Matching Efficiency vs p_{T};p_{T} [GeV];Efficiency",
						30, 1, 30);

  TProfile *p_matching_eff_vs_nhits = new TProfile("p_matching_eff_vs_nhits",
						   "Matching Efficiency vs nHGTDHits;nHGTDHits;Efficiency",
						   10, 0, 10);

  TH2D *h2_residual_vs_eta = new TH2D("h2_residual_vs_eta",
				      "Time Residual vs #eta;#eta;#Deltat [ps]",
				      50, -5, 5, 100, -500, 500);

  TH2D *h2_residual_vs_pt = new TH2D("h2_residual_vs_pt",
				     "Time Residual vs p_{T};p_{T} [GeV];#Deltat [ps]",
				     50, 0, 100, 100, -500, 500);

  TH2D *h2_residual_vs_nhits = new TH2D("h2_residual_vs_nhits",
					"Time Residual vs nHGTDHits;nHGTDHits;#Deltat [ps]",
					10, 0, 10, 100, -500, 500);

  TProfile *p_residual_rms_vs_nhits = new TProfile("p_residual_rms_vs_nhits",
						   "RMS(#Deltat) vs nHGTDHits;nHGTDHits;RMS(#Deltat) [ps]",
						   10, 0, 10, "s");  // "s" option for RMS

  // 1.3 Event-Level Matching Quality
  TH1D *h_avg_hs_residual = new TH1D("h_avg_hs_residual",
				     "Average HS Track Residual per Event;Avg |#Deltat| [ps];Events",
				     100, 0, 200);

  TProfile *p_resolution_vs_avg_residual = new TProfile("p_resolution_vs_avg_residual",
							"Vertex Resolution vs Avg HS Residual;Avg |#Deltat|_{HS} [ps];Vertex Resolution [ps]",
							50, 0, 200);

  // 1.4 Residual Cumulative Distribution
  TH1D *h_residual_abs_all = new TH1D("h_residual_abs_all",
				      "Cumulative |#Deltat| Distribution (All);|#Deltat| [ps];Cumulative Fraction",
				      200, 0, 500);

  TH1D *h_residual_abs_hs = new TH1D("h_residual_abs_hs",
				     "Cumulative |#Deltat| Distribution (HS);|#Deltat| [ps];Cumulative Fraction",
				     200, 0, 500);

  TH1D *h_residual_abs_pu = new TH1D("h_residual_abs_pu",
				     "Cumulative |#Deltat| Distribution (PU);|#Deltat| [ps];Cumulative Fraction",
				     200, 0, 500);

  // 1.5 Ideal Resolution Scenarios
  TH1D *h_matching_scenarios = new TH1D("h_matching_scenarios",
					"Event Outcomes by Matching Scenario;Scenario;Events Passing 60ps Cut",
					3, 0, 3);
  h_matching_scenarios->GetXaxis()->SetBinLabel(1, "HGTD Times");
  h_matching_scenarios->GetXaxis()->SetBinLabel(2, "Ideal Resolution");
  h_matching_scenarios->GetXaxis()->SetBinLabel(3, "Ideal Res+Eff");

  // 1.6 Failure Mode Breakdown
  TH1D *h_failure_modes = new TH1D("h_failure_modes",
				   "Failure Mode Breakdown;Mode;Events",
				   5, 0, 5);
  h_failure_modes->GetXaxis()->SetBinLabel(1, "Wrong Selection");
  h_failure_modes->GetXaxis()->SetBinLabel(2, "No Passing Cluster");
  h_failure_modes->GetXaxis()->SetBinLabel(3, "Fixable by Resolution");
  h_failure_modes->GetXaxis()->SetBinLabel(4, "Fixable by Efficiency");
  h_failure_modes->GetXaxis()->SetBinLabel(5, "Irredeemable");

  // Insufficient tracks histogram
  TH1D *h_insufficient_tracks = new TH1D("h_insufficient_tracks",
					 "Insufficient Tracks Metrics;Metric;Events",
					 2, 0, 2);
  h_insufficient_tracks->GetXaxis()->SetBinLabel(1, "Total with <N HS");
  h_insufficient_tracks->GetXaxis()->SetBinLabel(2, "Fail with <N HS");

  // Track events that fail HGTD but pass ideal scenarios
  int n_events_rescued_by_ideal_res = 0;
  int n_events_rescued_by_ideal_eff = 0;
  int n_events_total_processed = 0;

  // Track failure modes for events that fail time cut
  int n_fail_no_passing_cluster = 0;      // No cluster within 60ps (irredeemable by selection)
  int n_fail_wrong_selection = 0;         // At least one cluster within 60ps, but wrong one selected
  int n_fail_matching_only = 0;           // Would pass with ideal resolution
  int n_fail_efficiency_only = 0;         // Would pass with ideal efficiency but not resolution

  // Track insufficient tracks for events that fail
  int n_insufficient_tracks_total = 0;    // Total events with too few HS tracks with valid time
  int n_insufficient_tracks_fail = 0;     // Events with too few HS tracks that fail time cut
  const int INSUFFICIENT_TRACK_THRESHOLD = 4;  // Define "too few" as < 4 HS tracks with valid time

  // ========================================
  // ERROR SOURCE 2: WRONG CLUSTER SELECTION
  // ========================================

  // 2.1 Cluster-Level Metrics
  TH1D *h_n_clusters = new TH1D("h_n_clusters",
				"Number of Clusters per Event;N_{clusters};Events",
				20, 0, 20);

  TH1D *h_hs_cluster_rank = new TH1D("h_hs_cluster_rank",
				     "HS Cluster Rank by Score;Rank;Fraction of Events",
				     10, 0.5, 10.5);

  TH1D *h_score_ratio = new TH1D("h_score_ratio",
				 "Score Gap (TRKPTZ);(Score_{1st} - Score_{2nd}) / Score_{1st} [%];Fraction of Events",
				 50, 0, 100);

  // Score distributions for different methods
  std::map<Score, TH1D*> h_score_hs;
  std::map<Score, TH1D*> h_score_pu;

  h_score_hs[TRKPT]  = new TH1D("h_score_hs_trkpt",  "HS Cluster Score (TRKPT);Score;Clusters",  50, 0, 200);
  h_score_pu[TRKPT]  = new TH1D("h_score_pu_trkpt",  "PU Cluster Score (TRKPT);Score;Clusters",  50, 0, 200);
  h_score_hs[TRKPTZ] = new TH1D("h_score_hs_trkptz", "HS Cluster Score (TRKPTZ);Score;Clusters", 50, 0, 200);
  h_score_pu[TRKPTZ] = new TH1D("h_score_pu_trkptz", "PU Cluster Score (TRKPTZ);Score;Clusters", 50, 0, 200);
  h_score_hs[TESTML] = new TH1D("h_score_hs_dnn",    "HS Cluster Score (DNN);DNN Score;Clusters", 50, 0, 1);
  h_score_pu[TESTML] = new TH1D("h_score_pu_dnn",    "PU Cluster Score (DNN);DNN Score;Clusters", 50, 0, 1);

  // DNN-based HS cluster ranking and score ratio (parallel to existing TRKPTZ quantities)
  TH1D *h_hs_cluster_rank_dnn = new TH1D("h_hs_cluster_rank_dnn",
    "HS Cluster Rank by DNN Score;Rank;Fraction of Events",
    10, 0.5, 10.5);

  TH1D *h_score_ratio_dnn = new TH1D("h_score_ratio_dnn",
    "Score Gap (DNN);(Score_{1st} - Score_{2nd}) / Score_{1st} [%];Fraction of Events",
    50, 0, 100);

  // DNN time-based outcome (parallel to h_time_outcome which uses TRKPTZ)
  TH1D *h_time_outcome_dnn = new TH1D("h_time_outcome_dnn",
    "Time-Based Outcome (DNN);Outcome;Events",
    3, 0, 3);
  h_time_outcome_dnn->GetXaxis()->SetBinLabel(1, "Pass (<60ps)");
  h_time_outcome_dnn->GetXaxis()->SetBinLabel(2, "Fail (>60ps)");
  h_time_outcome_dnn->GetXaxis()->SetBinLabel(3, "No clusters");

  // DNN selection summary for extract_error_metrics.cxx
  // bin 1: DNN selects HS cluster correctly (pure selection check)
  // bin 2: DNN selects wrong cluster (HS exists but not chosen)
  // bin 3: events with wrong selection where HS ranked 1st by DNN
  TH1D *h_dnn_selection_outcome = new TH1D("h_dnn_selection_outcome",
    "DNN Cluster Selection Outcome;Outcome;Events",
    4, 0, 4);
  h_dnn_selection_outcome->GetXaxis()->SetBinLabel(1, "Correct HS");
  h_dnn_selection_outcome->GetXaxis()->SetBinLabel(2, "Wrong (HS exists)");
  h_dnn_selection_outcome->GetXaxis()->SetBinLabel(3, "No clusters");
  h_dnn_selection_outcome->GetXaxis()->SetBinLabel(4, "HS absent");

  // 2.2 Selection Confusion Matrix
  int n_tp = 0, n_fp = 0, n_fn = 0, n_tn = 0;  // True/False Positive/Negative

  TH1D *h_selection_outcome = new TH1D("h_selection_outcome",
				       "Cluster Selection Outcome;Outcome;Events",
				       4, 0, 4);
  h_selection_outcome->GetXaxis()->SetBinLabel(1, "Correct HS");
  h_selection_outcome->GetXaxis()->SetBinLabel(2, "Wrong (HS exists)");
  h_selection_outcome->GetXaxis()->SetBinLabel(3, "No clusters");
  h_selection_outcome->GetXaxis()->SetBinLabel(4, "HS absent");

  // Track time-based outcomes (60 ps criterion)
  TH1D *h_time_outcome = new TH1D("h_time_outcome",
				  "Time-Based Outcome;Outcome;Events",
				  3, 0, 3);
  h_time_outcome->GetXaxis()->SetBinLabel(1, "Pass (<60ps)");
  h_time_outcome->GetXaxis()->SetBinLabel(2, "Fail (>60ps)");
  h_time_outcome->GetXaxis()->SetBinLabel(3, "No clusters");

  // 2.3 Selection efficiency vs event properties
  TProfile *p_selection_eff_vs_nclusters = new TProfile("p_selection_eff_vs_nclusters",
							"Selection Efficiency vs N_{clusters};N_{clusters};Efficiency",
							15, 0, 15);

  TProfile *p_selection_eff_vs_njets = new TProfile("p_selection_eff_vs_njets",
						    "Selection Efficiency vs N_{forward jets};N_{jets};Efficiency",
						    10, 0, 10);

  TProfile *p_selection_eff_vs_nhs = new TProfile("p_selection_eff_vs_nhs",
						  "Selection Efficiency vs N_{HS tracks};N_{HS};Efficiency",
						  20, 0, 20);

  TProfile *p_selection_eff_vs_pufrac = new TProfile("p_selection_eff_vs_pufrac",
						     "Selection Efficiency vs PU Fraction;PU Fraction;Efficiency",
						     20, 0, 1);

  // 2.4 Resolution impact
  TH1D *h_resolution_correct = new TH1D("h_resolution_correct",
					"Vertex Resolution (Correct Selection);#Deltat_{vertex} [ps];Events",
					100, -200, 200);

  TH1D *h_resolution_wrong = new TH1D("h_resolution_wrong",
				      "Vertex Resolution (Wrong Selection);#Deltat_{vertex} [ps];Events",
				      100, -200, 200);

  TH1D *h_extra_error_from_selection = new TH1D("h_extra_error_from_selection",
						"Additional Error from Wrong Selection;Extra |#Deltat| [ps];Events",
						100, 0, 500);

  // ========================================
  // ERROR SOURCE 3: INSUFFICIENT TRACKS
  // ========================================

  // 3.1 Track Multiplicity
  TH1D *h_n_forward_tracks = new TH1D("h_n_forward_tracks",
				      "Forward Track Multiplicity;N_{forward tracks};Events",
				      50, 0, 50);

  TH1D *h_n_forward_hs = new TH1D("h_n_forward_hs",
				  "Forward HS Tracks;N_{HS tracks};Events",
				  30, 0, 30);

  TH1D *h_n_forward_hs_valid = new TH1D("h_n_forward_hs_valid",
					"Forward HS Tracks with Valid Time;N_{HS valid};Events",
					30, 0, 30);

  TH1D *h_n_forward_pu = new TH1D("h_n_forward_pu",
				  "Forward PU Tracks;N_{PU tracks};Events",
				  50, 0, 50);

  TProfile *p_time_validity_vs_eta = new TProfile("p_time_validity_vs_eta",
						  "Time Validity Fraction vs |#eta|;|#eta|;Valid Time Fraction",
						  30, 2.4, 4.0);

  // 3.2 Success rate vs track count (using TEfficiency like main algorithm)
  TH1D *h_success_vs_nhs_total = new TH1D("h_success_vs_nhs_total", "", 30, 0, 30);
  TH1D *h_success_vs_nhs_pass = new TH1D("h_success_vs_nhs_pass", "", 30, 0, 30);

  TH1D *h_success_vs_nhs_valid_total = new TH1D("h_success_vs_nhs_valid_total", "", 30, 0, 30);
  TH1D *h_success_vs_nhs_valid_pass = new TH1D("h_success_vs_nhs_valid_pass", "", 30, 0, 30);

  TH2D *h2_success_vs_hs_pu = new TH2D("h2_success_vs_hs_pu",
				       "Success Rate vs Track Counts;N_{HS valid};N_{PU valid}",
				       20, 0, 20, 30, 0, 30);

  // 3.3 Cluster size distributions
  TH1D *h_cluster_size_all = new TH1D("h_cluster_size_all",
				      "Cluster Size (All);N_{tracks};Clusters",
				      30, 0, 30);

  TH1D *h_cluster_size_hs = new TH1D("h_cluster_size_hs",
				     "Cluster Size (HS);N_{tracks};Clusters",
				     30, 0, 30);

  TH1D *h_cluster_size_pu = new TH1D("h_cluster_size_pu",
				     "Cluster Size (PU);N_{tracks};Clusters",
				     30, 0, 30);

  TH1D *h_n_hs_in_hs_cluster = new TH1D("h_n_hs_in_hs_cluster",
					"HS Tracks in HS Cluster;N_{HS in cluster};Clusters",
					30, 0, 30);

  // 3.4 Track loss waterfall
  TH1D *h_track_loss = new TH1D("h_track_loss",
				"Track Loss through Selection;Selection Stage;Avg N_{tracks}",
				8, 0, 8);
  h_track_loss->GetXaxis()->SetBinLabel(1, "All tracks");
  h_track_loss->GetXaxis()->SetBinLabel(2, "After |#eta| cut");
  h_track_loss->GetXaxis()->SetBinLabel(3, "After p_{T} cut");
  h_track_loss->GetXaxis()->SetBinLabel(4, "After quality");
  h_track_loss->GetXaxis()->SetBinLabel(5, "After vtx assoc");
  h_track_loss->GetXaxis()->SetBinLabel(6, "After PU removal");
  h_track_loss->GetXaxis()->SetBinLabel(7, "After time valid");
  h_track_loss->GetXaxis()->SetBinLabel(8, "HS tracks");

  // 3.5 HS track fragmentation
  TH1D *h_hs_time_spread = new TH1D("h_hs_time_spread",
				    "HS Track Time Spread (RMS);RMS_{t} [ps];Events",
				    50, 0, 200);

  TH1D *h_hs_z_spread = new TH1D("h_hs_z_spread",
				 "HS Track z_{0} Spread (RMS);RMS_{z} [mm];Events",
				 50, 0, 50);

  TProfile *p_success_vs_time_spread = new TProfile("p_success_vs_time_spread",
						    "Success vs Time Spread;RMS_{t} [ps];Success Rate",
						    40, 0, 200);

  TProfile *p_success_vs_z_spread = new TProfile("p_success_vs_z_spread",
						 "Success vs z Spread;RMS_{z} [mm];Success Rate",
						 40, 0, 50);

  // 3.6 Resolution vs track count
  TProfile *p_resolution_vs_nhs = new TProfile("p_resolution_vs_nhs",
					       "Resolution vs N_{HS in cluster};N_{HS};Resolution [ps]",
					       20, 0, 20);

  // ========================================
  // MISCLUSTERING STUDY: TESTML vs TEST_MISCL
  // ========================================
  // Quantifies the impact of misclustering on DNN efficiency by comparing:
  //   (1) TESTML: DNN score > 0.5 AND |Δt| < 60 ps  (overall DNN efficiency)
  //   (2) TEST_MISCL: same as (1) AND cluster purity > 0.75  (DNN efficiency with well-formed clusters)
  // The gap between the two curves = fraction of events misclustered but still passing time cut.

  // Summary counts (filled into h_misclustering for extract_error_metrics.cxx)
  //   bin 1: events where DNN passes time cut (TESTML)
  //   bin 2: events where DNN passes time cut AND purity > 0.75 (TEST_MISCL)
  //   bin 3: events where DNN selected cluster has purity > 0.75 (regardless of time)
  //   bin 4: total events entering this analysis (denominator)
  TH1D *h_misclustering = new TH1D("h_misclustering",
    "Misclustering Study;Metric;Events",
    4, 0, 4);
  h_misclustering->GetXaxis()->SetBinLabel(1, "DNN passes time");
  h_misclustering->GetXaxis()->SetBinLabel(2, "DNN+purity passes time");
  h_misclustering->GetXaxis()->SetBinLabel(3, "DNN selects pure cluster");
  h_misclustering->GetXaxis()->SetBinLabel(4, "Total events");

  // Purity distribution of the DNN-selected cluster
  TH1D *h_dnn_selected_purity = new TH1D("h_dnn_selected_purity",
    "Purity of DNN-Selected Cluster;Purity (HS p_{T} fraction);Events",
    50, 0, 1.05);

  // DNN score of the selected cluster (for both HS and PU selected clusters)
  TH1D *h_dnn_score_pure   = new TH1D("h_dnn_score_pure",
    "DNN Score (Pure Clusters, purity>0.75);DNN Score;Clusters",
    50, 0, 1);
  TH1D *h_dnn_score_impure = new TH1D("h_dnn_score_impure",
    "DNN Score (Impure Clusters, purity<0.75);DNN Score;Clusters",
    50, 0, 1);

  // Time resolution for DNN-selected clusters, split by purity
  TH1D *h_dnn_reso_pure   = new TH1D("h_dnn_reso_pure",
    "#Deltat_{vtx} (DNN, purity>0.75);#Deltat [ps];Events",
    100, -200, 200);
  TH1D *h_dnn_reso_impure = new TH1D("h_dnn_reso_impure",
    "#Deltat_{vtx} (DNN, purity<0.75);#Deltat [ps];Events",
    100, -200, 200);

  // Efficiency vs forward jet count for both DNN scenarios
  TH1D *h_dnn_eff_fjet_pass   = new TH1D("h_dnn_eff_fjet_pass",   "", 31, 0, 31);
  TH1D *h_dnn_eff_fjet_total  = new TH1D("h_dnn_eff_fjet_total",  "", 31, 0, 31);
  TH1D *h_miscl_eff_fjet_pass = new TH1D("h_miscl_eff_fjet_pass", "", 31, 0, 31);

  // Efficiency vs PU fraction for both DNN scenarios
  TH1D *h_dnn_eff_pufrac_pass   = new TH1D("h_dnn_eff_pufrac_pass",   "", 20, 0, 1);
  TH1D *h_dnn_eff_pufrac_total  = new TH1D("h_dnn_eff_pufrac_total",  "", 20, 0, 1);
  TH1D *h_miscl_eff_pufrac_pass = new TH1D("h_miscl_eff_pufrac_pass", "", 20, 0, 1);

  // ========================================
  // ERROR SOURCE 4: HS TRACK FRAGMENTATION
  // ========================================

  // 4.1 Cluster fragmentation metrics
  TH1D *h_n_clusters_with_hs = new TH1D("h_n_clusters_with_hs",
					"Number of Clusters Containing HS Tracks;N_{clusters with HS};Events",
					10, 0, 10);

  TH1D *h_hs_fragmentation_fraction = new TH1D("h_hs_fragmentation_fraction",
					       "Fraction of HS Tracks in Largest Cluster;Fraction;Events",
					       50, 0, 1.05);

  TH1D *h_n_clusters_passing = new TH1D("h_n_clusters_passing",
					"Number of Clusters Passing Efficiency Cut;N_{passing};Events",
					10, 0, 10);

  // 4.2 Fragmentation impact on efficiency
  TProfile *p_efficiency_vs_n_hs_clusters = new TProfile("p_efficiency_vs_n_hs_clusters",
							 "Efficiency vs N_{clusters with HS};N_{HS clusters};Efficiency",
							 10, 0, 10);

  TProfile *p_efficiency_vs_fragmentation = new TProfile("p_efficiency_vs_fragmentation",
							 "Efficiency vs HS Fragmentation;HS Fraction in Largest;Efficiency",
							 20, 0, 1.05);

  // 4.3 Fragmentation categorization
  TH1D *h_fragmentation_category = new TH1D("h_fragmentation_category",
					    "Fragmentation Category;Category;Events",
					    4, 0, 4);
  h_fragmentation_category->GetXaxis()->SetBinLabel(1, "Single cluster");
  h_fragmentation_category->GetXaxis()->SetBinLabel(2, "2 clusters");
  h_fragmentation_category->GetXaxis()->SetBinLabel(3, "3+ clusters");
  h_fragmentation_category->GetXaxis()->SetBinLabel(4, "No HS cluster");

  // Counters for fragmentation statistics
  int n_events_fragmented = 0;
  int n_events_single_hs_cluster = 0;
  int n_events_no_hs_cluster = 0;

  // ========================================
  // EXAMPLE EVENT TRACKING
  // ========================================

  struct ExampleEvent {
    TString file_num;
    Long64_t event_num;
    double extra_time;
    double metric;  // For sorting (e.g., avg_hs_residual, n_mismatched, etc.)
  };

  std::vector<ExampleEvent> examples_matching_errors;     // Error Source 1: High avg HS residual
  std::vector<ExampleEvent> examples_ideal_rescue;        // Error Source 1: Rescued by ideal scenarios
  std::vector<ExampleEvent> examples_pure_selection;      // Error Source 2: Pure selection failure
  std::vector<ExampleEvent> examples_no_passing_res;      // Error Source 2: Rescued by resolution
  std::vector<ExampleEvent> examples_no_passing_eff;      // Error Source 2: Rescued by efficiency
  std::vector<ExampleEvent> examples_irredeemable;        // Error Source 2: Truly irredeemable
  std::vector<ExampleEvent> examples_low_multiplicity;    // Error Source 3: Low HS track count
  std::vector<ExampleEvent> examples_fragmented;          // Error Source 4: Fragmented HS clusters

  // ========================================
  // EVENT LOOP
  // ========================================

  std::cout << "Processing " << chain.GetEntries() << " events..." << std::endl;

  Long64_t nEvents = 0;
  Long64_t maxEvents = chain.GetEntries();  // Process subset for speed

  // Counters for track loss waterfall
  double sum_all = 0, sum_eta = 0, sum_pt = 0, sum_qual = 0;
  double sum_assoc = 0, sum_pu = 0, sum_valid = 0, sum_hs = 0;

  while (reader.Next() && nEvents < maxEvents) {
    if (nEvents % 1000 == 0) {
      std::cout << "Event " << nEvents << " / " << maxEvents << "\r" << std::flush;
    }
    nEvents++;

    // Apply same basic cuts as main algorithm
    if (!branch.passBasicCuts()) continue;
    if (!branch.passJetPtCut()) continue;

    // Extract file number and event number for event display
    TString fileName = branch.reader.GetTree()->GetCurrentFile()->GetName();
    TString file_num = fileName(49, 6);
    Long64_t event_num = chain.GetReadEntry() - chain.GetChainOffset();

    // ========================================
    // ERROR SOURCE 2 & 3: CLUSTERING ANALYSIS
    // ========================================

    // Track selection matching event_processing.h exactly:
    //   - 3.0σ cut used for counting (nForwardTrack* stats, consistent with processEventData)
    //   - MAX_NSIGMA (2.0) cut used for clustering, no separate pileupRemoval step
    std::vector<int> tracks_count = getAssociatedTracks(&branch, MIN_TRACK_PT, MAX_TRACK_PT, 3.0);
    std::vector<int> tracks       = getAssociatedTracks(&branch, MIN_TRACK_PT, MAX_TRACK_PT, MAX_NSIGMA);

    // ========================================
    // ERROR SOURCE 1: TRACK TIME RESIDUALS
    // ========================================
    // ONLY analyze tracks that pass pileup removal and are used in clustering

    double sum_hs_residual = 0;
    double sum_hs_residual_sq = 0;  // For RMS calculation
    int n_hs_with_time = 0;

    // Track-level analysis - ONLY iterate over clustering tracks
    for (auto trk : tracks) {
      if (branch.trackTimeValid[trk] != 1) continue;

      // Get truth time
      double truth_time = -999;
      int particle_idx = branch.trackToParticle[trk];
      int truth_vtx = branch.trackToTruthvtx[trk];

      if (particle_idx != -1) {
	truth_time = branch.particleT[particle_idx];
      } else if (truth_vtx != -1) {
	truth_time = branch.truthVtxTime[truth_vtx];
      } else {
	continue;  // Skip unmatched pileup
      }

      double reco_time = branch.trackTime[trk];
      double residual = reco_time - truth_time;
      double time_res = branch.trackTimeRes[trk];
      double abs_residual = std::abs(residual);

      // Apply forward cuts for relevant tracks
      double eta = std::abs(branch.trackEta[trk]);
      double pt = branch.trackPt[trk];
      bool is_forward = (eta > MIN_ABS_ETA_TRACK && eta < MAX_ABS_ETA_TRACK &&
			 pt > MIN_TRACK_PT && pt < MAX_TRACK_PT);

      if (!is_forward) continue;

      bool is_hs = (truth_vtx == 0);
      bool is_well_matched = (abs_residual < 3.0 * time_res);

      // Fill histograms
      h_time_residual_all->Fill(residual);
      h_residual_abs_all->Fill(abs_residual);

      if (is_hs) {
	h_time_residual_hs->Fill(residual);
	h_residual_abs_hs->Fill(abs_residual);
	sum_hs_residual += abs_residual;
	sum_hs_residual_sq += residual * residual;  // Use signed residual for RMS
	n_hs_with_time++;
      } else {
	h_time_residual_pu->Fill(residual);
	h_residual_abs_pu->Fill(abs_residual);
      }

      if (is_well_matched) {
	h_time_residual_matched->Fill(residual);
      } else {
	h_time_residual_mismatched->Fill(residual);
      }

      // Matching efficiency vs kinematics
      p_matching_eff_vs_eta->Fill(eta, is_well_matched ? 1.0 : 0.0);
      p_matching_eff_vs_pt->Fill(pt, is_well_matched ? 1.0 : 0.0);

      int nhits = branch.trackHgtdHits[trk];
      p_matching_eff_vs_nhits->Fill(nhits, is_well_matched ? 1.0 : 0.0);

      // 2D distributions
      h2_residual_vs_eta->Fill(branch.trackEta[trk], residual);
      h2_residual_vs_pt->Fill(pt, residual);
      h2_residual_vs_nhits->Fill(nhits, residual);
      p_residual_rms_vs_nhits->Fill(nhits, residual);

      // Time validity
      p_time_validity_vs_eta->Fill(eta, 1.0);
    }

    // Event-level HS residual
    if (n_hs_with_time > 0) {
      double avg_hs_residual = sum_hs_residual / n_hs_with_time;
      h_avg_hs_residual->Fill(avg_hs_residual);
    }

    // Count forward tracks using the 3.0σ track set, matching processEventData's counting step
    int n_forward_total = 0, n_forward_hs_total = 0, n_forward_pu_total = 0;
    branch.countForwardTracks(n_forward_total, n_forward_hs_total, n_forward_pu_total, tracks_count, true);
    int n_forward_hs_valid_total = n_forward_hs_total;  // With checkValidTimes=true, these are the same

    // Track loss waterfall
    int stage_all = branch.trackPt.GetSize();
    int stage_eta = 0, stage_pt = 0, stage_qual = 0;
    int stage_assoc = 0, stage_pu = 0, stage_valid = 0, stage_hs = 0;

    for (int trk = 0; trk < branch.trackPt.GetSize(); trk++) {
      double eta = std::abs(branch.trackEta[trk]);
      double pt = branch.trackPt[trk];
      bool quality = branch.trackQuality[trk];
      bool time_valid = (branch.trackTimeValid[trk] == 1);
      bool is_forward = (eta > MIN_ABS_ETA_TRACK && eta < MAX_ABS_ETA_TRACK);
      bool pt_ok = (pt > MIN_TRACK_PT && pt < MAX_TRACK_PT);

      if (is_forward) stage_eta++;
      if (is_forward && pt_ok) stage_pt++;
      if (is_forward && pt_ok && quality) stage_qual++;

      bool passes_assoc = false;
      if (is_forward && pt_ok && quality) {
	passes_assoc = passTrackVertexAssociation(trk, 0, &branch, MAX_NSIGMA);
	if (passes_assoc) stage_assoc++;
      }

      // Check if in final selection
      bool in_final = std::find(tracks.begin(), tracks.end(), trk) != tracks.end();
      if (in_final) {
	stage_pu++;
	if (time_valid) {
	  stage_valid++;
	  if (branch.trackToTruthvtx[trk] == 0) stage_hs++;
	}
      }
    }

    // Fill track loss waterfall
    sum_all += stage_all;
    sum_eta += stage_eta;
    sum_pt += stage_pt;
    sum_qual += stage_qual;
    sum_assoc += stage_assoc;
    sum_pu += stage_pu;
    sum_valid += stage_valid;
    sum_hs += stage_hs;

    // Fill track multiplicity histograms
    h_n_forward_tracks->Fill(n_forward_total);
    h_n_forward_hs->Fill(n_forward_hs_total);
    h_n_forward_hs_valid->Fill(n_forward_hs_valid_total);
    h_n_forward_pu->Fill(n_forward_pu_total);

    // Calculate HS track spread
    std::vector<double> hs_times, hs_z0s;
    for (auto trk : tracks) {
      if (branch.trackTimeValid[trk] == 1 && branch.trackToTruthvtx[trk] == 0) {
	hs_times.push_back(branch.trackTime[trk]);
	hs_z0s.push_back(branch.trackZ0[trk]);
      }
    }

    double time_spread = 0, z_spread = 0;
    if (hs_times.size() > 1) {
      double t_mean = std::accumulate(hs_times.begin(), hs_times.end(), 0.0) / hs_times.size();
      double z_mean = std::accumulate(hs_z0s.begin(), hs_z0s.end(), 0.0) / hs_z0s.size();

      for (auto t : hs_times) time_spread += (t - t_mean) * (t - t_mean);
      for (auto z : hs_z0s) z_spread += (z - z_mean) * (z - z_mean);

      time_spread = std::sqrt(time_spread / hs_times.size());
      z_spread = std::sqrt(z_spread / hs_z0s.size());

      h_hs_time_spread->Fill(time_spread);
      h_hs_z_spread->Fill(z_spread);
    }

    // Perform clustering (matching main algorithm parameters for "HGTD Times" mode)
    std::vector<Cluster> clusters = clusterTracksInTime(
							tracks, &branch,
							3.0,      // distanceCut
							false,    // useSmearedTimes
							true,     // checkTimeValid (MUST match main algorithm!)
							30.0,     // smearRes
							true,     // useCone (main algorithm uses true)
							false     // usez0
							);

    int n_clusters = clusters.size();
    h_n_clusters->Fill(n_clusters);

    if (n_clusters == 0) {
      h_selection_outcome->Fill(2);  // No clusters
      h_success_vs_nhs_total->Fill(n_forward_hs_total);
      h_success_vs_nhs_valid_total->Fill(n_forward_hs_valid_total);
      if (time_spread > 0) {
	p_success_vs_time_spread->Fill(time_spread, 0);
	p_success_vs_z_spread->Fill(z_spread, 0);
      }
      continue;
    }

    // Match clusters to truth vertices
    int hs_cluster_idx = -1;
    for (int ic = 0; ic < clusters.size(); ic++) {
      auto& cluster = clusters[ic];

      // Calculate cluster purity
      cluster.calcPurity(&branch);

      // Identify which truth vertex this cluster matches
      std::map<int, double> vtx_pt;
      for (auto trk : cluster.trackIndices) {
	int truth_vtx = branch.trackToTruthvtx[trk];
	if (truth_vtx != -1) {
	  vtx_pt[truth_vtx] += branch.trackPt[trk];
	}
      }

      int best_vtx = -1;
      double max_pt = 0;
      for (auto [vtx, pt] : vtx_pt) {
	if (pt > max_pt) {
	  max_pt = pt;
	  best_vtx = vtx;
	}
      }

      bool is_hs_cluster = (best_vtx == 0);
      if (is_hs_cluster) {
	hs_cluster_idx = ic;

	// Count HS tracks in HS cluster
	int n_hs_in_cluster = 0;
	for (auto trk : cluster.trackIndices) {
	  if (branch.trackToTruthvtx[trk] == 0) n_hs_in_cluster++;
	}
	h_n_hs_in_hs_cluster->Fill(n_hs_in_cluster);
      }

      // Fill cluster size histograms
      h_cluster_size_all->Fill(cluster.trackIndices.size());
      if (is_hs_cluster) {
	h_cluster_size_hs->Fill(cluster.trackIndices.size());
      } else {
	h_cluster_size_pu->Fill(cluster.trackIndices.size());
      }

      // Fill score distributions (TRKPT, TRKPTZ, and DNN)
      double score_trkpt  = cluster.scores[TRKPT];
      double score_trkptz = cluster.scores[TRKPTZ];
      double score_dnn    = cluster.scores[TESTML];

      if (is_hs_cluster) {
	h_score_hs[TRKPT]->Fill(score_trkpt);
	h_score_hs[TRKPTZ]->Fill(score_trkptz);
	h_score_hs[TESTML]->Fill(score_dnn);
      } else {
	h_score_pu[TRKPT]->Fill(score_trkpt);
	h_score_pu[TRKPTZ]->Fill(score_trkptz);
	h_score_pu[TESTML]->Fill(score_dnn);
      }
    }

    // Select cluster using TRKPTZ score
    Cluster selected_cluster = chooseCluster(clusters, TRKPTZ);

    // Find which cluster was selected
    int selected_idx = -1;
    for (int ic = 0; ic < clusters.size(); ic++) {
      if (clusters[ic] == selected_cluster) {
	selected_idx = ic;
	break;
      }
    }

    // Determine selection outcome based on time resolution
    // An event passes if |selected_time - truth_time| <= 3 * PASS_SIGMA (60 ps)
    const double PASS_THRESHOLD = 3.0 * PASS_SIGMA;  // 60 ps
    double vtx_time_truth = branch.truthVtxTime[0];
    double selected_time = (selected_cluster.values.size() > 0) ? selected_cluster.values[0] : -999.0;
    double time_residual = selected_time - vtx_time_truth;
    bool passes_time_cut = (selected_cluster.values.size() > 0) && (std::abs(time_residual) <= PASS_THRESHOLD);

    bool correct_selection = (selected_idx == hs_cluster_idx);
    bool hs_exists = (hs_cluster_idx != -1);

    // Fill time-based outcome histogram (TRKPTZ)
    if (n_clusters == 0) {
      h_time_outcome->Fill(2);  // No clusters
    } else if (passes_time_cut) {
      h_time_outcome->Fill(0);  // Pass
    } else {
      h_time_outcome->Fill(1);  // Fail
    }

    // -------------------------------------------------------
    // DNN (TESTML) parallel selection for Error Source 2
    // -------------------------------------------------------
    Cluster dnn_selected_cluster = chooseCluster(clusters, TESTML);
    int dnn_selected_idx = -1;
    for (int ic = 0; ic < (int)clusters.size(); ic++) {
      if (clusters[ic] == dnn_selected_cluster) { dnn_selected_idx = ic; break; }
    }
    double dnn_sel_time     = (dnn_selected_cluster.values.size() > 0) ? dnn_selected_cluster.values[0] : -999.0;
    double dnn_sel_residual = dnn_sel_time - vtx_time_truth;
    bool   dnn_passes_time  = (dnn_selected_cluster.values.size() > 0) &&
                              (std::abs(dnn_sel_residual) <= PASS_THRESHOLD);
    bool   dnn_correct      = (dnn_selected_idx == hs_cluster_idx);

    // Time-based outcome for DNN
    if (n_clusters == 0) {
      h_time_outcome_dnn->Fill(2);
    } else if (dnn_passes_time) {
      h_time_outcome_dnn->Fill(0);
    } else {
      h_time_outcome_dnn->Fill(1);
    }

    // DNN selection outcome (HS identity check)
    if (hs_exists) {
      h_dnn_selection_outcome->Fill(dnn_correct ? 0 : 1);
    } else {
      h_dnn_selection_outcome->Fill(3);  // HS absent
    }

    // HS cluster rank by DNN score + score ratio
    if (hs_exists) {
      std::vector<std::pair<double, int>> dnn_ranking;
      for (int ic = 0; ic < (int)clusters.size(); ic++)
        dnn_ranking.push_back({clusters[ic].scores[TESTML], ic});
      std::sort(dnn_ranking.begin(), dnn_ranking.end(), std::greater<>());

      int dnn_rank = 1;
      for (auto [sc, idx] : dnn_ranking) {
        if (idx == hs_cluster_idx) { h_hs_cluster_rank_dnn->Fill(dnn_rank); break; }
        dnn_rank++;
      }
      if (dnn_ranking.size() >= 2) {
        double s1 = dnn_ranking[0].first;
        double s2 = dnn_ranking[1].first;
        double dnn_gap_pct = 100.0 * (s1 - s2) / (s1 + 1e-9);
        h_score_ratio_dnn->Fill(dnn_gap_pct);
      }
    }

    // ========================================
    // IDEAL RESOLUTION SCENARIOS
    // ========================================
    // Test if this event would pass with ideal resolution and/or efficiency

    bool passes_ideal_res = passes_time_cut;  // Start with HGTD result
    bool passes_ideal_eff = passes_time_cut;

    if (!passes_time_cut && n_clusters > 0) {
      // Scenario 2: Ideal Resolution (smear truth times, keep valid time requirement)
      std::vector<Cluster> clusters_ideal_res = clusterTracksInTime(
								    tracks, &branch,
								    3.0,      // distanceCut
								    true,     // useSmearedTimes
								    true,     // checkTimeValid
								    30.0,     // smearRes
								    true,     // useCone
								    false     // usez0
								    );

      if (clusters_ideal_res.size() > 0) {
	Cluster selected_ideal_res = chooseCluster(clusters_ideal_res, TRKPTZ);
	if (selected_ideal_res.values.size() > 0) {
	  double time_res_ideal_res = selected_ideal_res.values[0] - vtx_time_truth;
	  passes_ideal_res = (std::abs(time_res_ideal_res) <= PASS_THRESHOLD);
	}
      }

      // Scenario 3: Ideal Resolution + Efficiency (smear truth times, ignore valid time requirement)
      std::vector<Cluster> clusters_ideal_eff = clusterTracksInTime(
								    tracks, &branch,
								    3.0,      // distanceCut
								    true,     // useSmearedTimes
								    false,    // checkTimeValid (ignore validity!)
								    30.0,     // smearRes
								    true,     // useCone
								    false     // usez0
								    );

      if (clusters_ideal_eff.size() > 0) {
	Cluster selected_ideal_eff = chooseCluster(clusters_ideal_eff, TRKPTZ);
	if (selected_ideal_eff.values.size() > 0) {
	  double time_res_ideal_eff = selected_ideal_eff.values[0] - vtx_time_truth;
	  passes_ideal_eff = (std::abs(time_res_ideal_eff) <= PASS_THRESHOLD);
	}
      }
    }

    // Track events rescued by ideal scenarios
    n_events_total_processed++;
    if (!passes_time_cut && passes_ideal_res) n_events_rescued_by_ideal_res++;
    if (!passes_time_cut && !passes_ideal_res && passes_ideal_eff) n_events_rescued_by_ideal_eff++;

    // Fill scenario histogram
    if (passes_time_cut) h_matching_scenarios->Fill(0);           // HGTD passes
    if (passes_ideal_res) h_matching_scenarios->Fill(1);          // Ideal res passes
    if (passes_ideal_eff) h_matching_scenarios->Fill(2);          // Ideal eff passes

    // Check if ANY cluster (using HGTD times) passes the 60ps threshold
    bool any_cluster_passes = false;

    // Categorize failure modes for events that fail time cut
    if (!passes_time_cut && n_clusters > 0) {
      for (int ic = 0; ic < clusters.size(); ic++) {
	if (clusters[ic].values.size() > 0) {
	  double cluster_time = clusters[ic].values[0];
	  double cluster_residual = std::abs(cluster_time - vtx_time_truth);
	  if (cluster_residual <= PASS_THRESHOLD) {
	    any_cluster_passes = true;
	    break;
	  }
	}
      }

      if (any_cluster_passes) {
	// At least one cluster passes, but we selected the wrong one
	n_fail_wrong_selection++;
      } else {
	// No cluster passes with HGTD times
	n_fail_no_passing_cluster++;

	// Check which ideal scenario would rescue it
	if (passes_ideal_res) {
	  n_fail_matching_only++;
	} else if (passes_ideal_eff) {
	  n_fail_efficiency_only++;
	}
	// If neither ideal scenario rescues it, it's truly irredeemable
      }
    }

    // ========================================
    // MISCLUSTERING STUDY: TESTML vs TEST_MISCL
    // ========================================
    // Re-run cluster selection with the DNN score (TESTML) and evaluate purity of the
    // selected cluster.  This is the direct per-event analogue of the TEST_MISCL score
    // added to the main analysis: same cluster set, same DNN selection, purity check here.

    const double PURITY_THRESHOLD = 0.75;  // must match TEST_MISCL in event_processing.h
    const double DNN_SCORE_CUT    = 0.5;   // must match TESTML cut in event_processing.h

    int nForwardJetMiscl = 0;
    branch.countForwardJets(nForwardJetMiscl);
    double pu_frac_miscl = (n_forward_total > 0) ?
      (double)n_forward_pu_total / n_forward_total : 0.0;

    h_misclustering->Fill(3);   // denominator: every event that reaches this point

    if (n_clusters > 0) {
      // updateScores sets TESTML (and TEST_MISCL) on every cluster; calcPurity is
      // already called in the cluster-identity loop above.
      // chooseCluster with TESTML picks the cluster with the highest DNN score.
      Cluster dnn_cluster = chooseCluster(clusters, TESTML);

      double dnn_score    = dnn_cluster.scores.at(TESTML);
      double dnn_purity   = dnn_cluster.purity;
      double dnn_time     = (dnn_cluster.values.size() > 0) ? dnn_cluster.values[0] : -999.0;
      double dnn_diff     = dnn_time - vtx_time_truth;
      bool dnn_passes_time = (dnn_cluster.values.size() > 0) &&
                             (std::abs(dnn_diff) <= PASS_THRESHOLD) &&
                             (dnn_score > DNN_SCORE_CUT);
      bool dnn_passes_purity = dnn_purity > PURITY_THRESHOLD;

      // Fill purity and score distributions
      h_dnn_selected_purity->Fill(dnn_purity);
      if (dnn_passes_purity) {
        h_dnn_score_pure->Fill(dnn_score);
        h_dnn_reso_pure->Fill(dnn_diff);
        h_misclustering->Fill(2);   // DNN selects a pure cluster
      } else {
        h_dnn_score_impure->Fill(dnn_score);
        h_dnn_reso_impure->Fill(dnn_diff);
      }

      // Event-level efficiency counters
      h_dnn_eff_fjet_total->Fill(nForwardJetMiscl);
      h_dnn_eff_pufrac_total->Fill(pu_frac_miscl);

      if (dnn_passes_time) {
        h_misclustering->Fill(0);   // TESTML: DNN passes time cut
        h_dnn_eff_fjet_pass->Fill(nForwardJetMiscl);
        h_dnn_eff_pufrac_pass->Fill(pu_frac_miscl);

        if (dnn_passes_purity) {
          h_misclustering->Fill(1); // TEST_MISCL: DNN+purity passes time cut
          h_miscl_eff_fjet_pass->Fill(nForwardJetMiscl);
          h_miscl_eff_pufrac_pass->Fill(pu_frac_miscl);
        }
      }
    }

    // Track insufficient tracks for all events
    if (n_forward_hs_valid_total < INSUFFICIENT_TRACK_THRESHOLD) {
      n_insufficient_tracks_total++;
      if (!passes_time_cut) {
	n_insufficient_tracks_fail++;
      }
    }

    // ========================================
    // TRACK EXAMPLE EVENTS
    // ========================================

    // Error Source 1: Events with high HS matching errors or rescued by ideal scenarios
    // n_hs_with_time and sum_hs_residual_sq already calculated from clustering tracks only
    if (n_hs_with_time >= 4) {
      double rms_hs_residual = std::sqrt(sum_hs_residual_sq / n_hs_with_time);
      // Track events with high HS RMS (> 100 ps indicates poor HGTD matching for HS)
      if (rms_hs_residual > 100.0 && !passes_time_cut) {
	examples_matching_errors.push_back({file_num, event_num, selected_time, rms_hs_residual});
      }
    }

    if (!passes_time_cut && (passes_ideal_res || passes_ideal_eff)) {
      // Event rescued by ideal scenario
      double improvement = passes_ideal_eff ? 2.0 : 1.0;  // Prioritize ideal_eff rescues
      examples_ideal_rescue.push_back({file_num, event_num, selected_time, improvement});
    }

    // Error Source 2: Categorize failure types
    if (!passes_time_cut && n_clusters > 0) {
      if (any_cluster_passes) {
	// Pure selection failure
	examples_pure_selection.push_back({file_num, event_num, selected_time, std::abs(time_residual)});
      } else {
	// No passing cluster
	if (passes_ideal_res && !passes_ideal_eff) {
	  examples_no_passing_res.push_back({file_num, event_num, selected_time, std::abs(time_residual)});
	} else if (passes_ideal_eff) {
	  examples_no_passing_eff.push_back({file_num, event_num, selected_time, std::abs(time_residual)});
	} else {
	  examples_irredeemable.push_back({file_num, event_num, selected_time, std::abs(time_residual)});
	}
      }
    }

    // Error Source 3: Low multiplicity events that fail
    if (n_forward_hs_valid_total < INSUFFICIENT_TRACK_THRESHOLD && !passes_time_cut) {
      examples_low_multiplicity.push_back({file_num, event_num, selected_time, (double)n_forward_hs_valid_total});
    }

    int outcome = -1;
    if (hs_exists && correct_selection) {
      outcome = 0;  // Correct
      n_tp++;
    } else if (hs_exists && !correct_selection) {
      outcome = 1;  // Wrong
      n_fn++;
    } else if (!hs_exists) {
      outcome = 3;  // HS absent
    }

    if (outcome >= 0) h_selection_outcome->Fill(outcome);

    // Fill efficiency profiles
    int nForwardJet = 0;
    branch.countForwardJets(nForwardJet);
    double pu_frac = (n_forward_total > 0) ? (double)n_forward_pu_total / n_forward_total : 0;

    // Define success based on time resolution (matching actual algorithm)
    double success = passes_time_cut ? 1.0 : 0.0;
    p_selection_eff_vs_nclusters->Fill(n_clusters, success);
    p_selection_eff_vs_njets->Fill(nForwardJet, success);
    p_selection_eff_vs_nhs->Fill(n_forward_hs_valid_total, success);
    p_selection_eff_vs_pufrac->Fill(pu_frac, success);

    // ========================================
    // ERROR SOURCE 4: HS FRAGMENTATION ANALYSIS
    // ========================================

    // Count how many clusters contain HS tracks
    std::vector<int> hs_track_counts;  // HS tracks per cluster
    int n_clusters_with_hs_tracks = 0;
    int total_hs_in_clusters = 0;
    int max_hs_in_single_cluster = 0;
    int n_clusters_passing_eff = 0;

    for (int ic = 0; ic < clusters.size(); ic++) {
      int n_hs_in_this_cluster = 0;
      for (auto trk : clusters[ic].trackIndices) {
	if (branch.trackToTruthvtx[trk] == 0) {
	  n_hs_in_this_cluster++;
	}
      }

      if (n_hs_in_this_cluster > 0) {
	n_clusters_with_hs_tracks++;
	total_hs_in_clusters += n_hs_in_this_cluster;
	hs_track_counts.push_back(n_hs_in_this_cluster);
	if (n_hs_in_this_cluster > max_hs_in_single_cluster) {
	  max_hs_in_single_cluster = n_hs_in_this_cluster;
	}
      }

      // Check if this cluster passes efficiency cut
      if (clusters[ic].values.size() > 0) {
	double cluster_time = clusters[ic].values[0];
	double cluster_residual = std::abs(cluster_time - vtx_time_truth);
	if (cluster_residual <= 3.0 * PASS_SIGMA) {
	  n_clusters_passing_eff++;
	}
      }
    }

    // Fill fragmentation histograms
    h_n_clusters_with_hs->Fill(n_clusters_with_hs_tracks);
    h_n_clusters_passing->Fill(n_clusters_passing_eff);

    // Calculate fragmentation fraction (what fraction of HS tracks in largest cluster)
    double fragmentation_fraction = 0;
    if (total_hs_in_clusters > 0) {
      fragmentation_fraction = (double)max_hs_in_single_cluster / total_hs_in_clusters;
      h_hs_fragmentation_fraction->Fill(fragmentation_fraction);
    }

    // Categorize fragmentation
    int fragmentation_cat = -1;
    if (n_clusters_with_hs_tracks == 0) {
      fragmentation_cat = 3;  // No HS cluster
      n_events_no_hs_cluster++;
    } else if (n_clusters_with_hs_tracks == 1) {
      fragmentation_cat = 0;  // Single cluster (good)
      n_events_single_hs_cluster++;
    } else if (n_clusters_with_hs_tracks == 2) {
      fragmentation_cat = 1;  // 2 clusters (fragmented)
      n_events_fragmented++;
    } else {
      fragmentation_cat = 2;  // 3+ clusters (highly fragmented)
      n_events_fragmented++;
    }
    h_fragmentation_category->Fill(fragmentation_cat);

    // Error Source 4: Track fragmented events
    // Require truly fragmented events: multiple clusters with multiple HS tracks each
    if (n_clusters_with_hs_tracks >= 2 && fragmentation_fraction < 0.7) {
      // Sort to find second-largest cluster
      std::vector<int> sorted_counts = hs_track_counts;
      std::sort(sorted_counts.begin(), sorted_counts.end(), std::greater<int>());

      // Require at least 2 HS tracks in each of the top 2 clusters
      if (sorted_counts.size() >= 2 && sorted_counts[0] >= 2 && sorted_counts[1] >= 2) {
        // Store fragmentation fraction as metric (lower = more fragmented)
        examples_fragmented.push_back({file_num, event_num, selected_time, fragmentation_fraction});
      }
    }

    // Fill efficiency vs fragmentation profiles
    p_efficiency_vs_n_hs_clusters->Fill(n_clusters_with_hs_tracks, success);
    if (total_hs_in_clusters > 0) {
      p_efficiency_vs_fragmentation->Fill(fragmentation_fraction, success);
    }

    // Fill TEfficiency histograms (total always, pass only if passes)
    h_success_vs_nhs_total->Fill(n_forward_hs_total);
    h_success_vs_nhs_valid_total->Fill(n_forward_hs_valid_total);
    if (passes_time_cut) {
      h_success_vs_nhs_pass->Fill(n_forward_hs_total);
      h_success_vs_nhs_valid_pass->Fill(n_forward_hs_valid_total);
    }

    if (time_spread > 0) {
      p_success_vs_time_spread->Fill(time_spread, success);
      p_success_vs_z_spread->Fill(z_spread, success);
    }

    // Fill 2D success map
    h2_success_vs_hs_pu->Fill(n_forward_hs_valid_total, n_forward_pu_total, success);

    // Rank of HS cluster
    if (hs_exists) {
      std::vector<std::pair<double, int>> score_ranking;
      for (int ic = 0; ic < clusters.size(); ic++) {
	score_ranking.push_back({clusters[ic].scores[TRKPTZ], ic});
      }
      std::sort(score_ranking.begin(), score_ranking.end(), std::greater<>());

      int rank = 1;
      for (auto [score, idx] : score_ranking) {
	if (idx == hs_cluster_idx) {
	  h_hs_cluster_rank->Fill(rank);
	  break;
	}
	rank++;
      }

      // Percentage gap between 1st and 2nd score: 100*(S1-S2)/S1
      // High value means the top cluster is clearly dominant; low value means ambiguous
      if (score_ranking.size() >= 2) {
        double s1 = score_ranking[0].first;
        double s2 = score_ranking[1].first;
        double gap_pct = 100.0 * (s1 - s2) / (s1 + 1e-9);
        h_score_ratio->Fill(gap_pct);
      }
    }

    // Resolution impact
    if (selected_cluster.values.size() > 0) {
      // Reuse time_residual already calculated above

      // Vertex resolution correlation with matching quality
      if (n_hs_with_time > 0) {
	double avg_hs_residual = sum_hs_residual / n_hs_with_time;
	p_resolution_vs_avg_residual->Fill(avg_hs_residual, std::abs(time_residual));

	// Also track RMS for better understanding of spread
	double rms_hs_residual = std::sqrt(sum_hs_residual_sq / n_hs_with_time);
      }

      if (correct_selection) {
	h_resolution_correct->Fill(time_residual);

	// Resolution vs cluster size
	int n_hs_in_selected = 0;
	for (auto trk : selected_cluster.trackIndices) {
	  if (branch.trackToTruthvtx[trk] == 0) n_hs_in_selected++;
	}
	p_resolution_vs_nhs->Fill(n_hs_in_selected, std::abs(time_residual));
      } else {
	h_resolution_wrong->Fill(time_residual);

	// Extra error from wrong selection
	if (hs_exists) {
	  double hs_time = clusters[hs_cluster_idx].values[0];
	  double correct_residual = hs_time - vtx_time_truth;
	  double extra_error = std::abs(time_residual) - std::abs(correct_residual);
	  if (extra_error > 0) {
	    h_extra_error_from_selection->Fill(extra_error);
	  }
	}
      }
    }
  }

  std::cout << "\nProcessed " << nEvents << " events" << std::endl;

  // Create TEfficiency objects from pass/total histograms
  TEfficiency *eff_success_vs_nhs = new TEfficiency(*h_success_vs_nhs_pass, *h_success_vs_nhs_total);
  eff_success_vs_nhs->SetName("eff_success_vs_nhs");
  eff_success_vs_nhs->SetTitle("Algorithm Success vs N_{HS};N_{HS tracks};Efficiency");

  TEfficiency *eff_success_vs_nhs_valid = new TEfficiency(*h_success_vs_nhs_valid_pass, *h_success_vs_nhs_valid_total);
  eff_success_vs_nhs_valid->SetName("eff_success_vs_nhs_valid");
  eff_success_vs_nhs_valid->SetTitle("Algorithm Success vs N_{HS valid};N_{HS valid};Efficiency");

  // Fill track loss waterfall
  h_track_loss->SetBinContent(1, sum_all / nEvents);
  h_track_loss->SetBinContent(2, sum_eta / nEvents);
  h_track_loss->SetBinContent(3, sum_pt / nEvents);
  h_track_loss->SetBinContent(4, sum_qual / nEvents);
  h_track_loss->SetBinContent(5, sum_assoc / nEvents);
  h_track_loss->SetBinContent(6, sum_pu / nEvents);
  h_track_loss->SetBinContent(7, sum_valid / nEvents);
  h_track_loss->SetBinContent(8, sum_hs / nEvents);

  // Convert to cumulative distributions
  for (int i = 1; i <= h_residual_abs_all->GetNbinsX(); i++) {
    double integral_all = h_residual_abs_all->Integral(1, i);
    double integral_hs = h_residual_abs_hs->Integral(1, i);
    double integral_pu = h_residual_abs_pu->Integral(1, i);

    double total_all = h_residual_abs_all->Integral();
    double total_hs = h_residual_abs_hs->Integral();
    double total_pu = h_residual_abs_pu->Integral();

    if (total_all > 0) h_residual_abs_all->SetBinContent(i, integral_all / total_all);
    if (total_hs > 0) h_residual_abs_hs->SetBinContent(i, integral_hs / total_hs);
    if (total_pu > 0) h_residual_abs_pu->SetBinContent(i, integral_pu / total_pu);
  }

  // Print summary statistics
  std::cout << "\n========================================" << std::endl;
  std::cout << "         SUMMARY STATISTICS" << std::endl;
  std::cout << "========================================\n" << std::endl;

  std::cout << "Cluster Selection Performance:" << std::endl;
  std::cout << "  True Positives:  " << n_tp << std::endl;
  std::cout << "  False Negatives: " << n_fn << std::endl;
  if (n_tp + n_fn > 0) {
    std::cout << "  Efficiency:      " << 100.0 * n_tp / (n_tp + n_fn) << "%" << std::endl;
  }

  std::cout << "\nHS Track Fragmentation:" << std::endl;
  std::cout << "  Single HS Cluster:   " << n_events_single_hs_cluster
	    << " (" << 100.0 * n_events_single_hs_cluster / nEvents << "%)" << std::endl;
  std::cout << "  Fragmented (2+ clusters): " << n_events_fragmented
	    << " (" << 100.0 * n_events_fragmented / nEvents << "%)" << std::endl;
  std::cout << "  No HS Cluster:       " << n_events_no_hs_cluster
	    << " (" << 100.0 * n_events_no_hs_cluster / nEvents << "%)" << std::endl;

  std::cout << "\nIdeal Matching Scenarios:" << std::endl;
  double efficiency_hgtd = h_matching_scenarios->GetBinContent(1) / n_events_total_processed;
  double efficiency_ideal_res = h_matching_scenarios->GetBinContent(2) / n_events_total_processed;
  double efficiency_ideal_eff = h_matching_scenarios->GetBinContent(3) / n_events_total_processed;

  std::cout << "  HGTD Times:              " << 100.0 * efficiency_hgtd << "%" << std::endl;
  std::cout << "  Ideal Resolution:        " << 100.0 * efficiency_ideal_res << "%" << std::endl;
  std::cout << "  Ideal Resolution+Eff:    " << 100.0 * efficiency_ideal_eff << "%" << std::endl;

  std::cout << "\n  Events rescued by ideal resolution:     " << n_events_rescued_by_ideal_res
	    << " (" << 100.0 * n_events_rescued_by_ideal_res / n_events_total_processed << "%)" << std::endl;
  std::cout << "  Events rescued by ideal efficiency:      " << n_events_rescued_by_ideal_eff
	    << " (" << 100.0 * n_events_rescued_by_ideal_eff / n_events_total_processed << "%)" << std::endl;

  double gain_from_resolution = 100.0 * (efficiency_ideal_res - efficiency_hgtd);
  double gain_from_efficiency = 100.0 * (efficiency_ideal_eff - efficiency_ideal_res);
  double total_matching_gain = 100.0 * (efficiency_ideal_eff - efficiency_hgtd);

  std::cout << "\n  Potential gain from fixing resolution:  " << gain_from_resolution << "%" << std::endl;
  std::cout << "  Potential gain from fixing efficiency:  " << gain_from_efficiency << "%" << std::endl;
  std::cout << "  Total potential gain from matching:     " << total_matching_gain << "%" << std::endl;

  std::cout << "\nFailure Mode Breakdown (for events failing time cut):" << std::endl;
  int n_total_failures = n_events_total_processed - h_matching_scenarios->GetBinContent(1);
  std::cout << "  Total events failing time cut:          " << n_total_failures << std::endl;
  std::cout << "  Events with wrong selection:            " << n_fail_wrong_selection
	    << " (" << 100.0 * n_fail_wrong_selection / n_total_failures << "%)" << std::endl;
  std::cout << "  Events with no passing cluster:         " << n_fail_no_passing_cluster
	    << " (" << 100.0 * n_fail_no_passing_cluster / n_total_failures << "%)" << std::endl;
  std::cout << "    - Fixable by ideal resolution:        " << n_fail_matching_only
	    << " (" << 100.0 * n_fail_matching_only / n_total_failures << "%)" << std::endl;
  std::cout << "    - Fixable by ideal efficiency:        " << n_fail_efficiency_only
	    << " (" << 100.0 * n_fail_efficiency_only / n_total_failures << "%)" << std::endl;
  std::cout << "    - Truly irredeemable:                 " << (n_fail_no_passing_cluster - n_fail_matching_only - n_fail_efficiency_only)
	    << " (" << 100.0 * (n_fail_no_passing_cluster - n_fail_matching_only - n_fail_efficiency_only) / n_total_failures << "%)" << std::endl;

  // ========================================
  // PRINT EXAMPLE EVENTS FOR EVENT DISPLAY
  // ========================================

  std::cout << "\n========================================" << std::endl;
  std::cout << "   EXAMPLE EVENTS FOR VISUALIZATION" << std::endl;
  std::cout << "========================================\n" << std::endl;

  auto print_examples = [](const char* title, std::vector<ExampleEvent>& events, int max_print = 30) {
    std::cout << title << ":\n";
    if (events.empty()) {
      std::cout << "  (No examples found)\n";
      return;
    }

    // Sort by metric (descending for most severe cases)
    std::sort(events.begin(), events.end(),
              [](const ExampleEvent& a, const ExampleEvent& b) { return a.metric > b.metric; });

    int n_print = std::min(max_print, (int)events.size());
    for (int i = 0; i < n_print; i++) {
      printf("python3 event_display.py --file_num %s --event_num %lld --extra_time %.2f\n",
             events[i].file_num.Data(), events[i].event_num, events[i].extra_time);
    }
    std::cout << "\n";
  };

  std::cout << "ERROR SOURCE 1: WRONG HGTD MATCHING\n";
  std::cout << "-----------------------------------\n";
  print_examples("High HS matching errors (RMS of HS track residuals > 100 ps, ≥4 HS tracks)", examples_matching_errors, 20);
  print_examples("Events rescued by ideal resolution/efficiency", examples_ideal_rescue, 20);

  std::cout << "ERROR SOURCE 2: WRONG CLUSTER SELECTION\n";
  std::cout << "---------------------------------------\n";
  print_examples("Pure selection failures (≥1 cluster passes 60ps)", examples_pure_selection, 20);
  print_examples("Failures rescued by ideal resolution", examples_no_passing_res, 20);
  print_examples("Failures rescued by ideal efficiency", examples_no_passing_eff, 20);
  print_examples("Truly irredeemable failures", examples_irredeemable, 20);

  std::cout << "ERROR SOURCE 3: INSUFFICIENT TRACKS\n";
  std::cout << "-----------------------------------\n";
  print_examples("Low multiplicity failures (<4 HS tracks with valid time)", examples_low_multiplicity, 20);

  std::cout << "ERROR SOURCE 4: HS TRACK FRAGMENTATION\n";
  std::cout << "--------------------------------------\n";
  print_examples("Events with fragmented HS clusters (2+ clusters)", examples_fragmented, 20);

  // ========================================
  // PLOTTING
  // ========================================

  gStyle->SetOptStat(0);
  gStyle->SetPalette(kViridis);

  // Create output directory for plots
  gSystem->mkdir("error_analysis_plots", true);

  // Plot 1.1: Track Time Residual Distributions
  TCanvas *c1 = new TCanvas("c1", "Track Time Residuals", 1400, 900);
  c1->Divide(2, 2);

  c1->cd(1);
  h_time_residual_all->SetLineColor(kBlack);
  h_time_residual_all->SetLineWidth(2);
  h_time_residual_all->Scale(1.0 / h_time_residual_all->Integral());
  h_time_residual_all->Draw("HIST");
  h_time_residual_all->Fit("gaus", "Q", "", -100, 100);

  c1->cd(2);
  h_time_residual_hs->SetLineColor(C01);
  h_time_residual_hs->SetLineWidth(2);
  h_time_residual_hs->Scale(1.0 / h_time_residual_hs->Integral());
  h_time_residual_hs->Draw("HIST");

  h_time_residual_pu->SetLineColor(C02);
  h_time_residual_pu->SetLineWidth(2);
  h_time_residual_pu->Scale(1.0 / h_time_residual_pu->Integral());
  h_time_residual_pu->Draw("HIST SAME");

  TLegend *leg1 = new TLegend(0.6, 0.7, 0.88, 0.88);
  leg1->AddEntry(h_time_residual_hs, "HS tracks", "l");
  leg1->AddEntry(h_time_residual_pu, "PU tracks", "l");
  leg1->Draw();

  c1->cd(3);
  h_time_residual_matched->SetLineColor(C08);
  h_time_residual_matched->SetLineWidth(2);
  h_time_residual_matched->Scale(1.0 / h_time_residual_matched->Integral());
  h_time_residual_matched->Draw("HIST");

  h_time_residual_mismatched->SetLineColor(C07);
  h_time_residual_mismatched->SetLineWidth(2);
  h_time_residual_mismatched->Scale(1.0 / h_time_residual_mismatched->Integral());
  h_time_residual_mismatched->Draw("HIST SAME");

  TLegend *leg2 = new TLegend(0.5, 0.7, 0.88, 0.88);
  leg2->AddEntry(h_time_residual_matched, "Well-matched", "l");
  leg2->AddEntry(h_time_residual_mismatched, "Mis-matched", "l");
  leg2->Draw();

  c1->cd(4);
  h_residual_abs_all->SetLineColor(kBlack);
  h_residual_abs_all->SetLineWidth(2);
  h_residual_abs_all->GetYaxis()->SetTitle("Cumulative Fraction");
  h_residual_abs_all->Draw("HIST");

  h_residual_abs_hs->SetLineColor(C01);
  h_residual_abs_hs->SetLineWidth(2);
  h_residual_abs_hs->Draw("HIST SAME");

  h_residual_abs_pu->SetLineColor(C02);
  h_residual_abs_pu->SetLineWidth(2);
  h_residual_abs_pu->Draw("HIST SAME");

  TLegend *leg3 = new TLegend(0.6, 0.2, 0.88, 0.4);
  leg3->AddEntry(h_residual_abs_all, "All", "l");
  leg3->AddEntry(h_residual_abs_hs, "HS", "l");
  leg3->AddEntry(h_residual_abs_pu, "PU", "l");
  leg3->Draw();

  c1->SaveAs("error_analysis_plots/plot1_time_residuals.pdf");

  // Plot 1.2: Matching Efficiency vs Kinematics
  TCanvas *c2 = new TCanvas("c2", "Matching Efficiency", 1800, 600);
  c2->Divide(3, 1);

  c2->cd(1);
  p_matching_eff_vs_eta->SetMarkerStyle(20);
  p_matching_eff_vs_eta->SetMarkerColor(C01);
  p_matching_eff_vs_eta->SetLineColor(C01);
  p_matching_eff_vs_eta->SetMinimum(0);
  p_matching_eff_vs_eta->SetMaximum(1.1);
  p_matching_eff_vs_eta->Draw("PE");

  c2->cd(2);
  p_matching_eff_vs_pt->SetMarkerStyle(20);
  p_matching_eff_vs_pt->SetMarkerColor(C01);
  p_matching_eff_vs_pt->SetLineColor(C01);
  p_matching_eff_vs_pt->SetMinimum(0);
  p_matching_eff_vs_pt->SetMaximum(1.1);
  p_matching_eff_vs_pt->Draw("PE");

  c2->cd(3);
  p_matching_eff_vs_nhits->SetMarkerStyle(20);
  p_matching_eff_vs_nhits->SetMarkerColor(C01);
  p_matching_eff_vs_nhits->SetLineColor(C01);
  p_matching_eff_vs_nhits->SetMinimum(0);
  p_matching_eff_vs_nhits->SetMaximum(1.1);
  p_matching_eff_vs_nhits->Draw("PE");

  c2->SaveAs("error_analysis_plots/plot2_matching_efficiency.pdf");

  // Plot 1.3: Resolution vs Matching Quality
  TCanvas *c3 = new TCanvas("c3", "Resolution vs Matching", 1200, 600);
  c3->Divide(2, 1);

  c3->cd(1);
  p_residual_rms_vs_nhits->SetMarkerStyle(20);
  p_residual_rms_vs_nhits->SetMarkerColor(C01);
  p_residual_rms_vs_nhits->SetLineColor(C01);
  p_residual_rms_vs_nhits->Draw("PE");

  c3->cd(2);
  p_resolution_vs_avg_residual->SetMarkerStyle(20);
  p_resolution_vs_avg_residual->SetMarkerColor(C01);
  p_resolution_vs_avg_residual->SetLineColor(C01);
  p_resolution_vs_avg_residual->Draw("PE");

  c3->SaveAs("error_analysis_plots/plot3_resolution_vs_matching.pdf");

  // Plot 2.1: Cluster Score Distributions (TRKPT, TRKPTZ, DNN)
  TCanvas *c4 = new TCanvas("c4", "Cluster Scores", 1800, 600);
  c4->Divide(3, 1);

  c4->cd(1);
  gPad->SetLogy();
  h_score_hs[TRKPT]->SetLineColor(C01);
  h_score_hs[TRKPT]->SetFillColorAlpha(C01, 0.3);
  h_score_hs[TRKPT]->SetLineWidth(2);
  h_score_hs[TRKPT]->Draw("HIST");

  h_score_pu[TRKPT]->SetLineColor(C02);
  h_score_pu[TRKPT]->SetFillColorAlpha(C02, 0.3);
  h_score_pu[TRKPT]->SetLineWidth(2);
  h_score_pu[TRKPT]->Draw("HIST SAME");

  TLegend *leg4 = new TLegend(0.6, 0.7, 0.88, 0.88);
  leg4->AddEntry(h_score_hs[TRKPT], "HS clusters", "f");
  leg4->AddEntry(h_score_pu[TRKPT], "PU clusters", "f");
  leg4->Draw();

  c4->cd(2);
  gPad->SetLogy();
  h_score_hs[TRKPTZ]->SetLineColor(C01);
  h_score_hs[TRKPTZ]->SetFillColorAlpha(C01, 0.3);
  h_score_hs[TRKPTZ]->SetLineWidth(2);
  h_score_hs[TRKPTZ]->Draw("HIST");

  h_score_pu[TRKPTZ]->SetLineColor(C02);
  h_score_pu[TRKPTZ]->SetFillColorAlpha(C02, 0.3);
  h_score_pu[TRKPTZ]->SetLineWidth(2);
  h_score_pu[TRKPTZ]->Draw("HIST SAME");
  leg4->Draw();

  c4->cd(3);
  gPad->SetLogy();
  h_score_hs[TESTML]->SetLineColor(C01);
  h_score_hs[TESTML]->SetFillColorAlpha(C01, 0.3);
  h_score_hs[TESTML]->SetLineWidth(2);
  h_score_hs[TESTML]->Draw("HIST");

  h_score_pu[TESTML]->SetLineColor(C02);
  h_score_pu[TESTML]->SetFillColorAlpha(C02, 0.3);
  h_score_pu[TESTML]->SetLineWidth(2);
  h_score_pu[TESTML]->Draw("HIST SAME");
  leg4->Draw();

  c4->SaveAs("error_analysis_plots/plot4_cluster_scores.pdf");

  // Plot 2.1b: DNN vs TRKPTZ comparison — rank, score ratio, time outcome
  TCanvas *c4b = new TCanvas("c4b", "DNN vs TRKPTZ Selection", 1800, 600);
  c4b->Divide(3, 1);

  // Panel 1: HS cluster rank comparison (TRKPTZ vs DNN)
  c4b->cd(1);
  double rank_total_trkptz = h_hs_cluster_rank->Integral();
  double rank_total_dnn    = h_hs_cluster_rank_dnn->Integral();
  TH1D *h_hs_rank_trkptz_norm = (TH1D*)h_hs_cluster_rank->Clone("h_hs_rank_trkptz_norm");
  TH1D *h_hs_rank_dnn_norm    = (TH1D*)h_hs_cluster_rank_dnn->Clone("h_hs_rank_dnn_norm");
  if (rank_total_trkptz > 0) h_hs_rank_trkptz_norm->Scale(1.0 / rank_total_trkptz);
  if (rank_total_dnn    > 0) h_hs_rank_dnn_norm->Scale(1.0 / rank_total_dnn);
  h_hs_rank_trkptz_norm->SetLineColor(C01);
  h_hs_rank_trkptz_norm->SetLineWidth(2);
  h_hs_rank_trkptz_norm->SetTitle("HS Cluster Rank;Rank;Fraction of Events");
  h_hs_rank_trkptz_norm->Draw("HIST");
  h_hs_rank_dnn_norm->SetLineColor(C02);
  h_hs_rank_dnn_norm->SetLineWidth(2);
  h_hs_rank_dnn_norm->SetLineStyle(2);
  h_hs_rank_dnn_norm->Draw("HIST SAME");
  TLegend *leg4b = new TLegend(0.45, 0.65, 0.88, 0.88);
  leg4b->AddEntry(h_hs_rank_trkptz_norm, "TRKPTZ", "l");
  leg4b->AddEntry(h_hs_rank_dnn_norm,    "DNN",    "l");
  leg4b->Draw();

  // Panel 2: Score gap comparison — large gap means the selector is confident
  c4b->cd(2);
  if (h_score_ratio->Integral() > 0)     h_score_ratio->Scale(1.0 / h_score_ratio->Integral());
  if (h_score_ratio_dnn->Integral() > 0) h_score_ratio_dnn->Scale(1.0 / h_score_ratio_dnn->Integral());
  h_score_ratio->SetLineColor(C01);
  h_score_ratio->SetLineWidth(2);
  h_score_ratio->SetTitle("Score Gap (1st vs 2nd);(Score_{1st} - Score_{2nd}) / Score_{1st} [%];Fraction of Events");
  h_score_ratio->Draw("HIST");
  h_score_ratio_dnn->SetLineColor(C02);
  h_score_ratio_dnn->SetLineWidth(2);
  h_score_ratio_dnn->SetLineStyle(2);
  h_score_ratio_dnn->Draw("HIST SAME");
  leg4b->Draw();

  // Panel 3: Time-based outcome comparison (TRKPTZ vs DNN) — shows efficiency difference
  c4b->cd(3);
  TH1D *h_to_trkptz_norm = (TH1D*)h_time_outcome->Clone("h_to_trkptz_norm");
  TH1D *h_to_dnn_norm    = (TH1D*)h_time_outcome_dnn->Clone("h_to_dnn_norm");
  if (h_to_trkptz_norm->Integral() > 0) h_to_trkptz_norm->Scale(1.0 / h_to_trkptz_norm->Integral());
  if (h_to_dnn_norm->Integral()    > 0) h_to_dnn_norm->Scale(1.0 / h_to_dnn_norm->Integral());
  h_to_trkptz_norm->SetLineColor(C01);
  h_to_trkptz_norm->SetLineWidth(2);
  h_to_trkptz_norm->SetBarWidth(0.4);
  h_to_trkptz_norm->SetBarOffset(0.1);
  h_to_trkptz_norm->SetTitle("Time-Based Outcome (TRKPTZ vs DNN);Outcome;Fraction");
  h_to_trkptz_norm->SetFillColorAlpha(C01, 0.4);
  h_to_trkptz_norm->Draw("BAR");
  h_to_dnn_norm->SetLineColor(C02);
  h_to_dnn_norm->SetLineWidth(2);
  h_to_dnn_norm->SetBarWidth(0.4);
  h_to_dnn_norm->SetBarOffset(0.5);
  h_to_dnn_norm->SetFillColorAlpha(C02, 0.4);
  h_to_dnn_norm->Draw("BAR SAME");
  leg4b->Draw();

  c4b->SaveAs("error_analysis_plots/plot4b_dnn_vs_trkptz_selection.pdf");

  // Plot 2.2: Selection Efficiency vs Event Properties
  TCanvas *c5 = new TCanvas("c5", "Selection Efficiency", 1800, 1200);
  c5->Divide(2, 2);

  c5->cd(1);
  p_selection_eff_vs_nclusters->SetMarkerStyle(20);
  p_selection_eff_vs_nclusters->SetMarkerColor(C01);
  p_selection_eff_vs_nclusters->SetMinimum(0);
  p_selection_eff_vs_nclusters->SetMaximum(1.1);
  p_selection_eff_vs_nclusters->Draw("PE");

  c5->cd(2);
  p_selection_eff_vs_njets->SetMarkerStyle(20);
  p_selection_eff_vs_njets->SetMarkerColor(C01);
  p_selection_eff_vs_njets->SetMinimum(0);
  p_selection_eff_vs_njets->SetMaximum(1.1);
  p_selection_eff_vs_njets->Draw("PE");

  c5->cd(3);
  p_selection_eff_vs_nhs->SetMarkerStyle(20);
  p_selection_eff_vs_nhs->SetMarkerColor(C01);
  p_selection_eff_vs_nhs->SetMinimum(0);
  p_selection_eff_vs_nhs->SetMaximum(1.1);
  p_selection_eff_vs_nhs->Draw("PE");

  c5->cd(4);
  p_selection_eff_vs_pufrac->SetMarkerStyle(20);
  p_selection_eff_vs_pufrac->SetMarkerColor(C01);
  p_selection_eff_vs_pufrac->SetMinimum(0);
  p_selection_eff_vs_pufrac->SetMaximum(1.1);
  p_selection_eff_vs_pufrac->Draw("PE");

  c5->SaveAs("error_analysis_plots/plot5_selection_efficiency.pdf");

  // Plot 2.3: Resolution Impact
  TCanvas *c6 = new TCanvas("c6", "Resolution Impact", 1800, 600);
  c6->Divide(3, 1);

  c6->cd(1);
  h_resolution_correct->SetLineColor(C08);
  h_resolution_correct->SetFillColorAlpha(C08, 0.3);
  h_resolution_correct->SetLineWidth(2);
  h_resolution_correct->Draw("HIST");

  c6->cd(2);
  h_resolution_wrong->SetLineColor(C02);
  h_resolution_wrong->SetFillColorAlpha(C02, 0.3);
  h_resolution_wrong->SetLineWidth(2);
  h_resolution_wrong->Draw("HIST");

  c6->cd(3);
  h_extra_error_from_selection->SetLineColor(C07);
  h_extra_error_from_selection->SetFillColorAlpha(C07, 0.3);
  h_extra_error_from_selection->SetLineWidth(2);
  h_extra_error_from_selection->Draw("HIST");

  c6->SaveAs("error_analysis_plots/plot6_resolution_impact.pdf");

  // Plot 3.1: Track Multiplicity and Success Rate
  TCanvas *c7 = new TCanvas("c7", "Track Multiplicity", 1800, 1200);
  c7->Divide(2, 2);

  c7->cd(1);
  h_n_forward_hs_valid->SetLineColor(C01);
  h_n_forward_hs_valid->SetFillColorAlpha(C01, 0.3);
  h_n_forward_hs_valid->Draw("HIST");

  c7->cd(2);
  eff_success_vs_nhs_valid->SetMarkerStyle(20);
  eff_success_vs_nhs_valid->SetMarkerColor(C01);
  eff_success_vs_nhs_valid->SetLineColor(C01);
  eff_success_vs_nhs_valid->Draw("APE");
  gPad->Update();
  auto graph = eff_success_vs_nhs_valid->GetPaintedGraph();
  graph->SetMinimum(0);
  graph->SetMaximum(1.1);
  gPad->Update();

  c7->cd(3);
  h2_success_vs_hs_pu->Draw("COLZ");

  c7->cd(4);
  h_track_loss->SetLineColor(C01);
  h_track_loss->SetFillColorAlpha(C01, 0.5);
  h_track_loss->SetBarWidth(0.8);
  h_track_loss->Draw("BAR");
  gPad->SetGridy();

  c7->SaveAs("error_analysis_plots/plot7_track_multiplicity.pdf");

  // Plot 3.2: Cluster Sizes
  TCanvas *c8 = new TCanvas("c8", "Cluster Sizes", 1200, 600);
  c8->Divide(2, 1);

  c8->cd(1);
  h_cluster_size_hs->SetLineColor(C01);
  h_cluster_size_hs->SetFillColorAlpha(C01, 0.3);
  h_cluster_size_hs->SetLineWidth(2);
  h_cluster_size_hs->Draw("HIST");

  h_cluster_size_pu->SetLineColor(C02);
  h_cluster_size_pu->SetFillColorAlpha(C02, 0.3);
  h_cluster_size_pu->SetLineWidth(2);
  h_cluster_size_pu->Draw("HIST SAME");

  TLegend *leg5 = new TLegend(0.6, 0.7, 0.88, 0.88);
  leg5->AddEntry(h_cluster_size_hs, "HS clusters", "f");
  leg5->AddEntry(h_cluster_size_pu, "PU clusters", "f");
  leg5->Draw();

  c8->cd(2);
  h_n_hs_in_hs_cluster->SetLineColor(C01);
  h_n_hs_in_hs_cluster->SetFillColorAlpha(C01, 0.3);
  h_n_hs_in_hs_cluster->SetLineWidth(2);
  h_n_hs_in_hs_cluster->Draw("HIST");

  c8->SaveAs("error_analysis_plots/plot8_cluster_sizes.pdf");

  // Plot 3.3: HS Track Fragmentation
  TCanvas *c9 = new TCanvas("c9", "HS Fragmentation", 1800, 1200);
  c9->Divide(2, 2);

  c9->cd(1);
  h_hs_time_spread->SetLineColor(C01);
  h_hs_time_spread->SetFillColorAlpha(C01, 0.3);
  h_hs_time_spread->Draw("HIST");

  c9->cd(2);
  h_hs_z_spread->SetLineColor(C01);
  h_hs_z_spread->SetFillColorAlpha(C01, 0.3);
  h_hs_z_spread->Draw("HIST");

  c9->cd(3);
  p_success_vs_time_spread->SetMarkerStyle(20);
  p_success_vs_time_spread->SetMarkerColor(C01);
  p_success_vs_time_spread->SetMinimum(0);
  p_success_vs_time_spread->SetMaximum(1.1);
  p_success_vs_time_spread->Draw("PE");

  c9->cd(4);
  p_success_vs_z_spread->SetMarkerStyle(20);
  p_success_vs_z_spread->SetMarkerColor(C01);
  p_success_vs_z_spread->SetMinimum(0);
  p_success_vs_z_spread->SetMaximum(1.1);
  p_success_vs_z_spread->Draw("PE");

  c9->SaveAs("error_analysis_plots/plot9_fragmentation.pdf");

  // Plot 3.4: Resolution vs Track Count
  TCanvas *c10 = new TCanvas("c10", "Resolution vs Tracks", 800, 600);
  p_resolution_vs_nhs->SetMarkerStyle(20);
  p_resolution_vs_nhs->SetMarkerColor(C01);
  p_resolution_vs_nhs->SetLineColor(C01);
  p_resolution_vs_nhs->Draw("PE");
  c10->SaveAs("error_analysis_plots/plot10_resolution_vs_tracks.pdf");

  // Plot 3.5: Additional diagnostic plots
  TCanvas *c11 = new TCanvas("c11", "Additional Diagnostics", 1800, 600);
  c11->Divide(3, 1);

  c11->cd(1);
  h_selection_outcome->Scale(1/h_selection_outcome->Integral());
  h_selection_outcome->SetLineColor(C01);
  h_selection_outcome->SetFillColorAlpha(C01, 0.3);
  h_selection_outcome->SetBarWidth(0.8);
  h_selection_outcome->Draw("BAR");

  c11->cd(2);
  h_hs_cluster_rank->SetLineColor(C01);
  h_hs_cluster_rank->SetFillColorAlpha(C01, 0.3);
  h_hs_cluster_rank->Draw("HIST");

  c11->cd(3);
  h_score_ratio->SetTitle("Score Gap (TRKPTZ);(Score_{1st} - Score_{2nd}) / Score_{1st} [%];Fraction of Events");
  h_score_ratio->SetLineColor(C01);
  h_score_ratio->SetFillColorAlpha(C01, 0.3);
  h_score_ratio->Draw("HIST");

  c11->SaveAs("error_analysis_plots/plot11_diagnostics.pdf");

  // Plot 4.1: HS Track Fragmentation Analysis
  TCanvas *c12 = new TCanvas("c12", "HS Fragmentation", 1800, 1200);
  c12->Divide(2, 2);

  c12->cd(1);
  h_fragmentation_category->Scale(1.0 / h_fragmentation_category->Integral());
  h_fragmentation_category->SetLineColor(C01);
  h_fragmentation_category->SetFillColorAlpha(C01, 0.3);
  h_fragmentation_category->SetBarWidth(0.8);
  h_fragmentation_category->Draw("BAR");
  gPad->SetGridy();

  c12->cd(2);
  h_n_clusters_with_hs->SetLineColor(C01);
  h_n_clusters_with_hs->SetFillColorAlpha(C01, 0.3);
  h_n_clusters_with_hs->Draw("HIST");

  c12->cd(3);
  h_hs_fragmentation_fraction->SetLineColor(C01);
  h_hs_fragmentation_fraction->SetFillColorAlpha(C01, 0.3);
  h_hs_fragmentation_fraction->Draw("HIST");

  c12->cd(4);
  h_n_clusters_passing->SetLineColor(C01);
  h_n_clusters_passing->SetFillColorAlpha(C01, 0.3);
  h_n_clusters_passing->Draw("HIST");

  c12->SaveAs("error_analysis_plots/plot12_hs_fragmentation.pdf");

  // Plot 4.2: Efficiency vs Fragmentation
  TCanvas *c13 = new TCanvas("c13", "Efficiency vs Fragmentation", 1200, 600);
  c13->Divide(2, 1);

  c13->cd(1);
  p_efficiency_vs_n_hs_clusters->SetMarkerStyle(20);
  p_efficiency_vs_n_hs_clusters->SetMarkerColor(C01);
  p_efficiency_vs_n_hs_clusters->SetLineColor(C01);
  p_efficiency_vs_n_hs_clusters->SetMinimum(0);
  p_efficiency_vs_n_hs_clusters->SetMaximum(1.1);
  p_efficiency_vs_n_hs_clusters->Draw("PE");
  gPad->SetGridy();

  c13->cd(2);
  p_efficiency_vs_fragmentation->SetMarkerStyle(20);
  p_efficiency_vs_fragmentation->SetMarkerColor(C01);
  p_efficiency_vs_fragmentation->SetLineColor(C01);
  p_efficiency_vs_fragmentation->SetMinimum(0);
  p_efficiency_vs_fragmentation->SetMaximum(1.1);
  p_efficiency_vs_fragmentation->Draw("PE");
  gPad->SetGridy();

  c13->SaveAs("error_analysis_plots/plot13_efficiency_vs_fragmentation.pdf");

  // Plot MC: Misclustering Study (TESTML vs TEST_MISCL)
  // Plot MC.1: DNN-selected cluster purity and score distributions
  TCanvas *cMC1 = new TCanvas("cMC1", "Misclustering: Purity & DNN Score", 1800, 600);
  cMC1->Divide(3, 1);

  cMC1->cd(1);
  h_dnn_selected_purity->SetLineColor(C01);
  h_dnn_selected_purity->SetFillColorAlpha(C01, 0.3);
  h_dnn_selected_purity->SetLineWidth(2);
  h_dnn_selected_purity->GetXaxis()->SetTitle("Purity of DNN-selected cluster");
  h_dnn_selected_purity->Draw("HIST");
  // Draw purity threshold line
  TLine *purity_line = new TLine(0.75, 0, 0.75, h_dnn_selected_purity->GetMaximum());
  purity_line->SetLineColor(C02);
  purity_line->SetLineWidth(2);
  purity_line->SetLineStyle(2);
  purity_line->Draw("SAME");
  TLegend *legMC1 = new TLegend(0.15, 0.75, 0.55, 0.88);
  legMC1->AddEntry(h_dnn_selected_purity, "DNN-selected cluster", "f");
  legMC1->AddEntry(purity_line, "Purity threshold (0.75)", "l");
  legMC1->Draw();

  cMC1->cd(2);
  gPad->SetLogy();
  h_dnn_score_pure->SetLineColor(C08);
  h_dnn_score_pure->SetFillColorAlpha(C08, 0.3);
  h_dnn_score_pure->SetLineWidth(2);
  h_dnn_score_pure->Draw("HIST");
  h_dnn_score_impure->SetLineColor(C02);
  h_dnn_score_impure->SetFillColorAlpha(C02, 0.3);
  h_dnn_score_impure->SetLineWidth(2);
  h_dnn_score_impure->Draw("HIST SAME");
  TLegend *legMC2 = new TLegend(0.15, 0.75, 0.6, 0.88);
  legMC2->AddEntry(h_dnn_score_pure,   "Pure clusters (purity>0.75)",   "f");
  legMC2->AddEntry(h_dnn_score_impure, "Impure clusters (purity<0.75)", "f");
  legMC2->Draw();

  cMC1->cd(3);
  h_dnn_reso_pure->SetLineColor(C08);
  h_dnn_reso_pure->SetFillColorAlpha(C08, 0.3);
  h_dnn_reso_pure->SetLineWidth(2);
  h_dnn_reso_pure->Draw("HIST");
  h_dnn_reso_impure->SetLineColor(C02);
  h_dnn_reso_impure->SetFillColorAlpha(C02, 0.3);
  h_dnn_reso_impure->SetLineWidth(2);
  h_dnn_reso_impure->Draw("HIST SAME");
  legMC2->Draw();

  cMC1->SaveAs("error_analysis_plots/plotMC1_misclustering_purity.pdf");

  // Plot MC.2: Efficiency comparison (TESTML vs TEST_MISCL) vs forward jets and PU fraction
  TCanvas *cMC2 = new TCanvas("cMC2", "Misclustering: Efficiency Comparison", 1200, 600);
  cMC2->Divide(2, 1);

  // Build TEfficiency objects from pass/total histograms
  TEfficiency *eff_dnn_fjet   = new TEfficiency(*h_dnn_eff_fjet_pass,   *h_dnn_eff_fjet_total);
  TEfficiency *eff_miscl_fjet = new TEfficiency(*h_miscl_eff_fjet_pass, *h_dnn_eff_fjet_total);
  eff_dnn_fjet->SetName("eff_dnn_vs_fjet");
  eff_dnn_fjet->SetTitle("DNN Efficiency vs N_{forward jets};N_{forward jets};Efficiency");
  eff_miscl_fjet->SetName("eff_miscl_vs_fjet");
  eff_miscl_fjet->SetTitle("DNN+Purity Efficiency vs N_{forward jets};N_{forward jets};Efficiency");

  TEfficiency *eff_dnn_pufrac   = new TEfficiency(*h_dnn_eff_pufrac_pass,   *h_dnn_eff_pufrac_total);
  TEfficiency *eff_miscl_pufrac = new TEfficiency(*h_miscl_eff_pufrac_pass, *h_dnn_eff_pufrac_total);
  eff_dnn_pufrac->SetName("eff_dnn_vs_pufrac");
  eff_miscl_pufrac->SetName("eff_miscl_vs_pufrac");

  cMC2->cd(1);
  eff_dnn_fjet->SetMarkerStyle(20);
  eff_dnn_fjet->SetMarkerColor(C01);
  eff_dnn_fjet->SetLineColor(C01);
  eff_dnn_fjet->Draw("APE");
  gPad->Update();
  eff_dnn_fjet->GetPaintedGraph()->SetMinimum(0);
  eff_dnn_fjet->GetPaintedGraph()->SetMaximum(1.1);
  gPad->Update();
  eff_miscl_fjet->SetMarkerStyle(24);
  eff_miscl_fjet->SetMarkerColor(C02);
  eff_miscl_fjet->SetLineColor(C02);
  eff_miscl_fjet->Draw("PE SAME");
  gPad->SetGridy();
  TLegend *legMC3 = new TLegend(0.35, 0.15, 0.88, 0.35);
  legMC3->AddEntry(eff_dnn_fjet,   "DNN (TESTML)", "pe");
  legMC3->AddEntry(eff_miscl_fjet, "DNN + purity>0.75 (TEST_MISCL)", "pe");
  legMC3->Draw();

  cMC2->cd(2);
  eff_dnn_pufrac->SetMarkerStyle(20);
  eff_dnn_pufrac->SetMarkerColor(C01);
  eff_dnn_pufrac->SetLineColor(C01);
  eff_dnn_pufrac->Draw("APE");
  gPad->Update();
  eff_dnn_pufrac->GetPaintedGraph()->SetMinimum(0);
  eff_dnn_pufrac->GetPaintedGraph()->SetMaximum(1.1);
  gPad->Update();
  eff_miscl_pufrac->SetMarkerStyle(24);
  eff_miscl_pufrac->SetMarkerColor(C02);
  eff_miscl_pufrac->SetLineColor(C02);
  eff_miscl_pufrac->Draw("PE SAME");
  gPad->SetGridy();
  legMC3->Draw();

  cMC2->SaveAs("error_analysis_plots/plotMC2_misclustering_efficiency.pdf");

  // Fill failure mode histogram
  h_failure_modes->SetBinContent(1, n_fail_wrong_selection);
  h_failure_modes->SetBinContent(2, n_fail_no_passing_cluster);
  h_failure_modes->SetBinContent(3, n_fail_matching_only);
  h_failure_modes->SetBinContent(4, n_fail_efficiency_only);
  int n_fail_irredeemable = n_fail_no_passing_cluster - n_fail_matching_only - n_fail_efficiency_only;
  h_failure_modes->SetBinContent(5, n_fail_irredeemable);

  // Fill insufficient tracks histogram
  h_insufficient_tracks->SetBinContent(1, n_insufficient_tracks_total);
  h_insufficient_tracks->SetBinContent(2, n_insufficient_tracks_fail);

  // Write all histograms to file
  fout->Write();
  fout->Close();

  std::cout << "\n========================================" << std::endl;
  std::cout << "Analysis complete!" << std::endl;
  std::cout << "Output file: hgtd_matching_analysis.root" << std::endl;
  std::cout << "Plots saved to: error_analysis_plots/" << std::endl;
  std::cout << "========================================\n" << std::endl;

  return 0;
}
