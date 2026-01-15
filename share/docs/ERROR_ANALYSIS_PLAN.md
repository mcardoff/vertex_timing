# Error Analysis Plan: Quantifying Algorithm Failure Modes

## Executive Summary

This document outlines methods to quantify three major sources of error in the HGTD vertex timing algorithm:
1. Wrong HGTD matching (incorrect track time assignments)
2. Wrong cluster selection (incorrect hard scatter identification)
3. Insufficient tracks (too few tracks for reliable clustering)

Each section describes diagnostic metrics, proposed plots, and implementation strategies using existing truth information.

---

## Error Source 1: Wrong HGTD Matching

### Problem Statement
HGTD hits are incorrectly matched to tracks, causing wrong time measurements that propagate through clustering and vertex reconstruction.

### Available Truth Information
- `Track_truthVtx_idx`: True vertex assignment for each track
- `Track_truthPart_idx`: True particle assignment (-1 for pileup)
- `TruthVtx_time`: True vertex times
- `TruthPart_prodVtx_time`: True particle production times
- `Track_nHGTDHits`: Number of HGTD hits contributing to track time

### Quantification Metrics

#### 1.1 Track-Level Time Residuals
**Metric**: `Δt_track = Track_time - Truth_time`

For each track with valid HGTD time:
```cpp
double truth_time;
if (trackToParticle[idx] != -1) {
    truth_time = particleT[trackToParticle[idx]];
} else if (trackToTruthvtx[idx] != -1) {
    truth_time = truthVtxTime[trackToTruthvtx[idx]];
} else {
    // Unmatched pileup - skip or flag
}
double residual = trackTime[idx] - truth_time;
```

**Histograms**:
- `h_time_residual_all`: All tracks (-500 to 500 ps, 200 bins)
- `h_time_residual_vs_nhits`: 2D profile (nHGTDHits vs residual width)
- `h_time_residual_vs_eta`: Track |η| vs residual (identify acceptance issues)
- `h_time_residual_vs_pt`: Track pT vs residual

#### 1.2 Matching Quality Flag
**Metric**: Define "well-matched" as `|Δt_track| < 3σ_track_time`

```cpp
bool is_well_matched = std::abs(residual) < 3.0 * trackTimeRes[idx];
```

**Histograms**:
- `h_matching_efficiency_vs_eta`: Fraction of tracks with good matching vs η
- `h_matching_efficiency_vs_pt`: vs pT
- `h_matching_efficiency_vs_nhits`: vs number of HGTD hits

#### 1.3 Cluster Contamination from Mismatched Tracks
**Metric**: For each cluster, calculate fraction of tracks with large residuals

```cpp
int n_good = 0, n_bad = 0;
for (auto trk : cluster.trackIndices) {
    double residual = /* calculate as above */;
    if (std::abs(residual) < 3.0 * trackTimeRes[trk]) n_good++;
    else n_bad++;
}
double contamination = (double)n_bad / (n_good + n_bad);
```

**Histograms**:
- `h_cluster_contamination`: Contamination fraction (0 to 1, 50 bins)
- `h_resolution_vs_contamination`: 2D profile (contamination vs vertex time resolution)
- `h_purity_vs_contamination`: 2D profile (contamination vs cluster purity)

#### 1.4 Event-Level Matching Quality
**Metric**: Average track residual in event (only HS tracks)

```cpp
double sum_residual = 0;
int n_hs_tracks = 0;
for (int trk = 0; trk < trackTime.GetSize(); trk++) {
    if (trackToTruthvtx[trk] == 0 && trackTimeValid[trk] == 1) {
        sum_residual += std::abs(trackTime[trk] - truthVtxTime[0]);
        n_hs_tracks++;
    }
}
double avg_residual = (n_hs_tracks > 0) ? sum_residual / n_hs_tracks : -1;
```

**Histograms**:
- `h_avg_hs_residual`: Average residual per event (0 to 200 ps, 100 bins)
- `h_algo_success_vs_residual`: Algorithm pass rate vs avg residual
- `h_resolution_vs_avg_residual`: Achieved resolution vs avg residual

### Proposed Plots

#### Plot 1.1: Track Time Residual Distribution
- **Type**: Overlaid histograms
- **X-axis**: Time residual (ps)
- **Y-axis**: Fraction of tracks (normalized)
- **Curves**:
  - All tracks (black)
  - Tracks from HS (blue)
  - Tracks from PU (red)
  - Well-matched tracks (green dashed)
  - Mis-matched tracks (orange dashed)
- **Stats**: Show mean and RMS in legend
- **Purpose**: Quantify baseline matching quality

#### Plot 1.2: Matching Efficiency vs Kinematics
- **Type**: Multi-panel efficiency plot (2x2 grid)
- **Panels**:
  1. Efficiency vs |η| (2.4 to 4.0)
  2. Efficiency vs pT (1 to 30 GeV)
  3. Efficiency vs nHGTDHits (0 to 10)
  4. Efficiency vs vertex dz (0 to 500 mm)
- **Purpose**: Identify kinematic regions with poor matching

#### Plot 1.3: Algorithm Performance vs Matching Quality
- **Type**: 2D heatmap or profile
- **X-axis**: Average HS track residual (ps)
- **Y-axis**: Vertex time resolution (ps)
- **Color**: Event density or mean resolution
- **Overlays**:
  - Vertical line at 3σ expected (e.g., 75 ps)
  - Horizontal line at target resolution (e.g., 30 ps)
- **Purpose**: Directly correlate matching errors with resolution degradation

#### Plot 1.4: Residual Tails Analysis
- **Type**: Cumulative distribution
- **X-axis**: |Time residual| (ps)
- **Y-axis**: Fraction of tracks with residual < X
- **Curves**:
  - All tracks
  - HS tracks only
  - PU tracks only
- **Reference lines**: 1σ, 2σ, 3σ expected resolution
- **Purpose**: Quantify outlier rate

#### Plot 1.5: Resolution vs nHGTDHits
- **Type**: Profile plot with error bars
- **X-axis**: Number of HGTD hits
- **Y-axis**: RMS of time residual (ps)
- **Expected**: Should scale as `σ₀ / √(nHits)`
- **Overlay**: Fit to `A / √x + B` to extract baseline resolution
- **Purpose**: Verify resolution scaling and detect anomalies

---

## Error Source 2: Wrong Cluster Selection

### Problem Statement
The algorithm identifies clusters but selects the wrong one as the hard scatter due to insufficient discriminating power from pT and z information alone.

### Available Truth Information
- `TruthVtx_isHS`: Flags hard scatter vertex
- `Track_truthVtx_idx`: True vertex assignment
- Track pT and z0 for weighting

### Quantification Metrics

#### 2.1 Cluster-Level Truth Matching
**Metric**: Match each cluster to truth vertex with highest pT contribution

```cpp
// For each cluster
std::map<int, double> vtx_pt_map;  // vtxIdx -> sum pT
for (auto trk : cluster.trackIndices) {
    int truth_vtx = trackToTruthvtx[trk];
    if (truth_vtx == -1) continue;
    vtx_pt_map[truth_vtx] += trackPt[trk];
}

int matched_vtx = -1;
double max_pt = 0;
for (auto [vtx, pt] : vtx_pt_map) {
    if (pt > max_pt) {
        max_pt = pt;
        matched_vtx = vtx;
    }
}
bool cluster_is_hs = (matched_vtx == 0);  // Assuming HS is vertex 0
```

#### 2.2 Selection Confusion Matrix
**Metric**: Track algorithm selection decisions

For each event with ≥2 clusters:
```cpp
bool hs_cluster_exists = false;
bool hs_cluster_selected = false;
int selected_idx = /* from chooseCluster(...) */;

for (auto& cluster : clusters) {
    // Match to truth as above
    if (cluster_matches_hs) {
        hs_cluster_exists = true;
        if (&cluster == &clusters[selected_idx]) {
            hs_cluster_selected = true;
        }
    }
}

// Fill confusion matrix
if (hs_cluster_exists && hs_cluster_selected) n_true_positive++;
if (hs_cluster_exists && !hs_cluster_selected) n_false_negative++;
if (!hs_cluster_exists && hs_cluster_selected) n_false_positive++;
if (!hs_cluster_exists && !hs_cluster_selected) n_true_negative++;
```

**Metrics**:
- **Efficiency**: TP / (TP + FN) - fraction of HS clusters correctly selected
- **Purity**: TP / (TP + FP) - fraction of selections that are actually HS
- **False Selection Rate**: FN / (TP + FN) - how often we pick wrong cluster

#### 2.3 Score Discrimination Power
**Metric**: Separation between HS and PU cluster scores

For each scoring method (TRKPT, TRKPTZ, etc.):
```cpp
// Compute for all clusters in event
std::vector<double> hs_scores, pu_scores;
for (auto& cluster : clusters) {
    double score = cluster.scores[TRKPTZ];  // or other Score enum
    if (cluster_matches_hs) hs_scores.push_back(score);
    else pu_scores.push_back(score);
}

// Separation metric
double hs_max = *std::max_element(hs_scores.begin(), hs_scores.end());
double pu_max = *std::max_element(pu_scores.begin(), pu_scores.end());
double separation = (hs_max - pu_max) / hs_max;  // Relative separation
```

**Histograms**:
- `h_score_separation_{method}`: Separation distribution for each method
- `h_hs_score_{method}`: HS cluster scores
- `h_pu_score_{method}`: PU cluster scores

#### 2.4 Multi-Cluster Event Characterization
**Metric**: Properties of events with multiple clusters

```cpp
struct MultiClusterEvent {
    int n_clusters_total;
    int n_clusters_hs;
    int n_clusters_pu;
    double hs_cluster_rank;  // Rank by chosen score (1 = highest)
    double second_best_score_ratio;  // score[2nd] / score[1st]
};
```

**Histograms**:
- `h_n_clusters`: Number of clusters per event
- `h_hs_rank`: Rank of HS cluster by pT score (1, 2, 3, ...)
- `h_score_ratio`: Ratio of 2nd to 1st highest score (ambiguity measure)

#### 2.5 Wrong Selection Impact
**Metric**: Resolution degradation when wrong cluster selected

```cpp
double true_hs_time = truthVtxTime[0];
double selected_cluster_time = clusters[selected_idx].values[0];
double hs_cluster_time = /* time of HS cluster if exists */;

double error_wrong_selection = selected_cluster_time - true_hs_time;
double error_correct_selection = hs_cluster_time - true_hs_time;
double extra_error = std::abs(error_wrong_selection) - std::abs(error_correct_selection);
```

**Histograms**:
- `h_extra_error_from_selection`: Additional error from wrong selection
- `h_time_diff_hs_vs_selected`: |t_HS - t_selected| when HS exists but not selected

### Proposed Plots

#### Plot 2.1: Cluster Score Distributions
- **Type**: Overlaid histograms with fills
- **X-axis**: Cluster score (log scale if needed)
- **Y-axis**: Fraction of clusters (normalized)
- **Curves**:
  - HS clusters (blue, filled)
  - PU clusters (red, filled with transparency)
- **Separate panels** for each scoring method:
  - TRKPT
  - TRKPTZ
  - CALO60
  - CALO90
- **Stats**: Show separation power metric
- **Purpose**: Visualize discriminating power of each method

#### Plot 2.2: Selection Efficiency vs Event Properties
- **Type**: Profile plots (2x2 or 2x3 grid)
- **X-axes**:
  - Number of clusters in event
  - Number of forward jets
  - Pileup fraction
  - |ΔZ| between HS and nearest PU vertex
  - Number of HS tracks
  - Score ratio (2nd best / best)
- **Y-axis**: Selection efficiency (fraction correct)
- **Error bars**: Binomial errors
- **Purpose**: Identify conditions where selection fails

#### Plot 2.3: Confusion Matrix Heatmap
- **Type**: 2D heatmap (normalized)
- **Rows**: True state (HS exists, HS absent)
- **Columns**: Selected state (HS selected, PU selected, None selected)
- **Color**: Event fraction (0 to 1)
- **Annotations**: Print percentages in each cell
- **Separate matrices** for each scoring method
- **Purpose**: Quantify selection accuracy

#### Plot 2.4: ROC-like Curve for Cluster Selection
- **Type**: Scatter or line plot
- **X-axis**: False selection rate (wrong cluster chosen)
- **Y-axis**: True selection rate (correct cluster chosen)
- **Points**: Different score thresholds or methods
- **Diagonal line**: Random selection
- **Labels**: TRKPT, TRKPTZ, CALO60, etc.
- **Purpose**: Compare discriminating power of methods

#### Plot 2.5: Resolution vs Selection Quality
- **Type**: Box plots or violin plots
- **X-axis**: Selection outcome categories
  - "Correct HS"
  - "Wrong cluster (HS exists)"
  - "No clusters formed"
- **Y-axis**: Vertex time resolution (ps)
- **Purpose**: Quantify resolution penalty for wrong selection

#### Plot 2.6: HS Cluster Rank Distribution
- **Type**: Bar chart
- **X-axis**: Rank of HS cluster by score (1st, 2nd, 3rd, ...)
- **Y-axis**: Fraction of events
- **Separate bars**: Different scoring methods
- **Purpose**: How often is HS cluster the highest-ranked?

#### Plot 2.7: Score Separation vs ΔZ
- **Type**: 2D profile or heatmap
- **X-axis**: |Z_HS - Z_nearest_PU| (mm)
- **Y-axis**: Score separation (HS_score - max_PU_score) / HS_score
- **Color**: Mean or event density
- **Purpose**: Check if nearby PU vertices in Z make selection harder

---

## Error Source 3: Insufficient Tracks

### Problem Statement
Too few tracks in the event (especially forward HS tracks) lead to small, fragmented clusters that cannot be reliably identified as the hard scatter.

### Available Truth Information
- `Track_truthVtx_idx`: Count HS vs PU tracks
- Track kinematic cuts: |η| > 2.38, pT > 1 GeV, quality flag
- `Track_hasValidTime`: Valid HGTD measurement

### Quantification Metrics

#### 3.1 Track Multiplicity Distributions
**Metric**: Count tracks in various categories per event

```cpp
struct TrackCounts {
    int n_all_tracks;           // Total in event
    int n_forward_tracks;       // |η| > 2.38
    int n_forward_hs_tracks;    // Forward + from HS
    int n_forward_pu_tracks;    // Forward + from PU
    int n_forward_valid_time;   // Forward + valid HGTD time
    int n_forward_hs_valid;     // Forward + HS + valid time
};
```

**Histograms**:
- `h_n_forward_tracks`: Forward track multiplicity (0 to 50)
- `h_n_forward_hs_tracks`: Forward HS tracks (0 to 20)
- `h_n_forward_hs_valid`: Forward HS with valid time (0 to 20)
- `h_hs_time_efficiency`: Fraction of HS tracks with valid time

#### 3.2 Clustering Success vs Track Count
**Metric**: Algorithm outcome as function of track multiplicity

```cpp
enum AlgoOutcome {
    NO_CLUSTERS_FORMED,
    CLUSTERS_FORMED_HS_SELECTED,
    CLUSTERS_FORMED_WRONG_SELECTED,
    INSUFFICIENT_SEPARATION  // Multiple clusters, ambiguous scores
};

AlgoOutcome outcome = /* determine from analysis */;
```

**Histograms**:
- `h_outcome_vs_nhs`: 2D (n_hs_tracks vs outcome)
- `h_outcome_vs_nhs_valid`: 2D (n_hs_tracks_with_time vs outcome)
- `h_success_rate_vs_nhs`: Profile (n_hs_tracks vs success fraction)

#### 3.3 Cluster Formation Thresholds
**Metric**: Minimum track count for viable clustering

```cpp
// For events that pass/fail algorithm
int n_hs_pass_min = 1000, n_hs_fail_max = 0;
if (algo_passed) {
    n_hs_pass_min = std::min(n_hs_pass_min, n_forward_hs_valid);
} else {
    n_hs_fail_max = std::max(n_hs_fail_max, n_forward_hs_valid);
}
```

**Analysis**:
- Find 50% efficiency point: `n_hs_tracks` where success rate = 0.5
- Find 90% efficiency point
- Plot efficiency curve with threshold markers

#### 3.4 Cluster Size Distribution
**Metric**: Number of tracks per cluster

```cpp
for (auto& cluster : clusters) {
    int n_tracks = cluster.trackIndices.size();
    int n_hs_tracks = 0;
    for (auto trk : cluster.trackIndices) {
        if (trackToTruthvtx[trk] == 0) n_hs_tracks++;
    }
    // Fill histograms
}
```

**Histograms**:
- `h_cluster_size_all`: Tracks per cluster, all clusters
- `h_cluster_size_hs`: Tracks per cluster, HS cluster
- `h_cluster_size_pu`: Tracks per cluster, PU clusters
- `h_n_hs_in_hs_cluster`: HS tracks in the HS cluster (key metric!)

#### 3.5 Track Loss Along Selection Chain
**Metric**: Track attrition at each selection stage

```cpp
// Start with all tracks
int n_start = trackPt.GetSize();

// After forward cuts
int n_after_eta = /* count |η| > 2.38 */;

// After pT cut
int n_after_pt = /* count pT > 1 GeV */;

// After quality cut
int n_after_quality = /* count quality == true */;

// After vertex association
int n_after_assoc = /* count passPrimVertexAssociation */;

// After pileup removal
int n_after_pu = /* count passing pileup removal */;

// After time validity
int n_after_time = /* count hasValidTime == 1 */;

// After calo filtering (if applied)
int n_after_calo = /* count passing calo filter */;
```

**Visualization**: Sankey diagram or bar chart showing attrition

#### 3.6 HS Track Fragmentation
**Metric**: How spread out are HS tracks in time/space?

```cpp
// For HS tracks only
std::vector<double> hs_times, hs_z0s;
for (int trk = 0; trk < trackTime.GetSize(); trk++) {
    if (trackToTruthvtx[trk] == 0 && trackTimeValid[trk] == 1) {
        hs_times.push_back(trackTime[trk]);
        hs_z0s.push_back(trackZ0[trk]);
    }
}

// Compute spread
double time_rms = calculate_rms(hs_times);
double z_rms = calculate_rms(hs_z0s);

// Correlate with clustering success
```

**Histograms**:
- `h_hs_time_spread`: RMS of HS track times (0 to 200 ps)
- `h_hs_z_spread`: RMS of HS track z0 (0 to 50 mm)
- `h_success_vs_time_spread`: 2D profile
- `h_success_vs_z_spread`: 2D profile

### Proposed Plots

#### Plot 3.1: Track Multiplicity and Success Rate
- **Type**: Combined bar + line plot
- **X-axis**: Number of forward HS tracks with valid time
- **Y-axis (left)**: Number of events (bars)
- **Y-axis (right)**: Algorithm success rate (line with markers)
- **Vertical lines**: Threshold markers (50%, 90% efficiency)
- **Shaded regions**: "Insufficient tracks" zone
- **Purpose**: Show how track count affects performance

#### Plot 3.2: Track Loss Waterfall Chart
- **Type**: Horizontal bar chart (stacked or cascading)
- **Categories** (top to bottom):
  - All tracks in event
  - After |η| > 2.38 cut
  - After pT > 1 GeV cut
  - After quality cut
  - After vertex association
  - After pileup removal
  - After time validity requirement
  - After calo filter (if used)
  - HS tracks remaining
- **X-axis**: Average number of tracks
- **Colors**: HS tracks (blue), PU tracks (red), Lost tracks (gray)
- **Purpose**: Visualize where tracks are lost

#### Plot 3.3: Clustering Success vs Track Counts (2D)
- **Type**: 2D heatmap
- **X-axis**: Number of forward HS tracks with valid time
- **Y-axis**: Number of forward PU tracks with valid time
- **Color**: Algorithm success rate (0 to 1, diverging colormap)
- **Contours**: 0.5, 0.7, 0.9 success rate
- **Purpose**: Show interplay between HS and PU track counts

#### Plot 3.4: Cluster Size Distributions (Stacked)
- **Type**: Stacked histograms
- **X-axis**: Number of tracks in cluster
- **Y-axis**: Number of clusters
- **Stack layers**:
  - HS cluster (blue)
  - PU clusters (red)
- **Separate panels**:
  - Events that pass algorithm
  - Events that fail algorithm
- **Purpose**: Show typical cluster sizes for success/failure

#### Plot 3.5: HS Track Fragmentation Impact
- **Type**: 2D scatter with color
- **X-axis**: RMS of HS track times (ps)
- **Y-axis**: RMS of HS track z0 (mm)
- **Color**: Algorithm success (blue) vs failure (red)
- **Size**: Number of HS tracks
- **Purpose**: Check if spread in time/space causes failure

#### Plot 3.6: Efficiency Curves vs Cuts
- **Type**: Multiple efficiency curves
- **X-axis**: Cut threshold
- **Y-axis**: Algorithm efficiency (0 to 1)
- **Curves**:
  - Min HS tracks required (2, 3, 4, ...)
  - Max track time spread (50, 100, 150 ps)
  - Min cluster size (2, 3, 4, 5 tracks)
- **Purpose**: Optimize selection criteria for track requirements

#### Plot 3.7: Resolution vs Track Count
- **Type**: Profile plot with error bands
- **X-axis**: Number of HS tracks in selected cluster
- **Y-axis**: Vertex time resolution (ps)
- **Error bands**: Standard error of mean
- **Expected**: Resolution should improve with √N tracks
- **Overlay**: Fit to `A / √N + B` (statistical + systematic terms)
- **Purpose**: Quantify resolution dependence on statistics

---

## Implementation Strategy

### Phase 1: Instrumentation (Week 1)
1. Create new analysis executable: `error_source_analysis.cxx`
2. Add histogram booking in initialization
3. Implement per-event metrics calculation
4. Fill histograms during event loop
5. Truth-match clusters to vertices

### Phase 2: Validation (Week 1-2)
1. Run on small sample (1k events) to debug
2. Verify truth matching logic
3. Check histogram ranges and binning
4. Cross-check with existing `hgtd_matching.cxx` results

### Phase 3: Production (Week 2)
1. Run on full sample
2. Generate all plots with ROOT macros
3. Create summary statistics tables

### Phase 4: Interpretation (Week 2-3)
1. Analyze correlations between error sources
2. Identify dominant failure mode in different kinematic regions
3. Propose algorithmic improvements based on findings

### File Structure
```
error_analysis/
├── error_source_analysis.cxx        # Main analysis
├── error_analysis_plots.cxx         # Plotting macros
├── error_histograms.h               # Histogram definitions
└── error_metrics.h                  # Metric calculation functions
```

---

## Expected Outcomes

### Quantitative Results

1. **HGTD Matching**:
   - Baseline matching efficiency: __%
   - Residual RMS by η region: __ to __ ps
   - Fraction of events with ≥1 badly matched track: __%
   - Correlation: ΔResolution ∝ __ × AvgResidual

2. **Cluster Selection**:
   - Selection efficiency (HS present): __%
   - False selection rate (wrong cluster): __%
   - Score separation (HS vs PU): __ (TRKPT) vs __ (TRKPTZ)
   - Rank of HS cluster: __% are 1st, __% are 2nd, etc.

3. **Insufficient Tracks**:
   - 50% efficiency point: __ forward HS tracks
   - 90% efficiency point: __ forward HS tracks
   - Typical HS track loss: __% (forward selection) + __% (time validity)
   - Fragmentation: RMS spread = __ ps (time), __ mm (z)

### Qualitative Insights

- Which error source dominates in which kinematic regime?
- Are errors correlated? (e.g., mismatching worse when few tracks?)
- Where should algorithmic improvements focus?

---

## Next Steps

1. **Review this plan**: Confirm metrics and plots address your needs
2. **Prioritize**: Which error source to analyze first?
3. **Implementation**: Write `error_source_analysis.cxx` skeleton
4. **Iterate**: Refine based on initial results

---

## Notes

- All truth information is available in ntuples (verified in `clustering_structs.h`)
- Existing infrastructure (`Cluster` class, `AnalysisObj`, etc.) can be reused
- Plots should be publication-quality (ROOT TStyle, proper labels, legends)
- Consider separate analysis for different pileup scenarios (μ=200 vs μ=140)
