# Pileup Removal Analysis

## Overview

The script `analyze_pileup_removal.py` analyzes the effect of a pileup removal algorithm that:
1. **Pre-selects** forward tracks passing quality cuts:
   - Track quality flag = true
   - Valid timing measurement (`Track_hasValidTime == 1`)
   - 2.38 < |η| < 4.00
   - 1 < pT < 30 GeV
2. **Pre-filters** tracks to only those within **3σ of the reconstructed primary vertex** (RecoVtx[0])
3. For remaining tracks, finds their closest reconstructed vertex in nσ
4. **Removes** tracks whose closest vertex is NOT the primary vertex (index ≠ 0)

This is different from a naive closest-vertex approach because it enforces strict quality cuts and the 3σ cut first, ensuring we only consider high-quality tracks that could plausibly come from the primary interaction region.

## Key Changes

### Algorithm Implementation

The pileup removal has **three stages**:

```python
# Stage 0: Track quality and kinematic selection
if not track_quality or not track_hasValidTime:
    skip_track()
if not (2.38 < |eta| < 4.00) or not (1 < pT < 30):
    skip_track()

# Stage 1: 3σ cut to primary vertex
nsigma_to_primary = |z_primary - z0_track| / sqrt(var_z0_track)
if nsigma_to_primary > 3.0:
    reject_track()  # Excluded from analysis
    continue

# Stage 2: Closest vertex association (only for tracks passing stages 0 and 1)
closest_vtx = argmin(|z_vtx - z0_track| / sqrt(var_z0_track))
if closest_vtx != 0:
    remove_track()  # Classified as pileup
else:
    keep_track()    # Classified as hard scatter
```

### Track Categories

The analysis now tracks **three** categories of forward tracks:

1. **Kept tracks**: Within 3σ of primary AND closest to primary vertex
   - These are considered hard scatter candidates

2. **Removed tracks**: Within 3σ of primary BUT closest to a different vertex
   - These are classified as pileup by the algorithm

3. **Rejected tracks**: Beyond 3σ of primary vertex
   - These are excluded by the pre-filter
   - Not considered in efficiency calculations

## Metrics Provided

### Track-Level Analysis

For each category (kept, removed, rejected), the script provides:
- Truth vertex association (is it really HS or PU?)
- pT distribution
- z0 distribution
- η distribution
- nσ distance to primary vertex

### Event-Level Statistics

- Number of forward tracks per event
- Breakdown: kept / removed / rejected per event
- Fraction of tracks removed per event
- Number of reconstructed vertices per event

### Performance Metrics

**Within the 3σ window only:**
- **HS Efficiency**: % of HS tracks that are kept
- **PU Rejection**: % of PU tracks that are removed

**Overall:**
- % of HS tracks rejected by 3σ cut
- % of PU tracks rejected by 3σ cut

## Usage

```bash
# Run on your data files
python3 analyze_pileup_removal.py \
    --input '/path/to/data/*.root' \
    --output pileup_removal_analysis.pdf

# Test on subset of events
python3 analyze_pileup_removal.py \
    --input '/path/to/data/*.root' \
    --max-events 10000 \
    --output test.pdf
```

## Output

The script generates a multi-page PDF with:

### Page 1: Truth Vertex Association
- Truth vertex index distributions
- HS vs PU breakdown (3 categories: kept, removed, rejected)
- Efficiency metrics
- Which RecoVtx index removed tracks were associated to

### Page 2: pT Distributions
- All tracks, HS tracks, PU tracks separately
- Mean pT comparison

### Page 3: z0 Distributions
- Spatial distributions for all categories
- z0 spread (RMS) comparison

### Page 4: nσ to Primary Vertex
- Shows the 3σ cut threshold
- Distribution for kept, removed, and rejected tracks
- Statistics on signal loss

### Page 5: Event-Level Statistics
- Forward tracks per event
- Removed vs kept scatter plot
- Fraction removed distribution
- Number of reconstructed vertices

### Page 6: η Distributions
- Pseudorapidity distributions
- Summary statistics table

### Page 7: ROC Curves
**Important**: These ROC curves show **alternative simple threshold-based algorithms**, NOT the actual closest-vertex algorithm. They answer: "What if we used simple cuts instead of closest-vertex association?"

- **ROC 1**: nσ to primary vertex (scan 0-3σ)
  - Keep tracks with nσ < threshold
  - Shows performance of a simple nσ cut (no closest-vertex logic)
  - Current algorithm (red star) uses closest-vertex, so may not lie on this curve

- **ROC 2**: |Δz| to primary vertex (scan 0-5mm)
  - Keep tracks with |Δz| < threshold
  - Physical distance threshold instead of statistical

- **ROC 3**: Track pT threshold (scan 1-30 GeV)
  - Keep tracks with pT > threshold
  - Shows if kinematic cuts provide discrimination

- **ROC 4**: Comparison plot
  - Overlays all three alternative algorithms
  - Red star = actual closest-vertex algorithm performance
  - Diagonal line = random classifier
  - Shows which simple discriminant would work best
  - Gap between curves and red star indicates benefit of closest-vertex logic

## Key Questions Answered

1. **How many HS tracks are lost by the 3σ cut?**
   - Check "Rejected (>3σ)" category in Page 1

2. **Among tracks within 3σ, how well does closest-vertex work?**
   - HS Efficiency and PU Rejection metrics in Page 1

3. **Are there kinematic biases?**
   - Check pT and η distributions in Pages 2 and 6

4. **Why are HS tracks being removed?**
   - Check nσ distributions in Page 4
   - If removed HS tracks have small nσ to primary, it suggests multiple vertices are very close

5. **Is the algorithm consistent across events?**
   - Check event-level distributions in Page 5

6. **Is 3σ the optimal threshold?**
   - Check ROC curves in Page 7
   - Shows the full efficiency vs rejection trade-off curve
   - Compare where the current 3σ operating point sits

7. **Which variable provides the best discrimination?**
   - ROC comparison plot (Page 7, bottom right)
   - Curves closer to the top-left corner are better
   - nσ-based cuts are typically better than physical distance or pT cuts

## Algorithm Context

This matches the track association logic in your C++ code (see `src/clustering_structs.h:99-113`):

```cpp
bool passBasicCuts() {
    // ... checks for jets and vertex quality ...

    // HS vertex must be within 2mm (MAX_VTX_DZ = 5.0) of truth
    if(std::abs(this->truthVtxZ[0] - this->recoVtxZ[0]) > MAX_VTX_DZ) {
        return false;
    }
    return true;
}
```

And the track selection in `countForwardTracks()` at line 149-171, which applies:
- Quality cuts
- pT: 1-30 GeV
- η: 2.38-4.00
- Optional timing validity

The 3σ cut (MAX_NSIGMA = 3.0) is a standard track-to-vertex association threshold used throughout the analysis.
