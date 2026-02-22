# VBF Selection Cuts - Loosening for Timing Studies

## Summary

The VBF H→invisible selection cuts have been **loosened** from the original ATLAS analysis values to increase event statistics for timing studies. The overall selection efficiency increased from **2.8% to 6.5%** (2.3x improvement).

## Modified Cuts

| Cut Variable | ATLAS Value | This Analysis | Efficiency Change |
|--------------|-------------|---------------|-------------------|
| Leading jet pT | > 80 GeV | **> 60 GeV** | 34.7% → **58.8%** |
| Subleading jet pT | > 50 GeV | **> 40 GeV** | (combined above) |
| Δη_jj | > 3.8 | **> 3.0** | 39.4% → **52.0%** |
| m_jj | > 800 GeV | **> 600 GeV** | 26.3% → **33.9%** |

## Unchanged Cuts

The following cuts remain at ATLAS values:

| Cut | Value | Efficiency |
|-----|-------|------------|
| Jet multiplicity | 2-4 jets, pT > 25 GeV | 43.3% |
| Δφ_jj | < 2.0 | 68.5% |
| Opposite hemispheres | η_j1 × η_j2 < 0 | 64.3% |
| All jets pT sum | > 140 GeV | 96.7% |
| Forward jet | ≥1 in HGTD (2.38<\|η\|<4.0) | 52.5% |
| FSR compatibility | C_i < 0.6, m_rel < 0.05 for j3,j4 | 2.2% |

## Overall Impact

**Before (ATLAS cuts, no FSR):**
- Full VBF selection: 28/1000 events (2.8%)

**After (Loosened cuts, no FSR):**
- Full VBF selection: 65/1000 events (6.5%)
- **2.3x increase in statistics**

**After (Loosened cuts + FSR compatibility):**
- Full VBF selection: 4/1000 events (0.4%)
- **FSR cuts reduce efficiency by 16x** (65 → 4 events)

## Rationale

### Why Loosen These Cuts?

The cuts were loosened to:
1. **Increase statistics** for timing resolution studies
2. **Maintain VBF topology** characteristics:
   - Forward jets in opposite hemispheres
   - Large rapidity separation
   - High dijet mass
3. **Preserve HGTD acceptance** requirement for timing measurements

### Why These Specific Values?

**Leading/Subleading jet pT (80/50 → 60/40 GeV):**
- VBF jets are typically high-pT
- 60/40 GeV still selects hard scattering processes
- Maintains trigger acceptance
- Provides 24% efficiency gain (35% → 59%)

**Δη_jj (3.8 → 3.0):**
- VBF events have large rapidity separation
- Δη > 3.0 still ensures forward topology
- Provides 13% efficiency gain (39% → 52%)

**m_jj (800 → 600 GeV):**
- VBF process produces high-mass dijet systems
- 600 GeV maintains discrimination from QCD
- Provides 8% efficiency gain (26% → 34%)

## Physics Implications

### What We Keep:
- ✅ Forward jet topology (Δη > 3, opposite hemispheres)
- ✅ VBF kinematics (high m_jj, forward jets)
- ✅ HGTD acceptance (at least one jet in 2.38 < |η| < 4.0)
- ✅ Hard scattering (jet pT thresholds, Σ pT > 140 GeV)

### What Changes:
- Slightly lower jet pT thresholds
- Slightly relaxed rapidity separation
- Lower dijet mass threshold

### Still Good For:
- ✅ Timing resolution studies in VBF-like events
- ✅ Forward jet timing performance
- ✅ HGTD acceptance studies
- ✅ Pileup rejection in forward region

### Not Suitable For:
- ❌ ATLAS VBF H→invisible search (requires ATLAS cuts)
- ❌ Direct comparison to published limits
- ❌ Precise signal/background discrimination

## Implementation Details

### Modified Functions

All changes in `src/clustering_structs.h`:

1. **`passLeadingJetPt()`** (lines 248-257):
   ```cpp
   // Original: leadingPt > 80.0 && subleadingPt > 50.0
   // Updated: leadingPt > 60.0 && subleadingPt > 40.0
   ```

2. **`passLargeDeltaEta()`** (lines 280-289):
   ```cpp
   // Original: default minDeltaEta = 3.8
   // Updated: default minDeltaEta = 3.0
   ```

3. **`passLargeDijetMass()`** (lines 306-310):
   ```cpp
   // Original: default minMass = 800.0
   // Updated: default minMass = 600.0
   ```

### Added Functions

4. **`calcCentrality(int jetIdx)`** (lines 342-364):
   ```cpp
   // Calculate C_i = exp(-4/(Δη_jj)^2 * (η_i - η_center)^2)
   // Measures how central a jet is between the two leading jets
   ```

5. **`calcRelativeMass(int jetIdx)`** (lines 366-386):
   ```cpp
   // Calculate m_rel_i = min{m_j1i, m_j2i} / m_jj
   // Relative mass of jet i compared to leading dijet system
   ```

6. **`passFSRCompatibility()`** (lines 388-411):
   ```cpp
   // Check C_i < 0.6 and m_rel < 0.05 for jets 3 and 4
   // Rejects events with FSR-like additional jets
   ```

### Documentation Updates

- `VBF_SELECTION_README.md`: Added "Cut Modifications" section documenting all changes
- `test_vbf_selection.cxx`: Updated output labels to show loosened cut values
- Debug output updated to reflect new thresholds

## Testing

Tested on 1000 VBF H→inv events:

**Example passing events:**
- Event 9: m_jj = 1075 GeV, Δη = 5.57, pT(j1/j2) = 71/62 GeV
- Event 66: m_jj = 657 GeV, Δη = 3.77, pT(j1/j2) = 169/58 GeV
- Event 90: m_jj = 3248 GeV, Δη = 7.67, pT(j1/j2) = 116/43 GeV

All passing events show VBF-like topology with:
- Large rapidity separation (Δη > 3.0)
- High dijet mass (m_jj > 600 GeV)
- Forward jets (at least one in HGTD acceptance)

## Reverting to ATLAS Cuts

To use original ATLAS cuts, change in `src/clustering_structs.h`:

```cpp
// Line ~256: passLeadingJetPt()
return (leadingPt > 80.0 && subleadingPt > 50.0);

// Line ~281: passLargeDeltaEta()
bool passLargeDeltaEta(double minDeltaEta = 3.8) const {

// Line ~306: passLargeDijetMass()
bool passLargeDijetMass(double minMass = 800.0) const {
```

Then rebuild:
```bash
cd build && make
```

## Reference

**Original cuts from:**
"Search for invisible Higgs-boson decays in events with vector-boson fusion signatures using 139 fb⁻¹ of proton-proton data recorded by the ATLAS experiment"

ATLAS Collaboration, JHEP 08 (2022) 104

https://doi.org/10.1007/JHEP08(2022)104

---

**Date:** 2026-01-23

**Purpose:** Increase statistics for HGTD timing resolution studies while maintaining VBF event topology
