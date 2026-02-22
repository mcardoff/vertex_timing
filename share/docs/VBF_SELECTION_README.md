# VBF H→Invisible Signal Region Selection

## Overview

This document describes the VBF (Vector Boson Fusion) H→invisible signal region selection functions implemented in the vertex timing analysis project. The selection criteria are based on the ATLAS VBF H→invisible search published in **JHEP08(2022)104**, with **loosened cuts** to increase statistics for timing studies.

## Implementation

All VBF selection functions are implemented as member functions of the `BranchPointerWrapper` struct in `src/clustering_structs.h` (lines 223-362).

### Important Assumption

**Jets are pT-sorted in descending order** in the ROOT ntuples (verified with `check_jet_sorting.cxx` - 100% of 1000 events tested).
- Leading jet = `topoJetPt[0]` (highest pT)
- Subleading jet = `topoJetPt[1]` (second highest pT)

No additional sorting is required.

### Available Functions

#### Individual Selection Cuts

1. **`countJetsAbovePt(double ptThreshold = 25.0)`**
   - Counts jets with pT above the specified threshold
   - Default: 25 GeV

2. **`passJetMultiplicity()`**
   - Requires 2, 3, or 4 jets with pT > 25 GeV
   - Efficiency: ~43.3%

3. **`passLeadingJetPt()`**
   - Leading jet pT > 60 GeV (ATLAS: 80 GeV)
   - Subleading jet pT > 40 GeV (ATLAS: 50 GeV)
   - **Loosened for timing studies**
   - Efficiency: ~58.8% (was ~34.7% with ATLAS cuts)

4. **`passJetDeltaPhi()`**
   - Two leading jets not back-to-back: Δφ_jj < 2.0
   - Efficiency: ~68.5%

5. **`passOppositeHemispheres()`**
   - VBF topology: η_j1 × η_j2 < 0
   - Efficiency: ~64.3%

6. **`passLargeDeltaEta(double minDeltaEta = 3.0)`**
   - Large pseudorapidity separation: |Δη_jj| > 3.0 (ATLAS: 3.8)
   - **Loosened for timing studies**
   - Efficiency: ~52.0% (was ~39.4% with ATLAS cuts)

7. **`calcDijetMass()`**
   - Calculates invariant mass of two leading jets (in GeV)
   - Returns -1 if fewer than 2 jets

8. **`passLargeDijetMass(double minMass = 600.0)`**
   - Large dijet mass: m_jj > 600 GeV (ATLAS: 800 GeV)
   - **Loosened for timing studies**
   - Efficiency: ~33.9% (was ~26.3% with ATLAS cuts)

9. **`calcAllJetsPtSum()`**
   - Calculates scalar sum of all jet pT

10. **`passAllJetsPtSum(double minPtSum = 140.0)`**
    - Sum of all jet pT > 140 GeV
    - Efficiency: ~96.7%

11. **`passForwardJet()`**
    - At least one of the two leading jets in HGTD acceptance (2.38 < |η| < 4.0)
    - Ensures VBF events have forward jet activity for timing studies
    - Efficiency: ~52.5%

12. **`calcCentrality(int jetIdx)`**
    - Calculates centrality C_i for jet i between two leading jets
    - C_i = exp(-4/(Δη_jj)² × (η_i - η_center)²)
    - Measures how central a jet is in rapidity
    - Returns value between 0 (edge) and 1 (center)

13. **`calcRelativeMass(int jetIdx)`**
    - Calculates relative mass m_rel for jet i
    - m_rel_i = min{m_j1i, m_j2i} / m_jj
    - Compares dijet masses involving additional jet to leading dijet mass
    - Returns -1 if calculation fails

14. **`passFSRCompatibility()`**
    - Requires C_i < 0.6 and m_rel < 0.05 for jets 3 and 4 (if they exist)
    - Rejects events with FSR-like (Final State Radiation) additional jets
    - Passes automatically if only 2 jets in event
    - Efficiency: ~2.2%

#### Combined Selection

15. **`passVBFSignalRegion()`**
    - Applies all VBF cuts simultaneously (including forward jet and FSR requirements)
    - Uses **loosened cuts** compared to ATLAS analysis
    - Overall efficiency: **~0.4%** (was 6.5% without FSR, 2.8% with ATLAS cuts)
    - Prints detailed cut flow when `DEBUG = true`

## Cut Modifications for Timing Studies

The following cuts have been **loosened** from the original ATLAS VBF H→invisible analysis to increase statistics for timing studies while maintaining VBF-like topology:

| Cut | ATLAS Value | This Analysis | Reason |
|-----|-------------|---------------|--------|
| Leading jet pT | > 80 GeV | > 60 GeV | Increase acceptance |
| Subleading jet pT | > 50 GeV | > 40 GeV | Increase acceptance |
| Δη_jj | > 3.8 | > 3.0 | Increase acceptance while maintaining forward topology |
| m_jj | > 800 GeV | > 600 GeV | Increase acceptance while selecting VBF-like events |

**Unchanged cuts:**
- Jet multiplicity: 2-4 jets with pT > 25 GeV
- Δφ_jj: < 2.0
- Opposite hemispheres: η_j1 × η_j2 < 0
- All jets pT sum: > 140 GeV
- Forward jet requirement: ≥1 leading jet in HGTD (2.38 < |η| < 4.0)
- FSR compatibility: C_i < 0.6, m_rel < 0.05 for jets 3 and 4

## Excluded Requirements

The following cuts from JHEP08(2022)104 are **not implemented** due to unavailable variables in the ntuples:

- ❌ Lepton veto (no lepton candidates)
- ❌ Photon veto
- ❌ JVT (Jet Vertex Tagger) requirement
- ❌ fJVT (forward JVT) requirement
- ❌ b-tagging veto (≤1 b-tagged jet)
- ❌ FSR compatibility cuts (C_i < 0.6, m^rel_i < 0.05)
- ❌ Missing transverse energy: E^miss_T > 160 GeV
- ❌ Soft term requirement: soft-term E^miss_T < 20 GeV

## Usage Example

### In Analysis Code

```cpp
#include "src/clustering_functions.h"
#include "src/event_processing.h"

// Setup reader and branches
TTreeReader reader(&chain);
BranchPointerWrapper branch(reader);

// Event loop
while (reader.Next()) {

    // Check full VBF selection
    if (!branch.passVBFSignalRegion()) continue;

    // Or check individual cuts
    if (branch.passJetMultiplicity() &&
        branch.passLeadingJetPt() &&
        branch.calcDijetMass() > 1000.0) {

        // Process VBF-like event
        double mjj = branch.calcDijetMass();
        std::cout << "VBF event with m_jj = " << mjj << " GeV" << std::endl;
    }
}
```

### Testing

A dedicated test executable is provided:

```bash
cd /Users/mcard/project/vertex_timing

# Build
cd build
cmake ..
make test_vbf_selection

# Run on first 1000 events of file 000001
cd ..
./build/test_vbf_selection 000001 1000
```

Output example:
```
========================================
VBF H->Invisible Selection Test
========================================

Event 519 passes VBF selection:
  Number of jets: 2
  Leading jet pT: 92.7118 GeV
  Subleading jet pT: 47.9668 GeV
  Dijet mass: 763.892 GeV
  Delta eta_jj: 4.88193
  All jets pT sum: 140.679 GeV

========================================
VBF Selection Summary (N=1000 events)
========================================
Jet multiplicity (2-4, pT>25):  433 (43.3%)
Leading jet pT (>60,>40):       588 (58.8%)  [was 347 (34.7%) with ATLAS cuts]
Delta phi_jj (<2):              685 (68.5%)
Opposite hemispheres:           643 (64.3%)
Delta eta_jj (>3.0):            520 (52.0%)  [was 394 (39.4%) with ATLAS cuts]
Dijet mass (>600 GeV):          339 (33.9%)  [was 263 (26.3%) with ATLAS cuts]
All jets pT sum (>140):         967 (96.7%)
Forward jet (2.38<|eta|<4.0):   525 (52.5%)
FSR compatibility (C_i,m_rel):  22 (2.2%)
----------------------------------------
FULL VBF SELECTION:             4 (0.4%)     [was 65 (6.5%) without FSR, 28 (2.8%) with ATLAS cuts]
========================================
```

## File Locations

- **Implementation**: `src/clustering_structs.h` (lines 223-359)
- **Test program**: `test_vbf_selection.cxx`
- **Dependencies**:
  - `src/clustering_includes.h` (TLorentzVector added)
  - `CMakeLists.txt` (updated to build test executable)

## Physics Background

The VBF H→invisible analysis searches for Higgs bosons produced via vector boson fusion that decay to invisible particles (e.g., dark matter). The signature features:

- **Two high-pT jets** in opposite detector hemispheres (VBF topology)
- **Large rapidity separation** between jets (Δη > 3.8)
- **High dijet mass** (m_jj > 800 GeV)
- **Large missing transverse energy** from invisible decay products

The implemented cuts select events with VBF-like kinematics, which can be used for:
- Studying vertex timing performance in VBF events
- Testing clustering algorithms on forward-jet enriched samples
- Comparing timing resolution in different event topologies

## Reference

**"Search for invisible Higgs-boson decays in events with vector-boson fusion signatures using 139 fb⁻¹ of proton-proton data recorded by the ATLAS experiment"**
ATLAS Collaboration
JHEP 08 (2022) 104
https://doi.org/10.1007/JHEP08(2022)104
