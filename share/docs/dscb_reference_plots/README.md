# DSCB Reference Plots

Snapshot of the inclusive timing-residual fits produced with a **Gaussian + symmetric DSCB** model during the fit-method exploration on branch `claude/vibrant-antonelli-1a8701`.

## Model

Two-component sum, mean fixed at 0, full ±1000 ps fit range, bin-integration enabled (`"RQI"`):

- **Narrow Gaussian core**: parameters `Norm_g`, `Sigma_g ∈ [1, 30] ps`
- **Symmetric DSCB wider component**: parameters `Norm_d`, `Sigma_d ∈ [10, 500] ps`, `Alpha ∈ [0.1, 10]`, `N ∈ [1.01, 100]`

Defined in [src/plotting_utilities.h](../../../src/plotting_utilities.h) as `gausDSCBFunc` and `createGausDSCBFit`. These functions remain in source for reference but are **not** wired into the active fit pipeline — the production fits use the legacy double-Gaussian (`createDblFit`).

## Headline results (TRKPTZ, total stack)

| | σ_g | σ_d | α | n | χ²/ndf |
|---|---|---|---|---|---|
| Double Gaussian (production) | 13.23 | — | — | — | ≈4.81 |
| Gauss+DSCB (these plots)     | 10.22 | 24.37 | 2.07 | 1.06 | **3.56** |

The narrow σ_g comes out roughly 25% tighter and χ²/ndf improves by ~25% on TRKPTZ specifically. On HGTD-driven scores the gain is smaller. The `n` parameter consistently rails near its lower bound (1.01) — the data's deep tail is flatter than any finite-n DSCB allows, a known limitation accepted in lieu of adding a flat-pedestal component (which we tried; it didn't help TRKPTZ).

## Plot inventory

| File | Description |
|---|---|
| `inclusivereso_logscale.pdf` | All scores, log-y, full residual + per-purity-category pages |
| `inclusivereso_linscale.pdf` | Same on linear y |
| `inclusivereso_lowtrack_logscale.pdf` | Same again restricted to low-HS-track events |
| `inclusivereso_lowtrack_linscale.pdf` | Low-track, linear y |

Each PDF has the stacked total page plus three per-category (signal/mixed/background) pages per score.

## To re-enable the DSCB fit

Swap the three `createDblFit(...)` calls in `inclusivePlot` and `PlotObj::plotPostProcessing` for `createGausDSCBFit(...)` and use `FitParams::fillGausDSCB` instead of `fillEach`/`fillCoreGaus`. The slice-label code in `plotLogic` will also need a 6-parameter branch added.
