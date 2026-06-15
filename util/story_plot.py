#!/usr/bin/env python3
"""
story_plot.py

Five-page PDF, one panel per page, telling the story of the multi-Gaussian
timing structure in TRKPTZ cluster selection.

Page 1 — Baseline
          Standard TRKPTZ, all forward tracks.

Page 2 — Remove PU contamination
          HS-only (truth-matched) tracks, TRKPTZ selection.

Page 3 — Remove timing misassignments
          HS-only tracks, TRKPTZ selection, timing purity gate ≥ 95%.

Page 4 — Solve cluster selection
          HS-only tracks, oracle selection, timing purity gate ≥ 95%.

Page 5 — Split by N HS tracks in cluster  (≥ 9 shown)
          Single clean Gaussian — track-count variation was the root cause.

Pages 1–4 use a shared-mean mixture  Σ w_k · N(x; μ, σ_k).
Page 5 uses a single Gaussian + χ²/ndof.

Output: ../figs/story_plot.pdf
"""

import os, sys
import numpy as np
import pandas as pd
from scipy.stats import norm
from scipy.optimize import minimize
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.backends.backend_pdf as pdf_backend

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
CSV_PATH   = os.path.join(SCRIPT_DIR, "..", "hs_only_study.csv")
FIGS_DIR   = os.path.join(SCRIPT_DIR, "..", "figs")
PDF_PATH   = os.path.join(FIGS_DIR, "story_plot.pdf")

ZOOM = 120   # ps — shared x-range for all pages
BINS = 80
BW   = 2 * ZOOM / BINS
HS_TIMING_PURITY_CUT = 0.95

COMP_COLORS = ["#e377c2", "#ff7f0e", "#9467bd", "#8c564b"]

# ---------------------------------------------------------------------------
# Shared-mean Gaussian mixture
# ---------------------------------------------------------------------------
def fit_shared_mean_mixture(vals, n=2, n_restarts=6, seed=42,
                            sigma_inits=None, weight_reg=0.0):
    """
    Fit  Σ w_k · N(x; μ, σ_k)  with a shared mean.

    sigma_inits : optional list of length n — fixed σ starting point added
                  as an extra restart alongside the random ones.
    weight_reg  : Dirichlet-like regularisation strength. Adds
                  -weight_reg * Σ log(w_k) to the NLL, preventing any
                  component weight from collapsing to zero.  Values of
                  0.5-2.0 are typically enough to keep small components alive.
    """
    rng = np.random.default_rng(seed)
    N, mu0, s0 = len(vals), np.mean(vals), np.std(vals)

    def nll(params):
        mu     = params[0]
        sigmas = np.exp(params[1:1+n])
        if n == 1:
            weights = np.array([1.0])
        else:
            raw     = np.append(params[1+n:], 0.0)
            ew      = np.exp(raw - raw.max())
            weights = ew / ew.sum()
        pdf = sum(w * norm.pdf(vals, mu, s) for w, s in zip(weights, sigmas))
        reg = -weight_reg * np.sum(np.log(np.maximum(weights, 1e-300)))
        return -np.sum(np.log(np.maximum(pdf, 1e-300))) + reg

    bounds = [(-30, 30)] + [(-2, 5)] * n + [(-8, 8)] * max(n - 1, 0)
    best_nll, best_p = np.inf, None

    # Build list of starting points: random restarts + optional pinned init
    starts = []
    for _ in range(n_restarts):
        log_s = np.log(s0 * np.sort(rng.uniform(0.3, 2.0, n)))
        x0    = [mu0 + rng.uniform(-1, 1)] + list(log_s)
        if n > 1:
            x0 += list(rng.uniform(-1, 1, n - 1))
        starts.append(x0)
    if sigma_inits is not None:
        log_s_pinned = np.log(np.sort(sigma_inits))
        x0_pinned    = [mu0] + list(log_s_pinned)
        if n > 1:
            x0_pinned += [0.0] * (n - 1)
        starts.append(x0_pinned)

    for x0 in starts:
        res = minimize(nll, x0, method="L-BFGS-B", bounds=bounds,
                       options={"maxiter": 2000, "ftol": 1e-14})
        if res.fun < best_nll:
            best_nll, best_p = res.fun, res.x

    mu_fit     = best_p[0]
    sigmas_fit = np.exp(best_p[1:1+n])
    if n == 1:
        w_fit = np.array([1.0])
    else:
        raw   = np.append(best_p[1+n:], 0.0)
        ew    = np.exp(raw - raw.max())
        w_fit = ew / ew.sum()
    order = np.argsort(sigmas_fit)
    return mu_fit, sigmas_fit[order], w_fit[order]


def chi2_single_gauss(vals):
    mu, sigma = np.mean(vals), np.std(vals)
    edges = np.linspace(-ZOOM, ZOOM, BINS + 1)
    bw    = edges[1] - edges[0]
    obs, _ = np.histogram(vals, bins=edges)
    bc     = 0.5 * (edges[:-1] + edges[1:])
    exp    = norm.pdf(bc, mu, sigma) * len(vals) * bw
    mask   = exp >= 5
    if mask.sum() < 3:
        return mu, sigma, np.nan, 0
    chi2  = np.sum((obs[mask] - exp[mask])**2 / exp[mask])
    ndof  = int(mask.sum()) - 2
    return mu, sigma, chi2 / ndof if ndof > 0 else np.nan, ndof


# ---------------------------------------------------------------------------
# Page drawing
# ---------------------------------------------------------------------------
def compute_fit(vals, n_components, force_single, sigma_inits=None, weight_reg=0.0):
    """Return fit parameters without drawing anything."""
    if force_single:
        mu, sigma, chi2_ndf, ndof = chi2_single_gauss(vals)
        return ("single", mu, sigma, chi2_ndf, ndof)
    else:
        mu, sigmas, weights = fit_shared_mean_mixture(
            vals, n=n_components, sigma_inits=sigma_inits, weight_reg=weight_reg)
        return ("mixture", mu, sigmas, weights)


def draw_page(pdf, vals, color, page_title, subtitle,
              fit_result, ymax_linear, ymin_log=1e-5):
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))

    edges   = np.linspace(-ZOOM, ZOOM, BINS + 1)
    counts, _ = np.histogram(vals, bins=edges)
    density = counts / (len(vals) * BW)
    x       = np.linspace(-ZOOM, ZOOM, 1200)

    # Build fit curve(s)
    if fit_result[0] == "single":
        _, mu, sigma, chi2_ndf, ndof = fit_result
        fit_label = (f"Single Gaussian\n"
                     f"μ = {mu:+.1f} ps,  σ = {sigma:.1f} ps\n"
                     f"χ²/ndof = {chi2_ndf:.2f}  ({ndof} dof)")
        fit_total = norm.pdf(x, mu, sigma)
        fit_comps = []
    else:
        _, mu, sigmas, weights = fit_result
        fit_label = f"Shared-mean mixture\nμ = {mu:+.1f} ps"
        fit_total = sum(w * norm.pdf(x, mu, s) for w, s in zip(weights, sigmas))
        fit_comps = list(zip(sigmas, weights))

    for ax_idx, ax in enumerate(axes):
        log = (ax_idx == 1)

        ax.bar(edges[:-1], density, width=BW, align="edge",
               color=color, alpha=0.30, edgecolor="none")
        ax.step(np.append(edges[:-1], edges[-1]),
                np.append(density, 0),
                color=color, lw=1.8, where="post",
                label=f"N = {len(vals):,}")

        if fit_result[0] == "mixture":
            for k, (sig, w) in enumerate(fit_comps):
                y = w * norm.pdf(x, mu, sig)
                ax.plot(x, y, "--", color=COMP_COLORS[k % len(COMP_COLORS)],
                        lw=1.6, alpha=0.9,
                        label=f"σ = {sig:.1f} ps,  w = {w:.2f}")
        ax.plot(x, fit_total, color=color, lw=2.8, label=fit_label)

        ax.set_xlim(-ZOOM, ZOOM)
        ax.set_xlabel(r"$\Delta t$  [ps]", fontsize=11)
        ax.tick_params(labelsize=10)
        ax.legend(frameon=False, fontsize=9,
                  loc="upper right" if not log else "upper right")

        if log:
            ax.set_yscale("log")
            ax.set_ylim(ymin_log, ymax_linear * 3)
            ax.set_ylabel("Density  [ps⁻¹]  (log)", fontsize=11)
            ax.set_title("Log scale", fontsize=11)
        else:
            ax.set_ylim(0, ymax_linear)
            ax.set_ylabel("Density  [ps⁻¹]", fontsize=11)
            ax.set_title("Linear scale", fontsize=11)

    fig.suptitle(page_title, fontsize=14, fontweight="bold", y=1.01)
    fig.text(0.5, -0.04, subtitle, ha="center", va="top",
             fontsize=9, color="#444444", linespacing=1.5,
             transform=fig.transFigure)

    fig.tight_layout()
    pdf.savefig(fig, bbox_inches="tight")
    plt.close(fig)


def main():
    if not os.path.exists(CSV_PATH):
        print(f"ERROR: {CSV_PATH} not found.\nRun: cd build && ./hs_only_study",
              file=sys.stderr)
        sys.exit(1)

    df          = pd.read_csv(CSV_PATH)
    df_g        = df[df["hs_timing_purity"] >= HS_TIMING_PURITY_CUT]
    df_oracle_g = df_g[df_g["oracle_found"] == 1]
    df_hs9      = df_oracle_g[df_oracle_g["n_hs_in_cluster_oracle"] >= 9]

    print(f"Page 1 — Baseline:                {len(df):,}")
    print(f"Page 2 — Remove PU:               {len(df):,}")
    print(f"Page 3 — Remove timing misassign: {len(df_g):,}")
    print(f"Page 4 — Oracle selection:        {len(df_oracle_g):,}")
    print(f"Page 5 — n_hs ≥ 9:                {len(df_hs9):,}")

    pages = [
        dict(
            vals    = df["dt_std"].values,
            color   = "#1f77b4",
            title   = "Step 1 — Baseline",
            subtitle= (
                "Standard TRKPTZ selection, all forward tracks.\n"
                "Three components: narrow core (σ ≈ 11 ps), broad core (σ ≈ 25 ps),\n"
                "and a wide background tail."
            ),
            n_comp  = 3,
            single  = False,
        ),
        dict(
            vals        = df["dt_hs"].values,
            color       = "#2ca02c",
            n_comp      = 3,
            sigma_inits = [8.0, 14.0, 150.0],
            weight_reg  = 1.5,
            title       = "Step 2 — Remove PU contamination",
            subtitle= (
                "Cluster only truth-matched HS tracks; apply TRKPTZ selection.\n"
                "Overall σ: 53 ps → 16 ps.  Double-Gaussian persists.\n"
                "A broad tail (~0.4%) from HGTD timing misassignment on HS tracks remains."
            ),
            single  = False,
        ),
        dict(
            vals    = df_g["dt_hs"].values,
            color   = "#17becf",
            n_comp  = 2,
            title   = "Step 3 — Remove timing misassignments",
            subtitle= (
                "HS-only clustering + TRKPTZ, require HS timing purity ≥ 95%\n"
                "(≥ 95% of HS track pT within 3σ of per-particle truth time).\n"
                "σ: 16 ps → 11 ps.  Structure further reduced; 47% of events pass."
            ),
            single  = False,
        ),
        dict(
            vals    = df_oracle_g["dt_oracle"].values,
            color   = "#d62728",
            n_comp  = 2,
            title   = "Step 4 — Solve cluster selection",
            subtitle= (
                "HS-only clustering, oracle selection\n"
                "(best HS cluster within ±60 ps of truth), timing purity gate ≥ 95%.\n"
                "Minimal change from step 3 — cluster selection is not the remaining cause."
            ),
            single  = False,
        ),
        dict(
            vals    = df_hs9["dt_oracle"].values,
            color   = "#9467bd",
            n_comp  = 1,
            title   = "Step 5 — Split by N HS tracks in cluster  (≥ 9 shown)",
            subtitle= (
                "HS-only oracle selection, timing gate ≥ 95%,\n"
                "clusters containing ≥ 9 HS tracks.\n"
                "Single Gaussian — track-count variation σ ∝ 1/√n was the root cause."
            ),
            single  = True,
        ),
    ]

    os.makedirs(FIGS_DIR, exist_ok=True)

    # ── Pass 1: compute all fits and find global peak density ────────────────
    print("Computing fits…")
    for p in pages:
        p["fit"] = compute_fit(p["vals"], n_components=p["n_comp"],
                               force_single=p["single"],
                               sigma_inits=p.get("sigma_inits", None),
                               weight_reg=p.get("weight_reg", 0.0))
        edges  = np.linspace(-ZOOM, ZOOM, BINS + 1)
        counts, _ = np.histogram(p["vals"], bins=edges)
        p["peak"] = (counts / (len(p["vals"]) * BW)).max()
        print(f"  {p['title'][:40]:40s}  peak density = {p['peak']:.5f}")

    ymax = max(p["peak"] for p in pages) * 1.18

    # ── Pass 2: draw individual pages + final overlay page ───────────────────
    with pdf_backend.PdfPages(PDF_PATH) as pdf:
        for p in pages:
            draw_page(
                pdf,
                vals        = p["vals"],
                color       = p["color"],
                page_title  = p["title"],
                subtitle    = p["subtitle"],
                fit_result  = p["fit"],
                ymax_linear = ymax,
                ymin_log    = 1e-5,
            )

        # ── Page 6: all steps overlaid, linear + log ─────────────────────────
        fig, axes = plt.subplots(1, 2, figsize=(14, 6))
        edges = np.linspace(-ZOOM, ZOOM, BINS + 1)

        for p in pages:
            counts, _ = np.histogram(p["vals"], bins=edges)
            density   = counts / (len(p["vals"]) * BW)
            for ax in axes:
                ax.step(np.append(edges[:-1], edges[-1]),
                        np.append(density, 0),
                        color=p["color"], lw=1.8, where="post",
                        label=p["title"].replace(" — ", "\n", 1))
                ax.fill_between(edges[:-1], 0, density,
                                step="post", color=p["color"], alpha=0.12)

        axes[0].set_ylim(0, ymax)
        axes[0].set_ylabel("Density  [ps⁻¹]", fontsize=11)
        axes[0].set_title("Linear scale", fontsize=11)

        axes[1].set_yscale("log")
        axes[1].set_ylim(1e-5, ymax * 3)
        axes[1].set_ylabel("Density  [ps⁻¹]  (log)", fontsize=11)
        axes[1].set_title("Log scale", fontsize=11)

        for ax in axes:
            ax.set_xlim(-ZOOM, ZOOM)
            ax.set_xlabel(r"$\Delta t$  [ps]", fontsize=11)
            ax.legend(frameon=False, fontsize=8, loc="upper right")
            ax.tick_params(labelsize=10)

        fig.suptitle("All steps overlaid", fontsize=14, fontweight="bold", y=1.01)
        fig.tight_layout()
        pdf.savefig(fig, bbox_inches="tight")
        plt.close(fig)

    print(f"\nSaved → {PDF_PATH}")


if __name__ == "__main__":
    main()
