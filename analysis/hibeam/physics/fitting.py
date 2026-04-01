"""
hibeam.physics.fitting — Landau (Moyal) fitting engine
========================================================
The ONE place where Landau fitting lives.  Previously reimplemented in:
  - segmentation/final.py       (fit_landau — most complete version)
  - segmentation/segment.py     (landau_approx + curve_fit)
  - segmentation/simple.py      (landau_approx + curve_fit)
  - combined_analysis/combining_dedx.py  (fit_landau — simplified)
  - simulation_analysis/sim_dedx_plots.py
  - experimental_analysis/dx-graph.py

Now there is exactly one implementation.  All scripts call this.

References
----------
[1] Landau, L.D. (1944). J. Phys. USSR 8, 201.
[2] Moyal, J.E. (1955). Phil. Mag. 46, 263.
[3] Bichsel, H. (1988). Rev. Mod. Phys. 60, 663.
[4] ALICE Collaboration (2010). JINST 5, P09002.
[5] PDG (2022). Statistics chapter — χ² GOF, pulls, p-values.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any

import numpy as np
from scipy import optimize, stats


# ═══════════════════════════════════════════════════════════════════════════════
# Fit models
# ═══════════════════════════════════════════════════════════════════════════════

def moyal_scaled(x: np.ndarray, loc: float, scale: float, norm: float) -> np.ndarray:
    """Moyal PDF × normalisation constant.

    The Moyal distribution is the standard analytic approximation to the
    Landau.  The mode (MPV) is at x = loc; the width parameter 'scale'
    corresponds to the Landau width ξ.  norm = N_events × bin_width.
    """
    return norm * stats.moyal.pdf(x, loc=loc, scale=scale)


def moyal_gauss(x: np.ndarray, loc: float, scale: float,
                sigma: float, norm: float) -> np.ndarray:
    """Moyal ⊗ Gaussian convolution (Landau × detector resolution).

    The Gaussian component absorbs detector resolution, pad response,
    and kinetic-energy spread.  Ref: Meroli et al., JINST 6 (2011) P06013.
    """
    from scipy.signal import fftconvolve

    dx = x[1] - x[0] if len(x) > 1 else 1.0
    moyal_part = stats.moyal.pdf(x, loc=loc, scale=scale)
    gauss_kernel = np.exp(-0.5 * ((x - x.mean()) / max(sigma, 1e-6))**2)
    gauss_kernel /= gauss_kernel.sum()
    convolved = fftconvolve(moyal_part, gauss_kernel, mode="same")
    return norm * convolved


_MODELS = {
    "moyal":       (moyal_scaled, 3),       # (function, n_params)
    "moyal_gauss": (moyal_gauss,  4),
}


# ═══════════════════════════════════════════════════════════════════════════════
# Fit result container
# ═══════════════════════════════════════════════════════════════════════════════

@dataclass
class FitResult:
    """Container for all fit outputs — used by plotting and summary functions."""

    # Fit parameters
    loc:        float = 0.0
    loc_err:    float = 0.0
    scale:      float = 0.0
    scale_err:  float = 0.0
    norm:       float = 0.0
    norm_err:   float = 0.0
    mpv:        float = 0.0      # mode of Moyal = loc
    popt:       np.ndarray = field(default_factory=lambda: np.array([]))
    pcov:       np.ndarray = field(default_factory=lambda: np.array([]))

    # Fit-region goodness-of-fit
    chi2:       float = 0.0
    ndf:        int   = 0
    chi2_red:   float = 0.0
    p_value:    float = 0.0

    # Tail-region deviation
    chi2_tail:  float = float("nan")
    ndf_tail:   int   = 0
    p_tail:     float = float("nan")

    # Arrays for plotting
    cut_value:   float = 0.0
    bin_edges:   np.ndarray = field(default_factory=lambda: np.array([]))
    bin_centers: np.ndarray = field(default_factory=lambda: np.array([]))
    counts:      np.ndarray = field(default_factory=lambda: np.array([]))
    fit_curve:   np.ndarray = field(default_factory=lambda: np.array([]))
    fit_mask:    np.ndarray = field(default_factory=lambda: np.array([]))
    bin_width:   float = 0.0
    truncation:  float = 0.70
    model_name:  str   = "moyal"

    def as_dict(self) -> dict[str, Any]:
        """Convert to plain dict (backward-compatible with old code)."""
        return {k: getattr(self, k) for k in self.__dataclass_fields__}


# ═══════════════════════════════════════════════════════════════════════════════
# Main fit function
# ═══════════════════════════════════════════════════════════════════════════════

def fit_landau(
    values: np.ndarray,
    truncation: float = 0.70,
    n_bins: int = 100,
    model: str = "moyal",
    min_bin_counts: int = 5,
    max_iterations: int = 20_000,
) -> FitResult:
    """Fit a Landau (Moyal) distribution to energy-loss values.

    Uses binned Pearson χ² minimisation on the lower ``truncation``
    fraction of the data — the standard TPC dE/dx technique [Ref 4, 5].

    Parameters
    ----------
    values : np.ndarray
        Per-event energy loss values.
    truncation : float
        Fraction of data (from the low end) included in the fit.
    n_bins : int
        Number of histogram bins spanning the fit region.
    model : str
        Fit model: ``"moyal"`` or ``"moyal_gauss"``.
    min_bin_counts : int
        Minimum counts per bin for χ² validity.
    max_iterations : int
        Maximum iterations for curve_fit.

    Returns
    -------
    FitResult
        All fit parameters, uncertainties, GOF statistics, and
        arrays for plotting.
    """
    if model not in _MODELS:
        raise ValueError(f"Unknown model '{model}'.  Options: {list(_MODELS.keys())}")

    fit_func, n_params = _MODELS[model]

    if len(values) < 50:
        raise ValueError(f"Only {len(values)} events — need at least 50.")

    cut_value = np.percentile(values, truncation * 100.0)

    # ── Build histogram ──────────────────────────────────────────────────
    hist_min = float(values.min())
    hist_max = float(np.percentile(values, 99.0))
    fit_range = cut_value - hist_min
    if fit_range <= 0:
        raise ValueError("Truncation cut equals data minimum — degenerate data.")

    bin_width_target = fit_range / n_bins
    n_total = max(int(np.ceil((hist_max - hist_min) / bin_width_target)), n_bins + 1)
    counts, bin_edges = np.histogram(values, bins=n_total,
                                     range=(hist_min, hist_max))
    bin_centers = 0.5 * (bin_edges[:-1] + bin_edges[1:])
    bin_width   = float(bin_edges[1] - bin_edges[0])

    fit_mask = bin_centers <= cut_value
    valid = fit_mask & (counts >= min_bin_counts)
    if valid.sum() < n_params + 1:
        raise RuntimeError(
            f"Only {valid.sum()} bins with ≥{min_bin_counts} counts. "
            "Reduce n_bins or supply more events."
        )

    x_fit = bin_centers[valid]
    y_fit = counts[valid].astype(float)

    # ── Initial guesses ──────────────────────────────────────────────────
    peak_idx = np.argmax(counts[fit_mask])
    loc0   = float(bin_centers[fit_mask][peak_idx])
    scale0 = float(fit_range / 8.0)
    norm0  = float(counts.sum()) * bin_width

    if model == "moyal":
        p0 = [loc0, scale0, norm0]
    else:
        p0 = [loc0, scale0, scale0 * 0.5, norm0]

    # ── Fit ──────────────────────────────────────────────────────────────
    try:
        popt, pcov = optimize.curve_fit(
            fit_func, x_fit, y_fit,
            p0=p0,
            sigma=np.sqrt(np.maximum(y_fit, 1.0)),
            absolute_sigma=True,
            maxfev=max_iterations,
        )
    except RuntimeError as exc:
        raise RuntimeError(
            f"Fit did not converge: {exc}\n"
            "Try adjusting n_bins or check that the data looks Landau-like."
        ) from exc

    perr = np.sqrt(np.diag(pcov))
    loc, scale, norm = popt[0], popt[1], popt[-1]
    loc_err, scale_err, norm_err = perr[0], perr[1], perr[-1]

    # ── Fit curve ────────────────────────────────────────────────────────
    fit_curve = fit_func(bin_centers, *popt)

    # ── χ² — fit region ──────────────────────────────────────────────────
    exp_fit = fit_curve[valid]
    chi2    = float(np.sum((y_fit - exp_fit)**2 / np.maximum(exp_fit, 1e-6)))
    ndf     = int(valid.sum()) - n_params
    chi2_red = chi2 / max(ndf, 1)
    p_value  = float(stats.chi2.sf(chi2, max(ndf, 1)))

    # ── χ² — tail region ─────────────────────────────────────────────────
    tail_mask = (~fit_mask) & (counts >= 1)
    if tail_mask.sum() > 0:
        obs_tail = counts[tail_mask].astype(float)
        exp_tail = np.maximum(fit_curve[tail_mask], 1e-6)
        chi2_tail = float(np.sum((obs_tail - exp_tail)**2 / exp_tail))
        ndf_tail  = int(tail_mask.sum())
        p_tail    = float(stats.chi2.sf(chi2_tail, ndf_tail))
    else:
        chi2_tail, ndf_tail, p_tail = float("nan"), 0, float("nan")

    return FitResult(
        loc=loc, loc_err=loc_err,
        scale=scale, scale_err=scale_err,
        norm=norm, norm_err=norm_err,
        mpv=loc,
        popt=popt, pcov=pcov,
        chi2=chi2, ndf=ndf, chi2_red=chi2_red, p_value=p_value,
        chi2_tail=chi2_tail, ndf_tail=ndf_tail, p_tail=p_tail,
        cut_value=cut_value,
        bin_edges=bin_edges, bin_centers=bin_centers,
        counts=counts, fit_curve=fit_curve, fit_mask=fit_mask,
        bin_width=bin_width, truncation=truncation,
        model_name=model,
    )


# ═══════════════════════════════════════════════════════════════════════════════
# Fit uncertainty band (for plotting)
# ═══════════════════════════════════════════════════════════════════════════════

def fit_band(
    x: np.ndarray,
    result: FitResult,
) -> np.ndarray:
    """1σ uncertainty envelope via numerical Jacobian error propagation.

    Returns the standard deviation of the fit curve at each x-point.
    """
    fit_func, _ = _MODELS[result.model_name]
    popt = result.popt
    pcov = result.pcov
    eps = 1e-6

    jac = np.zeros((len(x), len(popt)))
    for i, p in enumerate(popt):
        step = max(abs(p) * eps, eps)
        p_up = popt.copy(); p_up[i] += step
        p_dn = popt.copy(); p_dn[i] -= step
        jac[:, i] = (fit_func(x, *p_up) - fit_func(x, *p_dn)) / (2.0 * step)

    variance = np.einsum("ij,jk,ik->i", jac, pcov, jac)
    return np.sqrt(np.maximum(variance, 0.0))
