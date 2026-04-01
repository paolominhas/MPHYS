"""
hibeam.plotting.overlays — Data vs simulation overlay figures
===============================================================
Peak-normalised or area-normalised overlay of experimental and
simulation dE/dx distributions with independent Landau fits.

Consolidates:
  - combined_analysis/dedx_data_comparison.py
  - combined_analysis/dedx_overlay_publication.py
  - combined_analysis/combination_try.py
  - combined_analysis/combining_dedx_simple.py
"""

from __future__ import annotations

from pathlib import Path
from typing import Optional

import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import moyal

from hibeam.physics.fitting import moyal_scaled, fit_landau
from hibeam.plotting.style import (
    SERIES_COLORS, add_logo, save_figure, get_colors,
)


def overlay_data_sim(
    data_values: np.ndarray,
    sim_values: np.ndarray,
    data_label: str = "Data",
    sim_label: str = "Simulation",
    sim_weights: np.ndarray | None = None,
    normalisation: str = "peak",
    truncation: float = 0.70,
    n_bins: int = 50,
    title: str = "",
    output: Optional[str | Path] = None,
    cfg: dict | None = None,
) -> plt.Figure:
    """Overlay two dE/dx distributions with independent Landau fits.

    Parameters
    ----------
    data_values, sim_values : np.ndarray
        Energy-loss values from experiment and simulation.
    data_label, sim_label : str
        Legend labels.
    sim_weights : np.ndarray or None
        Per-event weights for simulation (PrimaryWeight).
    normalisation : str
        ``"peak"`` — divide both by their histogram peak (dimensionless).
        ``"area"`` — normalise both to unit area.
    truncation : float
        Fit truncation fraction.
    n_bins : int
        Number of bins for the shared histogram.
    title : str
        Figure title.
    output : str or Path or None
        Save path (no extension — formats from config).
    cfg : dict or None
        Config for style overrides.

    Returns
    -------
    plt.Figure
    """
    colors = get_colors(cfg)
    c_data = "black"
    c_sim  = colors[0] if colors else "#4393c3"
    figsize = (cfg or {}).get("plotting", {}).get("figsize", [11, 9])

    if sim_weights is None:
        sim_weights = np.ones(len(sim_values))

    # ── Normalisation ────────────────────────────────────────────────────
    if normalisation == "peak":
        from hibeam.utils import histogram_peak
        data_peak = histogram_peak(data_values)
        sim_peak  = histogram_peak(sim_values, weights=sim_weights)
        data_norm = data_values / data_peak
        sim_norm  = sim_values / sim_peak
        x_label = r"$\mathrm{d}E/\mathrm{d}x$ [units of peak MPV]"
    else:
        data_norm = data_values
        sim_norm  = sim_values
        x_label = r"$\mathrm{d}E/\mathrm{d}x$"

    # ── Shared bin edges ─────────────────────────────────────────────────
    x_min = 0
    x_max = min(np.percentile(data_norm, 99.5),
                np.percentile(sim_norm, 99.5))
    edges   = np.linspace(x_min, x_max, n_bins + 1)
    centers = 0.5 * (edges[:-1] + edges[1:])
    bw      = edges[1] - edges[0]

    d_counts, _ = np.histogram(data_norm, bins=edges, density=True)
    s_counts, _ = np.histogram(sim_norm, bins=edges,
                                weights=sim_weights, density=True)
    d_err = np.sqrt(np.maximum(d_counts, 1e-9))
    s_err = np.sqrt(np.maximum(s_counts, 1e-9))

    # ── Independent fits on the bulk ─────────────────────────────────────
    cut = np.percentile(data_norm, truncation * 100)
    fm  = centers <= cut

    def _quick_fit(c, e, mask):
        p0 = [max(c[mask]), centers[mask][np.argmax(c[mask])],
              (cut - x_min) * 0.05]
        bounds = ([0, x_min, 1e-9], [np.inf, cut, cut])
        try:
            from scipy.optimize import curve_fit
            popt, _ = curve_fit(
                lambda x, a, l, s: a * moyal.pdf(x, l, s),
                centers[mask], c[mask], p0=p0, sigma=e[mask],
                bounds=bounds, maxfev=10000,
            )
            return popt
        except Exception:
            return p0

    d_popt = _quick_fit(d_counts, d_err, fm)
    s_popt = _quick_fit(s_counts, s_err, fm)
    x_smooth = np.linspace(x_min, x_max, 1000)

    # ── Figure ───────────────────────────────────────────────────────────
    fig, (ax_main, ax_res) = plt.subplots(
        2, 1, sharex=True, figsize=figsize,
        gridspec_kw={"height_ratios": [3, 1], "hspace": 0},
    )

    # Tail shading
    ax_main.axvspan(cut, x_max, alpha=0.08, color="gray",
                    label=f"Excluded tail ({100 - int(truncation*100)}%)")

    # Histograms
    try:
        import mplhep as hep
        hep.histplot(d_counts, bins=edges, yerr=d_err, ax=ax_main,
                     histtype="errorbar", color=c_data, label=data_label)
        hep.histplot(s_counts, bins=edges, yerr=s_err, ax=ax_main,
                     histtype="fill", color=c_sim, alpha=0.4,
                     label=sim_label)
    except ImportError:
        ax_main.bar(centers, d_counts, width=bw * 0.9, color=c_data,
                    alpha=0.6, label=data_label)
        ax_main.bar(centers, s_counts, width=bw * 0.7, color=c_sim,
                    alpha=0.4, label=sim_label)

    # Fit lines
    ax_main.plot(x_smooth, d_popt[0] * moyal.pdf(x_smooth, d_popt[1], d_popt[2]),
                 color="#E31A1C", lw=2,
                 label=f"{data_label} fit MPV={d_popt[1]:.3f}")
    ax_main.plot(x_smooth, s_popt[0] * moyal.pdf(x_smooth, s_popt[1], s_popt[2]),
                 color="darkorange", lw=2, ls="--",
                 label=f"{sim_label} fit MPV={s_popt[1]:.3f}")
    ax_main.axvline(cut, color="green", lw=1.5, ls=":")

    ax_main.set_yscale("log")
    ax_main.set_ylabel("Probability density")
    ax_main.set_ylim(bottom=1e-4)
    if title:
        ax_main.set_title(title, fontsize=11)
    ax_main.legend(loc="upper right", fontsize="small")
    add_logo(ax_main, cfg)

    # ── Residuals ────────────────────────────────────────────────────────
    d_fit_vals = d_popt[0] * moyal.pdf(centers, d_popt[1], d_popt[2])
    s_fit_vals = s_popt[0] * moyal.pdf(centers, s_popt[1], s_popt[2])
    d_res = np.where(fm, d_counts - d_fit_vals, 0)
    s_res = np.where(fm, s_counts - s_fit_vals, 0)

    ax_res.errorbar(centers[fm], d_res[fm], yerr=d_err[fm],
                    fmt="o", color=c_data, ms=4, label=f"{data_label} residuals")
    ax_res.errorbar(centers[fm], s_res[fm], yerr=s_err[fm],
                    fmt="s", color=c_sim, ms=4, alpha=0.7,
                    label=f"{sim_label} residuals")
    ax_res.axhline(0, color="gray", ls="--", lw=1)
    ax_res.axvline(cut, color="green", lw=1.5, ls=":")
    ax_res.set_xlabel(x_label)
    ax_res.set_ylabel("Data − Fit")
    ax_res.legend(loc="upper right", fontsize="small")

    plt.tight_layout()

    if output:
        formats = (cfg or {}).get("plotting", {}).get("formats", ["pdf", "png"])
        save_figure(fig, output, formats=formats)

    return fig
