"""
hibeam.plotting.histograms — Histogram and pull-panel figures
===============================================================
Publication-quality 1D histograms with Landau fit overlays, pull
residual panels, and annotation boxes showing fit statistics.

Consolidates plotting from:
  - segmentation/final.py    (_draw_panels, plot_dedx)
  - simulation_analysis/tpc_dedx_pub.py
  - simulation_analysis/sim_dedx_plots.py
  - experimental_analysis/npointsenergy_graph.py
"""

from __future__ import annotations

from pathlib import Path
from typing import Any, Optional

import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
import numpy as np

from hibeam.physics.fitting import FitResult, fit_band, moyal_scaled
from hibeam.plotting.style import (
    FIT_COLOR, TAIL_COLOR, add_logo, get_hist_colors, save_figure,
)


# ═══════════════════════════════════════════════════════════════════════════════
# Core draw routine
# ═══════════════════════════════════════════════════════════════════════════════

def draw_dedx_panels(
    ax_top: plt.Axes,
    ax_pull: plt.Axes,
    result: FitResult,
    x_label: str = r"$\sum E_\mathrm{dep}\ [\mathrm{MeV}]$",
    x_max: Optional[float] = None,
    log_y: bool = False,
    title: Optional[str] = None,
    cfg: dict | None = None,
) -> None:
    """Draw the main histogram + fit on ax_top and pulls on ax_pull.

    This is the core drawing function.  It does NOT create figures or
    save files — call :func:`plot_dedx` for the full pipeline.

    Parameters
    ----------
    ax_top, ax_pull : Axes
        Pre-created axes for the distribution and pull panels.
    result : FitResult
        Output of ``fitting.fit_landau()``.
    x_label : str
        Physical quantity label for the x-axis.
    x_max : float or None
        Hard x-axis upper limit.  None → use 95th percentile of bins.
    log_y : bool
        Logarithmic y-scale on the top panel.
    title : str or None
        Custom title.  None → default HIBEAM title.
    cfg : dict or None
        Config dict for color/style overrides.
    """
    bc  = result.bin_centers
    cnt = result.counts.astype(float)
    fc  = result.fit_curve
    fm  = result.fit_mask
    bw  = result.bin_width

    c_fit, c_tail = get_hist_colors(cfg)
    h_cfg = (cfg or {}).get("plotting", {}).get("histogram", {})
    fit_alpha  = h_cfg.get("fit_alpha", 0.55)
    tail_alpha = h_cfg.get("tail_alpha", 0.35)
    tail_hatch = h_cfg.get("tail_hatch", "////")

    if x_max is None:
        x_max = float(np.percentile(bc, 95) + bw / 2.0)

    trunc_pct = int(result.truncation * 100)

    # ── Uncertainties ────────────────────────────────────────────────────
    cnt_err  = np.sqrt(np.maximum(cnt, 1.0))
    denom    = np.sqrt(np.maximum(fc, 1.0))
    pulls    = (cnt - fc) / denom
    pull_err = cnt_err / denom

    # ── Step error band helper ───────────────────────────────────────────
    def _step_band(mask, color, alpha=0.25):
        x_s = np.repeat(bc[mask], 2)
        half = bw / 2.0
        x_edges = np.empty_like(x_s)
        x_edges[0::2] = bc[mask] - half
        x_edges[1::2] = bc[mask] + half
        y_lo = np.repeat(np.maximum(cnt[mask] - cnt_err[mask], 0.0), 2)
        y_hi = np.repeat(cnt[mask] + cnt_err[mask], 2)
        ax_top.fill_between(x_edges, y_lo, y_hi,
                            step=None, color=color, alpha=alpha, linewidth=0)

    # ── Top panel: fitted region ─────────────────────────────────────────
    ax_top.bar(bc[fm], cnt[fm], width=bw * 0.96, color=c_fit, alpha=fit_alpha,
               label=f"Data — fitted (bottom {trunc_pct}%)")
    _step_band(fm, c_fit, alpha=0.35)

    # ── Top panel: excluded tail ─────────────────────────────────────────
    ax_top.bar(bc[~fm], cnt[~fm], width=bw * 0.96,
               color=c_tail, alpha=tail_alpha, hatch=tail_hatch,
               edgecolor=c_tail,
               label=f"Excluded tail (top {100 - trunc_pct}%)")
    _step_band(~fm, c_tail, alpha=0.20)

    # ── Fit curve + uncertainty band ─────────────────────────────────────
    x_dense = np.linspace(float(bc[0]), float(x_max), 3000)
    y_dense = moyal_scaled(x_dense, result.loc, result.scale, result.norm)
    y_sigma = fit_band(x_dense, result)

    ax_top.plot(x_dense, y_dense, color="black", lw=2.0, zorder=5,
                label="Landau (Moyal) fit")
    ax_top.fill_between(x_dense, y_dense - y_sigma, y_dense + y_sigma,
                        color="black", alpha=0.15, zorder=4,
                        label=r"Fit $\pm1\sigma$")

    # ── Truncation line ──────────────────────────────────────────────────
    ax_top.axvline(result.cut_value, color="dimgray", lw=1.4, ls="--",
                   alpha=0.8, label=f"Truncation ({trunc_pct}th pctl)")

    # ── Annotation box ───────────────────────────────────────────────────
    chi2_tail_red = result.chi2_tail / max(result.ndf_tail, 1)
    ann = (
        rf"MPV $= {result.mpv:.4g}$ MeV" "\n"
        rf"$\xi = {result.scale:.4g} \pm {result.scale_err:.2g}$ MeV" "\n\n"
        rf"Fit: $\chi^2/\nu = {result.chi2_red:.2f}$,"
        rf" $p = {result.p_value:.3f}$" "\n\n"
        rf"Tail: $\chi^2/\nu = {chi2_tail_red:.2e}$"
    )
    ax_top.text(0.975, 0.965, ann, transform=ax_top.transAxes,
                va="top", ha="right", fontsize=9,
                bbox=dict(boxstyle="round,pad=0.5", facecolor="white",
                          edgecolor="#aaaaaa", alpha=0.93))

    y_label = "Events / bin" + (" (log)" if log_y else "")
    ax_top.set_ylabel(y_label, fontsize=12)
    if title:
        ax_top.set_title(title, fontsize=11)
    ax_top.legend(fontsize=8.5, loc="upper right",
                  framealpha=0.92, edgecolor="#aaaaaa")
    ax_top.set_xlim(float(bc[0]) - bw, x_max)
    if log_y:
        ax_top.set_yscale("log")
        ax_top.set_ylim(bottom=0.5)
    plt.setp(ax_top.get_xticklabels(), visible=False)

    # ── Pull panel ───────────────────────────────────────────────────────
    ax_pull.bar(bc[fm], pulls[fm], width=bw * 0.96, color=c_fit, alpha=fit_alpha)
    ax_pull.bar(bc[~fm], pulls[~fm], width=bw * 0.96,
                color=c_tail, alpha=tail_alpha, hatch=tail_hatch,
                edgecolor=c_tail)

    ax_pull.axhline(0, color="black", lw=1.0)
    ax_pull.axhline(+2, color="#888888", lw=1.0, ls="--")
    ax_pull.axhline(-2, color="#888888", lw=1.0, ls="--")
    ax_pull.axvline(result.cut_value, color="dimgray", lw=1.4, ls="--", alpha=0.8)

    pull_range = (cfg or {}).get("plotting", {}).get("pulls", {}).get("range", [-6.5, 6.5])
    ax_pull.set_ylabel(r"$(N-F)/\!\sqrt{F}$", fontsize=11)
    ax_pull.set_xlabel(x_label, fontsize=12)
    ax_pull.set_ylim(*pull_range)
    ax_pull.set_xlim(float(bc[0]) - bw, x_max)


# ═══════════════════════════════════════════════════════════════════════════════
# Full figure pipeline
# ═══════════════════════════════════════════════════════════════════════════════

def plot_dedx(
    result: FitResult,
    use_dedx: bool = False,
    output: Optional[str | Path] = None,
    title: Optional[str] = None,
    x_max: Optional[float] = None,
    cfg: dict | None = None,
) -> tuple[plt.Figure, plt.Figure]:
    """Produce linear + log figures of the energy-loss distribution.

    Parameters
    ----------
    result : FitResult
        Output of ``fitting.fit_landau()``.
    use_dedx : bool
        Label as dE/dx [MeV/mm] (True) or ΣEdep [MeV] (False).
    output : str or Path or None
        Base path.  Produces ``<stem>_linear.pdf`` and ``<stem>_log.pdf``.
    title : str or None
        Custom title.
    x_max : float or None
        Display x-limit.
    cfg : dict or None
        Config for style overrides.

    Returns
    -------
    (fig_linear, fig_log)
    """
    x_label = (
        r"$\mathrm{d}E/\mathrm{d}x\ [\mathrm{MeV}\,\mathrm{mm}^{-1}]$"
        if use_dedx
        else r"$\sum E_\mathrm{dep}\ [\mathrm{MeV}]$"
    )

    hr = (cfg or {}).get("plotting", {}).get("pulls", {}).get(
        "height_ratio", [3, 1])
    figsize = (cfg or {}).get("plotting", {}).get("figsize", [11, 7])
    formats = (cfg or {}).get("plotting", {}).get("formats", ["pdf", "png"])

    figures = {}
    for log_y, suffix in [(False, "linear"), (True, "log")]:
        fig = plt.figure(figsize=figsize)
        gs = gridspec.GridSpec(2, 1, height_ratios=hr, hspace=0.04,
                               left=0.10, right=0.97, top=0.93, bottom=0.10)
        ax_top = fig.add_subplot(gs[0])
        ax_pull = fig.add_subplot(gs[1], sharex=ax_top)

        draw_dedx_panels(ax_top, ax_pull, result,
                         x_label=x_label, x_max=x_max,
                         log_y=log_y, title=title, cfg=cfg)
        add_logo(ax_top, cfg)
        figures[suffix] = fig

    if output:
        output = Path(output)
        for suffix, fig in figures.items():
            save_figure(fig, output.parent / f"{output.stem}_{suffix}",
                        formats=formats, close=True)
    else:
        plt.show()

    return figures.get("linear"), figures.get("log")


# ═══════════════════════════════════════════════════════════════════════════════
# Simple histogram (no fit)
# ═══════════════════════════════════════════════════════════════════════════════

def plot_simple_histogram(
    values: np.ndarray,
    bins: int = 100,
    x_label: str = "",
    y_label: str = "Events / bin",
    title: str = "",
    output: Optional[str | Path] = None,
    cfg: dict | None = None,
    **hist_kwargs: Any,
) -> plt.Figure:
    """Quick histogram without fitting — for diagnostics and exploration."""
    figsize = (cfg or {}).get("plotting", {}).get("figsize", [10, 6])
    fig, ax = plt.subplots(figsize=figsize)

    color = get_hist_colors(cfg)[0]
    ax.hist(values, bins=bins, color=color, alpha=0.6,
            edgecolor="none", **hist_kwargs)
    ax.set_xlabel(x_label, fontsize=12)
    ax.set_ylabel(y_label, fontsize=11)
    if title:
        ax.set_title(title, fontsize=11)
    add_logo(ax, cfg)

    if output:
        formats = (cfg or {}).get("plotting", {}).get("formats", ["pdf", "png"])
        save_figure(fig, output, formats=formats)
    return fig
