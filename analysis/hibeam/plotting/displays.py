"""
hibeam.plotting.displays — 3D event displays and pad-plane diagnostics
========================================================================
Consolidates:
  - experimental_analysis/tpc_3d_tracks.py
  - experimental_analysis/plot_tracks.py
  - experimental_analysis/3Dplot.py, plot3d.py, newtracks.py
  - experimental_analysis/2Dhisto.py
  - simulation_analysis/3Dplot.py
  - simulation_analysis/pid_deltae_e.py (PID plot)
"""

from __future__ import annotations

from pathlib import Path
from typing import Optional

import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import numpy as np

from hibeam.plotting.style import add_logo, save_figure, get_colors


# ═══════════════════════════════════════════════════════════════════════════════
# 3D track scatter display
# ═══════════════════════════════════════════════════════════════════════════════

def event_display_3d(
    tracks: list[dict],
    event_ids: list[int] | None = None,
    max_tracks: int = 100,
    cmap: str = "viridis",
    title: str = "3D TPC track display",
    output: Optional[str | Path] = None,
    cfg: dict | None = None,
) -> plt.Figure:
    """3D scatter plot of TPC tracks coloured by charge/energy.

    Parameters
    ----------
    tracks : list[dict]
        Track dicts from ``exp_loader.load()["tracks"]``, each with
        keys ``x``, ``y``, ``z``, ``charge``.
    event_ids : list[int] or None
        Specific event indices to display.  None → first ``max_tracks``.
    max_tracks : int
        Maximum tracks to display if event_ids is None.
    cmap : str
        Matplotlib colourmap.
    title : str
        Figure title.
    output : str or Path or None
        Save path.
    cfg : dict or None
        Config dict.
    """
    fig = plt.figure(figsize=(10, 8))
    ax = fig.add_subplot(111, projection="3d")

    if event_ids is not None:
        display_tracks = [t for t in tracks if t.get("event") in event_ids]
    else:
        display_tracks = tracks[:max_tracks]

    all_energy = []
    for t in display_tracks:
        all_energy.extend(t["charge"])

    if not all_energy:
        ax.text(0.5, 0.5, 0.5, "No tracks to display",
                transform=ax.transAxes, ha="center")
        return fig

    vmin, vmax = np.percentile(all_energy, [2, 98])
    norm = mcolors.Normalize(vmin=vmin, vmax=vmax)

    for t in display_tracks:
        x = np.array(t["x"])
        y = np.array(t["y"])
        z = np.array(t["z"])
        e = np.array(t["charge"])
        sort_idx = np.argsort(z)
        sc = ax.scatter(x[sort_idx], z[sort_idx], y[sort_idx],
                        c=e[sort_idx], cmap=cmap, norm=norm,
                        s=8, alpha=0.7)

    fig.colorbar(sc, ax=ax, label="Charge [ADC]", shrink=0.6, pad=0.1)
    ax.set_xlabel("X [mm]")
    ax.set_ylabel("Z [drift]")
    ax.set_zlabel("Y [mm]")
    ax.set_title(title, fontsize=11)

    if output:
        formats = (cfg or {}).get("plotting", {}).get("formats", ["pdf", "png"])
        save_figure(fig, output, formats=formats)
    return fig


# ═══════════════════════════════════════════════════════════════════════════════
# Pad-plane heatmap
# ═══════════════════════════════════════════════════════════════════════════════

def pad_plane_heatmap(
    row: np.ndarray,
    col: np.ndarray,
    signal: np.ndarray,
    title: str = "Pad-plane occupancy",
    output: Optional[str | Path] = None,
    cfg: dict | None = None,
) -> plt.Figure:
    """2D heatmap of mean signal per (row, column) pad.

    Parameters
    ----------
    row, col, signal : np.ndarray
        Flat arrays of pad row, column, and signal per hit.
    """
    import awkward as ak

    row_flat = ak.to_numpy(ak.flatten(row)) if hasattr(row, "type") else row
    col_flat = ak.to_numpy(ak.flatten(col)) if hasattr(col, "type") else col
    sig_flat = ak.to_numpy(ak.flatten(signal)) if hasattr(signal, "type") else signal

    n_rows = int(row_flat.max()) + 1
    n_cols = int(col_flat.max()) + 1

    sum_grid   = np.zeros((n_rows, n_cols))
    count_grid = np.zeros((n_rows, n_cols))

    np.add.at(sum_grid, (row_flat.astype(int), col_flat.astype(int)), sig_flat)
    np.add.at(count_grid, (row_flat.astype(int), col_flat.astype(int)), 1)

    mean_grid = np.where(count_grid > 0, sum_grid / count_grid, 0)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))

    im1 = ax1.imshow(mean_grid, aspect="auto", origin="lower", cmap="inferno")
    fig.colorbar(im1, ax=ax1, label="Mean signal [ADC]")
    ax1.set_xlabel("Column")
    ax1.set_ylabel("Row")
    ax1.set_title("Mean signal per pad")

    im2 = ax2.imshow(count_grid, aspect="auto", origin="lower", cmap="viridis")
    fig.colorbar(im2, ax=ax2, label="Hit count")
    ax2.set_xlabel("Column")
    ax2.set_ylabel("Row")
    ax2.set_title("Hit occupancy")

    fig.suptitle(title, fontsize=12)
    fig.tight_layout()

    if output:
        formats = (cfg or {}).get("plotting", {}).get("formats", ["pdf", "png"])
        save_figure(fig, output, formats=formats)
    return fig


# ═══════════════════════════════════════════════════════════════════════════════
# PID ΔE-E 2D histogram
# ═══════════════════════════════════════════════════════════════════════════════

def pid_plot(
    delta_e: np.ndarray,
    e_residual: np.ndarray,
    title: str = r"$\Delta E$–$E$ particle identification",
    output: Optional[str | Path] = None,
    cfg: dict | None = None,
) -> plt.Figure:
    """2D histogram of ΔE (TPC) vs E (calorimeter) for PID.

    Parameters
    ----------
    delta_e : np.ndarray
        TPC deposit [MeV] for selected events.
    e_residual : np.ndarray
        Calorimeter deposit [MeV] for selected events.
    """
    fig, ax = plt.subplots(figsize=(9, 7))

    x_max = np.percentile(e_residual, 99)
    y_max = np.percentile(delta_e, 99)

    h = ax.hist2d(
        e_residual, delta_e,
        bins=[80, 80],
        range=[[0, x_max], [0, y_max]],
        cmap="magma",
        norm=mcolors.LogNorm(vmin=1),
    )
    fig.colorbar(h[3], ax=ax, label="Events")

    ax.set_xlabel(r"$E_\mathrm{residual}$ (calorimeter) [MeV]", fontsize=13)
    ax.set_ylabel(r"$\Delta E$ (TPC) [MeV]", fontsize=13)
    ax.set_title(title, fontsize=12)
    add_logo(ax, cfg)

    if output:
        formats = (cfg or {}).get("plotting", {}).get("formats", ["pdf", "png"])
        save_figure(fig, output, formats=formats)
    return fig
