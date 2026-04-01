"""
hibeam.plotting.segmentation_plots — Segmentation study figures
=================================================================
Multi-panel publication figures for the detector segmentation study.

Consolidates:
  - segmentation/segmentation_overlay.py
  - segmentation/overlay.py
  - segmentation/nhits.py, nhits_pub.py, nhits_vs_nsections.py
  - segmentation/segment.py
  - segmentation/simple.py
"""

from __future__ import annotations

from pathlib import Path
from typing import Any, Optional

import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
import numpy as np

from hibeam.physics.fitting import fit_landau, FitResult, moyal_scaled
from hibeam.plotting.style import (
    SERIES_COLORS, add_logo, save_figure, get_colors,
)


# ═══════════════════════════════════════════════════════════════════════════════
# Landau overlay — all segmentations on one figure
# ═══════════════════════════════════════════════════════════════════════════════

def overlay_segmentations(
    records: list[dict],
    truncation: float = 0.70,
    n_bins: int = 40,
    edep_floor_mev: float = 0.15,
    title: str = "Segmentation overlay",
    output: Optional[str | Path] = None,
    cfg: dict | None = None,
) -> tuple[plt.Figure, list[dict]]:
    """Overlay fitted Landau curves for all segmentations.

    Parameters
    ----------
    records : list[dict]
        Output of ``seg_loader.load()``.
    truncation : float
        Fit truncation fraction.
    n_bins : int
        Bins for fit (coarser for stability).
    edep_floor_mev : float
        Lower Edep cut in MeV.
    title : str
        Figure title.

    Returns
    -------
    (fig, fit_results) — figure and list of fit result dicts per nsec.
    """
    colors = get_colors(cfg)
    cmap = plt.cm.viridis(np.linspace(0.1, 0.9, len(records)))

    fig, (ax_main, ax_params) = plt.subplots(
        1, 2, figsize=(16, 7), gridspec_kw={"width_ratios": [2, 1]})

    fit_results = []
    nsec_list, mpv_list, xi_list = [], [], []

    for idx, rec in enumerate(records):
        nsec = rec["nsec"]
        edep = rec["edep_per_event"]

        if len(edep) < 50:
            print(f"    nSec={nsec}: too few events ({len(edep)}), skipping.")
            continue

        # Convert GeV → MeV if needed
        if edep.max() < 1.0:
            edep = edep * 1000.0

        # Apply floor cut
        if edep_floor_mev > 0:
            edep = edep[edep > edep_floor_mev]

        if len(edep) < 50:
            continue

        try:
            result = fit_landau(edep, truncation=truncation, n_bins=n_bins)
        except Exception as e:
            print(f"    nSec={nsec}: fit failed ({e})")
            continue

        # Plot fitted curve
        x = np.linspace(result.bin_centers[0], result.bin_centers[-1], 500)
        y = moyal_scaled(x, result.loc, result.scale, result.norm)
        color = cmap[idx]
        ax_main.plot(x, y, color=color, lw=1.5,
                     label=f"nSec={nsec}  MPV={result.mpv:.3f}")

        nsec_list.append(nsec)
        mpv_list.append(result.mpv)
        xi_list.append(result.scale)
        fit_results.append({"nsec": nsec, "result": result})

    ax_main.set_xlabel(r"$\sum E_\mathrm{dep}$ [MeV]", fontsize=12)
    ax_main.set_ylabel("Events / bin", fontsize=11)
    ax_main.set_title(title, fontsize=11)
    ax_main.legend(fontsize=7, loc="upper right", ncol=2)
    add_logo(ax_main, cfg)

    # ── Parameter panel ──────────────────────────────────────────────────
    if nsec_list:
        ax_params.plot(nsec_list, mpv_list, "o-", color=colors[0] if colors else "C0",
                       label="MPV [MeV]")
        ax_params.set_xlabel("nSections", fontsize=12)
        ax_params.set_ylabel("MPV [MeV]", fontsize=11, color=colors[0] if colors else "C0")

        ax2 = ax_params.twinx()
        ax2.plot(nsec_list, xi_list, "s--", color=colors[1] if len(colors) > 1 else "C1",
                 label=r"$\xi$ [MeV]")
        ax2.set_ylabel(r"Width $\xi$ [MeV]", fontsize=11,
                        color=colors[1] if len(colors) > 1 else "C1")
        ax_params.set_title("Fit parameters vs segmentation", fontsize=11)

    fig.tight_layout()

    if output:
        formats = (cfg or {}).get("plotting", {}).get("formats", ["pdf", "png"])
        save_figure(fig, output, formats=formats)

    return fig, fit_results


# ═══════════════════════════════════════════════════════════════════════════════
# nHits vs nSections
# ═══════════════════════════════════════════════════════════════════════════════

def plot_nhits_vs_nsec(
    records: list[dict],
    title: str = "nHits vs segmentation",
    output: Optional[str | Path] = None,
    cfg: dict | None = None,
) -> plt.Figure:
    """Two-panel figure: mean/mode nHits and zero-hit fraction vs nSec.

    Parameters
    ----------
    records : list[dict]
        Output of ``seg_loader.load()``.
    """
    colors = get_colors(cfg)
    nsec = [r["nsec"] for r in records]
    mean_all = [r["mean_nhits_all"] for r in records]
    mean_nz  = [r["mean_nhits_nonzero"] for r in records]
    mode_nz  = [r["mode_nhits"] for r in records]
    zero_frac = [r["zero_frac"] for r in records]

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5.5))

    # Left panel: mean / mode nHits
    ax1.plot(nsec, mean_all, "o-", color=colors[0] if colors else "C0",
             label="Mean (all events)")
    ax1.plot(nsec, mean_nz, "s-", color=colors[1] if len(colors) > 1 else "C1",
             label="Mean (non-zero)")
    ax1.plot(nsec, mode_nz, "^--", color=colors[2] if len(colors) > 2 else "C2",
             label="Mode (non-zero)")
    ax1.set_xlabel("nSections", fontsize=12)
    ax1.set_ylabel("nHits per event", fontsize=11)
    ax1.legend(fontsize=9)
    ax1.set_title(title, fontsize=11)

    # Right panel: zero-hit fraction
    ax2.plot(nsec, [f * 100 for f in zero_frac], "D-",
             color=colors[3] if len(colors) > 3 else "C3")
    ax2.set_xlabel("nSections", fontsize=12)
    ax2.set_ylabel("Zero-hit events [%]", fontsize=11)
    ax2.set_title("Zero-hit fraction vs segmentation", fontsize=11)

    fig.tight_layout()

    if output:
        formats = (cfg or {}).get("plotting", {}).get("formats", ["pdf", "png"])
        save_figure(fig, output, formats=formats)

    return fig


# ═══════════════════════════════════════════════════════════════════════════════
# Per-segmentation Edep histogram (publication quality)
# ═══════════════════════════════════════════════════════════════════════════════

def plot_nhits_distribution(
    records: list[dict],
    title: str = "nHits distribution",
    output: Optional[str | Path] = None,
    cfg: dict | None = None,
) -> plt.Figure:
    """Stacked histogram of nHits for selected segmentations."""
    colors = get_colors(cfg)
    cmap = plt.cm.viridis(np.linspace(0.1, 0.9, len(records)))

    fig, ax = plt.subplots(figsize=(10, 6))

    for idx, rec in enumerate(records):
        nhits = rec["nhits"]
        nhits_nz = nhits[nhits > 0]
        if len(nhits_nz) == 0:
            continue
        bins = np.arange(0, min(nhits_nz.max() + 2, 100))
        ax.hist(nhits_nz, bins=bins, histtype="step", lw=1.5,
                color=cmap[idx], label=f"nSec={rec['nsec']}")

    ax.set_xlabel("nHits per event (non-zero)", fontsize=12)
    ax.set_ylabel("Events", fontsize=11)
    ax.set_title(title, fontsize=11)
    ax.legend(fontsize=7, ncol=2, loc="upper right")
    add_logo(ax, cfg)

    if output:
        formats = (cfg or {}).get("plotting", {}).get("formats", ["pdf", "png"])
        save_figure(fig, output, formats=formats)

    return fig
