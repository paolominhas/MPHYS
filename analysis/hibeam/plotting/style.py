"""
hibeam.plotting.style — Central plot style management
=======================================================
ONE place for all matplotlib configuration.  Previously every file
set its own rcParams, colours, and fonts.

Usage
-----
    from hibeam.plotting import style
    style.apply(cfg)                   # from config.yaml
    style.apply_preset("ess")          # direct preset
    style.apply_preset("cms")          # mplhep CMS style

Colour constants are importable directly:
    from hibeam.plotting.style import ESS_CYAN, SERIES_COLORS
"""

from __future__ import annotations

import os
from pathlib import Path
from typing import Any

import matplotlib
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
from matplotlib.offsetbox import OffsetImage, AnnotationBbox

# ═══════════════════════════════════════════════════════════════════════════════
# ESS brand palette
# ═══════════════════════════════════════════════════════════════════════════════

ESS_NAVY   = "#00263A"
ESS_BLUE   = "#005B8E"
ESS_CYAN   = "#00A9CE"
ESS_SILVER = "#A8B2B8"
ESS_WHITE  = "#FFFFFF"
ESS_ORANGE = "#E5720F"
ESS_GREEN  = "#4CAF82"
ESS_PURPLE = "#C45AB3"
ESS_GOLD   = "#FFD166"

# Accessible series colours — distinct in print and greyscale
SERIES_COLORS = [ESS_CYAN, ESS_ORANGE, ESS_GREEN, ESS_PURPLE, ESS_GOLD]

# Histogram region colours
FIT_COLOR  = "#2166ac"   # strong blue — fitted region
TAIL_COLOR = "#d73027"   # vermilion — excluded tail


# ═══════════════════════════════════════════════════════════════════════════════
# Style presets
# ═══════════════════════════════════════════════════════════════════════════════

def _ess_params(cfg: dict | None = None) -> dict:
    """ESS / HIBEAM publication style."""
    p = cfg.get("plotting", {}) if cfg else {}
    return {
        "font.family":          p.get("font_family", "serif"),
        "font.serif":           p.get("font_serif",
                                      ["DejaVu Serif", "Georgia", "Times New Roman"]),
        "mathtext.fontset":     p.get("mathtext", "dejavuserif"),
        "font.size":            p.get("font_size", 10),
        "axes.titlesize":       p.get("title_size", 11),
        "axes.labelsize":       p.get("label_size", 10),
        "xtick.labelsize":      p.get("tick_size", 8.5),
        "ytick.labelsize":      p.get("tick_size", 8.5),
        "legend.fontsize":      p.get("legend_size", 8),
        "lines.linewidth":      1.6,
        "lines.markersize":     4.5,
        "axes.linewidth":       0.8,
        "xtick.major.width":    0.8,
        "ytick.major.width":    0.8,
        "xtick.minor.width":    0.5,
        "ytick.minor.width":    0.5,
        "xtick.direction":      "in",
        "ytick.direction":      "in",
        "xtick.top":            True,
        "ytick.right":          True,
        "xtick.minor.visible":  True,
        "ytick.minor.visible":  True,
        "figure.facecolor":     "white",
        "figure.figsize":       p.get("figsize", [11, 7]),
        "savefig.dpi":          p.get("dpi", 300),
        "savefig.bbox":         "tight",
    }


def _mplhep_style(name: str) -> None:
    """Apply an mplhep experiment style (CMS, ATLAS, ALICE)."""
    try:
        import mplhep as hep
        styles = {
            "cms":   hep.style.CMS,
            "atlas": "ATLAS",
            "alice": hep.style.ALICE if hasattr(hep.style, "ALICE") else hep.style.CMS,
        }
        style_obj = styles.get(name, hep.style.CMS)
        if isinstance(style_obj, str):
            hep.style.use(style_obj)
        else:
            plt.style.use(style_obj)
    except ImportError:
        pass  # fall back to ESS style


# ═══════════════════════════════════════════════════════════════════════════════
# Public API
# ═══════════════════════════════════════════════════════════════════════════════

def apply(cfg: dict | None = None) -> None:
    """Apply plot style from a config dict (or config.yaml DotDict).

    Reads ``cfg["plotting"]["style"]`` to choose the base preset,
    then overlays any custom rcParams from the config.
    """
    if cfg is None:
        apply_preset("ess")
        return

    p = cfg.get("plotting", {}) if isinstance(cfg, dict) else cfg
    style_name = p.get("style", "ess") if isinstance(p, dict) else "ess"
    apply_preset(style_name, cfg)


def apply_preset(name: str = "ess", cfg: dict | None = None) -> None:
    """Apply a named style preset.

    Parameters
    ----------
    name : str
        One of: ``"ess"``, ``"cms"``, ``"atlas"``, ``"alice"``, ``"minimal"``.
    cfg : dict, optional
        Full config dict — used to read custom overrides from the
        ``plotting:`` section.
    """
    # Use Agg backend for batch mode unless already interactive
    if not matplotlib.is_interactive():
        matplotlib.use("Agg")

    if name in ("cms", "atlas", "alice"):
        _mplhep_style(name)
        # Still apply ESS tweaks on top for tick direction etc.
        plt.rcParams.update(_ess_params(cfg))
    elif name == "minimal":
        plt.rcParams.update({
            "font.family": "sans-serif",
            "font.size": 11,
            "axes.linewidth": 0.5,
            "xtick.direction": "out",
            "ytick.direction": "out",
            "figure.facecolor": "white",
        })
    else:
        # Default: ESS
        plt.rcParams.update(_ess_params(cfg))


def get_colors(cfg: dict | None = None) -> list[str]:
    """Return the series colour list from config, or the default."""
    if cfg:
        p = cfg.get("plotting", {})
        return p.get("series_colors", SERIES_COLORS)
    return SERIES_COLORS


def get_hist_colors(cfg: dict | None = None) -> tuple[str, str]:
    """Return (fit_color, tail_color) from config."""
    if cfg:
        h = cfg.get("plotting", {}).get("histogram", {})
        return (h.get("fit_color", FIT_COLOR),
                h.get("tail_color", TAIL_COLOR))
    return FIT_COLOR, TAIL_COLOR


# ═══════════════════════════════════════════════════════════════════════════════
# Logo helper
# ═══════════════════════════════════════════════════════════════════════════════

def add_logo(ax: plt.Axes, cfg: dict | None = None) -> None:
    """Add the ESS logo to an axes if configured.

    Reads logo path, zoom, and position from ``cfg["plotting"]["logo"]``.
    Silently does nothing if the logo file is missing or disabled.
    """
    if cfg:
        logo_cfg = cfg.get("plotting", {}).get("logo", {})
        if not logo_cfg.get("show", True):
            return
        logo_path = cfg.get("paths", {}).get("logo", "assets/ess_logo.png")
        zoom = logo_cfg.get("zoom", 0.15)
        pos  = logo_cfg.get("position", [0.03, 0.95])
    else:
        logo_path = "assets/ess_logo.png"
        zoom = 0.15
        pos  = [0.03, 0.95]

    try:
        logo = mpimg.imread(logo_path)
        ab = AnnotationBbox(
            OffsetImage(logo, zoom=zoom),
            tuple(pos),
            xycoords="axes fraction",
            box_alignment=(0, 1),
            frameon=False,
        )
        ax.add_artist(ab)
    except (FileNotFoundError, OSError):
        pass


# ═══════════════════════════════════════════════════════════════════════════════
# Save helper
# ═══════════════════════════════════════════════════════════════════════════════

def save_figure(
    fig: plt.Figure,
    output_path: str | Path,
    formats: list[str] | None = None,
    dpi: int = 300,
    close: bool = True,
) -> list[Path]:
    """Save a figure in one or more formats.

    Parameters
    ----------
    fig : Figure
        Matplotlib figure to save.
    output_path : str or Path
        Base path (without extension).  Extensions are added per format.
    formats : list[str]
        File formats, e.g. ``["pdf", "png"]``.  Defaults to ``["pdf", "png"]``.
    dpi : int
        Resolution for raster formats.
    close : bool
        Close the figure after saving to free memory.

    Returns
    -------
    list[Path]
        Paths of saved files.
    """
    if formats is None:
        formats = ["pdf", "png"]

    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    saved = []
    stem = output_path.stem
    parent = output_path.parent

    for fmt in formats:
        path = parent / f"{stem}.{fmt}"
        fig.savefig(str(path), dpi=dpi, bbox_inches="tight")
        print(f"  Saved: {path}")
        saved.append(path)

    if close:
        plt.close(fig)

    return saved
