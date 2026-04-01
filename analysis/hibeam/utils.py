"""
hibeam.utils — Shared helper functions
=======================================
Small, reusable functions used across io, physics, and plotting modules.
"""

from __future__ import annotations

import glob
import os
import re
import warnings
from pathlib import Path
from typing import Optional

import numpy as np


# ── File discovery ────────────────────────────────────────────────────────────

def find_root_files(directory: str | Path, pattern: str = "*.root") -> list[Path]:
    """Return sorted list of ROOT files matching *pattern* in *directory*."""
    d = Path(directory)
    if not d.is_dir():
        raise FileNotFoundError(f"Directory not found: {d}")
    files = sorted(d.glob(pattern))
    if not files:
        warnings.warn(f"No files matching '{pattern}' in {d}", stacklevel=2)
    return files


def extract_nsec(filename: str | Path) -> Optional[int]:
    """Extract the nSections number from a segmentation filename.

    Matches patterns like ``nSec010``, ``NSec10``, ``nsec0010``.
    Returns None if no match is found.
    """
    m = re.search(r"[nN][sS]ec0*(\d+)", str(filename))
    return int(m.group(1)) if m else None


def make_short_label(path: str | Path, max_chars: int = 45) -> str:
    """Build a human-readable label from a ROOT filename.

    Strips common prefixes like ``tracks_centroids_`` and suffixes like
    ``_sorted``, then truncates to *max_chars*.
    """
    name = Path(path).stem
    for prefix in ("tracks_centroids_", "ana_geom_"):
        name = name.removeprefix(prefix)
    for suffix in ("_sorted", "_sorted_t0_100", "_sorted_t0"):
        name = name.removesuffix(suffix)
    return name[:max_chars]


# ── Array helpers ─────────────────────────────────────────────────────────────

def percentile_clip(arr: np.ndarray, lo: float = 0.0, hi: float = 95.0) -> np.ndarray:
    """Clip array to [lo, hi] percentile range.  Returns a copy."""
    lower = np.percentile(arr, lo) if lo > 0 else arr.min()
    upper = np.percentile(arr, hi)
    return arr[(arr >= lower) & (arr <= upper)].copy()


def safe_divide(a: np.ndarray, b: np.ndarray, fill: float = 0.0) -> np.ndarray:
    """Element-wise a/b, replacing division-by-zero with *fill*."""
    with np.errstate(divide="ignore", invalid="ignore"):
        result = np.where(b != 0, a / b, fill)
    return result


def histogram_peak(values: np.ndarray, bins: int = 200,
                   weights: np.ndarray | None = None) -> float:
    """Return the bin-centre of the histogram peak (coarse mode estimate)."""
    counts, edges = np.histogram(values, bins=bins, weights=weights)
    centres = 0.5 * (edges[:-1] + edges[1:])
    return float(centres[np.argmax(counts)])


# ── Branch discovery (used by multiple loaders) ──────────────────────────────

def find_branch(keys: list[str], candidates: list[str]) -> Optional[str]:
    """Return the first key in *keys* that matches any of *candidates*.

    Matching is done on the last component after ``/`` or ``.`` splitting,
    so ``"ProtoTPC/Edep"`` matches candidate ``"Edep"``.
    """
    key_set = set(keys)
    # First try exact match
    for c in candidates:
        if c in key_set:
            return c
    # Then try suffix match
    for c in candidates:
        for k in keys:
            leaf = k.split("/")[-1].split(".")[-1]
            if leaf == c:
                return k
    return None


# ── Printing helpers ──────────────────────────────────────────────────────────

def print_header(text: str, width: int = 60) -> None:
    """Print a section header with a decorative line."""
    print(f"\n{'─' * width}")
    print(f"  {text}")
    print(f"{'─' * width}")


def print_fit_summary(result: dict) -> None:
    """Pretty-print fit results from fitting.fit_landau()."""
    chi2_tail_red = result["chi2_tail"] / max(result["ndf_tail"], 1)
    print_header("FIT RESULTS")
    print(f"  Most probable value (MPV)  : {result['mpv']:.4g}")
    print(f"  Width  ξ (scale)           : {result['scale']:.4g}"
          f"  ±  {result['scale_err']:.3g}")
    print(f"  Fit region   χ²/ndf        : {result['chi2_red']:.3f}"
          f"    (p = {result['p_value']:.3f})")
    print(f"  Tail region  χ²/ndf        : {chi2_tail_red:.3f}"
          f"    (p = {result['p_tail']:.3f})")
    print("  (large χ²_tail expected: delta-ray enhancement)")
    print("─" * 60)
