"""
hibeam.physics.tracks — 3D track reconstruction and quality cuts
==================================================================
Functions for reconstructing track coordinates from raw TPC hits
and applying quality selections.

Consolidates logic from:
  - experimental_analysis/analysis2.py   (calculate_dEdx with truncated mean)
  - experimental_analysis/energyattempt.py  (cluster-based track finding)
"""

from __future__ import annotations

import numpy as np


def cluster_hits(
    row: np.ndarray,
    col: np.ndarray,
    time: np.ndarray,
    signal: np.ndarray,
    pitch_x: float = 1.0,
    pitch_y: float = 1.0,
    pitch_z: float = 1.0,
) -> dict:
    """Cluster raw TPC hits by pad row using charge-weighted centroids.

    For each unique pad row, computes the charge-weighted average
    position in (column, timestamp) space, then converts to physical
    coordinates using the pitch constants.

    Parameters
    ----------
    row, col, time, signal : np.ndarray
        Per-hit arrays (after noise filtering).
    pitch_x, pitch_y, pitch_z : float
        Physical pitch [mm] per column, row, and timestamp tick.

    Returns
    -------
    dict with arrays: x, y, z, energy (cluster centroids).
    """
    unique_rows = np.unique(row)
    cx, cy, cz, cE = [], [], [], []

    for r in unique_rows:
        idx = row == r
        hits_col = col[idx]
        hits_time = time[idx]
        hits_sig = signal[idx]

        total_E = np.sum(hits_sig)
        if total_E <= 0:
            continue

        avg_col  = np.sum(hits_col * hits_sig) / total_E
        avg_time = np.sum(hits_time * hits_sig) / total_E

        cx.append(avg_col * pitch_x)
        cy.append(float(r) * pitch_y)
        cz.append(avg_time * pitch_z)
        cE.append(total_E)

    return {
        "x": np.array(cx),
        "y": np.array(cy),
        "z": np.array(cz),
        "energy": np.array(cE),
    }


def track_length_3d(x: np.ndarray, y: np.ndarray, z: np.ndarray) -> float:
    """Total 3D path length of a track (sum of step lengths)."""
    if len(x) < 2:
        return 0.0
    dx = np.diff(x)
    dy = np.diff(y)
    dz = np.diff(z)
    return float(np.sum(np.sqrt(dx**2 + dy**2 + dz**2)))


def endpoint_length(x: np.ndarray, y: np.ndarray, z: np.ndarray) -> float:
    """Straight-line distance between first and last hits."""
    if len(x) < 2:
        return 0.0
    return float(np.sqrt(
        (x[-1] - x[0])**2 + (y[-1] - y[0])**2 + (z[-1] - z[0])**2
    ))


def chi2_quality(track: dict) -> float:
    """Compute χ²/ndf for a track from centroid residuals.

    Uses the linear fit parameters (slope_xy, intercept_xy, slope_zy,
    intercept_zy) stored in the track dict by exp_loader.
    """
    x = np.array(track["x"])
    y = np.array(track["y"])
    z = np.array(track["z"])
    n = len(x)

    if n < 3:
        return float("inf")

    # Fit x vs y and z vs y
    if "slope_xy" in track:
        slope_xy = track["slope_xy"]
        intercept_xy = track["intercept_xy"]
        slope_zy = track["slope_zy"]
        intercept_zy = track["intercept_zy"]
    else:
        # Compute from scratch
        slope_xy, intercept_xy = np.polyfit(y, x, 1)
        slope_zy, intercept_zy = np.polyfit(y, z, 1)

    res_x = x - (slope_xy * y + intercept_xy)
    res_z = z - (slope_zy * y + intercept_zy)

    sx = np.std(res_x) if np.std(res_x) > 0 else 1e-4
    sz = np.std(res_z) if np.std(res_z) > 0 else 1e-4

    chi2 = np.sum((res_x / sx)**2 + (res_z / sz)**2)
    ndf = max(2 * n - 4, 1)

    return chi2 / ndf
