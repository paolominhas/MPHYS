"""
hibeam.physics.dedx — Energy-loss computation
================================================
All methods for computing dE/dx from TPC data.  No plotting, no I/O.

Methods
-------
sum_edep       — Sum of Geant4 step energy deposits per event [MeV].
                 Proportional to dE/dx for a fixed-depth detector.
per_centroid   — ADC[j] / (track_length / nPoints) for experimental data.
truncated_mean — Sort dE/dx samples low→high, average the bottom fraction.

References
----------
[1] ALICE Collaboration (2010). JINST 5, P09002.
    Defines the standard 70% truncation fraction.
[2] Bichsel, H. (1988). Rev. Mod. Phys. 60, 663.
[3] STAR Collaboration (2003). NIM A 499, 659.
"""

from __future__ import annotations

from typing import Any, Optional

import awkward as ak
import numpy as np

from hibeam.config import TPCGeometry


# ═══════════════════════════════════════════════════════════════════════════════
# Track length estimation
# ═══════════════════════════════════════════════════════════════════════════════

def track_length_mm(
    pad_row: ak.Array,
    pad_col: ak.Array,
    timestamp: ak.Array,
    geom: TPCGeometry,
) -> float:
    """Estimate straight-track length [mm] from the 3D extent of pad hits.

    Physical coordinates:
        x = pad_column × pad_width   [mm]
        y = pad_row    × pad_height  [mm]
        z = timestamp  × drift_v × time_bin  [mm]
    """
    if len(pad_row) < 2:
        return float("nan")
    dx = (float(ak.max(pad_col))   - float(ak.min(pad_col)))   * geom.pad_width
    dy = (float(ak.max(pad_row))   - float(ak.min(pad_row)))   * geom.pad_height
    dz = (float(ak.max(timestamp)) - float(ak.min(timestamp))) * geom.drift_v * geom.time_bin
    return float(np.sqrt(dx**2 + dy**2 + dz**2))


# ═══════════════════════════════════════════════════════════════════════════════
# Main dE/dx dispatcher
# ═══════════════════════════════════════════════════════════════════════════════

def compute_dedx(
    data: dict[str, Any],
    method: str = "sum_edep",
    geom: Optional[TPCGeometry] = None,
    min_steps: int = 5,
    truncate_fraction: float = 0.3,
    low_cut_mev: float = 0.0,
) -> np.ndarray:
    """Compute per-event energy-loss values from TPC data.

    Parameters
    ----------
    data : dict
        Output of ``sim_loader.load()`` or ``exp_loader.load()``.
    method : str
        Computation method:
        - ``"sum_edep"``       — total Edep per event (simulation).
        - ``"per_centroid"``   — ADC / dx per centroid (experimental).
        - ``"truncated_mean"`` — truncated mean of dE/dx samples.
    geom : TPCGeometry or None
        If complete, enables true dE/dx [MeV/mm] for sum_edep method.
    min_steps : int
        Minimum G4 steps (or hits) per event to accept.
    truncate_fraction : float
        For ``"truncated_mean"``: discard the top this fraction.
    low_cut_mev : float
        Discard events with total Edep below this [MeV].

    Returns
    -------
    np.ndarray
        Per-event energy loss values.
    """
    if method == "sum_edep":
        return _sum_edep(data, geom=geom, min_steps=min_steps,
                         low_cut_mev=low_cut_mev)
    elif method == "per_centroid":
        # For experimental data, dedx is pre-computed by exp_loader
        return data["dedx"]
    elif method == "truncated_mean":
        return _truncated_mean_dedx(data, truncate_fraction=truncate_fraction,
                                    min_steps=min_steps)
    else:
        raise ValueError(
            f"Unknown method '{method}'.  "
            f"Options: sum_edep, per_centroid, truncated_mean."
        )


# ═══════════════════════════════════════════════════════════════════════════════
# Method implementations
# ═══════════════════════════════════════════════════════════════════════════════

def _sum_edep(
    data: dict,
    geom: Optional[TPCGeometry] = None,
    min_steps: int = 5,
    low_cut_mev: float = 0.0,
) -> np.ndarray:
    """Sum of G4-step energy deposits per event [MeV or MeV/mm].

    If a complete TPCGeometry is supplied AND pad branches are present,
    returns dE/dx [MeV/mm].  Otherwise returns sum(Edep) [MeV].
    """
    edep    = data["edep"]
    pad_row = data.get("pad_row")
    pad_col = data.get("pad_col")
    ts      = data.get("timestamp")

    use_geom = (
        geom is not None
        and geom.is_complete
        and pad_row is not None
        and pad_col is not None
        and ts is not None
    )

    result = []
    for i in range(data["n_events"]):
        ev_edep = edep[i]
        if len(ev_edep) < min_steps:
            continue
        e_total = float(ak.sum(ev_edep))

        if low_cut_mev > 0 and e_total < low_cut_mev / 1000.0:
            continue

        if not use_geom:
            result.append(e_total)
        else:
            length = track_length_mm(pad_row[i], pad_col[i], ts[i], geom)
            if np.isnan(length) or length <= 0.0:
                continue
            result.append(e_total / length)

    return np.asarray(result, dtype=np.float64)


def _truncated_mean_dedx(
    data: dict,
    truncate_fraction: float = 0.3,
    min_steps: int = 5,
) -> np.ndarray:
    """Truncated mean of per-hit dE/dx samples.

    For each event:
      1. Compute dr = 3D distance between consecutive hits.
      2. dE/dx = charge[i] / dr[i] for each pair.
      3. Sort ascending, average the bottom (1 - truncate_fraction).

    This is the standard method from ALICE JINST 5 (2010) P09002.
    """
    tracks = data.get("tracks", [])
    if not tracks:
        raise ValueError("No track data available for truncated_mean method.")

    results = []
    for track in tracks:
        x = np.array(track["x"])
        y = np.array(track["y"])
        z = np.array(track["z"])
        q = np.array(track["charge"])

        if len(x) < min_steps:
            continue

        # 3D path length between consecutive hits
        dr = np.sqrt(np.diff(x)**2 + np.diff(y)**2 + np.diff(z)**2)
        mask = dr > 0
        dr = dr[mask]
        dE = q[:-1][mask]

        if len(dr) < 3:
            continue

        dedx_samples = dE / dr
        sorted_samples = np.sort(dedx_samples)
        cutoff = int(len(sorted_samples) * (1 - truncate_fraction))
        results.append(np.mean(sorted_samples[:cutoff]))

    return np.asarray(results, dtype=np.float64)


# ═══════════════════════════════════════════════════════════════════════════════
# Peak-normalisation for data-vs-sim comparisons
# ═══════════════════════════════════════════════════════════════════════════════

def peak_normalise(values: np.ndarray, bins: int = 200,
                   weights: np.ndarray | None = None) -> tuple[np.ndarray, float]:
    """Divide values by their histogram peak position.

    Returns (normalised_values, peak_position).
    Makes the x-axis dimensionless (units of MPV) for comparing
    data and simulation without absolute calibration.
    """
    from hibeam.utils import histogram_peak
    peak = histogram_peak(values, bins=bins, weights=weights)
    return values / peak, peak
