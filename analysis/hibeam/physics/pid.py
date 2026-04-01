"""
hibeam.physics.pid — ΔE-E telescope particle identification
==============================================================
Computes the observables for a 2D ΔE-E plot using the Proto-TPC
(thin absorber, ΔE) and scintillator stack (thick absorber, E).

Different particle species trace distinct hyperbolic bands in the
ΔE vs E_residual plane, enabling Z and A identification.

References
----------
[1] Bethe, Ann. Phys. 5 (1930) 325.
[2] Carboni et al., NIM A 664 (2012) 251 — ΔE-E with FAZIA.
[3] Leo, Techniques for Nuclear and Particle Physics, Ch. 2.
"""

from __future__ import annotations

import numpy as np


def compute_pid_observables(
    tpc_mev: np.ndarray,
    cal_mev: np.ndarray,
    n_tpc: np.ndarray,
    n_cal: np.ndarray,
    min_tpc_hits: int = 1,
    min_scint_hits: int = 1,
) -> dict:
    """Select events with hits in both detectors and return ΔE vs E.

    Parameters
    ----------
    tpc_mev : np.ndarray
        Total TPC energy deposit per event [MeV].
    cal_mev : np.ndarray
        Total calorimeter deposit per event [MeV].
    n_tpc, n_cal : np.ndarray
        Hit counts per event.
    min_tpc_hits, min_scint_hits : int
        Minimum hit requirements.

    Returns
    -------
    dict
        delta_e   : np.ndarray — TPC deposit (ΔE) for selected events.
        e_residual: np.ndarray — calorimeter deposit (E) for selected events.
        mask      : np.ndarray[bool] — selection mask.
        n_selected: int
    """
    mask = (n_tpc >= min_tpc_hits) & (n_cal >= min_scint_hits)
    mask &= (tpc_mev > 0) & (cal_mev > 0)

    return {
        "delta_e":    tpc_mev[mask],
        "e_residual": cal_mev[mask],
        "mask":       mask,
        "n_selected": int(mask.sum()),
    }
