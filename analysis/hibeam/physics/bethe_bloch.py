"""
hibeam.physics.bethe_bloch — Bethe-Bloch mean energy-loss theory
==================================================================
Parameterisation of the mean energy loss per unit path length for
charged particles in a gaseous TPC medium.

References
----------
[BB]  Bethe & Bloch (1930/1933). Original stopping-power papers.
[PDG] PDG (2022). Section 34 — Passage of particles through matter.
[STI] Sternheimer (1971). Phys. Rev. B 3, 3681. Density effect.
"""

from __future__ import annotations

import numpy as np


# ── Physical constants ────────────────────────────────────────────────────────

K   = 0.307075       # Universal BB constant [MeV mol⁻¹ cm²]  (PDG Table 34.1)
M_E = 0.51099895     # Electron mass [MeV/c²]


# ── Default medium: Ar/CO₂ (90/10) — standard TPC fill gas ───────────────────

MEDIUM_AR_CO2 = {
    "name":    r"Ar/CO$_2$ 90/10",
    "Z_over_A": 0.45594,      # effective Z/A [mol/g]
    "I":        171.0e-6,     # mean excitation energy [MeV]
    "density":  1.662e-3,     # g/cm³ at STP
    "C_bar":    11.948,       # Sternheimer parameters
    "x0":       1.7635,
    "x1":       4.4855,
    "a":        0.19714,
    "k":        2.9618,
}


# ── Particle species ─────────────────────────────────────────────────────────

PARTICLES = {
    "electron": {"mass": M_E,      "z": 1, "color": "#e41a1c"},
    "pion":     {"mass": 139.570,   "z": 1, "color": "#377eb8"},
    "kaon":     {"mass": 493.677,   "z": 1, "color": "#4daf4a"},
    "proton":   {"mass": 938.272,   "z": 1, "color": "#984ea3"},
}


# ═══════════════════════════════════════════════════════════════════════════════
# Density effect correction
# ═══════════════════════════════════════════════════════════════════════════════

def density_effect(beta_gamma: np.ndarray, medium: dict) -> np.ndarray:
    """Sternheimer density-effect correction δ(βγ).

    Reduces energy loss at high momenta (Fermi plateau) by accounting
    for medium polarisation.  [STI, PDG eq. 34.6]
    """
    x     = np.log10(beta_gamma)
    C_bar = medium["C_bar"]
    x0, x1, a, k = medium["x0"], medium["x1"], medium["a"], medium["k"]

    safe_arg = np.maximum(x1 - x, 0.0)

    return np.where(
        x >= x1,
        2.0 * np.log(10.0) * x - C_bar,
        np.where(
            x >= x0,
            2.0 * np.log(10.0) * x - C_bar + a * safe_arg**k,
            0.0,
        ),
    )


# ═══════════════════════════════════════════════════════════════════════════════
# Bethe-Bloch formula
# ═══════════════════════════════════════════════════════════════════════════════

def bethe_bloch(
    momentum: np.ndarray,
    mass: float,
    z: float = 1.0,
    medium: dict | None = None,
) -> np.ndarray:
    """Mean energy loss ⟨-dE/dx⟩ [MeV/cm] from the Bethe-Bloch formula.

    Parameters
    ----------
    momentum : np.ndarray
        Particle momentum [MeV/c].
    mass : float
        Particle rest mass [MeV/c²].
    z : float
        Charge number (1 for π, K, p).
    medium : dict or None
        Absorber properties.  Defaults to Ar/CO₂ 90/10.

    Returns
    -------
    np.ndarray
        ⟨-dE/dx⟩ [MeV/cm] — always positive.
    """
    if medium is None:
        medium = MEDIUM_AR_CO2

    bg    = momentum / mass
    gamma = np.sqrt(1.0 + bg**2)
    beta  = bg / gamma
    beta2 = beta**2

    Z_A = medium["Z_over_A"]
    I   = medium["I"]

    T_max = (2.0 * M_E * bg**2) / (1.0 + 2.0 * gamma * M_E / mass + (M_E / mass)**2)
    delta = density_effect(bg, medium)

    prefactor = K * z**2 * Z_A / beta2
    bracket   = (0.5 * np.log(2.0 * M_E * bg**2 * T_max / I**2)
                 - beta2 - delta / 2.0)

    return prefactor * bracket


def bethe_bloch_all_particles(
    momentum: np.ndarray,
    medium: dict | None = None,
    path_length_cm: float = 1.0,
) -> dict[str, np.ndarray]:
    """Compute BB curves for all standard particles.

    Returns dict mapping particle name → dE/dx × path_length [MeV].
    """
    result = {}
    for name, props in PARTICLES.items():
        dedx_cm = bethe_bloch(momentum, props["mass"], props["z"], medium)
        valid = dedx_cm > 0
        dedx = np.where(valid, dedx_cm * path_length_cm, np.nan)
        result[name] = dedx
    return result
