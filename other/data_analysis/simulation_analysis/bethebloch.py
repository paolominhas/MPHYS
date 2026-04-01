#!/usr/bin/env python3
"""
bethe_bloch.py
==============
Bethe-Bloch mean energy-loss parameterisation for charged particles in
a gaseous TPC medium, with a simple overlay plot against the measured
dE/dx distribution from the HIBEAM simulation.

Imports load_tpc_data and compute_dedx from tpc_dedx_analysis.py so that
all data loading logic lives in one place.

Usage
-----
    python bethe_bloch.py /path/to/HIBEAMScatter.root [output.png]

References
----------
[BB]  Bethe, H. & Bloch, F. (1930/1933).  Original stopping-power papers.
[PDG] Particle Data Group (2022). Prog. Theor. Exp. Phys. 2022, 083C01.
      Section 34 — Passage of particles through matter.
      https://pdg.lbl.gov/2022/reviews/rpp2022-rev-passage-particles-matter.pdf
[STI] Sternheimer, R.M. (1971). Phys. Rev. B 3, 3681.
      Density effect parameterisation used here.
"""

import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

# Import data-loading utilities from the main analysis script.
# Make sure tpc_dedx_analysis.py is in the same directory (or on PYTHONPATH).
from other.data_analysis.simulation_analysis.final import load_tpc_data, compute_dedx


# ══════════════════════════════════════════════════════════════════════════════
# SECTION 1 — PHYSICAL CONSTANTS AND MATERIAL PROPERTIES
# ══════════════════════════════════════════════════════════════════════════════

# Universal Bethe-Bloch constant  K = 4πN_A r_e² m_e c²  [MeV mol⁻¹ cm²]
K = 0.307075          # PDG Table 34.1

# Electron mass [MeV/c²]
M_E = 0.51099895

# Speed of light — used only for unit clarity; β is dimensionless
C = 1.0


# Absorber medium: Ar/CO2 (90/10) — standard TPC fill gas
# Values from PDG Table 34.2 and NIST PSTAR database
MEDIUM = {
    "name"  : "Ar/CO$_2$ 90/10",
    "Z_over_A" : 0.45594,   # effective Z/A  [mol g⁻¹]
    "I"        : 171.0e-6,  # mean excitation energy  [MeV]  (Ar: 188 eV, CO2: 85 eV → mixture)
    "density"  : 1.662e-3,  # g/cm³ at STP
    # Sternheimer density-effect parameters [STI]
    "C_bar" : 11.948,
    "x0"    :  1.7635,
    "x1"    :  4.4855,
    "a"     :  0.19714,
    "k"     :  2.9618,
}


# Particle species: (name, mass [MeV/c²], charge z)
PARTICLES = {
    "electron" : (M_E,      1, "#e41a1c"),   # red
    "pion"     : (139.570,  1, "#377eb8"),   # blue
    "kaon"     : (493.677,  1, "#4daf4a"),   # green
    "proton"   : (938.272,  1, "#984ea3"),   # purple
}


# ══════════════════════════════════════════════════════════════════════════════
# SECTION 2 — BETHE-BLOCH FUNCTION
# ══════════════════════════════════════════════════════════════════════════════

def density_effect(beta_gamma: np.ndarray, medium: dict) -> np.ndarray:
    """
    Sternheimer density-effect correction δ(βγ) [STI, PDG eq. 34.6].

    This correction reduces the energy loss at high momenta (the
    'Fermi plateau') by accounting for the polarisation of the medium.
    Without it, Bethe-Bloch rises logarithmically forever, which is
    unphysical.
    """
    x     = np.log10(beta_gamma)
    C_bar = medium["C_bar"]
    x0    = medium["x0"]
    x1    = medium["x1"]
    a     = medium["a"]
    k     = medium["k"]

    # np.where evaluates ALL branches before selecting, so (x1-x)**k is
    # computed even where x > x1, giving a negative base to a fractional
    # power → NaN.  Clip to zero first so the arithmetic is always safe;
    # the np.where then discards those values anyway.
    safe_arg = np.maximum(x1 - x, 0.0)

    delta = np.where(
        x >= x1,
        2.0 * np.log(10.0) * x - C_bar,
        np.where(
            x >= x0,
            2.0 * np.log(10.0) * x - C_bar + a * safe_arg ** k,
            0.0,
        )
    )
    return delta


def bethe_bloch(
    momentum: np.ndarray,
    mass: float,
    z: float,
    medium: dict,
) -> np.ndarray:
    """
    Mean energy loss per unit path length  -<dE/dx>  [MeV/cm]
    from the Bethe-Bloch formula [PDG eq. 34.5].

    Parameters
    ----------
    momentum : np.ndarray
        Particle momentum  [MeV/c].
    mass : float
        Particle rest mass  [MeV/c²].
    z : float
        Particle charge number (1 for π, K, p).
    medium : dict
        Absorber properties (use the MEDIUM dict defined above).

    Returns
    -------
    np.ndarray
        -<dE/dx>  [MeV/cm] — always positive.
    """
    # Relativistic kinematics
    bg    = momentum / mass                           # β·γ
    gamma = np.sqrt(1.0 + bg**2)
    beta  = bg / gamma
    beta2 = beta**2

    Z_A = medium["Z_over_A"]
    I   = medium["I"]                                 # [MeV]

    # Maximum kinetic energy transfer to a free electron [PDG eq. 34.4]
    T_max = (2.0 * M_E * bg**2) / (1.0 + 2.0 * gamma * M_E / mass + (M_E / mass)**2)

    # Density effect
    delta = density_effect(bg, medium)

    # Bethe-Bloch [PDG eq. 34.5]
    prefactor = K * z**2 * Z_A / beta2
    bracket   = (0.5 * np.log(2.0 * M_E * bg**2 * T_max / I**2)
                 - beta2
                 - delta / 2.0)

    return prefactor * bracket


# ══════════════════════════════════════════════════════════════════════════════
# SECTION 3 — PLOTTING
# ══════════════════════════════════════════════════════════════════════════════

def plot_bethe_bloch(
    measured_values: np.ndarray,
    medium: dict = MEDIUM,
    momentum_range: tuple = (10.0, 1e5),   # MeV/c
    effective_track_length_cm: float = 0.6,
    output_path: str = None,
) -> plt.Figure:
    """
    Two-panel figure:

    Top panel — Bethe-Bloch curves for e, π, K, p vs momentum [MeV/c].
    The BB curves (MeV/cm) are multiplied by `effective_track_length_cm`
    to convert to MeV so they can be overlaid with the measured sum(Edep).
    Without full track-length reconstruction this conversion requires an
    assumed path length — one pad row (≈ 0.6 cm) is a reasonable default
    for a TPC.  The assumption is clearly labelled on the plot.

    A horizontal band shows measured mean ± 1σ for direct visual comparison.

    Bottom panel — measured sum(Edep) histogram.

    Parameters
    ----------
    measured_values : np.ndarray
        Per-event energy-loss values from compute_dedx() [MeV].
    medium : dict
        Absorber properties.
    momentum_range : tuple
        (p_min, p_max) in MeV/c for the Bethe-Bloch curves.
    effective_track_length_cm : float
        Path length used to convert BB [MeV/cm] → MeV.
        Default 0.6 cm ≈ one TPC pad row.  Adjust to match your geometry.
    output_path : str or None
        Save path; show interactively if None.
    """
    p = np.logspace(
        np.log10(momentum_range[0]),
        np.log10(momentum_range[1]),
        2000,
    )

    fig = plt.figure(figsize=(10, 7))
    gs  = gridspec.GridSpec(2, 1, height_ratios=[2, 1], hspace=0.36,
                            left=0.10, right=0.97, top=0.93, bottom=0.10)
    ax_bb   = fig.add_subplot(gs[0])
    ax_hist = fig.add_subplot(gs[1])

    # ── Top: Bethe-Bloch curves scaled to MeV ────────────────────────────────
    for name, (mass, z, colour) in PARTICLES.items():
        dedx_per_cm = bethe_bloch(p, mass, z, medium)          # MeV/cm
        dedx_mev    = dedx_per_cm * effective_track_length_cm  # MeV
        valid = dedx_mev > 0
        ax_bb.plot(p[valid], dedx_mev[valid],
                   lw=2.0, color=colour, label=name)

    # Measured mean ± 1σ band — now in the same units (MeV)
    mean_meas = float(np.mean(measured_values))
    std_meas  = float(np.std(measured_values))
    ax_bb.axhspan(
        mean_meas - std_meas, mean_meas + std_meas,
        color="gold", alpha=0.35,
        label=rf"Measured $\langle\Sigma E_{{\rm dep}}\rangle \pm 1\sigma$",
    )
    ax_bb.axhline(mean_meas, color="goldenrod", lw=1.2, ls="--")

    ax_bb.set_xscale("log")
    ax_bb.set_yscale("log")
    ax_bb.set_xlabel("Momentum  [MeV/$c$]", fontsize=12)
    ax_bb.set_ylabel(
        rf"$-\langle dE/dx\rangle \times {effective_track_length_cm}$ cm  [MeV]",
        fontsize=11,
    )
    ax_bb.set_title(
        f"Bethe-Bloch energy loss in {medium['name']}  "
        f"(path length = {effective_track_length_cm} cm assumed)\n"
        "[PDG 2022, Sec. 34;  Sternheimer 1971 density correction]",
        fontsize=10,
    )
    ax_bb.legend(fontsize=9, loc="upper right",
                 framealpha=0.92, edgecolor="#aaaaaa")
    ax_bb.set_xlim(*momentum_range)

    # ── Bottom: measured histogram ────────────────────────────────────────────
    counts, edges = np.histogram(measured_values, bins=120,
                                 range=(0, np.percentile(measured_values, 95)))
    centres = 0.5 * (edges[:-1] + edges[1:])
    bw      = edges[1] - edges[0]
    ax_hist.bar(centres, counts, width=bw * 0.92,
                color="#2166ac", alpha=0.60, label="Measured sum(Edep)")
    ax_hist.axvline(mean_meas, color="goldenrod", lw=1.5, ls="--",
                    label=rf"Mean = {mean_meas:.4g} MeV")
    ax_hist.set_xlabel(r"$\sum E_{\rm dep}$  [MeV]", fontsize=12)
    ax_hist.set_ylabel("Events / bin", fontsize=11)
    ax_hist.legend(fontsize=9, framealpha=0.92, edgecolor="#aaaaaa")

    if output_path:
        import os
        stem, ext = os.path.splitext(output_path)
        if not ext:
            ext = ".png"
        fig.savefig(f"{stem}{ext}", dpi=150, bbox_inches="tight")
        print(f"Figure saved: {stem}{ext}")
    else:
        plt.show()

    return fig


# ══════════════════════════════════════════════════════════════════════════════
# SECTION 4 — MAIN
# ══════════════════════════════════════════════════════════════════════════════

def main(filepath: str, output_path: str = None) -> None:
    print(f"\nLoading TPC data from: {filepath}")
    data   = load_tpc_data(filepath)
    values = compute_dedx(data, geom=None, min_steps=5)
    print(f"  {len(values):,} events accepted")

    print("\nGenerating Bethe-Bloch plot...")
    plot_bethe_bloch(values, output_path=output_path)


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python bethe_bloch.py <file.root> [output.png]")
        sys.exit(1)
    main(sys.argv[1], sys.argv[2] if len(sys.argv) > 2 else None)