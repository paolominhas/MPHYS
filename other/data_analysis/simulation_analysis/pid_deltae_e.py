#!/usr/bin/env python3
"""
pid_deltaE_E.py
===============
Particle identification via the ΔE-E telescope technique using the
HIBEAM Proto-TPC (ΔE) and scintillator stack (E).

Physics
-------
When a charged particle traverses two detector layers of different
thickness, the energy deposited in each depends on the particle's
charge Z, mass A, and velocity β through the Bethe-Bloch formula:

    -dE/dx ∝ Z² / β²  × [ ln(2mₑc²β²γ²/I) - β² ]

Different particle species therefore trace distinct loci (hyperbolic
bands) in a 2D plot of ΔE (thin absorber, here the TPC) vs E_residual
(thick absorber, here the scintillator stack).  Heavier or more-charged
particles deposit more energy in the thin layer for a given residual
energy, producing separated bands.

This is the standard ΔE-E telescope method used since the 1960s in
nuclear and particle physics for Z and A identification.

References
----------
[1] Bethe H., Ann. Phys. 5 (1930) 325  — original stopping-power theory
[2] Carboni S. et al., NIM A 664 (2012) 251  — ΔE-E + PSA with FAZIA Si telescopes
[3] Pouthas J. et al., NIM A 357 (1995) 418  — INDRA: 4π charged-particle array
[4] Leo W.R., Techniques for Nuclear and Particle Physics Experiments,
    Springer (1994) Ch. 2  — textbook treatment of ΔE-E method
[5] PDG, Phys. Rev. D 98 (2018) 030001, Sec. 34  — passage of particles through matter
[6] Grupen C. & Shwartz B., Particle Detectors, Cambridge (2008) Sec. 2.1
[7] Bichsel H., Rev. Mod. Phys. 60 (1988) 663  — stochastic energy loss theory

Data structure (from TFile::MakeProject headers)
-------------------------------------------------
Tree "hibeam" contains per-event branches:
  ProtoTPC/Edep   : vector<double>  — energy deposits per G4 step in TPC [GeV]
  ProtoTPC/SumEdep: double          — pre-summed TPC deposit [GeV]
  HRD/Edep        : vector<double>  — energy deposits in calorimeter (HRD) [GeV]

Usage
-----
    python3 pid_deltaE_E.py
    Run from the combined_analysis/ directory.
"""

import sys
from pathlib import Path

import awkward as ak
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.ticker as ticker
import mplhep as hep
import numpy as np
import uproot

# ══════════════════════════════════════════════════════════════════════════════
# CONFIG
# ══════════════════════════════════════════════════════════════════════════════
SIM_FILES = [
    {
        "file"  : "../simulation_data/KrakowScatter.root",
        "label" : r"Krakow $^2$H(p,p) elastic",
        "output": "pid_krakow",
    },
    {
        "file"  : "../simulation_data/MuonScatter_fixed.root",
        "label" : "Muon lab (MCPL)",
        "output": "pid_muon",
    },
]

MIN_TPC_HITS   = 1    # require at least 1 TPC hit
MIN_SCINT_HITS = 1    # require at least 1 scintillator hit
OUTDIR         = "."

# ══════════════════════════════════════════════════════════════════════════════
# STYLE — publication quality, matching the dE/dx comparison scripts
# ══════════════════════════════════════════════════════════════════════════════
plt.style.use(hep.style.CMS)
plt.rcParams.update({
    "font.size"       : 12,
    "axes.labelsize"  : 14,
    "xtick.labelsize" : 11,
    "ytick.labelsize" : 11,
    "xtick.direction" : "in",
    "ytick.direction" : "in",
    "xtick.top"       : True,
    "ytick.right"     : True,
})

# ══════════════════════════════════════════════════════════════════════════════
# DATA LOADING
# ══════════════════════════════════════════════════════════════════════════════
def find_branch(tree, candidates):
    """Return the first matching branch key from a list of candidates."""
    keys = set(tree.keys(recursive=True))
    for c in candidates:
        if c in keys:
            return c
    return None


def load_sim(filepath):
    """
    Load per-event TPC and calorimeter energy deposits.

    Unit convention: ProtoTPC/Edep is stored in GeV (Geant4 HEP default),
    while HRD/Edep may be in MeV.  We auto-detect by checking the raw
    magnitude: if max summed deposit > 500, it's already MeV.

    Returns
    -------
    tpc_mev  : 1D array — total energy deposited in TPC per event [MeV]
    cal_mev  : 1D array — total energy deposited in calorimeter per event [MeV]
    n_tpc    : 1D array — number of TPC hits per event
    n_cal    : 1D array — number of calorimeter hits per event
    """
    print(f"  Loading: {Path(filepath).name}")
    with uproot.open(filepath) as f:
        tree = f["hibeam"]

        # ── TPC energy deposit ────────────────────────────────────────────
        tpc_key = find_branch(tree, [
            "ProtoTPC/Edep", "ProtoTPC.Edep", "TPC/Edep", "TPC.Edep"])
        if tpc_key is None:
            raise KeyError(f"No TPC Edep branch. Available: "
                           f"{[k for k in tree.keys(recursive=True) if 'dep' in k.lower()]}")

        tpc_edep = tree[tpc_key].array(library="ak")
        n_tpc    = ak.to_numpy(ak.num(tpc_edep)).astype(int)
        tpc_raw  = ak.to_numpy(ak.sum(tpc_edep, axis=1)).astype(float)

        # ── Calorimeter / HRD energy deposit ──────────────────────────────
        cal_key = find_branch(tree, [
            "HRD/Edep", "HRD.Edep",
            "ScintHit/Edep", "ScintHit.Edep", "Scint/Edep", "Scint.Edep"])
        if cal_key is None:
            raise KeyError(f"No calorimeter Edep branch. Available: "
                           f"{[k for k in tree.keys(recursive=True)]}")

        cal_edep = tree[cal_key].array(library="ak")
        n_cal    = ak.to_numpy(ak.num(cal_edep)).astype(int)
        cal_raw  = ak.to_numpy(ak.sum(cal_edep, axis=1)).astype(float)

    # ── Auto-detect units: if max > 500 raw, assume already MeV ──────────
    def to_mev(arr, label):
        mx = float(arr[arr > 0].max()) if (arr > 0).any() else 0
        if mx > 500:
            print(f"    {label}: raw max = {mx:.1f} → assuming MeV (no conversion)")
            return arr
        else:
            print(f"    {label}: raw max = {mx:.4f} → assuming GeV, ×1000 → MeV")
            return arr * 1000.

    tpc_mev = to_mev(tpc_raw, f"TPC  ({tpc_key})")
    cal_mev = to_mev(cal_raw, f"Cal  ({cal_key})")

    print(f"    {len(tpc_mev):,} events total")
    m1 = n_tpc >= 1
    m2 = n_cal >= 1
    if m1.any():
        print(f"    TPC:  {m1.sum():,} with ≥1 hit,  "
              f"range {tpc_mev[m1].min():.3f}–{tpc_mev[m1].max():.1f} MeV")
    if m2.any():
        print(f"    Cal:  {m2.sum():,} with ≥1 hit,  "
              f"range {cal_mev[m2].min():.3f}–{cal_mev[m2].max():.1f} MeV")

    return tpc_mev, cal_mev, n_tpc, n_cal


# ══════════════════════════════════════════════════════════════════════════════
# PLOT
# ══════════════════════════════════════════════════════════════════════════════
def make_pid_plot(cfg, tpc_mev, cal_mev, n_tpc, n_cal):
    """
    2D histogram:  x = E_calorimeter (HRD total deposit)
                   y = E_TPC  (≈ ΔE in the thin gas layer)

    Standard ΔE-E presentation: the thin-detector signal on the y-axis,
    the residual/total energy on the x-axis [Leo 1994, Fig. 2.8].
    """
    # ── Event selection ───────────────────────────────────────────────────
    mask = (n_tpc >= MIN_TPC_HITS) & (n_cal >= MIN_SCINT_HITS) \
         & (tpc_mev > 0) & (cal_mev > 0)
    x = cal_mev[mask]
    y = tpc_mev[mask]
    print(f"  {mask.sum():,} events pass both-detector requirement")

    if mask.sum() < 50:
        print("  SKIP — too few events with hits in both detectors")
        return

    # ── Axis ranges: clip extreme tails ───────────────────────────────────
    x_hi = float(np.percentile(x, 99.5))
    y_hi = float(np.percentile(y, 99.5))

    # ── 2D histogram ──────────────────────────────────────────────────────
    fig, ax = plt.subplots(figsize=(7.5, 6.2))

    h = ax.hist2d(x, y, bins=[200, 200],
                  range=[[0, x_hi], [0, y_hi]],
                  norm=mcolors.LogNorm(vmin=1),
                  cmap="inferno", rasterized=True)

    cb = fig.colorbar(h[3], ax=ax, pad=0.015)
    cb.set_label("Events per bin", fontsize=12)

    # ── Annotations ───────────────────────────────────────────────────────
    ax.set_xlabel(r"$E_{\mathrm{HRD}}$  [MeV]  (calorimeter / scint. stack)", fontsize=14)
    ax.set_ylabel(r"$E_{\mathrm{TPC}}$  [MeV]  ($\approx\Delta E$)", fontsize=14)
    ax.set_title(f"{cfg['label']}  —  $\\Delta E$–$E$ particle ID",
                 fontsize=13, pad=6)

    # Physics note
    ax.text(0.03, 0.97,
            "$\\Delta E$–$E$ telescope method\n"
            "[Leo (1994) Ch. 2; Carboni et al.,\n"
            " NIM A 664 (2012) 251]\n"
            "Distinct bands $\\Leftrightarrow$ distinct species\n"
            "(separated by $Z$, $A$, $\\beta$)",
            transform=ax.transAxes, va="top", ha="left",
            fontsize=8, style="italic",
            bbox=dict(boxstyle="round,pad=0.4", facecolor="white",
                      edgecolor="#aaaaaa", alpha=0.90))

    # Event count
    ax.text(0.97, 0.03,
            f"{mask.sum():,} events",
            transform=ax.transAxes, va="bottom", ha="right",
            fontsize=9, color="white",
            bbox=dict(boxstyle="round,pad=0.3", facecolor="black",
                      alpha=0.5))

    # ESS branding
    try:
        lo = hep.cms.label(loc=0, data=False, label="Simulation",
                           rlabel="HIBEAM", ax=ax)
        lo[0].set_text("ESS")
    except Exception:
        pass

    plt.tight_layout()
    for ext in [".pdf", ".png"]:
        path = f"{OUTDIR}/{cfg['output']}{ext}"
        fig.savefig(path, dpi=300 if ext == ".png" else None,
                    bbox_inches="tight")
        print(f"    Saved: {path}")
    plt.close(fig)


# ══════════════════════════════════════════════════════════════════════════════
# MAIN
# ══════════════════════════════════════════════════════════════════════════════
for cfg in SIM_FILES:
    print(f"\n{'='*55}")
    print(f"  {cfg['label']}")
    print(f"{'='*55}")
    try:
        tpc, cal, nt, nc = load_sim(cfg["file"])
    except Exception as e:
        print(f"  SKIP — {e}")
        continue
    make_pid_plot(cfg, tpc, cal, nt, nc)

print("\nDone.")