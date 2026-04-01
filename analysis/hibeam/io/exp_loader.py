"""
hibeam.io.exp_loader — Load experimental TPC ROOT files
=========================================================
Reads reconstructed track data produced by the TPC reconstruction
chain.  Uses PyROOT with compiled C++ headers (TrackData.h) because
the ``tracks`` branch contains ``std::vector<TrackData>`` objects
that uproot cannot deserialise without the class dictionary.

Also provides an uproot-based loader for the raw ``tpc`` branch
(pad-level hits) used in pad-plane diagnostics and 3D displays.

Consolidates logic previously in:
  - combined_analysis/combining_dedx.py      (load_data)
  - combined_analysis/check.py
  - combined_analysis/combining_dedx_simple.py
  - experimental_analysis/full_data_ana.py
  - experimental_analysis/dx-graph.py
  - experimental_analysis/plot_tracks.py
  - experimental_analysis/testing.py
"""

from __future__ import annotations

import os
import warnings
from pathlib import Path
from typing import Any, Optional

import numpy as np


# ═══════════════════════════════════════════════════════════════════════════════
# PyROOT header management
# ═══════════════════════════════════════════════════════════════════════════════

_HEADERS_LOADED = False


def _ensure_headers(headers_dir: str | Path, headers_so: str | Path) -> None:
    """Load compiled C++ headers into PyROOT exactly once.

    The loading order is critical:
      1. AddIncludePath
      2. ProcessLine #include "TrackData.h"
      3. gSystem.Load the .so

    Getting this wrong causes silent data corruption (all-zero charge arrays).
    """
    global _HEADERS_LOADED
    if _HEADERS_LOADED:
        return

    import ROOT
    ROOT.gROOT.SetBatch(True)
    ROOT.gErrorIgnoreLevel = ROOT.kFatal

    headers_dir = Path(headers_dir).resolve()
    headers_so  = Path(headers_so).resolve()

    if not headers_so.exists():
        raise FileNotFoundError(
            f"Compiled headers not found: {headers_so}\n"
            "Run TFile::MakeProject on an experimental ROOT file, then "
            "compile with 'make -f MAKEP' to generate recovered_headers.so."
        )

    ROOT.gInterpreter.AddIncludePath(str(headers_dir))

    # Load TrackData.h if present; otherwise declare the struct inline
    hdr = headers_dir / "TrackData.h"
    if hdr.exists():
        ROOT.gInterpreter.ProcessLine(f'#include "{hdr}"')
    else:
        ROOT.gInterpreter.Declare("""
            #include <vector>
            class TrackData {
            public:
                int nPoints;
                float slope_xy, intercept_xy, slope_zy, intercept_zy;
                std::vector<int>   y;
                std::vector<float> x, z, charge, sigmas_x, sigmas_z, sigmas_y;
            };
        """)

    ROOT.gSystem.Load(str(headers_so))
    _HEADERS_LOADED = True


# ═══════════════════════════════════════════════════════════════════════════════
# Main track-based loader (PyROOT)
# ═══════════════════════════════════════════════════════════════════════════════

def load(
    filepath: str | Path,
    headers_dir: str | Path = "headers",
    headers_so: str | Path = "headers/recovered_headers.so",
    chi2_ndf_max: float = 25.0,
    min_track_points: int = 3,
    min_adc: float = 0.0,
    upper_percentile: float = 95.0,
) -> dict[str, Any]:
    """Load per-centroid dE/dx from the ``tracks`` branch via PyROOT.

    Applies chi²/ndf quality cuts and computes dE/dx as
    ADC[j] / (track_length / nPoints) for each centroid.

    Parameters
    ----------
    filepath : str or Path
        Experimental ROOT file with ``trackingData`` tree.
    headers_dir, headers_so : str or Path
        Location of compiled C++ headers.
    chi2_ndf_max : float
        Maximum chi²/ndf for track acceptance.
    min_track_points : int
        Minimum centroids per track.
    min_adc : float
        ADC noise floor — reject centroids below this.
    upper_percentile : float
        Clip the extreme high-ADC tail above this percentile.

    Returns
    -------
    dict
        dedx        : np.ndarray — per-centroid dE/dx [ADC/mm].
        tracks      : list[dict] — per-track metadata.
        n_tracks    : int
        n_events    : int
        source      : str
    """
    import ROOT

    filepath = Path(filepath)
    print(f"  Loading experimental: {filepath.name}")

    _ensure_headers(headers_dir, headers_so)

    tfile = ROOT.TFile.Open(str(filepath))
    if not tfile or tfile.IsZombie():
        raise IOError(f"Cannot open {filepath}")
    tree = tfile.Get("trackingData")
    if not tree:
        tfile.Close()
        raise KeyError("'trackingData' tree not found")

    tracks_vec = ROOT.std.vector("TrackData")()
    tree.SetBranchAddress("tracks", tracks_vec)

    all_dedx = []
    track_info = []

    n_events = tree.GetEntries()
    for i in range(n_events):
        tree.GetEntry(i)
        for track in tracks_vec:
            n = track.nPoints
            if n < min_track_points:
                continue

            # ── Chi²/ndf quality cut ────────────────────────────────────
            chi2 = 0.0
            for j in range(n):
                y_j = float(track.y[j])
                sx  = max(float(track.sigmas_x[j]), 1e-4)
                sz  = max(float(track.sigmas_z[j]), 1e-4)
                res_x = float(track.x[j]) - (track.slope_xy * y_j + track.intercept_xy)
                res_z = float(track.z[j]) - (track.slope_zy * y_j + track.intercept_zy)
                chi2 += (res_x / sx) ** 2 + (res_z / sz) ** 2
            ndf = max(2 * n - 4, 1)
            if chi2 / ndf > chi2_ndf_max:
                continue

            # ── 3D track length → step size per centroid ────────────────
            x0, x1 = float(track.x[0]), float(track.x[n - 1])
            y0, y1 = float(track.y[0]), float(track.y[n - 1])
            z0, z1 = float(track.z[0]), float(track.z[n - 1])
            length = np.sqrt((x1 - x0)**2 + (y1 - y0)**2 + (z1 - z0)**2)
            if length < 1e-3:
                continue
            dx = length / n

            # ── Per-centroid dE/dx ──────────────────────────────────────
            for j in range(n):
                adc = float(track.charge[j])
                if adc > min_adc:
                    all_dedx.append(adc / dx)

            track_info.append({
                "event":    i,
                "nPoints":  n,
                "chi2_ndf": chi2 / ndf,
                "length":   length,
                "x": [float(track.x[j]) for j in range(n)],
                "y": [float(track.y[j]) for j in range(n)],
                "z": [float(track.z[j]) for j in range(n)],
                "charge": [float(track.charge[j]) for j in range(n)],
            })

    tfile.Close()

    arr = np.array(all_dedx, dtype=np.float64)
    n_unique = len(np.unique(arr))
    print(f"    {len(arr):,} dE/dx values, {n_unique} unique, "
          f"{len(track_info)} tracks accepted")

    if n_unique < 20:
        raise RuntimeError(
            "Too few unique dE/dx values — charge fields may be corrupted. "
            "Check that TrackData.h matches the ROOT file being read."
        )

    # Clip extreme tail
    if upper_percentile < 100 and len(arr) > 0:
        upper = np.percentile(arr, upper_percentile)
        arr = arr[arr < upper]

    return {
        "dedx":     arr,
        "tracks":   track_info,
        "n_tracks": len(track_info),
        "n_events": n_events,
        "source":   filepath.stem,
    }


# ═══════════════════════════════════════════════════════════════════════════════
# Raw TPC pad-hit loader (uproot — no PyROOT needed)
# ═══════════════════════════════════════════════════════════════════════════════

def load_raw_hits(
    filepath: str | Path,
    tree_name: str = "trackingData",
    noise_sigma: float = 3.0,
    min_adc_signal: int = 10,
    entry_stop: Optional[int] = None,
) -> dict[str, Any]:
    """Load raw TPC pad hits for 3D displays and pad-plane diagnostics.

    Reads the ``tpc`` branch via uproot (fully split, no dictionary needed).

    Parameters
    ----------
    filepath : str or Path
        Experimental ROOT file.
    tree_name : str
        TTree name (default ``"trackingData"``).
    noise_sigma : float
        Reject hits where (val - pedestal) < noise_sigma × peddev.
    min_adc_signal : int
        Absolute minimum ADC above pedestal.
    entry_stop : int or None
        Load only the first N events (for quick inspection).

    Returns
    -------
    dict
        signal   : ak.Array — pedestal-subtracted ADC per hit.
        row      : ak.Array — pad row indices.
        col      : ak.Array — pad column indices.
        val      : ak.Array — raw ADC.
        pedestal : ak.Array
        peddev   : ak.Array or None
        n_events : int
        source   : str
    """
    import awkward as ak
    import uproot

    filepath = Path(filepath)
    print(f"  Loading raw hits: {filepath.name}")

    with uproot.open(str(filepath)) as f:
        tree = f[tree_name]

        val = tree["tpc/tpc.val"].array(library="ak", entry_stop=entry_stop)
        ped = tree["tpc/tpc.pedestal"].array(library="ak", entry_stop=entry_stop)
        row = tree["tpc/tpc.row"].array(library="ak", entry_stop=entry_stop)
        col = tree["tpc/tpc.column"].array(library="ak", entry_stop=entry_stop)

        # peddev is optional — some files don't have it
        try:
            dev = tree["tpc/tpc.peddev"].array(library="ak", entry_stop=entry_stop)
        except Exception:
            dev = None

    # Pedestal-subtracted signal
    sig = val - ped

    # Apply noise threshold
    if dev is not None:
        noise_floor = noise_sigma * dev
        good = sig > noise_floor
    else:
        good = sig > min_adc_signal

    return {
        "signal":   sig[good],
        "row":      row[good],
        "col":      col[good],
        "val_raw":  val,
        "pedestal": ped,
        "peddev":   dev,
        "noise_mask": good,
        "n_events": len(val),
        "source":   filepath.stem,
    }


# ═══════════════════════════════════════════════════════════════════════════════
# Branch inspector utility
# ═══════════════════════════════════════════════════════════════════════════════

def inspect_branches(filepath: str | Path, tree_name: str = "trackingData") -> None:
    """Print all branch names and types — run this to debug data issues."""
    import uproot

    filepath = Path(filepath)
    print(f"\nInspecting: {filepath.name}\n")

    with uproot.open(str(filepath)) as f:
        print("=== Top-level keys ===")
        for k in f.keys():
            print(f"  {k}")

        print(f"\n=== {tree_name} branches ===")
        try:
            tree = f[tree_name]
            for name, branch in tree.items():
                try:
                    dtype = branch.interpretation
                except Exception:
                    dtype = "?"
                print(f"  {name:<60}  {dtype}")
        except Exception as e:
            print(f"  Could not open {tree_name}: {e}")
