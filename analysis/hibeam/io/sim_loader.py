"""
hibeam.io.sim_loader — Load Geant4 simulation ROOT files
==========================================================
Reads processed HIBEAM ROOT ntuples produced by the hibeam_ana
analysis framework.  Uses uproot exclusively (no PyROOT dependency).

Consolidates logic previously scattered across:
  - segmentation/final.py          (load_tpc_data)
  - combined_analysis/combining_dedx.py  (load_sim)
  - simulation_analysis/prototpc_final.py
  - simulation_analysis/pid_deltae_e.py  (load_sim)
  - simulation_analysis/sim_prototpc_dedx.py

Every function returns a standardised dictionary so downstream code
never needs to know which file format it came from.
"""

from __future__ import annotations

import warnings
from pathlib import Path
from typing import Any, Optional

import awkward as ak
import numpy as np
import uproot

from hibeam.utils import find_branch


# ═══════════════════════════════════════════════════════════════════════════════
# Main TPC loader (for HIBEAMScatter.root and similar full-TPC files)
# ═══════════════════════════════════════════════════════════════════════════════

def load(
    filepath: str | Path,
    tree_name: str = "hibeam",
    branch_prefix: str = "TPC",
) -> dict[str, Any]:
    """Load TPC branch data from a processed HIBEAM ROOT ntuple.

    Branch discovery is dynamic: minor naming differences between ROOT
    split-class formats (``TPC/Edep`` vs ``TPC.Edep``) are handled
    automatically.

    Parameters
    ----------
    filepath : str or Path
        Path to ROOT file produced by hibeam_ana.
    tree_name : str
        TTree name inside the file (default ``"hibeam"``).
    branch_prefix : str
        Branch family to load.  ``"TPC"`` for the main TPC,
        ``"ProtoTPC"`` for the prototype.

    Returns
    -------
    dict
        edep       : ak.Array — per-event Edep values.
        pad_row    : ak.Array or None
        pad_col    : ak.Array or None
        timestamp  : ak.Array or None
        n_el       : ak.Array or None
        weights    : np.ndarray — PrimaryWeight (1.0 if absent).
        n_events   : int
        source     : str — filename stem for labelling.
        units      : str — "GeV" (raw Geant4) or "MeV" if auto-converted.
    """
    filepath = Path(filepath)
    print(f"  Loading simulation: {filepath.name}")

    with uproot.open(str(filepath)) as root_file:
        # ── Locate TTree ────────────────────────────────────────────────
        if tree_name not in root_file:
            available = list(root_file.keys())
            raise KeyError(
                f"TTree '{tree_name}' not found in {filepath}.\n"
                f"Available: {available}"
            )
        tree = root_file[tree_name]
        all_keys = tree.keys(recursive=True)

        # ── Discover branches belonging to the prefix ───────────────────
        if branch_prefix == "TPC":
            # Main TPC: include "TPC" keys but exclude "ProtoTPC" and "target"
            prefix_keys = [
                k for k in all_keys
                if "TPC" in k
                and not k.startswith("ProtoTPC")
                and not k.startswith("target")
            ]
        else:
            prefix_keys = [k for k in all_keys if k.startswith(branch_prefix)]

        if not prefix_keys:
            raise KeyError(
                f"No '{branch_prefix}' branches found.\n"
                f"Available: {all_keys}"
            )

        def _find(suffix: str) -> Optional[str]:
            """Find full branch key ending in suffix."""
            return find_branch(prefix_keys, [suffix])

        def _safe_array(suffix: str) -> Optional[ak.Array]:
            key = _find(suffix)
            if key is None:
                return None
            try:
                return tree[key].array(library="ak")
            except Exception as exc:
                warnings.warn(f"Could not read '{key}': {exc}", stacklevel=3)
                return None

        # ── Required: Edep ──────────────────────────────────────────────
        edep_key = _find("Edep")
        if edep_key is None:
            raise KeyError(
                f"'{branch_prefix}/Edep' not found.\n"
                f"Branches: {prefix_keys}"
            )
        edep = tree[edep_key].array(library="ak")

        # ── Optional: pad positions, electron counts ────────────────────
        pad_row   = _safe_array("padRow")
        pad_col   = _safe_array("padColumn")
        timestamp = _safe_array("timestamp")
        n_el      = _safe_array("nEl")
        nhits_arr = _safe_array("nHits")

        # ── Optional: PrimaryWeight (for weighted histograms) ───────────
        weights_key = find_branch(all_keys, ["PrimaryWeight"])
        if weights_key:
            w_raw = tree[weights_key].array(library="ak")
            weights = ak.to_numpy(ak.fill_none(ak.firsts(w_raw), 1.0))
        else:
            weights = np.ones(len(edep), dtype=np.float64)

    return {
        "edep":      edep,
        "pad_row":   pad_row,
        "pad_col":   pad_col,
        "timestamp": timestamp,
        "n_el":      n_el,
        "nhits":     nhits_arr,
        "weights":   weights,
        "n_events":  len(edep),
        "source":    filepath.stem,
        "units":     "GeV",
    }


# ═══════════════════════════════════════════════════════════════════════════════
# ProtoTPC convenience loader
# ═══════════════════════════════════════════════════════════════════════════════

def load_prototpc(
    filepath: str | Path,
    tree_name: str = "hibeam",
) -> dict[str, Any]:
    """Load ProtoTPC branches (used by Krakow/Muon setups).

    Thin wrapper around :func:`load` with ``branch_prefix="ProtoTPC"``.
    """
    return load(filepath, tree_name=tree_name, branch_prefix="ProtoTPC")


# ═══════════════════════════════════════════════════════════════════════════════
# PID loader (TPC + calorimeter for ΔE-E plots)
# ═══════════════════════════════════════════════════════════════════════════════

def load_pid(
    filepath: str | Path,
    tree_name: str = "hibeam",
) -> dict[str, Any]:
    """Load both TPC and calorimeter deposits for ΔE-E PID analysis.

    Returns
    -------
    dict
        tpc_mev   : np.ndarray — total TPC deposit per event [MeV].
        cal_mev   : np.ndarray — total calorimeter deposit per event [MeV].
        n_tpc     : np.ndarray — TPC hits per event.
        n_cal     : np.ndarray — calorimeter hits per event.
        n_events  : int
        source    : str
    """
    filepath = Path(filepath)
    print(f"  Loading PID data: {filepath.name}")

    with uproot.open(str(filepath)) as f:
        tree = f[tree_name]
        all_keys = tree.keys(recursive=True)

        # ── TPC deposit ─────────────────────────────────────────────────
        tpc_key = find_branch(all_keys, [
            "ProtoTPC/Edep", "ProtoTPC.Edep", "TPC/Edep", "TPC.Edep",
        ])
        if tpc_key is None:
            raise KeyError(f"No TPC Edep branch found in {filepath}")
        tpc_edep = tree[tpc_key].array(library="ak")
        n_tpc    = ak.to_numpy(ak.num(tpc_edep)).astype(int)
        tpc_raw  = ak.to_numpy(ak.sum(tpc_edep, axis=1)).astype(float)

        # ── Calorimeter deposit ─────────────────────────────────────────
        cal_key = find_branch(all_keys, [
            "HRD/Edep", "HRD.Edep",
            "ScintHit/Edep", "ScintHit.Edep", "Scint/Edep",
        ])
        if cal_key is None:
            raise KeyError(f"No calorimeter Edep branch found in {filepath}")
        cal_edep = tree[cal_key].array(library="ak")
        n_cal    = ak.to_numpy(ak.num(cal_edep)).astype(int)
        cal_raw  = ak.to_numpy(ak.sum(cal_edep, axis=1)).astype(float)

    # ── Auto-detect units ───────────────────────────────────────────────
    def _to_mev(arr: np.ndarray, label: str) -> np.ndarray:
        mx = float(arr[arr > 0].max()) if (arr > 0).any() else 0.0
        if mx > 500:
            return arr       # already MeV
        return arr * 1000.0  # GeV → MeV

    return {
        "tpc_mev":  _to_mev(tpc_raw, "TPC"),
        "cal_mev":  _to_mev(cal_raw, "Cal"),
        "n_tpc":    n_tpc,
        "n_cal":    n_cal,
        "n_events": len(tpc_raw),
        "source":   filepath.stem,
    }


# ═══════════════════════════════════════════════════════════════════════════════
# Helpers for downstream code
# ═══════════════════════════════════════════════════════════════════════════════

def edep_to_mev(data: dict) -> dict:
    """Convert Edep arrays from GeV to MeV in-place if needed.

    Checks the ``units`` key.  If already ``"MeV"``, does nothing.
    Returns the same dict for chaining.
    """
    if data.get("units") == "MeV":
        return data
    # Multiply the jagged Edep array by 1000
    data["edep"] = data["edep"] * 1000.0
    data["units"] = "MeV"
    return data
