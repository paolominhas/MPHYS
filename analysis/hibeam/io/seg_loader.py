"""
hibeam.io.seg_loader — Load segmentation study ROOT files
==========================================================
Reads ``ana_geom_nSec*.root`` files produced by the segmentation
analysis chain.  Each file contains a TTree with ProtoTPC energy
deposits at a different detector segmentation (nSections).

Consolidates logic previously in:
  - segmentation/overlay.py
  - segmentation/segmentation_overlay.py
  - segmentation/nhits.py, nhits_pub.py, nhits_vs_nsections.py
  - segmentation/analyse_segments.py
  - segmentation/segment.py, simple.py
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

import awkward as ak
import numpy as np
import uproot

from hibeam.utils import extract_nsec, find_root_files


# ═══════════════════════════════════════════════════════════════════════════════
# Per-file record
# ═══════════════════════════════════════════════════════════════════════════════

def load_single(
    filepath: str | Path,
    tree_name: str = "hibeam",
    edep_branch: str = "ProtoTPC/Edep",
    nhits_branch: str = "ProtoTPC/nHits",
    min_steps: int = 1,
    edep_floor_mev: float = 0.0,
) -> dict[str, Any]:
    """Load energy deposit and hit data from a single segmentation file.

    Parameters
    ----------
    filepath : str or Path
        Path to ``ana_geom_nSec*.root``.
    tree_name : str
        TTree name.
    edep_branch : str
        Per-hit Edep branch (jagged).
    nhits_branch : str
        Per-event hit count branch.
    min_steps : int
        Minimum G4 steps to accept an event.  Set to 1 to keep
        zero-hit events visible in nHits distributions.
    edep_floor_mev : float
        Hard lower-bound cut on per-hit Edep [MeV].  Removes the
        sub-threshold plateau that biases fits for nSec > 12.

    Returns
    -------
    dict with keys: nsec, edep_flat, edep_per_event, nhits, n_total,
    n_zero, zero_frac, mean_nhits_all, mean_nhits_nonzero, source.
    """
    filepath = Path(filepath)
    nsec = extract_nsec(filepath.name)

    with uproot.open(str(filepath)) as rf:
        tree = rf[tree_name]
        all_keys = tree.keys(recursive=True)

        # ── Edep (jagged) ───────────────────────────────────────────────
        edep_key = None
        for candidate in [edep_branch, edep_branch.replace("/", ".")]:
            if candidate in all_keys:
                edep_key = candidate
                break
        # Fallback: search by suffix
        if edep_key is None:
            for k in all_keys:
                if k.endswith("Edep") and "ProtoTPC" in k:
                    edep_key = k
                    break
        if edep_key is None:
            raise KeyError(
                f"Edep branch not found in {filepath.name}.\n"
                f"Available: {all_keys}"
            )

        edep_jagged = tree[edep_key].array(library="ak")
        edep_flat   = ak.to_numpy(ak.flatten(edep_jagged))

        # Per-event sum
        edep_sum = ak.to_numpy(ak.sum(edep_jagged, axis=1)).astype(float)

        # Step count per event
        n_steps = ak.to_numpy(ak.num(edep_jagged)).astype(int)

        # ── nHits (scalar per event) ────────────────────────────────────
        nhits_key = None
        for candidate in [nhits_branch, nhits_branch.replace("/", ".")]:
            if candidate in all_keys:
                nhits_key = candidate
                break
        if nhits_key:
            nhits = ak.to_numpy(tree[nhits_key].array(library="ak"))
        else:
            nhits = n_steps  # fallback: use step count

    # ── Apply cuts ──────────────────────────────────────────────────────
    # Per-hit floor cut
    if edep_floor_mev > 0:
        edep_flat = edep_flat[edep_flat > edep_floor_mev / 1000.0]  # GeV

    # Per-event step cut for Edep sums
    event_mask = n_steps >= min_steps
    edep_sum_cut = edep_sum[event_mask]

    # nHits statistics
    n_total = len(nhits)
    n_zero  = int((nhits == 0).sum())
    nz      = nhits[nhits > 0]

    return {
        "nsec":               nsec,
        "edep_flat":          edep_flat,
        "edep_per_event":     edep_sum_cut,
        "nhits":              nhits,
        "n_steps":            n_steps,
        "n_total":            n_total,
        "n_zero":             n_zero,
        "zero_frac":          n_zero / max(n_total, 1),
        "mean_nhits_all":     float(nhits.mean()) if n_total > 0 else 0.0,
        "mean_nhits_nonzero": float(nz.mean()) if len(nz) > 0 else 0.0,
        "mode_nhits":         int(np.bincount(nz).argmax()) if len(nz) > 0 else 0,
        "source":             filepath.stem,
    }


# ═══════════════════════════════════════════════════════════════════════════════
# Batch loader for an entire segmentation directory
# ═══════════════════════════════════════════════════════════════════════════════

def load(
    directory: str | Path,
    pattern: str = "ana_geom_*.root",
    tree_name: str = "hibeam",
    edep_branch: str = "ProtoTPC/Edep",
    nhits_branch: str = "ProtoTPC/nHits",
    min_steps: int = 1,
    edep_floor_mev: float = 0.0,
) -> list[dict[str, Any]]:
    """Load all segmentation files in a directory, sorted by nSections.

    Parameters
    ----------
    directory : str or Path
        Directory containing ``ana_geom_nSec*.root`` files.
    pattern : str
        Glob pattern for input files.
    (remaining parameters passed through to :func:`load_single`)

    Returns
    -------
    list[dict]
        One record per file, sorted ascending by ``nsec``.
        Each record has the keys documented in :func:`load_single`.
    """
    files = find_root_files(directory, pattern)
    files.sort(key=lambda p: extract_nsec(p.name) or 0)

    records = []
    for fpath in files:
        nsec = extract_nsec(fpath.name)
        if nsec is None:
            continue
        try:
            rec = load_single(
                fpath,
                tree_name=tree_name,
                edep_branch=edep_branch,
                nhits_branch=nhits_branch,
                min_steps=min_steps,
                edep_floor_mev=edep_floor_mev,
            )
            records.append(rec)
            print(f"    nSec={nsec:3d}: {rec['n_total']:,} events, "
                  f"{rec['n_zero']:,} zero-hit ({rec['zero_frac']:.1%}), "
                  f"mean nHits={rec['mean_nhits_nonzero']:.1f}")
        except Exception as e:
            print(f"    [SKIP] nSec={nsec}: {e}")

    if not records:
        raise RuntimeError(f"No valid segmentation files in {directory}")

    return records
