#!/usr/bin/env python3
"""
event_display_3d.py
===================
Quick 3D visualisation of TPC pad hits for a single HIBEAM event.
Pad hits are coloured by ADC value (green=low, red=high energy loss).

Uses load_tpc_data from tpc_dedx_analysis.py for data loading.

Usage
-----
    python event_display_3d.py <file.root> [event_number]

    event_number defaults to 0.
"""

import sys
import numpy as np
import awkward as ak
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors as mcolors

from final import load_tpc_data


def get_event_hits(data: dict, event: int) -> dict:
    """
    Extract pad hit coordinates and charges for one event.

    Returns dict with arrays: row, col, timestamp, edep, n_hits.
    Raises IndexError if event number is out of range.
    """
    n = data["n_events"]
    if event >= n or event < 0:
        raise IndexError(f"Event {event} out of range — file has {n} events.")

    pad_row = data.get("pad_row")
    pad_col = data.get("pad_col")
    ts      = data.get("timestamp")
    edep    = data["edep"]

    if any(x is None for x in [pad_row, pad_col, ts]):
        raise RuntimeError(
            "pad_row / pad_col / timestamp branches not found.\n"
            "These are produced by assign_pads() in the C++ filter — "
            "make sure you're using a processed file."
        )

    row  = ak.to_numpy(pad_row[event]).astype(float)
    col  = ak.to_numpy(pad_col[event]).astype(float)
    time = ak.to_numpy(ts[event]).astype(float)
    ev_edep = ak.to_numpy(edep[event]).astype(float)

    # Edep is per G4-step; broadcast to pad length if sizes differ.
    # The number of pads >> number of steps, so repeat each step's Edep
    # across its assigned pads (uniform approximation).
    if len(ev_edep) != len(row):
        # Repeat Edep values to match pad count
        repeats = int(np.ceil(len(row) / max(len(ev_edep), 1)))
        ev_edep = np.tile(ev_edep, repeats)[:len(row)]

    return {
        "row"    : row,
        "col"    : col,
        "time"   : time,
        "edep"   : ev_edep,
        "n_hits" : len(row),
    }


def plot_event_3d(hits: dict, event: int, output_path: str = None) -> plt.Figure:
    """
    3D scatter plot of pad hits coloured by energy deposition.

    Axes:
      x = pad column  (transverse)
      y = pad row     (transverse)
      z = timestamp   (drift / beam direction)

    Colour scale: green (low Edep) → yellow → red (high Edep).
    This is the standard 'hot track' convention used in ALICE and STAR
    event displays.
    """
    row   = hits["row"]
    col   = hits["col"]
    time  = hits["time"]
    edep  = hits["edep"]

    if hits["n_hits"] == 0:
        print(f"Event {event} has no pad hits — nothing to plot.")
        return None

    # Normalise Edep to [0, 1] for colour mapping
    e_min, e_max = edep.min(), edep.max()
    if e_max > e_min:
        e_norm = (edep - e_min) / (e_max - e_min)
    else:
        e_norm = np.zeros_like(edep)

    cmap   = cm.RdYlGn_r          # green=low, yellow=mid, red=high
    colours = cmap(e_norm)

    # Marker size proportional to charge (capped for clarity)
    sizes = 10.0 + 40.0 * e_norm

    fig = plt.figure(figsize=(10, 7))
    ax  = fig.add_subplot(111, projection="3d")

    sc = ax.scatter(
        col, row, time,
        c=e_norm, cmap=cmap,
        s=sizes, alpha=0.75,
        edgecolors="none",
        depthshade=True,
    )

    # Colour bar
    cbar = fig.colorbar(sc, ax=ax, pad=0.10, shrink=0.6)
    cbar.set_label(r"Normalised $E_\mathrm{dep}$  (green=low, red=high)",
                   fontsize=10)

    ax.set_xlabel("Pad column  (transverse)", fontsize=10, labelpad=8)
    ax.set_ylabel("Pad row  (transverse)",    fontsize=10, labelpad=8)
    ax.set_zlabel("Timestamp  (drift)",       fontsize=10, labelpad=8)
    ax.set_title(
        f"HIBEAM TPC Event Display  —  Event {event}\n"
        f"{hits['n_hits']} pad hits   "
        rf"$\Sigma E_\mathrm{{dep}}$ = {edep.sum():.4g} MeV",
        fontsize=11,
    )

    if output_path:
        fig.savefig(output_path, dpi=150, bbox_inches="tight")
        print(f"Figure saved: {output_path}")
    else:
        plt.show()

    return fig


def main(filepath: str, event: int = 0) -> None:
    print(f"\nLoading: {filepath}")
    data = load_tpc_data(filepath)
    print(f"  {data['n_events']:,} events in file")

    print(f"  Extracting event {event}...")
    hits = get_event_hits(data, event)
    print(f"  {hits['n_hits']} pad hits,  "
          f"sum(Edep) = {hits['edep'].sum():.4g} MeV")

    out = f"event_{event}_display.png"
    plot_event_3d(hits, event, output_path=out)


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python event_display_3d.py <file.root> [event_number]")
        sys.exit(1)
    filepath = sys.argv[1]
    event    = int(sys.argv[2]) if len(sys.argv) > 2 else 0
    main(filepath, event)