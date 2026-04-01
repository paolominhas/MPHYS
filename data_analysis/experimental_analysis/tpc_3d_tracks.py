"""
HIBEAM Prototype TPC — 3-D Track Visualisation
================================================
Reads the same experimental ROOT files as the pad-plane diagnostic
script (tree: trackingData, branch: tpc) and produces publication-
quality 3-D scatter plots of reconstructed TPC tracks.

Pad coordinates give (x, y) on the pad plane; the drift-time sample
gives the z-axis.  Colour encodes the pedestal-subtracted ADC signal,
which is proportional to the local ionisation energy deposit (dE/dx).

Figures are produced in both PDF (vector, for publication) and PNG
(raster, for slides / quick inspection).

Branch layout expected
----------------------
  tpc.val      — raw ADC value
  tpc.pedestal — electronic baseline per hit
  tpc.peddev   — pedestal RMS  (for noise threshold)
  tpc.row      — pad row   →  mapped to y
  tpc.column   — pad column →  mapped to x

If a time-sample branch exists (tpc.timebin or tpc.time), it is used
for the z-axis.  Otherwise the hit index within each event is used as
a proxy (suitable for single-track events where hits are time-ordered).

Usage
-----
    python tpc_3d_tracks.py                        # batch, all files
    python tpc_3d_tracks.py --file run042.root     # single file
    python tpc_3d_tracks.py --events 5,12,87       # specific events
"""

from __future__ import annotations
import argparse, glob, os, sys, textwrap
import numpy as np
import uproot
import awkward as ak

# ── Matplotlib setup — must precede any pyplot import ─────────────────────
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import cm, colors as mcolors
from mpl_toolkits.mplot3d import Axes3D          # noqa: F401  (registers projection)

# ── Publication style ─────────────────────────────────────────────────────
plt.rcParams.update({
    # Typography — use Computer Modern (TeX default) for consistency
    # with most HEP journals; falls back gracefully.
    "font.family":        "serif",
    "font.serif":         ["CMU Serif", "Computer Modern Roman",
                           "DejaVu Serif", "Times New Roman"],
    "mathtext.fontset":   "cm",
    "font.size":          10,
    "axes.titlesize":     11,
    "axes.labelsize":     10,
    "xtick.labelsize":    8,
    "ytick.labelsize":    8,
    "legend.fontsize":    8,
    # Lines / markers
    "lines.linewidth":    0.8,
    "lines.markersize":   3,
    # Figure defaults
    "figure.dpi":         150,
    "savefig.dpi":        300,
    "savefig.bbox":       "tight",
    "savefig.pad_inches": 0.05,
    # Axes
    "axes.linewidth":     0.6,
    "xtick.major.width":  0.6,
    "ytick.major.width":  0.6,
    "xtick.minor.width":  0.4,
    "ytick.minor.width":  0.4,
    "xtick.direction":    "in",
    "ytick.direction":    "in",
    "xtick.top":          True,
    "ytick.right":        True,
})


# ═══════════════════════════════════════════════════════════════════════════
#  Configuration
# ═══════════════════════════════════════════════════════════════════════════
INPUT_DIR       = "../experimental_data/"
OUTPUT_DIR      = "analysis_plots/"
NOISE_SIGMA     = 3.0      # noise rejection threshold (units of peddev)
MIN_HITS        = 8        # minimum hits for an event to be plotted
MAX_EVENTS_SHOW = 6        # how many events per multi-panel figure
PAD_PITCH_X     = 1.0      # mm — pad column pitch  (adjust to geometry)
PAD_PITCH_Y     = 1.0      # mm — pad row pitch
DRIFT_SPEED     = 1.0      # mm / time-bin  (adjust to detector config)
MARKER_SIZE     = 18        # scatter point size (pt²)
CMAP            = "inferno" # perceptually uniform, print-safe


# ═══════════════════════════════════════════════════════════════════════════
#  Data loading
# ═══════════════════════════════════════════════════════════════════════════
def load_events(path: str) -> list[dict]:
    """
    Load per-event TPC hits from a ROOT file.

    Returns a list of dicts, one per event that passes the noise cut
    and minimum-hit requirement.  Each dict contains numpy arrays:
        x, y, z, signal   (all in physical units where possible)
    """
    with uproot.open(path) as f:
        tree = f["trackingData"]
        keys = [k.split("/")[-1] for k in tree.keys()]

        val = tree["tpc/tpc.val"].array()
        ped = tree["tpc/tpc.pedestal"].array()
        dev = tree["tpc/tpc.peddev"].array()
        row = tree["tpc/tpc.row"].array()
        col = tree["tpc/tpc.column"].array()

        # Attempt to find a time-sample branch
        has_time = False
        for candidate in ("tpc.timebin", "tpc.time", "tpc.ts",
                          "tpc.timestamp", "tpc.sample"):
            full = f"tpc/{candidate}"
            if full in tree or candidate in keys:
                try:
                    tbin = tree[full].array()
                    has_time = True
                    break
                except Exception:
                    pass

    # Pedestal subtraction & noise mask
    sig  = val - ped
    good = sig > NOISE_SIGMA * dev

    events = []
    n_total = len(val)
    for i in range(n_total):
        mask_i = good[i]
        if ak.sum(mask_i) < MIN_HITS:
            continue

        s = ak.to_numpy(sig[i][mask_i]).astype(float)
        r = ak.to_numpy(row[i][mask_i]).astype(float)
        c = ak.to_numpy(col[i][mask_i]).astype(float)

        if has_time:
            t = ak.to_numpy(tbin[i][mask_i]).astype(float)
        else:
            # Fall back: use sequential index (adequate when hits are
            # stored in acquisition-time order within each event)
            t = np.arange(len(s), dtype=float)

        events.append({
            "event_id": i,
            "x":        c * PAD_PITCH_X,
            "y":        r * PAD_PITCH_Y,
            "z":        t * DRIFT_SPEED,
            "signal":   s,
        })

    print(f"  Events loaded  : {n_total}")
    print(f"  After cuts     : {len(events)}  "
          f"(≥ {MIN_HITS} hits, {NOISE_SIGMA}σ noise cut)")
    return events


# ═══════════════════════════════════════════════════════════════════════════
#  Single-event 3-D figure  (for highlight / supplemental material)
# ═══════════════════════════════════════════════════════════════════════════
def plot_single_event(evt: dict, filename: str, output_dir: str,
                      elev: float = 25, azim: float = -55) -> None:
    """
    One event → one figure with a single 3-D axis and a colour bar.
    """
    fig = plt.figure(figsize=(5.0, 4.5))     # single-column journal width
    ax  = fig.add_subplot(111, projection="3d")

    norm = mcolors.Normalize(vmin=evt["signal"].min(),
                             vmax=evt["signal"].max())
    mapper = cm.ScalarMappable(norm=norm, cmap=CMAP)

    sc = ax.scatter(evt["x"], evt["y"], evt["z"],
                    c=evt["signal"], cmap=CMAP, norm=norm,
                    s=MARKER_SIZE, edgecolors="none", alpha=0.85,
                    depthshade=True, rasterized=True)

    ax.set_xlabel(r"$x\;(\mathrm{column \times pitch})$ [mm]", labelpad=6)
    ax.set_ylabel(r"$y\;(\mathrm{row \times pitch})$ [mm]",    labelpad=6)
    ax.set_zlabel(r"$z\;(\mathrm{drift})$ [mm]",               labelpad=6)
    ax.view_init(elev=elev, azim=azim)

    # Colour bar
    cbar = fig.colorbar(mapper, ax=ax, shrink=0.55, pad=0.12, aspect=18)
    cbar.set_label(r"Signal  [ADC $-$ pedestal]", fontsize=9)
    cbar.ax.tick_params(labelsize=7)

    # Light grid styling
    ax.xaxis.pane.fill = False
    ax.yaxis.pane.fill = False
    ax.zaxis.pane.fill = False
    ax.xaxis.pane.set_edgecolor((0.85, 0.85, 0.85, 1.0))
    ax.yaxis.pane.set_edgecolor((0.85, 0.85, 0.85, 1.0))
    ax.zaxis.pane.set_edgecolor((0.85, 0.85, 0.85, 1.0))
    ax.grid(True, linewidth=0.3, alpha=0.5)

    eid = evt["event_id"]
    nhits = len(evt["signal"])
    ax.set_title(f"Event {eid}   ({nhits} hits)", fontsize=10, pad=10)

    stem = os.path.basename(filename).replace(".root", "")
    for ext in ("pdf", "png"):
        out = os.path.join(output_dir,
                           f"{stem}_event{eid:04d}_3d.{ext}")
        fig.savefig(out, dpi=300)
    plt.close(fig)
    print(f"    → Event {eid} saved")


# ═══════════════════════════════════════════════════════════════════════════
#  Multi-event panel figure  (main publication figure)
# ═══════════════════════════════════════════════════════════════════════════
def plot_multi_event(events: list[dict], filename: str,
                     output_dir: str, n_show: int = MAX_EVENTS_SHOW,
                     elev: float = 25, azim: float = -55) -> None:
    """
    Produce a grid of 3-D track views sharing a common colour scale,
    suitable for a full-width journal figure.
    """
    if len(events) == 0:
        return

    selected = events[:n_show]
    n = len(selected)
    ncols = min(n, 3)
    nrows = int(np.ceil(n / ncols))

    fig_w = 3.4 * ncols + 0.8          # ~3.4 in per panel + colorbar
    fig_h = 3.2 * nrows + 0.6
    fig = plt.figure(figsize=(fig_w, fig_h))

    # Global colour normalisation across all shown events
    all_sig = np.concatenate([e["signal"] for e in selected])
    vmin, vmax = np.percentile(all_sig, [1, 99])
    norm   = mcolors.Normalize(vmin=vmin, vmax=vmax)
    mapper = cm.ScalarMappable(norm=norm, cmap=CMAP)
    mapper.set_array([])

    axes = []
    for idx, evt in enumerate(selected):
        ax = fig.add_subplot(nrows, ncols, idx + 1, projection="3d")
        axes.append(ax)

        ax.scatter(evt["x"], evt["y"], evt["z"],
                   c=evt["signal"], cmap=CMAP, norm=norm,
                   s=MARKER_SIZE * 0.7, edgecolors="none",
                   alpha=0.85, depthshade=True, rasterized=True)

        ax.set_xlabel("$x$ [mm]", labelpad=4, fontsize=8)
        ax.set_ylabel("$y$ [mm]", labelpad=4, fontsize=8)
        ax.set_zlabel("$z$ [mm]", labelpad=4, fontsize=8)
        ax.tick_params(labelsize=6, pad=1)
        ax.view_init(elev=elev, azim=azim)

        ax.xaxis.pane.fill = False
        ax.yaxis.pane.fill = False
        ax.zaxis.pane.fill = False
        ax.xaxis.pane.set_edgecolor((0.85, 0.85, 0.85, 1.0))
        ax.yaxis.pane.set_edgecolor((0.85, 0.85, 0.85, 1.0))
        ax.zaxis.pane.set_edgecolor((0.85, 0.85, 0.85, 1.0))
        ax.grid(True, linewidth=0.25, alpha=0.4)

        eid   = evt["event_id"]
        nhits = len(evt["signal"])
        ax.set_title(f"Event {eid}  ({nhits} hits)",
                     fontsize=8, pad=4)

    # Shared colour bar
    fig.subplots_adjust(left=0.04, right=0.88,
                        bottom=0.06, top=0.91,
                        wspace=0.22, hspace=0.30)
    cax = fig.add_axes([0.91, 0.15, 0.015, 0.65])
    cbar = fig.colorbar(mapper, cax=cax)
    cbar.set_label(r"Signal  [ADC $-$ pedestal]", fontsize=9)
    cbar.ax.tick_params(labelsize=7)

    stem = os.path.basename(filename).replace(".root", "")
    fig.suptitle(
        f"Reconstructed TPC tracks — {stem}\n"
        rf"Noise cut: {NOISE_SIGMA}$\sigma$  |  "
        f"Min. hits: {MIN_HITS}  |  "
        f"Colour: pedestal-subtracted ADC",
        fontsize=10, y=0.97)

    for ext in ("pdf", "png"):
        out = os.path.join(output_dir, f"{stem}_tracks_3d.{ext}")
        fig.savefig(out, dpi=300)
    plt.close(fig)
    print(f"  ✔  Multi-panel saved → {stem}_tracks_3d.{{pdf,png}}")


# ═══════════════════════════════════════════════════════════════════════════
#  2-D projections panel  (complements the 3-D view)
# ═══════════════════════════════════════════════════════════════════════════
def plot_projections(evt: dict, filename: str, output_dir: str) -> None:
    """
    Three orthogonal projections (xy, xz, yz) of a single event.
    Useful where 3-D perspective may obscure structure.
    """
    fig, axes = plt.subplots(1, 3, figsize=(7.0, 2.6))

    pairs  = [("x", "y"), ("x", "z"), ("y", "z")]
    labels = [("$x$ [mm]", "$y$ [mm]"),
              ("$x$ [mm]", "$z$ [mm]"),
              ("$y$ [mm]", "$z$ [mm]")]

    norm = mcolors.Normalize(vmin=evt["signal"].min(),
                             vmax=evt["signal"].max())

    for ax, (kx, ky), (lx, ly) in zip(axes, pairs, labels):
        sc = ax.scatter(evt[kx], evt[ky], c=evt["signal"],
                        cmap=CMAP, norm=norm, s=MARKER_SIZE * 0.6,
                        edgecolors="none", alpha=0.85, rasterized=True)
        ax.set_xlabel(lx)
        ax.set_ylabel(ly)
        ax.set_aspect("auto")
        ax.grid(True, linewidth=0.25, alpha=0.3)

    cbar = fig.colorbar(sc, ax=axes.tolist(), shrink=0.85,
                        pad=0.03, aspect=25)
    cbar.set_label(r"Signal  [ADC $-$ pedestal]", fontsize=8)
    cbar.ax.tick_params(labelsize=7)

    eid = evt["event_id"]
    fig.suptitle(f"Event {eid} — orthogonal projections", fontsize=10,
                 y=1.02)

    stem = os.path.basename(filename).replace(".root", "")
    for ext in ("pdf", "png"):
        out = os.path.join(output_dir,
                           f"{stem}_event{eid:04d}_proj.{ext}")
        fig.savefig(out, dpi=300)
    plt.close(fig)


# ═══════════════════════════════════════════════════════════════════════════
#  Driver
# ═══════════════════════════════════════════════════════════════════════════
def analyse_file(path: str, output_dir: str,
                 event_indices: list[int] | None = None) -> None:
    filename = os.path.basename(path)
    print(f"\n{'─' * 64}\n  {filename}")
    events = load_events(path)
    if not events:
        print("  ⚠  No events survived cuts.")
        return

    # Sort by number of hits (descending) to show the richest tracks first
    events.sort(key=lambda e: len(e["signal"]), reverse=True)

    # If specific events were requested, filter
    if event_indices is not None:
        id_set = set(event_indices)
        events = [e for e in events if e["event_id"] in id_set]
        if not events:
            print(f"  ⚠  None of the requested events survived cuts.")
            return

    # Multi-panel overview
    plot_multi_event(events, path, output_dir)

    # Individual figures + projections for the top events
    for evt in events[:min(3, len(events))]:
        plot_single_event(evt, path, output_dir)
        plot_projections(evt, path, output_dir)


def main():
    parser = argparse.ArgumentParser(
        description="3-D TPC track visualisation for HIBEAM data.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=textwrap.dedent("""\
            examples:
              python tpc_3d_tracks.py
              python tpc_3d_tracks.py --file run042.root
              python tpc_3d_tracks.py --events 5,12,87
        """))
    parser.add_argument("--file", type=str, default=None,
                        help="Single ROOT file to process.")
    parser.add_argument("--events", type=str, default=None,
                        help="Comma-separated event indices to plot.")
    parser.add_argument("--input-dir", type=str, default=INPUT_DIR)
    parser.add_argument("--output-dir", type=str, default=OUTPUT_DIR)
    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)

    event_indices = None
    if args.events:
        event_indices = [int(x) for x in args.events.split(",")]

    if args.file:
        analyse_file(args.file, args.output_dir, event_indices)
    else:
        files = sorted(glob.glob(os.path.join(args.input_dir, "*.root")))
        if not files:
            print(f"No .root files found in {args.input_dir}")
            sys.exit(1)
        print(f"Found {len(files)} file(s)")
        for f in files:
            try:
                analyse_file(f, args.output_dir, event_indices)
            except Exception as exc:
                print(f"  ✗  {type(exc).__name__}: {exc}")

    print(f"\n{'═' * 64}\n  Done.")


if __name__ == "__main__":
    main()