#!/usr/bin/env python3
"""
nhits_publication.py
====================
Publication-quality nHits-per-event figure for journal submission.
A4 width, no legend, clean axes, PDF + PNG output.
"""

import argparse
import re
import sys
from pathlib import Path

import awkward as ak
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np
import uproot

parser = argparse.ArgumentParser()
parser.add_argument("--indir",   default="./ana_output_krakow")
parser.add_argument("--outdir",  default="./plots")
parser.add_argument("--pattern", default="ana_geom_*.root")
args = parser.parse_args()

indir  = Path(args.indir)
outdir = Path(args.outdir)
outdir.mkdir(parents=True, exist_ok=True)

def extract_nsec(fname):
    m = re.search(r'[nN][sS]ec0*(\d+)', fname)
    return int(m.group(1)) if m else 0

files = sorted(indir.glob(args.pattern), key=lambda p: extract_nsec(p.name))
if not files:
    sys.exit(f"ERROR: No files in {indir}")

records = []
for fpath in files:
    nsec = extract_nsec(fpath.name)
    try:
        with uproot.open(str(fpath)) as rf:
            nhits = ak.to_numpy(rf["hibeam"]["ProtoTPC/nHits"].array(library="ak"))
    except Exception as e:
        print(f"  [SKIP] nSec={nsec}: {e}")
        continue
    n_total = len(nhits)
    n_zero  = int((nhits == 0).sum())
    records.append({
        "nsec"     : nsec,
        "nhits"    : nhits,
        "n_total"  : n_total,
        "n_zero"   : n_zero,
        "zero_frac": n_zero / n_total,
        "mean_all" : float(nhits.mean()),
        "mean_nz"  : float(nhits[nhits > 0].mean()) if n_total > n_zero else 0.0,
    })

if not records:
    sys.exit("No data loaded.")

n_panels = len(records)
n_cols   = 4
n_rows   = (n_panels + n_cols - 1) // n_cols

# Shared x-limit
all_nz = np.concatenate([r["nhits"][r["nhits"] > 0] for r in records])
x_max  = int(np.percentile(all_nz, 99)) + 2

# Colours — simple, high contrast, colourblind-safe
C_ZERO = "#cf150b"   # red    — zero-hit bin
C_HIT  = "#1b20bf"   # blue   — hit bins
C_MEAN = "#010101"   # mean lines

# A4 width = 210 mm = 8.27 in
# Keep height proportional — 3 rows needs ~5.5 in
A4_W   = 8.27
PANEL_H = 1.9    # inches per row — keep panels compact
FIG_H   = PANEL_H * n_rows + 0.3  # small buffer

plt.rcParams.update({
    "font.family"      : "serif",
    "font.size"        : 9,
    "axes.titlesize"   : 10,
    "axes.labelsize"   : 9,
    "xtick.labelsize"  : 8,
    "ytick.labelsize"  : 8,
    "axes.linewidth"   : 0.7,
    "xtick.major.width": 0.7,
    "ytick.major.width": 0.7,
    "xtick.major.size" : 3.0,
    "ytick.major.size" : 3.0,
    "xtick.minor.size" : 1.5,
    "ytick.minor.size" : 1.5,
    "xtick.direction"  : "in",
    "ytick.direction"  : "in",
    "xtick.top"        : True,
    "ytick.right"      : True,
})

fig, axes_grid = plt.subplots(n_rows, n_cols,
                               figsize=(A4_W, FIG_H),
                               squeeze=False)

for idx, rec in enumerate(records):
    row = idx // n_cols
    col = idx  % n_cols
    ax  = axes_grid[row][col]

    nhits  = rec["nhits"]
    nsec   = rec["nsec"]

    bins   = np.arange(0, x_max + 2) - 0.5
    counts, edges = np.histogram(nhits[nhits <= x_max], bins=bins)
    bc = 0.5 * (edges[:-1] + edges[1:])
    nz = bc > 0.5

    # Non-zero bins
    ax.bar(bc[nz],  counts[nz],  width=0.85,
           color=C_HIT,  alpha=0.80, linewidth=0)
    # Zero bin
    ax.bar(bc[~nz], counts[~nz], width=0.85,
           color=C_ZERO, alpha=0.90, linewidth=0)

    # Mean lines — thin, unobtrusive
    ax.axvline(rec["mean_nz"],  color=C_MEAN, lw=0.9, ls="--")
    ax.axvline(rec["mean_all"], color=C_ZERO, lw=0.9, ls=":")

    # Panel label: nSec value top-right inside panel
    ax.text(0.97, 0.96, f"$N_{{\\rm seg}}={nsec}$",
            transform=ax.transAxes, ha="right", va="top",
            fontsize=11, color="black",
            bbox=dict(boxstyle="round,pad=0.2", facecolor="white", edgecolor="none", alpha=0.7))

    ax.set_yscale("log")
    ax.set_xlim(-1, x_max)
    ax.set_ylim(0.5, rec["n_total"] * 2.5)

    # Minimal ticks — 3 on x, log-spaced on y
    ax.xaxis.set_major_locator(ticker.MaxNLocator(integer=True, nbins=3))
    ax.yaxis.set_major_locator(ticker.LogLocator(base=10, numticks=4))
    ax.yaxis.set_major_formatter(ticker.LogFormatter(base=10, labelOnlyBase=True))
    ax.yaxis.set_minor_locator(ticker.NullLocator())

    # Axis labels only on outer panels
    is_left   = (col == 0)
    is_bottom = (row == n_rows - 1) or (idx >= n_panels - n_cols)

    if is_left:
        ax.set_ylabel("Events", fontsize=11, labelpad=2)
    else:
        ax.set_yticklabels([])

    if is_bottom:
        ax.set_xlabel("$n_{\\rm hits}$ / event", fontsize=11, labelpad=2)
    else:
        ax.set_xticklabels([])

# Hide unused axes
for j in range(n_panels, n_rows * n_cols):
    axes_grid[j // n_cols][j % n_cols].set_visible(False)

fig.tight_layout(pad=0.4, h_pad=0.5, w_pad=0.4)

# Save both PDF (vector) and PNG (raster)
for ext in [".pdf", ".png"]:
    path = outdir / f"nhits_segmentation{ext}"
    dpi  = 300 if ext == ".png" else None
    fig.savefig(path, dpi=dpi, bbox_inches="tight")
    print(f"Saved: {path}")