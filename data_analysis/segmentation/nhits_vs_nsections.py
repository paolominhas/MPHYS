#!/usr/bin/env python3
"""
nhits_vs_nsections.py
=====================
Two-panel publication figure: mean/mode nHits vs segmentation
for muon lab (left) and Krakow scattering (right).

Three lines per panel:
  - Mean nHits including zero-hit events
  - Mean nHits excluding zero-hit events
  - Mode (most frequent nHits) excluding zero-hit events

Single shared legend below both panels.
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
parser.add_argument("--muon_dir",   default="./ana_output_muon")
parser.add_argument("--krakow_dir", default="./ana_output_krakow")
parser.add_argument("--outdir",     default="./plots")
parser.add_argument("--pattern",    default="ana_geom_*.root")
args = parser.parse_args()

outdir = Path(args.outdir)
outdir.mkdir(parents=True, exist_ok=True)

def extract_nsec(fname):
    m = re.search(r'[nN][sS]ec0*(\d+)', fname)
    return int(m.group(1)) if m else 0

def load_records(directory, pattern):
    indir = Path(directory)
    files = sorted(indir.glob(pattern), key=lambda p: extract_nsec(p.name))
    if not files:
        print(f"  WARNING: No files in {indir}")
        return []
    records = []
    for fpath in files:
        nsec = extract_nsec(fpath.name)
        try:
            with uproot.open(str(fpath)) as rf:
                nhits = ak.to_numpy(
                    rf["hibeam"]["ProtoTPC/nHits"].array(library="ak"))
        except Exception as e:
            print(f"  [SKIP] nSec={nsec}: {e}")
            continue
        n_total  = len(nhits)
        nz       = nhits[nhits > 0]
        n_zero   = int((nhits == 0).sum())
        # Mode = most frequent value in non-zero hits
        if len(nz) > 0:
            vals, cnts = np.unique(nz, return_counts=True)
            mode_val   = int(vals[np.argmax(cnts)])
        else:
            mode_val = 0
        records.append({
            "nsec"     : nsec,
            "mean_all" : float(nhits.mean()),
            "mean_nz"  : float(nz.mean()) if len(nz) > 0 else 0.0,
            "mode_nz"  : float(mode_val),
            "zero_frac": n_zero / n_total,
        })
        print(f"  nSec={nsec:>3}  mean(all)={records[-1]['mean_all']:.2f}"
              f"  mean(nz)={records[-1]['mean_nz']:.2f}"
              f"  mode(nz)={mode_val}"
              f"  zero={records[-1]['zero_frac']*100:.1f}%")
    return records

print("\nLoading muon data...")
muon_rec = load_records(args.muon_dir,   args.pattern)
print("\nLoading Krakow data...")
krak_rec = load_records(args.krakow_dir, args.pattern)

if not muon_rec and not krak_rec:
    sys.exit("No data loaded.")

plt.rcParams.update({
    "font.family"       : "serif",
    "font.size"         : 12,
    "axes.labelsize"    : 13,
    "xtick.labelsize"   : 11,
    "ytick.labelsize"   : 11,
    "axes.linewidth"    : 0.8,
    "xtick.major.width" : 0.8,
    "ytick.major.width" : 0.8,
    "xtick.major.size"  : 4.0,
    "ytick.major.size"  : 4.0,
    "xtick.minor.size"  : 2.0,
    "ytick.minor.size"  : 2.0,
    "xtick.direction"   : "in",
    "ytick.direction"   : "in",
    "xtick.top"         : True,
    "ytick.right"       : True,
})

C_ALL  = "#888888"
C_NZ   = "#2166ac"
C_MODE = "#d73027"
LS_ALL  = "-"
LS_NZ   = "-"
LS_MODE = "--"
M_ALL  = "s"
M_NZ   = "o"
M_MODE = "^"
MS     = 6
LW     = 1.6

fig, (ax_mu, ax_kr) = plt.subplots(1, 2, figsize=(11.0, 4.5),
                                    sharey=False)
fig.subplots_adjust(left=0.08, right=0.97, bottom=0.22,
                    top=0.93, wspace=0.32)

def plot_panel(ax, records, title, label_suffix=""):
    if not records:
        ax.text(0.5, 0.5, "No data", transform=ax.transAxes,
                ha="center", va="center", fontsize=12, color="grey")
        ax.set_title(title); return [], []

    nsecs    = [r["nsec"]    for r in records]
    mean_all = [r["mean_all"] for r in records]
    mean_nz  = [r["mean_nz"]  for r in records]
    mode_nz  = [r["mode_nz"]  for r in records]

    h1, = ax.plot(nsecs, mean_all, marker=M_ALL, color=C_ALL,
                  markersize=MS, linewidth=LW, linestyle=LS_ALL, zorder=3)
    h2, = ax.plot(nsecs, mean_nz,  marker=M_NZ,  color=C_NZ,
                  markersize=MS, linewidth=LW, linestyle=LS_NZ,  zorder=4)
    h3, = ax.plot(nsecs, mode_nz,  marker=M_MODE, color=C_MODE,
                  markersize=MS, linewidth=LW, linestyle=LS_MODE, zorder=5)

    ax.set_xlabel("Number of segments $N_{\\rm seg}$")
    ax.set_ylabel(r"$n_{\rm hits}$ per event")
    ax.set_title(title, pad=6)
    ax.xaxis.set_major_locator(ticker.MaxNLocator(integer=True, nbins=6))
    ax.yaxis.set_major_locator(ticker.MaxNLocator(nbins=6))
    ax.set_xlim(nsecs[0] - 1, nsecs[-1] + 1)
    ax.set_ylim(bottom=0)

    # Panel label
    lbl = "(a)" if "uon" in title else "(b)"
    ax.text(0.03, 0.97, lbl, transform=ax.transAxes,
            ha="left", va="top", fontsize=13, fontweight="bold")

    return [h1, h2, h3], [
        r"Mean, all events",
        r"Mean, $n_{\rm hits}>0$ only",
        r"Mode, $n_{\rm hits}>0$ only",
    ]

handles_mu, labels_mu = plot_panel(ax_mu, muon_rec,
    "Muon lab  (MCPL cosmic muons)")
handles_kr, labels_kr = plot_panel(ax_kr, krak_rec,
    r"Krakow scattering  ($^{2}$H(p,p), 190 MeV)")

# Single shared legend centred below both panels
handles = handles_mu if handles_mu else handles_kr
labels  = labels_mu  if labels_mu  else labels_kr
if handles:
    fig.legend(handles, labels,
               loc="lower center",
               ncol=3,
               bbox_to_anchor=(0.5, 0.01),
               framealpha=0.92,
               edgecolor="#aaaaaa",
               fontsize=11,
               handlelength=2.0,
               columnspacing=1.5)

for ext in [".pdf", ".png"]:
    path = outdir / f"nhits_vs_nsections_comparison{ext}"
    fig.savefig(path, dpi=300 if ext == ".png" else None,
                bbox_inches="tight")
    print(f"\nSaved: {path}")

plt.close(fig)