#!/usr/bin/env python3
"""
nhits_and_overlay.py
====================
Two analyses built on segmentation_overlay.py:

  1. nHits per event vs nSections:
       - Per-file histogram of nHits showing the zero-hit spike
       - Mean nHits (excluding zero-hit events) vs nSections
       - Fraction of zero-hit events vs nSections

  2. Common-normalised Landau overlay:
       All Edep distributions divided by the SAME normalisation factor
       (total event count of the reference run, nSec=2) so the total
       area of each curve reflects the true hit rate — curves with more
       zero-hit events have less area.

Usage
-----
    python3 nhits_and_overlay.py --indir ./ana_output_krakow
"""

import argparse
import re
import sys
from pathlib import Path

import awkward as ak
import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
import numpy as np
import uproot

try:
    import final as tda
except ImportError:
    sys.exit("ERROR: tpc_dedx_analysis.py not found in the same directory.")

# ── Argument parsing ──────────────────────────────────────────────────────────
parser = argparse.ArgumentParser(description=__doc__,
    formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument("--indir",    default="./ana_output_krakow")
parser.add_argument("--outdir",   default="./seg_plots")
parser.add_argument("--pattern",  default="ana_geom_*.root")
parser.add_argument("--trunc",    default=0.70, type=float)
parser.add_argument("--nbins",    default=100,  type=int)
parser.add_argument("--minsteps", default=1,    type=int,
    help="Min G4 steps to accept (set to 1 so zero-hit events are visible)")
args = parser.parse_args()

indir  = Path(args.indir)
outdir = Path(args.outdir)
outdir.mkdir(parents=True, exist_ok=True)

def extract_nsec(fname):
    m = re.search(r'[nN][sS]ec0*(\d+)', fname)
    return int(m.group(1)) if m else 0

files = sorted(indir.glob(args.pattern), key=lambda p: extract_nsec(p.name))
if not files:
    sys.exit(f"ERROR: No files matching '{args.pattern}' in {indir}")

print(f"\nFound {len(files)} files in {indir}\n")
cmap = plt.cm.viridis(np.linspace(0.1, 0.9, len(files)))

# ── Load all data ─────────────────────────────────────────────────────────────
records = []   # list of dicts, one per file

for idx, fpath in enumerate(files):
    nsec = extract_nsec(fpath.name)
    try:
        with uproot.open(str(fpath)) as rf:
            tree    = rf["hibeam"]
            edep_jag = tree["ProtoTPC/Edep"].array(library="ak")
            nhits_arr= tree["ProtoTPC/nHits"].array(library="ak")
            layer_arr= tree["ProtoTPC/Layer"].array(library="ak")
            pad_row  = tree["ProtoTPC/padRow"].array(library="ak")
            pad_col  = tree["ProtoTPC/padColumn"].array(library="ak")
            ts       = tree["ProtoTPC/timestamp"].array(library="ak")
    except Exception as e:
        print(f"  [SKIP] nSec={nsec}: {e}")
        continue

    nhits      = ak.to_numpy(nhits_arr)
    # Unique layers hit per event — this is the true segmentation metric
    unique_layers = np.array([int(len(set(ak.to_list(layer_arr[i]))))
                              for i in range(len(layer_arr))])          # per-event hit count
    n_total    = len(nhits)
    zero_mask  = nhits == 0
    n_zero     = int(zero_mask.sum())
    zero_frac  = n_zero / n_total if n_total > 0 else 0.0

    nhits_nonzero = nhits[~zero_mask]
    mean_hits_all    = float(nhits.mean())
    mean_hits_nonzero= float(nhits_nonzero.mean()) if len(nhits_nonzero) > 0 else 0.0

    # Edep: sum hits per event (all events, including zero-hit ones)
    edep_per_event = ak.to_numpy(ak.sum(edep_jag, axis=1))

    # Build data dict for tpc_dedx_analysis (same structure as compute_dedx expects)
    data = {
        "edep"     : edep_jag,
        "pad_row"  : pad_row,
        "pad_col"  : pad_col,
        "timestamp": ts,
        "n_el"     : None,
        "n_events" : n_total,
    }

    # Fit using only hit events
    values = tda.compute_dedx(data, geom=None, min_steps=args.minsteps)
    values = values[values > 0]   # exclude zero-energy events

    # Remove the low-energy left bump — this secondary peak is caused by
    # partial tracks and delta rays and pulls the Moyal fit off the main
    # Landau peak. Cut below the 10th percentile of non-zero values.
    if len(values) > 50:
        low_cut = np.percentile(values, 10)
        values  = values[values > low_cut]

    fit_result = None
    if len(values) >= 50:
        try:
            fit_result = tda.fit_landau(values,
                                        truncation=args.trunc,
                                        n_bins=args.nbins)
        except Exception as e:
            print(f"    fit failed: {e}")

    records.append({
        "nsec"            : nsec,
        "idx"             : idx,
        "n_total"         : n_total,
        "nhits"           : nhits,
        "nhits_nonzero"   : nhits_nonzero,
        "unique_layers"   : unique_layers,
        "mean_hits_all"   : mean_hits_all,
        "mean_hits_nonzero": mean_hits_nonzero,
        "zero_frac"       : zero_frac,
        "n_zero"          : n_zero,
        "edep_per_event"  : edep_per_event,
        "values"          : values,
        "fit_result"      : fit_result,
        "fpath"           : fpath,
    })

    print(f"  nSec={nsec:>3}  total={n_total}  zero={n_zero} "
          f"({zero_frac*100:.1f}%)  "
          f"mean_hits(all)={mean_hits_all:.3f}  "
          f"mean_hits(nonzero)={mean_hits_nonzero:.3f}")

if not records:
    sys.exit("No files loaded successfully.")

# ══════════════════════════════════════════════════════════════════════════════
# FIGURE 1 — Per-file nHits histograms (showing zero spike)
#            One panel per segmentation, arranged in a grid
# ══════════════════════════════════════════════════════════════════════════════
n_files = len(records)
n_cols  = min(4, n_files)
n_rows  = (n_files + n_cols - 1) // n_cols

fig1, axes1 = plt.subplots(n_rows, n_cols,
                            figsize=(5 * n_cols, 4 * n_rows),
                            squeeze=False)
fig1.suptitle("nHits per event — all segmentations\n"
              "(zero-hit spike shown; dashed = mean excluding zeros)",
              fontsize=13, fontweight="bold")

for i, rec in enumerate(records):
    ax    = axes1[i // n_cols][i % n_cols]
    col   = cmap[rec["idx"]]
    nhits = rec["nhits"]

    # Clip x-axis at 99th percentile of non-zero hits so the bulk is visible
    nz    = rec["nhits_nonzero"]
    x_max = int(np.percentile(nz, 99)) + 2 if len(nz) > 0 else 50

    # Simple integer bins over visible range only
    bins = np.arange(0, x_max + 1) - 0.5
    ax.hist(nhits[nhits <= x_max], bins=bins,
            color=col, alpha=0.7, edgecolor="none")

    # Zero bar in red
    ax.bar(0, rec["n_zero"], width=0.96,
           color="crimson", alpha=0.8,
           label=f"Zero hits ({rec['zero_frac']*100:.1f}%)")

    # Two mean lines
    ax.axvline(rec["mean_hits_all"], color="black", lw=1.5, ls="--",
               label=f"Mean all = {rec['mean_hits_all']:.1f}")
    ax.axvline(rec["mean_hits_nonzero"], color="lime", lw=1.5, ls="-.",
               label=f"Mean non-zero = {rec['mean_hits_nonzero']:.1f}")

    ax.set_title(f"nSections = {rec['nsec']}", fontsize=10)
    ax.set_xlabel("nHits / event", fontsize=9)
    ax.set_ylabel("Events", fontsize=9)
    ax.set_xlim(-1, x_max)
    ax.legend(fontsize=7, loc="upper right")
    ax.grid(True, alpha=0.25)

# Hide unused subplots
for j in range(n_files, n_rows * n_cols):
    axes1[j // n_cols][j % n_cols].set_visible(False)

fig1.tight_layout()
p1 = outdir / "nhits_per_segmentation.png"
fig1.savefig(p1, dpi=150, bbox_inches="tight")
plt.close(fig1)
print(f"\n  Saved: {p1}")

# ══════════════════════════════════════════════════════════════════════════════
# FIGURE 2 — Mean nHits (excluding zeros) and zero-hit fraction vs nSections
# ══════════════════════════════════════════════════════════════════════════════
fig2, (ax2a, ax2b) = plt.subplots(1, 2, figsize=(13, 5))
fig2.suptitle("nHits vs Segmentation  (zero-hit events excluded from mean)",
              fontsize=13, fontweight="bold")

nsecs       = [r["nsec"]             for r in records]
mean_nz     = [r["mean_hits_nonzero"] for r in records]
mean_all    = [r["mean_hits_all"]    for r in records]
zero_fracs  = [r["zero_frac"] * 100  for r in records]

ax2a.plot(nsecs, mean_nz,  "o-", color="#2166ac", lw=2, ms=8,
          label="Mean nHits (zero-hit events excluded)")
ax2a.plot(nsecs, mean_all, "s--", color="#888888", lw=1.5, ms=6,
          label="Mean nHits (all events, including zeros)")

# Linear expectation from lowest point
if len(nsecs) > 1:
    slope = mean_nz[0] / nsecs[0]
    x_lin = np.linspace(nsecs[0], nsecs[-1], 200)
    ax2a.plot(x_lin, slope * x_lin, "k:", lw=1.5, label="Linear expectation")

ax2a.set_xlabel("nSections", fontsize=12)
ax2a.set_ylabel("Mean hits per event", fontsize=12)
ax2a.set_title("Mean nHits vs nSections", fontsize=11)
ax2a.legend(fontsize=9)
ax2a.grid(True, alpha=0.25)

ax2b.plot(nsecs, zero_fracs, "o-", color="#d73027", lw=2, ms=8)
ax2b.set_xlabel("nSections", fontsize=12)
ax2b.set_ylabel("Fraction of events with zero hits [%]", fontsize=12)
ax2b.set_title("Zero-hit event fraction vs nSections", fontsize=11)
ax2b.grid(True, alpha=0.25)

fig2.tight_layout()
p2 = outdir / "nhits_vs_nsections.png"
fig2.savefig(p2, dpi=150, bbox_inches="tight")
plt.close(fig2)
print(f"  Saved: {p2}")

# ══════════════════════════════════════════════════════════════════════════════
# FIGURE 2b — Mean UNIQUE LAYERS hit per event vs nSections
# This is the true segmentation metric — nHits counts G4 steps (~22 per
# segment), unique layers counts how many distinct readout segments a
# particle crossed.
# ══════════════════════════════════════════════════════════════════════════════
mean_ul_all = []
mean_ul_nz  = []
for rec in records:
    ul = rec["unique_layers"]
    mean_ul_all.append(float(ul.mean()))
    nz = ul[ul > 0]
    mean_ul_nz.append(float(nz.mean()) if len(nz) > 0 else 0.0)

fig2b, ax2c = plt.subplots(figsize=(9, 5))
ax2c.plot(nsecs, mean_ul_nz,  "o-", color="#2166ac", lw=2, ms=8,
          label="Mean unique layers (non-zero events only)")
ax2c.plot(nsecs, mean_ul_all, "s--", color="#888888", lw=1.5, ms=6,
          label="Mean unique layers (all events)")
ax2c.set_xlabel("nSections", fontsize=12)
ax2c.set_ylabel("Mean unique layers hit per event", fontsize=12)
ax2c.set_title("Unique Segments Hit vs Segmentation\n"
               "(nHits counts G4 steps ~22/segment; this counts distinct readout segments)",
               fontsize=11)
ax2c.legend(fontsize=10)
ax2c.grid(True, alpha=0.25)
fig2b.tight_layout()
p2b = outdir / "unique_layers_vs_nsections.png"
fig2b.savefig(p2b, dpi=150, bbox_inches="tight")
plt.close(fig2b)
print(f"  Saved: {p2b}")

# ══════════════════════════════════════════════════════════════════════════════
# FIGURE 3 — Common-normalised Landau overlay
#            All curves divided by the SAME factor = total events in the
#            reference run (nSec=2, first record).  This preserves relative
#            area so curves with more zero-hit events appear shorter.
# ══════════════════════════════════════════════════════════════════════════════
# Reference normalisation: total number of simulated events in first file
ref_n_total = records[0]["n_total"]
print(f"\n  Reference normalisation: {ref_n_total} events (nSec={records[0]['nsec']})")

fig3, ax3 = plt.subplots(figsize=(13, 7))
fig3.suptitle(
    "Landau fits — common normalisation\n"
    f"(all curves divided by N_ref = {ref_n_total:,} events from nSec={records[0]['nsec']})\n"
    "Reduced area indicates more zero-hit events",
    fontsize=12, fontweight="bold")

for rec, col in zip(records, cmap):
    fr = rec["fit_result"]
    if fr is None:
        continue

    bc   = fr["bin_centers"]
    cnt  = fr["counts"].astype(float)
    fm   = fr["fit_mask"]
    bw   = fr["bin_width"]

    # Common normalisation — divide by reference total, not by own peak
    norm_factor = ref_n_total * bw   # converts counts to probability density
    cnt_norm    = cnt / norm_factor

    # Fitted curve — normalised the same way
    x_dense  = np.linspace(float(bc[0]), float(bc[-1]), 2000)
    y_dense  = tda._moyal_scaled(x_dense, fr["loc"], fr["scale"], fr["norm"])
    y_norm   = y_dense / norm_factor

    # Draw fitted region as steps
    ax3.step(bc[fm], cnt_norm[fm], where="mid",
             color=col, linewidth=1.2, alpha=0.5)

    # Draw Landau fit curve
    ax3.plot(x_dense, y_norm, color=col, linewidth=2.0,
             label=(f"nSec={rec['nsec']:>3}  "
                    f"MPV={fr['mpv']:.4f}  "
                    f"hits={rec['zero_frac']*100:.1f}% zero"))

    # Cut threshold
    ax3.axvline(fr["cut_value"], color=col, lw=0.5, ls=":", alpha=0.4)

ax3.set_xlabel(r"$\sum E_\mathrm{dep}\ [\mathrm{MeV}]$", fontsize=13)
ax3.set_ylabel(f"Events / (bin × {ref_n_total:,})", fontsize=12)
ax3.legend(fontsize=8.5, ncol=3, loc="upper right",
           framealpha=0.92, edgecolor="#aaaaaa")
ax3.set_xlim(left=0)
ax3.grid(True, alpha=0.25)

fig3.tight_layout()
p3 = outdir / "overlay_common_normalised.png"
fig3.savefig(p3, dpi=150, bbox_inches="tight")
plt.close(fig3)
print(f"  Saved: {p3}")
print(f"\nAll plots saved to {outdir}/\n")