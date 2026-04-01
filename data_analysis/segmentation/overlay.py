#!/usr/bin/env python3
"""
segmentation_overlay.py
=======================
Wrapper around tpc_dedx_analysis.py (unchanged) that:

  1. Runs the full analysis pipeline on every ana_geom_*.root file in a
     directory, producing the individual linear + log fit plots for each.
  2. Produces an overlay figure comparing the fitted Landau curves and
     MPV / width parameters across all segmentations.

Usage
-----
    python3 segmentation_overlay.py --indir ./ana_output_krakow
    python3 segmentation_overlay.py --indir ./ana_output_krakow --outdir ./seg_plots
"""

import argparse
import re
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np

# ── Import the analysis module unchanged ─────────────────────────────────────
# tpc_dedx_analysis.py must be in the same directory or on PYTHONPATH.
try:
    import final as tda
except ImportError:
    sys.exit(
        "ERROR: tpc_dedx_analysis.py not found.\n"
        "Place it in the same directory as this script."
    )

# ── Argument parsing ──────────────────────────────────────────────────────────
parser = argparse.ArgumentParser(description=__doc__,
    formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument("--indir",   default="./ana_output_krakow",
    help="Directory containing ana_geom_*.root files")
parser.add_argument("--outdir",  default="./seg_plots",
    help="Output directory for all plots (default: ./seg_plots)")
parser.add_argument("--pattern", default="ana_geom_*.root",
    help="Glob pattern for input files (default: ana_geom_*.root)")
parser.add_argument("--branch",  default="ProtoTPC",
    help="TPC branch prefix used by load_tpc_data (default: ProtoTPC)")
parser.add_argument("--trunc",   default=0.70, type=float,
    help="Truncation fraction for Landau fit (default: 0.70)")
parser.add_argument("--nbins",   default=100, type=int,
    help="Histogram bins for fit (default: 100)")
parser.add_argument("--minsteps", default=5, type=int,
    help="Minimum G4 steps per event to accept (default: 5)")
args = parser.parse_args()

indir  = Path(args.indir)
outdir = Path(args.outdir)
outdir.mkdir(parents=True, exist_ok=True)

# ── Collect and sort files ────────────────────────────────────────────────────
def extract_nsec(fname):
    m = re.search(r'[nN][sS]ec0*(\d+)', fname)
    return int(m.group(1)) if m else 0

files = sorted(indir.glob(args.pattern), key=lambda p: extract_nsec(p.name))
if not files:
    sys.exit(f"ERROR: No files matching '{args.pattern}' in {indir}")

print(f"\nFound {len(files)} files in {indir}")
print(f"Output directory: {outdir}\n")

# Colour map — one colour per segmentation
cmap   = plt.cm.viridis(np.linspace(0.1, 0.9, len(files)))

# ── Per-file analysis (uses tpc_dedx_analysis unchanged) ─────────────────────
results = []   # list of (nsec, fit_result) pairs

for idx, fpath in enumerate(files):
    nsec = extract_nsec(fpath.name)
    print(f"  [{idx+1}/{len(files)}] nSections={nsec}  {fpath.name}")

    # ── Step 1: load data directly — bypasses load_tpc_data which
    # explicitly excludes ProtoTPC branches via "not k.startswith('ProtoTPC')"
    # We build the same dict that compute_dedx expects. ──────────────────────
    try:
        import uproot, awkward as ak
        with uproot.open(str(fpath)) as root_file:
            tree = root_file["hibeam"]
            edep     = tree["ProtoTPC/Edep"].array(library="ak")
            pad_row  = tree["ProtoTPC/padRow"].array(library="ak")
            pad_col  = tree["ProtoTPC/padColumn"].array(library="ak")
            ts       = tree["ProtoTPC/timestamp"].array(library="ak")
        data = {
            "edep"     : edep,
            "pad_row"  : pad_row,
            "pad_col"  : pad_col,
            "timestamp": ts,
            "n_el"     : None,
            "n_events" : len(edep),
        }
    except Exception as e:
        print(f"    SKIP — could not load: {e}")
        continue

    # ── Step 2: compute per-event energy loss ─────────────────────────────────
    values = tda.compute_dedx(data, geom=None, min_steps=args.minsteps)
    if len(values) == 0:
        print(f"    SKIP — no events passed minsteps cut")
        continue

    values = values[values > 0]

    # Remove the low-energy left bump (partial tracks / delta rays) that
    # pulls the Moyal fit away from the main Landau peak.
    # Cut below the 10th percentile of non-zero values.
    if len(values) > 50:
        low_cut = np.percentile(values, 10)
        values  = values[values > low_cut]

    print(f"    {len(values):,} events accepted")

    # ── Step 3: fit ───────────────────────────────────────────────────────────
    try:
        fit_result = tda.fit_landau(values,
                                    truncation=args.trunc,
                                    n_bins=args.nbins)
    except (ValueError, RuntimeError) as e:
        print(f"    SKIP — fit failed: {e}")
        continue

    print(f"    MPV={fit_result['mpv']:.5f}  "
          f"xi={fit_result['scale']:.5f}  "
          f"chi2/ndf={fit_result['chi2_red']:.2f}")

    # ── Step 4: individual plots using the module's plot_dedx exactly ─────────
    stem = outdir / f"edep_nSec{nsec:03d}"
    tda.plot_dedx(fit_result,
                  use_dedx=False,
                  output_path=str(stem))
    # plot_dedx saves <stem>_linear.png and <stem>_log.png

    results.append((nsec, fit_result))

if not results:
    sys.exit("No successful fits — nothing to overlay.")

# ── Overlay figure ────────────────────────────────────────────────────────────
# Shows all fitted Landau curves normalised to peak = 1 on one axis,
# plus two summary panels: MPV vs nSections and xi vs nSections.

fig_ov = plt.figure(figsize=(14, 10))
gs = gridspec.GridSpec(2, 2, figure=fig_ov,
                       hspace=0.35, wspace=0.32,
                       left=0.08, right=0.97,
                       top=0.93, bottom=0.08)

ax_curves = fig_ov.add_subplot(gs[0, :])   # full top row — Landau overlay
ax_mpv    = fig_ov.add_subplot(gs[1, 0])   # bottom left  — MPV vs nSec
ax_xi     = fig_ov.add_subplot(gs[1, 1])   # bottom right — xi  vs nSec

# ── Top panel: normalised Landau curves ──────────────────────────────────────
for (nsec, fr), col in zip(results, cmap):
    bc  = fr["bin_centers"]
    cnt = fr["counts"].astype(float)
    fm  = fr["fit_mask"]
    bw  = fr["bin_width"]

    # Normalise histogram to peak = 1
    peak = cnt[fm].max() if cnt[fm].max() > 0 else 1.0
    norm_cnt = cnt / peak

    # Draw fitted region as a step
    ax_curves.step(bc[fm], norm_cnt[fm], where="mid",
                   color=col, linewidth=1.5, alpha=0.6)

    # Draw the Landau fit curve normalised by the same factor
    x_dense = np.linspace(float(bc[0]), float(bc[-1]), 2000)
    y_dense = tda._moyal_scaled(x_dense, fr["loc"], fr["scale"], fr["norm"])
    y_dense_norm = y_dense / peak
    ax_curves.plot(x_dense, y_dense_norm,
                   color=col, linewidth=2.0,
                   label=f"nSec={nsec:>3}  MPV={fr['mpv']:.4f}")

    # Mark the cut threshold lightly
    ax_curves.axvline(fr["cut_value"], color=col,
                      lw=0.6, ls=":", alpha=0.5)

ax_curves.set_xlabel(r"$\sum E_\mathrm{dep}\ [\mathrm{MeV}]$", fontsize=12)
ax_curves.set_ylabel("Normalised yield  (peak = 1)", fontsize=12)
ax_curves.set_title(
    "Landau (Moyal) fits — all segmentations  "
    "(solid = fit, step = data, dotted = 70% cut)",
    fontsize=11)
ax_curves.legend(fontsize=8, ncol=4,
                 loc="upper right",
                 framealpha=0.92, edgecolor="#aaaaaa")
ax_curves.set_xlim(left=0)
ax_curves.grid(True, alpha=0.25)

# ── Bottom left: MPV vs nSections ─────────────────────────────────────────────
nsecs    = [r[0]           for r in results]
mpvs     = [r[1]["mpv"]    for r in results]
mpv_errs = [r[1]["loc_err"] for r in results]
xis      = [r[1]["scale"]   for r in results]
xi_errs  = [r[1]["scale_err"] for r in results]

ax_mpv.errorbar(nsecs, mpvs, yerr=mpv_errs,
                fmt="o-", color="#2166ac",
                capsize=4, linewidth=2, markersize=7,
                label="MPV ± σ")
mean_mpv = np.mean(mpvs)
ax_mpv.axhline(mean_mpv, color="dimgray", lw=1.2, ls="--",
               label=f"Mean = {mean_mpv:.4f} MeV")
var_mpv = np.std(mpvs) / mean_mpv * 100 if mean_mpv != 0 else 0
ax_mpv.text(0.97, 0.05,
            f"σ/μ = {var_mpv:.2f}%\n"
            f"({'STABLE' if var_mpv < 5 else 'VARIES'})",
            transform=ax_mpv.transAxes, ha="right", fontsize=10,
            color="green" if var_mpv < 5 else "red",
            bbox=dict(boxstyle="round", facecolor="white", alpha=0.8))
ax_mpv.set_xlabel("nSections", fontsize=12)
ax_mpv.set_ylabel("Landau MPV [MeV]", fontsize=12)
ax_mpv.set_title("MPV vs Segmentation", fontsize=11)
ax_mpv.legend(fontsize=9)
ax_mpv.grid(True, alpha=0.25)

# ── Bottom right: xi vs nSections ─────────────────────────────────────────────
ax_xi.errorbar(nsecs, xis, yerr=xi_errs,
               fmt="s-", color="#d73027",
               capsize=4, linewidth=2, markersize=7,
               label="ξ (width) ± σ")
mean_xi = np.mean(xis)
ax_xi.axhline(mean_xi, color="dimgray", lw=1.2, ls="--",
              label=f"Mean = {mean_xi:.4f} MeV")
var_xi = np.std(xis) / mean_xi * 100 if mean_xi != 0 else 0
ax_xi.text(0.97, 0.05,
           f"σ/μ = {var_xi:.2f}%\n"
           f"({'STABLE' if var_xi < 5 else 'VARIES'})",
           transform=ax_xi.transAxes, ha="right", fontsize=10,
           color="green" if var_xi < 5 else "red",
           bbox=dict(boxstyle="round", facecolor="white", alpha=0.8))
ax_xi.set_xlabel("nSections", fontsize=12)
ax_xi.set_ylabel("Landau ξ width [MeV]", fontsize=12)
ax_xi.set_title("Landau Width vs Segmentation", fontsize=11)
ax_xi.legend(fontsize=9)
ax_xi.grid(True, alpha=0.25)

fig_ov.suptitle(
    "ProtoTPC Segmentation Study — Landau Fit Parameter Comparison",
    fontsize=13, fontweight="bold")

ov_path = outdir / "overlay_segmentation.png"
fig_ov.savefig(ov_path, dpi=150, bbox_inches="tight")
plt.close(fig_ov)
print(f"\n  Overlay saved: {ov_path}")

# ── Text summary ──────────────────────────────────────────────────────────────
summary_path = outdir / "segmentation_summary.txt"
with open(summary_path, "w") as f:
    sep = "=" * 80
    f.write(sep + "\n")
    f.write("ProtoTPC Segmentation Study — Landau Fit Summary\n")
    f.write("Method: truncated mean (70%), Moyal approximation\n")
    f.write("Refs: Landau (1944); Moyal (1955); ALICE JINST 5 (2010) P09002\n")
    f.write(sep + "\n\n")
    hdr = (f"{'nSec':>5}  {'MPV':>10}  {'MPV_err':>9}  "
           f"{'xi':>10}  {'xi_err':>9}  {'chi2/ndf':>9}  {'p_val':>8}")
    f.write(hdr + "\n" + "-" * 80 + "\n")
    for nsec, fr in results:
        f.write(
            f"{nsec:>5}  {fr['mpv']:>10.6f}  {fr['loc_err']:>9.6f}  "
            f"{fr['scale']:>10.6f}  {fr['scale_err']:>9.6f}  "
            f"{fr['chi2_red']:>9.3f}  {fr['p_value']:>8.4f}\n")
    f.write("\n" + sep + "\n")
    f.write("Stability:\n")
    for label, vals, mean in [("MPV", mpvs, mean_mpv), ("xi", xis, mean_xi)]:
        s = np.std(vals)
        v = s / mean * 100 if mean != 0 else 0
        f.write(f"  {label:<5}: mean={mean:.6f}  std={s:.6f}  "
                f"var={v:.2f}%  "
                f"[{'STABLE <5%' if v < 5 else 'VARIES >5%'}]\n")

print(f"  Summary: {summary_path}")
print(f"\nDone. Individual plots in {outdir}/edep_nSec*.png\n")