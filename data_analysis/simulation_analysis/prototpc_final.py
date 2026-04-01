#!/usr/bin/env python3
"""
prototypc_dedx_analysis.py
==========================
Wrapper around tpc_dedx_analysis.py (unchanged) for experimental
arrangements that use the ProtoTPC rather than the main HIBEAM TPC:
  - Krakow scattering setup  (proton-deuteron elastic scattering)
  - Muon lab setup           (cosmic / MCPL muon beam)

The only difference from tpc_dedx_analysis.py is that load_tpc_data()
explicitly excludes ProtoTPC branches (via "not k.startswith('ProtoTPC')"),
so we bypass it and build the identical data dict directly from the
ProtoTPC sub-branches.  All downstream functions — compute_dedx,
fit_landau, plot_dedx — are called completely unchanged.

Scientific notes
----------------
- The 70% truncated-mean method is applied identically to the main TPC
  analysis [ALICE JINST 5 (2010) P09002].
- A low-energy pre-cut (bottom 10th percentile of non-zero events) is
  applied before fitting to suppress the secondary peak from partial tracks
  and delta-ray electrons that do not represent the primary dE/dx signal.
  This is standard in TPC dE/dx analyses where the active volume is short
  relative to typical track lengths [Bichsel, Rev.Mod.Phys. 60 (1988) 663].
- The ProtoTPC geometry constants (pad dimensions, drift velocity) can be
  supplied via --pad_width, --pad_height, --drift_v, --time_bin to enable
  true dE/dx [MeV/mm] output rather than sum(Edep) [MeV].

Usage
-----
    # Krakow scattering data (sum Edep, no geometry)
    python3 prototypc_dedx_analysis.py KrakowScatter.root

    # Muon lab with geometry → true dE/dx [MeV/mm]
    python3 prototypc_dedx_analysis.py muons.root --output muon_dedx.png \\
        --pad_width 1.0 --pad_height 1.0 --drift_v 0.077 --time_bin 100

    # Save figures
    python3 prototypc_dedx_analysis.py KrakowScatter.root --output krakow.png
"""

import argparse
import sys
import warnings
from pathlib import Path

import awkward as ak
import numpy as np
import uproot

# ── Import tpc_dedx_analysis unchanged ───────────────────────────────────────
try:
    import final as tda
except ImportError:
    sys.exit(
        "ERROR: tpc_dedx_analysis.py not found.\n"
        "Place it in the same directory as this script."
    )

# ── Argument parsing ──────────────────────────────────────────────────────────
parser = argparse.ArgumentParser(
    description=__doc__,
    formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument("filepath",
    help="Path to ROOT file containing ProtoTPC branch (hibeam tree)")
parser.add_argument("--output", default=None,
    help="Output figure path stem (two files written: _linear.png / _log.png)")
parser.add_argument("--tree",   default="hibeam",
    help="TTree name inside the ROOT file (default: hibeam)")
parser.add_argument("--trunc",  default=0.70, type=float,
    help="Truncated-mean fraction (default: 0.70)")
parser.add_argument("--nbins",  default=100, type=int,
    help="Histogram bins for fit (default: 100)")
parser.add_argument("--minsteps", default=5, type=int,
    help="Minimum G4 steps per event to accept (default: 5)")
parser.add_argument("--low_cut_pct", default=10.0, type=float,
    help="Low-energy pre-cut: remove bottom N%% of non-zero events "
         "to suppress partial-track / delta-ray peak (default: 10)")
# Optional geometry — enables true dE/dx [MeV/mm]
parser.add_argument("--pad_width",  default=None, type=float,
    help="Pad width [mm]  (enables dE/dx mode)")
parser.add_argument("--pad_height", default=None, type=float,
    help="Pad height [mm] (enables dE/dx mode)")
parser.add_argument("--drift_v",    default=None, type=float,
    help="Drift velocity [mm/ns]")
parser.add_argument("--time_bin",   default=None, type=float,
    help="Time bin width [ns]")
args = parser.parse_args()

# ══════════════════════════════════════════════════════════════════════════════
# STEP 1 — Load ProtoTPC branches directly
#          (bypasses load_tpc_data which excludes ProtoTPC)
# ══════════════════════════════════════════════════════════════════════════════
print(f"\nLoading ProtoTPC data from: {args.filepath}")

try:
    with uproot.open(args.filepath) as root_file:
        available = list(root_file.keys())
        if args.tree not in root_file:
            sys.exit(
                f"ERROR: TTree '{args.tree}' not found.\n"
                f"Available keys: {available}"
            )
        tree = root_file[args.tree]

        # Verify ProtoTPC branch exists
        all_keys = tree.keys(recursive=True)
        proto_keys = [k for k in all_keys if k.startswith("ProtoTPC")]
        if not proto_keys:
            sys.exit(
                f"ERROR: No ProtoTPC branches found in tree '{args.tree}'.\n"
                f"Available keys: {all_keys}"
            )

        print(f"  Found ProtoTPC sub-branches: {proto_keys}")

        # Load mandatory branch
        edep_jag = tree["ProtoTPC/Edep"].array(library="ak")

        # Load optional branches — warn if absent
        def safe_load(key):
            if key in all_keys:
                return tree[key].array(library="ak")
            warnings.warn(f"Branch '{key}' not found — track-length "
                          "normalisation unavailable.", stacklevel=2)
            return None

        pad_row = safe_load("ProtoTPC/padRow")
        pad_col = safe_load("ProtoTPC/padColumn")
        ts      = safe_load("ProtoTPC/timestamp")

        n_events = len(edep_jag)

except (OSError, KeyError) as e:
    sys.exit(f"ERROR opening file: {e}")

print(f"  Tree contains {n_events:,} events")

# Build the data dict expected by tda.compute_dedx
data = {
    "edep"      : edep_jag,
    "pad_row"   : pad_row,
    "pad_col"   : pad_col,
    "timestamp" : ts,
    "n_el"      : None,
    "n_events"  : n_events,
}

# ══════════════════════════════════════════════════════════════════════════════
# STEP 2 — Build geometry (optional)
# ══════════════════════════════════════════════════════════════════════════════
geom = None
if all(v is not None for v in
       [args.pad_width, args.pad_height, args.drift_v, args.time_bin]):
    geom = tda.TPCGeometry(
        pad_width  = args.pad_width,
        pad_height = args.pad_height,
        drift_v    = args.drift_v,
        time_bin   = args.time_bin,
    )
    print(f"  Geometry supplied → dE/dx [MeV/mm] mode")
    print(f"    pad_width={args.pad_width} mm  pad_height={args.pad_height} mm"
          f"  drift_v={args.drift_v} mm/ns  time_bin={args.time_bin} ns")
else:
    print("  No geometry supplied → sum(Edep) [MeV] mode")
    print("  (Supply --pad_width --pad_height --drift_v --time_bin for dE/dx)")

use_dedx = geom is not None and geom.is_complete

# ══════════════════════════════════════════════════════════════════════════════
# STEP 3 — Compute per-event energy loss (tda.compute_dedx unchanged)
# ══════════════════════════════════════════════════════════════════════════════
print(f"\nComputing per-event energy loss (min_steps={args.minsteps})...")
values = tda.compute_dedx(data, geom=geom, min_steps=args.minsteps)

# Remove zero-energy events
values    = values[values > 0]
n_hit     = len(values)
zero_frac = 1.0 - n_hit / n_events
print(f"  {n_hit:,} / {n_events:,} events with E_dep > 0  "
      f"({zero_frac*100:.1f}% zero-hit events excluded)")

if n_hit < 200:
    warnings.warn(f"Only {n_hit} hit events — fit may be unreliable.", stacklevel=1)

# ── GeV → MeV conversion ─────────────────────────────────────────────────────
# The hibeam_g4 simulation stores Edep in GeV (Geant4 internal units).
# Multiply by 1000 so all axes, fit parameters and annotations are in MeV.
values = values * 1000.0
print(f"  Energy converted GeV → MeV  "
      f"(range: {values.min():.3f} – {values.max():.3f} MeV)")

# ── Low-energy pre-cut ────────────────────────────────────────────────────────
# Removes partial-track / delta-ray secondary peak.
# Reference: Bichsel, Rev.Mod.Phys. 60 (1988) 663, Section IV.
if args.low_cut_pct > 0:
    low_cut  = np.percentile(values, args.low_cut_pct)
    n_before = len(values)
    values   = values[values > low_cut]
    print(f"  Low-energy pre-cut ({args.low_cut_pct:.0f}th pct, "
          f"{low_cut:.3f} MeV): {n_before - len(values):,} removed → "
          f"{len(values):,} remain")

# ══════════════════════════════════════════════════════════════════════════════
# STEP 4 — Fit Landau (tda.fit_landau unchanged)
# ══════════════════════════════════════════════════════════════════════════════
print(f"\nFitting Landau (Moyal) to bottom {args.trunc*100:.0f}% of values...")

try:
    fit_result = tda.fit_landau(
        values,
        truncation = args.trunc,
        n_bins     = args.nbins,
    )
except (ValueError, RuntimeError) as e:
    print(f"\nFit failed: {e}")
    print("  Try --nbins 50  or  --low_cut_pct 15  or  --minsteps 1")
    sys.exit(1)

# ══════════════════════════════════════════════════════════════════════════════
# STEP 5 — Print results
# ══════════════════════════════════════════════════════════════════════════════
chi2_tail_red = fit_result["chi2_tail"] / max(fit_result["ndf_tail"], 1)
units = "MeV/mm" if use_dedx else "MeV"
trunc_pct = int(args.trunc * 100)

print("\n" + "─" * 58)
print("  FIT RESULTS  (ProtoTPC)")
print("─" * 58)
print(f"  File              : {Path(args.filepath).name}")
print(f"  Mode              : {'dE/dx [MeV/mm]' if use_dedx else 'sum(Edep) [MeV]'}")
print(f"  Events (hit)      : {n_hit:,}  ({zero_frac*100:.1f}% zero-hit)")
print(f"  Events (fit)      : {len(values):,}  "
      f"(after {args.low_cut_pct:.0f}% low-cut + {trunc_pct}% truncation)")
print(f"  MPV               : {fit_result['mpv']:.4g} {units}")
print(f"  Width ξ           : {fit_result['scale']:.4g} ± "
      f"{fit_result['scale_err']:.3g} {units}")
print(f"  chi²/ndf (fit)    : {fit_result['chi2_red']:.3f}"
      f"  (p = {fit_result['p_value']:.4f})")
print(f"  chi²/ndf (tail)   : {chi2_tail_red:.3f}"
      f"  (p = {fit_result['p_tail']:.4f})")
print("─" * 58 + "\n")

# ══════════════════════════════════════════════════════════════════════════════
# STEP 6 — Publication-quality plot  (mplhep CMS/ATLAS style)
# ══════════════════════════════════════════════════════════════════════════════
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.ticker as ticker
import mplhep as hep
from scipy import stats as scipy_stats

# Use CMS style as the base — standard for HEP publications
plt.style.use(hep.style.CMS)
# Override a few settings for a cleaner single-panel figure
plt.rcParams.update({
    "font.size"       : 14,
    "axes.labelsize"  : 15,
    "xtick.labelsize" : 12,
    "ytick.labelsize" : 12,
    "legend.fontsize" : 11,
    "xtick.direction" : "in",
    "ytick.direction" : "in",
    "xtick.top"       : True,
    "ytick.right"     : True,
})

def _moyal(x, loc, scale, norm):
    return norm * scipy_stats.moyal.pdf(x, loc=loc, scale=scale)

def make_publication_plot(fit_result, use_dedx, output_path, units):
    bc_orig = fit_result["bin_centers"]
    bw_orig = fit_result["bin_width"]
    popt    = fit_result["popt"]
    mpv     = fit_result["mpv"]
    trunc_pct = int(args.trunc * 100)

    # ── Rebin to coarser bins appropriate for publication ────────────────────
    # The fit used fine bins (100) for precision — for display we merge every
    # 4 bins. This gives ~25 visible bars, which is standard for a Landau
    # distribution in a journal figure (e.g. ALICE JINST 5 P09002 Fig. 3).
    REBIN = 4
    n_new  = len(bc_orig) // REBIN
    # Trim to exact multiple
    cnt_raw = fit_result["counts"][:n_new * REBIN]
    fc_raw  = fit_result["fit_curve"][:n_new * REBIN]
    bc_raw  = bc_orig[:n_new * REBIN]

    cnt_new = cnt_raw.reshape(n_new, REBIN).sum(axis=1)
    fc_new  = fc_raw.reshape(n_new, REBIN).sum(axis=1)
    bc_new  = bc_raw.reshape(n_new, REBIN).mean(axis=1)
    bw_new  = bw_orig * REBIN
    edges_new = np.concatenate([bc_new - bw_new/2, [bc_new[-1] + bw_new/2]])

    # Fit mask on rebinned bins: bin is "fitted" if its centre ≤ cut value
    fm_new   = bc_new <= fit_result["cut_value"]

    cnt_err_new = np.sqrt(np.maximum(cnt_new, 1.0))
    fc_safe     = np.maximum(fc_new, 1e-6)
    pulls_new   = (cnt_new - fc_new) / np.sqrt(fc_safe)
    pull_err_new= cnt_err_new / np.sqrt(fc_safe)

    # ── Display range: rising edge to 6×MPV ──────────────────────────────────
    x_min   = max(float(bc_new[0]) - bw_new, 0.0)
    x_max   = min(float(bc_new[-1]) + bw_new, mpv * 6.0)
    x_dense = np.linspace(max(x_min, mpv - 4*fit_result["scale"]),
                          x_max, 3000)
    y_dense = _moyal(x_dense, *popt)

    C_FIT  = "#2166ac"
    C_TAIL = "#d73027"

    x_label = (
        r"$\mathrm{d}E/\mathrm{d}x\ [\mathrm{MeV}\,\mathrm{mm}^{-1}]$"
        if use_dedx else
        r"$\sum E_{\mathrm{dep}}\ [\mathrm{MeV}]$"
    )

    for log_y, suffix in [(False, "linear"), (True, "log")]:
        fig = plt.figure(figsize=(8, 7))
        gs  = gridspec.GridSpec(
            2, 1, height_ratios=[3, 1], hspace=0.0,
            left=0.12, right=0.97, top=0.92, bottom=0.11)
        ax_top = fig.add_subplot(gs[0])
        ax_res = fig.add_subplot(gs[1], sharex=ax_top)

        # ── Top panel: histograms ─────────────────────────────────────────────
        # Draw fitted and excluded regions separately using bar + hep errorbar
        ax_top.bar(bc_new[fm_new],  cnt_new[fm_new],  width=bw_new * 0.92,
                   color=C_FIT,  alpha=0.55, linewidth=0,
                   label=f"Data — fitted ({trunc_pct}% truncation)")
        ax_top.bar(bc_new[~fm_new], cnt_new[~fm_new], width=bw_new * 0.92,
                   color=C_TAIL, alpha=0.30, hatch="////",
                   edgecolor=C_TAIL, linewidth=0.5,
                   label="Excluded tail")

        # Error bars via hep.histplot errorbar mode for correct HEP style
        hep.histplot(cnt_new[fm_new],
                     bins=np.linspace(bc_new[fm_new][0] - bw_new/2,
                                      bc_new[fm_new][-1] + bw_new/2,
                                      fm_new.sum() + 1),
                     yerr=cnt_err_new[fm_new],
                     ax=ax_top, histtype="errorbar",
                     color="black", linewidth=0, markersize=3)

        # Fit curve
        ax_top.plot(x_dense, y_dense, color="black", lw=2.0,
                    label="Landau (Moyal) fit", zorder=5)
        ax_top.axvline(fit_result["cut_value"], color="dimgray",
                       lw=1.2, ls="--", alpha=0.8)

        # Annotation
        ann = "\n".join([
            f"MPV $= {fit_result['mpv']:.3g}$ {units}",
            f"$\\xi = {fit_result['scale']:.3g}"
            f"\\pm{fit_result['scale_err']:.2g}$ {units}",
            f"$\\chi^2/\\nu = {fit_result['chi2_red']:.2f}$,"
            f"  $p = {fit_result['p_value']:.3f}$",
        ])
        ax_top.text(0.97, 0.96, ann, transform=ax_top.transAxes,
                    va="top", ha="right", fontsize=11,
                    bbox=dict(boxstyle="round,pad=0.4", facecolor="white",
                              edgecolor="#aaaaaa", alpha=0.95))

        # ESS / HIBEAM branding via mplhep
        label_out = hep.cms.label(loc=0, data=False,
                                   label="Simulation", rlabel="HIBEAM",
                                   ax=ax_top)
        try:
            label_out[0].set_text("ESS")
        except (TypeError, IndexError):
            pass

        y_label = "Events / bin" + (" (log scale)" if log_y else "")
        ax_top.set_ylabel(y_label)
        ax_top.legend(loc="upper right", framealpha=0.92,
                      edgecolor="#aaaaaa", handlelength=1.8,
                      borderpad=0.5, labelspacing=0.3)
        ax_top.set_xlim(x_min, x_max)
        if log_y:
            ax_top.set_yscale("log")
            ax_top.set_ylim(bottom=0.5)
        else:
            ax_top.yaxis.set_major_locator(ticker.MaxNLocator(nbins=5))
        plt.setp(ax_top.get_xticklabels(), visible=False)

        # ── Residuals panel ───────────────────────────────────────────────────
        # Only plot the fitted region — tail residuals are by design large
        # and not informative about fit quality
        fm_idx = np.where(fm_new)[0]
        ax_res.errorbar(
            bc_new[fm_new],
            pulls_new[fm_new],
            yerr=pull_err_new[fm_new],
            fmt="o", color=C_FIT,
            markersize=4, linewidth=1.0, capsize=2,
            label=None)
        ax_res.axhline(0,  color="black",   lw=1.0)
        ax_res.axhline(+2, color="#888888", lw=1.0, ls="--",
                       label=r"$\pm 2\sigma$")
        ax_res.axhline(-2, color="#888888", lw=1.0, ls="--")
        ax_res.axvline(fit_result["cut_value"], color="dimgray",
                       lw=1.2, ls="--", alpha=0.8)

        # Set y-axis to include all pulls with ±2σ lines always visible
        pull_max = max(float(np.abs(pulls_new[fm_new]).max()), 2.5)
        ylim     = float(np.ceil(pull_max / 2.0) * 2.0) + 0.5
        ax_res.set_ylim(-ylim, ylim)
        ax_res.yaxis.set_major_locator(
            ticker.MultipleLocator(2.0))
        ax_res.legend(loc="upper right", fontsize=9,
                      framealpha=0.92, edgecolor="#aaaaaa")
        ax_res.set_ylabel(r"$(N-F)/\!\sqrt{F}$")
        ax_res.set_xlabel(x_label)
        ax_res.set_xlim(x_min, x_max)

        if output_path:
            import os
            stem, ext = os.path.splitext(output_path)
            if not ext:
                ext = ".png"
            for save_ext in [ext, ".pdf"]:
                path = f"{stem}_{suffix}{save_ext}"
                dpi  = 300 if save_ext == ".png" else None
                fig.savefig(path, dpi=dpi, bbox_inches="tight")
                print(f"  Saved: {path}")
        else:
            plt.show()
        plt.close(fig)

make_publication_plot(fit_result, use_dedx, args.output, units)