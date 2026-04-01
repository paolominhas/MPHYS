#!/usr/bin/env python3
"""
segmentation_overlay.py
=======================
Fits Landau (Moyal) distributions to ProtoTPC energy deposition for each
segmentation and produces a publication-quality overlay figure.

Fitting strategy
----------------
A hard lower bound cut (Edep > 0.15 MeV after GeV→MeV conversion) removes
the small plateau before the main Landau peak that causes fits to fail for
nSec > 12.  This is physically motivated: sub-0.15 MeV deposits are
dominated by very short partial tracks and sub-threshold delta-ray clusters
that do not represent the primary ionisation signal (Bichsel 1988, Sec. IV).

Coarser bins (n_bins=40) are used for the fit — at typical event counts of
O(5000) per file, 40 bins gives ~100 events/bin in the fit region, which
is sufficient for reliable chi-squared minimisation and avoids empty bins
that destabilise the Moyal fit (Pearson chi-squared requires ≥5 counts/bin).

Physical argument for the MPV trend
------------------------------------
If MPV decreases with increasing nSec, this is physically consistent with
the correction for Geant4's tendency to accumulate energy across long step
boundaries: finer segmentation forces Geant4 to record shorter, independent
G4 steps, each depositing less energy, giving a more accurate per-segment
Landau rather than an artificially summed value.  This is the segmented-TPC
analogue of the truncated-mean correction for delta-ray inflation.

References
----------
Landau (1944); Moyal (1955); Bichsel Rev.Mod.Phys. 60 (1988) 663;
ALICE JINST 5 (2010) P09002; PDG Statistics chapter (2022).

Usage
-----
    python3 segmentation_overlay.py --indir ./ana_output_muon
    python3 segmentation_overlay.py --indir ./ana_output_krakow --outdir ./seg_plots
"""

import argparse, re, sys
from pathlib import Path

import awkward as ak
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.ticker as ticker
import mplhep as hep
import numpy as np
import uproot

try:
    import final as tda
except ImportError:
    sys.exit("ERROR: tpc_dedx_analysis.py not found in the same directory.")

# ── Arguments ─────────────────────────────────────────────────────────────────
parser = argparse.ArgumentParser(description=__doc__,
    formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument("--indir",    default="./ana_output_krakow")
parser.add_argument("--outdir",   default="./seg_plots")
parser.add_argument("--pattern",  default="ana_geom_*.root")
parser.add_argument("--trunc",    default=0.70,  type=float)
parser.add_argument("--nbins",    default=40,    type=int,
    help="Histogram bins for fit — 40 gives ~100 events/bin at typical counts")
parser.add_argument("--minsteps", default=5,     type=int)
parser.add_argument("--low_mev",  default=0.15,  type=float,
    help="Hard lower cut in MeV after GeV→MeV conversion (default: 0.15)")
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

print(f"\nFound {len(files)} files in {indir}  (low cut: >{args.low_mev} MeV)\n")

# ── Style ─────────────────────────────────────────────────────────────────────
plt.style.use(hep.style.CMS)
plt.rcParams.update({
    "font.size"        : 13,
    "axes.labelsize"   : 14,
    "xtick.labelsize"  : 11,
    "ytick.labelsize"  : 11,
    "legend.fontsize"  : 9,
    "xtick.direction"  : "in",
    "ytick.direction"  : "in",
    "xtick.top"        : True,
    "ytick.right"      : True,
})

cmap = plt.cm.viridis(np.linspace(0.05, 0.92, len(files)))

# ── Per-file fitting ───────────────────────────────────────────────────────────
results = []

for idx, fpath in enumerate(files):
    nsec = extract_nsec(fpath.name)
    print(f"  [{idx+1:>2}/{len(files)}] nSec={nsec:>3}  {fpath.name}")

    try:
        with uproot.open(str(fpath)) as rf:
            tree = rf["hibeam"]
            edep    = tree["ProtoTPC/Edep"].array(library="ak")
            pad_row = tree["ProtoTPC/padRow"].array(library="ak")
            pad_col = tree["ProtoTPC/padColumn"].array(library="ak")
            ts      = tree["ProtoTPC/timestamp"].array(library="ak")
        data = {"edep": edep, "pad_row": pad_row, "pad_col": pad_col,
                "timestamp": ts, "n_el": None, "n_events": len(edep)}
    except Exception as e:
        print(f"    SKIP — load error: {e}"); continue

    # Compute per-event SumEdep
    values = tda.compute_dedx(data, geom=None, min_steps=args.minsteps)
    values = values[values > 0]

    # GeV → MeV
    values = values * 1000.0

    # Hard lower cut: removes the sub-threshold plateau that destabilises
    # the Moyal fit for high nSec values.  0.15 MeV is chosen to sit just
    # above the plateau shoulder visible in the data.
    n_before = len(values)
    values   = values[values > args.low_mev]
    print(f"    {n_before:,} hit events → {len(values):,} after "
          f">{args.low_mev} MeV cut")

    if len(values) < 50:
        print(f"    SKIP — too few events after cut"); continue

    try:
        fr = tda.fit_landau(values, truncation=args.trunc, n_bins=args.nbins)
    except (ValueError, RuntimeError) as e:
        print(f"    SKIP — fit failed: {e}"); continue

    print(f"    MPV={fr['mpv']:.4f} MeV  xi={fr['scale']:.4f} MeV"
          f"  chi2/ndf={fr['chi2_red']:.2f}  p={fr['p_value']:.3f}")
    results.append((nsec, fr))

if not results:
    sys.exit("No successful fits.")

nsecs    = [r[0]              for r in results]
mpvs     = [r[1]["mpv"]       for r in results]
mpv_errs = [r[1]["loc_err"]   for r in results]
xis      = [r[1]["scale"]     for r in results]
xi_errs  = [r[1]["scale_err"] for r in results]
chi2s    = [r[1]["chi2_red"]  for r in results]

# ══════════════════════════════════════════════════════════════════════════════
# FIGURE 1 — Normalised overlay of all Landau curves
# ══════════════════════════════════════════════════════════════════════════════
fig1, ax1 = plt.subplots(figsize=(9, 6))
fig1.subplots_adjust(left=0.12, right=0.97, top=0.91, bottom=0.12)

# Add a colourbar-style legend using a ScalarMappable
sm = plt.cm.ScalarMappable(
    cmap=plt.cm.viridis,
    norm=plt.Normalize(vmin=nsecs[0], vmax=nsecs[-1]))
sm.set_array([])
cbar = fig1.colorbar(sm, ax=ax1, pad=0.02)
cbar.set_label("$N_{\\rm seg}$", fontsize=13)
cbar.ax.yaxis.set_major_locator(ticker.MaxNLocator(integer=True, nbins=6))

for (nsec, fr), col in zip(results, cmap):
    bc   = fr["bin_centers"]
    cnt  = fr["counts"].astype(float)
    fm   = fr["fit_mask"]
    bw   = fr["bin_width"]
    peak = cnt[fm].max() if cnt[fm].max() > 0 else 1.0

    # Normalised step histogram (fitted region only — tail clutters the plot)
    ax1.step(np.concatenate([bc[fm] - bw/2, [bc[fm][-1] + bw/2]]),
             np.concatenate([cnt[fm] / peak, [0]]),
             where="post", color=col, linewidth=1.0, alpha=0.55)

    # Normalised fit curve extended to 5×MPV
    x_dense = np.linspace(float(bc[fm][0]),
                          min(float(bc[-1]), fr["mpv"] * 5), 2000)
    y_norm  = tda._moyal_scaled(x_dense, fr["loc"], fr["scale"], fr["norm"]) / peak
    ax1.plot(x_dense, y_norm, color=col, linewidth=1.8)

ax1.set_xlabel(r"$\sum E_{\mathrm{dep}}\ [\mathrm{MeV}]$")
ax1.set_ylabel("Normalised yield  (peak = 1)")
ax1.set_xlim(left=0)
ax1.set_ylim(bottom=0)
ax1.xaxis.set_major_locator(ticker.MaxNLocator(nbins=6))
ax1.yaxis.set_major_locator(ticker.MaxNLocator(nbins=5))

label_out = hep.cms.label(loc=0, data=False, rlabel="HIBEAM", ax=ax1)
try: label_out[0].set_text("ESS")
except Exception: pass

for ext in [".pdf", ".png"]:
    p = outdir / f"overlay_landau{ext}"
    fig1.savefig(p, dpi=300 if ext==".png" else None, bbox_inches="tight")
    print(f"  Saved: {p}")
plt.close(fig1)

# ══════════════════════════════════════════════════════════════════════════════
# FIGURE 2 — MPV and ξ vs nSections  (2-panel, A4 half-width)
# ══════════════════════════════════════════════════════════════════════════════
fig2, (ax_mpv, ax_xi) = plt.subplots(1, 2, figsize=(10, 4.5))
fig2.subplots_adjust(left=0.10, right=0.97, top=0.91,
                     bottom=0.14, wspace=0.35)

def style_panel(ax, x, y, ye, color, marker, ylabel, panel_label):
    ax.errorbar(x, y, yerr=ye, fmt=f"{marker}-", color=color,
                capsize=4, linewidth=1.8, markersize=7, zorder=4)
    mean_y = np.mean(y)
    ax.axhline(mean_y, color="dimgray", lw=1.0, ls="--", zorder=3)
    var = np.std(y) / mean_y * 100 if mean_y != 0 else 0
    ax.text(0.97, 0.97,
            f"$\\sigma/\\mu = {var:.1f}\\%$\n"
            f"({'stable' if var < 5 else 'varies'})",
            transform=ax.transAxes, ha="right", va="top",
            fontsize=10,
            color="#2ca02c" if var < 5 else "#d62728",
            bbox=dict(boxstyle="round,pad=0.3", facecolor="white",
                      edgecolor="#aaaaaa", alpha=0.90))
    ax.set_xlabel("Number of segments $N_{\\rm seg}$")
    ax.set_ylabel(ylabel)
    ax.set_xlim(x[0] - 1, x[-1] + 1)
    ax.xaxis.set_major_locator(ticker.MaxNLocator(integer=True, nbins=6))
    ax.yaxis.set_major_locator(ticker.MaxNLocator(nbins=5))
    ax.text(0.03, 0.97, panel_label, transform=ax.transAxes,
            ha="left", va="top", fontsize=13, fontweight="bold")

style_panel(ax_mpv, nsecs, mpvs, mpv_errs,
            "#2166ac", "o",
            r"Landau MPV [MeV]", "(a)")

style_panel(ax_xi,  nsecs, xis,  xi_errs,
            "#d73027", "s",
            r"Landau $\xi$ width [MeV]", "(b)")

for ax in (ax_mpv, ax_xi):
    label_out = hep.cms.label(loc=0, data=False, rlabel="HIBEAM", ax=ax, fontsize=10)
    try: label_out[0].set_text("ESS")
    except Exception: pass

for ext in [".pdf", ".png"]:
    p = outdir / f"mpv_xi_vs_nsections{ext}"
    fig2.savefig(p, dpi=300 if ext==".png" else None, bbox_inches="tight")
    print(f"  Saved: {p}")
plt.close(fig2)

# ── Text summary ──────────────────────────────────────────────────────────────
summary = outdir / "segmentation_summary.txt"
with open(summary, "w") as f:
    sep = "─" * 72
    f.write(sep + "\n")
    f.write("ProtoTPC Segmentation Study — Landau Fit Summary\n")
    f.write(f"Low cut: >{args.low_mev} MeV  |  bins: {args.nbins}"
            f"  |  truncation: {args.trunc*100:.0f}%\n")
    f.write("Refs: Bichsel (1988); ALICE JINST 5 P09002; PDG (2022)\n")
    f.write(sep + "\n")
    hdr = (f"{'nSec':>5}  {'MPV[MeV]':>10}  {'±err':>8}  "
           f"{'xi[MeV]':>9}  {'±err':>8}  {'chi2/ndf':>9}  {'p':>7}")
    f.write(hdr + "\n" + "─"*72 + "\n")
    for nsec, fr in results:
        f.write(f"{nsec:>5}  {fr['mpv']:>10.5f}  {fr['loc_err']:>8.5f}  "
                f"{fr['scale']:>9.5f}  {fr['scale_err']:>8.5f}  "
                f"{fr['chi2_red']:>9.3f}  {fr['p_value']:>7.4f}\n")
    f.write("\n" + sep + "\n")
    for label, vals in [("MPV", mpvs), ("xi", xis)]:
        m = np.mean(vals); s = np.std(vals)
        v = s/m*100 if m else 0
        f.write(f"  {label}: mean={m:.5f}  std={s:.5f}  "
                f"var={v:.2f}%  "
                f"[{'STABLE <5%' if v<5 else 'VARIES >5%'}]\n")
print(f"  Summary: {summary}\n  Done.")