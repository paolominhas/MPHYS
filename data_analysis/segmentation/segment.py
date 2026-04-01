#!/usr/bin/env python3
"""
segmentation_study.py
=====================
Reads all segmentation files, fits a Landau (Moyal) to each, and plots
the fit parameters (MPV, width, amplitude, chi2/ndf) vs nSections.

Usage:
    python3 segmentation_study.py --indir ./ana_output_krakow
    python3 segmentation_study.py --indir ./ana_output_krakow --branch ProtoTPC/Edep
"""

import sys, re, argparse
from pathlib import Path
import numpy as np
import uproot
import awkward as ak
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import mplhep as hep
from scipy.optimize import curve_fit
from scipy.stats import moyal, chi2

parser = argparse.ArgumentParser()
parser.add_argument("--indir",   default="./ana_output_muon")
parser.add_argument("--outdir",  default="./seg_study")
parser.add_argument("--branch",  default="ProtoTPC/Edep",
    help="Branch to read (default: ProtoTPC/Edep, alt: ProtoTPC_EDep)")
parser.add_argument("--trunc",   default=0.70, type=float)
parser.add_argument("--pattern", default="ana_geom_*.root")
args = parser.parse_args()

indir  = Path(args.indir)
outdir = Path(args.outdir)
outdir.mkdir(exist_ok=True)

def landau_approx(x, amp, loc, scale):
    return amp * moyal.pdf(x, loc=loc, scale=scale)

def extract_nsec(fname):
    m = re.search(r'[nN][sS]ec0*(\d+)', fname)
    return int(m.group(1)) if m else None

files = sorted(indir.glob(args.pattern),
               key=lambda p: extract_nsec(p.name) or 0)
if not files:
    sys.exit(f"No files matching '{args.pattern}' in {indir}")

print(f"\nFound {len(files)} files in {indir}\n")

# ── Per-file fit ───────────────────────────────────────────────────────────────
results = []
plt.style.use(hep.style.CMS)

# One overlay canvas for all fitted distributions
fig_ov, ax_ov = plt.subplots(figsize=(11, 7))
colors = plt.cm.viridis(np.linspace(0.1, 0.9, len(files)))

for idx, fpath in enumerate(files):
    nsec = extract_nsec(fpath.name)
    if nsec is None:
        continue

    try:
        tree        = uproot.open(f"{fpath}:hibeam")
        edep_jagged = tree[args.branch].array()
    except Exception as e:
        print(f"  [SKIP] {fpath.name}: {e}")
        continue

    # Sum hits per event to get total Edep per event
    energy       = ak.to_numpy(ak.sum(edep_jagged, axis=1))

    # Also try per-hit flat for nHits proxy
    edep_flat    = ak.to_numpy(ak.flatten(edep_jagged))
    nhits_mean   = len(edep_flat[edep_flat > 0]) / max(len(energy), 1)

    # Weights
    try:
        wjag   = tree["PrimaryWeight"].array()
        weights = ak.to_numpy(ak.fill_none(ak.firsts(wjag), 1.0))
    except Exception:
        weights = np.ones(len(energy))

    hit_mask     = energy > 0
    energy_hits  = energy[hit_mask]
    weights_hits = weights[hit_mask]

    if len(energy_hits) < 20:
        print(f"  [SKIP] nSec={nsec}: only {len(energy_hits)} hit events")
        continue

    cut_value  = np.percentile(energy_hits, args.trunc * 100)
    cut_mask   = energy_hits <= cut_value
    e_cut      = energy_hits[cut_mask]
    w_cut      = weights_hits[cut_mask]

    counts, edges = np.histogram(e_cut, bins=50, weights=w_cut)
    bin_centers   = (edges[:-1] + edges[1:]) / 2
    sumw2, _      = np.histogram(e_cut, bins=edges, weights=w_cut**2)
    errors        = np.sqrt(sumw2)
    errors[errors == 0] = 1

    p0 = [max(counts), bin_centers[np.argmax(counts)],
          (cut_value - e_cut.min()) * 0.05]
    try:
        popt, pcov = curve_fit(
            landau_approx, bin_centers, counts, p0=p0, sigma=errors,
            bounds=([0, e_cut.min(), 1e-9], [np.inf, cut_value, cut_value]),
            maxfev=20000)
        perr = np.sqrt(np.diag(pcov))

        # Chi2 / NDF
        y_fit  = landau_approx(bin_centers, *popt)
        nz     = errors > 1   # non-trivial bins
        chi2_v = np.sum(((counts[nz] - y_fit[nz]) / errors[nz])**2)
        ndf    = np.sum(nz) - 3
        chi2ndf = chi2_v / ndf if ndf > 0 else 999

        print(f"  nSec={nsec:>3}  MPV={popt[1]:.5f}±{perr[1]:.5f}  "
              f"xi={popt[2]:.5f}±{perr[2]:.5f}  "
              f"chi2/ndf={chi2ndf:.2f}  nHits/evt={nhits_mean:.1f}")

        results.append({
            "nsec"     : nsec,
            "mpv"      : popt[1], "mpv_err"  : perr[1],
            "xi"       : popt[2], "xi_err"   : perr[2],
            "amp"      : popt[0], "amp_err"  : perr[0],
            "chi2ndf"  : chi2ndf,
            "nhits"    : nhits_mean,
            "nevents"  : int(len(energy_hits)),
        })

    except RuntimeError:
        print(f"  nSec={nsec:>3}  fit failed")
        continue

    # Add to overlay
    norm_counts = counts / counts.max()
    x_sm  = np.linspace(edges[0], edges[-1], 500)
    y_sm  = landau_approx(x_sm, *popt)
    y_sm_n = y_sm / popt[0]   # normalise to peak=1

    ax_ov.step(bin_centers, norm_counts, where="mid",
               color=colors[idx], linewidth=1.5, alpha=0.7)
    ax_ov.plot(x_sm, y_sm_n, color=colors[idx], linewidth=2,
               linestyle="--", label=f"nSec={nsec}  MPV={popt[1]:.4f}")

ax_ov.set_xlabel("Energy Deposited [GeV]", fontsize=13)
ax_ov.set_ylabel("Normalised yield (peak=1)", fontsize=13)
ax_ov.set_title("Landau fits — all segmentations", fontweight="bold")
ax_ov.legend(fontsize=9, ncol=2)
ax_ov.grid(True, alpha=0.3)
fig_ov.tight_layout()
fig_ov.savefig(outdir / "overlay_all_fits.png", dpi=150, bbox_inches="tight")
plt.close(fig_ov)

if not results:
    sys.exit("No successful fits.")

results.sort(key=lambda r: r["nsec"])
nsecs    = [r["nsec"]     for r in results]
mpvs     = [r["mpv"]      for r in results]
mpv_errs = [r["mpv_err"]  for r in results]
xis      = [r["xi"]       for r in results]
xi_errs  = [r["xi_err"]   for r in results]
chi2ndfs = [r["chi2ndf"]  for r in results]
nhits    = [r["nhits"]    for r in results]

# ── Summary figure: 4 panels ──────────────────────────────────────────────────
plt.style.use(hep.style.CMS)
fig, axes = plt.subplots(2, 2, figsize=(14, 10))
fig.suptitle("ProtoTPC Segmentation Study — Landau Fit Parameters vs nSections",
             fontweight="bold", fontsize=14)

def plot_panel(ax, x, y, ye, ylabel, title, color):
    ax.errorbar(x, y, yerr=ye, fmt="o-", color=color,
                capsize=5, linewidth=2, markersize=8)
    # Reference line at mean
    mean_y = np.mean(y)
    ax.axhline(mean_y, color="grey", linestyle="--", linewidth=1,
               label=f"Mean = {mean_y:.5f}")
    ax.set_xlabel("nSections", fontsize=12)
    ax.set_ylabel(ylabel, fontsize=12)
    ax.set_title(title, fontsize=12)
    ax.legend(fontsize=10)
    ax.grid(True, alpha=0.3)
    # Annotate variation
    std_y  = np.std(y)
    var_pct = std_y / mean_y * 100 if mean_y != 0 else 0
    ax.text(0.98, 0.05,
            f"σ/μ = {var_pct:.2f}%\n({'STABLE' if var_pct < 5 else 'VARIES'})",
            transform=ax.transAxes, ha="right", fontsize=10,
            color="green" if var_pct < 5 else "red",
            bbox=dict(boxstyle="round", facecolor="white", alpha=0.7))

plot_panel(axes[0,0], nsecs, mpvs,  mpv_errs,
           "Landau MPV [GeV]", "MPV vs nSections", "steelblue")
plot_panel(axes[0,1], nsecs, xis,   xi_errs,
           "Landau ξ width [GeV]", "Width (ξ) vs nSections", "crimson")

# Chi2/ndf — no error bars
axes[1,0].plot(nsecs, chi2ndfs, "o-", color="seagreen",
               linewidth=2, markersize=8)
axes[1,0].axhline(1.0, color="grey", linestyle="--", label="χ²/ndf = 1 (ideal)")
axes[1,0].set_xlabel("nSections", fontsize=12)
axes[1,0].set_ylabel("χ²/ndf", fontsize=12)
axes[1,0].set_title("Fit quality vs nSections", fontsize=12)
axes[1,0].legend(fontsize=10); axes[1,0].grid(True, alpha=0.3)

# Mean nHits per event
axes[1,1].plot(nsecs, nhits, "o-", color="darkorange",
               linewidth=2, markersize=8)
# Overlay linear expectation from nSec=2 point
if len(nsecs) > 1:
    slope  = nhits[0] / nsecs[0]
    x_line = np.linspace(nsecs[0], nsecs[-1], 100)
    axes[1,1].plot(x_line, slope * x_line, "k--", linewidth=1.5,
                   label="Linear expectation")
axes[1,1].set_xlabel("nSections", fontsize=12)
axes[1,1].set_ylabel("Mean hits / event", fontsize=12)
axes[1,1].set_title("Mean nHits vs nSections", fontsize=12)
axes[1,1].legend(fontsize=10); axes[1,1].grid(True, alpha=0.3)

fig.tight_layout()
outpath = outdir / "segmentation_study.png"
fig.savefig(outpath, dpi=150, bbox_inches="tight")
plt.close(fig)
print(f"\n  Saved: {outpath}")
print(f"  Saved: {outdir}/overlay_all_fits.png")

# ── Text summary ──────────────────────────────────────────────────────────────
summary = outdir / "segmentation_summary.txt"
with open(summary, "w") as f:
    sep = "=" * 85
    f.write(sep + "\n")
    f.write("ProtoTPC Segmentation Study — Landau Fit Parameter Summary\n")
    f.write(sep + "\n")
    hdr = (f"{'nSec':>5}  {'MPV':>10}  {'MPV_err':>9}  {'xi':>10}  "
           f"{'xi_err':>9}  {'chi2/ndf':>9}  {'nHits/evt':>10}")
    f.write(hdr + "\n" + "-"*85 + "\n")
    for r in results:
        f.write(f"{r['nsec']:>5}  {r['mpv']:>10.6f}  {r['mpv_err']:>9.6f}  "
                f"{r['xi']:>10.6f}  {r['xi_err']:>9.6f}  "
                f"{r['chi2ndf']:>9.3f}  {r['nhits']:>10.2f}\n")
    f.write("\n" + sep + "\n")
    f.write("Stability (should be flat if segmentation has no physics effect):\n")
    for label, vals, errs in [
        ("MPV [GeV]",   mpvs, mpv_errs),
        ("Width xi [GeV]", xis, xi_errs)]:
        m = np.mean(vals); s = np.std(vals)
        vp = s/m*100 if m != 0 else 0
        f.write(f"  {label:<20}: mean={m:.6f}  std={s:.6f}  "
                f"var={vp:.2f}%  "
                f"[{'STABLE <5%' if vp<5 else 'VARIES >5%'}]\n")
print(f"  Saved: {summary}\n")