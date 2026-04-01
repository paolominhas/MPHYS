#!/usr/bin/env python3
"""
plot_edep.py
Fits per-hit ProtoTPC Edep with a Moyal (Landau) approximation + residuals.
One fit plot per segmentation, plus an overlay of all four.
"""

import sys, re
from pathlib import Path
import numpy as np
import uproot
import awkward as ak
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import mplhep as hep
from scipy.optimize import curve_fit
from scipy.stats import moyal

INDIR  = Path(sys.argv[1]) if len(sys.argv) > 1 else Path("./ana_output_muon")
OUTDIR = Path("edep_plots"); OUTDIR.mkdir(exist_ok=True)
TRUNC  = 0.70

colors  = ["steelblue", "crimson", "seagreen", "darkorange"]
markers = ["o", "s", "^", "D"]

def landau_approx(x, amp, loc, scale):
    return amp * moyal.pdf(x, loc=loc, scale=scale)

# ── Main loop ─────────────────────────────────────────────────────────────────
plt.style.use(hep.style.CMS)

# Collect everything first so the overlay has all 4
all_data   = []   # list of dicts per file
all_counts = []
all_edges  = []
all_popts  = []
all_labels = []

input_files = sorted(INDIR.glob("ana_geom_*.root"),
                     key=lambda p: int(re.search(r'[nN][sS]ec0*(\d+)', p.name).group(1)))

for idx, fpath in enumerate(input_files):
    m    = re.search(r'[nN][sS]ec0*(\d+)', fpath.name)
    nsec = int(m.group(1)) if m else idx

    print(f"\n  Processing: {fpath.name}  (nSections={nsec})")

    tree        = uproot.open(f"{fpath}:hibeam")
    edep_jagged = tree["ProtoTPC/Edep"].array()
    edep_flat   = ak.to_numpy(ak.flatten(edep_jagged))
    edep_hits   = edep_flat[edep_flat > 0]

    if len(edep_hits) == 0:
        print(f"    No hits — skipping.")
        continue

    # Truncation cut
    cut_value = np.percentile(edep_hits, TRUNC * 100)
    e_cut     = edep_hits[edep_hits <= cut_value]

    counts, edges = np.histogram(e_cut, bins=50)
    bin_centers   = (edges[:-1] + edges[1:]) / 2
    errors        = np.sqrt(counts)
    errors[errors == 0] = 1

    # Fit — use data-driven initial guess and wide bounds
    peak_idx = np.argmax(counts)
    p0 = [float(counts[peak_idx]),
          float(bin_centers[peak_idx]),
          float((cut_value - edep_hits.min()) * 0.05)]
    try:
        popt, _ = curve_fit(
            landau_approx, bin_centers, counts,
            p0=p0, sigma=errors,
            bounds=([0,        edep_hits.min(), 1e-9],
                    [np.inf,   cut_value,       cut_value]),
            maxfev=20000)
        print(f"    MPV={popt[1]:.5f} MeV  width={popt[2]:.5f} MeV  "
              f"entries={len(edep_hits)}  cut={cut_value:.5f} MeV")
    except RuntimeError:
        print(f"    Fit did not converge — using p0 estimate.")
        popt = p0

    all_data.append({
        "nsec": nsec, "idx": idx,
        "edep_hits": edep_hits, "e_cut": e_cut,
        "counts": counts, "edges": edges, "bin_centers": bin_centers,
        "errors": errors, "popt": popt, "cut_value": cut_value,
    })

# ── Per-file fit + residual plots ─────────────────────────────────────────────
for d in all_data:
    nsec        = d["nsec"]
    col         = colors[d["idx"] % len(colors)]
    popt        = d["popt"]
    edges       = d["edges"]
    counts      = d["counts"]
    bin_centers = d["bin_centers"]
    errors      = d["errors"]
    cut_value   = d["cut_value"]

    x_smooth   = np.linspace(edges[0], edges[-1], 500)
    y_smooth   = landau_approx(x_smooth, *popt)
    residuals  = counts - landau_approx(bin_centers, *popt)

    fig, (ax_main, ax_resid) = plt.subplots(
        2, 1, sharex=True,
        gridspec_kw={"height_ratios": [3, 1], "hspace": 0},
        figsize=(10, 8))
    fig.suptitle(f"ProtoTPC Per-hit Edep  |  nSections = {nsec}",
                 fontweight="bold", fontsize=14)

    # Main panel
    hep.histplot(counts, bins=edges, yerr=errors, ax=ax_main,
                 histtype="errorbar", color="black",
                 label=f"Per-hit Edep  (70% cut,  n={len(d['e_cut'])})")
    ax_main.plot(x_smooth, y_smooth, color=col, linewidth=2.5,
                 label=f"Landau (Moyal) fit\nMPV = {popt[1]:.5f} MeV")
    ax_main.axvline(cut_value, color="grey", linestyle=":",
                    linewidth=1.5, label=f"30% cut ({cut_value:.4f} MeV)")
    ax_main.set_ylabel("Hits / bin", fontsize=12)
    ax_main.legend(loc="upper right", fontsize=10)
    label_artists = hep.cms.label(loc=0, data=False, label="Simulation", ax=ax_main)
    label_artists[0].set_text("ESS")
    main_text.set_text("ESS")

    # Residuals panel
    ax_resid.errorbar(bin_centers, residuals, yerr=errors,
                      fmt="o", color="black", markersize=4)
    ax_resid.axhline(0, color="grey", linestyle="--")
    ax_resid.set_xlabel("Per-hit Edep [MeV]", fontsize=12)
    ax_resid.set_ylabel("Data − Fit", fontsize=11)
    lim = max(np.max(np.abs(residuals)), 1) * 1.5
    ax_resid.set_ylim(-lim, lim)

    fig.savefig(OUTDIR / f"edep_fit_nSec{nsec:03d}.png", dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"    Saved: edep_fit_nSec{nsec:03d}.png")

# ── Overlay — all segmentations ───────────────────────────────────────────────
fig_ov, ax_ov = plt.subplots(figsize=(11, 7))

for d in all_data:
    nsec        = d["nsec"]
    col         = colors[d["idx"] % len(colors)]
    popt        = d["popt"]
    edges       = d["edges"]
    counts      = d["counts"]
    bin_centers = d["bin_centers"]
    errors      = d["errors"]
    x_smooth    = np.linspace(edges[0], edges[-1], 500)
    y_smooth    = landau_approx(x_smooth, *popt)

    ax_ov.step(bin_centers, counts, where="mid",
           color=col, linewidth=2, alpha=0.8,
           label=f"nSections={nsec}  MPV={popt[1]:.4f} MeV")
    ax_ov.errorbar(bin_centers, counts, yerr=errors,
                   fmt="none", color=col, alpha=0.4, capsize=2)
    ax_ov.plot(x_smooth, y_smooth, color=col, linewidth=2,
           linestyle="--")

ax_ov.set_xlabel("Per-hit Edep [MeV]", fontsize=13)
ax_ov.set_ylabel("Hits / bin", fontsize=13)
ax_ov.set_title("ProtoTPC Per-hit Edep — all segmentations\n"
                "(solid = data, dashed = Landau fit)",
                fontsize=13, fontweight="bold")
ax_ov.legend(fontsize=11)
ax_ov.grid(True, alpha=0.3)
fig_ov.tight_layout()
fig_ov.savefig(OUTDIR / "overlay_edep.png", dpi=150, bbox_inches="tight")
plt.close(fig_ov)
print(f"\n  Overlay saved: overlay_edep.png")

# ── MPV vs nSections ──────────────────────────────────────────────────────────
if len(all_data) > 1:
    nsecs  = [d["nsec"]      for d in all_data]
    mpvs   = [d["popt"][1]   for d in all_data]
    widths = [d["popt"][2]   for d in all_data]

    fig_m, ax_m = plt.subplots(figsize=(8, 5))
    ax_m.errorbar(nsecs, mpvs, yerr=widths, fmt="o-", color="steelblue",
                  capsize=5, linewidth=2, markersize=8, label="MPV ± Landau width")
    ax_m.set_xlabel("nSections", fontsize=13)
    ax_m.set_ylabel("Landau MPV [MeV]", fontsize=13)
    ax_m.set_title("Landau MPV vs Segmentation", fontweight="bold", fontsize=14)
    ax_m.legend(); ax_m.grid(True, alpha=0.3)
    fig_m.tight_layout()
    fig_m.savefig(OUTDIR / "mpv_vs_nsections.png", dpi=150, bbox_inches="tight")
    plt.close(fig_m)

    print(f"\n  MPV summary:")
    for d in all_data:
        print(f"    nSec={d['nsec']:>3}  MPV={d['popt'][1]:.5f} MeV"
              f"  width={d['popt'][2]:.5f} MeV")

print(f"\n  All plots saved to {OUTDIR}/")