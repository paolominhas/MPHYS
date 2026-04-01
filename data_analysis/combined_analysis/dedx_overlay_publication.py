#!/usr/bin/env python3
"""
dedx_overlay_publication.py
============================
Produces 4 publication-quality dE/dx overlay figures:

  Krakow scattering:
    (1) run_0006  vs  KrakowScatter.root / 250proton.root
    (2) run_0042  vs  KrakowScatter.root / 250proton.root

  Muon lab:
    (3) run_0006  vs  MuonScatter_fixed.root
    (4) run_0042  vs  MuonScatter_fixed.root

Normalisation strategy
----------------------
Both distributions are normalised to a dimensionless axis x/MPV, where
MPV is located from the peak of a coarse histogram before fitting.
This makes both datasets sit on the same x-axis with no assumptions
about absolute calibration or unit conversion (ADC/mm vs MeV).

Data loading
------------
Experimental:  PyROOT + compiled TrackData.so — required because
  the 'tracks' branch holds std::vector<TrackData> objects that uproot
  cannot deserialise without the dictionary.  dE/dx is computed as
  ADC-per-centroid divided by the per-step path length (total 3D track
  length / number of centroids), exactly as in combining_dedx.py.

Simulation:    uproot — simple ProtoTPC/Edep summed per event, GeV→MeV.

Track quality cuts [literature]
--------------------------------
  chi²/ndf < 25  (computed from centroid residuals to the fitted line,
                  both x–y and z–y projections; ndf = 2n − 4)
                  Looser than ALICE (chi²/ndf < 4) because the prototype
                  TPC has fewer centroids per track than a full-scale TPC
                  [ALICE JINST 5 (2010) P09002; STAR NIM A 499 (2003) 659]
  nPoints >= 3   (minimum centroids)
  length  > 0    (rejects point-like clusters)

Fit
---
Moyal (Landau approximation) fitted to the bottom 70% of each
distribution independently on the dimensionless axis.
Residuals shown as (data − fit) / sqrt(|fit|) in a pull panel.

References
----------
Landau (1944); Moyal (1955); Bichsel Rev.Mod.Phys. 60 (1988) 663;
ALICE JINST 5 (2010) P09002; PDG Statistics (2022).

Usage
-----
    python3 dedx_overlay_publication.py
    Edit PAIRS and PATHS at the top for your file locations.
"""

import os
import sys
import warnings
from pathlib import Path

import ROOT
import awkward as ak
import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import mplhep as hep
import numpy as np
import uproot
from scipy import optimize, stats as scipy_stats

# ══════════════════════════════════════════════════════════════════════════════
# CONFIG — edit these paths
# ══════════════════════════════════════════════════════════════════════════════

# Directory containing recovered_headers/ subfolder with TrackData.so
HEADERS_DIR = Path("recovered_headers")

PAIRS = [
    {
        "title"   : r"Krakow $^2$H(p,p) elastic  —  run 0006",
        "exp_file": "../experimental_data/tracks_centroids_tpc_run_0006-sorted_t0_100.root",
        "sim_file": "../simulation_data/KrakowScatter.root",
        "output"  : "krakow_run0006",
        "x_label" : r"dE/dx  [units of peak MPV]",
    },
    {
        "title"   : r"Krakow $^2$H(p,p) elastic  —  run 0042",
        "exp_file": "../experimental_data/tracks_centroids_tpc_run_0042-sorted_t0.root",
        "sim_file": "../simulation_data/KrakowScatter.root",
        "output"  : "krakow_run0042",
        "x_label" : r"dE/dx  [units of peak MPV]",
    },
    {
        "title"   : "Muon lab (MCPL)  —  run 0006",
        "exp_file": "/home/s2289940/ESS/HIBEAM/data_analysis/experimental_data/tracks_centroids_chamber-fec01-debug-clk0-ped0-thr0-g0-s0-p1-1008-2_sorted.root",
        "sim_file": "../simulation_data/MuonScatter_fixed.root",
        "output"  : "muon_run1",
        "x_label" : r"dE/dx  [units of peak MPV]",
    },
    {
        "title"   : "Muon lab (MCPL)  —  run 0042",
        "exp_file": "/home/s2289940/ESS/HIBEAM/data_analysis/experimental_data/tracks_centroids_chamber-fec01-debug-clk0-ped0-thr0-g0-s0-p1-1008-4_sorted.root",
        "sim_file": "../simulation_data/MuonScatter_fixed.root",
        "output"  : "muon_run2",
        "x_label" : r"dE/dx  [units of peak MPV]",
    },
]

# Analysis parameters
TRUNCATION    = 0.70    # fit bottom 70% [ALICE JINST 5 (2010) P09002]
N_BINS        = 35      # coarse enough for ≥10 events/bin
CHI2_NDF_MAX  = 25.0    # track quality — looser than ALICE due to short tracks
MIN_POINTS    = 3       # minimum centroids per track
SIM_MIN_STEPS = 3       # minimum G4 steps per sim event

# ══════════════════════════════════════════════════════════════════════════════
# SECTION 1 — PyROOT setup  (done once at startup)
# ══════════════════════════════════════════════════════════════════════════════

def setup_pyroot():
    headers_so  = HEADERS_DIR / "recovered_headers.so"
    hdr_file    = HEADERS_DIR / "TrackData.h"

    if not headers_so.exists():
        sys.exit(f"ERROR: {headers_so} not found.\n"
                 "Run TFile::MakeProject to regenerate.")

    ROOT.gInterpreter.AddIncludePath(str(HEADERS_DIR))
    if hdr_file.exists():
        ROOT.gInterpreter.ProcessLine(f'#include "{hdr_file}"')
    ROOT.gSystem.Load(str(headers_so))
    print(f"  PyROOT: loaded {headers_so.name}")

# ══════════════════════════════════════════════════════════════════════════════
# SECTION 2 — EXPERIMENTAL DATA LOADING
# Per-centroid dE/dx:  adc[j] / (track_length / n_centroids)
# chi² recomputed from stored centroid residuals (both projections)
# ══════════════════════════════════════════════════════════════════════════════

def load_experimental(filepath):
    """
    Returns flat array of per-centroid dE/dx values [ADC/mm].

    Strategy: for each accepted track, the 3D length of the fitted line
    is divided equally among n centroids (delta_x = length/n), giving a
    uniform path length per centroid.  Each centroid's ADC is then divided
    by delta_x.  This is the most robust estimator given the short tracks
    typical of this prototype TPC.
    """
    print(f"  Loading exp: {Path(filepath).name}")
    tfile = ROOT.TFile.Open(str(filepath))
    if not tfile or tfile.IsZombie():
        raise IOError(f"Cannot open {filepath}")

    tree = tfile.Get("trackingData")
    if not tree:
        tfile.Close()
        raise KeyError("'trackingData' not found")

    tracks_vec = ROOT.std.vector("TrackData")()
    tree.SetBranchAddress("tracks", tracks_vec)

    dedx = []
    n_total = n_few = n_chi2 = n_short = 0

    for i in range(tree.GetEntries()):
        tree.GetEntry(i)
        for tr in tracks_vec:
            n = tr.nPoints
            n_total += 1

            # Cut 1: minimum centroids
            if n < MIN_POINTS:
                n_few += 1
                continue

            # Cut 2: chi²/ndf from centroid residuals to the fitted line
            # Both x–y and z–y projections contribute to chi²
            chi2 = 0.0
            for j in range(n):
                y_j = float(tr.y[j])
                sx  = max(float(tr.sigmas_x[j]), 1e-4)
                sz  = max(float(tr.sigmas_z[j]), 1e-4)
                res_x = float(tr.x[j]) - (tr.slope_xy * y_j + tr.intercept_xy)
                res_z = float(tr.z[j]) - (tr.slope_zy * y_j + tr.intercept_zy)
                chi2 += (res_x / sx)**2 + (res_z / sz)**2
            ndf = max(2 * n - 4, 1)
            if chi2 / ndf > CHI2_NDF_MAX:
                n_chi2 += 1
                continue

            # 3D track length from first to last centroid
            x0, x1 = float(tr.x[0]),   float(tr.x[n-1])
            y0, y1 = float(tr.y[0]),   float(tr.y[n-1])
            z0, z1 = float(tr.z[0]),   float(tr.z[n-1])
            length = float(np.sqrt((x1-x0)**2 + (y1-y0)**2 + (z1-z0)**2))
            if length < 1e-3:
                n_short += 1
                continue

            delta_x = length / n  # path length per centroid [mm]
            for j in range(n):
                adc = float(tr.charge[j])
                if adc > 0:
                    dedx.append(adc / delta_x)

    tfile.Close()
    arr = np.array(dedx, dtype=float)
    print(f"    {n_total} tracks → accepted {n_total-n_few-n_chi2-n_short}  "
          f"(few={n_few}, chi2={n_chi2}, short={n_short})  "
          f"→ {len(arr)} centroid dE/dx values")

    if len(np.unique(arr)) < 20:
        raise RuntimeError(
            "Too few unique dE/dx values — check TrackData.so matches the file.")

    # Remove extreme tail for display (top 5%) — does not affect fit range
    return arr[arr < np.percentile(arr, 95)]


# ══════════════════════════════════════════════════════════════════════════════
# SECTION 3 — SIMULATION DATA LOADING
# Sum Edep per event, convert GeV→MeV.
# ══════════════════════════════════════════════════════════════════════════════

def load_simulation(filepath):
    print(f"  Loading sim: {Path(filepath).name}")
    with uproot.open(str(filepath)) as f:
        tree  = f["hibeam"]
        bkeys = set(tree.keys(recursive=True))

        edep_key = next((k for k in
            ["ProtoTPC/Edep", "ProtoTPC.Edep", "TPC/Edep", "TPC.Edep"]
            if k in bkeys), None)
        if edep_key is None:
            raise KeyError(f"No Edep branch found. Keys: {list(bkeys)[:15]}")

        edep_jag = tree[edep_key].array(library="ak")

        weights = np.ones(len(edep_jag))
        if "PrimaryWeight" in bkeys:
            w_jag   = tree["PrimaryWeight"].array(library="ak")
            weights = ak.to_numpy(
                ak.fill_none(ak.firsts(w_jag), 1.0)).astype(float)

    n_steps = ak.to_numpy(ak.num(edep_jag)).astype(int)
    e_sum   = ak.to_numpy(ak.sum(edep_jag, axis=1)).astype(float) * 1000.0  # GeV→MeV
    mask    = (n_steps >= SIM_MIN_STEPS) & (e_sum > 0)
    vals    = e_sum[mask]
    wts     = weights[mask]
    print(f"    {len(e_sum)} events → {mask.sum()} accepted  "
          f"(range {vals.min():.2f}–{vals.max():.2f} MeV)")
    return vals, wts


# ══════════════════════════════════════════════════════════════════════════════
# SECTION 4 — PEAK LOCATION  (coarse, for normalisation)
# ══════════════════════════════════════════════════════════════════════════════

def find_peak(values, weights=None):
    """Locate histogram mode using 200 coarse bins."""
    counts, edges = np.histogram(values, bins=200, weights=weights)
    centers = 0.5 * (edges[:-1] + edges[1:])
    return float(centers[np.argmax(counts)])


# ══════════════════════════════════════════════════════════════════════════════
# SECTION 5 — HISTOGRAM ON SHARED DIMENSIONLESS AXIS
# ══════════════════════════════════════════════════════════════════════════════

def make_normalised_histogram(values, weights, edges):
    """
    Histogram on given edges, normalised to unit area (density=True).
    Errors: sqrt(sum(w²)) per bin, normalised by same factor as counts.
    """
    counts, _ = np.histogram(values, bins=edges, weights=weights)
    sumw2, _  = np.histogram(values, bins=edges,
                             weights=weights**2 if weights is not None
                             else np.ones(len(values)))
    bw   = float(edges[1] - edges[0])
    norm = float(counts.sum()) * bw
    if norm == 0:
        norm = 1.0
    density = counts / norm
    errors  = np.maximum(np.sqrt(sumw2) / norm, 1e-9)
    return density, errors


# ══════════════════════════════════════════════════════════════════════════════
# SECTION 6 — LANDAU (MOYAL) FIT
# ══════════════════════════════════════════════════════════════════════════════

def _moyal(x, loc, scale, amp):
    return amp * scipy_stats.moyal.pdf(x, loc=loc, scale=scale)


def fit_landau(centers, counts, errors, x_cut, label=""):
    """Fit Moyal to bins with centre ≤ x_cut (the truncation boundary)."""
    fm   = centers <= x_cut
    xf   = centers[fm]
    yf   = counts[fm]
    ye   = errors[fm]

    if fm.sum() < 4:
        warnings.warn(f"[{label}] too few bins for fit")
        return None

    loc0   = float(xf[np.argmax(yf)])
    scale0 = max(float((x_cut - xf[0]) / 8.0), 1e-9)
    amp0   = float(yf.max())

    try:
        popt, pcov = optimize.curve_fit(
            _moyal, xf, yf,
            p0=[loc0, scale0, amp0],
            sigma=ye, absolute_sigma=True,
            bounds=([xf[0], 1e-9, 0.0],
                    [x_cut,  x_cut, 100.0]),
            maxfev=40_000)
    except Exception as exc:
        warnings.warn(f"[{label}] fit failed: {exc}")
        return None

    perr     = np.sqrt(np.diag(pcov))
    y_pred   = _moyal(xf, *popt)
    chi2     = float(np.sum(((yf - y_pred) / ye)**2))
    ndf      = max(fm.sum() - 3, 1)
    chi2_red = chi2 / ndf
    p_val    = float(scipy_stats.chi2.sf(chi2, ndf))

    print(f"    [{label}] MPV={popt[0]:.4f}  xi={popt[1]:.4f}±{perr[1]:.4f}"
          f"  chi²/ndf={chi2_red:.2f}  p={p_val:.3f}")

    return {"popt": popt, "perr": perr,
            "mpv": popt[0], "xi": popt[1], "xi_err": perr[1],
            "chi2_red": chi2_red, "p_value": p_val}


# ══════════════════════════════════════════════════════════════════════════════
# SECTION 7 — PUBLICATION PLOT
# ══════════════════════════════════════════════════════════════════════════════

plt.style.use(hep.style.CMS)
plt.rcParams.update({
    "font.size"        : 13,
    "axes.labelsize"   : 14,
    "xtick.labelsize"  : 11,
    "ytick.labelsize"  : 11,
    "legend.fontsize"  : 10,
    "xtick.direction"  : "in",
    "ytick.direction"  : "in",
    "xtick.top"        : True,
    "ytick.right"      : True,
})

C_DATA = "#000000"   # black  — experimental data
C_SIM  = "#d73027"   # red    — simulation
C_TAIL = "#aaaaaa"   # grey   — excluded tail region


def make_plot(pair, exp_vals, exp_peak, d_counts, d_errors,
              sim_counts, sim_errors, d_fit, s_fit,
              centers, edges, x_cut, sim_peak):

    trunc_pct = int(TRUNCATION * 100)
    fm        = centers <= x_cut
    x_smooth  = np.linspace(float(centers[0]),
                             float(centers[fm][-1]), 2000)
    # CHANGED: full-range x for sim extrapolation
    x_full    = np.linspace(float(centers[0]),
                             float(centers[-1]), 2000)

    fig = plt.figure(figsize=(8.5, 7.5))
    gs  = gridspec.GridSpec(2, 1, height_ratios=[3, 1], hspace=0.0,
                            left=0.12, right=0.97, top=0.91, bottom=0.11)
    ax  = fig.add_subplot(gs[0])
    axr = fig.add_subplot(gs[1], sharex=ax)

    # ── Delta-ray tail region shading ─────────────────────────────────────────
    ax.axvspan(x_cut, float(edges[-1]), alpha=0.07, color=C_TAIL, zorder=0,
               label=f"Excluded tail  (top {100-trunc_pct}%,  $\\delta$-ray enhanced)")
    axr.axvspan(x_cut, float(edges[-1]), alpha=0.07, color=C_TAIL, zorder=0)

    # ── Histograms — data as errorbar, simulation as filled ───────────────────
    hep.histplot(sim_counts, bins=edges,
                 yerr=sim_errors, ax=ax,
                 histtype="fill", color=C_SIM, alpha=0.40,
                 label="Simulation (weighted, normalised)")

    hep.histplot(d_counts, bins=edges,
                 yerr=d_errors, ax=ax,
                 histtype="errorbar", color=C_DATA,
                 markersize=4, linewidth=1.2,
                 label="Data (normalised)")

    # ── Fit curves ────────────────────────────────────────────────────────────
    if d_fit is not None:
        ax.plot(x_smooth, _moyal(x_smooth, *d_fit["popt"]),
                color=C_DATA, lw=2.0, ls="-",
                label=f"Data Landau fit  MPV = {d_fit['mpv']:.3f}")
    # CHANGED: sim fit — solid in fitted region, dashed extrapolation
    if s_fit is not None:
        ax.plot(x_smooth, _moyal(x_smooth, *s_fit["popt"]),
                color=C_SIM,  lw=2.0, ls="-",
                label=f"Sim Landau fit   MPV = {s_fit['mpv']:.3f}")
        x_extrap = x_full[x_full > float(centers[fm][-1])]
        if len(x_extrap) > 0:
            ax.plot(x_extrap, _moyal(x_extrap, *s_fit["popt"]),
                    color=C_SIM, lw=1.8, ls="--",
                    label="Sim fit (extrapolated)")

    # Truncation line
    ax.axvline(x_cut, color="#4d9221", lw=1.2, ls=":",
               label=f"Truncation at {trunc_pct}th percentile")

    # CHANGED: annotation — values only, no ± errors
    d_mpv_raw = d_fit["mpv"] * exp_peak if d_fit else float("nan")
    s_mpv_raw = s_fit["mpv"] * sim_peak if s_fit else float("nan")
    ann_lines = []
    if d_fit is not None:
        ann_lines.append(f"Data:  MPV = {d_fit['mpv']:.3f}"
                         f"  $\\xi$ = {d_fit['xi']:.3f}"
                         f"  $\\chi^2$/ndf = {d_fit['chi2_red']:.2f}")
    if s_fit is not None:
        ann_lines.append(f"Sim:   MPV = {s_fit['mpv']:.3f}"
                         f"  $\\xi$ = {s_fit['xi']:.3f}"
                         f"  $\\chi^2$/ndf = {s_fit['chi2_red']:.2f}")
    ann_lines.append(f"Data MPV = {d_mpv_raw:.1f} ADC/mm")
    ann_lines.append(f"Sim  MPV = {s_mpv_raw:.3f} MeV")
    if d_mpv_raw and s_mpv_raw and d_mpv_raw != 0:
        ann_lines.append(f"Ratio = {s_mpv_raw/d_mpv_raw:.3e} MeV/(ADC/mm)")
    ann = "\n".join(ann_lines)
    ax.text(0.97, 0.50, ann, transform=ax.transAxes,
            va="top", ha="right", fontsize=8.5, family="monospace",
            bbox=dict(boxstyle="round,pad=0.35", facecolor="white",
                      edgecolor="#aaaaaa", alpha=0.93))

    ax.set_ylabel("Probability density  (unit area)")
    ax.set_title(pair["title"], pad=4, fontsize=13)
    ax.set_yscale("log")
    ax.set_ylim(bottom=1e-4)
    ax.legend(loc="upper right", fontsize=9,
              framealpha=0.92, edgecolor="#aaaaaa",
              handlelength=1.8, labelspacing=0.35)
    plt.setp(ax.get_xticklabels(), visible=False)

    # ESS branding
    try:
        lo = hep.cms.label(loc=0, data=True, label="Preliminary",
                           rlabel="HIBEAM", ax=ax)
        lo[0].set_text("ESS")
    except Exception:
        pass

    # ── Pull residuals: (data − fit) / sqrt|fit|  (fitted region only) ───────
    all_pulls = []
    for counts_arr, fit_res, color, marker in [
            (d_counts, d_fit, C_DATA, "o"),
            (sim_counts, s_fit, C_SIM, "s")]:
        if fit_res is None:
            continue
        yf    = _moyal(centers, *fit_res["popt"])
        denom = np.sqrt(np.maximum(np.abs(yf), 1e-10))
        pulls = (counts_arr - yf) / denom
        p_err = np.maximum(counts_arr, 1e-9)**0.5 / denom
        axr.errorbar(centers[fm], pulls[fm], yerr=p_err[fm],
                     fmt=marker, color=color, markersize=4,
                     linewidth=1.0, capsize=2)
        all_pulls.extend(pulls[fm].tolist())

    axr.axhline(0,  color="black",   lw=0.9)
    axr.axhline(+2, color="#888888", lw=0.9, ls="--",
                label=r"$\pm2\sigma$")
    axr.axhline(-2, color="#888888", lw=0.9, ls="--")
    axr.axvline(x_cut, color="#4d9221", lw=1.2, ls=":", alpha=0.7)

    pull_max = max(float(np.abs(all_pulls).max()) if all_pulls else 3.0, 2.5)
    ylim     = float(np.ceil(pull_max / 2.0) * 2.0) + 0.5
    axr.set_ylim(-ylim, ylim)
    axr.yaxis.set_major_locator(ticker.MultipleLocator(2.0))
    axr.legend(fontsize=9, loc="upper right",
               framealpha=0.92, edgecolor="#aaaaaa")
    axr.set_ylabel(r"$(N-F)/\!\sqrt{|F|}$", fontsize=11)
    axr.set_xlabel(pair["x_label"])

    for ext in [".pdf", ".png"]:
        path = Path(pair["output"] + ext)
        fig.savefig(path, dpi=300 if ext==".png" else None,
                    bbox_inches="tight")
        print(f"    Saved: {path}")
    plt.close(fig)


# ══════════════════════════════════════════════════════════════════════════════
# MAIN
# ══════════════════════════════════════════════════════════════════════════════

setup_pyroot()

for pair in PAIRS:
    print(f"\n{'='*60}")
    print(f"  {pair['title']}")
    print(f"{'='*60}")

    try:
        exp_vals = load_experimental(pair["exp_file"])
    except Exception as e:
        print(f"  SKIP — exp failed: {e}"); continue

    try:
        sim_vals, sim_wts = load_simulation(pair["sim_file"])
    except Exception as e:
        print(f"  SKIP — sim failed: {e}"); continue

    if len(exp_vals) < 100 or len(sim_vals) < 100:
        print(f"  SKIP — too few events"); continue

    # ── Peak-normalise to dimensionless axis x/MPV ────────────────────────────
    exp_peak = find_peak(exp_vals)
    sim_peak = find_peak(sim_vals, sim_wts)
    print(f"  Peaks:  data={exp_peak:.3f} ADC/mm   sim={sim_peak:.4f} MeV")

    exp_norm = exp_vals / exp_peak
    sim_norm = sim_vals / sim_peak

    # ── Shared bin edges: 0 to 99th pct of BOTH distributions ────────────────
    x_max = min(np.percentile(exp_norm, 99.5), np.percentile(sim_norm, 99.5))
    edges = np.linspace(0.0, x_max, N_BINS + 1)
    cents = 0.5 * (edges[:-1] + edges[1:])

    # Truncation boundary in normalised units (same for both)
    x_cut = float(np.percentile(
        np.concatenate([exp_norm, sim_norm]), TRUNCATION * 100.0))

    # ── Histograms ────────────────────────────────────────────────────────────
    d_cnt, d_err = make_normalised_histogram(
        exp_norm, np.ones(len(exp_norm)), edges)
    s_cnt, s_err = make_normalised_histogram(
        sim_norm, sim_wts, edges)

    # ── Fits ──────────────────────────────────────────────────────────────────
    print("  Fitting Landau distributions...")
    d_fit = fit_landau(cents, d_cnt, d_err, x_cut, "Data")
    s_fit = fit_landau(cents, s_cnt, s_err, x_cut, "Simulation")

    # ── Plot ──────────────────────────────────────────────────────────────────
    print("  Plotting...")
    make_plot(pair, exp_vals, exp_peak, d_cnt, d_err,
              s_cnt, s_err, d_fit, s_fit,
              cents, edges, x_cut, sim_peak)

print("\n  All done.")