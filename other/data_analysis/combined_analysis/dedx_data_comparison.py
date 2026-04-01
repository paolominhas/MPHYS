#!/usr/bin/env python3
"""
dedx_data_sim_comparison.py
============================
Overlays experimental TPC dE/dx data with HIBEAM Geant4 simulation on a
common peak-normalised axis, following the same style and fitting principles
as sim_dedx_plots.py.

Four figures produced:
  (1) Krakow run 0006  vs  KrakowScatter.root
  (2) Krakow run 0042  vs  KrakowScatter.root
  (3) Muon   run 0006  vs  MuonScatter_fixed.root
  (4) Muon   run 0042  vs  MuonScatter_fixed.root

Normalisation
-------------
Both distributions divided by their own histogram peak position so the
x-axis is dimensionless (units of MPV).  This removes the ADC/mm ↔ MeV
unit mismatch entirely without any affine rescaling assumptions.

Experimental loading
--------------------
PyROOT + TrackData.so (must be in recovered_headers/).
Per-centroid dE/dx = ADC[j] / (track_length / nPoints).
Quality cut: chi²/ndf < 25  (recomputed from stored centroid residuals).

Simulation loading
------------------
uproot: sum ProtoTPC/Edep per event, GeV→MeV.
Muon simulation: apply E > 15 MeV cut to remove sub-threshold bump
(delta-rays, clipping tracks, Michel electrons — Bichsel 1988 Sec IV).
Krakow simulation: no low-energy cut needed.

Fit
---
Moyal (Landau approximation) fitted independently to each distribution
on the bottom 70% of its own data [ALICE JINST 5 (2010) P09002].
Both fits plotted on the shared dimensionless axis.
Pull residuals shown for both in the bottom panel.

Usage
-----
    python3 dedx_data_sim_comparison.py
    Run from the combined_analysis/ directory (needs recovered_headers/).
"""

import os, sys, warnings
from pathlib import Path

import ROOT
import awkward as ak
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.ticker as ticker
import mplhep as hep
import numpy as np
import uproot
from scipy import optimize, stats as sc

# ══════════════════════════════════════════════════════════════════════════════
# CONFIG
# ══════════════════════════════════════════════════════════════════════════════
HEADERS_DIR = Path("recovered_headers")

PAIRS = [
    {
        "title"   : r"Krakow $^2$H(p,p)  —  run 0006",
        "exp_file": "../experimental_data/tracks_centroids_tpc_run_0006-sorted_t0_100.root",
        "sim_file": "../simulation_data/KrakowScatter.root",
        "sim_low_mev": 0.0,     # no cut for Krakow
        "output"  : "comparison_krakow_run0006",
    },
    {
        "title"   : r"Krakow $^2$H(p,p)  —  run 0042",
        "exp_file": "../experimental_data/tracks_centroids_tpc_run_0042-sorted_t0.root",
        "sim_file": "../simulation_data/KrakowScatter.root",
        "sim_low_mev": 0.0,
        "output"  : "comparison_krakow_run0042",
    },
    {
        "title"   : "Muon lab (MCPL)  —  run 0006",
        "exp_file": "../experimental_data/tracks_centroids_tpc_run_0006-sorted_t0_100.root",
        "sim_file": "../simulation_data/MuonScatter_fixed.root",
        "sim_low_mev": 15.0,    # remove sub-threshold bump
        "output"  : "comparison_muon_run0006",
    },
    {
        "title"   : "Muon lab (MCPL)  —  run 0042",
        "exp_file": "../experimental_data/tracks_centroids_tpc_run_0042-sorted_t0.root",
        "sim_file": "../simulation_data/MuonScatter_fixed.root",
        "sim_low_mev": 15.0,
        "output"  : "comparison_muon_run0042",
    },
]

TRUNCATION  = 0.70
N_BINS      = 35
MIN_STEPS   = 3
CHI2_MAX    = 25.0
MIN_POINTS  = 3
OUTDIR      = "."

# ══════════════════════════════════════════════════════════════════════════════
# STYLE
# ══════════════════════════════════════════════════════════════════════════════
plt.style.use(hep.style.CMS)
plt.rcParams.update({
    "font.size": 12, "axes.labelsize": 13,
    "xtick.labelsize": 10, "ytick.labelsize": 10,
    "legend.fontsize": 8.5,
    "xtick.direction": "in", "ytick.direction": "in",
    "xtick.top"  : False,   # no ticks on top axis
    "ytick.right": False,   # no ticks on right axis
})

C_DATA  = "#1f77b4"   # blue   — experimental data (standard matplotlib blue)
C_SIM   = "#d62728"   # red    — simulation
C_BUMP  = "#fff7ec"   # cream  — muon sub-threshold bump background

# ══════════════════════════════════════════════════════════════════════════════
# PYROOT SETUP
# ══════════════════════════════════════════════════════════════════════════════
def setup_pyroot():
    so  = HEADERS_DIR / "recovered_headers.so"
    hdr = HEADERS_DIR / "TrackData.h"
    if not so.exists():
        sys.exit(f"ERROR: {so} not found. Re-run TFile::MakeProject.")
    ROOT.gInterpreter.AddIncludePath(str(HEADERS_DIR))
    if hdr.exists():
        ROOT.gInterpreter.ProcessLine(f'#include "{hdr}"')
    ROOT.gSystem.Load(str(so))
    print(f"  PyROOT: loaded {so.name}")

# ══════════════════════════════════════════════════════════════════════════════
# EXPERIMENTAL LOADING  (PyROOT + TrackData)
# ══════════════════════════════════════════════════════════════════════════════
def load_experimental(filepath):
    """
    Per-centroid dE/dx [ADC/mm] from trackingData tree.
    dE/dx = ADC[j] / (total_3D_length / nPoints)
    Quality: chi²/ndf < CHI2_MAX,  nPoints >= MIN_POINTS
    """
    print(f"  Exp: {Path(filepath).name}")
    tfile = ROOT.TFile.Open(str(filepath))
    if not tfile or tfile.IsZombie():
        raise IOError(f"Cannot open {filepath}")
    tree = tfile.Get("trackingData")
    if not tree:
        tfile.Close()
        raise KeyError("'trackingData' not found")

    vec = ROOT.std.vector("TrackData")()
    tree.SetBranchAddress("tracks", vec)

    dedx = []
    n_total = n_few = n_chi2 = n_short = 0

    for i in range(tree.GetEntries()):
        tree.GetEntry(i)
        for tr in vec:
            n = tr.nPoints
            n_total += 1
            if n < MIN_POINTS:
                n_few += 1; continue

            # chi²/ndf from centroid residuals (both x-y and z-y projections)
            chi2 = 0.0
            for j in range(n):
                y_j = float(tr.y[j])
                sx  = max(float(tr.sigmas_x[j]), 1e-4)
                sz  = max(float(tr.sigmas_z[j]), 1e-4)
                chi2 += ((float(tr.x[j]) - (tr.slope_xy * y_j + tr.intercept_xy)) / sx)**2
                chi2 += ((float(tr.z[j]) - (tr.slope_zy * y_j + tr.intercept_zy)) / sz)**2
            ndf = max(2 * n - 4, 1)
            if chi2 / ndf > CHI2_MAX:
                n_chi2 += 1; continue

            # 3D track length → path length per centroid
            x0, x1 = float(tr.x[0]),   float(tr.x[n-1])
            y0, y1 = float(tr.y[0]),   float(tr.y[n-1])
            z0, z1 = float(tr.z[0]),   float(tr.z[n-1])
            length = float(np.sqrt((x1-x0)**2 + (y1-y0)**2 + (z1-z0)**2))
            if length < 1e-3:
                n_short += 1; continue

            dx = length / n
            for j in range(n):
                adc = float(tr.charge[j])
                if adc > 0:
                    dedx.append(adc / dx)

    tfile.Close()
    arr = np.array(dedx, dtype=float)
    print(f"    {n_total} tracks → {n_total-n_few-n_chi2-n_short} accepted  "
          f"→ {len(arr)} centroid values")
    if len(np.unique(arr)) < 20:
        raise RuntimeError("Too few unique dE/dx values — check TrackData.so")
    # Remove extreme tail (top 5%) for display
    return arr[arr < np.percentile(arr, 95)]

# ══════════════════════════════════════════════════════════════════════════════
# SIMULATION LOADING  (uproot)
# ══════════════════════════════════════════════════════════════════════════════
def load_simulation(filepath, low_mev=0.0):
    """Sum ProtoTPC/Edep per event [MeV]. Apply low_mev cut if > 0."""
    print(f"  Sim: {Path(filepath).name}")
    with uproot.open(str(filepath)) as f:
        tree  = f["hibeam"]
        bkeys = set(tree.keys(recursive=True))
        edep  = tree["ProtoTPC/Edep"].array(library="ak")
        n_s   = ak.to_numpy(ak.num(edep)).astype(int)
        e_s   = ak.to_numpy(ak.sum(edep, axis=1)).astype(float) * 1000.
        try:
            w = ak.to_numpy(ak.fill_none(
                ak.firsts(tree["PrimaryWeight"].array(library="ak")), 1.0))
        except Exception:
            w = np.ones(len(e_s))

    mask = (n_s >= MIN_STEPS) & (e_s > max(low_mev, 0.0))
    vals, wts = e_s[mask], w[mask]
    if len(vals) == 0:
        raise RuntimeError(f"No sim events passed cuts (low_mev={low_mev})")
    print(f"    {len(vals):,} events  ({vals.min():.2f}–{vals.max():.1f} MeV)"
          + (f"  [E>{low_mev:.0f} MeV cut]" if low_mev > 0 else ""))
    return vals, wts

# ══════════════════════════════════════════════════════════════════════════════
# FITTING  (shared with sim_dedx_plots.py logic)
# ══════════════════════════════════════════════════════════════════════════════
def _moyal(x, loc, scale, norm):
    return norm * sc.moyal.pdf(x, loc=loc, scale=scale)

def fit_on_edges(values, weights, edges, peak, label=""):
    """Fit Moyal to bottom TRUNCATION fraction on shared edges."""
    bc  = 0.5 * (edges[:-1] + edges[1:])
    bw  = float(edges[1] - edges[0])
    cut = float(np.percentile(values, TRUNCATION * 100))
    fm  = bc <= cut

    counts, _ = np.histogram(values, bins=edges, weights=weights)
    yn = counts / peak
    ye = np.maximum(np.sqrt(counts), 1.0) / peak

    if fm.sum() < 4:
        print(f"  [{label}] too few bins"); return None

    xf, yf, ef = bc[fm], yn[fm], ye[fm]
    loc0   = float(xf[np.argmax(yf)])
    scale0 = max((cut - float(bc[0])) / 8.0, 1e-9)
    norm0  = peak * bw
    try:
        popt, pcov = optimize.curve_fit(
            _moyal, xf, yf, p0=[loc0, scale0, norm0],
            sigma=ef, absolute_sigma=True,
            bounds=([float(bc[0]), 1e-9, 0], [cut, cut, norm0*20]),
            maxfev=40_000)
    except Exception as exc:
        print(f"  [{label}] fit failed: {exc}"); return None

    perr     = np.sqrt(np.diag(pcov))
    chi2_red = float(np.sum(((yf - _moyal(xf, *popt)) / ef)**2)) / max(fm.sum()-3, 1)
    print(f"  [{label}] MPV={popt[0]:.4g}  xi={popt[1]:.4g}±{perr[1]:.3g}"
          f"  chi²/ndf={chi2_red:.2f}")
    return {"popt": popt, "perr": perr, "mpv": popt[0], "xi": popt[1],
            "xi_err": perr[1], "chi2_red": chi2_red, "cut": cut,
            "yn": yn, "ye": ye, "fm": fm, "bc": bc, "bw": bw}

# ══════════════════════════════════════════════════════════════════════════════
# PLOT
# ══════════════════════════════════════════════════════════════════════════════
def make_comparison_plot(pair, exp_vals, sim_vals, sim_wts):
    low_mev  = pair["sim_low_mev"]
    is_muon  = low_mev > 0

    # ── Peak of each distribution (coarse histogram mode) ─────────────────────
    def coarse_peak(v, w=None):
        c, e = np.histogram(v, bins=200, weights=w)
        return float(0.5*(e[:-1]+e[1:])[np.argmax(c)])

    exp_pk_raw = coarse_peak(exp_vals)
    sim_pk_raw = coarse_peak(sim_vals, sim_wts)
    print(f"  Peaks:  data={exp_pk_raw:.3g} ADC/mm   sim={sim_pk_raw:.4g} MeV")

    # ── Dimensionless normalisation ────────────────────────────────────────────
    exp_norm = exp_vals / exp_pk_raw
    sim_norm = sim_vals / sim_pk_raw

    # ── Shared edges ──────────────────────────────────────────────────────────
    x_max = min(np.percentile(exp_norm, 99.5),
                np.percentile(sim_norm, 99.5), 8.0)
    edges = np.linspace(0, x_max, N_BINS + 1)
    bc    = 0.5*(edges[:-1]+edges[1:]);  bw = float(edges[1]-edges[0])

    # ── Histograms on shared edges ────────────────────────────────────────────
    exp_cnt, _ = np.histogram(exp_norm, bins=edges)
    sim_cnt, _ = np.histogram(sim_norm, bins=edges, weights=sim_wts)
    exp_pk_n   = float(exp_cnt.max());  sim_pk_n = float(sim_cnt.max())
    exp_yn = exp_cnt / exp_pk_n;  exp_ye = np.maximum(np.sqrt(exp_cnt),1.)/exp_pk_n
    sim_yn = sim_cnt / sim_pk_n;  sim_ye = np.maximum(np.sqrt(sim_cnt),1.)/sim_pk_n

    bump_edge = low_mev / sim_pk_raw if is_muon else 0.0
    bm = bc <= bump_edge

    # ── Fits ──────────────────────────────────────────────────────────────────
    print("  Fitting data...");  r_d = fit_on_edges(
        exp_norm, np.ones(len(exp_norm)), edges, exp_pk_n, "Data")
    print("  Fitting simulation...");  r_s = fit_on_edges(
        sim_norm, sim_wts, edges, sim_pk_n, "Simulation")

    cut_d = r_d["cut"] if r_d else x_max*TRUNCATION
    cut_s = r_s["cut"] if r_s else x_max*TRUNCATION
    cut_show = 0.5*(cut_d + cut_s)   # midpoint — single truncation line
    xs    = np.linspace(0, x_max, 3000)

    # Histogram region masks
    d_fit  = bc <= cut_d;          d_tail = ~d_fit
    s_fit  = (~bm) & (bc <= cut_s); s_tail = (~bm) & (~s_fit)

    # ── Figure ────────────────────────────────────────────────────────────────
    fig = plt.figure(figsize=(8.5, 6.8))
    gs  = gridspec.GridSpec(2, 1, height_ratios=[3,1], hspace=0.,
                            left=0.12, right=0.97, top=0.91, bottom=0.11)
    ax  = fig.add_subplot(gs[0])
    axr = fig.add_subplot(gs[1], sharex=ax)

    # Muon bump background
    if is_muon and bump_edge > 0:
        ax.axvspan(0, bump_edge, color=C_BUMP, alpha=1.0, zorder=0)
        axr.axvspan(0, bump_edge, color=C_BUMP, alpha=1.0, zorder=0)
        ax.axvline(bump_edge, color="#d6604d", lw=0.8, ls="-.", alpha=0.6, zorder=1)

    # ── Simulation: filled red histogram (background layer) ───────────────────
    ax.bar(bc[s_fit],  sim_yn[s_fit],  width=bw*0.88,
           color=C_SIM, alpha=0.45, linewidth=0, zorder=2)
    # Excluded tail — same colour but hatched AND much lower alpha to make
    # the contrast with the fitted region immediately obvious
    ax.bar(bc[s_tail], sim_yn[s_tail], width=bw*0.88,
           color=C_SIM, alpha=0.12, hatch="////",
           edgecolor=C_SIM, linewidth=0.5, zorder=2)
    if is_muon and bm.sum() > 0:
        ax.bar(bc[bm], sim_yn[bm], width=bw*0.88,
               color=C_SIM, alpha=0.35, edgecolor=C_SIM,
               linewidth=0.5, zorder=2)

    # Sim error band on fitted region
    if s_fit.sum() > 1:
        xe = np.concatenate([bc[s_fit]-bw/2,[bc[s_fit][-1]+bw/2]])
        ax.fill_between(xe,
            np.concatenate([sim_yn[s_fit]-sim_ye[s_fit],[sim_yn[s_fit][-1]-sim_ye[s_fit][-1]]]),
            np.concatenate([sim_yn[s_fit]+sim_ye[s_fit],[sim_yn[s_fit][-1]+sim_ye[s_fit][-1]]]),
            step="post", color=C_SIM, alpha=0.15, linewidth=0, zorder=2)

    # ── Data: blue error bars (foreground layer) ──────────────────────────────
    # Fitted region — solid, full opacity
    ax.errorbar(bc[d_fit],  exp_yn[d_fit],  yerr=exp_ye[d_fit],
                fmt="o", color=C_DATA, ms=4.5, lw=1.0, capsize=2,
                zorder=6)
    # Excluded tail — same markers but very faint
    ax.errorbar(bc[d_tail], exp_yn[d_tail], yerr=exp_ye[d_tail],
                fmt="o", color=C_DATA, alpha=0.18, ms=3, lw=0.7,
                capsize=1.5, zorder=6)

    # ── Fit curves ────────────────────────────────────────────────────────────
    if r_s:
        # Solid in fitted region
        xs_fit = np.linspace(0, cut_s, 2000)
        ax.plot(xs_fit, _moyal(xs_fit, *r_s["popt"]),
                color=C_SIM, lw=2.2, ls="--", zorder=7)
        # Dashed extrapolation beyond truncation
        xs_ext = np.linspace(cut_s, x_max, 1500)
        ax.plot(xs_ext, _moyal(xs_ext, *r_s["popt"]),
                color=C_SIM, lw=1.8, ls=":", alpha=0.7, zorder=7)
    if r_d:
        ax.plot(xs, _moyal(xs, *r_d["popt"]),
                color=C_DATA, lw=2.2, ls="-", zorder=8)

    # Truncation boundary
    ax.axvline(cut_show, color="dimgray", lw=1.0, ls=":", alpha=0.75, zorder=3)

    # ── Two separate legends stacked in upper-right ───────────────────────────
    # Legend 1: colour/type guide (what each colour/marker means)
    from matplotlib.lines  import Line2D
    from matplotlib.patches import Patch
    leg1_handles = [
        Line2D([0],[0], color=C_DATA, lw=2.2, ls="-",  marker="o",
               ms=4.5, label="Data  — Landau fit"),
        Line2D([0],[0], color=C_SIM,  lw=2.2, ls="--",
               label="Simulation  — Landau fit"),
        Patch(facecolor=C_DATA, alpha=0.55, label="Data  (fitted 70%)"),
        Patch(facecolor=C_DATA, alpha=0.15, label="Data  (excluded tail)"),
        Patch(facecolor=C_SIM,  alpha=0.45, label="Sim.  (fitted 70%)"),
        Patch(facecolor=C_SIM,  alpha=0.10, label="Sim.  (excluded tail)"),
    ]
    leg1 = ax.legend(handles=leg1_handles, loc="upper right",
                     framealpha=0.92, edgecolor="#aaaaaa",
                     handlelength=1.6, labelspacing=0.28,
                     borderpad=0.5, fontsize=8.5, ncol=2)
    ax.add_artist(leg1)   # keep leg1 when adding leg2

    # Legend 2: fit results — positioned just below leg1
    # Use bbox_to_anchor in axes fraction so spacing is consistent
    def fmt_result(label, r):
        if r is None: return f"{label}: fit failed"
        return (f"{label}:  "
                f"MPV={r['mpv']:.3f},  "
                f"$\\xi$={r['xi']:.3f},  "
                f"$\\chi^2/\\nu$={r['chi2_red']:.2f}")

    leg2_handles = [
        Line2D([0],[0], color=C_DATA, lw=0, marker="o", ms=0,
               label=fmt_result("Data", r_d)),
        Line2D([0],[0], color=C_SIM,  lw=0, marker="o", ms=0,
               label=fmt_result("Sim.", r_s)),
    ]
    # Place leg2 just below leg1 using bbox_to_anchor in display coords —
    # we render leg1 first to know its height, then position leg2 below it.
    # A fixed axes-fraction offset of 0.27 gives reliable spacing for this
    # figure size without needing a renderer call.
    leg2 = ax.legend(handles=leg2_handles, loc="upper right",
                     bbox_to_anchor=(1.0, 0.875),   # tight below leg1
                     bbox_transform=ax.transAxes,
                     framealpha=0.92, edgecolor="#aaaaaa",
                     handlelength=0, handletextpad=0,
                     labelspacing=0.30, borderpad=0.5,
                     fontsize=8.5, ncol=1)
    ax.add_artist(leg2)

    # ── Axes ──────────────────────────────────────────────────────────────────
    ax.set_yscale("log")
    pos = np.concatenate([exp_yn[exp_yn>0], sim_yn[sim_yn>0]])
    ax.set_ylim(pos.min()*0.3, 3.5)
    ax.set_xlim(0, x_max)
    ax.yaxis.set_major_locator(ticker.LogLocator(base=10,subs=(1.,),numticks=5))
    ax.yaxis.set_minor_locator(
        ticker.LogLocator(base=10,subs=np.arange(2,10)*0.1,numticks=20))
    ax.yaxis.set_minor_formatter(ticker.NullFormatter())
    ax.set_ylabel("Normalised yield  (peak = 1)")
    ax.set_title(pair["title"], pad=4, fontsize=12)
    plt.setp(ax.get_xticklabels(), visible=False)

    try:
        lo = hep.cms.label(loc=0, data=True,
                           rlabel="HIBEAM  Preliminary", ax=ax)
        lo[0].set_text("ESS")
    except Exception: pass

    # ── Pull panel ────────────────────────────────────────────────────────────
    for r, yn, ye, col, mk in [(r_d, exp_yn, exp_ye, C_DATA, "o"),
                                (r_s, sim_yn, sim_ye, C_SIM,  "s")]:
        if r is None: continue
        yf   = _moyal(bc, *r["popt"])
        den  = np.sqrt(np.maximum(np.abs(yf), 1e-10))
        pls  = np.clip((yn-yf)/den, -8, 8)
        perr = np.clip(ye/den, 0, 8)
        fm   = r["fm"]
        not_bm = ~bm if is_muon else np.ones(len(bc), bool)
        axr.errorbar(bc[fm], pls[fm], yerr=perr[fm],
                     fmt=mk, color=col, ms=4, lw=1., capsize=2, zorder=4)
        tail_m = ~fm & not_bm & (bc<=x_max)
        if tail_m.sum()>0:
            axr.errorbar(bc[tail_m], pls[tail_m], yerr=perr[tail_m],
                         fmt=mk, color=col, alpha=0.20, ms=3,
                         lw=0.7, capsize=1.5, zorder=3)

    axr.axhline(0,  color="black",   lw=0.9)
    axr.axhline(+2, color="#888888", lw=0.9, ls="--", label=r"$\pm2\sigma$")
    axr.axhline(-2, color="#888888", lw=0.9, ls="--")
    axr.axvline(cut_show, color="dimgray", lw=1., ls=":", alpha=0.7)
    axr.set_ylim(-4.5, 4.5)
    axr.yaxis.set_major_locator(ticker.MultipleLocator(2.))
    axr.legend(loc="upper right", fontsize=8, framealpha=0.92,
               edgecolor="#aaaaaa", handlelength=1.3, borderpad=0.4)
    axr.set_ylabel(r"$(N-F)/\!\sqrt{F}$", fontsize=10)
    axr.set_xlabel(r"dE/dx  [units of peak MPV]")

    for ext in [".pdf", ".png"]:
        fig.savefig(f"{OUTDIR}/{pair['output']}{ext}",
                    dpi=300 if ext==".png" else None, bbox_inches="tight")
        print(f"  Saved: {pair['output']}{ext}")
    plt.close(fig)
    return r_d, r_s, exp_pk_raw, sim_pk_raw

# ══════════════════════════════════════════════════════════════════════════════
# MAIN
# ══════════════════════════════════════════════════════════════════════════════
setup_pyroot()

results_log = []
for pair in PAIRS:
    print(f"\n{'='*58}\n  {pair['title']}\n{'='*58}")
    try:
        exp_vals = load_experimental(pair["exp_file"])
    except Exception as e:
        print(f"  SKIP — exp failed: {e}"); continue
    try:
        sim_vals, sim_wts = load_simulation(pair["sim_file"],
                                            pair["sim_low_mev"])
    except Exception as e:
        print(f"  SKIP — sim failed: {e}"); continue
    if len(exp_vals) < 50 or len(sim_vals) < 50:
        print("  SKIP — too few events"); continue
    r_d, r_s, epk, spk = make_comparison_plot(pair, exp_vals, sim_vals, sim_wts)
    results_log.append((pair, r_d, r_s, epk, spk))

# ── Write all fit results to a single text file ──────────────────────────────
with open(f"{OUTDIR}/fit_results_all_runs.txt", "w") as fout:
    fout.write("HIBEAM dE/dx  Data vs Simulation  —  Moyal (Landau) Fit Results\n")
    fout.write(f"Truncation: bottom {TRUNCATION*100:.0f}%\n")
    fout.write(f"Track cuts: chi2/ndf < {CHI2_MAX},  nPoints >= {MIN_POINTS}\n")
    fout.write(f"Sim cuts:   nSteps >= {MIN_STEPS}\n\n")
    for pair, r_d, r_s, epk, spk in results_log:
        fout.write(f"{'='*60}\n{pair['title']}\n{'='*60}\n")
        fout.write(f"Data peak (raw): {epk:.3f} ADC/mm\n")
        fout.write(f"Sim  peak (raw): {spk:.4f} MeV\n")
        if pair["sim_low_mev"] > 0:
            fout.write(f"Sim low-E cut:   {pair['sim_low_mev']:.0f} MeV\n")
        for tag, r in [("Data", r_d), ("Simulation", r_s)]:
            if r is None:
                fout.write(f"  {tag}: fit failed\n")
            else:
                fout.write(f"  {tag}:\n")
                fout.write(f"    MPV        = {r['mpv']:.4f}\n")
                fout.write(f"    xi         = {r['xi']:.4f}\n")
                fout.write(f"    xi_err     = {r['xi_err']:.4f}\n")
                fout.write(f"    chi2/ndf   = {r['chi2_red']:.2f}\n")
        fout.write("\n")
    print(f"\nFit results written to {OUTDIR}/fit_results_all_runs.txt")

print("\nAll done.")