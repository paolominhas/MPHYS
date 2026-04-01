#!/usr/bin/env python3
"""
sim_dedx_plots.py
=================
Publication-quality dE/dx simulation plots showing the effect of the
low-energy cut that removes the sub-threshold bump before the Landau peak.

For each dataset (Krakow / Muon) a single figure shows:
  - Full histogram of all accepted events on a LOG y-axis
  - Orange shading over the bump region (E ≤ LOW_MEV)
  - Dashed grey fit WITHOUT the cut  (biased — shown for comparison only)
  - Solid black fit WITH the cut     (correct Landau fit to MIP peak)
  - Pull residuals for the WITH-CUT fit (bottom panel)

The cut is a single line:
    mask = values > LOW_MEV    (currently 0.15 MeV)
  Applied after the Geant4 step count filter.

Low-energy bump physics (cosmic muons):
  1. Delta-ray electrons  — short secondary e⁻ from muon ionisation
  2. Clipping tracks      — muons crossing only part of the active volume
  3. Michel electrons     — from muon decay, E_max = 52 MeV
  All three deposit less than a full MIP, producing a secondary peak
  below the main Landau.  Ref: Bichsel Rev.Mod.Phys. 60 (1988) 663.

Usage:  python3 sim_dedx_plots.py
"""

import numpy as np
import uproot, awkward as ak
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.ticker as ticker
import mplhep as hep
from scipy import optimize, stats as sc

# ── Config ────────────────────────────────────────────────────────────────────
FILES = {
    "Krakow": {"path": "../simulation_data/KrakowScatter.root",
               "low_mev": 0.0},   # no low-energy cut
    "Muon"  : {"path": "../simulation_data/MuonScatter_fixed.root",
               "low_mev": 15.0},  # remove sub-threshold bump
}
OUTDIR     = "."
TRUNCATION = 0.70
N_BINS     = 35
MIN_STEPS  = 3
LOW_MEV    = 15.0   # cut threshold in MeV — removes low-energy bump

# ── Style ─────────────────────────────────────────────────────────────────────
plt.style.use(hep.style.CMS)
plt.rcParams.update({
    "font.size": 12, "axes.labelsize": 13,
    "xtick.labelsize": 10, "ytick.labelsize": 10,
    "legend.fontsize": 8.5,
    "xtick.direction": "in", "ytick.direction": "in",
    "xtick.top": True, "ytick.right": True,
})

C_HIST  = "#4393c3"   # blue   — main histogram
C_BUMP  = "#f4a582"   # orange — bump region
C_NOCUT = "#888888"   # grey   — fit without cut
C_CUT   = "#000000"   # black  — fit with cut

# ── Load ──────────────────────────────────────────────────────────────────────
def load(path):
    """Load all events with n_steps >= MIN_STEPS and Edep > 0 (no cut yet)."""
    with uproot.open(path) as f:
        tree = f["hibeam"]
        edep = tree["ProtoTPC/Edep"].array(library="ak")
        n    = ak.to_numpy(ak.num(edep)).astype(int)
        e    = ak.to_numpy(ak.sum(edep, axis=1)).astype(float) * 1000.
        try:
            w = ak.to_numpy(ak.fill_none(
                ak.firsts(tree["PrimaryWeight"].array(library="ak")), 1.0))
        except Exception:
            w = np.ones(len(e))

    mask = (n >= MIN_STEPS) & (e > 0)
    v, wt = e[mask], w[mask]
    n_bump = int((v <= LOW_MEV).sum())
    print(f"  {len(v):,} events with hits  "
          f"(range {v.min():.3f}–{v.max():.1f} MeV)")
    print(f"  Bump (E ≤ {LOW_MEV} MeV): {n_bump:,} "
          f"({100*n_bump/len(v):.1f}%)")
    return v, wt

# ── Moyal ─────────────────────────────────────────────────────────────────────
def _moyal(x, loc, scale, norm):
    return norm * sc.moyal.pdf(x, loc=loc, scale=scale)

# ── Fit a subset on shared edges ─────────────────────────────────────────────
def fit_subset(values, weights, edges, peak, label=""):
    bc  = 0.5 * (edges[:-1] + edges[1:])
    bw  = float(edges[1] - edges[0])
    cut = float(np.percentile(values, TRUNCATION * 100))
    fm  = bc <= cut

    counts, _ = np.histogram(values, bins=edges, weights=weights)
    yn = counts / peak
    ye = np.maximum(np.sqrt(counts), 1.0) / peak

    if fm.sum() < 4:
        print(f"  [{label}] too few bins — skip"); return None

    xf, yf, ef = bc[fm], yn[fm], ye[fm]
    loc0   = float(xf[np.argmax(yf)])
    scale0 = max((cut - float(bc[0])) / 8.0, 1e-9)
    norm0  = peak * bw
    try:
        popt, pcov = optimize.curve_fit(
            _moyal, xf, yf, p0=[loc0, scale0, norm0],
            sigma=ef, absolute_sigma=True,
            bounds=([float(bc[0]), 1e-9, 0], [cut, cut, norm0 * 20]),
            maxfev=40_000)
    except Exception as exc:
        print(f"  [{label}] fit failed: {exc}"); return None

    perr     = np.sqrt(np.diag(pcov))
    chi2_red = float(np.sum(((yf - _moyal(xf, *popt)) / ef)**2)) / max(fm.sum()-3, 1)
    print(f"  [{label}] MPV={popt[0]:.4g} MeV  xi={popt[1]:.4g}±{perr[1]:.3g}"
          f"  chi²/ndf={chi2_red:.2f}")
    return {"popt": popt, "perr": perr, "mpv": popt[0], "xi": popt[1],
            "xi_err": perr[1], "chi2_red": chi2_red, "cut": cut,
            "bc": bc, "bw": bw, "yn": yn, "ye": ye, "fm": fm}

# ── Main plot ─────────────────────────────────────────────────────────────────
def make_plot(label, values, weights, outname, low_mev=0.0):
    # Shared axis range
    c0, e0 = np.histogram(values, bins=200, weights=weights)
    mpv_g  = float(0.5 * (e0[:-1] + e0[1:])[np.argmax(c0)])
    x_max  = min(float(np.percentile(values, 99.5)), mpv_g * 7)
    edges  = np.linspace(0, x_max, N_BINS + 1)
    bc     = 0.5 * (edges[:-1] + edges[1:])
    bw     = float(edges[1] - edges[0])

    # Full histogram
    counts_all, _ = np.histogram(values, bins=edges, weights=weights)
    peak     = float(counts_all.max())
    yn_all   = counts_all / peak
    ye_all   = np.maximum(np.sqrt(counts_all), 1.0) / peak

    # Fit 1: no cut (always shown, even for muons, to show the bias)
    print(f"\n  Fit 1 — no low-energy cut:")
    r_nc = fit_subset(values, weights, edges, peak, "no cut") if low_mev > 0 else None

    # Fit 2: with cut (same as fit 1 for Krakow where low_mev=0)
    if low_mev > 0:
        v2, w2 = values[values > low_mev], weights[values > low_mev]
        print(f"\n  Fit 2 — E > {low_mev:.0f} MeV cut ({len(v2):,} events):")
        r_c = fit_subset(v2, w2, edges, peak, f"E>{low_mev:.0f}")
    else:
        print(f"\n  Single fit (no cut applied for {label}):")
        r_c = fit_subset(values, weights, edges, peak, label)
        r_nc = None  # only one fit for Krakow

    xs = np.linspace(0, x_max, 3000)

    # ── Figure ────────────────────────────────────────────────────────────────
    fig = plt.figure(figsize=(8.5, 7.2))
    gs  = gridspec.GridSpec(2, 1, height_ratios=[3, 1], hspace=0.0,
                            left=0.13, right=0.97, top=0.91, bottom=0.11)
    ax  = fig.add_subplot(gs[0])
    axr = fig.add_subplot(gs[1], sharex=ax)

    if low_mev > 0:
        ax.axvspan(0, low_mev, color=C_BUMP, alpha=0.30, zorder=1,
                   label=f"Sub-threshold bump  ($E \\leq {low_mev:.0f}$ MeV, excluded)")

    # Determine the 70% truncation cut for display (use WITH-cut fit if available,
    # else compute directly)
    cut_display = r_c["cut"] if r_c is not None else float(
        np.percentile(values[values > low_mev] if low_mev > 0 else values,
                      TRUNCATION * 100))

    bm      = bc <= low_mev  # bump region (empty for Krakow since low_mev=0)
    tail_m  = (~bm) & (bc > cut_display)
    fit_m   = (~bm) & (bc <= cut_display)

    ax.bar(bc[bm],     yn_all[bm],     width=bw*0.88, color=C_BUMP,
           alpha=0.80, linewidth=0, zorder=2)
    ax.bar(bc[fit_m],  yn_all[fit_m],  width=bw*0.88, color=C_HIST,
           alpha=0.55, linewidth=0, zorder=2,
           label=f"Simulation — fitted (bottom {int(TRUNCATION*100)}%)")
    ax.bar(bc[tail_m], yn_all[tail_m], width=bw*0.88, color="#d73027",
           alpha=0.28, hatch="////", edgecolor="#d73027", linewidth=0.5,
           zorder=2, label=r"Excluded tail  ($\delta$-ray enhanced)")

    # Statistical error band — fitted region only
    if fit_m.sum() > 1:
        xe = np.concatenate([bc[fit_m] - bw/2, [bc[fit_m][-1] + bw/2]])
        ax.fill_between(xe,
            np.concatenate([yn_all[fit_m] - ye_all[fit_m],
                            [yn_all[fit_m][-1] - ye_all[fit_m][-1]]]),
            np.concatenate([yn_all[fit_m] + ye_all[fit_m],
                            [yn_all[fit_m][-1] + ye_all[fit_m][-1]]]),
            step="post", color=C_HIST, alpha=0.18, linewidth=0, zorder=2)

    # Fit curves
    if r_nc is not None:
        ax.plot(xs, _moyal(xs, *r_nc["popt"]),
                color=C_NOCUT, lw=1.6, ls="--", zorder=4,
                label=(f"Fit — no cut  "
                       f"(MPV = {r_nc['mpv']:.3g} MeV,  "
                       f"$\\chi^2/\\nu$ = {r_nc['chi2_red']:.0f})"))

    if r_c is not None:
        cut_label = (f"Fit — $E > {low_mev:.0f}$ MeV  " if low_mev > 0
                     else "Landau (Moyal) fit  ")
        ax.plot(xs, _moyal(xs, *r_c["popt"]),
                color=C_CUT, lw=2.0, ls="-", zorder=5,
                label=(cut_label +
                       f"(MPV = {r_c['mpv']:.3g} MeV,  "
                       f"$\\chi^2/\\nu$ = {r_c['chi2_red']:.2f})"))
        ax.axvline(r_c["cut"], color="dimgray", lw=1.0, ls=":",
                   alpha=0.75, zorder=3,
                   label=f"70% truncation at {r_c['cut']:.1f} MeV")

    # Log scale and axes
    ax.set_yscale("log")
    ax.set_xlim(0, x_max)
    pos_min = yn_all[yn_all > 0].min()
    ax.set_ylim(bottom=pos_min * 0.3, top=3.0)
    ax.yaxis.set_major_locator(
        ticker.LogLocator(base=10, subs=(1.0,), numticks=6))
    ax.yaxis.set_minor_locator(
        ticker.LogLocator(base=10, subs=np.arange(2, 10)*0.1, numticks=20))
    ax.yaxis.set_minor_formatter(ticker.NullFormatter())

    # Annotation box — lower right, away from legend (upper right)
    if r_c is not None:
        cut_str = f"With cut ($E > {low_mev:.0f}$ MeV):" if low_mev > 0 else "Fit result:"
        ann = "\n".join([
            cut_str,
            f"  MPV $= {r_c['mpv']:.3g}$ MeV",
            f"  $\\xi = {r_c['xi']:.3g} \\pm {r_c['xi_err']:.2g}$ MeV",
            f"  $\\chi^2/\\nu = {r_c['chi2_red']:.2f}$",
        ])
        ax.text(0.97, 0.04, ann, transform=ax.transAxes,
                va="bottom", ha="right", fontsize=9,
                bbox=dict(boxstyle="round,pad=0.4", facecolor="white",
                          edgecolor="#aaaaaa", alpha=0.95, zorder=6))

    ax.set_ylabel("Normalised yield  (peak = 1,  log scale)")
    # Legend: upper right, small font to avoid overlap
    ax.legend(loc="upper right", framealpha=0.92, edgecolor="#aaaaaa",
              handlelength=1.8, labelspacing=0.30, borderpad=0.5)
    plt.setp(ax.get_xticklabels(), visible=False)

    try:
        lo = hep.cms.label(loc=0, data=False,
                           rlabel="HIBEAM  Simulation", ax=ax)
        lo[0].set_text("ESS")
    except Exception:
        pass

    # ── Pull panel ────────────────────────────────────────────────────────────
    if r_c is not None:
        fm   = r_c["fm"]
        yf   = _moyal(bc, *r_c["popt"])
        den  = np.sqrt(np.maximum(np.abs(yf), 1e-10))
        pls  = np.clip((yn_all - yf) / den, -8, 8)
        perr = np.clip(ye_all / den, 0, 8)

        axr.errorbar(bc[fit_m],  pls[fit_m],  yerr=perr[fit_m],
                     fmt="o", color=C_HIST,   ms=4, lw=1.0, capsize=2,
                     label="Fitted region")
        if tail_m.sum() > 0:
            axr.errorbar(bc[tail_m], pls[tail_m], yerr=perr[tail_m],
                         fmt="^", color="#d73027", alpha=0.70, ms=4,
                         lw=1.0, capsize=2, label="Excluded tail")
        if bm.sum() > 0:
            axr.errorbar(bc[bm], pls[bm], yerr=perr[bm],
                         fmt="s", color=C_BUMP, ms=4, lw=1.0, capsize=2,
                         label=f"Bump ($E\\leq{low_mev:.0f}$ MeV)")

    axr.axhline(0,  color="black",   lw=0.9)
    axr.axhline(+2, color="#888888", lw=0.9, ls="--")
    axr.axhline(-2, color="#888888", lw=0.9, ls="--",
                label=r"$\pm2\sigma$")
    if r_c is not None:
        axr.axvline(r_c["cut"], color="dimgray", lw=1.0, ls=":", alpha=0.7)
    if low_mev > 0:
        axr.axvspan(0, low_mev, color=C_BUMP, alpha=0.15)

    axr.set_ylim(-5.0, 5.0)
    axr.yaxis.set_major_locator(ticker.MultipleLocator(2.0))
    axr.legend(loc="upper left", fontsize=7.5, framealpha=0.92,
               edgecolor="#aaaaaa", handlelength=1.4,
               labelspacing=0.20, ncol=2, borderpad=0.4)
    axr.set_ylabel(r"$(N-F)/\!\sqrt{F}$", fontsize=10)
    axr.set_xlabel(r"$\sum E_{\mathrm{dep}}\ [\mathrm{MeV}]$")

    for ext in [".pdf", ".png"]:
        fig.savefig(f"{OUTDIR}/{outname}{ext}",
                    dpi=300 if ext == ".png" else None,
                    bbox_inches="tight")
        print(f"  Saved: {outname}{ext}")
    plt.close(fig)

# ── Run ───────────────────────────────────────────────────────────────────────
for label, cfg in FILES.items():
    print(f"\n{'='*55}\n  {label}\n{'='*55}")
    vals, wts = load(cfg["path"])
    make_plot(label, vals, wts, f"sim_{label.lower()}_dedx",
              low_mev=cfg["low_mev"])
print("\nDone.")