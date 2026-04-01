#!/usr/bin/env python3
"""
sim_protoTPC_dedx.py
====================
Produces publication-quality dE/dx plots for the HIBEAM Proto-TPC
simulation data, using ProtoTPC/Edep summed per event.

Runs on both Krakow (proton scattering) and Muon (MCPL) simulation
files.  Processing, binning, fitting and plotting are identical to
tpc_dedx_analysis.py.

Usage
-----
    python sim_protoTPC_dedx.py
    (run from combined_analysis/ directory)

References
----------
[1] Landau, L.D. (1944). J. Phys. USSR 8, 201.
[2] Moyal, J.E. (1955). Phil. Mag. 46, 263.
[3] Bichsel, H. (1988). Rev. Mod. Phys. 60, 663.
[4] ALICE Collaboration (2010). JINST 5, P09002.
[5] STAR Collaboration (2003). Nucl. Instrum. Meth. A 499, 659.
[6] PDG (2022). Prog. Theor. Exp. Phys. 2022, 083C01.
"""

import os
import sys
from pathlib import Path

import awkward as ak
import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import mplhep as hep
import numpy as np
import uproot
from scipy import optimize, stats


# ══════════════════════════════════════════════════════════════════════════════
# CONFIG
# ══════════════════════════════════════════════════════════════════════════════
SIMS = [
    {
        "file"      : "../simulation_data/KrakowScatter.root",
        "title"     : r"Proto-TPC $\mathrm{d}E/\mathrm{d}x$  —  Krakow $^2$H(p,p)",
        "output"    : "sim_dedx_krakow",
        "min_mev"   : 0.0,      # no low-energy cut
        "x_max"     : 40.0,     # MeV/mm
    },
    {
        "file"      : "../simulation_data/MuonScatter_fixed.root",
        "title"     : r"Proto-TPC $\mathrm{d}E/\mathrm{d}x$  —  Muon lab (MCPL)",
        "output"    : "sim_dedx_muon",
        "min_mev"   : 0.0,      # keep everything; fit handles truncation
        "x_max"     : 300.0,    # MeV/mm
    },
]

TRUNCATION = 0.70     # fit bottom 70 %  [ALICE JINST 5 (2010) P09002]
N_BINS     = 40       # coarse bins for publication
MIN_STEPS  = 5        # minimum G4 steps per event
OUTDIR     = "."

# ══════════════════════════════════════════════════════════════════════════════
# STYLE
# ══════════════════════════════════════════════════════════════════════════════
plt.style.use(hep.style.CMS)
plt.rcParams.update({
    "font.size"       : 12,
    "axes.labelsize"  : 14,
    "xtick.labelsize" : 11,
    "ytick.labelsize" : 11,
    "legend.fontsize" : 9,
    "xtick.direction" : "in",
    "ytick.direction" : "in",
    "xtick.top"       : True,
    "ytick.right"     : True,
})

C_FIT  = "#2166ac"
C_TAIL = "#d73027"


# ══════════════════════════════════════════════════════════════════════════════
# DATA LOADING  (uproot only — no PyROOT needed)
# ══════════════════════════════════════════════════════════════════════════════
def load_protoTPC(filepath, min_steps=5, min_mev=0.0):
    """
    Load ProtoTPC/Edep, sum per event, convert GeV → MeV.

    Returns
    -------
    values : np.ndarray — per-event sum(Edep) [MeV]
    """
    print(f"  Loading: {Path(filepath).name}")
    with uproot.open(filepath) as f:
        tree = f["hibeam"]
        keys = set(tree.keys(recursive=True))

        # Find branches
        def first_key(candidates):
            for k in candidates:
                if k in keys:
                    return k
            return None

        edep_key = first_key(["ProtoTPC/Edep", "ProtoTPC.Edep"])
        px_key   = first_key(["ProtoTPC/Pos_X", "ProtoTPC.Pos_X"])
        py_key   = first_key(["ProtoTPC/Pos_Y", "ProtoTPC.Pos_Y"])
        pz_key   = first_key(["ProtoTPC/Pos_Z", "ProtoTPC.Pos_Z"])

        if edep_key is None:
            raise KeyError(f"No ProtoTPC/Edep branch. Keys: {list(keys)[:20]}")

        edep = tree[edep_key].array(library="ak")

        # Position branches for path length
        has_pos = all(k is not None for k in [px_key, py_key, pz_key])
        if has_pos:
            pos_x = tree[px_key].array(library="ak")
            pos_y = tree[py_key].array(library="ak")
            pos_z = tree[pz_key].array(library="ak")
            print("    Position branches found → computing dE/dx [MeV/mm]")
        else:
            print("    WARNING: no Pos_X/Y/Z branches → falling back to sum(Edep)")

    n_steps = ak.to_numpy(ak.num(edep)).astype(int)
    e_sum   = ak.to_numpy(ak.sum(edep, axis=1)).astype(float) * 1000.0  # GeV → MeV

    if has_pos:
        # 3D path length from extent of hit positions [mm]
        # Geant4 positions are in mm by default
        dx = ak.to_numpy(ak.max(pos_x, axis=1) - ak.min(pos_x, axis=1)).astype(float)
        dy = ak.to_numpy(ak.max(pos_y, axis=1) - ak.min(pos_y, axis=1)).astype(float)
        dz = ak.to_numpy(ak.max(pos_z, axis=1) - ak.min(pos_z, axis=1)).astype(float)
        path_len = np.sqrt(dx**2 + dy**2 + dz**2)  # mm

        mask = (n_steps >= min_steps) & (e_sum > max(min_mev, 0.0)) & (path_len > 0.1)
        values = e_sum[mask] / path_len[mask]  # MeV/mm
        unit = "MeV/mm"
    else:
        mask = (n_steps >= min_steps) & (e_sum > max(min_mev, 0.0))
        values = e_sum[mask]
        unit = "MeV"

    print(f"    {len(e_sum):,} events total → {mask.sum():,} accepted "
          f"(min_steps>={min_steps}"
          + (f", E>{min_mev:.0f} MeV" if min_mev > 0 else "") + ")")
    print(f"    Range: {values.min():.4f} – {values.max():.2f} {unit}")

    return values


# ══════════════════════════════════════════════════════════════════════════════
# LANDAU FIT  (identical to tpc_dedx_analysis.py)
# ══════════════════════════════════════════════════════════════════════════════
def _moyal_scaled(x, loc, scale, norm):
    return norm * stats.moyal.pdf(x, loc=loc, scale=scale)


def fit_landau(values, truncation=TRUNCATION, n_bins=N_BINS):
    if len(values) < 50:
        raise ValueError(f"Only {len(values)} events — need >= 50.")

    cut_value = np.percentile(values, truncation * 100.0)

    hist_min = float(values.min())
    hist_max = float(np.percentile(values, 99.0))
    fit_range = cut_value - hist_min
    if fit_range <= 0:
        raise ValueError("Truncation cut equals data minimum.")

    bin_width_target = fit_range / n_bins
    n_total = max(int(np.ceil((hist_max - hist_min) / bin_width_target)),
                  n_bins + 1)
    counts, bin_edges = np.histogram(values, bins=n_total,
                                     range=(hist_min, hist_max))
    bin_centers = 0.5 * (bin_edges[:-1] + bin_edges[1:])
    bin_width   = float(bin_edges[1] - bin_edges[0])

    fit_mask = bin_centers <= cut_value
    valid_for_fit = fit_mask & (counts >= 5)
    if valid_for_fit.sum() < 4:
        raise RuntimeError(f"Only {valid_for_fit.sum()} usable bins.")

    x_fit = bin_centers[valid_for_fit]
    y_fit = counts[valid_for_fit].astype(float)

    peak_in_fit = np.argmax(counts[fit_mask])
    loc0   = float(bin_centers[fit_mask][peak_in_fit])
    scale0 = float((cut_value - bin_edges[0]) / 8.0)
    norm0  = float(counts.sum()) * bin_width

    try:
        popt, pcov = optimize.curve_fit(
            _moyal_scaled, x_fit, y_fit,
            p0=[loc0, scale0, norm0],
            sigma=np.sqrt(np.maximum(y_fit, 1.0)),
            absolute_sigma=True, maxfev=40_000,
            bounds=([float(x_fit[0]), 1e-9, 0],
                    [cut_value,       cut_value, norm0 * 20]))
    except RuntimeError as exc:
        raise RuntimeError(f"Fit did not converge: {exc}") from exc

    loc, scale, norm = popt
    loc_err, scale_err, norm_err = np.sqrt(np.diag(pcov))

    fit_curve = _moyal_scaled(bin_centers, *popt)

    exp_fit  = fit_curve[valid_for_fit]
    # chi2 using same sigma as passed to curve_fit: sqrt(max(observed, 1))
    sigma_fit = np.sqrt(np.maximum(y_fit, 1.0))
    chi2     = float(np.sum(((y_fit - exp_fit) / sigma_fit)**2))
    ndf      = int(valid_for_fit.sum()) - 3
    chi2_red = chi2 / max(ndf, 1)
    p_value  = float(stats.chi2.sf(chi2, max(ndf, 1)))

    print(f"    MPV = {loc:.4g} MeV/mm,  xi = {scale:.4g} MeV/mm,  "
          f"chi2/ndf = {chi2_red:.2f}  (p = {p_value:.3f})")

    return {
        "loc": loc, "loc_err": loc_err,
        "scale": scale, "scale_err": scale_err,
        "norm": norm, "norm_err": norm_err,
        "mpv": loc, "popt": popt, "pcov": pcov,
        "chi2": chi2, "ndf": ndf,
        "chi2_red": chi2_red, "p_value": p_value,
        "cut_value": cut_value,
        "bin_edges": bin_edges, "bin_centers": bin_centers,
        "counts": counts, "fit_curve": fit_curve,
        "fit_mask": fit_mask, "bin_width": bin_width,
        "truncation": truncation,
    }


# ══════════════════════════════════════════════════════════════════════════════
# PLOTTING  (identical style to tpc_dedx_analysis.py)
# ══════════════════════════════════════════════════════════════════════════════
def _draw_panels(ax_top, ax_res, fit_result, title, x_label,
                 x_max_display, log_y):
    bc  = fit_result["bin_centers"]
    cnt = fit_result["counts"].astype(float)
    fc  = fit_result["fit_curve"]
    fm  = fit_result["fit_mask"]
    bw  = fit_result["bin_width"]

    cnt_err  = np.sqrt(np.maximum(cnt, 1.0))
    denom    = np.sqrt(np.maximum(fc, 1.0))
    pulls    = (cnt - fc) / denom

    trunc_pct = int(fit_result["truncation"] * 100)

    # ── Histograms ────────────────────────────────────────────────────────
    ax_top.bar(bc[fm], cnt[fm], width=bw * 0.92,
               color=C_FIT, alpha=0.55,
               label=f"Data (bottom {trunc_pct}%)")
    ax_top.bar(bc[~fm], cnt[~fm], width=bw * 0.92,
               color=C_TAIL, alpha=0.30, hatch="////",
               edgecolor=C_TAIL, linewidth=0.5,
               label=f"Excluded tail (top {100 - trunc_pct}%)")

    # ── Error bands ───────────────────────────────────────────────────────
    for mask, color, alpha in [(fm, C_FIT, 0.30), (~fm, C_TAIL, 0.15)]:
        if mask.sum() == 0:
            continue
        x_e = np.empty(2 * mask.sum())
        x_e[0::2] = bc[mask] - bw / 2
        x_e[1::2] = bc[mask] + bw / 2
        y_lo = np.repeat(np.maximum(cnt[mask] - cnt_err[mask], 0.0), 2)
        y_hi = np.repeat(cnt[mask] + cnt_err[mask], 2)
        ax_top.fill_between(x_e, y_lo, y_hi, color=color,
                            alpha=alpha, linewidth=0)

    # ── Fit curve ─────────────────────────────────────────────────────────
    x_dense = np.linspace(float(bc[0]), float(x_max_display), 2000)
    y_dense = _moyal_scaled(x_dense, *fit_result["popt"])
    ax_top.plot(x_dense, y_dense, color="black", lw=2.0, zorder=5,
                label="Landau (Moyal) fit")

    # ── Truncation boundary ───────────────────────────────────────────────
    ax_top.axvline(fit_result["cut_value"], color="dimgray",
                   lw=1.2, ls="--", alpha=0.8,
                   label=f"{trunc_pct}th percentile")

    # ── Truncated mean (fit region only) ──────────────────────────────────
    cnt_fit_sum = float(np.sum(cnt[fm]))
    if cnt_fit_sum > 0:
        mean_fit = float(np.sum(bc[fm] * cnt[fm]) / cnt_fit_sum)
        ax_top.axvline(mean_fit, color="#4daf4a", lw=1.4, ls="-.",
                       label=rf"$\langle \Delta E \rangle_{{\mathrm{{fit}}}}$"
                             rf" = {mean_fit:.4g} MeV/mm")

    # ── MPV ───────────────────────────────────────────────────────────────
    ax_top.axvline(fit_result["mpv"], color="#ff7f00", lw=1.2, ls=":",
                   label=rf"MPV = {fit_result['mpv']:.4g} MeV/mm")

    # ── chi2 annotation ───────────────────────────────────────────────────
    ax_top.text(0.97, 0.50,
                rf"$\chi^2/\nu = {fit_result['chi2_red']:.2f}$",
                transform=ax_top.transAxes, va="top", ha="right",
                fontsize=11,
                bbox=dict(boxstyle="round,pad=0.4", facecolor="white",
                          edgecolor="#aaaaaa", alpha=0.92))

    # ── Axes ──────────────────────────────────────────────────────────────
    ax_top.set_ylabel("Events / bin" + ("  (log scale)" if log_y else ""),
                      fontsize=13)
    ax_top.set_title(title, fontsize=12, pad=4)
    ax_top.legend(loc="upper right", fontsize=8.5,
                  framealpha=0.92, edgecolor="#aaaaaa",
                  handlelength=1.6, labelspacing=0.30)
    ax_top.set_xlim(float(bc[0]) - bw, x_max_display)
    if log_y:
        ax_top.set_yscale("log")
        ax_top.set_ylim(bottom=0.5)
    plt.setp(ax_top.get_xticklabels(), visible=False)

    try:
        lo = hep.cms.label(loc=0, data=False, label="Simulation",
                           rlabel="HIBEAM", ax=ax_top)
        lo[0].set_text("ESS")
    except Exception:
        pass

    # ── Pull panel ────────────────────────────────────────────────────────
    ax_res.bar(bc[fm], pulls[fm], width=bw * 0.92,
               color=C_FIT, alpha=0.55)
    ax_res.bar(bc[~fm], pulls[~fm], width=bw * 0.92,
               color=C_TAIL, alpha=0.30, hatch="////", edgecolor=C_TAIL)

    ax_res.axhline(0,  color="black",   lw=0.9)
    ax_res.axhline(+2, color="#888888", lw=0.9, ls="--",
                   label=r"$\pm2\sigma$")
    ax_res.axhline(-2, color="#888888", lw=0.9, ls="--")
    ax_res.axvline(fit_result["cut_value"], color="dimgray",
                   lw=1.2, ls="--", alpha=0.8)
    ax_res.set_ylabel(r"$(N-F)/\!\sqrt{F}$", fontsize=11)
    ax_res.set_xlabel(x_label, fontsize=13)
    ax_res.set_ylim(-5.5, 5.5)
    ax_res.yaxis.set_major_locator(ticker.MultipleLocator(2.0))
    ax_res.set_xlim(float(bc[0]) - bw, x_max_display)
    ax_res.legend(fontsize=8, loc="upper right",
                  framealpha=0.92, edgecolor="#aaaaaa")


def make_plot(fit_result, title, output_base, x_max=300.0):
    """Produce linear + log PDF figures."""
    x_label = r"$\mathrm{d}E/\mathrm{d}x\ [\mathrm{MeV/mm}]$"

    for log_y, suffix in [(False, "linear"), (True, "log")]:
        fig = plt.figure(figsize=(8.5, 6.5))
        gs  = gridspec.GridSpec(
            2, 1, height_ratios=[3, 1], hspace=0.0,
            left=0.12, right=0.97, top=0.91, bottom=0.11)
        ax_top = fig.add_subplot(gs[0])
        ax_res = fig.add_subplot(gs[1], sharex=ax_top)
        _draw_panels(ax_top, ax_res, fit_result, title,
                     x_label, x_max, log_y)
        path = f"{OUTDIR}/{output_base}_{suffix}.pdf"
        fig.savefig(path, bbox_inches="tight")
        print(f"    Saved: {path}")
        plt.close(fig)


# ══════════════════════════════════════════════════════════════════════════════
# MAIN
# ══════════════════════════════════════════════════════════════════════════════
for cfg in SIMS:
    print(f"\n{'='*55}")
    print(f"  {cfg['title']}")
    print(f"{'='*55}")
    try:
        values = load_protoTPC(cfg["file"], min_steps=MIN_STEPS,
                               min_mev=cfg["min_mev"])
    except Exception as e:
        print(f"  SKIP — {e}")
        continue
    if len(values) < 50:
        print("  SKIP — too few events")
        continue

    print("  Fitting...")
    result = fit_landau(values)
    print("  Plotting...")
    make_plot(result, cfg["title"], cfg["output"], x_max=cfg.get("x_max", 300.0))

print("\nDone.")