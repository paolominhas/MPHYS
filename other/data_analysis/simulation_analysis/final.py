#!/usr/bin/env python3
"""
tpc_dedx_analysis.py
====================
Extracts TPC energy-loss data from a processed HIBEAM ROOT ntuple,
fits a Landau (Moyal) distribution to the bottom 70 % of values using
the truncated-mean method, and produces a publication-quality plot with
residuals and goodness-of-fit statistics.

The load_tpc_data() function is intentionally standalone so it can be
imported cleanly by other analysis scripts.

Usage
-----
    python tpc_dedx_analysis.py /path/to/HIBEAMScatter.root [figure.pdf]

References
----------
[1] Landau, L.D. (1944). J. Phys. USSR 8, 201.
[2] Moyal, J.E. (1955). Phil. Mag. 46, 263.
[3] Bichsel, H. (1988). Rev. Mod. Phys. 60, 663.
[4] ALICE Collaboration (2010). JINST 5, P09002.
[5] STAR Collaboration (2003). Nucl. Instrum. Meth. A 499, 659.
[6] Particle Data Group (2022). Prog. Theor. Exp. Phys. 2022, 083C01.
"""

import sys
import warnings
from dataclasses import dataclass
from typing import Dict, Optional

import awkward as ak
import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import mplhep as hep
import numpy as np
import uproot
from scipy import optimize, stats


# ── Publication style ─────────────────────────────────────────────────────────
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


# ══════════════════════════════════════════════════════════════════════════════
# SECTION 1 — TPC GEOMETRY
# ══════════════════════════════════════════════════════════════════════════════

@dataclass
class TPCGeometry:
    pad_width:  Optional[float] = None
    pad_height: Optional[float] = None
    drift_v:    Optional[float] = None
    time_bin:   Optional[float] = None

    @property
    def is_complete(self) -> bool:
        return all(v is not None for v in
                   [self.pad_width, self.pad_height, self.drift_v, self.time_bin])


# ══════════════════════════════════════════════════════════════════════════════
# SECTION 2 — DATA LOADING
# ══════════════════════════════════════════════════════════════════════════════

def load_tpc_data(filepath: str, tree_name: str = "hibeam") -> Dict:
    with uproot.open(filepath) as root_file:
        available = list(root_file.keys())
        if tree_name not in root_file:
            raise KeyError(
                f"TTree '{tree_name}' not found in {filepath}.\n"
                f"Available keys: {available}")
        tree = root_file[tree_name]
        all_keys = tree.keys(recursive=True)

        tpc_keys = [
            k for k in all_keys
            if ("TPC" in k)
            and not k.startswith("ProtoTPC")
            and not k.startswith("target")
        ]
        if not tpc_keys:
            raise KeyError(f"No TPC sub-branches found.\nAll keys: {all_keys}")

        def find_branch(suffix):
            for k in tpc_keys:
                if k.split("/")[-1] == suffix or k.split(".")[-1] == suffix:
                    return k
            return None

        edep_key = find_branch("Edep")
        if edep_key is None:
            raise KeyError(f"TPC Edep not found.\nTPC keys: {tpc_keys}")

        def safe_array(suffix):
            key = find_branch(suffix)
            if key is None:
                return None
            return tree[key].array(library="ak")

        edep = tree[edep_key].array(library="ak")
        return {
            "edep"      : edep,
            "pad_row"   : safe_array("padRow"),
            "pad_col"   : safe_array("padColumn"),
            "timestamp" : safe_array("timestamp"),
            "n_el"      : safe_array("nEl"),
            "n_events"  : len(edep),
        }


# ══════════════════════════════════════════════════════════════════════════════
# SECTION 3 — dE/dx COMPUTATION
# ══════════════════════════════════════════════════════════════════════════════

def _track_length_mm(pad_row, pad_col, timestamp, geom):
    if len(pad_row) < 2:
        return float("nan")
    dx = (float(ak.max(pad_col))   - float(ak.min(pad_col)))   * geom.pad_width
    dy = (float(ak.max(pad_row))   - float(ak.min(pad_row)))   * geom.pad_height
    dz = (float(ak.max(timestamp)) - float(ak.min(timestamp))) * geom.drift_v * geom.time_bin
    return float(np.sqrt(dx**2 + dy**2 + dz**2))


def compute_dedx(data, geom=None, min_steps=5):
    edep    = data["edep"]
    pad_row = data.get("pad_row")
    pad_col = data.get("pad_col")
    ts      = data.get("timestamp")

    use_geom = (geom is not None and geom.is_complete
                and pad_row is not None and pad_col is not None and ts is not None)

    result = []
    for i in range(data["n_events"]):
        ev_edep = edep[i]
        if len(ev_edep) < min_steps:
            continue
        e_total = float(ak.sum(ev_edep)) * 1000.0  # GeV → MeV
        if not use_geom:
            result.append(e_total)
        else:
            length = _track_length_mm(pad_row[i], pad_col[i], ts[i], geom)
            if np.isnan(length) or length <= 0.0:
                continue
            result.append(e_total / length)

    return np.asarray(result, dtype=np.float64)


# ══════════════════════════════════════════════════════════════════════════════
# SECTION 4 — LANDAU FIT
# ══════════════════════════════════════════════════════════════════════════════

def _moyal_scaled(x, loc, scale, norm):
    return norm * stats.moyal.pdf(x, loc=loc, scale=scale)


def fit_landau(values, truncation=0.70, n_bins=40):
    """
    Fit Moyal to the lower `truncation` fraction.
    n_bins=40 gives coarser binning suitable for publication figures.
    """
    if len(values) < 50:
        raise ValueError(f"Only {len(values)} events — need >= 50.")

    cut_value = np.percentile(values, truncation * 100.0)

    hist_min = float(values.min())
    hist_max = float(np.percentile(values, 99.0))
    fit_range = cut_value - hist_min
    if fit_range <= 0:
        raise ValueError("Truncation cut equals data minimum.")
    bin_width_target = fit_range / n_bins
    n_total = max(int(np.ceil((hist_max - hist_min) / bin_width_target)), n_bins + 1)
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
            absolute_sigma=True, maxfev=20_000)
    except RuntimeError as exc:
        raise RuntimeError(f"Fit did not converge: {exc}") from exc

    loc, scale, norm = popt
    loc_err, scale_err, norm_err = np.sqrt(np.diag(pcov))

    fit_curve = _moyal_scaled(bin_centers, loc, scale, norm)

    exp_fit = fit_curve[valid_for_fit]
    chi2    = float(np.sum((y_fit - exp_fit)**2 / np.maximum(exp_fit, 1e-6)))
    ndf     = int(valid_for_fit.sum()) - 3
    chi2_red = chi2 / max(ndf, 1)
    p_value  = float(stats.chi2.sf(chi2, max(ndf, 1)))

    tail_mask = (~fit_mask) & (counts >= 1)
    if tail_mask.sum() > 0:
        obs_tail = counts[tail_mask].astype(float)
        exp_tail = np.maximum(fit_curve[tail_mask], 1e-6)
        chi2_tail = float(np.sum((obs_tail - exp_tail)**2 / exp_tail))
        ndf_tail  = int(tail_mask.sum())
    else:
        chi2_tail, ndf_tail = float("nan"), 0

    return {
        "loc": loc, "loc_err": loc_err,
        "scale": scale, "scale_err": scale_err,
        "norm": norm, "norm_err": norm_err,
        "mpv": loc, "popt": popt, "pcov": pcov,
        "chi2": chi2, "ndf": ndf,
        "chi2_red": chi2_red, "p_value": p_value,
        "chi2_tail": chi2_tail, "ndf_tail": ndf_tail,
        "cut_value": cut_value,
        "bin_edges": bin_edges, "bin_centers": bin_centers,
        "counts": counts, "fit_curve": fit_curve,
        "fit_mask": fit_mask, "bin_width": bin_width,
        "truncation": truncation,
    }


# ══════════════════════════════════════════════════════════════════════════════
# SECTION 5 — PLOTTING
# ══════════════════════════════════════════════════════════════════════════════

C_FIT  = "#2166ac"   # blue  — fitted region
C_TAIL = "#d73027"   # red   — excluded tail


def _draw_panels(ax_top, ax_res, fit_result, x_label, x_max_display, log_y):
    bc  = fit_result["bin_centers"]
    cnt = fit_result["counts"].astype(float)
    fc  = fit_result["fit_curve"]
    fm  = fit_result["fit_mask"]
    bw  = fit_result["bin_width"]

    cnt_err  = np.sqrt(np.maximum(cnt, 1.0))
    denom    = np.sqrt(np.maximum(fc, 1.0))
    pulls    = (cnt - fc) / denom
    pull_err = cnt_err / denom

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
                             rf" = {mean_fit:.4g} MeV")

    # ── MPV ───────────────────────────────────────────────────────────────
    ax_top.axvline(fit_result["mpv"], color="#ff7f00", lw=1.2, ls=":",
                   label=rf"MPV = {fit_result['mpv']:.4g} MeV")

    # ── chi2 annotation (single value, clean) ─────────────────────────────
    ax_top.text(0.97, 0.50,
                rf"$\chi^2/\nu = {fit_result['chi2_red']:.2f}$",
                transform=ax_top.transAxes, va="top", ha="right",
                fontsize=11,
                bbox=dict(boxstyle="round,pad=0.4", facecolor="white",
                          edgecolor="#aaaaaa", alpha=0.92))

    # ── Axes ──────────────────────────────────────────────────────────────
    ax_top.set_ylabel("Events / bin" + ("  (log scale)" if log_y else ""),
                      fontsize=13)
    ax_top.set_title("TPC Energy-Loss Distribution  —  HIBEAM Simulation",
                     fontsize=12, pad=4)
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


def plot_dedx(fit_result, use_dedx=False, output_path=None):
    x_max_display = 250.0

    x_label = (
        r"$\mathrm{d}E/\mathrm{d}x\ [\mathrm{MeV}\,\mathrm{mm}^{-1}]$"
        if use_dedx
        else r"$\sum E_\mathrm{dep}\ [\mathrm{MeV}]$"
    )

    figures = {}
    for log_y, suffix in [(False, "linear"), (True, "log")]:
        fig = plt.figure(figsize=(8.5, 6.5))
        gs  = gridspec.GridSpec(
            2, 1, height_ratios=[3, 1], hspace=0.0,
            left=0.12, right=0.97, top=0.91, bottom=0.11)
        ax_top = fig.add_subplot(gs[0])
        ax_res = fig.add_subplot(gs[1], sharex=ax_top)
        _draw_panels(ax_top, ax_res, fit_result,
                     x_label, x_max_display, log_y)
        figures[suffix] = fig

    if output_path:
        import os
        stem, ext = os.path.splitext(output_path)
        if not ext:
            ext = ".pdf"
        for suffix, fig in figures.items():
            path = f"{stem}_{suffix}{ext}"
            fig.savefig(path, dpi=300 if ext == ".png" else None,
                        bbox_inches="tight")
            print(f"  Saved: {path}")
    else:
        plt.show()

    return figures["linear"], figures["log"]


# ══════════════════════════════════════════════════════════════════════════════
# MAIN
# ══════════════════════════════════════════════════════════════════════════════

def main(filepath, geom=None, truncation=0.70, n_bins=40,
         min_steps=5, output_path=None):
    print(f"\nLoading TPC data from:  {filepath}")
    data = load_tpc_data(filepath)
    print(f"  Tree contains {data['n_events']:,} events")

    use_dedx = geom is not None and geom.is_complete
    mode_str = "dE/dx [MeV/mm]" if use_dedx else "sum(Edep) [MeV]"
    print(f"  Energy-loss mode: {mode_str}")

    print("\nComputing per-event energy loss...")
    values = compute_dedx(data, geom=geom, min_steps=min_steps)
    print(f"  {len(values):,} events accepted (min_steps >= {min_steps})")

    print(f"\nFitting Landau (Moyal) to bottom {truncation*100:.0f}%...")
    result = fit_landau(values, truncation=truncation, n_bins=n_bins)

    print(f"\n{'─'*52}")
    print(f"  FIT RESULTS")
    print(f"{'─'*52}")
    print(f"  MPV                  : {result['mpv']:.4g} MeV")
    print(f"  Width xi             : {result['scale']:.4g} MeV")
    print(f"  Fit region chi2/ndf  : {result['chi2_red']:.2f}"
          f"    (p = {result['p_value']:.3f})")
    print(f"{'─'*52}\n")

    fig_lin, fig_log = plot_dedx(result, use_dedx=use_dedx,
                                  output_path=output_path)
    plt.close(fig_lin)
    plt.close(fig_log)


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python tpc_dedx_analysis.py <file.root> [output.pdf]")
        sys.exit(1)

    root_file  = sys.argv[1]
    fig_output = sys.argv[2] if len(sys.argv) > 2 else None
    
    main(root_file, geom=None, output_path=fig_output)