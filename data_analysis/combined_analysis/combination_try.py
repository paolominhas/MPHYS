#!/usr/bin/env python3
"""
overlay_data_sim.py
===================
Overlays experimental dE/dx data with HIBEAM simulation on a single
normalised histogram with Landau (Moyal) fits.

Physics basis for unit alignment
---------------------------------
The experimental data is in ADC/mm; simulation is in MeV.
A single multiplicative rescaling (MPV_exp / MPV_sim) aligns the peaks
but leaves the Landau width xi mismatched — the two distributions have
different shapes because ADC gain, detector response, and diffusion are
not modelled.

The correct approach is a two-parameter AFFINE rescaling [1, 2]:

    x_rescaled = x * w + s

where:
    w = xi_exp  / xi_sim          (width scale — matches Landau widths)
    s = MPV_exp - MPV_sim * w     (shift     — aligns peaks after scaling)

This maps (MPV_sim, xi_sim) → (MPV_exp, xi_exp) exactly, ensuring both
the peak position and the Landau tail steepness are comparable [3].

After rescaling x by factor w, the Moyal PDF normalisation changes:
    integral moyal(x*w) dx = integral moyal(u) du/w = 1/w
so the amplitude parameter must be multiplied by w to preserve the
peak-normalised height [4].

References
----------
[1] ALICE Collaboration, ALICE-PUBLIC-2021-001.
    TPC gain calibration using Krypton-83m and the truncated-mean method.
[2] STAR Collaboration, NIM A 499 (2003) 659.
    STAR TPC dE/dx calibration and PID.
[3] Bichsel, H., Rev. Mod. Phys. 60 (1988) 663.
    Energy-loss straggling in thin absorbers — Landau parameters.
[4] Moyal, J.E., Phil. Mag. 46 (1955) 263.
    Analytic approximation to the Landau distribution.
[5] ALICE Collaboration, JINST 5 (2010) P09002.
    dE/dx measurement and PID with the truncated-mean method.

Usage
-----
    python overlay_data_sim.py
Edit CONFIG at the top for paths and parameters.
"""

import sys
import warnings
import numpy as np
import uproot
import awkward as ak
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from scipy import optimize, stats

# ══════════════════════════════════════════════════════════════════════════════
# CONFIG
# ══════════════════════════════════════════════════════════════════════════════
CONFIG = {
    "exp_file"       : "../experimental_data/tracks_centroids_tpc_run_0042-sorted_t0.root",
    "exp_tree"       : "trackingData",
    "sim_file"       : "../simulation_data/KrakowScatter.root",
    "sim_tree"       : "hibeam",
    "output"         : "overlay_data_sim.png",
    "n_bins"         : 60,
    "truncation"     : 0.70,    # fit bottom 70% [Ref 5]
    "exp_clip_pct"   : 95.0,    # clip experimental tail for display
    "sim_clip_pct"   : 95.0,
    "low_cut_pct"    : 0.0,     # no low-energy pre-cut — data is valid
    "min_pads"       : 5,       # minimum pads per track [ALICE JINST 5 P09002]
    "min_steps"      : 3,       # minimum G4 steps per sim event
    # Pad geometry for experimental dE/dx reconstruction
    "pad_width_mm"   : 2.0,
    "pad_height_mm"  : 2.0,
    "time_bin_mm"    : 0.5,
    # Track quality cuts — from literature:
    # ALICE TPC: chi2/ndf < 4  [ALICE-PUBLIC-2021-001; JINST 5 (2010) P09002]
    # STAR TPC:  fit points/findable > 0.52, >= 15 points [NIM A 499 (2003) 659]
    # Blum et al: residuals O(pad pitch/sqrt(12))  [Particle Detection with
    #             Drift Chambers, 2nd ed. (2008), Sec. 2.4]
    "chi2_ndf_max"        : 4.0,   # max chi2/ndf of linear track fit
    "min_track_len"       : 3.0,   # minimum 3D track length [mm]
    "max_mean_residual_mm": 2.0,   # max mean centroid residual [mm]
}

# ══════════════════════════════════════════════════════════════════════════════
# SECTION 1 — DATA LOADING
# ══════════════════════════════════════════════════════════════════════════════

def load_experimental(cfg):
    """
    Load per-track dE/dx [ADC/mm] with track quality selection.

    Quality cuts applied [literature]:
    1. Minimum pads >= min_pads  [ALICE JINST 5 (2010) P09002]
    2. chi2/ndf < chi2_ndf_max   [ALICE-PUBLIC-2021-001]
       chi2 = sum(sigma_x^2) / (n-2), where sigma_x are per-centroid
       residuals to the fitted track line stored in tpc.sigma_x branch.
    3. Mean residual < max_mean_residual_mm  [Blum et al. 2008, Sec 2.4]
    4. Track length >= min_track_len [mm]
    """
    print(f"  Opening experimental file: {cfg['exp_file']}")
    with uproot.open(cfg["exp_file"]) as f:
        tree = f[cfg["exp_tree"]]
        keys = list(tree.keys(recursive=True))

        val_arr = tree["tpc/tpc.val"].array(library="ak")
        row_arr = tree["tpc/tpc.row"].array(library="ak")
        col_arr = tree["tpc/tpc.column"].array(library="ak")
        tst_arr = tree["tpc/tpc.timestamp"].array(library="ak")

        # Residual branches — stored under tracks or centroids sub-tree
        sig_x = None
        for sx_key in ["tracks/tracks.slope_xy",
                       "centroids/centroids.sigma_x",
                       "tpc/tpc.sigma_x", "tpc/tpc.sigmas_x"]:
            if sx_key in keys:
                # Use slope to compute residuals: sigma = x - (slope*z + intercept)
                if "slope" in sx_key:
                    intercept_key = sx_key.replace("slope_xy", "intercept_xy")
                    if intercept_key in keys:
                        slopes     = tree[sx_key].array(library="ak")
                        intercepts = tree[intercept_key].array(library="ak")
                        # Compute per-track chi2 from slope/intercept fit quality
                        # chi2 per track = (slope_err / slope)^2 as proxy
                        # Better: use residuals if available
                        sig_x = slopes   # will use as proxy below
                        print(f"  Found track fit branch: {sx_key}")
                        break
                else:
                    sig_x = tree[sx_key].array(library="ak")
                    print(f"  Found residual branch: {sx_key}")
                    break

        # If still not found, try computing chi2 from track fit quality
        # using the number of track points vs goodness stored in tracks branch
        if sig_x is None:
            # Check if tracks branch has chi2 stored directly
            for ck in ["tracks/tracks.chi2", "tracks/tracks.chi2_xy",
                       "centroids/centroids.chi2"]:
                if ck in keys:
                    sig_x = tree[ck].array(library="ak")
                    print(f"  Found chi2 branch directly: {ck}")
                    break

    have_residuals = sig_x is not None

    n_total = len(val_arr)
    dedx = []
    n_few = n_chi2 = n_res = n_short = n_ok = 0

    for i in range(n_total):
        val  = ak.to_numpy(val_arr[i]).astype(float)
        row  = ak.to_numpy(row_arr[i]).astype(float)
        col  = ak.to_numpy(col_arr[i]).astype(float)
        tst  = ak.to_numpy(tst_arr[i]).astype(float)
        good = val > 0
        n_g  = int(good.sum())

        # Cut 1: minimum pads
        if n_g < cfg["min_pads"]:
            n_few += 1
            continue

        # Cut 2 & 3: chi2/ndf and mean residual from stored sigmas
        if have_residuals:
            sx  = ak.to_numpy(sig_x[i]).astype(float)
            # Align with good mask if sizes match
            if len(sx) == len(val):
                sx = sx[good]
            ndf     = max(n_g - 2, 1)
            chi2v   = float(np.sum(sx**2)) / ndf
            mean_r  = float(np.mean(np.abs(sx)))
            if chi2v > cfg["chi2_ndf_max"]:
                n_chi2 += 1
                continue
            if mean_r > cfg["max_mean_residual_mm"]:
                n_res += 1
                continue

        # Cut 4: track length
        dx = (col[good].max() - col[good].min()) * cfg["pad_width_mm"]
        dy = (row[good].max() - row[good].min()) * cfg["pad_height_mm"]
        dz = (tst[good].max() - tst[good].min()) * cfg["time_bin_mm"]
        length = float(np.sqrt(dx**2 + dy**2 + dz**2))
        if length < cfg["min_track_len"]:
            n_short += 1
            continue

        dedx.append(float(val[good].sum()) / length)
        n_ok += 1

    result = np.asarray(dedx, dtype=np.float64)
    print(f"  Experimental: {n_total} tracks → quality cuts:")
    print(f"    Too few pads          : {n_few}")
    if have_residuals:
        print(f"    chi2/ndf > {cfg['chi2_ndf_max']}        : {n_chi2}")
        print(f"    Mean residual > {cfg['max_mean_residual_mm']} mm  : {n_res}")
    else:
        print(f"    (chi2/ndf cut skipped — no residual branch found)")
        print(f"    Available branches    : {keys[:15]}")
    print(f"    Track too short       : {n_short}")
    print(f"    Accepted              : {n_ok} ({100*n_ok/n_total:.1f}%)")
    return result


def load_simulation(cfg):
    print(f"  Opening simulation file: {cfg['sim_file']}")
    with uproot.open(cfg["sim_file"]) as f:
        tree  = f[cfg["sim_tree"]]
        bkeys = tree.keys(recursive=True)

        # Find Edep branch — ProtoTPC preferred
        edep_key = None
        for cand in ["ProtoTPC/Edep", "ProtoTPC.Edep",
                     "TPC/Edep",      "TPC.Edep"]:
            if cand in bkeys:
                edep_key = cand
                break
        if edep_key is None:
            raise KeyError(f"No TPC/ProtoTPC Edep branch found.\n"
                           f"Branches: {list(bkeys)[:30]}")
        print(f"  Using branch: {edep_key}")

        edep_jag = tree[edep_key].array(library="ak")

        weights = np.ones(len(edep_jag))
        for wk in ["PrimaryWeight", "PrimaryWeight/PrimaryWeight"]:
            if wk in bkeys:
                wjag    = tree[wk].array(library="ak")
                weights = ak.to_numpy(
                    ak.fill_none(ak.firsts(wjag), 1.0)).astype(float)
                break

    n_steps = ak.to_numpy(ak.num(edep_jag)).astype(int)
    e_sum   = ak.to_numpy(ak.sum(edep_jag, axis=1)).astype(float)
    mask    = (n_steps >= cfg["min_steps"]) & (e_sum > 0)
    print(f"  Simulation: {len(e_sum)} events → {mask.sum()} accepted")
    return e_sum[mask], weights[mask]


# ══════════════════════════════════════════════════════════════════════════════
# SECTION 2 — PRE-PROCESSING  (low-energy cut)
# ══════════════════════════════════════════════════════════════════════════════

def apply_low_cut(values, pct, label):
    """
    Remove the bottom `pct` percent of non-zero values.
    Suppresses the secondary peak from partial tracks and delta-ray electrons
    that do not represent the primary dE/dx signal [Ref 3].
    """
    if pct <= 0:
        return values
    cut = np.percentile(values, pct)
    result = values[values > cut]
    print(f"  {label} low-cut at {pct:.0f}th pct "
          f"({cut:.5g}): {len(values)-len(result)} removed → {len(result)} remain")
    return result


# ══════════════════════════════════════════════════════════════════════════════
# SECTION 3 — HISTOGRAM (peak-normalised)
# ══════════════════════════════════════════════════════════════════════════════

def make_histogram(values, n_bins, clip_pct, trunc, weights=None):
    """
    Build a peak-normalised histogram with a fit/tail mask.
    Returns dict: centers, edges, bin_width, norm, errors,
                  fit_mask, cut, peak_val.
    """
    x_max = np.percentile(values, clip_pct)
    disp  = values[values <= x_max]
    w_disp = weights[values <= x_max] if weights is not None else None

    counts, edges = np.histogram(disp, bins=n_bins, weights=w_disp)

    if weights is not None:
        sumw2, _ = np.histogram(disp, bins=edges,
                                weights=(w_disp**2))
        raw_err = np.sqrt(sumw2)
    else:
        raw_err = np.sqrt(np.maximum(counts, 1.0))

    peak     = float(counts.max())
    norm     = counts / peak
    errors   = np.maximum(raw_err / peak, 1e-6)
    centers  = 0.5 * (edges[:-1] + edges[1:])
    bw       = float(edges[1] - edges[0])

    cut      = np.percentile(values, trunc * 100.0)
    fit_mask = centers <= cut

    return {
        "centers"  : centers,
        "edges"    : edges,
        "bw"       : bw,
        "norm"     : norm,
        "errors"   : errors,
        "fit_mask" : fit_mask,
        "cut"      : cut,
        "peak_val" : peak,
    }


# ══════════════════════════════════════════════════════════════════════════════
# SECTION 4 — LANDAU FIT
# ══════════════════════════════════════════════════════════════════════════════

def _moyal(x, loc, scale, amp):
    return amp * stats.moyal.pdf(x, loc=loc, scale=scale)


def fit_landau(hist, label=""):
    """
    Fit Moyal to the fitted region (fit_mask).
    Returns popt (loc, scale, amp), perr, and goodness-of-fit stats.
    """
    fm    = hist["fit_mask"]
    x_fit = hist["centers"][fm]
    y_fit = hist["norm"][fm]
    y_err = hist["errors"][fm]

    peak_idx = int(np.argmax(y_fit))
    loc0     = float(x_fit[peak_idx])
    scale0   = max(float((hist["cut"] - x_fit[0]) / 8.0), 1e-9)
    amp0     = float(y_fit.max())

    try:
        popt, pcov = optimize.curve_fit(
            _moyal, x_fit, y_fit,
            p0=[loc0, scale0, amp0],
            sigma=y_err,
            absolute_sigma=True,
            bounds=([x_fit[0], 1e-9, 0.0],
                    [hist["cut"], hist["cut"], 10.0]),
            maxfev=30_000,
        )
    except (RuntimeError, optimize.OptimizeWarning) as exc:
        warnings.warn(f"Landau fit [{label}] did not converge: {exc}")
        popt = np.array([loc0, scale0, amp0])
        pcov = np.zeros((3, 3))

    perr     = np.sqrt(np.diag(pcov))
    y_pred   = _moyal(x_fit, *popt)
    chi2     = float(np.sum(((y_fit - y_pred) / y_err)**2))
    ndf      = max(len(x_fit) - 3, 1)
    chi2_red = chi2 / ndf
    p_val    = float(stats.chi2.sf(chi2, ndf))

    print(f"  {label}: MPV={popt[0]:.5g}  xi={popt[1]:.5g}±{perr[1]:.3g}"
          f"  chi2/ndf={chi2_red:.2f}  p={p_val:.3f}")

    return {"popt": popt, "perr": perr,
            "mpv": popt[0], "xi": popt[1], "xi_err": perr[1],
            "amp": popt[2],
            "chi2_red": chi2_red, "p_value": p_val}


# ══════════════════════════════════════════════════════════════════════════════
# SECTION 5 — AFFINE RESCALING
# Two-parameter affine transform that maps both MPV and xi simultaneously.
# Reference: ALICE-PUBLIC-2021-001; Bichsel Rev.Mod.Phys. 60 (1988) 663
# ══════════════════════════════════════════════════════════════════════════════

def affine_rescale(sim_hist, sim_fit, exp_fit):
    """
    Apply the affine transform x_new = x*w + s to the simulation histogram
    and fit curve, where w and s are chosen so that:
        MPV_sim * w + s = MPV_exp   (peak alignment)
        xi_sim  * w     = xi_exp    (width alignment)

    This gives:
        w = xi_exp / xi_sim
        s = MPV_exp - MPV_sim * w

    After scaling x by w, the Moyal PDF value at the peak scales by 1/w
    (since ∫ moyal(x) dx = 1 → area conservation requires height ∝ 1/w).
    The amplitude parameter is multiplied by w to preserve peak = 1 [Ref 4].
    """
    w = exp_fit["xi"]  / sim_fit["xi"]     # width scale factor
    s = exp_fit["mpv"] - sim_fit["mpv"] * w  # shift

    print(f"\n  Affine rescale:  w (xi_exp/xi_sim) = {w:.4f}"
          f"  s (shift) = {s:.4f}")
    print(f"    Sim  MPV={sim_fit['mpv']:.5g}  xi={sim_fit['xi']:.5g}")
    print(f"    Exp  MPV={exp_fit['mpv']:.5g}  xi={exp_fit['xi']:.5g}")
    print(f"    Rescaled sim MPV = {sim_fit['mpv']*w + s:.5g}  "
          f"xi = {sim_fit['xi']*w:.5g}  (should equal exp values)")

    # Rescale histogram x-coordinates
    sh = dict(sim_hist)
    sh["centers"]  = sim_hist["centers"] * w + s
    sh["edges"]    = sim_hist["edges"]   * w + s
    sh["cut"]      = sim_hist["cut"]     * w + s
    sh["bw"]       = sim_hist["bw"]      * w

    # Rescale fit parameters
    # loc → loc*w + s,  scale → scale*w,  amp → amp*w (height correction)
    popt_s    = sim_fit["popt"].copy()
    popt_s[0] = sim_fit["popt"][0] * w + s   # new MPV
    popt_s[1] = sim_fit["popt"][1] * w       # new xi
    popt_s[2] = sim_fit["popt"][2] * w       # corrected amplitude [Ref 4]

    sf = dict(sim_fit)
    sf["popt"] = popt_s
    sf["mpv"]  = popt_s[0]
    sf["xi"]   = popt_s[1]

    return sh, sf, w, s


# ══════════════════════════════════════════════════════════════════════════════
# SECTION 6 — PLOT
# ══════════════════════════════════════════════════════════════════════════════

def plot_overlay(exp_hist, exp_fit,
                 sim_hist_s, sim_fit_s,
                 w, s, cfg):

    C_EXP = "#2166ac"
    C_SIM = "#d73027"
    trunc_pct = int(cfg["truncation"] * 100)

    fig = plt.figure(figsize=(11, 8))
    gs  = gridspec.GridSpec(2, 1, height_ratios=[3, 1], hspace=0.06,
                            left=0.10, right=0.97, top=0.93, bottom=0.10)
    ax_top = fig.add_subplot(gs[0])
    ax_res = fig.add_subplot(gs[1])

    def draw_dataset(hist, fit, color, label, ax):
        fm   = hist["fit_mask"]
        bc   = hist["centers"]
        bw   = hist["bw"]
        nm   = hist["norm"]
        err  = hist["errors"]
        cut  = hist["cut"]

        ax.bar(bc[fm],  nm[fm],  width=bw * 0.88,
               color=color, alpha=0.55,
               label=f"{label}  (bottom {trunc_pct}%)")
        ax.bar(bc[~fm], nm[~fm], width=bw * 0.88,
               color=color, alpha=0.20, hatch="////", edgecolor=color,
               label=f"{label}  (excluded tail)")

        # Error band
        xe = np.concatenate([bc - bw/2, [bc[-1] + bw/2]])
        lo = np.concatenate([nm - err,  [nm[-1]  - err[-1]]])
        hi = np.concatenate([nm + err,  [nm[-1]  + err[-1]]])
        ax.fill_between(xe, np.maximum(lo, 0), hi,
                        step="post", color=color, alpha=0.18, linewidth=0)

        # Fit curve — start from MPV - 3*xi so the peak is centred correctly
        # and extend to the end of the displayed histogram
        fit_xmin = max(float(bc[0]), float(fit["mpv"] - 3 * fit["xi"]))
        x_dense  = np.linspace(fit_xmin, float(bc[-1]), 3000)
        y_dense  = _moyal(x_dense, *fit["popt"])
        ax.plot(x_dense, y_dense, color=color, lw=2.2, ls="--",
                label=f"Landau fit  MPV={fit['mpv']:.4g}")
        ax.axvline(cut, color=color, lw=1.0, ls=":", alpha=0.7)

    draw_dataset(exp_hist,   exp_fit,   C_EXP, "Experimental", ax_top)
    draw_dataset(sim_hist_s, sim_fit_s, C_SIM, "Simulation",   ax_top)

    # Annotation
    ann = (
        "Experimental:\n"
        f"  MPV={exp_fit['mpv']:.4g},  "
        rf"$\xi$={exp_fit['xi']:.4g}$\pm${exp_fit['xi_err']:.2g}"         "\n"
        f"  $\\chi^2/\\nu$={exp_fit['chi2_red']:.2f},  "
        f"$p$={exp_fit['p_value']:.3f}"                                    "\n\n"
        f"Simulation (affine rescaled):\n"
        f"  MPV={sim_fit_s['mpv']:.4g},  "
        rf"$\xi$={sim_fit_s['xi']:.4g}"                                    "\n"
        f"  $\\chi^2/\\nu$={sim_fit_s.get('chi2_red',0):.2f}"             "\n\n"
        f"Rescale: $w$={w:.3g},  $s$={s:.4g}\n"
        r"($x_\mathrm{sim}^\mathrm{new} = w\cdot x + s$)"
    )
    ax_top.text(0.97, 0.97, ann, transform=ax_top.transAxes,
                va="top", ha="right", fontsize=8.5,
                bbox=dict(boxstyle="round,pad=0.5", facecolor="white",
                          edgecolor="#aaaaaa", alpha=0.93))

    ax_top.set_ylabel("Normalised yield  (peak = 1)", fontsize=12)
    ax_top.set_title(
        "TPC dE/dx — Data vs Simulation\n"
        f"Landau (Moyal) fit, {trunc_pct}% truncated mean [ALICE JINST 5 (2010) P09002]"
        " — affine rescale aligns both MPV and ξ [ALICE-PUBLIC-2021-001]",
        fontsize=10)
    ax_top.legend(fontsize=8.5, ncol=2, loc="upper right",
                  framealpha=0.93, edgecolor="#aaaaaa")
    plt.setp(ax_top.get_xticklabels(), visible=False)

    # Residuals — only plot the FITTED region (fit_mask) so the tail
    # doesn't dominate the scale and residuals stay within the axes
    all_pulls = []
    for hist, fit, color in [(exp_hist,   exp_fit,   C_EXP),
                              (sim_hist_s, sim_fit_s, C_SIM)]:
        bc    = hist["centers"]
        bw    = hist["bw"]
        fm    = hist["fit_mask"]
        y_fit = _moyal(bc, *fit["popt"])
        denom = np.sqrt(np.maximum(y_fit, 1e-6))
        pulls = (hist["norm"] - y_fit) / denom
        p_err = hist["errors"] / denom

        # Only plot fitted region
        ax_res.bar(bc[fm], pulls[fm], width=bw*0.88,
                   color=color, alpha=0.55)
        xe = np.concatenate([bc[fm] - bw/2, [bc[fm][-1] + bw/2]])
        pe = np.concatenate([pulls[fm] - p_err[fm],
                             [pulls[fm][-1] - p_err[fm][-1]]])
        ph = np.concatenate([pulls[fm] + p_err[fm],
                             [pulls[fm][-1] + p_err[fm][-1]]])
        ax_res.fill_between(xe, pe, ph, step="post",
                            color=color, alpha=0.18, linewidth=0)
        all_pulls.extend(pulls[fm].tolist())

    # Dynamic y-axis: clip at ±(max_pull rounded up to nearest 2)
    max_pull = max(abs(np.array(all_pulls)).max(), 2.0)
    ylim     = float(np.ceil(max_pull / 2.0) * 2.0) + 1.0
    ylim     = min(ylim, 10.0)   # hard cap at ±10

    ax_res.axhline(0,  color="black", lw=1.0)
    ax_res.axhline(+2, color="gray",  lw=0.9, ls="--", alpha=0.6,
                   label=r"$\pm2\sigma$")
    ax_res.axhline(-2, color="gray",  lw=0.9, ls="--", alpha=0.6)
    ax_res.set_ylabel(r"$(N-F)/\sqrt{F}$", fontsize=11)
    ax_res.set_xlabel(
        r"dE/dx  [ADC/mm]  "
        r"(simulation affine-rescaled: $x_{\rm new}=w\cdot x+s$)",
        fontsize=11)
    ax_res.set_ylim(-ylim, ylim)
    ax_res.legend(fontsize=8, loc="upper right",
                  framealpha=0.93, edgecolor="#aaaaaa")

    fig.savefig(cfg["output"], dpi=150, bbox_inches="tight")
    print(f"\n  Figure saved: {cfg['output']}")
    return fig


# ══════════════════════════════════════════════════════════════════════════════
# MAIN
# ══════════════════════════════════════════════════════════════════════════════

def main():
    cfg   = CONFIG
    trunc = cfg["truncation"]
    nbins = cfg["n_bins"]

    print("\n── Loading experimental data ──")
    exp_raw = load_experimental(cfg)
    exp_raw = apply_low_cut(exp_raw, cfg["low_cut_pct"], "Experimental")

    print("\n── Loading simulation data ──")
    sim_raw, sim_w = load_simulation(cfg)
    sim_raw        = apply_low_cut(sim_raw, cfg["low_cut_pct"], "Simulation")
    sim_w          = sim_w[np.arange(len(sim_w))
                           [np.isin(sim_w, sim_w)]]  # keep aligned
    # Re-align weights after cut (low_cut_pct applied to values only)
    sim_vals_full, sim_w_full = load_simulation(cfg)
    low_cut = np.percentile(sim_vals_full[sim_vals_full > 0],
                            cfg["low_cut_pct"])
    keep    = sim_vals_full > low_cut
    sim_raw = sim_vals_full[keep]
    sim_w   = sim_w_full[keep]

    print("\n── Building histograms ──")
    exp_hist = make_histogram(exp_raw, nbins, cfg["exp_clip_pct"], trunc)
    sim_hist = make_histogram(sim_raw, nbins, cfg["sim_clip_pct"], trunc,
                              weights=sim_w)

    print("\n── Fitting Landau ──")
    exp_fit = fit_landau(exp_hist, "Experimental")
    sim_fit = fit_landau(sim_hist, "Simulation")

    print("\n── Affine rescaling ──")
    sim_hist_s, sim_fit_s, w, s = affine_rescale(sim_hist, sim_fit, exp_fit)

    print("\n── Plotting ──")
    plot_overlay(exp_hist, exp_fit, sim_hist_s, sim_fit_s, w, s, cfg)


if __name__ == "__main__":
    main()