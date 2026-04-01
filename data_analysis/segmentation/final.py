#!/usr/bin/env python3
"""
tpc_dedx_analysis.py
====================
Extracts TPC energy-loss data from a processed HIBEAM ROOT ntuple,
fits a Landau (Moyal) distribution to the bottom 70 % of values using
the truncated-mean method, and produces a plot with residuals and
goodness-of-fit statistics.

The load_tpc_data() function is intentionally standalone so it can be
imported cleanly by other analysis scripts.

Usage
-----
    python tpc_dedx_analysis.py /path/to/HIBEAMScatter.root [figure.png]

To enable true dE/dx [MeV/mm] rather than sum(Edep) [MeV], supply the
TPC geometry constants (see TPCGeometry and main() at the bottom).

References
----------
[1] Landau, L.D. (1944). J. Phys. USSR 8, 201.
    Original straggling distribution for thin absorbers.
[2] Moyal, J.E. (1955). Phil. Mag. 46, 263.
    Analytic approximation to the Landau distribution; implemented
    in scipy as scipy.stats.moyal.
[3] Bichsel, H. (1988). Rev. Mod. Phys. 60, 663.
    Definitive review of energy-loss straggling in detectors.
[4] ALICE Collaboration (2010). JINST 5, P09002.
    dE/dx measurement and PID with the truncated-mean method in a TPC.
    This paper defines the standard truncation fraction (~70 %).
[5] STAR Collaboration (2003). Nucl. Instrum. Meth. A 499, 659.
    Additional reference for TPC dE/dx truncated-mean methodology.
[6] Particle Data Group (2022). Prog. Theor. Exp. Phys. 2022, 083C01.
    Statistics chapter — chi-squared goodness-of-fit, pulls, p-values.
"""

import sys
import warnings
from dataclasses import dataclass
from typing import Dict, Optional

import awkward as ak
import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
import numpy as np
import uproot
from scipy import optimize, stats


# ══════════════════════════════════════════════════════════════════════════════
# SECTION 1 — TPC GEOMETRY
# Update these constants to match the values in TPCHit.hh.
# Leave as None to fall back to using sum(Edep) per event [MeV].
# ══════════════════════════════════════════════════════════════════════════════

@dataclass
class TPCGeometry:
    """
    Physical TPC parameters needed to convert pad indices → mm.
    All values must match those hard-coded in TPCHit.hh.
    If any field is None the geometry is considered incomplete and
    track-length normalisation is skipped.
    """
    pad_width:  Optional[float] = None   # mm  (padWidth  in TPCHit.hh)
    pad_height: Optional[float] = None   # mm  (padHeight in TPCHit.hh)
    drift_v:    Optional[float] = None   # mm/ns (driftV in TPCHit.hh)
    time_bin:   Optional[float] = None   # ns  (timeBin  in TPCHit.hh)

    @property
    def is_complete(self) -> bool:
        return all(v is not None for v in
                   [self.pad_width, self.pad_height, self.drift_v, self.time_bin])


# ══════════════════════════════════════════════════════════════════════════════
# SECTION 2 — DATA LOADING  (import this function in other scripts)
# ══════════════════════════════════════════════════════════════════════════════

def load_tpc_data(filepath: str, tree_name: str = "hibeam") -> Dict:
    """
    Load TPC branch data from a processed HIBEAM ROOT ntuple.

    Branch discovery is done dynamically so that minor naming differences
    between ROOT split-class formats are handled without hardcoding.

    Parameters
    ----------
    filepath : str
        Path to ROOT file produced by hibeam_ana.
    tree_name : str
        TTree name inside the file (default "hibeam").

    Returns
    -------
    dict
        'edep'      : ak.Array — per-event lists of G4-step Edep values [MeV]
        'pad_row'   : ak.Array or None — per-event pad row indices
        'pad_col'   : ak.Array or None — per-event pad column indices
        'timestamp' : ak.Array or None — per-event drift timestamps
        'n_el'      : ak.Array or None — per-event ionisation electron counts
        'n_events'  : int

    Raises
    ------
    KeyError
        If the TPC branch or its Edep sub-branch cannot be found.
    OSError
        If the ROOT file cannot be opened by uproot.
    """
    with uproot.open(filepath) as root_file:

        # ── Locate the TTree ────────────────────────────────────────────────
        available = list(root_file.keys())
        if tree_name not in root_file:
            raise KeyError(
                f"TTree '{tree_name}' not found in {filepath}.\n"
                f"Available keys: {available}"
            )
        tree = root_file[tree_name]

        # ── Discover TPC sub-branch keys ────────────────────────────────────
        # uproot represents split ROOT object branches as "Parent/Member"
        # or "Parent.Member" depending on ROOT version; we search for both.
        all_keys = tree.keys(recursive=True)

        # Collect keys belonging to the TPC branch but not ProtoTPC
        tpc_keys = [
            k for k in all_keys
            if ("TPC" in k)
            and not k.startswith("ProtoTPC")
            and not k.startswith("target")
        ]

        if not tpc_keys:
            raise KeyError(
                "No TPC sub-branches found in the tree.\n"
                f"All branch keys: {all_keys}"
            )

        def find_branch(suffix: str) -> Optional[str]:
            """Return the full key whose last component matches suffix."""
            for k in tpc_keys:
                if k.split("/")[-1] == suffix or k.split(".")[-1] == suffix:
                    return k
            return None

        edep_key = find_branch("Edep")
        if edep_key is None:
            raise KeyError(
                f"TPC Edep sub-branch not found.\nTPC keys found: {tpc_keys}"
            )

        # ── Load branches — warn rather than crash for optional ones ─────────
        def safe_array(suffix: str) -> Optional[ak.Array]:
            key = find_branch(suffix)
            if key is None:
                warnings.warn(
                    f"TPC sub-branch '{suffix}' not found; "
                    "track-length calculation will be unavailable.",
                    stacklevel=3,
                )
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

def _track_length_mm(
    pad_row, pad_col, timestamp, geom: TPCGeometry
) -> float:
    """
    Estimate straight-track length [mm] from the 3-D extent of pad hits.

    Uses the distance between extreme pad positions as a lower bound on
    the true path length. Valid for straight (non-curling) tracks through
    the TPC active volume.

    Physical coordinates are recovered as:
        x = pad_column * pad_width   [mm]
        y = pad_row    * pad_height  [mm]
        z = timestamp  * drift_v * time_bin  [mm]
    """
    if len(pad_row) < 2:
        return float("nan")
    dx = (float(ak.max(pad_col))   - float(ak.min(pad_col)))   * geom.pad_width
    dy = (float(ak.max(pad_row))   - float(ak.min(pad_row)))   * geom.pad_height
    dz = (float(ak.max(timestamp)) - float(ak.min(timestamp))) * geom.drift_v * geom.time_bin
    return float(np.sqrt(dx**2 + dy**2 + dz**2))


def compute_dedx(
    data: Dict,
    geom: Optional[TPCGeometry] = None,
    min_steps: int = 5,
) -> np.ndarray:
    """
    Compute per-event energy-loss values from TPC data.

    If a complete TPCGeometry is supplied AND pad position branches are
    present, returns dE/dx [MeV/mm].

    Otherwise returns total deposited energy sum(Edep) [MeV] per event.
    For a parallel beam traversing a fixed detector depth this is directly
    proportional to dE/dx [Ref 3].

    Parameters
    ----------
    data : dict
        Output of load_tpc_data().
    geom : TPCGeometry or None
        TPC geometry constants.  None → use sum(Edep) only.
    min_steps : int
        Minimum number of G4 steps required to accept an event.
        Rejects events with no or negligible TPC activity.

    Returns
    -------
    np.ndarray, shape (n_accepted_events,)
        Energy loss [MeV/mm or MeV] for each accepted event.
    """
    edep    = data["edep"]
    pad_row = data.get("pad_row")
    pad_col = data.get("pad_col")
    ts      = data.get("timestamp")

    use_geom = (
        geom is not None
        and geom.is_complete
        and pad_row is not None
        and pad_col is not None
        and ts is not None
    )

    result = []
    for i in range(data["n_events"]):
        ev_edep = edep[i]
        if len(ev_edep) < min_steps:
            continue
        e_total = float(ak.sum(ev_edep))

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

def _moyal_scaled(x: np.ndarray, loc: float, scale: float, norm: float) -> np.ndarray:
    """
    Moyal PDF multiplied by a free normalisation constant.

    The Moyal distribution [Ref 2] is the standard analytic approximation
    to the Landau distribution [Ref 1] used throughout HEP.  scipy.stats.moyal
    implements it exactly.  The mode (most probable value) of this
    parameterisation is at x = loc; the width parameter 'scale' corresponds
    to the Landau width ξ.

    norm = N_events * bin_width ensures the curve has the same integral
    as the histogram.
    """
    return norm * stats.moyal.pdf(x, loc=loc, scale=scale)


def fit_landau(
    values: np.ndarray,
    truncation: float = 0.70,
    n_bins: int = 100,
) -> Dict:
    """
    Fit a Landau (Moyal) distribution to the lower `truncation` fraction
    of energy-loss values using binned Pearson chi-squared minimisation.

    The truncated-mean method is the standard TPC dE/dx technique [Ref 4, 5]:
    the rising edge and peak of the Landau are well described by the
    single-collision model, while the high-energy Landau tail is inflated
    by delta-ray production and is intentionally excluded from the fit.

    Parameters
    ----------
    values : np.ndarray
        Per-event energy loss values.
    truncation : float
        Fraction of data (from the low end) included in the fit (default 0.70).
    n_bins : int
        Number of histogram bins spanning the full data range.

    Returns
    -------
    dict
        Fit parameters, uncertainties, goodness-of-fit statistics for both
        the fit region and the excluded tail, plus arrays for plotting.

    Notes on goodness-of-fit [Ref 6]
    ----------------------------------
    chi2     : Pearson chi-squared in the fit region, bins with ≥ 5 expected
    ndf      : degrees of freedom = (bins used) − 3 free parameters
    chi2_red : reduced chi-squared; values ~1 indicate a good fit
    p_value  : probability of observing a chi2 at least this large by chance
    chi2_tail / ndf_tail / p_tail :
               same quantities computed for the EXCLUDED tail region,
               quantifying how much the data there deviates from the fit
               extrapolated into that region.  A small p_tail is expected
               and physically meaningful (delta-ray enhancement) [Ref 3, 4].
    """
    if len(values) < 50:
        raise ValueError(
            f"Only {len(values)} events — need at least 50 for a reliable fit."
        )

    cut_value = np.percentile(values, truncation * 100.0)

    # ── Build histogram — guarantee ≥ n_bins bins inside the fit region ───────
    # The Landau tail is extreme: binning over the full range would compress
    # the entire peak into 1-2 bins.  Strategy:
    #   1. Cap display at the 99th percentile (tail above is unphysically long
    #      for a single-species sample and swamps the x-axis).
    #   2. Choose bin width so that exactly n_bins span [hist_min, cut_value].
    #      The same bin width extends uniformly to hist_max, guaranteeing
    #      density where the fit lives.
    hist_min = float(values.min())
    hist_max = float(np.percentile(values, 99.0))
    fit_range = cut_value - hist_min
    if fit_range <= 0:
        raise ValueError(
            "The truncation cut equals the data minimum — data may be degenerate."
        )
    bin_width_target = fit_range / n_bins
    n_total = max(int(np.ceil((hist_max - hist_min) / bin_width_target)), n_bins + 1)
    counts, bin_edges = np.histogram(values, bins=n_total,
                                     range=(hist_min, hist_max))
    bin_centers = 0.5 * (bin_edges[:-1] + bin_edges[1:])
    bin_width   = float(bin_edges[1] - bin_edges[0])

    # fit_mask: True for bins whose CENTRE lies in the fit region
    fit_mask = bin_centers <= cut_value

    # Only use bins with ≥ 5 counts for chi-squared validity [Ref 6]
    valid_for_fit = fit_mask & (counts >= 5)
    if valid_for_fit.sum() < 4:
        raise RuntimeError(
            f"Only {valid_for_fit.sum()} bins with ≥5 counts in the fit region. "
            "Reduce n_bins or supply more events."
        )

    x_fit = bin_centers[valid_for_fit]
    y_fit = counts[valid_for_fit].astype(float)

    # ── Initial parameter estimates ──────────────────────────────────────────
    # loc0  ~ position of histogram peak within the fit region
    # scale0 ~ rough width estimate
    # norm0  ~ total integral of the histogram × bin_width
    peak_in_fit = np.argmax(counts[fit_mask])
    loc0   = float(bin_centers[fit_mask][peak_in_fit])
    scale0 = float((cut_value - bin_edges[0]) / 8.0)
    norm0  = float(counts.sum()) * bin_width

    # ── Fit ──────────────────────────────────────────────────────────────────
    # sigma = sqrt(observed) is the standard Poisson approximation for
    # Pearson chi-squared minimisation [Ref 6].
    try:
        popt, pcov = optimize.curve_fit(
            _moyal_scaled,
            x_fit, y_fit,
            p0=[loc0, scale0, norm0],
            sigma=np.sqrt(np.maximum(y_fit, 1.0)),
            absolute_sigma=True,
            maxfev=20_000,
        )
    except RuntimeError as exc:
        raise RuntimeError(
            f"Landau fit did not converge: {exc}\n"
            "Try adjusting n_bins or check that the data looks Landau-like."
        ) from exc

    loc, scale, norm = popt
    loc_err, scale_err, norm_err = np.sqrt(np.diag(pcov))
    mpv = loc   # mode of scipy.stats.moyal is at x = loc

    # ── Full fitted curve at every bin centre ────────────────────────────────
    fit_curve = _moyal_scaled(bin_centers, loc, scale, norm)

    # ── Goodness of fit — fit region ─────────────────────────────────────────
    exp_fit = fit_curve[valid_for_fit]
    chi2    = float(np.sum((y_fit - exp_fit) ** 2 / np.maximum(exp_fit, 1e-6)))
    ndf     = int(valid_for_fit.sum()) - 3      # 3 free parameters
    chi2_red = chi2 / max(ndf, 1)
    p_value  = float(stats.chi2.sf(chi2, max(ndf, 1)))

    # ── Goodness of fit — tail region (top 1 - truncation) ───────────────────
    # A large chi2_tail (small p_tail) is physically expected: the Landau
    # tail is enhanced by delta-ray production relative to the fitted model.
    tail_mask = (~fit_mask) & (counts >= 1)
    if tail_mask.sum() > 0:
        obs_tail = counts[tail_mask].astype(float)
        exp_tail = np.maximum(fit_curve[tail_mask], 1e-6)
        chi2_tail  = float(np.sum((obs_tail - exp_tail) ** 2 / exp_tail))
        ndf_tail   = int(tail_mask.sum())
        p_tail     = float(stats.chi2.sf(chi2_tail, ndf_tail))
    else:
        chi2_tail, ndf_tail, p_tail = float("nan"), 0, float("nan")

    return {
        # Fit parameters
        "loc"       : loc,       "loc_err"   : loc_err,
        "scale"     : scale,     "scale_err" : scale_err,
        "norm"      : norm,      "norm_err"  : norm_err,
        "mpv"       : mpv,       "popt"      : popt,
        "pcov"      : pcov,
        # Fit-region goodness-of-fit
        "chi2"      : chi2,      "ndf"       : ndf,
        "chi2_red"  : chi2_red,  "p_value"   : p_value,
        # Tail-region deviation
        "chi2_tail" : chi2_tail, "ndf_tail"  : ndf_tail,
        "p_tail"    : p_tail,
        # Arrays for plotting
        "cut_value"  : cut_value,
        "bin_edges"  : bin_edges,
        "bin_centers": bin_centers,
        "counts"     : counts,
        "fit_curve"  : fit_curve,
        "fit_mask"   : fit_mask,
        "bin_width"  : bin_width,
        "truncation" : truncation,
    }


# ══════════════════════════════════════════════════════════════════════════════
# SECTION 5 — PLOTTING
# ══════════════════════════════════════════════════════════════════════════════

def _fit_band(x: np.ndarray, popt: np.ndarray, pcov: np.ndarray) -> np.ndarray:
    """
    1-sigma uncertainty envelope of the fitted curve via error propagation.

    Uses a central-difference numerical Jacobian:
        sigma_y(x)^2 = J(x) @ pcov @ J(x)^T

    where J_i(x) = d(f)/d(p_i) evaluated at popt.  This is the standard
    approach for propagating parameter uncertainties into a fit curve [Ref 6].
    """
    eps = 1e-6
    jac = np.zeros((len(x), len(popt)))
    for i, p in enumerate(popt):
        step = max(abs(p) * eps, eps)
        p_up = popt.copy(); p_up[i] += step
        p_dn = popt.copy(); p_dn[i] -= step
        jac[:, i] = (_moyal_scaled(x, *p_up) - _moyal_scaled(x, *p_dn)) / (2.0 * step)
    variance = np.einsum("ij,jk,ik->i", jac, pcov, jac)
    return np.sqrt(np.maximum(variance, 0.0))

def _draw_panels(
    ax_top: plt.Axes,
    ax_res: plt.Axes,
    fit_result: Dict,
    x_label: str,
    x_max_display: float,
    log_y: bool,
) -> None:
    """
    Core drawing routine shared by both the linear and log figures.

    Histogram uncertainties are shown as a continuous ±√N step band
    rather than per-bar error bars — this is far more readable at high
    bin density and is the standard approach in modern HEP publications
    (see e.g. the ATLAS / CMS style guides).

    Pull uncertainties are shown as a shaded band for the same reason.

    Parameters
    ----------
    ax_top, ax_res : Axes
        Top (distribution) and bottom (pull) panels.
    fit_result : dict
        Output of fit_landau().
    x_label : str
        Physical quantity label for the x-axis.
    x_max_display : float
        Hard x-axis upper limit — hides the extreme tail visually while
        keeping all data in the fit and statistics.
    log_y : bool
        Use logarithmic y-scale on the top panel.
    """
    bc  = fit_result["bin_centers"]
    cnt = fit_result["counts"].astype(float)
    fc  = fit_result["fit_curve"]
    fm  = fit_result["fit_mask"]
    bw  = fit_result["bin_width"]

    # ── Derived quantities ────────────────────────────────────────────────────
    cnt_err  = np.sqrt(np.maximum(cnt, 1.0))           # Poisson ±√N
    denom    = np.sqrt(np.maximum(fc, 1.0))
    pulls    = (cnt - fc) / denom
    pull_err = cnt_err / denom                          # σ_pull = √N/√F

    # Colours: accessible blue / vermilion pair (distinct even in greyscale)
    C_FIT  = "#2166ac"   # strong blue  — fitted region
    C_TAIL = "#d73027"   # vermilion    — excluded tail

    trunc_pct = int(fit_result["truncation"] * 100)

    # ── Helper: draw stepped error band around histogram ─────────────────────
    # Build the step outline: repeat each centre point by the two bin edges,
    # then fill between (cnt ± err).  This gives a clean band that follows
    # the histogram contour without cluttering each individual bar.
    def _step_band(mask, color, alpha_fill=0.25):
        x_s = np.repeat(bc[mask], 2)
        # shift to bin edges
        half = bw / 2.0
        x_edges = np.empty_like(x_s)
        x_edges[0::2] = bc[mask] - half
        x_edges[1::2] = bc[mask] + half
        y_lo = np.repeat(np.maximum(cnt[mask] - cnt_err[mask], 0.0), 2)
        y_hi = np.repeat(cnt[mask] + cnt_err[mask], 2)
        ax_top.fill_between(x_edges, y_lo, y_hi,
                            step=None, color=color, alpha=alpha_fill,
                            linewidth=0)

    # ── Top panel: fitted region (solid bars + error band) ───────────────────
    ax_top.bar(bc[fm], cnt[fm], width=bw * 0.96,
               color=C_FIT, alpha=0.55,
               label=f"Data — fitted region (bottom {trunc_pct}%)")
    _step_band(fm, C_FIT, alpha_fill=0.35)

    # ── Top panel: excluded tail (hatched bars + error band) ─────────────────
    # Hatching + different colour + explicit label makes the exclusion
    # immediately unambiguous.
    ax_top.bar(bc[~fm], cnt[~fm], width=bw * 0.96,
               color=C_TAIL, alpha=0.35, hatch="////",
               edgecolor=C_TAIL,
               label=f"NOT FITTED — excluded tail (top {100 - trunc_pct}%)")
    _step_band(~fm, C_TAIL, alpha_fill=0.20)

    # ── Fit curve with 1-sigma parameter uncertainty band ────────────────────
    x_dense = np.linspace(float(bc[0]), float(x_max_display), 3000)
    y_dense = _moyal_scaled(x_dense, fit_result["loc"],
                            fit_result["scale"], fit_result["norm"])
    y_sigma = _fit_band(x_dense, fit_result["popt"], fit_result["pcov"])
    ax_top.plot(x_dense, y_dense, color="black", lw=2.0, zorder=5,
                label="Landau (Moyal) fit")
    ax_top.fill_between(x_dense,
                        y_dense - y_sigma, y_dense + y_sigma,
                        color="black", alpha=0.15, zorder=4,
                        label=r"Fit $\pm1\sigma$ (param. uncertainty)")

    # ── Truncation boundary ───────────────────────────────────────────────────
    ax_top.axvline(fit_result["cut_value"], color="dimgray",
                   lw=1.4, ls="--", alpha=0.8,
                   label=f"Truncation at {trunc_pct}th percentile")

    # ── Annotation box ────────────────────────────────────────────────────────
    chi2_tail_red = fit_result["chi2_tail"] / max(fit_result["ndf_tail"], 1)
    ann = (
        rf"MPV $= {fit_result['mpv']:.4g}$ MeV"                                   "\n"
        rf"$\xi = {fit_result['scale']:.4g}"
        rf"\pm{fit_result['scale_err']:.2g}$ MeV"                                 "\n"
        "\n"
        "Fit region:\n"
        rf"$\chi^2/\nu = {fit_result['chi2_red']:.2f}$,"
        rf"  $p = {fit_result['p_value']:.3f}$"                                   "\n"
        "\n"
        "Excluded tail:\n"
        rf"$\chi^2/\nu = {chi2_tail_red:.2e}$,"
        rf"  $p \approx 0$"
    )
    ax_top.text(
        0.975, 0.965, ann,
        transform=ax_top.transAxes, va="top", ha="right", fontsize=9,
        bbox=dict(boxstyle="round,pad=0.5", facecolor="white",
                  edgecolor="#aaaaaa", alpha=0.93),
    )

    y_label = "Events / bin" + (" (log scale)" if log_y else "")
    ax_top.set_ylabel(y_label, fontsize=12)
    ax_top.set_title(
        "TPC Energy-Loss Distribution — HIBEAM Simulation\n"
        "Landau (Moyal) fit, truncated-mean method "
        "[ALICE JINST 5 (2010) P09002]",
        fontsize=11,
    )
    ax_top.legend(fontsize=8.5, loc="upper right",
                  framealpha=0.92, edgecolor="#aaaaaa")
    ax_top.set_xlim(float(bc[0]) - bw, x_max_display)
    if log_y:
        ax_top.set_yscale("log")
        ax_top.set_ylim(bottom=0.5)
    plt.setp(ax_top.get_xticklabels(), visible=False)

    # ── Pull panel: shaded band for ±1σ uncertainty ───────────────────────────
    # Colour each bar by its region, then overlay a grey shaded band for
    # ±pull_err.  This is much cleaner than per-bar error bars at this density.
    ax_res.bar(bc[fm],  pulls[fm],  width=bw * 0.96,
               color=C_FIT, alpha=0.55)
    ax_res.bar(bc[~fm], pulls[~fm], width=bw * 0.96,
               color=C_TAIL, alpha=0.35, hatch="////", edgecolor=C_TAIL)

    # Step-style ±1σ pull uncertainty band
    x_pull_edges = np.concatenate([bc - bw / 2.0, [bc[-1] + bw / 2.0]])
    pull_lo = np.concatenate([pulls - pull_err, [pulls[-1] - pull_err[-1]]])
    pull_hi = np.concatenate([pulls + pull_err, [pulls[-1] + pull_err[-1]]])
    ax_res.fill_between(x_pull_edges, pull_lo, pull_hi,
                        step="post", color="dimgray", alpha=0.30,
                        label=r"$\pm1\sigma$ (stat.)")

    ax_res.axhline(0,  color="black", lw=1.0)
    ax_res.axhline(+2, color="#888888", lw=1.0, ls="--",
                   label=r"$\pm2\sigma$")
    ax_res.axhline(-2, color="#888888", lw=1.0, ls="--")
    ax_res.axvline(fit_result["cut_value"], color="dimgray",
                   lw=1.4, ls="--", alpha=0.8)
    ax_res.set_ylabel(r"$(N-F)\,/\!\sqrt{F}$", fontsize=11)
    ax_res.set_xlabel(x_label, fontsize=12)
    ax_res.set_ylim(-6.5, 6.5)
    ax_res.set_xlim(float(bc[0]) - bw, x_max_display)
    ax_res.legend(fontsize=8, loc="upper right",
                  framealpha=0.92, edgecolor="#aaaaaa")


def plot_dedx(
    fit_result: Dict,
    use_dedx: bool = False,
    output_path: Optional[str] = None,
    x_display_percentile: float = 95.0,
) -> tuple:
    """
    Produce two publication-quality figures of the energy-loss distribution:
    one with a linear y-axis and one with a logarithmic y-axis.  The log
    version reveals structure in the tail that is invisible on a linear scale,
    which is the standard presentation in TPC PID papers [Ref 4, 5].

    The x-axis is clipped at `x_display_percentile` (default 95th percentile
    of the data) so the figure is not dominated by the extreme Landau tail.
    All bins remain in the fit and statistics; only the display range changes.

    Parameters
    ----------
    fit_result : dict
        Output of fit_landau().
    use_dedx : bool
        Label axis as dE/dx [MeV/mm] (True) or sum(Edep) [MeV] (False).
    output_path : str or None
        Base filename.  Two files are written:
          <stem>_linear.<ext>  and  <stem>_log.<ext>
        If None, both figures are displayed interactively.
    x_display_percentile : float
        Hard display cutoff (default 95).  Hides the extreme tail visually.

    Returns
    -------
    (fig_linear, fig_log) : tuple of matplotlib.figure.Figure
    """
    bc = fit_result["bin_centers"]

    # x upper limit for display — use the percentile stored in fit_result
    # (which is over the original values array) via the bin edges as a proxy.
    x_max_display = float(bc[bc <= np.percentile(bc, x_display_percentile)][-1]
                          + fit_result["bin_width"] / 2.0)

    x_label = (
        r"$\mathrm{d}E/\mathrm{d}x\ [\mathrm{MeV}\,\mathrm{mm}^{-1}]$"
        if use_dedx
        else r"$\sum E_\mathrm{dep}\ [\mathrm{MeV}]$"
    )

    figures = {}
    for log_y, suffix in [(False, "linear"), (True, "log")]:
        fig = plt.figure(figsize=(11, 7))
        gs  = gridspec.GridSpec(
            2, 1, height_ratios=[3, 1], hspace=0.04,
            left=0.10, right=0.97, top=0.93, bottom=0.10,
        )
        ax_top = fig.add_subplot(gs[0])
        ax_res = fig.add_subplot(gs[1], sharex=ax_top)
        _draw_panels(ax_top, ax_res, fit_result,
                     x_label, x_max_display, log_y)
        figures[suffix] = fig

    if output_path:
        import os
        stem, ext = os.path.splitext(output_path)
        if not ext:
            ext = ".png"
        for suffix, fig in figures.items():
            path = f"{stem}_{suffix}{ext}"
            fig.savefig(path, dpi=150, bbox_inches="tight")
            print(f"Figure saved: {path}")
    else:
        plt.show()

    return figures["linear"], figures["log"]


# ══════════════════════════════════════════════════════════════════════════════
# SECTION 6 — MAIN PIPELINE
# ══════════════════════════════════════════════════════════════════════════════

def main(
    filepath: str,
    geom: Optional[TPCGeometry] = None,
    truncation: float = 0.70,
    n_bins: int = 100,
    min_steps: int = 5,
    output_path: Optional[str] = None,
) -> None:
    """
    Full analysis pipeline: load → compute energy loss → fit → plot.

    Parameters
    ----------
    filepath : str
        Path to processed ROOT ntuple.
    geom : TPCGeometry or None
        Supply to get true dE/dx [MeV/mm]; None → use sum(Edep) [MeV].
    truncation : float
        Fraction for the truncated-mean fit (default 0.70, per [Ref 4]).
    n_bins : int
        Histogram bins (default 100; reduce if events are sparse).
    min_steps : int
        Minimum G4 steps per event to accept (default 5).
    output_path : str or None
        Save figure to this path; show interactively if None.
    """
    print(f"\nLoading TPC data from:  {filepath}")
    data = load_tpc_data(filepath)
    print(f"  Tree contains {data['n_events']:,} events")

    use_dedx = geom is not None and geom.is_complete
    mode_str = "dE/dx [MeV/mm]" if use_dedx else "sum(Edep) [MeV]  (no geometry supplied)"
    print(f"  Energy-loss mode: {mode_str}")

    print("\nComputing per-event energy loss...")
    values = compute_dedx(data, geom=geom, min_steps=min_steps)
    print(f"  {len(values):,} events accepted (min_steps ≥ {min_steps})")

    if len(values) < 200:
        warnings.warn(
            f"Only {len(values)} events passed the cut — fit may be unreliable.",
            stacklevel=1,
        )

    print(f"\nFitting Landau (Moyal) to bottom {truncation*100:.0f}% of values...")
    result = fit_landau(values, truncation=truncation, n_bins=n_bins)

    # ── Print summary ────────────────────────────────────────────────────────
    chi2_tail_red = result["chi2_tail"] / max(result["ndf_tail"], 1)
    print("\n" + "─" * 52)
    print("  FIT RESULTS")
    print("─" * 52)
    print(f"  Most probable value (MPV)  : {result['mpv']:.4g}")
    print(f"  Width  ξ (scale)           : {result['scale']:.4g}"
          f"  ±  {result['scale_err']:.3g}")
    print(f"  Fit region   chi²/ndf      : {result['chi2_red']:.3f}"
          f"    (p = {result['p_value']:.3f})")
    print(f"  Tail region  chi²/ndf      : {chi2_tail_red:.3f}"
          f"    (p = {result['p_tail']:.3f})")
    print("  (large chi²_tail is expected: delta-ray enhancement [Ref 3, 4])")
    print("─" * 52 + "\n")

    fig_lin, fig_log = plot_dedx(result, use_dedx=use_dedx, output_path=output_path)
    plt.close(fig_lin)
    plt.close(fig_log)


# ══════════════════════════════════════════════════════════════════════════════
# ENTRY POINT
# ══════════════════════════════════════════════════════════════════════════════

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python tpc_dedx_analysis.py <file.root> [output.png]")
        sys.exit(1)

    root_file  = sys.argv[1]
    fig_output = sys.argv[2] if len(sys.argv) > 2 else None

    # ── To enable true dE/dx [MeV/mm], fill in your TPCHit.hh constants: ────
    #
    # geometry = TPCGeometry(
    #     pad_width  = 3.0,    # mm
    #     pad_height = 6.0,    # mm
    #     drift_v    = 0.05,   # mm/ns  (e.g. 5 cm/μs for Ar/CO2)
    #     time_bin   = 100.0,  # ns
    # )
    # main(root_file, geom=geometry, output_path=fig_output)
    #
    # ─────────────────────────────────────────────────────────────────────────

    main(root_file, geom=None, output_path=fig_output)