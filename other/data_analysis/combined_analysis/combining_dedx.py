import ROOT
import os
import uproot
import awkward as ak
import numpy as np
import matplotlib
matplotlib.use('Agg')  # non-interactive backend, no display needed
import matplotlib.pyplot as plt
import mplhep as hep
import matplotlib.image as mpimg
from matplotlib.offsetbox import OffsetImage, AnnotationBbox
from scipy.optimize import curve_fit
from scipy.stats import moyal

# ==============================================================================
# CONFIGURATION
# ==============================================================================
DATA_FILE     = "../experimental_data/tracks_centroids_tpc_run_0042-sorted_t0.root"
SIM_FILE      = "../simulation_data/250proton.root"
CHI2_NDF_MAX  = 25.0   # track quality cut, shared with standalone script
MIN_TRACK_PTS = 7
SIM_PERCENTILE_CUT = 70  # Keep bottom N% of sim energy deposits
N_BINS        = 50
OUTPUT_PATH   = "combined_analysis.png"
FIT_TRUNCATION_PCT = 70
N_BINS             = 50
OUTPUT_PATH        = "combined_analysis.png"


# ==============================================================================
# STEP 1: LOAD REAL DATA  (PyROOT path, same logic as your original script)
# ==============================================================================
def load_data(file_path, min_adc_cut=0.0):
    """
    Reads dE/dx from the tracks branch via SetBranchAddress.
    Avoids centroids branch entirely (CheckByteCount corruption).
    Include path set before .so load — matches working diagnostic order.
    """
    script_dir  = os.path.dirname(os.path.abspath(__file__))
    headers_dir = os.path.join(script_dir, "recovered_headers")
    headers_so  = os.path.join(headers_dir, "recovered_headers.so")
    hdr_file    = os.path.join(headers_dir, "TrackData.h")

    if not os.path.exists(headers_so):
        raise FileNotFoundError(f"recovered_headers.so not found at {headers_so}")

    # ── Critical ordering: include path → header → .so ────────────────────────
    ROOT.gInterpreter.AddIncludePath(headers_dir)
    if os.path.exists(hdr_file):
        ROOT.gInterpreter.ProcessLine(f'#include "{hdr_file}"')
    else:
        raise FileNotFoundError(
            f"TrackData.h not found in {headers_dir}. "
            "Re-run TFile::MakeProject to regenerate headers.")
    ROOT.gSystem.Load(headers_so)

    tfile = ROOT.TFile.Open(file_path)
    if not tfile or tfile.IsZombie():
        raise IOError(f"Cannot open {file_path}")
    tree = tfile.Get("trackingData")
    if not tree:
        tfile.Close()
        raise KeyError("'trackingData' tree not found")

    tracks_vec = ROOT.std.vector("TrackData")()
    tree.SetBranchAddress("tracks", tracks_vec)

    CHI2_NDF_MAX  = 25.0
    MIN_TRACK_PTS = 3
    all_dedx = []

    for i in range(tree.GetEntries()):
        tree.GetEntry(i)
        for track in tracks_vec:
            n = track.nPoints
            if n < MIN_TRACK_PTS:
                continue

            # Chi²/ndf quality cut
            chi2 = 0.0
            for j in range(n):
                y_j = float(track.y[j])
                sx  = max(float(track.sigmas_x[j]), 1e-4)
                sz  = max(float(track.sigmas_z[j]), 1e-4)
                chi2 += ((float(track.x[j]) - (track.slope_xy * y_j + track.intercept_xy)) / sx) ** 2
                chi2 += ((float(track.z[j]) - (track.slope_zy * y_j + track.intercept_zy)) / sz) ** 2
            ndf = max(2 * n - 4, 1)
            if chi2 / ndf > CHI2_NDF_MAX:
                continue

            # 3-D track length → step size per centroid
            x0, x1 = float(track.x[0]), float(track.x[n-1])
            y0, y1 = float(track.y[0]), float(track.y[n-1])
            z0, z1 = float(track.z[0]), float(track.z[n-1])
            length  = np.sqrt((x1-x0)**2 + (y1-y0)**2 + (z1-z0)**2)
            if length < 1e-3:
                continue
            delta_x = length / n

            for j in range(n):
                adc = float(track.charge[j])
                if adc <= min_adc_cut:
                    continue
                all_dedx.append(adc / delta_x)

    tfile.Close()

    arr = np.array(all_dedx, dtype=float)
    print(f"  Data: {len(arr)} dE/dx values after chi²/ndf < {CHI2_NDF_MAX} cut")
    print(f"  Unique values: {len(np.unique(arr))}  (should be >> 10)")

    if len(np.unique(arr)) < 20:
        raise RuntimeError(
            "Too few unique dE/dx values — charge fields are likely still "
            "corrupted. Check that TrackData.h matches the file being read.")

    upper = np.percentile(arr, 95)
    return arr[arr < upper]

# ==============================================================================
# LOAD SIMULATION  -- return FULL array + a separate fit-range mask
# ==============================================================================
def load_sim(file_path, fit_truncation_pct):
    """
    Returns:
        energy_full  -- all events with hits (for display)
        weights_full -- corresponding weights
        fit_cut_val  -- the energy value at fit_truncation_pct (for fit masking)
    """
    tree = uproot.open(f"{file_path}:hibeam")

    edep_jagged    = tree["ProtoTPC/Edep"].array()
    energy         = ak.to_numpy(ak.sum(edep_jagged, axis=1))
    weights_jagged = tree["PrimaryWeight"].array()
    weights        = ak.to_numpy(ak.fill_none(ak.firsts(weights_jagged), 1.0))

    # Remove zero-deposit events (particle missed detector entirely)
    mask    = energy > 0
    energy  = energy[mask]
    weights = weights[mask]

    # Compute the truncation value but do NOT apply it here
    fit_cut_val = np.percentile(energy, fit_truncation_pct)

    print(f"  Sim: {len(energy)} events with hits")
    print(f"  Fit truncation at {fit_truncation_pct}th percentile = {fit_cut_val:.6f} GeV")
    print(f"  ({100-fit_truncation_pct}% of tail excluded from fit, included in plot)")

    return energy, weights, fit_cut_val

# ==============================================================================
# STEP 3: SHARED LANDAU FIT FUNCTION
# ==============================================================================
def landau(x, amp, loc, scale):
    return amp * moyal.pdf(x, loc=loc, scale=scale)

def fit_landau(bin_centers, counts, errors, x_range):
    p0     = [max(counts), bin_centers[np.argmax(counts)], (x_range[1] - x_range[0]) * 0.05]
    bounds = ([0, x_range[0], 1e-9], [np.inf, x_range[1], x_range[1]])
    try:
        popt, _ = curve_fit(landau, bin_centers, counts,
                            p0=p0, sigma=errors, bounds=bounds, maxfev=10000)
        return popt
    except RuntimeError:
        print("  Warning: fit did not converge, using initial guess.")
        return p0

# ==============================================================================
# PSEUDO-CALIBRATION
# ==============================================================================
def compute_calibration(data_arr, sim_energy, sim_w, fit_cut_val):
    """MPV-match data to simulation using only the fit region for both."""

    # Data: fit below 95th percentile (remove extreme ADC noise tail)
    d_fit    = data_arr[data_arr < np.percentile(data_arr, 95)]
    d_params = moyal.fit(d_fit)
    data_mpv = d_params[0]

    # Sim: fit below truncation cut only
    mask    = sim_energy <= fit_cut_val
    s_fit   = sim_energy[mask]
    s_w_fit = sim_w[mask]

    counts, edges = np.histogram(s_fit, bins=N_BINS, weights=s_w_fit)
    centers       = (edges[:-1] + edges[1:]) / 2
    sumw2, _      = np.histogram(s_fit, bins=edges, weights=s_w_fit**2)
    errors        = np.sqrt(sumw2)
    errors[errors == 0] = 1e-9

    p0     = [max(counts), centers[np.argmax(counts)], (edges[-1]-edges[0])*0.05]
    bounds = ([0, edges[0], 1e-9], [np.inf, edges[-1], edges[-1]])
    try:
        popt, _ = curve_fit(
            lambda x, a, l, s: a * moyal.pdf(x, l, s),
            centers, counts, p0=p0, sigma=errors, bounds=bounds, maxfev=10000
        )
        sim_mpv = popt[1]
    except RuntimeError:
        sim_mpv = centers[np.argmax(counts)]
        print("  Warning: sim fit did not converge, using histogram peak.")

    k = sim_mpv / data_mpv
    print(f"\n  Pseudo-calibration:")
    print(f"    Data MPV = {data_mpv:.4f} ADC/Z")
    print(f"    Sim  MPV = {sim_mpv:.6f} GeV")
    print(f"    k        = {k:.4e} GeV/(ADC/Z)\n")
    return k, data_mpv, sim_mpv



# ==============================================================================
# STEP 4: BUILD HISTOGRAMS ON A COMMON AXIS
# ==============================================================================
def make_histogram(values, weights, bins, density=True):
    """
    Returns (counts, edges, errors) on the supplied bin edges.
    density=True normalises so both datasets sit on the same vertical scale.
    """
    counts, edges = np.histogram(values, bins=bins, weights=weights, density=density)
    # Poisson-like error on weighted histogram
    sumw2, _  = np.histogram(values, bins=edges, weights=weights**2)
    bin_width = edges[1] - edges[0]
    norm      = np.sum(counts) * bin_width if density else 1.0
    errors    = np.sqrt(sumw2)
    if density:
        errors /= (np.sum(weights) * bin_width)  # same normalisation as density=True
    errors[errors == 0] = 1e-9  # avoid division-by-zero in curve_fit
    return counts, edges, errors

# ==============================================================================
# STEP 5: PLOT EVERYTHING
# ==============================================================================
# ==============================================================================
# PLOT
# ==============================================================================
def plot_combined(data_file, sim_file):

    print("Loading real data...")
    data_arr = load_data(data_file)

    print("Loading simulation...")
    sim_energy, sim_w, fit_cut_val = load_sim(sim_file, FIT_TRUNCATION_PCT)

    # ── Find each distribution's peak (mode) ──────────────────────────────────
    # Use a coarse histogram just to locate the peak reliably
    d_counts_coarse, d_edges_coarse = np.histogram(data_arr, bins=200)
    d_centers_coarse = (d_edges_coarse[:-1] + d_edges_coarse[1:]) / 2
    data_peak = d_centers_coarse[np.argmax(d_counts_coarse)]

    s_counts_coarse, s_edges_coarse = np.histogram(
        sim_energy[sim_energy <= fit_cut_val], bins=200,
        weights=sim_w[sim_energy <= fit_cut_val]
    )
    s_centers_coarse = (s_edges_coarse[:-1] + s_edges_coarse[1:]) / 2
    sim_peak = s_centers_coarse[np.argmax(s_counts_coarse)]

    print(f"  Data peak (raw):  {data_peak:.4f} ADC/Z")
    print(f"  Sim  peak (raw):  {sim_peak:.6f} GeV")

    # ── Normalise both x-axes to their own peak = 1 ───────────────────────────
    # This makes the x-axis dimensionless: "units of MPV"
    # No assumptions about absolute calibration needed
    data_norm = data_arr   / data_peak
    sim_norm  = sim_energy / sim_peak

    sim_cut_norm = fit_cut_val / sim_peak  # truncation in normalised units

    # ── Shared bin edges in normalised units ──────────────────────────────────
    x_min   = 0
    x_max   = min(
        np.percentile(data_norm, 99.5),
        np.percentile(sim_norm,  99.5)
    )
    edges   = np.linspace(x_min, x_max, N_BINS + 1)
    centers = (edges[:-1] + edges[1:]) / 2

    # ── Histograms normalised to unit area ────────────────────────────────────
    d_counts, _, d_errors = make_histogram(
        data_norm, np.ones(len(data_norm)), edges, density=True)
    s_counts, _, s_errors = make_histogram(
        sim_norm, sim_w, edges, density=True)

    # ── Fits in the bulk only (below truncation) ──────────────────────────────
    fit_mask    = centers <= sim_cut_norm
    x_fit_range = (x_min, sim_cut_norm)

    d_popt = fit_landau(centers[fit_mask], d_counts[fit_mask], d_errors[fit_mask], x_fit_range)
    s_popt = fit_landau(centers[fit_mask], s_counts[fit_mask], s_errors[fit_mask], x_fit_range)

    x_smooth = np.linspace(x_min, x_max, 1000)

    d_fit_bins              = landau(centers, *d_popt)
    s_fit_bins              = landau(centers, *s_popt)
    d_residuals             = np.zeros_like(d_counts)
    s_residuals             = np.zeros_like(s_counts)
    d_residuals[fit_mask]   = d_counts[fit_mask] - d_fit_bins[fit_mask]
    s_residuals[fit_mask]   = s_counts[fit_mask] - s_fit_bins[fit_mask]

    # ── Plot ──────────────────────────────────────────────────────────────────
    plt.style.use(hep.style.CMS)
    fig, (ax_main, ax_resid) = plt.subplots(
        2, 1, sharex=True,
        gridspec_kw={'height_ratios': [3, 1], 'hspace': 0},
        figsize=(11, 9)
    )

    # Delta-ray tail shading
    ax_main.axvspan(sim_cut_norm, x_max, alpha=0.08, color='gray',
                    label=f'Delta-ray tail (excluded from fit, {100-FIT_TRUNCATION_PCT}%)')
    ax_resid.axvspan(sim_cut_norm, x_max, alpha=0.08, color='gray')

    # Data and sim
    hep.histplot(d_counts, bins=edges, yerr=d_errors, ax=ax_main,
                 histtype='errorbar', color='black', label='Data')
    hep.histplot(s_counts, bins=edges, yerr=s_errors, ax=ax_main,
                 histtype='fill', color='steelblue', alpha=0.4,
                 label='Simulation (weighted)')

    # Fit lines
    ax_main.plot(x_smooth, landau(x_smooth, *d_popt), color='#E31A1C', lw=2,
                 label=f'Data Landau fit  MPV = {d_popt[1]:.3f} × peak')
    ax_main.plot(x_smooth, landau(x_smooth, *s_popt), color='darkorange', lw=2, ls='--',
                 label=f'Sim  Landau fit   MPV = {s_popt[1]:.3f} × peak')

    # Truncation line
    ax_main.axvline(sim_cut_norm, color='green', lw=1.5, ls=':',
                    label=f'Fit truncation ({FIT_TRUNCATION_PCT}th pct)')
    ax_resid.axvline(sim_cut_norm, color='green', lw=1.5, ls=':')

    ax_main.set_yscale('log')
    ax_main.set_ylabel("Probability Density (normalised)")
    ax_main.set_ylim(bottom=1e-4)
    ax_main.legend(loc='upper right', fontsize='small')

    # Annotation box showing the raw calibration values for the record
    info = (f"Data peak: {data_peak:.2f} ADC/Z\n"
            f"Sim  peak: {sim_peak:.4f} GeV\n"
            f"Ratio k:   {sim_peak/data_peak:.3e} GeV/(ADC/Z)")
    ax_main.text(0.97, 0.45, info, transform=ax_main.transAxes,
                 fontsize=8, va='top', ha='right',
                 bbox=dict(boxstyle='round', facecolor='white', alpha=0.7))

    main_text, _ = hep.cms.label(loc=0, data=False, label="Preliminary",
                                  rlabel="HIBEAM", ax=ax_main)
    main_text.set_text("ESS")

    try:
        logo = mpimg.imread("data/ess_logo.png")
        ab = AnnotationBbox(OffsetImage(logo, zoom=0.15), (0.03, 0.95),
                            xycoords='axes fraction', box_alignment=(0, 1), frameon=False)
        ax_main.add_artist(ab)
    except FileNotFoundError:
        pass

    # Residuals
    ax_resid.errorbar(centers[fit_mask], d_residuals[fit_mask], yerr=d_errors[fit_mask],
                      fmt='o', color='black', markersize=4, label='Data residuals')
    ax_resid.errorbar(centers[fit_mask], s_residuals[fit_mask], yerr=s_errors[fit_mask],
                      fmt='s', color='steelblue', markersize=4, alpha=0.7, label='Sim residuals')
    ax_resid.axhline(0, color='gray', ls='--', lw=1)

    limit = max(np.max(np.abs(np.concatenate(
        [d_residuals[fit_mask], s_residuals[fit_mask]]))), 1e-4) * 1.5
    ax_resid.set_ylim(-limit, limit)
    ax_resid.set_xlabel("dE/dx  [units of peak MPV]  — peaks aligned by normalisation")
    ax_resid.set_ylabel("Data − Fit")
    ax_resid.legend(loc='upper right', fontsize='small')

    plt.tight_layout()
    plt.savefig(OUTPUT_PATH, dpi=300)
    print(f"Saved → {OUTPUT_PATH}")

if __name__ == "__main__":
    plot_combined(DATA_FILE, SIM_FILE)