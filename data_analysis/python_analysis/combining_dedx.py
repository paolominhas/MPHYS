import ROOT
import os
import uproot
import awkward as ak
import numpy as np
import matplotlib.pyplot as plt
import mplhep as hep
import matplotlib.image as mpimg
from matplotlib.offsetbox import OffsetImage, AnnotationBbox
from scipy.optimize import curve_fit
from scipy.stats import moyal

# ==============================================================================
# CONFIGURATION
# ==============================================================================
DATA_FILE     = "../experimental_data/tracks_centroids_tpc_run_0006-sorted_t0_100.root"
SIM_FILE      = "../../initial-proj/final_sim_data/250proton.root"
MIN_ADC_CUT   = 500      # Noise floor for real data
SIM_PERCENTILE_CUT = 70  # Keep bottom N% of sim energy deposits
N_BINS        = 50
OUTPUT_PATH   = "combined_analysis.png"

# ==============================================================================
# STEP 1: LOAD REAL DATA  (PyROOT path, same logic as your original script)
# ==============================================================================
def load_data(file_path, min_adc_cut):
    """Returns a flat numpy array of dE/dx values from real detector data."""
    if os.path.exists("recovered_headers/recovered_headers.so"):
        ROOT.gSystem.Load("recovered_headers/recovered_headers.so")
    else:
        raise FileNotFoundError("recovered_headers.so not found.")

    f = ROOT.TFile.Open(file_path)
    tree = f.Get("trackingData")

    dedx_values = []
    for i in range(tree.GetEntries()):
        tree.GetEntry(i)
        for c in tree.centroids:
            if c.z > 0 and c.integrated_ADC_amplitude > min_adc_cut:
                dedx_values.append(c.integrated_ADC_amplitude / c.z)

    f.Close()
    arr = np.array(dedx_values)

    # Clip top 5% tail so the fit axis isn't dominated by outliers
    upper = np.percentile(arr, 95)
    return arr[arr < upper]

# ==============================================================================
# STEP 2: LOAD SIMULATION  (uproot path, same logic as your original script)
# ==============================================================================
def load_sim(file_path, percentile_cut):
    """Returns flat (energy, weights) numpy arrays from simulation."""
    tree = uproot.open(f"{file_path}:hibeam")

    edep_jagged   = tree["ProtoTPC/Edep"].array()
    energy        = ak.to_numpy(ak.sum(edep_jagged, axis=1))

    weights_jagged = tree["PrimaryWeight"].array()
    weights        = ak.to_numpy(ak.fill_none(ak.firsts(weights_jagged), 1.0))

    # Remove zero-deposit events
    mask    = energy > 0
    energy  = energy[mask]
    weights = weights[mask]

    # Apply the top-end percentile cut
    cut_val = np.percentile(energy, percentile_cut)
    mask2   = energy <= cut_val
    return energy[mask2], weights[mask2], cut_val

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
def plot_combined(data_file, sim_file):

    # --- Load ---
    print("Loading real data...")
    data_arr          = load_data(data_file, MIN_ADC_CUT)

    print("Loading simulation...")
    sim_energy, sim_w, cut_val = load_sim(sim_file, SIM_PERCENTILE_CUT)

    # --- Shared bin edges: use data range as the reference ---
    x_min  = min(data_arr.min(), sim_energy.min())
    x_max  = max(np.percentile(data_arr, 99), cut_val)
    edges  = np.linspace(x_min, x_max, N_BINS + 1)
    centers = (edges[:-1] + edges[1:]) / 2

    # --- Histograms (density-normalised so both live on the same y-axis) ---
    d_counts, _, d_errors = make_histogram(data_arr,   np.ones(len(data_arr)), edges, density=True)
    s_counts, _, s_errors = make_histogram(sim_energy, sim_w,                  edges, density=True)

    # --- Fits ---
    x_range = (x_min, x_max)
    d_popt  = fit_landau(centers, d_counts, d_errors, x_range)
    s_popt  = fit_landau(centers, s_counts, s_errors, x_range)

    x_smooth   = np.linspace(x_min, x_max, 1000)
    d_fit_bins = landau(centers, *d_popt)
    s_fit_bins = landau(centers, *s_popt)
    d_residuals = d_counts - d_fit_bins
    s_residuals = s_counts - s_fit_bins

    # --- Figure ---
    plt.style.use(hep.style.CMS)
    fig, (ax_main, ax_resid) = plt.subplots(
        2, 1,
        sharex=True,
        gridspec_kw={'height_ratios': [3, 1], 'hspace': 0},
        figsize=(11, 9)
    )

    # ── Main panel ─────────────────────────────────────────────────────────────

    # real data
    hep.histplot(d_counts, bins=edges, yerr=d_errors, ax=ax_main,
                 histtype='errorbar', color='black', label='Data (noise cut applied)')

    # simulation
    hep.histplot(s_counts, bins=edges, yerr=s_errors, ax=ax_main,
                 histtype='fill', color='steelblue', alpha=0.4, label='Simulation (weighted)')

    # fit
    ax_main.plot(x_smooth, landau(x_smooth, *d_popt), color='#E31A1C', lw=2,
                 label=f'Data Landau fit  MPV={d_popt[1]:.3f}')
    ax_main.plot(x_smooth, landau(x_smooth, *s_popt), color='darkorange', lw=2, linestyle='--',
                 label=f'Sim  Landau fit  MPV={s_popt[1]:.3f}')

    ax_main.set_yscale('log')
    ax_main.set_ylabel("Probability Density")
    ax_main.set_ylim(bottom=1e-4)
    ax_main.legend(loc='upper right')

    # ess brand
    main_text, _ = hep.cms.label(loc=0, data=False, label="Preliminary",
                                 rlabel="HIBEAM", ax=ax_main)
    main_text.set_text("ESS")

    try:
        logo = mpimg.imread("data/ess_logo.png")
        ab = AnnotationBbox(OffsetImage(logo, zoom=0.15), (0.03, 0.95),
                            xycoords='axes fraction', box_alignment=(0, 1), frameon=False)
        ax_main.add_artist(ab)
    except FileNotFoundError:
        print("ESS logo not found, skipping.")

    # residuals
    ax_resid.errorbar(centers, d_residuals, yerr=d_errors,
                      fmt='o', color='black', markersize=4, label='Data residuals')
    ax_resid.errorbar(centers, s_residuals, yerr=s_errors,
                      fmt='s', color='steelblue', markersize=4, alpha=0.7, label='Sim residuals')
    ax_resid.axhline(0, color='gray', linestyle='--', lw=1)

    limit = max(np.max(np.abs(np.concatenate([d_residuals, s_residuals]))), 1e-4) * 1.5
    ax_resid.set_ylim(-limit, limit)
    ax_resid.set_xlabel("dE/dx  [ADC / Z-distance]  |  Energy Deposited [GeV]")
    ax_resid.set_ylabel("Data − Fit")
    ax_resid.legend(loc='upper right', fontsize='small')

    # save
    plt.tight_layout()
    plt.savefig(OUTPUT_PATH, dpi=300)
    print(f"\nSaved → {OUTPUT_PATH}")
    plt.show()

if __name__ == "__main__":
    plot_combined(DATA_FILE, SIM_FILE)