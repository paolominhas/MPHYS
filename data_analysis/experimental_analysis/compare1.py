"""
Simulation vs Experimental Data Comparison
==========================================

Compares dE/dx distributions from:
  - Simulation: Geant4 output via uproot (energy deposits in GeV)
  - Experiment: TPC tracking data via PyROOT (ADC amplitudes)

IMPORTANT: The simulation records energy in GeV, while the experimental
data records integrated ADC amplitude / z-distance. These are different
units. We normalise both distributions to unit area so we can compare
SHAPES, but the x-axes are in different units. A proper comparison
requires either:
  (a) Converting ADC → energy using a calibration constant, or
  (b) Converting simulation energy → ADC using the inverse calibration.
Without calibration, we compare shapes only and report the MPV ratio
as an effective calibration factor.
"""

import os
import glob
import numpy as np
import matplotlib.pyplot as plt
import mplhep as hep
import ROOT
import uproot
import awkward as ak
from scipy.stats import moyal, ks_2samp
from scipy.optimize import curve_fit

plt.style.use(hep.style.CMS)

# ============================================================
# Configuration — edit these paths and parameters
# ============================================================

SIM_FILE = "../../initial-proj/final_sim_data/250proton.root"
SIM_TREE = "hibeam"
SIM_BRANCH_EDEP = "ProtoTPC/Edep"
SIM_BRANCH_WEIGHT = "PrimaryWeight"

EXP_DATA_DIR = "../experimental_data/"
EXP_HEADERS = "recovered_headers/recovered_headers.so"
ADC_NOISE_THRESHOLD = 500

N_BINS = 60
OUTPUT_DIR = "comparison_plots/"

# ============================================================
# 1. Load experimental data
# ============================================================

def load_headers():
    """Load headers by compiling at runtime instead of using .so"""
    header_dir = os.path.dirname(EXP_HEADERS)
    headers = ["CentroidData.h", "TrackData.h", "tpcData.h"]
    
    for h in headers:
        path = os.path.join(header_dir, h)
        if os.path.exists(path):
            ROOT.gInterpreter.Declare(open(path).read())
            print(f"  Loaded header: {h}")
        else:
            print(f"  WARNING: {h} not found")

def extract_dedx(file_path, min_adc_cut):
    """Extract dE/dx using TTree::Draw to avoid cppyy casting issues."""
    f = ROOT.TFile.Open(file_path)
    if not f or f.IsZombie():
        print(f"  WARNING: Cannot open {file_path}, skipping.")
        return np.array([])

    tree = f.Get("trackingData")
    if not tree:
        print(f"  WARNING: No 'trackingData' tree in {file_path}, skipping.")
        f.Close()
        return np.array([])

    # Use TTree::Draw into a histogram to avoid cppyy vector issues
    n_entries = tree.GetEntries()
    if n_entries == 0:
        f.Close()
        return np.array([])

    # Draw dE/dx with cuts directly — ROOT handles the vector iteration
    cut = f"centroids.z>0 && centroids.integrated_ADC_amplitude>{min_adc_cut}"
    expr = "centroids.integrated_ADC_amplitude/centroids.z"
    
    n_drawn = tree.Draw(f"{expr}>>htmp(10000, 0, 50000)", cut, "goff")
    
    if n_drawn <= 0:
        f.Close()
        return np.array([])

    # Extract drawn values from the tree's internal buffer
    buf = tree.GetV1()
    dedx_values = np.array([buf[i] for i in range(n_drawn)])

    f.Close()
    return dedx_values

def load_all_experimental(data_dir, min_adc_cut):
    """Loop over all .root files in a directory and extract dE/dx."""
    search = os.path.join(data_dir, "*.root")
    files = sorted(glob.glob(search))

    if not files:
        raise FileNotFoundError(f"No .root files found in {data_dir}")

    print(f"Found {len(files)} experimental files:")
    all_dedx = []
    for fpath in files:
        fname = os.path.basename(fpath)
        dedx = extract_dedx(fpath, min_adc_cut)
        print(f"  {fname}: {len(dedx)} hits")
        if len(dedx) > 0:
            all_dedx.append(dedx)

    if not all_dedx:
        raise ValueError("No valid dE/dx data extracted from any file.")

    combined = np.concatenate(all_dedx)
    print(f"Total experimental dE/dx entries: {len(combined)}")
    return combined

# ============================================================
# 2. Load simulation data
# ============================================================

def load_simulation(sim_file, tree_name, edep_branch, weight_branch):
    """
    Load simulation energy deposits from a ROOT file via uproot.
    
    Sums per-hit energy deposits into a total per event,
    and extracts the primary particle weights.
    """
    tree = uproot.open(f"{sim_file}:{tree_name}")

    edep_jagged = tree[edep_branch].array()
    energy = ak.to_numpy(ak.sum(edep_jagged, axis=1))

    weights = ak.to_numpy(
        ak.fill_none(ak.firsts(tree[weight_branch].array()), 1.0)
    )

    # Keep only events with actual energy deposits
    mask = energy > 0
    energy = energy[mask]
    weights = weights[mask]

    print(f"\nSimulation loaded from: {sim_file}")
    print(f"  Total events: {len(mask)}")
    print(f"  Events with deposits: {len(energy)}")
    print(f"  Energy range: {energy.min():.6f} - {energy.max():.6f} GeV")

    return energy, weights

# ============================================================
# 3. Fitting
# ============================================================

def landau_approx(x, amp, loc, scale):
    """Moyal (Landau approximation) scaled by amplitude."""
    return amp * moyal.pdf(x, loc=loc, scale=scale)

def fit_landau(bin_centres, counts, errors, x_max):
    """
    Fit a Moyal distribution to binned histogram data.
    Returns (parameters, uncertainties).
    """
    mask = counts > 0
    if np.sum(mask) < 5:
        print("  WARNING: Too few non-zero bins for fit.")
        return None, None

    p0 = [np.max(counts), bin_centres[np.argmax(counts)],
          (bin_centres[-1] - bin_centres[0]) * 0.05]

    errors_safe = errors.copy()
    errors_safe[errors_safe == 0] = np.max(errors_safe) * 0.1

    try:
        popt, pcov = curve_fit(
            landau_approx, bin_centres[mask], counts[mask],
            p0=p0, sigma=errors_safe[mask],
            bounds=([0, 0, 1e-8], [np.inf, x_max * 2, x_max]),
            maxfev=10000
        )
        perr = np.sqrt(np.diag(pcov))
        return popt, perr
    except RuntimeError:
        print("  WARNING: Fit did not converge.")
        return None, None

# ============================================================
# 4. Plotting
# ============================================================

def plot_single_distribution(data, weights, bin_edges, label, colour,
                             fit_params, ax):
    """Plot a single normalised histogram with its Landau fit."""
    counts, _ = np.histogram(data, bins=bin_edges, weights=weights)
    sumw2, _ = np.histogram(data, bins=bin_edges, weights=weights**2)
    errors = np.sqrt(sumw2)

    norm = np.sum(counts * np.diff(bin_edges))
    if norm == 0:
        return counts, errors, counts, errors

    counts_n = counts / norm
    errors_n = errors / norm

    hep.histplot(counts_n, bins=bin_edges, yerr=errors_n,
                 ax=ax, histtype="errorbar" if "Exp" in label else "step",
                 color=colour, label=label)

    if fit_params is not None:
        x_smooth = np.linspace(bin_edges[0], bin_edges[-1], 500)
        ax.plot(x_smooth, landau_approx(x_smooth, *fit_params),
                color=colour, linestyle="--", lw=1.5, alpha=0.6,
                label=f"Fit MPV={fit_params[1]:.4f}")

    return counts, errors, counts_n, errors_n

def plot_comparison(energy_exp, energy_sim, weights_sim, output_path):
    """
    Create the full comparison figure: overlay + ratio + residuals.
    """
    # --- Common binning on each dataset's own scale ---
    # Since units differ, we normalise each to unit area independently
    e_min_exp = 0
    e_max_exp = np.percentile(energy_exp, 99)
    bins_exp = np.linspace(e_min_exp, e_max_exp, N_BINS + 1)
    centres_exp = (bins_exp[:-1] + bins_exp[1:]) / 2

    e_min_sim = 0
    e_max_sim = np.percentile(energy_sim, 99)
    bins_sim = np.linspace(e_min_sim, e_max_sim, N_BINS + 1)
    centres_sim = (bins_sim[:-1] + bins_sim[1:]) / 2

    # Histogram and normalise
    counts_exp, _ = np.histogram(energy_exp, bins=bins_exp)
    err_exp = np.sqrt(counts_exp)
    norm_exp = np.sum(counts_exp * np.diff(bins_exp))
    counts_exp_n = counts_exp / norm_exp
    err_exp_n = err_exp / norm_exp

    counts_sim, _ = np.histogram(energy_sim, bins=bins_sim, weights=weights_sim)
    sumw2_sim, _ = np.histogram(energy_sim, bins=bins_sim, weights=weights_sim**2)
    err_sim = np.sqrt(sumw2_sim)
    norm_sim = np.sum(counts_sim * np.diff(bins_sim))
    counts_sim_n = counts_sim / norm_sim
    err_sim_n = err_sim / norm_sim

    # Fit both
    print("\nFitting experimental data...")
    popt_exp, perr_exp = fit_landau(centres_exp, counts_exp_n,
                                     err_exp_n, e_max_exp)

    print("Fitting simulation...")
    popt_sim, perr_sim = fit_landau(centres_sim, counts_sim_n,
                                     err_sim_n, e_max_sim)

    if popt_exp is not None and popt_sim is not None:
        print(f"\n{'='*55}")
        print(f"  Experimental: MPV = {popt_exp[1]:.4f} ± {perr_exp[1]:.4f}")
        print(f"  Simulation:   MPV = {popt_sim[1]:.4f} ± {perr_sim[1]:.4f}")
        print(f"  MPV ratio (exp/sim) = {popt_exp[1]/popt_sim[1]:.2f}")
        print(f"  (This ratio is an effective ADC-to-GeV calibration factor)")
        print(f"{'='*55}")

    # --- Figure: two panels side by side + residuals ---
    fig = plt.figure(figsize=(14, 10))

    # Top left: experimental
    ax_exp = fig.add_subplot(2, 2, 1)
    hep.histplot(counts_exp_n, bins=bins_exp, yerr=err_exp_n,
                 ax=ax_exp, histtype="errorbar", color="black",
                 label="Experimental data")
    if popt_exp is not None:
        x_sm = np.linspace(e_min_exp, e_max_exp, 500)
        ax_exp.plot(x_sm, landau_approx(x_sm, *popt_exp),
                    "b-", lw=2, label=f"Fit (MPV={popt_exp[1]:.2f})")
    ax_exp.set_xlabel("dE/dx (ADC / z)")
    ax_exp.set_ylabel("Normalised density")
    ax_exp.set_title("Experimental dE/dx")
    ax_exp.legend()

    # Top right: simulation
    ax_sim = fig.add_subplot(2, 2, 2)
    hep.histplot(counts_sim_n, bins=bins_sim, yerr=err_sim_n,
                 ax=ax_sim, histtype="step", color="red",
                 label="Simulation (weighted)")
    if popt_sim is not None:
        x_sm = np.linspace(e_min_sim, e_max_sim, 500)
        ax_sim.plot(x_sm, landau_approx(x_sm, *popt_sim),
                    "r--", lw=2, label=f"Fit (MPV={popt_sim[1]:.4f} GeV)")
    ax_sim.set_xlabel("Energy deposited (GeV)")
    ax_sim.set_ylabel("Normalised density")
    ax_sim.set_title("Simulated energy deposit")
    ax_sim.legend()

    # Bottom left: overlay using calibrated x-axis
    ax_overlay = fig.add_subplot(2, 2, 3)

    if popt_exp is not None and popt_sim is not None:
        # Scale experimental x-axis so MPVs align
        # This effectively converts ADC/z → GeV using MPV as calibration
        calib_factor = popt_sim[1] / popt_exp[1]

        energy_exp_calib = energy_exp * calib_factor
        e_max_calib = np.percentile(energy_exp_calib, 99)
        common_max = max(e_max_calib, e_max_sim)
        common_bins = np.linspace(0, common_max, N_BINS + 1)
        common_centres = (common_bins[:-1] + common_bins[1:]) / 2

        c_exp_c, _ = np.histogram(energy_exp_calib, bins=common_bins)
        e_exp_c = np.sqrt(c_exp_c)
        n_exp_c = np.sum(c_exp_c * np.diff(common_bins))
        c_exp_cn = c_exp_c / n_exp_c if n_exp_c > 0 else c_exp_c
        e_exp_cn = e_exp_c / n_exp_c if n_exp_c > 0 else e_exp_c

        c_sim_c, _ = np.histogram(energy_sim, bins=common_bins,
                                   weights=weights_sim)
        sw2_sim_c, _ = np.histogram(energy_sim, bins=common_bins,
                                     weights=weights_sim**2)
        e_sim_c = np.sqrt(sw2_sim_c)
        n_sim_c = np.sum(c_sim_c * np.diff(common_bins))
        c_sim_cn = c_sim_c / n_sim_c if n_sim_c > 0 else c_sim_c
        e_sim_cn = e_sim_c / n_sim_c if n_sim_c > 0 else e_sim_c

        hep.histplot(c_exp_cn, bins=common_bins, yerr=e_exp_cn,
                     ax=ax_overlay, histtype="errorbar", color="black",
                     label="Data (MPV-calibrated)")
        hep.histplot(c_sim_cn, bins=common_bins, yerr=e_sim_cn,
                     ax=ax_overlay, histtype="step", color="red",
                     label="Simulation")

        ax_overlay.set_xlabel("Energy deposited (GeV, calibrated)")
        ax_overlay.set_ylabel("Normalised density")
        ax_overlay.set_title("Overlay (data scaled by MPV ratio)")
        ax_overlay.legend()
    else:
        ax_overlay.text(0.5, 0.5, "Fit failed — cannot overlay",
                        transform=ax_overlay.transAxes, ha='center')

    # Bottom right: ratio
    ax_ratio = fig.add_subplot(2, 2, 4)

    if popt_exp is not None and popt_sim is not None:
        with np.errstate(divide="ignore", invalid="ignore"):
            ratio = np.where(c_sim_cn > 0, c_exp_cn / c_sim_cn, np.nan)
            ratio_err = np.where(c_sim_cn > 0, e_exp_cn / c_sim_cn, np.nan)

        ax_ratio.errorbar(common_centres, ratio, yerr=ratio_err,
                          fmt="o", color="black", markersize=3)
        ax_ratio.axhline(1.0, color="red", linestyle="--", lw=1)
        ax_ratio.set_xlabel("Energy deposited (GeV, calibrated)")
        ax_ratio.set_ylabel("Data / Simulation")
        ax_ratio.set_ylim(0.0, 2.5)
        ax_ratio.set_title("Ratio (calibrated data / simulation)")
        ax_ratio.grid(axis="y", alpha=0.3)
    else:
        ax_ratio.text(0.5, 0.5, "Fit failed — cannot compute ratio",
                      transform=ax_ratio.transAxes, ha='center')

    plt.tight_layout()
    plt.savefig(output_path, dpi=300)
    plt.show()
    print(f"\nPlot saved to: {output_path}")

    # --- Return fit results for metrics ---
    return (popt_exp, perr_exp, popt_sim, perr_sim,
            c_exp_cn if popt_exp is not None else None,
            c_sim_cn if popt_sim is not None else None,
            e_exp_cn if popt_exp is not None else None,
            e_sim_cn if popt_sim is not None else None,
            energy_exp_calib if popt_exp is not None else None)

# ============================================================
# 5. Quantitative agreement metrics
# ============================================================

def print_metrics(popt_exp, perr_exp, popt_sim, perr_sim,
                  c_exp_n, c_sim_n, e_exp_n, e_sim_n,
                  energy_exp_calib, energy_sim):
    """Print chi2, KS test, and MPV comparison."""
    print(f"\n{'='*55}")
    print("  QUANTITATIVE AGREEMENT METRICS")
    print(f"{'='*55}")

    # MPV comparison
    if popt_exp is not None and popt_sim is not None:
        mpv_diff = abs(popt_exp[1] - popt_sim[1])
        mpv_err = np.sqrt(perr_exp[1]**2 + perr_sim[1]**2)
        print(f"  Exp MPV:  {popt_exp[1]:.4f} ± {perr_exp[1]:.4f}")
        print(f"  Sim MPV:  {popt_sim[1]:.4f} ± {perr_sim[1]:.4f}")
        print(f"  Exp width: {popt_exp[2]:.4f} ± {perr_exp[2]:.4f}")
        print(f"  Sim width: {popt_sim[2]:.4f} ± {perr_sim[2]:.4f}")

    # Chi-squared on calibrated overlay
    if c_exp_n is not None and c_sim_n is not None:
        mask = (c_sim_n > 0) & (c_exp_n > 0)
        combined_err2 = e_sim_n**2 + e_exp_n**2
        combined_err2[combined_err2 == 0] = 1

        chi2 = np.sum(
            (c_exp_n[mask] - c_sim_n[mask])**2 / combined_err2[mask]
        )
        ndf = np.sum(mask) - 3  # subtract 3 fit parameters
        ndf = max(ndf, 1)
        print(f"\n  Chi2 / NDF = {chi2:.1f} / {ndf} = {chi2/ndf:.2f}")
        if chi2 / ndf < 2:
            print("  → Good agreement (chi2/NDF < 2)")
        elif chi2 / ndf < 5:
            print("  → Moderate agreement")
        else:
            print("  → Poor agreement — shapes differ significantly")

    # KS test on calibrated unbinned data
    if energy_exp_calib is not None:
        ks_stat, ks_pval = ks_2samp(energy_exp_calib, energy_sim)
        print(f"\n  KS statistic: {ks_stat:.4f}")
        print(f"  KS p-value:   {ks_pval:.4e}")
        if ks_pval > 0.05:
            print("  → Cannot reject null hypothesis (distributions compatible)")
        else:
            print("  → Distributions differ significantly (p < 0.05)")

    print(f"{'='*55}")

# ============================================================
# 6. Main execution
# ============================================================

if __name__ == "__main__":
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    # Load C++ headers for experimental data
    load_headers()

    # Load experimental data from all ROOT files
    print("Loading experimental data...")
    energy_exp = load_all_experimental(EXP_DATA_DIR, ADC_NOISE_THRESHOLD)

    # Load simulation
    print("\nLoading simulation...")
    energy_sim, weights_sim = load_simulation(
        SIM_FILE, SIM_TREE, SIM_BRANCH_EDEP, SIM_BRANCH_WEIGHT
    )

    # Plot and compare
    results = plot_comparison(
        energy_exp, energy_sim, weights_sim,
        os.path.join(OUTPUT_DIR, "sim_vs_data_comparison.png")
    )

    # Print metrics
    if results is not None:
        (popt_exp, perr_exp, popt_sim, perr_sim,
         c_exp_n, c_sim_n, e_exp_n, e_sim_n,
         energy_exp_calib) = results

        print_metrics(popt_exp, perr_exp, popt_sim, perr_sim,
                      c_exp_n, c_sim_n, e_exp_n, e_sim_n,
                      energy_exp_calib, energy_sim)