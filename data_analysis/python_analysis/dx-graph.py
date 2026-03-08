import ROOT
import os
import glob
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import moyal

# 1. Load the reconstructed C++ headers
if os.path.exists("recovered_headers/recovered_headers.so"):
    ROOT.gSystem.Load("recovered_headers/recovered_headers.so")
else:
    raise FileNotFoundError("recovered_headers.so not found. Please run MakeProject first.")

def extract_dedx(file_path, min_adc_cut):
    """Extracts dE/dx data from a single ROOT file, filtering out electronic noise."""
    f = ROOT.TFile.Open(file_path)
    tree = f.Get("trackingData")
    
    dedx_values = []
    
    for i in range(tree.GetEntries()):
        tree.GetEntry(i)
        for c in tree.centroids:
            # Physics Cut 1: Ignore exact zeros or negative Z values
            # Physics Cut 2: Mitigate electronic background with an ADC threshold
            if c.z > 0 and c.integrated_ADC_amplitude > min_adc_cut:
                dedx = c.integrated_ADC_amplitude / c.z
                dedx_values.append(dedx)
                
    f.Close()
    return np.array(dedx_values)

def analyze_and_plot(file_path, output_dir, min_adc_cut=500):
    """Calculates fit and residuals, then saves the combined plot."""
    filename = os.path.basename(file_path)
    print(f"Processing: {filename}...")
    
    # Extract data
    data = extract_dedx(file_path, min_adc_cut)
    
    if len(data) < 100:
        print(f"  -> Skipping {filename}: Not enough data after cuts ({len(data)} hits).")
        return

    # Mitigate the extreme high-end tail to help the fit converge properly
    upper_limit = np.percentile(data, 95)
    fit_data = data[data < upper_limit]

    # Fit the Moyal (Landau) distribution
    params = moyal.fit(fit_data)
    mpv, width = params[0], params[1]

    # --- Create the Combined Figure ---
    # We use gridspec to make the top plot larger than the residuals plot
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8), 
                                   gridspec_kw={'height_ratios': [3, 1]}, 
                                   sharex=True)
    
    # 1. Top Plot: Histogram and Fit
    counts, bins, _ = ax1.hist(fit_data, bins=100, density=True, 
                               alpha=0.6, color='skyblue', label='Data (Noise Cut Applied)')
    
    bin_centers = (bins[:-1] + bins[1:]) / 2
    fit_pdf = moyal.pdf(bin_centers, *params)
    
    ax1.plot(bin_centers, fit_pdf, 'r-', lw=2, 
             label=f'Landau Fit\nMPV: {mpv:.2f}\nWidth: {width:.2f}')
    
    ax1.set_title(f'dE/dx Energy Loss & Fit: {filename}')
    ax1.set_ylabel('Probability Density')
    ax1.legend()
    ax1.grid(axis='y', alpha=0.3)

    # 2. Bottom Plot: Residuals
    residuals = counts - fit_pdf
    ax2.bar(bin_centers, residuals, width=(bins[1]-bins[0]), color='purple', alpha=0.7)
    ax2.axhline(0, color='black', linestyle='--', lw=1)
    
    ax2.set_xlabel('dE/dx (ADC / Z-distance)')
    ax2.set_ylabel('Data - Fit')
    ax2.grid(axis='y', alpha=0.3)

    # Polish and Save
    plt.tight_layout()
    
    out_name = filename.replace(".root", "_analysis.png")
    out_path = os.path.join(output_dir, out_name)
    plt.savefig(out_path, dpi=300)
    plt.close(fig) # Close the figure to free up memory
    
    print(f"  -> Saved plot to {out_path}")

def run_batch(input_directory, output_directory, min_adc_cut=500):
    """Finds all .root files and processes them."""
    if not os.path.exists(output_directory):
        os.makedirs(output_directory)
        
    search_pattern = os.path.join(input_directory, "*.root")
    root_files = glob.glob(search_pattern)
    
    if not root_files:
        print(f"No .root files found in {input_directory}")
        return
        
    print(f"Found {len(root_files)} files. Starting batch process...\n")
    
    for file in root_files:
        analyze_and_plot(file, output_directory, min_adc_cut)
        
    print("\nBatch processing complete!")

# --- Execution ---
if __name__ == "__main__":
    # Define your folders here
    INPUT_DIR = "../experimental_data/" 
    OUTPUT_DIR = "analysis_plots/"
    
    # Change this threshold depending on where your electronic noise drops off
    # If the spike is still throwing off the fit, increase this number.
    ADC_NOISE_THRESHOLD = 500 
    
    run_batch(INPUT_DIR, OUTPUT_DIR, ADC_NOISE_THRESHOLD)