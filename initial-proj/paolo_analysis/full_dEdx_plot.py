import uproot
import awkward as ak
import matplotlib.pyplot as plt
import mplhep as hep
import numpy as np
from scipy.optimize import curve_fit
from scipy.stats import moyal
import matplotlib.image as mpimg
from matplotlib.offsetbox import OffsetImage, AnnotationBbox


root_file = "../final_sim_data/250proton.root"

try:
    # Open the ROOT file and access the 'hibeam' tree
    tree = uproot.open(f"{root_file}:hibeam")
    
    # Extract Energy: Using ProtoTPC.SumEdep as a baseline for total deposited energy.
    # (If you need individual hits instead, you'd pull "target_EDep" and use ak.flatten)
    # Pull the individual hit energies (this creates a jagged array of vectors)
    edep_jagged = tree["ProtoTPC/Edep"].array()

    # Sum the hits along axis 1 (this adds up all hits per event)
    energy = ak.sum(edep_jagged, axis=1)

    # Convert to a standard flat numpy array
    energy = ak.to_numpy(energy)
    
    # Extract Weights: PrimaryWeight is a vector<double>. We usually take the first element per event.
    weights_jagged = tree["PrimaryWeight"].array()
    weights = ak.firsts(weights_jagged) 
    
    # If any event has NO weight, default to 1.0 so we don't drop the data
    weights = ak.fill_none(weights, 1.0)

    # Convert awkward arrays to standard 1D numpy arrays for easy math
    energy = ak.to_numpy(energy)
    weights = ak.to_numpy(weights)

    if len(energy) == 0:
        raise ValueError("Tree is empty.")

except Exception as e:
    print(f"ROOT Loading Error: {e}")
    exit()

# Filter Zeros & Apply the 30% Top-End Cut


# Remove all events where the particle completely missed/deposited no energy
hit_mask = energy > 0
energy_hits = energy[hit_mask]
weights_hits = weights[hit_mask]

print(f"Total events simulated: {len(energy)}")
print(f"Events with actual hits (>0 GeV): {len(energy_hits)}")

# calculate the 70th percentile using ONLY the actual hits
cut_value = np.percentile(energy_hits, 70)

# Create a boolean mask to keep only data below the cut value
cut_mask = energy_hits <= cut_value
e_cut = energy_hits[cut_mask]
w_cut = weights_hits[cut_mask]

print(f"Events after 30% cut: {len(e_cut)} (Cut threshold: {cut_value:.6f} GeV)")

# The weights parameter applies the primary particle weights to the bin counts
counts, edges = np.histogram(e_cut, bins=50, weights=w_cut)
bin_centers = (edges[:-1] + edges[1:]) / 2

# error is the sqrt of the sum of the squares of weights
sumw2, _ = np.histogram(e_cut, bins=edges, weights=w_cut**2)
errors = np.sqrt(sumw2)

# Prevent division by zero in empty bins during the curve_fit
errors[errors == 0] = 1 

if np.sum(counts) == 0:
    print("!Warning! Histogram is empty after cuts. Cannot be fit.")
    exit()


# Landau Fit (Moyal Approximation)

def landau_approx(x, amp, loc, scale):
    return amp * moyal.pdf(x, loc=loc, scale=scale)

# Initial guess
p0 = [max(counts), bin_centers[np.argmax(counts)], 0.05]

# Bounds: amp > 0, loc between min/max of cut data, scale > 0
lower_bounds = [0, 0, 1e-6]
upper_bounds = [np.inf, cut_value, 1] 

try:
    popt, pcov = curve_fit(
        landau_approx, 
        bin_centers, 
        counts, 
        p0=p0, 
        sigma=errors, 
        bounds=(lower_bounds, upper_bounds)
    )
    print(f"Weighted Fit converged: MPV={popt[1]:.4f}, Width={popt[2]:.4f}")
except RuntimeError:
    print("Fit failed to converge. Plotting guess instead.")
    popt = p0 

x_smooth = np.linspace(0, cut_value, 1000)
y_smooth = landau_approx(x_smooth, *popt)
y_fit_bins = landau_approx(bin_centers, *popt)
residuals = counts - y_fit_bins

plt.style.use(hep.style.CMS)
fig, (ax_main, ax_resid) = plt.subplots(
    2, 1, 
    sharex=True, 
    gridspec_kw={'height_ratios': [3, 1], 'hspace': 0}, 
    figsize=(10, 8)
)









# 1. Plot the Data FIRST
hep.histplot(counts, bins=edges, yerr=errors, ax=ax_main, histtype='errorbar', color='black', label='Weighted Data (70% Cut)')

# 2. Plot the Fit Line SECOND
ax_main.plot(x_smooth, y_smooth, color='#E31A1C', linewidth=2.5, label=f'Landau Fit\nMPV={popt[1]:.4f} GeV')
ax_main.axvline(cut_value, color='blue', linestyle=':', label='30% Cut Threshold')

# 3. NOW call the legend (after the plots are drawn!)
ax_main.set_yscale('log')
ax_main.set_ylabel("Weighted Events / Bin")
ax_main.set_ylim(bottom=0.5) 
ax_main.legend(loc='upper right') # Forces the legend to render the labels we just made


# CUSTOM ESS / HIBEAM BRANDING

# Use the CMS labeler to get the perfect font, spacing, and layout
main_text, prelim_text = hep.cms.label(loc=0, data=False, label="Preliminary", rlabel="HIBEAM", ax=ax_main)

# Overwrite the "CMS" text with "ESS"
main_text.set_text("ESS")

# Inject the custom ESS Logo
try:
    logo_img = mpimg.imread("data/ess_logo.png") # <-- Update this path to your PNG!
    imagebox = OffsetImage(logo_img, zoom=0.15)  # <-- Adjust zoom to make it fit nicely
    
    # Place it dynamically (0.03 = 3% from left, 0.95 = 95% from bottom)
    ab = AnnotationBbox(
        imagebox, 
        (0.03, 0.95), 
        xycoords='axes fraction', 
        box_alignment=(0, 1), 
        frameon=False         
    )
    ax_main.add_artist(ab)
except FileNotFoundError:
    print("Logo PNG not found. Skipping image insertion.")
# ---------------------------------------------------------

# 4. Residuals Plot
ax_resid.errorbar(bin_centers, residuals, yerr=errors, fmt='o', color='black', markersize=4)
ax_resid.axhline(0, color='gray', linestyle='--')
ax_resid.set_xlabel("Energy Deposited [GeV]")
ax_resid.set_ylabel("Data - Fit")

limit = max(np.max(np.abs(residuals)), 1) * 1.5
ax_resid.set_ylim(-limit, limit)

plt.show()