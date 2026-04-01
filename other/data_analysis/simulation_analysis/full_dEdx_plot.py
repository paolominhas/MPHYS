import uproot
import awkward as ak
import matplotlib.pyplot as plt
import mplhep as hep
import numpy as np
from scipy.optimize import curve_fit
from scipy.stats import moyal
import matplotlib.image as mpimg
from matplotlib.offsetbox import OffsetImage, AnnotationBbox


root_file = "../simulation_data/mcpl_long_muon.root"

try:
    tree = uproot.open(f"{root_file}:hibeam")
    edep_jagged = tree["ProtoTPC/Edep"].array()
    energy = ak.to_numpy(ak.sum(edep_jagged, axis=1))
    weights_jagged = tree["PrimaryWeight"].array()
    weights = ak.firsts(weights_jagged)
    weights = ak.fill_none(weights, 1.0)
    energy = ak.to_numpy(energy)
    weights = ak.to_numpy(weights)
    if len(energy) == 0:
        raise ValueError("Tree is empty.")
except Exception as e:
    print(f"ROOT Loading Error: {e}")
    exit()

# Remove zero-energy events
hit_mask     = energy > 0
energy_hits  = energy[hit_mask]
weights_hits = weights[hit_mask]

print(f"Total events simulated: {len(energy)}")
print(f"Events with actual hits (>0 GeV): {len(energy_hits)}")

# Remove extreme outliers (top 1%) before doing anything — these are
# delta rays or secondary interactions, not the primary dE/dx signal.
# This is standard practice in TPC dE/dx analysis — see Allison & Cobb (1980)
# and the ALICE dE/dx procedure (ALICE-PUBLIC-2017-005).
outlier_threshold = np.percentile(energy_hits, 99)
clean_mask    = energy_hits <= outlier_threshold
energy_hits   = energy_hits[clean_mask]
weights_hits  = weights_hits[clean_mask]
print(f"After 1% outlier removal   : {len(energy_hits)}  (removed above {outlier_threshold:.6f} GeV)")

# 70th percentile cut — removes the top 30% highest energy deposits.
# This is the standard truncated mean method for dE/dx measurement:
# the Landau tail is dominated by rare hard delta-ray collisions which
# inflate the mean; removing the top 30% recovers a Gaussian-like
# estimator much closer to the MPV. Reference: Allison & Cobb (1980),
# also used by NA49, ALICE TPC, STAR TPC.
cut_value = np.percentile(energy_hits, 70)

cut_mask   = energy_hits <= cut_value
above_mask = energy_hits >  cut_value

e_cut      = energy_hits[cut_mask];   w_cut   = weights_hits[cut_mask]
e_above    = energy_hits[above_mask]; w_above = weights_hits[above_mask]

print(f"Events after 30% cut: {len(e_cut)} (Cut threshold: {cut_value:.6f} GeV)")

# Fitted region histogram — bin ONLY over the fitted data's natural range
# (same as original). This keeps the bin width fine enough for a good fit.
counts, edges = np.histogram(e_cut, bins=50, weights=w_cut)
bin_centers   = (edges[:-1] + edges[1:]) / 2
sumw2, _      = np.histogram(e_cut, bins=edges, weights=w_cut**2)
errors        = np.sqrt(sumw2)
errors[errors == 0] = 1

# Extend bin edges into the tail using the same bin width so red points align
bin_width      = edges[1] - edges[0]
tail_edges     = np.arange(edges[0], outlier_threshold + bin_width, bin_width)
tail_counts, _ = np.histogram(e_above, bins=tail_edges, weights=w_above)
sumw2_t, _     = np.histogram(e_above, bins=tail_edges, weights=w_above**2)
errors_tail    = np.sqrt(sumw2_t)
errors_tail[errors_tail == 0] = 1

if np.sum(counts) == 0:
    print("!Warning! Histogram is empty after cuts. Cannot be fit.")
    exit()

# Landau Fit (Moyal Approximation)
def landau_approx(x, amp, loc, scale):
    return amp * moyal.pdf(x, loc=loc, scale=scale)

p0           = [max(counts), bin_centers[np.argmax(counts)], 0.05]
lower_bounds = [0, 0, 1e-6]
upper_bounds = [np.inf, cut_value, 1]

try:
    popt, pcov = curve_fit(
        landau_approx, bin_centers, counts,
        p0=p0, sigma=errors,
        bounds=(lower_bounds, upper_bounds))
    print(f"Weighted Fit converged: MPV={popt[1]:.4f}, Width={popt[2]:.4f}")
except RuntimeError:
    print("Fit failed to converge. Plotting guess instead.")
    popt = p0

x_smooth   = np.linspace(0, outlier_threshold, 1000)
y_smooth   = landau_approx(x_smooth, *popt)
y_fit_bins = landau_approx(bin_centers, *popt)
residuals  = counts - y_fit_bins

plt.style.use(hep.style.CMS)
fig, (ax_main, ax_resid) = plt.subplots(
    2, 1,
    sharex=True,
    gridspec_kw={'height_ratios': [3, 1], 'hspace': 0},
    figsize=(10, 8))

# 1. Truncated (excluded) points in red — plotted first so black sits on top
hep.histplot(tail_counts, bins=tail_edges, yerr=errors_tail, ax=ax_main,
             histtype='errorbar', color='#E31A1C', alpha=0.7,
             label='Excluded data (top 30%)')

# 2. Fitted data in black
hep.histplot(counts, bins=edges, yerr=errors, ax=ax_main,
             histtype='errorbar', color='black',
             label='Weighted Data (70% Cut)')

# 3. Fit line
ax_main.plot(x_smooth, y_smooth, color='#E31A1C', linewidth=2.5,
             label=f'Landau Fit\nMPV={popt[1]:.4f} GeV')
ax_main.axvline(cut_value, color='blue', linestyle=':',
                label='30% Cut Threshold')

ax_main.set_yscale('log')
ax_main.set_ylabel("Weighted Events / Bin")
ax_main.set_ylim(bottom=0.5)
ax_main.set_xlim(0, outlier_threshold)
ax_main.legend(loc='upper right')

label_artists = hep.cms.label(loc=0, data=False,
                              text="Preliminary", rlabel="HIBEAM",
                              ax=ax_main)
label_artists[0].set_text("ESS")

try:
    logo_img = mpimg.imread("data/ess_logo.png")
    imagebox = OffsetImage(logo_img, zoom=0.15)
    ab = AnnotationBbox(imagebox, (0.03, 0.95),
                        xycoords='axes fraction',
                        box_alignment=(0, 1), frameon=False)
    ax_main.add_artist(ab)
except FileNotFoundError:
    print("Logo PNG not found. Skipping image insertion.")

# 4. Residuals — fitted region only
ax_resid.errorbar(bin_centers, residuals, yerr=errors,
                  fmt='o', color='black', markersize=4)
ax_resid.axhline(0, color='gray', linestyle='--')
ax_resid.set_xlabel("Energy Deposited [GeV]")
ax_resid.set_ylabel("Data - Fit")
limit = max(np.max(np.abs(residuals)), 1) * 1.5
ax_resid.set_ylim(-limit, limit)

plt.savefig("E_plot.png", dpi=150, bbox_inches="tight")
print("Saved → E_plot.png")
plt.show()