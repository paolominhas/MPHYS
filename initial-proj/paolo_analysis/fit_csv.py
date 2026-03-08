import matplotlib.pyplot as plt
import mplhep as hep
import numpy as np
import pandas as pd
from scipy.optimize import curve_fit
from scipy.stats import moyal

# (nice) try catch loop for the data errors if there is not numeric data, if no data etc
try:
    df = pd.read_csv("data/muon_long_output.csv")
    # everything is now a number just in case, no text no headers
    data = pd.to_numeric(df['total_energy_gev'], errors='coerce').dropna()
    if len(data) == 0:
        raise ValueError("No valid data found in CSV.")

except Exception as e:
    print(f"Data Loading Error: {e}")
    exit()

# histogram
counts, edges = np.histogram(data, bins=50, range=(0, 1))
bin_centers = (edges[:-1] + edges[1:]) / 2
errors = np.sqrt(counts)
errors[errors == 0] = 1 

# just in ase histogram is empty as it was
if np.sum(counts) == 0:
    print("!Warning! Histogram is empty (all counts are 0). Looks like it cannot be fit.")
    exit()

# the landau fit part in a function
def landau_approx(x, amp, loc, scale):
    return amp * moyal.pdf(x, loc=loc, scale=scale)

# Initial guess
# can assume the peak is about where the maximum count value is (MPV most probable value?)
p0 = [max(counts), bin_centers[np.argmax(counts)], 0.05]

# curve_fit bounds here:
# amplitude amp > 0
# location loc: between 0 and 1
# scale > 0 !IMPORTANT! -ve width causes crash ...
lower_bounds = [0, 0, 1e-6]       # amp, loc, scale
upper_bounds = [np.inf, 1, 1]     # amp, loc, scale

try:
    popt, pcov = curve_fit(
        landau_approx, 
        bin_centers, 
        counts, 
        p0=p0, 
        sigma=errors, 
        bounds=(lower_bounds, upper_bounds) # <--- This prevents the crash
    )
    print(f"Fit converged: MPV={popt[1]:.4f}, Width={popt[2]:.4f}")
except RuntimeError:
    print("Fit failed to converge. Plotting guess instead.")
    popt = p0 

# Calculate curves
x_smooth = np.linspace(0, 1, 1000)
y_smooth = landau_approx(x_smooth, *popt)
y_fit_bins = landau_approx(bin_centers, *popt)
residuals = counts - y_fit_bins

###########################
##          Plot!        ##
###########################

# Note: I need to find a way to use hep but without CERN collaboration labels all over everything
# This is from last meeting (9th Feb 2026)

plt.style.use(hep.style.CMS)
fig, (ax_main, ax_resid) = plt.subplots(2, 1, sharex=True, 
                                        gridspec_kw={'height_ratios': [3, 1], 'hspace': 0}, 
                                        figsize=(10, 8))

# plot
hep.histplot(counts, bins=edges, yerr=True, ax=ax_main, histtype='errorbar', color='black', label='Data')
ax_main.plot(x_smooth, y_smooth, color='#E31A1C', linewidth=2.5, label=f'Landau Fit\nMPV={popt[1]:.2f} GeV')
ax_main.set_yscale('log')
ax_main.set_ylabel("Events / Bin")
ax_main.legend()
ax_main.set_ylim(bottom=0.5) 
hep.cms.label(ax=ax_main, loc=0, data=False, label="Preliminary", rlabel="HIBEAM")

# residuals
ax_resid.errorbar(bin_centers, residuals, yerr=errors, fmt='o', color='black', markersize=4)
ax_resid.axhline(0, color='gray', linestyle='--')
ax_resid.set_xlabel("Energy Deposited [GeV]")
ax_resid.set_ylabel("Data - Fit")

# residual limits
limit = max(np.max(np.abs(residuals)), 1) * 1.5
ax_resid.set_ylim(-limit, limit)

plt.show()
