import ROOT
import os
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import linregress

# 1. Load the reconstructed library
if os.path.exists("recovered_headers/recovered_headers.so"):
    ROOT.gSystem.Load("recovered_headers/recovered_headers.so")
else:
    print("Error: Library not found. Please run MakeProject again.")
    exit()

# 2. Open file and tree
f = ROOT.TFile.Open("../experimental_data/muon_data.root")
tree = f.Get("trackingData")

# 3. Separate lists for different data types
centroid_list = []
track_list = []

print(f"Analyzing {tree.GetEntries()} events...")

# 4. Process the tree
for i in range(tree.GetEntries()):
    tree.GetEntry(i)
    
    # --- Process Centroids ---
    for c in tree.centroids:
        centroid_list.append({
            'event_id': i,
            'x': c.x, 'y': c.y, 'z': c.z,
            'sigma_x': c.sigma_x,
            'adc': c.integrated_ADC_amplitude
        })
        
    # --- Process Tracks ---
    # Assuming TrackData has common fields like theta, phi, length, or chi2
    for t in tree.tracks:
        # Note: You can use dir(t) to see exact field names if these differ
        track_list.append({
            'event_id': i,
            'track_id': len(track_list), # Unique ID for this track
            # Replace these with the actual fields found in your TrackData
            'chi2': getattr(t, 'chi2', 0), 
            'ndf': getattr(t, 'ndf', 0),
            'length': getattr(t, 'length', 0)
        })

# 5. Create DataFrames
df_centroids = pd.DataFrame(centroid_list)
df_tracks = pd.DataFrame(track_list)

# 6. Display Results
print("\n--- Centroids Summary ---")
print(df_centroids.head())
print(f"Total Centroids: {len(df_centroids)}")

print("\n--- Tracks Summary ---")
print(df_tracks.head())
print(f"Total Tracks: {len(df_tracks)}")

# 7. Save to HDF5 (Better than CSV for physics data)
# df_centroids.to_hdf('muon_analysis.h5', key='centroids')
# df_tracks.to_hdf('muon_analysis.h5', key='tracks')

# 1. Calculate dE/dx
# For a basic approach, we look at ADC (Energy) per unit Z (as a proxy for distance)
df_centroids['dEdZ'] = df_centroids['adc'] / df_centroids['z']

# 2. Linear Regression to find the "Expected" energy loss trend
# We'll filter out zeros or extreme outliers first
clean_df = df_centroids[(df_centroids['z'] > 0) & (df_centroids['adc'] > 100)]
slope, intercept, r_value, p_value, std_err = linregress(clean_df['z'], clean_df['adc'])

# 3. Calculate Residuals
# Residual = Observed ADC - Predicted ADC from the fit
clean_df['predicted_adc'] = intercept + slope * clean_df['z']
clean_df['residual'] = clean_df['adc'] - clean_df['predicted_adc']

# --- Plotting ---
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 12), sharex=True)

# Plot A: Energy vs Z with Fit
ax1.scatter(clean_df['z'], clean_df['adc'], alpha=0.1, s=1, label='Centroids')
ax1.plot(clean_df['z'], clean_df['predicted_adc'], color='red', label='Linear Fit')
ax1.set_ylabel('Energy (ADC Units)')
ax1.set_title('Energy Deposit vs. Drift Distance (Z)')
ax1.legend()

# Plot B: Residuals
ax2.scatter(clean_df['z'], clean_df['residual'], alpha=0.1, s=1, color='purple')
ax2.axhline(0, color='black', linestyle='--')
ax2.set_xlabel('Z (mm)')
ax2.set_ylabel('Residual (ADC)')
ax2.set_title('Residuals of Energy Fit')

plt.tight_layout()
plt.savefig("dedx_analysis.png")
print("Plots saved as dedx_analysis.png")

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import moyal

# 1. Prepare Data & Mitigate Backgrounds
# Apply a threshold cut to remove low-energy noise
threshold = 500  # Adjust based on your ADC distribution
analysis_df = df_centroids[df_centroids['adc'] > threshold].copy()

# Calculate dE/dx (Using Z as a path length proxy for now)
# Real dx would be sqrt(dx^2 + dy^2 + dz^2)
analysis_df['dedx'] = analysis_df['adc'] / analysis_df['z']

# Drop infinities or NaNs
data = analysis_df['dedx'].dropna()
data = data[data < np.percentile(data, 95)] # Cut extreme outliers for better fit

# 2. Fit to a Landau (Moyal) Distribution
# Moyal parameters: loc (peak position), scale (width)
params = moyal.fit(data)

# 3. Plotting
plt.figure(figsize=(10, 6))

# Histogram of data
count, bins, ignored = plt.hist(data, bins=100, density=True, 
                                alpha=0.6, color='skyblue', label='Experimental Data')

# Plot the Fit
x = np.linspace(min(bins), max(bins), 500)
plt.plot(x, moyal.pdf(x, *params), 'r-', lw=2, 
         label=f'Landau (Moyal) Fit\nMPV: {params[0]:.2f}\nWidth: {params[1]:.2f}')

plt.title('dE/dx Energy Loss Distribution (Muon Data)')
plt.xlabel('dE/dx (ADC/mm)')
plt.ylabel('Probability Density')
plt.legend()
plt.grid(axis='y', alpha=0.3)

plt.savefig("landau_fit.png")
print(f"Fit complete. Most Probable Value (MPV): {params[0]}")

# Calculate residuals between histogram and fit
bin_centers = (bins[:-1] + bins[1:]) / 2
fit_counts = moyal.pdf(bin_centers, *params)
# Normalize fit_counts to match histogram area if density=False was used
residuals = count - fit_counts

plt.figure(figsize=(10, 3))
plt.bar(bin_centers, residuals, width=(bins[1]-bins[0]), color='purple')
plt.axhline(0, color='black', linestyle='--')
plt.title('Fit Residuals')
plt.ylabel('Data - Fit')
plt.show()