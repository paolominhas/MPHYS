import uproot
import awkward as ak
import numpy as np
import matplotlib.pyplot as plt

# 1. Load the file and tree (grabbing the working cycle)
file = uproot.open("../initial-proj/hibeam_g4/data/TPC reconstruction ROOT files/tracks_centroids_chamber-fec01-debug-clk0-ped0-thr0-g0-s0-p1-1008-2_sorted.root")
tree = file["trackingData"] 

# 2. Extract the raw arrays
# We only need a few events to visualize, so we use entry_stop=20
print("Loading raw TPC data...")
data = tree.arrays([
    "tpc/tpc.row", 
    "tpc/tpc.column", 
    "tpc/tpc.timestamp", 
    "tpc/tpc.val", 
    "tpc/tpc.pedestal"
], entry_stop=20)

# 3. Setup the 3D plot
fig = plt.figure(figsize=(10, 8))
ax = fig.add_subplot(111, projection='3d')

print("Processing events...")
# 4. Loop through the events
for i, event in enumerate(data):
    # Extract arrays for this specific event
    row = event["tpc/tpc.row"]
    col = event["tpc/tpc.column"]
    time = event["tpc/tpc.timestamp"]
    val = event["tpc/tpc.val"]
    ped = event["tpc/tpc.pedestal"]
    
    # Calculate the actual signal (Energy = ADC value - Pedestal noise)
    signal = val - ped
    
    # Filter out the electronic noise (only keep hits with a clear signal)
    # 15 is a guess; you may need to increase or decrease this threshold!
    signal_mask = signal > 15 
    
    # Apply the mask to get clean hits
    clean_x = col[signal_mask]
    clean_y = row[signal_mask]
    clean_z = time[signal_mask]
    clean_energy = signal[signal_mask]
    
    # Skip events that are empty after filtering noise
    if len(clean_energy) == 0:
        continue
        
    # Plot the 3D hits. 
    # The color is weighted by the energy deposited.
    sc = ax.scatter(clean_x, clean_z, clean_y, 
                    c=clean_energy, cmap='turbo', 
                    s=15, alpha=0.8)

ax.set_xlabel('Column (X)')
ax.set_ylabel('Timestamp (Z / Drift Direction)')
ax.set_zlabel('Row (Y)')
ax.set_title('Raw 3D TPC Hits (Muon Beam)')

# Add colorbar only once
if 'sc' in locals():
    fig.colorbar(sc, ax=ax, label='Signal (ADC counts)')

plt.show()