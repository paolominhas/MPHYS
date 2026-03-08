import uproot
import awkward as ak
import numpy as np
import matplotlib.pyplot as plt


file = uproot.open("../initial-proj/hibeam_g4/data/muon_data.root")
tree = file["trackingData"] 
print("Loading raw TPC data...")
data = tree.arrays(["tpc/tpc.row", "tpc/tpc.column", "tpc/tpc.timestamp", 
                    "tpc/tpc.val", "tpc/tpc.pedestal"], entry_stop=10)


fig = plt.figure(figsize=(10, 8))
ax = fig.add_subplot(111, projection='3d')

for i, event in enumerate(data):
    # Extract arrays for this specific event
    row = event["tpc/tpc.row"]
    col = event["tpc/tpc.column"]
    time = event["tpc/tpc.timestamp"]
    val = event["tpc/tpc.val"]
    ped = event["tpc/tpc.pedestal"]
    
    signal = val - ped
    
    signal_mask = signal > 15 
    
    # Apply the mask to get clean hits
    clean_x = col[signal_mask]
    clean_y = row[signal_mask]
    clean_z = time[signal_mask]
    clean_energy = signal[signal_mask]
    
    if len(clean_energy) == 0:
        continue
        
    sc = ax.scatter(clean_x, clean_z, clean_y, 
                    c=clean_energy, cmap='viridis', 
                    s=clean_energy/5, alpha=0.7)

ax.set_xlabel('Column (X)')
ax.set_ylabel('Timestamp (Z / Drift Direction)')
ax.set_zlabel('Row (Y)')
ax.set_title('Raw 3D TPC Hits (Muon Beam)')
fig.colorbar(sc, ax=ax, label='Signal (val - pedestal)')

plt.show()