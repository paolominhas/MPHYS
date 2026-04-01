import uproot
import numpy as np
import matplotlib.pyplot as plt


PITCH_X = 1.0  # mm per column
PITCH_Y = 1.0  # mm per row
PITCH_Z = 1.0  # mm per timestamp tick (v_drift * sampling_time)


file = uproot.open("../initial-proj/hibeam_g4/data/TPC reconstruction ROOT files/tracks_centroids_tpc_run_0042-sorted_t0.root")
tree = file["trackingData"] # Use your working cycle

print("Loading raw TPC data...")
data = tree.arrays([
    "tpc/tpc.row", "tpc/tpc.column", "tpc/tpc.timestamp", 
    "tpc/tpc.val", "tpc/tpc.pedestal"
], entry_stop=5000) # Load 5000 events to get a good distribution

# Lists to store the final results for the histogram
all_dedx = []
all_track_lengths = []

print("Reconstructing tracks and calculating dE/dx...")

for event in data:
    row = event["tpc/tpc.row"].to_numpy()
    col = event["tpc/tpc.column"].to_numpy()
    time = event["tpc/tpc.timestamp"].to_numpy()
    signal = event["tpc/tpc.val"].to_numpy() - event["tpc/tpc.pedestal"].to_numpy()
    
    # Filter noise
    mask = signal > 15
    row, col, time, signal = row[mask], col[mask], time[mask], signal[mask]
    
    if len(row) < 10: 
        continue # Skip events with too few hits to form a track

    unique_rows = np.unique(row)
    
    cluster_x = []
    cluster_y = []
    cluster_z = []
    cluster_E = []
    
    for r in unique_rows:
        # Find all hits in this specific row
        idx = (row == r)
        hits_col = col[idx]
        hits_time = time[idx]
        hits_sig = signal[idx]
        
        total_E = np.sum(hits_sig)
        if total_E == 0: continue
        
        # Charge-weighted average to find the center of the cluster
        avg_col = np.sum(hits_col * hits_sig) / total_E
        avg_time = np.sum(hits_time * hits_sig) / total_E
        
        # Convert to physical units (mm)
        cluster_x.append(avg_col * PITCH_X)
        cluster_y.append(r * PITCH_Y)
        cluster_z.append(avg_time * PITCH_Z)
        cluster_E.append(total_E)
        
    cluster_x = np.array(cluster_x)
    cluster_y = np.array(cluster_y)
    cluster_z = np.array(cluster_z)
    cluster_E = np.array(cluster_E)
    
    # We need at least 5 clusters to get a meaningful truncated mean
    if len(cluster_E) < 5:
        continue

    # --- CALCULATE dE/dx ---
    # Calculate 3D distances between consecutive clusters
    dx = np.diff(cluster_x)
    dy = np.diff(cluster_y)
    dz = np.diff(cluster_z)
    
    # Distance formula: dx_3D = sqrt(dx^2 + dy^2 + dz^2)
    dr = np.sqrt(dx**2 + dy**2 + dz**2)
    
  
    dr_mask = dr > 0
    dr = dr[dr_mask]
    dE = cluster_E[:-1][dr_mask] 
    
    if len(dr) < 3: continue
    

    dedx_samples = dE / dr
    
    # Truncated Mean: Sort and drop the top 30%
    sorted_samples = np.sort(dedx_samples)
    cutoff = int(len(sorted_samples) * 0.70)
    truncated_mean = np.mean(sorted_samples[:cutoff])
    
    # Calculate total track length
    total_length = np.sum(dr)
    
    all_dedx.append(truncated_mean)
    all_track_lengths.append(total_length)

# ==========================================
# PLOTTING
# ==========================================
fig, ax = plt.subplots(figsize=(8, 6))

# Plot a 1D histogram of the dE/dx values
ax.hist(all_dedx, bins=50, color='royalblue', edgecolor='black', alpha=0.7)
ax.set_xlabel("Truncated dE/dx [ADC/mm]")
ax.set_ylabel("Number of Tracks")
ax.set_title("Muon Energy Loss (dE/dx) from Raw TPC Data")

plt.grid(axis='y', linestyle='--', alpha=0.7)
plt.show()