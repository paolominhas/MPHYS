import ROOT
import os
import glob
import random
import numpy as np
import matplotlib.pyplot as plt

# 1. SUPPRESS THE MASSIVE RED ERROR WALL
ROOT.gErrorIgnoreLevel = ROOT.kFatal 

# 2. Load the reconstructed C++ headers
if os.path.exists("recovered_headers/recovered_headers.so"):
    ROOT.gSystem.Load("recovered_headers/recovered_headers.so")
else:
    raise FileNotFoundError("recovered_headers.so not found. Run MakeProject first.")

def extract_random_tracks(file_path, num_tracks=100):
    """Extracts 3D coordinates and Energy from random events in the file."""
    f = ROOT.TFile.Open(file_path)
    tree = f.Get("trackingData")
    
    total_events = tree.GetEntries()
    if total_events == 0:
        return None
        
    # Pick random events
    sample_size = min(num_tracks, total_events)
    event_indices = random.sample(range(total_events), sample_size)
    
    track_data = []
    all_energy = [] # To calculate global color scale
    
    for i in event_indices:
        tree.GetEntry(i)
        
        x, y, z, e = [], [], [], []
        for c in tree.centroids:
            # Basic physical cuts: must have Z distance and positive energy
            if c.z > 0 and c.integrated_ADC_amplitude > 0:
                x.append(c.x)
                y.append(c.y)
                z.append(c.z)
                e.append(c.integrated_ADC_amplitude)
                all_energy.append(c.integrated_ADC_amplitude)
                
        if len(z) > 2: # Only keep tracks with at least 3 hits
            # Sort by Z-axis to draw continuous track lines naturally
            sort_idx = np.argsort(z)
            track_data.append({
                'x': np.array(x)[sort_idx],
                'y': np.array(y)[sort_idx],
                'z': np.array(z)[sort_idx],
                'e': np.array(e)[sort_idx]
            })
            
    f.Close()
    return track_data, all_energy

def plot_event_display(track_data, all_energy, filename, output_dir):
    """Creates the standard HEP 3D + 2D projection Event Display."""
    if not track_data:
        return
        
    # Calculate color scale limits (ignoring extreme outliers for a better visual)
    vmin = np.percentile(all_energy, 5)
    vmax = np.percentile(all_energy, 95)
    
    # Create a nice wide figure for the multi-view
    fig = plt.figure(figsize=(18, 8))
    fig.suptitle(f'Event Display: {filename} ({len(track_data)} Random Tracks)', fontsize=16)
    
    # --- Panel 1: 3D View ---
    ax_3d = fig.add_subplot(1, 2, 1, projection='3d')
    ax_3d.set_title("3D Track Topology")
    ax_3d.set_xlabel("Z (Drift/Beam Axis)")
    ax_3d.set_ylabel("X")
    ax_3d.set_zlabel("Y")
    
    # --- Panel 2: X-Z Top View ---
    ax_xz = fig.add_subplot(2, 2, 2)
    ax_xz.set_title("Top View (X vs Z)")
    ax_xz.set_ylabel("X Axis")
    ax_xz.grid(True, alpha=0.3)
    
    # --- Panel 3: Y-Z Side View ---
    ax_yz = fig.add_subplot(2, 2, 4, sharex=ax_xz) # Share Z axis with XZ view
    ax_yz.set_title("Side View (Y vs Z)")
    ax_yz.set_xlabel("Z Axis (mm)")
    ax_yz.set_ylabel("Y Axis")
    ax_yz.grid(True, alpha=0.3)
    
    # Plot each track
    # HEP standard uses 'turbo' or 'jet' for energy heatmaps
    cmap = 'turbo' 
    
    for trk in track_data:
        # 3D
        ax_3d.plot(trk['z'], trk['x'], trk['y'], color='gray', alpha=0.2, lw=0.8)
        ax_3d.scatter(trk['z'], trk['x'], trk['y'], c=trk['e'], cmap=cmap, vmin=vmin, vmax=vmax, s=15, alpha=0.9)
        
        # Top (X-Z)
        ax_xz.plot(trk['z'], trk['x'], color='gray', alpha=0.2, lw=0.8)
        sc = ax_xz.scatter(trk['z'], trk['x'], c=trk['e'], cmap=cmap, vmin=vmin, vmax=vmax, s=10, alpha=0.9)
        
        # Side (Y-Z)
        ax_yz.plot(trk['z'], trk['y'], color='gray', alpha=0.2, lw=0.8)
        ax_yz.scatter(trk['z'], trk['y'], c=trk['e'], cmap=cmap, vmin=vmin, vmax=vmax, s=10, alpha=0.9)

   
    cbar_ax = fig.add_axes([0.92, 0.15, 0.015, 0.7]) # [left, bottom, width, height]
    cbar = fig.colorbar(sc, cax=cbar_ax)
    cbar.set_label('Energy Loss (ADC Amplitude)', rotation=270, labelpad=20)
    
    
    plt.subplots_adjust(right=0.9)
    out_path = os.path.join(output_dir, filename.replace(".root", "_tracks.png"))
    plt.savefig(out_path, dpi=200)
    plt.close(fig)
    print(f"  -> Saved {out_path}")

def run_batch_displays(input_dir, output_dir):
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
        
    root_files = glob.glob(os.path.join(input_dir, "*.root"))
    print(f"Found {len(root_files)} files for track visualization.")
    
    for file in root_files:
        filename = os.path.basename(file)
        print(f"Processing: {filename}...")
        
        data = extract_random_tracks(file, num_tracks=100)
        if data is not None and len(data[0]) > 0:
            plot_event_display(data[0], data[1], filename, output_dir)
        else:
            print(f"  -> Skipped {filename}: No valid physical tracks found.")

if __name__ == "__main__":
    INPUT_DIR = "../experimental_data/"
    OUTPUT_DIR = "track_displays/"
    
    run_batch_displays(INPUT_DIR, OUTPUT_DIR)