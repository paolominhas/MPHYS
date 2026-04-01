import uproot
import awkward as ak
import numpy as np
import matplotlib.pyplot as plt
import mplhep as hep


hep.style.use(hep.style.ROOT)

def calculate_dEdx(x, y, z, q, truncate_fraction=0.3):
    """
    Calculates dE/dx for a track using the Truncated Mean method.
    """

    dx_vec = np.diff(x)
    dy_vec = np.diff(y)
    dz_vec = np.diff(z)
    
    # 3D path length between hits
    dr = np.sqrt(dx_vec**2 + dy_vec**2 + dz_vec**2)
    
    # Filter out zero distances (duplicate hits) to avoid division by zero
    mask = dr > 0
    dr = dr[mask]
    
    # The charge 'q' corresponds to the hit. We drop the last one to match the 'diff' array size
    dE = q[:-1][mask] 


    dedx_samples = dE / dr

    # 4. Truncated Mean Calculation
    # Sort samples from low to high
    sorted_samples = np.sort(dedx_samples)
    
    # Determine cutoff index (discard top X%)
    n_samples = len(sorted_samples)
    if n_samples < 3: return 0.0 # Too few hits to calculate
    
    cutoff_index = int(n_samples * (1 - truncate_fraction))
    
    # Average the lower fraction
    truncated_mean = np.mean(sorted_samples[:cutoff_index])
    
    return truncated_mean

# --- Main Analysis ---

# 1. Load the flattened data
print("Loading data...")
file = uproot.open("flat_data.root")
tree = file["TPCTree"]

# Read arrays lazily
arrays = tree.arrays(["x", "y", "z", "q"])

# 2. Loop over events to calculate dE/dx
# (We use a loop here because dE/dx logic is complex per track, 
# but this can be vectorized further if needed)
dedx_values = []
track_lengths = []

print("Analyzing tracks...")
for event in arrays:
    # Extract data as numpy arrays
    x = event.x.to_numpy()
    y = event.y.to_numpy()
    z = event.z.to_numpy()
    q = event.q.to_numpy()
    
    if len(x) < 5: continue # Skip empty or short events

    # Sort hits by Z (drift direction) to ensure we walk along the track
    # This assumes the track is roughly parallel to Z or simple.
    sort_idx = np.argsort(z)
    x, y, z, q = x[sort_idx], y[sort_idx], z[sort_idx], q[sort_idx]

    # Calculate Total Track Length
    total_length = np.sqrt((x[-1]-x[0])**2 + (y[-1]-y[0])**2 + (z[-1]-z[0])**2)
    
    # Calculate dE/dx
    val = calculate_dEdx(x, y, z, q)
    
    if val > 0:
        dedx_values.append(val)
        track_lengths.append(total_length)

# 3. Plotting
fig, ax = plt.subplots(figsize=(10, 8))

# 2D Histogram: dE/dx vs Track Length (Bragg Curve-ish look)
h = ax.hist2d(track_lengths, dedx_values, bins=[50, 50], cmap="turbo", cmin=1)

ax.set_xlabel("Track Length [mm]")
ax.set_ylabel("Truncated dE/dx [ADC/mm]")
ax.set_title("Particle Energy Loss Identification")
fig.colorbar(h[3], ax=ax, label="Counts")

plt.show()