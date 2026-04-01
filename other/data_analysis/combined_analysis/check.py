"""
Quick diagnostic — run this before the full combined script.
Plots a raw histogram of dE/dx from the tracks branch.
"""

import ROOT
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import os

# ── Config ────────────────────────────────────────────────────────────────────
DATA_FILE     = "../experimental_data/tracks_centroids_tpc_run_0006-sorted_t0_100.root"
HEADERS_SO    = "recovered_headers/recovered_headers.so"
HEADERS_H     = "recovered_headers/TrackData.h"
OUTPUT        = "diagnostic_dedx.png"
# ──────────────────────────────────────────────────────────────────────────────

ROOT.gInterpreter.AddIncludePath("recovered_headers/")

if os.path.exists(HEADERS_H):
    ROOT.gInterpreter.ProcessLine(f'#include "{HEADERS_H}"')
else:
    ROOT.gInterpreter.Declare("""
        #include <vector>
        class TrackData {
        public:
            int nPoints;
            float slope_xy, intercept_xy, slope_zy, intercept_zy;
            std::vector<int>   y;
            std::vector<float> x, z, charge, sigmas_x, sigmas_z, sigmas_y;
        };
    """)

ROOT.gSystem.Load(HEADERS_SO)

import os   # needed for path check above — move to top in real use

tfile = ROOT.TFile.Open(DATA_FILE)
tree  = tfile.Get("trackingData")

tracks_vec = ROOT.std.vector("TrackData")()
tree.SetBranchAddress("tracks", tracks_vec)

all_dedx     = []
raw_charge   = []   # before dividing by length — sanity check
raw_npoints  = []

for i in range(tree.GetEntries()):
    tree.GetEntry(i)
    for track in tracks_vec:
        n = track.nPoints
        if n < 2:
            continue

        x0, x1 = float(track.x[0]), float(track.x[n-1])
        y0, y1 = float(track.y[0]), float(track.y[n-1])
        z0, z1 = float(track.z[0]), float(track.z[n-1])
        length  = np.sqrt((x1-x0)**2 + (y1-y0)**2 + (z1-z0)**2)
        if length < 1e-3:
            continue
        delta_x = length / n

        for j in range(n):
            adc = float(track.charge[j])
            if adc <= 0:
                continue
            raw_charge.append(adc)
            all_dedx.append(adc / delta_x)

        raw_npoints.append(n)

tfile.Close()

all_dedx   = np.array(all_dedx)
raw_charge = np.array(raw_charge)
raw_npoints = np.array(raw_npoints)

# ── Print quick stats ─────────────────────────────────────────────────────────
print(f"Total dE/dx values : {len(all_dedx)}")
print(f"Unique dE/dx values: {len(np.unique(all_dedx))}  ← should be >> 10")
print(f"Min / Max          : {all_dedx.min():.3f} / {all_dedx.max():.3f}")
print(f"Mean / Median      : {all_dedx.mean():.3f} / {np.median(all_dedx):.3f}")
print(f"Std dev            : {all_dedx.std():.3f}")
print(f"Unique ADC values  : {len(np.unique(raw_charge))}  ← if small, charge is corrupted")
print(f"Unique nPoints     : {sorted(np.unique(raw_npoints))}")

# ── Plot ──────────────────────────────────────────────────────────────────────
fig, axes = plt.subplots(1, 3, figsize=(15, 4))

# 1. Raw dE/dx — if this looks like a few spikes it's still corrupted
axes[0].hist(all_dedx, bins=200, color="steelblue", alpha=0.8)
axes[0].set_title("dE/dx  (all values)")
axes[0].set_xlabel("ADC / mm")
axes[0].set_ylabel("Count")

# 2. Clip top 5% to see the bulk shape clearly
clipped = all_dedx[all_dedx < np.percentile(all_dedx, 95)]
axes[1].hist(clipped, bins=150, color="steelblue", alpha=0.8)
axes[1].set_title("dE/dx  (bottom 95%)")
axes[1].set_xlabel("ADC / mm")

# 3. Raw ADC charge — should be a smooth spectrum, not discrete spikes
axes[2].hist(raw_charge, bins=200, color="darkorange", alpha=0.8)
axes[2].set_title("Raw ADC charge per centroid")
axes[2].set_xlabel("ADC counts")

plt.suptitle("DATA DIAGNOSTIC — check these look like smooth Landau-shaped spectra",
             fontweight="bold")
plt.tight_layout()
plt.savefig(OUTPUT, dpi=150)
print(f"\nSaved → {OUTPUT}")