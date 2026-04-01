"""
nPoints vs ADC — using tpc branch (fully split, uproot-readable)
tpc.val = raw ADC value per hit
tpc.row = pad row per hit
nPoints per event = count of hits above noise floor
"""

import uproot
import awkward as ak
import numpy as np
import glob, os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

INPUT_DIR   = "../experimental_data/"
OUTPUT_PATH = "analysis_plots/npoints_vs_adc.png"
MIN_ADC     = 10   # noise floor on raw tpc.val

files = sorted(glob.glob(os.path.join(INPUT_DIR, "*.root")))
if not files:
    raise FileNotFoundError(f"No .root files in {INPUT_DIR}")

os.makedirs("analysis_plots", exist_ok=True)

COLORS = ["#2E86AB", "#E84855", "#F4A261", "#3BB273", "#7B2D8B"]

fig, ax = plt.subplots(figsize=(10, 7))
fig.patch.set_facecolor("white")

for idx, path in enumerate(files):
    label = os.path.basename(path).replace(".root","") \
                                  .replace("tracks_centroids_","") \
                                  .replace("_sorted","")[:45]
    try:
        with uproot.open(path) as f:
            tree = f["trackingData"]
            val  = tree["tpc/tpc.val"].array()
            row  = tree["tpc/tpc.row"].array()

        above   = val > MIN_ADC
        adc_sum = ak.to_numpy(ak.sum(val * above,  axis=1)).astype(float)
        n_pts   = ak.to_numpy(ak.sum(above,         axis=1)).astype(int)

        mask    = n_pts > 0
        adc_sum = adc_sum[mask]
        n_pts   = n_pts[mask]

        print(f"{label}: {mask.sum()} events, "
              f"max ADC={adc_sum.max():.0f}, max nPts={n_pts.max()}")

        c = COLORS[idx % len(COLORS)]
        ax.scatter(adc_sum, n_pts, alpha=0.20, s=5,
                   color=c, rasterized=True, label=label)

        # Profile: mean nPts in equal-occupancy ADC bins
        edges = np.unique(np.percentile(adc_sum, np.linspace(2, 98, 22)))
        px, py, pe = [], [], []
        for lo, hi in zip(edges[:-1], edges[1:]):
            sel = n_pts[(adc_sum >= lo) & (adc_sum < hi)]
            if len(sel) < 5:
                continue
            px.append(0.5*(lo+hi))
            py.append(sel.mean())
            pe.append(sel.std() / np.sqrt(len(sel)))
        ax.errorbar(px, py, yerr=pe, fmt="o-", color=c,
                    lw=2, ms=4, zorder=5)

    except Exception as e:
        print(f"  ✗ {label}: {e}")

ax.set_xlabel("Summed raw ADC per event  [ADC counts]", fontsize=12)
ax.set_xbound(0, 750000)
ax.set_ylabel("N hits above threshold per event", fontsize=12)
ax.set_ybound(0,4000)
ax.set_title(f"Hit Count vs. Total ADC  (all runs, ADC threshold = {MIN_ADC})",
             fontsize=12)
ax.legend(fontsize=7, loc="upper left", markerscale=3)
ax.grid(alpha=0.3)
plt.tight_layout()
plt.savefig(OUTPUT_PATH, dpi=180, bbox_inches="tight")
plt.close()
print(f"\n✔  Saved → {OUTPUT_PATH}")