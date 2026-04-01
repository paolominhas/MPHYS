import ROOT
import uproot
import awkward as ak
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import moyal

# ── Config ────────────────────────────────────────────────────────────────────
DATA_FILE = "../experimental_data/tracks_centroids_tpc_run_0042-sorted_t0.root"
SIM_FILE  = "../simulation_data/250proton.root"
OUTPUT    = "comparison.png"

# ── Load data ─────────────────────────────────────────────────────────────────
ROOT.gInterpreter.AddIncludePath("recovered_headers/")
ROOT.gInterpreter.ProcessLine('#include "TrackData.h"')
ROOT.gSystem.Load("recovered_headers/recovered_headers.so")

tfile = ROOT.TFile.Open(DATA_FILE)
tree  = tfile.Get("trackingData")
vec   = ROOT.std.vector("TrackData")()
tree.SetBranchAddress("tracks", vec)

data_dedx = []
for i in range(tree.GetEntries()):
    tree.GetEntry(i)
    for track in vec:
        n = track.nPoints
        if n < 3:
            continue
        length = np.sqrt((float(track.x[n-1]) - float(track.x[0]))**2 +
                         (float(track.y[n-1]) - float(track.y[0]))**2 +
                         (float(track.z[n-1]) - float(track.z[0]))**2)
        if length < 1e-3:
            continue
        dx = length / n
        for j in range(n):
            adc = float(track.charge[j])
            if adc > 0:
                data_dedx.append(adc / dx)
tfile.Close()

data_dedx = np.array(data_dedx)
print(f"Data: {len(data_dedx)} values, {len(np.unique(data_dedx))} unique")

# ── Load simulation ───────────────────────────────────────────────────────────
tree_sim   = uproot.open(f"{SIM_FILE}:hibeam")
edep       = ak.to_numpy(ak.sum(tree_sim["ProtoTPC/Edep"].array(), axis=1))
weights    = ak.to_numpy(ak.fill_none(ak.firsts(tree_sim["PrimaryWeight"].array()), 1.0))
sim_energy = edep[edep > 0]
sim_w      = weights[edep > 0]
print(f"Sim:  {len(sim_energy)} events with hits")

# ── Peak-normalise so both sit on the same x-axis ────────────────────────────
d_hist, d_edges = np.histogram(data_dedx, bins=200)
data_peak = ((d_edges[:-1] + d_edges[1:]) / 2)[np.argmax(d_hist)]

s_hist, s_edges = np.histogram(sim_energy, bins=200, weights=sim_w)
sim_peak  = ((s_edges[:-1] + s_edges[1:]) / 2)[np.argmax(s_hist)]

data_norm = data_dedx  / data_peak
sim_norm  = sim_energy / sim_peak

print(f"Data peak: {data_peak:.3f} ADC/mm")
print(f"Sim  peak: {sim_peak:.5f} GeV")

# ── Plot ──────────────────────────────────────────────────────────────────────
x_max = min(np.percentile(data_norm, 99), np.percentile(sim_norm, 99))
bins  = np.linspace(0, x_max, 80)

fig, ax = plt.subplots(figsize=(10, 6))

ax.hist(data_norm, bins=bins, density=True,
        histtype="step", color="black", lw=1.5, label="Data")

ax.hist(sim_norm,  bins=bins, density=True, weights=sim_w,
        histtype="stepfilled", color="steelblue", alpha=0.4, label="Simulation")

ax.set_xlabel("dE/dx  [units of peak]")
ax.set_ylabel("Probability Density")
ax.set_yscale("log")
ax.set_ylim(bottom=1e-3)
ax.legend()
ax.set_title(f"Data vs Simulation  |  data peak={data_peak:.1f} ADC/mm  |  sim peak={sim_peak:.4f} GeV")

plt.tight_layout()
plt.savefig(OUTPUT, dpi=150)
print(f"Saved → {OUTPUT}")