# Data directory

This directory is gitignored — ROOT files are too large for version control.

## How to populate it

Option A — **Symlink** your existing data directories:

```bash
ln -s /path/to/simulation/root/files    data/simulation/
ln -s /path/to/experimental/root/files  data/experimental/
ln -s /path/to/segmentation/outputs     data/segmentation/
```

Option B — **Copy** the files directly into these subdirectories.

## Expected layout

```
data/
├── simulation/
│   ├── KrakowScatter.root           # Krakow proton scattering (1M events)
│   ├── MuonScatter_fixed.root       # Muon lab MCPL beam
│   └── HIBEAMScatter.root           # Full HIBEAM simulation
│
├── experimental/
│   ├── tracks_centroids_tpc_run_0006-sorted_t0_100.root   # CCB run 6
│   ├── tracks_centroids_tpc_run_0042-sorted_t0.root       # CCB run 42
│   ├── tracks_centroids_chamber-...-1008-2_sorted.root    # Cosmics 500V
│   └── tracks_centroids_chamber-...-1008-4_sorted.root    # Cosmics 243V
│
└── segmentation/
    ├── ana_output_krakow/
    │   ├── ana_geom_axis2_nSec002.root
    │   ├── ana_geom_axis2_nSec004.root
    │   └── ...                      # nSec002 through nSec024
    └── ana_output_muon/
        ├── ana_geom_axis2_nSec002.root
        └── ...
```

If your files have different names, update the paths in `config.yaml`.
