# HIBEAM-MPHYS

<p align="center">
  <img src="analysis/assets/ess_logo.png" alt="ESS" width="200"/>
</p>

<h1 align="center">HIBEAM Prototype TPC</h1>

<p align="center">
  <em>Simulation, reconstruction, and analysis of the prototype Time Projection Chamber (and HIBEAM Annihilation Detector)<br/>
  for the HIBEAM experiment at the European Spallation Source</em>
</p>

<p align="center">
  <a href="#quick-start">Quick start</a> · 
  <a href="#project-structure">Structure</a> · 
  <a href="#pipeline">Pipeline</a> · 
  <a href="#results">Results</a> · 
  <a href="#references">References</a>
</p>

---

## Overview

This repository contains the full simulation and data analysis chain for the HIBEAM prototype TPC detector.  It covers everything from Geant4 geometry construction through Monte Carlo particle transport to publication-quality energy-loss analysis of both simulated and experimental beam-test data.

The project produced dE/dx measurements for two beam configurations at the IFJ Krakow CCB facility, validated against Geant4 simulation using Landau (Moyal) fitting with the standard 70% truncated-mean method.

---

## Project structure

```
HIBEAM/
│
├── initial-proj/                           Simulation chain (C++/Geant4)
│   ├── hibeam_g4_geobuilder/              ── Detector geometry builder
│   ├── hibeam_g4/                         ── Geant4 Monte Carlo simulation
│   └── hibeam_g4_analysis/                ── Raw → processed ROOT converter
│
└── analysis/                    Data analysis (Python)
    ├── config.yaml                        ── ALL tuneable parameters
    ├── hibeam/                            ── Importable Python package
    │   ├── io/                               Data loaders
    │   ├── physics/                          dE/dx, fitting, Bethe-Bloch
    │   └── plotting/                         Publication figures
    ├── scripts/                           ── CLI entry points
    ├── notebooks/                         ── Jupyter exploration
    └── output/                            ── Generated figures
```

---

## Pipeline

The data flows through five stages, from geometry definition to publication figures:

### Stage 1 — Geometry

`hibeam_g4_geobuilder` constructs the detector geometry using ROOT's TGeo classes.  It builds the complete HIBEAM annihilation setup: carbon foil, beam pipe, TPC, WASA calorimeter barrel, cosmic veto, and beam dump.

```bash
cd initial-proj/hibeam_g4_geobuilder_build/
./geobuilder --hibeam --save_as=geometry.root
```

### Stage 2 — Simulation

`hibeam_g4` runs the Geant4 Monte Carlo, transporting particles through the geometry and recording energy deposits in the TPC and calorimeter.  Supports GPS, MCPL, and proton-deuteron elastic scattering sources.

```bash
cd initial-proj/hibeam_g4_build/
./hibeam_g4 -m run.mac -c config.par output.root
```

### Stage 3 — Processing

`hibeam_g4_analysis` converts raw Geant4 output into processed ROOT ntuples with per-event energy sums, pad assignments, and hit counts.  Automatically detects which detector branches are present.

```bash
cd initial-proj/hibeam_g4_analysis_build/
./hibeam_ana --in=raw_output.root --out=processed.root
```

### Stage 4 — Analysis

The `hibeam` Python package reads both processed simulation files and experimental beam-test ROOT files, computes dE/dx using multiple methods, fits Landau distributions, and produces all figures.  Every parameter is controlled by a single `config.yaml`.

```bash
cd hibeam-tpc-analysis/
pip install -e .
python scripts/run_all.py
```

### Stage 5 — Figures

All output lands in `hibeam-tpc-analysis/output/` as both PDF (vector, for publication) and PNG (raster, for slides).

---

## Results

### Energy-loss distributions

Landau (Moyal) fits to the truncated dE/dx spectrum for simulation and experimental data, with pull residuals and goodness-of-fit statistics.

<p align="center">
  <img src="hibeam-tpc-analysis/output/sim/krakow_dedx_linear.png" width="48%"/>
  <img src="hibeam-tpc-analysis/output/sim/krakow_dedx_log.png" width="48%"/>
</p>
<p align="center"><em>Krakow proton scattering simulation — linear (left) and log (right) y-axes.<br/>
The Landau peak and truncation boundary are clearly visible; the excluded delta-ray tail is hatched.</em></p>

### Data vs simulation comparison

Experimental and simulated dE/dx distributions overlaid on a common peak-normalised axis, with independent Landau fits and residuals.

<p align="center">
  <img src="hibeam-tpc-analysis/output/comparison/run_0042_vs_krakow.png" width="60%"/>
</p>
<p align="center"><em>CCB run 42 (data, black) vs Krakow simulation (blue).  Both normalised to their respective peak MPV.</em></p>

### Segmentation study

Impact of detector segmentation on energy resolution, hit multiplicity, and Landau fit parameters.

<p align="center">
  <img src="hibeam-tpc-analysis/output/segmentation/krakow_overlay.png" width="48%"/>
  <img src="hibeam-tpc-analysis/output/segmentation/krakow_nhits_vs_nsec.png" width="48%"/>
</p>
<p align="center"><em>Left: Fitted Landau curves for each segmentation overlaid.  Right: Mean nHits and zero-hit fraction vs nSections.</em></p>

### Particle identification

ΔE-E telescope technique using the Proto-TPC (thin absorber) and scintillator stack (thick absorber).  Different particle species trace distinct hyperbolic bands.

<p align="center">
  <img src="hibeam-tpc-analysis/output/sim/krakow_pid.png" width="50%"/>
</p>
<p align="center"><em>2D ΔE–E distribution for Krakow proton scattering simulation.</em></p>

### 3D event displays

Reconstructed tracks in the prototype TPC, coloured by deposited charge.

<p align="center">
  <img src="hibeam-tpc-analysis/output/exp/run_0042_3d_display.png" width="48%"/>
  <img src="hibeam-tpc-analysis/output/exp/run_0042_padplane.png" width="48%"/>
</p>
<p align="center"><em>Left: 3D scatter plot of reconstructed tracks.  Right: Pad-plane hit occupancy and mean signal.</em></p>

> **Note:** The figure paths above assume you have run the full analysis pipeline.
> If figures are not yet generated, run `python scripts/run_all.py` from `hibeam-tpc-analysis/`.

---

## Quick start

### Prerequisites

**Simulation chain** (C++):
- CMake
- Geant4 (v11.0+)
- ROOT (6.14+)
- VGM (v5.2+)

**Analysis** (Python):
- Python 3.9+
- Dependencies installed via `pip install -e .` (numpy, scipy, matplotlib, uproot, awkward, mplhep, PyYAML)
- PyROOT (for experimental data only — simulation analysis uses uproot)

### Running the full chain

```bash
# 1. Build geometry
cd initial-proj/hibeam_g4_geobuilder_build/
./geobuilder --hibeam --save_as=geometry.root

# 2. Run simulation
cd ../hibeam_g4_build/
./hibeam_g4 -m macros/run_krakow.mac -c macros/config_krakow.par KrakowOutput.root

# 3. Process raw output
cd ../hibeam_g4_analysis_build/
./hibeam_ana --in=KrakowOutput.root --out=KrakowScatter.root

# 4. Run analysis (produces all figures)
cd ../../hibeam-tpc-analysis/
pip install -e .
python scripts/run_all.py
```

### Running just the analysis

If you already have processed ROOT files:

```bash
cd hibeam-tpc-analysis/

# Symlink your data
ln -s /path/to/your/simulation/files  data/simulation/
ln -s /path/to/your/experimental/files data/experimental/

# Edit paths if filenames differ
vim config.yaml

# Run everything
python scripts/run_all.py

# Or run individual steps
python scripts/run_sim_dedx.py
python scripts/run_exp_dedx.py
python scripts/run_comparison.py
python scripts/run_segmentation.py
```

### Changing analysis parameters

Everything is controlled by `config.yaml` — no Python editing needed:

```yaml
# Change the truncation fraction
fitting:
  truncation: 0.60          # was 0.70

# Switch to ATLAS plot style
plotting:
  style: atlas               # was ess

# Add a new experimental run
paths:
  experimental:
    run_0099: data/experimental/my_new_run.root

# Enable true dE/dx [MeV/mm]
geometry:
  pad_width: 3.0             # was null
  pad_height: 6.0
  drift_v: 0.05
  time_bin: 100.0
```

---

## Simulation tools reference

### hibeam_g4_geobuilder

Builds experimental geometries using ROOT TGeo classes for use with `hibeam_g4`.

```bash
./geobuilder [OPTIONS]
```

| Option | Description |
|---|---|
| `--hibeam` | Complete HIBEAM annihilation setup |
| `--sterile` | Sterile neutron annihilation setup |
| `--axion` | Axion search neutron detector |
| `--tpc` | HIBEAM TPC only |
| `--barrel` | WASA scintillator barrel only |
| `--save_as=FILE` | Save geometry to ROOT file |
| `--view` | Display in ROOT event viewer |

### hibeam_g4

Geant4 simulation of the HIBEAM detector and WASA calorimeter.

```bash
./hibeam_g4 [-g] [-t N] [-i input] [-m macro.mac] [-c config.par] output.root
```

Configuration file parameters include `Source` (gps/mcpl/scattering), `Geometry_Namefile`, `Detectors` (comma-separated sensitive volume list), and `WriteTree` (for full trajectory output).

### hibeam_g4_analysis

Converts raw Geant4 output to processed ROOT ntuples.  Auto-detects detector branches.

```bash
./hibeam_ana --in=input.root --out=output.root
```

Supported branches: Primary, TARGET, TPC, ProtoTPC, SECE (WASA), CV_bar (cosmic veto), HRDBar, Scintillator, TrigSci.

---

## Experimental data

Beam-test data was collected at the IFJ Krakow Cyclotron Centre Bronowice (CCB) with the prototype TPC:

| Run | Configuration | Description |
|---|---|---|
| 0006 | Baseline, 10×10 mm slits | CCB proton beam |
| 0042 | Baseline, 30×30 mm slits | CCB proton beam |
| 0925 | 500 V/cm, cosmics | Cosmic ray calibration |
| 0926 | 243 V/cm, cosmics | Cosmic ray calibration |

---

## References

**Energy loss and fitting:**
- L.D. Landau (1944). *J. Phys. USSR* 8, 201.
- J.E. Moyal (1955). *Phil. Mag.* 46, 263.
- H. Bichsel (1988). *Rev. Mod. Phys.* 60, 663.
- ALICE Collaboration (2010). *JINST* 5, P09002.
- Particle Data Group (2022). *PTEP* 2022, 083C01, Section 34.

**Detector techniques:**
- STAR Collaboration (2003). *NIM A* 499, 659.
- S. Carboni et al. (2012). *NIM A* 664, 251.
- S. Meroli et al. (2011). *JINST* 6, P06013.