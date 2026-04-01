# HIBEAM TPC Analysis

Clean, modular analysis framework for the HIBEAM prototype TPC detector.

## Quick start

```bash
# 1. Install dependencies
pip install -e .

# 2. Symlink your data (or adjust paths in config.yaml)
ln -s /path/to/your/root/files data/simulation/
ln -s /path/to/your/exp/files  data/experimental/

# 3. Run everything
python scripts/run_all.py
```

All figures appear in `output/`.

## How to change things

**Everything** is controlled by `config.yaml` — no Python editing needed for routine use.

| I want to...                         | Edit this in `config.yaml`           |
|--------------------------------------|--------------------------------------|
| Change truncation fraction           | `fitting: truncation: 0.60`          |
| Change plot colours                  | `plotting: series_colors:`           |
| Switch to ATLAS style                | `plotting: style: atlas`             |
| Add a new experimental run           | `paths: experimental: run_0099: ...` |
| Change the chi² quality cut          | `cuts: chi2_ndf_max: 10.0`           |
| Enable true dE/dx [MeV/mm]          | `geometry: pad_width: 3.0` (etc.)    |
| Change histogram bins                | `fitting: n_bins: 60`                |
| Change output format                 | `plotting: formats: [pdf]`           |
| Remove the ESS logo                  | `plotting: logo: show: false`        |

## Run individual analyses

```bash
python scripts/run_sim_dedx.py        # simulation dE/dx
python scripts/run_exp_dedx.py        # experimental dE/dx
python scripts/run_comparison.py      # data vs simulation overlays
python scripts/run_segmentation.py    # segmentation study
python scripts/run_all.py --only sim  # just one step
```

## Use as a library

```python
from hibeam import config
from hibeam.io import sim_loader
from hibeam.physics import dedx, fitting
from hibeam.plotting import style, histograms

cfg = config.load()
style.apply(cfg)

data   = sim_loader.load_prototpc("data/simulation/KrakowScatter.root")
sim_loader.edep_to_mev(data)
values = dedx.compute_dedx(data, min_steps=5)
result = fitting.fit_landau(values, truncation=0.70)
histograms.plot_dedx(result, output="output/my_plot")
```

## Architecture

```
hibeam/
├── io/          → Data loaders (sim, exp, segmentation, CSV)
├── physics/     → Pure computation (dE/dx, fitting, Bethe-Bloch, tracks, PID)
└── plotting/    → All figures (style, histograms, overlays, displays, segmentation)
```

Each layer depends only on the one above it.  `physics/` never imports `matplotlib`.
`io/` never imports from `physics/`.  This makes the code easy to test and reason about.

## Dependencies

- **uproot** + **awkward**: Read simulation ROOT files (no PyROOT needed).
- **PyROOT**: Read experimental data (requires compiled TrackData headers).
- **scipy**: Landau (Moyal) fitting.
- **matplotlib** + **mplhep**: Publication-quality figures.
- **pyyaml**: Config file parsing.
