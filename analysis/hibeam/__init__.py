"""
HIBEAM TPC Analysis Package
============================
Clean, modular analysis framework for the HIBEAM prototype TPC.

Subpackages
-----------
io       — Data loaders for simulation, experimental, and segmentation ROOT files.
physics  — Pure computation: dE/dx, Landau fitting, Bethe-Bloch, track reconstruction.
plotting — Publication-quality figures: histograms, overlays, event displays.

Quick start
-----------
    from hibeam import config
    cfg = config.load()                         # read config.yaml

    from hibeam.io import sim_loader
    data = sim_loader.load(cfg.paths.simulation.krakow)

    from hibeam.physics import dedx, fitting
    values = dedx.compute_dedx(data)
    result = fitting.fit_landau(values)

    from hibeam.plotting import style, histograms
    style.apply(cfg)
    histograms.plot_dedx(result, output="output/sim/krakow")
"""

__version__ = "1.0.0"
