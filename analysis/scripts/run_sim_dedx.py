#!/usr/bin/env python3
"""Run simulation dE/dx analysis for all configured datasets."""

import sys
from pathlib import Path
REPO_ROOT = Path(__file__).resolve().parent.parent
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))



from hibeam import config
from hibeam.io import sim_loader
from hibeam.physics import dedx, fitting
from hibeam.plotting import style, histograms
from hibeam.utils import print_fit_summary


def main():
    cfg = config.load()
    style.apply(cfg)

    sim_cfg = cfg.simulation
    fit_cfg = cfg.fitting

    for name, path in cfg.paths.simulation.items():
        resolved = config.resolve_path(cfg, path)
        if not resolved.exists():
            print(f"  [SKIP] {name}: {resolved} not found")
            continue

        print(f"\n{'═' * 60}")
        print(f"  Dataset: {name}")
        print(f"{'═' * 60}")

        # Load
        data = sim_loader.load_prototpc(str(resolved), tree_name=sim_cfg.tree_name)
        sim_loader.edep_to_mev(data)

        # Per-dataset overrides
        ds = getattr(sim_cfg, name, {})
        low_cut = ds.get("low_mev_cut", 0.0) if isinstance(ds, dict) else 0.0
        x_max   = ds.get("x_max", None) if isinstance(ds, dict) else None

        # Compute dE/dx
        values = dedx.compute_dedx(
            data, method="sum_edep",
            min_steps=sim_cfg.min_steps,
            low_cut_mev=low_cut,
        )
        print(f"  {len(values):,} events accepted")

        # Fit
        result = fitting.fit_landau(
            values,
            truncation=fit_cfg.truncation,
            n_bins=fit_cfg.n_bins,
            model=fit_cfg.model,
        )
        print_fit_summary(result.as_dict())

        # Plot
        output = config.resolve_path(cfg, f"{cfg.paths.output_dir}/sim/{name}_dedx")
        histograms.plot_dedx(
            result,
            output=output,
            title=f"Proto-TPC dE/dx — {name}",
            x_max=x_max,
            cfg=cfg,
        )


if __name__ == "__main__":
    main()
