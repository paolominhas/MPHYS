#!/usr/bin/env python3
"""Overlay experimental dE/dx with simulation for each run × sim pair."""

import sys
from pathlib import Path
REPO_ROOT = Path(__file__).resolve().parent.parent
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))



from hibeam import config
from hibeam.io import exp_loader, sim_loader
from hibeam.physics import dedx
from hibeam.plotting import style, overlays


def main():
    cfg = config.load()
    style.apply(cfg)

    cuts = cfg.cuts
    fit_cfg = cfg.fitting

    # Define comparison pairs: (exp_run, sim_dataset)
    pairs = [
        ("run_0006", "krakow"),
        ("run_0042", "krakow"),
        ("run_0006", "muon"),
        ("run_0042", "muon"),
    ]

    for exp_name, sim_name in pairs:
        exp_path = config.resolve_path(cfg, cfg.paths.experimental.get(exp_name, ""))
        sim_path = config.resolve_path(cfg, cfg.paths.simulation.get(sim_name, ""))

        if not exp_path.exists() or not sim_path.exists():
            print(f"  [SKIP] {exp_name} vs {sim_name}: file not found")
            continue

        print(f"\n{'═' * 60}")
        print(f"  {exp_name} vs {sim_name}")
        print(f"{'═' * 60}")

        # Load experimental
        exp_data = exp_loader.load(
            str(exp_path),
            headers_dir=str(config.resolve_path(cfg, cfg.paths.headers_dir)),
            headers_so=str(config.resolve_path(cfg, cfg.paths.headers_so)),
            chi2_ndf_max=cuts.chi2_ndf_max,
            min_track_points=cuts.min_track_points,
        )

        # Load simulation
        sim_data = sim_loader.load_prototpc(str(sim_path))
        sim_loader.edep_to_mev(sim_data)
        ds = getattr(cfg.simulation, sim_name, {})
        low_cut = ds.get("low_mev_cut", 0.0) if isinstance(ds, dict) else 0.0
        sim_values = dedx.compute_dedx(
            sim_data, method="sum_edep",
            min_steps=cfg.simulation.min_steps,
            low_cut_mev=low_cut,
        )

        output = config.resolve_path(
            cfg, f"{cfg.paths.output_dir}/comparison/{exp_name}_vs_{sim_name}")

        overlays.overlay_data_sim(
            data_values=exp_data["dedx"],
            sim_values=sim_values,
            sim_weights=sim_data["weights"],
            data_label=f"Data ({exp_name})",
            sim_label=f"Sim ({sim_name})",
            truncation=fit_cfg.truncation,
            title=f"dE/dx comparison: {exp_name} vs {sim_name}",
            output=output,
            cfg=cfg,
        )


if __name__ == "__main__":
    main()
