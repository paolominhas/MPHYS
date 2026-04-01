#!/usr/bin/env python3
"""Generate ΔE-E particle identification plots for all simulation datasets."""

from hibeam import config
from hibeam.io import sim_loader
from hibeam.physics import pid
from hibeam.plotting import style, displays


def main():
    cfg = config.load()
    style.apply(cfg)

    pid_cfg = cfg.pid
    outdir = config.resolve_path(cfg, cfg.paths.output_dir) / "sim"

    for name, path in cfg.paths.simulation.items():
        resolved = config.resolve_path(cfg, path)
        if not resolved.exists():
            print(f"  [SKIP] {name}: {resolved} not found")
            continue

        print(f"\n{'═' * 60}")
        print(f"  PID plot: {name}")
        print(f"{'═' * 60}")

        try:
            data = sim_loader.load_pid(str(resolved))

            obs = pid.compute_pid_observables(
                data["tpc_mev"], data["cal_mev"],
                data["n_tpc"], data["n_cal"],
                min_tpc_hits=pid_cfg.min_tpc_hits,
                min_scint_hits=pid_cfg.min_scint_hits,
            )

            print(f"  Selected events: {obs['n_selected']:,}")

            if obs["n_selected"] < 10:
                print("  [SKIP] Too few events for PID plot")
                continue

            displays.pid_plot(
                obs["delta_e"], obs["e_residual"],
                title=rf"$\Delta E$–$E$ PID — {name}",
                output=outdir / f"{name}_pid",
                cfg=cfg,
            )

        except Exception as e:
            print(f"  [SKIP] {e}")


if __name__ == "__main__":
    main()
