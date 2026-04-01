#!/usr/bin/env python3
"""Run experimental dE/dx analysis for all configured runs."""

from hibeam import config
from hibeam.io import exp_loader
from hibeam.physics import fitting
from hibeam.plotting import style, histograms
from hibeam.utils import print_fit_summary


def main():
    cfg = config.load()
    style.apply(cfg)

    cuts = cfg.cuts
    fit_cfg = cfg.fitting

    for name, path in cfg.paths.experimental.items():
        resolved = config.resolve_path(cfg, path)
        if not resolved.exists():
            print(f"  [SKIP] {name}: {resolved} not found")
            continue

        print(f"\n{'═' * 60}")
        print(f"  Run: {name}")
        print(f"{'═' * 60}")

        data = exp_loader.load(
            str(resolved),
            headers_dir=str(config.resolve_path(cfg, cfg.paths.headers_dir)),
            headers_so=str(config.resolve_path(cfg, cfg.paths.headers_so)),
            chi2_ndf_max=cuts.chi2_ndf_max,
            min_track_points=cuts.min_track_points,
            min_adc=cuts.min_adc,
            upper_percentile=cfg.experimental.upper_percentile,
        )

        values = data["dedx"]
        if len(values) < fit_cfg.min_events:
            print(f"  [SKIP] Only {len(values)} values — need {fit_cfg.min_events}")
            continue

        result = fitting.fit_landau(
            values,
            truncation=fit_cfg.truncation,
            n_bins=fit_cfg.n_bins,
            model=fit_cfg.model,
        )
        print_fit_summary(result.as_dict())

        output = config.resolve_path(cfg, f"{cfg.paths.output_dir}/exp/{name}_dedx")
        histograms.plot_dedx(
            result,
            output=output,
            title=f"Experimental dE/dx — {name}",
            cfg=cfg,
        )


if __name__ == "__main__":
    main()
