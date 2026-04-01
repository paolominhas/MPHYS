#!/usr/bin/env python3
"""Run the full segmentation study for all configured directories."""

from hibeam import config
from hibeam.io import seg_loader
from hibeam.plotting import style, segmentation_plots


def main():
    cfg = config.load()
    style.apply(cfg)

    seg = cfg.segmentation
    outdir = config.resolve_path(cfg, cfg.paths.output_dir) / "segmentation"

    for label, dir_key in [("krakow", "krakow_dir"), ("muon", "muon_dir")]:
        seg_dir = config.resolve_path(cfg, cfg.paths.segmentation.get(dir_key, ""))
        if not seg_dir.is_dir():
            print(f"  [SKIP] {label}: {seg_dir} not found")
            continue

        print(f"\n{'═' * 60}")
        print(f"  Segmentation study: {label}")
        print(f"{'═' * 60}")

        records = seg_loader.load(
            str(seg_dir),
            pattern=seg.pattern,
            tree_name=seg.tree_name,
            edep_branch=seg.branch,
            min_steps=seg.min_steps,
            edep_floor_mev=seg.edep_floor_mev,
        )

        # Landau overlay
        segmentation_plots.overlay_segmentations(
            records,
            truncation=cfg.fitting.truncation,
            n_bins=seg.fit_nbins,
            edep_floor_mev=seg.edep_floor_mev,
            title=f"Segmentation overlay — {label}",
            output=outdir / f"{label}_overlay",
            cfg=cfg,
        )

        # nHits vs nSections
        segmentation_plots.plot_nhits_vs_nsec(
            records,
            title=f"nHits vs segmentation — {label}",
            output=outdir / f"{label}_nhits_vs_nsec",
            cfg=cfg,
        )

        # nHits distribution
        segmentation_plots.plot_nhits_distribution(
            records,
            title=f"nHits distribution — {label}",
            output=outdir / f"{label}_nhits_dist",
            cfg=cfg,
        )


if __name__ == "__main__":
    main()
