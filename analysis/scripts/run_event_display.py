#!/usr/bin/env python3
"""Generate 3D event displays and pad-plane heatmaps for all experimental runs."""

import sys
from pathlib import Path
REPO_ROOT = Path(__file__).resolve().parent.parent
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))



from hibeam import config
from hibeam.io import exp_loader
from hibeam.plotting import style, displays


def main():
    cfg = config.load()
    style.apply(cfg)

    for name, path in cfg.paths.experimental.items():
        resolved = config.resolve_path(cfg, path)
        if not resolved.exists():
            print(f"  [SKIP] {name}: {resolved} not found")
            continue

        print(f"\n{'═' * 60}")
        print(f"  Event display: {name}")
        print(f"{'═' * 60}")

        outdir = config.resolve_path(cfg, cfg.paths.output_dir) / "exp"

        # ── 3D track display (needs PyROOT headers) ─────────────────────
        try:
            data = exp_loader.load(
                str(resolved),
                headers_dir=str(config.resolve_path(cfg, cfg.paths.headers_dir)),
                headers_so=str(config.resolve_path(cfg, cfg.paths.headers_so)),
                chi2_ndf_max=cfg.cuts.chi2_ndf_max,
                min_track_points=cfg.cuts.min_track_points,
            )

            if data["tracks"]:
                displays.event_display_3d(
                    data["tracks"],
                    max_tracks=100,
                    title=f"3D track display — {name}",
                    output=outdir / f"{name}_3d_display",
                    cfg=cfg,
                )
        except Exception as e:
            print(f"  [SKIP] Track display: {e}")

        # ── Pad-plane heatmap (uproot, no headers needed) ───────────────
        try:
            raw = exp_loader.load_raw_hits(
                str(resolved),
                noise_sigma=cfg.cuts.noise_sigma,
                min_adc_signal=cfg.cuts.min_adc_signal,
                entry_stop=5000,
            )

            displays.pad_plane_heatmap(
                raw["row"], raw["col"], raw["signal"],
                title=f"Pad-plane — {name}",
                output=outdir / f"{name}_padplane",
                cfg=cfg,
            )
        except Exception as e:
            print(f"  [SKIP] Pad-plane: {e}")


if __name__ == "__main__":
    main()
