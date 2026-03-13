"""
HIBEAM Prototype TPC — Pad-Plane Diagnostic Analysis
=====================================================
Uses tpc branch (fully split, uproot-readable):
  tpc.val      — raw ADC
  tpc.pedestal — per-hit electronic baseline
  tpc.peddev   — pedestal RMS (used for noise threshold)
  tpc.row      — pad row  (y-axis of pad plane)
  tpc.column   — pad column (x-axis of pad plane)

For each file produces a 4-panel figure:
  1. 2D heatmap: mean pedestal-subtracted ADC per (row, column) pad
  2. 2D heatmap: hit occupancy (how often each pad fires)
  3. Per-row profile: mean signal ADC, hit probability, flagged outliers
  4. Mean ADC/hit vs N-hits scatter (hexbin) — breaks the tautology
"""

import uproot
import awkward as ak
import numpy as np
import glob, os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.colors import LogNorm
import warnings
warnings.filterwarnings("ignore")

# ─── Configuration ────────────────────────────────────────────────────────────
INPUT_DIR      = "../experimental_data/"
OUTPUT_DIR     = "analysis_plots/"
NOISE_SIGMA    = 3.0    # reject hits where (val-ped) < NOISE_SIGMA * peddev
MIN_TRACK_HITS = 5      # events with fewer total hits are likely noise events
# ──────────────────────────────────────────────────────────────────────────────


def load_tpc(path):
    """
    Load tpc branch, subtract pedestal, apply noise threshold.
    Returns flat arrays of signal ADC, row, column, and per-event counts.
    """
    with uproot.open(path) as f:
        tree = f["trackingData"]
        val  = tree["tpc/tpc.val"].array()
        ped  = tree["tpc/tpc.pedestal"].array()
        dev  = tree["tpc/tpc.peddev"].array()
        row  = tree["tpc/tpc.row"].array()
        col  = tree["tpc/tpc.column"].array()

    # Pedestal-subtracted signal
    sig = val - ped

    # Noise mask: keep only hits > NOISE_SIGMA * peddev above pedestal
    noise_floor = NOISE_SIGMA * dev
    good = sig > noise_floor

    sig_good = sig[good]
    row_good = row[good]
    col_good = col[good]

    # Per-event quantities (for scatter panel)
    n_hits_per_event   = ak.to_numpy(ak.sum(good,     axis=1)).astype(int)
    adc_sum_per_event  = ak.to_numpy(ak.sum(sig_good if False   # placeholder
                                            else sig * good, axis=1)).astype(float)

    # Flatten to numpy
    sig_flat = ak.to_numpy(ak.flatten(sig_good))
    row_flat = ak.to_numpy(ak.flatten(row_good)).astype(int)
    col_flat = ak.to_numpy(ak.flatten(col_good)).astype(int)

    # Event-level: only keep events with enough hits
    event_mask = n_hits_per_event >= MIN_TRACK_HITS
    n_hits_evt = n_hits_per_event[event_mask]
    adc_sum_evt = adc_sum_per_event[event_mask]
    mean_adc_per_hit = np.where(n_hits_evt > 0,
                                adc_sum_evt / n_hits_evt, 0.0)

    n_total  = len(ak.flatten(val))
    n_signal = len(sig_flat)
    print(f"  Total raw hits    : {n_total:,}")
    print(f"  After noise cut   : {n_signal:,}  ({100*n_signal/n_total:.1f}%)")
    print(f"  Valid events      : {event_mask.sum():,}  "
          f"(>= {MIN_TRACK_HITS} hits/event)")

    return (sig_flat, row_flat, col_flat,
            n_hits_evt, mean_adc_per_hit)


def make_padplane_maps(sig, rows, cols):
    """Build 2D arrays (row × col) of mean ADC and occupancy."""
    if len(rows) == 0:
        return None, None, None, None

    row_min, row_max = rows.min(), rows.max()
    col_min, col_max = cols.min(), cols.max()
    n_rows = row_max - row_min + 1
    n_cols = col_max - col_min + 1

    adc_map  = np.zeros((n_rows, n_cols), dtype=float)
    cnt_map  = np.zeros((n_rows, n_cols), dtype=float)

    r_idx = rows - row_min
    c_idx = cols - col_min

    np.add.at(adc_map, (r_idx, c_idx), sig)
    np.add.at(cnt_map, (r_idx, c_idx), 1.0)

    mean_map = np.where(cnt_map > 0, adc_map / cnt_map, np.nan)

    row_labels = np.arange(row_min, row_max + 1)
    col_labels = np.arange(col_min, col_max + 1)
    return mean_map, cnt_map, row_labels, col_labels


def per_row_stats(sig, rows, n_hits_evt):
    """
    Per-row mean signal ADC and hit probability.
    Hit probability = fraction of valid events in which this row fired at all.
    """
    unique_rows = np.array(sorted(set(rows)), dtype=int)
    mean_adc  = np.array([sig[rows == r].mean() if np.sum(rows==r)>0
                          else np.nan for r in unique_rows])
    hit_count = np.array([np.sum(rows == r) for r in unique_rows], dtype=float)

    # Hit probability: we don't have per-event row info in flat arrays,
    # so use relative occupancy (hit_count / total hits) as a proxy
    hit_prob = hit_count / hit_count.sum()

    return unique_rows, mean_adc, hit_prob, hit_count


def analyse_file(path, output_dir):
    filename = os.path.basename(path)
    print(f"\n{'─'*64}\n  {filename}")

    sig, rows, cols, n_hits_evt, mean_adc_per_hit = load_tpc(path)
    if len(sig) < 100:
        print("  ⚠  Too few hits after cuts, skipping.")
        return

    mean_map, cnt_map, row_labels, col_labels = make_padplane_maps(sig, rows, cols)
    unique_rows, row_mean_adc, row_hit_prob, row_hit_count = \
        per_row_stats(sig, rows, n_hits_evt)

    # ── Figure ────────────────────────────────────────────────────────────
    fig = plt.figure(figsize=(16, 14))
    fig.patch.set_facecolor("white")
    gs = gridspec.GridSpec(2, 2, figure=fig,
                           hspace=0.40, wspace=0.35,
                           left=0.08, right=0.97, top=0.93, bottom=0.07)

    ax_adc   = fig.add_subplot(gs[0, 0])   # mean ADC heatmap
    ax_occ   = fig.add_subplot(gs[0, 1])   # occupancy heatmap
    ax_row   = fig.add_subplot(gs[1, 0])   # per-row bar chart
    ax_hex   = fig.add_subplot(gs[1, 1])   # hexbin scatter

    # ── Panel 1: Mean pedestal-subtracted ADC per pad ─────────────────────
    if mean_map is not None:
        extent = [col_labels[0]-0.5, col_labels[-1]+0.5,
                  row_labels[0]-0.5, row_labels[-1]+0.5]
        im1 = ax_adc.imshow(mean_map, aspect="auto", origin="lower",
                             extent=extent, cmap="plasma",
                             interpolation="nearest")
        plt.colorbar(im1, ax=ax_adc, label="Mean Signal ADC [ADC - pedestal]")

        # Overlay dead/cold pads: below 20% of global median
        median_adc = np.nanmedian(mean_map)
        cold_mask  = mean_map < 0.2 * median_adc
        r_cold, c_cold = np.where(cold_mask)
        if r_cold.size:
            ax_adc.scatter(col_labels[c_cold], row_labels[r_cold],
                           marker="x", color="cyan", s=20, lw=1,
                           label=f"Cold pads (<20% median, n={r_cold.size})")
            ax_adc.legend(fontsize=7, loc="upper right")

    ax_adc.set_xlabel("Column (pad strip)", fontsize=10)
    ax_adc.set_ylabel("Row (TPC layer)", fontsize=10)
    ax_adc.set_title("Mean Pedestal-Subtracted ADC\nper Pad Cell", fontsize=11)

    # ── Panel 2: Hit occupancy (log scale) ───────────────────────────────
    if cnt_map is not None:
        cnt_safe = np.where(cnt_map > 0, cnt_map, np.nan)
        im2 = ax_occ.imshow(cnt_safe, aspect="auto", origin="lower",
                             extent=extent,
                             norm=LogNorm(vmin=1,
                                          vmax=float(np.nanmax(cnt_safe))),
                             cmap="viridis", interpolation="nearest")
        plt.colorbar(im2, ax=ax_occ, label="Hit Count (log scale)")

        # Highlight dead pads (zero hits)
        dead_r, dead_c = np.where(cnt_map == 0)
        if dead_r.size:
            ax_occ.scatter(col_labels[dead_c], row_labels[dead_r],
                           marker="x", color="red", s=20, lw=1,
                           label=f"Dead pads (0 hits, n={dead_r.size})")
            ax_occ.legend(fontsize=7, loc="upper right")

    ax_occ.set_xlabel("Column (pad strip)", fontsize=10)
    ax_occ.set_ylabel("Row (TPC layer)", fontsize=10)
    ax_occ.set_title("Hit Occupancy per Pad Cell\n(log scale)", fontsize=11)

    # ── Panel 3: Per-row mean ADC with outlier flags ──────────────────────
    global_mean = np.nanmean(row_mean_adc)
    global_std  = np.nanstd(row_mean_adc)
    outlier     = np.abs(row_mean_adc - global_mean) > 2 * global_std
    low_occ     = row_hit_count < 0.5 * np.median(row_hit_count)

    bar_colors = ["#E84855" if (o or l) else "#2E86AB"
                  for o, l in zip(outlier, low_occ)]
    ax_row.bar(unique_rows, row_mean_adc, color=bar_colors,
               edgecolor="none", width=0.85)
    ax_row.axhline(global_mean, color="#3BB273", ls="--", lw=1.8,
                   label=f"Mean = {global_mean:.1f}")
    ax_row.axhspan(global_mean - 2*global_std,
                   global_mean + 2*global_std,
                   color="#3BB273", alpha=0.12, label="±2σ band")
    ax_row.set_xlabel("Row (TPC layer)", fontsize=10)
    ax_row.set_ylabel("Mean Signal ADC", fontsize=10)
    ax_row.set_title("Per-Row Mean Signal ADC\n"
                     "(red = outlier >2σ or <50% occupancy)", fontsize=11)
    ax_row.legend(fontsize=8)
    ax_row.grid(axis="y", alpha=0.3)

    # Print flagged rows
    bad_rows = unique_rows[outlier | low_occ]
    if len(bad_rows):
        print(f"  ⚠  Flagged rows: {list(bad_rows)}")

    # ── Panel 4: Hexbin — mean ADC/hit vs N hits per event ───────────────
    # This is the scientifically meaningful correlation:
    # high ionisation (large mean ADC/hit) should → fewer but stronger hits
    if len(n_hits_evt) > 10:
        hb = ax_hex.hexbin(n_hits_evt, mean_adc_per_hit,
                           gridsize=50, cmap="inferno",
                           mincnt=1, bins="log")
        plt.colorbar(hb, ax=ax_hex, label="log₁₀(event count)")

        # Profile
        edges = np.unique(np.percentile(n_hits_evt,
                                        np.linspace(2, 98, 18)))
        px, py, pe = [], [], []
        for lo, hi in zip(edges[:-1], edges[1:]):
            sel = mean_adc_per_hit[(n_hits_evt >= lo) & (n_hits_evt < hi)]
            if len(sel) < 5:
                continue
            px.append(0.5*(lo + hi))
            py.append(np.median(sel))
            pe.append(sel.std() / np.sqrt(len(sel)))
        ax_hex.errorbar(px, py, yerr=pe, fmt="o-", color="white",
                        lw=1.8, ms=4, zorder=5, label="Median profile")
        ax_hex.legend(fontsize=8)

    ax_hex.set_xlabel("N hits per event (above noise threshold)", fontsize=10)
    ax_hex.set_ylabel("Mean ADC per hit  [pedestal-subtracted]", fontsize=10)
    ax_hex.set_title("Mean Signal Strength vs. Track Multiplicity\n"
                     "(dense = typical; sparse = rare/large-angle)", fontsize=11)

    # ── Suptitle ──────────────────────────────────────────────────────────
    fig.suptitle(
        f"HIBEAM TPC Pad-Plane Diagnostics — {filename}\n"
        f"Signal hits: {len(sig):,}  |  "
        f"Noise threshold: {NOISE_SIGMA}σ  |  "
        f"Events: {len(n_hits_evt):,}",
        fontsize=12, fontweight="bold", y=0.97)

    out = os.path.join(output_dir,
                       filename.replace(".root", "_padplane.png"))
    plt.savefig(out, dpi=180, bbox_inches="tight")
    plt.close(fig)
    print(f"  ✔  Saved → {out}")


def run_batch(input_dir, output_dir):
    os.makedirs(output_dir, exist_ok=True)
    files = sorted(glob.glob(os.path.join(input_dir, "*.root")))
    if not files:
        print(f"No .root files in {input_dir}")
        return
    print(f"Found {len(files)} file(s)\n")
    for f in files:
        try:
            analyse_file(f, output_dir)
        except Exception as e:
            print(f"  ✗  {type(e).__name__}: {e}")
    print(f"\n{'═'*64}\n  Done.")


if __name__ == "__main__":
    run_batch(INPUT_DIR, OUTPUT_DIR)