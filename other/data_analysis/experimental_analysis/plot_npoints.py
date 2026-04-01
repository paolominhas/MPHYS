"""
plot_npoints.py
===============
Reads npoints_dump.csv produced by extract_npoints.C and plots:
  1. nPoints distribution across all files (stacked histogram)
  2. Per-file comparison
  3. Printed summary statistics
"""

import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import os

CSV_PATH    = "npoints_dump.csv"
OUTPUT_PATH = "analysis_plots/npoints_distribution.png"

if not os.path.exists(CSV_PATH):
    raise FileNotFoundError(
        f"'{CSV_PATH}' not found.\n"
        "Run first:  root -b -q extract_npoints.C"
    )

os.makedirs("analysis_plots", exist_ok=True)

df = pd.read_csv(CSV_PATH)
print(f"Loaded {len(df):,} track entries from {df['file'].nunique()} file(s)\n")

# ── Summary statistics ────────────────────────────────────────────────────────
print("=" * 55)
print(f"{'Statistic':<20} {'All files':>12}")
print("=" * 55)
for stat, val in [
    ("Total tracks",    f"{len(df):,}"),
    ("Mean nPoints",    f"{df['nPoints'].mean():.2f}"),
    ("Median nPoints",  f"{df['nPoints'].median():.0f}"),
    ("Std nPoints",     f"{df['nPoints'].std():.2f}"),
    ("Min nPoints",     f"{df['nPoints'].min()}"),
    ("Max nPoints",     f"{df['nPoints'].max()}"),
    ("Frac 7–10",       f"{((df['nPoints']>=7)&(df['nPoints']<=10)).mean()*100:.1f}%"),
]:
    print(f"  {stat:<18} {val:>12}")

print()
print("Per-file breakdown:")
print("-" * 55)
for fname, grp in df.groupby("file"):
    frac = ((grp['nPoints']>=7) & (grp['nPoints']<=10)).mean() * 100
    print(f"  {fname[:45]:<45}")
    print(f"    n={len(grp):,}  mean={grp['nPoints'].mean():.2f}  "
          f"median={grp['nPoints'].median():.0f}  "
          f"7–10: {frac:.1f}%")
print("=" * 55)

# ── Plot ──────────────────────────────────────────────────────────────────────
COLORS = ["#2E86AB", "#E84855", "#F4A261", "#3BB273", "#7B2D8B"]
files  = sorted(df["file"].unique())

fig = plt.figure(figsize=(14, 9))
fig.patch.set_facecolor("white")
gs  = gridspec.GridSpec(1, 2, figure=fig, wspace=0.35,
                        left=0.08, right=0.97, top=0.90, bottom=0.10)

ax_all = fig.add_subplot(gs[0])
ax_per = fig.add_subplot(gs[1])

# ── Left: combined distribution ───────────────────────────────────────────────
max_np = int(df["nPoints"].max())
bins   = np.arange(0, max_np + 2) - 0.5   # centred on integers

ax_all.hist(df["nPoints"], bins=bins, color="#2E86AB", edgecolor="white",
            linewidth=0.5, label=f"All tracks  (n={len(df):,})")

# Shade the expected 7–10 range
ax_all.axvspan(6.5, 10.5, color="#F4A261", alpha=0.20,
               label="Expected range (7–10)")
ax_all.axvline(df["nPoints"].mean(), color="#E84855", ls="--", lw=1.8,
               label=f"Mean = {df['nPoints'].mean():.2f}")
ax_all.axvline(df["nPoints"].median(), color="#3BB273", ls=":", lw=1.8,
               label=f"Median = {df['nPoints'].median():.0f}")

ax_all.set_xlabel("nPoints per track", fontsize=12)
ax_all.set_ylabel("Track count", fontsize=12)
ax_all.set_title("nPoints Distribution\n(all files combined)", fontsize=12)
ax_all.legend(fontsize=9)
ax_all.grid(axis="y", alpha=0.3)
ax_all.set_xticks(np.arange(0, max_np + 1, max(1, max_np // 15)))

# ── Right: per-file comparison ────────────────────────────────────────────────
bottom = np.zeros(max_np + 1)
x      = np.arange(0, max_np + 1)

for i, fname in enumerate(files):
    grp    = df[df["file"] == fname]
    counts = np.array([(grp["nPoints"] == v).sum() for v in x], dtype=float)
    label  = fname[:35] + "…" if len(fname) > 35 else fname
    ax_per.bar(x, counts, bottom=bottom, color=COLORS[i % len(COLORS)],
               edgecolor="white", linewidth=0.3, label=label, width=0.85)
    bottom += counts

ax_per.axvspan(6.5, 10.5, color="#F4A261", alpha=0.18, label="Expected (7–10)")
ax_per.set_xlabel("nPoints per track", fontsize=12)
ax_per.set_ylabel("Track count (stacked)", fontsize=12)
ax_per.set_title("nPoints Distribution\nper File (stacked)", fontsize=12)
ax_per.legend(fontsize=6, loc="upper right")
ax_per.grid(axis="y", alpha=0.3)
ax_per.set_xticks(np.arange(0, max_np + 1, max(1, max_np // 15)))

fig.suptitle("Track nPoints Verification  —  HIBEAM Prototype TPC",
             fontsize=13, fontweight="bold")

plt.savefig(OUTPUT_PATH, dpi=180, bbox_inches="tight")
plt.close()
print(f"\n✔  Plot saved → {OUTPUT_PATH}")