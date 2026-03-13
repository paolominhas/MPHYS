"""
npoints_vs_energy_cut.py  —  publication-quality figure
=========================================================
ESS / HIBEAM branding with clean scientific typography.

Error treatment:
  x  — quadrature peddev sum → ADC measurement uncertainty
  y  — SEM on mean nPoints
  survival — binomial  √(p(1−p)/N)
  histograms — Poisson  √N / N_total

Place your ESS logo PNG at:  data/ess_logo.png
(Any aspect ratio works; it will be letterboxed automatically.)

Run after:  root -b -q extract_npoints.C
"""

import uproot
import awkward as ak
import numpy as np
import pandas as pd
import os, glob
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.ticker as ticker
from matplotlib.offsetbox import OffsetImage, AnnotationBbox
import matplotlib.image as mpimg
from matplotlib.patches import FancyBboxPatch
import warnings
warnings.filterwarnings("ignore")

# ─── ESS Brand palette ────────────────────────────────────────────────────────
ESS_NAVY    = "#00263A"   # primary dark
ESS_BLUE    = "#005B8E"   # primary mid
ESS_CYAN    = "#00A9CE"   # primary bright
ESS_SILVER  = "#A8B2B8"
ESS_WHITE   = "#FFFFFF"
ESS_ORANGE  = "#E5720F"   # accent / warning
ESS_GREEN   = "#4CAF82"   # accent / good

# Data series colours — distinguishable at print and greyscale
SERIES_COLORS = [ESS_CYAN, ESS_ORANGE, ESS_GREEN, "#C45AB3", "#FFD166"]

# ─── Global matplotlib style ──────────────────────────────────────────────────
plt.rcParams.update({
    # Font
    "font.family"        : "serif",
    "font.serif"         : ["DejaVu Serif", "Georgia", "Times New Roman"],
    "mathtext.fontset"   : "dejavuserif",
    "font.size"          : 10,
    "axes.titlesize"     : 11,
    "axes.labelsize"     : 10,
    "xtick.labelsize"    : 8.5,
    "ytick.labelsize"    : 8.5,
    "legend.fontsize"    : 7.5,
    # Lines & markers
    "lines.linewidth"    : 1.6,
    "lines.markersize"   : 4.5,
    "errorbar.capsize"   : 2.5,
    # Axes
    "axes.spines.top"    : False,
    "axes.spines.right"  : False,
    "axes.linewidth"     : 0.8,
    "axes.edgecolor"     : ESS_NAVY,
    "axes.labelcolor"    : ESS_NAVY,
    "axes.titlecolor"    : ESS_NAVY,
    "axes.titlepad"      : 8,
    "axes.titleweight"   : "semibold",
    "axes.facecolor"     : "#F7F9FB",
    # Grid
    "axes.grid"          : True,
    "grid.color"         : "#D0D8E0",
    "grid.linewidth"     : 0.5,
    "grid.linestyle"     : "--",
    "grid.alpha"         : 0.7,
    # Ticks
    "xtick.color"        : ESS_NAVY,
    "ytick.color"        : ESS_NAVY,
    "xtick.direction"    : "in",
    "ytick.direction"    : "in",
    "xtick.minor.visible": True,
    "ytick.minor.visible": True,
    "xtick.major.width"  : 0.8,
    "ytick.major.width"  : 0.8,
    # Figure
    "figure.facecolor"   : ESS_WHITE,
    "figure.dpi"         : 150,
    "savefig.dpi"        : 200,
    "savefig.bbox"       : "tight",
    "savefig.facecolor"  : ESS_WHITE,
    # Legend
    "legend.framealpha"  : 0.92,
    "legend.edgecolor"   : "#C8D4DC",
    "legend.fancybox"    : False,
})

# ─── Config ───────────────────────────────────────────────────────────────────
CSV_PATH      = "npoints_dump.csv"
INPUT_DIR     = "../experimental_data/"
OUTPUT_PATH   = "analysis_plots/npoints_vs_energy_cut.png"
LOGO_PATH     = "data/ess_logo.png"
NOISE_SIGMA   = 3.0
N_THRESHOLDS  = 40
SNAPSHOT_PCTS = [0, 25, 50, 75, 90]
# ──────────────────────────────────────────────────────────────────────────────

os.makedirs("analysis_plots", exist_ok=True)

if not os.path.exists(CSV_PATH):
    raise FileNotFoundError(
        f"'{CSV_PATH}' not found — run:  root -b -q extract_npoints.C  first")

df_np = pd.read_csv(CSV_PATH)
print(f"Loaded {len(df_np):,} track rows  |  {df_np['file'].nunique()} file(s)\n")

FILE_MAP = {
    os.path.basename(p).replace(".root", ""): p
    for p in sorted(glob.glob(os.path.join(INPUT_DIR, "*.root")))
}


# ─── Data loading ─────────────────────────────────────────────────────────────

def load_event_quantities(path):
    with uproot.open(path) as f:
        tree  = f["trackingData"]
        val   = tree["tpc/tpc.val"].array()
        ped   = tree["tpc/tpc.pedestal"].array()
        dev   = tree["tpc/tpc.peddev"].array()

    sig      = val - ped
    good     = sig > (NOISE_SIGMA * dev)
    adc_evt  = ak.to_numpy(ak.sum(sig  * good, axis=1)).astype(float)
    sig_evt  = ak.to_numpy(np.sqrt(ak.sum((dev * good)**2, axis=1))).astype(float)
    return adc_evt, sig_evt


# ─── Per-file sweep ───────────────────────────────────────────────────────────

results = {}

for label, path in FILE_MAP.items():
    if label not in df_np["file"].values:
        continue
    print(f"Processing {label}…")
    try:
        adc_evt, sigma_evt = load_event_quantities(path)
    except Exception as e:
        print(f"  ✗ {e}"); continue

    df_f = df_np[df_np["file"] == label].copy()
    df_f["event_adc"]   = df_f["event"].apply(
        lambda e: float(adc_evt[e])   if e < len(adc_evt)   else np.nan)
    df_f["event_sigma"] = df_f["event"].apply(
        lambda e: float(sigma_evt[e]) if e < len(sigma_evt) else np.nan)
    df_f = df_f.dropna(subset=["event_adc"])

    adc_v   = df_f["event_adc"].values
    sig_v   = df_f["event_sigma"].values
    npts_v  = df_f["nPoints"].values
    total   = len(df_f)
    thrs    = np.linspace(0, np.percentile(adc_v, 95), N_THRESHOLDS)

    mean_n = np.full(N_THRESHOLDS, np.nan)
    sem_n  = np.full(N_THRESHOLDS, np.nan)
    surv   = np.full(N_THRESHOLDS, np.nan)
    surv_e = np.full(N_THRESHOLDS, np.nan)
    x_e    = np.full(N_THRESHOLDS, np.nan)

    for k, thr in enumerate(thrs):
        m = adc_v >= thr
        n = m.sum()
        if n < 2:
            continue
        mean_n[k] = npts_v[m].mean()
        sem_n[k]  = npts_v[m].std(ddof=1) / np.sqrt(n)
        p         = n / total
        surv[k]   = p * 100
        surv_e[k] = np.sqrt(p * (1 - p) / total) * 100
        x_e[k]    = sig_v[m].mean()

    snapshots = {}
    for pct in SNAPSHOT_PCTS:
        tv = np.percentile(adc_v, pct)
        snapshots[pct] = {"npts": npts_v[adc_v >= tv], "thr_val": tv}

    results[label] = dict(
        thresholds=thrs, mean_n=mean_n, sem_n=sem_n,
        surv=surv, surv_e=surv_e, x_e=x_e,
        snapshots=snapshots, npts_all=npts_v, adc_v=adc_v,
    )
    print(f"  nPoints {npts_v.min()}–{npts_v.max()}  "
          f"mean = {npts_v.mean():.2f} ± {npts_v.std(ddof=1)/np.sqrt(total):.3f}")
    print(f"  Median ADC σ per event = {np.median(sigma_evt):.1f} counts")


# ─── Figure layout ────────────────────────────────────────────────────────────

fig = plt.figure(figsize=(14, 12))
fig.patch.set_facecolor(ESS_WHITE)

# Outer margin for header bar
gs_outer = gridspec.GridSpec(2, 1, figure=fig,
                             height_ratios=[0.06, 1],
                             hspace=0.0,
                             left=0.0, right=1.0, top=1.0, bottom=0.0)

# Header bar axis (branding strip)
ax_header = fig.add_subplot(gs_outer[0])
ax_header.set_facecolor(ESS_NAVY)
for spine in ax_header.spines.values():
    spine.set_visible(False)
ax_header.set_xticks([]); ax_header.set_yticks([])

ax_header.text(0.50, 0.5,
               "HIBEAM Prototype TPC  ·  Track Length vs. Energy Cut",
               transform=ax_header.transAxes,
               ha="center", va="center",
               color=ESS_WHITE, fontsize=12, fontweight="bold",
               fontfamily="serif")

# Sub-label right side
ax_header.text(0.98, 0.5,
               "ESS / IFJ PAN Kraków",
               transform=ax_header.transAxes,
               ha="right", va="center",
               color=ESS_CYAN, fontsize=8, fontstyle="italic")

# Plot grid
gs_plots = gridspec.GridSpecFromSubplotSpec(
    2, 2, subplot_spec=gs_outer[1],
    hspace=0.42, wspace=0.30)

fig.subplots_adjust(left=0.09, right=0.97, bottom=0.08)

ax_mean = fig.add_subplot(gs_plots[0, 0])
ax_surv = fig.add_subplot(gs_plots[0, 1])
ax_snap = fig.add_subplot(gs_plots[1, 0])
ax_all  = fig.add_subplot(gs_plots[1, 1])

first_label = next(iter(results))


# ─── Panel helper: thin coloured top-border ──────────────────────────────────

def add_panel_accent(ax, color=ESS_CYAN, lw=3):
    """Draw a thin coloured line along the top of an axes."""
    ax.plot([0, 1], [1, 1], transform=ax.transAxes,
            color=color, lw=lw, clip_on=False, zorder=10,
            solid_capstyle="butt")


# ─── Plot panels ──────────────────────────────────────────────────────────────

for i, (label, res) in enumerate(results.items()):
    c     = SERIES_COLORS[i % len(SERIES_COLORS)]
    thr   = res["thresholds"]
    mn    = res["mean_n"]
    se    = res["sem_n"]
    sf    = res["surv"]
    sfe   = res["surv_e"]
    xe    = res["x_e"]
    short = (label[:26] + "…") if len(label) > 26 else label
    valid = ~np.isnan(mn)
    vs    = ~np.isnan(sf)

    # Panel 1 — mean nPoints ± SEM with x peddev error
    ax_mean.errorbar(thr[valid], mn[valid],
                     yerr=se[valid], xerr=xe[valid],
                     fmt="o", color=c, lw=1.5, ms=3.5,
                     capsize=2.5, capthick=1.0, elinewidth=0.8,
                     label=short, zorder=4)
    ax_mean.plot(thr[valid], mn[valid], "-", color=c, lw=1.3, alpha=0.6, zorder=3)

    # Panel 2 — survival fraction
    ax_surv.errorbar(thr[vs], sf[vs],
                     yerr=sfe[vs],
                     fmt="o", color=c, lw=1.5, ms=3.5,
                     capsize=2.5, capthick=1.0, elinewidth=0.8,
                     label=short, zorder=4)
    ax_surv.plot(thr[vs], sf[vs], "-", color=c, lw=1.3, alpha=0.6, zorder=3)

    # Panel 3 — snapshot histograms (first file)
    if label == first_label:
        snap_cm = [ESS_NAVY, ESS_BLUE, ESS_CYAN, ESS_GREEN, ESS_ORANGE]
        max_np  = int(res["npts_all"].max())
        bins    = np.arange(0, max_np + 2) - 0.5
        bin_c   = (bins[:-1] + bins[1:]) / 2
        for j, pct in enumerate(SNAPSHOT_PCTS):
            snap = res["snapshots"][pct]
            npts = snap["npts"]
            n    = len(npts)
            if n == 0: continue
            hc, _ = np.histogram(npts, bins=bins)
            hn    = hc / n
            he    = np.sqrt(hc) / n
            sc    = snap_cm[j]
            ax_snap.step(bins[:-1], hn, where="post", color=sc, lw=1.8,
                         label=f"$\\geq${snap['thr_val']:.0f} ADC  "
                               f"({pct}th pct, $n$={n:,})", zorder=4)
            ax_snap.errorbar(bin_c, hn, yerr=he, fmt="none",
                             color=sc, elinewidth=0.7, capsize=1.5, zorder=5)

    # Panel 4 — overall nPoints
    max_np  = int(res["npts_all"].max())
    bins    = np.arange(0, max_np + 2) - 0.5
    bin_c   = (bins[:-1] + bins[1:]) / 2
    n       = len(res["npts_all"])
    hc, _   = np.histogram(res["npts_all"], bins=bins)
    hn      = hc / n
    he      = np.sqrt(hc) / n
    ax_all.step(bins[:-1], hn, where="post", color=c, lw=1.8,
                label=short, zorder=4)
    ax_all.errorbar(bin_c, hn, yerr=he, fmt="none",
                    color=c, elinewidth=0.7, capsize=1.5, zorder=5)


# ─── Expected-range shading & reference lines ─────────────────────────────────

for ax in (ax_snap, ax_all):
    ax.axvspan(6.5, 10.5, color=ESS_CYAN, alpha=0.10, zorder=1,
               label="Expected 7–10 pts")
    ax.axvline(7,  color=ESS_CYAN, lw=0.7, ls=":", alpha=0.6, zorder=2)
    ax.axvline(10, color=ESS_CYAN, lw=0.7, ls=":", alpha=0.6, zorder=2)

ax_mean.axhspan(7, 10, color=ESS_CYAN, alpha=0.10, zorder=1,
                label="Expected 7–10 pts")
ax_mean.axhline(7,  color=ESS_CYAN, lw=0.7, ls=":", alpha=0.6)
ax_mean.axhline(10, color=ESS_CYAN, lw=0.7, ls=":", alpha=0.6)

ax_surv.axhline(50, color=ESS_SILVER, ls="--", lw=1.2,
                label="50% level", zorder=2)


# ─── Axis labels & titles ─────────────────────────────────────────────────────

ax_mean.set_xlabel(
    r"ADC energy threshold  [ped.-subtracted counts]" "\n"
    r"$x$-error $= \langle\sigma_\mathrm{ADC}\rangle$  "
    r"from quadrature $\sum\sigma_\mathrm{ped}$",
    fontsize=8.5)
ax_mean.set_ylabel(r"Mean $n_\mathrm{points}$ per track", fontsize=10)
ax_mean.set_title(r"Track Length vs. Energy Cut"
                  "\n"
                  r"$y$-error $= \mathrm{SEM} = \sigma/\sqrt{n}$", fontsize=10)
ax_mean.legend(loc="upper right")
add_panel_accent(ax_mean, ESS_CYAN)

ax_surv.set_xlabel("ADC energy threshold  [ped.-subtracted counts]", fontsize=9)
ax_surv.set_ylabel("Track survival  [%]", fontsize=10)
ax_surv.set_title("Survival Fraction vs. Energy Cut\n"
                  r"$y$-error $= \sqrt{p(1-p)/N}$  (binomial)", fontsize=10)
ax_surv.set_ylim(-2, 107)
ax_surv.yaxis.set_major_formatter(ticker.PercentFormatter())
ax_surv.legend(loc="upper right")
add_panel_accent(ax_surv, ESS_ORANGE)

ax_snap.set_xlabel(r"$n_\mathrm{points}$ per track", fontsize=10)
ax_snap.set_ylabel("Normalised count", fontsize=10)
ax_snap.set_title("$n_\\mathrm{points}$ at Selected Thresholds\n"
                  r"$y$-error $= \sqrt{N}/N_\mathrm{tot}$  (Poisson)  "
                  r"[first file]", fontsize=10)
ax_snap.legend(loc="upper right")
add_panel_accent(ax_snap, ESS_GREEN)

ax_all.set_xlabel(r"$n_\mathrm{points}$ per track", fontsize=10)
ax_all.set_ylabel("Normalised count", fontsize=10)
ax_all.set_title("Overall $n_\\mathrm{points}$  (no energy cut)\n"
                 r"$y$-error $= \sqrt{N}/N_\mathrm{tot}$  (Poisson)", fontsize=10)
ax_all.legend(loc="upper right")
add_panel_accent(ax_all, "#C45AB3")


# ─── Panel labels (a, b, c, d) ───────────────────────────────────────────────

for ax, lbl in zip([ax_mean, ax_surv, ax_snap, ax_all], list("abcd")):
    ax.text(-0.12, 1.04, f"({lbl})", transform=ax.transAxes,
            fontsize=12, fontweight="bold", color=ESS_NAVY, va="bottom")


# ─── ESS Logo ─────────────────────────────────────────────────────────────────

if os.path.exists(LOGO_PATH):
    try:
        logo = mpimg.imread(LOGO_PATH)
        # Place in header bar, left side
        imagebox = OffsetImage(logo, zoom=0.12)
        ab = AnnotationBbox(
            imagebox, (0.01, 0.5),
            xycoords=ax_header.transAxes,
            box_alignment=(0, 0.5),
            frameon=False)
        ax_header.add_artist(ab)
    except Exception as e:
        print(f"  Logo not loaded: {e}")
else:
    # Typographic fallback — styled ESS wordmark in the header
    ax_header.text(0.02, 0.5, "ESS",
                   transform=ax_header.transAxes,
                   ha="left", va="center",
                   color=ESS_CYAN, fontsize=14, fontweight="bold",
                   fontfamily="serif")


# ─── Footer strip ─────────────────────────────────────────────────────────────

fig.text(0.99, 0.005,
         f"Noise cut: {NOISE_SIGMA}σ  |  Threshold steps: {N_THRESHOLDS}  "
         f"|  HIBEAM TPC beam test, CCB IFJ PAN Kraków",
         ha="right", va="bottom",
         color=ESS_SILVER, fontsize=6.5, fontstyle="italic")


# ─── Save ─────────────────────────────────────────────────────────────────────

plt.savefig(OUTPUT_PATH, dpi=200, bbox_inches="tight",
            facecolor=ESS_WHITE)
plt.close()
print(f"\n✔  Saved → {OUTPUT_PATH}")