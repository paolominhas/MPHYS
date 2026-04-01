#!/usr/bin/env python3
"""
tpc_dedx_pub.py — Publication-quality TPC ΣEdep / dE/dx analysis
=================================================================
Fits a Landau (Moyal⊗Gaussian) distribution to the truncated (bottom
70 %) energy-loss spectrum from the HIBEAM TPC simulation, with the
Bethe–Bloch mean shown as a physics reference.

The Gaussian convolution absorbs detector resolution, pad response,
and the kinetic-energy spread of the pion spectrum [Refs 3, 6].

Plotting uses mplhep ATLAS style [1] from SFT/CERN (Scikit-HEP).

Usage:  python tpc_dedx_pub.py <file.root> [output_stem]

References
----------
[1] mplhep — Scikit-HEP, https://github.com/scikit-hep/mplhep
[2] L.D. Landau, J. Phys. USSR 8, 201 (1944)
[3] J.E. Moyal, Phil. Mag. 46, 263 (1955)
[4] ALICE Collaboration, JINST 5, P09002 (2010)
[5] PDG, Passage of Particles Through Matter (2024)
[6] S. Meroli et al., JINST 6, P06013 (2011)
"""

import sys, os
import numpy as np
import awkward as ak
import uproot
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.offsetbox import OffsetImage, AnnotationBbox
import mplhep as hep
from scipy import stats, optimize

# ── Style ─────────────────────────────────────────────────────────────────
hep.style.use("ATLAS")
plt.rcParams.update({
    "axes.linewidth":      0.8,
    "xtick.major.width":   0.8,
    "ytick.major.width":   0.8,
    "xtick.minor.width":   0.5,
    "ytick.minor.width":   0.5,
    "xtick.direction":     "in",
    "ytick.direction":     "in",
    "xtick.top":           True,
    "ytick.right":         True,
    "xtick.minor.visible": True,
    "ytick.minor.visible": True,
    "font.size":           13,
    "axes.labelsize":      15,
    "legend.fontsize":     10,
    "figure.facecolor":    "white",
})

# ── Physics constants ─────────────────────────────────────────────────────
M_E   = 0.511          # electron mass [MeV/c²]
M_PI  = 139.570        # charged pion mass [MeV/c²]
K_BB  = 0.307075       # [MeV cm² / mol]

# ── Ar/CO₂ 90:10 gas (Bragg additivity) ──────────────────────────────────
GAS_ZA  = 0.9*(18/39.948) + 0.1*(22/44.010)
GAS_I   = np.exp(0.9*np.log(188.0) + 0.1*np.log(85.0)) * 1e-6  # MeV
GAS_RHO = 1.662e-3  # g/cm³


# ═══════════════════════════════════════════════════════════════════════════
#  Bethe–Bloch  [PDG Eq. 34.5]
# ═══════════════════════════════════════════════════════════════════════════

def bethe_mean_MeV(T_MeV, path_cm, M=M_PI, rho=GAS_RHO):
    """⟨ΔE⟩ [MeV] for a pion of kinetic energy T traversing path_cm of gas."""
    g  = 1.0 + T_MeV / M
    b2 = 1.0 - 1.0/g**2
    bg2 = b2 * g**2
    Wm = 2*M_E*bg2 / (1 + 2*g*M_E/M + (M_E/M)**2)
    dEdx = K_BB * GAS_ZA / b2 * (0.5*np.log(2*M_E*bg2*Wm/GAS_I**2) - b2)
    return dEdx * rho * path_cm


# ═══════════════════════════════════════════════════════════════════════════
#  Data loading  (awkward v1/v2 compatible)
# ═══════════════════════════════════════════════════════════════════════════

def load_edep(filepath, tree="hibeam", min_hits=5):
    """Load per-event ΣEdep [MeV] from TPC branch. Returns 1-D numpy array."""
    with uproot.open(filepath) as f:
        t = f[tree]
        edep_key = next(
            (k for k in t.keys(recursive=True)
             if "Edep" in k and "TPC" in k
             and "Proto" not in k and "target" not in k), None)
        if edep_key is None:
            raise KeyError(f"TPC/Edep not found. Keys: {t.keys()}")
        edep = t[edep_key].array(library="ak")

    mask = ak.to_numpy(ak.num(edep) >= min_hits)
    sums = ak.to_numpy(ak.sum(edep[mask], axis=1)).astype(np.float64)
    print(f"  Events: {len(edep)}  →  {int(mask.sum())} after ≥{min_hits}-hit cut")
    return sums[sums > 0]


# ═══════════════════════════════════════════════════════════════════════════
#  Moyal ⊗ Gaussian fit  [Refs 2, 3, 6]
# ═══════════════════════════════════════════════════════════════════════════

def _moyal_gauss(x, loc, scale, sigma, norm):
    """
    Moyal (Landau approx.) convolved with Gaussian of width sigma.
    The convolution absorbs detector resolution and energy-spectrum
    smearing.  Standard parameterisation in TPC dE/dx [Refs 4, 6].

    Parameters: loc   — MPV of the underlying Moyal
                scale — Landau width ξ
                sigma — Gaussian resolution σ
                norm  — overall normalisation
    """
    dx = x[1] - x[0] if len(x) > 1 else 1e-4
    # Extend grid to avoid edge artefacts from convolution
    pad = max(int(6*sigma/dx), 30)
    x_ext = np.linspace(x[0] - pad*dx, x[-1] + pad*dx, len(x) + 2*pad)
    moyal = stats.moyal.pdf(x_ext, loc=loc, scale=scale)
    kern  = stats.norm.pdf(np.arange(-pad, pad+1)*dx, 0, max(sigma, 1e-8))
    kern /= kern.sum()
    conv  = np.convolve(moyal, kern, mode="same")
    return norm * np.interp(x, x_ext, conv)


def fit_landau(values, trunc=0.70, n_bins=80, x_max=None):
    """
    Fit Moyal⊗Gaussian to the bottom `trunc` fraction of data.
    Returns dict with fit results and arrays for plotting.
    """
    cut = np.percentile(values, trunc * 100)
    if x_max is None:
        x_max = np.percentile(values, 97)

    counts, edges = np.histogram(values, bins=n_bins, range=(0, x_max))
    bc = 0.5*(edges[:-1] + edges[1:])
    bw = edges[1] - edges[0]

    fm     = bc <= cut                    # fit-region mask
    usable = fm & (counts >= 5)           # Pearson chi² validity
    xf, yf = bc[usable], counts[usable].astype(float)

    # Initial guesses
    loc0   = float(bc[fm][np.argmax(counts[fm])])
    scale0 = bw * 3
    sigma0 = bw * 1.5                    # resolution ~ 1–2 bins
    norm0  = float(counts.sum() * bw)

    popt, pcov = optimize.curve_fit(
        _moyal_gauss, xf, yf,
        p0=[loc0, scale0, sigma0, norm0],
        sigma=np.sqrt(np.maximum(yf, 1.0)),
        absolute_sigma=True, maxfev=30000,
        bounds=([0, 1e-6, 1e-6, 0],       # all positive
                [np.inf, np.inf, np.inf, np.inf]),
    )
    perr = np.sqrt(np.diag(pcov))

    # Smooth curve for plotting
    xc = np.linspace(0, x_max, 2000)
    yc = _moyal_gauss(xc, *popt)

    # Chi²/ndf in fit region  (4 free parameters)
    y_exp = _moyal_gauss(xf, *popt)
    chi2  = float(np.sum((yf - y_exp)**2 / np.maximum(y_exp, 1e-6)))
    ndf   = len(xf) - 4

    # 1σ error band (numerical Jacobian)
    eps = 1e-6
    jac = np.zeros((len(xc), 4))
    for i in range(4):
        h = max(abs(popt[i])*eps, eps)
        p_up = popt.copy(); p_up[i] += h
        p_dn = popt.copy(); p_dn[i] -= h
        jac[:, i] = (_moyal_gauss(xc, *p_up) -
                      _moyal_gauss(xc, *p_dn)) / (2*h)
    band = np.sqrt(np.maximum(np.einsum("ij,jk,ik->i", jac, pcov, jac), 0))

    return dict(
        loc=popt[0], scale=popt[1], sigma=popt[2], norm=popt[3],
        loc_err=perr[0], scale_err=perr[1], sigma_err=perr[2],
        mpv=popt[0], popt=popt, pcov=pcov,
        chi2=chi2, ndf=ndf,
        p_value=float(stats.chi2.sf(chi2, max(ndf, 1))),
        cut=cut, trunc=trunc,
        centres=bc, counts=counts, edges=edges, bw=bw,
        fit_mask=fm, x_curve=xc, y_curve=yc, band=band,
    )


# ═══════════════════════════════════════════════════════════════════════════
#  ESS logo helper
# ═══════════════════════════════════════════════════════════════════════════

# Default path — override via ESS_LOGO env var or argument
_LOGO_PATH = os.environ.get(
    "ESS_LOGO",
    os.path.join(os.path.dirname(os.path.abspath(__file__)), "ess_logo.png")
)

def add_logo(ax, logo_path=_LOGO_PATH, height_inches=0.35, loc=(0.03, 0.03)):
    """
    Place the ESS logo in the lower-left corner, scaled to a fixed
    physical height (default 0.35 in ≈ 9 mm) regardless of the
    original image resolution.

    To use: place ess_logo.png next to this script, or set the
    ESS_LOGO environment variable, or pass logo_path explicitly.
    """
    if os.path.isfile(logo_path):
        img = plt.imread(logo_path)
        fig = ax.get_figure()
        dpi = fig.get_dpi()
        # Scale so the image renders at exactly height_inches tall
        zoom = (height_inches * dpi) / img.shape[0]
        imagebox = OffsetImage(img, zoom=zoom, alpha=0.8)
        imagebox.image.axes = ax
        ab = AnnotationBbox(
            imagebox, loc,
            xycoords="axes fraction",
            box_alignment=(0, 0),    # anchor lower-left
            frameon=False, pad=0,
        )
        ax.add_artist(ab)
    else:
        # Fallback: small italic text
        ax.text(loc[0], loc[1], r"$\it{ESS}$",
                transform=ax.transAxes, fontsize=11,
                color="#003366", alpha=0.6, va="bottom", ha="left")


# ═══════════════════════════════════════════════════════════════════════════
#  Publication figure
# ═══════════════════════════════════════════════════════════════════════════

def plot(fit, bethe_mean=None, output=None, logo_path=_LOGO_PATH):
    """
    Two-panel figure (distribution + pulls) in mplhep ATLAS style.
    Outputs both linear and log-y as vector PDF.
    """
    bc  = fit["centres"]
    cnt = fit["counts"].astype(float)
    fm  = fit["fit_mask"]
    bw  = fit["bw"]
    x_max = fit["x_curve"][-1]

    # Pulls using the convolved model
    fc_binned = _moyal_gauss(bc, *fit["popt"])
    denom     = np.sqrt(np.maximum(fc_binned, 1.0))
    pulls     = (cnt - fc_binned) / denom

    C_DATA, C_TAIL = "#1f77b4", "#d62728"
    trunc_pct = int(fit["trunc"] * 100)

    for log_y in [False, True]:
        fig = plt.figure(figsize=(8, 6.5))
        gs  = gridspec.GridSpec(2, 1, height_ratios=[3.5, 1], hspace=0.06,
                                left=0.13, right=0.96, top=0.92, bottom=0.11)
        ax  = fig.add_subplot(gs[0])
        axr = fig.add_subplot(gs[1], sharex=ax)

        # ── Histogram ─────────────────────────────────────────────────
        ax.bar(bc[fm],  cnt[fm],  width=bw*0.92,
               color=C_DATA, alpha=0.6, lw=0,
               label=f"Data (bottom {trunc_pct}%)")
        ax.bar(bc[~fm], cnt[~fm], width=bw*0.92,
               color=C_TAIL, alpha=0.3, lw=0,
               label=f"Excluded tail (top {100-trunc_pct}%)")

        # ── Fit curve + 1σ band ───────────────────────────────────────
        ax.plot(fit["x_curve"], fit["y_curve"],
                "k-", lw=1.8, label=r"Moyal $\otimes$ Gauss fit")
        ax.fill_between(fit["x_curve"],
                        fit["y_curve"] - fit["band"],
                        fit["y_curve"] + fit["band"],
                        color="black", alpha=0.10)

        # ── Truncation line ───────────────────────────────────────────
        ax.axvline(fit["cut"], color="grey", ls="--", lw=1.0, alpha=0.7)

        # ── Bethe–Bloch mean ──────────────────────────────────────────
        if bethe_mean is not None:
            ax.axvline(bethe_mean, color="#2ca02c", ls="-.", lw=1.5,
                       label=(r"$\langle\Delta E\rangle_{\rm Bethe}$"
                              f" = {bethe_mean:.3f} MeV"))

        # ── MPV ───────────────────────────────────────────────────────
        ax.axvline(fit["mpv"], color="#ff7f0e", ls=":", lw=1.5,
                   label=f"MPV = {fit['mpv']:.4f} MeV")

        # ── Axes & legend ─────────────────────────────────────────────
        ax.legend(loc="upper right", frameon=True, framealpha=0.9,
                  edgecolor="#cccccc", borderpad=0.5)
        ax.set_ylabel("Events / bin")
        ax.set_xlim(0, x_max)
        plt.setp(ax.get_xticklabels(), visible=False)
        if log_y:
            ax.set_yscale("log")
            ax.set_ylim(0.5, cnt.max() * 5)
        else:
            ax.set_ylim(0, cnt.max() * 1.35)

        # ── Fit quality box ───────────────────────────────────────────
        chi2r = fit["chi2"] / max(fit["ndf"], 1)
        txt = (
            rf"$\chi^2 / n_{{\rm df}} = {fit['chi2']:.1f}\,/\,"
            rf"{fit['ndf']} = {chi2r:.2f}$"  "\n"
            rf"$p = {fit['p_value']:.3f}$"  "\n"
            rf"$\xi = {fit['scale']:.4f}\pm{fit['scale_err']:.4f}$ MeV"  "\n"
            rf"$\sigma = {fit['sigma']:.4f}\pm{fit['sigma_err']:.4f}$ MeV"
        )
        ax.text(0.97, 0.58, txt,
                transform=ax.transAxes, ha="right", va="top", fontsize=9.5,
                bbox=dict(boxstyle="round,pad=0.4", fc="white",
                          ec="#bbbbbb", alpha=0.9))

        # ── Experiment label ──────────────────────────────────────────
        ax.text(0.05, 0.92, r"$\bf{HIBEAM}$ Simulation",
                transform=ax.transAxes, fontsize=14, va="top")

        # ── ESS logo (scaled to ~9 mm height) ────────────────────────
        add_logo(ax, logo_path=logo_path)

        # ── Pull panel ────────────────────────────────────────────────
        axr.bar(bc[fm],  pulls[fm],  width=bw*0.92,
                color=C_DATA, alpha=0.6, lw=0)
        axr.bar(bc[~fm], pulls[~fm], width=bw*0.92,
                color=C_TAIL, alpha=0.3, lw=0)
        axr.axhline(0,  color="black", lw=0.8)
        axr.axhline(+2, color="grey",  lw=0.7, ls="--")
        axr.axhline(-2, color="grey",  lw=0.7, ls="--")
        axr.axvline(fit["cut"], color="grey", ls="--", lw=1.0, alpha=0.7)
        axr.fill_between([0, x_max], -1, 1, color="grey", alpha=0.08)
        axr.set_ylabel(r"$(N\!-\!F)/\sqrt{F}$")
        axr.set_xlabel(r"$\sum E_{\mathrm{dep}}$ [MeV]")
        axr.set_ylim(-5, 5)
        axr.set_xlim(0, x_max)

        # ── Save ──────────────────────────────────────────────────────
        if output:
            tag = "log" if log_y else "linear"
            path = f"{output}_{tag}.pdf"
            fig.savefig(path, bbox_inches="tight")
            print(f"  ✔ {path}")
        plt.close(fig)


# ═══════════════════════════════════════════════════════════════════════════
#  Main
# ═══════════════════════════════════════════════════════════════════════════

def main():
    if len(sys.argv) < 2:
        print("Usage: python tpc_dedx_pub.py <file.root> [output_stem]")
        sys.exit(1)

    root_path = sys.argv[1]
    out_stem  = sys.argv[2] if len(sys.argv) > 2 else "tpc_edep"

    print(f"Loading: {root_path}")
    values = load_edep(root_path)
    print(f"  Accepted events: {len(values):,}")

    # ── Bethe–Bloch reference ─────────────────────────────────────────
    # Adjust T_KIN and PATH to match your TPC geometry in hibeam_g4
    T_KIN  = 200.0    # MeV  (typical annihilation pion)
    PATH   = 20.0     # cm   (effective gas path)
    bb     = bethe_mean_MeV(T_KIN, PATH)
    print(f"  Bethe–Bloch ⟨ΔE⟩ for {T_KIN:.0f} MeV π± "
          f"in {PATH:.0f} cm Ar/CO₂: {bb:.4f} MeV")

    # ── Fit (80 bins, Moyal⊗Gauss) ───────────────────────────────────
    result = fit_landau(values, trunc=0.70, n_bins=80, x_max=0.25)
    chi2r  = result["chi2"] / max(result["ndf"], 1)
    print(f"  MPV   = {result['mpv']:.4f} MeV")
    print(f"  ξ     = {result['scale']:.4f} ± {result['scale_err']:.4f} MeV")
    print(f"  σ_det = {result['sigma']:.4f} ± {result['sigma_err']:.4f} MeV")
    print(f"  χ²/ndf = {result['chi2']:.1f}/{result['ndf']}  "
          f"= {chi2r:.2f}  (p = {result['p_value']:.3f})")

    # ── Plot ──────────────────────────────────────────────────────────
    plot(result, bethe_mean=bb, output=out_stem)
    print("Done.")


if __name__ == "__main__":
    main()