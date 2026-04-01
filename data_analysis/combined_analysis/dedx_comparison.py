#!/usr/bin/env python3
"""
dedx_comparison.py — Proper per-step dE/dx from data and simulation
====================================================================
Computes dE/dx using consecutive-hit 3D distances (not total_length/N)
for both experimental TrackData and simulated ProtoTPCHit, then overlays
the two Landau spectra with a Moyal⊗Gaussian fit.

Data path  :  PyROOT  → TrackData (charge[j] / ds_j)  →  ADC/mm
Sim  path  :  uproot  → ProtoTPCHit (Edep[j] / ds_j)  →  MeV/mm

Usage:
    python dedx_comparison.py \
        --data  <tracks_centroids_*.root> \
        --sim   <simulation.root> \
        --headers recovered_headers/ \
        [--output dedx_proper]

Requires: ROOT (PyROOT), uproot, awkward, numpy, scipy, matplotlib, mplhep
"""

import os, sys, argparse
import numpy as np
import awkward as ak
import uproot
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import mplhep as hep
from scipy import stats, optimize


# ═══════════════════════════════════════════════════════════════════════════
#  Configuration  (single dict — threaded through every function)
# ═══════════════════════════════════════════════════════════════════════════
DEFAULT_CFG = dict(
    min_track_pts = 5,       # minimum centroids per track
    chi2_ndf_max  = 25.0,    # track quality cut
    min_ds_mm     = 0.5,     # reject unphysically short steps [mm]
    max_ds_mm     = 50.0,    # reject gap-jumps [mm]
    trunc_frac    = 0.70,    # truncated-mean fraction (bottom 70 %)
    n_bins        = 80,      # histogram bins
    sim_min_hits  = 5,       # minimum G4 hits per event
)


# ═══════════════════════════════════════════════════════════════════════════
#  1. EXPERIMENTAL DATA  (PyROOT — required for custom streamer classes)
# ═══════════════════════════════════════════════════════════════════════════

def load_data_dedx(filepath, headers_dir, cfg):
    """
    Per-centroid dE/dx [ADC/mm] using consecutive 3D distances.

    For each track:
      • apply chi²/ndf quality cut
      • compute ds_j = |r_{j} − r_{j-1}|  for j = 1 … N-1
      • dE/dx_j = charge[j] / ds_j

    Returns 1-D numpy array of all accepted dE/dx values.
    """
    import ROOT

    headers_so = os.path.join(headers_dir, "recovered_headers.so")
    hdr_file   = os.path.join(headers_dir, "TrackData.h")

    if not os.path.isfile(headers_so):
        raise FileNotFoundError(f"recovered_headers.so not found in {headers_dir}")
    if not os.path.isfile(hdr_file):
        raise FileNotFoundError(f"TrackData.h not found in {headers_dir}")

    ROOT.gInterpreter.AddIncludePath(headers_dir)
    ROOT.gInterpreter.ProcessLine(f'#include "{hdr_file}"')
    ROOT.gSystem.Load(headers_so)

    tfile = ROOT.TFile.Open(filepath)
    if not tfile or tfile.IsZombie():
        raise IOError(f"Cannot open {filepath}")
    tree = tfile.Get("trackingData")
    if not tree:
        tfile.Close()
        raise KeyError("'trackingData' tree not found")

    vec = ROOT.std.vector("TrackData")()
    tree.SetBranchAddress("tracks", vec)

    n_entries = tree.GetEntries()
    dedx_buf = np.empty(n_entries * 50, dtype=np.float64)
    idx = 0

    n_tracks_total = 0
    n_tracks_pass  = 0

    min_pts  = cfg["min_track_pts"]
    chi2_max = cfg["chi2_ndf_max"]
    ds_lo    = cfg["min_ds_mm"]
    ds_hi    = cfg["max_ds_mm"]

    for i in range(n_entries):
        tree.GetEntry(i)
        for track in vec:
            n = track.nPoints
            n_tracks_total += 1
            if n < min_pts:
                continue

            # ── Unpack coordinates once ───────────────────────────────
            xs = np.array([float(track.x[j]) for j in range(n)])
            ys = np.array([float(track.y[j]) for j in range(n)])
            zs = np.array([float(track.z[j]) for j in range(n)])
            qs = np.array([float(track.charge[j]) for j in range(n)])

            # ── Chi²/ndf quality cut ──────────────────────────────────
            sx = np.array([max(float(track.sigmas_x[j]), 1e-4) for j in range(n)])
            sz = np.array([max(float(track.sigmas_z[j]), 1e-4) for j in range(n)])
            res_x = (xs - (track.slope_xy * ys + track.intercept_xy)) / sx
            res_z = (zs - (track.slope_zy * ys + track.intercept_zy)) / sz
            chi2 = float(np.sum(res_x**2) + np.sum(res_z**2))
            ndf  = max(2 * n - 4, 1)
            if chi2 / ndf > chi2_max:
                continue

            n_tracks_pass += 1

            # ── Sort by padrow (y) so consecutive steps make sense ────
            order = np.argsort(ys)
            xs, ys, zs, qs = xs[order], ys[order], zs[order], qs[order]

            # ── Consecutive 3D distances ──────────────────────────────
            ds = np.sqrt(np.diff(xs)**2 + np.diff(ys)**2 + np.diff(zs)**2)

            # Use average of adjacent steps for interior points;
            # for first/last point, use the single adjacent step.
            ds_per_point = np.empty(n)
            ds_per_point[0]    = ds[0]
            ds_per_point[-1]   = ds[-1]
            if n > 2:
                ds_per_point[1:-1] = 0.5 * (ds[:-1] + ds[1:])

            # ── Filter and fill ───────────────────────────────────────
            good = ((qs > 0)
                    & (ds_per_point > ds_lo)
                    & (ds_per_point < ds_hi))
            vals = qs[good] / ds_per_point[good]

            need = idx + vals.size
            if need > dedx_buf.size:
                dedx_buf = np.resize(dedx_buf, need * 2)
            dedx_buf[idx:idx + vals.size] = vals
            idx += vals.size

    tfile.Close()
    result = dedx_buf[:idx].copy()
    print(f"  Data: {n_tracks_total} tracks → {n_tracks_pass} pass quality "
          f"→ {len(result):,} dE/dx values")
    return result


# ═══════════════════════════════════════════════════════════════════════════
#  2. SIMULATION  (uproot + awkward — fully vectorised)
# ═══════════════════════════════════════════════════════════════════════════

def load_sim_dedx(filepath, cfg, tree_name="hibeam",
                  sim_headers_dir=None):
    """
    Per-step dE/dx [MeV/mm] from Geant4 hit positions.

    Strategy:
      1. Try uproot (fast, vectorised).
      2. If uproot returns flat arrays (custom streamer class like
         ProtoTPCHit), fall back to PyROOT with simulation headers.

    Returns (1-D numpy array of dE/dx, weights or None).
    """
    ds_lo = cfg["min_ds_mm"]
    ds_hi = cfg["max_ds_mm"]
    min_hits = cfg["sim_min_hits"]

    # ── Attempt uproot first ──────────────────────────────────────────
    uproot_ok = False
    try:
        edep, px, py, pz, layer, weights = _load_sim_uproot(
            filepath, tree_name)
        # Check that arrays are actually jagged (not flat from a
        # custom streamer that uproot couldn't deserialize).
        # ak.num() raises ValueError on flat arrays — exactly the
        # error that triggers this fallback.
        try:
            ak.num(edep)
            uproot_ok = True
        except (ValueError, AttributeError):
            print("  uproot returned flat arrays (custom streamer) "
                  "→ falling back to PyROOT")
    except Exception as e:
        print(f"  uproot read failed ({e}) → falling back to PyROOT")

    if not uproot_ok:
        if sim_headers_dir is None:
            # Try auto-detecting next to the script
            guess = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                 "recovered_headers_simulation")
            if os.path.isdir(guess):
                sim_headers_dir = guess
            else:
                raise RuntimeError(
                    "uproot cannot read this simulation file's custom "
                    "classes.  Pass --sim-headers <dir> pointing to the "
                    "recovered_headers_simulation/ directory.")
        return _load_sim_pyroot(filepath, tree_name, sim_headers_dir, cfg)

    # ── From here: uproot succeeded with jagged arrays ────────────────
    nhits = ak.num(edep)
    mask  = nhits >= min_hits
    edep  = edep[mask]
    px, py, pz = px[mask], py[mask], pz[mask]
    if layer is not None:
        layer = layer[mask]
    if weights is not None:
        weights = weights[ak.to_numpy(mask)]

    if layer is not None:
        sort_idx = ak.argsort(layer, axis=1)
        edep = edep[sort_idx]
        px, py, pz = px[sort_idx], py[sort_idx], pz[sort_idx]

    result, n_ev = _compute_dedx_awkward(edep, px, py, pz, ds_lo, ds_hi)
    print(f"  Sim: {n_ev} events → {len(result):,} dE/dx values")
    return result, weights


def _load_sim_uproot(filepath, tree_name):
    """Read simulation branches with uproot. Returns awkward arrays."""
    with uproot.open(filepath) as f:
        t = f[tree_name]
        keys = t.keys(recursive=True)

        def find_key(pattern):
            for k in keys:
                if pattern in k and "target" not in k.lower():
                    return k
            return None

        edep_key  = find_key("Edep")
        posx_key  = find_key("Pos_X")
        posy_key  = find_key("Pos_Y")
        posz_key  = find_key("Pos_Z")
        layer_key = find_key("Layer")

        missing = []
        for name, key in [("Edep", edep_key), ("Pos_X", posx_key),
                          ("Pos_Y", posy_key), ("Pos_Z", posz_key)]:
            if key is None:
                missing.append(name)
        if missing:
            raise KeyError(f"Missing branches: {missing}")

        edep  = t[edep_key].array(library="ak")
        px    = t[posx_key].array(library="ak")
        py    = t[posy_key].array(library="ak")
        pz    = t[posz_key].array(library="ak")
        layer = t[layer_key].array(library="ak") if layer_key else None

    weights = None
    try:
        with uproot.open(filepath) as f:
            w = f[tree_name]["PrimaryWeight"].array(library="ak")
            weights = ak.to_numpy(ak.fill_none(ak.firsts(w), 1.0))
    except Exception:
        pass

    return edep, px, py, pz, layer, weights


def _load_sim_pyroot(filepath, tree_name, sim_headers_dir, cfg):
    """
    Fallback: read ProtoTPCHit via PyROOT + recovered simulation headers.
    Returns (1-D numpy array of dE/dx, weights or None).
    """
    import ROOT

    headers_so = os.path.join(sim_headers_dir, "recovered_headers.so")
    hdr_file   = os.path.join(sim_headers_dir, "ProtoTPCHit.h")

    if not os.path.isfile(headers_so):
        raise FileNotFoundError(
            f"recovered_headers.so not found in {sim_headers_dir}")

    ROOT.gInterpreter.AddIncludePath(sim_headers_dir)
    if os.path.isfile(hdr_file):
        ROOT.gInterpreter.ProcessLine(f'#include "{hdr_file}"')
    ROOT.gSystem.Load(headers_so)

    tfile = ROOT.TFile.Open(filepath)
    if not tfile or tfile.IsZombie():
        raise IOError(f"Cannot open {filepath}")
    tree = tfile.Get(tree_name)
    if not tree:
        tfile.Close()
        raise KeyError(f"'{tree_name}' tree not found")

    ds_lo    = cfg["min_ds_mm"]
    ds_hi    = cfg["max_ds_mm"]
    min_hits = cfg["sim_min_hits"]

    # Detect branch structure: could be a ProtoTPCHit object branch
    # or separate leaf branches — try the object route first.
    branch_names = [b.GetName() for b in tree.GetListOfBranches()]

    # Find the ProtoTPC branch (object holding vectors)
    proto_branch = None
    for bn in branch_names:
        if "ProtoTPC" in bn and "target" not in bn.lower():
            proto_branch = bn
            break

    n_entries = tree.GetEntries()
    dedx_buf  = np.empty(n_entries * 30, dtype=np.float64)
    idx       = 0
    n_events  = 0

    # Try reading PrimaryWeight
    weights = None
    try:
        w_arr = np.empty(n_entries, dtype=np.float64)
        has_weight = tree.GetBranch("PrimaryWeight") is not None
    except Exception:
        has_weight = False

    for i in range(n_entries):
        tree.GetEntry(i)

        # Get the ProtoTPCHit object
        hit = getattr(tree, proto_branch) if proto_branch else None
        if hit is None:
            continue

        n = hit.nHits
        if n < min_hits:
            continue

        # Unpack vectors
        edep  = np.array([float(hit.Edep[j])  for j in range(n)])
        pos_x = np.array([float(hit.Pos_X[j]) for j in range(n)])
        pos_y = np.array([float(hit.Pos_Y[j]) for j in range(n)])
        pos_z = np.array([float(hit.Pos_Z[j]) for j in range(n)])

        # Sort by Layer if available
        try:
            layers = np.array([int(hit.Layer[j]) for j in range(n)])
            order  = np.argsort(layers)
            edep, pos_x, pos_y, pos_z = (
                edep[order], pos_x[order], pos_y[order], pos_z[order])
        except Exception:
            pass

        # Consecutive 3D distances
        ds = np.sqrt(np.diff(pos_x)**2 + np.diff(pos_y)**2
                     + np.diff(pos_z)**2)
        ds_per = np.empty(n)
        ds_per[0]  = ds[0]
        ds_per[-1] = ds[-1]
        if n > 2:
            ds_per[1:-1] = 0.5 * (ds[:-1] + ds[1:])

        good = (edep > 0) & (ds_per > ds_lo) & (ds_per < ds_hi)
        vals = edep[good] / ds_per[good]

        need = idx + vals.size
        if need > dedx_buf.size:
            dedx_buf = np.resize(dedx_buf, need * 2)
        dedx_buf[idx:idx + vals.size] = vals
        idx += vals.size
        n_events += 1

    tfile.Close()
    result = dedx_buf[:idx].copy()
    result = result[np.isfinite(result) & (result > 0)]
    print(f"  Sim (PyROOT): {n_events} events → {len(result):,} dE/dx values")
    return result, weights


def _compute_dedx_awkward(edep, px, py, pz, ds_lo, ds_hi):
    """Vectorised dE/dx from jagged awkward arrays. Returns (flat array, n_events)."""
    dx = px[:, 1:] - px[:, :-1]
    dy = py[:, 1:] - py[:, :-1]
    dz = pz[:, 1:] - pz[:, :-1]
    ds = np.sqrt(dx**2 + dy**2 + dz**2)

    ds_left  = ak.concatenate([ds[:, :1], ds], axis=1)
    ds_right = ak.concatenate([ds, ds[:, -1:]], axis=1)
    ds_per_hit = 0.5 * (ds_left + ds_right)

    dedx    = edep / ds_per_hit
    step_ok = (ds_per_hit > ds_lo) & (ds_per_hit < ds_hi) & (edep > 0)

    result = ak.to_numpy(ak.flatten(dedx[step_ok])).astype(np.float64)
    result = result[np.isfinite(result) & (result > 0)]
    return result, len(edep)


# ═══════════════════════════════════════════════════════════════════════════
#  3. MOYAL ⊗ GAUSSIAN FIT
# ═══════════════════════════════════════════════════════════════════════════

def moyal_gauss(x, loc, scale, sigma, norm):
    """Moyal (Landau approx.) convolved with Gaussian resolution."""
    dx  = x[1] - x[0] if len(x) > 1 else 1e-4
    pad = max(int(6 * sigma / dx), 30)
    x_ext = np.linspace(x[0] - pad*dx, x[-1] + pad*dx, len(x) + 2*pad)
    moyal = stats.moyal.pdf(x_ext, loc=loc, scale=scale)
    kern  = stats.norm.pdf(np.arange(-pad, pad+1) * dx, 0, max(sigma, 1e-8))
    kern /= kern.sum()
    conv  = np.convolve(moyal, kern, mode="same")
    return norm * np.interp(x, x_ext, conv)


def fit_spectrum(values, cfg, label=""):
    """
    Fit Moyal⊗Gauss to the bottom `trunc_frac` fraction.
    Returns dict with everything needed for plotting.
    """
    trunc  = cfg["trunc_frac"]
    n_bins = cfg["n_bins"]

    cut   = np.percentile(values, trunc * 100)
    x_max = np.percentile(values, 97)

    counts, edges = np.histogram(values, bins=n_bins, range=(0, x_max))
    bc = 0.5 * (edges[:-1] + edges[1:])
    bw = edges[1] - edges[0]

    fm     = bc <= cut
    usable = fm & (counts >= 5)
    xf, yf = bc[usable], counts[usable].astype(float)

    loc0   = float(bc[fm][np.argmax(counts[fm])])
    scale0 = bw * 3
    sigma0 = bw * 1.5
    norm0  = float(counts.sum() * bw)

    try:
        popt, pcov = optimize.curve_fit(
            moyal_gauss, xf, yf,
            p0=[loc0, scale0, sigma0, norm0],
            sigma=np.sqrt(np.maximum(yf, 1.0)),
            absolute_sigma=True, maxfev=30_000,
            bounds=([0, 1e-6, 1e-6, 0],
                    [np.inf, np.inf, np.inf, np.inf]),
        )
        perr = np.sqrt(np.diag(pcov))
    except RuntimeError:
        print(f"  ⚠ {label} fit did not converge, using initial guess")
        popt = np.array([loc0, scale0, sigma0, norm0])
        pcov = np.zeros((4, 4))
        perr = np.zeros(4)

    xc = np.linspace(0, x_max, 2000)
    yc = moyal_gauss(xc, *popt)

    y_exp = moyal_gauss(xf, *popt)
    chi2  = float(np.sum((yf - y_exp)**2 / np.maximum(y_exp, 1e-6)))
    ndf   = max(len(xf) - 4, 1)

    # 1σ error band
    eps = 1e-6
    jac = np.zeros((len(xc), 4))
    for i in range(4):
        h = max(abs(popt[i]) * eps, eps)
        p_up = popt.copy(); p_up[i] += h
        p_dn = popt.copy(); p_dn[i] -= h
        jac[:, i] = (moyal_gauss(xc, *p_up) -
                     moyal_gauss(xc, *p_dn)) / (2 * h)
    band = np.sqrt(np.maximum(np.einsum("ij,jk,ik->i", jac, pcov, jac), 0))

    print(f"  {label}:  MPV = {popt[0]:.4f}   ξ = {popt[1]:.4f} ± {perr[1]:.4f}"
          f"   σ = {popt[2]:.4f} ± {perr[2]:.4f}"
          f"   χ²/ndf = {chi2:.1f}/{ndf} = {chi2/ndf:.2f}")

    return dict(
        label=label, popt=popt, pcov=pcov, perr=perr,
        mpv=popt[0], scale=popt[1], sigma=popt[2],
        chi2=chi2, ndf=ndf,
        p_value=float(stats.chi2.sf(chi2, max(ndf, 1))),
        cut=cut, trunc=trunc,
        centres=bc, counts=counts, edges=edges, bw=bw,
        fit_mask=fm, x_curve=xc, y_curve=yc, band=band,
    )


# ═══════════════════════════════════════════════════════════════════════════
#  4. CALIBRATION:  MPV-match data (ADC/mm) → sim (MeV/mm)
# ═══════════════════════════════════════════════════════════════════════════

def calibrate_to_sim(data_fit, sim_fit):
    """
    Returns scale factor k [MeV / ADC] so that  data_dedx * k  has
    the same MPV as sim_dedx.
    """
    k = sim_fit["mpv"] / data_fit["mpv"]
    print(f"\n  Pseudo-calibration:")
    print(f"    Data MPV = {data_fit['mpv']:.4f} ADC/mm")
    print(f"    Sim  MPV = {sim_fit['mpv']:.6f} MeV/mm")
    print(f"    k        = {k:.4e} MeV/ADC\n")
    return k


# ═══════════════════════════════════════════════════════════════════════════
#  5. PUBLICATION OVERLAY PLOT
# ═══════════════════════════════════════════════════════════════════════════

def plot_overlay(data_vals, sim_vals, data_fit, sim_fit,
                 calibration_k, cfg, output="dedx_proper"):
    """
    Two-panel overlay (density + ratio) with both distributions on a
    common MeV/mm axis after calibration.
    """
    n_bins    = cfg["n_bins"]
    trunc_pct = int(cfg["trunc_frac"] * 100)

    hep.style.use("ATLAS")
    plt.rcParams.update({
        "axes.linewidth": 0.8, "font.size": 13,
        "axes.labelsize": 15,  "legend.fontsize": 9.5,
        "figure.facecolor": "white",
        "xtick.direction": "in", "ytick.direction": "in",
        "xtick.top": True, "ytick.right": True,
        "xtick.minor.visible": True, "ytick.minor.visible": True,
    })

    # Calibrated data in MeV/mm
    data_cal = data_vals * calibration_k

    # Common x range
    x_max = min(np.percentile(data_cal, 98),
                np.percentile(sim_vals, 98))
    edges = np.linspace(0, x_max, n_bins + 1)
    bc    = 0.5 * (edges[:-1] + edges[1:])
    bw    = edges[1] - edges[0]

    # Density histograms
    d_cnt, _ = np.histogram(data_cal, bins=edges, density=True)
    s_cnt, _ = np.histogram(sim_vals, bins=edges, density=True)

    C_DATA = "#1f77b4"
    C_SIM  = "#d62728"

    for log_y in [False, True]:
        fig = plt.figure(figsize=(9, 7))
        gs  = gridspec.GridSpec(2, 1, height_ratios=[3.5, 1], hspace=0.06,
                                left=0.13, right=0.96, top=0.92, bottom=0.11)
        ax  = fig.add_subplot(gs[0])
        axr = fig.add_subplot(gs[1], sharex=ax)

        # ── Histograms ───────────────────────────────────────────────
        ax.bar(bc, d_cnt, width=bw*0.45, align="edge",
               color=C_DATA, alpha=0.55, lw=0, label="Data (calibrated)")
        ax.bar(bc + bw*0.45, s_cnt, width=bw*0.45, align="edge",
               color=C_SIM, alpha=0.45, lw=0, label="Simulation")

        # ── Fit curves (rescaled to density) ──────────────────────────
        d_refit = fit_spectrum(data_cal, cfg, label="Data (cal)")
        s_refit = fit_spectrum(sim_vals, cfg, label="Sim (refit)")

        ax.plot(d_refit["x_curve"],
                d_refit["y_curve"] / (len(data_cal) * d_refit["bw"]),
                color=C_DATA, ls="-", lw=1.8,
                label=(rf"Data fit: MPV = {d_refit['mpv']:.4f} MeV/mm"))
        ax.fill_between(d_refit["x_curve"],
                        (d_refit["y_curve"] - d_refit["band"]) / (len(data_cal) * d_refit["bw"]),
                        (d_refit["y_curve"] + d_refit["band"]) / (len(data_cal) * d_refit["bw"]),
                        color=C_DATA, alpha=0.10)

        ax.plot(s_refit["x_curve"],
                s_refit["y_curve"] / (len(sim_vals) * s_refit["bw"]),
                color=C_SIM, ls="--", lw=1.8,
                label=(rf"Sim fit: MPV = {s_refit['mpv']:.4f} MeV/mm"))
        ax.fill_between(s_refit["x_curve"],
                        (s_refit["y_curve"] - s_refit["band"]) / (len(sim_vals) * s_refit["bw"]),
                        (s_refit["y_curve"] + s_refit["band"]) / (len(sim_vals) * s_refit["bw"]),
                        color=C_SIM, alpha=0.10)

        # ── Truncation lines ─────────────────────────────────────────
        for fit, c in [(d_refit, C_DATA), (s_refit, C_SIM)]:
            ax.axvline(fit["cut"], color=c, ls=":", lw=0.8, alpha=0.5)

        ax.legend(loc="upper right", frameon=True, framealpha=0.9,
                  edgecolor="#cccccc")
        ax.set_ylabel("Probability density")
        ax.set_xlim(0, x_max)
        plt.setp(ax.get_xticklabels(), visible=False)

        if log_y:
            ax.set_yscale("log")
            ymax = max(d_cnt.max(), s_cnt.max())
            ax.set_ylim(ymax * 1e-4, ymax * 5)
        else:
            ax.set_ylim(0, max(d_cnt.max(), s_cnt.max()) * 1.35)

        # ── Info box ──────────────────────────────────────────────────
        txt = (
            rf"$k = {calibration_k:.3e}$ MeV/ADC" "\n"
            rf"Data $\xi = {d_refit['scale']:.4f}$ MeV/mm" "\n"
            rf"Sim  $\xi = {s_refit['scale']:.4f}$ MeV/mm" "\n"
            rf"Bottom {trunc_pct}% truncation"
        )
        ax.text(0.97, 0.55, txt,
                transform=ax.transAxes, ha="right", va="top", fontsize=9,
                bbox=dict(boxstyle="round,pad=0.4", fc="white",
                          ec="#bbbbbb", alpha=0.9))

        ax.text(0.05, 0.92, r"$\bf{HIBEAM}$ Simulation + Data",
                transform=ax.transAxes, fontsize=14, va="top")

        # ── Ratio panel ───────────────────────────────────────────────
        with np.errstate(divide="ignore", invalid="ignore"):
            ratio = np.where(s_cnt > 0, d_cnt / s_cnt, np.nan)
        axr.bar(bc, ratio - 1, width=bw*0.92, color="grey", alpha=0.5)
        axr.axhline(0, color="black", lw=0.8)
        axr.axhline(+0.2, color="grey", lw=0.6, ls="--")
        axr.axhline(-0.2, color="grey", lw=0.6, ls="--")
        axr.fill_between([0, x_max], -0.1, 0.1, color="grey", alpha=0.08)
        axr.set_ylabel("Data/Sim − 1")
        axr.set_xlabel(r"$\mathrm{d}E/\mathrm{d}x$  [MeV/mm]")
        axr.set_ylim(-1, 1)
        axr.set_xlim(0, x_max)

        tag  = "log" if log_y else "linear"
        path = f"{output}_{tag}.pdf"
        fig.savefig(path, bbox_inches="tight")
        print(f"  ✔ {path}")
        plt.close(fig)


# ═══════════════════════════════════════════════════════════════════════════
#  6. STANDALONE DIAGNOSTIC PANELS  (no calibration needed)
# ═══════════════════════════════════════════════════════════════════════════

def plot_diagnostics(data_vals, sim_vals, cfg, output="dedx_diagnostics.pdf"):
    """
    Side-by-side: raw data (ADC/mm) and raw sim (MeV/mm), each with
    its own Moyal⊗Gauss fit.
    """
    n_bins = cfg["n_bins"]
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))

    for ax, vals, unit, label, color in [
        (axes[0], data_vals, "ADC/mm", "Experimental", "#1f77b4"),
        (axes[1], sim_vals,  "MeV/mm", "Simulation",  "#d62728"),
    ]:
        fit = fit_spectrum(vals, cfg, label=label)
        clip = np.percentile(vals, 97)
        ax.hist(vals[vals < clip], bins=n_bins,
                color=color, alpha=0.5, density=True, label=label)
        ax.plot(fit["x_curve"], fit["y_curve"] / (len(vals) * fit["bw"]),
                "k-", lw=1.5, label=rf"Fit MPV = {fit['mpv']:.4f} {unit}")
        ax.axvline(fit["cut"], color="grey", ls="--", lw=0.8)
        ax.set_xlabel(f"dE/dx  [{unit}]")
        ax.set_ylabel("Probability density")
        ax.legend(fontsize=9)
        ax.set_title(f"{label} — per-step dE/dx")

    plt.tight_layout()
    fig.savefig(output, bbox_inches="tight")
    print(f"  ✔ {output}")
    plt.close(fig)


# ═══════════════════════════════════════════════════════════════════════════
#  Main
# ═══════════════════════════════════════════════════════════════════════════

def main():
    p = argparse.ArgumentParser(
        description="Per-step dE/dx comparison: data vs simulation")
    p.add_argument("--data", required=True,
                   help="Experimental ROOT file (trackingData tree)")
    p.add_argument("--sim", required=True,
                   help="Simulation ROOT file (hibeam tree)")
    p.add_argument("--headers", required=True,
                   help="Directory with recovered_headers.so + TrackData.h")
    p.add_argument("--sim-headers", default=None,
                   help="Directory with simulation recovered_headers.so + "
                        "ProtoTPCHit.h  (auto-detected if next to script)")
    p.add_argument("--sim-tree", default="hibeam",
                   help="Simulation tree name (default: hibeam)")
    p.add_argument("--output", default="dedx_proper",
                   help="Output stem (produces _linear.pdf, _log.pdf)")
    p.add_argument("--min-points", type=int, default=5)
    p.add_argument("--chi2-cut", type=float, default=25.0)
    p.add_argument("--trunc", type=float, default=0.70)
    p.add_argument("--n-bins", type=int, default=80)
    args = p.parse_args()

    # Build config from defaults + CLI overrides
    cfg = dict(DEFAULT_CFG)
    cfg["min_track_pts"] = args.min_points
    cfg["chi2_ndf_max"]  = args.chi2_cut
    cfg["trunc_frac"]    = args.trunc
    cfg["n_bins"]        = args.n_bins

    # ── Load ──────────────────────────────────────────────────────────
    print("Loading experimental data...")
    data_dedx = load_data_dedx(args.data, args.headers, cfg)

    print("Loading simulation...")
    sim_dedx, sim_weights = load_sim_dedx(args.sim, cfg,
                                          tree_name=args.sim_tree,
                                          sim_headers_dir=args.sim_headers)

    # ── Fit each independently ────────────────────────────────────────
    print("\nFitting data (ADC/mm)...")
    data_fit = fit_spectrum(data_dedx, cfg, label="Data")

    print("Fitting simulation (MeV/mm)...")
    sim_fit = fit_spectrum(sim_dedx, cfg, label="Sim")

    # ── Calibrate ─────────────────────────────────────────────────────
    k = calibrate_to_sim(data_fit, sim_fit)

    # ── Plot ──────────────────────────────────────────────────────────
    print("Plotting diagnostics...")
    plot_diagnostics(data_dedx, sim_dedx, cfg,
                     output=f"{args.output}_diagnostics.pdf")

    print("Plotting overlay...")
    plot_overlay(data_dedx, sim_dedx, data_fit, sim_fit, k, cfg,
                 output=args.output)

    print("\nDone.")


if __name__ == "__main__":
    main()