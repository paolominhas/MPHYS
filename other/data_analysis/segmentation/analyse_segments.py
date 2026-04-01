#!/usr/bin/env python3
"""
analyse_segmentation.py
=======================
Quantitative analysis of ProtoTPC detector segmentation impact on energy
deposition. Uses TTree::Draw to fill histograms — bypasses the missing
ProtoTPCHit/HRDSciHit class dictionaries entirely.

Branches read (all via Draw, no dictionary needed):
  SumEdep          — total Edep per event [MeV]
  Edep             — per-hit Edep (vector, Draw flattens automatically)
  Layer            — per-hit segment index
  nHits            — number of hits per event

Figures of merit:
  1. dE/dx resolution R = sigma/mean        (Blum, Riegler & Rolandi 2008)
  2. Landau MPV and width xi                (Landau 1944; Moyal 1955)
  3. Truncated mean resolution              (Allison & Cobb 1980)
  4. Per-layer Edep uniformity              (sigma/mean across layers)
  5. Mean hits per event vs nSections

Usage:
    python3 analyse_segmentation.py --indir ./ana_output_5_10_15_20_proton
"""

import argparse
import re
import sys
import ctypes
from pathlib import Path

parser = argparse.ArgumentParser(description=__doc__,
    formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument("--indir",    default="./ana_output_krakow",
    help="Directory containing ana_geom_nSec*.root files")
parser.add_argument("--outdir",   default="./seg_analysis",
    help="Output directory for plots and summary (default: ./seg_analysis)")
parser.add_argument("--pattern",  default="ana_geom_*.root",
    help="Glob pattern for input files (default: ana_geom_*.root)")
parser.add_argument("--trunc",    default=0.70, type=float,
    help="Truncation fraction for truncated mean (default: 0.70)")
parser.add_argument("--treename", default="hibeam",
    help="TTree name (default: hibeam)")
args = parser.parse_args()

try:
    import ROOT
    ROOT.gROOT.SetBatch(True)
    ROOT.gErrorIgnoreLevel = ROOT.kWarning
    ROOT.gStyle.SetOptStat(1111)
    ROOT.gStyle.SetOptFit(1111)
except ImportError:
    sys.exit("ERROR: PyROOT not available. Source thisroot.sh first.")

# ---------- helpers ----------------------------------------------------------
def extract_nsections(fname):
    m = re.search(r'[nN][sS]ec0*(\d+)', fname)
    return int(m.group(1)) if m else None

def extract_axis(fname):
    m = re.search(r'[aA]xis(\d+)', fname)
    return int(m.group(1)) if m else None

def draw_to_hist(tree, expr, cut, name, title, nbins, xlo, xhi):
    """
    Use TTree::Draw to fill a histogram — works without class dictionaries.
    Returns TH1D detached from the current directory.
    """
    h = ROOT.TH1D(name, title, nbins, xlo, xhi)
    h.SetDirectory(0)
    tree.Draw(f"{expr}>>{name}", cut, "goff")
    # TTree::Draw fills the named histo in gDirectory; retrieve it
    filled = ROOT.gDirectory.Get(name)
    if filled and filled != h:
        h = filled.Clone(name + "_clone")
        h.SetDirectory(0)
    return h

def probe_range(tree, expr, cut=""):
    """Quick Draw to find min/max of an expression."""
    tree.Draw(expr, cut, "goff")
    h = ROOT.gPad.GetPrimitive("htemp") if ROOT.gPad else None
    # Safer: use a wide temporary histogram and read back axis range
    tmp = ROOT.TH1D("_probe_tmp", "", 500, -1e9, 1e9)
    tmp.SetDirectory(0)
    tree.Draw(f"{expr}>>_probe_tmp", cut, "goff")
    if tmp.GetEntries() == 0:
        return None, None
    # Find actual filled range
    lo = tmp.GetXaxis().GetBinLowEdge(tmp.FindFirstBinAbove(0))
    hi = tmp.GetXaxis().GetBinUpEdge(tmp.FindLastBinAbove(0))
    ROOT.gDirectory.Delete("_probe_tmp;*")
    return lo, hi

def fit_landau(hist):
    if not hist or hist.GetEntries() < 20:
        return None
    hist.Fit("landau", "Q0")
    f = hist.GetFunction("landau")
    if not f:
        return None
    return {
        "mpv"       : f.GetParameter(1),
        "mpv_err"   : f.GetParError(1),
        "xi"        : f.GetParameter(2),
        "xi_err"    : f.GetParError(2),
        "chi2ndf"   : f.GetChisquare() / max(f.GetNDF(), 1),
    }

def fit_gauss(hist):
    if not hist or hist.GetEntries() < 10:
        return None
    hist.Fit("gaus", "Q0")
    f = hist.GetFunction("gaus")
    if not f:
        return None
    mean  = f.GetParameter(1)
    sigma = f.GetParameter(2)
    return {
        "mean"      : mean,
        "sigma"     : sigma,
        "resolution": abs(sigma / mean) if mean != 0 else 999,
        "chi2ndf"   : f.GetChisquare() / max(f.GetNDF(), 1),
    }

def make_graph(x, y, ye=None):
    n = len(x)
    xarr = (ctypes.c_double*n)(*[float(v) for v in x])
    yarr = (ctypes.c_double*n)(*[float(v) for v in y])
    if ye:
        earr = (ctypes.c_double*n)(*[float(v) for v in ye])
        zarr = (ctypes.c_double*n)(*([0.0]*n))
        return ROOT.TGraphErrors(n, xarr, yarr, zarr, earr)
    return ROOT.TGraph(n, xarr, yarr)

# ---------- collect files ----------------------------------------------------
indir  = Path(args.indir)
outdir = Path(args.outdir)
outdir.mkdir(parents=True, exist_ok=True)

files = sorted(indir.glob(args.pattern))
if not files:
    sys.exit(f"ERROR: No files matching '{args.pattern}' in {indir}")

print(f"\nFound {len(files)} file(s) in {indir}")
print(f"Output    : {outdir}\n")

results = []
canvas  = ROOT.TCanvas("c1", "c1", 1000, 700)
ROOT.gStyle.SetPadGridX(True)
ROOT.gStyle.SetPadGridY(True)

for fpath in files:
    nsec = extract_nsections(fpath.name)
    axis = extract_axis(fpath.name)
    if nsec is None:
        print(f"  [SKIP] Cannot parse nSections from: {fpath.name}")
        continue

    print(f"  Processing: {fpath.name}  (nSections={nsec}, axis={axis})")

    rfile = ROOT.TFile.Open(str(fpath))
    if not rfile or rfile.IsZombie():
        print(f"    ERROR: Cannot open file.")
        continue

    tree = rfile.Get(args.treename)
    if not tree:
        print(f"    ERROR: TTree '{args.treename}' not found.")
        rfile.Close()
        continue

    n_events = tree.GetEntries()
    print(f"    Events: {n_events}")

    tag = f"nSec{nsec:03d}_ax{axis}"

    # ------------------------------------------------------------------
    # 1. SumEdep per event  (scalar branch — Draw works fine)
    # ------------------------------------------------------------------
    lo, hi = probe_range(tree, "SumEdep", "SumEdep>0")
    if lo is None:
        print(f"    WARNING: SumEdep is all zero — particle never reached ProtoTPC.")
        print(f"             Check that --proto_tpc was set and geometry is correct.")
        rfile.Close()
        continue

    # Add 20% headroom
    xmax_sum = hi * 1.2

    h_sum = draw_to_hist(tree, "SumEdep", "SumEdep>0",
        f"hSum_{tag}",
        f"Total E_{{dep}} per event  |  nSections={nsec};"
        f"E_{{dep}} [MeV];Events",
        80, 0, xmax_sum)

    # ------------------------------------------------------------------
    # 2. Per-hit Edep  (vector branch — TTree::Draw flattens the vector)
    # ------------------------------------------------------------------
    lo2, hi2 = probe_range(tree, "Edep", "Edep>0")
    if lo2 is not None:
        h_hit = draw_to_hist(tree, "Edep", "Edep>0",
            f"hHit_{tag}",
            f"Per-hit E_{{dep}}  |  nSections={nsec};"
            f"E_{{dep}} [MeV];Hits",
            80, 0, hi2 * 1.2)
    else:
        h_hit = None

    # ------------------------------------------------------------------
    # 3. nHits per event
    # ------------------------------------------------------------------
    h_nhits = draw_to_hist(tree, "nHits", "",
        f"hNhits_{tag}",
        f"nHits per event  |  nSections={nsec};"
        f"nHits;Events",
        50, 0, 50)

    mean_nhits = h_nhits.GetMean() if h_nhits.GetEntries() > 0 else 0

    # ------------------------------------------------------------------
    # 4. Truncated mean per event
    #    Use Draw with a quantile cut: keep events below the trunc-th
    #    quantile of SumEdep as an approximation of truncated mean method
    # ------------------------------------------------------------------
    # Find the truncation threshold from quantile
    q_val = ctypes.c_double(0)
    q_prob = ctypes.c_double(args.trunc)
    h_sum.GetQuantiles(1, q_val, q_prob)
    trunc_threshold = q_val.value

    cut_trunc = f"SumEdep>0 && SumEdep<={trunc_threshold}"
    h_trunc = draw_to_hist(tree, "SumEdep", cut_trunc,
        f"hTrunc_{tag}",
        f"{int(args.trunc*100)}% Truncated E_{{dep}}  |  nSections={nsec};"
        f"E_{{dep}} [MeV];Events",
        60, 0, trunc_threshold * 1.1)

    # ------------------------------------------------------------------
    # 5. Per-layer mean Edep  (Draw with Layer as Y, Edep as X -> profile)
    # ------------------------------------------------------------------
    lo_l, hi_l = probe_range(tree, "Layer", "Edep>0")
    n_layers = int(hi_l) + 1 if hi_l is not None else nsec

    h_layer_prof = ROOT.TProfile(f"hLayerProf_{tag}",
        f"Mean E_{{dep}} per layer  |  nSections={nsec};"
        f"Layer index;Mean E_{{dep}} [MeV]",
        n_layers, -0.5, n_layers - 0.5)
    h_layer_prof.SetDirectory(0)
    tree.Draw(f"Edep:Layer>>{f'hLayerProf_{tag}'}", "Edep>0", "goff prof")

    # Layer uniformity: RMS of per-layer means / grand mean
    layer_means = [h_layer_prof.GetBinContent(b)
                   for b in range(1, h_layer_prof.GetNbinsX()+1)
                   if h_layer_prof.GetBinContent(b) > 0]
    if len(layer_means) > 1:
        gm = sum(layer_means)/len(layer_means)
        rms_layers = (sum((v-gm)**2 for v in layer_means)/len(layer_means))**0.5
        layer_uniformity = rms_layers / gm if gm > 0 else None
    else:
        layer_uniformity = None

    # ------------------------------------------------------------------
    # Fits
    # ------------------------------------------------------------------
    landau_sum  = fit_landau(h_sum)
    gauss_trunc = fit_gauss(h_trunc)

    full_mean = h_sum.GetMean()
    full_rms  = h_sum.GetRMS()
    full_res  = full_rms / full_mean if full_mean > 0 else 999

    row = {
        "file"            : fpath.name,
        "nsections"       : nsec,
        "axis"            : axis,
        "n_events"        : n_events,
        "mean_nhits"      : mean_nhits,
        "full_mean"       : full_mean,
        "full_rms"        : full_rms,
        "full_res"        : full_res,
        "n_active_layers" : len(layer_means),
        "layer_uniformity": layer_uniformity,
    }
    if landau_sum:
        row.update({
            "landau_mpv"    : landau_sum["mpv"],
            "landau_mpv_err": landau_sum["mpv_err"],
            "landau_xi"     : landau_sum["xi"],
            "landau_xi_err" : landau_sum["xi_err"],
            "landau_chi2ndf": landau_sum["chi2ndf"],
        })
    if gauss_trunc:
        row.update({
            "trunc_mean"    : gauss_trunc["mean"],
            "trunc_sigma"   : gauss_trunc["sigma"],
            "trunc_res"     : gauss_trunc["resolution"],
            "trunc_chi2ndf" : gauss_trunc["chi2ndf"],
        })

    results.append(row)

    # ------------------------------------------------------------------
    # Per-file 4-panel plot
    # ------------------------------------------------------------------
    canvas.Clear()
    canvas.Divide(2, 2)

    # Panel 1: SumEdep + Landau
    canvas.cd(1); ROOT.gPad.SetLeftMargin(0.14)
    h_sum.SetLineColor(ROOT.kBlue+1); h_sum.SetLineWidth(2)
    h_sum.Draw("HIST")
    if landau_sum and h_sum.GetFunction("landau"):
        h_sum.GetFunction("landau").SetLineColor(ROOT.kOrange+2)
        h_sum.GetFunction("landau").SetLineWidth(2)
        h_sum.GetFunction("landau").Draw("SAME")

    # Panel 2: Per-hit Edep + Landau
    canvas.cd(2); ROOT.gPad.SetLeftMargin(0.14)
    if h_hit and h_hit.GetEntries() > 0:
        lf_hit = fit_landau(h_hit)
        h_hit.SetLineColor(ROOT.kGreen+2); h_hit.SetLineWidth(2)
        h_hit.Draw("HIST")
        if lf_hit and h_hit.GetFunction("landau"):
            h_hit.GetFunction("landau").SetLineColor(ROOT.kOrange+2)
            h_hit.GetFunction("landau").SetLineWidth(2)
            h_hit.GetFunction("landau").Draw("SAME")

    # Panel 3: Truncated distribution + Gaussian
    canvas.cd(3); ROOT.gPad.SetLeftMargin(0.14)
    h_trunc.SetLineColor(ROOT.kRed+1); h_trunc.SetLineWidth(2)
    h_trunc.Draw("HIST")
    if gauss_trunc and h_trunc.GetFunction("gaus"):
        h_trunc.GetFunction("gaus").SetLineColor(ROOT.kOrange+2)
        h_trunc.GetFunction("gaus").SetLineWidth(2)
        h_trunc.GetFunction("gaus").Draw("SAME")

    # Panel 4: Per-layer mean Edep profile
    canvas.cd(4); ROOT.gPad.SetLeftMargin(0.14)
    h_layer_prof.SetLineColor(ROOT.kViolet+1)
    h_layer_prof.SetFillColorAlpha(ROOT.kViolet, 0.25)
    h_layer_prof.SetLineWidth(2)
    h_layer_prof.Draw("HISTE")

    plotfile = str(outdir / f"edep_nSec{nsec:03d}_axis{axis}.pdf")
    canvas.SaveAs(plotfile)
    print(f"    Saved: {plotfile}")

    rfile.Close()

if not results:
    sys.exit("\nERROR: No results collected. "
             "SumEdep was zero in all files.\n"
             "Possible causes:\n"
             "  - ProtoTPC not set as sensitive detector in config.par "
             "(add 'Detectors = ProtoTPC' or 'Detectors = ProtoTPC_0, ...')\n"
             "  - Particles not reaching the ProtoTPC — check geometry placement\n"
             "  - Wrong geometry file used in simulation")

# ---------- summary comparison plots ----------------------------------------
results.sort(key=lambda r: r["nsections"])

nsec_x      = [r["nsections"]                          for r in results]
full_res_y  = [r["full_res"]*100                       for r in results]
trunc_res_y = [r.get("trunc_res",0)*100                for r in results]
mpv_y       = [r.get("landau_mpv",0)                   for r in results]
mpv_ey      = [r.get("landau_mpv_err",0)               for r in results]
xi_y        = [r.get("landau_xi",0)                    for r in results]
xi_ey       = [r.get("landau_xi_err",0)                for r in results]
nhits_y     = [r["mean_nhits"]                         for r in results]
unifo_y     = [(r["layer_uniformity"]*100
                if r.get("layer_uniformity") else 0)   for r in results]

c2 = ROOT.TCanvas("c_summary","Summary",1200,900)
c2.Divide(2,3)

def style_graph(g, color, marker, title, xtitle, ytitle):
    g.SetTitle(title)
    g.SetMarkerStyle(marker); g.SetMarkerSize(1.2)
    g.SetMarkerColor(color);  g.SetLineColor(color); g.SetLineWidth(2)
    g.Draw("ALP")
    g.GetXaxis().SetTitle(xtitle)
    g.GetYaxis().SetTitle(ytitle)

c2.cd(1); ROOT.gPad.SetGrid()
style_graph(make_graph(nsec_x, full_res_y),
    ROOT.kBlue+1, 21,
    "Full dist. resolution #sigma/#mu vs nSections",
    "nSections", "Resolution #sigma/#mu [%]")

c2.cd(2); ROOT.gPad.SetGrid()
style_graph(make_graph(nsec_x, trunc_res_y),
    ROOT.kRed+1, 22,
    f"{int(args.trunc*100)}% Trunc. resolution vs nSections",
    "nSections", "Trunc. resolution [%]")

c2.cd(3); ROOT.gPad.SetGrid()
style_graph(make_graph(nsec_x, mpv_y, mpv_ey),
    ROOT.kGreen+2, 20,
    "Landau MPV vs nSections",
    "nSections", "Landau MPV [MeV]")

c2.cd(4); ROOT.gPad.SetGrid()
style_graph(make_graph(nsec_x, xi_y, xi_ey),
    ROOT.kOrange+2, 23,
    "Landau #xi (width) vs nSections",
    "nSections", "Landau #xi [MeV]")

c2.cd(5); ROOT.gPad.SetGrid()
style_graph(make_graph(nsec_x, nhits_y),
    ROOT.kViolet+1, 21,
    "Mean hits per event vs nSections",
    "nSections", "Mean nHits / event")

c2.cd(6); ROOT.gPad.SetGrid()
style_graph(make_graph(nsec_x, unifo_y),
    ROOT.kCyan+2, 22,
    "Layer-to-layer uniformity vs nSections",
    "nSections", "Layer RMS/mean [%]")

summary_plot = str(outdir / "summary_vs_nsections.pdf")
c2.SaveAs(summary_plot)
print(f"\n  Saved: {summary_plot}")

# ---------- text summary -----------------------------------------------------
summary_txt = outdir / "summary.txt"
with open(summary_txt, "w") as f:
    sep = "=" * 115
    f.write(sep + "\n")
    f.write("ProtoTPC Segmentation Impact Analysis\n")
    f.write("References: Blum et al. (2008); Allison & Cobb (1980); ALICE TDR (2000)\n")
    f.write(sep + "\n\n")
    hdr = (f"{'nSec':>5}  {'ax':>3}  {'events':>7}  {'mean_nHits':>11}  "
           f"{'full_mean':>10}  {'full_res%':>10}  {'trunc_res%':>10}  "
           f"{'L_MPV':>9}  {'L_xi':>8}  {'chi2/ndf':>9}  {'lyr_unif%':>10}")
    f.write(hdr + "\n" + "-"*115 + "\n")
    for r in results:
        f.write(
            f"{r['nsections']:>5}  {str(r['axis']):>3}  {r['n_events']:>7}  "
            f"{r['mean_nhits']:>11.2f}  {r['full_mean']:>10.5f}  "
            f"{r['full_res']*100:>10.3f}  "
            f"{r.get('trunc_res',0)*100:>10.3f}  "
            f"{r.get('landau_mpv',0):>9.5f}  "
            f"{r.get('landau_xi',0):>8.5f}  "
            f"{r.get('landau_chi2ndf',0):>9.3f}  "
            f"{(r['layer_uniformity']*100 if r.get('layer_uniformity') else 0):>10.3f}\n")
    f.write("\n" + sep + "\n")
    f.write("Stability check (theory: all metrics should be flat vs nSections)\n")
    f.write(sep + "\n")

    def stability(label, vals):
        import statistics
        clean = [v for v in vals if v > 0]
        if len(clean) < 2: return
        m = statistics.mean(clean)
        s = statistics.stdev(clean)
        vp = s/m*100
        flag = "STABLE (<5%)" if vp < 5 else "VARIES (>5%) — investigate"
        f.write(f"  {label:<38}: mean={m:.4f}  std={s:.4f}  "
                f"variation={vp:.2f}%  [{flag}]\n")

    stability("Full resolution sigma/mu [%]",   full_res_y)
    stability("Truncated resolution [%]",        trunc_res_y)
    stability("Landau MPV [MeV]",                mpv_y)
    stability("Landau xi [MeV]",                 xi_y)
    stability("Layer uniformity [%]",            unifo_y)

print(f"  Saved: {summary_txt}")
print("\nDone.\n")