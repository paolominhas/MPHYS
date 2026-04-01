"""
hibeam.io — Data loading modules
=================================
One loader per source format.  Each returns a standardised dictionary
so downstream physics/plotting code is source-agnostic.

Modules
-------
sim_loader  — Geant4 simulation ROOT files (uproot).
exp_loader  — Experimental TPC ROOT files (PyROOT + compiled headers).
seg_loader  — Segmentation study ROOT files (uproot).
csv_loader  — Auxiliary CSV / text data (pandas).
"""

from hibeam.io.sim_loader import load as load_sim
from hibeam.io.exp_loader import load as load_exp
from hibeam.io.seg_loader import load as load_seg
from hibeam.io.csv_loader import load_csv
