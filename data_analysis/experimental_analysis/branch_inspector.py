"""
Quick branch inspector — run this first to see the exact branch paths
in your trackingData tree so the analysis script uses the right names.
Usage:  python inspect_branches.py
"""
import uproot
import glob
import os

INPUT_DIR = "../experimental_data/"

files = sorted(glob.glob(os.path.join(INPUT_DIR, "*.root")))
if not files:
    print(f"No .root files found in {INPUT_DIR}")
    exit()

# Just inspect the first file
path = files[0]
print(f"Inspecting: {os.path.basename(path)}\n")

with uproot.open(path) as f:
    print("=== Top-level keys ===")
    for k in f.keys():
        print(f"  {k}")

    print("\n=== trackingData branches ===")
    try:
        tree = f["trackingData"]
        for name, branch in tree.items():
            try:
                dtype = branch.interpretation
            except Exception:
                dtype = "?"
            print(f"  {name:<60}  {dtype}")
    except Exception as e:
        print(f"  Could not open trackingData: {e}")