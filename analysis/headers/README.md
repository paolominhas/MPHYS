# Compiled C++ headers for experimental data

These files are produced by ROOT's `TFile::MakeProject` on the experimental
ROOT files.  They define the `TrackData`, `CentroidData`, and `tpcData`
classes needed to deserialise the `tracks` branch.

## Files

Copied from the existing `recovered_headers/` directory made with root:

```
headers/
├── TrackData.h                          # Track class definition
├── CentroidData.h                       # Centroid class definition
├── tpcData.h                            # TPC data container
├── recovered_headersProjectHeaders.h    # MakeProject master header
├── recovered_headersProjectSource.cxx   # MakeProject source
├── recovered_headersProjectInstances.h
├── recovered_headersProjectDict.cxx
├── recovered_headersProjectDict_rdict.pcm
├── recovered_headersLinkDef.h
├── MAKEP                                # Makefile from MakeProject
├── recovered_headersProjectSource.o     # Compiled object
├── recovered_headers.so                 # ← THE compiled shared library
└── LinkDef.h                            # (if present)
```

## No modifications needed

The headers are used exactly as-is.  The `exp_loader` module handles
the critical loading order (AddIncludePath → ProcessLine → gSystem.Load)
internally.

## If you need to regenerate them

If you get a new experimental ROOT file with a different class layout:

```python
import ROOT
f = ROOT.TFile.Open("your_new_file.root")
f.MakeProject("headers", "*", "RECREATE+")
f.Close()
```

Then compile:

```bash
cd headers/
make -f MAKEP
```

This produces a fresh `recovered_headers.so`.

## Troubleshooting

- **All-zero charge arrays**: The header loading order matters.
  `exp_loader` handles this, but if you call PyROOT manually, you must:
  1. `gInterpreter.AddIncludePath("headers/")`
  2. `gInterpreter.ProcessLine('#include "TrackData.h"')`
  3. `gSystem.Load("headers/recovered_headers.so")`

- **"CheckByteCount" errors on centroids branch**: This is a known
  issue with memberwise-serialised vectors.  `exp_loader` reads
  the `tracks` branch via `SetBranchAddress` instead, which avoids
  the problem entirely.
