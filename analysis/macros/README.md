# ROOT C++ macros

These are legacy ROOT macros from the original code I made.  Only
`analyse_points.C` is still needed — the others have been fully
superseded by the Python plotting modules.

## `analyse_points.C`

Extracts `nPoints` per track via `TTree::Draw("tracks.nPoints:Entry$")`,
producing `npoints_dump.csv`.  This approach uses no class dictionary
and is sometimes faster than PyROOT for bulk extraction. Still, this is
a bit old fashioned compared to the rest of the stuff.

```bash
cd data/experimental/
root -b -q ../../macros/analyse_points.C
# produces: npoints_dump.csv
```

Then load in Python:
```python
from hibeam.io import csv_loader
df = csv_loader.load_npoints("npoints_dump.csv")
```

## Superseded macros (kept for reference only)

| Macro | Replaced by | Notes |
|---|---|---|
| `colourEloss.c` | `hibeam.plotting.displays.event_display_3d()` | 3D event display |
| `paths.c` | `hibeam.plotting.displays.event_display_3d()` | 3D trajectory plot |
| `saveGraph.c` | Not needed | Incomplete spare code |
