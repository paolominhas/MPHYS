import ROOT

DATA_FILE = "../experimental_data/tracks_centroids_tpc_run_0006-sorted_t0_100.root"

f = ROOT.TFile.Open(DATA_FILE)
if not f or f.IsZombie():
    print("Could not open file!")
else:
    f.MakeProject("recovered_headers", "*", "new++")
    f.Close()
    print("Done.")