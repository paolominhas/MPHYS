import ROOT

# Point this at one of your real data files
DATA_FILE = "../experimental_data/proton250.root"

f = ROOT.TFile.Open(DATA_FILE)
if not f or f.IsZombie():
    print("Could not open file!")
else:
    # "new++" forces a full recompile and relink
    f.MakeProject("recovered_headers", "*", "new++")
    f.Close()
    print("Done — recovered_headers.so has been regenerated.")