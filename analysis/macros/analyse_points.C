/*
 * analyse_points.C
 * ================
 * Extracts nPoints per track from experimental ROOT files using
 * TTree::Draw — no class dictionary or .so required.
 *
 * Run from the repository root:
 *     root -b -q macros/analyse_points.C
 *
 * Produces: npoints_dump.csv
 * Load in Python: hibeam.io.csv_loader.load_npoints("npoints_dump.csv")
 */

#include <iostream>
#include <fstream>
#include <vector>
#include <string>

void analyse_points() {

    std::vector<std::string> files = {
        "data/experimental/tracks_centroids_tpc_run_0006-sorted_t0_100.root",
        "data/experimental/tracks_centroids_tpc_run_0042-sorted_t0.root",
        "data/experimental/tracks_centroids_chamber-fec01-debug-clk0-ped0-thr0-g0-s0-p1-1008-2_sorted.root",
        "data/experimental/tracks_centroids_chamber-fec01-debug-clk0-ped0-thr0-g0-s0-p1-1008-4_sorted.root"
    };

    std::ofstream csv("npoints_dump.csv");
    csv << "file,event,nPoints\n";

    for (auto& fname : files) {
        std::string label = fname.substr(fname.rfind('/') + 1);
        label = label.substr(0, label.rfind('.'));

        TFile* f = TFile::Open(fname.c_str());
        if (!f || f->IsZombie()) {
            std::cerr << "Cannot open " << fname << "\n";
            continue;
        }
        TTree* tree = (TTree*)f->Get("trackingData");
        if (!tree) {
            std::cerr << "No trackingData in " << fname << "\n";
            f->Close(); continue;
        }

        tree->SetEstimate(tree->GetEntries() * 20);

        Long64_t n = tree->Draw("tracks.nPoints:Entry$", "", "goff");
        std::cout << label << ": " << n << " track rows\n";

        Double_t* npts   = tree->GetV1();
        Double_t* evnums = tree->GetV2();

        for (Long64_t i = 0; i < n; i++) {
            csv << label << ","
                << (Long64_t)evnums[i] << ","
                << (int)npts[i] << "\n";
        }

        f->Close();
    }

    csv.close();
    std::cout << "\nWritten: npoints_dump.csv\n";
}
