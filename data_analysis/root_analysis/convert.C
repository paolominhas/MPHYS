#include <vector>
#include "TFile.h"
#include "TTree.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "TTreeReaderArray.h"

// STRUCTURAL ASSUMPTION:
// I am assuming the 'centroids' branch is a collection (std::vector or TClonesArray)
// and has members like .x, .y, .z, and .q (charge/energy).
// You may need to adjust "centroids.x" to whatever the actual member name is (e.g. "centroids.fX")
// Check this by running tree->Print() in your terminal first!

void convertToFlat(const char* inputFileName, const char* outputFileName = "flat_data.root") {
    
    // 1. Open Data
    TFile *inFile = TFile::Open(inputFileName);
    if (!inFile || inFile->IsZombie()) {
        printf("Error: Cannot open file %s\n", inputFileName);
        return;
    }

    // Access the tree (adjust name if strictly 'trackingData' without ;6)
    TTreeReader reader("trackingData;6", inFile); 

    // 2. Setup Input Readers (The Complex Data)
    // TPC Raw Hits (if you want raw pad data)
    TTreeReaderArray<Int_t>  raw_row(reader, "tpc.row");
    TTreeReaderArray<Int_t>  raw_col(reader, "tpc.column");
    TTreeReaderArray<Double_t> raw_val(reader, "tpc.val"); // Signal
    TTreeReaderArray<Double_t> raw_time(reader, "tpc.timestamp");

    // Reconstructed Centroids (The high level data)
    // NOTE: You must check exact variable names using tree->Print()! 
    // These are placeholders based on standard naming conventions.
    TTreeReaderArray<Double_t> cent_x(reader, "centroids.x"); 
    TTreeReaderArray<Double_t> cent_y(reader, "centroids.y");
    TTreeReaderArray<Double_t> cent_z(reader, "centroids.z");
    TTreeReaderArray<Double_t> cent_q(reader, "centroids.q"); // Charge/Energy

    // 3. Setup Output Tree (The Simple "Flat" Data)
    TFile *outFile = new TFile(outputFileName, "RECREATE");
    TTree *outTree = new TTree("TPCTree", "Flattened TPC Data for Python");

    // Variables for the new tree
    // We use std::vector because the number of hits changes every event
    std::vector<float> out_x, out_y, out_z, out_q;
    std::vector<int> out_row, out_col;
    std::vector<float> out_val, out_time;
    int event_id = 0;

    // Create Branches
    outTree->Branch("eventID", &event_id, "eventID/I");
    outTree->Branch("x", &out_x);
    outTree->Branch("y", &out_y);
    outTree->Branch("z", &out_z);
    outTree->Branch("q", &out_q); // Charge/Energy deposit
    
    // Raw branches (optional, good for debugging)
    outTree->Branch("raw_row", &out_row);
    outTree->Branch("raw_col", &out_col);
    outTree->Branch("raw_val", &out_val);

    // 4. Loop and Fill
    printf("Starting conversion...\n");
    while (reader.Next()) {
        // Clear vectors for new event
        out_x.clear(); out_y.clear(); out_z.clear(); out_q.clear();
        out_row.clear(); out_col.clear(); out_val.clear();

        // -- Process Reconstructed Centroids --
        for (int i = 0; i < cent_x.GetSize(); ++i) {
            out_x.push_back( static_cast<float>(cent_x[i]) );
            out_y.push_back( static_cast<float>(cent_y[i]) );
            out_z.push_back( static_cast<float>(cent_z[i]) );
            out_q.push_back( static_cast<float>(cent_q[i]) );
        }

        // -- Process Raw Hits (Optional) -- 
        for (int i = 0; i < raw_row.GetSize(); ++i) {
            out_row.push_back(raw_row[i]);
            out_col.push_back(raw_col[i]);
            out_val.push_back( static_cast<float>(raw_val[i]) );
        }

        outTree->Fill();
        event_id++;
        if (event_id % 100 == 0) printf("Processed event %d\r", event_id);
    }

    // 5. Write and Close
    outFile->Write();
    outFile->Close();
    inFile->Close();
    printf("\nConversion complete! Saved to %s\n", outputFileName);
}