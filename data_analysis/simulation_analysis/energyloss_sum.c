#include "TFile.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "TCanvas.h"
#include "TH1F.h"
#include <vector>

#include "TFile.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "TH1F.h"
#include <vector>
#include <fstream>  // Include for file writing at the bottom

void energyloss_sum() {
    TFile *f = new TFile("../final_sim_data/mcpl_long_muon.root", "READ");
    TTreeReader reader("hibeam", f);
    TTreeReaderValue<std::vector<double>> edep(reader, "ProtoTPC.Edep");

    // Open  CSV file
    std::ofstream csvFile("muon_long_output.csv");
    // Write Header (good for pandas df to be used in a second)
    csvFile << "total_energy_gev" << std::endl; 

    int eventCounter = 0; // Put this BEFORE the while loop

    while (reader.Next()) {
        eventCounter++;
        if (eventCounter % 100000 == 0) {
            printf("Processed %d events...\n", eventCounter);
        }
        double eventSum = 0.0;
        for (double hit : *edep) {
            eventSum += hit;
        }

        

        // Apply the filter (which is now heavily reduced to include that noise just for the moment)
        if (eventSum > 0.001) {
            // Write raw value of sum to CSV
            csvFile << eventSum << std::endl;
        }
    }
    
    csvFile.close(); // Close file remember!
    printf("Data saved to energy_output.csv\n");
}