#include <iostream>
#include <vector>
#include <numeric> // For std::accumulate
#include "TFile.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "TCanvas.h"
#include "TPolyLine3D.h"
#include "TPolyMarker3D.h"
#include "TH3F.h"
#include "TH1F.h"

void colourEloss() {
    TFile *f = new TFile("../hibeam_g4/outputs/MuonDataBronowice.root");
    TTreeReader reader("hibeam;1", f);
    TTreeReaderValue<std::vector<double>> posX(reader, "ProtoTPC_Position_X");
    TTreeReaderValue<std::vector<double>> posY(reader, "ProtoTPC_Position_Y");
    TTreeReaderValue<std::vector<double>> posZ(reader, "ProtoTPC_Position_Z");
    TTreeReaderValue<std::vector<double>> edep(reader, "ProtoTPC_EDep");
    TCanvas *c1 = new TCanvas("c1", "3D Tracks", 800, 600);
    TCanvas *c2 = new TCanvas("c2", "Energy Analysis", 800, 600);
    // 3D View Frame for the tracks
    c1->cd();
    TH3F *axes = new TH3F("axes", "TPC Tracks;X;Y;Z", 10, -50, 50, 10, -50, 50, 10, -50, 50);
    axes->SetStats(0);
    axes->Draw();

    // Total Energy Loss per Event histogram setup
    TH1F *hTotalEnergy = new TH1F("hTotalEnergy", "Total Energy Deposited in TPC;Energy (MeV);Count", 100, 0, 100); 
    int counter = 0;
    
    // Loop over Events
    while (reader.Next()) {
        int nHits = (*posX).size();
        if (nHits == 0) continue;

        // tracks
        c1->cd();
        TPolyMarker3D *markers = new TPolyMarker3D(nHits);
        
        // energy calc
        double total_event_energy = 0.0;
        for (int j = 0; j < nHits; ++j) {
            markers->SetPoint(j, (*posX)[j], (*posY)[j], (*posZ)[j]);
            
            // energy for this specific hit
            double energy_hit = (*edep)[counter];
            total_event_energy += energy_hit;
            // print energy of specific hits to console
            std::cout << "Event " << counter << " Hit " << j << ": " << energy_hit << " MeV" << std::endl;
        }

        markers->SetMarkerStyle(20);
        markers->Draw("SAME");

        // Fill histogram
        hTotalEnergy->Fill(total_event_energy);
        counter += 1;
    }
    c2->cd();
    hTotalEnergy->Draw();
}