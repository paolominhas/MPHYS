#include <iostream>
#include <vector>
#include "TFile.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "TCanvas.h"
#include "TPolyLine3D.h"
#include "TH3F.h"

void paths() {
    TFile *f = new TFile("../hibeam_g4/outputs/MCPLData.root");
    TTreeReader reader("hibeam;1", f);
    TTreeReaderValue<std::vector<double>> posX(reader, "ProtoTPC_GlobalPosition_X");
    TTreeReaderValue<std::vector<double>> posY(reader, "ProtoTPC_GlobalPosition_Y");
    TTreeReaderValue<std::vector<double>> posZ(reader, "ProtoTPC_GlobalPosition_Z");
    TCanvas *c1 = new TCanvas("c1", "Reconstructed Paths", 800, 600);
    TH3F *axes = new TH3F("axes", "Proto-TPC Trajectories;X;Y;Z",   10, -50, 50,      10, -50, 50,       10, -50, 50); 
    axes->SetStats(0); 
    axes->Draw();      
    int counter = 0;
    while (reader.Next()) {
        // Just 10 tracks here for now can be gotten rid of for all
        if (counter >= 1000) break;               
        counter++;
        int nHits = (*posX).size();                                 // access data by *dereferencing* the reader variable (*posX): gives the vector for current event
        if (nHits == 0) continue;                                   // if empty  skip=
        TPolyLine3D *path = new TPolyLine3D(nHits);

        for (int j = 0; j < nHits; ++j) {
            path->SetPoint(j, (*posX)[j], (*posY)[j], (*posZ)[j]);  // fill line with points from vectors
        }

        path->SetLineColor(kRed);   // red line
        path->SetLineWidth(2);      // thickness
        path->Draw("SAME");         // draw on top
    }
}