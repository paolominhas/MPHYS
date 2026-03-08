#include <iostream>
#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TGraph2D.h"
#include "TStyle.h"

void analyze_scatter() {
    TFile *file = TFile::Open("file.root", "READ");
    if (!file || file->IsZombie()) {
        std::cerr << "Error: Could not open file!" << std::endl;
        return;
    }

    TTree *tree = (TTree*)file->Get("hibeam;1");
    if (!tree) {
        std::cerr << "Error: Could not find TTree!" << std::endl;
        return;
    }

    // variables for data branch
    Double_t x_val, y_val, z_val;

    // Branch Addresses
    tree->SetBranchAddress("Pos_X", &x_val);
    tree->SetBranchAddress("Pos_Y", &y_val);
    tree->SetBranchAddress("Pos_Z", &z_val);
    
    // TGraph2D
    TGraph2D *graph = new TGraph2D();
    graph->SetTitle("3D Scatter Plot of x, y, z; X Axis; Y Axis; Z Axis");
    Long64_t nEntries = tree->GetEntries();
    std::cout << "Processing " << nEntries << " events..." << std::endl;

    int pointIndex = 0;
    for (Long64_t i = 0; i < nEntries; i++) {
        tree->GetEntry(i);
        graph->SetPoint(pointIndex, x_val, y_val, z_val);
        pointIndex++;
    }

    TCanvas *c1 = new TCanvas("c1", "3D Analysis", 800, 600);
    
    // 20 is a filled circle marker style and colour kBlue
    graph->SetMarkerStyle(20);
    graph->SetMarkerSize(0.6);
    graph->SetMarkerColor(kBlue);

    c1->Update();
    c1->Draw();
}