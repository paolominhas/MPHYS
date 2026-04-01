// export.c
#include <iostream>
#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TGraph2D.h"
#include "TStyle.h"


void export_geo() {
    // Open file
    TFile *f = TFile::Open("krakow_setup.root");
    
    // retrieve the TGeoManager object(s)
    TGeoManager *geo = (TGeoManager*)f->Get("geometry"); 
    
    if (geo) {
        // Export to .obj format
        geo->Export("detector.obj");
        cout << "Export successful!" << endl;
    } else {
        cout << "Could not find TGeoManager object." << endl;
    }
}