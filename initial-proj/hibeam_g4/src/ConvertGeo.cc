// ---------------------------------------------------------
// Implementation of the ConvertGeo class
// Created by C.Rappold (c.rappold@gsi.de)
//----------------------------------------------------------

#include "ConvertGeo.hh"

#include "Geant4GM/volumes/Factory.h"
#include "RootGM/volumes/Factory.h"
#include "TColor.h"
#include "TFile.h"
#include "TGeoManager.h"
#include "TROOT.h"

#include "G4SystemOfUnits.hh"
#include "G4UnitsTable.hh"

ConvertGeo::ConvertGeo(G4VPhysicalVolume* w) : physiWorld(w) {}
ConvertGeo::~ConvertGeo() {}

int ConvertGeo::Convert(const std::string& nameOut, const G4String& nameGeometry, const std::vector<G4String>& NameDetectorsSD, const Config& Config)
{
	// Import Geant4 geometry to VGM

	if(gGeoManager != nullptr)
	{
		std::cout << "!> ConvertGeo: gGeoManager not null !\n";
		gGeoManager->Print();
		TGeoManager* oldgeo1 = gGeoManager;
		gGeoManager          = nullptr;
		delete oldgeo1;
		oldgeo1 = nullptr;
	}

	Geant4GM::Factory g4Factory;
	g4Factory.Import(physiWorld);
	// where physiWorld is of G4VPhysicalVolume* type

	// Export VGM geometry to Root
	RootGM::Factory rtFactory;
	g4Factory.Export(&rtFactory);

	gGeoManager->CloseGeometry();

	TFile* f_outGeo = new TFile(nameOut.c_str(), "RECREATE");
	f_outGeo->cd();
	gGeoManager->Write();

	f_outGeo->WriteObjectAny(&NameDetectorsSD, "std::vector<std::string>", "nameDet");

	f_outGeo->Close();

	return 1;
}
