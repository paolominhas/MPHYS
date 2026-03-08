//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
// $Id: WasaDetectorConstruction.cc 77601 2013-11-26 17:08:44Z gcosmo $
//
/// \file WasaDetectorConstruction.cc
/// \brief Implementation of the WasaDetectorConstruction class

#include "WasaDetectorConstruction.hh"

#include "globals.hh"

#include "G4LogicalVolume.hh"
#include "G4Material.hh"
#include "G4NistManager.hh"
#include "G4PVPlacement.hh"
#include "G4PVReplica.hh"
#include "G4AutoDelete.hh"
#include "G4Colour.hh"
#include "G4GeometryManager.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4NistManager.hh"
#include "G4PhysicalConstants.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4VisAttributes.hh"

// VGM demo
#include "Geant4GM/volumes/Factory.h"
#include "RootGM/volumes/Factory.h"
#include "TGeoManager.h"
// end VGM demo

#include "G4SystemOfUnits.hh"

#include "G4AutoDelete.hh"
#include "G4ProductionCuts.hh"
#include "G4Region.hh"
#include "G4RegionStore.hh"
#include "G4SDManager.hh"
#include "SamplingD.hh"
#include "SensitiveD.hh"
#include "G4TransportationManager.hh"
#include "G4UserLimits.hh"

#include <iostream>

bool has_suffix(const std::string &str, const std::string &suffix)
{
    return str.size() >= suffix.size() &&
           str.compare(str.size() - suffix.size(), suffix.size(), suffix) == 0;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

WasaDetectorConstruction::WasaDetectorConstruction(Config& _par)
    : VDetectorConstruction(_par), experimentalHall_log(nullptr), experimentalHall_phys(nullptr)
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

WasaDetectorConstruction::~WasaDetectorConstruction() {}

G4VPhysicalVolume* WasaDetectorConstruction::Construct()
{
  //
  // Import geometry from Root or gdml
  //

  const std::string nameGeometry = Par.GetString("Geometry_Namefile");
  gGeoManager->Import(nameGeometry.c_str());
  const int checkOverlaps = Par.GetInt("CheckOverlaps");
  G4VPhysicalVolume* world;     

  if(has_suffix(nameGeometry,"gdml")){
	  G4GDMLParser parser;

	  if(checkOverlaps!=0) parser.SetOverlapCheck(true);
	  parser.Read(nameGeometry);

	  G4cout << std::endl;

	  const G4LogicalVolumeStore* lvs = G4LogicalVolumeStore::GetInstance();
	  std::vector<G4LogicalVolume*>::const_iterator lvciter;
	  for( lvciter = lvs->begin(); lvciter != lvs->end(); lvciter++ )
	  {
		  G4GDMLAuxListType auxInfo = parser.GetVolumeAuxiliaryInformation(*lvciter);

		  if (auxInfo.size()>0)
			  G4cout << "Auxiliary Information is found for Logical Volume :  "
				  << (*lvciter)->GetName() << G4endl;

		  print_aux(&auxInfo);
	  }

	  // now the 'global' auxiliary info
	  G4cout << std::endl;
	  G4cout << "Global auxiliary info:" << std::endl;
	  G4cout << std::endl;
	  print_aux(parser.GetAuxList());
	  G4cout << std::endl;

	  world = parser.GetWorldVolume();     
  }
  else if(has_suffix(nameGeometry,"root")){

	  // Import geometry from Root to VGM
	  RootGM::Factory rtFactory;
	  rtFactory.SetDebug(0);
	  rtFactory.Import(gGeoManager->GetTopNode());

	  // Export VGM geometry to Geant4
	  Geant4GM::Factory g4Factory;
	  g4Factory.SetDebug(0);
	  if(checkOverlaps!=0) g4Factory.SetSurfCheck(true);
	  rtFactory.Export(&g4Factory);

	  world = g4Factory.World();
  }
  experimentalHall_log  = world->GetLogicalVolume();
  experimentalHall_phys = world;

  const std::string allDetectors = Par.GetString("Detectors");
  std::istringstream streamDetectors{allDetectors};
  for(std::string detector{}; std::getline(streamDetectors, detector, ','); NameDetectorsSD.push_back(detector));
  NameDetectorsAll.insert(NameDetectorsAll.end(),NameDetectorsSD.begin(),NameDetectorsSD.end());

  const std::string allSamplingDetectors = Par.GetString("Sampling_Detectors");
  std::istringstream streamSamplingDetectors{allSamplingDetectors};
  for(std::string detector{}; std::getline(streamSamplingDetectors, detector, ','); NameDetectorsSample.push_back(detector));
  NameDetectorsAll.insert(NameDetectorsAll.end(),NameDetectorsSample.begin(),NameDetectorsSample.end());

  //G4VisAttributes VisDetectorSD(Blue);
  return world;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void WasaDetectorConstruction::DefineMaterials()
{
  // Dummy, as materials are imported via VGM
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* WasaDetectorConstruction::DefineVolumes()
{
  // Dummy, as geometry is imported via VGM

  return 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void WasaDetectorConstruction::ConstructSDandField()
{
  G4SDManager* SDman = G4SDManager::GetSDMpointer();

  for(auto& CurrentName : NameDetectorsSD)
    {
      G4LogicalVolume* Det = FindVolume(CurrentName);
      SD_Det* SD     = new SD_Det(CurrentName);
      // SD->Init();
      SDman->AddNewDetector(SD);
      Det->SetSensitiveDetector(SD);
    }

  std::cout << " Sensitive Detectors :" << std::endl;
  for(auto NameD : NameDetectorsSD)
    std::cout << NameD << std::endl;

  for(auto& CurrentName : NameDetectorsSample)
    {
      G4LogicalVolume* Det = FindVolume(CurrentName);
      Sample_Det* SD     = new Sample_Det(CurrentName);
      // SD->Init();
      SDman->AddNewDetector(SD);
      Det->SetSensitiveDetector(SD);
    }

  std::cout << " Sampling Detectors :" << std::endl;
  for(auto NameD : NameDetectorsSample)
    std::cout << NameD << std::endl;

  experimentalHall_log->SetUserLimits(new G4UserLimits(DBL_MAX, 2 * m, 10 * s, 0., 0.));
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* WasaDetectorConstruction::FindVolPhys(const G4String& name)
{
  G4PhysicalVolumeStore* pvStore = G4PhysicalVolumeStore::GetInstance();
  return pvStore->GetVolume(name);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void WasaDetectorConstruction::print_aux(const G4GDMLAuxListType* auxInfoList)
{
  G4String prepend="|";
  for(std::vector<G4GDMLAuxStructType>::const_iterator
      iaux = auxInfoList->begin(); iaux != auxInfoList->end(); iaux++ )
    {
      G4String str=iaux->type;
      G4String val=iaux->value;
      G4String unit=iaux->unit;

      G4cout << prepend << str << " : " << val  << " " << unit << G4endl;

      if (iaux->auxList) print_aux(iaux->auxList);
  }
  return;
}
