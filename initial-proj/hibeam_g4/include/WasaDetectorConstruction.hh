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
// $Id: WasaDetectorConstruction.hh 75215 2013-10-29 16:07:06Z gcosmo $
//
/// \file WasaDetectorConstruction.hh
/// \brief Definition of the WasaDetectorConstruction class

#ifndef WasaDetectorConstruction_h
#define WasaDetectorConstruction_h 1

#include "VDetectorConstruction.hh"
#include "G4Color.hh"
#include "Config.hh"
#include "G4VPhysicalVolume.hh"
#include "G4GDMLParser.hh"

class G4PVPlacement;
class G4VPhysicalVolume;
class G4Material;
class G4VSensitiveDetector;
class G4VisAttributes;

/// Detector construction class to define materials and geometry.

class WasaDetectorConstruction : public VDetectorConstruction
{
public:
  explicit WasaDetectorConstruction(Config& conf);
  virtual ~WasaDetectorConstruction();

  virtual G4VPhysicalVolume* Construct() override;
  virtual void ConstructSDandField() override;

private:
  G4LogicalVolume* experimentalHall_log;
  G4VPhysicalVolume* experimentalHall_phys;

  std::vector<G4PVPlacement*> AllPlacements;

  G4VPhysicalVolume* FindVolPhys(const G4String& name);
  void print_aux(const G4GDMLAuxListType* auxInfoList);
  // methods
  //
  void DefineMaterials();
  G4VPhysicalVolume* DefineVolumes();
};

#endif
