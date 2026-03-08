// ---------------------------------------------------------
// Definition of the VDetectorConstruction class
// Created by C.Rappold (c.rappold@gsi.de)
//----------------------------------------------------------

#ifndef VDetectorConstruction_h
#define VDetectorConstruction_h 1

#include "Config.hh"
#include "G4VPhysicalVolume.hh"
#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"

class G4VPhysicalVolume;

/// Detector construction class to define materials and geometry.

class VDetectorConstruction : public G4VUserDetectorConstruction
{
public:
  explicit VDetectorConstruction(Config& _par);
  virtual ~VDetectorConstruction();

public:
  virtual G4VPhysicalVolume* Construct() = 0;
  virtual void ConstructSDandField()     = 0;

  const std::vector<G4String>& GetNameDetectors() const { return NameDetectorsAll; };

  G4LogicalVolume* experimentalHall_logOutRoot;
  G4VPhysicalVolume* experimentalHall_physOutRoot;

protected:
  G4LogicalVolume* FindVolume(const G4String& name);
  Config& Par;
  std::vector<G4String> NameDetectorsSD;
  std::vector<G4String> NameDetectorsSample;
  std::vector<G4String> NameDetectorsAll;
};

#endif
