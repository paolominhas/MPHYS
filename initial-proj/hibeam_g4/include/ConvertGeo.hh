// ---------------------------------------------------------
// Definition of the ConvertGeo class
// Created by C.Rappold (c.rappold@gsi.de)
//----------------------------------------------------------

#ifndef ConvertGeo_h
#define ConvertGeo_h 1

#include "G4VPhysicalVolume.hh"
#include "Config.hh"

class ConvertGeo
{
  G4VPhysicalVolume* physiWorld;

public:
  ConvertGeo(G4VPhysicalVolume* w);
  ~ConvertGeo();
  int Convert(const std::string& nameOut, const G4String& nameGeometry, const std::vector<G4String>& NameDetectorsSD, const Config& Par);
};

#endif
