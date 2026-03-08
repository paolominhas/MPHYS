// -------------------------------------------------
// Definition of the SD_Det class
// Created by C.Rappold (c.rappold@gsi.de)
//--------------------------------------------------

#ifndef SENSITIVED_H
#define SENSITIVED_H

#include "SingleHit.hh"
#include "G4VSensitiveDetector.hh"

#include <unordered_map>

class G4Step;
class G4HCofThisEvent;
class G4TouchableHistory;

class SD_Det : public G4VSensitiveDetector
{
public:
  //
  SD_Det(const G4String& nameVol);
  virtual ~SD_Det();

  G4bool ProcessHits(G4Step* aStep, G4TouchableHistory* rohist);
  virtual void Initialize(G4HCofThisEvent* HCE);
  virtual void EndOfEvent(G4HCofThisEvent*);
  virtual void clear();

  SingleHitsCollection* fHitsCollection;
  G4int fHCID;

  std::unordered_map<int, std::unordered_map<int, int> > mapTrackID_Hits;
};

#endif // SENSITIVED_H
