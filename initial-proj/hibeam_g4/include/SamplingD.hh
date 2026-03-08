//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
/// -------------------------------------------------
// Definition of the Sample_Det class
// Created by C.Rappold (c.rappold@gsi.de)
//--------------------------------------------------

#ifndef SAMPLINGD_H
#define SAMPLINGD_H 1

#include "SingleHit.hh"
#include "G4VSensitiveDetector.hh"

#include <unordered_map>

class G4Step;
class G4HCofThisEvent;
class G4TouchableHistory;
class Sample_Det : public G4VSensitiveDetector
{
public:
    Sample_Det(const G4String& nameVol);
	~Sample_Det();

	G4bool ProcessHits(G4Step*aStep,G4TouchableHistory*ROhist);
	virtual void Initialize(G4HCofThisEvent*);
	virtual void EndOfEvent(G4HCofThisEvent*HCE);
	virtual void clear();

	SingleHitsCollection* fHitsCollection;
	G4int fHCID;

	std::unordered_map<int, std::unordered_map<int, int> > mapTrackID_Hits;
};
#endif


