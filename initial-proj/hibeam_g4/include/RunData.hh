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
// -------------------------------------------------
// Definition of the RunData class
// Created by C.Rappold (c.rappold@gsi.de)
//--------------------------------------------------

#ifndef RUNDATA_H
#define RUNDATA_H 1

#include "G4Event.hh"
#include "G4Run.hh"
#include "G4AnalysisManager.hh"
#include "Config.hh"
#include "PrimaryEvent.hh"
#include "MultiHits.hh"
#include "globals.hh"
#include "G4LogicalVolume.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

///  Run data class

class RunData : public G4Run
{
public:
  explicit RunData(const G4String& namefile);
  virtual ~RunData();

  void RecordEvent(const G4Event* event);
  void Merge(const G4Run*);

  void Reset();
  void InitTree(const std::vector<G4String>& nameDet, Config& config);
  void AccumulateFlux(const G4Step*);

private:
  const G4String& namefile;

  PrimaryEvent* fEvent;
  std::vector<MultiHits*> allHits;
  G4double WASAvol;
  std::map<std::string,int> eProc;
  std::map<std::string,int> gProc;
  std::map<std::string,int> volnames;

  bool fFactoryOn;
  bool writeTree;
  bool writeHistos;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif // RUNDATA_H
