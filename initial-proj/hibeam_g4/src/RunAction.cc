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
// -------------------------------------------------------
// Implementation of the RunAction class
// Created by C.Rappold (c.rappold@gsi.de)
//--------------------------------------------------------

#include "RunAction.hh"

#include "RunData.hh"

#include "G4Run.hh"
#include "G4RunManager.hh"
#include "G4SDManager.hh"
#include "G4AnalysisManager.hh"
#include "G4TScoreHistFiller.hh"
#include "G4SystemOfUnits.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunAction::RunAction(const G4String& name, const std::vector<G4String>& nameSD_Det, Config& config)
    : G4UserRunAction(), OutputFileName(name), NameDetectorsSD(nameSD_Det), Conf(config), run(nullptr)
{
  // set printing event number per each event
  // G4RunManager::GetRunManager()->SetPrintProgress(0);
  
  std::cout << "!> RunAction Ctr :" << OutputFileName << " " << NameDetectorsSD.size() << " " << nameSD_Det.size()
            << std::endl;
  // for(auto& nameD : NameDetectorsSD)
  //   std::cout<<"-> "<<nameD<<std::endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
RunAction::~RunAction()
{
}

G4Run* RunAction::GenerateRun() { return (new RunData(OutputFileName)); }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::BeginOfRunAction(const G4Run*)
{
  run = dynamic_cast<RunData*>(G4RunManager::GetRunManager()->GetNonConstCurrentRun());
  run->InitTree(NameDetectorsSD, Conf);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::AddFlux(const G4Step* step)
{
  run->AccumulateFlux(step);
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::EndOfRunAction(const G4Run*)
{
  G4cout << " End Run : Closing root file ";
  
  auto analysisManager = G4AnalysisManager::Instance();  
	analysisManager->Write();
	analysisManager->CloseFile();
  
  G4cout << " done !" << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
