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
// Implementation of the EventAction class
// Created by C.Rappold (c.rappold@gsi.de)
//--------------------------------------------------------

#include "EventAction.hh"

#include "RunAction.hh"
#include "G4Event.hh"
#include "G4RunManager.hh"
#include "G4SDManager.hh"
#include "RunData.hh"

//#include "Randomize.hh"
#include <iomanip>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

EventAction::EventAction(const std::vector<G4String>& name) : G4UserEventAction(), nameDetector(name)
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

EventAction::~EventAction() {}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::BeginOfEventAction(const G4Event* /*event*/)
{
  // for(auto&  array : list_Arrays)
  //   array.Clear("C");

  if(HCID.size() != nameDetector.size())
    {
      G4SDManager* sdManager = G4SDManager::GetSDMpointer();
      // std::cout<<"!> SDmanager :"<<std::endl;
      // sdManager->ListTree();
      for(auto name : nameDetector)
        {
          G4String tempName(name);
          tempName += "/Coll";
          HCID.push_back(sdManager->GetCollectionID(tempName));
        }
    }

  RunData* runData = dynamic_cast<RunData*>(G4RunManager::GetRunManager()->GetNonConstCurrentRun());
  runData->Reset();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::EndOfEventAction(const G4Event* event)
{
  RunData* runData = dynamic_cast<RunData*>(G4RunManager::GetRunManager()->GetNonConstCurrentRun());
  //runData->FillPerEvent(event);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
