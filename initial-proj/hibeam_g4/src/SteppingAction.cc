// -------------------------------------------------------
// Implementation of the SteppingAction class
// Created by L. Persson (selipe150@gmail.com)
//--------------------------------------------------------

#include "SteppingAction.hh"
#include "RunAction.hh"

#include "G4Step.hh"
#include "G4Event.hh"
#include "G4RunManager.hh"
#include "G4LogicalVolume.hh"
#include "G4ParticleDefinition.hh"


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SteppingAction::SteppingAction(RunAction* runAction, const G4String& OutputFile)
: fRunAction(runAction), fOutputFile(OutputFile)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SteppingAction::UserSteppingAction(const G4Step* step)
{
  fRunAction->AddFlux(step);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......