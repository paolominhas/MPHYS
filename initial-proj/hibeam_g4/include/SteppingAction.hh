#ifndef SteppingAction_h
#define SteppingAction_h 1

#include "G4UserSteppingAction.hh"
#include "globals.hh"

class G4LogicalVolume;
class RunAction;
class Config;

/// Stepping action class

class SteppingAction : public G4UserSteppingAction
{
  public:
    SteppingAction(RunAction* runAction, const G4String& OutputFile);
    ~SteppingAction() override = default;

    // method from the base class
    void UserSteppingAction(const G4Step*) override;

  private:
    RunAction* fRunAction = nullptr;
    const G4String& fOutputFile;
    G4LogicalVolume* fScoringVolume = nullptr;
};


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
