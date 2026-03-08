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
// Implementation of the ActionInitialization class
// Created by C.Rappold (c.rappold@gsi.de)
//--------------------------------------------------------

#include "ActionInitialization.hh"

#include "EventAction.hh"
#include "RunAction.hh"
#include "SteppingAction.hh"
#include "PrimaryGeneratorAction.hh"
#include "ScatteringGeneratorAction.hh"
#include "G4MCPLGenerator.hh"

#include <memory>
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

ActionInitialization::ActionInitialization(GeometryController* geoControl, Config& ConfFile)
    : G4VUserActionInitialization(), fGeoController(geoControl), Conf(ConfFile)
{
  OutputFile = Conf.GetString("Output_Namefile");
  std::cout << "ActionInit : done ! " << OutputFile << std::endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

ActionInitialization::~ActionInitialization() {}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ActionInitialization::BuildForMaster() const
{
  std::vector<G4String> nameD = fGeoController->GetNameDetectors();
  std::cout << "!> ActionInitialization BuildForMaster:" << nameD.size() << " "
            << fGeoController->GetNameDetectors().size() << std::endl;

  SetUserAction(new RunAction(OutputFile, fGeoController->GetNameDetectors(), Conf));
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ActionInitialization::Build() const
{
  std::vector<G4String> nameD = fGeoController->GetNameDetectors();
  std::cout << "!> ActionInitialization Build:" << nameD.size() << " " << fGeoController->GetNameDetectors().size()
            << std::endl;

  std::cout << "================================" << std::endl;
  const std::string source = Conf.GetString("Source");
  if(source=="mcpl"||source=="MCPL"){
  	std::cout << "Using MCPL source!" << std::endl;
  	SetUserAction(new G4MCPLGenerator(Conf));
  }
  else if(source=="scattering"||source=="Scattering"){
  	std::cout << "Using scattering source!" << std::endl;
	SetUserAction(new ScatteringGeneratorAction);
  }
  else if(source=="gps"||source=="GPS"){
  	std::cout << "Using general particle source!" << std::endl;
	SetUserAction(new PrimaryGeneratorAction);
  }
  std::cout << "================================" << std::endl;

  auto runAction = new RunAction(OutputFile, fGeoController->GetNameDetectors(), Conf);
  SetUserAction(runAction);

  EventAction* eventAction = new EventAction(fGeoController->GetNameDetectors());
  SetUserAction(eventAction);

  SetUserAction(new SteppingAction(runAction, OutputFile));
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
