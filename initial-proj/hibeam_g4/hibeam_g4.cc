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
//----------------------------------------------
// hibeam+g4 main of the G4 simulation
// Created by M.Holl (matthiasholl@ess.eu)
// Based on SolenoidSimple
// Created by C.Rappold (c.rappold@gsi.de)
//----------------------------------------------


#include <G4ProductionCuts.hh>
#include "G4SystemOfUnits.hh"

#include "VDetectorConstruction.hh"
#include "ActionInitialization.hh"
#include "Config.hh"
#include "ConvertGeo.hh"
#include "PhysicsList.hh"

#include "G4RunManagerFactory.hh"
#include "G4ScoringManager.hh"

#include "G4UImanager.hh"
#include "G4VisExecutive.hh"
#include "G4UIExecutive.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

int main(int argc,char** argv)
{

  std::cout<<" Config :"<<std::endl;
  Config config(argc,argv);
  if(config.ProperConf()!=0)
    return -1;
  config.CheckConfig();

  int guimode = config.GetInt("Gui"); 

  if(config.IsAvailable("HEPRand_Seed"))
    {
      std::cout<<"!> Set Seed for HepRandom :\n";
      //CLHEP::HepRandom::setTheEngine(new CLHEP::RanecuEngine);
      //CLHEP::HepRandom::setTheEngine(new CLHEP::MixMaxRng);
      CLHEP::HepRandom::setTheSeed(config.GetInt("HEPRand_Seed"));
      CLHEP::HepRandom::showEngineStatus();
    }
  else
    CLHEP::HepRandom::showEngineStatus();
  
  // Enable ui session if guimode is turned on
  G4UIExecutive* ui = nullptr;
  if ( guimode ) { ui = new G4UIExecutive(argc, argv); }
  
  // Construct the default run manager
  auto runManager = G4RunManagerFactory::CreateRunManager(G4RunManagerType::Default);
   if(config.IsAvailable("Threads"))
    {
      runManager->SetNumberOfThreads(config.GetInt("Threads"));
    }
  // Activate command-based scorer
  G4ScoringManager::GetScoringManager();
  
  std::string nameGeo = config.GetString("Geo");
  auto geometryController = new GeometryController(config);
  geometryController->SetGeometry(nameGeo);

  auto physicsList = new PhysicsList();
  G4ProductionCutsTable::GetProductionCutsTable()->SetEnergyRange(30.0*eV, 10.0*TeV);
  runManager->SetUserInitialization(physicsList);
	
  // User action initialization
  runManager->SetUserInitialization(new ActionInitialization(geometryController,config));
 
  // Visualization manager construction
  auto visManager = new G4VisExecutive;
  visManager->Initialize();
    
  // Get the pointer to the User Interface manager
  auto UImanager = G4UImanager::GetUIpointer();

  bool macro = config.IsAvailable("MacroFile");
  if ( !ui && macro)
    {
      // execute an argument macro file if exist
	  std::string command = "/control/execute ";
	  std::string fileName = config.GetString("MacroFile"); 
      UImanager->ApplyCommand(command+fileName);
    }
  else
    {
      UImanager->ApplyCommand("/control/execute init_vis.mac"); 
	    UImanager->ApplyCommand("/control/execute gui.mac");     
      // start interactive session
      ui->SessionStart();
      delete ui;
    }

  std::cout<<" Doing G4->ROOT convertion ? ";
  
  bool convert = config.IsAvailable("ConvertRoot");
  if(convert)
    {
      std::cout<<" yes"<<std::endl;
	  std::string NameConvert = config.GetString("ConvertRoot");
      geometryController->ConvertG4toRoot(NameConvert);
    }
  else
    std::cout<<" no"<<std::endl;

  // Job termination
  // Free the store: user actions, physics_list and detector_description are
  // owned and deleted by the run manager, so they should not be deleted 
  // in the main() program !

  delete visManager;
  delete runManager;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
