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
// Implementation of the RunData class
// Created by C.Rappold (c.rappold@gsi.de)
//--------------------------------------------------------

#include "RunData.hh"

#include "G4RunManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4UnitsTable.hh"
#include "SingleHit.hh"
#include "G4LogicalVolumeStore.hh"

#include <map>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunData::RunData(const G4String& name)
    : G4Run(), namefile(name), fFactoryOn(false), WASAvol(0), eProc(), gProc(), volnames()
{
  std::cout << "RunData Ctr :" << name << " " << namefile << std::endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunData::~RunData() { ; }

void RunData::InitTree(const std::vector<G4String>& nameDet, Config& config)
{
	
  	writeTree = !!config.GetInt("WriteTree");
  	writeHistos = !!config.GetInt("WriteHistograms");
	G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
	if ( ! fFactoryOn ) {
		analysisManager->SetDefaultFileType("root");
		analysisManager->SetVerboseLevel(1);
		// Only merge in MT mode to avoid warning when running in Sequential mode
 #ifdef G4MULTITHREADED
 		analysisManager->SetNtupleMerging(true);
 #endif
  		analysisManager->SetFileName(namefile.c_str());
		G4bool fileOpen = analysisManager->OpenFile();
		if (! fileOpen) {
			G4cerr << "\n---> RunData::InitTree(): cannot open "
				<< analysisManager->GetFileName() << G4endl;
			return;
		}

		// create tree and branches
		if(writeTree){
			analysisManager->CreateNtuple("hibeam", "HIBEAM simulation output"); 
			fEvent=new PrimaryEvent;
			analysisManager->CreateNtupleIColumn("PrimaryTrackID", fEvent->TrackID);
			analysisManager->CreateNtupleIColumn("PrimaryPDG", fEvent->PDG);
			analysisManager->CreateNtupleDColumn("PrimaryEkin", fEvent->Ekin);
			analysisManager->CreateNtupleDColumn("PrimaryTime", fEvent->Time);
			analysisManager->CreateNtupleDColumn("PrimaryPosX", fEvent->PosX);
			analysisManager->CreateNtupleDColumn("PrimaryPosY", fEvent->PosY);
			analysisManager->CreateNtupleDColumn("PrimaryPosZ", fEvent->PosZ);
			analysisManager->CreateNtupleDColumn("PrimaryMomX", fEvent->MomX);
			analysisManager->CreateNtupleDColumn("PrimaryMomY", fEvent->MomY);
			analysisManager->CreateNtupleDColumn("PrimaryMomZ", fEvent->MomZ);
			analysisManager->CreateNtupleDColumn("PrimaryWeight", fEvent->Weight);
			for(size_t iDet=0; iDet<nameDet.size(); iDet++)
			{
				allHits.push_back(new MultiHits);
				analysisManager->CreateNtupleIColumn(nameDet[iDet]+"_TrackID", allHits[iDet]->TrackID);
				analysisManager->CreateNtupleIColumn(nameDet[iDet]+"_LayerID", allHits[iDet]->LayerID);
				analysisManager->CreateNtupleIColumn(nameDet[iDet]+"_LayerID1", allHits[iDet]->LayerID1);
				analysisManager->CreateNtupleIColumn(nameDet[iDet]+"_LayerID2", allHits[iDet]->LayerID2);
				analysisManager->CreateNtupleIColumn(nameDet[iDet]+"_PDG", allHits[iDet]->Pdg);
				analysisManager->CreateNtupleDColumn(nameDet[iDet]+"_EDep", allHits[iDet]->Edep);
				analysisManager->CreateNtupleDColumn(nameDet[iDet]+"_Time", allHits[iDet]->Time);
				analysisManager->CreateNtupleDColumn(nameDet[iDet]+"_TrackLength", allHits[iDet]->TrackLength);
				analysisManager->CreateNtupleDColumn(nameDet[iDet]+"_Position_X", allHits[iDet]->HitPosX);
				analysisManager->CreateNtupleDColumn(nameDet[iDet]+"_Position_Y", allHits[iDet]->HitPosY);
				analysisManager->CreateNtupleDColumn(nameDet[iDet]+"_Position_Z", allHits[iDet]->HitPosZ);
				analysisManager->CreateNtupleDColumn(nameDet[iDet]+"_GlobalPosition_X", allHits[iDet]->GlobalPosX);
				analysisManager->CreateNtupleDColumn(nameDet[iDet]+"_GlobalPosition_Y", allHits[iDet]->GlobalPosY);
				analysisManager->CreateNtupleDColumn(nameDet[iDet]+"_GlobalPosition_Z", allHits[iDet]->GlobalPosZ);
				analysisManager->CreateNtupleDColumn(nameDet[iDet]+"_Momentum_X", allHits[iDet]->MomX);
				analysisManager->CreateNtupleDColumn(nameDet[iDet]+"_Momentum_Y", allHits[iDet]->MomY);
				analysisManager->CreateNtupleDColumn(nameDet[iDet]+"_Momentum_Z", allHits[iDet]->MomZ);

			}
			analysisManager->FinishNtuple();  
		}
		fFactoryOn = true;
	}

	if(writeHistos){
		analysisManager->CreateH1("Flux_TPC_n","Neutron flux in TPC", 323, 1e-11*MeV, 1.0659E+03*MeV,"MeV","none","log");
		analysisManager->CreateH1("Flux_TPC_n_linear","Neutron flux in TPC", 500, 0*MeV, 25*MeV,"MeV","none");
		analysisManager->CreateH1("Flux_TPC_n_time","Neutron flux in TPC", 208, 10*ns, 1.0e10*ns,"ns","none","log");
		analysisManager->CreateH2("Flux_TPC_n_Etime","Energy-time correlation", 100, 1e-11*MeV, 1.0659E+03*MeV, 100, 1e1*ns, 1e10*ns, "MeV", "ns", "none", "none", "log", "log");

		analysisManager->CreateH1("Flux_TPC_g","Gamma flux in TPC", 323, 1e-11*MeV, 1.0659E+03*MeV,"MeV","none","log");
		analysisManager->CreateH1("Flux_TPC_g_linear","Gamma flux in TPC", 500, 0*MeV, 25*MeV,"MeV","none");
		analysisManager->CreateH1("Flux_TPC_g_time","Gamma flux in TPC", 208, 10*ns, 1.0e10*ns,"ns","none","log");
		analysisManager->CreateH2("Flux_TPC_g_Etime","Energy-time correlation", 100, 1e-11*MeV, 1.0659E+03*MeV, 100, 1e1*ns, 1e10*ns, "MeV", "ns", "none", "none", "log", "log");

		analysisManager->CreateH1("Flux_TPC_e","Electron flux in TPC", 323, 1e-11*MeV, 1.0659E+03*MeV,"MeV","none","log");
		analysisManager->CreateH1("Flux_TPC_e_linear","Electron flux in TPC", 500, 0*MeV, 25*MeV,"MeV","none");
		analysisManager->CreateH1("Flux_TPC_e_time","Electron flux in TPC", 208, 10*ns, 1.0e10*ns,"ns","none","log");
		analysisManager->CreateH2("Flux_TPC_e_Etime","Energy-time correlation", 100, 1e-11*MeV, 1.0659E+03*MeV, 100, 1e1*ns, 1e10*ns, "MeV", "ns", "none", "none", "log", "log");

		analysisManager->CreateH1("Flux_TPC_p","Positron flux in TPC", 323, 1e-11*MeV, 1.0659E+03*MeV,"MeV","none","log");
		analysisManager->CreateH1("Flux_TPC_p_linear","Positron flux in TPC", 500, 0*MeV, 25*MeV,"MeV","none");
		analysisManager->CreateH1("Flux_TPC_p_time","Positron flux in TPC", 208, 10*ns, 1.0e10*ns,"ns","none","log");
		analysisManager->CreateH2("Flux_TPC_p_Etime","Energy-time correlation", 100, 1e-11*MeV, 1.0659E+03*MeV, 100, 1e1*ns, 1e10*ns, "MeV", "ns", "none", "none", "log", "log");

		analysisManager->CreateH1("Flux_WASA_n","Neutron flux in WASA", 323, 1e-11*MeV, 1.0659E+03*MeV,"MeV","none","log");
		analysisManager->CreateH1("Flux_WASA_n_linear","Neutron flux in WASA", 500, 0*MeV, 25*MeV,"MeV","none");
		analysisManager->CreateH1("Flux_WASA_n_time","Neutron flux in WASA", 208, 10*ns, 1.0e10*ns,"ns","none","log");
		analysisManager->CreateH2("Flux_WASA_n_Etime","Energy-time correlation", 100, 1e-11*MeV, 1.0659E+03*MeV, 100, 1e1*ns, 1e10*ns, "MeV", "ns", "none", "none", "log", "log");

		analysisManager->CreateH1("Flux_WASA_g","Gamma flux in WASA", 323, 1e-11*MeV, 1.0659E+03*MeV,"MeV","none","log");
		analysisManager->CreateH1("Flux_WASA_g_linear","Gamma flux in WASA", 500, 0*MeV, 25*MeV,"MeV","none");
		analysisManager->CreateH1("Flux_WASA_g_time","Gamma flux in WASA", 208, 10*ns, 1.0e10*ns,"ns","none","log");
		analysisManager->CreateH2("Flux_WASA_g_Etime","Energy-time correlation", 100, 1e-11*MeV, 1.0659E+03*MeV, 100, 1e1*ns, 1e10*ns, "MeV", "ns", "none", "none", "log", "log");

		analysisManager->CreateH1("Flux_WASA_e","Electron flux in WASA", 323, 1e-11*MeV, 1.0659E+03*MeV,"MeV","none","log");
		analysisManager->CreateH1("Flux_WASA_e_linear","Electron flux in WASA", 500, 0*MeV, 25*MeV,"MeV","none");
		analysisManager->CreateH1("Flux_WASA_e_time","Electron flux in WASA", 208, 10*ns, 1.0e10*ns,"ns","none","log");
		analysisManager->CreateH2("Flux_WASA_e_Etime","Energy-time correlation", 100, 1e-11*MeV, 1.0659E+03*MeV, 100, 1e1*ns, 1e10*ns, "MeV", "ns", "none", "none", "log", "log");

		analysisManager->CreateH1("Flux_WASA_p","Positron flux in WASA", 323, 1e-11*MeV, 1.0659E+03*MeV,"MeV","none","log");
		analysisManager->CreateH1("Flux_WASA_p_linear","Positron flux in WASA", 500, 0*MeV, 25*MeV,"MeV","none");
		analysisManager->CreateH1("Flux_WASA_p_time","Positron flux in WASA", 208, 10*ns, 1.0e10*ns,"ns","none","log");
		analysisManager->CreateH2("Flux_WASA_p_Etime","Energy-time correlation", 100, 1e-11*MeV, 1.0659E+03*MeV, 100, 1e1*ns, 1e10*ns, "MeV", "ns", "none", "none", "log", "log");

		analysisManager->CreateH1("Flux_TARGET_n","Neutron flux in target foil", 323, 1e-11*MeV, 1.0659E+03*MeV,"MeV","none","log");
		analysisManager->CreateH1("Flux_TARGET_n_linear","Neutron flux in target foil", 500, 0*MeV, 25*MeV,"MeV","none");
		analysisManager->CreateH1("Flux_TARGET_n_time","Neutron flux in target foil", 208, 10*ns, 1.0e10*ns,"ns","none","log");
		analysisManager->CreateH2("Flux_TARGET_n_Etime","Energy-time correlation", 100, 1e-11*MeV, 1.0659E+03*MeV, 100, 1e1*ns, 1e10*ns, "MeV", "ns", "none", "none", "log", "log");

		analysisManager->CreateH1("Flux_TARGET_g","Gamma flux in target foil", 323, 1e-11*MeV, 1.0659E+03*MeV,"MeV","none","log");
		analysisManager->CreateH1("Flux_TARGET_g_linear","Gamma flux in target foil", 500, 0*MeV, 25*MeV,"MeV","none");
		analysisManager->CreateH1("Flux_TARGET_g_time","Gamma flux in target foil", 208, 10*ns, 1.0e10*ns,"ns","none","log");
		analysisManager->CreateH2("Flux_TARGET_g_Etime","Energy-time correlation", 100, 1e-11*MeV, 1.0659E+03*MeV, 100, 1e1*ns, 1e10*ns, "MeV", "ns", "none", "none", "log", "log");

		analysisManager->CreateH1("Flux_TARGET_e","Electron flux in target foil", 323, 1e-11*MeV, 1.0659E+03*MeV,"MeV","none","log");
		analysisManager->CreateH1("Flux_TARGET_e_linear","Electron flux in target foil", 500, 0*MeV, 25*MeV,"MeV","none");
		analysisManager->CreateH1("Flux_TARGET_e_time","Electron flux in target foil", 208, 10*ns, 1.0e10*ns,"ns","none","log");
		analysisManager->CreateH2("Flux_TARGET_e_Etime","Energy-time correlation", 100, 1e-11*MeV, 1.0659E+03*MeV, 100, 1e1*ns, 1e10*ns, "MeV", "ns", "none", "none", "log", "log");

		analysisManager->CreateH1("Flux_TARGET_p","Positron flux in target foil", 323, 1e-11*MeV, 1.0659E+03*MeV,"MeV","none","log");
		analysisManager->CreateH1("Flux_TARGET_p_linear","Positron flux in target foil", 500, 0*MeV, 25*MeV,"MeV","none");
		analysisManager->CreateH1("Flux_TARGET_p_time","Positron flux in target foil", 208, 10*ns, 1.0e10*ns,"ns","none","log");  
		analysisManager->CreateH2("Flux_TARGET_p_Etime","Energy-time correlation", 100, 1e-11*MeV, 1.0659E+03*MeV, 100, 1e1*ns, 1e10*ns, "MeV", "ns", "none", "none", "log", "log");

		analysisManager->CreateH2("Flux_TPC_origin","Origin of late high-energy neutrons", 200, -25*cm, 25*cm, 200, -25*cm, 25*cm, "cm", "cm", "none", "none");
		analysisManager->CreateH1("Flux_TPC_triton","Spectrum of late high-energy neutrons", 40, 0*MeV, 20*MeV,"MeV","none","linear");

		analysisManager->CreateH2("Process_TPC_g","Originating process", 6, -0.5, 5.5, 10, -0.5, 9.5, "none", "none");
		analysisManager->CreateH2("Process_TPC_e","Originating process", 5, -0.5, 4.5, 10, -0.5, 9.5, "none", "none");
		analysisManager->CreateH2("Process_WASA_g","Originating process", 6, -0.5, 5.5, 10, -0.5, 9.5, "none", "none");
		analysisManager->CreateH2("Process_WASA_e","Originating process", 5, -0.5, 4.5, 10, -0.5, 9.5, "none", "none");

		analysisManager->CreateH1("Current_TPC_e_enter","Entering surface current", 323, 1e-11*MeV, 1.0659E+03*MeV,"MeV","none","log");
		analysisManager->CreateH1("Current_TPC_e_leave","Exiting surface current", 323, 1e-11*MeV, 1.0659E+03*MeV,"MeV","none","log");
		analysisManager->CreateH1("Current_TPC_g_enter","Entering surface current", 323, 1e-11*MeV, 1.0659E+03*MeV,"MeV","none","log");
		analysisManager->CreateH1("Current_TPC_g_leave","Exiting surface current", 323, 1e-11*MeV, 1.0659E+03*MeV,"MeV","none","log");

		analysisManager->CreateH1("Current_WASA_e_enter","Entering surface current", 323, 1e-11*MeV, 1.0659E+03*MeV,"MeV","none","log");
		analysisManager->CreateH1("Current_WASA_e_leave","Exiting surface current", 323, 1e-11*MeV, 1.0659E+03*MeV,"MeV","none","log");
		analysisManager->CreateH1("Current_WASA_g_enter","Entering surface current", 323, 1e-11*MeV, 1.0659E+03*MeV,"MeV","none","log");
		analysisManager->CreateH1("Current_WASA_g_leave","Exiting surface current", 323, 1e-11*MeV, 1.0659E+03*MeV,"MeV","none","log");
	}
  gProc = { {"nCapture", 0},
            {"annihil", 1},
            {"eBrem", 2},
            {"neutronInelastic", 3},
            {"photonNuclear", 4} };

  eProc = { {"eIoni", 0},
            {"compt", 1},
            {"phot", 2},
            {"conv", 3} };

  volnames = { {"TPC", 0},
               {"FLANGE", 1},
               {"PIPE", 1},
               {"YOKE", 3},
               {"SHIELD", 4},
               {"VETO", 5},
               {"TARGET", 6},
               {"COATING", 7},
               {"STOP", 8},
               {"LiF", 9} };

  for (int i = 0; i <= 16; i++) {
    std::string key = "SECE" + std::to_string(i);
    volnames[key] = 2;
  }

  G4LogicalVolumeStore* lvStore = G4LogicalVolumeStore::GetInstance();
  
  for(size_t iDet=0; iDet<nameDet.size(); iDet++) {
    if (nameDet[iDet].rfind("SECE", 0) == 0) {
      WASAvol += 48*lvStore->GetVolume(nameDet[iDet])->GetSolid()->GetCubicVolume(); //48 crystals per ring
      // std::cout << "Summed WASA volume " << WASAvol/cm3 <<std::endl;
    }
  }

  if(writeHistos){
	  std::cout << "Target volume: " << lvStore->GetVolume("TARGET")->GetSolid()->GetCubicVolume()/cm3 << std::endl;
	  std::cout << "TPC volume: " << lvStore->GetVolume("TPC")->GetSolid()->GetCubicVolume()/cm3 << std::endl;
	  std::cout << "WASA volume: " << WASAvol/cm3 << std::endl;
  }

  G4cout << "\n----> Output file is open in "
		<< analysisManager->GetFileName() << "."
		<< analysisManager->GetFileType() << G4endl;

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunData::RecordEvent(const G4Event* event)
{
	if (! fFactoryOn) return;
	if (! writeTree) return;
	G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();

	unsigned int Nprimary = event->GetNumberOfPrimaryVertex();
	if(Nprimary > 0)
	{
		// G4cout<<" Write Event str :"<<Nprimary<<G4endl;
		for(unsigned int iPrimary = 0; iPrimary < Nprimary; ++iPrimary)
		{
			G4PrimaryVertex* PrimVertex = event->GetPrimaryVertex(iPrimary);

			if(PrimVertex != nullptr)
			{
				unsigned int Nparticle     = PrimVertex->GetNumberOfParticle();
				//G4cout<< "Vtx#"<<iPrimary<<" Nparticle :"<<Nparticle<<G4endl;

				for(unsigned int iParticle = 0; iParticle < Nparticle; ++iParticle)
				{
					G4PrimaryParticle* PrimParticle = PrimVertex->GetPrimary(iParticle);

					if(PrimParticle != nullptr)
					{
						fEvent->PDG.push_back(PrimParticle->GetPDGcode());
						fEvent->Ekin.push_back(PrimParticle->GetKineticEnergy());
						fEvent->Time.push_back(PrimVertex->GetT0() / ms);
						fEvent->PosX.push_back(PrimVertex->GetX0() / cm);
						fEvent->PosY.push_back(PrimVertex->GetY0() / cm);
						fEvent->PosZ.push_back(PrimVertex->GetZ0() / cm);
						fEvent->MomX.push_back(PrimParticle->GetPx() / MeV);
						fEvent->MomY.push_back(PrimParticle->GetPy() / MeV);
						fEvent->MomZ.push_back(PrimParticle->GetPz() / MeV);
						fEvent->TrackID.push_back(PrimParticle->GetTrackID());
						fEvent->Weight.push_back(PrimParticle->GetWeight());
					}
				}
			}
		}
	}

	G4HCofThisEvent* hce = event->GetHCofThisEvent();
	if(!hce)
	{
		G4ExceptionDescription msg;
		msg << "No hits collection of this event found." << G4endl;
		G4Exception("RunData::RecordEvent()", "Code001", JustWarning, msg);
		return;
	}

	for(unsigned int idCol = 0; idCol < allHits.size(); ++idCol)
	{
		SingleHitsCollection* TempCol = dynamic_cast<SingleHitsCollection*>(hce->GetHC(idCol));
		if(TempCol != nullptr)
		{
			G4int nhits = TempCol->entries();
			// std::cout<<" Col "<<TempCol->GetName()<<" "<<nhits<<std::endl;
			for(G4int ihit = 0; ihit < nhits; ++ihit)
			{
				SingleHit* TempHit   = (*TempCol)[ihit];
				if(TempHit->Edep>0){
					allHits[idCol]->TrackID.push_back( TempHit->TrackID);
					allHits[idCol]->LayerID.push_back( TempHit->LayerID);
					allHits[idCol]->LayerID1.push_back( TempHit->LayerID1);
					allHits[idCol]->LayerID2.push_back( TempHit->LayerID2);
					allHits[idCol]->HitPosX.push_back( TempHit->HitPosX / cm);
					allHits[idCol]->HitPosY.push_back( TempHit->HitPosY / cm);
					allHits[idCol]->HitPosZ.push_back( TempHit->HitPosZ / cm);
					allHits[idCol]->GlobalPosX.push_back( TempHit->GlobalPosX / cm);
					allHits[idCol]->GlobalPosY.push_back( TempHit->GlobalPosY / cm);
					allHits[idCol]->GlobalPosZ.push_back( TempHit->GlobalPosZ / cm);
					allHits[idCol]->MomX.push_back( TempHit->MomX / GeV);
					allHits[idCol]->MomY.push_back( TempHit->MomY / GeV);
					allHits[idCol]->MomZ.push_back( TempHit->MomZ / GeV);
					allHits[idCol]->Edep.push_back( TempHit->Edep / MeV);
					allHits[idCol]->Time.push_back( TempHit->Time / ns);
					allHits[idCol]->TrackLength.push_back( TempHit->TrackLength / cm);
					allHits[idCol]->Pdg.push_back( TempHit->Pdg);
				}
			}
		}
		else
		{
			G4ExceptionDescription msg;
			msg << "Some of hits collections of this event not found." << G4endl;
			G4Exception("RunData::RecordEvent()", "Code002", JustWarning, msg);
			return;
		}
	}

	analysisManager->AddNtupleRow(0);	
}

void RunData::AccumulateFlux(const G4Step* step) 
{
	if (! writeHistos) return;
    // get volume of the current step
  G4LogicalVolume* volume
    = step->GetPreStepPoint()->GetTouchableHandle()
      ->GetVolume()->GetLogicalVolume();

  // collect relevant information about time step
  G4double stepLen = step->GetStepLength();
  G4double energy = step->GetPreStepPoint()->GetKineticEnergy();
  G4double initenergy = step->GetTrack()->GetVertexKineticEnergy();
  G4double time = (step->GetPreStepPoint()->GetGlobalTime()+step->GetPostStepPoint()->GetGlobalTime())/2;
  G4double weight = step->GetTrack()->GetWeight();
  G4String type = step->GetTrack()->GetParticleDefinition()->GetParticleName();
  G4ThreeVector vertex = step->GetTrack()->GetVertexPosition();

  auto analysisManager = G4AnalysisManager::Instance();  
  if (volume->GetName() == "TPC") {
    if (type == "neutron") {
      analysisManager->FillH1(0, energy, stepLen/cm*weight/(volume->GetSolid()->GetCubicVolume()/cm3));
      analysisManager->FillH1(1, energy, stepLen/cm*weight/(volume->GetSolid()->GetCubicVolume()/cm3));
      analysisManager->FillH1(2, time/ns, stepLen/cm*weight/(volume->GetSolid()->GetCubicVolume()/cm3)/(time/ns));
      analysisManager->FillH2(0, energy, time/ns, stepLen/cm*weight/(volume->GetSolid()->GetCubicVolume()/cm3));
    }
    if (type == "gamma") {
      analysisManager->FillH1(3, energy, stepLen/cm*weight/(volume->GetSolid()->GetCubicVolume()/cm3));
      analysisManager->FillH1(4, energy, stepLen/cm*weight/(volume->GetSolid()->GetCubicVolume()/cm3));
      analysisManager->FillH1(5, time/ns, stepLen/cm*weight/(volume->GetSolid()->GetCubicVolume()/cm3)/(time/ns));
      analysisManager->FillH2(1, energy, time/ns, stepLen/cm*weight/(volume->GetSolid()->GetCubicVolume()/cm3));
    }
    if (type == "e-") {
      analysisManager->FillH1(6, energy, stepLen/cm*weight/(volume->GetSolid()->GetCubicVolume()/cm3));
      analysisManager->FillH1(7, energy, stepLen/cm*weight/(volume->GetSolid()->GetCubicVolume()/cm3));
      analysisManager->FillH1(8, time/ns, stepLen/cm*weight/(volume->GetSolid()->GetCubicVolume()/cm3)/(time/ns));
      analysisManager->FillH2(2, energy, time/ns, stepLen/cm*weight/(volume->GetSolid()->GetCubicVolume()/cm3));
    }
    if (type == "e+") {
      analysisManager->FillH1(9, energy, stepLen/cm*weight/(volume->GetSolid()->GetCubicVolume()/cm3));
      analysisManager->FillH1(10, energy, stepLen/cm*weight/(volume->GetSolid()->GetCubicVolume()/cm3));
      analysisManager->FillH1(11, time/ns, stepLen/cm*weight/(volume->GetSolid()->GetCubicVolume()/cm3)/(time/ns));
      analysisManager->FillH2(3, energy, time/ns, stepLen/cm*weight/(volume->GetSolid()->GetCubicVolume()/cm3));
    }
  }
  else if (volume->GetName().rfind("SECE", 0) == 0) {
    if (type == "neutron") {
      analysisManager->FillH1(12, energy, stepLen/cm*weight/(WASAvol/cm3));
      analysisManager->FillH1(13, energy, stepLen/cm*weight/(WASAvol/cm3));
      analysisManager->FillH1(14, time/ns, stepLen/cm*weight/(WASAvol/cm3)/(time/ns));
      analysisManager->FillH2(4, energy, time/ns, stepLen/cm*weight/(WASAvol/cm3));
    }
    if (type == "gamma") {
      analysisManager->FillH1(15, energy, stepLen/cm*weight/(WASAvol/cm3));
      analysisManager->FillH1(16, energy, stepLen/cm*weight/(WASAvol/cm3));
      analysisManager->FillH1(17, time/ns, stepLen/cm*weight/(WASAvol/cm3)/(time/ns));
      analysisManager->FillH2(5, energy, time/ns, stepLen/cm*weight/(WASAvol/cm3));
    }
    if (type == "e-") {
      analysisManager->FillH1(18, energy, stepLen/cm*weight/(WASAvol/cm3));
      analysisManager->FillH1(19, energy, stepLen/cm*weight/(WASAvol/cm3));
      analysisManager->FillH1(20, time/ns, stepLen/cm*weight/(WASAvol/cm3)/(time/ns));
      analysisManager->FillH2(6, energy, time/ns, stepLen/cm*weight/(WASAvol/cm3));
    }
    if (type == "e+") {
      analysisManager->FillH1(21, energy, stepLen/cm*weight/(WASAvol/cm3));
      analysisManager->FillH1(22, energy, stepLen/cm*weight/(WASAvol/cm3));
      analysisManager->FillH1(23, time/ns, stepLen/cm*weight/(WASAvol/cm3)/(time/ns));
      analysisManager->FillH2(7, energy, time/ns, stepLen/cm*weight/(WASAvol/cm3));
    }
  } else if (volume->GetName() == "TARGET") {
    if (type == "neutron") {
      analysisManager->FillH1(24, energy, stepLen/cm*weight/(volume->GetSolid()->GetCubicVolume()/cm3));
      analysisManager->FillH1(25, energy, stepLen/cm*weight/(volume->GetSolid()->GetCubicVolume()/cm3));
      analysisManager->FillH1(26, time/ns, stepLen/cm*weight/(volume->GetSolid()->GetCubicVolume()/cm3)/(time/ns));
      analysisManager->FillH2(8, energy, time/ns, stepLen/cm*weight/(volume->GetSolid()->GetCubicVolume()/cm3));
    }
    if (type == "gamma") {
      analysisManager->FillH1(27, energy, stepLen/cm*weight/(volume->GetSolid()->GetCubicVolume()/cm3));
      analysisManager->FillH1(28, energy, stepLen/cm*weight/(volume->GetSolid()->GetCubicVolume()/cm3));
      analysisManager->FillH1(29, time/ns, stepLen/cm*weight/(volume->GetSolid()->GetCubicVolume()/cm3)/(time/ns));
      analysisManager->FillH2(9, energy, time/ns, stepLen/cm*weight/(volume->GetSolid()->GetCubicVolume()/cm3)); 
    } 
    if (type == "e-") {
      analysisManager->FillH1(30, energy, stepLen/cm*weight/(volume->GetSolid()->GetCubicVolume()/cm3));
      analysisManager->FillH1(31, energy, stepLen/cm*weight/(volume->GetSolid()->GetCubicVolume()/cm3));
      analysisManager->FillH1(32, time/ns, stepLen/cm*weight/(volume->GetSolid()->GetCubicVolume()/cm3)/(time/ns));
      analysisManager->FillH2(10, energy, time/ns, stepLen/cm*weight/(volume->GetSolid()->GetCubicVolume()/cm3));
    }
    if (type == "e+") {
      analysisManager->FillH1(33, energy, stepLen/cm*weight/(volume->GetSolid()->GetCubicVolume()/cm3));
      analysisManager->FillH1(34, energy, stepLen/cm*weight/(volume->GetSolid()->GetCubicVolume()/cm3));
      analysisManager->FillH1(35, time/ns, stepLen/cm*weight/(volume->GetSolid()->GetCubicVolume()/cm3)/(time/ns));
      analysisManager->FillH2(11, energy, time/ns, stepLen/cm*weight/(volume->GetSolid()->GetCubicVolume()/cm3));
    }
  }

  if (volume->GetName() == "TPC" && type == "neutron" && step->GetTrack()->GetCreatorProcess() && step->GetTrack()->GetCreatorProcess()->GetProcessName() == "tInelastic") {
    analysisManager->FillH2(12, vertex[0], vertex[1]);
    analysisManager->FillH1(36, energy);
  }

  if (volume->GetName() == "TPC" && type == "gamma" && step->GetTrack()->GetCreatorProcess()) {
    int volidx = volnames[step->GetTrack()->GetLogicalVolumeAtVertex()->GetName()];
    int procidx;
    if (gProc.find(step->GetTrack()->GetCreatorProcess()->GetProcessName()) == gProc.end()) {
      // std::cout << step->GetTrack()->GetCreatorProcess()->GetProcessName() << std::endl;
      procidx = 5;
    } else {
      procidx = gProc[step->GetTrack()->GetCreatorProcess()->GetProcessName()];
    }
    analysisManager->FillH2(13, procidx, volidx, stepLen/cm*weight/(volume->GetSolid()->GetCubicVolume()/cm3));
  }

  if (volume->GetName() == "TPC" && type == "e-" && step->GetTrack()->GetCreatorProcess()) {
    int volidx = volnames[step->GetTrack()->GetLogicalVolumeAtVertex()->GetName()];
    int procidx;
    if (eProc.find(step->GetTrack()->GetCreatorProcess()->GetProcessName()) == eProc.end()) {
      procidx = 4;
    } else {
      procidx = eProc[step->GetTrack()->GetCreatorProcess()->GetProcessName()];
    }
    analysisManager->FillH2(14, procidx, volidx, stepLen/cm*weight/(volume->GetSolid()->GetCubicVolume()/cm3));
  }

    if (volume->GetName().rfind("SECE", 0) == 0 && type == "gamma" && step->GetTrack()->GetCreatorProcess()) {
    int volidx = volnames[step->GetTrack()->GetLogicalVolumeAtVertex()->GetName()];
    int procidx;
    if (gProc.find(step->GetTrack()->GetCreatorProcess()->GetProcessName()) == gProc.end()) {
      procidx = 5;
    } else {
      procidx = gProc[step->GetTrack()->GetCreatorProcess()->GetProcessName()];
    }
    analysisManager->FillH2(15, procidx, volidx, stepLen/cm*weight/(WASAvol/cm3));
  }

  if (volume->GetName().rfind("SECE", 0) == 0 && type == "e-" && step->GetTrack()->GetCreatorProcess()) {
    int volidx = volnames[step->GetTrack()->GetLogicalVolumeAtVertex()->GetName()];
    int procidx;
    if (eProc.find(step->GetTrack()->GetCreatorProcess()->GetProcessName()) == eProc.end()) {
      procidx = 4;
    } else {
      procidx = eProc[step->GetTrack()->GetCreatorProcess()->GetProcessName()];
    }
    analysisManager->FillH2(16, procidx, volidx, stepLen/cm*weight/(WASAvol/cm3));
  }

  if (volume->GetName() == "TPC" && type == "e-" && step->GetPreStepPoint()->GetStepStatus() == fGeomBoundary) {
    analysisManager->FillH1(37, energy, weight);
  }
  if (volume->GetName() == "TPC" && type == "e-" && step->GetPostStepPoint()->GetStepStatus() == fGeomBoundary) {
    analysisManager->FillH1(38, energy, weight);
  }
  if (volume->GetName() == "TPC" && type == "gamma" && step->GetPreStepPoint()->GetStepStatus() == fGeomBoundary) {
    analysisManager->FillH1(39, energy, weight);
  } 
  if (volume->GetName() == "TPC" && type == "gamma" && step->GetPostStepPoint()->GetStepStatus() == fGeomBoundary) {
    analysisManager->FillH1(40, energy, weight);
  }

  if (volume->GetName().rfind("SECE", 0) == 0 && type == "e-" && step->GetPreStepPoint()->GetStepStatus() == fGeomBoundary) {
    analysisManager->FillH1(41, energy, weight);
  }
  if (volume->GetName().rfind("SECE", 0) == 0 && type == "e-" && step->GetPostStepPoint()->GetStepStatus() == fGeomBoundary) {
    analysisManager->FillH1(42, energy, weight);
  }
  if (volume->GetName().rfind("SECE", 0) == 0 && type == "gamma" && step->GetPreStepPoint()->GetStepStatus() == fGeomBoundary) {
    analysisManager->FillH1(43, energy, weight);
  } 
  if (volume->GetName().rfind("SECE", 0) == 0 && type == "gamma" && step->GetPostStepPoint()->GetStepStatus() == fGeomBoundary) {
    analysisManager->FillH1(44, energy, weight);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunData::Reset()
{
	if (! writeTree) return;
  for(size_t i = 0; i < allHits.size(); ++i) allHits[i]->Reset();
  fEvent->Zero();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void RunData::Merge(const G4Run* aRun)
{
  const RunData* localRun = static_cast<const RunData*>(aRun);
  G4Run::Merge(aRun);
} 
