// -------------------------------------------------------
// Implementation of the Sample_Det class
// Created by C.Rappold (c.rappold@gsi.de)
//--------------------------------------------------------

#include "SamplingD.hh"

#include "G4Material.hh"
#include "G4ParticleTypes.hh"
#include "G4SDManager.hh"
#include "G4Step.hh"
#include "G4VSolid.hh"
#include "G4VTouchable.hh"

#include <algorithm>
#include <sstream>
#include <string>
#include <vector>

Sample_Det::Sample_Det(const G4String& Vol_name)
    : G4VSensitiveDetector(Vol_name), fHitsCollection(nullptr), fHCID(-1)
{
  collectionName.insert("Coll");
}

Sample_Det::~Sample_Det() {}

void Sample_Det::Initialize(G4HCofThisEvent* hce)
{
  fHitsCollection = new SingleHitsCollection(SensitiveDetectorName, collectionName[0]);
  if(fHCID < 0)
    fHCID = G4SDManager::GetSDMpointer()->GetCollectionID(fHitsCollection);

  hce->AddHitsCollection(fHCID, fHitsCollection);
  if(mapTrackID_Hits.size() != 0)
    mapTrackID_Hits.clear();
}

void Sample_Det::clear() { mapTrackID_Hits.clear(); }

void Sample_Det::EndOfEvent(G4HCofThisEvent*)
{
  //   Hits->Clear("C");
  mapTrackID_Hits.clear();
}


//.....
G4bool Sample_Det::ProcessHits(G4Step* aStep, G4TouchableHistory* )
{

  	G4double energ_depos = aStep->GetTotalEnergyDeposit();
	
	G4String PhysName(aStep->GetPreStepPoint()->GetPhysicalVolume()->GetName());

	const G4TouchableHistory* touchable = dynamic_cast<const G4TouchableHistory*>(aStep->GetPreStepPoint()->GetTouchable());
	G4VPhysicalVolume* currentPhysical = touchable->GetVolume(0);
	G4int copyNo = currentPhysical->GetCopyNo();
	G4int copyNo1 = -1;
	if(touchable->GetHistoryDepth()>1) copyNo1 = touchable->GetCopyNumber(1);
	G4int copyNo2 = -1;
	if(touchable->GetHistoryDepth()>2) copyNo2 = touchable->GetCopyNumber(2);

	G4StepPoint* PreStep = aStep->GetPreStepPoint();
	G4StepPoint* PostStep = aStep->GetPostStepPoint();
	int CurrentTrack = aStep->GetTrack()->GetTrackID();

	auto it_layer = mapTrackID_Hits.find(copyNo);
	if(PostStep->GetStepStatus()==fGeomBoundary)
	{
		int IdHit        = fHitsCollection->GetSize();
		SingleHit* newHit = new SingleHit(fHCID);

		newHit->Pdg         = aStep->GetTrack()->GetDefinition()->GetPDGEncoding();
		newHit->TrackID     = aStep->GetTrack()->GetTrackID();
		newHit->Edep      = energ_depos;
		newHit->Time        = aStep->GetPreStepPoint()->GetGlobalTime();
		newHit->TrackLength = aStep->GetTrack()->GetTrackLength();
		newHit->HitPosX     = aStep->GetPreStepPoint()->GetPosition().x();
		newHit->HitPosY     = aStep->GetPreStepPoint()->GetPosition().y();
		newHit->HitPosZ     = aStep->GetPreStepPoint()->GetPosition().z();
		newHit->MomX        = aStep->GetTrack()->GetMomentum().x();
		newHit->MomY        = aStep->GetTrack()->GetMomentum().y();
		newHit->MomZ        = aStep->GetTrack()->GetMomentum().z();
		newHit->LayerID     = copyNo;
		newHit->LayerID1    = copyNo1;
		newHit->LayerID2    = copyNo2;

		mapTrackID_Hits.insert(std::pair<int, std::unordered_map<int, int> >(copyNo, {{CurrentTrack, IdHit}}));

		fHitsCollection->insert(newHit);
	}

	return true;
}

