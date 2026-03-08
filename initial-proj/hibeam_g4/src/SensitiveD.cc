// -------------------------------------------------------
// Implementation of the SD_Det class
// Created by C.Rappold (c.rappold@gsi.de)
//--------------------------------------------------------

#include "SensitiveD.hh"

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

SD_Det::SD_Det(const G4String& Vol_name)
    : G4VSensitiveDetector(Vol_name), fHitsCollection(nullptr), fHCID(-1)
{
  collectionName.insert("Coll");
}

SD_Det::~SD_Det() {}

void SD_Det::Initialize(G4HCofThisEvent* hce)
{
  fHitsCollection = new SingleHitsCollection(SensitiveDetectorName, collectionName[0]);
  if(fHCID < 0)
    fHCID = G4SDManager::GetSDMpointer()->GetCollectionID(fHitsCollection);

  hce->AddHitsCollection(fHCID, fHitsCollection);
  if(mapTrackID_Hits.size() != 0)
    mapTrackID_Hits.clear();
}

void SD_Det::clear() { mapTrackID_Hits.clear(); }

void SD_Det::EndOfEvent(G4HCofThisEvent*)
{
  //   Hits->Clear("C");
  mapTrackID_Hits.clear();
}

G4bool SD_Det::ProcessHits(G4Step* aStep, G4TouchableHistory*)
{

  G4double energ_depos = aStep->GetTotalEnergyDeposit();
  
   //std::cout<<" Current SD :"<<SensitiveDetectorName<<" "<<energ_depos<<std::endl;
  if(energ_depos < 1e-6)
    return true;

  G4StepPoint* preStepPoint = aStep->GetPreStepPoint();
  G4Track* track = aStep->GetTrack();
  //G4String PhysName(aStep->GetPreStepPoint()->GetPhysicalVolume()->GetName());

  G4TouchableHandle touchable = preStepPoint->GetTouchableHandle();
  G4ThreeVector worldPosition = preStepPoint->GetPosition();
  G4ThreeVector localPosition = touchable->GetHistory()->GetTopTransform().
  TransformPoint(worldPosition);

  //const G4TouchableHistory* touchable = dynamic_cast<const G4TouchableHistory*>(aStep->GetPreStepPoint()->GetTouchable());
  G4VPhysicalVolume* currentPhysical = touchable->GetVolume(0);
  G4int copyNo = touchable->GetCopyNumber();
  G4int copyNo1 = -1;
  if(touchable->GetHistoryDepth()>1) copyNo1 = touchable->GetCopyNumber(1);
  G4int copyNo2 = -1;
  if(touchable->GetHistoryDepth()>2) copyNo2 = touchable->GetCopyNumber(2);

  // G4int replicat1 = touchable->GetReplicaNumber(1);
  // G4int replicat = touchable->GetReplicaNumber(0);

  int CurrentTrack = track->GetTrackID();

  auto it_layer = mapTrackID_Hits.find(copyNo);
  if(it_layer == mapTrackID_Hits.end())
    {
      int IdHit        = fHitsCollection->GetSize();
      SingleHit* newHit = new SingleHit(fHCID);

      newHit->Pdg         = track->GetDefinition()->GetPDGEncoding();
      newHit->TrackID     = track->GetTrackID();
      newHit->Edep        = energ_depos;
      newHit->Time        = preStepPoint->GetGlobalTime();
      newHit->TrackLength = track->GetTrackLength();
      newHit->HitPosX     = localPosition.x();
      newHit->HitPosY     = localPosition.y();
      newHit->HitPosZ     = localPosition.z();
  	  newHit->GlobalPosX  = worldPosition.x();
      newHit->GlobalPosY  = worldPosition.y();
      newHit->GlobalPosZ  = worldPosition.z();
      newHit->MomX        = track->GetMomentum().x();
      newHit->MomY        = track->GetMomentum().y();
      newHit->MomZ        = track->GetMomentum().z();
      newHit->LayerID     = copyNo;
      newHit->LayerID1    = copyNo1;
      newHit->LayerID2    = copyNo2;

      mapTrackID_Hits.insert(std::pair<int, std::unordered_map<int, int> >(copyNo, {{CurrentTrack, IdHit}}));

      fHitsCollection->insert(newHit);
    }
  else
    {
      int IdHit = -1;

      auto it_track = it_layer->second.find(CurrentTrack);

      if(it_track != it_layer->second.end())
        IdHit = it_track->second;

      // G4cout<<" --> IdHit:"<<IdHit;

      if(IdHit == -1)
        {
          // G4cout<<" NewHit !"<<G4endl;
          IdHit            = fHitsCollection->GetSize();
          SingleHit* newHit = new SingleHit(fHCID);

          newHit->Pdg         = aStep->GetTrack()->GetDefinition()->GetPDGEncoding();
          newHit->TrackID     = aStep->GetTrack()->GetTrackID();
          newHit->Edep        = energ_depos;
          newHit->Time        = aStep->GetPreStepPoint()->GetGlobalTime();
		  newHit->TrackLength = aStep->GetTrack()->GetTrackLength();
		  newHit->HitPosX     = localPosition.x();
		  newHit->HitPosY     = localPosition.y();
		  newHit->HitPosZ     = localPosition.z();
		  newHit->GlobalPosX  = worldPosition.x();
		  newHit->GlobalPosY  = worldPosition.y();
		  newHit->GlobalPosZ  = worldPosition.z();
		  newHit->MomX        = aStep->GetTrack()->GetMomentum().x();
          newHit->MomY        = aStep->GetTrack()->GetMomentum().y();
          newHit->MomZ        = aStep->GetTrack()->GetMomentum().z();
          newHit->LayerID     = copyNo;
          newHit->LayerID1    = copyNo1;
          newHit->LayerID2    = copyNo2;
          it_layer->second.insert(std::pair<int, int>(CurrentTrack, IdHit));

          fHitsCollection->insert(newHit);
        }
      else
        {
          // G4cout<<" Hit#"<<IdHit<<G4endl;
          SingleHit* CurrentHit = dynamic_cast<SingleHit*>((*fHitsCollection)[IdHit]);
          if(CurrentHit != nullptr)
            {
              CurrentHit->Edep += energ_depos;
            }
        }
    }

  return true;
}
