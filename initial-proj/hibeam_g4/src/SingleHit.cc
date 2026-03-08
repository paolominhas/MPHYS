// -------------------------------------------------------
// Implementation of the SingleHit class
// Created by C.Rappold (c.rappold@gsi.de)
//--------------------------------------------------------

#include "SingleHit.hh"

#include "G4Circle.hh"
#include "G4Colour.hh"
#include "G4SystemOfUnits.hh"
#include "G4UnitsTable.hh"
#include "G4VVisManager.hh"
#include "G4VisAttributes.hh"
#include "G4ios.hh"

G4ThreadLocal G4Allocator<SingleHit>* SingleHitAllocator;

SingleHit::SingleHit()
    : TrackID(-1), HitPosX(-9999.), HitPosY(-9999.), HitPosZ(-9999.), GlobalPosX(-9999.), GlobalPosY(-9999.),
      GlobalPosZ(-9999.), MomX(-9999.), MomY(-9999.), MomZ(-9999.), Edep(-9999.), Time(-9999.),
      TrackLength(-9999.), LayerID(-1), LayerID1(-1), LayerID2(-1), Pdg(-999999)

{
}

SingleHit::SingleHit(G4int z)
    : TrackID(-1), HitPosX(-9999.), HitPosY(-9999.), HitPosZ(-9999.), GlobalPosX(-9999.), GlobalPosY(-9999.),
      GlobalPosZ(-9999.), MomX(-9999.), MomY(-9999.), MomZ(-9999.), Edep(-9999.), Time(-9999.),
      TrackLength(-9999.), LayerID(z), LayerID1(-1), LayerID2(-1), Pdg(-999999)

{
}

SingleHit::SingleHit(const SingleHit& hit) : G4VHit()
{
  TrackID     = hit.TrackID;
  HitPosX     = hit.HitPosX;
  HitPosY     = hit.HitPosY;
  HitPosZ     = hit.HitPosZ;
  GlobalPosX  = hit.GlobalPosX;
  GlobalPosY  = hit.GlobalPosY;
  GlobalPosZ  = hit.GlobalPosZ;
  MomX        = hit.MomX;
  MomY        = hit.MomY;
  MomZ        = hit.MomZ;
  Edep        = hit.Edep;
  Time        = hit.Time;
  TrackLength = hit.TrackLength;
  LayerID     = hit.LayerID;
  LayerID1    = hit.LayerID1;
  LayerID2    = hit.LayerID2;
  Pdg         = hit.Pdg;
}

const SingleHit& SingleHit::operator=(const SingleHit& hit)
{
  TrackID     = hit.TrackID;
  HitPosX     = hit.HitPosX;
  HitPosY     = hit.HitPosY;
  HitPosZ     = hit.HitPosZ;
  GlobalPosX  = hit.GlobalPosX;
  GlobalPosY  = hit.GlobalPosY;
  GlobalPosZ  = hit.GlobalPosZ;
  MomX        = hit.MomX;
  MomY        = hit.MomY;
  MomZ        = hit.MomZ;
  Edep        = hit.Edep;
  Time        = hit.Time;
  TrackLength = hit.TrackLength;
  LayerID     = hit.LayerID;
  LayerID1    = hit.LayerID1;
  LayerID2    = hit.LayerID2;
  Pdg         = hit.Pdg;

  return *this;
}

int SingleHit::operator==(const SingleHit& hit) const { return LayerID == hit.LayerID && TrackID == hit.TrackID; }

void SingleHit::Draw()
{
  G4VVisManager* pVVisManager = G4VVisManager::GetConcreteInstance();
  if(pVVisManager)
    {
      G4Circle circle(G4Point3D(HitPosX, HitPosY, HitPosZ));
      circle.SetScreenSize(2);
      circle.SetFillStyle(G4Circle::filled);
      G4Colour colour(1., 1., 0.);
      G4VisAttributes attribs(colour);
      circle.SetVisAttributes(attribs);
      pVVisManager->Draw(circle);
    }
}

SingleHit::~SingleHit() {}

void SingleHit::Print()
{
  std::cout << "The Hit #" << TrackID << std::endl;
  std::cout << "Edep:" << Edep << " Time:" << Time << std::endl;
  std::cout << "Position :" << HitPosX << " " << HitPosY << " " << HitPosZ << std::endl;
  std::cout << "Global Position :" << GlobalPosX << " " << GlobalPosY << " " << GlobalPosZ << std::endl;
  std::cout << "Momentum :" << MomX << " " << MomY << " " << MomZ << std::endl;
}
