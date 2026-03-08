// -------------------------------------------------
// Definition of the SingleHit class
// Created by C.Rappold (c.rappold@gsi.de)
//--------------------------------------------------

#ifndef SINGLEHIT_H
#define SINGLEHIT_H

#include "G4Allocator.hh"
#include "G4THitsCollection.hh"
#include "G4VHit.hh"

//#include <string>

class SingleHit : public G4VHit
{
public:
  SingleHit();
  explicit SingleHit(G4int z);
  SingleHit(const SingleHit& hit);
  virtual ~SingleHit();

  const SingleHit& operator=(const SingleHit& right);
  int operator==(const SingleHit& right) const;

  inline void* operator new(size_t);
  inline void operator delete(void* aHit);

  virtual void Draw();
  // virtual const std::map<G4String,G4AttDef>* GetAttDefs() const;
  // virtual std::vector<G4AttValue>* CreateAttValues() const;
  virtual void Print();

  G4int TrackID; ///< Track Id

  G4double HitPosX;
  G4double HitPosY;
  G4double HitPosZ;

  G4double GlobalPosX;
  G4double GlobalPosY;
  G4double GlobalPosZ;

  G4double MomX;
  G4double MomY;
  G4double MomZ;

  G4double Edep;
  G4double Time;
  G4double TrackLength;

  G4int LayerID;
  G4int LayerID1;
  G4int LayerID2;
  G4int Pdg;
};
typedef G4THitsCollection<SingleHit> SingleHitsCollection;

extern G4ThreadLocal G4Allocator<SingleHit>* SingleHitAllocator;

inline void* SingleHit::operator new(size_t)
{
  if(!SingleHitAllocator)
    SingleHitAllocator = new G4Allocator<SingleHit>;
  return (void*)SingleHitAllocator->MallocSingle();
}

inline void SingleHit::operator delete(void* aHit) { SingleHitAllocator->FreeSingle((SingleHit*)aHit); }

#endif // SINGLEHIT_H
