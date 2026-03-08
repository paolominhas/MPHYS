// -------------------------------------------------------------
// Implementation of the MultiHits class
// Created by C.Rappold (c.rappold@gsi.de)
//--------------------------------------------------------------

#include "MultiHits.hh"

#include <iostream>

MultiHits::MultiHits()
{
}

MultiHits::~MultiHits() {}
void MultiHits::Print() const
{
	for(int i=0; i<TrackID.size(); i++){
  std::cout << "The TMultiHits #" << TrackID[i] << " Layer:" << LayerID[i] << std::endl;
  std::cout << "Edep:" << Edep[i] << " Time:" << Time[i] << std::endl;
  std::cout << "Position :" << HitPosX[i] << " " << HitPosY[i] << " " << HitPosZ[i] << std::endl;
  std::cout << "Global Position :" << GlobalPosX[i] << " " << GlobalPosY[i] << " " << GlobalPosZ[i] << std::endl;
  std::cout << "Momentum :" << MomX[i] << " " << MomY[i] << " " << MomZ[i] << std::endl;
}
}

void MultiHits::Reset()
{
  TrackID.clear();
  
  LayerID.clear();
  LayerID1.clear();
  LayerID2.clear();
   
  HitPosX.clear();
  HitPosY.clear();
  HitPosZ.clear();
 
  GlobalPosX.clear();
  GlobalPosY.clear();
  GlobalPosZ.clear();

  MomX.clear();
  MomY.clear();
  MomZ.clear();

  Edep.clear();
  Time.clear();
  TrackLength.clear();

  Pdg.clear();
}
