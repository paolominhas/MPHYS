// -----------------------------------------------------
// Definition of the MultiHits class
// Created by C.Rappold (c.rappold@gsi.de)
//------------------------------------------------------

#ifndef MULTIHITS_H
#define MULTIHITS_H

#include <vector>
#include <string>

class MultiHits
{
public:
  MultiHits();
  ~MultiHits();

  void Print() const;
  void Reset();

  std::vector<int> TrackID;
  
  std::vector<int> LayerID;
  std::vector<int> LayerID1;
  std::vector<int> LayerID2;
   
  std::vector<double> HitPosX;
  std::vector<double> HitPosY;
  std::vector<double> HitPosZ;

  std::vector<double> GlobalPosX;
  std::vector<double> GlobalPosY;
  std::vector<double> GlobalPosZ;

  std::vector<double> MomX;
  std::vector<double> MomY;
  std::vector<double> MomZ;

  std::vector<double> Edep;
  std::vector<double> Time;
  std::vector<double> TrackLength;

  std::vector<int> Pdg;

};

#endif // MULTIHITS_H
