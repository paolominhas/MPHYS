// -----------------------------------------------------
// Definition of the PrimaryEvent class
// Created by C.Rappold (c.rappold@gsi.de)
//------------------------------------------------------

#ifndef PRIMARYEVENT_H
#define PRIMARYEVENT_H

#include <fstream>
#include <iostream>
#include <string>
#include <vector>

/********************************************************************/
class PrimaryEvent
{
public: // Primary decay vertex
  //
  // beam composition
  std::vector<int> PDG;
  std::vector<double> Ekin;
  std::vector<double> Time;
  std::vector<double> PosX;
  std::vector<double> PosY;
  std::vector<double> PosZ;
  std::vector<double> MomX;
  std::vector<double> MomY;
  std::vector<double> MomZ;
  std::vector<int> TrackID;
  std::vector<double> Weight;
  //***********************************************************

public:
  PrimaryEvent();
  ~PrimaryEvent();
  //
  void Display();
  void Zero();
  /*------------------------------------------------------------*/
};

#endif // PRIMARYEVENT_H
