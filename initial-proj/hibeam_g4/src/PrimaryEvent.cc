// -------------------------------------------------------------
// Implementation of the PrimaryEvent class
// Created by C.Rappold (c.rappold@gsi.de)
//--------------------------------------------------------------

#include "PrimaryEvent.hh"

/********************************************************************/
PrimaryEvent::PrimaryEvent() { Zero(); }

/********************************************************************/
PrimaryEvent::~PrimaryEvent() {}

/********************************************************************/
void PrimaryEvent::Zero()
{
  // Primary particle information
  PDG.resize(0);
  Ekin.resize(0);
  Time.resize(0);
  PosX.resize(0);
  PosY.resize(0);
  PosZ.resize(0);
  MomX.resize(0);
  MomY.resize(0);
  MomZ.resize(0);
  TrackID.resize(0);
  Weight.resize(0);
}

/********************************************************************/
void PrimaryEvent::Display()
{
  std::cout << "The PrimaryEvent class\n";
  //---------------------------------------------------------------
  // Primary deecay vertex
  std::cout << "Primary Event composition: \n";
  for(size_t i = 0; i < PDG.size(); i++)
    {
      std::cout << " Particle PDG Code: " << PDG[i] << "\n";
      std::cout << " Particle energy: " << Ekin[i] << "\n";
      std::cout << " Particle time: " << Time[i] << "\n";
      std::cout << " Particle position: (" << PosX[i] << ";";
      std::cout << PosY[i] << ";";
      std::cout << PosZ[i] << ")\n";
      std::cout << " Particle momentum: (" << MomX[i] << ";";
      std::cout << MomY[i] << ";";
      std::cout << MomZ[i] << ")\n";
      std::cout << " Particle weight: " << Weight[i] << "\n";
    }
}

/********************************************************************/
