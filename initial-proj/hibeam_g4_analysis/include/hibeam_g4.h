//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Fri May  3 11:28:41 2024 by ROOT version 6.30/02
// from TTree hibeam/HIBEAM simulation output
// found on file: test5.root
//////////////////////////////////////////////////////////

#ifndef hibeam_g4_h
#define hibeam_g4_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include <vector>

class hibeam_g4 {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   std::vector<int>     *TrackID;
   std::vector<int>     *LayerID;
   std::vector<int>     *LayerID1;
   std::vector<int>     *LayerID2;
   std::vector<int>     *PDG;
   std::vector<double>  *EDep;
   std::vector<double>  *Time;
   std::vector<double>  *TrackLength;
   std::vector<double>  *Position_X;
   std::vector<double>  *Position_Y;
   std::vector<double>  *Position_Z;
   std::vector<double>  *GlobalPosition_X;
   std::vector<double>  *GlobalPosition_Y;
   std::vector<double>  *GlobalPosition_Z;
   std::vector<double>  *Momentum_X;
   std::vector<double>  *Momentum_Y;
   std::vector<double>  *Momentum_Z;


   // List of branches
   TBranch        *b_TrackID;   //!
   TBranch        *b_LayerID;   //!
   TBranch        *b_LayerID1;   //!
   TBranch        *b_LayerID2;   //!
   TBranch        *b_PDG;   //!
   TBranch        *b_EDep;   //!
   TBranch        *b_Time;   //!
   TBranch        *b_TrackLength;   //!
   TBranch        *b_Position_X;   //!
   TBranch        *b_Position_Y;   //!
   TBranch        *b_Position_Z;   //!
   TBranch        *b_GlobalPosition_X;   //!
   TBranch        *b_GlobalPosition_Y;   //!
   TBranch        *b_GlobalPosition_Z;   //!
   TBranch        *b_Momentum_X;   //!
   TBranch        *b_Momentum_Y;   //!
   TBranch        *b_Momentum_Z;   //!


   hibeam_g4();
   hibeam_g4(TTree *tree, std::string branchName);
   ~hibeam_g4();
   void Init(TTree *tree, std::string branchName);
   void AddToOutput(TTree *tree, std::string branchName);
};

#endif

