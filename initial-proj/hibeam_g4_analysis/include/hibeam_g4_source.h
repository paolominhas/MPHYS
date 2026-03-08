//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Fri May  3 11:28:41 2024 by ROOT version 6.30/02
// from TTree hibeam/HIBEAM simulation output
// found on file: test5.root
//////////////////////////////////////////////////////////

#ifndef hibeam_g4_source_h
#define hibeam_g4_source_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include <vector>

class hibeam_g4_source {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   std::vector<int>     *PrimaryTrackID;
   std::vector<int>     *PrimaryPDG;
   std::vector<double>  *PrimaryEkin;
   std::vector<double>  *PrimaryTime;
   std::vector<double>  *PrimaryPosX;
   std::vector<double>  *PrimaryPosY;
   std::vector<double>  *PrimaryPosZ;
   std::vector<double>  *PrimaryMomX;
   std::vector<double>  *PrimaryMomY;
   std::vector<double>  *PrimaryMomZ;
   std::vector<double>  *PrimaryWeight;

   // List of branches
   TBranch        *b_PrimaryTrackID;   //!
   TBranch        *b_PrimaryPDG;   //!
   TBranch        *b_PrimaryEkin;   //!
   TBranch        *b_PrimaryTime;   //!
   TBranch        *b_PrimaryPosX;   //!
   TBranch        *b_PrimaryPosY;   //!
   TBranch        *b_PrimaryPosZ;   //!
   TBranch        *b_PrimaryMomX;   //!
   TBranch        *b_PrimaryMomY;   //!
   TBranch        *b_PrimaryMomZ;   //!
   TBranch        *b_PrimaryWeight;   //!


   hibeam_g4_source();
   hibeam_g4_source(TTree *tree);
   ~hibeam_g4_source();
   void Init(TTree *tree);
   void AddToOutput(TTree *tree);
};

#endif

