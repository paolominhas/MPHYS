#include "hibeam_g4_source.h"

hibeam_g4_source::hibeam_g4_source() : fChain(0) 
{
}


hibeam_g4_source::hibeam_g4_source(TTree *tree) : fChain(0) 
{
   Init(tree);
}


hibeam_g4_source::~hibeam_g4_source()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}


void hibeam_g4_source::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   PrimaryTrackID = 0;
   PrimaryPDG = 0;
   PrimaryEkin = 0;
   PrimaryTime = 0;
   PrimaryPosX = 0;
   PrimaryPosY = 0;
   PrimaryPosZ = 0;
   PrimaryMomX = 0;
   PrimaryMomY = 0;
   PrimaryMomZ = 0;
   PrimaryWeight = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("PrimaryTrackID", &PrimaryTrackID, &b_PrimaryTrackID);
   fChain->SetBranchAddress("PrimaryPDG", &PrimaryPDG, &b_PrimaryPDG);
   fChain->SetBranchAddress("PrimaryEkin", &PrimaryEkin, &b_PrimaryEkin);
   fChain->SetBranchAddress("PrimaryTime", &PrimaryTime, &b_PrimaryTime);
   fChain->SetBranchAddress("PrimaryPosX", &PrimaryPosX, &b_PrimaryPosX);
   fChain->SetBranchAddress("PrimaryPosY", &PrimaryPosY, &b_PrimaryPosY);
   fChain->SetBranchAddress("PrimaryPosZ", &PrimaryPosZ, &b_PrimaryPosZ);
   fChain->SetBranchAddress("PrimaryMomX", &PrimaryMomX, &b_PrimaryMomX);
   fChain->SetBranchAddress("PrimaryMomY", &PrimaryMomY, &b_PrimaryMomY);
   fChain->SetBranchAddress("PrimaryMomZ", &PrimaryMomZ, &b_PrimaryMomZ);
   fChain->SetBranchAddress("PrimaryWeight", &PrimaryWeight, &b_PrimaryWeight);
}


void hibeam_g4_source::AddToOutput(TTree *tree)
{
   tree->Branch("PrimaryTrackID", &PrimaryTrackID);
   tree->Branch("PrimaryPDG", &PrimaryPDG);
   tree->Branch("PrimaryEkin", &PrimaryEkin);
   tree->Branch("PrimaryTime", &PrimaryTime);
   tree->Branch("PrimaryPosX", &PrimaryPosX);
   tree->Branch("PrimaryPosY", &PrimaryPosY);
   tree->Branch("PrimaryPosZ", &PrimaryPosZ);
   tree->Branch("PrimaryMomX", &PrimaryMomX);
   tree->Branch("PrimaryMomY", &PrimaryMomY);
   tree->Branch("PrimaryMomZ", &PrimaryMomZ);
   tree->Branch("PrimaryWeight", &PrimaryWeight);
}
