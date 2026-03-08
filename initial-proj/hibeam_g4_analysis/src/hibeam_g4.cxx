#include "hibeam_g4.h"

hibeam_g4::hibeam_g4() : fChain(0) 
{
}


hibeam_g4::hibeam_g4(TTree *tree, std::string branchName) : fChain(0) 
{
	Init(tree,branchName);
}

hibeam_g4::~hibeam_g4()
{
   //if (!fChain) return;
   //delete fChain->GetCurrentFile();
}


void hibeam_g4::Init(TTree *tree, std::string branchName)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   TrackID = 0;
   LayerID = 0;
   LayerID1 = 0;
   LayerID2 = 0;
   PDG = 0;
   EDep = 0;
   Time = 0;
   TrackLength = 0;
   Position_X = 0;
   Position_Y = 0;
   Position_Z = 0;
   GlobalPosition_X = 0;
   GlobalPosition_Y = 0;
   GlobalPosition_Z = 0;
   Momentum_X = 0;
   Momentum_Y = 0;
   Momentum_Z = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress(Form("%s_TrackID",branchName.c_str()), &TrackID, &b_TrackID);
   fChain->SetBranchAddress(Form("%s_LayerID",branchName.c_str()), &LayerID, &b_LayerID);
   fChain->SetBranchAddress(Form("%s_LayerID1",branchName.c_str()), &LayerID1, &b_LayerID1);
   fChain->SetBranchAddress(Form("%s_LayerID2",branchName.c_str()), &LayerID2, &b_LayerID2);
   fChain->SetBranchAddress(Form("%s_PDG",branchName.c_str()), &PDG, &b_PDG);
   fChain->SetBranchAddress(Form("%s_EDep",branchName.c_str()), &EDep, &b_EDep);
   fChain->SetBranchAddress(Form("%s_Time",branchName.c_str()), &Time, &b_Time);
   fChain->SetBranchAddress(Form("%s_TrackLength",branchName.c_str()), &TrackLength, &b_TrackLength);
   fChain->SetBranchAddress(Form("%s_Position_X",branchName.c_str()), &Position_X, &b_Position_X);
   fChain->SetBranchAddress(Form("%s_Position_Y",branchName.c_str()), &Position_Y, &b_Position_Y);
   fChain->SetBranchAddress(Form("%s_Position_Z",branchName.c_str()), &Position_Z, &b_Position_Z);
   fChain->SetBranchAddress(Form("%s_GlobalPosition_X",branchName.c_str()), &GlobalPosition_X, &b_GlobalPosition_X);
   fChain->SetBranchAddress(Form("%s_GlobalPosition_Y",branchName.c_str()), &GlobalPosition_Y, &b_GlobalPosition_Y);
   fChain->SetBranchAddress(Form("%s_GlobalPosition_Z",branchName.c_str()), &GlobalPosition_Z, &b_GlobalPosition_Z);
   fChain->SetBranchAddress(Form("%s_Momentum_X",branchName.c_str()), &Momentum_X, &b_Momentum_X);
   fChain->SetBranchAddress(Form("%s_Momentum_Y",branchName.c_str()), &Momentum_Y, &b_Momentum_Y);
   fChain->SetBranchAddress(Form("%s_Momentum_Z",branchName.c_str()), &Momentum_Z, &b_Momentum_Z);
}


void hibeam_g4::AddToOutput(TTree *tree, std::string branchName)
{
   tree->Branch(Form("%s_TrackID",branchName.c_str()), &TrackID);
   tree->Branch(Form("%s_LayerID",branchName.c_str()), &LayerID);
   tree->Branch(Form("%s_LayerID1",branchName.c_str()), &LayerID1);
   tree->Branch(Form("%s_LayerID2",branchName.c_str()), &LayerID2);
   tree->Branch(Form("%s_PDG",branchName.c_str()), &PDG);
   tree->Branch(Form("%s_EDep",branchName.c_str()), &EDep);
   tree->Branch(Form("%s_Time",branchName.c_str()), &Time);
   tree->Branch(Form("%s_TrackLength",branchName.c_str()), &TrackLength);
   tree->Branch(Form("%s_Position_X",branchName.c_str()), &Position_X);
   tree->Branch(Form("%s_Position_Y",branchName.c_str()), &Position_Y);
   tree->Branch(Form("%s_Position_Z",branchName.c_str()), &Position_Z);
   tree->Branch(Form("%s_GlobalPosition_X",branchName.c_str()), &GlobalPosition_X);
   tree->Branch(Form("%s_GlobalPosition_Y",branchName.c_str()), &GlobalPosition_Y);
   tree->Branch(Form("%s_GlobalPosition_Z",branchName.c_str()), &GlobalPosition_Z);
   tree->Branch(Form("%s_Momentum_X",branchName.c_str()), &Momentum_X);
   tree->Branch(Form("%s_Momentum_Y",branchName.c_str()), &Momentum_Y);
   tree->Branch(Form("%s_Momentum_Z",branchName.c_str()), &Momentum_Z);
}
