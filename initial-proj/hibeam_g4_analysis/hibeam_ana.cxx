#include <stdio.h>
#include <iostream>
#include <vector>
#include "TFile.h"
#include "TTree.h"

#include "SECHit.hh"
#include "ProtoTPCHit.hh"
#include "TPCHit.hh"
#include "CVHit.hh"
#include "TrigSciHit.hh"
#include "GenSciHit.hh"
#include "HRDSciHit.hh"
#include "hibeam_g4.h"
#include "hibeam_g4_source.h"

int main(int argc, char *argv[])
{
	char *inname=NULL;
	char *outname=NULL;

	if (argc > 1){
		for(int i=0; i<argc; i++){
			if(strncmp(argv[i],"--in=",5)==0){
				inname = argv[i]+5;
			}
			else if(strncmp(argv[i],"--out=",6)==0){
				outname = argv[i]+6;
			}
		}
	}

	if(inname==NULL){
		printf("No input filename given!\n");
		printf("Usage: filter --in=/path/to/input.root --out=/path/to/output.root\n");
		return 0;
	}
	else{
		printf("Reading from file %s\n", inname);
	}
	if(outname==0){
		printf("No output filename given!\n");
		printf("Usage: filter --in=/path/to/input.root --out=/path/to/output.root\n");
		return 0;
	}
	else{
		printf("Writing to file %s\n", outname);
	}

	TFile *fin = TFile::Open(inname,"READ");
	TTree *tree = (TTree*) fin->Get("hibeam");

	hibeam_g4_source source; 
	hibeam_g4 target;
	hibeam_g4 tpc, protoTPC;
	hibeam_g4 cv, cvs;
	hibeam_g4 trigSci, genSci, hrdSci;
	std::vector<hibeam_g4> sec;
	sec.clear();

	if(tree->GetListOfBranches()->FindObject("PrimaryTrackID")){
		printf("The input file contains source data.\n");
		source.Init(tree);
	}
	if(tree->GetListOfBranches()->FindObject("TARGET_TrackID")){
		printf("The input file contains target data!\n");
		target.Init(tree,"TARGET");
	}
	if(tree->GetListOfBranches()->FindObject("TPC_TrackID")){
		printf("The input file contains TPC data.\n");
		tpc.Init(tree,"TPC");
	}
	if(tree->GetListOfBranches()->FindObject("ProtoTPC_TrackID")){
		printf("The input file contains TPC prototype data.\n");
		protoTPC.Init(tree,"ProtoTPC");
	}
	if(tree->GetListOfBranches()->FindObject("SECE1_TrackID")){
		printf("The input file contains SEC data.\n");
		for(int i=0; i<17; i++){
			hibeam_g4 secTmp(tree,Form("SECE%d",i));
			sec.push_back(secTmp);
		}
	}
	if(tree->GetListOfBranches()->FindObject("CV_bar_TrackID")){
		printf("The input file contains cosmic veto data.\n");
		cv.Init(tree,"CV_bar");
		cvs.Init(tree,"CV_shortbar");
	}
	if(tree->GetListOfBranches()->FindObject("trigSci_TrackID")){
		printf("The input file contains trigger scintillator data.\n");
		trigSci.Init(tree,"trigSci");
	}
	if(tree->GetListOfBranches()->FindObject("Scintillator_TrackID")){
		printf("The input file contains generic scintillator data.\n");
		genSci.Init(tree,"Scintillator");
	}
	if(tree->GetListOfBranches()->FindObject("Sci_bar_TrackID")){
		printf("The input file contains scintillator prototype data.\n");
		hrdSci.Init(tree,"Sci_bar");
	}
	else if(tree->GetListOfBranches()->FindObject("HRDBar_TrackID")){
		printf("The input file contains scintillator prototype data.\n");
		hrdSci.Init(tree,"HRDBar");
	}

	SECHit *hitSEC = nullptr;
	TPCHit *hitTPC = nullptr;
	ProtoTPCHit *hitProtoTPC = nullptr;
	CVHit *hitCV = nullptr;
	TrigSciHit *hitTrigSci = nullptr;
	GenSciHit *hitGenSci = nullptr;
	HRDSciHit *hitHRDSci = nullptr;
	
	TFile *fout = TFile::Open(outname,"RECREATE");
	TTree *newtree = new TTree("hibeam", "Processed HIBEAM simulation output");

	if(source.fChain!=0){
		source.AddToOutput(newtree);
	}
	if(target.fChain!=0){
		target.AddToOutput(newtree,"target");
	}
	if(sec.size()>0){
		newtree->Branch("SEC",&hitSEC,8000,1);
	}
	if(tpc.fChain!=0){
		newtree->Branch("TPC",&hitTPC,8000,1);
	}
	if(protoTPC.fChain!=0){
		//protoTPC.AddToOutput(newtree,"proto");
		newtree->Branch("ProtoTPC",&hitProtoTPC,8000,1);
	}
	if(cv.fChain!=0||cvs.fChain!=0){
		newtree->Branch("CV",&hitCV,8000,1);
	}
	if(trigSci.fChain!=0){
		newtree->Branch("trigSci",&hitTrigSci,8000,1);
	}
	if(genSci.fChain!=0){
		newtree->Branch("Scintillator",&hitGenSci,8000,1);
	}
	if(hrdSci.fChain!=0){
		newtree->Branch("HRD",&hitHRDSci,8000,1);
	}


	printf("Begin analysis.\n");
	int prev_trackid=-1;
	Long64_t nentries = tree->GetEntriesFast();
	for (Long64_t jentry=0; jentry<nentries;jentry++) {
		Long64_t ientry = tree->LoadTree(jentry);
		if (ientry < 0) break;
		tree->GetEntry(jentry);
		prev_trackid=-1;

		if(sec.size()>0){
			hitSEC->clear();
			for (int is=0; is<17; is++){
				for (int iv=0; iv<sec[is].TrackID->size(); iv++){
					hitSEC->push_back(is, sec[is].LayerID->at(iv), sec[is].EDep->at(iv), sec[is].Time->at(iv), sec[is].Position_X->at(iv), sec[is].Position_Y->at(iv), sec[is].Position_Z->at(iv));
				}
			}
			hitSEC->add_up();
			hitSEC->sort();
		}

		if(tpc.fChain!=0){
			hitTPC->clear();
			for (int iv=0; iv<tpc.TrackID->size(); iv++){
				hitTPC->assign_pads(tpc.EDep->at(iv), tpc.Time->at(iv), tpc.Position_X->at(iv), tpc.Position_Y->at(iv), tpc.Position_Z->at(iv));
			}
			hitTPC->add_up();
		}
		if(protoTPC.fChain!=0){
			hitProtoTPC->clear();
			for (int iv=0; iv<protoTPC.TrackID->size(); iv++){
				hitProtoTPC->push_back(protoTPC.EDep->at(iv), protoTPC.Time->at(iv), protoTPC.Position_X->at(iv), protoTPC.Position_Y->at(iv), protoTPC.Position_Z->at(iv), protoTPC.LayerID->at(iv));
			}
			hitProtoTPC->process();
		}
		if(cv.fChain!=0||cvs.fChain!=0){
			hitCV->clear();
			for (int iv=0; iv<cv.TrackID->size(); iv++){
				hitCV->push_back(cv.LayerID->at(iv), cv.LayerID1->at(iv), cv.LayerID2->at(iv), cv.EDep->at(iv), cv.Time->at(iv), cv.Position_X->at(iv), cv.Position_Y->at(iv), cv.Position_Z->at(iv));
			}
			for (int iv=0; iv<cvs.TrackID->size(); iv++){
				hitCV->push_back(cvs.LayerID->at(iv), cvs.LayerID1->at(iv), cvs.LayerID2->at(iv), cvs.EDep->at(iv), cvs.Time->at(iv), cvs.Position_X->at(iv), cvs.Position_Y->at(iv), cvs.Position_Z->at(iv));
			}
			hitCV->add_up();
			hitCV->sort();
		}

		if(trigSci.fChain!=0){
			hitTrigSci->clear();
			for (int iv=0; iv<trigSci.TrackID->size(); iv++){
				hitTrigSci->push_back(2*trigSci.LayerID->at(iv)+trigSci.LayerID1->at(iv), trigSci.EDep->at(iv), trigSci.Time->at(iv));
			}
			hitTrigSci->process();
		}

		if(genSci.fChain!=0){
			hitGenSci->clear();
			for (int iv=0; iv<genSci.TrackID->size(); iv++){
				hitGenSci->push_back(genSci.LayerID->at(iv), genSci.EDep->at(iv), genSci.Time->at(iv), genSci.Position_X->at(iv), genSci.Position_Y->at(iv), genSci.Position_Z->at(iv));
			}
			hitGenSci->process();
		}

		if(hrdSci.fChain!=0){
			hitHRDSci->clear();
			for (int iv=0; iv<hrdSci.TrackID->size(); iv++){
				hitHRDSci->push_back(hrdSci.TrackID->at(iv), hrdSci.LayerID2->at(iv), hrdSci.LayerID->at(iv), hrdSci.LayerID1->at(iv), hrdSci.EDep->at(iv), hrdSci.Time->at(iv), hrdSci.Position_X->at(iv), hrdSci.Position_Y->at(iv), hrdSci.Position_Z->at(iv));
			}
			hitHRDSci->process();
		}

		if((jentry%10000)==0){
			printf("%lld events read.\n",jentry);
		}

		newtree->Fill();
	}
	newtree->Print();
	newtree->AutoSave();

	printf("Finished. %lld events read in total.\n",nentries);
	fin->Close();
	fout->Close();

	delete fin;
	delete fout;
	return nentries;
}
