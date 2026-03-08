#include <stdio.h>
#include <iostream>
#include <vector>
#include "TFile.h"
#include "TTree.h"

#include "SECHit.hh"
#include "ProtoTPCHit.hh"
#include "TPCHit.hh"
#include "CVHit.hh"
#include "SimpleSciHit.hh"
#include "ScintHit.hh"
#include "ScintAna.hh"
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
	hibeam_g4 tpc, prototpc;
	hibeam_g4 cv, cvs;
	hibeam_g4 sci, scint;
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
		prototpc.Init(tree,"ProtoTPC");
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
	if(tree->GetListOfBranches()->FindObject("Scintillator_TrackID")){
		printf("The input file contains generic scintillator data.\n");
		sci.Init(tree,"Scintillator");
	}
	if(tree->GetListOfBranches()->FindObject("Sci_bar_TrackID")){
		printf("The input file contains scintillator prototype data.\n");
		scint.Init(tree,"Sci_bar");
	}


	SECHit *secH = nullptr;
	TPCHit *tpcH = nullptr;
	ProtoTPCHit *prototpcH = nullptr;
	CVHit *cvH = nullptr;
	SimpleSciHit *sciH = nullptr;
	ScintHit *scintH = nullptr;
	ScintAna *scintA = nullptr;
	
	TFile *fout = TFile::Open(outname,"RECREATE");
	TTree *newtree = new TTree("hibeam", "Processed HIBEAM simulation output");

	if(source.fChain!=0){
		source.AddToOutput(newtree);
	}
	if(target.fChain!=0){
		target.AddToOutput(newtree,"target");
	}
	if(sec.size()>0){
		newtree->Branch("SEC",&secH,8000,1);
	}
	if(tpc.fChain!=0){
		newtree->Branch("TPC",&tpcH,8000,1);
	}
	if(prototpc.fChain!=0){
		//prototpc.AddToOutput(newtree,"proto");
		newtree->Branch("ProtoTPC",&prototpcH,8000,1);
	}
	if(cv.fChain!=0||cvs.fChain!=0){
		newtree->Branch("CV",&cvH,8000,1);
	}
	if(sci.fChain!=0){
		newtree->Branch("Scintillator",&sciH,8000,1);
	}
	if(scint.fChain!=0){
		newtree->Branch("Scint",&scintH,8000,1);
		newtree->Branch("ScintAna",&scintA,8000,1);
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
			secH->clear();
			for (int is=0; is<17; is++){
				for (int iv=0; iv<sec[is].TrackID->size(); iv++){
					secH->push_back(is, sec[is].LayerID->at(iv), sec[is].EDep->at(iv), sec[is].Time->at(iv), sec[is].Position_X->at(iv), sec[is].Position_Y->at(iv), sec[is].Position_Z->at(iv));
				}
			}
			secH->add_up();
			secH->sort();
		}

		if(tpc.fChain!=0){
			tpcH->clear();
			for (int iv=0; iv<tpc.TrackID->size(); iv++){
				tpcH->push_back(tpc.EDep->at(iv), tpc.Time->at(iv), tpc.Position_X->at(iv), tpc.Position_Y->at(iv), tpc.Position_Z->at(iv));
			}
			tpcH->add_up();
		}
		if(prototpc.fChain!=0){
			prototpcH->clear();
			for (int iv=0; iv<prototpc.TrackID->size(); iv++){
				prototpcH->push_back(prototpc.EDep->at(iv), prototpc.Time->at(iv), prototpc.Position_X->at(iv), prototpc.Position_Y->at(iv), prototpc.Position_Z->at(iv), prototpc.LayerID->at(iv));
			}
			prototpcH->swap_xy();
			prototpcH->shift();
			prototpcH->assign_pads();
			prototpcH->sum_up();
			prototpcH->clean_up();
		}
		if(cv.fChain!=0||cvs.fChain!=0){
			cvH->clear();
			for (int iv=0; iv<cv.TrackID->size(); iv++){
				cvH->push_back(cv.LayerID->at(iv), cv.LayerID1->at(iv), cv.LayerID2->at(iv), cv.EDep->at(iv), cv.Time->at(iv), cv.Position_X->at(iv), cv.Position_Y->at(iv), cv.Position_Z->at(iv));
			}
			for (int iv=0; iv<cvs.TrackID->size(); iv++){
				cvH->push_back(cvs.LayerID->at(iv), cvs.LayerID1->at(iv), cvs.LayerID2->at(iv), cvs.EDep->at(iv), cvs.Time->at(iv), cvs.Position_X->at(iv), cvs.Position_Y->at(iv), cvs.Position_Z->at(iv));
			}
			cvH->add_up();
			cvH->sort();
		}

		if(sci.fChain!=0){
			sciH->clear();
			for (int iv=0; iv<sci.TrackID->size(); iv++){
				sciH->push_back(sci.LayerID->at(iv), sci.EDep->at(iv), sci.Time->at(iv), sci.Position_X->at(iv), sci.Position_Y->at(iv), sci.Position_Z->at(iv));
			}
			sciH->add_up();
			sciH->sort();
		}

		if(scint.fChain!=0){
			scintH->clear();
			scintA->clear();
			for (int iv=0; iv<scint.TrackID->size(); iv++){
				scintH->push_back(scint.TrackID->at(iv), scint.LayerID2->at(iv), scint.LayerID->at(iv), scint.LayerID1->at(iv), scint.EDep->at(iv), scint.Time->at(iv), scint.Position_X->at(iv), scint.Position_Y->at(iv), scint.Position_Z->at(iv));
			}
			scintH->add_up();
			scintH->sort();
			scintA->response(scintH);
			scintA->add_up();
			scintA->sort();
			scintA->sumE();
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
