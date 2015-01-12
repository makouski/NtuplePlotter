#include<iostream>

#include"EventTree.h"
#include"Selector.h"
#include"EventPick.h"

void fillCategory(EventTree* tree, TH1F* hist){
	int EleP = 0;
	int EleM = 0;
	int MuP = 0;
	int MuM = 0;
	int TauP = 0;
	int TauM = 0;
	for( int mcI = 0; mcI < tree->nMC_; ++mcI){
		if(abs(tree->mcMomPID->at(mcI))==24 && tree->mcParentage->at(mcI)==10){
			if( tree->mcPID->at(mcI) == 11 ) EleP = 1;
			if( tree->mcPID->at(mcI) == -11 ) EleM = 1;
			if( tree->mcPID->at(mcI) == 13 ) MuP = 1;
			if( tree->mcPID->at(mcI) == -13 ) MuM = 1;
			if( tree->mcPID->at(mcI) == 15) TauP = 1;
			if( tree->mcPID->at(mcI) == -15) TauM = 1;
		}
	}
	hist->Fill(1.0); // Total
	int nEle = EleP + EleM;
	int nMu = MuP + MuM;
	int nTau = TauP + TauM;
	if( nEle + nMu + nTau == 0) hist->Fill(2.0); // All Had
	if( nEle + nMu + nTau == 1) hist->Fill(3.0); // Single Lepton
	if( nEle + nMu + nTau == 2) hist->Fill(4.0); // Di Lepton
	return;
}

int main(int ac, char** av){
	if(ac < 2){
		std::cout << "usage: ./signalAcceptance inputFile[s]" << std::endl;
		return -1;
	}

	TH1F* allCategory = new TH1F("allCategory","all Category",11,0.5,11.5);
	TH1F* preselCategory = new TH1F("preselCategory","presel Category",11,0.5,11.5);
	TH1F* photonCategory = new TH1F("photonCategory"," photon Category",11,0.5,11.5);
	
	// object selector
	Selector* selectorLoose = new Selector();
	// create event selectors here
	EventPick* evtPickLoose = new EventPick("LoosePhotonID");
	
	EventTree* tree = new EventTree(ac-1, av+1);
	double PUweight = 1.0;
	
	Long64_t nEntr = tree->GetEntries();
	for(Long64_t entry=0; entry<nEntr; entry++){
		if(entry%10000 == 0) std::cout << "processing entry " << entry << " out of " << nEntr << std::endl;
		tree->GetEntry(entry);
		
		selectorLoose->process_objects(tree);
		evtPickLoose->process_event(tree, selectorLoose, PUweight);

		// fill the histograms
		fillCategory(tree, allCategory);
		if(evtPickLoose->passPreSel) fillCategory(tree, preselCategory);
		if(evtPickLoose->passAll) fillCategory(tree, photonCategory);
	}
	
	evtPickLoose->print_cutflow();

	// write histograms
	TFile outFile("signalAcc.root","RECREATE");
	allCategory->SetDirectory(outFile.GetDirectory(""));
	allCategory->Write();
	allCategory->SetDirectory(0);

	preselCategory->SetDirectory(outFile.GetDirectory(""));
	preselCategory->Write();
	preselCategory->SetDirectory(0);
	
	photonCategory->SetDirectory(outFile.GetDirectory(""));
	photonCategory->Write();
	photonCategory->SetDirectory(0);
	
	outFile.Close();

	delete tree;
	return 0;
}

