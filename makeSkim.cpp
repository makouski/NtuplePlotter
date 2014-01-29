#include<iostream>

#include"EventTree.h"
#include"Selector.h"
#include"EventPick.h"
#include<TFile.h>
#include<TTree.h>
#include<TDirectory.h>
#include<TObject.h>

int main(int ac, char** av){
	if(ac < 3){
		std::cout << "usage: ./makeSkim outputFileName inputFile[s]" << std::endl;
		return -1;
	}
	// input: dealing with TTree first
	EventTree* tree = new EventTree(ac-2, av+2);
	Selector* selector = new Selector();
	EventPick* evtPick = new EventPick("nominal");
	
	// no MET cut, for QCD fits
	evtPick->MET_cut = -1.0;

	// antiselection for QCD fit
	//selector->ele_Iso_MVA_invert = true;
	evtPick->NBjet_ge = 0;

	// add more branches to be saved
	//tree->chain->SetBranchStatus("*",1);

	TFile* outFile = new TFile( av[1] ,"RECREATE" );
	TDirectory* ggDir = outFile->mkdir("ggNtuplizer","ggNtuplizer");
	ggDir->cd();
	TTree* newTree = tree->chain->CloneTree(0);
	
	Long64_t nEntr = tree->GetEntries();
	for(Long64_t entry=0; entry<nEntr; entry++){
		if(entry%10000 == 0) std::cout << "processing entry " << entry << " out of " << nEntr << std::endl;
		//if(entry == 1000000) break;
		tree->GetEntry(entry);
		selector->process_objects(tree);
		evtPick->process_event(tree, selector);
		// make selection here
		if( evtPick->passPreSel )
			newTree->Fill();
	}
	newTree->Write();
	
	std::map<std::string, TH1F*> histMap;
	// copy histograms
	for(int fileInd = 2; fileInd < ac; ++fileInd){
		TFile* tempFile = new TFile(av[fileInd], "READ");
		TIter next(((TDirectory*)tempFile->Get("ggNtuplizer"))->GetListOfKeys());
		TObject* obj;
		while (obj = next()){
			std::string objName(obj->GetName());
			if( objName != "EventTree"){
				TH1F* hist = (TH1F*)tempFile->Get(("ggNtuplizer/"+objName).c_str());
				if( histMap.find(objName) != histMap.end() ){
					histMap[objName]->Add(hist);
				}
				else {
					hist->SetDirectory(0);
					histMap[objName] = hist;
				}
			}
		}
		tempFile->Close();
	}
	
	ggDir->cd();
	for(std::map<std::string, TH1F*>::iterator it = histMap.begin(); it!= histMap.end(); ++it){
		it->second->SetDirectory(ggDir);
		it->second->Write();
	}
	outFile->Close();
	
	return 0;
}
