#include<iostream>

#include"EventTree.h"
#include"Selector.h"
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
	EventTree* tree_ = new EventTree(ac-2, av+2);
	Selector* selector_ = new Selector("nominal");
	
	//selector_->jet_Pt_cuts.clear();
	//selector_->jet_Pt_cuts.push_back(30.0);
	
	selector_->ele_Pt_cut = 30.0;
	selector_->pho_Et_cut = 15.0;
	selector_->pho_noSigmaIEta_cut = true;
	selector_->pho_noChHadIso_cut = true;
	selector_->veto_jet_dR = 0.3;
	selector_->veto_pho_jet_dR = 0.5;
        selector_->veto_pho_lep_dR = 0.5;

	// add more branches to be saved
	tree_->chain->SetBranchStatus("*",1);

	//tree_->chain->SetBranchStatus("pho*",1);
	//tree_->chain->SetBranchStatus("phoCiC*",0);	
	// create output file
	TFile* outFile = new TFile( av[1] ,"RECREATE" );
	TDirectory* ggDir = outFile->mkdir("ggNtuplizer","ggNtuplizer");
	ggDir->cd();
	TTree* newTree = tree_->chain->CloneTree(0);
	
	Long64_t nEntr = tree_->GetEntries();
	for(Long64_t entry=0; entry<nEntr; entry++){
		if(entry%10000 == 0) std::cout << "processing entry " << entry << " out of " << nEntr << std::endl;
		//if(entry == 1000000) break;
		tree_->GetEntry(entry);
		selector_->process_objects(tree_);
		// make selection here
		if( selector_->Jets.size() >= 3 &&
		   selector_->Electrons.size() >= 1 &&
		   //selector_->MuonsLoose.size() == 0 &&
		   selector_->Photons.size() >= 1 )
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
	delete selector_;
	delete tree_;
	
	return 0;
}
