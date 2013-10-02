#include<iostream>

#include"EventTree.h"
#include"Selector.h"
#include"EventPick.h"
#include"Histogrammer.h"
#include"HistCollect.h"
#include"PUReweight.h"

void doJER(EventTree* tree);
double JERcorrection(double eta);
bool overlapWHIZARD(EventTree* tree);

int main(int ac, char** av){
	if(ac < 4){
		std::cout << "usage: ./makeTemplates sampleName outputDir inputFile[s]" << std::endl;
		return -1;
	}

	// book HistCollect
	HistCollect* looseCollect = new HistCollect("1pho",std::string("top_")+av[1]);
	// HistCollect for tight Photon ID
	//HistCollect* tightCollect = new HistCollect("1photight",std::string("top_")+av[1]);
	
	// object selectors
	Selector* selectorLoose = new Selector();
	//Selector* selectorTight = new Selector();
	// set up the parameters for object selectors here
	//selectorTight->pho_ID_ind = 2; // tight ID
	
	// create event selectors here
	EventPick* evtPickLoose = new EventPick("LoosePhotonID");
	//EventPick* evtPickTight = new EventPick("TightPhotonID");
	
	
	bool WHIZARD = false;
	if( std::string(av[1]).find("WHIZARD") != std::string::npos) WHIZARD = true;

	bool doOverlapRemoval = false;
	if( std::string(av[1]).find("TTJets") != std::string::npos || std::string(av[1]).find("TTgamma") != std::string::npos ) doOverlapRemoval = true;

	EventTree* tree = new EventTree(ac-3, av+3);
	double PUweight = 1.0;
	bool isMC = false;
	
	// initialize PU reweighting here
	PUReweight* PUweighter = new PUReweight(ac-3, av+3);
	
	Long64_t nEntr = tree->GetEntries();
	for(Long64_t entry=0; entry<nEntr; entry++){
		if(entry%10000 == 0) std::cout << "processing entry " << entry << " out of " << nEntr << std::endl;
		//if(entry==30000) break;
		tree->GetEntry(entry);
		isMC = !(tree->isData_);
		
		// JER smearing (not used now)
		//if(isMC) doJER(tree_);
		
		// apply PU reweighting
		PUweight = PUweighter->getWeight(tree->nPUInfo_, tree->puBX_, tree->nPU_);
		
		// do overlap removal here
		if( isMC && doOverlapRemoval && overlapWHIZARD(tree)){
			// overlapping part, not needed
			continue;
		}
		
		selectorLoose->process_objects(tree);
		//selectorTight->process_objects(tree);
		
		evtPickLoose->process_event(tree, selectorLoose, PUweight);
		//evtPickTight->process_event(tree, selectorTight, PUweight);
		
		// fill the histograms
		looseCollect->fill_histograms(selectorLoose, evtPickLoose, tree, isMC, PUweight);
	}
	
	looseCollect->write_histograms(evtPickLoose, isMC, av[2]);
	
	std::cout << "Average PU weight " << PUweighter->getAvgWeight() << std::endl;
	evtPickLoose->print_cutflow();
	
	delete tree;
	return 0;
}


// https://twiki.cern.ch/twiki/bin/view/CMS/JetResolution
void doJER(EventTree* tree){
	// scale
	for(int jetInd = 0; jetInd < tree->nJet_ ; ++jetInd){
		if(tree->jetGenJetIndex_->at(jetInd)>0){
			double oldPt = tree->jetPt_->at(jetInd);
			double genPt = tree->jetGenJetPt_->at(jetInd);
			double eta = tree->jetEta_->at(jetInd);
			tree->jetPt_->at(jetInd) = std::max(0.0, genPt + JERcorrection(eta)*(oldPt-genPt));
			//std::cout << "old " << oldPt << "  new " << tree->jetPt_->at(jetInd) << std::endl;
		}
	}
	// reorder if needed
	//for(int jetInd = 0; jetInd < tree->nJet_ ; ++jetInd){
	//	
	//}
}


double JERcorrection(double JetEta){
	double eta = TMath::Abs(JetEta);
	static const double corr[5] = {1.052, 1.057, 1.096, 1.134, 1.288};
	int region = 0;
	if( eta >= 0.5 ) region++;
	if( eta >= 1.1 ) region++;
	if( eta >= 1.7 ) region++;
	if( eta >= 2.3 ) region++;
	return corr[region];
}
