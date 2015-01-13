#include<iostream>

#include"EventTree.h"
#include"Selector.h"
#include"EventPick.h"

void doJER(EventTree* tree);
double JERcorrection(double eta);

void fillCategory(EventTree* tree, TH1F* hist, double weight){
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
	hist->Fill(1.0, weight); // Total
	int nEle = EleP + EleM;
	int nMu = MuP + MuM;
	int nTau = TauP + TauM;
	if( nEle + nMu + nTau == 0) hist->Fill(2.0, weight); // All Had
	if( nEle + nMu + nTau == 1) hist->Fill(3.0, weight); // Single Lepton
	if( nEle + nMu + nTau == 2) hist->Fill(4.0, weight); // Di Lepton
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
		
		doJER(tree);

		selectorLoose->process_objects(tree);
		evtPickLoose->process_event(tree, selectorLoose, PUweight);

		// fill the histograms
		fillCategory(tree, allCategory, PUweight);
		if(evtPickLoose->passPreSel) fillCategory(tree, preselCategory, PUweight);
		if(evtPickLoose->passAll) fillCategory(tree, photonCategory, PUweight);
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


// https://twiki.cern.ch/twiki/bin/view/CMS/JetResolution
void doJER(EventTree* tree){
	// correct MET when jets are smeared
	TLorentzVector tMET;
	tMET.SetPtEtaPhiM(tree->pfMET_,0.0,tree->pfMETPhi_,0.0);
	//std::cout << "before correction MET " << tree->pfMET_ << " " << tree->pfMETPhi_ << "    ";
	// scale jets
	for(int jetInd = 0; jetInd < tree->nJet_ ; ++jetInd){
		if(tree->jetPt_->at(jetInd) < 20) continue;
		if(tree->jetGenJetIndex_->at(jetInd)>0){
			TLorentzVector tjet;
			tjet.SetPtEtaPhiM(tree->jetPt_->at(jetInd), tree->jetEta_->at(jetInd), tree->jetPhi_->at(jetInd), 0.0);
			tMET+=tjet;
			double oldPt = tree->jetPt_->at(jetInd);
			double genPt = tree->jetGenJetPt_->at(jetInd);
			double eta = tree->jetEta_->at(jetInd);
			tree->jetPt_->at(jetInd) = std::max(0.0, genPt + JERcorrection(eta)*(oldPt-genPt));
			tjet.SetPtEtaPhiM(tree->jetPt_->at(jetInd), tree->jetEta_->at(jetInd), tree->jetPhi_->at(jetInd), 0.0);			
			tMET-=tjet;
			//std::cout << "old " << oldPt << "  new " << tree->jetPt_->at(jetInd) << std::endl;
		}
	}
	// save updated MET values
	tree->pfMET_ = tMET.Pt();
	tree->pfMETPhi_ = tMET.Phi();
	//std::cout << "after corrections " << tree->pfMET_ << " " << tree->pfMETPhi_ << std::endl;
}


double JERcorrection(double JetEta){
	double eta = TMath::Abs(JetEta);
	static const double corr[5] = {1.052, 1.057, 1.096, 1.134, 1.288};
	static const double corrDown[5] = {0.990, 1.001, 1.032, 1.042, 1.089};
	static const double corrUp[5] = {1.115, 1.114, 1.161, 1.228, 1.488};
	int region = 0;
	if( eta >= 0.5 ) region++;
	if( eta >= 1.1 ) region++;
	if( eta >= 1.7 ) region++;
	if( eta >= 2.3 ) region++;
	return corr[region];
	
	// should not get here
	return 1.0;
}

