#include<iostream>

#include"EventTree.h"
#include"Selector.h"
#include"EventPick.h"
#include"Histogrammer.h"
#include"HistCollect.h"
#include"PUReweight.h"

#include"TRandom3.h"
// temporary solution
#include"JECvariation.cpp"

int jecvar012_g = 1; // 0:down, 1:norm, 2:up
int jervar012_g = 1; // 0:down, 1:norm, 2:up
int eleeff012_g = 1; // 0:down, 1:norm, 2:up
int btagvar012_g = 1; // 0:down, 1:norm, 2:up
int phosmear01_g = 0; // 1 : do it, 0 : don't do it
int elesmear01_g = 0; // 1 : do it, 0 : don't do it

double getEleEff(EventTree* tree, EventPick* evt);
double getBtagSF(EventTree* tree, EventPick* evt);
void doEleSmearing(EventTree* tree);
void doPhoSmearing(EventTree* tree);
void doJER(EventTree* tree);
double JERcorrection(double eta);
bool overlapWHIZARD(EventTree* tree);

int main(int ac, char** av){
	if(ac < 4){
		std::cout << "usage: ./makeTemplates sampleName outputDir inputFile[s]" << std::endl;
		return -1;
	}

	std::cout << "JEC: " << jecvar012_g << "  JER: " << jervar012_g << "  EleEff: " << eleeff012_g << "  BtagVar: " << btagvar012_g << "  ";
	std::cout << "  PhoSmear: " << phosmear01_g << "  EleSmear: " << elesmear01_g << std::endl;
	// book HistCollect
	HistCollect* looseCollect = new HistCollect("1pho",std::string("top_")+av[1]);
	HistCollect* looseCollectNoMET = new HistCollect("1phoNoMET",std::string("top_")+av[1]);
	//HistCollect* fourjCollect = new HistCollect("1pho4j",std::string("top4j_")+av[1]);
	// HistCollect for tight Photon ID
	//HistCollect* tightCollect = new HistCollect("1photight",std::string("top_")+av[1]);
	
	// object selectors
	Selector* selectorLoose = new Selector();
	bool isQCD = false;
	if(std::string(av[3]).find("QCD") != std::string::npos){
		isQCD = true;
		selectorLoose->ele_MVA_range[0] = -1.0;
		selectorLoose->ele_MVA_range[1] = -0.1;
		selectorLoose->ele_RelIso_range[0] = 0.25;
		selectorLoose->ele_RelIso_range[1] = 1.0;
		//selectorLoose->ele_Iso_MVA_invert = true;

		// no MC categories for data-driven QCD needed
		looseCollect->fillRS = false;
		looseCollect->fillFE = false;
		looseCollect->fillFJRB = false;

		looseCollectNoMET->fillRS = false;
		looseCollectNoMET->fillFE = false;
		looseCollectNoMET->fillFJRB = false;
	}
	//selectorLoose->jet_Pt_cut = 15000;
	//selectorLoose->ele_Pt_cut = 25;
	//selectorLoose->pho_Et_cut = 20.0;
	//selectorLoose->pho_ID_ind = 2;  //pho ID	

	//Selector* selectorTight = new Selector();
	// set up the parameters for object selectors here
	//selectorTight->pho_ID_ind = 2; // tight ID
	
	// create event selectors here
	EventPick* evtPickLoose = new EventPick("LoosePhotonID");
	//evtPickLoose->veto_pho_jet_dR = 0.05;
	//evtPickLoose->Njet_ge = 4;
	EventPick* evtPickLooseNoMET = new EventPick("LoosePhotonID");
	evtPickLooseNoMET->MET_cut = -1.0;
	//evtPickLooseNoMET->veto_pho_jet_dR = 0.05;
	//evtPickLooseNoMET->Njet_ge = 4;
	//evtPickLoose->no_trigger = true;
	//evtPickLoose->NBjet_ge = 0;
	//evtPickLoose->Nele_eq = 1;
	//evtPickLoose->NlooseMuVeto_le = 99;
	
	//evtPickLoose->veto_jet_dR = 0.1;
	//evtPickLoose->veto_pho_jet_dR = 0.8;
	//evtPickLoose->veto_pho_lep_dR = 0.8;
	//EventPick* evtPickLoose4j = new EventPick("LoosePhotonID4j");
	//evtPickLoose4j->Njet_ge = 4;
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
	
	tree->GetEntry(0);
	isMC = !(tree->isData_);
	JECvariation* jecvar;
	jecvar = new JECvariation("./jecAK5PF/Summer12_V7", isMC);
	
	Long64_t nEntr = tree->GetEntries();
	for(Long64_t entry=0; entry<nEntr; entry++){
		if(entry%10000 == 0) std::cout << "processing entry " << entry << " out of " << nEntr << std::endl;
		//if(entry==5000000) break;
		tree->GetEntry(entry);
		isMC = !(tree->isData_);
		
		// apply PU reweighting
		if(isMC) PUweight = PUweighter->getWeight(tree->nPUInfo_, tree->puBX_, tree->nPU_);
		
		if(isMC && !isQCD){
			// JEC
			jecvar->applyJEC(tree, jecvar012_g); // 0:down, 1:norm, 2:up
			// JER smearing 
			doJER(tree);
			// photon energy smearing
			doPhoSmearing(tree);
			// electron energy smearing
			doEleSmearing(tree);
		}
		// do overlap removal here
		if( isMC && doOverlapRemoval && overlapWHIZARD(tree)){
			// overlapping part, not needed
			continue;
		}
		
		selectorLoose->process_objects(tree);
		//selectorTight->process_objects(tree);
		

		evtPickLoose->process_event(tree, selectorLoose, PUweight);
		evtPickLooseNoMET->process_event(tree, selectorLoose, PUweight);

		//evtPickLoose4j->process_event(tree, selectorLoose, PUweight);
		//evtPickTight->process_event(tree, selectorTight, PUweight);
		
		double evtWeight = PUweight;
		if(isMC && !isQCD){
			// electron trigger efficiency reweighting
			evtWeight *= getEleEff(tree, evtPickLoose);
			// b-tag SF reweighting
			evtWeight *= getBtagSF(tree, evtPickLoose);
		}
		// fill the histograms
		looseCollect->fill_histograms(selectorLoose, evtPickLoose, tree, isMC, evtWeight);
		looseCollectNoMET->fill_histograms(selectorLoose, evtPickLooseNoMET, tree, isMC, evtWeight);
		//fourjCollect->fill_histograms(selectorLoose, evtPickLoose4j, tree, isMC, evtWweight);
	}
	
	looseCollect->write_histograms(evtPickLoose, isMC, av[2]);
	looseCollectNoMET->write_histograms(evtPickLoose, isMC, av[2]);
	//fourjCollect->write_histograms(evtPickLoose4j, isMC, av[2]);

	std::cout << "Average PU weight " << PUweighter->getAvgWeight() << std::endl;
	evtPickLoose->print_cutflow();
	
	delete tree;
	return 0;
}

// AN2012_438_v10 page9
double getEleEff(EventTree* tree, EventPick* evt){
	static double trigEffSF_PtEta[3][3] = { {0.984, 0.967, 0.991}, {0.999, 0.983, 1.018}, {0.999, 0.988, 0.977} };
	static double trigEffSFerr_PtEta[3][3] = { {0.002, 0.002, 0.007}, {0.003, 0.002, 0.012}, {0.002, 0.003, 0.015} };
	static double idEffSF_PtEta[3][3] = { {0.950, 0.957, 0.922}, {0.966, 0.961, 0.941}, {0.961, 0.963, 0.971} };
	static double idEffSFerr_PtEta[3][3] = { {0.003, 0.002, 0.004}, {0.001, 0.002, 0.007}, {0.002, 0.003, 0.0} };
	
	if( evt->Electrons.size() < 1 ) return 1.0; // no electrons, no weight
	int eleInd = evt->Electrons[0];
	double pt = tree->elePt_->at(eleInd);
	double eta = tree->eleSCEta_->at(eleInd);
	int etaRegion = 0;
	if( eta > 0.80) etaRegion++;
	if( eta > 1.48) etaRegion++;

	int ptRegion = 0;
	if( pt > 40 ) ptRegion++;
	if( pt > 50 ) ptRegion++;

	if(eleeff012_g == 1) return trigEffSF_PtEta[ptRegion][etaRegion] * idEffSF_PtEta[ptRegion][etaRegion];
	if(eleeff012_g == 0) return (trigEffSF_PtEta[ptRegion][etaRegion] - trigEffSFerr_PtEta[ptRegion][etaRegion]) * (idEffSF_PtEta[ptRegion][etaRegion] - idEffSFerr_PtEta[ptRegion][etaRegion]);
	if(eleeff012_g == 2) return (trigEffSF_PtEta[ptRegion][etaRegion] + trigEffSFerr_PtEta[ptRegion][etaRegion]) * (idEffSF_PtEta[ptRegion][etaRegion] + idEffSFerr_PtEta[ptRegion][etaRegion]);	
}


void doEleSmearing(EventTree* tree){
	static TRandom3 rand;
	if(elesmear01_g == 0) return;
	for(int eleInd = 0; eleInd < tree->nEle_; ++eleInd){
		if(tree->elePt_->at(eleInd) < 15) continue;
		//std::cout << "electron Pt before " << tree->elePt_->at(eleInd) << "   "; 
		double factor = rand.Gaus(1.0,0.01);
		//std::cout << "factor " << factor << "  ";
		tree->elePt_->at(eleInd) *= factor;
		//std::cout << "electron Pt after " << tree->elePt_->at(eleInd) << std::endl;
	}
}

void doPhoSmearing(EventTree* tree){
	static TRandom3 rand;
	if(phosmear01_g == 0) return;
	for(int phoInd = 0; phoInd < tree->nPho_; ++phoInd){
		if(tree->phoEt_->at(phoInd) < 15) continue;
		//std::cout << "photon Et before " << tree->phoEt_->at(phoInd) << "   "; 
		double factor = rand.Gaus(1.0,0.01);
		//std::cout << "factor " << factor << "  ";
		tree->phoEt_->at(phoInd) *= factor;
		//std::cout << "photon Et after " << tree->phoEt_->at(phoInd) << std::endl;
	}
}

// https://twiki.cern.ch/twiki/bin/view/CMS/JetResolution
void doJER(EventTree* tree){
	// scale
	for(int jetInd = 0; jetInd < tree->nJet_ ; ++jetInd){
		if(tree->jetPt_->at(jetInd) < 20) continue;
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
	static const double corrDown[5] = {0.990, 1.001, 1.032, 1.042, 1.089};
	static const double corrUp[5] = {1.115, 1.114, 1.161, 1.228, 1.488};
	int region = 0;
	if( eta >= 0.5 ) region++;
	if( eta >= 1.1 ) region++;
	if( eta >= 1.7 ) region++;
	if( eta >= 2.3 ) region++;
	if(jervar012_g == 0) return corrDown[region];
	if(jervar012_g == 1) return corr[region];
	if(jervar012_g == 2) return corrUp[region];
}

// https://twiki.cern.ch/twiki/bin/viewauth/CMS/BTagSFMethods
// https://twiki.cern.ch/twiki/pub/CMS/BtagPOG/SFb-pt_NOttbar_payload_EPS13.txt
// weight for >=1 btag :  1 - prod(1-SFi) over all b-tagged jets

double getBtagSF(EventTree* tree, EventPick* evt){
	//Tagger: CSVM within 20 < pt < 800 GeV, abs(eta) < 2.4, x = pt
	static double ptmax[16] = {30, 40, 50, 60, 70, 80,100, 120, 160, 210, 260, 320, 400, 500, 600, 800};	
	static double SFb_error[16] = {0.0415694, 0.023429, 0.0261074, 0.0239251, 0.0232416, 0.0197251, 0.0217319, 0.0198108, 0.0193, 0.0276144, 0.020583, 0.026915, 0.0312739, 0.0415054, 0.074056, 0.0598311};
	
	double prod = 1.0;
	double jetpt;
	int ptInd=0;
	double SFb;
	for(std::vector<int>::const_iterator bjetInd = evt->bJets.begin(); bjetInd != evt->bJets.end(); bjetInd++){
		jetpt = tree->jetPt_->at(*bjetInd);
		SFb = (0.939158+(0.000158694*jetpt))+(-2.53962e-07*(jetpt*jetpt));
		ptInd = 0;
		if( btagvar012_g != 1){
			for(int ipt=0; ipt<15; ipt++) if(jetpt > ptmax[ipt]) ptInd++;

			if( btagvar012_g == 0 ) SFb -= SFb_error[ptInd];
			if( btagvar012_g == 2 ) SFb += SFb_error[ptInd];
		}
		prod *= 1.0 - SFb;
	}
	return 1.0 - prod;
}

