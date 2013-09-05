#include<iostream>

#include"EventTree.h"
#include"Selector.h"
#include"Histogrammer.h"

//#include"JetMETObjects/JetCorrectorParameters.h"
//#include"JetMETObjects/FactorizedJetCorrector.h"

//void applyJEC(EventTree* tree, FactorizedJetCorrector* JetCorrector);
void doJER(EventTree* tree);
double JERcorrection(double eta);
bool overlapWHIZARD(EventTree* tree);

int main(int ac, char** av){
	if(ac < 4){
		std::cout << "usage: ./makeTemplates sampleName outputDir inputFile[s]" << std::endl;
		return -1;
	}

	// everything
	Histogrammer* histnom_ = new Histogrammer((std::string("1pho_top_")+av[1]).c_str());
	Histogrammer* histnom_barrel = new Histogrammer((std::string("1pho_barrel_top_")+av[1]).c_str());
	Histogrammer* histnom_endcap = new Histogrammer((std::string("1pho_endcap_top_")+av[1]).c_str());

	// MC only: photon is:  real signal / real background / fake from e / fake from jet 
	// photon is real signal
	Histogrammer* histnom_rs = new Histogrammer((std::string("1pho_rs_top_")+av[1]).c_str());
	Histogrammer* histnom_barrel_rs = new Histogrammer((std::string("1pho_rs_barrel_top_")+av[1]).c_str());
	Histogrammer* histnom_endcap_rs = new Histogrammer((std::string("1pho_rs_endcap_top_")+av[1]).c_str());

	// photon is real but comes from background
	Histogrammer* histnom_rb = new Histogrammer((std::string("1pho_rb_top_")+av[1]).c_str());
	Histogrammer* histnom_barrel_rb = new Histogrammer((std::string("1pho_rb_barrel_top_")+av[1]).c_str());
	Histogrammer* histnom_endcap_rb = new Histogrammer((std::string("1pho_rb_endcap_top_")+av[1]).c_str());

	// photon is fake from e
	Histogrammer* histnom_fe = new Histogrammer((std::string("1pho_fe_top_")+av[1]).c_str());
	Histogrammer* histnom_barrel_fe = new Histogrammer((std::string("1pho_fe_barrel_top_")+av[1]).c_str());
	Histogrammer* histnom_endcap_fe = new Histogrammer((std::string("1pho_fe_endcap_top_")+av[1]).c_str());

	// photon is fake from jet
	Histogrammer* histnom_fj = new Histogrammer((std::string("1pho_fj_top_")+av[1]).c_str());
	Histogrammer* histnom_barrel_fj = new Histogrammer((std::string("1pho_fj_barrel_top_")+av[1]).c_str());
	Histogrammer* histnom_endcap_fj = new Histogrammer((std::string("1pho_fj_endcap_top_")+av[1]).c_str());
	
	// combination of fj+rb
	Histogrammer* histnom_fjrb = new Histogrammer((std::string("1pho_fjrb_top_")+av[1]).c_str());
	Histogrammer* histnom_barrel_fjrb = new Histogrammer((std::string("1pho_fjrb_barrel_top_")+av[1]).c_str());
	Histogrammer* histnom_endcap_fjrb = new Histogrammer((std::string("1pho_fjrb_endcap_top_")+av[1]).c_str());


	Selector* selectornom_ = new Selector("nominal");

	//selectornom_->doVeto = false;
	selectornom_->pho_noSigmaIEta_cut = true;
	//selectornom_->pho_noChHadIso_cut = true;
	//selectornom_->pho_noPhoIso_cut = true;
	//selectornom_->pho_noIso_cut = true;
	selectornom_->pho_Et_cut = 25;
	selectornom_->pho_ID_ind = 0; //2 = tight ID
	//selectornom_->MET_cut = 0;
	//selectornom_->NBjet_cut = 0; // no b-jet cut
	//selectornom_->veto_jet_dR = 0.3;
	
	//selectornom_->no_trigger = true;
	//selectornom_->NlooseEleVeto_cut = 10; // no Loose Electron Veto
	
	//selectornom_->jet_Pt_cuts.clear();
        //selectornom_->jet_Pt_cuts.push_back(35.0);

	//Histogrammer* hist_ = new Histogrammer((std::string("0pho_")+av[1]).c_str());
	//Selector* selector_ = new Selector("0pho");
	//selector_->pho_Et_cut = 2500;//no photons
	//selector_->Npho_cut = 0;
	//selector_->jet_Pt_cuts.clear();
        //selector_->jet_Pt_cuts.push_back(30.0);

	//selector_->pho_noSigmaIEta_cut = true;
	//selector_->pho_noIso_cut = true;
	//if( std::string(av[1]).find("WJets") == 0) selector_->pho_PID_cut = true;
	//std::cout << "pho_PID_cut " << selector_->pho_PID_cut << std::endl;
	
	//Histogrammer* hist1_ = new Histogrammer((std::string("0pho_4jets30_")+av[1]).c_str());
	//Selector* selector1_ = new Selector("0pho_4jets30");
	//selector1_->pho_Et_cut = 2500;//no photoins
	//selector1_->Npho_cut = 0;
	//selector1_->Njet_cut = 4;
	//selector1_->jet_Pt_cuts.clear();
	//selector1_->jet_Pt_cuts.push_back(20.0);
	
	bool WHIZARD = false;
	if( std::string(av[1]).find("WHIZARD") != std::string::npos) WHIZARD = true;

	bool doOverlapRemoval = false;
	if( std::string(av[1]).find("TTJets") != std::string::npos || std::string(av[1]).find("TTgamma") != std::string::npos ) doOverlapRemoval = true;

	EventTree* tree_ = new EventTree(ac-3, av+3);
	double PUweight = 1.0;
	bool isMC = false;
	TH1D* PUweightHist;
	
	Long64_t nEntr = tree_->GetEntries();
	if( nEntr>0 ){
		tree_->GetEntry(0);
		if( !tree_->isData_ ){
			TFile* pileupFile = new TFile("Pileup_Observed_69300.root","READ");
			PUweightHist = (TH1D*)pileupFile->Get("pileup");
			PUweightHist->SetDirectory(0);
			pileupFile->Close();
			isMC = true;
			double PUweightInt = PUweightHist->Integral();
			TH1F* mcPU=NULL;
			for(int nmcfile = 3; nmcfile<ac; nmcfile++){
				std::cout << "reading file " << av[nmcfile] << std::endl;
				TFile* mcFile = new TFile(av[nmcfile],"READ");
				if( mcPU==NULL) mcPU = (TH1F*)mcFile->Get("ggNtuplizer/hPU");
				else mcPU->Add((TH1F*)mcFile->Get("ggNtuplizer/hPU"));
				mcPU->SetDirectory(0);
				mcFile->Close();
			}
			mcPU->Scale(1.0/mcPU->Integral());
			PUweightHist->Divide(mcPU);
			PUweightHist->Scale(1.0/PUweightInt);
		}
	}

	// JEC stuff ===============================================================================================
//	JetCorrectorParameters *ResJetPar; 
//	JetCorrectorParameters *L3JetPar;
//	JetCorrectorParameters *L2JetPar;
//	JetCorrectorParameters *L1JetPar;
//	if( isMC ){
//		ResJetPar = new JetCorrectorParameters("./jecAK5PF/Summer12_V7_MC_L2L3Residual_AK5PF.txt");
//		L3JetPar  = new JetCorrectorParameters("./jecAK5PF/Summer12_V7_MC_L3Absolute_AK5PF.txt");
//		L2JetPar  = new JetCorrectorParameters("./jecAK5PF/Summer12_V7_MC_L2Relative_AK5PF.txt");
//		L1JetPar  = new JetCorrectorParameters("./jecAK5PF/Summer12_V7_MC_L1FastJet_AK5PF.txt");
//	}
//	else{
//		ResJetPar = new JetCorrectorParameters("./jecAK5PF/Summer12_V7_DATA_L2L3Residual_AK5PF.txt");
//		L3JetPar  = new JetCorrectorParameters("./jecAK5PF/Summer12_V7_DATA_L3Absolute_AK5PF.txt");
//		L2JetPar  = new JetCorrectorParameters("./jecAK5PF/Summer12_V7_DATA_L2Relative_AK5PF.txt");
//		L1JetPar  = new JetCorrectorParameters("./jecAK5PF/Summer12_V7_DATA_L1FastJet_AK5PF.txt");
//	}
//	std::vector<JetCorrectorParameters> vPar;
//	vPar.push_back(*L1JetPar);
//	vPar.push_back(*L2JetPar);
//	vPar.push_back(*L3JetPar);
//	if (!isMC) vPar.push_back(*ResJetPar);
//	FactorizedJetCorrector *JetCorrector = new FactorizedJetCorrector(vPar);
	// JEC stuff ===============================================================================================

	double sumPUweight=0;
	for(Long64_t entry=0; entry<nEntr; entry++){
		if(entry%10000 == 0) std::cout << "processing entry " << entry << " out of " << nEntr << std::endl;
		//if(entry==30000) break;
		tree_->GetEntry(entry);
		
		// JER smearing (not used now)
		//if(isMC) doJER(tree_);
		
		// JEC (was done in the ntuples)
		//applyJEC(tree_, JetCorrector);

		if(isMC) {
			for(int puInd=0; puInd<tree_->nPUInfo_; ++puInd){
				if( tree_->puBX_[puInd] == 0 ){
					PUweight = PUweightHist->GetBinContent(PUweightHist->GetXaxis()->FindBin(tree_->nPU_[puInd]));
					break;
				}
			}
		}
		sumPUweight+=PUweight;
		
		selectornom_->process_objects(tree_);
		//selector_->process_objects(tree_);
		//selector1_->process_objects(tree_);
		if(selectornom_->accept_counts()) {
			int phoind = selectornom_->Photons[0];
			// fill match and fake histograms first
			//if( tree_->phoGenIndex_[phoind] < 0 ) histnomfake_->fill(selectornom_,tree_, PUweight);
			//else {
			//	if( abs(tree_->phoGenMomPID[phoind])>23 ) histnomfake_->fill(selectornom_,tree_, PUweight);
			//	else histnommatch_->fill(selectornom_,tree_, PUweight);
			//}
			// filling analysis histograms
			if( isMC && doOverlapRemoval && overlapWHIZARD(tree_)){
				// overlapping part, not needed
			} else {
				// fill the histograms
				histnom_->fill(selectornom_,tree_, PUweight);
				bool barrel = fabs(tree_->phoEta_[phoind]) < 1.5;
				if(barrel) histnom_barrel->fill(selectornom_,tree_, PUweight);
				else histnom_endcap->fill(selectornom_,tree_, PUweight);
				// the following applies to MC only
				if(!isMC) continue;

				// split by photon origin
				bool rs = false, rb = false, fe = false, fj = false;
				// real photon coming from ISR or quarks/W
				int mcPhotonInd = -1;
				for(int mcInd=0; mcInd<tree_->nMC_; ++mcInd)
					if(dR(tree_->mcEta[mcInd],tree_->mcPhi[mcInd],tree_->phoEta_[phoind],tree_->phoPhi_[phoind]) < 0.2 && 
					fabs(tree_->phoEt_[phoind] - tree_->mcPt[mcInd]) / tree_->mcPt[mcInd] < 1.0 &&
					tree_->mcPID[mcInd] == 22){
						mcPhotonInd = mcInd;
						break;
					}

				if(mcPhotonInd >= 0){
					// signal: parents are quarks gluons or bosons
					if(tree_->mcParentage[mcPhotonInd]==2 || (tree_->mcParentage[mcPhotonInd]==10 /*&& abs(tree_->mcMomPID[mcPhotonInd])==24*/)) rs = true;
					else rb = true;
				}
				else{
					// no good matched Gen Photon found - our photon is fake
					//std::cout << "photon " << tree_->phoEt_[phoind] << " " << tree_->phoEta_[phoind] << " " << tree_->phoPhi_[phoind] << std::endl;
					for(int mcInd=0; mcInd<tree_->nMC_; ++mcInd){
						if(dR(tree_->mcEta[mcInd],tree_->mcPhi[mcInd],tree_->phoEta_[phoind],tree_->phoPhi_[phoind]) < 0.2){
							if(abs(tree_->mcPID[mcInd]) == 11 && fabs(tree_->phoEt_[phoind] - tree_->mcPt[mcInd]) / tree_->mcPt[mcInd] < 1.0) fe = true;
							if(fe) break;
							//std::cout << "ind " << mcInd << " PID " << tree_->mcPID[mcInd] << " genInd " << tree_->phoGenMomPID[phoind] << std::endl;
							//std::cout << "pt eta phi " << tree_->mcPt[mcInd] << " " << tree_->mcEta[mcInd] << " " << tree_->mcPhi[mcInd] << std::endl;
						}
					}
					if(!fe) fj = true; // call anything else a jet
				}
				// filling
				if(rs){
					histnom_rs->fill(selectornom_,tree_, PUweight);
					if(barrel) histnom_barrel_rs->fill(selectornom_,tree_, PUweight);
					else histnom_endcap_rs->fill(selectornom_,tree_, PUweight);
				}
				if(rb){
					histnom_rb->fill(selectornom_,tree_, PUweight);
					if(barrel) histnom_barrel_rb->fill(selectornom_,tree_, PUweight);
					else histnom_endcap_rb->fill(selectornom_,tree_, PUweight);
				}
				if(fe){
					histnom_fe->fill(selectornom_,tree_, PUweight);
					if(barrel) histnom_barrel_fe->fill(selectornom_,tree_, PUweight);
					else histnom_endcap_fe->fill(selectornom_,tree_, PUweight);
				}
				if(fj){
					histnom_fj->fill(selectornom_,tree_, PUweight);
					if(barrel) histnom_barrel_fj->fill(selectornom_,tree_, PUweight);
					else histnom_endcap_fj->fill(selectornom_,tree_, PUweight);
				}
				if(fj||rb){
					histnom_fjrb->fill(selectornom_,tree_, PUweight);
					if(barrel) histnom_barrel_fjrb->fill(selectornom_,tree_, PUweight);
					else histnom_endcap_fjrb->fill(selectornom_,tree_, PUweight);
				}



			}

		}
		//if(selector_->accept_counts()) hist_->fill(selector_,tree_, PUweight);
		//if(selector1_->accept_counts()) hist1_->fill(selector1_,tree_, PUweight);
	}
	sumPUweight/=nEntr;
	std::cout << "Average PU weight " << sumPUweight << std::endl;
	selectornom_->print_cutflow();
	//selector_->print_cutflow();
	//selector1_->print_cutflow();
	//hist_->write_histograms(av[2],selector_->histVector);
	std::vector<TH1F*> emptyVec;
	bool fillMCAll = true;	
	bool fillBarrel = true;
	bool fillEndcap = true;
	bool fillFE = true;
	bool fillFJ = true;
	bool fillRB = true;
	bool fillFJRB = true;

	histnom_->write_histograms(av[2],selectornom_->histVector);
	if(fillBarrel) histnom_barrel->write_histograms(av[2],emptyVec);
	if(fillEndcap) histnom_endcap->write_histograms(av[2],emptyVec);

	if(isMC){
		histnom_rs->write_histograms(av[2],selectornom_->histVector);
		if(fillBarrel) histnom_barrel_rs->write_histograms(av[2],emptyVec);
		if(fillEndcap) histnom_endcap_rs->write_histograms(av[2],emptyVec);
		
		if(fillRB){
			if(fillMCAll) histnom_rb->write_histograms(av[2],selectornom_->histVector);
			if(fillBarrel) histnom_barrel_rb->write_histograms(av[2],emptyVec);
			if(fillEndcap) histnom_endcap_rb->write_histograms(av[2],emptyVec);
		}
		if(fillFE){
			if(fillMCAll) histnom_fe->write_histograms(av[2],selectornom_->histVector);
			if(fillBarrel) histnom_barrel_fe->write_histograms(av[2],emptyVec);
			if(fillEndcap) histnom_endcap_fe->write_histograms(av[2],emptyVec);
		}
		if(fillFJ){
			if(fillMCAll) histnom_fj->write_histograms(av[2],selectornom_->histVector);
			if(fillBarrel) histnom_barrel_fj->write_histograms(av[2],emptyVec);
			if(fillEndcap) histnom_endcap_fj->write_histograms(av[2],emptyVec);
		}
		if(fillFJRB){
			if(fillMCAll) histnom_fjrb->write_histograms(av[2],selectornom_->histVector);
			if(fillBarrel) histnom_barrel_fjrb->write_histograms(av[2],emptyVec);
			if(fillEndcap) histnom_endcap_fjrb->write_histograms(av[2],emptyVec);
		}
	}
	//histnomfake_->write_histograms(av[2],emptyVec);
	//histnommatch_->write_histograms(av[2],emptyVec);
	//hist_->write_histograms(av[2],selector_->histVector);
	//hist1_->write_histograms(av[2],selector1_->histVector);
	
	delete tree_;
	return 0;
}

//void applyJEC(EventTree* tree, FactorizedJetCorrector* JetCorrector){
//	for(int jetInd = 0; jetInd < tree->nJet_ ; ++jetInd){
//		JetCorrector->setJetEta(tree->jetEta_[jetInd]);
//		JetCorrector->setJetPt(tree->jetRawPt_[jetInd]);
//		JetCorrector->setJetA(tree->jetArea_[jetInd]);
//		JetCorrector->setRho(tree->rho2012_); 
//		
//		double correction = JetCorrector->getCorrection();
//		//std::cout << tree->jetPt_[jetInd] << "  " << tree->jetRawPt_[jetInd] << "  " << tree->jetPt_[jetInd]/tree->jetRawPt_[jetInd] << "  " << correction << std::endl;
//		tree->jetPt_[jetInd] = tree->jetRawPt_[jetInd] * correction;
//		tree->jetEn_[jetInd] = tree->jetRawEn_[jetInd] * correction;
//	}
//}

bool overlapWHIZARD(EventTree* tree){
	const double Et_cut = 20;
	const double dR_cut = 0.1;
	// consider only mcParticles here
	bool haveOverlap = false;
	for(int phoInd=0; phoInd<tree->nMC_; ++phoInd){
		// find photons with Et above cut and that came from ISR or top
		if(tree->mcPID[phoInd]==22 && tree->mcPt[phoInd] > Et_cut && (tree->mcParentage[phoInd]==2 || (tree->mcParentage[phoInd]==10 && abs(tree->mcMomPID[phoInd])==24))){
			bool closeToLeg = false;
			bool haveLeg = false;
			// find "legs" (b-qark coming from top)
			for(int legInd=0; legInd<tree->nMC_; ++legInd){
				if(abs(tree->mcPID[legInd])==5 && abs(tree->mcMomPID[legInd])==6){
					haveLeg=true;
					if(dR(tree->mcEta[phoInd],tree->mcPhi[phoInd],tree->mcEta[legInd],tree->mcPhi[legInd]) < dR_cut) closeToLeg=true;
				}
			}
			if(haveLeg && !closeToLeg) haveOverlap = true;
		}
	}
	if (0 && !haveOverlap){
		std::cout << "no overlap" << std::endl;
		for(int legInd=0; legInd<tree->nMC_; ++legInd)
			std::cout << "PID " << tree->mcPID[legInd] 
			<< "  Pt " << tree->mcPt[legInd] 
			<< "  Eta " << tree->mcEta[legInd] 
			<< "  Phi " << tree->mcPhi[legInd] << std::endl;
	}
	return haveOverlap;
}

// https://twiki.cern.ch/twiki/bin/view/CMS/JetResolution
void doJER(EventTree* tree){
	// scale
	for(int jetInd = 0; jetInd < tree->nJet_ ; ++jetInd){
		if(tree->jetGenJetIndex_[jetInd]>0){
			double oldPt = tree->jetPt_[jetInd];
			double genPt = tree->jetGenJetPt_[jetInd];
			double eta = tree->jetEta_[jetInd];
			tree->jetPt_[jetInd] = std::max(0.0, genPt + JERcorrection(eta)*(oldPt-genPt));
			//std::cout << "old " << oldPt << "  new " << tree->jetPt_[jetInd] << std::endl;
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
