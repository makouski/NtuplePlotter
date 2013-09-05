#include"Histogrammer.h"

Histogrammer::Histogrammer(const char* titleIn){
	title.assign(titleIn);

	// 2d histograms
	hists2d["photon1_Sigma_ChIso"] = new TH2F("photon1_Sigma_ChIso","photon1 SigmaIetaIeta vs ChIso",160,0,0.04,300,-10,20);
	hists2d["photon1_Sigma_ChIso"]->SetDirectory(0);

	hists2d["photon1_Sigma_ChSCRIso"] = new TH2F("photon1_Sigma_ChSCRIso","photon1 SigmaIetaIeta vs ChSCRIso",160,0,0.04,300,-10,20);
	hists2d["photon1_Sigma_ChSCRIso"]->SetDirectory(0);

	hists2d["photon1_Sigma_PhoIso"] = new TH2F("photon1_Sigma_PhoIso","photon1 SigmaIetaIeta vs PhoIso",160,0,0.04,300,-10,20);
	hists2d["photon1_Sigma_PhoIso"]->SetDirectory(0);

	hists2d["photon1_Sigma_PhoSCRIso"] = new TH2F("photon1_Sigma_PhoSCRIso","photon1 SigmaIetaIeta vs PhoSCRIso",160,0,0.04,300,-10,20);
	hists2d["photon1_Sigma_PhoSCRIso"]->SetDirectory(0);

	hists2d["photon1_25_35_Sigma_ChIso"] = new TH2F("photon1_25_35_Sigma_ChIso","photon1 Et 25 to 35 SigmaIetaIeta vs ChIso",160,0,0.04,30,-10,20);
	hists2d["photon1_25_35_Sigma_ChIso"]->SetDirectory(0);
	//hists2d["photon1_25_35_Sigma_ChIso"]->Sumw2();

	hists2d["photon1_35_45_Sigma_ChIso"] = new TH2F("photon1_35_45_Sigma_ChIso","photon1 Et 35 to 45 SigmaIetaIeta vs ChIso",160,0,0.04,30,-10,20);
	hists2d["photon1_35_45_Sigma_ChIso"]->SetDirectory(0);
	//hists2d["photon1_35_45_Sigma_ChIso"]->Sumw2();

	hists2d["photon1_45_60_Sigma_ChIso"] = new TH2F("photon1_45_60_Sigma_ChIso","photon1 Et 45 to 60 SigmaIetaIeta vs ChIso",160,0,0.04,30,-10,20);
	hists2d["photon1_45_60_Sigma_ChIso"]->SetDirectory(0);

	hists2d["photon1_60_up_Sigma_ChIso"] = new TH2F("photon1_60_up_Sigma_ChIso","photon1 Et 60 up SigmaIetaIeta vs ChIso",160,0,0.04,30,-10,20);
	hists2d["photon1_60_up_Sigma_ChIso"]->SetDirectory(0);

	// creating histograms
	make_hist("ele1Pt","electron 1 Pt",30,0,300,"p_{T} (GeV)","Events / 10 GeV");
	make_hist("ele1Eta","electron 1 Eta",26,-2.6,2.6,"#eta","Events / 0.2");
	make_hist("ele1RelIso","electron 1 relative isolation",12,0,0.12,"RelIso","Events / 0.01");
	make_hist("ele1MVA","electron 1 MVA",80,0.5,1.3,"MVA Trig","Events / 0.01");
	make_hist("ele1D0","electron 1 Dxy_PV",80,-0.2,0.2,"Dxy_PV (cm)","Events / 0.005 cm");
	make_hist("ele1Dz","electron 1 Dz",75,-0.15,0.15,"D_{Z} (cm)","Events / 0.004 cm");
	make_hist("ele1EoverP","electron 1 EoverP",25,0,5,"E/P","Events / 0.2");
	make_hist("ele1sigmaIetaIeta","electron 1 sigmaIetaIeta",80,0,0.04,"#sigma_{i#etai#eta}","Events / 0.0005");
	make_hist("ele1MissHits","electron 1 missing hits",10,-0.5,9.5,"missing hits","Events");
	make_hist("ele1DrJet","dR electron 1 to closest jet",60,0,6,"#DeltaR(e,jet)","Events / 0.1");
	make_hist("ele1MotherID","electron 1 mother PDG ID",35,-0.5,34.5,"e mother PID","Events");
	make_hist("ele1GMotherID","electron 1 Gmother PDG ID",35,-0.5,34.5,"e Gmother PID","Events");
	
	make_hist("WtransMass","W transverse mass",20,0,200,"M(W_{T})(GeV)","Events / 10 GeV");
	make_hist("WtransMass_photon","W transverse mass with photon",20,0,200,"M_{T}(#gamma,#nu)(GeV)","Events / 10 GeV");
	make_hist("ele1pho1Mass","electron + photon mass",20,0,200,"M(e,#gamma)(GeV)","Events / 10 GeV");	
	
	make_hist("ele2Pt","electron 2 Pt",30,0,300,"p_{T} (GeV)","Events / 10 GeV");
	make_hist("ele2RelIso","electron 2 relative isolation",10,0,0.5,"RelIso","Events / 0.05");
	
	make_hist("looseEleDrGenPho","dR loose electron to Gen Photon",60,0,6,"#DeltaR(e_{loose},#gamma_{MC})","Events / 0.1");
	
	make_hist("photon1Et","photon 1 Et",20,0,200,"E_{T} (GeV)","Events / 10 GeV");
	make_hist("photon1Eta","photon 1 Eta",26,-2.6,2.6,"#eta","Events / 0.2");
	make_hist("photon1IsConv","photon 1 IsConv",2,-0.5,1.5,"","");
	make_hist("photon1HoverE","photon 1 HoverE",100,0,0.1,"H/E","Events / 0.001");
	make_hist("photon1SigmaIEtaIEta","photon 1 sigmaIetaIeta",80,0,0.04,"#sigma_{i#etai#eta}","Events / 0.0005");
	make_hist("photon1RelIso","photon 1 relative Isolation",30,-0.1,0.2,"RelIso","Events / 0.01");
	make_hist("photon1ChHadIso","photon 1 Charged Had Isolation",50,-2.5,10,"ChHadIso","Events / 0.25 GeV");
	make_hist("photon1NeuHadIso","photon 1 Neutral Had Isolation",60,-5,10,"NeuHadIso","Events / 0.25 GeV");
	make_hist("photon1PhoIso","photon 1 Photon Isolation",60,-5,10,"PhoIso","Events / 0.25 GeV");
	make_hist("photon1DrElectron","dR photon 1 to closest electron",60,0,6,"#DeltaR(#gamma,e)","Events / 0.1");
	make_hist("photon1DrJet","dR photon 1 to closest jet",60,0,6,"#DeltaR(#gamma,jet)","Events / 0.1");
	make_hist("photon1_j_inv_mass","invariant mass photon 1 and closest jet",50,0,200,"M(#gamma,jet)","Events / 4 GeV");
	make_hist("photon1_matchJetCSV","CSV of jet matched to Photon",400,-2,2,"CSV btag","");
	make_hist("photon1MotherID","photon 1 mother PDG ID",35,-0.5,34.5,"#gamma mother PID","Events");
	make_hist("photon1GMotherID","photon 1 Gmother PDG ID",35,-0.5,34.5,"#gamma Gmother PID","Events");
	make_hist("photon1DrMCbquark","dR photon 1 to gen level b",40,0,2,"#DeltaR(#gamma,b_{MC})","Events / 0.05");
	make_hist("GenPhotonEt","Et of matched Gen Photon",20,0,200,"E_{T} (GeV)","Events / 10 GeV");
	make_hist("nPhotons","number of photons",15,-0.5,14.5,"N_{#gamma}","Events");

	make_hist("jjg_inv_mass","j j gamma invariant mass",40,0,400,"GeV","");
	make_hist("jje_inv_mass","j j electron invariant mass",40,0,400,"GeV","");
	make_hist("bjjg_inv_mass","bj j gamma invariant mass",40,0,400,"GeV","");
	make_hist("bjje_inv_mass","bj j electron invariant mass",40,0,400,"GeV","");

	make_hist("ele_met_g_trans_mass","ele met gamma transverse mass",40,0,400,"","");	
	make_hist("Njets_vetoed","number of jets vetoed",10,-0.5,9.5,"","");

	make_hist("Ht","Ht",30,0,1500,"H_{T} (GeV)","Events / 50 GeV");
	make_hist("MET","Missing Transverse Momentum",30,0,300,"MET (GeV)","Events / 10 GeV");
	make_hist("nVtx","Number of Primary Vertices",50,0.5,50.5,"N_{PV}","Events");
	make_hist("nJets","number of jets",15,-0.5,14.5,"N_{jets}","Events");
	make_hist("jet1Pt","jet 1 pt",50,0,500,"p_{T} (GeV)","Events / 10 GeV");
	make_hist("jet1Eta","jet 1 Eta",26,-2.6,2.6,"#eta","Events / 0.2");
	make_hist("jet1CHF","jet1 CHF",10,0,1,"","");
	make_hist("jet1CEF","jet1 CEF",10,0,1,"","");
	make_hist("jet1NHF","jet1 NHF",10,0,1,"","");
	make_hist("jet1NEF","jet1 NEF",10,0,1,"","");
	make_hist("jet1NCharged","jet1 N charged",20,0,100,"","");
	make_hist("jet1NConstituents","jet1 N constituents",20,0,100,"","");

	make_hist("jet2Pt","jet 2 pt",30,0,300,"p_{T} (GeV)","Events / 10 GeV");
	make_hist("jet2Eta","jet 2 Eta",26,-2.6,2.6,"#eta","Events / 0.2");
	make_hist("jet3Pt","jet 3 pt",30,0,300,"p_{T} (GeV)","Events / 10 GeV");
	make_hist("jet3Eta","jet 3 Eta",26,-2.6,2.6,"#eta","Events / 0.2");
	make_hist("jet4Pt","jet 4 pt",15,0,150,"p_{T} (GeV)","Events / 10 GeV");
	make_hist("jet4Eta","jet 4 Eta",26,-2.6,2.6,"#eta","Events / 0.2");
}

void Histogrammer::fill(Selector* selector_, EventTree* tree, double weight){

        if(selector_->Jets.size() >= 2 && selector_->Photons.size() > 0){
		hists["jjg_inv_mass"]->Fill(selector_->JJGammaMassMin, weight);
		hists["bjjg_inv_mass"]->Fill(selector_->bJJGammaMassMin, weight);
	}
	if(selector_->Jets.size() >= 2 && selector_->Electrons.size() > 0){
		hists["jje_inv_mass"]->Fill(selector_->JJEMassMin, weight);
		hists["bjje_inv_mass"]->Fill(selector_->bJJEMassMin, weight);
	}
        if(selector_->Photons.size() > 0 && selector_->Electrons.size() > 0) 
		hists["ele_met_g_trans_mass"]->Fill(selector_->EGammaMETtransMassMin, weight);
        hists["Njets_vetoed"]->Fill(selector_->VetoedJets.size(), weight);

	// electrons
	if( selector_->Electrons.size() > 0 ){
		int ind = selector_->Electrons[0];
		hists["ele1Pt"]->Fill( tree->elePt_[ind], weight );
		hists["ele1Eta"]->Fill( tree->eleSCEta_[ind], weight );
		hists["ele1RelIso"]->Fill( selector_->Ele03RelIso[ind], weight );
		hists["ele1MVA"]->Fill( tree->eleIDMVATrig_[ind], weight );
		hists["ele1D0"]->Fill( tree->eleD0_[ind], weight );
		hists["ele1Dz"]->Fill( tree->eleDz_[ind], weight );
		hists["ele1EoverP"]->Fill( tree->eleEoverP_[ind], weight );
		hists["ele1sigmaIetaIeta"]->Fill( tree->eleSigmaIEtaIEta_[ind], weight );
		hists["ele1MissHits"]->Fill( tree->eleMissHits_[ind], weight );
		hists["ele1DrJet"]->Fill( minDr(tree->eleSCEta_[ind], tree->elePhi_[ind], selector_->Jets, tree->jetEta_, tree->jetPhi_), weight );
		hists["WtransMass"]->Fill( TMath::Sqrt(2*(tree->elePt_[ind])*(tree->pfMET_)*( 1.0 - TMath::Cos(dR(0.0,tree->elePhi_[ind],0.0,tree->pfMETPhi_)) )), weight );
		if( tree->isData_ == 0 ){
			if( tree->eleGenIndex_[ind] >= 0 ){
				hists["ele1MotherID"]->Fill( tree->eleGenMomPID_[ind], weight );
				if( TMath::Abs(tree->eleGenMomPID_[ind]) == 11 )
					hists["ele1GMotherID"]->Fill( tree->eleGenGMomPID_[ind], weight );
			}
			else hists["ele1MotherID"]->Fill( 0.0, weight );
		}
		if( selector_->Electrons.size() > 1 ){
			ind = selector_->Electrons[1];
			hists["ele2Pt"]->Fill( tree->elePt_[ind], weight );
			hists["ele2RelIso"]->Fill( selector_->Ele03RelIso[ind], weight );
		}
		if(selector_->Photons.size() > 0){
			TLorentzVector ele;
			TLorentzVector pho;
			int phoi = selector_->Photons[0];
			ele.SetPtEtaPhiM(tree->elePt_[ind], tree->eleSCEta_[ind], tree->elePhi_[ind], 0.0);
			pho.SetPtEtaPhiM(tree->phoEt_[phoi], tree->phoEta_[phoi], tree->phoPhi_[phoi], 0.0);
			hists["ele1pho1Mass"]->Fill( (ele+pho).M(), weight);
		}
	}
	
	// Loose Electrons (if any)
	if(selector_->ElectronsLoose.size() > 0 && tree->isData_ == 0){
		int eleInd = selector_->ElectronsLoose[0];
		double mindr = 999;
		for( int mcI = 0; mcI < tree->nMC_; ++mcI){
			if( tree->mcPID[mcI] == 22 ){
				double thisdr = dR(tree->mcEta[mcI], tree->mcPhi[mcI], tree->eleSCEta_[eleInd], tree->elePhi_[eleInd]);
				if( mindr > thisdr ) mindr = thisdr;
			}
		}
		hists["looseEleDrGenPho"]->Fill(mindr, weight);
	}
	
	// photons
	hists["nPhotons"]->Fill(selector_->Photons.size(), weight);
	if( selector_->Photons.size() > 0 ){
		int ind = selector_->Photons[0];
		
		hists2d["photon1_Sigma_ChIso"]->Fill(tree->phoSigmaIEtaIEta_[ind],selector_->Pho03ChHadIso[ind], weight);
		hists2d["photon1_Sigma_ChSCRIso"]->Fill(tree->phoSigmaIEtaIEta_[ind],selector_->Pho03ChHadSCRIso[ind], weight);
		hists2d["photon1_Sigma_PhoIso"]->Fill(tree->phoSigmaIEtaIEta_[ind],selector_->Pho03PhoIso[ind], weight);
		hists2d["photon1_Sigma_PhoSCRIso"]->Fill(tree->phoSigmaIEtaIEta_[ind],selector_->Pho03PhoSCRIso[ind], weight);

		if(tree->phoEt_[ind]>=25 && tree->phoEt_[ind]<35){
			hists2d["photon1_25_35_Sigma_ChIso"]->Fill(tree->phoSigmaIEtaIEta_[ind],selector_->Pho03ChHadIso[ind], weight);
		}
		if(tree->phoEt_[ind]>=35 && tree->phoEt_[ind]<45){
			hists2d["photon1_35_45_Sigma_ChIso"]->Fill(tree->phoSigmaIEtaIEta_[ind],selector_->Pho03ChHadIso[ind], weight);
		}
		if(tree->phoEt_[ind]>=45 && tree->phoEt_[ind]<60){
			hists2d["photon1_45_60_Sigma_ChIso"]->Fill(tree->phoSigmaIEtaIEta_[ind],selector_->Pho03ChHadIso[ind], weight);
		}
		if(tree->phoEt_[ind]>=60){
			hists2d["photon1_60_up_Sigma_ChIso"]->Fill(tree->phoSigmaIEtaIEta_[ind],selector_->Pho03ChHadIso[ind], weight);
		}

		hists["photon1Et"]->Fill( tree->phoEt_[ind], weight );
		hists["photon1Eta"]->Fill( tree->phoEta_[ind], weight );
		hists["photon1IsConv"]->Fill( tree->phoIsConv_[ind], weight );
		hists["photon1ChHadIso"]->Fill( selector_->Pho03ChHadIso[ind], weight );
		hists["photon1NeuHadIso"]->Fill( selector_->Pho03NeuHadIso[ind], weight );
		hists["photon1PhoIso"]->Fill( selector_->Pho03PhoIso[ind], weight );
		hists["photon1RelIso"]->Fill( 
			(selector_->Pho03ChHadIso[ind] + selector_->Pho03NeuHadIso[ind] + selector_->Pho03PhoIso[ind]) 
			/ tree->phoEt_[ind], weight );
		hists["photon1HoverE"]->Fill( tree->phoHoverE_[ind], weight );
		hists["photon1SigmaIEtaIEta"]->Fill( tree->phoSigmaIEtaIEta_[ind], weight );
		hists["photon1DrElectron"]->Fill( minDr(tree->phoEta_[ind], tree->phoPhi_[ind], selector_->Electrons, tree->eleSCEta_, tree->elePhi_), weight );
		hists["photon1DrJet"]->Fill( minDr(tree->phoEta_[ind], tree->phoPhi_[ind], selector_->Jets, tree->jetEta_, tree->jetPhi_), weight );
		hists["photon1_j_inv_mass"]->Fill( 
			phoJetmass(tree->phoEt_[ind], tree->phoEta_[ind], tree->phoPhi_[ind], selector_->Jets, tree->jetPt_, tree->jetEta_, tree->jetPhi_), 
			weight);
		hists["photon1_matchJetCSV"]->Fill(
			matchedJetBtag(tree->phoEta_[ind], tree->phoPhi_[ind], selector_->JetCandsNoDR, tree->jetEta_, tree->jetPhi_, tree->jetCombinedSecondaryVtxBJetTags_), 
			weight);
		hists["WtransMass_photon"]->Fill( TMath::Sqrt(2*(tree->phoEt_[ind])*(tree->pfMET_)*( 1.0 - TMath::Cos(dR(0.0,tree->phoPhi_[ind],0.0,tree->pfMETPhi_)) )), weight );

		if( tree->isData_ == 0 ){
			if( tree->phoGenIndex_[ind] >= 0 ){
				hists["photon1MotherID"]->Fill( tree->phoGenMomPID[ind], weight );
				if( TMath::Abs(tree->phoGenMomPID[ind]) == 22 ){
					hists["photon1GMotherID"]->Fill( tree->phoGenGMomPID[ind], weight );
				}
			}
			else {
				hists["photon1MotherID"]->Fill( 0.0, weight );
			}
			
			if(1){ //use alternative method for the time being
				// find the closest b-jet
				int phoGen=-1;
				double mindr = minDrPhoB(ind, tree);
				for( int mcI = 0; mcI < tree->nMC_; ++mcI){
					if( tree->mcIndex[mcI] == tree->phoGenIndex_[ind] ) 
						phoGen=mcI;
				}
				if( phoGen > 0)
					hists["GenPhotonEt"]->Fill(tree->mcPt[phoGen], weight);
				
				if(mindr<999) {
					hists["photon1DrMCbquark"]->Fill( mindr, weight );
					//if(mindr>0.2 && (tree->mcParentage[phoGen]==2 || (tree->mcParentage[phoGen]==10 && abs(tree->mcMomPID[phoGen])==24))){
					//	std::cout << "dr " << mindr << "  mcMomPID " << tree->mcMomPID[phoGen] 
					//		<< "  mcParentage " << tree->mcParentage[phoGen]
					//		<< "  mcPt " << tree->mcPt[phoGen] << std::endl;
					//}
				}
			}
			else{// dR that works for all ntuples
				double mindr = 999.0;
				for( int jetI = 0; jetI < tree->nJet_; ++jetI){
					if(abs(tree->jetPartonID_[jetI]) == 5){
						double dr = dR(tree->phoEta_[ind], tree->phoPhi_[ind], tree->jetGenJetEta_[jetI], tree->jetGenJetPhi_[jetI]);
						if( dr < mindr) mindr = dr;
					}
				}
				hists["photon1DrMCbquark"]->Fill( mindr, weight );
			}
		}
	}
	
	hists["Ht"]->Fill( calc_ht(selector_, tree), weight );
	hists["MET"]->Fill( tree->pfMET_, weight );
	hists["nVtx"]->Fill( tree->nVtx_, weight );
	hists["nJets"]->Fill( selector_->Jets.size(), weight );
	
	// jets
	if( selector_->Jets.size() > 0 ){
		int ind = selector_->Jets[0];
		hists["jet1Pt"]->Fill( tree->jetPt_[ind], weight );
		hists["jet1Eta"]->Fill( tree->jetEta_[ind], weight );
		hists["jet1CHF"]->Fill( tree->jetCHF_[ind], weight );
		hists["jet1CEF"]->Fill( tree->jetCEF_[ind], weight );
		hists["jet1NHF"]->Fill( tree->jetNHF_[ind], weight );
		hists["jet1NEF"]->Fill( tree->jetNEF_[ind], weight );
		hists["jet1NCharged"]->Fill( tree->jetNCharged_[ind], weight );
		hists["jet1NConstituents"]->Fill( tree->jetNConstituents_[ind], weight );
	}
	if( selector_->Jets.size() > 1 ){
		int ind = selector_->Jets[1];
		hists["jet2Pt"]->Fill( tree->jetPt_[ind], weight );
		hists["jet2Eta"]->Fill( tree->jetEta_[ind], weight );
	}
	if( selector_->Jets.size() > 2 ){
		int ind = selector_->Jets[2];
		hists["jet3Pt"]->Fill( tree->jetPt_[ind], weight );
		hists["jet3Eta"]->Fill( tree->jetEta_[ind], weight );
	}
	if( selector_->Jets.size() > 3 ){
		int ind = selector_->Jets[3];
		hists["jet4Pt"]->Fill( tree->jetPt_[ind], weight );
		hists["jet4Eta"]->Fill( tree->jetEta_[ind], weight );
	}
	
}

int Histogrammer::minDrIndex(double myEta, double myPhi, std::vector<int> Inds, Float_t* etas, Float_t* phis){
	double mindr = 999.0;
	double dr;
	int bestInd = -1;
	for( std::vector<int>::iterator it = Inds.begin(); it != Inds.end(); ++it){
		dr = dR(myEta, myPhi, etas[*it], phis[*it]);
		if( mindr > dr ) {
			mindr = dr;
			bestInd = *it;
		}
	}
	return bestInd;
}

double Histogrammer::minDr(double myEta, double myPhi, std::vector<int> Inds, Float_t* etas, Float_t* phis){
	int ind = minDrIndex(myEta, myPhi, Inds, etas, phis);
	if(ind>=0) return dR(myEta, myPhi, etas[ind], phis[ind]);
	else return 999.0;
}

double Histogrammer::matchedJetBtag(double phoEta, double phoPhi, std::vector<int> Inds, Float_t* jetEtas, Float_t* jetPhis, Float_t* jetCSV){
	int ind = minDrIndex(phoEta, phoPhi, Inds, jetEtas, jetPhis);
	if(ind>=0){
		if(jetCSV[ind]<-9) return -1.5;
		else return jetCSV[ind];
	}
	else return 1.9;
}

double Histogrammer::phoJetmass(double phoEt, double phoEta, double phoPhi, std::vector<int> Inds, Float_t* jetPts, Float_t* jetEtas, Float_t* jetPhis){
	double mindr = 999.0;
	double dr;
	int selected_jet=-1;
	TLorentzVector photonV;
	TLorentzVector jetV;
	photonV.SetPtEtaPhiM(phoEt, phoEta, phoPhi, 0.0);
	for( std::vector<int>::iterator it = Inds.begin(); it != Inds.end(); ++it){
		dr = dR(phoEta, phoPhi, jetEtas[*it], jetPhis[*it]);
		if( mindr > dr ) {
			mindr = dr;
			selected_jet = *it;
			jetV.SetPtEtaPhiM(jetPts[*it], jetEtas[*it], jetPhis[*it], 0.0);
		}
	}
	if(selected_jet>=0) return (jetV+photonV).M();
	else return 0.0;
}

double Histogrammer::minDrPhoB(int PhoInd, EventTree* tree){
	// find the closest b-jet
	TLorentzVector topVec;
	TLorentzVector WtopVec;
	TLorentzVector topBarVec;
	TLorentzVector WtopBarVec;
	TLorentzVector bTop;
	TLorentzVector bTopBar;
	int phoGen=-1;
	double mindr = 999.0;
	for( int mcI = 0; mcI < tree->nMC_; ++mcI){
		if( tree->mcIndex[mcI] == tree->phoGenIndex_[PhoInd] ) 
			phoGen=mcI;
		if( tree->mcPID[mcI] == 6) 
			topVec.SetPtEtaPhiM(tree->mcPt[mcI], tree->mcEta[mcI], tree->mcPhi[mcI], tree->mcMass[mcI]);
		if( tree->mcPID[mcI] == -6) 
			topBarVec.SetPtEtaPhiM(tree->mcPt[mcI], tree->mcEta[mcI], tree->mcPhi[mcI], tree->mcMass[mcI]);
		if( abs(tree->mcPID[mcI]) == 24 && tree->mcMomPID[mcI] == 6) 
			WtopVec.SetPtEtaPhiM(tree->mcPt[mcI], tree->mcEta[mcI], tree->mcPhi[mcI], tree->mcMass[mcI]);
		if( abs(tree->mcPID[mcI]) == 24 && tree->mcMomPID[mcI] == -6) 
			WtopBarVec.SetPtEtaPhiM(tree->mcPt[mcI], tree->mcEta[mcI], tree->mcPhi[mcI], tree->mcMass[mcI]);
	}
	if( phoGen > 0 && topVec.Pt() > 0.0001 && topBarVec.Pt() > 0.0001 && WtopVec.Pt() > 0.0001 && WtopBarVec.Pt() > 0.0001) {
		bTop = topVec - WtopVec;
		bTopBar = topBarVec - WtopBarVec;
		mindr = std::min(dR(tree->mcEta[phoGen], tree->mcPhi[phoGen], bTop.Eta(), bTop.Phi()),
				dR(tree->mcEta[phoGen], tree->mcPhi[phoGen], bTopBar.Eta(), bTopBar.Phi()));
	}
	return mindr;
}

double Histogrammer::calc_ht(Selector* selector_, EventTree* tree){
	double ht = 0.0;
	ht += tree->pfMET_;
	for( std::vector<int>::iterator it = selector_->Jets.begin(); it != selector_->Jets.end(); ++it)
		ht += tree->jetPt_[*it];
	for( std::vector<int>::iterator it = selector_->Electrons.begin(); it != selector_->Electrons.end(); ++it)
		ht += tree->elePt_[*it];
        for( std::vector<int>::iterator it = selector_->ElectronsLoose.begin(); it != selector_->ElectronsLoose.end(); ++it)
                ht += tree->elePt_[*it];
	for( std::vector<int>::iterator it = selector_->Muons.begin(); it != selector_->Muons.end(); ++it)
		ht += tree->muPt_[*it];	
        for( std::vector<int>::iterator it = selector_->MuonsLoose.begin(); it != selector_->MuonsLoose.end(); ++it)
                ht += tree->muPt_[*it];			
	for( std::vector<int>::iterator it = selector_->Photons.begin(); it != selector_->Photons.end(); ++it)
                ht += tree->phoEt_[*it];
	return ht;
}

void Histogrammer::make_hist(const char* hname, const char* htitle, int nbins, double xlow, double xhigh, const char* xlabel, const char* ylabel){
	TH1F* h = new TH1F(hname, htitle, nbins, xlow, xhigh);
	h->GetXaxis()->SetTitle(xlabel);
	h->GetYaxis()->SetTitle(ylabel);
	h->SetDirectory(0);
	//h->Sumw2();
	hists[hname] = h;
}

void Histogrammer::write_histograms(const char* folder, std::vector<TH1F*> histVector){
	std::string folderS(folder);

	TFile* outFile = new TFile((folderS+"/hist_"+title+".root").c_str(), "RECREATE");
	
	for( std::map< std::string, TH1F* >::iterator it = hists.begin(); it != hists.end(); ++it){
		it->second->SetDirectory(outFile->GetDirectory(""));
		it->second->Write();
		it->second->SetDirectory(0);
	}

	for( std::vector<TH1F*>::iterator it = histVector.begin(); it != histVector.end(); ++it){
		(*it)->SetDirectory(outFile->GetDirectory(""));
		(*it)->Write();
		(*it)->SetDirectory(0);
	}

	for( std::map< std::string, TH2F* >::iterator it = hists2d.begin(); it != hists2d.end(); ++it){
		it->second->SetDirectory(outFile->GetDirectory(""));
		it->second->Write();
		it->second->SetDirectory(0);
	}

	
	outFile->Close();
}

Histogrammer::~Histogrammer(){
}
