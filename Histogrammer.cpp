#include"Histogrammer.h"

Histogrammer::Histogrammer(std::string titleIn){
	title = titleIn;
	
	// 2d histograms
	make_hist2d("photon1_Sigma_ChIso","photon1 SigmaIetaIeta vs ChIso",160,0,0.04,300,-10,20);
	make_hist2d("photon1_Sigma_ChSCRIso","photon1 SigmaIetaIeta vs ChSCRIso",160,0,0.04,300,-10,20);
	make_hist2d("photon1_Sigma_PhoIso","photon1 SigmaIetaIeta vs PhoIso",160,0,0.04,300,-10,20);
	make_hist2d("photon1_Sigma_PhoSCRIso","photon1 SigmaIetaIeta vs PhoSCRIso",160,0,0.04,300,-10,20);
	make_hist2d("photon1_25_35_Sigma_ChSCRIso","photon1 Et 25 to 35 SigmaIetaIeta vs ChSCRIso",160,0,0.04,30,-10,20);
	make_hist2d("photon1_35_45_Sigma_ChSCRIso","photon1 Et 35 to 45 SigmaIetaIeta vs ChSCRIso",160,0,0.04,30,-10,20);
	make_hist2d("photon1_45_60_Sigma_ChSCRIso","photon1 Et 45 to 60 SigmaIetaIeta vs ChSCRIso",160,0,0.04,30,-10,20);
	make_hist2d("photon1_60_up_Sigma_ChSCRIso","photon1 Et 60 up SigmaIetaIeta vs ChSCRIso",160,0,0.04,30,-10,20);

	// creating histograms
	// electrons
	make_hist("ele1Pt","electron 1 Pt",30,0,300,"Electron p_{T} (GeV)","Events / 10 GeV");
	make_hist("ele1Eta","electron 1 Eta",26,-2.6,2.6,"Electron #eta","Events / 0.2");
	make_hist("ele1RelIso","electron 1 relative isolation",12,0,0.12,"Electron RelIso","Events / 0.01");
	make_hist("ele1MVA","electron 1 MVA",80,0.5,1.3,"Electron MVA Trig","Events / 0.01");
	make_hist("ele1D0","electron 1 Dxy_PV",80,-0.2,0.2,"Electron Dxy_PV (cm)","Events / 0.005 cm");
	make_hist("ele1Dz","electron 1 Dz",75,-0.15,0.15,"Electron D_{Z} (cm)","Events / 0.004 cm");
	make_hist("ele1EoverP","electron 1 EoverP",25,0,5,"Electron E/P","Events / 0.2");
	make_hist("ele1sigmaIetaIeta","electron 1 sigmaIetaIeta",80,0,0.04,"Electron #sigma_{i#etai#eta}","Events / 0.0005");
	make_hist("ele1MissHits","electron 1 missing hits",10,-0.5,9.5,"Electron missing hits","Events");
	make_hist("ele1DrJet","dR electron 1 to closest jet",60,0,6,"Electron #DeltaR(e,jet)","Events / 0.1");
	make_hist("ele1MotherID","electron 1 mother PDG ID",35,-0.5,34.5,"Electron mother PID","Events");
	make_hist("ele1GMotherID","electron 1 Gmother PDG ID",35,-0.5,34.5,"Electron Gmother PID","Events");
		
	make_hist("ele2Pt","electron 2 Pt",30,0,300,"p_{T} (GeV)","Events / 10 GeV");
	make_hist("ele2RelIso","electron 2 relative isolation",10,0,0.5,"RelIso","Events / 0.05");
	
	// photons
	make_hist("photon1Et","photon 1 Et",20,0,200,"Photon E_{T} (GeV)","Events / 10 GeV");
	make_hist("photon1Eta","photon 1 Eta",26,-2.6,2.6,"Photon #eta","Events / 0.2");
	make_hist("photon1IsConv","photon 1 IsConv",2,-0.5,1.5,"","");
	make_hist("photon1HoverE","photon 1 HoverE",100,0,0.1,"Photon H/E","Events / 0.001");
	make_hist("photon1SigmaIEtaIEta","photon 1 sigmaIetaIeta",80,0,0.04,"Photon #sigma_{i#etai#eta}","Events / 0.0005");
	make_hist("photon1RelIso","photon 1 relative Isolation",30,-0.1,0.2,"Photon RelIso","Events / 0.01");
	make_hist("photon1ChHadIso","photon 1 Charged Had Isolation",60,-2.0,10,"Photon ChHadIso","Events / 0.2 GeV");
	make_hist("photon1ChHadSCRIso","photon 1 Charged Had SCR Isolation",60,-2.0,10,"Photon ChHadSCRIso","Events / 0.2 GeV");
	make_hist("photon1ChHadRandIso","photon 1 Charged Had Rand Isolation",60,-2.0,10,"Photon ChHadRandIso","Events / 0.2 GeV");
	make_hist("photon1NeuHadIso","photon 1 Neutral Had Isolation",75,-5,10,"Photon NeuHadIso","Events / 0.2 GeV");
	make_hist("photon1PhoIso","photon 1 Photon Isolation",75,-5,10,"Photon PhoIso","Events / 0.2 GeV");
	make_hist("photon1PhoSCRIso","photon 1 Photon SCR Isolation",75,-5,10,"Photon PhoSCRIso","Events / 0.2 GeV");
	make_hist("photon1PhoRandIso","photon 1 Photon Rand Isolation",75,-5,10,"Photon PhoRandIso","Events / 0.2 GeV");
	make_hist("photon1DrElectron","dR photon 1 to closest electron",60,0,6,"#DeltaR(#gamma,e)","Events / 0.1");
	make_hist("photon1DrJet","dR photon 1 to closest jet",60,0,6,"#DeltaR(#gamma,jet)","Events / 0.1");
	make_hist("photon1MotherID","photon 1 mother PDG ID",35,-0.5,34.5,"Photon mother PID","Events");
	make_hist("photon1GMotherID","photon 1 Gmother PDG ID",35,-0.5,34.5,"Photon Gmother PID","Events");
	make_hist("photon1DrMCbquark","dR photon 1 to gen level b",40,0,2,"#DeltaR(#gamma,b_{MC})","Events / 0.05");
	make_hist("GenPhotonEt","Et of matched Gen Photon",20,0,200,"Gen Photon E_{T} (GeV)","Events / 10 GeV");
	make_hist("nPhotons","number of photons",15,-0.5,14.5,"N_{#gamma}","Events");

	// jets
	make_hist("jet1Pt","jet 1 pt",50,0,500,"p_{T} (GeV)","Events / 10 GeV");
	make_hist("jet1Eta","jet 1 Eta",26,-2.6,2.6,"#eta","Events / 0.2");
	make_hist("jet2Pt","jet 2 pt",30,0,300,"p_{T} (GeV)","Events / 10 GeV");
	make_hist("jet2Eta","jet 2 Eta",26,-2.6,2.6,"#eta","Events / 0.2");
	make_hist("jet3Pt","jet 3 pt",30,0,300,"p_{T} (GeV)","Events / 10 GeV");
	make_hist("jet3Eta","jet 3 Eta",26,-2.6,2.6,"#eta","Events / 0.2");
	make_hist("jet4Pt","jet 4 pt",15,0,150,"p_{T} (GeV)","Events / 10 GeV");
	make_hist("jet4Eta","jet 4 Eta",26,-2.6,2.6,"#eta","Events / 0.2");
	
	// event
	make_hist("looseEleDrGenPho","dR loose electron to Gen Photon",60,0,6,"#DeltaR(e_{loose},#gamma_{MC})","Events / 0.1");
	make_hist("WtransMass","W transverse mass",20,0,200,"M(W_{T})(GeV)","Events / 10 GeV");
	make_hist("ele1pho1Mass","electron + photon mass",20,0,200,"M(e,#gamma)(GeV)","Events / 10 GeV");	
	make_hist("Ht","Ht",30,0,1500,"H_{T} (GeV)","Events / 50 GeV");
	make_hist("MET","Missing Transverse Momentum",30,0,300,"MET (GeV)","Events / 10 GeV");
	make_hist("nVtx","Number of Primary Vertices",50,0.5,50.5,"N_{PV}","Events");
	make_hist("nJets","number of jets",15,-0.5,14.5,"N_{jets}","Events");
}

void Histogrammer::fill(Selector* selector, EventPick* selEvent, EventTree* tree, double weight){
	// 2d photon candidate histograms
	//std::cout << "here0" << std::endl;
	if(selEvent->PhotonsPresel.size()>0){
		int candArrInd = -1;
		int candInd = -1;
		for(int phoItmp = 0; phoItmp < selEvent->PhotonsPresel.size(); phoItmp++){
			if((int)selEvent->PhoPassChHadIso[phoItmp] + 
			   (int)selEvent->PhoPassPhoIso[phoItmp]+
			   (int)selEvent->PhoPassSih[phoItmp] >= 1){
				// at least 1 cut passed.
				candArrInd = selEvent->PhotonsPresel[phoItmp];
				candInd = phoItmp;
				break;
			}
		}
		//std::cout << "here01" << std::endl;
		if(candInd >= 0 && selEvent->PhoPassPhoIso[candInd]){
			hists2d["photon1_Sigma_ChIso"]->Fill(tree->phoSigmaIEtaIEta_->at(candArrInd),selector->Pho03ChHadIso[candArrInd], weight);
			hists2d["photon1_Sigma_ChSCRIso"]->Fill(tree->phoSigmaIEtaIEta_->at(candArrInd),selector->Pho03ChHadSCRIso[candArrInd], weight);
		
			if(tree->phoEt_->at(candArrInd)>=25 && tree->phoEt_->at(candArrInd)<35){
				hists2d["photon1_25_35_Sigma_ChSCRIso"]->Fill(tree->phoSigmaIEtaIEta_->at(candArrInd),selector->Pho03ChHadSCRIso[candArrInd], weight);
			}
			if(tree->phoEt_->at(candArrInd)>=35 && tree->phoEt_->at(candArrInd)<45){
				hists2d["photon1_35_45_Sigma_ChSCRIso"]->Fill(tree->phoSigmaIEtaIEta_->at(candArrInd),selector->Pho03ChHadSCRIso[candArrInd], weight);
			}
			if(tree->phoEt_->at(candArrInd)>=45 && tree->phoEt_->at(candArrInd)<60){
				hists2d["photon1_45_60_Sigma_ChSCRIso"]->Fill(tree->phoSigmaIEtaIEta_->at(candArrInd),selector->Pho03ChHadSCRIso[candArrInd], weight);
			}
			if(tree->phoEt_->at(candArrInd)>=60){
				hists2d["photon1_60_up_Sigma_ChSCRIso"]->Fill(tree->phoSigmaIEtaIEta_->at(candArrInd),selector->Pho03ChHadSCRIso[candArrInd], weight);
			}
		}
		if(candInd >= 0 && selEvent->PhoPassChHadIso[candInd]){
			hists2d["photon1_Sigma_PhoIso"]->Fill(tree->phoSigmaIEtaIEta_->at(candArrInd),selector->Pho03PhoIso[candArrInd], weight);
			hists2d["photon1_Sigma_PhoSCRIso"]->Fill(tree->phoSigmaIEtaIEta_->at(candArrInd),selector->Pho03PhoSCRIso[candArrInd], weight);
		}
		
	}
	//std::cout << "here1" << std::endl;
	// full event selection histograms
	
	// electrons
	if( selEvent->Electrons.size() > 0 ){
		int ind = selEvent->Electrons[0];
		hists["ele1Pt"]->Fill( tree->elePt_->at(ind), weight );
		hists["ele1Eta"]->Fill( tree->eleSCEta_->at(ind), weight );
		hists["ele1RelIso"]->Fill( selector->Ele03RelIso[ind], weight );
		hists["ele1MVA"]->Fill( tree->eleIDMVATrig_->at(ind), weight );
		hists["ele1D0"]->Fill( tree->eleD0_->at(ind), weight );
		hists["ele1Dz"]->Fill( tree->eleDz_->at(ind), weight );
		hists["ele1EoverP"]->Fill( tree->eleEoverP_->at(ind), weight );
		hists["ele1sigmaIetaIeta"]->Fill( tree->eleSigmaIEtaIEta_->at(ind), weight );
		hists["ele1MissHits"]->Fill( tree->eleMissHits_->at(ind), weight );
		hists["ele1DrJet"]->Fill( minDr(tree->eleSCEta_->at(ind), tree->elePhi_->at(ind), selEvent->Jets, tree->jetEta_, tree->jetPhi_), weight );
		hists["WtransMass"]->Fill( 
			TMath::Sqrt(2*(tree->elePt_->at(ind))*(tree->pfMET_)*( 1.0 - TMath::Cos(dR(0.0,tree->elePhi_->at(ind),0.0,tree->pfMETPhi_)) )), 
			weight );
		if( tree->isData_ == 0 ){
			if( tree->eleGenIndex_->at(ind) >= 0 ){
				hists["ele1MotherID"]->Fill( tree->eleGenMomPID_->at(ind), weight );
				if( TMath::Abs(tree->eleGenMomPID_->at(ind)) == 11 )
					hists["ele1GMotherID"]->Fill( tree->eleGenGMomPID_->at(ind), weight );
			}
			else hists["ele1MotherID"]->Fill( 0.0, weight );
		}
		if( selEvent->Electrons.size() > 1 ){
			ind = selEvent->Electrons[1];
			hists["ele2Pt"]->Fill( tree->elePt_->at(ind), weight );
			hists["ele2RelIso"]->Fill( selector->Ele03RelIso[ind], weight );
		}
		
		if(selEvent->Photons.size() > 0){
			TLorentzVector ele;
			TLorentzVector pho;
			int phoi = selEvent->Photons[0];
			ele.SetPtEtaPhiM(tree->elePt_->at(ind), tree->eleSCEta_->at(ind), tree->elePhi_->at(ind), 0.0);
			pho.SetPtEtaPhiM(tree->phoEt_->at(phoi), tree->phoEta_->at(phoi), tree->phoPhi_->at(phoi), 0.0);
			hists["ele1pho1Mass"]->Fill( (ele+pho).M(), weight);
		}
	}
	//std::cout << "here2" << std::endl;
	// Loose Electrons (if any)
	if(selEvent->ElectronsLoose.size() > 0 && tree->isData_ == 0){
		int eleInd = selEvent->ElectronsLoose[0];
		double mindr = 999;
		for( int mcI = 0; mcI < tree->nMC_; ++mcI){
			if( tree->mcPID->at(mcI) == 22 ){
				double thisdr = dR(tree->mcEta->at(mcI), tree->mcPhi->at(mcI), tree->eleSCEta_->at(eleInd), tree->elePhi_->at(eleInd));
				if( mindr > thisdr ) mindr = thisdr;
			}
		}
		hists["looseEleDrGenPho"]->Fill(mindr, weight);
	}
	//std::cout << "here3" << std::endl;
	// photons
	hists["nPhotons"]->Fill(selEvent->Photons.size(), weight);
	if( selEvent->Photons.size() > 0 ){
		int ind = selEvent->Photons[0];

		hists["photon1Et"]->Fill( tree->phoEt_->at(ind), weight );
		hists["photon1Eta"]->Fill( tree->phoEta_->at(ind), weight );
		hists["photon1IsConv"]->Fill( tree->phoIsConv_->at(ind), weight );
		
		hists["photon1ChHadIso"]->Fill( selector->Pho03ChHadIso[ind], weight );
		hists["photon1ChHadSCRIso"]->Fill( selector->Pho03ChHadSCRIso[ind], weight );
		hists["photon1ChHadRandIso"]->Fill( selector->Pho03RandChHadIso[ind], weight );
		
		hists["photon1NeuHadIso"]->Fill( selector->Pho03NeuHadIso[ind], weight );
		
		hists["photon1PhoIso"]->Fill( selector->Pho03PhoIso[ind], weight );
		hists["photon1PhoSCRIso"]->Fill( selector->Pho03PhoSCRIso[ind], weight );
		hists["photon1PhoRandIso"]->Fill( selector->Pho03RandPhoIso[ind], weight );
		
		hists["photon1RelIso"]->Fill( 
			(selector->Pho03ChHadIso[ind] + selector->Pho03NeuHadIso[ind] + selector->Pho03PhoIso[ind]) 
			/ tree->phoEt_->at(ind), weight );
		hists["photon1HoverE"]->Fill( tree->phoHoverE_->at(ind), weight );
		hists["photon1SigmaIEtaIEta"]->Fill( tree->phoSigmaIEtaIEta_->at(ind), weight );
		hists["photon1DrElectron"]->Fill( minDr(tree->phoEta_->at(ind), tree->phoPhi_->at(ind), selEvent->Electrons, tree->eleSCEta_, tree->elePhi_), weight );
		hists["photon1DrJet"]->Fill( minDr(tree->phoEta_->at(ind), tree->phoPhi_->at(ind), selector->Jets, tree->jetEta_, tree->jetPhi_), weight );

		if( tree->isData_ == 0 ){
			if( tree->phoGenIndex_->at(ind) >= 0 ){
				hists["photon1MotherID"]->Fill( tree->phoGenMomPID_->at(ind), weight );
				if( TMath::Abs(tree->phoGenMomPID_->at(ind)) == 22 ){
					hists["photon1GMotherID"]->Fill( tree->phoGenGMomPID_->at(ind), weight );
				}
			}
			else {
				hists["photon1MotherID"]->Fill( 0.0, weight );
			}
			
			// find the closest b-jet
			double mindr = minDrPhoB(ind, tree);
			int phoGen=-1;
			for( int mcI = 0; mcI < tree->nMC_; ++mcI){
				if( tree->mcIndex->at(mcI) == tree->phoGenIndex_->at(ind) ) 
					phoGen=mcI;
			}
			if( phoGen > 0)
				hists["GenPhotonEt"]->Fill(tree->mcPt->at(phoGen), weight);
				
			if(mindr<999) {
				hists["photon1DrMCbquark"]->Fill( mindr, weight );
			}
			
		}
	}
	//std::cout << "here4" << std::endl;
	hists["Ht"]->Fill( calc_ht(selEvent, tree), weight );
	hists["MET"]->Fill( tree->pfMET_, weight );
	hists["nVtx"]->Fill( tree->nVtx_, weight );
	hists["nJets"]->Fill( selEvent->Jets.size(), weight );
	//std::cout << "here5" << std::endl;
	// jets
	if( selEvent->Jets.size() > 0 ){
		int ind = selEvent->Jets[0];
		hists["jet1Pt"]->Fill( tree->jetPt_->at(ind), weight );
		hists["jet1Eta"]->Fill( tree->jetEta_->at(ind), weight );
	}
	if( selEvent->Jets.size() > 1 ){
		int ind = selEvent->Jets[1];
		hists["jet2Pt"]->Fill( tree->jetPt_->at(ind), weight );
		hists["jet2Eta"]->Fill( tree->jetEta_->at(ind), weight );
	}
	if( selEvent->Jets.size() > 2 ){
		int ind = selEvent->Jets[2];
		hists["jet3Pt"]->Fill( tree->jetPt_->at(ind), weight );
		hists["jet3Eta"]->Fill( tree->jetEta_->at(ind), weight );
	}
	if( selEvent->Jets.size() > 3 ){
		int ind = selEvent->Jets[3];
		hists["jet4Pt"]->Fill( tree->jetPt_->at(ind), weight );
		hists["jet4Eta"]->Fill( tree->jetEta_->at(ind), weight );
	}
	
}

int Histogrammer::minDrIndex(double myEta, double myPhi, std::vector<int> Inds, std::vector<float> *etas, std::vector<float> *phis){
	double mindr = 999.0;
	double dr;
	int bestInd = -1;
	for( std::vector<int>::iterator it = Inds.begin(); it != Inds.end(); ++it){
		dr = dR(myEta, myPhi, etas->at(*it), phis->at(*it));
		if( mindr > dr ) {
			mindr = dr;
			bestInd = *it;
		}
	}
	return bestInd;
}

double Histogrammer::minDr(double myEta, double myPhi, std::vector<int> Inds, std::vector<float> *etas, std::vector<float> *phis){
	int ind = minDrIndex(myEta, myPhi, Inds, etas, phis);
	if(ind>=0) return dR(myEta, myPhi, etas->at(ind), phis->at(ind));
	else return 999.0;
}


double Histogrammer::minDrPhoB(int PhoInd, EventTree* tree){
	// find the closest b-jet
	TLorentzVector b;
	TLorentzVector bBar;
	int phoGen=-1;
	double mindr = 999.0;
	for( int mcI = 0; mcI < tree->nMC_; ++mcI){
		if( tree->mcIndex->at(mcI) == tree->phoGenIndex_->at(PhoInd) ) 
			phoGen=mcI;
		if( tree->mcPID->at(mcI) == 5) 
			b.SetPtEtaPhiM(tree->mcPt->at(mcI), tree->mcEta->at(mcI), tree->mcPhi->at(mcI), tree->mcMass->at(mcI));
		if( tree->mcPID->at(mcI) == -5) 
			bBar.SetPtEtaPhiM(tree->mcPt->at(mcI), tree->mcEta->at(mcI), tree->mcPhi->at(mcI), tree->mcMass->at(mcI));
	}
	if( phoGen > 0 && b.Pt() > 0.0001 && bBar.Pt() > 0.0001 ) {
		mindr = std::min(dR(tree->mcEta->at(phoGen), tree->mcPhi->at(phoGen), b.Eta(), b.Phi()),
						 dR(tree->mcEta->at(phoGen), tree->mcPhi->at(phoGen), bBar.Eta(), bBar.Phi()));
	}
	return mindr;
}

double Histogrammer::calc_ht(EventPick* evtPick, EventTree* tree){
	double ht = 0.0;
	ht += tree->pfMET_;
	for( std::vector<int>::iterator it = evtPick->Jets.begin(); it != evtPick->Jets.end(); ++it)
		ht += tree->jetPt_->at(*it);
	for( std::vector<int>::iterator it = evtPick->Electrons.begin(); it != evtPick->Electrons.end(); ++it)
		ht += tree->elePt_->at(*it);
        for( std::vector<int>::iterator it = evtPick->ElectronsLoose.begin(); it != evtPick->ElectronsLoose.end(); ++it)
                ht += tree->elePt_->at(*it);
	//for( std::vector<int>::iterator it = evtPick->Muons.begin(); it != evtPick->Muons.end(); ++it)
	//	ht += tree->muPt_->at(*it);	
	for( std::vector<int>::iterator it = evtPick->MuonsLoose.begin(); it != evtPick->MuonsLoose.end(); ++it)
		ht += tree->muPt_->at(*it);			
	// photons are now also in jet collection
	//for( std::vector<int>::iterator it = evtPick->Photons.begin(); it != evtPick->Photons.end(); ++it)
	//	ht += tree->phoEt_->at(*it);
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

void Histogrammer::make_hist2d(const char* hname, const char* htitle, int nxbins, double xlow, double xhigh, int nybins, double ylow, double yhigh){
	TH2F* h2 = new TH2F(hname, htitle, nxbins, xlow, xhigh, nybins, ylow, yhigh);
	h2->SetDirectory(0);
	hists2d[hname] = h2;
}

void Histogrammer::write_histograms(std::string folderS, std::vector<TH1F*> histVector){

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
	// do not delete histograms
}
