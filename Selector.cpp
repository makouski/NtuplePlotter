#include"Selector.h"

double dR(double eta1, double phi1, double eta2, double phi2){
	double dphi = phi2 - phi1;
	double deta = eta2 - eta1;
	static const double pi = TMath::Pi();
	dphi = TMath::Abs( TMath::Abs(dphi) - pi ) - pi;
	return TMath::Sqrt( dphi*dphi + deta*deta );
}

Selector::Selector(const char* titleIn){
	title.assign(titleIn);
	
	// histograms
	//std::string histName(titleIn);
	Nele_before_cut = new TH1F("Nele_before_cut","number of electrons before cut",10,-0.5,9.5);
	Nele_before_cut->SetDirectory(0);
	histVector.push_back(Nele_before_cut);
	
	// cut-flow
	cutFlow = new TH1F("cut_flow","cut flow",10,-0.5,9.5);
	cutFlow->SetDirectory(0);
	set_cutflow_labels(); // keep the labels close to the cuts definitions (below)
	histVector.push_back(cutFlow);

	// assign cut values
	MET_cut = 20.0;
	no_trigger = false;
	
	// jets
	jet_Pt_cuts.push_back(30);
	btag_cut = 0.679;
	Njet_cut = 3;
	NBjet_cut = 1;

	// electrons
	ele_Pt_cut = 35.0;
	ele_PtLoose_cut = 20.0;
	ele_RelIso_range[0] = 0.0;
	ele_RelIso_range[1] = 0.1;
	ele_RelIsoLoose_cut = 0.2;
	ele_MVA_range[0] = 0.9;
	ele_MVA_range[1] = 1.0;
	ele_MVALoose_cut = 0.0;
	ele_Dxy_cut = 0.02;
	ele_MissInnHit_cut = 0;
	Nele_cut = 1;
	NlooseEleVeto_cut = 99; // no cut

	// photons
	pho_Et_cut = 25.0; 
	pho_ID_ind = 0; // 0 - Loose, 1 - Medium, 2 - Tight
	pho_noSigmaIEta_cut = false;
	pho_noIso_cut = false;
	pho_noChHadIso_cut = false;
	pho_noPhoIso_cut = false;
	pho_noPixelSeed_cut = false;
        pho_noEleVeto_cut = false;
	Npho_cut = 1;
	
	// muons
	mu_PtLoose_cut = 10.0;
	mu_RelIsoLoose_cut = 0.2;
	NlooseMuVeto_cut = 0;
	
	// delta R cuts
	doVeto = true;
	veto_jet_dR = 0.3;
	veto_lep_jet_dR = 0.5;// no effect (jets cleaned before leptons)
	veto_pho_jet_dR = 0.7;
	veto_pho_lep_dR = 0.7;
	W_mass_veto = 0.0; // turned off by default
	W_trans_mass_veto= 0.0; // turned off by default
}

void Selector::clear_vectors(){
	Photons.clear(); VetoedPhotons.clear(); WVetoedPhotons.clear();
	Electrons.clear(); ElectronsLoose.clear(); VetoedElectrons.clear();
	Muons.clear(); MuonsLoose.clear();
	JetCandsNoDR.clear(); JetCands.clear(); Jets.clear(); VetoedJets.clear();
	Ele03RelIso.clear();
	Mu04RelIso.clear();
	Pho03ChHadIso.clear(); Pho03ChHadSCRIso.clear(); Pho03NeuHadIso.clear(); Pho03PhoIso.clear(); Pho03PhoSCRIso.clear();
}

void Selector::filter_photons(){
	for(int phoInd = 0; phoInd < tree->nPho_; ++phoInd){
		double eta = tree->phoEta_[phoInd];
		double et = tree->phoEt_[phoInd];
		//Pho03ChHadIso.push_back(  std::max( tree->phoPFChIso_[phoInd]  - tree->rho2012_ * phoEffArea03ChHad(eta),  0.0) );
		//Pho03NeuHadIso.push_back( std::max( tree->phoPFNeuIso_[phoInd] - tree->rho2012_ * phoEffArea03NeuHad(eta), 0.0) );
		//Pho03PhoIso.push_back(    std::max( tree->phoPFPhoIso_[phoInd] - tree->rho2012_ * phoEffArea03Pho(eta),    0.0) );
		Pho03ChHadIso.push_back(    tree->phoPFChIso_[phoInd]  - tree->rho2012_ * phoEffArea03ChHad(eta) );
		Pho03ChHadSCRIso.push_back( tree->phoSCRChIso_[phoInd]  - tree->rho2012_ * phoEffArea03ChHad(eta) );
		Pho03NeuHadIso.push_back( tree->phoPFNeuIso_[phoInd] - tree->rho2012_ * phoEffArea03NeuHad(eta) );
		Pho03PhoIso.push_back(    tree->phoPFPhoIso_[phoInd] - tree->rho2012_ * phoEffArea03Pho(eta) );
		Pho03PhoSCRIso.push_back( tree->phoSCRPhoIso_[phoInd] - tree->rho2012_ * phoEffArea03Pho(eta) );

		// manual spike cleaning
		if (dR(tree->phoEta_[phoInd], tree->phoPhi_[phoInd], -1.76, 1.37) < 0.05) continue;
		if (dR(tree->phoEta_[phoInd], tree->phoPhi_[phoInd],  2.37, 2.69) < 0.05) continue;		

		int region = 0; //barrel
		if(TMath::Abs( eta )>1.5) region = 1; //endcap
		if( fidEtaPass( eta ) &&
		   et > pho_Et_cut &&
		   ( pho_noPixelSeed_cut || tree->phohasPixelSeed_[phoInd] == 0 ) &&
		   ( pho_noEleVeto_cut || tree->phoEleVeto_[phoInd] == 0 ) &&
		   tree->phoIsConv_[phoInd] == photonID_IsConv[region][pho_ID_ind] &&
		   tree->phoHoverE_[phoInd] < photonID_HoverE[region][pho_ID_ind] &&
		   ( pho_noSigmaIEta_cut || tree->phoSigmaIEtaIEta_[phoInd] < photonID_SigmaIEtaIEta[region][pho_ID_ind] ) &&
		   ( pho_noIso_cut || 
			( ((pho_noChHadIso_cut && Pho03ChHadIso[phoInd] < 10) || Pho03ChHadIso[phoInd] < photonID_RhoCorrR03ChHadIso[region][pho_ID_ind] ) &&
			 Pho03NeuHadIso[phoInd] < photonID_RhoCorrR03NeuHadIso_0[region][pho_ID_ind] + et * photonID_RhoCorrR03NeuHadIso_1[region][pho_ID_ind] &&
			 ( pho_noPhoIso_cut || Pho03PhoIso[phoInd] < photonID_RhoCorrR03PhoIso_0[region][pho_ID_ind] + et * photonID_RhoCorrR03PhoIso_1[region][pho_ID_ind]) ) )
		   ){
			Photons.push_back(phoInd);
		}
	}
}

void Selector::filter_electrons(){
	for(int eleInd = 0; eleInd < tree->nEle_; ++eleInd){
		double eta = tree->eleSCEta_[eleInd];
		double pt = tree->elePt_[eleInd];
		Ele03RelIso.push_back( (tree->elePFChIso03_[eleInd] + 
					std::max(tree->elePFNeuIso03_[eleInd] + 
					tree->elePFPhoIso03_[eleInd] - 
					std::max(0.0,(double)tree->rho2012_) * eleEffArea03(eta), 0.0)
					) / pt );
		if( fidEtaPass( eta ) &&
		   //fabs(eta) < 1.5 &&
		   pt > ele_Pt_cut &&
		   //pt>25 && pt<35 &&
		   ele_RelIso_range[0] <= Ele03RelIso[eleInd] &&
		   Ele03RelIso[eleInd] < ele_RelIso_range[1] &&
		   ele_MVA_range[0] < tree->eleIDMVATrig_[eleInd] &&
		   tree->eleIDMVATrig_[eleInd] <= ele_MVA_range[1] &&
		   tree->eleConvVtxFit_[eleInd] == 0 &&
		   TMath::Abs(tree->eleD0_[eleInd]) < ele_Dxy_cut &&
		   tree->eleMissHits_[eleInd] <= ele_MissInnHit_cut
		   ){
			Electrons.push_back(eleInd);
		}
		else if( fabs(eta) < 2.5 && 
			pt > ele_PtLoose_cut && 
			Ele03RelIso[eleInd] < ele_RelIsoLoose_cut && 
			tree->eleIDMVATrig_[eleInd] > ele_MVALoose_cut 
			){ 
				ElectronsLoose.push_back(eleInd);
		}
	}
}

void Selector::filter_muons(){
	for(int muInd = 0; muInd < tree->nMu_; ++muInd){
		double eta = tree->muEta_[muInd];
		double pt = tree->muPt_[muInd];
		Mu04RelIso.push_back( (tree->muPFIsoR04_CH_[muInd] +
					std::max(tree->muPFIsoR04_NH_[muInd] +
					tree->muPFIsoR04_Pho_[muInd] - 
					std::max(0.0,(double)tree->rho2012_) * muEffArea04(eta), 0.0)
					) / pt );
		if( TMath::Abs(eta) < 2.5 &&
		   pt > mu_PtLoose_cut &&
		   Mu04RelIso[muInd] < mu_RelIsoLoose_cut
		   ){
			MuonsLoose.push_back(muInd);
		}
	}
}


void Selector::vetoCollection(std::vector<int>* toClean, const Float_t* etaClean, const Float_t* phiClean, 
							   const std::vector<int>* toCheck, const Float_t* etaCheck, const Float_t* phiCheck, 
							   std::set<int>* vetoed, double dRcut){
	for(std::vector<int>::iterator cleanInd = toClean->begin(); cleanInd != toClean->end(); ++cleanInd)
		for(std::vector<int>::const_iterator checkInd = toCheck->begin() ; checkInd != toCheck->end(); ++checkInd)
			if( dR(etaClean[*cleanInd], phiClean[*cleanInd], etaCheck[*checkInd], phiCheck[*checkInd]) < dRcut )
				vetoed->insert(*cleanInd);
}

void Selector::cleanCollection(std::vector<int>* toClean, std::set<int>* vetoed){
	for(std::set<int>::iterator vInd = vetoed->begin(); vInd != vetoed->end(); ++vInd)
		for(std::vector<int>::iterator cleanInd = toClean->begin(); cleanInd != toClean->end(); ++cleanInd)
			if( *cleanInd == *vInd ){
				toClean->erase(cleanInd);
				break;
			}
}

// remove jets that are too close to lepton or photon (they are not "jets")
// remove leptons that are nearby jets
void Selector::make_dR_cuts(){
	// check jets close to lepton or photon (they are not "jets")
	vetoCollection(&JetCands, tree->jetEta_, tree->jetPhi_, 
					&Electrons, tree->eleSCEta_, tree->elePhi_, 
					&VetoedJets, veto_jet_dR);
	vetoCollection(&JetCands, tree->jetEta_, tree->jetPhi_, 
					&Photons, tree->phoEta_, tree->phoPhi_, 
					&VetoedJets, veto_jet_dR);
	cleanCollection(&JetCands, &VetoedJets);
	// check photons and electrons nearby jets
	vetoCollection(&Electrons, tree->eleSCEta_, tree->elePhi_,
				   &JetCands, tree->jetEta_, tree->jetPhi_,
				   &VetoedElectrons, veto_lep_jet_dR);
	vetoCollection(&Photons, tree->phoEta_, tree->phoPhi_,
				   &JetCands, tree->jetEta_, tree->jetPhi_,
				   &VetoedPhotons, veto_pho_jet_dR);
	
	// check photons nearby electrons
	vetoCollection(&Photons, tree->phoEta_, tree->phoPhi_,
				   &Electrons, tree->eleSCEta_, tree->elePhi_,
				   &VetoedPhotons, veto_pho_lep_dR);
	
	// clean Electrons and Photons collections
	cleanCollection(&Electrons, &VetoedElectrons);
	cleanCollection(&Photons, &VetoedPhotons);
}

// remove photons that can be electron or jet FSR
// remove photons that can be W FSR
void Selector::make_mass_cuts(){
	// iterate over the rest of the objects to find system of (jet, jet, photon) that has inv. mass below some threshold
	JJGammaMassMin = 9999.9;
	JJEMassMin = 9999.9;
	bJJGammaMassMin = 9999.9;
	bJJEMassMin = 9999.9;
	if(JetCands.size() >= 2){
		bool btag1;
		bool btag2;
		int btags;
		// tried mass of 2 highest Pt jets + ele or pho
		//if(Photons.size() > 0) JJGammaMassMin = invMass_jjp(JetCands[0], JetCands[1], Photons[0]);
		//if(Electrons.size() > 0) JJEMassMin = invMass_jje(JetCands[0], JetCands[1], Electrons[0]);
		for(std::vector<int>::iterator jetInd1 = JetCands.begin(); jetInd1 != JetCands.end()-1; ++jetInd1){
			btag1 =  tree->jetCombinedSecondaryVtxBJetTags_[*jetInd1] >= btag_cut;
			if(btag1) btags++;
			for(std::vector<int>::iterator jetInd2 = jetInd1+1; jetInd2 != JetCands.end(); ++jetInd2){
				btag2 = tree->jetCombinedSecondaryVtxBJetTags_[*jetInd2] >= btag_cut;
				if (btag2) btags++;
				//double jjmass = invMass_jj(*jetInd1, *jetInd2);
				for(std::vector<int>::iterator phoInd = Photons.begin(); phoInd != Photons.end(); ++phoInd){
					double ginvm = invMass_jjp(*jetInd1, *jetInd2, *phoInd);
					if( /*btags==0 &&*/ (!btag1 && !btag2) && JJGammaMassMin > ginvm) JJGammaMassMin = ginvm;
					if( (btag1||btag2) && bJJGammaMassMin > ginvm) bJJGammaMassMin = ginvm;
					if( (!btag1 && !btag2) && ginvm <= W_mass_veto ) WVetoedPhotons.insert(*phoInd);
				}
				for(std::vector<int>::iterator eleInd = Electrons.begin(); eleInd != Electrons.end(); ++eleInd){
					double einvm = invMass_jje(*jetInd1, *jetInd2, *eleInd);
					if( (!btag1 && !btag2) && JJEMassMin > einvm) JJEMassMin = einvm;
					if( (btag1||btag2) && bJJEMassMin > einvm) bJJEMassMin = einvm;
				}
			}
		}
	}
	// remove photons that can be lep W FSR
	EGammaMETtransMassMin = 9999.9;
	for(std::vector<int>::iterator phoInd = Photons.begin(); phoInd != Photons.end(); ++phoInd)
		for(std::vector<int>::iterator eleInd = Electrons.begin() ; eleInd != Electrons.end(); ++eleInd){
			double transm = clusterTransMass_pe(*phoInd, *eleInd);
			if( EGammaMETtransMassMin > transm ) EGammaMETtransMassMin = transm;
			if( transm < W_trans_mass_veto ) WVetoedPhotons.insert(*phoInd);
		}
	// finally remove vetoed Photons
	cleanCollection(&Photons, &WVetoedPhotons);
	// done!
}

double Selector::invMass_jj(int jetInd1, int jetInd2){
	TLorentzVector j1,j2;
	j1.SetPtEtaPhiM(tree->jetPt_[jetInd1], tree->jetEta_[jetInd1], tree->jetPhi_[jetInd1], 0.0);
	j2.SetPtEtaPhiM(tree->jetPt_[jetInd2], tree->jetEta_[jetInd2], tree->jetPhi_[jetInd2], 0.0);
	return (j1+j2).M();
}

double Selector::invMass_jje(int jetInd1, int jetInd2, int eleInd){
	TLorentzVector j1,j2,ele;
	j1.SetPtEtaPhiM(tree->jetPt_[jetInd1], tree->jetEta_[jetInd1], tree->jetPhi_[jetInd1], 0.0);
	j2.SetPtEtaPhiM(tree->jetPt_[jetInd2], tree->jetEta_[jetInd2], tree->jetPhi_[jetInd2], 0.0);
	ele.SetPtEtaPhiM(tree->elePt_[eleInd], tree->eleSCEta_[eleInd], tree->elePhi_[eleInd],  0.0);
	return (j1+j2+ele).M();
}

double Selector::invMass_jjp(int jetInd1, int jetInd2, int phoInd){
	TLorentzVector j1,j2,pho;
	j1.SetPtEtaPhiM(tree->jetPt_[jetInd1], tree->jetEta_[jetInd1], tree->jetPhi_[jetInd1], 0.0);
	j2.SetPtEtaPhiM(tree->jetPt_[jetInd2], tree->jetEta_[jetInd2], tree->jetPhi_[jetInd2], 0.0);
	pho.SetPtEtaPhiM(tree->phoEt_[phoInd], tree->phoEta_[phoInd],  tree->phoPhi_[phoInd],  0.0);
	return (j1+j2+pho).M();
}

// sqrt( ( sqrt(cluster.Pt^2 + cluster.M^2) + |MET.Pt| )^2 - (cluster.Pt + MET.Pt)^2 )
// here "cluster.Pt + MET.Pt" is vector sum of transverse momenta
double Selector::clusterTransMass_pe(int phoInd, int eleInd){
	TLorentzVector ele,pho,met;
	ele.SetPtEtaPhiM(tree->elePt_[eleInd], tree->eleSCEta_[eleInd], tree->elePhi_[eleInd], 0.0);
	pho.SetPtEtaPhiM(tree->phoEt_[phoInd], tree->phoEta_[phoInd],  tree->phoPhi_[phoInd],  0.0);
	met.SetPtEtaPhiM(tree->pfMET_, 0.0, tree->pfMETPhi_, 0.0);
	TLorentzVector sum = ele+pho;
	return TMath::Sqrt( TMath::Power( TMath::Sqrt(sum.Perp2() + sum.M2()) + met.Perp(), 2) - (met+sum).Perp2() );
}

// jet ID is not likely to be altered, so it is hardcoded
void Selector::filter_jets(){
	for(int jetInd = 0; jetInd < tree->nJet_; ++jetInd){
		if( TMath::Abs(tree->jetEta_[jetInd]) < 2.4  &&
		   tree->jetPt_[jetInd]         > 20.0 &&
		   tree->jetCHF_[jetInd]        > 0    &&
		   tree->jetNHF_[jetInd]        < 0.99 &&
		   tree->jetCEF_[jetInd]        < 0.99 &&
		   tree->jetNEF_[jetInd]        < 0.99 &&
		   tree->jetNCharged_[jetInd]   > 0    &&
		   tree->jetNConstituents_[jetInd] > 1){
				JetCands.push_back(jetInd);
				JetCandsNoDR.push_back(jetInd); // jet collection withour dR cuts and veto
			}
	}
	
	if( doVeto ) {
		make_dR_cuts();
	}
	make_mass_cuts();
	
	for(int candInd = 0; candInd < JetCands.size(); ++candInd){
		if( tree->jetPt_[JetCands[candInd]] > JetPtCut(candInd) ) 
			Jets.push_back(JetCands[candInd]);
		else break;
	}
}

double Selector::JetPtCut(int jetInd){
	int size = jet_Pt_cuts.size();
	if( jetInd < size ) return jet_Pt_cuts[jetInd];
	return jet_Pt_cuts[size-1];
}

int Selector::accept_counts(){
	cutFlow->Fill(0);
	if(Photons.size() < Npho_cut) return 0; 
	cutFlow->Fill(1);
	Nele_before_cut->Fill(Electrons.size());
	if(Electrons.size() != Nele_cut) return 0;
	//if(Electrons.size() < Nele_cut) return 0;
	cutFlow->Fill(2);
	if( MuonsLoose.size() > NlooseMuVeto_cut ) return 0;
	cutFlow->Fill(3);
	if( ElectronsLoose.size() > NlooseEleVeto_cut ) return 0;
	cutFlow->Fill(4);
	if(Jets.size() < Njet_cut) return 0;
	cutFlow->Fill(5);
	int bjets = 0;
	for(std::vector<int>::iterator jetInd = Jets.begin(); jetInd != Jets.end(); ++jetInd)
		if(tree->jetCombinedSecondaryVtxBJetTags_[*jetInd] > btag_cut) bjets++;
	if(bjets < NBjet_cut) return 0;
	cutFlow->Fill(6);
	if(tree->pfMET_ < MET_cut) return 0;
	cutFlow->Fill(7);
	if(!no_trigger && tree->HLT_[tree->HLTIndex_[17]]==0) return 0; // 17 --> HLT_Ele27_WP80_v
	cutFlow->Fill(8);

	return 1;
}

void Selector::print_cutflow(){
	std::cout << "Cut-Flow for the selector:   " << title << std::endl;
	std::cout << "Input Events                 " << cutFlow->GetBinContent(1) << std::endl;
	std::cout << "Events with >=" << Npho_cut << " photon     " << cutFlow->GetBinContent(2) << std::endl;
	//std::cout << "Events with >=" << Nele_cut << " electron   " << cutFlow->GetBinContent(3) << std::endl;
	std::cout << "Events with ==" << Nele_cut << " electron   " << cutFlow->GetBinContent(3) << std::endl;
	std::cout << "Events with <= " << NlooseMuVeto_cut << " loose muons " << cutFlow->GetBinContent(4) << std::endl;
	std::cout << "Events with <= " << NlooseEleVeto_cut << " loose electrons " << cutFlow->GetBinContent(5) << std::endl;
	std::cout << "Events with >= " << Njet_cut << " jets        " << cutFlow->GetBinContent(6) << std::endl;
	std::cout << "Events with >= " << NBjet_cut << " b-tag       " << cutFlow->GetBinContent(7) << std::endl;
	std::cout << "Events passing MET cut       " << cutFlow->GetBinContent(8) << std::endl;
	std::cout << "Passing Trigger              " << cutFlow->GetBinContent(9) << std::endl;
	std::cout << std::endl;
}

void Selector::set_cutflow_labels(){
	cutFlow->GetXaxis()->SetBinLabel(1,"Input");
        cutFlow->GetXaxis()->SetBinLabel(2,"Photon");
        cutFlow->GetXaxis()->SetBinLabel(3,"Electron");
        cutFlow->GetXaxis()->SetBinLabel(4,"Loose Mu");
        cutFlow->GetXaxis()->SetBinLabel(5,"Loose Ele");
        cutFlow->GetXaxis()->SetBinLabel(6,"N jets");
        cutFlow->GetXaxis()->SetBinLabel(7,"N b-tags");
        cutFlow->GetXaxis()->SetBinLabel(8,"MET");
        cutFlow->GetXaxis()->SetBinLabel(9,"Trigger");
        //cutFlow->GetXaxis()->SetBinLabel(1,"");

}


void Selector::process_objects(const EventTree* inp_tree){
	tree = inp_tree;
	clear_vectors();
	filter_photons();
	filter_electrons();
	filter_muons();
	filter_jets();
}

bool Selector::fidEtaPass(double Eta){
	double fabsEta = TMath::Abs(Eta);
	if( fabsEta > 2.5) return false;
	if( 1.4442 < fabsEta && fabsEta < 1.566) return false;
	return true;
}

// kEleGammaAndNeutralHadronIso03
// http://cmssw.cvs.cern.ch/cgi-bin/cmssw.cgi/UserCode/EGamma/EGammaAnalysisTools/interface/ElectronEffectiveArea.h
double Selector::eleEffArea03(double SCEta){
	double eta = TMath::Abs(SCEta);
	static const double area[7] = {0.130, 0.137, 0.067, 0.089, 0.107, 0.110, 0.138};
	int region = 0;
	if( eta >= 1.0 )   region++;
	if( eta >= 1.479 ) region++;
	if( eta >= 2.0 )   region++;
	if( eta >= 2.2 )   region++;
	if( eta >= 2.3 )   region++;
	if( eta >= 2.4 )   region++;
	return area[region];
}

// kMuGammaAndNeutralHadronIso04
// http://cmssw.cvs.cern.ch/cgi-bin/cmssw.cgi/UserCode/sixie/Muon/MuonAnalysisTools/interface/MuonEffectiveArea.h?revision=1.7&view=markup
double Selector::muEffArea04(double muEta){
	double eta = TMath::Abs(muEta);
	static const double area[6] = {0.674, 0.565, 0.442, 0.515, 0.821, 0.660};
	int region = 0;
	if( eta >= 1.0 )   region++;
	if( eta >= 1.479 ) region++;
	if( eta >= 2.0 )   region++;
	if( eta >= 2.2 )   region++;
	if( eta >= 2.3 )   region++;
	return area[region];
}

// https://twiki.cern.ch/twiki/bin/view/CMS/CutBasedPhotonID2012#Effective_Areas_for_rho_correcti
int Selector::phoRegion(double absEta){
	int region = 0;
	if( absEta >= 1.0  ) region++;
	if( absEta >= 1.479) region++;
	if( absEta >= 2.0  ) region++;
	if( absEta >= 2.2  ) region++;
	if( absEta >= 2.3  ) region++;
	if( absEta >= 2.4  ) region++;
	return region;
}
double Selector::phoEffArea03ChHad(double phoEta){
	double eta = TMath::Abs(phoEta);
	static const double area[7] = {0.012, 0.010, 0.014, 0.012, 0.016, 0.020, 0.012};
	return area[phoRegion(eta)];
}
double Selector::phoEffArea03NeuHad(double phoEta){
	double eta = TMath::Abs(phoEta);
	static const double area[7] = {0.030, 0.057, 0.039, 0.015, 0.024, 0.039, 0.072};
	return area[phoRegion(eta)];
}
double Selector::phoEffArea03Pho(double phoEta){
	double eta = TMath::Abs(phoEta);
	static const double area[7] = {0.148, 0.130, 0.112, 0.216, 0.262, 0.260, 0.266};
	return area[phoRegion(eta)];
}

Selector::~Selector(){
}
