#include"EventPick.h"

double secondMinDr(int myInd, std::vector<float> *etas, std::vector<float> *phis);

EventPick::EventPick(std::string titleIn){
	title = titleIn;

	cutFlow = new TH1F("cut_flow","cut flow",10,-0.5,9.5);
	cutFlow->SetDirectory(0);
	set_cutflow_labels(cutFlow); // keep the labels close to the cuts definitions (below)
	histVector.push_back(cutFlow);
	
	cutFlowWeight = new TH1F("cut_flow_weight","cut flow with PU weight",10,-0.5,9.5);
	cutFlowWeight->SetDirectory(0);
	set_cutflow_labels(cutFlowWeight);
	histVector.push_back(cutFlowWeight);

	genPhoRegionWeight = new TH1F("genPhoRegionWeight","GenPhoton passing fiducial cuts: barrel 0 or endcap 1",2,-0.5,1.5);
	genPhoRegionWeight->SetDirectory(0);
	histVector.push_back(genPhoRegionWeight);

	// assign cut values
	veto_jet_dR = 0.3;
	veto_lep_jet_dR = 0.5;
	veto_pho_jet_dR = 0.7;
	veto_pho_lep_dR = 0.7;
	MET_cut = 20.0;
	no_trigger = false;
	Njet_ge = 3;
	NBjet_ge = 1;
	Nele_eq = 1;
	NlooseEleVeto_le = 99; // no cut
	Npho_ge = 1;
	NlooseMuVeto_le = 0;
}

EventPick::~EventPick(){
}

void EventPick::process_event(const EventTree* inp_tree, const Selector* inp_selector, double weight){
	tree = inp_tree;
	selector = inp_selector;
	clear_vectors();

	passPreSel = false;
	passAll = false;
	
	// pre-selection: top ref selection
	// copy jet and electron collections, consiering overlap of jets with electrons, loose electrons, photons:
	// keep jets not close to electrons and photons (veto_jet_dR)
	for(std::vector<int>::const_iterator jetInd = selector->Jets.begin(); jetInd != selector->Jets.end(); jetInd++){
		bool goodJet = true;
		
		for(std::vector<int>::const_iterator eleInd = selector->Electrons.begin(); eleInd != selector->Electrons.end(); eleInd++)
			if(dR_jet_ele(*jetInd, *eleInd) < veto_jet_dR) goodJet = false;
		
		for(std::vector<int>::const_iterator eleInd = selector->ElectronsLoose.begin(); eleInd != selector->ElectronsLoose.end(); eleInd++)
			if(dR_jet_ele(*jetInd, *eleInd) < veto_jet_dR) goodJet = false;
		
		//for(int phoVi = 0; phoVi < selector->PhotonsPresel.size(); phoVi++){
		//	if(selector->PhoPassChHadIso.at(phoVi) && 
		//	selector->PhoPassPhoIso.at(phoVi) &&
		//	selector->PhoPassSih.at(phoVi) &&
		//	dR_jet_pho(*jetInd, selector->PhotonsPresel.at(phoVi)) < 0.1)
		//		goodJet = false;
		//}
 
		if(goodJet) Jets.push_back(*jetInd);
		
		// take care of bJet collection
		for(std::vector<int>::const_iterator bjetInd = selector->bJets.begin(); bjetInd != selector->bJets.end(); bjetInd++)
			if(*bjetInd == *jetInd && goodJet) bJets.push_back(*bjetInd);
	}
	
	// keep electrons that are not close to jets (veto_lep_jet_dR)
	for(std::vector<int>::const_iterator eleInd = selector->Electrons.begin(); eleInd != selector->Electrons.end(); eleInd++){
		bool goodEle = true;
		for(std::vector<int>::iterator jetInd = Jets.begin(); jetInd != Jets.end(); jetInd++)
			if(dR_jet_ele(*jetInd, *eleInd) < veto_lep_jet_dR) goodEle = false;
		
		if(goodEle) Electrons.push_back(*eleInd);
	}
	
	// loose electrons
	for(std::vector<int>::const_iterator eleInd = selector->ElectronsLoose.begin(); eleInd != selector->ElectronsLoose.end(); eleInd++){
		bool goodEle = true;
		for(std::vector<int>::iterator jetInd = Jets.begin(); jetInd != Jets.end(); jetInd++)
			if(dR_jet_ele(*jetInd, *eleInd) < veto_lep_jet_dR) goodEle = false;
		
		if(goodEle) ElectronsLoose.push_back(*eleInd);
	}
	
	// photon cleaning:
	for(int phoVi = 0; phoVi < selector->PhotonsPresel.size(); phoVi++){
		bool goodPhoton = true;
		// remove photons close to jets
		for(std::vector<int>::iterator jetInd = Jets.begin(); jetInd != Jets.end(); jetInd++){		
			double drjp = dR_jet_pho(*jetInd, selector->PhotonsPresel.at(phoVi));
			if(veto_jet_dR <= drjp && drjp < veto_pho_jet_dR) goodPhoton = false;
		}
		// and electrons
		for(std::vector<int>::iterator eleInd = Electrons.begin(); eleInd != Electrons.end(); eleInd++)
			if(dR_ele_pho(*eleInd, selector->PhotonsPresel.at(phoVi)) < veto_pho_lep_dR) goodPhoton = false;
		
		if(goodPhoton){
			PhotonsPresel.push_back(selector->PhotonsPresel.at(phoVi));
			PhoPassChHadIso.push_back(selector->PhoPassChHadIso.at(phoVi));
			PhoPassPhoIso.push_back(selector->PhoPassPhoIso.at(phoVi));
			PhoPassSih.push_back(selector->PhoPassSih.at(phoVi));
			if(PhoPassChHadIso.back() && PhoPassPhoIso.back() && PhoPassSih.back()) 
				Photons.push_back(selector->PhotonsPresel.at(phoVi));
		}
	}
	
	cutFlow->Fill(0.0); // Input events
	cutFlowWeight->Fill(0.0,weight);
	passPreSel = true;
	if(passPreSel && tree->IsVtxGood_>=0 && (no_trigger || tree->HLT_[tree->HLTIndex_[17]])) {cutFlow->Fill(1); cutFlowWeight->Fill(1,weight);}
	else passPreSel = false;
	if(passPreSel && Electrons.size() == Nele_eq) {cutFlow->Fill(2); cutFlowWeight->Fill(2,weight);}
	else passPreSel = false;
	if(passPreSel && selector->MuonsLoose.size() <= NlooseMuVeto_le) {cutFlow->Fill(3); cutFlowWeight->Fill(3,weight);}
	else passPreSel = false;
	if(passPreSel && selector->ElectronsLoose.size() <= NlooseEleVeto_le) {cutFlow->Fill(4); cutFlowWeight->Fill(4,weight);}
	else passPreSel = false;
	if(passPreSel && Jets.size() >= Njet_ge) {cutFlow->Fill(5); cutFlowWeight->Fill(5,weight);}
	else passPreSel = false;
	if(passPreSel && bJets.size() >= NBjet_ge) {cutFlow->Fill(6); cutFlowWeight->Fill(6,weight);}
	else passPreSel = false;
	if(passPreSel && tree->pfMET_ > MET_cut) {cutFlow->Fill(7); cutFlowWeight->Fill(7,weight);}
	else passPreSel = false;
	
	// require >=1 photon
	if(passPreSel && Photons.size() >= Npho_ge){
		cutFlow->Fill(8);
		cutFlowWeight->Fill(8,weight);
		passAll = true;
	}

	// saving information about Gen Level photons, if any
	// Save it only if PreSelection passed
	// Separate count for barrel and endcap (will be used separately anyway)
	
	bool foundGenPhotonBarrel = false;
	bool foundGenPhotonEndcap = false;
	if(passPreSel && !(tree->isData_)){
		for(int mcInd=0; mcInd<tree->nMC_; ++mcInd){
			if(tree->mcPID->at(mcInd) == 22 &&
			   (tree->mcParentage->at(mcInd)==2 || tree->mcParentage->at(mcInd)==10) &&
			   tree->mcPt->at(mcInd) > selector->pho_Et_cut){		
				if(secondMinDr(mcInd, tree->mcEta, tree->mcPhi) > 0.05){
					double fabsEta = TMath::Abs(tree->mcEta->at(mcInd));
					if(fabsEta < 1.4442) foundGenPhotonBarrel = true;
					if( 1.566 < fabsEta && fabsEta < 2.5) foundGenPhotonEndcap = true;
				}
			}
		}
	}
	if(foundGenPhotonBarrel) genPhoRegionWeight->Fill(0.0, weight);
	if(foundGenPhotonEndcap) genPhoRegionWeight->Fill(1.0, weight);
}

void EventPick::print_cutflow(){
	std::cout << "Cut-Flow for the event selector: " << title << std::endl;
	std::cout << "Input Events                 " << cutFlow->GetBinContent(1) << std::endl;
	std::cout << "Passing Trigger and PV       " << cutFlow->GetBinContent(2) << std::endl;
	std::cout << "Events with ==" << Nele_eq << " electron     " << cutFlow->GetBinContent(3) << std::endl;
	std::cout << "Events with <= " << NlooseMuVeto_le << " loose muons " << cutFlow->GetBinContent(4) << std::endl;
	std::cout << "Events with <= " << NlooseEleVeto_le << " loose electrons " << cutFlow->GetBinContent(5) << std::endl;
	std::cout << "Events with >= " << Njet_ge << " jets        " << cutFlow->GetBinContent(6) << std::endl;
	std::cout << "Events with >= " << NBjet_ge << " b-tag       " << cutFlow->GetBinContent(7) << std::endl;
	std::cout << "Events passing MET cut       " << cutFlow->GetBinContent(8) << std::endl;
	std::cout << "Events with >=" << Npho_ge << " photon       " << cutFlow->GetBinContent(9) << std::endl;
	std::cout << std::endl;
}

void EventPick::set_cutflow_labels(TH1F* hist){
	hist->GetXaxis()->SetBinLabel(1,"Input");
	hist->GetXaxis()->SetBinLabel(2,"Trigger");
	hist->GetXaxis()->SetBinLabel(3,"Electron");
	hist->GetXaxis()->SetBinLabel(4,"Loose Mu");
	hist->GetXaxis()->SetBinLabel(5,"Loose Ele");
	hist->GetXaxis()->SetBinLabel(6,"N jets");
	hist->GetXaxis()->SetBinLabel(7,"N b-tags");
	hist->GetXaxis()->SetBinLabel(8,"MET");
	hist->GetXaxis()->SetBinLabel(9,"Photon");
	//hist->GetXaxis()->SetBinLabel(1,"");
}

void EventPick::clear_vectors(){
	Electrons.clear();
	ElectronsLoose.clear();
	//Muons.clear();
	MuonsLoose.clear();
	Jets.clear();
	bJets.clear();
	Photons.clear();
	PhotonsPresel.clear();
	PhoPassChHadIso.clear();
	PhoPassPhoIso.clear();
	PhoPassSih.clear();
}

double EventPick::dR_jet_ele(int jetInd, int eleInd){
	return dR(tree->jetEta_->at(jetInd), tree->jetPhi_->at(jetInd), tree->eleSCEta_->at(eleInd), tree->elePhi_->at(eleInd));
}
double EventPick::dR_jet_pho(int jetInd, int phoInd){
	return dR(tree->jetEta_->at(jetInd), tree->jetPhi_->at(jetInd), tree->phoEta_->at(phoInd), tree->phoPhi_->at(phoInd));
}
double EventPick::dR_ele_pho(int eleInd, int phoInd){
	return dR(tree->eleSCEta_->at(eleInd), tree->elePhi_->at(eleInd), tree->phoEta_->at(phoInd), tree->phoPhi_->at(phoInd));
}
