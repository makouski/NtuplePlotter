#include<vector>
#include<set>
#include<iostream>
#include<TH1F.h>

#include"EventTree.h"
#include"Selector.h"

class EventPick{
public:
	EventPick(const char* titleIn);
	~EventPick();
	void process_event(const EventTree* inp_tree, const Selector* inp_selector);
	void print_cutflow();
	
	std::string title;
	
	// selected object indices
	std::vector<int> Electrons;
	std::vector<int> ElectronsLoose;
	//std::vector<int> Muons;
	std::vector<int> Jets;
	std::vector<int> bJets;
	// indices and selection cuts for photons
	std::vector<int> Photons;
	std::vector<int> PhotonsPresel;
	std::vector<bool> PhoPassChHadIso;
	std::vector<bool> PhoPassPhoIso;
	std::vector<bool> PhoPassSih;
	
	// delta R cuts
	double veto_jet_dR;
	double veto_lep_jet_dR;
	double veto_pho_jet_dR;
	double veto_pho_lep_dR;
	
	// cuts as parameters, to modify easily
	double MET_cut;
	bool no_trigger;
	
	int Njet_ge;
	int NBjet_ge;
	
	int Nele_eq;
	int NlooseEleVeto_le;

	int Npho_ge;
	int NlooseMuVeto_le;

	// variables showing passing or failing selections
	bool passPreSel; // passed preselection
	bool passAll; // single flag: event passed all cuts: preselection + photon
	
	// histograms
	std::vector<TH1F*> histVector;
	TH1F* cutFlow;
	
private:
	const EventTree* tree;
	const Selector* selector;
	
	void clear_vectors();
	void set_cutflow_labels();
	double dR_jet_ele(int jetInd, int eleInd);
	double dR_jet_pho(int jetInd, int phoInd);
	double dR_ele_pho(int eleInd, int phoInd);
}

EventPick::EventPick(const char* titleIn){
	title.assign(titleIn);

	cutFlow = new TH1F("cut_flow","cut flow",10,-0.5,9.5);
	cutFlow->SetDirectory(0);
	set_cutflow_labels(); // keep the labels close to the cuts definitions (below)
	histVector.push_back(cutFlow);
	
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

void EventPick::process_event(const EventTree* inp_tree, const Selector* inp_selector){
	tree = inp_tree;
	selector = inp_selector;
	clear_vectors();

	passPreSel = false;
	passAll = false;
	
	// pre-selection: top ref selection, no photons involved
	// copy jet and electron collections, consiering overlap of jets with electrons and loose electrons:
	// keep jets not close to electrons (veto_jet_dR)
	for(std::vector<int>::iterator jetInd = selector->Jets.begin(); jetInd != selector->Jets.end(); jetInd++){
		bool goodJet = true;
		
		for(std::vector<int>::iterator eleInd = selector->Electrons.begin(); eleInd != selector->Electrons.end(); eleInd++)
			if(dR_jet_ele(*jetInd, *eleInd) < veto_jet_dR) goodJet = false;
		
		for(std::vector<int>::iterator eleInd = selector->ElectronsLoose.begin(); eleInd != selector->ElectronsLoose.end(); eleInd++)
			if(dR_jet_ele(*jetInd, *eleInd) < veto_jet_dR) goodJet = false;
		
		if(goodJet) Jets.push_back(*jetInd);
		
		// take care of bJet collection
		for(std::vector<int>::iterator bjetInd = selector->bJets.begin(); bjetInd != selector->bJets.end(); bjetInd++)
			if(*bjetInd == *jetInd && goodJet) bJets.push_back(*bjetInd);
	}
		
	// keep electrons that are not close to jets (veto_lep_jet_dR)
	for(std::vector<int>::iterator eleInd = selector->Electrons.begin(); eleInd != selector->Electrons.end(); eleInd++){
		bool goodEle = true;
		for(std::vector<int>::iterator jetInd = Jets.begin(); jetInd != Jets.end(); jetInd++)
			if(dR_jet_ele(*jetInd, *eleInd) < veto_lep_jet_dR) goodEle = false;
		
		if(goodEle) Electrons.push_back(*eleInd);
	}
	
	cutFlow->Fill(0); // Input events
	passPreSel = true;
	if(passPreSel && Electrons.size() == Nele_eq) cutFlow->Fill(1);
	else passPreSel = false;
	if(passPreSel && selector->MuonsLoose.size() <= NlooseMuVeto_le) cutFlow->Fill(2);
	else passPreSel = false;
	if(passPreSel && passPreSel && selector->ElectronsLoose.size() <= NlooseEleVeto_le) cutFlow->Fill(3);
	else passPreSel = false;
	if(Jets.size() >= Njet_ge) cutFlow->Fill(4);
	else passPreSel = false;
	if(passPreSel && bJets.size() >= NBjet_ge) cutFlow->Fill(5);
	else passPreSel = false;
	if(passPreSel && tree->pfMET_ > MET_cut) cutFlow->Fill(6);
	else passPreSel = false;
	if(passPreSel && (no_trigger || tree->HLT_[tree->HLTIndex_[17]])) cutFlow->Fill(7);
	else passPreSel = false;
	
	// photon cleaning:
	for(int phoVi = 0; phoVi < selector->PhotonsPresel.size(); phoVi++){
		bool goodPhoton = true;
		// remove photons close (but not too close) to jets
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
	// require >=1 photon
	if(passPreSel && Photons.size() >= Npho_ge){
		cutFlow->Fill(8);
		passAll = true;
	}
}

void EventPick::print_cutflow(){
	std::cout << "Cut-Flow for the event selector: " << title << std::endl;
	std::cout << "Input Events                 " << cutFlow->GetBinContent(1) << std::endl;
	std::cout << "Events with ==" << Nele_eq << " electron   " << cutFlow->GetBinContent(2) << std::endl;
	std::cout << "Events with <= " << NlooseMuVeto_le << " loose muons " << cutFlow->GetBinContent(3) << std::endl;
	std::cout << "Events with <= " << NlooseEleVeto_le << " loose electrons " << cutFlow->GetBinContent(4) << std::endl;
	std::cout << "Events with >= " << Njet_ge << " jets        " << cutFlow->GetBinContent(5) << std::endl;
	std::cout << "Events with >= " << NBjet_ge << " b-tag       " << cutFlow->GetBinContent(6) << std::endl;
	std::cout << "Events passing MET cut       " << cutFlow->GetBinContent(7) << std::endl;
	std::cout << "Passing Trigger              " << cutFlow->GetBinContent(8) << std::endl;
	std::cout << "Events with >=" << Npho_ge << " photon     " << cutFlow->GetBinContent(9) << std::endl;
	std::cout << std::endl;
}

void EventPick::set_cutflow_labels(){
	cutFlow->GetXaxis()->SetBinLabel(1,"Input");
	cutFlow->GetXaxis()->SetBinLabel(2,"Electron");
	cutFlow->GetXaxis()->SetBinLabel(3,"Loose Mu");
	cutFlow->GetXaxis()->SetBinLabel(4,"Loose Ele");
	cutFlow->GetXaxis()->SetBinLabel(5,"N jets");
	cutFlow->GetXaxis()->SetBinLabel(6,"N b-tags");
	cutFlow->GetXaxis()->SetBinLabel(7,"MET");
	cutFlow->GetXaxis()->SetBinLabel(8,"Trigger");
	cutFlow->GetXaxis()->SetBinLabel(9,"Photon");
	//cutFlow->GetXaxis()->SetBinLabel(1,"");
}

void EventPick::clear_vectors(){
	Electrons.clear();
	ElectronsLoose.clear();
	//Muons.clear();
	Jets.clear();
	bJets.clear();
	Photons.clear();
	PhotonsPresel.clear();
	PhoPassChHadIso.clear();
	PhoPassPhoIso.clear();
	PhoPassSih.clear();
}

double EventPick::dR_jet_ele(int jetInd, int eleInd){
	return dR(tree->jetEta_[jetInd], tree->jetPhi_[jetInd], tree->eleSCEta_[eleInd], tree->elePhi_[eleInd]);
}
double EventPick::dR_jet_pho(int jetInd, int phoInd){
	return dR(tree->jetEta_[jetInd], tree->jetPhi_[jetInd], tree->phoEta_[phoInd], tree->phoPhi_[phoInd]);
}
double EventPick::dR_ele_pho(int eleInd, int phoInd){
	return dR(tree->eleSCEta_[eleInd], tree->elePhi_[eleInd], tree->phoEta_[phoInd], tree->phoPhi_[phoInd]);
}
