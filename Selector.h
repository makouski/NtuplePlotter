#ifndef SELECTOR_H
#define SELECTOR_H

#include<vector>
#include<set>
#include<iostream>
#include<algorithm>
#include<TH1F.h>
#include<TMath.h>
#include<TLorentzVector.h>
#include"EventTree.h"

// https://twiki.cern.ch/twiki/bin/viewauth/CMS/CutBasedPhotonID2012
// photon ID is not going to be changed every time this code runs
// barrel/endcap, Loose/Medium/Tight
const int    photonID_IsConv[2][3]                = { {0, 0, 0} ,             {0, 0, 0}             };
const double photonID_HoverE[2][3]                = { {0.05, 0.05, 0.05} ,    {0.05, 0.05, 0.05}    };
const double photonID_SigmaIEtaIEta[2][3]         = { {0.012, 0.011, 0.011} , {0.034, 0.033, 0.031} };
const double photonID_RhoCorrR03ChHadIso[2][3]    = { {2.6, 1.5, 0.7} ,       {2.3, 1.2, 0.5}       };
const double photonID_RhoCorrR03NeuHadIso_0[2][3] = { {3.5, 1.0, 0.4} ,       {2.9, 1.5, 1.5}       };
const double photonID_RhoCorrR03NeuHadIso_1[2][3] = { {0.04, 0.04, 0.04} ,    {0.04, 0.04, 0.04}    };
const double photonID_RhoCorrR03PhoIso_0[2][3]    = { {1.3, 0.7, 0.5} ,       {999, 1.0, 1.0}       };
const double photonID_RhoCorrR03PhoIso_1[2][3]    = { {0.005, 0.005, 0.005} , {0.005, 0.005, 0.005} };

double dR(double eta1, double phi1, double eta2, double phi2);

class Selector{
public:
	Selector(const char* titleIn);
	~Selector();
	void process_objects(const EventTree* inp_tree);
	int  accept_counts();
	void print_cutflow();
	
	std::string title;
	
	// selected object indices
	std::vector<int> Photons;
	std::set<int>    VetoedPhotons, WVetoedPhotons;
	std::vector<int> Electrons, ElectronsLoose;
	std::set<int>    VetoedElectrons;
	std::vector<int> Muons, MuonsLoose;
	std::vector<int> JetCands, JetCandsNoDR, Jets;
	std::set<int>    VetoedJets;
	
	// calculated rho corrected PF isolations
	std::vector<double> Ele03RelIso;
	std::vector<double> Mu04RelIso;
	std::vector<double> Pho03ChHadIso, Pho03ChHadSCRIso, Pho03NeuHadIso, Pho03PhoIso, Pho03PhoSCRIso;	
	
	// histograms
	TH1F* Nele_before_cut;
	std::vector<TH1F*> histVector;
	// some inv mass variables, for cutting or histogramminng
	double JJGammaMassMin;
	double JJEMassMin;
	double bJJGammaMassMin;
	double bJJEMassMin;
        
	double EGammaMETtransMassMin;

	// cut-flow
	TH1F* cutFlow;
	
	// cuts as parameters, to modify easily
	double MET_cut;
	bool no_trigger;
	
	// jets
	std::vector<double> jet_Pt_cuts;
	double btag_cut;
	int Njet_cut;
	int NBjet_cut;

	// electrons
	double ele_Pt_cut;
	double ele_PtLoose_cut;
	double ele_RelIso_range[2];
	double ele_RelIsoLoose_cut;
	double ele_MVA_range[2];
	double ele_MVALoose_cut;
	double ele_Dxy_cut;
	int    ele_MissInnHit_cut;
	int    Nele_cut, NlooseEleVeto_cut;
	
	// photons
	double pho_Et_cut;
	int    pho_ID_ind; // 0 - Loose, 1 - Medium, 2 - Tight
	bool   pho_noSigmaIEta_cut;
	bool   pho_noIso_cut;
	bool   pho_noChHadIso_cut;
	bool   pho_noPhoIso_cut;
	bool   pho_noPixelSeed_cut;
	bool   pho_noEleVeto_cut;
	int    Npho_cut;
	// muons
	double mu_PtLoose_cut;
	double mu_RelIsoLoose_cut;
	int    NlooseMuVeto_cut;

	// delta R cuts
	bool   doVeto;
	double veto_jet_dR;
	double veto_lep_jet_dR;
	double veto_pho_jet_dR;
	double veto_pho_lep_dR;
	double W_mass_veto;
	double W_trans_mass_veto;
	
private:
	const EventTree* tree;
	
	void clear_vectors();
	void set_cutflow_labels();
	void filter_photons();
	void filter_electrons();
	void filter_muons();
	void filter_jets();
	void make_dR_cuts();
	void make_mass_cuts();
	void vetoCollection(std::vector<int>* toClean, const Float_t* etaClean, const Float_t* phiClean, 
						const std::vector<int>* toCheck, const Float_t* etaCheck, const Float_t* phiCheck, 
						std::set<int>* vetoed, double dRcut);
	void cleanCollection(std::vector<int>* toClean, std::set<int>* vetoed);
	double invMass_jjp(int jetInd1, int jetInd2, int phoInd);
	double invMass_jje(int jetInd1, int jetInd2, int eleInd);
	double invMass_jj(int jetInd1, int jetInd2);
	double clusterTransMass_pe(int phoInd, int eleInd);
	
	double JetPtCut(int jetInd);
	bool fidEtaPass(double Eta);
	
	// effective areas, see Selector.cpp for more information
	double eleEffArea03(double SCEta);
	double muEffArea04(double muEta);
	double phoEffArea03ChHad(double phoEta);
	double phoEffArea03NeuHad(double phoEta);
	double phoEffArea03Pho(double phoEta);
	int phoRegion(double absEta);
	
};

#endif
