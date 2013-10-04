#include"EventTree.h"

EventTree::EventTree(int nFiles, char** fileNames){
	chain = new TChain("ggNtuplizer/EventTree");
	for(int fileI=0; fileI<nFiles; fileI++){
		chain->Add(fileNames[fileI]);
	}
	chain->SetBranchStatus("*",0);
	
	// keep some important branches
	chain->SetBranchStatus("nHLT",1);
	chain->SetBranchAddress("nHLT", &nHLT_);
	chain->SetBranchStatus("HLT",1);
	chain->SetBranchAddress("HLT", HLT_);
	chain->SetBranchStatus("HLTIndex",1);
	chain->SetBranchAddress("HLTIndex", HLTIndex_);
	chain->SetBranchStatus("bspotPos",1);
	chain->SetBranchAddress("bspotPos", bspotPos_);
	chain->SetBranchStatus("IsVtxGood",1);
	chain->SetBranchAddress("IsVtxGood", &IsVtxGood_);
	chain->SetBranchStatus("nPUInfo",1);
	chain->SetBranchAddress("nPUInfo", &nPUInfo_);
	chain->SetBranchStatus("nPU",1);
	chain->SetBranchAddress("nPU", &nPU_);
	chain->SetBranchStatus("puBX",1);
	chain->SetBranchAddress("puBX", &puBX_);
	chain->SetBranchStatus("puTrue",1);
	chain->SetBranchAddress("puTrue", &puTrue_);
	//chain->SetBranchStatus("",1);
	
	// event
	
	chain->SetBranchStatus("run",1);
	chain->SetBranchAddress("run", &run_);

	chain->SetBranchStatus("event",1);
	chain->SetBranchAddress("event", &event_);
	
	chain->SetBranchStatus("lumis",1);
	chain->SetBranchAddress("lumis", &lumis_);

	chain->SetBranchStatus("isData",1);
	chain->SetBranchAddress("isData", &isData_);

	chain->SetBranchStatus("nVtx",1);
	chain->SetBranchAddress("nVtx", &nVtx_);

	chain->SetBranchStatus("pfMET",1);
	chain->SetBranchAddress("pfMET", &pfMET_);

	chain->SetBranchStatus("pfMETPhi",1);
	chain->SetBranchAddress("pfMETPhi", &pfMETPhi_);
	
	// electrons	
	
	chain->SetBranchStatus("nEle",1);
	chain->SetBranchAddress("nEle", &nEle_);

	//chain->SetBranchStatus("eleEta",1);

	chain->SetBranchStatus("elePt",1);
	chain->SetBranchAddress("elePt", &elePt_);

	chain->SetBranchStatus("eleSCEta",1);
	chain->SetBranchAddress("eleSCEta", &eleSCEta_);

	chain->SetBranchStatus("elePhi",1);
	chain->SetBranchAddress("elePhi", &elePhi_);

	chain->SetBranchStatus("elePFChIso03",1);
	chain->SetBranchAddress("elePFChIso03", &elePFChIso03_);

	chain->SetBranchStatus("elePFNeuIso03",1);
	chain->SetBranchAddress("elePFNeuIso03", &elePFNeuIso03_);

	chain->SetBranchStatus("elePFPhoIso03",1);
	chain->SetBranchAddress("elePFPhoIso03", &elePFPhoIso03_);

	chain->SetBranchStatus("rho2012",1);
	chain->SetBranchAddress("rho2012", &rho2012_);

	chain->SetBranchStatus("eleIDMVATrig",1);
	chain->SetBranchAddress("eleIDMVATrig", &eleIDMVATrig_);

	chain->SetBranchStatus("eleD0",1);
	chain->SetBranchAddress("eleD0", &eleD0_);

	chain->SetBranchStatus("eleMissHits",1);
	chain->SetBranchAddress("eleMissHits", &eleMissHits_);

	chain->SetBranchStatus("eleConvVtxFit",1);
	chain->SetBranchAddress("eleConvVtxFit", &eleConvVtxFit_);

	chain->SetBranchStatus("eleDz",1);
	chain->SetBranchAddress("eleDz", &eleDz_);

	chain->SetBranchStatus("eleEoverP",1);
	chain->SetBranchAddress("eleEoverP", &eleEoverP_);

	chain->SetBranchStatus("eleSigmaIEtaIEta",1);
	chain->SetBranchAddress("eleSigmaIEtaIEta", &eleSigmaIEtaIEta_);

	chain->SetBranchStatus("eledEtaAtVtx",1);
	chain->SetBranchAddress("eledEtaAtVtx", &eledEtaAtVtx_);

	chain->SetBranchStatus("eledPhiAtVtx",1);
	chain->SetBranchAddress("eledPhiAtVtx", &eledPhiAtVtx_);

	chain->SetBranchStatus("eleEcalEn",1);
	chain->SetBranchAddress("eleEcalEn", &eleEcalEn_);

	chain->SetBranchStatus("eleHoverE",1);
	chain->SetBranchAddress("eleHoverE", &eleHoverE_);

	chain->SetBranchStatus("eleIsoTrkDR03",1);
	chain->SetBranchAddress("eleIsoTrkDR03", &eleIsoTrkDR03_);

	chain->SetBranchStatus("eleIsoEcalDR03",1);
	chain->SetBranchAddress("eleIsoEcalDR03", &eleIsoEcalDR03_);

	chain->SetBranchStatus("eleIsoHcalDR03",1);
	chain->SetBranchAddress("eleIsoHcalDR03", &eleIsoHcalDR03_);

	chain->SetBranchStatus("eleGenIndex",1);
	chain->SetBranchAddress("eleGenIndex", &eleGenIndex_);

	chain->SetBranchStatus("eleGenGMomPID",1);
	chain->SetBranchAddress("eleGenGMomPID", &eleGenGMomPID_);

	chain->SetBranchStatus("eleGenMomPID",1);
	chain->SetBranchAddress("eleGenMomPID", &eleGenMomPID_);

	// muons
	
	chain->SetBranchStatus("nMu",1);
	chain->SetBranchAddress("nMu", &nMu_);

	chain->SetBranchStatus("muPt",1);
	chain->SetBranchAddress("muPt", &muPt_);

	chain->SetBranchStatus("muEta",1);
	chain->SetBranchAddress("muEta", &muEta_);

	chain->SetBranchStatus("muPhi",1);
	chain->SetBranchAddress("muPhi", &muPhi_);
	
	chain->SetBranchStatus("muPFIsoR04_CH",1);
	chain->SetBranchAddress("muPFIsoR04_CH", &muPFIsoR04_CH_);
	
	chain->SetBranchStatus("muPFIsoR04_NH",1);
	chain->SetBranchAddress("muPFIsoR04_NH", &muPFIsoR04_NH_);
	
	chain->SetBranchStatus("muPFIsoR04_Pho",1);
	chain->SetBranchAddress("muPFIsoR04_Pho", &muPFIsoR04_Pho_);
	
	//chain->SetBranchStatus("rho25_muPFiso",1);
	//chain->SetBranchAddress("rho25_muPFiso", &rho25_muPFiso_);
	
	// jets
	
	chain->SetBranchStatus("nJet",1);
	chain->SetBranchAddress("nJet", &nJet_);

	chain->SetBranchStatus("jetPt",1);
	chain->SetBranchAddress("jetPt", &jetPt_);

	chain->SetBranchStatus("jetRawPt",1);
        chain->SetBranchAddress("jetRawPt", &jetRawPt_);
		
	chain->SetBranchStatus("jetEta",1);
	chain->SetBranchAddress("jetEta", &jetEta_);
	
	chain->SetBranchStatus("jetPhi",1);
	chain->SetBranchAddress("jetPhi", &jetPhi_);

	chain->SetBranchStatus("jetEn",1);
	chain->SetBranchAddress("jetEn", &jetEn_);

	chain->SetBranchStatus("jetArea",1);
	chain->SetBranchAddress("jetArea", &jetArea_);

	chain->SetBranchStatus("jetNConstituents",1);
	chain->SetBranchAddress("jetNConstituents", &jetNConstituents_);
	
	chain->SetBranchStatus("jetNCharged",1);
	chain->SetBranchAddress("jetNCharged", &jetNCharged_);
	
	chain->SetBranchStatus("jetCEF",1);
	chain->SetBranchAddress("jetCEF", &jetCEF_);

	chain->SetBranchStatus("jetNHF",1);
	chain->SetBranchAddress("jetNHF", &jetNHF_);
	
	chain->SetBranchStatus("jetNEF",1);
	chain->SetBranchAddress("jetNEF", &jetNEF_);
	
	chain->SetBranchStatus("jetCHF",1);
	chain->SetBranchAddress("jetCHF", &jetCHF_);

	chain->SetBranchStatus("jetCombinedSecondaryVtxBJetTags",1);
	chain->SetBranchAddress("jetCombinedSecondaryVtxBJetTags", &jetCombinedSecondaryVtxBJetTags_);
	
	chain->SetBranchStatus("jetCombinedSecondaryVtxMVABJetTags",1);
	chain->SetBranchAddress("jetCombinedSecondaryVtxMVABJetTags", &jetCombinedSecondaryVtxMVABJetTags_);

	//chain->SetBranchStatus("jetJetProbabilityBJetTags",1);
	//chain->SetBranchAddress("jetJetProbabilityBJetTags", &jetJetProbabilityBJetTags_);
	
	//chain->SetBranchStatus("jetJetBProbabilityBJetTags",1);
	//chain->SetBranchAddress("jetJetBProbabilityBJetTags", &jetJetBProbabilityBJetTags_);
	
	//chain->SetBranchStatus("jetTrackCountingHighPurBJetTags",1);
	//chain->SetBranchAddress("jetTrackCountingHighPurBJetTags", &jetTrackCountingHighPurBJetTags_);
	
	chain->SetBranchStatus("jetPartonID",1);
	chain->SetBranchAddress("jetPartonID", &jetPartonID_);
	
	chain->SetBranchStatus("jetGenPartonID",1);
	chain->SetBranchAddress("jetGenPartonID", &jetGenPartonID_);
	
	chain->SetBranchStatus("jetGenJetIndex",1);
	chain->SetBranchAddress("jetGenJetIndex", &jetGenJetIndex_);

	chain->SetBranchStatus("jetGenJetPt",1);
	chain->SetBranchAddress("jetGenJetPt", &jetGenJetPt_);

	chain->SetBranchStatus("jetGenPt",1);
	chain->SetBranchAddress("jetGenPt", &jetGenPt_);
	
	chain->SetBranchStatus("jetGenEta",1);
	chain->SetBranchAddress("jetGenEta", &jetGenEta_);
	
	chain->SetBranchStatus("jetGenPhi",1);
	chain->SetBranchAddress("jetGenPhi", &jetGenPhi_);

	// photons
	
	chain->SetBranchStatus("nPho",1);
	chain->SetBranchAddress("nPho", &nPho_);
	
	chain->SetBranchStatus("phoEt",1);
	chain->SetBranchAddress("phoEt", &phoEt_);
	
	chain->SetBranchStatus("phoEta",1);
	chain->SetBranchAddress("phoEta", &phoEta_);

	chain->SetBranchStatus("phoPhi",1);
	chain->SetBranchAddress("phoPhi", &phoPhi_);
	
	chain->SetBranchStatus("phoIsConv",1);
	chain->SetBranchAddress("phoIsConv", &phoIsConv_);
	
	chain->SetBranchStatus("phohasPixelSeed",1);
	chain->SetBranchAddress("phohasPixelSeed", &phohasPixelSeed_);

	chain->SetBranchStatus("phoEleVeto",1);
	chain->SetBranchAddress("phoEleVeto", &phoEleVeto_);

	chain->SetBranchStatus("phoHoverE",1);
	chain->SetBranchAddress("phoHoverE", &phoHoverE_);

	chain->SetBranchStatus("phoSigmaIEtaIEta",1);
	chain->SetBranchAddress("phoSigmaIEtaIEta", &phoSigmaIEtaIEta_);
	
	chain->SetBranchStatus("phoPFChIso",1);
	chain->SetBranchAddress("phoPFChIso", &phoPFChIso_);
	
	chain->SetBranchStatus("phoPFNeuIso",1);
	chain->SetBranchAddress("phoPFNeuIso", &phoPFNeuIso_);

	chain->SetBranchStatus("phoPFPhoIso",1);
	chain->SetBranchAddress("phoPFPhoIso", &phoPFPhoIso_);

	chain->SetBranchStatus("phoSCRPhoIso",1);
	chain->SetBranchAddress("phoSCRPhoIso", &phoSCRPhoIso_);

	chain->SetBranchStatus("phoSCRChIso",1);
	chain->SetBranchAddress("phoSCRChIso", &phoSCRChIso_);
	
	chain->SetBranchStatus("phoRandConePhoIso",1);
	chain->SetBranchAddress("phoRandConePhoIso", &phoRandConePhoIso_);
	
	chain->SetBranchStatus("phoRandConeChIso",1);
	chain->SetBranchAddress("phoRandConeChIso", &phoRandConeChIso_);

	//chain->SetBranchStatus("rho25",1);
	//chain->SetBranchAddress("rho25", &rho25_);
	
	//chain->SetBranchStatus("rho25_neu",1);
	//chain->SetBranchAddress("rho25_neu", &rho25_neu_);

	chain->SetBranchStatus("phoGenIndex",1);
	chain->SetBranchAddress("phoGenIndex", &phoGenIndex_);
	
	chain->SetBranchStatus("phoGenGMomPID",1);
	chain->SetBranchAddress("phoGenGMomPID", &phoGenGMomPID_);
	
	chain->SetBranchStatus("phoGenMomPID",1);
	chain->SetBranchAddress("phoGenMomPID", &phoGenMomPID_);
	
	// MC gen particles
	
	chain->SetBranchStatus("nMC",1);
	chain->SetBranchAddress("nMC", &nMC_);
	
	chain->SetBranchStatus("mcPt",1);
	chain->SetBranchAddress("mcPt", &mcPt);
	
	chain->SetBranchStatus("mcEta",1);
	chain->SetBranchAddress("mcEta", &mcEta);
	
	chain->SetBranchStatus("mcPhi",1);
	chain->SetBranchAddress("mcPhi", &mcPhi);
	
	chain->SetBranchStatus("mcMass",1);
	chain->SetBranchAddress("mcMass", &mcMass);
	
	chain->SetBranchStatus("mcPID",1);
	chain->SetBranchAddress("mcPID", &mcPID);
	
	chain->SetBranchStatus("mcMomPID",1);
	chain->SetBranchAddress("mcMomPID", &mcMomPID);
	
	chain->SetBranchStatus("mcGMomPID",1);
	chain->SetBranchAddress("mcGMomPID", &mcGMomPID);
	
	chain->SetBranchStatus("mcDecayType",1);
	chain->SetBranchAddress("mcDecayType", &mcDecayType);
	
	chain->SetBranchStatus("mcIndex",1);
	chain->SetBranchAddress("mcIndex", &mcIndex);
	
	chain->SetBranchStatus("mcMomPt",1);
	chain->SetBranchAddress("mcMomPt", &mcMomPt);
	
	chain->SetBranchStatus("mcMomEta",1);
	chain->SetBranchAddress("mcMomEta", &mcMomEta);
	
	chain->SetBranchStatus("mcMomPhi",1);
	chain->SetBranchAddress("mcMomPhi", &mcMomPhi);
	
	chain->SetBranchStatus("mcMomMass",1);
	chain->SetBranchAddress("mcMomMass", &mcMomMass);
	
	chain->SetBranchStatus("mcParentage",1);
	chain->SetBranchAddress("mcParentage", &mcParentage);
	
	//chain->SetBranchStatus("",1);
	//chain->SetBranchAddress("", _);
}

EventTree::~EventTree(){
	delete chain;
}

Long64_t EventTree::GetEntries(){
	return chain->GetEntries();
}

Int_t EventTree::GetEntry(Long64_t entry){
	chain->GetEntry(entry);
}
