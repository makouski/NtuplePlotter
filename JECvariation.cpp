#include"JetMETObjects/JetCorrectorParameters.h"
#include"JetMETObjects/FactorizedJetCorrector.h"

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


void applyJEC(EventTree* tree, FactorizedJetCorrector* JetCorrector){
	for(int jetInd = 0; jetInd < tree->nJet_ ; ++jetInd){
		JetCorrector->setJetEta(tree->jetEta_[jetInd]);
		JetCorrector->setJetPt(tree->jetRawPt_[jetInd]);
		JetCorrector->setJetA(tree->jetArea_[jetInd]);
		JetCorrector->setRho(tree->rho2012_);
		
		double correction = JetCorrector->getCorrection();
		//std::cout << tree->jetPt_[jetInd] << "  " << tree->jetRawPt_[jetInd] << "  " << tree->jetPt_[jetInd]/tree->jetRawPt_[jetInd] << "  " << correction << std::endl;
		tree->jetPt_[jetInd] = tree->jetRawPt_[jetInd] * correction;
		tree->jetEn_[jetInd] = tree->jetRawEn_[jetInd] * correction;
	}
	
}
