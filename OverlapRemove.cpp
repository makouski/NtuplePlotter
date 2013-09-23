#include"EventTree.h"
#include<iostream>
double dR(double eta1, double phi1, double eta2, double phi2);

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
