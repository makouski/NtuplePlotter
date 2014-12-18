#include"EventTree.h"
#include<iostream>
#include<cstdlib>
#include <math.h>

double dR(double eta1, double phi1, double eta2, double phi2);

bool overlapWHIZARD(EventTree* tree){
	const double Et_cut = 20;
	const double dR_cut = 0.1;
	// consider only mcParticles here
	bool haveOverlap = false;
	for(int phoInd=0; phoInd<tree->nMC_; ++phoInd){
		// find photons with Et above cut and that came from ISR or top
		if(tree->mcPID->at(phoInd)==22 && tree->mcPt->at(phoInd) > Et_cut && (tree->mcParentage->at(phoInd)==2 || (tree->mcParentage->at(phoInd)==10 && abs(tree->mcMomPID->at(phoInd))==24))){
			bool closeToLeg = false;
			bool haveLeg = false;
			// find "legs" (b-qark coming from top)
			for(int legInd=0; legInd<tree->nMC_; ++legInd){
				if(abs(tree->mcPID->at(legInd))==5 && abs(tree->mcMomPID->at(legInd))==6){
					haveLeg=true;
					if(dR(tree->mcEta->at(phoInd),tree->mcPhi->at(phoInd),tree->mcEta->at(legInd),tree->mcPhi->at(legInd)) < dR_cut) closeToLeg=true;
				}
			}
			if(haveLeg && !closeToLeg) haveOverlap = true;
		}
	}
	if (0 && !haveOverlap){
		std::cout << "no overlap" << std::endl;
		for(int legInd=0; legInd<tree->nMC_; ++legInd)
			std::cout << "PID " << tree->mcPID->at(legInd) 
			<< "  Pt " << tree->mcPt->at(legInd) 
			<< "  Eta " << tree->mcEta->at(legInd) 
			<< "  Phi " << tree->mcPhi->at(legInd) << std::endl;
	}
	return haveOverlap;
}

bool overlapMadGraph(EventTree* tree){
	const double Et_cut = 13;
	const double dR_cut = 0.3;
	const double Eta_cut = 3.0;
	bool haveOverlap = false;
	for(int phoInd=0; phoInd<tree->nMC_; ++phoInd){
		if(tree->mcIndex->at(phoInd) > 100) break;
		if(tree->mcPID->at(phoInd)==22 &&
			tree->mcPt->at(phoInd) > Et_cut &&
			fabs(tree->mcEta->at(phoInd)) < Eta_cut &&
			dR(tree->mcEta->at(phoInd),tree->mcPhi->at(phoInd),tree->mcMomEta->at(phoInd),tree->mcMomPhi->at(phoInd)) > dR_cut &&
			(tree->mcParentage->at(phoInd)==2 || tree->mcParentage->at(phoInd)==10 || tree->mcParentage->at(phoInd)==26)) haveOverlap = true;
	}
	return haveOverlap;
}

