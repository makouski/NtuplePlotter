#ifndef PUREWEIGHT_H
#define PUREWEIGHT_H

#include<TH1F.h>
#include<TFile.h>
#include<iostream>

#define PUFileName_Def "Pileup_Observed_69300.root"

class PUReweight{
public:
	PUReweight(int nFiles, char** fileNames);
	~PUReweight();
	double getWeight(int nPUInfo, int puBX[], int nPU[]);
	double getAvgWeight();
	
private:
	double PUweightSum;
	long int events;
	TH1D* PUweightHist;
};

#endif
