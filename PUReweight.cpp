#include"PUReweight.h"

PUReweight::PUReweight(int nFiles, char** fileNames){
	PUweightSum = 0.0;
	events = 0;
	TFile* pileupFile = new TFile(PUFileName_Def,"READ");
	PUweightHist = (TH1D*)pileupFile->Get("pileup");
	PUweightHist->SetDirectory(0);
	pileupFile->Close();
	double PUweightInt = PUweightHist->Integral();
	TH1F* mcPU=NULL;
	for(int nmcfile = 0; nmcfile<nFiles; nmcfile++){
		std::cout << "reading file " << std::string(fileNames[nmcfile]) << std::endl;
		TFile* mcFile = new TFile(fileNames[nmcfile],"READ");
		if( mcPU==NULL) mcPU = (TH1F*)mcFile->Get("ggNtuplizer/hPU");
		else mcPU->Add((TH1F*)mcFile->Get("ggNtuplizer/hPU"));
		mcPU->SetDirectory(0);
		mcFile->Close();
	}
	mcPU->Scale(1.0/mcPU->Integral());
	PUweightHist->Divide(mcPU);
	PUweightHist->Scale(1.0/PUweightInt);
	delete mcPU;
}

PUReweight::~PUReweight(){
	delete PUweightHist;
}

double PUReweight::getWeight(int nPUInfo, std::vector<int> *puBX, std::vector<int> *nPU){
	double PUweight=0.0;
	for(int puInd=0; puInd<nPUInfo; ++puInd){
		if( puBX->at(puInd) == 0 ){
			PUweight = PUweightHist->GetBinContent(PUweightHist->GetXaxis()->FindBin(nPU->at(puInd)));
			break;
		}
	}
	events++;
	PUweightSum+=PUweight;
	return PUweight;
}

double PUReweight::getAvgWeight(){
	if(events!=0) return PUweightSum/events;
	else return -1.0;
}