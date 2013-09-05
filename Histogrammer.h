#ifndef HISTOGRAMMER_H
#define HISTOGRAMMER_H

#include<iostream>
#include<vector>
#include<map>
#include<string>
#include"EventTree.h"
#include"Selector.h"
#include<TFile.h>
#include<TH1F.h>
#include<TH2F.h>
#include<TMath.h>
#include<TLorentzVector.h>

class Histogrammer{
public:
	Histogrammer(const char* titleIn);
	~Histogrammer(void);
	void fill(Selector* selector_, EventTree* tree, double weight);
	void write_histograms(const char* folder, std::vector<TH1F*> histVector);
	
private:
	std::string title;
	std::map< std::string, TH1F* > hists;
	std::map< std::string, TH2F* > hists2d;

	void make_hist(const char* hname, const char* htitle, int nbins, double xlow, double xhigh, const char* xlabel, const char* ylabel);
	int minDrIndex(double myEta, double myPhi, std::vector<int> Inds, Float_t* etas, Float_t* phis);
	double minDr(double myEta, double myPhi, std::vector<int> Inds, Float_t* etas, Float_t* phis);
	double matchedJetBtag(double phoEta, double phoPhi, std::vector<int> Inds, Float_t* jetEtas, Float_t* jetPhis, Float_t* jetCSV);
	double phoJetmass(double phoEt, double phoEta, double phoPhi, std::vector<int> Inds, Float_t* jetPts, Float_t* jetEtas, Float_t* jetPhis);
	double minDrPhoB(int PhoInd, EventTree* tree);
	double calc_ht(Selector* selector_, EventTree* tree);
};
#endif
