#ifndef PUREWEIGHT_H
#define PUREWEIGHT_H

#include<TH1F.h>
#include<TFile.h>
#include<iostream>
#include<vector>


class PUReweight{
public:
	PUReweight(TString fdata,TString fileNames);
	~PUReweight();
	double getWeight(int nPUInfo, std::vector<int> *puBX, std::vector<int> *nPU);
	double getAvgWeight();
	
private:
	double PUweightSum;
	long int events;
	TH1D* PUweightHist;
};

PUReweight::PUReweight(TString fdata,TString fileNames){
	PUweightSum = 0.0;
	events = 0;
        TFile* pileupFile = new TFile(fdata,"READ");
	PUweightHist = (TH1D*)pileupFile->Get("pileup");
	//	TCanvas* c1;
	//	c1->cd();
	//	PUweightHist->Draw();
	PUweightHist->SetDirectory(0);
	pileupFile->Close();
        TH1D* PUbackup;
if(PUweightHist->GetNbinsX() != 100){
std::cout << "Wrong number of bins in the pileup histogram" << std::endl;
PUbackup = new TH1D("pileup_new","pileup_new",100,0,100);
for(int ibin=1; ibin <= PUweightHist->GetNbinsX(); ibin++){
PUbackup->SetBinContent(ibin, PUweightHist->GetBinContent(ibin));
// assuming the same scale
}
PUweightHist = PUbackup;
}
	double PUweightInt = PUweightHist->Integral();
	TH1F* mcPU=NULL;
	//	std::cout << "reading file " << PUFileName_Def<<'\t' <<fileNames << std::endl;
	TFile* mcFile = new TFile(fileNames,"READ");
//	if( mcPU==NULL) mcPU = (TH1F*)mcFile->Get("ggNtuplizer/hPU");
//	else mcPU->Add((TH1F*)mcFile->Get("ggNtuplizer/hPU"));
	if( mcPU==NULL) mcPU = (TH1F*)mcFile->Get("hPU");
	else mcPU->Add((TH1F*)mcFile->Get("hPU"));
	mcPU->SetDirectory(0);
	mcFile->Close();
	
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

#endif
