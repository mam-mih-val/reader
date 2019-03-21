#pragma once
#include "TVector2.h"
#include "TH1F.h"
#include "TProfile.h"
#include "TCanvas.h"
#include "TPad.h"
#include "THStack.h"
#include "TLegend.h"
#include "TMath.h"
#include "TFile.h"

#include "HADES_constants.h"

#define  DATATREE_SHINE
#include "DataTreeEvent.h"

class Qvector
{
    private:
    unsigned int iNumberOfSE;
	TH1F* hCentralityPercentile;
    vector<TVector2> fQvector;	
	vector<TProfile*> hMeanQx;
	vector<TProfile*> hMeanQy;
	vector<TProfile*> hCorrelation;
	vector<TH1F*> hQx;
	vector<TH1F*> hQy;
	vector<TH1F*> hPsiEP; 

	void 	Estimate2SE(DataTreeEvent* fEvent);
	void 	Estimate3SE(DataTreeEvent* fEvent);
	float	GetCentralityClass(DataTreeEvent* fEvent);
	public:
    Qvector();
    ~Qvector() {};
	void		LoadCentralityPercentile(TString FileName);
    void        FillCorrections(DataTreeEvent* fEvent);
	void		SetNumberOfSE(unsigned int iNum) { iNumberOfSE = iNum; }
    void        Estimate(DataTreeEvent* fEvent);
	void		InitHistograms();
	void		SaveHistograms(TString sPicName);
	void		SaveHistogramsToROOTFile(TString sFileName);
	TVector2	GetQvector(unsigned int iSubEv=0) { return fQvector.at(iSubEv); }
    Float_t     GetComponent(int i, int j=0);
    Float_t     GetPsiEP(unsigned int iSubEv=0) { return fQvector.at(iSubEv).Phi(); }
};
