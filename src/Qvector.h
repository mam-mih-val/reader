#pragma once
#include "TVector2.h"
#include "TH1F.h"
#include "TProfile.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TPad.h"
#include "THStack.h"
#include "TMultiGraph.h"
#include "TLegend.h"
#include "TMath.h"
#include "TFile.h"

#include "HADES_constants.h"
#include "Centrality.h"

#define  DATATREE_SHINE
#include "DataTreeEvent.h"

using std::vector;
using std::cout;
using std::endl;

class Qvector
{
    private:
    unsigned int iNumberOfSE;
	Centrality* fCentrality;
	DataTreeEvent* fEvent;
    vector<TVector2> fQvector;
	vector<TProfile*> hMeanQx;
	vector<TProfile*> hMeanQy;
	vector<TProfile*> hCorrelation;
	vector<TProfile*> hMeanCosine;
	vector<TProfile*> hScalarProduct;
	vector<TH1F*> hQx;
	vector<TH1F*> hQy;
	vector<TH1F*> hPsiEP;
	vector<TH1F*> hDeltaPsiEP; 
	vector<TGraph*> hResolution;
	vector<TGraph*> hResolutionX;
	vector<TGraph*> hResolutionY;

    Qvector() {};
	void 	Estimate2SE();
	void 	Estimate3SE();
	void	FillCorrelations2SE();
	void	FillCorrelations3SE();
	void	FillMeanCosine2SE();
	void	FillScalarProduct3SE();
	void	FillResolutionProfile();
	void	Init2SEHistograms();
	void	Init3SEHistograms();
	void	EstimateResolution3SE();
	void	Recenter();
	
	public:
    Qvector(DataTreeEvent* _fEvent, Centrality* _centrality, unsigned int NumSE=2);
    ~Qvector();
    void        Estimate();
	void		EstimateResolution();
    void        FillCorrections();
    Float_t     GetComponent(int i, int j=0);
    Float_t     GetPsiEP(unsigned int iSubEv=0) { return fQvector.at(iSubEv).Phi(); }
	TVector2	GetQvector(unsigned int iSubEv=0) { return fQvector.at(iSubEv); }
	float		GetResolution();
	void		InitHistograms();
	void		LoadCentrality(Centrality* _centrality) { fCentrality = _centrality; }
	void		SaveHistograms(TString sPicName);
	void		SaveHistogramsToROOTFile(TString sFileName);
	void		SetNumberOfSE(unsigned int iNum) { iNumberOfSE = iNum; }
};
