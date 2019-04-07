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
	vector<TVector2> fResolution;
	vector<TProfile*> hMeanQx;
	vector<TProfile*> hMeanQy;
	vector<TProfile*> hCorrelation;
	vector<TH1F*> hQx;
	vector<TH1F*> hQy;
	vector<TH1F*> hPsiEP;
	vector<TH1F*> hDeltaPsiEP;
	vector<TGraph*> hResolutionX;
	vector<TGraph*> hResolutionY;
	TProfile*	hMeanCosine;
	TGraph* 	hResolutionEP;

    Qvector() {};
	void	ComputeResolution3SE();
	void 	Estimate2SE();
	void 	Estimate3SE();
	void	FillCorrelations2SE();
	void	FillCorrelations3SE();
	void	FillScalarProduct3SE();
	void	FillMeanCosine2SE();
	void	FillResolutionProfile();
	void	Init2SEHistograms();
	void	Init3SEHistograms();
	void	EstimateResolution2SE();
	void	EstimateResolution3SE();
	void	Recenter();
	
	public:
    Qvector(DataTreeEvent* _fEvent, Centrality* _centrality, unsigned int NumSE=2);
    ~Qvector();
	TVector2	At(int se=0){ return fQvector.at(se); }
    void        Estimate();
	void		Compute(); // One computes Qvector without filling any histograms
	void		ComputeResolution(); 
	void		EstimateResolution();
    void        FillCorrections();
    Float_t     GetComponent(int i, int j=0);
	float		X(int se=0) { return fQvector.at(se).X(); }
	float		Y(int se=0) { return fQvector.at(se).Y(); }
	float		Psi(int se=0) { return fQvector.at(se).Phi(); }
    Float_t     GetPsiEP(unsigned int iSubEv=0) { return fQvector.at(iSubEv).Phi(); }
	TVector2	GetQvector(unsigned int iSubEv=0) { return fQvector.at(iSubEv); }
	TVector2	Resolution(int se=0) { return fResolution.at(se); }
	void		InitHistograms();
	void		LoadCentrality(Centrality* _centrality) { fCentrality = _centrality; }
	void		SaveHistograms(TString sPicName);
	void		SaveHistogramsToROOTFile(TString sFileName);
	void		SetNumberOfSE(unsigned int iNum) { iNumberOfSE = iNum; }
};
