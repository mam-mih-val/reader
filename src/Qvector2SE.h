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
#include "Qvector.h"

#define  DATATREE_SHINE
#include "DataTreeEvent.h"

using std::vector;
using std::cout;
using std::endl;

class Qvector2SE : public Qvector
{
    private:
	float			fResolution;
	vector<TGraph*> hResolutionX;
	vector<TGraph*> hResolutionY;
	TProfile*		hMeanCosine;
	TGraph*			hResolutionEP;
    Qvector() {};
	void		EstimateResolution();
	public:
    Qvector2SE(DataTreeEvent* _fEvent, Centrality* _centrality, unsigned int NumSE=2);
    ~Qvector2SE();
	void		Compute(); // One computes Qvector without filling any histograms
    void        ComputeCorrections();
	void		ComputeCorrelations();
	void		ComputeResolution(); 
	void		Estimate();
	void		EsimateResolution();
	void		InitializeHistograms();
	float		Resolution() { return fResolution; }
	void		SaveHistograms(TString sFileName);
	void		SaveHistogramsToROOTFile(TString sFileName);
};