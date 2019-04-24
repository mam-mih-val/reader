#pragma once
#include "TVector2.h"
#include "TH1F.h"
#include "TProfile.h"
#include "TCanvas.h"
#include "TGraphErrors.h"
#include "TPad.h"
#include "THStack.h"
#include "TStyle.h"
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
using std::array;
using std::cout;
using std::endl;

class Qvector3SE : public Qvector
{
    private:
	vector<TVector2> fResolution;
	vector<TGraphErrors*> hResolutionX;
	vector<TGraphErrors*> hResolutionY;
	vector<TProfile*> hCorrMult;
	void		EstimateResolution();

	public:
    Qvector3SE(DataTreeEvent* _fEvent, Centrality* _centrality);
    ~Qvector3SE();
	void		Compute(); // One computes Qvector without filling any histograms
    void        ComputeCorrections();
	void		ComputeCorrelations();
	void		ComputeResolution(); 
	void		Estimate();
	void		InitializeHistograms();
	TVector2	Resolution(int se=0) { return fResolution.at(se); }
	void		SavePictures(TString sFileName);
	void		SaveHistogramsToROOTFile(TString sFileName);
};
