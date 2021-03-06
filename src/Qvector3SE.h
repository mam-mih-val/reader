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
#include "Qvector.h"
#include "Centrality.h"

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
    Qvector3SE(DataTreeEvent* _fEvent, Selector* _selector, Centrality* _centrality,bool _channelSelection=0,  TString _signal="adc", float _minSignal=0, float _maxSignal=9999, int _harm=1);
    ~Qvector3SE();
	void		Compute(); // One computes Qvector without filling any histograms
    void        ComputeCorrections();
	void		ComputeCorrelations();
	void		ComputeResolution(); 
	void		Estimate();
	void		InitializeHistograms();
	vector<TGraphErrors*> GetResolutionXGraph() { return hResolutionX; }
	vector<TGraphErrors*> GetResolutionYGraph() { return hResolutionY; }
	TVector2	Resolution(int se=0) { return fResolution.at(se); }
	void		SavePictures(TString sFileName);
	void		SaveHistogramsToROOTFile(TString sFileName);
};
