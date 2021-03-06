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
#include "Selector.h"

#define  DATATREE_SHINE
#include "DataTreeEvent.h"

using std::vector;
using std::cout;
using std::endl;

class Qvector
{
    protected:
    unsigned int iNumberOfSE;
	bool fChannelSelection;
	TString fSignal;
	float fMinSignal;
	float fMaxSignal;
	int	fHarm;
	Selector* fSelector;
	Centrality* fCentrality;
	DataTreeEvent* fEvent;
    vector<TVector2> fQvector;
	vector<TProfile*> hMeanQx;
	vector<TProfile*> hMeanQy;
	vector<TProfile*> hCorrelation;
	vector<TH1F*> hQx;
	vector<TH1F*> hQy;
	vector<TH1F*> hPsiEP;

    Qvector() {};
	void	Recenter();

	public:
    Qvector(DataTreeEvent* _fEvent, Centrality* _centrality, unsigned int NumSE=2);
    ~Qvector() {};
			TVector2	At(int se=0){ return fQvector.at(se); }
    virtual void        Estimate() {};
	virtual void		Compute() {}; // One computes Qvector without filling any histograms
	virtual void		ComputeResolution() {}; 
	virtual void		ComputeCorrelations() {};
    virtual void        ComputeCorrections() {};
	virtual void		EstimateResolution() {};
	vector<TProfile*>	GetCorrelationsHistograms() { return hCorrelation; }
			float		X(int se=0) { return fQvector.at(se).X(); }
			float		Y(int se=0) { return fQvector.at(se).Y(); }
			float		Psi(int se=0) { return fQvector.at(se).Phi(); }
			bool		ChannelSelection() { return fChannelSelection; }
			TString		Signal() { return fSignal; }
			float		MinSignal() { return fMinSignal; }
			float		MaxSignal() { return fMaxSignal; }
	virtual void		InitializeHistograms() {};
	virtual void		SavePictures(TString sFileName) {};
	virtual void		SaveHistogramsToROOTFile(TString sFileName);
};
