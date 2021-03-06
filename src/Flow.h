#pragma once
#include "TVector2.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TPad.h"
#include "THStack.h"
#include "TMath.h"
#include "TFile.h"
#include "TLorentzVector.h"

#include "HADES_constants.h"
#include "Centrality.h"
#include "Selector.h"
#include "Qvector3SE.h"
#include "TProfile2D.h"

#define  DATATREE_SHINE
#include "DataTreeEvent.h"

using std::vector;
using std::array;
using std::cout;
using std::endl;

class Flow
{
	protected:
	int fNumberOfSE;
	int fPid;
	int fHarm;
	enum eProfilesQa{
		kCosPhiYcm,
		kSinPhiYcm,
		kNumberOfQaProfiles
	};
	enum e1dQa{
		kCentrality = 0,
		kCentralityEstimator,
		kNumberOf1dQa
	};
	enum e2dQa{
		kM2VsP,
		kPtVsYforward,
		kPtVsYbackward,
		kFwAdcVsEstimator,
		kFwZVsEstimator,
		kFwAdcVsModuleId,
		kFwZVsModuleId,
		kFwTimeVsModuleId,
		kFwBetaVsModuleId,
		kNumberOf2dQa
	};
	vector<TVector2> fUvector;
	vector<TVector2> fFlow;
	DataTreeEvent* fEvent;
	Centrality* fCentrality;
	Selector* fSelector;
	
	array<TH1F*, kNumberOf1dQa> h1dQa;
	array<TProfile*, kNumberOfQaProfiles> hProfileQa;
	array<TH2F*, kNumberOf2dQa> h2dQa;
	array<TProfile2D*, 2> hMeanUn;
	
	
	virtual void InitializeQvectorCorrelations() {};
			void InitializeQaHistograms();
			void InitializeMeanUn();
	virtual void InitializeObservableFlow() {};
	virtual void FillPtDependence(int trackIdx) {};
			void FillQaHistograms(bool channelSelection=0,  TString signal="adc", float minSignal=0, float maxSignal=9999);
	virtual void FillYDependence(int trackIdx) {};
	virtual void SavePtDependence(TString sFileName) {};
	virtual void SaveYDependence(TString sFileName) {};
	Flow() {};
	public:
	~Flow() {};
	virtual void 		Estimate() {};
	virtual void		InitializeHistograms() {};
			void		FillMeanUn();
	virtual void		SavePictures(TString sFileName) {};
	virtual void		SaveHistogramsToRootFile(TString sFileName) {};
};