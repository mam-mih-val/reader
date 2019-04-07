#pragma once
#include "TVector2.h"
#include "TProfile3D.h"
#include "TCanvas.h"
#include "TPad.h"
#include "THStack.h"
#include "TMath.h"
#include "TFile.h"
#include "TLorentzVector.h"

#include "HADES_constants.h"
#include "Centrality.h"
#include "Selector.h"
#include "Qvector.h"

#define  DATATREE_SHINE
#include "DataTreeEvent.h"

using std::vector;
using std::cout;
using std::endl;

class Flow
{
	private:
	int fNumberOfSE;
	int fPid;
	vector<TVector2> fFlow;
	Centrality* fCentrality;
	DataTreeEvent* fEvent;
	TProfile* hRapidity;
    Qvector* fQvector;
	Selector* fSelector;
	vector<TProfile3D*> hFlowX;
	vector<TProfile3D*> hFlowY;
	
	Flow();
	void _Estimae3SE(int trackIdx);
	
	public:
	Flow(DataTreeEvent* _event, Centrality* _centrality, Qvector* _Qvector, Selector* _selector, int _pid, int numberOfSE=3);
	~Flow();
	void 		Estimate();
	void		InitializeHistograms();
	void		Estimae();
	void		SavePictures(TString sFileName);
	void		SaveHistogramsToRootFile(TString sFileName);
};