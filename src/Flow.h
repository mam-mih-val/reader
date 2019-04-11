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
#include "Qvector3SE.h"

#define  DATATREE_SHINE
#include "DataTreeEvent.h"

using std::vector;
using std::cout;
using std::endl;

class Flow
{
	protected:
	int fNumberOfSE;
	int fPid;
	vector<TVector2> fUvector;
	DataTreeEvent* fEvent;
	Centrality* fCentrality;
	Selector* fSelector;

	Flow() {};
	public:
	~Flow() {};
	virtual void 		Estimate() {};
	virtual void		InitializeHistograms() {};
	virtual void		SavePictures(TString sFileName) {};
	virtual void		SaveHistogramsToRootFile(TString sFileName) {};
};