#pragma once
#include "TVector2.h"
#include "TProfile3D.h"
#include "TCanvas.h"
#include "TPad.h"
#include "THStack.h"
#include "TStyle.h"
#include "TMath.h"
#include "TFile.h"
#include "TLorentzVector.h"

#include "HADES_constants.h"
#include "Centrality.h"
#include "Selector.h"
#include "Qvector.h"

#define  DATATREE_SHINE
#include "DataTreeEvent.h"
#include "Flow.h"

using std::vector;
using std::cout;
using std::endl;

class Flow3SE : public Flow
{
	protected:
	Qvector3SE* fQvector;
	array< vector<TProfile*>, 3 > xRapidity;
	array< vector<TProfile*>, 3 > yRapidity;
	array< vector<TProfile*>, 3 > xPt;
	array< vector<TProfile*>, 3 > yPt;
	Flow3SE() {};
	public:
	Flow3SE(DataTreeEvent* _event, Centrality* _centrality, Qvector3SE* _Qvector, Selector* _selector, int _pid=14);
	~Flow3SE();
	void 		Estimate();
	void		InitializeHistograms();
	void		SavePictures(TString sFileName);
	void		SaveHistogramsToRootFile(TString sFileName) {};
};