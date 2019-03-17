#pragma once
#include "TVector2.h"
#include "TH1F.h"
#include "TProfile.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TMath.h"

#include "HADES_constants.h"

#define  DATATREE_SHINE
#include "DataTreeEvent.h"

class Qvector
{
    private:
    Int_t iNumberOfSE;
	float fQ[2][2];
    std::vector<TVector2> fQvector(iNumberOfSE);
	enum eFlowProfile{
        meanQx1 = 0,
        meanQy1, // 1
		meanQx2, // 2
        meanQy2, // 3
        resolution, // 4
        NumOfFLowProfiles // 5
    };
     enum eQvectorDistribution{
        QxNotRecentred1=0,
		QxNotRecentred2, // 1
        QyNotRecentred1, // 2
		QyNotRecentred2, // 3
        QxRecentred1, // 4
        QxRecentred2, // 5
		QyRecentred1, // 6
		QyRecentred2,  // 7
		PsiEPNotRecentred1, // 8
		PsiEPNotRecentred2, // 9
		PsiEPRecentred1, // 10
		PsiEPRecentred2, // 11
        NumOfQvectorHistos // 12
    };
	enum eCanvasMap{
		MeanQvectors = 0,
		QvectorsDistribution,
		NumOfCanvases
	};
	TProfile* vProfile[NumOfFLowProfiles];
    TH1F*   vHisto1D[NumOfQvectorHistos];
	TCanvas* vCanvas[NumOfCanvases];
	
	vector<TProfile*> hMeanQx(iNumberOfSE);
	vector<TProfile*> hMeanQy(iNumberOfSE);
	vector<TH1F*> hQx(2*iNumberOfSE);
	vector<TH1F*> hQy(2*iNumberOfSE);
	vector<TH1F*> hPsiEP(2*iNumberOfSE); 

	void Estimate2SE(DataTreeEvent* fEvent);
	void Estimate3SE(DataTreeEvent* fEvent);
	public:
    Qvector();
    ~Qvector() {};
    void        FillCorrections(DataTreeEvent* fEvent);
	void		SetNumberOfSE(int iNum) { iNumberOfSE = iNum; }
    void        Estimate(DataTreeEvent* fEvent);
	void		InitHistograms();
	void		SaveHistograms(TString sPicName);
    Float_t     GetComponent(int i, int j=0);
    Float_t     GetPsiEP(int j=0);
};
