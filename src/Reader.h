#pragma once
#include <iostream>

#include "TROOT.h"
#include "TH2F.h"
#include "TProfile.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TChain.h"
#include "TMath.h"
#include "TVector3.h"

#include "Selector.h"
#include "Qvector.h"
#include "EventQA.h"
#include "TrackQA.h"

#include "HADES_constants.h"

class Reader
{
    private:
    TChain* fChain;
    TBranch* DTEvent;
    DataTreeEvent* fEvent;
    //Selector selector;
    
    enum eFlowProfile{
        meanQx = 0,
        meanQy,
        resolution,
        yMostCentral,
        yMidCentral,
        yPeripherial,
        ptMostCentral,
        ptMidCentral,
        ptPeripherial,
        NumOfFLowProfiles
    };
    enum eQvectorDistribution{
        QxNotRecentred=0,
        QyNotRecentred,
        QxRecentred,
        QyRecentred,
        NumOfQvectorHistos
    };
	enum eParticleMap{
		all = 0,
		electron,
		positron,
		pi_minus,
		pi_plus,
		proton,
		deuteron,
		helium3,
		helium4,
		NumOfParticles
	};
    TProfile* pFlowProfiles[NumOfFLowProfiles];
    TH1F*   vQvectorDistribution[4];
    public:
    Reader() {};
    Reader(char* cFileName);
    ~Reader();
    void            AddFile(char* cFileName);
    void            FillCorrectionHistos();
    DataTreeEvent*  GetEvent(int idx=0);
    void            GetFlow(int iNumHarm=1);
    void            BuildQAHistograms(TString sPicName);
    void            InitFlowHistos();
	void			BuildQvectorHistograms(TString sPicName);
    void            SaveFlowStatistics();
};