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
#include "Qvector3SE.h"
#include "EventQA.h"
#include "TrackQA.h"
#include "Centrality.h"
#include "Flow.h"
#include "Flow3SE.h"

#include "HADES_constants.h"

using std::vector;
using std::cout;
using std::endl;

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
    Reader(TString cFileName);
    ~Reader();
    void            AddFile(char* cFileName);
    void            FillCorrectionHistos();
    DataTreeEvent*  GetEvent(int idx=0);
    void            GetFlow(int iNumHarm=1);
    void            BuildEventQaHistograms(TString sPicName);
    void            BuildTrackQaHistograms(TString sPicName, int pid=14);
    void            InitFlowHistos();
	void			BuildQvector3SeHistograms(TString sPicName, bool channelSelection=0, TString signal="adc", float minSignal=0, float maxSignal=9999, int harm=1);
	void			BuildFlow3SeHistograms(TString sPicName, bool channelSelection=0, TString signal="adc", float minSignal=0, float maxSignal=9999, int pid=14, int harm=1);
    void            SaveFlowStatistics();
};