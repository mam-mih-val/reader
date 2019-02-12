#pragma once
#include <iostream>

#include "TROOT.h"
#include "TH2F.h"
#include "TProfile.h"
#include "TCanvas.h"
#include "TChain.h"
#include "TMath.h"

#include "Selector.h"
#include "Qvector.h"

#include "HADES_constants.h"

class Reader
{
    private:
    TChain* fChain;
    TBranch* DTEvent;
    DataTreeEvent* fEvent;
    Selector selector;
    Qvector fQ;
    enum eQAHisto1DMap{
        tracksMDC = 0,
        tracksMDC_selected,
        hitsTOF,
        hitsTOF_selected,
        chargeFW,
        chargeFW_selected,
        vertexZ,
        vertexZ_selected,
        ptMDC,
        ptMDC_selected,
        massTOF,
        massTOF_selected,
        rapidityMDC,
        rapidityMDC_selected,
        rapidityMDC_recentred,
        phiMDC,
        phiMDC_selected,
        betaTOF,
        betaTOF_selected,
        Num1DHistos
    };
    enum eQAHisto2DMap{
        tracks_hits = 0,
        tracks_hits_selected,
        tracks_charge,
        tracks_charge_selected,
        hits_charge,
        hits_charge_selected,
        vertexX_vertexY,
        vertexX_vertexY_selected,
        hitsFW_X_Y,
        hitsFW_X_Y_selected,
        Num2DHistos
    };
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
    
    TH1F* vHisto1D[Num1DHistos];
    TH2F* vHisto2D[Num2DHistos];
    TProfile* pFlowProfiles[NumOfFLowProfiles];
    TH1F*   vQvectorDistribution[4];
    public:
    Reader() {};
    Reader(char* cFileName);
    ~Reader();
    void            AddFile(char* cFileName);
    void            DrawQA1DHistos(char* cPictureName);
    void            DrawQA2DHistos(char* cPictureName);
    void            FillCorrectionHistos();
    DataTreeEvent*  GetEvent(int idx=0);
    void            GetEvent();
    void            GetFlow(int iNumHarm=1);
    void            GetQualityAssurance(Int_t iPT=HADES_constants::kPT2);
    void            InitQAHistos();
    void            InitFlowHistos();
    void            SaveFlowStatistics();
    void            SaveQAStatistics();
};
