#include "HADES_constants.h"

class Reader
{
    private:
    TChain* fChain;
    TBranch* DTEvent;
    DataTreeEvent* fEvent;
    Selector selector;
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
    void            FillCorrectionHistos();
    DataTreeEvent*  GetEvent(int idx=0);
    void            GetEvent();
    void            GetFlow(int iNumHarm=1);
    void            GetQualityAccurance();
    void            InitQAHistos();
    void            InitFlowHistos();
    void            SaveFlowStatistics();
    void            SaveQAStatistics();
};