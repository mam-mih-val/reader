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
    public:
    TH1F* vHisto1D[Num1DHistos];
    TH2F* vHisto2D[Num2DHistos];
    TProfile* pFlowProfiles[NumOfFLowProfiles];
    TH1F*   vQvectorDistribution[4];
    Reader() {};
    ~Reader();
    Reader(char* cFileName);
    void            Init();
    DataTreeEvent*  GetEvent(int idx);
    void            loop();
    void            GetEvent();
    void            SaveQAStatistics();
    void            GetQualityAccurance();
    void            InitQAHistos();
    void            InitFlowHistos();
    void            GetFlow(int iNumHarm=1);
    void            SaveFlowStatistics();
    void            GetQ();
};