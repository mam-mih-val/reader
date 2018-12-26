#include "TROOT.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "HADES_constants.h"
#include "DataTreeEvent.h"

class Reader
{
    private:
    TChain* fChain;
    TBranch* DTEvent;
    DataTreeEvent* fEvent;
    enum eHisto1DMap{
        tracksMDC = 0,
        tracksMDC_selected,
        hitsTOF,
        hitsTOF_selected,
        chargeFW,
        chargeFW_selected,
        Num1DHistos
    };
    enum eHisto2DMap{
        tracks_hits = 0,
        tracks_hits_selected,
        tracks_charge,
        tracks_charge_selected,
        hits_charge,
        hits_charge_selected,
        Num2DHistos
    };
    public:
    TH1F* vHisto1D[Num1DHistos];
    TH2F* vHisto2D[Num2DHistos];
    Reader() {};
    ~Reader();
    Reader(char* cFileName);
    void            Init();
    DataTreeEvent*  GetEvent(int idx);
    void            loop();
    void            GetEvent();
    void            SaveStatistics();
    void            GetQualityAccurance();
    void            InitHistos();
};