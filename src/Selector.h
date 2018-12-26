#include "TROOT.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "HADES_constants.h"
#include "DataTreeEvent.h"

class Selector
{
    private:
    DataTreeEvent*  fEvent;
    TH1F*           hRejected;
    enum cuts{
        cVeretexPositionZ = 0,
        cVeretexPositionXY,
        cVertexQuality,
        cTriggerVertexClust,
        cTriggerVertexCand,
        cTriggerGoodStart,
        cTriggerNoPileUp,
        cTriggerGoodStartVeto,
        cTriggerGoodStartMeta,
        cTriggerNoVeto
    };
    public:
    Selector();
    ~Selector();
    Bool_t IsCorrect(DataTreeEvent* _fEvent);
    void    SaveStatistics();
};