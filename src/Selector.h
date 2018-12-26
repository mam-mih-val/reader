#include "HADES_constants.h"

class Selector
{
    private:
    DataTreeEvent*  fEvent;
    TH1F*           hRejectedEvents;
    TH1F*           hRejectedTracks;
    enum event_cuts{
        cVeretexPositionZ = 0,
        cVeretexPositionXY,
        cVertexQuality,
        cTriggerVertexClust,
        cTriggerVertexCand,
        cTriggerGoodStart,
        cTriggerNoPileUp,
        cTriggerGoodStartVeto,
        cTriggerGoodStartMeta,
        cTriggerNoVeto,
        cNumOfEventCuts
    };
    enum track_cuts{
        cDCA = 0,
        cTrackHitMatchX,
        cTrackHitMatchY,
        cChi2,
        cBeta,
        cNumOfTrackCuts
    };
    public:
    Selector();
    ~Selector();
    Bool_t IsCorrectEvent(DataTreeEvent* _fEvent);
    Bool_t IsCorrectTrack(Int_t idx);
    void    SaveStatistics();
};