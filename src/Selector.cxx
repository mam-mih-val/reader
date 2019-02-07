#include "Selector.h"


Selector::Selector() 
{
    fEvent = new DataTreeEvent;
    hRejectedEvents = new TH1F("Amount of rejected events","",cNumOfEventCuts,0,8);
    hRejectedTracks = new TH1F("Amount of rejected tracks","",cNumOfTrackCuts,0,8);
}

Selector::~Selector()
{
    delete fEvent;
    delete hRejectedEvents;
    delete hRejectedTracks;
}

Bool_t Selector::IsCorrectEvent(DataTreeEvent* _fEvent, int iPT)
{
    fEvent = _fEvent;
    if( !fEvent->GetTrigger(iPT)->GetIsFired() )
    {
        return 0;
    }
    
    if (  fEvent->GetVertexPositionComponent(2) > 0 || fEvent->GetVertexPositionComponent(2) < -60 )
    {
        hRejectedEvents->Fill(cVeretexPositionZ);
        return 0;
    }
    Float_t Rx = fEvent->GetVertexPositionComponent(0), Ry = fEvent->GetVertexPositionComponent(1);
    if ( sqrt(Rx*Rx+Ry*Ry) > 3 )
    {
        hRejectedEvents->Fill(cVeretexPositionXY);
        return 0;
    }
    if ( fEvent->GetVertexQuality() < 0.5 || fEvent->GetVertexQuality() > 40 )
    {
        hRejectedEvents->Fill(cVertexQuality);
        return 0;
    }
    if ( !fEvent->GetTrigger(HADES_constants::kGoodVertexClust)->GetIsFired() ) 
    {
        hRejectedEvents->Fill(cTriggerVertexClust);
        return 0;
    }
    if ( ! fEvent->GetTrigger(HADES_constants::kGoodVertexCand)->GetIsFired() )
    {
        hRejectedEvents->Fill(cTriggerVertexCand);
        return 0;
    }
    if( !fEvent->GetTrigger(HADES_constants::kGoodSTART)->GetIsFired() )
    {
        hRejectedEvents->Fill(cTriggerGoodStart);
        return 0;
    }
    if( !fEvent->GetTrigger(HADES_constants::kNoPileUpSTART)->GetIsFired() )
    {
        hRejectedEvents->Fill(cTriggerNoPileUp);
        return 0;
       
    }
    if( !fEvent->GetTrigger(HADES_constants::kGoodSTARTVETO)->GetIsFired() )
    {
        hRejectedEvents->Fill(cTriggerGoodStartVeto);
        return 0;
    }
    if( !fEvent->GetTrigger(HADES_constants::kGoodSTARTMETA)->GetIsFired() )
    {
        hRejectedEvents->Fill(cTriggerGoodStartMeta);
        return 0;
    }
    if( !fEvent->GetTrigger(HADES_constants::kNoVETO)->GetIsFired() )
    {
        hRejectedEvents->Fill(cTriggerNoVeto);
        return 0;
    }
    return 1;
}

Bool_t Selector::IsCorrectTrack(Int_t idx)
{
    DataTreeTrack* fTrack = fEvent->GetVertexTrack(idx);
    DataTreeTOFHit* fHit = fEvent->GetTOFHit(idx);
    Float_t fTof = fHit->GetTime();
    Float_t fLen = fHit->GetPathLength();
    //cout << len/tof/299.792458 << endl;
    if ( fTrack->GetDCAComponent(1) > 15 )
    {
        hRejectedTracks->Fill(cDCA);
        return 0;
    }
    if ( fLen/fTof/299.792458 > 1 ) 
    {
        hRejectedTracks->Fill(cBeta);
        return 0;
    }
    if ( fHit->GetPositionComponent(0) < -5 || fHit->GetPositionComponent(0) > 5 )
    {
        hRejectedTracks->Fill(cTrackHitMatchX);
        return 0;
    }
    if ( fHit->GetPositionComponent(1) < -5 || fHit->GetPositionComponent(1) > 5 )
    {
        hRejectedTracks->Fill(cTrackHitMatchY);
        return 0;
    }    
    if ( fTrack->GetChi2() > 100 )
    {
        hRejectedTracks->Fill(cChi2);
        return 0;
    }    
    return 1;
}

void Selector::SaveStatistics() 
{  
    TFile* fFile = new TFile("../histograms/SelectorStatisticsHistos.root","recreate");
    hRejectedEvents->Write();
    hRejectedTracks->Write();
    fFile->Close();
}