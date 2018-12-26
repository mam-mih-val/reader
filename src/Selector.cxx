#include "Selector.h"


Selector::Selector() 
{
    fEvent = new DataTreeEvent;
    hRejected = new TH1F("Amount of rejected events","",8,0,8);
}

Selector::~Selector()
{
    delete fEvent;
    delete hRejected;
}

Bool_t Selector::IsCorrect(DataTreeEvent* _fEvent)
{
    fEvent = _fEvent;
    if (  fEvent->GetVertexPositionComponent(2) > 0 || fEvent->GetVertexPositionComponent(2) < -60 )
    {
        hRejected->Fill(cVeretexPositionZ);
        return 0;
    }
    Float_t Rx = fEvent->GetVertexPositionComponent(0), Ry = fEvent->GetVertexPositionComponent(1);
    if ( sqrt(Rx*Rx+Ry*Ry) > 3 )
    {
        hRejected->Fill(cVeretexPositionXY);
        return 0;
    }
    if ( fEvent->GetVertexQuality() < 0.5 || fEvent->GetVertexQuality() > 40 )
    {
        hRejected->Fill(cVertexQuality);
        return 0;
    }
    if ( !fEvent->GetTrigger(HADES_constants::kGoodVertexClust)->GetIsFired() ) 
    {
        hRejected->Fill(cTriggerVertexClust);
        return 0;
    }
    if ( ! fEvent->GetTrigger(HADES_constants::kGoodVertexCand)->GetIsFired() )
    {
        hRejected->Fill(cTriggerVertexCand);
        return 0;
    }
    if( !fEvent->GetTrigger(HADES_constants::kGoodSTART)->GetIsFired() )
    {
        hRejected->Fill(cTriggerGoodStart);
        return 0;
    }
    if( !fEvent->GetTrigger(HADES_constants::kNoPileUpSTART)->GetIsFired() )
    {
        hRejected->Fill(cTriggerNoPileUp);
        return 0;
       
    }
    if( !fEvent->GetTrigger(HADES_constants::kGoodSTARTVETO)->GetIsFired() )
    {
        hRejected->Fill(cTriggerGoodStartVeto);
        return 0;
    }
    if( !fEvent->GetTrigger(HADES_constants::kGoodSTARTMETA)->GetIsFired() )
    {
        hRejected->Fill(cTriggerGoodStartMeta);
        return 0;
    }
    if( !fEvent->GetTrigger(HADES_constants::kNoVETO)->GetIsFired() )
    {
        hRejected->Fill(cTriggerNoVeto);
        return 0;
    }
    
    return 1;
}

void Selector::SaveStatistics() 
{  
    TFile* fFile = new TFile("Selector_Statistics.root","recreate");
    hRejected->Write();
    fFile->Close();
}