#include "Selector.h"


Selector::Selector() 
{
    fEvent = new DataTreeEvent;
    hRejectedEvents = new TH1F("Amount of rejected events","",cNumOfEventCuts,0,8);
    hRejectedTracks = new TH1F("Amount of rejected tracks","",cNumOfTrackCuts,0,8);
    hIncorrectEvent = new TH1F("Amount of rejected events on this cut",";code of the cut;counts",cNumOfEventCuts,0,cNumOfEventCuts);
    hIncorrectTracks= new TH1F("Amount of rejected tracks on this cut",";code of the cut;counts",cNumOfTrackCuts,0,cNumOfTrackCuts);
}

Selector::~Selector()
{
    delete fEvent;
    delete hRejectedEvents;
    delete hRejectedTracks;
    delete hIncorrectEvent;
}

Bool_t Selector::IsCorrectEvent(DataTreeEvent* _fEvent, int iPT)
{
    fEvent = _fEvent;
    this->CheckEventCuts(fEvent);
/*    if( !fEvent->GetTrigger(iPT)->GetIsFired() )
    {
        return 0;
    }
*/    
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
    this->CheckTrackCuts(idx);
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

void Selector::CheckEventCuts(DataTreeEvent* _fEvent)
{
    if (  fEvent->GetVertexPositionComponent(2) > 0 || fEvent->GetVertexPositionComponent(2) < -60 )
    {
        hIncorrectEvent->Fill(cVeretexPositionZ);
    }
    Float_t Rx = fEvent->GetVertexPositionComponent(0), Ry = fEvent->GetVertexPositionComponent(1);
    if ( sqrt(Rx*Rx+Ry*Ry) > 3 )
    {
        hIncorrectEvent->Fill(cVeretexPositionXY);
    }
    if ( fEvent->GetVertexQuality() < 0.5 || fEvent->GetVertexQuality() > 40 )
    {
        hIncorrectEvent->Fill(cVertexQuality);
    }
    if ( !fEvent->GetTrigger(HADES_constants::kGoodVertexClust)->GetIsFired() ) 
    {
        hIncorrectEvent->Fill(cTriggerVertexClust);
    }
    if ( ! fEvent->GetTrigger(HADES_constants::kGoodVertexCand)->GetIsFired() )
    {
        hIncorrectEvent->Fill(cTriggerVertexCand);
    }
    if( !fEvent->GetTrigger(HADES_constants::kGoodSTART)->GetIsFired() )
    {
        hIncorrectEvent->Fill(cTriggerGoodStart);
    }
    if( !fEvent->GetTrigger(HADES_constants::kNoPileUpSTART)->GetIsFired() )
    {
        hIncorrectEvent->Fill(cTriggerNoPileUp);
    }
    if( !fEvent->GetTrigger(HADES_constants::kGoodSTARTVETO)->GetIsFired() )
    {
        hIncorrectEvent->Fill(cTriggerGoodStartVeto);
    }
    if( !fEvent->GetTrigger(HADES_constants::kGoodSTARTMETA)->GetIsFired() )
    {
        hIncorrectEvent->Fill(cTriggerGoodStartMeta);
    }
    if( !fEvent->GetTrigger(HADES_constants::kNoVETO)->GetIsFired() )
    {
        hIncorrectEvent->Fill(cTriggerNoVeto);
    }
    if( fEvent->GetNTOFHits()<5 )
    {
        hIncorrectEvent->Fill(cPT2);
    }
    if( !fEvent->GetTrigger(HADES_constants::kPT3)->GetIsFired() )
    {
        hIncorrectEvent->Fill(cPT3);
    }
}

void Selector::CheckTrackCuts(Int_t idx)
{
    DataTreeTrack* fTrack = fEvent->GetVertexTrack(idx);
    DataTreeTOFHit* fHit = fEvent->GetTOFHit(idx);
    Float_t fTof = fHit->GetTime();
    Float_t fLen = fHit->GetPathLength();
    //cout << len/tof/299.792458 << endl;
    if ( fTrack->GetDCAComponent(1) > 15 )
    {
        hIncorrectTracks->Fill(cDCA);
    }
    if ( fLen/fTof/299.792458 > 1 ) 
    {
        hIncorrectTracks->Fill(cBeta);
    }
    if ( fHit->GetPositionComponent(0) < -5 || fHit->GetPositionComponent(0) > 5 )
    {
        hIncorrectTracks->Fill(cTrackHitMatchX);
    }
    if ( fHit->GetPositionComponent(1) < -5 || fHit->GetPositionComponent(1) > 5 )
    {
        hIncorrectTracks->Fill(cTrackHitMatchY);
    }    
    if ( fTrack->GetChi2() > 100 )
    {
        hIncorrectTracks->Fill(cChi2);
    }    
}

void Selector::DrawStatistics()
{
    hIncorrectEvent->GetXaxis()->SetBinLabel(cVeretexPositionZ+1,"Vertex on Z");
    hIncorrectEvent->GetXaxis()->SetBinLabel(cVeretexPositionXY+1,"Vertex on X&Y");
    hIncorrectEvent->GetXaxis()->SetBinLabel(cVertexQuality+1,"Vertex Quality");
    hIncorrectEvent->GetXaxis()->SetBinLabel(cTriggerVertexClust+1,"Trigger VertexClust");
    hIncorrectEvent->GetXaxis()->SetBinLabel(cTriggerVertexCand+1,"Trigger VertexCand");
    hIncorrectEvent->GetXaxis()->SetBinLabel(cTriggerGoodStart+1,"Trigger GoodStart");
    hIncorrectEvent->GetXaxis()->SetBinLabel(cTriggerNoPileUp+1,"Trigger NoPileUp");
    hIncorrectEvent->GetXaxis()->SetBinLabel(cTriggerGoodStartVeto+1,"Trigger GoodStartVeto");
    hIncorrectEvent->GetXaxis()->SetBinLabel(cTriggerGoodStartMeta+1,"Trigger GoodStartMeta");
    hIncorrectEvent->GetXaxis()->SetBinLabel(cTriggerNoVeto+1,"Trigger NoVeto");
    hIncorrectEvent->GetXaxis()->SetBinLabel(cPT2+1,"Trigger PT2");
    hIncorrectEvent->GetXaxis()->SetBinLabel(cPT3+1,"Trigger PT3");
    hIncorrectEvent->GetXaxis()->SetTitle("");

    TCanvas* canv = new TCanvas("canv","Events",4000,3000);
    canv->cd();
    gStyle->SetOptStat(0);
    hIncorrectEvent->SetLineWidth(7);
    hIncorrectEvent->Draw();
    canv->SaveAs("../histograms/AmoutOfRejectedEvents.png");

    hIncorrectTracks->GetXaxis()->SetBinLabel(cDCA+1,"DCA");
    hIncorrectTracks->GetXaxis()->SetBinLabel(cTrackHitMatchX+1,"Track-Hit Match on X");
    hIncorrectTracks->GetXaxis()->SetBinLabel(cTrackHitMatchX+1,"Track-Hit Match on X");
    hIncorrectTracks->GetXaxis()->SetBinLabel(cTrackHitMatchY+1,"Track-Hit Match on Y");
    hIncorrectTracks->GetXaxis()->SetBinLabel(cChi2+1,"Chi2");
    hIncorrectTracks->GetXaxis()->SetBinLabel(cBeta+1,"Beta");
    hIncorrectTracks->GetXaxis()->SetTitle("");

    TCanvas* canv1 = new TCanvas("canv1","Tracks",4000,3000);
    canv1->cd();
    hIncorrectTracks->SetLineWidth(7);
    hIncorrectTracks->Draw();
    canv1->SaveAs("../histograms/AmoutOfRejectedTracks.png");
}

void Selector::SaveStatistics() 
{  
    this->DrawStatistics();
    TFile* fFile = new TFile("../histograms/SelectorStatisticsHistos.root","recreate");
    hRejectedEvents->Write();
    hRejectedTracks->Write();
    fFile->Close();
}