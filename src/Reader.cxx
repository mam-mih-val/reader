#include "Reader.h"

const double YCOR = 0.5*log(1.23*197+156.743) - 0.5*log(1.23*197-156.743);

Reader::Reader(char* cFileName)
{
    fChain = new TChain("DataTree");
    fChain->Add(cFileName);
    cout << fChain->GetEntries() << endl;
    fEvent = new DataTreeEvent;
    DTEvent = (TBranch*) fChain->GetBranch("DTEvent");
    DTEvent->SetAddress(&fEvent);
}

Reader::~Reader()
{
    delete fChain;
    delete fEvent;
}

void Reader::GetQualityAccurance()
{
    this->InitQAHistos();
    Long64_t lNEvents = fChain->GetEntries();
    Float_t fNHitsTOF;
    Float_t fNTracksMDC;
    Float_t fChargeFW;
    Float_t fVertexPosition[3];
    DataTreeTrack* fTrack;
    DataTreeTOFHit* fHit;
    for (long i=0; i<lNEvents; i++)
    {
        fChain->GetEntry(i);
        fNTracksMDC =    fEvent->GetCentralityEstimator(HADES_constants::kNtracks);
        fNHitsTOF =      fEvent->GetCentralityEstimator(HADES_constants::kNhitsTOF_RPC);
        fChargeFW =      fEvent->GetPSDEnergy();
        for(int j=0;j<3;j++)
            fVertexPosition[j] = fEvent->GetVertexPositionComponent(j);
        int iNTracks = fEvent->GetNVertexTracks();
        for (int j=0;j<iNTracks;j++)
        {
            fTrack = fEvent->GetVertexTrack(j);
            TLorentzVector fMomentum = fTrack->GetMomentum();
            fHit = fEvent->GetTOFHit(j);
            float fTof = fHit->GetTime(); 
            float fLen = fHit->GetPathLength();
            float fBeta = fLen/fTof/299.792458;
            float fMass2 = fMomentum.P()*fMomentum.P()*(1-fBeta*fBeta)/(fBeta*fBeta);
            vHisto1D[ptMDC]->Fill(fTrack->GetPt());
            vHisto1D[betaTOF]->Fill(fBeta);
            vHisto1D[massTOF]->Fill(fMass2);
            vHisto1D[rapidityMDC]->Fill(fMomentum.Rapidity());
            vHisto1D[phiMDC]->Fill(fMomentum.Phi());
        }
        vHisto1D[tracksMDC]->Fill(fNTracksMDC);
        vHisto1D[hitsTOF]->Fill(fNHitsTOF);
        vHisto1D[chargeFW]->Fill(fChargeFW);
        vHisto1D[vertexZ]->Fill(fVertexPosition[2]);
        vHisto2D[tracks_hits]->Fill(fNTracksMDC,fNHitsTOF);
        vHisto2D[tracks_charge]->Fill(fNTracksMDC,fChargeFW);
        vHisto2D[hits_charge]->Fill(fNHitsTOF,fChargeFW);
        vHisto2D[vertexX_vertexY]->Fill(fVertexPosition[0],fVertexPosition[1]);
        if(!selector.IsCorrectEvent(fEvent))
            continue;

        for (int j=0;j<iNTracks;j++)
        {
            
            if(!selector.IsCorrectTrack(j))
                continue;
            fTrack = fEvent->GetVertexTrack(j);
            TLorentzVector fMomentum = fTrack->GetMomentum();
            fHit = fEvent->GetTOFHit(j);
            float fTof = fHit->GetTime(); 
            float fLen = fHit->GetPathLength();
            float fBeta = fLen/fTof/299.792458;
            float fMass2 = fMomentum.P()*fMomentum.P()*(1-fBeta*fBeta)/(fBeta*fBeta);
            vHisto1D[ptMDC_selected]->Fill(fTrack->GetPt());
            vHisto1D[betaTOF_selected]->Fill(fBeta);
            vHisto1D[massTOF_selected]->Fill(fMass2);
            vHisto1D[rapidityMDC_selected]->Fill(fMomentum.Rapidity());
            vHisto1D[rapidityMDC_recentred]->Fill(fMomentum.Rapidity()-YCOR);
            vHisto1D[phiMDC_selected]->Fill(fMomentum.Phi());
        }
        vHisto1D[vertexZ_selected]->Fill(fVertexPosition[2]);
        vHisto1D[tracksMDC_selected]->Fill(fNTracksMDC);
        vHisto1D[hitsTOF_selected]->Fill(fNHitsTOF);
        vHisto1D[chargeFW_selected]->Fill(fChargeFW);
        vHisto2D[tracks_hits_selected]->Fill(fNTracksMDC,fNHitsTOF);
        vHisto2D[tracks_charge_selected]->Fill(fNTracksMDC,fChargeFW);
        vHisto2D[hits_charge_selected]->Fill(fNHitsTOF,fChargeFW);
        vHisto2D[vertexX_vertexY_selected]->Fill(fVertexPosition[0],fVertexPosition[1]);
    }
    selector.SaveStatistics();
    this->SaveQAStatistics();
    return;
}

DataTreeEvent* Reader::GetEvent(int idx)
{
    fChain->GetEntry(idx);
    return fEvent;
}

void Reader::InitQAHistos()
{
    vHisto1D[tracksMDC] =           new TH1F("tracksMDC","",100,0,100);
    vHisto1D[tracksMDC_selected] =  new TH1F("tracksMDC_selected","",100,0,100);
    vHisto1D[hitsTOF] =             new TH1F("hitsTOF","",200,0,200);
    vHisto1D[hitsTOF_selected] =    new TH1F("hitsTOF_selected","",100,0,200);
    vHisto1D[chargeFW] =            new TH1F("chargeFW","",100,0,8000);
    vHisto1D[chargeFW_selected] =   new TH1F("chargeFW_selected","",100,0,8000);
    vHisto1D[vertexZ] =             new TH1F("vertexZ","",100,-100,10);
    vHisto1D[vertexZ_selected] =    new TH1F("vertexZ_selected","",100,-100,10);
    vHisto1D[ptMDC] =               new TH1F("ptMDC","",100,0,2.5);
    vHisto1D[ptMDC_selected] =      new TH1F("ptMDC_selected","",100,0,2.5);
    vHisto1D[massTOF] =             new TH1F("massTOF","",100,0,4);
    vHisto1D[massTOF_selected] =    new TH1F("massTOF_selected","",100,0,4);
    vHisto1D[rapidityMDC] =         new TH1F("rapidityMDC","",100,-1,2);
    vHisto1D[rapidityMDC_selected]= new TH1F("rapidityMDC_selected","",100,-1,2);
    vHisto1D[rapidityMDC_recentred]=new TH1F("rapidityMDC_recentred","",100,-1,2);
    vHisto1D[phiMDC] =              new TH1F("phiMDC","",100,-3.1415,3.1415);
    vHisto1D[phiMDC_selected] =     new TH1F("phiMDC_selected","",100,-3.1415,3.1415);
    vHisto1D[betaTOF] =             new TH1F("betaTOF","",100,0,1.2);
    vHisto1D[betaTOF_selected] =     new TH1F("betaTOFSelected","",100,0,1.2);

    vHisto2D[tracks_hits] =         new TH2F("tracks&hits","",100,0,100,100,0,200);
    vHisto2D[tracks_hits_selected]= new TH2F("tracks&hits_selected","",100,0,100,100,0,200);
    vHisto2D[tracks_charge] =       new TH2F("tracks&charge","",100,0,100,100,0,8000);
    vHisto2D[tracks_charge_selected]=new TH2F("tracks&charge_selected","",100,0,100,100,0,8000);
    vHisto2D[hits_charge] =         new TH2F("hits&charge","",100,0,200,100,0,8000);
    vHisto2D[hits_charge_selected] =new TH2F("hits&charge_selected","",100,0,200,100,0,8000);
    vHisto2D[vertexX_vertexY] =     new TH2F("vertexX&vertexY","",100,-5,5,100,-5,5);
    vHisto2D[vertexX_vertexY_selected]=new TH2F("vertexX&vertexY_selected","",100,-5,5,100,-5,5);
}

void Reader::GetFlow(int iNumHarm=1)
{
    this->InitFlowHistos();
    Qvector fQ;
    Float_t fCentrality;
    Long64_t lNEvents = fChain->GetEntries();
    // Mean Q as function of Centrality loop
    for (int i=0;i<lNEvents;i++)
    {
        fChain->GetEntry(i);
        if( !selector.IsCorrectEvent(fEvent) )
            continue;
        fCentrality = fEvent->GetCentrality();
        fQ.Estimate(fEvent);
        pFlowProfiles[meanQx]->Fill(fCentrality,fQ.GetComponent(0));
        pFlowProfiles[meanQy]->Fill(fCentrality,fQ.GetComponent(1));
        vQvectorDistribution[QxNotRecentred]->Fill(fQ.GetComponent(0));
        vQvectorDistribution[QyNotRecentred]->Fill(fQ.GetComponent(1));
    }
    // Event resolution as function of centrality
    Float_t PsiEP[2];
    Float_t* fQCorection = new Float_t[2];
    for (int i=0;i<lNEvents;i++)
    {
        fChain->GetEntry(i);
        if(!selector.IsCorrectEvent(fEvent))
            continue;
        fCentrality = fEvent->GetCentrality();
        int bin = (fCentrality/5)+1;
        fQ.Estimate(fEvent,1);
        fQCorection[0] = pFlowProfiles[meanQx]->GetBinContent(bin);
        fQCorection[1] = pFlowProfiles[meanQy]->GetBinContent(bin);
        fQ.Recenter(fQCorection);
        PsiEP[0] = fQ.GetPsiEP(0);
        PsiEP[1] = fQ.GetPsiEP(1);
        if ( PsiEP[0] != PsiEP[0] || PsiEP[1] != PsiEP[1] )
            continue;
        Float_t fRes = 2*cos(PsiEP[0]-PsiEP[1]);
        pFlowProfiles[resolution]->Fill(fCentrality,fRes);
    }

    this->SaveFlowStatistics();
}

void Reader::SaveQAStatistics()
{
    TFile* fHistos = new TFile("../histograms/QualityAccuranceHistos.root","recreate");
    for(int i=0;i<Num1DHistos;i++)
    {
        vHisto1D[i]->Write();
    }
    for(int i=0;i<Num2DHistos;i++)
    {
        vHisto2D[i]->Write();
    }
}

void Reader::InitFlowHistos()
{
    pFlowProfiles[meanQx] =         new TProfile("MeanQx vs Centrality","",10,0,50);
    pFlowProfiles[meanQy] =         new TProfile("MeanQy vs Centrality","",10,0,50);
    pFlowProfiles[resolution] =     new TProfile("Resolution vs Centrality","",10,0,50);

    pFlowProfiles[yMostCentral] =   new TProfile("v1 vs y, most","",10,-0.8,0.8);
    pFlowProfiles[yMidCentral] =    new TProfile("v1 vs y, mid","",10,-0.8,0.8);
    pFlowProfiles[yPeripherial] =   new TProfile("v1 vs y, per","",10,-0.8,0.8);

    pFlowProfiles[ptMostCentral] =  new TProfile("v1 vs pt, most","",10,0.0,2.0);
    pFlowProfiles[ptMidCentral] =   new TProfile("v1 vs pt, mid","",10,0.0,2.0);
    pFlowProfiles[ptPeripherial] =  new TProfile("v1 vs pt, per","",10,0.0,2.0);

    vQvectorDistribution[QxNotRecentred] =       new TH1F("QxNotRecrntred","",100,-1,1);
    vQvectorDistribution[QyNotRecentred] =       new TH1F("QyNotRecentred","",100,-1,1);
    vQvectorDistribution[QxRecentred] =          new TH1F("QxRecentred","",100,-1,1);
    vQvectorDistribution[QyRecentred] =          new TH1F("QyRecrntred","",100,-1,1);
}

void Reader::SaveFlowStatistics()
{
    TFile* fFile = new TFile("../histograms/FlowStatistics.root","recreate");
    for(int i=0;i<NumOfFLowProfiles;i++)
    {
        pFlowProfiles[i]->Write();
    }
    fFile->Close();
}