#include "Reader.h"

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
    Long64_t lNEvents = fChain->GetEntries();
    Selector selector;
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
        if(!selector.IsCorrect(fEvent))
            continue;

        for (int j=0;j<iNTracks;j++)
        {
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
    return;
}

DataTreeEvent* Reader::GetEvent(int idx)
{
    fChain->GetEntry(idx);
    return fEvent;
}

void Reader::InitHistos()
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
    vHisto1D[rapidityMDC] =         new TH1F("rapidityMDC","",100,-2,4);
    vHisto1D[rapidityMDC_selected]= new TH1F("rapidityMDC_selected","",100,-2,4);
    vHisto1D[rapidityMDC_recentred]=new TH1F("rapidityMDC_recentred","",100,-2,4);
    vHisto1D[phiMDC] =              new TH1F("phiMDC","",100,-3.1415,3.1415);
    vHisto1D[phiMDC_selected] =     new TH1F("phiMDC_selected","",100,-3.1415,3.1415);
    vHisto1D[betaTOF] =             new TH1F("betaTOF","",100,0,1.2);
    vHisto1D[betaTOF_selected] =     new TH1F("betaTOFSelected","",100,0,1.2);

    vHisto2D[tracks_hits] =         new TH2F("tracks&hits","",100,0,100,100,0,200);
    vHisto2D[tracks_hits_selected]= new TH2F("tracks&hits_selected","",100,0,100,100,0,200);
    vHisto2D[tracks_charge] =       new TH2F("tracks&charge","",100,0,100,100,0,8000);
    vHisto2D[tracks_charge_selected]=new TH2F("tracks&charge_selected","",100,0,100,0,8000);
    vHisto2D[hits_charge] =         new TH2F("hits&charge","",100,0,200,100,0,8000);
    vHisto2D[hits_charge_selected] =new TH2F("hits&charge_selected","",100,0,200,100,0,8000);
    vHisto2D[vertexX_vertexY] =     new TH2F("vertexX&vertexY","",100,0,5,100,0,5);
    vHisto2D[vertexX_vertexY_selected]=new TH2F("vertexX&vertexY_selected","",100,0,5,100,0,5);
}

void Reader::SaveStatistics()
{
    TFile* fHistos = new TFile("../histograms/QualiryAccuranceHistos.root","recreate");
    for(int i=0;i<Num1DHistos;i++)
    {
        vHisto1D[i]->Write();
    }
    for(int i=0;i<Num2DHistos;i++)
    {
        vHisto2D[i]->Write();
    }
}