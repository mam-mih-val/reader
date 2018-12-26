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
    for (long i=0; i<lNEvents; i++)
    {
        fChain->GetEntry(i);
        fNTracksMDC =    fEvent->GetCentralityEstimator(HADES_constants::kNtracks);
        fNHitsTOF =      fEvent->GetCentralityEstimator(HADES_constants::kNhitsTOF_RPC);
        fChargeFW =      fEvent->GetPSDEnergy();
        vHisto1D[tracksMDC]->Fill(fNTracksMDC);
        vHisto1D[hitsTOF]->Fill(fNHitsTOF);
        vHisto1D[chargeFW]->Fill(fChargeFW);
        vHisto2D[tracks_hits]->Fill(fNTracksMDC,fNHitsTOF);
        vHisto2D[tracks_charge]->Fill(fNTracksMDC,fChargeFW);
        vHisto2D[hits_charge]->Fill(fNHitsTOF,fChargeFW);
        if(!selector.IsCorrect(fEvent))
            continue;
        vHisto1D[tracksMDC_selected]->Fill(fNTracksMDC);
        vHisto1D[hitsTOF_selected]->Fill(fNHitsTOF);
        vHisto1D[chargeFW_selected]->Fill(fChargeFW);
        vHisto2D[tracks_hits_selected]->Fill(fNTracksMDC,fNHitsTOF);
        vHisto2D[tracks_charge_selected]->Fill(fNTracksMDC,fChargeFW);
        vHisto2D[hits_charge_selected]->Fill(fNHitsTOF,fChargeFW);
    }
    selector.SaveStatistics();
    return;
}

DataTreeEvent* Reader::GetEvent(int idx)
{
    fChain->GetEntry(idx);
    return fEvent;
}

void Reader::SaveStatistics()
{
    TFile* fHistos = new TFile("histograms.root","recreate");
    for(int i=0;i<Num1DHistos;i++)
    {
        vHisto1D[i]->Write();
    }
    for(int i=0;i<Num2DHistos;i++)
    {
        vHisto2D[i]->Write();
    }
}

void Reader::InitHistos()
{
    vHisto1D[tracksMDC] =           new TH1F("tracksMDC","",100,0,100);
    vHisto1D[tracksMDC_selected] =  new TH1F("tracksMDC_selected","",100,0,100);
    vHisto1D[hitsTOF] =             new TH1F("hitsTOF","",200,0,200);
    vHisto1D[hitsTOF_selected] =    new TH1F("hitsTOF_selected","",100,0,200);
    vHisto1D[chargeFW] =            new TH1F("chargeFW","",100,0,8000);
    vHisto1D[chargeFW_selected] =   new TH1F("chargeFW_selected","",100,0,8000);

    vHisto2D[tracks_hits] =         new TH2F("tracks&hits","",100,0,100,100,0,200);
    vHisto2D[tracks_hits_selected]= new TH2F("tracks&hits_selected","",100,0,100,100,0,200);
    vHisto2D[tracks_charge] =       new TH2F("tracks&charge","",100,0,100,100,0,8000);
    vHisto2D[tracks_charge_selected]=new TH2F("tracks&charge_selected","",100,0,100,0,8000);
    vHisto2D[hits_charge] =         new TH2F("hits&charge","",100,0,200,100,0,8000);
    vHisto2D[hits_charge_selected] =new TH2F("hits&charge_selected","",100,0,200,100,0,8000);
}