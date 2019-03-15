#include "Reader.h"

const double YCOR = 0.5*log(1.23*197+156.743) - 0.5*log(1.23*197-156.743);

Reader::Reader(char* cFileName)
{
    fChain = new TChain("DataTree");
    fChain->Add(cFileName);
    cout << fChain->GetEntries() << "events" << endl;
    fEvent = new DataTreeEvent;
    fChain->SetBranchAddress("DTEvent", &fEvent);
}

Reader::~Reader()
{
    delete fChain;
    delete fEvent;
}

void Reader::AddFile(char* cFileName)
{
    fChain->Add(cFileName);
    cout << fChain->GetEntries() << " " << "events" << endl;
}

DataTreeEvent* Reader::GetEvent(int idx)
{
    fChain->GetEntry(idx);
    return fEvent;
}
/*
void Reader::InitFlowHistos()
{
    cout << "Initialization of flow histograms" << endl;
    pFlowProfiles[meanQx] =         new TProfile("MeanQx vs Centrality",";centrality;mean Qx",10,0,50);
    pFlowProfiles[meanQy] =         new TProfile("MeanQy vs Centrality",";centrality;mean Qy",10,0,50);
    pFlowProfiles[resolution] =     new TProfile("Resolution vs Centrality",";centrality;R_{1}",10,0,50);

    pFlowProfiles[yMostCentral] =   new TProfile("v1 vs y, most",";rapidity;v_{1}",10,-0.8,0.8);
    pFlowProfiles[yMidCentral] =    new TProfile("v1 vs y, mid",";rapidity;v_{1}",10,-0.8,0.8);
    pFlowProfiles[yPeripherial] =   new TProfile("v1 vs y, per",";rapidity;v_{1}",10,-0.8,0.8);

    pFlowProfiles[ptMostCentral] =  new TProfile("v1 vs pt, most",";pt, [#frac{GeV}{c}];v_{1}",10,0.0,2.0);
    pFlowProfiles[ptMidCentral] =   new TProfile("v1 vs pt, mid",";pt, [#frac{GeV}{c}];v_{1}",10,0.0,2.0);
    pFlowProfiles[ptPeripherial] =  new TProfile("v1 vs pt, per",";pt, [#frac{GeV}{c}];v_{1}",10,0.0,2.0);

    vQvectorDistribution[QxNotRecentred] =       new TH1F("QxNotRecrntred",";Qx;conts",100,-1,1);
    vQvectorDistribution[QyNotRecentred] =       new TH1F("QyNotRecentred",";Qy;counts",100,-1,1);
    vQvectorDistribution[QxRecentred] =          new TH1F("QxRecentred",";Qx recentred;couts",100,-1,1);
    vQvectorDistribution[QyRecentred] =          new TH1F("QyRecrntred",";Qy recentred;counts",100,-1,1);
}

void Reader::FillCorrectionHistos()
{
    this->InitFlowHistos();
    cout << "Mean Q-vector and EP resolution as funtions of centrality histograms are building" << endl;
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
        if( fQ.GetComponent(0)!=fQ.GetComponent(0) || fQ.GetComponent(1)!=fQ.GetComponent(1) )
            continue;
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
        Float_t fRes = cos( PsiEP[0]-PsiEP[1] );
        pFlowProfiles[resolution]->Fill(fCentrality,fRes);
    }
}

void Reader::GetFlow(int iNumHarm)
{
    this->FillCorrectionHistos();
    cout << "v1 as function of centrality, rapidity and pt histograms are building" << endl;
    Float_t     fCentrality;
    Float_t     fPsiEP;
    Float_t*    fQCorection = new Float_t[2];
    Long64_t    lNEvents = fChain->GetEntries();
    DataTreeTrack* fTrack;
    for(int i=0; i<lNEvents; i++)
    {
        fChain->GetEntry(i);
        if(!selector.IsCorrectEvent(fEvent))
            continue;
        fCentrality = fEvent->GetCentrality();
        int bin = (fCentrality/5)+1;
        fQ.Estimate(fEvent);
        fQCorection[0] = pFlowProfiles[meanQx]->GetBinContent(bin);
        fQCorection[1] = pFlowProfiles[meanQy]->GetBinContent(bin);
        fQ.Recenter(fQCorection);
        fPsiEP = fQ.GetPsiEP();
        if(fPsiEP != fPsiEP)
            continue;
        Float_t fRes = pFlowProfiles[resolution]->GetBinContent(bin);
        fRes = sqrt(fRes);
        Int_t iNumberOfTracks = fEvent->GetNVertexTracks();
        for(int j=0;j<iNumberOfTracks;j++)
        {
            fTrack =    fEvent->GetVertexTrack(j);
            Float_t fPt =       fTrack->GetPt();
            Float_t fPhi =      fTrack->GetPhi();
            Float_t fRapidity = fTrack->GetRapidity() - YCOR;
            Float_t vn =        cos(fPhi-fPsiEP)/fRes;
            if ( fCentrality > 0 && fCentrality < 20 )
            {
                pFlowProfiles[yMostCentral]->Fill(fRapidity,vn);
                pFlowProfiles[ptMostCentral]->Fill(fPt,vn);
                continue;
            }
            if ( fCentrality >= 20 && fCentrality < 30 )
            {
                pFlowProfiles[yMidCentral]->Fill(fRapidity,vn);
                pFlowProfiles[ptMidCentral]->Fill(fPt,vn);
                continue;
            }
            if ( fCentrality >= 30 && fCentrality < 50 )
            {
                pFlowProfiles[yPeripherial]->Fill(fRapidity,vn);
                pFlowProfiles[ptPeripherial]->Fill(fPt,vn);
                continue;
            }
        }
    }
    this->SaveFlowStatistics();
}

void Reader::SaveFlowStatistics()
{
    TFile* fFile = new TFile("../histograms/FlowStatistics.root","recreate");
    fFile->cd();
    for(int i=0;i<NumOfFLowProfiles;i++)
    {
        pFlowProfiles[i]->Write();
    }
    fFile->Close();
}
*/
void Reader::BuildQAHistograms(TString sPicName)
{
    EventQA* fEventQA = new EventQA;
	TrackQA* fTrackQA[NumOfParticles];
	fTrackQA[all] = 		new TrackQA(-1);
	fTrackQA[electron] = 	new TrackQA(3);
	fTrackQA[positron] = 	new TrackQA(2);
	fTrackQA[pi_minus] = 	new TrackQA(9);
	fTrackQA[pi_plus] = 	new TrackQA(8);
	fTrackQA[proton] = 		new TrackQA(14);
	fTrackQA[deuteron] =  	new TrackQA(45);
	fTrackQA[helium3] = 	new TrackQA(49);
	fTrackQA[helium4] = 	new TrackQA(47);
    Long64_t lNEvents = fChain->GetEntries();
    for(int i=0; i<lNEvents; i++)
    {
        fChain->GetEntry(i);
        fEventQA->FillHistograms(fEvent);
		for(int j=0; j<NumOfParticles;j++) 
			fTrackQA[j]->FillHistograms(fEvent);
    }
    fEventQA->SaveHistograms(sPicName);
	for(int j=0; j<NumOfParticles;j++) 
		fTrackQA[j]->SaveHistograms(sPicName);
}

void Reader::BuildQvectorHistograms(TString sPicName)
{
	Long64_t lNEvents = fChain->GetEntries();
	Qvector* fQ =  new Qvector;
    Selector* fSelector = new Selector; 
	for(int i=0; i<lNEvents; i++)
    {
        fChain->GetEntry(i);
		if( !fSelector->IsCorrectEvent(fEvent) )
			continue;
		fQ->FillCorrections(fEvent);
    }
	
	for(int i=0; i<lNEvents; i++)
    {
        fChain->GetEntry(i);
		if( !fSelector->IsCorrectEvent(fEvent) )
			continue;
		fQ->Estimate(fEvent);
    }
	fQ->SaveHistograms(sPicName);
}