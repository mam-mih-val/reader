#include "Reader.h"

const double YCOR = 0.5*log(1.23*197+156.743) - 0.5*log(1.23*197-156.743);
const double BETA = sqrt( 1.0 - 0.938*0.938/1.23/1.23 );

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

void Reader::DrawQA1DHistos(TString cPictureName)
{
    TCanvas* canv = new TCanvas("canv","QA",4500,2000);
    canv->Divide(3,2,0.005,0.0001);
    
    canv->cd(1)->SetLogy();
    vHisto1D[ptMDC]->SetLineWidth(5);
    vHisto1D[ptMDC]->Draw();
    
    canv->cd(4)->SetLogy();
    vHisto1D[ptMDC_selected]->SetLineWidth(5);
    vHisto1D[ptMDC_selected]->Draw();

    canv->cd(2)->SetLogy();
    vHisto1D[rapidityMDC]->SetLineWidth(5);
    vHisto1D[rapidityMDC]->Draw();
    
    canv->cd(5)->SetLogy();
    vHisto1D[rapidityMDC_selected]->SetLineWidth(5);
    vHisto1D[rapidityMDC_selected]->Draw();
    
    canv->cd(3);
    vHisto1D[phiMDC]->SetLineWidth(5);
    vHisto1D[phiMDC]->Draw();
    
    canv->cd(6);
    vHisto1D[phiMDC_selected]->SetLineWidth(5);
    vHisto1D[phiMDC_selected]->Draw();
    
    TString path = "../histograms/"+cPictureName+"_0.png";
    canv->SaveAs(path);

    TCanvas* canv1 = new TCanvas("canv1","QA1",4500,1500);
    canv1->Divide(3,1,0.005,0.0001);
    TLegend* legend = new TLegend(0.1,0.8,0.38,0.9);
    legend->AddEntry(vHisto1D[hitsTOF],"Unselected");
    legend->AddEntry(vHisto1D[hitsTOF_selected],"Selected");

    canv1->cd(1);
    vHisto1D[hitsTOF_selected]->SetLineWidth(5);
    vHisto1D[hitsTOF]->Draw();
    vHisto1D[hitsTOF_selected]->SetLineColor(1);
    vHisto1D[hitsTOF_selected]->SetLineWidth(5);
    vHisto1D[hitsTOF_selected]->Draw("same");
    legend->Draw();

    canv1->cd(2);
    vHisto1D[hitsTOF_uncuted]->SetLineWidth(5);
    vHisto1D[hitsTOF_uncuted]->Draw();
    vHisto1D[hitsTOF_uncuted_selected]->SetLineColor(1);
    vHisto1D[hitsTOF_uncuted_selected]->SetLineWidth(5);
    vHisto1D[hitsTOF_uncuted_selected]->Draw("same");
    legend->Draw();

    canv1->cd(3);
    vHisto1D[hitsTOF_matched]->SetLineWidth(5);
    vHisto1D[hitsTOF_matched]->Draw();
    vHisto1D[hitsTOF_matched_selected]->SetLineColor(1);
    vHisto1D[hitsTOF_matched_selected]->SetLineWidth(5);
    vHisto1D[hitsTOF_matched_selected]->Draw("same");
    legend->Draw();

    path = "../histograms/"+cPictureName+"_1.png";
    canv1->SaveAs(path);


}

void Reader::DrawQA2DHistos(TString cPictureName)
{
    TCanvas* canv = new TCanvas("canv","QA",4000,2500);
    TCanvas* canv1 = new TCanvas("canv1","QA",4000,2500);
    
    canv->Divide(2,2,0.005,0.0001);
    canv1->Divide(2,2,0.005,0.0001);
    
    canv->cd(1)->SetLogz();
    vHisto2D[hits_charge]->Draw("colz");
    
    canv->cd(3)->SetLogz();
    vHisto2D[hits_charge_selected]->Draw("colz");

    canv->cd(2)->SetLogz();
    vHisto2D[tracks_charge]->Draw("colz");
    
    canv->cd(4)->SetLogz();
    vHisto2D[tracks_charge_selected]->Draw("colz");
    
    canv1->cd(1)->SetLogz();
    vHisto2D[tracks_hits]->Draw("colz");
    
    canv1->cd(3)->SetLogz();
    vHisto2D[tracks_hits_selected]->Draw("colz");

    canv1->cd(2);
    vHisto2D[hitsFW_X_Y]->Draw("colz");
    
    canv1->cd(4);
    vHisto2D[hitsFW_X_Y_selected]->Draw("colz");

    TString path = "../histograms/"+cPictureName+"_0.png";
    TString path1 = "../histograms/"+cPictureName+"_1.png";
    
    canv->SaveAs(path);
    canv1->SaveAs(path1);

    TCanvas* canv2 = new TCanvas("canv2","QA",4000,2500);
    canv2->Divide(2,2,0.005,0.0001);
    
    canv2->cd(1);
    vHisto2D[phi_rapidity]->Draw("colz");

    canv2->cd(3);
    vHisto2D[phi_rapidity_selected]->Draw("colz");

    canv2->cd(2);
    vHisto2D[phi_pseudorapidity]->Draw("colz");

    canv2->cd(4);
    vHisto2D[phi_pseudorapidity_selected]->Draw("colz");
    path1 = "../histograms/"+cPictureName+"_2.png";
    canv2->SaveAs(path1);

    TCanvas* canv3 = new TCanvas("canv3","QA",4500,2000);
    canv3->Divide(3,2,0.005,0.0001);
    
    canv3->cd(1);
    vHisto2D[pt_rapidity]->Draw("colz");

    canv3->cd(4);
    vHisto2D[pt_rapidity_selected]->Draw("colz");

    canv3->cd(2);
    vHisto2D[pt_pseudorapidity]->Draw("colz");

    canv3->cd(5);
    vHisto2D[pt_pseudorapidity_selected]->Draw("colz");

    canv3->cd(3);
    vHisto2D[rapidity_pseudorapidity]->Draw("colz");

    canv3->cd(6);
    vHisto2D[rapidity_pseudorapidity_selected]->Draw("colz");

    path1 = "../histograms/"+cPictureName+"_3.png";
    canv3->SaveAs(path1);
}

void Reader::GetQualityAssurance(Int_t iPT)
{
    this->InitQAHistos();
    cout << "QA histograms are building" << endl;
    Long64_t lNEvents = fChain->GetEntries();
    Float_t fNHitsTOF;
    Float_t fNTracksMDC;
    Float_t fChargeFW;
    Float_t fVertexPosition[3];
    TVector3 b; b.SetXYZ(0,0,BETA/2);
    DataTreeTrack* fTrack;
    DataTreeTOFHit* fHit;
    for (long i=0; i<lNEvents; i++)
    {
        fChain->GetEntry(i);
        fNTracksMDC =    fEvent->GetCentralityEstimator(HADES_constants::kNselectedTracks);
        fNHitsTOF =      fEvent->GetCentralityEstimator(HADES_constants::kNhitsTOF_RPC_cut);
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
            float fMass2 = fMomentum.P()*fMomentum.P()*(1-fBeta*fBeta)/(fBeta*fBeta) * fHit->GetCharge();
            vHisto1D[ptMDC]->Fill(fTrack->GetPt());
            vHisto1D[betaTOF]->Fill(fBeta);
            vHisto1D[massTOF]->Fill(fMass2);
            fMomentum.Boost(b);
            vHisto1D[rapidityMDC]->Fill(fMomentum.Rapidity());
            vHisto1D[phiMDC]->Fill(fMomentum.Phi());
            vHisto2D[phi_rapidity]->Fill(fMomentum.Rapidity(),fMomentum.Phi());
            vHisto2D[phi_pt]->Fill(fMomentum.Pt(),fMomentum.Phi());
            vHisto2D[pt_rapidity]->Fill(fMomentum.Rapidity(),fMomentum.Pt());
            vHisto2D[phi_pseudorapidity]->Fill(fMomentum.PseudoRapidity(),fMomentum.Phi());
            vHisto2D[pt_pseudorapidity]->Fill(fMomentum.PseudoRapidity(),fMomentum.Pt());
            vHisto2D[rapidity_pseudorapidity]->Fill(fMomentum.PseudoRapidity(),fMomentum.Rapidity());
        }
        vHisto1D[hitsTOF_uncuted]->Fill( fEvent->GetCentralityEstimator(HADES_constants::kNhitsTOF_RPC) );
        vHisto1D[tracksMDC]->Fill(fNTracksMDC);
        vHisto1D[hitsTOF]->Fill(fNHitsTOF);
        vHisto1D[chargeFW]->Fill(fChargeFW);
        vHisto1D[vertexZ]->Fill(fVertexPosition[2]);
        vHisto2D[tracks_hits]->Fill(fNTracksMDC,fNHitsTOF);
        vHisto2D[tracks_charge]->Fill(fNTracksMDC,fChargeFW);
        vHisto2D[hits_charge]->Fill(fNHitsTOF,fChargeFW);
        vHisto2D[vertexX_vertexY]->Fill(fVertexPosition[0],fVertexPosition[1]);
        int iNPSDModules = fEvent->GetNPSDModules();
        DataTreePSDModule* fPSDModule;
        for(int j=0;j<iNPSDModules;j++)
        {
            fPSDModule = fEvent->GetPSDModule(j);
            vHisto2D[hitsFW_X_Y]->Fill(fPSDModule->GetPositionComponent(0),fPSDModule->GetPositionComponent(1),fPSDModule->GetEnergy());
        }
        if(!selector.IsCorrectEvent(fEvent, iPT))
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
            float fMass2 = fMomentum.P()*fMomentum.P()*(1-fBeta*fBeta)/(fBeta*fBeta) * fHit->GetCharge();
            vHisto1D[ptMDC_selected]->Fill(fTrack->GetPt());
            vHisto1D[betaTOF_selected]->Fill(fBeta);
            vHisto1D[massTOF_selected]->Fill(fMass2);
            fMomentum.Boost(b);
            vHisto1D[rapidityMDC_selected]->Fill(fMomentum.Rapidity());
            vHisto1D[phiMDC_selected]->Fill(fMomentum.Phi());
            vHisto2D[phi_rapidity_selected]->Fill(fMomentum.Rapidity(),fMomentum.Phi());
            vHisto2D[phi_pt_selected]->Fill(fMomentum.Pt(),fMomentum.Phi());
            vHisto2D[pt_rapidity_selected]->Fill(fMomentum.Rapidity(),fMomentum.Pt());
            vHisto2D[phi_pseudorapidity_selected]->Fill(fMomentum.PseudoRapidity(),fMomentum.Phi());
            vHisto2D[pt_pseudorapidity_selected]->Fill(fMomentum.PseudoRapidity(),fMomentum.Pt());
            vHisto2D[rapidity_pseudorapidity_selected]->Fill(fMomentum.PseudoRapidity(),fMomentum.Rapidity());
        }
        vHisto1D[hitsTOF_uncuted_selected]->Fill( fEvent->GetCentralityEstimator(HADES_constants::kNhitsTOF_RPC) );
        vHisto1D[vertexZ_selected]->Fill(fVertexPosition[2]);
        vHisto1D[tracksMDC_selected]->Fill(fNTracksMDC);
        vHisto1D[hitsTOF_selected]->Fill(fNHitsTOF);
        vHisto1D[chargeFW_selected]->Fill(fChargeFW);
        vHisto2D[tracks_hits_selected]->Fill(fNTracksMDC,fNHitsTOF);
        vHisto2D[tracks_charge_selected]->Fill(fNTracksMDC,fChargeFW);
        vHisto2D[hits_charge_selected]->Fill(fNHitsTOF,fChargeFW);
        vHisto2D[vertexX_vertexY_selected]->Fill(fVertexPosition[0],fVertexPosition[1]);
        iNPSDModules = fEvent->GetNPSDModules();
        for(int j=0;j<iNPSDModules;j++)
        {
            fPSDModule = fEvent->GetPSDModule(j);
            if ( fPSDModule->GetId() < 0 )
                continue;
            vHisto2D[hitsFW_X_Y_selected]->Fill(fPSDModule->GetPositionComponent(0),fPSDModule->GetPositionComponent(1),fPSDModule->GetEnergy());
        }
    }
    selector.SaveStatistics();
    this->SaveQAStatistics();
    return;
}

void Reader::GetPTChecking(Int_t iPT)
{
    this->InitQAHistos();
    Long64_t lNEvents = fChain->GetEntries();
    
    for (long i=0; i<lNEvents; i++)
    {
        fChain->GetEntry(i);
        vHisto1D[hitsTOF_uncuted]->Fill( fEvent->GetCentralityEstimator(HADES_constants::kNhitsTOF_RPC) );
        vHisto1D[hitsTOF]->Fill( fEvent->GetCentralityEstimator(HADES_constants::kNhitsTOF_RPC_cut) );
        vHisto1D[hitsTOF_matched]->Fill( fEvent->GetNTOFHits() );
    
        if( !fEvent->GetTrigger(iPT)->GetIsFired() ) 
            continue;
    
        vHisto1D[hitsTOF_uncuted_selected]->Fill( fEvent->GetCentralityEstimator(HADES_constants::kNhitsTOF_RPC) );
        vHisto1D[hitsTOF_selected]->Fill( fEvent->GetCentralityEstimator(HADES_constants::kNhitsTOF_RPC_cut) );
        vHisto1D[hitsTOF_matched_selected]->Fill( fEvent->GetNTOFHits() );
    }
}

DataTreeEvent* Reader::GetEvent(int idx)
{
    fChain->GetEntry(idx);
    return fEvent;
}

void Reader::InitQAHistos()
{
    cout << "Initialization of QA histograms" << endl;
    vHisto1D[tracksMDC] =           new TH1F("tracksMDC",";tracks MDC;counts",100,0,100);
    vHisto1D[tracksMDC_selected] =  new TH1F("tracksMDC_selected",";selected tracks MDC;counts",100,0,100);
    vHisto1D[hitsTOF] =             new TH1F("hitsTOF",";hits in TOF+RPC;counts",250,0,250);
    vHisto1D[hitsTOF_selected] =    new TH1F("hitsTOF_selected",";hits in TOF+RPC selected;counts",250,0,250);
    vHisto1D[chargeFW] =            new TH1F("chargeFW",";charge in FW;counts",100,0,8000);
    vHisto1D[chargeFW_selected] =   new TH1F("chargeFW_selected",";charge in FW selected;counts",100,0,8000);
    vHisto1D[vertexZ] =             new TH1F("vertexZ",";vertex on Z;conts",100,-100,10);
    vHisto1D[vertexZ_selected] =    new TH1F("vertexZ_selected",";vertex on Z, selected;counts",100,-100,10);
    vHisto1D[ptMDC] =               new TH1F("ptMDC",";pt, [#frac{GeV}{c}];counts",100,0,2.5);
    vHisto1D[ptMDC_selected] =      new TH1F("ptMDC_selected",";pt selected, [#frac{GeV}{c}];counts",100,0,2.5);
    vHisto1D[massTOF] =             new TH1F("massTOF",";m^{2}, [#frac{GeV}{c^{2}}];counts",100,-0.2,4);
    vHisto1D[massTOF_selected] =    new TH1F("massTOF_selected",";m^{2} selected, [#frac{GeV}{c^{2}}];counts",100,-0.2,4);
    vHisto1D[rapidityMDC] =         new TH1F("rapidityMDC",";rapidity, y;counts",100,-1,2);
    vHisto1D[rapidityMDC_selected]= new TH1F("rapidityMDC_selected",";rapidity selected, y;counts",100,-1,2);
    vHisto1D[rapidityMDC_recentred]=new TH1F("rapidityMDC_recentred",";rapidity recentred, y;counts",100,-1,2);
    vHisto1D[phiMDC] =              new TH1F("phiMDC",";phi;counts",100,-3.1415,3.1415);
    vHisto1D[phiMDC_selected] =     new TH1F("phiMDC_selected",";#phi selected;counts",100,-3.1415,3.1415);
    vHisto1D[betaTOF] =             new TH1F("betaTOF",";#beta;counts",100,0,1.2);
    vHisto1D[betaTOF_selected] =    new TH1F("betaTOFSelected",";#beta selected;counts",100,0,1.2);
    vHisto1D[hitsTOF_uncuted] =     new TH1F("hitsTOF_uncuted",";hits in TOF+RPC uncuted;counts",250,0,250);
    vHisto1D[hitsTOF_uncuted_selected]=new TH1F("hitsTOF_uncuted_selected",";hits in TOF+RPC uncuted selected;counts",250,0,250);
    vHisto1D[hitsTOF_matched] =     new TH1F("hitsTOF_matched",";hits in TOF+RPC matched;counts",100,0,100);
    vHisto1D[hitsTOF_matched_selected]=new TH1F("hitsTOF_matched_selected",";hits in TOF+RPC matched&selected;counts",100,0,100);

    vHisto2D[tracks_hits] =         new TH2F("tracks&hits",";tracks MDC;hits TOF+RPC",100,0,100,100,0,200);
    vHisto2D[tracks_hits_selected]= new TH2F("tracks&hits_selected",";selected tracks MDC;selected hits TOF+RPC",100,0,100,100,0,200);
    vHisto2D[tracks_charge] =       new TH2F("tracks&charge",";tracks MDC;charge FW",100,0,100,100,0,8000);
    vHisto2D[tracks_charge_selected]=new TH2F("tracks&charge_selected",";selected tracks MDC;selected charge FW",100,0,100,100,0,8000);
    vHisto2D[hits_charge] =         new TH2F("hits&charge",";hits TOF+RPC;charge FW;",100,0,200,100,0,8000);
    vHisto2D[hits_charge_selected] =new TH2F("hits&charge_selected",";selected hits TOF+RPC;selected charge FW",100,0,200,100,0,8000);
    vHisto2D[vertexX_vertexY] =     new TH2F("vertexX&vertexY",";vertex on X;vertex on Y",100,-5,5,100,-5,5);
    vHisto2D[vertexX_vertexY_selected]=new TH2F("vertexX&vertexY_selected",";selected vertex on X;selected vertex on Y",100,-5,5,100,-5,5);
    vHisto2D[hitsFW_X_Y]=           new TH2F("hits in FW coordinates",";X, [mm];Y, [mm]",49,-1000,1000,49,-1000,1000);
    vHisto2D[hitsFW_X_Y_selected]=  new TH2F("selected hits in FW coordinates",";X, [mm];Y, [mm]",49,-1000,1000,49,-1000,1000);
    vHisto2D[phi_rapidity]  =       new TH2F("phi&rapidity",";rapidity;phi, [rad]",100,-1,1,100,-3.14,3.14);
    vHisto2D[phi_rapidity_selected]=new TH2F("phi&rapidity_selected",";rapidity selected;phi selected, [rad]",100,-1,1,100,-3.14,3.14);
    vHisto2D[phi_pt]  =             new TH2F("phi&pt",";pt, [GeV/c];phi, [rad]",100,0,2,100,-3.14,3.14);
    vHisto2D[phi_pt_selected]  =    new TH2F("phi&pt_selected",";pt selected, [GeV/c];phi selected, [rad]",100,0,2,100,-3.14,3.14);
    vHisto2D[pt_rapidity]  =        new TH2F("pt&rapidity",";rapidity;pt, [GeV/c]",100,-1,1,100,0,2);
    vHisto2D[pt_rapidity_selected]= new TH2F("pt&rapidity_selected",";rapidity selected;pt selected, [GeV/c]",100,-1,1,100,0,2);
    vHisto2D[phi_pseudorapidity]  = new TH2F("phi&pseudorapidity",";pseudorapidity;phi, [rad]",100,-1,1,100,-3.14,3.14);
    vHisto2D[phi_pseudorapidity_selected]  = new TH2F("phi&pseudorapidity_selected",";pseudorapidity_selected;phi_selected, [rad]",100,-1,1,100,-3.14,3.14);
    vHisto2D[pt_pseudorapidity]  =  new TH2F("pt&pseudorapidity",";pseudorapidity;pt, [GeV/c]",100,-1,1,100,0,2);
    vHisto2D[pt_pseudorapidity_selected]=new TH2F("pt&pseudorapidity_selected",";pseudorapidity selected;pt selected, [GeV/c]",100,-1,1,100,0,2);
    vHisto2D[rapidity_pseudorapidity]=new TH2F("rapidity&pseudorapidity",";pseudorapidity;rapidity",100,-1,1,100,-1,1);
    vHisto2D[rapidity_pseudorapidity_selected]=new TH2F("rapidity&pseudorapidity_selected",";pseudorapidity selected;rapidity selected",100,-1,1,100,-1,1);
}

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

void Reader::SaveQAStatistics()
{
    TFile* fHistos = new TFile("../histograms/QualityAssuranceHistos.root","recreate");
    for(int i=0;i<Num1DHistos;i++)
    {
        vHisto1D[i]->Write();
    }
    for(int i=0;i<Num2DHistos;i++)
    {
        vHisto2D[i]->Write();
    }
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