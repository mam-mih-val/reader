#include "EventQA.h"

EventQA::EventQA(DataTreeEvent* _fEvent, Selector* _selector, Centrality* _centrality)
{
	fEvent = _fEvent;
	fSelector = _selector;
	fCentrality = _centrality;
    this->InitHistograms();
}

void EventQA::InitHistograms()
{
    cout << "Initialization of Event QA histograms" << endl;
    vHisto1D[tracksMDC] =           new TH1F("tracksMDC",";tracks MDC;counts",100,0,100);
    vHisto1D[tracksMDC_selected] =  new TH1F("tracksMDC_selected",";selected tracks MDC;counts",100,0,100);
    vHisto1D[hitsTOF] =             new TH1F("hitsTOF",";hits in TOF+RPC;counts",250,0,250);
    vHisto1D[hitsTOF_selected] =    new TH1F("hitsTOF_selected",";hits in TOF+RPC selected;counts",250,0,250);
    vHisto1D[chargeFW] =            new TH1F("chargeFW",";charge in FW;counts",110,0,110);
    vHisto1D[chargeFW_selected] =   new TH1F("chargeFW_selected",";charge in FW selected;counts",110,0,110);
    vHisto1D[vertexZ] =             new TH1F("vertexZ",";vertex on Z;conts",100,-100,10);
    vHisto1D[vertexZ_selected] =    new TH1F("vertexZ_selected",";vertex on Z, selected;counts",100,-100,10);
    vHisto1D[hitsTOF_uncuted] =     new TH1F("hitsTOF_uncuted",";hits in TOF+RPC uncuted;counts",250,0,250);
    vHisto1D[hitsTOF_uncuted_selected]=new TH1F("hitsTOF_uncuted_selected",";hits in TOF+RPC uncuted selected;counts",250,0,250);
    vHisto1D[hitsTOF_matched] =     		new TH1F("hitsTOF_matched",";hits in TOF+RPC matched;counts",100,0,100);
    vHisto1D[hitsTOF_matched_selected]=		new TH1F("hitsTOF_matched_selected",";hits in TOF+RPC matched&selected;counts",100,0,100);
	int nbins = fCentrality->GetNumClasses();
	vHisto1D[histo_centrality] =			new TH1F("centrality",";centrality class;counts",20,0,100);
	vHisto1D[histo_centrality_selected] =	new TH1F("centrality_selected",";centrality class;counts",20,0,100);

    vHisto2D[tracks_hits] =         new TH2F("tracks&hits",";tracks MDC;hits TOF+RPC",100,0,100,100,0,250);
    vHisto2D[tracks_hits_selected]= new TH2F("tracks&hits_selected",";selected tracks MDC;selected hits TOF+RPC",100,0,100,100,0,250);
    vHisto2D[tracks_charge] =       new TH2F("tracks&charge",";tracks MDC;charge FW",100,0,100,100,0,100);
    vHisto2D[tracks_charge_selected]=new TH2F("tracks&charge_selected",";selected tracks MDC;selected charge FW",100,0,100,100,0,100);
    vHisto2D[hits_charge] =         new TH2F("hits&charge",";hits TOF+RPC;charge FW;",100,0,250,100,0,100);
    vHisto2D[hits_charge_selected] =new TH2F("hits&charge_selected",";selected hits TOF+RPC;selected charge FW",100,0,250,100,0,100);
    vHisto2D[vertexX_vertexY] =     new TH2F("vertexX&vertexY",";vertex on X;vertex on Y",100,-5,5,100,-5,5);
    vHisto2D[vertexX_vertexY_selected]=new TH2F("vertexX&vertexY_selected",";selected vertex on X;selected vertex on Y",100,-5,5,100,-5,5);
    vHisto2D[hitsFW_X_Y]=           new TH2F("hits in FW coordinates",";X, [mm];Y, [mm]",200,-1000,1000,200,-1000,1000);
    vHisto2D[hitsFW_X_Y_selected]=  new TH2F("selected hits in FW coordinates",";X, [mm];Y, [mm]",200,-1000,1000,200,-1000,1000);

	vProfile[hits_centrality] = 		new TProfile("NumOfHits_centrality",";centrality class;Hits TOF+RPC",nbins,0,nbins);
	vProfile[hits_centrality_selected]=	new TProfile("NumOfHits_centrality_selected",";centrality class;Hits TOF+RPC",nbins,0,nbins);
}

void EventQA::FillHistograms()
{
	vProfile[hits_centrality]->Fill( fCentrality->GetCentralityClass(), fEvent->GetCentralityEstimator(HADES_constants::kNhitsTOF_RPC_cut) );
    vHisto1D[tracksMDC]->Fill( fEvent->GetCentralityEstimator(HADES_constants::kNselectedTracks) );
    vHisto1D[hitsTOF]->Fill( fEvent->GetCentralityEstimator(HADES_constants::kNhitsTOF_RPC_cut) );
    vHisto1D[chargeFW]->Fill( fEvent->GetPSDEnergy() );
    vHisto1D[vertexZ]->Fill( fEvent->GetVertexPositionComponent(2) );
    vHisto1D[hitsTOF_uncuted]->Fill( fEvent->GetCentralityEstimator(HADES_constants::kNhitsTOF_RPC) );
    vHisto1D[hitsTOF_matched]->Fill( fEvent->GetNTOFHits() );
	vHisto1D[histo_centrality]->Fill( fEvent->GetCentrality(HADES_constants::kNhitsTOF_RPC_cut) );

    vHisto2D[tracks_hits]->Fill( fEvent->GetCentralityEstimator(HADES_constants::kNselectedTracks), fEvent->GetCentralityEstimator(HADES_constants::kNhitsTOF_RPC_cut) );
    vHisto2D[tracks_charge]->Fill( fEvent->GetCentralityEstimator(HADES_constants::kNselectedTracks), fEvent->GetPSDEnergy() );
    vHisto2D[hits_charge]->Fill( fEvent->GetCentralityEstimator(HADES_constants::kNhitsTOF_RPC_cut), fEvent->GetPSDEnergy() );
    vHisto2D[vertexX_vertexY]->Fill( fEvent->GetVertexPositionComponent(0), fEvent->GetVertexPositionComponent(1) );
    vHisto2D[hitsFW_X_Y]->Fill( fEvent->GetVertexPositionComponent(0), fEvent->GetVertexPositionComponent(1) );
    int iNPSDModules = fEvent->GetNPSDModules();
    DataTreePSDModule* fPSDModule;
    for(int j=0;j<iNPSDModules;j++)
    {
        fPSDModule = fEvent->GetPSDModule(j);
        vHisto2D[hitsFW_X_Y]->Fill(fPSDModule->GetPositionComponent(0),fPSDModule->GetPositionComponent(1),fPSDModule->GetEnergy());
    
    }

    if ( fSelector->IsCorrectEvent(HADES_constants::kPT2) ) 
    {
		vProfile[hits_centrality_selected]->Fill( fCentrality->GetCentralityClass(), fEvent->GetCentralityEstimator(HADES_constants::kNhitsTOF_RPC_cut) );
        vHisto1D[tracksMDC_selected]->Fill( fEvent->GetCentralityEstimator(HADES_constants::kNselectedTracks) );
        vHisto1D[hitsTOF_selected]->Fill( fEvent->GetCentralityEstimator(HADES_constants::kNhitsTOF_RPC_cut) );
        vHisto1D[chargeFW_selected]->Fill( fEvent->GetPSDEnergy() );
        vHisto1D[vertexZ_selected]->Fill( fEvent->GetVertexPositionComponent(2) );
        vHisto1D[hitsTOF_uncuted_selected]->Fill( fEvent->GetCentralityEstimator(HADES_constants::kNhitsTOF_RPC_cut) );
        vHisto1D[hitsTOF_matched_selected]->Fill( fEvent->GetNTOFHits() );
		vHisto1D[histo_centrality_selected]->Fill( fEvent->GetCentrality(HADES_constants::kNhitsTOF_RPC_cut) );

        vHisto2D[tracks_hits_selected]->Fill( fEvent->GetCentralityEstimator(HADES_constants::kNselectedTracks), fEvent->GetCentralityEstimator(HADES_constants::kNhitsTOF_RPC_cut) );
        vHisto2D[tracks_charge_selected]->Fill( fEvent->GetCentralityEstimator(HADES_constants::kNselectedTracks), fEvent->GetPSDEnergy() );
        vHisto2D[hits_charge_selected]->Fill( fEvent->GetCentralityEstimator(HADES_constants::kNhitsTOF_RPC_cut), fEvent->GetPSDEnergy() );
        vHisto2D[vertexX_vertexY_selected]->Fill( fEvent->GetVertexPositionComponent(0), fEvent->GetVertexPositionComponent(1) );
        vHisto2D[hitsFW_X_Y_selected]->Fill( fEvent->GetVertexPositionComponent(0), fEvent->GetVertexPositionComponent(1) );
        int iNPSDModules = fEvent->GetNPSDModules();
        DataTreePSDModule* fPSDModule;
        for(int j=0;j<iNPSDModules;j++)
        {
            fPSDModule = fEvent->GetPSDModule(j);
            if( fPSDModule->HasPassedCuts() )
                vHisto2D[hitsFW_X_Y_selected]->Fill(fPSDModule->GetPositionComponent(0),fPSDModule->GetPositionComponent(1),fPSDModule->GetEnergy());
        }   
    }
}

void EventQA::SaveHistograms(TString PicName)
{
    vCanvas[multiplicity_vertex] = new TCanvas("multiplicity&vertex","QA",4000,2500);
    vCanvas[multiplicity_vertex]->Divide(2,2,0.005,0.0001);
    vCanvas[multiplicity_charge] = new TCanvas("multiplicity&charge","QA",4500,2000);
    vCanvas[multiplicity_charge]->Divide(3,2,0.005,0.0001);
    vCanvas[vertex_charge] = new TCanvas("vertex&charge","QA",4000,2500);
    vCanvas[vertex_charge]->Divide(2,2,0.005,0.0001);
	vCanvas[centrality] = new TCanvas("centrality&multiplicity","QA",4000,2500);

    gStyle->SetOptStat(0);
    TLegend* legend = new TLegend(0.1,0.8,0.38,0.9);
    
    vCanvas[multiplicity_vertex]->cd(1);
    vHisto1D[tracksMDC]->SetLineWidth(5);
    vHisto1D[tracksMDC]->Draw();
    vHisto1D[tracksMDC_selected]->SetLineWidth(5);
    vHisto1D[tracksMDC_selected]->SetLineColor(2);
    vHisto1D[tracksMDC_selected]->Draw("same");
    legend->AddEntry(vHisto1D[tracksMDC],"All Events");
    legend->AddEntry(vHisto1D[tracksMDC_selected],"Selected Events");
    legend->Draw();

    vCanvas[multiplicity_vertex]->cd(2);
    vHisto1D[hitsTOF]->SetLineWidth(5);
    vHisto1D[hitsTOF]->Draw();
    vHisto1D[hitsTOF_selected]->SetLineWidth(5);
    vHisto1D[hitsTOF_selected]->SetLineColor(2);
    vHisto1D[hitsTOF_selected]->Draw("same");
    legend->Draw();

    vCanvas[multiplicity_vertex]->cd(3);
    vHisto1D[chargeFW]->SetLineWidth(5);
    vHisto1D[chargeFW]->Draw();
    vHisto1D[chargeFW_selected]->SetLineWidth(5);
    vHisto1D[chargeFW_selected]->SetLineColor(2);
    vHisto1D[chargeFW_selected]->Draw("same");
    legend->Draw();

    vCanvas[multiplicity_vertex]->cd(4);
    vHisto1D[vertexZ]->SetLineWidth(5);
    vHisto1D[vertexZ]->Draw();
    vHisto1D[vertexZ_selected]->SetLineWidth(5);
    vHisto1D[vertexZ_selected]->SetLineColor(2);
    vHisto1D[vertexZ_selected]->Draw("same");
    legend->Draw();

    vCanvas[multiplicity_charge]->cd(1)->SetLogz();
    vHisto2D[tracks_hits]->Draw("colz");

    vCanvas[multiplicity_charge]->cd(4)->SetLogz();
    vHisto2D[tracks_hits_selected]->Draw("colz");

    vCanvas[multiplicity_charge]->cd(2)->SetLogz();
    vHisto2D[tracks_charge]->Draw("colz");
    
    vCanvas[multiplicity_charge]->cd(5)->SetLogz();
    vHisto2D[tracks_charge_selected]->Draw("colz");

    vCanvas[multiplicity_charge]->cd(3)->SetLogz();
    vHisto2D[hits_charge]->Draw("colz");

    vCanvas[multiplicity_charge]->cd(6)->SetLogz();
    vHisto2D[hits_charge_selected]->Draw("colz");

    vCanvas[vertex_charge]->cd(1)->SetLogz();
    vHisto2D[vertexX_vertexY]->Draw("colz");

    vCanvas[vertex_charge]->cd(3)->SetLogz();
    vHisto2D[vertexX_vertexY_selected]->Draw("colz");

    vCanvas[vertex_charge]->cd(2)->SetLogz();
    vHisto2D[hitsFW_X_Y]->Draw("colz");

    vCanvas[vertex_charge]->cd(4)->SetLogz();
    vHisto2D[hitsFW_X_Y_selected]->Draw("colz");
	
	vCanvas[centrality]->Divide(2,1);
	vCanvas[centrality]->cd(1);
	//vHisto1D[histo_centrality]->Draw();
	vHisto1D[histo_centrality]->SetLineWidth(6);
	vHisto1D[histo_centrality]->Draw();
	vCanvas[centrality]->cd(2);
	vProfile[hits_centrality_selected]->SetMarkerStyle(20);
	vProfile[hits_centrality_selected]->SetMarkerSize(5);
	vProfile[hits_centrality_selected]->SetMarkerColor(1);
	vProfile[hits_centrality_selected]->SetLineColor(1);
	vProfile[hits_centrality_selected]->SetLineWidth(6);
	vProfile[hits_centrality_selected]->Draw();

    for(int i=0; i<NumCanvases; i++)
    {
        TString sPath = "../histograms/Event_"+PicName+Form("_%i.png",i);
        vCanvas[i]->SaveAs(sPath);
    }
}

void EventQA::SaveHisogramsToROOTFile(TString FileName)
{
    TString sPath = "../histograms/Event"+FileName+".root";
    TFile* file = new TFile(sPath,"recreate");
    file->cd();
    for(int i=0; i<Num1DHistos; i++)
    {
        vHisto1D[i]->Write();
    }
    for(int i=0; i<Num2DHistos; i++)
    {
        vHisto2D[i]->Write();
    }
    file->Close();
}
