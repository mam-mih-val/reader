#include "TrackQA.h"

const double BETA = sqrt( 1.0 - 0.938*0.938/1.23/1.23 );

TrackQA::TrackQA(DataTreeEvent* _fEvent, Selector* _fSelector, int _iPid)
{
	fEvent = _fEvent;
	fSelector = _fSelector;
	iPid = _iPid;
	this->InitHistograms();
}

void TrackQA::InitHistograms()
{
	cout << "Initialization of Track QA histograms" << endl;
	vHisto1D[ptMDC] =					new TH1F( Form("ptMDC_%i",iPid),";pt, [#frac{GeV}{c}];counts",100,0,2.5);
	vHisto1D[ptMDC_selected] =			new TH1F( Form("ptMDC_%i_selected",iPid),";pt selected, [#frac{GeV}{c}];counts",100,0,2.5);
	vHisto1D[massTOF] =			 		new TH1F( Form("massTOF_%i",iPid),";m^{2}, [#frac{GeV^{2}}{c^{4}}];counts",100,-0.2,18);
	vHisto1D[massTOF_selected] =		new TH1F( Form("massTOF_%i_selected",iPid),";m^{2} selected, [#frac{GeV^{2}}{c^{4}}];counts",100,-0.2,18);
	vHisto1D[rapidityMDC] =				new TH1F( Form("rapidityMDC_%i",iPid),";rapidity, y;counts",100,-1,2);
	vHisto1D[rapidityMDC_selected]= 	new TH1F( Form("rapidityMDC_%i_selected",iPid),";rapidity selected, y;counts",100,-1,2);
	vHisto1D[pseudorapidityMDC] =		new TH1F( Form("pseudorapidityMDC_%i",iPid),";pseudorapidity;counts",100,0,3);
	vHisto1D[pseudorapidityMDC_selected]=new TH1F(Form("pseudorapidityMDC_%i_selected",iPid),";pseudorapidity_selected;counts",100,0,3);
	vHisto1D[phiMDC] =			  		new TH1F( Form("phiMDC_%i",iPid),";phi;counts",100,-3.1415,3.1415);
	vHisto1D[phiMDC_selected] =	 		new TH1F( Form("phiMDC_%i_selected",iPid),";#phi selected;counts",100,-3.1415,3.1415);
	vHisto1D[betaTOF] =			 		new TH1F( Form("betaTOF_%i",iPid),";#beta;counts",100,0,1.2);
	vHisto1D[betaTOF_selected] =		new TH1F( Form("betaTOF_%i_selected",iPid),";#beta selected;counts",100,0,1.2);
	
	vHisto2D[phi_rapidity]  =	   				new TH2F( Form("phi&rapidity_%i",iPid),";rapidity;phi, [rad]",100,-1,2,100,-3.14,3.14);
	vHisto2D[phi_rapidity_selected]=			new TH2F( Form("phi&rapidity_%i_selected",iPid),";rapidity selected;phi selected, [rad]",100,-1,2,100,-3.14,3.14);
	vHisto2D[phi_pt]  =			 				new TH2F( Form("phi&pti_%i",iPid),";pt, [GeV/c];phi, [rad]",100,0,2.5,100,-3.14,3.14);
	vHisto2D[phi_pt_selected]  =				new TH2F( Form("phi&pt_%i_selected",iPid),";pt selected, [GeV/c];phi selected, [rad]",100,0,2.5,100,-3.14,3.14);
	vHisto2D[pt_rapidity]  =					new TH2F( Form("pt&rapidity_%i",iPid),";rapidity;pt, [GeV/c]",100,-1,2,100,0,2.5);
	vHisto2D[pt_rapidity_selected]= 			new TH2F( Form("pt&rapidity_%i_selected",iPid),";rapidity selected;pt selected, [GeV/c]",100,-1,2,100,0,2.5);
	vHisto2D[phi_pseudorapidity]  = 			new TH2F( Form("phi&pseudorapidity_%i",iPid),";pseudorapidity;phi, [rad]",100,0,3,100,-3.14,3.14);
	vHisto2D[phi_pseudorapidity_selected]  = 	new TH2F( Form("phi&pseudorapidity_%i_selected",iPid),";pseudorapidity_selected;phi_selected, [rad]",100,0,3,100,-3.14,3.14);
	vHisto2D[pt_pseudorapidity]  =  			new TH2F( Form("pt&pseudorapidity_%i",iPid),";pseudorapidity;pt, [GeV/c]",100,0,3,100,0,2.5);
	vHisto2D[pt_pseudorapidity_selected]=		new TH2F( Form("pt&pseudorapidity_%i_selected",iPid),";pseudorapidity selected;pt selected, [GeV/c]",100,0,3,100,0,2.5);
	vHisto2D[rapidity_pseudorapidity]=			new TH2F( Form("rapidity&pseudorapidity_%i",iPid),";pseudorapidity;rapidity",100,0,3,100,-1,2);
	vHisto2D[rapidity_pseudorapidity_selected]= new TH2F( Form("rapidity&pseudorapidity_%i_selected",iPid),";pseudorapidity selected;rapidity selected",100,0,3,100,-1,2);
	vHisto2D[dEdXTOF_p] =			  			new TH2F( Form("TOF_dEdx&p_%i",iPid),";p, [GeV/c]; #frac{dE}{dx} TOF",100,-2.5,4.0,100,0,20);
	vHisto2D[dEdXTOF_p_selected] =	 			new TH2F( Form("TOF_dEdx&p_%i_selected",iPid),";p selected, [GeV/c]; #frac{dE}{dx} TOF selected",100,-2.5,4.0,100,0,20);
	vHisto2D[dEdXMDC_p] =			  			new TH2F( Form("MDC_dEdx&p_%i",iPid),";p, [GeV/c]; #frac{dE}{dx} MDC",100,-2.5,4.0,100,0,30);
	vHisto2D[dEdXMDC_p_selected] =	 			new TH2F( Form("MDC_dEdx&p_%i_selected",iPid),";p selected, [GeV/c]; #frac{dE}{dx} MDC selected",100,-2.5,4.0,100,0,30);
	vHisto2D[beta_p] =	 						new TH2F( Form("beta&p_%i",iPid),";p, [GeV/c]; #beta",100,-2.5,4.0,100,0.,1.1);
	vHisto2D[beta_p_selected] =	 				new TH2F( Form("beta&p_%i_selected",iPid),";p selected, [GeV/c]; #beta selected",100,-2.5,4.0,100,0.,1.1);
	vHisto2D[DCA_X_Y] = 						new TH2F( Form("DCA_X_Y_%i",iPid),"; DCA_X, [mm]; DCA_Y, [mm]",100,-30,30,100,30,30);
	vHisto2D[DCA_X_Y_selected] = 				new TH2F( Form("DCA_X_Y_%i_selected",iPid),"; DCA_X, selected, [mm]; DCA_Y, selected, [mm]",100,-30,30,100,30,30);
}

void TrackQA::FillHistograms()
{
	int iNTracks = fEvent->GetNVertexTracks();
	DataTreeTrack* fTrack;
	DataTreeTOFHit* fHit;
	TVector3 b; b.SetXYZ(0,0,-BETA);
	for (int i=0;i<iNTracks;i++)
	{
		fTrack = fEvent->GetVertexTrack(i);
		if ( iPid != -1 )
		{
			if( iPid != fTrack->GetPdgId() )
			{
				continue;
			}
		}
		TLorentzVector fMomentum = fTrack->GetMomentum();
		fHit = fEvent->GetTOFHit(i);
		float fTof = fHit->GetTime();
		float fCharge = fHit->GetCharge(); 
		float fLen = fHit->GetPathLength();
		float fBeta = fLen/fTof/299.792458;
		float fMass2 = fMomentum.P()*fMomentum.P()*(1-fBeta*fBeta)/(fBeta*fBeta) * fHit->GetCharge();
		vHisto2D[DCA_X_Y]->Fill( fTrack->GetDCAComponent(0), fTrack->GetDCAComponent(1) );
		vHisto1D[ptMDC]->Fill(fTrack->GetPt());
		vHisto1D[betaTOF]->Fill(fBeta);
		vHisto1D[massTOF]->Fill(fMass2*fCharge);
		vHisto1D[pseudorapidityMDC]->Fill(fMomentum.PseudoRapidity());
		vHisto2D[phi_pseudorapidity]->Fill(fMomentum.PseudoRapidity(),fMomentum.Phi());
		vHisto2D[pt_pseudorapidity]->Fill(fMomentum.PseudoRapidity(),fMomentum.Pt());
		vHisto2D[dEdXMDC_p]->Fill( fMomentum.P()*fCharge,fTrack->GetdEdx(HADES_constants::kMDC_all) );
		vHisto2D[dEdXTOF_p]->Fill( fMomentum.P()*fCharge,fTrack->GetdEdx(HADES_constants::kMETA) );
		vHisto2D[beta_p]->Fill(fMomentum.P()*fCharge,fBeta);
		float fPR = fMomentum.PseudoRapidity();
		fMomentum.Boost(b);
		vHisto1D[rapidityMDC]->Fill(fMomentum.Rapidity());
		vHisto1D[phiMDC]->Fill(fMomentum.Phi());
		vHisto2D[phi_rapidity]->Fill(fMomentum.Rapidity(),fMomentum.Phi());
		vHisto2D[phi_pt]->Fill(fMomentum.Pt(),fMomentum.Phi());
		vHisto2D[pt_rapidity]->Fill(fMomentum.Rapidity(),fMomentum.Pt());
		vHisto2D[rapidity_pseudorapidity]->Fill(fPR,fMomentum.Rapidity());

		if ( !fSelector->IsCorrectEvent() || !fSelector->IsCorrectTrack(i) )
			continue;
		
		vHisto2D[DCA_X_Y_selected]->Fill( fTrack->GetDCAComponent(0), fTrack->GetDCAComponent(1) );
		vHisto1D[ptMDC_selected]->Fill(fTrack->GetPt());
		vHisto1D[betaTOF_selected]->Fill(fBeta);
		vHisto1D[massTOF_selected]->Fill(fMass2*fCharge);
		fMomentum = fTrack->GetMomentum();
		vHisto1D[pseudorapidityMDC_selected]->Fill(fMomentum.PseudoRapidity());
		vHisto2D[phi_pseudorapidity_selected]->Fill(fMomentum.PseudoRapidity(),fMomentum.Phi());
		vHisto2D[pt_pseudorapidity_selected]->Fill(fMomentum.PseudoRapidity(),fMomentum.Pt());
		vHisto2D[dEdXMDC_p_selected]->Fill( fMomentum.P()*fCharge,fTrack->GetdEdx(HADES_constants::kMDC_all) );
		vHisto2D[dEdXTOF_p_selected]->Fill( fMomentum.P()*fCharge,fTrack->GetdEdx(HADES_constants::kMETA) );
		vHisto2D[beta_p_selected]->Fill(fMomentum.P()*fCharge,fBeta);
		fMomentum.Boost(b);
		vHisto1D[rapidityMDC_selected]->Fill(fMomentum.Rapidity());
		vHisto1D[phiMDC_selected]->Fill(fMomentum.Phi());
		vHisto2D[phi_rapidity_selected]->Fill(fMomentum.Rapidity(),fMomentum.Phi());
		vHisto2D[phi_pt_selected]->Fill(fMomentum.Pt(),fMomentum.Phi());
		vHisto2D[pt_rapidity_selected]->Fill(fMomentum.Rapidity(),fMomentum.Pt());
		vHisto2D[rapidity_pseudorapidity_selected]->Fill(fPR,fMomentum.Rapidity());
	}
}

void TrackQA::SaveHistograms(TString PicName)
{
	
	for(int i=0; i<NumCanvases; i++)
	{
		vCanvas[i] = new TCanvas(Form("Canvas_%i_%i",i,iPid),"QA",4500,2000);
		vCanvas[i]->Divide(3,2,0.005,0.0001);
	}

	gStyle->SetOptStat(0);
	TLegend* legend = new TLegend(0.1,0.8,0.38,0.9);

	vCanvas[OneDimHist]->cd(1)->SetLogy();
	vHisto1D[ptMDC]->SetLineWidth(5);
	vHisto1D[ptMDC]->Draw();
	vHisto1D[ptMDC_selected]->SetLineWidth(5);
	vHisto1D[ptMDC_selected]->SetLineColor(2);
	vHisto1D[ptMDC_selected]->Draw("same");
	legend->AddEntry(vHisto1D[ptMDC],"All Tracks");
	legend->AddEntry(vHisto1D[ptMDC_selected],"Selected Tracks");
	legend->Draw();

	vCanvas[OneDimHist]->cd(2)->SetLogy();
	vHisto1D[rapidityMDC]->SetLineWidth(5);
	vHisto1D[rapidityMDC]->Draw();
	vHisto1D[rapidityMDC_selected]->SetLineWidth(5);
	vHisto1D[rapidityMDC_selected]->SetLineColor(2);
	vHisto1D[rapidityMDC_selected]->Draw("same");
	legend->Draw();

	vCanvas[OneDimHist]->cd(3)->SetLogy();
	vHisto1D[pseudorapidityMDC]->SetLineWidth(5);
	vHisto1D[pseudorapidityMDC]->Draw();
	vHisto1D[pseudorapidityMDC_selected]->SetLineWidth(5);
	vHisto1D[pseudorapidityMDC_selected]->SetLineColor(2);
	vHisto1D[pseudorapidityMDC_selected]->Draw("same");
	legend->Draw();

	vCanvas[OneDimHist]->cd(4);
	vHisto1D[massTOF]->SetLineWidth(5);
	vHisto1D[massTOF]->Draw();
	vHisto1D[massTOF_selected]->SetLineWidth(5);
	vHisto1D[massTOF_selected]->SetLineColor(2);
	vHisto1D[massTOF_selected]->Draw("same");
	legend->Draw();

	vCanvas[OneDimHist]->cd(5);
	vHisto1D[phiMDC]->SetLineWidth(5);
	vHisto1D[phiMDC]->Draw();
	vHisto1D[phiMDC_selected]->SetLineWidth(5);
	vHisto1D[phiMDC_selected]->SetLineColor(2);
	vHisto1D[phiMDC_selected]->Draw("same");
	legend->Draw();

	vCanvas[OneDimHist]->cd(6)->SetLogy();
	vHisto1D[betaTOF]->SetLineWidth(5);
	vHisto1D[betaTOF]->Draw();
	vHisto1D[betaTOF_selected]->SetLineWidth(5);
	vHisto1D[betaTOF_selected]->SetLineColor(2);
	vHisto1D[betaTOF_selected]->Draw("same");
	legend->Draw();
//	***
	vCanvas[phi_kinematics]->cd(1)->SetLogz();
	vHisto2D[phi_rapidity]->Draw("colz");
	vCanvas[phi_kinematics]->cd(2)->SetLogz();
	vHisto2D[phi_pseudorapidity]->Draw("colz");
	vCanvas[phi_kinematics]->cd(3)->SetLogz();
	vHisto2D[phi_pt]->Draw("colz");
	vCanvas[phi_kinematics]->cd(4)->SetLogz();
	vHisto2D[phi_rapidity_selected]->Draw("colz");
	vCanvas[phi_kinematics]->cd(5)->SetLogz();
	vHisto2D[phi_pseudorapidity_selected]->Draw("colz");
	vCanvas[phi_kinematics]->cd(6)->SetLogz();
	vHisto2D[phi_pt_selected]->Draw("colz");
//	***
	vCanvas[pt_kinematics]->cd(1)->SetLogz();
	vHisto2D[pt_rapidity]->Draw("colz");
	vCanvas[pt_kinematics]->cd(2)->SetLogz();
	vHisto2D[pt_pseudorapidity]->Draw("colz");
	vCanvas[pt_kinematics]->cd(3)->SetLogz();
	vHisto2D[rapidity_pseudorapidity]->Draw("colz");
	vCanvas[pt_kinematics]->cd(4)->SetLogz();
	vHisto2D[pt_rapidity_selected]->Draw("colz");
	vCanvas[pt_kinematics]->cd(5)->SetLogz();
	vHisto2D[pt_pseudorapidity_selected]->Draw("colz");
	vCanvas[pt_kinematics]->cd(6)->SetLogz();
	vHisto2D[rapidity_pseudorapidity_selected]->Draw("colz");
//	***
	vCanvas[mass_qa]->cd(1)->SetLogz();
	vHisto2D[dEdXTOF_p]->Draw("colz");
	vCanvas[mass_qa]->cd(2)->SetLogz();
	vHisto2D[dEdXMDC_p]->Draw("colz");
	vCanvas[mass_qa]->cd(3)->SetLogz();
	vHisto2D[beta_p]->Draw("colz");
	vCanvas[mass_qa]->cd(4)->SetLogz();
	vHisto2D[dEdXTOF_p_selected]->Draw("colz");
	vCanvas[mass_qa]->cd(5)->SetLogz();
	vHisto2D[dEdXMDC_p_selected]->Draw("colz");
	vCanvas[mass_qa]->cd(6)->SetLogz();
	vHisto2D[beta_p_selected]->Draw("colz");
// ***
	for(int i=0; i<NumCanvases; i++)
    {
        TString sPath = "../histograms/Track_"+PicName+"_Pid_"+Form("%i_%i.png",iPid,i);
        vCanvas[i]->SaveAs(sPath);
    }
}

void TrackQA::SaveHisogramsToROOTFile(TString FileName)
{
    TString sPath = FileName+".root";
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