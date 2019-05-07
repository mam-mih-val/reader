#include "Flow3SE.h"

Flow3SE::Flow3SE(DataTreeEvent* _event, Centrality* _centrality, Qvector3SE* _Qvector, Selector* _selector, int _pid)
{
	fEvent = _event;
	fCentrality = _centrality;
	fQvector = _Qvector;
	fSelector = _selector;
	fPid = _pid;
	cout << "All was linked successfully" << endl;
	for(int i=0; i<3; i++)
		fUvector.push_back( TVector2{0.,0.} );
	this->InitializeHistograms();
}

Flow3SE::~Flow3SE()
{
	for(auto histo:xRapidity.at(0))
		delete histo;
	for(auto histo:xRapidity.at(1))
		delete histo;
	for(auto histo:xRapidity.at(2))
		delete histo;
	for(auto histo:yRapidity.at(0))
		delete histo;
	for(auto histo:yRapidity.at(1))
		delete histo;
	for(auto histo:yRapidity.at(2))
		delete histo;
	for(auto histo:xPt.at(0))
		delete histo;
	for(auto histo:xPt.at(1))
		delete histo;
	for(auto histo:xPt.at(2))
		delete histo;
	for(auto histo:yPt.at(0))
		delete histo;
	for(auto histo:yPt.at(1))
		delete histo;
	for(auto histo:yPt.at(2))
		delete histo;
}

void Flow3SE::InitializeHistograms()
{
	cout << "Flow histograms initialization" << endl;
	auto nbins = 20;
	for( int i=0; i<nbins; i++ )
	{
		xRapidity.at(0).push_back( new TProfile( Form( "v1_vs_y_on_X_SE1_centrality_%i", i ), ";rapidity, y; v1_{x}", 14, -0.7, 0.7) );
		xRapidity.at(1).push_back( new TProfile( Form( "v1_vs_y_on_X_SE2_centrality_%i", i ), ";rapidity, y; v1_{x}", 14, -0.7, 0.7) );
		xRapidity.at(2).push_back( new TProfile( Form( "v1_vs_y_on_X_SE3_centrality_%i", i ), ";rapidity, y; v1_{x}", 14, -0.7, 0.7) );
		yRapidity.at(0).push_back( new TProfile( Form( "v1_vs_y_on_Y_SE1_centrality_%i", i ), ";rapidity, y; v1_{y}", 14, -0.7, 0.7) );
		yRapidity.at(1).push_back( new TProfile( Form( "v1_vs_y_on_Y_SE2_centrality_%i", i ), ";rapidity, y; v1_{y}", 14, -0.7, 0.7) );
		yRapidity.at(2).push_back( new TProfile( Form( "v1_vs_y_on_Y_SE3_centrality_%i", i ), ";rapidity, y; v1_{y}", 14, -0.7, 0.7) );
		
		xPt.at(0).push_back( new TProfile( Form( "v1_vs_pt_on_X_SE1_centrality_%i", i ), ";pt, [#frac{GeV}{c}]; v1_{x}", 7, 0., 1.4) );
		xPt.at(1).push_back( new TProfile( Form( "v1_vs_pt_on_X_SE2_centrality_%i", i ), ";pt, [#frac{GeV}{c}]; v1_{x}", 7, 0., 1.4) );
		xPt.at(2).push_back( new TProfile( Form( "v1_vs_pt_on_X_SE3_centrality_%i", i ), ";pt, [#frac{GeV}{c}]; v1_{x}", 7, 0., 1.4) );
		yPt.at(0).push_back( new TProfile( Form( "v1_vs_pt_on_Y_SE1_centrality_%i", i ), ";pt, [#frac{GeV}{c}]; v1_{y}", 7, 0., 1.4) );
		yPt.at(1).push_back( new TProfile( Form( "v1_vs_pt_on_Y_SE2_centrality_%i", i ), ";pt, [#frac{GeV}{c}]; v1_{y}", 7, 0., 1.4) );
		yPt.at(2).push_back( new TProfile( Form( "v1_vs_pt_on_Y_SE3_centrality_%i", i ), ";pt, [#frac{GeV}{c}]; v1_{y}", 7, 0., 1.4) );
	}
}

void Flow3SE::Estimate()
{
	double BETA = sqrt( 1.0 - 0.938*0.938/1.23/1.23 );
	auto ntracks = fEvent->GetNVertexTracks();
	int centrality = fEvent->GetCentralityClass(HADES_constants::kNhitsTOF_RPC_cut);
	// cout << fCentrality->GetCentralityClass() << endl;
	TVector3 b; b.SetXYZ(0,0,-BETA);
	fQvector->Compute();
	for(int i=0; i<ntracks; i++)
	{
		DataTreeTrack* track = fEvent->GetVertexTrack(i);
		if( !fSelector->IsCorrectTrack(i) )
			continue;
		if( track->GetPdgId() != fPid )
			continue;
		TLorentzVector momentum = track->GetMomentum();
		momentum.Boost(b);
		for( auto &vector : fUvector )
			vector.SetMagPhi( 1., momentum.Phi() );
		for( int j=0; j<3; j++ )
		{
			// cout << fQvector->X(j) << endl;
			if( fQvector->X(j) < -990.0 )
				continue;
			auto flowX = fUvector.at(j).X()*fQvector->X(j) / fQvector->Resolution(j).X();
			auto flowY = fUvector.at(j).Y()*fQvector->Y(j) / fQvector->Resolution(j).Y();
			// cout << "flowX=" << flowX << " flowY=" << flowY << endl;
			auto w = momentum.Rapidity() > 0 ? 1 : -1;
			if( flowX == flowX )
			{
				xRapidity.at(j).at(centrality)->Fill( momentum.Rapidity(), flowX );
				xPt.at(j).at(centrality)->Fill( momentum.Pt(), w*flowX );
			}
			if( flowY == flowY )
			{
				yRapidity.at(j).at(centrality)->Fill( momentum.Rapidity(), flowY );
				yPt.at(j).at(centrality)->Fill( momentum.Pt(), w*flowY );
			}
		}
	}
}

void Flow3SE::SavePictures(TString sFileName)
{
	cout << "Saving pictures as png" << endl;
	vector<TCanvas*> canvas;
	vector<TLegend*> legend;
	array< vector<THStack*>, 4> stack; // 1, 2 - V vs y; 3, 4 - V vs Pt
	array< vector<TProfile*>, 3> hRapidityOnX;
	array< vector<TProfile*>, 3> hRapidityOnY;
	array< vector<TProfile*>, 3> hPtOnX;
	array< vector<TProfile*>, 3> hPtOnY;
	gStyle->SetErrorX(0);
	for( int se=0; se < 3; se++)
	{
		for(int i=0; i<4; i++)
		{
			hRapidityOnX.at(se).push_back( new TProfile( Form("v1_vs_y_on_X_cent_%i_SE%i", i, se), ";rapidity, y_{cm};v{1}", 14, -0.7, 0.7 ) );
			hRapidityOnX.at(se).back()->Sumw2();
			hRapidityOnY.at(se).push_back( new TProfile( Form("v1_vs_y_on_Y_cent_%i_SE%i", i, se), ";rapidity, y_{cm};v{1}", 14, -0.7, 0.7 ) );
			hRapidityOnY.at(se).back()->Sumw2();
			hPtOnX.at(se).push_back( new TProfile( Form("v1_vs_pt_on_X_cent_%i_SE%i", i, se), ";pt, [#frac{GeV}{c}];v{1}", 7, 0., 1.4 ) );
			hPtOnX.at(se).back()->Sumw2();
			hPtOnY.at(se).push_back( new TProfile( Form("v1_vs_pt_on_Y_cent_%i_SE%i", i, se), ";pt, [#frac{GeV}{c}];v{1}", 7, 0., 1.4 ) );
			hPtOnY.at(se).back()->Sumw2();		
		}
	}
	for( int se=0; se<3; se++ )
	{
		for( int i=0; i<4; i++ )
		{
			for( int j=i*2; j<(i+1)*2; j++ )
			{
				float wRapX = xRapidity.at(se).at(j)->GetEntries();
				float wRapY = yRapidity.at(se).at(j)->GetEntries();
				float wPtX = xPt.at(se).at(j)->GetEntries();
				float wPtY = yPt.at(se).at(j)->GetEntries();
				hRapidityOnX.at(se).at(i)->Add( xRapidity.at(se).at(j), wRapX );
				hRapidityOnY.at(se).at(i)->Add( yRapidity.at(se).at(j), wRapY );
				hPtOnX.at(se).at(i)->Add( xPt.at(se).at(j), wPtX );
				hPtOnY.at(se).at(i)->Add( yPt.at(se).at(j), wPtY );
			}
		}
	}
	for(int i=0; i<4; i++)
	{
		stack.at(0).push_back( new THStack( Form("v_{1}^{x} vs y_{cm}, cent=%i",i ), Form("%i-%i centrality;y_{cm};v_{1}^{x}",i*10,(i+1)*10) ) );
		stack.at(1).push_back( new THStack( Form("v_{1}^{y} vs y_{cm}, cent=%i",i ), Form("%i-%i centrality;y_{cm};v_{1}^{y}",i*10,(i+1)*10) ) );
		stack.at(2).push_back( new THStack( Form("v_{1}^{x} vs pt, cent=%i",i ), Form("%i-%i centrality;pt, [#frac{GeV}{c}];v_{1}^{x}", i*10, (i+1)*10) ) );
		stack.at(3).push_back( new THStack( Form("v_{1}^{y} vs pt, cent=%i",i ), Form("%i-%i centrality;pt, [#frac{GeV}{c}];v_{1}^{y}", i*10, (i+1)*10) ) );
		for( int se=0; se<3; se++ )
		{	
			stack.at(0).back()->Add( hRapidityOnX.at(se).at(i) );
			stack.at(1).back()->Add( hRapidityOnY.at(se).at(i) );
			stack.at(2).back()->Add( hPtOnX.at(se).at(i) );
			stack.at(3).back()->Add( hPtOnY.at(se).at(i) );
		}
	}
	for(int i=0; i<4; i++)
	{
		for(int se=0; se<3; se++)
		{
			hRapidityOnX.at(se).at(i)->SetLineColor(se+1);
			hRapidityOnX.at(se).at(i)->SetMarkerColor(se+1);
			hRapidityOnX.at(se).at(i)->SetMarkerStyle(20+se);
			hRapidityOnX.at(se).at(i)->SetMarkerSize(2);
			hRapidityOnX.at(se).at(i)->SetLineWidth(2);
			hRapidityOnX.at(se).at(i)->GetYaxis()->SetRangeUser(-0.5,0.5);

			hRapidityOnY.at(se).at(i)->SetLineColor(se+1);
			hRapidityOnY.at(se).at(i)->SetMarkerColor(se+1);
			hRapidityOnY.at(se).at(i)->SetMarkerStyle(20+se);
			hRapidityOnY.at(se).at(i)->SetMarkerSize(2);
			hRapidityOnY.at(se).at(i)->SetLineWidth(2);
			hRapidityOnY.at(se).at(i)->GetYaxis()->SetRangeUser(-0.5,0.5);
			
			hPtOnX.at(se).at(i)->SetLineColor(se+1);
			hPtOnX.at(se).at(i)->SetMarkerColor(se+1);
			hPtOnX.at(se).at(i)->SetMarkerStyle(20+se);
			hPtOnX.at(se).at(i)->SetMarkerSize(2);
			hPtOnX.at(se).at(i)->SetLineWidth(2);
			hPtOnX.at(se).at(i)->GetYaxis()->SetRangeUser(0.0,0.3);
			
			hPtOnY.at(se).at(i)->SetLineColor(se+1);
			hPtOnY.at(se).at(i)->SetMarkerColor(se+1);
			hPtOnY.at(se).at(i)->SetMarkerStyle(20+se);
			hPtOnY.at(se).at(i)->SetMarkerSize(2);
			hPtOnY.at(se).at(i)->SetLineWidth(2);
			hPtOnY.at(se).at(i)->GetYaxis()->SetRangeUser(0.0,0.3);

		}
	}
	legend.push_back( new TLegend(0.1,0.8,0.38,0.9) );
	legend.back()->AddEntry(hRapidityOnX.at(0).at(0), "SE 1");
	legend.back()->AddEntry(hRapidityOnX.at(1).at(0), "SE 2");
	legend.back()->AddEntry(hRapidityOnX.at(2).at(0), "SE 3");
	canvas.push_back( new TCanvas( "canv1", "", 2000, 750 ) );
	canvas.push_back( new TCanvas( "canv2", "", 2000, 750 ) );
	canvas.push_back( new TCanvas( "canv3", "", 2000, 750 ) );
	canvas.push_back( new TCanvas( "canv4", "", 2000, 750 ) );
	for( int i=0; i<4; i++ )
	{
		canvas.at(i)->Divide(4,1,0.01,0.1);
		for(int j=0; j<4; j++)
		{
			canvas.at(i)->cd(j+1);
			stack.at(i).at(j)->Draw();
			legend.back()->Draw();
		}
		canvas.at(i)->SaveAs( sFileName+Form("flow_%i",i)+".png" );
	}
}

void Flow3SE::SaveHistogramsToRootFile(TString sFileName)
{
	cout << "Saving pictures to ROOT file" << endl;
	auto file = new TFile( sFileName+".root", "recreate" );
	array< vector<TProfile*>, 3> hRapidityOnX;
	array< vector<TProfile*>, 3> hRapidityOnY;
	array< vector<TProfile*>, 3> hPtOnX;
	array< vector<TProfile*>, 3> hPtOnY;
	for( int se=0; se < 3; se++)
	{
		for(int i=0; i<4; i++)
		{
			hRapidityOnX.at(se).push_back( new TProfile( Form("v1_vs_y_on_X_cent_%i_SE%i", i, se), ";rapidity, y_{cm};v{1}", 14, -0.7, 0.7 ) );
			hRapidityOnX.at(se).back()->Sumw2();
			hRapidityOnY.at(se).push_back( new TProfile( Form("v1_vs_y_on_Y_cent_%i_SE%i", i, se), ";rapidity, y_{cm};v{1}", 14, -0.7, 0.7 ) );
			hRapidityOnY.at(se).back()->Sumw2();
			hPtOnX.at(se).push_back( new TProfile( Form("v1_vs_pt_on_X_cent_%i_SE%i", i, se), ";pt, [#frac{GeV}{c}];v{1}", 7, 0., 1.4 ) );
			hPtOnX.at(se).back()->Sumw2();
			hPtOnY.at(se).push_back( new TProfile( Form("v1_vs_pt_on_Y_cent_%i_SE%i", i, se), ";pt, [#frac{GeV}{c}];v{1}", 7, 0., 1.4 ) );
			hPtOnY.at(se).back()->Sumw2();		
		}
	}
	for( int se=0; se<3; se++ )
	{
		for( int i=0; i<4; i++ )
		{
			for( int j=i*2; j<(i+1)*2; j++ )
			{
				float wRapX = xRapidity.at(se).at(j)->GetEntries();
				float wRapY = yRapidity.at(se).at(j)->GetEntries();
				float wPtX = xPt.at(se).at(j)->GetEntries();
				float wPtY = yPt.at(se).at(j)->GetEntries();
				hRapidityOnX.at(se).at(i)->Add( xRapidity.at(se).at(j), wRapX );
				hRapidityOnY.at(se).at(i)->Add( yRapidity.at(se).at(j), wRapY );
				hPtOnX.at(se).at(i)->Add( xPt.at(se).at(j), wPtX );
				hPtOnY.at(se).at(i)->Add( yPt.at(se).at(j), wPtY );
			}
		}
	}
	file->cd();
	for( auto &vec : hRapidityOnX )
	{
		for( auto histo : vec )
			histo->Write();
	}
	for( auto &vec : hRapidityOnY )
	{
		for( auto histo : vec )
			histo->Write();
	}
	for( auto &vec : hPtOnX )
	{
		for( auto histo : vec )
			histo->Write();
	}
	for( auto &vec : hPtOnY )
	{
		for( auto histo : vec )
			histo->Write();
	}
	file->Close();
}