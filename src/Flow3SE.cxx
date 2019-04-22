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
	vector<TCanvas*> canvas;
	array< vector<THStack*>, 4> stack;
	array< vector<TProfile*>, 3> hRapidityOnX;
	array< vector<TProfile*>, 3> hRapidityOnY;
	array< vector<TProfile*>, 3> hPtOnX;
	array< vector<TProfile*>, 3> hPtOnY;
	gStyle->SetErrorX(0);
	for( int se=0; se < 3; i++)
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
				hRapidityOnX.at(se).at(i)->Add( xRapidity.at(k).at(j), wRapX );
				hRapidityOnY.at(se).at(i)->Add( yRapidity.at(k).at(j), wRapY );
				hPtOnX.at(se).at(i)->Add( xPt.at(k).at(j), wPtX );
				hPtOnY.at(se).at(i)->Add( yPt.at(k).at(j), wPtY );
			}
		}
	}
	stack.at(i).push_back( new THStack( "" ) )
	for(int i=0; i<4; i++)
	{
		hRapidityOnX.at(i)->SetLineColor(i+1);
		hRapidityOnX.at(i)->SetMarkerColor(i+1);
		hRapidityOnX.at(i)->SetMarkerStyle(20+i);
		hRapidityOnX.at(i)->SetMarkerSize(3);
		hRapidityOnX.at(i)->SetLineWidth(4);

		hRapidityOnY.at(i)->SetLineColor(i+1);
		hRapidityOnY.at(i)->SetMarkerColor(i+1);
		hRapidityOnY.at(i)->SetMarkerStyle(20+i);
		hRapidityOnY.at(i)->SetMarkerSize(3);
		hRapidityOnY.at(i)->SetLineWidth(4);
		
		hPtOnX.at(i)->SetLineColor(i+1);
		hPtOnX.at(i)->SetMarkerColor(i+1);
		hPtOnX.at(i)->SetMarkerStyle(20+i);
		hPtOnX.at(i)->SetMarkerSize(3);
		hPtOnX.at(i)->SetLineWidth(4);
		
		hPtOnY.at(i)->SetLineColor(i+1);
		hPtOnY.at(i)->SetMarkerColor(i+1);
		hPtOnY.at(i)->SetMarkerStyle(20+i);
		hPtOnY.at(i)->SetMarkerSize(3);
		hPtOnY.at(i)->SetLineWidth(4);
		
		stack.at(0)->Add( hRapidityOnX.at(i) );
		stack.at(1)->Add( hRapidityOnY.at(i) );
		stack.at(2)->Add( hPtOnX.at(i) );
		stack.at(3)->Add( hPtOnY.at(i) );
	}
	canvas.push_back( new TCanvas( "canv", "", 4000, 1500 ) );
	// canvas.back()->cd();
	// xRapidity.at(0).at(4)->Draw();
	canvas.back()->Divide(4,1);
	for(int i=0;i<4;i++)
	{
		canvas.back()->cd(i+1);
		stack.at(i)->Draw("NOSTACK");
		gPad->BuildLegend(0.1,0.75,0.38,0.9);
	}
	canvas.back()->SaveAs( "../histograms/"+sFileName+".png" );
}
