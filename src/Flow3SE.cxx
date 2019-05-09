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
	}
	for(int i=0; i<3; i++)
	{
		xPt.at(i).push_back( new TProfile( Form( "-0.05 < y_{cm} < 0.05, SE%i_{x}",  i ), ";pt, [#frac{GeV}{c}]; v1_{x}", 25, 0.0, 2.5) );
		xPt.at(i).push_back( new TProfile( Form( "-0.25 < y_{cm} < -0.15, SE%i_{x}", i ), ";pt, [#frac{GeV}{c}]; v1_{x}", 25, 0.0, 2.5) );
		xPt.at(i).push_back( new TProfile( Form( "-0.45 < y_{cm} < -0.35, SE%i_{x}", i ), ";pt, [#frac{GeV}{c}]; v1_{x}", 25, 0.0, 2.5) );
		xPt.at(i).push_back( new TProfile( Form( "-0.65 < y_{cm} < -0.55, SE%i_{x}", i ), ";pt, [#frac{GeV}{c}]; v1_{x}", 25, 0.0, 2.5) );
		xPt.at(i).push_back( new TProfile( Form( "-0.75 < y_{cm} < -0.65, SE%i_{x}", i ), ";pt, [#frac{GeV}{c}]; v1_{x}", 25, 0.0, 2.5) );
		
		yPt.at(i).push_back( new TProfile( Form( "-0.05 < y_{cm} < 0.05, SE%i_{y}",  i ), ";pt, [#frac{GeV}{c}]; v1_{x}", 25, 0.0, 2.5) );
		yPt.at(i).push_back( new TProfile( Form( "-0.25 < y_{cm} < -0.15, SE%i_{y}", i ), ";pt, [#frac{GeV}{c}]; v1_{x}", 25, 0.0, 2.5) );
		yPt.at(i).push_back( new TProfile( Form( "-0.45 < y_{cm} < -0.35, SE%i_{y}", i ), ";pt, [#frac{GeV}{c}]; v1_{x}", 25, 0.0, 2.5) );
		yPt.at(i).push_back( new TProfile( Form( "-0.65 < y_{cm} < -0.55, SE%i_{y}", i ), ";pt, [#frac{GeV}{c}]; v1_{x}", 25, 0.0, 2.5) );
		yPt.at(i).push_back( new TProfile( Form( "-0.75 < y_{cm} < -0.65, SE%i_{y}", i ), ";pt, [#frac{GeV}{c}]; v1_{x}", 25, 0.0, 2.5) );
	}
	for( auto &vector : xPt )
		for( auto histo : vector )
			histo->Sumw2();
	for( auto &vector : yPt )
		for( auto histo : vector )
			histo->Sumw2();
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
			if( flowX != flowX )
				continue;
			if( flowY != flowY )
				continue;
			xRapidity.at(j).at(centrality)->Fill( momentum.Rapidity(), flowX );
			yRapidity.at(j).at(centrality)->Fill( momentum.Rapidity(), flowY );
			if( fEvent->GetCentrality(HADES_constants::kNhitsTOF_RPC_cut) > 20 && fEvent->GetCentrality(HADES_constants::kNhitsTOF_RPC_cut) < 30 )
				continue;
			if( momentum.Rapidity() > -0.05 && momentum.Rapidity() < 0.05 )
			{
				xPt.at(j).at(0)->Fill( momentum.Pt(), flowX );
				yPt.at(j).at(0)->Fill( momentum.Pt(), flowX );
			}
			if( momentum.Rapidity() > -0.25 && momentum.Rapidity() < -0.15 )
			{
				xPt.at(j).at(1)->Fill( momentum.Pt(), flowX );
				yPt.at(j).at(1)->Fill( momentum.Pt(), flowX );
			}
			if( momentum.Rapidity() > -0.45 && momentum.Rapidity() < -0.35 )
			{
				xPt.at(j).at(2)->Fill( momentum.Pt(), flowX );
				yPt.at(j).at(2)->Fill( momentum.Pt(), flowX );
			}
			if( momentum.Rapidity() > -0.65 && momentum.Rapidity() < -0.55 )
			{
				xPt.at(j).at(3)->Fill( momentum.Pt(), flowX );
				yPt.at(j).at(3)->Fill( momentum.Pt(), flowX );
			}
			if( momentum.Rapidity() > -0.75 && momentum.Rapidity() < -0.65 )
			{
				xPt.at(j).at(4)->Fill( momentum.Pt(), flowX );
				yPt.at(j).at(4)->Fill( momentum.Pt(), flowX );
			}
		}
	}
}

void Flow3SE::SavePictures(TString sFileName)
{
	cout << "Saving pictures as png" << endl;
	vector<TCanvas*> canvas;
	vector<TLegend*> legend;
	array< vector<THStack*>, 5> stack; // 1, 2 - V vs y; 3, 4 - V vs Pt
	array< vector<TProfile*>, 3> hRapidityOnX;
	array< vector<TProfile*>, 3> hRapidityOnY;
	vector<TProfile*> hPt;
	gStyle->SetErrorX(0);
	for( int se=0; se < 3; se++)
	{
		for(int i=0; i<4; i++)
		{
			hRapidityOnX.at(se).push_back( new TProfile( Form("v1_vs_y_on_X_cent_%i_SE%i", i, se), ";rapidity, y_{cm};v{1}", 14, -0.7, 0.7 ) );
			hRapidityOnX.at(se).back()->Sumw2();
			hRapidityOnY.at(se).push_back( new TProfile( Form("v1_vs_y_on_Y_cent_%i_SE%i", i, se), ";rapidity, y_{cm};v{1}", 14, -0.7, 0.7 ) );
			hRapidityOnY.at(se).back()->Sumw2();
		}
	}
	hPt.push_back( new TProfile( "-0.05 < y_{cm} < 0.05", ";pt, [#frac{GeV}{c}]; v_{1}", 25, 0.0, 2.5) );	// 0
	hPt.push_back( new TProfile( "-0.25 < y_{cm} < -0.15", ";pt, [#frac{GeV}{c}]; v_{1}", 25, 0.0, 2.5) );	// 1
	hPt.push_back( new TProfile( "-0.45 < y_{cm} < -0.35", ";pt, [#frac{GeV}{c}]; v_{1}", 25, 0.0, 2.5) );	// 2
	hPt.push_back( new TProfile( "-0.65 < y_{cm} < -0.55", ";pt, [#frac{GeV}{c}]; v_{1}", 25, 0.0, 2.5) );	// 3
	hPt.push_back( new TProfile( "-0.75 < y_{cm} < -0.65", ";pt, [#frac{GeV}{c}]; v_{1}", 25, 0.0, 2.5) );	// 4
	
	for( int i=0; i< hPt.size(); i++ )
	{
		hPt.at(i)->Sumw2();
		for(int se=0; se<xPt.size(); se++)
		{
			float w1 = xPt.at(se).at(i)->GetEntries();
			hPt.at(i)->Add( xPt.at(se).at(i), w1 );
			float w2 = yPt.at(se).at(i)->GetEntries();
			hPt.at(i)->Add( yPt.at(se).at(i), w2 );
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
				hRapidityOnX.at(se).at(i)->Add( xRapidity.at(se).at(j), wRapX );
				hRapidityOnY.at(se).at(i)->Add( yRapidity.at(se).at(j), wRapY );
			}
		}
	}
	for(int i=0; i<4; i++)
	{
		stack.at(0).push_back( new THStack( Form("v_{1}^{x} vs y_{cm}, cent=%i",i ), Form("%i-%i centrality;y_{cm};v_{1}^{x}",i*10,(i+1)*10) ) );
		stack.at(1).push_back( new THStack( Form("v_{1}^{y} vs y_{cm}, cent=%i",i ), Form("%i-%i centrality;y_{cm};v_{1}^{y}",i*10,(i+1)*10) ) );
		for( int se=0; se<3; se++ )
		{	
			stack.at(0).back()->Add( hRapidityOnX.at(se).at(i) );
			stack.at(1).back()->Add( hRapidityOnY.at(se).at(i) );
		}
	}
	for( int se=0; se<3; se++ )
	{
		stack.at(2).push_back( new THStack( Form("v_{1}^{x} vs pt, SE%i",se ), Form("SE%i;pt, #frac{GeV}{c};v_{1}^{x}",se) ) );
		stack.at(3).push_back( new THStack( Form("v_{1}^{y} vs pt, SE%i",se ), Form("SE%i;pt, #frac{GeV}{c};v_{1}^{y}",se) ) );
		for( int i=0; i<xPt.at(se).size(); i++ )
		{
			stack.at(2).back()->Add( xPt.at(se).at(i) );
			stack.at(3).back()->Add( yPt.at(se).at(i) );
		}
	}
	stack.at(4).push_back( new THStack( "v_{1}^{y} vs pt", ";pt, #frac{GeV}{c};v_{1}") );
	for( int i=0; i<hPt.size(); i++ )
		stack.at(4).back()->Add( hPt.at(i) );
	for( int i=0; i<5; i++ )
	{
		for( int se=0; se<3; se++ )
		{
			xPt.at(se).at(i)->SetLineColor(i+1);
			xPt.at(se).at(i)->SetMarkerColor(i+1);
			xPt.at(se).at(i)->SetMarkerStyle(20+i);
			xPt.at(se).at(i)->SetMarkerSize(1);
			xPt.at(se).at(i)->SetLineWidth(2);
			xPt.at(se).at(i)->GetYaxis()->SetRangeUser(-0.6,0.4);

			yPt.at(se).at(i)->SetLineColor(i+1);
			yPt.at(se).at(i)->SetMarkerColor(i+1);
			yPt.at(se).at(i)->SetMarkerStyle(20+i);
			yPt.at(se).at(i)->SetMarkerSize(1);
			yPt.at(se).at(i)->SetLineWidth(2);
			yPt.at(se).at(i)->GetYaxis()->SetRangeUser(-0.6,0.4);
		}
		hPt.at(i)->SetLineColor(i+1);
		hPt.at(i)->SetMarkerColor(i+1);
		hPt.at(i)->SetMarkerStyle(20+i);
		hPt.at(i)->SetMarkerSize(1);
		hPt.at(i)->SetLineWidth(2);
		hPt.at(i)->GetYaxis()->SetRangeUser(-0.6,0.4);
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
		}
	}
	legend.push_back( new TLegend(0.1,0.8,0.38,0.9) );
	legend.back()->AddEntry(hRapidityOnX.at(0).at(0), "SE 1");
	legend.back()->AddEntry(hRapidityOnX.at(1).at(0), "SE 2");
	legend.back()->AddEntry(hRapidityOnX.at(2).at(0), "SE 3");
	canvas.push_back( new TCanvas( "canv1", "", 2000, 750 ) );
	canvas.push_back( new TCanvas( "canv2", "", 2000, 750 ) );
	canvas.push_back( new TCanvas( "canv3", "", 1500, 450 ) );
	canvas.push_back( new TCanvas( "canv4", "", 1500, 450 ) );
	for( int i=0; i<2; i++ )
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
	for( int i=2; i<4; i++ )
	{
		canvas.at(i)->Divide(3,1);
		for(int j=0; j<3; j++)
		{
			canvas.at(i)->cd(j+1);
			stack.at(i).at(j)->Draw();
			gPad->BuildLegend(0.1,0.65,0.5,0.9);
		}
		canvas.at(i)->SaveAs( sFileName+Form("flow_%i",i)+".png" );
	}
	canvas.push_back( new TCanvas( "canv5", "", 900, 750 ) );
	canvas.back()->cd();
	stack.back().back()->Draw("NOSTACK");
	gPad->BuildLegend(0.1,0.75,0.5,0.9);
	canvas.back()->SaveAs( sFileName+"flow_4.png" );
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