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
	for(int i=0; i<3; i++)
		fFlow.push_back( TVector2{0.,0.} );
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
	for( int se=0; se<3; se++ )
	{
		xRapidity.at(se).push_back( new TProfile( Form( "0.25 < pt < 0.3 SE%i_{x}", se ), ";y_{cm}; v_{1}^{x}", 14, -0.7, 0.7) ); // 0
		xRapidity.at(se).push_back( new TProfile( Form( "0.4 < pt < 0.45 SE%i_{x}", se ), ";y_{cm}; v_{1}^{x}", 14, -0.7, 0.7) ); // 1
		xRapidity.at(se).push_back( new TProfile( Form( "0.6 < pt < 0.65 SE%i_{x}", se ), ";y_{cm}; v_{1}^{x}", 14, -0.7, 0.7) ); // 2
		xRapidity.at(se).push_back( new TProfile( Form( "0.8 < pt < 0.85 SE%i_{x}", se ), ";y_{cm}; v_{1}^{x}", 14, -0.7, 0.7) ); // 3
		xRapidity.at(se).push_back( new TProfile( Form( "1.0 < pt < 1.05 SE%i_{x}", se ), ";y_{cm}; v_{1}^{x}", 14, -0.7, 0.7) ); // 4
		xRapidity.at(se).push_back( new TProfile( Form( "1.2 < pt < 1.25 SE%i_{x}", se ), ";y_{cm}; v_{1}^{x}", 14, -0.7, 0.7) ); // 5

		yRapidity.at(se).push_back( new TProfile( Form( "0.25 < pt < 0.3 SE%i_{y}", se ), ";y_{cm}; v_{1}^{y}", 14, -0.7, 0.7) );
		yRapidity.at(se).push_back( new TProfile( Form( "0.4 < pt < 0.45 SE%i_{y}", se ), ";y_{cm}; v_{1}^{y}", 14, -0.7, 0.7) );
		yRapidity.at(se).push_back( new TProfile( Form( "0.6 < pt < 0.65 SE%i_{y}", se ), ";y_{cm}; v_{1}^{y}", 14, -0.7, 0.7) );
		yRapidity.at(se).push_back( new TProfile( Form( "0.8 < pt < 0.85 SE%i_{y}", se ), ";y_{cm}; v_{1}^{y}", 14, -0.7, 0.7) );
		yRapidity.at(se).push_back( new TProfile( Form( "1.0 < pt < 1.05 SE%i_{y}", se ), ";y_{cm}; v_{1}^{y}", 14, -0.7, 0.7) );
		yRapidity.at(se).push_back( new TProfile( Form( "1.2 < pt < 1.25 SE%i_{y}", se ), ";y_{cm}; v_{1}^{y}", 14, -0.7, 0.7) );
	}
	for(int i=0; i<3; i++)
	{
		xPt.at(i).push_back( new TProfile( Form( "-0.05 < y_{cm} < 0.05, SE%i_{x}",  i ), ";pt, [#frac{GeV}{c}]; v1_{x}", 25, 0.0, 2.5) ); // 0
		xPt.at(i).push_back( new TProfile( Form( "-0.25 < y_{cm} < -0.15, SE%i_{x}", i ), ";pt, [#frac{GeV}{c}]; v1_{x}", 25, 0.0, 2.5) ); // 1
		xPt.at(i).push_back( new TProfile( Form( "-0.45 < y_{cm} < -0.35, SE%i_{x}", i ), ";pt, [#frac{GeV}{c}]; v1_{x}", 25, 0.0, 2.5) ); // 2
		xPt.at(i).push_back( new TProfile( Form( "-0.65 < y_{cm} < -0.55, SE%i_{x}", i ), ";pt, [#frac{GeV}{c}]; v1_{x}", 25, 0.0, 2.5) ); // 3
		xPt.at(i).push_back( new TProfile( Form( "-0.75 < y_{cm} < -0.65, SE%i_{x}", i ), ";pt, [#frac{GeV}{c}]; v1_{x}", 25, 0.0, 2.5) ); // 4
		
		yPt.at(i).push_back( new TProfile( Form( "-0.05 < y_{cm} < 0.05, SE%i_{y}",  i ), ";pt, [#frac{GeV}{c}]; v1_{x}", 25, 0.0, 2.5) );
		yPt.at(i).push_back( new TProfile( Form( "-0.25 < y_{cm} < -0.15, SE%i_{y}", i ), ";pt, [#frac{GeV}{c}]; v1_{x}", 25, 0.0, 2.5) );
		yPt.at(i).push_back( new TProfile( Form( "-0.45 < y_{cm} < -0.35, SE%i_{y}", i ), ";pt, [#frac{GeV}{c}]; v1_{x}", 25, 0.0, 2.5) );
		yPt.at(i).push_back( new TProfile( Form( "-0.65 < y_{cm} < -0.55, SE%i_{y}", i ), ";pt, [#frac{GeV}{c}]; v1_{x}", 25, 0.0, 2.5) );
		yPt.at(i).push_back( new TProfile( Form( "-0.75 < y_{cm} < -0.65, SE%i_{y}", i ), ";pt, [#frac{GeV}{c}]; v1_{x}", 25, 0.0, 2.5) );
	}
	for( int se=0; se<3; se++ )
	{
		hResolution.at(se).push_back( new TProfile( Form("resolution, SE%i_{x}", se), "resolution, SP;centrality;R_{1}^{x}", 20, 0, 100 ) );
		hResolution.at(se).push_back( new TProfile( Form("resolution, SE%i_{y}", se), "resolution, SP;centrality;R_{1}^{y}", 20, 0, 100 ) );
	}
	this->InitializeQvectorCorrelations();
	hKinematics.push_back( new TH2F( "m^{2} vs p", ";p, #frac{GeV}{c};m^{2}, #frac{GeV^{2}}{c^{4}}", 100, -3.0, 3.0, 100, 0.0, 18.0 ) );
	hKinematics.push_back( new TH2F( "pt vs y_{cm}", ";y_{cm};pt, #frac{GeV}{c}", 100, -1.0, 1.0, 100, 0.0, 3.0 ) );
	for( auto &vector : xPt )
		for( auto histo : vector )
			histo->Sumw2();
	for( auto &vector : yPt )
		for( auto histo : vector )
			histo->Sumw2();
}

void Flow3SE::InitializeQvectorCorrelations()
{
	hCorrelation.push_back( new TProfile("Qx_{a}Qx_{b}", ";Centrality;Qx_{a}Qx_{b}", 10, 0, 50) ); // 0
	hCorrelation.push_back( new TProfile("Qy_{a}Qy_{b}", ";Centrality;Qy_{a}Qy_{b}", 10, 0, 50) ); // 1
	hCorrelation.push_back( new TProfile("Qx_{a}Qy_{b}", ";Centrality;Qx_{a}Qy_{b}", 10, 0, 50) ); // 2
	hCorrelation.push_back( new TProfile("Qy_{a}Qx_{b}", ";Centrality;Qy_{a}Qx_{b}", 10, 0, 50) ); // 3

	hCorrelation.push_back( new TProfile("Qx_{b}Qx_{c}", ";Centrality;Qx_{b}Qx_{c}", 10, 0, 50) ); // 4
	hCorrelation.push_back( new TProfile("Qy_{b}Qy_{c}", ";Centrality;Qy_{b}Qy_{c}", 10, 0, 50) ); // 5
	hCorrelation.push_back( new TProfile("Qx_{b}Qy_{c}", ";Centrality;Qx_{b}Qx_{c}", 10, 0, 50) ); // 6
	hCorrelation.push_back( new TProfile("Qy_{b}Qx_{c}", ";Centrality;Qy_{b}Qx_{c}", 10, 0, 50) ); // 7

	hCorrelation.push_back( new TProfile("Qx_{a}Qx_{c}", ";Centrality;Qx_{a}Qx_{c}", 10, 0, 50) ); // 8
	hCorrelation.push_back( new TProfile("Qy_{a}Qy_{c}", ";Centrality;Qy_{a}Qy_{c}", 10, 0, 50) ); // 9
	hCorrelation.push_back( new TProfile("Qx_{a}Qy_{c}", ";Centrality;Qx_{a}Qx_{c}", 10, 0, 50) ); // 10
	hCorrelation.push_back( new TProfile("Qy_{a}Qx_{c}", ";Centrality;Qy_{a}Qx_{c}", 10, 0, 50) ); // 11
	this->InitializeQaHistograms();
}

void Flow3SE::InitializeObservableFlow()
{
	for( int se=0; se<3; se++ )
	{
		xObsRapidity.at(se).push_back( new TProfile( Form( "0.25 < pt < 0.3 SE%i_{x}", se ), ";y_{cm}; v_{1}^{x}", 14, -0.7, 0.7) ); // 0
		xObsRapidity.at(se).push_back( new TProfile( Form( "0.4 < pt < 0.45 SE%i_{x}", se ), ";y_{cm}; v_{1}^{x}", 14, -0.7, 0.7) ); // 1
		xObsRapidity.at(se).push_back( new TProfile( Form( "0.6 < pt < 0.65 SE%i_{x}", se ), ";y_{cm}; v_{1}^{x}", 14, -0.7, 0.7) ); // 2
		xObsRapidity.at(se).push_back( new TProfile( Form( "0.8 < pt < 0.85 SE%i_{x}", se ), ";y_{cm}; v_{1}^{x}", 14, -0.7, 0.7) ); // 3
		xObsRapidity.at(se).push_back( new TProfile( Form( "1.0 < pt < 1.05 SE%i_{x}", se ), ";y_{cm}; v_{1}^{x}", 14, -0.7, 0.7) ); // 4
		xObsRapidity.at(se).push_back( new TProfile( Form( "1.2 < pt < 1.25 SE%i_{x}", se ), ";y_{cm}; v_{1}^{x}", 14, -0.7, 0.7) ); // 5

		yObsRapidity.at(se).push_back( new TProfile( Form( "0.25 < pt < 0.3 SE%i_{y}", se ), ";y_{cm}; v_{1}^{y}", 14, -0.7, 0.7) );
		yObsRapidity.at(se).push_back( new TProfile( Form( "0.4 < pt < 0.45 SE%i_{y}", se ), ";y_{cm}; v_{1}^{y}", 14, -0.7, 0.7) );
		yObsRapidity.at(se).push_back( new TProfile( Form( "0.6 < pt < 0.65 SE%i_{y}", se ), ";y_{cm}; v_{1}^{y}", 14, -0.7, 0.7) );
		yObsRapidity.at(se).push_back( new TProfile( Form( "0.8 < pt < 0.85 SE%i_{y}", se ), ";y_{cm}; v_{1}^{y}", 14, -0.7, 0.7) );
		yObsRapidity.at(se).push_back( new TProfile( Form( "1.0 < pt < 1.05 SE%i_{y}", se ), ";y_{cm}; v_{1}^{y}", 14, -0.7, 0.7) );
		yObsRapidity.at(se).push_back( new TProfile( Form( "1.2 < pt < 1.25 SE%i_{y}", se ), ";y_{cm}; v_{1}^{y}", 14, -0.7, 0.7) );
	}
	for(int i=0; i<3; i++)
	{
		xObsPt.at(i).push_back( new TProfile( Form( "-0.05 < y_{cm} < 0.05, SE%i_{x}",  i ), ";pt, [#frac{GeV}{c}]; v1_{x}", 25, 0.0, 2.5) ); // 0
		xObsPt.at(i).push_back( new TProfile( Form( "-0.25 < y_{cm} < -0.15, SE%i_{x}", i ), ";pt, [#frac{GeV}{c}]; v1_{x}", 25, 0.0, 2.5) ); // 1
		xObsPt.at(i).push_back( new TProfile( Form( "-0.45 < y_{cm} < -0.35, SE%i_{x}", i ), ";pt, [#frac{GeV}{c}]; v1_{x}", 25, 0.0, 2.5) ); // 2
		xObsPt.at(i).push_back( new TProfile( Form( "-0.65 < y_{cm} < -0.55, SE%i_{x}", i ), ";pt, [#frac{GeV}{c}]; v1_{x}", 25, 0.0, 2.5) ); // 3
		xObsPt.at(i).push_back( new TProfile( Form( "-0.75 < y_{cm} < -0.65, SE%i_{x}", i ), ";pt, [#frac{GeV}{c}]; v1_{x}", 25, 0.0, 2.5) ); // 4
		
		yObsPt.at(i).push_back( new TProfile( Form( "-0.05 < y_{cm} < 0.05, SE%i_{y}",  i ), ";pt, [#frac{GeV}{c}]; v1_{x}", 25, 0.0, 2.5) );
		yObsPt.at(i).push_back( new TProfile( Form( "-0.25 < y_{cm} < -0.15, SE%i_{y}", i ), ";pt, [#frac{GeV}{c}]; v1_{x}", 25, 0.0, 2.5) );
		yObsPt.at(i).push_back( new TProfile( Form( "-0.45 < y_{cm} < -0.35, SE%i_{y}", i ), ";pt, [#frac{GeV}{c}]; v1_{x}", 25, 0.0, 2.5) );
		yObsPt.at(i).push_back( new TProfile( Form( "-0.65 < y_{cm} < -0.55, SE%i_{y}", i ), ";pt, [#frac{GeV}{c}]; v1_{x}", 25, 0.0, 2.5) );
		yObsPt.at(i).push_back( new TProfile( Form( "-0.75 < y_{cm} < -0.65, SE%i_{y}", i ), ";pt, [#frac{GeV}{c}]; v1_{x}", 25, 0.0, 2.5) );
	}
}

void Flow3SE::Estimate()
{
	this->FillQaHistograms();
	double BETA = sqrt( 1.0 - 0.938*0.938/1.23/1.23 );
	auto ntracks = fEvent->GetNVertexTracks();
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
		for( auto &vector : fFlow )
			vector.Set( 0., 0. );
		for( int se=0; se<3; se++ )
		{
			if( fQvector->X(se) < -990.0 )
				continue;
			fFlow.at(se).Set( fUvector.at(se).X()*fQvector->X(se) / fQvector->Resolution(se).X(), fUvector.at(se).Y()*fQvector->Y(se) / fQvector->Resolution(se).Y() );
		}
		this->FillPtDependence(i);
		this->FillYDependence(i);
	}
}

void Flow3SE::FillYDependence(int trackIdx)
{
	double BETA = sqrt( 1.0 - 0.938*0.938/1.23/1.23 );
	TVector3 b; b.SetXYZ(0,0,-BETA);
	DataTreeTrack* track = fEvent->GetVertexTrack(trackIdx);
	TLorentzVector momentum = track->GetMomentum();
	momentum.Boost(b);
	for( int se=0; se<3; se++ )
	{
		if( fFlow.at(se).X() != fFlow.at(se).X() )
			continue;
		if( fFlow.at(se).Y() != fFlow.at(se).Y() )
			continue;
		if( momentum.Pt() > 0.25 && momentum.Pt() < 0.3 )
		{
			xRapidity.at(se).at(0)->Fill( momentum.Rapidity(), fFlow.at(se).X() );
			yRapidity.at(se).at(0)->Fill( momentum.Rapidity(), fFlow.at(se).Y() );
		}
		if( momentum.Pt() > 0.4 && momentum.Pt() < 0.45 )
		{
			xRapidity.at(se).at(1)->Fill( momentum.Rapidity(), fFlow.at(se).X() );
			yRapidity.at(se).at(1)->Fill( momentum.Rapidity(), fFlow.at(se).Y() );
		}
		if( momentum.Pt() > 0.6 && momentum.Pt() < 0.65 )
		{
			xRapidity.at(se).at(2)->Fill( momentum.Rapidity(), fFlow.at(se).X() );
			yRapidity.at(se).at(2)->Fill( momentum.Rapidity(), fFlow.at(se).Y() );
		}
		if( momentum.Pt() > 0.8 && momentum.Pt() < 0.85 )
		{
			xRapidity.at(se).at(3)->Fill( momentum.Rapidity(), fFlow.at(se).X() );
			yRapidity.at(se).at(3)->Fill( momentum.Rapidity(), fFlow.at(se).Y() );
		}
		if( momentum.Pt() > 1.0 && momentum.Pt() < 1.05 )
		{
			xRapidity.at(se).at(4)->Fill( momentum.Rapidity(), fFlow.at(se).X() );
			yRapidity.at(se).at(4)->Fill( momentum.Rapidity(), fFlow.at(se).Y() );
		}
		if( momentum.Pt() > 1.2 && momentum.Pt() < 1.25 )
		{
			xRapidity.at(se).at(5)->Fill( momentum.Rapidity(), fFlow.at(se).X() );
			yRapidity.at(se).at(5)->Fill( momentum.Rapidity(), fFlow.at(se).Y() );
		}
	}

}

void Flow3SE::FillPtDependence(int trackIdx)
{
	double BETA = sqrt( 1.0 - 0.938*0.938/1.23/1.23 );
	TVector3 b; b.SetXYZ(0,0,-BETA);
	DataTreeTrack* track = fEvent->GetVertexTrack(trackIdx);
	TLorentzVector momentum = track->GetMomentum();
	momentum.Boost(b);
	for( int se=0; se<3; se++ )
	{
		if( fFlow.at(se).X() != fFlow.at(se).X() )
			continue;
		if( fFlow.at(se).Y() != fFlow.at(se).Y() )
			continue;
		if( momentum.Rapidity() > -0.05 && momentum.Rapidity() < 0.05 )
			{
				xPt.at(se).at(0)->Fill( momentum.Pt(), fFlow.at(se).X() );
				yPt.at(se).at(0)->Fill( momentum.Pt(), fFlow.at(se).Y() );
			}
			if( momentum.Rapidity() > -0.25 && momentum.Rapidity() < -0.15 )
			{
				xPt.at(se).at(1)->Fill( momentum.Pt(), fFlow.at(se).X() );
				yPt.at(se).at(1)->Fill( momentum.Pt(), fFlow.at(se).Y() );
			}
			if( momentum.Rapidity() > -0.45 && momentum.Rapidity() < -0.35 )
			{
				xPt.at(se).at(2)->Fill( momentum.Pt(), fFlow.at(se).X() );
				yPt.at(se).at(2)->Fill( momentum.Pt(), fFlow.at(se).Y() );
			}
			if( momentum.Rapidity() > -0.65 && momentum.Rapidity() < -0.55 )
			{
				xPt.at(se).at(3)->Fill( momentum.Pt(), fFlow.at(se).X() );
				yPt.at(se).at(3)->Fill( momentum.Pt(), fFlow.at(se).Y() );
			}
			if( momentum.Rapidity() > -0.75 && momentum.Rapidity() < -0.65 )
			{
				xPt.at(se).at(4)->Fill( momentum.Pt(), fFlow.at(se).X() );
				yPt.at(se).at(4)->Fill( momentum.Pt(), fFlow.at(se).Y() );
			}
	}
}

void Flow3SE::SavePictures(TString sFileName)
{
	cout << "Saving pictures as png" << endl;
	this->SavePtDependence(sFileName);
	this->SaveYDependence(sFileName);
}

void Flow3SE::SavePtDependence(TString sFileName)
{
	cout << "Saving pt-dependence" << endl;
	vector<TProfile*> hPt;
	vector<TCanvas*> canvas;
	array< vector<THStack*>, 3> stack;
	hPt.push_back( new TProfile( "-0.05 < y_{cm} < 0.05", ";pt, [#frac{GeV}{c}]; v_{1}", 25, 0.0, 2.5) );	// 0
	hPt.push_back( new TProfile( "-0.25 < y_{cm} < -0.15", ";pt, [#frac{GeV}{c}]; v_{1}", 25, 0.0, 2.5) );	// 1
	hPt.push_back( new TProfile( "-0.45 < y_{cm} < -0.35", ";pt, [#frac{GeV}{c}]; v_{1}", 25, 0.0, 2.5) );	// 2
	hPt.push_back( new TProfile( "-0.65 < y_{cm} < -0.55", ";pt, [#frac{GeV}{c}]; v_{1}", 25, 0.0, 2.5) );	// 3
	hPt.push_back( new TProfile( "-0.75 < y_{cm} < -0.65", ";pt, [#frac{GeV}{c}]; v_{1}", 25, 0.0, 2.5) );	// 4
	for( int i=0; i<hPt.size(); i++ )
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
		stack.at(0).push_back( new THStack( Form("v_{1}^{x} vs pt, SE%i",se ), Form("SE%i;pt, #frac{GeV}{c};v_{1}^{x}",se) ) );
		stack.at(1).push_back( new THStack( Form("v_{1}^{y} vs pt, SE%i",se ), Form("SE%i;pt, #frac{GeV}{c};v_{1}^{y}",se) ) );
		for( int i=0; i<xPt.at(se).size(); i++ )
		{
			stack.at(0).back()->Add( xPt.at(se).at(i) );
			stack.at(1).back()->Add( yPt.at(se).at(i) );
		}
	}
	stack.at(2).push_back( new THStack( "v_{1}^{y} vs pt", ";pt, #frac{GeV}{c};v_{1}") );
	for( int i=0; i<hPt.size(); i++ )
		stack.at(2).back()->Add( hPt.at(i));
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
	canvas.push_back( new TCanvas( "canv1", "", 1500, 400 ) );
	canvas.push_back( new TCanvas( "canv2", "", 1500, 400 ) );
	for( int i=0; i<2; i++ )
	{
		canvas.at(i)->Divide(3,1);
		for(int j=0; j<3; j++)
		{
			canvas.at(i)->cd(j+1);
			stack.at(i).at(j)->Draw();
			gPad->BuildLegend(0.1,0.65,0.5,0.9);
		}
		canvas.at(i)->SaveAs( sFileName+Form("_pt_%i",i)+".png" );
	}
	canvas.push_back( new TCanvas( "canv5", "", 900, 750 ) );
	canvas.back()->cd();
	stack.back().back()->Draw("NOSTACK");
	gPad->BuildLegend(0.1,0.75,0.5,0.9);
	canvas.back()->SaveAs( sFileName+"_pt_2.png" );
}

void Flow3SE::SaveYDependence(TString sFileName)
{
	cout << "Saving y-dependence" << endl;
	vector<TProfile*> hRapidity;
	vector<TCanvas*> canvas;
	array< vector<THStack*>, 3> stack;
	hRapidity.push_back( new TProfile( "0.25 < pt < 0.3", ";y_{cm};v_{1}", 14, -0.7, 0.7) );	// 0
	hRapidity.push_back( new TProfile( "0.4 < pt < 0.45", ";y_{cm};v_{1}", 14, -0.7, 0.7) );	// 1
	hRapidity.push_back( new TProfile( "0.6 < pt < 0.65", ";y_{cm};v_{1}", 14, -0.7, 0.7) );	// 2
	hRapidity.push_back( new TProfile( "0.8 < pt < 0.85", ";y_{cm};v_{1}", 14, -0.7, 0.7) );	// 3
	hRapidity.push_back( new TProfile( "1.0 < pt < 1.05", ";y_{cm};v_{1}", 14, -0.7, 0.7) );	// 4
	hRapidity.push_back( new TProfile( "1.2 < pt < 1.25", ";y_{cm};v_{1}", 14, -0.7, 0.7) );	// 5
	for( int i=0; i<hRapidity.size(); i++ )
	{
		hRapidity.at(i)->Sumw2();
		for(int se=0; se<xRapidity.size(); se++)
		{
			float w1 = xRapidity.at(se).at(i)->GetEntries();
			hRapidity.at(i)->Add( xRapidity.at(se).at(i), w1 );
			float w2 = yRapidity.at(se).at(i)->GetEntries();
			hRapidity.at(i)->Add( yRapidity.at(se).at(i), w2 );
		}
	}
	for( int se=0; se<3; se++ )
	{
		stack.at(0).push_back( new THStack( Form("v_{1}^{x} vs y, SE%i",se ), Form("SE%i;y_{cm};v_{1}^{x}",se) ) );
		stack.at(1).push_back( new THStack( Form("v_{1}^{y} vs y, SE%i",se ), Form("SE%i;y_{cm};v_{1}^{y}",se) ) );
		for( int i=0; i<xPt.at(se).size(); i++ )
		{
			stack.at(0).back()->Add( xRapidity.at(se).at(i) );
			stack.at(1).back()->Add( yRapidity.at(se).at(i) );
		}
	}
	stack.at(2).push_back( new THStack( "v_{1}^{y} vs y", ";y_{cm};v_{1}") );
	for( int i=0; i<hRapidity.size(); i++ )
		stack.at(2).back()->Add( hRapidity.at(i));
	for( int i=0; i<xRapidity.at(0).size(); i++ )
	{
		for( int se=0; se<3; se++ )
		{
			xRapidity.at(se).at(i)->SetLineColor(i+1);
			xRapidity.at(se).at(i)->SetMarkerColor(i+1);
			xRapidity.at(se).at(i)->SetMarkerStyle(20+i);
			xRapidity.at(se).at(i)->SetMarkerSize(1);
			xRapidity.at(se).at(i)->SetLineWidth(2);
			xRapidity.at(se).at(i)->GetYaxis()->SetRangeUser(-0.6,0.6);

			yRapidity.at(se).at(i)->SetLineColor(i+1);
			yRapidity.at(se).at(i)->SetMarkerColor(i+1);
			yRapidity.at(se).at(i)->SetMarkerStyle(20+i);
			yRapidity.at(se).at(i)->SetMarkerSize(1);
			yRapidity.at(se).at(i)->SetLineWidth(2);
			yRapidity.at(se).at(i)->GetYaxis()->SetRangeUser(-0.6,0.6);
		}
		hRapidity.at(i)->SetLineColor(i+1);
		hRapidity.at(i)->SetMarkerColor(i+1);
		hRapidity.at(i)->SetMarkerStyle(20+i);
		hRapidity.at(i)->SetMarkerSize(1);
		hRapidity.at(i)->SetLineWidth(2);
		hRapidity.at(i)->GetYaxis()->SetRangeUser(-0.6,0.4);
	}
	canvas.push_back( new TCanvas( "canv1", "", 1500, 400 ) );
	canvas.push_back( new TCanvas( "canv2", "", 1500, 400 ) );
	for( int i=0; i<2; i++ )
	{
		canvas.at(i)->Divide(3,1);
		for(int j=0; j<3; j++)
		{
			canvas.at(i)->cd(j+1);
			stack.at(i).at(j)->Draw();
			gPad->BuildLegend(0.1,0.65,0.45,0.9);
		}
		canvas.at(i)->SaveAs( sFileName+Form("_y_%i",i)+".png" );
	}
	canvas.push_back( new TCanvas( "canv5", "", 900, 750 ) );
	canvas.back()->cd();
	stack.back().back()->Draw("NOSTACK");
	gPad->BuildLegend(0.1,0.75,0.5,0.9);
	canvas.back()->SaveAs( sFileName+"_y_2.png" );
}

void Flow3SE::SaveHistogramsToRootFile(TString sFileName)
{
	auto file = new TFile(sFileName+".root");
	file->cd();
	for(int se=0; se<3; se++)
	{
		for( auto histo : xRapidity.at(se) )
			histo->Write();
		for( auto histo : yRapidity.at(se) )
			histo->Write();
		for( auto histo : xPt.at(se) )
			histo->Write();
		for( auto histo : yPt.at(se) )
			histo->Write();
	}
	file->Close();
}