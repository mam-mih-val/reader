#include "Flow.h"

Flow::Flow(DataTreeEvent* _event, Centrality* _centrality, Qvector* _Qvector, Selector* _selector, int _pid, int numberOfSE)
{
	fEvent = _event;
	fCentrality = _centrality;
	fQvector = _Qvector;
	fSelector = _selector;
	fNumberOfSE = numberOfSE;
	fPid = _pid;
	for(int i=0; i<fNumberOfSE; i++)
		fFlow.push_back( TVector2{0.,0.} );
	this->InitializeHistograms();
}

Flow::~Flow()
{
	for( auto histo : hFlowX )
		delete histo;
	for( auto histo : hFlowY )
		delete histo;
}

void Flow::InitializeHistograms()
{
	auto nbins = fCentrality->GetNumClasses();
	for(int i=0; i<fNumberOfSE; i++)
	{
		hFlowX.push_back( new TProfile3D( Form("FlowX%i",i),";centrality;rapidy, y;pt, [#frac{GeV}{c}];v1", nbins, 1, nbins, 10, -0.8, 0.8, 10, 0.0, 2.5 ) );
		hFlowY.push_back( new TProfile3D( Form("FlowY%i",i),";centrality;rapidy, y;pt, [#frac{GeV}{c}];v1", nbins, 1, nbins, 10, -0.8, 0.8, 10, 0.0, 2.5 ) );
	}
	hRapidity = new TProfile( "v1 vs y",";rapidity, y; v_{1}", 8,-0.8,0.8 );
}

void Flow::_Estimae3SE(int trackIdx)
{
	double BETA = sqrt( 1.0 - 0.938*0.938/1.23/1.23 );
	DataTreeTrack* track = fEvent->GetVertexTrack(trackIdx);
	TLorentzVector momentum = track->GetMomentum();
	TVector3 b; b.SetXYZ(0,0,-BETA);
	momentum.Boost(b);
	for( int j=0; j<fNumberOfSE; j++ )
	{
		fFlow.at(j).SetMagPhi( 1., momentum.Phi() );
		if( fQvector->X(j) < -990.0 )
			continue;
		auto flowX = fFlow.at(j).X()*fQvector->X(j) / fQvector->Resolution(j).X();
		// cout << "Vx = " << fFlow.at(j).X()*fQvector->X(j) << " Vy = " << fFlow.at(j).Y()*fQvector->Y(j) << endl;
		// auto flowX = fFlow.at(j).X()*fQvector->X(j);
		auto flowY = fFlow.at(j).Y()*fQvector->Y(j);
		if(flowX==flowX)
			hRapidity->Fill( momentum.Rapidity(), flowX );
		// if( flowX == flowX )
		// 	hFlowX.at(j)->Fill( fCentrality->GetCentralityClass(), momentum.Rapidity(), momentum.Pt(), flowX );
		// // if( flowY == flowY )
		// 	hFlowY.at(j)->Fill( fCentrality->GetCentralityClass(), momentum.Rapidity(), momentum.Pt(), flowY );
	}
}

void Flow::Estimate()
{
	int numberOfTracks = fEvent->GetNVertexTracks();
	DataTreeTrack* track;
	fQvector->Compute();
	for( int i=0; i<numberOfTracks; i++ )
	{
		track = fEvent->GetVertexTrack(i);
		if( !fSelector->IsCorrectTrack(i) )
			continue;
		if( track->GetPdgId() != fPid )
			continue;
		this->_Estimae3SE(i);
	}
}

void Flow::SavePictures(TString sFileName)
{
	cout << "Saving pictures as PNG" << endl;
	vector<TCanvas*> canvas;
	canvas.push_back( new TCanvas("canv","Flow",2000,2000) );
	vector<THStack*> stack;
	stack.push_back(new THStack("Stack",";Centrality;Correlations"));
	// stack.back()->Add(hRapidity);
	// stack.back()->Add(hFlowY.at(0)->ProjectionY());
	canvas.back()->cd();
	hRapidity->Draw();
	// stack.back()->Draw("PADS");
	// hFlowX.at(0)->Draw();
	int i=0;

	// canvas.back()->SaveAs( "../histograms/"+sFileName+".png" );
	for( auto canv : canvas )
	{
		canv->SaveAs( "../histograms/"+sFileName+Form("_%i_%i_SE.png",i,fNumberOfSE) );
		i++;
	}
}

void Flow::SaveHistogramsToRootFile(TString sFileName)
{
	auto file = new TFile("Flow.root","recreate");
	for( auto histo : hFlowX )
		histo->Write();
	for( auto histo : hFlowY )
		histo->Write();

	file->Close();	
}