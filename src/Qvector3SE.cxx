#include "Qvector3SE.h"

Qvector3SE::Qvector3SE(DataTreeEvent* _fEvent, Centrality* _centrality)
{
	fEvent = _fEvent;
	fCentrality = _centrality;
	iNumberOfSE = 3;
	for(unsigned int i=0;i<iNumberOfSE;i++)
		fQvector.push_back( TVector2(0.,0.) );
	for( unsigned int i=0;i<iNumberOfSE;i++ )
		fResolution.push_back( TVector2(0.,0.) );
	this->InitializeHistograms();
}

Qvector3SE::~Qvector3SE()
{
	for( auto histo : hMeanQx )
		delete histo;
	for( auto histo : hMeanQy )
		delete histo;
	for( auto histo : hCorrelation )
		delete histo;
	for( auto histo : hQx )
		delete histo;
	for( auto histo : hQy )
		delete histo;
	for( auto histo : hPsiEP )
		delete histo;
	for( auto histo : hResolutionX )
		delete histo;
	for( auto histo : hResolutionY )
		delete histo;
}

void Qvector3SE::InitializeHistograms()
{
	for(unsigned int i=0; i<iNumberOfSE; i++)
	{
		hMeanQx.push_back( new TProfile( Form("MeanQxSE%i",i),";Centrality;Qx", 8, 0, 40 ) );
		hMeanQy.push_back( new TProfile( Form("MeanQySE%i",i),";Centrality;Qy", 8, 0, 40 ) );
		hQx.push_back( new TH1F( Form("QxSE%i",i),";Qx;counts",100,-1.5,1.5) );
		hQy.push_back( new TH1F( Form("QySE%i",i),";Qy;counts",100,-1.5,1.5) );
		hPsiEP.push_back( new TH1F( Form("PsiEPSE%i",i),";PsiEP;counts",100,-3.15,3.15) );
	}
	for(unsigned int i=iNumberOfSE; i<2*iNumberOfSE; i++)
	{
		hQx.push_back( new TH1F( Form("QxRecentredSE%i",i-iNumberOfSE),";Qx;counts",100,-1.5,1.5) );
		hQy.push_back( new TH1F( Form("QyRecentredSE%i",i-iNumberOfSE),";Qy;counts",100,-1.5,1.5) );
		hPsiEP.push_back( new TH1F( Form("PsiEPRecentredSE%i",i-iNumberOfSE),";PsiEP;counts",100,0,6.3) );
	}
	hCorrelation.push_back( new TProfile("Qx_{a}Qx_{b}", ";Centrality;Qx_{a}Qx_{b}", 8, 0, 40) ); // 0
	hCorrelation.push_back( new TProfile("Qy_{a}Qy_{b}", ";Centrality;Qy_{a}Qy_{b}", 8, 0, 40) ); // 1
	hCorrelation.push_back( new TProfile("Qx_{a}Qy_{b}", ";Centrality;Qx_{a}Qy_{b}", 8, 0, 40) ); // 2
	hCorrelation.push_back( new TProfile("Qy_{a}Qx_{b}", ";Centrality;Qy_{a}Qx_{b}", 8, 0, 40) ); // 3

	hCorrelation.push_back( new TProfile("Qx_{b}Qx_{c}", ";Centrality;Qx_{b}Qx_{c}", 8, 0, 40) ); // 4
	hCorrelation.push_back( new TProfile("Qy_{b}Qy_{c}", ";Centrality;Qy_{b}Qy_{c}", 8, 0, 40) ); // 5
	hCorrelation.push_back( new TProfile("Qx_{b}Qy_{c}", ";Centrality;Qx_{b}Qx_{c}", 8, 0, 40) ); // 6
	hCorrelation.push_back( new TProfile("Qy_{b}Qx_{c}", ";Centrality;Qy_{b}Qx_{c}", 8, 0, 40) ); // 7

	hCorrelation.push_back( new TProfile("Qx_{a}Qx_{c}", ";Centrality;Qx_{a}Qx_{c}", 8, 0, 40) ); // 8
	hCorrelation.push_back( new TProfile("Qy_{a}Qy_{c}", ";Centrality;Qy_{a}Qy_{c}", 8, 0, 40) ); // 9
	hCorrelation.push_back( new TProfile("Qx_{a}Qy_{c}", ";Centrality;Qx_{a}Qx_{c}", 8, 0, 40) ); // 10
	hCorrelation.push_back( new TProfile("Qy_{a}Qx_{c}", ";Centrality;Qy_{a}Qx_{c}", 8, 0, 40) ); // 11
}

void Qvector3SE::Estimate()
{
	for( auto &vector : fQvector )
		vector.Set( 0., 0. );
	Int_t iNPSDModules = fEvent->GetNPSDModules();
    array<vector<DataTreePSDModule*>, 3> SubEvent;
	for(unsigned int i=0; i<iNPSDModules; i++)
	{
		if( !fEvent->GetPSDModule(i)->HasPassedCuts() )
			continue;
		if( fEvent->GetPSDModule(i)->GetRing() >= 0 && fEvent->GetPSDModule(i)->GetRing() < 5 )
			SubEvent.at(0).push_back( fEvent->GetPSDModule(i) );
		if( fEvent->GetPSDModule(i)->GetRing() >4 && fEvent->GetPSDModule(i)->GetRing() < 7 )
			SubEvent.at(1).push_back( fEvent->GetPSDModule(i) );
		if( fEvent->GetPSDModule(i)->GetRing() >= 7 && fEvent->GetPSDModule(i)->GetRing() <= 9 )
			SubEvent.at(2).push_back( fEvent->GetPSDModule(i) );
	}
	for( int i=0; i<SubEvent.size(); i++ )
	{
		if( SubEvent.at(i).size() < 3 )
		{
			fQvector.at(i).Set(-999.,-999.);
			continue;
		}
		float SumCharge=0;
		for( auto &module : SubEvent.at(i) )
		{
			double charge = module->GetEnergy();
			double phi = module->GetPhi();
			TVector2 add;
			add.SetMagPhi( charge, phi );
			fQvector.at(i)+=add;
			SumCharge+=charge;
		}
		fQvector.at(i)/=SumCharge;
		if( fQvector.at(i).Mod() > 1. )
		 	cout << "SE: " << i+1 << " Qx=" << fQvector.at(i).X() << " Qy=" << fQvector.at(i).Y() << 
			" |Q|=" << fQvector.at(i).Mod() << " charge of SE: " << SumCharge << endl;
	}
}

void Qvector3SE::Compute()
{
	this->Estimate();
	this->Recenter();
	this->EstimateResolution();
}

void Qvector3SE::ComputeCorrections()
{
	this->Estimate();
	for(int i=0;i<iNumberOfSE;i++)
	{
		if( fQvector.at(i).X() < -990. )
			continue;
		hMeanQx.at(i)->Fill( fEvent->GetCentrality(HADES_constants::kNhitsTOF_RPC_cut), fQvector.at(i).X() );
		hMeanQy.at(i)->Fill( fEvent->GetCentrality(HADES_constants::kNhitsTOF_RPC_cut), fQvector.at(i).Y() );
		hQx.at(i)->Fill( fQvector.at(i).X() );
		hQy.at(i)->Fill( fQvector.at(i).Y() );
		hPsiEP.at(i)->Fill( fQvector.at(i).Phi() );
	}
}

void Qvector3SE::ComputeCorrelations()
{
	this->Compute();
	for(unsigned int i=0; i<fQvector.size(); i++)
	{
		if( fQvector.at(i).X() < -990. )
			continue;
		hQx.at(iNumberOfSE+i)->Fill( fQvector.at(i).X() );
		hQy.at(iNumberOfSE+i)->Fill( fQvector.at(i).Y() );
		hPsiEP.at(iNumberOfSE+i)->Fill( fQvector.at(i).Phi() );
	}
	vector<TVector2> fQvector1;
	fQvector1.push_back( fQvector.at(1) );
	fQvector1.push_back( fQvector.at(2) );
	fQvector1.push_back( fQvector.at(0) );
	int bias = 4;
	for( int i=0; i<3; i++ )
	{
		if( fQvector.at(i).X() < -990.0 || fQvector1.at(i).X() < -990.0 )
			continue;
		hCorrelation.at(i*bias+0)->Fill( fEvent->GetCentrality(HADES_constants::kNhitsTOF_RPC_cut), ( fQvector.at(i).X() * fQvector1.at(i).X() ) );
		hCorrelation.at(i*bias+1)->Fill( fEvent->GetCentrality(HADES_constants::kNhitsTOF_RPC_cut), ( fQvector.at(i).Y() * fQvector1.at(i).Y() ) );
		hCorrelation.at(i*bias+2)->Fill( fEvent->GetCentrality(HADES_constants::kNhitsTOF_RPC_cut), ( fQvector.at(i).X() * fQvector1.at(i).Y() ) );
		hCorrelation.at(i*bias+3)->Fill( fEvent->GetCentrality(HADES_constants::kNhitsTOF_RPC_cut), ( fQvector.at(i).Y() * fQvector1.at(i).X() ) );
	}
}

void Qvector3SE::ComputeResolution()
{
	auto nbins = 8;
	for( int i=0; i<iNumberOfSE; i++ )
	{
		hResolutionX.push_back( new TGraph(nbins) );
		hResolutionY.push_back( new TGraph(nbins) );
		hResolutionX.back()->SetTitle(Form("R_{1,x}^{se%i}",i+1));
		hResolutionY.back()->SetTitle(Form("R_{1,y}^{se%i}",i+1));
	}
	
	for(int i=0; i<nbins;i++)
	{
		auto resX = hCorrelation.at(0)->GetBinContent(i+1) * hCorrelation.at(8)->GetBinContent(i+1) / hCorrelation.at(4)->GetBinContent(i+1);
		auto resY = hCorrelation.at(1)->GetBinContent(i+1) * hCorrelation.at(9)->GetBinContent(i+1) / hCorrelation.at(5)->GetBinContent(i+1);
		// cout << resX << endl;
		// cout << resY << endl;
		hResolutionX.at(0)->SetPoint( i, (float) (i+1)*5, sqrt(resX) );
		hResolutionY.at(0)->SetPoint( i, (float) (i+1)*5, sqrt(resY) );
	}
	for(int i=0; i<nbins;i++)
	{
		auto resX = hCorrelation.at(0)->GetBinContent(i+1) * hCorrelation.at(4)->GetBinContent(i+1) / hCorrelation.at(8)->GetBinContent(i+1);
		auto resY = hCorrelation.at(1)->GetBinContent(i+1) * hCorrelation.at(5)->GetBinContent(i+1) / hCorrelation.at(9)->GetBinContent(i+1);
		// cout << resX << endl;
		// cout << resY << endl;
		hResolutionX.at(1)->SetPoint( i, (float) (i+1)*5, sqrt(resX) );
		hResolutionY.at(1)->SetPoint( i, (float) (i+1)*5, sqrt(resY) );
	}
	for(int i=0; i<nbins;i++)
	{
		auto resX = hCorrelation.at(8)->GetBinContent(i+1) * hCorrelation.at(4)->GetBinContent(i+1) / hCorrelation.at(0)->GetBinContent(i+1);
		auto resY = hCorrelation.at(9)->GetBinContent(i+1) * hCorrelation.at(5)->GetBinContent(i+1) / hCorrelation.at(1)->GetBinContent(i+1);
		// cout << resX << endl;
		// cout << resY << endl;
		hResolutionX.at(2)->SetPoint( i, (float) (i+1)*5, sqrt(resX) );
		hResolutionY.at(2)->SetPoint( i, (float) (i+1)*5, sqrt(resY) );
	}
}

void Qvector3SE::SavePictures(TString sFileName)
{
	cout << "Drawing Qvector histograms" << endl;
	if( hResolutionX.size() == 0 )
		this->ComputeResolution();
	vector<TCanvas*> cCanvas;
	vector<TLegend*> legend;
	cCanvas.push_back( new TCanvas("canvas0","Qvectors",4000,2500) );
	cCanvas.back()->Divide(iNumberOfSE,2,0.005,0.0001);
	for(unsigned int i=0;i<iNumberOfSE;i++)
	{
		cCanvas.back()->cd(i+1);
		legend.push_back( new TLegend(0.1,0.8,0.38,0.9) );
		hQx.at(i)->SetLineWidth(5);
		hQx.at(i+iNumberOfSE)->SetLineColor(3);
		hQx.at(i+iNumberOfSE)->SetLineWidth(5);
		legend.back()->AddEntry(hQx.at(i),"Not Recentred");
		legend.back()->AddEntry(hQx.at(i+iNumberOfSE),"Recentred");
		hQx.at(i)->Draw();
		hQx.at(i+iNumberOfSE)->Draw("same");
		legend.back()->Draw();
		// Qy-distribution is drawing
		cCanvas.back()->cd(i+1+iNumberOfSE);
		hQy.at(i)->SetLineWidth(5);
		hQy.at(i+iNumberOfSE)->SetLineColor(3);
		hQy.at(i+iNumberOfSE)->SetLineWidth(5);
		hQy.at(i)->Draw();
		hQy.at(i+iNumberOfSE)->Draw("same");
		legend.back()->Draw();
	}
	vector<THStack*> hStack;
	int i=1;
	int bias = 4;
	cCanvas.push_back( new TCanvas("canvas2","Qvectors",3000,1000) ); // Q-vector correlations
	cCanvas.back()->Divide(3,1);
	for(int k=0; k<3; k++)
	{
		cCanvas.back()->cd(k+1);
		hStack.push_back( new THStack( Form("Stack_%i",k), ";Centrality;Correlations") );
		int n=1;
		for( int j=k*bias; j<bias*(k+1); j++ )
		{
			hStack.back()->Add( hCorrelation.at(j), "X0" );
			hCorrelation.at(j)->SetLineWidth(5);
			hCorrelation.at(j)->SetMarkerSize(2);
			hCorrelation.at(j)->SetLineColor(n);
			hCorrelation.at(j)->SetMarkerColor(n);
			hCorrelation.at(j)->SetMarkerStyle(20+n);
			n++;
		}
		hStack.back()->Draw("NOSTACK");
		gPad->BuildLegend(0.62,0.8,0.9,0.9);
	}
	cCanvas.push_back( new TCanvas("Resolution","canv",3000,3000) );
	auto hStack2 = new TMultiGraph("Stack2",";centrality;Resolution");
	i=1;
	cCanvas.back()->cd();
	for( auto &histo : hResolutionX )
	{
		histo->SetLineColor(0);
		histo->SetMarkerSize(6);
		histo->SetMarkerStyle(19+i);
		histo->SetMarkerColor(i);
		hStack2->Add( histo,"P" );
		i++;
	}
	i=1;
	for( auto &histo : hResolutionY )
	{
		histo->SetLineColor(0);
		histo->SetLineWidth(7);
		histo->SetMarkerSize(8);
		histo->SetMarkerStyle(23+i);
		histo->SetMarkerColor(i);
		hStack2->Add( histo,"P" );
		i++;
	}
	hStack2->Draw("A");
	gPad->BuildLegend(0.62,0.76,0.9,0.9);
	i=0;
	cout << "Saving Pictures as PNG" << endl;
	for( auto canvas : cCanvas )
	{
		canvas->SaveAs( "../histograms/"+sFileName+Form("_%i_%i_SE.png",i,iNumberOfSE) );
		i++;
	}
}

void Qvector3SE::EstimateResolution()
{
	auto cbin = fCentrality->GetCentralityClass();
	auto resX = hCorrelation.at(0)->GetBinContent(cbin) * hCorrelation.at(8)->GetBinContent(cbin) / hCorrelation.at(4)->GetBinContent(cbin);
	auto resY = hCorrelation.at(1)->GetBinContent(cbin) * hCorrelation.at(9)->GetBinContent(cbin) / hCorrelation.at(5)->GetBinContent(cbin);
	fResolution.at(0).Set( sqrt(resX), sqrt(resY) );
	resX = hCorrelation.at(0)->GetBinContent(cbin) * hCorrelation.at(4)->GetBinContent(cbin) / hCorrelation.at(8)->GetBinContent(cbin);
	resY = hCorrelation.at(1)->GetBinContent(cbin) * hCorrelation.at(5)->GetBinContent(cbin) / hCorrelation.at(9)->GetBinContent(cbin);
	fResolution.at(1).Set( sqrt(resX), sqrt(resY) );
	resX = hCorrelation.at(8)->GetBinContent(cbin) * hCorrelation.at(4)->GetBinContent(cbin) / hCorrelation.at(0)->GetBinContent(cbin);
	resY = hCorrelation.at(9)->GetBinContent(cbin) * hCorrelation.at(5)->GetBinContent(cbin) / hCorrelation.at(1)->GetBinContent(cbin);
	fResolution.at(2).Set( sqrt(resX), sqrt(resY) );
}