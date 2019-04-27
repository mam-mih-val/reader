#include "Qvector2SE.h"

Qvector2SE::Qvector2SE(DataTreeEvent* _fEvent, Centrality* _centrality)
{
	fEvent = _fEvent;
	fCentrality = _centrality;
	iNumberOfSE = 2;
	for(unsigned int i=0;i<iNumberOfSE;i++)
		fQvector.push_back( TVector2(0.,0.) );
	for( unsigned int i=0;i<iNumberOfSE;i++ )
		fResolution.push_back( TVector2(0.,0.) );
	this->InitHistograms();
}

Qvector2SE::~Qvector2SE()
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
}

void Qvector2SE::InitializeHistograms()
{
	for(unsigned int i=0; i<iNumberOfSE; i++)
	{
		hMeanQx.push_back( new TProfile( Form("MeanQxSE%i",i),";Centrality;Qx", nbins, 1, nbins+1 ) );
		hMeanQy.push_back( new TProfile( Form("MeanQySE%i",i),";Centrality;Qy", nbins, 1, nbins+1 ) );
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
	auto nbins = fCentrality->GetNumClasses();
	hCorrelation.push_back( new TProfile("Qx_{a}Qx_{b}", ";Centrality;Qx_{1}Qx_{2}", nbins, 1, nbins+1) ); // 0
	hCorrelation.push_back( new TProfile("Qy_{a}Qy_{b}", ";Centrality;Qy_{1}Qy_{2}", nbins, 1, nbins+1) ); // 1
	hCorrelation.push_back( new TProfile("Qx_{a}Qy_{b}", ";Centrality;Qx_{1}Qy_{2}", nbins, 1, nbins+1) ); // 2
	hCorrelation.push_back( new TProfile("Qy_{a}Qx_{b}", ";Centrality;Qy_{1}Qx_{2}", nbins, 1, nbins+1) ); // 3
	for(int i=0; i< nbins; i++)
	{
		hDeltaPsiEP.push_back( new TH1F( Form("PsiA_PsiB_centr_%i",i), Form("centrality %i;#Psi_{a}-#Psi_{b};counts",i), 100,0,3.1415) );
	}
	hMeanCosine = new TProfile( "MeanCosine",";centrality;<cos(#Psi_1^{a}-#Psi_{1}^{b})>", nbins, 1, nbins+1 );
}

void Qvector2SE::Estimate()
{
	for( auto &vector : fQvector )
		vector.Set( 0., 0. );
	Int_t iNPSDModules = fEvent->GetNPSDModules();
    DataTreePSDModule* fModule;
	vector<DataTreePSDModule*> vModules;
	for(unsigned int i=0; i<iNPSDModules; i++)
	{
		fModule = fEvent->GetPSDModule(i);
		if( fModule->GetId() <= 0 )
			continue;
		vModules.push_back(fModule);
	}
	random_shuffle( vModules.begin(),vModules.end() );
	vector<double> fSumCharge{0.,0.};
	vector<int> fHitsInSE{0,0};
	int p=1;
	for( auto &module : vModules )
	{
		TVector2 fAddition;
		fAddition.SetMagPhi( module->GetChargeZ(), module->GetPhi() );
		fQvector.at(p%2)+=fAddition;
		fSumCharge.at(p%2)+=module->GetChargeZ();
		fHitsInSE.at(p%2)+=1;
		p++;
	}
	for(unsigned int i=0; i<2; i++)
	{
		if( fSumCharge.at(i) < 1 )
		{
			fQvector.at(i).Set(-999.,-999.);
			continue;
		}
		fQvector.at(i) /= fSumCharge.at(i);
		//if( fQvector.at(i).Mod() >= 1 )
		//	cout << "SE: " << i+1 << " Qx=" << fQvector.at(i).X() << " Qy=" << fQvector.at(i).Y() << 
		//	" |Q|=" << fQvector.at(i).Mod() << " charge of SE: " << fSumCharge.at(i) << " hits in SE: " << fHitsInSE.at(i) << endl;
	}
}

void Qvector2SE::Compute()
{
	this->Estimate();
	this->Recenter();
	this->EstimateResolution();
}

void Qvector2SE::ComputeCorrections()
{
	this->Estimate();
	for(int i=0;i<iNumberOfSE;i++)
	{
		if( fQvector.at(i).X() < -990. )
			continue;
		hMeanQx.at(i)->Fill( fCentrality->GetCentralityClass(), fQvector.at(i).X() );
		hMeanQy.at(i)->Fill( fCentrality->GetCentralityClass(), fQvector.at(i).Y() );
		hQx.at(i)->Fill( fQvector.at(i).X() );
		hQy.at(i)->Fill( fQvector.at(i).Y() );
		hPsiEP.at(i)->Fill( fQvector.at(i).Phi() );
	}
}

void Qvector2SE::ComputeCorrelations()
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
	hMeanCosine->Fill( fCentrality->GetCentralityClass(), cos( fQvector.at(0).Phi() - fQvector.at(1).Phi() ) );
	if ( fQvector.at(0).X() > -990. && fQvector.at(1).X() > -990. )
	{
		hCorrelation.at(0)->Fill( fCentrality->GetCentralityClass(), ( fQvector.at(0).X() * fQvector.at(1).X() ) );
		hCorrelation.at(1)->Fill( fCentrality->GetCentralityClass(), ( fQvector.at(0).Y() * fQvector.at(1).Y() ) );
		hCorrelation.at(2)->Fill( fCentrality->GetCentralityClass(), ( fQvector.at(0).X() * fQvector.at(1).Y() ) );
		hCorrelation.at(3)->Fill( fCentrality->GetCentralityClass(), ( fQvector.at(0).Y() * fQvector.at(1).X() ) );
	}
}

void Qvector2SE::ComputeResolution()
{
	auto nbins = fCentrality->GetNumClasses();
	hResolutionEP = new TGraph(nbins-1);
	hResolutionEP->SetTitle(";centrality;R_{1}");
	for(int i=0; i<nbins-1;i++)
	{
		auto res = 2. * hMeanCosine->GetBinContent(i+1);
		hResolutionEP->SetPoint( i, (float) (i+1)*5, sqrt(res) );
	}
}

void Qvector2SE::SavePictures(TString sFileName)
{
	cout << "Drawing Qvector histograms" << endl;
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
	if( iNumberOfSE == 2 )
	{
		cCanvas.push_back( new TCanvas("canvas2","Qvectors",2700,2500) ); // Q-vector correlations
		cCanvas.back()->cd();
		hStack.push_back( new THStack("Stack",";Centrality;Correlations") );
		for( auto histo : hCorrelation )
		{
			hStack.back()->Add(histo);
			histo->SetLineWidth(7);
			histo->SetMarkerSize(4);
			histo->SetLineColor(i);
			histo->SetMarkerColor(i);
			histo->SetMarkerStyle(20+i);
			i++;
		}
		hStack.back()->Draw("NOSTACK");
		gPad->BuildLegend(0.1,0.8,0.38,0.9);
	}
	cCanvas.push_back( new TCanvas("MeanCosine","canv",3000,3000) );
	cCanvas.back()->cd();
	hResolutionEP->SetLineColor(0);
	hResolutionEP->SetLineWidth(7);
	hResolutionEP->SetMarkerSize(8);
	hResolutionEP->SetMarkerStyle(20);
	hResolutionEP->SetMarkerColor(1);
	hResolutionEP->Draw();
	i=0;
	cout << "Saving Pictures as PNG" << endl;
	for( auto canvas : cCanvas )
	{
		canvas->SaveAs( "../histograms/"+sPicName+Form("_%i_%i_SE.png",i,iNumberOfSE) );
		i++;
	}
}

void Qvector2SE::EsimateResolution()
{
	auto cbin = fCentrality->GetCentralityClass();
	fResolution = 2. * hMeanCosine->GetBinContent(i+1);
}