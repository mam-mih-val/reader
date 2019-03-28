#include "Qvector.h"

Qvector::Qvector(DataTreeEvent* _fEvent, Centrality* _centrality, unsigned int NumSE)
{
	fEvent = _fEvent;
	fCentrality = _centrality;
	iNumberOfSE = NumSE;
	for(unsigned int i=0;i<iNumberOfSE;i++)
		fQvector.push_back( TVector2(0.,0.) );
	this->InitHistograms();
}

Qvector::~Qvector()
{
	for( auto histo : hMeanQx )
		delete histo;
	for( auto histo : hMeanQy )
		delete histo;
	for( auto histo : hCorrelation )
		delete histo;
	for( auto histo : hResolution )
		delete histo;
	for( auto histo : hQx )
		delete histo;
	for( auto histo : hQy )
		delete histo;
	for( auto histo : hPsiEP )
		delete histo;
	for( auto histo : hDeltaPsiEP )
		delete histo;
}

void Qvector::InitHistograms()
{
	cout << "Initialization of Qvector histograms" << endl;
	auto nbins = fCentrality->GetNumClasses();
	for(unsigned int i=0; i<iNumberOfSE; i++)
	{
		hMeanQx.push_back( new TProfile( Form("MeanQxSE%i",i),";Centrality;Qx", nbins, 0, nbins ) );
		hMeanQy.push_back( new TProfile( Form("MeanQySE%i",i),";Centrality;Qy", nbins, 0, nbins ) );
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
	hResolution.push_back( new TProfile("MeanCosine",";centrality class;<cos(#Psi_{a}-#Psi_{b})>",nbins,0,nbins) ); // 0 
	if( iNumberOfSE == 2 )
	{
		hCorrelation.push_back( new TProfile("Qx_{a}Qx_{b}", ";Centrality;Qx_{1}Qx_{2}", nbins, 0, nbins) ); // 0
		hCorrelation.push_back( new TProfile("Qy_{a}Qy_{b}", ";Centrality;Qy_{1}Qy_{2}", nbins, 0, nbins) ); // 1
		hCorrelation.push_back( new TProfile("Qx_{a}Qy_{b}", ";Centrality;Qx_{1}Qy_{2}", nbins, 0, nbins) ); // 2
		hCorrelation.push_back( new TProfile("Qy_{a}Qx_{b}", ";Centrality;Qy_{1}Qx_{2}", nbins, 0, nbins) ); // 3

		for(int i=0; i< nbins; i++)
		{
			hDeltaPsiEP.push_back( new TH1F( Form("PsiA_PsiB_centr_%i",i), Form("centrality %i;#Psi_{a}-#Psi_{b};counts",i), 100,0,3.1415) );
		}
	}
	if( iNumberOfSE == 3 )
	{
		hCorrelation.push_back( new TProfile("Qx_{b}Qx_{c}", ";Centrality;Qx_{b}Qx_{c}", nbins, 0, nbins) ); // 0
		hCorrelation.push_back( new TProfile("Qy_{b}Qy_{c}", ";Centrality;Qy_{b}Qy_{c}", nbins, 0, nbins) ); // 1
		hCorrelation.push_back( new TProfile("Qx_{b}Qy_{c}", ";Centrality;Qx_{b}Qx_{c}", nbins, 0, nbins) ); // 2
		hCorrelation.push_back( new TProfile("Qy_{b}Qx_{c}", ";Centrality;Qy_{b}Qx_{c}", nbins, 0, nbins) ); // 3
	
		hCorrelation.push_back( new TProfile("Qx_{a}Qx_{c}", ";Centrality;Qx_{a}Qx_{c}", nbins, 0, nbins) ); // 4
		hCorrelation.push_back( new TProfile("Qy_{a}Qy_{c}", ";Centrality;Qy_{a}Qy_{c}", nbins, 0, nbins) ); // 5
		hCorrelation.push_back( new TProfile("Qx_{a}Qy_{c}", ";Centrality;Qx_{a}Qx_{c}", nbins, 0, nbins) ); // 6
		hCorrelation.push_back( new TProfile("Qy_{a}Qx_{c}", ";Centrality;Qy_{a}Qx_{c}", nbins, 0, nbins) ); // 7
	
		hCorrelation.push_back( new TProfile("Qx_{a}Qx_{b}", ";Centrality;Qx_{a}Qx_{b}", nbins, 0, nbins) ); // 8
		hCorrelation.push_back( new TProfile("Qy_{a}Qy_{b}", ";Centrality;Qy_{a}Qy_{b}", nbins, 0, nbins) ); // 9
		hCorrelation.push_back( new TProfile("Qx_{a}Qy_{b}", ";Centrality;Qx_{a}Qx_{b}", nbins, 0, nbins) ); // 10
		hCorrelation.push_back( new TProfile("Qy_{a}Qx_{b}", ";Centrality;Qy_{a}Qx_{b}", nbins, 0, nbins) ); // 11
		
		
	}
}

void Qvector::FillCorrections()
{
	if( iNumberOfSE == 2 )
		this->Estimate2SE();
	if( iNumberOfSE == 3 )
		this->Estimate3SE();
	
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

void Qvector::Estimate()
{
	if( iNumberOfSE == 2 )
		this->Estimate2SE();
	if( iNumberOfSE == 3 )
		this->Estimate3SE();

	//Float_t fCentrality = fEvent->GetCentrality();
	int iCentralityBin = (int) fCentrality->GetCentralityClass();
	for(unsigned int i=0;i<fQvector.size();i++)
	{
		if( fQvector.at(i).X() < -990. )
			continue;
		TVector2 fCorrVec;
		fCorrVec.Set( hMeanQx.at(i)->GetBinContent(iCentralityBin), hMeanQy.at(i)->GetBinContent(iCentralityBin) );
		fQvector.at(i)-=fCorrVec;
		hQx[iNumberOfSE+i]->Fill( fQvector.at(i).X() );
		hQy[iNumberOfSE+i]->Fill( fQvector.at(i).Y() );
		hPsiEP[iNumberOfSE+i]->Fill( fQvector.at(i).Phi() );
	}
	this->FillResolutionProfile();
	if ( iNumberOfSE == 2 && fQvector.at(0).X() > -990. && fQvector.at(1).X() > -990. )
	{
		hDeltaPsiEP.at( (int)fCentrality->GetCentralityClass() )->Fill( acos( cos(fQvector.at(0).Phi() - fQvector.at(1).Phi()) ) );
		//hPsiEP.back()->Fill( acos( cos(fQvector.at(0).Phi() - fQvector.at(1).Phi()) ) );
		hCorrelation.at(0)->Fill( fCentrality->GetCentralityClass(), ( fQvector.at(0).X() * fQvector.at(1).X() ) );
		hCorrelation.at(1)->Fill( fCentrality->GetCentralityClass(), ( fQvector.at(0).Y() * fQvector.at(1).Y() ) );
		hCorrelation.at(2)->Fill( fCentrality->GetCentralityClass(), ( fQvector.at(0).X() * fQvector.at(1).Y() ) );
		hCorrelation.at(3)->Fill( fCentrality->GetCentralityClass(), ( fQvector.at(0).Y() * fQvector.at(1).X() ) );
	}
	if ( iNumberOfSE == 3 && fQvector.at(1).X() > -990. && fQvector.at(2).X() > -990. )
	{
		hCorrelation.at(0)->Fill( fCentrality->GetCentralityClass(), ( fQvector.at(1).X() * fQvector.at(2).X() ) );
		hCorrelation.at(1)->Fill( fCentrality->GetCentralityClass(), ( fQvector.at(1).Y() * fQvector.at(2).Y() ) );
		hCorrelation.at(2)->Fill( fCentrality->GetCentralityClass(), ( fQvector.at(1).X() * fQvector.at(2).Y() ) );
		hCorrelation.at(3)->Fill( fCentrality->GetCentralityClass(), ( fQvector.at(1).Y() * fQvector.at(2).X() ) );
	}
	if ( iNumberOfSE == 3 && fQvector.at(0).X() > -990. && fQvector.at(2).X() > -990. )
	{
		hCorrelation.at(4)->Fill( fCentrality->GetCentralityClass(), ( fQvector.at(0).X() * fQvector.at(2).X() ) );
		hCorrelation.at(5)->Fill( fCentrality->GetCentralityClass(), ( fQvector.at(0).Y() * fQvector.at(2).Y() ) );
		hCorrelation.at(6)->Fill( fCentrality->GetCentralityClass(), ( fQvector.at(0).X() * fQvector.at(2).Y() ) );
		hCorrelation.at(7)->Fill( fCentrality->GetCentralityClass(), ( fQvector.at(0).Y() * fQvector.at(2).X() ) );
	}
	if ( iNumberOfSE == 3 && fQvector.at(0).X() > -990. && fQvector.at(1).X() > -990. )
	{
		hCorrelation.at(8)->Fill( fCentrality->GetCentralityClass(), ( fQvector.at(0).X() * fQvector.at(1).X() ) );
		hCorrelation.at(9)->Fill( fCentrality->GetCentralityClass(), ( fQvector.at(0).Y() * fQvector.at(1).Y() ) );
		hCorrelation.at(10)->Fill( fCentrality->GetCentralityClass(), ( fQvector.at(0).X() * fQvector.at(1).Y() ) );
		hCorrelation.at(11)->Fill( fCentrality->GetCentralityClass(), ( fQvector.at(0).Y() * fQvector.at(1).X() ) );
	}
}

void Qvector::SaveHistograms(TString sPicName)
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
		// Qy-distribution is draqwing
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
	if( iNumberOfSE == 3 )
	{
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
				hStack.back()->Add( hCorrelation.at(j) );
				hCorrelation.at(j)->SetLineWidth(5);
				hCorrelation.at(j)->SetMarkerSize(2);
				hCorrelation.at(j)->SetLineColor(n);
				hCorrelation.at(j)->SetMarkerColor(n);
				hCorrelation.at(j)->SetMarkerStyle(20+n);
				n++;
			}
			hStack.back()->Draw("NOSTACK");
			gPad->BuildLegend(0.1,0.62,0.92,0.9);
		}
	}

	cout << "Saving Pictures as PNG" << endl;
	cCanvas.push_back( new TCanvas("Psi","canv",4000,4000) );
	cCanvas.back()->cd();
	THStack* hStack1 = new THStack("Stack"," ");
	for( auto histo : hDeltaPsiEP )
	{
		histo->SetLineWidth(3);
		hStack1->Add(histo);
	}
	hStack1->Draw("PADS");
	cCanvas.push_back( new TCanvas("Resolution","canv",4000,3500) );
	cCanvas.back()->cd();
	hResolution.back()->SetLineWidth(7);
	hResolution.back()->SetMarkerSize(5);
	hResolution.back()->SetLineColor(9);
	hResolution.back()->SetMarkerColor(9);
	hResolution.back()->SetMarkerStyle(20);
	hResolution.back()->Draw();
	i=0;
	for( auto canvas : cCanvas )
	{
		canvas->SaveAs( "../histograms/"+sPicName+Form("_%i_%i_SE.png",i,iNumberOfSE) );
		i++;
	}
}

void Qvector::Estimate2SE()
{
	for( auto vector = begin(fQvector); vector != end(fQvector); ++vector )
		vector->Set( 0., 0. );
	Int_t iNPSDModules = fEvent->GetNPSDModules();
    DataTreePSDModule* fModule;
	vector<DataTreePSDModule*> vModules;
	for(int i=0; i<iNPSDModules;i++)
	{
		fModule = fEvent->GetPSDModule(i);
		if( fModule->GetId() <= 0 )
			continue;
		vModules.push_back(fModule);
	}
	TVector2 fAddition;
	random_shuffle( vModules.begin(),vModules.end() );
	vector<double> fSumCharge(0.);
	vector<int> fHitsInSE(0);
	int p=0;
	for( auto module = begin(vModules); module != end(vModules); ++module )
	{
		double charge = (*module)->GetEnergy();
		double phi = (*module)->GetPhi();
		fAddition.SetMagPhi( charge, phi );
		fQvector.at( p%2 )+=fAddition;
		fSumCharge.at(p%2)+=charge;
		fHitsInSE.at(p%2)++;
		p++;
	}
	for(unsigned int i=0; i<2; i++)
	{
		if( fSumCharge.at(i) < 70. )
		{
			fQvector.at(i).Set(-999.,-999.);
			continue;
		}
		if( fHitsInSE.at(i) < 3 )
		{
			fQvector.at(i).Set( -999., -999. );
			continue;
		}
		fQvector.at(i)*=( 1. / fSumCharge.at(i) );
		if( fQvector.at(i).Mod() >= 1 )
			cout << "SE: " << i+1 << " Qx=" << fQvector.at(i).X() << " Qy=" << fQvector.at(i).Y() << 
			" |Q|=" << fQvector.at(i).Mod() << " charge of SE: " << fSumCharge.at(i) << " hits in SE: " << fHitsInSE.at(i) << endl;
	}
}

void Qvector::SaveHistogramsToROOTFile(TString sFileName)
{
	TFile* file = new TFile("../histograms/"+sFileName+".root","recreate");
	file->cd();
	cout << "Saving histograms in root file" << endl;
	for( auto histo = begin(hMeanQx); histo != end(hMeanQx); ++histo )
		(*histo)->Write();
	for( auto histo = begin(hMeanQy); histo != end(hMeanQy); ++histo )
		(*histo)->Write();
	for( auto histo = begin(hQx); histo != end(hQx); ++histo )
		(*histo)->Write();
	for( auto histo = begin(hQy); histo != end(hQy); ++histo )
		(*histo)->Write();
	for( auto histo = begin(hPsiEP); histo != end(hPsiEP); ++histo )
		(*histo)->Write();
	for( auto histo : hCorrelation )
		histo->Write();
	file->Close();
}

void Qvector::Estimate3SE()
{
	for( auto &vector : fQvector )
		vector.Set( 0., 0. );
	Int_t iNPSDModules = fEvent->GetNPSDModules();
    array<vector<DataTreePSDModule*>, 3> SubEvent;
	for(int i=0; i<iNPSDModules; i++)
	{
		if( fEvent->GetPSDModule(i)->GetId() > 0 && fEvent->GetPSDModule(i)->GetId() < 6 )
			SubEvent.at(0).push_back( fEvent->GetPSDModule(i) );
		if( fEvent->GetPSDModule(i)->GetId() >= 6 && fEvent->GetPSDModule(i)->GetId() < 8 )
			SubEvent.at(1).push_back( fEvent->GetPSDModule(i) );
		if( fEvent->GetPSDModule(i)->GetId() >= 8 && fEvent->GetPSDModule(i)->GetId() < 11 )
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
		for( auto module : SubEvent.at(i) )
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

void Qvector::FillResolutionProfile()
{
	if( fQvector.at(0).X() > -990.0 && fQvector.at(1).X() > -990.0  )
		hResolution.back()->Fill( fCentrality->GetCentralityClass(), 2*cos( fQvector.at(0).Phi() - fQvector.at(1).Phi() ) );
}

float Qvector::GetResolution()
{
	auto bin = (int)hResolution.back()->FindBin( (float)fCentrality->GetCentralityClass() );
	return hResolution.back()->GetBinContent(bin);
}