#include "Qvector.h"

Qvector::Qvector()
{
	iNumberOfSE=2;
	for(unsigned int i=0;i<iNumberOfSE;i++)
		fQvector.push_back( TVector2(0.,0.) );
	this->LoadCentralityPercentile("~/centrality_epcorr_apr12_gen8_2018_07.root");
	this->InitHistograms();
}

void Qvector::InitHistograms()
{
	cout << "Initialization of Qvector histograms" << endl;
	auto nbins = hCentralityPercentile->GetNbinsX();
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
	hPsiEP.push_back( new TH1F("PsiEP",";#Psi_{a}-#Psi_{b};counts",100,0,3.1415) );
	if( iNumberOfSE == 2 )
	{
		hCorrelation.push_back( new TProfile("Qx_{a}Qx_{b}", ";Centrality;Qx_{1}Qx_{2}", nbins, 0, nbins) ); // 0
		hCorrelation.push_back( new TProfile("Qy_{a}Qy_{b}", ";Centrality;Qy_{1}Qy_{2}", nbins, 0, nbins) ); // 1
		hCorrelation.push_back( new TProfile("Qx_{a}Qy_{b}", ";Centrality;Qx_{1}Qy_{2}", nbins, 0, nbins) ); // 2
		hCorrelation.push_back( new TProfile("Qy_{a}Qx_{b}", ";Centrality;Qy_{1}Qx_{2}", nbins, 0, nbins) ); // 3
	}
	if( iNumberOfSE == 3 )
	{
		hCorrelation.push_back( new TProfile("Qx_{a}Qx_{b}", ";Centrality;Qx_{a}Qx_{b}", nbins, 0, nbins) ); // 0
		hCorrelation.push_back( new TProfile("Qx_{b}Qx_{c}", ";Centrality;Qx_{b}Qx_{c}", nbins, 0, nbins) ); // 1
		hCorrelation.push_back( new TProfile("Qx_{a}Qx_{c}", ";Centrality;Qx_{a}Qx_{c}", nbins, 0, nbins) ); // 2
		
		hCorrelation.push_back( new TProfile("Qy_{a}Qy_{b}", ";Centrality;Qy_{a}Qy_{b}", nbins, 0, nbins) ); // 3
		hCorrelation.push_back( new TProfile("Qy_{b}Qy_{c}", ";Centrality;Qy_{b}Qy_{c}", nbins, 0, nbins) ); // 4
		hCorrelation.push_back( new TProfile("Qy_{a}Qy_{c}", ";Centrality;Qy_{a}Qy_{c}", nbins, 0, nbins) ); // 5
		
		hCorrelation.push_back( new TProfile("Qx_{a}Qy_{b}", ";Centrality;Qx_{a}Qx_{b}", nbins, 0, nbins) ); // 6
		hCorrelation.push_back( new TProfile("Qy_{a}Qx_{b}", ";Centrality;Qy_{a}Qx_{b}", nbins, 0, nbins) ); // 7
		
		hCorrelation.push_back( new TProfile("Qx_{b}Qy_{c}", ";Centrality;Qx_{b}Qx_{c}", nbins, 0, nbins) ); // 8
		hCorrelation.push_back( new TProfile("Qy_{b}Qx_{c}", ";Centrality;Qy_{b}Qx_{c}", nbins, 0, nbins) ); // 9
		
		hCorrelation.push_back( new TProfile("Qx_{a}Qy_{c}", ";Centrality;Qx_{a}Qx_{c}", nbins, 0, nbins) ); // 10
		hCorrelation.push_back( new TProfile("Qy_{a}Qx_{c}", ";Centrality;Qy_{a}Qx_{c}", nbins, 0, nbins) ); // 11
	}
}

void Qvector::FillCorrections(DataTreeEvent* fEvent)
{
	if( iNumberOfSE == 2 )
		this->Estimate2SE(fEvent);
	if( iNumberOfSE == 3 )
		this->Estimate3SE(fEvent);
	
	for(int i=0;i<iNumberOfSE;i++)
	{
		if( fQvector.at(i).X() < -990. )
			continue;
		hMeanQx.at(i)->Fill( this->GetCentralityClass(fEvent), fQvector.at(i).X() );
		hMeanQy.at(i)->Fill( this->GetCentralityClass(fEvent), fQvector.at(i).Y() );
		hQx.at(i)->Fill( fQvector.at(i).X() );
		hQy.at(i)->Fill( fQvector.at(i).Y() );
		hPsiEP.at(i)->Fill( fQvector.at(i).Phi() );
	}
}

void Qvector::Estimate(DataTreeEvent* fEvent)
{
	if( iNumberOfSE == 2 )
		this->Estimate2SE(fEvent);
	if( iNumberOfSE == 3 )
		this->Estimate3SE(fEvent);

	//Float_t fCentrality = fEvent->GetCentrality();
	int iCentralityBin = (int) this->GetCentralityClass(fEvent);
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
	if ( iNumberOfSE == 2 && fQvector.at(0).X() > -990. && fQvector.at(1).X() > -990. )
	{
		hPsiEP.back()->Fill( acos( cos(fQvector.at(0).Phi() - fQvector.at(1).Phi()) ) );
		hCorrelation.at(0)->Fill( this->GetCentralityClass(fEvent), ( fQvector.at(0).X() * fQvector.at(1).X() ) );
		hCorrelation.at(1)->Fill( this->GetCentralityClass(fEvent), ( fQvector.at(0).Y() * fQvector.at(1).Y() ) );
		hCorrelation.at(2)->Fill( this->GetCentralityClass(fEvent), ( fQvector.at(0).X() * fQvector.at(1).Y() ) );
		hCorrelation.at(3)->Fill( this->GetCentralityClass(fEvent), ( fQvector.at(0).Y() * fQvector.at(1).X() ) );
	}
	if ( iNumberOfSE == 3 && fQvector.at(0).X() > -990. && fQvector.at(1).X() > -990. && fQvector.at(2).X() > -990. )
	{
		hCorrelation.at(0)->Fill( this->GetCentralityClass(fEvent), ( fQvector.at(0).X() * fQvector.at(1).X() ) );
		hCorrelation.at(1)->Fill( this->GetCentralityClass(fEvent), ( fQvector.at(1).X() * fQvector.at(2).X() ) );
		hCorrelation.at(2)->Fill( this->GetCentralityClass(fEvent), ( fQvector.at(0).X() * fQvector.at(2).X() ) );
		
		hCorrelation.at(3)->Fill( this->GetCentralityClass(fEvent), ( fQvector.at(0).Y() * fQvector.at(1).Y() ) );
		hCorrelation.at(4)->Fill( this->GetCentralityClass(fEvent), ( fQvector.at(1).Y() * fQvector.at(2).Y() ) );
		hCorrelation.at(5)->Fill( this->GetCentralityClass(fEvent), ( fQvector.at(0).Y() * fQvector.at(2).Y() ) );
		
		hCorrelation.at(6)->Fill( this->GetCentralityClass(fEvent), ( fQvector.at(0).X() * fQvector.at(1).Y() ) );
		hCorrelation.at(7)->Fill( this->GetCentralityClass(fEvent), ( fQvector.at(0).Y() * fQvector.at(1).X() ) );
		
		hCorrelation.at(8)->Fill( this->GetCentralityClass(fEvent), ( fQvector.at(1).X() * fQvector.at(2).Y() ) );
		hCorrelation.at(9)->Fill( this->GetCentralityClass(fEvent), ( fQvector.at(1).Y() * fQvector.at(2).X() ) );
		
		hCorrelation.at(10)->Fill( this->GetCentralityClass(fEvent), ( fQvector.at(0).X() * fQvector.at(2).Y() ) );
		hCorrelation.at(11)->Fill( this->GetCentralityClass(fEvent), ( fQvector.at(0).Y() * fQvector.at(2).X() ) );
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
	cCanvas.push_back( new TCanvas("canvas2","Qvectors",4000,2500) ); // Q-vector correlations
	cCanvas.back()->cd();
	int i=1;
	THStack* hStack = new THStack("Stack",";Centrality;Correlations");
	
	for( auto histo : hCorrelation )
	{
		hStack->Add(histo);
		histo->SetLineWidth(7);
		histo->SetMarkerSize(4);
		histo->SetLineColor(i);
		histo->SetMarkerColor(i);
		histo->SetMarkerStyle(20+i);
		//histo->GetYaxis()->SetRangeUser(-0.01, 0.07);
		i++;
	}
	hStack->Draw("NOSTACK");
	gPad->BuildLegend(0.1,0.8,0.38,0.9);
	cout << "Saving Pictures as PNG" << endl;
	cCanvas.push_back( new TCanvas("Psi","canv",1500,1000) );
	cCanvas.back()->cd();
	hPsiEP.back()->SetLineWidth(5);
	hPsiEP.back()->Draw();
	i=0;
	for( auto canvas : cCanvas )
	{
		canvas->SaveAs( "../histograms/"+sPicName+Form("_%i.png",i) );
		i++;
	}
}

void Qvector::Estimate2SE(DataTreeEvent* fEvent)
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
	int p=0;
	vector<double> fSumCharge;
	vector<int> fHitsInSE;
	for( int i=0; i<iNumberOfSE; i++ )
	{
		fSumCharge.push_back(0.);
		fHitsInSE.push_back(0);
	}
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
		if( fSumCharge.at(i) < 0.1 )
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

void Qvector::Estimate3SE(DataTreeEvent* fEvent)
{
	for( auto vector : fQvector )
		vector.Set( 0., 0. );
	Int_t iNPSDModules = fEvent->GetNPSDModules();
    DataTreePSDModule* fModule;
	vector<float> fSumCharge;
	vector<int> fHitsInSE;
	for( int i=0; i<iNumberOfSE; i++ )
	{
		fSumCharge.push_back(0.);
		fHitsInSE.push_back(0);
	}
	TVector2 fAddition;
	for(int i=0; i<iNPSDModules;i++)
	{
		fModule = fEvent->GetPSDModule(i);
		if( fModule->GetId() <= 0 )
			continue;
		double charge = fModule->GetEnergy();
		double phi = fModule->GetPhi();
		fAddition.SetMagPhi( charge, phi );
		if( fModule->GetId() <= 5 )
		{
			fQvector.at(0)+=fAddition;
			fSumCharge.at(0)+=charge;
			fHitsInSE.at(0)++;
		}
		if( fModule->GetId() == 6 || fModule->GetId() == 7 )
		{
			fQvector.at(1)+=fAddition;
			fSumCharge.at(1)+=charge;
			fHitsInSE.at(1)++;
		}
		if( fModule->GetId() >= 8 )
		{
			fQvector.at(2)+=fAddition;
			fSumCharge.at(2)+=charge;
			fHitsInSE.at(2)++;
		}
	}
	for( unsigned int i=0; i<fQvector.size(); i++ )
	{
		if( fSumCharge.at(i) < 0.1 )
		{
			fQvector.at(i).Set( -999., -999. );
			continue;
		}
		if( fHitsInSE.at(i) < 3 )
		{
			fQvector.at(i).Set( -999., -999. );
			continue;
		}
		fQvector.at(i) *= ( 1. / fSumCharge.at(i) );
		if( fQvector.at(i).Mod() >= 1 )
			cout << "SE: " << i+1 << " Qx=" << fQvector.at(i).X() << " Qy=" << fQvector.at(i).Y() << 
			" |Q|=" << fQvector.at(i).Mod() << " charge of SE: " << fSumCharge.at(i) << " hits in SE: " << fHitsInSE.at(i) << endl;
	}
}

float Qvector::GetCentralityClass(DataTreeEvent* fEvent)
{ 	
	auto TOFRPChits = fEvent->GetCentralityEstimator(HADES_constants::kNhitsTOF_cut) + fEvent->GetCentralityEstimator(HADES_constants::kNhitsRPC_cut);
	auto bin = hCentralityPercentile->FindBin(TOFRPChits);
	return hCentralityPercentile->GetBinContent(bin);
}

void Qvector::LoadCentralityPercentile(TString FileName)
{
	auto file = new TFile(FileName);
	if ( file->IsOpen() )
	{
		cout << "Centrality file loaded" << endl;
		file->cd("/Centrality/");
		hCentralityPercentile = (TH1F*) file->Get("/Centrality/TOFRPCtot_5pc_fixedCuts");
		cout << hCentralityPercentile->GetNbinsX() << " centrality classes" << endl;
	}
	else
	{
		cout << "Couldn't open the file" << endl;
		return;
	}
}