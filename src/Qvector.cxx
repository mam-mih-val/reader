#include "Qvector.h"

Qvector::Qvector()
{
	iNumberOfSE=2;
	for(int i=0;i<iNumberOfSE;i++)
		fQvector.push_back( TVector2(0.,0.) );
	this->InitHistograms();
}

void Qvector::InitHistograms()
{
	cout << "Initialization of Qvector histograms" << endl;
	
	for(unsigned int i=0; i<iNumberOfSE; i++)
	{
		hMeanQx.push_back( new TProfile( Form("MeanQxSE%i",i),";Centrality;Qx", 10,0,50 ) );
		hMeanQy.push_back(	new TProfile( Form("MeanQySE%i",i),";Centrality;Qy", 10,0,50 ) );
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
	hCorrelation.push_back( new TProfile("Qx_{1}Qx_{2}", ";Centrality;Qx_{1}Qx_{2}", 10, 0, 50) ); // 0
	hCorrelation.push_back( new TProfile("Qy_{1}Qy_{2}", ";Centrality;Qy_{1}Qy_{2}", 10, 0, 50) ); // 1
	hCorrelation.push_back( new TProfile("Qx_{1}Qy_{2}", ";Centrality;Qx_{1}Qy_{2}", 10, 0, 50) ); // 2
	hCorrelation.push_back( new TProfile("Qy_{1}Qx_{2}", ";Centrality;Qy_{1}Qx_{2}", 10, 0, 50) ); // 3
}

void Qvector::FillCorrections(DataTreeEvent* fEvent)
{
	this->Estimate2SE(fEvent);
	for(int i=0;i<iNumberOfSE;i++)
	{
		if( fQvector.at(i).X() == -999. )
			continue;
		hMeanQx.at(i)->Fill( fEvent->GetCentrality(), fQvector.at(i).X() );
		hMeanQy.at(i)->Fill( fEvent->GetCentrality(), fQvector.at(i).Y() );
		hQx.at(i)->Fill( fQvector.at(i).X() );
		hQy.at(i)->Fill( fQvector.at(i).Y() );
		hPsiEP.at(i)->Fill( fQvector.at(i).Phi() );
	}
}

void Qvector::Estimate(DataTreeEvent* fEvent)
{
	this->Estimate2SE(fEvent);
	Float_t fCentrality = fEvent->GetCentrality();
	int iCentralityBin = fCentrality/5+1;
	for(int i=0;i<fQvector.size();i++)
	{
		if( fQvector.at(i).X() == -999. )
			continue;
		TVector2 fCorrVec;
		fCorrVec.Set( hMeanQx.at(i)->GetBinContent(iCentralityBin), hMeanQy.at(i)->GetBinContent(iCentralityBin) );
		fQvector.at(i)-=fCorrVec;
		hQx[iNumberOfSE+i]->Fill( fQvector.at(i).X() );
		hQy[iNumberOfSE+i]->Fill( fQvector.at(i).Y() );
		hPsiEP[iNumberOfSE+i]->Fill( fQvector.at(i).Phi() );
	}
	if( fQvector.at(0).X() != -999. && fQvector.at(1).X() != -999. )
	{
		hCorrelation.at(0)->Fill( fCentrality, ( fQvector.at(0).X() * fQvector.at(1).X() ) );
		hCorrelation.at(1)->Fill( fCentrality, ( fQvector.at(0).Y() * fQvector.at(1).Y() ) );
		hCorrelation.at(2)->Fill( fCentrality, ( fQvector.at(0).X() * fQvector.at(1).Y() ) );
		hCorrelation.at(3)->Fill( fCentrality, ( fQvector.at(0).Y() * fQvector.at(1).X() ) );
	}
}

void Qvector::SaveHistograms(TString sPicName)
{
	cout << "Drawing Qvector histograms" << endl;
	vector<TCanvas*> cCanvas;
	
	vector<TLegend*> legend;
	cCanvas.push_back( new TCanvas("canvas0","Qvectors",4000,2500) );
	cCanvas.back()->Divide(iNumberOfSE,2,0.005,0.0001);
	for(int i=0;i<iNumberOfSE;i++)
	{
		cCanvas.back()->cd(i+1);
		legend.push_back( new TLegend(0.1,0.8,0.38,0.9) );
		hQx.at(i)->SetLineWidth(5);
		hQx.at(i+iNumberOfSE)->SetLineColor(3);
		hQx.at(i+iNumberOfSE)->SetLineWidth(5);
		legend.back()->AddEntry(hQx.at(i),"Not Recentred");
		legend.back()->AddEntry(hQx.at(i+iNumberOfSE),"Not Recentred");
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
	/*
	// Mean Qx vs Centrality is drawing
	cCanvas.push_back( new TCanvas("canvas1","Qvectors",4000,2500) ); // 1, Mean Q-vector component distribution
	cCanvas[1]->Divide(iNumberOfSE,2,0.005,0.0001);
		cCanvas.at(1)->cd(i+1);
		hMeanQx.at(i)->SetLineWidth(5);
		hMeanQx.at(i)->SetLineColor(1);
		hMeanQx.at(i)->SetMarkerSize(4);
		hMeanQx.at(i)->SetMarkerStyle(20);
		hMeanQx.at(i)->Draw();
		// Qy-distribution is draqwing
		cCanvas.at(1)->cd(i+iNumberOfSE+1);
		hMeanQy.at(i+iNumberOfSE)->SetLineWidth(5);
		hMeanQy.at(i+iNumberOfSE)->SetLineColor(1);
		hMeanQy.at(i+iNumberOfSE)->SetMarkerSize(4);
		hMeanQy.at(i+iNumberOfSE)->SetMarkerStyle(20);
		hMeanQy.at(i+iNumberOfSE)->Draw();
	*/
	cCanvas.push_back( new TCanvas("canvas2","Qvectors",4000,2500) ); // 2, Q-vector correlations
	cCanvas.back()->cd();
	int i=1;
	THStack* hStack = new THStack("Stack",";Centrality;Correlations");
	for( auto histo : hCorrelation )
	{
		histo->SetLineWidth(7);
		histo->SetMarkerSize(4);
		histo->SetLineColor(i);
		histo->SetMarkerColor(i);
		histo->SetMarkerStyle(20+i);
		hStack->Add(histo);
		i++;
	}
	hStack->Draw();
	gPad->BuildLegend();
	cout << "Saving Pictures as PNG" << endl;
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
	for(unsigned int i=0; i<iNPSDModules;i++)
	{
		fModule = fEvent->GetPSDModule(i);
		if( fModule->GetId() < 0 )
			continue;
		vModules.push_back(fModule);
	}
	TVector2 fAddition;
	random_shuffle( vModules.begin(),vModules.end() );
	int p=0;
	vector<float> fSumCharge(2);
	for( auto module = begin(vModules); module != end(vModules); ++module )
	{
		float charge = (*module)->GetEnergy();
		float phi = (*module)->GetPhi();
		fAddition.Set( charge*cos(phi), charge*sin(phi) );
		fQvector.at( p%2 )+=fAddition;
		fSumCharge.at(p%2)+=charge;
		p++;
	}
	for(unsigned int i=0; i<2; i++)
	{
		if( fSumCharge.at(i) == 0. )
		{
			fQvector.at(i).Set(-999.,-999.);
			continue;
		}
		fQvector.at(i)*=1/fSumCharge.at(i);
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