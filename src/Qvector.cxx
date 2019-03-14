#include "Qvector.h"

Qvector::Qvector()
{
    for (int i=0;i<2;i++)
    {
        for (int j=0;j<2;j++)
        {
            fQ[i][j]=0;
        }
    }
	this->InitHistograms();
}

void Qvector::InitHistograms()
{
	cout << "Initialization of Q-vectors histograms" << endl;
    vProfile[meanQx1] =         		new TProfile("MeanQx vs Centrality subevent 1",";centrality;mean Qx",10,0,50);
    vProfile[meanQy1] =         		new TProfile("MeanQy vs Centrality subevent 1",";centrality;mean Qy",10,0,50);
	vProfile[meanQx2] =         		new TProfile("MeanQx vs Centrality subevent 2",";centrality;mean Qx",10,0,50);
    vProfile[meanQy2] =         		new TProfile("MeanQy vs Centrality subevent 2",";centrality;mean Qy",10,0,50);
    vProfile[resolution] =     			new TProfile("Resolution vs Centrality",";centrality;R_{1}",10,0,50);

    vHisto1D[QxNotRecentred1] =       	new TH1F("QxNotRecrntred subevent 1",";Qx;counts",100,-1,1);
    vHisto1D[QyNotRecentred1] =       	new TH1F("QyNotRecentred subevent 1",";Qy;counts",100,-1,1);
    vHisto1D[QxRecentred1] =          	new TH1F("QxRecentred subevent 1",";Qx recentred;couts",100,-1,1);
    vHisto1D[QyRecentred1] =          	new TH1F("QyRecrntred subevent 1",";Qy recentred;counts",100,-1,1);
	vHisto1D[PsiEPNotRecentred1]=		new TH1F("PsiEPNotRecentred subevent 1",";#PsiEP [rad];counts",100,-3.15,3.15);
	vHisto1D[PsiEPRecentred1]=			new TH1F("PsiEPNotRecentred subevent 1",";#PsiEP recentred [rad];counts",100,-3.15,3.15);
	
	vHisto1D[QxNotRecentred2] =       	new TH1F("QxNotRecrntred subevent 2",";Qx;counts",100,-1,1);
    vHisto1D[QyNotRecentred2] =       	new TH1F("QyNotRecentred subevent 2",";Qy;counts",100,-1,1);
    vHisto1D[QxRecentred2] =          	new TH1F("QxRecentred subevent 2",";Qx recentred;couts",100,-1,1);
    vHisto1D[QyRecentred2] =          	new TH1F("QyRecrntred subevent 2",";Qy recentred;counts",100,-1,1);
	vHisto1D[PsiEPNotRecentred2]=		new TH1F("PsiEPNotRecentred subevent 2",";#PsiEP [rad];counts",100,-3.15,3.15);
	vHisto1D[PsiEPRecentred2]=			new TH1F("PsiEPNotRecentred subevent 2",";#PsiEP recentred [rad];counts",100,-3.15,3.15);
}

void Qvector::FillCorrections(DataTreeEvent* fEvent)
{
    for (int i=0;i<2;i++)
    {
        for (int j=0;j<2;j++)
        {
            fQ[i][j]=0;
        }
    }
    Int_t iNPSDModules = fEvent->GetNPSDModules();
    DataTreePSDModule* fModule;
    Float_t fPhi;
	Float_t fChargeSum[2];
	Float_t fChargeModule;
	fChargeSum[0]=0; fChargeSum[1]=0;
	vector<DataTreePSDModule*> vModules;
	for(unsigned int i=0; i<iNPSDModules;i++)
	{
		fModule = fEvent->GetPSDModule(i);
		if( fModule->GetId()<0 )
			continue;
		vModules.push_back(fModule);
	}
	random_shuffle( vModules.begin(),vModules.end() );
	for(unsigned int i=0;i<vModules.size();i++)
	{
		fChargeModule = vModules[i]->GetEnergy();
		fPhi = vModules[i]->GetPhi();
		int p = i%2;
		fQ[0][p]+=fChargeModule*cos(fPhi);
		fQ[1][p]+=fChargeModule*sin(fPhi);
		fChargeSum[p]+=fChargeModule;
	}
	for (int i=0;i<2;i++)
	{
		for (int j=0;j<2;j++)
		{
			fQ[i][j]/=fChargeSum[j];
		}
	}
	Float_t fCentrality = fEvent->GetCentrality();
	vProfile[meanQx1]->Fill(fCentrality,fQ[0][0]);
	vProfile[meanQx2]->Fill(fCentrality,fQ[0][1]);
	vProfile[meanQy1]->Fill(fCentrality,fQ[1][0]);
	vProfile[meanQy2]->Fill(fCentrality,fQ[1][1]);

	vHisto1D[QxNotRecentred1]->Fill(fQ[0][0]);
	vHisto1D[QyNotRecentred1]->Fill(fQ[1][0]);
	vHisto1D[QxNotRecentred2]->Fill(fQ[0][1]);
	vHisto1D[QyNotRecentred2]->Fill(fQ[1][1]);
	vHisto1D[PsiEPNotRecentred1]->Fill( TMath::ATan2( fQ[1][0],fQ[0][0]) );
	vHisto1D[PsiEPNotRecentred2]->Fill( TMath::ATan2( fQ[1][1],fQ[0][1]) );
}

void Qvector::Estimate(DataTreeEvent* fEvent)
{
	for (int i=0;i<2;i++)
    {
        for (int j=0;j<2;j++)
        {
            fQ[i][j]=0;
        }
    }
    Int_t iNPSDModules = fEvent->GetNPSDModules();
    DataTreePSDModule* fModule;
    Float_t fPhi;
	Float_t fChargeModule;
	Float_t fChargeSum[2];
	fChargeSum[0]=0; fChargeSum[1]=0;
	vector<DataTreePSDModule*> vModules;
	for(unsigned int i=0; i<iNPSDModules;i++)
	{
		fModule = fEvent->GetPSDModule(i);
		if( fModule->GetId()<0 )
			continue;
		vModules.push_back(fModule);
	}
	random_shuffle( vModules.begin(),vModules.end() );
	for(unsigned int i=0; i<vModules.size() ;i++)
	{
		fChargeModule = vModules[i]->GetEnergy();
		fPhi = vModules[i]->GetPhi();
		int p = i%2;
		fQ[0][p]+=fChargeModule*cos(fPhi);
		fQ[1][p]+=fChargeModule*sin(fPhi);
		fChargeSum[p]+=fChargeModule;
	}
	for (int i=0;i<2;i++)
	{
		for (int j=0;j<2;j++)
		{
			fQ[i][j]/=fChargeSum[j];
		}
	}
	Float_t fCentrality = fEvent->GetCentrality();
	//cout << fCentrality << endl;
	int iCentralityBin = fCentrality/5+1;
	fQ[0][0] -= vProfile[meanQx1]->GetBinContent( iCentralityBin );
	fQ[0][1] -= vProfile[meanQx2]->GetBinContent( iCentralityBin );
	fQ[1][0] -= vProfile[meanQy1]->GetBinContent( iCentralityBin );
	fQ[1][1] -= vProfile[meanQy2]->GetBinContent( iCentralityBin );

	vHisto1D[QxRecentred1]->Fill( fQ[0][0] );
	vHisto1D[QyRecentred1]->Fill( fQ[1][0] );
	vHisto1D[QxRecentred2]->Fill( fQ[0][1] );
	vHisto1D[QyRecentred2]->Fill( fQ[1][1] );
	vHisto1D[PsiEPRecentred1]->Fill( TMath::ATan2( fQ[1][0], fQ[0][0]) );
	vHisto1D[PsiEPRecentred2]->Fill( TMath::ATan2( fQ[1][1], fQ[0][1]) );
	vProfile[resolution]->Fill(fCentrality, cos( TMath::ATan2(fQ[1][0],fQ[0][0]) - TMath::ATan2(fQ[1][1],fQ[0][1]) ) );
}

Float_t Qvector::GetComponent(int i, int j)
{
    return fQ[i][j];
}

Float_t Qvector::GetPsiEP(int j)
{
    return atan2(fQ[1][j],fQ[0][j]);
}

void Qvector::SaveHistograms(TString sPicName)
{
	vCanvas[QvectorsDistribution] = new TCanvas("Canvas1","Qvectors",4500,2000);
	vCanvas[QvectorsDistribution]->Divide(3,2,0.005,0.0001);
	TLegend* legend = new TLegend(0.1,0.8,0.38,0.9);

	vCanvas[QvectorsDistribution]->cd(1);
	vHisto1D[QxNotRecentred1]->SetLineWidth(5);
	vHisto1D[QxRecentred1]->SetLineColor(1);
	vHisto1D[QxRecentred1]->SetLineWidth(5);
	legend->AddEntry(vHisto1D[QxNotRecentred1],"Not Recentred");
	legend->AddEntry(vHisto1D[QxRecentred1],"Recentred");
	vHisto1D[QxNotRecentred1]->Draw();
	vHisto1D[QxRecentred1]->Draw("same");
	legend->Draw();
	
	vCanvas[QvectorsDistribution]->cd(2);
	vHisto1D[QyNotRecentred1]->SetLineWidth(5);
	vHisto1D[QyRecentred1]->SetLineColor(1);
	vHisto1D[QyRecentred1]->SetLineWidth(5);
	vHisto1D[QyNotRecentred1]->Draw();
	vHisto1D[QyRecentred1]->Draw("same");
	legend->Draw();

	vCanvas[QvectorsDistribution]->cd(3);
	vHisto1D[PsiEPNotRecentred1]->SetLineWidth(5);
	vHisto1D[PsiEPRecentred1]->SetLineColor(1);
	vHisto1D[PsiEPRecentred1]->SetLineWidth(5);
	vHisto1D[PsiEPNotRecentred1]->Draw();
	vHisto1D[PsiEPRecentred1]->Draw("same");
	legend->Draw();

	vCanvas[QvectorsDistribution]->cd(4);
	vHisto1D[QxNotRecentred2]->SetLineWidth(5);
	vHisto1D[QxRecentred2]->SetLineColor(1);
	vHisto1D[QxRecentred2]->SetLineWidth(5);
	vHisto1D[QxNotRecentred2]->Draw();
	vHisto1D[QxRecentred2]->Draw("same");
	legend->Draw();
	
	vCanvas[QvectorsDistribution]->cd(5);
	vHisto1D[QyNotRecentred2]->SetLineWidth(5);
	vHisto1D[QyRecentred2]->SetLineColor(1);
	vHisto1D[QyRecentred2]->SetLineWidth(5);
	vHisto1D[QyNotRecentred2]->Draw();
	vHisto1D[QyRecentred2]->Draw("same");
	legend->Draw();

	vCanvas[QvectorsDistribution]->cd(6);
	vHisto1D[PsiEPNotRecentred2]->SetLineWidth(5);
	vHisto1D[PsiEPRecentred2]->SetLineColor(1);
	vHisto1D[PsiEPRecentred2]->SetLineWidth(5);
	vHisto1D[PsiEPNotRecentred2]->Draw();
	vHisto1D[PsiEPRecentred2]->Draw("same");
	legend->Draw();

	vCanvas[QvectorsDistribution]->SaveAs("../histograms/"+sPicName+".png");

	vCanvas[MeanQvectors] = new TCanvas("Canvas2","Qvectors",4000,2500);
	vCanvas[MeanQvectors]->Divide(3,2,0.005,0.0001);

	vCanvas[MeanQvectors]->cd(1);
	vProfile[meanQx1]->Draw();

	vCanvas[MeanQvectors]->cd(2);
	vProfile[meanQx2]->Draw();

	vCanvas[MeanQvectors]->cd(3);
	vProfile[meanQy1]->Draw();

	vCanvas[MeanQvectors]->cd(3);
	vProfile[meanQy2]->Draw();

	vCanvas[MeanQvectors]->SaveAs("../histograms/"+sPicName+"_1.png");
}