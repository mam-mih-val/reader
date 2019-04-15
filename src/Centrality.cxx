#include "Centrality.h"

Centrality::Centrality(DataTreeEvent* _fEvent, TString FileName)
{
	fEvent = _fEvent;
	this->LoadCentralityPercentile(FileName);
}

void Centrality::LoadCentralityPercentile(TString FileName)
{
	auto file = new TFile(FileName);
	if ( file->IsOpen() )
	{
		cout << "Centrality file loaded" << endl;
		file->cd("/Centrality/");
		hCentralityPercentile = (TH1F*) file->Get("/Centrality/TOFRPCtot_5pc_fixedCuts");
		cout << hCentralityPercentile->GetNbinsX() << " centrality classes" << endl;
		//file->Close();
	}
	else
		cout << "Couldn't open the file" << endl;
}

float Centrality::GetCentralityClass()
{
	auto TOFRPChits = fEvent->GetCentralityEstimator(HADES_constants::kNhitsTOF_RPC_cut);
	auto bin = hCentralityPercentile->FindBin(TOFRPChits);
	return hCentralityPercentile->GetBinContent(bin);
}

int	Centrality::GetNumClasses() 
{ 
	return hCentralityPercentile->GetNbinsX()-1; 
}
