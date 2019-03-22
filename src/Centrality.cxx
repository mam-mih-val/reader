#include "Centrality.h"

Centrality::Centrality(TString FileName)
{
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
	}
	else
		cout << "Couldn't open the file" << endl;
}

float Centrality::GetCentralityClass(DataTreeEvent* fEvent)
{
	auto TOFRPChits = fEvent->GetCentralityEstimator(HADES_constants::kNhitsTOF_cut) + fEvent->GetCentralityEstimator(HADES_constants::kNhitsRPC_cut);
	auto bin = hCentralityPercentile->FindBin(TOFRPChits);
	return hCentralityPercentile->GetBinContent(bin);
}