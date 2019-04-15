#include "Qvector.h"

Qvector::Qvector(DataTreeEvent* _fEvent, Centrality* _centrality, unsigned int NumSE)
{
	fEvent = _fEvent;
	fCentrality = _centrality;
	iNumberOfSE = NumSE;
	for(unsigned int i=0;i<iNumberOfSE;i++)
		fQvector.push_back( TVector2(0.,0.) );
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

void Qvector::Recenter()
{
	int bin = hMeanQx.at(0)->FindBin( fEvent->GetCentrality(HADES_constants::kNhitsTOF_RPC_cut) );
	for(unsigned int i=0;i<fQvector.size();i++)
	{
		if( fQvector.at(i).X() < -990. )
			continue;
		TVector2 fCorrVec;
		fCorrVec.Set( hMeanQx.at(i)->GetBinContent(bin), hMeanQy.at(i)->GetBinContent(bin) );
		fQvector.at(i)-=fCorrVec;
	}
}