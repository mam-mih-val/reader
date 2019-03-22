#include "Reader.h"

const double YCOR = 0.5*log(1.23*197+156.743) - 0.5*log(1.23*197-156.743);

Reader::Reader(char* cFileName)
{
    fChain = new TChain("DataTree");
    fChain->Add(cFileName);
    cout << fChain->GetEntries() << "events" << endl;
    fEvent = new DataTreeEvent;
    fChain->SetBranchAddress("DTEvent", &fEvent);
}

Reader::~Reader()
{
    delete fChain;
    delete fEvent;
}

void Reader::AddFile(char* cFileName)
{
    fChain->Add(cFileName);
    cout << fChain->GetEntries() << " " << "events" << endl;
}

DataTreeEvent* Reader::GetEvent(int idx)
{
    fChain->GetEntry(idx);
    return fEvent;
}

void Reader::BuildQAHistograms(TString sPicName)
{
    EventQA* fEventQA = new EventQA;
	TrackQA* fTrackQA[NumOfParticles];
	fTrackQA[all] = 		new TrackQA(-1);
	fTrackQA[electron] = 	new TrackQA(3);
	fTrackQA[positron] = 	new TrackQA(2);
	fTrackQA[pi_minus] = 	new TrackQA(9);
	fTrackQA[pi_plus] = 	new TrackQA(8);
	fTrackQA[proton] = 		new TrackQA(14);
	fTrackQA[deuteron] =  	new TrackQA(45);
	fTrackQA[helium3] = 	new TrackQA(49);
	fTrackQA[helium4] = 	new TrackQA(47);
    Long64_t lNEvents = fChain->GetEntries();
    for(int i=0; i<lNEvents; i++)
    {
        fChain->GetEntry(i);
        fEventQA->FillHistograms(fEvent);
		for(int j=0; j<NumOfParticles;j++) 
			fTrackQA[j]->FillHistograms(fEvent);
    }
    fEventQA->SaveHistograms(sPicName);
	for(int j=0; j<NumOfParticles;j++) 
		fTrackQA[j]->SaveHistograms(sPicName);
}

void Reader::BuildQvectorHistograms(TString sPicName)
{
	Long64_t lNEvents = fChain->GetEntries();
    Selector* fSelector = new Selector;
	Centrality* fCentrality = new Centrality("centrality_epcorr_apr12_gen8_2018_07.root");
	Qvector* fQ =  new Qvector(fCentrality,2);
	//fQ->LoadCentrality(fCentrality);
	cout << "Filling correction histograms" << endl;
	for(int i=0; i<lNEvents; i++)
    {
        fChain->GetEntry(i);
		if( !fSelector->IsCorrectEvent(fEvent) )
			continue;
		fQ->FillCorrections(fEvent);
    }
	cout << "Estimating Q-vectors" << endl;
	for(int i=0; i<lNEvents; i++)
    {
        fChain->GetEntry(i);
		if( !fSelector->IsCorrectEvent(fEvent) )
			continue;
		fQ->Estimate(fEvent);
    }
	fQ->SaveHistogramsToROOTFile(sPicName);
	fQ->SaveHistograms(sPicName);
}