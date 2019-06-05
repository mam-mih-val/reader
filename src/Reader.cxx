#include "Reader.h"

const double YCOR = 0.5*log(1.23*197+156.743) - 0.5*log(1.23*197-156.743);

Reader::Reader(TString cFileName)
{
	cout << "Reader initialization" << endl;
    fChain = new TChain("DataTree");
    fChain->Add(cFileName);
    cout << fChain->GetEntries() << " events found" << endl;
	fEvent = nullptr;
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

void Reader::BuildEventQaHistograms(TString sPicName)
{
	Selector* fSelector = new Selector(fEvent);
	Centrality* fCentrality = new Centrality(fEvent, "centrality_epcorr_apr12_gen8_2018_07.root");
    EventQA* fEventQA = new EventQA(fEvent, fSelector, fCentrality);
    Long64_t lNEvents = fChain->GetEntries();
    for(int i=0; i<lNEvents; i++)
    {
        fChain->GetEntry(i);
		if( !fEvent->GetTrigger(HADES_constants::kPT2)->GetIsFired() )
			continue;
        fEventQA->FillHistograms();
    }
    fEventQA->SaveHisogramsToROOTFile(sPicName);
}

void Reader::BuildTrackQaHistograms(TString sPicName)
{
	Selector* fSelector = new Selector(fEvent);
	Centrality* fCentrality = new Centrality(fEvent, "centrality_epcorr_apr12_gen8_2018_07.root");
	TrackQA* fTrackQA = new TrackQA(fEvent,fSelector,-1);
    Long64_t lNEvents = fChain->GetEntries();
    for(int i=0; i<lNEvents; i++)
    {
        fChain->GetEntry(i);
		fTrackQA->FillHistograms();
    }
	fTrackQA->SaveHisogramsToROOTFile(sPicName);
}
void Reader::BuildQvector3SeHistograms(TString sPicName, bool channelSelection, TString signal, float minSignal, float maxSignal)
{
	Long64_t lNEvents = fChain->GetEntries();
    Selector* fSelector = new Selector(fEvent);
	Centrality* fCentrality = new Centrality(fEvent,"centrality_epcorr_apr12_gen8_2018_07.root");
	Qvector3SE* fQ =  new Qvector3SE(fEvent, fSelector, fCentrality, channelSelection, signal, minSignal, maxSignal);
	cout << "Filling correction histograms" << endl;
	for(int i=0; i<lNEvents; i++)
    {
        fChain->GetEntry(i);
		if( !fSelector->IsCorrectEvent() )
			continue;
		fQ->ComputeCorrections();
    }
	cout << "Estimating Q-vectors" << endl;
	for(int i=0; i<lNEvents; i++)
    {
        fChain->GetEntry(i);
		if( !fSelector->IsCorrectEvent() )
			continue;
		fQ->ComputeCorrelations();
    }
	cout << "Estimating resolution" << endl;
	fQ->ComputeResolution();
	fQ->SavePictures(sPicName);
	fQ->SaveHistogramsToROOTFile(sPicName);
}

void Reader::BuildFlow3SeHistograms(TString sPicName, bool channelSelection, TString signal, float minSignal, float maxSignal)
{
	Long64_t lNEvents = fChain->GetEntries();
    Selector* fSelector = new Selector(fEvent);
	Centrality* fCentrality = new Centrality(fEvent,"centrality_epcorr_apr12_gen8_2018_07.root");
	Qvector3SE* fQ =  new Qvector3SE(fEvent, fSelector, fCentrality, channelSelection, signal, minSignal, maxSignal);
	cout << "Filling correction histograms" << endl;
	for(int i=0; i<lNEvents; i++)
    {
        fChain->GetEntry(i);
		if( !fSelector->IsCorrectEvent() )
			continue;
		if( fEvent->GetCentrality(HADES_constants::kNhitsTOF_RPC_cut) > 40. )
			continue;
		fQ->ComputeCorrections();
    }
	cout << "Estimating Q-vecor correlations" << endl;
	for(int i=0; i<lNEvents; i++)
    {
        fChain->GetEntry(i);
		if( !fSelector->IsCorrectEvent() )
			continue;
		fQ->ComputeCorrelations();
    }
	Flow3SE* flow = new Flow3SE(fEvent, fCentrality, fQ, fSelector, 14);
	cout << "Estimating flow" << endl;
	for(int i=0; i<lNEvents; i++)
    {
		fChain->GetEntry(i);
		if( fEvent->GetCentrality(HADES_constants::kNhitsTOF_RPC_cut) < 20 || fEvent->GetCentrality(HADES_constants::kNhitsTOF_RPC_cut) > 30 )
			continue;
		if( !fSelector->IsCorrectEvent() )
			continue;
		flow->Estimate();
    }
	// flow->SavePictures(sPicName);
	flow->SaveHistogramsToRootFile(sPicName);
}