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
    fEventQA->SaveHistograms(sPicName);
}

void Reader::BuildTrackQaHistograms(TString sPicName)
{
	Selector* fSelector = new Selector(fEvent);
	Centrality* fCentrality = new Centrality(fEvent, "centrality_epcorr_apr12_gen8_2018_07.root");
	TrackQA* fTrackQA[NumOfParticles];
	fTrackQA[all] = 		new TrackQA(fEvent,fSelector,-1);
	fTrackQA[electron] = 	new TrackQA(fEvent,fSelector,3);
	fTrackQA[positron] = 	new TrackQA(fEvent,fSelector,2);
	fTrackQA[pi_minus] = 	new TrackQA(fEvent,fSelector,9);
	fTrackQA[pi_plus] = 	new TrackQA(fEvent,fSelector,8);
	fTrackQA[proton] = 		new TrackQA(fEvent,fSelector,14);
	fTrackQA[deuteron] =  	new TrackQA(fEvent,fSelector,45);
	fTrackQA[helium3] = 	new TrackQA(fEvent,fSelector,49);
	fTrackQA[helium4] = 	new TrackQA(fEvent,fSelector,47);
    Long64_t lNEvents = fChain->GetEntries();
    for(int i=0; i<lNEvents; i++)
    {
        fChain->GetEntry(i);
		for(int j=0; j<NumOfParticles;j++) 
			fTrackQA[j]->FillHistograms();
    }
	for(int j=0; j<NumOfParticles;j++) 
		fTrackQA[j]->SaveHistograms(sPicName);
}
void Reader::BuildQvector3SeHistograms(TString sPicName, bool channelSelection, bool fwZ, bool protonSpectators)
{
	Long64_t lNEvents = fChain->GetEntries();
    Selector* fSelector = new Selector(fEvent);
	Centrality* fCentrality = new Centrality(fEvent,"centrality_epcorr_apr12_gen8_2018_07.root");
	Qvector3SE* fQ =  new Qvector3SE(fEvent, fSelector, fCentrality, channelSelection, fwZ, protonSpectators);
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

void Reader::BuildFlow3SeHistograms(TString sPicName, bool channelSelection, bool fwZ, bool protonSpectators)
{
	Long64_t lNEvents = fChain->GetEntries();
    Selector* fSelector = new Selector(fEvent);
	Centrality* fCentrality = new Centrality(fEvent,"centrality_epcorr_apr12_gen8_2018_07.root");
	Qvector3SE* fQ =  new Qvector3SE(fEvent, fSelector, fCentrality, channelSelection, fwZ, protonSpectators);
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