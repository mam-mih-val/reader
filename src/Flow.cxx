#include "Flow.h"

void Flow::InitializeQaHistograms()
{
    h1dQa.at(kCentrality) = new TH1F( "centrality", ";centrality;entries", 20, 0.0, 100.0 );
    h1dQa.at(kCentralityEstimator) = new TH1F( "centrality estimator", ";TOF+RPC hits;entries", 250, 0.0, 250.0 );
    
    h2dQa.at(kM2VsP) = new TH2F( "m^{2} vs p", ";p, #frac{GeV}{c};m^{2}, #frac{GeV^{2}}{c^{4}}", 100, -3.0, 3.0, 100, 0, 18.0 );
    h2dQa.at(kPtVsY) = new TH2F( "pt vs y_{cm}", ";y_{cm};pt, #frac{GeV}{c}", 100, -1.0, 1.0, 100, 0, 3.0 );
    h2dQa.at(kFwAdcVsEstimator) = new TH2F( "FW-ADC vs TOF+RPC-Hits", ";TOF+RPC-Hits;FW-ADC", 250, 0., 250.0, 100, 0.0, 8000.0 );
    h2dQa.at(kFwZVsEstimator) = new TH2F( "FW-Z vs TOF+RPC-Hits", ";TOF+RPC-Hits;FW-Z", 250, 0., 250.0, 100, 0.0, 100.0 );
    h2dQa.at(kFwAdcVsModuleId) = new TH2F( "FW-Module Id vs FW-ADC", ";FW-ADC;FW-Module Id", 100, 0.0, 1000.0, 304, 0., 304.0 );
    h2dQa.at(kFwZVsModuleId) = new TH2F( "FW-Module Id vs FW-Z", ";FW-Z;FW-Module Id", 20, 0.0, 20.0, 304, 0., 304.0 );
}

void Flow::FillQaHistograms(bool channelSelection, bool protonSpectators)
{
    h1dQa.at(kCentrality)->Fill( fEvent->GetCentrality(HADES_constants::kNhitsTOF_RPC_cut) );
    h1dQa.at(kCentralityEstimator)->Fill( fEvent->GetCentralityEstimator(HADES_constants::kNhitsTOF_RPC_cut) );

    double BETA = sqrt( 1.0 - 0.938*0.938/1.23/1.23 );
    TVector3 boost{0,0,-BETA};
    auto nTracks = fEvent->GetNVertexTracks();
    for(int i=0; i<nTracks; i++)
    {
        if( !fSelector->IsCorrectTrack(i) )
            continue;
        if( !fEvent->GetVertexTrack(i)->GetPdgId() )
            continue;
        auto momentum = fEvent->GetVertexTrack(i)->GetMomentum();
        auto tof = fEvent->GetTOFHit(i)->GetTime();
		auto charge = fEvent->GetTOFHit(i)->GetCharge(); 
		auto len = fEvent->GetTOFHit(i)->GetPathLength();
		auto beta = len/tof/299.792458;
		auto mass2 = momentum.P()*momentum.P()*(1-beta*beta)/(beta*beta);
        h2dQa.at(kM2VsP)->Fill(momentum.P(), mass2);
        momentum.Boost(boost);
        h2dQa.at(kPtVsY)->Fill(momentum.Rapidity(),momentum.Pt());
    }

    int nFwModules = fEvent->GetNPSDModules();
    int sumAdc=0, sumZ=0;
    for(int i=0; i<nFwModules; i++)
    {
        if( !fSelector->IsCorrectFwHit(i, channelSelection, protonSpectators) )
            continue;
        h2dQa.at(kFwAdcVsModuleId)->Fill( fEvent->GetPSDModule(i)->GetEnergy(), fEvent->GetPSDModule(i)->GetId() );
        h2dQa.at(kFwZVsModuleId)->Fill( fEvent->GetPSDModule(i)->GetChargeZ(), fEvent->GetPSDModule(i)->GetId() );
        sumAdc+=fEvent->GetPSDModule(i)->GetEnergy();
        sumZ+=fEvent->GetPSDModule(i)->GetChargeZ();
    }
    h2dQa.at(kFwAdcVsEstimator)->Fill(fEvent->GetCentralityEstimator(HADES_constants::kNhitsTOF_RPC_cut),sumAdc);
    h2dQa.at(kFwZVsEstimator)->Fill(fEvent->GetCentralityEstimator(HADES_constants::kNhitsTOF_RPC_cut),sumZ);
}