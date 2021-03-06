#include "Flow.h"

void Flow::InitializeQaHistograms()
{
    h1dQa.at(kCentrality) = new TH1F( "centrality", ";centrality;entries", 20, 0.0, 100.0 );
    h1dQa.at(kCentralityEstimator) = new TH1F( "centrality estimator", ";TOF+RPC hits;entries", 250, 0.0, 250.0 );
    hProfileQa.at(kCosPhiYcm) = new TProfile( "cos(#phi) vs y_{cm}", ";y_{cm}; <cos(#phi)>", 17, -0.85, 0.85 );
    hProfileQa.at(kSinPhiYcm) = new TProfile( "sin(#phi) vs y_{cm}", ";y_{cm}; <sin(#phi)>", 17, -0.85, 0.85 );
    
    h2dQa.at(kM2VsP) = new TH2F( "m^{2} vs p", ";p, #frac{GeV}{c};m^{2}, #frac{GeV^{2}}{c^{4}}", 100, -3.0, 3.0, 100, 0, 18.0 );
    h2dQa.at(kPtVsYforward) = new TH2F( "pt vs y_{cm} forward", ";y_{cm}^{fwrd};pt, #frac{GeV}{c}", 100, 0, 1.0, 300, 0, 3.0 );
    h2dQa.at(kPtVsYbackward) = new TH2F( "pt vs y_{cm} backward", ";y_{cm}^{bwrd};pt, #frac{GeV}{c}", 100, 0, 1.0, 300, 0, 3.0 );
    h2dQa.at(kFwAdcVsEstimator) = new TH2F( "FW-ADC vs TOF+RPC-Hits", ";TOF+RPC-Hits;FW-ADC", 250, 0., 250.0, 100, 0.0, 8000.0 );
    h2dQa.at(kFwZVsEstimator) = new TH2F( "FW-Z vs TOF+RPC-Hits", ";TOF+RPC-Hits;FW-Z", 250, 0., 250.0, 100, 0.0, 100.0 );
    h2dQa.at(kFwAdcVsModuleId) = new TH2F( "FW-Module Id vs FW-ADC", ";FW-ADC;FW-Module Id", 100, 0.0, 1000.0, 304, 0., 304.0 );
    h2dQa.at(kFwZVsModuleId) = new TH2F( "FW-Module Id vs FW-Z", ";FW-Z;FW-Module Id", 20, 0.0, 20.0, 304, 0., 304.0 );
    h2dQa.at(kFwTimeVsModuleId) = new TH2F( "FW-Module Id vs FW-Time", ";FW-Time;FW-Module Id", 200, 0.03, 0.05, 304, 0., 304.0 );
    h2dQa.at(kFwBetaVsModuleId) = new TH2F( "FW-Module Id vs FW-beta", ";FW-beta;FW-Module Id", 90, 0.7, 1.2, 304, 0., 304.0 );
}

void Flow::FillQaHistograms(bool channelSelection, TString signal, float minSignal, float maxSignal)
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
        if( fEvent->GetVertexTrack(i)->GetPdgId() != fPid )
            continue;
        auto momentum = fEvent->GetVertexTrack(i)->GetMomentum();
        auto tof = fEvent->GetTOFHit(i)->GetTime();
		auto charge = fEvent->GetTOFHit(i)->GetCharge(); 
		auto len = fEvent->GetTOFHit(i)->GetPathLength();
		auto beta = len/tof/299.792458;
		auto mass2 = momentum.P()*momentum.P()*(1-beta*beta)/(beta*beta);
        h2dQa.at(kM2VsP)->Fill(momentum.P(), mass2);
        momentum.Boost(boost);
        if( momentum.Rapidity() >= 0 )
            h2dQa.at(kPtVsYforward)->Fill( momentum.Rapidity(), momentum.Pt() );
        if( momentum.Rapidity() < 0 )
            h2dQa.at(kPtVsYbackward)->Fill( fabs(momentum.Rapidity()), momentum.Pt());
        hProfileQa.at(kCosPhiYcm)->Fill( momentum.Rapidity(), cos( momentum.Phi() ) );
        hProfileQa.at(kSinPhiYcm)->Fill( momentum.Rapidity(), sin( momentum.Phi() ) );
    }

    int nFwModules = fEvent->GetNPSDModules();
    int sumAdc=0, sumZ=0;
    for(int i=0; i<nFwModules; i++)
    {
        if( !fSelector->IsCorrectFwHit(i, channelSelection, signal, minSignal, maxSignal) )
            continue;
        h2dQa.at(kFwAdcVsModuleId)->Fill( fEvent->GetPSDModule(i)->GetEnergy(), fEvent->GetPSDModule(i)->GetId() );
        sumAdc+=fEvent->GetPSDModule(i)->GetEnergy();
        h2dQa.at(kFwZVsModuleId)->Fill( fEvent->GetPSDModule(i)->GetChargeZ(), fEvent->GetPSDModule(i)->GetId() );
        sumZ+=fEvent->GetPSDModule(i)->GetChargeZ();
        float beta = fEvent->GetPSDModule(i)->GetBeta();
        float x = fEvent->GetPSDModule(i)->GetPositionComponent(0);
        float y = fEvent->GetPSDModule(i)->GetPositionComponent(1);
        float z = fEvent->GetPSDModule(i)->GetPositionComponent(2);
        float distance = sqrt(x*x+y*y+z*z);
        float time = 299.792458*beta/distance;
        h2dQa.at(kFwBetaVsModuleId)->Fill( fEvent->GetPSDModule(i)->GetBeta(), fEvent->GetPSDModule(i)->GetId() );
        h2dQa.at(kFwTimeVsModuleId)->Fill( time, fEvent->GetPSDModule(i)->GetId() );
    }
    h2dQa.at(kFwAdcVsEstimator)->Fill(fEvent->GetCentralityEstimator(HADES_constants::kNhitsTOF_RPC_cut),sumAdc);
    h2dQa.at(kFwZVsEstimator)->Fill(fEvent->GetCentralityEstimator(HADES_constants::kNhitsTOF_RPC_cut),sumZ);
}

void Flow::InitializeMeanUn()
{
    hMeanUn.at(0) = new TProfile2D( "Mean_Unx", "<x_{n}>;y_{cm};pt, [GeV/c]", 16, -0.8, 0.8, 25, 0.0, 2.5 );
    hMeanUn.at(1) = new TProfile2D( "Mean_Uny", "<y_{n}>;y_{cm};pt, [GeV/c]", 16, -0.8, 0.8, 25, 0.0, 2.5 );
}


void Flow::FillMeanUn()
{
	int nTracks = fEvent->GetNVertexTracks();
    double BETA = sqrt( 1.0 - 0.938*0.938/1.23/1.23 );
    TVector3 boost{0,0,-BETA};
    for(int i=0; i<nTracks; i++)
    {
        auto track = fEvent->GetVertexTrack(i);
        if( !fSelector->IsCorrectTrack(i) )
            continue;
        auto p = track->GetMomentum();
        p.Boost(boost);
        TVector2 un;
        un.SetMagPhi(1.0, p.Phi());
        hMeanUn.at(0)->Fill(p.Rapidity(), p.Pt(), un.X() ); 
        hMeanUn.at(1)->Fill(p.Rapidity(), p.Pt(), un.Y() ); 
    }
}