#pragma once
#include <iostream>

#include "TROOT.h"
#include "TFile.h"
#include "TStyle.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TProfile.h"
#include "TCanvas.h"
#include "TLegend.h"

#include "Selector.h"
#include "Centrality.h"

#include "HADES_constants.h"

using std::vector;
using std::cout;
using std::endl;

class EventQA
{
    private:
    enum eQAHisto1DMap{
        tracksMDC = 0,
        tracksMDC_selected,
        hitsTOF,
        hitsTOF_selected,
        chargeFW,
        chargeFW_selected,
        vertexZ,
        vertexZ_selected,
        hitsTOF_uncuted,
        hitsTOF_uncuted_selected,
        hitsTOF_matched,
        hitsTOF_matched_selected,
		histo_centrality,
		histo_centrality_selected,
        Num1DHistos
    };
    enum eQAHisto2DMap{
        tracks_hits = 0,
        tracks_hits_selected,
        tracks_charge,
        tracks_charge_selected,
        hits_charge,
        hits_charge_selected,
        vertexX_vertexY,
        vertexX_vertexY_selected,
        hitsFW_X_Y,
        hitsFW_X_Y_selected,
        Num2DHistos
    };
	enum eQAProfileMap{
		hits_centrality = 0,
		hits_centrality_selected,
		NumProfiles
	};
    enum eCanvasMap{
        multiplicity_vertex = 0,
        multiplicity_charge,
        vertex_charge,
		centrality,
        NumCanvases
    };
    TH1F* vHisto1D[Num1DHistos];
    TH2F* vHisto2D[Num2DHistos];
	TProfile* vProfile[NumProfiles];
    TCanvas* vCanvas[NumCanvases];
    Selector* fSelector;
	Centrality* fCentrality;
	DataTreeEvent* fEvent;
    EventQA();
    public:
	EventQA(DataTreeEvent* _fEvent, Selector* _selector, Centrality* _centrality);
    ~EventQA();
    void    InitHistograms();
	void	LoadCentrality(Centrality* _centrality) { fCentrality = _centrality; }
    void    FillHistograms();
    void    SaveHistograms(TString PicName);
    void    SaveHisogramsToROOTFile(TString FileName);
};