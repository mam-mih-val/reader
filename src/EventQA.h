#pragma once
#include <iostream>

#include "TROOT.h"
#include "TFile.h"
#include "TStyle.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TLegend.h"

#include "Selector.h"

#include "HADES_constants.h"

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
    enum eCanvasMap{
        multiplicity_vertex = 0,
        multiplicity_charge,
        vertex_charge,
        NumCanvases
    };
    TH1F* vHisto1D[Num1DHistos];
    TH2F* vHisto2D[Num2DHistos];
    TCanvas* vCanvas[NumCanvases];
    Selector fSelector;
    public:
    EventQA();
    ~EventQA();
    void    InitHistograms();
    void    FillHistograms(DataTreeEvent* fEvent);
    void    SaveHistograms(TString PicName);
    void    SaveHisogramsToROOTFile(TString FileName);
};