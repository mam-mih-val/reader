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

using std::vector;
using std::cout;
using std::endl;

class TrackQA
{
    private:
    int iPid = -1; // Pid = -1 - all particle types
    enum eQAHisto1DMap{
        ptMDC = 0,
        ptMDC_selected, // 1
        massTOF, // 2
        massTOF_selected, // 3
        rapidityMDC, // 4
        rapidityMDC_selected, // 5
        phiMDC, // 6
        phiMDC_selected, // 7
        betaTOF, // 8
        betaTOF_selected, // 9
        pseudorapidityMDC, // 10
        pseudorapidityMDC_selected, // 11
        Num1DHistos // 12
    };
    enum eQAHisto2DMap{
        phi_rapidity = 0,
        phi_rapidity_selected, // 1
        phi_pt, // 2
        phi_pt_selected, // 3
        pt_rapidity, // 4
        pt_rapidity_selected, // 5
        phi_pseudorapidity, // 6 
        phi_pseudorapidity_selected, // 7 
        pt_pseudorapidity, // 8
        pt_pseudorapidity_selected, // 9
        rapidity_pseudorapidity, // 10
        rapidity_pseudorapidity_selected, // 11
        dEdXTOF_p, // 12
        dEdXTOF_p_selected, // 13
        dEdXMDC_p, // 14
        dEdXMDC_p_selected, // 15
        beta_p, // 16
        beta_p_selected, // 17
		DCA_X_Y, // 18
		DCA_X_Y_selected, // 19
        Num2DHistos // 20
    };
    enum eCanvasMap{
        OneDimHist = 0, // 
        phi_kinematics, // 1, phi vs pt, y, eta
        pt_kinematics, // 2, pt vs y, eta & y vs eta
        mass_qa, // 3
		DCA, //4
        NumCanvases // 5
    };
    TH1F* vHisto1D[Num1DHistos];
    TH2F* vHisto2D[Num2DHistos];
    TCanvas* vCanvas[NumCanvases];
    Selector fSelector;
    public:
    TrackQA() {};
	TrackQA(int iPid);
    ~TrackQA();
    void    InitHistograms();
    void    SetPid(int _Pid=-1) { iPid = _Pid; }
    void    FillHistograms(DataTreeEvent* fEvent);
    void    SaveHistograms(TString PicName);
    void    SaveHisogramsToROOTFile(TString FileName);
};