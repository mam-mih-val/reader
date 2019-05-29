vector<TGraphErrors*> YcmSystematics(TFile* file);
vector<TGraphErrors*> PtSystematics(TFile* file);

void ErrorPainter(TString fileName="000")
{
    gROOT->ProcessLine(".x SetStyle.C");
    auto file = new TFile(fileName+".root");
    vector<TGraphErrors*> hYcmSys = YcmSystematics(file);
    vector<TGraphErrors*> hPtSys = PtSystematics(file);
    
    vector<TMultiGraph*> stack{ 
    new TMultiGraph("Geometry ycm", "Geometry;y_{cm};Relative Error, [%]"),                 // 0
    new TMultiGraph("Sub-Events ycm", "Sub-Events;y_{cm};Relative Error, [%]"),             // 1
    new TMultiGraph("Geometry pt", "Geometry;pt, [#frac{Gev}{c}];Relative Error, [%]"),     // 2 
    new TMultiGraph("Sub-Events pt", "Sub-Events;pt, [#frac{Gev}{c}];Relative Error, [%]")  // 3
    };
    for(int i=0; i<2; i++)
    {
        hYcmSys.at(i)->SetMarkerStyle(20+i);
        hYcmSys.at(i)->SetMarkerSize(2);
        hYcmSys.at(i)->SetMarkerColor(i+1);
        hYcmSys.at(i)->SetLineColor(i+1);

        hPtSys.at(i)->SetMarkerStyle(20+i);
        hPtSys.at(i)->SetMarkerSize(2);
        hPtSys.at(i)->SetMarkerColor(i+1);
        hPtSys.at(i)->SetLineColor(i+1);
        
        stack.at(0)->Add(hYcmSys.at(i));
        stack.at(2)->Add(hPtSys.at(i));
    }
    for(int i=2; i<hYcmSys.size(); i++)
    {
        hYcmSys.at(i)->SetMarkerStyle(18+i);
        hYcmSys.at(i)->SetMarkerSize(2);
        hYcmSys.at(i)->SetMarkerColor(i-1);
        hYcmSys.at(i)->SetLineColor(i-1);

        hPtSys.at(i)->SetMarkerStyle(18+i);
        hPtSys.at(i)->SetMarkerSize(2);
        hPtSys.at(i)->SetMarkerColor(i-1);
        hPtSys.at(i)->SetLineColor(i-1);

        stack.at(1)->Add(hYcmSys.at(i));
        stack.at(3)->Add(hPtSys.at(i));
    }
    vector<TCanvas*> canvas{ new TCanvas("canv1", "", 2000, 800), new TCanvas("canv2", "", 2000, 800), };
    int bias = 0;
    for(int i=0; i<canvas.size(); i++)
    {
        canvas.at(i)->Divide(2,1);
        for(int j=0; j<2; j++)
        {
            canvas.at(i)->cd(j+1);
            stack.at(j+bias)->GetYaxis()->SetRangeUser(-5., 15.);
            stack.at(j+bias)->Draw("A");
            gPad->BuildLegend();
        }
        canvas.at(i)->SaveAs(fileName+Form("_%i_systematics.png", i) );
        bias+=2;
    }
}

vector<TGraphErrors*> YcmSystematics(TFile* file)
{

    enum eProfiles{
        kUnited=0,
        kX,
        kY,
        kA,
        kB,
        kC,
        kNumberOfProfiles
    };
    array<TProfile*, kNumberOfProfiles> hYcm;
    hYcm.at(kUnited) = new  TProfile( "v1_vs_y", "v_{1} vs y; y_{cm}; v_{1}", 14, -0.7, 0.7 );
    hYcm.at(kX) = new TProfile( "v1_vs_y_X", "v_{1} vs y; y_{cm}; v_{1}^{x}", 14, -0.7, 0.7 );
    hYcm.at(kY) = new TProfile( "v1_vs_y_Y", "v_{1} vs y; y_{cm}; v_{1}^{y}", 14, -0.7, 0.7 );
    hYcm.at(kA) = new TProfile( "v1_vs_y_A", "v_{1} vs y; y_{cm}; v_{1}^{a}", 14, -0.7, 0.7 );
    hYcm.at(kB) = new TProfile( "v1_vs_y_B", "v_{1} vs y; y_{cm}; v_{1}^{b}", 14, -0.7, 0.7 );
    hYcm.at(kC) = new TProfile( "v1_vs_y_C", "v_{1} vs y; y_{cm}; v_{1}^{c}", 14, -0.7, 0.7 );
    for( auto &histo : hYcm )
        histo->Sumw2();
    for(int se=0; se<3; se++)
    {
        int wX = ( (TProfile*) file->Get( Form( "/Flow/0.8 < pt < 0.85 SE%i_{x}", se) ) )->GetEntries();
        int wY = ( (TProfile*) file->Get( Form( "/Flow/0.8 < pt < 0.85 SE%i_{x}", se) ) )->GetEntries();
        hYcm.at(kUnited)->Add( (TProfile*) file->Get( Form( "/Flow/0.8 < pt < 0.85 SE%i_{x}", se) ), wX );
        hYcm.at(kUnited)->Add( (TProfile*) file->Get( Form( "/Flow/0.8 < pt < 0.85 SE%i_{y}", se) ), wY );
        hYcm.at(kX)->Add( (TProfile*) file->Get( Form( "/Flow/0.8 < pt < 0.85 SE%i_{x}", se) ), wX );
        hYcm.at(kY)->Add( (TProfile*) file->Get( Form( "/Flow/0.8 < pt < 0.85 SE%i_{y}", se) ), wY );
        hYcm.at(kA+se)->Add( (TProfile*) file->Get( Form( "/Flow/0.8 < pt < 0.85 SE%i_{x}", se) ), wX );
        hYcm.at(kA+se)->Add( (TProfile*) file->Get( Form( "/Flow/0.8 < pt < 0.85 SE%i_{y}", se) ), wY );
    }
    hYcm.at(kA)->Add( (TProfile*) file->Get( "/Flow/0.8 < pt < 0.85 SE0_{x}" ) );
    hYcm.at(kA)->Add( (TProfile*) file->Get( "/Flow/0.8 < pt < 0.85 SE0_{y}" ) );
    hYcm.at(kB)->Add( (TProfile*) file->Get( "/Flow/0.8 < pt < 0.85 SE1_{x}" ) );
    hYcm.at(kB)->Add( (TProfile*) file->Get( "/Flow/0.8 < pt < 0.85 SE1_{y}" ) );
    hYcm.at(kC)->Add( (TProfile*) file->Get( "/Flow/0.8 < pt < 0.85 SE2_{x}" ) );
    hYcm.at(kC)->Add( (TProfile*) file->Get( "/Flow/0.8 < pt < 0.85 SE2_{y}" ) );
    vector<TGraphErrors*> hYcmSys;
    int nBinsYcm = hYcm.at(kUnited)->GetNbinsX();
    for(int i=0; i<kNumberOfProfiles-1; i++)
    {
        hYcmSys.push_back( new TGraphErrors(nBinsYcm) );
    }
    hYcmSys.at(kX-1)->SetTitle("x-Component");
    hYcmSys.at(kY-1)->SetTitle("y-Component");
    hYcmSys.at(kA-1)->SetTitle("a-Sub-Event");
    hYcmSys.at(kB-1)->SetTitle("b-Sub-Event");
    hYcmSys.at(kC-1)->SetTitle("c-Sub-Event");
    for(int k=0; k<hYcmSys.size(); k++)
    {
        for(int i=0; i<nBinsYcm; i++)
        {
            float value = hYcm.at(k+1)->GetBinContent(i+1);
            float errValue = hYcm.at(k+1)->GetBinError(i+1)/value;
            float united = hYcm.at(kUnited)->GetBinContent(i+1);
            float errUnited = hYcm.at(kUnited)->GetBinError(i+1)/united;
            float sysErr = fabs(value-united)/fabs(value+united)*100; // percent
            float statErr = sqrt(errValue*errValue + errUnited*errUnited)*100; // percent
            float ycm = hYcm.at(kUnited)->GetXaxis()->GetBinCenter(i+1);
            // cout << statErr << endl;
            hYcmSys.at(k)->SetPoint(i, ycm, sysErr);
            hYcmSys.at(k)->SetPointError(i, 0, statErr);
        }
    }
    return hYcmSys;
}

vector<TGraphErrors*> PtSystematics(TFile* file)
{
    enum eProfiles{
        kUnited=0,
        kX,
        kY,
        kA,
        kB,
        kC,
        kNumberOfProfiles
    };
    array<TProfile*, kNumberOfProfiles> hPt;
    hPt.at(kUnited) = new  TProfile( "v1_vs_y", "v_{1} vs y; y_{cm}; v_{1}", 25, 0.0, 2.5 );
    hPt.at(kX) = new TProfile( "v1_vs_pt_X", "v_{1} vs pt; pt; v_{1}^{x}", 25, 0.0, 2.5 );
    hPt.at(kY) = new TProfile( "v1_vs_pt_Y", "v_{1} vs pt; pt; v_{1}^{y}", 25, 0.0, 2.5 );
    hPt.at(kA) = new TProfile( "v1_vs_pt_A", "v_{1} vs pt; pt; v_{1}^{a}", 25, 0.0, 2.5 );
    hPt.at(kB) = new TProfile( "v1_vs_pt_B", "v_{1} vs pt; pt; v_{1}^{b}", 25, 0.0, 2.5 );
    hPt.at(kC) = new TProfile( "v1_vs_pt_C", "v_{1} vs pt; pt; v_{1}^{c}", 25, 0.0, 2.5 );
    for( auto &histo : hPt )
        histo->Sumw2();
    for(int se=0; se<3; se++)
    {
        int wX = ( (TProfile*) file->Get( Form( "/Flow/-0.25 < y_{cm} < -0.15, SE%i_{x}", se) ) )->GetEntries();
        int wY = ( (TProfile*) file->Get( Form( "/Flow/-0.25 < y_{cm} < -0.15, SE%i_{y}", se) ) )->GetEntries();
        hPt.at(kUnited)->Add( (TProfile*) file->Get( Form( "/Flow/-0.25 < y_{cm} < -0.15, SE%i_{x}", se) ), wX );
        hPt.at(kUnited)->Add( (TProfile*) file->Get( Form( "/Flow/-0.25 < y_{cm} < -0.15, SE%i_{y}", se) ), wY );
        hPt.at(kX)->Add( (TProfile*) file->Get( Form( "/Flow/-0.25 < y_{cm} < -0.15, SE%i_{x}", se) ), wX );
        hPt.at(kY)->Add( (TProfile*) file->Get( Form( "/Flow/-0.25 < y_{cm} < -0.15, SE%i_{y}", se) ), wY );
        hPt.at(kA+se)->Add( (TProfile*) file->Get( Form( "/Flow/-0.25 < y_{cm} < -0.15, SE%i_{x}", se) ), wX );
        hPt.at(kA+se)->Add( (TProfile*) file->Get( Form( "/Flow/-0.25 < y_{cm} < -0.15, SE%i_{y}", se) ), wY );
    }
    hPt.at(kA)->Add( (TProfile*) file->Get( "/Flow/-0.25 < y_{cm} < -0.15, SE0_{x}" ) );
    hPt.at(kA)->Add( (TProfile*) file->Get( "/Flow/-0.25 < y_{cm} < -0.15, SE0_{y}" ) );
    hPt.at(kB)->Add( (TProfile*) file->Get( "/Flow/-0.25 < y_{cm} < -0.15, SE1_{x}" ) );
    hPt.at(kB)->Add( (TProfile*) file->Get( "/Flow/-0.25 < y_{cm} < -0.15, SE1_{y}" ) );
    hPt.at(kC)->Add( (TProfile*) file->Get( "/Flow/-0.25 < y_{cm} < -0.15, SE2_{x}" ) );
    hPt.at(kC)->Add( (TProfile*) file->Get( "/Flow/-0.25 < y_{cm} < -0.15, SE2_{y}" ) );
    vector<TGraphErrors*> hPtSys;
    int nBinsPt = hPt.at(kUnited)->GetNbinsX()-10;
    for(int i=0; i<kNumberOfProfiles-1; i++)
    {
        hPtSys.push_back( new TGraphErrors(nBinsPt) );
    }
    hPtSys.at(kX-1)->SetTitle("x-Component");
    hPtSys.at(kY-1)->SetTitle("y-Component");
    hPtSys.at(kA-1)->SetTitle("a-Sub-Event");
    hPtSys.at(kB-1)->SetTitle("b-Sub-Event");
    hPtSys.at(kC-1)->SetTitle("c-Sub-Event");
    for(int k=0; k<hPtSys.size(); k++)
    {
        for(int i=1; i<nBinsPt; i++)
        {
            float value = hPt.at(k+1)->GetBinContent(i+1);
            float errValue = hPt.at(k+1)->GetBinError(i+1)/value;
            float united = hPt.at(kUnited)->GetBinContent(i+1);
            float errUnited = hPt.at(kUnited)->GetBinError(i+1)/united;
            float sysErr = fabs(value-united)/fabs(value+united)*100; // percent
            float statErr = sqrt(errValue*errValue + errUnited*errUnited)*100; // percent
            float pt = hPt.at(kUnited)->GetXaxis()->GetBinCenter(i+1);
            hPtSys.at(k)->SetPoint(i, pt, sysErr);
            hPtSys.at(k)->SetPointError(i, 0, statErr);
        }
    }
    return hPtSys;
}