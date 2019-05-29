void ChannelErrorPainter()
{
    gROOT->ProcessLine(".x SetStyle.C");
    vector<TFile*> file{
        new TFile("000.root"),
        new TFile("001.root"),
        new TFile("gt120.root")
    };
    vector<TProfile*> hYcm{
        new TProfile( "80<FW-adc, ycm", ";y_{cm};v_{1}", 14, -0.7, 0.7 ),        //0
        new TProfile( "80<FW-adc<120, ycm", ";y_{cm};v_{1}", 14, -0.7, 0.7 ),    //1
        new TProfile( "120<FW-adc, ycm", ";y_{cm};v_{1}", 14, -0.7, 0.7 )        //2
    };
    vector<TProfile*> hPt{
        new TProfile( "80<FW-adc, pt", ";pt;v_{1}", 25, 0.0, 2.5 ),        //0
        new TProfile( "80<FW-adc<120, pt", ";pt;v_{1}", 25, 0.0, 2.5 ),    //1
        new TProfile( "120<FW-adc, pt", ";pt;v_{1}", 25, 0.0, 2.5 )        //2
    };
    for( auto &histo : hYcm )
        histo->Sumw2();
    for( auto &histo : hPt )
        histo->Sumw2();
    for( int i=0; i<file.size(); i++ )
    {
        for( int se=0; se<3; se++ )
        {
            int wX = ( (TProfile*) file.at(i)->Get(Form("/Flow/0.8 < pt < 0.85 SE%i_{x}", se)) )->GetEntries(); 
            int wY = ( (TProfile*) file.at(i)->Get(Form("/Flow/0.8 < pt < 0.85 SE%i_{y}", se)) )->GetEntries();
            hYcm.at(i)->Add( (TProfile*) file.at(i)->Get( Form("/Flow/0.8 < pt < 0.85 SE%i_{x}", se) ), wX ); 
            hYcm.at(i)->Add( (TProfile*) file.at(i)->Get( Form("/Flow/0.8 < pt < 0.85 SE%i_{y}", se) ), wY ); 

            wX = ( (TProfile*) file.at(i)->Get(Form("/Flow/-0.25 < y_{cm} < -0.15, SE%i_{x}", se)) )->GetEntries();
            wY = ( (TProfile*) file.at(i)->Get(Form("/Flow/-0.25 < y_{cm} < -0.15, SE%i_{x}", se)) )->GetEntries();
            hPt.at(i)->Add( (TProfile*) file.at(i)->Get( Form("/Flow/-0.25 < y_{cm} < -0.15, SE%i_{x}", se) ), wX ); 
            hPt.at(i)->Add( (TProfile*) file.at(i)->Get( Form("/Flow/-0.25 < y_{cm} < -0.15, SE%i_{y}", se) ), wY ); 
        }
    }
    int nbinsYcm = hYcm.at(0)->GetNbinsX();
    vector<TGraphErrors*> hYcmSys{
        new TGraphErrors(nbinsYcm),
        new TGraphErrors(nbinsYcm)
    };
    hYcmSys.at(0)->SetTitle( "80 < adc < 120" );
    hYcmSys.at(1)->SetTitle( "120 < adc" );
    for( int i=1; i<hYcm.size(); i++ )
    {
        for( int bin=0; bin<nbinsYcm; bin++ )
        {
            float binValue = hYcm.at(i)->GetBinContent(bin+1);
            float binError = hYcm.at(i)->GetBinError(bin+1)/binValue;
            float refValue = hYcm.at(0)->GetBinContent(bin+1);
            float refError = hYcm.at(0)->GetBinError(bin+1)/refValue;
            float sysErr = fabs(binValue-refValue)/fabs(binValue+refValue);
            float statErr = sqrt( binError*binError + refError*refError );
            float ycm = hYcm.at(0)->GetXaxis()->GetBinCenter(bin+1);
            hYcmSys.at(i-1)->SetPoint( bin, ycm, sysErr*100 );
            hYcmSys.at(i-1)->SetPointError( bin, 0, statErr*100 );
        }
    }
    int nbinsPt = hPt.at(0)->GetNbinsX() - 10;
    vector<TGraphErrors*> hPtSys{
        new TGraphErrors(nbinsPt),
        new TGraphErrors(nbinsPt)
    };
    hPtSys.at(0)->SetTitle( "80 < adc < 120" );
    hPtSys.at(1)->SetTitle( "120 < adc" );
    for( int i=1; i<hPt.size(); i++ )
    {
        for( int bin=1; bin<nbinsPt; bin++ )
        {
            float binValue = hPt.at(i)->GetBinContent(bin+1);
            float binError = hPt.at(i)->GetBinError(bin+1)/binValue;
            float refValue = hPt.at(0)->GetBinContent(bin+1);
            float refError = hPt.at(0)->GetBinError(bin+1)/refValue;
            float sysErr = fabs(binValue-refValue)/fabs(binValue+refValue);
            float statErr = sqrt( binError*binError + refError*refError );
            float pt = hPt.at(0)->GetXaxis()->GetBinCenter(bin+1);
            hPtSys.at(i-1)->SetPoint( bin, pt, sysErr*100 );
            hPtSys.at(i-1)->SetPointError( bin, 0, statErr*100 );
        }
    }
    vector<TMultiGraph*> stack{
        new TMultiGraph("y_{cm}", "Configuration Systematics;y_{cm};Relative Error [%]"),
        new TMultiGraph("pt", "Configuration Systematics;pt, [#frac{GeV}{c}];Relative Error [%]")
    };
    for(int i=0; i<hYcmSys.size(); i++)
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
        stack.at(1)->Add(hPtSys.at(i));
    }
    vector<THStack*> flow{
        new THStack("ycm-flow", ";y_{cm};v_{1}"),
        new THStack("pt-flow", ";pt, [#frac{GeV}{c}];v_{1}")
    };
    for(int i=0; i<hYcm.size(); i++)
    {
        hYcm.at(i)->GetYaxis()->SetRangeUser(-0.6, 0.6);
        hYcm.at(i)->SetMarkerStyle(20+i);
        hYcm.at(i)->SetMarkerSize(2);
        hYcm.at(i)->SetMarkerColor(i+1);
        hYcm.at(i)->SetLineColor(i+1);

        hPt.at(i)->GetYaxis()->SetRangeUser(-0.3, 0.05);
        hPt.at(i)->SetMarkerStyle(20+i);
        hPt.at(i)->SetMarkerSize(2);
        hPt.at(i)->SetMarkerColor(i+1);
        hPt.at(i)->SetLineColor(i+1);
        
        flow.at(0)->Add(hYcm.at(i));
        flow.at(1)->Add(hPt.at(i));
    }
    vector<TCanvas*> canvas{ new TCanvas( "errors", "", 2000, 800 ), new TCanvas( "flow", "", 2000, 800 ) };
    for( auto canv : canvas )
        canv->Divide(2,1);
    int i=1;
    for( auto graph : stack )
    {
        canvas.at(0)->cd(i);
        graph->GetYaxis()->SetRangeUser(-5.0, 15.0);
        graph->Draw("A");
        gPad->BuildLegend();
        i++;
    }
    i=1;
    for( auto histo : flow )
    {
        canvas.at(1)->cd(i);
        histo->Draw("NOSTACK");
        gPad->BuildLegend();
        i++;
    }
    canvas.at(0)->SaveAs("Configuration Systematics.png");
    canvas.at(1)->SaveAs("Comparasion Flow.png");
}