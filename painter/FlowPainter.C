void FlowPainter(TString fileName="000")
{
    gStyle->SetTitleSize(0.04,"X");
    gStyle->SetTitleSize(0.04,"Y");
    gStyle->SetTitleOffset(1.6,"Y");
    gStyle->SetTitleOffset(1.0,"X");
    gStyle->SetFrameLineWidth(2);
    gStyle->SetFrameFillColor(0);
    gStyle->SetPadColor(0);
    gStyle->SetLabelSize(0.04,"X");
    gStyle->SetLabelSize(0.04,"Y");
    gStyle->SetPadTopMargin(0.1);
    gStyle->SetPadBottomMargin(0.18);
    gStyle->SetPadLeftMargin(0.18);
    gStyle->SetPadTopMargin(0.1);
    gStyle->SetPadRightMargin(0.1);
    gStyle->SetMarkerSize(1.5);
    gStyle->SetErrorX(0);
    gStyle->SetOptStat(0);
    gStyle->SetCanvasColor(0);
    gStyle->SetTitleFillColor(0);
    gStyle->SetTitleBorderSize(0);
    gStyle->SetHistLineWidth(3);
    gStyle->SetLineWidth(2);
    gStyle->SetCanvasBorderMode(0);
    gStyle->SetLegendBorderSize(0);
    gStyle->SetPadBorderMode(0);


    auto file = new TFile(fileName+".root");
    vector<TProfile*> xRapidity;
    vector<TProfile*> yRapidity;
    vector<TProfile*> xPt;
    vector<TProfile*> yPt;
    vector<THStack*> stack;
    
    stack.push_back( new THStack("v1 vs ycm", "v_{1} vs y_{cm}; y_{cm}; v_{1}") );
    stack.push_back( new THStack("v1 vs pt", "v_{1} vs pt; pt, #frac{GeV}{c}; v_{1}") );
    
    for(int se=0; se<3; se++)
    {
        xRapidity.push_back( (TProfile*) file->Get( Form( "/Flow/0.8 < pt < 0.85 SE%i_{x}", se ) ) );
        stack.at(0)->Add(xRapidity.back());
        yRapidity.push_back( (TProfile*) file->Get( Form( "/Flow/0.8 < pt < 0.85 SE%i_{y}", se ) ) );
        stack.at(0)->Add(yRapidity.back());
        xPt.push_back( (TProfile*) file->Get( Form("/Flow/-0.25 < y_{cm} < -0.15, SE%i_{x}", se) ) );
        stack.at(1)->Add(xPt.back());
        yPt.push_back( (TProfile*) file->Get( Form("Flow/-0.25 < y_{cm} < -0.15, SE%i_{y}", se) ) );
        stack.at(1)->Add(yPt.back());
    }
    for(int i=0; i<3; i++)
    {
        xRapidity.at(i)->SetLineColor(i+1);
        xRapidity.at(i)->SetMarkerColor(i+1);
        xRapidity.at(i)->SetMarkerStyle(20+i);
        xRapidity.at(i)->SetMarkerSize(1);
        xRapidity.at(i)->SetLineWidth(2);
        xRapidity.at(i)->GetYaxis()->SetRangeUser(-0.6,0.6);

        yRapidity.at(i)->SetLineColor(i+1);
        yRapidity.at(i)->SetMarkerColor(i+1);
        yRapidity.at(i)->SetMarkerStyle(24+i);
        yRapidity.at(i)->SetMarkerSize(2);
        yRapidity.at(i)->SetLineWidth(2);
        yRapidity.at(i)->GetYaxis()->SetRangeUser(-0.6,0.4);

        xPt.at(i)->SetLineColor(i+1);
        xPt.at(i)->SetMarkerColor(i+1);
        xPt.at(i)->SetMarkerStyle(20+i);
        xPt.at(i)->SetMarkerSize(1);
        xPt.at(i)->SetLineWidth(2);
        xPt.at(i)->GetYaxis()->SetRangeUser(-0.27,-0.07);

        yPt.at(i)->SetLineColor(i+1);
        yPt.at(i)->SetMarkerColor(i+1);
        yPt.at(i)->SetMarkerStyle(24+i);
        yPt.at(i)->SetMarkerSize(2);
        yPt.at(i)->SetLineWidth(2);
        yPt.at(i)->GetYaxis()->SetRangeUser(-0.27,-0.07);
    }
    auto legend = new TLegend(0.2,0.65,0.5,0.8);
    for(int se=0; se<3; se++)
    {
        legend->AddEntry(xPt.at(se), Form("SE%i, X", se) );
        legend->AddEntry(yPt.at(se), Form("SE%i, Y", se) );
    }
    auto canvas = new TCanvas( "canvas", "", 1500, 700 );
    canvas->Divide(2,1);
    canvas->cd(1);
    stack.at(0)->Draw("NOSTACK");
    legend->Draw();
    canvas->cd(2);
    stack.at(1)->Draw("NOSTACK");
    legend->Draw();
    canvas->SaveAs(fileName+".png");
}