#include <TH2.h>
#include <TF1.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TRandom.h>
#include <TLegend.h>
#include <TLine.h>
#include <TArrow.h>
#include <TLatex.h>
#include <iostream>

using namespace std;

void hist(){
    TCanvas *canva = new TCanvas();
    canva->SetTickx();                     //covers the axis above the x axis with "scale"
    canva->SetTicky();                     //covers the axis beside the y axis with "scale"  
    canva->SetGridx();
    canva->SetGridy();

    TH1D *hist = new TH1D("hist","First Histogram",40,0,15);  //TH(histogram)1(Dimension)D(Doubles) Name Title NoofBins Begin and end //*No of bins dhould be square root of the number of entries
    for (int i=0; i<100000;i++){
        hist->Fill(gRandom->Gaus(5,1));
    };

    TF1 *fit= new TF1("fit","gaus",0,15);
    fit->SetParameter(0,40) ;      //Bin No
    fit->SetParameter(1,5);        //Mean
    fit->SetParameter(2,1) ;       //std dev ig
    fit->SetLineWidth(3);               //fit-line thickness
    fit->SetLineColor(kRed);               //fit-line thickness
    fit->SetLineStyle(3);               //1-solid ; 2- dotted

    double mean=fit->GetParameter(1);
    double std=fit->GetParameter(2);


    hist->SetStats(0);                   //*for removing the stat box on the top right corner
    hist->SetMinimum(0);
    hist->Draw("P*");
    hist->SetFillColor(kBlue-9);
    hist->Fit("fit","R");           //R stands for region of range :the range aplied here:LINK 5 Root Practice\1 Basics\4.1 Histograms\hist.C:21
                                    //Q is for quiet; "0" if for not drawing
    hist->SetTitle("First Histogram;A;B");               //!axis title can be set on the top
    hist->SetTitleSize(0.05,"X");                       //LINK https://youtu.be/-uugg0wshzc?list=PLLybgCU6QCGWLdDO4ZDaB0kLrO3maeYAe&t=342
    hist->SetTitleSize(0.05,"Y");    

    TLegend *leg=new TLegend(0.5,0.6,0.8,0.8);            //Consider the canvas as an cordinate region withe the top right corner as 1,1 and bottom left corner as 0,0;Cordinate (x1,y1,x2,y2)  
    leg->SetBorderSize(0);                                 //sets the margin to 0
    leg->AddEntry(hist,"Measured Data","f");               //L IS FOR LINE; F OR FILLING ; P FOR POINT
    leg->AddEntry(fit,"Fit Function","l");
    leg->Draw();
    

    TLine *line=new TLine(0,2000,15,2000);                    //Cordinate (x1,y1,x2,y2)
    line->SetLineWidth(2);
    line->SetLineColor(kBlack);
    line->Draw();

    double x0=6.0;
    int bin=hist->FindBin(x0);
    double y0=hist->GetBinContent(bin);

    TArrow  *arrow = new TArrow(10.0,4000,x0,y0);
    arrow->SetLineWidth(2);
    arrow->SetArrowSize(0.02);
    arrow->SetLineColor(kOrange);
    arrow->Draw();

    TLatex *text = new TLatex(10.0,4000,"lol");
    text->Draw();

    cout<< "Counts in each bin\n";
    for (int i=1; i<=hist->GetNbinsX(); i++){
        cout<<hist->GetBinCenter(i) <<"\t"<< hist->GetBinContent(i)<<"\n";  //tab is compulsory
    }
}