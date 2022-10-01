//*Creates Histogram  for px py pz in the same pad
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TMath.h>
#include <TLegend.h>
#include <iostream>
#include <fstream>
#include <TStopwatch.h>
#include <TProfile.h>
#include <cmath>

using namespace std;
void PHistGen(){
    TH1D *h1=new TH1D("h1","h1",40,-2,2);
    TH1D *h2=new TH1D("h2","h2", 40, -2, 2);
    TH1D *h3=new TH1D("h3","h3", 40, -2, 2);
    auto legend = new TLegend(0.1,0.7,0.48,0.9);

    ifstream px;
    ifstream py;
    ifstream pz;

    double valuepx;
    double valuepy;
    double valuepz;
    px.open("ofpx.txt");
    py.open("ofpy.txt");
    pz.open("ofpz.txt");
    while(!px.eof()){
        px >> valuepx;
        py >> valuepy;
        pz >> valuepz;
        h1->Fill(valuepx);
        h2->Fill(valuepy);
        h3->Fill(valuepz);

    };
    h1->Draw();
    h1->SetLineColor(kBlue);
    h1->SetTitle("p Distribution; Px; Counts");
    
    h2->Draw("SAME");
    h2->SetLineColor(kRed);
    h2->SetTitle("p Distribution; Py; Counts");

    h3->Draw("SAME");
    h3->SetLineColor(kGreen);
    h3->SetTitle("p Distribution; Pz; Counts");

    legend->SetHeader("The Legend Title","C"); // option "C" allows to center the header
    legend->AddEntry("h1","px","l");
    legend->AddEntry("h2","py","l");
    legend->AddEntry("h3","pz","l");
    legend->Draw();
};