#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>

void PlottingFromTreeInRootFile(){
    TFile *infile= new TFile("amptroot.root","READ");
    TTree *tree= (TTree*) infile->Get("tree");
    TH1D *hist= new TH1D("hist","Histogram",20,0,10);
    double x, y;

    tree->SetBranchAddress("x",&x);
    tree->SetBranchAddress("y",&y);

    int TotEntries= tree->GetEntries();
    int i=0;
    while(i<TotEntries){
        tree->GetEntry(i);
        hist->Fill(x);
        i++;
    };
    hist->Draw();
};