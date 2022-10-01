//*Plotting and saving histograms in root file;
#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TRandom.h>
void RootHisto(){
    TFile *file=new TFile("file.root","recreate"); //!This has to be before the hist...otherwise plot wont come
    TH1D *h = new TH1D("h","First Histogram",40,0,10);  //TH(histogram)1(Dimension)D(Doubles) Name Title NoofBins Begin and end
                                                  //No of bins dhould be square root of the number of entries  
    for (int i=0; i<100000;i++){
        h->Fill(gRandom->Gaus(5,1));
        h->SetMinimum(0);

    };
    file->Write();
    file->Close();

};