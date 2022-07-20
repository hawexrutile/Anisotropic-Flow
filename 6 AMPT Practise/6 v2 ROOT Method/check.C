#define check_cxx
#define BINpls 10
#include <algorithm>
#include "check.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TMath.h>
#include <iostream>
#include <TStopwatch.h>
#include <TProfile.h>
#include <cmath>


//?Why are theri two similar trees in root file?
using namespace std;
int MultBin=50;                  //Multiplicity Bin...   Total no of bins in the centrality distribution plot
const int N=8;
double lowcentral[N]={0.0};
double highcentral[N]={0.0};         //Global variable to be used i both functions

//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ien
void check::Loop(double lowerbc, double higherbc){

   if (fChain == 0) return;
   fChain->SetBranchAddress("Px",&Px);
   fChain->SetBranchAddress("Py",&Py);//?wHY DIISNT WE DO THIS FOR  PID
   Long64_t nentries = fChain->GetEntriesFast();
   Long64_t nbytes = 0, nb = 0;
   TCanvas *canva = new TCanvas("canva","Pt vs v2 and v3 for Proton, Pion and Kaon",1000,1000);

   //########################################################################
   //########################################################################
   

   
   //###########################################################################
   //###########################################################################

   // int canvasno=8;                                        //tOTAL nO OF plots
   // string canvaplotlist[8]={"g1","g2","g3","g4","g5","g6","g7","g8"};                 //List of plot ID's to use in giving plots positions
   // TCanvas * twntypcs=TCanvas("twntypcs","Twenty Pieces",1000,1000) ;                //To make 10 plots //?No idea why the error squigle fix it;
   // twntypcs->Divide(5,2);

   // for (int l;l<N;l++){                                       //For each Centrality Region


      canva->Divide(1,1);
      canva->cd(1);
      TProfile *g2 =new TProfile("g2", "Pt vs v2 ",BINpls,0,5);
      TProfile *g1 =new TProfile("g2", "Pt vs v3",BINpls,0,5);
      int TotEvents=nentries;// *TotEvents is for easy adjustament of event size 
      for (Long64_t jentry=0; jentry<TotEvents;jentry++) { //Total events=67295
         Long64_t ientry = LoadTree(jentry);                  //Load Tree apparently gives back the the respective possition in the tree
         if (ientry < 0) break;
         nb = fChain->GetEntry(jentry);   nbytes += nb;       //!Loads all the variables(branches) of the ith leaf
         // if (Cut(ientry) < 0) continue;
         double pt;
         double phi;
         // cout<<Mult<<"\n";
         if (jentry%5000==0) cout<< floor(((jentry*1.0)/TotEvents) *100) <<"% complete \n";
         
         if ((Mult<lowerbc) || (Mult>higherbc)) continue;
         for (int p =0;p < Mult ;p++){

            if(abs(PID[p])!=211 && abs(PID[p])!=321 && abs(PID[p])!=2212 ) continue;//We only need pion, proton or kaon otherwise you may skip....
            if ((Px[p]==0) && (Py[p]==0)) continue;        //?Do I need this"||"?---We need  0 pt so && should suffice...It diddnt make a difference when I changed from "||" to "&&" so it should be safe to change;
            pt=sqrt(pow(Px[p],2)+pow(Py[p],2));  //for that particular event        //Pt has units of GeV/c^2
            phi=(TMath::ATan2(Py[p],Px[p]));
            g2->Fill(pt,cos(2*(phi-Psi)));
            g1->Fill(pt,cos(3*(phi-Psi)));
         };


      };
      g1->SetTitle("Pt vs v2 and v3 for Proton, Pion and Kaon; pt; v2/v3");
      // g2->SetTitle("Pt vs v2 and v3 for Proton, Pion and Kaon; pt; v2/v3");
      g1-> SetMarkerStyle(23);       //TODO check if marker set style is working
      g2-> SetMarkerStyle(26);       //TODO check if marker set style is working
      g1->Draw();
      g2->Draw("SAME");
      canva->Print("lol.png");
      
};

void check::Multiplicity(){
   if (fChain == 0) return;
   fChain->SetBranchAddress("Px",&Px);
   fChain->SetBranchAddress("Py",&Py);
   Long64_t nentries = fChain->GetEntriesFast();
   Long64_t nbytes = 0, nb = 0;
   int TotEvents=nentries;
   int TotTracks=0;
   
   TFile *CenRegRoot=new TFile("CenRegRoot.root","recreate");
   TH1D *pion= new TH1D("pion","Pion Multiplicity",MultBin,395,5000);   // *int maxTrack=395; int minTrack=9000;
//######################################################################
   for (Long64_t jentry=0; jentry<TotEvents;jentry++) { //Total events=67295
      Long64_t ientry = LoadTree(jentry);                  //Load Tree apparently gives back the the respective possition in the tree
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      int   chrgi=0;
      TotTracks=TotTracks+Mult;                                     //Charged i entires
      for (int p =0;p < Mult ;p++){
         if(abs(PID[p])!=211 && abs(PID[p])!=321 && abs(PID[p])!=2212) continue;//We only need pion, proton or kaon otherwise you may skip....
         chrgi++ ;                                    //Charged i entires
      };
      if (jentry%5000==0) cout<< floor(((jentry*1.0)/TotEvents) *100) <<"% complete \n";
      pion->Fill(chrgi);//For countinf charged particles(the big 3s)
//########################################################################
   };

   CenRegRoot->Write();
   CenRegRoot->Close();

   // pion->Draw();
   // pion->SetLineColor(kBlue);
   // pion->SetTitle("Pion Multiplicity;Counts(N);Multiplicity");
   
   
};

void check::Printer(){
   TFile *CenRegRoot=new TFile("CenRegRoot.root","READ");
   TH1D *pion= (TH1D*) CenRegRoot->Get("pion");
   double_t *CenReg;
   CenReg= pion->GetIntegral();


	int w1=0,w2=0,w3=0,w4=0,w5=0,w6=0,w7=0,w8=0;
   for (int i=0;i<=MultBin;i++){           //*Why is i=2; ig i is the number of bins ; This loop is to set the starting and ending bin of each centrality region(eg 0%-5%)
		if (CenReg[i]<=1){			//The lowest if-block gets satisfied by the first bin
			highcentral[0]=i;
		};
		if (CenReg[i]>0.95){
			if (w1==0){
				w1=1;
				lowcentral[0]=i;
				highcentral[1]=i;
			};
		};
		if (CenReg[i]>0.9){
			if (w2==0){
				w2=1;
				lowcentral[1]=i;
				highcentral[2]=i;
			};
		};
		if (CenReg[i]>0.8){
			if (w3==0){
				w3=1;
				lowcentral[2]=i;
				highcentral[3]=i;
			};
		};
		if (CenReg[i]>0.7){
			if (w4==0){
				w4=1;
			lowcentral[3]=i;
			highcentral[4]=i;
			};
		};
		if (CenReg[i]>0.6){
			if (w5==0){
				w5=1;
				lowcentral[4]=i;
				highcentral[5]=i;
			};
		};
		if (CenReg[i]>0.5){
			if (w6==0){
				w6=1;
				lowcentral[5]=i;
				highcentral[6]=i;
			};
		};
		if (CenReg[i]>0.4){
			if (w7==0){
				w7=1;
				lowcentral[6]=i;
				highcentral[7]=i;
			};
		};
		if (CenReg[i]>0.2){
			if (w8==0){
				w8=1;
				lowcentral[7]=i;
			};
		};
	};
   for (int i=0;i<N;i++){
      // vector<string> canvaplotlist={"g1.png","g2.png","g3.png","g4.png","g5.png","g6.png","g7.png","g8.png"};                 //List of plot ID's to use in giving plots positions
      double lowerbc=pion->GetBinCenter(lowcentral[i]);                                           //For a given region gives the lowest bin's bin centre(ie the total no of trails)
      double higherbc=pion->GetBinCenter(highcentral[i]);   
      this->Loop(lowerbc,higherbc);
      // string lol=string(canvaplotlist[i]);
      
   };
};
