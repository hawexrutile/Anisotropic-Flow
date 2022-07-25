//*Code to create centrality-wise v2,v3 vs pt  distributin.
#define check_cxx
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
#include <TSystem.h>


//?Why are there two similar trees in root file?
using namespace std;             //Global variables have been defined bellow for use in multiple functions
const int MultBin=50;                  //Multiplicity Bin...   Total no of bins in the centrality distribution plot
const int N=5;                   //Total no of Centrality Regions
const int BINpls=10;             //Bins for TProfile

double lowcentral[N]={0.0};
double highcentral[N]={0.0};         //Global variable to be used i, both functions

TFile *CenRegRoot=new TFile("CenRegRoot.root","READ");             //!
TH1D *pion= (TH1D*) CenRegRoot->Get("pion");
TProfile *g3 =new TProfile("g3", "Centrality vs v3",N,0,100);      

//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ien
void check::Loop(int ij,TCanvas *canva ){           //SECTION //*Main loop function to calcualte and plot v2

   if (fChain == 0) return;
   Long64_t nentries = fChain->GetEntriesFast();             //total no of events        
   Long64_t nbytes = 0, nb = 0;
   
   canva->Divide(1,1);                                                //!
   canva->cd(1);
   
   TProfile *g2 =new TProfile("g2", "v2 vs Pt  ",BINpls,0,3.5);   //!dont specidfy y range here nstead here LINK 6 AMPT Practise\6 v2 ROOT Method\check.C:79
   TProfile *g1 =new TProfile("g1", "v3 vs Pt",BINpls,0,3.5);    //

   double lowerbc=pion->GetBinCenter(lowcentral[ij]);                                           //For a given region gives the lowest bin's bin centre(ie the total no of trails)
   double higherbc=pion->GetBinCenter(highcentral[ij]);

   int TotEvents=nentries;                                             // *TotEvents is for easy adjustament of event size 
   for (Long64_t jentry=0; jentry<TotEvents;jentry++) { 
//                                                                       Total events=67295
      Long64_t ientry = LoadTree(jentry);                              //Load Tree apparently gives back the the respective possition in the tree
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;                   //!Loads all the variables(branches) of the ith leaf
      // if (Cut(ientry) < 0) continue;
      double pt;
      double phi;
      // cout<<Mult<<"\n";
      if (jentry%10000==0) cout<< floor((ij*1.0+(jentry*1.0/TotEvents)) *100/(N*1.0)) <<"% complete \n";   //*Percentage finished indicator

      int   chrgi=0;                                     //Charged i entires
      for (int p =0;p < Mult ;p++){
         double eta= log(abs((Px[p]+Py[p])/(Px[p]+Py[p]+(2*Pz[p]))));
         if (!(((eta>2.8) && (eta<5.1)) || ((eta> -3.7) && (eta< -1.7)) ) ) continue;     // 2.8<eta<5.1 or -3.7<eta<-1.7 should be allowed.  //!Not-operator is important
         if(abs(PID[p])!=211 && abs(PID[p])!=321 && abs(PID[p])!=2212) continue;   //We only need pion, proton or kaon otherwise you may skip....
         chrgi++ ;                                    //Charged i entires
      };
      
      if ((chrgi<lowerbc) || (chrgi>higherbc)) continue;                      //to loop on the given centralities //! It should be charged particle not mult.
      for (int p =0;p < Mult ;p++){
         double eta= log(abs((Px[p]+Py[p])/(Px[p]+Py[p]+(2*Pz[p]))));
         if (!((eta> -0.8) && (eta<0.8)) ) continue;                         //CONFINING WITH ETA

         if(abs(PID[p])!=211 && abs(PID[p])!=321 && abs(PID[p])!=2212 ) continue;//We only need pion, proton or kaon otherwise you may skip....
         if ((Px[p]==0) && (Py[p]==0)) continue;        //We need  0 pt so && should suffice...It diddnt make a difference when I changed from "||" to "&&" so it should be safe to change;
         pt=sqrt(pow(Px[p],2)+pow(Py[p],2));  //for that particular event        //Pt has units of GeV/c^2
         phi=(TMath::ATan2(Py[p],Px[p]));
         if (phi<0.0)phi+=2*TMath::Pi();
         g2->Fill(pt,cos(2*(phi-Psi)));                             //*v2
         g1->Fill(pt,cos(3*(phi-Psi)));                             //*v3
         g1->GetYaxis()->SetRangeUser(-0.02, 0.3);                                //!Set Y range here
         g2->GetYaxis()->SetRangeUser(-0.02, 0.3);                                
         g3->Fill((ij*20)+10,cos(2*(phi-Psi)));

      };
   };
   g1->SetTitle("Pt vs v2 and v3; pt; v2/v3");
   g1-> SetMarkerStyle(23);      
   g2-> SetMarkerStyle(26);       
   g3-> SetMarkerStyle(26);       
   g1->Draw();
   g2->Draw("SAME");

};

void check::Multiplicity(){                                //ANCHOR Creates a multiplicity distribution and saves the hostpgram on a root file
   if (fChain == 0) return;
   fChain->SetBranchAddress("Px",&Px);
   fChain->SetBranchAddress("Py",&Py);
   Long64_t nentries = fChain->GetEntriesFast();
   Long64_t nbytes = 0, nb = 0;
   int TotEvents=nentries;   //Total events=67295           //?wHY tOTeVENTS
   
   TFile *CenRegRoot=new TFile("CenRegRoot.root","recreate");
   TH1D *pion= new TH1D("pion","Pion Multiplicity",MultBin,390,2000);   // !int maxTrack=395; int minTrack=2000 after eta confinemnt. it was 5000 bfr; 395 min is imp
//######################################################################
   for (Long64_t jentry=0; jentry<TotEvents;jentry++) {
      Long64_t ientry = LoadTree(jentry);                  //Load Tree apparently gives back the the respective possition in the tree
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      int   chrgi=0;                                     //Charged i entires
      for (int p =0;p < Mult ;p++){
         double eta= log(abs((Px[p]+Py[p])/(Px[p]+Py[p]+(2*Pz[p]))));
         if (!(((eta>2.8) && (eta<5.1)) || ((eta> -3.7) && (eta< -1.7)) ) ) continue;     // 2.8<eta<5.1 or -3.7<eta<-1.7 should be allowed.  //!Not-operator is important
         if(abs(PID[p])!=211 && abs(PID[p])!=321 && abs(PID[p])!=2212) continue;   //We only need pion, proton or kaon otherwise you may skip....
         chrgi++ ;                                    //Charged i entires
      };
      if (jentry%5000==0) cout<< floor(((jentry*1.0)/TotEvents) *100) <<"% complete \n";
      pion->Fill(chrgi);//For countinf charged particles(the big 3s)
      pion->Draw();
//########################################################################
   };
   CenRegRoot->Write();
   CenRegRoot->Close();
};

void check::Printer(){                                           //ANCHOR uses the multiplicity distribution and runs the 4-loop on Loop() function for all centralities...

   double_t *CenReg;
   CenReg= pion->GetIntegral();              //!cumulative sum of bins ...The resulting integral is normalized to 1
   // for (int i=0 ;i<MultBin;i++){
   //    cout<<"integral till bin "<<i<<" is "<<CenReg[i]<<"\n";
   // };

   int w1=0,w2=0,w3=0,w4=0,w5=0;         //Locks
   for (int i=0;i<=MultBin;i++){           //?Why is i=2; cuz in protay's code initial bins have to be skiped ; This loop is to set the starting and ending bin of each centrality region(eg 0%-5%)
                                          //!The lowest if-block gets satisfied by the first bin
      if (CenReg[i]<=1){		  //keeps on changing till the last bin
			highcentral[0]=i;      //Highest bin  of the highest centrality region(0-20%)
		};
		if (CenReg[i]>0.8){  
			if (w1==0){            //locks last
				w1=1;
				lowcentral[0]=i;    //Lowest Bin of the highest centrality region(0-20%)
				highcentral[1]=i;   //Highest bin  of the sevcond highest centrality region(20-40%)
            // cout<<1;
			};
		};
		if (CenReg[i]>0.6){
			if (w2==0){
				w2=1;
				lowcentral[1]=i;
				highcentral[2]=i;
            // cout<<2;
			};
		};
		if (CenReg[i]>0.4){
			if (w3==0){
				w3=1;
				lowcentral[2]=i;
				highcentral[3]=i;
            // cout<<3;
			};
		};
		if (CenReg[i]>0.2){
			if (w4==0){                    //locks second
				w4=1;
			   lowcentral[3]=i;               //Lowest bin  of the second  lowest centrality region(60-80%) 
			   highcentral[4]=i;              //Highest bin  of the lowest centrality region(80-100%) 
            // cout<<4;
			};
		};
		if (CenReg[i]>=0){                //Locks first   
			if (w5==0){
				w5=1;
				lowcentral[4]=i;            //Lowest bin  of the lowest centrality region(80-100%) 
            // cout<<5;
			};
		};
	};

   for (int i=0;i<N;i++){ 
      // cout <<lowcentral[i]<< "\t"<<highcentral[i]<<"\t"<<"\n";
                                  //*For ploting all the plots  
      TCanvas *canva = new TCanvas("canva","Pt vs v2 and v3 for Proton, Pion and Kaon",1000,1000); //I made a Tcanvas cuz to use the print function
      if (i==0){

         this->Loop(i,canva);
                        //For calling functions within the same class
         canva->Print("./Plot/g0.png");        //TODO remove memory leak warning
         if (canva){ 
            canva->Close(); 
            gSystem->ProcessEvents(); 
            delete canva; 
            canva = 0; 
         };
         
      }      
      else if (i==1){

         this->Loop(i,canva);                           //For calling functions within the same class
         canva->Print("./Plot/g1.png");        //!
         if (canva){ 
            canva->Close(); 
            gSystem->ProcessEvents(); 
            delete canva; 
            canva = 0; 
         };
      } 
      else if (i==2){ 
         this->Loop(i,canva);                           //For calling functions within the same class
         canva->Print("./Plot/g2.png");        //!
         if (canva){ 
            canva->Close(); 
            gSystem->ProcessEvents(); 
            delete canva; 
            canva = 0; 
         };
      }      
      else if (i==3){  
         this->Loop(i,canva);                           //For calling functions within the same class
         canva->Print("./Plot/g3.png");        //!
         if (canva){ 
            canva->Close(); 
            gSystem->ProcessEvents(); 
            delete canva; 
            canva = 0; 
         };
      }      
      else {
         this->Loop(i,canva);                           //For calling functions within the same class
         canva->Print("./Plot/g4.png");        //!
         g3->SetTitle("Centrality vs v2; Centrality; v2");
         g3->Draw();
         canva->Print("./Plot/Centvsv2.png");      //TODO do we need to close; does it rewrite files
         if (canva){ 
            canva->Close(); 
            gSystem->ProcessEvents(); 
            delete canva; 
            canva = 0; 
         };
      } ;       
   };

};
