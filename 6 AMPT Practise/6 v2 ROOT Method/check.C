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

//?Why are theri two similar trees in root file?
using namespace std;

//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ien
void check::Loop(){
   double binarray[BINpls];                    //! total no of bins has to be set here(dont use the bin as a variable cuz its an array instead use as a macro)
   double v2array[BINpls];                     //! same here

   if (fChain == 0) return;
   fChain->SetBranchAddress("Px",&Px);
   fChain->SetBranchAddress("Py",&Py);

   Long64_t nentries = fChain->GetEntriesFast();

   

   // int BINpls=10;                                                //set BIN value here; the number of points you want to plot
   double bw=5/(1.0*BINpls);                               //I subbed (Total width of plot) with 5 cuz it was 5 for my prev plot and other plots I saw online
   int i=0;
   double Bmin=0;
   double Bmax=0;

   while (i<BINpls){               //for each pt-BIN we run the loop on all the trails
   TStopwatch time;
   cout<<((i*1.0)/BINpls) *100 <<"% complete \n";
      double v2=0;
      double sum=0;
      int j=0;
      int lpos1=0;             //looping int fot calculationg Psi
      int lpos2=0;             //looping int for calculating tot
      int Tot=0;               //Variable to calculate v2
      Bmin=i*(bw);
      Bmax=(i+1)*(bw);
      double bc=((2.0*i+1.0)/2.0)*bw;    //BIN centre; converting to double is imp as otherwise it doesnt give the proper bin centre
      Long64_t nbytes = 0, nb = 0;
      for (Long64_t jentry=0; jentry<nentries-60000;jentry++) { //Tot events=67295
         // cout<<"event loop started \n";
         Long64_t ientry = LoadTree(jentry);                  //Load Tree apparently gives back the the respective possition in the tree
         if (ientry < 0) break;
         nb = fChain->GetEntry(jentry);   nbytes += nb;
         // if (Cut(ientry) < 0) continue;
         double pt;
         double phi;
         double SinSum=0;                                         //?Is the Psi defined by check.h actualy event plane?...yeah!!but ampt normaly keeps it 0 for this root file
         double CosSum=0;
         for (int p =0;p < Mult ;p++){   //This while loop is to calculate Psi. Here listtrail is a list of no events in each trail
            // pt=sqrt(pow(Px[p],2)+pow(Py[p],2));  //for that particular event
            if(abs(PID[p])!=211 && abs(PID[p])!=321 && abs(PID[p])!=2212) continue;//We only need pion, proton or kaon otherwise you may skip.This was important as it had a sihnificant change in the figure .....gave the maximum near 2 where as previously it was 1
            if ((Px[p]==0) && (Py[p]==0)) continue;
            phi=(TMath::ATan2(Py[p],Px[p]));
            // if (phi<0.0){         //To accomodate for inver trig anomalies
            //    phi+=2*TMath::Pi();
            // };
            SinSum=SinSum +sin(2*phi);        //summation of cos and sin phi
            CosSum=CosSum +cos(2*phi);        //?should Psi be for all the trails in the event or for the particular particle or for that particular ptbin?

         };
         // if ((SinSum==0)  (CosSum==0)) continue; //i kept these many "and if "statements cuz their was an inf/nan propagation....I fixed that now
         double Psi=(1/2)*TMath::ATan2(SinSum,CosSum);         //for each event we calcaulate Psi
         // if (Psi<0.0){         //To accomodate for inver trig anomalies
         //       Psi+=2*TMath::Pi();
         // };
         if (Psi!=0||Psi!=-0){
         cout<< Psi<<"\n";

         };
         // cout<< sizeof(Px)/sizeof(float_t)<<"\n";
         for (int p =0;p < Mult ;p++){
            if(abs(PID[p])!=211 && abs(PID[p])!=321 && abs(PID[p])!=2212) continue;//We only need pion, proton or kaon otherwise you may skip....
            if ((Px[p]==0) && (Py[p]==0)) continue;        //?Do i need this"||"?---We need  0 pt so && should suffice...It diddnt make a difference when I changed from "||" to "&&" so it should be safe to change;
            pt=sqrt(pow(Px[p],2)+pow(Py[p],2));  //for that particular event
            phi=(TMath::ATan2(Py[p],Px[p]));
            // if (phi<0.0){         //To accomodate for inverse trig anomalies
            //    phi+=2*TMath::Pi();
            // };
            // cout<<pt<<" ";
            if ((pt>=Bmin)&&(pt<Bmax)){     //if the ith particle's pt belongs to the BIN    //if phi is greater than 0 ; cuz we are averaging from 0-2pi
                  if (phi>TMath::Pi()){
                     cout<<phi<<"\t";
      
                  };
                  // cout<<"sumloop running\n";
                  sum=sum+cos(2*(phi-Psi));               
                  Tot++;
            };
         };

         // cout<<"event loop ended \n";
      };

      v2=sum/Tot;
      // cout<<sum<<"\t"<<Tot<<"\n";
      binarray[i]=bc;           //These arrays are used for plotting
      v2array[i]=v2;
      i++;
      cout<<time.RealTime()<<"s took for bin-"<< i<<"\n";//For productive purposes
      cout<<(BINpls-i)*time.RealTime()/60<<"m remaining  ###########################  ";
      // cout<<bc<<" "<<v2 <<"\n";
   };
   // for (double i:v2array ){
   //     cout<<bc<<" "<<v2 <<"\n";
   // };
   TGraph *g2= new TGraph(BINpls,binarray,v2array);
   g2-> SetMarkerStyle(22);
   g2->Draw("AP");
   g2->SetTitle("v2 vs pt ;Pt;v2");
};

