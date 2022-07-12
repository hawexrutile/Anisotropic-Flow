#define check_cxx
#define BINpls 100
#include <algorithm>
#include "check.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TGraph.h>

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

   

   // int BINpls=100;                                                //set BIN value here; the number of points you want to plot
   double bw=5/(1.0*BINpls);                               //I subbed (Total width of plot) with 5 cuz it was 5 for my prev plot and other plots I saw online
   int i=0;
   double Bmin=0;
   double Bmax=0;

   while (i<BINpls){               //for each pt-BIN we run the loop on all the trails
      double v2=0;
      double sum=0;
      int j=0;
      int lpos1=0;             //looping int fot calculationg Psi
      int lpos2=0;             //looping int for calculating tot
      int Tot=0;               //Variable to calculate v2
      Bmin=i*(bw);
      Bmax=(i+1)*(bw);
      double bc=((2*i+1)/2)*bw;    //BIN centre
      Long64_t nbytes = 0, nb = 0;
      for (Long64_t jentry=0; jentry<nentries;jentry++) {
         Long64_t ientry = LoadTree(jentry);                  //Load Tree apparently gives back the the respective possition in the tree
         if (ientry < 0) break;
         nb = fChain->GetEntry(jentry);   nbytes += nb;
         // if (Cut(ientry) < 0) continue;
         double pt;
         double phi;
            // double SinSum=0;                                         //?Is the Psi defined by check.h actualy event plane?
            // double CosSum=0;
            // while ((lpos1>=listtrail[j]) && (lpos1<listtrail[j+1])){   //This while loop is to calculate Psi. Here listtrail is a list of no events in each trail
            //          SinSum=SinSum +sin(2*listphi[lpos1]);        //summation of cos and sin phi
            //          CosSum=CosSum +cos(2*listphi[lpos1]);        //?should Psi be for all the trails in the event or for the particular particle or for that particular ptbin?
            //          lpos1++;

            // };
            // double Psi=(1/2)*atan(SinSum/CosSum);         //for each event we calcaulate Psi

         for (int p =0;p < sizeof(Px)/sizeof(float_t);p++){
            pt=sqrt(pow(Px[p],2)+pow(Py[p],2));  //for that particular event
            phi=(atan(Pz[p]/pt));
               if ((pt>=Bmin)&&(pt<Bmax && phi>=0)){     //if the ith particle's pt belongs to the BIN    //if phi is greater than 0 ; cuz we are averaging from 0-2pi
                     sum=sum+cos(phi-Psi);               
                     Tot++;
               };
         };
      };

      v2=sum/Tot;
      binarray[i]=bc;           //These arrays are used for plotting
      v2array[i]=v2;
      i++;
   };
   // for (double i:v2array ){
   //     cout<<i;
   // };
   TGraph *g2= new TGraph(BINpls,binarray,v2array);
   g2-> SetMarkerStyle(22);
   g2->Draw("AP");
   g2->SetTitle("v2 vs pt ;Pt;v2");
};

