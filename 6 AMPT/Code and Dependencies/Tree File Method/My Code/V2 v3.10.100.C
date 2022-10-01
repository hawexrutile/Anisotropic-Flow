//*Code to create centrality-wise v2,v3 vs pt  distributin.
#define check_cxx
#include <algorithm>
#include "check.h"
#include <TH2.h>
#include <TF1.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TLine.h>
#include <TLegend.h>
#include <TMath.h>
#include <iostream>
#include <TStopwatch.h>
#include <TProfile.h>
#include <cmath>
#include <TSystem.h>
//TODO Float_t vs float issue

//?Why are there two similar trees in root file?
using namespace std;             //Global variables have been defined bellow for use in multiple functions
const int MultBin=620;                  //Multiplicity Bin...   Total no of bins in the centrality distribution plot
const int N=10;                   //Total no of Centrality Regions
const int BINpls=10;             //Bins for TProfile
const int sub=0;

Float_t lowcentral[N]={0.0};
Float_t highcentral[N]={0.0};         //Global variable to be used i, both functions
Float_t res;

TFile *CenRegRoot=new TFile("CenRegRoot.root","READ");             //TODO add a if block here to check for file
TH1D *pion= (TH1D*) CenRegRoot->Get("pion");
TProfile *gcen =new TProfile("gcen", "Centrality vs v3",N,0,100);      
TProfile *gres =new TProfile("gres", "Centrality vs res",N,0,100);      


void check::Resolution(int ij){ 
   cout<<"running res, standby \n";
   if (fChain == 0) return;
   Long64_t nentries = fChain->GetEntriesFast();             //total no of events        
   Long64_t nbytes = 0, nb = 0;

   Float_t lowerbc=pion->GetBinCenter(lowcentral[ij]);        //For a given region gives the lowest bin's bin centre(ie the total no of trails)
   Float_t higherbc=pion->GetBinCenter(highcentral[ij]);
   Float_t  AvgNum;
   int totcen=0;
   int TotEvents=nentries-sub;                                             // *TotEvents is for easy adjustament of event size 
   for (Long64_t jentry=0; jentry<TotEvents;jentry++) {                                                                 //  Total events=67295
      Long64_t ientry = LoadTree(jentry);                              //Load Tree apparently gives back the the respective possition in the tree
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;                   //!Loads all the variables(branches) of the ith leaf
      // if (Cut(ientry) < 0) continue;
      Float_t pt;
      Float_t phi;
      Float_t eta;
      Float_t P;
      Double_t azimuthalpos,azimuthalneg,numpos,numneg,denopos,denoneg,azimuthal,Psi2pos,Psi2neg;
      numpos=0.0;//pos and neg is for eta values
      numneg=0.0;//we do this prob cuz to find resolution
      denopos=0.0;
      denoneg=0.0;
      azimuthalpos=0.0;
      azimuthalneg=0.0;
      // cout<<Mult<<"\n";
      // if (jentry%10000==0) cout<< floor((ij*1.0+(jentry*1.0/TotEvents)) *100/(N*1.0)) <<"% complete \n";   //*Percentage finished indicator

      int   chrgi=0;                                     //Charged i entires
      for (int p =0;p < Mult ;p++){
         if(abs(PID[p])!=211 && abs(PID[p])!=321 && abs(PID[p])!=2212) continue;   //We only need pion, proton or kaon otherwise you may skip....
         if ((Px[p]==0) && (Py[p]==0)) continue;
         P=sqrt(pow(Px[p],2)+pow(Py[p],2)+pow(Pz[p],2));
         eta=0.5* log((P+Pz[p])/(P-Pz[p]));
         // if (!((eta> -0.8) && (eta<0.8)) ) continue; 
         if ((eta> -0.5) && (eta< 0.5)) chrgi++ ;                  // 2.8<eta<5.1 or -3.7<eta< -1.7 should be allowed.  //!Not-operator is important
      };
      
      if ((chrgi<lowerbc) || (chrgi>higherbc)) continue;                      //to loop on the given centralities //! It should be charged particle not mult.
      totcen++;
      for (int p =0;p < Mult ;p++){
         if(abs(PID[p])!=211 && abs(PID[p])!=321 && abs(PID[p])!=2212) continue; 
         if ((Px[p]==0) && (Py[p]==0)) continue; 
         P=sqrt(pow(Px[p],2)+pow(Py[p],2)+pow(Pz[p],2));
         eta=0.5* log((P+Pz[p])/(P-Pz[p]));
         if((eta >=0) && (eta < 1)){         //pos eta  //if((eta >0.3) && (eta < 0.8)){  //
            azimuthalpos=TMath::ATan2(Py[p],Px[p]);  //ATan2 is like directly gives phi/azimuthal position
            if (azimuthalpos<0.0) azimuthalpos+=2*TMath::Pi();
            numpos+=TMath::Sin(2*azimuthalpos);
            denopos+=TMath::Cos(2*azimuthalpos);
         }

         else if((eta > -1) && (eta < 0)){  //neg eta  //else if((eta > -0.8) && (eta < -0.3)){  
            azimuthalneg=TMath::ATan2(Py[p],Px[p]);
            if (azimuthalneg<0.0) azimuthalneg+=2*TMath::Pi();
            numneg+=TMath::Sin(2*azimuthalneg);
            denoneg+=TMath::Cos(2*azimuthalneg);
         };
      };
      if (denopos==0||denoneg==0) continue;       //To avoid inf values in Psi
      Psi2pos=0.5*(TMath::ATan2(numpos,denopos)); //reactionplaneangle in positive eta
      Psi2neg=0.5*(TMath::ATan2(numneg,denoneg)); //reaction plane angle in negative eta
      if (Psi2pos<0.0)Psi2pos+=TMath::Pi();       //To account fr trig anomalies
      if (Psi2neg<0.0)Psi2neg+=TMath::Pi();     
      AvgNum+=TMath::Cos((2*Psi2pos)-(2*Psi2neg));//?resolution here                        //shopul it be 1/res or res
      // if (isnan(TMath::Cos(2*Psi2pos-2*Psi2neg))) cout<<2*Psi2pos-2*Psi2neg<<"\n";
   };
   res=AvgNum/totcen;
   // gres->Fill((ij*10)+5,res);

};


void check::Loop(int ij,TCanvas *canva ){           //SECTION //*Main loop function to calcualte and plot v2



   if (fChain == 0) return;
   Long64_t nentries = fChain->GetEntriesFast();             //total no of events        
   Long64_t nbytes = 0, nb = 0;
   

   Float_t xvalue[BINpls];
   Float_t xvaluer[BINpls];
   Float_t yvalue2[BINpls];
   Float_t yvaluer2[BINpls];
   Float_t yvalue3[BINpls];
   Float_t yvaluer3[BINpls];
   this->Resolution(ij);
   cout<<"res "<<res<<"\n";

   TProfile *e2 =new TProfile("e2", "v2 vs Pt  ",BINpls,0,5.0);   //!dont specidfy y range here nstead here LINK 6 AMPT Practise\6 v2 ROOT Method\check.C:79
   TProfile *e3 =new TProfile("e3", "v3 vs Pt",BINpls,0,5.0);    //


   Float_t lowerbc=pion->GetBinCenter(lowcentral[ij]);                                           //For a given region gives the lowest bin's bin centre(ie the total no of trails)
   cout<<"lowerbc "<< lowerbc<<"\n";
   Float_t higherbc=pion->GetBinCenter(highcentral[ij]);
   cout<<"higherbc "<< higherbc<<"\n";

   int TotEvents=nentries-sub;                                             // *TotEvents is for easy adjustament of event size 
   for (Long64_t jentry=0; jentry<TotEvents;jentry++) {                 //Eventloop
//                                                                       Total events=67295
      Long64_t ientry = LoadTree(jentry);                              //Load Tree apparently gives back the the respective possition in the tree
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;                   //!Loads all the variables(branches) of the ith leaf
      // if (Cut(ientry) < 0) continue;
      Float_t pt;
      Float_t phi;
      Float_t eta;
      Float_t P;
      // cout<<Mult<<"\n";
      if (jentry%10000==0) cout<< floor((ij*1.0+(jentry*1.0/TotEvents)) *100/(N*1.0)) <<"% complete \n";   //*Percentage finished indicator

      int   chrgi=0;                                     //Charged i entires
      for (int p =0;p < Mult ;p++){                       //*No of charged particle counter
         if(abs(PID[p])!=211 && abs(PID[p])!=321 && abs(PID[p])!=2212) continue;   //We only need pion, proton or kaon otherwise you may skip....
         if ((Px[p]==0) && (Py[p]==0)) continue;
         P=sqrt(pow(Px[p],2)+pow(Py[p],2)+pow(Pz[p],2));
         eta=0.5* log((P+Pz[p])/(P-Pz[p]));
         // if (!((eta> -0.8) && (eta<0.8)) ) continue; 
         if ((eta> -0.5) && (eta< 0.5) )   chrgi++ ;          // 2.8<eta<5.1 or -3.7<eta< -1.7 should be allowed.  //!Not-operator is important
      };

      if ((chrgi<lowerbc) || (chrgi>higherbc)) continue;                      //to loop on the given centralities //! It should be charged particle not mult.
      Float_t SinSum=0; 
      Float_t CosSum=0;
      for (int p =0;p < Mult ;p++){   //*calculate Psi2.
         if(abs(PID[p])!=211 && abs(PID[p])!=321 && abs(PID[p])!=2212 ) continue;//We only need pion, proton or kaon otherwise you may skip....
         if ((Px[p]==0) && (Py[p]==0)) continue;
         P=sqrt(pow(Px[p],2)+pow(Py[p],2)+pow(Pz[p],2));
         eta=0.5* log((P+Pz[p])/(P-Pz[p]));
         if (!((eta>= -1) && (eta< 0.0)) ) continue;       
         phi=(TMath::ATan2(Py[p],Px[p]));
         if (phi<0.0)phi+=2*TMath::Pi();
         SinSum+= sin(2*phi);        //summation of cos and sin phi
         CosSum+= cos(2*phi);        //?should Psi2 be for all the trails in the event or for the particular particle or for that particular ptbi
      };
      if (CosSum==0) continue;
      Psi=0.5*(TMath::ATan2(SinSum,CosSum));         //for each event we calcaulate Psi2
      if (Psi<0.0)  Psi+=2*TMath::Pi();                  //? SHould i normalise v0 or something
      for (int p =0;p < Mult ;p++){
         if(abs(PID[p])!=211 && abs(PID[p])!=321 && abs(PID[p])!=2212 ) continue;//We only need pion, proton or kaon otherwise you may skip....
         if ((Px[p]==0) && (Py[p]==0)) continue;        //We need  0 pt so && should suffice...It diddnt make a difference when I changed from "||" to "&&" so it should be safe to change;
         P=sqrt(pow(Px[p],2)+pow(Py[p],2)+pow(Pz[p],2));
         eta=0.5* log((P+Pz[p])/(P-Pz[p]));
         if (!((eta>= -1.0) && (eta< 0.0)) ) continue;  //if (!((eta> 0) && (eta< 0.8)) ) continue;                         //CONFINING WITH ETA
         pt=sqrt(pow(Px[p],2)+pow(Py[p],2));  //for that particular event        //Pt has units of GeV/c^2
         phi=(TMath::ATan2(Py[p],Px[p]));
         if (phi<0.0)phi+=2*TMath::Pi();
         if (fabs(cos(2*(phi-Psi))/res)>=1) continue;        //?why is it equal..?
         e2->Fill(pt,(cos(2*(phi-Psi)))/res);                             //*v2
         // e3->Fill(pt,cos(3*(phi-Psi)));                             //*v3     //TODO #1 Fix psi for v3
         if (pt>0.15 && pt<0.8)   gcen->Fill((ij*10)+5,(cos(2*(phi-Psi)))/res);
      };
   };
   cout<< gcen->GetBinContent(ij)<<"\t"<<ij<<"th gcen bin \n";


   e2->GetYaxis()->SetRangeUser(-1.0, 1.0);       
   e3->GetYaxis()->SetRangeUser(-1.0, 1.0);                                //!Set Y range here
   for (int ji=0;ji<N;ji++){
      xvalue[ji]=e2->GetBinCenter(ji);
      xvaluer[ji]=0;
      yvalue2[ji]=e2->GetBinContent(ji);
      yvaluer2[ji]=e2->GetBinError(ji);
      yvalue3[ji]=e3->GetBinContent(ji);
      yvaluer3[ji]=e3->GetBinError(ji);
   };

   TGraphErrors *g2 =new TGraphErrors(BINpls,xvalue,yvalue2,xvaluer,yvaluer2);   //!dont specidfy y range here nstead here LINK 6 AMPT Practise\6 v2 ROOT Method\check.C:79
   TGraphErrors *g3 =new TGraphErrors(BINpls,xvalue,yvalue3,xvaluer,yvaluer3);    
   // TF1 *fit2= new TF1("fit","pol2",0,4.0,"Q");            //"Q" is to supress terminal output
   // TF1 *fit3= new TF1("fit","pol2",0,4.0,"Q");
   // fit2->SetParNames("a","b","c");
   // fit2->SetLineWidth(3);               //fit-line thickness
   // fit2->SetLineColor(kBlue+3);               //fit-line thickness
   // fit2->SetLineStyle(1);               //1-solid ; 2- dotted
   // fit3->SetParNames("a","b","c");
   // fit3->SetLineWidth(3);               //fit-line thickness
   // fit3->SetLineColor(kMagenta+3);               //fit-line thickness
   // fit3->SetLineStyle(1);               //1-solid ; 2- dotted

   g2->GetXaxis()->SetRangeUser(0, 4.0); 
   g2->GetYaxis()->SetRangeUser(-0.2, 0.4);   
   g2->SetMarkerStyle(47);       
   g2->SetMarkerSize(3);      
   g2->Draw("AP ");
   g2->SetTitle("Pt vs v2; pt; v2");
   // g2->Fit(fit2,"R");
   g2->SetTitle(0);

   // g3->GetXaxis()->SetRangeUser(0, 5.0); 
   // g3->GetYaxis()->SetRangeUser(-1.0, 1.0);                                //!Set Y range here
   // g3->SetMarkerStyle(46);      
   // g3->SetMarkerSize(3);      
   // g2->SetTitle("Pt vs v2 and v3; pt; v2,v3");
   // g3->Draw("P SAME");
   // g3->Fit(fit3,"R");
   // g3->SetTitle(0);


   TLine *linex=new TLine(0,0,4.0,0);                    //Cordinate (x1,y1,x2,y2)
   linex->SetLineWidth(5);
   linex->SetLineColor(kBlack);
   linex->Draw();
   TLine *liney=new TLine(0,-0.2,0,0.4);                    //Cordinate (x1,y1,x2,y2)
   liney->SetLineWidth(5);
   liney->SetLineColor(kBlack);
   liney->Draw();

   TLegend *leg=new TLegend(0.55,0.1,0.85,0.2);         //Consider the canvas as an cordinate region withe the top right corner as 1,1 and bottom left corner as 0,0;Cordinate (x1,y1,x2,y2)  
   leg->SetBorderSize(0);                                 //sets the margin to 0
   leg->AddEntry(g2,"v2","p");               //L IS FOR LINE; F OR FILLING ; P FOR POINT
   // leg->AddEntry(g3,"v3","p");               //L IS FOR LINE; F OR FILLING ; P FOR POINT
   // leg->AddEntry(fit2,"v2 polinomial fit","l");
   // leg->AddEntry(fit3,"v3 polinomial fit","l");
   leg->Draw();

};

void check::Multiplicity(){                                //SECTION Creates a multiplicity distribution and saves the hostpgram on a root file
   if (fChain == 0) return;
   Long64_t nentries = fChain->GetEntriesFast();
   Long64_t nbytes = 0, nb = 0;
   int TotEvents=nentries;   //Total events=67295           //TotEvents has been defined, to make it easy to  change no of events, for quick bugfixing
   int min=500;
   int max=0;
   
   TFile *CenRegRoot=new TFile("CenRegRoot.root","recreate");
   TH1D *pion= new TH1D("pion","Pion Multiplicity",MultBin,0,700);   // !int maxTrack=395; int minTrack=2000 after eta confinemnt. it was 5000 bfr; 395 min is imp
//######################################################################
   for (Long64_t jentry=0; jentry<TotEvents;jentry++) {
      Long64_t ientry = LoadTree(jentry);                  //Load Tree apparently gives back the the respective possition in the tree
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      int   chrgi=0;    
      Float_t eta;                                 //Charged i entires
      Float_t P;
      for (int p =0;p < Mult ;p++){
         if(abs(PID[p])!=211 && abs(PID[p])!=321 && abs(PID[p])!=2212) continue;   //We only need pion, proton or kaon otherwise you may skip....
         if (Px[p]==0 && Py[p]==0) continue;
         P=sqrt(pow(Px[p],2)+pow(Py[p],2)+pow(Pz[p],2));
         eta=0.5 * log((P+Pz[p])/(P-Pz[p]));
         if (((eta> -0.5) && (eta< 0.5)) ) chrgi++;     // 2.8<eta<5.1 or -3.7<eta< -1.7 should be allowed.  //!Not-operator is important
      };
      if (jentry%5000==0) cout<< floor(((jentry*1.0)/TotEvents) *100) <<"% complete \n";
      // if (chrgi>max) max=chrgi;
      // if (chrgi<min) min=chrgi;
      pion->Fill(chrgi);//For counting charged particles(the big 3s)
      
//########################################################################
   };
      cout<<min<<"\t"<<max<<"\n";                       //;WHY IS NEG V2 EXACTLY AT 6; WHT GCEN HAS 0 ERROR;WHAT HAPPENS IF I REMOVE NEG V2;WHY IS V2 NEG IN THE FIRST PLACE
   pion->Draw();
   CenRegRoot->Write();
   CenRegRoot->Close();
};

// void check::Centralizer(){
//    map <string,int> locks;

// }
void check::Printer(){                                           //SECTION uses the multiplicity distribution and runs the 4-loop on Loop() function for all centralities...

   double_t *CenReg;
   CenReg= pion->GetIntegral();              //!cumulative sum of bins ...The resulting integral is normalized to 1

  
   int w1=0,w2=0,w3=0,w4=0,w5=0,w6=0,w7=0,w8=0,w9=0,w10=0;       //Locks
   for (int i=0;i<=MultBin;i++){           //?Why is i=2; cuz in protay's code initial bins have to be skiped ; This loop is to set the starting and ending bin of each centrality region(eg 0%-5%)
                                          //!The lowest if-block gets satisfied by the first bin
      if (CenReg[i]<=1){		  //keeps on changing till the last bin
			highcentral[0]=i;      //Highest bin  of the highest centrality region(0-20%)
		};
		if (CenReg[i]>0.9){  
			if (w1==0){            //locks last  //TODO fix the comments
				w1=1;
				lowcentral[0]=i;    //Lowest Bin of the highest centrality region(0-20%)
				highcentral[1]=i;   //Highest bin  of the sevcond highest centrality region(20-40%)
            // cout<<1;
			};
		};
		if (CenReg[i]>0.8){
			if (w2==0){
				w2=1;
				lowcentral[1]=i;
				highcentral[2]=i;
            // cout<<2;
			};
		};
		if (CenReg[i]>0.7){
			if (w3==0){
				w3=1;
				lowcentral[2]=i;
				highcentral[3]=i;
            // cout<<3;
			};
		};
		if (CenReg[i]>0.6){
			if (w4==0){                    //locks second
				w4=1;
			   lowcentral[3]=i;               //Lowest bin  of the second  lowest centrality region(60-80%) 
			   highcentral[4]=i;              //Highest bin  of the lowest centrality region(80-100%) 
            // cout<<4;
			};
		};
		if (CenReg[i]>0.5){
			if (w5==0){                    //locks second
				w5=1;
			   lowcentral[4]=i;               //Lowest bin  of the second  lowest centrality region(60-80%) 
			   highcentral[5]=i;              //Highest bin  of the lowest centrality region(80-100%) 
            // cout<<4;
			};
		};
		if (CenReg[i]>0.4){
			if (w6==0){                    //locks second
				w6=1;
			   lowcentral[5]=i;               //Lowest bin  of the second  lowest centrality region(60-80%) 
			   highcentral[6]=i;              //Highest bin  of the lowest centrality region(80-100%) 
            // cout<<4;
			};
		};
		if (CenReg[i]>0.3){
			if (w7==0){                    //locks second
				w7=1;
			   lowcentral[6]=i;               //Lowest bin  of the second  lowest centrality region(60-80%) 
			   highcentral[7]=i;              //Highest bin  of the lowest centrality region(80-100%) 
            // cout<<4;
			};
		};
		if (CenReg[i]>0.2){
			if (w8==0){                    //locks second
				w8=1;
			   lowcentral[7]=i;               //Lowest bin  of the second  lowest centrality region(60-80%) 
			   highcentral[8]=i;              //Highest bin  of the lowest centrality region(80-100%) 
            // cout<<4;
			};
		};
		if (CenReg[i]>0.1){
			if (w9==0){                    //locks second
				w9=1;
			   lowcentral[8]=i;               //Lowest bin  of the second  lowest centrality region(60-80%) 
			   highcentral[9]=i;              //Highest bin  of the lowest centrality region(80-100%) 
            // cout<<4;
			};
		};
		if (CenReg[i]>=0){                //Locks first   
			if (w10==0){
				w10=1;
				lowcentral[9]=i;            //Lowest bin  of the lowest centrality region(80-100%) 
            // cout<<5;
			};
		};
	};
   for (int i=0;i<N;i++){
		cout<<"low:"<<lowcentral[i]<<"\t"<<"high:"<<highcentral[i]<<endl;
	};
// //##############################################################################
   for (int i=0;i<N;i++){                                             //*For ploting all the plots  
      if (i==8 || i==9) continue;
      TCanvas *canva = new TCanvas("canva","Pt vs v2 and v3 for Proton, Pion and Kaon",1500,1000); //I made a Tcanvas cuz to use the print function  //dontuse TCanvas::Divide its harder to bring grid
      gPad->SetTickx();                     //covers the axis above the x axis with "scale"
      gPad->SetTicky();                     //covers the axis beside the y axis with "scale"  
      gPad->SetGridy();
      gPad->SetGridx();
      // canva->SetTickx();                     //covers the axis above the x axis with "scale"
      // canva->SetTicky();                     //covers the axis beside the y axis with "scale"  
      // canva->SetGridy();
      // canva->SetGridx();
      this->Loop(i,canva);                                          //For calling functions within the same class
      canva->Print(Form("./Plot/del%d.png",(i*10)+5));           //TODO remove memory leak warning
      if (canva){ 
         canva->Close(); 
         gSystem->ProcessEvents(); 
         delete canva; 
         canva = 0; 
      };
   };
//#######################################################################################
   TCanvas *canvacen = new TCanvas("canvacen","v2 vs Centrality for  Proton, Pion and Kaon",1500,1000);
   TF1 *fitcen= new TF1("fitcen","pol2",0,100);
   fitcen->SetLineWidth(5);                   //fit-line thickness
   fitcen->SetLineColor(kBlue);                   //fit-line thickness
   fitcen->SetLineStyle(3);
   gcen->SetStats(0);  
   canvacen->SetTickx();                     //covers the axis above the x axis with "scale"
   canvacen->SetTicky();                     //covers the axis beside the y axis with "scale"  
   canvacen->SetGridy();
   canvacen->SetGridx();
   gcen->GetYaxis()->SetRangeUser(0, 0.2);   
   gcen-> SetMarkerStyle(47);
   gcen->Fit("fitcen","R");  
   gcen-> SetMarkerSize(3); 
   gcen->SetTitle("v2 vs Centrality; Centrality; v2");
   gcen->Draw();
   TLegend *legcen=new TLegend(0.55,0.65,0.85,0.85);            //Consider the canvas as an cordinate region withe the top right corner as 1,1 and bottom left corner as 0,0;Cordinate (x1,y1,x2,y2)  
   legcen->SetBorderSize(0);                                 //sets the margin to 0
   legcen->AddEntry(gcen,"v2","p");               //L IS FOR LINE; F OR FILLING ; P FOR POINT
   legcen->AddEntry(fitcen,"v2 polinomial fit","l");
   legcen->Draw();
   canvacen->Print("./Plot/delcen.png"); 
   TCanvas *canvares = new TCanvas("canvares","v2 vs Centrality for  Proton, Pion and Kaon",1500,1000);
   gres->Draw();

};

void check::hist(){
   if (fChain == 0) return;
   Long64_t nentries = fChain->GetEntriesFast();
   Long64_t nbytes = 0, nb = 0;
   int TotEvents=nentries;   //Total events=67295           //TotEvents has been defined, to make it easy to  change no of events, for quick bugfixing
   TCanvas *canva = new TCanvas("canva","canva",2000,1000);
   canva->SetTickx();                     //covers the axis above the x axis with "scale"
   canva->SetTicky();                     //covers the axis beside the y axis with "scale"  
   canva->SetGridx();
   canva->SetGridy();
   TH1D *histeta = new TH1D("histeta","Pseudorapidity Distribution",1000,-8,8);  //TH(histogram)1(Dimension)D(Doubles) Name Title NoofBins Begin and end //*No of bins dhould be square root of the number of entries
   TH1D *histy = new TH1D("histy","Rapidity Distribution",1000,-8,8);  //TH(histogram)1(Dimension)D(Doubles) Name Title NoofBins Begin and end //*No of bins dhould be square root of the number of entries
   TH1D *histphi = new TH1D("histphi","Polar Angle Distribution",1000,-1.6,1.6);  //TH(histogram)1(Dimension)D(Doubles) Name Title NoofBins Begin and end //*No of bins dhould be square root of the number of entries
   TH1D *histpt = new TH1D("histpt","Transverse Momentum Distribution",1000,0,4);  //TH(histogram)1(Dimension)D(Doubles) Name Title NoofBins Begin and end //*No of bins dhould be square root of the number of entries
   
//######################################################################
   for (Long64_t jentry=0; jentry<TotEvents;jentry++) {
      Long64_t ientry = LoadTree(jentry);                  //Load Tree apparently gives back the the respective possition in the tree
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      Float_t eta;                                 //Charged i entires
      Float_t P;
      Float_t phi; 
      Float_t y;
      Float_t pt; 
      for (int p =0;p < Mult ;p++){
         if(abs(PID[p])!=211 && abs(PID[p])!=321 && abs(PID[p])!=2212) continue;   //We only need pion, proton or kaon otherwise you may skip....
         if (Px[p]==0 && Py[p]==0) continue;
         P=sqrt(pow(Px[p],2)+pow(Py[p],2)+pow(Pz[p],2));
         pt=sqrt(pow(Px[p],2)+pow(Py[p],2));
         eta=0.5 * log((P+Pz[p])/(P-Pz[p]));
         phi=atan(Py[p]/Px[p]);
         y=0.5 * log(abs((sqrt(Px[p]*Px[p]+Py[p]*Py[p]+Pz[p]*Pz[p]+Mass[p]*Mass[p])+Pz[p])/(sqrt(Px[p]*Px[p]+Py[p]*Py[p]+Pz[p]*Pz[p]+Mass[p]*Mass[p])-Pz[p])));
         histeta->Fill(eta);
         histphi->Fill(phi);
         histy->Fill(y);
         histpt->Fill(pt);
      };
      if (jentry%5000==0) cout<< floor(((jentry*1.0)/TotEvents) *100) <<"% complete \n";
      
//########################################################################
   };


   histeta->SetStats(0);                   //*for removing the stat box on the top right corner
   histeta->SetMinimum(0);
   histeta->SetFillColor(kBlue-9);
   histeta->SetTitle("Pseudoapidity Distribution;Pseudoapidity;Count"); 
   histeta->SetLineColor(kRed); 
   histeta->SetTitleSize(0.05,"X"); 
   histeta->SetTitleSize(0.05,"Y");    
   histeta->Draw();
   canva->Print("./Plot/eta.png");
   
   histy->SetStats(0);                   //*for removing the stat box on the top right corner
   histy->SetMinimum(0);
   histy->SetFillColor(kBlue-9);
   histy->SetLineColor(kRed); 
   histy->SetTitle("Rapidity Distribution;Rapidity;Count"); 
   histy->SetTitleSize(0.05,"X"); 
   histy->SetTitleSize(0.05,"Y");    
   histy->Draw();
   canva->Print("./Plot/y.png");
   
   histphi->SetStats(0);                   //*for removing th e stat box on the top right corner
   histphi->SetMinimum(0);
   histphi->SetFillColor(kBlue-9);
   histphi->SetLineColor(kRed); 
   histphi->SetTitle("Polar Angle Distribution;Polar Angle;Count"); 
   histphi->SetTitleSize(0.05,"X"); 
   histphi->SetTitleSize(0.05,"Y");    
   histphi->Draw();
   canva->Print("./Plot/phi.png");
   
   histpt->SetStats(0);                   //*for removing th e stat box on the top right corner
   histpt->SetMinimum(0);
   histpt->SetFillColor(kBlue-9);
   histpt->SetLineColor(kRed); 
   histpt->SetTitle("Transverse Momentum Distribution;Pt;Count"); 
   histpt->SetTitleSize(0.05,"X"); 
   histpt->SetTitleSize(0.05,"Y");    
   histpt->Draw();
   canva->Print("./Plot/pt.png");
   

};
