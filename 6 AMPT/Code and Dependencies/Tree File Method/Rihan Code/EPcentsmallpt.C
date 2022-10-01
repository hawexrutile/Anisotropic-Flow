#define ptyv2_cxx
#include "ptyv2.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TProfile.h>
#include <algorithm>
#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TMath.h>
#include <iostream>

using namespace std;

void ptyv2::Loop(){
	//   In a ROOT session, you can do:
	//      Root > .L ptyv2.C
	//      Root > ptyv2 t
	//      Root > t.GetEntry(12); // Fill t data members with entry number 12
	//      Root > t.Show();       // Show values of entry 12
	//      Root > t.Show(16);     // Read and show values of entry 16
	//      Root > t.Loop();       // Loop on all entries


	if (fChain == 0) return;
	Long64_t nentries = fChain->GetEntriesFast();
	Long64_t nbytes = 0, nb = 0;
	Int_t nEvTotal = fChain->GetEntries();



	//*Something about high and low centrality
	TFile *f=new TFile("CenRegRoot.root");//includes TH2F *myHist2; This distribution seems appropriately normalized to give actual actual percentage values in the get integral bin; CentralityDist---file name                                                                        
	TH1F *hi=(TH1F*)f->Get("pion");//

	TAxis *axis = hi->GetXaxis();
	int bmin = axis->FindBin(1); //in your case xmin=-1.5                                                                                        
	int bmax = axis->FindBin(620); //in your case xmax=0.8                                                                                       
	double integral = hi->Integral();
	Double_t *Integral;
	Integral=hi->GetIntegral();          //Get integral is the pointer to a aray of bin integrals;
	//  hi->Draw();                                                                                                                              

	const int N=8;                  //!const can be used to define N instead of "DEFINE"

	double lowcentral[N]={0.0};
	double highcentral[N]={0.0};

	int w1=0,w2=0,w3=0,w4=0,w5=0,w6=0,w7=0,w8=0;   //w's are used for locking the if-blocks; once locked the for loop will "continue"  till the next region is reached without setting of any other if-blocks...savy?

	for (int i=2;i<=620;i++){           //*Why is i=2; ig i is the number of bins ; This loop is to set the starting and ending bin of each centrality region(eg 0%-5%)
		if (Integral[i]<=1){			//The lowest if-block gets satisfied by the first bin
			highcentral[0]=i;
		};
		if (Integral[i]>0.95){
			if (w1==0){
				w1=1;
				lowcentral[0]=i;
				highcentral[1]=i;
			};
		};
		if (Integral[i]>0.9){
			if (w2==0){
				w2=1;
				lowcentral[1]=i;
				highcentral[2]=i;
			};
		};
		if (Integral[i]>0.8){
			if (w3==0){
				w3=1;
				lowcentral[2]=i;
				highcentral[3]=i;
			};
		};
		if (Integral[i]>0.7){
			if (w4==0){
				w4=1;
			lowcentral[3]=i;
			highcentral[4]=i;
			};
		};
		if (Integral[i]>0.6){
			if (w5==0){
				w5=1;
				lowcentral[4]=i;
				highcentral[5]=i;
				// lowcentral[1]=i;
				//highcentral[2]=i;
			};
		};
		if (Integral[i]>0.5){
			if (w6==0){
				w6=1;
				lowcentral[5]=i;
				highcentral[6]=i;
				// lowcentral[1]=i;
				//highcentral[2]=i;
			};
		};
		if (Integral[i]>0.4){
			if (w7==0){
				w7=1;
				lowcentral[6]=i;
				highcentral[7]=i;
				//lowcentral[1]=i;
				//highcentral[2]=i;
			};
		};

		if (Integral[i]>0.2){
			if (w8==0){
				w8=1;
				lowcentral[7]=i;
				//lowcentral[2]=i;
			};
		};

	};

	//*For listing the starting and ending bin of each centrality regions.
	for (int i=0;i<N;i++){
		lowcentral[i]=lowcentral[i]*1;
		highcentral[i]=highcentral[i]*1;
		cout<<"low:"<<lowcentral[i]<<"high:"<<highcentral[i]<<endl;
	};



	/////////////////////////////////////////////////////////////////////////


	Double_t pt,pp,eta,phi,Psi2pos,Psi2neg,px,py,pz,En,mass;
	Double_t impact;
	Int_t nTrack, pid, Ncharge;
	Int_t eventCount = 0;
	Double_t azimuthalpos,azimuthalneg,numpos,numneg,denopos,denoneg,azimuthal,v2=0.0,v=0.0;
	const int No=8; //10;

	TFile *fout = new TFile("EPv2vsPTAuAuCentRefSmallPt.root","RECREATE");
	//TFile *fout = new TFile("EPv2vsPTResCheckFineCent.root","RECREATE");


	//Define Histograms/Profiles here:

	Double_t PT[17]={0.2,0.4,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0,2.2,2.4,2.6,2.8,3.0,3.5,4.0};


	Int_t LowCentral[N]={0,5,10,20,30,40,50,60};
	Int_t HighCentral[N]={5,10,20,30,40,50,60,80};
	Double_t Cent[N+1]={0,5,10,20,30,40,50,60,80};
	Double_t Cento[No+1]={9,39,69,113,172,257,370,445,620};//{11,28,44,61,78,103,128,162,195,243,291,355,418,461,503,1252,2000};

	// TProfile *MPTvsCentPos=new TProfile("MeanptvsCentralPos","MeanptvsCentralPos",N,Cent);
	// TProfile *MPTvsCentNeg=new TProfile("MeanptvsCentralNeg","MeanptvsCentralNeg",N,Cent);
	// TProfile *MPTvsMultPos=new TProfile("MeanptvsMultPos","MeanptvsMultPos",No,Cento);
	// TProfile *MPTvsMultNeg=new TProfile("MeanptvsMultNeg","MeanptvsMultNeg",No,Cento);

	// TProfile *MPTvsCentPosAll=new TProfile("MeanptvsCentralPosAll","MeanptvsCentralPosAll",N,Cent);
	// TProfile *MPTvsCentNegAll=new TProfile("MeanptvsCentralNegAll","MeanptvsCentralNegAll",N,Cent);
	// TProfile *MPTvsMultPosAll=new TProfile("MeanptvsMultPosAll","MeanptvsMultPosAll",No,Cento);
	// TProfile *MPTvsMultNegAll=new TProfile("MeanptvsMultNegAll","MeanptvsMultNegAll",No,Cento);

	// MPTvsCentPos->Sumw2();
	// MPTvsCentNeg->Sumw2();
	// MPTvsMultPos->Sumw2();
	// MPTvsMultNeg->Sumw2();

	// MPTvsCentPosAll->Sumw2();
	// MPTvsCentNegAll->Sumw2();
	// MPTvsMultPosAll->Sumw2();
	// MPTvsMultNegAll->Sumw2();

	TProfile *hpRes[N];
	TProfile *hpPionPos[N];
	TProfile *cen=new TProfile("cen","v2",8,0,100);
	// TProfile *hpPionNeg[N];
	// TProfile *hpKaon[N];
	// TProfile *hpProton[N];

	char name[20];
	char name2[20];
	char name3[20];
	char name4[20];
	char name5[20];
	char name6[20];
	char name7[20];
	char name8[20];
	char name9[20];
	char name10[20];
	char name11[20];


	//I think this part is for calculating resolution improve v2 value
	for (int i=0;i<N;i++){
		sprintf(name,"hpRes_cent%d_%d",LowCentral[i],HighCentral[i]);
		hpRes[i]=new TProfile(name,"Resolution",1,-10,10,-2,2);

		sprintf(name2,"hpPionPos_cent%d_%d",LowCentral[i],HighCentral[i]);
		hpPionPos[i]=new TProfile(name2,"PionPos",16,PT);



		hpRes[i]->Sumw2();                  //?why did you call weights
		hpPionPos[i]->Sumw2();
		// hpPionNeg[i]->Sumw2();
		/*hpKaonPos[i]->Sumw2();
		hpKaonNeg[i]->Sumw2();
		hpProtonPos[i]->Sumw2();
		hpProtonNeg[i]->Sumw2();*/
		// hpKaon[i]->Sumw2();
		// hpProton[i]->Sumw2();
	};


	TH1D *heventcount=new TH1D("Events","Events",N,Cent);


	double multi;
	double countev[N]={0.0};
	double midcentral;


	//------------> Event loop starts
		cout<<"loop2-1"<<"\n";
	for(Long64_t jentry=0; jentry<nentries;jentry++){
	//    for(Long64_t jentry=0; jentry<2000;jentry++){
		Long64_t ientry = LoadTree(jentry);
		if (ientry < 0) break;
		nb = fChain->GetEntry(jentry);   nbytes += nb;
		// if (Cut(ientry) < 0) continue;
		impact = Imp;
		//cout<<"value of event loop:"<<jentry<<endl;
		multi=0.0;
		nTrack=Mult;
		for (int i=0;i<nTrack;i++){
			pid = PID[i];
			if(abs(pid)!=211 && abs(pid)!=321 && abs(pid)!=2212) continue;//"Continue" skips to the next iteration
			px = Px[i];
			py = Py[i];
			pz = Pz[i];
			pt = TMath::Sqrt(px*px+py*py);

			if(fabs(pt)==0) continue;

			pp = TMath::Sqrt(px*px+py*py+pz*pz);
			eta = 0.5*TMath::Log((pp+pz)/(pp-pz));
			if ((eta> -0.5 && eta <0.5)){
				multi+=1;
			};
		};

		//cout<<"Value:"<<multi<<endl;


		for (int c=0;c<N;c++){          //For each centrality region
			if (multi>lowcentral[c] && multi <highcentral[c]){
				midcentral=(LowCentral[c]+HighCentral[c])/2.0;
				nTrack = Mult; 
				Ncharge = 0;
				numpos=0.0;//pos and neg is for eta values
				numneg=0.0;//we do this prob cuz to find resolution
				denopos=0.0;
				denoneg=0.0;
				azimuthalpos=0.0;
				azimuthalneg=0.0;

				countev[c]+=1;
				heventcount->SetBinContent((c+1),countev[c]);

				for(int i=0; i<nTrack; i++){                //Loop for calculating Psi
					pid = PID[i];
					if(abs(pid)!=211 && abs(pid)!=321 && abs(pid)!=2212) continue;
					px = Px[i];
					py = Py[i];
					pt = TMath::Sqrt(px*px+py*py);
					if(fabs(pt)==0) continue;
					pz = Pz[i];
					pp = TMath::Sqrt(px*px+py*py+pz*pz);
					eta = 0.5*TMath::Log((pp+pz)/(pp-pz));  //to get a flat dist or something

					if((eta >0.3) && (eta < 0.8)){         //pos eta
						azimuthalpos=TMath::ATan2(py,px);  //ATan2 is like directly gives phi/azimuthal position
						if (azimuthalpos<0.0){
							azimuthalpos+=2*TMath::Pi();
						};
						numpos+=TMath::Sin(2*azimuthalpos);
						denopos+=TMath::Cos(2*azimuthalpos);
					}
					else if((eta > -0.8) && (eta < -0.3)){  //neg eta
						azimuthalneg=TMath::ATan2(py,px);
						if (azimuthalneg<0.0){
							azimuthalneg+=2*TMath::Pi();
						};
						numneg+=TMath::Sin(2*azimuthalneg);
						denoneg+=TMath::Cos(2*azimuthalneg);
					};
					Ncharge++; //do your thing: is used for cout
				};//trackloop #1

				if (denopos==0||denoneg==0) continue;     //To avoid inf values in Psi
				Psi2pos=0.5*(TMath::ATan2(numpos,denopos)); //reactionplaneangle in positive eta
				Psi2neg=0.5*(TMath::ATan2(numneg,denoneg));//reaction plane angle in negative eta
				if (Psi2pos<0.0)Psi2pos+=TMath::Pi();//To account fr trig anomalies
				if (Psi2neg<0.0)Psi2neg+=TMath::Pi();     
				//    cout<<"POS:"<<Psi2pos<<":NEG:"<<Psi2neg<<endl;
				//hpPsiPos[c]->Fill(Psi2pos);
				//hpPsiNeg[c]->Fill(Psi2neg);
				hpRes[c]->Fill(1,TMath::Cos(2*Psi2pos-2*Psi2neg));//*resolution here
				for(int i=0; i<nTrack; i++){ //trackloop #2
					pid = PID[i];
					if(abs(pid)!=211 && abs(pid)!=321 && abs(pid)!=2212) continue;
					px = Px[i];
					py = Py[i];
					pt = TMath::Sqrt(px*px+py*py);
					if(fabs(pt)==0) continue;
					pz = Pz[i];
					pp = TMath::Sqrt(px*px+py*py+pz*pz);
					eta = 0.5*TMath::Log((pp+pz)/(pp-pz));
					azimuthal=TMath::ATan2(py,px);
					if (azimuthal<0.0) azimuthal+=2*TMath::Pi();

					if (pt>0.15 && pt<0.8){
						if((eta >0.3) && (eta <0.8)){
							v2=TMath::Cos(2*azimuthal-2*Psi2neg);
							if (fabs(v2)<1.0){
								if (abs(pid)==2212||abs(pid)==321||abs(pid)==211){
									hpPionPos[c]->Fill(pt,v2);
									// cen->Fill(c,v2);
								};
							};
						}
						else if((eta > -0.8) && (eta < -0.3)){
							v2=TMath::Cos(2*azimuthal-2*Psi2pos);
							if (fabs(v2)<1.0){
								if (abs(pid)==2212||abs(pid)==321||abs(pid)==211){
									hpPionPos[c]->Fill(pt,v2);
									// cen->Fill(c,v2);
								};

							};
						};
					};
				};      //trackloop 2 ends
				eventCount++;
				if(eventCount%50000==0){
					cout<<" Event = "<<eventCount<<" b = "<<impact<<" Nch = "<<Ncharge<<endl;
				};
			};
		};
	};//Event loop
	fout->cd();
	fout->Write();
	fout->Close();
};//main functions Ends
