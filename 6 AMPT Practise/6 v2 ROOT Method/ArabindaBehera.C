#include<TTree.h>
#include<TRandom.h>
#include<iostream>
#include<math.h>
#include<stdio.h>
#include<fstream>
#include<map>
#include<TH2F.h>
#include<TCanvas.h>
#include<TLorentzVector.h>
#include<TROOT.h>

int makeAmptroot_50000(){

    Int_t evn,tn,nov,npp,npt,nesp,ninesp,nest,ninest,pid,counter=0;                     //input file variables:

    Float_t px,py,pz,mass,X,Y,Z,t,b;                                                    // for ampt.dat
    Int_t evn2, it, na, nb, nab, sr, stat, tpp;    
    Float_t nx, ny, zz, theta_p , phi_p , theta_t , phi_t , imp, psi;

    const Int_t mul = 90000;                                                            //Booked variables in tree:
    const Int_t nucl = 600;

    Int_t Event=0, Na, Nb ,Nab, Mult, Npartp, Npartt, Nesp,Ninesp, Nest, Ninest;        // event variables   
    Float_t Imp, Theta_p , Phi_p , Theta_t , Phi_t, Psi;

    Int_t Stat[nucl], PID[mul];                                                         //particle variables
    Float_t Nx[nucl], Ny[nucl], Nz[nucl], Px[mul], Py[mul], Pz[mul];
    Float_t XX[mul], YY[mul], ZZ[mul], TT[mul], Mass[mul];
    //particle variables
    Char_t outfile[100];
    //define a root file:----------------------------------------------
    //sprintf(outfile,"AuAu_200_SM_10mb_folder%d.root",folder);
    TFile *fampt = new TFile("AuAu_200_SM_for_50000_events.root",
    "recreate");
    TTree *tr = new TTree("tr","Reconst ntuple");
    //Define event branches:-------------------------------------------
    tr->Branch("Event", &Event,"Event/I");
    tr->Branch("Mult", &Mult, "Mult/I");
    // multiplicity = tracks
    tr->Branch("Npartp",&Npartp,"Npartp/I");
    tr->Branch("Npartt",&Npartt,"Npartt/I");
    tr->Branch("Nesp", &Nesp, "Nesp/I");
    tr->Branch("Ninesp",&Ninesp, "Ninesp/I");
    tr->Branch("Nest", &Nest, "Nest/I");
    tr->Branch("Ninest",&Ninest, "Ninest/I");
    tr->Branch("Imp",&Imp, "Imp/F");
    tr->Branch("Na", &Na, "Na/I"); //Na = Projectile mass no.
    tr->Branch("Nb", &Nb, "Nb/I"); //Nb = Target mass no.
    tr->Branch("Nab",&Nab, "Nab/I");
    tr->Branch("Psi",&Psi, "Psi/F");
    //particle branches:
    tr->Branch("Nx", &Nx, "Nx[Nab]/F"); // Nab = na+nb;
    tr->Branch("Ny", &Ny, "Ny[Nab]/F");
    tr->Branch("Nz", &Nz, "Nz[Nab]/F");
    tr->Branch("Stat",&Stat,"Stat[Nab]/I");
    tr->Branch("PID", &PID, "PID[Mult]/I");
    tr->Branch("Px", &Px, "Px[Mult]/F");
    tr->Branch("Py", &Py, "Py[Mult]/F");
    tr->Branch("Pz", &Pz, "Pz[Mult]/F");
    tr->Branch("Mass",&Mass,"Mass[Mult]/F");
    tr->Branch("XX", &XX, "XX[Mult]/F");
    tr->Branch("YY", &YY, "YY[Mult]/F");
    tr->Branch("ZZ", &ZZ, "ZZ[Mult]/F");
    tr->Branch("TT", &TT, "TT[Mult]/F");

    //**************************************************************************************
    
    
    cout<<"making .root file from .dat file..."<<endl;
    ifstream infile[3000];
    ifstream infile2[3000];
    // infile2 = npart-xy.dat
    char *fname = new char[100];
    char *fname2 = new char[100];
    for(int fid=2000; fid<2015; fid++){                  //Loop on folder;usually one root file for each folder
        sprintf(fname,"%d/ana/ampt.dat",fid);               //1 folder loop
        sprintf(fname2,"%d/ana/npart-xy.dat",fid);
        infile[fid].open(fname);
        infile2[fid].open(fname2);
        while(infile2[fid])
            {//2 event loop
            //infile[fid]>>evn>>tn>>nov>>b>>npp>>npt>>nesp>>ninesp>>nest>>ninest>>psi;
            infile[fid]>>evn>>tn>>nov>>b>>npp>>npt>>nesp>>ninesp>>nest>>ninest; // for AMPT without psi
            //infile2[fid]>>evn2>>it>>na>>theta_p>>phi_p>>nb>>theta_t>>phi_t>>imp;
            infile2[fid]>>evn2>>it>>na>>nb>>imp;

            if(infile2[fid].eof()) break;
            if(infile[fid].eof()) break;
            //cout<<"No of particles in "<<evn<<" event = "<<nov<<endl;
            Event = Event+1;
            Mult = nov;
            Imp = imp;
            Npartp = npp;
            Npartt = npt; tpp = npp+npt;
            Nesp = nesp;
            Ninesp = ninesp;
            Nest = nest;
            Ninest = ninest ;
            Na = na;
            Nb = nb; nab = na+nb;
            Nab = nab;
            Psi = psi;
            //for U+U only:
            Theta_p = theta_p;
            Phi_p = phi_p;
            Theta_t = theta_t;
            Phi_t = phi_t;
            //****************************ampt.dat particle loop********************************
            for(int j=0;j<nov;j++){                                  //particle loop
                infile[fid]>>pid>>px>>py>>pz>>mass>>X>>Y>>Z>>t;
                Px[j] = px;
                Py[j] = py;
                Pz[j] = pz;
                Mass[j] = mass;
                PID[j] = pid;
                XX[j] = X;
                YY[j] = Y;
                ZZ[j] = Z;
                TT[j] = t;
            }
            //***************************ampt.dat particle loop ends***************************
            //~~~~~~~~~~~ npart-xy.dat loop ~~~~~~~~~~~~~~~~~~~
            for(int k=0;k<nab;k++){
                infile2[fid]>>nx>>ny>>sr>>stat>>zz>>theta_p>>phi_p;
                Nx[k] = nx;
                Ny[k] = ny;
                Nz[k] = zz;
                Stat[k] = stat;
            }
            //~~~~~~~~~~~npart-xy.dat loop ends~~~~~~~~~~~~~~~~~~~~~
            if(evn!=evn2) {
                cout<<"..Error.Event mismatch... evn= "<<evn<<"\txy-ev= "
                <<evn2<<" f_id= "<<fid<<endl;
                break;
            }
            tr->Fill();
            if(Event%100==0){
                cout<<"am= "<<evn<<" xy-ev= "<<evn2<<" f_id= "<<fid<<" Mul= "
                <<Mult<<"\tnpart="<<tpp<<endl;
            }
        }//2 event loop end
    }//1 folder loop end
    fampt->cd();
    tr->Write();
    fampt->Close();
    cout<<" Tree written succesfully "<<endl;
    return 0;
}//main end