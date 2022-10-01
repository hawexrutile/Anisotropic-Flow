//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Thu Sep  5 10:44:15 2019 by ROOT version 5.34/30
// from TTree tr/Reconst ntuple
// found on file: AuAu200GeV_SM3mb_Alpha0p47_NT150_DecayOff_0_20fm_file1.root
//////////////////////////////////////////////////////////

#ifndef ptyv2_h
#define ptyv2_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <fstream>

// Header file for the classes stored in the TTree if any.

// Fixed size dimensions of array or collections stored in the TTree if any.

class ptyv2 {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   Int_t           Event;
   Int_t           Mult;
   Int_t           Npartp;
   Int_t           Npartt;
   Int_t           Nesp;
   Int_t           Ninesp;
   Int_t           Nest;
   Int_t           Ninest;
   Float_t         Imp;
   Int_t           Na;
   Int_t           Nb;
   Int_t           Nab;
   Float_t         Psi;
   Float_t         Theta_p;
   Float_t         Phi_p;
   Float_t         Theta_t;
   Float_t         Phi_t;
   Float_t         Nx[394];   //[Nab]
   Float_t         Ny[394];   //[Nab]
   Float_t         Nz[394];   //[Nab]
   Int_t         Stat[394];   //[Nab]
   Int_t          PID[10000];   //[Mult]
   Float_t         Px[10000];   //[Mult]
   Float_t         Py[10000];   //[Mult]
   Float_t         Pz[10000];   //[Mult]
   Float_t       Mass[10000];   //[Mult]
   Float_t         XX[10000];   //[Mult]
   Float_t         YY[10000];   //[Mult]
   Float_t         ZZ[10000];   //[Mult]
   Float_t         TT[10000];   //[Mult]
   Int_t           qMult;
   Int_t          qPID[12916];   //[qMult]
   Float_t         qPx[12916];   //[qMult]
   Float_t         qPy[12916];   //[qMult]
   Float_t         qPz[12916];   //[qMult]
   Float_t       qMass[12916];   //[qMult]
   Float_t         qXX[12916];   //[qMult]
   Float_t         qYY[12916];   //[qMult]
   Float_t         qZZ[12916];   //[qMult]
   Float_t         qTT[12916];   //[qMult]

   // List of branches
   TBranch        *b_Event;   //!
   TBranch        *b_Mult;   //!
   TBranch        *b_Npartp;   //!
   TBranch        *b_Npartt;   //!
   TBranch        *b_Nesp;   //!
   TBranch        *b_Ninesp;   //!
   TBranch        *b_Nest;   //!
   TBranch        *b_Ninest;   //!
   TBranch        *b_Imp;   //!
   TBranch        *b_Na;   //!
   TBranch        *b_Nb;   //!
   TBranch        *b_Nab;   //!
   TBranch        *b_Psi;   //!
   TBranch        *b_Theta_p;   //!
   TBranch        *b_Phi_p;   //!
   TBranch        *b_Theta_t;   //!
   TBranch        *b_Phi_t;   //!
   TBranch        *b_Nx;   //!
   TBranch        *b_Ny;   //!
   TBranch        *b_Nz;   //!
   TBranch        *b_Stat;   //!
   TBranch        *b_PID;   //!
   TBranch        *b_Px;   //!
   TBranch        *b_Py;   //!
   TBranch        *b_Pz;   //!
   TBranch        *b_Mass;   //!
   TBranch        *b_XX;   //!
   TBranch        *b_YY;   //!
   TBranch        *b_ZZ;   //!
   TBranch        *b_TT;   //!
   TBranch        *b_qMult;   //!
   TBranch        *b_qPID;   //!
   TBranch        *b_qPx;   //!
   TBranch        *b_qPy;   //!
   TBranch        *b_qPz;   //!
   TBranch        *b_qMass;   //!
   TBranch        *b_qXX;   //!
   TBranch        *b_qYY;   //!
   TBranch        *b_qZZ;   //!
   TBranch        *b_qTT;   //!

   ptyv2(TTree *tree=0);
   virtual ~ptyv2();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef ptyv2_cxx
ptyv2::ptyv2(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.


   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("AuAu200GeV_SM3mb_Alpha0p47_NT150_DecayOff_0_20fm_file10.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("AuAu200GeV_SM3mb_Alpha0p47_NT150_DecayOff_0_20fm_file10.root");
      }
      f->GetObject("tr",tree);

   }
   Init(tree);
   }


   // if(tree == 0){
   // TChain* chain = new TChain("tr");
   //  //    ifstream FILELIST("inputfiles.list");
   //    ifstream FILELIST("semfinal.list");                                    
   //  //ifstream FILELIST("prottay.list");                                    
   //  //ifstream FILELIST(inputfile);                                         
   //  //assert(FILELIST);                                                         
   // string filename;
   // while (getline(FILELIST,filename)){
   //    chain->Add(filename.c_str());
   //    cout<<" reading input file = "<<filename<<endl;
   // }
   // FILELIST.close();
   // tree = chain;
   // }
   // Init(tree);

// }


ptyv2::~ptyv2()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t ptyv2::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t ptyv2::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void ptyv2::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("Event", &Event, &b_Event);
   fChain->SetBranchAddress("Mult", &Mult, &b_Mult);
   fChain->SetBranchAddress("Npartp", &Npartp, &b_Npartp);
   fChain->SetBranchAddress("Npartt", &Npartt, &b_Npartt);
   fChain->SetBranchAddress("Nesp", &Nesp, &b_Nesp);
   fChain->SetBranchAddress("Ninesp", &Ninesp, &b_Ninesp);
   fChain->SetBranchAddress("Nest", &Nest, &b_Nest);
   fChain->SetBranchAddress("Ninest", &Ninest, &b_Ninest);
   fChain->SetBranchAddress("Imp", &Imp, &b_Imp);
   fChain->SetBranchAddress("Na", &Na, &b_Na);
   fChain->SetBranchAddress("Nb", &Nb, &b_Nb);
   fChain->SetBranchAddress("Nab", &Nab, &b_Nab);
   fChain->SetBranchAddress("Psi", &Psi, &b_Psi);
   fChain->SetBranchAddress("Theta_p", &Theta_p, &b_Theta_p);
   fChain->SetBranchAddress("Phi_p", &Phi_p, &b_Phi_p);
   fChain->SetBranchAddress("Theta_t", &Theta_t, &b_Theta_t);
   fChain->SetBranchAddress("Phi_t", &Phi_t, &b_Phi_t);
   fChain->SetBranchAddress("Nx", Nx, &b_Nx);
   fChain->SetBranchAddress("Ny", Ny, &b_Ny);
   fChain->SetBranchAddress("Nz", Nz, &b_Nz);
   fChain->SetBranchAddress("Stat", Stat, &b_Stat);
   fChain->SetBranchAddress("PID", PID, &b_PID);
   fChain->SetBranchAddress("Px", Px, &b_Px);
   fChain->SetBranchAddress("Py", Py, &b_Py);
   fChain->SetBranchAddress("Pz", Pz, &b_Pz);
   fChain->SetBranchAddress("Mass", Mass, &b_Mass);
   fChain->SetBranchAddress("XX", XX, &b_XX);
   fChain->SetBranchAddress("YY", YY, &b_YY);
   fChain->SetBranchAddress("ZZ", ZZ, &b_ZZ);
   fChain->SetBranchAddress("TT", TT, &b_TT);
   fChain->SetBranchAddress("qMult", &qMult, &b_qMult);
   fChain->SetBranchAddress("qPID", qPID, &b_qPID);
   fChain->SetBranchAddress("qPx", qPx, &b_qPx);
   fChain->SetBranchAddress("qPy", qPy, &b_qPy);
   fChain->SetBranchAddress("qPz", qPz, &b_qPz);
   fChain->SetBranchAddress("qMass", qMass, &b_qMass);
   fChain->SetBranchAddress("qXX", qXX, &b_qXX);
   fChain->SetBranchAddress("qYY", qYY, &b_qYY);
   fChain->SetBranchAddress("qZZ", qZZ, &b_qZZ);
   fChain->SetBranchAddress("qTT", qTT, &b_qTT);
   Notify();
}

Bool_t ptyv2::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void ptyv2::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t ptyv2::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef ptyv2_cxx
