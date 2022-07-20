
//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Wed Jul  6 19:52:18 2022 by ROOT version 5.34/30
// from TTree tr/Reconst ntuple
// found on file: AuAu200GeV_SM3mb_Alpha0p47_NT150_DecayOff_0_20fm_file10.root
//////////////////////////////////////////////////////////

#ifndef check_h
#define check_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.

// Fixed size dimensions of array or collections stored in the TTree if any.

class check {                                                                   //*Check Class
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
      Int_t           Stat[394];   //[Nab]
      Int_t           PID[8982];   //[Mult]
      Float_t         Px[8982];   //[Mult]
      Float_t         Py[8982];   //[Mult]
      Float_t         Pz[8982];   //[Mult]
      Float_t         Mass[8982];   //[Mult]
      Float_t         XX[8982];   //[Mult]
      Float_t         YY[8982];   //[Mult]
      Float_t         ZZ[8982];   //[Mult]
      Float_t         TT[8982];   //[Mult]
      Int_t           qMult;
      Int_t           qPID[12916];   //[qMult]
      Float_t         qPx[12916];   //[qMult]
      Float_t         qPy[12916];   //[qMult]
      Float_t         qPz[12916];   //[qMult]
      Float_t         qMass[12916];   //[qMult]
      Float_t         qXX[12916];   //[qMult]
      Float_t         qYY[12916];   //[qMult]
      Float_t         qZZ[12916];   //[qMult]
      Float_t         qTT[12916];   //[qMult]

      // List of branches
      TBranch        *b_Event;    //!Event
      TBranch        *b_Mult;     //! multiplicity = tracks
      TBranch        *b_Npartp;   //!No of Participant on projectile
      TBranch        *b_Npartt;   //!No of Participant on Target
      TBranch        *b_Nesp;     //!Elastic
      TBranch        *b_Ninesp;   //!Inelastic
      TBranch        *b_Nest;     //!Elastic
      TBranch        *b_Ninest;   //!Inelastic
      TBranch        *b_Imp;   //?Mass/impact parameter
      TBranch        *b_Na;   //!Projectile mass No.
      TBranch        *b_Nb;   //!Target mass No.
      TBranch        *b_Nab;   //!Nab = Na+Nb;
      TBranch        *b_Psi;   //!Event Plane Angle
      TBranch        *b_Theta_p;   //!projetileTheta
      TBranch        *b_Phi_p;   //!
      TBranch        *b_Theta_t;   //!Target Theta
      TBranch        *b_Phi_t;   //!
      TBranch        *b_Nx;   //!
      TBranch        *b_Ny;   //!
      TBranch        *b_Nz;   //!
      TBranch        *b_Stat;   //!
      TBranch        *b_PID;   //!Particle ID
      TBranch        *b_Px;   //!Momentum in x direction for particle i
      TBranch        *b_Py;   //!
      TBranch        *b_Pz;   //!
      TBranch        *b_Mass;   //!Mass of ith particle
      TBranch        *b_XX;   //!Position in x direction for particle i
      TBranch        *b_YY;   //!
      TBranch        *b_ZZ;   //!
      TBranch        *b_TT;   //!Time
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

      check(TTree *tree=0);                   //?IG this works like __init__ in python   ...for definign the drguments
      virtual ~check();
      virtual Int_t    Cut(Long64_t entry);
      virtual Int_t    GetEntry(Long64_t entry);
      virtual Long64_t LoadTree(Long64_t entry);
      virtual void     Init(TTree *tree);
      virtual void     Loop(double lowerbc, double higherbc);
      virtual void     Multiplicity();
      virtual void     Printer();
      virtual Bool_t   Notify();
      virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef check_cxx
check::check(TTree *tree) : fChain(0) {
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

check::~check(){                 //~ implies distructor
   if (!fChain) return;                 //fchain is tree/chain
   delete fChain->GetCurrentFile();
}

Int_t check::GetEntry(Long64_t entry){              //*Long 64 for huge int values 
                                       // Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t check::LoadTree(Long64_t entry){
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

void check::Init(TTree *tree)
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

Bool_t check::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void check::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t check::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef check_cxx