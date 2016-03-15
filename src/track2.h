//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Thu Nov 21 12:08:55 2013 by ROOT version 5.34/05
// from TTree track2/track 3
// found on file: nt_sbs_E11.0.root
//////////////////////////////////////////////////////////

#ifndef track2_h
#define track2_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.

// Fixed size dimensions of array or collections stored in the TTree if any.

class track2 {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   Int_t           Index;
   Int_t           PdgId;
   Int_t           TrackId;
   Int_t           TrackClass;
   Double_t        X0;
   Double_t        Y0;
   Double_t        Z0;
   Double_t        P0;
   Double_t        Theta0;
   Double_t        Phi0;
   Double_t        X0_tr;
   Double_t        Y0_tr;
   Double_t        Z0_tr;
   Double_t        Theta0_tr;
   Double_t        Phi0_tr;
   Double_t        Xtg_tr;
   Double_t        Ytg_tr;
   Double_t        Thetatg_tr;
   Double_t        Phitg_tr;
   Double_t        Xvb;
   Double_t        Yvb;
   Double_t        Zvb;
   Double_t        Pvb;
   Double_t        Thetavb;
   Double_t        Phivb;
   Double_t        Xvb_tr;
   Double_t        Yvb_tr;
   Double_t        Zvb_tr;
   Double_t        Thetavb_tr;
   Double_t        Phivb_tr;
   Double_t        Xfp_tr;
   Double_t        Yfp_tr;
   Double_t        Thetafp_tr;
   Double_t        Phifp_tr;
   Double_t        Xtg_rec_tr;
   Double_t        Ytg_rec_tr;
   Double_t        Thetatg_rec_tr;
   Double_t        Phitg_rec_tr;
   Double_t        X_rec_tr;
   Double_t        Y_rec_tr;
   Double_t        Z_rec_tr;
   Double_t        Theta_rec_tr;
   Double_t        Phi_rec_tr;
   Double_t        X_rec;
   Double_t        Y_rec;
   Double_t        Z_rec;
   Double_t        P_rec;
   Double_t        Theta_rec;
   Double_t        Phi_rec;
   Double_t        Delta;
   Double_t        Delta_rec;
   Double_t        Xtg_rec_db_tr;
   Double_t        Ytg_rec_db_tr;
   Double_t        Thetatg_rec_db_tr;
   Double_t        Phitg_rec_db_tr;
   Double_t        Delta_rec_db;
   Double_t        X_proj2tg_tr;
   Double_t        Y_proj2tg_tr;
   Double_t        X_rec2tg_tr;
   Double_t        Y_rec2tg_tr;
   Double_t        Theta_rec2tg_tr;
   Double_t        Phi_rec2tg_tr;
   Double_t        P_rec2tg;
   Double_t        X_proj2sl_tr;
   Double_t        Y_proj2sl_tr;
   Double_t        TrackRadlen;
   Double_t        Theta0Eff;
   Double_t        ElasXS;
   Double_t        Ei;
   Double_t        XS;
   Int_t           StepNum;
   Double_t        StepX[1024];   //[StepNum]
   Double_t        StepY[1024];   //[StepNum]
   Double_t        StepZ[1024];   //[StepNum]
   Double_t        StepdE[1024];   //[StepNum]
   Double_t        StepL[1024];   //[StepNum]
   Double_t        StepEkin[1024];   //[StepNum]
   Double_t        StepTL[1024];   //[StepNum]
   Double_t        StepRadlen[1024];   //[StepNum]
   Double_t        StepDsty[1024];   //[StepNum]
   Double_t        StepBx[1024];   //[StepNum]
   Double_t        StepBy[1024];   //[StepNum]
   Double_t        StepBz[1024];   //[StepNum]
   Double_t        TrackBdLx;
   Double_t        TrackBdLy;
   Double_t        TrackBdLz;
   Double_t        R0;
   Double_t        A0;
   Double_t        B0;

   // List of branches
   TBranch        *b_Index;   //!
   TBranch        *b_PdgId;   //!
   TBranch        *b_TrackId;   //!
   TBranch        *b_TrackClass;   //!
   TBranch        *b_X0;   //!
   TBranch        *b_Y0;   //!
   TBranch        *b_Z0;   //!
   TBranch        *b_P0;   //!
   TBranch        *b_Theta0;   //!
   TBranch        *b_Phi0;   //!
   TBranch        *b_X0_tr;   //!
   TBranch        *b_Y0_tr;   //!
   TBranch        *b_Z0_tr;   //!
   TBranch        *b_Theta0_tr;   //!
   TBranch        *b_Phi0_tr;   //!
   TBranch        *b_Xtg_tr;   //!
   TBranch        *b_Ytg_tr;   //!
   TBranch        *b_Thetatg_tr;   //!
   TBranch        *b_Phitg_tr;   //!
   TBranch        *b_Xvb;   //!
   TBranch        *b_Yvb;   //!
   TBranch        *b_Zvb;   //!
   TBranch        *b_Pvb;   //!
   TBranch        *b_Thetavb;   //!
   TBranch        *b_Phivb;   //!
   TBranch        *b_Xvb_tr;   //!
   TBranch        *b_Yvb_tr;   //!
   TBranch        *b_Zvb_tr;   //!
   TBranch        *b_Thetavb_tr;   //!
   TBranch        *b_Phivb_tr;   //!
   TBranch        *b_Xfp_tr;   //!
   TBranch        *b_Yfp_tr;   //!
   TBranch        *b_Thetafp_tr;   //!
   TBranch        *b_Phifp_tr;   //!
   TBranch        *b_Xtg_rec_tr;   //!
   TBranch        *b_Ytg_rec_tr;   //!
   TBranch        *b_Thetatg_rec_tr;   //!
   TBranch        *b_Phitg_rec_tr;   //!
   TBranch        *b_X_rec_tr;   //!
   TBranch        *b_Y_rec_tr;   //!
   TBranch        *b_Z_rec_tr;   //!
   TBranch        *b_Theta_rec_tr;   //!
   TBranch        *b_Phi_rec_tr;   //!
   TBranch        *b_X_rec;   //!
   TBranch        *b_Y_rec;   //!
   TBranch        *b_Z_rec;   //!
   TBranch        *b_P_rec;   //!
   TBranch        *b_Theta_rec;   //!
   TBranch        *b_Phi_rec;   //!
   TBranch        *b_Delta;   //!
   TBranch        *b_Delta_rec;   //!
   TBranch        *b_Xtg_rec_db_tr;   //!
   TBranch        *b_Ytg_rec_db_tr;   //!
   TBranch        *b_Thetatg_rec_db_tr;   //!
   TBranch        *b_Phitg_rec_db_tr;   //!
   TBranch        *b_Delta_rec_db;   //!
   TBranch        *b_X_proj2tg_tr;   //!
   TBranch        *b_Y_proj2tg_tr;   //!
   TBranch        *b_X_rec2tg_tr;   //!
   TBranch        *b_Y_rec2tg_tr;   //!
   TBranch        *b_Theta_rec2tg_tr;   //!
   TBranch        *b_Phi_rec2tg_tr;   //!
   TBranch        *b_P_rec2tg;   //!
   TBranch        *b_X_proj2sl_tr;   //!
   TBranch        *b_Y_proj2sl_tr;   //!
   TBranch        *b_TrackRadlen;   //!
   TBranch        *b_Theta0Eff;   //!
   TBranch        *b_ElasXS;   //!
   TBranch        *b_Ei;   //!
   TBranch        *b_XS;   //!
   TBranch        *b_StepNum;   //!
   TBranch        *b_StepX;   //!
   TBranch        *b_StepY;   //!
   TBranch        *b_StepZ;   //!
   TBranch        *b_StepdE;   //!
   TBranch        *b_StepL;   //!
   TBranch        *b_StepEkin;   //!
   TBranch        *b_StepTL;   //!
   TBranch        *b_StepRadlen;   //!
   TBranch        *b_StepDsty;   //!
   TBranch        *b_StepBx;   //!
   TBranch        *b_StepBy;   //!
   TBranch        *b_StepBz;   //!
   TBranch        *b_TrackBdLx;   //!
   TBranch        *b_TrackBdLy;   //!
   TBranch        *b_TrackBdLz;   //!
   TBranch        *b_R0;   //!
   TBranch        *b_A0;   //!
   TBranch        *b_B0;   //!

   track2(TTree *tree=0);
   virtual ~track2();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef track2_cxx
track2::track2(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
     tree=(TTree*)gDirectory->Get("track2");
   }
   Init(tree);
}

track2::~track2()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t track2::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t track2::LoadTree(Long64_t entry)
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

void track2::Init(TTree *tree)
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

   fChain->SetBranchAddress("Index", &Index, &b_Index);
   fChain->SetBranchAddress("PdgId", &PdgId, &b_PdgId);
   fChain->SetBranchAddress("TrackId", &TrackId, &b_TrackId);
   fChain->SetBranchAddress("TrackClass", &TrackClass, &b_TrackClass);
   fChain->SetBranchAddress("X0", &X0, &b_X0);
   fChain->SetBranchAddress("Y0", &Y0, &b_Y0);
   fChain->SetBranchAddress("Z0", &Z0, &b_Z0);
   fChain->SetBranchAddress("P0", &P0, &b_P0);
   fChain->SetBranchAddress("Theta0", &Theta0, &b_Theta0);
   fChain->SetBranchAddress("Phi0", &Phi0, &b_Phi0);
   fChain->SetBranchAddress("X0_tr", &X0_tr, &b_X0_tr);
   fChain->SetBranchAddress("Y0_tr", &Y0_tr, &b_Y0_tr);
   fChain->SetBranchAddress("Z0_tr", &Z0_tr, &b_Z0_tr);
   fChain->SetBranchAddress("Theta0_tr", &Theta0_tr, &b_Theta0_tr);
   fChain->SetBranchAddress("Phi0_tr", &Phi0_tr, &b_Phi0_tr);
   fChain->SetBranchAddress("Xtg_tr", &Xtg_tr, &b_Xtg_tr);
   fChain->SetBranchAddress("Ytg_tr", &Ytg_tr, &b_Ytg_tr);
   fChain->SetBranchAddress("Thetatg_tr", &Thetatg_tr, &b_Thetatg_tr);
   fChain->SetBranchAddress("Phitg_tr", &Phitg_tr, &b_Phitg_tr);
   fChain->SetBranchAddress("Xvb", &Xvb, &b_Xvb);
   fChain->SetBranchAddress("Yvb", &Yvb, &b_Yvb);
   fChain->SetBranchAddress("Zvb", &Zvb, &b_Zvb);
   fChain->SetBranchAddress("Pvb", &Pvb, &b_Pvb);
   fChain->SetBranchAddress("Thetavb", &Thetavb, &b_Thetavb);
   fChain->SetBranchAddress("Phivb", &Phivb, &b_Phivb);
   fChain->SetBranchAddress("Xvb_tr", &Xvb_tr, &b_Xvb_tr);
   fChain->SetBranchAddress("Yvb_tr", &Yvb_tr, &b_Yvb_tr);
   fChain->SetBranchAddress("Zvb_tr", &Zvb_tr, &b_Zvb_tr);
   fChain->SetBranchAddress("Thetavb_tr", &Thetavb_tr, &b_Thetavb_tr);
   fChain->SetBranchAddress("Phivb_tr", &Phivb_tr, &b_Phivb_tr);
   fChain->SetBranchAddress("Xfp_tr", &Xfp_tr, &b_Xfp_tr);
   fChain->SetBranchAddress("Yfp_tr", &Yfp_tr, &b_Yfp_tr);
   fChain->SetBranchAddress("Thetafp_tr", &Thetafp_tr, &b_Thetafp_tr);
   fChain->SetBranchAddress("Phifp_tr", &Phifp_tr, &b_Phifp_tr);
   fChain->SetBranchAddress("Xtg_rec_tr", &Xtg_rec_tr, &b_Xtg_rec_tr);
   fChain->SetBranchAddress("Ytg_rec_tr", &Ytg_rec_tr, &b_Ytg_rec_tr);
   fChain->SetBranchAddress("Thetatg_rec_tr", &Thetatg_rec_tr, &b_Thetatg_rec_tr);
   fChain->SetBranchAddress("Phitg_rec_tr", &Phitg_rec_tr, &b_Phitg_rec_tr);
   fChain->SetBranchAddress("X_rec_tr", &X_rec_tr, &b_X_rec_tr);
   fChain->SetBranchAddress("Y_rec_tr", &Y_rec_tr, &b_Y_rec_tr);
   fChain->SetBranchAddress("Z_rec_tr", &Z_rec_tr, &b_Z_rec_tr);
   fChain->SetBranchAddress("Theta_rec_tr", &Theta_rec_tr, &b_Theta_rec_tr);
   fChain->SetBranchAddress("Phi_rec_tr", &Phi_rec_tr, &b_Phi_rec_tr);
   fChain->SetBranchAddress("X_rec", &X_rec, &b_X_rec);
   fChain->SetBranchAddress("Y_rec", &Y_rec, &b_Y_rec);
   fChain->SetBranchAddress("Z_rec", &Z_rec, &b_Z_rec);
   fChain->SetBranchAddress("P_rec", &P_rec, &b_P_rec);
   fChain->SetBranchAddress("Theta_rec", &Theta_rec, &b_Theta_rec);
   fChain->SetBranchAddress("Phi_rec", &Phi_rec, &b_Phi_rec);
   fChain->SetBranchAddress("Delta", &Delta, &b_Delta);
   fChain->SetBranchAddress("Delta_rec", &Delta_rec, &b_Delta_rec);
   fChain->SetBranchAddress("Xtg_rec_db_tr", &Xtg_rec_db_tr, &b_Xtg_rec_db_tr);
   fChain->SetBranchAddress("Ytg_rec_db_tr", &Ytg_rec_db_tr, &b_Ytg_rec_db_tr);
   fChain->SetBranchAddress("Thetatg_rec_db_tr", &Thetatg_rec_db_tr, &b_Thetatg_rec_db_tr);
   fChain->SetBranchAddress("Phitg_rec_db_tr", &Phitg_rec_db_tr, &b_Phitg_rec_db_tr);
   fChain->SetBranchAddress("Delta_rec_db", &Delta_rec_db, &b_Delta_rec_db);
   fChain->SetBranchAddress("X_proj2tg_tr", &X_proj2tg_tr, &b_X_proj2tg_tr);
   fChain->SetBranchAddress("Y_proj2tg_tr", &Y_proj2tg_tr, &b_Y_proj2tg_tr);
   fChain->SetBranchAddress("X_rec2tg_tr", &X_rec2tg_tr, &b_X_rec2tg_tr);
   fChain->SetBranchAddress("Y_rec2tg_tr", &Y_rec2tg_tr, &b_Y_rec2tg_tr);
   fChain->SetBranchAddress("Theta_rec2tg_tr", &Theta_rec2tg_tr, &b_Theta_rec2tg_tr);
   fChain->SetBranchAddress("Phi_rec2tg_tr", &Phi_rec2tg_tr, &b_Phi_rec2tg_tr);
   fChain->SetBranchAddress("P_rec2tg", &P_rec2tg, &b_P_rec2tg);
   fChain->SetBranchAddress("X_proj2sl_tr", &X_proj2sl_tr, &b_X_proj2sl_tr);
   fChain->SetBranchAddress("Y_proj2sl_tr", &Y_proj2sl_tr, &b_Y_proj2sl_tr);
   fChain->SetBranchAddress("TrackRadlen", &TrackRadlen, &b_TrackRadlen);
   fChain->SetBranchAddress("Theta0Eff", &Theta0Eff, &b_Theta0Eff);
   fChain->SetBranchAddress("ElasXS", &ElasXS, &b_ElasXS);
   fChain->SetBranchAddress("Ei", &Ei, &b_Ei);
   fChain->SetBranchAddress("XS", &XS, &b_XS);
   fChain->SetBranchAddress("StepNum", &StepNum, &b_StepNum);
   fChain->SetBranchAddress("StepX", StepX, &b_StepX);
   fChain->SetBranchAddress("StepY", StepY, &b_StepY);
   fChain->SetBranchAddress("StepZ", StepZ, &b_StepZ);
   fChain->SetBranchAddress("StepdE", StepdE, &b_StepdE);
   fChain->SetBranchAddress("StepL", StepL, &b_StepL);
   fChain->SetBranchAddress("StepEkin", StepEkin, &b_StepEkin);
   fChain->SetBranchAddress("StepTL", StepTL, &b_StepTL);
   fChain->SetBranchAddress("StepRadlen", StepRadlen, &b_StepRadlen);
   fChain->SetBranchAddress("StepDsty", StepDsty, &b_StepDsty);
   fChain->SetBranchAddress("StepBx", StepBx, &b_StepBx);
   fChain->SetBranchAddress("StepBy", StepBy, &b_StepBy);
   fChain->SetBranchAddress("StepBz", StepBz, &b_StepBz);
   fChain->SetBranchAddress("TrackBdLx", &TrackBdLx, &b_TrackBdLx);
   fChain->SetBranchAddress("TrackBdLy", &TrackBdLy, &b_TrackBdLy);
   fChain->SetBranchAddress("TrackBdLz", &TrackBdLz, &b_TrackBdLz);
   fChain->SetBranchAddress("R0", &R0, &b_R0);
   fChain->SetBranchAddress("A0", &A0, &b_A0);
   fChain->SetBranchAddress("B0", &B0, &b_B0);
   Notify();
}

Bool_t track2::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void track2::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t track2::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef track2_cxx
