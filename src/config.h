//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Thu Nov 21 12:12:27 2013 by ROOT version 5.34/05
// from TTree config/run configuration
// found on file: nt_sbs_E11.0.root
//////////////////////////////////////////////////////////

#ifndef config_h
#define config_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.

// Fixed size dimensions of array or collections stored in the TTree if any.

class config {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   Int_t           Run;
   Int_t           SkimLevel;
   Int_t           BookTrees;
   Double_t        Beam;
   Double_t        BeamTiltedAngle;
   Double_t        TargetM;
   Double_t        TargetL;
   Double_t        TargetAtomicNumber;
   Double_t        TargetNeutronNumber;
   Int_t           WhereToStopRecon;
   Int_t           SnakeModel;
   Double_t        BPMYRes;
   Double_t        BPMXRes;
   Int_t           UseSeptumPlusStdHRS;
   Int_t           UseOpticsDB;
   Double_t        TargetXOffset;
   Double_t        TargetYOffset;
   Double_t        TargetZOffset;
   Double_t        PivotXOffset;
   Double_t        PivotYOffset;
   Double_t        PivotZOffset;
   Int_t           SetupLHRS;
   Int_t           SetupRHRS;
   Double_t        LHRSMomentum;
   Double_t        RHRSMomentum;
   Double_t        LHRSAngle;
   Double_t        RHRSAngle;
   Double_t        Pivot2LHRSVBFace;
   Double_t        Pivot2RHRSVBFace;
   Int_t           UseHelmField;
   Double_t        HelmXOffset;
   Double_t        HelmYOffset;
   Double_t        HelmZOffset;
   Double_t        HelmRotAxis1;
   Double_t        HelmRotAxis2;
   Double_t        HelmRotAxis3;
   Double_t        HelmRotAngle1;
   Double_t        HelmRotAngle2;
   Double_t        HelmRotAngle3;
   Double_t        HelmCurrentRatio;
   Int_t           UseSeptumField;
   Double_t        SeptumXOffset;
   Double_t        SeptumYOffset;
   Double_t        SeptumZOffset;
   Double_t        SeptumRotAxis1;
   Double_t        SeptumRotAxis2;
   Double_t        SeptumRotAxis3;
   Double_t        SeptumRotAngle1;
   Double_t        SeptumRotAngle2;
   Double_t        SeptumRotAngle3;
   Double_t        SeptumCurrentRatioL;
   Double_t        SeptumCurrentRatioR;
   Int_t           SetupG2PTarget;
   Int_t           TargetType;
   Double_t        ThirdArmAngle;
   Double_t        SetupThirdArmVD;
   Double_t        ThirdArmRotZAngle;
   Double_t        Pivot2ThirdArmFace;
   Int_t           SetupChicane;
   Int_t           SetupChicaneVD;
   Double_t        FZB1TiltedAngle;
   Double_t        FZB1PosX;
   Double_t        FZB1PosY;
   Double_t        FZB1PosZ;
   Double_t        FZB1Bx;
   Double_t        FZB1By;
   Double_t        FZB1Bz;
   Double_t        FZB2TiltedAngle;
   Double_t        FZB2PosX;
   Double_t        FZB2PosY;
   Double_t        FZB2PosZ;
   Double_t        FZB2Bx;
   Double_t        FZB2By;
   Double_t        FZB2Bz;
   Int_t           SetupVD;
   Double_t        VDAngle;
   Double_t        VDTiltAngle;
   Double_t        Pivot2VDFace;
   Int_t           SetupBigBite;
   Double_t        BigBiteAngle;
   Double_t        BigBiteTiltAngle;
   Double_t        Pivot2BigBiteFace;
   Int_t           SetupHAND;
   Double_t        Pivot2HANDLeadWall;
   Int_t           SetupSuperBigBite;
   Double_t        SuperBigBiteAngle;
   Double_t        Pivot2SuperBigBiteFace;
   Int_t           UseSBSField;
   Double_t        SBSFieldXOffset;
   Double_t        SBSFieldYOffset;
   Double_t        SBSFieldZOffset;
   Double_t        SBSFieldRotAxis1;
   Double_t        SBSFieldRotAxis2;
   Double_t        SBSFieldRotAxis3;
   Double_t        SBSFieldRotAngle1;
   Double_t        SBSFieldRotAngle2;
   Double_t        SBSFieldRotAngle3;
   Double_t        SBSFieldCurrentRatio;
   Int_t           SetupHMS;
   Double_t        HMSAngle;
   Double_t        Pivot2HMSFace;
   Int_t           SetupRTPC;
   Double_t        RatioHe2DME;
   Int_t           SDNum;
   Char_t          SDNameForSDID0[100];
   Char_t          SDNameForSDID1[100];

   // List of branches
   TBranch        *b_Run;   //!
   TBranch        *b_SkimLevel;   //!
   TBranch        *b_BookTrees;   //!
   TBranch        *b_Beam;   //!
   TBranch        *b_BeamTiltedAngle;   //!
   TBranch        *b_TargetM;   //!
   TBranch        *b_TargetL;   //!
   TBranch        *b_TargetAtomicNumber;   //!
   TBranch        *b_TargetNeutronNumber;   //!
   TBranch        *b_WhereToStopRecon;   //!
   TBranch        *b_SnakeModel;   //!
   TBranch        *b_BPMYRes;   //!
   TBranch        *b_BPMXRes;   //!
   TBranch        *b_UseSeptumPlusStdHRS;   //!
   TBranch        *b_UseOpticsDB;   //!
   TBranch        *b_TargetXOffset;   //!
   TBranch        *b_TargetYOffset;   //!
   TBranch        *b_TargetZOffset;   //!
   TBranch        *b_PivotXOffset;   //!
   TBranch        *b_PivotYOffset;   //!
   TBranch        *b_PivotZOffset;   //!
   TBranch        *b_SetupLHRS;   //!
   TBranch        *b_SetupRHRS;   //!
   TBranch        *b_LHRSMomentum;   //!
   TBranch        *b_RHRSMomentum;   //!
   TBranch        *b_LHRSAngle;   //!
   TBranch        *b_RHRSAngle;   //!
   TBranch        *b_Pivot2LHRSVBFace;   //!
   TBranch        *b_Pivot2RHRSVBFace;   //!
   TBranch        *b_UseHelmField;   //!
   TBranch        *b_HelmXOffset;   //!
   TBranch        *b_HelmYOffset;   //!
   TBranch        *b_HelmZOffset;   //!
   TBranch        *b_HelmRotAxis1;   //!
   TBranch        *b_HelmRotAxis2;   //!
   TBranch        *b_HelmRotAxis3;   //!
   TBranch        *b_HelmRotAngle1;   //!
   TBranch        *b_HelmRotAngle2;   //!
   TBranch        *b_HelmRotAngle3;   //!
   TBranch        *b_HelmCurrentRatio;   //!
   TBranch        *b_UseSeptumField;   //!
   TBranch        *b_SeptumXOffset;   //!
   TBranch        *b_SeptumYOffset;   //!
   TBranch        *b_SeptumZOffset;   //!
   TBranch        *b_SeptumRotAxis1;   //!
   TBranch        *b_SeptumRotAxis2;   //!
   TBranch        *b_SeptumRotAxis3;   //!
   TBranch        *b_SeptumRotAngle1;   //!
   TBranch        *b_SeptumRotAngle2;   //!
   TBranch        *b_SeptumRotAngle3;   //!
   TBranch        *b_SeptumCurrentRatioL;   //!
   TBranch        *b_SeptumCurrentRatioR;   //!
   TBranch        *b_SetupG2PTarget;   //!
   TBranch        *b_TargetType;   //!
   TBranch        *b_ThirdArmAngle;   //!
   TBranch        *b_SetupThirdArmVD;   //!
   TBranch        *b_ThirdArmRotZAngle;   //!
   TBranch        *b_Pivot2ThirdArmFace;   //!
   TBranch        *b_SetupChicane;   //!
   TBranch        *b_SetupChicaneVD;   //!
   TBranch        *b_FZB1TiltedAngle;   //!
   TBranch        *b_FZB1PosX;   //!
   TBranch        *b_FZB1PosY;   //!
   TBranch        *b_FZB1PosZ;   //!
   TBranch        *b_FZB1Bx;   //!
   TBranch        *b_FZB1By;   //!
   TBranch        *b_FZB1Bz;   //!
   TBranch        *b_FZB2TiltedAngle;   //!
   TBranch        *b_FZB2PosX;   //!
   TBranch        *b_FZB2PosY;   //!
   TBranch        *b_FZB2PosZ;   //!
   TBranch        *b_FZB2Bx;   //!
   TBranch        *b_FZB2By;   //!
   TBranch        *b_FZB2Bz;   //!
   TBranch        *b_SetupVD;   //!
   TBranch        *b_VDAngle;   //!
   TBranch        *b_VDTiltAngle;   //!
   TBranch        *b_Pivot2VDFace;   //!
   TBranch        *b_SetupBigBite;   //!
   TBranch        *b_BigBiteAngle;   //!
   TBranch        *b_BigBiteTiltAngle;   //!
   TBranch        *b_Pivot2BigBiteFace;   //!
   TBranch        *b_SetupHAND;   //!
   TBranch        *b_Pivot2HANDLeadWall;   //!
   TBranch        *b_SetupSuperBigBite;   //!
   TBranch        *b_SuperBigBiteAngle;   //!
   TBranch        *b_Pivot2SuperBigBiteFace;   //!
   TBranch        *b_UseSBSField;   //!
   TBranch        *b_SBSFieldXOffset;   //!
   TBranch        *b_SBSFieldYOffset;   //!
   TBranch        *b_SBSFieldZOffset;   //!
   TBranch        *b_SBSFieldRotAxis1;   //!
   TBranch        *b_SBSFieldRotAxis2;   //!
   TBranch        *b_SBSFieldRotAxis3;   //!
   TBranch        *b_SBSFieldRotAngle1;   //!
   TBranch        *b_SBSFieldRotAngle2;   //!
   TBranch        *b_SBSFieldRotAngle3;   //!
   TBranch        *b_SBSFieldCurrentRatio;   //!
   TBranch        *b_SetupHMS;   //!
   TBranch        *b_HMSAngle;   //!
   TBranch        *b_Pivot2HMSFace;   //!
   TBranch        *b_SetupRTPC;   //!
   TBranch        *b_RatioHe2DME;   //!
   TBranch        *b_SDNum;   //!
   TBranch        *b_SDNameForSDID0;   //!
   TBranch        *b_SDNameForSDID1;   //!

   config(TTree *tree=0);
   virtual ~config();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef config_cxx
config::config(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      tree=(TTree*)gDirectory->Get("config");
   }
   Init(tree);
}

config::~config()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t config::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t config::LoadTree(Long64_t entry)
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

void config::Init(TTree *tree)
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

   fChain->SetBranchAddress("Run", &Run, &b_Run);
   fChain->SetBranchAddress("SkimLevel", &SkimLevel, &b_SkimLevel);
   fChain->SetBranchAddress("BookTrees", &BookTrees, &b_BookTrees);
   fChain->SetBranchAddress("Beam", &Beam, &b_Beam);
   fChain->SetBranchAddress("BeamTiltedAngle", &BeamTiltedAngle, &b_BeamTiltedAngle);
   fChain->SetBranchAddress("TargetM", &TargetM, &b_TargetM);
   fChain->SetBranchAddress("TargetL", &TargetL, &b_TargetL);
   fChain->SetBranchAddress("TargetAtomicNumber", &TargetAtomicNumber, &b_TargetAtomicNumber);
   fChain->SetBranchAddress("TargetNeutronNumber", &TargetNeutronNumber, &b_TargetNeutronNumber);
   fChain->SetBranchAddress("WhereToStopRecon", &WhereToStopRecon, &b_WhereToStopRecon);
   fChain->SetBranchAddress("SnakeModel", &SnakeModel, &b_SnakeModel);
   fChain->SetBranchAddress("BPMYRes", &BPMYRes, &b_BPMYRes);
   fChain->SetBranchAddress("BPMXRes", &BPMXRes, &b_BPMXRes);
   fChain->SetBranchAddress("UseSeptumPlusStdHRS", &UseSeptumPlusStdHRS, &b_UseSeptumPlusStdHRS);
   fChain->SetBranchAddress("UseOpticsDB", &UseOpticsDB, &b_UseOpticsDB);
   fChain->SetBranchAddress("TargetXOffset", &TargetXOffset, &b_TargetXOffset);
   fChain->SetBranchAddress("TargetYOffset", &TargetYOffset, &b_TargetYOffset);
   fChain->SetBranchAddress("TargetZOffset", &TargetZOffset, &b_TargetZOffset);
   fChain->SetBranchAddress("PivotXOffset", &PivotXOffset, &b_PivotXOffset);
   fChain->SetBranchAddress("PivotYOffset", &PivotYOffset, &b_PivotYOffset);
   fChain->SetBranchAddress("PivotZOffset", &PivotZOffset, &b_PivotZOffset);
   fChain->SetBranchAddress("SetupLHRS", &SetupLHRS, &b_SetupLHRS);
   fChain->SetBranchAddress("SetupRHRS", &SetupRHRS, &b_SetupRHRS);
   fChain->SetBranchAddress("LHRSMomentum", &LHRSMomentum, &b_LHRSMomentum);
   fChain->SetBranchAddress("RHRSMomentum", &RHRSMomentum, &b_RHRSMomentum);
   fChain->SetBranchAddress("LHRSAngle", &LHRSAngle, &b_LHRSAngle);
   fChain->SetBranchAddress("RHRSAngle", &RHRSAngle, &b_RHRSAngle);
   fChain->SetBranchAddress("Pivot2LHRSVBFace", &Pivot2LHRSVBFace, &b_Pivot2LHRSVBFace);
   fChain->SetBranchAddress("Pivot2RHRSVBFace", &Pivot2RHRSVBFace, &b_Pivot2RHRSVBFace);
   fChain->SetBranchAddress("UseHelmField", &UseHelmField, &b_UseHelmField);
   fChain->SetBranchAddress("HelmXOffset", &HelmXOffset, &b_HelmXOffset);
   fChain->SetBranchAddress("HelmYOffset", &HelmYOffset, &b_HelmYOffset);
   fChain->SetBranchAddress("HelmZOffset", &HelmZOffset, &b_HelmZOffset);
   fChain->SetBranchAddress("HelmRotAxis1", &HelmRotAxis1, &b_HelmRotAxis1);
   fChain->SetBranchAddress("HelmRotAxis2", &HelmRotAxis2, &b_HelmRotAxis2);
   fChain->SetBranchAddress("HelmRotAxis3", &HelmRotAxis3, &b_HelmRotAxis3);
   fChain->SetBranchAddress("HelmRotAngle1", &HelmRotAngle1, &b_HelmRotAngle1);
   fChain->SetBranchAddress("HelmRotAngle2", &HelmRotAngle2, &b_HelmRotAngle2);
   fChain->SetBranchAddress("HelmRotAngle3", &HelmRotAngle3, &b_HelmRotAngle3);
   fChain->SetBranchAddress("HelmCurrentRatio", &HelmCurrentRatio, &b_HelmCurrentRatio);
   fChain->SetBranchAddress("UseSeptumField", &UseSeptumField, &b_UseSeptumField);
   fChain->SetBranchAddress("SeptumXOffset", &SeptumXOffset, &b_SeptumXOffset);
   fChain->SetBranchAddress("SeptumYOffset", &SeptumYOffset, &b_SeptumYOffset);
   fChain->SetBranchAddress("SeptumZOffset", &SeptumZOffset, &b_SeptumZOffset);
   fChain->SetBranchAddress("SeptumRotAxis1", &SeptumRotAxis1, &b_SeptumRotAxis1);
   fChain->SetBranchAddress("SeptumRotAxis2", &SeptumRotAxis2, &b_SeptumRotAxis2);
   fChain->SetBranchAddress("SeptumRotAxis3", &SeptumRotAxis3, &b_SeptumRotAxis3);
   fChain->SetBranchAddress("SeptumRotAngle1", &SeptumRotAngle1, &b_SeptumRotAngle1);
   fChain->SetBranchAddress("SeptumRotAngle2", &SeptumRotAngle2, &b_SeptumRotAngle2);
   fChain->SetBranchAddress("SeptumRotAngle3", &SeptumRotAngle3, &b_SeptumRotAngle3);
   fChain->SetBranchAddress("SeptumCurrentRatioL", &SeptumCurrentRatioL, &b_SeptumCurrentRatioL);
   fChain->SetBranchAddress("SeptumCurrentRatioR", &SeptumCurrentRatioR, &b_SeptumCurrentRatioR);
   fChain->SetBranchAddress("SetupG2PTarget", &SetupG2PTarget, &b_SetupG2PTarget);
   fChain->SetBranchAddress("TargetType", &TargetType, &b_TargetType);
   fChain->SetBranchAddress("ThirdArmAngle", &ThirdArmAngle, &b_ThirdArmAngle);
   fChain->SetBranchAddress("SetupThirdArmVD", &SetupThirdArmVD, &b_SetupThirdArmVD);
   fChain->SetBranchAddress("ThirdArmRotZAngle", &ThirdArmRotZAngle, &b_ThirdArmRotZAngle);
   fChain->SetBranchAddress("Pivot2ThirdArmFace", &Pivot2ThirdArmFace, &b_Pivot2ThirdArmFace);
   fChain->SetBranchAddress("SetupChicane", &SetupChicane, &b_SetupChicane);
   fChain->SetBranchAddress("SetupChicaneVD", &SetupChicaneVD, &b_SetupChicaneVD);
   fChain->SetBranchAddress("FZB1TiltedAngle", &FZB1TiltedAngle, &b_FZB1TiltedAngle);
   fChain->SetBranchAddress("FZB1PosX", &FZB1PosX, &b_FZB1PosX);
   fChain->SetBranchAddress("FZB1PosY", &FZB1PosY, &b_FZB1PosY);
   fChain->SetBranchAddress("FZB1PosZ", &FZB1PosZ, &b_FZB1PosZ);
   fChain->SetBranchAddress("FZB1Bx", &FZB1Bx, &b_FZB1Bx);
   fChain->SetBranchAddress("FZB1By", &FZB1By, &b_FZB1By);
   fChain->SetBranchAddress("FZB1Bz", &FZB1Bz, &b_FZB1Bz);
   fChain->SetBranchAddress("FZB2TiltedAngle", &FZB2TiltedAngle, &b_FZB2TiltedAngle);
   fChain->SetBranchAddress("FZB2PosX", &FZB2PosX, &b_FZB2PosX);
   fChain->SetBranchAddress("FZB2PosY", &FZB2PosY, &b_FZB2PosY);
   fChain->SetBranchAddress("FZB2PosZ", &FZB2PosZ, &b_FZB2PosZ);
   fChain->SetBranchAddress("FZB2Bx", &FZB2Bx, &b_FZB2Bx);
   fChain->SetBranchAddress("FZB2By", &FZB2By, &b_FZB2By);
   fChain->SetBranchAddress("FZB2Bz", &FZB2Bz, &b_FZB2Bz);
   fChain->SetBranchAddress("SetupVD", &SetupVD, &b_SetupVD);
   fChain->SetBranchAddress("VDAngle", &VDAngle, &b_VDAngle);
   fChain->SetBranchAddress("VDTiltAngle", &VDTiltAngle, &b_VDTiltAngle);
   fChain->SetBranchAddress("Pivot2VDFace", &Pivot2VDFace, &b_Pivot2VDFace);
   fChain->SetBranchAddress("SetupBigBite", &SetupBigBite, &b_SetupBigBite);
   fChain->SetBranchAddress("BigBiteAngle", &BigBiteAngle, &b_BigBiteAngle);
   fChain->SetBranchAddress("BigBiteTiltAngle", &BigBiteTiltAngle, &b_BigBiteTiltAngle);
   fChain->SetBranchAddress("Pivot2BigBiteFace", &Pivot2BigBiteFace, &b_Pivot2BigBiteFace);
   fChain->SetBranchAddress("SetupHAND", &SetupHAND, &b_SetupHAND);
   fChain->SetBranchAddress("Pivot2HANDLeadWall", &Pivot2HANDLeadWall, &b_Pivot2HANDLeadWall);
   fChain->SetBranchAddress("SetupSuperBigBite", &SetupSuperBigBite, &b_SetupSuperBigBite);
   fChain->SetBranchAddress("SuperBigBiteAngle", &SuperBigBiteAngle, &b_SuperBigBiteAngle);
   fChain->SetBranchAddress("Pivot2SuperBigBiteFace", &Pivot2SuperBigBiteFace, &b_Pivot2SuperBigBiteFace);
   fChain->SetBranchAddress("UseSBSField", &UseSBSField, &b_UseSBSField);
   fChain->SetBranchAddress("SBSFieldXOffset", &SBSFieldXOffset, &b_SBSFieldXOffset);
   fChain->SetBranchAddress("SBSFieldYOffset", &SBSFieldYOffset, &b_SBSFieldYOffset);
   fChain->SetBranchAddress("SBSFieldZOffset", &SBSFieldZOffset, &b_SBSFieldZOffset);
   fChain->SetBranchAddress("SBSFieldRotAxis1", &SBSFieldRotAxis1, &b_SBSFieldRotAxis1);
   fChain->SetBranchAddress("SBSFieldRotAxis2", &SBSFieldRotAxis2, &b_SBSFieldRotAxis2);
   fChain->SetBranchAddress("SBSFieldRotAxis3", &SBSFieldRotAxis3, &b_SBSFieldRotAxis3);
   fChain->SetBranchAddress("SBSFieldRotAngle1", &SBSFieldRotAngle1, &b_SBSFieldRotAngle1);
   fChain->SetBranchAddress("SBSFieldRotAngle2", &SBSFieldRotAngle2, &b_SBSFieldRotAngle2);
   fChain->SetBranchAddress("SBSFieldRotAngle3", &SBSFieldRotAngle3, &b_SBSFieldRotAngle3);
   fChain->SetBranchAddress("SBSFieldCurrentRatio", &SBSFieldCurrentRatio, &b_SBSFieldCurrentRatio);
   fChain->SetBranchAddress("SetupHMS", &SetupHMS, &b_SetupHMS);
   fChain->SetBranchAddress("HMSAngle", &HMSAngle, &b_HMSAngle);
   fChain->SetBranchAddress("Pivot2HMSFace", &Pivot2HMSFace, &b_Pivot2HMSFace);
   fChain->SetBranchAddress("SetupRTPC", &SetupRTPC, &b_SetupRTPC);
   fChain->SetBranchAddress("RatioHe2DME", &RatioHe2DME, &b_RatioHe2DME);
   fChain->SetBranchAddress("SDNum", &SDNum, &b_SDNum);
   fChain->SetBranchAddress("SDNameForSDID0", SDNameForSDID0, &b_SDNameForSDID0);
   fChain->SetBranchAddress("SDNameForSDID1", SDNameForSDID1, &b_SDNameForSDID1);
   Notify();
}

Bool_t config::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void config::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t config::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef config_cxx
