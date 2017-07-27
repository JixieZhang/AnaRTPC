//this script is used to reconstruct RTPC e-p events
//One can config BoNuS6 or BoNuS12 with either large pad or small pad
//The acceptance of RTPC and CLAS12 will be plot at the same time

#include "stdlib.h"
#include <iostream>
#include <fstream>
#include <iomanip>
#include "math.h"
#include "DriftEMagboltz.hh"
using namespace std;

#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TCanvas.h>
#include <TQObject.h>
#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TVector3.h"
#include "TString.h"
#include "TCut.h"
#include "TCutG.h"
#include "TPaveText.h"
#include "TText.h"
#include "TPad.h"
#include "TF1.h"
#include "TStyle.h"
#include "TGraph.h"
#include "TMultiGraph.h"
#include "TLorentzVector.h"
#include "TGraphErrors.h"

#include "track0.h"
#include "track1.h"
#include "track2.h"
#include "config.h"

#define MaxHit 200
#define MaxHit_rec 600

#define Store_Rec_Leaves 1

//#define RTPC_BoNuS6 1
//the following have been defined in AnaRTPC
extern double RTPC_Cathode_R; // = 30;    //mm
extern double RTPC_Anode_R; //   = 80;    //mm
extern double RTPC_ReadOut_R; // = 90;    //mm
extern double RTPC_TDC_Window; //= 25;    //ns
extern double RTPC_Pad_Z; // = 5.0;   //mm
extern double RTPC_Pad_W; // = 4.5;   //mm
extern double RTPC_Length; // = 400;            //mm

extern DriftEMagboltz *gEsim;

// estimate number of channels: assuming 3.0mm GEM1 glue area width
// Row# = int(2*PI*Readout_R/RTPC_Pad_W)
// Col# = int(RTPC_Length/RTPC_Pad_Z)
extern int GetNofCh(int padtype);

//config pad size:
//padtype:
//1) 4.5x5 mm, 5cm DriftRegion
//2) TDIS: 2x2 mm, 5<DriftRegion<15
//3) compass 2-D readout, 0.4x20 stripe, 0.4x0.4 mm equivalent, 5cm DriftRegion
//4) compass 2-D readout, 1.0x20 stripe, 1x1 mm equivalent, 5cm DriftRegion
//5) compass 2-D readout, 1.0x20 stripe, 1x1 mm equivalent, 4cm DriftRegion
//6) RTPC12 config: 2.8016x4.1016 mm, add 4mil gap between pads, R_in_readout=80.42
extern void ConfigPadSize(int padtype);

static track0 *Proton=0;
static track1 *Electron=0;
static config *Config=0;

//boxmuller gauss number generator
//input: mean m, standard deviation s
extern double fGaus(double m, double s);
extern void PrintAcc(TH3F *h3,TH3F *h3N, int Nmin=0, int TCS=0);

extern int FitATrack(int StepNum, double StepX[],double StepY[],double StepZ[],
                 double &R, double &A, double &B,
                 double &Phi_deg, double &Theta_deg, double &Z,
                 double RLimit_l=30.1, double RLimit_h=79.9);

extern void RTPC_Recon(double MeanBz, double R_rec, double Theta_rec, double Phi_rec,
                   double &P_corr, double &Theta_corr, double &Phi_corr);

extern const char *GetFileName(const char* str);

extern const char *GetBaseName(const char* str);


void CreateRTPCSuperEvent(const char *infile="nt.root", const char *outfile="nt_ep.root",
                    int padtype=2, int ntrack_per_event=25)
{
  ////////////////////////////////////////////////////////
  TFile *InFile = TFile::Open(infile);
  if(!InFile)  {
    cout<<"\nCan not open input file \""<<infile<<"\", I quit ...\n\n";
    return;
  }

  cout<<"infile="<<infile<<endl;
  ////////////////////////////////////////////////////////

  gStyle->SetOptStat(0);
  //set the draw-option "text" format
  gStyle->SetPaintTextFormat(".0f");

  ConfigPadSize(padtype);

  gEsim = new DriftEMagboltz(0.9, RTPC_Pad_W, RTPC_Pad_Z, RTPC_ReadOut_R, RTPC_Length);

  char pKey[100];
#ifdef RTPC_BoNuS6
  sprintf(pKey, "%s", "BoNuS6");
#else
  sprintf(pKey, "pad%.1fx%.1f", RTPC_Pad_W, RTPC_Pad_Z);
#endif
  system(Form("mkdir -p Graph_%s",pKey));
  system("rm -f Graph");
  system(Form("ln -sf Graph_%s Graph",pKey));

  const double deg = atan(1.0)/45.;

  Proton = new track0();
  if(gROOT->FindObject("track1"))
    {
      Electron = new track1();
    }
  Config = new config();

  //const double kMassPr=0.9383;

  /////////////////////////////////////////
  //tree variables
  int    ThrownIndex, Index;
  int    Pid;
  double Beam, Ei;

  double X0,Y0,Z0;

  //electron
  double P0_e,Theta0_e,Phi0_e;
  double Xvb_e,Yvb_e,Zvb_e;
  double Pvb_e,Thetavb_e,Phivb_e;

  double XS_e;

  double Theta0_rec_e,Phi0_rec_e,P0_rec_e;
  double X0_rec_e,Y0_rec_e,Z0_rec_e;


  //PROTON
  double P0_p,Theta0_p,Phi0_p;
  double Xvb_p,Yvb_p,Zvb_p;
  double Pvb_p,Thetavb_p,Phivb_p;

  double Theta0_rec_p,Phi0_rec_p,P0_rec_p;
  double X0_rec_p,Y0_rec_p,Z0_rec_p;

  double  Smin,Smax;
  double  MeanBz, dEdX;
  int     HitNum;
  double  StepX[MaxHit];
  double  StepY[MaxHit];
  double  StepZ[MaxHit];
  double  StepS[MaxHit];
  double  StepPhi[MaxHit];
  double  StepdE[MaxHit];
  double  StepL[MaxHit];
  double  R_sim, A_sim, B_sim, Theta_sim, Phi_sim, Z_sim, DCA_sim; //fit result to sim.track

  int     StepID[MaxHit];
  int     StepTDC[MaxHit];
  int     StepADC[MaxHit];
  int     ShiftTDC=0;

  int     HitNum_rec=0;
  double  StepX_rec[MaxHit_rec];
  double  StepY_rec[MaxHit_rec];
  double  StepZ_rec[MaxHit_rec];
  double  StepS_rec[MaxHit_rec];
  double  StepPhi_rec[MaxHit_rec];
  int     StepID_rec[MaxHit_rec];
  int     StepTDC_rec[MaxHit_rec];
  int     StepADC_rec[MaxHit_rec];

  int     HitNum_m=0;
  int     StepID_m[MaxHit];
  int     StepTDC_m[MaxHit];
  int     StepADC_m[MaxHit];
  double  StepX_rec_m[MaxHit];
  double  StepY_rec_m[MaxHit];
  double  StepZ_rec_m[MaxHit];
  double  StepS_rec_m[MaxHit];
  double  StepPhi_rec_m[MaxHit];
  double  R_rec, A_rec, B_rec, Theta_rec, Phi_rec, Z_rec, DCA_rec; //fit result to sim. track

  ///////////////////////////////////////////////////////////////////////
  //read config
  Config->GetEntry(0);
  Beam=Config->Beam;

  ///////////////////////////////////////////////////////////////////////

  char pOutAccFileName[255];
  sprintf(pOutAccFileName, "Graph/RTPCAcc.root");
  TFile *pFileAcc = new TFile(pOutAccFileName,"RECREATE");
  pFileAcc->cd();

  int NP=55, NTh_tr=35, NPh_tr=36;
  double P_min=0.0, P_max=11.0;
  double Th_min=5*deg, Th_max=40*deg;
  double Ph_min=-180*deg, Ph_max=180*deg;


  TH3F *h3N = new TH3F("h3N","Detected e^{-};#phi;#theta;P/GeV",
                   NPh_tr,Ph_min,Ph_max,NTh_tr,Th_min,Th_max,NP,P_min,P_max);
  TH3F *h3D = new TH3F("h3D","Thrown e^{-};#phi;#theta;P/GeV",
                   NPh_tr,Ph_min,Ph_max,NTh_tr,Th_min,Th_max,NP,P_min,P_max);

  TH2F *h2N_PT = new TH2F("h2N_PT","Detected e^{-};#theta;P/GeV",
                    NTh_tr,Th_min,Th_max,NP,P_min,P_max);
  TH2F *h2D_PT = new TH2F("h2D_PT","Thrown e^{-};#theta;P/GeV",
                    NTh_tr,Th_min,Th_max,NP,P_min,P_max);

  TH2F *h2N_TF = new TH2F("h2N_TF","Detected e^{-};#phi;#theta",
                    NPh_tr,Ph_min,Ph_max,NTh_tr,Th_min,Th_max);
  TH2F *h2D_TF = new TH2F("h2D_TF","Thrown e^{-};#phi;#theta",
                    NPh_tr,Ph_min,Ph_max,NTh_tr,Th_min,Th_max);

  TH2F *h2N_PF = new TH2F("h2N_PF","Detected e^{-};#phi;P/GeV",
                    NPh_tr,Ph_min,Ph_max,NP,P_min,P_max);
  TH2F *h2D_PF = new TH2F("h2D_PF","Thrown e^{-};#phi;P/GeV",
                    NPh_tr,Ph_min,Ph_max,NP,P_min,P_max);


  NP=30, NTh_tr=35, NPh_tr=36;
  P_min=0.05, P_max=0.35;
  Th_min=20*deg, Th_max=160*deg;
  Ph_min=-180*deg, Ph_max=180*deg;

  TH3F *h3N_p = new TH3F("h3N_p","Detected p_{RTPC};#phi;#theta;P/GeV",
                   NPh_tr,Ph_min,Ph_max,NTh_tr,Th_min,Th_max,NP,P_min,P_max);
  TH3F *h3D_p = new TH3F("h3D_p","Thrown p_{RTPC};#phi;#theta;P/GeV",
                   NPh_tr,Ph_min,Ph_max,NTh_tr,Th_min,Th_max,NP,P_min,P_max);

  TH2F *h2N_PT_p = new TH2F("h2N_PT_p","Detected p_{RTPC};#theta;P/GeV",
                      NTh_tr,Th_min,Th_max,NP,P_min,P_max);
  TH2F *h2D_PT_p = new TH2F("h2D_PT_p","Thrown p_{RTPC};#theta;P/GeV",
                      NTh_tr,Th_min,Th_max,NP,P_min,P_max);

  TH2F *h2N_TF_p = new TH2F("h2N_TF_p","Detected p_{RTPC};#phi;#theta",
                      NPh_tr,Ph_min,Ph_max,NTh_tr,Th_min,Th_max);
  TH2F *h2D_TF_p = new TH2F("h2D_TF_p","Thrown p_{RTPC};#phi;#theta",
                      NPh_tr,Ph_min,Ph_max,NTh_tr,Th_min,Th_max);

  TH2F *h2N_PF_p = new TH2F("h2N_PF_p","Detected p_{RTPC};#phi;P/GeV",
                      NPh_tr,Ph_min,Ph_max,NP,P_min,P_max);
  TH2F *h2D_PF_p = new TH2F("h2D_PF_p","Thrown p_{RTPC};#phi;P/GeV",
                      NPh_tr,Ph_min,Ph_max,NP,P_min,P_max);


  ///////////////////////////////////////////////////////////////////////

  char pOutFileName[255];
  if(strlen(outfile)<5)
    sprintf(pOutFileName, "%s_%s.root",GetBaseName(infile),pKey);
  else
    sprintf(pOutFileName, "%s",outfile);
  //sprintf(pOutFileName, "%s_%s.root",basename(outfile),pKey);

  cout<<"outfile="<<pOutFileName<<"  key="<<pKey<<endl;

  TFile *pFile=new TFile(pOutFileName,"RECREATE");
  pFile->cd();

  TH2F *h2Thre_C = new TH2F("h2Thre_C","Barely Reach Drift Region;#theta-90 (deg);P (MeV/c)",
                      180,-90,90,200,50,250);
  TH2F *h2Thre_A= new TH2F("h2Thre_A","Barely Reach 1st GEM Foil;#theta-90 (deg);P (MeV/c)",
                     180,-90,90,200,50,250);

  TTree *pTree=new TTree("ep","tagged DIS events");

  //////////////////////////////////////////////////////////////////////
  pTree->Branch("ThrownIndex",&ThrownIndex,"ThrownIndex/I");
  pTree->Branch("Index",&Index,"Index/I");
  pTree->Branch("Pid",&Pid,"Pid/I");

  pTree->Branch("Beam",&Beam,"Beam/D");
  pTree->Branch("Ei",&Ei,"Ei/D");

  pTree->Branch("X0",&X0,"X0/D");
  pTree->Branch("Y0",&Y0,"Y0/D");
  pTree->Branch("Z0",&Z0,"Z0/D");

  //electron which match to the proton,
  //this is used to study the electron background
  pTree->Branch("P0_e",&P0_e,"P0_e/D");
  pTree->Branch("Theta0_e",&Theta0_e,"Theta0_e/D");
  pTree->Branch("Phi0_e",&Phi0_e,"Phi0_e/D");
  pTree->Branch("Xvb_e",&Xvb_e,"Xvb_e/D");
  pTree->Branch("Yvb_e",&Yvb_e,"Yvb_e/D");
  pTree->Branch("Zvb_e",&Zvb_e,"Zvb_e/D");
  pTree->Branch("Pvb_e",&Pvb_e,"Pvb_e/D");
  pTree->Branch("Thetavb_e",&Thetavb_e,"Thetavb_e/D");
  pTree->Branch("Phivb_e",&Phivb_e,"Phivb_e/D");

  pTree->Branch("XS_e",&XS_e,"XS_e/D");

  pTree->Branch("Theta0_rec_e",&Theta0_rec_e,"Theta0_rec_e/D");
  pTree->Branch("Phi0_rec_e",&Phi0_rec_e,"Phi0_rec_e/D");
  pTree->Branch("P0_rec_e",&P0_rec_e,"P0_rec_e/D");

  pTree->Branch("X0_rec_e",&X0_rec_e,"X0_rec_e/D");
  pTree->Branch("Y0_rec_e",&Y0_rec_e,"Y0_rec_e/D");
  pTree->Branch("Z0_rec_e",&Z0_rec_e,"Z0_rec_e/D");


  //////////////////////////////////////////////////////////////////////
  //proton
  pTree->Branch("P0_p",&P0_p,"P0_p/D");
  pTree->Branch("Theta0_p",&Theta0_p,"Theta0_p/D");
  pTree->Branch("Phi0_p",&Phi0_p,"Phi0_p/D");
  pTree->Branch("Xvb_p",&Xvb_p,"Xvb_p/D");
  pTree->Branch("Yvb_p",&Yvb_p,"Yvb_p/D");
  pTree->Branch("Zvb_p",&Zvb_p,"Zvb_p/D");
  pTree->Branch("Pvb_p",&Pvb_p,"Pvb_p/D");
  pTree->Branch("Thetavb_p",&Thetavb_p,"Thetavb_p/D");
  pTree->Branch("Phivb_p",&Phivb_p,"Phivb_p/D");

  pTree->Branch("Theta0_rec_p",&Theta0_rec_p,"Theta0_rec_p/D");
  pTree->Branch("Phi0_rec_p",&Phi0_rec_p,"Phi0_rec_p/D");
  pTree->Branch("P0_rec_p",&P0_rec_p,"P0_rec_p/D");

  pTree->Branch("X0_rec_p",&X0_rec_p,"X0_rec_p/D");
  pTree->Branch("Y0_rec_p",&Y0_rec_p,"Y0_rec_p/D");
  pTree->Branch("Z0_rec_p",&Z0_rec_p,"Z0_rec_p/D");


  pTree->Branch("Smin",&Smin,"Smin/D");
  pTree->Branch("Smax",&Smax,"Smax/D");
  pTree->Branch("MeanBz",&MeanBz,"MeanBz/D");
  pTree->Branch("dEdX",&dEdX,"dEdX/D");

  pTree->Branch("HitNum",&HitNum,"HitNum/I");
  pTree->Branch("StepX",&StepX[0],"StepX[HitNum]/D");
  pTree->Branch("StepY",&StepY[0],"StepY[HitNum]/D");
  pTree->Branch("StepZ",&StepZ[0],"StepZ[HitNum]/D");
  pTree->Branch("StepS",&StepS[0],"StepS[HitNum]/D");
  pTree->Branch("StepPhi",&StepPhi[0],"StepPhi[HitNum]/D");
  pTree->Branch("StepdE",&StepdE[0],"StepdE[HitNum]/D");
  pTree->Branch("StepL",&StepL[0],"StepL[HitNum]/D");
  pTree->Branch("StepID",&StepID[0],"StepID[HitNum]/I");
  pTree->Branch("StepTDC",&StepTDC[0],"StepTDC[HitNum]/I");
  pTree->Branch("StepADC",&StepADC[0],"StepADC[HitNum]/I");
  pTree->Branch("ShiftTDC",&ShiftTDC,"ShiftTDC/I");

#if defined Store_Rec_Leaves && Store_Rec_Leaves>0
  //after hit recon
  pTree->Branch("HitNum_rec",&HitNum_rec,"HitNum_rec/I");
  pTree->Branch("StepX_rec",&StepX_rec[0],"StepX_rec[HitNum_rec]/D");
  pTree->Branch("StepY_rec",&StepY_rec[0],"StepY_rec[HitNum_rec]/D");
  pTree->Branch("StepZ_rec",&StepZ_rec[0],"StepZ_rec[HitNum_rec]/D");
  pTree->Branch("StepS_rec",&StepS_rec[0],"StepS_rec[HitNum_rec]/D");
  pTree->Branch("StepPhi_rec",&StepPhi_rec[0],"StepPhi_rec[HitNum_rec]/D");
  pTree->Branch("StepID_rec",&StepID_rec[0],"StepID_rec[HitNum_rec]/I");
  pTree->Branch("StepTDC_rec",&StepTDC_rec[0],"StepTDC_rec[HitNum_rec]/I");
  pTree->Branch("StepADC_rec",&StepADC_rec[0],"StepADC_rec[HitNum_rec]/I");
#endif

  pTree->Branch("R_sim",&R_sim,"R_sim/D");
  pTree->Branch("A_sim",&A_sim,"A_sim/D");
  pTree->Branch("B_sim",&B_sim,"B_sim/D");
  pTree->Branch("Theta_sim",&Theta_sim,"Theta_sim/D");
  pTree->Branch("Phi_sim",&Phi_sim,"Phi_sim/D");
  pTree->Branch("Z_sim",&Z_sim,"Z_sim/D");
  pTree->Branch("DCA_sim",&DCA_sim,"DCA_sim/D");

  //after hit merge
  pTree->Branch("HitNum_m",&HitNum_m,"HitNum_m/I");
  pTree->Branch("StepID_m",&StepID_m[0],"StepID_m[HitNum_m]/I");
  pTree->Branch("StepTDC_m",&StepTDC_m[0],"StepTDC_m[HitNum_m]/I");
  pTree->Branch("StepADC_m",&StepADC_m[0],"StepADC_m[HitNum_m]/I");
  pTree->Branch("StepX_rec_m",&StepX_rec_m[0],"StepX_rec_m[HitNum_m]/D");
  pTree->Branch("StepY_rec_m",&StepY_rec_m[0],"StepY_rec_m[HitNum_m]/D");
  pTree->Branch("StepZ_rec_m",&StepZ_rec_m[0],"StepZ_rec_m[HitNum_m]/D");
  pTree->Branch("StepS_rec_m",&StepS_rec_m[0],"StepS_rec_m[HitNum_m]/D");
  pTree->Branch("StepPhi_rec_m",&StepPhi_rec_m[0],"StepPhi_rec_m[HitNum_m]/D");
  pTree->Branch("R_rec",&R_rec,"R_rec/D");
  pTree->Branch("A_rec",&A_rec,"A_rec/D");
  pTree->Branch("B_rec",&B_rec,"B_rec/D");
  pTree->Branch("Theta_rec",&Theta_rec,"Theta_rec/D");
  pTree->Branch("Phi_rec",&Phi_rec,"Phi_rec/D");
  pTree->Branch("Z_rec",&Z_rec,"Z_rec/D");
  pTree->Branch("DCA_rec",&DCA_rec,"DCA_rec/D");


  //////////////////////////////////////////////////////////////////////
  //Resolution for NPS and LAC
  double pResBPM = 1.0; //mm
  //CLAS12 resolution: dp/p<1%, dtheta=1mrad, dphi=1mrad/sinTh
  double pResP_e=0.01;
  double pResTh_e=0.001;
  double pResPh_e=0.001;   //need to divided by sinTh
  double pResZ_e=4.0;

  //RTPC resolution
  //in Bonus6, the uncertainty of S is 0.35 mm, which is corresponding to 114/2 ns
  //in Bonus12, the time window is 25 ns, therefore dS=0.35/114*25=0.077mm
  //Assuming the readout pad is 4.5(phi)x5(z) mm, located at S=90mm,
  //then dPhi = RTPC_Pad_W/sqrt(12)/(2*PI*R) *(2*PI) = RTPC_Pad_W/sqrt(12)/R

  //I found that I need to make it a factor of 1.7 in order to match the resolution
  double pResS_p=1.7*0.35*RTPC_TDC_Window/114;
  double pResPh_p=RTPC_Pad_W/sqrt(12.)/RTPC_ReadOut_R;
  double pResZ_p=RTPC_Pad_Z/sqrt(12.);
  //////////////////////////////////////////////////////////////////////

  ofstream foutP;
  ofstream foutTh;
  ofstream foutPh;
  ofstream foutZ;
  char pMDFOutFileName[255];
  sprintf(pMDFOutFileName, "%s_%s_P_infile.dat",GetBaseName(infile),pKey);
  foutP.open(pMDFOutFileName);
  sprintf(pMDFOutFileName, "%s_%s_Th_infile.dat",GetBaseName(infile),pKey);
  foutTh.open(pMDFOutFileName);
  sprintf(pMDFOutFileName, "%s_%s_Ph_infile.dat",GetBaseName(infile),pKey);
  foutPh.open(pMDFOutFileName);
  sprintf(pMDFOutFileName, "%s_%s_Z_infile.dat",GetBaseName(infile),pKey);
  foutZ.open(pMDFOutFileName);

  foutP <<"     Pperp_tg       R_helix   cosTh_helix      Ph_helix       Z_helix     dca_helix        Weight"<<endl;
  foutP <<"Pperp_tg(GeV)   R_helix(mm) cos(Th_helix) Ph_helix(rad)   Z_helix(mm) dca_helix(mm)        Weight"<<endl;

  foutTh<<"     cosTh_tg       R_helix   cosTh_helix       Ph_helix      Z_helix     dca_helix        Weight"<<endl;
  foutTh<<"   cos(Th_tg)   R_helix(mm) cos(Th_helix)  Ph_helix(rad)  Z_helix(mm) dca_helix(mm)        Weight"<<endl;

  foutPh<<"        Ph_tg       R_helix   cosTh_helix      Ph_helix       Z_helix     dca_helix        Weight"<<endl;
  foutPh<<"   Ph_tg(rad)   R_helix(mm) cos(Th_helix)  Ph_helix(rad)  Z_helix(mm) dca_helix(mm)        Weight"<<endl;

  foutZ <<"         Z_tg       R_helix   cosTh_helix      Ph_helix       Z_helix     dca_helix        Weight"<<endl;
  foutZ <<"     Z_tg(mm)   R_helix(mm) cos(Th_helix)  Ph_helix(rad)  Z_helix(mm) dca_helix(mm)        Weight"<<endl;

  //foutP.setf(ios::scientific,ios::floatfield);
  //foutTh.setf(ios::scientific,ios::floatfield);
  //foutPh.setf(ios::scientific,ios::floatfield);
  //foutZ.setf(ios::scientific,ios::floatfield);
  //////////////////////////////////////////////////////////////////////

  Index=0;
  int pFoundGoodTrackInSuperEvent = 0;
  Long64_t nentries = Proton->fChain->GetEntriesFast();
  Long64_t nb0 = 0, nb1 = 0;
  for (Long64_t i=0; i<nentries;i++)
    //for (Long64_t i=0; i<5;i++)
  {
    if(!((i+1)%1000))
      printf("processing event %6d / %6d \r",int(i+1),int(nentries));

    //do proton first since its acceptance is smaller
    //read proton
    nb1 = Proton->fChain->GetEntry(i);
    if(nb1<=0) break;
    //apply cuts
    if ( Proton->Cut(i) < 0) continue;


    //read electron
    if(Electron)
    {
      nb0 = Electron->fChain->GetEntry(i);
      if(nb0<=0) break;
      //apply cuts
      if (Electron->Cut(i) < 0) continue;
      //apply trigger cuts
      if( Electron->Pvb < 0.3 )
      {
        //continue;
      }
    }

    ////////////////////////////////////////////////////////////
    //load variables
    X0=Proton->X0;
    Y0=Proton->Y0;
    Z0=Proton->Z0;
    Ei=Proton->Ei;

    ThrownIndex=Proton->Index;
    Pid=Proton->PdgId;

    //proton
    P0_p=Proton->P0;
    Theta0_p=Proton->Theta0;
    Phi0_p=Proton->Phi0;
    Xvb_p=Proton->Xvb;
    Yvb_p=Proton->Yvb;
    Zvb_p=Proton->Zvb;
    Pvb_p=Proton->Pvb;
    Thetavb_p=Proton->Thetavb;
    Phivb_p=Proton->Phivb;


    //by jixie: now I want to shift some tracks to create super events
    //one half will be shifted as 'early' the other half will be shifted 'late'
    //the ShiftTDC will be uniformly distributed within [1, pMaxTDC-1]
    if(!(Index%ntrack_per_event)) pFoundGoodTrackInSuperEvent = 0;
    if(ntrack_per_event>1)
    {
      //maximum TDC value
      int pMaxTDC = int(7000/NS_PER_TIC)-1;
      // how many TIC should be shifted in average
      double tmpN = pMaxTDC * 2.0 / ntrack_per_event;
      //how many extra TIC this track should be shifted? a value between 1 and tmpN,
      int tmpShift = rand() % int(ceil(tmpN)) + 1;

      int tureTrackIndex = int(ntrack_per_event/2.0) ;
      if((Index%ntrack_per_event)>=tureTrackIndex && Proton->StepNum>30 &&
         Proton->StepTL[Proton->StepNum-1]>RTPC_Anode_R-10.0 && !pFoundGoodTrackInSuperEvent )
      {
        ShiftTDC = 0;
        pFoundGoodTrackInSuperEvent = 1;
      }
      else ShiftTDC = -pMaxTDC + int((Index%ntrack_per_event) * tmpN) + tmpShift;
      if(ShiftTDC>pMaxTDC-1) ShiftTDC=pMaxTDC-1;
      cout<<"ThrownIndex="<<ThrownIndex<<"  Index="<<Index<<"  shiftTDC="<<ShiftTDC<<endl;
    }

    Smin=9999.;Smax=-9999.;
    HitNum=HitNum_rec=0;
    double SumBz=0, SumdEdX=0;
    for(int ss=0;ss<Proton->StepNum;ss++)
    {
      //only take Drift region
      double tmpS = sqrt(Proton->StepX[ss]*Proton->StepX[ss]+
                     Proton->StepY[ss]*Proton->StepY[ss]);
      if(tmpS>RTPC_Cathode_R && tmpS<RTPC_Anode_R)
      {
        if(Smin>tmpS) Smin=tmpS;
        else if(Smax<tmpS) Smax=tmpS;

        SumBz += Proton->StepBz[ss];
        if(Proton->StepL[ss]) SumdEdX += Proton->StepdE[ss]/Proton->StepL[ss];

        double tmpPhi=atan2(Proton->StepY[ss],Proton->StepX[ss]);
        StepX[HitNum]=Proton->StepX[ss];
        StepY[HitNum]=Proton->StepY[ss];
        StepZ[HitNum]=Proton->StepZ[ss];
        StepS[HitNum]=tmpS;
        StepPhi[HitNum]=tmpPhi;
        StepdE[HitNum]=Proton->StepdE[ss];
        StepL[HitNum]=Proton->StepL[ss];

        bool bUseMagboltz=true;
        if(!bUseMagboltz)
        {
          //Now mimic the digitization and reconstruction of RTPC track
          //by smearing real hit location
          double tmpS_rec = tmpS + fGaus(0,pResS_p);
          double tmpPhi_rec = tmpPhi + fGaus(0,pResPh_p);
          double tmpZ_rec = StepZ[HitNum] + fGaus(0,pResZ_p);

          StepX_rec[HitNum]=tmpS_rec*cos(tmpPhi_rec);
          StepY_rec[HitNum]=tmpS_rec*sin(tmpPhi_rec);
          StepZ_rec[HitNum]=tmpZ_rec;
          StepS_rec[HitNum]=tmpS_rec;
          StepPhi_rec[HitNum]=tmpPhi_rec;
        }
        else
        {
          //now use magboltz drift path
          double xi=StepX[HitNum], yi=StepY[HitNum], zi=StepZ[HitNum];
          double dE=StepdE[HitNum];

          int n_rec=0;
          double xo[4],yo[4],zo[4];
          int chan[4],tdc[4],adc[4];
          gEsim->DriftESim(xi,yi,zi,dE,n_rec,xo,yo,zo,chan,tdc,adc);

          StepID[HitNum]=chan[0];
          StepTDC[HitNum]=tdc[0];
          StepADC[HitNum]=adc[0];
          for(int t=1;t<n_rec;t++) StepADC[HitNum] += adc[t];

          for(int t=0;chan[t]>0 && t<n_rec;t++)
          {
            //by jixie: now I want to shift some tracks to create super events
            if(ntrack_per_event>1)
            {
              tdc[t] += ShiftTDC*NS_PER_TIC;
              if(tdc[t]<0) continue;
              gEsim->LookupXYZByIDTDC(chan[t],tdc[t],xo[t],yo[t],zo[t]);
              double tmpSS = sqrt(xo[t]*xo[t]+yo[t]*yo[t]);
              if(tmpSS<RTPC_Cathode_R-10.0 || tmpSS>RTPC_Anode_R+10.0) continue;
            }

            StepID_rec[HitNum_rec]=chan[t];
            StepTDC_rec[HitNum_rec]=tdc[t];
            StepADC_rec[HitNum_rec]=adc[t];
            StepX_rec[HitNum_rec]=xo[t];
            StepY_rec[HitNum_rec]=yo[t];
            StepZ_rec[HitNum_rec]=zo[t];
            StepPhi_rec[HitNum_rec]=atan2(yo[t],xo[t]);
            StepS_rec[HitNum_rec]=sqrt(xo[t]*xo[t]+yo[t]*yo[t]);

            HitNum_rec++;
            if(HitNum_rec>=MaxHit_rec) break;
          }
        }

        HitNum++;
        if(HitNum>=MaxHit) break;
      }
    }
    MeanBz = (HitNum>0) ? SumBz/HitNum : 0;
    dEdX = (HitNum>0) ? SumdEdX/HitNum : 0;

    ////////////////////////////////////////////////////////////////////////
    //merge step if ID and TDC both equal to previous step
    HitNum_m=0;
    for(int ss=0;StepID_rec[ss]>=0 && ss<HitNum_rec;ss++)
    {
      int found=0;
      for(int tt=HitNum_m-1;tt>=0;tt--)
      {
        if(StepID_rec[ss]==StepID_m[tt] && StepTDC_rec[ss]==StepTDC_m[tt])
        {
          StepADC_m[tt]+=StepADC_rec[ss];
          found=1;
          //cout<<"\t Merge ADC of _rec point="<<ss<<" to _m point="<<tt<<", HitNum_m="<<HitNum_m<<endl;
          break;
        }
      }
      if(found) continue;
      //keep only those valid reconstruction
      StepID_m[HitNum_m]=StepID_rec[ss];
      StepTDC_m[HitNum_m]=StepTDC_rec[ss];
      StepADC_m[HitNum_m]=StepADC_rec[ss];
      StepX_rec_m[HitNum_m]=StepX_rec[ss];
      StepY_rec_m[HitNum_m]=StepY_rec[ss];
      StepZ_rec_m[HitNum_m]=StepZ_rec[ss];
      StepPhi_rec_m[HitNum_m]=StepPhi_rec[ss];
      StepS_rec_m[HitNum_m]= StepS_rec[ss];

      HitNum_m++;
      if(HitNum_m>=MaxHit) break;
    }

    //do the helix fit and reconstruction
    if(HitNum_m>5)
    {
      FitATrack(HitNum,StepX,StepY,StepZ,R_sim,A_sim,B_sim,Phi_sim,
              Theta_sim,Z_sim,RTPC_Cathode_R+0.1,RTPC_Anode_R-0.1);
      DCA_sim=sqrt(A_sim*A_sim+B_sim*B_sim)-R_sim;
      FitATrack(HitNum_m,StepX_rec_m,StepY_rec_m,StepZ_rec_m,R_rec,A_rec,B_rec,
              Phi_rec,Theta_rec,Z_rec,RTPC_Cathode_R+0.1,RTPC_Anode_R-0.1);
      DCA_rec=sqrt(A_rec*A_rec+B_rec*B_rec)-R_rec;
    }

    double P_corr,Theta_corr,Phi_corr;
    RTPC_Recon( MeanBz, R_rec, Theta_rec, Phi_rec, P_corr, Theta_corr,Phi_corr);
    Theta0_rec_p=Theta_corr;
    Phi0_rec_p=Phi_corr;
    P0_rec_p=P_corr;
    X0_rec_p=X0 + fGaus(0.0,pResBPM);
    Y0_rec_p=Y0 + fGaus(0.0,pResBPM);
    Z0_rec_p=Z_rec;

    /////////////////////////////////////////////////
    //get thrown proton information

    h3D_p->Fill(Phi0_p,Theta0_p,P0_p);
    h2D_PT_p->Fill(Theta0_p,P0_p);
    h2D_TF_p->Fill(Phi0_p,Theta0_p);
    h2D_PF_p->Fill(Phi0_p,P0_p);
    //get reconstructed information
    if(fabs(1.0-P0_rec_p/P0_p)<0.1 && fabs(Z0-Z0_rec_p)<15 &&
     fabs(Theta0_p-Theta0_rec_p)<0.05 && fabs(Phi0_p-Phi0_rec_p)<0.05)
    {
      h3N_p->Fill(Phi0_p,Theta0_p,P0_p);
      h2N_PT_p->Fill(Theta0_p,P0_p);
      h2N_TF_p->Fill(Phi0_p,Theta0_p);
      h2N_PF_p->Fill(Phi0_p,P0_p);
    }

    //fill threshold histogram
    double Theta_p_min = atan((RTPC_Cathode_R+5)/(RTPC_Length/2-Z0));
    if(Smax<RTPC_Cathode_R+10 && Theta0_p>Theta_p_min && Smax>StepS[HitNum-1])
    {
      h2Thre_C->Fill(Theta0_p/deg-90,P0_p*1000);
    }
    if(Smax>RTPC_Anode_R-5 && Smax>StepS[HitNum-1])
    {
      h2Thre_A->Fill(Theta0_p/deg-90,P0_p*1000);
    }

    //////////////////////////////////////////////////////////
    //electron
    if(Electron)
    {
      P0_e=Electron->P0;
      Theta0_e=Electron->Theta0;
      Phi0_e=Electron->Phi0;
      Xvb_e=Electron->Xvb;
      Yvb_e=Electron->Yvb;
      Zvb_e=Electron->Zvb;
      Pvb_e=Electron->Pvb;
      Thetavb_e=Electron->Thetavb;
      Phivb_e=Electron->Phivb;

      XS_e = 0;
      if(Electron->Pvb>0.3)
      {
        //need to call P.Bosted XS
        XS_e = Electron->ElasXS;
      }

      //reconstruct the electron
      Theta0_rec_e=Theta0_e+fGaus(0.0,pResTh_e);
      Phi0_rec_e=Phi0_e+fGaus(0.0,pResPh_e/sin(Theta0_e));
      P0_rec_e=P0_e*fGaus(1.0,pResP_e);
      X0_rec_e = X0 + fGaus(0.0,pResBPM);
      Y0_rec_e = Y0 + fGaus(0.0,pResBPM);
      Z0_rec_e = Z0 + fGaus(0.0,pResZ_e);


      //get thrown information
      h3D->Fill(Phi0_e,Theta0_e,P0_e);
      h2D_PT->Fill(Theta0_e,P0_e);
      h2D_TF->Fill(Phi0_e,Theta0_e);
      h2D_PF->Fill(Phi0_e,P0_e);
      //get reconstructed information
      if(Electron->Pvb>0.3)
      {
        h3N->Fill(Phi0_e,Theta0_e,P0_e);
        h2N_PT->Fill(Theta0_e,P0_e);
        h2N_TF->Fill(Phi0_e,Theta0_e);
        h2N_PF->Fill(Phi0_e,P0_e);
      }
    }

    ///////////////////////////////////////////////
    //write out raw date file to run MulDFit
    if(HitNum>10 && fabs(1.0-P0_rec_p/P0_p)<0.1 && fabs(Z0-Z0_rec_p)<15 &&
     fabs(Theta0_p-Theta0_rec_p)<0.05 && fabs(Phi0_p-Phi0_rec_p)<0.05 &&
     DCA_sim<3)
    {
      double weight = 1.0/(fabs(sqrt(A_sim*A_sim+B_sim*B_sim)-R_sim)/0.2+0.1);
      weight=1.0;
      foutP<<setw(13)<<P0_p*sin(Theta0_p)<<" "
           <<setw(13)<<R_rec<<" "<<setw(13)<<cos(Theta_rec)<<" "<<setw(13)<<Phi_rec<<" "
           <<setw(13)<<Z_rec<<" "<<setw(13)<<DCA_rec<<" "<<setw(13)<<weight<<endl;

      foutTh<<setw(13)<<cos(Theta0_p)<<" "
          <<setw(13)<<R_rec<<" "<<setw(13)<<cos(Theta_rec)<<" "<<setw(13)<<Phi_rec<<" "
          <<setw(13)<<Z_rec<<" "<<setw(13)<<DCA_rec<<" "<<setw(13)<<weight<<endl;

      foutPh<<setw(13)<<Phi0_p<<" "
          <<setw(13)<<R_rec<<" "<<setw(13)<<cos(Theta_rec)<<" "<<setw(13)<<Phi_rec<<" "
          <<setw(13)<<Z_rec<<" "<<setw(13)<<DCA_rec<<" "<<setw(13)<<weight<<endl;

      foutZ<<setw(13)<<Z0<<" "
           <<setw(13)<<R_rec<<" "<<setw(13)<<cos(Theta_rec)<<" "<<setw(13)<<Phi_rec<<" "
           <<setw(13)<<Z_rec<<" "<<setw(13)<<DCA_rec<<" "<<setw(13)<<weight<<endl;
    }

    ///////////////////////////////////////////////

    //if(HitNum>5 && DCA_sim<0.2) {pTree->Fill(); Index++;}
    pTree->Fill(); Index++;
    /////////////////////////////////////////////////////////////

    //reset
    P0_e=Theta0_e=Phi0_e=-99.0;
    Xvb_e=Yvb_e=Zvb_e=-99.0;
    Pvb_e=Thetavb_e=Phivb_e=-99.0;

    P0_p=Theta0_p=Phi0_p=-99.0;
    Xvb_p=Yvb_p=Zvb_p=-99.0;
    Pvb_p=Thetavb_p=Phivb_p=-99.0;

  }

  TCanvas *c14 = new TCanvas("c14","RTPC Pr Threshold",800,600);
  c14->cd(0);
  h2Thre_C->Draw("");
  h2Thre_C->Fit("pol6");
  c14->Modified();
  c14->SaveAs(Form("Graph/PrThre_Cathode_%s.png","RTPC12"));

  TCanvas *c4 = new TCanvas("c4","RTPC Pr Threshold",800,600);
  c4->cd(0);
  h2Thre_A->Draw("*");
  h2Thre_A->Fit("pol6");
  c4->Modified();
  c4->SaveAs(Form("Graph/PrThre_Anode_%s.png","RTPC12"));


  pFile->Write("", TObject::kOverwrite);
  pFile->Close();

  foutP.close();
  foutTh.close();
  foutPh.close();
  foutZ.close();

  ////////////////////////////////////
  pFileAcc->cd();

  if(Electron)
    {
      TCanvas *c2 = new TCanvas("c2","e- acceptance",900,700);
      TH3F *h3Acc = (TH3F*) h3N->Clone("h3Acc_PTF");
      h3Acc->SetTitle("e^{-} Acceptance");
      h3Acc->Divide(h3D);
      h3Acc->Scale(100.0);

      TH2F *h2Acc_PT = (TH2F*) h2N_PT->Clone("h2Acc_PT");
      h2Acc_PT->SetTitle("e^{-} Acceptance");
      h2Acc_PT->Divide(h2D_PT);
      h2Acc_PT->Scale(100.0);

      TH2F *h2Acc_TF = (TH2F*) h2N_TF->Clone("h2Acc_TF");
      h2Acc_TF->SetTitle("e^{-} Acceptance");
      h2Acc_TF->Divide(h2D_TF);
      h2Acc_TF->Scale(100.0);

      TH2F *h2Acc_PF = (TH2F*) h2N_PF->Clone("h2Acc_PF");
      h2Acc_PF->SetTitle("e^{-} Acceptance");
      h2Acc_PF->Divide(h2D_PF);
      h2Acc_PF->Scale(100.0);

      c2->Divide(2,2);
      c2->cd(1);
      h2Acc_PT->Draw("colz text");
      c2->cd(2);
      h2Acc_TF->Draw("colz text");
      c2->cd(3);
      h2Acc_PF->Draw("colz text");

      c2->cd(4);
      h3Acc->Draw();

      c2->Modified();
      c2->SaveAs(Form("Graph/ElAcc_%s.png","CLAS12"));

      //DoAcc(pKey,pTree);
      PrintAcc(h3Acc,h3N,3,1);
      h3Acc->Write("", TObject::kOverwrite);
    }


  TCanvas *c3 = new TCanvas("c3","e- acceptance",900,700);
  TH3F *h3Acc_p = (TH3F*) h3N_p->Clone("h3Acc_PTF_p");
  h3Acc_p->SetTitle("p_{RTPC} Acceptance");
  h3Acc_p->Divide(h3D_p);
  h3Acc_p->Scale(100.0);

  TH2F *h2Acc_PT_p = (TH2F*) h2N_PT_p->Clone("h2Acc_PT_p");
  h2Acc_PT_p->SetTitle("p_{RTPC} Acceptance");
  h2Acc_PT_p->Divide(h2D_PT_p);
  h2Acc_PT_p->Scale(100.0);

  TH2F *h2Acc_TF_p = (TH2F*) h2N_TF_p->Clone("h2Acc_TF_p");
  h2Acc_TF_p->SetTitle("p_{RTPC} Acceptance");
  h2Acc_TF_p->Divide(h2D_TF_p);
  h2Acc_TF_p->Scale(100.0);

  TH2F *h2Acc_PF_p = (TH2F*) h2N_PF_p->Clone("h2Acc_PF_p");
  h2Acc_PF_p->SetTitle("p_{RTPC} Acceptance");
  h2Acc_PF_p->Divide(h2D_PF_p);
  h2Acc_PF_p->Scale(100.0);

  c3->Divide(2,2);
  c3->cd(1);
  h2Acc_PT_p->Draw("colz text");
  c3->cd(2);
  h2Acc_TF_p->Draw("colz text");
  c3->cd(3);
  h2Acc_PF_p->Draw("colz text");

  c3->cd(4);
  h3Acc_p->Draw();

  c3->Modified();
  c3->SaveAs(Form("Graph/PrAcc_%s.png","RTPC12"));

  //DoAcc(pKey,pTree);
  PrintAcc(h3Acc_p,h3N_p,3,1);
  h3Acc_p->Write("", TObject::kOverwrite);


  pFileAcc->Write("", TObject::kOverwrite);
  pFileAcc->Close();

  printf("\n\n");

}
