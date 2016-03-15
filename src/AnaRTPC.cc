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

//#define RTPC_BoNuS6 1

double RTPC_Cathode_R = 30;		//mm
double RTPC_Anode_R   = 80;		//mm
double RTPC_ReadOut_R = 90;		//mm
double RTPC_TDC_Window= 25;		//ns
double RTPC_Pad_Z = 5.0;		//mm
double RTPC_Pad_W = 4.5;		//mm
double RTPC_Length = 400;               //mm  

DriftEMagboltz *gEsim;

//estimate number of channel: assuming 30 degree missing in phi coverage 
// Row# = int(2*PI*Readout_R/RTPC_Pad_W)
// Col# = int(RTPC_Length/RTPC_Pad_Z)
int GetNofCh(int padtype)
{
  const double PI=asin(1.0)*2;
  int iRow = int(ceil(2*PI*RTPC_ReadOut_R/RTPC_Pad_W * 350./360.));
  int iCol = int(ceil(RTPC_Length/RTPC_Pad_Z));
  int N = iRow * iCol;
  if(padtype==3)
  {
    int iRow_zstrip = int(ceil(2*PI*RTPC_ReadOut_R/0.4 * 350./360.));
    int iCol_zstrip = int(ceil(RTPC_Length/20.0));
    int N_zstrip = iRow_zstrip*iCol_zstrip;

    int iRow_phistrip = int(ceil(2*PI*RTPC_ReadOut_R/20.0 * 350./360.));
    int iCol_phistrip = int(ceil(RTPC_Length/0.4));
    int N_phistrip = iRow_phistrip*iCol_phistrip;
    N = N_zstrip + N_phistrip;
  }
  else if(padtype==4 || padtype==5)
  {
    int iRow_zstrip = int(ceil(2*PI*RTPC_ReadOut_R/1.0 * 350./360.));
    int iCol_zstrip = int(ceil(RTPC_Length/20.0));
    int N_zstrip = iRow_zstrip*iCol_zstrip;

    int iRow_phistrip = int(ceil(2*PI*RTPC_ReadOut_R/20.0 * 350./360.));
    int iCol_phistrip = int(ceil(RTPC_Length/1.0));
    int N_phistrip = iRow_phistrip*iCol_phistrip;
    N = N_zstrip + N_phistrip;
  }
  cout<<"GetNofCh(padtype="<<padtype<<") return "<<N<<endl;
  return N;
}

//config pad size:  
//padtype = 1)4.5x5,  2) 2x2  3) compass 2-D readout, 0.4x0.4 equivalent
//4)compass 2-D readout, 1x20,1 x 1 mm equivalent, 5) 2.5x4 mm,Sebastian's
void ConfigPadSize(int padtype)
{
  //the following will be used to estimate RTPC resolution and fiting
#if defined RTPC_BoNuS6
  RTPC_Cathode_R = 30;		//mm
  RTPC_Anode_R   = 60;		//mm
  RTPC_ReadOut_R = 70;		//mm
  RTPC_TDC_Window= 114;		//ns
  RTPC_Pad_Z = 5.0;			//mm
  RTPC_Pad_W = 4.5;			//mm
  RTPC_Length = 200;			//mm
#else
  RTPC_Cathode_R  = 30;		//mm
  RTPC_Anode_R    = 80;		//mm
  RTPC_ReadOut_R  = 90;		//mm
  RTPC_TDC_Window = 340;		//ns
  RTPC_Length = 400;			//mm

  if(padtype == 1)
  {
    RTPC_Pad_Z = 5.0;		//mm
    RTPC_Pad_W = 4.5;		//mm
  }
  else if(padtype == 2)
  {
    RTPC_Pad_Z = 2.0;		//mm
    RTPC_Pad_W = 2.0;		//mm
  }
  else if(padtype == 3)
  {
    //20mm x 0.4mm z-strip and 0.4mm x 20mm phi-strip
    RTPC_Pad_Z = 0.4;		//mm
    RTPC_Pad_W = 0.4;		//mm
  }
  else if(padtype == 4)
  {
    //20mm x 1mm z-strip and 1mm x 20mm phi-strip
    RTPC_Pad_Z = 1.0;		//mm
    RTPC_Pad_W = 1.0;		//mm
  }
  else if(padtype == 5)
  {
    //20mm x 1mm z-strip and 1mm x 20mm phi-strip
    RTPC_Pad_Z = 1.0;		//mm
    RTPC_Pad_W = 1.0;		//mm
    RTPC_Pad_Z = 20.0;		//mm
    RTPC_Pad_W = 20.0;		//mm
    RTPC_Anode_R    = 70;	//mm
    RTPC_ReadOut_R  = 80;	//mm
  }
  else if(padtype == 6)
  {
    RTPC_Pad_Z = 4.0;		//mm
    RTPC_Pad_W = 2.8;		//mm
    RTPC_Anode_R    = 70;	//mm
    RTPC_ReadOut_R  = 80;	//mm
  }
#endif
  GetNofCh(padtype);
}
static track0 *Proton=0;
static track1 *Electron=0;
static config *Config=0;

//boxmuller gauss number generator
//input: mean m, standard deviation s 
double fGaus(double m, double s);
void PrintAcc(TH3F *h3,TH3F *h3N, int Nmin=0, int TCS=0);

int FitATrack(int StepNum, double StepX[],double StepY[],double StepZ[],
  double &R, double &A, double &B,
  double &Phi_deg, double &Theta_deg, double &Z,
  double RLimit_l=30.1, double RLimit_h=79.9);

void RTPC_Recon(double MeanBz, double R_rec, double Theta_rec, double Phi_rec,
  double &P_corr, double &Theta_corr, double &Phi_corr);

const char *GetFileName(const char* str)
{
  const char *name = strrchr(str,'/');
  //printf ("\"%s\": Last occurence of '/' found at %d \n",str,int(name-str+1));
  if(name==NULL) name=&str[0];
  else name++; 	
  return name;
}


const char *GetBaseName(const char* str)
{
  const char *name = strrchr(str,'/');
  if(name==NULL) name=&str[0];
  else name++; 	
  //printf ("\"%s\": Last occurence of '/' found at %d, name=%s\n",str,int(name-str),name);

  const char *surfix = strrchr(name,'.');	
  if(surfix==NULL)  surfix=&name[sizeof(name)];
  //printf ("\"%s\": Last occurence of '.' found at %d \n",name,int(surfix-name+1));

  int length = surfix-name;
  //this method have a memory leakage
  //char *basename = new char[length+1];	
  //strncpy(basename,name,length);basename[length]='\0';
  //printf("basename=%s\n",basename);
  //return basename;

  //I am using string to avoid memory leakage
  string basename;
  basename.append(name,length); 
  //printf("basename=%s\n",basename.c_str());
  return basename.c_str();
}


void AnaRTPC(const char *infile="nt.root", const char *outfile="nt_ep.root", int padtype=2)
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

  gEsim = new DriftEMagboltz(0.9,(float)RTPC_Pad_W,(float)RTPC_Pad_Z,(float)RTPC_ReadOut_R,(float)RTPC_Length);

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

  double  StepX_rec[MaxHit];
  double  StepY_rec[MaxHit];
  double  StepZ_rec[MaxHit];
  double  StepS_rec[MaxHit];
  double  StepPhi_rec[MaxHit];

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

  pTree->Branch("StepX_rec",&StepX_rec[0],"StepX_rec[HitNum]/D");
  pTree->Branch("StepY_rec",&StepY_rec[0],"StepY_rec[HitNum]/D");
  pTree->Branch("StepZ_rec",&StepZ_rec[0],"StepZ_rec[HitNum]/D");
  pTree->Branch("StepS_rec",&StepS_rec[0],"StepS_rec[HitNum]/D");
  pTree->Branch("StepPhi_rec",&StepPhi_rec[0],"StepPhi_rec[HitNum]/D");
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
  double pResS_p=0.35*RTPC_TDC_Window/114;  
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

    Smin=9999.;Smax=-9999.;
    HitNum=0;
    double SumBz=0, SumdEdX=0;
    for(int ss=0;ss<Proton->StepNum;ss++)
    {
      //only take Drift region
      double tmpS = sqrt(Proton->StepX[ss]*Proton->StepX[ss]+
	Proton->StepY[ss]*Proton->StepY[ss]);
      if(tmpS>RTPC_Cathode_R && tmpS<RTPC_Anode_R)
      {
	if(Smin>tmpS) Smin=tmpS;
	if(Smax<tmpS) Smax=tmpS;

	SumBz += Proton->StepBz[ss];
	SumdEdX += Proton->StepdE[ss]/Proton->StepL[ss];

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
	  float xi=StepX[HitNum], yi=StepY[HitNum], zi=StepZ[HitNum];
	  float dE=StepdE[HitNum];
	  float xo=0,yo=0,zo=0;
	  int chan=-1,tdc=-1,adc=-1;
	  gEsim->DriftESim(xi,yi,zi,dE,xo,yo,zo,chan,tdc,adc);

	  StepID[HitNum]=chan;
	  StepTDC[HitNum]=tdc;
	  StepADC[HitNum]=adc;
	  StepX_rec[HitNum]=xo;
	  StepY_rec[HitNum]=yo;
	  StepZ_rec[HitNum]=zo;
	  StepPhi_rec[HitNum]=atan2(yo,xo);
	  StepS_rec[HitNum]=sqrt(xo*xo+yo*yo);

	}


	HitNum++;
	if(HitNum>=200) break;
      }
    }
    MeanBz = (HitNum>0) ? SumBz/HitNum : 0;
    dEdX = (HitNum>0) ? SumdEdX/HitNum : 0;

    ////////////////////////////////////////////////////////////////////////
    //merge step if ID and TDC both equal to previous step
    HitNum_m=0;
    for(int ss=0;ss<HitNum;ss++)
    {
      if(ss>0 && StepID[ss]==StepID[ss-1] && StepTDC[ss]==StepTDC[ss-1])  
      {
	StepADC_m[HitNum_m-1]+=StepADC[ss];
	continue;
      }
      if(StepID[ss]>=0)
      {
	//keep only those valid reconstruction
	StepID_m[HitNum_m]=StepID[ss];
	StepTDC_m[HitNum_m]=StepTDC[ss];
	StepADC_m[HitNum_m]=StepADC[ss];
	StepX_rec_m[HitNum_m]=StepX_rec[ss];
	StepY_rec_m[HitNum_m]=StepY_rec[ss];
	StepZ_rec_m[HitNum_m]=StepZ_rec[ss];
	StepPhi_rec_m[HitNum_m]=StepPhi_rec[ss];
	StepS_rec_m[HitNum_m]= StepS_rec[ss];

	HitNum_m++;
      }
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
