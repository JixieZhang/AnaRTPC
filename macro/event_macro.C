#include <iostream>
#include "string.h"
#include "math.h"

#include "TROOT.h"
#include "TFile.h"
#include "TSystem.h"
#include "TTree.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TColor.h"
#include "TStyle.h"
#include "TGeoManager.h"
#include "TGeoTube.h"
#include "TGeoMaterial.h"
#include "TGeoVolume.h"
#include "TGeoMedium.h"
#include "TPolyLine3D.h"

using namespace std;
TGeoVolume *gWorld=0;

void set_color_env(){   
  //******************************************************************
  //code to improve the color palette
  //from the root web page and color codes from:
  //http://ultrahigh.org/2007/08/20/making-pretty-root-color-palettes/
  //******************************************************************
  
  const Int_t NRGBs = 5;
  const Int_t NCont = 255;
  
  Double_t stops[NRGBs] = { 0.00, 0.34, 0.61, 0.84, 1.00 };
  Double_t red[NRGBs]   = { 0.00, 0.00, 0.87, 1.00, 0.51 };
  Double_t green[NRGBs] = { 0.00, 0.81, 1.00, 0.20, 0.00 };
  Double_t blue[NRGBs]  = { 0.51, 1.00, 0.12, 0.00, 0.00 };
  TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
  gStyle->SetNumberContours(NCont);
  

  //  gStyle->SetNumberContours(255); 
  //******************************************************************
  

  /*/
const Int_t NRGBs = 6;
const Int_t NCont = 999;

  Double_t stops[NRGBs] = { 0.00, 0.1, 0.34, 0.61, 0.84, 1.00 };
  Double_t red[NRGBs]   = { 0.99, 0.0, 0.00, 0.87, 1.00, 0.51 };
  Double_t green[NRGBs] = { 0.00, 0.0, 0.81, 1.00, 0.20, 0.00 };
  Double_t blue[NRGBs]  = { 0.99, 0.0, 1.00, 0.12, 0.00, 0.00 };


  TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
  gStyle->SetNumberContours(NCont);

  */
}

void event3D()
{
  gSystem->Load("libGeom");
   TGeoManager *geom = new TGeoManager("simple1", "Simple geometry");

   //--- define some materials
   TGeoMaterial *matVacuum = new TGeoMaterial("Vacuum", 0,0,0);
   TGeoMaterial *matAl = new TGeoMaterial("Al", 26.98,13,2.7);
//   //--- define some media
   TGeoMedium *Vacuum = new TGeoMedium("Vacuum",1, matVacuum);
   TGeoMedium *Al = new TGeoMedium("Root Material",2, matAl);

     //--- make the top container volume
   Double_t worldx = 110.;
   Double_t worldy = 50.;
   Double_t worldz = 5.;

   TGeoVolume *top = geom->MakeBox("TOP", Vacuum, 100., 100., 100.);
   geom->SetTopVolume(top);

   gGeoManager->GetVolume("TOP")->InvisibleAll();

   //Measurements (in cm)
   
   Double_t target_rad = 0.3;
   Double_t target_len = 48.58;

   Double_t RTPC_len = 41.;
   
   Double_t foil1_R1 = 2.;
   Double_t foil1_R2 = 2.0000018;

   Double_t foil2_R1 = 3.;
   Double_t foil2_R2 = 3.0000018;

   Double_t GEM1_R1 = 7.;
   Double_t GEM1_R2 = 7.0005;

   Double_t GEM2_R1 = 7.3;
   Double_t GEM2_R2 = 7.3005;

   Double_t GEM3_R1 = 7.6;
   Double_t GEM3_R2 = 7.6005;

   Double_t Z_offset = 0;

   Int_t trans_lvl = 70;
   
   TGeoTube *target = new TGeoTube("target_tube",0, target_rad, target_len/2);
   
   TGeoVolume *target_vol = new TGeoVolume("target_vol",target, Al); //(*)

   top->AddNode(target_vol,1,new TGeoTranslation(0,0, Z_offset));
   gGeoManager->GetVolume("target_vol")->SetTransparency(0);   

   TGeoTube *foil1 = new TGeoTube("foil1_tube", foil1_R1, foil1_R2, RTPC_len/2);
   
   TGeoVolume *foil1_vol = new TGeoVolume("foil1_vol", foil1, Al); //(*)

   top->AddNode(foil1_vol,1,new TGeoTranslation(0,0, Z_offset));


   TGeoTube *foil2 = new TGeoTube("foil2_tube", foil2_R1, foil2_R2, RTPC_len/2);
   
   TGeoVolume *foil2_vol = new TGeoVolume("foil2_vol", foil2, Al); //(*)

   top->AddNode(foil2_vol,1,new TGeoTranslation(0,0, Z_offset));


   TGeoTube *gem1 = new TGeoTube("gem1_tube", GEM1_R1, GEM1_R2, RTPC_len/2);
   TGeoVolume *gem1_vol = new TGeoVolume("gem1_vol", gem1, Al); //(*)
   gem1_vol->SetLineColor(kOrange);
   top->AddNode(gem1_vol,1,new TGeoTranslation(0,0, Z_offset));

   TGeoTube *gem2 = new TGeoTube("gem2_tube", GEM2_R1, GEM2_R2, RTPC_len/2);
   TGeoVolume *gem2_vol = new TGeoVolume("gem2_vol", gem2, Al); //(*)
   gem2_vol->SetLineColor(kGreen);
   top->AddNode(gem2_vol,1,new TGeoTranslation(0,0, Z_offset));

   TGeoTube *gem3 = new TGeoTube("gem3_tube", GEM3_R1, GEM3_R2, RTPC_len/2);
   TGeoVolume *gem3_vol = new TGeoVolume("gem3_vol", gem3, Al); //(*)
   top->AddNode(gem3_vol,1,new TGeoTranslation(0,0, Z_offset));
   gem3_vol->SetLineColor(kRed);
   gGeoManager->GetVolume("gem3_vol")->SetTransparency(trans_lvl);   

   
   //--- draw the ROOT box.
   // by default the picture will appear in the standard ROOT TPad.
   //if you have activated the following line in system.rootrc,
   //it will appear in the GL viewer
   //#Viewer3D.DefaultDrawOption:   ogl
   
   geom->CloseGeometry();
   
   geom->SetVisLevel(4);
  
   //   top->SetVisibility(kFALSE);
   top->Draw("ogl");
   gWorld = top;
}

void event_macro(int ievent=0)
{
  
  set_color_env();
  //    gStyle->SetOptFit(1111);
  
  Bool_t VERBOSE = 0;
  
  gStyle->SetOptStat(0);
  
  Double_t rad2deg = 180/(4*atan(1.0));


  
  Int_t    steps;
  Int_t    rsteps;
  Int_t    steps_temp;//I need this variable 
  Int_t    fIndex;
  Int_t    fHit;

  Int_t    fSenPad[300] ;
  Int_t    fTDC[300] ;
  Int_t    fADC[300] ;
  Double_t fXRec[300], fYRec[300], fZRec[300];

  Double_t fX[300], fY[300], fZ[300];



  Double_t xrec[300],yrec[300],zrec[300];
  Double_t x[300],y[300],z[300];

  
  Double_t fEdep;

  Double_t fkineEne;
  Int_t    fPid;
  Double_t fTheta;
  Double_t fPhi;

  //   TFile *infile=new TFile("nt_out_nowire.root");

  TFile *infile=new TFile(" nt_out_1Ps24Pprime_1000.root");
  TTree *RTPCTree=(TTree*)infile->Get("ep");
  
  Int_t Entries = RTPCTree->GetEntries();
  cout<<"Entries: "<<Entries<<endl;
  
  RTPCTree ->SetBranchAddress("Index", &fIndex);
  RTPCTree ->SetBranchAddress("HitNum", &fHit);
  RTPCTree ->SetBranchAddress("Pid", &fPid);
  RTPCTree ->SetBranchAddress("StepID", &fSenPad[steps]);
  RTPCTree ->SetBranchAddress("StepTDC", &fTDC[steps]);
  RTPCTree ->SetBranchAddress("StepADC", &fADC[steps]);

  RTPCTree ->SetBranchAddress("StepX_rec", &fXRec[steps]);
  RTPCTree ->SetBranchAddress("StepY_rec", &fYRec[steps]);
  RTPCTree ->SetBranchAddress("StepZ_rec", &fZRec[steps]);


  RTPCTree ->SetBranchAddress("StepX", &fX[steps]);
  RTPCTree ->SetBranchAddress("StepY", &fY[steps]);
  RTPCTree ->SetBranchAddress("StepZ", &fZ[steps]);

  

  Int_t ncol = 100;
  Int_t nrow = 200;
  
  
  TH2F *hDisplay = new TH2F("hDisplay", "PAD Display", ncol+11, -5, ncol+5, nrow+11, -5, nrow+5);
  
  Int_t col, row;
    
  //Call the 3D geometry
  if(!gWorld)  event3D();
  else gWorld->Draw("ogl");
  
    hDisplay->SetXTitle("columns (Z direction)"); 
    hDisplay->SetYTitle("rows (Phi evolute)");
    hDisplay->SetZTitle("TDC time"); 

    Int_t p = 25*ievent;

   for(Int_t i = p; i <p+ 25; i++)
    {

      RTPCTree->GetEntry(i);

      steps_temp = fHit;

      //RTPCTree   ->Show(i,200);  
      //	 
      
          cout<<"Hits Num: "<<fHit<<endl;

      steps = 0;
      rsteps = 0;
      
      for (Int_t k=0;k<steps_temp;k++)
	{
	  //	      cout<<"Step: "<<fSenPad[k]<<endl;


	  if (fXRec[k]!=0) 
	    {
	      
	      cout<<"fRec[k]=("<<fXRec[k]<<", "<<fYRec[k]<<", "<<fZRec[k]<<")"<<endl;
	      xrec[steps] = fXRec[k]/10;// so we work in cm
	      yrec[steps] = fYRec[k]/10;
	      zrec[steps] = fZRec[k]/10;
	      steps++;
	    }

	  if (fX[k]!=0) 
	    {
	      x[rsteps] = fX[k]/10;// so we work in cm
	      y[rsteps] = fY[k]/10;
	      z[rsteps] = fZ[k]/10;
	      rsteps++;
	    } 



	  
	}

      for (Int_t k=0;k<steps;k++)
	{
	  row = fSenPad[k]/ncol;
	  col = fSenPad[k]%ncol; 
	  int time = fTDC[k];
	  //	  cout<<"row: "<<row<<" col: "<<col<<endl<<endl;
	  //hDisplay->Fill(col,row, time);
	  hDisplay->Fill(col,row);
	}
      
      //      cout<<"Entry: "<<fIndex<<endl;

      //cout<<steps<<endl;
      TPolyLine3D *track3D = new TPolyLine3D(steps,xrec,yrec,zrec);
      TPolyLine3D *track3Dr = new TPolyLine3D(rsteps,x,y,z);
   
      track3D->SetLineWidth(2);
      track3D->SetLineColor(kBlue);
      track3D->Draw("same");

      track3Dr->SetLineWidth(2);
      track3Dr->SetLineColor(kRed);
      track3Dr->Draw("same");

    }

   // gStyle->SetPalette(55);
    gStyle->SetCanvasPreferGL(kTRUE);

    TCanvas *c1a = new TCanvas("c1a","DISPLAY",1000,10,800,600);
    c1a->cd();
    hDisplay->Draw("COLZTEXT");    
  
    //   c1a->SetGrid();

    //    hDisplay->Draw("lego20");
    //    hDisplay->Draw("COLZ");

    //   DrawGrid();

    
}

void DrawGrid()
{

Int_t ncol = 100;
  Int_t nrow = 200;
  
  
   TPad *grid = new TPad("grid","",0,0,1,1);
   grid->Draw();
   grid->cd();
   grid->SetGrid();
   grid->SetFillStyle(4000);
   grid->SetFrameFillStyle(0);


   TH2 *hgrid = new TH2C("hgrid","", ncol+1, -5, ncol+5, nrow, -5, nrow-1+5);
   hgrid->Draw();
   hgrid->GetXaxis()->SetNdivisions(6,100);
   hgrid->GetYaxis()->SetNdivisions(200);
   hgrid->GetYaxis()->SetLabelOffset(999.);
   hgrid->GetXaxis()->SetLabelOffset(999.);
}

