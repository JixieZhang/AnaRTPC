//This script is to plot the trajectory of RTPC events
//
#include <iostream>
using namespace std;

#include <TQObject.h>
#include "TROOT.h"
#include <TStyle.h>
#include <TF1.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TTree.h>
#include <TFile.h>
#include <TVector3.h>
#include <TLorentzVector.h>
#include <TCanvas.h>
#include <TPad.h>
#include <TProfile.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TVirtualFitter.h>
#include <TMinuit.h>
#include <TMath.h>

#include "TLegend.h"
#include "TString.h"
#include "TCut.h"
#include "TCutG.h"
#include "TPaveText.h"
#include "TText.h"
#include "TSystem.h"
#include "TMultiGraph.h"
#include "TLorentzVector.h"
#include <TChain.h>
#include "TLine.h"

TH2F *h2yxframe = 0;
TH2F *h2xzframe = 0;
TH2F *h2yzframe = 0;
TH2F *h2rzframe = 0;
TCanvas *c1 = 0;
int  thisevent = 0;


void Init()
{
  c1 = new TCanvas("c1","RTPC track display",800,800);
  
  h2yxframe=new TH2F("yxframe","y-x",
		     200,-100,100,200,-100,100);

  h2yxframe->SetXTitle("X (mm)  ");
  h2yxframe->SetYTitle("Y (mm)  ");
  h2yxframe->SetMarkerStyle(1);
  h2yxframe->SetMarkerColor(3);

  for(int i=0;i<3600;i++)
    {
      double phi = i*2*3.14159/3600;
      h2yxframe->Fill(20*cos(phi),20*sin(phi));
      h2yxframe->Fill(30*cos(phi),30*sin(phi));
      h2yxframe->Fill(70*cos(phi),70*sin(phi));
      h2yxframe->Fill(80*cos(phi),80*sin(phi));
    }

  h2yzframe=new TH2F("yzframe","y-z",
		     400,-200,200,200,-100,100);

  h2yzframe->SetXTitle("Z (mm)  ");
  h2yzframe->SetYTitle("Y (mm)  ");
  h2yzframe->SetMarkerStyle(1);
  h2yzframe->SetMarkerColor(3);
  
  
  h2xzframe=new TH2F("xzframe","x-z",
		     400,-200,200,200,-100,100);

  h2xzframe->SetXTitle("Z (mm)  ");
  h2xzframe->SetYTitle("X (mm)  ");
  h2xzframe->SetMarkerStyle(1);
  h2xzframe->SetMarkerColor(3);

  h2rzframe=new TH2F("rzframe","r-z",
		     400,-200,200,100,0,100);

  h2rzframe->SetXTitle("Z (mm)  ");
  h2rzframe->SetYTitle("r (mm)  ");
  h2rzframe->SetMarkerStyle(1);
  h2rzframe->SetMarkerColor(3);

  for(int i=0;i<400;i++)
    {
      double z=-200+i*400/400;
      h2yzframe->Fill(z,-20);
      h2yzframe->Fill(z,-30);
      h2yzframe->Fill(z,-70);
      h2yzframe->Fill(z,-80);
      h2yzframe->Fill(z,20);
      h2yzframe->Fill(z,30);
      h2yzframe->Fill(z,70);
      h2yzframe->Fill(z,80);

      h2xzframe->Fill(z,-20);
      h2xzframe->Fill(z,-30);
      h2xzframe->Fill(z,-70);
      h2xzframe->Fill(z,-80);
      h2xzframe->Fill(z,20);
      h2xzframe->Fill(z,30);
      h2xzframe->Fill(z,70);
      h2xzframe->Fill(z,80);

      h2rzframe->Fill(z,20);
      h2rzframe->Fill(z,30);
      h2rzframe->Fill(z,70);
      h2rzframe->Fill(z,80);
    }

}

double GetVariable(const char *var,TCut cut)
{
  TTree *ep=(TTree*) gROOT->FindObject("ep");

  TH1F *htemp = (TH1F*) gROOT->FindObject("htemp");
  if(htemp) delete htemp;
  ep->Draw(Form("%s>>htemp",var),cut);
  htemp = (TH1F*) gROOT->FindObject("htemp");
  return double(htemp->GetMean());
}

int GetEntries(const char *var,TCut cut)
{
  TTree *ep=(TTree*) gROOT->FindObject("ep");

  TH1F *htemp = (TH1F*) gROOT->FindObject("htemp");
  if(htemp) delete htemp;
  ep->Draw(Form("%s>>htemp",var),cut);
  htemp = (TH1F*) gROOT->FindObject("htemp");
  //cout<<"Number of entries = "<<htemp->GetEntries()<<endl;
  return htemp->GetEntries();
}

//display the trajectory of an event
//if entry>0, display the specified event
//otherwise, display the current event
int  DisplayEvent( int entry=0, TCut extracut="")
{
  if(thisevent<=0) Init();

  TTree *ep=(TTree*) gROOT->FindObject("ep");


  int iEvent = (entry>0) ? entry : thisevent++;

  if(iEvent>= ep->GetEntries()) {
    cout<<"Reach the end of root file...\n";
    return -1;
  }
  
  char strcut[255]; sprintf(strcut,"Index==%d",iEvent);
  char strDCcut[255]; sprintf(strDCcut,"Index==%d && StepS>=30 && StepS<=80 && abs(StepZ)<210",iEvent);

  TCut cut = strcut;
  TCut DCcut = strDCcut;

  cut += extracut;
  DCcut += extracut;

  c1->Divide(2,2);
  c1->cd(1);
  double P0 = GetVariable("P0_p",cut);
  c1->cd(2);
  double Theta0 = GetVariable("Theta0_p*57.3",cut); 
  c1->cd(3);
  double Phi0 = GetVariable("Phi0_p*57.3",cut);
  c1->cd(4);
  int NHits = GetEntries("HitNum",DCcut);

  if(NHits<1) return 0;

  TPaveText *pt = new TPaveText(0.65,0.65,1-gStyle->GetPadRightMargin(),0.89,"brNDC");
  pt->SetFillColor(0);
  if(gROOT->IsBatch()) pt->SetFillStyle(4000);
  pt->SetBorderSize(0);
  //pt->AddText(Form("Event %d", iEvent));
  pt->AddText(Form("P=%.4f",P0));
  pt->AddText(Form("#theta=%.1f^{o}",Theta0));
  pt->AddText(Form("#phi=%.1f^{o}",Phi0));
  pt->AddText(Form("NHits=%d",NHits));
  h2yxframe->SetTitle(Form("Event %d: y-x", iEvent));

  TObject *obj=0;
  if(obj=gROOT->FindObject("gryx")) delete obj;
  if(obj=gROOT->FindObject("gryz")) delete obj;
  if(obj=gROOT->FindObject("grxz")) delete obj;
  if(obj=gROOT->FindObject("grrz")) delete obj;

  c1->Clear();
  //Get the TGraph object, in one canvas there can be only one 
  c1->cd(); h2yxframe->Draw();
   if(!gROOT->IsBatch()) h2yxframe->Draw("same");
  ep->Draw("StepY:StepX",DCcut,"same");
  TGraph *gryx=(TGraph*)(gROOT->FindObject("Graph")->Clone("gryx"));
  gryx->SetName("gryx");
  gryx->SetMarkerColor(4);gryx->SetMarkerStyle(4);

  c1->cd(); h2yzframe->Draw();
  ep->Draw("StepY:StepZ",DCcut,"same");
  TGraph *gryz=(TGraph*)(gROOT->FindObject("Graph")->Clone("gryz"));
  gryz->SetName("gryz");
  gryz->SetMarkerColor(4);gryz->SetMarkerStyle(4);

  c1->cd(); h2xzframe->Draw();
  ep->Draw("StepX:StepZ",DCcut,"same");
  TGraph *grxz=(TGraph*)(gROOT->FindObject("Graph")->Clone("grxz"));
  grxz->SetName("grxz");
  grxz->SetMarkerColor(4);grxz->SetMarkerStyle(4);
		
  c1->cd(); h2rzframe->Draw();
  ep->Draw("StepS:StepZ",DCcut,"same");
  TGraph *grrz=(TGraph*)(gROOT->FindObject("Graph")->Clone("grrz"));
  grrz->SetName("grrz");
  grrz->SetMarkerColor(4);grrz->SetMarkerStyle(4);



  c1->Clear();
  c1->Divide(2,2);

  c1->cd(1); h2yxframe->Draw();
  pt->Draw();
  if(!gROOT->IsBatch()) h2yxframe->Draw("same");
  gryx->Draw("same p");
  ep->Draw("StepY_rec:StepX_rec>>h2yx_rec",DCcut,"same*");

  c1->cd(2); h2yzframe->Draw();
  gryz->Draw("same p");
  ep->Draw("StepY_rec:StepZ_rec>>h2yz_rec",DCcut,"same*");

  c1->cd(3); h2xzframe->Draw();
  grxz->Draw("same p");
  ep->Draw("StepX_rec:StepZ_rec>>h2xz_rec",DCcut,"same*");

		
  c1->cd(4); h2rzframe->Draw();
  grrz->Draw("same p");
  ep->Draw("StepS_rec:StepZ_rec>>h2rz_rec",DCcut,"same*");

  c1->Update();
  c1->SaveAs(Form("Movie/RTPC_Event_%03d.png",iEvent));

  ep->Scan("P0_p:Theta0_p*57.3:Phi0_p*57.3:Z0:HitNum:Smax:StepS[HitNum-1]:StepPhi[HitNum-1]:StepZ[HitNum-1]",cut);

  return NHits;
}

//slice show n events, each last for specified second
//if istart<=0, start from current events, otherwise from specified event
void go(int n=1,int second=2, int istart=0)
{ 
  if(istart>0)  thisevent=istart;
  for(int i=0;i<n;i++)
    {
      int nHits = 0;
      while (nHits==0) 
	{
	  nHits = DisplayEvent(thisevent++);
	  if(nHits == -1) break;
	}
      c1->Update();
      if(nHits == -1) break;
      gSystem->Sleep(second*1000);
    }

}

void DisplayMyEvent(TCut cut="abs(sqrt(A_sim*A_sim+B_sim*B_sim)-R_sim)>1.0")
{ 
  int nHits = 0;
  while (nHits==0) 
    {
      nHits = DisplayEvent(0,cut);
      if(nHits == -1) break;
    }
  c1->Update();
}

//slice show n events, each last for specified second
//if istart<=0, start from current events, otherwise from specified event
void SliceShowMyEvent(int n=1, int second=2, int istart=0, 
		      TCut cut="abs(sqrt(A_sim*A_sim+B_sim*B_sim)-R_sim)<5.0")
{ 
  if(istart>0)  thisevent=istart;
  for(int i=0;i<n;i++)
    {
      int nHits = 0;
      while (nHits==0) 
	{
	  nHits = DisplayEvent(thisevent++,cut);
	  if(nHits == -1) break;
	}
      c1->Update();
      if(nHits == -1) break;
      gSystem->Sleep(second*1000);

    }
}
