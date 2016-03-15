//This script is to plot trajectory loss of 100ns tic vs 200ns tic
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

void trackloss(double P0_Mev, double Theta0_Deg)
{
  system("mkdir -p _Gr_kuhn");
  char key[100],file1[255],file2[255];
  sprintf(key,"P%.0fT%.0f",P0_Mev,Theta0_Deg);
  sprintf(file1,"nt_%s_100ns.root/ep",key);
  sprintf(file2,"nt_%s_200ns.root/ep",key);
  TChain *t1=new TChain("ep");t1->Add(file1);
  TChain *t2=new TChain("ep");t2->Add(file2);

  TCanvas *c11=new TCanvas("c11","track loss",600,500);

  t1->Draw("HitNum_m>>ht1(80,0,80)","HitNum_m>5",""); 
  TH1F* ht1=(TH1F*)gROOT->FindObject("ht1");
  ht1->SetTitle("Number of tracks for TIC=100ns(black) and 200ns(red)");
  ht1->SetXTitle("Number of TDC Hits");
  gPad->SetLogy(1);
  t2->Draw("HitNum_m>>ht2","HitNum_m>5","same");
  TH1F* ht2=(TH1F*)gROOT->FindObject("ht2");
  ht2->SetLineColor(2);
  TLegend *lg=new TLegend(0.65,0.72,1-gPad->GetRightMargin(),1.0-gPad->GetTopMargin(),"","brNDC");
  lg->SetFillColor(0);	lg->SetFillStyle(4000);
  lg->AddEntry("",Form("P=%.0f MeV, #theta=%.0f^{o}",P0_Mev,Theta0_Deg),"");
  lg->AddEntry(ht1,Form("100ns: %.0f",ht1->GetEntries()),"l");
  lg->AddEntry(ht2,Form("200ns: %.0f",ht2->GetEntries()),"l");
  lg->Draw();
  c11->SaveAs(Form("_Gr_kuhn/track_loss_%s.png",key));
}
