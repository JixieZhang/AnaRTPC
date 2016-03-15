
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TSystem.h>
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

void drawline(int istart=0, int iend=2)
{
  system("mkdir -p Graph/Gif");
  TTree *ep = (TTree*) gROOT->FindObject("ep");
  ep->SetAlias("CellID","StepID");
  ep->SetAlias("TDC","int(StepTDC+(Index%25)*2)");
  ep->SetAlias("TDCFold","int((StepTDC+(Index%25)*2)%80)");
  ep->SetAlias("Row","int(StepID/100)");
  ep->SetAlias("Col","int(StepID%100)");
  ep->SetAlias("EventID","int(Index/25)");

  TCanvas *c21=new TCanvas("c21","",900,500);
  TH2F* hFrame=new TH2F("hIDTDCFrame","",120,0,120,20000,0,20000);
  hFrame->SetXTitle("TDC (100ns per tic)");
  hFrame->SetYTitle("Cell ID");

  TH2F* hIDTDC=0;

  for(int ie=istart;ie<iend;ie++)
    {
      c21->cd(); c21->Clear(); 
      hFrame->SetTitle(Form("Event %d",ie));

      hFrame->Draw();
      for(int i=0;i<25;i++)
	{      
	  ep->Draw("CellID:TDC",Form("StepID>=0 && EventID==%d && int(Index%25)==%d",ie,i),"l*same");
	  //ep->Draw("StepID:StepTDC+(Index%25)*2",Form("StepID>=0 && (Index/25)==%d && int(Index%25)==%d",ie,i),"colrtext");
	  
	  gPad->Update();
	  c21->SaveAs(Form("Graph/Gif/SuperCell_Event%03d_%02d.png",ie,i));
	  
	  if(!gROOT->IsBatch()) {cout<<"sleep...\r";gSystem->Sleep(1000);}
	}
      c21->SaveAs(Form("Graph/SuperCell_Event%03d.png",ie));
      system(Form("convert -delay 50 -loop 0 Graph/Gif/SuperCell_Event%03d_??.png Graph/SuperCell_Event%03d.gif",ie,ie));
    }
}
void draw(int istart=1, int iend=20)
{
  TTree *ep = (TTree*) gROOT->FindObject("ep");

  ep->SetAlias("CellID","StepID");
  ep->SetAlias("TDC","int(StepTDC+(Index%25)*2)");
  ep->SetAlias("TDCFold","int((StepTDC+(Index%25)*2)%80)");
  ep->SetAlias("Row","int(StepID/100)");
  ep->SetAlias("Col","int(StepID%100)");
  ep->SetAlias("EventID","int(Index/25)");

  TCanvas *c11=new TCanvas("c11","",900,800);
  c11->SetRightMargin(0.10);
  c11->cd();

  ep->Draw("Row:Col","StepID>=0","*");
  c11->SaveAs("Graph/SuperCell.png");

  ep->Draw("CellID:TDC>>hIDTDC1(120,0,120,10,10,20)","StepID>=0 && EventID<100","colztext");
  c11->SaveAs("Graph/SuperCell_occupancy_10.png");

  TH1F *h1=0;
  TH2F *hIDTDC=0, *hIDTDCFold=0;
  TCanvas *c21=new TCanvas("c21","",900,500);
  c21->SetLeftMargin(0.06);
  c21->SetRightMargin(0.06);
  for(int i=istart;i<iend;i++)
    {
      hIDTDC = (TH2F*) gROOT->FindObject("hIDTDC");
      if(hIDTDC) delete hIDTDC;
      c21->cd();
      ep->Draw("CellID:TDC>>hIDTDC(120,0,120,20000,0,20000)",Form("StepID>=0 && EventID==%d",i),"text");
      //ep->Draw("StepID:StepTDC+(Index%25)*2>>hIDTDC(120,0,120,550,0,550)",Form("StepID>=0 && EventID==%d",i),"text");
      hIDTDC = (TH2F*) gROOT->FindObject("hIDTDC");
      if(hIDTDC) 
	{
	  hIDTDC->SetTitle(Form("Event %d; TDC/100ns; Cell ID",i));
	  hIDTDC->Draw("colrtext");
	}
      gPad->Update();
      c21->SaveAs(Form("Graph/SuperCell_occupancy_%03d.png",i));
      if(!gROOT->IsBatch()) {cout<<"sleep...\r";gSystem->Sleep(1000);}
    }
  cout<<""<<endl;
  
  TCanvas *c22=new TCanvas("c22","",900,500);
  c22->SetLeftMargin(0.06);
  c22->SetRightMargin(0.06);
  for(int i=istart;i<iend;i++)
    {
      hIDTDCFold = (TH2F*) gROOT->FindObject("hIDTDCFold");
      if(hIDTDCFold) delete hIDTDCFold;
      c22->cd();
      ep->Draw("CellID:TDCFold>>+hIDTDCFold(80,0,80,20000,0,20000)",Form("StepID>=0 && EventID==%d",i),"text");
      if(i>0) ep->Draw("CellID:TDCFold>>hIDTDCFold+",Form("StepID>=0 &&TDC>80 && EventID==%d",i-1),"text");
      hIDTDCFold = (TH2F*) gROOT->FindObject("hIDTDCFold");
      if(hIDTDCFold) 
	{
      hIDTDCFold->SetTitle(Form("Event %d; TDC/100ns; Cell ID",i));
      hIDTDCFold->Draw("colrtext");
      }
      gPad->Update();
      ep->Scan("(Index%100):HitNum:P0_p:Theta0_p*57.3",Form("HitNum>5 && EventID==%d",i));
      c22->SaveAs(Form("Graph/SuperCell_fold_occupancy_%03d.png",i));
      if(!gROOT->IsBatch()) {cout<<"sleep...\r";gSystem->Sleep(1000);}
    }
  cout<<""<<endl;


}
