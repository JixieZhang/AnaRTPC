//This script is to plot trajectory loss of 100ns tic vs 200ns tic
//
/*
TChain *t1=new TChain("ep");
TChain *t2=new TChain("ep");
t1->Add("nt_P150T90_100ns.root/ep");
t2->Add("nt_P150T90_200ns.root/ep");
t1->Draw("HitNum_m>>ht1(80,0,80)","HitNum_m>5","");
t2->Draw("HitNum_m>>ht2","HitNum_m>5","same");
ht2->SetLineColor(2);
*/

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


TF1* FitGaussian(TH1* h1, double &mean, double &sigma, double range_in_sigma=1.0)
{
	//cout<<"h1->GetEntries()="<<h1->GetEntries()<<endl;
	if(h1->GetEntries()<50) return NULL;

	double xmin=h1->GetMean()-1.0*h1->GetRMS();
	double xmax=h1->GetMean()+1.0*h1->GetRMS();
	h1->Fit("gaus","RQ","",xmin,xmax);
	TF1 *f=(TF1*)h1->GetListOfFunctions()->FindObject("gaus");
	mean=f->GetParameter(1);
	sigma=f->GetParameter(2);
	xmin=mean-range_in_sigma*sigma;
	xmax=mean+range_in_sigma*sigma;

	h1->Fit("gaus","RQ","",xmin,xmax);
	f=(TF1*)h1->GetListOfFunctions()->FindObject("gaus");
	mean=f->GetParameter(1);
	sigma=f->GetParameter(2);

	if(gStyle->GetOptFit()==0)
	{
		char str[100];
		TText *text=0;

		double xx=gStyle->GetPadLeftMargin()+0.03; 
		TPaveText *pt = new TPaveText(xx,0.20,xx+0.45,0.45,"brNDC");
		pt->SetBorderSize(0);
		pt->SetFillColor(0);
		sprintf(str,"Mean = %.3G",mean);
		text=pt->AddText(str);
		text->SetTextColor(2);
		sprintf(str,"Sigma = %.3G",sigma);
		text=pt->AddText(str);
		text->SetTextColor(2);
		pt->Draw("same");
	}

	return f;
}

TF1* FitGaussian(TH1* h1,double range_in_sigma=1.0)
{
	double mean,sigma;
	return FitGaussian(h1,mean,sigma,range_in_sigma);
}


void cmp_dPdPhi(double P0_Mev, double Theta0_Deg)
{
  gStyle->SetOptFit(0);
  system("mkdir -p _Gr_kuhn");
  char key[100],file1[255],file2[255];
  sprintf(key,"P%.0fT%.0f",P0_Mev,Theta0_Deg);
  sprintf(file1,"nt_%s_100ns.root/ep",key);
  sprintf(file2,"nt_%s_200ns.root/ep",key);
  TChain *t1=new TChain("ep");t1->Add(file1);
  TChain *t2=new TChain("ep");t2->Add(file2);

  double mean1,sigma1;
  double mean2,sigma2;
  TF1 *f1, *f2;

  double bincenter, binwidth, Prms=0.002;

  TCanvas *c11=new TCanvas("c11","dP",600,500);
  t1->Draw("P0_rec_p-P0_p>>hdP0(100,-0.05,0.05)","HitNum_m>5");
  TH1F* hdP0=(TH1F*)gROOT->FindObject("hdP0");
  Prms=hdP0->GetRMS();
  bincenter=(int(hdP0->GetMean()*1000))/1000.;
  binwidth=Prms/6.;
    
  char strPBin[200];
  sprintf(strPBin,"100,%.4f,%.4f",bincenter-50*binwidth,bincenter+50*binwidth);
  cout<<"strPBin=("<<strPBin<<")\n";

  /////////////////////////////////////////////////
  t1->Draw(Form("P0_rec_p-P0_p>>hdP1(%s)",strPBin),"HitNum_m>5");
  TH1F* hdP1=(TH1F*)gROOT->FindObject("hdP1");
  hdP1->SetTitle("P_{0}-P_{rec}");
  hdP1->SetXTitle("P_{0}-P_{rec} (GeV/c)");
  f1=FitGaussian(hdP1, mean1, sigma1);
  hdP1->SetLineColor(1);
  hdP1->SetMarkerColor(1);
  hdP1->SetMarkerStyle(4);
  f1->SetLineColor(1);
  //gPad->SetLogy(1);
  t2->Draw(Form("P0_rec_p-P0_p>>hdP2(%s)",strPBin),"HitNum_m>5","");
  TH1F* hdP2=(TH1F*)gROOT->FindObject("hdP2");
  hdP2->SetTitle("P_{0}-P_{rec}");
  hdP2->SetXTitle("P_{0}-P_{rec} (GeV/c)");
  f2=FitGaussian(hdP2, mean2, sigma2);
  hdP2->SetLineColor(2);
  hdP2->SetMarkerColor(2);
  hdP2->SetMarkerStyle(25);
  f2->SetLineColor(2);

  c11->Clear();
  if(hdP1->GetMaximum() > hdP2->GetMaximum())
    {
      hdP1->Draw("p");  
      hdP2->Draw("psame");
    }
  else
    {
      hdP2->Draw("p");
      hdP1->Draw("psame");
    }
  
  TLegend *lg=new TLegend(0.65,0.72,1-gPad->GetRightMargin(),1.0-gPad->GetTopMargin(),"","brNDC");
  lg->SetFillColor(0);	lg->SetFillStyle(4000);
  lg->AddEntry("",Form("P=%.0f MeV, #theta=%.0f^{o}",P0_Mev,Theta0_Deg),"");
  lg->AddEntry(hdP1,Form("100ns: #sigma=%.4f",sigma1),"pl");
  lg->AddEntry(hdP2,Form("200ns: #sigma=%.4f",sigma2),"pl");
  lg->Draw();
  c11->SaveAs(Form("_Gr_kuhn/cmp_dP_%s.png",key));


  double mean21,sigma21;
  double mean22,sigma22;
  TF1 *f21, *f22;
  TCanvas *c12=new TCanvas("c12","dPhi",600,500);
  c12->cd();
  t1->Draw("(Phi0_rec_p-Phi0_p)*57.3>>hdPhi0(200,-50,50)","HitNum_m>5");
  TH1F* hdPhi0=(TH1F*)gROOT->FindObject("hdPhi0");
  Prms=hdPhi0->GetRMS();
  bincenter=(int(hdPhi0->GetMean()*1000))/1000.;
  binwidth=Prms/6.;
  
  char strPhiBin[200];
  sprintf(strPhiBin,"100,%.0f,%.0f",bincenter-50*binwidth,
	  bincenter+50*binwidth);
  cout<<"strPhiBin=("<<strPhiBin<<")\n";

  ///////////////////////////////////////////////////////
  t1->Draw(Form("(Phi0_rec_p-Phi0_p)*57.3>>hdPhi1(%s)",strPhiBin),"HitNum_m>5");
  TH1F* hdPhi1=(TH1F*)gROOT->FindObject("hdPhi1");
  hdPhi1->SetTitle("#phi_{0}-#phi_{rec}");
  hdPhi1->SetXTitle("#phi_{0}-#phi_{rec} (deg)");
  f21=FitGaussian(hdPhi1, mean21, sigma21);
  hdPhi1->SetLineColor(1);
  hdPhi1->SetMarkerColor(1);
  hdPhi1->SetMarkerStyle(4);
  f21->SetLineColor(1);
  //gPad->SetLogy(1);
  t2->Draw(Form("(Phi0_rec_p-Phi0_p)*57.3>>hdPhi2(%s)",strPhiBin),"HitNum_m>5","");
  TH1F* hdPhi2=(TH1F*)gROOT->FindObject("hdPhi2");
  hdPhi2->SetTitle("#phi_{0}-#phi_{rec}");
  hdPhi2->SetXTitle("#phi_{0}-#phi_{rec} (deg)");
  f22=FitGaussian(hdPhi2, mean22, sigma22);
  hdPhi2->SetLineColor(2);
  hdPhi2->SetMarkerColor(2);
  hdPhi2->SetMarkerStyle(25);
  f22->SetLineColor(2);

  c12->Clear();
  if(hdPhi1->GetMaximum() > hdPhi2->GetMaximum())
    {
      hdPhi1->Draw("p");  
      hdPhi2->Draw("psame");
    }
  else
    {
      hdPhi2->Draw("p");
      hdPhi1->Draw("psame");
    }
  
  TLegend *lg2=new TLegend(0.65,0.72,1-gPad->GetRightMargin(),1.0-gPad->GetTopMargin(),"","brNDC");
  lg2->SetFillColor(0);	lg->SetFillStyle(4000);
  lg2->AddEntry("",Form("P=%.0f MeV, #theta=%.0f^{o}",P0_Mev,Theta0_Deg),"");
  lg2->AddEntry(hdPhi1,Form("100ns: #sigma=%.1f",sigma21),"pl");
  lg2->AddEntry(hdPhi2,Form("200ns: #sigma=%.1f",sigma22),"pl");
  lg2->Draw();
  c12->SaveAs(Form("_Gr_kuhn/cmp_dPhi_%s.png",key));
}
