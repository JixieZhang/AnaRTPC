//this script is used to plot resolutions as a function of
//other variables
#include "stdlib.h"
#include <iostream>
#include "math.h"
using namespace std;

#include "TLegend.h"
#include <TH1.h>
#include <TH2.h>
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
#include <TCanvas.h>
#include <TChain.h>


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


//study the pad resolutions
void PlotRes_2(const char *nt1="nt_ep_pad4.5x5.0.root",
			 const char *nt2="nt_ep_pad2.0x2.0.root")
{
	gStyle->SetOptFit(0);
	TChain *t1=new TChain("ep");
	TChain *t2=new TChain("ep");

	t1->Add(nt1);
	t2->Add(nt2);

	TCanvas *c4 = new TCanvas("c4","",800,600);
	t1->Draw("Theta0_p-Theta_sim>>h11(100,-0.1,0.1)");
	t1->Draw("Theta0_p-Theta_rec>>h12(100,-0.1,0.1)");
	t2->Draw("Theta0_p-Theta_sim>>h21(100,-0.1,0.1)");
	t2->Draw("Theta0_p-Theta_rec>>h22(100,-0.1,0.1)");

	TH1F *h11=(TH1F*)gROOT->FindObject("h11");
	TH1F *h12=(TH1F*)gROOT->FindObject("h12");
	TH1F *h21=(TH1F*)gROOT->FindObject("h21");
	TH1F *h22=(TH1F*)gROOT->FindObject("h22"); 

	c4->Clear();
	h11->SetTitle(";d#theta (mrad) ");
	h11->SetLineColor(1);h11->SetMarkerColor(1);
	h12->SetLineColor(2);h12->SetMarkerColor(2);
	h21->SetLineColor(3);h21->SetMarkerColor(3);
	h22->SetLineColor(4);h22->SetMarkerColor(4);

	double Mean,Sigma11,Sigma12,Sigma21,Sigma22;
	FitGaussian(h11, Mean, Sigma11);
	FitGaussian(h12, Mean, Sigma12);
	FitGaussian(h21, Mean, Sigma21);
	FitGaussian(h22, Mean, Sigma22);


	h11->Draw(); 
	h12->Draw("same");	
	h22->Draw("same");

	TLegend *lg=new TLegend(0.6,0.68,1-gPad->GetRightMargin(),1.0-gPad->GetTopMargin(),"","brNDC");
	lg->SetFillColor(0);	lg->SetFillStyle(4000);
	lg->AddEntry(h11,Form("Sim. Track: #sigma_{d#theta}=%.3f",Sigma11),"pl");
	lg->AddEntry(h12,Form("Pad: 4.5x5.0: #sigma_{d#theta}=%.3f",Sigma12),"pl");
	lg->AddEntry(h22,Form("Pad: 2.0x2.0: #sigma_{d#theta}=%.3f",Sigma22),"pl");
	lg->Draw();

	c4->SaveAs("Graph/dTheta_pad_resolution.png");
}


//study the pad resolutions
void PlotRes(const char *nt1="nt_ep_pad4.5x5.0.root",
			 const char *nt2="nt_ep_pad2.0x2.0.root",
			 const char *nt3="nt_ep_p250z0_2dreadout.root")
{
	gStyle->SetOptFit(0);
	TChain *t1=new TChain("ep");
	TChain *t2=new TChain("ep");
	TChain *t3=new TChain("ep");

	t1->Add(nt1);
	t2->Add(nt2);
	t3->Add(nt3);

	TCanvas *c4 = new TCanvas("c4","",800,600);
	t1->Draw("Theta0_p-Theta_sim>>h11(100,-0.1,0.1)");
	t1->Draw("Theta0_p-Theta_rec>>h12(100,-0.1,0.1)");
	t2->Draw("Theta0_p-Theta_sim>>h21(100,-0.1,0.1)");
	t2->Draw("Theta0_p-Theta_rec>>h22(100,-0.1,0.1)");
	t3->Draw("Theta0_p-Theta_sim>>h31(100,-0.1,0.1)");
	t3->Draw("Theta0_p-Theta_rec>>h32(100,-0.1,0.1)");

	TH1F *h11=(TH1F*)gROOT->FindObject("h11");
	TH1F *h12=(TH1F*)gROOT->FindObject("h12");
	TH1F *h21=(TH1F*)gROOT->FindObject("h21");
	TH1F *h22=(TH1F*)gROOT->FindObject("h22"); 
	TH1F *h31=(TH1F*)gROOT->FindObject("h31");
	TH1F *h32=(TH1F*)gROOT->FindObject("h32"); 

	c4->Clear();
	h11->SetTitle(";d#theta (rad) ");
	h11->SetLineColor(1);h11->SetMarkerColor(1);
	h12->SetLineColor(2);h12->SetMarkerColor(2);
	h21->SetLineColor(3);h21->SetMarkerColor(3);
	h22->SetLineColor(4);h22->SetMarkerColor(4);
	h31->SetLineColor(5);h31->SetMarkerColor(5);
	h32->SetLineColor(6);h32->SetMarkerColor(6);

	double Mean,Sigma11,Sigma12,Sigma21,Sigma22,Sigma31,Sigma32;
	FitGaussian(h11, Mean, Sigma11);
	FitGaussian(h12, Mean, Sigma12);
	FitGaussian(h21, Mean, Sigma21);
	FitGaussian(h22, Mean, Sigma22);
	FitGaussian(h31, Mean, Sigma31);
	FitGaussian(h32, Mean, Sigma32);


	h11->Draw(); 
	h12->Draw("same");	
	h22->Draw("same");
	h32->Draw("same");

	TLegend *lg=new TLegend(0.6,0.60,1-gPad->GetRightMargin(),1.0-gPad->GetTopMargin(),"","brNDC");
	lg->SetFillColor(0);	lg->SetFillStyle(4000);
	lg->AddEntry(h11,Form("Sim. Track:   #sigma_{d#theta}=%.3f",Sigma11),"pl");
	lg->AddEntry(h32,Form("2-D Readout:  #sigma_{d#theta}=%.3f",Sigma32),"pl");
	lg->AddEntry(h22,Form("Pad: 2.0x2.0: #sigma_{d#theta}=%.3f",Sigma22),"pl");
	lg->AddEntry(h12,Form("Pad: 4.5x5.0: #sigma_{d#theta}=%.3f",Sigma12),"pl");
	lg->Draw();

	c4->SaveAs("Graph/dTheta_pad_resolution.png");
}



//study the pid of rtpc
//plot P_rec vs P0, dEdX vs P0
void PlotPid(const char *infile="nt_ep.root")
{
	gStyle->SetOptFit(0);

	TFile *InFile = (TFile*)gROOT->GetListOfFiles()->FindObject(infile);
	if (!InFile) 
	{
		InFile = TFile::Open(infile);
		if(!InFile)  {
			cout<<"\nCan not open input file \""<<infile<<"\", I quit ...\n\n";
			exit(-1);
		}
	}
	InFile->cd();

	
	TTree *ep = (TTree*)gDirectory->Get("ep");

	TCanvas *c40 = new TCanvas("c40","",800,600);


	//ep->Draw("P0_rec_p:P0_p>>h(350,0,0.350,350,0,0.35)");
	//ep->Draw("dEdX:P0_p>>h(350,0,0.350,1000,0.0,100.0)");

	TH2F *PrecVsP0 = new TH2F("PrecVsP0","; P0 (GeV/c); P_rec (KeV/c)",
		350,0.0,0.35,350,0.0,0.35);
	TH2F *dEdXVsPrec = new TH2F("dEdXVsPrec","; Prec (GeV/c); dEdX (KeV/mm)",
		350,0.0,0.35,1000,0.0,100.0);
	TH2F *dEdXVsP0 = new TH2F("dEdXVsP0","; P0 (GeV/c); dEdX (KeV/mm)",
		350,0.0,0.35,1000,0.0,100.0);

	TH2F *frame[] = {PrecVsP0,dEdXVsPrec,dEdXVsP0};
	const char *strVar[] = {"P0_rec_p:P0_p","dEdX:P0_rec_p","dEdX:P0_p"};
	const char *strVarHe3[] = {"P0_rec_p*2:P0_p","dEdX:P0_rec_p*2","dEdX:P0_p"};
	const char *strBin[] = { "350,0.0,0.35,350,0.0,0.35",
		"350,0.0,0.35,1000,0.0,100.0", "350,0.0,0.35,1000,0.0,100.0"};

	//char strTg[255];
	TCut CutPi = "Smax>70 && Pid==211";
	TCut CutK  = "Smax>70 && Pid==321";
	TCut CutPr = "Smax>70 && Pid==2212";
	TCut CutHe3= "Smax>70 && Pid>1.0E6";


	for(int ii=0;ii<3;ii++)
	{
		c40->Clear();
		c40->cd();

		TH2F *hPi=(TH2F*)gROOT->FindObject("hPi");if(hPi) delete hPi;
		TH2F *hK=(TH2F*)gROOT->FindObject("hK"); if(hK) delete hK;
		TH2F *hPr=(TH2F*)gROOT->FindObject("hPr"); if(hPr) delete hPr;
		TH2F *hHe3=(TH2F*)gROOT->FindObject("hHe3"); if(hHe3) delete hHe3;

		frame[ii]->Draw();
		ep->Draw(Form("%s >> hPi(%s)",strVar[ii],strBin[ii]),CutPi,"");
		ep->Draw(Form("%s >> hK(%s)", strVar[ii],strBin[ii]),CutK,"");
		ep->Draw(Form("%s >> hPr(%s)",strVar[ii],strBin[ii]),CutPr,"");
		ep->Draw(Form("%s >> hHe3(%s)",strVarHe3[ii],strBin[ii]),CutHe3,"");

		hPi=(TH2F*)gROOT->FindObject("hPi");
		hK=(TH2F*)gROOT->FindObject("hK");
		hPr=(TH2F*)gROOT->FindObject("hPr");
		hHe3=(TH2F*)gROOT->FindObject("hHe3"); 

		hPi->SetLineColor(2);hPi->SetMarkerColor(2);hPi->SetMarkerStyle(2);
		hK->SetLineColor(3);hK->SetMarkerColor(3);hK->SetMarkerStyle(2);
		hPr->SetLineColor(4);hPr->SetMarkerColor(4);hPr->SetMarkerStyle(2);
		hHe3->SetLineColor(7);hHe3->SetMarkerColor(7);hHe3->SetMarkerStyle(2);
	
	
		if(ii>=1) c40->SetLogy(1);
		else c40->SetLogy(0);
		
		frame[ii]->Draw();
		
		hPi->Draw("same");	
		hK->Draw("same");
		hPr->Draw("same");
		hHe3->Draw("same");
		

		TLegend *lg=new TLegend(0.15,0.65,0.35,0.898,"","brNDC");
		lg->SetFillColor(0);	
		lg->SetBorderSize(0);
		lg->SetFillStyle(4000);
		//lg->AddEntry(hPi,"#color[2]{#pi^{+}}","plf");
		//lg->AddEntry(hK,"#color[3]{k^{+}}","plf");
		//lg->AddEntry(hPr,"#color[4]{proton}","plf");
		//lg->AddEntry(hHe3,"#color[7]{^{3}He}","plf");
		lg->AddEntry(hPi,"#color[1]{#pi^{+}}","plf");
		lg->AddEntry(hK,"#color[1]{k^{+}}","plf");
		lg->AddEntry(hPr,"#color[1]{proton}","plf");
		lg->AddEntry(hHe3,"#color[1]{^{3}He}","plf");


		lg->Draw();

		c40->SaveAs(Form("Graph/%s.png",frame[ii]->GetName()));
	}
}
