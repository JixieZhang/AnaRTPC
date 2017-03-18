//this script is used to fit pol7 to 2-D histo of 
//dPt(or dTh,dPh) vs Pt, Th, Ph
//for RTPC calibration
#include <fstream>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string>
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

/*
TF1* FitGaus(TH1* h1, double &mean, double &sigma, double range_in_sigma=1.0)
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

TF1* FitGaus(TH1* h1,double range_in_sigma=1.0)
{
	double mean,sigma;
	return FitGaus(h1,mean,sigma,range_in_sigma);
}
*/



//input: var="XS", "ALL" or "KLL" 
TF1* FitPol(TH1F *h1, const char *name, const char *title, double min, double max, int order=7)
{   
	char strF[200];
	ofstream fout;
	fout.open("RTPC_Calib_Para.inc",ios_base::app);

	char fcn[10];
	sprintf(fcn,"pol%d",order);
	h1->Fit(fcn,"RQ","",min,max);
	h1->Fit(fcn,"AR","",min,max);
	TF1 *f1 = (TF1*) h1->GetListOfFunctions()->FindObject(fcn);

	sprintf(strF,"\n//Pol%d fitted parameters for %s\n",order,title);
	printf(strF);   
	fout<<strF;
	sprintf(strF,"const double kPara_Pol%d_%s[%d] = {\n\t",order,name,order+1);
	printf(strF);   
	fout<<strF;
	for(int i=0;i<order;i++) 
	{
		sprintf(strF,"%.6E, ",f1->GetParameter(i)); 
		printf(strF);   
		fout<<strF;
		if (!((i+1)%5)) 
		{
			sprintf(strF,"\n\t");
			printf(strF);   
			fout<<strF;
		}
	}
	sprintf(strF,"%.6E \n};\n\n",f1->GetParameter(order));
	printf(strF);   
	fout<<strF;
	fout.close();

	return f1;
} 

//callibrate Theta0
void CalibTheta()
{
	//Fit Theta0_p-Theta_rec:Theta_rec over range of 20<R<160 degrees
	const double deg = atan(1.0)/45.;
	TTree *ep = (TTree*)gDirectory->Get("ep");

	char strName[100], strTitle[255], strVar[255];
	sprintf(strName,"dThVsTh");
	sprintf(strTitle,"Theta0_p-Theta_rec vs Theta_rec");
	sprintf(strVar,"Theta0_p-Theta_rec:Theta_rec");

	
	TCanvas *c10 = new TCanvas("c10","",0,10,800,600);
	ep->Draw(Form("%s >> %s",strVar,strName),"Smax>58","prof");
	TH1F *h1=(TH1F*) gROOT->FindObject(strName);
	h1->SetTitle(strTitle);
	h1->SetXTitle("#theta_{rec} (rad)");
	h1->SetYTitle("d#theta (rad)");

	//FitPol(h1, strName, strTitle, 20*deg, 160*deg, 7);
	FitPol(h1, strName, strTitle, 10*deg, 170*deg, 7);

	c10->cd();
	c10->SaveAs(Form("Graph/%s.png",strName));

}


//callibrate Phi0
void CalibPhi()
{
	//Fit Theta0_p-Theta_rec:Theta_rec over range of 20<R<160 degrees
	TTree *ep = (TTree*)gDirectory->Get("ep");

	char strName[100], strTitle[255], strVar[255];
	sprintf(strName,"dPhVsP");
	sprintf(strTitle,"Phi0_p-Phi_rec vs P0_p");
	sprintf(strVar,"Phi0_p-Phi_rec:P0_p");

	
	TCanvas *c20 = new TCanvas("c10","",0,10,800,600);
	ep->Draw(Form("%s >> %s",strVar,strName),"Smax>58 && abs(Phi0_p-Phi_rec)<0.1","prof");
	TH1F *h1=(TH1F*) gROOT->FindObject(strName);
	h1->SetTitle(strTitle);
	h1->SetXTitle("P (GeV/c)");
	h1->SetYTitle("d#phi(rad)");

	FitPol(h1, strName, strTitle,0.05, 0.3, 7);
	
	c20->cd();
	c20->SaveAs(Form("Graph/%s.png",strName));

}

//callibrate Phi0
void CalibPhi_by_Theta()
{
	//Fit Theta0_p-Theta_rec:Theta_rec over range of 20<R<160 degrees
	const double deg = atan(1.0)/45.;
	TTree *ep = (TTree*)gDirectory->Get("ep");

	char strName[100], strTitle[255], strVar[255];
	sprintf(strName,"dPhVsTh");
	sprintf(strTitle,"Phi0_p-Phi_rec vs Theta_rec");
	sprintf(strVar,"Phi0_p-Phi_rec:Theta_rec");

	
	TCanvas *c20 = new TCanvas("c10","",0,10,800,600);
	ep->Draw(Form("%s >> %s",strVar,strName),"Smax>58 && abs(Phi0_p-Phi_rec)<0.1","prof");
	TH1F *h1=(TH1F*) gROOT->FindObject(strName);
	h1->SetTitle(strTitle);
	h1->SetXTitle("#theta_{rec} (rad)");
	h1->SetYTitle("d#phi(rad)");

	//FitPol(h1, strName, strTitle, 20*deg, 160*deg, 7);
	FitPol(h1, strName, strTitle, 10*deg, 170*deg, 7);
	
	c20->cd();
	c20->SaveAs(Form("Graph/%s.png",strName));

}



//callibrate P0
void CalibP()
{
	TTree *ep = (TTree*)gDirectory->Get("ep");

	TCanvas *c30 = new TCanvas("c10","",600,850);
	c30->Divide(1,3);
	char strName[100], strTitle[255], strVar[255];

	
	//Fit P0*sin(Theta0_p)/abs(MeanBz):R_sim over range of 20<R<50
	c30->cd(1);
	sprintf(strName,"PperpOverBzVsR_1");
	sprintf(strTitle,"P0_p*sin(Theta0_p)/abs(MeanBz) vs R_sim");
	sprintf(strVar,"P0_p*sin(Theta0_p)/abs(MeanBz):R_sim");
	ep->Draw(Form("%s >> %s",strVar,strName),"Smax>40 && R_sim>15 && R_sim<60","prof");
	TH1F *h1=(TH1F*) gROOT->FindObject(strName);
	h1->SetTitle(strTitle);
	h1->SetXTitle("R_sim (mm)");
	h1->SetYTitle("Pperp/Bz");
	FitPol(h1, strName, strTitle, 19, 55, 7);
	
	//Fit P0*sin(Theta0_p)/abs(MeanBz):R_sim over range of 40<R<110
	c30->cd(2);
	sprintf(strName,"PperpOverBzVsR_2");
	sprintf(strTitle,"P0_p*sin(Theta0_p)/abs(MeanBz) vs R_sim");
	sprintf(strVar,"P0_p*sin(Theta0_p)/abs(MeanBz):R_sim");
	ep->Draw(Form("%s >> %s",strVar,strName),"Smax>40 && R_sim>30 && R_sim<120","prof");
	TH1F *h2=(TH1F*) gROOT->FindObject(strName);
	h2->SetTitle(strTitle);
	h2->SetXTitle("R_sim (mm)");
	h2->SetYTitle("Pperp/Bz");
	FitPol(h2, strName, strTitle, 45, 105, 3);


	//Fit P0*sin(Theta0_p)/abs(MeanBz):R_sim over range of 90<R<210
	c30->cd(3);
	sprintf(strName,"PperpOverBzVsR_3");
	sprintf(strTitle,"P0_p*sin(Theta0_p)/abs(MeanBz) vs R_sim");
	sprintf(strVar,"P0_p*sin(Theta0_p)/abs(MeanBz):R_sim");
	ep->Draw(Form("%s >> %s",strVar,strName),"Smax>40 && R_sim>100 && R_sim<300","prof");
	TH1F *h3=(TH1F*) gROOT->FindObject(strName);
	h3->SetTitle(strTitle);
	h3->SetXTitle("R_sim (mm)");
	h3->SetYTitle("Pperp/Bz");
	FitPol(h3, strName, strTitle, 100, 220, 1);

	c30->cd();
	c30->SaveAs("Graph/PperpOverBzVsR.png");
}

void Calib(const char *infile="nt_ep.root")
{
	TFile *InFile = (TFile*)gROOT->GetListOfFiles()->FindObject(infile);
	if (!InFile) 
	{
		InFile = TFile::Open(infile);
		if(!InFile)  {
			cout<<"\nCan not open input file \""<<infile<<"\", I quit ...\n\n";
			exit(-1);
		}
	}

	system("mkdir -p Graph");	
	system("rm -f RTPC_Calib_Para.inc");

	//callibrate Theta0
	CalibTheta();

	//callibrate Phi0
	CalibPhi();
	CalibPhi_by_Theta();

	//callibrate P0
	CalibP();

}
