//this script is used to plot 1-D DeltaX and 2-D DeltaX vs Var
//DeltaX = dP,dTheta,dPhi,dZ
//Var = P, Theta, Phi, Z
//This script is good for sanity check of optics|reconstruction

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


void DeltaXXX(const char* infile="nt_ep.root",const char* vipcut="(P0_p>0.07 && P0_p<0.1 && cos(Theta0_p)<-0.2)",
			  const char* cutname="Pr70to100", const char* treename="ep")
{
	
	///////////////////////////////////////////////////////////////
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
	TTree *ep = (TTree*)gDirectory->Get(treename);
	///////////////////////////////////////////////////////////////

	const double deg = atan(1.0)/45.;
	char key[200];              
	sprintf(key,"RTPC12_%s",cutname);


	///////////////////////////////////////////////////////////////
	//I have to put a wide range here otherwise I can not extract the sigma from the distribution
	TCut NoCut = "HitNum>5 && abs(P0_rec_p-P0_p)/P0_p<1.0 && abs(Theta0_rec_p-Theta0_p)<1.0 && abs(Phi0_rec_p-Phi0_p)<1.0 && abs(Z0_rec_p-Z0)<100";
	NoCut += vipcut;
	cout<<"TCut NoCut = \""<<NoCut<<"\""<<endl;

	gStyle->SetStatX(1.0-gStyle->GetPadRightMargin());  
	gStyle->SetStatH(0.060);
	gStyle->SetStatW(0.200);

	//////////////////////////////////////////////////////////
	//config the variables
	//plot 2-D dX, 4 pads:  
	//1) dX vs P0  
	//2) dX vs Theta0  
	//1) dX vs Phi0  
	//1) dX vs Z0
	const int nDelta=4;
	const int nVar=4;
	const char *strYVar[]={"(P0_rec_p-P0_p)*1000","(Theta0_rec_p-Theta0_p)*1000",
		"(Phi0_rec_p-Phi0_p)*1000","Z0_rec_p-Z0"};
	const char *strYTitle[]={"#deltaP (MeV/c)","#delta#theta (mrad)","#delta#phi (mrad)",
		"#deltaZ (mm)"};
	const char *strYName[]={"dP","dTheta","dPhi","dZ"};
	int    YBinNum[]={40,50,50,40};
	double YBinMin[]={-20,-50.,-50.,-10.};
	double YBinMax[]={ 20, 50., 50., 10.};

	double YSigmaFactor[nDelta]={5.0,5.0,5.0,5.0};   //sigma factor to cut on Y when plot 2-D 
	double YMean[nDelta],YSigma[nDelta];
	char   DeltaCut[512],strYCut[nDelta][255]; 
	DeltaCut[0]='\0'; //I have to give this terminal sign to this buffer, otherwise strlen() will not work
	sprintf(DeltaCut,"%s",vipcut);

	const char *strXVar[]={"P0_p*1000","Theta0_p","Phi0_p","Z0"};
	const char *strXTitle[]={"Thrown P (MeV/c)","Thrown #theta (rad)","Thrown #phi (rad)",
		"Thrown Z (mm)"};
	const char *strXName[]={"P","Theta","Phi","Z"};
	int    XBinNum[]={50,60,36,40};
	double XBinMin[]={ 50.,  0*deg,-180*deg,-200.};
	double XBinMax[]={300.,180*deg, 180*deg, 200.};

	char hTitle[255], strTg[255], hName[100], hName_1[100], hName_2[100];

	/////////////////////////////////////////////////////////
	TH1F *h1=0;
	TH2F *h2=0;
	TH1F *hMean=0, *hSigma=0,*hSigma_neg=0;
	

	//Plot 1-D 
	TCanvas *c31 = new TCanvas("c31","",800,700);
	c31->Clear();
	int nCol=int(ceil(nDelta/3.0));
	int nRow=int(ceil(double(nDelta)/double(nCol)));
	c31->Divide(nCol,nRow);
	for(int dd=0;dd<nDelta;dd++)
	{	
		c31->cd(dd+1);//gPad->SetRightMargin(0.12);
		sprintf(hName,"%s",strYName[dd]);
		sprintf(hTitle,"%s ;%s",strYName[dd],strYTitle[dd]);
		sprintf(strTg,"%s",strYVar[dd]);

		//cout<<"Name = "<<hName<<"    var = "<<strTg<<endl;
		//cout<<"Title = "<<hTitle<<endl;

		int UseGivenRange=0;
		if(UseGivenRange)
		{
			//if you know the range
			h1 = new TH1F(hName,hTitle,YBinNum[dd],YBinMin[dd],YBinMax[dd]);
			ep->Project(hName,strTg,NoCut);                                                                                  
			h1->Draw();
		}
		else
		{
			//if you do not know the range
			h1 = (TH1F*) gROOT->FindObject(hName);
			//iteration one, get initial mean and rms then plot again
			
			if(h1) delete h1;
			sprintf(strTg,"%s >> %s",strYVar[dd],hName);
			ep->Draw(strTg,NoCut);
			h1 = (TH1F*) gROOT->FindObject(hName);
			double Mean = h1->GetMean();
			double RMS = h1->GetRMS();
			delete h1;
			sprintf(strTg,"%s >> %s(120,%.5f,%.5f)",strYVar[dd],hName,Mean-3*RMS,Mean+3*RMS);
			ep->Draw(strTg,NoCut);
			h1 = (TH1F*) gROOT->FindObject(hName);
			double Sigma;
			FitGaus(h1,Mean,Sigma);
			double pCut=(YSigmaFactor[dd]+1.0)*Sigma;
			
			//cout<<hName<<":  Mean="<<Mean<<"  RMS="<<h1->GetRMS()<<"  Sigma="<<Sigma<<endl;
			//iteration two, plot again using sigmato get better range 
			delete h1;
			sprintf(strTg,"%s >> %s(100,%.5f,%.5f)",strYVar[dd],hName,Mean-pCut,Mean+pCut);
			ep->Draw(strTg,NoCut);
			h1 = (TH1F*) gROOT->FindObject(hName);

			h1->SetTitle(hTitle); 
		}
		
		h1->SetLineColor(4); 
		FitGaus(h1,YMean[dd],YSigma[dd]); 
		double pCut=YSigmaFactor[dd]*YSigma[dd];
		double pYmax=h1->GetMaximum()*0.7; 
		TLine *L1 = new TLine(-pCut,0,-pCut,pYmax);
		TLine *L2 = new TLine( pCut,0, pCut,pYmax);
		L1->SetLineWidth(2); L1->Draw();
		L2->SetLineWidth(2); L2->Draw();

		sprintf(strYCut[dd],"abs((%s)-(%.3f))<%.3f ",strYVar[dd],YMean[dd],pCut);
		if(strlen(DeltaCut)<1) strcpy(DeltaCut,strYCut[dd]);
		else 
		{
			strcat(DeltaCut," && ");
			strcat(DeltaCut,strYCut[dd]);
		}
	}
	c31->cd(0);
	c31->Modified();
	c31->SaveAs(Form("Graph/Delta_1D_%s.png",key));

	cout<<"\n/////////////////////////////////////////"<<endl;
	cout<<"TCut DeltaCut = \""<<DeltaCut<<"\""<<endl;
	cout<<"/////////////////////////////////////////\n"<<endl;


	/////////////////////////////////////////////////////////
	//plot 2-D
	for(int dd=0;dd<nDelta;dd++)
	{
		char tmpName[100];
		sprintf(tmpName,"c32_%s",strYName[dd]);
		TCanvas *c32 = new TCanvas(tmpName,strYName[dd],(1+dd)*30,(1+dd)*30,900,800);
		int pCol=int(ceil(nVar/3.0));
		int pRow=int(ceil(double(nVar)/double(pCol)));
		c32->Divide(pCol,pRow);
		for(int vv=0;vv<nVar;vv++)
		{
			c32->cd(vv+1);gPad->SetRightMargin(0.12);
			sprintf(hName,"%sVS%s",strYName[dd],strXName[vv]);
			sprintf(hTitle,"%s VS %s;%s ;%s ",strYName[dd],strXName[vv],strXTitle[vv],strYTitle[dd]);
			sprintf(strTg,"%s:%s",strYVar[dd],strXVar[vv]);
			
			//cout<<"Name = "<<hName<<"    var = "<<strTg<<endl;
			//cout<<"Title = "<<hTitle<<endl;

			h2 = new TH2F(hName,hTitle,
				XBinNum[vv],XBinMin[vv],XBinMax[vv],YBinNum[dd],YBinMin[dd],YBinMax[dd]);
			ep->Project(hName,strTg,DeltaCut);
			h2->Draw("colz");
			h2->FitSlicesY();
			sprintf(hName_1,"%s_1",hName);
			sprintf(hName_2,"%s_2",hName);
			hMean  = (TH1F*) gROOT->FindObject(hName_1);
			hSigma = (TH1F*) gROOT->FindObject(hName_2);

			hMean->SetMarkerStyle(20); hMean->SetMarkerColor(1);
			hSigma->SetMarkerStyle(22); hSigma->SetMarkerColor(2);hSigma->SetLineColor(2);
			
			hSigma_neg = (TH1F*) (hSigma->Clone("hSigma_neg"));
			hSigma_neg->Scale(-1.0);

			hMean->Draw("same");
			hSigma->Draw("lcsame");
			hSigma_neg->Draw("lcsame");
		}
		c32->cd(0);
		c32->Modified();
		c32->SaveAs(Form("Graph/%s_%s.png",strYName[dd],key));
	}
	return;
}
