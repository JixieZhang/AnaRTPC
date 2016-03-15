//////////////////////////////////////////////////////////////////////// 
//MulDFit.cc: implement of MulDFit class.
// Created by jixie zhang ,9/10/2014
//This class is trying to fit multiple dimention polynomial function
//  Y = Sum_abcde {K_abcde * A^a B^b C^c D^d E^e}
//    = Sum_bcde {Sum_a{K_abcde * A^a} * B^b C^c D^d E^e}
//in which self-variables A, B, C, D and E can go up to unlimited order
//
//One need to provide two input files: in_para.ini and in_file.dat
//The first file is the initial parameter file, which defines the fitting function 
//and the initial values for each matrix element. It also tells if a parameter is
//fixed in the fitting or not. Write1stParaFile() can create a template of in_para.ini
//The 2nd file privides the raw data. 
//
// 
//Format of the raw data file "in_file.dat": 
//The raw data file contains 7 collums
// var  self_var1  self_var2  self_var3  self_var4  self_var5 weight
// Line 1: name
// Line 2: Title
// From line 3 are data points
// Ror example:
// P_tg    X_fp    Theta_fp   Y_fp       Phi_fp     P_fp   Weight
//
//Format of the parameter file in_para.ini.  
//line 1-8: self-description
//line 9: MaxOrderA MaxOrderB MaxOrderC MaxOrderD MaxOrderE
//line 10: #lable of collum, which looks like the next line: 
//5     5     5     5     5
//Label ^b   ^c   ^d   ^e       K_0bcde       K_1bcde       K_2bcde       K_3bcde       K_4bcde       K_5bcde fix0 fix1 fix2 fix3 fix4 fix5
//B      0    0    0    0   7.32595e-06      0.996222    0.00165692     -0.011795    -0.0016533     0.0146345    0    0    0    0    0    0
//
//////////////////////////////////////////////////////////////////////
#include "MulDFit.h"

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

#include "TLegend.h"
#include "TString.h"
#include "TCut.h"
#include "TCutG.h"
#include "TPaveText.h"
#include "TText.h"
#include "TPad.h"
#include "TMultiGraph.h"
#include "TLorentzVector.h"
#include <TChain.h>
#include "TLine.h"


//////////////////////////////////////////////////
//global buffer begin
int gMaxOrderA=5;
int gMaxOrderB=5;
int gMaxOrderC=5;
int gMaxOrderD=5;
int gMaxOrderE=5;
int gMaxOrder=5;

//double mPara[kMaxParaA][kMaxParaB][kMaxParaC][kMaxParaD][kMaxParaE];
double ***** gPara=0;
double * gFitPara=0;
double * gFitParaErr=0;
int    * gFitParaFixedFlag=0;
//An index map that pointing 1-D parameter to 5-D array elements 
//mFitParaIndexMap[mParaNum][6], from 0 to 6 it is indexA...indexE and label
int  **  gFitParaIndexMap=0; 
int      gParaNum=0;  //number of parameter

//vector to store the raw data
vector <double> gVVar;  
vector <double> gVSelfVar[5]; 
vector <double> gVWeight;
//static buffer end
//////////////////////////////////////////////////

//I am not use this global function any longer
//user might need this routine
double myMulDPolN(double *x, double *par)
{
	//multi-dimention polN
	//  V = Sum_abcde {C_abcde * X^a T^b Y^C F^d P^e}
	//    = Sum_bcde {Sum_a{C_abcde * X^a} * T^b Y^C F^d P^e}

	//now need to set parameters values to the array element
	for(int i=0;i<gParaNum;i++)
	{
		gFitPara[i]=par[i];
		int ia = gFitParaIndexMap[i][0];
		int ib = gFitParaIndexMap[i][1];
		int ic = gFitParaIndexMap[i][2];
		int id = gFitParaIndexMap[i][3];
		int ie = gFitParaIndexMap[i][4];		
		gPara[ia][ib][ic][id][ie]=par[i];
#ifdef MulDFit_DEBUG 
		if(MulDFit_DEBUG>=5)
		{
			cout<<"myMulDPolN() set K_"<<ia<<ib<<ic<<id<<ie<<"="<<gPara[ia][ib][ic][id][ie]<<endl;
		}
#endif
	}


	double var=0;
	for(int ia=0;ia<=gMaxOrderA;ia++)
	{
		double itemA = pow(x[0],double(ia));
		for(int ib=0;ib<=gMaxOrderB;ib++)
		{
			double itemB = pow(x[1],double(ib));
			for(int ic=0;ic<=gMaxOrderC;ic++)
			{
				double itemC = pow(x[2],double(ic));
				for(int id=0;id<=gMaxOrderD;id++)
				{
					double itemD = pow(x[3],double(id));
					for(int ie=0;ie<=gMaxOrderE;ie++)
					{
						double itemE = pow(x[4],double(ie));
#ifdef MulDFit_DEBUG 
						if(MulDFit_DEBUG>=5)
						{
							cout<<"myMulDPolN():  K_"<<ia<<ib<<ic<<id<<ie<<"="<<(gPara)[ia][ib][ic][id][ie]<<endl;
						}
#endif
						var += gPara[ia][ib][ic][id][ie] * itemA * itemB *itemC *itemD *itemE;
					}
				}
			}
		}
	}
	return var;
}

//I am not use this global function any longer
//void myFcn(int &/*nPar*/, double */*grad*/, double &fval, double *par, int /*iflag */)
void myFcn(int &, double *, double &fval, double *par, int )
{
	double chi2 = 0;
	double var, var_cal;
	double weight = 1.0;  
	double x[5];

	for(unsigned int i=0;i<gVVar.size();i++)	 //calculate chisquare
	{	
		for(unsigned int iv=0;iv<5;iv++) x[iv]=gVSelfVar[iv][i];
		weight = gVWeight[i];
		var = gVVar[i];
		var_cal = myMulDPolN(x,par);
		chi2 += pow(var-var_cal,2.0)*weight*weight;
	}
	fval = chi2;
}


//////////////////////////////////////////////////////////////////////////////
MulDFit::MulDFit(const char* infile):ApplyMulDPolN(infile)
{
#ifdef MulDFit_DEBUG 
	if(MulDFit_DEBUG>0)
	{
		cout<<">>>>>MulDFit::MulDFit() started with debug level "<<MulDFit_DEBUG<<"<<<<<"<<endl;
	}
#endif

	//set the buffer to global pointer
	gParaNum = GetParaNum();
	gFitPara = GetFitPara();
	gFitParaErr = GetFitParaErr();
	gFitParaFixedFlag = GetFitParaFixedFlag();
	gFitParaIndexMap = GetFitParaIndexMap();
	gPara = GetMatrix(gMaxOrderA,gMaxOrderB,gMaxOrderC,gMaxOrderD,gMaxOrderE);
}

MulDFit::~MulDFit()
{
#ifdef MulDFit_DEBUG
	if(MulDFit_DEBUG>=3)
		cout<<"MulDFit::~MulDFit() "<<endl;
#endif
}


void MulDFit::LoadData()
{
#ifdef MulDFit_DEBUG
	if(MulDFit_DEBUG>=3)
		cout<<"MulDFit::LoadData()"<<endl;
#endif

	ifstream fin;
	fin.open("in_file.dat");
	char buf[512];
	//Line 1: variable names
	for(int i=0;i<7;i++) fin>>mDataVarName[i];
	fin.getline(buf,511);
	//Line 2: variable titles
	for(int i=0;i<7;i++) fin>>mDataVarTitle[i];
	fin.getline(buf,511);

	sprintf(mFitVarName,"%s_fit",this->GetBaseName(&mDataVarName[0][0]));
	sprintf(mFitVarTitle,"%s_fit",this->GetBaseName(&mDataVarTitle[0][0]));
	
#ifdef MulDFit_DEBUG
		if(MulDFit_DEBUG>=4)
		{
			for(int i=0;i<7;i++) cout<<setw(13)<<mDataVarName[i]<<" ";
			cout<<endl;
			for(int i=0;i<7;i++) cout<<setw(13)<<mDataVarTitle[i]<<" ";
			cout<<endl;
			cout<<"  Fitted_Variable:  name="<<mFitVarName<<"  title="<<mFitVarTitle<<endl;
		}
#endif


	double var,selfvar[5],weight;
	//now loop to the end
	while (!fin.eof())
	{
		fin>>var;
		for(int i=0;i<5;i++) fin>>selfvar[i];
		fin>>weight;

		//store them into the vector
		//sometimes fin will not check end of file, the last line is bad
		if(fabs(var-99999.0)>=1.0E-12 && fabs(weight-99999.0)>=1.0E-12)
		{
			if(!Cut(selfvar)) continue;

			gVVar.push_back(var);
			for(int i=0;i<5;i++) gVSelfVar[i].push_back(selfvar[i]);
			gVWeight.push_back(weight);
		}

#ifdef MulDFit_DEBUG
		if(MulDFit_DEBUG>=4)
		{
			cout<<setw(13)<<var<<" ";
			for(int i=0;i<5;i++) cout<<setw(13)<<selfvar[i]<<" ";
			cout<<setw(13)<<weight<<endl;
		}
#endif
		//reset the value
		var=weight=99999.0;
	}
	fin.close();

#ifdef MulDFit_DEBUG
	if(MulDFit_DEBUG>=1)
		cout<<"MulDFit::Load(): "<<gVVar.size()<<" data points are loaded"<<endl;
#endif
}



void MulDFit::WriteParaFile(const char* outfile)
{
#ifdef MulDFit_DEBUG
	if(MulDFit_DEBUG>=3)
		cout<<"MulDFit::WriteParaFile()"<<endl;
#endif

	ofstream fout;
	fout.open(outfile);
	//fout.setf(ios::scientific,ios::floatfield);

	//first 8 lines
	fout<<"#By Jixie Zhang, (jixie@jlab.org). Created by MulDFit.\n"
		<<"#This file defines the matrix which will be used by MulDPolN. \n"
		<<"#1) The first 8 lines is for comments and description. \n"
		<<"#2) The 9th line is to specify the maximum order for variables:\n#";
	for(int i=1;i<6;i++) fout<<setw(13)<<mDataVarName[i]<<" ";
	fout<<endl;
	fout<<"#3) The 10th line is the title of each collum of matrix elements at this file. \n"
		<<"#4) This matrix is fitted from a data set defined by next line:\n#";
	for(int i=0;i<6;i++) fout<<setw(13)<<mDataVarTitle[i]<<" ";
	fout<<endl;


	//line 9
	fout<<MaxOrderA
		<<" "<<setw(4)<<MaxOrderB
		<<" "<<setw(4)<<MaxOrderC
		<<" "<<setw(4)<<MaxOrderD
		<<" "<<setw(4)<<MaxOrderE
		<<endl;

	//line 10
	fout<<"Label ^b   ^c   ^d   ^e";
	for(int i=0;i<=MaxOrder;i++) fout<<"       K_"<<i<<"bcde";
	for(int i=0;i<=MaxOrder;i++) fout<<" fix"<<i;
	fout<<endl;


	int index=0;
	char label;
	int  ib=0,ic=0,id=0,ie=0;

	//now loop the parameter buffer to the end
	while(index<mParaNum) 
	{
		label=mFitParaIndexMap[index][5];

		int thisOrder=0;
		if(label=='A')  
		{
			thisOrder = MaxOrderA;
			ib=mFitParaIndexMap[index][1];
			ic=mFitParaIndexMap[index][2];
			id=mFitParaIndexMap[index][3];
			ie=mFitParaIndexMap[index][4];
		}
		else if(label=='B')  
		{
			thisOrder = MaxOrderB;
			ib=mFitParaIndexMap[index][0];
			ic=mFitParaIndexMap[index][2];
			id=mFitParaIndexMap[index][3];
			ie=mFitParaIndexMap[index][4];
		}
		else if(label=='C')  
		{
			thisOrder = MaxOrderC;
			ib=mFitParaIndexMap[index][0];
			ic=mFitParaIndexMap[index][1];
			id=mFitParaIndexMap[index][3];
			ie=mFitParaIndexMap[index][4];
		}
		else if(label=='D')  
		{
			thisOrder = MaxOrderD;
			ib=mFitParaIndexMap[index][0];
			ic=mFitParaIndexMap[index][1];
			id=mFitParaIndexMap[index][2];
			ie=mFitParaIndexMap[index][4];
		}
		else if(label=='E')  
		{
			thisOrder = MaxOrderE;
			ib=mFitParaIndexMap[index][0];
			ic=mFitParaIndexMap[index][1];
			id=mFitParaIndexMap[index][2];
			ie=mFitParaIndexMap[index][3];
		}


		fout<<char(label)<<"  "
			<<" "<<setw(4)<<ib
			<<" "<<setw(4)<<ic
			<<" "<<setw(4)<<id
			<<" "<<setw(4)<<ie;

		//write the fitted values
		for(int i=0;i<=MaxOrder;i++) 
		{
			if(i<=thisOrder) fout<<" "<<setw(13)<<mFitPara[index+i];
			else fout<<" "<<setw(13)<<" ";
		}
		//write the fixedflag
		for(int i=0;i<=MaxOrder;i++) 
		{
			if(i<=thisOrder) fout<<" "<<setw(4)<<mFitParaFixedFlag[index++];
			else fout<<" "<<setw(4)<<" ";
		} 
		fout<<endl;

	}

	fout.close();
}



void MulDFit::WriteParaFile_onelabel(ofstream &fout, char label, int orderA, 
									 int orderB,int orderC,int orderD, int orderE)
{
#ifdef MulDFit_DEBUG
	if(MulDFit_DEBUG>=3)
		cout<<"MulDFit::WriteParaFile_oneline()"<<endl;
#endif

	if(!fout.good())
	{
		cout<<"MulDFit::WriteParaFile_onelabel() can not open output file!"<<endl; 
		return;
	}

	int tmpMaxOrder = max(max(max(orderA,orderB),max(orderC,orderD)),orderE);

	if(orderA<1) return;

	//fout.setf(ios::scientific,ios::floatfield);

	//write A
	for(int ib=0;ib<=orderB;ib++)
	{
		for(int ic=0;ic<=orderC;ic++)
		{
			for(int id=0;id<=orderD;id++)
			{
				for(int ie=0;ie<=orderE;ie++)
				{
					if(ib+ic+id+ie>5) continue;
					fout<<label<<"  "
						<<" "<<setw(4)<<ib
						<<" "<<setw(4)<<ic
						<<" "<<setw(4)<<id
						<<" "<<setw(4)<<ie;

					//write the fitted values
					for(int ia=0;ia<=tmpMaxOrder;ia++) 
					{
						if(ia<=orderA)
						{
						//fout<<" "<<setw(13)<<0.0;
						double tmp = 0.0 +(((ia*(orderB+1) + ib)*(orderC+1) + ic)*(orderD+1) +id)*(orderE+1) + ie;
						fout<<" "<<setw(13)<<tmp*1.0E-3;
						}
						else fout<<" "<<setw(13)<<" ";
					}
					//write the fixedflag
					for(int ia=0;ia<=tmpMaxOrder;ia++) 
					{
						if(ia<=orderA) fout<<" "<<setw(4)<<0;
						else fout<<" "<<setw(4)<<" ";
					}
					fout<<endl;
				}
			}
		}
	}

}


void MulDFit::Write1stParaFile(char label, int orderA,int orderB,int orderC,int orderD, int orderE)
{
#ifdef MulDFit_DEBUG
	if(MulDFit_DEBUG>=3)
		cout<<"MulDFit::WriteParaFile()"<<endl;
#endif
	int tmpMaxOrder = max(max(max(orderA,orderB),max(orderC,orderD)),orderE);

	ofstream fout;
	fout.open("para_template.ini");
	cout<<"Creating input file template: para_template.ini "<<endl;
	//fout.setf(ios::scientific,ios::floatfield);

	
	const char *pDataVarName[6]={"Y","A","B","C","C","E"};
	const char *pDataVarTitle[6]={"Y","A","B","C","C","E"};

	//first 8 lines
	fout<<"#By Jixie Zhang, (jixie@jlab.org). Created by MulDFit.\n"
		<<"#This file defines the matrix which will be used by MulDPolN. \n"
		<<"#1) The first 8 lines is for comments and description. \n"
		<<"#2) The 9th line is to specify the maximum order for variables:\n#";
	for(int i=1;i<6;i++) fout<<setw(13)<<pDataVarName[i]<<" ";
	fout<<endl;
	fout<<"#3) The 10th line is the title of each collum of matrix elements at this file. \n"
		<<"#4) This matrix is fitted from a data set defined by next line:\n#";
	for(int i=0;i<6;i++) fout<<setw(13)<<pDataVarTitle[i]<<" ";
	fout<<endl;

	//lien 8
	fout<<orderA
		<<" "<<setw(4)<<orderB
		<<" "<<setw(4)<<orderC
		<<" "<<setw(4)<<orderD
		<<" "<<setw(4)<<orderE
		<<endl;

	//line 9
	fout<<"Label ^b   ^c   ^d   ^e";
	for(int i=0;i<=tmpMaxOrder;i++) fout<<"       K_"<<i<<"bcde";
	for(int i=0;i<=tmpMaxOrder;i++) fout<<" fix"<<i;
	fout<<endl;
	if(label=='A' || label=='Z') WriteParaFile_onelabel(fout, 'A', orderA, orderB, orderC, orderD, orderE);
	if(label=='B' || label=='Z') WriteParaFile_onelabel(fout, 'B', orderB, orderA, orderC, orderD, orderE);
	if(label=='C' || label=='Z') WriteParaFile_onelabel(fout, 'C', orderC, orderA, orderB, orderD, orderE);
	if(label=='D' || label=='Z') WriteParaFile_onelabel(fout, 'D', orderD, orderA, orderB, orderC, orderE);
	if(label=='E' || label=='Z') WriteParaFile_onelabel(fout, 'E', orderE, orderA, orderB, orderC, orderD);

	fout.close();
}


const char* MulDFit::GetBaseName(const char* str)
{
	const char *name = strrchr(str,'/');
	if(name==NULL) name=&str[0];
	else name++; 	
	//printf ("\"%s\": Last occurence of '/' found at %d, name=%s\n",str,int(name-str),name);

	const char *surfix = strrchr(name,'_');
	if(surfix==NULL)  surfix=&name[sizeof(name)];
	//printf ("\"%s\": Last occurence of '_' found at %d \n",name,int(surfix-name+1));
	
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


//Get lower and upper limit of TH1 which have continuous 4 bins above the given limit 
//if pYcut<0, then will use 2% of maximum height 
//if pYcut>0, then require this value
//void   GetTH1BoundaryValues(TH1 *h1,double &pXStart,double &pXEnd,double pYcut=-1.0);
void MulDFit::GetTH1BoundaryValues(TH1 *h1,double &pXStart,double &pXEnd,double pYcut)
{
  // Get the start bin and end bin 
  int pNXBin=h1->GetNbinsX();
  int pXBinStart=0,pXBinEnd=-1;

  if(pYcut<=0) pYcut=0.02*h1->GetMaximum(); 
  for(int i=1;i<=pNXBin;i++)
    {
      //continuous 4 bins none zero
      if(pXBinStart<=0 && i+3<=pNXBin)
	{
	  if (h1->GetBinContent(i)>pYcut   && h1->GetBinContent(i+1)>0 &&
	      h1->GetBinContent(i+2)>0 && h1->GetBinContent(i+3)>0 )
	    pXBinStart=i;
	}  

      if(pXBinEnd<=0 && pNXBin-i-2>=1)
	{
	  if (h1->GetBinContent(pNXBin-i+1)>pYcut && h1->GetBinContent(pNXBin-i)>0 &&
	      h1->GetBinContent(pNXBin-i-1)>0 && h1->GetBinContent(pNXBin-i-2)>0)
	    pXBinEnd=pNXBin-i+1;
	}
      if(pXBinStart>0 && pXBinEnd>0) break;
    }
  if(pXBinEnd<=0) pXBinEnd=pNXBin;

  pXStart=h1->GetBinLowEdge(pXBinStart);
  pXEnd=h1->GetBinLowEdge(pXBinEnd)+h1->GetBinWidth(pXBinEnd);
 
  double pYmax=h1->GetMaximum();	
  TLine *L1=new TLine(pXStart,0,pXStart,pYmax);
  L1->SetLineColor(2);
  h1->GetListOfFunctions()->Add(L1); 
  TLine *L2=new TLine(pXEnd,0,pXEnd,pYmax);
  h1->GetListOfFunctions()->Add(L2); 
  L2->SetLineColor(2);
}


TF1* MulDFit::FitGaus(TH1* h1, double &mean, double &sigma, double range_in_sigma)
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


void MulDFit::PlotDelta(TTree* tree)
{
#ifdef MulDFit_DEBUG
	if(MulDFit_DEBUG>=3)
		cout<<"MulDFit::PlotDelta()"<<endl;
#endif
	system("mkdir -p Graph");
	const double deg = atan(1.0)/45.;

	char var_basename[100];
	char key[100];
	sprintf(var_basename,"%s",this->GetBaseName(&mDataVarName[0][0]));
	sprintf(key,"RTPC12_%s",var_basename);

	////////////////////////////////////////////////////////////
	char  varName[100], varTitle[100];
	sprintf(varName,"d%s",var_basename);
	sprintf(varTitle,"#delta%s",var_basename);
	/////////////////////////////////////////////////////////////
	//I have to put a wide range here otherwise I can not extract the sigma from the distribution
	TCut NoCut = "abs(Y_fit-Y)/Y<1.0";
	cout<<"TCut NoCut = \""<<NoCut<<"\""<<endl;

	gStyle->SetStatX(1.0-gStyle->GetPadRightMargin());  
	gStyle->SetStatH(0.060);
	gStyle->SetStatW(0.200);

	/////////////////////////////////////////////////////////////
	//config the variables
	//plot 2-D dX, 4 pads:  
	//1) dX vs P0  
	//2) dX vs Theta0  
	//1) dX vs Phi0  
	//1) dX vs Z0
	const int nDelta=1;
	const int nVar=6;
	const char *strYVar[]={"(Y_fit-Y)"};
	const char *strYName[]={"dY"};
	const char *strYTitle[]={"#deltaY"};
	strYName[0]= &varName[0];
	strYTitle[0] = &varTitle[0];

	int    YBinNum[]={50};
	double YBinMin[]={-0.05};
	double YBinMax[]={ 0.05};

	double YSigmaFactor[nDelta]={5.0};   //sigma factor to cut on Y when plot 2-D 
	double YMean[nDelta],YSigma[nDelta];
	char   DeltaCut[512],strYCut[nDelta][255]; 
	DeltaCut[0]='\0'; //I have to give this terminal sign to this buffer, otherwise strlen() will not work

	const char *strXVar[]={"Y","A","B","C","D","E"};
	const char *strXName[]={"Y","A","B","C","D","E"};
	const char *strXTitle[]={"Y","A","B","C","D","E"};
	for(int i=0;i<6;i++) 
	{
		strXName[i] = &mDataVarName[i][0];
		strXTitle[i] = &mDataVarTitle[i][0];
	}
	int    XBinNum[]={  60,  50,     60,      36,   40,  50};
	double XBinMin[]={ 0.0,  0.,  0*deg,-180*deg,-200.,-100};
	double XBinMax[]={ 0.3,200.,180*deg, 180*deg, 200., 100};

	char hTitle[255], strTg[255], hName[100], hName_1[100], hName_2[100];

	/////////////////////////////////////////////////////////
	TH1F *h1=0;
	TH2F *h2=0;
	TH1F *hMean=0, *hSigma=0,*hSigma_neg=0;
	
	/////////////////////////////////////////////////////////
	//Plot 1-D X-Variables 
	//determine the range for all X-Variables
	int UseGivenRange=0;
	if(!UseGivenRange)
	{
		TCanvas *c30 = new TCanvas("c30","",900,800);
		c30->Clear();
		int nCol=int(ceil(nVar/3.0));
		int nRow=int(ceil(double(nVar)/double(nCol)));
		c30->Divide(nCol,nRow);
		for(int vv=0;vv<nVar;vv++)
		{	
			c30->cd(vv+1);//gPad->SetRightMargin(0.12);
			sprintf(hName,"%s",strXName[vv]);
			sprintf(hTitle,"%s ;%s",strXName[vv],strXTitle[vv]);
			sprintf(strTg,"%s",strXVar[vv]);

			//cout<<"Name = "<<hName<<"    var = "<<strTg<<endl;
			//cout<<"Title = "<<hTitle<<endl;

			//iteration one, get initial mean and rms then plot again

			h1 = (TH1F*) gROOT->FindObject(hName);
			if(h1) delete h1;
			sprintf(strTg,"%s >> %s",strXVar[vv],hName);
			tree->Draw(strTg,NoCut);
			h1 = (TH1F*) gROOT->FindObject(hName);
			h1->SetTitle(hTitle); 

			double pXStart=0., pXEnd=0.; 
			GetTH1BoundaryValues(h1,pXStart,pXEnd);
			h1->Draw();

			XBinMin[vv]=pXStart;
			XBinMax[vv]=pXEnd;
		}

		c30->cd(0);
		c30->Modified();
		c30->SaveAs(Form("Graph/Var_1D_%s.png",key));
	}


	/////////////////////////////////////////////////////////
	//Plot 1-D delta
	TCanvas *c31 = new TCanvas("c31","",600,500);
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

		if(UseGivenRange)
		{
			//if you know the range
			h1 = new TH1F(hName,hTitle,YBinNum[dd],YBinMin[dd],YBinMax[dd]);
			tree->Project(hName,strTg,NoCut);                                                                                  
			h1->Draw();
		}
		else
		{
			//if you do not know the range
			h1 = (TH1F*) gROOT->FindObject(hName);
			//iteration one, get initial mean and rms then plot again
			
			if(h1) delete h1;
			sprintf(strTg,"%s >> %s",strYVar[dd],hName);
			tree->Draw(strTg,NoCut);
			h1 = (TH1F*) gROOT->FindObject(hName);
			double Mean = h1->GetMean();
			double RMS = h1->GetRMS();
			delete h1;
			sprintf(strTg,"%s >> %s(120,%.5f,%.5f)",strYVar[dd],hName,Mean-3*RMS,Mean+3*RMS);
			tree->Draw(strTg,NoCut);
			h1 = (TH1F*) gROOT->FindObject(hName);
			double Sigma;
			FitGaus(h1,Mean,Sigma);
			double pCut=(YSigmaFactor[dd]+1.0)*Sigma;
			
			//cout<<hName<<":  Mean="<<Mean<<"  RMS="<<h1->GetRMS()<<"  Sigma="<<Sigma<<endl;
			//iteration two, plot again using sigmato get better range 
			delete h1;
			sprintf(strTg,"%s >> %s(100,%.5f,%.5f)",strYVar[dd],hName,Mean-pCut,Mean+pCut);
			tree->Draw(strTg,NoCut);
			h1 = (TH1F*) gROOT->FindObject(hName);

			h1->SetTitle(hTitle); 
		}
		
		h1->SetLineColor(4); 
		FitGaus(h1,YMean[dd],YSigma[dd]); 
		YBinNum[dd]=60;
		YBinMin[dd]=YMean[dd]-6.0*YSigma[dd];
		YBinMax[dd]=YMean[dd]+6.0*YSigma[dd];
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
			tree->Project(hName,strTg,DeltaCut);
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


void MulDFit::CheckResult()
{
#ifdef MulDFit_DEBUG
	if(MulDFit_DEBUG>=3)
		cout<<"MulDFit::CheckResult()"<<endl;
#endif

	//fill root tree then plot delta
	//tree buffer
	int    Index;
	double Y,A,B,C,D,E, Weight;
	double Y_fit;

	TFile *fFile = new TFile("muldfit_out.root","RECREATE");

	TTree *tree = new TTree("t","tree of MulDFit, Y=MulDPolN(x,par)");
	tree->Branch("Index",&Index,"Index/I");
	tree->Branch("Y",&Y,"Y/D");
	tree->Branch("A",&A,"A/D");
	tree->Branch("B",&B,"B/D");
	tree->Branch("C",&C,"C/D");
	tree->Branch("D",&D,"D/D");
	tree->Branch("E",&E,"E/D");
	tree->Branch("Y_fit",&Y_fit,"Y_fit/D");

	for(unsigned int i=0;i<gVVar.size();i++)
	{
		Index=i;
		Y=gVVar[i];
		A=gVSelfVar[0][i];
		B=gVSelfVar[1][i];
		C=gVSelfVar[2][i];
		D=gVSelfVar[3][i];
		E=gVSelfVar[4][i];
		Weight=gVWeight[i];

		double x[5]={A,B,C,D,E};		
		Y_fit = ApplyMulDPolN::Eval(x);

#ifdef MulDFit_DEBUG
	if(MulDFit_DEBUG>=4)
		cout<<"entry i="<<i<<"  Y="<<Y<<"  Y_fit="<<Y_fit<<endl;
#endif
		tree->Fill();
	}
	fFile->Write("",TObject::kOverwrite);

	PlotDelta(tree);

	fFile->Close();
}

//inline
bool MulDFit::Cut(double *)
{
	return true;
}

void MulDFit::Process()
{
#ifdef MulDFit_DEBUG
	if(MulDFit_DEBUG>=3)
	{
		cout<<"MulDFit::Process()\n";
	}
#endif

	LoadData();
	if(gVVar.size()<=(unsigned int)mParaNum)
	{
		cout<<"MulDFit::Process(): not enough data points to do Minuit fit!"<<endl;
		return;
	}

	for(int i=0;i<2;i++) 
	{
		DoMinuitFit();
		char tmpout[100];
		sprintf(tmpout,"_tmp_out_para_%d.ini",i);
		WriteParaFile(tmpout);
	}
	DoMinuitFit();
	WriteParaFile();

	CheckResult();
}


void MulDFit::DoMinuitFit()
{
#ifdef MulDFit_DEBUG
	if(MulDFit_DEBUG>=3)
	{
		cout<<"MulDFit::DoMinuitFit()\n";
	}
#endif

	//-----------------------------------------------------------------------//
	//initialize TMinuit with a maximum of mParaNum params
	TMinuit *pMinuit = new TMinuit(mParaNum);  
	pMinuit->SetFCN(myFcn);

	Double_t arglist[10];
	Int_t ierflg = 0;

	//-----------------------------------------------------------------------//
	/// set print level
	/*
	SET PRIntout  <level>
	Sets the print level, determining how much output will be
	produced. Allowed values and their meanings are displayed
	after a SHOw PRInt command, and are currently <level>=:
	[-1]  no output except from SHOW commands
	[0]  minimum output
	[1]  default value, normal output
	[2]  additional output giving intermediate results.
	[3]  maximum output, showing progress of minimizations.
	*/
	///set print level
	arglist[0] = -1;
#ifdef MulDFit_DEBUG
	if(MulDFit_DEBUG>=1)
	{
		arglist[0] = MulDFit_DEBUG-2;
	}
#endif
	pMinuit->mnexcm("SET PRIntout",arglist,1,ierflg);

	//-----------------------------------------------------------------------//
	/*
	SET NOWarnings
	Supresses Minuit warning messages.
	SET WARnings
	Instructs Minuit to output warning messages when suspicious
	conditions arise which may indicate unreliable results.
	This is the default.
	*/
	arglist[0] = 0;
	pMinuit->mnexcm("SET NOWarnings",arglist,0,ierflg);

#ifdef MulDFit_DEBUG
	if(MulDFit_DEBUG>=3)
	{
		arglist[0] = 1;
		pMinuit->mnexcm("SET WARnings",arglist,0,ierflg);
	}
#endif

	//-----------------------------------------------------------------------//
	/// Set starting values and step sizes for parameters
	/// use last fit result as the start point and use step size 0.000001
	for(int ii=0;ii<mParaNum;ii++)
	{
		char strParaName[100];
		sprintf(strParaName,"K_%d%d%d%d%d",mFitParaIndexMap[ii][0],mFitParaIndexMap[ii][1],
			mFitParaIndexMap[ii][2],mFitParaIndexMap[ii][3],mFitParaIndexMap[ii][4]);
		pMinuit->mnparm(ii, strParaName, mFitPara[ii], 0.001, 0,0,ierflg);
	}

	arglist[0] = 1;
	pMinuit->mnexcm("SET ERR", arglist ,1,ierflg);

	//-----------------------------------------------------------------------//
	/*
	//SET LIMits  0->[parno]  1->[lolim]  2->[uplim]
	//	  arglist[0] = 0;	   //this means set all parameters with the same limit
	//	  arglist[0] = 1;	   //this means set 1st parameters with the specified limit
	//	  arglist[0] = 2;	   //this means set 2nd parameters with the sepecified limit
	arglist[0] = 0;
	arglist[1] = 0.;
	arglist[2] = 0.;
	pMinuit->mnexcm("SET LIMits",arglist,3,ierflg);
	*/

	bool bFitWithLimit=false;
	if(bFitWithLimit)
	{
		for (int ii = 0; ii < mParaNum; ii++)
		{	
			//Set Min of Parameters, value==0 for No-Bound
			double lowlimit=0,uplimit=0;
			arglist[0] = 0;    
			arglist[1] = lowlimit;
			arglist[2] = uplimit;
			pMinuit->mnexcm("SET LIMits",arglist,3,ierflg);
		}
	}

	//-----------------------------------------------------------------------//
	//fix some parameters

	//Fix parameters in fortran kernal
	//arglist[0] = ii+1;   
	//pMinuit->mnexcm("FIX ", arglist ,1,ierflg);
	//or in C++ wrapper
	for (int ii = 0; ii < mParaNum; ii++)
	{	
		if(mFitParaFixedFlag[ii]==1) pMinuit->FixParameter(ii);
	}

	//-----------------------------------------------------------------------//
	// Now ready for minimization step
	arglist[0] = 2500;
	arglist[1] = 1.;// tolerance
	pMinuit->mnexcm("MIGRAD", arglist ,2,ierflg);

	//-----------------------------------------------------------------------//
	//get result
	for (int ii = 0; ii < mParaNum; ii++)
	{
		pMinuit->GetParameter(ii,mFitPara[ii],mFitParaErr[ii]);
	}


	//-----------------------------------------------------------------------//
	/*
	// Print results
	Double_t amin,edm,errdef;
	Int_t nvpar,nparx,icstat;
	pMinuit->mnstat(amin,edm,errdef,nvpar,nparx,icstat);
	pMinuit->mnprin(3,amin);
	*/


	//-----------------------------------------------------------------------//
	// free the memory
	delete pMinuit;
	//-----------------------------------------------------------------------//
}

