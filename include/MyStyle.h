#include <math.h>
#include <string>
#include <time.h>
#include <stdio.h>
#include <vector>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <TROOT.h>
#include <TStyle.h>
#include <TF1.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TFile.h>
#include <TDirectory.h>
#include <TObject.h>
#include <TVector3.h>
#include <TLorentzVector.h>
#include <TLorentzRotation.h>
#include <TCanvas.h>
#include <TPad.h>
#include <TProfile.h>
#include <TChain.h>
#include <TBenchmark.h>
#include <TAxis.h>
#include <TLegend.h>
#include <TLegendEntry.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TPaveStats.h>
using namespace std;

void SetMyFitStyle(bool bStatAt2ndLine=true);

void SetMyFitStyle(bool bStatAt2ndLine)
{
    TStyle *MyFitStyle=gStyle;
    MyFitStyle->SetPalette(1);

    MyFitStyle->SetPadLeftMargin(0.12);  
    //MyFitStyle->SetPadRightMargin(0.05); 
    MyFitStyle->SetPadBottomMargin(0.12); 

    //histo title
    //MyFitStyle->SetTitleX(MyFitStyle->GetPadLeftMargin()+0.04);
    //MyFitStyle->SetTitleY(0.99);
    if(bStatAt2ndLine) MyFitStyle->SetTitleW(0.85);  
    else MyFitStyle->SetTitleW(0.67);  //comment it to let it change automaticly
    MyFitStyle->SetTitleH(0.08);
    MyFitStyle->SetTitleBorderSize(0);
    MyFitStyle->SetTitleTextColor(1);
    //if(gROOT->IsBatch()) 
		MyFitStyle->SetTitleStyle(4000); 

    //stat
    //MyFitStyle->SetOptStat(10); //keep entries only
    MyFitStyle->SetOptStat(0);
    //MyFitStyle->SetOptFit(001);  //show x^2/ndf
    MyFitStyle->SetOptFit(0010); //not show x^2/ndf
    
	gStyle->SetStatX(1.0-gStyle->GetPadRightMargin());
    if(bStatAt2ndLine) 
    {
		gStyle->SetStatY(1.0-gStyle->GetPadTopMargin()); 
        MyFitStyle->SetStatW(0.17);
    }
    else 
    {
        MyFitStyle->SetStatY(0.99);
        MyFitStyle->SetStatW(0.16);
    }
    //MyFitStyle->SetStatH(MyFitStyle->GetTitleH()*1.33); //much better than 0.08==>0.1, 0.09==>0.12
    //MyFitStyle->SetStatFontSize(0.06);
	if(MyFitStyle->GetOptFit()==8 || MyFitStyle->GetOptFit()==9)  
		MyFitStyle->SetStatH(0.08);                                
    else   
		MyFitStyle->SetStatH(0.16); 
    MyFitStyle->SetStatTextColor(1);
    //if(gROOT->IsBatch()) 
		MyFitStyle->SetStatStyle(4000);

    //histo line and marker
    //1,6,7 small dot, medium dot and large dot
    //2 + ; 3 * ;4 o ;5 x
    //20 cylinder, 21 Square; 22 up_Triangle; 23 Down_triangle
    //29 Solid_5_Star; 30 Empty_5_Star

    //MyFitStyle->SetMarkerStyle(20);
    //MyFitStyle->SetMarkerColor(1);
    MyFitStyle->SetHistLineWidth(2);
    MyFitStyle->SetFuncWidth(3);   //set all func line width to 3
    //MyFitStyle->SetFuncColor(2); // //I may change this for peticular function line

    //pads and canvas
    MyFitStyle->SetPadGridX(1);
    MyFitStyle->SetPadGridY(1);
    MyFitStyle->SetPadBorderMode(0);
    MyFitStyle->SetFrameBorderMode(0);
    MyFitStyle->SetCanvasBorderMode(0);
    MyFitStyle->SetPadTickX(1);
    MyFitStyle->SetPadTickY(1);
    MyFitStyle->SetNdivisions(505,"X");
    MyFitStyle->SetNdivisions(505,"Y");


    //set the draw-option "text" format
    MyFitStyle->SetPaintTextFormat(".0f");

    //Fill area color
    MyFitStyle->SetFrameFillColor(0);
    MyFitStyle->SetTitleFillColor(0);
    MyFitStyle->SetStatColor(0);
    //MyFitStyle->SetHistFillColor(0);  //I may change this for peticular histogram
    MyFitStyle->SetCanvasColor(0);

    //all text in the canvas or histo
    MyFitStyle->SetTextColor(1);	//set the default text color to black 
    MyFitStyle->SetTextSize(0.06);	//don't know what it will do yet 


    //axis title
    MyFitStyle->SetTitleOffset(0.7,"X"); 
    MyFitStyle->SetTitleOffset(0.8,"Y"); 
    MyFitStyle->SetTitleSize(0.07,"X");
    MyFitStyle->SetTitleSize(0.07,"Y");
    MyFitStyle->SetTitleColor(1,"X");
    MyFitStyle->SetTitleColor(1,"Y");

    //axis label
    MyFitStyle->SetLabelOffset(0.001,"X");
    MyFitStyle->SetLabelOffset(0.001,"Y");
    MyFitStyle->SetLabelSize(0.06,"X");
    MyFitStyle->SetLabelSize(0.06,"Y");
    //MyFitStyle->SetLabelColor(1,"X"); 

    // The position of the date string can be controlled by:
    //  optdate = 10*format + mode
    //    mode = 1   (default) date is printed in the bottom/left corner.
    //    mode = 2   date is printed in the bottom/right corner.
    //    mode = 3   date is printed in the top/right corner.
    //    format = 0 (default) date has the format like: "Wed Sep 25 17:10:35 2002"
    //    format = 1 date has the format like: "2002-09-25"
    //    format = 2 date has the format like: "2002-09-25 17:10:35"
    //mutipad canvas will hide the date in the back
    //MyFitStyle->SetOptDate(11);
  
    gROOT->SetStyle(MyFitStyle->GetName());
}

void SavePad(char* name) 
{
	char name_eps[20];
	sprintf(name_eps,"%s.eps",name);
	char name_gif[20];
	sprintf(name_gif,"%s.gif",name);
	char name_root[20];
	sprintf(name_root,"%s.root",name);
	char name_C[20];
	sprintf(name_C,"%s.C",name);
	gPad->Print(name_eps);
	gPad->SaveAs(name_gif);
	//gPad->SaveAs(name_root);
	//gPad->SaveAs(name_C);
}

void SaveCanvas(char* name) 
{
	char name_eps[20];
	sprintf(name_eps,"%s.eps",name);
	char name_gif[20];
	sprintf(name_gif,"%s.gif",name);
	char name_root[20];
	sprintf(name_root,"%s.root",name);
	char name_C[20];
	sprintf(name_C,"%s.C",name);
	gPad->GetCanvas()->Print(name_eps);
	gPad->GetCanvas()->SaveAs(name_gif);
	//gPad->SaveAs(name_root);
	//gPad->SaveAs(name_C);
}

void SaveAs(char* name) 
{
	SaveCanvas(name); 
}
void Save(char* name) 
{
	SaveCanvas(name); 
}
void Print(char* name) 
{
	SaveCanvas(name); 
}

void SaveAll(char* name) 
{
	char name_eps[20];
	sprintf(name_eps,"%s.eps",name);
	char name_gif[20];
	sprintf(name_gif,"%s.gif",name);
	char name_root[20];
	sprintf(name_root,"%s.root",name);
	char name_C[20];
	sprintf(name_C,"%s.C",name);
	gPad->GetCanvas()->Print(name_eps);
	gPad->GetCanvas()->SaveAs(name_gif);
	gPad->SaveAs(name_root);
	gPad->SaveAs(name_C);
}