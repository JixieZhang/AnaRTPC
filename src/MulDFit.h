//////////////////////////////////////////////////////////////////////
// MulDFit.h: interface for the MulDFit class.
// Created by jixie zhang ,9/10/2014
//This class is trying to fit multiple dimention polenomial fuction
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

#ifndef _MulDFit_
#define _MulDFit_

#define MulDFit_DEBUG 2

#include <stdio.h>
#include <string>
#include <fstream>
#include <math.h>
#include <iostream>
#include <iomanip>      // std::setw
#include <algorithm>    // std::max
#include <vector>
using namespace std;

#include "ApplyMulDPolN.h"

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~`
class TF1;
class TH1;
class TTree;

class MulDFit:ApplyMulDPolN
{
public:
	MulDFit(const char* infile="in_para.ini");
	virtual ~MulDFit();

public:
	void   Process();
	void   LoadData();
	void   CheckResult();

	static double MulDPolN(double *x, double *par); //fit function 
	static void   WriteParaFile_onelabel(ofstream &fout, char lable, int orderA, 
		int orderB,int orderC,int orderD, int orderE);
	static void   Write1stParaFile(char label,int orderA,int orderB,int orderC,int orderD, int orderE);

private:
	bool   Cut(double *);
	void   WriteParaFile(const char* outfile="out_para.ini");
	void   DoMinuitFit();
	void   PlotDelta(TTree* tree);

private:
	//tools only, maybe I should move them to another class 
	//return the basename of a variable: remove whatever appends to last underscore
	//for example:  1) P_fit -->P, 2) P->P
	const char* GetBaseName(const char* str);

	//Get lower and upper limit of TH1 which have continuous 4 bins above the given limit 
	//if pYcut<0, then will use 2% of maximum height 
	//if pYcut>0, then require this value
	void   GetTH1BoundaryValues(TH1 *h1,double &pXStart,double &pXEnd,double pYcut=-1.0);

	//DO gaussian fit to histo and return its mean and sigma
	TF1*   FitGaus(TH1* h1, double &mean, double &sigma, double range_in_sigma=1.0);
	

private:

	char mDataVarName[7][100];
	char mDataVarTitle[7][100];
	char mFitVarName[100];
	char mFitVarTitle[100];

};
#endif //_MulDFit_

