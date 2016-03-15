// ApplyMulDPolN.h: interface for the ApplyMulDPolN class.
// Created by jixie zhang ,9/10/2014
//This class is trying apply the fited multiple dimention polenomial fuction
//  V = Sum_abcde {K_abcde * X^a T^b Y^C F^d P^e}
//    = Sum_bcde {Sum_a{K_abcde * X^a} * T^b Y^C F^d P^e}
//
//The initial value for all parameters should be given by a parameter file
//named in_para.ini. The format of this file is 
//line 1: MaxOrderA   MaxOrderB MaxOrderC MaxOrderD MaxOrderE
//line 2: #lable of collum, which looks like the next line: 
//5     5     5     5     5
//Label ^b   ^c   ^d   ^e   K_0bcde       K_1bcde       K_2bcde       K_3bcde       K_4bcde       K_5bcde     0fixed 1fixed 2fixed 3fixed 4fixed 5fixed
//A      0    0    0    0   0.000000e+00  0.000000e+00  0.000000e+00  0.000000e+00  0.000000e+00  0.000000e+00     0      0      0      0      0      0
//
//////////////////////////////////////////////////////////////////////

#ifndef _ApplyMulDPolN_
#define _ApplyMulDPolN_

#define ApplyMulDPolN_DEBUG 3

#include <stdio.h>
#include <string>
#include <fstream>
#include <math.h>
#include <iostream>
#include <iomanip>      // std::setw
#include <algorithm>    // std::max
#include <vector>
using namespace std;

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~`


class ApplyMulDPolN
{
public:
	ApplyMulDPolN(const char* infile="in_para.ini");
	virtual ~ApplyMulDPolN();

public:
	double Eval(double *x);

	template <typename T> static void FreeDynArray(T *****arr, int bin1, int bin2, int bin3,int bin4);
	template <typename T> static T ***** IniDynamicArray(int bin1, int bin2, int bin3, int bin4,int bin5, T defautval);
	template <typename T> static void DebugArray(T ***** arr, int bin1, int bin2, int bin3, int bin4,int bin5);

	static double MulDPolN(double *x, double *****matrix, int bin1, int bin2, int bin3, int bin4, int bin5); 
	static void   SetMatrix(double *par, int parnum, double **parmap, double *****matrix, int bin1, int bin2, int bin3, int bin4, int bin5);

public:
	void    SetMatrix(double *par);
	void    SetMatrix(double *****matrix,int bin1, int bin2, int bin3, int bin4, int bin5);

	int     GetParaNum() {return mParaNum;};
	double* GetFitPara() {return mFitPara;};
	double* GetFitParaErr() {return mFitParaErr;};
	int*    GetFitParaFixedFlag() {return mFitParaFixedFlag;};
	int**   GetFitParaIndexMap() {return mFitParaIndexMap;};

	double***** GetMatrix() {return mPara;};
	double***** GetMatrix(int &bin1, int &bin2, int &bin3, int &bin4, int &bin5);

private:
	double MulDPolN(double *x);  
	void   ReadParaFile(const char* infile="in_para.ini");

public:
	//since minuit can only access global variable, I have to make these variable static
	int MaxOrderA;
	int MaxOrderB;
	int MaxOrderC;
	int MaxOrderD;
	int MaxOrderE;
	int MaxOrder;  //store the maximum of all the above 5

	//double mPara[kMaxParaA][kMaxParaB][kMaxParaC][kMaxParaD][kMaxParaE];
	double *****mPara;
	double *mFitPara;
	double *mFitParaErr;
	int    *mFitParaFixedFlag;
	//An index map that pointing 1-D parameter to 5-D array elements 
	//mFitParaIndexMap[mParaNum][6], from 0 to 6 it is indexA...indexE and label
	int    **mFitParaIndexMap; 
	int    mParaNum;  //number of parameter

};
#endif //_ApplyMulDPolN_

