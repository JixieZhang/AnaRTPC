//apply RTPC reconstruction correction
#include "stdlib.h"
#include <iostream>
#include "math.h"
using namespace std;

#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
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

#include "BonusHelixFit.hh"
#include "RTPC_Calib_Para.inc"

//boxmuller gauss number generator
//input: mean m, standard deviation s 
double fGaus(double m, double s)	
{				        
	if(s==0.0) return m;

	double x1, x2, w, y1;
	static double y2;
	static int use_last = 0;

	if (use_last)		        /* use value from previous call */
	{
		y1 = y2;
		use_last = 0;
	}
	else
	{
		do {
			x1 = 2.0 * double(rand())/double(RAND_MAX) - 1.0;
			x2 = 2.0 * double(rand())/double(RAND_MAX) - 1.0;
			w = x1 * x1 + x2 * x2;
		} while ( w >= 1.0 );

		w = sqrt( (-2.0 * log( w ) ) / w );
		y1 = x1 * w;
		y2 = x2 * w;
		use_last = 1;
	}

	return( m + y1 * s );
}

void PrintAcc(TH3F *h3,TH3F *h3N, int Nmin=0, int TCS=1)
{
	h3->Print();

	int NBinX = h3->GetNbinsX();
	int NBinY = h3->GetNbinsY();
	int NBinZ = h3->GetNbinsZ(); 
	TAxis* XAxis = h3->GetXaxis();
	TAxis* YAxis = h3->GetYaxis();
	TAxis* ZAxis = h3->GetZaxis(); 
	double x,y,z,acc,err;
	int Ndet=0;
	char str[255];
	FILE *pFile = fopen(Form("%s.txt",h3->GetName()),"w");
	if(!TCS)sprintf(str,"%8s %8s %8s %8s %8s %8s\n","phi/deg","theta/deg","p/GeV","acc(%)","err(abs)","Ndet");
	else sprintf(str,"%8s %8s %8s %8s %8s %8s\n","phi/rad","theta/rad","p/GeV","acc(%)","err(abs)","Ndet");
	//cout<<str;
	fprintf(pFile,"%s",str);

	for(int i=1;i<NBinX;i++){
		x = XAxis->GetBinCenter(i);
		for(int j=1;j<NBinY;j++){
			y = YAxis->GetBinCenter(i);
			for(int k=1;k<NBinZ;k++){
				z = ZAxis->GetBinCenter(i);
				acc = h3->GetBinContent(i,j,k); 
				Ndet = h3N->GetBinContent(i,j,k);
				if(Ndet>=Nmin)
				{
					err = (Ndet>0) ? acc * pow(double(Ndet),-0.5) : 0.0 ;
					sprintf(str,"%8.2f %8.2f %8.3f %8.2f %8.2f %8d\n",x,y,z,acc,err,Ndet);
					//cout<<str;
					fprintf(pFile,"%s",str);
				}
			}
		}
	}
	fclose(pFile);
}


int FitATrack(int StepNum, double StepX[],double StepY[],double StepZ[],
			  double &R, double &A, double &B,
			  double &Phi, double &Theta, double &Z,
			  double RLimit_l=30.1, double RLimit_h=79.9)
{
	int pNpt=0;
	double pData[200][3];
	//helix can fit no more than 200 points
	//fill the vertex point
	pData[pNpt][0]=StepX[0];
	pData[pNpt][1]=StepY[0];
	pData[pNpt++][2]=StepZ[0];
	for (int j=1;j<StepNum;j++)
	{
		if (pNpt>=199) break;
		float tmpRR=sqrt(StepX[j]*StepX[j]+StepY[j]*StepY[j]);
		if (tmpRR<RLimit_l) continue; 
		else if (tmpRR>RLimit_h)  break; 
		else //process data points in the drift region
		{
			pData[pNpt][0]=StepX[j];
			pData[pNpt][1]=StepY[j];
			pData[pNpt++][2]=StepZ[j];
		}
	}

	//***********Do Helix Fit to simulated track*****************//
	if (pNpt<5) return 0;

	const double deg=atan(1.0)/45;
	double Phi_deg,Theta_deg;
	HelixFit(pNpt,pData,R,A,B,Phi_deg,Theta_deg,Z,1);
	//Phi located at [0,2pi], convert it to [-pi,pi]
	if(Phi_deg>180.0) Phi_deg-=360.0;
	Theta=Theta_deg*deg;
	Phi=Phi_deg*deg;

#ifdef BONUS_TREE_DEBUG
	if (Global_Debug_Level>=2)
	{
		printf("***********Do Helix Fit to simulated track %d*****************\n",i);
		printf("Fit Result:R=%.2f, A=%.2f, B=%.2f, Phi=%.2f, Theta=%.2f, Z=%.2f\n",
			tmpR,tmpA,tmpB,tmpPhi,tmpTheta,tmpZ);
	}
#endif

	return 1;
}


double CorrPhi_by_Theta(double Theta_corr, double Phi_rec)
{
	//****************************************
	//Minimizer is Linear                     
	//Chi2                      =      114.301
	//NDf                       =           78
	//p0                        =    -0.286451   +/-   0.0504813   
	//p1                        =     0.702482   +/-   0.316111    
	//p2                        =    -0.473434   +/-   0.786726    
	//p3                        =    -0.176363   +/-   1.01645     
	//p4                        =     0.423359   +/-   0.741797    
	//p5                        =    -0.242468   +/-   0.307939    
	//p6                        =    0.0623113   +/-   0.0677592   
	//p7                        =  -0.00631308   +/-   0.00613088  
	//
	////Pol7 fitted parameters for Phi0_p-Phi_rec vs Theta_rec
	//const double kPara_Pol7_dPhVsTh[8] = {                  
	//        -2.864512E-01, 7.024821E-01, -4.734338E-01, -1.763625E-01, 4.233588E-01, 
	//        -2.424678E-01, 6.231130E-02, -6.313080E-03                               
	//}

	double dPhi=0;
	for(int i=0;i<=7;i++) dPhi += kPara_Pol7_dPhVsTh[i]*pow(Theta_corr,double(i));

	//cout<<"Theta="<<Theta_corr<<"  dPhi="<<dPhi<<"  Phi="<<Phi_rec<<endl;
	double Phi_corr = Phi_rec + dPhi;

	return Phi_corr;
}



double CorrPhi(double P_corr, double Phi_rec)
{                  
	//Minimizer is Linear                     
	//Chi2                      =      64.6663
	//NDf                       =           77
	//p0                        =      7.93566   +/-   1.13459     
	//p1                        =     -350.845   +/-   57.8046     
	//p2                        =      6574.97   +/-   1227.71     
	//p3                        =     -67679.7   +/-   14102.5     
	//p4                        =       413046   +/-   94728.8     
	//p5                        = -1.49426e+06   +/-   372572      
	//p6                        =  2.96658e+06   +/-   795522      
	//p7                        = -2.49303e+06   +/-   712379      
	//
	////Pol7 fitted parameters for Phi0_p-Phi_rec vs P0_p
	//const double kPara_Pol7_dPhVsP[8] = {              
	//        7.935655E+00, -3.508450E+02, 6.574966E+03, -6.767971E+04, 4.130460E+05, 
	//        -1.494263E+06, 2.966578E+06, -2.493027E+06                              
	//};                                   
	double dPhi=0;
	for(int i=0;i<=7;i++) dPhi += kPara_Pol7_dPhVsP[i]*pow(P_corr,double(i));

	//cout<<"Theta="<<Theta_corr<<"  dPhi="<<dPhi<<"  Phi="<<Phi_rec<<endl;
	double Phi_corr = Phi_rec + dPhi;

	return Phi_corr;
}


double CorrTheta(double Theta_rec)
{
	//****************************************
	//Minimizer is Linear                     
	//Chi2                      =      102.661
	//NDf                       =           78
	//p0                        =    0.0103128   +/-   0.0312895   
	//p1                        =     0.164119   +/-   0.230239    
	//p2                        =    -0.649925   +/-   0.645826    
	//p3                        =      1.01193   +/-   0.908762    
	//p4                        =    -0.799845   +/-   0.702929    
	//p5                        =     0.337756   +/-   0.302911    
	//p6                        =   -0.0726412   +/-   0.0681149   
	//p7                        =   0.00624624   +/-   0.00622485  

	//Pol7 fitted parameters for Theta0_p-Theta_rec vs Theta_rec
	//const double kPara_Pol7_dThVsTh[8] = {                      
	//	1.031278E-02, 1.641190E-01, -6.499250E-01, 1.011926E+00, -7.998452E-01, 
	//	3.377556E-01, -7.264124E-02, 6.246237E-03                               
	//};	
	double dTheta=0;
	for(int i=0;i<=7;i++) dTheta += kPara_Pol7_dThVsTh[i]*pow(Theta_rec,double(i));

	double Theta_corr = Theta_rec + dTheta;

	//cout<<"Theta_rec="<<Theta_rec<<"  dSinTh="<<dSinTh<<"  Theta_corr="<<Theta_corr<<endl;
	return Theta_corr;
}

double GetPtanOverBzByR(double R)
{
	//****************************************
	//Minimizer is Linear                     
	//Chi2                      =      41.9928
	//NDf                       =           32
	//p0                        =      3.91575   +/-   1.61988     
	//p1                        =    -0.867099   +/-   0.347475    
	//p2                        =     0.081348   +/-   0.0315288   
	//p3                        =  -0.00417961   +/-   0.00156876  
	//p4                        =  0.000126991   +/-   4.62327e-05 
	//p5                        = -2.28036e-06   +/-   8.0721e-07  
	//p6                        =  2.23988e-08   +/-   7.73353e-09 
	//p7                        = -9.28288e-11   +/-   3.13752e-11 
	//
	////Pol7 fitted parameters for P0_p*sin(Theta0_p)/abs(MeanBz) vs R_sim
	//const double kPara_Pol7_PperpOverBzVsR_1[8] = {                     
	//        3.915747E+00, -8.670992E-01, 8.134804E-02, -4.179609E-03, 1.269908E-04, 
	//        -2.280364E-06, 2.239879E-08, -9.282882E-11                              
	//};             
	//****************************************
	//Minimizer is Linear                     
	//Chi2                      =      66.4879
	//NDf                       =           50
	//p0                        =    0.0163488   +/-   0.00211719  
	//p1                        = -0.000164565   +/-   8.7342e-05
	//p2                        =  4.57719e-06   +/-   1.17494e-06
	//p3                        = -1.54772e-08   +/-   5.15649e-09
	//
	////Pol3 fitted parameters for P0_p*sin(Theta0_p)/abs(MeanBz) vs R_sim
	//const double kPara_Pol3_PperpOverBzVsR_2[4] = {
	//        1.634879E-02, -1.645654E-04, 4.577192E-06, -1.547722E-08
	//};
	//
	//
	//****************************************
	//Minimizer is Linear
	//Chi2                      =      68.2429
	//NDf                       =           55
	//p0                        =   0.00122746   +/-   8.02035e-05
	//p1                        =  0.000289282   +/-   4.87419e-07
	//
	////Pol1 fitted parameters for P0_p*sin(Theta0_p)/abs(MeanBz) vs R_sim
	//const double kPara_Pol1_PperpOverBzVsR_3[2] = {
	//        1.227465E-03, 2.892817E-04
	//};



	double PtranOverBz=0;
	if(R<20) R=20;
	if(R<46)
	{
		for(int i=0;i<=7;i++) 
			PtranOverBz += kPara_Pol7_PperpOverBzVsR_1[i]*pow(R,double(i));
	}
	else if(R<100)
	{
		for(int i=0;i<=3;i++) 
			PtranOverBz += kPara_Pol3_PperpOverBzVsR_2[i]*pow(R,double(i));
	}
	else
	{
		for(int i=0;i<=1;i++) 
			PtranOverBz += kPara_Pol1_PperpOverBzVsR_3[i]*pow(R,double(i));
	}
	return PtranOverBz;
}



void RTPC_Recon(double MeanBz, double R_rec, double Theta_rec, double Phi_rec,
				double &P_corr, double &Theta_corr, double &Phi_corr)
{
	double PtranOverBz = GetPtanOverBzByR(R_rec);
	Theta_corr = CorrTheta(Theta_rec);
	//Theta_corr = Theta_rec;

	P_corr = -PtranOverBz*MeanBz/sin(Theta_corr);

	//Phi_corr = CorrPhi(P_corr,Phi_rec);
	//Phi_corr = CorrPhi_by_Theta(Theta_corr,Phi_rec);
	Phi_corr = Phi_rec;
}
