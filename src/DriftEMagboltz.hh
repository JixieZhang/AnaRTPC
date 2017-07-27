// source code for class  DriftEMagboltz
#ifndef _DriftEMagboltz_
#define _DriftEMagboltz_

#include <iostream>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "ChannelMap.hh"

using namespace std;

#define NS_PER_TIC 120    // nano-seconds per tic
#define NAL_SAMP 90       // Number of tic samples 

//#define NUM_PADS 18000    // Number of actual pads, defined in ChannelMap 

static const double PI=4.0*atan(1.0);
static const double rad2deg=180./PI;
//////////////////////////////////////////////////////////////////////////////////////
typedef struct {
  double phi;
  double z;
  double s;
} cylindrical_t; //this is the raw data type

typedef struct {
  double x;
  double y;
  double z;
} cartician_t; //this array is needed to measure the residuals

//////////////////////////////////////////////////////////////////////////////////////
class DriftEMagboltz
{
public:
  DriftEMagboltz(double R_He2DME=0.9, double pad_w=2.8016, double pad_l=4.1016, 
     double pad_s=80.42, double rtpc_l=400.);
  virtual ~DriftEMagboltz();

  //This routine will only be used by simulation
  //input:
  //Initial position x0,y0,z0 (in mm) and deltaE (in KeV)
  //output: x_r,y_r,z_r,chan,adc,tdc (in ns unit)
  //        tdc = (t_s2gem1+t_gem2pad)+tzero, tdc must be a multiple of NS_PER_TIC
  int DriftESim(double x0,double y0,double z0,double deltaE,
    double& x_r,double& y_r,double& z_r,int& chan,int& tdc,int& adc);

  //This routine will be used by simulation
  //One input position will create NofTIC of TDC: TDC_1st, TDC_2nd, TDC_3nd
  //where TDC_2nd=TDC_1st - 1 TDC_3nd=TDC_2st - 1
  //input:
  //Initial position x0,y0,z0 (in mm) and deltaE (in KeV)
  //        NofTIC is one ionization to how many tic, 
  //output: x_r,y_r,z_r,chan,adc,tdc (in tic unit)
  //        tdc = (t_s2gem1+t_gem2pad)+tzero, tdc must be a multiple of NS_PER_TIC
  //return NofTIC
  int DriftESim(double x0,double y0,double z0,double deltaE, int& NofTIC,
    double *x_r,double* y_r,double* z_r,int* chan,int* tdc,int* adc);

  //input: id and tdc in ns unit, this tdc should also include TPC_TZERO
  //output: (x_r,y_r,z_r) from look up table
  void LookupXYZByIDTDC(int chan,int tdc_ns,double& x_r,double& y_r,double& z_r);

  //input:
  //Initial position x0,y0,z0 (in mm) and deltaE (in KeV)
  //output: chan,adc,tdc (in tic unit)
  //        tdc = (t_s2gem1+t_gem2pad)+tzero, tdc must be a multiple of NS_PER_TIC
  int DriftEl2Pad(double x0,double y0,double z0,double dE_kev,int& chan,int& tdc,int& adc);


private:

  void InitElPathCell();
  //input:  channel id and tdc (in ns unit, already include TPC_TZERO)
  //output: reconstruncted ionization location x_r,y_r,z_r (in mm)
  void Reconstruct(int chan,int tdc_ns,double& x_r,double& y_r,double& z_r);
  
  //boxmuller gauss number generator
  //input: mean m, standard deviation s 
  double Gauss(double m, double s);

private:
  //the following coming from Magboltz

  //forward
  //return the time that electron drift from ionization location (s0) to GEM1, in ns
  double GetT_s2gem1(double s0_mm,double z0);
  //return the time that electron drift from GEM1 to PAD, in ns
  double GetT_gem2pad(double z0);
  //return the phi difference that electron drift from ionization location (s0) to GEM1, in rad
  double GetdPhi_s2gem1(double s0_mm,double z0);
  //return the phi kick that electron drift from GEM1 to PAD, in rad
  double GetdPhi_gem2pad(double z0);

  //backward
  //return S for given t_s2gem1
  double GetSByT(double t_s2gem1,double z0);
  //return the S_r for given pad_z and pad_phi
  double GetS_r(double z_pad, double phi_pad, double t_s2pad);
  double GetPhi_r(double z_pad, double phi_pad, double t_s2pad);

  //return the reconstructed S_r, Phi_r for given pad_z and pad_phi
  int   GetSPhi_r(double z_pad, double phi_pad, double t_s2pad, double &s_r, double &phi_r);

private:
  //determine the channel id by z(mm) and phi(rad)
  int   GetChanId(double z0,double phi_rad);
  //determine the channel z and phi by chan_id
  int   GetChanZPhi(int chan, double &z, double &phi);


private:
  //cylindrical_t rawCYL[NAL_SAMP][NUM_PADS];
  /*holds the signal data for the first linefit*/
  cartician_t rawXYZ[NAL_SAMP][NUM_PADS];
  /*holds the signal data for the first linefit*/
  double Ratio_He2DME;       //default He2DME Ratio=90:10
  int   TPC_TZERO;          //how long in time the DAQ will read ahead of trigger. in ns unit
  //if set to 1000.0, then 1 adc means 1 ev, if set to 20.0, then 1 adc means 50 ev
  double Kev2ADC;	       //default  coef=50.0; this will match real data better

  //the pad size is 2.79(phi)x4.0(z)  mm^2
  //static const double PAD_W = 2.79;   // PAD width in phi direction
  //static const double PAD_L = 4.0;    // PAD Length in z direction
  //static const double PAD_S = 80.0;   // The radius of the readout pad, 80 mm
  double PAD_W, PAD_L, PAD_S;
  double RTPC_L ;  //= 400.0;   // Length of the RTPC 

  ChannelMap *fChanMap;
};
typedef DriftEMagboltz RTPC12DriftEMagboltz;
#endif //_DriftEMagboltz_

/////////////////////////////////////////////////////////////
