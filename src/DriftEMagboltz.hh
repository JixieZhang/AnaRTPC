// source code for class  DriftEMagboltz
#ifndef _DriftEMagboltz_
#define _DriftEMagboltz_

#include <iostream>
#include <stdio.h>
#include <math.h>
#include <math.h>

using namespace std;

#define NS_PER_TIC 200    // nano-seconds per tic
#define NUM_PADS 18000    // Number of actual pads 
#define NAL_SAMP 55       // Number of tic samples 


static const double PI=4.0*atan(1.0);
static const double rad2deg=180./PI;
//////////////////////////////////////////////////////////////////////////////////////
typedef struct {
  float phi;
  float z;
  float s;
} cylindrical_t; //this is the raw data type

typedef struct {
  float x;
  float y;
  float z;
} cartician_t; //this array is needed to measure the residuals



//////////////////////////////////////////////////////////////////////////////////////
class DriftEMagboltz
{
public:
  DriftEMagboltz(float R_He2DME=0.9, float pad_w=2.79, float pad_l=4.0, float pad_s=80., 
    float rtpc_l=400.);
  virtual ~DriftEMagboltz();

  //This routine will only be used by simulation
  //input:
  //Initial position x0,y0,z0 (in mm) and deltaE (in KeV)
  //output: x_r,y_r,z_r,chan,adc,tdc (in ns unit)
  //        tdc = (t_s2gem1+t_gem2pad)+tzero, tdc must be a multiple of NS_PER_TIC
  int DriftESim(float x0,float y0,float z0,float deltaE,
    float& x_r,float& y_r,float& z_r,int& chan,int& tdc,int& adc);

  //This routine will be used by simulation
  //One input position will create NofTIC of TDC: TDC_1st, TDC_2nd, TDC_3nd
  //where TDC_2nd=TDC_1st - 1 TDC_3nd=TDC_2st - 1
  //input:
  //Initial position x0,y0,z0 (in mm) and deltaE (in KeV)
  //        NofTIC is one ionization to how many tic, 
  //output: x_r,y_r,z_r,chan,adc,tdc (in tic unit)
  //        tdc = (t_s2gem1+t_gem2pad)+tzero, tdc must be a multiple of NS_PER_TIC
  //return NofTIC
  int DriftESim(float x0,float y0,float z0,float deltaE, int& NofTIC,
    float *x_r,float* y_r,float* z_r,int* chan,int* tdc,int* adc);

  //input: id and tdc in ns unit, this tdc should also include TPC_TZERO
  //output: (x_r,y_r,z_r) from look up table
  void LookupXYZByIDTDC(int chan,int tdc_ns,float& x_r,float& y_r,float& z_r);

  //input:
  //Initial position x0,y0,z0 (in mm) and deltaE (in KeV)
  //output: chan,adc,tdc (in tic unit)
  //        tdc = (t_s2gem1+t_gem2pad)+tzero, tdc must be a multiple of NS_PER_TIC
  int DriftEl2Pad(float x0,float y0,float z0,float dE_kev,int& chan,int& tdc,int& adc);


private:

  void InitElPathCell();
  //input:  channel id and tdc (in ns unit, already include TPC_TZERO)
  //output: reconstruncted ionization location x_r,y_r,z_r (in mm)
  void Reconstruct(int chan,int tdc_ns,float& x_r,float& y_r,float& z_r);


private:
  //the following coming from Magboltz

  //forward
  //return the time that electron drift from ionization location (s0) to GEM1, in ns
  float GetT_s2gem1(float s0_mm,float z0);
  //return the time that electron drift from GEM1 to PAD, in ns
  float GetT_gem2pad(float z0);
  //return the phi difference that electron drift from ionization location (s0) to GEM1, in rad
  float GetdPhi_s2gem1(float s0_mm,float z0);
  //return the phi kick that electron drift from GEM1 to PAD, in rad
  float GetdPhi_gem2pad(float z0);

  //backward
  //return S for given t_s2gem1
  float GetSByT(float t_s2gem1,float z0);
  //return the S_r for given pad_z and pad_phi
  float GetS_r(float z_pad, float phi_pad, float t_s2pad);
  float GetPhi_r(float z_pad, float phi_pad, float t_s2pad);

  //return the reconstructed S_r, Phi_r for given pad_z and pad_phi
  int   GetSPhi_r(float z_pad, float phi_pad, float t_s2pad, float &s_r, float &phi_r);

private:
  //determine the channel id by z(mm) and phi(rad)
  int   GetChanId(float z0,float phi_rad);
  //determine the channel z and phi by chan_id
  int   GetChanZPhi(int chan, float &z, float &phi);


private:
  //cylindrical_t rawCYL[NAL_SAMP][NUM_PADS];
  /*holds the signal data for the first linefit*/
  cartician_t rawXYZ[NAL_SAMP][NUM_PADS];
  /*holds the signal data for the first linefit*/
  float Ratio_He2DME;       //default He2DME Ratio=90:10
  int   TPC_TZERO;          //how long in time the DAQ will read ahead of trigger. in ns unit
  //if set to 1000.0, then 1 adc means 1 ev, if set to 20.0, then 1 adc means 50 ev
  float Kev2ADC;	       //default  coef=50.0; this will match real data better

  //the pad size is 2.79(phi)x4.0(z)  mm^2
  //static const float PAD_W = 2.79;   // PAD width in phi direction
  //static const float PAD_L = 4.0;    // PAD Length in z direction
  //static const float PAD_S = 80.0;   // The radius of the readout pad, 80 mm
  float PAD_W, PAD_L, PAD_S;
  float RTPC_L ;  //= 400.0;   // Length of the RTPC 

};
typedef DriftEMagboltz BonusDriftEMagboltz;
#endif //_DriftEMagboltz_

/////////////////////////////////////////////////////////////
