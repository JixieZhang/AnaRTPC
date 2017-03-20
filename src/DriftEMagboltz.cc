// ********************************************************************
//
// $Id: DriftEMagboltz.cc,v 1.1, 2017/3/19 BONUS12 Exp $
// --------------------------------------------------------------
//
/*
This class is to simulate the drift path of ionizaton electron in the drift region
of the RTPC12 detector.  The TDC values have been digitalized in ns, it is a
multiple of NS_PER_TIC(200ns). The tdc include tpc_tzero, which is the extra 
time that the DAQ will read ahead of trigger, and t_gem2pad, the time that the 
ionization electron takes to drift from gem1 to pad. That says,
tdc = t_s2gem1 + t_gem2pad + tzero.

Magboltz simulation provide the following functions:
(note that z is one of the self variable because B field is not uniform)
1. t_s2gem1(s,z)     ==>drift time from the ionization location to gem1
2. t_gem2pad(z)      ==>the drift time from the first gem to pads
3. dPhi_s2gem1(s,z)  ==>phi deflection from cathode to R0
4. dPhi_gem2pad(z)   ==>phi deflection from the first gem to pads

5. GetS(z_pad,phi_pad,t_s2pad)   ==>inverse function of 1, to calculate the initial S
6. GetPhi(z_pad,phi_pad,t_s2pad) ==>inverse function of 3, to calculate the initial Phi

1,2,3,4 are used in digitization.
5,6 are used in reconstruction.
*/
////////////////////////////////////////////////////////////////////////////////////


////////////////////////
//#define DRIFTESIM_DEBUG 1
////////////////////////
#include <stdio.h>
#include "DriftEMagboltz.hh"


//////////////////////////////////////////////////////////////////////////////////////
DriftEMagboltz::DriftEMagboltz(float R_He2DME, float pad_w, float pad_l, float pad_s,
  float rtpc_l) : Ratio_He2DME(R_He2DME),PAD_W(pad_w),PAD_L(pad_l),PAD_S(pad_s),RTPC_L(rtpc_l)
{
  Kev2ADC=50.0;	 //set to 50.0, so 1 adc unit = 20 ev

  //TPC_TZERO is used to indicate how long in time the DAQ will read ahead of trigger. in ns unit
  //It is usually set when develop DAQ software, we do not know it now. 
  //I set it to 1600ns here.  Remember to update it when it is set
  TPC_TZERO = 1000; //1600;   
  //intialize the path cell for reconstruction
  InitElPathCell();
}
//////////////////////////////////////////////////////////////////////////////////////
DriftEMagboltz::~DriftEMagboltz()
{
  ;
}

//get drift time by s
float DriftEMagboltz::GetT_s2gem1(float s0_mm,float z0)
{
  z0+=0;
  float s0=s0_mm/10.;
  float t_us=0.0; 
  float a=-0.1642, b=-0.0947,c=8.8001;
  t_us = a*s0*s0+b*s0+c;
  return t_us*1000.;
}

//get time offset from 1st gem to pad
float DriftEMagboltz::GetT_gem2pad(float z0)
{
  //Assume A) gem1 to pad is 10mm,  B) the time to travel from gem1 to pad
  //equal to that it drifts 10cm to gem1
  return GetT_s2gem1(PAD_S-20.0,z0);
}

//get dphi by s
float DriftEMagboltz::GetdPhi_s2gem1(float s0_mm,float z0)
{
  z0+=0;
  float s0=s0_mm/10.;
  float a=0.0287, b=-0.5334, c=2.3475;
  float dphi=a*s0*s0+b*s0+c;
  return dphi;
}

//get dphi offset from 1st gem to pad
float DriftEMagboltz::GetdPhi_gem2pad(float z0)
{
  //Assume A) gem1 to pad is 10mm,  B) the phi deflection from gem1 to pad
  //equal to that to it drifts 10cm to gem1
  return GetdPhi_s2gem1(PAD_S-20.0,z0);
}

//reconstruct s by drifting time
//return negative s if error
float DriftEMagboltz::GetSByT(float t_s2gem1,float z0)
{
  float s_r=-10.0;
  z0+=0;
  float t_us=t_s2gem1/1000.;
  float a=-0.0335, b=-0.3208,c=6.9481;
  s_r = a*t_us*t_us+b*t_us+c;   //in cm

#ifdef DRIFTESIM_DEBUG
  if( DRIFTESIM_DEBUG>=4 )
  {
    printf("GetSByT(t_s2gem1=%.0f) ==> s_r=%.1fmm\n",t_s2gem1,s_r*10);
  }
#endif
  return s_r*10;   //turn s_r from cm to mm
}

//////////////////////////////////////////////////////////////////////////////////////
float DriftEMagboltz::GetS_r(float z_pad, float phi_pad, float t_s2pad)
{
  //input z_pad in mm, phi_pad in rad, t_s2pad in ns
  //output: s0_r, the radius in mm
  float s_r=-10.0;
  phi_pad += 0.0;  //to avoid warning

  //determine drift time from ionization location to gem1
  float t_s2gem1 = t_s2pad - GetT_gem2pad(z_pad);;

  //This part is the reverse function, in mm
  s_r = GetSByT(t_s2gem1,z_pad);  

#ifdef DRIFTESIM_DEBUG
  if( DRIFTESIM_DEBUG>=3 )
  {
    printf("GetS_r(z_pad=%.1fmm,phi_pad=%.1fdeg,t_s2pad=%.0f) ==> s_r=%.1fmm\n",
      z_pad, phi_pad*rad2deg, t_s2pad, s_r);
  }
#endif
  return s_r;  
}


//////////////////////////////////////////////////////////////////////////////////////
float DriftEMagboltz::GetPhi_r(float z_pad, float phi_pad, float t_s2pad)
{
  //input: z_pad in mm, phi_pad in rad, t_s2pad in ns
  //output: phi0_r in rad

  //since we do not have the function of dphi(t), we need to get S_r first
  //then use s the get total dphi
  float s_r = GetS_r(z_pad, phi_pad, t_s2pad);

  //calculate dphi
  float delta_phi = GetdPhi_s2gem1(s_r, z_pad) + GetdPhi_gem2pad(z_pad);

  float phi_r = phi_pad + delta_phi;
  if( phi_r < 0 ) phi_r += 2*PI;
  if( phi_r > 2*PI ) phi_r -= 2*PI;

#ifdef DRIFTESIM_DEBUG
  if( DRIFTESIM_DEBUG>=3 )
  {
    printf("GetPhi_r(z_pad=%.1fmm,phi_pad=%.1fdeg,t_s2pad=%.0f) ==> phi_r=%.1fdeg\n",
      z_pad, phi_pad*rad2deg, t_s2pad, phi_r*rad2deg);
  }
#endif

  return phi_r;
}

//////////////////////////////////////////////////////////////////////////////////////
int DriftEMagboltz::GetSPhi_r(float z_pad, float phi_pad, float t_s2pad, float &s_r, float &phi_r)
{
  //input: z_pad in mm, phi_pad in rad, t_s2pad in ns
  //output: phi0_r in rad

  //calculate S_r
  s_r = GetS_r(z_pad, phi_pad, t_s2pad);

  //calculate total dphi using S_r
  float dphi_gem2pad = GetdPhi_gem2pad(z_pad);
  float dphi_s2gem1=GetdPhi_s2gem1(s_r, z_pad);
  float delta_phi = dphi_s2gem1 + dphi_gem2pad;

  phi_r = phi_pad + delta_phi;
  if( phi_r < 0 ) phi_r += 2*PI;
  if( phi_r > 2*PI )	phi_r -= 2*PI;

#ifdef DRIFTESIM_DEBUG

  if( DRIFTESIM_DEBUG>=3 )
  {
    printf("GetSPhi_r(): dphi_s2gem1=%.1fdeg, dphi_gem2pad=%.1fdeg, delta_phi=%.1fdeg\n",
      dphi_s2gem1*rad2deg, dphi_gem2pad*rad2deg, delta_phi*rad2deg);
  }
  if( DRIFTESIM_DEBUG>=2 )
  {
    printf("GetSPhi_r(z_pad=%.1fmm,phi_pad=%.1fdeg,t_s2pad=%.0f) ==> s_r=%.1fmm,  phi_r=%.1fdeg\n",
      z_pad, phi_pad*rad2deg, t_s2pad, s_r, phi_r*rad2deg);
  }
#endif

  return 0;
}

//////////////////////////////////////////////////////////////////////////////////////
//input Phi angle in rad, from 0 to 2pi
int DriftEMagboltz::GetChanId(float z0, float phi_rad)
{
  //need to fill this routine when the readout pad pattern is ready
  //row shifting in z : 
  //row 0|1|2|3: 0|1|2|3 mm
  float phi_per_pad = PAD_W/PAD_S;
  int row = int(phi_rad/phi_per_pad);
  float z_shift = row%4;
  float z_start=z_shift-RTPC_L/2;
  float z_end=z_shift+RTPC_L/2;
  if(z0<z_start+0.01) return -10;
  else if(z0>z_end-0.01) return -9;
  int col = int((z0-z_shift+RTPC_L/2)/PAD_L);
  int Num_of_Col = int(ceil(RTPC_L/PAD_L));
  return row*Num_of_Col+col;
}

int DriftEMagboltz::GetChanZPhi(int chan, float &z, float &phi)
{
  //need to fill this routine when the readout pad pattern is ready
  float phi_per_pad = PAD_W/PAD_S;
  int Num_of_Col = int(ceil(RTPC_L/PAD_L));
  int row=chan/Num_of_Col;
  int col=chan%Num_of_Col; 
  float z_shift = row%4;
  z=(col+0.5)*PAD_L-RTPC_L/2+z_shift ;
  phi=(row+0.5)*phi_per_pad;
  return 0;
}

//////////////////////////////////////////////////////////////////////////////////////
int DriftEMagboltz::DriftEl2Pad(float x0,float y0,float z0,float dE_kev,
  int& chan,int& tdc,int& adc)
{
  //input:
  //Initial position x0,y0,z0 (in mm) and dE_kev (in KeV)
  //output: chan,adc,tdc (in nano second, not tic unit)
  //        tdc = int((t_s2gem1+t_gem2pad+tzero)/NS_PER_TIC) * NS_PER_TIC 

  //reset the values
  chan=-10;  tdc=-1;  adc=-1;

  float r0,phi0_rad;
  //convert (x0,y0,z0) into (r0,phi0,z0)
  r0=sqrt(x0*x0+y0*y0);  //in mm
  phi0_rad=atan2(y0,x0); //return (-Pi, + Pi)
  if( phi0_rad<0.)  phi0_rad+=2.0*PI;

#ifdef DRIFTESIM_DEBUG
  if( DRIFTESIM_DEBUG>=2 )
  {
    printf("\nDriftEl2Pad(x0=%.2f,y0=%.2f,z0=%.2f)",x0,y0,z0);
    printf(" ==> (r0,phi0_deg,z0)=(%.2f,%.2f,%.2f)\n",r0, phi0_rad*rad2deg,z0);
  }
#endif

  //determine drift_time
  float t_gem2pad = GetT_gem2pad(z0);
  float t_s2gem1 = GetT_s2gem1(r0, z0);
  float t_s2pad = t_s2gem1 + t_gem2pad;

#ifdef DRIFTESIM_DEBUG
  if( DRIFTESIM_DEBUG>=3 )
  {
    printf("\nDriftEl2Pad(): t_gem2pad=%.0f  t_s2gem1=%.0f  t_s2pad=%.0f\n",
      t_gem2pad,t_s2gem1,t_s2pad);
  }
#endif
  //////////////////////////////////////////////////////////
  //determine delta_phi, 
  float dphi_offset=GetdPhi_gem2pad(z0);
  float dphi_s2gem1=GetdPhi_s2gem1(r0,z0);
  float delta_phi = dphi_s2gem1 + dphi_offset;

#ifdef DRIFTESIM_DEBUG
  if( DRIFTESIM_DEBUG>=3 )
  {
    printf("DriftEl2Pad(): dphi_s2gem1=%.1fdeg, dphi_offset=%.1fdeg, delta_phi=%.1fdeg\n",
      dphi_s2gem1*rad2deg, dphi_offset*rad2deg, delta_phi*rad2deg);
  }
#endif
  //check the output//////////////////////////////////////////////////////////

  float phi_rad=float(phi0_rad-delta_phi);   //phi at pad pcb board
  //this is not the center of the pad yet!!!
  //convert the phi_deg into range [0,2*PI)
  if( phi_rad<0. )  phi_rad+=2.0*PI;

  //output////////////////////////////////////////////////////////////////
  //let z_pad=z0, ignore the motoin on z irection
  chan=GetChanId(z0,phi_rad);
  //if(chan<0) return chan;  //let it finish, do not stop here
  //if chan=-1 it means this channel is unreconstrutable
  tdc=int((t_s2pad+TPC_TZERO)/NS_PER_TIC)*NS_PER_TIC;
  adc=int(Kev2ADC*dE_kev);

#ifdef DRIFTESIM_DEBUG
  if( DRIFTESIM_DEBUG>=2 )
  {
    float z_pad_c,phi_pad_c;
    GetChanZPhi(chan,z_pad_c,phi_pad_c);
    if(phi_pad_c<0) phi_pad_c+=2.0*PI;
    printf("DriftEl2Pad(): output: (r,phi_deg,z)=(80,%.1f,%.1f) phi_pad=%.1f, z_pad=%.1f, chan=%5d, tdc=%5d, adc=%5d\n",
      phi_rad*rad2deg,z0,phi_pad_c*rad2deg,z_pad_c,chan,tdc,adc);
  }
#endif
  return chan;
}

//////////////////////////////////////////////////////////////////////////////////////
void DriftEMagboltz::Reconstruct(int chan,int tdc_ns,float& x_r,float& y_r,float& z_r)
{
  //input:  channel id and tdc (in ns unit, already include TPC_TZERO)
  //output: reconstruncted ionization location x_r,y_r,z_r (in mm)

  if( chan<0 || chan>=NUM_PADS )
  {
    printf("**Error! Invalid Channel, ID=%d",chan);
    return;
  }

  float z_mm=-210.0,phi_rad=-10.0;
  GetChanZPhi(chan,z_mm,phi_rad);

  //reset the values
  x_r=0.;y_r=0.;z_r=0.;
  float t_s2pad=tdc_ns-TPC_TZERO+0.5*NS_PER_TIC;
  float s_r_mm, phi_r_rad;
  GetSPhi_r(z_mm, phi_rad, t_s2pad, s_r_mm, phi_r_rad);

  if(s_r_mm>0.0) 
  {
    x_r=s_r_mm*cos(phi_r_rad);
    y_r=s_r_mm*sin(phi_r_rad);
    z_r=z_mm;
  }
  else
  {
    s_r_mm=0.0;
    return;
  }
#ifdef DRIFTESIM_DEBUG
  if( DRIFTESIM_DEBUG>=2 )
  {
    printf("Reconstruct(chan=%4d, tdc=%4d) ==> tic=%d, phi_pad_deg=%.2f, z_pad=%.2f\n",
      chan,tdc_ns,int(tdc_ns/NS_PER_TIC),phi_rad*rad2deg,z_mm);
    printf("==>Output(x_r,y_r,z_r)=(%.2f,%.2f,%.2f) ",x_r,y_r,z_r);
    printf(" or (s_r,phi_r_deg,z_r)=(%.2f,%.2f,%.2f)\n\n",
      s_r_mm, phi_r_rad*rad2deg,z_r);
  }    
#endif
  return;

}

////////////////////////////////////////////////////////////////////////////
//this routine will reconstruct each (id,tdc) before the program starts then 
//fill the result into a buffer. During analysis, the reconstruction will 
//go to this buffer to take the result.
//this will speed up the reconstruction 

void DriftEMagboltz::InitElPathCell()
{
  int chan, tdc, tic;
  float x_r, y_r, z_r, r_r;
  printf("DriftEMagboltz::InitElPathCell(): Initializing reconstruction path cells......\n");

  for( chan = 0; chan<NUM_PADS; chan++ )
  {
    //Kuhn asked to include ellary and later hits to study background
    for( tic=0;tic<NAL_SAMP;tic++ )    
    {
      tdc = tic*NS_PER_TIC;
      Reconstruct(chan,tdc, x_r, y_r, z_r);
      r_r=sqrt( x_r*x_r+y_r*y_r);
      if( r_r>=PAD_S || r_r<20. ) //set default values
      {
        rawXYZ[tic][chan].x = 0.;
        rawXYZ[tic][chan].y = 0.;
        rawXYZ[tic][chan].z = 0.;
        //to save time, jump out whem r_r<25.0 mm
        //if( r_r<25. ) break;
      }
      else
      {
        rawXYZ[tic][chan].x = x_r;
        rawXYZ[tic][chan].y = y_r;
        rawXYZ[tic][chan].z = z_r;
      }
    }//end loop over tic time bins
  }//end loop over pads
  printf("DriftEMagboltz::InitElPathCell(): Initializing reconstruction path cells done!\n");

#if defined DRIFTESIM_DEBUG 
  if ( DRIFTESIM_DEBUG>=1 )
  {
    FILE *pFile= fopen ("DriftPath_Jixie.txt" , "w");
    cout<<"\nJixie's Reconstruction map"<<endl;

    for( chan = 0; chan<=800; chan+=100 )
    {
      fprintf(pFile,"\n  ID  TIC    x(mm)    y(mm)    z(mm)    r(mm) phi(deg)\n");
      for( tic=0;tic<NAL_SAMP;tic++ )
      {
        float r_r=sqrt(rawXYZ[tic][chan].x*rawXYZ[tic][chan].x+
          rawXYZ[tic][chan].y*rawXYZ[tic][chan].y);
        float phi_deg=atan2(rawXYZ[tic][chan].y,rawXYZ[tic][chan].x)*rad2deg;
        if(phi_deg<0.) phi_deg+=360.;
        fprintf(pFile,"%4d %4d %8.2f %8.2f %8.2f %8.2f %8.2f\n",
          chan,tic,
          rawXYZ[tic][chan].x,
          rawXYZ[tic][chan].y,
          rawXYZ[tic][chan].z,
          r_r,phi_deg);
      }
    }
    fclose(pFile);
  }
#endif

}//end InitElPathCell

//////////////////////////////////////////////////////////////////////////////////////
//input: id and tdc_ns in ns unit
//output: (x_r,y_r,z_r) from look up table
void DriftEMagboltz::LookupXYZByIDTDC(int chan,int tdc_ns,
  float& x_r,float& y_r,float& z_r)
{
  x_r=y_r=z_r=0.0;
  int tic = int((tdc_ns-TPC_TZERO)/NS_PER_TIC);
  if(tic<0 || tic>=NAL_SAMP) return;
  //get (x_r, y_r, z_r) from vector rawXYZ
  x_r=this->rawXYZ[tic][chan].x;
  y_r=this->rawXYZ[tic][chan].y;
  z_r=this->rawXYZ[tic][chan].z;  
}

//////////////////////////////////////////////////////////////////////////////////////
//input: (x0,y0,z0) in mm and deltaE in KeV
//output: (x_r,y_r,z_r) and tdc in tic unit
//        tdc = int((t_s2gem1+t_gem2pad+tzero)/NS_PER_TIC) * NS_PER_TIC 
int DriftEMagboltz::DriftESim(float x0,float y0,float z0,float deltaE,
  float& x_r,float& y_r,float& z_r, int& chan,int& tdc,int& adc)
{
  //reset the values
  x_r=0.; y_r=0.; z_r=0.;  chan=-10;  tdc=-1;  adc=-1;

  float s0=sqrt(x0*x0+y0*y0);	//x0,y0,z0 in mm unit

  if( s0>PAD_S-10.0 || fabs(z0) > RTPC_L/2.+PAD_L )
  {
#ifdef DRIFTESIM_DEBUG
    if( DRIFTESIM_DEBUG>=1 )
    {
      printf("Magboltz:This inonized electron is out of Drift Volumn!!! \
             s0=%.1f > %.1f(first gem) or |z0|=%.1f > %.1f\n",
             s0,PAD_S-10.0,z0,RTPC_L/2.+PAD_L);
    }
#endif
    chan=-2;
    return -2;
  }

#ifdef DRIFTESIM_DEBUG
  if( DRIFTESIM_DEBUG>=1 )
  {
    float phi0_rad=atan2(y0,x0);
    if( phi0_rad<0. ) phi0_rad+=2.0*PI;
    printf("\nDriftESim(x0=%.2f,y0=%.2f,z0=%.2f)",x0,y0,z0);
    printf(" ==> (r0,phi0_deg,z0)=(%.2f,%.2f,%.2f)\n",s0, phi0_rad*rad2deg,z0);
  }
#endif

  //please note that the tdc is on tic unit
  int status=DriftEl2Pad( x0, y0, z0, deltaE, chan, tdc, adc);
  if( status<0 ) return status; 

  // dead region or unreconstructable
  int tic = int((tdc-TPC_TZERO)/NS_PER_TIC);
  if( chan<0 || chan>=NUM_PADS || tic>=NAL_SAMP || tic<0 ) return -1;

  //get (x_r, y_r, z_r) from vector rawXYZ
  x_r=this->rawXYZ[tic][chan].x;
  y_r=this->rawXYZ[tic][chan].y;
  z_r=this->rawXYZ[tic][chan].z;
#if defined DRIFTESIM_DEBUG 
  if ( DRIFTESIM_DEBUG>=1 )  
  {
    float phi0_rad, r_r, phi_r_rad, ds, dphi;
    phi0_rad=atan2(y0,x0);
    if( phi0_rad<0. ) phi0_rad+=2.0*PI;

    r_r=sqrt(x_r*x_r+y_r*y_r);
    ds=sqrt(x0*x0+y0*y0)-r_r;
    phi_r_rad=atan2(y_r,x_r);
    if( phi_r_rad<0. ) phi_r_rad+=2.0*PI;
    dphi=(phi0_rad-phi_r_rad)*rad2deg;
    printf("DriftESim(Jixie)==>Final output(x_r,y_r,z_r)=(%.2f,%.2f,%.2f)\n \
           -->(r,phi_deg,z)=(%.2f,%.2f,%.2f); dS=%.2f, dPhi=%.2f\n",
           x_r,y_r,z_r,r_r,phi_r_rad*rad2deg,z_r,ds,dphi);

  }
#endif

  return 0;
}


//////////////////////////////////////////////////////////////////////////////////////
//input: (x0,y0,z0) in mm and deltaE in KeV
//One input position will create NofTIC of TIC: TIC[NofTIC]
//if NofTIC==2, TIC[1]=TIC[0]+1  
//and  TIC[0] = tic from the original reconstruction code
//if NofTIC==3, TIC[2]=TIC[1]-1 and TIC[0]=TIC[1]+1 
//and  TIC[1] = tic from the original reconstruction code
//output: (x_r,y_r,z_r) and tdc in ns unit
//        tdc = int((t_s2gem1+t_gem2pad+tzero)/NS_PER_TIC) * NS_PER_TIC 
int DriftEMagboltz::DriftESim(float x0,float y0,float z0,float deltaE,int &NofTIC,
  float *x_r,float* y_r,float* z_r,int* chan,int* tdc,int* adc)
{
  //Suggested number of output tics for each ionization
  //the wave form is about 380 ns in width
  NofTIC = 380.0/NS_PER_TIC + 1;
  if(NofTIC>3) NofTIC=3;

  //reset the values 
  for(int t=0;t<NofTIC;t++) {
    chan[t]=-10;  
    tdc[t]=adc[t]=-1;
    x_r[t]=y_r[t]=z_r[t]=0.; 
  }

  float s0=sqrt(x0*x0+y0*y0);	//x0,y0,z0 in mm unit

  if( s0>PAD_S-10.0 || fabs(z0) > RTPC_L/2.+PAD_L )
  {
#ifdef DRIFTESIM_DEBUG
    if( DRIFTESIM_DEBUG>=1 )
    {
      printf("Magboltz:This inonized electron is out of Drift Volumn!!! \
             s0=%.1f > %.1f(first gem) or |z0|=%.1f > %.1f\n",
             s0,PAD_S-10.0,z0,RTPC_L/2.+PAD_L);
    }
#endif

    for(int t=0;t<NofTIC;t++) chan[t]=-2;
    return -2;
  }

#ifdef DRIFTESIM_DEBUG
  if( DRIFTESIM_DEBUG>=1 )
  {
    float phi0_rad=atan2(y0,x0);
    if( phi0_rad<0. ) phi0_rad+=2.0*PI;
    printf("\nDriftESim(x0=%.2f,y0=%.2f,z0=%.2f)",x0,y0,z0);
    printf(" ==> (r0,phi0_deg,z0)=(%.2f,%.2f,%.2f)\n",s0, phi0_rad*rad2deg,z0);
  }
#endif

  //please note that the tdc is in ns, nottic unit
  int status=DriftEl2Pad( x0, y0, z0, deltaE, chan[0], tdc[0], adc[0]);
  if ( status<0 ) return status; 

  // dead region or unreconstructable
  int tic = int((tdc[0]-TPC_TZERO)/NS_PER_TIC);
  if ( chan[0]<0 || chan[0]>=NUM_PADS || tic>=NAL_SAMP || tic<0 ) return -1;

  //Now determine how many reconstructed hit it might have
  //and determine how to distribute the ADC between them
  if (NofTIC>tic+1) NofTIC = tic+1;   //take care of tic==0 and tic==1 cases
  if (NofTIC>1 && adc[0]<3) NofTIC = 1;
  if (NofTIC>2 && adc[0]<7) NofTIC = 2;

  if (NofTIC == 2) 
  {//now split the energy by 2:1
    adc[1] = int(adc[0]*1.0/3.0);
    adc[0] -= adc[1];
    tdc[1] = tdc[0];
    tdc[0] += NS_PER_TIC;
    chan[1] = chan[0];
  }
  else if (NofTIC == 3) 
  {//now split the energy by 2:4:1
    adc[2] = int(adc[0]*1.0/7.0);
    adc[1] = int(adc[1]*4.0/7.0);
    adc[0] -= (adc[1] + adc[2]);
    tdc[2] = tdc[0]-NS_PER_TIC;
    tdc[1] = tdc[0];
    tdc[0] += NS_PER_TIC;
    chan[2]=chan[1]=chan[0];
  }

  //get (x_r, y_r, z_r) from vector rawXYZ
  for(int t=0;t<NofTIC;t++) 
  {
    tic = tdc[t]/NS_PER_TIC;
    x_r[t]=this->rawXYZ[tic][chan[t]].x;
    y_r[t]=this->rawXYZ[tic][chan[t]].y;
    z_r[t]=this->rawXYZ[tic][chan[t]].z;
  }


#if defined DRIFTESIM_DEBUG 
  if ( DRIFTESIM_DEBUG>=1 )  
  {
    float r_r, phi_r_rad, ds, dphi;
    float phi0_rad = atan2(y0,x0);
    if( phi0_rad<0. ) phi0_rad+=2.0*PI;
    printf("DriftESim(Jixie)==>Final output: %d points:\n",NofTIC);
    for(int t=0;t<NofTIC;t++) 
    {
      r_r=sqrt(x_r[t]*x_r[t]+y_r[t]*y_r[t]);
      ds=sqrt(x0*x0+y0*y0)-r_r;
      phi_r_rad=atan2(y_r[t],x_r[t]);
      if( phi_r_rad<0. ) phi_r_rad+=2.0*PI;
      dphi=(phi0_rad-phi_r_rad)*rad2deg;
      printf("\t %d: tdc=%d (x_r,y_r,z_r)=(%.2f,%.2f,%.2f) -->(r,phi_deg,z)=(%.2f,%.2f,%.2f); dS=%.2f, dPhi=%.2f\n",
        t,tdc[t],x_r[t],y_r[t],z_r[t],r_r,phi_r_rad*rad2deg,z_r[t],ds,dphi);
    }
  }
#endif

  return NofTIC;
}
