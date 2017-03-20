//this program is used to do the following:
//RTPC reconstruction: exe 11 infile.root 
//RTPC calibration: exe 12 
#include "stdlib.h"
#include <iostream>
#include <iomanip>
#include <stdio.h>
#include "math.h"
using namespace std;

#include "MyStyle.h"

#include "MulDFit.h"

void AnaRTPC(const char *infile,const char *outfile="", int readoutpad=2, int ntrack_per_event=1);

void DeltaXXX(const char* infile,const char* vipcut="",
  const char* cutname="Pall", const char* treename="ep");

void PlotPid(const char *infile="nt_ep.root");

void Calib(const char* infile="nt_ep.root");

int main(int argc, char** argv)
{

  if (argc<3) 
  {
    cout<<endl;
    cout<<"Usage: exe <job(10-19)> <infile=nt.root> [outfile=nt_ep.root] [readoutpad=6(1|2|3|4|5|6)] [ntrack_per_event=1]"<<endl;
    cout<<"       exe <job(20-29)> [in_para=in_para.ini] [in_file=in_file.dat]"<<endl;
    cout<<"\nFor 10<=job<=19 or job==99:"<<endl;
    cout<<"       exe <job(10-19)> <infile=nt.root> [outfile=nt_ep.root] [readoutpad=6(1|2|3|4|5|6)] [ntrack_per_event=1]"<<endl;
    cout<<"       readoutpad: 1)4.5x5,  2)2x2,  3)compass 2-D readout, 0.4x0.4 equivalent"<<endl;
    cout<<"       4)compass 2-D readout, 1.0x1.0 equivalent, 5cm drift distance,"<<endl;
    cout<<"       5)compass 2-D readout, 1.0x1.0 equivalent, 4cm drift distance, "<<endl;
    cout<<"       6)Final Configuration: 2.8x4 pad and 4cm drift distance."<<endl;
    cout<<"       ntrack_per_event is how many track in each super event. 1 means no no super evnet \\ \n"
        <<"       therefore no TDC shfited."<<endl;
    cout<<"       job=11  AnaRTPC, will create tree and files for MulDFit"<<endl;
    cout<<"       job=12  Get RTPC Calibration input file RTPC_Calib_Para.inc"<<endl;
    cout<<"       job=13  Plot Delta figures for RTPC"<<endl;
    cout<<"       job=14  Plot dEdX(PID) figures for RTPC\n"<<endl;
    cout<<"       job=18  AnaRTPC, create super event tree, combine ntrack in to one super event by \\ \n"
        <<"               shifting their hits"<<endl;
    cout<<"       job=99  Do all 11, 12, 13 and 14 jobs"<<endl;
    cout<<"example: "<<endl
      <<"       "<<argv[0]<<" 11  ../allz.root  nt_allz.root 6 1"<<endl
      <<"       "<<argv[0]<<" 12  ../allz.root  nt_allz.root "<<endl
      <<"       "<<argv[0]<<" 13  ../allz.root  nt_allz.root "<<endl
      <<"       "<<argv[0]<<" 14  ../allz.root  nt_allz.root "<<endl
      <<"       "<<argv[0]<<" 18  ../allz.root  nt_allz.root "<<endl
      <<"       "<<argv[0]<<" 99  ../allz.root  nt_allz.root "<<endl
      <<endl;
    cout<<"\nFor 20<=job<=29:"<<endl;
    cout<<"       exe <job> [in_para=in_para.ini] [in_file=in_file.dat] <label(A|B|C|D|E)> \\"<<endl
      <<"       [MaxOrderA=5] [MaxOrderB=5] [MaxOrderC=5] [MaxOrderD=5] [MaxOrderE=5] \n"<<endl;
    cout<<"       job=20:  create initial parameter file using the given label and MaxOrders. \\"<<endl;
    cout<<"                The label must be A, B, C, D or E."<<endl;
    cout<<"       job=21:  Run MulDFit to fit for a matrix and check the result"<<endl;
    cout<<"       job=22:  Check the result for a matrix"<<endl;
    cout<<endl;
    cout<<"example: "<<endl
      <<"       "<<argv[0]<<" 20  _in_para_p.ini  in_file_p.dat A 7 0 0 2 6"<<endl
      <<"       "<<argv[0]<<" 20 _in_para_th.ini in_file_th.dat B 0 2 0 2 3"<<endl
      <<"       "<<argv[0]<<" 20 _in_para_ph.ini in_file_ph.dat C 0 2 0 2 3"<<endl
      <<"       "<<argv[0]<<" 20  _in_para_z.ini  in_file_z.dat D 0 2 0 5 3"<<endl
      <<"       "<<argv[0]<<" 21 _in_para_th.ini in_file_th.dat"<<endl
      <<"       "<<argv[0]<<" 22 _in_para_th.ini in_file_th.dat"<<endl
      <<endl;
    exit(-1);
  }

  int job=11, readoutpad=6,ntrack_per_event=1;
  char infile[255], outfile[255];
  char label='A';
  int orderA=5, orderB=5, orderC=5, orderD=5, orderE=5;

  if(argc>1) job=atoi(argv[1]);
  if(argc>2) sprintf(infile,"%s",argv[2]);
  if(job==99 || (job>=10 && job<20))
  {
    if(argc>3) sprintf(outfile,"%s",argv[3]);
    else sprintf(outfile,"%s","nt_ep.root");
    if(argc>4) readoutpad=atoi(argv[4]);
    if(argc>5) ntrack_per_event=atoi(argv[5]);
  }
  if(job>=20 && job<29)
  {
    if(argc>2) system(Form("ln -sf %s in_para.ini",argv[2]));
    if(argc>3) system(Form("ln -sf %s in_file.dat",argv[3]));
    if(argc>4) label=argv[4][0];	
    if(argc>5) orderA=atoi(argv[5]);
    if(argc>6) orderB=atoi(argv[6]);
    if(argc>7) orderC=atoi(argv[7]);
    if(argc>8) orderD=atoi(argv[8]);
    if(argc>9) orderE=atoi(argv[9]);
  }

  SetMyFitStyle();
  if(job==11 || job==99) {AnaRTPC(infile,outfile,readoutpad,ntrack_per_event);}
  if(job==12 || job==99) Calib(outfile);
  if(job==13 || job==99) 
  {
    DeltaXXX(outfile,"P0_p>0.05 && P0_p<0.35","P050to350");
    DeltaXXX(outfile,"P0_p>0.07 && P0_p<0.10 && cos(Theta0_p)<-0.2","P070to100");
    DeltaXXX(outfile,"P0_p>0.17 && P0_p<0.22 && cos(Theta0_p)<-0.2","P170to220");
    DeltaXXX(outfile,"Z0>-10 && Z0<10","Z0000");
    DeltaXXX(outfile,"Z0>90 && Z0<110","Z+100");
    DeltaXXX(outfile,"Z0>-110 && Z0<-90","Z-100");
  }
  if(job==18) {AnaRTPC(infile,outfile,readoutpad,ntrack_per_event);}
  if(job==14 || job==99) PlotPid(outfile);

  if(job==20)
  {
    //prepare input configuration file
    //if(label=='A') MulDFit::Write1stParaFile('A',5,0,0,0,6); //for P
    //if(label=='B') MulDFit::Write1stParaFile('B',0,5,0,3,0); //for theta
    //if(label=='C') MulDFit::Write1stParaFile('C',0,5,1,0,5); //for phi
    //if(label=='D') MulDFit::Write1stParaFile('D',0,3,0,3,1); //for z
    MulDFit::Write1stParaFile(label,orderA,orderB,orderC,orderD,orderE);
    system(Form("cp -f para_template.ini %s",argv[2]));
  } 
  if(job==21)
  {
    //Use mulDFit to fit a matrix, will also check the result
    MulDFit pMulDFit;
    pMulDFit.Process();
  }
  if(job==22)
  {
    //check the result for given matrix and raw data file, do not fit the matrix
    MulDFit pMulDFit;
    pMulDFit.LoadData();
    pMulDFit.CheckResult();
  }

  return 0;

}


//ApplyMulDPolN *pMatrixT = new ApplyMulDPolN("matrix_th.ini");
//ApplyMulDPolN *pMatrixP = new ApplyMulDPolN("matrix_ph.ini");
//ApplyMulDPolN *pMatrixY = new ApplyMulDPolN("matrix_y.ini");
//ApplyMulDPolN *pMatrixD = new ApplyMulDPolN("matrix_d.ini");
//for(;;)
//{
//
//	double x[]={A,B,C,D,E};
//	th_tg = pMatrixT.Eval(x);
//	Ph_tg = pMatrixP.Eval(x);
//	Y_tg = pMatrixY.Eval(x);
//	d_tg = pMatrixD.Eval(x);
//}
