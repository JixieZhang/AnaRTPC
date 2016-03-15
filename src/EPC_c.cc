/*
This function is provided to calculate quasi-elastic and resonance cross section for different kinds of particles
Use EPC model for proton, neutron, pions

Meaning of parameters:
PID: praticle ID, following the PDG definition
=	2212	for p		;	2112	for n	;	211	for pi+	;
-211	for pi-		;	111		for pi0	;	11	for e-	;
22		for photon	;
Z,N: proton and neutron number of the nucleus.
Eb: incoming electron energy in GeV;
theta: scattering angle for outgoing particle in radian;
pf: outgoing particle momentum in GeV/c;
*/
#include <math.h>
namespace EPC
{

	extern "C" {
		void epc_(int *, int *, int *, double *, double *, double *, double *);
	}


	double GetXS(int PID, int Z, int N, double Eb, double theta, double pf)
	{
		double result=0.0;
		double deg=atan(1.0)/45.;
		int A;

		A=Z+N;

		theta=theta/deg;
		Eb=Eb*1000.0;
		pf=pf*1000.0;
	

		if((PID==2212)||(PID==2112)||(PID==211)||(PID==-211)||(PID==111))
			epc_(&PID,&Z,&N,&Eb,&pf,&theta,&result);
		else
			result=0;
		return result;
	}

}


