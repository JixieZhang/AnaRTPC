#ifndef _EPC_H
#define _EPC_H

//This function is provided to calculate quasi-elastic and resonance cross section for different kinds of particles
//Use EPC model for proton, neutron, pions
//
//Meaning of parameters:
//PID: praticle ID, following the PDG definition
//=	2212	for p		;	2112	for n	;	211	for pi+	;
//-211	for pi-		;	111		for pi0	;	11	for e-	;
//22		for photon	;
//Z,N: proton and neutron number of the nucleus.
//Eb: incoming electron energy in GeV;
//theta: scattering angle for outgoing particle in radian;
//pf: outgoing particle momentum in GeV/c;


namespace EPC
{
	double GetXS(int PID, int Z, int N, double Eb, double theta, double pf);
}


#endif
