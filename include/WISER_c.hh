#ifndef _WISER_H
#define _WISER_H

namespace WISER 
{
	//input:
	//       Z,N: proton and neutron number of the nucleus.
	//       PART: praticle type, 1(pi+), 2(pi-), 3(k+), 4(k-),5(pr);
	//Ei, Ef: incoming electron energy and outgoing pion momentum in GeV;
	//Theta: scattering angle for outgoing particle in radian;
	//RadLen: material thickness in unit of radiation length;
	double WISER(int PART, int Z, int N, double Ei, double Pf, double theta, double radlen);


	//input
	//		Z,N: proton and neutron number of the nucleus.
	//		Pid: praticle ID, following the PDG definition
	//       	2212	for p;	+/-211	for pi+/-;  +/-321	for k+/-;
	//Ei, Ef: incoming electron energy and outgoing pion momentum in GeV;
	//		Theta: scattering angle for outgoing particle in radian;
	//		RadLen: material thickness in unit of radiation length;
	double GetXS(int Pid, int Z, int N,  double Ei, double Pf, double Theta, double RadLen);
}; 

#endif

