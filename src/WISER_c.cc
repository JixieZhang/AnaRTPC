#include <stdio.h>
#include <math.h>

namespace WISER 
{

	extern "C" 
	{
		//! Z, N are target info
		//! PART is the particle whose cross section is to be calcualted
		//! E is beam energy in GeV
		//! P is scattered particle momentum in GeV
		//! TH is scattering angle in rad
		//! radlen is total radiation length
		//! xs is in mbarn/GeV/sr, per nucleon

		void wiser_(int* Z, int* N, int* PART, double* Ei, double* Pf, double* ang, 
			double* radlen, double* xs);
	}

	//input:
	//       Z,N: proton and neutron number of the nucleus.
	//       PART: praticle type, 1(pi+), 2(pi-), 3(k+), 4(k-),5(pr);
	//Ei, Ef: incoming electron energy and outgoing pion momentum in GeV;
	//Theta: scattering angle for outgoing particle in radian;
	//RadLen: material thickness in unit of radiation length;
	double WISER(int PART, int Z, int N, double Ei, double Pf, double theta, double radlen)
	{
		double XS;
		wiser_(&Z, &N, &PART, &Ei, &Pf, &theta, &radlen, &XS);
		return XS;
	}


	//input
	//       Z,N: proton and neutron number of the nucleus.
	//       Pid: praticle ID, following the PDG definition
	//       	2212	for p;	+/-211	for pi+/-;  +/-321	for k+/-;
	//Ei, Ef: incoming electron energy and outgoing pion momentum in GeV;
	//Theta: scattering angle for outgoing particle in radian;
	//RadLen: material thickness in unit of radiation length;
	double GetXS(int Pid, int Z, int N, double Ei, double Pf, double Theta, double RadLen)
	{
		if (fabs(RadLen)<1e-8) 
		{
			printf("Error: WISER::GetXS(): radiation length (%f) is not valid\n",RadLen);
			return -1;
		}
		int PartType=2; //pi-

		switch (Pid) 
		{
		case 211: // pi+
			PartType=1;
			break;
		case -211: // pi-
			PartType=2;
			break;
		case 321: // k+
			PartType=3;
			break;
		case -321: // k-
			PartType=4;
			break;
		case 2212: // p
			PartType=5;
			break;
		default:
			printf("Error: WISER::GetXS():Pid(%d) is not recognized\n",PartType);
			return -1;
			break;
		}

		return WISER(Z, N, PartType, Ei, Pf, Theta, RadLen);
	}
}
