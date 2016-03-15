//Header file for helix fit, c++ code!!!
#ifndef _HELIX_FIT_
#define _HELIX_FIT_
///////////////////////////////////////////////////////////////////////

/*------------------------------------------------------------------------\
//Function name: void helix_fit(int PointNum,double szPos[][3], 
//  double& Rho, double& A, double& B,double& Phi, double& Theta, 
//  double& X0, double& Y0,double& Z0, double &DCA, double &chi2);
//
//  Calculate the raidus of the trajectory
//
//Input parameters:
//  PointNum:       number of x-y points
//  szPos[any number>3][3]:  xyz array, 
//OutPut: 
//  Rho:       Radius, positive is clockwise, 
//  A,B:       helix center X,Y position
//  Theta,Phi: theta and phi angle at the initial step
//  X0,Y0,Z0:  vertex position (x,y,z) for inital step
//  DCA:       distance of closest approach to beamline
//  Chi2:      ch2ph+ch2z/(npt-5) 
\------------------------------------------------------------------------*/
void helix_fit(int PointNum,double szPos[][3], double& Rho, double& A, double& B,
  double& Phi, double& Theta, double& X0, double& Y0,double& Z0,  
  double& DCA, double& Chi2, int fit_track_to_beamline=1);

void HelixFit(int PointNum,double szPos[][3], double& R, double& A, double& B, 
  double& Phi_deg, double& Theta_deg, double& Z0, int fit_track_to_beamline=1);

////////////////////////////////////////////////////////////////////////
#endif //#ifdef _HELIX_FIT_
