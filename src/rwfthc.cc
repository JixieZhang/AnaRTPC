

#include "stdio.h"
#include "math.h"

// Version 6

#ifdef __cplusplus
extern "C" {
#endif

  // -------------------------------------------------------------------
  void rwsmav(float r[], float a[], float v[], int n)
  {
    //  Author: Martin Poppe. r[n] = a[n,n]*v[n]

    int i, k, ind;
    // Address in triangular matrix, row ii, column kk
    for(i=1; i<=n; i++)
    {
      r[i-1] = 0.0; 
      for(k=1; k<=n; k++)
      {
	if (i >= k) { ind = (i*i-i)/2 + k; r[i-1] += a[ind-1]*v[k-1]; }
	else        { ind = (k*k-k)/2 + i; r[i-1] += a[ind-1]*v[k-1]; }
      }
    }

  } // End of  void rwsmav(...)





  // -------------------------------------------------------------------
  void rwsmin(float v[], float b[], int n, int m, int* nrank)
  {
    // Author: V. Blobel

    // Obtain solution of linear equations V*X = B with symmetric matrix V
    // and inverse (for m=1) or matrix inversion only (for m=0).

    // V = Symmetric n-by-n matrix in symmetric storage mode 
    //   V(1) = V11,  V(2) = V12,  V(3) = V22,  V(4) = V13, ...
    //   replaced by inverse matrix
    // B = n-vector  (for m=0 use a dummy argument) 
    //   replaced by a solution vector
    // m = see above

    // Method of solution is by elimination selecting the pivot point on the 
    // diagonal each stage. The rank of the matrix is returned in nrank.

    // For nrank != n, all remaining rows and columns of the resulting matrix 
    // V and the corresponding elements of B are set to zero.


#define EPSS 1.0E-6 

    int    i, ii, ni, k, kk, j, jj, jl, jk, lk, l, ij;
    float  vkk, d, e; /* variable 'c' in original fortran file is not used */
    float  dr[2][200];

    // -------- Construct table -----------
    for(i=1; i<=n; i++) dr[0][i-1] = 1.0;
    ni = n; 
    ii = 0;
    for(i=1; i<=n; i++) { ii += i; dr[1][i-1] = fabs(v[ii-1]); }


    // --------- Loop begins ----------
    *nrank = n - ni;
    for (i=1; i<=ni; i++)
    {
      // --- Search for pivot and test for linearity and zero matrix
      k = kk = jj = 0; vkk = 0.0;
      for (j=1; j<=n; j++)
      {
	jj += j;
	if (dr[0][j-1] == 0.0)              break;
	if (fabs(v[jj-1]) <= vkk)           break;
	if (fabs(v[jj-1]) < EPSS*dr[1][j-1]) break;
	vkk = fabs(v[jj-1]);
	k=j; kk=jj;
      }
      if (k == 0) goto label_80;

      // --- Preparation for elimination ---
      *nrank = *nrank + 1;
      dr[0][k-1] = 0.0;
      d = 1.0/v[kk-1];
      v[kk-1] = -d;
      if (m == 1) b[k-1] *= d;
      jk = kk-k;
      jl = 0;

      // --- Elimination ---
      for (j=1; j<=n; j++)
      { 
	if        ( j-k < 0 ) goto label_24; 
	else { if ( j-k == 0) goto label_22;
	else           goto label_26; }

label_22:  jk = kk; jl += j; break;
label_24:  jk++;    goto label_28;
label_26:  jk = jk +j -1;

label_28:
	e = v[jk-1];
	v[jk-1] = d*e;

	if (m==1) b[j-1] -= b[k-1]*e;
	lk = kk-k;
	for (l=1; l<=j; l++)
	{
	  jl++;
	  if        (l-k < 0) goto label_34;
	  else { if (l==k)    goto label_32;
	  else         goto label_36; }
label_32: lk=kk;  break;
label_34: lk++;   goto label_38;
label_36: lk=lk + l - 1;
label_38: v[jl-1] -= v[lk-1]*e;
	}
      } // End of loop over j
    }   // End of loop over i


    // ----------- Change sign --------------
    ij=0;
    for (i=1; i<=n; i++)
      for (j=1; j<=i; j++)
      {
	ij++; v[ij-1] = -v[ij-1];
      }
      goto label_100;


      // --------- Clear rest of the matrix -------------
label_80:
      ij=0;
      for (i=1; i<=n; i++)
      {
	if(m == 1  &&  dr[0][i-1] != 0.0) b[i-1]=0.0;
	for (j=1; j<=i; j++)
	{
	  ij++; 
	  if (dr[0][i-1] + dr[0][j-1] != 0.0) v[ij-1]=0.0;
	  v[ij-1] = -v[ij-1];
	}
      }

label_100:

      return;
  } // End of void rwsmin(float v[], float b[], int n, int m, int* nrank



#define MAX_HITS_ON_CHAIN 200


  //-------------------------------------------------------------------------
  void rwfthc(int npt,                      float rf[MAX_HITS_ON_CHAIN],
    float pf[MAX_HITS_ON_CHAIN],  float wfi[MAX_HITS_ON_CHAIN],  
    float zf[MAX_HITS_ON_CHAIN],  float wzf[MAX_HITS_ON_CHAIN],  
    int iopt, 

    float vv0[],  float ee0[],    float* ch2ph,  float* ch2z,
    float del[MAX_HITS_ON_CHAIN],   float delz[MAX_HITS_ON_CHAIN])


    // ----- Function for fast helix fit. -----

    // A generalization of the TFTHEL routine to allow it to be called from a 
    // routine that contains any list of x and y values xf,yf for a set of npt 
    // points to be fitted.

    // ----- Input: -----

    // npt:         Number of 3-D points to be fit
    // xf[]:        Array of X-values of points to be fit
    // yf[]:        Array of Y-values of points to be fit
    // rf[]:        Array of R-values of points to be fit
    // pf[]:        Array of PHI-values of points to be fit
    // wfi[]:       Array of 1/(sig(rphi))**2 for each point
    // zf[]:        Array of Z-values of points to be fit
    // wzf[]:       Array of 1/(sig(z))**2 for each point
    // iopt:         0 -> Distance**2 =x**2 +y**2 minimized;
    //               1 -> Weighted with 1/SIMA(R*PHI)**2
    //               2 -> Error matrix calculated
    //               3 -> 3-Dimensional iteration

    // ----- Output: -----

    // vv0[5]:      [0] = 1/r*charge, positive if clockwise;
    //              [1] = tan(lambda)  { = dz/ds} tan(angle to (X,Y) PLANE);
    //              [2] = phi0         {0, 2*PI} angle to X-axis at r=d0;
    //              [3] = d0*sign      [cm] minimal distance to Z-axis,
    //                    +ve if axis encircled;
    //              [4] = z0           [cm]    z position at r=d0;
    // ee0[15]:     Inverse of error matrix in triangular form;
    // ch2ph:       chi squared = Sum(phi deviations / errors)^2;
    // ch2z:        chi squared = Sum(z   deviations / errors)^2;
    // del:         Unknown;
    // delz:        Unknown.

    // Note that the number of degrees of freedom = 2*npt-5

    // Based on subroutine CIRCLE.
    // Reference:  "Computer Physics Communications",  Volume 33, P. 329
    // Authors:  N. Chernov, G. Ososkov and M. Poppe

    // Modified by  Fred Weber, 8 Jun 1989.
    // Translated into C by Michael Ispiryan, 2006


  {
#define ITMAX  15       
#define IOWRIT 6
#define EPS    1.0e-16     
#define ONEPI  3.1415927    
#define PI     3.1415927  
#define TWOPI  6.2831854
#define PIBY2  1.57074635


    float   sp2[MAX_HITS_ON_CHAIN],  vv1[5];
    float   sxy[MAX_HITS_ON_CHAIN],   ss0[MAX_HITS_ON_CHAIN]; 
    float   eee[MAX_HITS_ON_CHAIN];
    float   grad[5], cov[15], dv[5];
    float   deln[MAX_HITS_ON_CHAIN],  delzn[MAX_HITS_ON_CHAIN];  

    double  xf[MAX_HITS_ON_CHAIN], yf[MAX_HITS_ON_CHAIN], wf[MAX_HITS_ON_CHAIN];

    double  alf,   a0,   a1,  a2,    a22,   bet,  cur,   dd,   den;
    double  det,   dy,    d2,   f,   fact,  fg,   f1,    g,    gam,   gam0; 
    double  g1,    h,     h2,   p2,  q2,    rm,   rn,    xa,   xb,    xd,   xi;
    double  xm,    xx,    xy,   x1,  x2,    den2, ya,    yb,   yd,    yi,   ym; 
    double  yy,    y1,    y2,   wn,  sa2b2, dd0,   phic;

    int     i,     n,     iter, nrank;

    float   chi2_here,  rr0,       asym, sst, ph0,   check;
    float   aa0,        ome,       gg0,  hh0, ff0,   sums, sumss, sumz, sumsz, sumw; 
    float   denom,      dzds_here, zz0,  eta, dfd,   dfo,  dpd,   dpo,  ggg,   dza;
    float   dzd,        dzo,       chi1;

    //nkb added these variables for the recalculation of dz/ds
    float   kangle, my_phi,     xc,        yc,   xdca;
    float   ydca, xbar, ybar, xpt, ypt, alpha, beta; 
    if (npt <= 2) 
    {
      fprintf(stderr,"BonusHelixFit::rwfthc(): Cannot fit less than 3 points; exiting..\n");  
      return;
    }
    if (npt > MAX_HITS_ON_CHAIN) 
    {
      fprintf(stderr,"BonusHelixFit::rwfthc(): Cannot fit more than %d points; exiting..\n", 
	MAX_HITS_ON_CHAIN);
      return;
    }
    for(i=0; i<npt; i++) 
    {
      xf[i] = rf[i]*cos(pf[i]); 
      yf[i] = rf[i]*sin(pf[i]);   
      wf[i] = wfi[i];
    }

    n = npt; 
    xm = 0.0; 
    ym = 0.0;
    for(i=0; i<15; i++) ee0[i]=0.0;
    for(i=0; i<5;  i++) { grad[i]=0.0; vv0[i]=0.0; }
    chi2_here = 0.0;
    *ch2ph    = 0.0;
    *ch2z     = 0.0;
    for(i=0; i<n; i++) sp2[i]=wf[i]*(rf[i]*rf[i]);

    if(iopt == 0)
    {
      for(i=0; i<n; i++)
      {
	wzf[i]=1.0;
	wf[i]=1.0; 
	xm += xf[i]; 
	ym += yf[i]; 
      }
      rn = 1.0/(double)(n);
    }

    else
    {
      wn=0.0;
      for(i=0; i<n; i++)
      {
	xm += xf[i]*wf[i]; 
	ym += yf[i]*wf[i];   
	wn += wf[i];
      }  
      rn = 1.0/(double)(wn);
    } // End of else


    xm *= rn; 
    ym *= rn; 
    x2=0.0;
    y2=0.0;
    xy=0.0;
    xd=0.0;
    yd=0.0;
    d2=0.0;

    for(i=0; i<n; i++)
    {
      xi  = xf[i] - xm;  
      yi  = yf[i] - ym;
      xx  = xi*xi;       
      yy  = yi*yi;
      x2 += xx*wf[i];
      y2 += yy*wf[i];
      xy += xi*yi*wf[i];
      dd  = xx + yy;
      xd += xi*dd*wf[i];
      yd += yi*dd*wf[i];  
      d2 += dd*dd*wf[i];
    }

    x2 *= rn;  y2 *= rn;  xy *= rn;  d2 *= rn;  xd *= rn;  yd *= rn; 
    f  = 3.0*x2 + y2;
    g  = 3.0*y2 + x2;
    fg = f*g;
    h  = xy + xy;     h2 = h*h;
    p2 = xd*xd;      q2 = yd*yd;   gam0 = x2 + y2; fact = gam0*gam0;
    a2 = (fg-h2-d2)/fact;
    fact *= gam0;
    a1 = (d2*(f+g) - 2.0*(p2+q2))/fact;
    fact *= gam0;
    a0 = (d2*(h2-fg) + 2.0*(p2*g + q2*f) - 4.0*xd*yd*h)/fact;
    a22 = a2 + a2; 
    yb = 1.0E+30; iter=0; xa = 1.0;


    // -------------------- Main iteration ----------------------------
label_103:
    ya = a0 + xa*(a1 + xa*(a2 + xa*(xa-4.0)));
    if (iter >= ITMAX) goto label_105;
    dy = a1 + xa*(a22 + xa*(4.0*xa - 12.0));
    xb = xa - ya/dy;
    if (fabs(ya)    >  fabs(yb)) xb = 0.5*(xb+xa);
    if (fabs(xa-xb) <  EPS)     goto label_105;
    xa = xb; yb = ya; iter++; 
    goto label_103;

label_105:
    gam = gam0*xb;
    f1 = f - gam;
    g1 = g - gam;
    x1 = xd*g1 - yd*h;
    y1 = yd*f1 - xd*h;
    det = f1*g1 - h2;   den2 = 1.0/(x1*x1 + y1*y1 + gam*det*det);
    if (den2 <= 0.0) goto label_999;
    den = sqrtf(den2); cur = det*den + 0.0000000001;
    alf = -(xm*det + x1)*den;
    bet = -(ym*det + y1)*den;
    rm = xm*xm + ym*ym;

    // -------  Calculation of standard circle parameters. NB: cur is 
    // -------  always positive.
    asym = bet*xm - alf*ym;
    sst = 1.0;
    if (asym<0.0) sst = -1.0;
    rr0 = sst*cur;
    if((alf*alf + bet*bet) <= 0.0) goto label_999;
    sa2b2 = 1.0/(sqrtf(alf*alf + bet*bet));
    dd0 = (1.0 - 1.0/sa2b2)/cur;
    phic = asin(alf*sa2b2) + PIBY2;
    if (bet > 0.0) phic = TWOPI - phic;
    ph0 = phic + PIBY2;

    if (rr0 <= 0.0)   ph0 -= ONEPI;
    if (ph0 >  TWOPI) ph0 -= TWOPI;
    if (ph0 <  0.0)   ph0 += TWOPI;

    vv0[0]=rr0;  vv0[2]=ph0;  vv0[3]=dd0;
    //printf("rr0,ph0,dd0 = %f %f %f\n",1/rr0,ph0,dd0);
    check = sst*rr0*dd0;
    if (check == 1.0) { dd0 -= 0.007; vv0[3] = dd0; }

    //  ------- Calculate phi distances to measured points 
    aa0=sst; ome=rr0; gg0=ome*dd0-aa0; hh0=1.0/gg0;
    for(i=0; i<n; i++)
    {
      asym = bet*xf[i] - alf*yf[i];   ss0[i] = 1.0;
      if (asym < 0.0) ss0[i] = -1.0;
      //By jixie: if rf[i] == 0, there is a problem on ff0. 
      if(!rf[i]) ff0 = 0; 
      else ff0 = ome*(rf[i]*rf[i] - dd0*dd0)/(2.0*rf[i]*gg0) + dd0/rf[i];

      if (ff0 < -1.0) ff0 = -1.0;
      if (ff0 >  1.0) ff0 =  1.0;

      del[i] = ph0 + (ss0[i]-aa0)*PIBY2 + ss0[i]*asin(ff0) - pf[i];
      if (del[i] >  ONEPI) del[i] -= TWOPI;
      if (del[i] < -ONEPI) del[i] += TWOPI;
    }

    
    // -------- Fit straight line in S-Z
    for(i=0; i<n; i++)
    {
      eee[i] = 0.5*vv0[0] * 
	       sqrtf(fabs( (rf[i]*rf[i] - vv0[3]*vv0[3]) / 
	                   (1.0-aa0*vv0[0]*vv0[3])         ));
      //printf("eee[%d] = %f\n",i,eee[i]);
      if (eee[i] >  0.9999) 
      { 
	//quiet = FALSE;
	//fprintf(stderr, "+Track circles too much for this code(eee=%f); bad dzds\n",eee[i]);
	//badarg = TRUE;//break;
	//printf("eee[%d] = %f\n",i,eee[i]);
	eee[i] =  0.9999;
      }
      if (eee[i] < -0.9999) 
      {
	//quiet = FALSE;
	//fprintf(stderr, "-Track circles too much for this code(eee=%f); bad dzds\n",eee[i]);
	//badarg = TRUE;//break;
	//printf("eee[%d] = %f\n",i,eee[i]);
	eee[i] = -0.9999;
      }

      sxy[i] = 2.0*asin(eee[i])/ome;
      //printf("original sxy[%d] = %f\n",i,sxy[i]);
    }
    
    //if(badarg)
    {
      /*
      for(i=0; i<n; i++)
      {
      //printf("rf[%d] = %f; pf[%d] = %f; zf[%d] = %f; wfi[%d] = %f;\n",
      //i, rf[i], i, pf[i], i, zf[i],i, wfi[i]); 

      printf("original sxy[%d] = %f, eee = %f\n",i,sxy[i], eee[i]);
      }*/

      //nate's attempt to use the points' arc distance from the dca as the parameter 's'
      //we only use this method if the argument of the arcsin is out of range
      my_phi = ph0 + PI;
      if (vv0[0]<0.0) my_phi+=PI;
      if(my_phi>2.0*PI) my_phi-=2.0*PI;
      xc   = -sin(my_phi)*((-vv0[3])+fabs(1.0/vv0[0]));
      yc   =  cos(my_phi)*((-vv0[3])+fabs(1.0/vv0[0]));
      xdca = -sin(my_phi)*(-vv0[3]);
      ydca =  cos(my_phi)*(-vv0[3]);
      xbar = xdca - xc;
      ybar = ydca - yc;
      //printf("xdca= %.1f ydca= %.1f xc= %.1f yc= %.1f xbar= %.1f ybar= %.1f\n",
      //	    xdca, ydca, xc, yc, xbar, ybar);
      for(i=0; i<n; i++)
      {
	/*//using law of cosines to determine s coordinate
	mydd =  (xf[i]-fabs(dd0)*cos(ph0))*(xf[i]-fabs(dd0)*cos(ph0)) 
	+ (yf[i]-fabs(dd0)*sin(ph0))*(yf[i]-fabs(dd0)*sin(ph0));
	eee[i] = 1 - mydd*rr0*rr0/2;
	if(fabs(eee[i]) > 1.0) 
	{
	//quiet = FALSE;
	//printf("eee[%d] = %f, rad = %f, mydd = %f\n",i,eee[i],1/rr0,sqrtf(mydd));
	//getchar();
	}
	sxy[i] = (1/(sst*rr0))*acos(eee[i]);
	//printf("law of cos sxy[%d] = %f\n",i,sxy[i]);
	*/
	//ksxy = sxy[i];
	xpt = xf[i] - xc;
	ypt = yf[i] - yc;
	alpha = atan2(ypt,xpt);
	beta  = atan2(ybar,xbar); //if(alpha > 2*PI)alpha -= 2*PI;
	//if(alpha < 0)   alpha += 2*PI;
	//if(beta > 2*PI) beta  -= 2*PI;
	//if(beta < 0)    beta  += 2*PI;
	//printf("alpha = %.2f beta = %.2f\n",alpha*180./PI,beta*180./PI);
	sxy[i] = beta - alpha;
	if(sxy[i] > PI) sxy[i] = sxy[i] - 2*PI;
	if(sxy[i] < -PI)sxy[i] = sxy[i] + 2*PI;
	//if(sxy[i] < 0)   sxy[i] += 2*PI;
	kangle = sxy[i];
	sxy[i] = (1/rr0)*sxy[i];

	//printf("[%d] xf= %.1f yf= %.1f xpt= %.1f ypt= %.1f\n",
	//i, xf[i], yf[i], xpt, ypt); 
	//HFILL(9916, ksxy - sxy[i], 0.0, 1.0);
	//printf("%f\n",ksxy-sxy[i]);
      }

    }


    sums = 0.0;
    sumss = 0.0;
    sumz =  0.0;
    sumsz =  0.0;
    sumw = 0.0;
    for(i=0; i<n; i++)
    {
      sumw   += wzf[i];
      sums   += sxy[i]*wzf[i];
      sumss  += sxy[i]*sxy[i]*wzf[i];
      sumz   += zf[i]*wzf[i];
      sumsz  += zf[i]*sxy[i]*wzf[i];
    }

    denom = sumw*sumss - sums*sums;
    if (fabs(denom) < 1.0E-6) 
    { 
      if (denom >= 0.0)  denom =  1.0E-6; 
      else               denom = -1.0E-6; 
    }

    dzds_here = (sumw*sumsz - sums*sumz)  / denom;
    zz0  = (sumss*sumz - sums*sumsz) / denom;
    vv0[1] = dzds_here;  vv0[4] = zz0;

    // --------- Calculation of chi**2
    for(i=0; i<n; i++)
    {
      delz[i]   = zz0 + dzds_here*sxy[i] - zf[i];
      *ch2ph   += sp2[i]*del[i]*del[i];
      *ch2z    += wzf[i]*delz[i]*delz[i];
      chi2_here = *ch2ph + *ch2z;
    }

    if (iopt < 2) return;


    // ----- Calculation of the error matrix -------
    for(i=0; i<n; i++)
    {
      ff0 = ome*(rf[i]*rf[i] - dd0*dd0) / (2.0*rf[i]*gg0) + dd0/rf[i];
      if (ff0 >  0.9999)  ff0 =  0.9999;
      if (ff0 < -0.9999)  ff0 = -0.9999; 
      eta = ss0[i] / sqrtf(fabs(1.0+ff0)*(1.0-ff0));
      dfd = (1.0 + hh0*hh0*(1.0-ome*ome*rf[i]*rf[i])) / (2.0*rf[i]);
      dfo = -aa0*(rf[i]*rf[i] - dd0*dd0)*hh0*hh0 / (2.0*rf[i]);
      dpd = eta*dfd;  dpo = eta*dfo; 
      // --- Derivatives of z component
      ggg = eee[i] / sqrtf(fabs( (1.0+eee[i])*(1.0-eee[i])));
      dza = sxy[i];
      check = rf[i]*rf[i] - vv0[3]*vv0[3];
      if (check == 0.0) check = 2.0*0.007;
      dzd = 2.0*(vv0[1]/vv0[0]) * fabs(ggg) * (0.5*aa0*vv0[0] /
	(1.0 - aa0*vv0[3]*vv0[0]) - vv0[3]/check);
      dzo = -vv0[1]*sxy[i]/vv0[0] + vv0[1]*ggg/(vv0[0]*vv0[0]) * 
	(2.0 + aa0*vv0[0]*vv0[3]/(1.0 - aa0*vv0[0]*vv0[3]));

      // ---- Error matrix
      ee0[0]  += sp2[i]*dpo*dpo + wzf[i]*dzo*dzo;
      ee0[1]  +=                wzf[i]*dza*dzo;
      ee0[2]  +=                wzf[i]*dza*dza;
      ee0[3]  += sp2[i]*dpo;
      ee0[5]  += sp2[i];
      ee0[6]  += sp2[i]*dpo*dpd + wzf[i]*dzo*dzd;
      ee0[8]  += sp2[i]*dpd; 
      ee0[9]  += sp2[i]*dpd*dpd + wzf[i]*dzd*dzd;
      ee0[10] +=                    wzf[i]*dzo;
      ee0[11] +=                    wzf[i]*dza;
      ee0[13] +=                    wzf[i]*dzd;
      ee0[14] +=                    wzf[i];

      // --- Gradient vector
      grad[0] += -del[i]*sp2[i]*dpo - delz[i]*wzf[i]*dzo; 
      grad[1] +=                  - delz[i]*wzf[i]*dza;
      grad[2] += -del[i]*sp2[i];
      grad[3] += -del[i]*sp2[i]*dpd - delz[i]*wzf[i]*dzd;
      grad[4] +=                  - delz[i]*wzf[i];   
    } // End of for(i...)


    if (iopt < 3) return;


    // --------------- Newton's next guess
    for(i=0; i<15; i++) cov[i] = ee0[i];

    rwsmin(cov, vv1, 5, 0, &nrank);
    rwsmav(dv,  cov, grad, 5);

    for(i=0; i<5; i++) vv1[i] = vv0[i] + dv[i];

    //------- New differences in phi and z
    gg0 = vv1[0]*vv1[3] - aa0;
    for(i=0; i<n; i++)
    {
      ff0 = vv1[0]*(rf[i]*rf[i] - vv1[3]*vv1[3]) / 
	(2.0*rf[i]*gg0) + vv1[3]/rf[i];

      if (ff0 >  1.0) ff0 =  1.0;
      if (ff0 < -1.0) ff0 = -1.0;

      deln[i] = vv1[2] + (ss0[i]-aa0)*PIBY2 + ss0[i]*asin(ff0) - pf[i];
      if (deln[i] >  ONEPI) deln[i] -= TWOPI;
      if (deln[i] < -ONEPI) deln[i] += TWOPI;
      eee[i] = 0.5*vv1[0]*sqrtf(fabs( (rf[i]*rf[i] - vv1[3]*vv1[3]) / 
	(1.0 - aa0*vv1[0]*vv1[3]) ));
      if (eee[i] >  0.9999) eee[i] =  0.9999;
      if (eee[i] < -0.9999) eee[i] = -0.9999;
      sxy[i]   = 2.0*asin(eee[i]) / vv1[0];
      delzn[i] = vv1[4] + vv1[1]*sxy[i] - zf[i];
    } 

    // ---------- Calculation of chi**2
    chi1 = *ch2ph = *ch2z = 0.0;
    for(i=0; i<n; i++)
    {
      chi1  += sp2[i]*deln[i]*deln[i] + wzf[i]*delzn[i]*delzn[i];
      *ch2ph += sp2[i]*deln[i]*deln[i];
      *ch2z  += wzf[i]*delzn[i]*delzn[i]; 
    }

    if (chi1 < chi2_here) { for(i=0; i<5; i++) vv0[i] = vv1[i]; }
    return;

    // --- Jump here if something goes crazy
label_999:
    *ch2ph = 1.0E+30; *ch2z = 1.0E+30;

    return;
  }


#ifdef __cplusplus
}
#endif
