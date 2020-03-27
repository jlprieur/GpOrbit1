/**********************************************************************
* "jlp_orbit1_cmpfit.cpp"
* Levenberg-Marquardt method
* Routines using the cminpack programs
* and in particular LMDER_
* to minimize the sum of squares of nonlinear functions 
* with user supplied Jacobian
* 
* Use cmpfit.a library:
* Original public domain version by B. Garbow, K. Hillstrom, J. More'
* (Argonne National Laboratory, MINPACK project, March 1980)
* Tranlation to C Language by S. Moshier (moshier.net) 
*
* JLP
* Version 14/07/2015
**********************************************************************/
#include <stdio.h>
#include <math.h>

#include "cmpfit.h"
#include "gpo_defs.h"  // MINI, MAXI, ...
#include "jlp_orbit1.h"
#include "jlp_kepler.h"

#define NMAX_MEAS 2048

/*
#define DEBUG
*/
/*
#define BOX_CONSTRAINTS
*/

static int static_is_initialized = 0;
static double elmt0[NELEMENTS], dt0;
static int oel_idx0[NELEMENTS]; 
static int n_free_elmt0;
static int nmeas1_8, nrv1_8, nrv2_8, orbit_type1_8;
static double epoch1_8[NMAX_MEAS], epoch_rv1_8[NMAX_MEAS]; 
static double epoch_rv2_8[NMAX_MEAS];
static double rho1_8[NMAX_MEAS], theta1_8[NMAX_MEAS]; 
static double rv1_8[NMAX_MEAS], rv2_8[NMAX_MEAS]; 
static double sigma_rho1_8[NMAX_MEAS];
static double sigma_rv1_8[NMAX_MEAS], sigma_rv2_8[NMAX_MEAS]; 

static int LSQ2_resid_and_gradient(int i_meas, double *resid0,
                                   double *fgrad0);

int fcn(int m0, int n0, double *xx, double *fvec, double **fjac, void *vars);

static int LSQ2_compute_dt0();

/**************************************************************************
* Init static arrays 
**************************************************************************/
int JLP_Orbit1::LSQ2_init_static_arrays(double *elmt, int *oel_idx,
                                        int n_free_elmt)
{
int i, j, m0, n0, maxfev, mode, nprint, info, nfev, njev;

//******************************************************************
// Save input values to static variables:
 n_free_elmt0 = n_free_elmt;
 for(i = 0; i < N10; i++) {
   elmt0[i] = elmt[i];
   oel_idx0[i] = oel_idx[i]; 
   }

//******************************************************************
// Save private variables to static variables:
  nmeas1_8 = nmeas1;
  nrv1_8 = nrv1;
  nrv2_8 = nrv2;
  orbit_type1_8 = orbit_type1;
  for(i = 0; i < nmeas1; i++) {
    epoch1_8[i] = epoch1[i];
    rho1_8[i] = rho1[i];
    theta1_8[i] = theta1[i];
    sigma_rho1_8[i] = sigma_rho1[i];
    }
  for(i = 0; i < nrv1; i++) {
    epoch_rv1_8[i] = epoch_rv1[i];
    rv1_8[i] = rv1[i];
    sigma_rv1_8[i] = sigma_rv1[i];
    }
  for(i = 0; i < nrv2; i++) {
    epoch_rv2_8[i] = epoch_rv2[i];
    rv2_8[i] = rv2[i];
    sigma_rv2_8[i] = sigma_rv2[i];
    }

// Compute dt0, time interval for gradient
  LSQ2_compute_dt0();

  static_is_initialized = 1234;

return(0);
}
/**************************************************************************
* Fit with cminfit routine
* Use Levenberg-Marquardt method
**************************************************************************/
int JLP_Orbit1::LSQ2_fit_with_cmpfit(double *elmt, double *elmt_err,
                                     int *oel_idx, int n_free_elmt,
                                     double *mean_theta_resid,
                                     double *mean_sigma_theta_resid,
                                     double *mean_rho_resid, 
                                     double *mean_sigma_rho_resid, 
                                     double *mean_sigma_rv1_resid,
                                     double *mean_sigma_rv2_resid,
                                     double *chisq2, wxString& Results_str) 
{
int i, j, m0, n0, nfev, njev;
double xx[N10], *fvec;
double perror[N10];                    //  Returned parameter errors
wxString sstr;
mp_result result;
mp_config m_config;
mp_par pars[N10];                        // Parameter constraints 

// {'P', 'T', 'e', 'a', 'W', 'w', 'i', 'K1', 'K2', 'V0'};
#ifdef BOX_CONSTRAINTS
static double xmin[N10] = {1.e-4, 1000., 0., 0.01, -360., -360., -90., 1.e-3, 1.e-3, -1.e+3};
static double xmax[N10] = {1.e+4, 3000., 1., 100., +360., +360., 90., 1.e+3, 1.e+3, 1.e+3};
#endif

memset(&m_config,0,sizeof(m_config));   /* Initialize constraint structure */
memset(pars,0,sizeof(pars));           /* Initialize constraint structure */
memset(&result, 0, sizeof(result));    // Zero results structure 
result.xerror = perror;

m_config.ftol = 10-10;
m_config.xtol = 10-10;
m_config.gtol = 10-10;
m_config.maxfev = 400;
m_config.maxiter = 1000;

for(i = 0; i < N10; i++) {
// Free parameters:
 pars[i].fixed = 0;
// Indicate that derivatives are provided in fcn():
 pars[i].side = 3;
}

#ifdef BOX_CONSTRAINTS
// Eccentricity between xmin[2] and xmax[2]:
pars[2].limited[0] = 1;
pars[2].limited[1] = 1;
pars[2].limits[0] = xmin[2];
pars[2].limits[1] = xmax[2];
#else
  /* No constraints */
#endif

// Save input values to static variables:
LSQ2_init_static_arrays(elmt, oel_idx, n_free_elmt);
 
// m0: total nber of measurements:
  m0 = 2 * nmeas1 + nrv1 + nrv2;

// n0: nber of elements to be determined by the fit: 
  n0 = n_free_elmt0;

// Allocate memory:
  fvec = new double[m0];

// Initial solution, starting values for the fit: 
  for(i = 0; i < n_free_elmt0; i++) { 
    j = oel_idx0[i],
    xx[i] = elmt0[j];
    }

// Launch fit routine:
//  mpfit(fcn, m0, n0, xx, pars, &m_config, NULL, &result); 
  mpfit(fcn, m0, n0, xx, 0, &m_config, NULL, &result); 

// String to be displayed in logbook:
*chisq2 = result.bestnorm;
Results_str.Printf("Successful LSquares2 fit. m0=%d n0=%d niter=%d nfev=%d chisq2=%.2f\n",
                   m0, n0, result.niter, result.nfev, *chisq2);

// Generate full solution
  for(i = 0; i < n_free_elmt0; i++) { 
    j = oel_idx0[i],
    elmt0[j] = xx[i];
    elmt_err[j] = result.xerror[i];
    }

// Compute residuals with that full solution:
LSQ2_compute_residuals(elmt0, oel_idx0, n_free_elmt0, mean_theta_resid,
                       mean_sigma_theta_resid, mean_rho_resid, 
                       mean_sigma_rho_resid, mean_sigma_rv1_resid, 
                       mean_sigma_rv2_resid, chisq2);

sstr.Printf(" chisq=%.2f mean residuals: %.2f (theta) %.3f (rho) %.3f (rv1) %.3f (rv2)\n",
             *chisq2, *mean_sigma_theta_resid, *mean_sigma_rho_resid,
             *mean_sigma_rv1_resid, *mean_sigma_rv2_resid);

Results_str.Append(sstr);

//******************************************************************
// Save to output parameters: 
 for(i = 0; i < N10; i++) {
   elmt[i] = elmt0[i];
   }

return 0;
}
/****************************************************************************
* LSQ2_compute_residuals
* Compute residuals and chisq2
****************************************************************************/
int JLP_Orbit1::LSQ2_compute_residuals(double *elmt, int *oel_idx, 
                                       int n_free_elmt,
                                       double *mean_theta_resid,
                                       double *mean_sigma_theta_resid,
                                       double *mean_rho_resid,
                                       double *mean_sigma_rho_resid,
                                       double *mean_sigma_rv1_resid,
                                       double *mean_sigma_rv2_resid, 
                                       double *chisq2)
{
double sum_resi, sumsq_resi, sum_weight, weight; 
double resid0, fgrad0[N10], wmean;
double K1, K2;
int i, i0; 

// Save input values to static variables:
if(static_is_initialized != 1234) {
  LSQ2_init_static_arrays(elmt, oel_idx, n_free_elmt);
  }
 
*mean_theta_resid = 0.;
*mean_rho_resid = 0.;
*mean_sigma_theta_resid = 0.;
*mean_sigma_rho_resid = 0.;
*mean_sigma_rv1_resid = 0.;
*mean_sigma_rv2_resid = 0.;
*chisq2 = 0.;

// {'P', 'T', 'e', 'a', 'W', 'w', 'i', 'K1', 'K2', 'V0'}
 K1 = elmt0[7];
 K2 = elmt0[8];

// Rho:
if(nmeas1_8 > 0) {
 sum_resi = 0.;
 sumsq_resi = 0.;
 sum_weight = 0.;
 for(i = 0; i < nmeas1_8; i++) {
    LSQ2_resid_and_gradient(i, &resid0, fgrad0);
    i0 = i;
    if(rho1_8[i0] > 0 && sigma_rho1_8[i0] > 0) {
      weight = rho1_8[i0] / sigma_rho1_8[i0];
      sum_resi += resid0;
      sumsq_resi += SQUARE(resid0) / weight;
      *chisq2 += SQUARE(resid0);
      sum_weight += weight;
      }
   }
 if(sum_weight > 0) {
   *mean_rho_resid = sum_resi / sum_weight; 
   *mean_sigma_rho_resid = sqrt(sumsq_resi / sum_weight 
                                - SQUARE(*mean_rho_resid));
   }
 }

// Theta:
if(nmeas1_8 > 0) {
 sumsq_resi = 0.;
 sum_resi = 0.;
 sum_weight = 0.;
 for(i = nmeas1_8; i < 2 * nmeas1_8; i++) {
    LSQ2_resid_and_gradient(i, &resid0, fgrad0);
    i0 = i - nmeas1_8;
    if(sigma_rho1_8[i0] > 0) {
      weight = ABS(rho1_8[i0]) / sigma_rho1_8[i0];
      sum_resi += resid0;
      sumsq_resi += SQUARE(resid0) / weight;
      *chisq2 += SQUARE(resid0);
      sum_weight += weight;
      }
   }
 if(sum_weight > 0) {
   *mean_theta_resid = sum_resi / sum_weight; 
   *mean_sigma_theta_resid = sqrt(sumsq_resi / sum_weight 
                                  - SQUARE(*mean_theta_resid));
   }
 }
// Rv1:
if(nrv1_8 > 0) {
 sumsq_resi = 0.;
 sum_resi = 0.;
 sum_weight = 0.;
 for(i = 2 * nmeas1_8; i < 2 * nmeas1_8 + nrv1_8; i++) {
    LSQ2_resid_and_gradient(i, &resid0, fgrad0);
    i0 = i - 2 * nmeas1_8;
    if(sigma_rv1_8[i0] > 0) {
      weight = K1 / sigma_rv1_8[i0];
      sum_resi += resid0;
      sumsq_resi += SQUARE(resid0) / weight;
      *chisq2 += SQUARE(resid0);
      sum_weight += weight;
      }
   }
 if(sum_weight > 0) {
   wmean = sum_resi / sum_weight;
   *mean_sigma_rv1_resid = sqrt(sumsq_resi / sum_weight - wmean * wmean);
   }
 }
// Rv2:
if(nrv2_8 > 0) {
 sumsq_resi = 0.;
 sum_resi = 0.;
 sum_weight = 0.;
 for(i = 2 * nmeas1_8 + nrv1_8; i < 2 * nmeas1_8 + nrv1_8 + nrv2_8; i++) {
    LSQ2_resid_and_gradient(i, &resid0, fgrad0);
    i0 = i - 2 * nmeas1_8 - nrv1_8;
    if(sigma_rv2_8[i0] > 0) {
      weight = K2 / sigma_rv2_8[i0];
      sum_resi += resid0;
      sumsq_resi += SQUARE(resid0) / weight;
      *chisq2 += SQUARE(resid0);
      sum_weight += weight;
      }
   }
 if(sum_weight > 0) {
   wmean = sum_resi / sum_weight;
   *mean_sigma_rv2_resid = sqrt(sumsq_resi / sum_weight - wmean * wmean);
   }
 }

*chisq2 /= (double)(2 * nmeas1_8 + nrv1_8 + nrv2_8);

return(0);
}
/****************************************************************************
* fcn()
* subroutine fcn for mpfit routine 
*
* INPUT:
*  m0: nber of data points
*  n0: nber of parameters to be computed
*  xx[n0]: array of n0 parameters
*  vars - private data (struct vars_struct *)
* 
*
* OUTPUT:
*  fvec[m0]: deviates
*  fjac[m0]: jacobian 
****************************************************************************/
int fcn(int m0, int n0, double *xx, double *fvec, double **fjac, void *vars)
{      
int i, j;
double resid0, fgrad0[NELEMENTS];
#ifdef DEBUG_1
double chi2 = 0.;
#endif

/* Load current solution: */
 for(i = 0; i < n_free_elmt0; i++) { 
   j = oel_idx0[i],
   elmt0[j] = xx[i];
   }

// Deviates and Jacobian:
 for(i = 0; i < m0; i++) {
   LSQ2_resid_and_gradient(i, &resid0, fgrad0);
   fvec[i] = resid0;
#ifdef DEBUG_1
   chi2 += resid0 * resid0;
#endif

// Only when fjac is non zero, are the user-computed derivatives needed !!!
// Build the Jacobian:
// Derivative of the ith deviate with respect to the jth parameter:
// fjac[j][i] = d(deviate[i]) / d(xx[j]):
   if(fjac != NULL) for(j = 0; j < n0; j++) fjac[j][i] = fgrad0[j];

 }

/* Print statements for debugging */
#ifdef DEBUG_1
      printf("fcn() called: chisq2=.3f\n", chi2);
      for(i = 0; i < n0; i++) printf("%.3f ", xx[i]);
      printf("\n");
#endif

return(0);
}
/*************************************************************************
* LSQ2_resid_and_gradient
* calculates deviates and their derivatives with respect 
* to the free orbital elements
* in each data measurement data point (theta, rho rv1, and rv2)
*
* {'P', 'T', 'e', 'a', 'W', 'w', 'i', 'K1', 'K2', 'V0'};
*
* INPUT:
* elmt0: temporary (full) solution of the orbital elements 
* oel_idx0: index array to go from oelmnt2 to elmt0
*
* i_meas: index of the measurement to be used 
*
* OUTPUT:
* resid0: O-C residual
* fgrad_rho: derivative of the rho residuals (between elmt0 and emt0+d_elmt) 
*           relative to each of the free orbital elements
*************************************************************************/
static int LSQ2_resid_and_gradient(int i_meas, double *resid0, double *fgrad0)
{
double elmt1[NELEMENTS], d_elmt[NELEMENTS];
double elmt_var, eper_var, wtheta, wweight;
double epoch_o, theta_c, rho_c, rv1_c, rv2_c, E_anom;
double K1, K2, fvalue;
int i, j, k, i0;

*resid0 = 0.;
fvalue = 0.;
wweight = 0;

// ckeyword[NELEMENTS] = {'P', 'T', 'e', 'a', 'W', 'w', 'i', 'K1', 'K2', 'V0'};

// 1. Setup for derivative calculation
// Relative element variation for derivative calculation
 elmt_var = 0.01;

// {'P', 'T', 'e', 'a', 'W', 'w', 'i', 'K1', 'K2', 'V0'}
 K1 = elmt0[7];
 K2 = elmt0[8];

// Element derivative:
 d_elmt[0] = elmt_var * elmt0[0];
 d_elmt[1] = elmt_var * elmt0[0];
 d_elmt[2] = elmt_var;
 d_elmt[3] = elmt_var * elmt0[3];
 d_elmt[4] = 1.; 
 d_elmt[5] = 1.; 
 d_elmt[6] = 1.; 
 d_elmt[7] = elmt_var * elmt0[7];
 d_elmt[8] = elmt_var * elmt0[8];
 d_elmt[9] = elmt_var * elmt0[7];

// Relative period variation may be smaller if long time interval:
 eper_var = 0.05 * elmt0[0] / MAXI(dt0, 1.);
 if(eper_var < 0.01) d_elmt[0] = elmt0[0] * eper_var;

// Calling Kepler_Init_External: initialize Kepler parameters with elmt0
  Kepler_Init_External(elmt0, orbit_type1_8);

// Ephemerids and residuals:

//************* Rho: *********************************************
  if(i_meas < nmeas1_8) {
    i0 = i_meas;
    epoch_o = epoch1_8[i0];
    Kepler_Ephemerid1(epoch_o, &theta_c, &rho_c, &E_anom);
    fvalue = rho_c;
    if(sigma_rho1_8[i0] > 0. && rho1_8[i0] > 0.) {
      wweight = rho1_8[i0] / sigma_rho1_8[i0];
      *resid0 = (rho1_8[i0] - rho_c) * wweight;
      }

//************* Theta: *********************************************
    } else if(i_meas >= nmeas1_8 && i_meas < 2 * nmeas1_8) {
    i0 = i_meas - nmeas1_8;
    epoch_o = epoch1_8[i0];
    Kepler_Ephemerid1(epoch_o, &theta_c, &rho_c, &E_anom);
    fvalue = theta_c;
    if(sigma_rho1_8[i0] > 0.) {
      wtheta = theta1_8[i0] - theta_c;
// Adjust quadrant
      while(ABS(wtheta) > 90.) {
       if(wtheta > 0) wtheta -= 180.;
          else wtheta += 180.;
       }
      wweight = ABS(rho1_8[i0]) / sigma_rho1_8[i0];
      *resid0 = wtheta * wweight;
      }

//************* Rv1: *********************************************
    } else if(i_meas >= 2* nmeas1_8 && i_meas < 2 * nmeas1_8 + nrv1_8) {
    i0 = i_meas - 2 * nmeas1_8;
    epoch_o = epoch_rv1_8[i0];
    Kepler_Ephemerid2(epoch_o, &rv1_c, &rv2_c, &E_anom);
    fvalue = rv1_c;
    if(sigma_rv1_8[i0] > 0.) {
       wweight = K1 / sigma_rv1_8[i0];
       *resid0 = (rv1_8[i0] - rv1_c) * wweight;
       }

//************* Rv2: *********************************************
    } else if((i_meas >= 2 * nmeas1_8 + nrv1_8) 
              && (i_meas < 2 * nmeas1_8 + nrv1_8 + nrv2_8)) {
    i0 = i_meas - 2 * nmeas1_8 - nrv1_8;
    epoch_o = epoch_rv2_8[i0];
    Kepler_Ephemerid2(epoch_o, &rv1_c, &rv2_c, &E_anom);
    fvalue = rv2_c;
    if(sigma_rv2_8[i0] > 0.) {
       wweight = K2 / sigma_rv2_8[i0];
       *resid0 = (rv2_8[i0] - rv2_c) * wweight;
       }
    }

// Gradients:

  for(k = 0; k < N10; k++) fgrad0[k] = 0.; 

  if(i_meas < nmeas1_8) {
    i0 = i_meas;
    epoch_o = epoch1_8[i0];
    for(i = 0; i < n_free_elmt0; i++) {
      for(k = 0; k < N10; k++) elmt1[k] = elmt0[k]; 
      j = oel_idx0[i];
      elmt1[j] += d_elmt[j];
      Kepler_Init_External(elmt1, orbit_type1_8);
      Kepler_Ephemerid1(epoch_o, &theta_c, &rho_c, &E_anom);
      fgrad0[i] = wweight * (fvalue - rho_c) / d_elmt[j];
      }
    } else if((i_meas >= nmeas1_8) && (i_meas < 2 * nmeas1_8)) {
    i0 = i_meas - nmeas1_8;
    epoch_o = epoch1_8[i0];
    for(i = 0; i < n_free_elmt0; i++) {
      for(k = 0; k < N10; k++) elmt1[k] = elmt0[k]; 
      j = oel_idx0[i];
      elmt1[j] += d_elmt[j];
      Kepler_Init_External(elmt1, orbit_type1_8);
      Kepler_Ephemerid1(epoch_o, &theta_c, &rho_c, &E_anom);
      wtheta = wweight * (fvalue - theta_c) / d_elmt[j];
// Adjust quadrant
      while(ABS(wtheta) > 90.) {
        if(wtheta > 0) wtheta -= 180.;
         else wtheta += 180.;
         }
      fgrad0[i] = wtheta;
      }
    } else if((i_meas >= 2* nmeas1_8) 
              && (i_meas < 2 * nmeas1_8 + nrv1_8)) {
    i0 = i_meas - 2 * nmeas1_8;
    epoch_o = epoch_rv1_8[i0];
    for(i = 0; i < n_free_elmt0; i++) {
      for(k = 0; k < N10; k++) elmt1[k] = elmt0[k]; 
      j = oel_idx0[i];
      elmt1[j] += d_elmt[j];
      Kepler_Init_External(elmt1, orbit_type1_8);
      Kepler_Ephemerid2(epoch_o, &rv1_c, &rv2_c, &E_anom);
      fgrad0[i] = wweight * (fvalue - rv1_c) / d_elmt[j];
      }
    } else if((i_meas >= 2 * nmeas1_8 + nrv1_8) 
             && (i_meas < 2 * nmeas1_8 + nrv1_8 + nrv2_8)) {
    i0 = i_meas - 2 * nmeas1_8 - nrv1_8;
    epoch_o = epoch_rv2_8[i0];
    for(i = 0; i < n_free_elmt0; i++) {
      for(k = 0; k < N10; k++) elmt1[k] = elmt0[k]; 
      j = oel_idx0[i];
      elmt1[j] += d_elmt[j];
      Kepler_Init_External(elmt1, orbit_type1_8);
      Kepler_Ephemerid2(epoch_o, &rv1_c, &rv2_c, &E_anom);
      fgrad0[i] = wweight * (fvalue - rv2_c) / d_elmt[j];
      }
    }
return(0);
}
/************************************************************************
* Compute dt0, time interval for gradient
************************************************************************/
static int LSQ2_compute_dt0()
{
int i;
double tmin, tmax;
// dt0 : 5 percent of observing time:
 tmin = epoch_rv1_8[0];
 tmax = tmin;
 for(i = 0; i < nrv1_8; i++) {
   tmin = MINI(epoch_rv1_8[i], tmin);
   tmax = MAXI(epoch_rv1_8[i], tmax);
   }
 for(i = 0; i < nrv2_8; i++) {
   tmin = MINI(epoch_rv2_8[i], tmin);
   tmax = MAXI(epoch_rv2_8[i], tmax);
   }
 for(i = 0; i < nmeas1_8; i++) {
   tmin = MINI(epoch1_8[i], tmin);
   tmax = MAXI(epoch1_8[i], tmax);
   }
 dt0 = tmax - tmin;
#ifdef DEBUG
 printf("LSQ2 tmin=%f tmax=%f dt0=%f\n", tmin, tmax, dt0);
#endif

return(0);
}
