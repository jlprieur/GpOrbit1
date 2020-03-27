/****************************************************************************
* Name: jlp_orbit1_lsqfit.cpp
* Fit orbital elements 
* for mixed orbits (visual and radial-velocity measurements)
* From Tokovinin's ORBIT.FOR program
* 
* JLP
* Version 22/06/2015
****************************************************************************/
#include <stdlib.h>   // exit()
#include <math.h>     // cos(), sin(), etc
#include "jlp_orbit1.h"
#include "jlp_kepler.h"
#include "gpo_rw_files.h" // set_idx_free_flags() 

/*
#define DEBUG
*/

static int covar_sum(int n_free_oelmnt2, double alpha[NELEMENTS][NELEMENTS], 
                     double *beta, double *B_vector, int *oel_idx2, double dy, 
                     double *chisq2, double sig2i);

// #define ITER_MAXI 1024
#define ITER_MAXI 4096
 
/*****************************************************************************
* Least-squares fit 
* Fit orbital elements 
* for mixed orbits (visual and radial-velocity measurements)
* From Tokovinin's ORBIT.FOR program
* 
* INPUT:
* max_noelements2: maximum nber of orbital elements
*
* OUTPUT:
* oelements2: orbital elements
* oelements2_err: errors of orbital elements
* chisq2: Chi square
* results_str: string containing information on fit results
*
*****************************************************************************/
int JLP_Orbit1::LSquaresFit1(double *oelements2, double *oelements2_err, 
                             const int max_noelements2, 
                             double *mean_theta_resid, 
                             double *mean_sigma_theta_resid, 
                             double *mean_rho_resid,
                             double *mean_sigma_rho_resid,
                             double *mean_sigma_rv1_resid, 
                             double *mean_sigma_rv2_resid,
                             double *chisq2, wxString &results_str)
{
int i, j, k, iter, status = 1;
double oelmnt2[NELEMENTS], oelmnt2_err[NELEMENTS], scale_oelt2[NELEMENTS];
double alpha[NELEMENTS][NELEMENTS], beta[NELEMENTS];
// kappa: speed convergence, used as a damping factor 
// (kappa = 0.1 for safe convergence, 1.0 for high speed)
double sigma, Ochisq, kappa = 0.1;
// oel_idx2 : index list of elements to be computed
int oel_idx2[NELEMENTS];
// n_free_oelmnt2: nber of orbital elements to be fitted
int n_free_oelmnt2;
/*
* {'P', 'T', 'e', 'a', 'W', 'w', 'i', 'K1', 'K2', 'V0'};
* orbit_type = 1 visual & spectro SB1 (9 elements idx1= 1 1 1 1 1 1 1 1 0 1)
* orbit_type = 2 visual & spectro SB2 (10 elements idx1= 1 1 1 1 1 1 1 1 1 1)
* orbit_type = 3 visual (7 elements: idx2= 1 1 1 1 1 1 1 0 0 0)
* orbit_type = 4 spectroscopic SB1 (6 elements: idx3= 1 1 1 0 0 1 0 1 0 1)
* orbit_type = 5 spectroscopic SB2 (7 elements: idx3= 1 1 1 0 0 1 0 1 1 1)
*/
int idx_free[5][N10];


 if(initialized != 1234 || max_noelements2 != N10) {
   fprintf(stderr, "LSquaresFit1/Error: initialized=%d max_noelements2=%d\n",
           initialized, max_noelements2); 
   return(-1);
   }

if(orbit_type1 != orbit_type_from_elements1) {
  fprintf(stderr, "LSquaresFit1/Error: orbit_type_from_elements1=%d orbit_type1=%d\n",
         orbit_type_from_elements1, orbit_type1);
  results_str = wxT("LSquaresFit1/Error: inconsistant orbit type\n");
  return(-1);
  }

if(orbit_type1 < 1 || orbit_type1 > 5) {
     results_str = wxT("LSquaresFit1/Error: unsupported orbit type\n");
     return(-1);
  }

// Array used to scale orbital elements:
for(i = 0; i < N10; i++) scale_oelt2[i] = 1.;
// From radians to degrees:
for(i = 4; i < 7; i++) scale_oelt2[i] = 180. / PI;

// Initialize oelmnt2 with internal orbital elements:
for(i = 0; i < N10; i++) {
  oelmnt2[i] = oelements1[i];
  }

// Initialize errors:
for(i = 0; i < N10; i++) {
  oelmnt2_err[i] = 0.;
  }

// Index of possibly free parameters:
set_idx_free_flags(idx_free);

// ckeyword[NELEMENTS] = {'P', 'T', 'e', 'a', 'W', 'w', 'i', 'K1', 'K2', 'V0'};
// Build the index list of free elements to be computed:
k = 0;
  for(i = 0; i < N10; i++) {
    if((fixed_flags1[i] == 0) && (idx_free[orbit_type1 - 1][i] == 1)) {
       oel_idx2[k] = i;
       k++;
      }
    }
  n_free_oelmnt2 = k;

// Dummy value to start first iteration:
Ochisq = 1.e+10;

//**************************************************************************
// Main loop:
//**************************************************************************
// nmeas1: nber of measurements 
// n_free_oelmnt2: nber of orbital elements to be fitted
for(iter = 0; iter < ITER_MAXI; iter++) {

#ifdef DEBUG
  printf("LSquaresFit, iter=%d n_free_oelmnt2=%d nmeas1=%d\n", 
          iter, n_free_oelmnt2, nmeas1);
  for(i = 0; i < n_free_oelmnt2; i++) {
    j = oel_idx2[i];
    printf("current: oelmnt2[%d]=%f err=%f\n", j, oelmnt2[j], oelmnt2_err[j]); 
    }
#endif

// Initialize alpha covariance matrix and beta matrix
  for(i = 0; i < n_free_oelmnt2; i++) beta[i] = 0.;

  for(i = 0; i < n_free_oelmnt2; i++) {
    for(k = 0; k < i; k++) {
      alpha[i][k] = 0.;
      }
    }

// Calling Kepler_Init_External: initialize Kepler parameters with oelmnt2
  Kepler_Init_External(oelmnt2, orbit_type1);
  *chisq2 = 0.;

// Build alpha covariance matrix and beta matrix 
// from rv measurements:
  status = LSQ1_covar_from_rv_meas(alpha, beta, oel_idx2, oelmnt2, 
                                  n_free_oelmnt2, chisq2, mean_sigma_rv1_resid,                                  mean_sigma_rv2_resid);

// Build alpha covariance matrix and beta matrix 
// from visual measurements:
  status = LSQ1_covar_from_visual_meas(alpha, beta, oel_idx2, n_free_oelmnt2, 
                                      chisq2, mean_theta_resid, 
                                      mean_sigma_theta_resid, 
                                      mean_rho_resid, mean_sigma_rho_resid);

#ifdef DEBUG
  printf(" iter=%d chisq=%.2f mean residuals: %.3f (theta) %.3f (rho) %.3f (rv1) %.3f (rv2)\n",
          iter, *chisq2, *mean_sigma_theta_resid, *mean_sigma_rho_resid, 
          *mean_sigma_rv1_resid, *mean_sigma_rv2_resid);
#endif

// Stop criterium when chisq2 has increased (relative to previous iteration) 
      if(iter > 1 && *chisq2 >= 0.99999999 * Ochisq) {
#ifdef DEBUG
         printf("New chisq2=%.2f (previous iteration: %.2f)\n", 
                 *chisq2, Ochisq);
         printf("Convergence has thus been reached: exit from loop.\n");
#endif
         status = 0;
         break;
         }

// Covariance matrix has been built in the bottom left corner
// Make covariance matrix symmetric in the upper right corner:
   for(i = 1; i < n_free_oelmnt2; i++) {
     for(k = 0; k < i; k++) {
      alpha[k][i] = alpha[i][k];
      }
     }

// Debug: print diagonal:
#ifdef DEBUG1
   for(i = 0; i < n_free_oelmnt2; i++) {
     j = oel_idx2[i];
     printf("LSQ1/DEBUG: before inversion j=%d alpha[%d][%d]=%f beta[%d]=%f\n",
             j, i, i, alpha[i][i], i, beta[i]);
     }
#endif


   status = LSQ1_GaussInversion(alpha, n_free_oelmnt2, beta);
   if(status) break;
   scale_oelt2[0] = - SQUARE(oelmnt2[0]) / (2. * PI);

// Debug: print diagonal:
#ifdef DEBUG1
   for(i = 0; i < n_free_oelmnt2; i++) {
     j = oel_idx2[i];
     printf("LSQ1/DEBUG: after inversion j=%d alpha[%d][%d]=%f beta[%d]=%f\n",
             j, i, i, alpha[i][i], i, beta[i]);
     }
#endif

// Prepare next iteration: moves towards current solution (contained in beta array)
    Ochisq = *chisq2;

    if(iter == 0) 
       kappa = 0.;
    else
       kappa = 0.1;

     sigma = sqrt(*chisq2 / ( 2. * nmeas1 + nrv1 + nrv2 - n_free_oelmnt2));
     for(i = 0; i < n_free_oelmnt2; i++) {
        j = oel_idx2[i];
// Scale orbital elemnts (angles were obtained in radians,
// but W, w and i are stored in degrees):
        oelmnt2[j] += kappa * scale_oelt2[j] * beta[i];
        oelmnt2_err[j] = ABS(scale_oelt2[j]) * sigma *
                             sqrt(alpha[i][i]); 
        }                           

  } // EOF loop on iter

if(status == 0) {

#ifdef DEBUG
  printf("OK: successful fit final values:\n");
  for(j = 0; j < N10; j++) printf("%d value: %f err: %f\n", 
                                       j, oelmnt2[j], oelmnt2_err[j]); 
#endif

/*
* {'P', 'T', 'e', 'a', 'W', 'w', 'i', 'K1', 'K2', 'V0'};
* orbit_type = 1 visual & spectro SB1 (9 elements idx= 1 1 1 1 1 1 1 1 0 1)
* orbit_type = 2 visual & spectro SB2 (10 elements idx= 1 1 1 1 1 1 1 1 1 1)
* orbit_type = 3 visual (7 elements: idx= 1 1 1 1 1 1 1 0 0 0)
* orbit_type = 4 spectroscopic SB1 (6 elements: idx= 1 1 1 0 0 1 0 1 0 1)
* orbit_type = 5 spectroscopic SB2 (7 elements: idx= 1 1 1 0 0 1 0 1 1 1)
*/

if(orbit_type1 == 3) {
  results_str.Printf("Successful LSquares1 fit. nmeas=%d (iter=%d)\n\
 chisq2=%.2f mean residuals: %.2f+/-%.2f (theta) %.3f+/-%.3f (rho)\n",
                      nmeas1, iter, *chisq2, *mean_theta_resid,
                      *mean_sigma_theta_resid, *mean_rho_resid,
                      *mean_sigma_rho_resid);
  } else {
  results_str.Printf("Successful LSquares1 fit. nmeas=%d nrv1=%d nrv2=%d (iter=%d)\n\
 chisq2=%.2f mean residuals: %.2f+/-%.2f (theta) %.3f+/-%.3f (rho)\n\
 %.3f (rv1) %.3f (rv2)\n",
                      nmeas1, nrv1, nrv2, iter, *chisq2, 
                      *mean_theta_resid, *mean_sigma_theta_resid, 
                      *mean_rho_resid, *mean_sigma_rho_resid,
                      *mean_sigma_rv1_resid, *mean_sigma_rv2_resid);
  } // EOF orbit_type2 == 3 

// Save results before exiting:
 for(i = 0; i < N10; i++) {
  oelements2[i] = oelmnt2[i];
  oelements2_err[i] = oelmnt2_err[i];
  } 

} // EOF status == 0 

return(status);
}
/*************************************************************************
* Build the covariance matrices alpha and beta, 
* and update chisq2
*************************************************************************/
static int covar_sum(int n_free_oelmnt2, double alpha[NELEMENTS][NELEMENTS], 
                     double *beta, double *B_vector, int *oel_idx2, double dy, 
                     double *chisq2, double sig2i)
{ 
double ww;
int i, j, k;

  for(i = 0; i < n_free_oelmnt2; i++) {
    j = oel_idx2[i];
    ww = B_vector[j] * sig2i;
// beta += dy * B_vector * sig2i;
    beta[i] += dy * ww;
    for(k = 0; k <= i; k++) {
      j = oel_idx2[k];
// alpha += B_vector * B_vector * sig2i;
      alpha[i][k] += ww * B_vector[j];
      }
  }
// chisq2
  *chisq2 += dy * dy * sig2i;
return(0);
}
/*************************************************************************
* Build alpha covariance matrix and beta matrix 
* from visual measurements
*
* INPUT:
*  alpha, beta: covariance matrices
*  oel_idx2: index matrix (used when parameters are not free parameters)
*  n_free_oelmnt2: nber of orbital elements to be fitted
*
* OUTPUT:
*  chisq2: new value of chisq2 (integration starting from input value)
*  mean_sigma_theta_resid, mean_sigma_rho_resid: mean residuals for rho and theta
*************************************************************************/
int JLP_Orbit1::LSQ1_covar_from_visual_meas(double alpha[NELEMENTS][NELEMENTS], 
                                           double *beta, int *oel_idx2,
                                           const int n_free_oelmnt2, 
                                           double *chisq2, 
                                           double *mean_theta_resid,
                                           double *mean_sigma_theta_resid,
                                           double *mean_rho_resid,
                                           double *mean_sigma_rho_resid)
{
double wtheta, wrho, sumsq_wtheta, sumsq_wrho; 
double sum_wrho, sum_wtheta;
double sum_weight_rho, sum_weight_theta, wweight;
double sig2i, dy; 
double epoch_o, theta_c, rho_c, E_anom;
double B_vector[NELEMENTS];
int i;

 *mean_theta_resid = 0.;
 *mean_sigma_theta_resid = 0.;
 *mean_rho_resid = 0.;
 *mean_sigma_rho_resid = 0.;

// Return if no data:
 if(nmeas1 == 0) return(1);

 sum_wtheta = 0.;
 sum_wrho = 0.;
 sumsq_wtheta = 0.;
 sumsq_wrho = 0.;
 sum_weight_theta = 0.;
 sum_weight_rho = 0.;

// Data on visual binary orbit
for(i = 0; i < nmeas1; i++) {

// WEIGHTS TAKEN INTO ACCOUNT HERE:
if(sigma_rho1[i] > 0.) {

// Weight by default in Tokovinin's method:
  wweight = ABS(rho1[i]) / sigma_rho1[i];
  sig2i = SQUARE(wweight);

// Compute residuals
 epoch_o = epoch1[i]; 
 Kepler_Ephemerid1(epoch_o, &theta_c, &rho_c, &E_anom);

//*********** Process theta
 wtheta = theta1[i] - theta_c;

// Adjust quadrant
 while(ABS(wtheta) > 90.) {
   if(wtheta > 0) wtheta -= 180.;
   else wtheta += 180.;
   }

 sum_wtheta += wweight * wtheta;
 sumsq_wtheta += wweight * SQUARE(wtheta);
 sum_weight_theta += wweight;

// 6 - Derivatives relative to theta 
// -> B vector: B = d_theta/d_element
  Kepler_Dtheta(epoch_o, theta_c, rho_c, E_anom, B_vector, N7);

// Conversion to radians:
  dy = wtheta * PI / 180.;

// Build the covariance matrices alpha and beta, and update chisq:
  covar_sum(n_free_oelmnt2, alpha, beta, B_vector, oel_idx2, dy, chisq2, sig2i);

//*********** Process Rho
  if(rho1[i] > 0) {
  wrho = rho1[i] - rho_c;
  sum_wrho += wweight * wrho;
  sumsq_wrho += wweight * SQUARE(wrho);
  sum_weight_rho += wweight;

/* DEBUGG
  printf("LSquaresFit, i=%d epoch_o=%f rho_c=%f theta_c=%f\n", 
          i, epoch1[i], rho_c, theta_c);
*/

// 7 - Derivatives relative to rho: 
// -> B vector: B = d_theta/d_element
  Kepler_Drho(epoch_o, theta_c, rho_c, E_anom, B_vector, N7);
  dy = wrho / rho1[i];
// Build the covariance matrices alpha and beta, and update chisq:
  covar_sum(n_free_oelmnt2, alpha, beta, B_vector, oel_idx2, dy, chisq2, sig2i);
  } // EOF rho1[i] > 0

} // EOF if sigma_rho1[i] > 0 .... 

} // EOF loop on i=0,nmeas1

// Mean residuals:
  if(sum_weight_theta > 0) {
    *mean_theta_resid = sum_wtheta / sum_weight_theta;
    *mean_sigma_theta_resid = sqrt(sumsq_wtheta/ sum_weight_theta 
                                   - SQUARE(*mean_theta_resid));
    }
  if(sum_weight_rho > 0) {
    *mean_rho_resid = sum_wrho / sum_weight_rho;
    *mean_sigma_rho_resid = sqrt(sumsq_wrho/ sum_weight_rho
                                 - SQUARE(*mean_rho_resid));
    }

return(0);
}
/*************************************************************************
* Build alpha covariance matrix and beta matrix 
* from radial velocity measurements
*
* INPUT:
*  alpha, beta: covariance matrices
*  oel_idx2: index matrix (used when parameters are not free parameters)
*  n_free_oelmnt2: nber of orbital elements to be fitted
*
* OUTPUT:
*  chisq2: new value of chisq2 (integration starting from input value)
*  mean_sigma_rv1_resid, mean_sigma_rv2_resid: mean residuals for rv1 and rv2 
*************************************************************************/
int JLP_Orbit1::LSQ1_covar_from_rv_meas(double alpha[NELEMENTS][NELEMENTS], 
                                        double beta[NELEMENTS], int *oel_idx2,
                                        double oelmnt2[NELEMENTS],
                                        const int n_free_oelmnt2, 
                                        double *chisq2, 
                                        double *mean_sigma_rv1_resid,
                                        double *mean_sigma_rv2_resid)
{
double wrv1, wrv2, sumsq_wrv1, sumsq_wrv2, ww_rv1, ww_rv2;
double sum_wrv1, sum_wrv2, sum_weight_rv1, sum_weight_rv2; 
double wmean, wweight;
double sig2i, dy; 
double epoch_o, rv1_c, rv2_c, E_anom;
double B_vector[NELEMENTS], K1, K2;
int i;

 *mean_sigma_rv1_resid = 0.;
 *mean_sigma_rv2_resid = 0.;

// Return if no data:
 if(nrv1 == 0 && nrv2 == 0) return(1);

 K1 = oelmnt2[7];
 K2 = oelmnt2[8];

 sum_wrv1 = 0.;
 sum_wrv2 = 0.;
 sumsq_wrv1 = 0.;
 sumsq_wrv2 = 0.;
 sum_weight_rv1 = 0.;
 sum_weight_rv2 = 0.;

/************************************************************************/
// Radial velocity measurements of component 1:
for(i = 0; i < nrv1; i++) {

// WEIGHTS TAKEN INTO ACCOUNT HERE:
 if(sigma_rv1[i] > 0.) {
  wweight = K1 / sigma_rv1[i];
  sig2i = SQUARE(wweight);

// Compute residuals
 epoch_o = epoch_rv1[i]; 
 Kepler_Ephemerid2(epoch_o, &rv1_c, &rv2_c, &E_anom);

 wrv1 = rv1[i] - rv1_c;
 sum_wrv1 += wweight * wrv1;
 sumsq_wrv1 += wweight * SQUARE(wrv1);
 sum_weight_rv1 += wweight;

// Derivatives relative to V1 
// -> B vector: B = 1/K1 * d_V1/d_element
  Kepler_DV1(epoch_o, E_anom, B_vector, N10, &ww_rv1);

/* Debug only to check they are equal:
  printf("DEBUG: rv1_c=%f ww_rv1=%f\n", rv1_c, ww_rv1);
*/

  dy = wrv1 / K1;

// Build the covariance matrices alpha and beta, and update chisq:
  covar_sum(n_free_oelmnt2, alpha, beta, B_vector, oel_idx2, dy, chisq2, sig2i);

  } // EOF sigma_rv1[i] > 0
} // EOF loop on i=0,nrv1

/************************************************************************/
// Radial velocity measurements of component 2:
for(i = 0; i < nrv2; i++) {

// WEIGHTS TAKEN INTO ACCOUNT HERE:
 if(sigma_rv2[i] > 0.) {
  wweight = K2 / sigma_rv2[i];
  sig2i = SQUARE(wweight);

// Compute residuals
 epoch_o = epoch_rv2[i]; 
 Kepler_Ephemerid2(epoch_o, &rv1_c, &rv2_c, &E_anom);

 wrv2 = rv2[i] - rv2_c;
 sum_wrv2 += wweight * wrv2;
 sumsq_wrv2 += wweight * SQUARE(wrv2);
 sum_weight_rv2 += wweight;

// Derivatives relative to V2 
// -> B vector: B = 1/K2 * d_V2/d_element
  Kepler_DV2(epoch_o, E_anom, B_vector, N10, &ww_rv2);

/* Debug only to check they are equal:
  printf("DEBUG: rv2_c=%f ww_rv2=%f\n", rv2_c, ww_rv2);
*/

  dy = wrv2 / K2;

// Build the covariance matrices alpha and beta, and update chisq:
  covar_sum(n_free_oelmnt2, alpha, beta, B_vector, oel_idx2, dy, chisq2, sig2i);

  } // EOF sigma_rv2[i] > 0

} // EOF loop on i=0,nrv2

// Mean residuals:
  if(sum_weight_rv1 > 0) {
    wmean = sum_wrv1 / sum_weight_rv1;
    *mean_sigma_rv1_resid = sqrt(sumsq_wrv1/ sum_weight_rv1 
                                 - SQUARE(wmean) );
    }
  if(sum_weight_rv2 > 0) {
    wmean = sum_wrv2 / sum_weight_rv2;
    *mean_sigma_rv2_resid = sqrt(sumsq_wrv2/ sum_weight_rv2 
                                 - SQUARE(wmean) );
    }

return(0);
}
