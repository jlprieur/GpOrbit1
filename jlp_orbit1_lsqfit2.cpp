/****************************************************************************
* Name: jlp_orbit1_lsqfit2.cpp
* Fit orbital elements by Levenberg-Marquardt method
* for mixed orbits (visual and radial-velocity measurements)
* From Tokovinin's IDL orbit.pro procedure
* 
* JLP
* Version 22/06/2015
****************************************************************************/
#include <stdlib.h>   // exit()
#include <math.h>     // cos(), sin(), etc

#include "gpo_defs.h"  // MINI, MAXI, ...
#include "gpo_rw_files.h" // set_idx_free_flags()
#include "jlp_orbit1.h"

#define DEBUG

/*****************************************************************************
* Least-squares fit 
* Fit orbital elements by Levenberg-Marquardt method
* for mixed orbits (visual and radial-velocity measurements)
* From Tokovinin's IDL orbit.pro procedure
* 
* ckeyword[NELEMENTS] = {'P', 'T', 'e', 'a', 'W', 'w', 'i', 'K1', 'K2', 'V0'};
*
* INPUT:
* max_nelmt0: max nber of orbital elements
*
* OUTPUT:
* elmt0[N10]: orbital elements
* elmt0_err[N10]: errors of orbital elements
* chisq2: Chi square
* results_str: string containing information on fit results
*
*****************************************************************************/
int JLP_Orbit1::LSquaresFit2(double *elmt0, double *elmt0_err,
                             const int max_nelmt0, 
                             double *mean_theta_resid,      
                             double *mean_sigma_theta_resid,
                             double *mean_rho_resid,
                             double *mean_sigma_rho_resid,
                             double *mean_sigma_rv1_resid,
                             double *mean_sigma_rv2_resid,
                             double *chisq2, wxString &results_str)
{
int i, j, k, iter, status;
double scale_oelt2[NELEMENTS];
int idx_free[5][N10];
// oel_idx0 : index list of elements to be computed
int oel_idx0[NELEMENTS];
// n_free_elmt0: nber of orbital elements to be fitted
int n_free_elmt0;

 if(initialized != 1234 || max_nelmt0 != N10) {
   fprintf(stderr, "LSquaresFit2/Error: initialized=%d max_nelmt0=%d\n",
           initialized, max_nelmt0);
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

// Initialize elmt0 with internal orbital elements:
for(i = 0; i < N10; i++) {
  elmt0[i] = oelements1[i];
  }

// Initialize errors:
for(i = 0; i < N10; i++) {
  elmt0_err[i] = 0.;
  }

/*
* {'P', 'T', 'e', 'a', 'W', 'w', 'i', 'K1', 'K2', 'V0'};
* orbit_type = 1 visual & spectro SB1 (9 elements idx1= 1 1 1 1 1 1 1 1 0 1)
* orbit_type = 2 visual & spectro SB2 (10 elements idx1= 1 1 1 1 1 1 1 1 1 1)
* orbit_type = 3 visual (7 elements: idx2= 1 1 1 1 1 1 1 0 0 0)
* orbit_type = 4 spectroscopic SB1 (6 elements: idx3= 1 1 1 0 0 1 0 1 0 1)
* orbit_type = 5 spectroscopic SB2 (7 elements: idx3= 1 1 1 0 0 1 0 1 1 1)
*/
// Index of possibly free parameters:
set_idx_free_flags(idx_free);

// ckeyword[NELEMENTS] = {'P', 'T', 'e', 'a', 'W', 'w', 'i', 'K1', 'K2', 'V0'};
// Build the index list of free elements to be computed:
k = 0;
  for(i = 0; i < N10; i++) {
    if((fixed_flags1[i] == 0) && (idx_free[orbit_type1 - 1][i] == 1)) {
       oel_idx0[k] = i;
       k++;
      }
    }
  n_free_elmt0 = k;

// Array used to scale orbital elements:
for(i = 0; i < N10; i++) scale_oelt2[i] = 1.; 
for(i = 4; i < 7; i++) scale_oelt2[i] = 180. / PI; 

// Fit elements to data with Levenberg-Marquardt method:
 status = LSQ2_fit_with_cmpfit(elmt0, elmt0_err, oel_idx0, n_free_elmt0, 
                               mean_theta_resid, mean_sigma_theta_resid, 
                               mean_rho_resid, mean_sigma_rho_resid,
                               mean_sigma_rv1_resid, mean_sigma_rv2_resid,
                               chisq2, results_str);

return(status);
}
