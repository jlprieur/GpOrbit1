/****************************************************************************
* Name: jlp_orbit1.cpp
* (for visual orbits)
* 
* JLP
* Version 08/06/2015
****************************************************************************/
#include <stdlib.h>   // exit()
#include <math.h>     // cos(), sin(), etc
#include "jlp_orbit1.h"
#include "jlp_kepler.h"  // Kepler_Init_External

/*
#define DEBUG
*/
/*****************************************************************************
* Constructor without any arguments:
*****************************************************************************/
JLP_Orbit1::JLP_Orbit1()
{
int i;

 orbit_type_from_elements1 = 0;

 ClearMeasurements();

 for(i = 0; i < N10; i++) fixed_flags1[i] = 0;

 initialized = 1234;

return;
}
/*****************************************************************************
* Constructor with measurements and orbital elements 
*****************************************************************************/
JLP_Orbit1::JLP_Orbit1(double *epoch0, double *theta0, double *rho0, 
                       double *sigma_rho0, int nmeas0, 
                       double *oelements0, double *oelements_err0,
                       const int max_noelements0, int orbit_type0)
{
int i;

 orbit_type_from_elements1 = 0;

 ClearMeasurements();

 for(i = 0; i < N10; i++) fixed_flags1[i] = 0;

 initialized = 1234;

// Load measurements:
 LoadMeasurements(epoch0, theta0, rho0, sigma_rho0, nmeas0);

// Set orbital elements
 SetOrbitalElements(oelements0, oelements_err0, max_noelements0, orbit_type0);

return;
}
/*****************************************************************************
* Clear all measurements
*****************************************************************************/
void JLP_Orbit1::ClearMeasurements()
{
// nmeas1: number of (rho, theta) measurements read from input file
// nrv1, nrv2: number of radial velocity measurements read from input file
 nmeas1 = 0;
 nrv1 = 0;
 nrv2 = 0;
 orbit_type1 = 0;
}
/*****************************************************************************
* Load measurements
*****************************************************************************/
int JLP_Orbit1::LoadMeasurements(double *epoch0, double *theta0, double *rho0,
                                 double *sigma_rho0, const int nmeas0)
{
int i;

nmeas1 = 0;

// Possibility of not loading measurements with an "empty" call:
// LoadMeasurements(NULL,NULL,NULL,NULL,0)
 if(nmeas0 <= 0) return(1);

 if(nmeas0 > NMEAS_MAX){
   fprintf(stderr, "LoadMeasurements/Fatal error: nmeas0=%d > NMEAS_MAX=%d\n",
           nmeas0, NMEAS_MAX); 
   exit(-1);
   }

 nmeas1 = nmeas0;

// Transfer data:
 for(i = 0; i < nmeas1; i++) {
   epoch1[i] = epoch0[i];
   theta1[i] = theta0[i];
   rho1[i] = rho0[i];
   sigma_rho1[i] = sigma_rho0[i];
   }

// Set orbit_type1
 set_orbit_type_from_data();

return(0);
}
/*****************************************************************************
* Load radial velocities 
*****************************************************************************/
int JLP_Orbit1::LoadRadialVelocities(double *epoch_rv1_0, double *epoch_rv2_0, 
                                     double *rv1_0, double *rv2_0, 
                                     double *sigma_rv1_0, double *sigma_rv2_0,
                                     const int nrv1_0, const int nrv2_0)
{
int i, nmax;

nrv1 = 0;
nrv2 = 0;

 if(nrv1_0 > NMEAS_MAX || nrv2_0 > NMEAS_MAX){
   fprintf(stderr, "LoadMeasurements/Fatal error: nrv1=%d, nrv2=%d > NMEAS_MAX=%d\n",
           nrv1_0, nrv2_0, NMEAS_MAX); 
   exit(-1);
   }

nrv1 = nrv1_0;
nrv2 = nrv2_0;

if(nrv1_0 != 0){

// Transfer data:
 for(i = 0; i < nrv1; i++) {
   rv1[i] = rv1_0[i];
   epoch_rv1[i] = epoch_rv1_0[i];
   sigma_rv1[i] = sigma_rv1_0[i];
   }
}

if(nrv2_0 != 0) { 

// Transfer data:
 for(i = 0; i < nrv2; i++) {
   rv2[i] = rv2_0[i];
   epoch_rv2[i] = epoch_rv2_0[i];
   sigma_rv2[i] = sigma_rv2_0[i];
   }
}

// Set orbit type:
   set_orbit_type_from_data();

return(0);
}
/*****************************************************************************
* {'P', 'T', 'e', 'a', 'W', 'w', 'i', 'K1', 'K2', 'V0'};
* orbit_type = 1 visual & spectro SB1 (9 elements idx1= 1 1 1 1 1 1 1 1 0 1)
* orbit_type = 2 visual & spectro SB2 (10 elements idx1= 1 1 1 1 1 1 1 1 1 1)
* orbit_type = 3 visual (7 elements: idx2= 1 1 1 1 1 1 1 0 0 0)
* orbit_type = 4 spectroscopic SB1 (6 elements: idx3= 1 1 1 0 0 1 0 1 0 1)
* orbit_type = 5 spectroscopic SB2 (7 elements: idx3= 1 1 1 0 0 1 0 1 1 1)
*****************************************************************************/
int JLP_Orbit1::set_orbit_type_from_data()
{
int status = 0;

orbit_type1 = 0;
if(nmeas1 > 0 && nrv1 > 0 && nrv2 == 0) orbit_type1 = 1;
if(nmeas1 > 0 && nrv1 > 0 && nrv2 > 0) orbit_type1 = 2;
if(nmeas1 > 0 && nrv1 == 0 && nrv2 == 0) orbit_type1 = 3;
if(nmeas1 == 0 && nrv1 > 0 && nrv2 == 0) orbit_type1 = 4;
if(nmeas1 == 0 && nrv1 > 0 && nrv2 > 1) orbit_type1 = 5;
if(orbit_type1 == 0) {
  fprintf(stderr, "set_orbit_type_from_data/Error: nmeas1=%d nrv1=%d nrv2=%d\n", 
          nmeas1, nrv1, nrv2);
  status = -1;
  }

#ifdef DEBUG
printf("set_orbit_type_from_data: orbit_type=%d\n", orbit_type1);
#endif
return(status);
}
/*****************************************************************************
* Set orbital elements
* {'P', 'T', 'e', 'a', 'W', 'w', 'i', 'K1', 'K2', 'V0'};
* orbit_type = 1 visual & spectro SB1 (9 elements idx1= 1 1 1 1 1 1 1 1 0 1)
* orbit_type = 2 visual & spectro SB2 (10 elements idx1= 1 1 1 1 1 1 1 1 1 1)
* orbit_type = 3 visual (7 elements: idx1= 1 1 1 1 1 1 1 0 0 0)
* orbit_type = 4 spectroscopic SB1 (6 elements: idx1= 1 1 1 0 0 1 0 1 0 1)
* orbit_type = 5 spectroscopic SB2 (7 elements: idx1= 1 1 1 0 0 1 0 1 1 1)
*****************************************************************************/
int JLP_Orbit1::SetOrbitalElements(double *oelements0, double *oelements_err0, 
                                   const int max_noelements0,
                                   const int orbit_type0)
{
int i;

// Erase all arrays first to handle error cases:
 for(i = 0; i < N10; i++) {
   oelements1[i] = 0.;
   oelements_err1[i] = 0.;
   }

 orbit_type_from_elements1 = 0;

// Possibility of resetting elements with an "empty" call:
// SetOrbitalElements(NULL,0,0)
 if(orbit_type0 <= 0) return(1);

 if(max_noelements0 != N10) {
   orbit_type_from_elements1 = 0;
   fprintf(stderr, "SetOrbitalElements/Error: max_noelements0 = %d\n",
           max_noelements0);
   exit(-1);
   }

 for(i = 0; i < N10; i++) {
   oelements1[i] = oelements0[i];
   oelements_err1[i] = oelements_err0[i];
   }

 orbit_type_from_elements1 = orbit_type0;

#ifdef DEBUG
printf("SetOrbitalElements: orbit_type=%d\n", orbit_type0);
#endif

return(0);
}
/***********************************************************************
* Get values of orbital elements:
***********************************************************************/
int JLP_Orbit1::GetOrbitalElements(double *oelements0, double *oelements_err0,
                                   const int max_nelmts, int *orbit_type0)
{
int i;

// If no elements have been loaded yet, return error flag:
      if((initialized != 1234) || (max_nelmts != N10)) {
    fprintf(stderr, "GetOrbitalElements/Fatal Error: max_nelmts=%d initialized=%d\n",
             max_nelmts, initialized);
    exit(-1);
    }

  for(i = 0; i < N10; i++) {
     oelements0[i] = oelements1[i];
     oelements_err0[i] = oelements_err1[i];
     }
  *orbit_type0 = orbit_type1;

 return(0);
}
/*************************************************************
* Compute the ephemerids corresponding to the observation epoch:
*
* INPUT:
*  Omega_node (radians), omega_peri (radians), i_incl (radians),
*  e_eccent, T_periastron (years), Period (years),
*  a_smaxis (arcseconds) =  orbital elements
*  epoch_o = Epoch of the ephemerid
*
* OUTPUT:
*  theta_c = ANGOLO DI POSIZIONE CALCOLATO (degrees)
*  rho_c= DISTANZA ANGOLARE CALCOLATA
*
**************************************************************/
int JLP_Orbit1::ComputeEphemerid(float epoch_o, float *theta_c, float *rho_c)
{
double Omega_node, omega_peri, i_incl, e_eccent, T_periastron; 
double Period, a_smaxis, mean_motion, c_tolerance; 
double daa, aa, ab, pp, theta;
double ee, eps, cc, mean_anomaly, eccentric_anomaly, true_anomaly;

// ckeyword[NELEMENTS] = {'P', 'T', 'e', 'a', 'W', 'w', 'i', 'K1', 'K2', 'V0'};
Period = oelements1[0];
T_periastron = oelements1[1];
e_eccent = oelements1[2];
a_smaxis = oelements1[3];
Omega_node = oelements1[4] * PI / 180.;
omega_peri = oelements1[5] * PI / 180.;
i_incl = oelements1[6] * PI / 180.;

mean_motion = 2. * PI / Period;

/*  
* c_tolerance = smallest increment allowed in the iterative process
*                used for solving Kepler's equation
*/
c_tolerance = ABS(1.5E-5 * cos(i_incl)
                      / sqrt((1.0 + e_eccent)/(1.0 - e_eccent)));

/* mean_anomaly = ANOMALIA MEDIA
* true_anomaly = ANOMALIA VERA
* eccentric_anomaly = ANOMALIA ECCENTRICA
*/
   daa = epoch_o - T_periastron;
/*   mean_anomaly = mean_motion * AMOD(daa, Period);
* real function that returns the value of daa modulo Period
* i.e.:
*/
   daa = daa - Period * (float)((int)(daa / Period));
   mean_anomaly = mean_motion * daa;
/* Initialize ee: */
   ee = mean_anomaly + e_eccent*sin(mean_anomaly)
             / sqrt(1.0 + e_eccent * (e_eccent - 2.0 * cos(mean_anomaly)));
/* Iterative resolution of Kepler's equation: */
   do {
    eps = (mean_anomaly - ee) + e_eccent * sin(ee);
    cc = eps / (1.0 - e_eccent * cos(ee));
/* Correction of the current value of ee: */
    ee += cc;
   } while(ABS(cc) > c_tolerance);

   eccentric_anomaly = ee * 0.5;
   aa = sqrt((1.0 + e_eccent) / (1.0 - e_eccent)) * sin(eccentric_anomaly);
/* In FORTRAN:
* ATAN2(y,x) returns the argument alpha of the complex number x + i y,
* expressed in radians, in the range [-PI, +PI],
* such that: tan(alpha) = y / x
* if x > 0  alpha = arctan(y/x)
* if x < 0  alpha = +PI + arctan(y/x) if y > 0
*       or  alpha = -PI + arctan(y/x) if y < 0
* if x = 0  alpha = PI/2 if y > 0 or alpha = -PI/2 if y < 0
*
* Same syntax in C: atan2(y,x) = arctan(y/x) = ATAN2(y,x)
*/
   true_anomaly = 2.0 * atan2(aa, cos(eccentric_anomaly));
   pp = true_anomaly + omega_peri;
   if(pp > 2.*PI) pp -= 2.*PI;
   ab = cos(i_incl) * sin(pp);
   theta = atan2(ab, cos(pp));
   if(theta < 0.0) theta += 2.*PI;

/* Computed position angle: */
   *theta_c = theta + Omega_node;

/* Conversion to degrees: */
 *theta_c *= (180./PI);
/* Put the result into the interval [0,360] */
 if(*theta_c < 0.0) *theta_c += 360.0;
 if(*theta_c > 360.0) *theta_c -= 360.0;

/* Computed separation angle: */
   *rho_c = a_smaxis * (1.0 - e_eccent * cos(ee)) * cos(pp) / cos(theta);
/* DEBUG:
printf("theta=%f theta_c=%f\n", theta, *theta_c);
printf("pp=%f ee=%f e_ccent=%f a_smaxis=%f rho_c=%f\n", pp, ee, e_eccent,
       a_smaxis, *rho_c);
*/

return(0);
}
/*****************************************************************************
* Initialize private Kepler variables with internal orbital elements
*
* double k_P, k_PMU, k_T0, k_SF, k_CF2, k_CF, k_CF3, k_EC;
* double k_A, k_WW, k_CWW, k_SWW, k_W, k_CW, k_SW, k_SI, k_CI, k_TI;
* double k_AA, k_BB, k_FF, k_GG;
*
* ckeyword[NELEMENTS] = {'P', 'T', 'e', 'a', 'W', 'w', 'i', 'K1', 'K2', 'V0'};
*****************************************************************************/
int JLP_Orbit1::Kepler_Init_Internal()
{ int status = -1;
 if(initialized == 1234 && orbit_type1 != 0) {
 status = Kepler_Init_External(oelements1, orbit_type1);
 }
return(status);
}
