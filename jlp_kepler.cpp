/****************************************************************************
* Name: jlp_orbit1_kepler.cpp
* (for visual orbits)
* From A. Tokovinin (version of 1991) 
* 
* JLP
* Version 08/06/2015
****************************************************************************/
#include <stdio.h>  
#include <stdlib.h>   // exit()
#include <math.h>     // cos(), sin(), etc
#include "gpo_defs.h"

/* Contained here and declared in "jlp_orbit1.h"
    int Kepler_Init_Internal(const int orbit_type0);
    int Kepler_Init_External(double *oelmts1, const int orbit_type0);
    int Kepler_Ephemerid1(double epoch_o, double *theta_c, double *rho_c);
    int Kepler_Ephemerid2(const double epoch_o, double *rv1_c, 
                          double *rv2_c, double *E_anom_c)
    int Kepler_Dtheta(double *B_vector);
    int Kepler_Drho(double *B_vector);
*/

// Kepler private parameters:
static double k_P, k_PMU, k_T0, k_SF, k_CF2, k_CF, k_CF3, k_EC;
static double k_A, k_WW, k_CWW, k_SWW, k_W, k_CW, k_SW, k_SI, k_CI, k_TI;
static double k_AA, k_BB, k_FF, k_GG;
static double k_K1, k_K2, k_V0;

/*****************************************************************************
* Initialize private Kepler variables with external orbital elements 
*
* double k_P, k_PMU, k_T0, k_SF, k_CF2, k_CF, k_CF3, k_EC; 
* double k_A, k_WW, k_CWW, k_SWW, k_W, k_CW, k_SW, k_SI, k_CI, k_TI;
* double k_AA, k_BB, k_FF, k_GG;
*
* ckeyword[NELEMENTS] = {'P', 'T', 'e', 'a', 'W', 'w', 'i', 'K1', 'K2', 'V0'};
* {'P', 'T', 'e', 'a', 'W', 'w', 'i', 'K1', 'K2', 'V0'};
* orbit_type = 1 visual & spectro SB1 (9 elements idx1= 1 1 1 1 1 1 1 1 0 1)
* orbit_type = 2 visual & spectro SB2 (10 elements idx1= 1 1 1 1 1 1 1 1 1 1)
* orbit_type = 3 visual (7 elements: idx2= 1 1 1 1 1 1 1 0 0 0)
* orbit_type = 4 spectroscopic SB1 (6 elements: idx3= 1 1 1 0 0 1 0 1 0 1)
* orbit_type = 5 spectroscopic SB2 (7 elements: idx3= 1 1 1 0 0 1 0 1 1 1)
*
*****************************************************************************/
int Kepler_Init_External(double *oelmts1, const int orbit_type0)
{
double i_radians;

if((orbit_type0 <= 0) || (orbit_type0 > 5)) {
  fprintf(stderr, "Kepler_Init_External/Fatal error: orbit_type0 = %d\n", 
          orbit_type0);
  exit(-1);
  }

  k_P = oelmts1[0];  // period
  k_PMU = 2. * PI / k_P;

  k_T0 = oelmts1[1]; // T_periastron

  k_SF = oelmts1[2];  // eccentricity (e)
  k_CF2 = 1. - k_SF * k_SF; // 1 - e^2
  k_CF = sqrt(k_CF2);       // sqrt(1 - e^2)
  k_CF3 = k_CF * k_CF2;     // (1 - e^2)^(3/2)
  k_EC = sqrt( (1. + k_SF) / (1. - k_SF) );  // sqrt((1 + e)(1 - e))

  k_A = oelmts1[3];               // semi major-axis (a)

  k_WW = oelmts1[4] * PI / 180.;  // Omega (in radians)
  k_CWW = cos(k_WW);              // cos(Omega)
  k_SWW = sin(k_WW);              // sin(Omega)

  k_W = oelmts1[5] * PI / 180.;   // omega (in radians)
  k_CW = cos(k_W);                // cos(omega)
  k_SW = sin(k_W);                // sin(omega)

  i_radians = oelmts1[6] * PI /180.; // inclination (in radians)
  k_SI = sin(i_radians);             // sin(inclination)
  k_CI = cos(i_radians);             // cos(inclination)
  if(k_CI == 0.) k_CI = 1.e-9;
  k_TI = k_SI / k_CI;                // tan(inclination)

// Thiele-van-den-Bos elements
  k_AA = k_A * (k_CW * k_CWW - k_SW * k_SWW * k_CI);
  k_BB = k_A * (k_CW * k_SWW + k_SW * k_CWW * k_CI);
  k_FF = k_A * (-k_SW * k_CWW - k_CW * k_SWW * k_CI);
  k_GG = k_A * (-k_SW * k_SWW + k_CW * k_CWW * k_CI);

// For radial velocities:
/*
* {'P', 'T', 'e', 'a', 'W', 'w', 'i', 'K1', 'K2', 'V0'};
* orbit_type = 1 visual & spectro SB1 (9 elements idx1= 1 1 1 1 1 1 1 1 0 1)
* orbit_type = 2 visual & spectro SB2 (10 elements idx1= 1 1 1 1 1 1 1 1 1 1)
* orbit_type = 3 visual (7 elements: idx2= 1 1 1 1 1 1 1 0 0 0)
* orbit_type = 4 spectroscopic SB1 (6 elements: idx3= 1 1 1 0 0 1 0 1 0 1)
* orbit_type = 5 spectroscopic SB2 (7 elements: idx3= 1 1 1 0 0 1 0 1 1 1)
*/
  if(orbit_type0 != 0 && orbit_type0 != 3) {
    k_K1 = oelmts1[7];
    k_K2 = oelmts1[8];
    k_V0 = oelmts1[9];
    } else {
    k_K1 = 0.;
    k_K2 = 0.;
    k_V0 = 0.;
    } 

return(0);
}

/*****************************************************************************
* Sove Kepler equation 
* and return eccentric and true anomaly:
*****************************************************************************/
int Solve_Kepler(const double epoch_o, double *E_anomaly,
                             double *V_anomaly) 
{
double phase, dt, ww, M_anom, E_anom1, E_anom, V_anom;
int i, isafe;

// Phase for the epoch:
  dt = epoch_o - k_T0;
  ww = dt / k_P;
  phase = ww - (int)ww;
  if(phase < 0) phase += 1.;

// Solve Kepler equation
// E - e sin(E) = M
// with M = n (t - t0) : mean anomaly
// E: eccentric anomaly
  M_anom = phase * 2. * PI;
  E_anom = M_anom;

// NB: k_SF is the eccentricity
  isafe = 256;
  for(i = 0; i < isafe; i++) {
//DDEBUG: printf("epoch=%f E_anom = %f\n", epoch_o, E_anom);
    E_anom1 = E_anom 
           + (M_anom + k_SF * sin(E_anom) - E_anom)/(1.0 - k_SF * cos(E_anom));
    if(ABS(E_anom1 - E_anom) < 1.e-6) break; 
    E_anom = E_anom1;
    }
 
// Eccentric anomaly:
 *E_anomaly = E_anom;

// True anomaly: v such that tan(v/2) = sqrt((1 + e)/(1-e)) * tan(E/2) 
 *V_anomaly = 2. * atan(k_EC * tan(E_anom/2.));

return(0);
}

/*****************************************************************************
* Compute ephemerid rho and theta for a given epoch
*****************************************************************************/
int Kepler_Ephemerid1(const double epoch_o, double *theta_c, 
                                  double *rho_c, double *E_anom_c)
{
double E_anom, V_anom;
double U, CU, TMW, R;

// Sove Kepler equation and return eccentric and true anomaly:
 Solve_Kepler(epoch_o, &E_anom, &V_anom); 

 U = V_anom + k_W;
 CU = cos(U);

// k_CI: cos(i)
 TMW = atan( k_CI * tan(U) );
 if(CU < 0.) TMW += PI;

 R = k_A * (1. - k_SF * cos(E_anom));
 *rho_c = R * CU / cos(TMW);
 *theta_c = (TMW + k_WW) * 180. / PI;
 if(*theta_c > 360.) *theta_c -= 360.;
 if(*theta_c < 0.) *theta_c += 360.;

// Save eccentric anomaly (used for gradients...)
 *E_anom_c = E_anom;

return(0);
}
/*****************************************************************************
* Compute ephemerid rv1 and rv2 for a given epoch 
*****************************************************************************/
int Kepler_Ephemerid2(const double epoch_o, double *rv1_c, 
                                  double *rv2_c, double *E_anom_c)
{
double E_anom, V_anom;
double U, CU, A1;

// Sove Kepler equation and return eccentric and true anomaly:
 Solve_Kepler(epoch_o, &E_anom, &V_anom); 

 U = V_anom + k_W;
 CU = cos(U);
 A1 = k_SF * k_CW + CU;

// Radial velocities for components 1 and 2:

 *rv1_c = k_V0 + k_K1 * A1;

 *rv2_c = k_V0 - k_K2 * A1;

// Save eccentric anomaly (used for gradients...)
 *E_anom_c = E_anom;

return(0);
}
/*****************************************************************************
* Compute B vector (derivative of theta relative to orbital elements 
* B(j) = dTETA/dEL(j), RES1=TETA
*
* Private parameters:
* double k_P, k_PMU, k_T0, k_SF, k_CF2, k_CF, k_CF3, k_EC; 
* double k_A, k_WW, k_CWW, k_SWW, k_W, k_CW, k_SW, k_SI, k_CI, k_TI;
* double k_AA, k_BB, k_FF, k_GG;
*****************************************************************************/
int Kepler_Dtheta(const double epoch_o, const double theta_c, 
                              const double rho_c, const double E_anom,
                              double *B, const int n_B) 
{
double R, TMW, dt, V_anom;
double A1, A2, A3, A4;

if(n_B < N7) {
  fprintf(stderr, "Kepler_Dtheta/Fatal error: n_B=%d < 7 !\n", n_B);
  exit(-1);
  }


  A1 = SQUARE(k_A / rho_c);
  R = k_A * (1. - k_SF * cos(E_anom));
  A2 = R / k_A;
  A3 = A1 * k_CF * k_CI;
  TMW = theta_c * PI / 180. - k_WW;
  A4 = cos(TMW) * sin(TMW);

  dt = epoch_o - k_T0;
  B[0] = A3 * dt;
  B[1] = -A3 * k_PMU;
  B[2] = A1 * k_CI * sin(E_anom) * (A2 + k_CF2) / k_CF;
  B[3] = 0.;
  B[4] = 1.;
// True anomaly: v such that tan(v/2) = sqrt((1 + e)/(1-e)) * tan(E/2) 
  V_anom = 2. * atan(k_EC * tan(E_anom/2.));
  B[5] = k_CI * SQUARE(cos(TMW) / cos(V_anom + k_W));
  B[6] = -k_TI * A4;
  B[7] = 0.;
  B[8] = 0.;
  B[9] = 0.;

return(0);
}
/*****************************************************************************
* Compute B vector (derivative of drho/rho relative to orbital elements 
* B(j)=1/RO*dRO/dEL(j)
*
* Private parameters:
* double k_P, k_PMU, k_T0, k_SF, k_CF2, k_CF, k_CF3, k_EC; 
* double k_A, k_WW, k_CWW, k_SWW, k_W, k_CW, k_SW, k_SI, k_CI, k_TI;
* double k_AA, k_BB, k_FF, k_GG;
*****************************************************************************/
int Kepler_Drho(const double epoch_o, const double theta_c, 
                            const double rho_c, const double E_anom,
                            double *B, const int n_B) 
{
double R, TMW, dt, A1, A2, A3, A4, A5;

if(n_B < N7) {
  fprintf(stderr, "Kepler_Drho/Fatal error: n_B=%d < 7 !n", n_B);
  exit(-1);
  }

  R = k_A * (1. - k_SF * cos(E_anom));
  A1 = SQUARE(k_A / R);
  A2 = R / k_A;
  TMW = theta_c * PI / 180. - k_WW;
  A4 = cos(TMW) * sin(TMW);
  A5 = A4 * k_TI * k_SI;
  A3 = A1 * (k_SF * sin(E_anom) - A5 * k_CF);

  dt = epoch_o - k_T0;
  B[0] = A3 * dt;
  B[1] = -A3 * k_PMU;
  B[2] = -A1 * ((cos(E_anom) - k_SF) * k_CF 
               + A4 * (A2 + k_CF2) * sin(E_anom)) / k_CF;
  B[3] = 1. / k_A;
  B[4] = 0.;
  B[5] = -A5;
  B[6] = -k_TI * SQUARE(sin(TMW));
  B[7] = 0.;
  B[8] = 0.;
  B[9] = 0.;

return(0);
}
/*****************************************************************************
* Compute B vector (derivative of dV1/K1 relative to orbital elements 
* B(j)=1/K1 * dV1/dEL(j)
* and also compute rv1_c
*
* Private parameters:
* double k_P, k_PMU, k_T0, k_SF, k_CF2, k_CF, k_CF3, k_EC; 
* double k_A, k_WW, k_CWW, k_SWW, k_W, k_CW, k_SW, k_SI, k_CI, k_TI;
* double k_AA, k_BB, k_FF, k_GG;
* double k_V0, k_K1, k_K2;
*****************************************************************************/
int Kepler_DV1(const double epoch_o, const double E_anom,
                           double *B, const int n_B, double *rv1_c) 
{
double V_anom, CV, SV, U, CU, SU, A1, dt; 

if(n_B != N10) {
  fprintf(stderr, "Kepler_DV1/Fatal error: n_B=%d != 10 !n", n_B);
  exit(-1);
  }

// True anomaly: v such that tan(v/2) = sqrt((1 + e)/(1-e)) * tan(E/2) 
  V_anom = 2. * atan(k_EC * tan(E_anom/2.));
  CV = cos(V_anom);
  SV = sin(V_anom);
  U = V_anom + k_W;
  CU = cos(U);
  SU = sin(U);
  A1 = SQUARE(1. + k_SF * CV) * SU / k_CF3;

  dt = epoch_o - k_T0;
  B[0] = - A1 * dt;
  B[1] = A1 * k_PMU;
  B[2] = k_CW - SU * SV * (2. + k_SF * CV) / k_CF2; 
  B[3] = 0.;
  B[4] = 0.;
  B[5] = -SU - k_SF * k_SW;
  B[6] = 0.;
  B[7] = CU + k_SF * k_CW;
  *rv1_c = k_V0 + k_K1 * B[7]; 
  B[7] /= k_K1;
  B[8] = 0.;
  B[9] = 1. / k_K1;

return(0);
}
/*****************************************************************************
* Compute B vector (derivative of dV2/K2 relative to orbital elements 
* B(j)=1/K2 * dV2/dEL(j)
* and also compute rv2_c
*
* Private parameters:
* double k_P, k_PMU, k_T0, k_SF, k_CF2, k_CF, k_CF3, k_EC; 
* double k_A, k_WW, k_CWW, k_SWW, k_W, k_CW, k_SW, k_SI, k_CI, k_TI;
* double k_AA, k_BB, k_FF, k_GG;
* double k_V0, k_K1, k_K2;
*****************************************************************************/
int Kepler_DV2(const double epoch_o, const double E_anom,
                           double *B, const int n_B, double *rv2_c) 
{
double V_anom, CV, SV, U, CU, SU, A1, dt; 

if(n_B != N10) {
  fprintf(stderr, "Kepler_DV1/Fatal error: n_B=%d != 10 !n", n_B);
  exit(-1);
  }

// True anomaly: v such that tan(v/2) = sqrt((1 + e)/(1-e)) * tan(E/2) 
  V_anom = 2. * atan(k_EC * tan(E_anom/2.));
  CV = cos(V_anom);
  SV = sin(V_anom);
  U = V_anom + k_W;
  CU = cos(U);
  SU = sin(U);
  A1 = SQUARE(1. + k_SF * CV) * SU / k_CF3;

  dt = epoch_o - k_T0;
  B[0] = A1 * dt;
  B[1] = -A1 * k_PMU;
  B[2] = -k_CW + SU * SV * (2. + k_SF * CV) / k_CF2; 
  B[3] = 0.;
  B[4] = 0.;
  B[5] = SU + k_SF * k_SW;
  B[6] = 0.;
  B[7] = -CU - k_SF * k_CW;
  *rv2_c = k_V0 + k_K2 * B[7]; 
  B[8] = B[7] / k_K2;
  B[7] = 0.;
  B[9] = 1. / k_K2;

return(0);
}
