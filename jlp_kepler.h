/****************************************************************************
* Name: jlp_kepler.h
* 
* JLP
* Version 08/06/2015
****************************************************************************/
#ifndef _jlp_kepler_
#define _jlp_kepler_

// Routines contained in "jlp_orbit1_kepler.cpp":
  int Kepler_Init_External(double *oelmts0, const int orbit_type0);
  int Solve_Kepler(const double epoch_o, double *E_anomaly,
                   double *V_anomaly);
  int Kepler_Ephemerid1(const double epoch_o, double *theta_c, double *rho_c,
                        double *E_anom_c);
  int Kepler_Ephemerid2(const double epoch_o, double *rv1_c,
                        double *rv2_c, double *E_anom_c);
  int Kepler_Dtheta(const double epoch_o, const double theta_c,
                    const double rho_c, const double E_anom,
                    double *B, const int n_B);
  int Kepler_Drho(const double epoch_o, const double theta_c,
                  const double rho_c, const double E_anom,
                  double *B, const int n_B);
  int Kepler_DV1(const double epoch_o, const double E_anom,
                 double *B, const int n_B, double *rv1_c);
  int Kepler_DV2(const double epoch_o, const double E_anom,
                 double *B, const int n_B, double *rv2_c);

// Kepler private parameters:
static double k_P, k_PMU, k_T0, k_SF, k_CF2, k_CF, k_CF3, k_EC;
static double k_A, k_WW, k_CWW, k_SWW, k_W, k_CW, k_SW, k_SI, k_CI, k_TI;
static double k_AA, k_BB, k_FF, k_GG;
static double k_K1, k_K2, k_V0;

#endif
