/****************************************************************************
* Name: jlp_orbit1.h
* For visual orbits only
* 
* JLP
* Version 08/06/2015
****************************************************************************/
#ifndef _jlp_orbit1_
#define _jlp_orbit1_

#include <stdio.h>
// For compilers that support precompilation, includes "wx/wx.h".
#include "wx/wxprec.h"

#ifndef WX_PRECOMP
    #include "wx/wx.h"
#endif

#include "wx/tglbtn.h"
#include "wx/bookctrl.h"
#include "wx/imaglist.h"
#include "wx/cshelp.h"

#if wxUSE_TOOLTIPS
    #include "wx/tooltip.h"
#endif

#ifndef ABS
#define ABS(a) ((a) < 0.0  ? (-(a)) : (a))
#endif

#ifndef PI
#define PI 3.14159265
#endif

#ifndef SQUARE
#define SQUARE(a) ((a) * (a))
#endif

#include "gpo_defs.h" // NELEMENTS, NMEAS_MAX

//----------------------------------------------------------------------
// class definitions
//----------------------------------------------------------------------

class JLP_Orbit1
{
public:

// Constructors in "jlp_orbit1.cpp":
  JLP_Orbit1();

  JLP_Orbit1(double *epoch0, double *theta0, double *rho0, double *sigma_rho0, 
             int nmeas0, double *oelements0, double *oelements_err0, 
             const int max_noelements0, int orbit_type0);

// Destructor:
  ~JLP_Orbit1(){
    return;
    };

// Routines contained in "jlp_orbit1.cpp":
  void ClearMeasurements();
  int LoadMeasurements(double *epoch0, double *theta0, double *rho0,
                       double *sigma_rho0, const int nmeas0);
  int LoadRadialVelocities(double *epoch_rv1_0, double *epoch_rv2_0,
                           double *rv1_0, double *rv2_0,
                           double *sigma_rv1_0, double *sigma_rv2_0,
                           const int nrv1_0, const int nrv2_0);
  int set_orbit_type_from_data();
  int SetOrbitalElements(double *oelements0, double *oelements_err0,
                         const int max_noelements0, const int orbit_type0);
  int GetOrbitalElements(double *oelements0, double *oelements_err0,
                         const int max_nelmts, int *orbit_type0);

  int ComputeEphemerid(float epoch_o, float *theta_c, float *rho_c);
  int Kepler_Init_Internal();

// Routines contained in "jlp_orbit1_lsqfit1.cpp":
  int LSquaresFit1(double *oelements2, double *oelements2_err,
                   const int max_noelements2, 
                   double *mean_theta_resid, double *mean_sigma_theta_resid,
                   double *mean_rho_resid, double *mean_sigma_rho_resid, 
                   double *mean_sigma_rv1_resid, double *mean_sigma_rv2_resid, 
                   double *chisq2, wxString &results_str);

  int LSQ1_covar_from_visual_meas(double alpha[NELEMENTS][NELEMENTS], 
                                  double *beta, int *oel_idx2,
                                  const int noel2, double *chisq2,  
                                  double *mean_theta_resid,
                                  double *mean_sigma_theta_resid,
                                  double *mean_rho_resid,
                                  double *mean_sigma_rho_resid);
  int LSQ1_covar_from_rv_meas(double alpha[NELEMENTS][NELEMENTS], 
                              double beta[NELEMENTS], int *oel_idx2,
                              double oelmnt2[NELEMENTS],
                              const int noel2, double *chisq2,  
                              double *mean_rv1_resid, double *mean_rv2_resid);

// Routines contained in "jlp_orbit1_gauss.cpp":
  int LSQ1_GaussInversion(double alpha[NELEMENTS][NELEMENTS],
                          const int n_free_oelements, double beta[NELEMENTS]);

// Routines contained in "jlp_orbit1_lsqfit2.cpp":
  int LSquaresFit2(double *oelements2, double *oelements2_err,
                   const int max_noelements2, 
                   double *mean_theta_resid, double *mean_sigma_theta_resid,
                   double *mean_rho_resid, double *mean_sigma_rho_resid, 
                   double *mean_sigma_rv1_resid, double *mean_sigma_rv2_resid, 
                   double *chisq2, wxString &results_str);

// Routine contained in "jlp_orbit1_cmpfit.cpp":
  int LSQ2_fit_with_cmpfit(double *elmt, double *elmt_err, int *oel_idx, 
                           int n_free_elmt, 
                           double *mean_theta_resid,
                           double *mean_sigma_theta_resid,
                           double *mean_rho_resid,
                           double *mean_sigma_rho_resid,
                           double *mean_sigma_rv1_resid,
                           double *mean_sigma_rv2_resid,
                           double *chisq2, wxString& Results_str);
  int LSQ2_compute_residuals(double *elmt, int *oel_idx, int n_free_elmt,
                             double *mean_theta_resid,
                             double *mean_sigma_theta_resid,
                             double *mean_rho_resid,
                             double *mean_sigma_rho_resid,
                             double *mean_sigma_rv1_resid,
                             double *mean_sigma_rv2_resid,
                             double *chisq2);
  int LSQ2_init_static_arrays(double *elmt, int *oel_idx, int n_free_elmt);

// Accessors:

// nmeas1: number of (rho, theta) measurements
// nrv1: number of radial velocity measurements of the primary
// nrv2: number of radial velocity measurements of the secondary
   int Get_nmeas1() {return(nmeas1);}
   int Get_nrv1() {return(nrv1);}
   int Get_nrv2() {return(nrv2);}
   int GetOrbitType() {return(orbit_type1);}
   int GetOrbitTypeFromElements() {return(orbit_type_from_elements1);}
// 'P', 'T', 'e', 'a', 'W', 'w', 'i', 'K1', 'K2', 'V0'
   double GetPeriod() {return(oelements1[0]);}
   double GetT0() {return(oelements1[1]);}
   double GetEccent() {return(oelements1[2]);}
   double GetSemiAxis() {return(oelements1[3]);}
   double GetOmega() {return(oelements1[4]);}
   double Getomega() {return(oelements1[5]);}
   double GetIncl() {return(oelements1[6]);}
   double GetK1() {return(oelements1[7]);}
   double GetK2() {return(oelements1[8]);}
   double GetV0() {return(oelements1[9]);}

/***********************************************************************
* Get values of fixed flags (used by the fitting procedure) 
***********************************************************************/
   int GetFixedFlags(int *fflags0, const int max_nelmts0)
     { int i;
      
      if((initialized != 1234) || (max_nelmts0 != N10)) { 
        fprintf(stderr, " GetFixedFlags/Error: nelmts0=%d initialized=%d\n",
                max_nelmts0, initialized);
        return(-1);
        }

      for(i = 0; i < N10; i++) fflags0[i] = fixed_flags1[i];
            
     return(0);
     }

/***********************************************************************
* Set values of fixed flags (used by the fitting procedure) 
***********************************************************************/
   int SetFixedFlags(int *fflags0, const int max_nelmts0)
     { int i;
      
      if((initialized != 1234) || (max_nelmts0 != N10)) { 
        fprintf(stderr, " SetFixedFlags/Error: max_nelmts0=%d initialized=%d\n",
                max_nelmts0, initialized);
        return(-1);
        }

      for(i = 0; i < N10; i++) fixed_flags1[i] = fflags0[i];

     return(0);
     }

private:
  int initialized;
// Orbit type: 1=visual&spectro... 
  int orbit_type1, orbit_type_from_elements1;
  double oelements1[NELEMENTS], oelements_err1[NELEMENTS], chisq2;
  double epoch1[NMEAS_MAX], theta1[NMEAS_MAX], rho1[NMEAS_MAX];
  double sigma_rho1[NMEAS_MAX];
  double epoch_rv1[NMEAS_MAX], rv1[NMEAS_MAX]; 
  double sigma_rv1[NMEAS_MAX], sigma_rv2[NMEAS_MAX];
  double epoch_rv2[NMEAS_MAX], rv2[NMEAS_MAX]; 

// nmeas1: number of (rho, theta) measurements read from input file
// nrv1, nrv2: number of radial velocity measurements read from input file
  int nmeas1, nrv1, nrv2;

// fixed_flags1: flags of the elements whose value will remain fixed in the fit 
  int fixed_flags1[NELEMENTS];
  int idx1[NELEMENTS];

};
#endif
