/****************************************************************************
* Name: gpo_rw_files.h
* 
* JLP
* Version 08/06/2015
****************************************************************************/
#ifndef _gpo_rw_files__ 
#define _gpo_rw_files__ 

#include "gpo_defs.h"         // NMEAS_MAX, NELEMENTS, N7, ...     

// In "gpo_rw_tokofiles.cpp":
int read_toko_visual_datafile(const char *in_data_file, double *epoch0,
                       double *theta0, double *rho0, double *sigma_rho0,
                       char notes_nmeas0[NMEAS_MAX][128],
                       const int nmeas_max, int *nmeas0);
int read_toko_rv_datafile(const char *in_data_file, 
                          double *epoch_rv10, double *rv10, double *sigma_rv10,
                          char notes_rv10[NMEAS_MAX][128],
                          double *epoch_rv20, double *rv20, double *sigma_rv20,
                          char notes_rv20[NMEAS_MAX][128],
                          const int nmeas_max, int *nrv10, int *nrv20);
int rd_oelements_toko(const char *orbit_infile, double *orbital_elements,
                      double *orbit_equinox0, char *object_name0, 
                      double *right_ascension0, double *declination0, 
                      int *orbit_type0);

// In "gpo_rw_datafiles.cpp":
int read_WDS_visual_datafile(const char *in_data_file, double *epoch0,
                             double *theta0, double *rho0, double *weight0,
                             char notes_nmeas[NMEAS_MAX][128],
                             const int nmeas_max, int *nmeas0);
int read_mini_visual_datafile(const char *in_data_file, double *epoch0,
                              double *theta0, double *rho0, double *weight0,
                              char notes_nmeas[NMEAS_MAX][128],
                              const int nmeas_max, int *nmeas0);
int precession_correction(double *dtheta_precess, double alpha,       
                          double delta, double epoch_o, double orbit_equinox);

// In "gpo_rw_orbitfiles.cpp":
int rd_visual_oelements_OC6(const char *orbit_infile, 
                            double *orbital_elements, 
                            double *orbit_equinox, char *object_name,  
                            double *right_ascension0, double *declination0,
                            int *orbit_type0);
int rd_visual_oelements_scardia(const char *orbit_infile, 
                                double *orbital_elements, double *orbit_equinox,
                                char *object_name, double *right_ascension0, 
                                double *declination0, int *orbit_type0);
void set_idx_free_flags(int idx_free[5][N10]);

// In "gpo_rw_coravel.cpp":
int read_carquillat_file(const char *infile, double *orbital_elements,
                         int *orbit_type0,
                         double *epoch_rv10, double *rv10, double *sigma_rv10,
                         double *weight_rv10, char notes_rv10[NMEAS_MAX][128],
                         double *epoch_rv20, double *rv20, double *sigma_rv20,
                         double *weight_rv20, char notes_rv20[NMEAS_MAX][128],
                         const int nmeas_max, int *nrv10, int *nrv20);
int read_CORAVEL_rv_datafile(const char *in_data_file,
                          double *epoch_rv10, double *rv10, double *sigma_rv10,
                          double *weight_rv10, char notes_rv10[NMEAS_MAX][128],
                          double *epoch_rv20, double *rv20, double *sigma_rv20,
                          double *weight_rv20, char notes_rv20[NMEAS_MAX][128],
                          const int nmeas_max, int *nrv10, int *nrv20);
void besselian_yr_to_julian_day(double epoch_b, double *epoch_jd);
void julian_day_to_besselian_yr(double epoch_jd, double *epoch_b);
int ConvertElementToSpectro(double *oelmts0, double *oelmts_err0,
                            double *spectro_oelmts, double *spectro_oelmts_err,
                            const int nmax_oelmts, const int orbit_type0);
int ConvertElementFromSpectro(double *spectro_oelmts, 
                              double *spectro_oelmts_err,
                              double *oelmts0, double *oelmts_err0,
                              const int nmax_oelmts, const int orbit_type0);


#endif
