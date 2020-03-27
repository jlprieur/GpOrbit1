/*************************************************************************
* gpo_rw_coravel.cpp
* To read/write spectroscopic binary data (radial velocitiesi in CORAVEL format)
*
* Contains:
*
* JLP
* Version 03/07/2015
**************************************************************************/
#include <stdio.h> 
#include <stdlib.h>    // exit() 
#include <ctype.h>     // isdigit() 
#include <string.h>    // strcpy() 
#include <math.h>      // sin() 
#include "gpo_rw_files.h"
#include "gpo_defs.h"  // ABS()

/*
#define DEBUG
*/

#define STOP_ERROR { fprintf(stderr, "gpo_rw_coravel/Fatal error reading input file\n"); fclose(fp); exit(-1);}

/* Prototypes of functions defined here:
*/
static int decode_data_line_coravel(char *in_line, double *wepoch,
                                   double *wrv1, double *wdrv1, double *wrv2,
                                   double *wdrv2, double *weight, 
                                   int *is_rv1, int *is_rv2);


/**********************************************************************
* Read RV orbit and data file in Carquillat's format
*
*
SB1: (e.g. HD32835)
288.2     45400.85  3.14      0.1       30.       15.
 51      15
45338.504 +4.9
45340.435 +3.3
48265.346 -14.8
48292.50  -14.0  0.5
48641.48  +17.0  0.5
48672.362 +32.6
48739.34  +36.1  0.5

SB2: (e.g. HD151604)
19.6986   53063.36  5.17      0.565     71.4      72.7      -13.7    
33       5
52821.44  -55.6  1.          +26.9  1.
52821.57                     +28.1  1.
52822.42                     +33.1  1.
52822.44  -61.4  1.
52824.46  -67.9  1.          +40.2  1.
52824.60  -67.8  1.
53087.62  +19.0  1.          -48.0  1.
53088.60   +7.4  1.          -37.6  1.

SB2: (e.g. HD7119)
6.761504  48945.9   0.68      0.028     42.48     46.81     -9.
1687.000  50684.4   0.00000   0.00000    5.87      0.       -2.66    2
1687.000  50460.84  0.095     0.01000    3.00      0.0      -1.66
60       5
48940.474 -21.6  1.            7.0  0.5
48966.423  19.3  1.          -39.0  0.5
48967.430 -18.2  1.            4.2  0.5
50418.429  21.9  1.          -43.4  0.5
50419.363  33.4  1.          -57.0  0.5
50419.504  34.3  1.            0.0  0.0
50420.382  10.9  1.          -33.4  0.5
50420.534   5.7  1.          -22.2  0.5
***********************************************************/
int read_carquillat_file(const char *infile, double *orbital_elements, 
                         int *orbit_type0,
                         double *epoch_rv10, double *rv10, double *sigma_rv10, 
                         double *weight_rv10, char notes_rv10[NMEAS_MAX][128],
                         double *epoch_rv20, double *rv20, double *sigma_rv20,
                         double *weight_rv20, char notes_rv20[NMEAS_MAX][128],
                         const int nmeas_max, int *nrv10, int *nrv20)
{
FILE *fp_in;
float ww0, ww1, ww2, ww3, ww4, ww5, ww6, ww7, ww8; 
double sum_nval, bessel_year;
double omega_peri, k1, k2, v0;
double e_eccent, T_periastron, Period;
int nval, i, iline, is_sb2 = 0;
char in_line[80], *pc;

*orbit_type0 = 0;
*nrv10 = 0;
*nrv20 = 0;
k2 = 0;

//*********************************************************************
/* 1. Read all lines to determine if is an SB1 or an SB2:: */
//*********************************************************************

/****************** Open file for the first time *******************/
if((fp_in = fopen(infile, "r")) == NULL) {
  fprintf(stderr, "Fatal error opening input file >%s<\n", infile);
  exit(-1);
  }

sum_nval = 0.;
iline = 0;
  while(!feof(fp_in)) {
    if(fgets(in_line, 80, fp_in)) {
      if(in_line[0] != '%' && in_line[0] != '#') {
        nval = sscanf(in_line, "%f %f %f %f %f", &ww0, &ww1, &ww2, &ww3, &ww4);
        if(nval > 1) {
          sum_nval += nval;
          iline++;
          }
       }
    }
  }
// Mean nval gives an indication on SB1/SB2 cases:
if(iline == 0) {
  fprintf(stderr, "read_carquillat_file/Error: empty file iline=%d\n", iline); 
  return(-1);
  }
sum_nval /= (double)iline;
if(sum_nval > 3.) {
  is_sb2 = 1;
#ifdef DEBUG
  printf("read_carquillat_file/SB2 data detected in file\n");
#endif
  } else {
  is_sb2 = 0;
#ifdef DEBUG
  printf("read_carquillat_file/SB1 data detected in file\n");
#endif
  }

fclose(fp_in);

/* Format of the preliminary orbital parameters:
C Period (days)
C Periastron passage or passage to the ascending node (if e=0)
C omega angle (argument of periatron in radians) or 0 if e=0
C eccentricity
C k1 (semi amplitude of the rv curve of the primary in km/s)
C V0 radial velocity of the system gravity center in km/s
C Optionnal parameters: 
C X1=0: free period
C X1=1: fixed periode fixe
C Y1=0: free eccentricity 
C Y1=1: fixed eccentricity
C Y1=2: circular orbit
C Example of SB1:
C 14.2081   48682.860 0.        0.        42.5      2.8      1     2
C Number of radial velocities (in this file)
C Nomber of iterations (that were needed for convergence)
C 42       5
C Julian date (-2000000), radial velocity of the primary, poids de la mesure
C vitesse radiale de la secondaire, poids de la mesure
C
C For SB1:
C P T0 w e K1 V0 optionnal parameters
C Fortran format for SB1:
1001  FORMAT(F10.6,F10.3,2F10.5,F10.2,F6.2,1X,2I1)
      READ(LU_IN,1001) PP,TZERO,OMEGA,EE,K1,VZERO,X1,Y1
C Data fortran format of SB1:
      READ(LU_IN,1201)TOBS(I),VIT(I),POI(I)
1201  FORMAT(F9.3,F7.1,F6.2)
C
C For SB2:
C P T0 w e K1 K2 V0 optionnal parameters
C Fortran format for SB2:
1000  FORMAT(F10.7,F10.3,2F10.4,2F10.2,F6.2,2X,2I1)
      READ(LU_IN,1000) PP,TZERO,OMEGA,EE,K1,K2,VZERO,X1,Y1
C Data fortran format of SB2:
1200  FORMAT((F9.3,2(F7.1,F6.2,7X)))
      READ(LU_IN,1200)TOBS(I),VIT1(I),POI1(I),VIT2(I),POI2(I)
*/

//*********************************************************************
// 2. Read orbital elements
//*********************************************************************
/****************** Open file for the second time *******************/
if((fp_in = fopen(infile, "r")) == NULL) {
  fprintf(stderr, "Fatal error opening input file >%s<\n", infile);
  exit(-1);
  }

if(fgets(in_line, 80, fp_in)) {
    if(in_line[0] != '%' && in_line[0] != '#') {
      nval = sscanf(in_line, "%f %f %f %f %f %f %f %f %f", 
                    &ww0, &ww1, &ww2, &ww3, &ww4, &ww5, &ww6, &ww7, &ww8);
      if((is_sb2 && nval < 7) || (nval < 6)) {
        fprintf(stderr,"read_carquillat_file/Error: bad format in first line\n");
        return(-1);
        }
      Period = ww0;
// Conversion of the epoch from Julian days to Besselian years:
      julian_day_to_besselian_yr(ww1, &bessel_year);
      T_periastron = bessel_year;
      omega_peri = ww2 * 180. / PI;
      e_eccent = ww3;
      k1 = ww4;
      if(is_sb2) {
        k2 = ww5;
        v0 = ww6;
        } else {
        v0 = ww5;
        }
     }
} else {
  fprintf(stderr,"read_carquillat_file/Error reading first uncommented line\n");
  return(-1);
  }

// Orbital elements: P T0 e a W w i K1 K2 V0
/*****************************************************************************
* Set orbital elements
* {'P', 'T', 'e', 'a', 'W', 'w', 'i', 'K1', 'K2', 'V0'};
* orbit_type = 1 visual & spectro SB1 (9 elements idx1= 1 1 1 1 1 1 1 1 0 1)
* orbit_type = 2 visual & spectro SB2 (10 elements idx1= 1 1 1 1 1 1 1 1 1 1)
* orbit_type = 3 visual (7 elements: idx1= 1 1 1 1 1 1 1 0 0 0)
* orbit_type = 4 spectroscopic SB1 (6 elements: idx1= 1 1 1 0 0 1 0 1 0 1)
* orbit_type = 5 spectroscopic SB2 (7 elements: idx1= 1 1 1 0 0 1 0 1 1 1)
*****************************************************************************/

// Conversion fom days to years:
orbital_elements[0] = Period / 365.242198781; 
orbital_elements[1] = T_periastron;
orbital_elements[2] = e_eccent; 
orbital_elements[3] = 0.; 
orbital_elements[4] = 0.; 
orbital_elements[5] = omega_peri; 
orbital_elements[6] = 0.; 
orbital_elements[7] = k1; 
orbital_elements[8] = k2; 
orbital_elements[9] = v0; 

if(k2 == 0) *orbit_type0 = 4;
 else *orbit_type0 = 5;

#ifdef DEBUG
for(i = 0; i < N10; i++) {
  printf("i=%d elmnt=%f\n", i, orbital_elements[i]);
  }
#endif

//*********************************************************************
// 2. Read radial velocities 
//*********************************************************************
*nrv10 = 0.;
*nrv20 = 0.;
while(!feof(fp_in)) {
  if(fgets(in_line, 80, fp_in)) {
    if(in_line[0] != '%' && in_line[0] != '#') {
      nval = sscanf(in_line, "%f", &ww0);
// Confirm that it is a radial velocity measurement with the Julian date:
     if(ww0 > 40000. && nval == 1) {
// Radial velocity:
// Conversion of the epoch from Julian days to Besselian years:
       julian_day_to_besselian_yr(ww0, &bessel_year);
/*
C Data fortran format of SB1:
      READ(LU_IN,1201)TOBS(I),VIT(I),POI(I)
1201  FORMAT(F9.3,F7.1,F6.2)
C Data fortran format of SB2:
1200  FORMAT((F9.3,2(F7.1,F6.2,7X)))
      READ(LU_IN,1200)TOBS(I),VIT1(I),POI1(I),VIT2(I),POI2(I)
*/

// Load rv1 and weight1 for SB1 and SB2:
       nval = sscanf(in_line, "%f %f %f", &ww0, &ww1, &ww2); 
       if(nval > 1) {
         epoch_rv10[*nrv10] = bessel_year;
         rv10[*nrv10] = ww1;
         if(nval == 3) {
           weight_rv10[*nrv10] = ww2;
           } else {
           weight_rv10[*nrv10] = 1.0;
           }
// Assume error of 0.3 km/s for Coravel
         sigma_rv10[*nrv10] = 0.3;
// Fill notes with empty string:
         strcpy(notes_rv10[*nrv10], "");
// Cancel this measurement if weight was set to zero in input file:
         if(weight_rv10[*nrv10] > 0) (*nrv10)++;
         } // EOF nval > 1

// Load rv2 and weight2 for SB2:
      if(is_sb2) {
         nval = sscanf(&in_line[28], "%f %f", &ww0, &ww1); 
         if(nval > 0) {
         epoch_rv20[*nrv20] = bessel_year;
         rv20[*nrv20] = ww0;
         if(nval == 2) {
           weight_rv20[*nrv20] = ww1;
           } else {
           weight_rv20[*nrv20] = 1.0;
           }
// Assume error of 0.3 km/s for Coravel
         sigma_rv20[*nrv20] = 0.3;
// Fill notes with empty string:
         strcpy(notes_rv20[*nrv20], "");
// Cancel this measurement if weight was set to zero in input file:
         if(weight_rv20[*nrv20] > 0) (*nrv20)++;
         } // EOF nval > 1
       }
     } // EOF nval == 1
   } // EOF if line is not commented
  } // EOF if fgets
} // EOF while !feof()

fclose(fp_in);

return(0);
}

/**********************************************************************
* Read radial-velocity data file in CORAVEL format
*
* INPUT:
* orbit_infile: name of the file
*
* OUTPUT:
* epoch_rv10, rv10, sigma_rv10, notes_rv10
* epoch_rv20, rv20, sigma_rv20, notes_rv20
* nrv10, nrv20
*
**********************************************************************/
int read_CORAVEL_rv_datafile(const char *in_data_file, 
                          double *epoch_rv10, double *rv10, double *sigma_rv10, 
                          double *weight_rv10, char notes_rv10[NMEAS_MAX][128],
                          double *epoch_rv20, double *rv20, double *sigma_rv20,
                          double *weight_rv20, char notes_rv20[NMEAS_MAX][128],
                          const int nmeas_max, int *nrv10, int *nrv20)
{
double wepoch, wrv1, wdrv1, wrv2, wdrv2, weight;
char in_line[80], buffer[80], *pc;
FILE *fp_in;
int status, i, is_rv1, is_rv2, is_theta_rho;
char notes[128];

if((fp_in = fopen(in_data_file, "r")) == NULL) {
  fprintf(stderr, "read_coravel_rv_datafile/Fatal error opening input file >%s<\n",
          in_data_file);
  return(-1);
  }

*nrv10 = 0;
*nrv20 = 0;
while(!feof(fp_in)) {
  if(fgets(buffer, 80, fp_in)) {

// Remove headings blanks if present:
    pc = buffer;
    if(*pc == ' ')  {
      while(*pc == ' ') pc++;
      }

// Check if first character is a digit: */
    if(*pc && isdigit(*pc))  {

      strcpy(in_line, pc);
      status = decode_data_line_coravel(in_line, &wepoch, &wrv1, &wdrv1, &wrv2,
                                        &wdrv2, &weight, &is_rv1, &is_rv2);
      if(status != 0){
        fprintf(stderr, "read_coravel_rv_datafile/Error reading line >%s<\n ichar=%d Last line or End of File ?)\n",
                in_line, (int)in_line[0]);
      } else {

// Possibility of rv1 and rv2 in same line:
        if(is_rv1) {
          epoch_rv10[*nrv10] = wepoch;
          rv10[*nrv10] = wrv1;
          sigma_rv10[*nrv10] = wdrv1;
          weight_rv10[*nrv10] = weight;
          strcpy(notes_rv10[*nrv10], notes);
          if(*nrv10 >= nmeas_max) {
           fprintf(stderr, "read_coravel_rv_datafile/Fatal error: nrv10=%d reaches limit storage\n", *nrv10);
           exit(1);
           }
          (*nrv10)++;
          }
        if(is_rv2) {
          epoch_rv20[*nrv20] = wepoch;
          rv20[*nrv20] = wrv2;
          sigma_rv20[*nrv20] = wdrv2;
          weight_rv20[*nrv20] = weight;
          strcpy(notes_rv20[*nrv20], notes);
          if(*nrv20 >= nmeas_max) {
           fprintf(stderr, "read_coravel_rv_datafile/Fatal error: nrv20=%d reaches limit storage\n", *nrv20);
           exit(1);
           }
          (*nrv20)++;
          }
     } // EOF status == 0
   } /* EOF in_line != % */
 } /* EOF fgets */
} /* EOF while !feof*/

fclose(fp_in);

#ifdef DEBUG
printf("read_coravel_rv_datafile/measurements successfully read : nrv10=%d nrv20=%d\n",
       *nrv10, *nrv20);
#endif

return(0);
}
/***************************************************************************
* Read a measuremement line from data file in CORAVEL format
*
***************************************************************************/
static int decode_data_line_coravel(char *in_line, double *wepoch,
                                   double *wrv1, double *wdrv1, double *wrv2,
                                   double *wdrv2, double *weight, 
                                   int *is_rv1, int *is_rv2)
{
int i, nval, is_SB2, status = -1;
double ww0, ww1, ww2, ww3, ww4, ww5, bessel_year;

/**************** Example: ***************************************
Julian_date rv1 sigma_rv1 [weight]

48940.406 -24.3
48966.368   6.0
48967.343  28.5
48969.337 -6.3
48970.250 -31.2

or:
Julian_date rv1 sigma_rv1 rv2 sigma_rv2 [weight]

46632.364   -6.37 0.44 176 12.32
46632.364   -0.20 0.80 176  6.78
47012.343   -4.26 0.31 565 22.42
******************************************************************/

  *wepoch = 0.;
  *wrv1 = 0.;
  *wdrv1 = 0.;
  *wrv2 = 0.;
  *wdrv2 = 0.;
  *is_rv1 = 0;
  *is_rv2 = 0;
  *weight = 1.;

// Decode epoch:
  nval = sscanf(in_line, "%lf", &ww0);
  if(nval != 1) {
    fprintf(stderr, "decode_data_line_coravel/Error reading >%s<\n", in_line);
    return(1);
    }

// Confirm that it is a radial velocity measurement with the Julian date:
  if(ww0 > 40000.) {
// Radial velocity:
// Conversion of the epoch from Julian days to Besselian years:
    julian_day_to_besselian_yr(ww0, &bessel_year);
    *wepoch = bessel_year;

/*
* The  radial velocity  data  can be  in  two forms:
*
* epoch RV1 err1 for SB1 (one radial velocity)
* or:
* epoch RV1, err1, RV2, Err2  for SB2 (both velocities).
*/
    nval = sscanf(in_line,"%lf %lf %lf %lf %lf %f", 
                  &ww0, &ww1, &ww2, &ww3, &ww4, &ww5);
    switch(nval) {
// epoch rv1 err1 rv2 err2 weight
     case 6:
       *wrv1 = ww1;
       *wdrv1 = ww2;
       *wrv2 = ww3;
       *wdrv2 = ww4;
       *weight = ww5;
       *is_rv1 = 1;
       *is_rv2 = 1;
       status = 0;
       break;
// epoch rv1 err1 rv2 err2
     case 5:
       *wrv1 = ww1;
       *wdrv1 = ww2;
       *wrv2 = ww3;
       *wdrv2 = ww4;
       *is_rv1 = 1;
       *is_rv2 = 1;
       status = 0;
       break;
// epoch rv1 rv2 weight 
     case 4:
       *wrv1 = ww1;
       *wdrv1 = ww2;
       *weight = ww3;
       *is_rv1 = 1;
       status = 0;
       break;
// epoch rv1 rv2
// or:
// epoch rv1 weight 
     case 3:
       *wrv1 = ww1;
       *is_rv1 = 1;
// If 3rd value is small, assume it is the error 
// otherwise, assume it is rv2:
       if(ABS(ww2) < 2.) { 
          *wdrv1 = ww2;
         } else {
          *wrv2 = ww2;
          *is_rv2 = 1;
         }
       status = 0;
       break;
// epoch rv1 
     case 2:
       *wrv1 = ww1;
       *is_rv1 = 1;
       status = 0;
       break;
     default: 
       fprintf(stderr, "Error reading line >%s< nval=%d\n", in_line, nval);
       status = -1;
       break;
     } // EOF switch(nval)
  } // EOF ww0 > 40000

return(status);
}
/***********************************************************************
* Conversion from Besselian years to Julian days
***********************************************************************/
void besselian_yr_to_julian_day(double epoch_b, double *epoch_jd)
{
*epoch_jd = 365.242198781 *(epoch_b - 1900.) + 15020.31352;
}

/***********************************************************************
* Conversion from Julian days to Besselian years
***********************************************************************/
void julian_day_to_besselian_yr(double epoch_jd, double *epoch_b)
{
*epoch_b = 1900.0 + (epoch_jd - 15020.31352) / 365.242198781;
}
/********************************************************************
* Conversion to Julian days, days etc, according to orbit_type 
* {'P', 'T', 'e', 'a', 'W', 'w', 'i', 'K1', 'K2', 'V0'};
********************************************************************/
int ConvertElementToSpectro(double *oelmts0, double *oelmts_err0,
                            double *spectro_oelmts, double *spectro_oelmts_err,
                            const int nmax_oelmts, const int orbit_type0)
{
double JDay, bessel_year;
int i;

 if(nmax_oelmts != N10) {
   fprintf(stderr, "ConvertElementToSpectro/Error: nmax_elmts=%d\n",
           nmax_oelmts);
   return(-1);
   }

for(i = 0; i < N10; i++) {
   spectro_oelmts[i] = oelmts0[i];
   spectro_oelmts_err[i] = oelmts_err0[i];
   }

// Conversion if needed:
if(orbit_type0 != 3) {

  spectro_oelmts[0] = oelmts0[0] * YEAR_TO_DAYS;
  spectro_oelmts_err[0] = oelmts_err0[0] * YEAR_TO_DAYS;

// Conversion of the epoch from Julian days to Besselian years:
// julian_day_to_besselian_yr(JDay, &bessel_year);
// Conversion of the epoch from Besselian years to Julian days:
  bessel_year = oelmts0[1];
  besselian_yr_to_julian_day(bessel_year, &JDay);
  spectro_oelmts[1] = JDay;
  spectro_oelmts_err[1] = oelmts_err0[1] * YEAR_TO_DAYS;
  }

return(0);
}
/********************************************************************
* Conversion from Julian days, days etc, according to orbit_type 
* {'P', 'T', 'e', 'a', 'W', 'w', 'i', 'K1', 'K2', 'V0'};
********************************************************************/
int ConvertElementFromSpectro(double *spectro_oelmts, 
                              double *spectro_oelmts_err,
                              double *oelmts0, double *oelmts_err0,
                              const int nmax_oelmts, const int orbit_type0)
{
double JDay, bessel_year;
int i;

 if(nmax_oelmts != N10) {
   fprintf(stderr, "ConvertElementToSpectro/Error: nmax_elmts=%d\n",
           nmax_oelmts);
   return(-1);
   }

for(i = 0; i < N10; i++) {
   oelmts0[i] = spectro_oelmts[i];
   oelmts_err0[i] = spectro_oelmts_err[i];
   }

// Conversion if needed:
if(orbit_type0 != 3) {

  oelmts0[0] = spectro_oelmts[0] / YEAR_TO_DAYS;
  oelmts_err0[0] = spectro_oelmts_err[0] / YEAR_TO_DAYS;

// Conversion of the epoch from Julian days to Besselian years:
  JDay = spectro_oelmts[1];
  julian_day_to_besselian_yr(JDay, &bessel_year);
  oelmts0[1] = bessel_year;
  oelmts_err0[1] = spectro_oelmts_err[1] / YEAR_TO_DAYS;
  }

return(0);
}
