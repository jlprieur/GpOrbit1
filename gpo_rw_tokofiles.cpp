/*************************************************************************
* gpo_rw_tokofiles.cpp
* from "orbit/orbit_plot_utils.c"
* To read/write orbit data files in Tokovinin's format
*
* Contains:
* read_toko_rv_datafile()
* read_toko_visual_datafile()
* rd_oelements_toko()
* ComputeResiduals()
*
* JLP
* Version 13/06/2015
**************************************************************************/
#include <stdio.h> 
#include <stdlib.h>    //exit() 
#include <ctype.h>     //isalpha()

#include "gpo_frame.h"
#include "gpo_rw_files.h" // set_idx_free_flags()
#include "jlp_kepler.h"

/* 
#define DEBUG
*/

static int decode_oelement_line_toko1(char *in_line, double *orbital_elements, 
                                    int nelmts); 
static int read_int_keyword(char *in_line, char *keyword, int *ivalue);
static int read_float_keyword(char *in_line, char *keyword, float *fvalue);
static int extract_end_notes_from_dataline(char *in_line, char *notes);
static int find_component_in_rv_notes(char *notes, int *is_rv1, int *is_rv2,
                                      int *is_SB2);
static int orbit_output_curves(float *xplot, float *yplot, int npts_max, 
                               int *npts, int ncurves, char *outfile);
static int decode_data_line_toko1(char *in_line, double *wepoch, 
                                  double *wtheta, double *wrho, double *wdrho,
                                  double *rv1, double *drv1, double *rv2, 
                                  double *drv2, char *notes,
                                  int *is_rv1, int *is_rv2, int *is_theta_rho);

#define UNDETERMINED -12345.

/******************************************************
* Format 1:  Tokovinin 1991 (without radial velocities)
*******************************************************
Object: COU14i GII2010                  
R.A.:   21.502                           
Decl:   17.18                           
N_VR:   0                           
N_OB:   47                           
P:        26.132                  
T0:     1963.887                  
e:         0.239                  
a:         0.3664                  
W:       231.80                  
w:       252.08                  
i:        70.30                  
K1:    *                *         
K2:    *                *          
V0:    *                *         
1976.8594  275.9  0.200  0.003 n1 Mca4m
1977.4819  286.1  0.179  0.003 n1 Mca4m
1977.9189  294.2  0.171  0.003 n1 Mca2m
*******************************************************
* Format 2:  Tokovinin 2015 (without radial velocities) 
*******************************************************
Object:   FIN379   
RA:       2.4412
Dec:     -25.30
P          6.703
T       2008.8426
e           0.506 
a           0.0996
W           4.6 
w          368.9 
i           42.4 
*K1         17.76
*K2           0.0              0.00
*V0         -22.27
      1963.0500    142.70   0.114000  0.05   I1 Fin1963a 26J  0
      1964.0341    171.40   0.140000  0.05   I1 Fin1964a 26J  1
      1965.0460    185.40   0.138000  0.05   I1 Fin1965a 26J  1
*******************************************************
* Format 3:  Tokovinin 2015 (with radial velocities) 
*******************************************************
Object: Mlr 224 = GL 765.2
R.A.:   19.4040
Dec:   76.1812
P       11.769
T       1993.513
e       0.224
a       0.225
W       106.34
w        89.40
i        82.56
K1       7.54
K2       6.96
V0      -3.91
45533.4644  -10.69    0.51 Va COR
45543.4416  -11.25    0.51 Va COR
49547.4899    0.12    0.63 Vb COR
49592.4857    0.94    0.57 Vb COR
C RVM data
46632.364   -6.37 0.44 176 12.32  3.18 Va
46632.364   -0.20 0.80 176  6.78  3.18 Vb
47012.343   -4.26 0.31 565 22.42  2.48 Va+b
C Visual and interferometric data
1971.6   275.6 0.21 0.04 I1 M1
1971.57  275.6 0.17 0.04 I1 M1
*******************************************************/
/*************************************************************************
* read_toko_rv_datafile (Tokovinin's format) 
*
* INPUT:
* in_data_file: name of the file
*
* OUTPUT:
* epoch_rv10, rv10, sigma_rv10, notes_rv10
* epoch_rv20, rv20, sigma_rv20, notes_rv20
* nrv10, nrv20
*************************************************************************/
int read_toko_rv_datafile(const char *in_data_file, 
                          double *epoch_rv10, double *rv10, double *sigma_rv10,
                          char notes_rv10[NMEAS_MAX][128],
                          double *epoch_rv20, double *rv20, double *sigma_rv20,
                          char notes_rv20[NMEAS_MAX][128],
                          const int nmeas_max, int *nrv10, int *nrv20)
{
double wepoch, wtheta, wrho, wdrho, wrv1, wdrv1, wrv2, wdrv2; 
char in_line[80], buffer[80], *pc;
FILE *fp_in;
int status, i, is_rv1, is_rv2, is_theta_rho;
char notes[128];

if((fp_in = fopen(in_data_file, "r")) == NULL) {
  fprintf(stderr, "read_toko_rv_datafile/Fatal error opening input file >%s<\n",
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
      status = decode_data_line_toko1(in_line, &wepoch, &wtheta, &wrho, 
                                      &wdrho, &wrv1, &wdrv1, &wrv2, &wdrv2,
                                      notes, &is_rv1, &is_rv2, &is_theta_rho); 
      if(status != 0){
        fprintf(stderr, "read_toko_rv_datafile/Error reading line >%s<\n ichar=%d Last line or End of File ?)\n",
                in_line, (int)in_line[0]);
      } else {

// Possibility of rv1 and rv2 in same line:
        if(is_rv1) {
          epoch_rv10[*nrv10] = wepoch;
          rv10[*nrv10] = wrv1;
          sigma_rv10[*nrv10] = wdrv1;
          strcpy(notes_rv10[*nrv10], notes);
          if(*nrv10 >= nmeas_max) {
           fprintf(stderr, "read_toko_rv_datafile/Fatal error: nrv10=%d reaches limit storage\n", *nrv10);
           exit(1);
           }
          (*nrv10)++;
          }
        if(is_rv2) {
          epoch_rv20[*nrv20] = wepoch;
          rv20[*nrv20] = wrv2;
          sigma_rv20[*nrv20] = wdrv2;
          strcpy(notes_rv20[*nrv20], notes);
          if(*nrv20 >= nmeas_max) {
           fprintf(stderr, "read_toko_rv_datafile/Fatal error: nrv20=%d reaches limit storage\n", *nrv20);
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
printf("read_toko_rv_datafile/measurements successfully read : nrv10=%d nrv20=%d\n", 
       *nrv10, *nrv20);
#endif

return(0);
}
/*************************************************************************
* read_toko_visual_datafile (Tokovinin's format) 
*
* INPUT:
* in_data_file: name of the file
*
* OUTPUT:
* epoch0, theta0, rho0, sigma_rho0: measurement values 
* notes_meas0: notes read in measurement line
* nmeas0: number of measurements that could be read in input file 
*************************************************************************/
int read_toko_visual_datafile(const char *in_data_file, double *epoch0,
                       double *theta0, double *rho0, double *sigma_rho0,
                       char notes_meas0[NMEAS_MAX][128],
                       const int nmeas_max, int *nmeas0)
{
double wepoch, wtheta, wrho, wdrho, wrv1, wdrv1, wrv2, wdrv2; 
char in_line[80], buffer[80], *pc;
FILE *fp_in;
int status, i, is_rv1, is_rv2, is_theta_rho;
char notes[128];

if((fp_in = fopen(in_data_file, "r")) == NULL) {
  fprintf(stderr, "read_toko_visual_datafile/Fatal error opening input file >%s<\n",
          in_data_file);
  return(-1);
  }

*nmeas0 = 0;
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

      status = decode_data_line_toko1(in_line, &wepoch, &wtheta, &wrho, 
                                      &wdrho, &wrv1, &wdrv1, &wrv2, &wdrv2,
                                      notes, &is_rv1, &is_rv2, &is_theta_rho); 
      if(status != 0){
        fprintf(stderr, "read_toko_visual_datafile/Error reading line >%s<\n ichar=%d Last line or End of File ?)\n",
                in_line, (int)in_line[0]);
      } else {
      if(is_theta_rho) {
        epoch0[*nmeas0] = wepoch; 
        theta0[*nmeas0] = wtheta; 
        rho0[*nmeas0] = wrho; 
        sigma_rho0[*nmeas0] = wdrho; 
        strcpy(notes_meas0[*nmeas0], notes);
        if(*nmeas0 >= nmeas_max) {
         fprintf(stderr, "read_toko_visual_datafile/Fatal error: *nmeas0=%d reaches limit storage\n %s", *nmeas0);
         exit(1);
         }
        (*nmeas0)++;
       } // EOF is_theta_rho 
     } // EOF status == 0
   } /* EOF in_line != % */
 } /* EOF fgets */
} /* EOF while !feof*/

fclose(fp_in);

#ifdef DEBUG
printf("read_toko_visual_datafile/measurements successfully read : nmeas0=%d\n", 
       *nmeas0);
#endif

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

/*************************************************************************
* rd_oelements_toko (Tokovinin's format) 
*
* Orbital elements: P T0 e a W w i K1 K2 V0
*
* INPUT:
* orbit_infile: name of the file
*
* OUTPUT:
* orbital_elements: orbital elements
* orbit_type0: type of orbit (visual, spectroscopic, etc)
*
*************************************************************************/
int rd_oelements_toko(const char *orbit_infile, double *oelements0,
                      double *orbit_equinox0,
                      char *object_name0, double *right_ascension0,
                      double *declination0, int *orbit_type0)
{
double oelements[N10], oelements_err[N10], oelements_err0[N10];
double bessel_year, ww;
char in_line[80], buffer[80], *pc;
FILE *fp_in;
int status;
int i, nval, nelmts, header_is_read;
char notes[128];

strcpy(object_name0, "");
*right_ascension0 = 0.;
*declination0 = 0.;
*orbit_equinox0 = 0.;

nelmts = NELEMENTS;

for(i = 0; i < nelmts; i++) oelements[i] = UNDETERMINED;

if((fp_in = fopen(orbit_infile, "r")) == NULL) {
  fprintf(stderr, "rd_oelements_toko/Fatal error opening input file >%s<\n",
          orbit_infile);
  return(-1);
  }

header_is_read = 0;
while(!feof(fp_in) && !header_is_read) {
  if(fgets(buffer, 80, fp_in)) {

// Remove headings blanks if present:
    pc = buffer;
    if(*pc == ' ')  {
      while(*pc == ' ') pc++;
      }

// Check if first character is alpha, not a comment, or '*': 
    if((isalpha(*pc) && *pc != 'C') || *pc == '*') {

     strcpy(in_line, pc);
// Special treatment of the first lines:
// to read object name, R.A. and Decl:
     if(!strncmp(in_line, "Object", 6)) {
       pc = in_line;
// Skip Object or Object:
       while(isalpha(*pc) || *pc == ':') pc++;
       while(*pc == ' ') pc++;
       strcpy(object_name0, pc);
// Cleans strings from possible "\r" at the end:
       pc = object_name0;
       while(*pc && *pc != '\r' && *pc != '\n') pc++;
       *pc = '\0';
     } else if(!strncmp(in_line, "RA", 2)
              || !strncmp(in_line, "R.A.", 4)) {
       pc = in_line;
// Skip RA RA: R.A. or R.A.:
       while(isalpha(*pc) || *pc == '.' || *pc == ':') pc++;
       while(*pc == ' ') pc++;
       nval = sscanf(pc, "%lf", &ww);
       if(nval == 1) *right_ascension0 = ww;
// Declination:
     } else if(!strncmp(in_line, "Dec", 3)) {
       pc = in_line;
// Skip Dec Dec: Decl or Decl:
       while(isalpha(*pc) || *pc == '.' || *pc == ':') pc++;
       while(*pc == ' ') pc++;
       nval = sscanf(pc, "%lf", &ww);
       if(nval == 1) *declination0 = ww;
// Equinox:
     } else if(!strncmp(in_line, "Equinox", 8)) {
       pc = in_line;
// Skip Equinox or Equinox:
       while(isalpha(*pc) || *pc == '.' || *pc == ':') pc++;
       while(*pc == ' ') pc++;
       nval = sscanf(pc, "%lf", &ww);
       if(nval == 1) *orbit_equinox0 = ww;
     } else {
      decode_oelement_line_toko1(in_line, oelements, nelmts);
     }
// Last header line is 'V0' or '*V0':
    if(in_line[0] == 'V' || in_line[1] == 'V') header_is_read = 1;
   } /* EOF in_line != % */
 } /* EOF fgets */
} /* EOF while !feof*/

fclose(fp_in);

/*****************************************************************************
* {'P', 'T', 'e', 'a', 'W', 'w', 'i', 'K1', 'K2', 'V0'};
* orbit_type = 1 visual & spectro SB1 (9 elements idx1= 1 1 1 1 1 1 1 1 0 1)
* orbit_type = 2 visual & spectro SB2 (10 elements idx1= 1 1 1 1 1 1 1 1 1 1)
* orbit_type = 3 visual (7 elements: idx2= 1 1 1 1 1 1 1 0 0 0)
* orbit_type = 4 spectroscopic SB1 (6 elements: idx3= 1 1 1 0 0 1 0 1 0 1)
* orbit_type = 5 spectroscopic SB2 (7 elements: idx3= 1 1 1 0 0 1 0 1 1 1)
*****************************************************************************/

// Check if all elements have been retrieved from file:
nval = 0;
for(i = 0; i < N10; i++) {
  if(oelements[i] != UNDETERMINED) nval++; 
  }

*orbit_type0 = 0;
if(nval == 6 || nval == 7 || nval == 9 || nval == 10) {
  if(nval == 7 && oelements[7] == UNDETERMINED) {
   *orbit_type0 = 3;
  } else if(nval == 7 && oelements[3] == UNDETERMINED) {
    if(oelements[8] != UNDETERMINED) {
     *orbit_type0 = 5;
    } else {
     *orbit_type0 = 4;
    }
  } else if( nval == 10) {
   if(oelements[8] == 0) *orbit_type0 = 1;
     else *orbit_type0 = 2;
  }
} // EOF if nval == 6 ...


// JLP2015: correction for spectroscopic orbits: 
for(i = 0; i < nelmts; i++) oelements_err[i] = 0.;
ConvertElementFromSpectro(oelements, oelements_err, 
                          oelements0, oelements_err0, 
                          N10, *orbit_type0);

#ifdef DEBUG
printf("obj=%s \n ra=%.3f \n dec=%.3f \n", 
       object_name0, right_ascension0, declination0);
printf("rd_oelements_toko/ orbit_type0=%d\n", *orbit_type0);

for(i = 0; i < N10; i++) {
   printf(" oelements0[%d] = %f \n", i, oelements0[i]);
   }
#endif

return(0);
}
/***************************************************************************
* Read a oelement line from data file in Tokovinin's format
*
***************************************************************************/
static int decode_oelement_line_toko1(char *in_line, double *orbital_elements, 
                                    int nelmts) 
{
int status, n_fkeys = NELEMENTS;
char fkeyword[NELEMENTS][8], ffkeyword[NELEMENTS][8];
int i, ivalue;
float fvalue;

if(nelmts != n_fkeys) {
  fprintf(stderr,"decode_oelement_line_toko1/Fatal error, wrong size nelmts=%d n_fkeys=%d !\n",
         nelmts, n_fkeys);
 }

// Float keywords: P, T0, e, a, W, w, i, K1, K2 and V0 
// or: P, T, e, a, o, w, i, K1, K2 and V0
strcpy(fkeyword[0], "P"); 
strcpy(fkeyword[1], "T0"); 
strcpy(fkeyword[2], "e"); 
strcpy(fkeyword[3], "a"); 
strcpy(fkeyword[4], "W"); 
strcpy(fkeyword[5], "w"); 
strcpy(fkeyword[6], "i"); 
strcpy(fkeyword[7], "K1"); 
strcpy(fkeyword[8], "K2"); 
strcpy(fkeyword[9], "V0"); 
//
strcpy(ffkeyword[0], "P"); 
strcpy(ffkeyword[1], "T"); 
strcpy(ffkeyword[2], "e"); 
strcpy(ffkeyword[3], "a"); 
strcpy(ffkeyword[4], "o"); 
strcpy(ffkeyword[5], "w"); 
strcpy(ffkeyword[6], "i"); 
strcpy(ffkeyword[7], "K1"); 
strcpy(ffkeyword[8], "K2"); 
strcpy(ffkeyword[9], "V0"); 
//

for(i = 0; i < n_fkeys; i++) {
// First try:
  status = read_float_keyword(in_line, fkeyword[i], &fvalue);
// Second try:
  if(status) {
     status = read_float_keyword(in_line, ffkeyword[i], &fvalue);
// If '*' in front of keyword, non used parameter:
     if(status == 0 && ffkeyword[i][0] == '*') {
       return(0);
       }
     }
  if(status == 0) {
    orbital_elements[i] = fvalue; 
// Return if one float keyword has been found in current line:
    return(0);
    }
}

return(-1);
}
/***********************************************************************
* Read keyword assumin syntax with column separator:
* for instance: 
* N_OB:   47                        (keyword: "N_OB")
* N_VR=45                           (keyword: "N_VR")
* NO LONGER USED !
***********************************************************************/
static int read_int_keyword(char *in_line, char *keyword, int *ivalue)
{
 int klen, status = -1, ival;
 char *pc;

 *ivalue = UNDETERMINED; 
 klen = strlen(keyword);

// Skip ':' or '='
 if(in_line[klen + 1] == ':' || in_line[klen + 1] == '=') klen++; 

// Test if line starts with keyword:
 if(strncmp(in_line,keyword, klen) == 0) { 
// Good line:
   if(in_line[klen + 1] != '\0') { 
     pc = &in_line[klen + 1]; 
     if(sscanf(pc, "%d", &ival) == 1) {*ivalue = ival; status = 0;}
   }
 }
return(status);
}
/***********************************************************************
* Read keyword assumin syntax with column separator:
* for instance: 
* P:   115.8767                    (keyword: "P")
* K1    *                *         (keyword: "K1")
***********************************************************************/
static int read_float_keyword(char *in_line, char *keyword, float *fvalue)
{
 int klen, status = -1;
 float fval;
 char *pc;

 *fvalue = UNDETERMINED; 
 klen = strlen(keyword);

// Skip ':' or '='
 if(in_line[klen + 1] == ':' || in_line[klen + 1] == '=') klen++; 

// Test if line starts with keyword:
 if(strncmp(in_line,keyword, klen) == 0) { 
// Good line:
   if(in_line[klen + 1] != '\0') { 
     pc = &in_line[klen + 1]; 
     if(sscanf(pc, "%f", &fval) == 1) {
       *fvalue = fval; 
        status = 0;
       }
   }
 }
return(status);
}
/**************************************************************************
* Routine to output one or two curves (yplot versus xplot) to file
*
* INPUT:
* xplot[ncurves*npts_max], yplot[ncurves*npts_max]: arrays to be plotted
* npts_max: maximum size of curve to be plotted
* npts[]: number of points of the curves contained in arrays xplot, yplot
* ncurves: number of curves
* outfile: output file name
*
**************************************************************************/
static int orbit_output_curves(float *xplot, float *yplot, int npts_max, 
                               int *npts, int ncurves, char *outfile)
{
FILE *fp_out;
register int i;

if((fp_out = fopen(outfile, "w")) == NULL) {
  fprintf(stderr, "orbit_output_curves/Fatal error opening output file >%s<\n",
          outfile);
  return(-1);
  }

 fprintf(fp_out,"%%ncurves=%d npts1=%d, npts2=%d\n", ncurves, npts[0], npts[1]);

/* First curve: */
 fprintf(fp_out,"%% Curve #1\n");
 fprintf(fp_out,"%d \n", npts[0]);
for(i = 0; i < npts[0]; i++) 
   fprintf(fp_out, "%f %f\n", xplot[i], yplot[i]);

/* Second curve: */
if(ncurves > 1) {
  fprintf(fp_out,"%% Curve #2\n");
  fprintf(fp_out,"%d \n", npts[1]);
  for(i = 0; i < npts[1]; i++) 
   fprintf(fp_out, "%f %f\n", xplot[npts_max + i], yplot[npts_max + i]);
  }

fclose(fp_out);
return(0);
}
/***************************************************************************
*  Extract notes from a measuremement line in Tokovinin's format
***************************************************************************/
static int extract_end_notes_from_dataline(char *in_line, char *notes)
{
int status = -1;
char *pc;

*notes = '\0';

// Go to first alphbetic character (' ', '.', '1', ... are not alphabetic)
pc = in_line;
while(*pc) {
  if(!isalpha(*pc)) pc++;
  else break;
  }

// Copy input line from the first alphabetic character:
if(*pc && isalpha(*pc)) {
   status = 0;
   strcpy(notes, pc);
   }

// Cleans notes from possible *\r" at the end:
pc = notes;
while(*pc && *pc != '\r' && *pc != '\n') pc++;
*pc = '\0';

return(status);
}
/***************************************************************************
* Read a measuremement line from data file in Tokovinin's format
*
***************************************************************************/
static int decode_data_line_toko1(char *in_line, double *wepoch, 
                                  double *wtheta, double *wrho, double *wdrho,
                                  double *wrv1, double *wdrv1, double *wrv2, 
                                  double *wdrv2, char *notes,
                                  int *is_rv1, int *is_rv2, int *is_theta_rho) 
{
int i, nval, is_SB2, status = -1;
double ww0, ww1, ww2, ww3, ww4, bessel_year;

/**************** Example: ***************************************
45533.4644  -10.69    0.51 Va COR
45543.4416  -11.25    0.51 Va COR
49547.4899    0.12    0.63 Vb COR
49592.4857    0.94    0.57 Vb COR
C RVM data
46632.364   -6.37 0.44 176 12.32  3.18 Va
46632.364   -0.20 0.80 176  6.78  3.18 Vb
47012.343   -4.26 0.31 565 22.42  2.48 Va+b
C Visual and interferometric data
      1963.0500    142.70   0.114000  0.05   I1 Fin1963a 26J  0
      1964.0341    171.40   0.140000  0.05   I1 Fin1964a 26J  1
      1965.0460    185.40   0.138000  0.05   I1 Fin1965a 26J  1
1971.6   275.6 0.21 0.04 I1 M1
1971.57  275.6 0.17 0.04 I1 M1
******************************************************************/

  *wepoch = 0.;
  *wtheta = 0.;
  *wrho = 0.;
  *wdrho = 0.; 
  *wrv1 = 0.; 
  *wdrv1 = 0.; 
  *wrv2 = 0.; 
  *wdrv2 = 0.; 
  *notes = '\0';
  *is_theta_rho = 0;
  *is_rv1 = 0;
  *is_rv2 = 0;
  strcpy(notes, "");

// Decode epoch:
  nval = sscanf(in_line, "%lf", &ww0);
  if(nval != 1) {
    fprintf(stderr, "decode_data_line_toko1/Error reading >%s<\n", in_line);
    return(1);
    }

// Extract notes
  extract_end_notes_from_dataline(in_line, notes);

// Check if interferometric measurement with date:
// epoch, theta, rho, drho measurements;
  if(ww0 < 2100.) {
    nval = sscanf(in_line,"%lf %lf %lf %lf", &ww0, &ww1, &ww2, &ww3);
    if(nval == 4) {
      *is_theta_rho = 1;
      *wepoch = ww0; 
      *wtheta = ww1; 
      *wrho = ww2; 
      *wdrho = ww3; 
      status = 0;
      }
  } else {
// Radial velocity:
// Conversion of the epoch from Julian days to Besselian years:
    julian_day_to_besselian_yr(ww0, &bessel_year);
    *wepoch = bessel_year; 

/*
* The  radial velocity  data  can be  in  two forms:
*
* (time, RV, error,  ...  Va, ...)  where Va  marks the primary component, 
* Vb marks the secondary. Lines with the V: tag are ignored.
*
* (time, RV1, err1, RV2, Err2, ... V2, ...) for SB2 (both velocities).
*/
    find_component_in_rv_notes(notes, is_rv1, is_rv2, &is_SB2);
    if(is_SB2) { 
       nval = sscanf(in_line,"%lf %lf %lf %lf %lf", &ww0, &ww1, &ww2, &ww3, &ww4);
       if(nval == 5) {
         *wrv1 = ww1; 
         *wdrv1 = ww2; 
         *wrv2 = ww3; 
         *wdrv2 = ww4; 
         status = 0;
         } else {
           fprintf(stderr, "Error/SB2: >%s< nval=%d\n", in_line, nval);
           status = 2;
         }
    } else { 
       nval = sscanf(in_line,"%lf %lf %lf", &ww0, &ww1, &ww2);
       if(nval == 3) {
        if(is_rv1) {
          *wrv1 = ww1; 
          *wdrv1 = ww2; 
          }
        if(is_rv2) {
          *wrv2 = ww1; 
          *wdrv2 = ww2; 
          }
        status = 0;
         } else {
           fprintf(stderr, "Error/Va or Vb: >%s< nval=%d\n", in_line, nval);
           status = 3;
         }
    }
  } // EOF ww0 >= 2100

return(status);
}
/**********************************************************************
* Look for Va V1 in notes:
* Va Vb Va+b V1 V2 can be found in radial velocity measurement lines
*
* WARNING: V2 means SB2 object with time RV1 err1 RV2 err2
**********************************************************************/
static int find_component_in_rv_notes(char *notes, int *is_rv1, int *is_rv2,
                                      int *is_SB2) 
{ 
int istart, found, i0;
char *pc, buffer[128];

 *is_rv1 = 0;
 *is_rv2 = 0;
 *is_SB2 = 0;

pc = notes;
istart = 0;
found = 0;
while(*pc) {
 if(*pc == 'V') break; 
 pc++;
 istart++;
 }

if(*pc == 'V') {
  pc++; 
  istart++;
// Case: 'Va' or 'V1'
  if((*pc == 'a') || (*pc == '1')) {
    *is_rv1 = 1; 
    i0 = istart + 1; found = 1;
    pc++; 
    istart++;
    if(*pc == '+') {
// Case: 'Va+b'
     pc++; 
     istart++;
     if(*pc == 'b') {
       *is_rv2 = 1;
       i0 = istart + 1; found = 1;
       }
    }
// Case: 'Vb'
  } else if(*pc == 'b') {
     *is_rv2 = 1; 
     i0 = istart + 1; found = 1;
// Case: 'V2'
  } else if(*pc == '2') {
    *is_rv1 = 1;
    *is_rv2 = 1;
    *is_SB2 = 1; 
    i0 = istart + 1; found = 1;
  }
}

// Truncate notes if component number was found:
 if(found) {
   strcpy(buffer, &notes[i0]);
   strcpy(notes, buffer);
   }

return(0);
}
/*****************************************************************************
* UpdatePanelText 
* prepare various wxString *_txt, ready for output to file
* 
* INPUT:
* panel_type: 0: orbit and data 
*            1: orbit only 
*            2: data only 
*            3: compute residuals 
******************************************************************************/
int GpFrame::UpdatePanelText(const int textpanel_type)
{
wxString Buffer_txt;

// Initialize to empty string:
  Buffer_txt = wxT("");

if((initialized != 1234) || (m_jlp_orbit1 == NULL)) return(-1);

// Save panel type to private variable:
m_textpanel_type = textpanel_type;

switch(m_textpanel_type) {
// panel_type=0: orbit and data
   case 0:
     PrepareOrbitText(Buffer_txt);
     PrepareDataText(Buffer_txt);
     ShowOrbitDataButton->SetValue(true); 
     break;
// panel_type=1: orbit only 
   case 1:
     PrepareOrbitText(Buffer_txt);
     ShowOrbitButton->SetValue(true); 
     break;
// panel_type=1: data only 
   case 2:
     PrepareDataText(Buffer_txt);
     ShowDataButton->SetValue(true); 
     break;
// panel_type=3: residuals 
   case 3:
     PrepareOrbitText(Buffer_txt);
     PrepareResidualsText(Buffer_txt);
     ShowResidualsButton->SetValue(true); 
     break;
  }

TextPanel_txt[m_textpanel_type] = Buffer_txt;
OrbitDataTextCtrl->SetValue(TextPanel_txt[m_textpanel_type]);

return(0);
}
/*****************************************************************************
* PrepareOrbitText
* prepare wxString *_txt, ready for output to file
* 
* INPUT/OUTPUT:
*   Buffer_txt
******************************************************************************/
int GpFrame::PrepareOrbitText(wxString &Buffer_txt)
{
wxString sstr;
double elmts[NELEMENTS], elmts_err[NELEMENTS];
double oelmnts0[NELEMENTS], oelmnts_err0[NELEMENTS];
int i0, orbit_type0;
int idx_free[5][N10];

if((initialized != 1234) || (m_jlp_orbit1 == NULL)) return(-1);

m_jlp_orbit1->GetOrbitalElements(elmts, elmts_err, N10, &orbit_type0);

ConvertElementToSpectro(elmts, elmts_err, oelmnts0, oelmnts_err0,
                        N10, orbit_type0);

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

/************************************************************/
// header:
/************************************************************/
sstr.Printf("Object: %s\n", object_name1);
Buffer_txt.Append(sstr);
if(right_ascension1 != 0. || declination1 != 0.) {
  sstr.Printf("R.A.:  %.3f\n", right_ascension1);
  Buffer_txt.Append(sstr);
  sstr.Printf("Dec:   %.3f\n", declination1);
  Buffer_txt.Append(sstr);
  }
if(orbit_equinox2 != 0.) {
  sstr.Printf("Equinox:  %.1f\n", orbit_equinox2);
  Buffer_txt.Append(sstr);
  }
sstr.Printf("P:     %.8g\n", oelmnts0[0]);
Buffer_txt.Append(sstr);
sstr.Printf("T:     %.8g\n", oelmnts0[1]);
Buffer_txt.Append(sstr);
sstr.Printf("e:     %.8g\n", oelmnts0[2]);
Buffer_txt.Append(sstr);

i0 = orbit_type0 - 1;
if(idx_free[i0][3] == 1)
  sstr.Printf("a:     %.8g\n", oelmnts0[3]);
else
  sstr.Printf("*a:         \n");
Buffer_txt.Append(sstr);

if(idx_free[i0][4] == 1)
  sstr.Printf("W:     %.8g\n", oelmnts0[4]);
else
  sstr.Printf("*W:         \n");
Buffer_txt.Append(sstr);

if(idx_free[i0][5] == 1)
  sstr.Printf("w:     %.8g\n", oelmnts0[5]);
else
  sstr.Printf("*w:         \n");
Buffer_txt.Append(sstr);

if(idx_free[i0][6] == 1)
  sstr.Printf("i:     %.8g\n", oelmnts0[6]);
else
  sstr.Printf("*i:         \n");
Buffer_txt.Append(sstr);

if(idx_free[i0][7] == 1)
  sstr.Printf("K1:    %.8g\n", oelmnts0[7]);
else
  sstr.Printf("*K1:         \n");
Buffer_txt.Append(sstr);

if(idx_free[i0][8] == 1)
  sstr.Printf("K2:    %.8g\n", oelmnts0[8]);
else
  sstr.Printf("*K2:         \n");
Buffer_txt.Append(sstr);

if(idx_free[i0][9] == 1)
  sstr.Printf("V0:    %.8g\n", oelmnts0[9]);
else
  sstr.Printf("*V0:         \n");
Buffer_txt.Append(sstr);

return(0);
}
/*****************************************************************************
* PrepareDataText
* prepare wxString *_txt, ready for output to file
* 
* INPUT/OUTPUT:
*   Buffer_txt
******************************************************************************/
int GpFrame::PrepareDataText(wxString &Buffer_txt)
{
wxString sstr;
double oelmnts0[NELEMENTS], oelmnts_err0[NELEMENTS];
double epoch_jd;
int i, orbit_type0;

if((initialized != 1234) || (m_jlp_orbit1 == NULL)) return(-1);

m_jlp_orbit1->GetOrbitalElements(oelmnts0, oelmnts_err0, N10, &orbit_type0);

//*************************************************
// Spectrocopic data: rv1
//*************************************************
  if(nrv1 > 0) {
  sstr.Printf("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n");
  Buffer_txt.Append(sstr);
  sstr.Printf("%% RVa:\n%% epoch(JD) rv sig_rv notes\n");
  Buffer_txt.Append(sstr);
  for(i = 0; i < nrv1; i++) {
// Conversion of the epoch to Julian days since radial velocities:
    besselian_yr_to_julian_day(epoch_rv1[i], &epoch_jd);
    sstr.Printf("%.4f %+.2f %.2f Va %s\n",
                epoch_jd, rv1[i], sigma_rv1[i], notes_rv1[i]);
    Buffer_txt.Append(sstr);
    }
  } // EOF nrv1 != 0
//*************************************************
// Spectrocopic data: rv2
//*************************************************
  if(nrv2 > 0) {
  sstr.Printf("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n");
  Buffer_txt.Append(sstr);
  sstr.Printf("%% RVb:\n%% epoch(JD) rv sig_rv notes\n");
  Buffer_txt.Append(sstr);
  for(i = 0; i < nrv2; i++) {
// Conversion of the epoch to Julian days since radial velocities:
    besselian_yr_to_julian_day(epoch_rv2[i], &epoch_jd);
    sstr.Printf("%.4f %+.2f %.2f Vb %s\n",
                epoch_jd, rv2[i], sigma_rv2[i], notes_rv2[i]);
    Buffer_txt.Append(sstr);
    }
  } // EOF nrv2 != 0
//*************************************************
// Visual data: rho, theta
//*************************************************
  if(nmeas1 > 0) {
  sstr.Printf("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n");
  Buffer_txt.Append(sstr);
  sstr.Printf("%% epoch theta rho sig_rho notes\n");
  Buffer_txt.Append(sstr);
  sstr.Printf("%%\n");
  Buffer_txt.Append(sstr);
  for(i = 0; i < nmeas1; i++) {
    sstr.Printf("%.4f %+.2f %.3f %.3f %s\n",
                epoch1[i], theta1[i], rho1[i], sigma_rho1[i],
                notes_meas1[i]);
    Buffer_txt.Append(sstr);
    }
  } // EOF nmeas1 != 0

return(0);
}
/***********************************************************************
* PrepareResidualsText for text panel
* Compute Residuals
* and load results to Buffer_txt (ready for output to file)
************************************************************************/
int GpFrame::PrepareResidualsText(wxString &Buffer_txt)
{
wxString sstr;
double epoch_jd;
double epoch_o, theta_c, rho_c, E_anom, wtheta, wrho; 
double rv1_c, rv2_c, wrv1, wrv2;
double sum_weight_rv1, sum_weight_rv2, wweight;
double sum_weight_rho, sum_weight_theta;
double sumsq_wrv1, sum_wrv1, sumsq_wrv2, sum_wrv2; 
double sumsq_wrho, sum_wrho, sumsq_wtheta, sum_wtheta, ww;
double K1, K2;
int i, orbit_type0, orbit_type_from_elements0;

if((initialized != 1234) || (m_jlp_orbit1 == NULL)) return(-1);

  mean_rho_resid1 = 0.;
  mean_theta_resid1 = 0.;
  mean_sigma_rho_resid1 = 0.;
  mean_sigma_theta_resid1 = 0.;
  mean_sigma_rv1_resid1 = 0.;
  mean_sigma_rv2_resid1 = 0.;

  orbit_type_from_elements0 = m_jlp_orbit1->GetOrbitTypeFromElements();

if(orbit_type_from_elements0 == 0) {
  sstr.Printf("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n"); 
  Buffer_txt.Append(sstr);
  sstr.Printf("No orbit loaded yet...\n"); 
  Buffer_txt.Append(sstr);
  return(-1);
  }

orbit_type0 = m_jlp_orbit1->GetOrbitType();
if(orbit_type0 != orbit_type_from_elements0) {
  sstr.Printf("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n"); 
  Buffer_txt.Append(sstr);
  sstr.Printf("Wrong orbit loaded: its type is incompatible with loaded data.\n"); 
  Buffer_txt.Append(sstr);
  return(-1);
  }

// Initialize Kepler parameters:
m_jlp_orbit1->Kepler_Init_Internal();

/************************************************************/
// Radial velocities
/************************************************************/
sum_wrv1 = 0.;
sumsq_wrv1 = 0.;
sum_wrv2 = 0.;
sumsq_wrv2 = 0.;
sum_weight_rv1 = 0.;
sum_weight_rv2 = 0.;
K1 = m_jlp_orbit1->GetK1();
K2 = m_jlp_orbit1->GetK2();

if(nrv1 > 0) {

sstr.Printf("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n"); 
Buffer_txt.Append(sstr);
sstr.Printf("%% RVa:\n%% epoch(JD) rv sig_rv rv_O-C notes\n"); 
Buffer_txt.Append(sstr);
for(i = 0; i < nrv1; i++) {

// Compute residuals
  epoch_o = epoch_rv1[i];
  Kepler_Ephemerid2(epoch_o, &rv1_c, &rv2_c, &E_anom);
  wrv1 = rv1[i] - rv1_c;

// Cumulative sums for computing the mean residuals:
  if(sigma_rv1[i] > 0) {
    wweight = K1 / sigma_rv1[i];
    sum_wrv1 += wrv1 * wweight;
    sumsq_wrv1 += SQUARE(wrv1) * wweight;
    sum_weight_rv1 += wweight;
    }

// Conversion of the epoch to Julian days since radial velocities:
  besselian_yr_to_julian_day(epoch_rv1[i], &epoch_jd);
  sstr.Printf("%.4f %+.2f %.2f %+.2f Va %s\n",
              epoch_jd, rv1[i], sigma_rv1[i], wrv1, 
              notes_rv1[i]);
  Buffer_txt.Append(sstr);
  }

// Mean residuals:
if(sum_weight_rv1 > 0.) {
  ww = sum_wrv1 / sum_weight_rv1;
  mean_sigma_rv1_resid1 = sqrt(sumsq_wrv1 / sum_weight_rv1 - ww * ww);
  sstr.Printf("%%\n%% RVa residuals: mean=%.2f sigma=%.2f\n%%\n", 
               ww, mean_sigma_rv1_resid1);
  Buffer_txt.Append(sstr);
  }

} // EOF nrv1 != 0

if(nrv2 > 0) {
sstr.Printf("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n"); 
Buffer_txt.Append(sstr);
sstr.Printf("%% RVb:\n%% epoch(JD) rv sig_rv weight rv_O-C notes\n");
Buffer_txt.Append(sstr);
sstr.Printf("%%\n"); 
Buffer_txt.Append(sstr);
for(i = 0; i < nrv2; i++) {

// Compute residuals
  epoch_o = epoch_rv2[i];
  Kepler_Ephemerid2(epoch_o, &rv1_c, &rv2_c, &E_anom);
  wrv2 = rv2[i] - rv2_c;

// Cumulative sums for computing the mean residuals:
  if(sigma_rv2[i] > 0) {
    wweight = K2 / sigma_rv2[i];
    sum_wrv2 += wrv2 * wweight;
    sumsq_wrv2 += SQUARE(wrv2) * wweight;
    sum_weight_rv2 += wweight;
    }

// Conversion of the epoch to Julian days 
  besselian_yr_to_julian_day(epoch_rv2[i], &epoch_jd);
  sstr.Printf("%.4f %+.2f %.2f %+.2f Vb %s\n", 
              epoch_jd, rv2[i], sigma_rv2[i], wrv2, 
              notes_rv2[i]);
  Buffer_txt.Append(sstr);
  }

if(sum_weight_rv2 > 0.) {
  ww = sum_wrv2 / sum_weight_rv2;
  mean_sigma_rv2_resid1 = sqrt(sumsq_wrv2 / sum_weight_rv2 - ww * ww);
  sstr.Printf("%%\n%% RVb residuals: mean=%.2f sigma=%.2f\n%%\n", 
               ww, mean_sigma_rv2_resid1);
  Buffer_txt.Append(sstr);
  }

} // EOF nrv2 != 0

/************************************************************/
// Theta rho measurements 
/************************************************************/

if(nmeas1 > 0) {

sum_wrho = 0.;
sumsq_wrho = 0.;
sum_wtheta = 0.;
sumsq_wtheta = 0.;
sum_weight_rho = 0.;
sum_weight_theta = 0.;

sstr.Printf("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n"); 
Buffer_txt.Append(sstr);
sstr.Printf("%% epoch theta rho sig_rho theta_O-C rho_O-C notes\n");
Buffer_txt.Append(sstr);
sstr.Printf("%%\n"); 
Buffer_txt.Append(sstr);
for(i = 0; i < nmeas1; i++) {
// Compute residuals
   epoch_o = epoch1[i];
   Kepler_Ephemerid1(epoch_o, &theta_c, &rho_c, &E_anom);

   wtheta = theta1[i] - theta_c;

// Will set residuals in the -90, +90 degrees range
   if(wtheta < -90.) wtheta += 360.;
// Adjust quadrant
   while(ABS(wtheta) > 90.) {
     if(wtheta > 0) wtheta -= 180.;
     else wtheta += 180.;
     }

// Cumulative sums for computing the mean residuals:
  if(sigma_rho1[i] > 0) {
     wweight = ABS(rho1[i]) / sigma_rho1[i];
     sum_wtheta += wtheta * wweight;
     sumsq_wtheta += SQUARE(wtheta) * wweight;
     sum_weight_theta += wweight;
     }

// Rho:
  wrho = rho1[i] - rho_c;

// Cumulative sums for computing the mean residuals:
  if(rho1[i] > 0 && sigma_rho1[i] > 0) {
     wweight = rho1[i] / sigma_rho1[i];
     sum_wrho += wrho * wweight;
     sumsq_wrho += SQUARE(wrho) * wweight;
     sum_weight_rho += wweight;
     }

  sstr.Printf("%.4f %+.2f %.3f %.3f %+.2f %+.3f %s\n", 
              epoch1[i], theta1[i], rho1[i], sigma_rho1[i], 
              wtheta, wrho, notes_meas1[i]);
  Buffer_txt.Append(sstr);
  }

// Mean standard deviation: 
if(sum_weight_rho > 0.) {
  mean_rho_resid1 = sum_wrho / sum_weight_rho;
  mean_sigma_rho_resid1 = sqrt(sumsq_wrho / sum_weight_rho 
                               - SQUARE(mean_rho_resid1));
  sstr.Printf("%%\n%% Rho residuals: mean=%.3f sigma=%.3f\n", 
               mean_rho_resid1, mean_sigma_rho_resid1);
  Buffer_txt.Append(sstr);
  }
if(sum_weight_theta > 0.) {
  mean_theta_resid1 = sum_wtheta / sum_weight_theta;
  mean_sigma_theta_resid1 = sqrt(sumsq_wtheta / sum_weight_theta 
                               - SQUARE(mean_theta_resid1));
  sstr.Printf("%% theta residuals: mean=%.2f sigma=%.2f\n%%\n", 
               mean_theta_resid1, mean_sigma_theta_resid1);
  Buffer_txt.Append(sstr);
  }

} // EOF nmeas1 != 0

return(0);
}
