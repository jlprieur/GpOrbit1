/*************************************************************************
* gpo_rw_datafiles.cpp
* To read/write visual data files
*
* Contains:
*
* JLP
* Version 29/06/2015
**************************************************************************/
#include <stdio.h> 
#include <stdlib.h>    // exit() 
#include <ctype.h>     // isdigit() 
#include <string.h>    // strcpy() 
#include <math.h>      // sin() 
#include "gpo_rw_files.h"

#define DEGTORAD 3.14159/180.

static int extract_end_notes_from_dataline(char *in_line, char *notes);
static int decode_mini_dataline(char *in_line, double *wepoch,
                                double *wtheta, double *wrho, double *wweight,
                                char *notes);
static int decode_WDS_dataline(char *in_line, double *wepoch,
                                double *wtheta, double *wrho, double *wweight,
                                char *notes);

/**********************************************************************
* Read visual data file in WDS format
*
* INPUT:
* in_data_file: name of the file
*
* OUTPUT:
* epoch0, theta0, rho0, sigma_rho0, weight0: measurement values                
* notes_meas: notes read in measurement line
* nmeas0: number of measurements that could be read in input file
*
**********************************************************************/
int read_WDS_visual_datafile(const char *in_data_file, double *epoch0,
                             double *theta0, double *rho0, double *weight0,
                             char notes_meas0[NMEAS_MAX][128],
                             const int nmeas_max, int *nmeas0)
{
double wepoch, wtheta, wrho, wweight;
char in_line[80], buffer[80], *pc;
FILE *fp_in;
int status, i, nval;
char notes[128];

if((fp_in = fopen(in_data_file, "r")) == NULL) {
  fprintf(stderr, "read_WDS_visual_datafile/Fatal error opening input file >%s<\n",
          in_data_file);
  return(-1);
  }

*nmeas0 = 0;
while(!feof(fp_in)) {
  if(fgets(buffer, 80, fp_in)) {

// Remove headings blanks if present:
    pc = buffer;
    while(*pc == ' ') pc++;

// Check if first character is a digit: */
    if(*pc && isdigit(*pc))  {

    strcpy(in_line, pc);
    status = decode_WDS_dataline(in_line, &wepoch, &wtheta, &wrho, &wweight,
                                  notes);
      if(status != 0){
        fprintf(stderr, "read_WDS_visual_datafile/Error reading line >%s<\n ichar=%d Last line or End of File ?)\n",
                in_line, (int)in_line[0]);
      } else {
      epoch0[*nmeas0] = wepoch;
      theta0[*nmeas0] = wtheta;
      rho0[*nmeas0] = wrho;
      weight0[*nmeas0] = wweight;
      strcpy(notes_meas0[*nmeas0], notes);
      if(*nmeas0 >= nmeas_max) {
        fprintf(stderr, "read_WDS_visual_datafile/Fatal error: nmeas0=%d reaches limit storage\n %s", *nmeas0);
        exit(1);
        }
      (*nmeas0)++;
      } // EOF status == 0
    } /* EOF in_line is_digit */
  } /* EOF fgets */
} /* EOF while !feof*/

fclose(fp_in);

return(0);
}
/**********************************************************************
* Read visual data file in minimum format
*
* INPUT:
* in_data_file: name of the file
*
* OUTPUT:
* epoch0, theta0, rho0, weight0: measurement values                
* notes_meas: notes read in measurement line
* nmeas0: number of measurements that could be read in input file
*
**********************************************************************/
int read_mini_visual_datafile(const char *in_data_file, double *epoch0,
                              double *theta0, double *rho0, double *weight0,
                              char notes_meas0[NMEAS_MAX][128],
                              const int nmeas_max, int *nmeas0)
{
double wepoch, wtheta, wrho, wweight;
char in_line[80], buffer[80], *pc;
FILE *fp_in;
int status, i, nval;
char notes[128];

if((fp_in = fopen(in_data_file, "r")) == NULL) {
  fprintf(stderr, "read_mini_visual_datafile/Fatal error opening input file >%s<\n",
          in_data_file);
  return(-1);
  }

*nmeas0 = 0;
while(!feof(fp_in)) {
  if(fgets(buffer, 80, fp_in)) {

// Remove headings blanks if present:
    pc = buffer;
    while(*pc == ' ') pc++;

// Check if first character is a digit: */
    if(*pc && isdigit(*pc))  {

    strcpy(in_line, pc);
    status = decode_mini_dataline(in_line, &wepoch, &wtheta, &wrho, &wweight,
                                  notes);
      if(status != 0){
        fprintf(stderr, "read_mini_visual_datafile/Error reading line >%s<\n ichar=%d Last line or End of File ?)\n",
                in_line, (int)in_line[0]);
      } else {
      epoch0[*nmeas0] = wepoch;
      theta0[*nmeas0] = wtheta;
      rho0[*nmeas0] = wrho;
      weight0[*nmeas0] = wweight;
      strcpy(notes_meas0[*nmeas0], notes);
      if(*nmeas0 >= nmeas_max) {
        fprintf(stderr, "read_mini_visual_datafile/Fatal error: nmeas0=%d reaches limit storage\n %s", *nmeas0);
        exit(1);
        }
      (*nmeas0)++;
      } // EOF status == 0
    } /* EOF in_line is_digit */
  } /* EOF fgets */
} /* EOF while !feof*/

fclose(fp_in);
 
return(0);
}
/***************************************************************************
* Read a visual measuremement line from data file in minimum format
*
* Format:
* epoch rho theta weight notes
***************************************************************************/
static int decode_mini_dataline(char *in_line, double *wepoch,
                                double *wtheta, double *wrho, double *wweight,
                                char *notes)
{
int i, nval, status = -1;
double ww0, ww1, ww2, ww3;

  *wepoch = 0.;
  *wtheta = 0.;
  *wrho = 0.;
  *wweight = 0.;
  strcpy(notes, "");

// Decode epoch, theta, rho, drho measurements;
  nval = sscanf(in_line,"%lf %lf %lf %lf", &ww0, &ww1, &ww2, &ww3);
  if(nval != 4) {
    fprintf(stderr, "decode_mini_dataline/Error reading >%s<\n", in_line);
    return(1);
    } else {
      *wepoch = ww0;
      *wtheta = ww1;
      *wrho = ww2;
      *wweight = ww3;
// Extract notes
      extract_end_notes_from_dataline(in_line, notes);
      status = 0;
      }

return(0);
}
/***************************************************************************
* Read a visual measuremement line from data file in WDS format
* or extended WDS format
*
* Format:
* epoch rho theta nights observer diameter [weight]
*
* Example of extended WDS format:
%% WDS extended format
%% epoch rho_arcsec theta_deg nights observer telescope_diameter_cm Marco's_weight
 1781.790   -5.120  338.018  1 H    20   .4
 1802.080   -1.000  333.560  1 H    20   .4
 1821.930   -5.428  336.004  2 SHJ  10   .6
***************************************************************************/
static int decode_WDS_dataline(char *in_line, double *wepoch,
                                double *wtheta, double *wrho, double *wweight,
                                char *notes)
{
int i, nval, status = -1, diam, diameter;
double ww0;
char *pc;

  *wepoch = 0.;
  *wtheta = 0.;
  *wrho = 0.;
  *wweight = 0.;
  *notes = '\0';

//***********************************************************************
// 1. Read epoch:
// Remove headings blanks if present:
  pc = in_line;
  while(*pc == ' ') pc++;

  nval = 0;
  if(isdigit(*pc)) nval = sscanf(pc, "%lf", &ww0);
  if(nval != 1) {
    fprintf(stderr, "decode_WDS_dataline/Error reading epoch >%s<\n", in_line);
    return(1);
    }
  *wepoch = ww0;

//***********************************************************************
// 2. Read rho:
// Skip digits:
  while(isdigit(*pc) || *pc == '.' || *pc == '-') pc++;
// Skip blanks:
  while(*pc == ' ') pc++;
  nval = 0;
  if(isdigit(*pc) || *pc == '-') nval = sscanf(pc, "%lf", &ww0);
  if(nval != 1) {
    fprintf(stderr, "decode_WDS_dataline/Error reading rho >%s<\n", in_line);
    return(1);
    }
  *wrho = ww0;

//***********************************************************************
// 3. Read theta:
// Skip digits:
  while(isdigit(*pc) || *pc == '.' || *pc == '-') pc++;
// Skip blanks:
  while(*pc == ' ') pc++;
  nval = 0;
  if(isdigit(*pc) || *pc == '-') nval = sscanf(pc, "%lf", &ww0);
  if(nval != 1) {
    fprintf(stderr, "decode_WDS_dataline/Error reading theta >%s<\n", in_line);
    return(1);
    }
  *wtheta = ww0;

//***********************************************************************
// 4. Read notes containing observer name:

// Skip blanks and nber of nights:
  while(isdigit(*pc) || *pc == '.' || *pc == '-' || *pc == ' ') pc++;
  while(*pc && !isalpha(*pc)) pc++;

// Copy notes:
  i = 0;
  while(*pc && *pc != '\r' && *pc != '\n' && *pc != ' ') {
    notes[i++] = *pc; 
    pc++;
    }
  notes[i] = '\0';

// End of line, so I return from here:
  if(*notes == '\0' || *pc != ' ') return(0);

//***********************************************************************
// 5. Read diameter (optional parameter):
  while(*pc == ' ') pc++;
  nval = 0;
  if(isdigit(*pc)) nval = sscanf(pc, "%d", &diam);
  if(nval != 1) return(0); 
  diameter = diam;

//***********************************************************************
// 6. Read weight (optional parameter):
// Skip digits:
  while(isdigit(*pc) || *pc == '.' || *pc == '-') pc++;
// Skip blanks:
  while(*pc == ' ') pc++;
  nval = 0;
  if(isdigit(*pc)) nval = sscanf(pc, "%lf", &ww0);
  if(nval != 1) return(0); 
  *wweight = ww0;

return(0);
}
/***************************************************************************
*  Extract notes from a measuremement line in minimum format
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
*  Extract notes from a measuremement line in Tokovinin's format
***************************************************************************/
static int toko_extract_notes_from_dataline(char *in_line, char *notes)
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
/**********************************************************************
* Computing the precession correction to be applied to positions angles
* using the formula of Armellini (1931) (or Couteau, 1978)
*
* INPUT :
* theta: before correction (in radians)
* alpha, delta: coordinates of object (in radians)
* epoch_o: epoch of observation
* orbit_equinox: equinox used as a reference for computing the orbit
*
* OUTPUT :
* dtheta_precess: correction for precession (in radians)
***********************************************************************/
int precession_correction(double *dtheta_precess, double alpha, 
                          double delta, double epoch_o, double orbit_equinox)
{

/* Precession correction of theta in arcseconds
*/
  *dtheta_precess = -20.0 * (epoch_o - orbit_equinox) * sin(alpha) / cos(delta);

#ifdef DEBUG_1
  printf(" epoch=%f equinox=%f sin(alpha)=%f cos(delta)=%f\n",
         epoch_o, orbit_equinox, sin(alpha), cos(delta));
  printf(" correction for precession: %f (arcsec) or %f (degrees)\n",
         *dtheta_precess,  *dtheta_precess/3600.);
#endif

/* Conversion to radians: */
  *dtheta_precess *= DEGTORAD/3600.;

return(0);
}
