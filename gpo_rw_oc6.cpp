/*************************************************************************
* gpo_rw_oc6.cpp
* To read/write OC6 visual orbit files
* from "jlp/src/orbits/OC6_catalog_utils.c"
*
* Contains:
*
* JLP
* Version 28/06/2015
**************************************************************************/
#include <stdio.h> 
#include <stdlib.h>    // exit() 
#include <string.h>    // strcpy() 
#include <ctype.h>                   /* isprint... */
#include <math.h>
#include <time.h>                    /* date */

#include "gpo_defs.h"  // MINI, MAXI, etc
#include "gpo_frame.h"
#include "jlp_string.h"    // compact_string()

#include "gpo_defs.h"      // YEAR_TO_DAYS

/* The prototypes of routines included here
* are defined in "gpo_rw_oc6.h":
*/
#include "gpo_rw_oc6.h"

/*
#define DEBUG
*/

#define DEGTORAD   (PI/180.00)

static int get_coordinates_from_WDS_name(char *WDS_name, 
                                         double *right_ascension0, 
                                         double *declination0);

/*****************************************************************************
* {'P', 'T', 'e', 'a', 'W', 'w', 'i', 'K1', 'K2', 'V0'};
* orbit_type = 1 visual & spectro SB1 (9 elements idx1= 1 1 1 1 1 1 1 1 0 1)
* orbit_type = 2 visual & spectro SB2 (10 elements idx1= 1 1 1 1 1 1 1 1 1 1)
* orbit_type = 3 visual (7 elements: idx2= 1 1 1 1 1 1 1 0 0 0)
* orbit_type = 4 spectroscopic SB1 (6 elements: idx3= 1 1 1 0 0 1 0 1 0 1)
* orbit_type = 5 spectroscopic SB2 (7 elements: idx3= 1 1 1 0 0 1 0 1 1 1)
*****************************************************************************/

/**********************************************************************
* Read visual orbit file in OC6 format
*
**********************************************************************/
int rd_visual_oelements_OC6(const char *orbit_infile, 
                            double *oelemnt0, double *orbit_equinox0,
                            char *object_name0, double *right_ascension0,
                            double *declination0, int *orbit_type0)
{
float Omega_node, omega_peri, i_incl, e_eccent, T_periastron, orbit_equinox;
float mean_motion, a_smaxis, Period;
int iline, status, is_master_file, line_length;
char object_name[60], discov_name[40], comp_name[10], WDS_name[40];
char ADS_name[40], author[60], *pc;
/* Maximum line seems to be 265 for OC6 catalog... */
char in_line1[300], in_line2[300], buffer[128];
FILE *fp_in;

strcpy(object_name0, "");
*orbit_type0 = 0;
*right_ascension0 = 0.;
*declination0 = 0.;

/* Open input OC6-formatted file containing the orbital parameters: */
if((fp_in = fopen(orbit_infile, "r")) == NULL) {
   fprintf(stderr, "rd_visual_oelements_OC6/Fatal error opening input file: %s\n",
           orbit_infile);
    return(-1);
  }

iline = 0;
while(!feof(fp_in)) {
  if(fgets(in_line1,300,fp_in)) {
    line_length = (int)strlen(in_line1);
/* strlen(line1) = 279 if master file
*  strlen(line1) = 265 if OC6 file
*/
    is_master_file = (line_length > 270) ? 1 : 0;
    iline++;
/* Commented lines can start with % or # : */
    if(in_line1[0] != '%' && in_line1[0] != '#') {

/* Check if it is a very short line
*/
    if(line_length < 10) {
      printf("WARNING line #%d is very short (length=%d): >%s<\n",
              iline, line_length, in_line1);
    } else {
    status = get_orbit_from_OC6_list(in_line1, iline, is_master_file, WDS_name, 
                            ADS_name, discov_name, comp_name, object_name, 
                            author, &Omega_node, &omega_peri, &i_incl, 
                            &e_eccent, &T_periastron, &Period, 
                            &a_smaxis, &mean_motion, &orbit_equinox);
    if(status) {
     fprintf(stderr, "rd_visual_oelements_0C6/Error reading orbital parameters in line #%d (status=%d)\n", iline, status);
      } 
// Break at first uncommented long line:
     break;
    } // EOF line_length >= 10
  } // EOF not commented line
} // EOF fgets
} // EOF !feof
 
fclose(fp_in);

// ckeyword[NELEMENTS] = {'P', 'T', 'e', 'a', 'W', 'w', 'i', 'K1', 'K2', 'V0'};

if(status == 0) {
// The name of the object is retrieved from the orbit filename
strcpy(buffer, orbit_infile);

// Remove extension:
pc = buffer;
while(*pc && *pc != '.') pc++;
*pc = '\0';

// Remove directory:
pc--;
while(*pc != buffer[0] && *pc != '/' && *pc != '\\') pc--;
if(*pc != '/' || *pc != '\\') pc++;

strcpy(object_name0, pc);

// Get approximative coordinates from WDS name:
get_coordinates_from_WDS_name(WDS_name, right_ascension0, declination0);

#ifdef DEBUG
printf("WDS_name = %s  ADS_name = %s\n", WDS_name, ADS_name);
printf("discov_name = %s  comp_name = %s\n", discov_name, comp_name);
printf("author = %s\n", author);
#endif

oelemnt0[0] = Period; 
oelemnt0[1] = T_periastron; 
oelemnt0[2] = e_eccent; 
oelemnt0[3] = a_smaxis; 
oelemnt0[4] = Omega_node; 
oelemnt0[5] = omega_peri; 
oelemnt0[6] = i_incl; 
*orbit_type0 = 3;
*orbit_equinox0 = orbit_equinox;

#ifdef DEBUG
for(int i = 0; i < N7; i++) printf("element[%d]=%f\n", i, oelemnt0[i]);
#endif
  }

return(status);
}
/***************************************************************************
* get_orbit_from_OC6_list
* Read input line from a list of selected orbits in OC6 format 
* and retrieve the orbital parameters and the measurements
* in OC6 ("Sixth Orbital Catalog") format (i.e., iformat = 2 or -2)
*
* Input lines contain
*  - name of object and orbital parameters in OC6 format (if iformat = 2 or -2) 
*  - the measurements (if iformat = -2)
*    otherwise, (if format = 2) the program will retrieve the measurements 
*      from the calibrated table.
*
* Example:
*
* With iformat = 2:
*
* 000000.00+000000.0 00093+7943 STF   2          102    431    760   6.68   6.89    540.      y    .         0.995  a   .      110.1       .     171.2        .      1887.5     y    .       0.715     .       333.7       .  
   2000      3 n Hei1997  abcdefghikj.png
* 092059.40+381117.9 09210+3811 STF1338AB       7307  80441  45858   6.72   7.08    444.27    y    .         1.624  a   .       33.4       .     177.4        .      1983.69    y    .       0.247     .        83.6       .          1999 2002-03   3.4   n Sca2002b wds09210+3811b.png
*
*
* INPUT:
* in_line: line from the list extracted from the OC6 catalog,
*           corresponding to the object
* iline: number of the corresponding line in input file
* is_master_file: flag set to 1 if master file ("orb6.master", 
*                          to 0 if OC6 file ("orb6orbits.txt")
* 
* OUTPUT:
* object designation, orbital parameters and arrays with the measurements
*
***************************************************************************/
int get_orbit_from_OC6_list(char *in_line, int iline, int is_master_file, 
                       char *WDS_name, char *ADS_name, char *discov_name, 
                       char *comp_name, char *object_name, char *author, 
                       float *Omega_node, float *omega_peri, float *i_incl, 
                       float *e_eccent, float *T_periastron, float *Period, 
                       float *a_smaxis, float *mean_motion, 
                       float *orbit_equinox)
{
char buffer[80];
int nval;
double ww;

#ifdef DEBUG
printf("get_orbit_from_OC6_list/iline=%d\n %s", iline, in_line);
#endif

/* Retrieve the object designation and the author of the orbit
*/
strncpy(WDS_name, &in_line[19], 10);
WDS_name[10] = '\0';

/* Discover's name in fixed format with 7 characters (like with the WDS) */
strncpy(discov_name, &in_line[30], 7);
discov_name[7] = '\0';

strncpy(comp_name, &in_line[37], 7);
comp_name[7] = '\0';
/* Remove all blanks: */
compact_string(comp_name, 10);

strncpy(ADS_name, &in_line[45], 5);
ADS_name[5] = '\0';
/* Remove all blanks: */
compact_string(ADS_name, 40);

/* Period: */
strncpy(buffer, &in_line[80], 11);
buffer[11] = '\0';
if(sscanf(buffer, "%f", Period) != 1) {
  fprintf(stderr, "WDS=%s Error reading period: buffer=%s\n", WDS_name, buffer);
  return(-1);
  }

// Conversion of period to years:
/* Input file units: minutes, hours, days, centuries, or years */
if(in_line[92] == 'm') 
  *Period /= (YEAR_TO_DAYS * 24. * 60.);
else if(in_line[92] == 'h') 
  *Period /= (YEAR_TO_DAYS * 24.);
else if(in_line[92] == 'd') 
  *Period /= YEAR_TO_DAYS;
else if(in_line[92] == 'c') 
  *Period *= 100.;
else if(in_line[92] != 'y') 
  {
  fprintf(stderr, "Error reading period of %s: unknown unit (P=%.20s unit=%c)\n", 
          WDS_name, &in_line[81], in_line[92]);
  return(-1);
  }

/* Semi-major axis: */
strncpy(buffer, &in_line[105], 9);
buffer[9] = '\0';
*a_smaxis = 0.;
if((nval = sscanf(buffer, "%f", a_smaxis)) != 1) {
  fprintf(stderr, "Error reading semi-major axis: nval=%d buffer=%s\n", 
          nval, buffer);
  fprintf(stderr, "iline=%d a_smaxis=%f inline=%s\n", iline, *a_smaxis, in_line);
  return(-1);
  }
/* Units: milliarcseconds or arcseconds */
if(in_line[114] == 'm') 
  *a_smaxis *= 1000.;
else if(in_line[114] != 'a') 
  {
  fprintf(stderr, "Error reading semi-major axis: unknown unit (a=%.20s unit=%c)\n", 
          &in_line[105], in_line[114]);
  return(-1);
  }

/* Inclination (degrees) */
strncpy(buffer, &in_line[125], 8);
buffer[8] = '\0';
if(sscanf(buffer, "%f", i_incl) != 1) {
  fprintf(stderr, "Error reading inclination: buffer=%s\n", buffer);
  return(-1);
  }

/* Node, Omega (degrees) */
strncpy(buffer, &in_line[143], 8);
buffer[8] = '\0';
if(sscanf(buffer, "%f", Omega_node) != 1) {
  fprintf(stderr, "Error reading Omega (node): buffer=%s\n", buffer);
  return(-1);
  }

/* T of periastron passage */
strncpy(buffer, &in_line[162], 12);
buffer[12] = '\0';
if(sscanf(buffer, "%f", T_periastron) != 1) {
  fprintf(stderr, "WDS_name=%s /error reading T periastron: buffer=%s\n", 
          WDS_name, buffer);
  return(-1);
  }
/* Units: modified Julian date or fractionnal Besselian year */
/* The time of periastron passage (T0) and code for units:
                            d = Julian date (-2,400,000 days)
                            m = modified Julian date (MJD = JD-2,400,000.5 days)
                            y = fractional Besselian year
*/

#ifdef DEBUG
printf("ADS_name=%s raw periastron = %f yr (unit=%c)\n", 
        ADS_name, *T_periastron, in_line[174]);
#endif

/* Modified JD: JD - 2400000
*  reduced JD:  JD - 2400000.5
*/
if(in_line[174] == 'm')
  {
/* JLP 2015... */
  *T_periastron = 1900.0 + (*T_periastron - 15020.31352 - 0.5) / 365.242198781; 
  }
else if(in_line[174] == 'd')
  {
/* d = modified Julian date (JulianDate - 2,400,000 days) */
/* Bessel epoch = 1900.0 + (JulianDate - 2415020.31352) / 365.242198781 */
  *T_periastron = 1900.0 + (*T_periastron - 15020.31352) / 365.242198781; 
  }
else if(in_line[174] != 'y') 
  {
  fprintf(stderr, "Error reading periastron: unknown unit (P=%.20s unit=%c)\n", 
          &in_line[162], in_line[174]);
  return(-1);
  }

/* Eccentricity */
strncpy(buffer, &in_line[187], 8);
buffer[8] = '\0';
if(sscanf(buffer, "%f", e_eccent) != 1) {
  fprintf(stderr, "Error reading eccentricity: buffer=%s\n", buffer);
  return(-1);
  }

/* omega (longitude of periastron) */
strncpy(buffer, &in_line[205], 8);
buffer[8] = '\0';
if(sscanf(buffer, "%f", omega_peri) != 1) {
  fprintf(stderr, "Error reading omega_periastron: buffer=%s\n", buffer);
  return(-1);
  }

/* Equinox */
strncpy(buffer, &in_line[223], 4);
buffer[4] = '\0';
/* Default value is 2000.0 */
if(sscanf(buffer, "%f", orbit_equinox) != 1) {
  *orbit_equinox = 2000.0;
}

/* Reference of orbit (author) */
if(is_master_file) 
   strncpy(author, &in_line[251], 8);
  else
   strncpy(author, &in_line[237], 8);
author[8] = '\0';
/* Remove heading and trailing blanks: */
trim_string(author, 60);

/* Restriction of the object name to ADS name or discoverer name */
if(ADS_name[0] != '\0' && ADS_name[0] != '.') 
      sprintf(object_name, "ADS %s", ADS_name);
else if (discov_name[0] != '\0') strcpy(object_name, discov_name);

#ifdef DEBUG_1
printf("Object=%s WDS=%s ADS=%s discov=%s comp=%s author=%s\n",
        object_name, WDS_name, ADS_name, discov_name, comp_name, author);

 printf("Omega_node=%.3f omega_peri=%.3f incl=%.3f e=%.4f T=%.3f P=%.3f a=%.5f Equinox=%.3f\n", 
        *Omega_node, *omega_peri, *i_incl, *e_eccent, *T_periastron, 
        *Period, *a_smaxis, *orbit_equinox);
#endif

 *mean_motion = (360.0 / *Period);

return(0);
}
/***************************************************************************
* line_extraction_from_OC6_catalog ("orb6.master" in November 2009)
* Return the line in OC6 that corresponds to the input object name 
*
* INPUT:
* OC6_fname: name of master file or Sixth Orbit catalog 
* is_master_file: flag set to 1 if master file ("orb6.master", 
*                          to 0 if OC6 file ("orb6orbits.txt")
* ads_name[60]: ADS name of object ('\0' if not in ADS)
* discov_name[60]: discovery name of object 
* comp_name: name of the companion (e.g., AB, Aa, etc) 
* norbits_per_object: maximum number of orbits to be loaded to output list
*                     fror each object
* 
* OUTPUT:
* fp_out: pointer to the output ASCII file 
* found: flag set to one if at least one orbit was found for object
*
***************************************************************************/
int line_extraction_from_OC6_catalog(char *OC6_fname, int is_master_file,
                                     char *ads_name, 
                                     char *discov_name, char *comp_name, 
                                     FILE *fp_out, int *found, 
                                     int *candidate_found, 
                                     int norbits_per_object)
{
/* Line length is 278 + "\n" for OC6 catalog (orb6.master, november 2009)... */
/* Line length is 264 + "\n" for OC6 catalog (orb6orbits.txt, november 2009)... */
int iline, status, max_norbits = 1024, norbits, imin, nlines_in_header;
int OC6_comp_is_AB, comp_is_AB, discov_name_only;
char OC6_comp_really_compacted[40], comp_really_compacted[40];
char OC6_ads_name[60], OC6_comp_name[40], OC6_discov_name[40];
char line_buffer[300], compacted_ads_name[60], compacted_comp_name[40];
char compacted_discov_name[40];
/* Assume that the number of orbits for one object is always less than 1024: */
char object_orbits[300*1024];
FILE *fp_in;
register int i;
// To print information on screen:
int italk = 0;

*found = 0;
*candidate_found = 0;
norbits = 0;

strcpy(compacted_ads_name, ads_name);
compact_string(compacted_ads_name, 60);

if(compacted_ads_name[0] == '\0') discov_name_only = 1;
else discov_name_only = 0; 

strcpy(compacted_comp_name, comp_name);
compact_string(compacted_comp_name, 40);
really_compact_companion(compacted_comp_name, comp_really_compacted, 40);

strcpy(compacted_discov_name, discov_name);
compact_string(compacted_discov_name, 40);

#ifdef DEBUG
 printf("CURRENT OBJECT: ads_name=%s comp_name=%s discov_name=%s \n", 
ads_name, comp_name, discov_name);
#endif

/* Open OC6 catalog: */
if((fp_in = fopen(OC6_fname,"r")) == NULL) {
  fprintf(stderr, "line_extraction_from_OC6_catalog/Fatal error opening %s\n",
          OC6_fname);
  exit(-1);
  }

iline = 0;

if(is_master_file)
  nlines_in_header = 4;
else
  nlines_in_header = 8;

while(!feof(fp_in)) {
/* Should read more than 
* 278 characters in master OC6 file
* 264 characters in non-master OC6 file
* in order to be sure to copy the complete line: */
  if(fgets(line_buffer, 280, fp_in)) {
    iline++;
    if(*found && !strncmp(line_buffer,"    ",4)) break; 
/* Skip the header and empty lines (that
* generally indicate the end of a given object):*/
    if(iline > nlines_in_header && strncmp(line_buffer,"    ",4)) {
/* Get the object name of each line */
     status = get_name_from_OC6_line(line_buffer, OC6_ads_name,
                                     OC6_discov_name, OC6_comp_name);
     if(status) {
     fprintf(stderr, "line_extraction_from_OC6_catalog/Error processing line #%d\n", iline); 
     return(-1);
     }
#ifdef DEBUG_1
if(iline < nlines_in_header+3) 
     printf("OK1found=%d iline=%d: >%s<\n ads=%s< disc=%s< comp=%s<\n", 
             *found, iline, line_buffer, OC6_ads_name, OC6_discov_name, 
             OC6_comp_name);
#endif
/* If "ads_name" was found, copy the full line containing 
* the orbital elements: */
     compact_string(OC6_ads_name, 60);
/* CASE 1 */
    if(discov_name_only){
    if(!strcmp(OC6_discov_name, compacted_discov_name)) {
        for(i = 0; i < 300; i++) 
           object_orbits[i + 300*norbits] = line_buffer[i];
        *found = 1;
        norbits++;
        if(norbits >= max_norbits) {
            fprintf(stderr, "line_extraction_from_OC6_catalog/Fatal error: %d orbits found for a single object!\n", norbits);
            fprintf(stderr, "(OC6_discov_name=%s)\n", OC6_discov_name);
            exit(-1);
            }
      }
/* CASE 2 */
    } else if(!strcmp(compacted_ads_name, OC6_ads_name)
         && !strcmp(OC6_discov_name, compacted_discov_name)) {
/* Test on the companion names if present in the object name: 
*/
     compact_string(OC6_comp_name, 40);
     really_compact_companion(OC6_comp_name, OC6_comp_really_compacted, 40);
#ifdef DEBUG
printf("DEBUG/OC6_object=%s OC6_companion=%s| comp=%s| (really compacted: OC6 comp=%s comp=%s)\n", 
        OC6_ads_name, OC6_comp_name, comp_name, OC6_comp_really_compacted, 
        comp_really_compacted);
#endif
        comp_is_AB = 0;
        if((compacted_comp_name[0] == '\0') 
           || !strcmp(compacted_comp_name,"AB")
           || !strncmp(compacted_comp_name,"Aa-B",4)) comp_is_AB = 1;
        OC6_comp_is_AB = 0;
        if((OC6_comp_name[0] == '\0') || !strcmp(OC6_comp_name,"AB")
           || !strncmp(OC6_comp_name,"Aa-B",4)) OC6_comp_is_AB = 1;
           
        if((*compacted_comp_name != '\0' 
             && !strcmp(OC6_comp_name, compacted_comp_name))
/* If not mentionned in object name, should
* be either not mentioned in OC6 or equal to AB: */
           || (comp_is_AB && OC6_comp_is_AB)
           || (*comp_really_compacted != '\0' && 
               !strcmp(comp_really_compacted, OC6_comp_really_compacted))){
        for(i = 0; i < 300; i++) 
           object_orbits[i + 300*norbits] = line_buffer[i];
        *found = 1;
        if(*found && italk > 0) {
          printf("From_OC6_cat/Object found now in OC6: >%s< >%s< >%s<\n", 
                  OC6_ads_name, OC6_discov_name, OC6_comp_name); 
          }
        norbits++;
        if(norbits >= max_norbits) {
            fprintf(stderr, "line_extraction_from_OC6_catalog/Fatal error: %d orbits found for a single object!\n", norbits);
            fprintf(stderr, "(OC6_ads_name=%s OC6_discov_name=%s)\n", 
                    OC6_ads_name, OC6_discov_name);
            exit(-1);
            }
        } else if(!(*found)) {
        if(italk > 0) {
        printf("CURRENT OBJECT: ads_name=%s comp_name=%s discov_name=%s \n",
	 ads_name, comp_name, discov_name);
        printf("From_OC6_cat/Not yet found, possible candidate in OC6: >%s< >%s< >%s< (companion names look different though...)\n", 
               OC6_ads_name, OC6_discov_name, OC6_comp_name); 
        }
         *candidate_found = 1;
        }
#ifdef DEBUG
  printf("DEBUG/iline=%d OC6_ads_name=>%s< OC6_comp_name=>%s<\n", 
          iline, OC6_ads_name, OC6_comp_name);
#endif
     } /* EOF if(!strcmp ads_name ...) */
    } /* EOF if line > nlines_header */
  } /* EOF if fgets */ 
 }

#ifdef DEBUG
printf("line_extraction_from_OC6_catalog: %d lines read and %d orbits found for current object\n", 
        iline, norbits);
#endif

/* number of orbits per object: 
* 0=all 1=last 2=last two orbits, etc.
*/
 
if(*found) {
if(norbits_per_object > 0) 
  imin = MAXI(norbits - norbits_per_object, 0);
else
  imin = 0;

for(i = norbits - 1; i >= imin; i--)
   fprintf(fp_out, "%s", &object_orbits[i*300]);
}

fclose(fp_in);
return(0);
}
/***************************************************************************
* get_name_from_OC6_line
* Extract the object and companion names from the line in OC6 
* that corresponds to the input object name 
*
*
* Example:
*
* 000000.00+000000.0 00093+7943 STF   2          102    431    760   6.68   6.89    540.      y    .         0.995  a   .      110.1       .     171.2        .      1887.5     y    .       0.715     .       333.7       .  
   2000      3 n Hei1997  abcdefghikj.png
*
*
* INPUT:
* in_line: full line of OC6 corresponding to the object
* 
* OUTPUT:
* OC6_ads_name: ADS name of object (NULL if not in ADS catalog) 
* OC6_comp_name: name of the companion (e.g., AB, Aa, etc) 
*
***************************************************************************/
int get_name_from_OC6_line(char *in_line, char *OC6_ads_name,
                           char *OC6_discov_name, char *OC6_comp_name) 
{
char WDS_name[40], ads_name[40], discov_name[40], comp_name[40];

/* Retrieve the object designation and the author of the orbit
*/
strncpy(WDS_name, &in_line[19], 10);
WDS_name[10] = '\0';

/* Discover's name in fixed format with 7 characters */
strncpy(discov_name, &in_line[30], 7);
discov_name[7] = '\0';
/* Remove all blanks: */
compact_string(discov_name, 10);
strcpy(OC6_discov_name, discov_name);

/* Companion name */
strncpy(comp_name, &in_line[37], 7);
comp_name[7] = '\0';
/* Remove all blanks: */
compact_string(comp_name, 10);
strcpy(OC6_comp_name, comp_name);

/* ADS name */
OC6_ads_name[0] = '\0';
strncpy(ads_name, &in_line[45], 5);
ads_name[5] = '\0';
/* Remove all blanks: */
compact_string(ads_name, 40);
if(ads_name[0] == '.') ads_name[0] = '\0'; 
if(ads_name[0] != '\0') sprintf(OC6_ads_name, "ADS %s", ads_name);

#ifdef DEBUG_1
  printf("WDS=%s OC6_ADS=%s OC6_discov=%s OC6_comp=%s\n",
          WDS_name, OC6_ads_name, OC6_discov_name, OC6_comp_name);
#endif

return(0);
}
/******************************************************************
* Convert companion name to avoid problems when finding orbits
*
* e.g.  Bb-C  -> BC
*       Aa-Bb -> AB
*       AaBb -> AB
*       Cc-D -> CD
*       Cc,D -> CD
*       AB-C -> AC
*       AC   -> AC
*****************************************************************/
int really_compact_companion(char *name_in, char *name_out, int length)
{
char *pc;
int k, separation_found, main_found;

name_in[length-1] = '\0';

/******* Look for separation */
pc = name_in;
separation_found = 0;
k = 0;
while(*pc) {
  if(*pc == '-' || *pc == ',') {
  separation_found = 1;
  } else if (isupper(*pc)) {
   name_out[k++] = *pc;
  }
  pc++;
  }
name_out[k] = '\0';
/* Job done if no separation: */
if(!separation_found) return(0);

/******* Look for main components if separation: */
pc = name_in;
main_found = 0;
k = 0;
while(*pc) {
  if(*pc == '-' || *pc == ',') {
   main_found = 0;
  } else if(!main_found && isupper(*pc)) {
   name_out[k++] = *pc;
   main_found = 1;
   }
  pc++;
  }
name_out[k] = '\0';

return(0);
}
/************************************************************************
*
*************************************************************************/
static int get_coordinates_from_WDS_name(char *WDS_name, 
                                         double *right_ascension0, 
                                         double *declination0)
{
int nval, i0, i1, i2, i3, sign;
double ww;
char *pc, buffer[80];

*right_ascension0 = 0.;
*declination0 = 0.;

strcpy(buffer, WDS_name);
nval = sscanf(buffer, "%02d%02d%02d", &i0, &i1, &i2);
if(nval != 3) {
  fprintf(stderr, "get_coordinates_from_WDS_name/Error reading %s\n", WDS_name);
  return(-1);
  }
printf("%d %d %d \n", i0, i1, i2);

pc = buffer;
while(isdigit(*pc)) pc++;
if(*pc == '+') sign = 1;
else if(*pc == '-') sign = -1; 
else {
  fprintf(stderr, "get_coordinates_from_WDS_name/Error reading %s\n", WDS_name);
  return(-1);
  }
if(*pc) pc++;
strcpy(buffer, pc);

return(0);
}
