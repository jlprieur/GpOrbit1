/*************************************************************************
* gpo_rw_orbitfiles.cpp
* To read/write visual orbit files
*
* Contains:
*
* JLP
* Version 28/06/2015
**************************************************************************/
#include <stdio.h> 
#include <stdlib.h>    // exit() 
#include <string.h>    // strcpy() 
#include "gpo_rw_files.h"
#include "gpo_frame.h"

#define DEGTORAD 3.14159/180.

/*
#define DEBUG
*/

/*****************************************************************************
* {'P', 'T', 'e', 'a', 'W', 'w', 'i', 'K1', 'K2', 'V0'};
* orbit_type = 1 visual & spectro SB1 (9 elements idx1= 1 1 1 1 1 1 1 1 0 1)
* orbit_type = 2 visual & spectro SB2 (10 elements idx1= 1 1 1 1 1 1 1 1 1 1)
* orbit_type = 3 visual (7 elements: idx2= 1 1 1 1 1 1 1 0 0 0)
* orbit_type = 4 spectroscopic SB1 (6 elements: idx3= 1 1 1 0 0 1 0 1 0 1)
* orbit_type = 5 spectroscopic SB2 (7 elements: idx3= 1 1 1 0 0 1 0 1 1 1)
*****************************************************************************/

/***************************************************************************
* rd_visual_oelements_scardia (Scardia's format)
* Read visual orbital elements in file and retrieve the orbital parameters
*
* INPUT:
* fp_in: pointer to the input file containing the measurements and the
*        orbital elements
*
* OUTPUT:
* arrays filled with the orbital parameters
// Orbital elements: P T0 e a W w i K1 K2 V0
* double orbital_elements[NELEMENTS];
*
*
***************************************************************************/
int rd_visual_oelements_scardia(const char *orbit_infile, 
                                double *orbital_elements, double *orbit_equinox,
                                char *object_name, double *right_ascension0,
                                double *declination0, int *orbit_type0)
{
FILE *fp_in;
float Omega_node, omega_peri, i_incl;
float e_eccent, T_periastron, Period;
float a_smaxis, mean_motion, equinox;
int nval, i;
char in_line[80], buffer[128], *pc;

*orbit_type0 = 0;

// Coordinates are not available yet (not necessary since Marco makes 
// precession corrections himself in data files) 
*right_ascension0 = 0.;
*declination0 = 0.;

if((fp_in = fopen(orbit_infile, "r")) == NULL) {
  fprintf(stderr, "Fatal error opening input file >%s<\n", orbit_infile);
  exit(-1);
  }

/* Read first uncommented line: */
  while(!feof(fp_in)) {
    if(fgets(in_line, 80, fp_in)) {
      if(in_line[0] != '%' && in_line[0] != '#') break;
    }
  }
fclose(fp_in);

/* Read orbital elements from input file: */
/* Marco's format */
// Format: W w i e T P a equinox
   nval = sscanf(in_line, "%f %f %f %f %f %f %f %f",
                 &Omega_node, &omega_peri, &i_incl, &e_eccent,
                 &T_periastron, &Period, &a_smaxis, &equinox);
   if(nval != 8) {
     fprintf(stderr, "read_orbital_elements_from_file/Fatal error reading input file (nval=%d) \n %s \n",
              nval, in_line);
     return(-2);
    }

#ifdef DEBUG
 printf("nval=%d Omega_node=%.3f omega_peri=%.3f incl=%.3f e=%.3f T=%.3f P=%.3f a=%.3f Equinox=%.3f\n",
        nval, Omega_node, omega_peri, i_incl, e_eccent,
        T_periastron, Period, a_smaxis, equinox);
#endif

  mean_motion = 360.0 / Period;

// Orbital elements: P T0 e a W w i K1 K2 V0
 orbital_elements[0] = Period;
 orbital_elements[1] = T_periastron;
 orbital_elements[2] = e_eccent;
 orbital_elements[3] = a_smaxis;
 orbital_elements[4] = Omega_node;
 orbital_elements[5] = omega_peri;
 orbital_elements[6] = i_incl;
 *orbit_equinox = equinox;
 *orbit_type0 = 3;

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

strcpy(object_name, pc);

return(0);
}
/*****************************************************************************
*
* {'P', 'T', 'e', 'a', 'W', 'w', 'i', 'K1', 'K2', 'V0'};
* orbit_type = 1 visual & spectro SB1 (9 elements idx1= 1 1 1 1 1 1 1 1 0 1)
* orbit_type = 2 visual & spectro SB2 (10 elements idx1= 1 1 1 1 1 1 1 1 1 1)
* orbit_type = 3 visual (7 elements: idx2= 1 1 1 1 1 1 1 0 0 0)
* orbit_type = 4 spectroscopic SB1 (6 elements: idx3= 1 1 1 0 0 1 0 1 0 1)
* orbit_type = 5 spectroscopic SB2 (7 elements: idx3= 1 1 1 0 0 1 0 1 1 1)
*
*****************************************************************************/
void set_idx_free_flags(int idx_free[5][N10])
{
int i, j;
int idx_free_flags[5][N10] = { {1, 1, 1, 1, 1, 1, 1, 1, 0, 1},
                               {1, 1, 1, 1, 1, 1, 1, 1, 1, 1},
                               {1, 1, 1, 1, 1, 1, 1, 0, 0, 0},
                               {1, 1, 1, 0, 0, 1, 0, 1, 0, 1},
                               {1, 1, 1, 0, 0, 1, 0, 1, 1, 1} };
for(i = 0; i < 5; i++)
  for(j = 0; j < N10; j++)
   idx_free[i][j] = idx_free_flags[i][j];

}
