/************************************************************************
* "gpo_rw_OC6.h"
* Set of routines used to read OC6 catalog
*
* JLP 
* Version 21/06/2015
*************************************************************************/
#ifndef _gpo_rw_OC6_h /* BOF sentry */
#define _gpo_rw_OC6_h

/* Declaring linkage specification to have "correct names"
* that can be linked with C programs */

int get_orbit_from_OC6_list(char *in_line, int iline, int is_master_file,
                       char *WDS_name, char *ADS_name, char *discov_name,
                       char *comp_name, char *object_name, char *author,
                       float *Omega_node, float *omega_peri, float *i_incl, 
                       float *e_eccent, float *T_periastron, float *Period, 
                       float *a_smaxis, float *mean_motion, 
                       float *orbit_equinox);
int line_extraction_from_OC6_catalog(char *OC6_fname, int is_master_file,
                                     char *ads_name, char *discov_name,
                                     char *comp_name, FILE *fp_out, int *found, 
                                     int *candidate_found, 
                                     int norbits_per_object);
int get_name_from_OC6_line(char *in_line, char *OC6_object_name,
                           char *OC6_discov_name, char *OC6_comp_name);
int really_compact_companion(char *name_in, char *name_out, int length);

#endif
