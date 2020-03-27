/************************************************************************
* "jlp_string.cpp"
* String utilities
*
* JLP 
* Version 05/07/2015
*************************************************************************/
#include <stdio.h>
#include <stdlib.h>    // exit()
#include <string.h>    // strcpy()
#include <ctype.h>                   /* isprint... */

/* The prototypes of routines included here
*/
#include "jlp_string.h"

/************************************************************************
* Removes the trailing blanks 
* and the extra-blanks (more than one successive blanks) contained in string
*
************************************************************************/
int trim_string(char *str1, int len1)
{
char str1_copy[256];
char *pc1;
int i;

if(len1 > 256) {
  fprintf(stderr,"trim_string/Fatal error; string is too large: len=%d\n",
          len1);
  exit(-1);
  }
str1[len1 - 1] = '\0';

/* Remove the non printable characters (tab, eof line, cr, etc): */
pc1 = str1; 
i = 0;
while(*pc1) {
 if(isprint(*pc1)) str1_copy[i++] = *pc1;
 pc1++;
 }
str1_copy[i] = '\0';


/* Remove the first blanks: */
pc1 = str1_copy; 
while(*pc1 == ' ') pc1++;
strcpy(str1, pc1);

/* Remove the trailing blanks: */
pc1 = str1; 
i = 0;
while(*pc1) {pc1++; i++;}
i--;
while(i > 0 && str1[i] == ' ') str1[i--] = '\0';
strcpy(str1_copy, str1);
strcpy(str1, str1_copy);

/* Remove the contiguous blanks contained in string : */
pc1 = str1_copy;
for(i = 0; i < len1 && str1[i] != '\0'; i++)
  if(str1[i] != ' ' || str1_copy[i+1] != ' ') *(pc1++) = str1[i];
*pc1 = '\0';

strcpy(str1, str1_copy);

return (0);
}
/************************************************************************
* Removes all the blanks contained in string
*
************************************************************************/
int compact_string(char *str1, int len1)
{
char str1_copy[256];
char *pc1;
int i;

if(len1 > 256) {
  fprintf(stderr,"compact_string/Fatal error; string is too large: len=%d (max=256)\n",
          len1);
  exit(-1);
  }
str1[len1 - 1] = '\0';

/* Copy the printable characters, but not the blanks: */
pc1 = str1; 
i = 0;
while(*pc1) {
 if(isprint(*pc1) && (*pc1 != ' ')) str1_copy[i++] = *pc1;
 pc1++;
 }
str1_copy[i] = '\0';

strcpy(str1, str1_copy);

return(0);
}
/************************************************************************
* Removes the '\l' '\n' and all non-printable characters in a string
*
************************************************************************/
int cleanup_string(char *str1, int len1)
{
char str1_copy[256];
char *pc1;
int i;

if(len1 > 256) {
  fprintf(stderr,"cleanup_string/Fatal error; string is too large: len=%d\n",
          len1);
  exit(-1);
  }
str1[len1 - 1] = '\0';

/* Only copy the printable characters ( not the tab, eof line, cr, etc): */
pc1 = str1; 
i = 0;
while(*pc1) {
 if(isprint(*pc1)) str1_copy[i++] = *pc1;
 pc1++;
 }
str1_copy[i] = '\0';

strcpy(str1, str1_copy);

return(0);
}
