/****************************************************************************
* Name: gpo_defs.h
* 
* JLP
* Version 08/06/2015
****************************************************************************/
#ifndef _gpo_defs__ 
#define _gpo_defs__ 

// Maximum number of points (for curves) 
#define NPLOT_MAX 4096

// Maximum number of observations:
#define NMEAS_MAX 4096

// Maximum number of results (in gpo_widget_panel.cpp)
#define NRESULTS_MAX 20 

// Panel text types
#define NTYPES_TEXTPANEL 5

// Language: English, French, Italian, Spanish, German
#define NLANG 5
#define NMAX_MESSAGES 1024 

// Number of elements:
#define NELEMENTS 10 
#define N7 7 
#define N10 10 

#ifndef ABS
#define ABS(a) ((a) < 0.0  ? (-(a)) : (a))
#endif

#ifndef SQUARE
#define SQUARE(a) ((a) * (a))
#endif

#ifndef PI
#define PI 3.14159265
#endif

#ifndef MINI
#define MINI(a,b) ((a) < (b) ? (a) : (b))
#endif

#ifndef MAXI
#define MAXI(a,b) ((a) < (b) ? (b) : (a))
#endif

#define YEAR_TO_DAYS 365.242198781

#endif
