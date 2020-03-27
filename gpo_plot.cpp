/****************************************************************************
* Name: gpo_plot.cpp 
* 
* JLP
* Version 25/06/2015
****************************************************************************/
#include "gpo_frame.h"
#include "gpo_defs.h"
#include "jlp_kepler.h"

/********************************************************************
* Plot rho versus epoch 
********************************************************************/
int GpFrame::plot_rho_vs_epoch()
{
int i, k, xgrid_is_wanted, ygrid_is_wanted, jlp_axes_are_wanted;
char xlabel[40], ylabel[40], title[80];
char nchar_type[4], pcolor[32], plot_fname[128];
double step, theta_c, rho_c, E_anom;
int orbit_type_from_elements0;

if(initialized != 1234) return(1);

// Clear display:
  m_GraphicPanel->wxGP_ClearDrawingDisplay();

if(nmeas1 == 0){ 
  fprintf(stderr, "plot_rho_vs_epoch/Error: nmeas1=0 !\n");
  return(1);
  }

k = 0;
  for(i = 0; i < nmeas1; i++) {
/* Discard all negative rho's */
    if(rho1[i] > 0) {
      xplot1[k] = epoch1[i];
      yplot1[k] = rho1[i]; 
      k++;
      }
    }
nplot1 = k;

if(nplot1 == 0) {
  fprintf(stderr, "plot_rho_vs_epoch/all rho values are negative !\n");
  return(1);
  }

// Load arrays (reset_first=1):
  strcpy(nchar_type, "82");
  strcpy(pcolor, "Black");
  m_GraphicPanel->wxGP_LoadPlotData(xplot1, yplot1, errorx1, errory1,
                                    nplot1, nchar_type, pcolor, plot_fname, 1);

// Setting orbit_type_from_elements to zero is equivalent 
// to saying no orbit has been entered yet
  orbit_type_from_elements0 = m_jlp_orbit1->GetOrbitTypeFromElements();

  if(orbit_type_from_elements0 != 0) {
//**************************************************************************
// Computing the C curve with the orbit:
// Calling Kepler_Init_Internal: initialize elements and functions of elements
  m_jlp_orbit1->Kepler_Init_Internal();

  nplot1 = 512;
  step = (epoch1[nmeas1 - 1] - epoch1[0] )/ (double)(nplot1 - 1); 
  for(i = 0; i < nplot1; i++) {
    xplot1[i] = epoch1[0] + step * (double)i;
// Compute ephemerids (theta in degrees, rho in arcseconds
    Kepler_Ephemerid1(xplot1[i], &theta_c, &rho_c, &E_anom);
    yplot1[i] = rho_c;
    }

// Load computed arrays (reset_first=0):
  strcpy(nchar_type, "L0");
  strcpy(pcolor, "Black");
  m_GraphicPanel->wxGP_LoadPlotData(xplot1, yplot1, errorx1, errory1,
                                    nplot1, nchar_type, pcolor, plot_fname, 0);
//**************************************************************************
  }

// no grid, but JLP_axes:
  strcpy(xlabel, "Epoch [year]");
  strcpy(ylabel, "Rho [arcsec]");
  strcpy(title, "");

// JLP2015: grid does not seem to be working...
  xgrid_is_wanted = 0;
  ygrid_is_wanted = 0;
// JLP axes are better for rho and theta
  jlp_axes_are_wanted = 1;

// iplan=0 x1=0 x2=0 y1=0 y2=0
  m_GraphicPanel->wxGP_LoadPlotSettings(xlabel, ylabel, title, xgrid_is_wanted, 
                                   ygrid_is_wanted, jlp_axes_are_wanted,
                                   0, 0, 0, 0, 0);

// Call Newplot():
  m_GraphicPanel->wxGP_PlotToDrawingDisplay();

return(0);
}
/********************************************************************
* Plot rho residuals versus epoch 
********************************************************************/
int GpFrame::plot_rho_resid_vs_epoch()
{
int i, k, xgrid_is_wanted, ygrid_is_wanted, jlp_axes_are_wanted;
char xlabel[40], ylabel[40], title[80];
char nchar_type[4], pcolor[32];
double theta_c, rho_c, E_anom;
int orbit_type_from_elements0;

if(initialized != 1234 || m_jlp_orbit1 == NULL) return(1);

// Clear display:
  m_GraphicPanel->wxGP_ClearDrawingDisplay();

orbit_type_from_elements0 = m_jlp_orbit1->GetOrbitTypeFromElements();

if((nmeas1 == 0) || (orbit_type_from_elements0 == 0) 
    || (orbit_type_from_elements0 > 3)){ 
  fprintf(stderr, "plot_rho_vs_epoch/Error: nmeas1=%d orbit_type=%d !\n",
          nmeas1, orbit_type_from_elements0);
  return(1);
  }

// Calling Kepler_Init_Internal: initialize elements and functions of elements
m_jlp_orbit1->Kepler_Init_Internal();

k = 0;
  for(i = 0; i < nmeas1; i++) {
/* Discard all negative rho's */
    if(rho1[i] > 0) {
      xplot1[k] = epoch1[i];
// Compute ephemerids (theta in degrees, rho in arcseconds
      Kepler_Ephemerid1(epoch1[i], &theta_c, &rho_c, &E_anom);
      yplot1[k] = rho1[i] - rho_c; 
      k++;
      }
    }
nplot1 = k;

if(nplot1 == 0) {
  fprintf(stderr, "plot_rho_resid_vs_epoch/all rho values are negative !\n");
  return(1);
  }

// Load arrays (reset_first=1):
  strcpy(nchar_type, "82");
  strcpy(pcolor, "Black");
  m_GraphicPanel->wxGP_LoadPlotData(xplot1, yplot1, errorx1, errory1,
                                    nplot1, nchar_type, pcolor, plot_fname, 1);

// no grid, but JLP_axes:
  strcpy(xlabel, "Epoch [year]");
  strcpy(ylabel, "Rho(O-C) [arcsec]");
  strcpy(title, "");

// JLP2015: grid does not seem to be working...
  xgrid_is_wanted = 0;
  ygrid_is_wanted = 0;
// JLP axes are better for rho and theta
  jlp_axes_are_wanted = 1;

// iplan=0 x1=0 x2=0 y1=0 y2=0
  m_GraphicPanel->wxGP_LoadPlotSettings(xlabel, ylabel, title, xgrid_is_wanted, 
                                   ygrid_is_wanted, jlp_axes_are_wanted,
                                   0, 0, 0, 0, 0);

// Call Newplot():
  m_GraphicPanel->wxGP_PlotToDrawingDisplay();

return(0);
}
/********************************************************************
* Plot theta versus epoch 
********************************************************************/
int GpFrame::plot_theta_resid_vs_epoch()
{
double theta_c, rho_c, E_anom, wtheta;
int i, xgrid_is_wanted, ygrid_is_wanted, jlp_axes_are_wanted;
char xlabel[40], ylabel[40], title[80];
char nchar_type[4], pcolor[32];
int orbit_type_from_elements0;

if(initialized != 1234 || m_jlp_orbit1 == NULL) return(1);

// Clear display:
  m_GraphicPanel->wxGP_ClearDrawingDisplay();

orbit_type_from_elements0 = m_jlp_orbit1->GetOrbitTypeFromElements();

if((nmeas1 == 0) || (orbit_type_from_elements0 == 0) 
    || (orbit_type_from_elements0 > 3)){
  fprintf(stderr, "plot_theta_vs_epoch/Error: nmeas1=%d orbit_type=%d !\n",
          nmeas1, orbit_type_from_elements0);
  return(1);
  }

// Calling Kepler_Init_Internal: initialize elements and functions of elements
  m_jlp_orbit1->Kepler_Init_Internal();

nplot1 = nmeas1;
  for(i = 0; i < nplot1; i++) {
    xplot1[i] = epoch1[i];
// Compute ephemerids (theta in degrees, rho in arcseconds
    Kepler_Ephemerid1(epoch1[i], &theta_c, &rho_c, &E_anom);
    wtheta = theta1[i] - theta_c;
// Will set residuals in the -90, +90 degrees range
    if(wtheta < -90.) wtheta += 360.;
// Adjust quadrant
    while(ABS(wtheta) > 90.) {
      if(wtheta > 0) wtheta -= 180.;
      else wtheta += 180.;
      }
    yplot1[i] = wtheta; 
    }

// Load data arrays (reset_first=1):
  strcpy(nchar_type, "82");
  strcpy(pcolor, "Black");
  m_GraphicPanel->wxGP_LoadPlotData(xplot1, yplot1, errorx1, errory1,
                                    nplot1, nchar_type, pcolor, plot_fname, 1);

// no grid, but JLP_axes:
  strcpy(xlabel, "Epoch [year]");
  strcpy(ylabel, "Theta(O-C) [degr]");
  strcpy(title, "");

// JLP2015: grid does not seem to be working...
  xgrid_is_wanted = 0;
  ygrid_is_wanted = 0;
// JLP axes are better for rho and theta
  jlp_axes_are_wanted = 1;

// iplan=0 x1=0 x2=0 y1=0 y2=0
  m_GraphicPanel->wxGP_LoadPlotSettings(xlabel, ylabel, title, xgrid_is_wanted, 
                                   ygrid_is_wanted, jlp_axes_are_wanted,
                                   0, 0, 0, 0, 0);

// Call Newplot():
  m_GraphicPanel->wxGP_PlotToDrawingDisplay();

return(0);
}
/********************************************************************
* Plot theta residuals versus epoch 
********************************************************************/
int GpFrame::plot_theta_vs_epoch()
{
double step, theta_c, rho_c, E_anom;
int i, xgrid_is_wanted, ygrid_is_wanted, jlp_axes_are_wanted;
char xlabel[40], ylabel[40], title[80];
char nchar_type[4], pcolor[32];
int orbit_type_from_elements0;

if(initialized != 1234) return(1);

// Clear display:
  m_GraphicPanel->wxGP_ClearDrawingDisplay();

if(nmeas1 == 0){ 
  fprintf(stderr, "plot_theta_vs_epoch/Error: nmeas1=0 !\n");
  return(1);
  }

nplot1 = nmeas1;
  for(i = 0; i < nplot1; i++) {
    xplot1[i] = epoch1[i];
    yplot1[i] = theta1[i]; 
    }

// Load data arrays (reset_first=1):
  strcpy(nchar_type, "82");
  strcpy(pcolor, "Black");
  m_GraphicPanel->wxGP_LoadPlotData(xplot1, yplot1, errorx1, errory1,
                                    nplot1, nchar_type, pcolor, plot_fname, 1);

// Setting orbit_type_from_elements to zero is equivalent
// to saying no orbit has been entered yet
  orbit_type_from_elements0 = m_jlp_orbit1->GetOrbitTypeFromElements();

  if(orbit_type_from_elements0 != 0) {
//**************************************************************************
// Computing the C curve with the orbit:
// Calling Kepler_Init_Internal: initialize elements and functions of elements
  m_jlp_orbit1->Kepler_Init_Internal();

  nplot1 = 512;
  step = (epoch1[nmeas1 - 1] - epoch1[0] )/ (double)(nplot1 - 1); 
  for(i = 0; i < nplot1; i++) {
    xplot1[i] = epoch1[0] + step * (double)i;
// Compute ephemerids (theta in degrees, rho in arcseconds
    Kepler_Ephemerid1(xplot1[i], &theta_c, &rho_c, &E_anom);
    yplot1[i] = theta_c;
    }

// Load computed arrays (reset_first=0):
  strcpy(nchar_type, "L0");
  strcpy(pcolor, "Black");
  m_GraphicPanel->wxGP_LoadPlotData(xplot1, yplot1, errorx1, errory1,
                                    nplot1, nchar_type, pcolor, plot_fname, 0);
//**************************************************************************
  }

// no grid, but JLP_axes:
  strcpy(xlabel, "Epoch [year]");
  strcpy(ylabel, "Theta [degr]");
  strcpy(title, "");

// JLP2015: grid does not seem to be working...
  xgrid_is_wanted = 0;
  ygrid_is_wanted = 0;
// JLP axes are better for rho and theta
  jlp_axes_are_wanted = 1;

// iplan=0 x1=0 x2=0 y1=0 y2=0
  m_GraphicPanel->wxGP_LoadPlotSettings(xlabel, ylabel, title, xgrid_is_wanted, 
                                   ygrid_is_wanted, jlp_axes_are_wanted,
                                   0, 0, 0, 0, 0);

// Call Newplot():
  m_GraphicPanel->wxGP_PlotToDrawingDisplay();

return(0);
}
/********************************************************************
* Plot orbit in sky plane. 
********************************************************************/
int GpFrame::plot_orbit_skyplane()
{
int i, k, xgrid_is_wanted, ygrid_is_wanted, jlp_axes_are_wanted;
double T_periastron, Period, ww;
double oelmnts0[NELEMENTS], oelmnts_err0[NELEMENTS];
double step, theta_c, rho_c, E_anom;
double rho, theta, cross_width, epoch_o;
char xlabel[40], ylabel[40], title[80];
char nchar_type[4], pcolor[32];
int orbit_type_from_elements0, orbit_type0;

 if(initialized != 1234 || m_jlp_orbit1 == NULL) {
    return(1);
   }

// Clear display:
  m_GraphicPanel->wxGP_ClearDrawingDisplay();

if(nmeas1 == 0){ 
  fprintf(stderr, "plot_orbit_skyplane/Error: nmeas1=0 !\n");
  return(1);
  }

/* Transformation with a change of origin for theta (zero at the bottom of the
* plot)
*/
  k = 0;
  for(i = 0; i < nmeas1; i++) {
    rho = rho1[i];
    theta = theta1[i];
    theta -= 90.;
/* Conversion to radians: */
    theta *= (PI/180.);
/* Discard all negative rho's */
      if(rho > 0) {
        xplot1[k] = rho * cos(theta);
        yplot1[k] = rho * sin(theta);
        k++;
      }
    } /* EOF loop on i */

nplot1 = k;
if(nplot1 == 0) {
  fprintf(stderr, "plot_robit_skyplane/all rho values are negative !\n");
  return(1);
  }

// Load arrays (reset_first=1):
// 42: +  52: x 82: o
  strcpy(nchar_type, "82");
  strcpy(pcolor, "Black");
  m_GraphicPanel->wxGP_LoadPlotData(xplot1, yplot1, errorx1, errory1,
                                    nplot1, nchar_type, pcolor, plot_fname, 1);

// Setting orbit_type_from_elements to zero is equivalent
// to saying no orbit has been entered yet
  orbit_type_from_elements0 = m_jlp_orbit1->GetOrbitTypeFromElements();

  if(orbit_type_from_elements0 != 0) {
//**************************************************************************
// Computing the orbit:
// ckeyword[NELEMENTS] = {'P', 'T', 'e', 'a', 'W', 'w', 'i', 'K1', 'K2', 'V0'};
  m_jlp_orbit1->GetOrbitalElements(oelmnts0, oelmnts_err0, N10, &orbit_type0);
  Period = oelmnts0[0];
  T_periastron = oelmnts0[1];

// Calling Kepler_Init_Internal: initialize elements and functions of elements
  m_jlp_orbit1->Kepler_Init_Internal();

  nplot1 = 512;
  step = Period / (double)(nplot1 -1);
  for(i = 0; i < nplot1; i++) {
    epoch_o = T_periastron + step * (double)i;
// Compute ephemerids (theta in degrees, rho in arcseconds
    Kepler_Ephemerid1(epoch_o, &theta_c, &rho_c, &E_anom);
/* Transformation with a change of origin for theta (zero at the bottom of the
* plot)
*/
    theta_c -= 90.;
    ww = theta_c * PI / 180.;
    xplot1[i] = rho_c * cos(ww);
    yplot1[i] = rho_c * sin(ww);
    }

// Load computed arrays (reset_first=0):
  strcpy(nchar_type, "L0");
  strcpy(pcolor, "Black");
  m_GraphicPanel->wxGP_LoadPlotData(xplot1, yplot1, errorx1, errory1,
                                    nplot1, nchar_type, pcolor, plot_fname, 0);
//**************************************************************************
  }

// Draw central cross and North-East label with width of a/10.: 
// olmnts0[3] = a 
  cross_width = oelmnts0[3]/10.;
  draw_central_cross(cross_width);
  draw_north_east_label(cross_width);

/* Draw line of absids (Epoch of periastron, and opposite at +Period/2) */
  draw_line_of_absids();

// no grid, but JLP_axes:
  strcpy(xlabel, "X [arcsec]");
  strcpy(ylabel, "Y [arcsec]");
  strcpy(title, "");

// JLP2015: grid does not seem to be working...
  xgrid_is_wanted = 0;
  ygrid_is_wanted = 0;
// JLP axes are not needed here but better for reading values 
  jlp_axes_are_wanted = 1;

// iplan=1 x1=0 x2=0 y1=0 y2=0
  m_GraphicPanel->wxGP_LoadPlotSettings(xlabel, ylabel, title, xgrid_is_wanted, 
                                        ygrid_is_wanted, jlp_axes_are_wanted,
                                        1, 0, 0, 0, 0);

// Call Newplot():
  m_GraphicPanel->wxGP_PlotToDrawingDisplay();

return(0);
}
/******************************************************************
* Draw line of absids
*
********************************************************************/
int GpFrame::draw_line_of_absids()
{
double oelmnts0[NELEMENTS], oelmnts_err0[NELEMENTS];
double Period, T_periastron, ww;
float epoch_o, theta_c, rho_c;
char nchar_type[4], pcolor[32];
int orbit_type0;

// ckeyword[NELEMENTS] = {'P', 'T', 'e', 'a', 'W', 'w', 'i', 'K1', 'K2', 'V0'};
  m_jlp_orbit1->GetOrbitalElements(oelmnts0, oelmnts_err0, N10, &orbit_type0);
  Period = oelmnts0[0];
  T_periastron = oelmnts0[1];

  epoch_o = T_periastron;
  m_jlp_orbit1->ComputeEphemerid(epoch_o, &theta_c, &rho_c);
  theta_c -= 90.;
  ww = theta_c * PI / 180.;
  xplot1[0] = rho_c * cos(ww);
  yplot1[0] = rho_c * sin(ww);

  epoch_o = T_periastron + Period / 2.;
  m_jlp_orbit1->ComputeEphemerid(epoch_o, &theta_c, &rho_c);
  theta_c -= 90.;
  ww = theta_c * PI / 180.;
  xplot1[1] = rho_c * cos(ww);
  yplot1[1] = rho_c * sin(ww);
  nplot1 = 2;

// Load computed arrays (reset_first=0):
// Gray dash line
  strcpy(nchar_type, "L1");
  strcpy(pcolor, "Gray");
  m_GraphicPanel->wxGP_LoadPlotData(xplot1, yplot1, errorx1, errory1,
                                    nplot1, nchar_type, pcolor, plot_fname, 0);

return(0);
}
/******************************************************************
* Draw central cross 
*
* cross_width : in user coordinates
********************************************************************/
int GpFrame::draw_central_cross(double cross_width)
{
double x1, y1, x2, y2;
int lwidth = 0;
char pcolor[64];

strcpy(pcolor, "Black");

/* User coordinates (arcseconds with orbit center at 0,0) */
  x1 = -cross_width/2.; 
  x2 = cross_width/2.; 
  y1 = 0.; 
  y2 = 0.;
  m_GraphicPanel->wxGP_PlotLine1(x1, y1, x2, y2, lwidth, pcolor); 
  x1 = 0.; 
  x2 = 0.; 
  y1 = -cross_width/2.; 
  y2 = cross_width/2.;
  m_GraphicPanel->wxGP_PlotLine1(x1, y1, x2, y2, lwidth, pcolor); 

return(0);
}
/******************************************************************
* Determine the sense of motion of the companion
*
********************************************************************/
int GpFrame::compute_sense_of_motion(int *north_to_east)
{
double oelmnts0[NELEMENTS], oelmnts_err0[NELEMENTS];
double T_periastron;
float epoch_o1, epoch_o2, theta_c1, theta_c2, rho_c1, rho_c2;
int orbit_type0;

// ckeyword[NELEMENTS] = {'P', 'T', 'e', 'a', 'W', 'w', 'i', 'K1', 'K2', 'V0'};
  m_jlp_orbit1->GetOrbitalElements(oelmnts0, oelmnts_err0, N10, &orbit_type0);
  T_periastron = oelmnts0[1];
  epoch_o1 = T_periastron;
  m_jlp_orbit1->ComputeEphemerid(epoch_o1, &theta_c1, &rho_c1);

// Small difference of epochs to have a small difference in the angles too..
  epoch_o2 = T_periastron + 0.001;
  m_jlp_orbit1->ComputeEphemerid(epoch_o2, &theta_c2, &rho_c2);

// Look for the sense of motion of the companion:
  if(theta_c2 > theta_c1)
    *north_to_east = 1;
  else
    *north_to_east = 0;
return(0);
}
/******************************************************************
* Draw North-East label: 
*
* INPUT:
* cross_width (in arcseconds)
********************************************************************/
int GpFrame::draw_north_east_label(const double cross_width)
{
char pcolor[32];
double x1, y1, x2, y2;
int lwidth = 0, north_to_east;

// Determine the sense of motion
 compute_sense_of_motion(&north_to_east);

/* Draw central cross (to do the same as "2binpl.for"): */
/* User coordinates (arcseconds) */

  strcpy(pcolor, "Black");
  nplot1 = 2;

  x1 = -cross_width/2.; 
  y1 = 0.; 
  x2 = cross_width/2.;
  y2 = 0.; 
  m_GraphicPanel->wxGP_PlotLine1(x1, y1, x2, y2, lwidth, pcolor); 

  x1 = 0.; 
  y1 = -cross_width/2.; 
  x2 = 0.; 
  y2 = cross_width/2.;
  m_GraphicPanel->wxGP_PlotLine1(x1, y1, x2, y2, lwidth, pcolor); 

#ifdef TTTT
INT4 ix1, ix2, iy1, iy2, max_length, idrawit;
float x1, x2, y1, y2;
float angle, expand, length;
char xlabel[20];

/* MGO coordinates */
ix1 = 25000; iy1 = 9000;
JLP_RELOC(&ix1, &iy1, &idv);
ix2 = ix1 + 2500;
iy2 = iy1;
JLP_DRAW(&ix2, &iy2, &idv);
/*
void JLP_SPLABEL(char *xlabel, INT4 *max_length, INT4 *ix, INT4 *iy,
                 float *angle, float *expand, INT4 *idrawit, float *length,
                 INT4 *idv1)
*/
max_length = 1;
strcpy(xlabel, "E");
ix2 += 50;
iy2 -= 300;
angle = 0.;
expand = 1.2;
idrawit = 1;
JLP_SPLABEL(xlabel, &max_length, &ix2, &iy2, &angle, &expand, &idrawit,
            &length, &idv);

JLP_RELOC(&ix1, &iy1, &idv);
ix2 = ix1;
iy2 = iy1 - 2500;
JLP_DRAW(&ix2, &iy2, &idv);
ix2 -= 500;
iy2 -= 800;
strcpy(xlabel, "N");
JLP_SPLABEL(xlabel, &max_length, &ix2, &iy2, &angle, &expand, &idrawit,
            &length, &idv);

#endif
return(0);
}
/********************************************************************
* Plot radial velocities versus epoch 
********************************************************************/
int GpFrame::plot_rv_vs_epoch()
{
int i, xgrid_is_wanted, ygrid_is_wanted, jlp_axes_are_wanted;
char xlabel[40], ylabel[40], title[80];
char nchar_type[4], pcolor[32];
double rv1_c, rv2_c, E_anom;
double step, epoch_start, epoch_end;
int orbit_type_from_elements0;

if(initialized != 1234) return(1);

// Clear display:
  m_GraphicPanel->wxGP_ClearDrawingDisplay();

if(nrv1 == 0 && nrv2 == 0){
  fprintf(stderr, "plot_rv_vs_epoch/Error: nrv1=0 and nrv2=0 !\n");
  return(1);
  }

epoch_start = 1e+6;
epoch_end = -1;

if(nrv1 != 0) {
  nplot1 = nrv1;
  for(i = 0; i < nrv1; i++) {
    xplot1[i] = epoch_rv1[i];
    epoch_start = MINI(epoch_start, xplot1[i]);
    epoch_end = MAXI(epoch_end, xplot1[i]);
    yplot1[i] = rv1[i]; 
    }

// Load arrays (reset_first=1):
  strcpy(nchar_type, "92");
  strcpy(pcolor, "Black");
  m_GraphicPanel->wxGP_LoadPlotData(xplot1, yplot1, errorx1, errory1,
                                    nplot1, nchar_type, pcolor, plot_fname, 1);
}

if(nrv2 != 0) {
  nplot1 = nrv2;
  for(i = 0; i < nrv2; i++) {
    xplot1[i] = epoch_rv2[i];
    epoch_start = MINI(epoch_start, xplot1[i]);
    epoch_end = MAXI(epoch_end, xplot1[i]);
    yplot1[i] = rv2[i]; 
    }

  strcpy(nchar_type, "82");
  strcpy(pcolor, "Black");
  m_GraphicPanel->wxGP_LoadPlotData(xplot1, yplot1, errorx1, errory1,
                                    nplot1, nchar_type, pcolor, plot_fname, 0);
}

// Setting orbit_type_from_elements to zero is equivalent
// to saying no orbit has been entered yet
  orbit_type_from_elements0 = m_jlp_orbit1->GetOrbitTypeFromElements();

  if(orbit_type_from_elements0 != 0) {
//**************************************************************************
// Computing the C curve with the orbit:
// Calling Kepler_Init_Internal: initialize elements and functions of elements
  m_jlp_orbit1->Kepler_Init_Internal();

  nplot1 = 512;
  step = (epoch_end - epoch_start)/ (double)(nplot1 - 1); 
  for(i = 0; i < nplot1; i++) {
    xplot1[i] = epoch_start + step * (double)i;
// Compute ephemerids for a given epoch (rv1, rv2 in km/s)
    Kepler_Ephemerid2(xplot1[i], &rv1_c, &rv2_c, &E_anom);
    yplot1[i] = rv1_c;
    yplot2[i] = rv2_c;
    }

// Load computed arrays (reset_first=0):
  strcpy(nchar_type, "L0");
  strcpy(pcolor, "Black");
  m_GraphicPanel->wxGP_LoadPlotData(xplot1, yplot1, errorx1, errory1,
                                    nplot1, nchar_type, pcolor, plot_fname, 0);
  strcpy(nchar_type, "L1");
  strcpy(pcolor, "Black");
  m_GraphicPanel->wxGP_LoadPlotData(xplot1, yplot1, errorx1, errory1,
                                    nplot1, nchar_type, pcolor, plot_fname, 0);
//****************************************************************************
  }

// no grid, but JLP_axes:
  strcpy(xlabel, "Epoch [year]");
  strcpy(ylabel, "Radial velocity [km/s]");
  strcpy(title, "");

// JLP2015: grid does not seem to be working...
  xgrid_is_wanted = 0;
  ygrid_is_wanted = 0;
// JLP axes are better for rho and theta
  jlp_axes_are_wanted = 1;

// iplan=0 x1=0 x2=0 y1=0 y2=0
  m_GraphicPanel->wxGP_LoadPlotSettings(xlabel, ylabel, title, xgrid_is_wanted, 
                                   ygrid_is_wanted, jlp_axes_are_wanted,
                                   0, 0, 0, 0, 0);

// Call Newplot():
  m_GraphicPanel->wxGP_PlotToDrawingDisplay();

return(0);
}
/********************************************************************
* Plot residuals of radial velocities versus epoch 
********************************************************************/
int GpFrame::plot_rv_resid_vs_epoch()
{
int i, xgrid_is_wanted, ygrid_is_wanted, jlp_axes_are_wanted;
char xlabel[40], ylabel[40], title[80];
char nchar_type[4], pcolor[32];
double rv1_c, rv2_c, E_anom;
int orbit_type_from_elements0;

if(initialized != 1234 || m_jlp_orbit1 == NULL) return(1);

// Clear display:
  m_GraphicPanel->wxGP_ClearDrawingDisplay();

orbit_type_from_elements0 = m_jlp_orbit1->GetOrbitTypeFromElements();

if((nrv1 == 0) && (nrv2 == 0) || (orbit_type_from_elements0 == 0) 
   || (orbit_type_from_elements0 == 3)){
  fprintf(stderr, "plot_rv_resid_vs_epoch/Error: nrv1=%d nrv2=%d orbit_type=%d !\n",
          nrv1, nrv2, orbit_type_from_elements0);
  return(1);
  }

// Calling Kepler_Init_Internal: initialize elements and functions of elements
m_jlp_orbit1->Kepler_Init_Internal();

if(nrv1 != 0) {
  nplot1 = nrv1;
  for(i = 0; i < nrv1; i++) {
    xplot1[i] = epoch_rv1[i];
// Compute ephemerids for a given epoch (rv1, rv2 in km/s)
    Kepler_Ephemerid2(epoch_rv1[i], &rv1_c, &rv2_c, &E_anom);
    yplot1[i] = rv1[i] - rv1_c; 
    }

// Load arrays (reset_first=1):
  strcpy(nchar_type, "92");
  strcpy(pcolor, "Black");
  m_GraphicPanel->wxGP_LoadPlotData(xplot1, yplot1, errorx1, errory1,
                                    nplot1, nchar_type, pcolor, plot_fname, 1);
}

if(nrv2 != 0) {
  nplot1 = nrv2;
  for(i = 0; i < nrv2; i++) {
    xplot1[i] = epoch_rv2[i];
// Compute ephemerids for a given epoch (rv1, rv2 in km/s)
    Kepler_Ephemerid2(epoch_rv2[i], &rv1_c, &rv2_c, &E_anom);
    yplot1[i] = rv2[i] - rv2_c; 
    }

  strcpy(nchar_type, "82");
  strcpy(pcolor, "Black");
  m_GraphicPanel->wxGP_LoadPlotData(xplot1, yplot1, errorx1, errory1,
                                    nplot1, nchar_type, pcolor, plot_fname, 0);
}

// no grid, but JLP_axes:
  strcpy(xlabel, "Epoch [year]");
  strcpy(ylabel, "Radial velocity(O-C) [km/s]");
  strcpy(title, "");

// JLP2015: grid does not seem to be working...
  xgrid_is_wanted = 0;
  ygrid_is_wanted = 0;
// JLP axes are better for rho and theta
  jlp_axes_are_wanted = 1;

// iplan=0 x1=0 x2=0 y1=0 y2=0
  m_GraphicPanel->wxGP_LoadPlotSettings(xlabel, ylabel, title, xgrid_is_wanted, 
                                   ygrid_is_wanted, jlp_axes_are_wanted,
                                   0, 0, 0, 0, 0);

// Call Newplot():
  m_GraphicPanel->wxGP_PlotToDrawingDisplay();

return(0);
}
/********************************************************************
* Plot radial velocities versus phase 
********************************************************************/
int GpFrame::plot_rv_vs_phase()
{
double oelmnts0[NELEMENTS], oelmnts_err0[NELEMENTS];
int i, xgrid_is_wanted, ygrid_is_wanted, jlp_axes_are_wanted;
char xlabel[40], ylabel[40], title[80];
char nchar_type[4], pcolor[32];
double rv1_c, rv2_c, E_anom;
double step, T0, Period, epoch_o, dt, ww;
int orbit_type0, orbit_type_from_elements0;

if(initialized != 1234 || m_jlp_orbit1 == NULL) return(1);

// Clear display:
  m_GraphicPanel->wxGP_ClearDrawingDisplay();

orbit_type_from_elements0 = m_jlp_orbit1->GetOrbitTypeFromElements();

if((nrv1 == 0) && (nrv2 == 0) || (orbit_type_from_elements0 == 0) 
   || (orbit_type_from_elements0 == 3)){
  fprintf(stderr, "plot_rv_vs_phase/Error: nrv1=%d nrv2=%d orbit_type=%d !\n",
          nrv1, nrv2, orbit_type_from_elements0);
  return(1);
  }

 m_jlp_orbit1->GetOrbitalElements(oelmnts0, oelmnts_err0, N10, &orbit_type0);
 Period = oelmnts0[0];
 T0 = oelmnts0[1];

if(nrv1 != 0) {
  nplot1 = nrv1;
  for(i = 0; i < nrv1; i++) {
// Compute the phase:
    ww = (epoch_rv1[i] - T0) / Period;
    xplot1[i] = ww - (int)ww;
    if(xplot1[i] < 0.) xplot1[i] += 1.;
    yplot1[i] = rv1[i]; 
    }

// Load arrays (reset_first=1):
  strcpy(nchar_type, "92");
  strcpy(pcolor, "Black");
  m_GraphicPanel->wxGP_LoadPlotData(xplot1, yplot1, errorx1, errory1,
                                    nplot1, nchar_type, pcolor, plot_fname, 1);
}

if(nrv2 != 0) {
  nplot1 = nrv2;
  for(i = 0; i < nrv2; i++) {
// Compute the phase:
    ww = (epoch_rv2[i] - T0) / Period;
    xplot1[i] = ww - (int)ww;
    if(xplot1[i] < 0.) xplot1[i] += 1.;
    yplot1[i] = rv2[i]; 
    }

  strcpy(nchar_type, "82");
  strcpy(pcolor, "Black");
  m_GraphicPanel->wxGP_LoadPlotData(xplot1, yplot1, errorx1, errory1,
                                    nplot1, nchar_type, pcolor, plot_fname, 0);
}

// Computing the C curve with the orbit:
// Calling Kepler_Init_Internal: initialize elements and functions of elements
  m_jlp_orbit1->Kepler_Init_Internal();

  nplot1 = 512;
  step = 1.0/ (double)(nplot1 - 1); 
  for(i = 0; i < nplot1; i++) {
    xplot1[i] = step * (double)i;
    epoch_o = T0 + xplot1[i] * Period; 
// Compute ephemerids for a given epoch (rv1, rv2 in km/s)
    Kepler_Ephemerid2(epoch_o, &rv1_c, &rv2_c, &E_anom);
    yplot1[i] = rv1_c;
    yplot2[i] = rv2_c;
    }

// Load computed arrays (reset_first=0):
  strcpy(nchar_type, "L0");
  strcpy(pcolor, "Black");
  m_GraphicPanel->wxGP_LoadPlotData(xplot1, yplot1, errorx1, errory1,
                                    nplot1, nchar_type, pcolor, plot_fname, 0);
  strcpy(nchar_type, "L1");
  strcpy(pcolor, "Black");
  m_GraphicPanel->wxGP_LoadPlotData(xplot1, yplot1, errorx1, errory1,
                                    nplot1, nchar_type, pcolor, plot_fname, 0);

// no grid, but JLP_axes:
  strcpy(xlabel, "Phase");
  strcpy(ylabel, "Radial velocity [km/s]");
  strcpy(title, "");

// JLP2015: grid does not seem to be working...
  xgrid_is_wanted = 0;
  ygrid_is_wanted = 0;
// JLP axes are not needed here but better for reading values 
  jlp_axes_are_wanted = 1;

// iplan=0 x1=0 x2=1 y1=0 y2=0
  m_GraphicPanel->wxGP_LoadPlotSettings(xlabel, ylabel, title, xgrid_is_wanted, 
                                   ygrid_is_wanted, jlp_axes_are_wanted,
                                   0, 0., 1., 0, 0);

// Call Newplot():
  m_GraphicPanel->wxGP_PlotToDrawingDisplay();

return(0);
}
