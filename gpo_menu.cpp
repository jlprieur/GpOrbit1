/******************************************************************************
* Name:        gpo_menu (GpFrame class)
* Purpose:     handling menu events of GpFrame class
*
* Author:      JLP 
* Version:     08/06/2015 
******************************************************************************/
#include "gpo_frame.h"
#include "gpo_frame_id.h"  // Menu identifiers
#include "gpo_rw_files.h"  // rd_orbit_datafile ... 

/************************************************************************
* Plot on the graphic panel
*************************************************************************/
void GpFrame::OnPlot( wxCommandEvent &event )
{

if(initialized != 1234) return;

// Select graphic panel:
iPage = 1;
m_notebook->SetSelection(iPage);

switch(event.GetId()) {
// Clear the graphic panel
 default :
 case ID_PLOT_IDLE:
    m_GraphicPanel->wxGP_ClearDrawingDisplay();
    break;
// Display rho vs epoch
 case ID_PLOT_RHO:
    plot_rho_vs_epoch();
    break;
// Display rho residuals vs epoch
 case ID_PLOT_RHO_RESID:
    plot_rho_resid_vs_epoch();
    break;
// Display theta vs epoch
 case ID_PLOT_THETA:
    plot_theta_vs_epoch();
    break;
// Display theta residuals vs epoch
 case ID_PLOT_THETA_RESID:
    plot_theta_resid_vs_epoch();
    break;
// Display sky plane orbit 
 case ID_PLOT_ORBIT:
    plot_orbit_skyplane();
    break;
// Display radial velocity vs epoch
 case ID_PLOT_RV_EPOCH:
    plot_rv_vs_epoch();
    break;
// Display radial velocity residuals vs epoch
 case ID_PLOT_RV_RESID:
    plot_rv_resid_vs_epoch();
    break;
// Display radial velocity vs phase
 case ID_PLOT_RV_PHASE:
    plot_rv_vs_phase();
    break;
}
}
/************************************************************************
** Input ascii orbit data file in Tokovinin's format
*************************************************************************/
void GpFrame::OnLoadTokoFile( wxCommandEvent &WXUNUSED(event) )
{
wxString fname;

fname = wxFileSelector( wxT("Load Orbit and Data File (in Tokovinin's format)"),
            wxT(""), wxT(""), wxT("dat|DAT|inp"),
            wxT("Files (*.dat;*.DAT;*.inp)|*.dat;*.DAT;*.inp"));

if(!fname.empty()) LoadTokoFile(fname);

}
/************************************************************************
** Input ascii orbit data file in Tokovinin's format
*************************************************************************/
void GpFrame::LoadTokoFile(wxString &fname)
{
wxString str1, path, extension, short_fname;
double oelmnts0[NELEMENTS], oelmnts_err0[NELEMENTS], orbit_equinox0;
double epoch0[NMEAS_MAX], theta0[NMEAS_MAX], rho0[NMEAS_MAX];
double sigma_rho0[NMEAS_MAX], weight0[NMEAS_MAX];
double epoch_rv10[NMEAS_MAX], rv10[NMEAS_MAX], sigma_rv10[NMEAS_MAX];
double epoch_rv20[NMEAS_MAX], rv20[NMEAS_MAX], sigma_rv20[NMEAS_MAX];
double weight_rv10[NMEAS_MAX], weight_rv20[NMEAS_MAX];
char notes_meas0[NMEAS_MAX][128];
char notes_rv10[NMEAS_MAX][128], notes_rv20[NMEAS_MAX][128];
char object_name0[128]; 
double right_ascension0, declination0;
char in_orbit_datafile[256];
int status_orbit, status_data, status, i, orbit_type0;
int nmeas0, nrv10, nrv20;

// Clean the screen and make some initialization:
 ResetBeforeLoadingDataFile();
 ResetBeforeLoadingOrbitFile();

// Initialize elements and errors:
 for(i = 0; i < N10; i++) {
   oelmnts0[i] = 0.;
   oelmnts_err0[i] = 0.;
   }
 orbit_type0 = 0;

//**********************************************************************//
// Read orbital parameters:
//**********************************************************************//

  status_orbit = rd_oelements_toko(fname.mb_str(), oelmnts0,
                             &orbit_equinox0,
                             object_name0, &right_ascension0,
                             &declination0, &orbit_type0);
// Transfer to internal parameters;
  if(status_orbit == 0) {
    orbit_equinox2 = orbit_equinox0;

    strcpy(object_name1, object_name0);
    right_ascension1 = right_ascension0;
    declination1 = declination0;

    UpdateElementsAfterFileIsLoaded(oelmnts0, oelmnts_err0, N10, orbit_type0);
    }


//**********************************************************************//
// Read visual-measurement data
//**********************************************************************//
  status_data = read_toko_visual_datafile(fname.mb_str(), epoch0, theta0, rho0,
                                          sigma_rho0, notes_meas0, 
                                          NMEAS_MAX, &nmeas0);

// If successful, transfer data to private arrays:
  if(status_data == 0)  {
    nmeas1 = nmeas0;
    for(i = 0; i < nmeas1; i++) {
      epoch1[i] = epoch0[i];
      theta1[i] = theta0[i];
      rho1[i] = rho0[i];
      sigma_rho1[i] = sigma_rho0[i];
      strcpy(notes_meas1[i], notes_meas0[i]);
      }
    }

//**********************************************************************//
// Read radial velocity data
//**********************************************************************//
  status = read_toko_rv_datafile(fname.mb_str(), epoch_rv10, rv10,
                                 sigma_rv10, notes_rv10, 
                                 epoch_rv20, rv20, sigma_rv20, 
                                 notes_rv20, NMEAS_MAX, &nrv10, &nrv20);

// If successful, transfer data to private arrays:
  if(status == 0)  {
    nrv1 = nrv10;
    nrv2 = nrv20;
    for(i = 0; i < nrv1; i++) {
      epoch_rv1[i] = epoch_rv10[i];
      rv1[i] = rv10[i];
      sigma_rv1[i] = sigma_rv10[i];
      strcpy(notes_rv1[i], notes_rv10[i]);
      }
    for(i = 0; i < nrv2; i++) {
      epoch_rv2[i] = epoch_rv20[i];
      rv2[i] = rv20[i];
      sigma_rv2[i] = sigma_rv20[i];
      strcpy(notes_rv2[i], notes_rv20[i]);
      }
    }
   status_data *= status;

// Compute residuals and enable/disable menu items:
   UpdateDataAfterFileIsLoaded();
   UpdateMenuAfterFileIsLoaded();

//******************
// Write to logbook:
// Logbook:
  if(status_orbit == 0 || status_data == 0) {

  LogPanel->Clear();

// Display the content of the logfile on TextPanel:
  OrbitDataTextCtrl->LoadFile(fname);

// Removes the directory name (since the full path is generally too long...)
  wxFileName::SplitPath(fname, &path, &short_fname, &extension);

  if(status_data == 0) {
    Data_fname = fname; 
    short_Data_fname = short_fname + _T(".") + extension;
    str1 = _T("%% Input data file: ") + short_Data_fname + _T("\n");
    WriteToLogbook(str1, true);
    } 
  if(status_orbit == 0) {
    Orbit_fname = fname; 
    short_Orbit_fname = short_fname + _T(".") + extension;
    str1 = _T("%% Input orbit file: ") + short_Orbit_fname + _T("\n");
    WriteToLogbook(str1, true);
    }
//******************************************************************

    } // EOF status_orbit... == 0

// Select Data file in TextPanel:
 ShowDataFileButton->SetValue(true);
 TextCancelChangesButton->Show(false);
 TextValidChangesButton->Show(false);
 TextSaveToFileButton->Show(false);

return;
}
/************************************************************************
** Input ascii orbit data file in Carquillat's format
*************************************************************************/
void GpFrame::OnLoadCarquiFile( wxCommandEvent &WXUNUSED(event) )
{
wxString str1, path, extension, fname, short_fname;
double oelmnts0[NELEMENTS], oelmnts_err0[NELEMENTS];
double epoch0[NMEAS_MAX], theta0[NMEAS_MAX], rho0[NMEAS_MAX];
double sigma_rho0[NMEAS_MAX], weight0[NMEAS_MAX];
double epoch_rv10[NMEAS_MAX], rv10[NMEAS_MAX], sigma_rv10[NMEAS_MAX];
double epoch_rv20[NMEAS_MAX], rv20[NMEAS_MAX], sigma_rv20[NMEAS_MAX];
double weight_rv10[NMEAS_MAX], weight_rv20[NMEAS_MAX];
char notes_meas0[NMEAS_MAX][128];
char notes_rv10[NMEAS_MAX][128], notes_rv20[NMEAS_MAX][128];
char object_name0[128];
double right_ascension0, declination0;
char in_orbit_datafile[256];
int status_orbit, status_data, status, i, orbit_type0;
int nmeas0, nrv10, nrv20;

 fname = wxFileSelector( wxT("Load Orbit and Data File (in Carquillat's format)"));

 if(fname.empty()) return;

// Initialize elements and errors:
 for(i = 0; i < N10; i++) {
   oelmnts0[i] = 0.;
   oelmnts_err0[i] = 0.;
   }
 orbit_type0 = 0;

// Clean the screen and make some initialization:
 ResetBeforeLoadingDataFile();
 ResetBeforeLoadingOrbitFile();

 status = read_carquillat_file(fname.mb_str(), oelmnts0, &orbit_type0,
                               epoch_rv10, rv10, sigma_rv10, weight_rv10,
                               notes_rv10, epoch_rv20, rv20, sigma_rv20,
                               weight_rv20, notes_rv20, NMEAS_MAX,
                               &nrv10, &nrv20);
 if(status) return; 

// Transfer orbitl elements to internal parameters;
 orbit_equinox2 = 0.;
 right_ascension1 = 0.;
 declination1 = 0.;
 UpdateElementsAfterFileIsLoaded(oelmnts0, oelmnts_err0, N10, orbit_type0);

// Transfer data to private arrays:
 nrv1 = nrv10;
 nrv2 = nrv20;
 for(i = 0; i < nrv1; i++) {
   epoch_rv1[i] = epoch_rv10[i];
   rv1[i] = rv10[i];
   sigma_rv1[i] = sigma_rv10[i];
   strcpy(notes_rv1[i], notes_rv10[i]);
   }
 for(i = 0; i < nrv2; i++) {
   epoch_rv2[i] = epoch_rv20[i];
   rv2[i] = rv20[i];
   sigma_rv2[i] = sigma_rv20[i];
   strcpy(notes_rv2[i], notes_rv20[i]);
   }

// Compute residuals and enable/disable menu items:
  UpdateDataAfterFileIsLoaded();
  UpdateMenuAfterFileIsLoaded();

//******************
// Write to logbook:
// Logbook:
  LogPanel->Clear();

// Display the content of the logfile on TextPanel:
  OrbitDataTextCtrl->LoadFile(fname);

// Removes the directory name (since the full path is generally too long...)
  wxFileName::SplitPath(fname, &path, &short_fname, &extension);

// Object name is retrieved from filename:
  strcpy(object_name1, short_fname.mb_str());
// Update and display it in Orbit panel:
  WidgetPanel_Update_OElements();

  Data_fname = fname; 
  short_Data_fname = short_fname + _T(".") + extension;
  str1 = _T("%% Input data file: ") + short_Data_fname + _T("\n");
  WriteToLogbook(str1, true);

  Orbit_fname = fname; 
  short_Orbit_fname = short_fname + _T(".") + extension;
  str1 = _T("%% Input orbit file: ") + short_Orbit_fname + _T("\n");
  WriteToLogbook(str1, true);

// Select Data file in TextPanel:
 ShowDataFileButton->SetValue(true);
 TextCancelChangesButton->Show(false);
 TextValidChangesButton->Show(false);
 TextSaveToFileButton->Show(false);

return;
}
/************************************************************************
* Clean the screen and make some initialization:
*
*************************************************************************/
void GpFrame::ResetBeforeLoadingDataFile()
{

// Erase graphic screen:
  m_GraphicPanel->wxGP_ClearDrawingDisplay();

// Init nber of measurements and sigma residual values:
  nmeas1 = 0;
  nrv1 = 0;
  nrv2 = 0;
  precession_corrected = 0.;
  mean_sigma_rho_resid1 = 0;
  mean_sigma_theta_resid1 = 0;
  mean_sigma_rv1_resid1 = 0;
  mean_sigma_rv2_resid1 = 0;

// Initialize also m_jlp_orbit1 internal parameters: nrv1, .. orbit_type1, ...
  m_jlp_orbit1->ClearMeasurements();

// Select "idle" item in plot menu
  menuPlot->Check(ID_PLOT_IDLE, true);
  m_GraphicPanel->wxGP_ClearDrawingDisplay();

// Clean the result display on widget panel:
 WidgetPanel_Update_Results();
}
/************************************************************************
* Clean the screen and make some initialization:
*
*************************************************************************/
void GpFrame::ResetBeforeLoadingOrbitFile()
{
double oelmnts0[NELEMENTS], oelmnts_err0[NELEMENTS];
int oelmnts_fflags0[NELEMENTS]; 
int i, orbit_type0;

orbit_equinox2 = 0.;
right_ascension1 = 0.;
declination1 = 0.;

// Erase errors and elements:
 for(i = 0; i < N10; i++) {
   oelmnts0[i] = 0.;
   oelmnts_err0[i] = 0.;
   oelmnts_fflags0[i] = 0;
   }

// This flag is also used to mean that no orbit has been entered yet: 
orbit_type0 = 0;

 m_jlp_orbit1->SetOrbitalElements(oelmnts0, oelmnts_err0, N10, orbit_type0);
 m_jlp_orbit1->SetFixedFlags(oelmnts_fflags0, N10);

// Clean the elements display on widget panel:
 WidgetPanel_Update_OElements();

}
/************************************************************************
*
*************************************************************************/
void GpFrame::UpdateElementsAfterFileIsLoaded(double *oelmnts0,
                                            double *oelmnts_err0,
                                            const int max_noelmnts0,
                                            int orbit_type0)
{
int i;
int oelmnts_fflags0[NELEMENTS];

if(initialized != 1234 || m_jlp_orbit1 == NULL) return;

// Do nothing if no orbit has been loaded:
if(orbit_type0 == 0 || max_noelmnts0 != N10) {
  fprintf(stderr, "UpdateElementsAfterFileIsLoaded/Error: orbit_type=%d noelmnts=%d\n",
          orbit_type0, max_noelmnts0);
  return;
  }

// Load elements to JLP_Orbits1 object:
m_jlp_orbit1->SetOrbitalElements(oelmnts0, oelmnts_err0, N10, orbit_type0);

// Set all fixed flags to zero:
for(i = 0; i < N10; i++) oelmnts_fflags0[i] = 0;
m_jlp_orbit1->SetFixedFlags(oelmnts_fflags0, N10);

  if(precession_corrected == -1 && nmeas1 > 0) {
    CorrectForPrecession();
  }

// Hide/show parameters :
UpdateOrbitPanelAfterFileIsLoaded();

WidgetPanel_Update_OElements();

return;
}
/************************************************************************
* Enable/disable menu items:
*
*************************************************************************/
void GpFrame::UpdateMenuAfterFileIsLoaded()
{
int orbit_type0;

if(initialized != 1234 || m_jlp_orbit1 == NULL) return;

orbit_type0 = m_jlp_orbit1->GetOrbitType();

// Case when visual measurements are absent:
 if(nmeas1 == 0) {

// Disable visual plot items in plot menu
   menuPlot->Enable(ID_PLOT_RHO, false);
   menuPlot->Enable(ID_PLOT_RHO, false);
   menuPlot->Enable(ID_PLOT_THETA, false);
   menuPlot->Enable(ID_PLOT_ORBIT, false);
   menuPlot->Enable(ID_PLOT_RHO_RESID, false);
   menuPlot->Enable(ID_PLOT_THETA_RESID, false);

// Case when visual measurements are present:
   } else {

// Enable plot items in plot menu
    menuPlot->Enable(ID_PLOT_RHO, true);
    menuPlot->Enable(ID_PLOT_THETA, true);
   menuPlot->Enable(ID_PLOT_ORBIT, true);
    if(orbit_type0 != 0) {
      menuPlot->Enable(ID_PLOT_RHO_RESID, true);
      menuPlot->Enable(ID_PLOT_THETA_RESID, true);
      }
    }

// Case when rv measurements are absent:
 if(nrv1 == 0 && nrv2 == 0) {

// Disable rv plot items in plot menu
   menuPlot->Enable(ID_PLOT_RV_EPOCH, false);
   menuPlot->Enable(ID_PLOT_RV_PHASE, false);
   menuPlot->Enable(ID_PLOT_RV_RESID, false);

// Case when rv measurements are present:
   } else {

// Enable rv plot items in plot menu
    menuPlot->Enable(ID_PLOT_RV_EPOCH, true);
    menuPlot->Enable(ID_PLOT_RV_PHASE, true);
    if(orbit_type0 != 0) {
      menuPlot->Enable(ID_PLOT_RV_RESID, true);
      } else {
      menuPlot->Enable(ID_PLOT_RV_RESID, false);
      }
    }

return;
}
/************************************************************************
* Compute residuals and enable/disable menu items:
*
*************************************************************************/
void GpFrame::UpdateDataAfterFileIsLoaded()
{
double epoch_o, dtheta_precess;
int i, orbit_type_from_elements0;

// Precession correction:
orbit_type_from_elements0 = m_jlp_orbit1->GetOrbitTypeFromElements();
  if(orbit_type_from_elements0 == 0) {
    precession_corrected = -1;
  } else if(nmeas1 > 0) {
    CorrectForPrecession();
  }

// Load data to JLP_Orbits1 object:
   if(nmeas1 > 0) m_jlp_orbit1->LoadMeasurements(epoch1, theta1, rho1, 
                                                 sigma_rho1, nmeas1);
   if(nrv1 > 0 || nrv2 > 0) 
            m_jlp_orbit1->LoadRadialVelocities(epoch_rv1, epoch_rv2, rv1, rv2,
                                   sigma_rv1, sigma_rv2, nrv1, nrv2);

/* Compute new residuals in Residuals_txt 
* panel_type: 0: orbit and data, 1: orbit only, 2: data only 3: residuals 
*/
   UpdatePanelText(3);

   WidgetPanel_Update_Results();

// Hide/show useful items in data fit panel:
   UpdateDataFitPanelAfterFileIsLoaded();

return;
}
/************************************************************************
* Input visual/rv data file with visual measurements
*************************************************************************/
void GpFrame::OnLoadDataFile( wxCommandEvent &event )
{
wxString fname, str1;
int iformat = 1;

switch(event.GetId()) {
 default :
 case ID_LOAD_WDS_DATAFILE : 
   str1 = wxT("Select WDS visual data file");
   iformat = 1;
   break;
 case ID_LOAD_MINI_DATAFILE : 
   str1 = wxT("Select mini. visual data file");
   iformat = 2;
   break;
 case ID_LOAD_TOKO_DATAFILE : 
   str1 = wxT("Select visual & RV data file (Tokovinin's format)");
   iformat = 3;
   break;
 case ID_LOAD_CORAVEL_DATAFILE : 
   str1 = wxT("Select RV data file (CORAVEL format)");
   iformat = 4;
   break;
}

 fname = wxFileSelector( str1);

 if(fname.empty()) return;

 LoadDataFile(iformat, fname);

// Select Data file in TextPanel: 
 ShowDataFileButton->SetValue(true);
 TextCancelChangesButton->Show(false);
 TextValidChangesButton->Show(false);
 TextSaveToFileButton->Show(false);

}
/************************************************************************
* Input data file with visual/rv measurements
*
* Format:
*   1 = WDS data file
*   2 = mini data file
*   3 = Tokovinin's formatted  data file
*   4 = CORAVEL  data file
*************************************************************************/
void GpFrame::LoadDataFile(const int iformat, wxString &fname)
{
wxString str1, path, extension, short_fname;
char in_orbit_datafile[256];
double epoch0[NMEAS_MAX], theta0[NMEAS_MAX], rho0[NMEAS_MAX];
double sigma_rho0[NMEAS_MAX], weight0[NMEAS_MAX];
double weight_rv10[NMEAS_MAX], weight_rv20[NMEAS_MAX];
double epoch_rv10[NMEAS_MAX], rv10[NMEAS_MAX], sigma_rv10[NMEAS_MAX];
double epoch_rv20[NMEAS_MAX], rv20[NMEAS_MAX], sigma_rv20[NMEAS_MAX];
char notes_meas0[NMEAS_MAX][128]; 
char notes_rv10[NMEAS_MAX][128], notes_rv20[NMEAS_MAX][128];
int nmeas0, nrv10, nrv20;
int status_rv, status_v, i;

// Do some initialization:
 ResetBeforeLoadingDataFile();

// Read file and decode orbital parameters:

 status_v = 1;
 status_rv = 1;
 nmeas1 = 0;
 nrv1 = 0;
 nrv2 = 0;

 switch(iformat) {
   default:

// WDS format
   case 1:
     status_v = read_WDS_visual_datafile(fname.mb_str(), epoch0, theta0, rho0, 
                                       weight0, notes_meas0, NMEAS_MAX, 
                                       &nmeas0);

// Fake values for sigma since they are missing in WDS files:
// weight = rho / drho in Tokovinin's LSquares1
     for(i = 0; i < nmeas0; i++) {
       if(weight0[i] > 0.) sigma_rho0[i] = ABS(rho0[i]) / weight0[i];
        else sigma_rho0[i] = 10.;
       }
     break;

// Minimum format
   case 2:
     status_v = read_mini_visual_datafile(fname.mb_str(), epoch0, theta0, rho0, 
                                        weight0, notes_meas0, NMEAS_MAX, 
                                        &nmeas0);

// Fake values for sigma since they are missing in those files:
// weight = rho / drho in Tokovinin's LSquares1
     for(i = 0; i < nmeas0; i++) {
       if(weight0[i] > 0.) sigma_rho0[i] = ABS(rho0[i]) / weight0[i];
        else sigma_rho0[i] = 10.;
       }
     break;

// Tokovinin's format
   case 3:

// Read visual-measurement data
     status_v = read_toko_visual_datafile(fname.mb_str(), epoch0, theta0, rho0, 
                                          sigma_rho0, notes_meas0, 
                                          NMEAS_MAX, &nmeas0);

// Read radial velocity data
     status_rv = read_toko_rv_datafile(fname.mb_str(), epoch_rv10, rv10, 
                                    sigma_rv10, notes_rv10, 
                                    epoch_rv20, rv20, sigma_rv20, 
                                    notes_rv20, NMEAS_MAX, 
                                    &nrv10, &nrv20); 

     break;

// CORAVEL format
   case 4:

     status_rv = read_CORAVEL_rv_datafile(fname.mb_str(), epoch_rv10, rv10, 
                                    sigma_rv10, weight_rv10, notes_rv10, 
                                    epoch_rv20, rv20, sigma_rv20, weight_rv20, 
                                    notes_rv20, NMEAS_MAX, 
                                    &nrv10, &nrv20); 
     break;
 }

// If successful, transfer data to private arrays:
   if(status_v == 0)  {
     nmeas1 = nmeas0;
     for(i = 0; i < nmeas1; i++) {
       epoch1[i] = epoch0[i];
       theta1[i] = theta0[i];
       rho1[i] = rho0[i];
       sigma_rho1[i] = sigma_rho0[i];
       strcpy(notes_meas1[i], notes_meas0[i]);
       }
     }

// If successful, transfer data to private arrays:
   if(status_rv == 0)  {
     nrv1 = nrv10;
     nrv2 = nrv20;
     for(i = 0; i < nrv1; i++) {
       epoch_rv1[i] = epoch_rv10[i];
       rv1[i] = rv10[i];
       sigma_rv1[i] = sigma_rv10[i];
       strcpy(notes_rv1[i], notes_rv10[i]);
       }
     for(i = 0; i < nrv2; i++) {
       epoch_rv2[i] = epoch_rv20[i];
       rv2[i] = rv20[i];
       sigma_rv2[i] = sigma_rv20[i];
       strcpy(notes_rv2[i], notes_rv20[i]);
       }
     }

// Compute residuals and enable/disable menu items:
UpdateDataAfterFileIsLoaded();
UpdateMenuAfterFileIsLoaded();

if(nmeas1 != 0 || nrv1 != 0 || nrv2 != 0) {

// Display the content of the file on TextPanel:
  OrbitDataTextCtrl->LoadFile(fname);

//******************
// Write to logbook:
// Logbook:

  Data_fname = fname; 
// Removes the directory name (since the full path is generally too long...)
  wxFileName::SplitPath(fname, &path, &short_fname, &extension);
  short_Data_fname = short_fname; 

  str1 = _T("%% Input data file: ") + short_Data_fname + _T("\n");
  WriteToLogbook(str1, true);
//******************************************************************
  }

return;
}
/************************************************************************
* Input OC6 visual orbit file
*************************************************************************/
void GpFrame::OnLoadOrbitFile( wxCommandEvent &event )
{
int iformat = 1;
wxString str1, fname;

switch(event.GetId()) {
 default :
 case ID_LOAD_SCARDIA_ORBITFILE :
   str1 = wxT("Select Visual Orbit File (Scardia's format)");
   iformat = 1;
   break;
 case ID_LOAD_OC6_ORBITFILE :
   str1 = wxT("Select Visual Orbit File (OC6 format)");
   iformat = 2;
   break;
 case ID_LOAD_TOKO_ORBITFILE :
   str1 = wxT("Select Visual & RV Orbit File (Tokovinin's format)");
   iformat = 3;
   break;
}

// Display File Selector Dialog: 
 fname = wxFileSelector( str1);

 if(fname.empty()) return;

 LoadOrbitFile(iformat, fname);

// Select Orbit file in TextPanel: 
 ShowOrbitFileButton->SetValue(true);
 TextCancelChangesButton->Show(false);
 TextValidChangesButton->Show(false);
 TextSaveToFileButton->Show(false);

}
/************************************************************************
* Input ascii orbit file
*
* INPUT:
* iformat:
*   1 = Scardia
*   2 = OC6 
*   3 = Tokovinin 
*************************************************************************/
void GpFrame::LoadOrbitFile(const int iformat, wxString &fname)
{
wxString str1, path, extension, short_fname;
char in_orbit_datafile[256];
double oelmnts0[NELEMENTS], oelmnts_err0[NELEMENTS], orbit_equinox0;
int status, i, orbit_type0;
char object_name0[128];
double right_ascension0, declination0;

// Initialize elements and errors:
 for(i = 0; i < N10; i++) {
   oelmnts0[i] = 0.;
   oelmnts_err0[i] = 0.;
   }
 orbit_type0 = 0;

// Do some initialization:
 ResetBeforeLoadingOrbitFile();

// Read file and decode orbital parameters:
 switch(iformat) {
   default:
   case 1:
     status = rd_visual_oelements_scardia(fname.mb_str(), oelmnts0, 
                                          &orbit_equinox0, object_name0,
                                          &right_ascension0, &declination0,
                                          &orbit_type0);
     break;
   case 2:
     status = rd_visual_oelements_OC6(fname.mb_str(), oelmnts0, 
                                      &orbit_equinox0, object_name0,
                                      &right_ascension0, &declination0,
                                      &orbit_type0);
     break;
   case 3:
     status = rd_oelements_toko(fname.mb_str(), oelmnts0, 
                                &orbit_equinox0, object_name0, 
                                &right_ascension0, &declination0, 
                                &orbit_type0);
     break;
   }

// Make some initialization:
if(status == 0) {
  orbit_equinox2 = orbit_equinox0;
  strcpy(object_name1, object_name0);
  right_ascension1 = right_ascension0;
  declination1 =  declination0;
  UpdateElementsAfterFileIsLoaded(oelmnts0, oelmnts_err0, N10, orbit_type0);
  Orbit_fname = fname; 

// Display the content of the logfile on TextPanel:
  OrbitDataTextCtrl->LoadFile(fname);
  m_textpanel_type = -1;

// Compute residuals and enable/disable menu items:
  UpdateMenuAfterFileIsLoaded();

//******************
// Write to logbook:
// Logbook:

// Removes the directory name (since the full path is generally too long...)
  wxFileName::SplitPath(fname, &path, &short_fname, &extension);
  short_Orbit_fname = short_fname + _T(".") + extension;

  str1 = _T("%%\n%% Input orbit file: ") + short_Orbit_fname + _T("\n");
  WriteToLogbook(str1, true);
//******************************************************************
  }

return;
}
/*************************************************************************
* Save orbit and data to Tokovinin's formatted file
*
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
*
*
*************************************************************************/
void GpFrame::OnSaveOrbitAndData( wxCommandEvent &WXUNUSED(event) )
{
wxString savefilename;

savefilename = wxFileSelector( 
            wxT("Save orbit and data to file (Tokovinin's format)"),
            wxT(""), wxT(""), wxT("dat|DAT|inp"),
            wxT("Files (*.dat;*.DAT;*.inp)|*.dat;*.DAT;*.inp"),
            wxFD_SAVE, this);

   if ( savefilename.empty() ) return;

/* Prepare orbit and data text: 
* panel_type: 0: orbit and data, 1: orbit only, 2: data only 3: residuals 
*/
 UpdatePanelText(0);

// Save to file:
 OrbitDataTextCtrl->SaveFile(savefilename);

return;
}
/************************************************************************
* Save data to Tokovinin's formatted file
*************************************************************************/
void GpFrame::OnSaveData( wxCommandEvent &WXUNUSED(event) )
{
wxString savefilename;

savefilename = wxFileSelector( 
            wxT("Save data to file (Tokovinin's format)"),
            wxT(""), wxT(""), wxT("dat|DAT|inp"),
            wxT("Files (*.dat;*.DAT;*.inp)|*.dat;*.DAT;*.inp"),
            wxFD_SAVE, this);

   if ( savefilename.empty() ) return;

/* Prepare data text: 
* panel_type: 0: orbit and data, 1: orbit only, 2: data only 3: residuals 
*/
 UpdatePanelText(2);

// Save to file:
 OrbitDataTextCtrl->SaveFile(savefilename);

return;
}
/************************************************************************
* Save orbit to Tokovinin's formatted file
*************************************************************************/
void GpFrame::OnSaveOrbit( wxCommandEvent &WXUNUSED(event) )
{
wxString savefilename;

savefilename = wxFileSelector( 
            wxT("Save orbit to file (Tokovinin's format)"),
            wxT(""), wxT(""), wxT("dat|DAT|inp"),
            wxT("Files (*.dat;*.DAT;*.inp)|*.dat;*.DAT;*.inp"),
            wxFD_SAVE, this);

   if ( savefilename.empty() ) return;

/* Prepare orbit text: 
* panel_type: 0: orbit and data, 1: orbit only, 2: data only 3: residuals 
*/
 UpdatePanelText(1);

// Save to file:
 OrbitDataTextCtrl->SaveFile(savefilename);

return;
}
/************************************************************************
* Save residuals to file
*************************************************************************/
void GpFrame::OnSaveResiduals( wxCommandEvent &WXUNUSED(event) )
{
int status;
wxString savefilename;

   savefilename = wxFileSelector( 
            wxT("Save residuals to File (Tokovinin's format))"),
            wxT(""), wxT(""), wxT("dat|DAT|inp"),
            wxT("Files (*.dat;*.DAT;*.inp)|*.dat;*.DAT;*.inp"),
            wxFD_SAVE, this);

   if ( savefilename.empty() ) return;

/* Prepare orbit text: 
* panel_type: 0: orbit and data, 1: orbit only, 2: data only 3: residuals 
*/
 UpdatePanelText(3);

// Save to file:
 OrbitDataTextCtrl->SaveFile(savefilename);

return;
}
/************************************************************************
* Select notebook page 
*************************************************************************/
void GpFrame::OnSelectPage( wxBookCtrlEvent &event )
{
switch(event.GetId()) {
 default :
 case ID_TEXT_PAGE :
   iPage = 0;
   break;
 case ID_GRAPHIC_PAGE :
   iPage = 1;
   break;
 case ID_WIDGET_PAGE :
   iPage = 2;
   break;
}

return;
}
