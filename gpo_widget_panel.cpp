/****************************************************************************
* Name: gpo_widget_panel.cpp 
* 
* JLP
* Version 17/07/2015
****************************************************************************/
#include "gpo_frame.h"
#include "gpo_frame_id.h"

#include "gpo_defs.h"      // YEAR_TO_DAYS
#include "gpo_rw_files.h"  // besselian_to_julian() 

/********************************************************************
* Widget panel with orbital elements, results, etc. 
********************************************************************/
int GpFrame::WidgetPanel_Setup()
{
wxBoxSizer *w_topsizer, *w_hsizer2;
int status;

  w_topsizer = new wxBoxSizer(wxVERTICAL);

  w_hsizer2 = new wxBoxSizer( wxHORIZONTAL );

// Create oelmt_sizer for displaying orbital elements:
  WidgetPanel_OElmnts_Setup();
  w_hsizer2->Add(oelmt_sizer, 0, wxALIGN_CENTRE_VERTICAL | wxALL, 6);

// Create results_sizer for displaying the data fit results 
  WidgetPanel_Results_Setup();
  w_hsizer2->Add(results_sizer, 0, wxALIGN_CENTRE_VERTICAL | wxALL, 10);

  w_topsizer->Add(w_hsizer2, 0, wxALIGN_CENTER | wxTOP, 10);

  m_WidgetPanel->SetSizer(w_topsizer);

  Centre();

// Select the first data fit method:
  data_fit_method1 = 1;
  DataFitButton1->SetValue(true);

return(0);
}
/********************************************************************
* Select labels for displaying the orbital elements
********************************************************************/
int GpFrame::WidgetPanel_OElmnts_str_Setup(const int orbit_type0)
{
if(orbit_type0 == 3) {
// i=251 "period  (yr) :"
  oelmt_str[0] = Str0[iLang][251];
// i=252 "T_peri  (yr) :"
  oelmt_str[1] = Str0[iLang][252];
} else {
// i=253 "period  (d) :"
  oelmt_str[0] = Str0[iLang][253];
// i=254 "T_peri  (JD) :"
  oelmt_str[1] = Str0[iLang][254];
}

// i=255 "eccentricity :"
  oelmt_str[2] = Str0[iLang][255]; 
// i=256 "a  (arcsec)  :"
  oelmt_str[3] = Str0[iLang][256]; 
  oelmt_str[4] = _T("Omega  (deg) :");
  oelmt_str[5] = _T("omega  (deg) :");
  oelmt_str[6] = _T("incl.  (deg) :");

  oelmt_str[7] = _T("K1  (km/s)  :");
  oelmt_str[8] = _T("K2  (km/s)  :");
  oelmt_str[9] = _T("V0  (km/s)  :");
}
/********************************************************************
* Sub-panel with the orbital elements 
// Create oelmt_sizer for displaying orbital elements:
********************************************************************/
int GpFrame::WidgetPanel_OElmnts_Setup()
{
wxBoxSizer *w_hsizer0, *w_hsizer1, *orb_elements_sizer; 
int i, irows, icols, vgap, hgap = 12;

// Setup with orbit type = 3 (visual orbit) by default:
WidgetPanel_OElmnts_str_Setup(3);

// Setting the flexible grid sizer at maximum size:
  irows = 1 + N10;
  icols = 4;
  vgap = 12;
  fgs1 = new wxFlexGridSizer(irows, icols, vgap, hgap);

// i=257 "  Name"
  wxStaticText *label_txt1 = new wxStaticText(m_WidgetPanel, -1, 
                                              Str0[iLang][257]);
  fgs1->Add(label_txt1);

// i=258 "   Value"
  wxStaticText *label_txt2 = new wxStaticText(m_WidgetPanel, -1, 
                                              Str0[iLang][258], 
                                              wxPoint(-1, -1), wxSize(80, 28));
  fgs1->Add(label_txt2);

// i=259 "Fixed"
  wxStaticText *label_txt3 = new wxStaticText(m_WidgetPanel, -1, 
                                              Str0[iLang][259]);
  fgs1->Add(label_txt3);

// i=260 "Error"
  wxStaticText *label_txt4 = new wxStaticText(m_WidgetPanel, -1,
                                              Str0[iLang][260]);
  fgs1->Add(label_txt4);

// Loading all the sizer items:
/* Orbital elements:
* oelements {'P', 'T', 'e', 'a', 'W', 'w', 'i', 'K1', 'K2', 'V0'};
* oelmnt[0] : P (years)
* oelmnt[1] : T_peri (years)
* oelmnt[2] : e=eccentricity 
* oelmnt[3] : a=semi-major axis  (arcsec) 
* oelmnt[4] : W=Omega, node arg. (deg) 
* oelmnt[5] : w=omega, periastron argument (deg) 
* oelmnt[6] : i=inclination (deg) 
* oelmnt[7] : K1 (km/s)
* oelmnt[8] : K2 (km/s)
* oelmnt[9] : V0 (km/s)
*/
  for(i = 0; i < N10; i++) {
    oelmt_name_static[i] = new wxStaticText(m_WidgetPanel, -1, oelmt_str[i]);
    oelmt_txt_ctrl[i] = new wxTextCtrl(m_WidgetPanel, -1, wxT(""),
                                   wxPoint(-1, -1), wxSize(120, 28));
    oelmt_fflag_checkbox[i] = new wxCheckBox(m_WidgetPanel, wxID_ANY, 
                                       wxT(""));
    oelmt_err_static[i] = new wxStaticText(m_WidgetPanel, -1, wxT("0."),
                                   wxPoint(-1, -1), wxSize(100, 28));
      fgs1->Add(oelmt_name_static[i]);
      fgs1->Add(oelmt_txt_ctrl[i]);
      fgs1->Add(oelmt_fflag_checkbox[i]);
      fgs1->Add(oelmt_err_static[i]);
    }

////////////////////////////////////////////////////////////////////////
// Sizer surrounded with a rectangle, with a title on top:
// i=261 "Orbital elements"
  oelmt_sizer = new wxStaticBoxSizer(wxVERTICAL, m_WidgetPanel,
                                     Str0[iLang][261]);

// horizontal box to put some space (10 px) around the table:
  orb_elements_sizer = new wxBoxSizer(wxHORIZONTAL);
  orb_elements_sizer->Add(fgs1, 0, wxALL, 10);
  oelmt_sizer->Add(orb_elements_sizer, 0, wxALL, 10);

///////////////////////////////////////////////////////////////////////
// Create two buttons: Valid and Cancel Changes:
 w_hsizer1 = new wxBoxSizer( wxHORIZONTAL );

// i=262 "Valid changes" 
 ValidOrbitChangesButton = new wxButton(m_WidgetPanel, ID_VALID_ORBIT_CHANGES,
                                        Str0[iLang][262]);
// i=263 "Cancel changes" 
 CancelOrbitChangesButton = new wxButton(m_WidgetPanel, ID_CANCEL_ORBIT_CHANGES,
                                         Str0[iLang][263]);

// Add buttons, horizontally unstretchable, with minimal size:
 w_hsizer1->Add( ValidOrbitChangesButton, 0);
 w_hsizer1->Add( CancelOrbitChangesButton, 0,
                    wxALIGN_RIGHT |wxLEFT | wxRIGHT, 20);

// Add button sizer with minimal size:
  oelmt_sizer->Add(w_hsizer1, 0, wxALIGN_CENTER | wxALL, 10);

// Border at bottom:
// oelmt_sizer->Add(-1, 25);

return(0);
}
/********************************************************************
* Sub-panel with the results of the data fit 
* Create datafit_sizer for displaying results of data fit
********************************************************************/
int GpFrame::WidgetPanel_Results_Setup()
{
wxBoxSizer *w_hsizer1, *w_hsizer2, *results1_sizer; 
int i, irows, icols, vgap, hgap = 12;
wxString result_name_str[NRESULTS_MAX];
wxString result_units_str[NRESULTS_MAX];

nresults1 = 13;

// Setting the flexible grid sizer at maximum size:
  irows = nresults1;
  icols = 3;
  vgap = 8;
  fgs2 = new wxFlexGridSizer(irows, icols, vgap, hgap);

// nmeas1, nrv1, nrv2, mean_sigma_rho_resid1, mean_sigma_theta_resid1,
// mean_sigma_rv1_resid1, mean_sigma_rv2_resid1
// i=264 "nber of visual meas. (meas1):"
  result_name_str[0] = Str0[iLang][264];
  result_units_str[0] = _T("");
// i=265 "nber of primary rv meas. (nrv1) :"
  result_name_str[1] = Str0[iLang][265]; 
  result_units_str[1] = _T("");
// i=266 "nber of secondary rv meas. (nrv2) :"
  result_name_str[2] = Str0[iLang][266]; 
  result_units_str[2] = _T("");
// i=267 "theta corrected for precession :"
  result_name_str[3] = Str0[iLang][267];
  result_units_str[3] = _T("");
// i=268 "mean rho residual :"
  result_name_str[4] = Str0[iLang][268];
// i=270 "(arcsec)"
  result_units_str[4] = Str0[iLang][270];
// i=269 "mean sigma rho residual :"
  result_name_str[5] = Str0[iLang][269];
// i=270 "(arcsec)"
  result_units_str[5] = Str0[iLang][270];
// i=271 "mean theta residual :"
  result_name_str[6] = Str0[iLang][271];
  result_units_str[6] = _T("(deg.)");
// i=272 "mean sigma theta residual :"
  result_name_str[7] = Str0[iLang][272];
  result_units_str[7] = _T("(deg.)");
// i=273 "mean sigma rv1 residual :"
  result_name_str[8] = Str0[iLang][273]; 
  result_units_str[8] = _T("(km/s)");
// i=274 "mean sigma rv2 residual :"
  result_name_str[9] = Str0[iLang][274]; 
  result_units_str[9] = _T("(km/s)");
// i=277 "Parallax :"
  result_name_str[10] = Str0[iLang][277]; 
  result_units_str[10] = _T("(mas)");
// i=278 "Semi major axis :"
  result_name_str[11] = Str0[iLang][278]; 
  result_units_str[11] = _T("(au)");
// i=279 "Total mass :"
  result_name_str[12] = Str0[iLang][279]; 
  result_units_str[12] = _T("(Msol)");

// Loading all the results items:
  for(i = 0; i < nresults1; i++) {
    result_name_static[i] = new wxStaticText(m_WidgetPanel, -1, 
                                             result_name_str[i]);
    result_value_static[i] = new wxStaticText(m_WidgetPanel, -1, wxT("0."),
                                   wxPoint(-1, -1), wxSize(120, 28));
    result_units_static[i] = new wxStaticText(m_WidgetPanel, -1, 
                                             result_units_str[i]);
    fgs2->Add(result_name_static[i]);
    fgs2->Add(result_value_static[i]);
    fgs2->Add(result_units_static[i]);
    }

////////////////////////////////////////////////////////////////////////
// Sizer surrounded with a rectangle, with a title on top:
// i=280 "Data fit"
  results_sizer = new wxStaticBoxSizer(wxVERTICAL, m_WidgetPanel,
                                       Str0[iLang][280]);

// horizontal box to put some space (10 px) around the table:
  results1_sizer = new wxBoxSizer(wxHORIZONTAL);
  results1_sizer->Add(fgs2, 0, wxALL, 10);
  results_sizer->Add(results1_sizer, 0, wxALL, 10);

///////////////////////////////////////////////////////////////////////
// Create new button:  "New fit"
// Button for starting new fit computation:
  w_hsizer1 = new wxBoxSizer( wxHORIZONTAL );
// i=281 "New orbit computation"
  NewFitButton = new wxButton(m_WidgetPanel, ID_NEW_ORBIT_FIT,
                              Str0[iLang][281]);
// Add button, horizontally unstretchable, with minimal size:
  w_hsizer1->Add( NewFitButton, 0);
  results_sizer->Add(w_hsizer1, 0, wxALIGN_CENTER | wxALL, 10);

///////////////////////////////////////////////////////////////////////
// Create radio buttons 
 DataFitButton1 = new wxRadioButton(m_WidgetPanel, ID_DATA_FIT1,
                                    _T("Gauss-Newton (Tokovinin 1992)"), wxPoint(-1,-1),
                                    wxSize(-1,-1), wxRB_GROUP);
 DataFitButton2 = new wxRadioButton(m_WidgetPanel, ID_DATA_FIT2,
                                    _T("Levenberg-Marquardt"));

// Add buttons with minimal size:
 w_hsizer2 = new wxBoxSizer( wxHORIZONTAL );
 w_hsizer2->Add( DataFitButton1, 0);
 w_hsizer2->Add( DataFitButton2, 0, wxALIGN_RIGHT |wxLEFT | wxRIGHT, 20);

// Add button sizer with minimal size:
// with an all-around border with a width of 10 and implicit top alignment
 results_sizer->Add(w_hsizer2, 0, wxALIGN_CENTER | wxALL, 10);

// Border at bottom:
// results_sizer->Add(-1, 25);

return(0);
}
/********************************************************************
* Update data fit results values 
* Refresh screen by displaying internal parameters values
********************************************************************/
int GpFrame::WidgetPanel_Update_Results()
{
int i;
wxString sstr;

 if(initialized != 1234) return(1); 

// Static text:
// nmeas1, nrv1, nrv2, mean_sigma_rho_resid1, mean_sigma_theta_resid1,
// mean_sigma_rv1_resid1, mean_sigma_rv2_resid1

// result_value_static[0] : Str0[iLang][264], nber of visual meas. 
  sstr.Printf(_T("%d"), nmeas1);
  result_value_static[0]->SetLabel(sstr);

// result_value_static[1] : Str0[iLang][265], nber of primary rv meas. 
  sstr.Printf(_T("%d"), nrv1);
  result_value_static[1]->SetLabel(sstr);

// result_value_static[1] : Str0[iLang][266], nber of secondary rv meas. 
  sstr.Printf(_T("%d"), nrv2);
  result_value_static[2]->SetLabel(sstr);

// Precession correction:
// i=275 "yes"
  if(precession_corrected == 1) sstr = Str0[iLang][275];
// i=276 "no"
    else sstr = Str0[iLang][276];
// result_value_static[3] : yes/no (theta corrected for precession)
  result_value_static[3]->SetLabel(sstr);

// result_value_static[4] : Str0[iLang][268], mean rho residual
// Only 3 non-zero figures are enough for display, hence "%.3g" :
  sstr.Printf(_T("%.3g"), mean_rho_resid1);
  result_value_static[4]->SetLabel(sstr);

// result_value_static[5] : Str0[iLang][269], mean sigma rho residual
  sstr.Printf(_T("%.3g"), mean_sigma_rho_resid1);
  result_value_static[5]->SetLabel(sstr);

// result_value_static[6] : Str0[iLang][271], mean theta residual
  sstr.Printf(_T("%.3g"), mean_theta_resid1);
  result_value_static[6]->SetLabel(sstr);

// result_value_static[7] : Str0[iLang][272], mean sigma theta residual
  sstr.Printf(_T("%.3g"), mean_sigma_theta_resid1);
  result_value_static[7]->SetLabel(sstr);

// result_value_static[8] : Str0[iLang][273], mean rv1 residual
  sstr.Printf(_T("%.3g"), mean_sigma_rv1_resid1);
  result_value_static[8]->SetLabel(sstr);

// result_value_static[9] : Str0[iLang][274], mean rv2 residual
  sstr.Printf(_T("%.3g"), mean_sigma_rv2_resid1);
  result_value_static[9]->SetLabel(sstr);

// result_value_static[10] : Str0[iLang][277], Parallax 
  sstr.Printf(_T("%.3g +/- %.1f"), parallax1, parallax_err1);
  result_value_static[10]->SetLabel(sstr);
// result_value_static[11] : Str0[iLang][278], Semi major axis 
  sstr.Printf(_T("%.2f +/- %.2f"), semi_majax1, semi_majax_err1);
  result_value_static[11]->SetLabel(sstr);
// result_value_static[12] : Str0[iLang][279], Total mass 
  sstr.Printf(_T("%.2f +/- %.2f"), total_mass1, total_mass_err1);
  result_value_static[12]->SetLabel(sstr);

return(0);
}
/********************************************************************
* Compute the semi major axis in AU and the sum of masses M1+M2
*
* Orbital elements:
* oelements {'P', 'T', 'e', 'a', 'W', 'w', 'i', 'K1', 'K2', 'V0'};
* oelmnt[0] : P (years)
* oelmnt[1] : T_peri (years)
* oelmnt[2] : e=eccentricity 
* oelmnt[3] : a=semi-major axis  (arcsec) 
* oelmnt[4] : W=Omega, node arg. (deg) 
* oelmnt[5] : w=omega, periastron argument (deg) 
* oelmnt[6] : i=inclination (deg) 
* oelmnt[7] : K1 (km/s)  
* oelmnt[8] : K2 (km/s)
* oelmnt[9] : V0 (km/s)
*
********************************************************************/
int GpFrame::ComputeMassAndSemiAxis(double *oelmnts0, double *oelmnts_err0)
{
int orbit_type0, status;

 if(initialized != 1234 || m_jlp_orbit1 == NULL) {
    return(1);
   }

semi_majax1 = 0.;
semi_majax_err1 = 0.;
total_mass1 = 0.;
total_mass_err1 = 0.;

// To avoid errors with null values: 
if((parallax1 <= 0.) || (oelmnts0[0] <= 0.) || (oelmnts0[3] <= 0.)) return(-1); 


// Semi major axis: conversion to astronomical units (AU) 
// a_AU = a / parallax_in_arcsec:
// a_AU = a_smaxis[0] / (parallax / 1000.);
// err_a_AU = a_AU * sqrt(SQUARE(err_a_smaxis[0] / a_smaxis[0])
//                       + SQUARE(err_parallax / parallax));
semi_majax1 = oelmnts0[3] / (parallax1 / 1000.);
semi_majax_err1 = semi_majax1 * sqrt(SQUARE(oelmnts_err0[3] / oelmnts0[3])
                       + SQUARE(parallax_err1 / parallax1));
// Sum of the masses: computation of value and error 
// total_masses = pow(a_AU, 3) / (Period[0] * Period[0]);
// err_total_masses = total_masses * sqrt( 9.* SQUARE(err_a_smaxis[0] / a_smaxis[0])
//                               + 9.* SQUARE(err_parallax / parallax)
//                                + 4. * SQUARE(err_Period[0] / Period[0]));
 total_mass1 = pow(semi_majax1, 3) / (oelmnts0[0] * oelmnts0[0]);
 total_mass_err1 = total_mass1 * sqrt( 9.* SQUARE(semi_majax_err1 / semi_majax1)
                               + 9.* SQUARE(parallax_err1 / parallax1)
                               + 4. * SQUARE(oelmnts_err0[0] / oelmnts0[0]));
return(0);
}
/********************************************************************
* Update orbital elements values 
* Refresh screen by displaying internal parameters values
* ckeyword[NELEMENTS] = {'P', 'T', 'e', 'a', 'W', 'w', 'i', 'K1', 'K2', 'V0'};
********************************************************************/
int GpFrame::WidgetPanel_Update_OElements()
{
double oelmnts0[NELEMENTS], oelmnts_err0[NELEMENTS];
int oelmnts_fflags0[NELEMENTS];
int orbit_type0, status;

 if(initialized != 1234 || m_jlp_orbit1 == NULL) {
    return(1); 
   } 

 m_jlp_orbit1->GetOrbitalElements(oelmnts0, oelmnts_err0, N10, &orbit_type0);
 m_jlp_orbit1->GetFixedFlags(oelmnts_fflags0, N10);

 status = WidgetPanel_Display_OElements(oelmnts0, oelmnts_err0, 
                                        oelmnts_fflags0, N10);

return(status);
}
/********************************************************************
* Update orbital elements values 
* Refresh screen by displaying internal parameters values
* ckeyword[NELEMENTS] = {'P', 'T', 'e', 'a', 'W', 'w', 'i', 'K1', 'K2', 'V0'};
********************************************************************/
int GpFrame::WidgetPanel_Display_OElements(double *oelmnts0, 
                                           double *oelmnts_err0, 
                                           int *oelmnts_fflags0, 
                                           const int noelmnts_max) 
{
double elmts[NELEMENTS], elmts_err[NELEMENTS];
int i, orbit_type_from_elements0;
wxString sstr, sstr1, sstr2, sstr3, sstr4;

 if(initialized != 1234 || m_jlp_orbit1 == NULL) {
    return(1); 
   } 

// i=282 "Orbital elements"
// i=283 "Orbital elements of"
  if(*object_name1) { 
    sstr1 = Str0[iLang][283];
    sstr2.Printf(wxT(" %s"), object_name1);
  } else {
    sstr1 = Str0[iLang][282];
    sstr2 = _T("");
  }

// i=284 "(equinox"
  if(orbit_equinox2 > 100.) {
    sstr3 = wxT(" ") + Str0[iLang][284];
    sstr4.Printf(wxT(" %.1f)"), orbit_equinox2);
  } else {
    sstr3 = _T("");
    sstr4 = _T("");
  } 

// Title:
  sstr = sstr1 + sstr2 + sstr3 + sstr4;

// Update title:
  oelmt_sizer->GetStaticBox()->SetLabel(sstr);

// Names: 
// Setup with current orbit type:
  orbit_type_from_elements0 = m_jlp_orbit1->GetOrbitTypeFromElements();
  WidgetPanel_OElmnts_str_Setup(orbit_type_from_elements0);
  for(i = 0; i < N10; i++) oelmt_name_static[i]->SetLabel(oelmt_str[i]);

  ConvertElementToSpectro(oelmnts0, oelmnts_err0, elmts, elmts_err,
                          N10, orbit_type_from_elements0);

// Text ctrl and static text:
for(i = 0; i < N10; i++) {
// 7 figures are needed for correct display, hence "%.7g" :
  sstr.Printf(_T("%.7g"), (float)elmts[i]);
  oelmt_txt_ctrl[i]->SetValue(sstr);
// Only 3 non-zero figures are enough for display, hence "%.3g" :
  sstr.Printf(_T("%.3g"), (float)elmts_err[i]);
  oelmt_err_static[i]->SetLabel(sstr);
  }

// Check boxes: 
 for(i = 0; i < N10; i++) {
   if(oelmnts_fflags0[i]) { 
     oelmt_fflag_checkbox[i]->SetValue(true);
    } else {
     oelmt_fflag_checkbox[i]->SetValue(false);
    }
   }

return(0);
}
/**************************************************************************
* Handle "CancelOrbitChanges" button:
**************************************************************************/
void GpFrame::OnCancelOrbitChanges( wxCommandEvent& WXUNUSED(event) )
{
double oelmnts0[NELEMENTS], oelmnts_err0[NELEMENTS];
int orbit_type0;

// Refresh screen by displaying internal parameters values
 WidgetPanel_Update_OElements();

// Compute new residuals to update Results Panel
//
// Compute new residuals in Residuals_txt and load then in OrbitDataText
// panel_type = 3 for residuals
 UpdatePanelText(3);

// Compute new value of mass and semi-axis in au with internal values:
 m_jlp_orbit1->GetOrbitalElements(oelmnts0, oelmnts_err0, N10, &orbit_type0);
 ComputeMassAndSemiAxis(oelmnts0, oelmnts_err0);

 WidgetPanel_Update_Results();

}
/**************************************************************************
* Handle "ValidOrbitChanges" button:
**************************************************************************/
void GpFrame::OnValidOrbitChanges( wxCommandEvent& WXUNUSED(event) )
{
double oelmnts[NELEMENTS], oelmnts_err[NELEMENTS];
double oelmnts0[NELEMENTS], oelmnts_err0[NELEMENTS];
int oelmnts_fflags[NELEMENTS];
double ww;
int i, orbit_type0;
wxString sstr, results_str;

 if((initialized != 1234) || (m_jlp_orbit1 == NULL)) return; 

// Get info from internal orbit object: 
 m_jlp_orbit1->GetOrbitalElements(oelmnts, oelmnts_err, N10, &orbit_type0);
 m_jlp_orbit1->GetFixedFlags(oelmnts_fflags, N10);

// Read displayed values:
 for(i = 0; i < N10; i++) {
     if(oelmt_txt_ctrl[i]->GetValue().ToDouble(&ww)) {
       oelmnts[i] = ww;
     }
     if(oelmt_fflag_checkbox[i]->IsChecked()) {
       oelmnts_fflags[i] = 1;
       } else {
       oelmnts_fflags[i] = 0;
       }
   }

// Conversion if needed
 ConvertElementFromSpectro(oelmnts, oelmnts_err,
                           oelmnts0, oelmnts_err0, N10, orbit_type0);

// Save values to internal parameters
 m_jlp_orbit1->SetOrbitalElements(oelmnts0, oelmnts_err0, N10, orbit_type0);
 m_jlp_orbit1->SetFixedFlags(oelmnts_fflags, N10);

// Compute new residuals in Residuals_txt and load then in OrbitDataText
// panel_type = 3 for residuals
 UpdatePanelText(3);

// Compute new value of mass and semi-axis in au:
 ComputeMassAndSemiAxis(oelmnts, oelmnts_err);
 WidgetPanel_Update_Results();

// JLP2018:  Write results to logbook:
 results_str.Printf("OK: nmeas=%d mean residuals: %.2f+/-%.2f (theta) %.3f+/-%.3f (rho)\n",
                      nmeas1, mean_theta_resid1,
                      mean_sigma_theta_resid1, mean_rho_resid1,
                      mean_sigma_rho_resid1);
 WriteToLogbook(results_str, false);


}
/**************************************************************************
* Handle "NewOrbitComputation" button:
**************************************************************************/
void GpFrame::OnNewOrbitFit( wxCommandEvent& WXUNUSED(event) )
{
int status, i;
double oelmnts[NELEMENTS], oelmnts_err[NELEMENTS];
double mean_theta_resid, mean_rho_resid; 
double mean_sigma_theta_resid, mean_sigma_rho_resid; 
double mean_sigma_rv1_resid, mean_sigma_rv2_resid, chisq2;
int oelmnts_fflags[NELEMENTS];
wxString fit_results_str;

if(initialized != 1234 || m_jlp_orbit1 == NULL) return;
 
 if( data_fit_method1 == 1) {
  status = m_jlp_orbit1->LSquaresFit1(oelmnts, oelmnts_err, N10, 
                                      &mean_theta_resid,
                                      &mean_sigma_theta_resid, 
                                      &mean_rho_resid,
                                      &mean_sigma_rho_resid, 
                                      &mean_sigma_rv1_resid,
                                      &mean_sigma_rv2_resid,
                                      &chisq2,
                                      fit_results_str);
 } else {
  status = m_jlp_orbit1->LSquaresFit2(oelmnts, oelmnts_err, N10, 
                                      &mean_theta_resid,
                                      &mean_sigma_theta_resid, 
                                      &mean_rho_resid,
                                      &mean_sigma_rho_resid, 
                                      &mean_sigma_rv1_resid,
                                      &mean_sigma_rv2_resid,
                                      &chisq2,
                                      fit_results_str);
 }

   if(status == 0) {
    m_jlp_orbit1->GetFixedFlags(oelmnts_fflags, N10);

// Display new orbital elements and their errors:
    WidgetPanel_Display_OElements(oelmnts, oelmnts_err,
                                  oelmnts_fflags, N10);

// Overwrite sigma residuals, since they are slightly different 
// with LSquaresFit (that uses the measurement errors instead of the weights):
    mean_theta_resid1 = mean_theta_resid;
    mean_rho_resid1 = mean_rho_resid;
    mean_sigma_theta_resid1 = mean_sigma_theta_resid;
    mean_sigma_rho_resid1 = mean_sigma_rho_resid;
    mean_sigma_rv1_resid1 = mean_sigma_rv1_resid;
    mean_sigma_rv2_resid1 = mean_sigma_rv2_resid;

// Compute new value of mass and semi-axis in au:
    ComputeMassAndSemiAxis(oelmnts, oelmnts_err);

    WidgetPanel_Update_Results();

// Clear graphic panel:
    m_GraphicPanel->wxGP_ClearDrawingDisplay();
    menuPlot->Check(ID_PLOT_IDLE, true);
    }

// Write comments of results to logbook: 
  WriteToLogbook(fit_results_str, false);
 
}
/************************************************************************
* Hide/Show unused/used orbital elements in orbit panel:
*
* {'P', 'T', 'e', 'a', 'W', 'w', 'i', 'K1', 'K2', 'V0'};
* orbit_type = 1 visual & spectro SB1 (9 elements idx1= 1 1 1 1 1 1 1 1 0 1)
* orbit_type = 2 visual & spectro SB2 (10 elements idx1= 1 1 1 1 1 1 1 1 1 1)
* orbit_type = 3 visual (7 elements: idx2= 1 1 1 1 1 1 1 0 0 0)
* orbit_type = 4 spectroscopic SB1 (6 elements: idx3= 1 1 1 0 0 1 0 1 0 1)
* orbit_type = 5 spectroscopic SB2 (7 elements: idx3= 1 1 1 0 0 1 0 1 1 1)
*
*************************************************************************/
void GpFrame::UpdateOrbitPanelAfterFileIsLoaded()
{
bool show_flag;
int i, orbit_type0;

if(initialized != 1234 || m_jlp_orbit1 == NULL) return;

orbit_type0 = m_jlp_orbit1->GetOrbitTypeFromElements();

// Show all items
   show_flag = true;
   for(i = 0; i < 10; i++) {
     oelmt_name_static[i]->Show(show_flag);
     oelmt_txt_ctrl[i]->Show(show_flag);
     oelmt_fflag_checkbox[i]->Show(show_flag);
     oelmt_err_static[i]->Show(show_flag);
     }

// Hide some items if needed:
   show_flag = false;

// Visual orbit: hide 'K1', 'K2' and 'V0' items
   if(orbit_type0 == 3) {
     for(i = 7; i < 10; i++) {
       oelmt_name_static[i]->Show(show_flag);
       oelmt_txt_ctrl[i]->Show(show_flag);
       oelmt_fflag_checkbox[i]->Show(show_flag);
       oelmt_err_static[i]->Show(show_flag);
       }
     } else if ((orbit_type0 == 4)
              || (orbit_type0 == 5)) {
//
// Spectroscopic orbit: hide 'a', 'W' and i'' items
       i =  3;
       oelmt_name_static[i]->Show(show_flag);
       oelmt_txt_ctrl[i]->Show(show_flag);
       oelmt_fflag_checkbox[i]->Show(show_flag);
       oelmt_err_static[i]->Show(show_flag);
       i =  4;
       oelmt_name_static[i]->Show(show_flag);
       oelmt_txt_ctrl[i]->Show(show_flag);
       oelmt_fflag_checkbox[i]->Show(show_flag);
       oelmt_err_static[i]->Show(show_flag);
       i =  6;
       oelmt_name_static[i]->Show(show_flag);
       oelmt_txt_ctrl[i]->Show(show_flag);
       oelmt_fflag_checkbox[i]->Show(show_flag);
       oelmt_err_static[i]->Show(show_flag);
       if(orbit_type0 == 4) {
         i =  8;
         oelmt_name_static[i]->Show(show_flag);
         oelmt_txt_ctrl[i]->Show(show_flag);
         oelmt_fflag_checkbox[i]->Show(show_flag);
         oelmt_err_static[i]->Show(show_flag);
         }
     } else if(orbit_type0 == 1) {
         i =  8;
         oelmt_name_static[i]->Show(show_flag);
         oelmt_txt_ctrl[i]->Show(show_flag);
         oelmt_fflag_checkbox[i]->Show(show_flag);
         oelmt_err_static[i]->Show(show_flag);
     } 

// Update the screen with those settings:
   m_WidgetPanel->Layout();
}
/************************************************************************
* Hide/Show unused/used items in data fit panel:
*
* {'P', 'T', 'e', 'a', 'W', 'w', 'i', 'K1', 'K2', 'V0'};
* orbit_type = 1 visual & spectro SB1 (9 elements idx1= 1 1 1 1 1 1 1 1 0 1)
* orbit_type = 2 visual & spectro SB2 (10 elements idx1= 1 1 1 1 1 1 1 1 1 1)
* orbit_type = 3 visual (7 elements: idx2= 1 1 1 1 1 1 1 0 0 0)
* orbit_type = 4 spectroscopic SB1 (6 elements: idx3= 1 1 1 0 0 1 0 1 0 1)
* orbit_type = 5 spectroscopic SB2 (7 elements: idx3= 1 1 1 0 0 1 0 1 1 1)
*
*************************************************************************/
void GpFrame::UpdateDataFitPanelAfterFileIsLoaded()
{
double oelmnts0[NELEMENTS], oelmnts_err0[NELEMENTS];
bool show_flag;
int i, orbit_type0;

if(initialized != 1234 || m_jlp_orbit1 == NULL) return;

orbit_type0 = m_jlp_orbit1->GetOrbitType();

// Show all items
  show_flag = true;
  for(i = 0; i < nresults1; i++) {
    result_name_static[i]->Show(show_flag);
    result_value_static[i]->Show(show_flag);
    result_units_static[i]->Show(show_flag);
    }

// Hide some items if needed:
   show_flag = false;

// Visual orbit: 
   if(orbit_type0 == 3) {
// Hide nrv1, nrv2 items:
     for(i = 1; i <= 2; i++) {
       result_name_static[i]->Show(show_flag);
       result_value_static[i]->Show(show_flag);
       result_units_static[i]->Show(show_flag);
       }
// Hide mean_sigma_nrv1, mean_sigma_nrv2 items:
     for(i = 8; i <= 9; i++) {
       result_name_static[i]->Show(show_flag);
       result_value_static[i]->Show(show_flag);
       result_units_static[i]->Show(show_flag);
       }
     } else if ((orbit_type0 == 4)
              || (orbit_type0 == 5)) {
//
// Spectroscopic orbit: 
// Hide nmeas1 item
       i =  0;
       result_name_static[i]->Show(show_flag);
       result_value_static[i]->Show(show_flag);
       result_units_static[i]->Show(show_flag);
// Hide precession_corrected, mean_sigma_rho, mean_sigma_theta items:
     for(i = 3; i <= 7; i++) {
       result_name_static[i]->Show(show_flag);
       result_value_static[i]->Show(show_flag);
       result_units_static[i]->Show(show_flag);
       }
     } 

// Update the screen with those settings:
   m_WidgetPanel->Layout();

// JLP2018: Update results:
// Compute new value of mass and semi-axis in au with internal values:
 m_jlp_orbit1->GetOrbitalElements(oelmnts0, oelmnts_err0, N10, &orbit_type0);
 ComputeMassAndSemiAxis(oelmnts0, oelmnts_err0);
 WidgetPanel_Update_Results();

return;
}
/***********************************************************************
* Handle radio button to select data fit method 
************************************************************************/
void GpFrame::OnSelectDataFitMethod(wxCommandEvent &event)
{

switch(event.GetId()) {
 default :
 case ID_DATA_FIT1:
   data_fit_method1 = 1;
   break;
 case ID_DATA_FIT2:
   data_fit_method1 = 2;
   break;
}

return;
}
