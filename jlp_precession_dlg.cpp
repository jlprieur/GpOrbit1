/******************************************************************************
* jlp_precession_dlg.cpp
* Dialog box used when performing cosmetic corrections on images
*
* Author:      JLP 
* Version:     28/01/2009
******************************************************************************/
#include "jlp_precession_dlg.h"

//*************************************************************************
enum
{
   ID_PRECESS_OBJECT    = 800,
   ID_PRECESS_RA,
   ID_PRECESS_DEC,
   ID_PRECESS_EQUINOX,
   ID_PRECESS_OK,
   ID_PRECESS_CANCEL,
};

BEGIN_EVENT_TABLE(JLP_Precession_Dlg, wxDialog)
EVT_BUTTON  (ID_PRECESS_OK, JLP_Precession_Dlg::OnOKButton)
EVT_BUTTON  (ID_PRECESS_CANCEL, JLP_Precession_Dlg::OnCancelButton)
EVT_TEXT    (ID_PRECESS_OBJECT, JLP_Precession_Dlg::OnChangeParam)
EVT_TEXT    (ID_PRECESS_RA, JLP_Precession_Dlg::OnChangeParam)
EVT_TEXT    (ID_PRECESS_DEC, JLP_Precession_Dlg::OnChangeParam)
EVT_TEXT    (ID_PRECESS_EQUINOX, JLP_Precession_Dlg::OnChangeParam)
END_EVENT_TABLE()

/********************************************************************
* Constructor:
********************************************************************/
JLP_Precession_Dlg::JLP_Precession_Dlg(wxFrame *parent, double RA_0, 
                                       double DEC_0, double Equinox_0, 
                                       char *object_name0,
                                       const wxString &title)
        : wxDialog(parent, -1, title, wxPoint(400,100), wxDefaultSize,
                   wxDEFAULT_DIALOG_STYLE | wxRESIZE_BORDER)
{
wxString Object_str, RA_str, DEC_str, Equinox_str;
wxBoxSizer *topsizer;
wxFlexGridSizer *fgs1;
int nrows, ncols, vgap = 12, hgap = 12;

// Save input parameters to private variables:
  strcpy(object_name1, object_name0);
  RA_1 = RA_0;
  DEC_1 = DEC_0;
  Equinox_1 = Equinox_0;

// Initialize the text strings:
  Object_str = wxString(object_name1);
  RA_str.Printf(_T("%.4f"), RA_1);
  DEC_str.Printf(_T("%.4f"), DEC_1);
  Equinox_str.Printf(_T("%.2f"), Equinox_1);

// To avoid initialization problems with Windows:
// (An event is sent to "ChangeParameters"
//  as soon as a Text control is created...)
  initialized = 0; 

// Flexible grid sizer:
  nrows = 4;
  ncols = 2;
  fgs1 = new wxFlexGridSizer(nrows, ncols, vgap, hgap);

// Create the text controls: 
  TextCtrl_Object = new wxTextCtrl( this, ID_PRECESS_OBJECT, Object_str,
                                    wxPoint(-1, -1), wxSize(250, 30));
  TextCtrl_RA = new wxTextCtrl( this, ID_PRECESS_RA, RA_str,
                                    wxPoint(-1, -1), wxSize(150, 30));
  TextCtrl_DEC = new wxTextCtrl( this, ID_PRECESS_DEC, DEC_str,
                                    wxPoint(-1, -1), wxSize(150, 30));
  TextCtrl_Equinox = new wxTextCtrl( this, ID_PRECESS_EQUINOX, Equinox_str,
                                    wxPoint(-1, -1), wxSize(150, 30));
 
topsizer = new wxBoxSizer( wxVERTICAL );


// Sizer surrounded with a rectangle, with a title on top:
wxStaticBoxSizer *Precession_sizer = new wxStaticBoxSizer(wxVERTICAL, this, 
                                   _T(" Object/orbit parameters"));
    
  fgs1->Add( new wxStaticText( this, wxID_ANY, _T("Object name:") ));
  fgs1->Add(TextCtrl_Object); 
  fgs1->Add( new wxStaticText( this, wxID_ANY, _T("R.A.:") ));
  fgs1->Add(TextCtrl_RA, 0, wxRIGHT, 0); 
  fgs1->Add( new wxStaticText( this, wxID_ANY, _T("DEC:") ));
  fgs1->Add(TextCtrl_DEC, 0, wxRIGHT, 0); 
  fgs1->Add( new wxStaticText( this, wxID_ANY, _T("Equinox:") ));
  fgs1->Add(TextCtrl_Equinox, 0, wxRIGHT, 0); 

  Precession_sizer->Add(fgs1, 0, wxALIGN_CENTER|wxALL, 20);
  topsizer->Add(Precession_sizer, 0, wxALIGN_CENTER|wxALL, 20);

wxBoxSizer *button_sizer = new wxBoxSizer( wxHORIZONTAL );

//create two buttons that are horizontally unstretchable, 
  // with an all-around border with a width of 10 and implicit top alignment
 button_sizer->Add(
    new wxButton(this, ID_PRECESS_OK, _T("OK") ), 0, wxALIGN_LEFT|wxALL, 10);

 button_sizer->Add(
   new wxButton(this, ID_PRECESS_CANCEL, _T("Cancel") ), 0, wxALIGN_CENTER|wxALL, 10);

  //create a sizer with no border and centered horizontally
  topsizer->Add(button_sizer, 0, wxALIGN_CENTER);

  SetSizer(topsizer);      // use the sizer for layout

  topsizer->SetSizeHints( this );   // set size hints to honour minimum size

  initialized = 1234; 
return;
}
/**************************************************************************
* Handle "OK" button:
**************************************************************************/
void JLP_Precession_Dlg::OnOKButton( wxCommandEvent& WXUNUSED(event) )
{
// Close dialog and return status = 0:
  EndModal(0); 
}
/**************************************************************************
* Handle "Cancel" button:
**************************************************************************/
void JLP_Precession_Dlg::OnCancelButton( wxCommandEvent& WXUNUSED(event) )
{
// Close dialog and return status = 1:
  EndModal(1); 
}
/**************************************************************************
* Handle text editing 
*
    TextCtrl_Object
    TextCtrl_RA
    TextCtrl_DEC_
    TextCtrl_Equinox
*
**************************************************************************/
void JLP_Precession_Dlg::OnChangeParam( wxCommandEvent& event )
{
double old_value, new_value;
int status = 0;
wxString w_str;

// First check that all text controls are created:
 if(initialized != 1234) return;

  switch (event.GetId())
  {
   case ID_PRECESS_OBJECT:
    {
// Get new value
      w_str = TextCtrl_Object->GetValue();
      strcpy(object_name1, w_str.mb_str());
      break;
    }
   case ID_PRECESS_RA:
    {
// Get new value
    if(TextCtrl_RA->GetValue().ToDouble(&new_value)) {
      old_value = RA_1;
      RA_1 = new_value; 
      if(DataIsOK() == false) status = -1;
      }
// If bad new value restore previous value:
    if(status) { 
      RA_1 = old_value; 
      wxMessageBox(_T("Bad RA value!"), _T("Object/orbit parameters"), 
                   wxICON_ERROR);
      w_str.Printf(_T("%.4f"), RA_1); 
      TextCtrl_RA->SetValue(w_str);
      }
      break;
    }
   case ID_PRECESS_DEC:
    {
// Get new value
    if(TextCtrl_DEC->GetValue().ToDouble(&new_value)) {
      old_value = DEC_1;
      DEC_1 = new_value; 
      if(DataIsOK() == false) status = -1;
      }
// If bad new value restore previous value:
    if(status) { 
      DEC_1 = old_value;
      wxMessageBox(_T("Bad DEC value!"), _T("Object/orbit parameters"), 
                   wxICON_ERROR);
      w_str.Printf(_T("%.4f"), DEC_1); 
      TextCtrl_DEC->SetValue(w_str);
      }
      break;
    }
   case ID_PRECESS_EQUINOX:
    {
// Get new value
    if(TextCtrl_Equinox->GetValue().ToDouble(&new_value)) {
      old_value = Equinox_1;
      Equinox_1 = new_value; 
      if(DataIsOK() == false) status = -1;
      }
// If bad new value restore previous value:
    if(status) { 
      Equinox_1 = old_value;
      wxMessageBox(_T("Bad Equinox value!"), _T("Object/orbit parameters"), 
                   wxICON_ERROR);
      w_str.Printf(_T("%.2f"), Equinox_1); 
      TextCtrl_Equinox->SetValue(w_str);
      }
      break;
    }
  }  // EOF switch
return;
}
