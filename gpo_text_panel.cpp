/****************************************************************************
* Name: gpo_text_panel.cpp 
* 
* JLP
* Version 08/06/2015
****************************************************************************/
#include "gpo_frame.h"
#include "gpo_frame_id.h"
#include "jlp_precession_dlg.h"  // JLP_Precession_Dlg class...
#include "gpo_rw_files.h"        // precession_correction()

/********************************************************************
* Dislay/edit window
*
********************************************************************/
int GpFrame::TextPanel_Setup()
{
wxBoxSizer *topsizer, *hsizer1, *hsizer2, *hsizer3;

topsizer = new wxBoxSizer( wxVERTICAL );

// Create the text control (read-only here !) 
 OrbitDataTextCtrl = new wxTextCtrl(m_TextPanel, wxID_ANY, _T(""), 
                         wxPoint(-1,-1), wxSize(-1, -1),
                         wxTE_MULTILINE );
//                         wxTE_MULTILINE | wxTE_READONLY );

// 1 = Fill all space horizontally 
 hsizer1 = new wxBoxSizer( wxHORIZONTAL );
 hsizer1->Add(OrbitDataTextCtrl, 1, wxEXPAND);

// 1 = Fill all space vertically 
// border = 10 pix : left, right and top
 topsizer->Add(hsizer1, 1, wxEXPAND | wxLEFT | wxRIGHT | wxTOP, 10);

// Create radio buttons 
//i=220 "Original data file"
 ShowDataFileButton = new wxRadioButton(m_TextPanel, ID_SHOW_DATAFILE,
                                        Str0[iLang][220],
                                        wxPoint(-1,-1), wxSize(-1,-1), 
                                        wxRB_GROUP);
//i=221 "Original orbit file"
 ShowOrbitFileButton = new wxRadioButton(m_TextPanel, ID_SHOW_ORBITFILE,
                                        Str0[iLang][221]);
//i=222 "Data"
 ShowDataButton = new wxRadioButton(m_TextPanel, ID_SHOW_DATA,
                                        Str0[iLang][222]);
//i=223 "Orbit"
 ShowOrbitButton = new wxRadioButton(m_TextPanel, ID_SHOW_ORBIT,
                                        Str0[iLang][223]);
//i=224 "Orbit and data"
 ShowOrbitDataButton = new wxRadioButton(m_TextPanel, ID_SHOW_ORBITDATA,
                                        Str0[iLang][224]);
//i=225 "Residuals"
 ShowResidualsButton = new wxRadioButton(m_TextPanel, ID_SHOW_RESIDUALS,
                                        Str0[iLang][225]);

// Add buttons with minimal size: 
 hsizer2 = new wxBoxSizer( wxHORIZONTAL );
 hsizer2->Add( ShowDataFileButton, 0);
 hsizer2->Add( ShowOrbitFileButton, 0, wxALIGN_RIGHT |wxLEFT | wxRIGHT, 20);
 hsizer2->Add( ShowDataButton, 0, wxALIGN_RIGHT |wxLEFT | wxRIGHT, 20);
 hsizer2->Add( ShowOrbitButton, 0, wxALIGN_RIGHT |wxLEFT | wxRIGHT, 20);
 hsizer2->Add( ShowOrbitDataButton, 0, wxALIGN_RIGHT |wxLEFT | wxRIGHT, 20);
 hsizer2->Add( ShowResidualsButton, 0, wxALIGN_RIGHT |wxLEFT | wxRIGHT, 20);

// Add button sizer with minimal size: 
// with an all-around border with a width of 20 and implicit top alignment
 topsizer->Add(hsizer2, 0, wxALIGN_RIGHT | wxALL, 20);

///////////////////////////////////////////////////////////////////////
// Create two buttons: Valid and Cancel Changes:
 hsizer3 = new wxBoxSizer( wxHORIZONTAL );

// i=226 "Cancel text changes"
 TextCancelChangesButton = new wxButton(m_TextPanel, ID_TEXT_CANCEL_CHANGES,
                                        Str0[iLang][226]);
// i=227 "Valid text changes"
 TextValidChangesButton = new wxButton(m_TextPanel, ID_TEXT_VALID_CHANGES,
                                       Str0[iLang][227]);
// i=228 "Save to file"
 TextSaveToFileButton = new wxButton(m_TextPanel, ID_TEXT_SAVETOFILE,
                                         Str0[iLang][228]);

// Add buttons, horizontally unstretchable, with minimal size:
 hsizer3->Add( TextCancelChangesButton, 0);
 hsizer3->Add( TextValidChangesButton, 0, wxALIGN_RIGHT |wxLEFT | wxRIGHT, 20);
 hsizer3->Add( TextSaveToFileButton, 0, wxALIGN_RIGHT |wxLEFT | wxRIGHT, 20);

// Add button sizer with minimal size:
 topsizer->Add(hsizer3, 0, wxALIGN_RIGHT | wxALL, 20);

// Border at bottom:
 topsizer->Add(-1, 25);

// Sizer implementation on the panel:
  m_TextPanel->SetSizer(topsizer);      // use the sizer for layout

return(0);
}
/***********************************************************************
* Handle radio button to display Data File on OrbitDataTextCtrl 
************************************************************************/
void GpFrame::OnShowToTextPanel(wxCommandEvent &event)
{
bool show_buttons = true;

if(initialized != 1234) return;

 switch(event.GetId()){
// Display Data File
  case ID_SHOW_DATAFILE:
    show_buttons = false;
    m_textpanel_type = -2;
    if(!Data_fname.IsEmpty()) 
      OrbitDataTextCtrl->LoadFile(Data_fname);
    else
      OrbitDataTextCtrl->Clear();
    break;
// Display Orbit File
  case ID_SHOW_ORBITFILE:
    show_buttons = false;
    m_textpanel_type = -1;
    if(!Orbit_fname.IsEmpty()) 
      OrbitDataTextCtrl->LoadFile(Orbit_fname);
    else
      OrbitDataTextCtrl->Clear();
    break;
/*
* panel_type: 0: orbit and data  1: orbit only 2: data only 3: residuals 
*/
// Display orbit and data:
  case ID_SHOW_ORBITDATA:
    UpdatePanelText(0);
    break;
// Display orbit:
  case ID_SHOW_ORBIT:
    UpdatePanelText(1);
    break;
// Display data:
  case ID_SHOW_DATA:
    UpdatePanelText(2);
    break;
// Display residuals:
  case ID_SHOW_RESIDUALS:
    UpdatePanelText(3);
    break;
  }

// Show/Hide relevant buttons:
 TextCancelChangesButton->Show(show_buttons);
 TextValidChangesButton->Show(show_buttons);
 TextSaveToFileButton->Show(show_buttons);
//    OrbitDataTextCtrl->SetDefaultStyle(wxTextAttr(*wxRED));

}
/***********************************************************************
* 
************************************************************************/
void GpFrame::CorrectForPrecession()
{
JLP_Precession_Dlg *PrecessDlg;
double epoch_o, dtheta_precess;
int status, i;

  precession_corrected = 0;

  PrecessDlg = new JLP_Precession_Dlg(this, right_ascension1, declination1,
                                      orbit_equinox2, object_name1, 
                                      wxT("Correction for precession"));
  status = PrecessDlg->ShowModal();

// Retrieve the object/orbit parameters:
  PrecessDlg->RetrieveData(&right_ascension1, &declination1,
                           &orbit_equinox2, object_name1);

// If OK, return 0:
  if(status == 0 && orbit_equinox2 != 0. 
     && (right_ascension1 != 0. || declination1 != 0.)) {
   precession_corrected = 1;
   for(i = 0; i < nmeas1; i++) {
     epoch_o = epoch1[i];
     precession_correction(&dtheta_precess, right_ascension1, declination1,
                           epoch_o, orbit_equinox2);
     theta1[i] += dtheta_precess;
     }
   }

// Update object name, equinox, etc on OrbitPanel: 
  WidgetPanel_Update_OElements();

delete PrecessDlg;
}
/**************************************************************************
* Handle "TextCancelChanges" button:
**************************************************************************/
void GpFrame::OnTextCancelChanges( wxCommandEvent& WXUNUSED(event) )
{
 if(m_textpanel_type >= 0 && m_textpanel_type <= NTYPES_TEXTPANEL)
    UpdatePanelText(m_textpanel_type);
}
/**************************************************************************
* Handle "ValidOrbitChanges" button:
**************************************************************************/
void GpFrame::OnTextValidChanges( wxCommandEvent& WXUNUSED(event) )
{
wxString savefilename;

if(m_textpanel_type < 0 || m_textpanel_type > 2) return;

savefilename = wxFileSelector(
            wxT("Save to file (Tokovinin's format)"),
            wxT(""), wxT(""), wxT("dat|DAT|inp"),
            wxT("Files (*.dat;*.DAT;*.inp)|*.dat;*.DAT;*.inp"),
            wxFD_SAVE, this);

 if ( savefilename.empty() ) return;

// Save to file:
  OrbitDataTextCtrl->SaveFile(savefilename);

/*
* panel_type: 0: orbit and data  1: orbit only 2: data only 3: residuals 
*/
 switch(m_textpanel_type) {
   case 0:
     LoadTokoFile(savefilename);
     break;
   case 1:
// Format: 1 = Scardia, 2 = OC6, 3 = Tokovinin
     LoadOrbitFile(3, savefilename);
     break;
   case 2:
// Format: 1 = WDS, 2 = minimum format, 3 = Tokovinin, 4 = CORAVEL
     LoadDataFile(3, savefilename);
     break;
 }
}
/**************************************************************************
* Handle "TextSaveToFile" button:
**************************************************************************/
void GpFrame::OnTextSaveToFile( wxCommandEvent& WXUNUSED(event) )
{
wxString savefilename;

savefilename = wxFileSelector(
            wxT("Save to file (Tokovinin's format)"),
            wxT(""), wxT(""), wxT("dat|DAT|inp"),
            wxT("Files (*.dat;*.DAT;*.inp)|*.dat;*.DAT;*.inp"),
            wxFD_SAVE, this);

 if ( savefilename.empty() ) return;

// Save to file:
 OrbitDataTextCtrl->SaveFile(savefilename);

}
