/****************************************************************************
* Name: GpOrbit1.cpp
* From Gplot1.cpp
* 
* JLP
* Version 23/06/2015
****************************************************************************/
#include <stdlib.h>   // exit() 
#include "time.h"

#if defined(__WXGTK__) || defined(__WXMOTIF__) || defined(__WXMAC__) || defined(__WXMGL__) || defined(__WXX11__)
    #define USE_XPM
#endif

#ifdef USE_XPM
    #include "mondrian.xpm"
#endif

#include "wx/progdlg.h"

#if !wxUSE_TOGGLEBTN
    #define wxToggleButton wxCheckBox
    #define EVT_TOGGLEBUTTON EVT_CHECKBOX
#endif

// JLP routines:
#include "gpo_frame.h"
#include "gpo_frame_id.h"  // My ID numbers (ID_xxx)
#include "jlp_wx_gpanel.h"
#include "gpo_language.h" // EACUTE, etc...

/*
#DEFINE DEBUG
*/

BEGIN_EVENT_TABLE(GpFrame, wxFrame)

// Notebook:
  EVT_NOTEBOOK_PAGE_CHANGED (ID_TEXT_PAGE, GpFrame::OnSelectPage)
  EVT_NOTEBOOK_PAGE_CHANGED (ID_GRAPHIC_PAGE, GpFrame::OnSelectPage)
  EVT_NOTEBOOK_PAGE_CHANGED (ID_WIDGET_PAGE, GpFrame::OnSelectPage)

// Text panel:
  EVT_RADIOBUTTON  (ID_SHOW_DATAFILE, GpFrame::OnShowToTextPanel)
  EVT_RADIOBUTTON  (ID_SHOW_ORBITFILE, GpFrame::OnShowToTextPanel)
  EVT_RADIOBUTTON  (ID_SHOW_DATA, GpFrame::OnShowToTextPanel)
  EVT_RADIOBUTTON  (ID_SHOW_ORBIT, GpFrame::OnShowToTextPanel)
  EVT_RADIOBUTTON  (ID_SHOW_ORBITDATA, GpFrame::OnShowToTextPanel)
  EVT_RADIOBUTTON  (ID_SHOW_RESIDUALS, GpFrame::OnShowToTextPanel)
  EVT_BUTTON  (ID_TEXT_CANCEL_CHANGES, GpFrame::OnTextCancelChanges)
  EVT_BUTTON  (ID_TEXT_VALID_CHANGES, GpFrame::OnTextValidChanges)
  EVT_BUTTON  (ID_TEXT_SAVETOFILE, GpFrame::OnTextSaveToFile)

// Widget panel:
  EVT_BUTTON  (ID_VALID_ORBIT_CHANGES, GpFrame::OnValidOrbitChanges)
  EVT_BUTTON  (ID_CANCEL_ORBIT_CHANGES, GpFrame::OnCancelOrbitChanges)
  EVT_BUTTON  (ID_NEW_ORBIT_FIT, GpFrame::OnNewOrbitFit)
  EVT_RADIOBUTTON (ID_DATA_FIT1, GpFrame::OnSelectDataFitMethod)
  EVT_RADIOBUTTON (ID_DATA_FIT2, GpFrame::OnSelectDataFitMethod)

// Menu/File:
  EVT_MENU(ID_LOAD_TOKO_FILE, GpFrame::OnLoadTokoFile)
  EVT_MENU(ID_LOAD_CARQUI_FILE, GpFrame::OnLoadCarquiFile)
// Load data
  EVT_MENU(ID_LOAD_TOKO_DATAFILE, GpFrame::OnLoadDataFile)
  EVT_MENU(ID_LOAD_CORAVEL_DATAFILE, GpFrame::OnLoadDataFile)
  EVT_MENU(ID_LOAD_WDS_DATAFILE, GpFrame::OnLoadDataFile)
  EVT_MENU(ID_LOAD_MINI_DATAFILE, GpFrame::OnLoadDataFile)
// Load orbit
  EVT_MENU(ID_LOAD_TOKO_ORBITFILE, GpFrame::OnLoadOrbitFile)
  EVT_MENU(ID_LOAD_OC6_ORBITFILE, GpFrame::OnLoadOrbitFile)
  EVT_MENU(ID_LOAD_SCARDIA_ORBITFILE, GpFrame::OnLoadOrbitFile)
// Save 
  EVT_MENU(ID_SAVE_ORBIT_DATA, GpFrame::OnSaveOrbitAndData)
  EVT_MENU(ID_SAVE_ORBIT, GpFrame::OnSaveOrbit)
  EVT_MENU(ID_SAVE_DATA, GpFrame::OnSaveData)
  EVT_MENU(ID_SAVE_RESIDUALS, GpFrame::OnSaveResiduals)
  EVT_MENU(ID_QUIT,           GpFrame::OnQuit)

// Menu/Plot
  EVT_MENU(ID_PLOT_IDLE, GpFrame::OnPlot)
  EVT_MENU(ID_PLOT_RHO, GpFrame::OnPlot)
  EVT_MENU(ID_PLOT_RHO_RESID, GpFrame::OnPlot)
  EVT_MENU(ID_PLOT_THETA, GpFrame::OnPlot)
  EVT_MENU(ID_PLOT_THETA_RESID, GpFrame::OnPlot)
  EVT_MENU(ID_PLOT_ORBIT, GpFrame::OnPlot)
  EVT_MENU(ID_PLOT_RV_EPOCH, GpFrame::OnPlot)
  EVT_MENU(ID_PLOT_RV_RESID, GpFrame::OnPlot)
  EVT_MENU(ID_PLOT_RV_PHASE, GpFrame::OnPlot)

// Menu/Logbook
  EVT_MENU(ID_LOGBOOK_SHOW, GpFrame::OnViewLogbook)
  EVT_MENU(ID_LOGBOOK_HIDE, GpFrame::OnViewLogbook)
  EVT_MENU(ID_LOGBOOK_CLEAN, GpFrame::OnCleanLogbook)
  EVT_MENU(ID_LOGBOOK_CLEAR, GpFrame::OnClearLogbook)
  EVT_MENU(ID_LOGBOOK_SAVE, GpFrame::OnSaveLogbook)

// Menu/Language
  EVT_MENU(ID_LANG_EN, GpFrame::OnSelectLanguage)
  EVT_MENU(ID_LANG_FR, GpFrame::OnSelectLanguage)
  EVT_MENU(ID_LANG_IT, GpFrame::OnSelectLanguage)
  EVT_MENU(ID_LANG_SP, GpFrame::OnSelectLanguage)

// Menu/Miscellaneous:
  EVT_MENU(ID_CONTEXT_HELP,   GpFrame::OnContextHelp)
  EVT_MENU(ID_ABOUT,          GpFrame::OnAbout)
  EVT_MENU(ID_HELP,           GpFrame::OnHelp)

END_EVENT_TABLE()

//----------------------------------------------------------------------
// MyApp
//----------------------------------------------------------------------

class MyApp: public wxApp
{
public:
   bool OnInit();
};

IMPLEMENT_APP(MyApp)

bool MyApp::OnInit()
{
// Transform coma into point for numbers:
setlocale(LC_NUMERIC, "C");

    // use standard command line handling:
    if ( !wxApp::OnInit() )
        return false;

    // parse the cmd line
    int iwidth = 800, iheight = 800;
    if ( argc == 3 )
    {
        wxSscanf(wxString(argv[1]), wxT("%d"), &iwidth);
        wxSscanf(wxString(argv[2]), wxT("%d"), &iheight);
    }

#if wxUSE_HELP
    wxHelpProvider::Set( new wxSimpleHelpProvider );
#endif // wxUSE_HELP

// Create the main frame window
    GpFrame *frame = new GpFrame(_T("GpOrbit1"), iwidth, iheight);

// Give it an icon
// The wxICON() macros loads an icon from a resource under Windows
// and uses an #included XPM image under GTK+ and Motif

#ifdef USE_XPM
    frame->SetIcon( wxICON(mondrian) );
#endif

    frame->Show(true);

    return true;
}

/**********************************************************************
* GpFrame constructor
*
* INPUT:
*   iwidth, iheight : size of created window
*
***********************************************************************/
GpFrame::GpFrame(const wxChar *title, int iwidth, int iheight)
       : wxFrame(NULL, wxID_ANY, title, wxPoint(-1, -1), wxSize(iwidth, iheight))
{
int i, iwidth1, iheight1;
wxSize m_size1;
wxString str1;
wxBoxSizer *w_hsizer0, *w_vsizer0;

// Initialize private variables:
  Data_fname = wxT("");
  Orbit_fname = wxT("");
  short_Orbit_fname = wxT("");
  short_Data_fname = wxT("");
  for(i = 0; i < NTYPES_TEXTPANEL; i++) TextPanel_txt[i] = wxT("");
// Original input data type as default:
  m_textpanel_type = -2;

// English language as default:
  iLang = 0;
// Load menu messages to Str0:
  LoadMenuMessages();

// Text panel as default in notebook:
  m_notebook = NULL;
  iPage = 0;

// nmeas1: number of (rho, theta) measurements read from input file
// nrv1, nrv2: number of radial velocity measurements read from input file
  nmeas1 = 0;
  nrv1 = 0;
  nrv2 = 0;
  precession_corrected = 0;

// Residuals:
  mean_sigma_rho_resid1 = 0.;
  mean_sigma_theta_resid1 = 0.;
  mean_sigma_rv1_resid1 = 0.;
  mean_sigma_rv2_resid1 = 0.;

// Status bar:
// Create a status bar with two fields at the bottom:
  m_StatusBar = CreateStatusBar(2);
// First field has a variable length, second has a fixed length:
  int widths[2];
  widths[0] = -1;
  widths[1] = 200;
  SetStatusWidths( 2, widths );

// Create an "empty" orbit oject:
  m_jlp_orbit1 = new JLP_Orbit1();
  nplot1 = 0;

//
  m_size1 = this->GetClientSize();
#ifdef DEBUG
  printf("GpFrame/Size: width=%d height=%d iwidth=%d iheight=%d\n", 
         m_size1.x, m_size1.y, iwidth, iheight);
#endif
  m_size1.x -= 20;
  m_size1.y -= 20;

// Create topsizer to locate panels and log window
  m_topsizer = new wxBoxSizer( wxVERTICAL );
// Create book control (multi-panels):
//  m_notebook = new wxBookCtrl(this, ID_BOOK_CTRL);
  m_notebook = new wxNotebook(this, ID_NOTEBOOK);

// Create Logbook panel first:
  str1 = wxString("");
  iwidth1 = m_size1.x;
  iheight1 = (int)((double)m_size1.y / 6.);
  LogPanel = new JLP_wxLogbook(this, str1, iwidth1, iheight1);

  wxLog::SetActiveTarget(new wxLogTextCtrl(LogPanel));

// Create text panel:
  m_TextPanel = new wxPanel(m_notebook, ID_TEXT_PAGE);
  TextPanel_Setup();
  m_notebook->AddPage(m_TextPanel, _T("Text Panel"));

// Then create the graphic panel:
  iwidth = m_size1.x;
  iheight = (int)((double)m_size1.y * 5. )/ 6.;
  m_GraphicPanel = new JLP_wxGraphicPanel(m_notebook, ID_GRAPHIC_PAGE,
                                          m_StatusBar, LogPanel, 20, 20, 
                                          iwidth, iheight);
  m_notebook->AddPage(m_GraphicPanel, _T("Graphic Panel"));

// Create widget panel:
  m_WidgetPanel = new wxPanel(m_notebook, ID_WIDGET_PAGE);
  WidgetPanel_Setup();
  m_notebook->AddPage(m_WidgetPanel, _T("Widget Panel"));

// The initial size of m_scrolled1 is interpreted as the minimal size:
// 1 : make vertically stretchable
// wxEXPAND : make horizontally stretchable, and the item will be expanded
// to fill the space assigned to the item.
// Proportion set to 5, i.e., graphic panel will be 5/6 of the window 
  m_topsizer->Add(m_notebook, 5, wxEXPAND | wxALL);

// Proportion set to 1, i.e., log window will be 1/6 of the window 
  m_topsizer->Add(LogPanel, 1, wxEXPAND);

// Sizer implementation on the panel:
SetSizerAndFit(m_topsizer);

// Create a menu on top of the window:
  Gp_SetupMenu();

initialized = 1234;

return;
}
/********************************************************************
* Setup the menu on top of main frame
********************************************************************/
void GpFrame::Gp_SetupMenu()
{
wxString sstr1[NLANG], sstr2[NLANG], help_str, item_str;

SetHelpText( _T("Program to plot curves from data contained in ASCII files") );

  menu_bar = new wxMenuBar;

// ***************** File menu **********************************
  wxMenu *file_menu = new wxMenu;

  menu_bar->Append(file_menu, _T("File"));

  wxMenu *loadfile_menu = new wxMenu;
  loadfile_menu->Append(ID_LOAD_TOKO_FILE, _T("Orbit and data (Tokovinin's fmt)"),
                    _T("Orbital elements, visual and radial velocity data"));
  loadfile_menu->Append(ID_LOAD_CARQUI_FILE, _T("Orbit and data (Carquillat's fmt)"),
                    _T("Orbital elements and radial velocity data"));
  loadfile_menu->AppendSeparator();
  loadfile_menu->Append(ID_LOAD_TOKO_DATAFILE, _T("Visual and RV data (Tokovinin's fmt)"),
                    _T("Visual and radial velocity data"));
  loadfile_menu->Append(ID_LOAD_CORAVEL_DATAFILE, _T("RV data (CORAVEL fmt)"),
                    _T("Radial velocity data"));
  loadfile_menu->Append(ID_LOAD_WDS_DATAFILE, _T("Visual data (WDS fmt)"),
  _T("WDS format: epoch rho theta nights observer diameter + weight_if_known"));
  loadfile_menu->Append(ID_LOAD_MINI_DATAFILE, _T("Visual data (mini. fmt)"),
                    _T("Format: epoch rho theta weight notes"));
  loadfile_menu->AppendSeparator();
  sstr1[0] = _T("Visual and RV orbit (Tokovinin's fmt)"); 
  sstr1[1] = _T("Orbite visuelle et spectro. (Format de Tokovinin)"); 
/*
  sstr1[2]= _T("Orbita visuale ") + EGRAVE + _T(" espettro. (Formato di Tokovinin)"); 
*/
  sstr1[2]= _T("Orbita visuale e espettro. (Formato di Tokovinin)"); 
  sstr1[3] = _T("Orbita visuale y espectro. (Formato de Tokovinin)"); 
  sstr2[0] = _T("in Tokovinin's format");
  sstr2[1] = _T("Format de Tokovinin");
  sstr2[2] = _T("Formato di Tokovinin");
  sstr2[3] = _T("Formato de Tokovinin");
  help_str = sstr2[0];
  loadfile_menu->Append(ID_LOAD_TOKO_ORBITFILE, sstr1[0], sstr2[0]); 
  loadfile_menu->Append(ID_LOAD_SCARDIA_ORBITFILE, 
                    _T("Visual orbit (Scardia)"),
                    _T("Format: W_Node w_peri i e T_peri P a equinox"));
  loadfile_menu->Append(ID_LOAD_OC6_ORBITFILE, 
                    _T("Visual orbit (OC6)"),
                    _T("OC6 format (WDS 6th Orbit Catalog)"));
  file_menu->Append(wxID_ANY, wxT("Load file"), loadfile_menu,
                    _T("Load data and orbit"));
  wxMenu *savefile_menu = new wxMenu;
  savefile_menu->Append(ID_SAVE_RESIDUALS, _T("Residuals"),
                    _T("Output file in Tokovinin's format"));
  savefile_menu->Append(ID_SAVE_ORBIT_DATA, _T("Orbit"),
                    _T("Output file in Tokovinin's format"));
  savefile_menu->Append(ID_SAVE_ORBIT_DATA, _T("Data"),
                    _T("Output file in Tokovinin's format"));
  savefile_menu->Append(ID_SAVE_ORBIT_DATA, _T("Orbit and data"),
                    _T("Output file in Tokovinin's format"));
  file_menu->Append(wxID_ANY, wxT("Save to file"), savefile_menu,
                    _T("Load data and orbit"));

  file_menu->Append(ID_QUIT, _T("Exit"), _T("Quit program"));

// ***************** Plot menu ******************************
  menuPlot = new wxMenu;
  menuPlot->Append( ID_PLOT_IDLE, _T("Idle"),
                       wxT("Clear graphic panel"), wxITEM_RADIO);
  menuPlot->Append( ID_PLOT_RHO, _T("rho vs epoch"),
                       wxT("Plot rho vs epoch"), wxITEM_RADIO);
  menuPlot->Append( ID_PLOT_THETA, _T("theta vs epoch"),
                       wxT("Plot theta vs epoch"), wxITEM_RADIO);
  menuPlot->Append( ID_PLOT_ORBIT, _T("orbit on sky plane"),
                       wxT("Plot orbit on sky plane"), wxITEM_RADIO);
  menuPlot->Append( ID_PLOT_RV_EPOCH, _T("rv vs epoch"),
                       wxT("Radial velocity vs epoch"), wxITEM_RADIO);
  menuPlot->Append( ID_PLOT_RV_PHASE, _T("rv vs phase"),
                       wxT("Radial velocity vs phase"), wxITEM_RADIO);
  menuPlot->Append( ID_PLOT_RHO_RESID, _T("rho residuals vs epoch"),
                       wxT("rho residuals vs epoch"), wxITEM_RADIO);
  menuPlot->Append( ID_PLOT_THETA_RESID, _T("theta residuals vs epoch"),
                       wxT("theta residuals vs epoch"), wxITEM_RADIO);
  menuPlot->Append( ID_PLOT_RV_RESID, _T("rv residuals vs epoch"),
                       wxT("Radial velocity residuals vs epoch"), wxITEM_RADIO);
  menu_bar->Append(menuPlot, _T("Plot"));

// ***************** Logbook menu ******************************
  menuLog = new wxMenu;
  menuLog->Append( ID_LOGBOOK_SHOW, _T("Show logbook"),
                       wxT("Display the logbook window"), wxITEM_RADIO);
  menuLog->Append( ID_LOGBOOK_HIDE, _T("Hide logbook"),
                       wxT("Hide the logbook window"), wxITEM_RADIO);
  menuLog->Append( ID_LOGBOOK_CLEAR, _T("Clear the logbook"),
                       wxT("Clear the logbook content"));
  menuLog->Append( ID_LOGBOOK_CLEAN, _T("Clean the logbook"),
                       wxT("Clean the logbook content"));
  menuLog->Append( ID_LOGBOOK_SAVE, _T("Save cleaned logbook"),
                       wxT("Save a selection from the logbook content"));
  menu_bar->Append(menuLog, _T("Logbook"));

// ***************** Language menu ******************************
  wxMenu *language_menu = new wxMenu;
  language_menu->Append(ID_LANG_EN, _T("English"), wxT(""), wxITEM_RADIO);
// ccedil: 00E7
// ntilde: 00F1
  language_menu->Append(ID_LANG_FR, _T("Fran\u00E7ais"), wxT(""), wxITEM_RADIO);
  language_menu->Append(ID_LANG_IT, _T("Italiano"), wxT(""), wxITEM_RADIO);
  language_menu->Append(ID_LANG_SP, _T("Espa\u00F1ol"), wxT(""), wxITEM_RADIO);
  menu_bar->Append(language_menu, _T("Language"));

// ***************** Help menu ******************************
  wxMenu *help_menu = new wxMenu;
  help_menu->Append(ID_HELP, _T("&Help"));
  help_menu->Append(ID_CONTEXT_HELP, _T("&Context help...\tCtrl-H"),
                     _T("Get context help for a control"));
  help_menu->Append(ID_ABOUT, _T("&About\tF1"));
  menu_bar->Append(help_menu, _T("&Help"));

  SetMenuBar(menu_bar);

return;
}

void GpFrame::OnQuit (wxCommandEvent& WXUNUSED(event) )
{
    Close(true);
}

/*****************************************************************
* Help 
*****************************************************************/
void GpFrame::OnHelp( wxCommandEvent& WXUNUSED(event) )
{
 (void)wxMessageBox(_T("Sorry: \"Help\" is not implemented yet\n") 
                    _T("Current version: June 2015"),
                    _T("GpOrbit1"),
                     wxICON_INFORMATION | wxOK );
}
/*****************************************************************
* About
*****************************************************************/
void GpFrame::OnAbout( wxCommandEvent& WXUNUSED(event) )
{
 (void)wxMessageBox( _T("GpOrbit1\n")
                     _T("Jean-Louis Prieur (c) 2015\n")
                     _T("Created with wxWidgets"), _T("GpOrbit1"), 
                     wxICON_INFORMATION | wxOK );
}
/*****************************************************************
* Context help
*****************************************************************/
void GpFrame::OnContextHelp(wxCommandEvent& WXUNUSED(event))
{
    // starts a local event loop
    wxContextHelp chelp(this);
}
/************************************************************************
** Display text in status bar 
*************************************************************************/
void GpFrame::SetText_to_StatusBar(wxString str1, const int icol)
{
// Update the first field (since 2nd argument is 0 here) of the status bar:
  if(m_StatusBar != NULL) m_StatusBar->SetStatusText(str1, icol);
}
