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
#include "jlp_language_dlg.h"
// #include "gpo_language.h" // EACUTE, etc...

#define DEBUG
/*
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
    int iwidth = 800, iheight = 500;
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
int i, iwidth1, iheight1, status;
wxSize m_size1;
wxString str1;
wxBoxSizer *w_hsizer0, *w_vsizer0;
wxString buffer;

// Initialize private variables:
  Data_fname = wxT("");
  Orbit_fname = wxT("");
  short_Orbit_fname = wxT("");
  short_Data_fname = wxT("");
  for(i = 0; i < NTYPES_TEXTPANEL; i++) TextPanel_txt[i] = wxT("");
// Original input data type as default:
  m_textpanel_type = -2;

// iLang: 0=English 1=French 2=Italian 3=Spanish 4=German
// Language as default:
// French
  iLang = 1;

// Prompt the user for a new value:
  SelectLanguageSetup();

// Load menu messages to Str0:
  status = LoadMenuMessages();
  if(status) {
    buffer.Printf(wxT("Error loading file GpOrbit1_messages.txt"));
    wxMessageBox(buffer, wxT("LoadmenuMessages"), wxOK | wxICON_ERROR);
    Close();
    }

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

  parallax1 = 127.;
  parallax_err1 = 2.;
  semi_majax1 = 0.; 
  semi_majax_err1 = 0.;
  total_mass1 = 0.; 
  total_mass_err1 = 0.;

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

// Initialize number of points:
  nplot1 = 0;

  plot_fname[0] = '\0';

// Initialize error values:
  for(i = 0; i < NPLOT_MAX; i++) {
    errorx1[i] = 0.;
    errory1[i] = 0.;
    }

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
// i=219 "Text Panel"
  m_notebook->AddPage(m_TextPanel, Str0[iLang][219]);

// Then create the graphic panel:
  iwidth = m_size1.x;
  iheight = (int)((double)m_size1.y * 5. )/ 6.;
/* OLD
  m_GraphicPanel = new JLP_wxGraphicPanel(m_notebook, ID_GRAPHIC_PAGE,
                                          m_StatusBar, LogPanel, 20, 20,
*/
  m_GraphicPanel = new JLP_wxGraphicPanel((wxFrame *)m_notebook, 
                                          ID_GRAPHIC_PAGE,
                                          LogPanel, 20, 20,
                                          iwidth, iheight);
// i=230 "Graphic Panel"
  m_notebook->AddPage(m_GraphicPanel, Str0[iLang][230]);

// Create widget panel:
  m_WidgetPanel = new wxPanel(m_notebook, ID_WIDGET_PAGE);
  WidgetPanel_Setup();

// i=250 "Widget Panel"
  m_notebook->AddPage(m_WidgetPanel, Str0[iLang][250]);

// The initial size of m_scrolled1 is interpreted as the minimal size:
// 1 : make vertically stretchable
// wxEXPAND : make horizontally stretchable, and the item will be expanded
// to fill the space assigned to the item.
// Proportion set to 4, i.e., graphic panel will be 5/6 of the window
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
* iLang: 0=English 1=French 2=Italian 3=Spanish 4=German
********************************************************************/
void GpFrame::Gp_SetupMenu()
{

SetHelpText( _T("Program to plot curves from data contained in ASCII files") );

  menu_bar = new wxMenuBar;

// ***************** File menu **********************************
  wxMenu *file_menu = new wxMenu;

// i=0: "File"
  menu_bar->Append(file_menu, Str0[iLang][0]);

  wxMenu *loadfile_menu = new wxMenu;
// i=1 "Orbit and data (Tokovinin's fmt)"
// i=2 "Orbital elements, visual and radial velocity data"
  loadfile_menu->Append(ID_LOAD_TOKO_FILE, Str0[iLang][1], Str0[iLang][2]);
// i=3 "Orbit and data (Carquillat's fmt)"
// i=4 "Orbital elements and radial velocity data"
  loadfile_menu->Append(ID_LOAD_CARQUI_FILE, Str0[iLang][3], Str0[iLang][4]);
  loadfile_menu->AppendSeparator();
// i=5 "Visual and RV data (Tokovinin's fmt)"
// i=6 "Visual and radial velocity data"
  loadfile_menu->Append(ID_LOAD_TOKO_DATAFILE, Str0[iLang][5], Str0[iLang][6]);
// i=7 "RV data (CORAVEL fmt)"
// i=8 "Radial velocity data"
  loadfile_menu->Append(ID_LOAD_CORAVEL_DATAFILE, Str0[iLang][7],
                        Str0[iLang][8]);
// i=9 "Visual data (WDS fmt)"
// i=10 "WDS format: epoch rho theta nights observer diameter [weight]"
  loadfile_menu->Append(ID_LOAD_WDS_DATAFILE, Str0[iLang][9], Str0[iLang][10]);
// i=11 "Visual data (mini. fmt)"
// i=12 "Format: epoch rho theta weight notes"
  loadfile_menu->Append(ID_LOAD_MINI_DATAFILE, Str0[iLang][11],
                        Str0[iLang][12]);
  loadfile_menu->AppendSeparator();
// i=13 "Visual and RV orbit (Tokovinin's fmt)"
// i=14 "Tokovinin's format"
  loadfile_menu->Append(ID_LOAD_TOKO_ORBITFILE, Str0[iLang][13],
                        Str0[iLang][14]);
// i=15 "Visual orbit (Scardia)"
// i=16 "Format: W_Node w_peri i e T_peri P a equinox"
  loadfile_menu->Append(ID_LOAD_SCARDIA_ORBITFILE, Str0[iLang][15],
                        Str0[iLang][16]);
// i=17 "Visual orbit (OC6)"
// i=18 "OC6 format (WDS 6th Orbit Catalog)"
  loadfile_menu->Append(ID_LOAD_OC6_ORBITFILE, Str0[iLang][17],
                        Str0[iLang][18]);
// i=19 "Load file"
// i=20 "Load data and orbit"
  file_menu->Append(wxID_ANY, Str0[iLang][19], loadfile_menu,
                    Str0[iLang][19]);
  wxMenu *savefile_menu = new wxMenu;
// i=21 "Residuals"
// i=22 "Output file in Tokovinin's format"
  savefile_menu->Append(ID_SAVE_RESIDUALS, Str0[iLang][21],
                    Str0[iLang][22]);
// i=23 "Orbit"
// i=24 "Output file in Tokovinin's format"
  savefile_menu->Append(ID_SAVE_ORBIT_DATA, Str0[iLang][23],
                        Str0[iLang][24]);
// i=25 "Data"
// i=26 "Output file in Tokovinin's format"
  savefile_menu->Append(ID_SAVE_ORBIT_DATA, Str0[iLang][25],
                        Str0[iLang][26]);
// i=27 "Orbit and Data"
// i=28 "Output file in Tokovinin's format"
  savefile_menu->Append(ID_SAVE_ORBIT_DATA, Str0[iLang][27],
                        Str0[iLang][28]);
// i=29 "Save to file"
// i=30 "Load data and orbit"
  file_menu->Append(wxID_ANY, Str0[iLang][29], savefile_menu,
                    Str0[iLang][30]);
// i=31 "Exit"
// i=32 "Quit program"
  file_menu->Append(ID_QUIT, Str0[iLang][31], Str0[iLang][32]);

// ***************** Plot menu ******************************
  menuPlot = new wxMenu;
// i=40 "Idle"
// i=41 "Clear graphic panel"
  menuPlot->Append( ID_PLOT_IDLE, Str0[iLang][40], Str0[iLang][41],
                    wxITEM_RADIO);
// i=42 "rho vs epoch"
// i=43 "Plot rho vs epoch"
  menuPlot->Append( ID_PLOT_RHO, Str0[iLang][42], Str0[iLang][43],
                    wxITEM_RADIO);
// i=44 "theta vs epoch"
// i=45 "Plot theta vs epoch"
  menuPlot->Append( ID_PLOT_THETA, Str0[iLang][44], Str0[iLang][45],
                    wxITEM_RADIO);
// i=46 "orbit on sky plane"
// i=47 "Plot orbit on sky plane"
  menuPlot->Append( ID_PLOT_ORBIT, Str0[iLang][46], Str0[iLang][47],
                    wxITEM_RADIO);
// i=48 "rv vs epoch"
// i=49 "Radial velocity vs epoch"
  menuPlot->Append( ID_PLOT_RV_EPOCH, Str0[iLang][48], Str0[iLang][49],
                    wxITEM_RADIO);
// i=50 "rv vs phase"
// i=51 "Radial velocity vs phase"
  menuPlot->Append( ID_PLOT_RV_PHASE, Str0[iLang][50], Str0[iLang][51],
                    wxITEM_RADIO);
// i=52 "rho residuals vs epoch"
// i=53 "rho residuals vs epoch"
  menuPlot->Append( ID_PLOT_RHO_RESID, Str0[iLang][52], Str0[iLang][53],
                    wxITEM_RADIO);
// i=54 "theta residuals vs epoch"
// i=55 "theta residuals vs epoch"
  menuPlot->Append( ID_PLOT_THETA_RESID, Str0[iLang][54], Str0[iLang][55],
                    wxITEM_RADIO);
// i=56 "rv residuals vs epoch"
// i=57 "Radial velocity residuals vs epoch"
  menuPlot->Append( ID_PLOT_RV_RESID, Str0[iLang][56], Str0[iLang][57],
                    wxITEM_RADIO);
// i=58 "Plot"
  menu_bar->Append(menuPlot, Str0[iLang][58]);

// ***************** Logbook menu ******************************
  menuLog = new wxMenu;
// i=100 "Show logbook"
// i=101 "Display the logbook window"
  menuLog->Append( ID_LOGBOOK_SHOW, Str0[iLang][100], Str0[iLang][101],
                    wxITEM_RADIO);
// i=102 "Hide logbook"
// i=103 "Hide the logbook window"
  menuLog->Append( ID_LOGBOOK_HIDE, Str0[iLang][102], Str0[iLang][103],
                    wxITEM_RADIO);
// i=104 "Clear the logbook"
// i=105 "Clear the logbook content"
  menuLog->Append( ID_LOGBOOK_CLEAR, Str0[iLang][104], Str0[iLang][105],
                    wxITEM_RADIO);
// i=106 "Clean the logbook"
// i=107 "Clean the logbook content"
  menuLog->Append( ID_LOGBOOK_CLEAN, Str0[iLang][106], Str0[iLang][107],
                    wxITEM_RADIO);
// i=108 "Save cleaned logbook"
// i=109 "Save a selection from the logbook content"
  menuLog->Append( ID_LOGBOOK_SAVE, Str0[iLang][108], Str0[iLang][109],
                    wxITEM_RADIO);
// i=110 "Logbook"
  menu_bar->Append(menuLog, Str0[iLang][110]);

// ***************** Help menu ******************************
// i=200 "Help"
// i=201 "Help about this program"
  wxMenu *help_menu = new wxMenu;
  help_menu->Append(ID_HELP, Str0[iLang][200], Str0[iLang][201]);
// i=202 "&Context help...\tCtrl-H"
// i=203 "Get context help for a control"
  help_menu->Append(ID_CONTEXT_HELP, Str0[iLang][202], Str0[iLang][203]);
// i=204 "About"
// i=205 "About this program"
  help_menu->Append(ID_ABOUT, Str0[iLang][204], Str0[iLang][205]);
// i=206 "Help"
  menu_bar->Append(help_menu, Str0[iLang][206]);

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
wxString sstr;

// i=207 "Sorry: \"Help\" is not implemented yet\n"
// i=208 "Current version: June 2015"
 sstr = Str0[iLang][207] + "\n" + Str0[iLang][208];
 (void)wxMessageBox( sstr, _T("GpOrbit1"), wxICON_INFORMATION | wxOK );
}
/*****************************************************************
* About
*****************************************************************/
void GpFrame::OnAbout( wxCommandEvent& WXUNUSED(event) )
{
wxString sstr;

// i=209 "Created with wxWidgets"
 sstr = _T("GpOrbit1\nJean-Louis Prieur (c) 2015\n") + Str0[iLang][207];
 (void)wxMessageBox( sstr, _T("GpOrbit1"), wxICON_INFORMATION | wxOK );
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
/**************************************************************************
*
***************************************************************************/
void GpFrame::SelectLanguageSetup()
{
JLP_Language_Dlg *LanguageDlg;
int status, i, i_lang;

  LanguageDlg = new JLP_Language_Dlg(this, wxT("Language Selection"));

  status = LanguageDlg->ShowModal();

// Retrieve the object/orbit parameters:
  LanguageDlg->RetrieveData(&i_lang);

// Set private variable iLang if status is OK:
  if(status == 0 && i_lang >= 0 && i_lang < NLANG) {
   iLang = i_lang;
   }

delete LanguageDlg;

// Exit if Cancel
if(status) {
 exit(-1);
}

return;
}
