/****************************************************************************
* Name: gpo_frame.h
* 
* JLP
* Version 08/06/2015
****************************************************************************/
#ifndef _gpo_frame__ 
#define _gpo_frame__ 

// For compilers that support precompilation, includes "wx/wx.h".
#include "wx/wxprec.h"

#ifndef WX_PRECOMP
    #include "wx/wx.h"
#endif

#include "wx/tglbtn.h"
#include "wx/bookctrl.h"
#include "wx/imaglist.h"
#include "wx/cshelp.h"

#if wxUSE_TOOLTIPS
    #include "wx/tooltip.h"
#endif

#include "gpo_defs.h"         // NMEAS_MAX, NELEMENTS, N7, ... 
#include "jlp_wxlogbook.h"    // JLP_wxLogbook class
#include "jlp_wx_gpanel.h"    // JLP_wxGraphicPanel class
#include "jlp_orbit1.h"       // JLP_Orbit1 class

//----------------------------------------------------------------------
// class definitions
//----------------------------------------------------------------------

class GpFrame: public wxFrame
{
public:
    GpFrame(const wxChar *title, int x, int y);
    ~GpFrame() {return;};

    void Gp_SetupMenu();
    void OnQuit(wxCommandEvent& event);
    void OnAbout(wxCommandEvent& event);
    void OnHelp(wxCommandEvent& event);
    void SelectLanguageSetup();


    void SetText_to_StatusBar(wxString text, const int icol);

// In "gpo_menu.cpp":
    void OnPlot(wxCommandEvent& event);
    void OnLoadTokoFile(wxCommandEvent& event);
    void LoadTokoFile(wxString &fname);
    void OnLoadCarquiFile(wxCommandEvent& event);
    void OnLoadDataFile(wxCommandEvent& event);
    void LoadDataFile(int iformat, wxString &fname);
    void OnLoadOrbitFile(wxCommandEvent& event);
    void LoadOrbitFile(int iformat, wxString &fname);
    void OnSaveOrbit(wxCommandEvent& event);
    void OnSaveData(wxCommandEvent& event);
    void OnSaveOrbitAndData(wxCommandEvent& event);
    void OnSaveResiduals(wxCommandEvent& event);
    void ResetBeforeLoadingDataFile();
    void ResetBeforeLoadingOrbitFile();
    void UpdateMenuAfterFileIsLoaded();
    void UpdateDataAfterFileIsLoaded();
    void UpdateElementsAfterFileIsLoaded(double *oelmnts0, double *oelmnts_err0,
                                       const int max_noelmnts0,
                                      int orbit_type0);
    void UpdateOrbitPanelAfterFileIsLoaded();
    void UpdateDataFitPanelAfterFileIsLoaded();
    int SetOrbitType();

// In "gpo_rw_tokofiles.cpp":
    int save_to_toko_file(const char *out_orbit_data_file);
    int UpdatePanelText(const int textpanel_type);
    int PrepareOrbitDataText(wxString &Buffer_txt);
    int PrepareOrbitText(wxString &Buffer_txt);
    int PrepareDataText(wxString &Buffer_txt);
    int PrepareResidualsText(wxString &Buffer_txt);

// In "gp_frame_menu.cpp":
    void OnViewLogbook(wxCommandEvent& event);
    void OnSaveLogbook(wxCommandEvent& event);
    void OnClearLogbook(wxCommandEvent& event);
    void OnCleanLogbook(wxCommandEvent& event);

// In "gpo_text_panel.cpp":
    int TextPanel_Setup(); 
    void OnShowToTextPanel(wxCommandEvent& event);
    void CorrectForPrecession();

// In "gpo_widget_panel.cpp":
    int WidgetPanel_Setup(); 
    int WidgetPanel_OElmnts_str_Setup(const int orbit_type0);
    int WidgetPanel_OElmnts_Setup();
    int WidgetPanel_Results_Setup();
    int WidgetPanel_Update_OElements();
    int WidgetPanel_Update_Results();
    int ComputeMassAndSemiAxis(double *oelmnts0, double *oelmnts_err0);
    int WidgetPanel_Display_OElements(double *oelmnts0, double *oelmnts_err0,
                                      int *oelmnts_fflags0, 
                                      const int noelmnts_max);

// In "gpo_language.cpp":
    void OnSelectLanguage(wxCommandEvent& event);
    int LoadMenuMessages();

// In "gpo_plot.cpp":
    int plot_rv_vs_epoch();
    int plot_rv_resid_vs_epoch();
    int plot_rv_vs_phase();
    int plot_rho_vs_epoch();
    int plot_rho_resid_vs_epoch();
    int plot_theta_vs_epoch();
    int plot_theta_resid_vs_epoch();
    int plot_orbit_skyplane();
    int compute_sense_of_motion(int *north_to_east);
    int draw_line_of_absids();
    int draw_central_cross(double user_cross_width);
    int draw_north_east_label(double cross_width);

// in "gpo_logbook.cpp":
    int SaveLogbook(wxString save_filename);
    void ShowLogbook();
    void HideLogbook();
    void ClearLogbook();
    void CleanLogbook();
    int WriteToLogbook(wxString str1, bool SaveToFile);
    int AddNewPointToLogbook(double xx, double yy, double value);

#if wxUSE_TOOLTIPS
    void OnToggleTooltips(wxCommandEvent& event);
#endif // wxUSE_TOOLTIPS

    void OnContextHelp(wxCommandEvent& event);

private:
  void OnResize(wxSizeEvent &event);

// Notebook:
  void OnSelectPage(wxBookCtrlEvent& event);

// Text panel:
  void OnTextCancelChanges(wxCommandEvent& event);
  void OnTextValidChanges(wxCommandEvent& event);
  void OnTextSaveToFile(wxCommandEvent& event);

// Widget panel:
  void OnCancelOrbitChanges(wxCommandEvent& event);
  void OnValidOrbitChanges(wxCommandEvent& event);
  void OnNewOrbitFit(wxCommandEvent& event);
  void OnSelectDataFitMethod(wxCommandEvent& event);

  JLP_wxGraphicPanel *m_GraphicPanel;
  JLP_Orbit1 *m_jlp_orbit1;

  int initialized, iLang, iPage;

  wxString Orbit_fname, short_Orbit_fname;
  wxString Data_fname, short_Data_fname;
  wxString TextPanel_txt[NTYPES_TEXTPANEL]; 

// Private arrays: 
// epoch1, theta1, rho1, sigma_rho1: arrays filled with the measurement values
  double epoch1[NMEAS_MAX], theta1[NMEAS_MAX], rho1[NMEAS_MAX]; 
  double sigma_rho1[NMEAS_MAX];
  double epoch_rv1[NMEAS_MAX], epoch_rv2[NMEAS_MAX];
  double rv1[NMEAS_MAX], rv2[NMEAS_MAX];
  double sigma_rv1[NMEAS_MAX], sigma_rv2[NMEAS_MAX];
  char notes_meas1[NMEAS_MAX][128]; 
  char notes_rv1[NMEAS_MAX][128], notes_rv2[NMEAS_MAX][128];
  char object_name1[128], plot_fname[128]; 
  double right_ascension1, declination1;

// nmeas1: number of (rho, theta) measurements read from input file
// nrv1, nrv2: number of radial velocity measurements read from input file
  int nmeas1, nrv1, nrv2;

  double xplot1[NPLOT_MAX], yplot1[NPLOT_MAX], yplot2[NPLOT_MAX];
  double errorx1[NPLOT_MAX], errory1[NPLOT_MAX];
  int nplot1;
  double mean_sigma_rv1_resid1, mean_sigma_rv2_resid1; 
  double mean_rho_resid1, mean_theta_resid1;
  double mean_sigma_rho_resid1, mean_sigma_theta_resid1;
  double parallax1, parallax_err1, semi_majax1, semi_majax_err1; 
  double total_mass1, total_mass_err1;

// Orbital parameters: 
  double orbit_equinox2;
  int precession_corrected;

// Menus:
  wxMenuBar *menu_bar;
  wxMenu *menuFile, *menuLog, *menuPlot;
  wxBoxSizer  *m_topsizer;
  wxStatusBar *m_StatusBar;
  wxString Str0[NLANG][NMAX_MESSAGES];
 
// Notebook:
  wxBookCtrl *m_notebook;

// Logbook:
  wxString    m_Logbook;
  JLP_wxLogbook *LogPanel;

// TextPanel:
  wxPanel *m_TextPanel;
  wxButton *TextValidChangesButton, *TextCancelChangesButton;
  wxButton *TextSaveToFileButton;
  wxTextCtrl *OrbitDataTextCtrl;
  wxRadioButton *ShowOrbitDataButton, *ShowDataFileButton, *ShowOrbitFileButton;
  wxRadioButton *ShowDataButton, *ShowOrbitButton, *ShowResidualsButton; 
  int m_textpanel_type;

// WidgetPanel:
  wxPanel *m_WidgetPanel;
  wxButton *ValidOrbitChangesButton, *CancelOrbitChangesButton;
  wxButton *NewFitButton;
  wxFlexGridSizer *fgs1, *fgs2;
  wxString oelmt_str[NELEMENTS];
  wxStaticBoxSizer *oelmt_sizer, *results_sizer;
  wxStaticText *oelmt_name_static[NELEMENTS], *oelmt_err_static[NELEMENTS];
  wxStaticText *result_name_static[NRESULTS_MAX]; 
  wxStaticText *result_value_static[NRESULTS_MAX];
  wxStaticText *result_units_static[NRESULTS_MAX];
  wxTextCtrl *oelmt_txt_ctrl[NELEMENTS];
  wxCheckBox *oelmt_fflag_checkbox[NELEMENTS];
  wxRadioButton *DataFitButton1, *DataFitButton2; 
  int nresults1, data_fit_method1;

  DECLARE_EVENT_TABLE()
};

#endif
