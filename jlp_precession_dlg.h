/******************************************************************************
* jlp_precession_dlg.h
* To apply precession correction on data 
*
* Author:      JLP 
* Version:     28/07/2015
******************************************************************************/
#ifndef jlp_precession_dlg_h    // sentry 
#define jlp_precession_dlg_h

// For compilers that support precompilation, includes "wx/wx.h".
#include "wx/wxprec.h"

#ifndef WX_PRECOMP
    #include "wx/wx.h"
#endif

#include "wx/filename.h"

/********************************************************************
* Class JLP_Precession_Dlg
*********************************************************************/

class JLP_Precession_Dlg: public wxDialog
{
public:

// Constructor:
     JLP_Precession_Dlg(wxFrame *parent, double RA_0, double DEC_0,
                   double Equinox_0, char *object_name0, 
                   const wxString &title);
// Destructor: 
    ~JLP_Precession_Dlg(){
       };

// Handling events:
     void OnOKButton( wxCommandEvent &event );
     void OnCancelButton( wxCommandEvent &event );
     void OnChangeParam( wxCommandEvent& event );

// Accessors:
   int RetrieveData(double *RA_2, double *DEC_2, double *Equinox_2, 
                    char *object_name2) {
       *RA_2 = RA_1; 
       *DEC_2 = DEC_1; 
       *Equinox_2 = Equinox_1; 
       strcpy(object_name2, object_name1);
       return(0);
       }

protected:
   bool DataIsOK() {
       if(RA_1 < 0. || RA_1 > 24. 
          || DEC_1 < -90. || DEC_1 > 90.) return(false); 
       if(Equinox_1 < 0. || Equinox_1 > 3000.) return(false); 
       return(true);
       }

private:
    wxTextCtrl *TextCtrl_Object, *TextCtrl_RA, *TextCtrl_DEC;
    wxTextCtrl *TextCtrl_Equinox; 
    char object_name1[128];
    double RA_1, DEC_1, Equinox_1;
    int initialized;

    DECLARE_EVENT_TABLE()
};

#endif               // EOF sentry
