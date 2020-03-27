/****************************************************************************
* Name: gpo_frame_id.h  (GpFrame class)
* list of ID used by GpFrame class
* 
* JLP
* Version 09/06/2015
****************************************************************************/
#ifndef _gpo_frame_id_h_
#define _gpo_frame_id_h_

// For compilers that support precompilation, includes "wx/wx.h".
#include "wx/wxprec.h"

#ifndef WX_PRECOMP
    #include "wx/wx.h"
#endif

//---------------------------------------------------------------------
//---------------------------------------------------------------------
enum{
//  ID_QUIT         = wxID_EXIT,
  ID_QUIT         = 1050, 

// Menu/File:
  ID_SAVE_TO_PST  = 1051,

  ID_LOAD_TOKO_FILE,
  ID_LOAD_CARQUI_FILE,
  ID_LOAD_TOKO_DATAFILE,
  ID_LOAD_CORAVEL_DATAFILE,
  ID_LOAD_WDS_DATAFILE,
  ID_LOAD_MINI_DATAFILE,

  ID_LOAD_TOKO_ORBITFILE,
  ID_LOAD_OC6_ORBITFILE,
  ID_LOAD_SCARDIA_ORBITFILE,

  ID_SAVE_ORBIT,
  ID_SAVE_DATA,
  ID_SAVE_ORBIT_DATA,
  ID_SAVE_RESIDUALS,

// Menu/Help context
  ID_CONTEXT_HELP,

// Menu/Notebook:
  ID_NOTEBOOK,
  ID_TEXT_PAGE,
  ID_GRAPHIC_PAGE,
  ID_WIDGET_PAGE,

// Text panel: 
  ID_SHOW_DATAFILE,
  ID_SHOW_ORBITFILE,
  ID_SHOW_DATA,
  ID_SHOW_ORBIT,
  ID_SHOW_ORBITDATA,
  ID_SHOW_RESIDUALS,
  ID_TEXT_VALID_CHANGES,
  ID_TEXT_CANCEL_CHANGES,
  ID_TEXT_SAVETOFILE,

// Widget panel: 
  ID_VALID_ORBIT_CHANGES,
  ID_CANCEL_ORBIT_CHANGES,
  ID_NEW_ORBIT_FIT,
  ID_DATA_FIT1,
  ID_DATA_FIT2,

// Widget panel:
  ID_TEXT_OEL0,
  ID_TEXT_OEL1,
  ID_TEXT_OEL2,
  ID_TEXT_OEL3,
  ID_TEXT_OEL4,
  ID_TEXT_OEL5,
  ID_TEXT_OEL6,

// Plot
  ID_PLOT_IDLE,
  ID_PLOT_RHO,
  ID_PLOT_RHO_RESID,
  ID_PLOT_THETA,
  ID_PLOT_THETA_RESID,
  ID_PLOT_ORBIT,
  ID_PLOT_RV_EPOCH,
  ID_PLOT_RV_PHASE,
  ID_PLOT_RV_RESID,
 
// Logbook
  ID_LOGBOOK_SHOW,
  ID_LOGBOOK_HIDE,
  ID_LOGBOOK_CLEAR,
  ID_LOGBOOK_CLEAN,
  ID_LOGBOOK_SAVE,

// Help:
  ID_ABOUT          = wxID_ABOUT,
  ID_HELP           = wxID_HELP
};

#endif
