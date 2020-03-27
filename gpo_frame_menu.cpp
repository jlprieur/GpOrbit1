/******************************************************************************
* Name:        gpo_frame_menu (GpFrame class)
* Purpose:     handling menu events of GpFrame class
*
* Author:      JLP 
* Version:     08/06/2015 
******************************************************************************/
#include "gpo_frame.h"
#include "gpo_frame_id.h"  // Menu identifiers

/* Declared in gpo_frame.h

void OnSaveToPostscript(wxCommandEvent &event);
void OnViewLogbook(wxCommandEvent& event);
void OnSaveLogbook(wxCommandEvent& WXUNUSED(event));
void OnClearLogbook(wxCommandEvent& event);
void OnCleanLogbook(wxCommandEvent& event);

*/

/************************************************************************
* Showing/hiding logbook panel 
************************************************************************/
void GpFrame::OnViewLogbook(wxCommandEvent& event)
{
  switch (event.GetId())
  {
   case ID_LOGBOOK_SHOW:
     ShowLogbook();
     break;
   case ID_LOGBOOK_HIDE:
     HideLogbook();
     break;
   }
}
/************************************************************************
* Save useful content of logbook to file 
************************************************************************/
void GpFrame::OnSaveLogbook(wxCommandEvent& WXUNUSED(event))
{
wxString save_filename;

// Select name for output logbook file:
wxFileDialog
saveFileDialog(this, wxT("Save logbook to file"), wxT(""), wxT(""), 
               wxT("Logbook files (*.log;*.txt)|*.log;*.txt"), 
               wxFD_SAVE|wxFD_OVERWRITE_PROMPT);

 if (saveFileDialog.ShowModal() == wxID_CANCEL) return;

save_filename = saveFileDialog.GetFilename();

SaveLogbook(save_filename);

return;
}
/*******************************************************************
* Clear the logbook: erase all its content
********************************************************************/
void GpFrame::OnClearLogbook(wxCommandEvent& event)
{
ClearLogbook();
return;
}
/*******************************************************************
* Clean the logbook: only keep its useful content
********************************************************************/
void GpFrame::OnCleanLogbook(wxCommandEvent& event)
{
CleanLogbook();
}
