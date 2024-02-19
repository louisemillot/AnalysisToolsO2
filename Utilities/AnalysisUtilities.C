#include "AnalysisUtilities.h"
#include "../Settings/GlobalSettings.h"

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////// Analysis Utilities ////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void NormaliseYieldToNEntries(TH1D* histogram) { 
  histogram->Scale(1./histogram->GetEntries(),"width"); // If option contains "width" the bin contents and errors are divided by the bin width.
}

void NormaliseYieldToNEvents(TH1D* histogram, int nEvents) { 
  histogram->Scale(1./nEvents,"width"); // If option contains "width" the bin contents and errors are divided by the bin width.
}

int GetNEventsSel8(TFile* file_O2Analysis) {
  return ((TH1I*)file_O2Analysis->Get("event-selection-task/hColCounterAcc"))->GetEntries();
}

int GetNEventsGen(TFile* file_O2Analysis) {
  return ((TH1I*)file_O2Analysis->Get("event-selection-task/hColCounterAcc"))->GetEntries();
}

int GetNEventsSel8Centrality(TFile* file_O2Analysis, float centralityLow, float centralityHigh) {
  TH1D* H1D_Centrality_FT0C= (TH1D*)file_O2Analysis->Get("centrality-qa/hCentFT0C");
  int iBinCent_low = H1D_Centrality_FT0C->GetXaxis()->FindBin(centralityLow + GLOBAL_epsilon);
  int iBinCent_high = H1D_Centrality_FT0C->GetXaxis()->FindBin(centralityHigh - GLOBAL_epsilon);
  return H1D_Centrality_FT0C->Integral(iBinCent_low, iBinCent_high);
}
