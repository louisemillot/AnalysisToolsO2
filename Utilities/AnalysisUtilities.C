#include "AnalysisUtilities.h"
#include "../Settings/GlobalSettings.h"

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////// Analysis Utilities ////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void NormaliseYieldToNEntries(TH1D* histogram) { // also takes care of the transformation of raw count to dCount/dQuantity (like dN/dpT), thanks to option "width"
  histogram->Scale(1./histogram->GetEntries(),"width"); // If option contains "width" the bin contents and errors are divided by the bin width.
  cout << "takes overflow/underflow into account; am I sure that is what I want?" << endl;
}

void NormaliseYieldToIntegral(TH1D* histogram) { // also takes care of the transformation of raw count to dCount/dQuantity (like dN/dpT), thanks to option "width"
  histogram->Scale(1./histogram->Integral("width"),"width"); // If option contains "width" the bin contents and errors are divided by the bin width.
}


void NormaliseYieldToNEvents(TH1D* histogram, int nEvents) { // also takes care of the transformation of raw count to dCount/dQuantity (like dN/dpT), thanks to option "width"
  histogram->Scale(1./nEvents,"width"); // If option contains "width" the bin contents and errors are divided by the bin width.
}
void NormaliseRawHistToNEvents(TH1D* histogram, int nEvents) { // also takes care of the transformation of raw count to dCount/dQuantity (like dN/dpT), thanks to option "width"
  histogram->Scale(1./nEvents,""); // If option contains "width" the bin contents and errors are divided by the bin width.
}

void TransformRawHistToYield(TH1D* histogram){
  histogram->Scale(1.,"width"); // If option contains "width" the bin contents and errors are divided by the bin width.
}

int GetNEventsSel8(TFile* file_O2Analysis) {
  return ((TH1I*)file_O2Analysis->Get("event-selection-task/hColCounterAcc"))->GetEntries(); //this is only sel8 (no sel8Full for example) and doesn't exclude collisions cut by the vertexZ pos cut
}

int GetNEventsSelected_JetFramework(TFile* file_O2Analysis) {
  return ((TH1I*)file_O2Analysis->Get("jet-finder-charged-qa/h_collisions"))->GetBinContent(2); //this should be the actual selection AND take vertexZ into account; sadly only works for jets
}
int GetNEventsSelected_TrackEffWorkflow(TFile* file_O2Analysis) {
  return ((TH1I*)file_O2Analysis->Get("track-efficiency/h_collisions"))->GetBinContent(2); //this should be the actual selection AND take vertexZ into account; sadly only works for jets
}

int GetNEventsSelected_JetFramework_weighted(TFile* file_O2Analysis) {
  return ((TH1I*)file_O2Analysis->Get("jet-finder-charged-qa/h_collisions_weighted"))->GetBinContent(2); //this should be the actual selection AND take vertexZ into account; sadly only works for jets
}


int GetNEventsGen(TFile* file_O2Analysis) {
  return ((TH1I*)file_O2Analysis->Get("event-selection-task/hColCounterAcc"))->GetEntries();
}

int GetNEventsSel8Centrality(TFile* file_O2Analysis, float centralityLow, float centralityHigh) { // probably same issue as getneventssel8
  TH1D* H1D_Centrality_FT0C= (TH1D*)file_O2Analysis->Get("centrality-qa/hCentFT0C");
  int iBinCent_low = H1D_Centrality_FT0C->GetXaxis()->FindBin(centralityLow + GLOBAL_epsilon);
  int iBinCent_high = H1D_Centrality_FT0C->GetXaxis()->FindBin(centralityHigh - GLOBAL_epsilon);
  return H1D_Centrality_FT0C->Integral(iBinCent_low, iBinCent_high);
}


int GetNEventsSelectedCentrality_JetFramework(TFile* file_O2Analysis, float centralityLow, float centralityHigh, const char trainId[]) { // should check it gives the correct number of coll (is posZ taken into account, fullsel8 etc)
  TH2D* H2D_Centrality_FT0C= (TH2D*)file_O2Analysis->Get("jet-finder-charged-qa"+(TString)trainId+"/h2_centrality_collisions");
  TH1D* H1D_Centrality_FT0C= (TH1D*)H2D_Centrality_FT0C->ProjectionX("H1D_Centrality_FT0C", 2, 2);

  int iBinCent_low = H1D_Centrality_FT0C->GetXaxis()->FindBin(centralityLow + GLOBAL_epsilon);
  int iBinCent_high = H1D_Centrality_FT0C->GetXaxis()->FindBin(centralityHigh - GLOBAL_epsilon);
  return H1D_Centrality_FT0C->Integral(iBinCent_low, iBinCent_high);
}