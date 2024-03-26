#ifndef ANALYSIS_UTILITIES_H
#define ANALYSIS_UTILITIES_H

// Analysis Utilities
void NormaliseYieldToNEntries(TH1D* histogram);
void NormaliseYieldToNEvents(TH1D* histogram, int nEvents);
void NormaliseYieldToIntegral(TH1D* histogram);
int GetNEvents(TFile* file_O2Analysis, TString analysisWorkflow);
int GetNEventsSel8Centrality(TFile* file_O2Analysis, float centralityLow, float centralityHigh);

#endif
