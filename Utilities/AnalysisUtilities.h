#ifndef ANALYSIS_UTILITIES_H
#define ANALYSIS_UTILITIES_H

// Analysis Utilities
void NormaliseYieldToNEntries(TH1D* histogram);
void NormaliseYieldToNEvents(TH1D* histogram, int nEvents);
void NormaliseRawHistToNEvents(TH1D* histogram, int nEvents);
void NormaliseYieldToIntegral(TH1D* histogram);
int GetNEvents(TFile* file_O2Analysis, TString analysisWorkflow);
int GetNEventsGen(TFile* file_O2Analysis);
int GetNEventsSelected_JetFramework(TFile* file_O2Analysis);
int GetNEventsSelected_JetFramework_weighted(TFile* file_O2Analysis);
int GetNEventsSel8Centrality(TFile* file_O2Analysis, float centralityLow, float centralityHigh);
int GetNEventsSelectedCentrality_JetFramework(TFile* file_O2Analysis, float centralityLow, float centralityHigh, const char trainId[]);

void TransformRawHistToYield(TH1D* histogram);

#endif
