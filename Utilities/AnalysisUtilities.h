#ifndef ANALYSIS_UTILITIES_H
#define ANALYSIS_UTILITIES_H

// Analysis Utilities
void NormaliseYieldToNEntries(TH1D* histogram);
void NormaliseRawHistToNEntries(TH1D* histogram);
void NormaliseYieldToNEvents(TH1D* histogram, double nEvents);
void NormaliseRawHistToNEvents(TH1D* histogram, double nEvents);
void NormaliseYieldToIntegral(TH1D* histogram);
void NormaliseRawHistToIntegral(TH1D* histogram);
int GetNEvents(TFile* file_O2Analysis, TString analysisWorkflow);
int GetNEventsGen(TFile* file_O2Analysis);
int GetNEventsSelected_JetFramework(TFile* file_O2Analysis, const char analysisWorkflow[]);
double GetNEventsSelected_JetFramework_weighted(TFile* file_O2Analysis, const char analysisWorkflow[]);
int GetNEventsSelected_TrackEffWorkflow(TFile* file_O2Analysis, const char analysisWorkflow[]);
int GetNEventsSelected_TrackEffWorkflow_gen(TFile* file_O2Analysis, const char analysisWorkflow[]);
double GetNEventsSelected_TrackEffWorkflow_weighted(TFile* file_O2Analysis, const char analysisWorkflow[]);
double GetNEventsSelected_TrackEffWorkflow_gen_weighted(TFile* file_O2Analysis, const char analysisWorkflow[]);
int GetNEventsSel8Centrality(TFile* file_O2Analysis, float centralityLow, float centralityHigh);
int GetNEventsSelectedCentrality_JetFramework(TFile* file_O2Analysis, float centralityLow, float centralityHigh, const char wagonId[]);

void TransformRawHistToYield(TH1D* histogram);

#endif
