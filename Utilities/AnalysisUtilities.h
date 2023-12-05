#ifndef ANALYSIS_UTILITIES_H
#define ANALYSIS_UTILITIES_H

// Analysis Utilities
void NormaliseYieldToNJets(TH1D* histogram);
void NormaliseYieldToNEvents(TH1D* histogram, int nEvents);
int GetNEvents(TFile* file_O2Analysis, TString analysisWorkflow);

#endif