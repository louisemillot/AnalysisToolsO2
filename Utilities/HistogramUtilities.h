#ifndef HISTOGRAM_UTILITIES_H
#define HISTOGRAM_UTILITIES_H

#include "TPolyLine.h"

// Histogram
float findMinFloat(float* array, int length);
float findMaxFloat(float* array, int length);

std::vector<double> GetTH1Bins(TH1 H1_histo);
TH2D RebinVariableBins2D(TH2D* H2D_hist, int nBinsX, int nBinsY, float* binsX, float* binsY);

TString contextCustom(TString* mainContext, TString* secondaryContext, TString* tertiaryContext, const char options[]);
TString contextCustom(TString* mainContext, TString* secondaryContext, const char options[]);
TString contextCustom(TString* mainContext, const char options[]);
TString contextEtaRange(float* EtaRange);
TString contextPtRange(float* PtRange);
TString contextJetRadius(float jetRadius);

void Draw_TH1_Histograms_in_one(TH1D** histograms_collection, const TString* legendList_string, int collectionSize, const TString Context, TString* pdfName, TString* &texXtitle, TString* &texYtitle, TString* texCollisionDataInfo, std::array<std::array<float, 2>, 2> drawnWindow, const char options[], TF1** optionalFitCollection);
void Draw_TH1_Histograms_in_one(TH1D** histograms_collection, const TString* legendList_string, int collectionSize, const TString Context, TString* pdfName, TString* &texXtitle, TString* &texYtitle, TString* texCollisionDataInfo, std::array<std::array<float, 2>, 2> drawnWindow, const char options[]);
void Draw_TH1_Histogram(TH1D* H1D_Sigma_asFunctionOf_Centrality, TString Context, TString* pdfName, TString* &texXtitle, TString* &texYtitle, TString* texCollisionDataInfo, std::array<std::array<float, 2>, 2> drawnWindow, const char options[]);
void Draw_TH2_Histograms(TH2D** histograms_collection, const TString* legendList_string, int collectionSize, const TString Context, TString* pdfName, TString* &texXtitle, TString* &texYtitle, TString* texCollisionDataInfo, std::array<std::array<float, 2>, 2> drawnWindow, const char options[]);
void Draw_TH2_Histograms(TH2D** histograms_collection, const TString* legendList_string, int collectionSize, TString Context, TString* pdfName, TString* &texXtitle, TString* &texYtitle, TString* texCollisionDataInfo, std::array<std::array<float, 2>, 2> drawnWindow, const char options[], TPolyLine** optionalLine) ;

// Preferred colors and markers
// const int fillColors[] = {kGray+1,  kRed-10, kBlue-9, kGreen-8, kMagenta-9, kOrange-9,kCyan-8,kYellow-7,}; // for syst bands
const int colors[]     = {kRed+1, kBlack, kBlue+1, kGreen+3, kMagenta+1, kOrange-1, kCyan+2, kYellow+2, kGray+1};
const int markers[]    = {kFullSquare, kFullCircle, kFullDiamond, kFullTriangleUp, kFullTriangleDown, kFullCross, kFullStar, kFullCrossX, kFullFourTrianglesX, kFullDoubleDiamond};

#endif
