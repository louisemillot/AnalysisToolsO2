#ifndef HISTOGRAM_UTILITIES_H
#define HISTOGRAM_UTILITIES_H

#include "TPolyLine.h"

// Histogram
float findMinFloat(float* array, int length);
float findMaxFloat(float* array, int length);

std::vector<double> MakeVariableBinning_twoWidths(double xMin, int nLeft, double xMiddle, double xMax, int nRight);
std::vector<double> GetTH1Bins(TH1 H1_histo);

TH2D RebinVariableBins2D(TH2D* H2D_hist, int nBinsX, int nBinsY, double* binsX, double* binsY);
TH2D RebinVariableBins2D_PriorWeightedBinMerging(TH2D* H2D_hist, int nBinsX, int nBinsY, double* binsX, double* binsY, TH1D* H1D_priorSpectrum);
TH2D RebinVariableBins2D(TH2D* H2D_hist, int nBinsX, int nBinsY, double* binsX, double* binsY, bool debug);
TH2D RebinVariableBins2D_PriorWeightedBinMerging(TH2D* H2D_hist, int nBinsX, int nBinsY, double* binsX, double* binsY, TH1D* H1D_priorSpectrum, bool debug);
void NormaliseYSlicesToOne(TH2D* H2D_hist);
void NormaliseXSlicesToOne(TH2D* H2D_hist);
void NormaliseXSlicesToOneNoUnderOverFlows(TH2D* H2D_hist);
void WeightMatrixWithPrior(TH2D* H2D_hist, TH1D* priorSpectrum);
void TransformRawResponseToYieldResponse(TH2D* H2D_hist);


TH2* NormalizeResponsMatrixYaxisWithPrior(TH2 *h2RM, TH1 *hPrior);
TH2D* RebinVariableBins2D_aliPhysics(TH2D* hRMFine, int nBinsX, int nBinsY, double* binsX, double* binsY, bool useFunctionWeight);

// weird stuff happening here with optional debug thing, maybe remove it completely
TH2D GetTransposeHistogram(TH2D* inputHist);
TH2D GetMatrixProductTH2xTH2(TH2D* histA, TH2D* histB);
TH1D GetMatrixVectorProductTH2xTH1(TH2D* histA, TH1D* histU);

TString contextCustomThreeFields(TString mainContext, TString secondaryContext, TString tertiaryContext, std::string options);
TString contextCustomTwoFields(TString mainContext, TString secondaryContext, std::string options);
TString contextCustomOneField(TString mainContext, std::string options);
TString contextEtaRange(float* EtaRange);
TString contextPtRange(float* PtRange);
TString contextJetRadius(float jetRadius);

TString contextDatasetRadiusCompAndVarRange(TString mainContext, int iDataset, float* variableRange, std::string options);
TString contextDatasetCompAndRadiusAndVarRange(TString mainContext, float jetRadius, float* variableRange, std::string options);
TString contextDatasetCompAndRadius(TString mainContext, float jetRadius, std::string options);
TString contextDatasetComp(TString mainContext, std::string options);

void CentralityLegend(TString* centralityLegend, const float** arrayCentralityIntervals, int nCentralityBins);
// void IterationLegend(TString* iterationLegend, int nIterationmax);
void IterationLegend(TString* iterationLegend, int unfoldIterationMin, int unfoldIterationMax, int step);


void Draw_TH1_Histograms_in_one(TH1D** histograms_collection, const TString* legendList_string, int collectionSize, const TString Context, TString* pdfName, TString* &texXtitle, TString* &texYtitle, TString* texCollisionDataInfo, std::array<std::array<float, 2>, 2> drawnWindow, std::array<std::array<float, 2>, 2> legendPlacement, std::string options, TF1** optionalFitCollection);
void Draw_TH1_Histograms_in_one(TH1D** histograms_collection, const TString* legendList_string, int collectionSize, const TString Context, TString* pdfName, TString* &texXtitle, TString* &texYtitle, TString* texCollisionDataInfo, std::array<std::array<float, 2>, 2> drawnWindow, std::array<std::array<float, 2>, 2> legendPlacement, std::string options);
void Draw_TH1_Histogram(TH1D* H1D_Sigma_asFunctionOf_Centrality, TString Context, TString* pdfName, TString* &texXtitle, TString* &texYtitle, TString* texCollisionDataInfo, std::array<std::array<float, 2>, 2> drawnWindow, std::array<std::array<float, 2>, 2> legendPlacement, std::string options);
void Draw_TH2_Histograms(TH2D** histograms_collection, const TString* legendList_string, int collectionSize, const TString Context, TString* pdfName, TString* &texXtitle, TString* &texYtitle, TString* texCollisionDataInfo, std::array<std::array<float, 2>, 2> drawnWindow, std::string options);
void Draw_TH2_Histograms(TH2D** histograms_collection, const TString* legendList_string, int collectionSize, TString Context, TString* pdfName, TString* &texXtitle, TString* &texYtitle, TString* texCollisionDataInfo, std::array<std::array<float, 2>, 2> drawnWindow, std::string options, TPolyLine** optionalLine) ;
void Draw_TH2_Histogram(TH2D* histogram, TString Context, TString* pdfName, TString* &texXtitle, TString* &texYtitle, TString* texCollisionDataInfo, std::array<std::array<float, 2>, 2> drawnWindow, std::string options);

// Preferred colors and markers
// const int fillColors[] = {kGray+1,  kRed-10, kBlue-9, kGreen-8, kMagenta-9, kOrange-9,kCyan-8,kYellow-7,}; // for syst bands
const int colors[]               = {kRed+1, kBlack, kBlue+1, kGreen+3, kMagenta+1, kOrange-1, kCyan+2, kYellow+2, kGray+1};
const int markers[]              = {kFullSquare, kFullCircle, kFullDiamond, kFullTriangleUp, kFullTriangleDown, kFullCross, kFullStar, kFullCrossX, kFullFourTrianglesX, kFullDoubleDiamond};
const int markersColorPairs[]    = {kFullSquare, kOpenSquare, kFullCircle, kOpenCircle, kFullDiamond, kOpenDiamond, kFullTriangleUp, kOpenTriangleUp, kFullTriangleDown, kOpenTriangleDown};

#endif
