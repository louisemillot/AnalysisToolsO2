#ifndef HISTOGRAM_UTILITIES_H
#define HISTOGRAM_UTILITIES_H

// Histogram
float findMinFloat(float* array, int length);
float findMaxFloat(float* array, int length);
void Draw_TH1_Histograms_in_one(TH1D** histograms_collection, const TString* legendList_string, Int_t collectionSize, const TString Context, TString* pdfName, TString* &texXtitle, TString* &texYtitle, TString* texCollisionDataInfo, const char options[]);
void Draw_TH2_Histograms(TH1D** histograms_collection, const TString* legendList_string, Int_t collectionSize, const TString Context, TString* pdfName, TString* &texXtitle, TString* &texYtitle, TString* texCollisionDataInfo, const char options[]);

// Preferred colors and markers
// const Int_t fillColors[] = {kGray+1,  kRed-10, kBlue-9, kGreen-8, kMagenta-9, kOrange-9,kCyan-8,kYellow-7,}; // for syst bands
const Int_t colors[]     = {kRed+1, kBlack, kBlue+1, kGreen+3, kMagenta+1, kOrange-1, kCyan+2, kYellow+2, kGray+1};
const Int_t markers[]    = {kFullSquare, kFullCircle, kFullDiamond, kFullTriangleUp, kFullTriangleDown, kFullCross, kFullStar, kFullCrossX, kFullFourTrianglesX, kFullDoubleDiamond};

#endif
