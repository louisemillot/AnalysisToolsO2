#include "TStyle.h"
#include "TFile.h"
#include "TFitResult.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TF1.h"
#include "TLatex.h"
#include "TLine.h"
#include "TRatioPlot.h"
#include "TLegend.h"
#include "TH1.h"
// #include <RooUnfold.h>
// #include "RooUnfoldResponse.h"
// #include "RooUnfoldBayes.h"
// #include "RooUnfoldBinByBin.h"

//My Libraries
#include "./TrackQC_settings.h"
#include "../Settings/AxisTitles.h"
#include "../Utilities/AnalysisUtilities.h"
#include "../Utilities/HistogramUtilities.h"
#include "../Utilities/AnalysisUtilities.C" // bizarre but only including the .h fils doesn't work for the standard 'root macro.C+' method, I need to include the .C as well
#include "../Utilities/HistogramUtilities.C" // bizarre but only including the .h fils doesn't work for the standard 'root macro.C+' method, I need to include the .C as well

#include<array>
#include <iomanip>
#include <sstream>
#include <string.h>
using namespace std;

// Misc utilities
void SetStyle(Bool_t graypalette=kFALSE);
void LoadLibs();

// Plot Utilities
TString contextTrackDatasetComp(const char options[]);

//////////// QC plot functions
// Dataset comparison
void Draw_Pt_DatasetComparison();
void Draw_Eta_DatasetComparison();
void Draw_Phi_DatasetComparison();


/////////////////////////////////////////////////////
///////////////////// Main Macro ////////////////////
/////////////////////////////////////////////////////

void TrackQC() {
  // Load necessary libraries
  LoadLibs();
  // Set the default style
  SetStyle();

  // TString* SaveAs_Title = new TString("");
  TString* texXtitle = new TString("");
  TString* texYtitle = new TString("");
  // TString* Extra = new TString("");


  Draw_Pt_DatasetComparison();
  Draw_Eta_DatasetComparison();
  Draw_Phi_DatasetComparison();

}

/////////////////////////////////////////////////////
/////////////////// Misc utilities //////////////////
/////////////////////////////////////////////////////

void LoadLibs() {
  // gSystem->Load("libCore.so");  
  // gSystem->Load("libGeom.so");
  // gSystem->Load("libPhysics.so");
  // gSystem->Load("libVMC");
  // gSystem->Load("libTree");
  // gSystem->Load("libMinuit");
  // gSystem->Load("libSTEERBase");
  // gSystem->Load("libESD");
  // gSystem->Load("libAOD");
  // gSystem->Load("libANALYSIS");
  // gSystem->Load("libANALYSISalice");
  // gSystem->Load("libCORRFW");
  // gSystem->Load("libPWGTools");
}

void SetStyle(Bool_t graypalette) {
  cout << "Setting style!" << endl;
  
  gStyle->Reset("Plain");
  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);
  if(graypalette) gStyle->SetPalette(8,0);
  else gStyle->SetPalette(1);
  gStyle->SetCanvasColor(10);
  gStyle->SetCanvasBorderMode(0);
  gStyle->SetFrameLineWidth(1);
  gStyle->SetFrameFillColor(kWhite);
  gStyle->SetPadColor(10);
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);
  gStyle->SetPadBottomMargin(0.15);
  gStyle->SetPadLeftMargin(0.15);
  gStyle->SetHistLineWidth(1);
  gStyle->SetHistLineColor(kRed);
  gStyle->SetFuncWidth(2);
  gStyle->SetFuncColor(kGreen);
  gStyle->SetLineWidth(2);
  gStyle->SetLineScalePS(1);
  gStyle->SetLabelSize(0.045,"xyz");
  gStyle->SetLabelOffset(0.01,"y");
  gStyle->SetLabelOffset(0.01,"x");
  gStyle->SetLabelColor(kBlack,"xyz");
  gStyle->SetTitleSize(0.05,"xyz");
  gStyle->SetTitleOffset(1.4,"y");
  gStyle->SetTitleOffset(1.2,"x");
  gStyle->SetTitleFillColor(kWhite);
  gStyle->SetTextSizePixels(26);
  gStyle->SetTextFont(42);
  //  gStyle->SetTickLength(0.04,"X");  gStyle->SetTickLength(0.04,"Y"); 

  gStyle->SetLegendBorderSize(0);
  gStyle->SetLegendFillColor(kWhite);
  //  gStyle->SetFillColor(kWhite);
  gStyle->SetLegendFont(42);
}




//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////// Context Utilities /////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

TString contextTrackDatasetComp(const char options[]){
  TString texcontextDatasetComp;
  //  if (strstr(options, "track") == NULL) {
  texcontextDatasetComp = *texDatasetsComparisonCommonDenominator;
  return texcontextDatasetComp;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////// QC  plot functions /////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void Draw_Pt_DatasetComparison() {

  TH1D* H1D_trackPt[nDatasets];
  TH1D* H1D_trackPt_rebinned[nDatasets];
  
  TH1D* H1D_trackPt_rebinned_ratios[nDatasets];

  bool divideSuccess = false;

  for(int iDataset = 0; iDataset < nDatasets; iDataset++){

    H1D_trackPt[iDataset] = (TH1D*)((TH1D*)file_O2Analysis_list[iDataset]->Get(analysisWorkflow+"/h_track_pt"))->Clone("Draw_Pt_DatasetComparison"+Datasets[iDataset]);

    H1D_trackPt_rebinned[iDataset] = (TH1D*)H1D_trackPt[iDataset]->Rebin(1.,"trackPt_rebinned_"+Datasets[iDataset]);

    // NormaliseYieldToNJets(H1D_trackPt_rebinned[iDataset]);
    NormaliseYieldToNEvents(H1D_trackPt_rebinned[iDataset], GetNEventsSel8(file_O2Analysis_list[iDataset]));

    H1D_trackPt_rebinned_ratios[iDataset] = (TH1D*)H1D_trackPt_rebinned[iDataset]->Clone("trackPt_rebinned_ratios"+Datasets[iDataset]);
    H1D_trackPt_rebinned_ratios[iDataset]->Reset("M");
    divideSuccess = H1D_trackPt_rebinned_ratios[iDataset]->Divide(H1D_trackPt_rebinned[iDataset], H1D_trackPt_rebinned[0]);
  }

  TString* pdfName = new TString("track_Pt_DataComp");
  TString* pdfName_ratio = new TString("track_Pt_DataComp_ratio");

  TString textContext(contextTrackDatasetComp(""));

  Draw_TH1_Histograms_in_one(H1D_trackPt_rebinned, DatasetsNames, nDatasets, textContext, pdfName, texPtX, texTrackNormPtYield, texCollisionDataInfo, "logy");
  if (divideSuccess == true) {
    Draw_TH1_Histograms_in_one(H1D_trackPt_rebinned_ratios, DatasetsNames, nDatasets, textContext, pdfName_ratio, texPtX, texRatioDatasets, texCollisionDataInfo, "autoratio,avoidFirst");
  }
  else {
    cout << "Divide failed in Draw_Pt_DatasetComparison" << endl;
  }
}


void Draw_Eta_DatasetComparison() {

  TH1D* H1D_trackEta[nDatasets];
  TH1D* H1D_trackEta_rebinned[nDatasets];
  
  TH1D* H1D_trackEta_rebinned_ratios[nDatasets];

  bool divideSuccess = false;

  for(int iDataset = 0; iDataset < nDatasets; iDataset++){

    H1D_trackEta[iDataset] = (TH1D*)((TH1D*)file_O2Analysis_list[iDataset]->Get(analysisWorkflow+"/h_track_eta"))->Clone("Draw_Eta_DatasetComparison"+Datasets[iDataset]);

    H1D_trackEta_rebinned[iDataset] = (TH1D*)H1D_trackEta[iDataset]->Rebin(1.,"trackEta_rebinned_"+Datasets[iDataset]);

    // NormaliseYieldToNJets(H1D_trackEta_rebinned[iDataset]);
    NormaliseYieldToNEvents(H1D_trackEta_rebinned[iDataset], GetNEventsSel8(file_O2Analysis_list[iDataset]));

    H1D_trackEta_rebinned_ratios[iDataset] = (TH1D*)H1D_trackEta_rebinned[iDataset]->Clone("trackEta_rebinned_ratios"+Datasets[iDataset]);
    H1D_trackEta_rebinned_ratios[iDataset]->Reset("M");
    divideSuccess = H1D_trackEta_rebinned_ratios[iDataset]->Divide(H1D_trackEta_rebinned[iDataset], H1D_trackEta_rebinned[0]);
  }

  TString* pdfName = new TString("track_Eta_DataComp");
  TString* pdfName_ratio = new TString("track_Eta_DataComp_ratio");

  TString textContext(contextTrackDatasetComp(""));

  Draw_TH1_Histograms_in_one(H1D_trackEta_rebinned, DatasetsNames, nDatasets, textContext, pdfName, texEtaX, texTrackNormEtaYield, texCollisionDataInfo, "");
  if (divideSuccess == true) {
    Draw_TH1_Histograms_in_one(H1D_trackEta_rebinned_ratios, DatasetsNames, nDatasets, textContext, pdfName_ratio, texEtaX, texRatioDatasets, texCollisionDataInfo, "autoratio,avoidFirst");
  }
  else {
    cout << "Divide failed in Draw_Eta_DatasetComparison" << endl;
  }
}

void Draw_Phi_DatasetComparison() {

  TH1D* H1D_trackPhi[nDatasets];
  TH1D* H1D_trackPhi_rebinned[nDatasets];
  
  TH1D* H1D_trackPhi_rebinned_ratios[nDatasets];

  bool divideSuccess = false;

  for(int iDataset = 0; iDataset < nDatasets; iDataset++){

    H1D_trackPhi[iDataset] = (TH1D*)((TH1D*)file_O2Analysis_list[iDataset]->Get(analysisWorkflow+"/h_track_phi"))->Clone("Draw_Phi_DatasetComparison"+Datasets[iDataset]);

    H1D_trackPhi_rebinned[iDataset] = (TH1D*)H1D_trackPhi[iDataset]->Rebin(1.,"trackPhi_rebinned_"+Datasets[iDataset]);

    // NormaliseYieldToNJets(H1D_trackPhi_rebinned[iDataset]);
    NormaliseYieldToNEvents(H1D_trackPhi_rebinned[iDataset], GetNEventsSel8(file_O2Analysis_list[iDataset]));

    H1D_trackPhi_rebinned_ratios[iDataset] = (TH1D*)H1D_trackPhi_rebinned[iDataset]->Clone("trackPhi_rebinned_ratios"+Datasets[iDataset]);
    H1D_trackPhi_rebinned_ratios[iDataset]->Reset("M");
    divideSuccess = H1D_trackPhi_rebinned_ratios[iDataset]->Divide(H1D_trackPhi_rebinned[iDataset], H1D_trackPhi_rebinned[0]);
  }

  TString* pdfName = new TString("track_Phi_DataComp");
  TString* pdfName_ratio = new TString("track_Phi_DataComp_ratio");

  TString textContext(contextTrackDatasetComp(""));

  Draw_TH1_Histograms_in_one(H1D_trackPhi_rebinned, DatasetsNames, nDatasets, textContext, pdfName, texPhiX, texTrackNormPhiYield, texCollisionDataInfo, "");
  if (divideSuccess == true) {
    Draw_TH1_Histograms_in_one(H1D_trackPhi_rebinned_ratios, DatasetsNames, nDatasets, textContext, pdfName_ratio, texPhiX, texRatioDatasets, texCollisionDataInfo, "autoratio,avoidFirst");
  }
  else {
    cout << "Divide failed in Draw_Phi_DatasetComparison" << endl;
  }
}
