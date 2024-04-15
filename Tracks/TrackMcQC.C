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
#include "./TrackMcQC_settings.h"
#include "../Settings/AxisTitles.h"
#include "../Settings/GlobalSettings.h"
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

// // Plot Utilities
// TString contextDatasetComp(const char options[]);

//////////// QC plot functions
// Dataset comparison

void Draw_Efficiency_Pt_DatasetComparison(float* etaRange);
void Draw_Efficiency_Eta_DatasetComparison(float* PtRange);
void Draw_Efficiency_Phi_DatasetComparison(float* ptRange, float* etaRange);
void Draw_Efficiency_Eta_PtRangeComparison(float* ptRange, int nPtRanges, int iDataset, const char options[]);
void Draw_Efficiency_Phi_PtRangeComparison(float* ptRange, int nPtRanges, float* etaRange, int iDataset, const char options[]);

/////////////////////////////////////////////////////
///////////////////// Main Macro ////////////////////
/////////////////////////////////////////////////////

void TrackMcQC() {
  // Load necessary libraries
  LoadLibs();
  // Set the default style
  SetStyle();

  // TString* SaveAs_Title = new TString("");
  TString* texXtitle = new TString("");
  TString* texYtitle = new TString("");
  // TString* Extra = new TString("");

  float etaRange[2] = {-0.8, 0.8};
  Draw_Efficiency_Pt_DatasetComparison(etaRange);
  float ptRange1[2] = {0.15, 100};
  Draw_Efficiency_Eta_DatasetComparison(ptRange1);
  Draw_Efficiency_Phi_DatasetComparison(ptRange1, etaRange);

  int nPtRanges = 4;
  float ptRange[5] = {0, 0.15, 0.5, 1, 100};
  for(int iDataset = 0; iDataset < nDatasets; iDataset++){
    Draw_Efficiency_Eta_PtRangeComparison(ptRange, nPtRanges, iDataset, "");
    Draw_Efficiency_Eta_PtRangeComparison(ptRange, nPtRanges, iDataset, "scaled");
    Draw_Efficiency_Phi_PtRangeComparison(ptRange, nPtRanges, etaRange, iDataset, "");
    Draw_Efficiency_Phi_PtRangeComparison(ptRange, nPtRanges, etaRange, iDataset, "scaled");
  }
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




// //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// ////////////////////////////////////////////////////////////////////////////// Context Utilities /////////////////////////////////////////////////////////////////////////////////////
// //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// TString contextDatasetComp(const char options[]){
//   TString texcontextDatasetCompAndRadiusAndVarRange;
//   //  if (strstr(options, "track") == NULL) {
//   texcontextDatasetCompAndRadiusAndVarRange = *texDatasetsComparisonCommonDenominator;
//   return texcontextDatasetCompAndRadiusAndVarRange;
// }

// TString contextPtRange(float* PtRange){
//   std::stringstream ss;
//   ss << PtRange[0] << " < #it{p}_{T} < " << PtRange[1];
//   TString textContext((TString)ss.str());
//   // TString texDataset(Form("%.0f", PtRange[0])+" < #it{p}_{T} < "+Form("%.0f", PtRange[1]));
//   return textContext;
// }

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////// QC  plot functions /////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


void Draw_Efficiency_Pt_DatasetComparison(float* etaRange) {

  TH3D* H3D_trackPtEtaPhi_mcparticles[nDatasets];
  TH3D* H3D_trackPtEtaPhi_assoctracks[nDatasets];
  TH1D* H1D_trackPt_mcparticles[nDatasets];
  TH1D* H1D_trackPt_assoctracks[nDatasets];
  
  TH1D* H1D_trackPt_efficiency[nDatasets];

  bool divideSuccess = false;

  for(int iDataset = 0; iDataset < nDatasets; iDataset++){

    H3D_trackPtEtaPhi_mcparticles[iDataset] = (TH3D*)((TH3D*)file_O2Analysis_list[iDataset]->Get(analysisWorkflow+"/h3_track_pt_track_eta_track_phi_mcparticles"))->Clone("Draw_Efficiency_Pt_DatasetComparison_mcparticles"+Datasets[iDataset]);
    H3D_trackPtEtaPhi_assoctracks[iDataset] = (TH3D*)((TH3D*)file_O2Analysis_list[iDataset]->Get(analysisWorkflow+"/h3_track_pt_track_eta_track_phi_associatedtrackSelColl"))->Clone("Draw_Efficiency_Pt_DatasetComparison_associatedtrackSelColl"+Datasets[iDataset]);

    int ibinEta_low_mcpart = H3D_trackPtEtaPhi_mcparticles[iDataset]->GetYaxis()->FindBin(etaRange[0]+deltaEtaMcVsTrackEfficiency);
    int ibinEta_high_mcpart = H3D_trackPtEtaPhi_mcparticles[iDataset]->GetYaxis()->FindBin(etaRange[1]-deltaEtaMcVsTrackEfficiency);
    int ibinEta_low_track = H3D_trackPtEtaPhi_mcparticles[iDataset]->GetYaxis()->FindBin(etaRange[0]);
    int ibinEta_high_track = H3D_trackPtEtaPhi_mcparticles[iDataset]->GetYaxis()->FindBin(etaRange[1]);

    H1D_trackPt_mcparticles[iDataset] = (TH1D*)H3D_trackPtEtaPhi_mcparticles[iDataset]->ProjectionX("trackPt_mcparticles"+Datasets[iDataset], ibinEta_low_mcpart, ibinEta_high_mcpart, 0, -1, "e");
    H1D_trackPt_assoctracks[iDataset] = (TH1D*)H3D_trackPtEtaPhi_assoctracks[iDataset]->ProjectionX("trackPt_assoctracks"+Datasets[iDataset], ibinEta_low_track, ibinEta_high_track, 0, -1, "e");
    

    H1D_trackPt_efficiency[iDataset] = (TH1D*)H1D_trackPt_mcparticles[iDataset]->Clone("trackPt_efficiency"+Datasets[iDataset]);
    H1D_trackPt_efficiency[iDataset]->Reset("M");
    cout << "tracks count = " << H1D_trackPt_assoctracks[iDataset]->GetEntries() << ", part count = " << H1D_trackPt_mcparticles[iDataset]->GetEntries() << endl;
    divideSuccess = H1D_trackPt_efficiency[iDataset]->Divide(H1D_trackPt_assoctracks[iDataset], H1D_trackPt_mcparticles[iDataset]);
  }

  TString* pdfNameEntriesNorm = new TString("track_Pt_efficiency"+dummyName[0]+"_@eta["+Form("%.1f", etaRange[0])+","+Form("%.1f", etaRange[1])+"]");

  TString textContext(contextCustomTwoFields(*texDatasetsComparisonCommonDenominator, contextEtaRange(etaRange), ""));

  if (divideSuccess == true) {
    Draw_TH1_Histograms_in_one(H1D_trackPt_efficiency, DatasetsNames, nDatasets, textContext, pdfNameEntriesNorm, texPtX, texTrackEfficiency, texCollisionDataInfo, "logx");
  }
  else {
    cout << "Divide failed in Draw_Efficiency_Eta_DatasetComparison" << endl;
  }
}


void Draw_Efficiency_Eta_DatasetComparison(float* ptRange) {

  TH3D* H3D_trackPtEtaPhi_mcparticles[nDatasets];
  TH3D* H3D_trackPtEtaPhi_assoctracks[nDatasets];
  TH1D* H1D_trackEta_mcparticles[nDatasets];
  TH1D* H1D_trackEta_assoctracks[nDatasets];
  
  TH1D* H1D_trackEta_efficiency[nDatasets];

  bool divideSuccess = false;

  for(int iDataset = 0; iDataset < nDatasets; iDataset++){

    H3D_trackPtEtaPhi_mcparticles[iDataset] = (TH3D*)((TH3D*)file_O2Analysis_list[iDataset]->Get(analysisWorkflow+"/h3_track_pt_track_eta_track_phi_mcparticles"))->Clone("Draw_Efficiency_Eta_DatasetComparison_mcparticles"+Datasets[iDataset]);
    H3D_trackPtEtaPhi_assoctracks[iDataset] = (TH3D*)((TH3D*)file_O2Analysis_list[iDataset]->Get(analysisWorkflow+"/h3_track_pt_track_eta_track_phi_associatedtrackSelColl"))->Clone("Draw_Efficiency_Eta_DatasetComparison_associatedtrackSelColl"+Datasets[iDataset]);
    // int ibincutlow = H1D_trackEta[iDataset]->GetXaxis()->FindBin(-0.8);
    // int ibincuthigh = H1D_trackEta[iDataset]->GetXaxis()->FindBin(0.8);
    // for (int ibin = 1; ibin < ibincutlow; ibin++){
    //   H1D_trackEta[iDataset]->SetBinContent(ibin, 0);
    // }
    // for (int ibin = ibincuthigh; ibin < H1D_trackEta[iDataset]->GetNbinsX(); ibin++){
    //   H1D_trackEta[iDataset]->SetBinContent(ibin, 0);
    // }


    int ibinPt_low = H3D_trackPtEtaPhi_mcparticles[iDataset]->GetXaxis()->FindBin(ptRange[0]);
    int ibinPt_high = H3D_trackPtEtaPhi_mcparticles[iDataset]->GetXaxis()->FindBin(ptRange[1]);
    H3D_trackPtEtaPhi_mcparticles[iDataset]->GetXaxis()->SetRange(ibinPt_low,ibinPt_high);
    H3D_trackPtEtaPhi_assoctracks[iDataset]->GetXaxis()->SetRange(ibinPt_low,ibinPt_high);

    H1D_trackEta_mcparticles[iDataset] = (TH1D*)H3D_trackPtEtaPhi_mcparticles[iDataset]->ProjectionY("trackEta_mcparticles"+Datasets[iDataset]+"_@pt["+Form("%.2f", ptRange[0])+","+Form("%.1f", ptRange[1])+"]", ibinPt_low, ibinPt_high, 0, -1, "e");
    H1D_trackEta_assoctracks[iDataset] = (TH1D*)H3D_trackPtEtaPhi_assoctracks[iDataset]->ProjectionY("trackEta_assoctracks"+Datasets[iDataset]+"_@pt["+Form("%.2f", ptRange[0])+","+Form("%.1f", ptRange[1])+"]", ibinPt_low, ibinPt_high, 0, -1, "e");
    // H1D_trackEta_rebinned[iDataset] = (TH1D*)H1D_trackEta[iDataset]->Rebin(1.,"trackEta_rebinned"+Datasets[iDataset]);

    // NormaliseYieldToNEntries(H1D_trackEta_rebinned[iDataset]);
    // H1D_trackEta_rebinned[iDataset]->Scale(1./H1D_trackEta_rebinned[iDataset]->Integral(1, H1D_trackEta_rebinned[iDataset]->GetNbinsX()),"width"); // If option contains "width" the bin contents and errors are divided by the bin width.
    // H1D_trackEta_rebinned[iDataset]->Scale(1./H1D_trackEta_rebinned[iDataset]->Integral(ibincutlow, ibincuthigh),"width"); // If option contains "width" the bin contents and errors are divided by the bin width.

    H1D_trackEta_efficiency[iDataset] = (TH1D*)H1D_trackEta_mcparticles[iDataset]->Clone("trackEta_efficiency"+Datasets[iDataset]);
    H1D_trackEta_efficiency[iDataset]->Reset("M");
    divideSuccess = H1D_trackEta_efficiency[iDataset]->Divide(H1D_trackEta_assoctracks[iDataset], H1D_trackEta_mcparticles[iDataset]);
  }

  TString* pdfNameEntriesNorm = new TString("track_Eta_efficiency"+dummyName[0]+"_@pt["+Form("%.2f", ptRange[0])+","+Form("%.1f", ptRange[1])+"]");

  // TString textContext(contextDatasetComp(""));
  TString textContext(contextCustomTwoFields(*texDatasetsComparisonCommonDenominator, contextPtRange(ptRange), ""));

  if (divideSuccess == true) {
    Draw_TH1_Histograms_in_one(H1D_trackEta_efficiency, DatasetsNames, nDatasets, textContext, pdfNameEntriesNorm, texEtaX, texTrackEfficiency, texCollisionDataInfo, "");
  }
  else {
    cout << "Divide failed in Draw_Efficiency_Eta_DatasetComparison" << endl;
  }
}

void Draw_Efficiency_Phi_DatasetComparison(float* ptRange, float* etaRange) {

  TH3D* H3D_trackPtEtaPhi_mcparticles[nDatasets];
  TH3D* H3D_trackPtEtaPhi_assoctracks[nDatasets];
  TH1D* H1D_trackPhi_mcparticles[nDatasets];
  TH1D* H1D_trackPhi_assoctracks[nDatasets];
  
  TH1D* H1D_trackPhi_efficiency[nDatasets];

  bool divideSuccess = false;

  for(int iDataset = 0; iDataset < nDatasets; iDataset++){

    H3D_trackPtEtaPhi_mcparticles[iDataset] = (TH3D*)((TH3D*)file_O2Analysis_list[iDataset]->Get(analysisWorkflow+"/h3_track_pt_track_eta_track_phi_mcparticles"))->Clone("Draw_Efficiency_Phi_DatasetComparison_mcparticles"+Datasets[iDataset]);
    H3D_trackPtEtaPhi_assoctracks[iDataset] = (TH3D*)((TH3D*)file_O2Analysis_list[iDataset]->Get(analysisWorkflow+"/h3_track_pt_track_eta_track_phi_associatedtrackSelColl"))->Clone("Draw_Efficiency_Phi_DatasetComparison_associatedtrackSelColl"+Datasets[iDataset]);
    // int ibincutlow = H1D_trackPhi[iDataset]->GetXaxis()->FindBin(-0.8);
    // int ibincuthigh = H1D_trackPhi[iDataset]->GetXaxis()->FindBin(0.8);
    // for (int ibin = 1; ibin < ibincutlow; ibin++){
    //   H1D_trackPhi[iDataset]->SetBinContent(ibin, 0);
    // }
    // for (int ibin = ibincuthigh; ibin < H1D_trackPhi[iDataset]->GetNbinsX(); ibin++){
    //   H1D_trackPhi[iDataset]->SetBinContent(ibin, 0);
    // }


    int ibinPt_low = H3D_trackPtEtaPhi_mcparticles[iDataset]->GetXaxis()->FindBin(ptRange[0]);
    int ibinPt_high = H3D_trackPtEtaPhi_mcparticles[iDataset]->GetXaxis()->FindBin(ptRange[1]);
    H3D_trackPtEtaPhi_mcparticles[iDataset]->GetXaxis()->SetRange(ibinPt_low,ibinPt_high);
    H3D_trackPtEtaPhi_assoctracks[iDataset]->GetXaxis()->SetRange(ibinPt_low,ibinPt_high);
    int ibinEta_low_mcpart = H3D_trackPtEtaPhi_mcparticles[iDataset]->GetYaxis()->FindBin(etaRange[0]+deltaEtaMcVsTrackEfficiency);
    int ibinEta_high_mcpart = H3D_trackPtEtaPhi_mcparticles[iDataset]->GetYaxis()->FindBin(etaRange[1]-deltaEtaMcVsTrackEfficiency);
    int ibinEta_low_track = H3D_trackPtEtaPhi_mcparticles[iDataset]->GetYaxis()->FindBin(etaRange[0]);
    int ibinEta_high_track = H3D_trackPtEtaPhi_mcparticles[iDataset]->GetYaxis()->FindBin(etaRange[1]);
    H3D_trackPtEtaPhi_mcparticles[iDataset]->GetXaxis()->SetRange(ibinEta_low_mcpart,ibinEta_high_mcpart);
    H3D_trackPtEtaPhi_assoctracks[iDataset]->GetXaxis()->SetRange(ibinEta_low_track,ibinEta_high_track);

    H1D_trackPhi_mcparticles[iDataset] = (TH1D*)H3D_trackPtEtaPhi_mcparticles[iDataset]->ProjectionZ("trackPhi_mcparticles"+Datasets[iDataset]+"_@pt["+Form("%.2f", ptRange[0])+","+Form("%.1f", ptRange[1])+"]", ibinPt_low, ibinPt_high, ibinEta_low_mcpart, ibinEta_high_mcpart, "e");
    H1D_trackPhi_assoctracks[iDataset] = (TH1D*)H3D_trackPtEtaPhi_assoctracks[iDataset]->ProjectionZ("trackPhi_assoctracks"+Datasets[iDataset]+"_@pt["+Form("%.2f", ptRange[0])+","+Form("%.1f", ptRange[1])+"]", ibinPt_low, ibinPt_high, ibinEta_low_track, ibinEta_high_track, "e");
    // H1D_trackPhi_rebinned[iDataset] = (TH1D*)H1D_trackPhi[iDataset]->Rebin(1.,"trackPhi_rebinned"+Datasets[iDataset]);

    // NormaliseYieldToNEntries(H1D_trackPhi_rebinned[iDataset]);
    // H1D_trackPhi_rebinned[iDataset]->Scale(1./H1D_trackPhi_rebinned[iDataset]->Integral(1, H1D_trackPhi_rebinned[iDataset]->GetNbinsX()),"width"); // If option contains "width" the bin contents and errors are divided by the bin width.
    // H1D_trackPhi_rebinned[iDataset]->Scale(1./H1D_trackPhi_rebinned[iDataset]->Integral(ibincutlow, ibincuthigh),"width"); // If option contains "width" the bin contents and errors are divided by the bin width.

    H1D_trackPhi_efficiency[iDataset] = (TH1D*)H1D_trackPhi_mcparticles[iDataset]->Clone("trackPhi_efficiency"+Datasets[iDataset]);
    H1D_trackPhi_efficiency[iDataset]->Reset("M");
    divideSuccess = H1D_trackPhi_efficiency[iDataset]->Divide(H1D_trackPhi_assoctracks[iDataset], H1D_trackPhi_mcparticles[iDataset]);
  }

  TString* pdfNameEntriesNorm = new TString("track_Phi_efficiency"+dummyName[0]+"_@pt["+Form("%.2f", ptRange[0])+","+Form("%.1f", ptRange[1])+"]"+"_@eta["+Form("%.1f", etaRange[0])+","+Form("%.1f", etaRange[1])+"]");

  // TString textContext(contextDatasetComp(""));
  TString textContext(contextCustomTwoFields(*texDatasetsComparisonCommonDenominator, "#splitline{"+contextPtRange(ptRange)+"}{"+contextEtaRange(etaRange)+"}", ""));

  if (divideSuccess == true) {
    Draw_TH1_Histograms_in_one(H1D_trackPhi_efficiency, DatasetsNames, nDatasets, textContext, pdfNameEntriesNorm, texPhiX, texTrackEfficiency, texCollisionDataInfo, "");
  }
  else {
    cout << "Divide failed in Draw_Efficiency_Phi_DatasetComparison" << endl;
  }
}

void Draw_Efficiency_Eta_PtRangeComparison(float* ptRange, int nPtRanges, int iDataset, const char options[]) {

  TH3D* H3D_trackPtEtaPhi_mcparticles;
  TH3D* H3D_trackPtEtaPhi_assoctracks;

  H3D_trackPtEtaPhi_mcparticles = (TH3D*)((TH3D*)file_O2Analysis_list[iDataset]->Get(analysisWorkflow+"/h3_track_pt_track_eta_track_phi_mcparticles"))->Clone("Draw_Efficiency_Eta_PtRangeComparison_mcparticles"+Datasets[iDataset]);
  H3D_trackPtEtaPhi_assoctracks = (TH3D*)((TH3D*)file_O2Analysis_list[iDataset]->Get(analysisWorkflow+"/h3_track_pt_track_eta_track_phi_associatedtrackSelColl"))->Clone("Draw_Efficiency_Eta_PtRangeComparison_associatedtrackSelColl"+Datasets[iDataset]);
  // H3D_jetRjetPtjetEta = (TH3D*)((TH3D*)file_O2Analysis_list[iDataset]->Get(analysisWorkflow+"/h3_jet_r_jet_pt_jet_eta"+jetFinderQaHistType[iJetFinderQaType]))->Clone("Draw_Eta_PtCutComparison"+Datasets[iDataset]);

  TH1D* H1D_trackEta_mcparticles[nPtRanges];
  TH1D* H1D_trackEta_assoctracks[nPtRanges];
  
  TH1D* H1D_trackEta_efficiency[nPtRanges];

  bool divideSuccess = false;


  TString PtCutsLegend[nPtRanges];
  std::stringstream ss;

  TString* yAxisLabel;

  for(int iBinPt = 0; iBinPt < nPtRanges; iBinPt++){
    // int ibincutlow = H1D_trackEta[iDataset]->GetXaxis()->FindBin(-0.8);
    // int ibincuthigh = H1D_trackEta[iDataset]->GetXaxis()->FindBin(0.8);
    // for (int ibin = 1; ibin < ibincutlow; ibin++){
    //   H1D_trackEta[iDataset]->SetBinContent(ibin, 0);
    // }
    // for (int ibin = ibincuthigh; ibin < H1D_trackEta[iDataset]->GetNbinsX(); ibin++){
    //   H1D_trackEta[iDataset]->SetBinContent(ibin, 0);
    // }


    int ibinPt_low = H3D_trackPtEtaPhi_mcparticles->GetXaxis()->FindBin(ptRange[iBinPt]);
    int ibinPt_high = H3D_trackPtEtaPhi_mcparticles->GetXaxis()->FindBin(ptRange[iBinPt+1]);
    H3D_trackPtEtaPhi_mcparticles->GetXaxis()->SetRange(ibinPt_low,ibinPt_high);
    H3D_trackPtEtaPhi_assoctracks->GetXaxis()->SetRange(ibinPt_low,ibinPt_high);

    H1D_trackEta_mcparticles[iBinPt] = (TH1D*)H3D_trackPtEtaPhi_mcparticles->ProjectionY("trackEta_mcparticles"+Datasets[iDataset]+"_@pt["+Form("%.2f", ptRange[iBinPt])+","+Form("%.1f", ptRange[iBinPt+1])+"]", ibinPt_low, ibinPt_high, 0, -1, "e");
    H1D_trackEta_assoctracks[iBinPt] = (TH1D*)H3D_trackPtEtaPhi_assoctracks->ProjectionY("trackEta_assoctracks"+Datasets[iDataset]+"_@pt["+Form("%.2f", ptRange[iBinPt])+","+Form("%.1f", ptRange[iBinPt+1])+"]", ibinPt_low, ibinPt_high, 0, -1, "e");
    // H1D_trackEta_rebinned[iDataset] = (TH1D*)H1D_trackEta[iDataset]->Rebin(1.,"trackEta_rebinned"+Datasets[iDataset]);

    // NormaliseYieldToNEntries(H1D_trackEta_rebinned[iDataset]);
    // H1D_trackEta_rebinned[iDataset]->Scale(1./H1D_trackEta_rebinned[iDataset]->Integral(1, H1D_trackEta_rebinned[iDataset]->GetNbinsX()),"width"); // If option contains "width" the bin contents and errors are divided by the bin width.
    // H1D_trackEta_rebinned[iDataset]->Scale(1./H1D_trackEta_rebinned[iDataset]->Integral(ibincutlow, ibincuthigh),"width"); // If option contains "width" the bin contents and errors are divided by the bin width.

    H1D_trackEta_efficiency[iBinPt] = (TH1D*)H1D_trackEta_mcparticles[iBinPt]->Clone("trackEta_efficiency"+Datasets[iBinPt]);
    H1D_trackEta_efficiency[iBinPt]->Reset("M");
    divideSuccess = H1D_trackEta_efficiency[iBinPt]->Divide(H1D_trackEta_assoctracks[iBinPt], H1D_trackEta_mcparticles[iBinPt]);
    
    if (strstr(options, "scaled") != NULL) {
        NormaliseYieldToIntegral(H1D_trackEta_efficiency[iBinPt]);
    }

    ss << "pt[" << ptRange[iBinPt] << ", " << ptRange[iBinPt+1] << "]";
    PtCutsLegend[iBinPt] = (TString)ss.str();
    ss.str("");
    ss.clear();
  }
  // for(int iBinPt = 0; iBinPt < nPtRanges; iBinPt++){
  //   H1D_trackEta_efficiency[iBinPt]->Scale(H1D_trackEta_efficiency[nPtRanges-1]->GetEntries()/H1D_trackEta_efficiency[iBinPt]->GetEntries(),"width"); // If option contains "width" the bin contents and errors are divided by the bin width.
  // }

  TString optDelimiter;
  strcmp(options, "scaled") == 0 ? optDelimiter = "_" : optDelimiter = "";

  TString* pdfNameEntriesNorm = new TString("track_Eta_efficiency_ptRangeComp"+Datasets[iDataset]+optDelimiter+options);

  // TString textContext(contextDatasetComp(""));
  TString textContext(contextCustomTwoFields(*texDatasetsComparisonCommonDenominator, DatasetsNames[iDataset],  ""));

  if (divideSuccess == true) {
    Draw_TH1_Histograms_in_one(H1D_trackEta_efficiency, PtCutsLegend, nPtRanges, textContext, pdfNameEntriesNorm, texEtaX, texTrackEfficiency, texCollisionDataInfo, "");
  }
  else {
    cout << "Divide failed in Draw_Efficiency_Eta_DatasetComparison" << endl;
  }
}


void Draw_Efficiency_Phi_PtRangeComparison(float* ptRange, int nPtRanges, float* etaRange, int iDataset, const char options[]) {

  TH3D* H3D_trackPtEtaPhi_mcparticles;
  TH3D* H3D_trackPtEtaPhi_assoctracks;

  H3D_trackPtEtaPhi_mcparticles = (TH3D*)((TH3D*)file_O2Analysis_list[iDataset]->Get(analysisWorkflow+"/h3_track_pt_track_eta_track_phi_mcparticles"))->Clone("Draw_Efficiency_Phi_PtRangeComparison_mcparticles"+Datasets[iDataset]);
  H3D_trackPtEtaPhi_assoctracks = (TH3D*)((TH3D*)file_O2Analysis_list[iDataset]->Get(analysisWorkflow+"/h3_track_pt_track_eta_track_phi_associatedtrackSelColl"))->Clone("Draw_Efficiency_Phi_PtRangeComparison_associatedtrackSelColl"+Datasets[iDataset]);
  // H3D_jetRjetPtjetEta = (TH3D*)((TH3D*)file_O2Analysis_list[iDataset]->Get(analysisWorkflow+"/h3_jet_r_jet_pt_jet_eta"+jetFinderQaHistType[iJetFinderQaType]))->Clone("Draw_Phi_PtCutComparison"+Datasets[iDataset]);

  TH1D* H1D_trackPhi_mcparticles[nPtRanges];
  TH1D* H1D_trackPhi_assoctracks[nPtRanges];
  
  TH1D* H1D_trackPhi_efficiency[nPtRanges];

  bool divideSuccess = false;


  TString PtCutsLegend[nPtRanges];
  std::stringstream ss;

  TString* yAxisLabel;

  for(int iBinPt = 0; iBinPt < nPtRanges; iBinPt++){
    // int ibincutlow = H1D_trackPhi[iDataset]->GetXaxis()->FindBin(-0.8);
    // int ibincuthigh = H1D_trackPhi[iDataset]->GetXaxis()->FindBin(0.8);
    // for (int ibin = 1; ibin < ibincutlow; ibin++){
    //   H1D_trackPhi[iDataset]->SetBinContent(ibin, 0);
    // }
    // for (int ibin = ibincuthigh; ibin < H1D_trackPhi[iDataset]->GetNbinsX(); ibin++){
    //   H1D_trackPhi[iDataset]->SetBinContent(ibin, 0);
    // }

    int ibinPt_low = H3D_trackPtEtaPhi_mcparticles->GetXaxis()->FindBin(ptRange[iBinPt]);
    int ibinPt_high = H3D_trackPtEtaPhi_mcparticles->GetXaxis()->FindBin(ptRange[iBinPt+1]);
    H3D_trackPtEtaPhi_mcparticles->GetXaxis()->SetRange(ibinPt_low,ibinPt_high);
    H3D_trackPtEtaPhi_assoctracks->GetXaxis()->SetRange(ibinPt_low,ibinPt_high);
    int ibinEta_low_mcpart = H3D_trackPtEtaPhi_mcparticles->GetYaxis()->FindBin(etaRange[0]+deltaEtaMcVsTrackEfficiency);
    int ibinEta_high_mcpart = H3D_trackPtEtaPhi_mcparticles->GetYaxis()->FindBin(etaRange[1]-deltaEtaMcVsTrackEfficiency);
    int ibinEta_low_track = H3D_trackPtEtaPhi_mcparticles->GetYaxis()->FindBin(etaRange[0]);
    int ibinEta_high_track = H3D_trackPtEtaPhi_mcparticles->GetYaxis()->FindBin(etaRange[1]);
    H3D_trackPtEtaPhi_mcparticles->GetXaxis()->SetRange(ibinEta_low_mcpart,ibinEta_high_mcpart);
    H3D_trackPtEtaPhi_assoctracks->GetXaxis()->SetRange(ibinEta_low_track,ibinEta_high_track);

    H1D_trackPhi_mcparticles[iBinPt] = (TH1D*)H3D_trackPtEtaPhi_mcparticles->ProjectionZ("trackPhi_mcparticles"+Datasets[iDataset]+"_@pt["+Form("%.2f", ptRange[iBinPt])+","+Form("%.1f", ptRange[iBinPt+1])+"]", ibinPt_low, ibinPt_high, ibinEta_low_mcpart, ibinEta_high_mcpart, "e");
    H1D_trackPhi_assoctracks[iBinPt] = (TH1D*)H3D_trackPtEtaPhi_assoctracks->ProjectionZ("trackPhi_assoctracks"+Datasets[iDataset]+"_@pt["+Form("%.2f", ptRange[iBinPt])+","+Form("%.1f", ptRange[iBinPt+1])+"]", ibinPt_low, ibinPt_high, ibinEta_low_track, ibinEta_high_track, "e");
    // H1D_trackPhi_rebinned[iDataset] = (TH1D*)H1D_trackPhi[iDataset]->Rebin(1.,"trackPhi_rebinned"+Datasets[iDataset]);

    // NormaliseYieldToNEntries(H1D_trackPhi_rebinned[iDataset]);
    // H1D_trackPhi_rebinned[iDataset]->Scale(1./H1D_trackPhi_rebinned[iDataset]->Integral(1, H1D_trackPhi_rebinned[iDataset]->GetNbinsX()),"width"); // If option contains "width" the bin contents and errors are divided by the bin width.
    // H1D_trackPhi_rebinned[iDataset]->Scale(1./H1D_trackPhi_rebinned[iDataset]->Integral(ibincutlow, ibincuthigh),"width"); // If option contains "width" the bin contents and errors are divided by the bin width.

    H1D_trackPhi_efficiency[iBinPt] = (TH1D*)H1D_trackPhi_mcparticles[iBinPt]->Clone("trackPhi_efficiency"+Datasets[iBinPt]);
    H1D_trackPhi_efficiency[iBinPt]->Reset("M");
    divideSuccess = H1D_trackPhi_efficiency[iBinPt]->Divide(H1D_trackPhi_assoctracks[iBinPt], H1D_trackPhi_mcparticles[iBinPt]);
    
    if (strstr(options, "scaled") != NULL) {
        NormaliseYieldToIntegral(H1D_trackPhi_efficiency[iBinPt]);
    }

    ss << "pt[" << ptRange[iBinPt] << ", " << ptRange[iBinPt+1] << "]";
    PtCutsLegend[iBinPt] = (TString)ss.str();
    ss.str("");
    ss.clear();
  }
//   for(int iBinPt = 0; iBinPt < nPtRanges; iBinPt++){
//     H1D_trackPhi_efficiency[iBinPt]->Scale(H1D_trackPhi_efficiency[nPtRanges-1]->GetEntries()/H1D_trackPhi_efficiency[iBinPt]->GetEntries(),"width"); // If option contains "width" the bin contents and errors are divided by the bin width.
//   }

  TString optDelimiter;
  strcmp(options, "scaled") == 0 ? optDelimiter = "_" : optDelimiter = "";

  TString* pdfNameEntriesNorm = new TString("track_Phi_efficiency_ptRangeComp"+Datasets[iDataset]+"_@eta["+Form("%.1f", etaRange[0])+","+Form("%.1f", etaRange[1])+"]"+optDelimiter+options);

  // TString textContext(contextDatasetComp(""));
  TString textContext(contextCustomThreeFields(*texDatasetsComparisonCommonDenominator, DatasetsNames[iDataset], contextEtaRange(etaRange), ""));

  if (divideSuccess == true) {
    Draw_TH1_Histograms_in_one(H1D_trackPhi_efficiency, PtCutsLegend, nPtRanges, textContext, pdfNameEntriesNorm, texPhiX, texTrackEfficiency, texCollisionDataInfo, "");
  }
  else {
    cout << "Divide failed in Draw_Efficiency_Phi_DatasetComparison" << endl;
  }
}


