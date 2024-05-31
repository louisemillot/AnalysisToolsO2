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
void Draw_Efficiency_Phi_DatasetComparison_finerPhi(float* ptRange, float* etaRange);

void Draw_Efficiency_Pt_DatasetComparison_legacy(float* etaRange);
void Draw_Efficiency_Eta_DatasetComparison_legacy(float* PtRange);
void Draw_Efficiency_Phi_DatasetComparison_legacy(float* ptRange, float* etaRange);
void Draw_Efficiency_Eta_PtRangeComparison_legacy(float* ptRange, int nPtRanges, int iDataset, const char options[]);
void Draw_Efficiency_Phi_PtRangeComparison_legacy(float* ptRange, int nPtRanges, float* etaRange, int iDataset, const char options[]);
void Draw_Efficiency_Phi_DatasetComparison_finerPhi_legacy(float* ptRange, float* etaRange);

void Draw_mcPt_DatasetComparison(int particleStatusOption);
void Draw_Vx_DatasetComparison(int particleStatusOption);
void Draw_Vy_DatasetComparison(int particleStatusOption);
void Draw_Vz_DatasetComparison(int particleStatusOption);
void Draw_y_iu_DatasetComparison();
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

  // float etaRange[2] = {-0.9, 0.9};
  // Draw_Efficiency_Pt_DatasetComparison(etaRange);
  // float ptRange1[2] = {0.15, 100};
  // Draw_Efficiency_Eta_DatasetComparison(ptRange1);
  // Draw_Efficiency_Phi_DatasetComparison(ptRange1, etaRange);
  // // Draw_Efficiency_Phi_DatasetComparison_finerPhi(ptRange1, etaRange); // only works with very specific datasets created locally

  // int nPtRanges = 4;
  // float ptRange[5] = {0, 0.15, 0.5, 1, 100};
  // for(int iDataset = 0; iDataset < nDatasets; iDataset++){
  //   Draw_Efficiency_Eta_PtRangeComparison(ptRange, nPtRanges, iDataset, "");
  //   Draw_Efficiency_Eta_PtRangeComparison(ptRange, nPtRanges, iDataset, "scaled");
  //   Draw_Efficiency_Phi_PtRangeComparison(ptRange, nPtRanges, etaRange, iDataset, "");
  //   Draw_Efficiency_Phi_PtRangeComparison(ptRange, nPtRanges, etaRange, iDataset, "scaled");
  // }
  int particleStatusOption = 2;
  Draw_mcPt_DatasetComparison(particleStatusOption);
  Draw_Vx_DatasetComparison(particleStatusOption);
  Draw_Vy_DatasetComparison(particleStatusOption);
  Draw_Vz_DatasetComparison(particleStatusOption);
  Draw_y_iu_DatasetComparison();

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


void Draw_Efficiency_Pt_DatasetComparison_legacy(float* etaRange) {
  cout << "test0" << endl;
  TH3D* H3D_trackPtEtaPhi_mcparticles[nDatasets];
  TH3D* H3D_trackPtEtaPhi_assoctracks[nDatasets];
  TH3D* H3D_trackPtEtaPhi_assoctracks_split[nDatasets];
  TH3D* H3D_trackPtEtaPhi_assoctracks_nonSplitSecondary[nDatasets];

  TH1D* H1D_trackPt_mcparticles[nDatasets];

  TH1D* H1D_trackPt_assoctracks[nDatasets];
  TH1D* H1D_trackPt_assoctracks_split[nDatasets];
  TH1D* H1D_trackPt_assoctracks_nonSplitSecondary[nDatasets];
  
  TH1D* H1D_trackPt_efficiency[nDatasets];
  // TH1D* H1D_trackPt_efficiency_split[nDatasets];
  TH1D* H1D_trackPt_efficiency_splitCorrected[nDatasets];
  TH1D* H1D_trackPt_efficiency_splitAndSecondaryCorrected[nDatasets];

  bool divideSuccess = false;
  bool divideSuccess_split = false;

  cout << "test1" << endl;

  for(int iDataset = 0; iDataset < nDatasets; iDataset++){

    H3D_trackPtEtaPhi_mcparticles[iDataset] = (TH3D*)((TH3D*)file_O2Analysis_list[iDataset]->Get(analysisWorkflow[iDataset]+"/h3_track_pt_track_eta_track_phi_mcparticles"))->Clone("Draw_Efficiency_Pt_DatasetComparison_mcparticles"+Datasets[iDataset]);
    cout << "test2.1" << endl;
    H3D_trackPtEtaPhi_assoctracks[iDataset] = (TH3D*)((TH3D*)file_O2Analysis_list[iDataset]->Get(analysisWorkflow[iDataset]+"/h3_track_pt_track_eta_track_phi_associatedtrackSelColl"))->Clone("Draw_Efficiency_Pt_DatasetComparison_associatedtrackSelColl"+Datasets[iDataset]);
    cout << "test2.2" << endl;
    H3D_trackPtEtaPhi_assoctracks_split[iDataset] = (TH3D*)((TH3D*)file_O2Analysis_list[iDataset]->Get(analysisWorkflow[iDataset]+"/h3_track_pt_track_eta_track_phi_associatedtrackSelCollSplit"))->Clone("Draw_Efficiency_Pt_DatasetComparison_associatedtrackSelCollSplit"+Datasets[iDataset]);
    // H3D_trackPtEtaPhi_assoctracks_nonSplitSecondary[iDataset] = (TH3D*)((TH3D*)file_O2Analysis_list[iDataset]->Get(analysisWorkflow[iDataset]+"/h3_track_pt_track_eta_track_phi_associatedtrackSelCollNonSplitNonPrimary"))->Clone("Draw_Efficiency_Pt_DatasetComparison_associatedtrackSelCollNonSplitSecondary"+Datasets[iDataset]);
    cout << "test2.3" << endl;

    int ibinEta_low_mcpart = H3D_trackPtEtaPhi_mcparticles[iDataset]->GetYaxis()->FindBin(etaRange[0]+deltaEtaMcVsTrackEfficiency);
    int ibinEta_high_mcpart = H3D_trackPtEtaPhi_mcparticles[iDataset]->GetYaxis()->FindBin(etaRange[1]-deltaEtaMcVsTrackEfficiency);
    int ibinEta_low_track = H3D_trackPtEtaPhi_mcparticles[iDataset]->GetYaxis()->FindBin(etaRange[0]);
    int ibinEta_high_track = H3D_trackPtEtaPhi_mcparticles[iDataset]->GetYaxis()->FindBin(etaRange[1]);

    // cout << "H3D_trackPtEtaPhi_assoctracks -  preSetRange bin count = " << H3D_trackPtEtaPhi_assoctracks[iDataset]->GetNbinsY() << endl;
    // H3D_trackPtEtaPhi_mcparticles[iDataset]->GetYaxis()->SetRange(ibinEta_low_mcpart,ibinEta_high_mcpart);
    // H3D_trackPtEtaPhi_assoctracks[iDataset]->GetYaxis()->SetRange(ibinEta_low_track,ibinEta_high_track);
    // H3D_trackPtEtaPhi_assoctracks_split[iDataset]->GetYaxis()->SetRange(ibinEta_low_track,ibinEta_high_track);
    // H3D_trackPtEtaPhi_assoctracks_nonSplitSecondary[iDataset]->GetYaxis()->SetRange(ibinEta_low_track,ibinEta_high_track);
    // cout << "H3D_trackPtEtaPhi_assoctracks - postSetRange bin count = " << H3D_trackPtEtaPhi_assoctracks[iDataset]->GetNbinsY() << endl;

    H1D_trackPt_mcparticles[iDataset] = (TH1D*)H3D_trackPtEtaPhi_mcparticles[iDataset]->ProjectionX("trackPt_mcparticles"+Datasets[iDataset], ibinEta_low_mcpart, ibinEta_high_mcpart, 0, -1, "e");
    cout << "test2.31" << endl;
    H1D_trackPt_assoctracks[iDataset] = (TH1D*)H3D_trackPtEtaPhi_assoctracks[iDataset]->ProjectionX("trackPt_tracks_split"+Datasets[iDataset], ibinEta_low_track, ibinEta_high_track, 0, -1, "e");
    cout << "test2.32" << endl;
    H1D_trackPt_assoctracks_split[iDataset] = (TH1D*)H3D_trackPtEtaPhi_assoctracks_split[iDataset]->ProjectionX("trackPt_assoctracks_split"+Datasets[iDataset], ibinEta_low_track, ibinEta_high_track, 0, -1, "e");
    cout << "test2.33" << endl;
    // H1D_trackPt_assoctracks_nonSplitSecondary[iDataset] = (TH1D*)H3D_trackPtEtaPhi_assoctracks_nonSplitSecondary[iDataset]->ProjectionX("trackPt_assoctracks_nonSplitSecondary"+Datasets[iDataset], ibinEta_low_track, ibinEta_high_track, 0, -1, "e");
    // cout << "test2.4" << endl;


    // naming the histogram and resetting it to have a chosen name
    H1D_trackPt_efficiency[iDataset] = (TH1D*)H1D_trackPt_mcparticles[iDataset]->Clone("trackPt_efficiency"+Datasets[iDataset]);
    H1D_trackPt_efficiency[iDataset]->Reset("M");
    cout << "  - " << DatasetsNames[iDataset] << ":" << endl;
    cout << "      tracks count = " << H1D_trackPt_assoctracks[iDataset]->GetEntries() << ", part count = " << H1D_trackPt_mcparticles[iDataset]->GetEntries() << endl;
    divideSuccess = H1D_trackPt_efficiency[iDataset]->Divide(H1D_trackPt_assoctracks[iDataset], H1D_trackPt_mcparticles[iDataset]);
    cout << "test2.5" << endl;

    // H1D_trackPt_efficiency_split[iDataset] = (TH1D*)H1D_trackPt_mcparticles[iDataset]->Clone("trackPt_efficiency_split"+Datasets[iDataset]);
    // H1D_trackPt_efficiency_split[iDataset]->Reset("M");
    // cout << "split tracks count = " << H1D_trackPt_assoctracks_split[iDataset]->GetEntries() << ", part count = " << H1D_trackPt_mcparticles[iDataset]->GetEntries() << endl;
    // divideSuccess_split = H1D_trackPt_efficiency_split[iDataset]->Divide(H1D_trackPt_assoctracks_split[iDataset], H1D_trackPt_mcparticles[iDataset]);

    // H1D_trackPt_efficiency_splitCorrected[iDataset] = (TH1D*)H1D_trackPt_efficiency[iDataset]->Clone("trackEta_efficiency_splitCorrected"+Datasets[iDataset]);
    // H1D_trackPt_efficiency_splitCorrected[iDataset]->Add(H1D_trackPt_efficiency_split[iDataset],-1.);

    // first naming the histogram and resetting it to have a chosen unique name
    H1D_trackPt_efficiency_splitCorrected[iDataset] = (TH1D*)H1D_trackPt_mcparticles[iDataset]->Clone("trackPt_efficiency_split"+Datasets[iDataset]);
    H1D_trackPt_efficiency_splitCorrected[iDataset]->Reset("M");
    divideSuccess_split = H1D_trackPt_efficiency_splitCorrected[iDataset]->Divide(H1D_trackPt_assoctracks_split[iDataset], H1D_trackPt_mcparticles[iDataset]);
    H1D_trackPt_efficiency_splitCorrected[iDataset]->Scale(-1.);
    H1D_trackPt_efficiency_splitCorrected[iDataset]->Add(H1D_trackPt_efficiency[iDataset],1.);
    cout << "test2.6" << endl;

    // H1D_trackPt_efficiency_splitAndSecondaryCorrected[iDataset] = (TH1D*)H1D_trackPt_mcparticles[iDataset]->Clone("trackPt_efficiency_splitAndSecondaryCorrected"+Datasets[iDataset]);
    // H1D_trackPt_efficiency_splitAndSecondaryCorrected[iDataset]->Reset("M");
    // divideSuccess_split = H1D_trackPt_efficiency_splitAndSecondaryCorrected[iDataset]->Divide(H1D_trackPt_assoctracks_nonSplitSecondary[iDataset], H1D_trackPt_mcparticles[iDataset]);
    // H1D_trackPt_efficiency_splitAndSecondaryCorrected[iDataset]->Scale(-1.);
    // H1D_trackPt_efficiency_splitAndSecondaryCorrected[iDataset]->Add(H1D_trackPt_efficiency_splitCorrected[iDataset],1.);
    cout << "test2.7" << endl;
  }

  TString* pdfNameEntriesNorm = new TString("track_Pt_efficiency"+dummyName[0]+"_@eta["+Form("%.1f", etaRange[0])+","+Form("%.1f", etaRange[1])+"]");
  TString* pdfNameEntriesNorm_split = new TString("track_Pt_efficiency"+dummyName[0]+"_@eta["+Form("%.1f", etaRange[0])+","+Form("%.1f", etaRange[1])+"]_splitTracks");
  TString* pdfNameEntriesNorm_splitCorrected = new TString("track_Pt_efficiency"+dummyName[0]+"_@eta["+Form("%.1f", etaRange[0])+","+Form("%.1f", etaRange[1])+"]_splitCorrected");
  // TString* pdfNameEntriesNorm_splitAndSecondaryCorrected = new TString("track_Pt_efficiency"+dummyName[0]+"_@eta["+Form("%.1f", etaRange[0])+","+Form("%.1f", etaRange[1])+"]_splitAndSecondaryCorrected");

  TString textContext(contextCustomTwoFields(*texDatasetsComparisonCommonDenominator, contextEtaRange(etaRange), ""));

  if (divideSuccess == true) {
    Draw_TH1_Histograms_in_one(H1D_trackPt_efficiency, DatasetsNames, nDatasets, textContext, pdfNameEntriesNorm, texPtMC, texTrackEfficiency, texCollisionDataInfo, drawnWindowAuto, "logx,efficiency,150MevLine");
  }
  else {
    cout << "Divide failed in Draw_Efficiency_Pt_DatasetComparison" << endl;
  }
  Draw_TH1_Histograms_in_one(H1D_trackPt_efficiency_splitCorrected, DatasetsNames, nDatasets, textContext, pdfNameEntriesNorm_splitCorrected, texPtMC, texTrackEfficiency, texCollisionDataInfo, drawnWindowAuto, "logx,efficiency,150MevLine");
  // Draw_TH1_Histograms_in_one(H1D_trackPt_efficiency_splitAndSecondaryCorrected, DatasetsNames, nDatasets, textContext, pdfNameEntriesNorm_splitAndSecondaryCorrected, texPtMC, texTrackEfficiency, texCollisionDataInfo, drawnWindowAuto, "logx,efficiency,150MevLine");
}


void Draw_Efficiency_Eta_DatasetComparison_legacy(float* ptRange) {

  TH3D* H3D_trackPtEtaPhi_mcparticles[nDatasets];
  TH3D* H3D_trackPtEtaPhi_assoctracks[nDatasets];
  TH3D* H3D_trackPtEtaPhi_assoctracks_split[nDatasets];
  TH3D* H3D_trackPtEtaPhi_assoctracks_nonSplitSecondary[nDatasets];

  TH1D* H1D_trackEta_mcparticles[nDatasets];

  TH1D* H1D_trackEta_assoctracks[nDatasets];
  TH1D* H1D_trackEta_assoctracks_split[nDatasets];
  TH1D* H1D_trackEta_assoctracks_nonSplitSecondary[nDatasets];

  TH1D* H1D_trackEta_efficiency[nDatasets];
  // TH1D* H1D_trackEta_efficiency_split[nDatasets];
  TH1D* H1D_trackEta_efficiency_splitCorrected[nDatasets];
  TH1D* H1D_trackEta_efficiency_splitAndSecondaryCorrected[nDatasets];

  bool divideSuccess = false;
  bool divideSuccess_split = false;

  for(int iDataset = 0; iDataset < nDatasets; iDataset++){

    H3D_trackPtEtaPhi_mcparticles[iDataset] = (TH3D*)((TH3D*)file_O2Analysis_list[iDataset]->Get(analysisWorkflow[iDataset]+"/h3_track_pt_track_eta_track_phi_mcparticles"))->Clone("Draw_Efficiency_Eta_DatasetComparison_mcparticles"+Datasets[iDataset]);
    H3D_trackPtEtaPhi_assoctracks[iDataset] = (TH3D*)((TH3D*)file_O2Analysis_list[iDataset]->Get(analysisWorkflow[iDataset]+"/h3_track_pt_track_eta_track_phi_associatedtrackSelColl"))->Clone("Draw_Efficiency_Eta_DatasetComparison_associatedtrackSelColl"+Datasets[iDataset]);
    H3D_trackPtEtaPhi_assoctracks_split[iDataset] = (TH3D*)((TH3D*)file_O2Analysis_list[iDataset]->Get(analysisWorkflow[iDataset]+"/h3_track_pt_track_eta_track_phi_associatedtrackSelCollSplit"))->Clone("Draw_Efficiency_Eta_DatasetComparison_associatedtrackSelCollSplit"+Datasets[iDataset]);
    // H3D_trackPtEtaPhi_assoctracks_nonSplitSecondary[iDataset] = (TH3D*)((TH3D*)file_O2Analysis_list[iDataset]->Get(analysisWorkflow[iDataset]+"/h3_track_pt_track_eta_track_phi_associatedtrackSelCollNonSplitNonPrimary"))->Clone("Draw_Efficiency_Eta_DatasetComparison_associatedtrackSelCollNonSplitSecondary"+Datasets[iDataset]);
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
    // H3D_trackPtEtaPhi_mcparticles[iDataset]->GetXaxis()->SetRange(ibinPt_low,ibinPt_high);
    // H3D_trackPtEtaPhi_assoctracks[iDataset]->GetXaxis()->SetRange(ibinPt_low,ibinPt_high);
    // H3D_trackPtEtaPhi_assoctracks_split[iDataset]->GetXaxis()->SetRange(ibinPt_low,ibinPt_high);
    // H3D_trackPtEtaPhi_assoctracks_nonSplitSecondary[iDataset]->GetXaxis()->SetRange(ibinPt_low,ibinPt_high);

    H1D_trackEta_mcparticles[iDataset] = (TH1D*)H3D_trackPtEtaPhi_mcparticles[iDataset]->ProjectionY("trackEta_mcparticles"+Datasets[iDataset]+"_@pt["+Form("%.2f", ptRange[0])+","+Form("%.1f", ptRange[1])+"]", ibinPt_low, ibinPt_high, 0, -1, "e");
    H1D_trackEta_assoctracks[iDataset] = (TH1D*)H3D_trackPtEtaPhi_assoctracks[iDataset]->ProjectionY("trackEta_assoctracks"+Datasets[iDataset]+"_@pt["+Form("%.2f", ptRange[0])+","+Form("%.1f", ptRange[1])+"]", ibinPt_low, ibinPt_high, 0, -1, "e");
    H1D_trackEta_assoctracks_split[iDataset] = (TH1D*)H3D_trackPtEtaPhi_assoctracks_split[iDataset]->ProjectionY("trackEta_assoctracks_split"+Datasets[iDataset]+"_@pt["+Form("%.2f", ptRange[0])+","+Form("%.1f", ptRange[1])+"]", ibinPt_low, ibinPt_high, 0, -1, "e");
    // H1D_trackEta_assoctracks_nonSplitSecondary[iDataset] = (TH1D*)H3D_trackPtEtaPhi_assoctracks_nonSplitSecondary[iDataset]->ProjectionY("trackEta_assoctracks_nonSplitSecondary"+Datasets[iDataset]+"_@pt["+Form("%.2f", ptRange[0])+","+Form("%.1f", ptRange[1])+"]", ibinPt_low, ibinPt_high, 0, -1, "e");
    // H1D_trackEta_rebinned[iDataset] = (TH1D*)H1D_trackEta[iDataset]->Rebin(1.,"trackEta_rebinned"+Datasets[iDataset]);

    // NormaliseYieldToNEntries(H1D_trackEta_rebinned[iDataset]);
    // H1D_trackEta_rebinned[iDataset]->Scale(1./H1D_trackEta_rebinned[iDataset]->Integral(1, H1D_trackEta_rebinned[iDataset]->GetNbinsX()),"width"); // If option contains "width" the bin contents and errors are divided by the bin width.
    // H1D_trackEta_rebinned[iDataset]->Scale(1./H1D_trackEta_rebinned[iDataset]->Integral(ibincutlow, ibincuthigh),"width"); // If option contains "width" the bin contents and errors are divided by the bin width.

    H1D_trackEta_efficiency[iDataset] = (TH1D*)H1D_trackEta_mcparticles[iDataset]->Clone("trackEta_efficiency"+Datasets[iDataset]);
    H1D_trackEta_efficiency[iDataset]->Reset("M");
    divideSuccess = H1D_trackEta_efficiency[iDataset]->Divide(H1D_trackEta_assoctracks[iDataset], H1D_trackEta_mcparticles[iDataset]);

    // H1D_trackEta_efficiency_split[iDataset] = (TH1D*)H1D_trackEta_mcparticles[iDataset]->Clone("trackEta_efficiency_split"+Datasets[iDataset]);
    // H1D_trackEta_efficiency_split[iDataset]->Reset("M");
    // divideSuccess_split = H1D_trackEta_efficiency_split[iDataset]->Divide(H1D_trackEta_assoctracks_split[iDataset], H1D_trackEta_mcparticles[iDataset]);

    // H1D_trackEta_efficiency_splitCorrected[iDataset] = (TH1D*)H1D_trackEta_efficiency[iDataset]->Clone("trackEta_efficiency_splitCorrected"+Datasets[iDataset]);
    // H1D_trackEta_efficiency_splitCorrected[iDataset]->Add(H1D_trackEta_efficiency_split[iDataset],-1.);


    H1D_trackEta_efficiency_splitCorrected[iDataset] = (TH1D*)H1D_trackEta_mcparticles[iDataset]->Clone("trackEta_efficiency_split"+Datasets[iDataset]);
    H1D_trackEta_efficiency_splitCorrected[iDataset]->Reset("M");
    divideSuccess_split = H1D_trackEta_efficiency_splitCorrected[iDataset]->Divide(H1D_trackEta_assoctracks_split[iDataset], H1D_trackEta_mcparticles[iDataset]);
    H1D_trackEta_efficiency_splitCorrected[iDataset]->Scale(-1.);
    H1D_trackEta_efficiency_splitCorrected[iDataset]->Add(H1D_trackEta_efficiency[iDataset],1.);

    // H1D_trackEta_efficiency_splitAndSecondaryCorrected[iDataset] = (TH1D*)H1D_trackEta_mcparticles[iDataset]->Clone("trackEta_efficiency_splitAndSecondaryCorrected"+Datasets[iDataset]);
    // H1D_trackEta_efficiency_splitAndSecondaryCorrected[iDataset]->Reset("M");
    // divideSuccess_split = H1D_trackEta_efficiency_splitAndSecondaryCorrected[iDataset]->Divide(H1D_trackEta_assoctracks_nonSplitSecondary[iDataset], H1D_trackEta_mcparticles[iDataset]);
    // H1D_trackEta_efficiency_splitAndSecondaryCorrected[iDataset]->Scale(-1.);
    // H1D_trackEta_efficiency_splitAndSecondaryCorrected[iDataset]->Add(H1D_trackEta_efficiency_splitCorrected[iDataset],1.);
  }

  TString* pdfNameEntriesNorm = new TString("track_Eta_efficiency"+dummyName[0]+"_@pt["+Form("%.2f", ptRange[0])+","+Form("%.1f", ptRange[1])+"]");
  TString* pdfNameEntriesNorm_split = new TString("track_Eta_efficiency"+dummyName[0]+"_@pt["+Form("%.2f", ptRange[0])+","+Form("%.1f", ptRange[1])+"]_splitTracks");
  TString* pdfNameEntriesNorm_splitCorrected = new TString("track_Eta_efficiency"+dummyName[0]+"_@pt["+Form("%.2f", ptRange[0])+","+Form("%.1f", ptRange[1])+"]_splitCorrected");
  // TString* pdfNameEntriesNorm_splitAndSecondaryCorrected = new TString("track_Eta_efficiency"+dummyName[0]+"_@pt["+Form("%.2f", ptRange[0])+","+Form("%.1f", ptRange[1])+"]_splitAndSecondaryCorrected");

  // TString textContext(contextDatasetComp(""));
  TString textContext(contextCustomTwoFields(*texDatasetsComparisonCommonDenominator, contextPtRange(ptRange), ""));

  if (divideSuccess == true) {
    Draw_TH1_Histograms_in_one(H1D_trackEta_efficiency, DatasetsNames, nDatasets, textContext, pdfNameEntriesNorm, texEtaX, texTrackEfficiency, texCollisionDataInfo, drawnWindowAuto, "efficiency");
  }
  else {
    cout << "Divide failed in Draw_Efficiency_Eta_DatasetComparison" << endl;
  }
  Draw_TH1_Histograms_in_one(H1D_trackEta_efficiency_splitCorrected, DatasetsNames, nDatasets, textContext, pdfNameEntriesNorm_splitCorrected, texEtaX, texTrackEfficiency, texCollisionDataInfo, drawnWindowAuto, "efficiency");
  // Draw_TH1_Histograms_in_one(H1D_trackEta_efficiency_splitAndSecondaryCorrected, DatasetsNames, nDatasets, textContext, pdfNameEntriesNorm_splitAndSecondaryCorrected, texEtaX, texTrackEfficiency, texCollisionDataInfo, drawnWindowAuto, "efficiency");
}

void Draw_Efficiency_Phi_DatasetComparison_legacy(float* ptRange, float* etaRange) {

  TH3D* H3D_trackPtEtaPhi_mcparticles[nDatasets];
  TH3D* H3D_trackPtEtaPhi_assoctracks[nDatasets];
  TH3D* H3D_trackPtEtaPhi_assoctracks_split[nDatasets];
  TH3D* H3D_trackPtEtaPhi_assoctracks_nonSplitSecondary[nDatasets];

  TH1D* H1D_trackPhi_mcparticles[nDatasets];

  TH1D* H1D_trackPhi_assoctracks[nDatasets];
  TH1D* H1D_trackPhi_assoctracks_split[nDatasets];
  TH1D* H1D_trackPhi_assoctracks_nonSplitSecondary[nDatasets];
  
  TH1D* H1D_trackPhi_efficiency[nDatasets];
  // TH1D* H1D_trackPhi_efficiency_split[nDatasets];
  TH1D* H1D_trackPhi_efficiency_splitCorrected[nDatasets];
  TH1D* H1D_trackPhi_efficiency_splitAndSecondaryCorrected[nDatasets];

  bool divideSuccess = false;
  bool divideSuccess_split = false;

  for(int iDataset = 0; iDataset < nDatasets; iDataset++){

    H3D_trackPtEtaPhi_mcparticles[iDataset] = (TH3D*)((TH3D*)file_O2Analysis_list[iDataset]->Get(analysisWorkflow[iDataset]+"/h3_track_pt_track_eta_track_phi_mcparticles"))->Clone("Draw_Efficiency_Phi_DatasetComparison_mcparticles"+Datasets[iDataset]);
    H3D_trackPtEtaPhi_assoctracks[iDataset] = (TH3D*)((TH3D*)file_O2Analysis_list[iDataset]->Get(analysisWorkflow[iDataset]+"/h3_track_pt_track_eta_track_phi_associatedtrackSelColl"))->Clone("Draw_Efficiency_Phi_DatasetComparison_associatedtrackSelColl"+Datasets[iDataset]);
    H3D_trackPtEtaPhi_assoctracks_split[iDataset] = (TH3D*)((TH3D*)file_O2Analysis_list[iDataset]->Get(analysisWorkflow[iDataset]+"/h3_track_pt_track_eta_track_phi_associatedtrackSelCollSplit"))->Clone("Draw_Efficiency_Phi_DatasetComparison_associatedtrackSelCollSplit"+Datasets[iDataset]);
    // H3D_trackPtEtaPhi_assoctracks_nonSplitSecondary[iDataset] = (TH3D*)((TH3D*)file_O2Analysis_list[iDataset]->Get(analysisWorkflow[iDataset]+"/h3_track_pt_track_eta_track_phi_associatedtrackSelCollNonSplitNonPrimary"))->Clone("Draw_Efficiency_Phi_DatasetComparison_associatedtrackSelCollNonSplitSecondary"+Datasets[iDataset]);
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
    // H3D_trackPtEtaPhi_mcparticles[iDataset]->GetXaxis()->SetRange(ibinPt_low,ibinPt_high);
    // H3D_trackPtEtaPhi_assoctracks[iDataset]->GetXaxis()->SetRange(ibinPt_low,ibinPt_high);
    // H3D_trackPtEtaPhi_assoctracks_split[iDataset]->GetXaxis()->SetRange(ibinPt_low,ibinPt_high);
    // H3D_trackPtEtaPhi_assoctracks_nonSplitSecondary[iDataset]->GetXaxis()->SetRange(ibinPt_low,ibinPt_high);
    int ibinEta_low_mcpart = H3D_trackPtEtaPhi_mcparticles[iDataset]->GetYaxis()->FindBin(etaRange[0]+deltaEtaMcVsTrackEfficiency);
    int ibinEta_high_mcpart = H3D_trackPtEtaPhi_mcparticles[iDataset]->GetYaxis()->FindBin(etaRange[1]-deltaEtaMcVsTrackEfficiency);
    int ibinEta_low_track = H3D_trackPtEtaPhi_mcparticles[iDataset]->GetYaxis()->FindBin(etaRange[0]);
    int ibinEta_high_track = H3D_trackPtEtaPhi_mcparticles[iDataset]->GetYaxis()->FindBin(etaRange[1]);
    // H3D_trackPtEtaPhi_mcparticles[iDataset]->GetXaxis()->SetRange(ibinEta_low_mcpart,ibinEta_high_mcpart);
    // H3D_trackPtEtaPhi_assoctracks[iDataset]->GetXaxis()->SetRange(ibinEta_low_track,ibinEta_high_track);
    // H3D_trackPtEtaPhi_assoctracks_split[iDataset]->GetXaxis()->SetRange(ibinEta_low_track,ibinEta_high_track);
    // H3D_trackPtEtaPhi_assoctracks_nonSplitSecondary[iDataset]->GetXaxis()->SetRange(ibinEta_low_track,ibinEta_high_track);

    H1D_trackPhi_mcparticles[iDataset] = (TH1D*)H3D_trackPtEtaPhi_mcparticles[iDataset]->ProjectionZ("trackPhi_mcparticles"+Datasets[iDataset]+"_@pt["+Form("%.2f", ptRange[0])+","+Form("%.1f", ptRange[1])+"]", ibinPt_low, ibinPt_high, ibinEta_low_mcpart, ibinEta_high_mcpart, "e");
    H1D_trackPhi_assoctracks[iDataset] = (TH1D*)H3D_trackPtEtaPhi_assoctracks[iDataset]->ProjectionZ("trackPhi_assoctracks"+Datasets[iDataset]+"_@pt["+Form("%.2f", ptRange[0])+","+Form("%.1f", ptRange[1])+"]", ibinPt_low, ibinPt_high, ibinEta_low_track, ibinEta_high_track, "e");
    H1D_trackPhi_assoctracks_split[iDataset] = (TH1D*)H3D_trackPtEtaPhi_assoctracks_split[iDataset]->ProjectionZ("trackPhi_assoctracks_split"+Datasets[iDataset]+"_@pt["+Form("%.2f", ptRange[0])+","+Form("%.1f", ptRange[1])+"]", ibinPt_low, ibinPt_high, ibinEta_low_track, ibinEta_high_track, "e");
    // H1D_trackPhi_assoctracks_nonSplitSecondary[iDataset] = (TH1D*)H3D_trackPtEtaPhi_assoctracks_nonSplitSecondary[iDataset]->ProjectionZ("trackPhi_assoctracks_nonSplitSecondary"+Datasets[iDataset]+"_@pt["+Form("%.2f", ptRange[0])+","+Form("%.1f", ptRange[1])+"]", ibinPt_low, ibinPt_high, ibinEta_low_track, ibinEta_high_track, "e");
    // H1D_trackPhi_rebinned[iDataset] = (TH1D*)H1D_trackPhi[iDataset]->Rebin(1.,"trackPhi_rebinned"+Datasets[iDataset]);

    // NormaliseYieldToNEntries(H1D_trackPhi_rebinned[iDataset]);
    // H1D_trackPhi_rebinned[iDataset]->Scale(1./H1D_trackPhi_rebinned[iDataset]->Integral(1, H1D_trackPhi_rebinned[iDataset]->GetNbinsX()),"width"); // If option contains "width" the bin contents and errors are divided by the bin width.
    // H1D_trackPhi_rebinned[iDataset]->Scale(1./H1D_trackPhi_rebinned[iDataset]->Integral(ibincutlow, ibincuthigh),"width"); // If option contains "width" the bin contents and errors are divided by the bin width.

    H1D_trackPhi_efficiency[iDataset] = (TH1D*)H1D_trackPhi_mcparticles[iDataset]->Clone("trackPhi_efficiency"+Datasets[iDataset]);
    H1D_trackPhi_efficiency[iDataset]->Reset("M");
    divideSuccess = H1D_trackPhi_efficiency[iDataset]->Divide(H1D_trackPhi_assoctracks[iDataset], H1D_trackPhi_mcparticles[iDataset]);

    // H1D_trackPhi_efficiency_split[iDataset] = (TH1D*)H1D_trackPhi_mcparticles[iDataset]->Clone("trackPhi_efficiency_split"+Datasets[iDataset]);
    // H1D_trackPhi_efficiency_split[iDataset]->Reset("M");
    // divideSuccess = H1D_trackPhi_efficiency_split[iDataset]->Divide(H1D_trackPhi_assoctracks_split[iDataset], H1D_trackPhi_mcparticles[iDataset]);

    // H1D_trackPhi_efficiency_splitCorrected[iDataset] = (TH1D*)H1D_trackPhi_efficiency[iDataset]->Clone("trackPhi_efficiency_splitCorrected"+Datasets[iDataset]);
    // H1D_trackPhi_efficiency_splitCorrected[iDataset]->Add(H1D_trackPhi_efficiency_split[iDataset],-1.);

    H1D_trackPhi_efficiency_splitCorrected[iDataset] = (TH1D*)H1D_trackPhi_mcparticles[iDataset]->Clone("trackPhi_efficiency_split"+Datasets[iDataset]);
    H1D_trackPhi_efficiency_splitCorrected[iDataset]->Reset("M");
    divideSuccess_split = H1D_trackPhi_efficiency_splitCorrected[iDataset]->Divide(H1D_trackPhi_assoctracks_split[iDataset], H1D_trackPhi_mcparticles[iDataset]);
    H1D_trackPhi_efficiency_splitCorrected[iDataset]->Scale(-1.);
    H1D_trackPhi_efficiency_splitCorrected[iDataset]->Add(H1D_trackPhi_efficiency[iDataset],1.);

    // H1D_trackPhi_efficiency_splitAndSecondaryCorrected[iDataset] = (TH1D*)H1D_trackPhi_mcparticles[iDataset]->Clone("trackPhi_efficiency_splitAndSecondaryCorrected"+Datasets[iDataset]);
    // H1D_trackPhi_efficiency_splitAndSecondaryCorrected[iDataset]->Reset("M");
    // divideSuccess_split = H1D_trackPhi_efficiency_splitAndSecondaryCorrected[iDataset]->Divide(H1D_trackPhi_assoctracks_nonSplitSecondary[iDataset], H1D_trackPhi_mcparticles[iDataset]);
    // H1D_trackPhi_efficiency_splitAndSecondaryCorrected[iDataset]->Scale(-1.);
    // H1D_trackPhi_efficiency_splitAndSecondaryCorrected[iDataset]->Add(H1D_trackPhi_efficiency_splitCorrected[iDataset],1.);
  }

  TString* pdfNameEntriesNorm = new TString("track_Phi_efficiency"+dummyName[0]+"_@pt["+Form("%.2f", ptRange[0])+","+Form("%.1f", ptRange[1])+"]"+"_@eta["+Form("%.1f", etaRange[0])+","+Form("%.1f", etaRange[1])+"]");
  TString* pdfNameEntriesNorm_split = new TString("track_Phi_efficiency"+dummyName[0]+"_@pt["+Form("%.2f", ptRange[0])+","+Form("%.1f", ptRange[1])+"]"+"_@eta["+Form("%.1f", etaRange[0])+","+Form("%.1f", etaRange[1])+"]_split");
  TString* pdfNameEntriesNorm_splitCorrected = new TString("track_Phi_efficiency"+dummyName[0]+"_@pt["+Form("%.2f", ptRange[0])+","+Form("%.1f", ptRange[1])+"]"+"_@eta["+Form("%.1f", etaRange[0])+","+Form("%.1f", etaRange[1])+"]_splitCorrected");
  // TString* pdfNameEntriesNorm_splitAndSecondaryCorrected = new TString("track_Phi_efficiency"+dummyName[0]+"_@pt["+Form("%.2f", ptRange[0])+","+Form("%.1f", ptRange[1])+"]"+"_@eta["+Form("%.1f", etaRange[0])+","+Form("%.1f", etaRange[1])+"]_splitAndSecondaryCorrected");

  // TString textContext(contextDatasetComp(""));
  TString textContext(contextCustomTwoFields(*texDatasetsComparisonCommonDenominator, "#splitline{"+contextPtRange(ptRange)+"}{"+contextEtaRange(etaRange)+"}", ""));

  if (divideSuccess == true) {
    Draw_TH1_Histograms_in_one(H1D_trackPhi_efficiency, DatasetsNames, nDatasets, textContext, pdfNameEntriesNorm, texPhiX, texTrackEfficiency, texCollisionDataInfo, drawnWindowAuto, "efficiency");
  }
  else {
    cout << "Divide failed in Draw_Efficiency_Phi_DatasetComparison" << endl;
  }
  Draw_TH1_Histograms_in_one(H1D_trackPhi_efficiency_splitCorrected, DatasetsNames, nDatasets, textContext, pdfNameEntriesNorm_splitCorrected, texPhiX, texTrackEfficiency, texCollisionDataInfo, drawnWindowAuto, "efficiency");
  // Draw_TH1_Histograms_in_one(H1D_trackPhi_efficiency_splitAndSecondaryCorrected, DatasetsNames, nDatasets, textContext, pdfNameEntriesNorm_splitAndSecondaryCorrected, texPhiX, texTrackEfficiency, texCollisionDataInfo, drawnWindowAuto, "efficiency");
}


void Draw_Efficiency_Phi_DatasetComparison_finerPhi_legacy(float* ptRange, float* etaRange) {

  TH3D* H3D_trackPtEtaPhi_mcparticles[nDatasets];
  TH3D* H3D_trackPtEtaPhi_assoctracks[nDatasets];
  TH3D* H3D_trackPtEtaPhi_assoctracks_split[nDatasets];
  TH3D* H3D_trackPtEtaPhi_assoctracks_nonSplitSecondary[nDatasets];

  TH1D* H1D_trackPhi_mcparticles[nDatasets];

  TH1D* H1D_trackPhi_assoctracks[nDatasets];
  TH1D* H1D_trackPhi_assoctracks_split[nDatasets];
  TH1D* H1D_trackPhi_assoctracks_nonSplitSecondary[nDatasets];
  
  TH1D* H1D_trackPhi_efficiency[nDatasets];
  // TH1D* H1D_trackPhi_efficiency_split[nDatasets];
  TH1D* H1D_trackPhi_efficiency_splitCorrected[nDatasets];
  TH1D* H1D_trackPhi_efficiency_splitAndSecondaryCorrected[nDatasets];

  bool divideSuccess = false;
  bool divideSuccess_split = false;

  for(int iDataset = 0; iDataset < nDatasets; iDataset++){

    H3D_trackPtEtaPhi_mcparticles[iDataset] = (TH3D*)((TH3D*)file_O2Analysis_list[iDataset]->Get(analysisWorkflow[iDataset]+"/h3_track_pt_track_eta_track_phi_mcparticles_finerPhi"))->Clone("Draw_Efficiency_Phi_DatasetComparison_mcparticles_finerPhi"+Datasets[iDataset]);
    H3D_trackPtEtaPhi_assoctracks[iDataset] = (TH3D*)((TH3D*)file_O2Analysis_list[iDataset]->Get(analysisWorkflow[iDataset]+"/h3_track_pt_track_eta_track_phi_associatedtrackSelColl_finerPhi"))->Clone("Draw_Efficiency_Phi_DatasetComparison_associatedtrackSelColl_finerPhi"+Datasets[iDataset]);
    H3D_trackPtEtaPhi_assoctracks_split[iDataset] = (TH3D*)((TH3D*)file_O2Analysis_list[iDataset]->Get(analysisWorkflow[iDataset]+"/h3_track_pt_track_eta_track_phi_associatedtrackSelCollSplit_finerPhi"))->Clone("Draw_Efficiency_Phi_DatasetComparison_associatedtrackSelCollSplit_finerPhi"+Datasets[iDataset]);
    // H3D_trackPtEtaPhi_assoctracks_nonSplitSecondary[iDataset] = (TH3D*)((TH3D*)file_O2Analysis_list[iDataset]->Get(analysisWorkflow[iDataset]+"/h3_track_pt_track_eta_track_phi_associatedtrackSelCollNonSplitNonPrimary"))->Clone("Draw_Efficiency_Phi_DatasetComparison_associatedtrackSelCollNonSplitSecondary"+Datasets[iDataset]);
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
    // H3D_trackPtEtaPhi_mcparticles[iDataset]->GetXaxis()->SetRange(ibinPt_low,ibinPt_high);
    // H3D_trackPtEtaPhi_assoctracks[iDataset]->GetXaxis()->SetRange(ibinPt_low,ibinPt_high);
    // H3D_trackPtEtaPhi_assoctracks_split[iDataset]->GetXaxis()->SetRange(ibinPt_low,ibinPt_high);
    // H3D_trackPtEtaPhi_assoctracks_nonSplitSecondary[iDataset]->GetXaxis()->SetRange(ibinPt_low,ibinPt_high);
    int ibinEta_low_mcpart = H3D_trackPtEtaPhi_mcparticles[iDataset]->GetYaxis()->FindBin(etaRange[0]+deltaEtaMcVsTrackEfficiency);
    int ibinEta_high_mcpart = H3D_trackPtEtaPhi_mcparticles[iDataset]->GetYaxis()->FindBin(etaRange[1]-deltaEtaMcVsTrackEfficiency);
    int ibinEta_low_track = H3D_trackPtEtaPhi_mcparticles[iDataset]->GetYaxis()->FindBin(etaRange[0]);
    int ibinEta_high_track = H3D_trackPtEtaPhi_mcparticles[iDataset]->GetYaxis()->FindBin(etaRange[1]);
    // H3D_trackPtEtaPhi_mcparticles[iDataset]->GetXaxis()->SetRange(ibinEta_low_mcpart,ibinEta_high_mcpart);
    // H3D_trackPtEtaPhi_assoctracks[iDataset]->GetXaxis()->SetRange(ibinEta_low_track,ibinEta_high_track);
    // H3D_trackPtEtaPhi_assoctracks_split[iDataset]->GetXaxis()->SetRange(ibinEta_low_track,ibinEta_high_track);
    // H3D_trackPtEtaPhi_assoctracks_nonSplitSecondary[iDataset]->GetXaxis()->SetRange(ibinEta_low_track,ibinEta_high_track);

    H1D_trackPhi_mcparticles[iDataset] = (TH1D*)H3D_trackPtEtaPhi_mcparticles[iDataset]->ProjectionZ("trackPhi_mcparticles_finerPhi"+Datasets[iDataset]+"_@pt["+Form("%.2f", ptRange[0])+","+Form("%.1f", ptRange[1])+"]", ibinPt_low, ibinPt_high, ibinEta_low_mcpart, ibinEta_high_mcpart, "e");
    H1D_trackPhi_assoctracks[iDataset] = (TH1D*)H3D_trackPtEtaPhi_assoctracks[iDataset]->ProjectionZ("trackPhi_assoctracks_finerPhi"+Datasets[iDataset]+"_@pt["+Form("%.2f", ptRange[0])+","+Form("%.1f", ptRange[1])+"]", ibinPt_low, ibinPt_high, ibinEta_low_track, ibinEta_high_track, "e");
    H1D_trackPhi_assoctracks_split[iDataset] = (TH1D*)H3D_trackPtEtaPhi_assoctracks_split[iDataset]->ProjectionZ("trackPhi_assoctracks_nonSplitSecondary_finerPhi"+Datasets[iDataset]+"_@pt["+Form("%.2f", ptRange[0])+","+Form("%.1f", ptRange[1])+"]", ibinPt_low, ibinPt_high, ibinEta_low_track, ibinEta_high_track, "e");
    // H1D_trackPhi_assoctracks_nonSplitSecondary[iDataset] = (TH1D*)H3D_trackPtEtaPhi_assoctracks_nonSplitSecondary[iDataset]->ProjectionZ("trackPhi_assoctracks_split"+Datasets[iDataset]+"_@pt["+Form("%.2f", ptRange[0])+","+Form("%.1f", ptRange[1])+"]", ibinPt_low, ibinPt_high, ibinEta_low_track, ibinEta_high_track, "e");
    // H1D_trackPhi_rebinned[iDataset] = (TH1D*)H1D_trackPhi[iDataset]->Rebin(1.,"trackPhi_rebinned"+Datasets[iDataset]);

    // NormaliseYieldToNEntries(H1D_trackPhi_rebinned[iDataset]);
    // H1D_trackPhi_rebinned[iDataset]->Scale(1./H1D_trackPhi_rebinned[iDataset]->Integral(1, H1D_trackPhi_rebinned[iDataset]->GetNbinsX()),"width"); // If option contains "width" the bin contents and errors are divided by the bin width.
    // H1D_trackPhi_rebinned[iDataset]->Scale(1./H1D_trackPhi_rebinned[iDataset]->Integral(ibincutlow, ibincuthigh),"width"); // If option contains "width" the bin contents and errors are divided by the bin width.

    H1D_trackPhi_efficiency[iDataset] = (TH1D*)H1D_trackPhi_mcparticles[iDataset]->Clone("trackPhi_efficiency"+Datasets[iDataset]);
    H1D_trackPhi_efficiency[iDataset]->Reset("M");
    divideSuccess = H1D_trackPhi_efficiency[iDataset]->Divide(H1D_trackPhi_assoctracks[iDataset], H1D_trackPhi_mcparticles[iDataset]);

    // H1D_trackPhi_efficiency_split[iDataset] = (TH1D*)H1D_trackPhi_mcparticles[iDataset]->Clone("trackPhi_efficiency_split"+Datasets[iDataset]);
    // H1D_trackPhi_efficiency_split[iDataset]->Reset("M");
    // divideSuccess = H1D_trackPhi_efficiency_split[iDataset]->Divide(H1D_trackPhi_assoctracks_split[iDataset], H1D_trackPhi_mcparticles[iDataset]);

    // H1D_trackPhi_efficiency_splitCorrected[iDataset] = (TH1D*)H1D_trackPhi_efficiency[iDataset]->Clone("trackPhi_efficiency_splitCorrected"+Datasets[iDataset]);
    // H1D_trackPhi_efficiency_splitCorrected[iDataset]->Add(H1D_trackPhi_efficiency_split[iDataset],-1.);

    H1D_trackPhi_efficiency_splitCorrected[iDataset] = (TH1D*)H1D_trackPhi_mcparticles[iDataset]->Clone("trackPhi_efficiency_split"+Datasets[iDataset]);
    H1D_trackPhi_efficiency_splitCorrected[iDataset]->Reset("M");
    divideSuccess_split = H1D_trackPhi_efficiency_splitCorrected[iDataset]->Divide(H1D_trackPhi_assoctracks_split[iDataset], H1D_trackPhi_mcparticles[iDataset]);
    H1D_trackPhi_efficiency_splitCorrected[iDataset]->Scale(-1.);
    H1D_trackPhi_efficiency_splitCorrected[iDataset]->Add(H1D_trackPhi_efficiency[iDataset],1.);

    // H1D_trackPhi_efficiency_splitAndSecondaryCorrected[iDataset] = (TH1D*)H1D_trackPhi_mcparticles[iDataset]->Clone("trackPhi_efficiency_splitAndSecondaryCorrected"+Datasets[iDataset]);
    // H1D_trackPhi_efficiency_splitAndSecondaryCorrected[iDataset]->Reset("M");
    // divideSuccess_split = H1D_trackPhi_efficiency_splitAndSecondaryCorrected[iDataset]->Divide(H1D_trackPhi_assoctracks_nonSplitSecondary[iDataset], H1D_trackPhi_mcparticles[iDataset]);
    // H1D_trackPhi_efficiency_splitAndSecondaryCorrected[iDataset]->Scale(-1.);
    // H1D_trackPhi_efficiency_splitAndSecondaryCorrected[iDataset]->Add(H1D_trackPhi_efficiency_splitCorrected[iDataset],1.);
  }

  TString* pdfNameEntriesNorm = new TString("track_Phi_efficiency"+dummyName[0]+"_@pt["+Form("%.2f", ptRange[0])+","+Form("%.1f", ptRange[1])+"]"+"_@eta["+Form("%.1f", etaRange[0])+","+Form("%.1f", etaRange[1])+"]_finerPhi");
  TString* pdfNameEntriesNorm_split = new TString("track_Phi_efficiency"+dummyName[0]+"_@pt["+Form("%.2f", ptRange[0])+","+Form("%.1f", ptRange[1])+"]"+"_@eta["+Form("%.1f", etaRange[0])+","+Form("%.1f", etaRange[1])+"]_split_finerPhi");
  TString* pdfNameEntriesNorm_splitCorrected = new TString("track_Phi_efficiency"+dummyName[0]+"_@pt["+Form("%.2f", ptRange[0])+","+Form("%.1f", ptRange[1])+"]"+"_@eta["+Form("%.1f", etaRange[0])+","+Form("%.1f", etaRange[1])+"]_splitCorrected_finerPhi");
  // TString* pdfNameEntriesNorm_splitAndSecondaryCorrected = new TString("track_Phi_efficiency"+dummyName[0]+"_@pt["+Form("%.2f", ptRange[0])+","+Form("%.1f", ptRange[1])+"]"+"_@eta["+Form("%.1f", etaRange[0])+","+Form("%.1f", etaRange[1])+"]_splitAndSecondaryCorrected");

  // TString textContext(contextDatasetComp(""));
  TString textContext(contextCustomTwoFields(*texDatasetsComparisonCommonDenominator, "#splitline{"+contextPtRange(ptRange)+"}{"+contextEtaRange(etaRange)+"}", ""));

  if (divideSuccess == true) {
    Draw_TH1_Histograms_in_one(H1D_trackPhi_efficiency, DatasetsNames, nDatasets, textContext, pdfNameEntriesNorm, texPhiX, texTrackEfficiency, texCollisionDataInfo, drawnWindowAuto, "efficiency");
  }
  else {
    cout << "Divide failed in Draw_Efficiency_Phi_DatasetComparison" << endl;
  }
  Draw_TH1_Histograms_in_one(H1D_trackPhi_efficiency_splitCorrected, DatasetsNames, nDatasets, textContext, pdfNameEntriesNorm_splitCorrected, texPhiX, texTrackEfficiency, texCollisionDataInfo, drawnWindowAuto, "efficiency");
  // Draw_TH1_Histograms_in_one(H1D_trackPhi_efficiency_splitAndSecondaryCorrected, DatasetsNames, nDatasets, textContext, pdfNameEntriesNorm_splitAndSecondaryCorrected, texPhiX, texTrackEfficiency, texCollisionDataInfo, drawnWindowAuto, "efficiency");
}

void Draw_Efficiency_Eta_PtRangeComparison_legacy(float* ptRange, int nPtRanges, int iDataset, const char options[]) {

  TH3D* H3D_trackPtEtaPhi_mcparticles;
  TH3D* H3D_trackPtEtaPhi_assoctracks;

  H3D_trackPtEtaPhi_mcparticles = (TH3D*)((TH3D*)file_O2Analysis_list[iDataset]->Get(analysisWorkflow[iDataset]+"/h3_track_pt_track_eta_track_phi_mcparticles"))->Clone("Draw_Efficiency_Eta_PtRangeComparison_mcparticles"+Datasets[iDataset]);
  H3D_trackPtEtaPhi_assoctracks = (TH3D*)((TH3D*)file_O2Analysis_list[iDataset]->Get(analysisWorkflow[iDataset]+"/h3_track_pt_track_eta_track_phi_associatedtrackSelColl"))->Clone("Draw_Efficiency_Eta_PtRangeComparison_associatedtrackSelColl"+Datasets[iDataset]);
  // H3D_jetRjetPtjetEta = (TH3D*)((TH3D*)file_O2Analysis_list[iDataset]->Get(analysisWorkflow[iDataset]+"/h3_jet_r_jet_pt_jet_eta"+jetFinderQaHistType[iJetFinderQaType]))->Clone("Draw_Eta_PtCutComparison"+Datasets[iDataset]);

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
    Draw_TH1_Histograms_in_one(H1D_trackEta_efficiency, PtCutsLegend, nPtRanges, textContext, pdfNameEntriesNorm, texEtaX, texTrackEfficiency, texCollisionDataInfo, drawnWindowAuto, "");
  }
  else {
    cout << "Divide failed in Draw_Efficiency_Eta_DatasetComparison" << endl;
  }
}


void Draw_Efficiency_Phi_PtRangeComparison_legacy(float* ptRange, int nPtRanges, float* etaRange, int iDataset, const char options[]) {

  TH3D* H3D_trackPtEtaPhi_mcparticles;
  TH3D* H3D_trackPtEtaPhi_assoctracks;

  H3D_trackPtEtaPhi_mcparticles = (TH3D*)((TH3D*)file_O2Analysis_list[iDataset]->Get(analysisWorkflow[iDataset]+"/h3_track_pt_track_eta_track_phi_mcparticles"))->Clone("Draw_Efficiency_Phi_PtRangeComparison_mcparticles"+Datasets[iDataset]);
  H3D_trackPtEtaPhi_assoctracks = (TH3D*)((TH3D*)file_O2Analysis_list[iDataset]->Get(analysisWorkflow[iDataset]+"/h3_track_pt_track_eta_track_phi_associatedtrackSelColl"))->Clone("Draw_Efficiency_Phi_PtRangeComparison_associatedtrackSelColl"+Datasets[iDataset]);
  // H3D_jetRjetPtjetEta = (TH3D*)((TH3D*)file_O2Analysis_list[iDataset]->Get(analysisWorkflow[iDataset]+"/h3_jet_r_jet_pt_jet_eta"+jetFinderQaHistType[iJetFinderQaType]))->Clone("Draw_Phi_PtCutComparison"+Datasets[iDataset]);

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
    Draw_TH1_Histograms_in_one(H1D_trackPhi_efficiency, PtCutsLegend, nPtRanges, textContext, pdfNameEntriesNorm, texPhiX, texTrackEfficiency, texCollisionDataInfo, drawnWindowAuto, "");
  }
  else {
    cout << "Divide failed in Draw_Efficiency_Phi_DatasetComparison" << endl;
  }
}


void Draw_Efficiency_Pt_DatasetComparison(float* etaRange) {
  cout << "test0" << endl;
  TH3D* H3D_trackPtEtaPhi_mcparticles[nDatasets];
  TH3D* H3D_trackPtEtaPhi_assoctracks[nDatasets];
  TH3D* H3D_trackPtEtaPhi_assoctracks_split[nDatasets];
  TH3D* H3D_trackPtEtaPhi_assoctracks_nonSplitSecondary[nDatasets];

  TH1D* H1D_trackPt_mcparticles[nDatasets];

  TH1D* H1D_trackPt_assoctracks[nDatasets];
  TH1D* H1D_trackPt_assoctracks_split[nDatasets];
  TH1D* H1D_trackPt_assoctracks_nonSplitSecondary[nDatasets];
  
  TH1D* H1D_trackPt_efficiency[nDatasets];
  // TH1D* H1D_trackPt_efficiency_split[nDatasets];
  TH1D* H1D_trackPt_efficiency_splitCorrected[nDatasets];
  TH1D* H1D_trackPt_efficiency_splitAndSecondaryCorrected[nDatasets];

  bool divideSuccess = false;
  bool divideSuccess_split = false;

  cout << "test1" << endl;

  for(int iDataset = 0; iDataset < nDatasets; iDataset++){

    H3D_trackPtEtaPhi_mcparticles[iDataset] = (TH3D*)((TH3D*)file_O2Analysis_list[iDataset]->Get(analysisWorkflow[iDataset]+"/h3_particle_pt_particle_eta_particle_phi_mcpartofinterest"))->Clone("Draw_Efficiency_Pt_DatasetComparison_mcparticles"+Datasets[iDataset]);
    cout << "test2.1" << endl;
    H3D_trackPtEtaPhi_assoctracks[iDataset] = (TH3D*)((TH3D*)file_O2Analysis_list[iDataset]->Get(analysisWorkflow[iDataset]+"/h3_particle_pt_particle_eta_particle_phi_associatedtrack_primary"))->Clone("Draw_Efficiency_Pt_DatasetComparison_associatedtrackSelColl"+Datasets[iDataset]);
    cout << "test2.2" << endl;
    H3D_trackPtEtaPhi_assoctracks_split[iDataset] = (TH3D*)((TH3D*)file_O2Analysis_list[iDataset]->Get(analysisWorkflow[iDataset]+"/h3_particle_pt_particle_eta_particle_phi_associatedtrack_split"))->Clone("Draw_Efficiency_Pt_DatasetComparison_associatedtrackSelCollSplit"+Datasets[iDataset]);
    // H3D_trackPtEtaPhi_assoctracks_nonSplitSecondary[iDataset] = (TH3D*)((TH3D*)file_O2Analysis_list[iDataset]->Get(analysisWorkflow[iDataset]+"/h3_track_pt_track_eta_track_phi_associatedtrackSelCollNonSplitNonPrimary"))->Clone("Draw_Efficiency_Pt_DatasetComparison_associatedtrackSelCollNonSplitSecondary"+Datasets[iDataset]);
    cout << "test2.3" << endl;

    int ibinEta_low_mcpart = H3D_trackPtEtaPhi_mcparticles[iDataset]->GetYaxis()->FindBin(etaRange[0]+deltaEtaMcVsTrackEfficiency);
    int ibinEta_high_mcpart = H3D_trackPtEtaPhi_mcparticles[iDataset]->GetYaxis()->FindBin(etaRange[1]-deltaEtaMcVsTrackEfficiency);
    int ibinEta_low_track = H3D_trackPtEtaPhi_mcparticles[iDataset]->GetYaxis()->FindBin(etaRange[0]);
    int ibinEta_high_track = H3D_trackPtEtaPhi_mcparticles[iDataset]->GetYaxis()->FindBin(etaRange[1]);

    // cout << "H3D_trackPtEtaPhi_assoctracks -  preSetRange bin count = " << H3D_trackPtEtaPhi_assoctracks[iDataset]->GetNbinsY() << endl;
    // H3D_trackPtEtaPhi_mcparticles[iDataset]->GetYaxis()->SetRange(ibinEta_low_mcpart,ibinEta_high_mcpart);
    // H3D_trackPtEtaPhi_assoctracks[iDataset]->GetYaxis()->SetRange(ibinEta_low_track,ibinEta_high_track);
    // H3D_trackPtEtaPhi_assoctracks_split[iDataset]->GetYaxis()->SetRange(ibinEta_low_track,ibinEta_high_track);
    // H3D_trackPtEtaPhi_assoctracks_nonSplitSecondary[iDataset]->GetYaxis()->SetRange(ibinEta_low_track,ibinEta_high_track);
    // cout << "H3D_trackPtEtaPhi_assoctracks - postSetRange bin count = " << H3D_trackPtEtaPhi_assoctracks[iDataset]->GetNbinsY() << endl;

    H1D_trackPt_mcparticles[iDataset] = (TH1D*)H3D_trackPtEtaPhi_mcparticles[iDataset]->ProjectionX("trackPt_mcparticles"+Datasets[iDataset], ibinEta_low_mcpart, ibinEta_high_mcpart, 0, -1, "e");
    cout << "test2.31" << endl;
    H1D_trackPt_assoctracks[iDataset] = (TH1D*)H3D_trackPtEtaPhi_assoctracks[iDataset]->ProjectionX("trackPt_tracks_split"+Datasets[iDataset], ibinEta_low_track, ibinEta_high_track, 0, -1, "e");
    cout << "test2.32" << endl;
    H1D_trackPt_assoctracks_split[iDataset] = (TH1D*)H3D_trackPtEtaPhi_assoctracks_split[iDataset]->ProjectionX("trackPt_assoctracks_split"+Datasets[iDataset], ibinEta_low_track, ibinEta_high_track, 0, -1, "e");
    cout << "test2.33" << endl;
    // H1D_trackPt_assoctracks_nonSplitSecondary[iDataset] = (TH1D*)H3D_trackPtEtaPhi_assoctracks_nonSplitSecondary[iDataset]->ProjectionX("trackPt_assoctracks_nonSplitSecondary"+Datasets[iDataset], ibinEta_low_track, ibinEta_high_track, 0, -1, "e");
    // cout << "test2.4" << endl;


    // naming the histogram and resetting it to have a chosen name
    H1D_trackPt_efficiency[iDataset] = (TH1D*)H1D_trackPt_mcparticles[iDataset]->Clone("trackPt_efficiency"+Datasets[iDataset]);
    H1D_trackPt_efficiency[iDataset]->Reset("M");
    cout << "  - " << DatasetsNames[iDataset] << ":" << endl;
    cout << "      tracks count = " << H1D_trackPt_assoctracks[iDataset]->GetEntries() << ", part count = " << H1D_trackPt_mcparticles[iDataset]->GetEntries() << endl;
    divideSuccess = H1D_trackPt_efficiency[iDataset]->Divide(H1D_trackPt_assoctracks[iDataset], H1D_trackPt_mcparticles[iDataset]);
    cout << "test2.5" << endl;

    // H1D_trackPt_efficiency_split[iDataset] = (TH1D*)H1D_trackPt_mcparticles[iDataset]->Clone("trackPt_efficiency_split"+Datasets[iDataset]);
    // H1D_trackPt_efficiency_split[iDataset]->Reset("M");
    // cout << "split tracks count = " << H1D_trackPt_assoctracks_split[iDataset]->GetEntries() << ", part count = " << H1D_trackPt_mcparticles[iDataset]->GetEntries() << endl;
    // divideSuccess_split = H1D_trackPt_efficiency_split[iDataset]->Divide(H1D_trackPt_assoctracks_split[iDataset], H1D_trackPt_mcparticles[iDataset]);

    // H1D_trackPt_efficiency_splitCorrected[iDataset] = (TH1D*)H1D_trackPt_efficiency[iDataset]->Clone("trackEta_efficiency_splitCorrected"+Datasets[iDataset]);
    // H1D_trackPt_efficiency_splitCorrected[iDataset]->Add(H1D_trackPt_efficiency_split[iDataset],-1.);

    // first naming the histogram and resetting it to have a chosen unique name
    H1D_trackPt_efficiency_splitCorrected[iDataset] = (TH1D*)H1D_trackPt_mcparticles[iDataset]->Clone("trackPt_efficiency_split"+Datasets[iDataset]);
    H1D_trackPt_efficiency_splitCorrected[iDataset]->Reset("M");
    divideSuccess_split = H1D_trackPt_efficiency_splitCorrected[iDataset]->Divide(H1D_trackPt_assoctracks_split[iDataset], H1D_trackPt_mcparticles[iDataset]);
    H1D_trackPt_efficiency_splitCorrected[iDataset]->Scale(-1.);
    H1D_trackPt_efficiency_splitCorrected[iDataset]->Add(H1D_trackPt_efficiency[iDataset],1.);
    cout << "test2.6" << endl;

    // H1D_trackPt_efficiency_splitAndSecondaryCorrected[iDataset] = (TH1D*)H1D_trackPt_mcparticles[iDataset]->Clone("trackPt_efficiency_splitAndSecondaryCorrected"+Datasets[iDataset]);
    // H1D_trackPt_efficiency_splitAndSecondaryCorrected[iDataset]->Reset("M");
    // divideSuccess_split = H1D_trackPt_efficiency_splitAndSecondaryCorrected[iDataset]->Divide(H1D_trackPt_assoctracks_nonSplitSecondary[iDataset], H1D_trackPt_mcparticles[iDataset]);
    // H1D_trackPt_efficiency_splitAndSecondaryCorrected[iDataset]->Scale(-1.);
    // H1D_trackPt_efficiency_splitAndSecondaryCorrected[iDataset]->Add(H1D_trackPt_efficiency_splitCorrected[iDataset],1.);
    cout << "test2.7" << endl;
  }

  TString* pdfNameEntriesNorm = new TString("track_Pt_efficiency"+dummyName[0]+"_@eta["+Form("%.1f", etaRange[0])+","+Form("%.1f", etaRange[1])+"]");
  TString* pdfNameEntriesNorm_split = new TString("track_Pt_efficiency"+dummyName[0]+"_@eta["+Form("%.1f", etaRange[0])+","+Form("%.1f", etaRange[1])+"]_splitTracks");
  TString* pdfNameEntriesNorm_splitCorrected = new TString("track_Pt_efficiency"+dummyName[0]+"_@eta["+Form("%.1f", etaRange[0])+","+Form("%.1f", etaRange[1])+"]_splitCorrected");
  // TString* pdfNameEntriesNorm_splitAndSecondaryCorrected = new TString("track_Pt_efficiency"+dummyName[0]+"_@eta["+Form("%.1f", etaRange[0])+","+Form("%.1f", etaRange[1])+"]_splitAndSecondaryCorrected");
  cout << "test3" << endl;

  TString textContext(contextCustomTwoFields(*texDatasetsComparisonCommonDenominator, contextEtaRange(etaRange), ""));
  cout << "test4" << endl;

  if (divideSuccess == true) {
    Draw_TH1_Histograms_in_one(H1D_trackPt_efficiency, DatasetsNames, nDatasets, textContext, pdfNameEntriesNorm, texPtMC, texTrackEfficiency, texCollisionDataInfo, drawnWindowAuto, "logx,efficiency,150MevLine");
  }
  else {
    cout << "Divide failed in Draw_Efficiency_Pt_DatasetComparison" << endl;
  }
  cout << "test5" << endl;
  Draw_TH1_Histograms_in_one(H1D_trackPt_efficiency_splitCorrected, DatasetsNames, nDatasets, textContext, pdfNameEntriesNorm_splitCorrected, texPtMC, texTrackEfficiency, texCollisionDataInfo, drawnWindowAuto, "logx,efficiency,150MevLine");
  // Draw_TH1_Histograms_in_one(H1D_trackPt_efficiency_splitAndSecondaryCorrected, DatasetsNames, nDatasets, textContext, pdfNameEntriesNorm_splitAndSecondaryCorrected, texPtMC, texTrackEfficiency, texCollisionDataInfo, drawnWindowAuto, "logx,efficiency,150MevLine");
}


void Draw_Efficiency_Eta_DatasetComparison(float* ptRange) {

  TH3D* H3D_trackPtEtaPhi_mcparticles[nDatasets];
  TH3D* H3D_trackPtEtaPhi_assoctracks[nDatasets];
  TH3D* H3D_trackPtEtaPhi_assoctracks_split[nDatasets];
  TH3D* H3D_trackPtEtaPhi_assoctracks_nonSplitSecondary[nDatasets];

  TH1D* H1D_trackEta_mcparticles[nDatasets];

  TH1D* H1D_trackEta_assoctracks[nDatasets];
  TH1D* H1D_trackEta_assoctracks_split[nDatasets];
  TH1D* H1D_trackEta_assoctracks_nonSplitSecondary[nDatasets];

  TH1D* H1D_trackEta_efficiency[nDatasets];
  // TH1D* H1D_trackEta_efficiency_split[nDatasets];
  TH1D* H1D_trackEta_efficiency_splitCorrected[nDatasets];
  TH1D* H1D_trackEta_efficiency_splitAndSecondaryCorrected[nDatasets];

  bool divideSuccess = false;
  bool divideSuccess_split = false;

  for(int iDataset = 0; iDataset < nDatasets; iDataset++){

    H3D_trackPtEtaPhi_mcparticles[iDataset] = (TH3D*)((TH3D*)file_O2Analysis_list[iDataset]->Get(analysisWorkflow[iDataset]+"/h3_particle_pt_particle_eta_particle_phi_mcpartofinterest"))->Clone("Draw_Efficiency_Eta_DatasetComparison_mcparticles"+Datasets[iDataset]);
    H3D_trackPtEtaPhi_assoctracks[iDataset] = (TH3D*)((TH3D*)file_O2Analysis_list[iDataset]->Get(analysisWorkflow[iDataset]+"/h3_particle_pt_particle_eta_particle_phi_associatedtrack_primary"))->Clone("Draw_Efficiency_Eta_DatasetComparison_associatedtrackSelColl"+Datasets[iDataset]);
    H3D_trackPtEtaPhi_assoctracks_split[iDataset] = (TH3D*)((TH3D*)file_O2Analysis_list[iDataset]->Get(analysisWorkflow[iDataset]+"/h3_particle_pt_particle_eta_particle_phi_associatedtrack_split"))->Clone("Draw_Efficiency_Eta_DatasetComparison_associatedtrackSelCollSplit"+Datasets[iDataset]);
    // H3D_trackPtEtaPhi_assoctracks_nonSplitSecondary[iDataset] = (TH3D*)((TH3D*)file_O2Analysis_list[iDataset]->Get(analysisWorkflow[iDataset]+"/h3_track_pt_track_eta_track_phi_associatedtrackSelCollNonSplitNonPrimary"))->Clone("Draw_Efficiency_Eta_DatasetComparison_associatedtrackSelCollNonSplitSecondary"+Datasets[iDataset]);
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
    // H3D_trackPtEtaPhi_mcparticles[iDataset]->GetXaxis()->SetRange(ibinPt_low,ibinPt_high);
    // H3D_trackPtEtaPhi_assoctracks[iDataset]->GetXaxis()->SetRange(ibinPt_low,ibinPt_high);
    // H3D_trackPtEtaPhi_assoctracks_split[iDataset]->GetXaxis()->SetRange(ibinPt_low,ibinPt_high);
    // H3D_trackPtEtaPhi_assoctracks_nonSplitSecondary[iDataset]->GetXaxis()->SetRange(ibinPt_low,ibinPt_high);

    H1D_trackEta_mcparticles[iDataset] = (TH1D*)H3D_trackPtEtaPhi_mcparticles[iDataset]->ProjectionY("trackEta_mcparticles"+Datasets[iDataset]+"_@pt["+Form("%.2f", ptRange[0])+","+Form("%.1f", ptRange[1])+"]", ibinPt_low, ibinPt_high, 0, -1, "e");
    H1D_trackEta_assoctracks[iDataset] = (TH1D*)H3D_trackPtEtaPhi_assoctracks[iDataset]->ProjectionY("trackEta_assoctracks"+Datasets[iDataset]+"_@pt["+Form("%.2f", ptRange[0])+","+Form("%.1f", ptRange[1])+"]", ibinPt_low, ibinPt_high, 0, -1, "e");
    H1D_trackEta_assoctracks_split[iDataset] = (TH1D*)H3D_trackPtEtaPhi_assoctracks_split[iDataset]->ProjectionY("trackEta_assoctracks_split"+Datasets[iDataset]+"_@pt["+Form("%.2f", ptRange[0])+","+Form("%.1f", ptRange[1])+"]", ibinPt_low, ibinPt_high, 0, -1, "e");
    // H1D_trackEta_assoctracks_nonSplitSecondary[iDataset] = (TH1D*)H3D_trackPtEtaPhi_assoctracks_nonSplitSecondary[iDataset]->ProjectionY("trackEta_assoctracks_nonSplitSecondary"+Datasets[iDataset]+"_@pt["+Form("%.2f", ptRange[0])+","+Form("%.1f", ptRange[1])+"]", ibinPt_low, ibinPt_high, 0, -1, "e");
    // H1D_trackEta_rebinned[iDataset] = (TH1D*)H1D_trackEta[iDataset]->Rebin(1.,"trackEta_rebinned"+Datasets[iDataset]);

    // NormaliseYieldToNEntries(H1D_trackEta_rebinned[iDataset]);
    // H1D_trackEta_rebinned[iDataset]->Scale(1./H1D_trackEta_rebinned[iDataset]->Integral(1, H1D_trackEta_rebinned[iDataset]->GetNbinsX()),"width"); // If option contains "width" the bin contents and errors are divided by the bin width.
    // H1D_trackEta_rebinned[iDataset]->Scale(1./H1D_trackEta_rebinned[iDataset]->Integral(ibincutlow, ibincuthigh),"width"); // If option contains "width" the bin contents and errors are divided by the bin width.

    H1D_trackEta_efficiency[iDataset] = (TH1D*)H1D_trackEta_mcparticles[iDataset]->Clone("trackEta_efficiency"+Datasets[iDataset]);
    H1D_trackEta_efficiency[iDataset]->Reset("M");
    divideSuccess = H1D_trackEta_efficiency[iDataset]->Divide(H1D_trackEta_assoctracks[iDataset], H1D_trackEta_mcparticles[iDataset]);

    // H1D_trackEta_efficiency_split[iDataset] = (TH1D*)H1D_trackEta_mcparticles[iDataset]->Clone("trackEta_efficiency_split"+Datasets[iDataset]);
    // H1D_trackEta_efficiency_split[iDataset]->Reset("M");
    // divideSuccess_split = H1D_trackEta_efficiency_split[iDataset]->Divide(H1D_trackEta_assoctracks_split[iDataset], H1D_trackEta_mcparticles[iDataset]);

    // H1D_trackEta_efficiency_splitCorrected[iDataset] = (TH1D*)H1D_trackEta_efficiency[iDataset]->Clone("trackEta_efficiency_splitCorrected"+Datasets[iDataset]);
    // H1D_trackEta_efficiency_splitCorrected[iDataset]->Add(H1D_trackEta_efficiency_split[iDataset],-1.);


    H1D_trackEta_efficiency_splitCorrected[iDataset] = (TH1D*)H1D_trackEta_mcparticles[iDataset]->Clone("trackEta_efficiency_split"+Datasets[iDataset]);
    H1D_trackEta_efficiency_splitCorrected[iDataset]->Reset("M");
    divideSuccess_split = H1D_trackEta_efficiency_splitCorrected[iDataset]->Divide(H1D_trackEta_assoctracks_split[iDataset], H1D_trackEta_mcparticles[iDataset]);
    H1D_trackEta_efficiency_splitCorrected[iDataset]->Scale(-1.);
    H1D_trackEta_efficiency_splitCorrected[iDataset]->Add(H1D_trackEta_efficiency[iDataset],1.);

    // H1D_trackEta_efficiency_splitAndSecondaryCorrected[iDataset] = (TH1D*)H1D_trackEta_mcparticles[iDataset]->Clone("trackEta_efficiency_splitAndSecondaryCorrected"+Datasets[iDataset]);
    // H1D_trackEta_efficiency_splitAndSecondaryCorrected[iDataset]->Reset("M");
    // divideSuccess_split = H1D_trackEta_efficiency_splitAndSecondaryCorrected[iDataset]->Divide(H1D_trackEta_assoctracks_nonSplitSecondary[iDataset], H1D_trackEta_mcparticles[iDataset]);
    // H1D_trackEta_efficiency_splitAndSecondaryCorrected[iDataset]->Scale(-1.);
    // H1D_trackEta_efficiency_splitAndSecondaryCorrected[iDataset]->Add(H1D_trackEta_efficiency_splitCorrected[iDataset],1.);
  }

  TString* pdfNameEntriesNorm = new TString("track_Eta_efficiency"+dummyName[0]+"_@pt["+Form("%.2f", ptRange[0])+","+Form("%.1f", ptRange[1])+"]");
  TString* pdfNameEntriesNorm_split = new TString("track_Eta_efficiency"+dummyName[0]+"_@pt["+Form("%.2f", ptRange[0])+","+Form("%.1f", ptRange[1])+"]_splitTracks");
  TString* pdfNameEntriesNorm_splitCorrected = new TString("track_Eta_efficiency"+dummyName[0]+"_@pt["+Form("%.2f", ptRange[0])+","+Form("%.1f", ptRange[1])+"]_splitCorrected");
  // TString* pdfNameEntriesNorm_splitAndSecondaryCorrected = new TString("track_Eta_efficiency"+dummyName[0]+"_@pt["+Form("%.2f", ptRange[0])+","+Form("%.1f", ptRange[1])+"]_splitAndSecondaryCorrected");

  // TString textContext(contextDatasetComp(""));
  TString textContext(contextCustomTwoFields(*texDatasetsComparisonCommonDenominator, contextPtRange(ptRange), ""));

  if (divideSuccess == true) {
    Draw_TH1_Histograms_in_one(H1D_trackEta_efficiency, DatasetsNames, nDatasets, textContext, pdfNameEntriesNorm, texEtaX, texTrackEfficiency, texCollisionDataInfo, drawnWindowAuto, "efficiency");
  }
  else {
    cout << "Divide failed in Draw_Efficiency_Eta_DatasetComparison" << endl;
  }
  Draw_TH1_Histograms_in_one(H1D_trackEta_efficiency_splitCorrected, DatasetsNames, nDatasets, textContext, pdfNameEntriesNorm_splitCorrected, texEtaX, texTrackEfficiency, texCollisionDataInfo, drawnWindowAuto, "efficiency");
  // Draw_TH1_Histograms_in_one(H1D_trackEta_efficiency_splitAndSecondaryCorrected, DatasetsNames, nDatasets, textContext, pdfNameEntriesNorm_splitAndSecondaryCorrected, texEtaX, texTrackEfficiency, texCollisionDataInfo, drawnWindowAuto, "efficiency");
}

void Draw_Efficiency_Phi_DatasetComparison(float* ptRange, float* etaRange) {

  TH3D* H3D_trackPtEtaPhi_mcparticles[nDatasets];
  TH3D* H3D_trackPtEtaPhi_assoctracks[nDatasets];
  TH3D* H3D_trackPtEtaPhi_assoctracks_split[nDatasets];
  TH3D* H3D_trackPtEtaPhi_assoctracks_nonSplitSecondary[nDatasets];

  TH1D* H1D_trackPhi_mcparticles[nDatasets];

  TH1D* H1D_trackPhi_assoctracks[nDatasets];
  TH1D* H1D_trackPhi_assoctracks_split[nDatasets];
  TH1D* H1D_trackPhi_assoctracks_nonSplitSecondary[nDatasets];
  
  TH1D* H1D_trackPhi_efficiency[nDatasets];
  // TH1D* H1D_trackPhi_efficiency_split[nDatasets];
  TH1D* H1D_trackPhi_efficiency_splitCorrected[nDatasets];
  TH1D* H1D_trackPhi_efficiency_splitAndSecondaryCorrected[nDatasets];

  bool divideSuccess = false;
  bool divideSuccess_split = false;

  for(int iDataset = 0; iDataset < nDatasets; iDataset++){

    H3D_trackPtEtaPhi_mcparticles[iDataset] = (TH3D*)((TH3D*)file_O2Analysis_list[iDataset]->Get(analysisWorkflow[iDataset]+"/h3_particle_pt_particle_eta_particle_phi_mcpartofinterest"))->Clone("Draw_Efficiency_Phi_DatasetComparison_mcparticles"+Datasets[iDataset]);
    H3D_trackPtEtaPhi_assoctracks[iDataset] = (TH3D*)((TH3D*)file_O2Analysis_list[iDataset]->Get(analysisWorkflow[iDataset]+"/h3_particle_pt_particle_eta_particle_phi_associatedtrack_primary"))->Clone("Draw_Efficiency_Phi_DatasetComparison_associatedtrackSelColl"+Datasets[iDataset]);
    H3D_trackPtEtaPhi_assoctracks_split[iDataset] = (TH3D*)((TH3D*)file_O2Analysis_list[iDataset]->Get(analysisWorkflow[iDataset]+"/h3_particle_pt_particle_eta_particle_phi_associatedtrack_split"))->Clone("Draw_Efficiency_Phi_DatasetComparison_associatedtrackSelCollSplit"+Datasets[iDataset]);
    // H3D_trackPtEtaPhi_assoctracks_nonSplitSecondary[iDataset] = (TH3D*)((TH3D*)file_O2Analysis_list[iDataset]->Get(analysisWorkflow[iDataset]+"/h3_track_pt_track_eta_track_phi_associatedtrackSelCollNonSplitNonPrimary"))->Clone("Draw_Efficiency_Phi_DatasetComparison_associatedtrackSelCollNonSplitSecondary"+Datasets[iDataset]);
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
    // H3D_trackPtEtaPhi_mcparticles[iDataset]->GetXaxis()->SetRange(ibinPt_low,ibinPt_high);
    // H3D_trackPtEtaPhi_assoctracks[iDataset]->GetXaxis()->SetRange(ibinPt_low,ibinPt_high);
    // H3D_trackPtEtaPhi_assoctracks_split[iDataset]->GetXaxis()->SetRange(ibinPt_low,ibinPt_high);
    // H3D_trackPtEtaPhi_assoctracks_nonSplitSecondary[iDataset]->GetXaxis()->SetRange(ibinPt_low,ibinPt_high);
    int ibinEta_low_mcpart = H3D_trackPtEtaPhi_mcparticles[iDataset]->GetYaxis()->FindBin(etaRange[0]+deltaEtaMcVsTrackEfficiency);
    int ibinEta_high_mcpart = H3D_trackPtEtaPhi_mcparticles[iDataset]->GetYaxis()->FindBin(etaRange[1]-deltaEtaMcVsTrackEfficiency);
    int ibinEta_low_track = H3D_trackPtEtaPhi_mcparticles[iDataset]->GetYaxis()->FindBin(etaRange[0]);
    int ibinEta_high_track = H3D_trackPtEtaPhi_mcparticles[iDataset]->GetYaxis()->FindBin(etaRange[1]);
    // H3D_trackPtEtaPhi_mcparticles[iDataset]->GetXaxis()->SetRange(ibinEta_low_mcpart,ibinEta_high_mcpart);
    // H3D_trackPtEtaPhi_assoctracks[iDataset]->GetXaxis()->SetRange(ibinEta_low_track,ibinEta_high_track);
    // H3D_trackPtEtaPhi_assoctracks_split[iDataset]->GetXaxis()->SetRange(ibinEta_low_track,ibinEta_high_track);
    // H3D_trackPtEtaPhi_assoctracks_nonSplitSecondary[iDataset]->GetXaxis()->SetRange(ibinEta_low_track,ibinEta_high_track);

    H1D_trackPhi_mcparticles[iDataset] = (TH1D*)H3D_trackPtEtaPhi_mcparticles[iDataset]->ProjectionZ("trackPhi_mcparticles"+Datasets[iDataset]+"_@pt["+Form("%.2f", ptRange[0])+","+Form("%.1f", ptRange[1])+"]", ibinPt_low, ibinPt_high, ibinEta_low_mcpart, ibinEta_high_mcpart, "e");
    H1D_trackPhi_assoctracks[iDataset] = (TH1D*)H3D_trackPtEtaPhi_assoctracks[iDataset]->ProjectionZ("trackPhi_assoctracks"+Datasets[iDataset]+"_@pt["+Form("%.2f", ptRange[0])+","+Form("%.1f", ptRange[1])+"]", ibinPt_low, ibinPt_high, ibinEta_low_track, ibinEta_high_track, "e");
    H1D_trackPhi_assoctracks_split[iDataset] = (TH1D*)H3D_trackPtEtaPhi_assoctracks_split[iDataset]->ProjectionZ("trackPhi_assoctracks_split"+Datasets[iDataset]+"_@pt["+Form("%.2f", ptRange[0])+","+Form("%.1f", ptRange[1])+"]", ibinPt_low, ibinPt_high, ibinEta_low_track, ibinEta_high_track, "e");
    // H1D_trackPhi_assoctracks_nonSplitSecondary[iDataset] = (TH1D*)H3D_trackPtEtaPhi_assoctracks_nonSplitSecondary[iDataset]->ProjectionZ("trackPhi_assoctracks_nonSplitSecondary"+Datasets[iDataset]+"_@pt["+Form("%.2f", ptRange[0])+","+Form("%.1f", ptRange[1])+"]", ibinPt_low, ibinPt_high, ibinEta_low_track, ibinEta_high_track, "e");
    // H1D_trackPhi_rebinned[iDataset] = (TH1D*)H1D_trackPhi[iDataset]->Rebin(1.,"trackPhi_rebinned"+Datasets[iDataset]);

    // NormaliseYieldToNEntries(H1D_trackPhi_rebinned[iDataset]);
    // H1D_trackPhi_rebinned[iDataset]->Scale(1./H1D_trackPhi_rebinned[iDataset]->Integral(1, H1D_trackPhi_rebinned[iDataset]->GetNbinsX()),"width"); // If option contains "width" the bin contents and errors are divided by the bin width.
    // H1D_trackPhi_rebinned[iDataset]->Scale(1./H1D_trackPhi_rebinned[iDataset]->Integral(ibincutlow, ibincuthigh),"width"); // If option contains "width" the bin contents and errors are divided by the bin width.

    H1D_trackPhi_efficiency[iDataset] = (TH1D*)H1D_trackPhi_mcparticles[iDataset]->Clone("trackPhi_efficiency"+Datasets[iDataset]);
    H1D_trackPhi_efficiency[iDataset]->Reset("M");
    divideSuccess = H1D_trackPhi_efficiency[iDataset]->Divide(H1D_trackPhi_assoctracks[iDataset], H1D_trackPhi_mcparticles[iDataset]);

    // H1D_trackPhi_efficiency_split[iDataset] = (TH1D*)H1D_trackPhi_mcparticles[iDataset]->Clone("trackPhi_efficiency_split"+Datasets[iDataset]);
    // H1D_trackPhi_efficiency_split[iDataset]->Reset("M");
    // divideSuccess = H1D_trackPhi_efficiency_split[iDataset]->Divide(H1D_trackPhi_assoctracks_split[iDataset], H1D_trackPhi_mcparticles[iDataset]);

    // H1D_trackPhi_efficiency_splitCorrected[iDataset] = (TH1D*)H1D_trackPhi_efficiency[iDataset]->Clone("trackPhi_efficiency_splitCorrected"+Datasets[iDataset]);
    // H1D_trackPhi_efficiency_splitCorrected[iDataset]->Add(H1D_trackPhi_efficiency_split[iDataset],-1.);

    H1D_trackPhi_efficiency_splitCorrected[iDataset] = (TH1D*)H1D_trackPhi_mcparticles[iDataset]->Clone("trackPhi_efficiency_split"+Datasets[iDataset]);
    H1D_trackPhi_efficiency_splitCorrected[iDataset]->Reset("M");
    divideSuccess_split = H1D_trackPhi_efficiency_splitCorrected[iDataset]->Divide(H1D_trackPhi_assoctracks_split[iDataset], H1D_trackPhi_mcparticles[iDataset]);
    H1D_trackPhi_efficiency_splitCorrected[iDataset]->Scale(-1.);
    H1D_trackPhi_efficiency_splitCorrected[iDataset]->Add(H1D_trackPhi_efficiency[iDataset],1.);

    // H1D_trackPhi_efficiency_splitAndSecondaryCorrected[iDataset] = (TH1D*)H1D_trackPhi_mcparticles[iDataset]->Clone("trackPhi_efficiency_splitAndSecondaryCorrected"+Datasets[iDataset]);
    // H1D_trackPhi_efficiency_splitAndSecondaryCorrected[iDataset]->Reset("M");
    // divideSuccess_split = H1D_trackPhi_efficiency_splitAndSecondaryCorrected[iDataset]->Divide(H1D_trackPhi_assoctracks_nonSplitSecondary[iDataset], H1D_trackPhi_mcparticles[iDataset]);
    // H1D_trackPhi_efficiency_splitAndSecondaryCorrected[iDataset]->Scale(-1.);
    // H1D_trackPhi_efficiency_splitAndSecondaryCorrected[iDataset]->Add(H1D_trackPhi_efficiency_splitCorrected[iDataset],1.);
  }

  TString* pdfNameEntriesNorm = new TString("track_Phi_efficiency"+dummyName[0]+"_@pt["+Form("%.2f", ptRange[0])+","+Form("%.1f", ptRange[1])+"]"+"_@eta["+Form("%.1f", etaRange[0])+","+Form("%.1f", etaRange[1])+"]");
  TString* pdfNameEntriesNorm_split = new TString("track_Phi_efficiency"+dummyName[0]+"_@pt["+Form("%.2f", ptRange[0])+","+Form("%.1f", ptRange[1])+"]"+"_@eta["+Form("%.1f", etaRange[0])+","+Form("%.1f", etaRange[1])+"]_split");
  TString* pdfNameEntriesNorm_splitCorrected = new TString("track_Phi_efficiency"+dummyName[0]+"_@pt["+Form("%.2f", ptRange[0])+","+Form("%.1f", ptRange[1])+"]"+"_@eta["+Form("%.1f", etaRange[0])+","+Form("%.1f", etaRange[1])+"]_splitCorrected");
  // TString* pdfNameEntriesNorm_splitAndSecondaryCorrected = new TString("track_Phi_efficiency"+dummyName[0]+"_@pt["+Form("%.2f", ptRange[0])+","+Form("%.1f", ptRange[1])+"]"+"_@eta["+Form("%.1f", etaRange[0])+","+Form("%.1f", etaRange[1])+"]_splitAndSecondaryCorrected");

  // TString textContext(contextDatasetComp(""));
  TString textContext(contextCustomTwoFields(*texDatasetsComparisonCommonDenominator, "#splitline{"+contextPtRange(ptRange)+"}{"+contextEtaRange(etaRange)+"}", ""));

  if (divideSuccess == true) {
    Draw_TH1_Histograms_in_one(H1D_trackPhi_efficiency, DatasetsNames, nDatasets, textContext, pdfNameEntriesNorm, texPhiX, texTrackEfficiency, texCollisionDataInfo, drawnWindowAuto, "efficiency");
  }
  else {
    cout << "Divide failed in Draw_Efficiency_Phi_DatasetComparison" << endl;
  }
  Draw_TH1_Histograms_in_one(H1D_trackPhi_efficiency_splitCorrected, DatasetsNames, nDatasets, textContext, pdfNameEntriesNorm_splitCorrected, texPhiX, texTrackEfficiency, texCollisionDataInfo, drawnWindowAuto, "efficiency");
  // Draw_TH1_Histograms_in_one(H1D_trackPhi_efficiency_splitAndSecondaryCorrected, DatasetsNames, nDatasets, textContext, pdfNameEntriesNorm_splitAndSecondaryCorrected, texPhiX, texTrackEfficiency, texCollisionDataInfo, drawnWindowAuto, "efficiency");
}


void Draw_Efficiency_Eta_PtRangeComparison(float* ptRange, int nPtRanges, int iDataset, const char options[]) {

  TH3D* H3D_trackPtEtaPhi_mcparticles;
  TH3D* H3D_trackPtEtaPhi_assoctracks;

  H3D_trackPtEtaPhi_mcparticles = (TH3D*)((TH3D*)file_O2Analysis_list[iDataset]->Get(analysisWorkflow[iDataset]+"/h3_particle_pt_particle_eta_particle_phi_mcpartofinterest"))->Clone("Draw_Efficiency_Eta_PtRangeComparison_mcparticles"+Datasets[iDataset]);
  H3D_trackPtEtaPhi_assoctracks = (TH3D*)((TH3D*)file_O2Analysis_list[iDataset]->Get(analysisWorkflow[iDataset]+"/h3_particle_pt_particle_eta_particle_phi_associatedtrack_primary"))->Clone("Draw_Efficiency_Eta_PtRangeComparison_associatedtrackSelColl"+Datasets[iDataset]);
  // H3D_jetRjetPtjetEta = (TH3D*)((TH3D*)file_O2Analysis_list[iDataset]->Get(analysisWorkflow[iDataset]+"/h3_jet_r_jet_pt_jet_eta"+jetFinderQaHistType[iJetFinderQaType]))->Clone("Draw_Eta_PtCutComparison"+Datasets[iDataset]);

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
    Draw_TH1_Histograms_in_one(H1D_trackEta_efficiency, PtCutsLegend, nPtRanges, textContext, pdfNameEntriesNorm, texEtaX, texTrackEfficiency, texCollisionDataInfo, drawnWindowAuto, "");
  }
  else {
    cout << "Divide failed in Draw_Efficiency_Eta_DatasetComparison" << endl;
  }
}


void Draw_Efficiency_Phi_PtRangeComparison(float* ptRange, int nPtRanges, float* etaRange, int iDataset, const char options[]) {

  TH3D* H3D_trackPtEtaPhi_mcparticles;
  TH3D* H3D_trackPtEtaPhi_assoctracks;

  H3D_trackPtEtaPhi_mcparticles = (TH3D*)((TH3D*)file_O2Analysis_list[iDataset]->Get(analysisWorkflow[iDataset]+"/h3_particle_pt_particle_eta_particle_phi_mcpartofinterest"))->Clone("Draw_Efficiency_Phi_PtRangeComparison_mcparticles"+Datasets[iDataset]);
  H3D_trackPtEtaPhi_assoctracks = (TH3D*)((TH3D*)file_O2Analysis_list[iDataset]->Get(analysisWorkflow[iDataset]+"/h3_particle_pt_particle_eta_particle_phi_associatedtrack_primary"))->Clone("Draw_Efficiency_Phi_PtRangeComparison_associatedtrackSelColl"+Datasets[iDataset]);
  // H3D_jetRjetPtjetEta = (TH3D*)((TH3D*)file_O2Analysis_list[iDataset]->Get(analysisWorkflow[iDataset]+"/h3_jet_r_jet_pt_jet_eta"+jetFinderQaHistType[iJetFinderQaType]))->Clone("Draw_Phi_PtCutComparison"+Datasets[iDataset]);

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
    Draw_TH1_Histograms_in_one(H1D_trackPhi_efficiency, PtCutsLegend, nPtRanges, textContext, pdfNameEntriesNorm, texPhiX, texTrackEfficiency, texCollisionDataInfo, drawnWindowAuto, "");
  }
  else {
    cout << "Divide failed in Draw_Efficiency_Phi_DatasetComparison" << endl;
  }
}



void Draw_mcPt_DatasetComparison(int particleStatusOption) {
  TH1D* H1D_vx[nDatasets];
  TH1D* H1D_vx_rebinned[nDatasets];

  TString Folder;
  TString Name;
  if (particleStatusOption == 0) {
    Folder = "";
    Name = "";
  }
  if (particleStatusOption == 1) {
    Folder = "Primaries/";
    Name = "Primaries_";
  }
  if (particleStatusOption == 2) {
    Folder = "TransportSecondaries/";
    Name = "TransportSecondaries_";
  }

  for(int iDataset = 0; iDataset < nDatasets; iDataset++){

    H1D_vx[iDataset] = (TH1D*)((TH1D*)file_O2Analysis_list[iDataset]->Get(analysisWorkflow[iDataset]+"/Kine/"+Folder+"ptmc"))->Clone("Draw_mcPt_DatasetComparison"+Datasets[iDataset]);
    H1D_vx_rebinned[iDataset] = (TH1D*)H1D_vx[iDataset]->Rebin(1.,"ptmc_rebinned"+Datasets[iDataset]);

    NormaliseYieldToNEntries(H1D_vx_rebinned[iDataset]);
  }

  TString* pdfName = new TString("track_"+Name+"Pt_DataComp");

  TString textContext("");

  Draw_TH1_Histograms_in_one(H1D_vx_rebinned, DatasetsNames, nDatasets, textContext, pdfName, texPtMC, texTrackPtYield_EntriesNorm, texCollisionDataInfo, drawnWindowAuto, "logx,logy");
}

void Draw_Vx_DatasetComparison(int particleStatusOption) {
  TH1D* H1D_vx[nDatasets];
  TH1D* H1D_vx_rebinned[nDatasets];

  TH1D* dummy = new TH1D("h1_2", "h1 title_2", 5000, -1, 1);

  TString Folder;
  TString Name;
  if (particleStatusOption == 0) {
    Folder = "";
    Name = "";
  }
  if (particleStatusOption == 1) {
    Folder = "Primaries/";
    Name = "Primaries_";
  }
  if (particleStatusOption == 2) {
    Folder = "TransportSecondaries/";
    Name = "TransportSecondaries_";
  }

  for(int iDataset = 0; iDataset < nDatasets; iDataset++){

    H1D_vx[iDataset] = (TH1D*)((TH1D*)file_O2Analysis_list[iDataset]->Get(analysisWorkflow[iDataset]+"/Kine/"+Folder+"vxmc"))->Clone("Draw_Vx_DatasetComparison"+Datasets[iDataset]);
    if (particleStatusOption == 1) {
      std::vector<double> O2H1DYbinsVector = GetTH1Bins(dummy);
      double* O2yBins = &O2H1DYbinsVector[0];
      H1D_vx_rebinned[iDataset] = (TH1D*)H1D_vx[iDataset]->Rebin(5000,"vxmc_rebinned"+Datasets[iDataset], O2yBins);
    }
    else {
      H1D_vx_rebinned[iDataset] = (TH1D*)H1D_vx[iDataset]->Rebin(50.,"vxmc_rebinned"+Datasets[iDataset]);
    }

    NormaliseYieldToNEntries(H1D_vx_rebinned[iDataset]);
  }

  TString* pdfName = new TString("track_"+Name+"Vx_DataComp");

  TString textContext("");

  Draw_TH1_Histograms_in_one(H1D_vx_rebinned, DatasetsNames, nDatasets, textContext, pdfName, texVx, texPartVxYield_EntriesNorm, texCollisionDataInfo, drawnWindowAuto, "logy");
}

void Draw_Vy_DatasetComparison(int particleStatusOption) {
  TH1D* H1D_vy[nDatasets];
  TH1D* H1D_vy_rebinned[nDatasets];

  TH1D* dummy = new TH1D("h1_2", "h1 title_2", 5000, -1, 1);

  TString Folder;
  TString Name;
  if (particleStatusOption == 0) {
    Folder = "";
    Name = "";
  }
  if (particleStatusOption == 1) {
    Folder = "Primaries/";
    Name = "Primaries_";
  }
  if (particleStatusOption == 2) {
    Folder = "TransportSecondaries/";
    Name = "TransportSecondaries_";
  }

  for(int iDataset = 0; iDataset < nDatasets; iDataset++){

    H1D_vy[iDataset] = (TH1D*)((TH1D*)file_O2Analysis_list[iDataset]->Get(analysisWorkflow[iDataset]+"/Kine/"+Folder+"vymc"))->Clone("Draw_Vy_DatasetComparison"+Datasets[iDataset]);
    if (particleStatusOption == 1) {
      std::vector<double> O2H1DYbinsVector = GetTH1Bins(dummy);
      double* O2yBins = &O2H1DYbinsVector[0];
      H1D_vy_rebinned[iDataset] = (TH1D*)H1D_vy[iDataset]->Rebin(5000,"vymc_rebinned"+Datasets[iDataset], O2yBins);
    }
    else {
      H1D_vy_rebinned[iDataset] = (TH1D*)H1D_vy[iDataset]->Rebin(50.,"vymc_rebinned"+Datasets[iDataset]);
    }

    NormaliseYieldToNEntries(H1D_vy_rebinned[iDataset]);
  }

  TString* pdfName = new TString("track_"+Name+"Vy_DataComp");

  TString textContext("");

  Draw_TH1_Histograms_in_one(H1D_vy_rebinned, DatasetsNames, nDatasets, textContext, pdfName, texVy, texPartVyYield_EntriesNorm, texCollisionDataInfo, drawnWindowAuto, "logy");
}

void Draw_Vz_DatasetComparison(int particleStatusOption) {
  TH1D* H1D_vz[nDatasets];
  TH1D* H1D_vz_rebinned[nDatasets];

  TString Folder;
  TString Name;
  if (particleStatusOption == 0) {
    Folder = "";
    Name = "";
  }
  if (particleStatusOption == 1) {
    Folder = "Primaries/";
    Name = "Primaries_";
  }
  if (particleStatusOption == 2) {
    Folder = "TransportSecondaries/";
    Name = "TransportSecondaries_";
  }

  for(int iDataset = 0; iDataset < nDatasets; iDataset++){

    H1D_vz[iDataset] = (TH1D*)((TH1D*)file_O2Analysis_list[iDataset]->Get(analysisWorkflow[iDataset]+"/Kine/"+Folder+"vzmc"))->Clone("Draw_Vz_DatasetComparison"+Datasets[iDataset]);
    H1D_vz_rebinned[iDataset] = (TH1D*)H1D_vz[iDataset]->Rebin(50.,"vzmc_rebinned"+Datasets[iDataset]);

    NormaliseYieldToNEntries(H1D_vz_rebinned[iDataset]);
  }

  TString* pdfName = new TString("track_"+Name+"Vz_DataComp");

  TString textContext("");

  Draw_TH1_Histograms_in_one(H1D_vz_rebinned, DatasetsNames, nDatasets, textContext, pdfName, texVz, texPartVzYield_EntriesNorm, texCollisionDataInfo, drawnWindowAuto, "logy");
}


void Draw_y_iu_DatasetComparison() {
  TH1D* H1D_IUy[nDatasets];
  TH1D* H1D_IUy_rebinned[nDatasets];

  TH1D* dummy = new TH1D("h1", "h1 title", 500, -20, 20);

  for(int iDataset = 0; iDataset < nDatasets; iDataset++){
    H1D_IUy[iDataset] = (TH1D*)((TH1D*)file_O2Analysis_list[iDataset]->Get(analysisWorkflow[iDataset]+"/IU/y"))->Clone("Draw_Y_IU_DatasetComparison"+Datasets[iDataset]);
    // H1D_IUy[iDataset]->GetXaxis()->SetRangeUser(-20,20);
    std::vector<double> O2H1DYbinsVector = GetTH1Bins(dummy);
    double* O2yBins = &O2H1DYbinsVector[0];
    H1D_IUy_rebinned[iDataset] = (TH1D*)H1D_IUy[iDataset]->Rebin(500, "track_y_iu_rebinned"+Datasets[iDataset], O2yBins);
    // H1D_IUy_rebinned[iDataset]->GetXaxis()->SetRangeUser(-20,20);

    // int ibinYlow = H1D_IUy_rebinned[iDataset]->GetXaxis()->FindBin(-20);
    // int ibinYhigh = H1D_IUy_rebinned[iDataset]->GetXaxis()->FindBin(20);

    NormaliseYieldToNEntries(H1D_IUy_rebinned[iDataset]);
  }

  TString* pdfName = new TString("track_IU_y_DataComp");

  TString textContext("");

  Draw_TH1_Histograms_in_one(H1D_IUy_rebinned, DatasetsNames, nDatasets, textContext, pdfName, texIUy, texPartYYield_EntriesNorm, texCollisionDataInfo, drawnWindowAuto, "");
}
