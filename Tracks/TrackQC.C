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
#include "THnSparse.h"
// #include <RooUnfold.h>
// #include "RooUnfoldResponse.h"
// #include "RooUnfoldBayes.h"
// #include "RooUnfoldBinByBin.h"

//My Libraries
#include "./TrackQC_settings.h"
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

// Plot Utilities
TString contextTrackDatasetComp(const char options[]);

//////////// QC plot functions
// Dataset comparison
void Draw_Pt_DatasetComparison();
void Draw_Eta_DatasetComparison();
void Draw_Phi_DatasetComparison();
void Draw_Eta_DatasetComparison_EntriesNorm();

void Draw_Eta_DatasetComparison_trackSelComp();
void Draw_Phi_DatasetComparison_trackSelComp();


void Draw_Pt_CentralityComparison(int iDataset);
void Draw_Eta_CentralityComparison(int iDataset);
void Draw_Phi_CentralityComparison(int iDataset);

void Draw_Pt_Run2Run3Comparison_0010Cent(int iDataset);
void Draw_Eta_Run2Run3Comparison_0010Cent(int iDataset);
void Draw_Phi_Run2Run3Comparison_0010Cent(int iDataset);

void Draw_Phi_DatasetComparison_PtRange(float* PtRange); //works only on modified jetfinderQA for h_track_pt_track_eta_track_phi
void Draw_Eta_DatasetComparison_PtRange(float* PtRange); //works only on modified jetfinderQA for h_track_pt_track_eta_track_phi

void Draw_Sigmapt_vs_pt_DatasetComp();

void Draw_Sigmapt_vs_pt_DatasetComp_legacyTH3();
void Draw_Sigmapt_nonGlobal_uniformTracks_legacyTH3();
void Draw_Sigmapt_nonGlobal_uniformTracks_fromSubtraction_legacyTH3();
void Draw_Sigmapt_nonGlobal_uniformTracks_centralEta_legacyTH3();
void Draw_Sigmapt_vs_pt_DatasetComp_centralEta_legacyTH3();

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


  // Draw_Pt_DatasetComparison();
  // Draw_Eta_DatasetComparison();
  // Draw_Phi_DatasetComparison();

  Draw_Eta_DatasetComparison_trackSelComp();
  Draw_Phi_DatasetComparison_trackSelComp();
  // // Draw_Eta_DatasetComparison_EntriesNorm();


  // for(int iDataset = 0; iDataset < nDatasets; iDataset++){
  //   Draw_Pt_CentralityComparison(iDataset);
  //   Draw_Eta_CentralityComparison(iDataset);
  //   Draw_Phi_CentralityComparison(iDataset);
  // }
  // for(int iDataset = 0; iDataset < nDatasets; iDataset++){
  //   Draw_Pt_Run2Run3Comparison_0010Cent(iDataset);
  //   Draw_Eta_Run2Run3Comparison_0010Cent(iDataset);
  //   Draw_Phi_Run2Run3Comparison_0010Cent(iDataset);
  // }

  // float ptRange1[2] = {0.15, 100};
  // Draw_Phi_DatasetComparison_PtRange(ptRange1); //works only on modified jetfinderQA for h_track_pt_track_eta_track_phi
  // Draw_Eta_DatasetComparison_PtRange(ptRange1); //works only on modified jetfinderQA for h_track_pt_track_eta_track_phi
  // float ptRange2[2] = {2, 100};
  // Draw_Phi_DatasetComparison_PtRange(ptRange2); //works only on modified jetfinderQA for h_track_pt_track_eta_track_phi
  // Draw_Eta_DatasetComparison_PtRange(ptRange2); //works only on modified jetfinderQA for h_track_pt_track_eta_track_phi
  // float ptRange3[2] = {3, 100};
  // Draw_Phi_DatasetComparison_PtRange(ptRange3); //works only on modified jetfinderQA for h_track_pt_track_eta_track_phi
  // Draw_Eta_DatasetComparison_PtRange(ptRange3); //works only on modified jetfinderQA for h_track_pt_track_eta_track_phi
  // float ptRange4[2] = {4, 100};
  // Draw_Phi_DatasetComparison_PtRange(ptRange4); //works only on modified jetfinderQA for h_track_pt_track_eta_track_phi
  // Draw_Eta_DatasetComparison_PtRange(ptRange4); //works only on modified jetfinderQA for h_track_pt_track_eta_track_phi

  Draw_Sigmapt_vs_pt_DatasetComp();
  // Draw_Sigmapt_nonGlobal_uniformTracks();
  // Draw_Sigmapt_nonGlobal_uniformTracks_fromSubtraction();
  // Draw_Sigmapt_nonGlobal_uniformTracks_centralEta();
  // Draw_Sigmapt_vs_pt_DatasetComp_centralEta();
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
  TString texcontextDatasetCompAndRadiusAndVarRange;
  //  if (strstr(options, "track") == NULL) {
  texcontextDatasetCompAndRadiusAndVarRange = *texDatasetsComparisonCommonDenominator;
  return texcontextDatasetCompAndRadiusAndVarRange;
}

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

void Draw_Pt_DatasetComparison() {
  TH2D* H2D_centrality_track[nDatasets];

  TH1D* H1D_trackPt[nDatasets];
  TH1D* H1D_trackPt_rebinned[nDatasets];
  
  TH1D* H1D_trackPt_rebinned_ratios[nDatasets];

  bool divideSuccess = false;

  for(int iDataset = 0; iDataset < nDatasets; iDataset++){

    H2D_centrality_track[iDataset] = (TH2D*)((TH2D*)file_O2Analysis_list[iDataset]->Get(analysisWorkflow[iDataset]+"/h2_centrality_track_pt"))->Clone("Draw_Pt_DatasetComparison"+Datasets[iDataset]);
    H1D_trackPt[iDataset] = (TH1D*)H2D_centrality_track[iDataset]->ProjectionY("trackPt_"+Datasets[iDataset], 1, H2D_centrality_track[iDataset]->GetNbinsX(), "e");


    H1D_trackPt_rebinned[iDataset] = (TH1D*)H1D_trackPt[iDataset]->Rebin(1.,"trackPt_rebinned_"+Datasets[iDataset]);

    // NormaliseYieldToNEntries(H1D_trackPt_rebinned[iDataset]);
    NormaliseYieldToNEvents(H1D_trackPt_rebinned[iDataset], GetNEventsSelected_TrackEffWorkflow(file_O2Analysis_list[iDataset]));

    H1D_trackPt_rebinned_ratios[iDataset] = (TH1D*)H1D_trackPt_rebinned[iDataset]->Clone("trackPt_rebinned_ratios"+Datasets[iDataset]);
    H1D_trackPt_rebinned_ratios[iDataset]->Reset("M");
    divideSuccess = H1D_trackPt_rebinned_ratios[iDataset]->Divide(H1D_trackPt_rebinned[iDataset], H1D_trackPt_rebinned[0]);
  }

  TString* pdfName = new TString("track_Pt_DataComp");
  TString* pdfName_ratio = new TString("track_Pt_DataComp_ratio");

  TString textContext(contextTrackDatasetComp(""));

  Draw_TH1_Histograms_in_one(H1D_trackPt_rebinned, DatasetsNames, nDatasets, textContext, pdfName, texPtX, texTrackPtYield_EventNorm, texCollisionDataInfo, drawnWindowAuto, "logx,logy");
  if (divideSuccess == true) {
    Draw_TH1_Histograms_in_one(H1D_trackPt_rebinned_ratios, DatasetsNames, nDatasets, textContext, pdfName_ratio, texPtX, texRatioDatasets, texCollisionDataInfo, drawnWindowAuto, "autoratio,avoidFirst");
  }
  else {
    cout << "Divide failed in Draw_Pt_DatasetComparison" << endl;
  }
}


void Draw_Eta_DatasetComparison() {
  TH2D* H2D_centrality_track[nDatasets];

  TH1D* H1D_trackEta[nDatasets];
  TH1D* H1D_trackEta_rebinned[nDatasets];
  
  TH1D* H1D_trackEta_rebinned_ratios[nDatasets];

  bool divideSuccess = false;

  for(int iDataset = 0; iDataset < nDatasets; iDataset++){


    H2D_centrality_track[iDataset] = (TH2D*)((TH2D*)file_O2Analysis_list[iDataset]->Get(analysisWorkflow[iDataset]+"/h2_centrality_track_eta"))->Clone("Draw_Eta_DatasetComparison"+Datasets[iDataset]);
    H1D_trackEta[iDataset] = (TH1D*)H2D_centrality_track[iDataset]->ProjectionY("trackEta_"+Datasets[iDataset], 1, H2D_centrality_track[iDataset]->GetNbinsX(), "e");

    H1D_trackEta_rebinned[iDataset] = (TH1D*)H1D_trackEta[iDataset]->Rebin(1.,"trackEta_rebinned"+Datasets[iDataset]);

    // NormaliseYieldToNEntries(H1D_trackEta_rebinned[iDataset]);
    NormaliseYieldToNEvents(H1D_trackEta_rebinned[iDataset], GetNEventsSelected_TrackEffWorkflow(file_O2Analysis_list[iDataset]));

    H1D_trackEta_rebinned_ratios[iDataset] = (TH1D*)H1D_trackEta_rebinned[iDataset]->Clone("trackEta_rebinned_ratios"+Datasets[iDataset]);
    H1D_trackEta_rebinned_ratios[iDataset]->Reset("M");
    divideSuccess = H1D_trackEta_rebinned_ratios[iDataset]->Divide(H1D_trackEta_rebinned[iDataset], H1D_trackEta_rebinned[0]);
  }

  TString* pdfNameEventNorm = new TString("track_Eta_DataComp_EventNorm");
  TString* pdfNameEventNorm_ratio = new TString("track_Eta_DataComp_EventNorm_ratio");

  TString textContext(contextTrackDatasetComp(""));

  Draw_TH1_Histograms_in_one(H1D_trackEta_rebinned, DatasetsNames, nDatasets, textContext, pdfNameEventNorm, texEtaX, texTrackEtaYield_EventNorm, texCollisionDataInfo, drawnWindowAuto, "");
  if (divideSuccess == true) {
    Draw_TH1_Histograms_in_one(H1D_trackEta_rebinned_ratios, DatasetsNames, nDatasets, textContext, pdfNameEventNorm_ratio, texEtaX, texRatioDatasets, texCollisionDataInfo, drawnWindowAuto, "autoratio,avoidFirst");
  }
  else {
    cout << "Divide failed in Draw_Eta_DatasetComparison" << endl;
  }
}

void Draw_Eta_DatasetComparison_EntriesNorm() {
  TH2D* H2D_centrality_track[nDatasets];

  TH1D* H1D_trackEta[nDatasets];
  TH1D* H1D_trackEta_rebinned[nDatasets];
  
  TH1D* H1D_trackEta_rebinned_ratios[nDatasets];

  bool divideSuccess = false;

  for(int iDataset = 0; iDataset < nDatasets; iDataset++){

    H2D_centrality_track[iDataset] = (TH2D*)((TH2D*)file_O2Analysis_list[iDataset]->Get(analysisWorkflow[iDataset]+"/h2_centrality_track_eta"))->Clone("Draw_Eta_DatasetComparison"+Datasets[iDataset]);
    H1D_trackEta[iDataset] = (TH1D*)H2D_centrality_track[iDataset]->ProjectionY("trackEta_"+Datasets[iDataset], 1, H2D_centrality_track[iDataset]->GetNbinsX(), "e");

    H1D_trackEta_rebinned[iDataset] = (TH1D*)H1D_trackEta[iDataset]->Rebin(1.,"trackEta_rebinned_EntriesNorm"+Datasets[iDataset]);

    NormaliseYieldToNEntries(H1D_trackEta_rebinned[iDataset]);
    // H1D_trackEta_rebinned[iDataset]->Scale(1./H1D_trackEta_rebinned[iDataset]->Integral(1, H1D_trackEta_rebinned[iDataset]->GetNbinsX()),"width"); // If option contains "width" the bin contents and errors are divided by the bin width.
    // H1D_trackEta_rebinned[iDataset]->Scale(1./H1D_trackEta_rebinned[iDataset]->Integral(ibincutlow, ibincuthigh),"width"); // If option contains "width" the bin contents and errors are divided by the bin width.

    H1D_trackEta_rebinned_ratios[iDataset] = (TH1D*)H1D_trackEta_rebinned[iDataset]->Clone("trackEta_rebinned_ratios_EntriesNorm"+Datasets[iDataset]);
    H1D_trackEta_rebinned_ratios[iDataset]->Reset("M");
    divideSuccess = H1D_trackEta_rebinned_ratios[iDataset]->Divide(H1D_trackEta_rebinned[iDataset], H1D_trackEta_rebinned[0]);
  }

  TString* pdfNameEntriesNorm = new TString("track_Eta_DataComp_EntriesNorm");
  TString* pdfNameEntriesNorm_ratio = new TString("track_Eta_DataComp_EntriesNorm_ratio");

  TString textContext(contextTrackDatasetComp(""));

  Draw_TH1_Histograms_in_one(H1D_trackEta_rebinned, DatasetsNames, nDatasets, textContext, pdfNameEntriesNorm, texEtaX, texTrackEtaYield_EntriesNorm, texCollisionDataInfo, drawnWindowAuto, "minYnotZero");
  if (divideSuccess == true) {
    Draw_TH1_Histograms_in_one(H1D_trackEta_rebinned_ratios, DatasetsNames, nDatasets, textContext, pdfNameEntriesNorm_ratio, texEtaX, texRatioDatasets, texCollisionDataInfo, drawnWindowAuto, "autoratio,avoidFirst");
  }
  else {
    cout << "Divide failed in Draw_Eta_DatasetComparison" << endl;
  }
}

void Draw_Phi_DatasetComparison() { 
  TH2D* H2D_centrality_track[nDatasets];

  TH1D* H1D_trackPhi[nDatasets];
  TH1D* H1D_trackPhi_rebinned[nDatasets];
  
  TH1D* H1D_trackPhi_rebinned_ratios[nDatasets];

  bool divideSuccess = false;

  for(int iDataset = 0; iDataset < nDatasets; iDataset++){

    H2D_centrality_track[iDataset] = (TH2D*)((TH2D*)file_O2Analysis_list[iDataset]->Get(analysisWorkflow[iDataset]+"/h2_centrality_track_phi"))->Clone("Draw_Phi_DatasetComparison"+Datasets[iDataset]);
    H1D_trackPhi[iDataset] = (TH1D*)H2D_centrality_track[iDataset]->ProjectionY("trackPhi_"+Datasets[iDataset], 1, H2D_centrality_track[iDataset]->GetNbinsX(), "e");

    H1D_trackPhi_rebinned[iDataset] = (TH1D*)H1D_trackPhi[iDataset]->Rebin(1.,"trackPhi_rebinned_"+Datasets[iDataset]);

    // NormaliseYieldToNEntries(H1D_trackPhi_rebinned[iDataset]);
    NormaliseYieldToNEvents(H1D_trackPhi_rebinned[iDataset], GetNEventsSelected_TrackEffWorkflow(file_O2Analysis_list[iDataset]));

    H1D_trackPhi_rebinned_ratios[iDataset] = (TH1D*)H1D_trackPhi_rebinned[iDataset]->Clone("trackPhi_rebinned_ratios"+Datasets[iDataset]);
    H1D_trackPhi_rebinned_ratios[iDataset]->Reset("M");
    divideSuccess = H1D_trackPhi_rebinned_ratios[iDataset]->Divide(H1D_trackPhi_rebinned[iDataset], H1D_trackPhi_rebinned[0]);
  }

  TString* pdfName = new TString("track_Phi_DataComp");
  TString* pdfName_ratio = new TString("track_Phi_DataComp_ratio");

  TString textContext(contextTrackDatasetComp(""));

  Draw_TH1_Histograms_in_one(H1D_trackPhi_rebinned, DatasetsNames, nDatasets, textContext, pdfName, texPhiX, texTrackPhiYield_EventNorm, texCollisionDataInfo, drawnWindowAuto, "");
  if (divideSuccess == true) {
    Draw_TH1_Histograms_in_one(H1D_trackPhi_rebinned_ratios, DatasetsNames, nDatasets, textContext, pdfName_ratio, texPhiX, texRatioDatasets, texCollisionDataInfo, drawnWindowAuto, "autoratio,avoidFirst");
  }
  else {
    cout << "Divide failed in Draw_Phi_DatasetComparison" << endl;
  }
}


void Draw_Eta_DatasetComparison_PtRange(float* PtRange) { //works only on modified jetfinderQA for h_track_pt_track_eta_track_phi

  TH3D* H3D_trackPtEtaPhi[nDatasets];
  TH1D* H1D_trackEta[nDatasets];
  TH1D* H1D_trackEta_rebinned[nDatasets];
  
  TH1D* H1D_trackEta_rebinned_ratios[nDatasets];

  bool divideSuccess = false;

  float PtCutLow =0;
  float PtCutHigh = 0;
  for(int iDataset = 0; iDataset < nDatasets; iDataset++){
    H3D_trackPtEtaPhi[iDataset] = (TH3D*)((TH3D*)file_O2Analysis_list[iDataset]->Get(analysisWorkflow[iDataset]+"/h_track_pt_track_eta_track_phi"))->Clone("Draw_Eta_DatasetComparison_PtRange"+Datasets[iDataset]+Form("%.1f", PtRange[0])+"<pt<"+Form("%.1f", PtRange[1]));

    PtCutLow = PtRange[0];
    PtCutHigh = PtRange[1];
    int ibinPt_low = H3D_trackPtEtaPhi[iDataset]->GetXaxis()->FindBin(PtCutLow);
    int ibinPt_high = H3D_trackPtEtaPhi[iDataset]->GetXaxis()->FindBin(PtCutHigh);
    H3D_trackPtEtaPhi[iDataset]->GetXaxis()->SetRange(ibinPt_low,ibinPt_high);


    H1D_trackEta[iDataset] = (TH1D*)H3D_trackPtEtaPhi[iDataset]->ProjectionY("trackEta_"+Datasets[iDataset]+Form("%.2f", PtCutLow)+"<pt<"+Form("%.1f", PtCutHigh), ibinPt_low, ibinPt_high, 0, -1, "e");
    if (PtCutLow > 3){
      H1D_trackEta_rebinned[iDataset] = (TH1D*)H1D_trackEta[iDataset]->Rebin(2.,"trackEta_rebinned_"+Datasets[iDataset]+Form("%.2f", PtCutLow)+"<pt<"+Form("%.1f", PtCutHigh));
    }
    else {
      H1D_trackEta_rebinned[iDataset] = (TH1D*)H1D_trackEta[iDataset]->Rebin(1.,"trackEta_rebinned_"+Datasets[iDataset]+Form("%.2f", PtCutLow)+"<pt<"+Form("%.1f", PtCutHigh));
    }
    // NormaliseYieldToNEntries(H1D_trackEta_rebinned[iDataset]);
    NormaliseYieldToNEvents(H1D_trackEta_rebinned[iDataset], GetNEventsSelected_TrackEffWorkflow(file_O2Analysis_list[iDataset]));

    H1D_trackEta_rebinned_ratios[iDataset] = (TH1D*)H1D_trackEta_rebinned[iDataset]->Clone("trackEta_rebinned_ratios"+Datasets[iDataset]);
    H1D_trackEta_rebinned_ratios[iDataset]->Reset("M");
    divideSuccess = H1D_trackEta_rebinned_ratios[iDataset]->Divide(H1D_trackEta_rebinned[iDataset], H1D_trackEta_rebinned[0]);
  }

  TString* pdfName = new TString("track_Eta_DataComp"+dummyName[0]+"_@pt["+Form("%.2f", PtCutLow)+","+Form("%.1f", PtCutHigh)+"]");
  TString* pdfName_ratio = new TString("track_Eta_DataComp_ratio"+dummyName[0]+"_@pt["+Form("%.2f", PtCutLow)+","+Form("%.1f", PtCutHigh)+"]");

  TString textContext(contextPtRange(PtRange));

  Draw_TH1_Histograms_in_one(H1D_trackEta_rebinned, DatasetsNames, nDatasets, textContext, pdfName, texEtaX, texTrackEtaYield_EventNorm, texCollisionDataInfo, drawnWindowAuto, "");
  if (divideSuccess == true) {
    Draw_TH1_Histograms_in_one(H1D_trackEta_rebinned_ratios, DatasetsNames, nDatasets, textContext, pdfName_ratio, texEtaX, texRatioDatasets, texCollisionDataInfo, drawnWindowAuto, "autoratio,avoidFirst");
  }
  else {
    cout << "Divide failed in Draw_Eta_DatasetComparison" << endl;
  }
}

void Draw_Phi_DatasetComparison_PtRange(float* PtRange) { //works only on modified jetfinderQA for h_track_pt_track_eta_track_phi

  TH3D* H3D_trackPtEtaPhi[nDatasets];
  TH1D* H1D_trackPhi[nDatasets];
  TH1D* H1D_trackPhi_rebinned[nDatasets];
  
  TH1D* H1D_trackPhi_rebinned_ratios[nDatasets];

  bool divideSuccess = false;

  float PtCutLow =0;
  float PtCutHigh = 0;
  for(int iDataset = 0; iDataset < nDatasets; iDataset++){
    H3D_trackPtEtaPhi[iDataset] = (TH3D*)((TH3D*)file_O2Analysis_list[iDataset]->Get(analysisWorkflow[iDataset]+"/h_track_pt_track_eta_track_phi"))->Clone("Draw_Phi_DatasetComparison_PtRange"+Datasets[iDataset]+Form("%.1f", PtRange[0])+"<pt<"+Form("%.1f", PtRange[1]));

    PtCutLow = PtRange[0];
    PtCutHigh = PtRange[1];
    int ibinPt_low = H3D_trackPtEtaPhi[iDataset]->GetXaxis()->FindBin(PtCutLow);
    int ibinPt_high = H3D_trackPtEtaPhi[iDataset]->GetXaxis()->FindBin(PtCutHigh);
    H3D_trackPtEtaPhi[iDataset]->GetXaxis()->SetRange(ibinPt_low,ibinPt_high);


    H1D_trackPhi[iDataset] = (TH1D*)H3D_trackPtEtaPhi[iDataset]->ProjectionZ("trackPhi_"+Datasets[iDataset]+Form("%.2f", PtCutLow)+"<pt<"+Form("%.1f", PtCutHigh), ibinPt_low, ibinPt_high, 0, -1, "e");
    H1D_trackPhi_rebinned[iDataset] = (TH1D*)H1D_trackPhi[iDataset]->Rebin(1.,"trackPhi_rebinned_"+Datasets[iDataset]+Form("%.2f", PtCutLow)+"<pt<"+Form("%.1f", PtCutHigh));
    

    // NormaliseYieldToNEntries(H1D_trackPhi_rebinned[iDataset]);
    NormaliseYieldToNEvents(H1D_trackPhi_rebinned[iDataset], GetNEventsSelected_TrackEffWorkflow(file_O2Analysis_list[iDataset]));

    H1D_trackPhi_rebinned_ratios[iDataset] = (TH1D*)H1D_trackPhi_rebinned[iDataset]->Clone("trackPhi_rebinned_ratios"+Datasets[iDataset]);
    H1D_trackPhi_rebinned_ratios[iDataset]->Reset("M");
    divideSuccess = H1D_trackPhi_rebinned_ratios[iDataset]->Divide(H1D_trackPhi_rebinned[iDataset], H1D_trackPhi_rebinned[0]);
  }

  TString* pdfName = new TString("track_Phi_DataComp"+dummyName[0]+"_@pt["+Form("%.2f", PtCutLow)+","+Form("%.1f", PtCutHigh)+"]");
  TString* pdfName_ratio = new TString("track_Phi_DataComp_ratio"+dummyName[0]+"_@pt["+Form("%.2f", PtCutLow)+","+Form("%.1f", PtCutHigh)+"]");

  TString textContext(contextPtRange(PtRange));

  Draw_TH1_Histograms_in_one(H1D_trackPhi_rebinned, DatasetsNames, nDatasets, textContext, pdfName, texPhiX, texTrackPhiYield_EventNorm, texCollisionDataInfo, drawnWindowAuto, "");
  if (divideSuccess == true) {
    Draw_TH1_Histograms_in_one(H1D_trackPhi_rebinned_ratios, DatasetsNames, nDatasets, textContext, pdfName_ratio, texPhiX, texRatioDatasets, texCollisionDataInfo, drawnWindowAuto, "autoratio,avoidFirst");
  }
  else {
    cout << "Divide failed in Draw_Phi_DatasetComparison" << endl;
  }
}



void Draw_Pt_CentralityComparison(int iDataset) {

  TH2D* H3D_trackPttrackCent;
  TH1D* H1D_trackPt[nCentralityBins];
  TH1D* H1D_trackPt_rebinned[nCentralityBins];
  
  H3D_trackPttrackCent = (TH2D*)((TH2D*)file_O2Analysis_list[iDataset]->Get(analysisWorkflow[iDataset]+"/h2_centrality_track_pt"))->Clone("Draw_Pt_CentralityComparison"+Datasets[iDataset]);
  H3D_trackPttrackCent->Sumw2();

  int ibinCent_low, ibinCent_high;
  TString CentralityLegend[nCentralityBins];
  std::stringstream ss;
  for(int iCentralityBin = 0; iCentralityBin < nCentralityBins; iCentralityBin++){

    ibinCent_low = H3D_trackPttrackCent->GetXaxis()->FindBin(arrayCentralityBinning[iCentralityBin]);
    ibinCent_high = H3D_trackPttrackCent->GetXaxis()->FindBin(arrayCentralityBinning[iCentralityBin+1])-1;
    H1D_trackPt[iCentralityBin] = (TH1D*)H3D_trackPttrackCent->ProjectionY("trackPt_"+Datasets[iDataset]+"_@cent["+Form("%.1d", ibinCent_low)+","+Form("%.1d", ibinCent_high)+"]", ibinCent_low, ibinCent_high, "e");
    H1D_trackPt_rebinned[iCentralityBin] = (TH1D*)H1D_trackPt[iCentralityBin]->Rebin(1.,"trackPt_rebinned_"+Datasets[iDataset]+"_@cent["+Form("%.1d", ibinCent_low)+","+Form("%.1d", ibinCent_high)+"]");

    NormaliseYieldToNEvents(H1D_trackPt_rebinned[iCentralityBin], GetNEventsSelectedCentrality_JetFramework(file_O2Analysis_list[iDataset], arrayCentralityBinning[iCentralityBin], arrayCentralityBinning[iCentralityBin+1], trainId));

    ss << "Cent " << arrayCentralityBinning[iCentralityBin] << " - " << arrayCentralityBinning[iCentralityBin+1] << " ";
    CentralityLegend[iCentralityBin] = (TString)ss.str();
    ss.str("");
    ss.clear();
  }

  TString* pdfName = new TString("track_CentralityComp_"+Datasets[iDataset]+"_Pt");

  TString textContext(contextCustomOneField(*texDatasetsComparisonCommonDenominator, ""));

  Draw_TH1_Histograms_in_one(H1D_trackPt_rebinned, CentralityLegend, nCentralityBins, textContext, pdfName, texPtX, texTrackPtYield_EventNorm_CentWindow, texCollisionDataInfo, drawnWindowAuto, "logy");
  TString* pdfName2 = new TString("track_CentralityComp_"+Datasets[iDataset]+"_Pt_logx");
  Draw_TH1_Histograms_in_one(H1D_trackPt_rebinned, CentralityLegend, nCentralityBins, textContext, pdfName2, texPtX, texTrackPtYield_EventNorm_CentWindow, texCollisionDataInfo, drawnWindowAuto, "logx,logy");
}



void Draw_Eta_CentralityComparison(int iDataset) {

  TH2D* H3D_trackEtatrackCent;
  TH1D* H1D_trackEta[nCentralityBins];
  TH1D* H1D_trackEta_rebinned[nCentralityBins];
  
  H3D_trackEtatrackCent = (TH2D*)((TH2D*)file_O2Analysis_list[iDataset]->Get(analysisWorkflow[iDataset]+"/h2_centrality_track_eta"))->Clone("Draw_Eta_CentralityComparison"+Datasets[iDataset]);
  H3D_trackEtatrackCent->Sumw2();

  int ibinCent_low, ibinCent_high;
  TString CentralityLegend[nCentralityBins];
  std::stringstream ss;
  for(int iCentralityBin = 0; iCentralityBin < nCentralityBins; iCentralityBin++){

    ibinCent_low = H3D_trackEtatrackCent->GetXaxis()->FindBin(arrayCentralityBinning[iCentralityBin]);
    ibinCent_high = H3D_trackEtatrackCent->GetXaxis()->FindBin(arrayCentralityBinning[iCentralityBin+1])-1;
    H1D_trackEta[iCentralityBin] = (TH1D*)H3D_trackEtatrackCent->ProjectionY("trackEta_"+Datasets[iDataset]+"_@cent["+Form("%.1d", ibinCent_low)+","+Form("%.1d", ibinCent_high)+"]", ibinCent_low, ibinCent_high, "e");
    H1D_trackEta_rebinned[iCentralityBin] = (TH1D*)H1D_trackEta[iCentralityBin]->Rebin(1.,"trackEta_rebinned_"+Datasets[iDataset]+"_@cent["+Form("%.1d", ibinCent_low)+","+Form("%.1d", ibinCent_high)+"]");

    NormaliseYieldToNEvents(H1D_trackEta_rebinned[iCentralityBin], GetNEventsSelectedCentrality_JetFramework(file_O2Analysis_list[iDataset], arrayCentralityBinning[iCentralityBin], arrayCentralityBinning[iCentralityBin+1], trainId));

    ss << "Cent " << arrayCentralityBinning[iCentralityBin] << " - " << arrayCentralityBinning[iCentralityBin+1] << " ";
    CentralityLegend[iCentralityBin] = (TString)ss.str();
    ss.str("");
    ss.clear();
  }

  TString* pdfName = new TString("track_CentralityComp_"+Datasets[iDataset]+"_Eta");

  TString textContext(contextCustomOneField(*texDatasetsComparisonCommonDenominator, ""));

  Draw_TH1_Histograms_in_one(H1D_trackEta_rebinned, CentralityLegend, nCentralityBins, textContext, pdfName, texPtX, texTrackEtaYield_EventNorm_CentWindow, texCollisionDataInfo, drawnWindowAuto, "");
}

void Draw_Phi_CentralityComparison(int iDataset) {

  TH2D* H3D_trackPhitrackCent;
  TH1D* H1D_trackPhi[nCentralityBins];
  TH1D* H1D_trackPhi_rebinned[nCentralityBins];
  
  H3D_trackPhitrackCent = (TH2D*)((TH2D*)file_O2Analysis_list[iDataset]->Get(analysisWorkflow[iDataset]+"/h2_centrality_track_phi"))->Clone("Draw_Phi_CentralityComparison"+Datasets[iDataset]);
  H3D_trackPhitrackCent->Sumw2();

  int ibinCent_low, ibinCent_high;
  TString CentralityLegend[nCentralityBins];
  std::stringstream ss;
  for(int iCentralityBin = 0; iCentralityBin < nCentralityBins; iCentralityBin++){

    ibinCent_low = H3D_trackPhitrackCent->GetXaxis()->FindBin(arrayCentralityBinning[iCentralityBin]);
    ibinCent_high = H3D_trackPhitrackCent->GetXaxis()->FindBin(arrayCentralityBinning[iCentralityBin+1])-1;
    H1D_trackPhi[iCentralityBin] = (TH1D*)H3D_trackPhitrackCent->ProjectionY("trackPhi_"+Datasets[iDataset]+"_@cent["+Form("%.1d", ibinCent_low)+","+Form("%.1d", ibinCent_high)+"]", ibinCent_low, ibinCent_high, "e");
    H1D_trackPhi_rebinned[iCentralityBin] = (TH1D*)H1D_trackPhi[iCentralityBin]->Rebin(1.,"trackPhi_rebinned_"+Datasets[iDataset]+"_@cent["+Form("%.1d", ibinCent_low)+","+Form("%.1d", ibinCent_high)+"]");

    NormaliseYieldToNEvents(H1D_trackPhi_rebinned[iCentralityBin], GetNEventsSelectedCentrality_JetFramework(file_O2Analysis_list[iDataset], arrayCentralityBinning[iCentralityBin], arrayCentralityBinning[iCentralityBin+1], trainId));

    ss << "Cent " << arrayCentralityBinning[iCentralityBin] << " - " << arrayCentralityBinning[iCentralityBin+1] << " ";
    CentralityLegend[iCentralityBin] = (TString)ss.str();
    ss.str("");
    ss.clear();
  }

  TString* pdfName = new TString("track_CentralityComp_"+Datasets[iDataset]+"_Phi");

  TString textContext(contextCustomOneField(*texDatasetsComparisonCommonDenominator, ""));

  Draw_TH1_Histograms_in_one(H1D_trackPhi_rebinned, CentralityLegend, nCentralityBins, textContext, pdfName, texPtX, texTrackPhiYield_EventNorm_CentWindow, texCollisionDataInfo, drawnWindowAuto, "");
}



void Draw_Pt_Run2Run3Comparison_0010Cent(int iDataset) {

  int ibinCent_low, ibinCent_high;

  cout << "o2physics" << endl;
  // O2Physics Run 3

  TH2D* H3D_trackPttrackCent;
  TH1D* H1D_trackPt;
  TH1D* H1D_trackPt_rebinned;

  H3D_trackPttrackCent = (TH2D*)((TH2D*)file_O2Analysis_list[iDataset]->Get(analysisWorkflow[iDataset]+"/h2_centrality_track_pt"))->Clone("Draw_Pt_CentralityComparison"+DatasetsNames[iDataset]);
  H3D_trackPttrackCent->Sumw2();

  ibinCent_low = H3D_trackPttrackCent->GetXaxis()->FindBin((double)0 + GLOBAL_epsilon);
  ibinCent_high = H3D_trackPttrackCent->GetXaxis()->FindBin((double)10 - GLOBAL_epsilon);
  H1D_trackPt = (TH1D*)H3D_trackPttrackCent->ProjectionY("trackPt_"+DatasetsNames[iDataset]+"_@cent["+Form("%.1d", ibinCent_low)+","+Form("%.1d", ibinCent_high)+"]", ibinCent_low, ibinCent_high, "e");
  H1D_trackPt_rebinned = (TH1D*)H1D_trackPt->Rebin(1.,"trackPt_rebinned_"+DatasetsNames[iDataset]+"_@cent["+Form("%.1d", ibinCent_low)+","+Form("%.1d", ibinCent_high)+"]");

  NormaliseYieldToNEvents(H1D_trackPt_rebinned, GetNEventsSelectedCentrality_JetFramework(file_O2Analysis_list[iDataset], 0, 10, trainId));

  cout << "aliPhysics" << endl;
  // AliPhysics Run 2

  TH1D* H1D_run2_trackPt;
  TH1D* H1D_run2_trackPt_rebinned;
  H1D_run2_trackPt = (TH1D*)((TObject*)((TObject*)file_AliAnalysis->Get("JetCore_JetCore_57_2050_04_histos"))->FindObject("tracks"))->FindObject("histTrackPt_0");
  std::vector<double> O2H1DPtbinsVector = GetTH1Bins(H1D_trackPt_rebinned);
  double* O2ptBins = &O2H1DPtbinsVector[0];
  H1D_run2_trackPt_rebinned = (TH1D*)H1D_run2_trackPt->Rebin(H1D_trackPt_rebinned->GetNbinsX(), "H1D_run2_trackPt_rebinned", O2ptBins);

  int nEvents_CentWindow;
  TH1D* H1D_Centrality_Run2= (TH1D*)((TObject*)file_AliAnalysis->Get("JetCore_JetCore_57_2050_04_histos"))->FindObject("fHistCentrality");
  int iBinCent_low = H1D_Centrality_Run2->GetXaxis()->FindBin(0 + GLOBAL_epsilon);
  int iBinCent_high = H1D_Centrality_Run2->GetXaxis()->FindBin(10 - GLOBAL_epsilon);
  nEvents_CentWindow = H1D_Centrality_Run2->Integral(iBinCent_low, iBinCent_high);

  H1D_run2_trackPt_rebinned->Scale(1./nEvents_CentWindow,"width"); // If option contains "width" the bin contents and errors are divided by the bin width.

  cout << "ratio" << endl;
  // ratio

  TH1D* H1D_trackVar_rebinned_ratio;
  bool divideSuccess = false;

  H1D_trackVar_rebinned_ratio = (TH1D*)H1D_trackPt_rebinned->Clone("trackVar_rebinned_ratios"+DatasetsNames[iDataset]);
  H1D_trackVar_rebinned_ratio->Reset("M");
  divideSuccess = H1D_trackVar_rebinned_ratio->Divide(H1D_trackPt_rebinned, H1D_run2_trackPt_rebinned);

  cout << "plotting" << endl;
  // plotting
  TH1D* hist_list[2] = {H1D_trackPt_rebinned, H1D_run2_trackPt_rebinned};
  TString RunCompLegend[2];
  std::stringstream ss;

  ss << "run 3";
  RunCompLegend[0] = (TString)ss.str();
  ss.str("");
  ss.clear();
  ss << "run 2";
  RunCompLegend[1] = (TString)ss.str();
  ss.str("");
  ss.clear();

  TString* pdfName = new TString("track_Run2_vs_Run3_"+DatasetsNames[iDataset]+"_@cent[00,10]_Pt");
  TString* pdfName_ratio = new TString("track_Run2_vs_Run3_"+DatasetsNames[iDataset]+"_@cent[00,10]_Pt_ratio");

  TString textContext(contextCustomOneField(RunCompLegend[0]+" vs "+RunCompLegend[1], ""));

  Draw_TH1_Histograms_in_one(hist_list, RunCompLegend, 2, textContext, pdfName, texPtX, texTrackPtYield_EventNorm_CentWindow, texCollisionDataInfo, drawnWindowAuto, "logy");
  TString* pdfName2 = new TString("track_Run2_vs_Run3_"+DatasetsNames[iDataset]+"_@cent[00,10]_Pt_logx");
  Draw_TH1_Histograms_in_one(hist_list, RunCompLegend, 2, textContext, pdfName2, texPtX, texTrackPtYield_EventNorm_CentWindow, texCollisionDataInfo, drawnWindowAuto, "logx,logy");
  if (divideSuccess == true) {
    Draw_TH1_Histogram(H1D_trackVar_rebinned_ratio, textContext, pdfName_ratio, texPtX, texRatioRun3Run2, texCollisionDataInfo, drawnWindowAuto, "autoratio, logx, ratioLine");
  }
  else {
    cout << "Divide failed in Draw_Pt_Run2Run3Comparison_0010Cent" << endl;
  }
}

void Draw_Eta_Run2Run3Comparison_0010Cent(int iDataset) {

  int ibinCent_low, ibinCent_high;


  // O2Physics Run 3

  TH2D* H3D_trackEtatrackCent;
  TH1D* H1D_trackEta;
  TH1D* H1D_trackEta_rebinned;

  H3D_trackEtatrackCent = (TH2D*)((TH2D*)file_O2Analysis_list[iDataset]->Get(analysisWorkflow[iDataset]+"/h2_centrality_track_eta"))->Clone("Draw_Eta_CentralityComparison"+DatasetsNames[iDataset]);
  H3D_trackEtatrackCent->Sumw2();

  ibinCent_low = H3D_trackEtatrackCent->GetXaxis()->FindBin((double)0 + GLOBAL_epsilon);
  ibinCent_high = H3D_trackEtatrackCent->GetXaxis()->FindBin((double)10 - GLOBAL_epsilon);
  H1D_trackEta = (TH1D*)H3D_trackEtatrackCent->ProjectionY("trackEta_"+DatasetsNames[iDataset]+"_@cent["+Form("%.1d", ibinCent_low)+","+Form("%.1d", ibinCent_high)+"]", ibinCent_low, ibinCent_high, "e");
  H1D_trackEta_rebinned = (TH1D*)H1D_trackEta->Rebin(1.,"trackEta_rebinned_"+DatasetsNames[iDataset]+"_@cent["+Form("%.1d", ibinCent_low)+","+Form("%.1d", ibinCent_high)+"]");

  NormaliseYieldToNEvents(H1D_trackEta_rebinned, GetNEventsSelectedCentrality_JetFramework(file_O2Analysis_list[iDataset], 0, 10, trainId));

  // AliPhysics Run 2

  TH1D* H1D_run2_trackEta;
  TH1D* H1D_run2_trackEta_rebinned;
  H1D_run2_trackEta = (TH1D*)((TObject*)((TObject*)file_AliAnalysis->Get("JetCore_JetCore_57_2050_04_histos"))->FindObject("tracks"))->FindObject("histTrackEta_0");
  std::vector<double> O2H1DEtabinsVector = GetTH1Bins(H1D_trackEta_rebinned);
  double* O2etaBins = &O2H1DEtabinsVector[0];
  H1D_run2_trackEta_rebinned = (TH1D*)H1D_run2_trackEta->Rebin(H1D_trackEta_rebinned->GetNbinsX(), "H1D_run2_trackEta_rebinned", O2etaBins);

  int nEvents_CentWindow;
  TH1D* H1D_Centrality_Run2= (TH1D*)((TObject*)file_AliAnalysis->Get("JetCore_JetCore_57_2050_04_histos"))->FindObject("fHistCentrality");
  int iBinCent_low = H1D_Centrality_Run2->GetXaxis()->FindBin(0 + GLOBAL_epsilon);
  int iBinCent_high = H1D_Centrality_Run2->GetXaxis()->FindBin(10 - GLOBAL_epsilon);
  nEvents_CentWindow = H1D_Centrality_Run2->Integral(iBinCent_low, iBinCent_high);

  H1D_run2_trackEta_rebinned->Scale(1./nEvents_CentWindow,"width"); // If option contains "width" the bin contents and errors are divided by the bin width.

  // plotting
  TH1D* hist_list[2] = {H1D_trackEta_rebinned, H1D_run2_trackEta_rebinned};
  TString RunCompLegend[2];
  std::stringstream ss;

  ss << "run 3";
  RunCompLegend[0] = (TString)ss.str();
  ss.str("");
  ss.clear();
  ss << "run 2";
  RunCompLegend[1] = (TString)ss.str();
  ss.str("");
  ss.clear();



  TString* pdfName = new TString("track_Run2_vs_Run3_"+DatasetsNames[iDataset]+"_@cent[00,10]_Eta");

  TString textContext(contextCustomOneField(RunCompLegend[0]+" vs "+RunCompLegend[1], ""));

  Draw_TH1_Histograms_in_one(hist_list, RunCompLegend, 2, textContext, pdfName, texEtaX, texTrackEtaYield_EventNorm_CentWindow, texCollisionDataInfo, drawnWindowAuto, "");
}

void Draw_Phi_Run2Run3Comparison_0010Cent(int iDataset) {

  int ibinCent_low, ibinCent_high;


  // O2Physics Run 3

  TH2D* H3D_trackPhitrackCent;
  TH1D* H1D_trackPhi;
  TH1D* H1D_trackPhi_rebinned;

  H3D_trackPhitrackCent = (TH2D*)((TH2D*)file_O2Analysis_list[iDataset]->Get(analysisWorkflow[iDataset]+"/h2_centrality_track_phi"))->Clone("Draw_Phi_CentralityComparison"+DatasetsNames[iDataset]);
  H3D_trackPhitrackCent->Sumw2();

  ibinCent_low = H3D_trackPhitrackCent->GetXaxis()->FindBin((double)0 + GLOBAL_epsilon);
  ibinCent_high = H3D_trackPhitrackCent->GetXaxis()->FindBin((double)10 - GLOBAL_epsilon);
  H1D_trackPhi = (TH1D*)H3D_trackPhitrackCent->ProjectionY("trackPhi_"+DatasetsNames[iDataset]+"_@cent["+Form("%.1d", ibinCent_low)+","+Form("%.1d", ibinCent_high)+"]", ibinCent_low, ibinCent_high, "e");
  H1D_trackPhi_rebinned = (TH1D*)H1D_trackPhi->Rebin(1.,"trackPhi_rebinned_"+DatasetsNames[iDataset]+"_@cent["+Form("%.1d", ibinCent_low)+","+Form("%.1d", ibinCent_high)+"]");

  NormaliseYieldToNEvents(H1D_trackPhi_rebinned, GetNEventsSelectedCentrality_JetFramework(file_O2Analysis_list[iDataset], 0, 10, trainId));

  // AliPhysics Run 2

  TH1D* H1D_run2_trackPhi;
  TH1D* H1D_run2_trackPhi_rebinned;
  H1D_run2_trackPhi = (TH1D*)((TObject*)((TObject*)file_AliAnalysis->Get("JetCore_JetCore_57_2050_04_histos"))->FindObject("tracks"))->FindObject("histTrackPhi_0");
  std::vector<double> O2H1DPhibinsVector = GetTH1Bins(H1D_trackPhi_rebinned);
  double* O2phiBins = &O2H1DPhibinsVector[0];
  H1D_run2_trackPhi_rebinned = (TH1D*)H1D_run2_trackPhi->Rebin(H1D_trackPhi_rebinned->GetNbinsX(), "H1D_run2_trackPhi_rebinned", O2phiBins);

  int nEvents_CentWindow;
  TH1D* H1D_Centrality_Run2= (TH1D*)((TObject*)file_AliAnalysis->Get("JetCore_JetCore_57_2050_04_histos"))->FindObject("fHistCentrality");
  int iBinCent_low = H1D_Centrality_Run2->GetXaxis()->FindBin(0 + GLOBAL_epsilon);
  int iBinCent_high = H1D_Centrality_Run2->GetXaxis()->FindBin(10 - GLOBAL_epsilon);
  nEvents_CentWindow = H1D_Centrality_Run2->Integral(iBinCent_low, iBinCent_high);

  H1D_run2_trackPhi_rebinned->Scale(1./nEvents_CentWindow,"width"); // If option contains "width" the bin contents and errors are divided by the bin width.

  // plotting
  TH1D* hist_list[2] = {H1D_trackPhi_rebinned, H1D_run2_trackPhi_rebinned};
  TString RunCompLegend[2];
  std::stringstream ss;

  ss << "run 3";
  RunCompLegend[0] = (TString)ss.str();
  ss.str("");
  ss.clear();
  ss << "run 2";
  RunCompLegend[1] = (TString)ss.str();
  ss.str("");
  ss.clear();



  TString* pdfName = new TString("track_Run2_vs_Run3_"+DatasetsNames[iDataset]+"_@cent[00,10]_Phi");

  TString textContext(contextCustomOneField(RunCompLegend[0]+" vs "+RunCompLegend[1], ""));

  Draw_TH1_Histograms_in_one(hist_list, RunCompLegend, 2, textContext, pdfName, texPhiX, texTrackPhiYield_EventNorm_CentWindow, texCollisionDataInfo, drawnWindowAuto, "");
}




void Draw_Sigmapt_vs_pt_DatasetComp_legacyTH3() {
  TH3D* HDsparse_cent_sigmapt_pt[nDatasets];
  TH3D* H3D_cent_sigmapt_pt[nDatasets];

  TH2D* H2D_sigmapt_pt[nDatasets];
  TH2D* H2D_sigmapt_pt_rebinned[nDatasets];

  TH1D* H1D_sigmapt_pt_X_forMedian[nDatasets];

  TH1D* H1D_sigmapt_pt_mean[nDatasets];
  TH1D* H1D_sigmapt_pt_mean_withProfile[nDatasets];
  TH1D* H1D_sigmapt_pt_median[nDatasets];

  for(int iDataset = 0; iDataset < nDatasets; iDataset++){
    HDsparse_cent_sigmapt_pt[iDataset] = (TH3D*)((TH3D*)file_O2Analysis_list[iDataset]->Get(analysisWorkflow[iDataset]+"/h3_centrality_track_pt_track_sigmapt"))->Clone("Draw_Sigmapt_vs_pt_DatasetComp"+Datasets[iDataset]);
    H2D_sigmapt_pt[iDataset] = (TH2D*)HDsparse_cent_sigmapt_pt[iDataset]->Project3D(dummyName[0]+Form("%d", iDataset)+"_sigmapt_e_zy"); //can't use letter D in this or it seems to replace the histogram in current pad (see documentation of ProjectionX function. Isn't mentioned in project3D sadly)
    H2D_sigmapt_pt[iDataset]->Sumw2();

    // H2D_sigmapt_pt_rebinned[iDataset] = (TH2D*)H2D_sigmapt_pt[iDataset]->RebinX(1.,"H2D_sigmapt_pt_rebinned"+Datasets[iDataset]);
    // H2D_sigmapt_pt_rebinned[iDataset]->GetXaxis()->SetRange(0,H2D_sigmapt_pt_rebinned[iDataset]->GetNbinsX()-10);
    // H2D_rhoCentrality_rebinned[iDataset]->GetYaxis()->SetRange(1, H2D_sigmapt_pt_rebinned[iDataset]->FindLastBinAbove(1, 2)); //(asks for the last bin on the y axis (axis number 2) to have strictly more than 1 entry)

    H1D_sigmapt_pt_mean[iDataset] = (TH1D*)HDsparse_cent_sigmapt_pt[iDataset]->ProjectionY("H1D_sigmapt_pt_mean"+Datasets[iDataset], 0, -1, 0, -1, "e");
    H1D_sigmapt_pt_median[iDataset] = (TH1D*)HDsparse_cent_sigmapt_pt[iDataset]->ProjectionY("H1D_sigmapt_pt_median"+Datasets[iDataset], 0, -1, 0, -1, "e");
    
    H1D_sigmapt_pt_mean_withProfile[iDataset] = (TH1D*)H2D_sigmapt_pt[iDataset]->ProfileX("H1D_sigmapt_pt_mean_withProfile"+Datasets[iDataset], 0, -1, "e");



    // get mean of sigmapt
    for(int iBin = 1; iBin <= H1D_sigmapt_pt_mean[iDataset]->GetNbinsX(); iBin++){
      H2D_sigmapt_pt[iDataset]->GetXaxis()->SetRange(iBin,iBin);

      H1D_sigmapt_pt_mean[iDataset]->SetBinContent(iBin, H2D_sigmapt_pt[iDataset]->GetMean(2));
      H1D_sigmapt_pt_mean[iDataset]->SetBinError(iBin, H2D_sigmapt_pt[iDataset]->GetMeanError(2));
      // use TH1::GetMean(axisNumber 1/2/3), GetMeanError() after using setrange
    }

    // median instead of mean:
    double x, q;
    q = 0.5;
    for(int iBin = 1; iBin <= H1D_sigmapt_pt_median[iDataset]->GetNbinsX(); iBin++){
      // H2D_sigmapt_pt[iDataset]->GetXaxis()->SetRange(iBin,iBin);
      H1D_sigmapt_pt_X_forMedian[iDataset] = (TH1D*)H2D_sigmapt_pt[iDataset]->ProjectionY("H2D_sigmapt_pt_X"+Datasets[iDataset], iBin, iBin, "e");

      H1D_sigmapt_pt_X_forMedian[iDataset]->GetQuantiles(1, &x, &q);
      H1D_sigmapt_pt_median[iDataset]->SetBinContent(iBin, x);
      H1D_sigmapt_pt_median[iDataset]->SetBinError(iBin, 0.0001); // no idea how to get the error on the median calculation
    }

  }

  TString* pdfName = new TString("track_sigmapt_vs_pt_DataComp");
  TString* pdfName_logy = new TString("track_sigmapt_vs_pt_DataComp_logy");

  TString* pdfName_mean = new TString("track_sigmapt_mean_vs_pt_DataComp");
  TString* pdfName_mean_logy = new TString("track_sigmapt_mean_vs_pt_DataComp_logy");

  TString* pdfName_mean_withProfile = new TString("track_sigmapt_mean_vs_pt_DataComp_withProfile");
  TString* pdfName_mean_withProfile_logy = new TString("track_sigmapt_mean_vs_pt_DataComp_withProfile_logy");

  TString* pdfName_median = new TString("track_sigmapt_median_vs_pt_DataComp");
  TString* pdfName_median_logy = new TString("track_sigmapt_median_vs_pt_DataComp_logy");

  // TString textContext(contextDatasetComp(""));
  TString textContext(contextCustomOneField(*texDatasetsComparisonCommonDenominator, ""));

  Draw_TH2_Histograms(H2D_sigmapt_pt, DatasetsNames, nDatasets, textContext, pdfName, texPtX, texSigmaPt, texCollisionDataInfo, drawnWindowAuto, "logx,logz,autoRangeSame"); // ?
  Draw_TH2_Histograms(H2D_sigmapt_pt, DatasetsNames, nDatasets, textContext, pdfName_logy, texPtX, texSigmaPt, texCollisionDataInfo, drawnWindowAuto, "logx,logy,logz,autoRangeSame"); // ?

  Draw_TH1_Histograms_in_one(H1D_sigmapt_pt_mean, DatasetsNames, nDatasets, textContext, pdfName_mean, texPtX, texSigmaPtMean, texCollisionDataInfo, drawnWindowAuto, "logx");
  Draw_TH1_Histograms_in_one(H1D_sigmapt_pt_mean, DatasetsNames, nDatasets, textContext, pdfName_mean_logy, texPtX, texSigmaPtMean, texCollisionDataInfo, drawnWindowAuto, "logx,logy");

  Draw_TH1_Histograms_in_one(H1D_sigmapt_pt_mean_withProfile, DatasetsNames, nDatasets, textContext, pdfName_mean_withProfile, texPtX, texSigmaPtMean, texCollisionDataInfo, drawnWindowAuto, "logx");
  Draw_TH1_Histograms_in_one(H1D_sigmapt_pt_mean_withProfile, DatasetsNames, nDatasets, textContext, pdfName_mean_withProfile_logy, texPtX, texSigmaPtMean, texCollisionDataInfo, drawnWindowAuto, "logx,logy");

  Draw_TH1_Histograms_in_one(H1D_sigmapt_pt_median, DatasetsNames, nDatasets, textContext, pdfName_median, texPtX, texSigmaPtMedian, texCollisionDataInfo, drawnWindowAuto, "logx");
  Draw_TH1_Histograms_in_one(H1D_sigmapt_pt_median, DatasetsNames, nDatasets, textContext, pdfName_median_logy, texPtX, texSigmaPtMedian, texCollisionDataInfo, drawnWindowAuto, "logx,logy");
}

void Draw_Sigmapt_nonGlobal_uniformTracks_legacyTH3() {
  TH3D* HDsparse_cent_sigmapt_pt[nDatasets];
  TH3D* H3D_cent_sigmapt_pt[nDatasets];

  TH2D* H2D_sigmapt_pt[nDatasets];
  TH2D* H2D_sigmapt_pt_rebinned[nDatasets];

  for(int iDataset = 0; iDataset < nDatasets; iDataset++){
    HDsparse_cent_sigmapt_pt[iDataset] = (TH3D*)((TH3D*)file_O2Analysis_list[iDataset]->Get(analysisWorkflow[iDataset]+"/h3_centrality_track_pt_track_sigmapt_nonGlobal_uniformTracks"))->Clone("Draw_Sigmapt_vs_pt_DatasetComp"+Datasets[iDataset]);
    H2D_sigmapt_pt[iDataset] = (TH2D*)HDsparse_cent_sigmapt_pt[iDataset]->Project3D(dummyName[0]+Form("%d", iDataset)+"_sigmapt_e_zy"); //can't use letter D in this or it seems to replace the histogram in current pad (see documentation of ProjectionX function. Isn't mentioned in project3D sadly)
    H2D_sigmapt_pt[iDataset]->Sumw2();

    // H2D_sigmapt_pt_rebinned[iDataset] = (TH2D*)H2D_sigmapt_pt[iDataset]->RebinX(1.,"H2D_sigmapt_pt_rebinned"+Datasets[iDataset]);
    // H2D_sigmapt_pt_rebinned[iDataset]->GetXaxis()->SetRange(0,H2D_sigmapt_pt_rebinned[iDataset]->GetNbinsX()-10);
    // H2D_rhoCentrality_rebinned[iDataset]->GetYaxis()->SetRange(1, H2D_sigmapt_pt_rebinned[iDataset]->FindLastBinAbove(1, 2)); //(asks for the last bin on the y axis (axis number 2) to have strictly more than 1 entry)
  }

  TString* pdfName = new TString("track_sigmapt_nonGlobal_uniformTracks_fromWorkflow");

  // TString textContext(contextDatasetComp(""));
  TString textContext(contextCustomOneField(*texDatasetsComparisonCommonDenominator, ""));

  Draw_TH2_Histograms(H2D_sigmapt_pt, DatasetsNames, nDatasets, textContext, pdfName, texPtX, texSigmaPt, texCollisionDataInfo, drawnWindowAuto, "logx,logz,autoRangeSame"); // ?
}

void Draw_Sigmapt_nonGlobal_uniformTracks_fromSubtraction_legacyTH3() {
  // function hardcoded for the case: dataset1, dataset2 = globalTracks, uniformTracks
  if (nDatasets != 2) {
    cout << "Draw_Sigmapt_nonGlobal_uniformTracks should only be run with nDatasets = 2" << endl;
  }
  TH3D* HDsparse_cent_sigmapt_pt[nDatasets];
  TH3D* H3D_cent_sigmapt_pt[nDatasets];

  TH2D* H2D_sigmapt_pt[nDatasets];
  TH2D* H2D_sigmapt_pt_rebinned[nDatasets];

  TH2D* H2D_sigmapt_pt_nonGlobal_uniformTracks[1];

  for(int iDataset = 0; iDataset < nDatasets; iDataset++){
    HDsparse_cent_sigmapt_pt[iDataset] = (TH3D*)((TH3D*)file_O2Analysis_list[iDataset]->Get(analysisWorkflow[iDataset]+"/h3_centrality_track_pt_track_sigmapt"))->Clone("Draw_Sigmapt_nonGlobal_uniformTracks_subtract"+Datasets[iDataset]);
    H2D_sigmapt_pt[iDataset] = (TH2D*)HDsparse_cent_sigmapt_pt[iDataset]->Project3D(dummyName[0]+Form("%d", iDataset)+"_sigmapt_e_zy"); //can't use letter D in this or it seems to replace the histogram in current pad (see documentation of ProjectionX function. Isn't mentioned in project3D sadly)
    H2D_sigmapt_pt[iDataset]->Sumw2();

    // H2D_sigmapt_pt_rebinned[iDataset] = (TH2D*)H2D_sigmapt_pt[iDataset]->RebinX(1.,"H2D_sigmapt_pt_rebinned"+Datasets[iDataset]);
    // H2D_sigmapt_pt_rebinned[iDataset]->GetXaxis()->SetRange(0,H2D_sigmapt_pt_rebinned[iDataset]->GetNbinsX()-10);
    // H2D_rhoCentrality_rebinned[iDataset]->GetYaxis()->SetRange(1, H2D_sigmapt_pt_rebinned[iDataset]->FindLastBinAbove(1, 2)); //(asks for the last bin on the y axis (axis number 2) to have strictly more than 1 entry)
  }

  H2D_sigmapt_pt_nonGlobal_uniformTracks[0] = (TH2D*)H2D_sigmapt_pt[1]->Clone("H2D_sigmapt_pt_nonGlobal_uniformTracks");
  // H2D_sigmapt_pt_nonGlobal_uniformTracks[0]->Scale(-1.);
  H2D_sigmapt_pt_nonGlobal_uniformTracks[0]->Add(H2D_sigmapt_pt[0],-1.);

  TString* pdfName = new TString("track_sigmapt_nonGlobal_uniformTracks_fromSubtraction");

  // TString textContext(contextDatasetComp(""));
  TString textContext(contextCustomOneField(*texDatasetsComparisonCommonDenominator, ""));
  TString DatasetsName[1] = {"nonGlobal uniformTracks"};

  Draw_TH2_Histograms(H2D_sigmapt_pt_nonGlobal_uniformTracks, DatasetsName, 1, textContext, pdfName, texPtX, texSigmaPt, texCollisionDataInfo, drawnWindowAuto, "logx, logz,autoRangeSame"); // ?
}


void Draw_Sigmapt_nonGlobal_uniformTracks_centralEta_legacyTH3() {
  TH3D* HDsparse_cent_sigmapt_pt[nDatasets];
  TH3D* H3D_cent_sigmapt_pt[nDatasets];

  TH2D* H2D_sigmapt_pt[nDatasets];
  TH2D* H2D_sigmapt_pt_rebinned[nDatasets];

  int ibinEta_low, ibinEta_high;
  for(int iDataset = 0; iDataset < nDatasets; iDataset++){
    HDsparse_cent_sigmapt_pt[iDataset] = (TH3D*)((TH3D*)file_O2Analysis_list[iDataset]->Get(analysisWorkflow[iDataset]+"/h3_track_eta_track_pt_track_sigmapt_nonGlobal_uniformTracks"))->Clone("Draw_Sigmapt_nonGlobal_uniformTracks_centralEta"+Datasets[iDataset]);

    ibinEta_low = HDsparse_cent_sigmapt_pt[iDataset]->GetXaxis()->FindBin((double)(-0.1) + GLOBAL_epsilon);
    ibinEta_high = HDsparse_cent_sigmapt_pt[iDataset]->GetXaxis()->FindBin((double)0.1 - GLOBAL_epsilon);
    HDsparse_cent_sigmapt_pt[iDataset]->GetXaxis()->SetRange(ibinEta_low, ibinEta_high);
    H2D_sigmapt_pt[iDataset] = (TH2D*)HDsparse_cent_sigmapt_pt[iDataset]->Project3D(dummyName[0]+Form("%d", iDataset)+"_sigmapt_e_zy"); //can't use letter D in this or it seems to replace the histogram in current pad (see documentation of ProjectionX function. Isn't mentioned in project3D sadly)
    H2D_sigmapt_pt[iDataset]->Sumw2();

    // H2D_sigmapt_pt_rebinned[iDataset] = (TH2D*)H2D_sigmapt_pt[iDataset]->RebinX(1.,"H2D_sigmapt_pt_rebinned"+Datasets[iDataset]);
    // H2D_sigmapt_pt_rebinned[iDataset]->GetXaxis()->SetRange(0,H2D_sigmapt_pt_rebinned[iDataset]->GetNbinsX()-10);
    // H2D_rhoCentrality_rebinned[iDataset]->GetYaxis()->SetRange(1, H2D_sigmapt_pt_rebinned[iDataset]->FindLastBinAbove(1, 2)); //(asks for the last bin on the y axis (axis number 2) to have strictly more than 1 entry)
  }

  TString* pdfName = new TString("track_sigmapt_nonGlobal_uniformTracks_@eta[-0.1;0.1]");

  // TString textContext(contextDatasetComp(""));
  TString textContext(contextCustomOneField(*texDatasetsComparisonCommonDenominator, ""));

  Draw_TH2_Histograms(H2D_sigmapt_pt, DatasetsNames, nDatasets, textContext, pdfName, texPtX, texSigmaPt, texCollisionDataInfo, drawnWindowAuto, "logx,logz,autoRangeSame"); // ?
}

void Draw_Sigmapt_vs_pt_DatasetComp_centralEta_legacyTH3() {
  TH3D* HDsparse_cent_sigmapt_pt[nDatasets];
  TH3D* H3D_cent_sigmapt_pt[nDatasets];

  TH2D* H2D_sigmapt_pt[nDatasets];
  TH2D* H2D_sigmapt_pt_rebinned[nDatasets];

  int ibinEta_low, ibinEta_high;
  for(int iDataset = 0; iDataset < nDatasets; iDataset++){
    HDsparse_cent_sigmapt_pt[iDataset] = (TH3D*)((TH3D*)file_O2Analysis_list[iDataset]->Get(analysisWorkflow[iDataset]+"/h3_track_eta_track_pt_track_sigmapt"))->Clone("Draw_Sigmapt_vs_pt_DatasetComp_centralEta"+Datasets[iDataset]);

    ibinEta_low = HDsparse_cent_sigmapt_pt[iDataset]->GetXaxis()->FindBin((double)(-0.1) + GLOBAL_epsilon);
    ibinEta_high = HDsparse_cent_sigmapt_pt[iDataset]->GetXaxis()->FindBin((double)0.1 - GLOBAL_epsilon);
    HDsparse_cent_sigmapt_pt[iDataset]->GetXaxis()->SetRange(ibinEta_low, ibinEta_high);
    H2D_sigmapt_pt[iDataset] = (TH2D*)HDsparse_cent_sigmapt_pt[iDataset]->Project3D(dummyName[0]+Form("%d", iDataset)+"_sigmapt_e_zy"); //can't use letter D in this or it seems to replace the histogram in current pad (see documentation of ProjectionX function. Isn't mentioned in project3D sadly)
    H2D_sigmapt_pt[iDataset]->Sumw2();
  }

  TString* pdfName = new TString("track_sigmapt_vs_pt_DataComp_@eta[-0.1;0.1]");

  // TString textContext(contextDatasetComp(""));
  TString textContext(contextCustomOneField(*texDatasetsComparisonCommonDenominator, ""));

  Draw_TH2_Histograms(H2D_sigmapt_pt, DatasetsNames, nDatasets, textContext, pdfName, texPtX, texSigmaPt, texCollisionDataInfo, drawnWindowAuto, "logx,logz,autoRangeSame"); // ?
}



void Draw_Sigmapt_vs_pt_DatasetComp() {
  TH3D* HDsparse_cent_sigmapt_pt[nDatasets];
  TH3D* H3D_cent_sigmapt_pt[nDatasets];

  TH2D* H2D_sigmapt_pt_concatenated[nDatasets];
  TH2D* H2D_sigmapt_pt[nDatasets];
  TH2D* H2D_sigmapt_pt_high[nDatasets];

  TH1D* H1D_sigmapt_pt_X_forMedian[nDatasets];

  TH1D* H1D_sigmapt_pt_mean_withProfile[nDatasets];
  TH1D* H1D_sigmapt_pt_high_mean_withProfile[nDatasets];
  TH1D* H1D_sigmapt_pt_median[nDatasets];

  TH1D* H1D_sigmapt_pt_mean_withProfile_concatenated[nDatasets];

  double rebin_xMin = 0;
  double rebin_xMiddle = 10;
  double rebin_xMax = 100;
  int rebin_nLeft = 50;
  int rebin_nRight = 10;
  int nbinsX;
  for(int iDataset = 0; iDataset < nDatasets; iDataset++){
    H2D_sigmapt_pt[iDataset] = (TH2D*)((TH2D*)file_O2Analysis_list[iDataset]->Get(analysisWorkflow[iDataset]+"/h2_track_pt_track_sigmapt"))->Clone("Draw_Sigmapt_vs_pt_DatasetComp_left"+Datasets[iDataset]);
    H2D_sigmapt_pt_high[iDataset] = (TH2D*)((TH2D*)file_O2Analysis_list[iDataset]->Get(analysisWorkflow[iDataset]+"/h2_track_pt_high_track_sigmapt"))->Clone("Draw_Sigmapt_vs_pt_high_DatasetComp_right"+Datasets[iDataset]);
    // H2D_sigmapt_pt[iDataset] = (TH2D*)HDsparse_cent_sigmapt_pt[iDataset]->Project3D(dummyName[0]+Form("%d", iDataset)+"_sigmapt_e_zy"); //can't use letter D in this or it seems to replace the histogram in current pad (see documentation of ProjectionX function. Isn't mentioned in project3D sadly)
    // H2D_sigmapt_pt[iDataset]->Sumw2();

    // H2D_sigmapt_pt_rebinned_simple[iDataset] = (TH2D*)H2D_sigmapt_pt[iDataset]->RebinX(1.,"H2D_sigmapt_pt_rebinnedSimple"+Datasets[iDataset]);
    // H2D_sigmapt_pt_rebinned[iDataset]->GetXaxis()->SetRange(0,H2D_sigmapt_pt_rebinned[iDataset]->GetNbinsX()-10);
    // H2D_rhoCentrality_rebinned[iDataset]->GetYaxis()->SetRange(1, H2D_sigmapt_pt_rebinned[iDataset]->FindLastBinAbove(1, 2)); //(asks for the last bin on the y axis (axis number 2) to have strictly more than 1 entry)

    // H1D_sigmapt_pt_median[iDataset] = (TH1D*)H2D_sigmapt_pt[iDataset]->ProjectionX("H1D_sigmapt_pt_median"+Datasets[iDataset], 0, -1, "e");
    // H1D_sigmapt_pt_median[iDataset]->Reset("M");
    // cout << "testC" << endl;
 
    // x-axis
    std::vector<double> xbinsVectorLeft = GetTH1Bins((TH1D*)H2D_sigmapt_pt[iDataset]->ProjectionX("H1D_sigmapt_pt_left"+Datasets[iDataset], 0, -1, "e"));
    std::vector<double> xbinsVectorRight = GetTH1Bins((TH1D*)H2D_sigmapt_pt_high[iDataset]->ProjectionX("H1D_sigmapt_pt_right"+Datasets[iDataset], 0, -1, "e"));
    xbinsVectorRight.erase(xbinsVectorRight.begin());
    cout << "xbinsVectorLeft.front() = " << xbinsVectorLeft.front() << ", xbinsVectorLeft.back() = " << xbinsVectorLeft.back() << ", xbinsVectorRight.front() = " << xbinsVectorRight.front() << ", xbinsVectorRight.back() = " << xbinsVectorRight.back() << endl;
    std::vector<double> xbinsVectorCombination = xbinsVectorLeft;
    xbinsVectorCombination.insert( xbinsVectorCombination.end(), xbinsVectorRight.begin(), xbinsVectorRight.end() );
    double* xbins_new = &xbinsVectorCombination[0];
    cout << "xbinsVectorCombination.size() = " << xbinsVectorCombination.size() << endl;



    H1D_sigmapt_pt_mean_withProfile[iDataset] = (TH1D*)H2D_sigmapt_pt[iDataset]->ProfileX("H1D_sigmapt_pt_rebinned_mean_withProfile"+Datasets[iDataset], 0, -1, "e");
    H1D_sigmapt_pt_high_mean_withProfile[iDataset] = (TH1D*)H2D_sigmapt_pt_high[iDataset]->ProfileX("H1D_sigmapt_pt_high_rebinned_mean_withProfile"+Datasets[iDataset], 0, -1, "e");

    TH1D H1D_sigmapt_pt_mean_withProfile_concatenated_temp("H1D_sigmapt_pt_mean_withProfile_concatenated_temp"+Datasets[iDataset], "H1D_sigmapt_pt_mean_withProfile_concatenated_temp"+Datasets[iDataset], xbinsVectorCombination.size()-1, xbins_new);
    H1D_sigmapt_pt_mean_withProfile_concatenated[iDataset] = (TH1D*)H1D_sigmapt_pt_mean_withProfile_concatenated_temp.Clone("H1D_sigmapt_pt_mean_withProfile_concatenated"+Datasets[iDataset]);

    // y-axis
    std::vector<double> ybinsVector = GetTH1Bins((TH1D*)H2D_sigmapt_pt[iDataset]->ProjectionY("H1D_sigmapt_pt_y"+Datasets[iDataset], 0, -1, "e"));
    double* ybins_new = &ybinsVector[0];

    TH2D H2D_sigmapt_pt_concatenated_temp("H2D_sigmapt_pt_concatenated_temp"+Datasets[iDataset], "H2D_sigmapt_pt_concatenated_temp"+Datasets[iDataset], xbinsVectorCombination.size()-1, xbins_new, ybinsVector.size()-1, ybins_new);
    H2D_sigmapt_pt_concatenated[iDataset] = (TH2D*)H2D_sigmapt_pt_concatenated_temp.Clone("H2D_sigmapt_pt_concatenated"+Datasets[iDataset]);


    double binContent, binError;
    int iBinLowY_OldHist, iBinHighY_OldHist;
    for(int iBinX = 1; iBinX <= H2D_sigmapt_pt[iDataset]->GetNbinsX(); iBinX++){
      // H1D_sigmapt_pt_mean_withProfile_concatenated[iDataset]->SetBinContent(iBinX, H1D_sigmapt_pt_mean_withProfile[iDataset]->GetBinContent(iBinX));
      // H1D_sigmapt_pt_mean_withProfile_concatenated[iDataset]->SetBinError(iBinX, H1D_sigmapt_pt_mean_withProfile[iDataset]->GetBinError(iBinX));
      for(int iBinY = 1; iBinY <= H2D_sigmapt_pt[iDataset]->GetNbinsY(); iBinY++){
        H2D_sigmapt_pt_concatenated[iDataset]->SetBinContent(iBinX, iBinY, H2D_sigmapt_pt[iDataset]->GetBinContent(iBinX, iBinY));
        H2D_sigmapt_pt_concatenated[iDataset]->SetBinError(iBinX, iBinY, H2D_sigmapt_pt[iDataset]->GetBinError(iBinX, iBinY));
        if (iDataset ==3){
          H2D_sigmapt_pt_concatenated[iDataset]->SetBinContent(iBinX, iBinY, H2D_sigmapt_pt[iDataset]->GetBinContent(iBinX, iBinY) - H2D_sigmapt_pt[0]->GetBinContent(iBinX, iBinY));
          H2D_sigmapt_pt_concatenated[iDataset]->SetBinError(iBinX, iBinY, sqrt(H2D_sigmapt_pt_concatenated[iDataset]->GetBinContent(iBinX, iBinY)));
        } 
      }
    }

    for(int iBinX = 1; iBinX <= H2D_sigmapt_pt_high[iDataset]->GetNbinsX(); iBinX++){
      // H1D_sigmapt_pt_mean_withProfile_concatenated[iDataset]->SetBinContent(H2D_sigmapt_pt[iDataset]->GetNbinsX()+iBinX, H1D_sigmapt_pt_high_mean_withProfile[iDataset]->GetBinContent(iBinX));
      // H1D_sigmapt_pt_mean_withProfile_concatenated[iDataset]->SetBinError(H2D_sigmapt_pt[iDataset]->GetNbinsX()+iBinX, H1D_sigmapt_pt_high_mean_withProfile[iDataset]->GetBinError(iBinX));
      for(int iBinY = 1; iBinY <= H2D_sigmapt_pt[iDataset]->GetNbinsY(); iBinY++){
        H2D_sigmapt_pt_concatenated[iDataset]->SetBinContent(H2D_sigmapt_pt[iDataset]->GetNbinsX()+iBinX, iBinY, H2D_sigmapt_pt_high[iDataset]->GetBinContent(iBinX, iBinY));
        H2D_sigmapt_pt_concatenated[iDataset]->SetBinError(H2D_sigmapt_pt[iDataset]->GetNbinsX()+iBinX, iBinY, H2D_sigmapt_pt_high[iDataset]->GetBinError(iBinX, iBinY));
        if (iDataset ==3){
          H2D_sigmapt_pt_concatenated[iDataset]->SetBinContent(H2D_sigmapt_pt[iDataset]->GetNbinsX()+iBinX, iBinY, H2D_sigmapt_pt_high[iDataset]->GetBinContent(iBinX, iBinY) - H2D_sigmapt_pt_high[0]->GetBinContent(iBinX, iBinY));
          H2D_sigmapt_pt_concatenated[iDataset]->SetBinError(H2D_sigmapt_pt[iDataset]->GetNbinsX()+iBinX, iBinY, sqrt(H2D_sigmapt_pt_concatenated[iDataset]->GetBinContent(H2D_sigmapt_pt[iDataset]->GetNbinsX()+iBinX, iBinY)));
        } 
      }
    }

    H1D_sigmapt_pt_mean_withProfile_concatenated[iDataset] = (TH1D*)H2D_sigmapt_pt_concatenated[iDataset]->ProfileX("H1D_sigmapt_pt_concatenated_rebinned_mean_withProfile"+Datasets[iDataset], 0, -1, "e");



    // // median instead of mean:
    // double x, q;
    // q = 0.5;
    // for(int iBin = 1; iBin <= H1D_sigmapt_pt_median[iDataset]->GetNbinsX(); iBin++){
    //   // H2D_sigmapt_pt[iDataset]->GetXaxis()->SetRange(iBin,iBin);
    //   H1D_sigmapt_pt_X_forMedian[iDataset] = (TH1D*)H2D_sigmapt_pt[iDataset]->ProjectionY("H2D_sigmapt_pt_X"+Datasets[iDataset], iBin, iBin, "e");

    //   H1D_sigmapt_pt_X_forMedian[iDataset]->GetQuantiles(1, &x, &q);
    //   H1D_sigmapt_pt_median[iDataset]->SetBinContent(iBin, x);
    //   H1D_sigmapt_pt_median[iDataset]->SetBinError(iBin, 0.0001); // no idea how to get the error on the median calculation
    // }
  }

  TString* pdfName = new TString("track_sigmapt_vs_pt_DataComp");
  TString* pdfName_logy = new TString("track_sigmapt_vs_pt_DataComp_logy");

  TString* pdfName_mean_withProfile = new TString("track_sigmapt_mean_vs_pt_DataComp_withProfile");
  TString* pdfName_mean_withProfile_logy = new TString("track_sigmapt_mean_vs_pt_DataComp_withProfile_logy");

  // TString* pdfName_median = new TString("track_sigmapt_median_vs_pt_DataComp");
  // TString* pdfName_median_logy = new TString("track_sigmapt_median_vs_pt_DataComp_logy");

  // TString textContext(contextDatasetComp(""));
  TString textContext(contextCustomOneField(*texDatasetsComparisonCommonDenominator, ""));

  std::array<std::array<float, 2>, 2> drawnWindowSigma = {{{0.1, 100}, {0.001, 100}}}; // {{xmin, xmax}, {ymin, ymax}}
  // Draw_TH2_Histograms(H2D_sigmapt_pt_concatenated, DatasetsNames, nDatasets, textContext, pdfName, texPtX, texSigmaPt, texCollisionDataInfo, drawnWindowSigma, "logx,logz,autoRangeSame"); // ?
  Draw_TH2_Histograms(H2D_sigmapt_pt_concatenated, DatasetsNames, nDatasets, textContext, pdfName_logy, texPtX, texSigmaPt, texCollisionDataInfo, drawnWindowSigma, "logx,logy,logz,autoRangeSame"); // ?

  // Draw_TH1_Histograms_in_one(H1D_sigmapt_pt_mean_withProfile_concatenated, DatasetsNames, nDatasets, textContext, pdfName_mean_withProfile, texPtX, texSigmaPtMean, texCollisionDataInfo, drawnWindowAuto, "logx");
  Draw_TH1_Histograms_in_one(H1D_sigmapt_pt_mean_withProfile_concatenated, DatasetsNames, nDatasets, textContext, pdfName_mean_withProfile_logy, texPtX, texSigmaPtMean, texCollisionDataInfo, drawnWindowAuto, "logx,logy");

  // Draw_TH1_Histograms_in_one(H1D_sigmapt_pt_median, DatasetsNames, nDatasets, textContext, pdfName_median, texPtX, texSigmaPtMedian, texCollisionDataInfo, drawnWindowAuto, "logx");
  // Draw_TH1_Histograms_in_one(H1D_sigmapt_pt_median, DatasetsNames, nDatasets, textContext, pdfName_median_logy, texPtX, texSigmaPtMedian, texCollisionDataInfo, drawnWindowAuto, "logx,logy");
}




void Draw_Eta_DatasetComparison_trackSelComp() {
  TH2D* H2D_centrality_track[nDatasets];
  TH1D* H1D_trackEta[nDatasets];
  TH1D* H1D_trackEta_rebinned[nDatasets];
  TH1D* H1D_trackEta_rebinned_ratios[nDatasets];

  TH1D* H1D_trackEta_cutVariationComp[nDatasets];

  bool divideSuccess = false;

  int iGlobal = 0;
  int iITSOnly = 2;
  int iGlobalWithBadTPCCrossedRows = 3;
  int iUniform = 1;

  H2D_centrality_track[iGlobal] = (TH2D*)((TH2D*)file_O2Analysis_list[iGlobal]->Get(analysisWorkflow[iGlobal]+"/h2_centrality_track_eta"))->Clone("Draw_Eta_DatasetComparison"+Datasets[iGlobal]);
  H1D_trackEta[iGlobal] = (TH1D*)H2D_centrality_track[iGlobal]->ProjectionY("trackEta_"+Datasets[iGlobal], 1, H2D_centrality_track[iGlobal]->GetNbinsX(), "e");
  H1D_trackEta_rebinned[iGlobal] = (TH1D*)H1D_trackEta[iGlobal]->Rebin(1.,"trackEta_rebinned"+Datasets[iGlobal]);
  // NormaliseYieldToNEvents(H1D_trackEta_rebinned[iGlobal], GetNEventsSelected_TrackEffWorkflow(file_O2Analysis_list[iGlobal]));
  // NormaliseYieldToNEntries(H1D_trackEta_rebinned[iGlobal]);

  H2D_centrality_track[iITSOnly] = (TH2D*)((TH2D*)file_O2Analysis_list[iITSOnly]->Get(analysisWorkflow[iITSOnly]+"/h2_centrality_track_eta"))->Clone("Draw_Eta_DatasetComparison"+Datasets[iITSOnly]);
  H1D_trackEta[iITSOnly] = (TH1D*)H2D_centrality_track[iITSOnly]->ProjectionY("trackEta_"+Datasets[iITSOnly], 1, H2D_centrality_track[iITSOnly]->GetNbinsX(), "e");
  H1D_trackEta_rebinned[iITSOnly] = (TH1D*)H1D_trackEta[iITSOnly]->Rebin(1.,"trackEta_rebinned"+Datasets[iITSOnly]);
  // NormaliseYieldToNEvents(H1D_trackEta_rebinned[iITSOnly], GetNEventsSelected_TrackEffWorkflow(file_O2Analysis_list[iITSOnly]));
  // NormaliseYieldToNEntries(H1D_trackEta_rebinned[iITSOnly]);

  H2D_centrality_track[iGlobalWithBadTPCCrossedRows] = (TH2D*)((TH2D*)file_O2Analysis_list[iGlobalWithBadTPCCrossedRows]->Get(analysisWorkflow[iGlobalWithBadTPCCrossedRows]+"/h2_centrality_track_eta"))->Clone("Draw_Eta_DatasetComparison"+Datasets[iGlobalWithBadTPCCrossedRows]);
  H1D_trackEta[iGlobalWithBadTPCCrossedRows] = (TH1D*)H2D_centrality_track[iGlobalWithBadTPCCrossedRows]->ProjectionY("trackEta_"+Datasets[iGlobalWithBadTPCCrossedRows], 1, H2D_centrality_track[iGlobalWithBadTPCCrossedRows]->GetNbinsX(), "e");
  H1D_trackEta_rebinned[iGlobalWithBadTPCCrossedRows] = (TH1D*)H1D_trackEta[iGlobalWithBadTPCCrossedRows]->Rebin(1.,"trackEta_rebinned"+Datasets[iGlobalWithBadTPCCrossedRows]);
  H1D_trackEta_rebinned[iGlobalWithBadTPCCrossedRows]->Add(H1D_trackEta_rebinned[iGlobal], -1.);
  // NormaliseYieldToNEvents(H1D_trackEta_rebinned[iGlobalWithBadTPCCrossedRows], GetNEventsSelected_TrackEffWorkflow(file_O2Analysis_list[iGlobalWithBadTPCCrossedRows]));
  // NormaliseYieldToNEntries(H1D_trackEta_rebinned[iGlobalWithBadTPCCrossedRows]);

  H2D_centrality_track[iUniform] = (TH2D*)((TH2D*)file_O2Analysis_list[iUniform]->Get(analysisWorkflow[iUniform]+"/h2_centrality_track_eta"))->Clone("Draw_Eta_DatasetComparison"+Datasets[iUniform]);
  H1D_trackEta[iUniform] = (TH1D*)H2D_centrality_track[iUniform]->ProjectionY("trackEta_"+Datasets[iUniform], 1, H2D_centrality_track[iUniform]->GetNbinsX(), "e");
  H1D_trackEta_rebinned[iUniform] = (TH1D*)H1D_trackEta[iUniform]->Rebin(1.,"trackEta_rebinned"+Datasets[iUniform]);
  // NormaliseYieldToNEvents(H1D_trackEta_rebinned[iUniform], GetNEventsSelected_TrackEffWorkflow(file_O2Analysis_list[iUniform]));
  // NormaliseYieldToNEntries(H1D_trackEta_rebinned[iUniform]);

  H1D_trackEta_cutVariationComp[0] = (TH1D*)H1D_trackEta_rebinned[iGlobal]->Clone("H1D_trackEta_interm1_global_inArray");
  NormaliseYieldToIntegral(H1D_trackEta_cutVariationComp[0]);

  TH1D* H1D_trackEta_interm1_globalAndItsOnly = (TH1D*)H1D_trackEta_rebinned[iGlobal]->Clone("H1D_trackEta_interm1_globalAndItsOnly");
  H1D_trackEta_interm1_globalAndItsOnly->Add(H1D_trackEta_rebinned[iITSOnly]);
  H1D_trackEta_cutVariationComp[1] = (TH1D*)H1D_trackEta_interm1_globalAndItsOnly->Clone("H1D_trackEta_interm1_globalAndItsOnly_inArray");
  NormaliseYieldToIntegral(H1D_trackEta_cutVariationComp[1]);

  H1D_trackEta_cutVariationComp[2] = (TH1D*)H1D_trackEta_rebinned[iGlobalWithBadTPCCrossedRows]->Clone("H1D_trackEta_interm1_globalAndTpcCrossedRowsRelaxed_inArray");
  NormaliseYieldToIntegral(H1D_trackEta_cutVariationComp[2]);

  H1D_trackEta_cutVariationComp[3] = (TH1D*)H1D_trackEta_rebinned[iUniform]->Clone("H1D_trackEta_interm1_uniform_inArray");
  NormaliseYieldToIntegral(H1D_trackEta_cutVariationComp[3]);







  TString* pdfNameEventNorm = new TString("track_Eta_cutVariationComp_EventNorm");
  TString* pdfNameEventNorm_ratio = new TString("track_Eta_cutVariationComp_EventNorm_ratio");

  TString textContext(contextTrackDatasetComp(""));

  TString cutVariationComp_names[4] = {"global", "global + itsOnly", "globalWithBadTPCCrossedRows", "uniform"};
  Draw_TH1_Histograms_in_one(H1D_trackEta_cutVariationComp, cutVariationComp_names, 4, textContext, pdfNameEventNorm, texEtaX, texTrackEtaYield_EntriesNorm, texCollisionDataInfo, drawnWindowAuto, "");
}

void Draw_Phi_DatasetComparison_trackSelComp() {
  TH2D* H2D_centrality_track[nDatasets];
  TH1D* H1D_trackPhi[nDatasets];
  TH1D* H1D_trackPhi_rebinned[nDatasets];
  TH1D* H1D_trackPhi_rebinned_ratios[nDatasets];

  TH1D* H1D_trackPhi_cutVariationComp[10];

  bool divideSuccess = false;

  int iGlobal = 0;
  int iITSOnly = 2;
  int iGlobalWithBadTPCCrossedRows = 3;
  int iUniform = 1;

  H2D_centrality_track[iGlobal] = (TH2D*)((TH2D*)file_O2Analysis_list[iGlobal]->Get(analysisWorkflow[iGlobal]+"/h2_centrality_track_phi"))->Clone("Draw_Phi_DatasetComparison"+Datasets[iGlobal]);
  H1D_trackPhi[iGlobal] = (TH1D*)H2D_centrality_track[iGlobal]->ProjectionY("trackPhi_"+Datasets[iGlobal], 1, H2D_centrality_track[iGlobal]->GetNbinsX(), "e");
  H1D_trackPhi_rebinned[iGlobal] = (TH1D*)H1D_trackPhi[iGlobal]->Rebin(1.,"trackPhi_rebinned"+Datasets[iGlobal]);
  // NormaliseYieldToNEvents(H1D_trackPhi_rebinned[iGlobal], GetNEventsSelected_TrackEffWorkflow(file_O2Analysis_list[iGlobal]));
  // NormaliseYieldToNEntries(H1D_trackPhi_rebinned[iGlobal]);

  H2D_centrality_track[iITSOnly] = (TH2D*)((TH2D*)file_O2Analysis_list[iITSOnly]->Get(analysisWorkflow[iITSOnly]+"/h2_centrality_track_phi"))->Clone("Draw_Phi_DatasetComparison"+Datasets[iITSOnly]);
  H1D_trackPhi[iITSOnly] = (TH1D*)H2D_centrality_track[iITSOnly]->ProjectionY("trackPhi_"+Datasets[iITSOnly], 1, H2D_centrality_track[iITSOnly]->GetNbinsX(), "e");
  H1D_trackPhi_rebinned[iITSOnly] = (TH1D*)H1D_trackPhi[iITSOnly]->Rebin(1.,"trackPhi_rebinned"+Datasets[iITSOnly]);
  // NormaliseYieldToNEvents(H1D_trackPhi_rebinned[iITSOnly], GetNEventsSelected_TrackEffWorkflow(file_O2Analysis_list[iITSOnly]));
  // NormaliseYieldToNEntries(H1D_trackPhi_rebinned[iITSOnly]);

  H2D_centrality_track[iGlobalWithBadTPCCrossedRows] = (TH2D*)((TH2D*)file_O2Analysis_list[iGlobalWithBadTPCCrossedRows]->Get(analysisWorkflow[iGlobalWithBadTPCCrossedRows]+"/h2_centrality_track_phi"))->Clone("Draw_Phi_DatasetComparison"+Datasets[iGlobalWithBadTPCCrossedRows]);
  H1D_trackPhi[iGlobalWithBadTPCCrossedRows] = (TH1D*)H2D_centrality_track[iGlobalWithBadTPCCrossedRows]->ProjectionY("trackPhi_"+Datasets[iGlobalWithBadTPCCrossedRows], 1, H2D_centrality_track[iGlobalWithBadTPCCrossedRows]->GetNbinsX(), "e");
  H1D_trackPhi_rebinned[iGlobalWithBadTPCCrossedRows] = (TH1D*)H1D_trackPhi[iGlobalWithBadTPCCrossedRows]->Rebin(1.,"trackPhi_rebinned"+Datasets[iGlobalWithBadTPCCrossedRows]);
  // NormaliseYieldToNEvents(H1D_trackPhi_rebinned[iGlobalWithBadTPCCrossedRows], GetNEventsSelected_TrackEffWorkflow(file_O2Analysis_list[iGlobalWithBadTPCCrossedRows]));
  // NormaliseYieldToNEntries(H1D_trackPhi_rebinned[iGlobalWithBadTPCCrossedRows]);
  H1D_trackPhi_rebinned[iGlobalWithBadTPCCrossedRows]->Add(H1D_trackPhi_rebinned[iGlobal], -1.);

  H2D_centrality_track[iUniform] = (TH2D*)((TH2D*)file_O2Analysis_list[iUniform]->Get(analysisWorkflow[iUniform]+"/h2_centrality_track_phi"))->Clone("Draw_Phi_DatasetComparison"+Datasets[iUniform]);
  H1D_trackPhi[iUniform] = (TH1D*)H2D_centrality_track[iUniform]->ProjectionY("trackPhi_"+Datasets[iUniform], 1, H2D_centrality_track[iUniform]->GetNbinsX(), "e");
  H1D_trackPhi_rebinned[iUniform] = (TH1D*)H1D_trackPhi[iUniform]->Rebin(1.,"trackPhi_rebinned"+Datasets[iUniform]);
  // NormaliseYieldToNEvents(H1D_trackPhi_rebinned[iUniform], GetNEventsSelected_TrackEffWorkflow(file_O2Analysis_list[iUniform]));
  // NormaliseYieldToNEntries(H1D_trackPhi_rebinned[iUniform]);

  H1D_trackPhi_cutVariationComp[0] = (TH1D*)H1D_trackPhi_rebinned[iGlobal]->Clone("H1D_trackPhi_global_inArray");
  NormaliseYieldToNEvents(H1D_trackPhi_cutVariationComp[0], GetNEventsSelected_TrackEffWorkflow(file_O2Analysis_list[iGlobal]));
  // NormaliseYieldToIntegral(H1D_trackPhi_cutVariationComp[0]);

  TH1D* H1D_trackPhi_interm1_globalAndItsOnly = (TH1D*)H1D_trackPhi_rebinned[iGlobal]->Clone("H1D_trackPhi_interm1_globalAndItsOnly");
  H1D_trackPhi_interm1_globalAndItsOnly->Add(H1D_trackPhi_rebinned[iITSOnly]);
  H1D_trackPhi_cutVariationComp[1] = (TH1D*)H1D_trackPhi_interm1_globalAndItsOnly->Clone("H1D_trackPhi_interm1_globalAndItsOnly_inArray");
  NormaliseYieldToNEvents(H1D_trackPhi_cutVariationComp[1], GetNEventsSelected_TrackEffWorkflow(file_O2Analysis_list[iGlobal]));
  // NormaliseYieldToIntegral(H1D_trackPhi_cutVariationComp[1]);

  TH1D* H1D_trackPhi_interm1_globalTpcCrossedRowsRelaxed = (TH1D*)H1D_trackPhi_rebinned[iGlobalWithBadTPCCrossedRows]->Clone("H1D_trackPhi_interm1_globalTpcCrossedRowsRelaxed");
  H1D_trackPhi_cutVariationComp[2] = (TH1D*)H1D_trackPhi_interm1_globalTpcCrossedRowsRelaxed->Clone("H1D_trackPhi_interm1_globalTpcCrossedRowsRelaxed_inArray");
  NormaliseYieldToNEvents(H1D_trackPhi_cutVariationComp[2], GetNEventsSelected_TrackEffWorkflow(file_O2Analysis_list[iGlobal]));
  // NormaliseYieldToIntegral(H1D_trackPhi_cutVariationComp[2]);

  H1D_trackPhi_cutVariationComp[3] = (TH1D*)H1D_trackPhi_rebinned[iUniform]->Clone("H1D_trackPhi_uniform_inArray");
  NormaliseYieldToNEvents(H1D_trackPhi_cutVariationComp[3], GetNEventsSelected_TrackEffWorkflow(file_O2Analysis_list[iGlobal]));
  // NormaliseYieldToIntegral(H1D_trackPhi_cutVariationComp[3]);


  // TH1D* H1D_trackPhi_interm1_ItsOnly = (TH1D*)H1D_trackPhi_rebinned[iITSOnly]->Clone("H1D_trackPhi_interm1_ItsOnly");
  // H1D_trackPhi_interm1_globalAndItsOnly->Add(H1D_trackPhi_rebinned[iGlobal]);
  // // H1D_trackPhi_interm1_globalAndItsOnly->Add(H1D_trackPhi_rebinned[iITSOnly]);
  // H1D_trackPhi_cutVariationComp[4] = H1D_trackPhi_interm1_ItsOnly;
  // NormaliseYieldToNEvents(H1D_trackPhi_cutVariationComp[4], GetNEventsSelected_TrackEffWorkflow(file_O2Analysis_list[iGlobal]));

  // cout << "Nevents = " << GetNEventsSelected_TrackEffWorkflow(file_O2Analysis_list[iGlobal]) << " globalTracks" << endl;
  // cout << "Nevents = " << GetNEventsSelected_TrackEffWorkflow(file_O2Analysis_list[iITSOnly]) << " global + itsOnly" << endl;
  // cout << "Nevents = " << GetNEventsSelected_TrackEffWorkflow(file_O2Analysis_list[iGlobalWithBadTPCCrossedRows]) << " global + noCrossedRowsCutTracks" << endl;
  // cout << "Nevents = " << GetNEventsSelected_TrackEffWorkflow(file_O2Analysis_list[iUniform]) << " uniform" << endl;



  // NormaliseYieldToNEntries(H1D_trackEta_rebinned[iDataset]);


  TString* pdfNameEventNorm = new TString("track_Phi_cutVariationComp_EventNorm");
  TString* pdfNameEventNorm_ratio = new TString("track_Phi_cutVariationComp_EventNorm_ratio");

  TString textContext(contextTrackDatasetComp(""));

  TString cutVariationComp_names[4] = {"global", "global + itsOnly", "globalWithBadTPCCrossedRows", "uniform"};
  Draw_TH1_Histograms_in_one(H1D_trackPhi_cutVariationComp, cutVariationComp_names, 4, textContext, pdfNameEventNorm, texPhiX, texTrackPhiYield_EventNorm, texCollisionDataInfo, drawnWindowAuto, "");
}
