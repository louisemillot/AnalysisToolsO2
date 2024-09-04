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
#include "TPolyLine.h"

// #include <RooUnfold.h>
// #include "RooUnfoldResponse.h"
// #include "RooUnfoldBayes.h"
// #include "RooUnfoldBinByBin.h"

//My Libraries
#include "./JetQC_settings.h"
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
// TString contextDataset1D(int iDataset, float* variableRange, const char options[]);
// TString contextDatasetCompAndRadiusAndVarRange(float jetRadius, float* variableRange, const char options[]);
// TString contextPtRange(float* PtRange);
// TString contextEtaRange(float* PtRange);
// TString contextJetRadius(float jetRadius);

//////////// QC plot functions
// Radius comparison
void Draw_Pt_RadiusComparison(int iDataset, float* etaRange);
void Draw_Eta_RadiusComparison(int iDataset, float* PtRange);
void Draw_Phi_RadiusComparison(int iDataset, float* PtRange);
void Draw_jetNTracks_RadiusComparison_withPtRange(int iDataset, float* PtRange);
void Draw_JetArea_vs_JetPt_RadiusComparison(int iDataset, float* PtRange);
void Draw_JetNTracks_vs_JetPt_RadiusComparison(int iDataset, float* PtRange);
void Draw_Pt_ratio_etaNeg_etaPos_RadiusComparison(int iDataset, float* etaRange);
void Draw_JetPhi_vs_JetEta_RadiusComparison(int iDataset);
void Draw_JetEta_vs_JetPt_RadiusComparison(int iDataset, float* PtRange);
void Draw_JetPhi_vs_JetPt_RadiusComparison(int iDataset, float* PtRange);
void Draw_JetTRDratio_vs_JetEta(int iDataset);
void Draw_JetTRDcount_vs_JetEta(int iDataset);
void Draw_Pt_ratio_etaNeg_etaPos_TRDonly_vs_noTRD(int iDataset, float* etaRange);

void Draw_Pt_RadiusComparison_mcp(int iDataset, float* etaRange);


// Dataset comparison
void Draw_Pt_DatasetComparison(float jetRadius, float* etaRange, const char options[]);
void Draw_Eta_DatasetComparison(float jetRadius, float* PtRange, const char options[]);
void Draw_Phi_DatasetComparison(float jetRadius, float* PtRange, const char options[]);
void Draw_Pt_ratio_etaNeg_etaPos_DatasetComparison(float jetRadius, float* etaRange);
void Draw_Area_PtIntegrated_DatasetComparison(float jetRadius, float* PtRange);
void Draw_Rho_vs_Centrality_DatasetComp();
void Draw_Rho_vs_SelectedMultiplicity_DatasetComp();
void Draw_Rho_vs_LeadingJetPt_DatasetComp();
void Draw_BkgFluctuations_vs_Centrality_DatasetComp(std::array<std::array<float, 2>, 2> drawnWindow);
void Draw_SelectedMultiplicity_vs_Centrality_DatasetComp();
void Draw_Rho_vs_SelectedMultiplicity_DatasetCompRatio();
void Draw_Rho_vs_SelectedMultiplicity_DatasetComp_withCutDemarcation();
void Draw_Ncoll_vs_centrality(const char options[]);


// Rebin comparison
void Draw_Area_PtIntegrated_BinningComparison(int iDataset, float jetRadius, float* PtRange);

// Centrality comparison
void Draw_Pt_CentralityComparison(float jetRadius, int iDataset);
void Draw_BkgFluctuations_CentralityProjection(int iDataset, std::array<std::array<float, 2>, 2> drawnWindowZoom, const char options[]);
void Draw_Eta_CentralityComparison(float jetRadius, int iDataset);
void Draw_Phi_CentralityComparison(float jetRadius, int iDataset);
void Draw_BkgFluctuations_withFit_CentralityProjection(int iDataset, std::array<std::array<float, 2>, 2> drawnWindowZoom);
void Draw_Rho_CentralityProjection(int iDataset, const char options[]);
void Draw_Rho_CentralityProjection_DatasetComp(float* centRange, const char options[]);

void Draw_RhoMean_asFunctionOf_Centrality(int iDataset, const char options[]);
// NTracks comp
void Draw_Rho_withFit_NTracksProjection(int iDataset);

// Pt cut comparison
void Draw_Eta_PtCutComparison(float jetRadius, int iDataset, float* PtCuts, int nPtCut, const char options[]);
void Draw_Phi_PtCutComparison(float jetRadius, int iDataset, float* PtCuts, int nPtCut, const char options[]);

// Run 2 vs Run 3 comparison
void Draw_Pt_Run2Run3Comparison_0010Cent_R040(int iDataset, const char options[]);

// MC jet pt resolution
void Draw_jet_resolution_MC_PtRangeComparison(int iDataset, double jetRadius, const char options[]);

void Count_Nevents_perCentClass(int iDataset, const char options[]);



void Draw_Pt_PbPbToPPComparison_HARDCODED(float jetRadius, float* etaRange, const char options[]);

/////////////////////////////////////////////////////
///////////////////// Main Macro ////////////////////
/////////////////////////////////////////////////////

void JetQC() {
  // Load necessary libraries
  LoadLibs();
  // Set the default style
  SetStyle();

  // TString* SaveAs_Title = new TString("");
  TString* texXtitle = new TString("");
  TString* texYtitle = new TString("");
  // TString* Extra = new TString("");

  // int iDataset = 0;

  float etaRangeSym[2] = {-0.5, 0.5};
  float etaRangeNeg[2] = {-0.5, 0};
  float etaRangePos[2] = {0, 0.5};

  float jetRadiusForDataComp = 0.2;
  float jetR02 = 0.2;
  float jetR06 = 0.6;
  

  // const int nPtMinCuts = 7;
  // float jetPtMinCut;
  // float jetPtMinCutArray[nPtMinCuts] = {-999, 0, 5, 10, 15, 20, 40.};
  const int nPtMinCuts = 3;
  float jetPtMinCut;
  float jetPtMinCutArray[nPtMinCuts] = {0, 30, 60};

  for(int iPtMinCut = 0; iPtMinCut < nPtMinCuts; iPtMinCut++){
    jetPtMinCut = jetPtMinCutArray[iPtMinCut];

    float ptRange[2] = {jetPtMinCut, 200};
    float PtRangeZoom0[2] = {jetPtMinCut, 100};
  //   float PtRangeZoom020[2] = {jetPtMinCut, 20};
  //   float PtRangeZoom2030[2] = {20, 30};
  //   float PtRangeZoom3040[2] = {30, 40};
  //   float PtRangeZoom4050[2] = {40, 50};
  //   float PtRangeZoom5060[2] = {50, 60};
  //   float PtRangeZoom8090[2] = {80, 90};

    Draw_Eta_DatasetComparison(jetRadiusForDataComp, ptRange, "normEvents");
    Draw_Phi_DatasetComparison(jetRadiusForDataComp, ptRange, "normEvents");
    Draw_Pt_DatasetComparison(jetRadiusForDataComp, etaRangeSym, "normEvents");

    // Draw_Eta_DatasetComparison(jetR02, ptRange, "normEvents");
    // Draw_Eta_DatasetComparison(jetR06, ptRange, "normEvents");
    // Draw_Area_PtIntegrated_DatasetComparison(jetRadiusForDataComp, ptRange);
    // Draw_Area_PtIntegrated_DatasetComparison(jetR02, ptRange);
    // Draw_Area_PtIntegrated_DatasetComparison(jetR06, ptRange);
 
    for(int iDataset = 0; iDataset < nDatasets; iDataset++){
      // Draw_Eta_RadiusComparison(iDataset, ptRange);
      // Draw_Phi_RadiusComparison(iDataset, ptRange);
      // Draw_jetNTracks_RadiusComparison_withPtRange(iDataset, PtRangeZoom0);
      // Draw_jetNTracks_RadiusComparison_withPtRange(iDataset, PtRangeZoom2030);
      // Draw_jetNTracks_RadiusComparison_withPtRange(iDataset, PtRangeZoom3040);
      // Draw_jetNTracks_RadiusComparison_withPtRange(iDataset, PtRangeZoom4050);
      // Draw_jetNTracks_RadiusComparison_withPtRange(iDataset, PtRangeZoom5060);
      // Draw_jetNTracks_RadiusComparison_withPtRange(iDataset, PtRangeZoom8090);

      // // 2D plots
      // Draw_JetArea_vs_JetPt_RadiusComparison(iDataset, PtRangeZoom0);
      // Draw_JetArea_vs_JetPt_RadiusComparison(iDataset, ptRange);
      // Draw_JetNTracks_vs_JetPt_RadiusComparison(iDataset, PtRangeZoom0);
      // Draw_JetPhi_vs_JetEta_RadiusComparison(iDataset);
      // Draw_JetEta_vs_JetPt_RadiusComparison(iDataset, PtRangeZoom0);
      // Draw_JetPhi_vs_JetPt_RadiusComparison(iDataset, PtRangeZoom0);
    }
  
  //   // Draw_Area_PtIntegrated_BinningComparison(iDataset, 0.6, ptRange);
  }

  // Draw_Pt_PbPbToPPComparison_HARDCODED(jetRadiusForDataComp, etaRangeSym, "normEvents");
  // Draw_Pt_ratio_etaNeg_etaPos_DatasetComparison(jetRadiusForDataComp, etaRangeSym);
  Draw_Ncoll_vs_centrality("");

  for(int iDataset = 0; iDataset < nDatasets; iDataset++){
  //   // Draw_Pt_RadiusComparison(iDataset, etaRangeSym);
  //   // Draw_Pt_RadiusComparison_mcp(iDataset, etaRangeSym);
  //   // Draw_Pt_RadiusComparison(iDataset, etaRangeNeg);
  //   // Draw_Pt_RadiusComparison(iDataset, etaRangePos);
    // Draw_Pt_ratio_etaNeg_etaPos_RadiusComparison(iDataset, etaRangeSym);
    // Draw_Eta_PtCutComparison(jetR02, iDataset, jetPtMinCutArray, nPtMinCuts, "normEntries");
  //   Draw_Eta_PtCutComparison(jetR06, iDataset, jetPtMinCutArray, nPtMinCuts, "normEntries");
    // Draw_Eta_PtCutComparison(jetRadiusForDataComp, iDataset, jetPtMinCutArray, nPtMinCuts, "normEntries");

    // Draw_Phi_PtCutComparison(jetR02, iDataset, jetPtMinCutArray, nPtMinCuts, "normEntries");
  //   Draw_Phi_PtCutComparison(jetR06, iDataset, jetPtMinCutArray, nPtMinCuts, "normEntries");
  //   Draw_Phi_PtCutComparison(jetRadiusForDataComp, iDataset, jetPtMinCutArray, nPtMinCuts, "normEntries");

    // Draw_Pt_CentralityComparison(jetRadiusForDataComp, iDataset);
    // Draw_Eta_CentralityComparison(jetRadiusForDataComp, iDataset);
    // Draw_Phi_CentralityComparison(jetRadiusForDataComp, iDataset);

  //   // Draw_Pt_Run2Run3Comparison_0010Cent_R040(iDataset, "jetCorrected");
  //   // Draw_Pt_Run2Run3Comparison_0010Cent_R040(iDataset, "jetUncorrected");

  // Draw_jet_resolution_MC_PtRangeComparison(iDataset, jetRadiusForDataComp, "");
  // Count_Nevents_perCentClass(iDataset, "");
  }


  // ////// Background //////
  Draw_Rho_vs_Centrality_DatasetComp();

  Draw_Rho_vs_SelectedMultiplicity_DatasetComp();
  // Draw_Rho_vs_SelectedMultiplicity_DatasetCompRatio();

  const std::array<std::array<float, 2>, 2> drawnWindowBkgFluctVsArea = {{{-999, -999}, {-50, 150}}}; // {{xmin, xmax}, {ymin, ymax}}

  // Draw_BkgFluctuations_vs_Centrality_DatasetComp(drawnWindowBkgFluctVsArea);
  // Draw_SelectedMultiplicity_vs_Centrality_DatasetComp();

  std::array<std::array<float, 2>, 2> drawnWindowBkgFluctZoom = {{{-10, 20}, 
                                                                  {0.001, 12}}}; // {{xmin, xmax}, {ymin, ymax}}
  for(int iDataset = 0; iDataset < nDatasets; iDataset++){
    Draw_Rho_CentralityProjection(iDataset, "normEntries");
    Draw_RhoMean_asFunctionOf_Centrality(iDataset,"");
    Draw_BkgFluctuations_withFit_CentralityProjection(iDataset, drawnWindowBkgFluctZoom);
  }
  for(int iCentBin = 0; iCentBin < nCentralityBins; iCentBin++){
    float centRange[2] = {arrayCentralityBinning[iCentBin], arrayCentralityBinning[iCentBin+1]};
      Draw_Rho_CentralityProjection_DatasetComp(centRange, "normEntries");
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




//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////// Context Utilities /////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// TString contextDatasetRadiusCompAndVarRange(int iDataset, float* variableRange, const char options[]){
//   TString texcontextDatasetRadiusCompAndVarRange;
//   if (strstr(options, "pt") != NULL) { //  || strstr(options, "ratio") != NULL not sure why I had this here
//     texcontextDatasetRadiusCompAndVarRange = "#splitline{"+*texDatasetsComparisonCommonDenominator+" "+DatasetsNames[iDataset]+"}{#splitline{2023 QC}{"+contextPtRange(variableRange)+"}}";
//   }
//   if (strstr(options, "eta") != NULL) { //  || strstr(options, "ratio") != NULL not sure why I had this here
//     texcontextDatasetRadiusCompAndVarRange = "#splitline{"+*texDatasetsComparisonCommonDenominator+" "+DatasetsNames[iDataset]+"}{#splitline{2023 QC}{"+contextEtaRange(variableRange)+"}}";
//   }

//   return texcontextDatasetRadiusCompAndVarRange;
// }

// TString contextDatasetCompAndRadiusAndVarRange(float jetRadius, float* variableRange, const char options[]){
//   TString texcontextDatasetCompAndRadiusAndVarRange;
//   if (strstr(options, "pt") != NULL) { //  || strstr(options, "ratio") != NULL not sure why I had this here
//     texcontextDatasetCompAndRadiusAndVarRange = "#splitline{"+*texDatasetsComparisonCommonDenominator+"}{#splitline{"+contextJetRadius(jetRadius)+"}{"+contextPtRange(variableRange)+"}}";
//   }
//   if (strstr(options, "eta") != NULL) { //  || strstr(options, "ratio") != NULL not sure why I had this here
//     texcontextDatasetCompAndRadiusAndVarRange = "#splitline{"+*texDatasetsComparisonCommonDenominator+"}{#splitline{"+contextJetRadius(jetRadius)+"}{"+contextEtaRange(variableRange)+"}}";
//   }

//   return texcontextDatasetCompAndRadiusAndVarRange;
// }

// TString contextDatasetCompAndRadius(float jetRadius, const char options[]){
//   TString texcontextDatasetCompAndRadius;
//   texcontextDatasetCompAndRadius = "#splitline{"+*texDatasetsComparisonCommonDenominator+"}{"+contextJetRadius(jetRadius)+"}";

//   return texcontextDatasetCompAndRadius;
// }

// TString contextDatasetComp(const char options[]){
//   TString texcontextDatasetComp;
//   texcontextDatasetComp = *texDatasetsComparisonCommonDenominator;

//   return texcontextDatasetComp;
// }

// TString contextPtRange(float* PtRange){
//   std::stringstream ss;
//   ss << PtRange[0] << " < #it{p}_{T} < " << PtRange[1];
//   TString textContext((TString)ss.str());
//   // TString texDataset(Form("%.0f", PtRange[0])+" < #it{p}_{T} < "+Form("%.0f", PtRange[1]));
//   return textContext;
// }

// TString contextEtaRange(float* EtaRange){
//   std::stringstream ss;
//   ss << EtaRange[0] << " < #eta < " << EtaRange[1];
//   TString textContext((TString)ss.str());
//   return textContext;
// }

// TString contextJetRadius(float jetRadius){
//   std::stringstream ss;
//   ss << " R = " << jetRadius;
//   TString textContext((TString)ss.str());
//   // TString texDataset(Form("%.0f", PtRange[0])+" < #it{p}_{T} < "+Form("%.0f", PtRange[1]));
//   return textContext;
// }

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////// QC  plot functions /////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void Draw_Pt_RadiusComparison(int iDataset, float* etaRange) {

  TH3D* H3D_jetRjetPtjetEta;
  TH1D* H1D_jetPt[nRadius];
  TH1D* H1D_jetPt_rebinned[nRadius];
  
  H3D_jetRjetPtjetEta = (TH3D*)((TH3D*)file_O2Analysis_list[iDataset]->Get(analysisWorkflow[iDataset]+"/h3_jet_r_jet_pt_jet_eta"+jetFinderQaHistType[iJetFinderQaType]))->Clone("Draw_Pt_RadiusComparison"+Datasets[iDataset]+Form("%.1f", etaRange[0])+"<eta<"+Form("%.1f", etaRange[1]));
  H3D_jetRjetPtjetEta->Sumw2();

  float EtaCutLow = etaRange[0];
  float EtaCutHigh = etaRange[1];
  int ibinEta_low = H3D_jetRjetPtjetEta->GetZaxis()->FindBin(EtaCutLow);
  int ibinEta_high = H3D_jetRjetPtjetEta->GetZaxis()->FindBin(EtaCutHigh);
  if (ibinEta_low == 0) 
    cout << "WARNING: Pt_RadiusComparison is counting the underflow with the chosen etaRange" << endl;
  if (ibinEta_high == H3D_jetRjetPtjetEta->GetZaxis()->GetNbins()+1) 
    cout << "WARNING: Pt_RadiusComparison is counting the overflow with the chosen etaRange" << endl;
  int ibinJetRadius = 0;

  for(int iRadius = 0; iRadius < nRadius; iRadius++){
    ibinJetRadius = H3D_jetRjetPtjetEta->GetXaxis()->FindBin(arrayRadius[iRadius]+GLOBAL_epsilon);

    H1D_jetPt[iRadius] = (TH1D*)H3D_jetRjetPtjetEta->ProjectionY("jetPt_"+RadiusLegend[iRadius]+Form("%.1f", EtaCutLow)+"<eta<"+Form("%.1f", EtaCutHigh), ibinJetRadius, ibinJetRadius, ibinEta_low, ibinEta_high, "e");
    H1D_jetPt_rebinned[iRadius] = (TH1D*)H1D_jetPt[iRadius]->Rebin(1.,"jetPt_rebinned_"+RadiusLegend[iRadius]+Form("%.1f", EtaCutLow)+"<eta<"+Form("%.1f", EtaCutHigh));

    // NormaliseYieldToNEntries(H1D_jetPt_rebinned[iRadius]);
    NormaliseYieldToNEvents(H1D_jetPt_rebinned[iRadius], GetNEventsSelected_JetFramework(file_O2Analysis_list[iDataset]));
  }

  TString* pdfName = new TString("jet_"+jetType[iJetType]+"_"+jetLevel[iJetLevel]+"_"+DatasetsNames[iDataset]+"_Pt_@eta["+Form("%.1f", EtaCutLow)+","+Form("%.1f", EtaCutHigh)+"]"+jetFinderQaHistType[iJetFinderQaType]);

  // TString textContext(contextDatasetRadiusCompAndVarRange(iDataset, etaRange, "eta"));
  TString textContext(contextCustomThreeFields(*texDatasetsComparisonCommonDenominator, Datasets[iDataset], contextEtaRange(etaRange), ""));

  Draw_TH1_Histograms_in_one(H1D_jetPt_rebinned, RadiusLegend, nRadius, textContext, pdfName, texPtX, texJetPtYield_EventNorm, texCollisionDataInfo, drawnWindowAuto, "logy");
}

void Draw_Eta_RadiusComparison(int iDataset, float* PtRange) {

  TH3D* H3D_jetRjetPtjetEta;
  TH1D* H1D_jetEta[nRadius];
  TH1D* H1D_jetEta_rebinned[nRadius];
  
  H3D_jetRjetPtjetEta = (TH3D*)((TH3D*)file_O2Analysis_list[iDataset]->Get(analysisWorkflow[iDataset]+"/h3_jet_r_jet_pt_jet_eta"+jetFinderQaHistType[iJetFinderQaType]))->Clone("Draw_Eta_RadiusComparison"+Datasets[iDataset]+Form("%.1f", PtRange[0])+"<pt<"+Form("%.1f", PtRange[1]));
  H3D_jetRjetPtjetEta->Sumw2();

  float PtCutLow = PtRange[0];
  float PtCutHigh = PtRange[1];
  int ibinPt_low = H3D_jetRjetPtjetEta->GetYaxis()->FindBin(PtCutLow);
  int ibinPt_high = H3D_jetRjetPtjetEta->GetYaxis()->FindBin(PtCutHigh);
  if (ibinPt_low == 0) 
    cout << "WARNING: Eta_RadiusComparison is counting the underflow with the chosen PtRange" << endl;
  if (ibinPt_high == H3D_jetRjetPtjetEta->GetYaxis()->GetNbins()+1) 
    cout << "WARNING: Eta_RadiusComparison is counting the overflow with the chosen PtRange" << endl;
  int ibinJetRadius = 0;

  for(int iRadius = 0; iRadius < nRadius; iRadius++){
    ibinJetRadius = H3D_jetRjetPtjetEta->GetXaxis()->FindBin(arrayRadius[iRadius]+GLOBAL_epsilon);
    H1D_jetEta[iRadius] = (TH1D*)H3D_jetRjetPtjetEta->ProjectionZ("jetEta_"+RadiusLegend[iRadius]+Form("%.1f", PtCutLow)+"<pt<"+Form("%.1f", PtCutHigh), ibinJetRadius, ibinJetRadius, ibinPt_low, ibinPt_high, "e");
    H1D_jetEta_rebinned[iRadius] = (TH1D*)H1D_jetEta[iRadius]->Rebin(1.,"jetEta_rebinned_"+RadiusLegend[iRadius]+Form("%.1f", PtCutLow)+"<pt<"+Form("%.1f", PtCutHigh));

    // NormaliseYieldToNEntries(H1D_jetEta_rebinned[iRadius]);
    NormaliseYieldToNEvents(H1D_jetEta_rebinned[iRadius], GetNEventsSelected_JetFramework(file_O2Analysis_list[iDataset]));
  }

  TString* pdfName = new TString("jet_"+jetType[iJetType]+"_"+jetLevel[iJetLevel]+"_"+DatasetsNames[iDataset]+"_Eta_@pt["+Form("%03.0f", PtCutLow)+","+Form("%03.0f", PtCutHigh)+"]"+jetFinderQaHistType[iJetFinderQaType]);

  // TString textContext(contextDatasetRadiusCompAndVarRange(iDataset, PtRange, "pt"));
  TString textContext(contextCustomThreeFields(*texDatasetsComparisonCommonDenominator, Datasets[iDataset], contextPtRange(PtRange), ""));

  const std::array<std::array<float, 2>, 2> drawnWindow = {{{-1, 1}, {-999, -999}}}; // {{xmin, xmax}, {ymin, ymax}}

  Draw_TH1_Histograms_in_one(H1D_jetEta_rebinned, RadiusLegend, nRadius, textContext, pdfName, texEtaX, texJetEtaYield_EventNorm, texCollisionDataInfo, drawnWindow, "");
}

void Draw_Phi_RadiusComparison(int iDataset, float* PtRange) {

  TH3D* H3D_jetRjetPtjetPhi;
  TH1D* H1D_jetPhi[nRadius];
  TH1D* H1D_jetPhi_rebinned[nRadius];
  
  H3D_jetRjetPtjetPhi = (TH3D*)((TH3D*)file_O2Analysis_list[iDataset]->Get(analysisWorkflow[iDataset]+"/h3_jet_r_jet_pt_jet_phi"+jetFinderQaHistType[iJetFinderQaType]))->Clone("Draw_Phi_RadiusComparison"+Datasets[iDataset]+Form("%.1f", PtRange[0])+"<pt<"+Form("%.1f", PtRange[1]));
  H3D_jetRjetPtjetPhi->Sumw2();

  float PtCutLow = PtRange[0];
  float PtCutHigh = PtRange[1];
  int ibinPt_low = H3D_jetRjetPtjetPhi->GetYaxis()->FindBin(PtCutLow);
  int ibinPt_high = H3D_jetRjetPtjetPhi->GetYaxis()->FindBin(PtCutHigh);
  if (ibinPt_low == 0) 
    cout << "WARNING: Phi_RadiusComparison is counting the underflow with the chosen PtRange" << endl;
  if (ibinPt_high == H3D_jetRjetPtjetPhi->GetYaxis()->GetNbins()) 
    cout << "WARNING: Phi_RadiusComparison is counting the overflow with the chosen PtRange" << endl;
  int ibinJetRadius = 0;

  for(int iRadius = 0; iRadius < nRadius; iRadius++){
    ibinJetRadius = H3D_jetRjetPtjetPhi->GetXaxis()->FindBin(arrayRadius[iRadius]+GLOBAL_epsilon);
    H1D_jetPhi[iRadius] = (TH1D*)H3D_jetRjetPtjetPhi->ProjectionZ("jetPhi_"+RadiusLegend[iRadius]+Form("%.1f", PtCutLow)+"<pt<"+Form("%.1f", PtCutHigh), ibinJetRadius, ibinJetRadius, ibinPt_low, ibinPt_high, "e");
    H1D_jetPhi_rebinned[iRadius] = (TH1D*)H1D_jetPhi[iRadius]->Rebin(1.,"jetPhi_rebinned_"+RadiusLegend[iRadius]+Form("%.1f", PtCutLow)+"<pt<"+Form("%.1f", PtCutHigh));

    // NormaliseYieldToNEntries(H1D_jetPhi_rebinned[iRadius]);
    NormaliseYieldToNEvents(H1D_jetPhi_rebinned[iRadius], GetNEventsSelected_JetFramework(file_O2Analysis_list[iDataset]));
  }
 
  TString* pdfName = new TString("jet_"+jetType[iJetType]+"_"+jetLevel[iJetLevel]+"_"+DatasetsNames[iDataset]+"_Phi_@pt["+Form("%03.0f", PtCutLow)+","+Form("%03.0f", PtCutHigh)+"]"+jetFinderQaHistType[iJetFinderQaType]);

  // TString textContext(contextDatasetRadiusCompAndVarRange(iDataset, PtRange, "pt"));
  TString textContext(contextCustomThreeFields(*texDatasetsComparisonCommonDenominator, Datasets[iDataset], contextPtRange(PtRange), ""));

  Draw_TH1_Histograms_in_one(H1D_jetPhi_rebinned, RadiusLegend, nRadius, textContext, pdfName, texPhiX, texJetPhiYield_EventNorm, texCollisionDataInfo, drawnWindowAuto, "histWithLine");
}

void Draw_jetNTracks_RadiusComparison_withPtRange(int iDataset, float* PtRange) {

  TH3D* H3D_jetRjetPtjetNTracks;
  TH1D* H1D_jetNTracks[nRadius];
  TH1D* H1D_jetNTracks_rebinned[nRadius];
  
  H3D_jetRjetPtjetNTracks = (TH3D*)((TH3D*)file_O2Analysis_list[iDataset]->Get(analysisWorkflow[iDataset]+"/h3_jet_r_jet_pt_jet_ntracks"+jetFinderQaHistType[iJetFinderQaType]))->Clone("Draw_jetNTracks_RadiusComparison_withPtRange"+Datasets[iDataset]+Form("%.1f", PtRange[0])+"<pt<"+Form("%.1f", PtRange[1]));
  H3D_jetRjetPtjetNTracks->Sumw2();

  float PtCutLow = PtRange[0];
  float PtCutHigh = PtRange[1];
  int ibinPt_low = H3D_jetRjetPtjetNTracks->GetYaxis()->FindBin(PtCutLow);
  int ibinPt_high = H3D_jetRjetPtjetNTracks->GetYaxis()->FindBin(PtCutHigh);
  if (ibinPt_low == 0) 
    cout << "WARNING: NTracks_RadiusComparison is counting the underflow with the chosen PtRange" << endl;
  if (ibinPt_high == H3D_jetRjetPtjetNTracks->GetYaxis()->GetNbins()+1) 
    cout << "WARNING: NTracks_RadiusComparison is counting the overflow with the chosen PtRange" << endl;
  int ibinJetRadius = 0;

  for(int iRadius = 0; iRadius < nRadius; iRadius++){
    ibinJetRadius = H3D_jetRjetPtjetNTracks->GetXaxis()->FindBin(arrayRadius[iRadius]+GLOBAL_epsilon);
    H1D_jetNTracks[iRadius] = (TH1D*)H3D_jetRjetPtjetNTracks->ProjectionZ("jetNTracks_"+RadiusLegend[iRadius]+Form("%.1f", PtCutLow)+"<pt<"+Form("%.1f", PtCutHigh), ibinJetRadius, ibinJetRadius, ibinPt_low, ibinPt_high, "e");
    H1D_jetNTracks_rebinned[iRadius] = (TH1D*)H1D_jetNTracks[iRadius]->Rebin(1.,"jetNTracks_rebinned_"+RadiusLegend[iRadius]+Form("%.1f", PtCutLow)+"<pt<"+Form("%.1f", PtCutHigh));

    // NormaliseYieldToNEntries(H1D_jetNTracks_rebinned[iRadius]);
    NormaliseYieldToNEvents(H1D_jetNTracks_rebinned[iRadius], GetNEventsSelected_JetFramework(file_O2Analysis_list[iDataset]));
  }

  TString* pdfName = new TString("jet_"+jetType[iJetType]+"_"+jetLevel[iJetLevel]+"_"+DatasetsNames[iDataset]+"_NTracks_@pT["+Form("%03.0f", PtCutLow)+","+Form("%03.0f", PtCutHigh)+"]"+jetFinderQaHistType[iJetFinderQaType]);

  // TString textContext(contextDatasetRadiusCompAndVarRange(iDataset, PtRange, "pt"));
  TString textContext(contextCustomThreeFields(*texDatasetsComparisonCommonDenominator, Datasets[iDataset], contextPtRange(PtRange), ""));

  Draw_TH1_Histograms_in_one(H1D_jetNTracks_rebinned, RadiusLegend, nRadius, textContext, pdfName, texNTracksX, texJetNormNTracksYield, texCollisionDataInfo, drawnWindowAuto, "");
}



void Draw_JetArea_vs_JetPt_RadiusComparison(int iDataset, float* PtRange) {

  TH3D* H3D_jetRjetPtjetArea;
  TH2D* H2D_jetArea[nRadius];
  TH2D* H2D_jetArea_rebinned[nRadius];
  
  H3D_jetRjetPtjetArea = (TH3D*)((TH3D*)file_O2Analysis_list[iDataset]->Get(analysisWorkflow[iDataset]+"/h3_jet_r_jet_pt_jet_area"+jetFinderQaHistType[iJetFinderQaType]))->Clone("Draw_JetArea_vs_JetPt_RadiusComparison"+Datasets[iDataset]+Form("%.1f", PtRange[0])+"<pt<"+Form("%.1f", PtRange[1]));
  H3D_jetRjetPtjetArea->Sumw2();

  float PtCutLow = PtRange[0];
  float PtCutHigh = PtRange[1];
  int ibinPt_low = H3D_jetRjetPtjetArea->GetYaxis()->FindBin(PtCutLow);
  int ibinPt_high = H3D_jetRjetPtjetArea->GetYaxis()->FindBin(PtCutHigh);

  int ibinJetRadius;

  for(int iRadius = 0; iRadius < nRadius; iRadius++){
    ibinJetRadius = H3D_jetRjetPtjetArea->GetXaxis()->FindBin(arrayRadius[iRadius]+GLOBAL_epsilon);
    H3D_jetRjetPtjetArea->GetXaxis()->SetRange(ibinJetRadius, ibinJetRadius);
    H3D_jetRjetPtjetArea->GetYaxis()->SetRange(ibinPt_low,ibinPt_high);
    H3D_jetRjetPtjetArea->GetZaxis()->SetRange(1,H3D_jetRjetPtjetArea->GetZaxis()->FindBin(areaDisplayMax[iRadius]-GLOBAL_epsilon));

    H2D_jetArea[iRadius] = (TH2D*)H3D_jetRjetPtjetArea->Project3D(RadiusLegend[iRadius]+Form("%d", iDataset)+"_jetArea_e_zy"); //can't use letter D in this or it seems to replace the histogram in current pad (see documentation of ProjectionX function. Isn't mentioned in project3D sadly)
    // H2D_jetArea_rebinned[iRadius] = (TH2D*)H2D_jetArea[iRadius]->RebinY(4.,"H1D_jetArea_rebinned_"+RadiusLegend[iRadius]);
    H2D_jetArea[iRadius]->Scale(1./H2D_jetArea[iRadius]->Integral(1, H2D_jetArea[iRadius]->GetNbinsX(), 1, H2D_jetArea[iRadius]->GetNbinsY(), "width"));
  }

  TString* pdfName = new TString("jet_"+jetType[iJetType]+"_"+jetLevel[iJetLevel]+"_"+DatasetsNames[iDataset]+"_JetArea-vs-Pt_@pT["+Form("%03.0f", PtCutLow)+","+Form("%03.0f", PtCutHigh)+"]"+jetFinderQaHistType[iJetFinderQaType]);

  // TString textContext(contextDatasetRadiusCompAndVarRange(iDataset, PtRange, "pt"));
  // TString textContext(contextCustomThreeFields(*texDatasetsComparisonCommonDenominator, Datasets[iDataset], contextPtRange(PtRange), ""));
  TString textContext(contextCustomOneField(*texDatasetsComparisonCommonDenominator, ""));

  // Draw_TH2_Histograms(H2D_jetArea, RadiusLegend, nRadius, textContext, pdfName, texPtX, texJetArea, texCollisionDataInfo, drawnWindowAuto, "");
  TString* pdfNamelogz = new TString(*pdfName + "_logz");
  Draw_TH2_Histograms(H2D_jetArea, RadiusLegend, nRadius, textContext, pdfNamelogz, texPtX, texJetArea, texCollisionDataInfo, drawnWindowAuto, "logz");
  // TString* pdfNamelogy = new TString(*pdfName + "_logx");
  // Draw_TH2_Histograms(H2D_jetArea, RadiusLegend, nRadius, textContext, pdfNamelogy, texPtX, texJetArea, texCollisionDataInfo, drawnWindowAuto, "logx");
  // TString* pdfNamelogyz = new TString(*pdfName + "_logxz");
  // Draw_TH2_Histograms(H2D_jetArea, RadiusLegend, nRadius, textContext, pdfNamelogyz, texPtX, texJetArea, texCollisionDataInfo, drawnWindowAuto, "logxlogz");
}

void Draw_JetNTracks_vs_JetPt_RadiusComparison(int iDataset, float* PtRange) {

  TH3D* H3D_jetRjetPtjetNTracks;
  TH2D* H2D_jetNTracks[nRadius];
  TH2D* H2D_jetNTracks_rebinned[nRadius];
  
  H3D_jetRjetPtjetNTracks = (TH3D*)((TH3D*)file_O2Analysis_list[iDataset]->Get(analysisWorkflow[iDataset]+"/h3_jet_r_jet_pt_jet_ntracks"+jetFinderQaHistType[iJetFinderQaType]))->Clone("Draw_JetNTracks_vs_JetPt_RadiusComparison"+Datasets[iDataset]+Form("%.1f", PtRange[0])+"<pt<"+Form("%.1f", PtRange[1]));
  H3D_jetRjetPtjetNTracks->Sumw2();

  float PtCutLow = PtRange[0];
  float PtCutHigh = PtRange[1];
  int ibinPt_low = H3D_jetRjetPtjetNTracks->GetYaxis()->FindBin(PtCutLow);
  int ibinPt_high = H3D_jetRjetPtjetNTracks->GetYaxis()->FindBin(PtCutHigh);
  int ibinJetRadius;

  for(int iRadius = 0; iRadius < nRadius; iRadius++){
    ibinJetRadius = H3D_jetRjetPtjetNTracks->GetXaxis()->FindBin(arrayRadius[iRadius]+GLOBAL_epsilon);
    H3D_jetRjetPtjetNTracks->GetXaxis()->SetRange(ibinJetRadius, ibinJetRadius);
    H3D_jetRjetPtjetNTracks->GetYaxis()->SetRange(ibinPt_low,ibinPt_high);
    H2D_jetNTracks[iRadius] = (TH2D*)H3D_jetRjetPtjetNTracks->Project3D(RadiusLegend[iRadius]+Form("%d", iDataset)+"_jetNTracks_e_zy"); //can't use letter D in this or it seems to replace the histogram in current pad (see documentation of ProjectionX function. Isn't mentioned in project3D sadly)
    // H2D_jetNTracks_rebinned[iRadius] = (TH2D*)H2D_jetNTracks[iRadius]->Rebin(1.,"jetNTracks_rebinned_"+RadiusLegend[iRadius]);
  }

  TString* pdfName = new TString("jet_"+jetType[iJetType]+"_"+jetLevel[iJetLevel]+"_"+DatasetsNames[iDataset]+"_JetNTracks-vs-Pt_@pT["+Form("%03.0f", PtCutLow)+","+Form("%03.0f", PtCutHigh)+"]"+jetFinderQaHistType[iJetFinderQaType]);

  // TString textContext(contextDatasetRadiusCompAndVarRange(iDataset, PtRange, "pt"));
  TString textContext(contextCustomThreeFields(*texDatasetsComparisonCommonDenominator, Datasets[iDataset], contextPtRange(PtRange), ""));

  Draw_TH2_Histograms(H2D_jetNTracks, RadiusLegend, nRadius, textContext, pdfName, texPtX, texJetNTracks, texCollisionDataInfo, drawnWindowAuto, "");
  TString* pdfNamelogz = new TString(*pdfName + "_logz");
  Draw_TH2_Histograms(H2D_jetNTracks, RadiusLegend, nRadius, textContext, pdfNamelogz, texPtX, texJetNTracks, texCollisionDataInfo, drawnWindowAuto, "logz");
}

void Draw_JetEta_vs_JetPt_RadiusComparison(int iDataset, float* PtRange) {

  TH3D* H3D_jetRjetPtjetEta;
  TH2D* H2D_jetEta[nRadius];
  TH2D* H2D_jetEta_rebinned[nRadius];
  
  H3D_jetRjetPtjetEta = (TH3D*)((TH3D*)file_O2Analysis_list[iDataset]->Get(analysisWorkflow[iDataset]+"/h3_jet_r_jet_pt_jet_eta"+jetFinderQaHistType[iJetFinderQaType]))->Clone("Draw_JetEta_vs_JetPt_RadiusComparison"+Datasets[iDataset]+Form("%.1f", PtRange[0])+"<pt<"+Form("%.1f", PtRange[1]));
  H3D_jetRjetPtjetEta->Sumw2();

  float PtCutLow = PtRange[0];
  float PtCutHigh = PtRange[1];
  int ibinPt_low = H3D_jetRjetPtjetEta->GetYaxis()->FindBin(PtCutLow);
  int ibinPt_high = H3D_jetRjetPtjetEta->GetYaxis()->FindBin(PtCutHigh);
  int ibinJetRadius;

  for(int iRadius = 0; iRadius < nRadius; iRadius++){
    ibinJetRadius = H3D_jetRjetPtjetEta->GetXaxis()->FindBin(arrayRadius[iRadius]+GLOBAL_epsilon);
    H3D_jetRjetPtjetEta->GetXaxis()->SetRange(ibinJetRadius, ibinJetRadius);
    H3D_jetRjetPtjetEta->GetYaxis()->SetRange(ibinPt_low,ibinPt_high);
    H2D_jetEta[iRadius] = (TH2D*)H3D_jetRjetPtjetEta->Project3D(RadiusLegend[iRadius]+Form("%d", iDataset)+Form("%.1f", PtRange[0])+"<pt<"+Form("%.1f", PtRange[1])+"_jetEta_jetPt_e_zy"); //can't use letter D in this or it seems to replace the histogram in current pad (see documentation of ProjectionX function. Isn't mentioned in project3D sadly)
    // H2D_jetEta_rebinned[iRadius] = (TH2D*)H2D_jetEta[iRadius]->RebinY(1.,"H1D_jetEta_rebinned_"+RadiusLegend[iRadius]);
  }

  TString* pdfName = new TString("jet_"+jetType[iJetType]+"_"+jetLevel[iJetLevel]+"_"+Form("%d", iDataset)+"_JetEta-vs-Pt_@pT["+Form("%03.0f", PtCutLow)+","+Form("%03.0f", PtCutHigh)+"]"+jetFinderQaHistType[iJetFinderQaType]);

  // TString textContext(contextDatasetRadiusCompAndVarRange(iDataset, PtRange, "pt"));
  TString textContext(contextCustomThreeFields(*texDatasetsComparisonCommonDenominator, Datasets[iDataset], contextPtRange(PtRange), ""));

  // Draw_TH2_Histograms(H2D_jetEta, RadiusLegend, nRadius, textContext, pdfName, texPtX, texEtaX, texCollisionDataInfo, drawnWindowAuto, "");
  TString* pdfNamelogz = new TString(*pdfName + "_logz");
  Draw_TH2_Histograms(H2D_jetEta, RadiusLegend, nRadius, textContext, pdfNamelogz, texPtX, texEtaX, texCollisionDataInfo, drawnWindowAuto, "logz");
  // TString* pdfNamelogy = new TString(*pdfName + "_logx");
  // Draw_TH2_Histograms(H2D_jetEta, RadiusLegend, nRadius, textContext, pdfNamelogy, texPtX, texEtaX, texCollisionDataInfo, drawnWindowAuto, "logx");
  // TString* pdfNamelogyz = new TString(*pdfName + "_logxz");
  // Draw_TH2_Histograms(H2D_jetEta, RadiusLegend, nRadius, textContext, pdfNamelogyz, texPtX, texEtaX, texCollisionDataInfo, drawnWindowAuto, "logxlogz");
}

void Draw_JetPhi_vs_JetPt_RadiusComparison(int iDataset, float* PtRange) {

  TH3D* H3D_jetRjetPtjetPhi;
  TH2D* H2D_jetPhi[nRadius];
  TH2D* H2D_jetPhi_rebinned[nRadius];
  
  H3D_jetRjetPtjetPhi = (TH3D*)((TH3D*)file_O2Analysis_list[iDataset]->Get(analysisWorkflow[iDataset]+"/h3_jet_r_jet_pt_jet_phi"+jetFinderQaHistType[iJetFinderQaType]))->Clone("Draw_JetPhi_vs_JetPt_RadiusComparison"+Datasets[iDataset]+Form("%.1f", PtRange[0])+"<pt<"+Form("%.1f", PtRange[1]));
  H3D_jetRjetPtjetPhi->Sumw2();

  float PtCutLow = PtRange[0];
  float PtCutHigh = PtRange[1];
  int ibinPt_low = H3D_jetRjetPtjetPhi->GetYaxis()->FindBin(PtCutLow);
  int ibinPt_high = H3D_jetRjetPtjetPhi->GetYaxis()->FindBin(PtCutHigh);
  int ibinJetRadius;

  for(int iRadius = 0; iRadius < nRadius; iRadius++){
    ibinJetRadius = H3D_jetRjetPtjetPhi->GetXaxis()->FindBin(arrayRadius[iRadius]+GLOBAL_epsilon);
    H3D_jetRjetPtjetPhi->GetXaxis()->SetRange(ibinJetRadius, ibinJetRadius);
    H3D_jetRjetPtjetPhi->GetYaxis()->SetRange(ibinPt_low,ibinPt_high);
    H2D_jetPhi[iRadius] = (TH2D*)H3D_jetRjetPtjetPhi->Project3D(RadiusLegend[iRadius]+Form("%d", iDataset)+Form("%.1f", PtRange[0])+"<pt<"+Form("%.1f", PtRange[1])+"_jetPhi_jetPt_e_zy"); //can't use letter D in this or it seems to replace the histogram in current pad (see documentation of ProjectionX function. Isn't mentioned in project3D sadly)
    // H2D_jetPhi_rebinned[iRadius] = (TH2D*)H2D_jetPhi[iRadius]->RebinY(1.,"H1D_jetPhi_rebinned_"+RadiusLegend[iRadius]);
  }

  TString* pdfName = new TString("jet_"+jetType[iJetType]+"_"+jetLevel[iJetLevel]+"_"+DatasetsNames[iDataset]+"_JetPhi-vs-Pt_@pT["+Form("%03.0f", PtCutLow)+","+Form("%03.0f", PtCutHigh)+"]"+jetFinderQaHistType[iJetFinderQaType]);

  // TString textContext(contextDatasetRadiusCompAndVarRange(iDataset, PtRange, "pt"));
  TString textContext(contextCustomThreeFields(*texDatasetsComparisonCommonDenominator, Datasets[iDataset], contextPtRange(PtRange), ""));

  // Draw_TH2_Histograms(H2D_jetPhi, RadiusLegend, nRadius, textContext, pdfName, texPtX, texPhiX, texCollisionDataInfo, drawnWindowAuto, "");
  TString* pdfNamelogz = new TString(*pdfName + "_logz");
  Draw_TH2_Histograms(H2D_jetPhi, RadiusLegend, nRadius, textContext, pdfNamelogz, texPtX, texPhiX, texCollisionDataInfo, drawnWindowAuto, "logz");
  // TString* pdfNamelogy = new TString(*pdfName + "_logx");
  // Draw_TH2_Histograms(H2D_jetPhi, RadiusLegend, nRadius, textContext, pdfNamelogy, texPtX, texPhiX, texCollisionDataInfo, drawnWindowAuto, "logx");
  // TString* pdfNamelogyz = new TString(*pdfName + "_logxz");
  // Draw_TH2_Histograms(H2D_jetPhi, RadiusLegend, nRadius, textContext, pdfNamelogyz, texPtX, texPhiX, texCollisionDataInfo, drawnWindowAuto, "logxlogz");
}

void Draw_JetPhi_vs_JetEta_RadiusComparison(int iDataset) {

  TH3D* H3D_jetRjetPhijetEta;
  TH2D* H2D_jetPhijetEta[nRadius];
  TH2D* H2D_jetPhijetEta_rebinned[nRadius];
  
  H3D_jetRjetPhijetEta = (TH3D*)((TH3D*)file_O2Analysis_list[iDataset]->Get(analysisWorkflow[iDataset]+"/h3_jet_r_jet_eta_jet_phi"+jetFinderQaHistType[iJetFinderQaType]))->Clone("Draw_JetPhi_vs_JetEta_RadiusComparison"+Datasets[iDataset]);
  H3D_jetRjetPhijetEta->Sumw2();

  int ibinJetRadius;

  for(int iRadius = 0; iRadius < nRadius; iRadius++){
    ibinJetRadius = H3D_jetRjetPhijetEta->GetXaxis()->FindBin(arrayRadius[iRadius]+GLOBAL_epsilon);
    H3D_jetRjetPhijetEta->GetXaxis()->SetRange(ibinJetRadius, ibinJetRadius);
    H2D_jetPhijetEta[iRadius] = (TH2D*)H3D_jetRjetPhijetEta->Project3D(RadiusLegend[iRadius]+Form("%d", iDataset)+"_jetPhi_jetEta_e_zy"); //can't use letter D in this or it seems to replace the histogram in current pad (see documentation of ProjectionX function. Isn't mentioned in project3D sadly)
    // H2D_jetPhijetEta_rebinned[iRadius] = (TH2D*)H2D_jetArea[iRadius]->RebinY(2.,"H1D_jetPhi_jetEta_rebinned_"+RadiusLegend[iRadius]);
  }

  TString* pdfName = new TString("jet_"+jetType[iJetType]+"_"+jetLevel[iJetLevel]+"_"+DatasetsNames[iDataset]+"_JetPhi-vs-Eta"+jetFinderQaHistType[iJetFinderQaType]);

  TString textContext("#splitline{"+*texDatasetsComparisonCommonDenominator+" "+DatasetsNames[iDataset]+"}{2023 QC}");

  // Draw_TH2_Histograms(H2D_jetPhijetEta, RadiusLegend, nRadius, textContext, pdfName, texEtaX, texPhiX, texCollisionDataInfo, drawnWindowAuto, "");
  TString* pdfNamelogz = new TString(*pdfName + "_logz");
  Draw_TH2_Histograms(H2D_jetPhijetEta, RadiusLegend, nRadius, textContext, pdfNamelogz, texPtX, texJetArea, texCollisionDataInfo, drawnWindowAuto, "logz");
  // TString* pdfNamelogy = new TString(*pdfName + "_logx");
  // Draw_TH2_Histograms(H2D_jetPhijetEta, RadiusLegend, nRadius, textContext, pdfNamelogy, texPtX, texJetArea, texCollisionDataInfo, drawnWindowAuto, "logx");
  // TString* pdfNamelogyz = new TString(*pdfName + "_logxz");
  // Draw_TH2_Histograms(H2D_jetPhijetEta, RadiusLegend, nRadius, textContext, pdfNamelogyz, texPtX, texJetArea, texCollisionDataInfo, drawnWindowAuto, "logxlogz");
  
}

void Draw_Pt_DatasetComparison(float jetRadius, float* etaRange, const char options[]) {

  TH3D* H3D_jetRjetPtjetEta[nDatasets];
  TH1D* H1D_jetPt[nDatasets];
  TH1D* H1D_jetPt_rebinned[nDatasets];
  
  TH1D* H1D_jetPt_rebinned_ratios[nDatasets];

  float EtaCutLow = etaRange[0];
  float EtaCutHigh = etaRange[1];
  int ibinJetRadius = 0;

  bool divideSuccess = false;

  TString* yAxisLabel;
  for(int iDataset = 0; iDataset < nDatasets; iDataset++){

    H3D_jetRjetPtjetEta[iDataset] = (TH3D*)((TH3D*)file_O2Analysis_list[iDataset]->Get(analysisWorkflow[iDataset]+"/h3_jet_r_jet_pt_jet_eta"+jetFinderQaHistType[iJetFinderQaType]))->Clone("Draw_Pt_DatasetComparison"+Datasets[iDataset]+"Radius"+Form("%.1f",jetRadius)+jetRadius+Form("%.1f", etaRange[0])+"<eta<"+Form("%.1f", etaRange[1]));

    int ibinEta_low = H3D_jetRjetPtjetEta[iDataset]->GetZaxis()->FindBin(EtaCutLow);
    int ibinEta_high = H3D_jetRjetPtjetEta[iDataset]->GetZaxis()->FindBin(EtaCutHigh);
    if (ibinEta_low == 0) 
      cout << "WARNING: Pt_DatasetComparison is counting the underflow with the chosen etaRange" << endl;
    if (ibinEta_high == H3D_jetRjetPtjetEta[iDataset]->GetZaxis()->GetNbins()+1) 
      cout << "WARNING: Pt_DatasetComparison is counting the overflow with the chosen etaRange" << endl;
    ibinJetRadius = H3D_jetRjetPtjetEta[iDataset]->GetXaxis()->FindBin(jetRadius+GLOBAL_epsilon);
 
    H1D_jetPt[iDataset] = (TH1D*)H3D_jetRjetPtjetEta[iDataset]->ProjectionY("jetPt_"+Datasets[iDataset]+"Radius"+Form("%.1f",jetRadius)+Form("%.1f", EtaCutLow)+"<eta<"+Form("%.1f", EtaCutHigh), ibinJetRadius, ibinJetRadius, ibinEta_low, ibinEta_high, "e");
    H1D_jetPt_rebinned[iDataset] = (TH1D*)H1D_jetPt[iDataset]->Rebin(5.,"jetPt_rebinned_"+Datasets[iDataset]+"Radius"+Form("%.1f",jetRadius)+Form("%.1f", EtaCutLow)+"<eta<"+Form("%.1f", EtaCutHigh));

    if (strstr(options, "normEntries") != NULL) {
      NormaliseYieldToNEntries(H1D_jetPt_rebinned[iDataset]);
      yAxisLabel = texJetPtYield_EntriesNorm;
    }
    int Nevents;
    if (strstr(options, "normEvents") != NULL) {
      if (isDatasetWeighted[iDataset]) {
        Nevents = GetNEventsSelected_JetFramework_weighted(file_O2Analysis_list[iDataset]);
      } else {
        Nevents = GetNEventsSelected_JetFramework(file_O2Analysis_list[iDataset]);
      }
      NormaliseYieldToNEvents(H1D_jetPt_rebinned[iDataset], Nevents);
      yAxisLabel = texJetPtYield_EventNorm;
    }

    H1D_jetPt_rebinned_ratios[iDataset] = (TH1D*)H1D_jetPt_rebinned[iDataset]->Clone("jetPt_rebinned_ratios"+Datasets[iDataset]+"Radius"+Form("%.1f",jetRadius)+Form("%.1f", EtaCutLow)+"<eta<"+Form("%.1f", EtaCutHigh));
    H1D_jetPt_rebinned_ratios[iDataset]->Reset("M");
    divideSuccess = H1D_jetPt_rebinned_ratios[iDataset]->Divide(H1D_jetPt_rebinned[iDataset], H1D_jetPt_rebinned[0]);
  }

  TString* pdfName = new TString("jet_"+jetType[iJetType]+"_"+jetLevel[iJetLevel]+"_DataComp_R="+Form("%.1f", jetRadius)+"_Pt_@eta["+Form("%.1f", EtaCutLow)+","+Form("%.1f", EtaCutHigh)+"]"+jetFinderQaHistType[iJetFinderQaType]);
  TString* pdfName_ratio = new TString("jet_"+jetType[iJetType]+"_"+jetLevel[iJetLevel]+"_DataComp_R="+Form("%.1f", jetRadius)+"_Pt_@eta["+Form("%.1f", EtaCutLow)+","+Form("%.1f", EtaCutHigh)+"]"+jetFinderQaHistType[iJetFinderQaType]+"_ratio");

  // TString textContext(contextDatasetCompAndRadiusAndVarRange(jetRadius, etaRange, "eta"));
  TString textContext(contextCustomThreeFields(*texDatasetsComparisonCommonDenominator, "", "#splitline{"+contextJetRadius(jetRadius)+"}{"+contextEtaRange(etaRange)+"}", ""));

  Draw_TH1_Histograms_in_one(H1D_jetPt_rebinned, DatasetsNames, nDatasets, textContext, pdfName, texPtX, texJetPtYield_EventNorm, texCollisionDataInfo, drawnWindowAuto, "logy");
  if (divideSuccess == true) {
    Draw_TH1_Histograms_in_one(H1D_jetPt_rebinned_ratios, DatasetsNames, nDatasets, textContext, pdfName_ratio, texPtX, texRatioDatasets, texCollisionDataInfo, drawnWindowAuto, "standardratio,avoidFirst");
  }
  else {
    cout << "Divide failed in Draw_Pt_DatasetComparison" << endl;
  }
}

void Draw_Eta_DatasetComparison(float jetRadius, float* PtRange, const char options[]) {

  TH3D* H3D_jetRjetPtjetEta[nDatasets];
  TH1D* H1D_jetEta[nDatasets];
  TH1D* H1D_jetEta_rebinned[nDatasets];
  
  TH1D* H1D_jetEta_rebinned_ratios[nDatasets];

  float PtCutLow = PtRange[0];
  float PtCutHigh = PtRange[1];
  int ibinJetRadius = 0;

  bool divideSuccess = false;

  TString* yAxisLabel;

  for(int iDataset = 0; iDataset < nDatasets; iDataset++){

    H3D_jetRjetPtjetEta[iDataset] = (TH3D*)((TH3D*)file_O2Analysis_list[iDataset]->Get(analysisWorkflow[iDataset]+"/h3_jet_r_jet_pt_jet_eta"+jetFinderQaHistType[iJetFinderQaType]))->Clone("Draw_Pt_DatasetComparison"+Datasets[iDataset]+"Radius"+Form("%.1f",jetRadius)+Form("%.1f", PtRange[0])+"<pt<"+Form("%.1f", PtRange[1]));

    int ibinPt_low = H3D_jetRjetPtjetEta[iDataset]->GetYaxis()->FindBin(PtCutLow);
    int ibinPt_high = H3D_jetRjetPtjetEta[iDataset]->GetYaxis()->FindBin(PtCutHigh);
    if (ibinPt_low == 0) 
      cout << "WARNING: Eta_DatasetComparison is counting the underflow with the chosen PtRange" << endl;
    if (ibinPt_high == H3D_jetRjetPtjetEta[iDataset]->GetYaxis()->GetNbins()+1) 
      cout << "WARNING: Eta_DatasetComparison is counting the overflow with the chosen PtRange" << endl;
    ibinJetRadius = H3D_jetRjetPtjetEta[iDataset]->GetXaxis()->FindBin(jetRadius+GLOBAL_epsilon);

    H1D_jetEta[iDataset] = (TH1D*)H3D_jetRjetPtjetEta[iDataset]->ProjectionZ("jetEta_"+Datasets[iDataset]+"Radius"+Form("%.1f",jetRadius)+Form("%.1f", PtCutLow)+"<pt<"+Form("%.1f", PtCutHigh), ibinJetRadius, ibinJetRadius, ibinPt_low, ibinPt_high, "e");
    H1D_jetEta_rebinned[iDataset] = (TH1D*)H1D_jetEta[iDataset]->Rebin(1.,"jetEta_rebinned_"+Datasets[iDataset]+"Radius"+Form("%.1f",jetRadius)+Form("%.1f", PtCutLow)+"<pt<"+Form("%.1f", PtCutHigh));

    if (strstr(options, "normEntries") != NULL) {
      NormaliseYieldToNEntries(H1D_jetEta_rebinned[iDataset]);
      yAxisLabel = texJetEtaYield_EntriesNorm;
    }
    int Nevents;
    if (strstr(options, "normEvents") != NULL) {
      if (isDatasetWeighted[iDataset]) {
        Nevents = GetNEventsSelected_JetFramework_weighted(file_O2Analysis_list[iDataset]);
      } else {
        Nevents = GetNEventsSelected_JetFramework(file_O2Analysis_list[iDataset]);
      }
      NormaliseYieldToNEvents(H1D_jetEta_rebinned[iDataset], Nevents);
      yAxisLabel = texJetEtaYield_EventNorm;
    }
  
    H1D_jetEta_rebinned_ratios[iDataset] = (TH1D*)H1D_jetEta_rebinned[iDataset]->Clone("jetEta_rebinned_ratios"+Datasets[iDataset]+"Radius"+Form("%.1f",jetRadius)+Form("%.1f", PtCutLow)+"<pt<"+Form("%.1f", PtCutHigh));
    H1D_jetEta_rebinned_ratios[iDataset]->Reset("M");
    divideSuccess = H1D_jetEta_rebinned_ratios[iDataset]->Divide(H1D_jetEta_rebinned[iDataset], H1D_jetEta_rebinned[0]);
  }

  TString* pdfName = new TString("jet_"+jetType[iJetType]+"_"+jetLevel[iJetLevel]+"_DataComp_R="+Form("%.1f", jetRadius)+"_Eta_@pT["+Form("%03.0f", PtCutLow)+","+Form("%03.0f", PtCutHigh)+"]"+jetFinderQaHistType[iJetFinderQaType]);
  TString* pdfName_ratio = new TString("jet_"+jetType[iJetType]+"_"+jetLevel[iJetLevel]+"_DataComp_R="+Form("%.1f", jetRadius)+"_Eta_@pT["+Form("%03.0f", PtCutLow)+","+Form("%03.0f", PtCutHigh)+"]"+jetFinderQaHistType[iJetFinderQaType]+"_ratio");

  // TString textContext(contextDatasetCompAndRadiusAndVarRange(jetRadius, PtRange, "pt"));
  TString textContext(contextCustomThreeFields(*texDatasetsComparisonCommonDenominator, "", "#splitline{"+contextJetRadius(jetRadius)+"}{"+contextPtRange(PtRange)+"}", ""));

  const std::array<std::array<float, 2>, 2> drawnWindow = {{{-1, 1}, {-999, -999}}}; // {{xmin, xmax}, {ymin, ymax}}

  Draw_TH1_Histograms_in_one(H1D_jetEta_rebinned, DatasetsNames, nDatasets, textContext, pdfName, texEtaX, yAxisLabel, texCollisionDataInfo, drawnWindow, "");

  if (divideSuccess == true) {
    Draw_TH1_Histograms_in_one(H1D_jetEta_rebinned_ratios, DatasetsNames, nDatasets, textContext, pdfName_ratio, texEtaX, texRatioDatasets, texCollisionDataInfo, drawnWindow, "standardratio,avoidFirst");
  }
  else {
    cout << "Divide failed in Draw_Eta_DatasetComparison" << endl;
  }
}

void Draw_Phi_DatasetComparison(float jetRadius, float* PtRange, const char options[]) {

  TH3D* H3D_jetRjetPtjetPhi[nDatasets];
  TH1D* H1D_jetPhi[nDatasets];
  TH1D* H1D_jetPhi_rebinned[nDatasets];
  
  TH1D* H1D_jetPhi_rebinned_ratios[nDatasets];

  float PtCutLow = PtRange[0];
  float PtCutHigh = PtRange[1];
  int ibinJetRadius = 0;

  bool divideSuccess = false;

  TString* yAxisLabel;

  for(int iDataset = 0; iDataset < nDatasets; iDataset++){

    H3D_jetRjetPtjetPhi[iDataset] = (TH3D*)((TH3D*)file_O2Analysis_list[iDataset]->Get(analysisWorkflow[iDataset]+"/h3_jet_r_jet_pt_jet_phi"+jetFinderQaHistType[iJetFinderQaType]))->Clone("Draw_Phi_DatasetComparison"+Datasets[iDataset]+"Radius"+Form("%.1f",jetRadius)+Form("%.1f", PtRange[0])+"<pt<"+Form("%.1f", PtRange[1]));

    int ibinPt_low = H3D_jetRjetPtjetPhi[iDataset]->GetYaxis()->FindBin(PtCutLow);
    int ibinPt_high = H3D_jetRjetPtjetPhi[iDataset]->GetYaxis()->FindBin(PtCutHigh);
    if (ibinPt_low == 0) 
      cout << "WARNING: Phi_DatasetComparison is counting the underflow with the chosen PtRange" << endl;
    if (ibinPt_high == H3D_jetRjetPtjetPhi[iDataset]->GetYaxis()->GetNbins()+1) 
      cout << "WARNING: Phi_DatasetComparison is counting the overflow with the chosen PtRange" << endl;
    ibinJetRadius = H3D_jetRjetPtjetPhi[iDataset]->GetXaxis()->FindBin(jetRadius+GLOBAL_epsilon);

    H1D_jetPhi[iDataset] = (TH1D*)H3D_jetRjetPtjetPhi[iDataset]->ProjectionZ("jetPhi_"+Datasets[iDataset]+"Radius"+Form("%.1f",jetRadius)+Form("%.1f", PtCutLow)+"<pt<"+Form("%.1f", PtCutHigh), ibinJetRadius,ibinJetRadius, ibinPt_low, ibinPt_high, "e");
    H1D_jetPhi_rebinned[iDataset] = (TH1D*)H1D_jetPhi[iDataset]->Rebin(5.,"jetPhi_rebinned_"+Datasets[iDataset]+"Radius"+Form("%.1f",jetRadius)+Form("%.1f", PtCutLow)+"<pt<"+Form("%.1f", PtCutHigh));

    if (strstr(options, "normEntries") != NULL) {
      NormaliseYieldToNEntries(H1D_jetPhi_rebinned[iDataset]);
      yAxisLabel = texJetPhiYield_EntriesNorm;
    }
    int Nevents;
    if (strstr(options, "normEvents") != NULL) {
      if (isDatasetWeighted[iDataset]) {
        Nevents = GetNEventsSelected_JetFramework_weighted(file_O2Analysis_list[iDataset]);
      } else {
        Nevents = GetNEventsSelected_JetFramework(file_O2Analysis_list[iDataset]);
      }
      NormaliseYieldToNEvents(H1D_jetPhi_rebinned[iDataset], Nevents);
      yAxisLabel = texJetPhiYield_EventNorm;
    }

    H1D_jetPhi_rebinned_ratios[iDataset] = (TH1D*)H1D_jetPhi_rebinned[iDataset]->Clone("jetPhi_rebinned_ratios"+Datasets[iDataset]+"Radius"+Form("%.1f",jetRadius)+Form("%.1f", PtCutLow)+"<pt<"+Form("%.1f", PtCutHigh));
    H1D_jetPhi_rebinned_ratios[iDataset]->Reset("M");
    divideSuccess = H1D_jetPhi_rebinned_ratios[iDataset]->Divide(H1D_jetPhi_rebinned[iDataset], H1D_jetPhi_rebinned[0]);
  }

  TString* pdfName = new TString("jet_"+jetType[iJetType]+"_"+jetLevel[iJetLevel]+"_DataComp_R="+Form("%.1f", jetRadius)+"_Phi_@pT["+Form("%03.0f", PtCutLow)+","+Form("%03.0f", PtCutHigh)+"]"+jetFinderQaHistType[iJetFinderQaType]);
  TString* pdfName_ratio = new TString("jet_"+jetType[iJetType]+"_"+jetLevel[iJetLevel]+"_DataComp_R="+Form("%.1f", jetRadius)+"_Phi_@pT["+Form("%03.0f", PtCutLow)+","+Form("%03.0f", PtCutHigh)+"]"+jetFinderQaHistType[iJetFinderQaType]+"_ratio");

  // TString textContext(contextDatasetCompAndRadiusAndVarRange(jetRadius, PtRange, "pt"));
  TString textContext(contextCustomThreeFields(*texDatasetsComparisonCommonDenominator, "", "#splitline{"+contextJetRadius(jetRadius)+"}{"+contextPtRange(PtRange)+"}", ""));

  Draw_TH1_Histograms_in_one(H1D_jetPhi_rebinned, DatasetsNames, nDatasets, textContext, pdfName, texPhiX, texJetPhiYield_EventNorm, texCollisionDataInfo, drawnWindowAuto, "histWithLine");

  if (divideSuccess == true) {
    Draw_TH1_Histograms_in_one(H1D_jetPhi_rebinned_ratios, DatasetsNames, nDatasets, textContext, pdfName_ratio, texPhiX, texRatioDatasets, texCollisionDataInfo, drawnWindowAuto, "standardratio,avoidFirst");
  }
  else {
    cout << "Divide failed in Draw_Phi_DatasetComparison" << endl;
  }
}

void Draw_Pt_ratio_etaNeg_etaPos_RadiusComparison(int iDataset, float* etaRange) {

  TH3D* H3D_jetRjetPtjetEta;
  TH1D* H1D_jetPt_left[nRadius];
  TH1D* H1D_jetPt_right[nRadius];
  TH1D* H1D_jetPt_left_rebinned[nRadius];
  TH1D* H1D_jetPt_right_rebinned[nRadius];
  TH1D* H1D_jetPt_rebinned_ratios[nRadius];

  H3D_jetRjetPtjetEta = (TH3D*)((TH3D*)file_O2Analysis_list[iDataset]->Get(analysisWorkflow[iDataset]+"/h3_jet_r_jet_pt_jet_eta"+jetFinderQaHistType[iJetFinderQaType]))->Clone("Draw_Pt_ratio_etaNeg_etaPos_RadiusComparison"+Datasets[iDataset]+Form("%.1f", etaRange[0])+"<eta<"+Form("%.1f", etaRange[1]));
  H3D_jetRjetPtjetEta->Sumw2();


  float EtaCutLow = etaRange[0];
  float EtaCutHigh = etaRange[1];
  if (EtaCutLow > 0 || EtaCutHigh < 0 ) {
    cout << "eta=0 should be within the [EtaCutLow,EtaCutLow] range for the Draw_Pt_ratio_etaNeg_etaPos function" << endl;
    return;
  }
  int ibinEta_low = H3D_jetRjetPtjetEta->GetZaxis()->FindBin(EtaCutLow);
  int ibinEta_zero = H3D_jetRjetPtjetEta->GetZaxis()->FindBin(0.);
  int ibinEta_high = H3D_jetRjetPtjetEta->GetZaxis()->FindBin(EtaCutHigh);
    if (ibinEta_low == 0) 
      cout << "WARNING: Pt_ratio_etaNeg_etaPos_RadiusComparison is counting the underflow with the chosen etaRange" << endl;
    if (ibinEta_high == H3D_jetRjetPtjetEta->GetZaxis()->GetNbins()+1) 
      cout << "WARNING: Pt_ratio_etaNeg_etaPos_RadiusComparison is counting the overflow with the chosen etaRange" << endl;
  int ibinJetRadius = 0;

  bool divideSuccess = false;

  for(int iRadius = 0; iRadius < nRadius; iRadius++){
    ibinJetRadius = H3D_jetRjetPtjetEta->GetXaxis()->FindBin(arrayRadius[iRadius]+GLOBAL_epsilon);
    H1D_jetPt_left[iRadius] = (TH1D*)H3D_jetRjetPtjetEta->ProjectionY("jetPt_left_"+Datasets[iDataset]+RadiusLegend[iRadius], ibinJetRadius, ibinJetRadius, ibinEta_low, ibinEta_zero-1, "e");
    H1D_jetPt_right[iRadius] = (TH1D*)H3D_jetRjetPtjetEta->ProjectionY("jetPt_right_"+Datasets[iDataset]+RadiusLegend[iRadius], ibinJetRadius, ibinJetRadius, ibinEta_zero, ibinEta_high, "e");
    H1D_jetPt_left_rebinned[iRadius] = (TH1D*)H1D_jetPt_left[iRadius]->Rebin(1.,"jetPt_left_rebinned_"+Datasets[iDataset]+RadiusLegend[iRadius]+Form("%.1f",arrayRadius[iRadius])+Form("%.1f", etaRange[0])+"<eta<"+Form("%.1f", etaRange[1]));
    H1D_jetPt_right_rebinned[iRadius] = (TH1D*)H1D_jetPt_right[iRadius]->Rebin(1.,"jetPt_right_rebinned_"+Datasets[iDataset]+RadiusLegend[iRadius]+Form("%.1f",arrayRadius[iRadius])+Form("%.1f", etaRange[0])+"<eta<"+Form("%.1f", etaRange[1]));

    // NormaliseYieldToNEntries(H1D_jetPt_left_rebinned[iRadius]);
    // NormaliseYieldToNEntries(H1D_jetPt_right_rebinned[iRadius]);
    NormaliseYieldToNEvents(H1D_jetPt_left_rebinned[iRadius], GetNEventsSelected_JetFramework(file_O2Analysis_list[iDataset]));
    NormaliseYieldToNEvents(H1D_jetPt_right_rebinned[iRadius], GetNEventsSelected_JetFramework(file_O2Analysis_list[iDataset]));

    H1D_jetPt_rebinned_ratios[iRadius] = (TH1D*)H1D_jetPt_left_rebinned[iRadius]->Clone("jetPt_rebinned_ratios"+Datasets[iRadius]+RadiusLegend[iRadius]+Form("%.1f",arrayRadius[iRadius])+Form("%.1f", etaRange[0])+"<eta<"+Form("%.1f", etaRange[1]));
    H1D_jetPt_rebinned_ratios[iRadius]->Reset("M");
    divideSuccess = H1D_jetPt_rebinned_ratios[iRadius]->Divide(H1D_jetPt_right_rebinned[iRadius], H1D_jetPt_left_rebinned[iRadius]);
  }

  TString* pdfName = new TString("jet_"+jetType[iJetType]+"_"+jetLevel[iJetLevel]+"_"+DatasetsNames[iDataset]+"_pT_etaRightLeftRatio_@eta["+Form("%.1f", EtaCutLow)+","+Form("%.1f", EtaCutHigh)+"]"+jetFinderQaHistType[iJetFinderQaType]);

  // TString textContext(contextDatasetRadiusCompAndVarRangecontextDatasetRadiusCompAndVarRange(iDataset, etaRange, "eta"));
  TString textContext(contextCustomThreeFields(*texDatasetsComparisonCommonDenominator, Datasets[iDataset], contextEtaRange(etaRange), ""));

  if (divideSuccess == true) {
    Draw_TH1_Histograms_in_one(H1D_jetPt_rebinned_ratios, RadiusLegend, nRadius, textContext, pdfName, texPtX, texJetRatioEtaComparison, texCollisionDataInfo, drawnWindowAuto, "standardratio");
  }
  else {
    cout << "Divide failed in Draw_Pt_ratio_etaNeg_etaPos_RadiusComparison" << endl;
  }
}

void Draw_Pt_ratio_etaNeg_etaPos_DatasetComparison(float jetRadius, float* etaRange) {

  TH3D* H3D_jetRjetPtjetEta[nDatasets];
  TH1D* H1D_jetPt_left[nDatasets];
  TH1D* H1D_jetPt_right[nDatasets];
  TH1D* H1D_jetPt_left_rebinned[nDatasets];
  TH1D* H1D_jetPt_right_rebinned[nDatasets];
  TH1D* H1D_jetPt_rebinned_ratios[nDatasets];

  float EtaCutLow = etaRange[0];
  float EtaCutHigh = etaRange[1];
  if (EtaCutLow > 0 || EtaCutHigh < 0 ) {
    cout << "eta=0 should be within the [EtaCutLow,EtaCutLow] range for the Draw_Pt_ratio_etaNeg_etaPos function" << endl;
    return;
  }
  int ibinJetRadius = 0;
  int nEvents;

  bool divideSuccess = false;

  for(int iDataset = 0; iDataset < nDatasets; iDataset++){
    H3D_jetRjetPtjetEta[iDataset] = (TH3D*)((TH3D*)file_O2Analysis_list[iDataset]->Get(analysisWorkflow[iDataset]+"/h3_jet_r_jet_pt_jet_eta"+jetFinderQaHistType[iJetFinderQaType]))->Clone("Draw_Pt_ratio_etaNeg_etaPos_DatasetComparison"+Datasets[iDataset]+"Radius"+Form("%.1f",jetRadius)+Form("%.1f", etaRange[0])+"<eta<"+Form("%.1f", etaRange[1]));
    H3D_jetRjetPtjetEta[iDataset]->Sumw2();

    int ibinEta_low = H3D_jetRjetPtjetEta[iDataset]->GetZaxis()->FindBin(EtaCutLow);
    int ibinEta_zero = H3D_jetRjetPtjetEta[iDataset]->GetZaxis()->FindBin(0.);
    int ibinEta_high = H3D_jetRjetPtjetEta[iDataset]->GetZaxis()->FindBin(EtaCutHigh);
    if (ibinEta_low == 0) 
      cout << "WARNING: Pt_ratio_etaNeg_etaPos_DatasetComparison is counting the underflow with the chosen etaRange" << endl;
    if (ibinEta_high == H3D_jetRjetPtjetEta[iDataset]->GetZaxis()->GetNbins()+1) 
      cout << "WARNING: Pt_ratio_etaNeg_etaPos_DatasetComparison is counting the overflow with the chosen etaRange" << endl;

    ibinJetRadius = H3D_jetRjetPtjetEta[iDataset]->GetXaxis()->FindBin(jetRadius+GLOBAL_epsilon);

    H1D_jetPt_left[iDataset] = (TH1D*)H3D_jetRjetPtjetEta[iDataset]->ProjectionY("jetPt_left_"+Datasets[iDataset]+"Radius"+Form("%.1f",jetRadius)+Form("%.1f", etaRange[0])+"<eta<"+Form("%.1f", etaRange[1]),ibinJetRadius,ibinJetRadius, ibinEta_low, ibinEta_zero-1, "e");
    H1D_jetPt_right[iDataset] = (TH1D*)H3D_jetRjetPtjetEta[iDataset]->ProjectionY("jetPt_right_"+Datasets[iDataset]+"Radius"+Form("%.1f",jetRadius)+Form("%.1f", etaRange[0])+"<eta<"+Form("%.1f", etaRange[1]),ibinJetRadius,ibinJetRadius, ibinEta_zero, ibinEta_high, "e");
    H1D_jetPt_left_rebinned[iDataset] = (TH1D*)H1D_jetPt_left[iDataset]->Rebin(1.,"jetPt_left_rebinned_"+Datasets[iDataset]+"Radius"+Form("%.1f",jetRadius)+Form("%.1f", etaRange[0])+"<eta<"+Form("%.1f", etaRange[1]));
    H1D_jetPt_right_rebinned[iDataset] = (TH1D*)H1D_jetPt_right[iDataset]->Rebin(1.,"jetPt_right_rebinned_"+Datasets[iDataset]+"Radius"+Form("%.1f",jetRadius)+Form("%.1f", etaRange[0])+"<eta<"+Form("%.1f", etaRange[1]));
    
    // NormaliseYieldToNEntries(H1D_jetPt_left_rebinned[iDataset]);
    // NormaliseYieldToNEntries(H1D_jetPt_right_rebinned[iDataset]);
    NormaliseYieldToNEvents(H1D_jetPt_left_rebinned[iDataset], GetNEventsSelected_JetFramework(file_O2Analysis_list[iDataset]));
    NormaliseYieldToNEvents(H1D_jetPt_right_rebinned[iDataset],  GetNEventsSelected_JetFramework(file_O2Analysis_list[iDataset]));

    H1D_jetPt_rebinned_ratios[iDataset] = (TH1D*)H1D_jetPt_left_rebinned[iDataset]->Clone("jetPt_rebinned_ratios"+Datasets[iDataset]+"Radius"+Form("%.1f",jetRadius)+Form("%.1f", etaRange[0])+"<eta<"+Form("%.1f", etaRange[1]));
    H1D_jetPt_rebinned_ratios[iDataset]->Reset("M");
    divideSuccess = H1D_jetPt_rebinned_ratios[iDataset]->Divide(H1D_jetPt_right_rebinned[iDataset], H1D_jetPt_left_rebinned[iDataset]);
  }
 
  TString* pdfName = new TString("jet_"+jetType[iJetType]+"_"+jetLevel[iJetLevel]+"_DataComp_R="+Form("%.1f", jetRadius)+"_pT_etaRightLeftRatio_@eta["+Form("%.1f", EtaCutLow)+","+Form("%.1f", EtaCutHigh)+"]"+jetFinderQaHistType[iJetFinderQaType]);

  // TString textContext(contextDatasetCompAndRadiusAndVarRange(jetRadius, etaRange, "eta"));
  TString textContext(contextCustomThreeFields(*texDatasetsComparisonCommonDenominator, "", "#splitline{"+contextJetRadius(jetRadius)+"}{"+contextEtaRange(etaRange)+"}", ""));

  if (divideSuccess == true) {
    Draw_TH1_Histograms_in_one(H1D_jetPt_rebinned_ratios, DatasetsNames, nDatasets, textContext, pdfName, texPtX, texJetRatioEtaComparison, texCollisionDataInfo, drawnWindowAuto, "standardratio");
  }
  else {
    cout << "Divide failed in Draw_Pt_ratio_etaNeg_etaPos_DatasetComparison" << endl;
  }
}


void Draw_Area_PtIntegrated_DatasetComparison(float jetRadius, float* PtRange) {

  TH3D* H3D_jetRjetPtjetArea[nDatasets];
  TH1D* H1D_jetArea[nDatasets];
  TH1D* H1D_jetArea_rebinned[nDatasets];
  
  TH1D* H1D_jetArea_rebinned_ratios[nDatasets];

  float PtCutLow = PtRange[0];
  float PtCutHigh = PtRange[1];
  int ibinJetRadius = 0;

  bool divideSuccess = false;

  for(int iDataset = 0; iDataset < nDatasets; iDataset++){

    H3D_jetRjetPtjetArea[iDataset] = (TH3D*)((TH3D*)file_O2Analysis_list[iDataset]->Get(analysisWorkflow[iDataset]+"/h3_jet_r_jet_pt_jet_area"+jetFinderQaHistType[iJetFinderQaType]))->Clone("Draw_Area_PtIntegrated_DatasetComparison"+Datasets[iDataset]+"Radius"+Form("%.1f",jetRadius)+Form("%.1f", PtRange[0])+"<pt<"+Form("%.1f", PtRange[1]));

    int ibinPt_low = H3D_jetRjetPtjetArea[iDataset]->GetYaxis()->FindBin(PtCutLow);
    int ibinPt_high = H3D_jetRjetPtjetArea[iDataset]->GetYaxis()->FindBin(PtCutHigh);
    if (ibinPt_low == 0) 
      cout << "WARNING: Area_DatasetComparison is counting the underflow with the chosen PtRange" << endl;
    if (ibinPt_high == H3D_jetRjetPtjetArea[iDataset]->GetYaxis()->GetNbins()+1) 
      cout << "WARNING: Area_DatasetComparison is counting the overflow with the chosen PtRange" << endl;
    ibinJetRadius = H3D_jetRjetPtjetArea[iDataset]->GetXaxis()->FindBin(jetRadius+GLOBAL_epsilon);

    H1D_jetArea[iDataset] = (TH1D*)H3D_jetRjetPtjetArea[iDataset]->ProjectionZ("jetArea_"+Datasets[iDataset]+"Radius"+Form("%.1f",jetRadius)+Form("%.1f", PtCutLow)+"<pt<"+Form("%.1f", PtCutHigh), ibinJetRadius,ibinJetRadius, ibinPt_low, ibinPt_high, "e");
    H1D_jetArea[iDataset]->Sumw2();
    H1D_jetArea_rebinned[iDataset] = (TH1D*)H1D_jetArea[iDataset]->Rebin(1.,"jetArea_rebinned_"+Datasets[iDataset]+"Radius"+Form("%.1f",jetRadius)+Form("%.1f", PtCutLow)+"<pt<"+Form("%.1f", PtCutHigh));

    // NormaliseYieldToNEntries(H1D_jetArea_rebinned[iDataset]);
    NormaliseYieldToNEvents(H1D_jetArea_rebinned[iDataset], GetNEventsSelected_JetFramework(file_O2Analysis_list[iDataset]));

    H1D_jetArea_rebinned_ratios[iDataset] = (TH1D*)H1D_jetArea_rebinned[iDataset]->Clone("jetArea_rebinned_ratios"+Datasets[iDataset]+"Radius"+Form("%.1f",jetRadius)+Form("%.1f", PtCutLow)+"<pt<"+Form("%.1f", PtCutHigh));
    H1D_jetArea_rebinned_ratios[iDataset]->Reset("M");
    divideSuccess = H1D_jetArea_rebinned_ratios[iDataset]->Divide(H1D_jetArea_rebinned[iDataset], H1D_jetArea_rebinned[0]);
  }

  TString* pdfName = new TString("jet_"+jetType[iJetType]+"_"+jetLevel[iJetLevel]+"_DataComp_R="+Form("%.1f", jetRadius)+"_Area_@pT["+Form("%03.0f", PtCutLow)+","+Form("%03.0f", PtCutHigh)+"]"+jetFinderQaHistType[iJetFinderQaType]);
  TString* pdfName_ratio = new TString("jet_"+jetType[iJetType]+"_"+jetLevel[iJetLevel]+"_DataComp_R="+Form("%.1f", jetRadius)+"_Area_@pT["+Form("%03.0f", PtCutLow)+","+Form("%03.0f", PtCutHigh)+"]"+jetFinderQaHistType[iJetFinderQaType]+"_ratio");

  // TString textContext(contextDatasetCompAndRadiusAndVarRange(jetRadius, PtRange, "pt"));
  TString textContext(contextCustomThreeFields(*texDatasetsComparisonCommonDenominator, "", "#splitline{"+contextJetRadius(jetRadius)+"}{"+contextPtRange(PtRange)+"}", ""));

  Draw_TH1_Histograms_in_one(H1D_jetArea_rebinned, DatasetsNames, nDatasets, textContext, pdfName, texJetArea, texJetNormAreaYield, texCollisionDataInfo, drawnWindowAuto, "");

  if (divideSuccess == true) {
    Draw_TH1_Histograms_in_one(H1D_jetArea_rebinned_ratios, DatasetsNames, nDatasets, textContext, pdfName_ratio, texJetArea, texRatioDatasets, texCollisionDataInfo, drawnWindowAuto, "standardratio,avoidFirst");
  }
  else {
    cout << "Divide failed in Draw_Area_PtIntegrated_DatasetComparison" << endl;
  }
}


void Draw_JetTRDratio_vs_JetEta(int iDataset) {
  const int nHistCollection = 1;
  const TString histCollectionLegend[nHistCollection] = {""};

  TH3D* H3D_jetPtjetEtajetTRDratio;
  TH2D* H2D_jetEtajetTRDratio[nHistCollection];
  TH2D* H2D_jetEtajetTRDratio_rebinned[nHistCollection];
  
  H3D_jetPtjetEtajetTRDratio = (TH3D*)((TH3D*)file_O2Analysis_list[iDataset]->Get(analysisWorkflow[iDataset]+"/h3_jet_pt_jet_eta_jet_nTRDtracks_ratio"+jetFinderQaHistType[iJetFinderQaType]))->Clone("Draw_JetTRDratio_vs_JetEta"+Datasets[iDataset]);
  H3D_jetPtjetEtajetTRDratio->Sumw2();

  H2D_jetEtajetTRDratio[0] = (TH2D*)H3D_jetPtjetEtajetTRDratio->Project3D((TString)Form("%d", iDataset)+"_jetEta_jetTRratio_e_zy"); //can't use letter D in this or it seems to replace the histogram in current pad (see documentation of ProjectionX function. Isn't mentioned in project3D sadly)
  // H2D_jetEtajetTRDratio_rebinned = (TH2D*)H2D_jetEtajetTRDratio->RebinY(2.,"H1D_jetEta_jetTRD_rebinned_"+RadiusLegend[iRadius]);

  TString* pdfName = new TString("jet_"+jetType[iJetType]+"_"+jetLevel[iJetLevel]+"_"+DatasetsNames[iDataset]+"_JetTRDratio-vs-Eta"+jetFinderQaHistType[iJetFinderQaType]);

  TString textContext("#splitline{"+*texDatasetsComparisonCommonDenominator+" "+DatasetsNames[iDataset]+"}{2023 QC}");

  Draw_TH2_Histograms(H2D_jetEtajetTRDratio, histCollectionLegend, nHistCollection, textContext, pdfName, texEtaX, texJetNTrdTracksRatio, texCollisionDataInfo, drawnWindowAuto, "");
}

void Draw_JetTRDcount_vs_JetEta(int iDataset) {
  const int nHistCollection = 1;
  const TString histCollectionLegend[nHistCollection] = {""};

  TH3D* H3D_jetPtjetEtajetTRD;
  TH2D* H2D_jetEtajetTRD[nHistCollection];
  TH2D* H2D_jetEtajetTRD_rebinned[nHistCollection];
  
  H3D_jetPtjetEtajetTRD = (TH3D*)((TH3D*)file_O2Analysis_list[iDataset]->Get(analysisWorkflow[iDataset]+"/h3_jet_pt_jet_eta_jet_nTRDtracks"+jetFinderQaHistType[iJetFinderQaType]))->Clone("Draw_JetTRDcount_vs_JetEta"+Datasets[iDataset]);
  H3D_jetPtjetEtajetTRD->Sumw2();

  H2D_jetEtajetTRD[0] = (TH2D*)H3D_jetPtjetEtajetTRD->Project3D((TString)Form("%d", iDataset)+"_jetEta_jetTR_e_zy"); //can't use letter D in this or it seems to replace the histogram in current pad (see documentation of ProjectionX function. Isn't mentioned in project3D sadly)
  // // H2D_jetEtajetTRD_rebinned = (TH2D*)H2D_jetEtajetTRD->RebinY(2.,"H1D_jetEta_jetTRD_rebinned_"+RadiusLegend[iRadius]);

  TString* pdfName = new TString("jet_"+jetType[iJetType]+"_"+jetLevel[iJetLevel]+"_"+DatasetsNames[iDataset]+"_JetTRD-vs-Eta"+jetFinderQaHistType[iJetFinderQaType]);

  TString textContext("#splitline{"+*texDatasetsComparisonCommonDenominator+" "+DatasetsNames[iDataset]+"}{2023 QC}");

  Draw_TH2_Histograms(H2D_jetEtajetTRD, histCollectionLegend, nHistCollection, textContext, pdfName, texEtaX, texJetNTrdTracks, texCollisionDataInfo, drawnWindowAuto, "");
}

void Draw_Pt_ratio_etaNeg_etaPos_TRDonly_vs_noTRD(int iDataset, float* etaRange) {
  const int nHistCollection = 1;
  const TString histCollectionLegend[nHistCollection] = {""};

  TH3D* H3D_jetPtjetEtajetTRD;
  TH1D* H1D_jetPt_left_TRDonly[nHistCollection];
  TH1D* H1D_jetPt_right_TRDonly[nHistCollection];
  TH1D* H1D_jetPt_left_TRDonly_rebinned[nHistCollection];
  TH1D* H1D_jetPt_right_TRDonly_rebinned[nHistCollection];
  TH1D* H1D_jetPt_TRDonly_rebinned_ratios[nHistCollection];

  H3D_jetPtjetEtajetTRD = (TH3D*)((TH3D*)file_O2Analysis_list[iDataset]->Get(analysisWorkflow[iDataset]+"/h3_jet_pt_jet_eta_jet_nTRDtracks"+jetFinderQaHistType[iJetFinderQaType]))->Clone("Draw_Pt_ratio_etaNeg_etaPos_JetTRDonly"+Datasets[iDataset]+Form("%.1f", etaRange[0])+"<eta<"+Form("%.1f", etaRange[1]));
  H3D_jetPtjetEtajetTRD->Sumw2();

  TH1D* H1D_jetPt_left_noTRD[nHistCollection];
  TH1D* H1D_jetPt_right_noTRD[nHistCollection];
  TH1D* H1D_jetPt_left_noTRD_rebinned[nHistCollection];
  TH1D* H1D_jetPt_right_noTRD_rebinned[nHistCollection];
  TH1D* H1D_jetPt_noTRD_rebinned_ratios[nHistCollection];

  H3D_jetPtjetEtajetTRD = (TH3D*)((TH3D*)file_O2Analysis_list[iDataset]->Get(analysisWorkflow[iDataset]+"/h3_jet_pt_jet_eta_jet_nTRDtracks"+jetFinderQaHistType[iJetFinderQaType]))->Clone("Draw_Pt_ratio_etaNeg_etaPos_JetNoTRD"+Datasets[iDataset]+Form("%.1f", etaRange[0])+"<eta<"+Form("%.1f", etaRange[1]));
  H3D_jetPtjetEtajetTRD->Sumw2();

  TH1D* H1D_jetPt_noTRD_vs_onlyTRD_ratios[nHistCollection];

  float EtaCutLow = etaRange[0];
  float EtaCutHigh = etaRange[1];
  if (EtaCutLow > 0 || EtaCutHigh < 0 ) {
    cout << "eta=0 should be within the [EtaCutLow,EtaCutLow] range for the Draw_Pt_ratio_etaNeg_etaPos function" << endl;
    return;
  }
  int ibinEta_low = H3D_jetPtjetEtajetTRD->GetYaxis()->FindBin(EtaCutLow);
  int ibinEta_zero = H3D_jetPtjetEtajetTRD->GetYaxis()->FindBin(0.);
  int ibinEta_high = H3D_jetPtjetEtajetTRD->GetYaxis()->FindBin(EtaCutHigh);
    if (ibinEta_low == 0) 
      cout << "WARNING: Pt_ratio_etaNeg_etaPos_RadiusComparison is counting the underflow with the chosen etaRange" << endl;
    if (ibinEta_high == H3D_jetPtjetEtajetTRD->GetYaxis()->GetNbins()+1) 
      cout << "WARNING: Pt_ratio_etaNeg_etaPos_RadiusComparison is counting the overflow with the chosen etaRange" << endl;
  int ibinJetRadius = 0;

  bool divideSuccess = false;


  //// jets with TRDmatched tracks only 
  H1D_jetPt_left_TRDonly[0] = (TH1D*)H3D_jetPtjetEtajetTRD->ProjectionX("jetPt_TRDonly_left"+Datasets[iDataset]+Form("%.1f", etaRange[0])+"<eta<"+Form("%.1f", etaRange[1]), ibinEta_low, ibinEta_zero-1, 2, -1, "e"); // asking bin number >= 2 avoids TRD match == 0; it's good to also keep overflow bin with -1 for upper limit as they do have TRD matches
  H1D_jetPt_right_TRDonly[0] = (TH1D*)H3D_jetPtjetEtajetTRD->ProjectionX("jetPt_TRDonly_right"+Datasets[iDataset]+Form("%.1f", etaRange[0])+"<eta<"+Form("%.1f", etaRange[1]), ibinEta_zero, ibinEta_high, 2, -1, "e"); // asking bin number >= 2 avoids TRD match == 0; it's good to also keep overflow bin with -1 for upper limit as they do have TRD matches
  H1D_jetPt_left_TRDonly_rebinned[0] = (TH1D*)H1D_jetPt_left_TRDonly[0]->Rebin(1.,"jetPt_TRDonly_left_rebinned"+Datasets[iDataset]+Form("%.1f", etaRange[0])+"<eta<"+Form("%.1f", etaRange[1]));
  H1D_jetPt_right_TRDonly_rebinned[0] = (TH1D*)H1D_jetPt_right_TRDonly[0]->Rebin(1.,"jetPt_TRDonly_right_rebinned"+Datasets[iDataset]+Form("%.1f", etaRange[0])+"<eta<"+Form("%.1f", etaRange[1]));

  NormaliseYieldToNEvents(H1D_jetPt_left_TRDonly_rebinned[0], GetNEventsSelected_JetFramework(file_O2Analysis_list[iDataset]));
  NormaliseYieldToNEvents(H1D_jetPt_right_TRDonly_rebinned[0], GetNEventsSelected_JetFramework(file_O2Analysis_list[iDataset]));

  H1D_jetPt_TRDonly_rebinned_ratios[0] = (TH1D*)H1D_jetPt_left_TRDonly_rebinned[0]->Clone("jetPt_TRDonly_rebinned_ratios"+Datasets[iDataset]+Form("%.1f", etaRange[0])+"<eta<"+Form("%.1f", etaRange[1]));
  H1D_jetPt_TRDonly_rebinned_ratios[0]->Reset("M");
  divideSuccess = H1D_jetPt_TRDonly_rebinned_ratios[0]->Divide(H1D_jetPt_right_TRDonly_rebinned[0], H1D_jetPt_left_TRDonly_rebinned[0]);

  TString* pdfNameA = new TString("jet_"+jetType[iJetType]+"_"+jetLevel[iJetLevel]+"_"+DatasetsNames[iDataset]+"_pT_etaRightLeftRatio_JetTRDonly"+jetFinderQaHistType[iJetFinderQaType]);

  // TString textContextA(contextDatasetRadiusCompAndVarRange(iDataset, etaRange, "eta"));
  TString textContextA(contextCustomThreeFields(*texDatasetsComparisonCommonDenominator, Datasets[iDataset], contextEtaRange(etaRange), ""));

  if (divideSuccess == true) {
    Draw_TH1_Histograms_in_one(H1D_jetPt_TRDonly_rebinned_ratios, histCollectionLegend, nHistCollection, textContextA, pdfNameA, texPtX, texJetRatioEtaComparison, texCollisionDataInfo, drawnWindowAuto, "standardratio");
  }
  else {
    cout << "Divide failed in Draw_Pt_ratio_etaNeg_etaPos_RadiusComparison 1" << endl;
  }


  divideSuccess = false;
  //// jets without TRDmatched tracks only 

  H1D_jetPt_left_noTRD[0] = (TH1D*)H3D_jetPtjetEtajetTRD->ProjectionX("jetPt_noTRD_left"+Datasets[iDataset]+Form("%.1f", etaRange[0])+"<eta<"+Form("%.1f", etaRange[1]), ibinEta_low, ibinEta_zero-1, 0, 1, "e"); // asking bin number >= 2 avoids TRD match == 0; it's good to also keep overflow bin with -1 for upper limit as they do have TRD matches
  H1D_jetPt_right_noTRD[0] = (TH1D*)H3D_jetPtjetEtajetTRD->ProjectionX("jetPt_noTRD_right"+Datasets[iDataset]+Form("%.1f", etaRange[0])+"<eta<"+Form("%.1f", etaRange[1]), ibinEta_zero, ibinEta_high, 0, 1, "e"); // asking bin number >= 2 avoids TRD match == 0; it's good to also keep overflow bin with -1 for upper limit as they do have TRD matches
  H1D_jetPt_left_noTRD_rebinned[0] = (TH1D*)H1D_jetPt_left_noTRD[0]->Rebin(1.,"jetPt_noTRD_left_rebinned"+Datasets[iDataset]+Form("%.1f", etaRange[0])+"<eta<"+Form("%.1f", etaRange[1]));
  H1D_jetPt_right_noTRD_rebinned[0] = (TH1D*)H1D_jetPt_right_noTRD[0]->Rebin(1.,"jetPt_noTRD_right_rebinned"+Datasets[iDataset]+Form("%.1f", etaRange[0])+"<eta<"+Form("%.1f", etaRange[1]));

  NormaliseYieldToNEvents(H1D_jetPt_left_noTRD_rebinned[0], GetNEventsSelected_JetFramework(file_O2Analysis_list[iDataset]));
  NormaliseYieldToNEvents(H1D_jetPt_right_noTRD_rebinned[0], GetNEventsSelected_JetFramework(file_O2Analysis_list[iDataset]));

  H1D_jetPt_noTRD_rebinned_ratios[0] = (TH1D*)H1D_jetPt_left_noTRD_rebinned[0]->Clone("jetPt_noTRD_rebinned_ratios"+Datasets[iDataset]+Form("%.1f", etaRange[0])+"<eta<"+Form("%.1f", etaRange[1]));
  H1D_jetPt_noTRD_rebinned_ratios[0]->Reset("M");
  divideSuccess = H1D_jetPt_noTRD_rebinned_ratios[0]->Divide(H1D_jetPt_right_noTRD_rebinned[0], H1D_jetPt_left_noTRD_rebinned[0]);

  TString* pdfNameB = new TString("jet_"+jetType[iJetType]+"_"+jetLevel[iJetLevel]+"_"+DatasetsNames[iDataset]+"_pT_etaRightLeftRatio_JetNoTRD"+jetFinderQaHistType[iJetFinderQaType]);

  // TString textContextB(contextDatasetRadiusCompAndVarRange(iDataset, etaRange, "eta"));
  TString textContextB(contextCustomThreeFields(*texDatasetsComparisonCommonDenominator, Datasets[iDataset], contextEtaRange(etaRange), ""));

  if (divideSuccess == true) {
    Draw_TH1_Histograms_in_one(H1D_jetPt_noTRD_rebinned_ratios, histCollectionLegend, nHistCollection, textContextB, pdfNameB, texPtX, texJetRatioEtaComparison, texCollisionDataInfo, drawnWindowAuto, "standardratio");
  }
  else {
    cout << "Divide failed in Draw_Pt_ratio_etaNeg_etaPos_RadiusComparison 2" << endl;
  }


  //// ratio between no trd and only trd
  H1D_jetPt_noTRD_vs_onlyTRD_ratios[0] = (TH1D*)H1D_jetPt_left_noTRD_rebinned[0]->Clone("jetPt_noTRD_vs_onlyTRD_ratios"+Datasets[iDataset]+Form("%.1f", etaRange[0])+"<eta<"+Form("%.1f", etaRange[1]));
  H1D_jetPt_noTRD_vs_onlyTRD_ratios[0]->Reset("M");
  divideSuccess = H1D_jetPt_noTRD_vs_onlyTRD_ratios[0]->Divide(H1D_jetPt_TRDonly_rebinned_ratios[0], H1D_jetPt_noTRD_rebinned_ratios[0]);

  TString* pdfNameC = new TString("jet_"+jetType[iJetType]+"_"+jetLevel[iJetLevel]+"_"+DatasetsNames[iDataset]+"_pT_etaRightLeftRatio_Jet_NoTRD_vs_OnlyTRD_ratio"+jetFinderQaHistType[iJetFinderQaType]);

  // TString textContextC(contextDatasetRadiusCompAndVarRange(iDataset, etaRange, "eta"));
  TString textContextC(contextCustomThreeFields(*texDatasetsComparisonCommonDenominator, Datasets[iDataset], contextEtaRange(etaRange), ""));

  if (divideSuccess == true) {
    Draw_TH1_Histograms_in_one(H1D_jetPt_noTRD_vs_onlyTRD_ratios, histCollectionLegend, nHistCollection, textContextC, pdfNameC, texPtX, texJetRatioTRDvsNoTRD, texCollisionDataInfo, drawnWindowAuto, "standardratio");
  }
  else {
    cout << "Divide failed in Draw_Pt_ratio_etaNeg_etaPos_RadiusComparison 3" << endl;
  }

}


// To-do list:
// - implement carolina's macro of the gpad thing to automate the division of canvas based on how many plots one wants ?
// - 1D graphs: 2023 QC displayed on top of texCollisionDataInfo; but 2D graphs fine --> doesnt like the nested splitline


// option e for projection computes errors

void Draw_Area_PtIntegrated_BinningComparison(int iDataset, float jetRadius, float* PtRange) {


  // int nRebinValues = 2;
  // int rebinValue[2] = {1, 20};
  // const TString rebinValuesNames[2] = {"rebin(1)", "rebin(20)"};
  // const TString rebinValuesTitle = {"rebin(1,20)"};

  int nRebinValues = 1;
  int rebinValue[1] = {1};
  const TString rebinValuesNames[1] = {"rebin(1)"};
  const TString rebinValuesTitle = {"rebin(1)"};

  TH3D* H3D_jetRjetPtjetArea;
  TH1D* H1D_jetArea;
  TH1D* H1D_jetArea_rebinned[nRebinValues];
  
  TH1D* H1D_jetArea_rebinned_ratios[nRebinValues];

  float PtCutLow = PtRange[0];
  float PtCutHigh = PtRange[1];
  int ibinJetRadius = 0;

  bool divideSuccess = false;

  H3D_jetRjetPtjetArea = (TH3D*)((TH3D*)file_O2Analysis_list[iDataset]->Get(analysisWorkflow[iDataset]+"/h3_jet_r_jet_pt_jet_area"+jetFinderQaHistType[iJetFinderQaType]))->Clone("Draw_Area_PtIntegrated_DatasetComparison"+Datasets[iDataset]+"Radius"+Form("%.1f",jetRadius)+Form("%.1f", PtRange[0])+"<pt<"+Form("%.1f", PtRange[1]));

  int ibinPt_low = H3D_jetRjetPtjetArea->GetYaxis()->FindBin(PtCutLow);
  int ibinPt_high = H3D_jetRjetPtjetArea->GetYaxis()->FindBin(PtCutHigh);
  if (ibinPt_low == 0) 
    cout << "WARNING: Area_DatasetComparison is counting the underflow with the chosen PtRange" << endl;
  if (ibinPt_high == H3D_jetRjetPtjetArea->GetYaxis()->GetNbins()+1) 
    cout << "WARNING: Area_DatasetComparison is counting the overflow with the chosen PtRange" << endl;
  ibinJetRadius = H3D_jetRjetPtjetArea->GetXaxis()->FindBin(jetRadius+GLOBAL_epsilon);

  H1D_jetArea = (TH1D*)H3D_jetRjetPtjetArea->ProjectionZ("jetArea_"+Datasets[iDataset]+"Radius"+Form("%.1f",jetRadius)+Form("%.1f", PtCutLow)+"<pt<"+Form("%.1f", PtCutHigh), ibinJetRadius,ibinJetRadius, ibinPt_low, ibinPt_high, "e");
  H1D_jetArea->Sumw2();

  for(int iRebinValue = 0; iRebinValue < nRebinValues; iRebinValue++){
    H1D_jetArea_rebinned[iRebinValue] = (TH1D*)H1D_jetArea->Rebin(rebinValue[iRebinValue],"jetArea_rebinned_"+Datasets[iDataset]+"rebin"+rebinValue[iRebinValue]+"Radius"+Form("%.1f",jetRadius)+Form("%.1f", PtCutLow)+"<pt<"+Form("%.1f", PtCutHigh));

    NormaliseYieldToNEvents(H1D_jetArea_rebinned[iRebinValue], GetNEventsSelected_JetFramework(file_O2Analysis_list[iDataset]));

    H1D_jetArea_rebinned_ratios[iRebinValue] = (TH1D*)H1D_jetArea_rebinned[iRebinValue]->Clone("jetArea_rebinned1_ratios"+Datasets[iDataset]+"rebin"+rebinValue[iRebinValue]+"Radius"+Form("%.1f",jetRadius)+Form("%.1f", PtCutLow)+"<pt<"+Form("%.1f", PtCutHigh));
    H1D_jetArea_rebinned_ratios[iRebinValue]->Reset("M");
    divideSuccess = H1D_jetArea_rebinned_ratios[iRebinValue]->Divide(H1D_jetArea_rebinned[iRebinValue], H1D_jetArea_rebinned[0]);
  }
  // for (int i = 0; i<300; i++){
  //   H1D_jetArea_rebinned[nRebinValues+1]->Fill(i*0.005, 1);
  // }

  TString* pdfName = new TString("jet_"+jetType[iJetType]+"_"+jetLevel[iJetLevel]+"_"+rebinValuesTitle+"_R="+Form("%.1f", jetRadius)+"_Area_@pT["+Form("%03.0f", PtCutLow)+","+Form("%03.0f", PtCutHigh)+"]"+jetFinderQaHistType[iJetFinderQaType]);
  TString* pdfName_ratio = new TString("jet_"+jetType[iJetType]+"_"+jetLevel[iJetLevel]+"_"+rebinValuesTitle+"_R="+Form("%.1f", jetRadius)+"_Area_@pT["+Form("%03.0f", PtCutLow)+","+Form("%03.0f", PtCutHigh)+"]"+jetFinderQaHistType[iJetFinderQaType]+"_ratio");

  // TString textContext(contextDatasetCompAndRadiusAndVarRange(jetRadius, PtRange, "pt"));
  TString textContext(contextCustomThreeFields(*texDatasetsComparisonCommonDenominator, "", "#splitline{"+contextJetRadius(jetRadius)+"}{"+contextPtRange(PtRange)+"}", ""));

  Draw_TH1_Histograms_in_one(H1D_jetArea_rebinned, rebinValuesNames, nRebinValues, textContext, pdfName, texJetArea, texJetNormAreaYield, texCollisionDataInfo, drawnWindowAuto, "");

  if (divideSuccess == true) {
    Draw_TH1_Histograms_in_one(H1D_jetArea_rebinned_ratios, rebinValuesNames, nRebinValues, textContext, pdfName_ratio, texJetArea, texRatioDatasets, texCollisionDataInfo, drawnWindowAuto, "autoratio,avoidFirst");
  }
  else {
    cout << "Divide failed in Draw_Area_PtIntegrated_BinningComparison" << endl;
  }



  TH1D* H1D_JetAreaBinSpacing[1] = {new TH1D("H1D_JetAreaBinSpacing", "H1D_JetAreaBinSpacing", 200, 0.004, 0.006)};

  int nonEmptyBinsCount = 0;
  double previousBinCenter = 0;
  for (int iAreaBin = 0; iAreaBin < H1D_jetArea->GetNbinsX(); iAreaBin++){
    if (H1D_jetArea->GetBinContent(iAreaBin)>0) {
      nonEmptyBinsCount++;
      cout << "AIMERIC ------------------- Non-empty Area bins: " << H1D_jetArea->GetBinCenter(iAreaBin) << endl;
      H1D_JetAreaBinSpacing[0]->Fill(H1D_jetArea->GetBinCenter(iAreaBin) - previousBinCenter);
      previousBinCenter = H1D_jetArea->GetBinCenter(iAreaBin);
    }
  }


  int jetAreaVsNGhost_nCurves = 2;
  const TString jetAreaVsNGhost_CurveNames[2] = {"O2Physics distribution", "nGhost_x_GhostArea"};
  TH1D* H1D_jetArea_vs_nGhost[2] = {new TH1D("H1D_jetArea_vs_nGhost1", "H1D_jetArea_vs_nGhost1", nonEmptyBinsCount, 0, nonEmptyBinsCount), new TH1D("H1D_jetArea_vs_nGhost2", "H1D_jetArea_vs_nGhost2", nonEmptyBinsCount, 0, nonEmptyBinsCount)};

  nonEmptyBinsCount = 0; //reset for the loop
  float ghostArea = 0.005;
  for (int iAreaBin = 0; iAreaBin < H1D_jetArea->GetNbinsX(); iAreaBin++){
    if (H1D_jetArea->GetBinContent(iAreaBin)>0) {
      nonEmptyBinsCount++;
      H1D_jetArea_vs_nGhost[0]->SetBinContent(nonEmptyBinsCount, H1D_jetArea->GetBinCenter(iAreaBin));
      H1D_jetArea_vs_nGhost[1]->SetBinContent(nonEmptyBinsCount, (nonEmptyBinsCount+3) * ghostArea);
    }
  }


  TString* pdfName_JetAreaBinSpacing = new TString("jet_"+jetType[iJetType]+"_"+jetLevel[iJetLevel]+"_R="+Form("%.1f", jetRadius)+"_AreaBinSpacing_@pT["+Form("%03.0f", PtCutLow)+","+Form("%03.0f", PtCutHigh)+"]"+jetFinderQaHistType[iJetFinderQaType]);
  TString* pdfName_jetArea_vs_nGhost = new TString("jet_"+jetType[iJetType]+"_"+jetLevel[iJetLevel]+"_R="+Form("%.1f", jetRadius)+"_Area_vs_nGhost_@pT["+Form("%03.0f", PtCutLow)+","+Form("%03.0f", PtCutHigh)+"]"+jetFinderQaHistType[iJetFinderQaType]);

  Draw_TH1_Histograms_in_one(H1D_JetAreaBinSpacing, rebinValuesNames, nRebinValues, textContext, pdfName_JetAreaBinSpacing, texBinSpacing, texCount, texCollisionDataInfo, drawnWindowAuto, "");
  Draw_TH1_Histograms_in_one(H1D_jetArea_vs_nGhost, jetAreaVsNGhost_CurveNames, jetAreaVsNGhost_nCurves, textContext, pdfName_jetArea_vs_nGhost, texNGhost, texJetArea, texCollisionDataInfo, drawnWindowAuto, "");


  TH1D* H1D_jetArea_vs_nGhost_ratios[2] = {new TH1D("H1D_jetArea_vs_nGhost1_ratio", "H1D_jetArea_vs_nGhost1_ratio", nonEmptyBinsCount, 0, nonEmptyBinsCount), new TH1D("H1D_jetArea_vs_nGhost2_ratio", "H1D_jetArea_vs_nGhost2_ratio", nonEmptyBinsCount, 0, nonEmptyBinsCount)};

  H1D_jetArea_vs_nGhost_ratios[0] = (TH1D*)H1D_jetArea_vs_nGhost[0]->Clone("jetArea_rebinned0_ratios"+Datasets[iDataset]+"Radius"+Form("%.1f",jetRadius)+Form("%.1f", PtCutLow)+"<pt<"+Form("%.1f", PtCutHigh));
  H1D_jetArea_vs_nGhost_ratios[0]->Reset("M");
  H1D_jetArea_vs_nGhost_ratios[1] = (TH1D*)H1D_jetArea_vs_nGhost[1]->Clone("jetArea_rebinned1_ratios"+Datasets[iDataset]+"Radius"+Form("%.1f",jetRadius)+Form("%.1f", PtCutLow)+"<pt<"+Form("%.1f", PtCutHigh));
  H1D_jetArea_vs_nGhost_ratios[1]->Reset("M");
  divideSuccess = H1D_jetArea_vs_nGhost_ratios[0]->Divide(H1D_jetArea_vs_nGhost[0], H1D_jetArea_vs_nGhost[0]);
  divideSuccess = H1D_jetArea_vs_nGhost_ratios[1]->Divide(H1D_jetArea_vs_nGhost[1], H1D_jetArea_vs_nGhost[0]);

  TString* pdfName_jetArea_vs_nGhost_ratio = new TString("jet_"+jetType[iJetType]+"_"+jetLevel[iJetLevel]+"_R="+Form("%.1f", jetRadius)+"_Area_vs_nGhost_@pT["+Form("%03.0f", PtCutLow)+","+Form("%03.0f", PtCutHigh)+"]"+jetFinderQaHistType[iJetFinderQaType]+"_ratio");

  if (divideSuccess == true) {
    Draw_TH1_Histograms_in_one(H1D_jetArea_vs_nGhost_ratios, jetAreaVsNGhost_CurveNames, 2, textContext, pdfName_ratio, texNGhost, texJetRatioAreaJetVsNGhost, texCollisionDataInfo, drawnWindowAuto, "autoratio,avoidFirst");
  }
  else {
    cout << "Divide failed in Draw_Area_PtIntegrated_BinningComparison" << endl;
  }
}











void Draw_Rho_vs_Centrality_DatasetComp() {

  TH2D* H2D_rhoCentrality[nDatasets];
  TH2D* H2D_rhoCentrality_rebinned[nDatasets];

  for(int iDataset = 0; iDataset < nDatasets; iDataset++){
    H2D_rhoCentrality[iDataset] = (TH2D*)((TH2D*)file_O2Analysis_list[iDataset]->Get(analysisWorkflow[iDataset]+"/h2_centrality_rho"))->Clone("Draw_Rho_vs_Centrality"+Datasets[iDataset]);
    H2D_rhoCentrality[iDataset]->Sumw2();

    H2D_rhoCentrality_rebinned[iDataset] = (TH2D*)H2D_rhoCentrality[iDataset]->RebinX(10,"H2D_rhoCentrality_rebinned"+Datasets[iDataset]);
    H2D_rhoCentrality_rebinned[iDataset]->GetXaxis()->SetRange(0,H2D_rhoCentrality_rebinned[iDataset]->GetNbinsX()-10);
    // H2D_rhoCentrality_rebinned[iDataset]->GetYaxis()->SetRange(1, H2D_rhoCentrality_rebinned[iDataset]->FindLastBinAbove(1, 2)); //(asks for the last bin on the y axis (axis number 2) to have strictly more than 1 entry)
    H2D_rhoCentrality_rebinned[iDataset]->Scale(1./H2D_rhoCentrality_rebinned[iDataset]->Integral(1, H2D_rhoCentrality_rebinned[iDataset]->GetNbinsX(), 1, H2D_rhoCentrality_rebinned[iDataset]->GetNbinsY(), "width"));
  }

  TString* pdfName = new TString("jet_"+jetType[iJetType]+"_"+jetLevel[iJetLevel]+"_DataComp_rho_vs_centrality");

  // TString textContext(contextDatasetComp(""));
  TString textContext(contextCustomOneField(*texDatasetsComparisonCommonDenominator, ""));

  Draw_TH2_Histograms(H2D_rhoCentrality_rebinned, DatasetsNames, nDatasets, textContext, pdfName, texCentrality, texRho, texCollisionDataInfo, drawnWindowAuto, "logz,autoRangeSame");
}


void Draw_SelectedMultiplicity_vs_Centrality_DatasetComp() {

  TH2D* H2D_multCentrality[nDatasets];
  TH2D* H2D_multCentrality_rebinned[nDatasets];

  for(int iDataset = 0; iDataset < nDatasets; iDataset++){
    H2D_multCentrality[iDataset] = (TH2D*)((TH2D*)file_O2Analysis_list[iDataset]->Get(analysisWorkflow[iDataset]+"/h2_centrality_ntracks"))->Clone("Draw_SelectedMultiplicity_vs_Centrality_DatasetComp"+Datasets[iDataset]);
    H2D_multCentrality[iDataset]->Sumw2();

    H2D_multCentrality_rebinned[iDataset] = (TH2D*)H2D_multCentrality[iDataset]->RebinX(10,"H2D_multCentrality_rebinned"+Datasets[iDataset]);
    H2D_multCentrality_rebinned[iDataset]->GetXaxis()->SetRange(0,H2D_multCentrality_rebinned[iDataset]->GetNbinsX()-10);
    // H2D_multCentrality_rebinned[iDataset]->GetYaxis()->SetRange(1, H2D_multCentrality_rebinned[iDataset]->FindLastBinAbove(1, 2)); //(asks for the last bin on the y axis (axis number 2) to have strictly more than 1 entry)
  }

  TString* pdfName = new TString("jet_"+jetType[iJetType]+"_"+jetLevel[iJetLevel]+"_DataComp_selMult_vs_centrality");

  TString textContext(contextCustomOneField(*texDatasetsComparisonCommonDenominator, ""));

  Draw_TH2_Histograms(H2D_multCentrality_rebinned, DatasetsNames, nDatasets, textContext, pdfName, texCentrality, texSelectedMultiplicity, texCollisionDataInfo, drawnWindowAuto, "logz,autoRangeSame");
}

void Draw_Rho_vs_SelectedMultiplicity_DatasetComp() {

  TH2D* H2D_rhoMult[nDatasets];
  TH2D* H2D_rhoMult_rebinned[nDatasets];

  for(int iDataset = 0; iDataset < nDatasets; iDataset++){
    H2D_rhoMult[iDataset] = (TH2D*)((TH2D*)file_O2Analysis_list[iDataset]->Get(analysisWorkflow[iDataset]+"/h2_ntracks_rho"))->Clone("Draw_Rho_vs_Mult"+Datasets[iDataset]);
    H2D_rhoMult[iDataset]->Sumw2();

    H2D_rhoMult_rebinned[iDataset] = (TH2D*)H2D_rhoMult[iDataset]->RebinX(1.,"H2D_rhoMult_rebinned"+Datasets[iDataset]);
    // H2D_rhoMult_rebinned[iDataset]->GetXaxis()->SetRange(0,H2D_rhoMult_rebinned[iDataset]->FindLastBinAbove(1, 1)); //(asks for the last bin on the x axis (axis number 1) to have strictly more than 1 entry)
    // H2D_rhoMult_rebinned[iDataset]->GetYaxis()->SetRange(1, H2D_rhoMult_rebinned[iDataset]->FindLastBinAbove(1, 2));
    H2D_rhoMult_rebinned[iDataset]->Scale(1./H2D_rhoMult_rebinned[iDataset]->Integral(1, H2D_rhoMult_rebinned[iDataset]->GetNbinsX(), 1, H2D_rhoMult_rebinned[iDataset]->GetNbinsY(), "width"));
  }

  TString* pdfName = new TString("jet_"+jetType[iJetType]+"_"+jetLevel[iJetLevel]+"_DataComp_rho_vs_multiplicitySelected");

  TString textContext(contextCustomOneField(*texDatasetsComparisonCommonDenominator, ""));

  Draw_TH2_Histograms(H2D_rhoMult_rebinned, DatasetsNames, nDatasets, textContext, pdfName, texSelectedMultiplicity, texRho, texCollisionDataInfo, drawnWindowAuto, "logz,autoRangeSame");
}

void Draw_Rho_vs_SelectedMultiplicity_DatasetComp_withCutDemarcation() {

  TH2D* H2D_rhoMult[nDatasets];
  TH2D* H2D_rhoMult_rebinned[nDatasets];

  for(int iDataset = 0; iDataset < nDatasets; iDataset++){
    H2D_rhoMult[iDataset] = (TH2D*)((TH2D*)file_O2Analysis_list[iDataset]->Get(analysisWorkflow[iDataset]+"/h2_ntracks_rho"))->Clone("Draw_Rho_vs_Mult"+Datasets[iDataset]);
    H2D_rhoMult[iDataset]->Sumw2();

    H2D_rhoMult_rebinned[iDataset] = (TH2D*)H2D_rhoMult[iDataset]->RebinX(1.,"H2D_rhoMult_rebinned"+Datasets[iDataset]);
    // H2D_rhoMult_rebinned[iDataset]->GetXaxis()->SetRange(0,H2D_rhoMult_rebinned[iDataset]->FindLastBinAbove(1, 1)); //(asks for the last bin on the x axis (axis number 1) to have strictly more than 1 entry)
    // H2D_rhoMult_rebinned[iDataset]->GetYaxis()->SetRange(1, H2D_rhoMult_rebinned[iDataset]->FindLastBinAbove(1, 2));
  }

  TString* pdfName = new TString("jet_"+jetType[iJetType]+"_"+jetLevel[iJetLevel]+"_DataComp_rho_vs_multiplicitySelected");

  TString textContext(contextCustomOneField(*texDatasetsComparisonCommonDenominator, ""));

  float lineEdgesX[4] = {0, 500, 1200, 1200};
  float lineEdgesY[4] = {10, 35, 84, 100};
  TPolyLine* demarcation = new TPolyLine(4, lineEdgesX, lineEdgesY);
  Draw_TH2_Histograms(H2D_rhoMult_rebinned, DatasetsNames, nDatasets, textContext, pdfName, texSelectedMultiplicity, texRho, texCollisionDataInfo, drawnWindowAuto, "logz,autoRangeSame,drawLines", demarcation);
}


void Draw_Rho_vs_LeadingJetPt_DatasetComp() {

  TH2D* H2D_rhoLeadJet[nDatasets];
  TH2D* H2D_rhoLeadJet_rebinned[nDatasets];

  for(int iDataset = 0; iDataset < nDatasets; iDataset++){
    H2D_rhoLeadJet[iDataset] = (TH2D*)((TH2D*)file_O2Analysis_list[iDataset]->Get(analysisWorkflow[iDataset]+"/h2_leadingjetpt_rho"))->Clone("Draw_Rho_vs_LeadingJetPt"+Datasets[iDataset]);
    H2D_rhoLeadJet[iDataset]->Sumw2();

    H2D_rhoLeadJet_rebinned[iDataset] = (TH2D*)H2D_rhoLeadJet[iDataset]->RebinX(1.,"H2D_rhoLeadJet_rebinned"+Datasets[iDataset]);
    // H2D_rhoLeadJet_rebinned[iDataset]->GetXaxis()->SetRange(0,H2D_rhoLeadJet_rebinned[iDataset]->FindLastBinAbove(1, 1)); //(asks for the last bin on the x axis (axis number 1) to have strictly more than 1 entry)
    // H2D_rhoLeadJet_rebinned[iDataset]->GetYaxis()->SetRange(1, H2D_rhoLeadJet_rebinned[iDataset]->FindLastBinAbove(1, 2));
  }

  TString* pdfName = new TString("jet_"+jetType[iJetType]+"_"+jetLevel[iJetLevel]+"_DataComp_rho_vs_leadingJetPt");

  TString textContext(contextCustomOneField(*texDatasetsComparisonCommonDenominator, ""));

  Draw_TH2_Histograms(H2D_rhoLeadJet_rebinned, DatasetsNames, nDatasets, textContext, pdfName, texLeadJetPt, texRho, texCollisionDataInfo, drawnWindowAuto, "logz,autoRangeSame");
}

void Draw_BkgFluctuations_vs_Centrality_DatasetComp(std::array<std::array<float, 2>, 2> drawnWindow) {

  TH2D* H2D_temp[nDatasets];
  TH2D* H2D_fluctuations[nDatasets];
  TH2D* H2D_fluctuations_rebinned[nDatasets];

  for(int iDataset = 0; iDataset < nDatasets; iDataset++){
    H2D_temp[iDataset] = (TH2D*)file_O2Analysis_list[iDataset]->Get(analysisWorkflow[iDataset]+"/h2_centrality_rhorandomcone");
    if (!H2D_temp[iDataset])
    {
      cout << "h2_centrality_rhorandomcone not found in file; trying h2_centrality_rhoRandomCone" << endl;
      H2D_temp[iDataset] = (TH2D*)file_O2Analysis_list[iDataset]->Get(analysisWorkflow[iDataset]+"/h2_centrality_rhoRandomCone");
      if (!H2D_temp[iDataset]) {
        cout << "h2_centrality_rhoRandomCone not found in file either; aborting Draw_BkgFluctuations_vs_Centrality_DatasetComp()" << endl;
        return;
      }
    }
    H2D_fluctuations[iDataset] = (TH2D*)(H2D_temp[iDataset])->Clone("Draw_BkgFluctuations_vs_Centrality"+Datasets[iDataset]);
    H2D_fluctuations[iDataset]->Sumw2();

    float FluctuLow = -60;
    float FluctuHigh = 60;
    int ibinFluctu_low = H2D_fluctuations[iDataset]->GetYaxis()->FindBin(FluctuLow);
    int ibinFluctu_high = H2D_fluctuations[iDataset]->GetYaxis()->FindBin(FluctuHigh);

    H2D_fluctuations_rebinned[iDataset] = (TH2D*)H2D_fluctuations[iDataset]->RebinX(10,"H2D_rhoCentrality_rebinned"+Datasets[iDataset]);
    H2D_fluctuations_rebinned[iDataset]->GetXaxis()->SetRange(0,H2D_fluctuations_rebinned[iDataset]->GetNbinsX()-10);
    // int symBinLimit = min(H2D_fluctuations_rebinned[iDataset]->FindFirstBinAbove(1, 2), abs(H2D_fluctuations_rebinned[iDataset]->GetNbinsY() - H2D_fluctuations_rebinned[iDataset]->FindLastBinAbove(1, 2))); //(asks for the first/last bin on the y axis (axis number 2) to have strictly more than 1 entry)
    // H2D_fluctuations_rebinned[iDataset]->GetYaxis()->SetRange(symBinLimit, H2D_fluctuations_rebinned[iDataset]->GetNbinsY() - symBinLimit); //getting symmetric window around 0 on Y axis
    H2D_fluctuations_rebinned[iDataset]->Scale(1./H2D_fluctuations_rebinned[iDataset]->Integral(1, H2D_fluctuations_rebinned[iDataset]->GetNbinsX(), 1, H2D_fluctuations_rebinned[iDataset]->GetNbinsY(), "width"));
  }

  TString* pdfName = new TString("jet_"+jetType[iJetType]+"_"+jetLevel[iJetLevel]+"_BkgFluctuations_vs_centrality");
  TString placeHolderLegend[1] = {""};

  TString textContext(contextCustomOneField(*texDatasetsComparisonCommonDenominator, ""));
  Draw_TH2_Histograms(H2D_fluctuations_rebinned, DatasetsNames, nDatasets, textContext, pdfName, texCentrality, texBkgFluctuationRandomCone, texCollisionDataInfo, drawnWindow, "logz,autoRangeSameSym");
}




void Draw_Pt_CentralityComparison(float jetRadius, int iDataset) {

  TH3D* H3D_jetRjetPtjetCent;
  TH1D* H1D_jetPt[nCentralityBins];
  TH1D* H1D_jetPt_rebinned[nCentralityBins];
  
  H3D_jetRjetPtjetCent = (TH3D*)((TH3D*)file_O2Analysis_list[iDataset]->Get(analysisWorkflow[iDataset]+"/h3_jet_r_jet_pt_centrality"+jetFinderQaHistType[iJetFinderQaType]))->Clone("Draw_Pt_CentralityComparison"+Datasets[iDataset]+Form("%.1f", jetRadius));
  H3D_jetRjetPtjetCent->Sumw2();

  int ibinJetRadius = H3D_jetRjetPtjetCent->GetXaxis()->FindBin(jetRadius+GLOBAL_epsilon);
  int ibinCent_low, ibinCent_high;
  TString CentralityLegend[nCentralityBins];
  std::stringstream ss;
  for(int iCentralityBin = 0; iCentralityBin < nCentralityBins; iCentralityBin++){

    ibinCent_low = H3D_jetRjetPtjetCent->GetZaxis()->FindBin(arrayCentralityBinning[iCentralityBin]);
    ibinCent_high = H3D_jetRjetPtjetCent->GetZaxis()->FindBin(arrayCentralityBinning[iCentralityBin+1])-1;
    H1D_jetPt[iCentralityBin] = (TH1D*)H3D_jetRjetPtjetCent->ProjectionY("jetPt_"+Datasets[iDataset]+Form("%.1f", jetRadius)+"_@cent["+Form("%.1d", ibinCent_low)+","+Form("%.1d", ibinCent_high)+"]", ibinJetRadius, ibinJetRadius, ibinCent_low, ibinCent_high, "e");
    H1D_jetPt_rebinned[iCentralityBin] = (TH1D*)H1D_jetPt[iCentralityBin]->Rebin(1.,"jetPt_rebinned_"+Datasets[iDataset]+Form("%.1f", jetRadius)+"_@cent["+Form("%.1d", ibinCent_low)+","+Form("%.1d", ibinCent_high)+"]");

    // NormaliseYieldToNEntries(H1D_jetPt_rebinned[iRadius]);
    // NormaliseYieldToNEvents(H1D_jetPt_rebinned[iCentralityBin], GetNEventsSelected_JetFramework(file_O2Analysis_list[iDataset]));
    NormaliseYieldToNEvents(H1D_jetPt_rebinned[iCentralityBin], GetNEventsSelectedCentrality_JetFramework(file_O2Analysis_list[iDataset], arrayCentralityBinning[iCentralityBin], arrayCentralityBinning[iCentralityBin+1], trainId));
    // NormaliseYieldToNEvents(H1D_jetPt_rebinned[iCentralityBin], GetNEventsSelected_JetFramework(file_O2Analysis_list[iDataset]));

    ss << "Cent " << arrayCentralityBinning[iCentralityBin] << " - " << arrayCentralityBinning[iCentralityBin+1] << " ";
    CentralityLegend[iCentralityBin] = (TString)ss.str();
    ss.str("");
    ss.clear();
  }

  TString* pdfName = new TString("jet_"+jetType[iJetType]+"_"+jetLevel[iJetLevel]+"_CentralityComp_R="+Form("%.1f", jetRadius)+"_"+DatasetsNames[iDataset]+"_Pt"+jetFinderQaHistType[iJetFinderQaType]);

  // TString textContext(contextDatasetCompAndRadius(jetRadius, ""));
  TString textContext(contextCustomTwoFields(*texDatasetsComparisonCommonDenominator, contextJetRadius(jetRadius), ""));

  std::array<std::array<float, 2>, 2> drawnWindow = {{{0, 200}, {-999, -999}}}; // {{xmin, xmax}, {ymin, ymax}}
  if(iJetFinderQaType == 1) {
    drawnWindow = {{{-20, 200}, {-999, -999}}};
  }
  Draw_TH1_Histograms_in_one(H1D_jetPt_rebinned, CentralityLegend, nCentralityBins, textContext, pdfName, texPtX, texJetPtYield_EventNorm_CentWindow, texCollisionDataInfo, drawnWindow, "logy");
  // Draw_TH1_Histograms_in_one(H1D_jetPt_rebinned, CentralityLegend, nCentralityBins, textContext, pdfName, texPtX, texJetPtYield_EventNorm, texCollisionDataInfo, drawnWindowAuto, "logy");
  TString* pdfName2 = new TString("jet_"+jetType[iJetType]+"_"+jetLevel[iJetLevel]+"_CentralityComp_R="+Form("%.1f", jetRadius)+"_"+DatasetsNames[iDataset]+"_Pt"+jetFinderQaHistType[iJetFinderQaType]+"_logx");
  const std::array<std::array<float, 2>, 2> drawnWindowLog = {{{0.4, 200}, {-999, -999}}}; // {{xmin, xmax}, {ymin, ymax}}

  Draw_TH1_Histograms_in_one(H1D_jetPt_rebinned, CentralityLegend, nCentralityBins, textContext, pdfName2, texPtX, texJetPtYield_EventNorm_CentWindow, texCollisionDataInfo, drawnWindowLog, "logx,logy");

}



void Draw_Eta_CentralityComparison(float jetRadius, int iDataset) { //for now only works with one Radius, the one defined in the jetfinderQA as selectedJetsRadius

  TH2D* H2D_jetCentjetEta;
  TH1D* H1D_jetEta[nCentralityBins];
  TH1D* H1D_jetEta_rebinned[nCentralityBins];
  
  H2D_jetCentjetEta = (TH2D*)((TH2D*)file_O2Analysis_list[iDataset]->Get(analysisWorkflow[iDataset]+"/h2_centrality_jet_eta"+jetFinderQaHistType[iJetFinderQaType]))->Clone("Draw_Eta_CentralityComparison"+Datasets[iDataset]+Form("%.1f", jetRadius));
  H2D_jetCentjetEta->Sumw2();

  // int ibinJetRadius = H3D_jetRjetEtajetCent->GetXaxis()->FindBin(jetRadius+GLOBAL_epsilon);
  int ibinCent_low, ibinCent_high;
  TString CentralityLegend[nCentralityBins];
  std::stringstream ss;
  for(int iCentralityBin = 0; iCentralityBin < nCentralityBins; iCentralityBin++){

    ibinCent_low = H2D_jetCentjetEta->GetXaxis()->FindBin(arrayCentralityBinning[iCentralityBin]);
    ibinCent_high = H2D_jetCentjetEta->GetXaxis()->FindBin(arrayCentralityBinning[iCentralityBin+1])-1;
    // H1D_jetEta[iCentralityBin] = (TH1D*)H3D_jetRjetEtajetCent->ProjectionY("jetEta_"+Datasets[iDataset]+Form("%.1f", jetRadius)+"_@cent["+Form("%.1d", ibinCent_low)+","+Form("%.1d", ibinCent_high)+"]", ibinJetRadius, ibinJetRadius, ibinCent_low, ibinCent_high, "e");
    H1D_jetEta[iCentralityBin] = (TH1D*)H2D_jetCentjetEta->ProjectionY("jetEta_"+Datasets[iDataset]+Form("%.1f", jetRadius)+"_@cent["+Form("%.1d", ibinCent_low)+","+Form("%.1d", ibinCent_high)+"]", ibinCent_low, ibinCent_high, "e");
    H1D_jetEta_rebinned[iCentralityBin] = (TH1D*)H1D_jetEta[iCentralityBin]->Rebin(1.,"jetEta_rebinned_"+Datasets[iDataset]+Form("%.1f", jetRadius)+"_@cent["+Form("%.1d", ibinCent_low)+","+Form("%.1d", ibinCent_high)+"]");

    // NormaliseYieldToNEntries(H1D_jetEta_rebinned[iRadius]);
    NormaliseYieldToNEvents(H1D_jetEta_rebinned[iCentralityBin], GetNEventsSelectedCentrality_JetFramework(file_O2Analysis_list[iDataset], arrayCentralityBinning[iCentralityBin], arrayCentralityBinning[iCentralityBin+1], trainId));

    ss << "Cent " << arrayCentralityBinning[iCentralityBin] << " - " << arrayCentralityBinning[iCentralityBin+1] << " ";
    CentralityLegend[iCentralityBin] = (TString)ss.str();
    ss.str("");
    ss.clear();
  }

  TString* pdfName = new TString("jet_"+jetType[iJetType]+"_"+jetLevel[iJetLevel]+"_CentralityComp_R="+Form("%.1f", jetRadius)+"_"+DatasetsNames[iDataset]+"_Eta"+jetFinderQaHistType[iJetFinderQaType]);

  // TString textContext(contextDatasetCompAndRadius(jetRadius, ""));
  TString textContext(contextCustomTwoFields(*texDatasetsComparisonCommonDenominator, contextJetRadius(jetRadius), ""));

  const std::array<std::array<float, 2>, 2> drawnWindow = {{{-1, 1}, {-999, -999}}}; // {{xmin, xmax}, {ymin, ymax}}

  Draw_TH1_Histograms_in_one(H1D_jetEta_rebinned, CentralityLegend, nCentralityBins, textContext, pdfName, texEtaX, texJetEtaYield_EventNorm_CentWindow, texCollisionDataInfo, drawnWindow, "");
  // TString* pdfName2 = new TString("jet_"+jetType[iJetType]+"_"+jetLevel[iJetLevel]+"_CentralityComp_R="+Form("%.1f", jetRadius)+"_"+Datasets[iDataset]+"_Eta"+jetFinderQaHistType[iJetFinderQaType]+"_logx");
  // Draw_TH1_Histograms_in_one(H1D_jetEta_rebinned, CentralityLegend, nCentralityBins, textContext, pdfName2, texEtaX, texJetEtaYield_EventNorm_CentWindow, texCollisionDataInfo, drawnWindowAuto, "logx,logy");
}

void Draw_Phi_CentralityComparison(float jetRadius, int iDataset) { //for now only works with one Radius, the one defined in the jetfinderQA as selectedJetsRadius

  TH2D* H2D_jetCentjetPhi;
  TH1D* H1D_jetPhi[nCentralityBins];
  TH1D* H1D_jetPhi_rebinned[nCentralityBins];
  
  H2D_jetCentjetPhi = (TH2D*)((TH2D*)file_O2Analysis_list[iDataset]->Get(analysisWorkflow[iDataset]+"/h2_centrality_jet_phi"+jetFinderQaHistType[iJetFinderQaType]))->Clone("Draw_Phi_CentralityComparison"+Datasets[iDataset]+Form("%.1f", jetRadius));
  H2D_jetCentjetPhi->Sumw2();

  // int ibinJetRadius = H3D_jetRjetPhijetCent->GetXaxis()->FindBin(jetRadius+GLOBAL_epsilon);
  int ibinCent_low, ibinCent_high;
  TString CentralityLegend[nCentralityBins];
  std::stringstream ss;
  for(int iCentralityBin = 0; iCentralityBin < nCentralityBins; iCentralityBin++){

    ibinCent_low = H2D_jetCentjetPhi->GetXaxis()->FindBin(arrayCentralityBinning[iCentralityBin]);
    ibinCent_high = H2D_jetCentjetPhi->GetXaxis()->FindBin(arrayCentralityBinning[iCentralityBin+1])-1;
    // H1D_jetPhi[iCentralityBin] = (TH1D*)H3D_jetRjetPhijetCent->ProjectionY("jetPhi_"+Datasets[iDataset]+Form("%.1f", jetRadius)+"_@cent["+Form("%.1d", ibinCent_low)+","+Form("%.1d", ibinCent_high)+"]", ibinJetRadius, ibinJetRadius, ibinCent_low, ibinCent_high, "e");
    H1D_jetPhi[iCentralityBin] = (TH1D*)H2D_jetCentjetPhi->ProjectionY("jetPhi_"+Datasets[iDataset]+Form("%.1f", jetRadius)+"_@cent["+Form("%.1d", ibinCent_low)+","+Form("%.1d", ibinCent_high)+"]", ibinCent_low, ibinCent_high, "e");
    H1D_jetPhi_rebinned[iCentralityBin] = (TH1D*)H1D_jetPhi[iCentralityBin]->Rebin(1.,"jetPhi_rebinned_"+Datasets[iDataset]+Form("%.1f", jetRadius)+"_@cent["+Form("%.1d", ibinCent_low)+","+Form("%.1d", ibinCent_high)+"]");

    // NormaliseYieldToNEntries(H1D_jetPhi_rebinned[iRadius]);
    NormaliseYieldToNEvents(H1D_jetPhi_rebinned[iCentralityBin], GetNEventsSelectedCentrality_JetFramework(file_O2Analysis_list[iDataset], arrayCentralityBinning[iCentralityBin], arrayCentralityBinning[iCentralityBin+1], trainId));


    ss << "Cent " << arrayCentralityBinning[iCentralityBin] << " - " << arrayCentralityBinning[iCentralityBin+1] << " ";
    CentralityLegend[iCentralityBin] = (TString)ss.str();
    ss.str("");
    ss.clear();
  }

  TString* pdfName = new TString("jet_"+jetType[iJetType]+"_"+jetLevel[iJetLevel]+"_CentralityComp_R="+Form("%.1f", jetRadius)+"_"+DatasetsNames[iDataset]+"_Phi"+jetFinderQaHistType[iJetFinderQaType]);

  // TString textContext(contextDatasetCompAndRadius(jetRadius, ""));
  TString textContext(contextCustomTwoFields(*texDatasetsComparisonCommonDenominator, contextJetRadius(jetRadius), ""));

  Draw_TH1_Histograms_in_one(H1D_jetPhi_rebinned, CentralityLegend, nCentralityBins, textContext, pdfName, texPhiX, texJetPhiYield_EventNorm_CentWindow, texCollisionDataInfo, drawnWindowAuto, "histWithLine");
  // TString* pdfName2 = new TString("jet_"+jetType[iJetType]+"_"+jetLevel[iJetLevel]+"_CentralityComp_R="+Form("%.1f", jetRadius)+"_"+Datasets[iDataset]+"_Phi"+jetFinderQaHistType[iJetFinderQaType]+"_logx");
  // Draw_TH1_Histograms_in_one(H1D_jetPhi_rebinned, CentralityLegend, nCentralityBins, textContext, pdfName2, texPhiX, texJetPhiYield_EventNorm_CentWindow, texCollisionDataInfo, drawnWindowAuto, "logx,logy");
}

void Draw_BkgFluctuations_CentralityProjection(int iDataset, std::array<std::array<float, 2>, 2> drawnWindowZoom, const char options[]) {

  TH2D* H2D_fluctuations_centrality;

  TH1D* H1D_fluctuations[nCentralityBins];
  TH1D* H1D_fluctuations_rebinned[nCentralityBins];

  H2D_fluctuations_centrality = (TH2D*)((TH2D*)file_O2Analysis_list[iDataset]->Get(analysisWorkflow[iDataset]+"/h2_centrality_rhorandomcone"))->Clone("Draw_BkgFluctuations_CentralityProjection"+Datasets[iDataset]);
  H2D_fluctuations_centrality->Sumw2();

  float FluctuLow = -60;
  float FluctuHigh = 60;
  int ibinFluctu_low  = H2D_fluctuations_centrality->GetYaxis()->FindBin(FluctuLow);
  int ibinFluctu_high = H2D_fluctuations_centrality->GetYaxis()->FindBin(FluctuHigh);
  H2D_fluctuations_centrality->GetYaxis()->SetRange(ibinFluctu_low, ibinFluctu_high);

  int ibinCent_low, ibinCent_high;
  TString CentralityLegend[nCentralityBins];
  std::stringstream ss;

  TString* yAxisLabel;

  for(int iCentralityBin = 0; iCentralityBin < nCentralityBins; iCentralityBin++){
    ibinCent_low = H2D_fluctuations_centrality->GetXaxis()->FindBin(arrayCentralityBinning[iCentralityBin]);
    ibinCent_high = H2D_fluctuations_centrality->GetXaxis()->FindBin(arrayCentralityBinning[iCentralityBin+1])-1;
    // cout << "ibinCent_low = " << ibinCent_low << ", ibinCent_high = " << ibinCent_high << endl;
    H1D_fluctuations[iCentralityBin] = (TH1D*)H2D_fluctuations_centrality->ProjectionY("bkgFluctuationCentrality_"+Datasets[iDataset]+"_@cent["+Form("%.1d", ibinCent_low)+","+Form("%.1d", ibinCent_high)+"]", ibinCent_low, ibinCent_high, "e");

    H1D_fluctuations_rebinned[iCentralityBin] = (TH1D*)H1D_fluctuations[iCentralityBin]->Rebin(1.,"bkgFluctuationCentrality_rebinned_"+Datasets[iDataset]+"_@cent["+Form("%.1d", ibinCent_low)+","+Form("%.1d", ibinCent_high)+"]");

    if (strstr(options, "normEntries") != NULL) {
      NormaliseYieldToNEntries(H1D_fluctuations_rebinned[iCentralityBin]);
      yAxisLabel = texEntriesNorm_BkgFluctuationYield;
    }
    if (strstr(options, "normEvents") != NULL) {
      NormaliseYieldToNEvents(H1D_fluctuations_rebinned[iCentralityBin], GetNEventsSelected_JetFramework(file_O2Analysis_list[iDataset]));
      yAxisLabel = texCollNorm_BkgFluctuationYield;
    }
    if (strstr(options, "normEventsCentrality") != NULL) {
      NormaliseYieldToNEvents(H1D_fluctuations_rebinned[iCentralityBin], GetNEventsSelectedCentrality_JetFramework(file_O2Analysis_list[iDataset], arrayCentralityBinning[iCentralityBin], arrayCentralityBinning[iCentralityBin+1], trainId));
      yAxisLabel = texCollNorm_BkgFluctuationYield_CentWindow;
    }
    
    ss << "Cent " << arrayCentralityBinning[iCentralityBin] << " - " << arrayCentralityBinning[iCentralityBin+1] << " ";
    CentralityLegend[iCentralityBin] = (TString)ss.str();
    ss.str("");
    ss.clear();
  }

  TString* pdfName = new TString("jet_"+jetType[iJetType]+"_"+jetLevel[iJetLevel]+"_CentralityComp_"+DatasetsNames[iDataset]+"_BkgFluctuations_logy");

  TString textContext(contextCustomOneField(*texDatasetsComparisonCommonDenominator, ""));

  Draw_TH1_Histograms_in_one(H1D_fluctuations_rebinned, CentralityLegend, nCentralityBins, textContext, pdfName, texBkgFluctuationRandomCone, yAxisLabel, texCollisionDataInfo, drawnWindowAuto, "logy");
  TString* pdfName2 = new TString("jet_"+jetType[iJetType]+"_"+jetLevel[iJetLevel]+"_CentralityComp_"+DatasetsNames[iDataset]+"_BkgFluctuations");
  Draw_TH1_Histograms_in_one(H1D_fluctuations_rebinned, CentralityLegend, nCentralityBins, textContext, pdfName2, texBkgFluctuationRandomCone, yAxisLabel, texCollisionDataInfo, drawnWindowAuto, "");

  // zoom around 0
  // for(int iCentralityBin = 0; iCentralityBin < nCentralityBins; iCentralityBin++){
  //   H1D_fluctuations_rebinned[iCentralityBin]->GetXaxis()->SetRangeUser(zoomX[0], zoomX[1]);
  // }
  TString* pdfName_zoom = new TString("jet_"+jetType[iJetType]+"_"+jetLevel[iJetLevel]+"_CentralityComp_"+DatasetsNames[iDataset]+"_BkgFluctuations_logy_zoom");
  Draw_TH1_Histograms_in_one(H1D_fluctuations_rebinned, CentralityLegend, nCentralityBins, textContext, pdfName_zoom, texBkgFluctuationRandomCone, yAxisLabel, texCollisionDataInfo, drawnWindowZoom, "logy");
}

void Draw_Rho_CentralityProjection(int iDataset, const char options[]) {

  TH2D* H2D_rhoCentrality;

  TH1D* H1D_rho[nCentralityBins];
  TH1D* H1D_rho_rebinned[nCentralityBins];

  H2D_rhoCentrality = (TH2D*)((TH2D*)file_O2Analysis_list[iDataset]->Get(analysisWorkflow[iDataset]+"/h2_centrality_rho"))->Clone("Draw_Rho_CentralityProjection"+Datasets[iDataset]);
  H2D_rhoCentrality->Sumw2();

  // float FluctuLow = -60;
  // float FluctuHigh = 60;
  // int ibinFluctu_low  = H2D_rhoCentrality->GetYaxis()->FindBin(FluctuLow);
  // int ibinFluctu_high = H2D_rhoCentrality->GetYaxis()->FindBin(FluctuHigh);
  // H2D_rhoCentrality->GetYaxis()->SetRange(ibinFluctu_low, ibinFluctu_high);

  int ibinCent_low, ibinCent_high;
  TString CentralityLegend[nCentralityBins];
  std::stringstream ss;

  TString* yAxisLabel;

  for(int iCentralityBin = 0; iCentralityBin < nCentralityBins; iCentralityBin++){
    ibinCent_low = H2D_rhoCentrality->GetXaxis()->FindBin(arrayCentralityBinning[iCentralityBin]);
    ibinCent_high = H2D_rhoCentrality->GetXaxis()->FindBin(arrayCentralityBinning[iCentralityBin+1])-1;
    // cout << "ibinCent_low = " << ibinCent_low << ", ibinCent_high = " << ibinCent_high << endl;
    H1D_rho[iCentralityBin] = (TH1D*)H2D_rhoCentrality->ProjectionY("bkgRhoCentrality_"+Datasets[iDataset]+"_@cent["+Form("%.1d", ibinCent_low)+","+Form("%.1d", ibinCent_high)+"]", ibinCent_low, ibinCent_high, "e");

    H1D_rho_rebinned[iCentralityBin] = (TH1D*)H1D_rho[iCentralityBin]->Rebin(1.,"bkgRhoCentrality_rebinned_"+Datasets[iDataset]+"_@cent["+Form("%.1d", ibinCent_low)+","+Form("%.1d", ibinCent_high)+"]");

    if (strstr(options, "normEntries") != NULL) {
      NormaliseYieldToNEntries(H1D_rho_rebinned[iCentralityBin]);
      yAxisLabel = texEntriesNormRho;
    }
    if (strstr(options, "normEvents") != NULL) {
      NormaliseYieldToNEvents(H1D_rho_rebinned[iCentralityBin], GetNEventsSelected_JetFramework(file_O2Analysis_list[iDataset]));
      yAxisLabel = texCollNorm_RhoYield;
    }
    if (strstr(options, "normEventsCentrality") != NULL) {
      NormaliseYieldToNEvents(H1D_rho_rebinned[iCentralityBin], GetNEventsSelectedCentrality_JetFramework(file_O2Analysis_list[iDataset], arrayCentralityBinning[iCentralityBin], arrayCentralityBinning[iCentralityBin+1], trainId));
      yAxisLabel = texCollNorm_RhoYield_CentWindow;
    }
    
    ss << "Cent " << arrayCentralityBinning[iCentralityBin] << " - " << arrayCentralityBinning[iCentralityBin+1] << " ";
    CentralityLegend[iCentralityBin] = (TString)ss.str();
    ss.str("");
    ss.clear();
  }

  TString* pdfName = new TString("jet_"+jetType[iJetType]+"_"+jetLevel[iJetLevel]+"_CentralityComp_"+DatasetsNames[iDataset]+"_BkgRho_logy");

  TString textContext(contextCustomOneField(*texDatasetsComparisonCommonDenominator, ""));

  Draw_TH1_Histograms_in_one(H1D_rho_rebinned, CentralityLegend, nCentralityBins, textContext, pdfName, texRho, yAxisLabel, texCollisionDataInfo, drawnWindowAuto, "autoXrange,logy");
  TString* pdfName2 = new TString("jet_"+jetType[iJetType]+"_"+jetLevel[iJetLevel]+"_CentralityComp_"+DatasetsNames[iDataset]+"_BkgRho");
  Draw_TH1_Histograms_in_one(H1D_rho_rebinned, CentralityLegend, nCentralityBins, textContext, pdfName2, texRho, yAxisLabel, texCollisionDataInfo, drawnWindowAuto, "autoXrange");
}

void Draw_Rho_CentralityProjection_DatasetComp(float* centRange, const char options[]) {

  TH2D* H2D_rhoCentrality[nDatasets];

  TH1D* H1D_rho[nDatasets];
  TH1D* H1D_rho_rebinned[nDatasets];


  // float FluctuLow = -60;
  // float FluctuHigh = 60;
  // int ibinFluctu_low  = H2D_rhoCentrality->GetYaxis()->FindBin(FluctuLow);
  // int ibinFluctu_high = H2D_rhoCentrality->GetYaxis()->FindBin(FluctuHigh);
  // H2D_rhoCentrality->GetYaxis()->SetRange(ibinFluctu_low, ibinFluctu_high);

  int ibinCent_low, ibinCent_high;

  TString CentralityLegend;
  std::stringstream ss;
  TString* yAxisLabel;

  for(int iDataset = 0; iDataset < nDatasets; iDataset++){
    H2D_rhoCentrality[iDataset] = (TH2D*)((TH2D*)file_O2Analysis_list[iDataset]->Get(analysisWorkflow[iDataset]+"/h2_centrality_rho"))->Clone("Draw_Rho_CentralityProjection"+Datasets[iDataset]);
    H2D_rhoCentrality[iDataset]->Sumw2();
    ibinCent_low = H2D_rhoCentrality[iDataset]->GetXaxis()->FindBin(centRange[0]);
    ibinCent_high = H2D_rhoCentrality[iDataset]->GetXaxis()->FindBin(centRange[1])-1;

    H1D_rho[iDataset] = (TH1D*)H2D_rhoCentrality[iDataset]->ProjectionY("bkgRhoCentrality_"+Datasets[iDataset]+"_@cent["+Form("%.1d", ibinCent_low)+","+Form("%.1d", ibinCent_high)+"]", ibinCent_low, ibinCent_high, "e");

    H1D_rho_rebinned[iDataset] = (TH1D*)H1D_rho[iDataset]->Rebin(1.,"bkgRhoCentrality_rebinned_"+Datasets[iDataset]+"_@cent["+Form("%.1d", ibinCent_low)+","+Form("%.1d", ibinCent_high)+"]");

    if (strstr(options, "normEntries") != NULL) {
      NormaliseYieldToNEntries(H1D_rho_rebinned[iDataset]);
      yAxisLabel = texEntriesNormRho;
    }
    if (strstr(options, "normEvents") != NULL) {
      NormaliseYieldToNEvents(H1D_rho_rebinned[iDataset], GetNEventsSelected_JetFramework(file_O2Analysis_list[iDataset]));
      yAxisLabel = texCollNorm_RhoYield;
    }
    if (strstr(options, "normEventsCentrality") != NULL) {
      NormaliseYieldToNEvents(H1D_rho_rebinned[iDataset], GetNEventsSelectedCentrality_JetFramework(file_O2Analysis_list[iDataset], centRange[0], centRange[1], trainId));
      yAxisLabel = texCollNorm_RhoYield_CentWindow;
    }
    
    ss << "Cent " << centRange[0] << "% - " << centRange[1] << "% ";
    CentralityLegend = (TString)ss.str();
    ss.str("");
    ss.clear();
  }

  TString* pdfName = new TString("jet_"+jetType[iJetType]+"_"+jetLevel[iJetLevel]+"_DatasetComp_BkgRho"+"_@cent["+Form("%.1d", ibinCent_low)+","+Form("%.1d", ibinCent_high)+"]_logy");

  TString textContext(contextCustomTwoFields(*texDatasetsComparisonCommonDenominator, CentralityLegend, ""));

  Draw_TH1_Histograms_in_one(H1D_rho_rebinned, DatasetsNames, nDatasets, textContext, pdfName, texRho, yAxisLabel, texCollisionDataInfo, drawnWindowAuto, "autoXrange,logy");
  TString* pdfName2 = new TString("jet_"+jetType[iJetType]+"_"+jetLevel[iJetLevel]+"_DatasetComp_BkgRho"+"_@cent["+Form("%.1d", ibinCent_low)+","+Form("%.1d", ibinCent_high)+"]");
  Draw_TH1_Histograms_in_one(H1D_rho_rebinned, DatasetsNames, nDatasets, textContext, pdfName2, texRho, yAxisLabel, texCollisionDataInfo, drawnWindowAuto, "autoXrange");
}

void Draw_RhoMean_asFunctionOf_Centrality(int iDataset, const char options[]) {

  TH2D* H2D_rhoCentrality;

  TH1D* H1D_rho[nCentralityBins];
  TH1D* H1D_rho_rebinned[nCentralityBins];

  H2D_rhoCentrality = (TH2D*)((TH2D*)file_O2Analysis_list[iDataset]->Get(analysisWorkflow[iDataset]+"/h2_centrality_rho"))->Clone("Draw_Rho_CentralityProjection"+Datasets[iDataset]);
  H2D_rhoCentrality->Sumw2();

  // float FluctuLow = -60;
  // float FluctuHigh = 60;
  // int ibinFluctu_low  = H2D_rhoCentrality->GetYaxis()->FindBin(FluctuLow);
  // int ibinFluctu_high = H2D_rhoCentrality->GetYaxis()->FindBin(FluctuHigh);
  // H2D_rhoCentrality->GetYaxis()->SetRange(ibinFluctu_low, ibinFluctu_high);

  double RhoMean[nCentralityBins], RhoMeanError[nCentralityBins];

  int ibinCent_low, ibinCent_high;
  for(int iCentralityBin = 0; iCentralityBin < nCentralityBins; iCentralityBin++){
    ibinCent_low = H2D_rhoCentrality->GetXaxis()->FindBin(arrayCentralityBinning[iCentralityBin]);
    ibinCent_high = H2D_rhoCentrality->GetXaxis()->FindBin(arrayCentralityBinning[iCentralityBin+1])-1;
    // cout << "ibinCent_low = " << ibinCent_low << ", ibinCent_high = " << ibinCent_high << endl;
    H1D_rho[iCentralityBin] = (TH1D*)H2D_rhoCentrality->ProjectionY("bkgRhoCentrality_"+Datasets[iDataset]+"_@cent["+Form("%.1d", ibinCent_low)+","+Form("%.1d", ibinCent_high)+"]", ibinCent_low, ibinCent_high, "e");

    // H1D_rho_rebinned[iCentralityBin] = (TH1D*)H1D_rho[iCentralityBin]->Rebin(1.,"bkgRhoCentrality_rebinned_"+Datasets[iDataset]+"_@cent["+Form("%.1d", ibinCent_low)+","+Form("%.1d", ibinCent_high)+"]");

    // H1D_rho_rebinned[iCentralityBin]

    RhoMean[iCentralityBin] = H1D_rho[iCentralityBin]->GetMean();
    RhoMeanError[iCentralityBin] = H1D_rho[iCentralityBin]->GetMeanError();

    // if (strstr(options, "normEntries") != NULL) {
    //   NormaliseYieldToNEntries(H1D_rho_rebinned[iCentralityBin]);
    // }
    // if (strstr(options, "normEvents") != NULL) {
    //   NormaliseYieldToNEvents(H1D_rho_rebinned[iCentralityBin], GetNEventsSelected_JetFramework(file_O2Analysis_list[iDataset]));
    // }
    // if (strstr(options, "normEventsCentrality") != NULL) {
    //   NormaliseYieldToNEvents(H1D_rho_rebinned[iCentralityBin], GetNEventsSelectedCentrality_JetFramework(file_O2Analysis_list[iDataset], arrayCentralityBinning[iCentralityBin], arrayCentralityBinning[iCentralityBin+1], trainId));
    // }
  }

  TH1D* H1D_MeanRho_asFunctionOf_Centrality = new TH1D("H1D_MeanRho_asFunctionOf_Centrality", "H1D_MeanRho_asFunctionOf_Centrality", nCentralityBins, arrayCentralityBinning);
  H1D_MeanRho_asFunctionOf_Centrality->Sumw2();

  for(int iCentralityBin = 0; iCentralityBin < nCentralityBins; iCentralityBin++){
    H1D_MeanRho_asFunctionOf_Centrality->SetBinContent(iCentralityBin+1, RhoMean[iCentralityBin]);
    H1D_MeanRho_asFunctionOf_Centrality->SetBinError(iCentralityBin+1, RhoMeanError[iCentralityBin]);
  }

  TString textContext(contextCustomOneField(*texDatasetsComparisonCommonDenominator, ""));

  TString* pdfName = new TString("jet_"+jetType[iJetType]+"_"+jetLevel[iJetLevel]+"_CentralityComp_"+DatasetsNames[iDataset]+"_RhoMean_asFunctionOf_Centrality");
  Draw_TH1_Histogram(H1D_MeanRho_asFunctionOf_Centrality, textContext, pdfName, texCentrality, texRhoMean, texCollisionDataInfo, drawnWindowAuto, "");
}

void Draw_Eta_PtCutComparison(float jetRadius, int iDataset, float* PtCuts, int nPtCut, const char options[]) {

  TH3D* H3D_jetRjetPtjetEta;
  TH1D* H1D_jetEta[nPtCut];
  TH1D* H1D_jetEta_rebinned[nPtCut];
  
  H3D_jetRjetPtjetEta = (TH3D*)((TH3D*)file_O2Analysis_list[iDataset]->Get(analysisWorkflow[iDataset]+"/h3_jet_r_jet_pt_jet_eta"+jetFinderQaHistType[iJetFinderQaType]))->Clone("Draw_Eta_PtCutComparison"+Datasets[iDataset]);
  H3D_jetRjetPtjetEta->Sumw2();

  TString PtCutsLegend[nPtCut];
  std::stringstream ss;

  int ibinJetRadius = H3D_jetRjetPtjetEta->GetXaxis()->FindBin(jetRadius+GLOBAL_epsilon);

  TString* yAxisLabel;

  for(int iBinPt = 0; iBinPt < nPtCut; iBinPt++){

    float PtCutLow = PtCuts[iBinPt];
    int ibinPt_low = H3D_jetRjetPtjetEta->GetYaxis()->FindBin(PtCutLow);
    if (ibinPt_low == 0) 
      cout << "WARNING: Eta_RadiusComparison is counting the underflow with the chosen PtRange" << endl;

    H1D_jetEta[iBinPt] = (TH1D*)H3D_jetRjetPtjetEta->ProjectionZ("jetEta_"+DatasetsNames[iDataset]+Form("%.1f", jetRadius)+Form("%.1f", PtCutLow)+"<pt", ibinJetRadius, ibinJetRadius, ibinPt_low, H3D_jetRjetPtjetEta->GetYaxis()->GetNbins(), "e");
    H1D_jetEta_rebinned[iBinPt] = (TH1D*)H1D_jetEta[iBinPt]->Rebin(2.,"jetEta_rebinned_"+DatasetsNames[iDataset]+Form("%.1f", jetRadius)+Form("%.1f", PtCutLow)+"<pt");

    if (strstr(options, "normEntries") != NULL) {
      NormaliseYieldToNEntries(H1D_jetEta_rebinned[iBinPt]);
      yAxisLabel = texJetEtaYield_EntriesNorm;
    }
    if (strstr(options, "normEvents") != NULL) {
      NormaliseYieldToNEvents(H1D_jetEta_rebinned[iBinPt], GetNEventsSelected_JetFramework(file_O2Analysis_list[iDataset]));
      yAxisLabel = texJetEtaYield_EventNorm;
    }


    ss << "Pt > " << PtCuts[iBinPt];
    PtCutsLegend[iBinPt] = (TString)ss.str();
    ss.str("");
    ss.clear();
  }

  TString* pdfName = new TString("jet_"+jetType[iJetType]+"_"+jetLevel[iJetLevel]+"_"+DatasetsNames[iDataset]+"_R="+Form("%.1f", jetRadius)+"_PtCutComp_Eta"+jetFinderQaHistType[iJetFinderQaType]);

  // TString textContext(contextDatasetCompAndRadius(jetRadius, ""));
  TString textContext(contextCustomTwoFields(*texDatasetsComparisonCommonDenominator, contextJetRadius(jetRadius), ""));

  const std::array<std::array<float, 2>, 2> drawnWindow = {{{-1, 1}, {-999, -999}}}; // {{xmin, xmax}, {ymin, ymax}}

  Draw_TH1_Histograms_in_one(H1D_jetEta_rebinned, PtCutsLegend, nPtCut, textContext, pdfName, texEtaX, yAxisLabel, texCollisionDataInfo, drawnWindow, "");
}


void Draw_Phi_PtCutComparison(float jetRadius, int iDataset, float* PtCuts, int nPtCut, const char options[]) {

  TH3D* H3D_jetRjetPtjetPhi;
  TH1D* H1D_jetPhi[nPtCut];
  TH1D* H1D_jetPhi_rebinned[nPtCut];
  
  H3D_jetRjetPtjetPhi = (TH3D*)((TH3D*)file_O2Analysis_list[iDataset]->Get(analysisWorkflow[iDataset]+"/h3_jet_r_jet_pt_jet_phi"+jetFinderQaHistType[iJetFinderQaType]))->Clone("Draw_Phi_PtCutComparison"+Datasets[iDataset]);
  H3D_jetRjetPtjetPhi->Sumw2();

  TString PtCutsLegend[nPtCut];
  std::stringstream ss;

  int ibinJetRadius = H3D_jetRjetPtjetPhi->GetXaxis()->FindBin(jetRadius+GLOBAL_epsilon);

  TString* yAxisLabel;

  for(int iBinPt = 0; iBinPt < nPtCut; iBinPt++){

    float PtCutLow = PtCuts[iBinPt];
    int ibinPt_low = H3D_jetRjetPtjetPhi->GetYaxis()->FindBin(PtCutLow);
    if (ibinPt_low == 0) 
      cout << "WARNING: Phi_RadiusComparison is counting the underflow with the chosen PtRange" << endl;

    H1D_jetPhi[iBinPt] = (TH1D*)H3D_jetRjetPtjetPhi->ProjectionZ("jetPhi_"+DatasetsNames[iDataset]+Form("%.1f", jetRadius)+Form("%.1f", PtCutLow)+"<pt", ibinJetRadius, ibinJetRadius, ibinPt_low, H3D_jetRjetPtjetPhi->GetYaxis()->GetNbins(), "e");
    H1D_jetPhi_rebinned[iBinPt] = (TH1D*)H1D_jetPhi[iBinPt]->Rebin(2.,"jetPhi_rebinned_"+DatasetsNames[iDataset]+Form("%.1f", jetRadius)+Form("%.1f", PtCutLow)+"<pt");

    if (strstr(options, "normEntries") != NULL) {
      NormaliseYieldToNEntries(H1D_jetPhi_rebinned[iBinPt]);
      yAxisLabel = texJetPhiYield_EntriesNorm;
    }
    if (strstr(options, "normEvents") != NULL) {
      NormaliseYieldToNEvents(H1D_jetPhi_rebinned[iBinPt], GetNEventsSelected_JetFramework(file_O2Analysis_list[iDataset]));
      yAxisLabel = texJetPhiYield_EventNorm;
    }    


    ss << "Pt > " << PtCuts[iBinPt];
    PtCutsLegend[iBinPt] = (TString)ss.str();
    ss.str("");
    ss.clear();
  }

  TString* pdfName = new TString("jet_"+jetType[iJetType]+"_"+jetLevel[iJetLevel]+"_"+DatasetsNames[iDataset]+"_R="+Form("%.1f", jetRadius)+"_PtCutComp_Phi"+jetFinderQaHistType[iJetFinderQaType]);

  // TString textContext(contextDatasetCompAndRadius(jetRadius, ""));
  TString textContext(contextCustomTwoFields(*texDatasetsComparisonCommonDenominator, contextJetRadius(jetRadius), ""));

  Draw_TH1_Histograms_in_one(H1D_jetPhi_rebinned, PtCutsLegend, nPtCut, textContext, pdfName, texPhiX, yAxisLabel, texCollisionDataInfo, drawnWindowAuto, "histWithLine");
}


void Draw_BkgFluctuations_withFit_CentralityProjection(int iDataset, std::array<std::array<float, 2>, 2> drawnWindowZoom) {

  TH2D* H2D_fluctuations_centrality;

  TH1D* H1D_fluctuations[nCentralityBins];
  TH1D* H1D_fluctuations_rebinned[nCentralityBins];

  H2D_fluctuations_centrality = (TH2D*)((TH2D*)file_O2Analysis_list[iDataset]->Get(analysisWorkflow[iDataset]+"/h2_centrality_rhorandomcone"))->Clone("Draw_BkgFluctuations_withFit_CentralityProjection"+Datasets[iDataset]);
  H2D_fluctuations_centrality->Sumw2();

  float FluctuLow = -60;
  float FluctuHigh = 60;
  int ibinFluctu_low  = H2D_fluctuations_centrality->GetYaxis()->FindBin(FluctuLow);
  int ibinFluctu_high = H2D_fluctuations_centrality->GetYaxis()->FindBin(FluctuHigh);
  H2D_fluctuations_centrality->GetYaxis()->SetRange(ibinFluctu_low, ibinFluctu_high);

  int ibinCent_low, ibinCent_high;
  TString CentralityLegend[nCentralityBins];
  std::stringstream ss;
  ss.precision(2);
  TString* yAxisLabel;
  
  ////////////////////////////////// Fit initialisation //////////////////////////////////
  //Fit tools initialisation
  TF1 *gaussInit[nCentralityBins];
  TF1 *gaussFinal[nCentralityBins];
  TF1 *gaussDrawn[nCentralityBins]; // drawn over the full range
  // TF1 *GaussPlusPolynom[nCentralityBins];
  // TF1 *bkgparab[nCentralityBins];
  TFitResultPtr fFitResult[nCentralityBins];
  // Double_t parGaussParab[nbinpT][6]; //parab backround
  double parGaussInit[nCentralityBins][5]; //linear background
  double parGaussFinal[nCentralityBins][5]; //linear background
  double SignalMean[nCentralityBins], SignalStandardDev[nCentralityBins];
  double SignalMeanError[nCentralityBins], SignalStandardDevError[nCentralityBins];

  for(int iCentralityBin = 0; iCentralityBin < nCentralityBins; iCentralityBin++){
    ibinCent_low = H2D_fluctuations_centrality->GetXaxis()->FindBin(arrayCentralityBinning[iCentralityBin]);
    ibinCent_high = H2D_fluctuations_centrality->GetXaxis()->FindBin(arrayCentralityBinning[iCentralityBin+1])-1;
    // cout << "ibinCent_low = " << ibinCent_low << ", ibinCent_high = " << ibinCent_high << endl;
    H1D_fluctuations[iCentralityBin] = (TH1D*)H2D_fluctuations_centrality->ProjectionY("bkgFluctuationCentrality_"+Datasets[iDataset]+"_@cent["+Form("%.1d", ibinCent_low)+","+Form("%.1d", ibinCent_high)+"]", ibinCent_low, ibinCent_high, "e"); // e computes the new errors

    H1D_fluctuations_rebinned[iCentralityBin] = (TH1D*)H1D_fluctuations[iCentralityBin]->Rebin(1.,"bkgFluctuationCentrality_rebinned_"+Datasets[iDataset]+"_@cent["+Form("%.1d", ibinCent_low)+","+Form("%.1d", ibinCent_high)+"]");

    // if (strstr(options, "normEntries") != NULL) {
    NormaliseYieldToNEntries(H1D_fluctuations_rebinned[iCentralityBin]);
    yAxisLabel = texEntriesNorm_BkgFluctuationYield;
    // }
    // if (strstr(options, "normEvents") != NULL) {
    //   NormaliseYieldToNEvents(H1D_fluctuations_rebinned[iCentralityBin], GetNEventsSelected_JetFramework(file_O2Analysis_list[iDataset]));
    //   yAxisLabel = texCollNorm_BkgFluctuationYield;
    // }
    // if (strstr(options, "normEventsCentrality") != NULL) {
    //   NormaliseYieldToNEvents(H1D_fluctuations_rebinned[iCentralityBin], GetNEventsSelectedCentrality_JetFramework(file_O2Analysis_list[iDataset], arrayCentralityBinning[iCentralityBin], arrayCentralityBinning[iCentralityBin+1], trainId));
    //   yAxisLabel = texCollNorm_BkgFluctuationYield_CentWindow;
    // }

    ////////////////////////////////////////////////////////////////////
    //////////////////////////// Fit start /////////////////////////////
    ////////////////////////////////////////////////////////////////////

    double xHistMax = H1D_fluctuations_rebinned[iCentralityBin]->GetBinContent(H1D_fluctuations_rebinned[iCentralityBin]->GetMaximumBin());
    double yHistMax = H1D_fluctuations_rebinned[iCentralityBin]->GetBinContent(H1D_fluctuations_rebinned[iCentralityBin]->GetMaximumBin());
    
    // gaussInit[iCentralityBin] = new TF1("gaussInit_"+Datasets[iDataset]+"_@cent["+Form("%.1d", ibinCent_low)+","+Form("%.1d", ibinCent_high)+"]", "gaus", H1D_fluctuations_rebinned[iCentralityBin]->GetBinContent(H1D_fluctuations_rebinned[iCentralityBin]->GetMinimumBin()), H1D_fluctuations_rebinned[iCentralityBin]->GetBinContent(H1D_fluctuations_rebinned[iCentralityBin]->GetMaximumBin())); // change -40 60 to values actually taken from histogram automatically
    gaussInit[iCentralityBin] = new TF1("gaussInit_"+Datasets[iDataset]+"_@cent["+Form("%.1d", ibinCent_low)+","+Form("%.1d", ibinCent_high)+"]", "gaus", H1D_fluctuations_rebinned[iCentralityBin]->GetXaxis()->GetXmin(), H1D_fluctuations_rebinned[iCentralityBin]->GetXaxis()->GetXmax());
    gaussInit[iCentralityBin]->SetParName(0, "norm");
    gaussInit[iCentralityBin]->SetParName(1, "mean");
    gaussInit[iCentralityBin]->SetParName(2, "sigma");
    gaussInit[iCentralityBin]->SetParameters(1., 0, 20);
    // gaussInit[iCentralityBin]->SetParLimits(0, 0., 1.1 * yHistMax);
    // gaussInit[iCentralityBin]->SetParLimits(1, -10, 10);
    // gaussInit[iCentralityBin]->SetParLimits(2, 0.1, 100);

    fFitResult[iCentralityBin] = H1D_fluctuations_rebinned[iCentralityBin]->Fit(gaussInit[iCentralityBin], "R0QL"); // P: Use Pearson chi-square method, using expected errors instead of the observed one given by TH1::GetBinError (default case). The expected error is instead estimated from the square-root of the bin function value. (WL for weithged likelihood is currently bugged in root, the fit crashes)

    gaussInit[iCentralityBin]->GetParameters(&parGaussInit[iCentralityBin][0]);


    gaussFinal[iCentralityBin] = new TF1("gaussFinal_"+Datasets[iDataset]+"_@cent["+Form("%.1d", ibinCent_low)+","+Form("%.1d", ibinCent_high)+"]", "gaus", -3*parGaussInit[iCentralityBin][2], parGaussInit[iCentralityBin][1]+0.5*parGaussInit[iCentralityBin][2]);
    gaussFinal[iCentralityBin]->SetParName(0, "norm");
    gaussFinal[iCentralityBin]->SetParName(1, "mean");
    gaussFinal[iCentralityBin]->SetParName(2, "sigma");
    gaussFinal[iCentralityBin]->SetParameters(parGaussInit[iCentralityBin][0], parGaussInit[iCentralityBin][1], parGaussInit[iCentralityBin][2]);
    // gaussFinal[iCentralityBin]->SetParLimits(0, 0., 1.1*yHistMax);
    // gaussFinal[iCentralityBin]->SetParLimits(1, -10, 10);
    // gaussFinal[iCentralityBin]->SetParLimits(2, 0.1, 100);

    fFitResult[iCentralityBin] = H1D_fluctuations_rebinned[iCentralityBin]->Fit(gaussFinal[iCentralityBin], "R0QP"); // P: Use Pearson chi-square method, using expected errors instead of the observed one given by TH1::GetBinError (default case). The expected error is instead estimated from the square-root of the bin function value. (WL for weithged likelihood is currently bugged in root, the fit crashes)
    // gauss[iCentralityBin]->Draw("same");

    gaussFinal[iCentralityBin]->GetParameters(&parGaussFinal[iCentralityBin][0]);



    SignalMean[iCentralityBin] = parGaussFinal[iCentralityBin][1];
    SignalStandardDev[iCentralityBin] = parGaussFinal[iCentralityBin][2];

    SignalMeanError[iCentralityBin] = gaussFinal[iCentralityBin]->GetParError(1);
    SignalStandardDevError[iCentralityBin] = gaussFinal[iCentralityBin]->GetParError(2);

    // cout << "iCentralityBin = " << iCentralityBin << ", SignalMean = " << SignalMean[iCentralityBin] << ", SignalStandardDev = " << SignalStandardDev[iCentralityBin] << ", SignalAmplitude = " << gaussFinal[iCentralityBin]->GetParameter(0) << endl;
    // cout << "                 " << "  SignalMeanError = " << SignalMeanError[iCentralityBin] << ", SignalStandardDevError = " << SignalStandardDevError[iCentralityBin] << ", SignalAmplitudeError = " << gaussFinal[iCentralityBin]->GetParError(0) << endl;

    gaussDrawn[iCentralityBin] = new TF1("gaussDrawn_"+Datasets[iDataset]+"_@cent["+Form("%.1d", ibinCent_low)+","+Form("%.1d", ibinCent_high)+"]", "gaus", H1D_fluctuations_rebinned[iCentralityBin]->GetXaxis()->GetXmin(), H1D_fluctuations_rebinned[iCentralityBin]->GetXaxis()->GetXmax()); // change -40 60 to values actually taken from histogram automatically
    gaussDrawn[iCentralityBin]->SetParameters(parGaussFinal[iCentralityBin][0], parGaussFinal[iCentralityBin][1], parGaussFinal[iCentralityBin][2]);


    ss << "" << arrayCentralityBinning[iCentralityBin] << "-" << arrayCentralityBinning[iCentralityBin+1] << "%; #sigma = " << SignalStandardDev[iCentralityBin] << ", <#delta #it{p}_{T}> = " << SignalMean[iCentralityBin];
    CentralityLegend[iCentralityBin] = (TString)ss.str();
    ss.str("");
    ss.clear();
  }

  TString textContext(contextCustomOneField(*texDatasetsComparisonCommonDenominator, ""));

  // ------ plot the bkg fluctuation distributions with the fits ------ //

  TString* pdfName2 = new TString("jet_"+jetType[iJetType]+"_"+jetLevel[iJetLevel]+"_CentralityComp_"+DatasetsNames[iDataset]+"_BkgFluctuations");
  Draw_TH1_Histograms_in_one(H1D_fluctuations_rebinned, CentralityLegend, nCentralityBins, textContext, pdfName2, texBkgFluctuationRandomCone, yAxisLabel, texCollisionDataInfo, drawnWindowAuto, "fit", gaussDrawn); 
  TString* pdfName = new TString("jet_"+jetType[iJetType]+"_"+jetLevel[iJetLevel]+"_CentralityComp_"+DatasetsNames[iDataset]+"_BkgFluctuations_logy");
  Draw_TH1_Histograms_in_one(H1D_fluctuations_rebinned, CentralityLegend, nCentralityBins, textContext, pdfName, texBkgFluctuationRandomCone, yAxisLabel, texCollisionDataInfo, drawnWindowAuto, "logy,fit", gaussDrawn);

  // zoom around 0
  // for(int iCentralityBin = 0; iCentralityBin < nCentralityBins; iCentralityBin++){
  //   H1D_fluctuations_rebinned[iCentralityBin]->GetXaxis()->SetRangeUser(zoomX[0], zoomX[1]);
  // }
  TString* pdfName2_zoom = new TString("jet_"+jetType[iJetType]+"_"+jetLevel[iJetLevel]+"_CentralityComp_"+DatasetsNames[iDataset]+"_BkgFluctuations_logy_zoom");
  Draw_TH1_Histograms_in_one(H1D_fluctuations_rebinned, CentralityLegend, nCentralityBins, textContext, pdfName2_zoom, texBkgFluctuationRandomCone, yAxisLabel, texCollisionDataInfo, drawnWindowZoom, "logy,fit", gaussDrawn); // replace ", texCollisionDataInfo, drawnWindowAuto, " with ", drawnWindowAuto, " with drawnWindowAuto defined as something like {-999, -999}, and in hist drawing functino I chek for it, and if it's set I draw as set,else I do as usual


  // ------ plot the std deviation as f° of centrality ------ //

  TH1D* H1D_SigmaBkgFluct_asFunctionOf_Centrality = new TH1D("H1D_SigmaBkgFluct_asFunctionOf_Centrality", "H1D_SigmaBkgFluct_asFunctionOf_Centrality", nCentralityBins, arrayCentralityBinning);
  TH1D* H1D_MeanBkgFluct_asFunctionOf_Centrality = new TH1D("H1D_MeanBkgFluct_asFunctionOf_Centrality", "H1D_MeanBkgFluct_asFunctionOf_Centrality", nCentralityBins, arrayCentralityBinning);
  H1D_SigmaBkgFluct_asFunctionOf_Centrality->Sumw2();
  H1D_MeanBkgFluct_asFunctionOf_Centrality->Sumw2();

  for(int iCentralityBin = 0; iCentralityBin < nCentralityBins; iCentralityBin++){
    H1D_SigmaBkgFluct_asFunctionOf_Centrality->SetBinContent(iCentralityBin+1, SignalStandardDev[iCentralityBin]);
    H1D_SigmaBkgFluct_asFunctionOf_Centrality->SetBinError(iCentralityBin+1, SignalStandardDevError[iCentralityBin]);
    H1D_MeanBkgFluct_asFunctionOf_Centrality->SetBinContent(iCentralityBin+1, SignalMean[iCentralityBin]);
    H1D_MeanBkgFluct_asFunctionOf_Centrality->SetBinError(iCentralityBin+1, SignalMeanError[iCentralityBin]);
  }

  TString* pdfName3 = new TString("jet_"+jetType[iJetType]+"_"+jetLevel[iJetLevel]+"_CentralityComp_"+DatasetsNames[iDataset]+"_BkgFluctuationsSigma_asFunctionOf_Centrality");
  Draw_TH1_Histogram(H1D_SigmaBkgFluct_asFunctionOf_Centrality, textContext, pdfName3, texCentrality, texSigmaBkgFluctFit, texCollisionDataInfo, drawnWindowAuto, "");
  TString* pdfName4 = new TString("jet_"+jetType[iJetType]+"_"+jetLevel[iJetLevel]+"_CentralityComp_"+DatasetsNames[iDataset]+"_BkgFluctuationsMean_asFunctionOf_Centrality");
  Draw_TH1_Histogram(H1D_MeanBkgFluct_asFunctionOf_Centrality, textContext, pdfName4, texCentrality, texMeanBkgFluctFit, texCollisionDataInfo, drawnWindowAuto, "");

  cout << "Draw_BkgFluctuations_withFit_CentralityProjection -------- replace P and L in fits with WL once the installed root version has the released bugfix for it" << endl;

}

void Draw_Rho_withFit_NTracksProjection(int iDataset) { /// should be bkgfluct vs ntracks 

  TH2D* H2D_fluctuations_centrality;

  TH1D* H1D_fluctuations[nTracksBins];
  TH1D* H1D_fluctuations_rebinned[nTracksBins];

  H2D_fluctuations_centrality = (TH2D*)((TH2D*)file_O2Analysis_list[iDataset]->Get(analysisWorkflow[iDataset]+"/h2_ntracks_rho"))->Clone("Draw_Rho_withFit_NTracksProjection"+Datasets[iDataset]);
  H2D_fluctuations_centrality->Sumw2();

  // float FluctuLow = -60;
  // float FluctuHigh = 60;
  // int ibinFluctu_low  = H2D_fluctuations_centrality->GetYaxis()->FindBin(FluctuLow);
  // int ibinFluctu_high = H2D_fluctuations_centrality->GetYaxis()->FindBin(FluctuHigh);
  // H2D_fluctuations_centrality->GetYaxis()->SetRange(ibinFluctu_low, ibinFluctu_high);

  int ibinCent_low, ibinCent_high;
  TString CentralityLegend[nTracksBins];
  std::stringstream ss;
  ss.precision(2);
  TString* yAxisLabel;

  ////////////////////////////////// Fit initialisation //////////////////////////////////
  //Fit tools initialisation
  TF1 *gaussInit[nTracksBins];
  TF1 *gaussFinal[nTracksBins];
  TF1 *gaussDrawn[nTracksBins]; // drawn over the full range
  // TF1 *GaussPlusPolynom[nTracksBins];
  // TF1 *bkgparab[nTracksBins];
  TFitResultPtr fFitResult[nTracksBins];
  // Double_t parGaussParab[nbinpT][6]; //parab backround
  double parGaussInit[nTracksBins][5]; //linear background
  double parGaussFinal[nTracksBins][5]; //linear background
  double SignalMean[nTracksBins], SignalStandardDev[nTracksBins];
  double SignalMeanError[nTracksBins], SignalStandardDevError[nTracksBins];

  for(int iTracksBin = 0; iTracksBin < nTracksBins; iTracksBin++){
    ibinCent_low = H2D_fluctuations_centrality->GetXaxis()->FindBin(arrayNTracksBinning[iTracksBin]);
    ibinCent_high = H2D_fluctuations_centrality->GetXaxis()->FindBin(arrayNTracksBinning[iTracksBin+1])-1;
    // cout << "ibinCent_low = " << ibinCent_low << ", ibinCent_high = " << ibinCent_high << endl;
    H1D_fluctuations[iTracksBin] = (TH1D*)H2D_fluctuations_centrality->ProjectionY("bkgFluctuationCentrality_"+Datasets[iDataset]+"_@cent["+Form("%.1d", ibinCent_low)+","+Form("%.1d", ibinCent_high)+"]", ibinCent_low, ibinCent_high, "e"); // e computes the new errors

    H1D_fluctuations_rebinned[iTracksBin] = (TH1D*)H1D_fluctuations[iTracksBin]->Rebin(1.,"bkgFluctuationCentrality_rebinned_"+Datasets[iDataset]+"_@cent["+Form("%.1d", ibinCent_low)+","+Form("%.1d", ibinCent_high)+"]");

    // if (strstr(options, "normEntries") != NULL) {
    NormaliseYieldToNEntries(H1D_fluctuations_rebinned[iTracksBin]);
    yAxisLabel = texEntriesNorm_BkgFluctuationYield;
    // }
    // if (strstr(options, "normEvents") != NULL) {
    //   NormaliseYieldToNEvents(H1D_fluctuations_rebinned[iTracksBin], GetNEventsSelected_JetFramework(file_O2Analysis_list[iDataset]));
    //   yAxisLabel = texCollNorm_BkgFluctuationYield;
    // }
    // if (strstr(options, "normEventsCentrality") != NULL) {
    //   NormaliseYieldToNEvents(H1D_fluctuations_rebinned[iTracksBin], GetNEventsSelectedCentrality_JetFramework(file_O2Analysis_list[iDataset], arrayCentralityBinning[iTracksBin], arrayCentralityBinning[iTracksBin+1], trainId));
    //   yAxisLabel = texCollNorm_BkgFluctuationYield_CentWindow;
    // }

    ////////////////////////////////////////////////////////////////////
    //////////////////////////// Fit start /////////////////////////////
    ////////////////////////////////////////////////////////////////////

    double xHistMax = H1D_fluctuations_rebinned[iTracksBin]->GetBinContent(H1D_fluctuations_rebinned[iTracksBin]->GetMaximumBin());
    double yHistMax = H1D_fluctuations_rebinned[iTracksBin]->GetBinContent(H1D_fluctuations_rebinned[iTracksBin]->GetMaximumBin());
    
    // gaussInit[iTracksBin] = new TF1("gaussInit_"+Datasets[iDataset]+"_@cent["+Form("%.1d", ibinCent_low)+","+Form("%.1d", ibinCent_high)+"]", "gaus", H1D_fluctuations_rebinned[iTracksBin]->GetBinContent(H1D_fluctuations_rebinned[iTracksBin]->GetMinimumBin()), H1D_fluctuations_rebinned[iTracksBin]->GetBinContent(H1D_fluctuations_rebinned[iCentralityBin]->GetMaximumBin())); // change -40 60 to values actually taken from histogram automatically
    gaussInit[iTracksBin] = new TF1("gaussInit_"+Datasets[iDataset]+"_@cent["+Form("%.1d", ibinCent_low)+","+Form("%.1d", ibinCent_high)+"]", "gaus", H1D_fluctuations_rebinned[iTracksBin]->GetXaxis()->GetXmin(), H1D_fluctuations_rebinned[iTracksBin]->GetXaxis()->GetXmax());
    gaussInit[iTracksBin]->SetParName(0, "norm");
    gaussInit[iTracksBin]->SetParName(1, "mean");
    gaussInit[iTracksBin]->SetParName(2, "sigma");
    gaussInit[iTracksBin]->SetParameters(1., 0, 20);
    // gaussInit[iTracksBin]->SetParLimits(0, 0., 1.1 * yHistMax);
    // gaussInit[iTracksBin]->SetParLimits(1, -10, 10);
    // gaussInit[iTracksBin]->SetParLimits(2, 0.1, 100);

    fFitResult[iTracksBin] = H1D_fluctuations_rebinned[iTracksBin]->Fit(gaussInit[iTracksBin], "R0QL"); // P: Use Pearson chi-square method, using expected errors instead of the observed one given by TH1::GetBinError (default case). The expected error is instead estimated from the square-root of the bin function value. (WL for weithged likelihood is currently bugged in root, the fit crashes)

    gaussInit[iTracksBin]->GetParameters(&parGaussInit[iTracksBin][0]);


    gaussFinal[iTracksBin] = new TF1("gaussFinal_"+Datasets[iDataset]+"_@cent["+Form("%.1d", ibinCent_low)+","+Form("%.1d", ibinCent_high)+"]", "gaus", -3*parGaussInit[iTracksBin][2], parGaussInit[iTracksBin][1]+0.5*parGaussInit[iTracksBin][2]);
    gaussFinal[iTracksBin]->SetParName(0, "norm");
    gaussFinal[iTracksBin]->SetParName(1, "mean");
    gaussFinal[iTracksBin]->SetParName(2, "sigma");
    gaussFinal[iTracksBin]->SetParameters(parGaussInit[iTracksBin][0], parGaussInit[iTracksBin][1], parGaussInit[iTracksBin][2]);
    // gaussFinal[iTracksBin]->SetParLimits(0, 0., 1.1*yHistMax);
    // gaussFinal[iTracksBin]->SetParLimits(1, -10, 10);
    // gaussFinal[iTracksBin]->SetParLimits(2, 0.1, 100);

    fFitResult[iTracksBin] = H1D_fluctuations_rebinned[iTracksBin]->Fit(gaussFinal[iTracksBin], "R0QP"); // P: Use Pearson chi-square method, using expected errors instead of the observed one given by TH1::GetBinError (default case). The expected error is instead estimated from the square-root of the bin function value. (WL for weithged likelihood is currently bugged in root, the fit crashes)
    // gauss[iTracksBin]->Draw("same");

    gaussFinal[iTracksBin]->GetParameters(&parGaussFinal[iTracksBin][0]);



    SignalMean[iTracksBin] = parGaussFinal[iTracksBin][1];
    SignalStandardDev[iTracksBin] = parGaussFinal[iTracksBin][2];

    SignalMeanError[iTracksBin] = gaussFinal[iTracksBin]->GetParError(1);
    SignalStandardDevError[iTracksBin] = gaussFinal[iTracksBin]->GetParError(2);

    cout << "iTracksBin = " << iTracksBin << ", SignalMean = " << SignalMean[iTracksBin] << ", SignalStandardDev = " << SignalStandardDev[iTracksBin] << ", SignalAmplitude = " << gaussFinal[iTracksBin]->GetParameter(0) << endl;
    cout << "                 " << "  SignalMeanError = " << SignalMeanError[iTracksBin] << ", SignalStandardDevError = " << SignalStandardDevError[iTracksBin] << ", SignalAmplitudeError = " << gaussFinal[iTracksBin]->GetParError(0) << endl;

    gaussDrawn[iTracksBin] = new TF1("gaussDrawn_"+Datasets[iDataset]+"_@cent["+Form("%.1d", ibinCent_low)+","+Form("%.1d", ibinCent_high)+"]", "gaus", H1D_fluctuations_rebinned[iTracksBin]->GetXaxis()->GetXmin(), H1D_fluctuations_rebinned[iTracksBin]->GetXaxis()->GetXmax()); // change -40 60 to values actually taken from histogram automatically
    gaussDrawn[iTracksBin]->SetParameters(parGaussFinal[iTracksBin][0], parGaussFinal[iTracksBin][1], parGaussFinal[iTracksBin][2]);


    ss << "" << arrayNTracksBinning[iTracksBin] << "-" << arrayNTracksBinning[iTracksBin+1] << "%; #sigma = " << SignalStandardDev[iTracksBin] << ", <#delta #it{p}_{T}> = " << SignalMean[iTracksBin];
    CentralityLegend[iTracksBin] = (TString)ss.str();
    ss.str("");
    ss.clear();
  }

  TString textContext(contextCustomOneField(*texDatasetsComparisonCommonDenominator, ""));

  // ------ plot the bkg fluctuation distributions with the fits ------ //
  TString* pdfName2 = new TString("jet_"+jetType[iJetType]+"_"+jetLevel[iJetLevel]+"_nTracksComp_"+DatasetsNames[iDataset]+"_Rho");
  Draw_TH1_Histograms_in_one(H1D_fluctuations_rebinned, CentralityLegend, nTracksBins, textContext, pdfName2, texBkgFluctuationRandomCone, yAxisLabel, texCollisionDataInfo, drawnWindowAuto, "fit", gaussDrawn);
  // TString* pdfName = new TString("jet_"+jetType[iJetType]+"_"+jetLevel[iJetLevel]+"_nTracksComp_"+Datasets[iDataset]+"_Rho_logy");
  // Draw_TH1_Histograms_in_one(H1D_fluctuations_rebinned, CentralityLegend, nTracksBins, textContext, pdfName, texBkgFluctuationRandomCone, yAxisLabel, texCollisionDataInfo, drawnWindowAuto, "logy,fit", gaussDrawn);

  // // ------ plot the std deviation as f° of centrality ------ //

  // TH1D* H1D_SigmaBkgFluct_asFunctionOf_Centrality = new TH1D("H1D_SigmaBkgFluct_asFunctionOf_Centrality", "H1D_SigmaBkgFluct_asFunctionOf_Centrality", nTracksBins, arrayNTracksBinning);
  // TH1D* H1D_MeanBkgFluct_asFunctionOf_Centrality = new TH1D("H1D_MeanBkgFluct_asFunctionOf_Centrality", "H1D_MeanBkgFluct_asFunctionOf_Centrality", nTracksBins, arrayNTracksBinning);
  // H1D_SigmaBkgFluct_asFunctionOf_Centrality->Sumw2();
  // H1D_MeanBkgFluct_asFunctionOf_Centrality->Sumw2();

  // for(int iTracksBin = 0; iTracksBin < nTracksBins; iTracksBin++){
  //   H1D_SigmaBkgFluct_asFunctionOf_Centrality->SetBinContent(iTracksBin+1, SignalStandardDev[iTracksBin]);
  //   H1D_SigmaBkgFluct_asFunctionOf_Centrality->SetBinError(iTracksBin+1, SignalStandardDevError[iTracksBin]);
  //   H1D_MeanBkgFluct_asFunctionOf_Centrality->SetBinContent(iTracksBin+1, SignalMean[iTracksBin]);
  //   H1D_MeanBkgFluct_asFunctionOf_Centrality->SetBinError(iTracksBin+1, SignalMeanError[iTracksBin]);
  // }

  // TString* pdfName3 = new TString("jet_"+jetType[iJetType]+"_"+jetLevel[iJetLevel]+"_nTracksComp_"+Datasets[iDataset]+"_RhoSigma_asFunctionOf_Centrality");
  // Draw_TH1_Histogram(H1D_SigmaBkgFluct_asFunctionOf_Centrality, textContext, pdfName3, texCentrality, texSigmaBkgFluctFit, texCollisionDataInfo, drawnWindowAuto, "");
  // TString* pdfName4 = new TString("jet_"+jetType[iJetType]+"_"+jetLevel[iJetLevel]+"_nTracksComp_"+Datasets[iDataset]+"_RhoMean_asFunctionOf_Centrality");
  // Draw_TH1_Histogram(H1D_MeanBkgFluct_asFunctionOf_Centrality, textContext, pdfName4, texCentrality, texMeanBkgFluctFit, texCollisionDataInfo, drawnWindowAuto, "");

  // cout << "Draw_BkgFluctuations_withFit_CentralityProjection -------- replace P and L in fits with WL once the installed root version has the released bugfix for it" << endl;


}

void Draw_Pt_Run2Run3Comparison_0010Cent_R040(int iDataset, const char options[]) {

  TString VarRun3;
  TString VarRun2;
  TString correctedStatus;
  if (strstr(options, "jetUncorrected") != NULL){
    VarRun3 = "h3_jet_r_jet_pt_centrality";
    VarRun2 = "histJetPt_0";
    correctedStatus = "";
  }
  if (strstr(options, "jetCorrected") != NULL){
    VarRun3 = "h3_jet_r_jet_pt_jet_eta_rhoareasubtracted";
    VarRun2 = "histJetCorrPt_0";
    correctedStatus = "_corrected";
  }

  int ibinCent_low, ibinCent_high;
  float jetRadius = 0.4;

  // O2Physics Run 3

  TH3D* H3D_jetRjetVarCent;
  TH1D* H1D_jetVar;
  TH1D* H1D_jetVar_rebinned;

  H3D_jetRjetVarCent = (TH3D*)((TH3D*)file_O2Analysis_list[iDataset]->Get(analysisWorkflow[iDataset]+"/"+VarRun3))->Clone("Draw_Pt_CentralityComparison"+Datasets[iDataset]+Form("%.1f", jetRadius));
  H3D_jetRjetVarCent->Sumw2();

  int ibinJetRadius = H3D_jetRjetVarCent->GetXaxis()->FindBin(jetRadius+GLOBAL_epsilon);

  ibinCent_low = H3D_jetRjetVarCent->GetZaxis()->FindBin((double)0 + GLOBAL_epsilon);
  ibinCent_high = H3D_jetRjetVarCent->GetZaxis()->FindBin((double)10 - GLOBAL_epsilon)-1;
  H1D_jetVar = (TH1D*)H3D_jetRjetVarCent->ProjectionY("jetPt_"+Datasets[iDataset]+Form("%.1f", jetRadius)+"_@cent["+Form("%.1d", ibinCent_low)+","+Form("%.1d", ibinCent_high)+"]", ibinJetRadius, ibinJetRadius, ibinCent_low, ibinCent_high, "e");
  H1D_jetVar_rebinned = (TH1D*)H1D_jetVar->Rebin(1.,"jetPt_rebinned_"+Datasets[iDataset]+Form("%.1f", jetRadius)+"_@cent["+Form("%.1d", ibinCent_low)+","+Form("%.1d", ibinCent_high)+"]");

  // NormaliseYieldToNEntries(H1D_jetVar_rebinned[iRadius]);
  // NormaliseYieldToNEvents(H1D_jetVar_rebinned[iCentralityBin], GetNEventsSelected_JetFramework(file_O2Analysis_list[iDataset]));
  // NormaliseYieldToNEvents(H1D_jetVar_rebinned, GetNEventsSelectedCentrality_JetFramework(file_O2Analysis_list[iDataset], arrayCentralityBinning, arrayCentralityBinning, trainId));
  // NormaliseYieldToNEvents(H1D_jetVar_rebinned[iCentralityBin], GetNEventsSelected_JetFramework(file_O2Analysis_list[iDataset]));

  NormaliseYieldToNEvents(H1D_jetVar_rebinned, GetNEventsSelectedCentrality_JetFramework(file_O2Analysis_list[iDataset], 0, 10, trainId));

  // AliPhysics Run 2

  TH1D* H1D_run2_jetVar;
  TH1D* H1D_run2_jetVar_rebinned;
  H1D_run2_jetVar = (TH1D*)((TObject*)((TObject*)file_AliAnalysis->Get("JetCore_JetCore_57_2050_04_histos"))->FindObject("Jet_AKTChargedR040_tracks_pT0150_pt_scheme"))->FindObject(VarRun2);
  std::vector<double> O2H1DPtbinsVector = GetTH1Bins(H1D_jetVar_rebinned);
  double* O2ptBins = &O2H1DPtbinsVector[0];
  H1D_run2_jetVar_rebinned = (TH1D*)H1D_run2_jetVar->Rebin(H1D_jetVar_rebinned->GetNbinsX(), "H1D_run2_jetVar_rebinned", O2ptBins);

  int nEvents_CentWindow;
  TH1D* H1D_Centrality_Run2= (TH1D*)((TObject*)file_AliAnalysis->Get("JetCore_JetCore_57_2050_04_histos"))->FindObject("fHistCentrality");
  int iBinCent_low = H1D_Centrality_Run2->GetXaxis()->FindBin(0 + GLOBAL_epsilon);
  int iBinCent_high = H1D_Centrality_Run2->GetXaxis()->FindBin(10 - GLOBAL_epsilon);
  nEvents_CentWindow = H1D_Centrality_Run2->Integral(iBinCent_low, iBinCent_high);

  H1D_run2_jetVar_rebinned->Scale(1./nEvents_CentWindow,"width"); // If option contains "width" the bin contents and errors are divided by the bin width.

  // plotting
  TH1D* hist_list[2] = {H1D_jetVar_rebinned, H1D_run2_jetVar_rebinned};
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



  TString* pdfName = new TString("jet_Run2_vs_Run3_"+DatasetsNames[iDataset]+"_Pt"+correctedStatus);

  TString textContext(contextCustomOneField(RunCompLegend[0]+" vs "+RunCompLegend[1], ""));

  Draw_TH1_Histograms_in_one(hist_list, RunCompLegend, 2, textContext, pdfName, texPtX, texJetPtYield_EventNorm_CentWindow, texCollisionDataInfo, drawnWindowAuto, "logy");
  TString* pdfName2 = new TString("jet_Run2_vs_Run3_"+DatasetsNames[iDataset]+"_Pt_logx");
  Draw_TH1_Histograms_in_one(hist_list, RunCompLegend, 2, textContext, pdfName2, texPtX, texJetPtYield_EventNorm_CentWindow, texCollisionDataInfo, drawnWindowAuto, "logx,logy");
}

void Draw_Rho_vs_SelectedMultiplicity_DatasetCompRatio() {

  TH2D* H2D_rhoMult[nDatasets-1];
  TH2D* H2D_rhoMult_rebinned[nDatasets-1];
  TH2D* H2D_rhoMult_ratios[nDatasets-1];

  TH2D* H2D_rhoMult_ref = (TH2D*)((TH2D*)file_O2Analysis_list[0]->Get(analysisWorkflow[0]+"/h2_ntracks_rho"))->Clone("Draw_Rho_vs_Mult_ref"+Datasets[0]);
  H2D_rhoMult_ref->Sumw2();

  TH2D* H2D_rhoMult_ref_rebinned = (TH2D*)H2D_rhoMult_ref->RebinX(1.,"H2D_rhoMult_ref_rebinned"+Datasets[0]);

  bool divideSuccess = false;
  for(int iDataset = 0; iDataset < nDatasets-1; iDataset++){
    H2D_rhoMult[iDataset] = (TH2D*)((TH2D*)file_O2Analysis_list[iDataset+1]->Get(analysisWorkflow[iDataset+1]+"/h2_ntracks_rho"))->Clone("Draw_Rho_vs_SelectedMultiplicity_DatasetCompRatio"+Datasets[iDataset]);
    H2D_rhoMult[iDataset]->Sumw2();

    H2D_rhoMult_rebinned[iDataset] = (TH2D*)H2D_rhoMult[iDataset]->RebinX(1.,"H2D_rhoMult_ratio_rebinned"+Datasets[iDataset]);
    // H2D_rhoMult_rebinned[iDataset]->GetXaxis()->SetRange(0,H2D_rhoMult_rebinned[iDataset]->FindLastBinAbove(1, 1)); //(asks for the last bin on the x axis (axis number 1) to have strictly more than 1 entry)
    // H2D_rhoMult_rebinned[iDataset]->GetYaxis()->SetRange(1, H2D_rhoMult_rebinned[iDataset]->FindLastBinAbove(1, 2));

    divideSuccess = H2D_rhoMult_ratios[iDataset]->Divide(H2D_rhoMult_rebinned[iDataset], H2D_rhoMult_ref_rebinned);
    cout << "division of H2 for dataset " << iDataset << " is successful? " << divideSuccess << endl;
  }

  TString* pdfName = new TString("jet_"+jetType[iJetType]+"_"+jetLevel[iJetLevel]+"_DataComp_rho_vs_multiplicitySelected_ratio");

  TString textContext(contextCustomOneField(*texDatasetsComparisonCommonDenominator, ""));

  Draw_TH2_Histograms(H2D_rhoMult_ratios, DatasetsNames, nDatasets-1, textContext, pdfName, texSelectedMultiplicity, texRho, texCollisionDataInfo, drawnWindowAuto, "logz,autoRangeSame");
}

// void Draw_Rho_vs_SelectedMultiplicity_CentralityComp(int iDataset) {

//   TH3D* H3D_rhoMult = (TH3D*)((TH3D*)file_O2Analysis_list[iDataset]->Get(analysisWorkflow[iDataset]+"/h3_centrality_ntracks_rho"))->Clone("Draw_Rho_vs_Mult"+Datasets[iDataset]);
//   H3D_rhoMult->Sumw2();

//   TH2D* H2D_rhoMult[nCentralityBins];
//   TH2D* H2D_rhoMult_rebinned[nCentralityBins];
  
//   int ibinCent_low, ibinCent_high;
//   TString CentralityLegend[nCentralityBins];
//   std::stringstream ss;
//   for(int iCentralityBin = 0; iCentralityBin < nCentralityBins; iCentralityBin++){

//     ibinCent_low = H3D_rhoMult->GetXaxis()->FindBin(arrayCentralityBinning[iCentralityBin]);
//     ibinCent_high = H3D_rhoMult->GetXaxis()->FindBin(arrayCentralityBinning[iCentralityBin+1])-1;

//     H2D_rhoMult[iCentralityBin] = (TH2D*)H3D_rhoMult->ProjectionX("H2D_rhoMult_"+Datasets[iDataset]+Form("%.1f", jetRadius)+"_@cent["+Form("%.1d", ibinCent_low)+","+Form("%.1d", ibinCent_high)+"]", ibinCent_low, ibinCent_high, "e");

//     H2D_rhoMult_rebinned[iCentralityBin] = (TH2D*)H2D_rhoMult[iCentralityBin]->RebinX(1.,"H2D_rhoMult_rebinned"+Datasets[iDataset]+"_@cent["+Form("%.1d", ibinCent_low)+","+Form("%.1d", ibinCent_high)+"]");
//     // H2D_rhoMult_rebinned[iDataset]->GetXaxis()->SetRange(0,H2D_rhoMult_rebinned[iDataset]->FindLastBinAbove(1, 1)); //(asks for the last bin on the x axis (axis number 1) to have strictly more than 1 entry)
//     // H2D_rhoMult_rebinned[iDataset]->GetYaxis()->SetRange(1, H2D_rhoMult_rebinned[iDataset]->FindLastBinAbove(1, 2));

//     ss << "Cent " << arrayCentralityBinning[iCentralityBin] << " - " << arrayCentralityBinning[iCentralityBin+1] << " ";
//     CentralityLegend[iCentralityBin] = (TString)ss.str();
//     ss.str("");
//     ss.clear();
//   }

//   TString* pdfName = new TString("jet_"+jetType[iJetType]+"_"+jetLevel[iJetLevel]+"_CentralityComp_R="+Form("%.1f", jetRadius)+"_"+Datasets[iDataset]+"rho_vs_multiplicitySelected");

//   TString textContext(contextCustomOneField(*texDatasetsComparisonCommonDenominator, ""));

//   Draw_TH2_Histograms(H2D_rhoMult_rebinned, DatasetsNames, nCentralityBins, textContext, pdfName, texSelectedMultiplicity, texRho, texCollisionDataInfo, drawnWindowAuto, "logz,autoRangeSame");
// }


// void Draw_Pt_CentralityAndDatasetComparison(float jetRadius) {

//   TH3D* H3D_jetRjetPtjetCent;
//   TH1D* H1D_jetPt[nCentralityBins];
//   TH1D* H1D_jetPt_rebinned[nCentralityBins];
  
//   H3D_jetRjetPtjetCent = (TH3D*)((TH3D*)file_O2Analysis_list[iDataset]->Get(analysisWorkflow[iDataset]+"/h3_jet_r_jet_pt_centrality"+jetFinderQaHistType[iJetFinderQaType]))->Clone("Draw_Pt_CentralityComparison"+Datasets[iDataset]+Form("%.1f", jetRadius));
//   H3D_jetRjetPtjetCent->Sumw2();

//   int ibinJetRadius = H3D_jetRjetPtjetCent->GetXaxis()->FindBin(jetRadius+GLOBAL_epsilon);
//   int ibinCent_low, ibinCent_high;
//   TString CentralityLegend[nCentralityBins];
//   std::stringstream ss;
//   for(int iCentralityBin = 0; iCentralityBin < nCentralityBins; iCentralityBin++){

//     ibinCent_low = H3D_jetRjetPtjetCent->GetZaxis()->FindBin(arrayCentralityBinning[iCentralityBin]);
//     ibinCent_high = H3D_jetRjetPtjetCent->GetZaxis()->FindBin(arrayCentralityBinning[iCentralityBin+1])-1;
//     H1D_jetPt[iCentralityBin] = (TH1D*)H3D_jetRjetPtjetCent->ProjectionY("jetPt_"+Datasets[iDataset]+Form("%.1f", jetRadius)+"_@cent["+Form("%.1d", ibinCent_low)+","+Form("%.1d", ibinCent_high)+"]", ibinJetRadius, ibinJetRadius, ibinCent_low, ibinCent_high, "e");
//     H1D_jetPt_rebinned[iCentralityBin] = (TH1D*)H1D_jetPt[iCentralityBin]->Rebin(1.,"jetPt_rebinned_"+Datasets[iDataset]+Form("%.1f", jetRadius)+"_@cent["+Form("%.1d", ibinCent_low)+","+Form("%.1d", ibinCent_high)+"]");

//     // NormaliseYieldToNEntries(H1D_jetPt_rebinned[iRadius]);
//     NormaliseYieldToNEvents(H1D_jetPt_rebinned[iCentralityBin], GetNEventsSelectedCentrality_JetFramework(file_O2Analysis_list[iDataset], arrayCentralityBinning[iCentralityBin], arrayCentralityBinning[iCentralityBin+1], trainId));
//     // NormaliseYieldToNEvents(H1D_jetPt_rebinned[iCentralityBin], GetNEventsSelected_JetFramework(file_O2Analysis_list[iDataset]));

//     ss << "Cent " << arrayCentralityBinning[iCentralityBin] << " - " << arrayCentralityBinning[iCentralityBin+1] << " ";
//     CentralityLegend[iCentralityBin] = (TString)ss.str();
//     ss.str("");
//     ss.clear();
//   }

//   TString* pdfName = new TString("jet_"+jetType[iJetType]+"_"+jetLevel[iJetLevel]+"_CentralityComp_R="+Form("%.1f", jetRadius)+"_"+Datasets[iDataset]+"_Pt"+jetFinderQaHistType[iJetFinderQaType]);

//   // TString textContext(contextDatasetCompAndRadius(jetRadius, ""));
//   TString textContext(contextCustomTwoFields(*texDatasetsComparisonCommonDenominator, contextJetRadius(jetRadius), ""));

//   Draw_TH1_Histograms_in_one(H1D_jetPt_rebinned, CentralityLegend, nCentralityBins, textContext, pdfName, texPtX, texJetPtYield_EventNorm_CentWindow, texCollisionDataInfo, drawnWindowAuto, "logy");
//   // Draw_TH1_Histograms_in_one(H1D_jetPt_rebinned, CentralityLegend, nCentralityBins, textContext, pdfName, texPtX, texJetPtYield_EventNorm, texCollisionDataInfo, drawnWindowAuto, "logy");
//   if (iJetFinderQaType == 0) {
//     TString* pdfName2 = new TString("jet_"+jetType[iJetType]+"_"+jetLevel[iJetLevel]+"_CentralityComp_R="+Form("%.1f", jetRadius)+"_"+Datasets[iDataset]+"_Pt"+jetFinderQaHistType[iJetFinderQaType]+"_logx");
//     Draw_TH1_Histograms_in_one(H1D_jetPt_rebinned, CentralityLegend, nCentralityBins, textContext, pdfName2, texPtX, texJetPtYield_EventNorm_CentWindow, texCollisionDataInfo, drawnWindowAuto, "logx,logy");
//   }
// }

void Draw_Pt_RadiusComparison_mcp(int iDataset, float* etaRange) {

  TH3D* H3D_jetRjetPtjetEta;
  TH1D* H1D_jetPt[nRadius];
  TH1D* H1D_jetPt_rebinned[nRadius];
  
  H3D_jetRjetPtjetEta = (TH3D*)((TH3D*)file_O2Analysis_list[iDataset]->Get(analysisWorkflow[iDataset]+"/h3_jet_r_part_jet_pt_part_jet_eta_part"+jetFinderQaHistType[iJetFinderQaType]))->Clone("Draw_Pt_RadiusComparison"+Datasets[iDataset]+Form("%.1f", etaRange[0])+"<eta<"+Form("%.1f", etaRange[1]));
  H3D_jetRjetPtjetEta->Sumw2();

  float EtaCutLow = etaRange[0];
  float EtaCutHigh = etaRange[1];
  int ibinEta_low = H3D_jetRjetPtjetEta->GetZaxis()->FindBin(EtaCutLow);
  int ibinEta_high = H3D_jetRjetPtjetEta->GetZaxis()->FindBin(EtaCutHigh);
  if (ibinEta_low == 0) 
    cout << "WARNING: Pt_RadiusComparison is counting the underflow with the chosen etaRange" << endl;
  if (ibinEta_high == H3D_jetRjetPtjetEta->GetZaxis()->GetNbins()+1) 
    cout << "WARNING: Pt_RadiusComparison is counting the overflow with the chosen etaRange" << endl;
  int ibinJetRadius = 0;

  for(int iRadius = 0; iRadius < nRadius; iRadius++){
    ibinJetRadius = H3D_jetRjetPtjetEta->GetXaxis()->FindBin(arrayRadius[iRadius]+GLOBAL_epsilon);

    H1D_jetPt[iRadius] = (TH1D*)H3D_jetRjetPtjetEta->ProjectionY("jetPt_"+RadiusLegend[iRadius]+Form("%.1f", EtaCutLow)+"<eta<"+Form("%.1f", EtaCutHigh), ibinJetRadius, ibinJetRadius, ibinEta_low, ibinEta_high, "e");
    H1D_jetPt_rebinned[iRadius] = (TH1D*)H1D_jetPt[iRadius]->Rebin(1.,"jetPt_rebinned_"+RadiusLegend[iRadius]+Form("%.1f", EtaCutLow)+"<eta<"+Form("%.1f", EtaCutHigh));

    // NormaliseYieldToNEntries(H1D_jetPt_rebinned[iRadius]);
    NormaliseYieldToNEvents(H1D_jetPt_rebinned[iRadius], GetNEventsSelected_JetFramework(file_O2Analysis_list[iDataset]));
  }

  TString* pdfName = new TString("jet_"+jetType[iJetType]+"_"+jetLevel[iJetLevel]+"_"+DatasetsNames[iDataset]+"_Pt_@eta["+Form("%.1f", EtaCutLow)+","+Form("%.1f", EtaCutHigh)+"]"+jetFinderQaHistType[iJetFinderQaType]+"_mcp");

  // TString textContext(contextDatasetRadiusCompAndVarRange(iDataset, etaRange, "eta"));
  TString textContext(contextCustomThreeFields(*texDatasetsComparisonCommonDenominator, Datasets[iDataset], contextEtaRange(etaRange), ""));

  Draw_TH1_Histograms_in_one(H1D_jetPt_rebinned, RadiusLegend, nRadius, textContext, pdfName, texPtX, texJetPtYield_EventNorm, texCollisionDataInfo, drawnWindowAuto, "logy");
}





void Draw_Ncoll_vs_centrality(const char options[]) {

  TH2D* H2D_centralityCollisions[nDatasets];
  TH1D* H1D_Ncoll[nDatasets];
  TH1D* H1D_Ncoll_rebinned[nDatasets];
  
  for(int iDataset = 0; iDataset < nDatasets; iDataset++){

    H2D_centralityCollisions[iDataset] = (TH2D*)((TH2D*)file_O2Analysis_list[iDataset]->Get(analysisWorkflow[iDataset]+"/h2_centrality_collisions"))->Clone("Draw_Ncoll_vs_centrality"+Datasets[iDataset]);

    H1D_Ncoll[iDataset] = (TH1D*)H2D_centralityCollisions[iDataset]->ProjectionX("Ncoll_"+Datasets[iDataset], 3, 3, "e");
    H1D_Ncoll_rebinned[iDataset] = (TH1D*)H1D_Ncoll[iDataset]->Rebin(10.,"Ncoll_rebinned_"+Datasets[iDataset]);

    H1D_Ncoll_rebinned[iDataset]->GetXaxis()->SetRange(1,H1D_Ncoll_rebinned[iDataset]->GetXaxis()->FindBin(98));
    H1D_Ncoll_rebinned[iDataset]->Scale(1./H1D_Ncoll_rebinned[iDataset]->Integral(1,100));
    // NormaliseYieldToNEntries(H1D_Ncoll_rebinned[iDataset]);
  }

  TString* pdfName = new TString("nCollSelected_vs_centrality_DataComp");

  TString textContext(contextCustomOneField(*texDatasetsComparisonCommonDenominator, ""));

  Draw_TH1_Histograms_in_one(H1D_Ncoll_rebinned, DatasetsNames, nDatasets, textContext, pdfName, texCentrality, texCount, texCollisionDataInfo, drawnWindowAuto, "");

}



void Draw_jet_resolution_MC_PtRangeComparison(int iDataset, double jetRadius, const char options[]) {

  const int nPtRange = 4;
  double ptRange[nPtRange+1] = {10, 20, 40, 80, 200};

  TH3D* H3D_jetPtjetRes = (TH3D*)((TH3D*)file_O2Analysis_list[iDataset]->Get(analysisWorkflow[iDataset]+"/h3_jet_r_jet_pt_tag_jet_pt_base_diff_matchedgeo"))->Clone("Draw_jet_resolution_MC_DatasetComparison");
  TH1D* H1D_jetRes[nPtRange];
  TH1D* H1D_jetRes_rebinned[nPtRange];
  
  TString* yAxisLabel;

  for(int iPtRange = 0; iPtRange < nPtRange; iPtRange++){
    cout << "test" << endl;
    float PtCutLow = ptRange[iPtRange];
    float PtCutHigh = ptRange[iPtRange+1];

    int ibinPt_low = H3D_jetPtjetRes->GetYaxis()->FindBin(PtCutLow);
    int ibinPt_high = H3D_jetPtjetRes->GetYaxis()->FindBin(PtCutHigh);
    int ibinJetRadius = H3D_jetPtjetRes->GetXaxis()->FindBin(jetRadius+GLOBAL_epsilon);

    H1D_jetRes[iPtRange] = (TH1D*)H3D_jetPtjetRes->ProjectionZ("jetRes_"+(TString)Form("%.1f", PtCutLow)+"<eta<"+(TString)Form("%.1f", PtCutHigh), ibinJetRadius, ibinJetRadius, ibinPt_low, ibinPt_high, "e");
    H1D_jetRes_rebinned[iPtRange] = (TH1D*)H1D_jetRes[iPtRange]->Rebin(10.,"jetRes_rebinned_"+(TString)Form("%.1f", PtCutLow)+"<eta<"+(TString)Form("%.1f", PtCutHigh));

    NormaliseYieldToIntegral(H1D_jetRes_rebinned[iPtRange]);
    yAxisLabel = texJetPtYield_EntriesNorm;
    // if (strstr(options, "normEvents") != NULL) {
    //   NormaliseYieldToNEvents(H1D_jetRes_rebinned[iDataset], GetNEventsSelected_JetFramework(file_O2Analysis_list[iDataset]));
    //   yAxisLabel = texJetPtYield_EventNorm;
    // }

    // int ibinJetRes_low = H1D_jetRes_rebinned[iDataset]->GetXaxis()->FindBin(-2);
    // int ibinJetRes_high = H1D_jetRes_rebinned[iDataset]->GetXaxis()->FindBin(1.1);
    // H1D_jetRes_rebinned[iDataset]->GetXaxis()->SetRangeUser(ibinJetRes_low, ibinJetRes_high);
  }
  TString* pdfName = new TString("jet_"+jetType[iJetType]+"_"+jetLevel[iJetLevel]+"_Resolution_PtRangeComp_"+DatasetsNames[iDataset]);

  TString textContext(contextCustomOneField(*texDatasetsComparisonCommonDenominator, ""));

  TString PtRangeNames[nPtRange] = {(TString)"10 < pt < 20 GeV", (TString)"20 < pt < 40 GeV", (TString)"40 < pt < 80 GeV", (TString)"80 < pt < 200 GeV"};

  std::array<std::array<float, 2>, 2> drawnWindow = {{{-2, 1.5}, {-999, -999}}}; // {{xmin, xmax}, {ymin, ymax}}

  Draw_TH1_Histograms_in_one(H1D_jetRes_rebinned, PtRangeNames, nPtRange, textContext, pdfName, texResMC, texJetPtYield_EntriesNorm, texCollisionDataInfo, drawnWindow, "logy");
}


void Count_Nevents_perCentClass(int iDataset, const char options[]) {

  TH2D* H2D_centrality_collisions = (TH2D*)((TH2D*)file_O2Analysis_list[iDataset]->Get(analysisWorkflow[iDataset]+"/h2_centrality_collisions"))->Clone("Count_Nevents_perCentClass");
  
  TString* yAxisLabel;

  const int nCent = 7;
  float centRange[nCent+1] = {-100, 00, 10, 20, 40, 60, 80, 140};

  cout << "dataset " << DatasetsNames[iDataset]<< endl;
  for(int iCent = 0; iCent < nCent; iCent++){
    float centLow = centRange[iCent];
    float centHigh = centRange[iCent+1];

    int ibinCent_low = H2D_centrality_collisions->GetXaxis()->FindBin(centLow);
    int ibinCent_high = H2D_centrality_collisions->GetXaxis()->FindBin(centHigh);

    TH1D* H1D_collisions = H2D_centrality_collisions->ProjectionY("H1D_collisions"+(TString)Datasets[iDataset]+"_@cent["+Form("%.1f", centRange[0])+","+Form("%.1f", centRange[1])+"]", ibinCent_low, ibinCent_high);

    cout << "Centrality " << centLow << "-" << centHigh << "% : nColl total = " << (long int)H1D_collisions->GetBinContent(1) << ", nColl selected = " << (long int)H1D_collisions->GetBinContent(2) << endl;
  }
}





void Draw_Pt_PbPbToPPComparison_HARDCODED(float jetRadius, float* etaRange, const char options[]) {

  TH3D* H3D_jetRjetPtjetEta[nDatasets];
  TH1D* H1D_jetPt[nDatasets];
  TH1D* H1D_jetPt_rebinned[nDatasets];
  
  TH1D* H1D_jetPt_rebinned_ratios[nDatasets];

  float EtaCutLow = etaRange[0];
  float EtaCutHigh = etaRange[1];
  int ibinJetRadius = 0;

  bool divideSuccess = false;

  TString* yAxisLabel;
  // cout << "test0" << endl;
  for(int iDataset = 0; iDataset < nDatasets; iDataset++){
    // cout << "test0.1" << endl;

    if (DatasetsNames[iDataset].EqualTo((TString)"LHC22o pass6")) { // hardcoded here "LHC22o pass7"
      H3D_jetRjetPtjetEta[iDataset] = (TH3D*)((TH3D*)file_O2Analysis_list[iDataset]->Get(analysisWorkflow[iDataset]+"/h3_jet_r_jet_pt_jet_eta"+jetFinderQaHistType[0]))->Clone("Draw_Pt_DatasetComparison"+Datasets[iDataset]+"Radius"+Form("%.1f",jetRadius)+jetRadius+Form("%.1f", etaRange[0])+"<eta<"+Form("%.1f", etaRange[1]));
    } else {
      H3D_jetRjetPtjetEta[iDataset] = (TH3D*)((TH3D*)file_O2Analysis_list[iDataset]->Get(analysisWorkflow[iDataset]+"/h3_jet_r_jet_pt_jet_eta"+jetFinderQaHistType[iJetFinderQaType]))->Clone("Draw_Pt_DatasetComparison"+Datasets[iDataset]+"Radius"+Form("%.1f",jetRadius)+jetRadius+Form("%.1f", etaRange[0])+"<eta<"+Form("%.1f", etaRange[1]));
    }
    // cout << "test0.2" << endl;

    int ibinEta_low = H3D_jetRjetPtjetEta[iDataset]->GetZaxis()->FindBin(EtaCutLow);
    int ibinEta_high = H3D_jetRjetPtjetEta[iDataset]->GetZaxis()->FindBin(EtaCutHigh);
    if (ibinEta_low == 0) 
      cout << "WARNING: Pt_DatasetComparison is counting the underflow with the chosen etaRange" << endl;
    if (ibinEta_high == H3D_jetRjetPtjetEta[iDataset]->GetZaxis()->GetNbins()+1) 
      cout << "WARNING: Pt_DatasetComparison is counting the overflow with the chosen etaRange" << endl;
    ibinJetRadius = H3D_jetRjetPtjetEta[iDataset]->GetXaxis()->FindBin(jetRadius+GLOBAL_epsilon);
    // cout << "test0.3" << endl;
 
    H1D_jetPt[iDataset] = (TH1D*)H3D_jetRjetPtjetEta[iDataset]->ProjectionY("jetPt_"+Datasets[iDataset]+"Radius"+Form("%.1f",jetRadius)+Form("%.1f", EtaCutLow)+"<eta<"+Form("%.1f", EtaCutHigh), ibinJetRadius, ibinJetRadius, ibinEta_low, ibinEta_high, "e");

    if (DatasetsNames[iDataset].EqualTo((TString)"LHC22o pass6")) {
      std::vector<double> O2H1DPtbinsVector = GetTH1Bins(H1D_jetPt_rebinned[0]); // hardcoded here H1D_jetPt_rebinned[0]
      double* O2ptBins = &O2H1DPtbinsVector[0];
      H1D_jetPt_rebinned[iDataset] = (TH1D*)H1D_jetPt[iDataset]->Rebin(H1D_jetPt_rebinned[0]->GetNbinsX(), "jetPt_rebinned_"+Datasets[iDataset]+"Radius"+Form("%.1f",jetRadius)+Form("%.1f", EtaCutLow)+"<eta<"+Form("%.1f", EtaCutHigh), O2ptBins);  // hardcoded here H1D_jetPt_rebinned[0]
    } else {
      H1D_jetPt_rebinned[iDataset] = (TH1D*)H1D_jetPt[iDataset]->Rebin(5.,"jetPt_rebinned_"+Datasets[iDataset]+"Radius"+Form("%.1f",jetRadius)+Form("%.1f", EtaCutLow)+"<eta<"+Form("%.1f", EtaCutHigh));
    }
    
    // cout << "test0.35" << endl;

    if (strstr(options, "normEntries") != NULL) {
      NormaliseYieldToNEntries(H1D_jetPt_rebinned[iDataset]);
      yAxisLabel = texJetPtYield_EntriesNorm;
    }
    int Nevents;
    if (strstr(options, "normEvents") != NULL) {
      if (isDatasetWeighted[iDataset]) {
        Nevents = GetNEventsSelected_JetFramework_weighted(file_O2Analysis_list[iDataset]);
      } else {
        Nevents = GetNEventsSelected_JetFramework(file_O2Analysis_list[iDataset]);
      }
      NormaliseYieldToNEvents(H1D_jetPt_rebinned[iDataset], Nevents);
      yAxisLabel = texJetPtYield_EventNorm;
    }
    // cout << "test0.4" << endl;

    H1D_jetPt_rebinned_ratios[iDataset] = (TH1D*)H1D_jetPt_rebinned[iDataset]->Clone("jetPt_rebinned_ratios"+Datasets[iDataset]+"Radius"+Form("%.1f",jetRadius)+Form("%.1f", EtaCutLow)+"<eta<"+Form("%.1f", EtaCutHigh));
    H1D_jetPt_rebinned_ratios[iDataset]->Reset("M");
    divideSuccess = H1D_jetPt_rebinned_ratios[iDataset]->Divide(H1D_jetPt_rebinned[iDataset], H1D_jetPt_rebinned[0]);
  }
  // cout << "test1" << endl;

  TString* pdfName = new TString("jet_"+jetType[iJetType]+"_"+jetLevel[iJetLevel]+"_DataComp_R="+Form("%.1f", jetRadius)+"_Pt_@eta["+Form("%.1f", EtaCutLow)+","+Form("%.1f", EtaCutHigh)+"]"+jetFinderQaHistType[iJetFinderQaType]);
  TString* pdfName_ratio = new TString("jet_"+jetType[iJetType]+"_"+jetLevel[iJetLevel]+"_DataComp_R="+Form("%.1f", jetRadius)+"_Pt_@eta["+Form("%.1f", EtaCutLow)+","+Form("%.1f", EtaCutHigh)+"]"+jetFinderQaHistType[iJetFinderQaType]+"_ratio");

  // TString textContext(contextDatasetCompAndRadiusAndVarRange(jetRadius, etaRange, "eta"));
  TString textContext(contextCustomThreeFields(*texDatasetsComparisonCommonDenominator, "", "#splitline{"+contextJetRadius(jetRadius)+"}{"+contextEtaRange(etaRange)+"}", ""));

  Draw_TH1_Histograms_in_one(H1D_jetPt_rebinned, DatasetsNames, nDatasets, textContext, pdfName, texPtX, texJetPtYield_EventNorm, texCollisionDataInfo, drawnWindowAuto, "logy");
  if (divideSuccess == true) {
    Draw_TH1_Histograms_in_one(H1D_jetPt_rebinned_ratios, DatasetsNames, nDatasets, textContext, pdfName_ratio, texPtX, texRatioDatasets, texCollisionDataInfo, drawnWindowAuto, "autoratio,avoidFirst");
  }
  else {
    cout << "Divide failed in Draw_Pt_DatasetComparison" << endl;
  }
}



// to do list: 
// - add option normEntries vs normEvents (vs normCent if needed) for all functions

