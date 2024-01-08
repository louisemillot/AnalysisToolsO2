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
#include "./JetQC_settings.h"
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
TString contextDataset1D(int iDataset, float* variableRange, const char options[]);
TString contextDatasetComp(float jetRadius, float* variableRange, const char options[]);
TString contextPtRange(float* PtRange);
TString contextEtaRange(float* PtRange);
TString contextJetRadius(float jetRadius);

//////////// QC plot functions
// Radius comparison
void Draw_Pt_RadiusComparison(int iDataset, float* etaRange);
void Draw_Eta_RadiusComparison(int iDataset, float* PtRange);
void Draw_Phi_RadiusComparison(int iDataset, float* PtRange);
void Draw_NTracks_RadiusComparison_withPtRange(int iDataset, float* PtRange);
void Draw_LeadingTrackPt_vs_JetPt_RadiusComparison(int iDataset);
void Draw_JetArea_vs_JetPt_RadiusComparison(int iDataset, float* PtRange);
void Draw_JetNTracks_vs_JetPt_RadiusComparison(int iDataset, float* PtRange);
void Draw_Pt_ratio_etaNeg_etaPos_RadiusComparison(int iDataset, float* etaRange);
void Draw_JetPhi_vs_JetEta_RadiusComparison(int iDataset);
void Draw_JetEta_vs_JetPt_RadiusComparison(int iDataset, float* PtRange);
void Draw_JetPhi_vs_JetPt_RadiusComparison(int iDataset, float* PtRange);
void Draw_JetTRDratio_vs_JetEta(int iDataset);
void Draw_JetTRDcount_vs_JetEta(int iDataset);
void Draw_Pt_ratio_etaNeg_etaPos_TRDonly_vs_noTRD(int iDataset, float* etaRange);

// Dataset comparison
void Draw_Pt_DatasetComparison(float jetRadius, float* etaRange);
void Draw_Eta_DatasetComparison(float jetRadius, float* PtRange);
void Draw_Phi_DatasetComparison(float jetRadius, float* PtRange);
void Draw_Pt_ratio_etaNeg_etaPos_DatasetComparison(float jetRadius, float* etaRange);
void Draw_Area_PtIntegrated_DatasetComparison(float jetRadius, float* PtRange);


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

  int iDataset = 0;

  float etaRangeSym[2] = {-0.5, 0.5};
  float etaRangeNeg[2] = {-0.5, 0};
  float etaRangePos[2] = {0, 0.5};

  float jetRadiusForDataComp = 0.4;
  float jetR02 = 0.2;
  float jetR06 = 0.6;
  

  const int nPtMinCuts = 2;
  float jetPtMinCut;
  float jetPtMinCutArray[nPtMinCuts] = {1., 10.};

  for(int iPtMinCut = 0; iPtMinCut < nPtMinCuts; iPtMinCut++){
    jetPtMinCut = jetPtMinCutArray[iPtMinCut];

    float ptRange[2] = {jetPtMinCut, 200};
    float PtRangeZoom0[2] = {jetPtMinCut, 100};
    float PtRangeZoom020[2] = {jetPtMinCut, 20};
    float PtRangeZoom2030[2] = {20, 30};
    float PtRangeZoom3040[2] = {30, 40};
    float PtRangeZoom4050[2] = {40, 50};
    float PtRangeZoom5060[2] = {50, 60};
    float PtRangeZoom8090[2] = {80, 90};

    Draw_Eta_RadiusComparison(iDataset, ptRange);
    Draw_Phi_RadiusComparison(iDataset, ptRange);
    Draw_NTracks_RadiusComparison_withPtRange(iDataset, PtRangeZoom0);
    Draw_NTracks_RadiusComparison_withPtRange(iDataset, PtRangeZoom2030);
    Draw_NTracks_RadiusComparison_withPtRange(iDataset, PtRangeZoom3040);
    Draw_NTracks_RadiusComparison_withPtRange(iDataset, PtRangeZoom4050);
    Draw_NTracks_RadiusComparison_withPtRange(iDataset, PtRangeZoom5060);
    Draw_NTracks_RadiusComparison_withPtRange(iDataset, PtRangeZoom8090);
    // Draw_LeadingTrackPt_vs_JetPt_RadiusComparison(iDataset); leading pT not implemented yet
    Draw_JetArea_vs_JetPt_RadiusComparison(iDataset, PtRangeZoom0);
    Draw_JetNTracks_vs_JetPt_RadiusComparison(iDataset, PtRangeZoom0);
    // Draw_JetArea_vs_JetPt_RadiusComparison(iDataset, PtRangeZoom020);
    // Draw_JetNTracks_vs_JetPt_RadiusComparison(iDataset, PtRangeZoom020);
    Draw_JetPhi_vs_JetEta_RadiusComparison(iDataset);
    Draw_JetEta_vs_JetPt_RadiusComparison(iDataset, PtRangeZoom0);
    Draw_JetPhi_vs_JetPt_RadiusComparison(iDataset, PtRangeZoom0);

    Draw_Eta_DatasetComparison(jetRadiusForDataComp, ptRange);
    Draw_Eta_DatasetComparison(jetR02, ptRange);
    Draw_Eta_DatasetComparison(jetR06, ptRange);
    Draw_Phi_DatasetComparison(jetRadiusForDataComp, ptRange);
    Draw_Area_PtIntegrated_DatasetComparison(jetRadiusForDataComp, ptRange);
    Draw_Area_PtIntegrated_DatasetComparison(jetR02, ptRange);
    Draw_Area_PtIntegrated_DatasetComparison(jetR06, ptRange);
  }

  Draw_Pt_RadiusComparison(iDataset, etaRangeSym);
  Draw_Pt_RadiusComparison(iDataset, etaRangeNeg);
  Draw_Pt_RadiusComparison(iDataset, etaRangePos);
  Draw_Pt_DatasetComparison(jetRadiusForDataComp, etaRangeSym);
  Draw_Pt_ratio_etaNeg_etaPos_RadiusComparison(iDataset, etaRangeSym);
  Draw_Pt_ratio_etaNeg_etaPos_DatasetComparison(jetRadiusForDataComp, etaRangeSym);
  Draw_JetTRDcount_vs_JetEta(iDataset);
  Draw_JetTRDratio_vs_JetEta(iDataset);
  Draw_Pt_ratio_etaNeg_etaPos_TRDonly_vs_noTRD(iDataset, etaRangeSym);
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

TString contextDatasetRadiusComp(int iDataset, float* variableRange, const char options[]){
  TString texContextDatasetRadiusComp;
  if (strstr(options, "pt") != NULL) { //  || strstr(options, "ratio") != NULL not sure why I had this here
    texContextDatasetRadiusComp = "#splitline{"+*texDatasetsComparisonCommonDenominator+" "+DatasetsNames[iDataset]+"}{#splitline{2023 QC}{"+contextPtRange(variableRange)+"}}";
  }
  if (strstr(options, "eta") != NULL) { //  || strstr(options, "ratio") != NULL not sure why I had this here
    texContextDatasetRadiusComp = "#splitline{"+*texDatasetsComparisonCommonDenominator+" "+DatasetsNames[iDataset]+"}{#splitline{2023 QC}{"+contextEtaRange(variableRange)+"}}";
  }

  return texContextDatasetRadiusComp;
}

TString contextDatasetComp(float jetRadius, float* variableRange, const char options[]){
  TString texcontextDatasetComp;
  if (strstr(options, "pt") != NULL) { //  || strstr(options, "ratio") != NULL not sure why I had this here
    texcontextDatasetComp = "#splitline{"+*texDatasetsComparisonCommonDenominator+"}{#splitline{"+contextJetRadius(jetRadius)+"}{"+contextPtRange(variableRange)+"}}";
  }
  if (strstr(options, "eta") != NULL) { //  || strstr(options, "ratio") != NULL not sure why I had this here
    texcontextDatasetComp = "#splitline{"+*texDatasetsComparisonCommonDenominator+"}{#splitline{"+contextJetRadius(jetRadius)+"}{"+contextEtaRange(variableRange)+"}}";
  }

  return texcontextDatasetComp;
}

TString contextPtRange(float* PtRange){
  std::stringstream ss;
  ss << PtRange[0] << " < #it{p}_{T} < " << PtRange[1];
  TString textContext((TString)ss.str());
  // TString texDataset(Form("%.0f", PtRange[0])+" < #it{p}_{T} < "+Form("%.0f", PtRange[1]));
  return textContext;
}

TString contextEtaRange(float* EtaRange){
  std::stringstream ss;
  ss << EtaRange[0] << " < #eta < " << EtaRange[1];
  TString textContext((TString)ss.str());
  return textContext;
}

TString contextJetRadius(float jetRadius){
  std::stringstream ss;
  ss << " R = " << jetRadius;
  TString textContext((TString)ss.str());
  // TString texDataset(Form("%.0f", PtRange[0])+" < #it{p}_{T} < "+Form("%.0f", PtRange[1]));
  return textContext;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////// QC  plot functions /////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void Draw_Pt_RadiusComparison(int iDataset, float* etaRange) {

  TH3D* H3D_jetRjetPtjetEta;
  TH1D* H1D_jetPt[nRadius];
  TH1D* H1D_jetPt_rebinned[nRadius];
  
  H3D_jetRjetPtjetEta = (TH3D*)((TH3D*)file_O2Analysis_list[iDataset]->Get(analysisWorkflow+"/h3_jet_r_jet_pt_jet_eta"))->Clone("Draw_Pt_RadiusComparison"+Datasets[iDataset]+Form("%.1f", etaRange[0])+"<eta<"+Form("%.1f", etaRange[1]));
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
    ibinJetRadius = H3D_jetRjetPtjetEta->GetXaxis()->FindBin(arrayRadius[iRadius]+e_binEdge);

    H1D_jetPt[iRadius] = (TH1D*)H3D_jetRjetPtjetEta->ProjectionY("jetPt_"+RadiusLegend[iRadius]+Form("%.1f", EtaCutLow)+"<eta<"+Form("%.1f", EtaCutHigh), ibinJetRadius, ibinJetRadius, ibinEta_low, ibinEta_high);
    H1D_jetPt_rebinned[iRadius] = (TH1D*)H1D_jetPt[iRadius]->Rebin(1.,"jetPt_rebinned_"+RadiusLegend[iRadius]+Form("%.1f", EtaCutLow)+"<eta<"+Form("%.1f", EtaCutHigh));

    // NormaliseYieldToNJets(H1D_jetPt_rebinned[iRadius]);
    NormaliseYieldToNEvents(H1D_jetPt_rebinned[iRadius], GetNEventsSel8(file_O2Analysis_list[iDataset]));
  }

  TString* pdfName = new TString("jet_"+jetType[iJetType]+"_"+jetLevel[iJetLevel]+"_"+Datasets[iDataset]+"_Pt_@eta["+Form("%.1f", EtaCutLow)+","+Form("%.1f", EtaCutHigh)+"]");

  TString textContext(contextDatasetRadiusComp(iDataset, etaRange, "eta"));
  
  Draw_TH1_Histograms_in_one(H1D_jetPt_rebinned, RadiusLegend, nRadius, textContext, pdfName, texPtX, texJetNormPtYield, texCollisionDataInfo, "logy");
}

void Draw_Eta_RadiusComparison(int iDataset, float* PtRange) {

  TH3D* H3D_jetRjetPtjetEta;
  TH1D* H1D_jetEta[nRadius];
  TH1D* H1D_jetEta_rebinned[nRadius];
  
  H3D_jetRjetPtjetEta = (TH3D*)((TH3D*)file_O2Analysis_list[iDataset]->Get(analysisWorkflow+"/h3_jet_r_jet_pt_jet_eta"))->Clone("Draw_Eta_RadiusComparison"+Datasets[iDataset]+Form("%.1f", PtRange[0])+"<pt<"+Form("%.1f", PtRange[1]));
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
    ibinJetRadius = H3D_jetRjetPtjetEta->GetXaxis()->FindBin(arrayRadius[iRadius]+e_binEdge);
    H1D_jetEta[iRadius] = (TH1D*)H3D_jetRjetPtjetEta->ProjectionZ("jetEta_"+RadiusLegend[iRadius]+Form("%.1f", PtCutLow)+"<pt<"+Form("%.1f", PtCutHigh), ibinJetRadius, ibinJetRadius, ibinPt_low, ibinPt_high);
    H1D_jetEta_rebinned[iRadius] = (TH1D*)H1D_jetEta[iRadius]->Rebin(1.,"jetEta_rebinned_"+RadiusLegend[iRadius]+Form("%.1f", PtCutLow)+"<pt<"+Form("%.1f", PtCutHigh));

    // NormaliseYieldToNJets(H1D_jetEta_rebinned[iRadius]);
    NormaliseYieldToNEvents(H1D_jetEta_rebinned[iRadius], GetNEventsSel8(file_O2Analysis_list[iDataset]));
  }

  TString* pdfName = new TString("jet_"+jetType[iJetType]+"_"+jetLevel[iJetLevel]+"_"+Datasets[iDataset]+"_Eta_@pt["+Form("%03.0f", PtCutLow)+","+Form("%03.0f", PtCutHigh)+"]");

  TString textContext(contextDatasetRadiusComp(iDataset, PtRange, "pt"));

  Draw_TH1_Histograms_in_one(H1D_jetEta_rebinned, RadiusLegend, nRadius, textContext, pdfName, texEtaX, texJetNormEtaYield, texCollisionDataInfo, "");
}

void Draw_Phi_RadiusComparison(int iDataset, float* PtRange) {

  TH3D* H3D_jetRjetPtjetPhi;
  TH1D* H1D_jetPhi[nRadius];
  TH1D* H1D_jetPhi_rebinned[nRadius];
  
  H3D_jetRjetPtjetPhi = (TH3D*)((TH3D*)file_O2Analysis_list[iDataset]->Get(analysisWorkflow+"/h3_jet_r_jet_pt_jet_phi"))->Clone("Draw_Phi_RadiusComparison"+Datasets[iDataset]+Form("%.1f", PtRange[0])+"<pt<"+Form("%.1f", PtRange[1]));
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
    ibinJetRadius = H3D_jetRjetPtjetPhi->GetXaxis()->FindBin(arrayRadius[iRadius]+e_binEdge);
    H1D_jetPhi[iRadius] = (TH1D*)H3D_jetRjetPtjetPhi->ProjectionZ("jetPhi_"+RadiusLegend[iRadius]+Form("%.1f", PtCutLow)+"<pt<"+Form("%.1f", PtCutHigh), ibinJetRadius, ibinJetRadius, ibinPt_low, ibinPt_high);
    H1D_jetPhi_rebinned[iRadius] = (TH1D*)H1D_jetPhi[iRadius]->Rebin(1.,"jetPhi_rebinned_"+RadiusLegend[iRadius]+Form("%.1f", PtCutLow)+"<pt<"+Form("%.1f", PtCutHigh));

    // NormaliseYieldToNJets(H1D_jetPhi_rebinned[iRadius]);
    NormaliseYieldToNEvents(H1D_jetPhi_rebinned[iRadius], GetNEventsSel8(file_O2Analysis_list[iDataset]));
  }
 
  TString* pdfName = new TString("jet_"+jetType[iJetType]+"_"+jetLevel[iJetLevel]+"_"+Datasets[iDataset]+"_Phi_@pt["+Form("%03.0f", PtCutLow)+","+Form("%03.0f", PtCutHigh)+"]");

  TString textContext(contextDatasetRadiusComp(iDataset, PtRange, "pt"));

  Draw_TH1_Histograms_in_one(H1D_jetPhi_rebinned, RadiusLegend, nRadius, textContext, pdfName, texPhiX, texJetNormPhiYield, texCollisionDataInfo, "");
}

void Draw_NTracks_RadiusComparison_withPtRange(int iDataset, float* PtRange) {

  TH3D* H3D_jetRjetPtjetNTracks;
  TH1D* H1D_jetNTracks[nRadius];
  TH1D* H1D_jetNTracks_rebinned[nRadius];
  
  H3D_jetRjetPtjetNTracks = (TH3D*)((TH3D*)file_O2Analysis_list[iDataset]->Get(analysisWorkflow+"/h3_jet_r_jet_pt_jet_ntracks"))->Clone("Draw_NTracks_RadiusComparison_withPtRange"+Datasets[iDataset]+Form("%.1f", PtRange[0])+"<pt<"+Form("%.1f", PtRange[1]));
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
    ibinJetRadius = H3D_jetRjetPtjetNTracks->GetXaxis()->FindBin(arrayRadius[iRadius]+e_binEdge);
    H1D_jetNTracks[iRadius] = (TH1D*)H3D_jetRjetPtjetNTracks->ProjectionZ("jetNTracks_"+RadiusLegend[iRadius]+Form("%.1f", PtCutLow)+"<pt<"+Form("%.1f", PtCutHigh), ibinJetRadius, ibinJetRadius, ibinPt_low, ibinPt_high);
    H1D_jetNTracks_rebinned[iRadius] = (TH1D*)H1D_jetNTracks[iRadius]->Rebin(1.,"jetNTracks_rebinned_"+RadiusLegend[iRadius]+Form("%.1f", PtCutLow)+"<pt<"+Form("%.1f", PtCutHigh));

    // NormaliseYieldToNJets(H1D_jetNTracks_rebinned[iRadius]);
    NormaliseYieldToNEvents(H1D_jetNTracks_rebinned[iRadius], GetNEventsSel8(file_O2Analysis_list[iDataset]));
  }

  TString* pdfName = new TString("jet_"+jetType[iJetType]+"_"+jetLevel[iJetLevel]+"_"+Datasets[iDataset]+"_NTracks_@pT["+Form("%03.0f", PtCutLow)+","+Form("%03.0f", PtCutHigh)+"]");

  TString textContext(contextDatasetRadiusComp(iDataset, PtRange, "pt"));

  Draw_TH1_Histograms_in_one(H1D_jetNTracks_rebinned, RadiusLegend, nRadius, textContext, pdfName, texNTracksX, texJetNormNTracksYield, texCollisionDataInfo, "");
}

void Draw_JetArea_vs_JetPt_RadiusComparison(int iDataset, float* PtRange) {

  TH3D* H3D_jetRjetPtjetArea;
  TH2D* H2D_jetArea[nRadius];
  TH2D* H2D_jetArea_rebinned[nRadius];
  
  H3D_jetRjetPtjetArea = (TH3D*)((TH3D*)file_O2Analysis_list[iDataset]->Get(analysisWorkflow+"/h3_jet_r_jet_pt_jet_area"))->Clone("Draw_JetArea_vs_JetPt_RadiusComparison"+Datasets[iDataset]+Form("%.1f", PtRange[0])+"<pt<"+Form("%.1f", PtRange[1]));
  H3D_jetRjetPtjetArea->Sumw2();

  float PtCutLow = PtRange[0];
  float PtCutHigh = PtRange[1];
  int ibinPt_low = H3D_jetRjetPtjetArea->GetYaxis()->FindBin(PtCutLow);
  int ibinPt_high = H3D_jetRjetPtjetArea->GetYaxis()->FindBin(PtCutHigh);

  for(int iRadius = 0; iRadius < nRadius; iRadius++){
    H3D_jetRjetPtjetArea->GetXaxis()->SetRange(iRadius+1,iRadius+1);
    H3D_jetRjetPtjetArea->GetYaxis()->SetRange(ibinPt_low,ibinPt_high);
    H2D_jetArea[iRadius] = (TH2D*)H3D_jetRjetPtjetArea->Project3D(RadiusLegend[iRadius]+"_jetArea_zy");
    // H2D_jetArea_rebinned[iRadius] = (TH2D*)H2D_jetArea[iRadius]->RebinY(2.,"H1D_jetArea_rebinned_"+RadiusLegend[iRadius]);
  }

  TString* pdfName = new TString("jet_"+jetType[iJetType]+"_"+jetLevel[iJetLevel]+"_"+Datasets[iDataset]+"_JetArea-vs-Pt_@pT["+Form("%03.0f", PtCutLow)+","+Form("%03.0f", PtCutHigh)+"]");

  TString textContext(contextDatasetRadiusComp(iDataset, PtRange, "pt"));

  Draw_TH2_Histograms(H2D_jetArea, RadiusLegend, nRadius, textContext, pdfName, texPtX, texJetArea, texCollisionDataInfo, "");
  TString* pdfNamelogz = new TString(*pdfName + "_logz");
  Draw_TH2_Histograms(H2D_jetArea, RadiusLegend, nRadius, textContext, pdfNamelogz, texPtX, texJetArea, texCollisionDataInfo, "logz");
  TString* pdfNamelogy = new TString(*pdfName + "_logx");
  Draw_TH2_Histograms(H2D_jetArea, RadiusLegend, nRadius, textContext, pdfNamelogy, texPtX, texJetArea, texCollisionDataInfo, "logx");
  TString* pdfNamelogyz = new TString(*pdfName + "_logxz");
  Draw_TH2_Histograms(H2D_jetArea, RadiusLegend, nRadius, textContext, pdfNamelogyz, texPtX, texJetArea, texCollisionDataInfo, "logxlogz");
}

void Draw_JetNTracks_vs_JetPt_RadiusComparison(int iDataset, float* PtRange) {

  TH3D* H3D_jetRjetPtjetNTracks;
  TH2D* H2D_jetNTracks[nRadius];
  TH2D* H2D_jetNTracks_rebinned[nRadius];
  
  H3D_jetRjetPtjetNTracks = (TH3D*)((TH3D*)file_O2Analysis_list[iDataset]->Get(analysisWorkflow+"/h3_jet_r_jet_pt_jet_ntracks"))->Clone("Draw_JetNTracks_vs_JetPt_RadiusComparison"+Datasets[iDataset]+Form("%.1f", PtRange[0])+"<pt<"+Form("%.1f", PtRange[1]));
  H3D_jetRjetPtjetNTracks->Sumw2();

  float PtCutLow = PtRange[0];
  float PtCutHigh = PtRange[1];
  int ibinPt_low = H3D_jetRjetPtjetNTracks->GetYaxis()->FindBin(PtCutLow);
  int ibinPt_high = H3D_jetRjetPtjetNTracks->GetYaxis()->FindBin(PtCutHigh);

  for(int iRadius = 0; iRadius < nRadius; iRadius++){
    H3D_jetRjetPtjetNTracks->GetXaxis()->SetRange(iRadius+1,iRadius+1);
    H3D_jetRjetPtjetNTracks->GetYaxis()->SetRange(ibinPt_low,ibinPt_high);
    H2D_jetNTracks[iRadius] = (TH2D*)H3D_jetRjetPtjetNTracks->Project3D(RadiusLegend[iRadius]+"_jetNTracks_zy");
    // H2D_jetNTracks_rebinned[iRadius] = (TH2D*)H2D_jetNTracks[iRadius]->Rebin(1.,"jetNTracks_rebinned_"+RadiusLegend[iRadius]);
  }

  TString* pdfName = new TString("jet_"+jetType[iJetType]+"_"+jetLevel[iJetLevel]+"_"+Datasets[iDataset]+"_JetNTracks-vs-Pt_@pT["+Form("%03.0f", PtCutLow)+","+Form("%03.0f", PtCutHigh)+"]");

  TString textContext(contextDatasetRadiusComp(iDataset, PtRange, "pt"));

  Draw_TH2_Histograms(H2D_jetNTracks, RadiusLegend, nRadius, textContext, pdfName, texPtX, texJetNTracks, texCollisionDataInfo, "");
  TString* pdfNamelogz = new TString(*pdfName + "_logz");
  Draw_TH2_Histograms(H2D_jetNTracks, RadiusLegend, nRadius, textContext, pdfNamelogz, texPtX, texJetNTracks, texCollisionDataInfo, "logz");
}

void Draw_JetEta_vs_JetPt_RadiusComparison(int iDataset, float* PtRange) {

  TH3D* H3D_jetRjetPtjetEta;
  TH2D* H2D_jetEta[nRadius];
  TH2D* H2D_jetEta_rebinned[nRadius];
  
  H3D_jetRjetPtjetEta = (TH3D*)((TH3D*)file_O2Analysis_list[iDataset]->Get(analysisWorkflow+"/h3_jet_r_jet_pt_jet_eta"))->Clone("Draw_JetEta_vs_JetPt_RadiusComparison"+Datasets[iDataset]+Form("%.1f", PtRange[0])+"<pt<"+Form("%.1f", PtRange[1]));
  H3D_jetRjetPtjetEta->Sumw2();

  float PtCutLow = PtRange[0];
  float PtCutHigh = PtRange[1];
  int ibinPt_low = H3D_jetRjetPtjetEta->GetYaxis()->FindBin(PtCutLow);
  int ibinPt_high = H3D_jetRjetPtjetEta->GetYaxis()->FindBin(PtCutHigh);

  for(int iRadius = 0; iRadius < nRadius; iRadius++){
    H3D_jetRjetPtjetEta->GetXaxis()->SetRange(iRadius+1,iRadius+1);
    H3D_jetRjetPtjetEta->GetYaxis()->SetRange(ibinPt_low,ibinPt_high);
    H2D_jetEta[iRadius] = (TH2D*)H3D_jetRjetPtjetEta->Project3D(RadiusLegend[iRadius]+"_jetEta_zy");
    // H2D_jetEta_rebinned[iRadius] = (TH2D*)H2D_jetEta[iRadius]->RebinY(2.,"H1D_jetEta_rebinned_"+RadiusLegend[iRadius]);
  }

  TString* pdfName = new TString("jet_"+jetType[iJetType]+"_"+jetLevel[iJetLevel]+"_"+Datasets[iDataset]+"_JetEta-vs-Pt_@pT["+Form("%03.0f", PtCutLow)+","+Form("%03.0f", PtCutHigh)+"]");

  TString textContext(contextDatasetRadiusComp(iDataset, PtRange, "pt"));

  Draw_TH2_Histograms(H2D_jetEta, RadiusLegend, nRadius, textContext, pdfName, texPtX, texEtaX, texCollisionDataInfo, "");
  TString* pdfNamelogz = new TString(*pdfName + "_logz");
  Draw_TH2_Histograms(H2D_jetEta, RadiusLegend, nRadius, textContext, pdfNamelogz, texPtX, texEtaX, texCollisionDataInfo, "logz");
  TString* pdfNamelogy = new TString(*pdfName + "_logx");
  Draw_TH2_Histograms(H2D_jetEta, RadiusLegend, nRadius, textContext, pdfNamelogy, texPtX, texEtaX, texCollisionDataInfo, "logx");
  TString* pdfNamelogyz = new TString(*pdfName + "_logxz");
  Draw_TH2_Histograms(H2D_jetEta, RadiusLegend, nRadius, textContext, pdfNamelogyz, texPtX, texEtaX, texCollisionDataInfo, "logxlogz");
}

void Draw_JetPhi_vs_JetPt_RadiusComparison(int iDataset, float* PtRange) {

  TH3D* H3D_jetRjetPtjetPhi;
  TH2D* H2D_jetPhi[nRadius];
  TH2D* H2D_jetPhi_rebinned[nRadius];
  
  H3D_jetRjetPtjetPhi = (TH3D*)((TH3D*)file_O2Analysis_list[iDataset]->Get(analysisWorkflow+"/h3_jet_r_jet_pt_jet_phi"))->Clone("Draw_JetPhi_vs_JetPt_RadiusComparison"+Datasets[iDataset]+Form("%.1f", PtRange[0])+"<pt<"+Form("%.1f", PtRange[1]));
  H3D_jetRjetPtjetPhi->Sumw2();

  float PtCutLow = PtRange[0];
  float PtCutHigh = PtRange[1];
  int ibinPt_low = H3D_jetRjetPtjetPhi->GetYaxis()->FindBin(PtCutLow);
  int ibinPt_high = H3D_jetRjetPtjetPhi->GetYaxis()->FindBin(PtCutHigh);

  for(int iRadius = 0; iRadius < nRadius; iRadius++){
    H3D_jetRjetPtjetPhi->GetXaxis()->SetRange(iRadius+1,iRadius+1);
    H3D_jetRjetPtjetPhi->GetYaxis()->SetRange(ibinPt_low,ibinPt_high);
    H2D_jetPhi[iRadius] = (TH2D*)H3D_jetRjetPtjetPhi->Project3D(RadiusLegend[iRadius]+Datasets[iDataset]+Form("%.1f", PtRange[0])+"<pt<"+Form("%.1f", PtRange[1])+"_jetPhi_jetPt_zy");
    // H2D_jetPhi_rebinned[iRadius] = (TH2D*)H2D_jetPhi[iRadius]->RebinY(2.,"H1D_jetPhi_rebinned_"+RadiusLegend[iRadius]);
  }

  TString* pdfName = new TString("jet_"+jetType[iJetType]+"_"+jetLevel[iJetLevel]+"_"+Datasets[iDataset]+"_JetPhi-vs-Pt_@pT["+Form("%03.0f", PtCutLow)+","+Form("%03.0f", PtCutHigh)+"]");

  TString textContext(contextDatasetRadiusComp(iDataset, PtRange, "pt"));

  Draw_TH2_Histograms(H2D_jetPhi, RadiusLegend, nRadius, textContext, pdfName, texPtX, texPhiX, texCollisionDataInfo, "");
  TString* pdfNamelogz = new TString(*pdfName + "_logz");
  Draw_TH2_Histograms(H2D_jetPhi, RadiusLegend, nRadius, textContext, pdfNamelogz, texPtX, texPhiX, texCollisionDataInfo, "logz");
  TString* pdfNamelogy = new TString(*pdfName + "_logx");
  Draw_TH2_Histograms(H2D_jetPhi, RadiusLegend, nRadius, textContext, pdfNamelogy, texPtX, texPhiX, texCollisionDataInfo, "logx");
  TString* pdfNamelogyz = new TString(*pdfName + "_logxz");
  Draw_TH2_Histograms(H2D_jetPhi, RadiusLegend, nRadius, textContext, pdfNamelogyz, texPtX, texPhiX, texCollisionDataInfo, "logxlogz");
}

void Draw_JetPhi_vs_JetEta_RadiusComparison(int iDataset) {

  TH3D* H3D_jetRjetPhijetEta;
  TH2D* H2D_jetPhijetEta[nRadius];
  TH2D* H2D_jetPhijetEta_rebinned[nRadius];
  
  H3D_jetRjetPhijetEta = (TH3D*)((TH3D*)file_O2Analysis_list[iDataset]->Get(analysisWorkflow+"/h3_jet_r_jet_eta_jet_phi"))->Clone("Draw_JetPhi_vs_JetEta_RadiusComparison"+Datasets[iDataset]);
  H3D_jetRjetPhijetEta->Sumw2();

  for(int iRadius = 0; iRadius < nRadius; iRadius++){
    H2D_jetPhijetEta[iRadius] = (TH2D*)H3D_jetRjetPhijetEta->Project3D(RadiusLegend[iRadius]+Datasets[iDataset]+"_jetPhi_jetEta_zy");
    // H2D_jetPhijetEta_rebinned[iRadius] = (TH2D*)H2D_jetArea[iRadius]->RebinY(2.,"H1D_jetPhi_jetEta_rebinned_"+RadiusLegend[iRadius]);
  }

  TString* pdfName = new TString("jet_"+jetType[iJetType]+"_"+jetLevel[iJetLevel]+"_"+Datasets[iDataset]+"_JetPhi-vs-Eta");

  TString textContext("#splitline{"+*texDatasetsComparisonCommonDenominator+" "+DatasetsNames[iDataset]+"}{2023 QC}");

  Draw_TH2_Histograms(H2D_jetPhijetEta, RadiusLegend, nRadius, textContext, pdfName, texEtaX, texPhiX, texCollisionDataInfo, "");
  // TString* pdfNamelogz = new TString(*pdfName + "_logz");
  // Draw_TH2_Histograms(H2D_jetPhijetEta, RadiusLegend, nRadius, textContext, pdfNamelogz, texPtX, texJetArea, texCollisionDataInfo, "logz");
  // TString* pdfNamelogy = new TString(*pdfName + "_logx");
  // Draw_TH2_Histograms(H2D_jetPhijetEta, RadiusLegend, nRadius, textContext, pdfNamelogy, texPtX, texJetArea, texCollisionDataInfo, "logx");
  // TString* pdfNamelogyz = new TString(*pdfName + "_logxz");
  // Draw_TH2_Histograms(H2D_jetPhijetEta, RadiusLegend, nRadius, textContext, pdfNamelogyz, texPtX, texJetArea, texCollisionDataInfo, "logxlogz");
  
}

void Draw_Pt_DatasetComparison(float jetRadius, float* etaRange) {

  TH3D* H3D_jetRjetPtjetEta[nDatasets];
  TH1D* H1D_jetPt[nDatasets];
  TH1D* H1D_jetPt_rebinned[nDatasets];
  
  TH1D* H1D_jetPt_rebinned_ratios[nDatasets];

  float EtaCutLow = etaRange[0];
  float EtaCutHigh = etaRange[1];
  int ibinJetRadius = 0;

  bool divideSuccess = false;

  for(int iDataset = 0; iDataset < nDatasets; iDataset++){

    H3D_jetRjetPtjetEta[iDataset] = (TH3D*)((TH3D*)file_O2Analysis_list[iDataset]->Get(analysisWorkflow+"/h3_jet_r_jet_pt_jet_eta"))->Clone("Draw_Pt_DatasetComparison"+Datasets[iDataset]+"Radius"+Form("%.1f",jetRadius)+jetRadius+Form("%.1f", etaRange[0])+"<eta<"+Form("%.1f", etaRange[1]));

    int ibinEta_low = H3D_jetRjetPtjetEta[iDataset]->GetZaxis()->FindBin(EtaCutLow);
    int ibinEta_high = H3D_jetRjetPtjetEta[iDataset]->GetZaxis()->FindBin(EtaCutHigh);
    if (ibinEta_low == 0) 
      cout << "WARNING: Pt_DatasetComparison is counting the underflow with the chosen etaRange" << endl;
    if (ibinEta_high == H3D_jetRjetPtjetEta[iDataset]->GetZaxis()->GetNbins()+1) 
      cout << "WARNING: Pt_DatasetComparison is counting the overflow with the chosen etaRange" << endl;
    ibinJetRadius = H3D_jetRjetPtjetEta[iDataset]->GetXaxis()->FindBin(jetRadius+e_binEdge);
    
    H1D_jetPt[iDataset] = (TH1D*)H3D_jetRjetPtjetEta[iDataset]->ProjectionY("jetPt_"+Datasets[iDataset]+"Radius"+Form("%.1f",jetRadius)+Form("%.1f", EtaCutLow)+"<eta<"+Form("%.1f", EtaCutHigh), ibinJetRadius, ibinJetRadius, ibinEta_low, ibinEta_high);
    H1D_jetPt_rebinned[iDataset] = (TH1D*)H1D_jetPt[iDataset]->Rebin(1.,"jetPt_rebinned_"+Datasets[iDataset]+"Radius"+Form("%.1f",jetRadius)+Form("%.1f", EtaCutLow)+"<eta<"+Form("%.1f", EtaCutHigh));

    // NormaliseYieldToNJets(H1D_jetPt_rebinned[iDataset]);
    NormaliseYieldToNEvents(H1D_jetPt_rebinned[iDataset], GetNEventsSel8(file_O2Analysis_list[iDataset]));

    H1D_jetPt_rebinned_ratios[iDataset] = (TH1D*)H1D_jetPt_rebinned[iDataset]->Clone("jetPt_rebinned_ratios"+Datasets[iDataset]+"Radius"+Form("%.1f",jetRadius)+Form("%.1f", EtaCutLow)+"<eta<"+Form("%.1f", EtaCutHigh));
    H1D_jetPt_rebinned_ratios[iDataset]->Reset("M");
    divideSuccess = H1D_jetPt_rebinned_ratios[iDataset]->Divide(H1D_jetPt_rebinned[iDataset], H1D_jetPt_rebinned[0]);
  }

  TString* pdfName = new TString("jet_"+jetType[iJetType]+"_"+jetLevel[iJetLevel]+"_DataComp_R="+Form("%.1f", jetRadius)+"_Pt_@eta["+Form("%.1f", EtaCutLow)+","+Form("%.1f", EtaCutHigh)+"]");
  TString* pdfName_ratio = new TString("jet_"+jetType[iJetType]+"_"+jetLevel[iJetLevel]+"_DataComp_R="+Form("%.1f", jetRadius)+"_Pt_@eta["+Form("%.1f", EtaCutLow)+","+Form("%.1f", EtaCutHigh)+"]_ratio");

  TString textContext(contextDatasetComp(jetRadius, etaRange, "eta"));

  Draw_TH1_Histograms_in_one(H1D_jetPt_rebinned, DatasetsNames, nDatasets, textContext, pdfName, texPtX, texJetNormPtYield, texCollisionDataInfo, "logy");
  if (divideSuccess == true) {
    Draw_TH1_Histograms_in_one(H1D_jetPt_rebinned_ratios, DatasetsNames, nDatasets, textContext, pdfName_ratio, texPtX, texRatioDatasets, texCollisionDataInfo, "autoratio,avoidFirst");
  }
  else {
    cout << "Divide failed in Draw_Pt_DatasetComparison" << endl;
  }
}

void Draw_Eta_DatasetComparison(float jetRadius, float* PtRange) {

  TH3D* H3D_jetRjetPtjetEta[nDatasets];
  TH1D* H1D_jetEta[nDatasets];
  TH1D* H1D_jetEta_rebinned[nDatasets];
  
  TH1D* H1D_jetEta_rebinned_ratios[nDatasets];

  float PtCutLow = PtRange[0];
  float PtCutHigh = PtRange[1];
  int ibinJetRadius = 0;

  bool divideSuccess = false;

  for(int iDataset = 0; iDataset < nDatasets; iDataset++){

    H3D_jetRjetPtjetEta[iDataset] = (TH3D*)((TH3D*)file_O2Analysis_list[iDataset]->Get(analysisWorkflow+"/h3_jet_r_jet_pt_jet_eta"))->Clone("Draw_Pt_DatasetComparison"+Datasets[iDataset]+"Radius"+Form("%.1f",jetRadius)+Form("%.1f", PtRange[0])+"<pt<"+Form("%.1f", PtRange[1]));

    int ibinPt_low = H3D_jetRjetPtjetEta[iDataset]->GetYaxis()->FindBin(PtCutLow);
    int ibinPt_high = H3D_jetRjetPtjetEta[iDataset]->GetYaxis()->FindBin(PtCutHigh);
    if (ibinPt_low == 0) 
      cout << "WARNING: Eta_DatasetComparison is counting the underflow with the chosen PtRange" << endl;
    if (ibinPt_high == H3D_jetRjetPtjetEta[iDataset]->GetYaxis()->GetNbins()+1) 
      cout << "WARNING: Eta_DatasetComparison is counting the overflow with the chosen PtRange" << endl;
    ibinJetRadius = H3D_jetRjetPtjetEta[iDataset]->GetXaxis()->FindBin(jetRadius+e_binEdge);

    H1D_jetEta[iDataset] = (TH1D*)H3D_jetRjetPtjetEta[iDataset]->ProjectionZ("jetEta_"+Datasets[iDataset]+"Radius"+Form("%.1f",jetRadius)+Form("%.1f", PtCutLow)+"<pt<"+Form("%.1f", PtCutHigh), ibinJetRadius, ibinJetRadius, ibinPt_low, ibinPt_high);
    H1D_jetEta_rebinned[iDataset] = (TH1D*)H1D_jetEta[iDataset]->Rebin(1.,"jetEta_rebinned_"+Datasets[iDataset]+"Radius"+Form("%.1f",jetRadius)+Form("%.1f", PtCutLow)+"<pt<"+Form("%.1f", PtCutHigh));

    // NormaliseYieldToNJets(H1D_jetEta_rebinned[iDataset]);
    NormaliseYieldToNEvents(H1D_jetEta_rebinned[iDataset], GetNEventsSel8(file_O2Analysis_list[iDataset]));

    H1D_jetEta_rebinned_ratios[iDataset] = (TH1D*)H1D_jetEta_rebinned[iDataset]->Clone("jetEta_rebinned_ratios"+Datasets[iDataset]+"Radius"+Form("%.1f",jetRadius)+Form("%.1f", PtCutLow)+"<pt<"+Form("%.1f", PtCutHigh));
    H1D_jetEta_rebinned_ratios[iDataset]->Reset("M");
    divideSuccess = H1D_jetEta_rebinned_ratios[iDataset]->Divide(H1D_jetEta_rebinned[iDataset], H1D_jetEta_rebinned[0]);
  }

  TString* pdfName = new TString("jet_"+jetType[iJetType]+"_"+jetLevel[iJetLevel]+"_DataComp_R="+Form("%.1f", jetRadius)+"_Eta_@pT["+Form("%03.0f", PtCutLow)+","+Form("%03.0f", PtCutHigh)+"]");
  TString* pdfName_ratio = new TString("jet_"+jetType[iJetType]+"_"+jetLevel[iJetLevel]+"_DataComp_R="+Form("%.1f", jetRadius)+"_Eta_@pT["+Form("%03.0f", PtCutLow)+","+Form("%03.0f", PtCutHigh)+"]_ratio");

  TString textContext(contextDatasetComp(jetRadius, PtRange, "pt"));

  Draw_TH1_Histograms_in_one(H1D_jetEta_rebinned, DatasetsNames, nDatasets, textContext, pdfName, texEtaX, texJetNormEtaYield, texCollisionDataInfo, "");

  if (divideSuccess == true) {
    Draw_TH1_Histograms_in_one(H1D_jetEta_rebinned_ratios, DatasetsNames, nDatasets, textContext, pdfName_ratio, texEtaX, texRatioDatasets, texCollisionDataInfo, "autoratio,avoidFirst");
  }
  else {
    cout << "Divide failed in Draw_Eta_DatasetComparison" << endl;
  }
}

void Draw_Phi_DatasetComparison(float jetRadius, float* PtRange) {

  TH3D* H3D_jetRjetPtjetPhi[nDatasets];
  TH1D* H1D_jetPhi[nDatasets];
  TH1D* H1D_jetPhi_rebinned[nDatasets];
  
  TH1D* H1D_jetPhi_rebinned_ratios[nDatasets];

  float PtCutLow = PtRange[0];
  float PtCutHigh = PtRange[1];
  int ibinJetRadius = 0;

  bool divideSuccess = false;

  for(int iDataset = 0; iDataset < nDatasets; iDataset++){

    H3D_jetRjetPtjetPhi[iDataset] = (TH3D*)((TH3D*)file_O2Analysis_list[iDataset]->Get(analysisWorkflow+"/h3_jet_r_jet_pt_jet_phi"))->Clone("Draw_Phi_DatasetComparison"+Datasets[iDataset]+"Radius"+Form("%.1f",jetRadius)+Form("%.1f", PtRange[0])+"<pt<"+Form("%.1f", PtRange[1]));

    int ibinPt_low = H3D_jetRjetPtjetPhi[iDataset]->GetYaxis()->FindBin(PtCutLow);
    int ibinPt_high = H3D_jetRjetPtjetPhi[iDataset]->GetYaxis()->FindBin(PtCutHigh);
    if (ibinPt_low == 0) 
      cout << "WARNING: Phi_DatasetComparison is counting the underflow with the chosen PtRange" << endl;
    if (ibinPt_high == H3D_jetRjetPtjetPhi[iDataset]->GetYaxis()->GetNbins()+1) 
      cout << "WARNING: Phi_DatasetComparison is counting the overflow with the chosen PtRange" << endl;
    ibinJetRadius = H3D_jetRjetPtjetPhi[iDataset]->GetXaxis()->FindBin(jetRadius+e_binEdge);

    H1D_jetPhi[iDataset] = (TH1D*)H3D_jetRjetPtjetPhi[iDataset]->ProjectionZ("jetPhi_"+Datasets[iDataset]+"Radius"+Form("%.1f",jetRadius)+Form("%.1f", PtCutLow)+"<pt<"+Form("%.1f", PtCutHigh), ibinJetRadius,ibinJetRadius, ibinPt_low, ibinPt_high);
    H1D_jetPhi_rebinned[iDataset] = (TH1D*)H1D_jetPhi[iDataset]->Rebin(1.,"jetPhi_rebinned_"+Datasets[iDataset]+"Radius"+Form("%.1f",jetRadius)+Form("%.1f", PtCutLow)+"<pt<"+Form("%.1f", PtCutHigh));

    // NormaliseYieldToNJets(H1D_jetPhi_rebinned[iDataset]);
    NormaliseYieldToNEvents(H1D_jetPhi_rebinned[iDataset], GetNEventsSel8(file_O2Analysis_list[iDataset]));

    H1D_jetPhi_rebinned_ratios[iDataset] = (TH1D*)H1D_jetPhi_rebinned[iDataset]->Clone("jetPhi_rebinned_ratios"+Datasets[iDataset]+"Radius"+Form("%.1f",jetRadius)+Form("%.1f", PtCutLow)+"<pt<"+Form("%.1f", PtCutHigh));
    H1D_jetPhi_rebinned_ratios[iDataset]->Reset("M");
    divideSuccess = H1D_jetPhi_rebinned_ratios[iDataset]->Divide(H1D_jetPhi_rebinned[iDataset], H1D_jetPhi_rebinned[0]);
  }

  TString* pdfName = new TString("jet_"+jetType[iJetType]+"_"+jetLevel[iJetLevel]+"_DataComp_R="+Form("%.1f", jetRadius)+"_Phi_@pT["+Form("%03.0f", PtCutLow)+","+Form("%03.0f", PtCutHigh)+"]");
  TString* pdfName_ratio = new TString("jet_"+jetType[iJetType]+"_"+jetLevel[iJetLevel]+"_DataComp_R="+Form("%.1f", jetRadius)+"_Phi_@pT["+Form("%03.0f", PtCutLow)+","+Form("%03.0f", PtCutHigh)+"]_ratio");

  TString textContext(contextDatasetComp(jetRadius, PtRange, "pt"));

  Draw_TH1_Histograms_in_one(H1D_jetPhi_rebinned, DatasetsNames, nDatasets, textContext, pdfName, texPhiX, texJetNormPhiYield, texCollisionDataInfo, "");

  if (divideSuccess == true) {
    Draw_TH1_Histograms_in_one(H1D_jetPhi_rebinned_ratios, DatasetsNames, nDatasets, textContext, pdfName_ratio, texPhiX, texRatioDatasets, texCollisionDataInfo, "autoratio,avoidFirst");
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

  H3D_jetRjetPtjetEta = (TH3D*)((TH3D*)file_O2Analysis_list[iDataset]->Get(analysisWorkflow+"/h3_jet_r_jet_pt_jet_eta"))->Clone("Draw_Pt_ratio_etaNeg_etaPos_RadiusComparison"+Datasets[iDataset]+Form("%.1f", etaRange[0])+"<eta<"+Form("%.1f", etaRange[1]));
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
    ibinJetRadius = H3D_jetRjetPtjetEta->GetXaxis()->FindBin(arrayRadius[iRadius]+e_binEdge);
    H1D_jetPt_left[iRadius] = (TH1D*)H3D_jetRjetPtjetEta->ProjectionY("jetPt_left_"+Datasets[iDataset]+RadiusLegend[iRadius], ibinJetRadius, ibinJetRadius, ibinEta_low, ibinEta_zero-1);
    H1D_jetPt_right[iRadius] = (TH1D*)H3D_jetRjetPtjetEta->ProjectionY("jetPt_right_"+Datasets[iDataset]+RadiusLegend[iRadius], ibinJetRadius, ibinJetRadius, ibinEta_zero, ibinEta_high);
    H1D_jetPt_left_rebinned[iRadius] = (TH1D*)H1D_jetPt_left[iRadius]->Rebin(1.,"jetPt_left_rebinned_"+Datasets[iDataset]+RadiusLegend[iRadius]+Form("%.1f",arrayRadius[iRadius])+Form("%.1f", etaRange[0])+"<eta<"+Form("%.1f", etaRange[1]));
    H1D_jetPt_right_rebinned[iRadius] = (TH1D*)H1D_jetPt_right[iRadius]->Rebin(1.,"jetPt_right_rebinned_"+Datasets[iDataset]+RadiusLegend[iRadius]+Form("%.1f",arrayRadius[iRadius])+Form("%.1f", etaRange[0])+"<eta<"+Form("%.1f", etaRange[1]));

    // NormaliseYieldToNJets(H1D_jetPt_left_rebinned[iRadius]);
    // NormaliseYieldToNJets(H1D_jetPt_right_rebinned[iRadius]);
    NormaliseYieldToNEvents(H1D_jetPt_left_rebinned[iRadius], GetNEventsSel8(file_O2Analysis_list[iDataset]));
    NormaliseYieldToNEvents(H1D_jetPt_right_rebinned[iRadius], GetNEventsSel8(file_O2Analysis_list[iDataset]));

    H1D_jetPt_rebinned_ratios[iRadius] = (TH1D*)H1D_jetPt_left_rebinned[iRadius]->Clone("jetPt_rebinned_ratios"+Datasets[iRadius]+RadiusLegend[iRadius]+Form("%.1f",arrayRadius[iRadius])+Form("%.1f", etaRange[0])+"<eta<"+Form("%.1f", etaRange[1]));
    H1D_jetPt_rebinned_ratios[iRadius]->Reset("M");
    divideSuccess = H1D_jetPt_rebinned_ratios[iRadius]->Divide(H1D_jetPt_right_rebinned[iRadius], H1D_jetPt_left_rebinned[iRadius]);
  }

  TString* pdfName = new TString("jet_"+jetType[iJetType]+"_"+jetLevel[iJetLevel]+"_"+Datasets[iDataset]+"_pT_etaRightLeftRatio_@eta["+Form("%.1f", EtaCutLow)+","+Form("%.1f", EtaCutHigh)+"]");

  TString textContext(contextDatasetRadiusComp(iDataset, etaRange, "eta"));

  if (divideSuccess == true) {
    Draw_TH1_Histograms_in_one(H1D_jetPt_rebinned_ratios, RadiusLegend, nRadius, textContext, pdfName, texPtX, texJetRatioEtaComparison, texCollisionDataInfo, "standardratio");
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
    H3D_jetRjetPtjetEta[iDataset] = (TH3D*)((TH3D*)file_O2Analysis_list[iDataset]->Get(analysisWorkflow+"/h3_jet_r_jet_pt_jet_eta"))->Clone("Draw_Pt_ratio_etaNeg_etaPos_DatasetComparison"+Datasets[iDataset]+"Radius"+Form("%.1f",jetRadius)+Form("%.1f", etaRange[0])+"<eta<"+Form("%.1f", etaRange[1]));
    H3D_jetRjetPtjetEta[iDataset]->Sumw2();

    int ibinEta_low = H3D_jetRjetPtjetEta[iDataset]->GetZaxis()->FindBin(EtaCutLow);
    int ibinEta_zero = H3D_jetRjetPtjetEta[iDataset]->GetZaxis()->FindBin(0.);
    int ibinEta_high = H3D_jetRjetPtjetEta[iDataset]->GetZaxis()->FindBin(EtaCutHigh);
    if (ibinEta_low == 0) 
      cout << "WARNING: Pt_ratio_etaNeg_etaPos_DatasetComparison is counting the underflow with the chosen etaRange" << endl;
    if (ibinEta_high == H3D_jetRjetPtjetEta[iDataset]->GetZaxis()->GetNbins()+1) 
      cout << "WARNING: Pt_ratio_etaNeg_etaPos_DatasetComparison is counting the overflow with the chosen etaRange" << endl;

    ibinJetRadius = H3D_jetRjetPtjetEta[iDataset]->GetXaxis()->FindBin(jetRadius+e_binEdge);

    H1D_jetPt_left[iDataset] = (TH1D*)H3D_jetRjetPtjetEta[iDataset]->ProjectionY("jetPt_left_"+Datasets[iDataset]+"Radius"+Form("%.1f",jetRadius)+Form("%.1f", etaRange[0])+"<eta<"+Form("%.1f", etaRange[1]),ibinJetRadius,ibinJetRadius, ibinEta_low, ibinEta_zero-1);
    H1D_jetPt_right[iDataset] = (TH1D*)H3D_jetRjetPtjetEta[iDataset]->ProjectionY("jetPt_right_"+Datasets[iDataset]+"Radius"+Form("%.1f",jetRadius)+Form("%.1f", etaRange[0])+"<eta<"+Form("%.1f", etaRange[1]),ibinJetRadius,ibinJetRadius, ibinEta_zero, ibinEta_high);
    H1D_jetPt_left_rebinned[iDataset] = (TH1D*)H1D_jetPt_left[iDataset]->Rebin(1.,"jetPt_left_rebinned_"+Datasets[iDataset]+"Radius"+Form("%.1f",jetRadius)+Form("%.1f", etaRange[0])+"<eta<"+Form("%.1f", etaRange[1]));
    H1D_jetPt_right_rebinned[iDataset] = (TH1D*)H1D_jetPt_right[iDataset]->Rebin(1.,"jetPt_right_rebinned_"+Datasets[iDataset]+"Radius"+Form("%.1f",jetRadius)+Form("%.1f", etaRange[0])+"<eta<"+Form("%.1f", etaRange[1]));
    
    // NormaliseYieldToNJets(H1D_jetPt_left_rebinned[iDataset]);
    // NormaliseYieldToNJets(H1D_jetPt_right_rebinned[iDataset]);
    NormaliseYieldToNEvents(H1D_jetPt_left_rebinned[iDataset], GetNEventsSel8(file_O2Analysis_list[iDataset]));
    NormaliseYieldToNEvents(H1D_jetPt_right_rebinned[iDataset],  GetNEventsSel8(file_O2Analysis_list[iDataset]));

    H1D_jetPt_rebinned_ratios[iDataset] = (TH1D*)H1D_jetPt_left_rebinned[iDataset]->Clone("jetPt_rebinned_ratios"+Datasets[iDataset]+"Radius"+Form("%.1f",jetRadius)+Form("%.1f", etaRange[0])+"<eta<"+Form("%.1f", etaRange[1]));
    H1D_jetPt_rebinned_ratios[iDataset]->Reset("M");
    divideSuccess = H1D_jetPt_rebinned_ratios[iDataset]->Divide(H1D_jetPt_right_rebinned[iDataset], H1D_jetPt_left_rebinned[iDataset]);
  }
 
  TString* pdfName = new TString("jet_"+jetType[iJetType]+"_"+jetLevel[iJetLevel]+"_DataComp_R="+Form("%.1f", jetRadius)+"_pT_etaRightLeftRatio_@eta["+Form("%.1f", EtaCutLow)+","+Form("%.1f", EtaCutHigh)+"]");

  TString textContext(contextDatasetComp(jetRadius, etaRange, "eta"));

  if (divideSuccess == true) {
    Draw_TH1_Histograms_in_one(H1D_jetPt_rebinned_ratios, DatasetsNames, nDatasets, textContext, pdfName, texPtX, texJetRatioEtaComparison, texCollisionDataInfo, "standardratio");
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

    H3D_jetRjetPtjetArea[iDataset] = (TH3D*)((TH3D*)file_O2Analysis_list[iDataset]->Get(analysisWorkflow+"/h3_jet_r_jet_pt_jet_area"))->Clone("Draw_Area_PtIntegrated_DatasetComparison"+Datasets[iDataset]+"Radius"+Form("%.1f",jetRadius)+Form("%.1f", PtRange[0])+"<pt<"+Form("%.1f", PtRange[1]));

    int ibinPt_low = H3D_jetRjetPtjetArea[iDataset]->GetYaxis()->FindBin(PtCutLow);
    int ibinPt_high = H3D_jetRjetPtjetArea[iDataset]->GetYaxis()->FindBin(PtCutHigh);
    if (ibinPt_low == 0) 
      cout << "WARNING: Area_DatasetComparison is counting the underflow with the chosen PtRange" << endl;
    if (ibinPt_high == H3D_jetRjetPtjetArea[iDataset]->GetYaxis()->GetNbins()+1) 
      cout << "WARNING: Area_DatasetComparison is counting the overflow with the chosen PtRange" << endl;
    ibinJetRadius = H3D_jetRjetPtjetArea[iDataset]->GetXaxis()->FindBin(jetRadius+e_binEdge);

    H1D_jetArea[iDataset] = (TH1D*)H3D_jetRjetPtjetArea[iDataset]->ProjectionZ("jetArea_"+Datasets[iDataset]+"Radius"+Form("%.1f",jetRadius)+Form("%.1f", PtCutLow)+"<pt<"+Form("%.1f", PtCutHigh), ibinJetRadius,ibinJetRadius, ibinPt_low, ibinPt_high);
    H1D_jetArea[iDataset]->Sumw2();
    H1D_jetArea_rebinned[iDataset] = (TH1D*)H1D_jetArea[iDataset]->Rebin(2.,"jetArea_rebinned_"+Datasets[iDataset]+"Radius"+Form("%.1f",jetRadius)+Form("%.1f", PtCutLow)+"<pt<"+Form("%.1f", PtCutHigh));

    // NormaliseYieldToNJets(H1D_jetArea_rebinned[iDataset]);
    NormaliseYieldToNEvents(H1D_jetArea_rebinned[iDataset], GetNEventsSel8(file_O2Analysis_list[iDataset]));

    H1D_jetArea_rebinned_ratios[iDataset] = (TH1D*)H1D_jetArea_rebinned[iDataset]->Clone("jetArea_rebinned_ratios"+Datasets[iDataset]+"Radius"+Form("%.1f",jetRadius)+Form("%.1f", PtCutLow)+"<pt<"+Form("%.1f", PtCutHigh));
    H1D_jetArea_rebinned_ratios[iDataset]->Reset("M");
    divideSuccess = H1D_jetArea_rebinned_ratios[iDataset]->Divide(H1D_jetArea_rebinned[iDataset], H1D_jetArea_rebinned[0]);
  }

  TString* pdfName = new TString("jet_"+jetType[iJetType]+"_"+jetLevel[iJetLevel]+"_DataComp_R="+Form("%.1f", jetRadius)+"_Area_@pT["+Form("%03.0f", PtCutLow)+","+Form("%03.0f", PtCutHigh)+"]");
  TString* pdfName_ratio = new TString("jet_"+jetType[iJetType]+"_"+jetLevel[iJetLevel]+"_DataComp_R="+Form("%.1f", jetRadius)+"_Area_@pT["+Form("%03.0f", PtCutLow)+","+Form("%03.0f", PtCutHigh)+"]_ratio");

  TString textContext(contextDatasetComp(jetRadius, PtRange, "pt"));

  Draw_TH1_Histograms_in_one(H1D_jetArea_rebinned, DatasetsNames, nDatasets, textContext, pdfName, texJetArea, texJetNormAreaYield, texCollisionDataInfo, "");

  if (divideSuccess == true) {
    Draw_TH1_Histograms_in_one(H1D_jetArea_rebinned_ratios, DatasetsNames, nDatasets, textContext, pdfName_ratio, texJetArea, texRatioDatasets, texCollisionDataInfo, "autoratio,avoidFirst");
  }
  else {
    cout << "Divide failed in Draw_Area_PtIntegrated_DatasetComparison" << endl;
  }
}


void Draw_JetTRDratio_vs_JetEta(int iDataset) {
  const Int_t nHistCollection = 1;
  const TString histCollectionLegend[nHistCollection] = {""};

  TH3D* H3D_jetPtjetEtajetTRDratio;
  TH2D* H2D_jetEtajetTRDratio[nHistCollection];
  TH2D* H2D_jetEtajetTRDratio_rebinned[nHistCollection];
  
  H3D_jetPtjetEtajetTRDratio = (TH3D*)((TH3D*)file_O2Analysis_list[iDataset]->Get(analysisWorkflow+"/h3_jet_pt_jet_eta_jet_nTRDtracks_ratio"))->Clone("Draw_JetTRDratio_vs_JetEta"+Datasets[iDataset]);
  H3D_jetPtjetEtajetTRDratio->Sumw2();

  H2D_jetEtajetTRDratio[0] = (TH2D*)H3D_jetPtjetEtajetTRDratio->Project3D(Datasets[iDataset]+"_jetEta_jetTRratio_zy"); //can't use letter D in this or it seems to replace the histogram in current pad (see documentation of ProjectionX function. Isn't mentioned in project3D sadly)
  // H2D_jetEtajetTRDratio_rebinned = (TH2D*)H2D_jetEtajetTRDratio->RebinY(2.,"H1D_jetEta_jetTRD_rebinned_"+RadiusLegend[iRadius]);

  TString* pdfName = new TString("jet_"+jetType[iJetType]+"_"+jetLevel[iJetLevel]+"_"+Datasets[iDataset]+"_JetTRDratio-vs-Eta");

  TString textContext("#splitline{"+*texDatasetsComparisonCommonDenominator+" "+DatasetsNames[iDataset]+"}{2023 QC}");

  Draw_TH2_Histograms(H2D_jetEtajetTRDratio, histCollectionLegend, nHistCollection, textContext, pdfName, texEtaX, texJetNTrdTracksRatio, texCollisionDataInfo, "");
}

void Draw_JetTRDcount_vs_JetEta(int iDataset) {
  const Int_t nHistCollection = 1;
  const TString histCollectionLegend[nHistCollection] = {""};

  TH3D* H3D_jetPtjetEtajetTRD;
  TH2D* H2D_jetEtajetTRD[nHistCollection];
  TH2D* H2D_jetEtajetTRD_rebinned[nHistCollection];
  
  H3D_jetPtjetEtajetTRD = (TH3D*)((TH3D*)file_O2Analysis_list[iDataset]->Get(analysisWorkflow+"/h3_jet_pt_jet_eta_jet_nTRDtracks"))->Clone("Draw_JetTRDcount_vs_JetEta"+Datasets[iDataset]);
  H3D_jetPtjetEtajetTRD->Sumw2();

  H2D_jetEtajetTRD[0] = (TH2D*)H3D_jetPtjetEtajetTRD->Project3D(Datasets[iDataset]+"_jetEta_jetTR_zy");
  // // H2D_jetEtajetTRD_rebinned = (TH2D*)H2D_jetEtajetTRD->RebinY(2.,"H1D_jetEta_jetTRD_rebinned_"+RadiusLegend[iRadius]);

  TString* pdfName = new TString("jet_"+jetType[iJetType]+"_"+jetLevel[iJetLevel]+"_"+Datasets[iDataset]+"_JetTRD-vs-Eta");

  TString textContext("#splitline{"+*texDatasetsComparisonCommonDenominator+" "+DatasetsNames[iDataset]+"}{2023 QC}");

  Draw_TH2_Histograms(H2D_jetEtajetTRD, histCollectionLegend, nHistCollection, textContext, pdfName, texEtaX, texJetNTrdTracks, texCollisionDataInfo, "");
}

void Draw_Pt_ratio_etaNeg_etaPos_TRDonly_vs_noTRD(int iDataset, float* etaRange) {
  const Int_t nHistCollection = 1;
  const TString histCollectionLegend[nHistCollection] = {""};

  TH3D* H3D_jetPtjetEtajetTRD;
  TH1D* H1D_jetPt_left_TRDonly[nHistCollection];
  TH1D* H1D_jetPt_right_TRDonly[nHistCollection];
  TH1D* H1D_jetPt_left_TRDonly_rebinned[nHistCollection];
  TH1D* H1D_jetPt_right_TRDonly_rebinned[nHistCollection];
  TH1D* H1D_jetPt_TRDonly_rebinned_ratios[nHistCollection];

  H3D_jetPtjetEtajetTRD = (TH3D*)((TH3D*)file_O2Analysis_list[iDataset]->Get(analysisWorkflow+"/h3_jet_pt_jet_eta_jet_nTRDtracks"))->Clone("Draw_Pt_ratio_etaNeg_etaPos_JetTRDonly"+Datasets[iDataset]+Form("%.1f", etaRange[0])+"<eta<"+Form("%.1f", etaRange[1]));
  H3D_jetPtjetEtajetTRD->Sumw2();

  TH1D* H1D_jetPt_left_noTRD[nHistCollection];
  TH1D* H1D_jetPt_right_noTRD[nHistCollection];
  TH1D* H1D_jetPt_left_noTRD_rebinned[nHistCollection];
  TH1D* H1D_jetPt_right_noTRD_rebinned[nHistCollection];
  TH1D* H1D_jetPt_noTRD_rebinned_ratios[nHistCollection];

  H3D_jetPtjetEtajetTRD = (TH3D*)((TH3D*)file_O2Analysis_list[iDataset]->Get(analysisWorkflow+"/h3_jet_pt_jet_eta_jet_nTRDtracks"))->Clone("Draw_Pt_ratio_etaNeg_etaPos_JetNoTRD"+Datasets[iDataset]+Form("%.1f", etaRange[0])+"<eta<"+Form("%.1f", etaRange[1]));
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
  H1D_jetPt_left_TRDonly[0] = (TH1D*)H3D_jetPtjetEtajetTRD->ProjectionX("jetPt_TRDonly_left"+Datasets[iDataset]+Form("%.1f", etaRange[0])+"<eta<"+Form("%.1f", etaRange[1]), ibinEta_low, ibinEta_zero-1, 2, -1); // asking bin number >= 2 avoids TRD match == 0; it's good to also keep overflow bin with -1 for upper limit as they do have TRD matches
  H1D_jetPt_right_TRDonly[0] = (TH1D*)H3D_jetPtjetEtajetTRD->ProjectionX("jetPt_TRDonly_right"+Datasets[iDataset]+Form("%.1f", etaRange[0])+"<eta<"+Form("%.1f", etaRange[1]), ibinEta_zero, ibinEta_high, 2, -1); // asking bin number >= 2 avoids TRD match == 0; it's good to also keep overflow bin with -1 for upper limit as they do have TRD matches
  H1D_jetPt_left_TRDonly_rebinned[0] = (TH1D*)H1D_jetPt_left_TRDonly[0]->Rebin(1.,"jetPt_TRDonly_left_rebinned"+Datasets[iDataset]+Form("%.1f", etaRange[0])+"<eta<"+Form("%.1f", etaRange[1]));
  H1D_jetPt_right_TRDonly_rebinned[0] = (TH1D*)H1D_jetPt_right_TRDonly[0]->Rebin(1.,"jetPt_TRDonly_right_rebinned"+Datasets[iDataset]+Form("%.1f", etaRange[0])+"<eta<"+Form("%.1f", etaRange[1]));

  NormaliseYieldToNEvents(H1D_jetPt_left_TRDonly_rebinned[0], GetNEventsSel8(file_O2Analysis_list[iDataset]));
  NormaliseYieldToNEvents(H1D_jetPt_right_TRDonly_rebinned[0], GetNEventsSel8(file_O2Analysis_list[iDataset]));

  H1D_jetPt_TRDonly_rebinned_ratios[0] = (TH1D*)H1D_jetPt_left_TRDonly_rebinned[0]->Clone("jetPt_TRDonly_rebinned_ratios"+Datasets[iDataset]+Form("%.1f", etaRange[0])+"<eta<"+Form("%.1f", etaRange[1]));
  H1D_jetPt_TRDonly_rebinned_ratios[0]->Reset("M");
  divideSuccess = H1D_jetPt_TRDonly_rebinned_ratios[0]->Divide(H1D_jetPt_right_TRDonly_rebinned[0], H1D_jetPt_left_TRDonly_rebinned[0]);

  TString* pdfNameA = new TString("jet_"+jetType[iJetType]+"_"+jetLevel[iJetLevel]+"_"+Datasets[iDataset]+"_pT_etaRightLeftRatio_JetTRDonly");

  TString textContextA(contextDatasetRadiusComp(iDataset, etaRange, "eta"));

  if (divideSuccess == true) {
    Draw_TH1_Histograms_in_one(H1D_jetPt_TRDonly_rebinned_ratios, histCollectionLegend, nHistCollection, textContextA, pdfNameA, texPtX, texJetRatioEtaComparison, texCollisionDataInfo, "standardratio");
  }
  else {
    cout << "Divide failed in Draw_Pt_ratio_etaNeg_etaPos_RadiusComparison 1" << endl;
  }


  divideSuccess = false;
  //// jets without TRDmatched tracks only 

  H1D_jetPt_left_noTRD[0] = (TH1D*)H3D_jetPtjetEtajetTRD->ProjectionX("jetPt_noTRD_left"+Datasets[iDataset]+Form("%.1f", etaRange[0])+"<eta<"+Form("%.1f", etaRange[1]), ibinEta_low, ibinEta_zero-1, 0, 1); // asking bin number >= 2 avoids TRD match == 0; it's good to also keep overflow bin with -1 for upper limit as they do have TRD matches
  H1D_jetPt_right_noTRD[0] = (TH1D*)H3D_jetPtjetEtajetTRD->ProjectionX("jetPt_noTRD_right"+Datasets[iDataset]+Form("%.1f", etaRange[0])+"<eta<"+Form("%.1f", etaRange[1]), ibinEta_zero, ibinEta_high, 0, 1); // asking bin number >= 2 avoids TRD match == 0; it's good to also keep overflow bin with -1 for upper limit as they do have TRD matches
  H1D_jetPt_left_noTRD_rebinned[0] = (TH1D*)H1D_jetPt_left_noTRD[0]->Rebin(1.,"jetPt_noTRD_left_rebinned"+Datasets[iDataset]+Form("%.1f", etaRange[0])+"<eta<"+Form("%.1f", etaRange[1]));
  H1D_jetPt_right_noTRD_rebinned[0] = (TH1D*)H1D_jetPt_right_noTRD[0]->Rebin(1.,"jetPt_noTRD_right_rebinned"+Datasets[iDataset]+Form("%.1f", etaRange[0])+"<eta<"+Form("%.1f", etaRange[1]));

  NormaliseYieldToNEvents(H1D_jetPt_left_noTRD_rebinned[0], GetNEventsSel8(file_O2Analysis_list[iDataset]));
  NormaliseYieldToNEvents(H1D_jetPt_right_noTRD_rebinned[0], GetNEventsSel8(file_O2Analysis_list[iDataset]));

  H1D_jetPt_noTRD_rebinned_ratios[0] = (TH1D*)H1D_jetPt_left_noTRD_rebinned[0]->Clone("jetPt_noTRD_rebinned_ratios"+Datasets[iDataset]+Form("%.1f", etaRange[0])+"<eta<"+Form("%.1f", etaRange[1]));
  H1D_jetPt_noTRD_rebinned_ratios[0]->Reset("M");
  divideSuccess = H1D_jetPt_noTRD_rebinned_ratios[0]->Divide(H1D_jetPt_right_noTRD_rebinned[0], H1D_jetPt_left_noTRD_rebinned[0]);

  TString* pdfNameB = new TString("jet_"+jetType[iJetType]+"_"+jetLevel[iJetLevel]+"_"+Datasets[iDataset]+"_pT_etaRightLeftRatio_JetNoTRD");

  TString textContextB(contextDatasetRadiusComp(iDataset, etaRange, "eta"));

  if (divideSuccess == true) {
    Draw_TH1_Histograms_in_one(H1D_jetPt_noTRD_rebinned_ratios, histCollectionLegend, nHistCollection, textContextB, pdfNameB, texPtX, texJetRatioEtaComparison, texCollisionDataInfo, "standardratio");
  }
  else {
    cout << "Divide failed in Draw_Pt_ratio_etaNeg_etaPos_RadiusComparison 2" << endl;
  }


  //// ratio between no trd and only trd
  H1D_jetPt_noTRD_vs_onlyTRD_ratios[0] = (TH1D*)H1D_jetPt_left_noTRD_rebinned[0]->Clone("jetPt_noTRD_vs_onlyTRD_ratios"+Datasets[iDataset]+Form("%.1f", etaRange[0])+"<eta<"+Form("%.1f", etaRange[1]));
  H1D_jetPt_noTRD_vs_onlyTRD_ratios[0]->Reset("M");
  divideSuccess = H1D_jetPt_noTRD_vs_onlyTRD_ratios[0]->Divide(H1D_jetPt_TRDonly_rebinned_ratios[0], H1D_jetPt_noTRD_rebinned_ratios[0]);

  TString* pdfNameC = new TString("jet_"+jetType[iJetType]+"_"+jetLevel[iJetLevel]+"_"+Datasets[iDataset]+"_pT_etaRightLeftRatio_Jet_NoTRD_vs_OnlyTRD_ratio");

  TString textContextC(contextDatasetRadiusComp(iDataset, etaRange, "eta"));

  if (divideSuccess == true) {
    Draw_TH1_Histograms_in_one(H1D_jetPt_noTRD_vs_onlyTRD_ratios, histCollectionLegend, nHistCollection, textContextC, pdfNameC, texPtX, texJetRatioTRDvsNoTRD, texCollisionDataInfo, "standardratio");
  }
  else {
    cout << "Divide failed in Draw_Pt_ratio_etaNeg_etaPos_RadiusComparison 3" << endl;
  }

}


// To-do list:
// - implement carolina's macro of the gpad thing to automate the division of canvas based on how many plots one wants ?
// - 1D graphs: 2023 QC displayed on top of texCollisionDataInfo; but 2D graphs fine --> doesnt like the nested splitline


// option e for projection computes errors

