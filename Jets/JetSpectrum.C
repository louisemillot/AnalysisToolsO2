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
#include <RooUnfold.h>
#include "RooUnfoldResponse.h"
#include "RooUnfoldBayes.h"
#include "RooUnfoldBinByBin.h"
#include "RooUnfoldSvd.h"

//My Libraries
#include "./JetSpectrum_settings.h"
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
TString contextDataset1D(int iDataset, float* variableRange, const char options[]);
TString contextDatasetCompAndRadiusAndVarRange(float jetRadius, float* variableRange, const char options[]);
TString contextPtRange(float* PtRange);
TString contextEtaRange(float* PtRange);
TString contextJetRadius(float jetRadius);
TString contextCentRange(float* centRange);

//////////// Pt Spectrum analysis functions
void Get_Pt_spectrum_bkgCorrected(TH1D* &H1D_jetPt_raw, int iDataset, int iRadius, float* centRange, const char options[]);
void Get_Pt_spectrum_bkgCorrected_genBinning(TH1D* &H1D_jetPt_raw, int iDataset, int iRadius, float* centRange, const char options[]);
void Get_Pt_spectrum_mcp(TH1D* &H1D_jetPt_mcp, int iDataset, int iRadius, const char options[]);
void Get_Pt_spectrum_mcp_recBinning(TH1D* &H1D_jetPt_mcp, int iDataset, int iRadius, const char options[]);
void Get_Pt_spectrum_mcdMatched(TH1D* &H1D_jetPt_mcd, int iDataset, int iRadius, const char options[]);
void Get_Pt_spectrum_mcpMatched(TH1D* &H1D_jetPt_mcp, int iDataset, int iRadius, const char options[]);
void Get_Pt_spectrum_unfolded(TH1D* &H1D_jetPt_unfolded, int iDataset, int iRadius, float* centRange, const char options[]);
void Get_Pt_spectrum_bkgCorrected_preWidthScaling(TH1D* &H1D_jetPt_raw, int iDataset, int iRadius, float* centRange, const char options[]);
void Get_Pt_spectrum_bkgCorrected_genBinning_preWidthScaling(TH1D* &H1D_jetPt_raw, int iDataset, int iRadius, float* centRange, const char options[]);
void Get_Pt_spectrum_mcp_preWidthScaling(TH1D* &H1D_jetPt_mcp, int iDataset, int iRadius, const char options[]);
void Get_Pt_spectrum_mcp_recBinning_preWidthScaling(TH1D* &H1D_jetPt_mcp, int iDataset, int iRadius, const char options[]);
void Get_Pt_spectrum_mcdMatched_preWidthScaling(TH1D* &H1D_jetPt_mcd, int iDataset, int iRadius, const char options[]);
void Get_Pt_spectrum_mcpMatched_preWidthScaling(TH1D* &H1D_jetPt_mcpMatched, int iDataset, int iRadius, const char options[]);
void Get_Pt_spectrum_unfolded_preWidthScaling(TH1D* &H1D_jetPt_unfolded, int iDataset, int iRadius, float* centRange, const char options[]);
void Get_Pt_spectrum_mcp_folded_preWidthScaling(TH1D* &H1D_jetPt_mcp_folded, int iDataset, int iRadius, float* centRange, const char options[]);
void Get_Pt_spectrum_mcp_folded(TH1D* &H1D_jetPt_mcp_folded, int iDataset, int iRadius, float* centRange, const char options[]);
void Get_Pt_spectrum_mcp_folded_genBinning_preWidthScaling(TH1D* &H1D_jetPt_mcp_folded, int iDataset, int iRadius, float* centRange, const char options[]);
void Get_Pt_spectrum_mcp_folded_genBinning(TH1D* &H1D_jetPt_mcp_folded, int iDataset, int iRadius, float* centRange, const char options[]);

void Get_PtResponseMatrix_Fluctuations(TH2D* &H2D_jetPtResponseMatrix_fluctuations, int iDataset, int iRadius, float* centRange);
void Get_PtResponseMatrix_detectorResponse(TH2D* &H2D_jetPtResponseMatrix_detectorResponse, int iDataset, int iRadius);
void Get_PtResponseMatrix_DetectorAndFluctuationsCombined(TH2D* &H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined, TH2D* H2D_jetPtResponseMatrix_detectorResponse, TH2D* H2D_jetPtResponseMatrix_fluctuations, int iDataset, int iRadius, float* centRange);

void Draw_ResponseMatrices_Fluctuations(int iDataset, int iRadius);
void Draw_ResponseMatrices_detectorResponse(int iDataset, int iRadius);
void Draw_ResponseMatrices_DetectorAndFluctuationsCombined(int iDataset, int iRadius);

void Draw_Pt_spectrum_raw(int iDataset, int iRadius, const char options[]);
void Draw_Pt_spectrum_mcp(int iDataset, int iRadius, const char options[]);
void Draw_Pt_spectrum_mcdMatched(int iDataset, int iRadius, const char options[]);
void Draw_Pt_spectrum_unfolded_FluctResponseOnly(int iDataset, int iRadius, const char options[]);
void Draw_Pt_spectrum_unfolded(int iDataset, int iRadius, const char options[]);
void Draw_Pt_TestSpectrum_unfolded(int iDataset, int iRadius, const char options[]);

void Draw_Pt_efficiency_jets(int iDataset, int iRadius, const char options[]);


/////////////////////////////////////////////////////
///////////////////// Main Macro ////////////////////
/////////////////////////////////////////////////////

void JetSpectrum() {
  // Load necessary libraries
  LoadLibs();
  // Set the default style
  SetStyle();

  // TString* SaveAs_Title = new TString("");
  TString* texXtitle = new TString("");
  TString* texYtitle = new TString("");
  // TString* Extra = new TString("");

  int iDataset = 0;

  // float etaRangeSym[2] = {-0.5, 0.5};
  // float etaRangeNeg[2] = {-0.5, 0};
  // float etaRangePos[2] = {0, 0.5};

  float jetRadiusForDataComp = 0.4;
  float jetR02 = 0.2;
  float jetR06 = 0.6;
  
  // float ptRange[2] = {1, 200};

  int iRadius = 1;

  // Draw_Pt_spectrum_raw(iDataset, iRadius, "evtNorm");

  Draw_ResponseMatrices_Fluctuations(iDataset, iRadius);
  Draw_ResponseMatrices_detectorResponse(iDataset, iRadius);
  Draw_ResponseMatrices_DetectorAndFluctuationsCombined(iDataset, iRadius);

  // // Draw_Pt_spectrum_unfolded_FluctResponseOnly(iDataset, iRadius, "evtNorm"); // NOT FIXED YET - result meaningless

  // Draw_Pt_spectrum_mcp(iDataset, iRadius, "evtNorm");
  // Draw_Pt_spectrum_mcdMatched(iDataset, iRadius, "evtNorm");

  // Draw_Pt_efficiency_jets(iDataset, iRadius, "");

  Draw_Pt_spectrum_unfolded(iDataset, iRadius, "evtNorm");
  // Draw_Pt_TestSpectrum_unfolded(iDataset, iRadius, "evtNorm");


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

// TString contextDatasetRadiusCompAndVarRange(int iDataset, float* variableRange, const char options[]){
//   TString texcontextDatasetRadiusCompAndVarRange;
//   if (strstr(options, "pt") != NULL) { //  || strstr(options, "ratio") != NULL not sure why I had this here
//     texcontextDatasetRadiusCompAndVarRange = "#splitline{"+*texDatasetsComparisonCommonDenominator+" "+DatasetsNames[iDataset]+"}{#splitline{2023 QC}{"+contextPtRange(variableRange)+"}}";
//   }
//   if (strstr(options, "eta") != NULL) { //  || strstr(options, "ratio") != NULL not sure why I had this here
//     texcontextDatasetRadiusCompAndVarRange = "#splitline{"+*texDatasetsComparisonCommonDenominator+" "+DatasetsNames[iDataset]+"}{#splitline{2023 QC}{"+contextEtaRange(variableRange)+"}}";
//   }
//   if (strstr(options, "cent") != NULL) { //  || strstr(options, "ratio") != NULL not sure why I had this here
//     texcontextDatasetRadiusCompAndVarRange = "#splitline{"+*texDatasetsComparisonCommonDenominator+" "+DatasetsNames[iDataset]+"}{#splitline{2023 QC}{"+contextCentRange(variableRange)+"}}";
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

// TString contextCentRange(float* centRange){
//   std::stringstream ss;
//   ss << centRange[0] << " < #cent < " << centRange[1];
//   TString textContext((TString)ss.str());
//   return textContext;
// }

// void CentralityLegend(TString* centralityLegend, const float* arrayCentralityBinning, int nCentralityBins){
//   std::stringstream ss;
//   ss.precision(2);
//   for(int iCentralityBin = 0; iCentralityBin < nCentralityBins; iCentralityBin++){
//     ss << "" << arrayCentralityBinning[iCentralityBin] << "-" << arrayCentralityBinning[iCentralityBin+1] << "%";
//     centralityLegend[iCentralityBin] = (TString)ss.str();
//     ss.str("");
//     ss.clear();
//   }
// }

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////// Spectrum getting functions ///////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


void Get_Pt_spectrum_bkgCorrected_preWidthScaling(TH1D* &H1D_jetPt, int iDataset, int iRadius, float* centRange, const char options[]) {
  TH3D* H3D_jetRjetPtcolCent;
  TH1D* H1D_jetPt_defaultBin;
  // TH1D* H1D_jetPt_raw[nRadius];

  H3D_jetRjetPtcolCent = (TH3D*)((TH3D*)file_O2Analysis_list[iDataset]->Get(analysisWorkflowData+"/h3_jet_r_jet_pt_centrality_rhoareasubtracted"))->Clone("Get_Pt_spectrum_bkgCorrected"+Datasets[iDataset]+Form("%.1f", centRange[0])+"<cent<"+Form("%.1f", centRange[1]));
  H3D_jetRjetPtcolCent->Sumw2();
  int binCentEdges[2] = {H3D_jetRjetPtcolCent->GetZaxis()->FindBin(centRange[0]), H3D_jetRjetPtcolCent->GetZaxis()->FindBin(centRange[1])};
  if (binCentEdges[0] == 0) 
    cout << "WARNING: Get_Pt_spectrum_bkgCorrected is counting the underflow with the chosen centRange" << endl;
  if (binCentEdges[1] == H3D_jetRjetPtcolCent->GetZaxis()->GetNbins()+1) 
    cout << "WARNING: Get_Pt_spectrum_bkgCorrected is counting the overflow with the chosen centRange" << endl;

  int ibinJetRadius = H3D_jetRjetPtcolCent->GetXaxis()->FindBin(arrayRadius[iRadius]+GLOBAL_epsilon);

  H1D_jetPt_defaultBin = (TH1D*)H3D_jetRjetPtcolCent->ProjectionY("jetPt_bkgCorrected_"+RadiusLegend[iRadius]+Form("%.1f", centRange[0])+"<cent<"+Form("%.1f", centRange[1]), ibinJetRadius, ibinJetRadius, binCentEdges[0], binCentEdges[1]);
  H1D_jetPt = (TH1D*)H1D_jetPt_defaultBin->Rebin(nBinPtJetsRec[iRadius],"jetPt_bkgCorrected_rebinned_"+RadiusLegend[iRadius]+Form("%.1f", centRange[0])+"<cent<"+Form("%.1f", centRange[1]), ptBinsJetsRec[iRadius]);

  if (strstr(options, "evtNorm") != NULL) {
    // NormaliseYieldToNEvents(H1D_jetPt, GetNEventsSelectedCentrality_JetFramework(file_O2Analysis_list[iDataset], centRange[0], centRange[1], trainId));
    NormaliseRawHistToNEvents(H1D_jetPt, GetNEventsSelectedCentrality_JetFramework(file_O2Analysis_list[iDataset], centRange[0], centRange[1], trainId));
  }
  cout << "GetNEventsSelectedCentrality_JetFramework(file_O2Analysis_list[iDataset]) = " << GetNEventsSelectedCentrality_JetFramework(file_O2Analysis_list[iDataset], centRange[0], centRange[1], trainId) << endl;
  if (strstr(options, "entriesNorm") != NULL) {
    NormaliseYieldToNEntries(H1D_jetPt);
  }
}


void Get_Pt_spectrum_bkgCorrected_genBinning_preWidthScaling(TH1D* &H1D_jetPt, int iDataset, int iRadius, float* centRange, const char options[]) {
  TH3D* H3D_jetRjetPtcolCent;
  TH1D* H1D_jetPt_defaultBin;
  // TH1D* H1D_jetPt_raw[nRadius];

  H3D_jetRjetPtcolCent = (TH3D*)((TH3D*)file_O2Analysis_list[iDataset]->Get(analysisWorkflowData+"/h3_jet_r_jet_pt_centrality_rhoareasubtracted"))->Clone("Get_Pt_spectrum_bkgCorrected_genBinning"+Datasets[iDataset]+Form("%.1f", centRange[0])+"<cent<"+Form("%.1f", centRange[1]));
  H3D_jetRjetPtcolCent->Sumw2();
  int binCentEdges[2] = {H3D_jetRjetPtcolCent->GetZaxis()->FindBin(centRange[0]), H3D_jetRjetPtcolCent->GetZaxis()->FindBin(centRange[1])};
  if (binCentEdges[0] == 0) 
    cout << "WARNING: Get_Pt_spectrum_bkgCorrected_genBinning is counting the underflow with the chosen centRange" << endl;
  if (binCentEdges[1] == H3D_jetRjetPtcolCent->GetZaxis()->GetNbins()+1) 
    cout << "WARNING: Get_Pt_spectrum_bkgCorrected_genBinning is counting the overflow with the chosen centRange" << endl;

  int ibinJetRadius = H3D_jetRjetPtcolCent->GetXaxis()->FindBin(arrayRadius[iRadius]+GLOBAL_epsilon);

  H1D_jetPt_defaultBin = (TH1D*)H3D_jetRjetPtcolCent->ProjectionY("jetPt_bkgCorrected_genBinning_"+RadiusLegend[iRadius]+Form("%.1f", centRange[0])+"<cent<"+Form("%.1f", centRange[1]), ibinJetRadius, ibinJetRadius, binCentEdges[0], binCentEdges[1]);
  H1D_jetPt = (TH1D*)H1D_jetPt_defaultBin->Rebin(nBinPtJetsGen[iRadius],"jetPt_bkgCorrected_rebinned_genBinning_"+RadiusLegend[iRadius]+Form("%.1f", centRange[0])+"<cent<"+Form("%.1f", centRange[1]), ptBinsJetsGen[iRadius]);

  if (strstr(options, "evtNorm") != NULL) {
    // NormaliseYieldToNEvents(H1D_jetPt, GetNEventsSelectedCentrality_JetFramework(file_O2Analysis_list[iDataset], centRange[0], centRange[1], trainId));
    NormaliseRawHistToNEvents(H1D_jetPt, GetNEventsSelectedCentrality_JetFramework(file_O2Analysis_list[iDataset], centRange[0], centRange[1], trainId));
  }
  cout << "GetNEventsSelectedCentrality_JetFramework(file_O2Analysis_list[iDataset]) = " << GetNEventsSelectedCentrality_JetFramework(file_O2Analysis_list[iDataset], centRange[0], centRange[1], trainId) << endl;
  if (strstr(options, "entriesNorm") != NULL) {
    NormaliseYieldToNEntries(H1D_jetPt);
  }
}

void Get_Pt_spectrum_mcp_preWidthScaling(TH1D* &H1D_jetPt_mcp, int iDataset, int iRadius, const char options[]) {
  TH3D* H3D_jetRjetPtcolCent;
  TH1D* H1D_jetPt_mcp_defaultBin;
  // TH1D* H1D_jetPt_raw[nRadius];

  // //////// with centrality; sadly we don't have a mcp version that gives centrality info
  // H3D_jetRjetPtcolCent = (TH3D*)((TH3D*)file_O2Analysis_ppSimDetectorEffect[iCent]->Get(analysisWorkflow+"/h3_jet_r_jet_pt_centrality"))->Clone("Get_Pt_spectrum_mcp"+Datasets[iDataset]+Form("%.1f", centRange[0])+"<cent<"+Form("%.1f", centRange[1]));
  // H3D_jetRjetPtcolCent->Sumw2();

  // int binCentEdges[2] = {H3D_jetRjetPtcolCent->GetZaxis()->FindBin(centRange[0]), H3D_jetRjetPtcolCent->GetZaxis()->FindBin(centRange[1])};
  // if (binCentEdges[0] == 0) 
  //   cout << "WARNING: Get_Pt_spectrum_bkgCorrected is counting the underflow with the chosen centRange" << endl;
  // if (binCentEdges[1] == H3D_jetRjetPtcolCent->GetZaxis()->GetNbins()+1)
  //   cout << "WARNING: Get_Pt_spectrum_bkgCorrected is counting the overflow with the chosen centRange" << endl;
  

  // int ibinJetRadius = H3D_jetRjetPtcolCent->GetXaxis()->FindBin(arrayRadius[iRadius]+GLOBAL_epsilon);

  // H1D_jetPt_mcp_defaultBin = (TH1D*)H3D_jetRjetPtcolCent->ProjectionY("jetPt_mcp_"+RadiusLegend[iRadius]+Form("%.1f", centRange[0])+"<cent<"+Form("%.1f", centRange[1]), ibinJetRadius, ibinJetRadius, binCentEdges[0], binCentEdges[1]);
  // // H1D_jetPt_mcp_defaultBin = (TH1D*)H3D_jetRjetPtcolCent->ProjectionY("jetPt_mcp_"+RadiusLegend[iRadius]+Form("%.1f", centRange[0])+"<cent<"+Form("%.1f", centRange[1]), ibinJetRadius, ibinJetRadius, 0, -1);

  // H1D_jetPt_mcp = (TH1D*)H1D_jetPt_mcp_defaultBin->Rebin(nBinPtJetsGen[iRadius],"jetPt_mcp_rebinned_"+RadiusLegend[iRadius]+Form("%.1f", centRange[0])+"<cent<"+Form("%.1f", centRange[1]), ptBinsJetsGen[iRadius]);




  //////// without centrality; sadly we don't have a mcp version that gives centrality info, so gotta use one with all centralities (though one can edit the cent window in the jetfinderqa options)
  H3D_jetRjetPtcolCent = (TH3D*)((TH3D*)file_O2Analysis_ppSimDetectorEffect->Get(analysisWorkflowMC+"/h3_jet_r_part_jet_pt_part_jet_eta_part"))->Clone("Get_Pt_spectrum_mcp"+Datasets[iDataset]);
  H3D_jetRjetPtcolCent->Sumw2();

  int ibinJetRadius = H3D_jetRjetPtcolCent->GetXaxis()->FindBin(arrayRadius[iRadius]+GLOBAL_epsilon);

  H1D_jetPt_mcp_defaultBin = (TH1D*)H3D_jetRjetPtcolCent->ProjectionY("jetPt_mcp_"+RadiusLegend[iRadius]+Datasets[iDataset], ibinJetRadius, ibinJetRadius, 0, -1);

  H1D_jetPt_mcp = (TH1D*)H1D_jetPt_mcp_defaultBin->Rebin(nBinPtJetsGen[iRadius],"jetPt_mcp_rebinned_"+RadiusLegend[iRadius]+Datasets[iDataset], ptBinsJetsGen[iRadius]);



  if (strstr(options, "evtNorm") != NULL) {
    // NormaliseYieldToNEvents(H1D_jetPt_mcp, GetNEventsGen(file_O2Analysis_list[iDataset]));
    // NormaliseYieldToNEvents(H1D_jetPt_mcp, GetNEventsSelected_JetFramework_weighted(file_O2Analysis_ppSimDetectorEffect));
    NormaliseRawHistToNEvents(H1D_jetPt_mcp, GetNEventsSelected_JetFramework_weighted(file_O2Analysis_ppSimDetectorEffect));
    
    // NormaliseYieldToNEntries(H1D_jetPt_mcp);
    
  }
  if (strstr(options, "entriesNorm") != NULL) {
    NormaliseYieldToNEntries(H1D_jetPt_mcp);
  }
}


void Get_Pt_spectrum_mcp_recBinning_preWidthScaling(TH1D* &H1D_jetPt_mcp, int iDataset, int iRadius, const char options[]) {
  TH3D* H3D_jetRjetPtcolCent;
  TH1D* H1D_jetPt_mcp_defaultBin;
  // TH1D* H1D_jetPt_raw[nRadius];

  // //////// with centrality; sadly we don't have a mcp version that gives centrality info
  // H3D_jetRjetPtcolCent = (TH3D*)((TH3D*)file_O2Analysis_ppSimDetectorEffect[iCent]->Get(analysisWorkflow+"/h3_jet_r_jet_pt_centrality"))->Clone("Get_Pt_spectrum_mcp"+Datasets[iDataset]+Form("%.1f", centRange[0])+"<cent<"+Form("%.1f", centRange[1]));
  // H3D_jetRjetPtcolCent->Sumw2();

  // int binCentEdges[2] = {H3D_jetRjetPtcolCent->GetZaxis()->FindBin(centRange[0]), H3D_jetRjetPtcolCent->GetZaxis()->FindBin(centRange[1])};
  // if (binCentEdges[0] == 0) 
  //   cout << "WARNING: Get_Pt_spectrum_bkgCorrected is counting the underflow with the chosen centRange" << endl;
  // if (binCentEdges[1] == H3D_jetRjetPtcolCent->GetZaxis()->GetNbins()+1)
  //   cout << "WARNING: Get_Pt_spectrum_bkgCorrected is counting the overflow with the chosen centRange" << endl;
  

  // int ibinJetRadius = H3D_jetRjetPtcolCent->GetXaxis()->FindBin(arrayRadius[iRadius]+GLOBAL_epsilon);

  // H1D_jetPt_mcp_defaultBin = (TH1D*)H3D_jetRjetPtcolCent->ProjectionY("jetPt_mcp_"+RadiusLegend[iRadius]+Form("%.1f", centRange[0])+"<cent<"+Form("%.1f", centRange[1]), ibinJetRadius, ibinJetRadius, binCentEdges[0], binCentEdges[1]);
  // // H1D_jetPt_mcp_defaultBin = (TH1D*)H3D_jetRjetPtcolCent->ProjectionY("jetPt_mcp_"+RadiusLegend[iRadius]+Form("%.1f", centRange[0])+"<cent<"+Form("%.1f", centRange[1]), ibinJetRadius, ibinJetRadius, 0, -1);

  // H1D_jetPt_mcp = (TH1D*)H1D_jetPt_mcp_defaultBin->Rebin(nBinPtJetsGen[iRadius],"jetPt_mcp_rebinned_"+RadiusLegend[iRadius]+Form("%.1f", centRange[0])+"<cent<"+Form("%.1f", centRange[1]), ptBinsJetsGen[iRadius]);




  //////// without centrality; sadly we don't have a mcp version that gives centrality info, so gotta use one with all centralities (though one can edit the cent window in the jetfinderqa options)
  H3D_jetRjetPtcolCent = (TH3D*)((TH3D*)file_O2Analysis_ppSimDetectorEffect->Get(analysisWorkflowMC+"/h3_jet_r_part_jet_pt_part_jet_eta_part"))->Clone("Get_Pt_spectrum_mcp_recBinning_preWidthScaling"+Datasets[iDataset]);
  H3D_jetRjetPtcolCent->Sumw2();

  int ibinJetRadius = H3D_jetRjetPtcolCent->GetXaxis()->FindBin(arrayRadius[iRadius]+GLOBAL_epsilon);

  H1D_jetPt_mcp_defaultBin = (TH1D*)H3D_jetRjetPtcolCent->ProjectionY("jetPt_mcp_rebinRec_"+RadiusLegend[iRadius]+Datasets[iDataset], ibinJetRadius, ibinJetRadius, 0, -1);

  H1D_jetPt_mcp = (TH1D*)H1D_jetPt_mcp_defaultBin->Rebin(nBinPtJetsRec[iRadius],"jetPt_mcp_rebinnedRec_"+RadiusLegend[iRadius]+Datasets[iDataset], ptBinsJetsRec[iRadius]);



  if (strstr(options, "evtNorm") != NULL) {
    // NormaliseYieldToNEvents(H1D_jetPt_mcp, GetNEventsGen(file_O2Analysis_list[iDataset]));
    // NormaliseYieldToNEvents(H1D_jetPt_mcp, GetNEventsSelected_JetFramework_weighted(file_O2Analysis_ppSimDetectorEffect));
    NormaliseRawHistToNEvents(H1D_jetPt_mcp, GetNEventsSelected_JetFramework_weighted(file_O2Analysis_ppSimDetectorEffect));
    
    // NormaliseYieldToNEntries(H1D_jetPt_mcp);
    
  }
  if (strstr(options, "entriesNorm") != NULL) {
    NormaliseYieldToNEntries(H1D_jetPt_mcp);
  }
}

void Get_Pt_spectrum_mcpMatched_preWidthScaling(TH1D* &H1D_jetPt_mcpMatched, int iDataset, int iRadius, const char options[]) {
  TH3D* H3D_jetRjetPtcolCent;
  TH1D* H1D_jetPt_mcpMatched_defaultBin;
  // TH1D* H1D_jetPt_raw[nRadius];

  // //////// with centrality; sadly we don't have a mcp version that gives centrality info
  // H3D_jetRjetPtcolCent = (TH3D*)((TH3D*)file_O2Analysis_ppSimDetectorEffect[iCent]->Get(analysisWorkflow+"/h3_jet_r_jet_pt_centrality"))->Clone("Get_Pt_spectrum_mcp"+Datasets[iDataset]+Form("%.1f", centRange[0])+"<cent<"+Form("%.1f", centRange[1]));
  // H3D_jetRjetPtcolCent->Sumw2();

  // int binCentEdges[2] = {H3D_jetRjetPtcolCent->GetZaxis()->FindBin(centRange[0]), H3D_jetRjetPtcolCent->GetZaxis()->FindBin(centRange[1])};
  // if (binCentEdges[0] == 0) 
  //   cout << "WARNING: Get_Pt_spectrum_bkgCorrected is counting the underflow with the chosen centRange" << endl;
  // if (binCentEdges[1] == H3D_jetRjetPtcolCent->GetZaxis()->GetNbins()+1)
  //   cout << "WARNING: Get_Pt_spectrum_bkgCorrected is counting the overflow with the chosen centRange" << endl;
  

  // int ibinJetRadius = H3D_jetRjetPtcolCent->GetXaxis()->FindBin(arrayRadius[iRadius]+GLOBAL_epsilon);

  // H1D_jetPt_mcpMatched_defaultBin = (TH1D*)H3D_jetRjetPtcolCent->ProjectionY("jetPt_mcpMatched_"+RadiusLegend[iRadius]+Form("%.1f", centRange[0])+"<cent<"+Form("%.1f", centRange[1]), ibinJetRadius, ibinJetRadius, binCentEdges[0], binCentEdges[1]);
  // // H1D_jetPt_mcpMatched_defaultBin = (TH1D*)H3D_jetRjetPtcolCent->ProjectionY("jetPt_mcpMatched_"+RadiusLegend[iRadius]+Form("%.1f", centRange[0])+"<cent<"+Form("%.1f", centRange[1]), ibinJetRadius, ibinJetRadius, 0, -1);

  // H1D_jetPt_mcpMatched = (TH1D*)H1D_jetPt_mcpMatched_defaultBin->Rebin(nBinPtJetsGen[iRadius],"jetPt_mcpMatched_rebinned_"+RadiusLegend[iRadius]+Form("%.1f", centRange[0])+"<cent<"+Form("%.1f", centRange[1]), ptBinsJetsGen[iRadius]);






  //////// without centrality; sadly we don't have a mcp version that gives centrality info, so gotta use one with all centralities (though one can edit the cent window in the jetfinderqa options)
  // H3D_jetRjetPtcolCent = (TH3D*)((TH3D*)file_O2Analysis_ppSimDetectorEffect->Get(analysisWorkflow+"/h3_jet_r_jet_pt_jet_eta"))->Clone("Get_Pt_spectrum_mcpMatched"+Datasets[iDataset]);
  // H3D_jetRjetPtcolCent->Sumw2();
  H3D_jetRjetPtcolCent = (TH3D*)((TH3D*)file_O2Analysis_ppSimDetectorEffect->Get(analysisWorkflowMC+"/h3_jet_r_jet_pt_tag_jet_pt_base_matchedgeo"))->Clone("Get_Pt_spectrum_mcpMatched"+Datasets[iDataset]);// base is mcd in jetfinderQA as of 06/2024, thus tag is mcp, and so hist is (x=r, y=mcp, z=mcd)
  H3D_jetRjetPtcolCent->Sumw2();


  int ibinJetRadius = H3D_jetRjetPtcolCent->GetXaxis()->FindBin(arrayRadius[iRadius]+GLOBAL_epsilon);

  H1D_jetPt_mcpMatched_defaultBin = (TH1D*)H3D_jetRjetPtcolCent->ProjectionY("jetPt_mcpMatched_"+RadiusLegend[iRadius]+Datasets[iDataset], ibinJetRadius, ibinJetRadius, 0, -1);

  // should be rebinned using rec bins, because in the unfolding it'll be given to RooUnfoldResponse(mcpMatched, mcp, response) to get the efficiency 
  H1D_jetPt_mcpMatched = (TH1D*)H1D_jetPt_mcpMatched_defaultBin->Rebin(nBinPtJetsGen[iRadius],"jetPt_mcpMatched_rebinned_"+RadiusLegend[iRadius]+Datasets[iDataset], ptBinsJetsGen[iRadius]);







  if (strstr(options, "evtNorm") != NULL) {
    // NormaliseYieldToNEvents(H1D_jetPt_mcpMatched, GetNEventsGen(file_O2Analysis_list[iDataset]));
    // NormaliseYieldToNEvents(H1D_jetPt_mcpMatched, GetNEventsSelected_JetFramework_weighted(file_O2Analysis_ppSimDetectorEffect));
    NormaliseRawHistToNEvents(H1D_jetPt_mcpMatched, GetNEventsSelected_JetFramework_weighted(file_O2Analysis_ppSimDetectorEffect));
    // NormaliseYieldToNEntries(H1D_jetPt_mcpMatched);
  }
  if (strstr(options, "entriesNorm") != NULL) {
    NormaliseYieldToNEntries(H1D_jetPt_mcpMatched);
  }

}

void Get_Pt_spectrum_mcdMatched_preWidthScaling(TH1D* &H1D_jetPt_mcdMatched, int iDataset, int iRadius, const char options[]) {
  TH3D* H3D_jetRjetPtcolCent;
  TH1D* H1D_jetPt_mcdMatched_defaultBin;
  // TH1D* H1D_jetPt_raw[nRadius];

  // //////// with centrality; sadly we don't have a mcd version that gives centrality info
  // H3D_jetRjetPtcolCent = (TH3D*)((TH3D*)file_O2Analysis_ppSimDetectorEffect[iCent]->Get(analysisWorkflow+"/h3_jet_r_jet_pt_centrality"))->Clone("Get_Pt_spectrum_mcd"+Datasets[iDataset]+Form("%.1f", centRange[0])+"<cent<"+Form("%.1f", centRange[1]));
  // H3D_jetRjetPtcolCent->Sumw2();

  // int binCentEdges[2] = {H3D_jetRjetPtcolCent->GetZaxis()->FindBin(centRange[0]), H3D_jetRjetPtcolCent->GetZaxis()->FindBin(centRange[1])};
  // if (binCentEdges[0] == 0) 
  //   cout << "WARNING: Get_Pt_spectrum_bkgCorrected is counting the underflow with the chosen centRange" << endl;
  // if (binCentEdges[1] == H3D_jetRjetPtcolCent->GetZaxis()->GetNbins()+1)
  //   cout << "WARNING: Get_Pt_spectrum_bkgCorrected is counting the overflow with the chosen centRange" << endl;
  

  // int ibinJetRadius = H3D_jetRjetPtcolCent->GetXaxis()->FindBin(arrayRadius[iRadius]+GLOBAL_epsilon);

  // H1D_jetPt_mcdMatched_defaultBin = (TH1D*)H3D_jetRjetPtcolCent->ProjectionY("jetPt_mcdMatched_"+RadiusLegend[iRadius]+Form("%.1f", centRange[0])+"<cent<"+Form("%.1f", centRange[1]), ibinJetRadius, ibinJetRadius, binCentEdges[0], binCentEdges[1]);
  // // H1D_jetPt_mcdMatched_defaultBin = (TH1D*)H3D_jetRjetPtcolCent->ProjectionY("jetPt_mcdMatched_"+RadiusLegend[iRadius]+Form("%.1f", centRange[0])+"<cent<"+Form("%.1f", centRange[1]), ibinJetRadius, ibinJetRadius, 0, -1);

  // H1D_jetPt_mcdMatched = (TH1D*)H1D_jetPt_mcdMatched_defaultBin->Rebin(nBinPtJetsGen[iRadius],"jetPt_mcdMatched_rebinned_"+RadiusLegend[iRadius]+Form("%.1f", centRange[0])+"<cent<"+Form("%.1f", centRange[1]), ptBinsJetsGen[iRadius]);






  //////// without centrality; sadly we don't have a mcd version that gives centrality info, so gotta use one with all centralities (though one can edit the cent window in the jetfinderqa options)
  // H3D_jetRjetPtcolCent = (TH3D*)((TH3D*)file_O2Analysis_ppSimDetectorEffect->Get(analysisWorkflow+"/h3_jet_r_jet_pt_jet_eta"))->Clone("Get_Pt_spectrum_mcdMatched"+Datasets[iDataset]);
  // H3D_jetRjetPtcolCent->Sumw2();
  H3D_jetRjetPtcolCent = (TH3D*)((TH3D*)file_O2Analysis_ppSimDetectorEffect->Get(analysisWorkflowMC+"/h3_jet_r_jet_pt_tag_jet_pt_base_matchedgeo"))->Clone("Get_Pt_spectrum_mcdMatched"+Datasets[iDataset]);// base is mcd in jetfinderQA as of 06/2024, thus tag is mcp, and so hist is (x=r, y=mcp, z=mcd)
  H3D_jetRjetPtcolCent->Sumw2();


  int ibinJetRadius = H3D_jetRjetPtcolCent->GetXaxis()->FindBin(arrayRadius[iRadius]+GLOBAL_epsilon);

  H1D_jetPt_mcdMatched_defaultBin = (TH1D*)H3D_jetRjetPtcolCent->ProjectionZ("jetPt_mcdMatched_"+RadiusLegend[iRadius]+Datasets[iDataset], ibinJetRadius, ibinJetRadius, 0, -1);

  // should be rebinned using rec bins, because in the unfolding it'll be given to RooUnfoldResponse(mcdMatched, mcp, response) to get the efficiency 
  H1D_jetPt_mcdMatched = (TH1D*)H1D_jetPt_mcdMatched_defaultBin->Rebin(nBinPtJetsRec[iRadius],"jetPt_mcdMatched_rebinned_"+RadiusLegend[iRadius]+Datasets[iDataset], ptBinsJetsRec[iRadius]);







  if (strstr(options, "evtNorm") != NULL) {
    // NormaliseYieldToNEvents(H1D_jetPt_mcdMatched, GetNEventsGen(file_O2Analysis_list[iDataset]));
    // NormaliseYieldToNEvents(H1D_jetPt_mcdMatched, GetNEventsSelected_JetFramework_weighted(file_O2Analysis_ppSimDetectorEffect));
    NormaliseRawHistToNEvents(H1D_jetPt_mcdMatched, GetNEventsSelected_JetFramework_weighted(file_O2Analysis_ppSimDetectorEffect));
    // NormaliseYieldToNEntries(H1D_jetPt_mcdMatched);
  }
  if (strstr(options, "entriesNorm") != NULL) {
    NormaliseYieldToNEntries(H1D_jetPt_mcdMatched);
  }

}

void Get_Pt_spectrum_bkgCorrected(TH1D* &H1D_jetPt, int iDataset, int iRadius, float* centRange, const char options[]) {
  Get_Pt_spectrum_bkgCorrected_preWidthScaling(H1D_jetPt, iDataset, iRadius, centRange, options);
  TransformRawHistToYield(H1D_jetPt);
}

void Get_Pt_spectrum_mcp(TH1D* &H1D_jetPt_mcp, int iDataset, int iRadius, const char options[]) {
  Get_Pt_spectrum_mcp_preWidthScaling(H1D_jetPt_mcp, iDataset, iRadius, options);
  TransformRawHistToYield(H1D_jetPt_mcp);
}
void Get_Pt_spectrum_mcdMatched(TH1D* &H1D_jetPt_mcdMatched, int iDataset, int iRadius, const char options[]) {
  Get_Pt_spectrum_mcdMatched_preWidthScaling(H1D_jetPt_mcdMatched, iDataset, iRadius, options);
  TransformRawHistToYield(H1D_jetPt_mcdMatched);
}
void Get_Pt_spectrum_mcpMatched(TH1D* &H1D_jetPt_mcpMatched, int iDataset, int iRadius, const char options[]) {
  Get_Pt_spectrum_mcpMatched_preWidthScaling(H1D_jetPt_mcpMatched, iDataset, iRadius, options);
  TransformRawHistToYield(H1D_jetPt_mcpMatched);
}

void Get_Pt_spectrum_mcp_recBinning(TH1D* &H1D_jetPt_mcp, int iDataset, int iRadius, const char options[]) {
  Get_Pt_spectrum_mcp_recBinning_preWidthScaling(H1D_jetPt_mcp, iDataset, iRadius, options);
  TransformRawHistToYield(H1D_jetPt_mcp);
}

void Get_Pt_spectrum_bkgCorrected_genBinning(TH1D* &H1D_jetPt, int iDataset, int iRadius, float* centRange, const char options[]) {
  Get_Pt_spectrum_bkgCorrected_genBinning_preWidthScaling(H1D_jetPt, iDataset, iRadius, centRange, options);
  TransformRawHistToYield(H1D_jetPt);
}


void Get_Pt_spectrum_mcp_folded_preWidthScaling(TH1D* &H1D_jetPt_mcp_folded, int iDataset, int iRadius, float* centRange, const char options[]) {
  TH3D* H3D_jetRjetPtcolCent;
  TH1D* H1D_jetPt_mcp_defaultBin;
  TH1D* H1D_jetPt_mcp_fineBin;
  // TH1D* H1D_jetPt_raw[nRadius];
  TString partialUniqueSpecifier = Datasets[iDataset]+"_R="+Form("%.1f",arrayRadius[iRadius]);

  // //////// with centrality; sadly we don't have a mcp version that gives centrality info
  // H3D_jetRjetPtcolCent = (TH3D*)((TH3D*)file_O2Analysis_ppSimDetectorEffect[iCent]->Get(analysisWorkflow+"/h3_jet_r_jet_pt_centrality"))->Clone("Get_Pt_spectrum_mcp"+Datasets[iDataset]+Form("%.1f", centRange[0])+"<cent<"+Form("%.1f", centRange[1]));
  // H3D_jetRjetPtcolCent->Sumw2();

  // int binCentEdges[2] = {H3D_jetRjetPtcolCent->GetZaxis()->FindBin(centRange[0]), H3D_jetRjetPtcolCent->GetZaxis()->FindBin(centRange[1])};
  // if (binCentEdges[0] == 0) 
  //   cout << "WARNING: Get_Pt_spectrum_bkgCorrected is counting the underflow with the chosen centRange" << endl;
  // if (binCentEdges[1] == H3D_jetRjetPtcolCent->GetZaxis()->GetNbins()+1)
  //   cout << "WARNING: Get_Pt_spectrum_bkgCorrected is counting the overflow with the chosen centRange" << endl;
  

  // int ibinJetRadius = H3D_jetRjetPtcolCent->GetXaxis()->FindBin(arrayRadius[iRadius]+GLOBAL_epsilon);

  // H1D_jetPt_mcp_defaultBin = (TH1D*)H3D_jetRjetPtcolCent->ProjectionY("jetPt_mcp_"+RadiusLegend[iRadius]+Form("%.1f", centRange[0])+"<cent<"+Form("%.1f", centRange[1]), ibinJetRadius, ibinJetRadius, binCentEdges[0], binCentEdges[1]);
  // // H1D_jetPt_mcp_defaultBin = (TH1D*)H3D_jetRjetPtcolCent->ProjectionY("jetPt_mcp_"+RadiusLegend[iRadius]+Form("%.1f", centRange[0])+"<cent<"+Form("%.1f", centRange[1]), ibinJetRadius, ibinJetRadius, 0, -1);

  // H1D_jetPt_mcp = (TH1D*)H1D_jetPt_mcp_defaultBin->Rebin(nBinPtJetsGen[iRadius],"jetPt_mcp_rebinned_"+RadiusLegend[iRadius]+Form("%.1f", centRange[0])+"<cent<"+Form("%.1f", centRange[1]), ptBinsJetsGen[iRadius]);



  //////// without centrality; sadly we don't have a mcp version that gives centrality info, so gotta use one with all centralities (though one can edit the cent window in the jetfinderqa options)
  H3D_jetRjetPtcolCent = (TH3D*)((TH3D*)file_O2Analysis_ppSimDetectorEffect->Get(analysisWorkflowMC+"/h3_jet_r_part_jet_pt_part_jet_eta_part"))->Clone("Get_Pt_spectrum_mcp_folded_preWidthScaling_"+partialUniqueSpecifier);
  H3D_jetRjetPtcolCent->Sumw2();

  int ibinJetRadius = H3D_jetRjetPtcolCent->GetXaxis()->FindBin(arrayRadius[iRadius]+GLOBAL_epsilon);

  H1D_jetPt_mcp_defaultBin = (TH1D*)H3D_jetRjetPtcolCent->ProjectionY("H1D_jetPt_mcp_defaultBin_"+partialUniqueSpecifier, ibinJetRadius, ibinJetRadius, 0, -1);
  H1D_jetPt_mcp_fineBin = (TH1D*)H1D_jetPt_mcp_defaultBin->Rebin(nBinPtJetsFine[iRadius],"H1D_jetPt_mcp_fineBin_"+partialUniqueSpecifier, ptBinsJetsFine[iRadius]);

  TH2D* H2D_jetPtResponseMatrix_fluctuations;
  Get_PtResponseMatrix_Fluctuations(H2D_jetPtResponseMatrix_fluctuations, iDataset, iRadius, centRange);

  TH1D* H1D_jetPt_mcp_folded_fineBin = (TH1D*)H1D_jetPt_mcp_fineBin->Clone("H1D_jetPt_mcp_folded_fineBin_"+partialUniqueSpecifier); // OK because I have them both have a fine binning, same for both, before I multiply them
  H1D_jetPt_mcp_folded_fineBin->Reset("M");
  // matrix product of fluct response times truth; assumes the two are of the same size binning wise
  // Careful: xy of hist and ij of Resp(i,j) are inverted ! hist(j,i) = matrix(i,j) and so if vector(i)=SUM(matrixA(i,k)vector(k)) then hist(i)=SUM(histA(k,i)histB(k)), and if we replace j,i by gen,rec we get hist(iX)=SUM(histA(k,iX)histB(k))
  double matrixElementSum_iX, matrixElementSum_error_iX;
  for(int iBinX = 1; iBinX <= H2D_jetPtResponseMatrix_fluctuations->GetNbinsY(); iBinX++){ // 0 and n+1 take underflow and overflow into account
    matrixElementSum_iX = 0;
    matrixElementSum_error_iX = 0;
    for (int iMatrix = 1; iMatrix <= H2D_jetPtResponseMatrix_fluctuations->GetNbinsX(); iMatrix++){ // 0 and n+1 take underflow and overflow into account
      matrixElementSum_iX += H2D_jetPtResponseMatrix_fluctuations->GetBinContent(iMatrix, iBinX) * H1D_jetPt_mcp_fineBin->GetBinContent(iMatrix); 
      matrixElementSum_error_iX += pow(H2D_jetPtResponseMatrix_fluctuations->GetBinError(iMatrix, iBinX), 2) + pow(H1D_jetPt_mcp_fineBin->GetBinError(iMatrix), 2); // simple sigma_ab = a2sigma_a2 + b2sigma_b2 ; that assumes there are no correlations; here it s background fluct from PbPB sim, and detector effects from a pp sim, so we can maybe say theyre not correlated
    }
    H1D_jetPt_mcp_folded_fineBin->SetBinContent(iBinX, matrixElementSum_iX);
    H1D_jetPt_mcp_folded_fineBin->SetBinError(iBinX, sqrt(matrixElementSum_error_iX));
  }









  H1D_jetPt_mcp_folded = (TH1D*)H1D_jetPt_mcp_folded_fineBin->Rebin(nBinPtJetsRec[iRadius],"H1D_jetPt_mcp_folded_"+RadiusLegend[iRadius]+Datasets[iDataset], ptBinsJetsRec[iRadius]);



  if (strstr(options, "evtNorm") != NULL) {
    // NormaliseYieldToNEvents(H1D_jetPt_mcp_folded, GetNEventsGen(file_O2Analysis_list[iDataset]));
    // NormaliseYieldToNEvents(H1D_jetPt_mcp_folded, GetNEventsSelected_JetFramework_weighted(file_O2Analysis_ppSimDetectorEffect));
    NormaliseRawHistToNEvents(H1D_jetPt_mcp_folded, GetNEventsSelected_JetFramework_weighted(file_O2Analysis_ppSimDetectorEffect));
    
    // NormaliseYieldToNEntries(H1D_jetPt_mcp_folded);
    
  }
  if (strstr(options, "entriesNorm") != NULL) {
    NormaliseYieldToNEntries(H1D_jetPt_mcp_folded);
  }
}

void Get_Pt_spectrum_mcp_folded_genBinning_preWidthScaling(TH1D* &H1D_jetPt_mcp_folded, int iDataset, int iRadius, float* centRange, const char options[]) {
  TH3D* H3D_jetRjetPtcolCent;
  TH1D* H1D_jetPt_mcp_defaultBin;
  TH1D* H1D_jetPt_mcp_fineBin;
  // TH1D* H1D_jetPt_raw[nRadius];
  TString partialUniqueSpecifier = Datasets[iDataset]+"_R="+Form("%.1f",arrayRadius[iRadius]);

  // //////// with centrality; sadly we don't have a mcp version that gives centrality info
  // H3D_jetRjetPtcolCent = (TH3D*)((TH3D*)file_O2Analysis_ppSimDetectorEffect[iCent]->Get(analysisWorkflow+"/h3_jet_r_jet_pt_centrality"))->Clone("Get_Pt_spectrum_mcp"+Datasets[iDataset]+Form("%.1f", centRange[0])+"<cent<"+Form("%.1f", centRange[1]));
  // H3D_jetRjetPtcolCent->Sumw2();

  // int binCentEdges[2] = {H3D_jetRjetPtcolCent->GetZaxis()->FindBin(centRange[0]), H3D_jetRjetPtcolCent->GetZaxis()->FindBin(centRange[1])};
  // if (binCentEdges[0] == 0) 
  //   cout << "WARNING: Get_Pt_spectrum_bkgCorrected is counting the underflow with the chosen centRange" << endl;
  // if (binCentEdges[1] == H3D_jetRjetPtcolCent->GetZaxis()->GetNbins()+1)
  //   cout << "WARNING: Get_Pt_spectrum_bkgCorrected is counting the overflow with the chosen centRange" << endl;
  

  // int ibinJetRadius = H3D_jetRjetPtcolCent->GetXaxis()->FindBin(arrayRadius[iRadius]+GLOBAL_epsilon);

  // H1D_jetPt_mcp_defaultBin = (TH1D*)H3D_jetRjetPtcolCent->ProjectionY("jetPt_mcp_"+RadiusLegend[iRadius]+Form("%.1f", centRange[0])+"<cent<"+Form("%.1f", centRange[1]), ibinJetRadius, ibinJetRadius, binCentEdges[0], binCentEdges[1]);
  // // H1D_jetPt_mcp_defaultBin = (TH1D*)H3D_jetRjetPtcolCent->ProjectionY("jetPt_mcp_"+RadiusLegend[iRadius]+Form("%.1f", centRange[0])+"<cent<"+Form("%.1f", centRange[1]), ibinJetRadius, ibinJetRadius, 0, -1);

  // H1D_jetPt_mcp = (TH1D*)H1D_jetPt_mcp_defaultBin->Rebin(nBinPtJetsGen[iRadius],"jetPt_mcp_rebinned_"+RadiusLegend[iRadius]+Form("%.1f", centRange[0])+"<cent<"+Form("%.1f", centRange[1]), ptBinsJetsGen[iRadius]);



  //////// without centrality; sadly we don't have a mcp version that gives centrality info, so gotta use one with all centralities (though one can edit the cent window in the jetfinderqa options)
  H3D_jetRjetPtcolCent = (TH3D*)((TH3D*)file_O2Analysis_ppSimDetectorEffect->Get(analysisWorkflowMC+"/h3_jet_r_part_jet_pt_part_jet_eta_part"))->Clone("Get_Pt_spectrum_mcp_folded_preWidthScaling_"+partialUniqueSpecifier);
  H3D_jetRjetPtcolCent->Sumw2();

  int ibinJetRadius = H3D_jetRjetPtcolCent->GetXaxis()->FindBin(arrayRadius[iRadius]+GLOBAL_epsilon);

  H1D_jetPt_mcp_defaultBin = (TH1D*)H3D_jetRjetPtcolCent->ProjectionY("H1D_jetPt_mcp_defaultBin_"+partialUniqueSpecifier, ibinJetRadius, ibinJetRadius, 0, -1);
  H1D_jetPt_mcp_fineBin = (TH1D*)H1D_jetPt_mcp_defaultBin->Rebin(nBinPtJetsFine[iRadius],"H1D_jetPt_mcp_fineBin_"+partialUniqueSpecifier, ptBinsJetsFine[iRadius]);

  TH2D* H2D_jetPtResponseMatrix_fluctuations;
  Get_PtResponseMatrix_Fluctuations(H2D_jetPtResponseMatrix_fluctuations, iDataset, iRadius, centRange);

  TH1D* H1D_jetPt_mcp_folded_fineBin = (TH1D*)H1D_jetPt_mcp_fineBin->Clone("H1D_jetPt_mcp_folded_fineBin_"+partialUniqueSpecifier); // OK because I have them both have a fine binning, same for both, before I multiply them
  H1D_jetPt_mcp_folded_fineBin->Reset("M");
  // matrix product of fluct response times truth; assumes the two are of the same size binning wise
  // Careful: xy of hist and ij of Resp(i,j) are inverted ! hist(j,i) = matrix(i,j) and so if vector(i)=SUM(matrixA(i,k)vector(k)) then hist(i)=SUM(histA(k,i)histB(k)), and if we replace j,i by gen,rec we get hist(iX)=SUM(histA(k,iX)histB(k))
  double matrixElementSum_iX, matrixElementSum_error_iX;
  for(int iBinX = 1; iBinX <= H2D_jetPtResponseMatrix_fluctuations->GetNbinsY(); iBinX++){ // 0 and n+1 take underflow and overflow into account
    matrixElementSum_iX = 0;
    matrixElementSum_error_iX = 0;
    for (int iMatrix = 1; iMatrix <= H2D_jetPtResponseMatrix_fluctuations->GetNbinsX(); iMatrix++){ // 0 and n+1 take underflow and overflow into account
      matrixElementSum_iX += H2D_jetPtResponseMatrix_fluctuations->GetBinContent(iMatrix, iBinX) * H1D_jetPt_mcp_fineBin->GetBinContent(iMatrix); 
      matrixElementSum_error_iX += pow(H2D_jetPtResponseMatrix_fluctuations->GetBinError(iMatrix, iBinX), 2) + pow(H1D_jetPt_mcp_fineBin->GetBinError(iMatrix), 2); // simple sigma_ab = a2sigma_a2 + b2sigma_b2 ; that assumes there are no correlations; here it s background fluct from PbPB sim, and detector effects from a pp sim, so we can maybe say theyre not correlated
    }
    H1D_jetPt_mcp_folded_fineBin->SetBinContent(iBinX, matrixElementSum_iX);
    H1D_jetPt_mcp_folded_fineBin->SetBinError(iBinX, sqrt(matrixElementSum_error_iX));
  }









  H1D_jetPt_mcp_folded = (TH1D*)H1D_jetPt_mcp_folded_fineBin->Rebin(nBinPtJetsGen[iRadius],"H1D_jetPt_mcp_folded_gen_binning_"+RadiusLegend[iRadius]+Datasets[iDataset], ptBinsJetsGen[iRadius]);



  if (strstr(options, "evtNorm") != NULL) {
    // NormaliseYieldToNEvents(H1D_jetPt_mcp_folded, GetNEventsGen(file_O2Analysis_list[iDataset]));
    // NormaliseYieldToNEvents(H1D_jetPt_mcp_folded, GetNEventsSelected_JetFramework_weighted(file_O2Analysis_ppSimDetectorEffect));
    NormaliseRawHistToNEvents(H1D_jetPt_mcp_folded, GetNEventsSelected_JetFramework_weighted(file_O2Analysis_ppSimDetectorEffect));
    
    // NormaliseYieldToNEntries(H1D_jetPt_mcp_folded);
    
  }
  if (strstr(options, "entriesNorm") != NULL) {
    NormaliseYieldToNEntries(H1D_jetPt_mcp_folded);
  }
}

void Get_Pt_spectrum_mcp_folded(TH1D* &H1D_jetPt_mcp_folded, int iDataset, int iRadius, float* centRange, const char options[]) {
  Get_Pt_spectrum_mcp_folded_preWidthScaling(H1D_jetPt_mcp_folded, iDataset, iRadius, centRange, options);
  TransformRawHistToYield(H1D_jetPt_mcp_folded);
}

void Get_Pt_spectrum_mcp_folded_genBinning(TH1D* &H1D_jetPt, int iDataset, int iRadius, float* centRange, const char options[]) {
  Get_Pt_spectrum_mcp_folded_genBinning_preWidthScaling(H1D_jetPt, iDataset, iRadius, centRange, options);
  TransformRawHistToYield(H1D_jetPt);
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////// Response matrix functions ///////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


void Get_PtResponseMatrix_DetectorAndFluctuationsCombined_coarseBinningBeforeMult_obsolete(TH2D* &H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined, TH2D* H2D_jetPtResponseMatrix_detectorResponse, TH2D* H2D_jetPtResponseMatrix_fluctuations, int iDataset, int iRadius, float* centRange) { // to be created once the matching is ready
  // if (!(H2D_jetPtResponseMatrix_fluctuations->GetNbinsX() == H2D_jetPtResponseMatrix_detectorResponse->GetNbinsY())){
  //   cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl;
  //   cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl;
  //   cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl;
  //   cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl;
  //   cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl;
  //   cout << "!!!!!!!!!!!!!!!!!!!!!!                                                                                 !!!!!!!!!!!!!!!!!!!!!!!" << endl;
  //   cout << "!!!!!!!!!!!!!!!!!!!!!!                                                                                 !!!!!!!!!!!!!!!!!!!!!!!" << endl;
  //   cout << "!!!!!!!!!!!!!!!!!!!!!!                                                                                 !!!!!!!!!!!!!!!!!!!!!!!" << endl;
  //   cout << "!!!!!!!!!!!!!!!!!!!!!!                         matrix multiplication impossible                        !!!!!!!!!!!!!!!!!!!!!!!" << endl;
  //   cout << "!!!!!!!!!!!!!!!!!!!!!!                     dimensions of the matrices incompatible                     !!!!!!!!!!!!!!!!!!!!!!!" << endl;
  //   cout << "!!!!!!!!!!!!!!!!!!!!!!                                                                                 !!!!!!!!!!!!!!!!!!!!!!!" << endl;
  //   cout << "!!!!!!!!!!!!!!!!!!!!!!                                                                                 !!!!!!!!!!!!!!!!!!!!!!!" << endl;
  //   cout << "!!!!!!!!!!!!!!!!!!!!!!                                                                                 !!!!!!!!!!!!!!!!!!!!!!!" << endl;
  //   cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl;
  //   cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl;
  //   cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl;
  //   cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl;
  //   cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl;
  // }


  TString partialUniqueSpecifier = Datasets[iDataset]+"_R="+Form("%.1f",arrayRadius[iRadius])+"_@cent["+Form("%.1f", centRange[0])+","+Form("%.1f", centRange[1])+"]";

  // TH2D *H2D_PtResponseMatrix = new TH2D("H2D_PtResponseMatrix", "H2D_PtResponseMatrix", nbinpT, pTbins, nbinpT, pTbins);
  // si je fais une multiplication de fluct(mcpBins, mcdBins) x det(mcpBins, mcdBins) il faut que mcpBins == mcdBins ou la multiplication de matrice fonctionne pas
  H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined = (TH2D*)H2D_jetPtResponseMatrix_fluctuations->Clone("Get_PtResponseMatrix_DetectorAndFluctuationsCombined"+partialUniqueSpecifier);
  H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined->Reset("M");

  // matrix product of fluct response times det response; assumes the two are of the same size binning wise
  double matrixElementSum_iRec_iGen, matrixElementSum_error_iRec_iGen;
  for(int iBinRec = 1; iBinRec <= H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined->GetNbinsX(); iBinRec++){ // 0 and n+1 take underflow and overflow into account
    for(int iBinGen = 1; iBinGen <= H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined->GetNbinsY(); iBinGen++){ // 0 and n+1 take underflow and overflow into account
      matrixElementSum_iRec_iGen = 0;
      matrixElementSum_error_iRec_iGen = 0;
      // matrix multiplication: detector response matrix times background fluctuation response matrix
      for (int iMatrix = 0; iMatrix <= H2D_jetPtResponseMatrix_fluctuations->GetNbinsX()+1; iMatrix++){ // 0 and n+1 take underflow and overflow into account
        matrixElementSum_iRec_iGen += H2D_jetPtResponseMatrix_fluctuations->GetBinContent(iBinRec, iMatrix) * H2D_jetPtResponseMatrix_detectorResponse->GetBinContent(iMatrix, iBinGen);
        matrixElementSum_error_iRec_iGen += pow(H2D_jetPtResponseMatrix_fluctuations->GetBinError(iBinRec, iMatrix), 2) + pow(H2D_jetPtResponseMatrix_detectorResponse->GetBinError(iMatrix, iBinGen), 2); // simple sigma_ab = a2sigma_a2 + b2sigma_b2 ; that assumes there are no correlations; here it s background fluct from PbPB sim, and detector effects from a pp sim, so we can maybe say theyre not correlated
      }
      // cout << "iBinRec, iBinGen, matrixElementSum_iRec_iGen = " << iBinRec << ", " << iBinGen << ", " << matrixElementSum_iRec_iGen << endl;
      H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined->SetBinContent(iBinRec, iBinGen, matrixElementSum_iRec_iGen);
      H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined->SetBinError(iBinRec, iBinGen, sqrt(matrixElementSum_error_iRec_iGen));
    }
  }


  if (doNormalisation) {
    // Normalisation of the combined response matrix: y-axis is weighted by the prior, the y-projection of the detector response matrix, the generator-level jet spectrum
    TH1D* priorSpectrum = (TH1D*)H2D_jetPtResponseMatrix_detectorResponse->ProjectionY("priorSpectrum"+partialUniqueSpecifier, 1, H2D_jetPtResponseMatrix_detectorResponse->GetNbinsX());
    cout << "prior weighting of combined response matrix: should I keep under/overflows in the projection? Should I normalise the over/underflow bins as well?" << endl;
    cout << "prior and combined resp are correlated -> I should put that into the error calculation" << endl;
    double combinedResponseContent, priorSpectrumContent;
    double combinedResponseError, priorSpectrumError;
    for(int iBinY = 1; iBinY <= H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined->GetNbinsY(); iBinY++){
      priorSpectrumContent = priorSpectrum->GetBinContent(iBinY);
      priorSpectrumError = priorSpectrum->GetBinError(iBinY);
      cout << "is ibinY = " << iBinY << " overflow? " << priorSpectrum->IsBinOverflow(iBinY) << endl;
      for(int iBinX = 1; iBinX <= H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined->GetNbinsX(); iBinX++){
        combinedResponseContent = H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined->GetBinContent(iBinX, iBinY);
        combinedResponseError = H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined->GetBinError(iBinX, iBinY);
        H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined->SetBinContent(iBinX, iBinY, combinedResponseContent * 1./priorSpectrumContent);
        H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined->SetBinError(iBinX, iBinY, sqrt( combinedResponseError*combinedResponseError / (combinedResponseContent*combinedResponseContent) + priorSpectrumError*priorSpectrumError / (priorSpectrumContent*priorSpectrumContent) ));

        // cout << "is H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined(" << iBinX << ", " << iBinY << ") overflow? " << H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined->IsBinOverflow(iBinX, iBinY) << endl;
      }
    }
  }
}

void Get_PtResponseMatrix_detectorResponse_coarseBinningBeforeMult_obsolete(TH2D* &H2D_jetPtResponseMatrix_detectorResponse, int iDataset, int iRadius) {
  TH3D* H3D_jetRpartPtdetPt;
  TH2D* H2D_gen_det_geoMatched;
  TH2D* H2D_gen_det_geoMatched_rebinned;

  TString partialUniqueSpecifier = Datasets[iDataset]+"_R="+Form("%.1f",arrayRadius[iRadius]);

  H3D_jetRpartPtdetPt = (TH3D*)((TH3D*)file_O2Analysis_ppSimDetectorEffect->Get(analysisWorkflowMC+"/h3_jet_r_jet_pt_tag_jet_pt_base_matchedgeo"))->Clone("Get_PtResponseMatrix_detectorResponse"+partialUniqueSpecifier);// base is mcd in jetfinderQA as of 06/2024, thus tag is mcp, and so hist is (x=r, y=mcp, z=mcd)
  H3D_jetRpartPtdetPt->Sumw2();

  int ibinJetRadius = H3D_jetRpartPtdetPt->GetXaxis()->FindBin(arrayRadius[iRadius]+GLOBAL_epsilon);
  H3D_jetRpartPtdetPt->GetXaxis()->SetRange(ibinJetRadius, ibinJetRadius);

  // project H3D onto a H2D, option "yz" means y goes on y-axis while z goes on x-axis, and so H2D_gen_det_geoMatched will be (x=mcd, y=mcp)
  H2D_gen_det_geoMatched = (TH2D*)H3D_jetRpartPtdetPt->Project3D(partialUniqueSpecifier+"_genrec_e_yz"); //can't use letter D in this or it seems to replace the histogram in current pad (see documentation of ProjectionX function. Isn't mentioned in project3D sadly)

  // keep (gen, gen) for the bins; rec will be introduced in the fluctuation response, and by multiplication will stay in the combined matrix
  TH2D* H2D_response = (TH2D*)RebinVariableBins2D(H2D_gen_det_geoMatched, nBinPtJetsGen[iRadius], nBinPtJetsGen[iRadius], ptBinsJetsGen[iRadius], ptBinsJetsGen[iRadius]).Clone("Get_PtResponseMatrix_detectorResponse_rebinned"+partialUniqueSpecifier);


  // Normalisation of the response matrix: each pt gen slice is normalised to unity
  double genSliceNorm = 1;
  for(int iBinY = 1; iBinY <= H2D_response->GetNbinsY(); iBinY++){
    genSliceNorm = H2D_response->Integral(1, H2D_response->GetNbinsX(), iBinY, iBinY);
    for(int iBinX = 1; iBinX <= H2D_response->GetNbinsX(); iBinX++){
      H2D_response->SetBinContent(iBinX, iBinY, H2D_response->GetBinContent(iBinX, iBinY) * 1./genSliceNorm);
      H2D_response->SetBinError(iBinX, iBinY, H2D_response->GetBinError(iBinX, iBinY) * 1./genSliceNorm);
      // cout << "H2D_response(" << iBinX << ", " << iBinY << ") = " << H2D_response->GetBinContent(iBinX, iBinY) << endl;
    }
  }

  // H2D_jetPtResponseMatrix_detectorResponse = (TH2D*)H2D_gen_det_geoMatched_rebinned->Clone("H2D_jetPtResponseMatrix_detectorResponse"+partialUniqueSpecifier); // should be using the one below; this one is just a test
  H2D_jetPtResponseMatrix_detectorResponse = (TH2D*)H2D_response->Clone("H2D_jetPtResponseMatrix_detectorResponse"+partialUniqueSpecifier); 

  cout << "############ detectorResponse matrix ########################" << endl;
  cout << "############ H3D_jetRpartPtdetPt:" << endl;
  cout << "############ nBins_x = " << H3D_jetRpartPtdetPt->GetNbinsX() << ", nBins_y = " << H3D_jetRpartPtdetPt->GetNbinsY() << ", nBins_z = " << H3D_jetRpartPtdetPt->GetNbinsZ() << endl;
  cout << "############ nEntries = " << H3D_jetRpartPtdetPt->GetEntries() << endl;
  cout << "############ nEffectiveEntries = " << H3D_jetRpartPtdetPt->GetEffectiveEntries() << endl;
  cout << "############ H2D_gen_det_geoMatched:" << endl;
  cout << "############ nBins_x = " << H2D_gen_det_geoMatched->GetNbinsX() << ", nBins_y = " << H2D_gen_det_geoMatched->GetNbinsY() << endl;
  cout << "############ nEntries = " << H2D_gen_det_geoMatched->GetEntries() << endl;
  cout << "############ H2D_response:" << endl;
  cout << "############ nBins_x = " << H2D_response->GetNbinsX() << ", nBins_y = " << H2D_response->GetNbinsY() << endl;
  cout << "############ nEntries = " << H2D_response->GetEntries() << endl;
  cout << "############ H2D_jetPtResponseMatrix_detectorResponse:" << endl;
  cout << "############ nBins_x = " << H2D_jetPtResponseMatrix_detectorResponse->GetNbinsX() << ", nBins_y = " << H2D_jetPtResponseMatrix_detectorResponse->GetNbinsY() << endl;
  cout << "############ nEntries = " << H2D_jetPtResponseMatrix_detectorResponse->GetEntries() << endl;
  cout << "############ ----test---- H2D_gen_det_geoMatched->Integral(100, 200, 100, 200) = " << H2D_gen_det_geoMatched->Integral(100, 200, 100, 200) << endl;
  cout << "############ ----test---- H3D_jetRpartPtdetPt->Integral(1, 1, 100, 200, 100, 200) = " << H3D_jetRpartPtdetPt->Integral(1, 1, 100, 200, 100, 200) << endl;

  // its better but the matrix is pretty empty at high pt; why???
  // gotta "normalise the projection on the pgen T -axis to unity" for each pt gen
}

void Get_PtResponseMatrix_Fluctuations_coarseBinningBeforeMult_obsolete(TH2D* &H2D_jetPtResponseMatrix_fluctuations, int iDataset, int iRadius, float* centRange) { 
  // see Hiroki Yokoyama thesis
  // iRadius is for chosing the pT binning

  TString partialUniqueSpecifier = Datasets[iDataset]+"_R="+Form("%.1f",arrayRadius[iRadius])+"_@cent["+Form("%.1f", centRange[0])+","+Form("%.1f", centRange[1])+"]";

  TH2D* H2D_fluctuations_centrality;
  TH1D* H1D_fluctuations;

  H2D_fluctuations_centrality = (TH2D*)((TH2D*)file_O2Analysis_list[iDataset]->Get(analysisWorkflowData+"/h2_centrality_rhorandomconewithoutleadingjet"))->Clone("Get_PtResponseMatrix_Fluctuations"+Datasets[iDataset]);
  H2D_fluctuations_centrality->Sumw2();

  int ibinCent_low, ibinCent_high;

  ibinCent_low = H2D_fluctuations_centrality->GetXaxis()->FindBin(centRange[0]);
  ibinCent_high = H2D_fluctuations_centrality->GetXaxis()->FindBin(centRange[1]);
  H1D_fluctuations = (TH1D*)H2D_fluctuations_centrality->ProjectionY("bkgFluctuationCentrality_highRes_"+partialUniqueSpecifier, ibinCent_low, ibinCent_high);

  NormaliseYieldToNEntries(H1D_fluctuations); // normalising fluctuations to 1
    
  TH2D H2D_response = TH2D("H2D_response_"+partialUniqueSpecifier, "H2D_response_"+partialUniqueSpecifier, nBinPtJetsRec[iRadius], ptBinsJetsRec[iRadius], nBinPtJetsGen[iRadius], ptBinsJetsGen[iRadius]); // actually doesn't work if original histogram has fixed bin size

  //==================== Build response matrix: shift deltaPt by pT gen along the pT rec axis ====================//
  int ibinZeroFluct= H1D_fluctuations->FindBin(0+GLOBAL_epsilon);
  
  for(int iBinRec = 1; iBinRec <= H2D_response.GetNbinsX(); iBinRec++){
    for(int iBinGen = 1; iBinGen <= H2D_response.GetNbinsY(); iBinGen++){
      double ptGen = H2D_response.GetYaxis()->GetBinCenter(iBinGen);
      double ptRec_low = H2D_response.GetXaxis()->GetBinLowEdge(iBinRec);
      double ptRec_up = H2D_response.GetXaxis()->GetBinLowEdge(iBinRec+1);
      // double xPtRecWidth = H2D_response->GetXaxis()->GetBinWidth(iBinRec);
      // if (ibinZero + (iBinRec - iBinGen) <= H1D_fluctuations_highRes->GetNbinsX()) { // make sure it's within bin range
      H2D_response.SetBinContent(iBinRec, iBinGen, H1D_fluctuations->Integral(H1D_fluctuations->GetXaxis()->FindBin(ptRec_low - ptGen + GLOBAL_epsilon), H1D_fluctuations->GetXaxis()->FindBin(ptRec_up - ptGen + GLOBAL_epsilon))); 
      // cout << "FluctResp(" << iBinRec << ", " << iBinGen << ") = " << H2D_response.GetBinContent(iBinRec, iBinGen) << endl;
        // H2D_response_highRes->SetBinError(iBinRec, iBinGen, something); 
    }
  }
  // }
  //========================================= Build response matrix end =========================================//

  H2D_jetPtResponseMatrix_fluctuations = (TH2D*)H2D_response.Clone("H2D_jetPtResponseMatrix_fluctuations"+partialUniqueSpecifier); // is replaced below with pT binning from analysis; rebin 2D with variable bins doesn't work, need to use custom function RebinVariableBins2D
  // TH2D H2D_jetPtResponseMatrix_fluctuations_rebinned = RebinVariableBins2D(H2D_jetPtResponseMatrix_fluctuations, nBinPtJetsRec[iRadius], nBinPtJetsGen[iRadius], ptBinsJetsRec[iRadius], ptBinsJetsGen[iRadius]);

  // cout << "H2D_jetPtResponseMatrix_fluctuations nbins X =" << H2D_jetPtResponseMatrix_fluctuations->GetNbinsX() << endl;
  // cout << "H2D_jetPtResponseMatrix_fluctuations_rebinned nbins X =" << H2D_jetPtResponseMatrix_fluctuations_rebinned.GetNbinsX() << endl;
  // H2D_jetPtResponseMatrix_fluctuations = (TH2D*)RebinVariableBins2D(H2D_jetPtResponseMatrix_fluctuations, nBinPtJetsRec[iRadius], nBinPtJetsGen[iRadius], ptBinsJetsRec[iRadius], ptBinsJetsGen[iRadius]).Clone("H2D_jetPtResponseMatrix_fluctuations_rebinned_"+partialUniqueSpecifier);
}


void Get_PtResponseMatrix_DetectorAndFluctuationsCombined(TH2D* &H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined, TH2D* H2D_jetPtResponseMatrix_detectorResponse, TH2D* H2D_jetPtResponseMatrix_fluctuations, int iDataset, int iRadius, float* centRange) { // to be created once the matching is ready


  TString partialUniqueSpecifier = Datasets[iDataset]+"_R="+Form("%.1f",arrayRadius[iRadius])+"_@cent["+Form("%.1f", centRange[0])+","+Form("%.1f", centRange[1])+"]";

  // TH2D *H2D_PtResponseMatrix = new TH2D("H2D_PtResponseMatrix", "H2D_PtResponseMatrix", nbinpT, pTbins, nbinpT, pTbins);
  // si je fais une multiplication de fluct(mcpBins, mcdBins) x det(mcpBins, mcdBins) il faut que mcpBins == mcdBins ou la multiplication de matrice fonctionne pas
  TH2D* H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined_preRebin = (TH2D*)H2D_jetPtResponseMatrix_fluctuations->Clone("Get_PtResponseMatrix_DetectorAndFluctuationsCombined_preRebin"+partialUniqueSpecifier); // OK because I have them both have a fine binning, same for both, before I multiply them
  H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined_preRebin->Reset("M");

  // // matrix product of fluct response times det response; assumes the two are of the same size binning wise
  // double matrixElementSum_iRec_iGen, matrixElementSum_error_iRec_iGen;
  // for(int iBinRec = 0; iBinRec <= H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined_preRebin->GetNbinsX()+1; iBinRec++){ // 0 and n+1 take underflow and overflow into account
  //   for(int iBinGen = 0; iBinGen <= H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined_preRebin->GetNbinsY()+1; iBinGen++){ // 0 and n+1 take underflow and overflow into account
  //     matrixElementSum_iRec_iGen = 0;
  //     matrixElementSum_error_iRec_iGen = 0;
  //     // matrix multiplication: detector response matrix times background fluctuation response matrix
  //     for (int iMatrix = 0; iMatrix <= H2D_jetPtResponseMatrix_fluctuations->GetNbinsX()+1; iMatrix++){ // 0 and n+1 take underflow and overflow into account
  //       matrixElementSum_iRec_iGen += H2D_jetPtResponseMatrix_fluctuations->GetBinContent(iBinRec, iMatrix) * H2D_jetPtResponseMatrix_detectorResponse->GetBinContent(iMatrix, iBinGen);
  //       matrixElementSum_error_iRec_iGen += pow(H2D_jetPtResponseMatrix_fluctuations->GetBinError(iBinRec, iMatrix), 2) + pow(H2D_jetPtResponseMatrix_detectorResponse->GetBinError(iMatrix, iBinGen), 2); // simple sigma_ab = a2sigma_a2 + b2sigma_b2 ; that assumes there are no correlations; here it s background fluct from PbPB sim, and detector effects from a pp sim, so we can maybe say theyre not correlated
  //     }
  //     // cout << "iBinRec, iBinGen, matrixElementSum_iRec_iGen = " << iBinRec << ", " << iBinGen << ", " << matrixElementSum_iRec_iGen << endl;
  //     H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined_preRebin->SetBinContent(iBinGen, iBinRec, matrixElementSum_iRec_iGen);
  //     H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined_preRebin->SetBinError(iBinGen, iBinRec, sqrt(matrixElementSum_error_iRec_iGen));
  //   }
  // }

  // matrix product of fluct response times det response; assumes the two are of the same size binning wise
  // Careful: xy of hist and ij of Resp(i,j) are inverted ! hist(j,i) = matrix(i,j) and so if matrix(i,j)=SUM(matrixA(i,k)matrixB(k,j)) then hist(j,i)=SUM(histA(k,i)histB(j,k)), and if we replace j,i by gen,rec we get hist(gen,rec)=SUM(histA(k,rec)histB(gen,k))
  double matrixElementSum_iGen_iRec, matrixElementSum_error_iGen_iRec;
  for(int iBinGen = 1; iBinGen <= H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined_preRebin->GetNbinsX(); iBinGen++){ // 0 and n+1 take underflow and overflow into account
    for(int iBinRec = 1; iBinRec <= H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined_preRebin->GetNbinsY(); iBinRec++){ // 0 and n+1 take underflow and overflow into account
      matrixElementSum_iGen_iRec = 0;
      matrixElementSum_error_iGen_iRec = 0;
      for (int iMatrix = 1; iMatrix <= H2D_jetPtResponseMatrix_fluctuations->GetNbinsX(); iMatrix++){ // 0 and n+1 take underflow and overflow into account
        matrixElementSum_iGen_iRec += H2D_jetPtResponseMatrix_fluctuations->GetBinContent(iMatrix, iBinRec) * H2D_jetPtResponseMatrix_detectorResponse->GetBinContent(iBinGen, iMatrix); 
        matrixElementSum_error_iGen_iRec += pow(H2D_jetPtResponseMatrix_fluctuations->GetBinError(iMatrix, iBinRec), 2) + pow(H2D_jetPtResponseMatrix_detectorResponse->GetBinError(iBinGen, iMatrix), 2); // simple sigma_ab = a2sigma_a2 + b2sigma_b2 ; that assumes there are no correlations; here it s background fluct from PbPB sim, and detector effects from a pp sim, so we can maybe say theyre not correlated
      }
      H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined_preRebin->SetBinContent(iBinGen, iBinRec, matrixElementSum_iGen_iRec);
      H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined_preRebin->SetBinError(iBinGen, iBinRec, sqrt(matrixElementSum_error_iGen_iRec));
    }
  }

  H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined = (TH2D*)RebinVariableBins2D(H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined_preRebin, nBinPtJetsGen[iRadius], nBinPtJetsRec[iRadius], ptBinsJetsGen[iRadius], ptBinsJetsRec[iRadius]).Clone("Get_PtResponseMatrix_DetectorAndFluctuationsCombined"+partialUniqueSpecifier);


  if (doNormalisation) { // should I not take the bin width into account? else some bins after normalisation will still have values above 1? eh why though, shouldnt be true

    //  Hiroki's version

    // Normalisation of the combined response matrix: y-axis is weighted by the prior, the generator-level jet spectrum, at step 3 of detResp construction in hiroki's thesis (ie before matching to fill response and normalisation of pt-gen axis to unity)
    
    // old // TH2D* H2D_jetPtResponseMatrix_detectorResponse_rebinned = (TH2D*)RebinVariableBins2D(H2D_jetPtResponseMatrix_detectorResponse, nBinPtJetsGen[iRadius], nBinPtJetsGen[iRadius], ptBinsJetsGen[iRadius], ptBinsJetsGen[iRadius]).Clone("Get_PtResponseMatrix_DetectorAndFluctuationsCombined"+partialUniqueSpecifier);
    // old // TH1D* priorSpectrum = (TH1D*)H2D_jetPtResponseMatrix_detectorResponse_rebinned->ProjectionY("priorSpectrum"+partialUniqueSpecifier, 1, H2D_jetPtResponseMatrix_detectorResponse_rebinned->GetNbinsX());

    TH1D* priorSpectrum;
    Get_Pt_spectrum_mcp(priorSpectrum, iDataset, iRadius, "evtNorm"); // _recBinning if needed ; 
    cout << "maybe need to have priorSpectrum normalised to number of events? or maybe I shouldn't?" << endl;

    cout << "prior weighting of combined response matrix: should I keep under/overflows in the projection? Should I normalise the over/underflow bins as well?" << endl;
    cout << "prior and combined resp are correlated -> I should put that into the error calculation" << endl;
    double combinedResponseContent, priorSpectrumContent;
    double combinedResponseError, priorSpectrumError;
    for(int iBinX = 1; iBinX <= H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined->GetNbinsX(); iBinX++){
      priorSpectrumContent = priorSpectrum->GetBinContent(iBinX);
      priorSpectrumError = priorSpectrum->GetBinError(iBinX);
      cout << "is iBinX = " << iBinX << " overflow? " << priorSpectrum->IsBinOverflow(iBinX) << endl;
      for(int iBinY = 1; iBinY <= H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined->GetNbinsY(); iBinY++){
        combinedResponseContent = H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined->GetBinContent(iBinX, iBinY);
        combinedResponseError = H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined->GetBinError(iBinX, iBinY);
        H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined->SetBinContent(iBinX, iBinY, combinedResponseContent * 1./priorSpectrumContent);
        H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined->SetBinError(iBinX, iBinY, sqrt( combinedResponseError*combinedResponseError / (combinedResponseContent*combinedResponseContent) + priorSpectrumError*priorSpectrumError / (priorSpectrumContent*priorSpectrumContent) ));

        // cout << "is H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined(" << iBinX << ", " << iBinY << ") overflow? " << H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined->IsBinOverflow(iBinX, iBinY) << endl;
      }
    }
  }
  
}

void Get_PtResponseMatrix_detectorResponse(TH2D* &H2D_jetPtResponseMatrix_detectorResponse, int iDataset, int iRadius) {
  TH3D* H3D_jetRpartPtdetPt;
  TH2D* H2D_gen_det_geoMatched;
  TH2D* H2D_gen_det_geoMatched_rebinned;

  TString partialUniqueSpecifier = Datasets[iDataset]+"_R="+Form("%.1f",arrayRadius[iRadius]);

  H3D_jetRpartPtdetPt = (TH3D*)((TH3D*)file_O2Analysis_ppSimDetectorEffect->Get(analysisWorkflowMC+"/h3_jet_r_jet_pt_tag_jet_pt_base_matchedgeo"))->Clone("Get_PtResponseMatrix_detectorResponse"+partialUniqueSpecifier);// base is mcd in jetfinderQA as of 06/2024, thus tag is mcp, and so hist is (x=r, y=mcp, z=mcd)
  H3D_jetRpartPtdetPt->Sumw2();

  int ibinJetRadius = H3D_jetRpartPtdetPt->GetXaxis()->FindBin(arrayRadius[iRadius]+GLOBAL_epsilon);
  H3D_jetRpartPtdetPt->GetXaxis()->SetRange(ibinJetRadius, ibinJetRadius);

  // project H3D onto a H2D, option "zy" means z goes on y-axis while y goes on x-axis, and so H2D_gen_det_geoMatched will be (x=mcp, y=mcd)
  H2D_gen_det_geoMatched = (TH2D*)H3D_jetRpartPtdetPt->Project3D(partialUniqueSpecifier+"_genrec_e_zy"); //can't use letter D in this or it seems to replace the histogram in current pad (see documentation of ProjectionX function. Isn't mentioned in project3D sadly)

  // keep (gen, gen) for the bins; rec will be introduced in the fluctuation response, and by multiplication will stay in the combined matrix
  TH2D* H2D_response = (TH2D*)RebinVariableBins2D(H2D_gen_det_geoMatched, nBinPtJetsFine[iRadius], nBinPtJetsFine[iRadius], ptBinsJetsFine[iRadius], ptBinsJetsFine[iRadius]).Clone("Get_PtResponseMatrix_detectorResponse_rebinned"+partialUniqueSpecifier);


  // Normalisation of the response matrix: each pt gen slice is normalised to unity
  double genSliceNorm = 1;
  for(int iBinX = 1; iBinX <= H2D_response->GetNbinsX(); iBinX++){
    genSliceNorm = H2D_response->Integral(1, H2D_response->GetNbinsY(), iBinX, iBinX);
    for(int iBinY = 1; iBinY <= H2D_response->GetNbinsY(); iBinY++){
      H2D_response->SetBinContent(iBinX, iBinY, H2D_response->GetBinContent(iBinX, iBinY) * 1./genSliceNorm);
      H2D_response->SetBinError(iBinX, iBinY, H2D_response->GetBinError(iBinX, iBinY) * 1./genSliceNorm);
      // cout << "H2D_response(" << iBinX << ", " << iBinY << ") = " << H2D_response->GetBinContent(iBinX, iBinY) << endl;
    }
  }

  // H2D_jetPtResponseMatrix_detectorResponse = (TH2D*)H2D_gen_det_geoMatched_rebinned->Clone("H2D_jetPtResponseMatrix_detectorResponse"+partialUniqueSpecifier); // should be using the one below; this one is just a test
  H2D_jetPtResponseMatrix_detectorResponse = (TH2D*)H2D_response->Clone("H2D_jetPtResponseMatrix_detectorResponse"+partialUniqueSpecifier); 

  cout << "############ detectorResponse matrix ########################" << endl;
  cout << "############ H3D_jetRpartPtdetPt:" << endl;
  cout << "############ nBins_x = " << H3D_jetRpartPtdetPt->GetNbinsX() << ", nBins_y = " << H3D_jetRpartPtdetPt->GetNbinsY() << ", nBins_z = " << H3D_jetRpartPtdetPt->GetNbinsZ() << endl;
  cout << "############ nEntries = " << H3D_jetRpartPtdetPt->GetEntries() << endl;
  cout << "############ nEffectiveEntries = " << H3D_jetRpartPtdetPt->GetEffectiveEntries() << endl;
  cout << "############ H2D_gen_det_geoMatched:" << endl;
  cout << "############ nBins_x = " << H2D_gen_det_geoMatched->GetNbinsX() << ", nBins_y = " << H2D_gen_det_geoMatched->GetNbinsY() << endl;
  cout << "############ nEntries = " << H2D_gen_det_geoMatched->GetEntries() << endl;
  cout << "############ H2D_response:" << endl;
  cout << "############ nBins_x = " << H2D_response->GetNbinsX() << ", nBins_y = " << H2D_response->GetNbinsY() << endl;
  cout << "############ nEntries = " << H2D_response->GetEntries() << endl;
  cout << "############ H2D_jetPtResponseMatrix_detectorResponse:" << endl;
  cout << "############ nBins_x = " << H2D_jetPtResponseMatrix_detectorResponse->GetNbinsX() << ", nBins_y = " << H2D_jetPtResponseMatrix_detectorResponse->GetNbinsY() << endl;
  cout << "############ nEntries = " << H2D_jetPtResponseMatrix_detectorResponse->GetEntries() << endl;
  cout << "############ ----test---- H2D_gen_det_geoMatched->Integral(100, 200, 100, 200) = " << H2D_gen_det_geoMatched->Integral(100, 200, 100, 200) << endl;
  cout << "############ ----test---- H3D_jetRpartPtdetPt->Integral(1, 1, 100, 200, 100, 200) = " << H3D_jetRpartPtdetPt->Integral(1, 1, 100, 200, 100, 200) << endl;

  // its better but the matrix is pretty empty at high pt; why???
  // gotta "normalise the projection on the pgen T -axis to unity" for each pt gen
}

void Get_PtResponseMatrix_Fluctuations(TH2D* &H2D_jetPtResponseMatrix_fluctuations, int iDataset, int iRadius, float* centRange) { 
  // see Hiroki Yokoyama thesis
  // iRadius is for chosing the pT binning

  TString partialUniqueSpecifier = Datasets[iDataset]+"_R="+Form("%.1f",arrayRadius[iRadius])+"_@cent["+Form("%.1f", centRange[0])+","+Form("%.1f", centRange[1])+"]";

  TH2D* H2D_fluctuations_centrality;
  TH1D* H1D_fluctuations;

  H2D_fluctuations_centrality = (TH2D*)((TH2D*)file_O2Analysis_list[iDataset]->Get(analysisWorkflowData+"/h2_centrality_rhorandomconewithoutleadingjet"))->Clone("Get_PtResponseMatrix_Fluctuations"+Datasets[iDataset]);
  H2D_fluctuations_centrality->Sumw2();

  int ibinCent_low, ibinCent_high;

  ibinCent_low = H2D_fluctuations_centrality->GetXaxis()->FindBin(centRange[0]);
  ibinCent_high = H2D_fluctuations_centrality->GetXaxis()->FindBin(centRange[1]);
  H1D_fluctuations = (TH1D*)H2D_fluctuations_centrality->ProjectionY("bkgFluctuationCentrality_highRes_"+partialUniqueSpecifier, ibinCent_low, ibinCent_high);

  NormaliseYieldToNEntries(H1D_fluctuations); // normalising fluctuations to 1
    
  TH2D H2D_response = TH2D("H2D_response_"+partialUniqueSpecifier, "H2D_response_"+partialUniqueSpecifier, nBinPtJetsFine[iRadius], ptBinsJetsFine[iRadius], nBinPtJetsFine[iRadius], ptBinsJetsFine[iRadius]); // actually doesn't work if original histogram has fixed bin size

  //==================== Build response matrix: shift deltaPt by pT gen along the pT rec axis ====================//
  int ibinZeroFluct= H1D_fluctuations->FindBin(0+GLOBAL_epsilon);
  
  for(int iBinGen = 1; iBinGen <= H2D_response.GetNbinsX(); iBinGen++){
    for(int iBinRec = 1; iBinRec <= H2D_response.GetNbinsY(); iBinRec++){
      double ptGen = H2D_response.GetXaxis()->GetBinCenter(iBinGen);
      double ptRec_low = H2D_response.GetYaxis()->GetBinLowEdge(iBinRec);
      double ptRec_up = H2D_response.GetYaxis()->GetBinLowEdge(iBinRec+1);
      // double xPtRecWidth = H2D_response->GetXaxis()->GetBinWidth(iBinRec);
      // if (ibinZero + (iBinRec - iBinGen) <= H1D_fluctuations_highRes->GetNbinsX()) { // make sure it's within bin range
      H2D_response.SetBinContent(iBinGen, iBinRec, H1D_fluctuations->Integral(H1D_fluctuations->GetXaxis()->FindBin(ptRec_low - ptGen + GLOBAL_epsilon), H1D_fluctuations->GetXaxis()->FindBin(ptRec_up - ptGen + GLOBAL_epsilon))); 
      // cout << "FluctResp(" << iBinRec << ", " << iBinGen << ") = " << H2D_response.GetBinContent(iBinRec, iBinGen) << endl;
        // H2D_response_highRes->SetBinError(iBinRec, iBinGen, something); 
    }
  }

  // }
  //========================================= Build response matrix end =========================================//

  H2D_jetPtResponseMatrix_fluctuations = (TH2D*)H2D_response.Clone("H2D_jetPtResponseMatrix_fluctuations"+partialUniqueSpecifier); // is replaced below with pT binning from analysis; rebin 2D with variable bins doesn't work, need to use custom function RebinVariableBins2D
  // TH2D H2D_jetPtResponseMatrix_fluctuations_rebinned = RebinVariableBins2D(H2D_jetPtResponseMatrix_fluctuations, nBinPtJetsRec[iRadius], nBinPtJetsGen[iRadius], ptBinsJetsRec[iRadius], ptBinsJetsGen[iRadius]);

  // cout << "H2D_jetPtResponseMatrix_fluctuations nbins X =" << H2D_jetPtResponseMatrix_fluctuations->GetNbinsX() << endl;
  // cout << "H2D_jetPtResponseMatrix_fluctuations_rebinned nbins X =" << H2D_jetPtResponseMatrix_fluctuations_rebinned.GetNbinsX() << endl;
  // H2D_jetPtResponseMatrix_fluctuations = (TH2D*)RebinVariableBins2D(H2D_jetPtResponseMatrix_fluctuations, nBinPtJetsRec[iRadius], nBinPtJetsGen[iRadius], ptBinsJetsRec[iRadius], ptBinsJetsGen[iRadius]).Clone("H2D_jetPtResponseMatrix_fluctuations_rebinned_"+partialUniqueSpecifier);

}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////// Spectrum Unfolding functions ///////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


void Get_Pt_spectrum_unfolded_preWidthScaling(TH1D* &H1D_jetPt_unfolded, int iDataset, int iRadius, float* centRange, const char options[]) {
  //for now makes the assumption gen and rec have the same pT binning

  // one sees that if using width scaling before putting it in the unfolding, width changes in the bin array will lead to bad effects (if width is twice larger, the bins will have a twice higher multiplication factor)
  // actually only seems to be an issue if I then scale the unfolded hist by binwidth
  // but I'm still worried about bin migration accounting if width have different width, so might as well do it at the end to be sure
  // TH1D* measured;
  // Get_Pt_spectrum_bkgCorrected_preWidthScaling(measured, iDataset, iRadius, centRange, options);
  // TH1D* mcp;
  // Get_Pt_spectrum_mcp_preWidthScaling(mcp, iDataset, iRadius, options);
  // TH1D* mcdMatched;
  // Get_Pt_spectrum_mcdMatched_preWidthScaling(mcdMatched, iDataset, iRadius, options);

    TH1D* measured;
  Get_Pt_spectrum_bkgCorrected_preWidthScaling(measured, iDataset, iRadius, centRange, "");
  TH1D* mcp;
  Get_Pt_spectrum_mcp_preWidthScaling(mcp, iDataset, iRadius, "");
  TH1D* mcdMatched;
  Get_Pt_spectrum_mcdMatched_preWidthScaling(mcdMatched, iDataset, iRadius, "");
  // Get_Pt_spectrum_mcpMatched_preWidthScaling(mcdMatched, iDataset, iRadius, "");


  // TH1D* measured;
  // Get_Pt_spectrum_bkgCorrected(measured, iDataset, iRadius, centRange, options);
  // TH1D* mcp;
  // Get_Pt_spectrum_mcp(mcp, iDataset, iRadius, options);
  // TH1D* mcdMatched;
  // Get_Pt_spectrum_mcdMatched(mcdMatched, iDataset, iRadius, options);

  TH2D* H2D_jetPtResponseMatrix_fluctuations;
  TH2D* H2D_jetPtResponseMatrix_detectorResponse;
  TH2D* H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined;
  
  Get_PtResponseMatrix_Fluctuations(H2D_jetPtResponseMatrix_fluctuations, iDataset, iRadius, centRange);
  Get_PtResponseMatrix_detectorResponse(H2D_jetPtResponseMatrix_detectorResponse, iDataset, iRadius);
  // compute matrixFluctuations times matrixDetector
  Get_PtResponseMatrix_DetectorAndFluctuationsCombined(H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined, H2D_jetPtResponseMatrix_detectorResponse, H2D_jetPtResponseMatrix_fluctuations, iDataset, iRadius, centRange);

  // // test
  // Get_PtResponseMatrix_detectorResponse(H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined, iDataset, iRadius);


  //  = (TH1D*)H1D_jetPt->Clone("mcp_jetPt_"+RadiusLegend[iRadius]+Form("%.1f", centRange[0])+"<cent<"+Form("%.1f", centRange[1]));
  // TH1D* mcp_rebinned = (TH1D*)H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined->Rebin2D(1., "mcp_rebinned_jetPt_"+RadiusLegend[iRadius]+Form("%.1f", centRange[0])+"<cent<"+Form("%.1f", centRange[1])); // rebin 2D here hasn't been parmetered well; this is the template for rebin() for 1D only
  
  // // TH1D* projX = (TH1D*)H2D_jetPtResponseMatrix_fluctuations[iRadius]->ProjectionX("projX",0,-1);
  // // TH1D* projY = (TH1D*)H2D_jetPtResponseMatrix_fluctuations[iRadius]->ProjectionY("projY",0,-1);
  // cout << "Debug log for unfolding - start" << endl;
  // for(int ibinPt = 1; ibinPt <= nBinPtJetsRec[iRadius]; ibinPt++){
  //   cout << "measured(" << ibinPt << ") = " << measured->GetBinContent(ibinPt) << endl;
  //   cout << "mcp_rebinned(" << ibinPt << ") = " << mcp_rebinned->GetBinContent(ibinPt) << endl;
  //   for(int jbinPt = 1; jbinPt <= nBinPtJetsRec[iRadius]; jbinPt++){
  //     cout << "H2D_jetPtResponseMatrix_fluctuations(" << ibinPt << "," << jbinPt << ") = " << H2D_jetPtResponseMatrix_fluctuations->GetBinContent(ibinPt,jbinPt) << endl;
  //   }
  // }
  // // cout << "Debug log for unfolding - end" << endl;
  // TH1D* responseProjRec = (TH1D*) H2D_jetPtResponseMatrix_detectorResponse->ProjectionX("responseProjRec"+RadiusLegend[iRadius]+Form("%.1f", centRange[0])+"<cent<"+Form("%.1f", centRange[1]), 0, -1);
  // TH1D* responseProjGen = (TH1D*) H2D_jetPtResponseMatrix_detectorResponse->ProjectionY("responseProjGen"+RadiusLegend[iRadius]+Form("%.1f", centRange[0])+"<cent<"+Form("%.1f", centRange[1]), 0, -1);
  // RooUnfoldResponse* response = new RooUnfoldResponse(responseProjRec, responseProjGen, H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined); // measured and mcp_rebinned are here to take inneficiencies and fakes into account
  // actually re-reading the bit of code in unfoldSpec_statVariations, the projections are just used to go from thnsparse (with only 1 dim anyway) to th1D format. No actual projection is done
  // and in the code the gen "projection" is actually the generated jets (mcp jets) before matching (so if it's not matched to a mcd jet the "projection" histogram still get an entry) . I.E. I should just use measured and mcp like below (https://github.com/alisw/AliPhysics/blob/686b64904ff2ab14cf84c21a7f6a9334cd01b76a/PWGJE/UserTasks/AliAnalysisTaskIDFFTCF.cxx#L1963)

  ////////////////////////////////////////
  //// setting up the response matrix ////
  ////////////////////////////////////////
  //    RooUnfoldResponse constructor - create from already-filled histograms
  //  "response" gives the response matrix, measured X truth.
  //  "measured" and "truth" give the projections of "response" onto the X-axis and Y-axis respectively,
  //  but with additional entries in "measured" for measurements with no corresponding truth (fakes/background) and
  //  in "truth" for unmeasured events (inefficiency).
  //  "measured" and/or "truth" can be specified as 0 (1D case only) or an empty histograms (no entries) as a shortcut
  //  to indicate, respectively, no fakes and/or no inefficiency.


  // based on https://github.com/alisw/AliPhysics/blob/master/PWGLF/FORWARD/analysis2/scripts/UnfoldMult.C
  TH1D* responseProjY = H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined->ProjectionY(H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined->GetName()+(TString)"_projY", 1, H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined->GetNbinsY());
  TH1D* responseProjX = H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined->ProjectionX(H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined->GetName()+(TString)"_projX", 1, H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined->GetNbinsX());
  TH2D* responseTranspose = (TH2D*)GetTransposeHistogram(H2D_jetPtResponseMatrix_fluctuations).Clone(H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined->GetName()+(TString)"_transposed");

  // RooUnfoldResponse* response = new RooUnfoldResponse(responseProjY, responseProjX, responseTranspose); // measured and mcp_rebinned are here to take inneficiencies and fakes into account; or is it really what's happening? 'Alternatively, the response matrix can be constructed from a pre-existing TH2D 2-dimensional histogram (with truth and measured distribution TH1D histograms for normalisation).' from https://hepunx.rl.ac.uk/~adye/software/unfold/RooUnfold.html, and I'm already normalising so maybe tehre's no need for more normalisation
  RooUnfoldResponse* response = new RooUnfoldResponse(mcdMatched, mcp, responseTranspose); // measured and mcp_rebinned are here to take inneficiencies and fakes into account; or is it really what's happening? 'Alternatively, the response matrix can be constructed from a pre-existing TH2D 2-dimensional histogram (with truth and measured distribution TH1D histograms for normalisation).' from https://hepunx.rl.ac.uk/~adye/software/unfold/RooUnfold.html, and I'm already normalising so maybe tehre's no need for more normalisation

  // // meh not sure it's true ; first two arguments of RooUnfoldResponse should be the reco and truth events (here jets) used to when filling the response matrix
  // RooUnfoldResponse* response = new RooUnfoldResponse(mcdMatched, mcp, H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined); // measured and mcp_rebinned are here to take inneficiencies and fakes into account; or is it really what's happening? 'Alternatively, the response matrix can be constructed from a pre-existing TH2D 2-dimensional histogram (with truth and measured distribution TH1D histograms for normalisation).' from https://hepunx.rl.ac.uk/~adye/software/unfold/RooUnfold.html, and I'm already normalising so maybe tehre's no need for more normalisation
  // RooUnfoldResponse* response = new RooUnfoldResponse(mcdMatched, 0, H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined); // measured and mcp_rebinned are here to take inneficiencies and fakes into account
  // RooUnfoldResponse* response = new RooUnfoldResponse(projX, projY, H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined, "", ""); // measured and mcp_rebinned are here to take inneficiencies and fakes into account
  cout << "Get_Pt_spectrum_unfolded(): should I use response->UseOverflow() ? using it gives a ratio unfolded/mcp much higher than without using it" << endl;
  // response->UseOverflow();

  RooUnfold* unfold = new RooUnfoldBayes(response, measured, 4);
  // RooUnfold* unfold = new RooUnfoldBinByBin(response, measured);
  // RooUnfold* unfold = new RooUnfoldSvd(response, measured, 8); // that's how it is done in hiroki yokoyama's thesis
  // unfold->IncludeSystematics(RooUnfolding::kAll);
  // TH1D* hist_unfold = (TH1D*)unfold.Hunfold();
  // TH1D* hist_unfold = (TH1D*)unfold->Hreco()->Clone("hist_unfold");

  // works but too simplistic; doesn't account for bin migration
  // RooUnfold* unfold = new RooUnfoldBinByBin(response, measured);
  TH1D* hist_unfold = static_cast<TH1D*>(unfold->Hreco(RooUnfold::kCovariance));

  // RooUnfoldBinByBin unfold(*H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined, measured);
  // TH1D* hist_unfold = static_cast<TH1D*>(unfold.Hreco(RooUnfold::kCovariance));


  TString partialUniqueSpecifier = Datasets[iDataset]+"_R="+Form("%.1f",arrayRadius[iRadius])+"_@cent["+Form("%.1f", centRange[0])+","+Form("%.1f", centRange[1])+"]";

  H1D_jetPt_unfolded = (TH1D*)hist_unfold->Clone("H1D_jetPt_unfolded"+partialUniqueSpecifier);

  cout << "mcp NbinsX = " << mcp->GetNbinsX() << endl;
  cout << "unfolded NbinsX = " << H1D_jetPt_unfolded->GetNbinsX() << endl;

  // only do this norm if Get_Pt_spectrum_bkgCorrected_preWidthScaling is given "" as options instead of "evtNorm" 
  NormaliseRawHistToNEvents(H1D_jetPt_unfolded, GetNEventsSelectedCentrality_JetFramework(file_O2Analysis_list[iDataset], centRange[0], centRange[1], trainId));

}

void Get_Pt_spectrum_unfolded(TH1D* &H1D_jetPt_unfolded, int iDataset, int iRadius, float* centRange, const char options[]) {
  Get_Pt_spectrum_unfolded_preWidthScaling(H1D_jetPt_unfolded, iDataset, iRadius, centRange, options);
  TransformRawHistToYield(H1D_jetPt_unfolded); // see width comments at beginning of Get_Pt_spectrum_unfolded
}

void Get_Pt_TestSpectrum_unfolded_preWidthScaling(TH1D* &H1D_jetPt_unfolded, int iDataset, int iRadius, float* centRange, const char options[]) {
  //for now makes the assumption gen and rec have the same pT binning

  // one sees that if using width scaling before putting it in the unfolding, width changes in the bin array will lead to bad effects (if width is twice larger, the bins will have a twice higher multiplication factor)
  // actually only seems to be an issue if I then scale the unfolded hist by binwidth
  // but I'm still worried about bin migration accounting if width have different width, so might as well do it at the end to be sure
  // TH1D* measured;
  // Get_Pt_spectrum_bkgCorrected_preWidthScaling(measured, iDataset, iRadius, centRange, options);
  // TH1D* mcp;
  // Get_Pt_spectrum_mcp_preWidthScaling(mcp, iDataset, iRadius, options);
  // TH1D* mcdMatched;
  // Get_Pt_spectrum_mcdMatched_preWidthScaling(mcdMatched, iDataset, iRadius, options);

  TH1D* measured;
  Get_Pt_spectrum_mcp_folded_preWidthScaling(measured, iDataset, iRadius, centRange, "");
  TH1D* mcp;
  Get_Pt_spectrum_mcp_preWidthScaling(mcp, iDataset, iRadius, "");
  TH1D* mcdMatched;
  Get_Pt_spectrum_mcdMatched_preWidthScaling(mcdMatched, iDataset, iRadius, "");
  // Get_Pt_spectrum_mcpMatched_preWidthScaling(mcdMatched, iDataset, iRadius, "");


  // TH1D* measured;
  // Get_Pt_spectrum_bkgCorrected(measured, iDataset, iRadius, centRange, options);
  // TH1D* mcp;
  // Get_Pt_spectrum_mcp(mcp, iDataset, iRadius, options);
  // TH1D* mcdMatched;
  // Get_Pt_spectrum_mcdMatched(mcdMatched, iDataset, iRadius, options);

  TH2D* H2D_jetPtResponseMatrix_fluctuations;
  TH2D* H2D_jetPtResponseMatrix_detectorResponse;
  TH2D* H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined;
  Get_PtResponseMatrix_Fluctuations(H2D_jetPtResponseMatrix_fluctuations, iDataset, iRadius, centRange);
  Get_PtResponseMatrix_detectorResponse(H2D_jetPtResponseMatrix_detectorResponse, iDataset, iRadius);
  // compute matrixFluctuations times matrixDetector
  Get_PtResponseMatrix_DetectorAndFluctuationsCombined(H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined, H2D_jetPtResponseMatrix_detectorResponse, H2D_jetPtResponseMatrix_fluctuations, iDataset, iRadius, centRange);

  // // test
  // Get_PtResponseMatrix_detectorResponse(H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined, iDataset, iRadius);


  //  = (TH1D*)H1D_jetPt->Clone("mcp_jetPt_"+RadiusLegend[iRadius]+Form("%.1f", centRange[0])+"<cent<"+Form("%.1f", centRange[1]));
  // TH1D* mcp_rebinned = (TH1D*)H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined->Rebin2D(1., "mcp_rebinned_jetPt_"+RadiusLegend[iRadius]+Form("%.1f", centRange[0])+"<cent<"+Form("%.1f", centRange[1])); // rebin 2D here hasn't been parmetered well; this is the template for rebin() for 1D only
  
  // // TH1D* projX = (TH1D*)H2D_jetPtResponseMatrix_fluctuations[iRadius]->ProjectionX("projX",0,-1);
  // // TH1D* projY = (TH1D*)H2D_jetPtResponseMatrix_fluctuations[iRadius]->ProjectionY("projY",0,-1);
  // cout << "Debug log for unfolding - start" << endl;
  // for(int ibinPt = 1; ibinPt <= nBinPtJetsRec[iRadius]; ibinPt++){
  //   cout << "measured(" << ibinPt << ") = " << measured->GetBinContent(ibinPt) << endl;
  //   cout << "mcp_rebinned(" << ibinPt << ") = " << mcp_rebinned->GetBinContent(ibinPt) << endl;
  //   for(int jbinPt = 1; jbinPt <= nBinPtJetsRec[iRadius]; jbinPt++){
  //     cout << "H2D_jetPtResponseMatrix_fluctuations(" << ibinPt << "," << jbinPt << ") = " << H2D_jetPtResponseMatrix_fluctuations->GetBinContent(ibinPt,jbinPt) << endl;
  //   }
  // }
  // // cout << "Debug log for unfolding - end" << endl;
  // TH1D* responseProjRec = (TH1D*) H2D_jetPtResponseMatrix_detectorResponse->ProjectionX("responseProjRec"+RadiusLegend[iRadius]+Form("%.1f", centRange[0])+"<cent<"+Form("%.1f", centRange[1]), 0, -1);
  // TH1D* responseProjGen = (TH1D*) H2D_jetPtResponseMatrix_detectorResponse->ProjectionY("responseProjGen"+RadiusLegend[iRadius]+Form("%.1f", centRange[0])+"<cent<"+Form("%.1f", centRange[1]), 0, -1);
  // RooUnfoldResponse* response = new RooUnfoldResponse(responseProjRec, responseProjGen, H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined); // measured and mcp_rebinned are here to take inneficiencies and fakes into account
  // actually re-reading the bit of code in unfoldSpec_statVariations, the projections are just used to go from thnsparse (with only 1 dim anyway) to th1D format. No actual projection is done
  // and in the code the gen "projection" is actually the generated jets (mcp jets) before matching (so if it's not matched to a mcd jet the "projection" histogram still get an entry) . I.E. I should just use measured and mcp like below (https://github.com/alisw/AliPhysics/blob/686b64904ff2ab14cf84c21a7f6a9334cd01b76a/PWGJE/UserTasks/AliAnalysisTaskIDFFTCF.cxx#L1963)

  ////////////////////////////////////////
  //// setting up the response matrix ////
  ////////////////////////////////////////
  //    RooUnfoldResponse constructor - create from already-filled histograms
  //  "response" gives the response matrix, measured X truth.
  //  "measured" and "truth" give the projections of "response" onto the X-axis and Y-axis respectively,
  //  but with additional entries in "measured" for measurements with no corresponding truth (fakes/background) and
  //  in "truth" for unmeasured events (inefficiency).
  //  "measured" and/or "truth" can be specified as 0 (1D case only) or an empty histograms (no entries) as a shortcut
  //  to indicate, respectively, no fakes and/or no inefficiency.

  // based on https://github.com/alisw/AliPhysics/blob/master/PWGLF/FORWARD/analysis2/scripts/UnfoldMult.C
  TH1D* responseProjY = H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined->ProjectionY(H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined->GetName()+(TString)"_projY", 1, H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined->GetNbinsY());
  TH1D* responseProjX = H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined->ProjectionX(H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined->GetName()+(TString)"_projX", 1, H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined->GetNbinsX());
  TH2D* responseTranspose = (TH2D*)GetTransposeHistogram(H2D_jetPtResponseMatrix_fluctuations).Clone(H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined->GetName()+(TString)"_transposed");

  RooUnfoldResponse* response = new RooUnfoldResponse(responseProjY, responseProjX, responseTranspose); // measured and mcp_rebinned are here to take inneficiencies and fakes into account; or is it really what's happening? 'Alternatively, the response matrix can be constructed from a pre-existing TH2D 2-dimensional histogram (with truth and measured distribution TH1D histograms for normalisation).' from https://hepunx.rl.ac.uk/~adye/software/unfold/RooUnfold.html, and I'm already normalising so maybe tehre's no need for more normalisation
  // RooUnfoldResponse* response = new RooUnfoldResponse(mcdMatched, 0, H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined); // measured and mcp_rebinned are here to take inneficiencies and fakes into account
  // RooUnfoldResponse* response = new RooUnfoldResponse(projX, projY, H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined, "", ""); // measured and mcp_rebinned are here to take inneficiencies and fakes into account
  cout << "Get_Pt_spectrum_unfolded(): should I use response->UseOverflow() ? using it gives a ratio unfolded/mcp much higher than without using it" << endl;
  // response->UseOverflow();

  RooUnfold* unfold = new RooUnfoldBayes(response, measured, 4);
  // RooUnfold* unfold = new RooUnfoldBinByBin(response, measured);
  // RooUnfold* unfold = new RooUnfoldSvd(response, measured, 8); // that's how it is done in hiroki yokoyama's thesis
  // unfold->IncludeSystematics(RooUnfolding::kAll);
  // TH1D* hist_unfold = (TH1D*)unfold.Hunfold();
  // TH1D* hist_unfold = (TH1D*)unfold->Hreco()->Clone("hist_unfold");

  // works but too simplistic; doesn't account for bin migration
  // RooUnfold* unfold = new RooUnfoldBinByBin(response, measured);
  TH1D* hist_unfold = static_cast<TH1D*>(unfold->Hreco(RooUnfold::kCovariance));

  // RooUnfoldBinByBin unfold(*H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined, measured);
  // TH1D* hist_unfold = static_cast<TH1D*>(unfold.Hreco(RooUnfold::kCovariance));


  TString partialUniqueSpecifier = Datasets[iDataset]+"_R="+Form("%.1f",arrayRadius[iRadius])+"_@cent["+Form("%.1f", centRange[0])+","+Form("%.1f", centRange[1])+"]";

  H1D_jetPt_unfolded = (TH1D*)hist_unfold->Clone("H1D_jetPt_unfolded"+partialUniqueSpecifier);

  cout << "mcp NbinsX = " << mcp->GetNbinsX() << endl;
  cout << "unfolded NbinsX = " << H1D_jetPt_unfolded->GetNbinsX() << endl;

  // only do this norm if Get_Pt_spectrum_bkgCorrected_preWidthScaling is given "" as options instead of "evtNorm" 
  // NormaliseRawHistToNEvents(H1D_jetPt_unfolded, GetNEventsSelectedCentrality_JetFramework(file_O2Analysis_list[iDataset], centRange[0], centRange[1], trainId));
  NormaliseRawHistToNEvents(H1D_jetPt_unfolded, GetNEventsSelected_JetFramework_weighted(file_O2Analysis_ppSimDetectorEffect));
}

void Get_Pt_TestSpectrum_unfolded(TH1D* &H1D_jetPt_unfolded, int iDataset, int iRadius, float* centRange, const char options[]) {
  Get_Pt_TestSpectrum_unfolded_preWidthScaling(H1D_jetPt_unfolded, iDataset, iRadius, centRange, options);
  TransformRawHistToYield(H1D_jetPt_unfolded); // see width comments at beginning of Get_Pt_spectrum_unfolded
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////// Spectrum plotting functions ///////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


void Draw_Pt_spectrum_raw(int iDataset, int iRadius, const char options[]) {

  TH1D* H1D_jetPt_raw[nCentralityBins];
  float centRange[2];
  for(int iCent = 0; iCent < nCentralityBins; iCent++){
    // centRange[0] = arrayCentralityBinning[iCent];
    // centRange[1] = arrayCentralityBinning[iCent+1];
    centRange[0] = arrayCentralityIntervals[iCent][0];
    centRange[1] = arrayCentralityIntervals[iCent][1];

    Get_Pt_spectrum_bkgCorrected(H1D_jetPt_raw[iCent], iDataset, iRadius, centRange, options);
  }


  TString* pdfName = new TString("jet_"+jetType[iJetType]+"_"+jetLevel[iJetLevel]+"_"+Datasets[iDataset]+"_R="+Form("%.1f", arrayRadius[iRadius])+"_Pt_raw");

  TString textContext(contextCustomOneField(*texDatasetsComparisonCommonDenominator, ""));

  TString centralityLegend[nCentralityBins]; CentralityLegend(centralityLegend, arrayCentralityIntervals, nCentralityBins);

  TString* yAxisLabel;
  if (strstr(options, "evtNorm") != NULL) {
    yAxisLabel = texJetPtYield_EventNorm;
  }
  if (strstr(options, "entriesNorm") != NULL) {
    yAxisLabel = texJetPtYield_EntriesNorm;
  }

  // Draw_TH1_Histograms_in_one(H1D_jetPt_raw, centralityLegend, nCentralityBins, textContext, pdfName, texPtX, yAxisLabel, texCollisionDataInfo, "logy");
  Draw_TH1_Histograms_in_one(H1D_jetPt_raw, centralityLegend, nCentralityBins, textContext, pdfName, texPtX, yAxisLabel, texCollisionDataInfo, drawnWindowAuto, "logy");
}

void Draw_Pt_spectrum_mcp(int iDataset, int iRadius, const char options[]) {

  TH1D* H1D_jetPt_mcp_genBinning;
  TH1D* H1D_jetPt_mcp_recBinning;
  TH1D* H1D_jetPt_mcp_collection[2];
  Get_Pt_spectrum_mcp(H1D_jetPt_mcp_genBinning, iDataset, iRadius, options);
  Get_Pt_spectrum_mcp_recBinning(H1D_jetPt_mcp_recBinning, iDataset, iRadius, options);
  H1D_jetPt_mcp_collection[0] = H1D_jetPt_mcp_genBinning;
  H1D_jetPt_mcp_collection[1] = H1D_jetPt_mcp_recBinning;


  TString* pdfName = new TString("jet_"+jetType[iJetType]+"_"+jetLevel[iJetLevel]+"_"+Datasets[iDataset]+"_R="+Form("%.1f", arrayRadius[iRadius])+"_Pt_mcp");

  TString textContext(contextCustomOneField(*texDatasetsComparisonCommonDenominator, ""));

  TString* yAxisLabel;
  if (strstr(options, "evtNorm") != NULL) {
    yAxisLabel = texJetPtYield_EventNorm;
  }
  if (strstr(options, "entriesNorm") != NULL) {
    yAxisLabel = texJetPtYield_EntriesNorm;
  }
  TString genVsRecBinningLegend[2] = {"gen binning", "rec binning"};

  Draw_TH1_Histograms_in_one(H1D_jetPt_mcp_collection, genVsRecBinningLegend, 2, textContext, pdfName, texPtX, yAxisLabel, texCollisionDataInfo, drawnWindowAuto, "logy"); 
}

void Draw_Pt_spectrum_mcdMatched(int iDataset, int iRadius, const char options[]) {

  TH1D* H1D_jetPt_mcp;
  Get_Pt_spectrum_mcdMatched(H1D_jetPt_mcp, iDataset, iRadius, options);


  TString* pdfName = new TString("jet_"+jetType[iJetType]+"_"+jetLevel[iJetLevel]+"_"+Datasets[iDataset]+"_R="+Form("%.1f", arrayRadius[iRadius])+"_Pt_mcdMatched");

  TString textContext(contextCustomOneField(*texDatasetsComparisonCommonDenominator, ""));

  TString* yAxisLabel;
  if (strstr(options, "evtNorm") != NULL) {
    yAxisLabel = texJetPtYield_EventNorm;
  }
  if (strstr(options, "entriesNorm") != NULL) {
    yAxisLabel = texJetPtYield_EntriesNorm;
  }

  Draw_TH1_Histogram(H1D_jetPt_mcp, textContext, pdfName, texPtX, yAxisLabel, texCollisionDataInfo, drawnWindowAuto, "logy");
}


void Draw_Pt_efficiency_jets(int iDataset, int iRadius, const char options[]) {
  TH1D* H1D_jetPt_mcp;
  TH1D* H1D_jetPt_mcpMatched;
  TH1D* H1D_jetEfficiency;
  bool divideSuccess;
  Get_Pt_spectrum_mcp(H1D_jetPt_mcp, iDataset, iRadius, options);
  Get_Pt_spectrum_mcpMatched(H1D_jetPt_mcpMatched, iDataset, iRadius, options);

  H1D_jetEfficiency = (TH1D*)H1D_jetPt_mcpMatched->Clone("H1D_jetEfficiency"+Datasets[iDataset]+"_R="+Form("%.1f", arrayRadius[iRadius]));
  divideSuccess = H1D_jetEfficiency->Divide(H1D_jetPt_mcp);

  TString* pdfName = new TString("jet_"+jetType[iJetType]+"_"+jetLevel[iJetLevel]+"_"+Datasets[iDataset]+"_R="+Form("%.1f", arrayRadius[iRadius])+"_Pt_efficiency");

  TString textContext(contextCustomOneField(*texDatasetsComparisonCommonDenominator, ""));

  if (divideSuccess){
    Draw_TH1_Histogram(H1D_jetEfficiency, textContext, pdfName, texPtX, texJetEfficiency, texCollisionDataInfo, drawnWindowAuto, "efficiency");
  }
}

void Draw_ResponseMatrices_Fluctuations(int iDataset, int iRadius) {

  TH2D* H2D_jetPtResponseMatrix_fluctuations[nCentralityBins];

  float centRange[2];
  for(int iCent = 0; iCent < nCentralityBins; iCent++){
    // centRange[0] = arrayCentralityBinning[iCent];
    // centRange[1] = arrayCentralityBinning[iCent+1];
    centRange[0] = arrayCentralityIntervals[iCent][0];
    centRange[1] = arrayCentralityIntervals[iCent][1];
    Get_PtResponseMatrix_Fluctuations(H2D_jetPtResponseMatrix_fluctuations[iCent], iDataset, iRadius, centRange);
  }
  TString* pdfName_logz = new TString("background_responseMatrices_"+jetType[iJetType]+"_"+Datasets[iDataset]+"_centralityComp_logz");
  TString* pdfNameFullRes_logz = new TString("background_responseMatrices_"+jetType[iJetType]+"_"+Datasets[iDataset]+"_centralityComp_FullRes_logz");
  TString* pdfName = new TString("background_responseMatrices_"+jetType[iJetType]+"_"+Datasets[iDataset]+"_centralityComp");
  TString* pdfNameFullRes = new TString("background_responseMatrices_"+jetType[iJetType]+"_"+Datasets[iDataset]+"_centralityComp_FullRes");

  TString textContext(contextCustomOneField(*texDatasetsComparisonCommonDenominator, ""));

  TString centralityLegend[nCentralityBins]; CentralityLegend(centralityLegend, arrayCentralityIntervals, nCentralityBins);

  // Draw_TH2_Histograms(H2D_jetPtResponseMatrix_fluctuations, centralityLegend, nCentralityBins, textContext, pdfName, texPtJetRecX, texPtJetGenX, texCollisionDataInfo, "logz");

  Draw_TH2_Histograms(H2D_jetPtResponseMatrix_fluctuations, centralityLegend, nCentralityBins, textContext, pdfName_logz, texPtJetGenX, texPtJetRecX, texCollisionDataInfo, drawnWindowAuto, "logz");
  Draw_TH2_Histograms(H2D_jetPtResponseMatrix_fluctuations, centralityLegend, nCentralityBins, textContext, pdfName, texPtJetGenX, texPtJetRecX, texCollisionDataInfo, drawnWindowAuto, "");
}

void Draw_ResponseMatrices_detectorResponse(int iDataset, int iRadius) {

  TH2D* H2D_jetPtResponseMatrix_detectorResponse;

  Get_PtResponseMatrix_detectorResponse(H2D_jetPtResponseMatrix_detectorResponse, iDataset, iRadius);

  TString* pdfName = new TString("detectorResponse_responseMatrices_"+jetType[iJetType]+"_"+Datasets[iDataset]);
  TString* pdfName_logz = new TString("detectorResponse_responseMatrices_"+jetType[iJetType]+"_"+Datasets[iDataset]+"_logz");

  TString textContext(contextCustomOneField(*texDatasetsComparisonCommonDenominator, ""));

  // Draw_TH2_Histograms(H2D_jetPtResponseMatrix_detectorResponse, centralityLegend, nCentralityBins, textContext, pdfName, texPtJetRecX, texPtJetGenX, texCollisionDataInfo, "logz");

  Draw_TH2_Histogram(H2D_jetPtResponseMatrix_detectorResponse, textContext, pdfName, texPtJetGenX, texPtJetRecX, texCollisionDataInfo, drawnWindowAuto, "");
  Draw_TH2_Histogram(H2D_jetPtResponseMatrix_detectorResponse, textContext, pdfName_logz, texPtJetGenX, texPtJetRecX, texCollisionDataInfo, drawnWindowAuto, "logz");
}

void Draw_ResponseMatrices_DetectorAndFluctuationsCombined(int iDataset, int iRadius) {

  TH2D* H2D_jetPtResponseMatrix_fluctuations[nCentralityBins];
  TH2D* H2D_jetPtResponseMatrix_detectorResponse[nCentralityBins];
  TH2D* H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined[nCentralityBins];

  float centRange[2];
  for(int iCent = 0; iCent < nCentralityBins; iCent++){
    // centRange[0] = arrayCentralityBinning[iCent];
    // centRange[1] = arrayCentralityBinning[iCent+1];
    centRange[0] = arrayCentralityIntervals[iCent][0];
    centRange[1] = arrayCentralityIntervals[iCent][1];

    Get_PtResponseMatrix_Fluctuations(H2D_jetPtResponseMatrix_fluctuations[iCent], iDataset, iRadius, centRange);
    Get_PtResponseMatrix_detectorResponse(H2D_jetPtResponseMatrix_detectorResponse[iCent], iDataset, iRadius);
    Get_PtResponseMatrix_DetectorAndFluctuationsCombined(H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined[iCent], H2D_jetPtResponseMatrix_detectorResponse[iCent], H2D_jetPtResponseMatrix_fluctuations[iCent], iDataset, iRadius, centRange);
  }
  TString* pdfName = new TString("detectorAndFluctuationsCombined_responseMatrices_"+jetType[iJetType]+"_"+Datasets[iDataset]+"_centralityComp");
  TString* pdfName_logz = new TString("detectorAndFluctuationsCombined_responseMatrices_"+jetType[iJetType]+"_"+Datasets[iDataset]+"_centralityComp_logz");

  TString textContext(contextCustomOneField(*texDatasetsComparisonCommonDenominator, ""));

  TString centralityLegend[nCentralityBins]; CentralityLegend(centralityLegend, arrayCentralityIntervals, nCentralityBins);

  // Draw_TH2_Histograms(H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined, centralityLegend, nCentralityBins, textContext, pdfName, texPtJetRecX, texPtJetGenX, texCollisionDataInfo, "logz");

  Draw_TH2_Histograms(H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined, centralityLegend, nCentralityBins, textContext, pdfName, texPtJetGenX, texPtJetRecX, texCollisionDataInfo, drawnWindowAuto, "");
  Draw_TH2_Histograms(H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined, centralityLegend, nCentralityBins, textContext, pdfName_logz, texPtJetGenX, texPtJetRecX, texCollisionDataInfo, drawnWindowAuto, "logz");
}

void Draw_Pt_spectrum_unfolded(int iDataset, int iRadius, const char options[]) {

  TH1D* H1D_jetPt_measured[nCentralityBins];
  TH1D* H1D_jetPt_unfolded[nCentralityBins];
  TH1D* H1D_jetPt_unfolded_mcpComp[nCentralityBins][2];
  TH1D* H1D_jetPt_unfolded_measuredComp[nCentralityBins][2];
  TH1D* H1D_jetPt_mcp;
  TH1D* H1D_jetPt_ratio_mcp[nCentralityBins];
  TH1D* H1D_jetPt_ratio_measured[nCentralityBins];
  bool divideSuccessMcp[nCentralityBins];
  bool divideSuccessMeasured[nCentralityBins];
  float centRange[2];
  TString partialUniqueSpecifier;

  Get_Pt_spectrum_mcp(H1D_jetPt_mcp, iDataset, iRadius, options);
  // H1D_jetPt_mcp = (TH1D*)H1D_jetPt_mcp_defaultBin->Rebin(nBinPtJetsGen[iRadius],"jetPt_mcp_rebinned_"+RadiusLegend[iRadius]+Datasets[iDataset], ptBinsJetsGen[iRadius]);

  for(int iCent = 0; iCent < nCentralityBins; iCent++){
    centRange[0] = arrayCentralityIntervals[iCent][0];
    centRange[1] = arrayCentralityIntervals[iCent][1];
    partialUniqueSpecifier = Datasets[iDataset]+"_R="+Form("%.1f",arrayRadius[iRadius])+"_@cent["+Form("%.1f", centRange[0])+","+Form("%.1f", centRange[1])+"]";
    
    Get_Pt_spectrum_unfolded(H1D_jetPt_unfolded[iCent], iDataset, iRadius, centRange, options);

    // comparison with raw measured
    Get_Pt_spectrum_bkgCorrected_genBinning(H1D_jetPt_measured[iCent], iDataset, iRadius, centRange, options);
    H1D_jetPt_unfolded_measuredComp[iCent][0] = (TH1D*)H1D_jetPt_unfolded[iCent]->Clone("H1D_jetPt_unfolded_measuredComp"+partialUniqueSpecifier);
    H1D_jetPt_unfolded_measuredComp[iCent][1] = (TH1D*)H1D_jetPt_measured[iCent]->Clone("H1D_jetPt_measured_measuredComp"+partialUniqueSpecifier);
    H1D_jetPt_ratio_measured[iCent] = (TH1D*)H1D_jetPt_measured[iCent]->Clone("H1D_jetPt_ratio_measured"+partialUniqueSpecifier);
    divideSuccessMeasured[iCent] = H1D_jetPt_ratio_measured[iCent]->Divide(H1D_jetPt_unfolded[iCent]);
    for(int iBin = 1; iBin <= H1D_jetPt_ratio_measured[iCent]->GetNbinsX(); iBin++){
      cout << "ratio measured/unfolded bin contents: " << H1D_jetPt_ratio_measured[iCent]->GetBinContent(iBin) << endl;
      cout << "nevents rec, mcd, mcp = " << GetNEventsSelectedCentrality_JetFramework(file_O2Analysis_list[iDataset], centRange[0], centRange[1], trainId) << ", " << GetNEventsSelected_JetFramework_weighted(file_O2Analysis_ppSimDetectorEffect) << ", " << GetNEventsSelected_JetFramework_weighted(file_O2Analysis_ppSimDetectorEffect)  << endl;
    }

    // comparison with mcp truth
    H1D_jetPt_unfolded_mcpComp[iCent][0] = (TH1D*)H1D_jetPt_unfolded[iCent]->Clone("H1D_jetPt_unfolded_mcpComp"+partialUniqueSpecifier);
    H1D_jetPt_unfolded_mcpComp[iCent][1] = (TH1D*)H1D_jetPt_mcp->Clone("H1D_jetPt_unfolded_mcpComp"+partialUniqueSpecifier);
    H1D_jetPt_ratio_mcp[iCent] = (TH1D*)H1D_jetPt_mcp->Clone("H1D_jetPt_ratio_mcp"+partialUniqueSpecifier);
    divideSuccessMcp[iCent] = H1D_jetPt_ratio_mcp[iCent]->Divide(H1D_jetPt_unfolded[iCent]);
    for(int iBin = 1; iBin <= H1D_jetPt_ratio_mcp[iCent]->GetNbinsX(); iBin++){
      cout << "ratio mcp/unfolded bin contents: " << H1D_jetPt_ratio_mcp[iCent]->GetBinContent(iBin) << endl;
      cout << "nevents rec, mcd, mcp = " << GetNEventsSelectedCentrality_JetFramework(file_O2Analysis_list[iDataset], centRange[0], centRange[1], trainId) << ", " << GetNEventsSelected_JetFramework_weighted(file_O2Analysis_ppSimDetectorEffect) << ", " << GetNEventsSelected_JetFramework_weighted(file_O2Analysis_ppSimDetectorEffect)  << endl;
    }
  }

  TString* pdfName = new TString("jet_"+jetType[iJetType]+"_"+jetLevel[iJetLevel]+"_"+Datasets[iDataset]+"_R="+Form("%.1f", arrayRadius[iRadius])+"_Pt_unfolded");

  TString textContext(contextCustomOneField(*texDatasetsComparisonCommonDenominator, ""));

  TString centralityLegend[nCentralityBins]; CentralityLegend(centralityLegend, arrayCentralityIntervals, nCentralityBins);

  TString* yAxisLabel;
  if (strstr(options, "evtNorm") != NULL) {
    yAxisLabel = texJetPtYield_EventNorm;
  }
  if (strstr(options, "entriesNorm") != NULL) {
    yAxisLabel = texJetPtYield_EntriesNorm;
  }
  // Draw_TH1_Histograms_in_one(H1D_jetPt_unfolded, centralityLegend, nCentralityBins, textContext, pdfName, texPtX, yAxisLabel, texCollisionDataInfo, "logy");
  Draw_TH1_Histograms_in_one(H1D_jetPt_unfolded, centralityLegend, nCentralityBins, textContext, pdfName, texPtX, yAxisLabel, texCollisionDataInfo, drawnWindowAuto, "logy");


    // comparison with raw measured
  TString unfoldedMeasuredCompLegend[2] = {"unfolded data", "measured raw"};
  for(int iCent = 0; iCent < nCentralityBins; iCent++){
    TString* pdfName_measuredComp = new TString("jet_"+jetType[iJetType]+"_"+jetLevel[iJetLevel]+"_"+Datasets[iDataset]+"_R="+Form("%.1f", arrayRadius[iRadius])+"_@cent["+Form("%.1f", centRange[0])+","+Form("%.1f", centRange[1])+"]"+"_Pt_unfolded_measuredComp");
    Draw_TH1_Histograms_in_one(H1D_jetPt_unfolded_measuredComp[iCent], unfoldedMeasuredCompLegend, 2, textContext, pdfName_measuredComp, texPtX, yAxisLabel, texCollisionDataInfo, drawnWindowAuto, "logy");
  }

  bool divideSuccessMeasured_boolsum = true;
  for(int iCent = 0; iCent < nCentralityBins; iCent++){
    if (!(divideSuccessMeasured[iCent])) {
      divideSuccessMeasured_boolsum = false;
    }
  }
  if (divideSuccessMeasured_boolsum){
    TString* pdfName_ratio_measured = new TString("jet_"+jetType[iJetType]+"_"+jetLevel[iJetLevel]+"_"+Datasets[iDataset]+"_R="+Form("%.1f", arrayRadius[iRadius])+"_Pt_unfolded_ratioMeasured");
    Draw_TH1_Histograms_in_one(H1D_jetPt_ratio_measured, centralityLegend, nCentralityBins, textContext, pdfName_ratio_measured, texPtX, texRatioMeasuredUnfolded, texCollisionDataInfo, drawnWindowAuto, "standardratio,ratioLine");
  }


    // comparison with mcp truth
  TString unfoldedTruthCompLegend[2] = {"unfolded data", "mcp truth"};
  for(int iCent = 0; iCent < nCentralityBins; iCent++){
    TString* pdfName_mcpComp = new TString("jet_"+jetType[iJetType]+"_"+jetLevel[iJetLevel]+"_"+Datasets[iDataset]+"_R="+Form("%.1f", arrayRadius[iRadius])+"_@cent["+Form("%.1f", centRange[0])+","+Form("%.1f", centRange[1])+"]"+"_Pt_unfolded_mcpComp");
    Draw_TH1_Histograms_in_one(H1D_jetPt_unfolded_mcpComp[iCent], unfoldedTruthCompLegend, 2, textContext, pdfName_mcpComp, texPtX, yAxisLabel, texCollisionDataInfo, drawnWindowAuto, "logy");
  }

  bool divideSuccessMcp_boolsum = true;
  for(int iCent = 0; iCent < nCentralityBins; iCent++){
    if (!(divideSuccessMcp[iCent])) {
      divideSuccessMcp_boolsum = false;
    }
  }
  if (divideSuccessMcp_boolsum){
    TString* pdfName_ratio_mcp = new TString("jet_"+jetType[iJetType]+"_"+jetLevel[iJetLevel]+"_"+Datasets[iDataset]+"_R="+Form("%.1f", arrayRadius[iRadius])+"_Pt_unfolded_ratioMcp");
    Draw_TH1_Histograms_in_one(H1D_jetPt_ratio_mcp, centralityLegend, nCentralityBins, textContext, pdfName_ratio_mcp, texPtX, texRatioMcpUnfolded, texCollisionDataInfo, drawnWindowAuto, "standardratio,ratioLine");
  }
}

void Draw_Pt_TestSpectrum_unfolded(int iDataset, int iRadius, const char options[]) {

  TH1D* H1D_jetPt_measured[nCentralityBins];
  TH1D* H1D_jetPt_unfolded[nCentralityBins];
  TH1D* H1D_jetPt_unfolded_mcpComp[nCentralityBins][2];
  TH1D* H1D_jetPt_unfolded_measuredComp[nCentralityBins][2];
  TH1D* H1D_jetPt_mcp;
  TH1D* H1D_jetPt_ratio_mcp[nCentralityBins];
  TH1D* H1D_jetPt_ratio_measured[nCentralityBins];
  bool divideSuccessMcp[nCentralityBins];
  bool divideSuccessMeasured[nCentralityBins];
  float centRange[2];
  TString partialUniqueSpecifier;

  Get_Pt_spectrum_mcp(H1D_jetPt_mcp, iDataset, iRadius, options);
  // H1D_jetPt_mcp = (TH1D*)H1D_jetPt_mcp_defaultBin->Rebin(nBinPtJetsGen[iRadius],"jetPt_mcp_rebinned_"+RadiusLegend[iRadius]+Datasets[iDataset], ptBinsJetsGen[iRadius]);

  for(int iCent = 0; iCent < nCentralityBins; iCent++){
    centRange[0] = arrayCentralityIntervals[iCent][0];
    centRange[1] = arrayCentralityIntervals[iCent][1];
    partialUniqueSpecifier = Datasets[iDataset]+"_R="+Form("%.1f",arrayRadius[iRadius])+"_@cent["+Form("%.1f", centRange[0])+","+Form("%.1f", centRange[1])+"]";
    
    Get_Pt_TestSpectrum_unfolded(H1D_jetPt_unfolded[iCent], iDataset, iRadius, centRange, options);
    // comparison with raw measured
    Get_Pt_spectrum_mcp_folded_genBinning(H1D_jetPt_measured[iCent], iDataset, iRadius, centRange, options);

    H1D_jetPt_unfolded_measuredComp[iCent][0] = (TH1D*)H1D_jetPt_unfolded[iCent]->Clone("H1D_jetPt_unfolded_measuredComp"+partialUniqueSpecifier);
    H1D_jetPt_unfolded_measuredComp[iCent][1] = (TH1D*)H1D_jetPt_measured[iCent]->Clone("H1D_jetPt_measured_measuredComp"+partialUniqueSpecifier);
    H1D_jetPt_ratio_measured[iCent] = (TH1D*)H1D_jetPt_measured[iCent]->Clone("H1D_jetPt_ratio_measured"+partialUniqueSpecifier);
    divideSuccessMeasured[iCent] = H1D_jetPt_ratio_measured[iCent]->Divide(H1D_jetPt_unfolded[iCent]);
    for(int iBin = 1; iBin <= H1D_jetPt_ratio_measured[iCent]->GetNbinsX(); iBin++){
      cout << "ratio folded/unfolded bin contents: " << H1D_jetPt_ratio_measured[iCent]->GetBinContent(iBin) << endl;
      cout << "nevents rec, mcd, mcp = " << GetNEventsSelectedCentrality_JetFramework(file_O2Analysis_list[iDataset], centRange[0], centRange[1], trainId) << ", " << GetNEventsSelected_JetFramework_weighted(file_O2Analysis_ppSimDetectorEffect) << ", " << GetNEventsSelected_JetFramework_weighted(file_O2Analysis_ppSimDetectorEffect)  << endl;
    }

    // comparison with mcp truth
    H1D_jetPt_unfolded_mcpComp[iCent][0] = (TH1D*)H1D_jetPt_unfolded[iCent]->Clone("H1D_jetPt_unfolded_mcpComp"+partialUniqueSpecifier);
    H1D_jetPt_unfolded_mcpComp[iCent][1] = (TH1D*)H1D_jetPt_mcp->Clone("H1D_jetPt_unfolded_mcpComp"+partialUniqueSpecifier);
    H1D_jetPt_ratio_mcp[iCent] = (TH1D*)H1D_jetPt_mcp->Clone("H1D_jetPt_ratio_mcp"+partialUniqueSpecifier);
    divideSuccessMcp[iCent] = H1D_jetPt_ratio_mcp[iCent]->Divide(H1D_jetPt_unfolded[iCent]);
    for(int iBin = 1; iBin <= H1D_jetPt_ratio_mcp[iCent]->GetNbinsX(); iBin++){
      cout << "ratio mcp/unfolded bin contents: " << H1D_jetPt_ratio_mcp[iCent]->GetBinContent(iBin) << endl;
      cout << "nevents rec, mcd, mcp = " << GetNEventsSelectedCentrality_JetFramework(file_O2Analysis_list[iDataset], centRange[0], centRange[1], trainId) << ", " << GetNEventsSelected_JetFramework_weighted(file_O2Analysis_ppSimDetectorEffect) << ", " << GetNEventsSelected_JetFramework_weighted(file_O2Analysis_ppSimDetectorEffect)  << endl;
    }
  }

  TString* pdfName = new TString("jet_"+jetType[iJetType]+"_"+jetLevel[iJetLevel]+"_"+Datasets[iDataset]+"_R="+Form("%.1f", arrayRadius[iRadius])+"_Pt_TestRefold_unfolded");

  TString textContext(contextCustomOneField(*texDatasetsComparisonCommonDenominator, ""));

  TString centralityLegend[nCentralityBins]; CentralityLegend(centralityLegend, arrayCentralityIntervals, nCentralityBins);

  TString* yAxisLabel;
  if (strstr(options, "evtNorm") != NULL) {
    yAxisLabel = texJetPtYield_EventNorm;
  }
  if (strstr(options, "entriesNorm") != NULL) {
    yAxisLabel = texJetPtYield_EntriesNorm;
  }
  // Draw_TH1_Histograms_in_one(H1D_jetPt_unfolded, centralityLegend, nCentralityBins, textContext, pdfName, texPtX, yAxisLabel, texCollisionDataInfo, "logy");
  Draw_TH1_Histograms_in_one(H1D_jetPt_unfolded, centralityLegend, nCentralityBins, textContext, pdfName, texPtX, yAxisLabel, texCollisionDataInfo, drawnWindowAuto, "logy");


    // comparison with raw measured
  TString unfoldedMeasuredCompLegend[2] = {"unfolded data", "folded truth"};
  for(int iCent = 0; iCent < nCentralityBins; iCent++){
    TString* pdfName_measuredComp = new TString("jet_"+jetType[iJetType]+"_"+jetLevel[iJetLevel]+"_"+Datasets[iDataset]+"_R="+Form("%.1f", arrayRadius[iRadius])+"_@cent["+Form("%.1f", centRange[0])+","+Form("%.1f", centRange[1])+"]"+"_Pt_TestRefold_measuredComp");
    Draw_TH1_Histograms_in_one(H1D_jetPt_unfolded_measuredComp[iCent], unfoldedMeasuredCompLegend, 2, textContext, pdfName_measuredComp, texPtX, yAxisLabel, texCollisionDataInfo, drawnWindowAuto, "logy");
  }

  bool divideSuccessMeasured_boolsum = true;
  for(int iCent = 0; iCent < nCentralityBins; iCent++){
    if (!(divideSuccessMeasured[iCent])) {
      divideSuccessMeasured_boolsum = false;
    }
  }
  if (divideSuccessMeasured_boolsum){
    TString* pdfName_ratio_measured = new TString("jet_"+jetType[iJetType]+"_"+jetLevel[iJetLevel]+"_"+Datasets[iDataset]+"_R="+Form("%.1f", arrayRadius[iRadius])+"_Pt_TestRefold_ratioMeasured");
    Draw_TH1_Histograms_in_one(H1D_jetPt_ratio_measured, centralityLegend, nCentralityBins, textContext, pdfName_ratio_measured, texPtX, texRatioMeasuredUnfolded, texCollisionDataInfo, drawnWindowAuto, "standardratio,ratioLine");
  }


    // comparison with mcp truth
  TString unfoldedTruthCompLegend[2] = {"unfolded data", "mcp truth"};
  for(int iCent = 0; iCent < nCentralityBins; iCent++){
    TString* pdfName_mcpComp = new TString("jet_"+jetType[iJetType]+"_"+jetLevel[iJetLevel]+"_"+Datasets[iDataset]+"_R="+Form("%.1f", arrayRadius[iRadius])+"_@cent["+Form("%.1f", centRange[0])+","+Form("%.1f", centRange[1])+"]"+"_Pt_TestRefold_mcpComp");
    Draw_TH1_Histograms_in_one(H1D_jetPt_unfolded_mcpComp[iCent], unfoldedTruthCompLegend, 2, textContext, pdfName_mcpComp, texPtX, yAxisLabel, texCollisionDataInfo, drawnWindowAuto, "logy");
  }

  bool divideSuccessMcp_boolsum = true;
  for(int iCent = 0; iCent < nCentralityBins; iCent++){
    if (!(divideSuccessMcp[iCent])) {
      divideSuccessMcp_boolsum = false;
    }
  }
  if (divideSuccessMcp_boolsum){
    TString* pdfName_ratio_mcp = new TString("jet_"+jetType[iJetType]+"_"+jetLevel[iJetLevel]+"_"+Datasets[iDataset]+"_R="+Form("%.1f", arrayRadius[iRadius])+"_Pt_TestRefold_ratioMcp");
    Draw_TH1_Histograms_in_one(H1D_jetPt_ratio_mcp, centralityLegend, nCentralityBins, textContext, pdfName_ratio_mcp, texPtX, texRatioMcpUnfolded, texCollisionDataInfo, drawnWindowAuto, "standardratio,ratioLine");
  }
}

// // // PhD fuction for V0s
// void Get_V0_RawYield_vsPt_FeeddownCorrected_withUnfolding(TH1D* &hRawYield_vsPt, TH3D* H3D_detectedV0s_DataLike, TH1I* H1I_SelectedEventCount, TFile *file_O2Analysis, double* SignalMean, double* SignalStandardDev, int ipart, double* pTbins, int nbinpT, int ibinXaxisCut_low,  int ibinXaxisCut_high, int warning_cutArry_ID, double SideBandSizeMultiplierModifier, int SignalExtractionType, int InvMassRebinFactor) {
//   float SelectedEventCount = H1I_SelectedEventCount->GetEntries();

//   TH1D* mcp = (TH1D*)file_O2Analysis->Get("lambdakzero-particle-count-mc/h"+NamePart[ipart]+"Count_PtDiff");
//   TH1D* mcp_rebinned = (TH1D*)mcp->Rebin(nbinpT,"mcp_rebinned",pTbins);
  
//   Get_RawYield_vsPt_FeeddownCorrected_RawCount(hRawYield_vsPt, H3D_detectedV0s_DataLike, H1I_SelectedEventCount, file_O2Analysis, SignalMean, SignalStandardDev, ipart, pTbins, nbinpT, ibinXaxisCut_low, ibinXaxisCut_high, warning_cutArry_ID, SideBandSizeMultiplierModifier, SignalExtractionType, InvMassRebinFactor);
//   TH2D *H2D_PtResponseMatrix = new TH2D("H2D_PtResponseMatrix", "H2D_PtResponseMatrix", nbinpT, pTbins, nbinpT, pTbins);
//   Get_PtResponseMatrix(H2D_PtResponseMatrix, file_O2Analysis, pTbins, nbinpT);

//   TH1D* measured = (TH1D*)hRawYield_vsPt->Clone("measured");

//   // TH1D* projX = (TH1D*)H2D_PtResponseMatrix->ProjectionX("projX",0,-1);
//   // TH1D* projY = (TH1D*)H2D_PtResponseMatrix->ProjectionY("projY",0,-1);
//   cout << "--------------------------------------------------------------------- PRE UNFOLDING TEST" << endl;
//   for(int ibinPt = 1; ibinPt <= nbinpT; ibinPt++){
//     cout << "measured(" << ibinPt << ") = " << measured->GetBinContent(ibinPt) << endl;
//     cout << "mcp_rebinned(" << ibinPt << ") = " << mcp_rebinned->GetBinContent(ibinPt) << endl;
//     for(int jbinPt = 1; jbinPt <= nbinpT; jbinPt++){
//     cout << "H2D_PtResponseMatrix(" << ibinPt << "," << jbinPt << ") = " << H2D_PtResponseMatrix->GetBinContent(ibinPt,jbinPt) << endl;
//     }
//   }

//   RooUnfoldResponse* response = new RooUnfoldResponse(measured, mcp_rebinned, H2D_PtResponseMatrix);
//   // RooUnfoldResponse* response = new RooUnfoldResponse(projX, projY, H2D_PtResponseMatrix, "", "");
//   // RooUnfold* unfold = new RooUnfoldBayes(response, measured, 4);
//   RooUnfold* unfold = new RooUnfoldBinByBin(response, measured);

//   // TH1D* hist_unfold = (TH1D*)unfold.Hunfold();
//   // TH1D* hist_unfold = (TH1D*)unfold->Hreco()->Clone("hist_unfold");
//   TH1D* hist_unfold = static_cast<TH1D*>(unfold->Hreco());


//   hRawYield_vsPt = (TH1D*)hist_unfold->Clone("hRawYield_vsPt");

//   for(int ibinPt = 1; ibinPt <= nbinpT; ibinPt++){
//     double dpT = hRawYield_vsPt->GetXaxis()->GetBinWidth(ibinPt);
//     double drapidity = 1.5; //1.5
//     double d2N_dpTdy = hRawYield_vsPt->GetBinContent(ibinPt) *1./dpT*1./drapidity;
//     hRawYield_vsPt->SetBinContent(ibinPt,d2N_dpTdy);
//     hRawYield_vsPt->SetBinError(ibinPt,hRawYield_vsPt->GetBinError(ibinPt) *1./dpT*1./drapidity);  // error on d2N_dpTdy 
//   }    

//   hRawYield_vsPt->Scale(1./SelectedEventCount);

//   // Lambda and AntiLambda Feeddown
//   TH1D *H1D_Feeddown_Correction = new TH1D("H1D_Feeddown_Correction", "H1D_Feeddown_Correction", nbinpT, pTbins);
//   if (ipart == 1 || ipart == 2) {
//     Get_FeedDown(H1D_Feeddown_Correction, ipart, file_O2Analysis, pTbins, nbinpT);
//     H1D_Feeddown_Correction->Scale(-1.);
//     hRawYield_vsPt->Add(H1D_Feeddown_Correction);
//     // cout << "test00000" << endl;
//     for(int ibinPt = 1; ibinPt <= nbinpT; ibinPt++){
//     cout << "H1D_Feeddown_Correction(" << ibinPt << ") = " << H1D_Feeddown_Correction->GetBinContent(ibinPt) << endl;
//     }
//   }
//   delete H1D_Feeddown_Correction;


//   cout << "--------------------------------------------------------------------- POST UNFOLDING TEST" << hist_unfold->GetEntries() << endl;
//   for(int ibinPt = 1; ibinPt <= nbinpT; ibinPt++){
//     cout << "truth_rebinned(" << ibinPt << ") = " << truth_rebinned->GetBinContent(ibinPt) << endl;
//     cout << "H2D_PtResponseMatrix(" << ibinPt << ") = " << H2D_PtResponseMatrix->GetBinContent(ibinPt,ibinPt) << endl;
//     cout << "hist_unfold(" << ibinPt << ") = " << hist_unfold->GetBinContent(ibinPt) << endl;
//     cout << "measured(" << ibinPt << ") = " << measured->GetBinContent(ibinPt) << endl;
//     cout << "measuredPostUnfoldPostCorrections(" << ibinPt << ") = " << hRawYield_vsPt->GetBinContent(ibinPt) << endl;
//   }
// }


// unfolding: do the refolding to check

// think about addding to o2physics the efficiecny based removal of tracks before jet finding


//   Marta thesis really detailed;
//   give details about choice of bins:

//   each slice of gen pt in det and fluct matrices are both normalised to 1 (proba distrib)

//   detector, both axes have same binning

//   combined matrix is rebinned to available data binning

//   generated axis should be bigger binning! ie fewer bins
//   marta does bin=5GeV for rec and 10GeV for gen; because this 10GeV will be the binning of the unfolded spectrum

//   The coarse response matrix is obtained by merging slices of the fine response matrix. Each slice in generated pT of the fine response matrix is weighted as function of a steeply falling spectrum.



//   TODOLIST : all the cent stuff for file_O2Analysis_ppSimDetectorEffect should be removed: it's pp, we're not looking at centrality




// ratio bin contents: 252.506
// ratio bin contents: 252.507
// ratio bin contents: 252.507
// ratio bin contents: 252.507
// ratio bin contents: 252.507
// ratio bin contents: 252.507
// ratio bin contents: 252.507
// ratio bin contents: 252.507
// ratio bin contents: 252.507
// ratio bin contents: 252.507


// whyyyyyy


// matrix multiplication: think about whether I really want to include over/underflows


// FUCK 
// Warning: RooUnfoldResponse measured X truth is 15 X 10, but matrix is 10 X 15
// except my TH2 has 15 lines and 10 columns ... so I guess that's not what roounfold wants

