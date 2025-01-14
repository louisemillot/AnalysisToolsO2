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
#include "TSVDUnfold.h"

//My Libraries
#include "./JetSpectrum_settings.h"
#include "../Settings/AxisTitles.h"
#include "../Settings/GlobalSettings.h"
#include "../Utilities/AnalysisUtilities.h"
#include "../Utilities/HistogramUtilities.h"
#include "../Utilities/HistogramPlotting.h"
#include "../Utilities/AnalysisUtilities.C" // bizarre but only including the .h fils doesn't work for the standard 'root macro.C+' method, I need to include the .C as well
#include "../Utilities/HistogramUtilities.C" // bizarre but only including the .h fils doesn't work for the standard 'root macro.C+' method, I need to include the .C as well
#include "../Utilities/HistogramPlotting.C" // bizarre but only including the .h fils doesn't work for the standard 'root macro.C+' method, I need to include the .C as well

#include<array>
#include <iomanip>
#include <sstream>
#include <string.h>
using namespace std;

// Misc utilities
void SetStyle(Bool_t graypalette=kFALSE);
void LoadLibs();

// Plot Utilities
TString contextDataset1D(int iDataset, float* variableRange, std::string options);
TString contextDatasetCompAndRadiusAndVarRange(float jetRadius, float* variableRange, std::string options);
TString contextPtRange(float* PtRange);
TString contextEtaRange(float* PtRange);
TString contextJetRadius(float jetRadius);
TString contextCentRange(float* centRange);

//////////// Pt Spectrum analysis functions
void Get_Pt_spectrum_bkgCorrected(TH1D* &H1D_jetPt_raw, int iDataset, int iRadius, float* centRange, std::string options);
void Get_Pt_spectrum_bkgCorrected_genBinning(TH1D* &H1D_jetPt_raw, int iDataset, int iRadius, float* centRange, std::string options);
void Get_Pt_spectrum_bkgCorrected_fineBinning(TH1D* &H1D_jetPt_raw, int iDataset, int iRadius, float* centRange, std::string options);
void Get_Pt_spectrum_mcp(TH1D* &H1D_jetPt_mcp, int iDataset, int iRadius, bool controlMC, std::string options);
void Get_Pt_spectrum_mcp_recBinning(TH1D* &H1D_jetPt_mcp, int iDataset, int iRadius, bool controlMC, std::string options);
void Get_Pt_spectrum_mcdMatched(TH1D* &H1D_jetPt_mcd, int iDataset, int iRadius, std::string options);
void Get_Pt_spectrum_mcpMatched(TH1D* &H1D_jetPt_mcp, int iDataset, int iRadius, std::string options);
std::pair<int, RooUnfold*> Get_Pt_spectrum_unfolded(TH1D* &H1D_jetPt_unfolded, int iDataset, int iRadius, float* centRange, int unfoldParameterBayes, std::string options);
void Get_Pt_spectrum_bkgCorrected_preWidthScalingAndEvtNorm(TH1D* &H1D_jetPt_raw, int iDataset, int iRadius, float* centRange, std::string options);
void Get_Pt_spectrum_bkgCorrected_genBinning_preWidthScalingAndEvtNorm(TH1D* &H1D_jetPt_raw, int iDataset, int iRadius, float* centRange, std::string options);
void Get_Pt_spectrum_bkgCorrected_fineBinning_preWidthScalingAndEvtNorm(TH1D* &H1D_jetPt_raw, int iDataset, int iRadius, float* centRange, std::string options);
void Get_Pt_spectrum_mcp_preWidthScalingAndEvtNorm(TH1D* &H1D_jetPt_mcp, int iDataset, int iRadius, bool controlMC, std::string options);
void Get_Pt_spectrum_mcp_recBinning_preWidthScalingAndEvtNorm(TH1D* &H1D_jetPt_mcp, int iDataset, int iRadius, bool controlMC, std::string options);
void Get_Pt_spectrum_mcdMatched_preWidthScalingAndEvtNorm(TH1D* &H1D_jetPt_mcd, int iDataset, int iRadius, std::string options);
void Get_Pt_spectrum_mcpMatched_preWidthScalingAndEvtNorm(TH1D* &H1D_jetPt_mcpMatched, int iDataset, int iRadius, std::string options);
void Get_Pt_spectrum_mcp_fineBinning_preWidthScalingAndEvtNorm(TH1D* &H1D_jetPt_mcp, int iDataset, int iRadius, bool controlMC, std::string options);
void Get_Pt_spectrum_mcd_fineBinning_preWidthScalingAndEvtNorm(TH1D* &H1D_jetPt_mcp, int iDataset, int iRadius, std::string options);
void Get_Pt_spectrum_mcd_recBinning_preWidthScalingAndEvtNorm(TH1D* &H1D_jetPt_mcp, int iDataset, int iRadius, std::string options);
void Get_Pt_spectrum_mcp_fineBinning(TH1D* &H1D_jetPt_mcp, int iDataset, int iRadius, bool controlMC, std::string options);
void Get_Pt_spectrum_mcd_fineBinning(TH1D* &H1D_jetPt_mcp, int iDataset, int iRadius, std::string options);
void Get_Pt_spectrum_mcd_recBinning(TH1D* &H1D_jetPt_mcp, int iDataset, int iRadius, std::string options);
std::pair<int, RooUnfold*> Get_Pt_spectrum_unfolded_preWidthScalingAndEvtNorm(TH1D* &H1D_jetPt_unfolded, int iDataset, int iRadius, float* centRange, int unfoldParameterBayes, std::string options);
void Get_Pt_spectrum_unfoldedThenRefolded_preWidthScalingAndEvtNorm(TH1D* &H1D_jetPt_mcp_unfoldedThenRefolded, int iDataset, int iRadius, float* centRange, int unfoldParameterBayes, std::string options);
void Get_Pt_spectrum_unfoldedThenRefolded(TH1D* &H1D_jetPt_mcp_unfoldedThenRefolded, int iDataset, int iRadius, float* centRange, int unfoldParameterBayes, std::string options);

void Get_PtResponseMatrix_Fluctuations(TH2D* &H2D_jetPtResponseMatrix_fluctuations, int iDataset, int iRadius, float* centRange);
void Get_PtResponseMatrix_detectorResponse(TH2D* &H2D_jetPtResponseMatrix_detectorResponse, int iDataset, int iRadius);
void Get_PtResponseMatrix_DetectorAndFluctuationsCombined(TH2D* &H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined, TH2D* H2D_jetPtResponseMatrix_detectorResponse, TH2D* H2D_jetPtResponseMatrix_fluctuations, int iDataset, int iRadius, float* centRange, std::string options);
void Get_PtResponseMatrix_DetectorAndFluctuationsCombined_fineBinning(TH2D* &H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined_fineBinning, TH2D* H2D_jetPtResponseMatrix_detectorResponse, TH2D* H2D_jetPtResponseMatrix_fluctuations, int iDataset, int iRadius, float* centRange, std::string options);
void ReweightResponseMatrixWithPrior(TH2D* &H2D_jetPtResponseMatrix, int iDataset, int iRadius, float* centRange, std::string options);


bool Get_Pt_JetEfficiency(TH1D* &H1D_jetEfficiency, bool* divideSuccess, int iDataset, int iRadius, float* centRange, std::string options);
void Get_ResponseMatrix_Pt_KinematicEffiency(TH1D* &H1D_kinematicEfficiency, TH2D* H2D_jetPtResponseMatrix, TString name_H1D_kinematicEfficiency, int iRadius);

void Draw_ResponseMatrices_Fluctuations(int iDataset, int iRadius);
void Draw_ResponseMatrices_detectorResponse(int iDataset, int iRadius);
void Draw_ResponseMatrices_DetectorAndFluctuationsCombined(int iDataset, int iRadius, std::string options);

void Draw_Pt_spectrum_raw(int iDataset, int iRadius, std::string options);
void Draw_Pt_spectrum_mcp(int iDataset, int iRadius, std::string options);
void Draw_Pt_spectrum_mcdMatched(int iDataset, int iRadius, std::string options);
void Draw_Pt_spectrum_unfolded_FluctResponseOnly(int iDataset, int iRadius, std::string options);
void Draw_Pt_spectrum_unfolded(int iDataset, int iRadius, int unfoldParameterBayes, std::string options);
void Draw_Pt_TestSpectrum_unfolded(int iDataset, int iRadius, std::string options);
void Draw_Pt_spectrum_unfolded_parameterVariation(int iDataset, int iRadius, int unfoldIterationMin, int unfoldIterationMax, int step, std::string options);

void Draw_Pt_efficiency_jets(int iDataset, int iRadius, std::string options);
void Draw_kinematicEfficiency(int iDataset, int iRadius, std::string options);
void Draw_FakeRatio(int iDataset, int iRadius, std::string options);

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

  // gathers the analysis options in a single char[]
  snprintf(optionsAnalysis, sizeof(optionsAnalysis), "%s,%s,%s,%s", mergingPrior, unfoldingPrior, unfoldingMethod, normMethod);
  cout << "Analysis options are: " << optionsAnalysis << endl;

  int iDataset = 0;
  int iRadius = 1;

  // find a way to input mcpPrior/mcdPrior and bayes/svd as a variables rather than typed out like this

  Draw_ResponseMatrices_Fluctuations(iDataset, iRadius);
  Draw_ResponseMatrices_detectorResponse(iDataset, iRadius);
  Draw_ResponseMatrices_DetectorAndFluctuationsCombined(iDataset, iRadius, optionsAnalysis);

  // // Draw_Pt_spectrum_unfolded_FluctResponseOnly(iDataset, iRadius, optionsAnalysis); // NOT FIXED YET - result meaningless

  // Draw_Pt_spectrum_raw(iDataset, iRadius, optionsAnalysis);
  // Draw_Pt_spectrum_mcp(iDataset, iRadius, optionsAnalysis);
  // Draw_Pt_spectrum_mcdMatched(iDataset, iRadius, optionsAnalysis);

  Draw_Pt_efficiency_jets(iDataset, iRadius, optionsAnalysis);

  Draw_kinematicEfficiency(iDataset, iRadius, optionsAnalysis);
  Draw_FakeRatio(iDataset, iRadius, optionsAnalysis);

  // int unfoldParameterBayes = 70;
  // Draw_Pt_spectrum_unfolded(iDataset, iRadius, unfoldParameterBayes, optionsAnalysis); // "evtNorm"

  int unfoldParameterBayesMin = 0;
  int unfoldParameterBayesMax = 100;
  int unfoldParameterBayesStep = 10;
  Draw_Pt_spectrum_unfolded_parameterVariation(iDataset, iRadius, unfoldParameterBayesMin, unfoldParameterBayesMax, unfoldParameterBayesStep, optionsAnalysis); // "evtNorm"
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

// TString contextDatasetRadiusCompAndVarRange(int iDataset, float* variableRange, std::string options){
//   TString texcontextDatasetRadiusCompAndVarRange;
//   if (options.find("pt") != std::string::npos) { //  || options.find("ratio") != NULL not sure why I had this here
//     texcontextDatasetRadiusCompAndVarRange = "#splitline{"+*texDatasetsComparisonCommonDenominator+" "+DatasetsNames[iDataset]+"}{#splitline{2023 QC}{"+contextPtRange(variableRange)+"}}";
//   }
//   if (options.find("eta") != std::string::npos) { //  || options.find("ratio") != NULL not sure why I had this here
//     texcontextDatasetRadiusCompAndVarRange = "#splitline{"+*texDatasetsComparisonCommonDenominator+" "+DatasetsNames[iDataset]+"}{#splitline{2023 QC}{"+contextEtaRange(variableRange)+"}}";
//   }
//   if (options.find("cent") != std::string::npos) { //  || options.find("ratio") != NULL not sure why I had this here
//     texcontextDatasetRadiusCompAndVarRange = "#splitline{"+*texDatasetsComparisonCommonDenominator+" "+DatasetsNames[iDataset]+"}{#splitline{2023 QC}{"+contextCentRange(variableRange)+"}}";
//   }

//   return texcontextDatasetRadiusCompAndVarRange;
// }

// TString contextDatasetCompAndRadiusAndVarRange(float jetRadius, float* variableRange, std::string options){
//   TString texcontextDatasetCompAndRadiusAndVarRange;
//   if (options.find("pt") != std::string::npos) { //  || options.find("ratio") != NULL not sure why I had this here
//     texcontextDatasetCompAndRadiusAndVarRange = "#splitline{"+*texDatasetsComparisonCommonDenominator+"}{#splitline{"+contextJetRadius(jetRadius)+"}{"+contextPtRange(variableRange)+"}}";
//   }
//   if (options.find("eta") != std::string::npos) { //  || options.find("ratio") != NULL not sure why I had this here
//     texcontextDatasetCompAndRadiusAndVarRange = "#splitline{"+*texDatasetsComparisonCommonDenominator+"}{#splitline{"+contextJetRadius(jetRadius)+"}{"+contextEtaRange(variableRange)+"}}";
//   }

//   return texcontextDatasetCompAndRadiusAndVarRange;
// }

// TString contextDatasetCompAndRadius(float jetRadius, std::string options){
//   TString texcontextDatasetCompAndRadius;
//   texcontextDatasetCompAndRadius = "#splitline{"+*texDatasetsComparisonCommonDenominator+"}{"+contextJetRadius(jetRadius)+"}";

//   return texcontextDatasetCompAndRadius;
// }

// TString contextDatasetComp(std::string options){
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


void Get_Pt_spectrum_bkgCorrected_preWidthScalingAndEvtNorm(TH1D* &H1D_jetPt, int iDataset, int iRadius, float* centRange, __attribute__ ((unused)) std::string options) {
  TH3D* H3D_jetRjetPtcolCent;
  TH1D* H1D_jetPt_defaultBin;
  // TH1D* H1D_jetPt_raw[nRadius];

  if (isPbPb == true){
  H3D_jetRjetPtcolCent = (TH3D*)((TH3D*)file_O2Analysis_list[iDataset]->Get(analysisWorkflowData+"/h3_jet_r_jet_pt_centrality_rhoareasubtracted"))->Clone("Get_Pt_spectrum_bkgCorrected"+Datasets[iDataset]+Form("%.1f", centRange[0])+"<cent<"+Form("%.1f", centRange[1]));
  } else {
    if (!controlMC) {
      H3D_jetRjetPtcolCent = (TH3D*)((TH3D*)file_O2Analysis_list[iDataset]->Get(analysisWorkflowData+"/h3_jet_r_jet_pt_centrality"))->Clone("Get_Pt_spectrum_bkgCorrected"+Datasets[iDataset]+Form("%.1f", centRange[0])+"<cent<"+Form("%.1f", centRange[1]));
    } else {
      H3D_jetRjetPtcolCent = (TH3D*)((TH3D*)file_O2Analysis_ppSimDetectorEffect_unfoldingControl->Get(analysisWorkflowData+"/h3_jet_r_jet_pt_centrality"))->Clone("Get_Pt_spectrum_bkgCorrected"+Datasets[iDataset]+Form("%.1f", centRange[0])+"<cent<"+Form("%.1f", centRange[1]));
    }
  }
  H3D_jetRjetPtcolCent->Sumw2();
  int binCentEdges[2] = {H3D_jetRjetPtcolCent->GetZaxis()->FindBin(centRange[0]), H3D_jetRjetPtcolCent->GetZaxis()->FindBin(centRange[1])};
  if (binCentEdges[0] == 0) 
    cout << "WARNING: Get_Pt_spectrum_bkgCorrected is counting the underflow with the chosen centRange" << endl;
  if (binCentEdges[1] == H3D_jetRjetPtcolCent->GetZaxis()->GetNbins()+1) 
    cout << "WARNING: Get_Pt_spectrum_bkgCorrected is counting the overflow with the chosen centRange" << endl;
  if (isPbPb == false){
    binCentEdges[0] = 0;
    binCentEdges[1] = -1;
  }

  int ibinJetRadius = H3D_jetRjetPtcolCent->GetXaxis()->FindBin(arrayRadius[iRadius]+GLOBAL_epsilon);

  H1D_jetPt_defaultBin = (TH1D*)H3D_jetRjetPtcolCent->ProjectionY("jetPt_bkgCorrected_"+RadiusLegend[iRadius]+Form("%.1f", centRange[0])+"<cent<"+Form("%.1f", centRange[1]), ibinJetRadius, ibinJetRadius, binCentEdges[0], binCentEdges[1], "e");
  H1D_jetPt = (TH1D*)H1D_jetPt_defaultBin->Rebin(nBinPtJetsRec[iRadius],"jetPt_bkgCorrected_rebinned_"+RadiusLegend[iRadius]+Form("%.1f", centRange[0])+"<cent<"+Form("%.1f", centRange[1]), ptBinsJetsRec[iRadius]);

  if (doWidthScalingEarly) {
    TransformRawHistToYield(H1D_jetPt);
  }

}


void Get_Pt_spectrum_bkgCorrected_genBinning_preWidthScalingAndEvtNorm(TH1D* &H1D_jetPt, int iDataset, int iRadius, float* centRange, __attribute__ ((unused)) std::string options) {
  TH3D* H3D_jetRjetPtcolCent;
  TH1D* H1D_jetPt_defaultBin;
  // TH1D* H1D_jetPt_raw[nRadius];
  if (isPbPb == true){
      H3D_jetRjetPtcolCent = (TH3D*)((TH3D*)file_O2Analysis_list[iDataset]->Get(analysisWorkflowData+"/h3_jet_r_jet_pt_centrality_rhoareasubtracted"))->Clone("Get_Pt_spectrum_bkgCorrected_genBinning"+Datasets[iDataset]+Form("%.1f", centRange[0])+"<cent<"+Form("%.1f", centRange[1]));
  } else {
    if (!controlMC) {
      H3D_jetRjetPtcolCent = (TH3D*)((TH3D*)file_O2Analysis_list[iDataset]->Get(analysisWorkflowData+"/h3_jet_r_jet_pt_centrality"))->Clone("Get_Pt_spectrum_bkgCorrected"+Datasets[iDataset]+Form("%.1f", centRange[0])+"<cent<"+Form("%.1f", centRange[1]));
    } else {
      H3D_jetRjetPtcolCent = (TH3D*)((TH3D*)file_O2Analysis_ppSimDetectorEffect_unfoldingControl->Get(analysisWorkflowData+"/h3_jet_r_jet_pt_centrality"))->Clone("Get_Pt_spectrum_bkgCorrected"+Datasets[iDataset]+Form("%.1f", centRange[0])+"<cent<"+Form("%.1f", centRange[1]));
    }
  }
  H3D_jetRjetPtcolCent->Sumw2();
  int binCentEdges[2] = {H3D_jetRjetPtcolCent->GetZaxis()->FindBin(centRange[0]), H3D_jetRjetPtcolCent->GetZaxis()->FindBin(centRange[1])};
  if (binCentEdges[0] == 0) 
    cout << "WARNING: Get_Pt_spectrum_bkgCorrected_genBinning is counting the underflow with the chosen centRange" << endl;
  if (binCentEdges[1] == H3D_jetRjetPtcolCent->GetZaxis()->GetNbins()+1) 
    cout << "WARNING: Get_Pt_spectrum_bkgCorrected_genBinning is counting the overflow with the chosen centRange" << endl;
  if (isPbPb == false){
    binCentEdges[0] = 0;
    binCentEdges[1] = -1;
  }

  int ibinJetRadius = H3D_jetRjetPtcolCent->GetXaxis()->FindBin(arrayRadius[iRadius]+GLOBAL_epsilon);

  H1D_jetPt_defaultBin = (TH1D*)H3D_jetRjetPtcolCent->ProjectionY("jetPt_bkgCorrected_genBinning_"+RadiusLegend[iRadius]+Form("%.1f", centRange[0])+"<cent<"+Form("%.1f", centRange[1]), ibinJetRadius, ibinJetRadius, binCentEdges[0], binCentEdges[1], "e");
  H1D_jetPt = (TH1D*)H1D_jetPt_defaultBin->Rebin(nBinPtJetsGen[iRadius],"jetPt_bkgCorrected_rebinned_genBinning_"+RadiusLegend[iRadius]+Form("%.1f", centRange[0])+"<cent<"+Form("%.1f", centRange[1]), ptBinsJetsGen[iRadius]);

  if (doWidthScalingEarly) {
    TransformRawHistToYield(H1D_jetPt);
  }

}

void Get_Pt_spectrum_bkgCorrected_fineBinning_preWidthScalingAndEvtNorm(TH1D* &H1D_jetPt, int iDataset, int iRadius, float* centRange, __attribute__ ((unused)) std::string options) {
  TH3D* H3D_jetRjetPtcolCent;
  TH1D* H1D_jetPt_defaultBin;
  // TH1D* H1D_jetPt_raw[nRadius];

  if (isPbPb == true){
      H3D_jetRjetPtcolCent = (TH3D*)((TH3D*)file_O2Analysis_list[iDataset]->Get(analysisWorkflowData+"/h3_jet_r_jet_pt_centrality_rhoareasubtracted"))->Clone("Get_Pt_spectrum_bkgCorrected_fineBinning"+Datasets[iDataset]+Form("%.1f", centRange[0])+"<cent<"+Form("%.1f", centRange[1]));
  } else {
    if (!controlMC) {
      H3D_jetRjetPtcolCent = (TH3D*)((TH3D*)file_O2Analysis_list[iDataset]->Get(analysisWorkflowData+"/h3_jet_r_jet_pt_centrality"))->Clone("Get_Pt_spectrum_bkgCorrected"+Datasets[iDataset]+Form("%.1f", centRange[0])+"<cent<"+Form("%.1f", centRange[1]));
    } else {
      H3D_jetRjetPtcolCent = (TH3D*)((TH3D*)file_O2Analysis_ppSimDetectorEffect_unfoldingControl->Get(analysisWorkflowData+"/h3_jet_r_jet_pt_centrality"))->Clone("Get_Pt_spectrum_bkgCorrected"+Datasets[iDataset]+Form("%.1f", centRange[0])+"<cent<"+Form("%.1f", centRange[1]));
    }
  }
  H3D_jetRjetPtcolCent->Sumw2();
  int binCentEdges[2] = {H3D_jetRjetPtcolCent->GetZaxis()->FindBin(centRange[0]), H3D_jetRjetPtcolCent->GetZaxis()->FindBin(centRange[1])};
  if (binCentEdges[0] == 0) 
    cout << "WARNING: Get_Pt_spectrum_bkgCorrected_fineBinning is counting the underflow with the chosen centRange" << endl;
  if (binCentEdges[1] == H3D_jetRjetPtcolCent->GetZaxis()->GetNbins()+1) 
    cout << "WARNING: Get_Pt_spectrum_bkgCorrected_fineBinning is counting the overflow with the chosen centRange" << endl;
  if (isPbPb == false){
    binCentEdges[0] = 0;
    binCentEdges[1] = -1;
  }

  int ibinJetRadius = H3D_jetRjetPtcolCent->GetXaxis()->FindBin(arrayRadius[iRadius]+GLOBAL_epsilon);

  H1D_jetPt_defaultBin = (TH1D*)H3D_jetRjetPtcolCent->ProjectionY("jetPt_bkgCorrected_fineBinning_"+RadiusLegend[iRadius]+Form("%.1f", centRange[0])+"<cent<"+Form("%.1f", centRange[1]), ibinJetRadius, ibinJetRadius, binCentEdges[0], binCentEdges[1], "e");
  H1D_jetPt = (TH1D*)H1D_jetPt_defaultBin->Rebin(nBinPtJetsFine[iRadius],"jetPt_bkgCorrected_rebinned_fineBinning_"+RadiusLegend[iRadius]+Form("%.1f", centRange[0])+"<cent<"+Form("%.1f", centRange[1]), ptBinsJetsFine[iRadius]);

  if (doWidthScalingEarly) {
    TransformRawHistToYield(H1D_jetPt);
  }

}

void Get_Pt_spectrum_mcp_preWidthScalingAndEvtNorm(TH1D* &H1D_jetPt_mcp, int iDataset, int iRadius, bool fcontrolMC, __attribute__ ((unused)) std::string options) {
  TH3D* H3D_jetRjetPtjetEta;
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
  if (!fcontrolMC){
    H3D_jetRjetPtjetEta = (TH3D*)((TH3D*)file_O2Analysis_ppSimDetectorEffect->Get(analysisWorkflowMC+"/h3_jet_r_part_jet_pt_part_jet_eta_part"))->Clone("Get_Pt_spectrum_mcp"+Datasets[iDataset]);
  } else {
    H3D_jetRjetPtjetEta = (TH3D*)((TH3D*)file_O2Analysis_ppSimDetectorEffect_unfoldingControl->Get(analysisWorkflowMC+"/h3_jet_r_part_jet_pt_part_jet_eta_part"))->Clone("Get_Pt_spectrum_mcp_unfoldingControl"+Datasets[iDataset]);
  }
  H3D_jetRjetPtjetEta->Sumw2();

  int ibinJetRadius = H3D_jetRjetPtjetEta->GetXaxis()->FindBin(arrayRadius[iRadius]+GLOBAL_epsilon);
  int binEtaEdges[2] = {H3D_jetRjetPtjetEta->GetZaxis()->FindBin(etaAnalysisRange[0]), H3D_jetRjetPtjetEta->GetZaxis()->FindBin(etaAnalysisRange[1])-1};

  H1D_jetPt_mcp_defaultBin = (TH1D*)H3D_jetRjetPtjetEta->ProjectionY("jetPt_mcp_"+RadiusLegend[iRadius]+Datasets[iDataset], ibinJetRadius, ibinJetRadius, binEtaEdges[0], binEtaEdges[1], "e");

  H1D_jetPt_mcp = (TH1D*)H1D_jetPt_mcp_defaultBin->Rebin(nBinPtJetsGen[iRadius],"jetPt_mcp_rebinned_"+RadiusLegend[iRadius]+Datasets[iDataset], ptBinsJetsGen[iRadius]);

  if (doWidthScalingEarly) {
    TransformRawHistToYield(H1D_jetPt_mcp);
  }

}

void Get_Pt_spectrum_mcp_fineBinning_preWidthScalingAndEvtNorm(TH1D* &H1D_jetPt_mcp, int iDataset, int iRadius, bool fcontrolMC, __attribute__ ((unused)) std::string options) {
  TH3D* H3D_jetRjetPtjetEta;
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
  if (!fcontrolMC){
    H3D_jetRjetPtjetEta = (TH3D*)((TH3D*)file_O2Analysis_ppSimDetectorEffect->Get(analysisWorkflowMC+"/h3_jet_r_part_jet_pt_part_jet_eta_part"))->Clone("Get_Pt_spectrum_mcp_fineBinning_preWidthScalingAndEvtNorm"+Datasets[iDataset]);
  } else {
    H3D_jetRjetPtjetEta = (TH3D*)((TH3D*)file_O2Analysis_ppSimDetectorEffect_unfoldingControl->Get(analysisWorkflowMC+"/h3_jet_r_part_jet_pt_part_jet_eta_part"))->Clone("Get_Pt_spectrum_mcp_fineBinning_preWidthScalingAndEvtNorm_unfoldingControl"+Datasets[iDataset]);
  }
  H3D_jetRjetPtjetEta->Sumw2();

  int ibinJetRadius = H3D_jetRjetPtjetEta->GetXaxis()->FindBin(arrayRadius[iRadius]+GLOBAL_epsilon);
  int binEtaEdges[2] = {H3D_jetRjetPtjetEta->GetZaxis()->FindBin(etaAnalysisRange[0]), H3D_jetRjetPtjetEta->GetZaxis()->FindBin(etaAnalysisRange[1])-1};

  H1D_jetPt_mcp_defaultBin = (TH1D*)H3D_jetRjetPtjetEta->ProjectionY("jetPt_mcp_fineBinning_"+RadiusLegend[iRadius]+Datasets[iDataset], ibinJetRadius, ibinJetRadius, binEtaEdges[0], binEtaEdges[1], "e");

  H1D_jetPt_mcp = (TH1D*)H1D_jetPt_mcp_defaultBin->Rebin(nBinPtJetsFine[iRadius],"jetPt_mcp_rebinned_fineBinning_"+RadiusLegend[iRadius]+Datasets[iDataset], ptBinsJetsFine[iRadius]);

  if (doWidthScalingEarly) {
    TransformRawHistToYield(H1D_jetPt_mcp);
  }
}

void Get_Pt_spectrum_mcd_fineBinning_preWidthScalingAndEvtNorm(TH1D* &H1D_jetPt_mcd, int iDataset, int iRadius, __attribute__ ((unused)) std::string options) {
  TH3D* H3D_jetRjetPtjetEta;
  TH1D* H1D_jetPt_mcd_defaultBin;
  // TH1D* H1D_jetPt_raw[nRadius];

  //////// without centrality; sadly we don't have a mcd version that gives centrality info, so gotta use one with all centralities (though one can edit the cent window in the jetfinderqa options)
  H3D_jetRjetPtjetEta = (TH3D*)((TH3D*)file_O2Analysis_ppSimDetectorEffect->Get(analysisWorkflowMC+"/h3_jet_r_jet_pt_jet_eta"))->Clone("Get_Pt_spectrum_mcd_fineBinning_preWidthScalingAndEvtNorm"+Datasets[iDataset]);
  H3D_jetRjetPtjetEta->Sumw2();

  int ibinJetRadius = H3D_jetRjetPtjetEta->GetXaxis()->FindBin(arrayRadius[iRadius]+GLOBAL_epsilon);
  int binEtaEdges[2] = {H3D_jetRjetPtjetEta->GetZaxis()->FindBin(etaAnalysisRange[0]), H3D_jetRjetPtjetEta->GetZaxis()->FindBin(etaAnalysisRange[1])-1};

  H1D_jetPt_mcd_defaultBin = (TH1D*)H3D_jetRjetPtjetEta->ProjectionY("jetPt_mcd_fineBinning_"+RadiusLegend[iRadius]+Datasets[iDataset], ibinJetRadius, ibinJetRadius, binEtaEdges[0], binEtaEdges[1], "e");

  H1D_jetPt_mcd = (TH1D*)H1D_jetPt_mcd_defaultBin->Rebin(nBinPtJetsFine[iRadius],"jetPt_mcd_rebinned_fineBinning_"+RadiusLegend[iRadius]+Datasets[iDataset], ptBinsJetsFine[iRadius]);

  if (doWidthScalingEarly) {
    TransformRawHistToYield(H1D_jetPt_mcd);
  }

}

void Get_Pt_spectrum_mcd_recBinning_preWidthScalingAndEvtNorm(TH1D* &H1D_jetPt_mcd, int iDataset, int iRadius, __attribute__ ((unused)) std::string options) {
  TH3D* H3D_jetRjetPtjetEta;
  TH1D* H1D_jetPt_mcd_defaultBin;
  // TH1D* H1D_jetPt_raw[nRadius];

  //////// without centrality; sadly we don't have a mcd version that gives centrality info, so gotta use one with all centralities (though one can edit the cent window in the jetfinderqa options)
  H3D_jetRjetPtjetEta = (TH3D*)((TH3D*)file_O2Analysis_ppSimDetectorEffect->Get(analysisWorkflowMC+"/h3_jet_r_jet_pt_jet_eta"))->Clone("Get_Pt_spectrum_mcd_recBinning_preWidthScalingAndEvtNorm"+Datasets[iDataset]);
  H3D_jetRjetPtjetEta->Sumw2();

  int ibinJetRadius = H3D_jetRjetPtjetEta->GetXaxis()->FindBin(arrayRadius[iRadius]+GLOBAL_epsilon);
  int binEtaEdges[2] = {H3D_jetRjetPtjetEta->GetZaxis()->FindBin(etaAnalysisRange[0]), H3D_jetRjetPtjetEta->GetZaxis()->FindBin(etaAnalysisRange[1])-1};

  H1D_jetPt_mcd_defaultBin = (TH1D*)H3D_jetRjetPtjetEta->ProjectionY("jetPt_mcd_recBinning_"+RadiusLegend[iRadius]+Datasets[iDataset], ibinJetRadius, ibinJetRadius, binEtaEdges[0], binEtaEdges[1], "e");

  H1D_jetPt_mcd = (TH1D*)H1D_jetPt_mcd_defaultBin->Rebin(nBinPtJetsRec[iRadius],"jetPt_mcd_rebinned_recBinning_"+RadiusLegend[iRadius]+Datasets[iDataset], ptBinsJetsRec[iRadius]);

  if (doWidthScalingEarly) {
    TransformRawHistToYield(H1D_jetPt_mcd);
  }

}

void Get_Pt_spectrum_mcd_preWidthScalingAndEvtNorm(TH1D* &H1D_jetPt_mcd, int iDataset, int iRadius, __attribute__ ((unused)) std::string options) {
  TH3D* H3D_jetRjetPtjetEta;
  TH1D* H1D_jetPt_mcd_defaultBin;
  // TH1D* H1D_jetPt_raw[nRadius];

  //////// without centrality; sadly we don't have a mcd version that gives centrality info, so gotta use one with all centralities (though one can edit the cent window in the jetfinderqa options)
  H3D_jetRjetPtjetEta = (TH3D*)((TH3D*)file_O2Analysis_ppSimDetectorEffect->Get(analysisWorkflowMC+"/h3_jet_r_jet_pt_jet_eta"))->Clone("Get_Pt_spectrum_mcd_preWidthScalingAndEvtNorm"+Datasets[iDataset]);
  H3D_jetRjetPtjetEta->Sumw2();

  int ibinJetRadius = H3D_jetRjetPtjetEta->GetXaxis()->FindBin(arrayRadius[iRadius]+GLOBAL_epsilon);
  int binEtaEdges[2] = {H3D_jetRjetPtjetEta->GetZaxis()->FindBin(etaAnalysisRange[0]), H3D_jetRjetPtjetEta->GetZaxis()->FindBin(etaAnalysisRange[1])-1};

  H1D_jetPt_mcd_defaultBin = (TH1D*)H3D_jetRjetPtjetEta->ProjectionY("jetPt_mcd_genBinning_"+RadiusLegend[iRadius]+Datasets[iDataset], ibinJetRadius, ibinJetRadius, binEtaEdges[0], binEtaEdges[1], "e");

  H1D_jetPt_mcd = (TH1D*)H1D_jetPt_mcd_defaultBin->Rebin(nBinPtJetsGen[iRadius],"jetPt_mcd_rebinned_genBinning_"+RadiusLegend[iRadius]+Datasets[iDataset], ptBinsJetsGen[iRadius]);

  if (doWidthScalingEarly) {
    TransformRawHistToYield(H1D_jetPt_mcd);
  }
}


void Get_Pt_spectrum_mcp_recBinning_preWidthScalingAndEvtNorm(TH1D* &H1D_jetPt_mcp, int iDataset, int iRadius, bool fcontrolMC, __attribute__ ((unused)) std::string options) {
  TH3D* H3D_jetRjetPtjetEta;
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
  if (!fcontrolMC){
    H3D_jetRjetPtjetEta = (TH3D*)((TH3D*)file_O2Analysis_ppSimDetectorEffect->Get(analysisWorkflowMC+"/h3_jet_r_part_jet_pt_part_jet_eta_part"))->Clone("Get_Pt_spectrum_mcp_recBinning_preWidthScalingAndEvtNorm"+Datasets[iDataset]);
  } else {
    H3D_jetRjetPtjetEta = (TH3D*)((TH3D*)file_O2Analysis_ppSimDetectorEffect_unfoldingControl->Get(analysisWorkflowMC+"/h3_jet_r_part_jet_pt_part_jet_eta_part"))->Clone("Get_Pt_spectrum_mcp_recBinning_preWidthScalingAndEvtNorm_unfoldingControl"+Datasets[iDataset]);
  }
  H3D_jetRjetPtjetEta->Sumw2();

  int ibinJetRadius = H3D_jetRjetPtjetEta->GetXaxis()->FindBin(arrayRadius[iRadius]+GLOBAL_epsilon);
  int binEtaEdges[2] = {H3D_jetRjetPtjetEta->GetZaxis()->FindBin(etaAnalysisRange[0]), H3D_jetRjetPtjetEta->GetZaxis()->FindBin(etaAnalysisRange[1])-1};

  H1D_jetPt_mcp_defaultBin = (TH1D*)H3D_jetRjetPtjetEta->ProjectionY("jetPt_mcp_rebinRec_"+RadiusLegend[iRadius]+Datasets[iDataset], ibinJetRadius, ibinJetRadius, binEtaEdges[0], binEtaEdges[1], "e");

  H1D_jetPt_mcp = (TH1D*)H1D_jetPt_mcp_defaultBin->Rebin(nBinPtJetsRec[iRadius],"jetPt_mcp_rebinnedRec_"+RadiusLegend[iRadius]+Datasets[iDataset], ptBinsJetsRec[iRadius]);


  if (doWidthScalingEarly) {
    TransformRawHistToYield(H1D_jetPt_mcp);
  }

}

void Get_Pt_spectrum_mcpMatched_preWidthScalingAndEvtNorm(TH1D* &H1D_jetPt_mcpMatched, int iDataset, int iRadius, __attribute__ ((unused)) std::string options) {
  TH3D* H3D_jetRjetPtTagjetPtBase;
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
  H3D_jetRjetPtTagjetPtBase = (TH3D*)((TH3D*)file_O2Analysis_ppSimDetectorEffect->Get(analysisWorkflowMC+"/h3_jet_r_jet_pt_tag_jet_pt_base_matchedgeo"))->Clone("Get_Pt_spectrum_mcpMatched"+Datasets[iDataset]);// base is mcd in jetfinderQA as of 06/2024, thus tag is mcp, and so hist is (x=r, y=mcp, z=mcd)
  H3D_jetRjetPtTagjetPtBase->Sumw2();


  int ibinJetRadius = H3D_jetRjetPtTagjetPtBase->GetXaxis()->FindBin(arrayRadius[iRadius]+GLOBAL_epsilon);

  H1D_jetPt_mcpMatched_defaultBin = (TH1D*)H3D_jetRjetPtTagjetPtBase->ProjectionY("jetPt_mcpMatched_"+RadiusLegend[iRadius]+Datasets[iDataset], ibinJetRadius, ibinJetRadius, 1, H3D_jetRjetPtTagjetPtBase->GetNbinsZ(), "e");

  // should be rebinned using rec bins, because in the unfolding it'll be given to RooUnfoldResponse(mcpMatched, mcp, response) to get the efficiency 
  H1D_jetPt_mcpMatched = (TH1D*)H1D_jetPt_mcpMatched_defaultBin->Rebin(nBinPtJetsGen[iRadius],"jetPt_mcpMatched_rebinned_"+RadiusLegend[iRadius]+Datasets[iDataset], ptBinsJetsGen[iRadius]);





  if (doWidthScalingEarly) {
    TransformRawHistToYield(H1D_jetPt_mcpMatched);
  }



}

void Get_Pt_spectrum_mcpMatched_fineBinning_preWidthScalingAndEvtNorm(TH1D* &H1D_jetPt_mcpMatched, int iDataset, int iRadius, __attribute__ ((unused)) std::string options) {
  TH3D* H3D_jetRjetPtTagjetPtBase;
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
  H3D_jetRjetPtTagjetPtBase = (TH3D*)((TH3D*)file_O2Analysis_ppSimDetectorEffect->Get(analysisWorkflowMC+"/h3_jet_r_jet_pt_tag_jet_pt_base_matchedgeo"))->Clone("Get_Pt_spectrum_mcpMatched"+Datasets[iDataset]);// base is mcd in jetfinderQA as of 06/2024, thus tag is mcp, and so hist is (x=r, y=mcp, z=mcd)
  H3D_jetRjetPtTagjetPtBase->Sumw2();


  int ibinJetRadius = H3D_jetRjetPtTagjetPtBase->GetXaxis()->FindBin(arrayRadius[iRadius]+GLOBAL_epsilon);

  H1D_jetPt_mcpMatched_defaultBin = (TH1D*)H3D_jetRjetPtTagjetPtBase->ProjectionY("jetPt_mcpMatched_"+RadiusLegend[iRadius]+Datasets[iDataset], ibinJetRadius, ibinJetRadius, 1, H3D_jetRjetPtTagjetPtBase->GetNbinsZ(), "e");

  // should be rebinned using rec bins, because in the unfolding it'll be given to RooUnfoldResponse(mcpMatched, mcp, response) to get the efficiency 
  H1D_jetPt_mcpMatched = (TH1D*)H1D_jetPt_mcpMatched_defaultBin->Rebin(nBinPtJetsFine[iRadius],"jetPt_mcpMatched_rebinned_"+RadiusLegend[iRadius]+Datasets[iDataset], ptBinsJetsFine[iRadius]);

  if (doWidthScalingEarly) {
    TransformRawHistToYield(H1D_jetPt_mcpMatched);
  }

}



void Get_Pt_spectrum_mcdMatched_preWidthScalingAndEvtNorm(TH1D* &H1D_jetPt_mcdMatched, int iDataset, int iRadius, __attribute__ ((unused)) std::string options) {
  TH3D* H3D_jetRjetPtTagjetPtBase;
  TH1D* H1D_jetPt_mcdMatched_defaultBin;

  H3D_jetRjetPtTagjetPtBase = (TH3D*)((TH3D*)file_O2Analysis_ppSimDetectorEffect->Get(analysisWorkflowMC+"/h3_jet_r_jet_pt_tag_jet_pt_base_matchedgeo"))->Clone("Get_Pt_spectrum_mcdMatched"+Datasets[iDataset]);// base is mcd in jetfinderQA as of 06/2024, thus tag is mcp, and so hist is (x=r, y=mcp, z=mcd)
  H3D_jetRjetPtTagjetPtBase->Sumw2();

  int ibinJetRadius = H3D_jetRjetPtTagjetPtBase->GetXaxis()->FindBin(arrayRadius[iRadius]+GLOBAL_epsilon);

  H1D_jetPt_mcdMatched_defaultBin = (TH1D*)H3D_jetRjetPtTagjetPtBase->ProjectionZ("jetPt_mcdMatched_"+RadiusLegend[iRadius]+Datasets[iDataset], ibinJetRadius, ibinJetRadius, 1, H3D_jetRjetPtTagjetPtBase->GetNbinsZ(), "e");

  H1D_jetPt_mcdMatched = (TH1D*)H1D_jetPt_mcdMatched_defaultBin->Rebin(nBinPtJetsGen[iRadius],"jetPt_mcdMatched_rebinned_"+RadiusLegend[iRadius]+Datasets[iDataset], ptBinsJetsGen[iRadius]);

  if (doWidthScalingEarly) {
    TransformRawHistToYield(H1D_jetPt_mcdMatched);
  }
}

void Get_Pt_spectrum_mcdMatched_recBinning_preWidthScalingAndEvtNorm(TH1D* &H1D_jetPt_mcdMatched, int iDataset, int iRadius, __attribute__ ((unused)) std::string options) {
  TH3D* H3D_jetRjetPtTagjetPtBase;
  TH1D* H1D_jetPt_mcdMatched_defaultBin;

  H3D_jetRjetPtTagjetPtBase = (TH3D*)((TH3D*)file_O2Analysis_ppSimDetectorEffect->Get(analysisWorkflowMC+"/h3_jet_r_jet_pt_tag_jet_pt_base_matchedgeo"))->Clone("Get_Pt_spectrum_mcdMatched_recBinning"+Datasets[iDataset]);// base is mcd in jetfinderQA as of 06/2024, thus tag is mcp, and so hist is (x=r, y=mcp, z=mcd)
  H3D_jetRjetPtTagjetPtBase->Sumw2();

  int ibinJetRadius = H3D_jetRjetPtTagjetPtBase->GetXaxis()->FindBin(arrayRadius[iRadius]+GLOBAL_epsilon);

  H1D_jetPt_mcdMatched_defaultBin = (TH1D*)H3D_jetRjetPtTagjetPtBase->ProjectionZ("jetPt_mcdMatched_recBinning"+RadiusLegend[iRadius]+Datasets[iDataset], ibinJetRadius, ibinJetRadius, 1, H3D_jetRjetPtTagjetPtBase->GetNbinsZ(), "e");

  H1D_jetPt_mcdMatched = (TH1D*)H1D_jetPt_mcdMatched_defaultBin->Rebin(nBinPtJetsRec[iRadius],"jetPt_mcdMatched_recBinning_rebinned_"+RadiusLegend[iRadius]+Datasets[iDataset], ptBinsJetsRec[iRadius]);

  if (doWidthScalingEarly) {
    TransformRawHistToYield(H1D_jetPt_mcdMatched);
  }
}

void Get_Pt_spectrum_mcdMatched_fineBinning_preWidthScalingAndEvtNorm(TH1D* &H1D_jetPt_mcdMatched, int iDataset, int iRadius, __attribute__ ((unused)) std::string options) {
  TH3D* H3D_jetRjetPtTagjetPtBase;
  TH1D* H1D_jetPt_mcdMatched_defaultBin;

  H3D_jetRjetPtTagjetPtBase = (TH3D*)((TH3D*)file_O2Analysis_ppSimDetectorEffect->Get(analysisWorkflowMC+"/h3_jet_r_jet_pt_tag_jet_pt_base_matchedgeo"))->Clone("Get_Pt_spectrum_mcdMatched_fineBinning"+Datasets[iDataset]);// base is mcd in jetfinderQA as of 06/2024, thus tag is mcp, and so hist is (x=r, y=mcp, z=mcd)
  H3D_jetRjetPtTagjetPtBase->Sumw2();

  int ibinJetRadius = H3D_jetRjetPtTagjetPtBase->GetXaxis()->FindBin(arrayRadius[iRadius]+GLOBAL_epsilon);

  H1D_jetPt_mcdMatched_defaultBin = (TH1D*)H3D_jetRjetPtTagjetPtBase->ProjectionZ("jetPt_mcdMatched_fineBinning_"+RadiusLegend[iRadius]+Datasets[iDataset], ibinJetRadius, ibinJetRadius, 1, H3D_jetRjetPtTagjetPtBase->GetNbinsZ(), "e");

  H1D_jetPt_mcdMatched = (TH1D*)H1D_jetPt_mcdMatched_defaultBin->Rebin(nBinPtJetsFine[iRadius],"jetPt_mcdMatched_fineBinning_rebinned_"+RadiusLegend[iRadius]+Datasets[iDataset], ptBinsJetsFine[iRadius]);


  if (doWidthScalingEarly) {
    TransformRawHistToYield(H1D_jetPt_mcdMatched);
  }
}

void Get_Pt_spectrum_bkgCorrected_preWidthScaling(TH1D* &H1D_jetPt, int iDataset, int iRadius, float* centRange, std::string options) {
  Get_Pt_spectrum_bkgCorrected_preWidthScalingAndEvtNorm(H1D_jetPt, iDataset, iRadius, centRange, options);


  if (options.find("evtNorm") != std::string::npos && doEvtNorm) {
    if (isPbPb) {
      NormaliseRawHistToNEvents(H1D_jetPt, GetNEventsSelectedCentrality_JetFramework(file_O2Analysis_list[iDataset], centRange[0], centRange[1], trainId));
    } else {
      if (!controlMC) {
        NormaliseRawHistToNEvents(H1D_jetPt, GetNEventsSelectedCentrality_JetFramework(file_O2Analysis_list[iDataset], -999, 999, trainId));
      } else {
        NormaliseRawHistToNEvents(H1D_jetPt, GetNEventsSelectedCentrality_JetFramework(file_O2Analysis_ppSimDetectorEffect_unfoldingControl, -999, 999, trainId));
      }
    }
  }
  if (options.find("entriesNorm") != std::string::npos) {
    NormaliseYieldToNEntries(H1D_jetPt);
  }
}
void Get_Pt_spectrum_bkgCorrected(TH1D* &H1D_jetPt, int iDataset, int iRadius, float* centRange, std::string options) {
  Get_Pt_spectrum_bkgCorrected_preWidthScaling(H1D_jetPt, iDataset, iRadius, centRange, options);

  if (doWidthScalingAtEnd) {
    TransformRawHistToYield(H1D_jetPt);
  }
}

void Get_Pt_spectrum_mcp_preWidthScaling(TH1D* &H1D_jetPt_mcp, int iDataset, int iRadius, bool fcontrolMC, std::string options) {
  Get_Pt_spectrum_mcp_preWidthScalingAndEvtNorm(H1D_jetPt_mcp, iDataset, iRadius, fcontrolMC, options);

  if (options.find("evtNorm") != std::string::npos && doEvtNorm) {
    if (isPbPb == true || ppMcIsWeighted) {
      NormaliseRawHistToNEvents(H1D_jetPt_mcp, GetNEventsSelected_JetFramework_weighted(file_O2Analysis_ppSimDetectorEffect));
    } else {
      if (!fcontrolMC) {
        NormaliseRawHistToNEvents(H1D_jetPt_mcp, GetNEventsSelected_JetFramework(file_O2Analysis_ppSimDetectorEffect));
      } else {
        NormaliseRawHistToNEvents(H1D_jetPt_mcp, GetNEventsSelected_JetFramework(file_O2Analysis_ppSimDetectorEffect_unfoldingControl));
      }
    }
  }
  if (options.find("entriesNorm") != std::string::npos) {
    NormaliseYieldToNEntries(H1D_jetPt_mcp);
  }
}
void Get_Pt_spectrum_mcp(TH1D* &H1D_jetPt_mcp, int iDataset, int iRadius, bool fcontrolMC, std::string options) {
  Get_Pt_spectrum_mcp_preWidthScaling(H1D_jetPt_mcp, iDataset, iRadius, fcontrolMC, options);

  if (doWidthScalingAtEnd) {
    TransformRawHistToYield(H1D_jetPt_mcp);
  }
}

void Get_Pt_spectrum_mcp_fineBinning_preWidthScaling(TH1D* &H1D_jetPt_mcp, int iDataset, int iRadius, bool fcontrolMC, std::string options) {
  Get_Pt_spectrum_mcp_fineBinning_preWidthScalingAndEvtNorm(H1D_jetPt_mcp, iDataset, iRadius, fcontrolMC, options);

  if (options.find("evtNorm") != std::string::npos && doEvtNorm) {
    if (isPbPb == true || ppMcIsWeighted) {
      NormaliseRawHistToNEvents(H1D_jetPt_mcp, GetNEventsSelected_JetFramework_weighted(file_O2Analysis_ppSimDetectorEffect));
    } else {
      if (!fcontrolMC) {
        NormaliseRawHistToNEvents(H1D_jetPt_mcp, GetNEventsSelected_JetFramework(file_O2Analysis_ppSimDetectorEffect));
      } else {
        NormaliseRawHistToNEvents(H1D_jetPt_mcp, GetNEventsSelected_JetFramework(file_O2Analysis_ppSimDetectorEffect_unfoldingControl));
      }
    }
  }
  if (options.find("entriesNorm") != std::string::npos) {
    NormaliseYieldToNEntries(H1D_jetPt_mcp);
  }
}
void Get_Pt_spectrum_mcp_fineBinning(TH1D* &H1D_jetPt_mcp, int iDataset, int iRadius, bool fcontrolMC, std::string options) {
  Get_Pt_spectrum_mcp_fineBinning_preWidthScaling(H1D_jetPt_mcp, iDataset, iRadius, fcontrolMC, options);

  if (doWidthScalingAtEnd) {
    TransformRawHistToYield(H1D_jetPt_mcp);
  }
}

void Get_Pt_spectrum_mcd_fineBinning_preWidthScaling(TH1D* &H1D_jetPt_mcd, int iDataset, int iRadius, std::string options) {
  Get_Pt_spectrum_mcd_fineBinning_preWidthScalingAndEvtNorm(H1D_jetPt_mcd, iDataset, iRadius, options);

  if (options.find("evtNorm") != std::string::npos && doEvtNorm) {
    if (isPbPb == true || ppMcIsWeighted) {
      NormaliseRawHistToNEvents(H1D_jetPt_mcd, GetNEventsSelected_JetFramework_weighted(file_O2Analysis_ppSimDetectorEffect));
    } else {
      NormaliseRawHistToNEvents(H1D_jetPt_mcd, GetNEventsSelected_JetFramework(file_O2Analysis_ppSimDetectorEffect));
    }
  }
  if (options.find("entriesNorm") != std::string::npos) {
    NormaliseYieldToNEntries(H1D_jetPt_mcd);
  }
}
void Get_Pt_spectrum_mcd_fineBinning(TH1D* &H1D_jetPt_mcd, int iDataset, int iRadius, std::string options) {
  Get_Pt_spectrum_mcd_fineBinning_preWidthScaling(H1D_jetPt_mcd, iDataset, iRadius, options);

  if (doWidthScalingAtEnd) {
    TransformRawHistToYield(H1D_jetPt_mcd);
  }
}

void Get_Pt_spectrum_mcd_recBinning_preWidthScaling(TH1D* &H1D_jetPt_mcd, int iDataset, int iRadius, std::string options) {
  Get_Pt_spectrum_mcd_recBinning_preWidthScalingAndEvtNorm(H1D_jetPt_mcd, iDataset, iRadius, options);

  if (options.find("evtNorm") != std::string::npos && doEvtNorm) {
    if (isPbPb == true || ppMcIsWeighted) {
      NormaliseRawHistToNEvents(H1D_jetPt_mcd, GetNEventsSelected_JetFramework_weighted(file_O2Analysis_ppSimDetectorEffect));
    } else {
      NormaliseRawHistToNEvents(H1D_jetPt_mcd, GetNEventsSelected_JetFramework(file_O2Analysis_ppSimDetectorEffect));
    }
  }
  if (options.find("entriesNorm") != std::string::npos) {
    NormaliseYieldToNEntries(H1D_jetPt_mcd);
  }
}
void Get_Pt_spectrum_mcd_recBinning(TH1D* &H1D_jetPt_mcd, int iDataset, int iRadius, std::string options) {
  Get_Pt_spectrum_mcd_recBinning_preWidthScaling(H1D_jetPt_mcd, iDataset, iRadius, options);

  if (doWidthScalingAtEnd) {
    TransformRawHistToYield(H1D_jetPt_mcd);
  }
}

void Get_Pt_spectrum_mcd_preWidthScaling(TH1D* &H1D_jetPt_mcd, int iDataset, int iRadius, std::string options) {
  Get_Pt_spectrum_mcd_preWidthScalingAndEvtNorm(H1D_jetPt_mcd, iDataset, iRadius, options);

  if (options.find("evtNorm") != std::string::npos && doEvtNorm) {
    if (isPbPb == true || ppMcIsWeighted) {
      NormaliseRawHistToNEvents(H1D_jetPt_mcd, GetNEventsSelected_JetFramework_weighted(file_O2Analysis_ppSimDetectorEffect));
    } else {
      NormaliseRawHistToNEvents(H1D_jetPt_mcd, GetNEventsSelected_JetFramework(file_O2Analysis_ppSimDetectorEffect));
    }
  }
  if (options.find("entriesNorm") != std::string::npos) {
    NormaliseYieldToNEntries(H1D_jetPt_mcd);
  }
}
void Get_Pt_spectrum_mcd(TH1D* &H1D_jetPt_mcd, int iDataset, int iRadius, std::string options) {
  Get_Pt_spectrum_mcd_preWidthScaling(H1D_jetPt_mcd, iDataset, iRadius, options);

  if (doWidthScalingAtEnd) {
    TransformRawHistToYield(H1D_jetPt_mcd);
  }
}

void Get_Pt_spectrum_mcdMatched_preWidthScaling(TH1D* &H1D_jetPt_mcdMatched, int iDataset, int iRadius, std::string options) {
  Get_Pt_spectrum_mcdMatched_preWidthScalingAndEvtNorm(H1D_jetPt_mcdMatched, iDataset, iRadius, options);

  if (options.find("evtNorm") != std::string::npos && doEvtNorm) {
    if (isPbPb == true || ppMcIsWeighted) {
      NormaliseRawHistToNEvents(H1D_jetPt_mcdMatched, GetNEventsSelected_JetFramework_weighted(file_O2Analysis_ppSimDetectorEffect));
    } else {
      NormaliseRawHistToNEvents(H1D_jetPt_mcdMatched, GetNEventsSelected_JetFramework(file_O2Analysis_ppSimDetectorEffect));
    }
  }
  if (options.find("entriesNorm") != std::string::npos) {
    NormaliseYieldToNEntries(H1D_jetPt_mcdMatched);
  }
}
void Get_Pt_spectrum_mcdMatched(TH1D* &H1D_jetPt_mcdMatched, int iDataset, int iRadius, std::string options) {
  Get_Pt_spectrum_mcdMatched_preWidthScaling(H1D_jetPt_mcdMatched, iDataset, iRadius, options);

  if (doWidthScalingAtEnd) {
    TransformRawHistToYield(H1D_jetPt_mcdMatched);
  }
}

void Get_Pt_spectrum_mcdMatched_recBinning_preWidthScaling(TH1D* &H1D_jetPt_mcdMatched, int iDataset, int iRadius, std::string options) {
  Get_Pt_spectrum_mcdMatched_recBinning_preWidthScalingAndEvtNorm(H1D_jetPt_mcdMatched, iDataset, iRadius, options);

  if (options.find("evtNorm") != std::string::npos && doEvtNorm) {
    if (isPbPb == true || ppMcIsWeighted) {
      NormaliseRawHistToNEvents(H1D_jetPt_mcdMatched, GetNEventsSelected_JetFramework_weighted(file_O2Analysis_ppSimDetectorEffect));
    } else {
      NormaliseRawHistToNEvents(H1D_jetPt_mcdMatched, GetNEventsSelected_JetFramework(file_O2Analysis_ppSimDetectorEffect));
    }
  }
  if (options.find("entriesNorm") != std::string::npos) {
    NormaliseYieldToNEntries(H1D_jetPt_mcdMatched);
  }
}
void Get_Pt_spectrum_mcdMatched_recBinning(TH1D* &H1D_jetPt_mcdMatched, int iDataset, int iRadius, std::string options) {
  Get_Pt_spectrum_mcdMatched_recBinning_preWidthScaling(H1D_jetPt_mcdMatched, iDataset, iRadius, options);

  if (doWidthScalingAtEnd) {
    TransformRawHistToYield(H1D_jetPt_mcdMatched);
  }
}

void Get_Pt_spectrum_mcdMatched_fineBinning_preWidthScaling(TH1D* &H1D_jetPt_mcdMatched, int iDataset, int iRadius, std::string options) {
  Get_Pt_spectrum_mcdMatched_fineBinning_preWidthScalingAndEvtNorm(H1D_jetPt_mcdMatched, iDataset, iRadius, options);

  if (options.find("evtNorm") != std::string::npos && doEvtNorm) {
    if (isPbPb == true || ppMcIsWeighted) {
      NormaliseRawHistToNEvents(H1D_jetPt_mcdMatched, GetNEventsSelected_JetFramework_weighted(file_O2Analysis_ppSimDetectorEffect));
    } else {
      NormaliseRawHistToNEvents(H1D_jetPt_mcdMatched, GetNEventsSelected_JetFramework(file_O2Analysis_ppSimDetectorEffect));
    }
  }
  if (options.find("entriesNorm") != std::string::npos) {
    NormaliseYieldToNEntries(H1D_jetPt_mcdMatched);
  }
}
void Get_Pt_spectrum_mcdMatched_fineBinning(TH1D* &H1D_jetPt_mcdMatched, int iDataset, int iRadius, std::string options) {
  Get_Pt_spectrum_mcdMatched_fineBinning_preWidthScaling(H1D_jetPt_mcdMatched, iDataset, iRadius, options);

  if (doWidthScalingAtEnd) {
    TransformRawHistToYield(H1D_jetPt_mcdMatched);
  }
}

void Get_Pt_spectrum_mcpMatched_preWidthScaling(TH1D* &H1D_jetPt_mcpMatched, int iDataset, int iRadius, std::string options) {
  Get_Pt_spectrum_mcpMatched_preWidthScalingAndEvtNorm(H1D_jetPt_mcpMatched, iDataset, iRadius, options);

  if (options.find("evtNorm") != std::string::npos && doEvtNorm) {
    if (isPbPb == true || ppMcIsWeighted) {
      NormaliseRawHistToNEvents(H1D_jetPt_mcpMatched, GetNEventsSelected_JetFramework_weighted(file_O2Analysis_ppSimDetectorEffect));
    } else {
      NormaliseRawHistToNEvents(H1D_jetPt_mcpMatched, GetNEventsSelected_JetFramework(file_O2Analysis_ppSimDetectorEffect));
    }
  }
  if (options.find("entriesNorm") != std::string::npos) {
    NormaliseYieldToNEntries(H1D_jetPt_mcpMatched);
  }
}
void Get_Pt_spectrum_mcpMatched(TH1D* &H1D_jetPt_mcpMatched, int iDataset, int iRadius, std::string options) {
  Get_Pt_spectrum_mcpMatched_preWidthScaling(H1D_jetPt_mcpMatched, iDataset, iRadius, options);

  if (doWidthScalingAtEnd) {
    TransformRawHistToYield(H1D_jetPt_mcpMatched);
  }
}

void Get_Pt_spectrum_mcpMatched_fineBinning_preWidthScaling(TH1D* &H1D_jetPt_mcpMatched, int iDataset, int iRadius, std::string options) {
  Get_Pt_spectrum_mcpMatched_fineBinning_preWidthScalingAndEvtNorm(H1D_jetPt_mcpMatched, iDataset, iRadius, options);



  if (options.find("evtNorm") != std::string::npos && doEvtNorm) {
    if (isPbPb == true || ppMcIsWeighted) {
      NormaliseRawHistToNEvents(H1D_jetPt_mcpMatched, GetNEventsSelected_JetFramework_weighted(file_O2Analysis_ppSimDetectorEffect));
    } else {
      NormaliseRawHistToNEvents(H1D_jetPt_mcpMatched, GetNEventsSelected_JetFramework(file_O2Analysis_ppSimDetectorEffect));
    }
  }
  if (options.find("entriesNorm") != std::string::npos) {
    NormaliseYieldToNEntries(H1D_jetPt_mcpMatched);
  }
}
void Get_Pt_spectrum_mcpMatched_fineBinning(TH1D* &H1D_jetPt_mcpMatched, int iDataset, int iRadius, std::string options) {
  Get_Pt_spectrum_mcpMatched_fineBinning_preWidthScaling(H1D_jetPt_mcpMatched, iDataset, iRadius, options);

  if (doWidthScalingAtEnd) {
    TransformRawHistToYield(H1D_jetPt_mcpMatched);
  }
}

void Get_Pt_spectrum_mcp_recBinning_preWidthScaling(TH1D* &H1D_jetPt_mcp, int iDataset, int iRadius, bool controlMC, std::string options) {
  Get_Pt_spectrum_mcp_recBinning_preWidthScalingAndEvtNorm(H1D_jetPt_mcp, iDataset, iRadius, controlMC, options);

  if (options.find("evtNorm") != std::string::npos && doEvtNorm) {
    if (isPbPb == true || ppMcIsWeighted) {
      NormaliseRawHistToNEvents(H1D_jetPt_mcp, GetNEventsSelected_JetFramework_weighted(file_O2Analysis_ppSimDetectorEffect));
    } else {
      if (!controlMC) {
        NormaliseRawHistToNEvents(H1D_jetPt_mcp, GetNEventsSelected_JetFramework(file_O2Analysis_ppSimDetectorEffect));
      } else {
        NormaliseRawHistToNEvents(H1D_jetPt_mcp, GetNEventsSelected_JetFramework(file_O2Analysis_ppSimDetectorEffect_unfoldingControl));
      }
    }
  }
  if (options.find("entriesNorm") != std::string::npos) {
    NormaliseYieldToNEntries(H1D_jetPt_mcp);
  }
}
void Get_Pt_spectrum_mcp_recBinning(TH1D* &H1D_jetPt_mcp, int iDataset, int iRadius, bool fcontrolMC, std::string options) {
  Get_Pt_spectrum_mcp_recBinning_preWidthScaling(H1D_jetPt_mcp, iDataset, iRadius, fcontrolMC, options);

  if (doWidthScalingAtEnd) {
    TransformRawHistToYield(H1D_jetPt_mcp);
  }
}

void Get_Pt_spectrum_bkgCorrected_genBinning_preWidthScaling(TH1D* &H1D_jetPt, int iDataset, int iRadius, float* centRange, std::string options) {
  Get_Pt_spectrum_bkgCorrected_genBinning_preWidthScalingAndEvtNorm(H1D_jetPt, iDataset, iRadius, centRange, options);


  if (options.find("evtNorm") != std::string::npos && doEvtNorm) {
    if (isPbPb) {
      NormaliseRawHistToNEvents(H1D_jetPt, GetNEventsSelectedCentrality_JetFramework(file_O2Analysis_list[iDataset], centRange[0], centRange[1], trainId));
    } else {
      if (!controlMC) {
        NormaliseRawHistToNEvents(H1D_jetPt, GetNEventsSelectedCentrality_JetFramework(file_O2Analysis_list[iDataset], -999, 999, trainId));
      } else {
        NormaliseRawHistToNEvents(H1D_jetPt, GetNEventsSelectedCentrality_JetFramework(file_O2Analysis_ppSimDetectorEffect_unfoldingControl, -999, 999, trainId));
      }
    }
  }
  if (options.find("entriesNorm") != std::string::npos) {
    NormaliseYieldToNEntries(H1D_jetPt);
  }
}
void Get_Pt_spectrum_bkgCorrected_genBinning(TH1D* &H1D_jetPt, int iDataset, int iRadius, float* centRange, std::string options) {
  Get_Pt_spectrum_bkgCorrected_genBinning_preWidthScaling(H1D_jetPt, iDataset, iRadius, centRange, options);

  if (doWidthScalingAtEnd) {
    TransformRawHistToYield(H1D_jetPt);
  }
}

void Get_Pt_spectrum_bkgCorrected_fineBinning_preWidthScaling(TH1D* &H1D_jetPt, int iDataset, int iRadius, float* centRange, std::string options) {
  Get_Pt_spectrum_bkgCorrected_fineBinning_preWidthScalingAndEvtNorm(H1D_jetPt, iDataset, iRadius, centRange, options);

  if (options.find("evtNorm") != std::string::npos && doEvtNorm) {
    if (isPbPb) {
      NormaliseRawHistToNEvents(H1D_jetPt, GetNEventsSelectedCentrality_JetFramework(file_O2Analysis_list[iDataset], centRange[0], centRange[1], trainId));
    } else {
      if (!controlMC) {
        NormaliseRawHistToNEvents(H1D_jetPt, GetNEventsSelectedCentrality_JetFramework(file_O2Analysis_list[iDataset], -999, 999, trainId));
      } else {
        NormaliseRawHistToNEvents(H1D_jetPt, GetNEventsSelectedCentrality_JetFramework(file_O2Analysis_ppSimDetectorEffect_unfoldingControl, -999, 999, trainId));
      }
    }
  }
  if (options.find("entriesNorm") != std::string::npos) {
    NormaliseYieldToNEntries(H1D_jetPt);
  }
}
void Get_Pt_spectrum_bkgCorrected_fineBinning(TH1D* &H1D_jetPt, int iDataset, int iRadius, float* centRange, std::string options) {
  Get_Pt_spectrum_bkgCorrected_fineBinning_preWidthScaling(H1D_jetPt, iDataset, iRadius, centRange, options);

  if (doWidthScalingAtEnd) {
    TransformRawHistToYield(H1D_jetPt);
  }
}


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////// Response matrix functions ///////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void Get_PtResponseMatrix_DetectorAndFluctuationsCombined_fineBinning(TH2D* &H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined_fineBinning, TH2D* H2D_jetPtResponseMatrix_detectorResponse, TH2D* H2D_jetPtResponseMatrix_fluctuations, int iDataset, int iRadius, float* centRange, __attribute__ ((unused)) std::string options) {
  // https://github.com/alisw/AliPhysics/blob/master/PWGJE/PWGJE/AliAnaChargedJetResponseMaker.cxx for ann example that works, by marta verveij

  TString partialUniqueSpecifier = Datasets[iDataset]+"_R="+Form("%.1f",arrayRadius[iRadius])+"_@cent["+Form("%.1f", centRange[0])+","+Form("%.1f", centRange[1])+"]";

  // matrix product of fluct response times det response; assumes the two are of the same size binning wise
  // Careful: xy of hist and ij of Resp(i,j) are inverted ! hist(j,i) = matrix(i,j) and so if matrix(i,j)=SUM(matrixA(i,k)matrixB(k,j)) then hist(j,i)=SUM(histA(k,i)histB(j,k)), and if we replace j,i by gen,rec we get hist(gen,rec)=SUM(histA(k,rec)histB(gen,k))
  H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined_fineBinning = (TH2D*)GetMatrixProductTH2xTH2(H2D_jetPtResponseMatrix_fluctuations, H2D_jetPtResponseMatrix_detectorResponse).Clone("Get_PtResponseMatrix_DetectorAndFluctuationsCombined_fineBinning"+partialUniqueSpecifier);


  if (drawIntermediateResponseMatrices) {
    TString* pdfName_preRebin = new TString("responseMatrix_combined_preRebinning"+jetType[iJetType]+"_"+Datasets[iDataset]);
    TString* pdfName_preRebin_logz = new TString("responseMatrix_combined_preRebinning"+(TString)"_R="+Form("%.1f",arrayRadius[iRadius])+"_"+Datasets[iDataset]+"_logz");

    TString textContext_preRebin(contextCustomOneField(*texDatasetsComparisonCommonDenominator, ""));
    
    Draw_TH2_Histogram(H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined_fineBinning, textContext_preRebin, pdfName_preRebin, texPtJetRecX, texPtJetGenX, texCollisionDataInfo, drawnWindowAuto, legendPlacementAuto, contextPlacementAuto, "");
    Draw_TH2_Histogram(H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined_fineBinning, textContext_preRebin, pdfName_preRebin_logz, texPtJetRecX, texPtJetGenX, texCollisionDataInfo, drawnWindowAuto, legendPlacementAuto, contextPlacementAuto, "logz");
  }
  // cout << "bin(topleft 1,N) = " << H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined_fineBinning->GetBinContent(1,H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined_fineBinning->GetNbinsY()) << endl;
  // cout << "bin(bottom left 1,1) = " << H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined_fineBinning->GetBinContent(1,1) << endl;
}


void Get_PtResponseMatrix_DetectorAndFluctuationsCombined(TH2D* &H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined, TH2D* H2D_jetPtResponseMatrix_detectorResponse, TH2D* H2D_jetPtResponseMatrix_fluctuations, int iDataset, int iRadius, float* centRange, std::string options) {
  // https://github.com/alisw/AliPhysics/blob/master/PWGJE/PWGJE/AliAnaChargedJetResponseMaker.cxx for ann example that works, by marta verveij

  TString partialUniqueSpecifier = Datasets[iDataset]+"_R="+Form("%.1f",arrayRadius[iRadius])+"_@cent["+Form("%.1f", centRange[0])+","+Form("%.1f", centRange[1])+"]";

  TH2D* H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined_preRebin;
  Get_PtResponseMatrix_DetectorAndFluctuationsCombined_fineBinning(H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined_preRebin, H2D_jetPtResponseMatrix_detectorResponse, H2D_jetPtResponseMatrix_fluctuations, iDataset, iRadius, centRange, options);

  // Merge bins of the combined response matrix with fine binning to get the coarse one
  TH1D* priorSpectrumMerging;
  bool debugBool = false;
  Get_Pt_spectrum_mcp_fineBinning(priorSpectrumMerging, iDataset, iRadius, false, options); //take mcp as prior by default
  if (options.find("mcpPriorMerging") != std::string::npos) {
    priorSpectrumMerging->Reset("M");
    Get_Pt_spectrum_mcp_fineBinning(priorSpectrumMerging, iDataset, iRadius, false, options); 
    H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined = (TH2D*)RebinVariableBins2D_PriorWeightedBinMerging(H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined_preRebin, nBinPtJetsRec[iRadius], nBinPtJetsGen[iRadius], ptBinsJetsRec[iRadius], ptBinsJetsGen[iRadius], priorSpectrumMerging, debugBool).Clone("Get_PtResponseMatrix_DetectorAndFluctuationsCombined"+partialUniqueSpecifier);
  }
  if (options.find("mcdPriorMerging") != std::string::npos) {
    priorSpectrumMerging->Reset("M");
    Get_Pt_spectrum_mcd_fineBinning(priorSpectrumMerging, iDataset, iRadius, options);
    H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined = (TH2D*)RebinVariableBins2D_PriorWeightedBinMerging(H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined_preRebin, nBinPtJetsRec[iRadius], nBinPtJetsGen[iRadius], ptBinsJetsRec[iRadius], ptBinsJetsGen[iRadius], priorSpectrumMerging, debugBool).Clone("Get_PtResponseMatrix_DetectorAndFluctuationsCombined"+partialUniqueSpecifier);
  }
  if (options.find("measuredPriorMerging") != std::string::npos) {
    priorSpectrumMerging->Reset("M");
    Get_Pt_spectrum_bkgCorrected_fineBinning(priorSpectrumMerging, iDataset, iRadius, centRange, options);
    H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined = (TH2D*)RebinVariableBins2D_PriorWeightedBinMerging(H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined_preRebin, nBinPtJetsRec[iRadius], nBinPtJetsGen[iRadius], ptBinsJetsRec[iRadius], ptBinsJetsGen[iRadius], priorSpectrumMerging, debugBool).Clone("Get_PtResponseMatrix_DetectorAndFluctuationsCombined"+partialUniqueSpecifier);
  }
  if (options.find("noMergingPrior") != std::string::npos) {
    H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined = (TH2D*)RebinVariableBins2D(H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined_preRebin, nBinPtJetsRec[iRadius], nBinPtJetsGen[iRadius], ptBinsJetsRec[iRadius], ptBinsJetsGen[iRadius], debugBool).Clone("Get_PtResponseMatrix_DetectorAndFluctuationsCombined"+partialUniqueSpecifier);
  }
  if (options.find("testAliPhysics") != std::string::npos) {
    H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined = (TH2D*)RebinVariableBins2D_aliPhysics(H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined_preRebin, nBinPtJetsRec[iRadius], nBinPtJetsGen[iRadius], ptBinsJetsRec[iRadius], ptBinsJetsGen[iRadius], debugBool)->Clone("Get_PtResponseMatrix_DetectorAndFluctuationsCombined"+partialUniqueSpecifier);
  }
  // normalising priorSpectrum with evtNorm doesn't change anything as the weighting does prior_content(i)/prior_integral()
  // dividing priorSpectrum by binwidth doesn't change anything for the same reason



  // H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined = (TH2D*)RebinVariableBins2D(H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined_preRebin, nBinPtJetsRec[iRadius], nBinPtJetsGen[iRadius], ptBinsJetsRec[iRadius], ptBinsJetsGen[iRadius], true).Clone("Get_PtResponseMatrix_DetectorAndFluctuationsCombined"+partialUniqueSpecifier);


  // When looking at combined response matrix before normalisation, large bins in y will look strange, and out of place compared to other bin slices of same size; this is because it potentially merges A LOT of bins together; it'll look a lot better after normalisation:


  // debug
  // for(int iBinY = 0; iBinY <= H2D_jetPtResponseMatrix_detectorResponse->GetNbinsY()+1; iBinY++){ // 0 and n+1 take underflow and overflow into account
  //   for(int iBinX = 0; iBinX <= H2D_jetPtResponseMatrix_detectorResponse->GetNbinsX()+1; iBinX++){ // 0 and n+1 take underflow and overflow into account
  //     if (iBinX == 0) {
  //       cout << "iBinX = 0 --> H2D_jetPtResponseMatrix_detectorResponse->GetBinContent(0,"<< iBinY << ") = " << H2D_jetPtResponseMatrix_detectorResponse->GetBinContent(0,iBinY) << endl;
  //       cout << "iBinX = 0 --> H2D_jetPtResponseMatrix_fluctuations->GetBinContent(0,"<< iBinY << ") = " << H2D_jetPtResponseMatrix_fluctuations->GetBinContent(0,iBinY) << endl;
  //     }
  //   }
  // }
  if (drawIntermediateResponseMatrices) {
    TH2D* H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined_preNorm = (TH2D*)H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined->Clone("H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined_preNorm"+partialUniqueSpecifier);

    TString* pdfName = new TString("responseMatrix_combined_preNorm"+jetType[iJetType]+"_"+Datasets[iDataset]);
    TString* pdfName_logz = new TString("responseMatrix_combined_preNorm"+(TString)"_R="+Form("%.1f",arrayRadius[iRadius])+"_"+Datasets[iDataset]+"_logz");

    TString textContext(contextCustomOneField(*texDatasetsComparisonCommonDenominator, ""));

    Draw_TH2_Histogram(H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined_preNorm, textContext, pdfName, texPtJetRecX, texPtJetGenX, texCollisionDataInfo, drawnWindowAuto, legendPlacementAuto, contextPlacementAuto, "");
    Draw_TH2_Histogram(H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined_preNorm, textContext, pdfName_logz, texPtJetRecX, texPtJetGenX, texCollisionDataInfo, drawnWindowAuto, legendPlacementAuto, contextPlacementAuto, "logz");
  }

  // TransformRawResponseToYieldResponse(H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined);
  // cout << "Should I normalise the combined matrix?" << endl;
  // cout << "     Marta doesn't do it, but it looks like I need it due to merging (some bins are very large)" << endl;
  // cout << "     it's actually done in AliAnaChargedJetResponseMaker::MakeResponseMatrixRebin" << endl;


  if (useYSliceNorm) {
    NormaliseYSlicesToOne(H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined); // Marta doesn't do it, but it looks like I need it due to merging (some bins are very large) ; actually marta probably uses it: AliAnaChargedJetResponseMaker::MakeResponseMatrixRebin does it by default inside the rebinning function, and it takes into account the whole Fine range (1fine, Nfine)
  } 
  if (scaleRespByWidth) {
    H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined->Scale(1., "width");
  }

  if (drawIntermediateResponseMatrices) {
    TH2D* H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined_preWeighting = (TH2D*)H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined->Clone("H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined_preWeighting"+partialUniqueSpecifier);

    TString* pdfName = new TString("responseMatrix_combined_preWeighting"+jetType[iJetType]+"_"+Datasets[iDataset]);
    TString* pdfName_logz = new TString("responseMatrix_combined_preWeighting"+(TString)"_R="+Form("%.1f",arrayRadius[iRadius])+"_"+Datasets[iDataset]+"_logz");

    TString textContext(contextCustomOneField(*texDatasetsComparisonCommonDenominator, ""));

    Draw_TH2_Histogram(H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined_preWeighting, textContext, pdfName, texPtJetRecX, texPtJetGenX, texCollisionDataInfo, drawnWindowAuto, legendPlacementAuto, contextPlacementAuto, "");
    Draw_TH2_Histogram(H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined_preWeighting, textContext, pdfName_logz, texPtJetRecX, texPtJetGenX, texCollisionDataInfo, drawnWindowAuto, legendPlacementAuto, contextPlacementAuto, "logz");
  }
}

void ReweightResponseMatrixWithPrior(TH2D* &H2D_jetPtResponseMatrix, int iDataset, int iRadius, float* centRange, std::string options) {
  
  // before this, all y-slices (ie pt gen slices) have been normalised to 1;means each pt gen slice has a proba of 1
  // withthis function, we give each slice a weight so that they have different normalisation values, corresponding to the prior 

  // prior choice; none by default (flat)
  TH1D* priorSpectrumWeighting;
  if (options.find("mcpPriorUnfolding") != std::string::npos) {
    if (!normGenAndMeasByNEvts) {
      Get_Pt_spectrum_mcp_preWidthScalingAndEvtNorm(priorSpectrumWeighting, iDataset, iRadius, false, options); 
    } else {
      Get_Pt_spectrum_mcp_preWidthScaling(priorSpectrumWeighting, iDataset, iRadius, false, options); 
    }
    WeightMatrixWithPrior(H2D_jetPtResponseMatrix, priorSpectrumWeighting);
    // for (int i = 1; i < priorSpectrumWeighting->GetNbinsX(); i++)
    // {
    //   cout << "prior(" << i << ")" << priorSpectrumWeighting->GetBinContent(i)<< endl;
    //   cout << "responseIntegralLine(" << i << ")" << H2D_jetPtResponseMatrix->Integral(1,H2D_jetPtResponseMatrix->GetNbinsX(), i, i)<< endl;
    // }
    
  }
  if (options.find("mcdPriorUnfolding") != std::string::npos) {
    if (!normGenAndMeasByNEvts) {
      Get_Pt_spectrum_mcd_preWidthScalingAndEvtNorm(priorSpectrumWeighting, iDataset, iRadius, options); 
    } else {
      Get_Pt_spectrum_mcd_preWidthScaling(priorSpectrumWeighting, iDataset, iRadius, options); 
    }
    WeightMatrixWithPrior(H2D_jetPtResponseMatrix, priorSpectrumWeighting);
  }
  if (options.find("measuredPriorUnfolding") != std::string::npos) {
    if (!normGenAndMeasByNEvts) {
      Get_Pt_spectrum_bkgCorrected_genBinning_preWidthScalingAndEvtNorm(priorSpectrumWeighting, iDataset, iRadius, centRange, options);
    } else {
      Get_Pt_spectrum_bkgCorrected_genBinning_preWidthScaling(priorSpectrumWeighting, iDataset, iRadius, centRange, options); 
    }
    WeightMatrixWithPrior(H2D_jetPtResponseMatrix, priorSpectrumWeighting);
  }
  if (options.find("testAliPhysics") != std::string::npos) {
    if (!normGenAndMeasByNEvts) {
      Get_Pt_spectrum_bkgCorrected_genBinning_preWidthScalingAndEvtNorm(priorSpectrumWeighting, iDataset, iRadius, centRange, options); 
    } else {
      Get_Pt_spectrum_bkgCorrected_genBinning_preWidthScaling(priorSpectrumWeighting, iDataset, iRadius, centRange, options); 
    }
    H2D_jetPtResponseMatrix = (TH2D*)NormalizeResponsMatrixYaxisWithPrior(H2D_jetPtResponseMatrix, priorSpectrumWeighting)->Clone(H2D_jetPtResponseMatrix->GetName()+(TString)"_testAliPhysics");
  }
  // cout << "((((((((((((((((((((((((()))))))))))))))))))))))))" << endl;
  // cout << "pre norm that shouldn't be" << endl;
  // cout << "H2D_jetPtResponseMatrix->Integral(1, N, 1, 1)" << H2D_jetPtResponseMatrix->Integral(1, H2D_jetPtResponseMatrix->GetNbinsX(), 1, 1) << endl;
  // cout << "H2D_jetPtResponseMatrix->Integral(0, -1, 1, 1)" << H2D_jetPtResponseMatrix->Integral(0, -1, 1, 1) << endl;
  // // NormaliseXSlicesToOneNoUnderOverFlows(H2D_jetPtResponseMatrix);
  // // NormaliseYSlicesToOneNoUnderOverFlows(H2D_jetPtResponseMatrix);
  // cout << "post norm that shouldn't be" << endl;
  // cout << "H2D_jetPtResponseMatrix->Integral(1, N, 1, 1)" << H2D_jetPtResponseMatrix->Integral(1, H2D_jetPtResponseMatrix->GetNbinsX(), 1, 1) << endl;
  // cout << "H2D_jetPtResponseMatrix->Integral(0, -1, 1, 1)" << H2D_jetPtResponseMatrix->Integral(0, -1, 1, 1) << endl;
  // cout << "((((((((((((((((((((((((()))))))))))))))))))))))))" << endl;
}


void ReweightResponseMatrixWithPrior_fineBinning(TH2D* &H2D_jetPtResponseMatrix, int iDataset, int iRadius, float* centRange, std::string options) {
  
  // before this, all y-slices (ie pt gen slices) have been normalised to 1;means each pt gen slice has a proba of 1
  // withthis function, we give each slice a weight so that they have different normalisation values, corresponding to the prior 

  // prior choice; none by default (flat)
  TH1D* priorSpectrumWeighting;
  if (options.find("mcpPriorUnfolding") != std::string::npos) {
    if (!normGenAndMeasByNEvts) {
      Get_Pt_spectrum_mcp_fineBinning_preWidthScalingAndEvtNorm(priorSpectrumWeighting, iDataset, iRadius, false, options); 
    } else {
      Get_Pt_spectrum_mcp_fineBinning_preWidthScaling(priorSpectrumWeighting, iDataset, iRadius, false, options); 
    }
    WeightMatrixWithPrior(H2D_jetPtResponseMatrix, priorSpectrumWeighting);
  }
  if (options.find("mcdPriorUnfolding") != std::string::npos) {
    if (!normGenAndMeasByNEvts) {
      Get_Pt_spectrum_mcd_fineBinning_preWidthScalingAndEvtNorm(priorSpectrumWeighting, iDataset, iRadius, options); 
    } else {
      Get_Pt_spectrum_mcd_fineBinning_preWidthScaling(priorSpectrumWeighting, iDataset, iRadius, options); 
    }
    WeightMatrixWithPrior(H2D_jetPtResponseMatrix, priorSpectrumWeighting);
  }
  if (options.find("measuredPriorUnfolding") != std::string::npos) {
    if (!normGenAndMeasByNEvts) {
      Get_Pt_spectrum_bkgCorrected_fineBinning_preWidthScalingAndEvtNorm(priorSpectrumWeighting, iDataset, iRadius, centRange, options); 
    } else {
      Get_Pt_spectrum_bkgCorrected_fineBinning_preWidthScaling(priorSpectrumWeighting, iDataset, iRadius, centRange, options); 
    }
    WeightMatrixWithPrior(H2D_jetPtResponseMatrix, priorSpectrumWeighting);
  }
  if (options.find("testAliPhysics") != std::string::npos) {
    if (!normGenAndMeasByNEvts) {
      Get_Pt_spectrum_bkgCorrected_fineBinning_preWidthScalingAndEvtNorm(priorSpectrumWeighting, iDataset, iRadius, centRange, options); 
    } else {
      Get_Pt_spectrum_bkgCorrected_fineBinning_preWidthScaling(priorSpectrumWeighting, iDataset, iRadius, centRange, options); 
    }
    H2D_jetPtResponseMatrix = (TH2D*)NormalizeResponsMatrixYaxisWithPrior(H2D_jetPtResponseMatrix, priorSpectrumWeighting)->Clone(H2D_jetPtResponseMatrix->GetName()+(TString)"_testAliPhysics");
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
  // project H3D onto a H2D, option "yz" means y goes on y-axis while z goes on x-axis, and so H2D_gen_det_geoMatched will be (x=mcd, y=mcp)
  H2D_gen_det_geoMatched = (TH2D*)H3D_jetRpartPtdetPt->Project3D(partialUniqueSpecifier+"_genrec_e_yz"); //can't use letter D in this or it seems to replace the histogram in current pad (see documentation of ProjectionX function. Isn't mentioned in project3D sadly)

  // keep (gen, gen) for the bins; rec will be introduced in the fluctuation response, and by multiplication will stay in the combined matrix
  TH2D* H2D_response = (TH2D*)RebinVariableBins2D(H2D_gen_det_geoMatched, nBinPtJetsFine[iRadius], nBinPtJetsFine[iRadius], ptBinsJetsFine[iRadius], ptBinsJetsFine[iRadius]).Clone("Get_PtResponseMatrix_detectorResponse_rebinned"+partialUniqueSpecifier);

  if (useYSliceNorm) {
    NormaliseYSlicesToOne(H2D_response);
  }
  if (normDetRespByNEvts) {
    if (isPbPb || ppMcIsWeighted) {
      H2D_response->Scale(1./GetNEventsSelected_JetFramework_weighted(file_O2Analysis_ppSimDetectorEffect));
    } else {
      int Nevents = GetNEventsSelected_JetFramework(file_O2Analysis_ppSimDetectorEffect);
      for (int iBinX = 0; iBinX < H2D_response->GetNbinsX(); iBinX++) {
        for (int iBinY = 0; iBinY < H2D_response->GetNbinsY(); iBinY++) {
          H2D_response->SetBinContent(iBinX, iBinY, H2D_response->GetBinContent(iBinX, iBinY) * 1./Nevents);
          H2D_response->SetBinError(iBinX, iBinY, H2D_response->GetBinError(iBinX, iBinY) * 1./Nevents);
        }
      }
      // H2D_response->Scale(1./GetNEventsSelected_JetFramework(file_O2Analysis_ppSimDetectorEffect));
    }
  }
  cout << "Detector response building: errors here should probably be reduced to take into account correlations, as the normalisation factor is built from same matrix" << endl;

  // H2D_jetPtResponseMatrix_detectorResponse = (TH2D*)H2D_gen_det_geoMatched_rebinned->Clone("H2D_jetPtResponseMatrix_detectorResponse"+partialUniqueSpecifier); // should be using the one below; this one is just a test
  H2D_jetPtResponseMatrix_detectorResponse = (TH2D*)H2D_response->Clone("H2D_jetPtResponseMatrix_detectorResponse"+partialUniqueSpecifier); 

  // cout << "............................................................." <<endl;
  // cout << "............................................................." <<endl;
  // for(int iBinX = 1; iBinX <= H2D_jetPtResponseMatrix_detectorResponse->GetNbinsX(); iBinX++){ // 0 and n+1 would take underflow and overflow into account, don't want that
  //   for(int iBinY = 1; iBinY <= H2D_jetPtResponseMatrix_detectorResponse->GetNbinsY(); iBinY++){ // 0 and n+1 would take underflow and overflow into account, don't want that
  //     cout << "iBinX = " << iBinX << ", iBinY = " << iBinY << "         --------          detResponseContent = " << H2D_jetPtResponseMatrix_detectorResponse->GetBinContent(iBinX, iBinY) << ", detResponseError = " << H2D_jetPtResponseMatrix_detectorResponse->GetBinError(iBinX, iBinY) << endl;
  //   }
  // }
  // cout << "............................................................." <<endl;
  // cout << "............................................................." <<endl;
}

void Get_PtResponseMatrix_Fluctuations(TH2D* &H2D_jetPtResponseMatrix_fluctuations, int iDataset, int iRadius, float* centRange) { 
  // see Hiroki Yokoyama thesis
  // iRadius is for chosing the pT binning
  TString partialUniqueSpecifier = Datasets[iDataset]+"_R="+Form("%.1f",arrayRadius[iRadius])+"_@cent["+Form("%.1f", centRange[0])+","+Form("%.1f", centRange[1])+"]";

  if (isPbPb == false){
    TH2D H2D_identity = TH2D("H2D_response_"+partialUniqueSpecifier, "H2D_response_"+partialUniqueSpecifier, nBinPtJetsFine[iRadius], ptBinsJetsFine[iRadius], nBinPtJetsFine[iRadius], ptBinsJetsFine[iRadius]);
    for(int iBinRec = 0; iBinRec <= H2D_identity.GetNbinsX()+1; iBinRec++){
      H2D_identity.SetBinContent(iBinRec, iBinRec, 1);
      H2D_identity.SetBinError(iBinRec, iBinRec, 0);
    }
    H2D_jetPtResponseMatrix_fluctuations = (TH2D*)H2D_identity.Clone("H2D_jetPtResponseMatrix_fluctuations"+partialUniqueSpecifier);
  } else {
    cout << "I should check that the average of each ptGen slice is as displaced to the diagonal as the randomCone distrib is; ie should I use GetBinLowEdge or GetBinLowEdge+width for ptGen" << endl;


    TH2D* H2D_fluctuations_centrality;
    TH1D* H1D_fluctuations;

    H2D_fluctuations_centrality = (TH2D*)((TH2D*)file_O2Analysis_list[iDataset]->Get(analysisWorkflowData+"/h2_centrality_rhorandomconewithoutleadingjet"))->Clone("Get_PtResponseMatrix_Fluctuations"+Datasets[iDataset]);
    H2D_fluctuations_centrality->Sumw2();

    int ibinCent_low, ibinCent_high;

    ibinCent_low = H2D_fluctuations_centrality->GetXaxis()->FindBin(centRange[0]);
    ibinCent_high = H2D_fluctuations_centrality->GetXaxis()->FindBin(centRange[1]);
    H1D_fluctuations = (TH1D*)H2D_fluctuations_centrality->ProjectionY("bkgFluctuationCentrality_highRes_"+partialUniqueSpecifier, ibinCent_low, ibinCent_high, "e");

    // NormaliseYieldToNEntries(H1D_fluctuations); // normalising fluctuations to 1
    NormaliseRawHistToIntegral(H1D_fluctuations); // normalising fluctuations to 1
      
    TH2D H2D_response = TH2D("H2D_response_"+partialUniqueSpecifier, "H2D_response_"+partialUniqueSpecifier, nBinPtJetsFine[iRadius], ptBinsJetsFine[iRadius], nBinPtJetsFine[iRadius], ptBinsJetsFine[iRadius]); // actually doesn't work if original histogram has fixed bin size

    //==================== Build response matrix: shift deltaPt by pT gen along the pT rec axis ====================//
    int ibinZeroFluct= H1D_fluctuations->FindBin(0+GLOBAL_epsilon);
    double integralError;
    // for(int iBinRec = 1; iBinRec <= H2D_response.GetNbinsX(); iBinRec++){
    for(int iBinRec = 0; iBinRec <= H2D_response.GetNbinsX()+1; iBinRec++){
      // for(int iBinGen = 1; iBinGen <= H2D_response.GetNbinsY(); iBinGen++){
      for(int iBinGen = 0; iBinGen <= H2D_response.GetNbinsY()+1; iBinGen++){
        double ptGen = H2D_response.GetYaxis()->GetBinLowEdge(iBinGen); // was bincenter before but then it'd give .5 values of GeV, and 
        double ptRec_low = H2D_response.GetXaxis()->GetBinLowEdge(iBinRec);
        double ptRec_up = H2D_response.GetXaxis()->GetBinLowEdge(iBinRec+1);
        // double xPtRecWidth = H2D_response->GetXaxis()->GetBinWidth(iBinRec);
        // if (ibinZero + (iBinRec - iBinGen) <= H1D_fluctuations_highRes->GetNbinsX()) { // make sure it's within bin range
        int iBin_fluct_low = H1D_fluctuations->GetXaxis()->FindBin(ptRec_low - ptGen + GLOBAL_epsilon);
        int iBin_fluct_high = H1D_fluctuations->GetXaxis()->FindBin(ptRec_up - ptGen - GLOBAL_epsilon);
        H2D_response.SetBinContent(iBinRec, iBinGen, H1D_fluctuations->IntegralAndError(iBin_fluct_low, iBin_fluct_high, integralError)); 
        // cout << "FluctResp(" << iBinRec << ", " << iBinGen << ") = " << H2D_response.GetBinContent(iBinRec, iBinGen) << endl;
        H2D_response.SetBinError(iBinRec, iBinGen, integralError); 
        // if (iBinRec == 0 && iBinGen == 0){
        //   cout << "iBinRec " << iBinRec << ", iBinGen "<< iBinGen << ": iBin_fluct_low = " << iBin_fluct_low << ", iBin_fluct_high = " << iBin_fluct_high << "             - ptRec_low-ptGen = " << ptRec_low - ptGen + GLOBAL_epsilon<< ", ptRec_up-ptGen = " << ptRec_up - ptGen + GLOBAL_epsilon<< endl;
        // }
      }
    }

    //========================================= Build response matrix end =========================================//

    H2D_jetPtResponseMatrix_fluctuations = (TH2D*)H2D_response.Clone("H2D_jetPtResponseMatrix_fluctuations"+partialUniqueSpecifier);
  }
}




//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////// Efficiency functions //////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


void  Get_ResponseMatrix_Pt_KinematicEffiency(TH1D* &H1D_kinematicEfficiency, TH2D* H2D_jetPtResponseMatrix_fineBinning, TString name_H1D_kinematicEfficiency, int iRadius){
  // assumes the response matrix has the fine binning, and will get the kinematic efficiency for rec axis binning equals to ptBinsJetsRec
  // cout << "Get_ResponseMatrix_Pt_KinematicEffiency" << endl; 
  
  int ibinRec_min = H2D_jetPtResponseMatrix_fineBinning->GetXaxis()->FindBin(ptBinsJetsRec[iRadius][0]);
  int ibinRec_max = H2D_jetPtResponseMatrix_fineBinning->GetXaxis()->FindBin(ptBinsJetsRec[iRadius][nBinPtJetsRec[iRadius]])-1;

  // cout << "ibinRec_min = " << ibinRec_min << ", ibinRec_max = " << ibinRec_max << ", ptmin = " << ptBinsJetsRec[iRadius][0] << ", ptmax = " << ptBinsJetsRec[iRadius][nBinPtJetsRec[iRadius]] << endl;

  TH1D* H1D_kinematicEfficiency_preRebin = H2D_jetPtResponseMatrix_fineBinning->ProjectionY("H1D_kinematicEfficiency_preRebin"+name_H1D_kinematicEfficiency, ibinRec_min, ibinRec_max, "e");
  
  H1D_kinematicEfficiency = (TH1D*)H1D_kinematicEfficiency_preRebin->Rebin(nBinPtJetsGen[iRadius],"H1D_kinematicEfficiency"+name_H1D_kinematicEfficiency+RadiusLegend[iRadius], ptBinsJetsGen[iRadius]);

  double integralOfResponse_iBinGen, integralOfResponse_iBinGen_error;
  double binContent, binError, binErrorA, binErrorB;
  for(int iBinGen = 1; iBinGen <= nBinPtJetsGen[iRadius]; iBinGen++){
    int ibinGen_low = H2D_jetPtResponseMatrix_fineBinning->GetYaxis()->FindBin(ptBinsJetsGen[iRadius][iBinGen-1]);
    int ibinGen_high = H2D_jetPtResponseMatrix_fineBinning->GetYaxis()->FindBin(ptBinsJetsGen[iRadius][iBinGen])-1;
    integralOfResponse_iBinGen = H2D_jetPtResponseMatrix_fineBinning->IntegralAndError( 1, nBinPtJetsFine[iRadius], ibinGen_low, ibinGen_high, integralOfResponse_iBinGen_error);


    H1D_kinematicEfficiency->GetBinContent(iBinGen) == 0 ? binErrorB = 0 : binErrorB = H1D_kinematicEfficiency->GetBinError(iBinGen)*H1D_kinematicEfficiency->GetBinError(iBinGen) / (H1D_kinematicEfficiency->GetBinContent(iBinGen)*H1D_kinematicEfficiency->GetBinContent(iBinGen));
    integralOfResponse_iBinGen == 0                                    ? binErrorA = 0 : binErrorA = integralOfResponse_iBinGen_error*integralOfResponse_iBinGen_error / (integralOfResponse_iBinGen*integralOfResponse_iBinGen);
    integralOfResponse_iBinGen == 0                                    ? binContent = 0 : binContent = H1D_kinematicEfficiency->GetBinContent(iBinGen) *1./integralOfResponse_iBinGen; // do I really give the value 0 if denominator is 0 ? 


    // H1D_kinematicEfficiency->SetBinContent(iBinGen, H1D_kinematicEfficiency->GetBinContent(iBinGen) * 1./H2D_jetPtResponseMatrix_fineBinning->Integral( 0, -1, ibinGen_low, ibinGen_high));
    // cout << "H1D_kinematicEfficiency_numerator(" << iBinGen << ") = " << H1D_kinematicEfficiency->GetBinContent(iBinGen) << ", error= = " << H1D_kinematicEfficiency->GetBinError(iBinGen) << ", ibinGen_low = " << ibinGen_low << ", ibinGen_high = " << ibinGen_high << "Integral divider = " << integralOfResponse_iBinGen << "Integral error = " << integralOfResponse_iBinGen_error << endl;
    H1D_kinematicEfficiency->SetBinContent(iBinGen, H1D_kinematicEfficiency->GetBinContent(iBinGen) * 1./integralOfResponse_iBinGen);
    H1D_kinematicEfficiency->SetBinError(iBinGen, sqrt(binContent*binContent * (binErrorA + binErrorB))); // sigma(A/B)2 / (A/B) = sigma(A)2 /A2 + sigma(B)2 /B2
    // cout << "H1D_kinematicEfficiency(" << iBinGen << ") = " << binContent << ", error= = " << H1D_kinematicEfficiency->GetBinError(iBinGen) << endl;

  }
}




bool  Get_Pt_JetEfficiency(TH1D* &H1D_jetEfficiency, int iDataset, int iRadius, float* centRange, std::string options){
  TH1D* H1D_jetPt_mcp;
  TH1D* H1D_jetPt_mcpMatched;
  bool divideSuccess;

  Get_Pt_spectrum_mcp(H1D_jetPt_mcp, iDataset, iRadius, false, options);
  Get_Pt_spectrum_mcpMatched(H1D_jetPt_mcpMatched, iDataset, iRadius, options);

  H1D_jetEfficiency = (TH1D*)H1D_jetPt_mcpMatched->Clone("H1D_jetEfficiency"+Datasets[iDataset]+"_R="+Form("%.1f", arrayRadius[iRadius])+"_@cent["+Form("%.1f", centRange[0])+","+Form("%.1f", centRange[1])+"]");
  divideSuccess = H1D_jetEfficiency->Divide(H1D_jetPt_mcp);
  if (!divideSuccess){
    cout << "################## Get_Pt_JetEfficiency FAILED!!!!! ##################" << endl;
  }
  return divideSuccess;
}
bool  Get_Pt_JetEfficiency_fineBinning(TH1D* &H1D_jetEfficiency, int iDataset, int iRadius, float* centRange, std::string options){
  TH1D* H1D_jetPt_mcp;
  TH1D* H1D_jetPt_mcpMatched;
  bool divideSuccess;

  Get_Pt_spectrum_mcp_fineBinning(H1D_jetPt_mcp, iDataset, iRadius, false, options);
  Get_Pt_spectrum_mcpMatched_fineBinning(H1D_jetPt_mcpMatched, iDataset, iRadius, options);

  H1D_jetEfficiency = (TH1D*)H1D_jetPt_mcpMatched->Clone("H1D_jetEfficiency_fineBinning"+Datasets[iDataset]+"_R="+Form("%.1f", arrayRadius[iRadius])+"_@cent["+Form("%.1f", centRange[0])+","+Form("%.1f", centRange[1])+"]");
  divideSuccess = H1D_jetEfficiency->Divide(H1D_jetPt_mcp);
  if (!divideSuccess){
    cout << "################## Get_Pt_JetEfficiency FAILED!!!!! ##################" << endl;
  }
  return divideSuccess;
}

bool Get_Pt_JetFakes(TH1D* &H1D_jetFakes, int iDataset, int iRadius, float* centRange, std::string options){
  TH1D* H1D_jetPt_mcd;
  TH1D* H1D_jetPt_mcdMatched;
  bool divideSuccess;

  Get_Pt_spectrum_mcd_recBinning(H1D_jetPt_mcd, iDataset, iRadius, options);
  Get_Pt_spectrum_mcdMatched_recBinning(H1D_jetPt_mcdMatched, iDataset, iRadius, options);


  H1D_jetFakes = (TH1D*)H1D_jetPt_mcdMatched->Clone("H1D_jetFakes"+Datasets[iDataset]+"_R="+Form("%.1f", arrayRadius[iRadius])+"_@cent["+Form("%.1f", centRange[0])+","+Form("%.1f", centRange[1])+"]");
  divideSuccess = H1D_jetFakes->Divide(H1D_jetPt_mcd);
  if (!divideSuccess){
    cout << "################## Get_Pt_JetFakes FAILED!!!!! ##################" << endl;
  }

  return divideSuccess;
}
bool Get_Pt_JetFakes_fineBinning(TH1D* &H1D_jetFakes, int iDataset, int iRadius, float* centRange, std::string options){
  TH1D* H1D_jetPt_mcd;
  TH1D* H1D_jetPt_mcdMatched;
  bool divideSuccess;

  Get_Pt_spectrum_mcd_fineBinning(H1D_jetPt_mcd, iDataset, iRadius, options);
  Get_Pt_spectrum_mcdMatched_fineBinning(H1D_jetPt_mcdMatched, iDataset, iRadius, options);

  H1D_jetFakes = (TH1D*)H1D_jetPt_mcdMatched->Clone("H1D_jetFakes_fineBinning"+Datasets[iDataset]+"_R="+Form("%.1f", arrayRadius[iRadius])+"_@cent["+Form("%.1f", centRange[0])+","+Form("%.1f", centRange[1])+"]");
  divideSuccess = H1D_jetFakes->Divide(H1D_jetPt_mcd);
  if (!divideSuccess){
    cout << "################## Get_Pt_JetFakes FAILED!!!!! ##################" << endl;
  }
  return divideSuccess;
}





//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////// RooUnfold Custom Utilities ///////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



int GetSvdBestRegularisationParameter(TSVDUnfold* unfoldTSvd){
  TH1D* H1D_D = unfoldTSvd->GetD();
  double minimum = 999999;
  double tempContent;
  int k = -99;
  for (int iBinX = 1; iBinX < H1D_D->GetNbinsX(); iBinX++) {
    tempContent = abs(abs(H1D_D->GetBinContent(iBinX)) - 1); // one wants |k| closest to 1 as possible
    if (minimum > tempContent) {
      minimum = tempContent;
      k = iBinX - 1; // or is it iBinX?
    }
  }
  cout << "k = " << k << ", |H1D-1| minimum = " << minimum << endl;
  return k;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////// Spectrum Unfolding functions ///////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


std::pair<int, RooUnfold*> Get_Pt_spectrum_unfolded_preWidthScalingAndEvtNorm(TH1D* &H1D_jetPt_unfolded, int iDataset, int iRadius, float* centRange, int unfoldParameterBayes, std::string options) {
  //for now makes the assumption gen and rec have the same pT binning

  TString partialUniqueSpecifier = Datasets[iDataset]+"_R="+Form("%.1f",arrayRadius[iRadius])+"_@cent["+Form("%.1f", centRange[0])+","+Form("%.1f", centRange[1])+"]";

  TH1D* measured;
  TH1D* mcp;
  TH1D* mcd;

  if (!normGenAndMeasByNEvts) {
    Get_Pt_spectrum_bkgCorrected_preWidthScalingAndEvtNorm(measured, iDataset, iRadius, centRange, options);
    Get_Pt_spectrum_mcp_preWidthScalingAndEvtNorm(mcp, iDataset, iRadius, false, options);
    Get_Pt_spectrum_mcd_recBinning_preWidthScalingAndEvtNorm(mcd, iDataset, iRadius, options);

    if (useFineBinningTest) {
      Get_Pt_spectrum_bkgCorrected_fineBinning_preWidthScalingAndEvtNorm(measured, iDataset, iRadius, centRange, options);
      Get_Pt_spectrum_mcp_fineBinning_preWidthScalingAndEvtNorm(mcp, iDataset, iRadius, false, options);
      Get_Pt_spectrum_mcd_fineBinning_preWidthScalingAndEvtNorm(mcd, iDataset, iRadius, options);
    }
  } else{
    Get_Pt_spectrum_bkgCorrected_preWidthScaling(measured, iDataset, iRadius, centRange, options);
    Get_Pt_spectrum_mcp_preWidthScaling(mcp, iDataset, iRadius, false, options);
    Get_Pt_spectrum_mcd_recBinning_preWidthScaling(mcd, iDataset, iRadius, options);

    if (useFineBinningTest) {
      Get_Pt_spectrum_bkgCorrected_fineBinning_preWidthScaling(measured, iDataset, iRadius, centRange, options);
      Get_Pt_spectrum_mcp_fineBinning_preWidthScaling(mcp, iDataset, iRadius, false, options);
      Get_Pt_spectrum_mcd_fineBinning_preWidthScaling(mcd, iDataset, iRadius, options);
    }
  } 

  bool divideSuccessFakes;
  TH1D* H1D_jetFakes;
  if (useInitialResponseMethod) {
    if (applyFakes) {
      if (!useFineBinningTest){ 
        divideSuccessFakes = Get_Pt_JetFakes(H1D_jetFakes, iDataset, iRadius, centRange, options);
      } else {
        divideSuccessFakes = Get_Pt_JetFakes_fineBinning(H1D_jetFakes, iDataset, iRadius, centRange, options);
      }
      if (divideSuccessFakes){
        measured->Multiply(H1D_jetFakes);
      } else {
        cout << "################## measured->Multiply(H1D_jetFakes) failed because Get_Pt_JetFakes() FAILED!!!!! ##################" << endl;
      }
    }
  }


  // TH1D* measured;
  // Get_Pt_spectrum_bkgCorrected(measured, iDataset, iRadius, centRange, options);
  // TH1D* mcp;
  // Get_Pt_spectrum_mcp(mcp, iDataset, iRadius, options);
  // TH1D* mcdMatched;
  // Get_Pt_spectrum_mcdMatched(mcdMatched, iDataset, iRadius, options);

  TH1D* H1D_kinematicEfficiency;
  TH2D* H2D_jetPtResponseMatrix_fluctuations;
  TH2D* H2D_jetPtResponseMatrix_detectorResponse;
  TH2D* H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined;
  TH2D* H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined_fineBinning;
  
  Get_PtResponseMatrix_Fluctuations(H2D_jetPtResponseMatrix_fluctuations, iDataset, iRadius, centRange);
  Get_PtResponseMatrix_detectorResponse(H2D_jetPtResponseMatrix_detectorResponse, iDataset, iRadius);
  // compute matrixFluctuations times matrixDetector

  Get_PtResponseMatrix_DetectorAndFluctuationsCombined_fineBinning(H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined_fineBinning, H2D_jetPtResponseMatrix_detectorResponse, H2D_jetPtResponseMatrix_fluctuations, iDataset, iRadius, centRange, options);
  Get_ResponseMatrix_Pt_KinematicEffiency(H1D_kinematicEfficiency, H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined_fineBinning, partialUniqueSpecifier, iRadius); // I want the efficiency before the reweighting and normalisation
  if (useFineBinningTest) {
    ReweightResponseMatrixWithPrior_fineBinning(H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined_fineBinning, iDataset, iRadius, centRange, options);
  }
  
  Get_PtResponseMatrix_DetectorAndFluctuationsCombined(H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined, H2D_jetPtResponseMatrix_detectorResponse, H2D_jetPtResponseMatrix_fluctuations, iDataset, iRadius, centRange, options);
  ReweightResponseMatrixWithPrior(H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined, iDataset, iRadius, centRange, options);

  if (drawIntermediateResponseMatrices) {
    TH2D* H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined_postWeighting = (TH2D*)H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined->Clone("H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined_postWeighting"+partialUniqueSpecifier);

    TString* pdfNamePost = new TString("responseMatrix_combined_postWeighting"+jetType[iJetType]+"_"+Datasets[iDataset]);
    TString* pdfNamePost_logz = new TString("responseMatrix_combined_postWeighting"+(TString)"_R="+Form("%.1f",arrayRadius[iRadius])+"_"+Datasets[iDataset]+"_logz");

    TString textContextPost(contextCustomOneField(*texDatasetsComparisonCommonDenominator, ""));

    Draw_TH2_Histogram(H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined_postWeighting, textContextPost, pdfNamePost, texPtJetRecX, texPtJetGenX, texCollisionDataInfo, drawnWindowAuto, legendPlacementAuto, contextPlacementAuto, "");
    Draw_TH2_Histogram(H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined_postWeighting, textContextPost, pdfNamePost_logz, texPtJetRecX, texPtJetGenX, texCollisionDataInfo, drawnWindowAuto, legendPlacementAuto, contextPlacementAuto, "logz");
  }

  cout << "RooUnfoldResponse setting - start" << endl;
  RooUnfoldResponse* response;
  if (useInitialResponseMethod) {

    // // based on Marta's work: https://twiki.cern.ch/twiki/bin/viewauth/ALICE/JEJetSpectrumUnfolding

    if (useFineBinningTest) {
      response = new RooUnfoldResponse(0, 0, H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined_fineBinning); // measured and mcp_rebinned are here to take inneficiencies and fakes into account; or is it really what's happening? 'Alternatively, the response matrix can be constructed from a pre-existing TH2D 2-dimensional histogram (with truth and measured distribution TH1D histograms for normalisation).' from https://hepunx.rl.ac.uk/~adye/software/unfold/RooUnfold.html, and I'm already normalising so maybe tehre's no need for more normalisation
    } else {
      response = new RooUnfoldResponse(0, 0, H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined); // measured and mcp_rebinned are here to take inneficiencies and fakes into account; or is it really what's happening? 'Alternatively, the response matrix can be constructed from a pre-existing TH2D 2-dimensional histogram (with truth and measured distribution TH1D histograms for normalisation).' from https://hepunx.rl.ac.uk/~adye/software/unfold/RooUnfold.html, and I'm already normalising so maybe tehre's no need for more normalisation
    }

    cout << "Get_Pt_spectrum_unfolded(): should I use response->UseOverflow() ? using it gives a ratio unfolded/mcp much higher than without using it" << endl;
    // response->UseOverflow();

  } else {
    TH2D *Respt;
    if (useFineBinningTest){
      Respt = H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined_fineBinning;
    } else {
      Respt = H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined;
    }
    TH1D *mcdMatched = (TH1D*)Respt->ProjectionX("Respt_projX", 1, Respt->GetNbinsY());
    TH1D *mcpMatched = (TH1D*)Respt->ProjectionY("Respt_projY", 1, Respt->GetNbinsX());
    TH1D *fake = (TH1D*)mcd->Clone();
    if (normDetRespByNEvts) {
      if (isPbPb || ppMcIsWeighted) {
        fake->Scale(1./GetNEventsSelected_JetFramework_weighted(file_O2Analysis_ppSimDetectorEffect));
      } else {
        fake->Scale(1./GetNEventsSelected_JetFramework(file_O2Analysis_ppSimDetectorEffect));
      }
    }
    for (auto i = 1; i <= fake->GetNbinsX(); i++) { 
      cout << "Fake i = " << i << " content: "<< fake->GetBinContent(i) << " - " << mcdMatched->GetBinContent(i) << " = " << fake->GetBinContent(i) - mcdMatched->GetBinContent(i) << " ------ error of mcd = " << fake->GetBinError(i) << ", error of mcdMatched = " << mcdMatched->GetBinError(i) << endl;
    }
    fake->Add(mcdMatched, -1);
    TH1D *miss = (TH1D*)mcp->Clone(); // for each bin of gen-level jet distribution, how many are not matched to a rec-level jet
    if (normDetRespByNEvts) {
      if (isPbPb || ppMcIsWeighted) {
        miss->Scale(1./GetNEventsSelected_JetFramework_weighted(file_O2Analysis_ppSimDetectorEffect));
      } else {
        miss->Scale(1./GetNEventsSelected_JetFramework(file_O2Analysis_ppSimDetectorEffect));
      }
    }
    for (auto i = 1; i <= fake->GetNbinsX(); i++) { 
      cout << "Miss i = " << i << " content: "<< miss->GetBinContent(i) << " - " << mcpMatched->GetBinContent(i) << " = " << miss->GetBinContent(i) - mcpMatched->GetBinContent(i) << " ------ error of mcp = " << miss->GetBinError(i) << ", error of mcpMatched = " << mcpMatched->GetBinError(i) << endl;
    }
    miss->Add(mcpMatched, -1); // for each bin of rec-level jet distribution, how many are not matched to a gen-level jet

    response = new RooUnfoldResponse(mcd, mcp);
    for (auto i = 1; i <= Respt->GetNbinsX(); i++) {
      for (auto j = 1; j <= Respt->GetNbinsY(); j++) { // ptpair
        Double_t bincenx = Respt->GetXaxis()->GetBinCenter(i);
        Double_t binceny = Respt->GetYaxis()->GetBinCenter(j);
        Double_t bincont = Respt->GetBinContent(i, j);
        response->Fill(bincenx, binceny, bincont);
        // cout << "response i,j = " << i << "," << j << " content: "<< bincont << endl;
      }
    }
    for (auto i = 1; i <= miss->GetNbinsX(); i++) { 
      Double_t bincenx = miss->GetXaxis()->GetBinCenter(i);
      Double_t bincont = miss->GetBinContent(i);
      response->Miss(bincenx, bincont);
      // cout << "Miss i = " << i << " content: "<< bincont << endl;
    }
    for (auto i = 1; i <= fake->GetNbinsX(); i++) {
      Double_t bincenx = fake->GetXaxis()->GetBinCenter(i);
      Double_t bincont = fake->GetBinContent(i);
      response->Fake(bincenx, bincont);
      // cout << "Fake i = " << i << " content: "<< bincont << endl;
    }


    // Get_ResponseMatrix_Pt_KinematicEffiency(H1D_kinematicEfficiency, (TH2D*)response->Hresponse(), partialUniqueSpecifier, iRadius); // I want the efficiency before the reweighting and normalisation

  }
  cout << "RooUnfoldResponse setting - end" << endl;




  // note about unfolding: if SUM(Mkj over k)=1 then it's also true for the inverse of M (can easily demosntrate it, just write it);
  RooUnfoldBayes* unfoldBayes = new RooUnfoldBayes(response, measured, unfoldParameterBayes);
  int unfoldParameterSvdInitial = 1;
  RooUnfoldSvd* unfoldSvdInitialiser = new RooUnfoldSvd(response, measured, unfoldParameterSvdInitial); // instance of RooUnfoldSvd only used to find the best regularisation parameter; the unfolded spectrum returned by it is not retrived
  RooUnfold* unfold = unfoldBayes; // default Bayes
  TH1D* hist_unfold;

  int unfoldParameter;
  if (options.find("Svd") != std::string::npos) {
    unfoldSvdInitialiser->Hreco(); // necessary to have GetD() give a meaningful output
    TSVDUnfold *tsvdUnfold = (TSVDUnfold*)unfoldSvdInitialiser->Impl();
    unfoldParameter = GetSvdBestRegularisationParameter(tsvdUnfold);
    RooUnfoldSvd* unfoldSvd = new RooUnfoldSvd(response, measured, unfoldParameter); // the RooUnfoldSvd instance that is actually used to unfold, with the best regularisation parameter

    // unfoldSvd->SetRegParm(-200);  // doesnt seem to change anything, even if I put 0 here ... annoying; but accoring to mailing list of roounfold it doesn't work well
    hist_unfold = (TH1D*)(unfoldSvd->Hreco());
    unfold = unfoldSvd;

    // plot svd d distribution
    TH1D* H1D_D = tsvdUnfold->GetD();
    TString* pdfName_regparam = new TString("Svd_regularisationd_distribution_"+(TString)"_R="+Form("%.1f",arrayRadius[iRadius])+"_"+Datasets[iDataset]);
    TString textContext(contextCustomOneField(*texDatasetsComparisonCommonDenominator, ""));
    Draw_TH1_Histogram(H1D_D, textContext, pdfName_regparam, texSvdDvector, texSvdK, texCollisionDataInfo, drawnWindowAuto, legendPlacementAuto, contextPlacementAuto, "logy");
  } else if (options.find("Bayes") != std::string::npos) {
    unfoldParameter = unfoldParameterBayes;
    hist_unfold = (TH1D*)(unfoldBayes->Hreco());
    unfold = unfoldBayes;
  }  


  H1D_jetPt_unfolded = (TH1D*)hist_unfold->Clone("H1D_jetPt_unfolded"+partialUniqueSpecifier);

  bool divideSuccessEff;
  TH1D* H1D_jetEfficiency;
  if (!useFineBinningTest){ 
    divideSuccessEff = Get_Pt_JetEfficiency(H1D_jetEfficiency, iDataset, iRadius, centRange, options);
  } else {
    divideSuccessEff = Get_Pt_JetEfficiency_fineBinning(H1D_jetEfficiency, iDataset, iRadius, centRange, options);
  }

  if (useInitialResponseMethod) {
    if (applyEfficiencies > 0) {
      if (applyEfficiencies == 2 || applyEfficiencies == 3) {
        if (divideSuccessEff){
          H1D_jetPt_unfolded->Divide(H1D_jetEfficiency);
        } else {
          cout << "################## Get_Pt_JetEfficiency FAILED!!!!! in Get_Pt_spectrum_unfolded_preWidthScalingAndEvtNorm ##################" << endl;
        }
      }
      if ((applyEfficiencies == 1 || applyEfficiencies == 3) && !useFineBinningTest) {
        H1D_jetPt_unfolded->Divide(H1D_kinematicEfficiency);
      }
    }
  }

  if (drawIntermediateResponseMatrices) {
    TH2D* H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined_postUnfolding = (TH2D*)unfold->response()->Hresponse()->Clone("H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined_postUnfolding"+partialUniqueSpecifier);

    TString* pdfName = new TString("responseMatrix_combined_postUnfolding"+jetType[iJetType]+"_"+Datasets[iDataset]);
    TString* pdfName_logz = new TString("responseMatrix_combined_postUnfolding"+(TString)"_R="+Form("%.1f",arrayRadius[iRadius])+"_"+Datasets[iDataset]+"_logz");

    TString textContext(contextCustomOneField(*texDatasetsComparisonCommonDenominator, ""));

    Draw_TH2_Histogram(H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined_postUnfolding, textContext, pdfName, texPtJetRecX, texPtJetGenX, texCollisionDataInfo, drawnWindowAuto, legendPlacementAuto, contextPlacementAuto, "");
    Draw_TH2_Histogram(H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined_postUnfolding, textContext, pdfName_logz, texPtJetRecX, texPtJetGenX, texCollisionDataInfo, drawnWindowAuto, legendPlacementAuto, contextPlacementAuto, "logz");
  }


  // if (doWidthScalingEarly) {
  //   TransformRawHistToYield(H1D_jetPt_unfolded);
  // }

  std::pair<int, RooUnfold*> unfoldInfo(unfoldParameter, unfold);
  return unfoldInfo;
}


std::pair<int, RooUnfold*> Get_Pt_spectrum_unfolded_preWidthScaling(TH1D* &H1D_jetPt_unfolded, int iDataset, int iRadius, float* centRange, int unfoldParameterBayes, std::string options) {
  std::pair<int, RooUnfold*> unfoldInfo = Get_Pt_spectrum_unfolded_preWidthScalingAndEvtNorm(H1D_jetPt_unfolded, iDataset, iRadius, centRange, unfoldParameterBayes, options);

  if (!normunfoldedByNEvts){
    if (options.find("evtNorm") != std::string::npos && doEvtNorm) {
      if (isPbPb) {
        NormaliseRawHistToNEvents(H1D_jetPt_unfolded, GetNEventsSelectedCentrality_JetFramework(file_O2Analysis_list[iDataset], centRange[0], centRange[1], trainId));
      } else {
        if (!controlMC) {
          NormaliseRawHistToNEvents(H1D_jetPt_unfolded, GetNEventsSelectedCentrality_JetFramework(file_O2Analysis_list[iDataset], -999, 999, trainId));
        } else {
          NormaliseRawHistToNEvents(H1D_jetPt_unfolded, GetNEventsSelectedCentrality_JetFramework(file_O2Analysis_ppSimDetectorEffect_unfoldingControl, -999, 999, trainId));
        }
      }
    }
    if (options.find("entriesNorm") != std::string::npos) {
      NormaliseYieldToNEntries(H1D_jetPt_unfolded);
    }
  }
  
  return unfoldInfo;
}
std::pair<int, RooUnfold*> Get_Pt_spectrum_unfolded(TH1D* &H1D_jetPt_unfolded, int iDataset, int iRadius, float* centRange, int unfoldParameterBayes, std::string options) {
  std::pair<int, RooUnfold*> unfoldInfo = Get_Pt_spectrum_unfolded_preWidthScaling(H1D_jetPt_unfolded, iDataset, iRadius, centRange, unfoldParameterBayes, options);

  if (doWidthScalingAtEnd) {
    TransformRawHistToYield(H1D_jetPt_unfolded);
  }
  return unfoldInfo;
}



void Get_Pt_spectrum_unfoldedThenRefolded_preWidthScalingAndEvtNorm(TH1D* &H1D_jetPt_unfoldedThenRefolded, int iDataset, int iRadius, float* centRange, int unfoldParameterBayes, std::string options) {
  TH1D* H1D_jetPt_unfolded;
  // TH1D* H1D_jetPt_raw[nRadius];
  TString partialUniqueSpecifier = Datasets[iDataset]+"_R="+Form("%.1f",arrayRadius[iRadius])+"_@cent["+Form("%.1f", centRange[0])+","+Form("%.1f", centRange[1])+"]";


  RooUnfold* unfold = Get_Pt_spectrum_unfolded_preWidthScalingAndEvtNorm(H1D_jetPt_unfolded, iDataset, iRadius, centRange, unfoldParameterBayes, options).second;
  // Get_Pt_spectrum_unfolded(H1D_jetPt_unfolded, iDataset, iRadius, centRange, unfoldParameterBayes, options);

  // cout << "((((((((((((((((()))))))))))))))))" << endl;
  // cout << "REFOLDING TEST: " << endl;
  // for(int iBinX = 0; iBinX <= H1D_jetPt_unfolded->GetNbinsX()+1; iBinX++){
  //   cout << "H1D_jetPt_unfolded(" << iBinX << ") = " << H1D_jetPt_unfolded->GetBinContent(iBinX) << ", error = "<< H1D_jetPt_unfolded->GetBinError(iBinX) << endl;
  // }

  TH2D* H2D_jetPtResponseMatrix_fluctuations;
  TH2D* H2D_jetPtResponseMatrix_detectorResponse;
  TH2D* H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined;
  TH2D* H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined_fineBinning;
  TH1D* H1D_kinematicEfficiency;

  Get_PtResponseMatrix_Fluctuations(H2D_jetPtResponseMatrix_fluctuations, iDataset, iRadius, centRange);
  Get_PtResponseMatrix_detectorResponse(H2D_jetPtResponseMatrix_detectorResponse, iDataset, iRadius);



  Get_PtResponseMatrix_DetectorAndFluctuationsCombined_fineBinning(H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined_fineBinning, H2D_jetPtResponseMatrix_detectorResponse, H2D_jetPtResponseMatrix_fluctuations, iDataset, iRadius, centRange, options);
  Get_ResponseMatrix_Pt_KinematicEffiency(H1D_kinematicEfficiency, H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined_fineBinning, partialUniqueSpecifier, iRadius); // I want the efficiency before the reweighting and normalisation
  if (useFineBinningTest) {
    ReweightResponseMatrixWithPrior_fineBinning(H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined_fineBinning, iDataset, iRadius, centRange, options);
  }

  Get_PtResponseMatrix_DetectorAndFluctuationsCombined(H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined, H2D_jetPtResponseMatrix_detectorResponse, H2D_jetPtResponseMatrix_fluctuations, iDataset, iRadius, centRange, options);
  ReweightResponseMatrixWithPrior(H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined, iDataset, iRadius, centRange, options);
  

  if (drawIntermediateResponseMatrices) {
    TH2D* H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined_duringRefolding = (TH2D*)H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined->Clone("H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined_duringRefolding"+partialUniqueSpecifier);

    TString* pdfName = new TString("responseMatrix_combined_duringRefolding"+jetType[iJetType]+"_"+Datasets[iDataset]);
    TString* pdfName_logz = new TString("responseMatrix_combined_duringRefolding"+(TString)"_R="+Form("%.1f",arrayRadius[iRadius])+"_"+Datasets[iDataset]+"_logz");

    TString textContext(contextCustomOneField(*texDatasetsComparisonCommonDenominator, ""));

    Draw_TH2_Histogram(H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined_duringRefolding, textContext, pdfName, texPtJetRecX, texPtJetGenX, texCollisionDataInfo, drawnWindowAuto, legendPlacementAuto, contextPlacementAuto, "");
    Draw_TH2_Histogram(H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined_duringRefolding, textContext, pdfName_logz, texPtJetRecX, texPtJetGenX, texCollisionDataInfo, drawnWindowAuto, legendPlacementAuto, contextPlacementAuto, "logz");
  }

  cout << "((((((((((((((((()))))))))))))))))" << endl;
  cout << "REFOLDING TEST: pre efficiency" << endl;
  for(int iBinX = 0; iBinX <= H1D_jetPt_unfolded->GetNbinsX()+1; iBinX++){
    cout << "H1D_jetPt_unfolded(" << iBinX << ") = " << H1D_jetPt_unfolded->GetBinContent(iBinX) << endl;
  }
  

  bool divideSuccessEff;
  TH1D* H1D_jetEfficiency;
  if (!useFineBinningTest){ 
    divideSuccessEff = Get_Pt_JetEfficiency(H1D_jetEfficiency, iDataset, iRadius, centRange, options);
  } else {
    divideSuccessEff = Get_Pt_JetEfficiency_fineBinning(H1D_jetEfficiency, iDataset, iRadius, centRange, options);
  }

  if (applyEfficiencies > 0) {
    if (applyEfficiencies == 2 || applyEfficiencies == 3) {
      if (divideSuccessEff){
        H1D_jetPt_unfolded->Multiply(H1D_jetEfficiency);
      } else {
        cout << "################## H1D_jetPt_unfolded->Multiply(H1D_jetEfficiency) failed because Get_Pt_JetEfficiency() FAILED!!!!! ##################" << endl;
      }
    }
    if ((applyEfficiencies == 1 || applyEfficiencies == 3) && !useFineBinningTest) {
      H1D_jetPt_unfolded->Multiply(H1D_kinematicEfficiency);
    }
  }

  TH2D* refoldingResponseMatrix;
  if (useInitialResponseMethod) {
    if (useFineBinningTest) {
      refoldingResponseMatrix = (TH2D*)H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined_fineBinning->Clone("refoldingResponseMatrix"+partialUniqueSpecifier);
    } else {
      refoldingResponseMatrix = (TH2D*)H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined->Clone("refoldingResponseMatrix"+partialUniqueSpecifier);
    }
  } else {
    refoldingResponseMatrix = (TH2D*)unfold->response()->HresponseNoOverflow();
    // refoldingResponseMatrix = (TH2D*)unfold->response()->Hresponse();
    cout << "refold why use HresponseNoOverflow instead of Hresponse?" << endl;
  }

  if (normaliseRespYSliceForRefold){
    NormaliseYSlicesToOne(refoldingResponseMatrix);
  }
  // H1D_jetPt_unfoldedThenRefolded = (TH1D*)GetMatrixVectorProductTH2xTH1((TH2D*)unfold->response()->Hresponse(), H1D_jetPt_unfolded).Clone("Get_Pt_spectrum_unfoldedThenRefolded_preWidthScalingAndEvtNorm"+partialUniqueSpecifier); //gives the same
  // H1D_jetPt_unfoldedThenRefolded = (TH1D*)GetMatrixVectorProductTH2xTH1(H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined, H1D_jetPt_unfolded).Clone("Get_Pt_spectrum_unfoldedThenRefolded_preWidthScalingAndEvtNorm"+partialUniqueSpecifier); //gives the same
  H1D_jetPt_unfoldedThenRefolded = (TH1D*)GetMatrixVectorProductTH2xTH1(refoldingResponseMatrix, H1D_jetPt_unfolded).Clone("Get_Pt_spectrum_unfoldedThenRefolded_preWidthScalingAndEvtNorm"+partialUniqueSpecifier); //gives the same


  bool divideSuccessFakes;
  TH1D* H1D_jetFakes;
  if (applyFakes) { // for useInitialResponseMethod set to false, it will give a slightly different result as Get_Pt_JetFakes doesn't get fakes from the same histogram as the one used to encode fakes in the response matrix
    if (!useFineBinningTest){ 
      divideSuccessFakes = Get_Pt_JetFakes(H1D_jetFakes, iDataset, iRadius, centRange, options);
    } else {
      divideSuccessFakes = Get_Pt_JetFakes_fineBinning(H1D_jetFakes, iDataset, iRadius, centRange, options);
    }
    if (divideSuccessFakes){
      H1D_jetPt_unfoldedThenRefolded->Divide(H1D_jetFakes);
    } else {
      cout << "################## H1D_jetPt_unfoldedThenRefolded->Divide(H1D_jetFakes) failed because Get_Pt_JetFakes() FAILED!!!!! ##################" << endl;
    }
  }

  // cout << "((((((((((((((((()))))))))))))))))" << endl;
  // cout << "REFOLDING TEST: post efficiency" << endl;
  // for(int iBinX = 0; iBinX <= H1D_jetPt_unfoldedThenRefolded->GetNbinsX()+1; iBinX++){
  //   cout << "H1D_jetPt_unfoldedThenRefolded(" << iBinX << ") = " << H1D_jetPt_unfoldedThenRefolded->GetBinContent(iBinX) << ", error = "<< H1D_jetPt_unfoldedThenRefolded->GetBinError(iBinX) << endl;
  // }
  

  // cout << "Refolding check of Combined matrix = " << endl; 
  // for(int iBinY = 1; iBinY <= H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined->GetNbinsY(); iBinY++){
  //   cout << "combinedMatrix integral of slice iBinY = " << iBinY << " is: " << H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined->Integral(0, -1, iBinY, iBinY) << endl;
  // }


  // if (useFineBinningTest) {
  //   H1D_jetPt_unfoldedThenRefolded = (TH1D*)GetMatrixVectorProductTH2xTH1(H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined_fineBinning, H1D_jetPt_unfolded).Clone("Get_Pt_spectrum_unfoldedThenRefolded_preWidthScalingAndEvtNorm"+partialUniqueSpecifier);
  // } else {
  //   H1D_jetPt_unfoldedThenRefolded = (TH1D*)GetMatrixVectorProductTH2xTH1(H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined, H1D_jetPt_unfolded).Clone("Get_Pt_spectrum_unfoldedThenRefolded_preWidthScalingAndEvtNorm"+partialUniqueSpecifier);
  // }
  

  // if (doWidthScalingEarly) {
  //   TransformRawHistToYield(H1D_jetPt_unfoldedThenRefolded);
  // }

  // cout << "still not giving back the measured used as input to the unfolding; got an issue somewhere" << endl;
}

void Get_Pt_spectrum_unfoldedThenRefolded_preWidthScaling(TH1D* &H1D_jetPt_unfoldedThenRefolded, int iDataset, int iRadius, float* centRange, int unfoldParameterBayes, std::string options) {
  Get_Pt_spectrum_unfoldedThenRefolded_preWidthScalingAndEvtNorm(H1D_jetPt_unfoldedThenRefolded, iDataset, iRadius, centRange, unfoldParameterBayes, options);

  if (!normunfoldedByNEvts){
    if (options.find("evtNorm") != std::string::npos && doEvtNorm) {
        if (isPbPb) {
        NormaliseRawHistToNEvents(H1D_jetPt_unfoldedThenRefolded, GetNEventsSelectedCentrality_JetFramework(file_O2Analysis_list[iDataset], centRange[0], centRange[1], trainId));
      } else {
        // NormaliseRawHistToNEvents(H1D_jetPt_unfoldedThenRefolded, GetNEventsSelectedCentrality_JetFramework(file_O2Analysis_list[iDataset], 0, -1, trainId));
        if (!controlMC) {
          NormaliseRawHistToNEvents(H1D_jetPt_unfoldedThenRefolded, GetNEventsSelectedCentrality_JetFramework(file_O2Analysis_list[iDataset], -999, 999, trainId));
        } else {
          NormaliseRawHistToNEvents(H1D_jetPt_unfoldedThenRefolded, GetNEventsSelectedCentrality_JetFramework(file_O2Analysis_ppSimDetectorEffect_unfoldingControl, -999, 999, trainId));
        }
      }
    }
    if (options.find("entriesNorm") != std::string::npos) {
      NormaliseYieldToNEntries(H1D_jetPt_unfoldedThenRefolded);
    }
  }
}
void Get_Pt_spectrum_unfoldedThenRefolded(TH1D* &H1D_jetPt_unfoldedThenRefolded, int iDataset, int iRadius, float* centRange, int unfoldParameterBayes, std::string options) {
  Get_Pt_spectrum_unfoldedThenRefolded_preWidthScaling(H1D_jetPt_unfoldedThenRefolded, iDataset, iRadius, centRange, unfoldParameterBayes, options);

  if (doWidthScalingAtEnd) {
    TransformRawHistToYield(H1D_jetPt_unfoldedThenRefolded);
  }
}


void Get_Pt_spectrum_unfoldedThenRefolded_RooUnfoldMethod_preWidthScalingAndEvtNorm(TH1D* &H1D_jetPt_unfoldedThenRefolded, int iDataset, int iRadius, float* centRange, int unfoldParameterBayes, std::string options) {
  // Matches exactly with the manual method IF NO PRIOR
  // if I have a non flat prior, then the roounfold method gives me a good closure test, but not the manual method!

  TH1D* H1D_jetPt_unfolded;
  // TH1D* H1D_jetPt_raw[nRadius];
  TString partialUniqueSpecifier = Datasets[iDataset]+"_R="+Form("%.1f",arrayRadius[iRadius])+"_@cent["+Form("%.1f", centRange[0])+","+Form("%.1f", centRange[1])+"]";


  RooUnfold* unfold = Get_Pt_spectrum_unfolded_preWidthScalingAndEvtNorm(H1D_jetPt_unfolded, iDataset, iRadius, centRange, unfoldParameterBayes, options).second;

  TH2D* H2D_jetPtResponseMatrix_fluctuations;
  TH2D* H2D_jetPtResponseMatrix_detectorResponse;
  TH2D* H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined;
  TH2D* H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined_fineBinning;
  TH1D* H1D_kinematicEfficiency;

  Get_PtResponseMatrix_Fluctuations(H2D_jetPtResponseMatrix_fluctuations, iDataset, iRadius, centRange);
  Get_PtResponseMatrix_detectorResponse(H2D_jetPtResponseMatrix_detectorResponse, iDataset, iRadius);

  Get_PtResponseMatrix_DetectorAndFluctuationsCombined_fineBinning(H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined_fineBinning, H2D_jetPtResponseMatrix_detectorResponse, H2D_jetPtResponseMatrix_fluctuations, iDataset, iRadius, centRange, options);
  Get_ResponseMatrix_Pt_KinematicEffiency(H1D_kinematicEfficiency, H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined_fineBinning, partialUniqueSpecifier, iRadius); // I want the efficiency before the reweighting and normalisation

  Get_PtResponseMatrix_DetectorAndFluctuationsCombined(H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined, H2D_jetPtResponseMatrix_detectorResponse, H2D_jetPtResponseMatrix_fluctuations, iDataset, iRadius, centRange, options);
  ReweightResponseMatrixWithPrior(H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined, iDataset, iRadius, centRange, options);

  bool divideSuccessEff;
  TH1D* H1D_jetEfficiency;
  if (!useFineBinningTest){ 
    divideSuccessEff = Get_Pt_JetEfficiency(H1D_jetEfficiency, iDataset, iRadius, centRange, options);
  } else {
    divideSuccessEff = Get_Pt_JetEfficiency_fineBinning(H1D_jetEfficiency, iDataset, iRadius, centRange, options);
  }

  if (useInitialResponseMethod) { // this condition isn't present in the manual refolding method as in the manual case the response matrix object is just a th2 and not a roounfoldresponse object and doesn't have the efficiencies/fakes encoded in it
    if (applyEfficiencies > 0) {
      if (applyEfficiencies == 2 || applyEfficiencies == 3) {
        if (divideSuccessEff){
          H1D_jetPt_unfolded->Multiply(H1D_jetEfficiency);
        } else {
          cout << "################## H1D_jetPt_unfolded->Multiply(H1D_jetEfficiency) failed because Get_Pt_JetEfficiency() FAILED!!!!! ##################" << endl;
        }
      }
      if ((applyEfficiencies == 1 || applyEfficiencies == 3) && !useFineBinningTest) {
        H1D_jetPt_unfolded->Multiply(H1D_kinematicEfficiency);
      }
    }
  }

  // cout << "Refolding check of Combined matrix = " << endl; 
  // for(int iBinY = 1; iBinY <= H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined->GetNbinsY(); iBinY++){
  //   cout << "combinedMatrix integral of slice iBinY = " << iBinY << " is: " << H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined->Integral(0, -1, iBinY, iBinY) << endl;
  // }


  // RooUnfoldResponse* rooresponse_postUnfold = unfold->response();
  // H1D_jetPt_unfoldedThenRefolded = (TH1D*)unfold->response()->ApplyToTruth(H1D_jetPt_unfolded, "H1D_jetPt_unfoldedThenRefolded_withApplyToTruth")->Clone("Get_Pt_spectrum_unfoldedThenRefolded_preWidthScalingAndEvtNorm"+partialUniqueSpecifier);

  // RooUnfoldResponse* response = new RooUnfoldResponse(0, 0, H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined); // measured and mcp_rebinned are here to take inneficiencies and fakes into account; or is it really what's happening? 'Alternatively, the response matrix can be constructed from a pre-existing TH2D 2-dimensional histogram (with truth and measured distribution TH1D histograms for normalisation).' from https://hepunx.rl.ac.uk/~adye/software/unfold/RooUnfold.html, and I'm already normalising so maybe tehre's no need for more normalisation
  H1D_jetPt_unfoldedThenRefolded = (TH1D*)unfold->response()->ApplyToTruth(H1D_jetPt_unfolded, "H1D_jetPt_unfoldedThenRefolded_withApplyToTruth")->Clone("Get_Pt_spectrum_unfoldedThenRefolded_preWidthScalingAndEvtNorm"+partialUniqueSpecifier);


  bool divideSuccessFakes;
  TH1D* H1D_jetFakes;
  if (useInitialResponseMethod) { // this condition isn't present in the manual refolding method as in the manual case the response matrix object is just a th2 and not a roounfoldresponse object and doesn't have the efficiencies/fakes encoded in it
    if (applyFakes) {
      if (!useFineBinningTest){ 
        divideSuccessFakes = Get_Pt_JetFakes(H1D_jetFakes, iDataset, iRadius, centRange, options);
      } else {
        divideSuccessFakes = Get_Pt_JetFakes_fineBinning(H1D_jetFakes, iDataset, iRadius, centRange, options);
      }
      if (divideSuccessFakes){
        H1D_jetPt_unfoldedThenRefolded->Divide(H1D_jetFakes);
      } else {
        cout << "################## H1D_jetPt_unfoldedThenRefolded->Divide(H1D_jetFakes) failed because Get_Pt_JetFakes() FAILED!!!!! ##################" << endl;
      }
    }
  }

  //not sure I should normalise --> probably not as H1D_jetPt_unfolded is already normalised in Get_Pt_spectrum_unfolded_preWidthScalingAndEvtNorm

  // cout << "still not giving back the measured used as input to the unfolding; got an issue somewhere" << endl;
}

void Get_Pt_spectrum_unfoldedThenRefolded_RooUnfoldMethod_preWidthScaling(TH1D* &H1D_jetPt_unfoldedThenRefolded, int iDataset, int iRadius, float* centRange, int unfoldParameterBayes, std::string options) {
  Get_Pt_spectrum_unfoldedThenRefolded_RooUnfoldMethod_preWidthScalingAndEvtNorm(H1D_jetPt_unfoldedThenRefolded, iDataset, iRadius, centRange, unfoldParameterBayes, options);

  if (!normunfoldedByNEvts){
    if (options.find("evtNorm") != std::string::npos && doEvtNorm) {
        if (isPbPb) {
        NormaliseRawHistToNEvents(H1D_jetPt_unfoldedThenRefolded, GetNEventsSelectedCentrality_JetFramework(file_O2Analysis_list[iDataset], centRange[0], centRange[1], trainId));
      } else {
        // NormaliseRawHistToNEvents(H1D_jetPt_unfoldedThenRefolded, GetNEventsSelectedCentrality_JetFramework(file_O2Analysis_list[iDataset], -999, 999, trainId));
        if (!controlMC) {
          NormaliseRawHistToNEvents(H1D_jetPt_unfoldedThenRefolded, GetNEventsSelectedCentrality_JetFramework(file_O2Analysis_list[iDataset], -999, 999, trainId));
        } else {
          NormaliseRawHistToNEvents(H1D_jetPt_unfoldedThenRefolded, GetNEventsSelectedCentrality_JetFramework(file_O2Analysis_ppSimDetectorEffect_unfoldingControl, -999, 999, trainId));
        }
      }
    }
    if (options.find("entriesNorm") != std::string::npos) {
      NormaliseYieldToNEntries(H1D_jetPt_unfoldedThenRefolded);
    }
  }
}
void Get_Pt_spectrum_unfoldedThenRefolded_RooUnfoldMethod(TH1D* &H1D_jetPt_unfoldedThenRefolded, int iDataset, int iRadius, float* centRange, int unfoldParameterBayes, std::string options) {
  Get_Pt_spectrum_unfoldedThenRefolded_RooUnfoldMethod_preWidthScaling(H1D_jetPt_unfoldedThenRefolded, iDataset, iRadius, centRange, unfoldParameterBayes, options);

  if (doWidthScalingAtEnd) {
    TransformRawHistToYield(H1D_jetPt_unfoldedThenRefolded);
  }
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////// Spectrum plotting functions ///////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


void Draw_Pt_spectrum_raw(int iDataset, int iRadius, std::string options) {

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
  yAxisLabel = texCount;
  if (options.find("evtNorm") != std::string::npos && doEvtNorm) {
    yAxisLabel = texJetPtYield_EventNorm;
  }
  if (options.find("entriesNorm") != std::string::npos) {
    yAxisLabel = texJetPtYield_EntriesNorm;
  }

  // Draw_TH1_Histograms_in_one(H1D_jetPt_raw, centralityLegend, nCentralityBins, textContext, pdfName, texPtX, yAxisLabel, texCollisionDataInfo, "logy");
  Draw_TH1_Histograms_in_one(H1D_jetPt_raw, centralityLegend, nCentralityBins, textContext, pdfName, texPtX, yAxisLabel, texCollisionDataInfo, drawnWindowAuto, legendPlacementAuto, contextPlacementAuto, "logy");
}

void Draw_Pt_spectrum_mcp(int iDataset, int iRadius, std::string options) {

  TH1D* H1D_jetPt_mcp_genBinning;
  TH1D* H1D_jetPt_mcp_recBinning;
  TH1D* H1D_jetPt_mcp_collection[2];
  Get_Pt_spectrum_mcp(H1D_jetPt_mcp_genBinning, iDataset, iRadius, controlMC, options);
  Get_Pt_spectrum_mcp_recBinning(H1D_jetPt_mcp_recBinning, iDataset, iRadius, controlMC, options);
  H1D_jetPt_mcp_collection[0] = H1D_jetPt_mcp_genBinning;
  H1D_jetPt_mcp_collection[1] = H1D_jetPt_mcp_recBinning;


  TString* pdfName = new TString("jet_"+jetType[iJetType]+"_"+jetLevel[iJetLevel]+"_"+Datasets[iDataset]+"_R="+Form("%.1f", arrayRadius[iRadius])+"_Pt_mcp");

  TString textContext(contextCustomOneField(*texDatasetsComparisonCommonDenominator, ""));

  TString* yAxisLabel;
  yAxisLabel = texCount;
  if (options.find("evtNorm") != std::string::npos && doEvtNorm) {
    yAxisLabel = texJetPtYield_EventNorm;
  }
  if (options.find("entriesNorm") != std::string::npos) {
    yAxisLabel = texJetPtYield_EntriesNorm;
  }
  TString genVsRecBinningLegend[2] = {"gen binning", "rec binning"};

  Draw_TH1_Histograms_in_one(H1D_jetPt_mcp_collection, genVsRecBinningLegend, 2, textContext, pdfName, texPtX, yAxisLabel, texCollisionDataInfo, drawnWindowAuto, legendPlacementAuto, contextPlacementAuto, "logy"); 
}

void Draw_Pt_spectrum_mcdMatched(int iDataset, int iRadius, std::string options) {

  TH1D* H1D_jetPt_mcdMatched;
  Get_Pt_spectrum_mcdMatched(H1D_jetPt_mcdMatched, iDataset, iRadius, options);


  TString* pdfName = new TString("jet_"+jetType[iJetType]+"_"+jetLevel[iJetLevel]+"_"+Datasets[iDataset]+"_R="+Form("%.1f", arrayRadius[iRadius])+"_Pt_mcdMatched");

  TString textContext(contextCustomOneField(*texDatasetsComparisonCommonDenominator, ""));

  TString* yAxisLabel;
  yAxisLabel = texCount;
  if (options.find("evtNorm") != std::string::npos && doEvtNorm) {
    yAxisLabel = texJetPtYield_EventNorm;
  }
  if (options.find("entriesNorm") != std::string::npos) {
    yAxisLabel = texJetPtYield_EntriesNorm;
  }

  Draw_TH1_Histogram(H1D_jetPt_mcdMatched, textContext, pdfName, texPtX, yAxisLabel, texCollisionDataInfo, drawnWindowAuto, legendPlacementAuto, contextPlacementAuto, "logy");
}


void Draw_Pt_efficiency_jets(int iDataset, int iRadius, std::string options) {
  TH1D* H1D_jetEfficiency[nCentralityBins];
  bool divideSuccess[nCentralityBins];

  float centRange[2];
  for(int iCent = 0; iCent < nCentralityBins; iCent++){
    centRange[0] = arrayCentralityIntervals[iCent][0];
    centRange[1] = arrayCentralityIntervals[iCent][1];
    divideSuccess[iCent] = Get_Pt_JetEfficiency(H1D_jetEfficiency[iCent], iDataset, iRadius, centRange, options);
  }
  TString* pdfName = new TString("jet_"+jetType[iJetType]+"_"+jetLevel[iJetLevel]+"_"+Datasets[iDataset]+"_R="+Form("%.1f", arrayRadius[iRadius])+"_Pt_efficiency");

  TString textContext(contextCustomOneField(*texDatasetsComparisonCommonDenominator, ""));

  TString centralityLegend[nCentralityBins]; CentralityLegend(centralityLegend, arrayCentralityIntervals, nCentralityBins);


  bool divideSuccess_boolsum = true;
  for(int iCent = 0; iCent < nCentralityBins; iCent++){
    if (!(divideSuccess[iCent])) {
      divideSuccess_boolsum = false;
    }
  }
  if (divideSuccess_boolsum){
    Draw_TH1_Histograms_in_one(H1D_jetEfficiency, centralityLegend, nCentralityBins, textContext, pdfName, texPtJetGenX, texJetEfficiency, texCollisionDataInfo, drawnWindowAuto, legendPlacementAuto, contextPlacementAuto, "efficiency");
  }
}

void Draw_kinematicEfficiency(int iDataset, int iRadius, std::string options) {

  TH2D* H2D_jetPtResponseMatrix_fluctuations[nCentralityBins];
  TH2D* H2D_jetPtResponseMatrix_detectorResponse;
  TH2D* H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined_fineBinning[nCentralityBins];
  TH1D* H1D_kinematicEfficiency[nCentralityBins];

  TString partialUniqueSpecifier = Datasets[iDataset]+"_R="+Form("%.1f",arrayRadius[iRadius]);
  TString name_H1D_kinematicEfficiency = Datasets[iDataset]+"_R="+Form("%.1f",arrayRadius[iRadius]);

  float centRange[2];
  Get_PtResponseMatrix_detectorResponse(H2D_jetPtResponseMatrix_detectorResponse, iDataset, iRadius);
  for(int iCent = 0; iCent < nCentralityBins; iCent++){
    centRange[0] = arrayCentralityIntervals[iCent][0];
    centRange[1] = arrayCentralityIntervals[iCent][1];

    Get_PtResponseMatrix_Fluctuations(H2D_jetPtResponseMatrix_fluctuations[iCent], iDataset, iRadius, centRange);
    Get_PtResponseMatrix_DetectorAndFluctuationsCombined_fineBinning(H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined_fineBinning[iCent], H2D_jetPtResponseMatrix_detectorResponse, H2D_jetPtResponseMatrix_fluctuations[iCent], iDataset, iRadius, centRange, options);

    name_H1D_kinematicEfficiency += (TString)"_@cent["+Form("%.1f", centRange[0])+","+Form("%.1f", centRange[1])+"]";

    Get_ResponseMatrix_Pt_KinematicEffiency(H1D_kinematicEfficiency[iCent], H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined_fineBinning[iCent], name_H1D_kinematicEfficiency, iRadius);
  }
  TString* pdfName = new TString("kinematicEfficiency_"+partialUniqueSpecifier+"_centralityComp");
  // TString* pdfName_logz = new TString("kinematicEfficiency_"+(TString)"_R="+Form("%.1f",arrayRadius[iRadius])+"_"+Datasets[iDataset]+"_centralityComp_logz");

  TString textContext(contextCustomOneField(*texDatasetsComparisonCommonDenominator, ""));

  TString centralityLegend[nCentralityBins]; CentralityLegend(centralityLegend, arrayCentralityIntervals, nCentralityBins);

  // Draw_TH2_Histograms(H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined, centralityLegend, nCentralityBins, textContext, pdfName, texPtJetRecX, texPtJetGenX, texCollisionDataInfo, "logz");

  Draw_TH1_Histograms_in_one(H1D_kinematicEfficiency, centralityLegend, nCentralityBins, textContext, pdfName, texPtJetGenX, texJetKinematicEfficiency, texCollisionDataInfo, drawnWindowAuto, legendPlacementAuto, contextPlacementAuto, "");
  // Draw_TH2_Histograms(H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined, centralityLegend, nCentralityBins, textContext, pdfName_logz, texPtJetRecX, texPtJetGenX, texCollisionDataInfo, drawnWindowAuto, legendPlacementAuto, contextPlacementAuto, "logz");
}


void Draw_FakeRatio(int iDataset, int iRadius, std::string options) {
  TH1D* H1D_fakeRatio[nCentralityBins];

  TString partialUniqueSpecifier = Datasets[iDataset]+"_R="+Form("%.1f",arrayRadius[iRadius]);
  TString name_H1D_fakeRatio = Datasets[iDataset]+"_R="+Form("%.1f",arrayRadius[iRadius]);

  float centRange[2];
  for(int iCent = 0; iCent < nCentralityBins; iCent++){
    centRange[0] = arrayCentralityIntervals[iCent][0];
    centRange[1] = arrayCentralityIntervals[iCent][1];

    name_H1D_fakeRatio += (TString)"_@cent["+Form("%.1f", centRange[0])+","+Form("%.1f", centRange[1])+"]";

    Get_Pt_JetFakes(H1D_fakeRatio[iCent], iDataset, iRadius, centRange, options);
  }
  TString* pdfName = new TString("fakeRatio_"+partialUniqueSpecifier+"_centralityComp");
  // TString* pdfName_logz = new TString("kinematicEfficiency_"+(TString)"_R="+Form("%.1f",arrayRadius[iRadius])+"_"+Datasets[iDataset]+"_centralityComp_logz");

  TString textContext(contextCustomOneField(*texDatasetsComparisonCommonDenominator, ""));

  TString centralityLegend[nCentralityBins]; CentralityLegend(centralityLegend, arrayCentralityIntervals, nCentralityBins);

  // Draw_TH2_Histograms(H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined, centralityLegend, nCentralityBins, textContext, pdfName, texPtJetRecX, texPtJetGenX, texCollisionDataInfo, "logz");

  Draw_TH1_Histograms_in_one(H1D_fakeRatio, centralityLegend, nCentralityBins, textContext, pdfName, texPtJetRecX, texFakeRatio, texCollisionDataInfo, drawnWindowAuto, legendPlacementAuto, contextPlacementAuto, "");
  // Draw_TH2_Histograms(H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined, centralityLegend, nCentralityBins, textContext, pdfName_logz, texPtJetRecX, texPtJetGenX, texCollisionDataInfo, drawnWindowAuto, legendPlacementAuto, contextPlacementAuto, "logz");
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
  TString* pdfName_logz = new TString("responseMatrix_fluctuationsBackground_"+(TString)"_R="+Form("%.1f",arrayRadius[iRadius])+"_"+Datasets[iDataset]+"_centralityComp_logz");
  TString* pdfNameFullRes_logz = new TString("responseMatrix_fluctuationsBackground_"+(TString)"_R="+Form("%.1f",arrayRadius[iRadius])+"_"+Datasets[iDataset]+"_centralityComp_FullRes_logz");
  TString* pdfName = new TString("responseMatrix_fluctuationsBackground_"+(TString)"_R="+Form("%.1f",arrayRadius[iRadius])+"_"+Datasets[iDataset]+"_centralityComp");
  TString* pdfNameFullRes = new TString("responseMatrix_fluctuationsBackground_"+(TString)"_R="+Form("%.1f",arrayRadius[iRadius])+"_"+Datasets[iDataset]+"_centralityComp_FullRes");

  TString textContext(contextCustomOneField(*texDatasetsComparisonCommonDenominator, ""));

  TString centralityLegend[nCentralityBins]; CentralityLegend(centralityLegend, arrayCentralityIntervals, nCentralityBins);

  // Draw_TH2_Histograms(H2D_jetPtResponseMatrix_fluctuations, centralityLegend, nCentralityBins, textContext, pdfName, texPtJetRecX, texPtJetGenX, texCollisionDataInfo, "logz");

  Draw_TH2_Histograms(H2D_jetPtResponseMatrix_fluctuations, centralityLegend, nCentralityBins, textContext, pdfName_logz, texPtJetRecX, texPtJetGenX, texCollisionDataInfo, drawnWindowAuto, legendPlacementAuto, contextPlacementAuto, "logz");
  Draw_TH2_Histograms(H2D_jetPtResponseMatrix_fluctuations, centralityLegend, nCentralityBins, textContext, pdfName, texPtJetRecX, texPtJetGenX, texCollisionDataInfo, drawnWindowAuto, legendPlacementAuto, contextPlacementAuto, "");
}

void Draw_ResponseMatrices_detectorResponse(int iDataset, int iRadius) {

  TH2D* H2D_jetPtResponseMatrix_detectorResponse;

  Get_PtResponseMatrix_detectorResponse(H2D_jetPtResponseMatrix_detectorResponse, iDataset, iRadius);

  TString* pdfName = new TString("responseMatrix_detectorEffects_"+jetType[iJetType]+"_"+Datasets[iDataset]);
  TString* pdfName_logz = new TString("responseMatrix_detectorEffects_"+(TString)"_R="+Form("%.1f",arrayRadius[iRadius])+"_"+Datasets[iDataset]+"_logz");

  TString textContext(contextCustomOneField(*texDatasetsComparisonCommonDenominator, ""));

  // Draw_TH2_Histograms(H2D_jetPtResponseMatrix_detectorResponse, centralityLegend, nCentralityBins, textContext, pdfName, texPtJetRecX, texPtJetGenX, texCollisionDataInfo, "logz");

  Draw_TH2_Histogram(H2D_jetPtResponseMatrix_detectorResponse, textContext, pdfName, texPtJetRecX, texPtJetGenX, texCollisionDataInfo, drawnWindowAuto, legendPlacementAuto, contextPlacementAuto, "");
  Draw_TH2_Histogram(H2D_jetPtResponseMatrix_detectorResponse, textContext, pdfName_logz, texPtJetRecX, texPtJetGenX, texCollisionDataInfo, drawnWindowAuto, legendPlacementAuto, contextPlacementAuto, "logz");
}

void Draw_ResponseMatrices_DetectorAndFluctuationsCombined(int iDataset, int iRadius, std::string options) {

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
    Get_PtResponseMatrix_DetectorAndFluctuationsCombined(H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined[iCent], H2D_jetPtResponseMatrix_detectorResponse[iCent], H2D_jetPtResponseMatrix_fluctuations[iCent], iDataset, iRadius, centRange, options);
    ReweightResponseMatrixWithPrior(H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined[iCent], iDataset, iRadius, centRange, options);
  }
  TString* pdfName = new TString("responseMatrix_combined"+(TString)"_R="+Form("%.1f",arrayRadius[iRadius])+"_"+Datasets[iDataset]+"_centralityComp");
  TString* pdfName_logz = new TString("responseMatrix_combined"+(TString)"_R="+Form("%.1f",arrayRadius[iRadius])+"_"+Datasets[iDataset]+"_centralityComp_logz");

  TString textContext(contextCustomOneField(*texDatasetsComparisonCommonDenominator, ""));

  TString centralityLegend[nCentralityBins]; CentralityLegend(centralityLegend, arrayCentralityIntervals, nCentralityBins);

  // Draw_TH2_Histograms(H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined, centralityLegend, nCentralityBins, textContext, pdfName, texPtJetRecX, texPtJetGenX, texCollisionDataInfo, "logz");

  Draw_TH2_Histograms(H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined, centralityLegend, nCentralityBins, textContext, pdfName, texPtJetRecX, texPtJetGenX, texCollisionDataInfo, drawnWindowAuto, legendPlacementAuto, contextPlacementAuto, "");
  Draw_TH2_Histograms(H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined, centralityLegend, nCentralityBins, textContext, pdfName_logz, texPtJetRecX, texPtJetGenX, texCollisionDataInfo, drawnWindowAuto, legendPlacementAuto, contextPlacementAuto, "logz");
}

void Draw_Pt_spectrum_unfolded(int iDataset, int iRadius, int unfoldParameterBayes, std::string options) {

  TH1D* H1D_jetPt_measured[nCentralityBins];
  TH1D* H1D_jetPt_measured_genBinning[nCentralityBins];
  TH1D* H1D_jetPt_unfolded[nCentralityBins];
  TH1D* H1D_jetPt_unfoldedThenRefolded[nCentralityBins];
  TH1D* H1D_jetPt_unfoldedThenRefolded_RooUnfoldMethod[nCentralityBins];
  TH1D* H1D_jetPt_unfolded_mcpComp[nCentralityBins][2];
  TH1D* H1D_jetPt_unfolded_run2Comp[nCentralityBins][2];
  TH1D* H1D_jetPt_unfolded_measuredComp[nCentralityBins][2];
  TH1D* H1D_jetPt_unfolded_refoldedComp[nCentralityBins][3];
  TH1D* H1D_jetPt_mcp;
  TH1D* H1D_jetPt_run2;
  TH1D* H1D_jetPt_ratio_mcp[nCentralityBins];
  TH1D* H1D_jetPt_ratio_run2[nCentralityBins];
  TH1D* H1D_jetPt_ratio_measured[nCentralityBins];
  TH1D* H1D_jetPt_ratio_measuredRefolded[nCentralityBins][2];

  if (isPbPb) {
    H1D_jetPt_run2 = (TH1D*)((TH1D*)file_O2Analysis_run2ComparisonFile->Get("Bayesian_Unfoldediter15"))->Clone("H1D_jetPt_run2");
    int NcollRun2 = 4619963; // central (see Laura discussion mattermost) 
    H1D_jetPt_run2->Scale(1./NcollRun2);
  }

  bool divideSuccessMcp[nCentralityBins];
  bool divideSuccessRun2[nCentralityBins];
  bool divideSuccessMeasured[nCentralityBins];
  bool divideSuccessMeasuredRefolded[nCentralityBins][2];
  float centRange[2];
  TString partialUniqueSpecifier;

  if (!useFineBinningTest) {
    Get_Pt_spectrum_mcp(H1D_jetPt_mcp, iDataset, iRadius, controlMC, options);
  } else {
    Get_Pt_spectrum_mcp_fineBinning(H1D_jetPt_mcp, iDataset, iRadius, controlMC, options);
  }
  // H1D_jetPt_mcp = (TH1D*)H1D_jetPt_mcp_defaultBin->Rebin(nBinPtJetsGen[iRadius],"jetPt_mcp_rebinned_"+RadiusLegend[iRadius]+Datasets[iDataset], ptBinsJetsGen[iRadius]);

  int unfoldParameter[nCentralityBins];
  for(int iCent = 0; iCent < nCentralityBins; iCent++){
    centRange[0] = arrayCentralityIntervals[iCent][0];
    centRange[1] = arrayCentralityIntervals[iCent][1];
    partialUniqueSpecifier = Datasets[iDataset]+"_R="+Form("%.1f",arrayRadius[iRadius])+"_@cent["+Form("%.1f", centRange[0])+","+Form("%.1f", centRange[1])+"]";
    
    unfoldParameter[iCent] = Get_Pt_spectrum_unfolded(H1D_jetPt_unfolded[iCent], iDataset, iRadius, centRange, unfoldParameterBayes, options).first;
    
    // comparison with raw measured
    if (!useFineBinningTest) {
      Get_Pt_spectrum_bkgCorrected_genBinning(H1D_jetPt_measured_genBinning[iCent], iDataset, iRadius, centRange, options);
    } else {
      Get_Pt_spectrum_bkgCorrected_fineBinning(H1D_jetPt_measured_genBinning[iCent], iDataset, iRadius, centRange, options);
    }

    H1D_jetPt_unfolded_measuredComp[iCent][0] = (TH1D*)H1D_jetPt_unfolded[iCent]->Clone("H1D_jetPt_unfolded_measuredComp"+partialUniqueSpecifier);
    H1D_jetPt_unfolded_measuredComp[iCent][1] = (TH1D*)H1D_jetPt_measured_genBinning[iCent]->Clone("H1D_jetPt_measured_genBinning_measuredComp"+partialUniqueSpecifier);
    H1D_jetPt_ratio_measured[iCent] = (TH1D*)H1D_jetPt_measured_genBinning[iCent]->Clone("H1D_jetPt_ratio_measured"+partialUniqueSpecifier);
    divideSuccessMeasured[iCent] = H1D_jetPt_ratio_measured[iCent]->Divide(H1D_jetPt_unfolded[iCent]);

    // comparison with mcp truth
    H1D_jetPt_unfolded_mcpComp[iCent][0] = (TH1D*)H1D_jetPt_unfolded[iCent]->Clone("H1D_jetPt_unfolded_mcpComp"+partialUniqueSpecifier);
    H1D_jetPt_unfolded_mcpComp[iCent][1] = (TH1D*)H1D_jetPt_mcp->Clone("H1D_jetPt_unfolded_mcpComp"+partialUniqueSpecifier);
    H1D_jetPt_ratio_mcp[iCent] = (TH1D*)H1D_jetPt_mcp->Clone("H1D_jetPt_ratio_mcp"+partialUniqueSpecifier);
    divideSuccessMcp[iCent] = H1D_jetPt_ratio_mcp[iCent]->Divide(H1D_jetPt_unfolded[iCent]);

    // comparison with run2 
    if (isPbPb) {
      H1D_jetPt_unfolded_run2Comp[iCent][0] = (TH1D*)H1D_jetPt_unfolded[iCent]->Clone("H1D_jetPt_unfolded_run2Comp"+partialUniqueSpecifier);
      H1D_jetPt_unfolded_run2Comp[iCent][1] = (TH1D*)H1D_jetPt_run2->Clone("H1D_jetPt_unfolded_run2Comp"+partialUniqueSpecifier);
      H1D_jetPt_ratio_run2[iCent] = (TH1D*)H1D_jetPt_run2->Clone("H1D_jetPt_ratio_run2"+partialUniqueSpecifier);
      divideSuccessRun2[iCent] = H1D_jetPt_ratio_run2[iCent]->Divide(H1D_jetPt_unfolded[iCent]);
    }

    // comparison with refolded
    if (!useFineBinningTest) {
      Get_Pt_spectrum_bkgCorrected(H1D_jetPt_measured[iCent], iDataset, iRadius, centRange, options);
    } else {
      Get_Pt_spectrum_bkgCorrected_fineBinning(H1D_jetPt_measured[iCent], iDataset, iRadius, centRange, options);
    }
    Get_Pt_spectrum_unfoldedThenRefolded(H1D_jetPt_unfoldedThenRefolded[iCent], iDataset, iRadius, centRange, unfoldParameterBayes, options);
    Get_Pt_spectrum_unfoldedThenRefolded_RooUnfoldMethod(H1D_jetPt_unfoldedThenRefolded_RooUnfoldMethod[iCent], iDataset, iRadius, centRange, unfoldParameterBayes, options);
    H1D_jetPt_unfolded_refoldedComp[iCent][0] = (TH1D*)H1D_jetPt_unfoldedThenRefolded[iCent]->Clone("H1D_jetPt_refolded_refoldedComp"+partialUniqueSpecifier);
    H1D_jetPt_unfolded_refoldedComp[iCent][1] = (TH1D*)H1D_jetPt_unfoldedThenRefolded_RooUnfoldMethod[iCent]->Clone("H1D_jetPt_unfoldedThenRefolded_RooUnfoldMethod"+partialUniqueSpecifier);
    H1D_jetPt_unfolded_refoldedComp[iCent][2] = (TH1D*)H1D_jetPt_measured[iCent]->Clone("H1D_jetPt_measured_refoldedComp"+partialUniqueSpecifier);
    H1D_jetPt_ratio_measuredRefolded[iCent][0] = (TH1D*)H1D_jetPt_unfoldedThenRefolded[iCent]->Clone("H1D_jetPt_ratio_refoldedComp"+partialUniqueSpecifier);
    H1D_jetPt_ratio_measuredRefolded[iCent][1] = (TH1D*)H1D_jetPt_unfoldedThenRefolded_RooUnfoldMethod[iCent]->Clone("H1D_jetPt_ratio_refoldedComp_RooUnfoldMethod"+partialUniqueSpecifier);
    divideSuccessMeasuredRefolded[iCent][0] = H1D_jetPt_ratio_measuredRefolded[iCent][0]->Divide(H1D_jetPt_measured[iCent]);
    divideSuccessMeasuredRefolded[iCent][1] = H1D_jetPt_ratio_measuredRefolded[iCent][1]->Divide(H1D_jetPt_measured[iCent]);


    // H1D_jetPt_unfolded_refoldedComp_RooUnfoldMethod[iCent][0] = (TH1D*)H1D_jetPt_unfoldedThenRefolded_RooUnfoldMethod[iCent]->Clone("H1D_jetPt_refolded_refoldedComp_RooUnfoldMethod"+partialUniqueSpecifier);
    // H1D_jetPt_unfolded_refoldedComp_RooUnfoldMethod[iCent][1] = (TH1D*)H1D_jetPt_measured_RooUnfoldMethod[iCent]->Clone("H1D_jetPt_measured_refoldedComp_RooUnfoldMethod"+partialUniqueSpecifier);
    // H1D_jetPt_ratio_measuredRefolded_RooUnfoldMethod[iCent] = (TH1D*)H1D_jetPt_unfoldedThenRefolded_RooUnfoldMethod[iCent]->Clone("H1D_jetPt_ratio_refoldedComp_RooUnfoldMethod"+partialUniqueSpecifier);
    // divideSuccessMeasuredRefolded_RooUnfoldMethod[iCent] = H1D_jetPt_ratio_measuredRefolded_RooUnfoldMethod[iCent]->Divide(H1D_jetPt_measured[iCent]);
  }

  TString unfoldingCode;
  if (useInitialResponseMethod){
    unfoldingCode = "myUnfold";
  } else {
    unfoldingCode = "joUnfold";
  }
  if (controlMC){
    unfoldingCode += "_controlMC";
  }
  TString unfoldingInfo[nCentralityBins];
  for(int iCent = 0; iCent < nCentralityBins; iCent++){
    unfoldingInfo[iCent] = (TString)unfoldingMethod+"-k="+Form("%i", unfoldParameter[iCent])+"-"+(TString)mergingPrior+"-"+(TString)unfoldingPrior+"-"+unfoldingCode;
  }

  struct stat st1{};
  if (stat("pdfFolder/IterationsDump", &st1) == -1) {
      mkdir("pdfFolder/IterationsDump", 0700);
  }
  struct stat st2{};
  if (stat("pngFolder/IterationsDump", &st2) == -1) {
      mkdir("pngFolder/IterationsDump", 0700);
  }

  TString textContext[nCentralityBins];
  for(int iCent = 0; iCent < nCentralityBins; iCent++){    
    centRange[0] = arrayCentralityIntervals[iCent][0];
    centRange[1] = arrayCentralityIntervals[iCent][1];
    textContext[iCent] = contextCustomThreeFields(*texDatasetsComparisonCommonDenominator, contextJetRadius(arrayRadius[iRadius]), contextCentRange(centRange), "");
  }

  TString dummyLegend[1] = {""};

  TString* yAxisLabel = texCount;
  if (options.find("evtNorm") != std::string::npos && doEvtNorm) {
    yAxisLabel = texJetPtYield_EventNorm;
  }
  if (options.find("entriesNorm") != std::string::npos) {
    yAxisLabel = texJetPtYield_EntriesNorm;
  }

    // comparison with raw measured
  TString unfoldedMeasuredCompLegend[2] = {"unfolded data", "measured raw (gen binning)"};
  for(int iCent = 0; iCent < nCentralityBins; iCent++){
    TString* pdfName_measuredComp = new TString("IterationsDump/jet_"+Datasets[iDataset]+"_R="+Form("%.1f", arrayRadius[iRadius])+"_@cent["+Form("%.0f", centRange[0])+","+Form("%.0f", centRange[1])+"]"+"_Pt_unfolded_"+unfoldingInfo[iCent]+"_measuredComp");
    Draw_TH1_Histograms_in_one(H1D_jetPt_unfolded_measuredComp[iCent], unfoldedMeasuredCompLegend, 2, textContext[iCent], pdfName_measuredComp, texPtX, yAxisLabel, texCollisionDataInfo, drawnWindowAuto, legendPlacementAuto, contextPlacementAuto, "logy");
    if (divideSuccessMeasured[iCent]){
      TString* pdfName_ratio_measured = new TString("IterationsDump/jet_"+Datasets[iDataset]+"_R="+Form("%.1f", arrayRadius[iRadius])+"_@cent["+Form("%.0f", centRange[0])+","+Form("%.0f", centRange[1])+"]"+"_Pt_unfolded_"+unfoldingInfo[iCent]+"_ratioMeasured");
      Draw_TH1_Histogram(H1D_jetPt_ratio_measured[iCent], textContext[iCent], pdfName_ratio_measured, texPtX, texRatioMeasuredUnfolded, texCollisionDataInfo, drawnWindowAuto, legendPlacementAuto, contextPlacementAuto, "standardratio,ratioLine");
    }
  }

    // comparison with mcp truth
  TString unfoldedTruthCompLegend[2] = {"unfolded data", "mcp truth"};
  for(int iCent = 0; iCent < nCentralityBins; iCent++){
    TString* pdfName_mcpComp = new TString("IterationsDump/jet_"+Datasets[iDataset]+"_R="+Form("%.1f", arrayRadius[iRadius])+"_@cent["+Form("%.0f", centRange[0])+","+Form("%.0f", centRange[1])+"]"+"_Pt_unfolded_"+unfoldingInfo[iCent]+"_mcpComp");
    Draw_TH1_Histograms_in_one(H1D_jetPt_unfolded_mcpComp[iCent], unfoldedTruthCompLegend, 2, textContext[iCent], pdfName_mcpComp, texPtX, yAxisLabel, texCollisionDataInfo, drawnWindowAuto, legendPlacementAuto, contextPlacementAuto, "logy");
    if (divideSuccessMcp[iCent]){
      TString* pdfName_ratio_mcp = new TString("IterationsDump/jet_"+Datasets[iDataset]+"_R="+Form("%.1f", arrayRadius[iRadius])+"_@cent["+Form("%.0f", centRange[0])+","+Form("%.0f", centRange[1])+"]"+"_Pt_unfolded_"+unfoldingInfo[iCent]+"_ratioMcp");
      Draw_TH1_Histogram(H1D_jetPt_ratio_mcp[iCent], textContext[iCent], pdfName_ratio_mcp, texPtX, texRatioMcpUnfolded, texCollisionDataInfo, drawnWindowAuto, legendPlacementAuto, contextPlacementAuto, "standardratio,ratioLine");
    }
  }


  // comparison with refolded
  TString unfoldedRefoldedCompLegend[3] = {"refolded manually", "refolded roounfold", "measured"};
  for(int iCent = 0; iCent < nCentralityBins; iCent++){
    TString* pdfName_refoldedComp = new TString("IterationsDump/jet_"+Datasets[iDataset]+"_R="+Form("%.1f", arrayRadius[iRadius])+"_@cent["+Form("%.0f", centRange[0])+","+Form("%.0f", centRange[1])+"]"+"_Pt_unfolded_"+unfoldingInfo[iCent]+"_RefoldedComp");
    Draw_TH1_Histograms_in_one(H1D_jetPt_unfolded_refoldedComp[iCent], unfoldedRefoldedCompLegend, 3, textContext[iCent], pdfName_refoldedComp, texPtX, yAxisLabel, texCollisionDataInfo, drawnWindowAuto, legendPlacementAuto, contextPlacementAuto, "logy");
    if (divideSuccessMeasuredRefolded[iCent]) {
      TString* pdfName_ratio_refoldedComp = new TString("IterationsDump/jet_"+Datasets[iDataset]+"_R="+Form("%.1f", arrayRadius[iRadius])+"_@cent["+Form("%.0f", centRange[0])+","+Form("%.0f", centRange[1])+"]"+"_Pt_unfolded_"+unfoldingInfo[iCent]+"_ratioRefoldedUnfolded");
      Draw_TH1_Histograms_in_one(H1D_jetPt_ratio_measuredRefolded[iCent], unfoldedRefoldedCompLegend, 2, textContext[iCent], pdfName_ratio_refoldedComp, texPtX, texRatioRefoldedMeasured, texCollisionDataInfo, drawnWindowAuto, legendPlacementAuto, contextPlacementAuto, "standardratio,ratioLine");
    }
  }


  // comparison with Run 2
  if (isPbPb) {
    TString unfoldedRun2CompLegend[2] = {"unfolded Run3", "unfolded Run2 0-10%"};
    for(int iCent = 0; iCent < nCentralityBins; iCent++){
      TString* pdfName_run2Comp = new TString("IterationsDump/jet_"+Datasets[iDataset]+"_R="+Form("%.1f", arrayRadius[iRadius])+"_@cent["+Form("%.0f", centRange[0])+","+Form("%.0f", centRange[1])+"]"+"_Pt_unfolded_"+unfoldingInfo[iCent]+"_run2Comp");
      Draw_TH1_Histograms_in_one(H1D_jetPt_unfolded_run2Comp[iCent], unfoldedRun2CompLegend, 2, textContext[iCent], pdfName_run2Comp, texPtX, yAxisLabel, texCollisionDataInfo, drawnWindowAuto, legendPlacementAuto, contextPlacementAuto, "logy");
      if (divideSuccessRun2[iCent]){
        TString* pdfName_ratio_run2 = new TString("IterationsDump/jet_"+Datasets[iDataset]+"_R="+Form("%.1f", arrayRadius[iRadius])+"_@cent["+Form("%.0f", centRange[0])+","+Form("%.0f", centRange[1])+"]"+"_Pt_unfolded_"+unfoldingInfo[iCent]+"_ratioRun2");
        Draw_TH1_Histogram(H1D_jetPt_ratio_run2[iCent], textContext[iCent], pdfName_ratio_run2, texPtX, texRatioRun2Unfolded, texCollisionDataInfo, drawnWindowAuto, legendPlacementAuto, contextPlacementAuto, "standardratio,ratioLine");
      }
    }
  }
}

void Draw_Pt_spectrum_unfolded_parameterVariation(int iDataset, int iRadius, int unfoldIterationMin, int unfoldIterationMax, int step, std::string options) {

  const int nUnfoldIteration = std::floor((unfoldIterationMax - unfoldIterationMin + 1)/step);

  TH1D* H1D_jetPt_measured[nCentralityBins];
  TH1D* H1D_jetPt_measured_genBinning[nCentralityBins];
  TH1D* H1D_jetPt_unfolded[nCentralityBins][nUnfoldIteration];
  TH1D* H1D_jetPt_unfoldedThenRefolded[nCentralityBins][nUnfoldIteration];
  TH1D* H1D_jetPt_unfolded_mcpComp[nCentralityBins][nUnfoldIteration+1];
  TH1D* H1D_jetPt_unfolded_measuredComp[nCentralityBins][nUnfoldIteration+1];
  TH1D* H1D_jetPt_unfolded_refoldedComp[nCentralityBins][nUnfoldIteration+1];
  TH1D* H1D_jetPt_mcp;
  TH1D* H1D_jetPt_ratio_mcp[nCentralityBins][nUnfoldIteration];
  TH1D* H1D_jetPt_ratio_measured[nCentralityBins][nUnfoldIteration];
  TH1D* H1D_jetPt_ratio_measuredRefolded[nCentralityBins][nUnfoldIteration];

  bool divideSuccessMcp[nCentralityBins][nUnfoldIteration];
  bool divideSuccessMeasured[nCentralityBins][nUnfoldIteration];
  bool divideSuccessMeasuredRefolded[nCentralityBins][nUnfoldIteration];
  float centRange[2];
  TString partialUniqueSpecifier;

  Get_Pt_spectrum_mcp(H1D_jetPt_mcp, iDataset, iRadius, controlMC, options);
  // H1D_jetPt_mcp = (TH1D*)H1D_jetPt_mcp_defaultBin->Rebin(nBinPtJetsGen[iRadius],"jetPt_mcp_rebinned_"+RadiusLegend[iRadius]+Datasets[iDataset], ptBinsJetsGen[iRadius]);

  for(int iCent = 0; iCent < nCentralityBins; iCent++){
    centRange[0] = arrayCentralityIntervals[iCent][0];
    centRange[1] = arrayCentralityIntervals[iCent][1];
    partialUniqueSpecifier = Datasets[iDataset]+"_R="+Form("%.1f",arrayRadius[iRadius])+"_@cent["+Form("%.1f", centRange[0])+","+Form("%.1f", centRange[1])+"]";
    Get_Pt_spectrum_bkgCorrected(H1D_jetPt_measured[iCent], iDataset, iRadius, centRange, options);
    Get_Pt_spectrum_bkgCorrected_genBinning(H1D_jetPt_measured_genBinning[iCent], iDataset, iRadius, centRange, options);

    int unfoldIterations = 0;

    for(int iUnfoldIteration = 0; iUnfoldIteration < nUnfoldIteration; iUnfoldIteration++){

      unfoldIterations = iUnfoldIteration * step + unfoldIterationMin; 
      // unfoldingIterationLegend[iUnfoldIteration] = unfoldIterations;

      cout << "((((((((((((()))))))))))))" << endl;
      cout << "Iteration "<< iUnfoldIteration << endl;
      cout << "((((((((((((()))))))))))))" << endl;
      Get_Pt_spectrum_unfolded(H1D_jetPt_unfolded[iCent][iUnfoldIteration], iDataset, iRadius, centRange, unfoldIterations, options);

      // comparison with raw measured
      H1D_jetPt_unfolded_measuredComp[iCent][iUnfoldIteration] = (TH1D*)H1D_jetPt_unfolded[iCent][iUnfoldIteration]->Clone("H1D_jetPt_unfolded_measuredComp"+partialUniqueSpecifier);
      H1D_jetPt_ratio_measured[iCent][iUnfoldIteration] = (TH1D*)H1D_jetPt_measured_genBinning[iCent]->Clone("H1D_jetPt_ratio_measured"+partialUniqueSpecifier);
      divideSuccessMeasured[iCent][iUnfoldIteration] = H1D_jetPt_ratio_measured[iCent][iUnfoldIteration]->Divide(H1D_jetPt_unfolded[iCent][iUnfoldIteration]);

      // comparison with mcp truth
      H1D_jetPt_unfolded_mcpComp[iCent][iUnfoldIteration] = (TH1D*)H1D_jetPt_unfolded[iCent][iUnfoldIteration]->Clone("H1D_jetPt_unfolded_mcpComp"+partialUniqueSpecifier);
      H1D_jetPt_ratio_mcp[iCent][iUnfoldIteration] = (TH1D*)H1D_jetPt_mcp->Clone("H1D_jetPt_ratio_mcp"+partialUniqueSpecifier);
      divideSuccessMcp[iCent][iUnfoldIteration] = H1D_jetPt_ratio_mcp[iCent][iUnfoldIteration]->Divide(H1D_jetPt_unfolded[iCent][iUnfoldIteration]);


      // comparison with refolded
      Get_Pt_spectrum_unfoldedThenRefolded(H1D_jetPt_unfoldedThenRefolded[iCent][iUnfoldIteration], iDataset, iRadius, centRange, unfoldIterations, options);
      H1D_jetPt_unfolded_refoldedComp[iCent][iUnfoldIteration] = (TH1D*)H1D_jetPt_unfoldedThenRefolded[iCent][iUnfoldIteration]->Clone("H1D_jetPt_refolded_refoldedComp"+partialUniqueSpecifier);
      H1D_jetPt_ratio_measuredRefolded[iCent][iUnfoldIteration] = (TH1D*)H1D_jetPt_unfoldedThenRefolded[iCent][iUnfoldIteration]->Clone("H1D_jetPt_ratio_refoldedComp"+partialUniqueSpecifier);
      divideSuccessMeasuredRefolded[iCent][iUnfoldIteration] = H1D_jetPt_ratio_measuredRefolded[iCent][iUnfoldIteration]->Divide(H1D_jetPt_measured[iCent]);
    }
    H1D_jetPt_unfolded_measuredComp[iCent][nUnfoldIteration] = (TH1D*)H1D_jetPt_measured_genBinning[iCent]->Clone("H1D_jetPt_measured_genBinning_measuredComp"+partialUniqueSpecifier);
    H1D_jetPt_unfolded_mcpComp[iCent][nUnfoldIteration] = (TH1D*)H1D_jetPt_mcp->Clone("H1D_jetPt_unfolded_mcpComp"+partialUniqueSpecifier);
    H1D_jetPt_unfolded_refoldedComp[iCent][nUnfoldIteration] = (TH1D*)H1D_jetPt_measured[iCent]->Clone("H1D_jetPt_measured_refoldedComp"+partialUniqueSpecifier);
  }

  TString unfoldingInfo = (TString)unfoldingMethod+"_"+(TString)unfoldingPrior+"_kmax="+Form("%i", unfoldIterationMax);

  TString unfoldingIterationLegend[nUnfoldIteration+1]; IterationLegend(unfoldingIterationLegend, unfoldIterationMin, unfoldIterationMax, step);

  

  for(int iCent = 0; iCent < nCentralityBins; iCent++){
    partialUniqueSpecifier = Datasets[iDataset]+"_R="+Form("%.1f",arrayRadius[iRadius])+"_@cent["+Form("%.1f", centRange[0])+","+Form("%.1f", centRange[1])+"]";

    TString textContext(contextCustomTwoFields(*texDatasetsComparisonCommonDenominator, (TString)Form("%.1f", centRange[0])+"-"+Form("%.1f", centRange[1])+"%", ""));

    TString* pdfName = new TString("jet_"+jetType[iJetType]+"_"+jetLevel[iJetLevel]+"_"+partialUniqueSpecifier+"_Pt_unfolded_"+unfoldingInfo);

    TString* yAxisLabel = texCount;
    if (options.find("evtNorm") != std::string::npos && doEvtNorm) {
      yAxisLabel = texJetPtYield_EventNorm;
    }
    if (options.find("entriesNorm") != std::string::npos) {
      yAxisLabel = texJetPtYield_EntriesNorm;
    }
    Draw_TH1_Histograms_in_one(H1D_jetPt_unfolded[iCent], unfoldingIterationLegend, nUnfoldIteration, textContext, pdfName, texPtX, yAxisLabel, texCollisionDataInfo, drawnWindowAuto, legendPlacementAuto, contextPlacementAuto, "logy");


      // comparison with raw measured
    // TString unfoldedMeasuredCompLegend[2] = {"unfolded data", "measured raw (gen binning)"};
    unfoldingIterationLegend[nUnfoldIteration] = (TString)"raw measured";
    TString* pdfName_measuredComp = new TString("jet_"+jetType[iJetType]+"_"+jetLevel[iJetLevel]+"_"+partialUniqueSpecifier+"_Pt_unfolded_"+unfoldingInfo+"_measuredComp");
    Draw_TH1_Histograms_in_one(H1D_jetPt_unfolded_measuredComp[iCent], unfoldingIterationLegend, nUnfoldIteration+1, textContext, pdfName_measuredComp, texPtX, yAxisLabel, texCollisionDataInfo, drawnWindowAuto, legendPlacementAuto, contextPlacementAuto, "logy");

    bool divideSuccessMeasured_boolsum = true;
    for(int iUnfoldIteration = 0; iUnfoldIteration < nUnfoldIteration; iUnfoldIteration++){
      if (!(divideSuccessMeasured[iUnfoldIteration])) {
        divideSuccessMeasured_boolsum = false;
      }
    }
    if (divideSuccessMeasured_boolsum){
      TString* pdfName_ratio_measured = new TString("jet_"+jetType[iJetType]+"_"+jetLevel[iJetLevel]+"_"+partialUniqueSpecifier+"_Pt_unfolded_"+unfoldingInfo+"_ratioMeasured");
      Draw_TH1_Histograms_in_one(H1D_jetPt_ratio_measured[iCent], unfoldingIterationLegend, nUnfoldIteration, textContext, pdfName_ratio_measured, texPtX, texRatioMeasuredUnfolded, texCollisionDataInfo, drawnWindowAuto, legendPlacementAuto, contextPlacementAuto, "standardratio,ratioLine");
    }


      // comparison with mcp truth
    // TString unfoldedTruthCompLegend[2] = {"unfolded data", "mcp truth"};
    unfoldingIterationLegend[nUnfoldIteration] = (TString)"mcp";
    TString* pdfName_mcpComp = new TString("jet_"+jetType[iJetType]+"_"+jetLevel[iJetLevel]+"_"+partialUniqueSpecifier+"_@cent["+Form("%.1f", centRange[0])+","+Form("%.1f", centRange[1])+"]"+"_Pt_unfolded_"+unfoldingInfo+"_mcpComp");
    Draw_TH1_Histograms_in_one(H1D_jetPt_unfolded_mcpComp[iCent], unfoldingIterationLegend, nUnfoldIteration+1, textContext, pdfName_mcpComp, texPtX, yAxisLabel, texCollisionDataInfo, drawnWindowAuto, legendPlacementAuto, contextPlacementAuto, "logy");

    bool divideSuccessMcp_boolsum = true;
    for(int iUnfoldIteration = 0; iUnfoldIteration < nUnfoldIteration; iUnfoldIteration++){
      if (!(divideSuccessMcp[iUnfoldIteration])) {
        divideSuccessMcp_boolsum = false;
      }
    }
    if (divideSuccessMcp_boolsum){
      TString* pdfName_ratio_mcp = new TString("jet_"+jetType[iJetType]+"_"+jetLevel[iJetLevel]+"_"+partialUniqueSpecifier+"_Pt_unfolded_"+unfoldingInfo+"_ratioMcp");
      Draw_TH1_Histograms_in_one(H1D_jetPt_ratio_mcp[iCent], unfoldingIterationLegend, nUnfoldIteration, textContext, pdfName_ratio_mcp, texPtX, texRatioMcpUnfolded, texCollisionDataInfo, drawnWindowAuto, legendPlacementAuto, contextPlacementAuto, "standardratio,ratioLine");
    }

    // comparison with refolded
    // TString unfoldedRefoldedCompLegend[2] = {"refolded", "measured"};
    unfoldingIterationLegend[nUnfoldIteration] = (TString)"measured";
    for(int iIteration = 0; iIteration < nUnfoldIteration; iIteration++){
      unfoldingIterationLegend[iIteration] += (TString)" refolded";
    }
    TString* pdfName_refoldedComp = new TString("jet_"+jetType[iJetType]+"_"+jetLevel[iJetLevel]+"_"+partialUniqueSpecifier+"_@cent["+Form("%.1f", centRange[0])+","+Form("%.1f", centRange[1])+"]"+"_Pt_unfolded_"+unfoldingInfo+"_RefoldedComp");
    Draw_TH1_Histograms_in_one(H1D_jetPt_unfolded_refoldedComp[iCent], unfoldingIterationLegend, nUnfoldIteration+1, textContext, pdfName_refoldedComp, texPtX, yAxisLabel, texCollisionDataInfo, drawnWindowAuto, legendPlacementAuto, contextPlacementAuto, "logy");

    bool divideSuccessRefoldedComp_boolsum = true;
    for(int iUnfoldIteration = 0; iUnfoldIteration < nUnfoldIteration; iUnfoldIteration++){
      if (!(divideSuccessMeasuredRefolded[iUnfoldIteration])) {
        divideSuccessRefoldedComp_boolsum = false;
      }
    }
    if (divideSuccessRefoldedComp_boolsum){
      TString* pdfName_ratio_refoldedComp = new TString("jet_"+jetType[iJetType]+"_"+jetLevel[iJetLevel]+"_"+partialUniqueSpecifier+"_Pt_unfolded_"+unfoldingInfo+"_ratioRefoldedUnfolded");
      Draw_TH1_Histograms_in_one(H1D_jetPt_ratio_measuredRefolded[iCent], unfoldingIterationLegend, nUnfoldIteration, textContext, pdfName_ratio_refoldedComp, texPtX, texRatioRefoldedMeasured, texCollisionDataInfo, drawnWindowAuto, legendPlacementAuto, contextPlacementAuto, "standardratio,ratioLine");
    }
  }
}




// rename refoldedUnfolded as closure test?
// and try and spend 15 min to clean hist names for the spectrum analysis



// WARNING FOR EFFICIENCIES I SHOULD REREAD THIS BELOW!!
// hMcEfficiency_vsPt->Divide(hMcSignalCount_vsPt,TrueV0PtSpectrum_AnalysisBins, 1., 1., "b"); // option b for binomial because efficiency: https://twiki.cern.ch/twiki/bin/view/ALICE/PWGLFPAGSTRANGENESSEfficiency
