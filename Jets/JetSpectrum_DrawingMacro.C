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
#include <RooUnfold.h> // one should likely do `aliBuild build RooUnfold` then `alienv enter RooUnfold/latest` as alidist roounfold version is usually quite old
#include "RooUnfoldResponse.h"
#include "RooUnfoldBayes.h"
#include "RooUnfoldBinByBin.h"
#include "RooUnfoldSvd.h"
#include "TSVDUnfold.h"

//My Libraries
#include "./JetSpectrum_settings.h"
#include "./JetSpectrum_inputs.h"

#include "./JetSpectrum_ResponseMatrixFunctions.h"
#include "./JetSpectrum_ResponseMatrixFunctions.C"
#include "./JetSpectrum_SpectraGetters.h"
#include "./JetSpectrum_SpectraGetters.C"
#include "./JetSpectrum_Unfolding.h"
#include "./JetSpectrum_Unfolding.C"
#include "./JetSpectrum_EfficiencyPurityGetters.h"
#include "./JetSpectrum_EfficiencyPurityGetters.C"

#include "../Settings/AxisTitles.h"
#include "../Settings/GlobalSettings.h"
#include "../Utilities/AnalysisUtilities.h"
#include "../Utilities/HistogramUtilities.h"
#include "../Utilities/HistogramPlotting.h"
#include "../Utilities/AnalysisUtilities.C" 
#include "../Utilities/HistogramUtilities.C"
#include "../Utilities/HistogramPlotting.C" 

#include<array>
#include <iomanip>
#include <sstream>
#include <string.h>
#include <stdlib.h>     /* abort, NULL */
using namespace std;

// Misc utilities
void SetStyle(Bool_t graypalette=kFALSE);
void LoadLibs();

void Draw_ResponseMatrices_Fluctuations(int iDataset, int iRadius);
void Draw_ResponseMatrices_detectorResponse(int iDataset, int iRadius);
void Draw_ResponseMatrices_DetectorAndFluctuationsCombined(int iDataset, int iRadius, std::string options);

void Draw_Pt_spectrum_raw(int iDataset, int iRadius, std::string options);
void Draw_Pt_spectrum_mcp(int iDataset, int iRadius, std::string options);
void Draw_Pt_spectrum_mcdMatched(int iDataset, int iRadius, std::string options);
void Draw_Pt_spectrum_unfolded_FluctResponseOnly(int iDataset, int iRadius, std::string options);
void Draw_Pt_spectrum_unfolded(int iDataset, int iRadius, int unfoldParameterInput, std::string options);
void Draw_Pt_TestSpectrum_unfolded(int iDataset, int iRadius, std::string options);
void Draw_Pt_spectrum_unfolded_parameterVariation(int iDataset, int iRadius, int unfoldIterationMin, int unfoldIterationMax, int step, std::string options);

void Draw_Pt_efficiency_jets(int iDataset, int iRadius, std::string options);
void Draw_kinematicEfficiency(int iDataset, int iRadius, std::string options);
void Draw_FakeRatio(int iDataset, int iRadius, std::string options);

/////////////////////////////////////////////////////
///////////////////// Main Macro ////////////////////
/////////////////////////////////////////////////////

void JetSpectrum_DrawingMacro() {
  // Load necessary libraries
  LoadLibs();
  // Set the default style
  SetStyle();

  // TString* SaveAs_Title = new TString("");
  TString* texXtitle = new TString("");
  TString* texYtitle = new TString("");
  // TString* Extra = new TString("");

  // gathers the analysis options in a single char[]
  snprintf(optionsAnalysis, sizeof(optionsAnalysis), "%s,%s,%s", mergingPrior, unfoldingPrior, unfoldingMethod);
  cout << "Analysis options are: " << optionsAnalysis << endl;

  int iDataset = 0;
  int iRadius = 0;

  // // find a way to input mcpPrior/mcdPrior and bayes/svd as a variables rather than typed out like this
  // Draw_ResponseMatrices_Fluctuations(iDataset, iRadius);
  // Draw_ResponseMatrices_detectorResponse(iDataset, iRadius);
  Draw_ResponseMatrices_DetectorAndFluctuationsCombined(iDataset, iRadius, optionsAnalysis);

  // Draw_Pt_spectrum_unfolded_FluctResponseOnly(iDataset, iRadius, optionsAnalysis); // NOT FIXED YET - result meaningless
  // Draw_Pt_spectrum_raw(iDataset, iRadius, optionsAnalysis);
  // Draw_Pt_spectrum_mcp(iDataset, iRadius, optionsAnalysis);
  // Draw_Pt_spectrum_mcdMatched(iDataset, iRadius, optionsAnalysis);

  // Draw_Pt_efficiency_jets(iDataset, iRadius, optionsAnalysis);
  // Draw_kinematicEfficiency(iDataset, iRadius, optionsAnalysis);
  // Draw_FakeRatio(iDataset, iRadius, optionsAnalysis);

  // int unfoldParameterInput = 6;
  // Draw_Pt_spectrum_unfolded(iDataset, iRadius, unfoldParameterInput, optionsAnalysis);
  // int unfoldParameterInput2 = 8;
  // Draw_Pt_spectrum_unfolded(iDataset, iRadius, unfoldParameterInput2, optionsAnalysis);
  // int unfoldParameterInput3 = 10;
  // Draw_Pt_spectrum_unfolded(iDataset, iRadius, unfoldParameterInput3, optionsAnalysis);

  // int unfoldParameterInputMin = 4;
  // int unfoldParameterInputMax = 13;
  // int unfoldParameterInputStep = 2;
  // Draw_Pt_spectrum_unfolded_parameterVariation(iDataset, iRadius, unfoldParameterInputMin, unfoldParameterInputMax, unfoldParameterInputStep, optionsAnalysis);
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
////////////////////////////////////////////////////////////////////////// Spectrum plotting functions ///////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


void Draw_Pt_spectrum_raw(int iDataset, int iRadius, std::string options) {

  TH1D* H1D_jetPt_raw;

  Get_Pt_spectrum_bkgCorrected_recBinning(H1D_jetPt_raw, iDataset, iRadius, options);

  TString* pdfName = new TString("jet_"+jetType[iJetType]+"_"+jetLevel[iJetLevel]+"_"+Datasets[iDataset]+DatasetsNames[iDataset]+"_R="+Form("%.1f", arrayRadius[iRadius])+"_Pt_raw");

  TString textContext(contextCustomOneField(*texDatasetsComparisonCommonDenominator, ""));

  TString* yAxisLabel;
  yAxisLabel = texCount;
  if (normaliseDistribsBeforeUnfolding) {
    yAxisLabel = texJetPtYield_EventNorm;
  }

  Draw_TH1_Histogram(H1D_jetPt_raw, textContext, pdfName, texPtX, yAxisLabel, texCollisionDataInfo, drawnWindowUnfoldedMeasurement, legendPlacementAuto, contextPlacementAuto, "logy");
}

void Draw_Pt_spectrum_mcp(int iDataset, int iRadius, std::string options) {

  TH1D* H1D_jetPt_mcp_genBinning;
  TH1D* H1D_jetPt_mcp_recBinning;
  TH1D* H1D_jetPt_mcp_collection[2];
  Get_Pt_spectrum_mcp_genBinning(H1D_jetPt_mcp_genBinning, iDataset, iRadius, controlMC, options);
  Get_Pt_spectrum_mcp_recBinning(H1D_jetPt_mcp_recBinning, iDataset, iRadius, controlMC, options);
  H1D_jetPt_mcp_collection[0] = H1D_jetPt_mcp_genBinning;
  H1D_jetPt_mcp_collection[1] = H1D_jetPt_mcp_recBinning;


  TString* pdfName = new TString("jet_"+jetType[iJetType]+"_"+jetLevel[iJetLevel]+"_"+Datasets[iDataset]+DatasetsNames[iDataset]+"_R="+Form("%.1f", arrayRadius[iRadius])+"_Pt_mcp");

  TString textContext(contextCustomOneField(*texDatasetsComparisonCommonDenominator, ""));

  TString* yAxisLabel;
  yAxisLabel = texCount;
  if (normaliseDistribsBeforeUnfolding) {
    yAxisLabel = texJetPtYield_EventNorm;
  }
  TString genVsRecBinningLegend[2] = {"gen binning", "rec binning"};

  Draw_TH1_Histograms(H1D_jetPt_mcp_collection, genVsRecBinningLegend, 2, textContext, pdfName, texPtX, yAxisLabel, texCollisionDataInfo, drawnWindowAuto, legendPlacementAuto, contextPlacementAuto, "logy"); 
}

void Draw_Pt_spectrum_mcdMatched(int iDataset, int iRadius, std::string options) {

  TH1D* H1D_jetPt_mcdMatched;
  Get_Pt_spectrum_mcdMatched_genBinning(H1D_jetPt_mcdMatched, iDataset, iRadius, options);


  TString* pdfName = new TString("jet_"+jetType[iJetType]+"_"+jetLevel[iJetLevel]+"_"+Datasets[iDataset]+DatasetsNames[iDataset]+"_R="+Form("%.1f", arrayRadius[iRadius])+"_Pt_mcdMatched");

  TString textContext(contextCustomOneField(*texDatasetsComparisonCommonDenominator, ""));

  TString* yAxisLabel;
  yAxisLabel = texCount;
  if (normaliseDistribsBeforeUnfolding) {
    yAxisLabel = texJetPtYield_EventNorm;
  }

  Draw_TH1_Histogram(H1D_jetPt_mcdMatched, textContext, pdfName, texPtX, yAxisLabel, texCollisionDataInfo, drawnWindowUnfoldedMeasurement, legendPlacementAuto, contextPlacementAuto, "logy");
}


void Draw_Pt_efficiency_jets(int iDataset, int iRadius, std::string options) {
  TH1D* H1D_jetEfficiency;
  bool divideSuccess;

  divideSuccess = Get_Pt_JetEfficiency(H1D_jetEfficiency, iDataset, iRadius, options);
  
  TString* pdfName = new TString("jet_"+jetType[iJetType]+"_"+jetLevel[iJetLevel]+"_"+Datasets[iDataset]+DatasetsNames[iDataset]+"_R="+Form("%.1f", arrayRadius[iRadius])+"_Pt_efficiency");

  TString textContext(contextCustomOneField(*texDatasetsComparisonCommonDenominator, ""));

  if (divideSuccess){
    Draw_TH1_Histogram(H1D_jetEfficiency, textContext, pdfName, texPtJetGenX, texJetEfficiency, texCollisionDataInfo, drawnWindowAuto, legendPlacementAuto, contextPlacementAuto, "efficiency");
  }
}

void Draw_kinematicEfficiency(int iDataset, int iRadius, std::string options) {

  TH2D* H2D_jetPtResponseMatrix_fluctuations;
  TH2D* H2D_jetPtResponseMatrix_detectorResponse;
  TH2D* H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined_fineBinning;
  TH1D* H1D_kinematicEfficiency;

  TString partialUniqueSpecifier = Datasets[iDataset]+"_R="+Form("%.1f",arrayRadius[iRadius]);
  TString name_H1D_kinematicEfficiency = Datasets[iDataset]+"_R="+Form("%.1f",arrayRadius[iRadius]);

  Get_PtResponseMatrix_detectorResponse(H2D_jetPtResponseMatrix_detectorResponse, iDataset, iRadius);

  Get_PtResponseMatrix_Fluctuations(H2D_jetPtResponseMatrix_fluctuations, iDataset, iRadius);
  Get_PtResponseMatrix_DetectorAndFluctuationsCombined_fineBinning(H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined_fineBinning, H2D_jetPtResponseMatrix_detectorResponse, H2D_jetPtResponseMatrix_fluctuations, iDataset, iRadius, options);

  Get_ResponseMatrix_Pt_KinematicEffiency(H1D_kinematicEfficiency, H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined_fineBinning, name_H1D_kinematicEfficiency, iRadius);

  TString priorInfo = (TString)(TString)mergingPrior+"-"+(TString)unfoldingPrior;

  TString* pdfName = new TString("kinematicEfficiency_"+partialUniqueSpecifier+priorInfo);

  TString textContext(contextCustomOneField(*texDatasetsComparisonCommonDenominator, ""));

  Draw_TH1_Histogram(H1D_kinematicEfficiency, textContext, pdfName, texPtJetGenX, texJetKinematicEfficiency, texCollisionDataInfo, drawnWindowAuto, legendPlacementAuto, contextPlacementAuto, "");
}


void Draw_FakeRatio(int iDataset, int iRadius, std::string options) {
  TH1D* H1D_fakeRatio;

  TString partialUniqueSpecifier = Datasets[iDataset]+"_R="+Form("%.1f",arrayRadius[iRadius]);
  TString name_H1D_fakeRatio = Datasets[iDataset]+"_R="+Form("%.1f",arrayRadius[iRadius]);


  Get_Pt_JetFakes(H1D_fakeRatio, iDataset, iRadius, options);

  TString* pdfName = new TString("fakeRatio_"+partialUniqueSpecifier);

  TString textContext(contextCustomOneField(*texDatasetsComparisonCommonDenominator, ""));

  Draw_TH1_Histogram(H1D_fakeRatio, textContext, pdfName, texPtJetRecX, texFakeRatio, texCollisionDataInfo, drawnWindowUnfoldedMeasurement, legendPlacementAuto, contextPlacementAuto, "");
}

void Draw_ResponseMatrices_Fluctuations(int iDataset, int iRadius) {

  TH2D* H2D_jetPtResponseMatrix_fluctuations;

  Get_PtResponseMatrix_Fluctuations(H2D_jetPtResponseMatrix_fluctuations, iDataset, iRadius);

  TString priorInfo = (TString)(TString)mergingPrior+"-"+(TString)unfoldingPrior;

  struct stat st1{};
  if (stat("pdfFolder/ResponseMatrices", &st1) == -1) {
      mkdir("pdfFolder/ResponseMatrices", 0700);
  }
  struct stat st2{};
  if (stat("pngFolder/ResponseMatrices", &st2) == -1) {
      mkdir("pngFolder/ResponseMatrices", 0700);
  }

  TString* pdfName_logz = new TString("ResponseMatrices/responseMatrix_fluctuationsBackground_"+(TString)"_R="+Form("%.1f",arrayRadius[iRadius])+"_"+Datasets[iDataset]+DatasetsNames[iDataset]+"_"+priorInfo+"_logz");
  TString* pdfNameFullRes_logz = new TString("ResponseMatrices/responseMatrix_fluctuationsBackground_"+(TString)"_R="+Form("%.1f",arrayRadius[iRadius])+"_"+Datasets[iDataset]+DatasetsNames[iDataset]+"_"+priorInfo+"FullRes_logz");
  TString* pdfName = new TString("ResponseMatrices/responseMatrix_fluctuationsBackground_"+(TString)"_R="+Form("%.1f",arrayRadius[iRadius])+"_"+Datasets[iDataset]+DatasetsNames[iDataset]+"_"+priorInfo);
  TString* pdfNameFullRes = new TString("ResponseMatrices/responseMatrix_fluctuationsBackground_"+(TString)"_R="+Form("%.1f",arrayRadius[iRadius])+"_"+Datasets[iDataset]+DatasetsNames[iDataset]+"_"+priorInfo+"_FullRes");

  TString textContext = contextCustomTwoFields(*texDatasetsComparisonCommonDenominator, contextJetRadius(arrayRadius[iRadius]), "");


  double th2ContourCustom[1] = {0.000001}; // hardcoded at 10-6 for now
  int contourNumberCustom = 1;

  // Draw_TH2_Histogram(H2D_jetPtResponseMatrix_fluctuations, textContext, pdfName_logz, texPtJetBkgCorrX, texPtJetBkgFreeX, texCollisionDataInfo, drawnWindowAuto, th2ContourCustom, contourNumberCustom, "logz");
  Draw_TH2_Histogram(H2D_jetPtResponseMatrix_fluctuations, textContext, pdfName_logz, texPtJetBkgCorrX, texPtJetBkgFreeX, texCollisionDataInfo, drawnWindowAuto, th2ContoursNone, contourNumberNone, "logz");
  Draw_TH2_Histogram(H2D_jetPtResponseMatrix_fluctuations, textContext, pdfName, texPtJetBkgCorrX, texPtJetBkgFreeX, texCollisionDataInfo, drawnWindowAuto, th2ContoursNone, contourNumberNone, "");
}

void Draw_ResponseMatrices_detectorResponse(int iDataset, int iRadius) {

  TH2D* H2D_jetPtResponseMatrix_detectorResponse;
  cout << "Draw_ResponseMatrices_detectorResponse 1" << endl;
  Get_PtResponseMatrix_detectorResponse(H2D_jetPtResponseMatrix_detectorResponse, iDataset, iRadius);
  cout << "Draw_ResponseMatrices_detectorResponse 2" << endl;

  TString priorInfo = (TString)(TString)mergingPrior+"-"+(TString)unfoldingPrior;

  struct stat st1{};
  if (stat("pdfFolder/ResponseMatrices", &st1) == -1) {
      mkdir("pdfFolder/ResponseMatrices", 0700);
  }
  struct stat st2{};
  if (stat("pngFolder/ResponseMatrices", &st2) == -1) {
      mkdir("pngFolder/ResponseMatrices", 0700);
  }

  TString* pdfName = new TString("ResponseMatrices/responseMatrix_detectorEffects_"+jetType[iJetType]+"_"+Datasets[iDataset]+DatasetsNames[iDataset]+"_"+priorInfo);
  TString* pdfName_logz = new TString("ResponseMatrices/responseMatrix_detectorEffects_"+(TString)"_R="+Form("%.1f",arrayRadius[iRadius])+"_"+Datasets[iDataset]+DatasetsNames[iDataset]+"_"+priorInfo+"_logz");

  TString textContext = contextCustomTwoFields(*texDatasetsComparisonCommonDenominator, contextJetRadius(arrayRadius[iRadius]), "");

  Draw_TH2_Histogram(H2D_jetPtResponseMatrix_detectorResponse, textContext, pdfName, texPtJetRecX, texPtJetGenX, texCollisionMCInfo, drawnWindowAuto, th2ContoursNone, contourNumberNone, "");
  Draw_TH2_Histogram(H2D_jetPtResponseMatrix_detectorResponse, textContext, pdfName_logz, texPtJetRecX, texPtJetGenX, texCollisionMCInfo, drawnWindowAuto, th2ContoursNone, contourNumberNone, "logz");
}

void Draw_ResponseMatrices_DetectorAndFluctuationsCombined(int iDataset, int iRadius, std::string options) {

  TH2D* H2D_jetPtResponseMatrix_fluctuations;
  TH2D* H2D_jetPtResponseMatrix_detectorResponse;
  TH2D* H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined;


  Get_PtResponseMatrix_Fluctuations(H2D_jetPtResponseMatrix_fluctuations, iDataset, iRadius);
  Get_PtResponseMatrix_detectorResponse(H2D_jetPtResponseMatrix_detectorResponse, iDataset, iRadius);
  Get_PtResponseMatrix_DetectorAndFluctuationsCombined(H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined, H2D_jetPtResponseMatrix_detectorResponse, H2D_jetPtResponseMatrix_fluctuations, iDataset, iRadius, options);
  // FinaliseResponseMatrix(H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined, iDataset, iRadius, options);

  TString priorInfo = (TString)(TString)mergingPrior+"-"+(TString)unfoldingPrior;

  struct stat st1{};
  if (stat("pdfFolder/ResponseMatrices", &st1) == -1) {
      mkdir("pdfFolder/ResponseMatrices", 0700);
  }
  struct stat st2{};
  if (stat("pngFolder/ResponseMatrices", &st2) == -1) {
      mkdir("pngFolder/ResponseMatrices", 0700);
  }

  TString* pdfName = new TString("ResponseMatrices/responseMatrix_combined"+(TString)"_R="+Form("%.1f",arrayRadius[iRadius])+"_"+Datasets[iDataset]+DatasetsNames[iDataset]+"_"+priorInfo);
  TString* pdfName_logz = new TString("ResponseMatrices/responseMatrix_combined"+(TString)"_R="+Form("%.1f",arrayRadius[iRadius])+"_"+Datasets[iDataset]+DatasetsNames[iDataset]+"_"+priorInfo+"_logz");

  TString texCombinedMatrix = contextCustomOneField((TString)"Combined matrix - "+(TString)*texEnergy, "");
  TString textContextMatrixCombined = contextCustomFourFields((TString)"Detector response: "+(TString)*texCollisionMCType, "", (TString)"Fluctuations response: "+*texCollisionDataType, contextJetRadius(arrayRadius[iRadius]), "");

  TH2D* MatrixResponse = (TH2D*)GetTransposeHistogram(H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined).Clone("Draw_ResponseMatrices_DetectorAndFluctuationsCombined"+(TString)"_R="+Form("%.1f",arrayRadius[iRadius])+"_"+Datasets[iDataset]+DatasetsNames[iDataset]+"_"+priorInfo);

  Draw_TH2_Histogram(MatrixResponse, textContextMatrixCombined, pdfName, texPtJetGenX, texPtJetRecX, &texCombinedMatrix, drawnWindowAuto, th2ContoursNone, contourNumberNone, "");
  Draw_TH2_Histogram(MatrixResponse, textContextMatrixCombined, pdfName_logz, texPtJetGenX, texPtJetRecX, &texCombinedMatrix, drawnWindowAuto, th2ContoursNone, contourNumberNone, "logz");
}

void Draw_Pt_spectrum_unfolded(int iDataset, int iRadius, int unfoldParameterInput, std::string options) {

  TH1D* H1D_jetPt_measured;
  TH1D* H1D_jetPt_measured_genBinning;
  TH1D* H1D_jetPt_unfolded;
  TH1D* H1D_jetPt_unfoldedThenRefolded;
  TH1D* H1D_jetPt_unfoldedThenRefolded_RooUnfoldMethod;
  TH1D* H1D_jetPt_unfolded_mcpComp[2];
  TH1D* H1D_jetPt_unfolded_run2Comp[2];
  TH1D* H1D_jetPt_unfolded_measuredComp[2];
  TH1D* H1D_jetPt_unfolded_refoldedComp[3];
  TH1D* H1D_jetPt_mcp;
  TH1D* H1D_jetPt_run2;
  TH1D* H1D_jetPt_ratio_mcp;
  TH1D* H1D_jetPt_ratio_run2;
  TH1D* H1D_jetPt_ratio_measured;
  TH1D* H1D_jetPt_ratio_measuredRefolded[2];

  if (isDataPbPb && comparePbPbWithRun2) {
    H1D_jetPt_run2 = (TH1D*)((TH1D*)file_O2Analysis_run2ComparisonFile->Get("Bayesian_Unfoldediter15"))->Clone("H1D_jetPt_run2");
    int NcollRun2 = 4619963; // central (see Laura discussion mattermost) 
    H1D_jetPt_run2->Scale(1./NcollRun2);
  }

  bool divideSuccessMcp;
  bool divideSuccessRun2;
  bool divideSuccessMeasured;
  bool divideSuccessMeasuredRefolded[2];
  TString partialUniqueSpecifier;

  if (!useFineBinningTest) {
    Get_Pt_spectrum_mcp_genBinning(H1D_jetPt_mcp, iDataset, iRadius, controlMC, options);
  } else {
    Get_Pt_spectrum_mcp_fineBinning(H1D_jetPt_mcp, iDataset, iRadius, controlMC, options);
  }
  // H1D_jetPt_mcp = (TH1D*)H1D_jetPt_mcp_defaultBin->Rebin(nBinPtJetsGen[iRadius],"jetPt_mcp_rebinned_"+RadiusLegend[iRadius]+Datasets[iDataset]+DatasetsNames[iDataset], ptBinsJetsGen[iRadius]);

  int unfoldParameter;

  partialUniqueSpecifier = Datasets[iDataset]+"_R="+Form("%.1f",arrayRadius[iRadius]);
  

  unfoldParameter = Get_Pt_spectrum_unfolded(H1D_jetPt_unfolded, iDataset, iRadius, unfoldParameterInput, options).first;

  // comparison with raw measured
  if (!useFineBinningTest) {
    Get_Pt_spectrum_bkgCorrected_genBinning(H1D_jetPt_measured_genBinning, iDataset, iRadius, options);
  } else {
    Get_Pt_spectrum_bkgCorrected_fineBinning(H1D_jetPt_measured_genBinning, iDataset, iRadius, options);
  }

  H1D_jetPt_unfolded_measuredComp[0] = (TH1D*)H1D_jetPt_unfolded->Clone("H1D_jetPt_unfolded_measuredComp"+partialUniqueSpecifier);
  H1D_jetPt_unfolded_measuredComp[1] = (TH1D*)H1D_jetPt_measured_genBinning->Clone("H1D_jetPt_measured_genBinning_measuredComp"+partialUniqueSpecifier);
  H1D_jetPt_ratio_measured = (TH1D*)H1D_jetPt_measured_genBinning->Clone("H1D_jetPt_ratio_measured"+partialUniqueSpecifier);
  divideSuccessMeasured = H1D_jetPt_ratio_measured->Divide(H1D_jetPt_unfolded);

  // comparison with mcp truth
  H1D_jetPt_unfolded_mcpComp[0] = (TH1D*)H1D_jetPt_unfolded->Clone("H1D_jetPt_unfolded_mcpComp"+partialUniqueSpecifier);
  H1D_jetPt_unfolded_mcpComp[1] = (TH1D*)H1D_jetPt_mcp->Clone("H1D_jetPt_unfolded_mcpComp"+partialUniqueSpecifier);
  H1D_jetPt_ratio_mcp = (TH1D*)H1D_jetPt_mcp->Clone("H1D_jetPt_ratio_mcp"+partialUniqueSpecifier);
  divideSuccessMcp = H1D_jetPt_ratio_mcp->Divide(H1D_jetPt_unfolded);

  // comparison with run2 
  if (isDataPbPb && comparePbPbWithRun2) {
    H1D_jetPt_unfolded_run2Comp[0] = (TH1D*)H1D_jetPt_unfolded->Clone("H1D_jetPt_unfolded_run2Comp"+partialUniqueSpecifier);
    H1D_jetPt_unfolded_run2Comp[1] = (TH1D*)H1D_jetPt_run2->Clone("H1D_jetPt_unfolded_run2Comp"+partialUniqueSpecifier);
    H1D_jetPt_ratio_run2 = (TH1D*)H1D_jetPt_run2->Clone("H1D_jetPt_ratio_run2"+partialUniqueSpecifier);
    divideSuccessRun2 = H1D_jetPt_ratio_run2->Divide(H1D_jetPt_unfolded);
  }

  // comparison with refolded
  if (!useFineBinningTest) {
    Get_Pt_spectrum_bkgCorrected_recBinning(H1D_jetPt_measured, iDataset, iRadius, options);
  } else {
    Get_Pt_spectrum_bkgCorrected_fineBinning(H1D_jetPt_measured, iDataset, iRadius, options);
  }
  Get_Pt_spectrum_dataUnfoldedThenRefolded(H1D_jetPt_unfoldedThenRefolded, iDataset, iRadius, unfoldParameterInput, options);
  Get_Pt_spectrum_dataUnfoldedThenRefolded_RooUnfoldMethod(H1D_jetPt_unfoldedThenRefolded_RooUnfoldMethod, iDataset, iRadius, unfoldParameterInput, options);
  H1D_jetPt_unfolded_refoldedComp[0] = (TH1D*)H1D_jetPt_unfoldedThenRefolded->Clone("H1D_jetPt_refolded_refoldedComp"+partialUniqueSpecifier);
  H1D_jetPt_unfolded_refoldedComp[1] = (TH1D*)H1D_jetPt_unfoldedThenRefolded_RooUnfoldMethod->Clone("H1D_jetPt_unfoldedThenRefolded_RooUnfoldMethod"+partialUniqueSpecifier);
  H1D_jetPt_unfolded_refoldedComp[2] = (TH1D*)H1D_jetPt_measured->Clone("H1D_jetPt_measured_refoldedComp"+partialUniqueSpecifier);
  H1D_jetPt_ratio_measuredRefolded[0] = (TH1D*)H1D_jetPt_unfoldedThenRefolded->Clone("H1D_jetPt_ratio_refoldedComp"+partialUniqueSpecifier);
  H1D_jetPt_ratio_measuredRefolded[1] = (TH1D*)H1D_jetPt_unfoldedThenRefolded_RooUnfoldMethod->Clone("H1D_jetPt_ratio_refoldedComp_RooUnfoldMethod"+partialUniqueSpecifier);
  divideSuccessMeasuredRefolded[0] = H1D_jetPt_ratio_measuredRefolded[0]->Divide(H1D_jetPt_measured);
  divideSuccessMeasuredRefolded[1] = H1D_jetPt_ratio_measuredRefolded[1]->Divide(H1D_jetPt_measured);


  // H1D_jetPt_unfolded_refoldedComp_RooUnfoldMethod[0] = (TH1D*)H1D_jetPt_unfoldedThenRefolded_RooUnfoldMethod->Clone("H1D_jetPt_refolded_refoldedComp_RooUnfoldMethod"+partialUniqueSpecifier);
  // H1D_jetPt_unfolded_refoldedComp_RooUnfoldMethod[1] = (TH1D*)H1D_jetPt_measured_RooUnfoldMethod->Clone("H1D_jetPt_measured_refoldedComp_RooUnfoldMethod"+partialUniqueSpecifier);
  // H1D_jetPt_ratio_measuredRefolded_RooUnfoldMethod = (TH1D*)H1D_jetPt_unfoldedThenRefolded_RooUnfoldMethod->Clone("H1D_jetPt_ratio_refoldedComp_RooUnfoldMethod"+partialUniqueSpecifier);
  // divideSuccessMeasuredRefolded_RooUnfoldMethod = H1D_jetPt_ratio_measuredRefolded_RooUnfoldMethod->Divide(H1D_jetPt_measured);

  TString unfoldingCode;
  if (useManualRespMatrixSettingMethod){
    unfoldingCode = "myUnfold";
  } else {
    unfoldingCode = "joUnfold";
  }
  if (controlMC){
    unfoldingCode += "_controlMC";
  }
  TString unfoldingInfo = (TString)unfoldingMethod+"-k="+Form("%i", unfoldParameter)+"-"+(TString)mergingPrior+"-"+(TString)unfoldingPrior+"-"+unfoldingCode;

  struct stat st1{};
  if (stat("pdfFolder/IterationsDump", &st1) == -1) {
      mkdir("pdfFolder/IterationsDump", 0700);
  }
  struct stat st2{};
  if (stat("pngFolder/IterationsDump", &st2) == -1) {
      mkdir("pngFolder/IterationsDump", 0700);
  }

  TString textContext = contextCustomTwoFields(*texDatasetsComparisonCommonDenominator, contextJetRadius(arrayRadius[iRadius]), "");

  TString dummyLegend[1] = {""};

  TString* yAxisLabel = texCount;
  if (normaliseDistribsBeforeUnfolding || normaliseDistribsAfterUnfolding) { //should probably check if having both on doesn't lead to double normalisation
    yAxisLabel = texJetPtYield_EventNorm;
  }

  TString pdfTitleBase = (TString)"IterationsDump/jet_";//+Datasets[iDataset]+DatasetsNames[iDataset]+"_R="+Form("%.1f", arrayRadius[iRadius])+"_Pt_unfolded_";
  // std::array<std::array<float, 2>, 2> drawnWindow = {{{ptWindowDisplay[0], ptWindowDisplay[1]}, {-999, -999}}}; // {{xmin, xmax}, {ymin, ymax}}

    // comparison with raw measured
  TString unfoldedMeasuredCompLegend[2] = {"unfolded data", "measured raw (gen binning)"};
  TString* pdfName_measuredComp = new TString(pdfTitleBase+unfoldingInfo+"_measuredComp");
  Draw_TH1_Histograms(H1D_jetPt_unfolded_measuredComp, unfoldedMeasuredCompLegend, 2, textContext, pdfName_measuredComp, texPtX, yAxisLabel, texCollisionDataInfo, drawnWindowUnfoldedMeasurement, legendPlacementAuto, contextPlacementAuto, "logy");
  if (divideSuccessMeasured){
    TString* pdfName_ratio_measured = new TString(pdfTitleBase+unfoldingInfo+"_ratioMeasured");
    Draw_TH1_Histogram(H1D_jetPt_ratio_measured, textContext, pdfName_ratio_measured, texPtX, texRatioMeasuredUnfolded, texCollisionDataInfo, drawnWindowUnfoldedMeasurement, legendPlacementAuto, contextPlacementAuto, "zoomToOneLarge,ratioLine");
  }

    // comparison with mcp truth
  TString unfoldedTruthCompLegend[2] = {"unfolded data", "mcp truth"};
  TString* pdfName_mcpComp = new TString(pdfTitleBase+unfoldingInfo+"_mcpComp");
  Draw_TH1_Histograms(H1D_jetPt_unfolded_mcpComp, unfoldedTruthCompLegend, 2, textContext, pdfName_mcpComp, texPtX, yAxisLabel, texCollisionDataInfo, drawnWindowUnfoldedMeasurement, legendPlacementAuto, contextPlacementAuto, "logy");
  if (divideSuccessMcp){
    TString* pdfName_ratio_mcp = new TString(pdfTitleBase+unfoldingInfo+"_ratioMcp");
    Draw_TH1_Histogram(H1D_jetPt_ratio_mcp, textContext, pdfName_ratio_mcp, texPtX, texRatioMcpUnfolded, texCollisionDataInfo, drawnWindowUnfoldedMeasurement, legendPlacementAuto, contextPlacementAuto, "zoomToOneLarge,ratioLine");
  }


  // comparison with refolded
  TString unfoldedRefoldedCompLegend[3] = {"refolded manually", "refolded roounfold (noErrors)", "measured"};
  TString* pdfName_refoldedComp = new TString(pdfTitleBase+unfoldingInfo+"_RefoldedComp");
  Draw_TH1_Histograms(H1D_jetPt_unfolded_refoldedComp, unfoldedRefoldedCompLegend, 3, textContext, pdfName_refoldedComp, texPtX, yAxisLabel, texCollisionDataInfo, drawnWindowUnfoldedMeasurement, legendPlacementAuto, contextPlacementAuto, "logy");
  if (divideSuccessMeasuredRefolded[0] && divideSuccessMeasuredRefolded[1]) {
    TString* pdfName_ratio_refoldedComp = new TString(pdfTitleBase+unfoldingInfo+"_ratioRefoldedUnfolded");
    Draw_TH1_Histograms(H1D_jetPt_ratio_measuredRefolded, unfoldedRefoldedCompLegend, 2, textContext, pdfName_ratio_refoldedComp, texPtX, texRatioRefoldedMeasured, texCollisionDataInfo, drawnWindowUnfoldedMeasurement, legendPlacementAuto, contextPlacementAuto, "zoomToOneLarge,ratioLine");
  }


  // comparison with Run 2
  if (isDataPbPb && comparePbPbWithRun2) {
  TString unfoldedRun2CompLegend[2] = {"unfolded Run3", "unfolded Run2 0-10%"};
    TString* pdfName_run2Comp = new TString(pdfTitleBase+unfoldingInfo+"_run2Comp");
    Draw_TH1_Histograms(H1D_jetPt_unfolded_run2Comp, unfoldedRun2CompLegend, 2, textContext, pdfName_run2Comp, texPtX, yAxisLabel, texCollisionDataInfo, drawnWindowUnfoldedMeasurement, legendPlacementAuto, contextPlacementAuto, "logy");
    if (divideSuccessRun2){
      TString* pdfName_ratio_run2 = new TString(pdfTitleBase+unfoldingInfo+"_ratioRun2");
      Draw_TH1_Histogram(H1D_jetPt_ratio_run2, textContext, pdfName_ratio_run2, texPtX, texRatioRun2Unfolded, texCollisionDataInfo, drawnWindowUnfoldedMeasurement, legendPlacementAuto, contextPlacementAuto, "zoomToOneLarge,ratioLine");
    }
  }

}

void Draw_Pt_spectrum_unfolded_parameterVariation(int iDataset, int iRadius, int unfoldIterationMin, int unfoldIterationMax, int step, std::string options) {

  const int nUnfoldIteration = std::floor((unfoldIterationMax - unfoldIterationMin + 1)/step);

  TH1D* H1D_jetPt_measured;
  TH1D* H1D_jetPt_measured_genBinning;
  TH1D* H1D_jetPt_unfolded[nUnfoldIteration];
  TH1D* H1D_jetPt_unfoldedThenRefolded[nUnfoldIteration];
  TH1D* H1D_jetPt_unfolded_mcpComp[nUnfoldIteration+1];
  TH1D* H1D_jetPt_unfolded_measuredComp[nUnfoldIteration+1];
  TH1D* H1D_jetPt_unfolded_refoldedComp[nUnfoldIteration+1];
  TH1D* H1D_jetPt_mcp;
  TH1D* H1D_jetPt_ratio_mcp[nUnfoldIteration];
  TH1D* H1D_jetPt_ratio_measured[nUnfoldIteration];
  TH1D* H1D_jetPt_ratio_measuredRefolded[nUnfoldIteration];

  bool divideSuccessMcp[nUnfoldIteration];
  bool divideSuccessMeasured[nUnfoldIteration];
  bool divideSuccessMeasuredRefolded[nUnfoldIteration];
  TString partialUniqueSpecifier;

  Get_Pt_spectrum_mcp_genBinning(H1D_jetPt_mcp, iDataset, iRadius, controlMC, options);
  // H1D_jetPt_mcp = (TH1D*)H1D_jetPt_mcp_defaultBin->Rebin(nBinPtJetsGen[iRadius],"jetPt_mcp_rebinned_"+RadiusLegend[iRadius]+Datasets[iDataset]+DatasetsNames[iDataset], ptBinsJetsGen[iRadius]);

  partialUniqueSpecifier = Datasets[iDataset]+"_R="+Form("%.1f",arrayRadius[iRadius]);
  Get_Pt_spectrum_bkgCorrected_recBinning(H1D_jetPt_measured, iDataset, iRadius, options);
  Get_Pt_spectrum_bkgCorrected_genBinning(H1D_jetPt_measured_genBinning, iDataset, iRadius, options);

  int unfoldIterations = 0;

  for(int iUnfoldIteration = 0; iUnfoldIteration < nUnfoldIteration; iUnfoldIteration++){

    unfoldIterations = iUnfoldIteration * step + unfoldIterationMin; 
    // unfoldingIterationLegend[iUnfoldIteration] = unfoldIterations;

    cout << "((((((((((((()))))))))))))" << endl;
    cout << "Iteration "<< iUnfoldIteration << endl;
    cout << "((((((((((((()))))))))))))" << endl;
    Get_Pt_spectrum_unfolded(H1D_jetPt_unfolded[iUnfoldIteration], iDataset, iRadius, unfoldIterations, options);

    // comparison with raw measured
    H1D_jetPt_unfolded_measuredComp[iUnfoldIteration] = (TH1D*)H1D_jetPt_unfolded[iUnfoldIteration]->Clone("H1D_jetPt_unfolded_measuredComp"+partialUniqueSpecifier);
    H1D_jetPt_ratio_measured[iUnfoldIteration] = (TH1D*)H1D_jetPt_measured_genBinning->Clone("H1D_jetPt_ratio_measured"+partialUniqueSpecifier);
    divideSuccessMeasured[iUnfoldIteration] = H1D_jetPt_ratio_measured[iUnfoldIteration]->Divide(H1D_jetPt_unfolded[iUnfoldIteration]);

    // comparison with mcp truth
    H1D_jetPt_unfolded_mcpComp[iUnfoldIteration] = (TH1D*)H1D_jetPt_unfolded[iUnfoldIteration]->Clone("H1D_jetPt_unfolded_mcpComp"+partialUniqueSpecifier);
    H1D_jetPt_ratio_mcp[iUnfoldIteration] = (TH1D*)H1D_jetPt_mcp->Clone("H1D_jetPt_ratio_mcp"+partialUniqueSpecifier);
    divideSuccessMcp[iUnfoldIteration] = H1D_jetPt_ratio_mcp[iUnfoldIteration]->Divide(H1D_jetPt_unfolded[iUnfoldIteration]);


    // comparison with refolded
    Get_Pt_spectrum_dataUnfoldedThenRefolded(H1D_jetPt_unfoldedThenRefolded[iUnfoldIteration], iDataset, iRadius, unfoldIterations, options);
    H1D_jetPt_unfolded_refoldedComp[iUnfoldIteration] = (TH1D*)H1D_jetPt_unfoldedThenRefolded[iUnfoldIteration]->Clone("H1D_jetPt_refolded_refoldedComp"+partialUniqueSpecifier);
    H1D_jetPt_ratio_measuredRefolded[iUnfoldIteration] = (TH1D*)H1D_jetPt_unfoldedThenRefolded[iUnfoldIteration]->Clone("H1D_jetPt_ratio_refoldedComp"+partialUniqueSpecifier);
    divideSuccessMeasuredRefolded[iUnfoldIteration] = H1D_jetPt_ratio_measuredRefolded[iUnfoldIteration]->Divide(H1D_jetPt_measured);
  }
  H1D_jetPt_unfolded_measuredComp[nUnfoldIteration] = (TH1D*)H1D_jetPt_measured_genBinning->Clone("H1D_jetPt_measured_genBinning_measuredComp"+partialUniqueSpecifier);
  H1D_jetPt_unfolded_mcpComp[nUnfoldIteration] = (TH1D*)H1D_jetPt_mcp->Clone("H1D_jetPt_unfolded_mcpComp"+partialUniqueSpecifier);
  H1D_jetPt_unfolded_refoldedComp[nUnfoldIteration] = (TH1D*)H1D_jetPt_measured->Clone("H1D_jetPt_measured_refoldedComp"+partialUniqueSpecifier);


  TString unfoldingInfo = (TString)unfoldingMethod+"_"+(TString)unfoldingPrior+"_kmax="+Form("%i", unfoldIterationMax);

  TString unfoldingIterationLegend[nUnfoldIteration+1]; IterationLegend(unfoldingIterationLegend, unfoldIterationMin, unfoldIterationMax, step);

  
  partialUniqueSpecifier = Datasets[iDataset]+"_R="+Form("%.1f",arrayRadius[iRadius]);

  TString textContext = contextCustomTwoFields(*texDatasetsComparisonCommonDenominator, contextJetRadius(arrayRadius[iRadius]), "");

  TString* pdfName = new TString("jet_"+jetType[iJetType]+"_"+jetLevel[iJetLevel]+"_"+partialUniqueSpecifier+"_Pt_unfolded_"+unfoldingInfo);

  TString* yAxisLabel = texCount;
  if (normaliseDistribsBeforeUnfolding || normaliseDistribsAfterUnfolding) { //should probably check if having both on doesn't lead to double normalisation
    yAxisLabel = texJetPtYield_EventNorm;
  }
  Draw_TH1_Histograms(H1D_jetPt_unfolded, unfoldingIterationLegend, nUnfoldIteration, textContext, pdfName, texPtX, yAxisLabel, texCollisionDataInfo, drawnWindowUnfoldedMeasurement, legendPlacementAuto, contextPlacementAuto, "logy");


    // comparison with raw measured
  // TString unfoldedMeasuredCompLegend[2] = {"unfolded data", "measured raw (gen binning)"};
  unfoldingIterationLegend[nUnfoldIteration] = (TString)"raw measured";
  TString* pdfName_measuredComp = new TString("jet_"+jetType[iJetType]+"_"+jetLevel[iJetLevel]+"_"+partialUniqueSpecifier+"_Pt_unfolded_"+unfoldingInfo+"_measuredComp");
  Draw_TH1_Histograms(H1D_jetPt_unfolded_measuredComp, unfoldingIterationLegend, nUnfoldIteration+1, textContext, pdfName_measuredComp, texPtX, yAxisLabel, texCollisionDataInfo, drawnWindowUnfoldedMeasurement, legendPlacementAuto, contextPlacementAuto, "logy");

  bool divideSuccessMeasured_boolsum = true;
  for(int iUnfoldIteration = 0; iUnfoldIteration < nUnfoldIteration; iUnfoldIteration++){
    if (!(divideSuccessMeasured[iUnfoldIteration])) {
      divideSuccessMeasured_boolsum = false;
    }
  }
  if (divideSuccessMeasured_boolsum){
    TString* pdfName_ratio_measured = new TString("jet_"+jetType[iJetType]+"_"+jetLevel[iJetLevel]+"_"+partialUniqueSpecifier+"_Pt_unfolded_"+unfoldingInfo+"_ratioMeasured");
    Draw_TH1_Histograms(H1D_jetPt_ratio_measured, unfoldingIterationLegend, nUnfoldIteration, textContext, pdfName_ratio_measured, texPtX, texRatioMeasuredUnfolded, texCollisionDataInfo, drawnWindowUnfoldedMeasurement, legendPlacementAuto, contextPlacementAuto, "zoomToOneLarge,ratioLine");
  }


    // comparison with mcp truth
  // TString unfoldedTruthCompLegend[2] = {"unfolded data", "mcp truth"};
  unfoldingIterationLegend[nUnfoldIteration] = (TString)"mcp";
  TString* pdfName_mcpComp = new TString("jet_"+jetType[iJetType]+"_"+jetLevel[iJetLevel]+"_"+partialUniqueSpecifier+"_Pt_unfolded_"+unfoldingInfo+"_mcpComp");
  Draw_TH1_Histograms(H1D_jetPt_unfolded_mcpComp, unfoldingIterationLegend, nUnfoldIteration+1, textContext, pdfName_mcpComp, texPtX, yAxisLabel, texCollisionDataInfo, drawnWindowUnfoldedMeasurement, legendPlacementAuto, contextPlacementAuto, "logy");

  bool divideSuccessMcp_boolsum = true;
  for(int iUnfoldIteration = 0; iUnfoldIteration < nUnfoldIteration; iUnfoldIteration++){
    if (!(divideSuccessMcp[iUnfoldIteration])) {
      divideSuccessMcp_boolsum = false;
    }
  }
  if (divideSuccessMcp_boolsum){
    TString* pdfName_ratio_mcp = new TString("jet_"+jetType[iJetType]+"_"+jetLevel[iJetLevel]+"_"+partialUniqueSpecifier+"_Pt_unfolded_"+unfoldingInfo+"_ratioMcp");
    Draw_TH1_Histograms(H1D_jetPt_ratio_mcp, unfoldingIterationLegend, nUnfoldIteration, textContext, pdfName_ratio_mcp, texPtX, texRatioMcpUnfolded, texCollisionDataInfo, drawnWindowUnfoldedMeasurement, legendPlacementAuto, contextPlacementAuto, "zoomToOneLarge,ratioLine");
  }

  // comparison with refolded
  // TString unfoldedRefoldedCompLegend[2] = {"refolded", "measured"};
  unfoldingIterationLegend[nUnfoldIteration] = (TString)"measured";
  for(int iIteration = 0; iIteration < nUnfoldIteration; iIteration++){
    unfoldingIterationLegend[iIteration] += (TString)" refolded";
  }
  TString* pdfName_refoldedComp = new TString("jet_"+jetType[iJetType]+"_"+jetLevel[iJetLevel]+"_"+partialUniqueSpecifier+"_Pt_unfolded_"+unfoldingInfo+"_RefoldedComp");
  Draw_TH1_Histograms(H1D_jetPt_unfolded_refoldedComp, unfoldingIterationLegend, nUnfoldIteration+1, textContext, pdfName_refoldedComp, texPtX, yAxisLabel, texCollisionDataInfo, drawnWindowUnfoldedMeasurement, legendPlacementAuto, contextPlacementAuto, "logy");

  bool divideSuccessRefoldedComp_boolsum = true;
  for(int iUnfoldIteration = 0; iUnfoldIteration < nUnfoldIteration; iUnfoldIteration++){
    if (!(divideSuccessMeasuredRefolded[iUnfoldIteration])) {
      divideSuccessRefoldedComp_boolsum = false;
    }
  }
  if (divideSuccessRefoldedComp_boolsum){
    TString* pdfName_ratio_refoldedComp = new TString("jet_"+jetType[iJetType]+"_"+jetLevel[iJetLevel]+"_"+partialUniqueSpecifier+"_Pt_unfolded_"+unfoldingInfo+"_ratioRefoldedUnfolded");
    Draw_TH1_Histograms(H1D_jetPt_ratio_measuredRefolded, unfoldingIterationLegend, nUnfoldIteration, textContext, pdfName_ratio_refoldedComp, texPtX, texRatioRefoldedMeasured, texCollisionDataInfo, drawnWindowUnfoldedMeasurement, legendPlacementAuto, contextPlacementAuto, "zoomToOneLarge,ratioLine");
  }
}




// rename refoldedUnfolded as closure test?
// and try and spend 15 min to clean hist names for the spectrum analysis



// WARNING FOR EFFICIENCIES I SHOULD REREAD THIS BELOW!!
// hMcEfficiency_vsPt->Divide(hMcSignalCount_vsPt,TrueV0PtSpectrum_AnalysisBins, 1., 1., "b"); // option b for binomial because efficiency: https://twiki.cern.ch/twiki/bin/view/ALICE/PWGLFPAGSTRANGENESSEfficiency
