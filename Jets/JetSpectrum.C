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
#include <stdlib.h>     /* abort, NULL */
using namespace std;

// Misc utilities
void SetStyle(Bool_t graypalette=kFALSE);
void LoadLibs();

//////////// Pt Spectrum analysis functions
void Get_Pt_spectrum_bkgCorrected_recBinning(TH1D* &H1D_jetPt_raw, int iDataset, int iRadius, std::string options);
void Get_Pt_spectrum_bkgCorrected_genBinning(TH1D* &H1D_jetPt_raw, int iDataset, int iRadius, std::string options);
void Get_Pt_spectrum_bkgCorrected_fineBinning(TH1D* &H1D_jetPt_raw, int iDataset, int iRadius, std::string options);
void Get_Pt_spectrum_mcp_genBinning(TH1D* &H1D_jetPt_mcp, int iDataset, int iRadius, bool fcontrolMC, std::string options);
void Get_Pt_spectrum_mcp_recBinning(TH1D* &H1D_jetPt_mcp, int iDataset, int iRadius, bool fcontrolMC, std::string options);
void Get_Pt_spectrum_mcdMatched_genBinning(TH1D* &H1D_jetPt_mcd, int iDataset, int iRadius, std::string options);
void Get_Pt_spectrum_mcpMatched_genBinning(TH1D* &H1D_jetPt_mcp, int iDataset, int iRadius, std::string options);
std::pair<int, RooUnfold*> Get_Pt_spectrum_unfolded(TH1D* &H1D_jetPt_unfolded, int iDataset, int iRadius, int unfoldParameterInput, std::string options);
void Get_Pt_spectrum_bkgCorrected_recBinning_preWidthScalingAndEvtNorm(TH1D* &H1D_jetPt_raw, int iDataset, int iRadius, std::string options);
void Get_Pt_spectrum_bkgCorrected_genBinning_preWidthScalingAndEvtNorm(TH1D* &H1D_jetPt_raw, int iDataset, int iRadius, std::string options);
void Get_Pt_spectrum_bkgCorrected_fineBinning_preWidthScalingAndEvtNorm(TH1D* &H1D_jetPt_raw, int iDataset, int iRadius, std::string options);
void Get_Pt_spectrum_mcp_genBinning_preWidthScalingAndEvtNorm(TH1D* &H1D_jetPt_mcp, int iDataset, int iRadius, bool fcontrolMC, std::string options);
void Get_Pt_spectrum_mcp_recBinning_preWidthScalingAndEvtNorm(TH1D* &H1D_jetPt_mcp, int iDataset, int iRadius, bool fcontrolMC, std::string options);
void Get_Pt_spectrum_mcdMatched_genBinning_preWidthScalingAndEvtNorm(TH1D* &H1D_jetPt_mcd, int iDataset, int iRadius, std::string options);
void Get_Pt_spectrum_mcpMatched_genBinning_preWidthScalingAndEvtNorm(TH1D* &H1D_jetPt_mcpMatched, int iDataset, int iRadius, std::string options);
void Get_Pt_spectrum_mcp_fineBinning_preWidthScalingAndEvtNorm(TH1D* &H1D_jetPt_mcp, int iDataset, int iRadius, bool fcontrolMC, std::string options);
void Get_Pt_spectrum_mcd_fineBinning_preWidthScalingAndEvtNorm(TH1D* &H1D_jetPt_mcp, int iDataset, int iRadius, std::string options);
void Get_Pt_spectrum_mcd_recBinning_preWidthScalingAndEvtNorm(TH1D* &H1D_jetPt_mcp, int iDataset, int iRadius, std::string options);
void Get_Pt_spectrum_mcp_fineBinning(TH1D* &H1D_jetPt_mcp, int iDataset, int iRadius, bool fcontrolMC, std::string options);
void Get_Pt_spectrum_mcd_fineBinning(TH1D* &H1D_jetPt_mcp, int iDataset, int iRadius, std::string options);
void Get_Pt_spectrum_mcd_recBinning(TH1D* &H1D_jetPt_mcp, int iDataset, int iRadius, std::string options);
std::pair<int, RooUnfold*> Get_Pt_spectrum_unfolded_preWidthScalingAndEvtNorm(TH1D* &H1D_jetPt_unfolded, int iDataset, int iRadius, int unfoldParameterInput, std::string options);
void Get_Pt_spectrum_unfoldedThenRefolded_preWidthScalingAndEvtNorm(TH1D* &H1D_jetPt_mcp_unfoldedThenRefolded, int iDataset, int iRadius, int unfoldParameterInput, std::string options);
void Get_Pt_spectrum_unfoldedThenRefolded(TH1D* &H1D_jetPt_mcp_unfoldedThenRefolded, int iDataset, int iRadius, int unfoldParameterInput, std::string options);

void Get_PtResponseMatrix_Fluctuations(TH2D* &H2D_jetPtResponseMatrix_fluctuations, int iDataset, int iRadius);
void Get_PtResponseMatrix_detectorResponse(TH2D* &H2D_jetPtResponseMatrix_detectorResponse, int iDataset, int iRadius);
void Get_PtResponseMatrix_DetectorAndFluctuationsCombined(TH2D* &H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined, TH2D* H2D_jetPtResponseMatrix_detectorResponse, TH2D* H2D_jetPtResponseMatrix_fluctuations, int iDataset, int iRadius, std::string options);
void Get_PtResponseMatrix_DetectorAndFluctuationsCombined_fineBinning(TH2D* &H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined_fineBinning, TH2D* H2D_jetPtResponseMatrix_detectorResponse, TH2D* H2D_jetPtResponseMatrix_fluctuations, int iDataset, int iRadius, std::string options);
void ReweightResponseMatrixWithPrior(TH2D* &H2D_jetPtResponseMatrix, int iDataset, int iRadius, std::string options);


bool Get_Pt_JetEfficiency(TH1D* &H1D_jetEfficiency, bool* divideSuccess, int iDataset, int iRadius, std::string options);
void Get_ResponseMatrix_Pt_KinematicEffiency(TH1D* &H1D_kinematicEfficiency, TH2D* H2D_jetPtResponseMatrix, TString name_H1D_kinematicEfficiency, int iRadius);

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
  snprintf(optionsAnalysis, sizeof(optionsAnalysis), "%s,%s,%s", mergingPrior, unfoldingPrior, unfoldingMethod);
  cout << "Analysis options are: " << optionsAnalysis << endl;

  int iDataset = 0;
  int iRadius = 0;

  // // find a way to input mcpPrior/mcdPrior and bayes/svd as a variables rather than typed out like this
  // Draw_ResponseMatrices_Fluctuations(iDataset, iRadius);
  // Draw_ResponseMatrices_detectorResponse(iDataset, iRadius);
  // Draw_ResponseMatrices_DetectorAndFluctuationsCombined(iDataset, iRadius, optionsAnalysis);

  // // // Draw_Pt_spectrum_unfolded_FluctResponseOnly(iDataset, iRadius, optionsAnalysis); // NOT FIXED YET - result meaningless
  // Draw_Pt_spectrum_raw(iDataset, iRadius, optionsAnalysis);
  // Draw_Pt_spectrum_mcp(iDataset, iRadius, optionsAnalysis);
  // Draw_Pt_spectrum_mcdMatched(iDataset, iRadius, optionsAnalysis);

  for (int iDataset = 0; iDataset < nDatasets; iDataset++) {
    // Draw_Pt_efficiency_jets(iDataset, iRadius, optionsAnalysis);
  }
  

  // Draw_kinematicEfficiency(iDataset, iRadius, optionsAnalysis);
  // Draw_FakeRatio(iDataset, iRadius, optionsAnalysis);

  // int unfoldParameterInput = 20;
  // Draw_Pt_spectrum_unfolded(iDataset, iRadius, unfoldParameterInput, optionsAnalysis);

  int unfoldParameterInputMin = 6;
  int unfoldParameterInputMax = 15;
  int unfoldParameterInputStep = 2;
  Draw_Pt_spectrum_unfolded_parameterVariation(iDataset, iRadius, unfoldParameterInputMin, unfoldParameterInputMax, unfoldParameterInputStep, optionsAnalysis);
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
/////////////////////////////////////////////////////////////////////////// Spectrum getting functions ///////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


void Get_Pt_spectrum_bkgCorrected_recBinning_preWidthScalingAndEvtNorm(TH1D* &H1D_jetPt, int iDataset, int iRadius, __attribute__ ((unused)) std::string options) {
  TH1D* H1D_jetPt_defaultBin;

  TString histogramName = "";
  if (doBkgSubtractionInData) {
    histogramName = "h_jet_pt_rhoareasubtracted";
  } else {
    histogramName = "h_jet_pt";
  }

  if (!controlMC) {
    H1D_jetPt_defaultBin = (TH1D*)((TH1D*)file_O2Analysis_list[iDataset]->Get(analysisWorkflowData+"/"+histogramName))->Clone("Get_Pt_spectrum_bkgCorrected_genBinning"+Datasets[iDataset]+DatasetsNames[iDataset]);
  } else {
    H1D_jetPt_defaultBin = (TH1D*)((TH1D*)file_O2Analysis_ppSimDetectorEffect_unfoldingControl->Get(analysisWorkflowData+"/"+histogramName))->Clone("Get_Pt_spectrum_bkgCorrected_genBinning_unfoldingControl"+Datasets[iDataset]+DatasetsNames[iDataset]);
  }
  H1D_jetPt_defaultBin->Sumw2();
  
  H1D_jetPt = (TH1D*)H1D_jetPt_defaultBin->Rebin(nBinPtJetsRec[iRadius],"jetPt_bkgCorrected_rebinned_recBinning_"+Datasets[iDataset]+DatasetsNames[iDataset]+RadiusLegend[iRadius], ptBinsJetsRec[iRadius]);

  if (doWidthScalingEarly) {
    TransformRawHistToYield(H1D_jetPt);
  }
}


void Get_Pt_spectrum_bkgCorrected_genBinning_preWidthScalingAndEvtNorm(TH1D* &H1D_jetPt, int iDataset, int iRadius, __attribute__ ((unused)) std::string options) {
  TH1D* H1D_jetPt_defaultBin;

  TString histogramName = "";
  if (doBkgSubtractionInData) {
    histogramName = "h_jet_pt_rhoareasubtracted";
  } else {
    histogramName = "h_jet_pt";
  }

  if (!controlMC) {
    H1D_jetPt_defaultBin = (TH1D*)((TH1D*)file_O2Analysis_list[iDataset]->Get(analysisWorkflowData+"/"+histogramName))->Clone("Get_Pt_spectrum_bkgCorrected_genBinning"+Datasets[iDataset]+DatasetsNames[iDataset]);
  } else {
    H1D_jetPt_defaultBin = (TH1D*)((TH1D*)file_O2Analysis_ppSimDetectorEffect_unfoldingControl->Get(analysisWorkflowData+"/"+histogramName))->Clone("Get_Pt_spectrum_bkgCorrected_genBinning_unfoldingControl"+Datasets[iDataset]+DatasetsNames[iDataset]);
  }
  H1D_jetPt_defaultBin->Sumw2();

  H1D_jetPt = (TH1D*)H1D_jetPt_defaultBin->Rebin(nBinPtJetsGen[iRadius],"jetPt_bkgCorrected_rebinned_genBinning_"+Datasets[iDataset]+DatasetsNames[iDataset]+RadiusLegend[iRadius], ptBinsJetsGen[iRadius]);

  if (doWidthScalingEarly) {
    TransformRawHistToYield(H1D_jetPt);
  }
}

void Get_Pt_spectrum_bkgCorrected_fineBinning_preWidthScalingAndEvtNorm(TH1D* &H1D_jetPt, int iDataset, int iRadius, __attribute__ ((unused)) std::string options) {
  TH1D* H1D_jetPt_defaultBin;

  TString histogramName = "";
  if (doBkgSubtractionInData) {
    histogramName = "h_jet_pt_rhoareasubtracted";
  } else {
    histogramName = "h_jet_pt";
  }

  if (!controlMC) {
    H1D_jetPt_defaultBin = (TH1D*)((TH1D*)file_O2Analysis_list[iDataset]->Get(analysisWorkflowData+"/"+histogramName))->Clone("Get_Pt_spectrum_bkgCorrected_genBinning"+Datasets[iDataset]+DatasetsNames[iDataset]);
  } else {
    H1D_jetPt_defaultBin = (TH1D*)((TH1D*)file_O2Analysis_ppSimDetectorEffect_unfoldingControl->Get(analysisWorkflowData+"/"+histogramName))->Clone("Get_Pt_spectrum_bkgCorrected_genBinning_unfoldingControl"+Datasets[iDataset]+DatasetsNames[iDataset]);
  }
  H1D_jetPt_defaultBin->Sumw2();

  H1D_jetPt = (TH1D*)H1D_jetPt_defaultBin->Rebin(nBinPtJetsFine[iRadius],"jetPt_bkgCorrected_rebinned_fineBinning_"+Datasets[iDataset]+DatasetsNames[iDataset]+RadiusLegend[iRadius], ptBinsJetsFine[iRadius]);

  if (doWidthScalingEarly) {
    TransformRawHistToYield(H1D_jetPt);
  }
}

void Get_Pt_spectrum_mcp_genBinning_preWidthScalingAndEvtNorm(TH1D* &H1D_jetPt_mcp, int iDataset, int iRadius, bool fcontrolMC, __attribute__ ((unused)) std::string options) {
  TH1D* H1D_jetPt_defaultBin;

  TString histogramName = "";
  if (doBkgSubtractionInMC) {
    histogramName = "h_jet_pt_part_rhoareasubtracted";
  } else {
    histogramName = "h_jet_pt_part";
  }

  if (!fcontrolMC){
    H1D_jetPt_defaultBin = (TH1D*)((TH1D*)file_O2Analysis_MCfileForMatrix->Get(analysisWorkflowMC+"/"+histogramName))->Clone("Get_Pt_spectrum_mcp_genBinning"+Datasets[iDataset]+DatasetsNames[iDataset]);
  } else {
    H1D_jetPt_defaultBin = (TH1D*)((TH1D*)file_O2Analysis_ppSimDetectorEffect_unfoldingControl->Get(analysisWorkflowMC+"/"+histogramName))->Clone("Get_Pt_spectrum_mcp_genBinning_unfoldingControl"+Datasets[iDataset]+DatasetsNames[iDataset]);
  }
  H1D_jetPt_defaultBin->Sumw2();

  H1D_jetPt_mcp = (TH1D*)H1D_jetPt_defaultBin->Rebin(nBinPtJetsGen[iRadius],"jetPt_mcp_rebinned_genBinning_"+Datasets[iDataset]+DatasetsNames[iDataset]+RadiusLegend[iRadius]+Datasets[iDataset]+DatasetsNames[iDataset], ptBinsJetsGen[iRadius]);

  if (doWidthScalingEarly) {
    TransformRawHistToYield(H1D_jetPt_mcp);
  }
}

void Get_Pt_spectrum_mcp_fineBinning_preWidthScalingAndEvtNorm(TH1D* &H1D_jetPt_mcp, int iDataset, int iRadius, bool fcontrolMC, __attribute__ ((unused)) std::string options) {
  TH1D* H1D_jetPt_defaultBin;

  TString histogramName = "";
  if (doBkgSubtractionInMC) {
    histogramName = "h_jet_pt_part_rhoareasubtracted";
  } else {
    histogramName = "h_jet_pt_part";
  }

  if (!fcontrolMC){
    H1D_jetPt_defaultBin = (TH1D*)((TH1D*)file_O2Analysis_MCfileForMatrix->Get(analysisWorkflowMC+"/"+histogramName))->Clone("Get_Pt_spectrum_mcp_fineBinning"+Datasets[iDataset]+DatasetsNames[iDataset]);
  } else {
    H1D_jetPt_defaultBin = (TH1D*)((TH1D*)file_O2Analysis_ppSimDetectorEffect_unfoldingControl->Get(analysisWorkflowMC+"/"+histogramName))->Clone("Get_Pt_spectrum_mcp_fineBinning_unfoldingControl"+Datasets[iDataset]+DatasetsNames[iDataset]);
  }
  H1D_jetPt_defaultBin->Sumw2();
  cout << "test Get_Pt_spectrum_mcp_fineBinning_preWidthScalingAndEvtNorm 1" << endl;

  H1D_jetPt_mcp = (TH1D*)H1D_jetPt_defaultBin->Rebin(nBinPtJetsFine[iRadius],"jetPt_mcp_rebinned_fineBinning_"+Datasets[iDataset]+DatasetsNames[iDataset]+RadiusLegend[iRadius]+Datasets[iDataset]+DatasetsNames[iDataset], ptBinsJetsFine[iRadius]);
  cout << "test Get_Pt_spectrum_mcp_fineBinning_preWidthScalingAndEvtNorm 2" << endl;

  if (doWidthScalingEarly) {
    TransformRawHistToYield(H1D_jetPt_mcp);
  }
}

void Get_Pt_spectrum_mcp_recBinning_preWidthScalingAndEvtNorm(TH1D* &H1D_jetPt_mcp, int iDataset, int iRadius, bool fcontrolMC, __attribute__ ((unused)) std::string options) {
  TH1D* H1D_jetPt_defaultBin;

  TString histogramName = "";
  if (doBkgSubtractionInMC) {
    histogramName = "h_jet_pt_part_rhoareasubtracted";
  } else {
    histogramName = "h_jet_pt_part";
  }

  if (!fcontrolMC){
    H1D_jetPt_defaultBin = (TH1D*)((TH1D*)file_O2Analysis_MCfileForMatrix->Get(analysisWorkflowMC+"/"+histogramName))->Clone("Get_Pt_spectrum_mcp_recBinning"+Datasets[iDataset]+DatasetsNames[iDataset]);
  } else {
    H1D_jetPt_defaultBin = (TH1D*)((TH1D*)file_O2Analysis_ppSimDetectorEffect_unfoldingControl->Get(analysisWorkflowMC+"/"+histogramName))->Clone("Get_Pt_spectrum_mcp_recBinning_unfoldingControl"+Datasets[iDataset]+DatasetsNames[iDataset]);
  }
  H1D_jetPt_defaultBin->Sumw2();

  H1D_jetPt_mcp = (TH1D*)H1D_jetPt_defaultBin->Rebin(nBinPtJetsGen[iRadius],"jetPt_mcp_rebinned_recBinning"+RadiusLegend[iRadius]+Datasets[iDataset]+DatasetsNames[iDataset], ptBinsJetsGen[iRadius]);

  if (doWidthScalingEarly) {
    TransformRawHistToYield(H1D_jetPt_mcp);
  }
}

void Get_Pt_spectrum_mcd_fineBinning_preWidthScalingAndEvtNorm(TH1D* &H1D_jetPt_mcd, int iDataset, int iRadius, __attribute__ ((unused)) std::string options) {
  TH1D* H1D_jetPt_defaultBin;

  TString histogramName = "";
  if (doBkgSubtractionInMC) {
    histogramName = "h_jet_pt_rhoareasubtracted";
  } else {
    histogramName = "h_jet_pt";
  }

  H1D_jetPt_defaultBin = (TH1D*)((TH1D*)file_O2Analysis_MCfileForMatrix->Get(analysisWorkflowMC+"/"+histogramName))->Clone("Get_Pt_spectrum_mcd_fineBinning"+Datasets[iDataset]+DatasetsNames[iDataset]);
  H1D_jetPt_defaultBin->Sumw2();

  H1D_jetPt_mcd = (TH1D*)H1D_jetPt_defaultBin->Rebin(nBinPtJetsFine[iRadius],"jetPt_mcd_rebinned_fineBinning_"+Datasets[iDataset]+DatasetsNames[iDataset]+RadiusLegend[iRadius]+Datasets[iDataset]+DatasetsNames[iDataset], ptBinsJetsFine[iRadius]);

  if (doWidthScalingEarly) {
    TransformRawHistToYield(H1D_jetPt_mcd);
  }
}

void Get_Pt_spectrum_mcd_recBinning_preWidthScalingAndEvtNorm(TH1D* &H1D_jetPt_mcd, int iDataset, int iRadius, __attribute__ ((unused)) std::string options) {
  TH1D* H1D_jetPt_defaultBin;

  TString histogramName = "";
  if (doBkgSubtractionInMC) {
    histogramName = "h_jet_pt_rhoareasubtracted";
  } else {
    histogramName = "h_jet_pt";
  }

  H1D_jetPt_defaultBin = (TH1D*)((TH1D*)file_O2Analysis_MCfileForMatrix->Get(analysisWorkflowMC+"/"+histogramName))->Clone("Get_Pt_spectrum_mcd_recBinning"+Datasets[iDataset]+DatasetsNames[iDataset]);
  H1D_jetPt_defaultBin->Sumw2();

  H1D_jetPt_mcd = (TH1D*)H1D_jetPt_defaultBin->Rebin(nBinPtJetsRec[iRadius],"jetPt_mcd_rebinned_recBinning_"+RadiusLegend[iRadius]+Datasets[iDataset]+DatasetsNames[iDataset], ptBinsJetsRec[iRadius]);

  if (doWidthScalingEarly) {
    TransformRawHistToYield(H1D_jetPt_mcd);
  }
}

void Get_Pt_spectrum_mcd_genBinning_preWidthScalingAndEvtNorm(TH1D* &H1D_jetPt_mcd, int iDataset, int iRadius, __attribute__ ((unused)) std::string options) {
  TH1D* H1D_jetPt_defaultBin;

  TString histogramName = "";
  if (doBkgSubtractionInMC) {
    histogramName = "h_jet_pt_rhoareasubtracted";
  } else {
    histogramName = "h_jet_pt";
  }

    H1D_jetPt_defaultBin = (TH1D*)((TH1D*)file_O2Analysis_MCfileForMatrix->Get(analysisWorkflowMC+"/"+histogramName))->Clone("Get_Pt_spectrum_mcd_genBinning"+Datasets[iDataset]+DatasetsNames[iDataset]);
  H1D_jetPt_defaultBin->Sumw2();

  H1D_jetPt_mcd = (TH1D*)H1D_jetPt_defaultBin->Rebin(nBinPtJetsGen[iRadius],"jetPt_mcd_rebinned_genBinning_"+RadiusLegend[iRadius]+Datasets[iDataset]+DatasetsNames[iDataset], ptBinsJetsGen[iRadius]);

  if (doWidthScalingEarly) {
    TransformRawHistToYield(H1D_jetPt_mcd);
  }
}

void Get_Pt_spectrum_mcpMatched_genBinning_preWidthScalingAndEvtNorm(TH1D* &H1D_jetPt_mcpMatched, int iDataset, int iRadius, __attribute__ ((unused)) std::string options) {
  TH2D* H2D_jetPtMcdjetPtMcp;
  TH3D* H3D_jetRjetPtTagjetPtBase;

  TH1D* H1D_jetPt_mcpMatched_defaultBin;
  if (analysisWorkflowMC.Contains("jet-spectra-charged") == true) {
    if (doBkgSubtractionInMC) {
      H2D_jetPtMcdjetPtMcp = (TH2D*)((TH2D*)file_O2Analysis_MCfileForMatrix->Get(analysisWorkflowMC+"/h2_jet_pt_mcd_jet_pt_mcp_matchedgeo_rhoareasubtracted"))->Clone("Get_Pt_spectrum_mcpMatched_genBinning"+Datasets[iDataset]+DatasetsNames[iDataset]);
    } else {
      H2D_jetPtMcdjetPtMcp = (TH2D*)((TH2D*)file_O2Analysis_MCfileForMatrix->Get(analysisWorkflowMC+"/h2_jet_pt_mcd_jet_pt_mcp_matchedgeo"))->Clone("Get_Pt_spectrum_mcpMatched_genBinning"+Datasets[iDataset]+DatasetsNames[iDataset]);
    }
    H2D_jetPtMcdjetPtMcp->Sumw2();

    H1D_jetPt_mcpMatched_defaultBin = (TH1D*)H2D_jetPtMcdjetPtMcp->ProjectionY("jetPt_mcpMatched_genBinning_"+RadiusLegend[iRadius]+Datasets[iDataset]+DatasetsNames[iDataset], 1, H2D_jetPtMcdjetPtMcp->GetNbinsX(), "e");

  } else if (analysisWorkflowMC.Contains("jet-finder-charged-qa") == true) {
    H3D_jetRjetPtTagjetPtBase = (TH3D*)((TH3D*)file_O2Analysis_MCfileForMatrix->Get(analysisWorkflowMC+"/h3_jet_r_jet_pt_tag_jet_pt_base_matchedgeo"))->Clone("Get_Pt_spectrum_mcpMatched_genBinning"+Datasets[iDataset]+DatasetsNames[iDataset]);// base is mcd in jetfinderQA as of 06/2024, thus tag is mcp, and so hist is (x=r, y=mcp, z=mcd)
    H3D_jetRjetPtTagjetPtBase->Sumw2();

    int ibinJetRadius = H3D_jetRjetPtTagjetPtBase->GetXaxis()->FindBin(arrayRadius[iRadius]+GLOBAL_epsilon);
    H1D_jetPt_mcpMatched_defaultBin = (TH1D*)H3D_jetRjetPtTagjetPtBase->ProjectionY("jetPt_mcpMatched_genBinning_"+RadiusLegend[iRadius]+Datasets[iDataset]+DatasetsNames[iDataset], ibinJetRadius, ibinJetRadius, 1, H3D_jetRjetPtTagjetPtBase->GetNbinsZ(), "e");
  } else {
    cout << "analysisWorkflowMC name is wrong, need to be corrected" << endl;
    abort();
  }

  H1D_jetPt_mcpMatched = (TH1D*)H1D_jetPt_mcpMatched_defaultBin->Rebin(nBinPtJetsGen[iRadius],"jetPt_mcpMatched_genBinning_rebinned_"+RadiusLegend[iRadius]+Datasets[iDataset]+DatasetsNames[iDataset], ptBinsJetsGen[iRadius]);

  if (doWidthScalingEarly) {
    TransformRawHistToYield(H1D_jetPt_mcpMatched);
  }
}

void Get_Pt_spectrum_mcpMatched_fineBinning_preWidthScalingAndEvtNorm(TH1D* &H1D_jetPt_mcpMatched, int iDataset, int iRadius, __attribute__ ((unused)) std::string options) {
  TH2D* H2D_jetPtMcdjetPtMcp;
  TH3D* H3D_jetRjetPtTagjetPtBase;

  TH1D* H1D_jetPt_mcpMatched_defaultBin;

  if (analysisWorkflowMC.Contains("jet-spectra-charged") == true) {
    if (doBkgSubtractionInMC) {
      H2D_jetPtMcdjetPtMcp = (TH2D*)((TH2D*)file_O2Analysis_MCfileForMatrix->Get(analysisWorkflowMC+"/h2_jet_pt_mcd_jet_pt_mcp_matchedgeo_rhoareasubtracted"))->Clone("Get_Pt_spectrum_mcpMatched_fineBinning"+Datasets[iDataset]+DatasetsNames[iDataset]);
    } else {
      H2D_jetPtMcdjetPtMcp = (TH2D*)((TH2D*)file_O2Analysis_MCfileForMatrix->Get(analysisWorkflowMC+"/h2_jet_pt_mcd_jet_pt_mcp_matchedgeo"))->Clone("Get_Pt_spectrum_mcpMatched_fineBinning"+Datasets[iDataset]+DatasetsNames[iDataset]);
    }
    H2D_jetPtMcdjetPtMcp->Sumw2();

    H1D_jetPt_mcpMatched_defaultBin = (TH1D*)H2D_jetPtMcdjetPtMcp->ProjectionY("jetPt_mcpMatched_fineBinning_"+RadiusLegend[iRadius]+Datasets[iDataset]+DatasetsNames[iDataset], 1, H2D_jetPtMcdjetPtMcp->GetNbinsX(), "e");
  
  } else if (analysisWorkflowMC.Contains("jet-finder-charged-qa") == true) {
    H3D_jetRjetPtTagjetPtBase = (TH3D*)((TH3D*)file_O2Analysis_MCfileForMatrix->Get(analysisWorkflowMC+"/h3_jet_r_jet_pt_tag_jet_pt_base_matchedgeo"))->Clone("Get_Pt_spectrum_mcpMatched_fineBinning"+Datasets[iDataset]+DatasetsNames[iDataset]);// base is mcd in jetfinderQA as of 06/2024, thus tag is mcp, and so hist is (x=r, y=mcp, z=mcd)
    H3D_jetRjetPtTagjetPtBase->Sumw2();

    int ibinJetRadius = H3D_jetRjetPtTagjetPtBase->GetXaxis()->FindBin(arrayRadius[iRadius]+GLOBAL_epsilon);
    H1D_jetPt_mcpMatched_defaultBin = (TH1D*)H3D_jetRjetPtTagjetPtBase->ProjectionY("jetPt_mcpMatched_fineBinning_"+RadiusLegend[iRadius]+Datasets[iDataset]+DatasetsNames[iDataset], ibinJetRadius, ibinJetRadius, 1, H3D_jetRjetPtTagjetPtBase->GetNbinsZ(), "e");
  } else {
    cout << "analysisWorkflowMC name is wrong, need to be corrected" << endl;
    abort();
  }

  H1D_jetPt_mcpMatched = (TH1D*)H1D_jetPt_mcpMatched_defaultBin->Rebin(nBinPtJetsFine[iRadius],"jetPt_mcpMatched_fineBinning_rebinned"+RadiusLegend[iRadius]+Datasets[iDataset]+DatasetsNames[iDataset], ptBinsJetsFine[iRadius]);

  if (doWidthScalingEarly) {
    TransformRawHistToYield(H1D_jetPt_mcpMatched);
  }
}

void Get_Pt_spectrum_mcdMatched_genBinning_preWidthScalingAndEvtNorm(TH1D* &H1D_jetPt_mcdMatched, int iDataset, int iRadius, __attribute__ ((unused)) std::string options) {
  TH2D* H2D_jetPtMcdjetPtMcd;
  TH3D* H3D_jetRjetPtTagjetPtBase;

  TH1D* H1D_jetPt_mcdMatched_defaultBin;
  if (analysisWorkflowMC.Contains("jet-spectra-charged") == true) {
    if (doBkgSubtractionInMC) {
      H2D_jetPtMcdjetPtMcd = (TH2D*)((TH2D*)file_O2Analysis_MCfileForMatrix->Get(analysisWorkflowMC+"/h2_jet_pt_mcd_jet_pt_mcp_matchedgeo_rhoareasubtracted"))->Clone("Get_Pt_spectrum_mcdMatched_genBinning"+Datasets[iDataset]+DatasetsNames[iDataset]);
    } else {
      H2D_jetPtMcdjetPtMcd = (TH2D*)((TH2D*)file_O2Analysis_MCfileForMatrix->Get(analysisWorkflowMC+"/h2_jet_pt_mcd_jet_pt_mcp_matchedgeo"))->Clone("Get_Pt_spectrum_mcdMatched_genBinning"+Datasets[iDataset]+DatasetsNames[iDataset]);
    }
    H2D_jetPtMcdjetPtMcd->Sumw2();

    H1D_jetPt_mcdMatched_defaultBin = (TH1D*)H2D_jetPtMcdjetPtMcd->ProjectionX("jetPt_mcdMatched_genBinning_"+RadiusLegend[iRadius]+Datasets[iDataset]+DatasetsNames[iDataset], 1, H2D_jetPtMcdjetPtMcd->GetNbinsY(), "e");

  } else if (analysisWorkflowMC.Contains("jet-finder-charged-qa") == true) {
    H3D_jetRjetPtTagjetPtBase = (TH3D*)((TH3D*)file_O2Analysis_MCfileForMatrix->Get(analysisWorkflowMC+"/h3_jet_r_jet_pt_tag_jet_pt_base_matchedgeo"))->Clone("Get_Pt_spectrum_mcdMatched_genBinning"+Datasets[iDataset]+DatasetsNames[iDataset]);// base is mcd in jetfinderQA as of 06/2024, thus tag is mcp, and so hist is (x=r, y=mcp, z=mcd)
    H3D_jetRjetPtTagjetPtBase->Sumw2();

    int ibinJetRadius = H3D_jetRjetPtTagjetPtBase->GetXaxis()->FindBin(arrayRadius[iRadius]+GLOBAL_epsilon);
    H1D_jetPt_mcdMatched_defaultBin = (TH1D*)H3D_jetRjetPtTagjetPtBase->ProjectionZ("jetPt_mcdMatched_genBinning_"+RadiusLegend[iRadius]+Datasets[iDataset]+DatasetsNames[iDataset], ibinJetRadius, ibinJetRadius, 1, H3D_jetRjetPtTagjetPtBase->GetNbinsZ(), "e");
  } else {
    cout << "analysisWorkflowMC name is wrong, need to be corrected" << endl;
    abort();
  }

  H1D_jetPt_mcdMatched = (TH1D*)H1D_jetPt_mcdMatched_defaultBin->Rebin(nBinPtJetsGen[iRadius],"jetPt_mcdMatched_genBinning_rebinned_"+RadiusLegend[iRadius]+Datasets[iDataset]+DatasetsNames[iDataset], ptBinsJetsGen[iRadius]);

  if (doWidthScalingEarly) {
    TransformRawHistToYield(H1D_jetPt_mcdMatched);
  }
}

void Get_Pt_spectrum_mcdMatched_recBinning_preWidthScalingAndEvtNorm(TH1D* &H1D_jetPt_mcdMatched, int iDataset, int iRadius, __attribute__ ((unused)) std::string options) {
  TH2D* H2D_jetPtMcdjetPtMcd;
  TH3D* H3D_jetRjetPtTagjetPtBase;

  TH1D* H1D_jetPt_mcdMatched_defaultBin;

  if (analysisWorkflowMC.Contains("jet-spectra-charged") == true) {
    if (doBkgSubtractionInMC) {
      H2D_jetPtMcdjetPtMcd = (TH2D*)((TH2D*)file_O2Analysis_MCfileForMatrix->Get(analysisWorkflowMC+"/h2_jet_pt_mcd_jet_pt_mcp_matchedgeo_rhoareasubtracted"))->Clone("Get_Pt_spectrum_mcdMatched_recBinning"+Datasets[iDataset]+DatasetsNames[iDataset]);
    } else {
      H2D_jetPtMcdjetPtMcd = (TH2D*)((TH2D*)file_O2Analysis_MCfileForMatrix->Get(analysisWorkflowMC+"/h2_jet_pt_mcd_jet_pt_mcp_matchedgeo"))->Clone("Get_Pt_spectrum_mcdMatched_recBinning"+Datasets[iDataset]+DatasetsNames[iDataset]);
    }
    H2D_jetPtMcdjetPtMcd->Sumw2();

    H1D_jetPt_mcdMatched_defaultBin = (TH1D*)H2D_jetPtMcdjetPtMcd->ProjectionX("jetPt_mcdMatched_recBinning_"+RadiusLegend[iRadius]+Datasets[iDataset]+DatasetsNames[iDataset], 1, H2D_jetPtMcdjetPtMcd->GetNbinsY(), "e");
  
  } else if (analysisWorkflowMC.Contains("jet-finder-charged-qa") == true) {
    H3D_jetRjetPtTagjetPtBase = (TH3D*)((TH3D*)file_O2Analysis_MCfileForMatrix->Get(analysisWorkflowMC+"/h3_jet_r_jet_pt_tag_jet_pt_base_matchedgeo"))->Clone("Get_Pt_spectrum_mcdMatched_recBinning"+Datasets[iDataset]+DatasetsNames[iDataset]);// base is mcd in jetfinderQA as of 06/2024, thus tag is mcp, and so hist is (x=r, y=mcp, z=mcd)
    H3D_jetRjetPtTagjetPtBase->Sumw2();

    int ibinJetRadius = H3D_jetRjetPtTagjetPtBase->GetXaxis()->FindBin(arrayRadius[iRadius]+GLOBAL_epsilon);
    H1D_jetPt_mcdMatched_defaultBin = (TH1D*)H3D_jetRjetPtTagjetPtBase->ProjectionZ("jetPt_mcdMatched_recBinning_"+RadiusLegend[iRadius]+Datasets[iDataset]+DatasetsNames[iDataset], ibinJetRadius, ibinJetRadius, 1, H3D_jetRjetPtTagjetPtBase->GetNbinsZ(), "e");
  } else {
    cout << "analysisWorkflowMC name is wrong, need to be corrected" << endl;
    abort();
  }

  H1D_jetPt_mcdMatched = (TH1D*)H1D_jetPt_mcdMatched_defaultBin->Rebin(nBinPtJetsRec[iRadius],"jetPt_mcdMatched_recBinning_rebinned_"+RadiusLegend[iRadius]+Datasets[iDataset]+DatasetsNames[iDataset], ptBinsJetsRec[iRadius]);

  if (doWidthScalingEarly) {
    TransformRawHistToYield(H1D_jetPt_mcdMatched);
  }
}

void Get_Pt_spectrum_mcdMatched_fineBinning_preWidthScalingAndEvtNorm(TH1D* &H1D_jetPt_mcdMatched, int iDataset, int iRadius, __attribute__ ((unused)) std::string options) {
  TH2D* H2D_jetPtMcdjetPtMcd;
  TH3D* H3D_jetRjetPtTagjetPtBase;

  TH1D* H1D_jetPt_mcdMatched_defaultBin;

  if (analysisWorkflowMC.Contains("jet-spectra-charged") == true) {
    if (doBkgSubtractionInMC) {
      H2D_jetPtMcdjetPtMcd = (TH2D*)((TH2D*)file_O2Analysis_MCfileForMatrix->Get(analysisWorkflowMC+"/h2_jet_pt_mcd_jet_pt_mcp_matchedgeo_rhoareasubtracted"))->Clone("Get_Pt_spectrum_mcdMatched_fineBinning"+Datasets[iDataset]+DatasetsNames[iDataset]);
    } else {
      H2D_jetPtMcdjetPtMcd = (TH2D*)((TH2D*)file_O2Analysis_MCfileForMatrix->Get(analysisWorkflowMC+"/h2_jet_pt_mcd_jet_pt_mcp_matchedgeo"))->Clone("Get_Pt_spectrum_mcdMatched_fineBinning"+Datasets[iDataset]+DatasetsNames[iDataset]);
    }
    H2D_jetPtMcdjetPtMcd->Sumw2();

    H1D_jetPt_mcdMatched_defaultBin = (TH1D*)H2D_jetPtMcdjetPtMcd->ProjectionX("jetPt_mcdMatched_fineBinning_"+RadiusLegend[iRadius]+Datasets[iDataset]+DatasetsNames[iDataset], 1, H2D_jetPtMcdjetPtMcd->GetNbinsY(), "e");
  
  } else if (analysisWorkflowMC.Contains("jet-finder-charged-qa") == true) {
    H3D_jetRjetPtTagjetPtBase = (TH3D*)((TH3D*)file_O2Analysis_MCfileForMatrix->Get(analysisWorkflowMC+"/h3_jet_r_jet_pt_tag_jet_pt_base_matchedgeo"))->Clone("Get_Pt_spectrum_mcdMatched_fineBinning"+Datasets[iDataset]+DatasetsNames[iDataset]);// base is mcd in jetfinderQA as of 06/2024, thus tag is mcp, and so hist is (x=r, y=mcp, z=mcd)
    H3D_jetRjetPtTagjetPtBase->Sumw2();

    int ibinJetRadius = H3D_jetRjetPtTagjetPtBase->GetXaxis()->FindBin(arrayRadius[iRadius]+GLOBAL_epsilon);
    H1D_jetPt_mcdMatched_defaultBin = (TH1D*)H3D_jetRjetPtTagjetPtBase->ProjectionZ("jetPt_mcdMatched_fineBinning_"+RadiusLegend[iRadius]+Datasets[iDataset]+DatasetsNames[iDataset], ibinJetRadius, ibinJetRadius, 1, H3D_jetRjetPtTagjetPtBase->GetNbinsZ(), "e");
  } else {
    cout << "analysisWorkflowMC name is wrong, need to be corrected" << endl;
    abort();
  }

  H1D_jetPt_mcdMatched = (TH1D*)H1D_jetPt_mcdMatched_defaultBin->Rebin(nBinPtJetsFine[iRadius],"jetPt_mcdMatched_fineBinning_rebinned_"+RadiusLegend[iRadius]+Datasets[iDataset]+DatasetsNames[iDataset], ptBinsJetsFine[iRadius]);

  if (doWidthScalingEarly) {
    TransformRawHistToYield(H1D_jetPt_mcdMatched);
  }
}

void Get_Pt_spectrum_bkgCorrected_recBinning_preWidthScaling(TH1D* &H1D_jetPt, int iDataset, int iRadius, std::string options) {
  Get_Pt_spectrum_bkgCorrected_recBinning_preWidthScalingAndEvtNorm(H1D_jetPt, iDataset, iRadius, options);


  if (normaliseDistribsBeforeUnfolding) {
    if (!controlMC) {
      NormaliseRawHistToNEvents(H1D_jetPt, GetNEventsSelected_JetFramework(file_O2Analysis_list[iDataset], analysisWorkflowData));
    } else {
      NormaliseRawHistToNEvents(H1D_jetPt, GetNEventsSelected_JetFramework(file_O2Analysis_ppSimDetectorEffect_unfoldingControl, analysisWorkflow_unfoldingControl));
    }
  }
}
void Get_Pt_spectrum_bkgCorrected_recBinning(TH1D* &H1D_jetPt, int iDataset, int iRadius, std::string options) {
  Get_Pt_spectrum_bkgCorrected_recBinning_preWidthScaling(H1D_jetPt, iDataset, iRadius, options);

  if (doWidthScalingAtEnd) {
    TransformRawHistToYield(H1D_jetPt);
  }
}
void Get_Pt_spectrum_bkgCorrected_genBinning_preWidthScaling(TH1D* &H1D_jetPt, int iDataset, int iRadius, std::string options) {
  Get_Pt_spectrum_bkgCorrected_genBinning_preWidthScalingAndEvtNorm(H1D_jetPt, iDataset, iRadius, options);


  if (normaliseDistribsBeforeUnfolding) {
    if (!controlMC) {
      NormaliseRawHistToNEvents(H1D_jetPt, GetNEventsSelected_JetFramework(file_O2Analysis_list[iDataset], analysisWorkflowData));
    } else {
      NormaliseRawHistToNEvents(H1D_jetPt, GetNEventsSelected_JetFramework(file_O2Analysis_ppSimDetectorEffect_unfoldingControl, analysisWorkflow_unfoldingControl));
    }
  }
}

void Get_Pt_spectrum_bkgCorrected_genBinning(TH1D* &H1D_jetPt, int iDataset, int iRadius, std::string options) {
  Get_Pt_spectrum_bkgCorrected_genBinning_preWidthScaling(H1D_jetPt, iDataset, iRadius, options);

  if (doWidthScalingAtEnd) {
    TransformRawHistToYield(H1D_jetPt);
  }
}

void Get_Pt_spectrum_bkgCorrected_fineBinning_preWidthScaling(TH1D* &H1D_jetPt, int iDataset, int iRadius, std::string options) {
  Get_Pt_spectrum_bkgCorrected_fineBinning_preWidthScalingAndEvtNorm(H1D_jetPt, iDataset, iRadius, options);

  if (normaliseDistribsBeforeUnfolding) {
    if (!controlMC) {
      NormaliseRawHistToNEvents(H1D_jetPt, GetNEventsSelected_JetFramework(file_O2Analysis_list[iDataset], analysisWorkflowData));
    } else {
      NormaliseRawHistToNEvents(H1D_jetPt, GetNEventsSelected_JetFramework(file_O2Analysis_ppSimDetectorEffect_unfoldingControl, analysisWorkflow_unfoldingControl));
    }
  }
}
void Get_Pt_spectrum_bkgCorrected_fineBinning(TH1D* &H1D_jetPt, int iDataset, int iRadius, std::string options) {
  Get_Pt_spectrum_bkgCorrected_fineBinning_preWidthScaling(H1D_jetPt, iDataset, iRadius, options);

  if (doWidthScalingAtEnd) {
    TransformRawHistToYield(H1D_jetPt);
  }
}


void Get_Pt_spectrum_mcp_genBinning_preWidthScaling(TH1D* &H1D_jetPt_mcp, int iDataset, int iRadius, bool fcontrolMC, std::string options) {
  Get_Pt_spectrum_mcp_genBinning_preWidthScalingAndEvtNorm(H1D_jetPt_mcp, iDataset, iRadius, fcontrolMC, options);

  if (normaliseDistribsBeforeUnfolding) {
    if (mcIsWeighted) {
      if (!fcontrolMC) {
        NormaliseRawHistToNEvents(H1D_jetPt_mcp, GetNEventsSelected_JetFramework_weighted(file_O2Analysis_MCfileForMatrix, analysisWorkflowMC));
      } else {
        NormaliseRawHistToNEvents(H1D_jetPt_mcp, GetNEventsSelected_JetFramework_weighted(file_O2Analysis_ppSimDetectorEffect_unfoldingControl, analysisWorkflow_unfoldingControl));
      }
    } else {
      if (!fcontrolMC) {
        NormaliseRawHistToNEvents(H1D_jetPt_mcp, GetNEventsSelected_JetFramework(file_O2Analysis_MCfileForMatrix, analysisWorkflowMC));
      } else {
        NormaliseRawHistToNEvents(H1D_jetPt_mcp, GetNEventsSelected_JetFramework(file_O2Analysis_ppSimDetectorEffect_unfoldingControl, analysisWorkflow_unfoldingControl));
      }
    }
  }
}
void Get_Pt_spectrum_mcp_genBinning(TH1D* &H1D_jetPt_mcp, int iDataset, int iRadius, bool fcontrolMC, std::string options) {
  Get_Pt_spectrum_mcp_genBinning_preWidthScaling(H1D_jetPt_mcp, iDataset, iRadius, fcontrolMC, options);

  if (doWidthScalingAtEnd) {
    TransformRawHistToYield(H1D_jetPt_mcp);
  }
}

void Get_Pt_spectrum_mcp_fineBinning_preWidthScaling(TH1D* &H1D_jetPt_mcp, int iDataset, int iRadius, bool fcontrolMC, std::string options) {
  Get_Pt_spectrum_mcp_fineBinning_preWidthScalingAndEvtNorm(H1D_jetPt_mcp, iDataset, iRadius, fcontrolMC, options);

  if (normaliseDistribsBeforeUnfolding) {
    if (mcIsWeighted) {
      if (!fcontrolMC) {
        NormaliseRawHistToNEvents(H1D_jetPt_mcp, GetNEventsSelected_JetFramework_weighted(file_O2Analysis_MCfileForMatrix, analysisWorkflowMC));
      } else {
        NormaliseRawHistToNEvents(H1D_jetPt_mcp, GetNEventsSelected_JetFramework_weighted(file_O2Analysis_ppSimDetectorEffect_unfoldingControl, analysisWorkflow_unfoldingControl));
      }
    } else {
      if (!fcontrolMC) {
        NormaliseRawHistToNEvents(H1D_jetPt_mcp, GetNEventsSelected_JetFramework(file_O2Analysis_MCfileForMatrix, analysisWorkflowMC));
      } else {
        NormaliseRawHistToNEvents(H1D_jetPt_mcp, GetNEventsSelected_JetFramework(file_O2Analysis_ppSimDetectorEffect_unfoldingControl, analysisWorkflow_unfoldingControl));
      }
    }
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

  if (normaliseDistribsBeforeUnfolding) {
    if (mcIsWeighted) {
      NormaliseRawHistToNEvents(H1D_jetPt_mcd, GetNEventsSelected_JetFramework_weighted(file_O2Analysis_MCfileForMatrix, analysisWorkflowMC));
    } else {
      NormaliseRawHistToNEvents(H1D_jetPt_mcd, GetNEventsSelected_JetFramework(file_O2Analysis_MCfileForMatrix, analysisWorkflowMC));
    }
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

  if (normaliseDistribsBeforeUnfolding) {
    if (mcIsWeighted) {
      NormaliseRawHistToNEvents(H1D_jetPt_mcd, GetNEventsSelected_JetFramework_weighted(file_O2Analysis_MCfileForMatrix, analysisWorkflowMC));
    } else {
      NormaliseRawHistToNEvents(H1D_jetPt_mcd, GetNEventsSelected_JetFramework(file_O2Analysis_MCfileForMatrix, analysisWorkflowMC));
    }
  }
}
void Get_Pt_spectrum_mcd_recBinning(TH1D* &H1D_jetPt_mcd, int iDataset, int iRadius, std::string options) {
  Get_Pt_spectrum_mcd_recBinning_preWidthScaling(H1D_jetPt_mcd, iDataset, iRadius, options);

  if (doWidthScalingAtEnd) {
    TransformRawHistToYield(H1D_jetPt_mcd);
  }
}

void Get_Pt_spectrum_mcd_genBinning_preWidthScaling(TH1D* &H1D_jetPt_mcd, int iDataset, int iRadius, std::string options) {
  Get_Pt_spectrum_mcd_genBinning_preWidthScalingAndEvtNorm(H1D_jetPt_mcd, iDataset, iRadius, options);

  if (normaliseDistribsBeforeUnfolding) {
    if (mcIsWeighted) {
      NormaliseRawHistToNEvents(H1D_jetPt_mcd, GetNEventsSelected_JetFramework_weighted(file_O2Analysis_MCfileForMatrix, analysisWorkflowMC));
    } else {
      NormaliseRawHistToNEvents(H1D_jetPt_mcd, GetNEventsSelected_JetFramework(file_O2Analysis_MCfileForMatrix, analysisWorkflowMC));
    }
  }
}
void Get_Pt_spectrum_mcd_genBinning(TH1D* &H1D_jetPt_mcd, int iDataset, int iRadius, std::string options) {
  Get_Pt_spectrum_mcd_genBinning_preWidthScaling(H1D_jetPt_mcd, iDataset, iRadius, options);

  if (doWidthScalingAtEnd) {
    TransformRawHistToYield(H1D_jetPt_mcd);
  }
}

void Get_Pt_spectrum_mcdMatched_genBinning_preWidthScaling(TH1D* &H1D_jetPt_mcdMatched, int iDataset, int iRadius, std::string options) {
  Get_Pt_spectrum_mcdMatched_genBinning_preWidthScalingAndEvtNorm(H1D_jetPt_mcdMatched, iDataset, iRadius, options);

  if (normaliseDistribsBeforeUnfolding) {
    if (mcIsWeighted) {
      NormaliseRawHistToNEvents(H1D_jetPt_mcdMatched, GetNEventsSelected_JetFramework_weighted(file_O2Analysis_MCfileForMatrix, analysisWorkflowMC));
    } else {
      NormaliseRawHistToNEvents(H1D_jetPt_mcdMatched, GetNEventsSelected_JetFramework(file_O2Analysis_MCfileForMatrix, analysisWorkflowMC));
    }
  }
}
void Get_Pt_spectrum_mcdMatched_genBinning(TH1D* &H1D_jetPt_mcdMatched, int iDataset, int iRadius, std::string options) {
  Get_Pt_spectrum_mcdMatched_genBinning_preWidthScaling(H1D_jetPt_mcdMatched, iDataset, iRadius, options);

  if (doWidthScalingAtEnd) {
    TransformRawHistToYield(H1D_jetPt_mcdMatched);
  }
}

void Get_Pt_spectrum_mcdMatched_recBinning_preWidthScaling(TH1D* &H1D_jetPt_mcdMatched, int iDataset, int iRadius, std::string options) {
  Get_Pt_spectrum_mcdMatched_recBinning_preWidthScalingAndEvtNorm(H1D_jetPt_mcdMatched, iDataset, iRadius, options);

  if (normaliseDistribsBeforeUnfolding) {
    if (mcIsWeighted) {
      NormaliseRawHistToNEvents(H1D_jetPt_mcdMatched, GetNEventsSelected_JetFramework_weighted(file_O2Analysis_MCfileForMatrix, analysisWorkflowMC));
    } else {
      NormaliseRawHistToNEvents(H1D_jetPt_mcdMatched, GetNEventsSelected_JetFramework(file_O2Analysis_MCfileForMatrix, analysisWorkflowMC));
    }
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

  if (normaliseDistribsBeforeUnfolding) {
    if (mcIsWeighted) {
      NormaliseRawHistToNEvents(H1D_jetPt_mcdMatched, GetNEventsSelected_JetFramework_weighted(file_O2Analysis_MCfileForMatrix, analysisWorkflowMC));
    } else {
      NormaliseRawHistToNEvents(H1D_jetPt_mcdMatched, GetNEventsSelected_JetFramework(file_O2Analysis_MCfileForMatrix, analysisWorkflowMC));
    }
  }
}
void Get_Pt_spectrum_mcdMatched_fineBinning(TH1D* &H1D_jetPt_mcdMatched, int iDataset, int iRadius, std::string options) {
  Get_Pt_spectrum_mcdMatched_fineBinning_preWidthScaling(H1D_jetPt_mcdMatched, iDataset, iRadius, options);

  if (doWidthScalingAtEnd) {
    TransformRawHistToYield(H1D_jetPt_mcdMatched);
  }
}

void Get_Pt_spectrum_mcpMatched_genBinning_preWidthScaling(TH1D* &H1D_jetPt_mcpMatched, int iDataset, int iRadius, std::string options) {
  Get_Pt_spectrum_mcpMatched_genBinning_preWidthScalingAndEvtNorm(H1D_jetPt_mcpMatched, iDataset, iRadius, options);

  if (normaliseDistribsBeforeUnfolding) {
    if (mcIsWeighted) {
      NormaliseRawHistToNEvents(H1D_jetPt_mcpMatched, GetNEventsSelected_JetFramework_weighted(file_O2Analysis_MCfileForMatrix, analysisWorkflowMC));
    } else {
      NormaliseRawHistToNEvents(H1D_jetPt_mcpMatched, GetNEventsSelected_JetFramework(file_O2Analysis_MCfileForMatrix, analysisWorkflowMC));
    }
  }
}
void Get_Pt_spectrum_mcpMatched_genBinning(TH1D* &H1D_jetPt_mcpMatched, int iDataset, int iRadius, std::string options) {
  Get_Pt_spectrum_mcpMatched_genBinning_preWidthScaling(H1D_jetPt_mcpMatched, iDataset, iRadius, options);

  if (doWidthScalingAtEnd) {
    TransformRawHistToYield(H1D_jetPt_mcpMatched);
  }
}

void Get_Pt_spectrum_mcpMatched_fineBinning_preWidthScaling(TH1D* &H1D_jetPt_mcpMatched, int iDataset, int iRadius, std::string options) {
  Get_Pt_spectrum_mcpMatched_fineBinning_preWidthScalingAndEvtNorm(H1D_jetPt_mcpMatched, iDataset, iRadius, options);

  if (normaliseDistribsBeforeUnfolding) {
    if (mcIsWeighted) {
      NormaliseRawHistToNEvents(H1D_jetPt_mcpMatched, GetNEventsSelected_JetFramework_weighted(file_O2Analysis_MCfileForMatrix, analysisWorkflowMC));
    } else {
      NormaliseRawHistToNEvents(H1D_jetPt_mcpMatched, GetNEventsSelected_JetFramework(file_O2Analysis_MCfileForMatrix, analysisWorkflowMC));
    }
  }
}
void Get_Pt_spectrum_mcpMatched_fineBinning(TH1D* &H1D_jetPt_mcpMatched, int iDataset, int iRadius, std::string options) {
  Get_Pt_spectrum_mcpMatched_fineBinning_preWidthScaling(H1D_jetPt_mcpMatched, iDataset, iRadius, options);

  if (doWidthScalingAtEnd) {
    TransformRawHistToYield(H1D_jetPt_mcpMatched);
  }
}

void Get_Pt_spectrum_mcp_recBinning_preWidthScaling(TH1D* &H1D_jetPt_mcp, int iDataset, int iRadius, bool fcontrolMC, std::string options) {
  Get_Pt_spectrum_mcp_recBinning_preWidthScalingAndEvtNorm(H1D_jetPt_mcp, iDataset, iRadius, fcontrolMC, options);

  if (normaliseDistribsBeforeUnfolding) {
    if (mcIsWeighted) {
      if (!fcontrolMC) {
        NormaliseRawHistToNEvents(H1D_jetPt_mcp, GetNEventsSelected_JetFramework_weighted(file_O2Analysis_MCfileForMatrix, analysisWorkflowMC));
      } else {
        NormaliseRawHistToNEvents(H1D_jetPt_mcp, GetNEventsSelected_JetFramework_weighted(file_O2Analysis_ppSimDetectorEffect_unfoldingControl, analysisWorkflow_unfoldingControl));
      }
    } else {
      if (!fcontrolMC) {
        NormaliseRawHistToNEvents(H1D_jetPt_mcp, GetNEventsSelected_JetFramework(file_O2Analysis_MCfileForMatrix, analysisWorkflowMC));
      } else {
        NormaliseRawHistToNEvents(H1D_jetPt_mcp, GetNEventsSelected_JetFramework(file_O2Analysis_ppSimDetectorEffect_unfoldingControl, analysisWorkflow_unfoldingControl));
      }
    }
  }
}
void Get_Pt_spectrum_mcp_recBinning(TH1D* &H1D_jetPt_mcp, int iDataset, int iRadius, bool fcontrolMC, std::string options) {
  Get_Pt_spectrum_mcp_recBinning_preWidthScaling(H1D_jetPt_mcp, iDataset, iRadius, fcontrolMC, options);

  if (doWidthScalingAtEnd) {
    TransformRawHistToYield(H1D_jetPt_mcp);
  }
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////// Response matrix functions ///////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void Get_PtResponseMatrix_DetectorAndFluctuationsCombined_fineBinning(TH2D* &H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined_fineBinning, TH2D* H2D_jetPtResponseMatrix_detectorResponse, TH2D* H2D_jetPtResponseMatrix_fluctuations, int iDataset, int iRadius, __attribute__ ((unused)) std::string options) {
  // https://github.com/alisw/AliPhysics/blob/master/PWGJE/PWGJE/AliAnaChargedJetResponseMaker.cxx for ann example that works, by marta verveij

  TString partialUniqueSpecifier = Datasets[iDataset]+"_R="+Form("%.1f",arrayRadius[iRadius]);

  // matrix product of fluct response times det response; assumes the two are of the same size binning wise
  // Careful: xy of hist and ij of Resp(i,j) are inverted ! hist(j,i) = matrix(i,j) and so if matrix(i,j)=SUM(matrixA(i,k)matrixB(k,j)) then hist(j,i)=SUM(histA(k,i)histB(j,k)), and if we replace j,i by gen,rec we get hist(gen,rec)=SUM(histA(k,rec)histB(gen,k))
  H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined_fineBinning = (TH2D*)GetMatrixProductTH2xTH2(H2D_jetPtResponseMatrix_fluctuations, H2D_jetPtResponseMatrix_detectorResponse).Clone("Get_PtResponseMatrix_DetectorAndFluctuationsCombined_fineBinning"+partialUniqueSpecifier);

  if (drawIntermediateResponseMatrices) {
    TString* pdfName_preRebin = new TString("responseMatrix_combined_preRebinning"+jetType[iJetType]+"_"+Datasets[iDataset]+DatasetsNames[iDataset]);
    TString* pdfName_preRebin_logz = new TString("responseMatrix_combined_preRebinning"+(TString)"_R="+Form("%.1f",arrayRadius[iRadius])+"_"+Datasets[iDataset]+DatasetsNames[iDataset]+"_logz");

    TString textContext_preRebin(contextCustomOneField(*texDatasetsComparisonCommonDenominator, ""));

    Draw_TH2_Histogram(H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined_fineBinning, textContext_preRebin, pdfName_preRebin, texPtJetRecX, texPtJetGenX, texCollisionDataInfo, drawnWindowAuto, th2ContoursNone, contourNumberNone, "");
    Draw_TH2_Histogram(H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined_fineBinning, textContext_preRebin, pdfName_preRebin_logz, texPtJetRecX, texPtJetGenX, texCollisionDataInfo, drawnWindowAuto, th2ContoursNone, contourNumberNone, "logz");
  }
  // cout << "bin(topleft 1,N) = " << H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined_fineBinning->GetBinContent(1,H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined_fineBinning->GetNbinsY()) << endl;
  // cout << "bin(bottom left 1,1) = " << H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined_fineBinning->GetBinContent(1,1) << endl;
}


void Get_PtResponseMatrix_DetectorAndFluctuationsCombined(TH2D* &H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined, TH2D* H2D_jetPtResponseMatrix_detectorResponse, TH2D* H2D_jetPtResponseMatrix_fluctuations, int iDataset, int iRadius, std::string options) {
  // https://github.com/alisw/AliPhysics/blob/master/PWGJE/PWGJE/AliAnaChargedJetResponseMaker.cxx for ann example that works, by marta verveij

  TString partialUniqueSpecifier = Datasets[iDataset]+"_R="+Form("%.1f",arrayRadius[iRadius]);
  cout << "test Get_PtResponseMatrix_DetectorAndFluctuationsCombined 0" << endl;

  TH2D* H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined_preRebin;
  Get_PtResponseMatrix_DetectorAndFluctuationsCombined_fineBinning(H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined_preRebin, H2D_jetPtResponseMatrix_detectorResponse, H2D_jetPtResponseMatrix_fluctuations, iDataset, iRadius, options);
  cout << "test Get_PtResponseMatrix_DetectorAndFluctuationsCombined 1" << endl;

  // Merge bins of the combined response matrix with fine binning to get the coarse one
  TH1D* priorSpectrumMerging;
  bool debugBool = false;
  Get_Pt_spectrum_mcp_fineBinning(priorSpectrumMerging, iDataset, iRadius, false, options); //take mcp as prior by default
  cout << "test Get_PtResponseMatrix_DetectorAndFluctuationsCombined 2" << endl;
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
    Get_Pt_spectrum_bkgCorrected_fineBinning(priorSpectrumMerging, iDataset, iRadius, options);
    H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined = (TH2D*)RebinVariableBins2D_PriorWeightedBinMerging(H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined_preRebin, nBinPtJetsRec[iRadius], nBinPtJetsGen[iRadius], ptBinsJetsRec[iRadius], ptBinsJetsGen[iRadius], priorSpectrumMerging, debugBool).Clone("Get_PtResponseMatrix_DetectorAndFluctuationsCombined"+partialUniqueSpecifier);
  }
  if (options.find("noMergingPrior") != std::string::npos) {
    cout << "test Get_PtResponseMatrix_DetectorAndFluctuationsCombined 2.1" << endl;
    H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined = (TH2D*)RebinVariableBins2D(H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined_preRebin, nBinPtJetsRec[iRadius], nBinPtJetsGen[iRadius], ptBinsJetsRec[iRadius], ptBinsJetsGen[iRadius], debugBool).Clone("Get_PtResponseMatrix_DetectorAndFluctuationsCombined"+partialUniqueSpecifier);
    cout << "test Get_PtResponseMatrix_DetectorAndFluctuationsCombined 2.2" << endl;
  }
  if (options.find("testAliPhysics") != std::string::npos) {
    H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined = (TH2D*)RebinVariableBins2D_aliPhysics(H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined_preRebin, nBinPtJetsRec[iRadius], nBinPtJetsGen[iRadius], ptBinsJetsRec[iRadius], ptBinsJetsGen[iRadius], debugBool)->Clone("Get_PtResponseMatrix_DetectorAndFluctuationsCombined"+partialUniqueSpecifier);
  }
  // normalising priorSpectrum with evtNorm doesn't change anything as the weighting does prior_content(i)/prior_integral()
  // dividing priorSpectrum by binwidth doesn't change anything for the same reason
  cout << "test Get_PtResponseMatrix_DetectorAndFluctuationsCombined 3" << endl;


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

    TString* pdfName = new TString("responseMatrix_combined_preNorm"+jetType[iJetType]+"_"+Datasets[iDataset]+DatasetsNames[iDataset]);
    TString* pdfName_logz = new TString("responseMatrix_combined_preNorm"+(TString)"_R="+Form("%.1f",arrayRadius[iRadius])+"_"+Datasets[iDataset]+DatasetsNames[iDataset]+"_logz");

    TString textContext(contextCustomOneField(*texDatasetsComparisonCommonDenominator, ""));

    Draw_TH2_Histogram(H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined_preNorm, textContext, pdfName, texPtJetRecX, texPtJetGenX, texCollisionDataInfo, drawnWindowAuto, th2ContoursNone, contourNumberNone, "");
    Draw_TH2_Histogram(H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined_preNorm, textContext, pdfName_logz, texPtJetRecX, texPtJetGenX, texCollisionDataInfo, drawnWindowAuto, th2ContoursNone, contourNumberNone, "logz");
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

    TString* pdfName = new TString("responseMatrix_combined_preWeighting"+jetType[iJetType]+"_"+Datasets[iDataset]+DatasetsNames[iDataset]);
    TString* pdfName_logz = new TString("responseMatrix_combined_preWeighting"+(TString)"_R="+Form("%.1f",arrayRadius[iRadius])+"_"+Datasets[iDataset]+DatasetsNames[iDataset]+"_logz");

    TString textContext(contextCustomOneField(*texDatasetsComparisonCommonDenominator, ""));

    Draw_TH2_Histogram(H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined_preWeighting, textContext, pdfName, texPtJetRecX, texPtJetGenX, texCollisionDataInfo, drawnWindowAuto, th2ContoursNone, contourNumberNone, "");
    Draw_TH2_Histogram(H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined_preWeighting, textContext, pdfName_logz, texPtJetRecX, texPtJetGenX, texCollisionDataInfo, drawnWindowAuto, th2ContoursNone, contourNumberNone, "logz");
  }
}

void ReweightResponseMatrixWithPrior(TH2D* &H2D_jetPtResponseMatrix, int iDataset, int iRadius, std::string options) {
  
  // before this, all y-slices (ie pt gen slices) have been normalised to 1;means each pt gen slice has a proba of 1
  // withthis function, we give each slice a weight so that they have different normalisation values, corresponding to the prior 

  // prior choice; none by default (flat)
  TH1D* priorSpectrumWeighting;
  if (options.find("mcpPriorUnfolding") != std::string::npos) {
    if (!normGenAndMeasByNEvts) {
      Get_Pt_spectrum_mcp_genBinning_preWidthScalingAndEvtNorm(priorSpectrumWeighting, iDataset, iRadius, false, options); 
    } else {
      Get_Pt_spectrum_mcp_genBinning_preWidthScaling(priorSpectrumWeighting, iDataset, iRadius, false, options); 
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
      Get_Pt_spectrum_mcd_genBinning_preWidthScalingAndEvtNorm(priorSpectrumWeighting, iDataset, iRadius, options); 
    } else {
      Get_Pt_spectrum_mcd_genBinning_preWidthScaling(priorSpectrumWeighting, iDataset, iRadius, options); 
    }
    WeightMatrixWithPrior(H2D_jetPtResponseMatrix, priorSpectrumWeighting);
  }
  if (options.find("measuredPriorUnfolding") != std::string::npos) {
    if (!normGenAndMeasByNEvts) {
      Get_Pt_spectrum_bkgCorrected_genBinning_preWidthScalingAndEvtNorm(priorSpectrumWeighting, iDataset, iRadius, options);
    } else {
      Get_Pt_spectrum_bkgCorrected_genBinning_preWidthScaling(priorSpectrumWeighting, iDataset, iRadius, options); 
    }
    WeightMatrixWithPrior(H2D_jetPtResponseMatrix, priorSpectrumWeighting);
  }
  if (options.find("testAliPhysics") != std::string::npos) {
    if (!normGenAndMeasByNEvts) {
      Get_Pt_spectrum_bkgCorrected_genBinning_preWidthScalingAndEvtNorm(priorSpectrumWeighting, iDataset, iRadius, options); 
    } else {
      Get_Pt_spectrum_bkgCorrected_genBinning_preWidthScaling(priorSpectrumWeighting, iDataset, iRadius, options); 
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


void ReweightResponseMatrixWithPrior_fineBinning(TH2D* &H2D_jetPtResponseMatrix, int iDataset, int iRadius, std::string options) {
  
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
      Get_Pt_spectrum_bkgCorrected_fineBinning_preWidthScalingAndEvtNorm(priorSpectrumWeighting, iDataset, iRadius, options); 
    } else {
      Get_Pt_spectrum_bkgCorrected_fineBinning_preWidthScaling(priorSpectrumWeighting, iDataset, iRadius, options); 
    }
    WeightMatrixWithPrior(H2D_jetPtResponseMatrix, priorSpectrumWeighting);
  }
  if (options.find("testAliPhysics") != std::string::npos) {
    if (!normGenAndMeasByNEvts) {
      Get_Pt_spectrum_bkgCorrected_fineBinning_preWidthScalingAndEvtNorm(priorSpectrumWeighting, iDataset, iRadius, options); 
    } else {
      Get_Pt_spectrum_bkgCorrected_fineBinning_preWidthScaling(priorSpectrumWeighting, iDataset, iRadius, options); 
    }
    H2D_jetPtResponseMatrix = (TH2D*)NormalizeResponsMatrixYaxisWithPrior(H2D_jetPtResponseMatrix, priorSpectrumWeighting)->Clone(H2D_jetPtResponseMatrix->GetName()+(TString)"_testAliPhysics");
  }
}

void Get_PtResponseMatrix_detectorResponse(TH2D* &H2D_jetPtResponseMatrix_detectorResponse, int iDataset, int iRadius) {
  TH3D* H3D_jetRpartPtdetPt;
  TH2D* H2D_jetPtMcdjetPtMcd;

  TH2D* H2D_gen_det_geoMatched;
  TH2D* H2D_gen_det_geoMatched_rebinned;
  cout << "Get_PtResponseMatrix_detectorResponse 1" << endl;
  TString partialUniqueSpecifier = Datasets[iDataset]+"_R="+Form("%.1f",arrayRadius[iRadius]);
  if (analysisWorkflowMC.Contains("jet-finder-charged-qa") == true) {
    H3D_jetRpartPtdetPt = (TH3D*)((TH3D*)file_O2Analysis_MCfileForMatrix->Get(analysisWorkflowMC+"/h3_jet_r_jet_pt_tag_jet_pt_base_matchedgeo"))->Clone("Get_PtResponseMatrix_detectorResponse"+partialUniqueSpecifier);// base is mcd in jetfinderQA as of 06/2024, thus tag is mcp, and so hist is (x=r, y=mcp, z=mcd)
    H3D_jetRpartPtdetPt->Sumw2();

    int ibinJetRadius = H3D_jetRpartPtdetPt->GetXaxis()->FindBin(arrayRadius[iRadius]+GLOBAL_epsilon);
    H3D_jetRpartPtdetPt->GetXaxis()->SetRange(ibinJetRadius, ibinJetRadius);
    // project H3D onto a H2D, option "yz" means y goes on y-axis while z goes on x-axis, and so H2D_gen_det_geoMatched will be (x=mcd, y=mcp)
    H2D_gen_det_geoMatched = (TH2D*)H3D_jetRpartPtdetPt->Project3D(partialUniqueSpecifier+"_genrec_e_yz"); //can't use letter D in this or it seems to replace the histogram in current pad (see documentation of ProjectionX function. Isn't mentioned in project3D sadly)
  } else if (analysisWorkflowMC.Contains("jet-spectra-charged") == true) {
    if (doBkgSubtractionInMC) {
      H2D_jetPtMcdjetPtMcd = (TH2D*)((TH2D*)file_O2Analysis_MCfileForMatrix->Get(analysisWorkflowMC+"/h2_jet_pt_mcd_jet_pt_mcp_matchedgeo_rhoareasubtracted"))->Clone("Get_PtResponseMatrix_detectorResponse"+Datasets[iDataset]+DatasetsNames[iDataset]);
    } else {
      H2D_jetPtMcdjetPtMcd = (TH2D*)((TH2D*)file_O2Analysis_MCfileForMatrix->Get(analysisWorkflowMC+"/h2_jet_pt_mcd_jet_pt_mcp_matchedgeo"))->Clone("Get_PtResponseMatrix_detectorResponse"+Datasets[iDataset]+DatasetsNames[iDataset]);
    }
    H2D_jetPtMcdjetPtMcd->Sumw2();

    // H2D_gen_det_geoMatched = (TH2D*)GetTransposeHistogram(H2D_jetPtMcdjetPtMcd).Clone(partialUniqueSpecifier+"_genrec");
    H2D_gen_det_geoMatched = (TH2D*)H2D_jetPtMcdjetPtMcd->Clone(partialUniqueSpecifier+"_genrec");
  }

  // keep (gen, gen) for the bins; rec will be introduced in the fluctuation response, and by multiplication will stay in the combined matrix
  TH2D* H2D_response = (TH2D*)RebinVariableBins2D(H2D_gen_det_geoMatched, nBinPtJetsFine[iRadius], nBinPtJetsFine[iRadius], ptBinsJetsFine[iRadius], ptBinsJetsFine[iRadius]).Clone("Get_PtResponseMatrix_detectorResponse_rebinned"+partialUniqueSpecifier);

  if (useYSliceNorm) {
    NormaliseYSlicesToOne(H2D_response);
  }
  if (normDetRespByNEvts) {
    if (mcIsWeighted) {
      H2D_response->Scale(1./GetNEventsSelected_JetFramework_weighted(file_O2Analysis_MCfileForMatrix, analysisWorkflowMC));
    } else {
      int Nevents = GetNEventsSelected_JetFramework(file_O2Analysis_MCfileForMatrix, analysisWorkflowMC);
      for (int iBinX = 0; iBinX < H2D_response->GetNbinsX(); iBinX++) {
        for (int iBinY = 0; iBinY < H2D_response->GetNbinsY(); iBinY++) {
          H2D_response->SetBinContent(iBinX, iBinY, H2D_response->GetBinContent(iBinX, iBinY) * 1./Nevents);
          H2D_response->SetBinError(iBinX, iBinY, H2D_response->GetBinError(iBinX, iBinY) * 1./Nevents);
        }
      }
      // H2D_response->Scale(1./GetNEventsSelected_JetFramework(file_O2Analysis_MCfileForMatrix, analysisWorkflowMC));
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

void Get_PtResponseMatrix_Fluctuations(TH2D* &H2D_jetPtResponseMatrix_fluctuations, int iDataset, int iRadius) { 
  // see Hiroki Yokoyama thesis
  // iRadius is for chosing the pT binning
  TString partialUniqueSpecifier = Datasets[iDataset]+"_R="+Form("%.1f",arrayRadius[iRadius]);

  if (useFactorisedMatrix == false){
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

    H2D_fluctuations_centrality = (TH2D*)((TH2D*)file_O2Analysis_list[iDataset]->Get(analysisWorkflowBkg+"/h2_centrality_rhorandomcone"+randomConeTypeList[randomConeType]))->Clone("Get_PtResponseMatrix_Fluctuations"+Datasets[iDataset]+DatasetsNames[iDataset]);
    H2D_fluctuations_centrality->Sumw2();


    int ibinCent_low = H2D_fluctuations_centrality->GetXaxis()->FindBin(centralityRange[0]);
    int ibinCent_high = H2D_fluctuations_centrality->GetXaxis()->FindBin(centralityRange[1]);
    H1D_fluctuations = (TH1D*)H2D_fluctuations_centrality->ProjectionY("bkgFluctuationCentrality_highRes_"+partialUniqueSpecifier, ibinCent_low, ibinCent_high, "e");

    NormaliseRawHistToIntegral(H1D_fluctuations); // normalising fluctuations to 1

    TH2D H2D_response = TH2D("H2D_response_"+partialUniqueSpecifier, "H2D_response_"+partialUniqueSpecifier, nBinPtJetsFine[iRadius], ptBinsJetsFine[iRadius], nBinPtJetsFine[iRadius], ptBinsJetsFine[iRadius]); // actually doesn't work if original histogram has fixed bin size

    //==================== Build response matrix: shift deltaPt by pT gen along the pT rec axis ====================//
    int ibinZeroFluct= H1D_fluctuations->FindBin(0+GLOBAL_epsilon);
    double integralError;
    for(int iBinRec = 0; iBinRec <= H2D_response.GetNbinsX()+1; iBinRec++){
      for(int iBinGen = 0; iBinGen <= H2D_response.GetNbinsY()+1; iBinGen++){
        double ptGen = H2D_response.GetYaxis()->GetBinLowEdge(iBinGen); // was bincenter before but then it'd give .5 values of GeV, and 
        double ptRec_low = H2D_response.GetXaxis()->GetBinLowEdge(iBinRec);
        double ptRec_up = H2D_response.GetXaxis()->GetBinLowEdge(iBinRec+1);
        int iBin_fluct_low = H1D_fluctuations->GetXaxis()->FindBin(ptRec_low - ptGen + GLOBAL_epsilon);
        int iBin_fluct_high = H1D_fluctuations->GetXaxis()->FindBin(ptRec_up - ptGen - GLOBAL_epsilon);
        H2D_response.SetBinContent(iBinRec, iBinGen, H1D_fluctuations->IntegralAndError(iBin_fluct_low, iBin_fluct_high, integralError)); 
        H2D_response.SetBinError(iBinRec, iBinGen, integralError); 
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




bool  Get_Pt_JetEfficiency(TH1D* &H1D_jetEfficiency, int iDataset, int iRadius, std::string options){
  TH1D* H1D_jetPt_mcp;
  TH1D* H1D_jetPt_mcpMatched;
  bool divideSuccess;

  Get_Pt_spectrum_mcp_genBinning(H1D_jetPt_mcp, iDataset, iRadius, false, options);
  Get_Pt_spectrum_mcpMatched_genBinning(H1D_jetPt_mcpMatched, iDataset, iRadius, options);

  H1D_jetEfficiency = (TH1D*)H1D_jetPt_mcpMatched->Clone("H1D_jetEfficiency"+Datasets[iDataset]+DatasetsNames[iDataset]+"_R="+Form("%.1f", arrayRadius[iRadius]));
  divideSuccess = H1D_jetEfficiency->Divide(H1D_jetEfficiency, H1D_jetPt_mcp, 1., 1., "b");
  if (!divideSuccess){
    cout << "################## Get_Pt_JetEfficiency FAILED!!!!! ##################" << endl;
  }
  return divideSuccess;
}
bool  Get_Pt_JetEfficiency_fineBinning(TH1D* &H1D_jetEfficiency, int iDataset, int iRadius, std::string options){
  TH1D* H1D_jetPt_mcp;
  TH1D* H1D_jetPt_mcpMatched;
  bool divideSuccess;

  Get_Pt_spectrum_mcp_fineBinning(H1D_jetPt_mcp, iDataset, iRadius, false, options);
  Get_Pt_spectrum_mcpMatched_fineBinning(H1D_jetPt_mcpMatched, iDataset, iRadius, options);

  H1D_jetEfficiency = (TH1D*)H1D_jetPt_mcpMatched->Clone("H1D_jetEfficiency_fineBinning"+Datasets[iDataset]+DatasetsNames[iDataset]+"_R="+Form("%.1f", arrayRadius[iRadius]));
  divideSuccess = H1D_jetEfficiency->Divide(H1D_jetEfficiency, H1D_jetPt_mcp, 1., 1., "b");
  if (!divideSuccess){
    cout << "################## Get_Pt_JetEfficiency FAILED!!!!! ##################" << endl;
  }
  return divideSuccess;
}

bool Get_Pt_JetFakes(TH1D* &H1D_jetFakes, int iDataset, int iRadius, std::string options){
  TH1D* H1D_jetPt_mcd;
  TH1D* H1D_jetPt_mcdMatched;
  bool divideSuccess;

    cout << "test Get_Pt_JetFakes 1" << endl;

  Get_Pt_spectrum_mcd_recBinning(H1D_jetPt_mcd, iDataset, iRadius, options);
    cout << "test Get_Pt_JetFakes 2" << endl;
  Get_Pt_spectrum_mcdMatched_recBinning(H1D_jetPt_mcdMatched, iDataset, iRadius, options);
    cout << "test Get_Pt_JetFakes 3" << endl;


  H1D_jetFakes = (TH1D*)H1D_jetPt_mcdMatched->Clone("H1D_jetFakes"+Datasets[iDataset]+DatasetsNames[iDataset]+"_R="+Form("%.1f", arrayRadius[iRadius]));
  divideSuccess = H1D_jetFakes->Divide(H1D_jetPt_mcd);
  if (!divideSuccess){
    cout << "################## Get_Pt_JetFakes FAILED!!!!! ##################" << endl;
  }

  return divideSuccess;
}
bool Get_Pt_JetFakes_fineBinning(TH1D* &H1D_jetFakes, int iDataset, int iRadius, std::string options){
  TH1D* H1D_jetPt_mcd;
  TH1D* H1D_jetPt_mcdMatched;
  bool divideSuccess;

  Get_Pt_spectrum_mcd_fineBinning(H1D_jetPt_mcd, iDataset, iRadius, options);
  Get_Pt_spectrum_mcdMatched_fineBinning(H1D_jetPt_mcdMatched, iDataset, iRadius, options);

  H1D_jetFakes = (TH1D*)H1D_jetPt_mcdMatched->Clone("H1D_jetFakes_fineBinning"+Datasets[iDataset]+DatasetsNames[iDataset]+"_R="+Form("%.1f", arrayRadius[iRadius]));
  divideSuccess = H1D_jetFakes->Divide(H1D_jetPt_mcd);
  if (!divideSuccess){
    cout << "################## Get_Pt_JetFakes FAILED!!!!! ##################" << endl;
  }
  return divideSuccess;
}





//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////// RooUnfold Custom Utilities ///////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



int GetSvdBestRegularisationParameter_notYetSatisfying(TSVDUnfold* unfoldTSvd){
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


std::pair<int, RooUnfold*> Get_Pt_spectrum_unfolded_preWidthScalingAndEvtNorm(TH1D* &H1D_jetPt_unfolded, int iDataset, int iRadius, int unfoldParameterInput, std::string options) {
  //for now makes the assumption gen and rec have the same pT binning

  TString partialUniqueSpecifier = Datasets[iDataset]+"_R="+Form("%.1f",arrayRadius[iRadius]);

  TH1D* measured;
  TH1D* mcp;
  TH1D* mcd;

  if (!normGenAndMeasByNEvts) {
    Get_Pt_spectrum_bkgCorrected_recBinning_preWidthScalingAndEvtNorm(measured, iDataset, iRadius, options);
    Get_Pt_spectrum_mcp_genBinning_preWidthScalingAndEvtNorm(mcp, iDataset, iRadius, false, options);
    Get_Pt_spectrum_mcd_recBinning_preWidthScalingAndEvtNorm(mcd, iDataset, iRadius, options);

    if (useFineBinningTest) {
      Get_Pt_spectrum_bkgCorrected_fineBinning_preWidthScalingAndEvtNorm(measured, iDataset, iRadius, options);
      Get_Pt_spectrum_mcp_fineBinning_preWidthScalingAndEvtNorm(mcp, iDataset, iRadius, false, options);
      Get_Pt_spectrum_mcd_fineBinning_preWidthScalingAndEvtNorm(mcd, iDataset, iRadius, options);
    }
  } else{
    Get_Pt_spectrum_bkgCorrected_recBinning_preWidthScaling(measured, iDataset, iRadius, options);
    Get_Pt_spectrum_mcp_genBinning_preWidthScaling(mcp, iDataset, iRadius, false, options);
    Get_Pt_spectrum_mcd_recBinning_preWidthScaling(mcd, iDataset, iRadius, options);

    if (useFineBinningTest) {
      Get_Pt_spectrum_bkgCorrected_fineBinning_preWidthScaling(measured, iDataset, iRadius, options);
      Get_Pt_spectrum_mcp_fineBinning_preWidthScaling(mcp, iDataset, iRadius, false, options);
      Get_Pt_spectrum_mcd_fineBinning_preWidthScaling(mcd, iDataset, iRadius, options);
    }
  } 

  bool divideSuccessFakes;
  TH1D* H1D_jetFakes;
  if (useManualRespMatrixSettingMethod) {
    if (applyFakes) {
      if (!useFineBinningTest){ 
        divideSuccessFakes = Get_Pt_JetFakes(H1D_jetFakes, iDataset, iRadius, options);
      } else {
        divideSuccessFakes = Get_Pt_JetFakes_fineBinning(H1D_jetFakes, iDataset, iRadius, options);
      }
      if (divideSuccessFakes){
        measured->Multiply(H1D_jetFakes);
      } else {
        cout << "################## measured->Multiply(H1D_jetFakes) failed because Get_Pt_JetFakes() FAILED!!!!! ##################" << endl;
      }
    }
  }


  // TH1D* measured;
  // Get_Pt_spectrum_bkgCorrected_recBinning(measured, iDataset, iRadius, options);
  // TH1D* mcp;
  // Get_Pt_spectrum_mcp_genBinning(mcp, iDataset, iRadius, options);
  // TH1D* mcdMatched;
  // Get_Pt_spectrum_mcdMatched_genBinning(mcdMatched, iDataset, iRadius, options);

  TH1D* H1D_kinematicEfficiency;
  TH2D* H2D_jetPtResponseMatrix_fluctuations;
  TH2D* H2D_jetPtResponseMatrix_detectorResponse;
  TH2D* H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined;
  TH2D* H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined_fineBinning;
  
  Get_PtResponseMatrix_Fluctuations(H2D_jetPtResponseMatrix_fluctuations, iDataset, iRadius);
  Get_PtResponseMatrix_detectorResponse(H2D_jetPtResponseMatrix_detectorResponse, iDataset, iRadius);
  // compute matrixFluctuations times matrixDetector

  Get_PtResponseMatrix_DetectorAndFluctuationsCombined_fineBinning(H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined_fineBinning, H2D_jetPtResponseMatrix_detectorResponse, H2D_jetPtResponseMatrix_fluctuations, iDataset, iRadius, options);
  Get_ResponseMatrix_Pt_KinematicEffiency(H1D_kinematicEfficiency, H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined_fineBinning, partialUniqueSpecifier, iRadius); // I want the efficiency before the reweighting and normalisation
  if (useFineBinningTest) {
    ReweightResponseMatrixWithPrior_fineBinning(H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined_fineBinning, iDataset, iRadius, options);
  }
  
  Get_PtResponseMatrix_DetectorAndFluctuationsCombined(H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined, H2D_jetPtResponseMatrix_detectorResponse, H2D_jetPtResponseMatrix_fluctuations, iDataset, iRadius, options);
  ReweightResponseMatrixWithPrior(H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined, iDataset, iRadius, options);

  if (drawIntermediateResponseMatrices) {
    TH2D* H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined_postWeighting = (TH2D*)H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined->Clone("H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined_postWeighting"+partialUniqueSpecifier);

    TString* pdfNamePost = new TString("responseMatrix_combined_postWeighting"+jetType[iJetType]+"_"+Datasets[iDataset]+DatasetsNames[iDataset]);
    TString* pdfNamePost_logz = new TString("responseMatrix_combined_postWeighting"+(TString)"_R="+Form("%.1f",arrayRadius[iRadius])+"_"+Datasets[iDataset]+DatasetsNames[iDataset]+"_logz");

    TString textContextPost(contextCustomOneField(*texDatasetsComparisonCommonDenominator, ""));

    Draw_TH2_Histogram(H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined_postWeighting, textContextPost, pdfNamePost, texPtJetRecX, texPtJetGenX, texCollisionDataInfo, drawnWindowAuto, th2ContoursNone, contourNumberNone, "");
    Draw_TH2_Histogram(H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined_postWeighting, textContextPost, pdfNamePost_logz, texPtJetRecX, texPtJetGenX, texCollisionDataInfo, drawnWindowAuto, th2ContoursNone, contourNumberNone, "logz");
  }

  cout << "RooUnfoldResponse setting - start" << endl;
  RooUnfoldResponse* response;
  if (useManualRespMatrixSettingMethod) {

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
      if (mcIsWeighted) {
        fake->Scale(1./GetNEventsSelected_JetFramework_weighted(file_O2Analysis_MCfileForMatrix, analysisWorkflowMC));
      } else {
        fake->Scale(1./GetNEventsSelected_JetFramework(file_O2Analysis_MCfileForMatrix, analysisWorkflowMC));
      }
    }
    for (auto i = 1; i <= fake->GetNbinsX(); i++) { 
      cout << "Fake i = " << i << " content: "<< fake->GetBinContent(i) << " - " << mcdMatched->GetBinContent(i) << " = " << fake->GetBinContent(i) - mcdMatched->GetBinContent(i) << " ------ error of mcd = " << fake->GetBinError(i) << ", error of mcdMatched = " << mcdMatched->GetBinError(i) << endl;
    }
    fake->Add(mcdMatched, -1);
    TH1D *miss = (TH1D*)mcp->Clone(); // for each bin of gen-level jet distribution, how many are not matched to a rec-level jet
    if (normDetRespByNEvts) {
      if (mcIsWeighted) {
        miss->Scale(1./GetNEventsSelected_JetFramework_weighted(file_O2Analysis_MCfileForMatrix, analysisWorkflowMC));
      } else {
        miss->Scale(1./GetNEventsSelected_JetFramework(file_O2Analysis_MCfileForMatrix, analysisWorkflowMC));
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
  RooUnfoldBayes* unfoldBayes = new RooUnfoldBayes(response, measured, unfoldParameterInput);
  int unfoldParameterSvdInitial = 1;
  RooUnfoldSvd* unfoldSvdInitialiser = new RooUnfoldSvd(response, measured, unfoldParameterSvdInitial); // instance of RooUnfoldSvd only used to find the best regularisation parameter; the unfolded spectrum returned by it is not retrived
  RooUnfold* unfold = unfoldBayes; // default Bayes
  TH1D* hist_unfold;
  cout << "Get_Pt_spectrum_unfolded_preWidthScalingAndEvtNorm test - 1" << endl;

  int unfoldParameter;
  if (options.find("Svd") != std::string::npos) {
    unfoldSvdInitialiser->Hreco(); // necessary to have GetD() give a meaningful output
    TSVDUnfold *tsvdUnfold = (TSVDUnfold*)unfoldSvdInitialiser->Impl();
    if (automaticBestSvdParameter) {
      unfoldParameter = GetSvdBestRegularisationParameter_notYetSatisfying(tsvdUnfold);
    } else {
      unfoldParameter = unfoldParameterInput;
    }
    RooUnfoldSvd* unfoldSvd = new RooUnfoldSvd(response, measured, unfoldParameter); // the RooUnfoldSvd instance that is actually used to unfold, with the best regularisation parameter

    // unfoldSvd->SetRegParm(-200);  // doesnt seem to change anything, even if I put 0 here ... annoying; but accoring to mailing list of roounfold it doesn't work well
    hist_unfold = (TH1D*)(unfoldSvd->Hreco());
    unfold = unfoldSvd;

    // plot svd d distribution
    TH1D* H1D_D = tsvdUnfold->GetD();
    TString* pdfName_regparam = new TString("Svd_regularisationd_distribution_"+(TString)"_R="+Form("%.1f",arrayRadius[iRadius])+"_"+Datasets[iDataset]+DatasetsNames[iDataset]);
    TString textContext(contextCustomOneField(*texDatasetsComparisonCommonDenominator, ""));
    std::array<std::array<float, 2>, 2> drawnWindowSvdParam = {{{0, 30}, {0.01, 10000}}}; // {{xmin, xmax}, {ymin, ymax}}
    Draw_TH1_Histogram(H1D_D, textContext, pdfName_regparam, texSvdDvector, texSvdK, texCollisionDataInfo, drawnWindowSvdParam, legendPlacementAuto, contextPlacementAuto, "logy");

  } else if (options.find("Bayes") != std::string::npos) {
    unfoldParameter = unfoldParameterInput;
    hist_unfold = (TH1D*)(unfoldBayes->Hreco());
    unfold = unfoldBayes;
  }  

  cout << "Get_Pt_spectrum_unfolded_preWidthScalingAndEvtNorm test - 2" << endl;

  H1D_jetPt_unfolded = (TH1D*)hist_unfold->Clone("H1D_jetPt_unfolded"+partialUniqueSpecifier);

  bool divideSuccessEff;
  TH1D* H1D_jetEfficiency;
  if (!useFineBinningTest){ 
    divideSuccessEff = Get_Pt_JetEfficiency(H1D_jetEfficiency, iDataset, iRadius, options);
  } else {
    divideSuccessEff = Get_Pt_JetEfficiency_fineBinning(H1D_jetEfficiency, iDataset, iRadius, options);
  }
  cout << "Get_Pt_spectrum_unfolded_preWidthScalingAndEvtNorm test - 3" << endl;

  if (useManualRespMatrixSettingMethod) {
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
  cout << "Get_Pt_spectrum_unfolded_preWidthScalingAndEvtNorm test - 4" << endl;

  if (drawIntermediateResponseMatrices) {
    TH2D* H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined_postUnfolding = (TH2D*)unfold->response()->Hresponse()->Clone("H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined_postUnfolding"+partialUniqueSpecifier);

    TString* pdfName = new TString("responseMatrix_combined_postUnfolding"+jetType[iJetType]+"_"+Datasets[iDataset]+DatasetsNames[iDataset]);
    TString* pdfName_logz = new TString("responseMatrix_combined_postUnfolding"+(TString)"_R="+Form("%.1f",arrayRadius[iRadius])+"_"+Datasets[iDataset]+DatasetsNames[iDataset]+"_logz");

    TString textContext(contextCustomOneField(*texDatasetsComparisonCommonDenominator, ""));

    Draw_TH2_Histogram(H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined_postUnfolding, textContext, pdfName, texPtJetRecX, texPtJetGenX, texCollisionDataInfo, drawnWindowAuto, th2ContoursNone, contourNumberNone, "");
    Draw_TH2_Histogram(H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined_postUnfolding, textContext, pdfName_logz, texPtJetRecX, texPtJetGenX, texCollisionDataInfo, drawnWindowAuto, th2ContoursNone, contourNumberNone, "logz");
  }

  cout << "Get_Pt_spectrum_unfolded_preWidthScalingAndEvtNorm test - 5" << endl;

  // if (doWidthScalingEarly) {
  //   TransformRawHistToYield(H1D_jetPt_unfolded);
  // }

  std::pair<int, RooUnfold*> unfoldInfo(unfoldParameter, unfold);
  return unfoldInfo;
}


std::pair<int, RooUnfold*> Get_Pt_spectrum_unfolded_preWidthScaling(TH1D* &H1D_jetPt_unfolded, int iDataset, int iRadius, int unfoldParameterInput, std::string options) {
  std::pair<int, RooUnfold*> unfoldInfo = Get_Pt_spectrum_unfolded_preWidthScalingAndEvtNorm(H1D_jetPt_unfolded, iDataset, iRadius, unfoldParameterInput, options);
  cout << "Get_Pt_spectrum_unfolded_preWidthScaling test 1" << endl;
  if (normaliseDistribsAfterUnfolding){
    if (!controlMC) {
      NormaliseRawHistToNEvents(H1D_jetPt_unfolded, GetNEventsSelected_JetFramework(file_O2Analysis_list[iDataset], analysisWorkflowData));
    } else {
      if (!mcIsWeighted) {
        NormaliseRawHistToNEvents(H1D_jetPt_unfolded, GetNEventsSelected_JetFramework_weighted(file_O2Analysis_ppSimDetectorEffect_unfoldingControl, analysisWorkflow_unfoldingControl));
      } else {
        NormaliseRawHistToNEvents(H1D_jetPt_unfolded, GetNEventsSelected_JetFramework(file_O2Analysis_ppSimDetectorEffect_unfoldingControl, analysisWorkflow_unfoldingControl));
      }
    }
  }
   cout << "Get_Pt_spectrum_unfolded_preWidthScaling test 2" << endl;

  return unfoldInfo;
}
std::pair<int, RooUnfold*> Get_Pt_spectrum_unfolded(TH1D* &H1D_jetPt_unfolded, int iDataset, int iRadius, int unfoldParameterInput, std::string options) {
  std::pair<int, RooUnfold*> unfoldInfo = Get_Pt_spectrum_unfolded_preWidthScaling(H1D_jetPt_unfolded, iDataset, iRadius, unfoldParameterInput, options);

  if (doWidthScalingAtEnd) {
    TransformRawHistToYield(H1D_jetPt_unfolded);
  }
  return unfoldInfo;
}



void Get_Pt_spectrum_unfoldedThenRefolded_preWidthScalingAndEvtNorm(TH1D* &H1D_jetPt_unfoldedThenRefolded, int iDataset, int iRadius, int unfoldParameterInput, std::string options) {
  TH1D* H1D_jetPt_unfolded;
  // TH1D* H1D_jetPt_raw[nRadius];
  TString partialUniqueSpecifier = Datasets[iDataset]+"_R="+Form("%.1f",arrayRadius[iRadius]);


  RooUnfold* unfold = Get_Pt_spectrum_unfolded_preWidthScalingAndEvtNorm(H1D_jetPt_unfolded, iDataset, iRadius, unfoldParameterInput, options).second;
  // Get_Pt_spectrum_unfolded(H1D_jetPt_unfolded, iDataset, iRadius, unfoldParameterInput, options);

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

  Get_PtResponseMatrix_Fluctuations(H2D_jetPtResponseMatrix_fluctuations, iDataset, iRadius);
  Get_PtResponseMatrix_detectorResponse(H2D_jetPtResponseMatrix_detectorResponse, iDataset, iRadius);

  Get_PtResponseMatrix_DetectorAndFluctuationsCombined_fineBinning(H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined_fineBinning, H2D_jetPtResponseMatrix_detectorResponse, H2D_jetPtResponseMatrix_fluctuations, iDataset, iRadius, options);
  Get_ResponseMatrix_Pt_KinematicEffiency(H1D_kinematicEfficiency, H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined_fineBinning, partialUniqueSpecifier, iRadius); // I want the efficiency before the reweighting and normalisation
  if (useFineBinningTest) {
    ReweightResponseMatrixWithPrior_fineBinning(H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined_fineBinning, iDataset, iRadius, options);
  }

  Get_PtResponseMatrix_DetectorAndFluctuationsCombined(H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined, H2D_jetPtResponseMatrix_detectorResponse, H2D_jetPtResponseMatrix_fluctuations, iDataset, iRadius, options);
  ReweightResponseMatrixWithPrior(H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined, iDataset, iRadius, options);
  

  if (drawIntermediateResponseMatrices) {
    TH2D* H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined_duringRefolding = (TH2D*)H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined->Clone("H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined_duringRefolding"+partialUniqueSpecifier);

    TString* pdfName = new TString("responseMatrix_combined_duringRefolding"+jetType[iJetType]+"_"+Datasets[iDataset]+DatasetsNames[iDataset]);
    TString* pdfName_logz = new TString("responseMatrix_combined_duringRefolding"+(TString)"_R="+Form("%.1f",arrayRadius[iRadius])+"_"+Datasets[iDataset]+DatasetsNames[iDataset]+"_logz");

    TString textContext(contextCustomOneField(*texDatasetsComparisonCommonDenominator, ""));

    Draw_TH2_Histogram(H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined_duringRefolding, textContext, pdfName, texPtJetRecX, texPtJetGenX, texCollisionDataInfo, drawnWindowAuto, th2ContoursNone, contourNumberNone, "");
    Draw_TH2_Histogram(H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined_duringRefolding, textContext, pdfName_logz, texPtJetRecX, texPtJetGenX, texCollisionDataInfo, drawnWindowAuto, th2ContoursNone, contourNumberNone, "logz");
  }

  cout << "((((((((((((((((()))))))))))))))))" << endl;
  cout << "REFOLDING TEST: pre efficiency" << endl;
  for(int iBinX = 0; iBinX <= H1D_jetPt_unfolded->GetNbinsX()+1; iBinX++){
    cout << "H1D_jetPt_unfolded(" << iBinX << ") = " << H1D_jetPt_unfolded->GetBinContent(iBinX) << endl;
  }
  

  bool divideSuccessEff;
  TH1D* H1D_jetEfficiency;
  if (!useFineBinningTest){ 
    divideSuccessEff = Get_Pt_JetEfficiency(H1D_jetEfficiency, iDataset, iRadius, options);
  } else {
    divideSuccessEff = Get_Pt_JetEfficiency_fineBinning(H1D_jetEfficiency, iDataset, iRadius, options);
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
  if (useManualRespMatrixSettingMethod) {
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
  if (applyFakes) { // for useManualRespMatrixSettingMethod set to false, it will give a slightly different result as Get_Pt_JetFakes doesn't get fakes from the same histogram as the one used to encode fakes in the response matrix
    if (!useFineBinningTest){ 
      divideSuccessFakes = Get_Pt_JetFakes(H1D_jetFakes, iDataset, iRadius, options);
    } else {
      divideSuccessFakes = Get_Pt_JetFakes_fineBinning(H1D_jetFakes, iDataset, iRadius, options);
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

void Get_Pt_spectrum_unfoldedThenRefolded_preWidthScaling(TH1D* &H1D_jetPt_unfoldedThenRefolded, int iDataset, int iRadius, int unfoldParameterInput, std::string options) {
  Get_Pt_spectrum_unfoldedThenRefolded_preWidthScalingAndEvtNorm(H1D_jetPt_unfoldedThenRefolded, iDataset, iRadius, unfoldParameterInput, options);

  if (normaliseDistribsAfterUnfolding){
    if (!controlMC) {
      NormaliseRawHistToNEvents(H1D_jetPt_unfoldedThenRefolded, GetNEventsSelected_JetFramework(file_O2Analysis_list[iDataset], analysisWorkflowData));
    } else {
      if (!mcIsWeighted) {
        NormaliseRawHistToNEvents(H1D_jetPt_unfoldedThenRefolded, GetNEventsSelected_JetFramework_weighted(file_O2Analysis_ppSimDetectorEffect_unfoldingControl, analysisWorkflow_unfoldingControl));
      } else {
        NormaliseRawHistToNEvents(H1D_jetPt_unfoldedThenRefolded, GetNEventsSelected_JetFramework(file_O2Analysis_ppSimDetectorEffect_unfoldingControl, analysisWorkflow_unfoldingControl));
      }
    }
  }
}
void Get_Pt_spectrum_unfoldedThenRefolded(TH1D* &H1D_jetPt_unfoldedThenRefolded, int iDataset, int iRadius, int unfoldParameterInput, std::string options) {
  Get_Pt_spectrum_unfoldedThenRefolded_preWidthScaling(H1D_jetPt_unfoldedThenRefolded, iDataset, iRadius, unfoldParameterInput, options);

  if (doWidthScalingAtEnd) {
    TransformRawHistToYield(H1D_jetPt_unfoldedThenRefolded);
  }
}


void Get_Pt_spectrum_unfoldedThenRefolded_RooUnfoldMethod_preWidthScalingAndEvtNorm(TH1D* &H1D_jetPt_unfoldedThenRefolded, int iDataset, int iRadius, int unfoldParameterInput, std::string options) {
  // Matches exactly with the manual method IF NO PRIOR
  // if I have a non flat prior, then the roounfold method gives me a good closure test, but not the manual method!

  TH1D* H1D_jetPt_unfolded;
  // TH1D* H1D_jetPt_raw[nRadius];
  TString partialUniqueSpecifier = Datasets[iDataset]+"_R="+Form("%.1f",arrayRadius[iRadius]);


  RooUnfold* unfold = Get_Pt_spectrum_unfolded_preWidthScalingAndEvtNorm(H1D_jetPt_unfolded, iDataset, iRadius, unfoldParameterInput, options).second;

  TH2D* H2D_jetPtResponseMatrix_fluctuations;
  TH2D* H2D_jetPtResponseMatrix_detectorResponse;
  TH2D* H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined;
  TH2D* H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined_fineBinning;
  TH1D* H1D_kinematicEfficiency;

  Get_PtResponseMatrix_Fluctuations(H2D_jetPtResponseMatrix_fluctuations, iDataset, iRadius);
  Get_PtResponseMatrix_detectorResponse(H2D_jetPtResponseMatrix_detectorResponse, iDataset, iRadius);

  Get_PtResponseMatrix_DetectorAndFluctuationsCombined_fineBinning(H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined_fineBinning, H2D_jetPtResponseMatrix_detectorResponse, H2D_jetPtResponseMatrix_fluctuations, iDataset, iRadius, options);
  Get_ResponseMatrix_Pt_KinematicEffiency(H1D_kinematicEfficiency, H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined_fineBinning, partialUniqueSpecifier, iRadius); // I want the efficiency before the reweighting and normalisation

  Get_PtResponseMatrix_DetectorAndFluctuationsCombined(H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined, H2D_jetPtResponseMatrix_detectorResponse, H2D_jetPtResponseMatrix_fluctuations, iDataset, iRadius, options);
  ReweightResponseMatrixWithPrior(H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined, iDataset, iRadius, options);

  bool divideSuccessEff;
  TH1D* H1D_jetEfficiency;
  if (!useFineBinningTest){ 
    divideSuccessEff = Get_Pt_JetEfficiency(H1D_jetEfficiency, iDataset, iRadius, options);
  } else {
    divideSuccessEff = Get_Pt_JetEfficiency_fineBinning(H1D_jetEfficiency, iDataset, iRadius, options);
  }

  if (useManualRespMatrixSettingMethod) { // this condition isn't present in the manual refolding method as in the manual case the response matrix object is just a th2 and not a roounfoldresponse object and doesn't have the efficiencies/fakes encoded in it
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
  if (useManualRespMatrixSettingMethod) { // this condition isn't present in the manual refolding method as in the manual case the response matrix object is just a th2 and not a roounfoldresponse object and doesn't have the efficiencies/fakes encoded in it
    if (applyFakes) {
      if (!useFineBinningTest){ 
        divideSuccessFakes = Get_Pt_JetFakes(H1D_jetFakes, iDataset, iRadius, options);
      } else {
        divideSuccessFakes = Get_Pt_JetFakes_fineBinning(H1D_jetFakes, iDataset, iRadius, options);
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

void Get_Pt_spectrum_unfoldedThenRefolded_RooUnfoldMethod_preWidthScaling(TH1D* &H1D_jetPt_unfoldedThenRefolded, int iDataset, int iRadius, int unfoldParameterInput, std::string options) {
  Get_Pt_spectrum_unfoldedThenRefolded_RooUnfoldMethod_preWidthScalingAndEvtNorm(H1D_jetPt_unfoldedThenRefolded, iDataset, iRadius, unfoldParameterInput, options);

  if (normaliseDistribsAfterUnfolding){
    if (!controlMC) {
      NormaliseRawHistToNEvents(H1D_jetPt_unfoldedThenRefolded, GetNEventsSelected_JetFramework(file_O2Analysis_list[iDataset], analysisWorkflowData));
    } else {
      if (!mcIsWeighted) {
        NormaliseRawHistToNEvents(H1D_jetPt_unfoldedThenRefolded, GetNEventsSelected_JetFramework_weighted(file_O2Analysis_ppSimDetectorEffect_unfoldingControl, analysisWorkflow_unfoldingControl));
      } else {
        NormaliseRawHistToNEvents(H1D_jetPt_unfoldedThenRefolded, GetNEventsSelected_JetFramework(file_O2Analysis_ppSimDetectorEffect_unfoldingControl, analysisWorkflow_unfoldingControl));
      }
    }
  }
}
void Get_Pt_spectrum_unfoldedThenRefolded_RooUnfoldMethod(TH1D* &H1D_jetPt_unfoldedThenRefolded, int iDataset, int iRadius, int unfoldParameterInput, std::string options) {
  Get_Pt_spectrum_unfoldedThenRefolded_RooUnfoldMethod_preWidthScaling(H1D_jetPt_unfoldedThenRefolded, iDataset, iRadius, unfoldParameterInput, options);

  if (doWidthScalingAtEnd) {
    TransformRawHistToYield(H1D_jetPt_unfoldedThenRefolded);
  }
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

  Draw_TH1_Histogram(H1D_jetPt_raw, textContext, pdfName, texPtX, yAxisLabel, texCollisionDataInfo, drawnWindow, legendPlacementAuto, contextPlacementAuto, "logy");
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

  Draw_TH1_Histograms(H1D_jetPt_mcp_collection, genVsRecBinningLegend, 2, textContext, pdfName, texPtX, yAxisLabel, texCollisionDataInfo, drawnWindow, legendPlacementAuto, contextPlacementAuto, "logy"); 
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

  Draw_TH1_Histogram(H1D_jetPt_mcdMatched, textContext, pdfName, texPtX, yAxisLabel, texCollisionDataInfo, drawnWindow, legendPlacementAuto, contextPlacementAuto, "logy");
}


void Draw_Pt_efficiency_jets(int iDataset, int iRadius, std::string options) {
  TH1D* H1D_jetEfficiency;
  bool divideSuccess;

  divideSuccess = Get_Pt_JetEfficiency(H1D_jetEfficiency, iDataset, iRadius, options);
  
  TString* pdfName = new TString("jet_"+jetType[iJetType]+"_"+jetLevel[iJetLevel]+"_"+Datasets[iDataset]+DatasetsNames[iDataset]+"_R="+Form("%.1f", arrayRadius[iRadius])+"_Pt_efficiency");

  TString textContext(contextCustomOneField(*texDatasetsComparisonCommonDenominator, ""));

  if (divideSuccess){
    Draw_TH1_Histogram(H1D_jetEfficiency, textContext, pdfName, texPtJetGenX, texJetEfficiency, texCollisionDataInfo, drawnWindow, legendPlacementAuto, contextPlacementAuto, "efficiency");
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

  TString* pdfName = new TString("kinematicEfficiency_"+partialUniqueSpecifier);

  TString textContext(contextCustomOneField(*texDatasetsComparisonCommonDenominator, ""));

  Draw_TH1_Histogram(H1D_kinematicEfficiency, textContext, pdfName, texPtJetGenX, texJetKinematicEfficiency, texCollisionDataInfo, drawnWindow, legendPlacementAuto, contextPlacementAuto, "");
}


void Draw_FakeRatio(int iDataset, int iRadius, std::string options) {
  TH1D* H1D_fakeRatio;

  TString partialUniqueSpecifier = Datasets[iDataset]+"_R="+Form("%.1f",arrayRadius[iRadius]);
  TString name_H1D_fakeRatio = Datasets[iDataset]+"_R="+Form("%.1f",arrayRadius[iRadius]);


  Get_Pt_JetFakes(H1D_fakeRatio, iDataset, iRadius, options);

  TString* pdfName = new TString("fakeRatio_"+partialUniqueSpecifier);

  TString textContext(contextCustomOneField(*texDatasetsComparisonCommonDenominator, ""));

  Draw_TH1_Histogram(H1D_fakeRatio, textContext, pdfName, texPtJetRecX, texFakeRatio, texCollisionDataInfo, drawnWindow, legendPlacementAuto, contextPlacementAuto, "");
}

void Draw_ResponseMatrices_Fluctuations(int iDataset, int iRadius) {

  TH2D* H2D_jetPtResponseMatrix_fluctuations;

  Get_PtResponseMatrix_Fluctuations(H2D_jetPtResponseMatrix_fluctuations, iDataset, iRadius);

  TString* pdfName_logz = new TString("responseMatrix_fluctuationsBackground_"+(TString)"_R="+Form("%.1f",arrayRadius[iRadius])+"_"+Datasets[iDataset]+DatasetsNames[iDataset]+"_logz");
  TString* pdfNameFullRes_logz = new TString("responseMatrix_fluctuationsBackground_"+(TString)"_R="+Form("%.1f",arrayRadius[iRadius])+"_"+Datasets[iDataset]+DatasetsNames[iDataset]+"FullRes_logz");
  TString* pdfName = new TString("responseMatrix_fluctuationsBackground_"+(TString)"_R="+Form("%.1f",arrayRadius[iRadius])+"_"+Datasets[iDataset]+DatasetsNames[iDataset]);
  TString* pdfNameFullRes = new TString("responseMatrix_fluctuationsBackground_"+(TString)"_R="+Form("%.1f",arrayRadius[iRadius])+"_"+Datasets[iDataset]+DatasetsNames[iDataset]+"_FullRes");

  TString textContext(contextCustomOneField(*texDatasetsComparisonCommonDenominator, ""));


  double th2ContourCustom[1] = {0.000001}; // hardcoded at 10-6 for now
  int contourNumberCustom = 1;

  Draw_TH2_Histogram(H2D_jetPtResponseMatrix_fluctuations, textContext, pdfName_logz, texPtJetBkgCorrX, texPtJetBkgFreeX, texCollisionDataInfo, drawnWindowAuto, th2ContourCustom, contourNumberCustom, "logz");
  Draw_TH2_Histogram(H2D_jetPtResponseMatrix_fluctuations, textContext, pdfName, texPtJetBkgCorrX, texPtJetBkgFreeX, texCollisionDataInfo, drawnWindowAuto, th2ContoursNone, contourNumberNone, "");
}

void Draw_ResponseMatrices_detectorResponse(int iDataset, int iRadius) {

  TH2D* H2D_jetPtResponseMatrix_detectorResponse;
  cout << "Draw_ResponseMatrices_detectorResponse 1" << endl;
  Get_PtResponseMatrix_detectorResponse(H2D_jetPtResponseMatrix_detectorResponse, iDataset, iRadius);
  cout << "Draw_ResponseMatrices_detectorResponse 2" << endl;

  TString* pdfName = new TString("responseMatrix_detectorEffects_"+jetType[iJetType]+"_"+Datasets[iDataset]+DatasetsNames[iDataset]);
  TString* pdfName_logz = new TString("responseMatrix_detectorEffects_"+(TString)"_R="+Form("%.1f",arrayRadius[iRadius])+"_"+Datasets[iDataset]+DatasetsNames[iDataset]+"_logz");

  TString textContext(contextCustomOneField(*texDatasetsComparisonCommonDenominator, ""));

  Draw_TH2_Histogram(H2D_jetPtResponseMatrix_detectorResponse, textContext, pdfName, texPtJetRecX, texPtJetGenX, texCollisionDataInfo, drawnWindowAuto, th2ContoursNone, contourNumberNone, "");
  Draw_TH2_Histogram(H2D_jetPtResponseMatrix_detectorResponse, textContext, pdfName_logz, texPtJetRecX, texPtJetGenX, texCollisionDataInfo, drawnWindowAuto, th2ContoursNone, contourNumberNone, "logz");
}

void Draw_ResponseMatrices_DetectorAndFluctuationsCombined(int iDataset, int iRadius, std::string options) {

  TH2D* H2D_jetPtResponseMatrix_fluctuations;
  TH2D* H2D_jetPtResponseMatrix_detectorResponse;
  TH2D* H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined;


  Get_PtResponseMatrix_Fluctuations(H2D_jetPtResponseMatrix_fluctuations, iDataset, iRadius);
  Get_PtResponseMatrix_detectorResponse(H2D_jetPtResponseMatrix_detectorResponse, iDataset, iRadius);
  Get_PtResponseMatrix_DetectorAndFluctuationsCombined(H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined, H2D_jetPtResponseMatrix_detectorResponse, H2D_jetPtResponseMatrix_fluctuations, iDataset, iRadius, options);
  ReweightResponseMatrixWithPrior(H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined, iDataset, iRadius, options);

  TString* pdfName = new TString("responseMatrix_combined"+(TString)"_R="+Form("%.1f",arrayRadius[iRadius])+"_"+Datasets[iDataset]+DatasetsNames[iDataset]);
  TString* pdfName_logz = new TString("responseMatrix_combined"+(TString)"_R="+Form("%.1f",arrayRadius[iRadius])+"_"+Datasets[iDataset]+DatasetsNames[iDataset]+"_logz");

  TString textContext(contextCustomOneField(*texDatasetsComparisonCommonDenominator, ""));

  Draw_TH2_Histogram(H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined, textContext, pdfName, texPtJetRecX, texPtJetGenX, texCollisionDataInfo, drawnWindowAuto, th2ContoursNone, contourNumberNone, "");
  Draw_TH2_Histogram(H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined, textContext, pdfName_logz, texPtJetRecX, texPtJetGenX, texCollisionDataInfo, drawnWindowAuto, th2ContoursNone, contourNumberNone, "logz");
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
  
  cout << "Draw_Pt_spectrum_unfolded test - -2" << endl;

  unfoldParameter = Get_Pt_spectrum_unfolded(H1D_jetPt_unfolded, iDataset, iRadius, unfoldParameterInput, options).first;
  cout << "Draw_Pt_spectrum_unfolded test - -1" << endl;

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
  cout << "Draw_Pt_spectrum_unfolded test - 0" << endl;

  // comparison with refolded
  if (!useFineBinningTest) {
    Get_Pt_spectrum_bkgCorrected_recBinning(H1D_jetPt_measured, iDataset, iRadius, options);
  } else {
    Get_Pt_spectrum_bkgCorrected_fineBinning(H1D_jetPt_measured, iDataset, iRadius, options);
  }
  Get_Pt_spectrum_unfoldedThenRefolded(H1D_jetPt_unfoldedThenRefolded, iDataset, iRadius, unfoldParameterInput, options);
  Get_Pt_spectrum_unfoldedThenRefolded_RooUnfoldMethod(H1D_jetPt_unfoldedThenRefolded_RooUnfoldMethod, iDataset, iRadius, unfoldParameterInput, options);
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
  cout << "Draw_Pt_spectrum_unfolded test - 1" << endl;

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
  cout << "Draw_Pt_spectrum_unfolded test - 2" << endl;

  TString pdfTitleBase = (TString)"IterationsDump/jet_";//+Datasets[iDataset]+DatasetsNames[iDataset]+"_R="+Form("%.1f", arrayRadius[iRadius])+"_Pt_unfolded_";
  // std::array<std::array<float, 2>, 2> drawnWindow = {{{ptWindowDisplay[0], ptWindowDisplay[1]}, {-999, -999}}}; // {{xmin, xmax}, {ymin, ymax}}

    // comparison with raw measured
  TString unfoldedMeasuredCompLegend[2] = {"unfolded data", "measured raw (gen binning)"};
  TString* pdfName_measuredComp = new TString(pdfTitleBase+unfoldingInfo+"_measuredComp");
  Draw_TH1_Histograms(H1D_jetPt_unfolded_measuredComp, unfoldedMeasuredCompLegend, 2, textContext, pdfName_measuredComp, texPtX, yAxisLabel, texCollisionDataInfo, drawnWindow, legendPlacementAuto, contextPlacementAuto, "logy");
  if (divideSuccessMeasured){
    TString* pdfName_ratio_measured = new TString(pdfTitleBase+unfoldingInfo+"_ratioMeasured");
    Draw_TH1_Histogram(H1D_jetPt_ratio_measured, textContext, pdfName_ratio_measured, texPtX, texRatioMeasuredUnfolded, texCollisionDataInfo, drawnWindow, legendPlacementAuto, contextPlacementAuto, "zoomToOneLarge,ratioLine");
  }
  cout << "Draw_Pt_spectrum_unfolded test - 3" << endl;

    // comparison with mcp truth
  TString unfoldedTruthCompLegend[2] = {"unfolded data", "mcp truth"};
  TString* pdfName_mcpComp = new TString(pdfTitleBase+unfoldingInfo+"_mcpComp");
  Draw_TH1_Histograms(H1D_jetPt_unfolded_mcpComp, unfoldedTruthCompLegend, 2, textContext, pdfName_mcpComp, texPtX, yAxisLabel, texCollisionDataInfo, drawnWindow, legendPlacementAuto, contextPlacementAuto, "logy");
  if (divideSuccessMcp){
    TString* pdfName_ratio_mcp = new TString(pdfTitleBase+unfoldingInfo+"_ratioMcp");
    Draw_TH1_Histogram(H1D_jetPt_ratio_mcp, textContext, pdfName_ratio_mcp, texPtX, texRatioMcpUnfolded, texCollisionDataInfo, drawnWindow, legendPlacementAuto, contextPlacementAuto, "zoomToOneLarge,ratioLine");
  }

  cout << "Draw_Pt_spectrum_unfolded test - 4" << endl;

  // comparison with refolded
  TString unfoldedRefoldedCompLegend[3] = {"refolded manually", "refolded roounfold", "measured"};
  TString* pdfName_refoldedComp = new TString(pdfTitleBase+unfoldingInfo+"_RefoldedComp");
  Draw_TH1_Histograms(H1D_jetPt_unfolded_refoldedComp, unfoldedRefoldedCompLegend, 3, textContext, pdfName_refoldedComp, texPtX, yAxisLabel, texCollisionDataInfo, drawnWindow, legendPlacementAuto, contextPlacementAuto, "logy");
  if (divideSuccessMeasuredRefolded[0] && divideSuccessMeasuredRefolded[1]) {
    TString* pdfName_ratio_refoldedComp = new TString(pdfTitleBase+unfoldingInfo+"_ratioRefoldedUnfolded");
    Draw_TH1_Histograms(H1D_jetPt_ratio_measuredRefolded, unfoldedRefoldedCompLegend, 2, textContext, pdfName_ratio_refoldedComp, texPtX, texRatioRefoldedMeasured, texCollisionDataInfo, drawnWindow, legendPlacementAuto, contextPlacementAuto, "zoomToOneLarge,ratioLine");
  }

  cout << "Draw_Pt_spectrum_unfolded test - 5" << endl;

  // comparison with Run 2
  if (isDataPbPb && comparePbPbWithRun2) {
  TString unfoldedRun2CompLegend[2] = {"unfolded Run3", "unfolded Run2 0-10%"};
    TString* pdfName_run2Comp = new TString(pdfTitleBase+unfoldingInfo+"_run2Comp");
    Draw_TH1_Histograms(H1D_jetPt_unfolded_run2Comp, unfoldedRun2CompLegend, 2, textContext, pdfName_run2Comp, texPtX, yAxisLabel, texCollisionDataInfo, drawnWindow, legendPlacementAuto, contextPlacementAuto, "logy");
    if (divideSuccessRun2){
      TString* pdfName_ratio_run2 = new TString(pdfTitleBase+unfoldingInfo+"_ratioRun2");
      Draw_TH1_Histogram(H1D_jetPt_ratio_run2, textContext, pdfName_ratio_run2, texPtX, texRatioRun2Unfolded, texCollisionDataInfo, drawnWindow, legendPlacementAuto, contextPlacementAuto, "zoomToOneLarge,ratioLine");
    }
  }
  cout << "Draw_Pt_spectrum_unfolded test - 6" << endl;

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
    Get_Pt_spectrum_unfoldedThenRefolded(H1D_jetPt_unfoldedThenRefolded[iUnfoldIteration], iDataset, iRadius, unfoldIterations, options);
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

  TString textContext(contextCustomOneField(*texDatasetsComparisonCommonDenominator, ""));

  TString* pdfName = new TString("jet_"+jetType[iJetType]+"_"+jetLevel[iJetLevel]+"_"+partialUniqueSpecifier+"_Pt_unfolded_"+unfoldingInfo);

  TString* yAxisLabel = texCount;
  if (normaliseDistribsBeforeUnfolding || normaliseDistribsAfterUnfolding) { //should probably check if having both on doesn't lead to double normalisation
    yAxisLabel = texJetPtYield_EventNorm;
  }
  Draw_TH1_Histograms(H1D_jetPt_unfolded, unfoldingIterationLegend, nUnfoldIteration, textContext, pdfName, texPtX, yAxisLabel, texCollisionDataInfo, drawnWindow, legendPlacementAuto, contextPlacementAuto, "logy");


    // comparison with raw measured
  // TString unfoldedMeasuredCompLegend[2] = {"unfolded data", "measured raw (gen binning)"};
  unfoldingIterationLegend[nUnfoldIteration] = (TString)"raw measured";
  TString* pdfName_measuredComp = new TString("jet_"+jetType[iJetType]+"_"+jetLevel[iJetLevel]+"_"+partialUniqueSpecifier+"_Pt_unfolded_"+unfoldingInfo+"_measuredComp");
  Draw_TH1_Histograms(H1D_jetPt_unfolded_measuredComp, unfoldingIterationLegend, nUnfoldIteration+1, textContext, pdfName_measuredComp, texPtX, yAxisLabel, texCollisionDataInfo, drawnWindow, legendPlacementAuto, contextPlacementAuto, "logy");

  bool divideSuccessMeasured_boolsum = true;
  for(int iUnfoldIteration = 0; iUnfoldIteration < nUnfoldIteration; iUnfoldIteration++){
    if (!(divideSuccessMeasured[iUnfoldIteration])) {
      divideSuccessMeasured_boolsum = false;
    }
  }
  if (divideSuccessMeasured_boolsum){
    TString* pdfName_ratio_measured = new TString("jet_"+jetType[iJetType]+"_"+jetLevel[iJetLevel]+"_"+partialUniqueSpecifier+"_Pt_unfolded_"+unfoldingInfo+"_ratioMeasured");
    Draw_TH1_Histograms(H1D_jetPt_ratio_measured, unfoldingIterationLegend, nUnfoldIteration, textContext, pdfName_ratio_measured, texPtX, texRatioMeasuredUnfolded, texCollisionDataInfo, drawnWindow, legendPlacementAuto, contextPlacementAuto, "zoomToOneLarge,ratioLine");
  }


    // comparison with mcp truth
  // TString unfoldedTruthCompLegend[2] = {"unfolded data", "mcp truth"};
  unfoldingIterationLegend[nUnfoldIteration] = (TString)"mcp";
  TString* pdfName_mcpComp = new TString("jet_"+jetType[iJetType]+"_"+jetLevel[iJetLevel]+"_"+partialUniqueSpecifier+"_Pt_unfolded_"+unfoldingInfo+"_mcpComp");
  Draw_TH1_Histograms(H1D_jetPt_unfolded_mcpComp, unfoldingIterationLegend, nUnfoldIteration+1, textContext, pdfName_mcpComp, texPtX, yAxisLabel, texCollisionDataInfo, drawnWindow, legendPlacementAuto, contextPlacementAuto, "logy");

  bool divideSuccessMcp_boolsum = true;
  for(int iUnfoldIteration = 0; iUnfoldIteration < nUnfoldIteration; iUnfoldIteration++){
    if (!(divideSuccessMcp[iUnfoldIteration])) {
      divideSuccessMcp_boolsum = false;
    }
  }
  if (divideSuccessMcp_boolsum){
    TString* pdfName_ratio_mcp = new TString("jet_"+jetType[iJetType]+"_"+jetLevel[iJetLevel]+"_"+partialUniqueSpecifier+"_Pt_unfolded_"+unfoldingInfo+"_ratioMcp");
    Draw_TH1_Histograms(H1D_jetPt_ratio_mcp, unfoldingIterationLegend, nUnfoldIteration, textContext, pdfName_ratio_mcp, texPtX, texRatioMcpUnfolded, texCollisionDataInfo, drawnWindow, legendPlacementAuto, contextPlacementAuto, "zoomToOneLarge,ratioLine");
  }

  // comparison with refolded
  // TString unfoldedRefoldedCompLegend[2] = {"refolded", "measured"};
  unfoldingIterationLegend[nUnfoldIteration] = (TString)"measured";
  for(int iIteration = 0; iIteration < nUnfoldIteration; iIteration++){
    unfoldingIterationLegend[iIteration] += (TString)" refolded";
  }
  TString* pdfName_refoldedComp = new TString("jet_"+jetType[iJetType]+"_"+jetLevel[iJetLevel]+"_"+partialUniqueSpecifier+"_Pt_unfolded_"+unfoldingInfo+"_RefoldedComp");
  Draw_TH1_Histograms(H1D_jetPt_unfolded_refoldedComp, unfoldingIterationLegend, nUnfoldIteration+1, textContext, pdfName_refoldedComp, texPtX, yAxisLabel, texCollisionDataInfo, drawnWindow, legendPlacementAuto, contextPlacementAuto, "logy");

  bool divideSuccessRefoldedComp_boolsum = true;
  for(int iUnfoldIteration = 0; iUnfoldIteration < nUnfoldIteration; iUnfoldIteration++){
    if (!(divideSuccessMeasuredRefolded[iUnfoldIteration])) {
      divideSuccessRefoldedComp_boolsum = false;
    }
  }
  if (divideSuccessRefoldedComp_boolsum){
    TString* pdfName_ratio_refoldedComp = new TString("jet_"+jetType[iJetType]+"_"+jetLevel[iJetLevel]+"_"+partialUniqueSpecifier+"_Pt_unfolded_"+unfoldingInfo+"_ratioRefoldedUnfolded");
    Draw_TH1_Histograms(H1D_jetPt_ratio_measuredRefolded, unfoldingIterationLegend, nUnfoldIteration, textContext, pdfName_ratio_refoldedComp, texPtX, texRatioRefoldedMeasured, texCollisionDataInfo, drawnWindow, legendPlacementAuto, contextPlacementAuto, "zoomToOneLarge,ratioLine");
  }
}




// rename refoldedUnfolded as closure test?
// and try and spend 15 min to clean hist names for the spectrum analysis



// WARNING FOR EFFICIENCIES I SHOULD REREAD THIS BELOW!!
// hMcEfficiency_vsPt->Divide(hMcSignalCount_vsPt,TrueV0PtSpectrum_AnalysisBins, 1., 1., "b"); // option b for binomial because efficiency: https://twiki.cern.ch/twiki/bin/view/ALICE/PWGLFPAGSTRANGENESSEfficiency
