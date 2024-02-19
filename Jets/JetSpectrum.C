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
void Draw_Pt_spectrum_raw(int iDataset, int iRadius, const char options[]);
void Get_PtResponseMatrix_Fluctuations(TH2D* &H2D_jetPtResponseMatrix_fluctuations, int iDataset, int iRadius, float* centRange);
void Get_PtResponseMatrix_Fluctuations_FullResForCheck(TH2D* &H2D_jetPtResponseMatrix_fluctuations, int iDataset, int iRadius, float* centRange);
void Draw_ResponseMatrices_Fluctuations(int iDataset, int iRadius);
void Get_Pt_spectrum_truth(TH1D* &H1D_jetPt_truth, int iDataset, int iRadius, float* centRange, const char options[]);
void Draw_Pt_spectrum_truth(int iDataset, int iRadius, const char options[]);
void Get_Pt_spectrum_unfolded(TH1D* &H1D_jetPt_unfolded, int iDataset, int iRadius, float* centRange, const char options[]);
void Draw_Pt_spectrum_unfolded_FluctResponseOnly(int iDataset, int iRadius, const char options[]);

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

  Draw_Pt_spectrum_raw(iDataset, iRadius, "entriesNorm");
  Draw_ResponseMatrices_Fluctuations(iDataset, iRadius);
  Draw_Pt_spectrum_unfolded_FluctResponseOnly(iDataset, iRadius, "entriesNorm");
  Draw_Pt_spectrum_truth(iDataset, iRadius, "entriesNorm");

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

TString contextDatasetRadiusCompAndVarRange(int iDataset, float* variableRange, const char options[]){
  TString texcontextDatasetRadiusCompAndVarRange;
  if (strstr(options, "pt") != NULL) { //  || strstr(options, "ratio") != NULL not sure why I had this here
    texcontextDatasetRadiusCompAndVarRange = "#splitline{"+*texDatasetsComparisonCommonDenominator+" "+DatasetsNames[iDataset]+"}{#splitline{2023 QC}{"+contextPtRange(variableRange)+"}}";
  }
  if (strstr(options, "eta") != NULL) { //  || strstr(options, "ratio") != NULL not sure why I had this here
    texcontextDatasetRadiusCompAndVarRange = "#splitline{"+*texDatasetsComparisonCommonDenominator+" "+DatasetsNames[iDataset]+"}{#splitline{2023 QC}{"+contextEtaRange(variableRange)+"}}";
  }
  if (strstr(options, "cent") != NULL) { //  || strstr(options, "ratio") != NULL not sure why I had this here
    texcontextDatasetRadiusCompAndVarRange = "#splitline{"+*texDatasetsComparisonCommonDenominator+" "+DatasetsNames[iDataset]+"}{#splitline{2023 QC}{"+contextCentRange(variableRange)+"}}";
  }

  return texcontextDatasetRadiusCompAndVarRange;
}

TString contextDatasetCompAndRadiusAndVarRange(float jetRadius, float* variableRange, const char options[]){
  TString texcontextDatasetCompAndRadiusAndVarRange;
  if (strstr(options, "pt") != NULL) { //  || strstr(options, "ratio") != NULL not sure why I had this here
    texcontextDatasetCompAndRadiusAndVarRange = "#splitline{"+*texDatasetsComparisonCommonDenominator+"}{#splitline{"+contextJetRadius(jetRadius)+"}{"+contextPtRange(variableRange)+"}}";
  }
  if (strstr(options, "eta") != NULL) { //  || strstr(options, "ratio") != NULL not sure why I had this here
    texcontextDatasetCompAndRadiusAndVarRange = "#splitline{"+*texDatasetsComparisonCommonDenominator+"}{#splitline{"+contextJetRadius(jetRadius)+"}{"+contextEtaRange(variableRange)+"}}";
  }

  return texcontextDatasetCompAndRadiusAndVarRange;
}

TString contextDatasetCompAndRadius(float jetRadius, const char options[]){
  TString texcontextDatasetCompAndRadius;
  texcontextDatasetCompAndRadius = "#splitline{"+*texDatasetsComparisonCommonDenominator+"}{"+contextJetRadius(jetRadius)+"}";

  return texcontextDatasetCompAndRadius;
}

TString contextDatasetComp(const char options[]){
  TString texcontextDatasetComp;
  texcontextDatasetComp = *texDatasetsComparisonCommonDenominator;

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

TString contextCentRange(float* centRange){
  std::stringstream ss;
  ss << centRange[0] << " < #cent < " << centRange[1];
  TString textContext((TString)ss.str());
  return textContext;
}

void CentralityLegend(TString* centralityLegend, const float* arrayCentralityBinning, int nCentralityBins){
  std::stringstream ss;
  ss.precision(2);
  for(int iCentralityBin = 0; iCentralityBin < nCentralityBins; iCentralityBin++){
    ss << "" << arrayCentralityBinning[iCentralityBin] << "-" << arrayCentralityBinning[iCentralityBin+1] << "%";
    centralityLegend[iCentralityBin] = (TString)ss.str();
    ss.str("");
    ss.clear();
  }
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////// Spectrum analysis functions ///////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


void Get_Pt_spectrum_bkgCorrected(TH1D* &H1D_jetPt, int iDataset, int iRadius, float* centRange, const char options[]) {
  TH3D* H3D_jetRjetPtcolCent;
  TH1D* H1D_jetPt_defaultBin;
  // TH1D* H1D_jetPt_raw[nRadius];

  H3D_jetRjetPtcolCent = (TH3D*)((TH3D*)file_O2Analysis_list[iDataset]->Get(analysisWorkflow+"/h3_jet_r_jet_pt_centrality_rhoareasubtracted"))->Clone("Get_Pt_spectrum_bkgCorrected"+Datasets[iDataset]+Form("%.1f", centRange[0])+"<cent<"+Form("%.1f", centRange[1]));
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
    NormaliseYieldToNEvents(H1D_jetPt, GetNEventsSel8(file_O2Analysis_list[iDataset]));
  }
  if (strstr(options, "entriesNorm") != NULL) {
    NormaliseYieldToNEntries(H1D_jetPt);
  }

}

void Get_Pt_spectrum_truth(TH1D* &H1D_jetPt_truth, int iDataset, int iRadius, float* centRange, const char options[]) {
  TH3D* H3D_jetRjetPtcolCent;
  TH1D* H1D_jetPt_truth_defaultBin;
  // TH1D* H1D_jetPt_raw[nRadius];

  H3D_jetRjetPtcolCent = (TH3D*)((TH3D*)file_O2Analysis_list[iDataset]->Get(analysisWorkflow+"/h3_jet_r_jet_pt_centrality_rhoareasubtracted"))->Clone("Get_Pt_spectrum_truth"+Datasets[iDataset]+Form("%.1f", centRange[0])+"<cent<"+Form("%.1f", centRange[1]));
  H3D_jetRjetPtcolCent->Sumw2();
  int binCentEdges[2] = {H3D_jetRjetPtcolCent->GetZaxis()->FindBin(centRange[0]), H3D_jetRjetPtcolCent->GetZaxis()->FindBin(centRange[1])};
  if (binCentEdges[0] == 0) 
    cout << "WARNING: Get_Pt_spectrum_bkgCorrected is counting the underflow with the chosen centRange" << endl;
  if (binCentEdges[1] == H3D_jetRjetPtcolCent->GetZaxis()->GetNbins()+1) 
    cout << "WARNING: Get_Pt_spectrum_bkgCorrected is counting the overflow with the chosen centRange" << endl;

  int ibinJetRadius = H3D_jetRjetPtcolCent->GetXaxis()->FindBin(arrayRadius[iRadius-1]+GLOBAL_epsilon);

  H1D_jetPt_truth_defaultBin = (TH1D*)H3D_jetRjetPtcolCent->ProjectionY("jetPt_truth_"+RadiusLegend[iRadius]+Form("%.1f", centRange[0])+"<cent<"+Form("%.1f", centRange[1]), ibinJetRadius, ibinJetRadius, binCentEdges[0], binCentEdges[1]);
  H1D_jetPt_truth = (TH1D*)H1D_jetPt_truth_defaultBin->Rebin(nBinPtJetsRec[iRadius],"jetPt_truth_rebinned_"+RadiusLegend[iRadius]+Form("%.1f", centRange[0])+"<cent<"+Form("%.1f", centRange[1]), ptBinsJetsRec[iRadius]);

  if (strstr(options, "evtNorm") != NULL) {
    // NormaliseYieldToNEvents(H1D_jetPt_truth, GetNEventsGen(file_O2Analysis_list[iDataset]));
    cout << "Normalisation to number of events is not implemented for MC truth data yet. Switching to normalisation to the number of entries (#jets)." << endl;
    NormaliseYieldToNEntries(H1D_jetPt_truth);
  }
  if (strstr(options, "entriesNorm") != NULL) {
    NormaliseYieldToNEntries(H1D_jetPt_truth);
  }

}

// void Get_PtResponseMatrix_DetectorEfficiency(TH2D* H2D_jetPtResponseMatrix_detectorEfficiency, int iDataset, float* centRange, const char options[]) { // to be created once the matching is ready
// }

// void Get_PtResponseMatrix_FluctuationsLegacyWrong(TH2D* &H2D_jetPtResponseMatrix_fluctuations, int iDataset, int iRadius, float* centRange) { 
//   // see Hiroki Yokoyama thesis
//   // iRadius is for chosing the pT binning

//   TString partialUniqueSpecifier = Datasets[iDataset]+"_R="+Form("%.1f",arrayRadius[iRadius])+"_@cent["+Form("%.1f", centRange[0])+","+Form("%.1f", centRange[1])+"]";

//   TH2D* H2D_fluctuations_centrality;
//   TH1D* H1D_fluctuations_highRes;
//   TH1D* H1D_pt_forBinningCopying;

//   H2D_fluctuations_centrality = (TH2D*)((TH2D*)file_O2Analysis_list[iDataset]->Get(analysisWorkflow+"/h2_centrality_rhoRandomCone"))->Clone("Get_PtResponseMatrix_Fluctuations"+Datasets[iDataset]);
//   H2D_fluctuations_centrality->Sumw2();

//   H1D_pt_forBinningCopying = (TH1D*)((TH3D*)file_O2Analysis_list[iDataset]->Get(analysisWorkflow+"/h3_jet_r_jet_pt_centrality"))->ProjectionY("jetPt_forBinningCopying_"+partialUniqueSpecifier, 0, -1, 0, -1);

//   int ibinCent_low, ibinCent_high;

//   ibinCent_low = H2D_fluctuations_centrality->GetXaxis()->FindBin(centRange[0]);
//   ibinCent_high = H2D_fluctuations_centrality->GetXaxis()->FindBin(centRange[1]);
//   H1D_fluctuations_highRes = (TH1D*)H2D_fluctuations_centrality->ProjectionY("bkgFluctuationCentrality_highRes_"+partialUniqueSpecifier, ibinCent_low, ibinCent_high);

//   NormaliseYieldToNEntries(H1D_fluctuations_highRes); // normalising fluctuations to 1
    
//   TH2D H2D_response_highRes = TH2D("H2D_response_highRes_"+partialUniqueSpecifier, "H2D_response_highRes_"+partialUniqueSpecifier, H1D_pt_forBinningCopying->GetNbinsX(), H1D_pt_forBinningCopying->GetXaxis()->GetXmin(), H1D_pt_forBinningCopying->GetXaxis()->GetXmax(), H1D_pt_forBinningCopying->GetNbinsX(), H1D_pt_forBinningCopying->GetXaxis()->GetXmin(), H1D_pt_forBinningCopying->GetXaxis()->GetXmax()); // actually doesn't work if original histogram has fixed bin size
//   TH2D H2D_response_highRes = TH2D("H2D_response_highRes_"+partialUniqueSpecifier, "H2D_response_highRes_"+partialUniqueSpecifier, H1D_pt_forBinningCopying->GetNbinsX(), H1D_pt_forBinningCopying->GetXaxis()->GetXmin(), H1D_pt_forBinningCopying->GetXaxis()->GetXmax(), H1D_pt_forBinningCopying->GetNbinsX(), H1D_pt_forBinningCopying->GetXaxis()->GetXmin(), H1D_pt_forBinningCopying->GetXaxis()->GetXmax()); // actually doesn't work if original histogram has fixed bin size
//   // TH2D* H2D_response_rebinned; 

//   //==================== rebin fluctuations to have same binning as AnalysisResults rec pt ====================//
//   std::vector<double> ptBinsVector = GetTH1Bins(H1D_pt_forBinningCopying);
//   double* ptBins = &ptBinsVector[0];
//   H1D_fluctuations_highRes->Rebin(H1D_pt_forBinningCopying->GetNbinsX(), "bkgFluctuationCentrality_highRes_AnalysisResutltsPtBinning"+partialUniqueSpecifier, ptBins);
//   //====================================== rebin end ==========================================//

  
//   // //==================== Build response matrix: shift deltaPt by pT gen along the pT rec axis ====================//
//   // int ibinZero= H1D_fluctuations_highRes->FindBin(0+GLOBAL_epsilon);
//   // for(int iBinGen = 1; iBinGen <= H1D_fluctuations_highRes->GetNbinsX(); iBinGen++){
//   //   for(int iBinRec = 1; iBinRec <= H1D_fluctuations_highRes->GetNbinsX(); iBinRec++){
//   //     if (ibinZero + (iBinRec - iBinGen) <= H1D_fluctuations_highRes->GetNbinsX()) { // make sure it's within bin range
//   //       H2D_response_highRes.SetBinContent(iBinRec, iBinGen, H1D_fluctuations_highRes->GetBinContent(ibinZero + (iBinRec - iBinGen)));
//   //       // H2D_response_highRes->SetBinError(iBinRec, iBinGen, something); 
//   //     }
//   //   }
//   // }
//   // //========================================= Build response matrix end =========================================//

//   //==================== Build response matrix: shift deltaPt by pT gen along the pT rec axis ====================//
//   int ibinZero= H1D_fluctuations_highRes->FindBin(0+GLOBAL_epsilon);
//   for(int iBinGen = 1; iBinGen <= H1D_fluctuations_highRes->GetNbinsX(); iBinGen++){
//     for(int iBinRec = 1; iBinRec <= H1D_fluctuations_highRes->GetNbinsX(); iBinRec++){
//       if (ibinZero + (iBinRec - iBinGen) <= H1D_fluctuations_highRes->GetNbinsX()) { // make sure it's within bin range
//         H2D_response_highRes.SetBinContent(iBinRec, iBinGen, H1D_fluctuations_highRes->GetBinContent(ibinZero + (iBinRec - iBinGen)));
//         // H2D_response_highRes->SetBinError(iBinRec, iBinGen, something); 
//       }
//     }
//   }
//   //========================================= Build response matrix end =========================================//

//   H2D_jetPtResponseMatrix_fluctuations = (TH2D*)H2D_response_highRes.Clone("H2D_jetPtResponseMatrix_fluctuations"+partialUniqueSpecifier); // to be replaced with pT binning from analysis; rebin 2D with variable bins doesn't work, gonna have to make a custom function that does it
//   // TH2D H2D_jetPtResponseMatrix_fluctuations_rebinned = RebinVariableBins2D(H2D_jetPtResponseMatrix_fluctuations, nBinPtJetsRec[iRadius], nBinPtJetsGen[iRadius], ptBinsJetsRec[iRadius], ptBinsJetsGen[iRadius]);

//   // cout << "H2D_jetPtResponseMatrix_fluctuations nbins X =" << H2D_jetPtResponseMatrix_fluctuations->GetNbinsX() << endl;
//   // cout << "H2D_jetPtResponseMatrix_fluctuations_rebinned nbins X =" << H2D_jetPtResponseMatrix_fluctuations_rebinned.GetNbinsX() << endl;
//   H2D_jetPtResponseMatrix_fluctuations = (TH2D*)RebinVariableBins2D(H2D_jetPtResponseMatrix_fluctuations, nBinPtJetsRec[iRadius], nBinPtJetsGen[iRadius], ptBinsJetsRec[iRadius], ptBinsJetsGen[iRadius]).Clone("H2D_jetPtResponseMatrix_fluctuations_rebinned_"+partialUniqueSpecifier);

//   cout << "normalisation looks wrong; maybe need to normalise each ptRec line to unity? else I have bins with values 10" << endl;
// }


void Get_PtResponseMatrix_Fluctuations(TH2D* &H2D_jetPtResponseMatrix_fluctuations, int iDataset, int iRadius, float* centRange) { 
  // see Hiroki Yokoyama thesis
  // iRadius is for chosing the pT binning

  TString partialUniqueSpecifier = Datasets[iDataset]+"_R="+Form("%.1f",arrayRadius[iRadius])+"_@cent["+Form("%.1f", centRange[0])+","+Form("%.1f", centRange[1])+"]";

  TH2D* H2D_fluctuations_centrality;
  TH1D* H1D_fluctuations;

  H2D_fluctuations_centrality = (TH2D*)((TH2D*)file_O2Analysis_list[iDataset]->Get(analysisWorkflow+"/h2_centrality_rhoRandomCone"))->Clone("Get_PtResponseMatrix_Fluctuations"+Datasets[iDataset]);
  H2D_fluctuations_centrality->Sumw2();

  int ibinCent_low, ibinCent_high;

  ibinCent_low = H2D_fluctuations_centrality->GetXaxis()->FindBin(centRange[0]);
  ibinCent_high = H2D_fluctuations_centrality->GetXaxis()->FindBin(centRange[1]);
  H1D_fluctuations = (TH1D*)H2D_fluctuations_centrality->ProjectionY("bkgFluctuationCentrality_highRes_"+partialUniqueSpecifier, ibinCent_low, ibinCent_high);

  NormaliseYieldToNEntries(H1D_fluctuations); // normalising fluctuations to 1
    
  TH2D H2D_response = TH2D("H2D_response_"+partialUniqueSpecifier, "H2D_response_"+partialUniqueSpecifier, nBinPtJetsRec[iRadius], ptBinsJetsRec[iRadius], nBinPtJetsGen[iRadius], ptBinsJetsGen[iRadius]); // actually doesn't work if original histogram has fixed bin size

  //==================== Build response matrix: shift deltaPt by pT gen along the pT rec axis ====================//
  int ibinZeroFluct= H1D_fluctuations->FindBin(0+GLOBAL_epsilon);
  
  for(int iBinGen = 1; iBinGen <= H2D_response.GetNbinsX(); iBinGen++){
    for(int iBinRec = 1; iBinRec <= H2D_response.GetNbinsX(); iBinRec++){
      double ptGen = H2D_response.GetYaxis()->GetBinCenter(iBinGen);
      double ptRec_low = H2D_response.GetXaxis()->GetBinLowEdge(iBinRec);
      double ptRec_up = H2D_response.GetXaxis()->GetBinLowEdge(iBinRec+1);
      // double xPtRecWidth = H2D_response.GetXaxis()->GetBinWidth(iBinRec);
      // if (ibinZero + (iBinRec - iBinGen) <= H1D_fluctuations_highRes->GetNbinsX()) { // make sure it's within bin range
      H2D_response.SetBinContent(iBinRec, iBinGen, H1D_fluctuations->Integral(H1D_fluctuations->GetXaxis()->FindBin(ptRec_low - ptGen + GLOBAL_epsilon), H1D_fluctuations->GetXaxis()->FindBin(ptRec_up - ptGen + GLOBAL_epsilon))); 
        // H2D_response_highRes->SetBinError(iBinRec, iBinGen, something); 
    }
  }
  // }
  //========================================= Build response matrix end =========================================//

  H2D_jetPtResponseMatrix_fluctuations = (TH2D*)H2D_response.Clone("H2D_jetPtResponseMatrix_fluctuations"+partialUniqueSpecifier); // to be replaced with pT binning from analysis; rebin 2D with variable bins doesn't work, gonna have to make a custom function that does it
  // TH2D H2D_jetPtResponseMatrix_fluctuations_rebinned = RebinVariableBins2D(H2D_jetPtResponseMatrix_fluctuations, nBinPtJetsRec[iRadius], nBinPtJetsGen[iRadius], ptBinsJetsRec[iRadius], ptBinsJetsGen[iRadius]);

  // cout << "H2D_jetPtResponseMatrix_fluctuations nbins X =" << H2D_jetPtResponseMatrix_fluctuations->GetNbinsX() << endl;
  // cout << "H2D_jetPtResponseMatrix_fluctuations_rebinned nbins X =" << H2D_jetPtResponseMatrix_fluctuations_rebinned.GetNbinsX() << endl;
  H2D_jetPtResponseMatrix_fluctuations = (TH2D*)RebinVariableBins2D(H2D_jetPtResponseMatrix_fluctuations, nBinPtJetsRec[iRadius], nBinPtJetsGen[iRadius], ptBinsJetsRec[iRadius], ptBinsJetsGen[iRadius]).Clone("H2D_jetPtResponseMatrix_fluctuations_rebinned_"+partialUniqueSpecifier);

  cout << "normalisation looks wrong; maybe need to normalise each ptRec line to unity? else I have bins with values 10" << endl;
}

void Get_PtResponseMatrix_Fluctuations_FullResForCheck(TH2D* &H2D_jetPtResponseMatrix_fluctuations, int iDataset, int iRadius, float* centRange) { 
  // see Hiroki Yokoyama thesis
  // iRadius is for chosing the pT binning

  TString partialUniqueSpecifier = Datasets[iDataset]+"_R="+Form("%.1f",arrayRadius[iRadius])+"_@cent["+Form("%.1f", centRange[0])+","+Form("%.1f", centRange[1])+"]";

  TH2D* H2D_fluctuations_centrality;
  TH1D* H1D_fluctuations_highRes;
  TH1D* H1D_pt_forBinningCopying;

  H2D_fluctuations_centrality = (TH2D*)((TH2D*)file_O2Analysis_list[iDataset]->Get(analysisWorkflow+"/h2_centrality_rhoRandomCone"))->Clone("Get_PtResponseMatrix_Fluctuations"+Datasets[iDataset]);
  H2D_fluctuations_centrality->Sumw2();

  H1D_pt_forBinningCopying = (TH1D*)((TH3D*)file_O2Analysis_list[iDataset]->Get(analysisWorkflow+"/h3_jet_r_jet_pt_centrality"))->ProjectionY("jetPt_forBinningCopying_"+partialUniqueSpecifier, 0, -1, 0, -1);

  int ibinCent_low, ibinCent_high;

  ibinCent_low = H2D_fluctuations_centrality->GetXaxis()->FindBin(centRange[0]);
  ibinCent_high = H2D_fluctuations_centrality->GetXaxis()->FindBin(centRange[1]);
  H1D_fluctuations_highRes = (TH1D*)H2D_fluctuations_centrality->ProjectionY("bkgFluctuationCentrality_highRes_"+partialUniqueSpecifier, ibinCent_low, ibinCent_high);

  NormaliseYieldToNEntries(H1D_fluctuations_highRes); // normalising fluctuations to 1
    
  TH2D H2D_response_highRes = TH2D("H2D_response_highRes_"+partialUniqueSpecifier, "H2D_response_highRes_"+partialUniqueSpecifier, H1D_pt_forBinningCopying->GetNbinsX(), H1D_pt_forBinningCopying->GetXaxis()->GetXmin(), H1D_pt_forBinningCopying->GetXaxis()->GetXmax(), H1D_pt_forBinningCopying->GetNbinsX(), H1D_pt_forBinningCopying->GetXaxis()->GetXmin(), H1D_pt_forBinningCopying->GetXaxis()->GetXmax()); // actually doesn't work if original histogram has fixed bin size
  // TH2D* H2D_response_rebinned; 

  //==================== rebin fluctuations to have same binning as AnalysisResults rec pt ====================//
  std::vector<double> ptBinsVector = GetTH1Bins(H1D_pt_forBinningCopying);
  double* ptBins = &ptBinsVector[0];
  H1D_fluctuations_highRes->Rebin(H1D_pt_forBinningCopying->GetNbinsX(), "bkgFluctuationCentrality_highRes_AnalysisResutltsPtBinning"+partialUniqueSpecifier, ptBins);
  //====================================== rebin end ==========================================//

  
  //==================== Build response matrix: shift deltaPt by pT gen along the pT rec axis ====================//
  int ibinZero= H1D_fluctuations_highRes->FindBin(0+GLOBAL_epsilon);
  for(int iBinGen = 1; iBinGen <= H1D_fluctuations_highRes->GetNbinsX(); iBinGen++){
    for(int iBinRec = 1; iBinRec <= H1D_fluctuations_highRes->GetNbinsX(); iBinRec++){
      if (ibinZero + (iBinRec - iBinGen) <= H1D_fluctuations_highRes->GetNbinsX()) { // make sure it's within bin range
        H2D_response_highRes.SetBinContent(iBinRec, iBinGen, H1D_fluctuations_highRes->GetBinContent(ibinZero + (iBinRec - iBinGen)));
        // H2D_response_highRes->SetBinError(iBinRec, iBinGen, something); 
      }
    }
  }
  //========================================= Build response matrix end =========================================//

  H2D_jetPtResponseMatrix_fluctuations = (TH2D*)H2D_response_highRes.Clone("H2D_jetPtResponseMatrix_fluctuations"+partialUniqueSpecifier); // to be replaced with pT binning from analysis; rebin 2D with variable bins doesn't work, gonna have to make a custom function that does it
  // TH2D H2D_jetPtResponseMatrix_fluctuations_rebinned = RebinVariableBins2D(H2D_jetPtResponseMatrix_fluctuations, nBinPtJetsRec[iRadius], nBinPtJetsGen[iRadius], ptBinsJetsRec[iRadius], ptBinsJetsGen[iRadius]);

  // cout << "H2D_jetPtResponseMatrix_fluctuations nbins X =" << H2D_jetPtResponseMatrix_fluctuations->GetNbinsX() << endl;
  // cout << "H2D_jetPtResponseMatrix_fluctuations_rebinned nbins X =" << H2D_jetPtResponseMatrix_fluctuations_rebinned.GetNbinsX() << endl;
  // H2D_jetPtResponseMatrix_fluctuations = (TH2D*)RebinVariableBins2D(H2D_jetPtResponseMatrix_fluctuations, nBinPtJetsRec[iRadius], nBinPtJetsGen[iRadius], ptBinsJetsRec[iRadius], ptBinsJetsGen[iRadius]).Clone("H2D_jetPtResponseMatrix_fluctuations_rebinned_"+partialUniqueSpecifier);
}

void Get_Pt_spectrum_unfolded(TH1D* &H1D_jetPt_unfolded, int iDataset, int iRadius, float* centRange, const char options[]) { // to be tested
  //for now makes the assumption gen and rec have the same pT binning

  TH1D* measured;
  Get_Pt_spectrum_bkgCorrected(measured, iDataset, iRadius, centRange, options);
  TH1D* truth;
  Get_Pt_spectrum_truth(truth, iDataset, iRadius, centRange, options);

  TH2D* H2D_jetPtResponseMatrix;
  Get_PtResponseMatrix_Fluctuations(H2D_jetPtResponseMatrix, iDataset, iRadius, centRange);
  //Get_PtResponseMatrix_Detector()
  // compute matrixFluctuations times matrixDetector

  //  = (TH1D*)H1D_jetPt->Clone("truth_jetPt_"+RadiusLegend[iRadius]+Form("%.1f", centRange[0])+"<cent<"+Form("%.1f", centRange[1]));
  TH1D* truth_rebinned = (TH1D*)truth->Rebin(1., "truth_rebinned_jetPt_"+RadiusLegend[iRadius]+Form("%.1f", centRange[0])+"<cent<"+Form("%.1f", centRange[1]));
  
  // // TH1D* projX = (TH1D*)H2D_jetPtResponseMatrix[iRadius]->ProjectionX("projX",0,-1);
  // // TH1D* projY = (TH1D*)H2D_jetPtResponseMatrix[iRadius]->ProjectionY("projY",0,-1);
  // cout << "Debug log for unfolding - start" << endl;
  // for(int ibinPt = 1; ibinPt <= nBinPtJetsRec[iRadius]; ibinPt++){
  //   cout << "measured(" << ibinPt << ") = " << measured->GetBinContent(ibinPt) << endl;
  //   cout << "truth_rebinned(" << ibinPt << ") = " << truth_rebinned->GetBinContent(ibinPt) << endl;
  //   for(int jbinPt = 1; jbinPt <= nBinPtJetsRec[iRadius]; jbinPt++){
  //     cout << "H2D_jetPtResponseMatrix(" << ibinPt << "," << jbinPt << ") = " << H2D_jetPtResponseMatrix->GetBinContent(ibinPt,jbinPt) << endl;
  //   }
  // }
  // // cout << "Debug log for unfolding - end" << endl;

  RooUnfoldResponse* response = new RooUnfoldResponse(measured, truth_rebinned, H2D_jetPtResponseMatrix);
  // RooUnfoldResponse* response = new RooUnfoldResponse(projX, projY, H2D_jetPtResponseMatrix, "", "");
  // RooUnfold* unfold = new RooUnfoldBayes(response, measured[iRadius], 4);
  RooUnfold* unfold = new RooUnfoldBinByBin(response, measured);
  // unfold->IncludeSystematics(RooUnfolding::kAll);
  // TH1D* hist_unfold = (TH1D*)unfold.Hunfold();
  // TH1D* hist_unfold = (TH1D*)unfold->Hreco()->Clone("hist_unfold");
  TH1D* hist_unfold = static_cast<TH1D*>(unfold->Hreco());


  H1D_jetPt_unfolded = (TH1D*)hist_unfold->Clone("hRawYield_vsPt");
}


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////// Spectrum plotting functions ///////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


void Draw_Pt_spectrum_raw(int iDataset, int iRadius, const char options[]) {

  TH1D* H1D_jetPt_raw[nCentralityBins];
  float centRange[2];
  for(int iCent = 0; iCent < nCentralityBins; iCent++){
    centRange[0] = arrayCentralityBinning[iCent];
    centRange[1] = arrayCentralityBinning[iCent+1];

    Get_Pt_spectrum_bkgCorrected(H1D_jetPt_raw[iCent], iDataset, iRadius, centRange, options);
  }


  TString* pdfName = new TString("jet_"+jetType[iJetType]+"_"+jetLevel[iJetLevel]+"_"+Datasets[iDataset]+"_R="+Form("%.1f", arrayRadius[iRadius])+"_Pt_raw");

  TString textContext(contextDatasetComp("")); // placeholder
  
  TString centralityLegend[nCentralityBins]; CentralityLegend(centralityLegend, arrayCentralityBinning, nCentralityBins);

  TString* yAxisLabel;
  if (strstr(options, "evtNorm") != NULL) {
    yAxisLabel = texJetPtYield_EventNorm;
  }
  if (strstr(options, "entriesNorm") != NULL) {
    yAxisLabel = texJetPtYield_EntriesNorm;
  }

  Draw_TH1_Histograms_in_one(H1D_jetPt_raw, centralityLegend, nCentralityBins, textContext, pdfName, texPtX, yAxisLabel, texCollisionDataInfo, "logy");
}

void Draw_Pt_spectrum_truth(int iDataset, int iRadius, const char options[]) {

  TH1D* H1D_jetPt_truth[nCentralityBins];
  float centRange[2];
  for(int iCent = 0; iCent < nCentralityBins; iCent++){
    centRange[0] = arrayCentralityBinning[iCent];
    centRange[1] = arrayCentralityBinning[iCent+1];

    Get_Pt_spectrum_truth(H1D_jetPt_truth[iCent], iDataset, iRadius, centRange, options);
  }


  TString* pdfName = new TString("jet_"+jetType[iJetType]+"_"+jetLevel[iJetLevel]+"_"+Datasets[iDataset]+"_R="+Form("%.1f", arrayRadius[iRadius])+"_Pt_truth");

  TString textContext(contextDatasetComp("")); // placeholder
  
  TString centralityLegend[nCentralityBins]; CentralityLegend(centralityLegend, arrayCentralityBinning, nCentralityBins);

  TString* yAxisLabel;
  if (strstr(options, "evtNorm") != NULL) {
    yAxisLabel = texJetPtYield_EventNorm;
  }
  if (strstr(options, "entriesNorm") != NULL) {
    yAxisLabel = texJetPtYield_EntriesNorm;
  }

  Draw_TH1_Histograms_in_one(H1D_jetPt_truth, centralityLegend, nCentralityBins, textContext, pdfName, texPtX, yAxisLabel, texCollisionDataInfo, "logy");
}

void Draw_ResponseMatrices_Fluctuations(int iDataset, int iRadius) {

  TH2D* H2D_jetPtResponseMatrix_fluctuations[nCentralityBins];
  TH2D* H2D_jetPtResponseMatrix_fluctuations_fullRes[nCentralityBins];

  float centRange[2];
  for(int iCent = 0; iCent < nCentralityBins; iCent++){
    centRange[0] = arrayCentralityBinning[iCent];
    centRange[1] = arrayCentralityBinning[iCent+1];
    Get_PtResponseMatrix_Fluctuations(H2D_jetPtResponseMatrix_fluctuations[iCent], iDataset, iRadius, centRange);
    Get_PtResponseMatrix_Fluctuations_FullResForCheck(H2D_jetPtResponseMatrix_fluctuations_fullRes[iCent], iDataset, iRadius, centRange);
  }
  TString* pdfName = new TString("background_responseMatrices_"+jetType[iJetType]+"_"+Datasets[iDataset]+"_centralityComp");
  TString* pdfNameFullRes = new TString("background_responseMatrices_"+jetType[iJetType]+"_"+Datasets[iDataset]+"_centralityComp_FullRes");

  TString textContext(contextDatasetComp("")); // placeholder
  
  TString centralityLegend[nCentralityBins]; CentralityLegend(centralityLegend, arrayCentralityBinning, nCentralityBins);

  Draw_TH2_Histograms(H2D_jetPtResponseMatrix_fluctuations, centralityLegend, nCentralityBins, textContext, pdfName, texPtJetRecX, texPtJetGenX, texCollisionDataInfo, "logz");
  Draw_TH2_Histograms(H2D_jetPtResponseMatrix_fluctuations_fullRes, centralityLegend, nCentralityBins, textContext, pdfNameFullRes, texPtJetRecX, texPtJetGenX, texCollisionDataInfo, "logz");
}

void Draw_Pt_spectrum_unfolded_FluctResponseOnly(int iDataset, int iRadius, const char options[]) {

  TH1D* H1D_jetPt_unfolded[nCentralityBins];
  float centRange[2];
  for(int iCent = 0; iCent < nCentralityBins; iCent++){
    centRange[0] = arrayCentralityBinning[iCent];
    centRange[1] = arrayCentralityBinning[iCent+1];

    Get_Pt_spectrum_unfolded(H1D_jetPt_unfolded[iCent], iDataset, iRadius, centRange, options);
  }


  TString* pdfName = new TString("jet_"+jetType[iJetType]+"_"+jetLevel[iJetLevel]+"_"+Datasets[iDataset]+"_R="+Form("%.1f", arrayRadius[iRadius])+"_Pt_unfolded");

  TString textContext(contextDatasetComp("")); // placeholder
  
  TString centralityLegend[nCentralityBins]; CentralityLegend(centralityLegend, arrayCentralityBinning, nCentralityBins);

  TString* yAxisLabel;
  if (strstr(options, "evtNorm") != NULL) {
    yAxisLabel = texJetPtYield_EventNorm;
  }
  if (strstr(options, "entriesNorm") != NULL) {
    yAxisLabel = texJetPtYield_EntriesNorm;
  }
  Draw_TH1_Histograms_in_one(H1D_jetPt_unfolded, centralityLegend, nCentralityBins, textContext, pdfName, texPtX, yAxisLabel, texCollisionDataInfo, "logy");
}

// // // PhD fuction for V0s
// void Get_V0_RawYield_vsPt_FeeddownCorrected_withUnfolding(TH1D* &hRawYield_vsPt, TH3D* H3D_detectedV0s_DataLike, TH1I* H1I_SelectedEventCount, TFile *file_O2Analysis, double* SignalMean, double* SignalStandardDev, int ipart, double* pTbins, int nbinpT, int ibinXaxisCut_low,  int ibinXaxisCut_high, int warning_cutArry_ID, double SideBandSizeMultiplierModifier, int SignalExtractionType, int InvMassRebinFactor) {
//   float SelectedEventCount = H1I_SelectedEventCount->GetEntries();

//   TH1D* truth = (TH1D*)file_O2Analysis->Get("lambdakzero-particle-count-mc/h"+NamePart[ipart]+"Count_PtDiff");
//   TH1D* truth_rebinned = (TH1D*)truth->Rebin(nbinpT,"truth_rebinned",pTbins);
  
//   Get_RawYield_vsPt_FeeddownCorrected_RawCount(hRawYield_vsPt, H3D_detectedV0s_DataLike, H1I_SelectedEventCount, file_O2Analysis, SignalMean, SignalStandardDev, ipart, pTbins, nbinpT, ibinXaxisCut_low, ibinXaxisCut_high, warning_cutArry_ID, SideBandSizeMultiplierModifier, SignalExtractionType, InvMassRebinFactor);
//   TH2D *H2D_PtResponseMatrix = new TH2D("H2D_PtResponseMatrix", "H2D_PtResponseMatrix", nbinpT, pTbins, nbinpT, pTbins);
//   Get_PtResponseMatrix(H2D_PtResponseMatrix, file_O2Analysis, pTbins, nbinpT);

//   // Get_RawYield_vsPt_FeeddownCorrected(hRawYield_vsPt, H3D_detectedV0s_DataLike, H1I_SelectedEventCount, file_O2Analysis, SignalMean, SignalStandardDev, ipart, pTbins, nbinpT, ibinXaxisCut_low, ibinXaxisCut_high, warning_cutArry_ID, SideBandSizeMultiplierModifier, SignalExtractionType, InvMassRebinFactor);
//   // TH2D *H2D_PtResponseMatrix_Density = new TH2D("H2D_PtResponseMatrix", "H2D_PtResponseMatrix", nbinpT, pTbins, nbinpT, pTbins);
//   // Get_PtResponseMatrix_Density(H2D_PtResponseMatrix, file_O2Analysis, pTbins, nbinpT);
//   // for(int ibinPt = 1; ibinPt <= nbinpT; ibinPt++){
//   //   double dpT = truth_rebinned->GetXaxis()->GetBinWidth(ibinPt);
//   //   double drapidity = 1.5; //1.5
//   //   double d2N_dpTdy = truth_rebinned->GetBinContent(ibinPt) *1./dpT*1./drapidity;
//   //   truth_rebinned->SetBinContent(ibinPt,d2N_dpTdy);
//   //   truth_rebinned->SetBinError(ibinPt,truth_rebinned->GetBinError(ibinPt) *1./dpT*1./drapidity);  // error on d2N_dpTdy 
//   // }    
//   // truth_rebinned->Scale(1./SelectedEventCount);

//   TH1D* measured = (TH1D*)hRawYield_vsPt->Clone("measured");

//   // TH1D* projX = (TH1D*)H2D_PtResponseMatrix->ProjectionX("projX",0,-1);
//   // TH1D* projY = (TH1D*)H2D_PtResponseMatrix->ProjectionY("projY",0,-1);
//   cout << "--------------------------------------------------------------------- PRE UNFOLDING TEST" << endl;
//   for(int ibinPt = 1; ibinPt <= nbinpT; ibinPt++){
//     cout << "measured(" << ibinPt << ") = " << measured->GetBinContent(ibinPt) << endl;
//     cout << "truth_rebinned(" << ibinPt << ") = " << truth_rebinned->GetBinContent(ibinPt) << endl;
//     for(int jbinPt = 1; jbinPt <= nbinpT; jbinPt++){
//     cout << "H2D_PtResponseMatrix(" << ibinPt << "," << jbinPt << ") = " << H2D_PtResponseMatrix->GetBinContent(ibinPt,jbinPt) << endl;
//     }
//   }

//   RooUnfoldResponse* response = new RooUnfoldResponse(measured, truth_rebinned, H2D_PtResponseMatrix);
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

//   // TUnfold unfold(H2D_PtResponseMatrix, TUnfold::kHistMapOutputHoriz, TUnfold::kRegModeNone); // kRegModeNone regularisation method choice, here none, will see what's best later
//   // // set input distribution and bias scale (=0)
//   // cout << "AMIERIC DEBUG unfolding 4" << endl;
//   // if(unfold.SetInput(hRawYield_vsPt, 0.0)>=10000) {
//   //   std::cout<<"Unfolding result may be wrong\n";
//   // }
//   // cout << "AMIERIC DEBUG unfolding 5" << endl;

//   // // do the unfolding here
//   // double tauMin=0.0;
//   // double tauMax=0.0;
//   // int nScan=30;
//   // int iBest;
//   // TSpline *logTauX,*logTauY;
//   // TGraph *lCurve;
//   // // this method scans the parameter tau
//   // // finally, the unfolding is done for the "best" choice of tau
//   // iBest=unfold.ScanLcurve(nScan,tauMin,tauMax,&lCurve,&logTauX,&logTauY);
//   // std::cout<<"tau="<<unfold.GetTau()<<"\n";
//   // std::cout<<"chi**2="<<unfold.GetChi2A()<<"+"<<unfold.GetChi2L()
//   //          <<" / "<<unfold.GetNdf()<<"\n";
//   // // save point corresponding to the kink in the L curve as TGraph
//   // double t[1],x[1],y[1];
//   // logTauX->GetKnot(iBest,t[0],x[0]);
//   // logTauY->GetKnot(iBest,t[0],y[0]);
//   // TGraph *bestLcurve=new TGraph(1,x,y);
//   // TGraph *bestLogTauX=new TGraph(1,t,x);

//   // //============================================================
//   // // extract unfolding results into histograms
//   // // set up a bin map, excluding underflow and overflow bins
//   // // the binMap relates the the output of the unfolding to the final
//   // // histogram bins
//   // int *binMap=new int[nbinpT+2];
//   // for(int i=1;i<=nbinpT;i++) binMap[i]=i;
//   // binMap[0]=-1;
//   // binMap[nbinpT+1]=-1;
//   // TH1D *histPtUnfold=new TH1D("Unfolded",";pT(gen)",nbinpT,pTbins);
//   // unfold.GetOutput(histPtUnfold,binMap);
//   // TH1D *histPtDetFold=new TH1D("FoldedBack","pT(det)",nbinpT,pTbins);
//   // unfold.GetFoldedOutput(histPtDetFold); // Aimeric: probably for comparison? some kind of closure test
//   // // store global correlation coefficients
//   // TH1D *histRhoi=new TH1D("rho_I","mass",nbinpT,pTbins);
//   // unfold.GetRhoI(histRhoi,binMap);
//   // delete[] binMap;
//   // binMap=0;
// }
