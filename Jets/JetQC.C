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
TString contextDataset(int iRun);
TString contextPtRange(float* PtRange);
TString contextJetRadius(float jetRadius);

//////////// QC plot functions
// Radius comparison
void Draw_Pt_RadiusComparison(int iRun, float* etaRange);
void Draw_Eta_RadiusComparison(int iRun, float* PtRange);
void Draw_Phi_RadiusComparison(int iRun, float* PtRange);
void Draw_NTracks_RadiusComparison_withPtRange(int iRun, float* PtRange);
void Draw_LeadingTrackPt_vs_JetPt_RadiusComparison(int iRun);
void Draw_JetArea_vs_JetPt_RadiusComparison(int iRun, float* PtRange);
void Draw_JetNTracks_vs_JetPt_RadiusComparison(int iRun, float* PtRange);
void Draw_Pt_ratio_etaNeg_etaPos_RadiusComparison(int iRun, float* etaRange);

// Run comparison
void Draw_Pt_RunComparison(float jetRadius, float* etaRange);
void Draw_Eta_RunComparison(float jetRadius, float* PtRange);
void Draw_Phi_RunComparison(float jetRadius, float* PtRange);
void Draw_Pt_ratio_etaNeg_etaPos_RunComparison(float jetRadius, float* etaRange);


// Options to be set:
//////// -------- Full Analysis -------- ////////
TString* texCollisionDataInfo = new TString("Pb-Pb #sqrt{#it{s}} = 5.36 TeV"); 
const Int_t nRuns = 9;
const TString Runs[nRuns] = {"LHC23zzi_cpass0", "LHC23zzf_cpass0", "LHC23zzg_cpass0", "LHC23zx_cpass0", "LHC23zy_cpass0", "LHC23zz_cpass0", "LHC23zza_cpass0", "LHC23zzb_cpass0", "LHC23zze_cpass0"};
TFile* file_O2Analysis_list[nRuns] = {new TFile("Datasets/"+Runs[0]+"/AnalysisResults.root"),
                                      new TFile("Datasets/"+Runs[1]+"/AnalysisResults.root"),
                                      new TFile("Datasets/"+Runs[2]+"/AnalysisResults.root"),
                                      new TFile("Datasets/"+Runs[3]+"/AnalysisResults.root"),
                                      new TFile("Datasets/"+Runs[4]+"/AnalysisResults.root"),
                                      new TFile("Datasets/"+Runs[5]+"/AnalysisResults.root"),
                                      new TFile("Datasets/"+Runs[6]+"/AnalysisResults.root"),
                                      new TFile("Datasets/"+Runs[7]+"/AnalysisResults.root"),
                                      new TFile("Datasets/"+Runs[8]+"/AnalysisResults.root")
                                      };
const TString analysisWorkflow = "jet-finder-charged-qa";
// //////// -------- Flat Phi Periods -------- ////////
// TString* texCollisionDataInfo = new TString("Pb-Pb #sqrt{#it{s}} = 5.36 TeV"); 
// const Int_t nRuns = 3;
// const TString Runs[nRuns] = {"LHC23zzf_cpass0", "LHC23zzg_cpass0", "LHC23zzi_cpass0"};
// TFile* file_O2Analysis_list[nRuns] = {new TFile("Datasets/"+Runs[0]+"/AnalysisResults.root"),
//                                       new TFile("Datasets/"+Runs[1]+"/AnalysisResults.root"),
//                                       new TFile("Datasets/"+Runs[2]+"/AnalysisResults.root")
//                                       };
// const TString analysisWorkflow = "jet-finder-charged-qa";
// //////// -------- Cpass2 test -------- ////////
// TString* texCollisionDataInfo = new TString("Pb-Pb #sqrt{#it{s}} = 5.36 TeV"); 
// const Int_t nRuns = 1;
// const TString Runs[nRuns] = {"LHC23zzh_cpass2"};
// TFile* file_O2Analysis_list[nRuns] = {new TFile("Datasets/"+Runs[0]+"/AnalysisResults.root")
//                                       };
// const TString analysisWorkflow = "jet-finder-charged-qa";
// //////// -------- Cpass1 test -------- ////////
// TString* texCollisionDataInfo = new TString("Pb-Pb #sqrt{#it{s}} = 5.36 TeV"); 
// const Int_t nRuns = 1;
// const TString Runs[nRuns] = {"LHC23zzh_cpass1"};
// TFile* file_O2Analysis_list[nRuns] = {new TFile("Datasets/"+Runs[0]+"/AnalysisResults.root")
//                                       };
// const TString analysisWorkflow = "jet-finder-charged-qa";
////// -------- pre renaming of jet-finder-qa to jet-finder-charged-qa -------- ////////
// TString* texCollisionDataInfo = new TString("Pb-Pb #sqrt{#it{s}} = 5.36 TeV"); 
// const Int_t nRuns = 4;
// const TString Runs[nRuns] = {"LHC23zx_cpass0", "LHC23zy_cpass0", "LHC23zz_cpass0", "LHC23zzb_cpass0"};
// TFile* file_O2Analysis_list[nRuns] = {new TFile("Datasets/"+Runs[0]+"/AnalysisResults.root"),
//                                       new TFile("Datasets/"+Runs[1]+"/AnalysisResults.root"),
//                                       new TFile("Datasets/"+Runs[2]+"/AnalysisResults.root"),
//                                       new TFile("Datasets/"+Runs[3]+"/AnalysisResults.root")
//                                       };
// const TString analysisWorkflow = "jet-finder-qa";

// Analysis settings
const Int_t nJetType = 3;
const TString jetType[nJetType] = {"charged", "neutral", "full"};
const Int_t nJetLevel = 3;
const TString jetLevel[nJetLevel] = {"data", "mcd", "mcp"};
const Int_t nRadius = 3;
const TString RadiusLegend[nRadius] = {"R = 0.2", "R = 0.4", "R = 0.6"};
Float_t arrayRadius[nRadius] = {0.2, 0.4, 0.6};
// const Int_t nRadius = 1;
// const TString RadiusLegend[nRadius] = {"R = 0.4"};
// Float_t arrayRadius[nRadius] = {0.4};

// Choice of jet type (charged, neutral, full) and level (data, detector level, particle level)
const Int_t iJetType = 0;
const Int_t iJetLevel = 0;

// Commonly used titles
TString* texRapidity = new TString("y");
TString* texPtX = new TString("#it{p}_{T} (GeV/#it{c})");
TString* texEtaX = new TString("#it{#eta}");
TString* texPhiX = new TString("#it{#phi} (rad)");
TString* texNTracksX = new TString("#it{N}_{tracks}");
TString* texLeadPt = new TString("#it{p}_{T,track}^{leading}");
TString* texJetArea = new TString("#it{A}_{jet}");
TString* texJetNTracks = new TString("#it{N}_{tracks}");

TString* texPtDifferentialYield = new TString("1/#it{N}_{ev} d#it{N}/d#it{p}_{T} (GeV/#it{c})^{-1}");

TString* texCount = new TString("count");

TString* texPtMeasured = new TString("#it{p}_{T}^{measured} (GeV/#it{c})");
TString* texPtMC = new TString("#it{p}_{T}^{MC} (GeV/#it{c})");

TString* texNSubJettinessRatio = new TString("#it{#tau_{2}/#tau_{1}}");
TString* texDeltaR = new TString("#it{#Delta R}");
TString* texdN_dsubratio = new TString("1/N^{ev} d#it{N}/d(#it{#tau_{2}/#tau_{1}})");
TString* texdN_dDeltaR = new TString("1/N^{ev} d#it{N}/d(#it{#Delta R})");

TString* texNormPtYield = new TString("1/#it{N}_{ev} d#it{N}_{jet}/d#it{p}_{T}");
TString* texNormEtaYield = new TString("1/#it{N}_{ev} d#it{N}_{jet}/d#it{#eta}");
TString* texNormPhiYield = new TString("1/#it{N}_{ev} d#it{N}_{jet}/d#it{#phi}");
TString* texNormNTracksYield = new TString("1/#it{N}_{ev} d#it{N}_{jet}/d#it{N}_{tracks}");


TString* texRatioRuns = new TString("Run / Run "+Runs[0]);
TString* texRatioEtaComparison = new TString("#it{N}_{jet}^{ #eta > 0} / #it{N}_{jet}^{ #eta < 0}");


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

  int iRun = 0;

  float etaRangeSym[2] = {-0.5, 0.5};
  float etaRangeNeg[2] = {-0.5, 0};
  float etaRangePos[2] = {0, 0.5};

  float jetRadiusForRunComp = 0.4;

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

    // Draw_Eta_RadiusComparison(iRun, ptRange);
    // Draw_Phi_RadiusComparison(iRun, ptRange);
    // Draw_NTracks_RadiusComparison_withPtRange(iRun, PtRangeZoom0);
    // Draw_NTracks_RadiusComparison_withPtRange(iRun, PtRangeZoom2030);
    // Draw_NTracks_RadiusComparison_withPtRange(iRun, PtRangeZoom3040);
    // Draw_NTracks_RadiusComparison_withPtRange(iRun, PtRangeZoom4050);
    // Draw_NTracks_RadiusComparison_withPtRange(iRun, PtRangeZoom5060);
    // Draw_NTracks_RadiusComparison_withPtRange(iRun, PtRangeZoom8090);
    // // Draw_LeadingTrackPt_vs_JetPt_RadiusComparison(iRun); leading pT not implemented yet
    Draw_JetArea_vs_JetPt_RadiusComparison(iRun, PtRangeZoom0);
    Draw_JetNTracks_vs_JetPt_RadiusComparison(iRun, PtRangeZoom0);
    Draw_JetArea_vs_JetPt_RadiusComparison(iRun, PtRangeZoom020);
    Draw_JetNTracks_vs_JetPt_RadiusComparison(iRun, PtRangeZoom020);

    // Draw_Eta_RunComparison(jetRadiusForRunComp, ptRange);
    // Draw_Phi_RunComparison(jetRadiusForRunComp, ptRange);
  }

    // Draw_Pt_RadiusComparison(iRun, etaRangeSym);
    // Draw_Pt_RadiusComparison(iRun, etaRangeNeg);
    // Draw_Pt_RadiusComparison(iRun, etaRangePos);
    // Draw_Pt_RunComparison(jetRadiusForRunComp, etaRangeSym);
    // Draw_Pt_ratio_etaNeg_etaPos_RadiusComparison(iRun, etaRangeSym);
    // Draw_Pt_ratio_etaNeg_etaPos_RunComparison(jetRadiusForRunComp, etaRangeSym);

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

TString contextDataset(int iRun){
  TString texDataset("2023 QC "+Runs[iRun]);
  return texDataset;
}

TString contextPtRange(float* PtRange){

  std::stringstream ss;
  ss << PtRange[0] << " < #it{p}_{T} < " << PtRange[1];
  TString textContext((TString)ss.str());
  // TString texDataset(Form("%.0f", PtRange[0])+" < #it{p}_{T} < "+Form("%.0f", PtRange[1]));
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

void Draw_Pt_RadiusComparison(int iRun, float* etaRange) {

  TH3D* H3D_jetRjetPtjetEta;
  TH1D* H1D_jetPt[nRadius];
  TH1D* H1D_jetPt_rebinned[nRadius];
  
  H3D_jetRjetPtjetEta = (TH3D*)file_O2Analysis_list[iRun]->Get(analysisWorkflow+"/h3_jet_r_jet_pt_jet_eta");
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
    ibinJetRadius = H3D_jetRjetPtjetEta->GetXaxis()->FindBin(arrayRadius[iRadius]);

    H1D_jetPt[iRadius] = (TH1D*)H3D_jetRjetPtjetEta->ProjectionY("jetPt_"+RadiusLegend[iRadius]+Form("%.1f", EtaCutLow)+"<eta<"+Form("%.1f", EtaCutHigh), ibinJetRadius, ibinJetRadius, ibinEta_low, ibinEta_high);
    H1D_jetPt_rebinned[iRadius] = (TH1D*)H1D_jetPt[iRadius]->Rebin(1.,"H1D_jetPt_rebinned_"+RadiusLegend[iRadius]);

    // NormaliseYieldToNJets(H1D_jetPt_rebinned[iRadius]);
    NormaliseYieldToNEvents(H1D_jetPt_rebinned[iRadius], GetNEventsSel8(file_O2Analysis_list[iRun]));
  }

  TString* pdfName = new TString("jet_"+jetType[iJetType]+"_"+jetLevel[iJetLevel]+"_"+Runs[iRun]+"_Pt_@eta["+Form("%.1f", EtaCutLow)+","+Form("%.1f", EtaCutHigh)+"]");

  std::stringstream ss;
  ss << EtaCutLow << " < #eta < " << EtaCutHigh;
  TString textContext("#splitline{"+contextDataset(iRun)+"}{"+(TString)ss.str()+"}");

  Draw_TH1_Histograms_in_one(H1D_jetPt_rebinned, RadiusLegend, nRadius, textContext, pdfName, texPtX, texNormPtYield, texCollisionDataInfo, "logy");
}

void Draw_Eta_RadiusComparison(int iRun, float* PtRange) {

  TH3D* H3D_jetRjetPtjetEta;
  TH1D* H1D_jetEta[nRadius];
  TH1D* H1D_jetEta_rebinned[nRadius];
  
  H3D_jetRjetPtjetEta = (TH3D*)file_O2Analysis_list[iRun]->Get(analysisWorkflow+"/h3_jet_r_jet_pt_jet_eta");
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
    ibinJetRadius = H3D_jetRjetPtjetEta->GetXaxis()->FindBin(arrayRadius[iRadius]);
    H1D_jetEta[iRadius] = (TH1D*)H3D_jetRjetPtjetEta->ProjectionZ("jetEta_"+RadiusLegend[iRadius]+Form("%.1f", PtCutLow)+"<pt<"+Form("%.1f", PtCutHigh), ibinJetRadius, ibinJetRadius, ibinPt_low, ibinPt_high);
    H1D_jetEta_rebinned[iRadius] = (TH1D*)H1D_jetEta[iRadius]->Rebin(1.,"H1D_jetEta_rebinned_"+RadiusLegend[iRadius]);

    // NormaliseYieldToNJets(H1D_jetEta_rebinned[iRadius]);
    NormaliseYieldToNEvents(H1D_jetEta_rebinned[iRadius], GetNEventsSel8(file_O2Analysis_list[iRun]));
  }

  TString* pdfName = new TString("jet_"+jetType[iJetType]+"_"+jetLevel[iJetLevel]+"_"+Runs[iRun]+"_Eta_@pt["+Form("%03.0f", PtCutLow)+","+Form("%03.0f", PtCutHigh)+"]");

  TString textContext("#splitline{"+contextDataset(iRun)+"}{"+contextPtRange(PtRange)+"}");

  Draw_TH1_Histograms_in_one(H1D_jetEta_rebinned, RadiusLegend, nRadius, textContext, pdfName, texEtaX, texNormEtaYield, texCollisionDataInfo, "");
}

void Draw_Phi_RadiusComparison(int iRun, float* PtRange) {

  TH3D* H3D_jetRjetPtjetPhi;
  TH1D* H1D_jetPhi[nRadius];
  TH1D* H1D_jetPhi_rebinned[nRadius];
  
  H3D_jetRjetPtjetPhi = (TH3D*)file_O2Analysis_list[iRun]->Get(analysisWorkflow+"/h3_jet_r_jet_pt_jet_phi");
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
    ibinJetRadius = H3D_jetRjetPtjetPhi->GetXaxis()->FindBin(arrayRadius[iRadius]);
    H1D_jetPhi[iRadius] = (TH1D*)H3D_jetRjetPtjetPhi->ProjectionZ("jetPhi_"+RadiusLegend[iRadius]+Form("%.1f", PtCutLow)+"<pt<"+Form("%.1f", PtCutHigh), ibinJetRadius, ibinJetRadius, ibinPt_low, ibinPt_high);
    H1D_jetPhi_rebinned[iRadius] = (TH1D*)H1D_jetPhi[iRadius]->Rebin(1.,"H1D_jetPhi_rebinned_"+RadiusLegend[iRadius]);

    // NormaliseYieldToNJets(H1D_jetPhi_rebinned[iRadius]);
    NormaliseYieldToNEvents(H1D_jetPhi_rebinned[iRadius], GetNEventsSel8(file_O2Analysis_list[iRun]));
  }
 
  TString* pdfName = new TString("jet_"+jetType[iJetType]+"_"+jetLevel[iJetLevel]+"_"+Runs[iRun]+"_Phi_@pt["+Form("%03.0f", PtCutLow)+","+Form("%03.0f", PtCutHigh)+"]");

  TString textContext("#splitline{"+contextDataset(iRun)+"}{"+contextPtRange(PtRange)+"}");

  Draw_TH1_Histograms_in_one(H1D_jetPhi_rebinned, RadiusLegend, nRadius, textContext, pdfName, texPhiX, texNormPhiYield, texCollisionDataInfo, "");
}

void Draw_NTracks_RadiusComparison_withPtRange(int iRun, float* PtRange) {

  TH3D* H3D_jetRjetPtjetNTracks;
  TH1D* H1D_jetNTracks[nRadius];
  TH1D* H1D_jetNTracks_rebinned[nRadius];
  
  H3D_jetRjetPtjetNTracks = (TH3D*)file_O2Analysis_list[iRun]->Get(analysisWorkflow+"/h3_jet_r_jet_pt_jet_ntracks");
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
    ibinJetRadius = H3D_jetRjetPtjetNTracks->GetXaxis()->FindBin(arrayRadius[iRadius]);
    H1D_jetNTracks[iRadius] = (TH1D*)H3D_jetRjetPtjetNTracks->ProjectionZ("jetNTracks_"+RadiusLegend[iRadius], ibinJetRadius, ibinJetRadius, ibinPt_low, ibinPt_high);
    H1D_jetNTracks_rebinned[iRadius] = (TH1D*)H1D_jetNTracks[iRadius]->Rebin(1.,"H1D_jetNTracks_rebinned_"+RadiusLegend[iRadius]);

    // NormaliseYieldToNJets(H1D_jetNTracks_rebinned[iRadius]);
    NormaliseYieldToNEvents(H1D_jetNTracks_rebinned[iRadius], GetNEventsSel8(file_O2Analysis_list[iRun]));
  }

  TString* pdfName = new TString("jet_"+jetType[iJetType]+"_"+jetLevel[iJetLevel]+"_"+Runs[iRun]+"_NTracks_@pT["+Form("%03.0f", PtCutLow)+","+Form("%03.0f", PtCutHigh)+"]");

  TString textContext("#splitline{"+contextDataset(iRun)+"}{"+contextPtRange(PtRange)+"}");

  Draw_TH1_Histograms_in_one(H1D_jetNTracks_rebinned, RadiusLegend, nRadius, textContext, pdfName, texNTracksX, texNormNTracksYield, texCollisionDataInfo, "");
}


void Draw_LeadingTrackPt_vs_JetPt_RadiusComparison(int iRun) {

  TH3D* H3D_jetRjetPtjetLeadPt;
  TH2D* H2D_jetLeadPt[nRadius];
  TH2D* H2D_jetLeadPt_rebinned[nRadius];
  
  H3D_jetRjetPtjetLeadPt = (TH3D*)file_O2Analysis_list[iRun]->Get(analysisWorkflow+"/h3_jet_r_jet_pt_leadingtrack_pt");
  H3D_jetRjetPtjetLeadPt->Sumw2();

  for(int iRadius = 0; iRadius < nRadius; iRadius++){
    H3D_jetRjetPtjetLeadPt->GetXaxis()->SetRange(iRadius+1,iRadius+1);
    H2D_jetLeadPt[iRadius] = (TH2D*)H3D_jetRjetPtjetLeadPt->Project3D(RadiusLegend[iRadius]+"_zy");
    // H2D_jetLeadPt_rebinned[iRadius] = (TH2D*)H2D_jetLeadPt[iRadius]->Rebin(1.,"H1D_jetLeadPt_rebinned_"+RadiusLegend[iRadius]);
  }

  TString* pdfName = new TString("jet_"+jetType[iJetType]+"_"+jetLevel[iJetLevel]+"_"+Runs[iRun]+"_LeadPt");
  // std::stringstream ss;
  // ss << PtCutLow << " < #it{p}_{T} < " << PtCutHigh;
  // TString textContext(ss.str());
  TString textContext(contextDataset(iRun));

  Draw_TH2_Histograms(H2D_jetLeadPt, RadiusLegend, nRadius, textContext, pdfName, texPtX, texLeadPt, texCollisionDataInfo, "");
}

void Draw_JetArea_vs_JetPt_RadiusComparison(int iRun, float* PtRange) {

  TH3D* H3D_jetRjetPtjetArea;
  TH2D* H2D_jetArea[nRadius];
  TH2D* H2D_jetArea_rebinned[nRadius];
  
  H3D_jetRjetPtjetArea = (TH3D*)file_O2Analysis_list[iRun]->Get(analysisWorkflow+"/h3_jet_r_jet_pt_jet_area");
  H3D_jetRjetPtjetArea->Sumw2();

  float PtCutLow = PtRange[0];
  float PtCutHigh = PtRange[1];
  int ibinPt_low = H3D_jetRjetPtjetArea->GetYaxis()->FindBin(PtCutLow);
  int ibinPt_high = H3D_jetRjetPtjetArea->GetYaxis()->FindBin(PtCutHigh);

  for(int iRadius = 0; iRadius < nRadius; iRadius++){
    H3D_jetRjetPtjetArea->GetXaxis()->SetRange(iRadius+1,iRadius+1);
    H3D_jetRjetPtjetArea->GetYaxis()->SetRange(ibinPt_low,ibinPt_high);
    H2D_jetArea[iRadius] = (TH2D*)H3D_jetRjetPtjetArea->Project3D(RadiusLegend[iRadius]+"_jetArea_zy");
    // H2D_jetArea_rebinned[iRadius] = (TH2D*)H2D_jetArea[iRadius]->Rebin(1.,"H1D_jetArea_rebinned_"+RadiusLegend[iRadius]);
  }

  TString* pdfName = new TString("jet_"+jetType[iJetType]+"_"+jetLevel[iJetLevel]+"_"+Runs[iRun]+"_JetArea-vs-Pt_@pT["+Form("%03.0f", PtCutLow)+","+Form("%03.0f", PtCutHigh)+"]");
  // std::stringstream ss;
  // ss << PtCutLow << " < #it{p}_{T} < " << PtCutHigh;
  // TString textContext(ss.str());
  TString textContext(contextDataset(iRun));

  Draw_TH2_Histograms(H2D_jetArea, RadiusLegend, nRadius, textContext, pdfName, texPtX, texJetArea, texCollisionDataInfo, "");
  TString* pdfNamelogz = new TString(*pdfName + "_logz");
  Draw_TH2_Histograms(H2D_jetArea, RadiusLegend, nRadius, textContext, pdfNamelogz, texPtX, texJetArea, texCollisionDataInfo, "logz");
  TString* pdfNamelogy = new TString(*pdfName + "_logx");
  Draw_TH2_Histograms(H2D_jetArea, RadiusLegend, nRadius, textContext, pdfNamelogy, texPtX, texJetArea, texCollisionDataInfo, "logx");
  TString* pdfNamelogyz = new TString(*pdfName + "_logxz");
  Draw_TH2_Histograms(H2D_jetArea, RadiusLegend, nRadius, textContext, pdfNamelogyz, texPtX, texJetArea, texCollisionDataInfo, "logxlogz");
}

void Draw_JetNTracks_vs_JetPt_RadiusComparison(int iRun, float* PtRange) {

  TH3D* H3D_jetRjetPtjetNTracks;
  TH2D* H2D_jetNTracks[nRadius];
  TH2D* H2D_jetNTracks_rebinned[nRadius];
  
  H3D_jetRjetPtjetNTracks = (TH3D*)file_O2Analysis_list[iRun]->Get(analysisWorkflow+"/h3_jet_r_jet_pt_jet_ntracks");
  H3D_jetRjetPtjetNTracks->Sumw2();

  float PtCutLow = PtRange[0];
  float PtCutHigh = PtRange[1];
  int ibinPt_low = H3D_jetRjetPtjetNTracks->GetYaxis()->FindBin(PtCutLow);
  int ibinPt_high = H3D_jetRjetPtjetNTracks->GetYaxis()->FindBin(PtCutHigh);

  for(int iRadius = 0; iRadius < nRadius; iRadius++){
    H3D_jetRjetPtjetNTracks->GetXaxis()->SetRange(iRadius+1,iRadius+1);
    H3D_jetRjetPtjetNTracks->GetYaxis()->SetRange(ibinPt_low,ibinPt_high);
    H2D_jetNTracks[iRadius] = (TH2D*)H3D_jetRjetPtjetNTracks->Project3D(RadiusLegend[iRadius]+"_jetNTracks_zy");
    // H2D_jetNTracks_rebinned[iRadius] = (TH2D*)H2D_jetNTracks[iRadius]->Rebin(1.,"H1D_jetNTracks_rebinned_"+RadiusLegend[iRadius]);
  }

  TString* pdfName = new TString("jet_"+jetType[iJetType]+"_"+jetLevel[iJetLevel]+"_"+Runs[iRun]+"_JetNTracks-vs-Pt_@pT["+Form("%03.0f", PtCutLow)+","+Form("%03.0f", PtCutHigh)+"]");

  TString textContext(contextDataset(iRun));

  Draw_TH2_Histograms(H2D_jetNTracks, RadiusLegend, nRadius, textContext, pdfName, texPtX, texJetNTracks, texCollisionDataInfo, "");
  TString* pdfNamelogz = new TString(*pdfName + "_logz");
  Draw_TH2_Histograms(H2D_jetNTracks, RadiusLegend, nRadius, textContext, pdfNamelogz, texPtX, texJetNTracks, texCollisionDataInfo, "logz");
}

void Draw_Pt_RunComparison(float jetRadius, float* etaRange) {

  TH3D* H3D_jetRjetPtjetEta[nRuns];
  TH1D* H1D_jetPt[nRuns];
  TH1D* H1D_jetPt_rebinned[nRuns];
  
  TH1D* H1D_jetPt_rebinned_ratios[nRuns];

  float EtaCutLow = etaRange[0];
  float EtaCutHigh = etaRange[1];
  int ibinJetRadius = 0;

  bool divideSuccess = false;

  for(int iRun = 0; iRun < nRuns; iRun++){

    H3D_jetRjetPtjetEta[iRun] = (TH3D*)file_O2Analysis_list[iRun]->Get(analysisWorkflow+"/h3_jet_r_jet_pt_jet_eta");

    int ibinEta_low = H3D_jetRjetPtjetEta[iRun]->GetZaxis()->FindBin(EtaCutLow);
    int ibinEta_high = H3D_jetRjetPtjetEta[iRun]->GetZaxis()->FindBin(EtaCutHigh);
    if (ibinEta_low == 0) 
      cout << "WARNING: Pt_RunComparison is counting the underflow with the chosen etaRange" << endl;
    if (ibinEta_high == H3D_jetRjetPtjetEta[iRun]->GetZaxis()->GetNbins()+1) 
      cout << "WARNING: Pt_RunComparison is counting the overflow with the chosen etaRange" << endl;
    ibinJetRadius = H3D_jetRjetPtjetEta[iRun]->GetXaxis()->FindBin(jetRadius);
    
    H1D_jetPt[iRun] = (TH1D*)H3D_jetRjetPtjetEta[iRun]->ProjectionY("jetPt_"+Runs[iRun]+Form("%.1f", EtaCutLow)+"<eta<"+Form("%.1f", EtaCutHigh), ibinJetRadius, ibinJetRadius, ibinEta_low, ibinEta_high);
    H1D_jetPt_rebinned[iRun] = (TH1D*)H1D_jetPt[iRun]->Rebin(1.,"H1D_jetPt_rebinned_"+Runs[iRun]);

    // NormaliseYieldToNJets(H1D_jetPt_rebinned[iRun]);
    NormaliseYieldToNEvents(H1D_jetPt_rebinned[iRun], GetNEventsSel8(file_O2Analysis_list[iRun]));

    H1D_jetPt_rebinned_ratios[iRun] = (TH1D*)H1D_jetPt_rebinned[iRun]->Clone("H1D_jetPt_rebinned_ratios"+Runs[iRun]);
    H1D_jetPt_rebinned_ratios[iRun]->Reset("M");
    divideSuccess = H1D_jetPt_rebinned_ratios[iRun]->Divide(H1D_jetPt_rebinned[iRun], H1D_jetPt_rebinned[0]);
  }

  TString* pdfName = new TString("jet_"+jetType[iJetType]+"_"+jetLevel[iJetLevel]+"_RunComp_R="+Form("%.1f", jetRadius)+"_Pt_@eta["+Form("%.1f", EtaCutLow)+","+Form("%.1f", EtaCutHigh)+"]");
  TString* pdfName_ratio = new TString("jet_"+jetType[iJetType]+"_"+jetLevel[iJetLevel]+"_RunComp_R="+Form("%.1f", jetRadius)+"_Pt_@eta["+Form("%.1f", EtaCutLow)+","+Form("%.1f", EtaCutHigh)+"]_ratio");

  std::stringstream ss;
  ss << EtaCutLow << " < #eta < " << EtaCutHigh;
  TString textContext("#splitline{"+contextJetRadius(jetRadius)+"}{"+(TString)ss.str()+"}");

  Draw_TH1_Histograms_in_one(H1D_jetPt_rebinned, Runs, nRuns, textContext, pdfName, texPtX, texNormPtYield, texCollisionDataInfo, "logy");
  if (divideSuccess == true) {
    Draw_TH1_Histograms_in_one(H1D_jetPt_rebinned_ratios, Runs, nRuns, textContext, pdfName_ratio, texPtX, texRatioRuns, texCollisionDataInfo, "standardratio,avoidFirst");
  }
  else {
    cout << "Divide failed in Draw_Pt_RunComparison" << endl;
  }
}

void Draw_Eta_RunComparison(float jetRadius, float* PtRange) {

  TH3D* H3D_jetRjetPtjetEta[nRuns];
  TH1D* H1D_jetEta[nRuns];
  TH1D* H1D_jetEta_rebinned[nRuns];
  
  TH1D* H1D_jetEta_rebinned_ratios[nRuns];

  float PtCutLow = PtRange[0];
  float PtCutHigh = PtRange[1];
  int ibinJetRadius = 0;

  bool divideSuccess = false;

  for(int iRun = 0; iRun < nRuns; iRun++){

    H3D_jetRjetPtjetEta[iRun] = (TH3D*)file_O2Analysis_list[iRun]->Get(analysisWorkflow+"/h3_jet_r_jet_pt_jet_eta");

    int ibinPt_low = H3D_jetRjetPtjetEta[iRun]->GetYaxis()->FindBin(PtCutLow);
    int ibinPt_high = H3D_jetRjetPtjetEta[iRun]->GetYaxis()->FindBin(PtCutHigh);
    if (ibinPt_low == 0) 
      cout << "WARNING: Eta_RunComparison is counting the underflow with the chosen PtRange" << endl;
    if (ibinPt_high == H3D_jetRjetPtjetEta[iRun]->GetYaxis()->GetNbins()+1) 
      cout << "WARNING: Eta_RunComparison is counting the overflow with the chosen PtRange" << endl;
    ibinJetRadius = H3D_jetRjetPtjetEta[iRun]->GetXaxis()->FindBin(jetRadius);

    H1D_jetEta[iRun] = (TH1D*)H3D_jetRjetPtjetEta[iRun]->ProjectionZ("jetEta_"+Runs[iRun]+Form("%.1f", PtCutLow)+"<pt<"+Form("%.1f", PtCutHigh), ibinJetRadius, ibinJetRadius, ibinPt_low, ibinPt_high);
    H1D_jetEta_rebinned[iRun] = (TH1D*)H1D_jetEta[iRun]->Rebin(1.,"H1D_jetEta_rebinned_"+Runs[iRun]);

    // NormaliseYieldToNJets(H1D_jetEta_rebinned[iRun]);
    NormaliseYieldToNEvents(H1D_jetEta_rebinned[iRun], GetNEventsSel8(file_O2Analysis_list[iRun]));

    H1D_jetEta_rebinned_ratios[iRun] = (TH1D*)H1D_jetEta_rebinned[iRun]->Clone("H1D_jetEta_rebinned_ratios"+Runs[iRun]);
    H1D_jetEta_rebinned_ratios[iRun]->Reset("M");
    divideSuccess = H1D_jetEta_rebinned_ratios[iRun]->Divide(H1D_jetEta_rebinned[iRun], H1D_jetEta_rebinned[0]);
  }

  TString* pdfName = new TString("jet_"+jetType[iJetType]+"_"+jetLevel[iJetLevel]+"_RunComp_R="+Form("%.1f", jetRadius)+"_Eta_@pT["+Form("%03.0f", PtCutLow)+","+Form("%03.0f", PtCutHigh)+"]");
  TString* pdfName_ratio = new TString("jet_"+jetType[iJetType]+"_"+jetLevel[iJetLevel]+"_RunComp_R="+Form("%.1f", jetRadius)+"_Eta_@pT["+Form("%03.0f", PtCutLow)+","+Form("%03.0f", PtCutHigh)+"]_ratio");

  TString textContext("#splitline{"+contextJetRadius(jetRadius)+"}{"+contextPtRange(PtRange)+"}");

  Draw_TH1_Histograms_in_one(H1D_jetEta_rebinned, Runs, nRuns, textContext, pdfName, texEtaX, texNormEtaYield, texCollisionDataInfo, "");

  if (divideSuccess == true) {
    Draw_TH1_Histograms_in_one(H1D_jetEta_rebinned_ratios, Runs, nRuns, textContext, pdfName_ratio, texEtaX, texRatioRuns, texCollisionDataInfo, "autoratio,avoidFirst");
  }
  else {
    cout << "Divide failed in Draw_Eta_RunComparison" << endl;
  }
}

void Draw_Phi_RunComparison(float jetRadius, float* PtRange) {

  TH3D* H3D_jetRjetPtjetPhi[nRuns];
  TH1D* H1D_jetPhi[nRuns];
  TH1D* H1D_jetPhi_rebinned[nRuns];
  
  TH1D* H1D_jetPhi_rebinned_ratios[nRuns];

  float PtCutLow = PtRange[0];
  float PtCutHigh = PtRange[1];
  int ibinJetRadius = 0;

  bool divideSuccess = false;

  for(int iRun = 0; iRun < nRuns; iRun++){

    H3D_jetRjetPtjetPhi[iRun] = (TH3D*)file_O2Analysis_list[iRun]->Get(analysisWorkflow+"/h3_jet_r_jet_pt_jet_phi");

    int ibinPt_low = H3D_jetRjetPtjetPhi[iRun]->GetYaxis()->FindBin(PtCutLow);
    int ibinPt_high = H3D_jetRjetPtjetPhi[iRun]->GetYaxis()->FindBin(PtCutHigh);
    if (ibinPt_low == 0) 
      cout << "WARNING: Phi_RunComparison is counting the underflow with the chosen PtRange" << endl;
    if (ibinPt_high == H3D_jetRjetPtjetPhi[iRun]->GetYaxis()->GetNbins()+1) 
      cout << "WARNING: Phi_RunComparison is counting the overflow with the chosen PtRange" << endl;
    ibinJetRadius = H3D_jetRjetPtjetPhi[iRun]->GetXaxis()->FindBin(jetRadius);

    H1D_jetPhi[iRun] = (TH1D*)H3D_jetRjetPtjetPhi[iRun]->ProjectionZ("jetPhi_"+Runs[iRun]+Form("%.1f", PtCutLow)+"<pt<"+Form("%.1f", PtCutHigh), ibinJetRadius,ibinJetRadius, ibinPt_low, ibinPt_high);
    H1D_jetPhi_rebinned[iRun] = (TH1D*)H1D_jetPhi[iRun]->Rebin(1.,"H1D_jetPhi_rebinned_"+Runs[iRun]);

    // NormaliseYieldToNJets(H1D_jetPhi_rebinned[iRun]);
    NormaliseYieldToNEvents(H1D_jetPhi_rebinned[iRun], GetNEventsSel8(file_O2Analysis_list[iRun]));

    H1D_jetPhi_rebinned_ratios[iRun] = (TH1D*)H1D_jetPhi_rebinned[iRun]->Clone("H1D_jetPhi_rebinned_ratios"+Runs[iRun]);
    H1D_jetPhi_rebinned_ratios[iRun]->Reset("M");
    divideSuccess = H1D_jetPhi_rebinned_ratios[iRun]->Divide(H1D_jetPhi_rebinned[iRun], H1D_jetPhi_rebinned[0]);
  }

  TString* pdfName = new TString("jet_"+jetType[iJetType]+"_"+jetLevel[iJetLevel]+"_RunComp_R="+Form("%.1f", jetRadius)+"_Phi_@pT["+Form("%03.0f", PtCutLow)+","+Form("%03.0f", PtCutHigh)+"]");
  TString* pdfName_ratio = new TString("jet_"+jetType[iJetType]+"_"+jetLevel[iJetLevel]+"_RunComp_R="+Form("%.1f", jetRadius)+"_Phi_@pT["+Form("%03.0f", PtCutLow)+","+Form("%03.0f", PtCutHigh)+"]_ratio");

  TString textContext("#splitline{"+contextJetRadius(jetRadius)+"}{"+contextPtRange(PtRange)+"}");

  Draw_TH1_Histograms_in_one(H1D_jetPhi_rebinned, Runs, nRuns, textContext, pdfName, texPhiX, texNormPhiYield, texCollisionDataInfo, "");

  if (divideSuccess == true) {
    Draw_TH1_Histograms_in_one(H1D_jetPhi_rebinned_ratios, Runs, nRuns, textContext, pdfName_ratio, texPhiX, texRatioRuns, texCollisionDataInfo, "standardratio,avoidFirst");
  }
  else {
    cout << "Divide failed in Draw_Phi_RunComparison" << endl;
  }
}

void Draw_Pt_ratio_etaNeg_etaPos_RadiusComparison(int iRun, float* etaRange) {

  TH3D* H3D_jetRjetPtjetEta;
  TH1D* H1D_jetPt_left[nRadius];
  TH1D* H1D_jetPt_right[nRadius];
  TH1D* H1D_jetPt_left_rebinned[nRadius];
  TH1D* H1D_jetPt_right_rebinned[nRadius];
  TH1D* H1D_jetPt_rebinned_ratios[nRadius];

  H3D_jetRjetPtjetEta = (TH3D*)file_O2Analysis_list[iRun]->Get(analysisWorkflow+"/h3_jet_r_jet_pt_jet_eta");
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
    ibinJetRadius = H3D_jetRjetPtjetEta->GetXaxis()->FindBin(arrayRadius[iRadius]);
    H1D_jetPt_left[iRadius] = (TH1D*)H3D_jetRjetPtjetEta->ProjectionY("jetPt_left_"+RadiusLegend[iRadius], ibinJetRadius, ibinJetRadius, ibinEta_low, ibinEta_zero-1);
    H1D_jetPt_right[iRadius] = (TH1D*)H3D_jetRjetPtjetEta->ProjectionY("jetPt_right_"+RadiusLegend[iRadius], ibinJetRadius, ibinJetRadius, ibinEta_zero, ibinEta_high);
    H1D_jetPt_left_rebinned[iRadius] = (TH1D*)H1D_jetPt_left[iRadius]->Rebin(1.,"H1D_jetPt_left_rebinned_"+RadiusLegend[iRadius]);
    H1D_jetPt_right_rebinned[iRadius] = (TH1D*)H1D_jetPt_right[iRadius]->Rebin(1.,"H1D_jetPt_right_rebinned_"+RadiusLegend[iRadius]);

    // NormaliseYieldToNJets(H1D_jetPt_left_rebinned[iRadius]);
    // NormaliseYieldToNJets(H1D_jetPt_right_rebinned[iRadius]);
    NormaliseYieldToNEvents(H1D_jetPt_left_rebinned[iRadius], GetNEventsSel8(file_O2Analysis_list[iRun]));
    NormaliseYieldToNEvents(H1D_jetPt_right_rebinned[iRadius], GetNEventsSel8(file_O2Analysis_list[iRun]));

    H1D_jetPt_rebinned_ratios[iRadius] = (TH1D*)H1D_jetPt_left_rebinned[iRadius]->Clone("H1D_jetPt_rebinned_ratios"+Runs[iRadius]);
    H1D_jetPt_rebinned_ratios[iRadius]->Reset("M");
    divideSuccess = H1D_jetPt_rebinned_ratios[iRadius]->Divide(H1D_jetPt_right_rebinned[iRadius], H1D_jetPt_left_rebinned[iRadius]);
  }
 
  TString* pdfName = new TString("jet_"+jetType[iJetType]+"_"+jetLevel[iJetLevel]+"_"+Runs[iRun]+"_pT_etaRightLeftRatio");
  std::stringstream ss;
  ss << EtaCutLow << " < #eta < " << EtaCutHigh;
  TString textContext("#splitline{"+contextDataset(iRun)+"}{"+(TString)ss.str()+"}");

  if (divideSuccess == true) {
    Draw_TH1_Histograms_in_one(H1D_jetPt_rebinned_ratios, RadiusLegend, nRadius, textContext, pdfName, texPtX, texRatioEtaComparison, texCollisionDataInfo, "standardratio");
  }
  else {
    cout << "Divide failed in Draw_Pt_ratio_etaPos_etaNeg" << endl;
  }
}

void Draw_Pt_ratio_etaNeg_etaPos_RunComparison(float jetRadius, float* etaRange) {

  TH3D* H3D_jetRjetPtjetEta;
  TH1D* H1D_jetPt_left[nRuns];
  TH1D* H1D_jetPt_right[nRuns];
  TH1D* H1D_jetPt_left_rebinned[nRuns];
  TH1D* H1D_jetPt_right_rebinned[nRuns];
  TH1D* H1D_jetPt_rebinned_ratios[nRuns];

  float EtaCutLow = etaRange[0];
  float EtaCutHigh = etaRange[1];
  if (EtaCutLow > 0 || EtaCutHigh < 0 ) {
    cout << "eta=0 should be within the [EtaCutLow,EtaCutLow] range for the Draw_Pt_ratio_etaNeg_etaPos function" << endl;
    return;
  }
  int ibinJetRadius = 0;
  int nEvents;

  bool divideSuccess = false;

  for(int iRun = 0; iRun < nRuns; iRun++){
    H3D_jetRjetPtjetEta = (TH3D*)file_O2Analysis_list[iRun]->Get(analysisWorkflow+"/h3_jet_r_jet_pt_jet_eta");
    H3D_jetRjetPtjetEta->Sumw2();

    int ibinEta_low = H3D_jetRjetPtjetEta->GetZaxis()->FindBin(EtaCutLow);
    int ibinEta_zero = H3D_jetRjetPtjetEta->GetZaxis()->FindBin(0.);
    int ibinEta_high = H3D_jetRjetPtjetEta->GetZaxis()->FindBin(EtaCutHigh);
    if (ibinEta_low == 0) 
      cout << "WARNING: Pt_ratio_etaNeg_etaPos_RunComparison is counting the underflow with the chosen etaRange" << endl;
    if (ibinEta_high == H3D_jetRjetPtjetEta->GetZaxis()->GetNbins()+1) 
      cout << "WARNING: Pt_ratio_etaNeg_etaPos_RunComparison is counting the overflow with the chosen etaRange" << endl;

    ibinJetRadius = H3D_jetRjetPtjetEta->GetXaxis()->FindBin(jetRadius);

    H1D_jetPt_left[iRun] = (TH1D*)H3D_jetRjetPtjetEta->ProjectionY("jetPt_left_"+Runs[iRun],ibinJetRadius,ibinJetRadius, ibinEta_low, ibinEta_zero-1);
    H1D_jetPt_right[iRun] = (TH1D*)H3D_jetRjetPtjetEta->ProjectionY("jetPt_right_"+Runs[iRun],ibinJetRadius,ibinJetRadius, ibinEta_zero, ibinEta_high);
    H1D_jetPt_left_rebinned[iRun] = (TH1D*)H1D_jetPt_left[iRun]->Rebin(1.,"H1D_jetPt_left_rebinned_"+Runs[iRun]);
    H1D_jetPt_right_rebinned[iRun] = (TH1D*)H1D_jetPt_right[iRun]->Rebin(1.,"H1D_jetPt_right_rebinned_"+Runs[iRun]);
    
    // NormaliseYieldToNJets(H1D_jetPt_left_rebinned[iRun]);
    // NormaliseYieldToNJets(H1D_jetPt_right_rebinned[iRun]);
    NormaliseYieldToNEvents(H1D_jetPt_left_rebinned[iRun], GetNEventsSel8(file_O2Analysis_list[iRun]));
    NormaliseYieldToNEvents(H1D_jetPt_right_rebinned[iRun],  GetNEventsSel8(file_O2Analysis_list[iRun]));

    H1D_jetPt_rebinned_ratios[iRun] = (TH1D*)H1D_jetPt_left_rebinned[iRun]->Clone("H1D_jetPt_rebinned_ratios"+Runs[iRun]);
    H1D_jetPt_rebinned_ratios[iRun]->Reset("M");
    divideSuccess = H1D_jetPt_rebinned_ratios[iRun]->Divide(H1D_jetPt_right_rebinned[iRun], H1D_jetPt_left_rebinned[iRun]);
  }
 
  TString* pdfName = new TString("jet_"+jetType[iJetType]+"_"+jetLevel[iJetLevel]+"_RunComp_R="+Form("%.1f", jetRadius)+"_pT_etaRightLeftRatio");

  TString textContext(contextJetRadius(jetRadius));

  if (divideSuccess == true) {
    Draw_TH1_Histograms_in_one(H1D_jetPt_rebinned_ratios, Runs, nRuns, textContext, pdfName, texPtX, texRatioEtaComparison, texCollisionDataInfo, "standardratio");
  }
  else {
    cout << "Divide failed in Draw_Pt_ratio_etaPos_etaNeg" << endl;
  }
}

// To-do list:
// - implement carolina's macro of the gpad thing to automate the division of canvas based on how many plots one wants
// - change color and markers to make run comparison better
// - add info to 2D plots, like R=0.2 or R=0.4 or R=0.6 etc ; and clean them in general, the color scale is truncated and can't read it
// - changer nom run -> period
// - eta+ vs eta- comparison: make sure the separation bin is chosen so that there's no artificial imbalance




// option e for projection computes errors