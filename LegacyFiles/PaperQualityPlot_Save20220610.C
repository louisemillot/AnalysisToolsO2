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

void SetStyle(Bool_t graypalette=kFALSE);
void myPadSetUp(TPad *currentPad, float currentLeft=0.11, float currentTop=0.04, float currentRight=0.04, float currentBottom=0.15);
void DrawLogo (Int_t logo=0, Double_t xmin =  0.28, Double_t ymin= 0.68) ;
void FakeHistosOnlyForExample(TH1*&hstat, TH1*&hsyst, TH1*&hsystCorr);
void LoadLibs();

void PseudoEfficiency_HEPcomparison_allCountSignalRegion(TH1D* &hstat, TH1D* &hsyst, TH1D*&hsystCorr, Int_t ipart, TFile* file_O2Analysis, TString* SaveAs_Title, TString* &texXtitle, TString* &texYtitle, Double_t* pTbins, Int_t nbinpT);
void PseudoEfficiency_HEPcomparison_allCountSignalRegionOLD_withBugs(TH1D* &hstat, TH1D* &hsyst, TH1D*&hsystCorr, Int_t ipart, TFile* file_O2Analysis, TString* SaveAs_Title, Double_t* pTbins, Int_t nbinpT);
void RawSpectrum_O2data_allCountSignalRegion(TH1D* &hstat, TH1D* &hsyst, TH1D*&hsystCorr, Int_t ipart, TFile* file_O2Analysis, TString* SaveAs_Title, TString* &texXtitle, TString* &texYtitle, Double_t* pTbins, Int_t nbinpT);
void InvMass_Plot(TH1D* &hstat, TH1D* &hsyst, TH1D*&hsystCorr, Int_t ipart, TFile* file_O2Analysis, TString* SaveAs_Title, TString* &texXtitle, TString* &texYtitle, bool isMC);
void PtDistribution_O2data(TH1D* &hstat, TH1D* &hsyst, TH1D*&hsystCorr, Int_t ipart, TFile* file_O2Analysis, TString* SaveAs_Title, TString* &texXtitle, TString* &texYtitle);
void DcaHistoProcessing_O2data(TH1D* &hstat, TH1D* &hsyst, TH1D*&hsystCorr, Int_t ipart, TFile* file_O2Analysis);
void Efficiency_O2data_allCountSignalRegion(TH1D* &hstat, TH1D* &hsyst, TH1D*&hsystCorr, Int_t ipart, TFile* file_O2Analysis, TString* SaveAs_Title, TString* &texXtitle, TString* &texYtitle, Double_t* pTbins, Int_t nbinpT);
void InvMassDistributions_withFit(TH1D* &hstat, TH1D* &hsyst, TH1D* &hsystCorr, Int_t ipart, TFile* file_O2Analysis, TString* SaveAs_Title, TString* &texXtitle, TString* &texYtitle, Double_t* pTbins, Int_t nbinpT);
void InvMassDistributions_withFit_MC_datalike(TH1D* &hstat, TH1D* &hsyst, TH1D* &hsystCorr, Int_t ipart, TFile* file_O2Analysis, TString* SaveAs_Title, TString* &texXtitle, TString* &texYtitle, Double_t* pTbins, Int_t nbinpT);


void PtDifferential_SigmaOfFit(TH1D* &hstat, TH1D* &hsyst, TH1D*&hsystCorr, Int_t ipart, TFile* file_O2Analysis, TString* SaveAs_Title, TString* &texXtitle, TString* &texYtitle, Double_t* pTbins, Int_t nbinpT, bool isMC);
void RawSpectrum_O2data_withFit(TH1D* &hstat, TH1D* &hsyst, TH1D*&hsystCorr, Int_t ipart, TFile* file_O2Analysis, TString* SaveAs_Title, TString* &texXtitle, TString* &texYtitle, Double_t* pTbins, Int_t nbinpT);
void PseudoEfficiency_HEPcomparison_withFit(TH1D* &hstat, TH1D* &hsyst, TH1D*&hsystCorr, Int_t ipart, TFile* file_O2Analysis, TString* SaveAs_Title, TString* &texXtitle, TString* &texYtitle, Double_t* pTbins, Int_t nbinpT);
void Efficiency_O2MCdata_withFit(TH1D* &hstat, TH1D* &hsyst, TH1D*&hsystCorr, Int_t ipart, TFile* file_O2Analysis, TString* SaveAs_Title, TString* &texXtitle, TString* &texYtitle, Double_t* pTbins, Int_t nbinpT);
void Efficiency_O2data_truePt_trueV0s(TH1D* &hstat, TH1D* &hsyst, TH1D*&hsystCorr, Int_t ipart, TFile* file_O2Analysis, TString* SaveAs_Title, TString* &texXtitle, TString* &texYtitle, Double_t* pTbins, Int_t nbinpT);
void RawSpectrum_O2data_truePt_trueV0s(TH1D* &hstat, TH1D* &hsyst, TH1D*&hsystCorr, Int_t ipart, TFile* file_O2Analysis, TString* SaveAs_Title, TString* &texXtitle, TString* &texYtitle, Double_t* pTbins, Int_t nbinpT);

void V0Daughters_TrackingEfficiencies(TH1D* &hstat, TH1D* &hsyst, TH1D*&hsystCorr, Int_t ipart, TFile* file_O2Analysis, TString* SaveAs_Title, Double_t* pTbins, Int_t nbinpT);
void V0Daughters_TrackingEfficiencies_HistogramProcessing(TH1D* &daughter1, TH1D* &daughter2, Int_t ipart, TFile* file_O2Analysis, Double_t* pTbinsProtons, Int_t nbinpTProtons, Double_t* pTbinsPions, Int_t nbinpTPions);
void V0Daughters_TrackingEfficiencies_Daughter1(TH1D* &hstat, TH1D* &hsyst, TH1D*&hsystCorr, Int_t ipart, TFile* file_O2Analysis, TString* SaveAs_Title, TString* &texXtitle, TString* &texYtitle, TString* &Extra, Double_t* pTbinsProtons, Int_t nbinpTProtons, Double_t* pTbinsPions, Int_t nbinpTPions);
void V0Daughters_TrackingEfficiencies_Daughter2(TH1D* &hstat, TH1D* &hsyst, TH1D*&hsystCorr, Int_t ipart, TFile* file_O2Analysis, TString* SaveAs_Title, TString* &texXtitle, TString* &texYtitle, TString* &Extra, Double_t* pTbinsProtons, Int_t nbinpTProtons, Double_t* pTbinsPions, Int_t nbinpTPions);

void pT_Spectrum_postAnalyserCuts(TH1D* &hstat, TH1D* &hsyst, TH1D*&hsystCorr, Int_t ipart, TFile* file_O2Analysis, TString* SaveAs_Title, TString* &texXtitle, TString* &texYtitle, Double_t* pTbins, Int_t nbinpT, bool isMC);
void pT_Spectrum_preAnalyserCutsK0S(TH1D* &hstat, TH1D* &hsyst, TH1D*&hsystCorr, Int_t ipart, TFile* file_O2Analysis, TString* SaveAs_Title, TString* &texXtitle, TString* &texYtitle, Double_t* pTbins, Int_t nbinpT, bool isMC);

void Efficiency_DcaScan(TH1D* &hstat, TH1D* &hsyst, TH1D*&hsystCorr, Int_t ipart, TFile* file_O2Analysis, TString* SaveAs_Title, TString* &texXtitle, TString* &texYtitle, Double_t* pTbins, Int_t nbinpT, Double_t maxDCAcut);


bool reject = false ;
float fparab(float *x, float *par);
float fline(float *x, float *par);
// Histogram processing functions; 

// returns pt diff raw spectrum in hPtDiff_RawSpectrum_O2
// void PtDifferential_Fit(TH1D* &hPtDiff_RawSpectrum_O2, TH3D* H3D_DetectedV0s_O2, TH3D* H3D_DetectedV0s_O2_rebinnedZ, Int_t ipart, Double_t* pTbins, Int_t nbinpT, Double_t* parGaussianParab_Sigma);
void PtDifferential_Fit(TH1D* &hPtDiff_RawSpectrum_O2, TH3D* H3D_DetectedV0s_O2, TH3D* H3D_DetectedV0s_O2_rebinnedZ, Int_t ipart, Double_t* pTbins, Int_t nbinpT, Double_t* parGaussianParab_Sigma, Double_t* parGaussianParab_Sigma_error);

//Custom Get functions
void Get_McEfficiency_vsPt(TH1D* &hMcEfficiency_vsPt, TH3D* H3D_detectedV0s_TruePt, TH1D* TrueV0PtSpectrum, Int_t ipart, Double_t* pTbins, Int_t nbinpT, Int_t ibinXaxisCut);
void Get_McSignalBackgroundRatio_vsPt(TH1D* &hSignalBackgroundRatio_vsPt, TH3D* H3D_detectedV0s_DataLike, Int_t ipart, Double_t* pTbins, Int_t nbinpT, Int_t ibinXaxisCut);

// Preferred colors and markers
const Int_t fillColors[] = {kGray+1,  kRed-10, kBlue-9, kGreen-8, kMagenta-9, kOrange-9,kCyan-8,kYellow-7}; // for syst bands
const Int_t colors[]     = {kBlack, kRed+1 , kBlue+1, kGreen+3, kMagenta+1, kOrange-1,kCyan+2,kYellow+2};
const Int_t markers[]    = {kFullCircle, kFullSquare,kOpenCircle,kOpenSquare,kOpenDiamond,kOpenCross,kFullCross,kFullDiamond,kFullStar,kOpenStar};

//Options:
TFile* file_O2Analysis = new TFile("AnalysisResults.root");
// TFile* file_O2Analysis = new TFile("AnalysisResults_noMC.root");
// TFile* file_O2Analysis = new TFile("AnalysisResultsMC.root");

const Int_t numPart = 3;
const Int_t ipart = 0;
const TString NamehistoInvMass[numPart] = {"InvMassK0S", "InvMassLambda", "InvMassAntiLambda"};//,"InvMassXiPlus", "InvMassXiMinus", "InvMassOmegaPlus", "InvMassOmegaMinus"};
const TString NamePart[numPart] = {"K0Short", "Lambda", "AntiLambda"};//, "XiPlu", "XiMin", "OmPlu", "OmMin"};
const TString NamePart_Latex[numPart] = {"K^{0}_{S}", "#Lambda", "#bar{#Lambda}"};//, "XiPlu", "XiMin", "OmPlu", "OmMin"};
const Float_t LowerLimitDisplay_Mass[numPart] = {0.455, 1.085, 1.085};
const Float_t UpperLimitDisplay_Mass[numPart] = {0.525, 1.145, 1.145};

const Float_t MassPart[numPart] = {0.497611, 1.115683, 1.115683};

const Float_t min_range_signal[numPart] = {0.47, 1.105, 1.105};
const Float_t max_range_signal[numPart] = {0.51, 1.1205, 1.1205};
const Float_t liminf[numPart] = {0.45, 1.08, 1.08};
const Float_t limsup[numPart] = {0.55, 1.18, 1.18};

// pT binning for K0S, Lambda and Antilambda
// Double_t pTbins[numPart][20] = {{0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.6, 2.0, 2.4, 3.0},{0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 2.0, 2.4, 3.0, 3.5},{0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 2.0, 2.4, 3.0, 3.5}};
// Int_t nbinpT[numPart] = {16,9,9};
// Double_t pTbins[numPart][20] = {{0.0, 0.2, 0.4, 0.6, 0.8, 1., 1.2, 1.6, 2., 2.4, 3.0},{0.0, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 2.0, 2.4, 3.5},{0.0, 0.6, 0.8, 1.2, 1.6, 2.0, 2.4, 3.5}};
// Int_t nbinpT[numPart] = {10,9,7};
// Double_t pTbins[numPart][20] = {{0.0, 0.4, 0.8, 1.2, 1.6, 3.0},{0.0, 0.6, 1.0, 1.4, 2.0, 3.5},{0.0, 0.6, 1.0, 1.4, 2.0, 3.5}};
// Int_t nbinpT[numPart] = {5,5,5};
// Double_t pTbins[numPart][20] = {{0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.6, 2.0, 2.4, 3.0},{0.0, 0.6, 1.0, 1.4, 2.0, 3.0},{0.0, 0.6, 1.0, 1.4, 2.0, 3.0}};
// Int_t nbinpT[numPart] = {16,5,5};
Double_t pTbins[numPart][4] = {{0., 10.},{0., 3.},{0., 3.}};
Int_t nbinpT[numPart] = {1,1,1};

//For V0 daughter tracking efficiencies
Int_t daughterSpecies = 0; //0 pion, 1 proton
Double_t pTbinsPions[20] = {0.0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.4, 0.5, 0.6, 0.8, 1.0, 1.4, 1.8, 2.2, 2.6};
Int_t nbinpTPions = 15;
Double_t pTbinsProtons[20] = {0.0, 0.1, 0.2, 0.3, 0.5, 0.7, 0.9, 1.2, 1.4, 1.6, 2.0, 3.0, 4.0};
Int_t nbinpTProtons = 12;
Double_t daughterWindow[2] = {pTbinsPions[nbinpTPions], pTbinsProtons[nbinpTProtons]};

// const Float_t MassPart[numPart] = {0.497611, 1.115683, 1.115683};// 1.32171, 1.32171, 1.67245, 1.67245};
// const Float_t WidthPartLimitsFit[2][numPart] = {{0.003, 0.001, 0.001}, {0.01, 0.005, 0.005}};//{{0.003, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001}, {0.01, 0.005, 0.005, 0.005, 0.005, 0.008, 0.008}};

// === Commonly used x/ titles: ===
// pt invariant yields
TString* texPtX = new TString("#it{p}_{T} (GeV/#it{c})");
TString* texPtY = new TString("1/#it{N}_{ev} 1/(2#pi#it{p}_{T}) d#it{N}/(d#it{p}_{T}d#it{y}) ((GeV/#it{c})^{-2})");
// mt invariant yields
TString* texMtX = new TString("#it{m}_{T} (GeV/#it{c}^{2})");
TString* texMtY = new TString("1/#it{N}_{ev} 1/(2#pi#it{m}_{T}) d#it{N}/(d#it{m}_{T}d#it{y}) ((GeV/#it{c}^{2})^{-2})"); 
// Invariant mass with decay products K and pi
TString* texMassX = new TString("#it{M}_{K#pi} (GeV/#it{c}^{2})");
TString* texMassY = new TString("d#it{N}/(d#it{M}_{K#pi})");
// Invariant mass with decay products pi+ and pi-
TString* texMassPiPiX = new TString("#it{M}_{#it{#pi}^{+}#it{#pi}^{-}} (GeV/#it{c}^{2})");
TString* texMassPiPiY = new TString("d#it{N}/(d#it{M}_{#it{#pi}^{+}#it{#pi}^{-}})");
// <pt>, npart
TString* texMeanPt = new TString("#LT#it{p}_{T}#GT (GeV/#it{c})");
TString* texMeanNpart = new TString("#LT#it{N}_{part}#GT");
//AIMERIC TILTES
TString* texptDifferentialYield = new TString("1/#it{N}_{ev} d#it{N}/d#it{p}_{T} (GeV/#it{c})^{-1}");
TString* texptDifferentialYield_HEPratio_pseudoEff = new TString("O2/HEP ratio of 1/#it{N}_{ev} d#it{N}/d#it{p}_{T} (GeV/#it{c})^{-1}");
// TString* texRawYield = new TString("raw yield 1/#it{N}_{ev} d#it{N}/d#it{p}_{T}");
TString* texRawYield = new TString("d#it{N}/d#it{p}_{T}");
TString* texCount = new TString("count");
TString* texSgimaGaussFit = new TString("#sigma of Gaussian fit (GeV/#it{c}^{2})");
TString* texEfficiency = new TString("Efficiency #it{V0}_{detected}/#it{V0}_{MC} (%)");
TString* texDaughterRecoEfficiency = new TString("Track Reco Efficiency #it{N}_{reco}/#it{N}_{genMC} (%)");

TString* texK0S_daughter_PiPlus = new TString("#pi^{+}");
TString* texLambda_daughter_PiMinus = new TString("#pi^{-}");
TString* texAntiLambda_daughter_PiPlus = new TString("#pi^{+}");

TString* texK0S_daughter_PiMinus = new TString("#pi^{-}");
TString* texLambda_daughter_protonPlus = new TString("p");
TString* texAntiLambda_daughter_protonMinus = new TString("#bar{p}");

const TString* texInvMass_K0ShortDecay = new TString("M_{#pi^{+}#pi^{-}} (GeV/#it{c}^{2})");
const TString* texInvMass_LambdaDecay = new TString("M_{p#pi^{-}} (GeV/#it{c}^{2})");
const TString* texInvMass_AntiLambdaDecay = new TString("M_{#bar{p}#pi^{+}} (GeV/#it{c}^{2})");

const TString* texInvMassDecays_titles[numPart] = {texInvMass_K0ShortDecay,texInvMass_LambdaDecay,texInvMass_AntiLambdaDecay};

void PaperQualityPlot_Original() {
  // Load necessary libraries
  LoadLibs();
  // Set the default style
  SetStyle();

  //histograms definition and processing
  TH1D* hstat, *hsyst, *hsystCorr;
  Int_t icolor=0;

  TString* SaveAs_Title = new TString("");
  TString* texXtitle = new TString("");
  TString* texYtitle = new TString("");
  TString* Extra = new TString("");

  // PtDistribution_O2data(hstat, hsyst, hsystCorr, ipart, file_O2Analysis, SaveAs_Title);
  // PseudoEfficiency_HEPcomparison_withFit(hstat, hsyst, hsystCorr, ipart, file_O2Analysis, SaveAs_Title, texXtitle, texYtitle, pTbins[ipart], nbinpT[ipart]);
  // InvMassDistributions_withFit(hstat, hsyst, hsystCorr, ipart, file_O2Analysis, SaveAs_Title, texXtitle, texYtitle, pTbins[ipart], nbinpT[ipart]);
  InvMassDistributions_withFit_MC_datalike(hstat, hsyst, hsystCorr, ipart, file_O2Analysis, SaveAs_Title, texXtitle, texYtitle, pTbins[ipart], nbinpT[ipart]);

  // PseudoEfficiency_HEPcomparison_allCountSignalRegion(hstat, hsyst, hsystCorr, ipart, file_O2Analysis, SaveAs_Title, texXtitle, texYtitle, pTbins[ipart], nbinpT[ipart]);
  // PseudoEfficiency_HEPcomparison_allCountSignalRegionOLD_withBugs(hstat, hsyst, hsystCorr, ipart, file_O2Analysis, SaveAs_Title, pTbins[ipart], nbinpT[ipart]);
  // RawSpectrum_O2data_withFit(hstat, hsyst, hsystCorr, ipart, file_O2Analysis, SaveAs_Title, texXtitle, texYtitle, pTbins[ipart], nbinpT[ipart]);
  // Efficiency_O2MCdata_withFit(hstat, hsyst, hsystCorr, ipart, file_O2Analysis, SaveAs_Title, texXtitle, texYtitle, pTbins[ipart], nbinpT[ipart]);
  // Efficiency_O2data_allCountSignalRegion(hstat, hsyst, hsystCorr, ipart, file_O2Analysis, SaveAs_Title, texXtitle, texYtitle, pTbins[ipart], nbinpT[ipart]);
  // Efficiency_O2data_truePt_trueV0s(hstat, hsyst, hsystCorr, ipart, file_O2Analysis, SaveAs_Title, texXtitle, texYtitle, pTbins[ipart], nbinpT[ipart]);
  // RawSpectrum_O2data_truePt_trueV0s(hstat, hsyst, hsystCorr, ipart, file_O2Analysis, SaveAs_Title, texXtitle, texYtitle, pTbins[ipart], nbinpT[ipart]);
  // V0Daughters_TrackingEfficiencies_Daughter2(hstat, hsyst, hsystCorr, ipart, file_O2Analysis, SaveAs_Title, texXtitle, texYtitle, Extra, pTbinsProtons, nbinpTProtons, pTbinsPions, nbinpTPions);

  Double_t maxDCAcut = 1;
  TH3D* H3D_detectedV0s_TruePt_DcaHist = (TH3D*)file_O2Analysis->Get("lambdakzero-analysis-mc/h3dMass"+NamePart[ipart]+"Dca");
  Double_t* DCAbins_fullFileRange = (Double_t*)H3D_detectedV0s_TruePt_DcaHist->GetXaxis();
  // Efficiency_DcaScan(hstat, hsyst, hsystCorr, ipart, file_O2Analysis, SaveAs_Title, texXtitle, texYtitle, pTbins[ipart], nbinpT[ipart], maxDCAcut);

  bool isMC = true;
  // PtDifferential_SigmaOfFit(hstat, hsyst, hsystCorr, ipart, file_O2Analysis, SaveAs_Title, texXtitle, texYtitle, pTbins[ipart] ,nbinpT[ipart], isMC);
  // InvMass_Plot(hstat, hsyst, hsystCorr, ipart, file_O2Analysis, SaveAs_Title, texXtitle, texYtitle, isMC);
  // pT_Spectrum_postAnalyserCuts(hstat, hsyst, hsystCorr, ipart, file_O2Analysis, SaveAs_Title, texXtitle, texYtitle, pTbins[ipart], nbinpT[ipart], isMC);
  // pT_Spectrum_preAnalyserCutsK0S(hstat, hsyst, hsystCorr, ipart, file_O2Analysis, SaveAs_Title, texXtitle, texYtitle, pTbins[ipart], nbinpT[ipart], isMC);

 // Prepare Figure, please stick to the default canvas size(s) unless absolutely necessary in your case
  // Rectangular
  TCanvas *cfig = new TCanvas("cfig", "Alice Figure Template", 800, 600); 
  // Square
  //TCanvas *cfig = new TCanvas("cfig", "Alice Figure Template", 800, 800); 
  // cfig->SetLogy();

  cfig->cd();
  // Set Titles etc..
  // DrawFrame(xmin, ymin, xmax, ymax)
  // TH1 * h = cfig->DrawFrame(0,0.0001,5,0.02);// RawSpectrum_O2data_allCountSignalRegion LOGY scale
  // TH1 * h = cfig->DrawFrame(0,0,3,0.1); // RawSpectrum_O2data_K0S per event
  // TH1 * h = cfig->DrawFrame(0,0,3,4000); // RawSpectrum_O2data_K0S
  // TH1 * h = cfig->DrawFrame(0,0,3.5,0.01); // RawSpectrum_O2data_Lambda per event
  // TH1 * h = cfig->DrawFrame(0,0,3.5,500); // RawSpectrum_O2data_Lambda
  // TH1 * h = cfig->DrawFrame(0.45,0,0.55,1000); // PseudoEfficiency_HEPcomparison_withFit
  TH1 * h = cfig->DrawFrame(0,0,3,0.017); // PtDifferential_SigmaOfFit
  // TH1 * h = cfig->DrawFrame(0,0,3,0.15); // PseudoEfficiency_HEPcomparison_allCountSignalRegion
  // TH1 * h = cfig->DrawFrame(0,0,3,0.06); // Efficiency_K0S
  // TH1 * h = cfig->DrawFrame(0,0,3,0.7); // Efficiency_K0S new MC
  // TH1 * h = cfig->DrawFrame(0,0,3,0.64); // Efficiency_Lambda
  // TH1 * h = cfig->DrawFrame(0,0,10,0.05); // Momentum Distribution
  // TH1 * h = cfig->DrawFrame(0.45,0,0.55,200); // InvMass K0S
  // TH1 * h = cfig->DrawFrame(LowerLimitDisplay_Mass[ipart],0,UpperLimitDisplay_Mass[ipart],1.1*hstat->GetMaximum()); // InvMass Lambda
  // TH1 * h = cfig->DrawFrame(0,0,0.1,500); // DcaHistoProcessing_O2data
  // TH1 * h = cfig->DrawFrame(0,0,daughterWindow[daughterSpecies],0.2); // DaughterTrackingEfficiencies

  // if (0)
  //   {
  //     TLegend *legt = new TLegend( 0.46, 0.6, 0.72, 0.8, "Test Labels Legend");
  //     legt->AddEntry((TObject*)0, texMtX, "");
  //     legt->AddEntry((TObject*)0, texMtY, "");
  //     legt->AddEntry((TObject*)0, texMassX, "");
  //     legt->AddEntry((TObject*)0, texMeanPt, "");
  //     legt->AddEntry((TObject*)0, texMeanNpart, "");
  //     legt->SetFillColor(0);
  //     legt->SetTextSize(gStyle->GetTextSize()*0.6);
  //     legt->Draw();
  //   }

  // Set titles
  // h->SetXTitle(texPtX);
  cout << "title: " << texXtitle->Data() <<endl;
  h->SetXTitle(texXtitle->Data());
  // Please be consistent on the y label
  // h->SetYTitle(texptDifferentialYield);//(texPtY);
  h->SetYTitle(texYtitle->Data());//(texPtY);

//   h->GetXaxis()->SetTitle(texXtitle);
//   h->GetYaxis()->SetTitle(texYtitle);

  // Draw your histos here:
  hsystCorr->SetFillColor(fillColors[icolor]);
  hsyst    ->SetFillColor(fillColors[icolor]);
  // hsystCorr->Draw("E3,same"); //SHOULD BE PLOTTED EVENTUALLY
  hsyst->SetFillStyle(0); // To draw empty boxes
  hsyst->SetLineColor(colors[icolor]); // To draw empty boxes
  hsystCorr->SetLineColor(colors[icolor]); // To draw empty boxes
  // hsyst->Draw("E2,same"); //SHOULD BE PLOTTED EVENTUALLY
  hstat->Draw("E,same");
  hstat->SetMarkerStyle(markers[0]);
  // use the same color for markers and lines
  hstat->SetMarkerColor(colors [icolor]);
  hstat->SetLineColor  (colors [icolor]);

  // Draw the logo   
  //  0: Just "ALICE" (for final data), to be added only if ALICE does not appear otherwise (e.g. in the legend)
  //  >0: ALICE Preliminary
  // DrawLogo(1, 0.59, 0.81);

  // You should always specify the colliding system
  // NOTATION: pp, p-Pb, Pb-Pb. 
  // Don't forget to use #sqrt{s_{NN}} for p-Pb and Pb-Pb
  // You can change the position of this with
  TLatex * textContext = new TLatex (0.18,0.82,"#bf{ALICE Run 3 Performance}"); //BOLD
  textContext->SetTextSize(0.05);
  textContext->SetNDC(kTRUE); //remove if I want x,y in TLatex to be in the coordinate system of the histogram
  textContext->Draw();
  TLatex * textColl = new TLatex (0.18,0.75,"pp #sqrt{#it{s}} = 900 GeV, #it{B_{z}} = 0.2 T");
  textColl->SetTextSize(0.04);
  textColl->SetNDC(kTRUE); //remove if I want x,y in TLatex to be in the coordinate system of the histogram
  textColl->Draw();
  TLatex * text_part = new TLatex (0.18,0.70,NamePart_Latex[ipart]);
  text_part->SetTextSize(0.04);
  text_part->SetNDC(kTRUE); //remove if I want x,y in TLatex to be in the coordinate system of the histogram
  text_part->Draw();
  TLatex * text_extra = new TLatex (0.65,0.55,*Extra);
  text_extra->SetNDC(kTRUE); //remove if I want x,y in TLatex to be in the coordinate system of the histogram
  text_extra->Draw();
  // TLatex * text2 = new TLatex (0.55,55,"V0A Multiplicity Classes (Pb-Side)");
  // text2->SetTextSizePixels(24);
  // text2->Draw();

  // //Legend, if needed
  // TLegend * leg = new TLegend(  0.19,  0.19,  0.57, 0.42);
  // leg->AddEntry(hstat,     "0-5\%, stat errors",   "LPE");
  // leg->AddEntry(hsyst,     "syst error (Uncorrelated)",  "F");
  // leg->AddEntry(hsystCorr, "syst error (Correlated)",    "F" );
  // leg->SetFillColor(0);
  // leg->SetTextSize(gStyle->GetTextSize()*0.8);
  // leg->Draw();

  // // Save to HEP data

  // AliHEPDataParser * hepParser = new AliHEPDataParser(hstat, hsyst);
  // hepParser->SetTitle("pt distribution of pi+-, arXiv:XXXX.YYYY");
  // hepParser->SetName("1/Nev 1/p_T 1/2pi d^2N/(dp_Tdy) (GeV/c)^{-1}"); 
  // hepParser->SetXaxisName("PT IN GEV/c");
  // hepParser->SetReaction("RE: P PB --> PI + X");
  // hepParser->SetEnergy("SQRT(SNN) : 5020.0 GeV");
  // hepParser->SetRapidityRange("YRAP : -0.5 - +0.5");
  // hepParser->SaveHEPDataFile("figTemplateHEPData.txt");    // it must be specified explicity if graphs are to be used



  // cfig->SaveAs("HEPcomp_ratio_RawSpectrum_PtDifferential_K0sCount_WithSel8.pdf","pdf"); //HEPcomp_Ratio_
  // cfig->SaveAs("RawSpectrumLog_PtDifferential_K0sCount_WithSel8.pdf","pdf"); //HEPcomp_Ratio_
  // cfig->SaveAs("DCA_K0S_WithSel8.pdf","pdf"); 
  // cfig->SaveAs(strcat(SaveAs_Title, ".pdf"),"pdf"); 
  cfig->SaveAs(*SaveAs_Title+"_"+NamePart[ipart]+".pdf","pdf"); 

}

//________________________________
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

void DrawLogo (Int_t logo, Double_t xmin, Double_t ymin) {

  // Logo is not needed anymore, now we only write alice preliminary
  // Logo:
  // 0: Justr writes "ALICE" (for final data)
  // Anything eles: writes "ALICE Preliminary"

  TLatex *   tex = new TLatex(xmin,ymin, logo ? "ALICE Preliminary" : "ALICE");
  tex->SetNDC();
  tex->SetTextFont(42);
  tex->Draw();

  // OLD logo
  //  TPad * currentPad = gPad;
  // Double_t AliLogo_LowX =xmin;
  // Double_t AliLogo_LowY = ymin;
  // Double_t AliLogo_Height = size;
  // //ALICE logo is a  file that is 821x798 pixels->should be wider than a square
  // Double_t AliLogo_Width  = (821./798.) * AliLogo_Height * gPad->GetWh() / gPad->GetWw();
  
  // TPad *myPadLogo = new TPad("myPadLogo", "Pad for ALICE Logo",AliLogo_LowX,AliLogo_LowY,AliLogo_LowX+AliLogo_Width,AliLogo_LowY+AliLogo_Height);
  // myPadSetUp(myPadLogo,0,0,0,0);
  // //myPadLogo->SetFixedAspectRatio(1);
  // myPadLogo->Draw();
  // myPadLogo->cd();
  // if (logo == 0) {
  //   myPadLogo->SetFillColor(2); // color to first figure out where is the pad then comment !
  // } else if (logo == 1){
  //   TASImage *myAliceLogo = new TASImage(performanceLogoPath);
  //   myAliceLogo->Draw();
  // } else if (logo == 2) {
  //   TASImage *myAliceLogo = new TASImage(preliminaryLogoPath);
  //   myAliceLogo->Draw();
  // }
  // // go back to the old pad
  // currentPad->cd();

}

void myPadSetUp(TPad *currentPad, float currentLeft, float currentTop, float currentRight, float currentBottom){
  currentPad->SetLeftMargin(currentLeft);
  currentPad->SetTopMargin(currentTop);
  currentPad->SetRightMargin(currentRight);
  currentPad->SetBottomMargin(currentBottom);
  return;
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
  gStyle->SetLabelSize(0.045,"xyz");
  gStyle->SetLabelOffset(0.01,"y");
  gStyle->SetLabelOffset(0.01,"x");
  gStyle->SetLabelColor(kBlack,"xyz");
  gStyle->SetTitleSize(0.05,"xyz");
  gStyle->SetTitleOffset(1.25,"y");
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

void FakeHistosOnlyForExample(TH1* &hstat, TH1* &hsyst, TH1*&hsystCorr) {

  TF1 * fExpo = new TF1 ("fExpo", "expo");
  fExpo->SetParameters(10, -0.3);
  hstat     = new TH1F("hstat", "hstat", 100, 0, 10);
  hsyst     = new TH1F("hsyst", "hsyst", 100, 0, 10);
  hsystCorr = new TH1F("hsystCorr", "hsystCorr", 100, 0, 10);
  hstat     ->Sumw2(); 
  hsyst     ->Sumw2();
  hsystCorr ->Sumw2(); 
  hstat->FillRandom("fExpo",20000);
  Int_t nbinx = hstat->GetNbinsX();
  
  for(Int_t ibinx = 1; ibinx <= nbinx; ibinx++){
    hsyst->SetBinContent (ibinx, hstat->GetBinContent(ibinx));
    hsyst->SetBinError   (ibinx, hstat->GetBinContent(ibinx)*0.08);
    hsystCorr->SetBinContent (ibinx, hstat->GetBinContent(ibinx));
    hsystCorr->SetBinError   (ibinx, hstat->GetBinContent(ibinx)*0.15);
  }
  

}

void RawSpectrum_O2data_allCountSignalRegion(TH1D* &hstat, TH1D* &hsyst, TH1D*&hsystCorr, Int_t ipart, TFile* file_O2Analysis, TString* SaveAs_Title, TString* &texXtitle, TString* &texYtitle, Double_t* pTbins, Int_t nbinpT) {

  //O2 file
  // TFile* file_O2Analysis = new TFile("AnalysisResults_WithSel8Cut.root");

  TH1I* H1I_SelectedEventCount = (TH1I*)file_O2Analysis->Get("lambdakzero-analysis/hSelectedEventCounter");
  float SelectedEventCount = H1I_SelectedEventCount->GetEntries();

  TH3D* H3D_DetectedV0s_O2 = (TH3D*)file_O2Analysis->Get("lambdakzero-analysis/h3dMass"+NamePart[ipart]);
  TH1D* H1D_DetectedV0s_PtAxis_O2 = (TH1D*)H3D_DetectedV0s_O2->ProjectionY("InvMass_K0s",0,-1,0,-1);
  // TH1D* hPtDiff_RawSpectrum_O2 = (TH1D*)H1D_DetectedV0s_PtAxis_O2->Clone("hPtDiff_RawSpectrum");
  // hPtDiff_RawSpectrum_O2->Reset("M"); //should be empty, only used that to get same Pt binning as used in analysis
  
  //Mass window in which we count K0S
  float MassSignalInterval_lowEdge = 0.48;
  float MassSignalInterval_upEdge = 0.5;
  int MassSignalInterval_lowEdgeBin = H3D_DetectedV0s_O2->GetZaxis()->FindBin(MassSignalInterval_lowEdge);
  int MassSignalInterval_upEdgeBin = H3D_DetectedV0s_O2->GetZaxis()->FindBin(MassSignalInterval_upEdge);

  // Double_t xbins[9] = {0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.6, 2.0, 3.0};
  // Int_t nbinpT = 8;
  TH1D *H1D_DetectedV0s_O2_SmallpTinterval[nbinpT];
  TH1D *hPtDiff_RawSpectrum_O2 = new TH1D("hPtDiff_RawSpectrum_O2","hPtDiff_RawSpectrum_O2",nbinpT,pTbins);

  for(int ibinPt = 1; ibinPt <= nbinpT; ibinPt++){
    int iPt = ibinPt-1; //different convention for bins that start at 1 and arrays that start at 0

    int ibinPt_originalAxis_low = H3D_DetectedV0s_O2->GetYaxis()->FindBin(pTbins[iPt]);
    int ibinPt_originalAxis_up = H3D_DetectedV0s_O2->GetYaxis()->FindBin(pTbins[iPt+1]) - 1;
    cout << "ibinPt: " << ibinPt << "; originalAxis [low;up] = [" << ibinPt_originalAxis_low << ";" << ibinPt_originalAxis_up << "]" << endl;

    H1D_DetectedV0s_O2_SmallpTinterval[iPt] = (TH1D*)H3D_DetectedV0s_O2->ProjectionZ("InvMass_"+NamePart[ipart]+Form("_Pt%i", ibinPt),0,-1,ibinPt_originalAxis_low,ibinPt_originalAxis_up);

    float SignalCountRough = H1D_DetectedV0s_O2_SmallpTinterval[iPt]->Integral(MassSignalInterval_lowEdgeBin,MassSignalInterval_upEdgeBin);

    float dpT = hPtDiff_RawSpectrum_O2->GetXaxis()->GetBinWidth(ibinPt);
    float dN_dpT = SignalCountRough *1./dpT;
    hPtDiff_RawSpectrum_O2->SetBinContent(ibinPt,dN_dpT);
    hPtDiff_RawSpectrum_O2->SetBinError(ibinPt,sqrt(SignalCountRough*1./dpT));
  }

  // H3D_DetectedV0s_O2->RebinY(10);
  // hPtDiff_RawSpectrum_O2->Rebin(10);

  // Int_t nbinPt = hPtDiff_RawSpectrum_O2->GetNbinsX();
  // for(Int_t ibinPt = 1; ibinPt <= nbinPt; ibinPt++){
  //   float SignalCountRough = H3D_DetectedV0s_O2->Integral(0,-1,ibinPt,ibinPt+1,MassSignalInterval_lowEdgeBin,MassSignalInterval_upEdgeBin);

  //   float dpT = hPtDiff_RawSpectrum_O2->GetXaxis()->GetBinWidth(ibinPt);
  //   float dN_dpT = SignalCountRough *1./dpT;
  //   hPtDiff_RawSpectrum_O2->SetBinContent(ibinPt,dN_dpT);
  //   hPtDiff_RawSpectrum_O2->SetBinError(ibinPt,sqrt(SignalCountRough*1./dpT));
  // }

  hPtDiff_RawSpectrum_O2->Scale(1./SelectedEventCount);
  hstat = (TH1D*)hPtDiff_RawSpectrum_O2->Clone("hstat");
  // hstat->Rebin(10); probably shouldn't rebin after we've divided by dpt


  //////Error Bars///////
  hsyst     = new TH1D("hsyst", "hsyst", 100, 0, 10);
  hsystCorr = new TH1D("hsystCorr", "hsystCorr", 100, 0, 10);
  hstat     ->Sumw2(); 
  hsyst     ->Sumw2();
  hsystCorr ->Sumw2(); 
  Int_t nbinx = hstat->GetNbinsX();
  
  for(Int_t ibinx = 1; ibinx <= nbinx; ibinx++){
    hsyst->SetBinContent (ibinx, hstat->GetBinContent(ibinx));
    hsyst->SetBinError   (ibinx, hstat->GetBinContent(ibinx)*0.08);
    hsystCorr->SetBinContent (ibinx, hstat->GetBinContent(ibinx));
    hsystCorr->SetBinError   (ibinx, hstat->GetBinContent(ibinx)*0.15);
  } 

  *SaveAs_Title += "RawSpectrum_O2data_allCountSignalRegion";
  texXtitle = texPtX;
  texYtitle = texRawYield;
}

void PtDistribution_O2data(TH1D* &hstat, TH1D* &hsyst, TH1D*&hsystCorr, Int_t ipart, TFile* file_O2Analysis, TString* SaveAs_Title, TString* &texXtitle, TString* &texYtitle) {

  //O2 file
  // TFile* file_O2Analysis = new TFile("AnalysisResults.root");

  TH1I* H1I_SelectedEventCount = (TH1I*)file_O2Analysis->Get("lambdakzero-analysis/hSelectedEventCounter");
  float SelectedEventCount = H1I_SelectedEventCount->GetEntries();

  TH3D* H3D_DetectedV0s_O2 = (TH3D*)file_O2Analysis->Get("lambdakzero-analysis/h3dMass"+NamePart[ipart]);

  float DCAcut = 1;
  int DCAcut_Bin = H3D_DetectedV0s_O2->GetXaxis()->FindBin(DCAcut);
  cout << "DCA cut bin: " << DCAcut_Bin << endl;

  TH1D* H1D_DetectedV0s_pT_O2 = (TH1D*)H3D_DetectedV0s_O2->ProjectionY(NamehistoInvMass[ipart],1,-1,1,-1);
  H1D_DetectedV0s_pT_O2->Sumw2();
  H1D_DetectedV0s_pT_O2->Scale(1./SelectedEventCount); 
  hstat = (TH1D*)H1D_DetectedV0s_pT_O2->Clone("hstat");

  hstat->Rebin(10);   

  //////Error Bars///////
  hsyst     = new TH1D("hsyst", "hsyst", 100, 0, 10);
  hsystCorr = new TH1D("hsystCorr", "hsystCorr", 100, 0, 10);
  hstat     ->Sumw2(); 
  hsyst     ->Sumw2();
  hsystCorr ->Sumw2(); 
  Int_t nbinx = hstat->GetNbinsX();
  
  for(Int_t ibinx = 1; ibinx <= nbinx; ibinx++){
    hsyst->SetBinContent (ibinx, hstat->GetBinContent(ibinx));
    hsyst->SetBinError   (ibinx, hstat->GetBinContent(ibinx)*0.08);
    hsystCorr->SetBinContent (ibinx, hstat->GetBinContent(ibinx));
    hsystCorr->SetBinError   (ibinx, hstat->GetBinContent(ibinx)*0.15);
  } 

  *SaveAs_Title += "Pt_Distrib_O2data";
  texXtitle = texPtX;
  texYtitle = texRawYield;
}

void InvMass_Plot(TH1D* &hstat, TH1D* &hsyst, TH1D*&hsystCorr, Int_t ipart, TFile* file_O2Analysis, TString* SaveAs_Title, TString* &texXtitle, TString* &texYtitle, bool isMC) {

  TH1I* H1I_SelectedEventCount;
  TH3D* H3D_DetectedV0s_O2;

  if (isMC) {
    H3D_DetectedV0s_O2 = (TH3D*)file_O2Analysis->Get("lambdakzero-analysis-mc/h3dMass"+NamePart[ipart]);
    H1I_SelectedEventCount = (TH1I*)file_O2Analysis->Get("lambdakzero-analysis-mc/hSelectedEventCounter");
  }
  else {
    H3D_DetectedV0s_O2 = (TH3D*)file_O2Analysis->Get("lambdakzero-analysis/h3dMass"+NamePart[ipart]);
    H1I_SelectedEventCount = (TH1I*)file_O2Analysis->Get("lambdakzero-analysis/hSelectedEventCounter");
  }

  float SelectedEventCount = H1I_SelectedEventCount->GetEntries();

  // float DCAcut = 1;
  // int DCAcut_Bin = H3D_DetectedV0s_O2->GetXaxis()->FindBin(DCAcut);
  // cout << "DCA cut bin: " << DCAcut_Bin << endl;

  TH1D* H1D_DetectedV0s_InvMass_O2 = (TH1D*)H3D_DetectedV0s_O2->ProjectionZ(NamehistoInvMass[ipart],0,-1,0,-1);
  H1D_DetectedV0s_InvMass_O2->Sumw2();
  hstat = (TH1D*)H1D_DetectedV0s_InvMass_O2->Clone("hstat");
  // float_t ScalingFactor = SelectedEventCount;
  hstat->Rebin(10);   

  // float_t ScalingFactor = hstat->GetMaximum();  
  // hstat->Scale(1./ScalingFactor);

  //choice of number of bin
  // float nEntries = hstat->GetEntries();
  // float kbin0 = hstat->GetNbinsX();
  // // float ratio_binning = kbin0*1./(sqrt(nEntries)); //sturge
  // float ratio_binning = kbin0*1./(log(nEntries)/log(2)+1); //square root
  // hstat->Rebin(ratio_binning); //sturge  

  //////Error Bars///////
  hsyst     = new TH1D("hsyst", "hsyst", 100, 0, 10);
  hsystCorr = new TH1D("hsystCorr", "hsystCorr", 100, 0, 10);
  hstat     ->Sumw2(); 
  hsyst     ->Sumw2();
  hsystCorr ->Sumw2(); 
  Int_t nbinx = hstat->GetNbinsX();
  
  for(Int_t ibinx = 1; ibinx <= nbinx; ibinx++){
    hsyst->SetBinContent (ibinx, hstat->GetBinContent(ibinx));
    hsyst->SetBinError   (ibinx, hstat->GetBinContent(ibinx)*0.08);
    hsystCorr->SetBinContent (ibinx, hstat->GetBinContent(ibinx));
    hsystCorr->SetBinError   (ibinx, hstat->GetBinContent(ibinx)*0.15);
  } 

  *SaveAs_Title += "InvMass";
  texXtitle = texPtX;
  texYtitle = texCount;
}

void DcaHistoProcessing_O2data(TH1D* &hstat, TH1D* &hsyst, TH1D*&hsystCorr, Int_t ipart, TFile* file_O2Analysis) {

  //O2 file
  // TFile* file_O2Analysis = new TFile("AnalysisResults.root");

  TH1I* H1I_SelectedEventCount = (TH1I*)file_O2Analysis->Get("lambdakzero-analysis/hSelectedEventCounter");
  float SelectedEventCount = H1I_SelectedEventCount->GetEntries();

  TH3D* H3D_DetectedV0s_O2 = (TH3D*)file_O2Analysis->Get("lambdakzero-analysis/h3dMass"+NamePart[ipart]+"Dca");

  // float DCAcut = aa1;
  // int DCAcut_Bin = H3D_DetectedV0s_O2->GetZaxis()->FindBin(DCAcut);
  // cout << "DCA cut bin: " << DCAcut_Bin << endl;

  TH1D* H1D_DetectedV0s_InvMass_O2 = (TH1D*)H3D_DetectedV0s_O2->ProjectionX(NamehistoInvMass[ipart],0,-1,0,-1);
  H1D_DetectedV0s_InvMass_O2->Sumw2();
  hstat = (TH1D*)H1D_DetectedV0s_InvMass_O2->Clone("hstat");
  hstat->Rebin(10);

  //////Error Bars///////
  hsyst     = new TH1D("hsyst", "hsyst", 100, 0, 10);
  hsystCorr = new TH1D("hsystCorr", "hsystCorr", 100, 0, 10);
  hstat     ->Sumw2(); 
  hsyst     ->Sumw2();
  hsystCorr ->Sumw2(); 
  Int_t nbinx = hstat->GetNbinsX();
  
  for(Int_t ibinx = 1; ibinx <= nbinx; ibinx++){
    hsyst->SetBinContent (ibinx, hstat->GetBinContent(ibinx));
    hsyst->SetBinError   (ibinx, hstat->GetBinContent(ibinx)*0.08);
    hsystCorr->SetBinContent (ibinx, hstat->GetBinContent(ibinx));
    hsystCorr->SetBinError   (ibinx, hstat->GetBinContent(ibinx)*0.15);
  } 

}

void PseudoEfficiency_HEPcomparison_allCountSignalRegionOLD_withBugs(TH1D* &hstat, TH1D* &hsyst, TH1D*&hsystCorr, Int_t ipart, TFile* file_O2Analysis, TString* SaveAs_Title, Double_t* pTbins, Int_t nbinpT) {
//////Error Bars initialisation///////
  hsyst     = new TH1D("hsyst", "hsyst", 100, 0, 10);
  hsystCorr = new TH1D("hsystCorr", "hsystCorr", 100, 0, 10);
  hsyst     ->Sumw2(); 
  hsystCorr ->Sumw2(); 

  //HEP file
  TFile* file_HEPAnalysis = new TFile("HEPData-ins881474-v1-root.root");
  TH1D* hPtDiff_RawSpectrum_HEP[numPart]; //= (TH1D*)file_HEPAnalysis->Get("Table 1/Hist1D_y1");
  TH1D* hPtDiff_RawSpectrum_HEP_statErrors[numPart]; // = (TH1D*)file_HEPAnalysis->Get("Table 1/Hist1D_y1_e1");
  hPtDiff_RawSpectrum_HEP[0] = (TH1D*)file_HEPAnalysis->Get("Table 1/Hist1D_y1"); //K0S
  hPtDiff_RawSpectrum_HEP_statErrors[0] = (TH1D*)file_HEPAnalysis->Get("Table 1/Hist1D_y1_e1");
  hPtDiff_RawSpectrum_HEP[1] = (TH1D*)file_HEPAnalysis->Get("Table 2/Hist1D_y1"); //Lambda
  hPtDiff_RawSpectrum_HEP_statErrors[1] = (TH1D*)file_HEPAnalysis->Get("Table 2/Hist1D_y1_e1");
  hPtDiff_RawSpectrum_HEP[2] = (TH1D*)file_HEPAnalysis->Get("Table 3/Hist1D_y1"); //AntiLambda
  hPtDiff_RawSpectrum_HEP_statErrors[2] = (TH1D*)file_HEPAnalysis->Get("Table 3/Hist1D_y1_e1");


  //O2 file
  // TFile* file_O2Analysis = new TFile("AnalysisResults_WithSel8Cut.root");

  TH1I* H1I_SelectedEventCount = (TH1I*)file_O2Analysis->Get("lambdakzero-analysis/hSelectedEventCounter");
  float SelectedEventCount = H1I_SelectedEventCount->GetEntries();

  TH3D* H3D_DetectedV0s_O2 = (TH3D*)file_O2Analysis->Get("lambdakzero-analysis/h3dMass"+NamePart[ipart]);
  TH1D* H1D_DetectedV0s_PtAxis_O2 = (TH1D*)H3D_DetectedV0s_O2->ProjectionY(NamehistoInvMass[ipart],0,-1,0,-1);
  
  TH1D* hPtDiff_RawSpectrum_O2 = (TH1D*)H1D_DetectedV0s_PtAxis_O2->Clone("hPtDiff_RawSpectrum_O2"); 
  hPtDiff_RawSpectrum_O2->Reset("M"); //should be empty, only used that to get same Pt binning as used in analysis
  
  //Mass window in which we count K0S
  float MassSignalInterval_lowEdge = 0.48;
  float MassSignalInterval_upEdge = 0.5;
  int MassSignalInterval_lowEdgeBin = H3D_DetectedV0s_O2->GetZaxis()->FindBin(MassSignalInterval_lowEdge);
  int MassSignalInterval_upEdgeBin = H3D_DetectedV0s_O2->GetZaxis()->FindBin(MassSignalInterval_upEdge);

  Int_t nbinPt = hPtDiff_RawSpectrum_O2->GetNbinsX();
  for(Int_t ibinPt = 1; ibinPt <= nbinPt; ibinPt++){
    float SignalCountRough = H3D_DetectedV0s_O2->Integral(0,-1,ibinPt,ibinPt+1,MassSignalInterval_lowEdgeBin,MassSignalInterval_upEdgeBin);

    hPtDiff_RawSpectrum_O2->SetBinContent(ibinPt,SignalCountRough); //hPtDiff_RawSpectrum_O2 is actually a count, not the 1/Nev*dN/dpT
  }

  Double_t xbins[17] = {0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.6, 2.0, 2.4, 3.0};
  TH1D* hPtDiff_RawSpectrum_O2_rebinned = (TH1D*)hPtDiff_RawSpectrum_O2->Rebin(16,"hPtDiff_RawSpectrum_O2_rebinned",xbins);
  hPtDiff_RawSpectrum_O2_rebinned->Sumw2();

  double placeholder = pTbins[0];
  int placeholder2 = nbinpT;

  // ->Divide(hPtDiff_RawSpectrum_HEP);
  Int_t nbinpT_2 = hPtDiff_RawSpectrum_O2_rebinned->GetNbinsX();
  for(Int_t ibinPt = 1; ibinPt <= nbinpT_2; ibinPt++){
    float SignalCount_ibinPt = hPtDiff_RawSpectrum_O2_rebinned->GetBinContent(ibinPt);
    float dpT = hPtDiff_RawSpectrum_O2_rebinned->GetXaxis()->GetBinWidth(ibinPt);
    // float ratio_HEP_O2 = (SignalCount_ibinPt*1./dpT)*1./hPtDiff_RawSpectrum_HEP->GetBinContent(ibinPt);// = dN/dpT_O2 / dN/dpT_HEP)
    hPtDiff_RawSpectrum_O2_rebinned->SetBinContent(ibinPt,SignalCount_ibinPt*1./dpT);

    hPtDiff_RawSpectrum_O2_rebinned->SetBinError(ibinPt,sqrt(SignalCount_ibinPt*1./dpT));//*1./hPtDiff_RawSpectrum_HEP->GetBinContent(ibinPt)
    hPtDiff_RawSpectrum_HEP[ipart]->SetBinError(ibinPt,hPtDiff_RawSpectrum_HEP_statErrors[ipart]->GetBinContent(ibinPt));
  }

// 'MIGHT BE MISSING a division by rapidity interval for HEP DATA because it is dN**2/dpT/dy and rapidity interval is [-0.75,0.75]'

  hstat = (TH1D*)hPtDiff_RawSpectrum_O2->Rebin(16,"hstat",xbins);
  hPtDiff_RawSpectrum_O2_rebinned->Scale(1./SelectedEventCount);
  hstat->Reset("M");
  hstat->Divide(hPtDiff_RawSpectrum_O2_rebinned,hPtDiff_RawSpectrum_HEP[ipart]);

  //Error Bars Processing
  Int_t nbinx = hstat->GetNbinsX();

  for(Int_t ibinx = 1; ibinx <= nbinx; ibinx++){
    hsyst->SetBinContent (ibinx, hstat->GetBinContent(ibinx));
    hsyst->SetBinError   (ibinx, hstat->GetBinContent(ibinx)*0.08);
    hsystCorr->SetBinContent (ibinx, hstat->GetBinContent(ibinx));
    hsystCorr->SetBinError   (ibinx, hstat->GetBinContent(ibinx)*0.15);
  } 

  *SaveAs_Title += "PseudoEfficiency_HEPcomparison_allCountSignalRegion";
}

void PseudoEfficiency_HEPcomparison_allCountSignalRegion(TH1D* &hstat, TH1D* &hsyst, TH1D*&hsystCorr, Int_t ipart, TFile* file_O2Analysis, TString* SaveAs_Title, TString* &texXtitle, TString* &texYtitle, Double_t* pTbins, Int_t nbinpT) {
//////Error Bars initialisation///////
  hsyst     = new TH1D("hsyst", "hsyst", 100, 0, 10);
  hsystCorr = new TH1D("hsystCorr", "hsystCorr", 100, 0, 10);
  hsyst     ->Sumw2(); 
  hsystCorr ->Sumw2(); 

  //HEP file
  TFile* file_HEPAnalysis = new TFile("HEPData-ins881474-v1-root.root");
  TH1D* hPtDiff_RawSpectrum_HEP[numPart]; //= (TH1D*)file_HEPAnalysis->Get("Table 1/Hist1D_y1");
  TH1D* hPtDiff_RawSpectrum_HEP_statErrors[numPart]; // = (TH1D*)file_HEPAnalysis->Get("Table 1/Hist1D_y1_e1");
  hPtDiff_RawSpectrum_HEP[0] = (TH1D*)file_HEPAnalysis->Get("Table 1/Hist1D_y1"); //K0S
  hPtDiff_RawSpectrum_HEP_statErrors[0] = (TH1D*)file_HEPAnalysis->Get("Table 1/Hist1D_y1_e1");
  hPtDiff_RawSpectrum_HEP[1] = (TH1D*)file_HEPAnalysis->Get("Table 2/Hist1D_y1"); //Lambda
  hPtDiff_RawSpectrum_HEP_statErrors[1] = (TH1D*)file_HEPAnalysis->Get("Table 2/Hist1D_y1_e1");
  hPtDiff_RawSpectrum_HEP[2] = (TH1D*)file_HEPAnalysis->Get("Table 3/Hist1D_y1"); //AntiLambda
  hPtDiff_RawSpectrum_HEP_statErrors[2] = (TH1D*)file_HEPAnalysis->Get("Table 3/Hist1D_y1_e1");


  //O2 file
  // TFile* file_O2Analysis = new TFile("AnalysisResults_WithSel8Cut.root");

  TH1I* H1I_SelectedEventCount = (TH1I*)file_O2Analysis->Get("lambdakzero-analysis/hSelectedEventCounter");
  float SelectedEventCount = H1I_SelectedEventCount->GetEntries();

  TH3D* H3D_DetectedV0s_O2 = (TH3D*)file_O2Analysis->Get("lambdakzero-analysis/h3dMass"+NamePart[ipart]);
  TH1D* H1D_DetectedV0s_PtAxis_O2 = (TH1D*)H3D_DetectedV0s_O2->ProjectionY(NamehistoInvMass[ipart],0,-1,0,-1);
  
  // TH1D* hPtDiff_RawSpectrum_O2 = (TH1D*)H1D_DetectedV0s_PtAxis_O2->Clone("hPtDiff_RawSpectrum_O2"); 
  // hPtDiff_RawSpectrum_O2->Reset("M"); //should be empty, only used that to get same Pt binning as used in analysis
  
  //Mass window in which we count K0S
  float MassSignalInterval_lowEdge = 0.48;//45
  float MassSignalInterval_upEdge = 0.5;//55 gives a sudden much higher yield; is there overflow/underflow at work there?
  int MassSignalInterval_lowEdgeBin = H3D_DetectedV0s_O2->GetZaxis()->FindBin(MassSignalInterval_lowEdge);
  int MassSignalInterval_upEdgeBin = H3D_DetectedV0s_O2->GetZaxis()->FindBin(MassSignalInterval_upEdge);

  // aaa TH1D* H1D_DetectedV0s_O2_ProjZ = (TH1D*)H3D_DetectedV0s_O2_rebinnedZ->ProjectionZ("InvMass_"+NamePart[ipart]+Form("_Pt%i", ibinPt),0,nx,ibinPt_originalAxis_low,ibinPt_originalAxis_up);

  TH1D *H1D_DetectedV0s_O2_SmallpTinterval[nbinpT];
  TH1D *hPtDiff_RawSpectrum_O2 = new TH1D("hPtDiff_RawSpectrum_O2","hPtDiff_RawSpectrum_O2",nbinpT,pTbins);
  // Int_t nbinPt = hPtDiff_RawSpectrum_O2->GetNbinsX();

  for(int ibinPt = 1; ibinPt <= nbinpT; ibinPt++){
    int iPt = ibinPt-1; //different convention for bins that start at 1 and arrays that start at 0

    int ibinPt_originalAxis_low = H3D_DetectedV0s_O2->GetYaxis()->FindBin(pTbins[iPt]);
    int ibinPt_originalAxis_up = H3D_DetectedV0s_O2->GetYaxis()->FindBin(pTbins[iPt+1]) - 1;
    cout << "ibinPt: " << ibinPt << "; originalAxis [low;up] = [" << ibinPt_originalAxis_low << ";" << ibinPt_originalAxis_up << "]" << endl;

    H1D_DetectedV0s_O2_SmallpTinterval[iPt] = (TH1D*)H3D_DetectedV0s_O2->ProjectionZ("InvMass_"+NamePart[ipart]+Form("_Pt%i", ibinPt),0,-1,ibinPt_originalAxis_low,ibinPt_originalAxis_up);

    // for(Int_t ibinPt = 1; ibinPt <= nbinPt; ibinPt++){

    float SignalCountRough = H1D_DetectedV0s_O2_SmallpTinterval[iPt]->Integral(MassSignalInterval_lowEdgeBin,MassSignalInterval_upEdgeBin);

    hPtDiff_RawSpectrum_O2->SetBinContent(ibinPt,SignalCountRough); //hPtDiff_RawSpectrum_O2 is actually a count, not the 1/Nev*dN/dpT

    double dpT = hPtDiff_RawSpectrum_O2->GetXaxis()->GetBinWidth(ibinPt);
    float drapidity = 1.5;//1.5
    double dN_dpT = SignalCountRough *1./dpT*1./drapidity;
    hPtDiff_RawSpectrum_O2->SetBinContent(ibinPt,dN_dpT);
    hPtDiff_RawSpectrum_O2->SetBinError(ibinPt,sqrt(SignalCountRough*1./dpT*1./drapidity)); //Issue with errors; this error of the O2 data spectrum is false, not right formula with fit
    cout << "(ibinPt,dN_dpT) = (" << ibinPt << "," << dN_dpT*1./SelectedEventCount << ")" << endl;

    hPtDiff_RawSpectrum_HEP[ipart]->SetBinError(ibinPt,hPtDiff_RawSpectrum_HEP_statErrors[ipart]->GetBinContent(ibinPt));
  }

  // Double_t xbins[17] = {0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.6, 2.0, 2.4, 3.0};
  // TH1D* hPtDiff_RawSpectrum_O2_rebinned = (TH1D*)hPtDiff_RawSpectrum_O2->Rebin(16,"hPtDiff_RawSpectrum_O2_rebinned",xbins);
  // hPtDiff_RawSpectrum_O2_rebinned->Sumw2();

  // double placeholder = pTbins[0];
  // int placeholder2 = nbinpT;

  // // ->Divide(hPtDiff_RawSpectrum_HEP);
  // Int_t nbinpT_2 = hPtDiff_RawSpectrum_O2_rebinned->GetNbinsX();
  // for(Int_t ibinPt = 1; ibinPt <= nbinpT_2; ibinPt++){
  //   float SignalCount_ibinPt = hPtDiff_RawSpectrum_O2_rebinned->GetBinContent(ibinPt);
  //   float dpT = hPtDiff_RawSpectrum_O2_rebinned->GetXaxis()->GetBinWidth(ibinPt);
  //   // float ratio_HEP_O2 = (SignalCount_ibinPt*1./dpT)*1./hPtDiff_RawSpectrum_HEP->GetBinContent(ibinPt);// = dN/dpT_O2 / dN/dpT_HEP)
  //   hPtDiff_RawSpectrum_O2_rebinned->SetBinContent(ibinPt,SignalCount_ibinPt*1./dpT);

  //   hPtDiff_RawSpectrum_O2_rebinned->SetBinError(ibinPt,sqrt(SignalCount_ibinPt*1./dpT));//*1./hPtDiff_RawSpectrum_HEP->GetBinContent(ibinPt)
  //   hPtDiff_RawSpectrum_HEP[ipart]->SetBinError(ibinPt,hPtDiff_RawSpectrum_HEP_statErrors[ipart]->GetBinContent(ibinPt));
  // }

// 'MIGHT BE MISSING a division by rapidity interval for HEP DATA because it is dN**2/dpT/dy and rapidity interval is [-0.75,0.75]'

  hstat = (TH1D*)hPtDiff_RawSpectrum_O2->Clone("hstat");
  hPtDiff_RawSpectrum_O2->Scale(1./SelectedEventCount);
  hstat->Reset("M");
  hstat->Divide(hPtDiff_RawSpectrum_O2,hPtDiff_RawSpectrum_HEP[ipart]);

  //Error Bars Processing
  Int_t nbinx = hstat->GetNbinsX();

  for(Int_t ibinx = 1; ibinx <= nbinx; ibinx++){
    hsyst->SetBinContent (ibinx, hstat->GetBinContent(ibinx));
    hsyst->SetBinError   (ibinx, hstat->GetBinContent(ibinx)*0.08);
    hsystCorr->SetBinContent (ibinx, hstat->GetBinContent(ibinx));
    hsystCorr->SetBinError   (ibinx, hstat->GetBinContent(ibinx)*0.15);
  } 

  *SaveAs_Title += "PseudoEfficiency_HEPcomparison_allCountSignalRegion";
  texXtitle = texPtX;
  texYtitle = texptDifferentialYield_HEPratio_pseudoEff;
}



Double_t fparab(Double_t *x, Double_t *par) {
  const Int_t numPart=7;
  // Float_t liminf[numPart]={0.48, 1.11, 1.11, 1.31,  1.31,  1.665, 1.665};
  // Float_t limsup[numPart]={0.54, 1.13, 1.13, 1.335, 1.335, 1.685, 1.685}; Original
  // Float_t limsup[numPart]={0.515, 1.13, 1.13, 1.335, 1.335, 1.685, 1.685}; //Tweaked
  Int_t part=par[3];
  if (reject && x[0] > min_range_signal[part] && x[0] < max_range_signal[part]) { //Original was using this duplicate but different definition of limsup and liminf
    TF1::RejectPoint();
    return 0;
  }
  return par[0] + par[1]*x[0] + par[2]*x[0]*x[0];
}
Double_t fline(Double_t *x, Double_t *par) {
  const Int_t numPart=7;
  // Float_t liminf[numPart]={0.48, 1.11, 1.11, 1.31,  1.31,  1.665, 1.665};
  // Float_t limsup[numPart]={0.54, 1.13, 1.13, 1.335, 1.335, 1.685, 1.685}; Original
  // Float_t limsup[numPart]={0.515, 1.13, 1.13, 1.335, 1.335, 1.685, 1.685}; //Tweaked
  Int_t part=par[2];
  if (reject && x[0] > min_range_signal[part] && x[0] < max_range_signal[part]) { //Original was using this duplicate but different definition of limsup and liminf
    TF1::RejectPoint();
    return 0;
  }
  return par[0] + par[1]*x[0];
}

void RawSpectrum_O2data_withFit(TH1D* &hstat, TH1D* &hsyst, TH1D*&hsystCorr, Int_t ipart, TFile* file_O2Analysis, TString* SaveAs_Title, TString* &texXtitle, TString* &texYtitle, Double_t* pTbins, Int_t nbinpT) {

  //O2 file
  TH1I* H1I_SelectedEventCount = (TH1I*)file_O2Analysis->Get("lambdakzero-analysis/hSelectedEventCounter");
  float SelectedEventCount = H1I_SelectedEventCount->GetEntries();

  TH3D* H3D_DetectedV0s_O2 = (TH3D*)file_O2Analysis->Get("lambdakzero-analysis/h3dMass"+NamePart[ipart]);
  TH1D* H1D_DetectedV0s_PtAxis_O2 = (TH1D*)H3D_DetectedV0s_O2->ProjectionY(NamehistoInvMass[ipart],0,-1,0,-1);

  //Rebinning InvMass axis if necessary
  TH3D* H3D_DetectedV0s_O2_rebinnedZ = (TH3D*)H3D_DetectedV0s_O2->RebinZ(1,"H3D_DetectedV0s_O2_rebinnedZ");

  //Fit tools initialisation
  TF1 *gauss[nbinpT];
  TF1 *GaussPlusPolynom[nbinpT];
  TF1 *bkgparab[nbinpT];
  Float_t fFitResult[nbinpT];
  Double_t parGaussParab[nbinpT][6];  

  TH1D *H1D_DetectedV0s_O2_SmallpTinterval[nbinpT];
  TH1D *hPtDiff_RawSpectrum_O2 = new TH1D("hPtDiff_RawSpectrum_O2","hPtDiff_RawSpectrum_O2",nbinpT,pTbins);

  TCanvas *canvasMassvsPt = new TCanvas ("canvasMassPt", NamehistoInvMass[ipart], 1200, 800);
  canvasMassvsPt->Divide(4,4);

  for(int ibinPt = 1; ibinPt <= nbinpT; ibinPt++){
    int iPt = ibinPt-1; //different convention for bins that start at 1 and arrays that start at 0

    canvasMassvsPt->cd(ibinPt); //ATTENTION: +1 because canvas divide starts counting from 1

    int ibinPt_originalAxis_low = H3D_DetectedV0s_O2->GetYaxis()->FindBin(pTbins[iPt]);
    int ibinPt_originalAxis_up = H3D_DetectedV0s_O2->GetYaxis()->FindBin(pTbins[iPt+1]) - 1;
    cout << "ibinPt: " << ibinPt << "; originalAxis [low;up] = [" << ibinPt_originalAxis_low << ";" << ibinPt_originalAxis_up << "]" << endl;

    H1D_DetectedV0s_O2_SmallpTinterval[iPt] = (TH1D*)H3D_DetectedV0s_O2_rebinnedZ->ProjectionZ("InvMass_"+NamePart[ipart]+Form("_Pt%i", ibinPt),0,-1,ibinPt_originalAxis_low,ibinPt_originalAxis_up);

    gauss[iPt] = new TF1("gauss_"+NamePart[ipart]+Form("_Pt%i", ibinPt), "gaus", min_range_signal[ipart], max_range_signal[ipart]);
    gauss[iPt]->SetLineColor(kRed);   
    gauss[iPt]->SetParName(0, "norm");
    gauss[iPt]->SetParName(1, "mean");
    gauss[iPt]->SetParName(2, "sigma");
    gauss[iPt]->SetParameters(1., MassPart[ipart], 0.0002);
    gauss[iPt]->SetParLimits(0, 0., 1.1*H1D_DetectedV0s_O2_SmallpTinterval[iPt]->GetBinContent(H1D_DetectedV0s_O2_SmallpTinterval[iPt]->GetMaximumBin()));
    gauss[iPt]->SetParLimits(1, min_range_signal[ipart], max_range_signal[ipart]);
    gauss[iPt]->SetParLimits(2, 0.0001, 0.05);

    bkgparab[iPt] = new TF1("parab_"+NamePart[ipart]+Form("_Pt%i", ibinPt), fparab, liminf[ipart], limsup[ipart], 4);
    bkgparab[iPt]->SetLineColor(kGreen);
    bkgparab[iPt]->FixParameter(3, ipart);

    GaussPlusPolynom[iPt] = new TF1("GaussPlusPolynom_"+NamePart[ipart]+Form("_Pt%i", ibinPt), "gaus(0)+pol2(3)", liminf[ipart], limsup[ipart]);
    GaussPlusPolynom[iPt]->SetLineColor(kAzure); 
    GaussPlusPolynom[iPt]->SetParName(0, "norm");
    GaussPlusPolynom[iPt]->SetParName(1, "mean");
    GaussPlusPolynom[iPt]->SetParName(2, "sigma");
    GaussPlusPolynom[iPt]->SetParName(3, "p0");
    GaussPlusPolynom[iPt]->SetParName(4, "p1");
    GaussPlusPolynom[iPt]->SetParName(5, "p2");

    GaussPlusPolynom[iPt]->SetParameters(1., MassPart[ipart], 0.0002);
    GaussPlusPolynom[iPt]->SetParLimits(0, 0., 1.1*H1D_DetectedV0s_O2_SmallpTinterval[iPt]->GetBinContent(H1D_DetectedV0s_O2_SmallpTinterval[iPt]->GetMaximumBin()));
    GaussPlusPolynom[iPt]->SetParLimits(1, min_range_signal[ipart], max_range_signal[ipart]);
    GaussPlusPolynom[iPt]->SetParLimits(2, 0.0001, 0.05);

    H1D_DetectedV0s_O2_SmallpTinterval[iPt]->Fit(gauss[iPt], "R0QL");
    reject = true; // for background fit: ignore signal range when fitting
    H1D_DetectedV0s_O2_SmallpTinterval[iPt]->Fit(bkgparab[iPt], "R0QL");
    reject = false;

    gauss[iPt]->GetParameters(&parGaussParab[iPt][0]);
    bkgparab[iPt]->GetParameters(&parGaussParab[iPt][3]);
    GaussPlusPolynom[iPt]->SetParameters(parGaussParab[iPt]);

    fFitResult[iPt] = H1D_DetectedV0s_O2_SmallpTinterval[iPt]->Fit(GaussPlusPolynom[iPt],"SRB+QL"); //L option is log likelyhood: would supposedly be better in case of counts

    //give background the final parameters from that last signal+background fit
    GaussPlusPolynom[iPt]->GetParameters(&parGaussParab[iPt][0]);
    bkgparab[iPt]->SetParameters(&parGaussParab[iPt][3]);

    //plot Inv Mass histogram and the background and signal fits for each iPt bin
    H1D_DetectedV0s_O2_SmallpTinterval[iPt]->GetXaxis()->SetTitle("inv. mass (GeV/c^{2})");
    H1D_DetectedV0s_O2_SmallpTinterval[iPt]->GetYaxis()->SetTitle("Counts");
    gPad->DrawFrame(LowerLimitDisplay_Mass[ipart],0,UpperLimitDisplay_Mass[ipart],1.2*H1D_DetectedV0s_O2_SmallpTinterval[iPt]->GetMaximum());
    H1D_DetectedV0s_O2_SmallpTinterval[iPt]->Draw("same");
    bkgparab[iPt]->Draw("same");
    // gauss[iPt]->Draw("same");
    TLegend *legt = new TLegend( 0.62, 0.6, 0.87, 0.8, Form("pT [%.1f;%.1f]", pTbins[iPt], pTbins[iPt+1]));
    legt->Draw();


    cout << "xmin_InvMass= " << liminf[ipart] << ", xmax_InvMass= " << limsup[ipart] << endl;
    int ibin_liminf = H1D_DetectedV0s_O2_SmallpTinterval[iPt]->GetXaxis()->FindBin(liminf[ipart]);
    int ibin_limsup = H1D_DetectedV0s_O2_SmallpTinterval[iPt]->GetXaxis()->FindBin(limsup[ipart]);
    int nbinInvMass = ibin_limsup - ibin_liminf;
    cout << "ibin_liminf= " << ibin_liminf << ", ibin_limsup= " << ibin_limsup << endl;
    cout << "nbinInvMass = " << nbinInvMass << endl;

    int totalCount = H1D_DetectedV0s_O2_SmallpTinterval[iPt]->Integral(ibin_liminf+1,ibin_limsup-1); //integral of a hist counts the number of bins, does not care for bin width: not true integral
    double backgroundCount = bkgparab[iPt]->Integral(liminf[ipart],limsup[ipart])*1./((limsup[ipart]-liminf[ipart])/nbinInvMass); //bkgparab is not an histogram but a function; here we need to convert it to a count

    float SignalCountFit = totalCount - backgroundCount;
    cout << "totalCount= " << totalCount << ", backgroundCount= " << backgroundCount << endl;
    cout << "SignalCountFit= " << SignalCountFit << endl;

    float dpT = hPtDiff_RawSpectrum_O2->GetXaxis()->GetBinWidth(ibinPt);
    float dN_dpT = SignalCountFit *1./dpT;
    hPtDiff_RawSpectrum_O2->SetBinContent(ibinPt,dN_dpT);
    hPtDiff_RawSpectrum_O2->SetBinError(ibinPt,sqrt(SignalCountFit*1./dpT));
    cout << "(ibinPt,dN_dpT) = (" << ibinPt << "," << dN_dpT*1./SelectedEventCount << ")" << endl;
  }
  canvasMassvsPt->SaveAs("Fits_invMass_ptDiff_"+NamePart[ipart]+".pdf");

  // hPtDiff_RawSpectrum_O2->Scale(1./SelectedEventCount);
  hstat = (TH1D*)hPtDiff_RawSpectrum_O2->Clone("hstat");

  //////Error Bars///////
  hsyst     = new TH1D("hsyst", "hsyst", 100, 0, 10);
  hsystCorr = new TH1D("hsystCorr", "hsystCorr", 100, 0, 10);
  hstat     ->Sumw2(); 
  hsyst     ->Sumw2();
  hsystCorr ->Sumw2(); 
  Int_t nbinx = hstat->GetNbinsX();
  
  for(Int_t ibinx = 1; ibinx <= nbinx; ibinx++){
    hsyst->SetBinContent (ibinx, hstat->GetBinContent(ibinx));
    hsyst->SetBinError   (ibinx, hstat->GetBinContent(ibinx)*0.08);
    hsystCorr->SetBinContent (ibinx, hstat->GetBinContent(ibinx));
    hsystCorr->SetBinError   (ibinx, hstat->GetBinContent(ibinx)*0.15);
  } 

  *SaveAs_Title += "RawSpectrum_O2data_withFit";
  texXtitle = texPtX;
  texYtitle = texRawYield;
}

void PseudoEfficiency_HEPcomparison_withFit(TH1D* &hstat, TH1D* &hsyst, TH1D* &hsystCorr, Int_t ipart, TFile* file_O2Analysis, TString* SaveAs_Title, TString* &texXtitle, TString* &texYtitle, Double_t* pTbins, Int_t nbinpT) {

  //HEP file
  TFile* file_HEPAnalysis = new TFile("HEPData-ins881474-v1-root.root");
  TH1D* hPtDiff_RawSpectrum_HEP[numPart]; //= (TH1D*)file_HEPAnalysis->Get("Table 1/Hist1D_y1");
  TH1D* hPtDiff_RawSpectrum_HEP_statErrors[numPart]; // = (TH1D*)file_HEPAnalysis->Get("Table 1/Hist1D_y1_e1");
  hPtDiff_RawSpectrum_HEP[0] = (TH1D*)file_HEPAnalysis->Get("Table 1/Hist1D_y1"); //K0S
  hPtDiff_RawSpectrum_HEP_statErrors[0] = (TH1D*)file_HEPAnalysis->Get("Table 1/Hist1D_y1_e1");
  hPtDiff_RawSpectrum_HEP[1] = (TH1D*)file_HEPAnalysis->Get("Table 2/Hist1D_y1"); //Lambda
  hPtDiff_RawSpectrum_HEP_statErrors[1] = (TH1D*)file_HEPAnalysis->Get("Table 2/Hist1D_y1_e1");
  hPtDiff_RawSpectrum_HEP[2] = (TH1D*)file_HEPAnalysis->Get("Table 3/Hist1D_y1"); //AntiLambda
  hPtDiff_RawSpectrum_HEP_statErrors[2] = (TH1D*)file_HEPAnalysis->Get("Table 3/Hist1D_y1_e1");

  //O2 file
  TH1I* H1I_SelectedEventCount = (TH1I*)file_O2Analysis->Get("lambdakzero-analysis/hSelectedEventCounter");
  int SelectedEventCount = H1I_SelectedEventCount->GetEntries();

  TH3D* H3D_DetectedV0s_O2 = (TH3D*)file_O2Analysis->Get("lambdakzero-analysis/h3dMass"+NamePart[ipart]);

  TH1D* H1D_DetectedV0s_PtAxis_O2 = (TH1D*)H3D_DetectedV0s_O2->ProjectionY(NamehistoInvMass[ipart],0,-1,0,-1);

  //Rebinning InvMass axis if necessary
  TH3D* H3D_DetectedV0s_O2_rebinnedZ = (TH3D*)H3D_DetectedV0s_O2->RebinZ(10,"H3D_DetectedV0s_O2_rebinnedZ");

  //Fit tools initialisation
  TF1 *gauss[nbinpT];
  TF1 *GaussPlusPolynom[nbinpT];
  TF1 *bkgparab[nbinpT];
  
  Float_t fFitResult[nbinpT];
  Double_t parGaussParab[nbinpT][6];  
  Double_t parGaussParabError_mu[nbinpT];  
  Double_t parGaussParabError_sigma[nbinpT];  

  TH1D *H1D_DetectedV0s_O2_SmallpTinterval[nbinpT];
  TH1D *hPtDiff_RawSpectrum_O2 = new TH1D("hPtDiff_RawSpectrum_O2","hPtDiff_RawSpectrum_O2",nbinpT,pTbins);

  // TCanvas *canvasMassvsPt = new TCanvas ("canvasMassPt", NamehistoInvMass[ipart], 800, 600);
  TCanvas *canvasMassvsPt = new TCanvas ("canvasMassPt", NamehistoInvMass[ipart], 1600, 1200);
  canvasMassvsPt->Divide(std::ceil(std::sqrt(nbinpT)),std::ceil(std::sqrt(nbinpT)));
  TH1 *hFrame[nbinpT];

  for(int ibinPt = 1; ibinPt <= nbinpT; ibinPt++){
    int iPt = ibinPt-1; //different convention for bins that start at 1 and arrays that start at 0

    canvasMassvsPt->cd(ibinPt); //ATTENTION: +1 because canvas divide starts counting from 1

    int ibinPt_originalAxis_low = H3D_DetectedV0s_O2->GetYaxis()->FindBin(pTbins[iPt]);
    int ibinPt_originalAxis_up = H3D_DetectedV0s_O2->GetYaxis()->FindBin(pTbins[iPt+1]) - 1;
    cout << "ibinPt: " << ibinPt << "; originalAxis [low;up] = [" << ibinPt_originalAxis_low << ";" << ibinPt_originalAxis_up << "]" << endl;

    H1D_DetectedV0s_O2_SmallpTinterval[iPt] = (TH1D*)H3D_DetectedV0s_O2_rebinnedZ->ProjectionZ("InvMass_"+NamePart[ipart]+Form("_Pt%i", ibinPt),0,-1,ibinPt_originalAxis_low,ibinPt_originalAxis_up);
    hFrame[ibinPt] = canvasMassvsPt->DrawFrame(LowerLimitDisplay_Mass[ipart],0,UpperLimitDisplay_Mass[ipart],1.45*H1D_DetectedV0s_O2_SmallpTinterval[iPt]->GetMaximum());
    hFrame[ibinPt]->SetXTitle(texInvMassDecays_titles[ipart]->Data());
    hFrame[ibinPt]->SetYTitle("Counts");
    H1D_DetectedV0s_O2_SmallpTinterval[iPt]->Draw("E,same");

    gauss[iPt] = new TF1("gauss_"+NamePart[ipart]+Form("_Pt%i", ibinPt), "gaus", min_range_signal[ipart], max_range_signal[ipart]);
    gauss[iPt]->SetLineColor(kBlue);   
    gauss[iPt]->SetParName(0, "norm");
    gauss[iPt]->SetParName(1, "mean");
    gauss[iPt]->SetParName(2, "sigma");
    gauss[iPt]->SetParameters(1., MassPart[ipart], 0.0002);
    gauss[iPt]->SetParLimits(0, 0., 1.1*H1D_DetectedV0s_O2_SmallpTinterval[iPt]->GetBinContent(H1D_DetectedV0s_O2_SmallpTinterval[iPt]->GetMaximumBin()));
    gauss[iPt]->SetParLimits(1, min_range_signal[ipart], max_range_signal[ipart]);
    gauss[iPt]->SetParLimits(2, 0.0001, 0.05);

    bkgparab[iPt] = new TF1("parab_"+NamePart[ipart]+Form("_Pt%i", ibinPt), fparab, liminf[ipart], limsup[ipart], 4);
    bkgparab[iPt]->SetLineColor(kGreen);
    bkgparab[iPt]->FixParameter(3, ipart);

    GaussPlusPolynom[iPt] = new TF1("GaussPlusPolynom_"+NamePart[ipart]+Form("_Pt%i", ibinPt), "gaus(0)+pol2(3)", liminf[ipart], limsup[ipart]);
    GaussPlusPolynom[iPt]->SetLineColor(kRed); 
    GaussPlusPolynom[iPt]->SetParName(0, "norm");
    GaussPlusPolynom[iPt]->SetParName(1, "mean");
    GaussPlusPolynom[iPt]->SetParName(2, "sigma");
    GaussPlusPolynom[iPt]->SetParName(3, "p0");
    GaussPlusPolynom[iPt]->SetParName(4, "p1");
    GaussPlusPolynom[iPt]->SetParName(5, "p2");

    GaussPlusPolynom[iPt]->SetParameters(1., MassPart[ipart], 0.0002);
    GaussPlusPolynom[iPt]->SetParLimits(0, 0., 1.1*H1D_DetectedV0s_O2_SmallpTinterval[iPt]->GetBinContent(H1D_DetectedV0s_O2_SmallpTinterval[iPt]->GetMaximumBin()));
    GaussPlusPolynom[iPt]->SetParLimits(1, min_range_signal[ipart], max_range_signal[ipart]);
    GaussPlusPolynom[iPt]->SetParLimits(2, 0.0001, 0.05);

    H1D_DetectedV0s_O2_SmallpTinterval[iPt]->Fit(gauss[iPt], "R0QL");
    reject = true; // for background fit: ignore signal range when fitting
    H1D_DetectedV0s_O2_SmallpTinterval[iPt]->Fit(bkgparab[iPt], "R0QL");
    reject = false; // to allow plotting it without the gap

    gauss[iPt]->GetParameters(&parGaussParab[iPt][0]);
    bkgparab[iPt]->GetParameters(&parGaussParab[iPt][3]);
    GaussPlusPolynom[iPt]->SetParameters(parGaussParab[iPt]);

    fFitResult[iPt] = H1D_DetectedV0s_O2_SmallpTinterval[iPt]->Fit(GaussPlusPolynom[iPt],"SRB+QL"); //L option is log likelyhood: would supposedly be better in case of counts

    //give background the final parameters from that last signal+background fit
    GaussPlusPolynom[iPt]->GetParameters(&parGaussParab[iPt][0]);
    bkgparab[iPt]->SetParameters(&parGaussParab[iPt][3]);

    parGaussParabError_mu[iPt] = GaussPlusPolynom[iPt]->GetParError(0);
    parGaussParabError_sigma[iPt] = GaussPlusPolynom[iPt]->GetParError(1);

    //plot Inv Mass histogram and the background and signal fits for each iPt bin
    // gPad->DrawFrame(LowerLimitDisplay_Mass[ipart],0,UpperLimitDisplay_Mass[ipart],1.45*H1D_DetectedV0s_O2_SmallpTinterval[iPt]->GetMaximum());
    // gPad->SetXTitle(texInvMassDecays_titles[ipart]->Data());
    // gPad->SetYTitle("Counts");
    // H1D_DetectedV0s_O2_SmallpTinterval[iPt]->Draw("E,same");
    H1D_DetectedV0s_O2_SmallpTinterval[iPt]->GetXaxis()->SetTitle(texInvMassDecays_titles[ipart]->Data());
    H1D_DetectedV0s_O2_SmallpTinterval[iPt]->GetYaxis()->SetTitle("Counts");
    gauss[iPt]->GetHistogram()->GetYaxis()->SetTitle("efafesf");
    GaussPlusPolynom[iPt]->GetHistogram()->GetYaxis()->SetTitle("del_Phi");
    Int_t icolor=0;
    H1D_DetectedV0s_O2_SmallpTinterval[iPt]->SetMarkerStyle(markers[0]);
    H1D_DetectedV0s_O2_SmallpTinterval[iPt]->SetMarkerColor(colors [icolor]);
    H1D_DetectedV0s_O2_SmallpTinterval[iPt]->SetLineColor  (colors [icolor]);
    // bkgparab[iPt]->Draw("same");
    // gauss[iPt]->Draw("same");


    TLatex * textContext = new TLatex (0.18,0.82,NamePart_Latex[ipart]); //BOLD
    textContext->SetTextSize(0.05);
    textContext->SetNDC(kTRUE); //remove if I want x,y in TLatex to be in the coordinate system of the histogram
    textContext->Draw();
    TLatex * textColl = new TLatex (0.18,0.75,"pp #sqrt{#it{s}} = 900 GeV, #it{B_{z}} = 0.2 T");
    textColl->SetTextSize(0.04);
    textColl->SetNDC(kTRUE); //remove if I want x,y in TLatex to be in the coordinate system of the histogram
    textColl->Draw();
    TLatex * text_pTbin = new TLatex (0.18,0.70,Form("pT [%.1f;%.1f] GeV", pTbins[iPt], pTbins[iPt+1]));
    text_pTbin->SetTextSize(0.04);
    text_pTbin->SetNDC(kTRUE); //remove if I want x,y in TLatex to be in the coordinate system of the histogram
    text_pTbin->Draw();
    TLatex * text_fit = new TLatex (0.18,0.63,Form("Gaussian Fit:"));
    text_fit->SetTextSize(0.035);
    text_fit->SetNDC(kTRUE);
    text_fit->Draw();
    TLatex * text_fitParam_mu = new TLatex (0.18,0.59,Form("#it{#mu} = %.0f #pm %.2f MeV/#it{c}^{2}", parGaussParab[iPt][1]*1000., parGaussParabError_mu[iPt]*1000.));
    text_fitParam_mu->SetTextSize(0.035);
    text_fitParam_mu->SetNDC(kTRUE);
    text_fitParam_mu->Draw();
    TLatex * text_fitParam_sigma = new TLatex (0.18,0.55,Form("#it{#sigma} = %.1f #pm %.2f MeV/#it{c}^{2}", parGaussParab[iPt][2]*1000., parGaussParabError_sigma[iPt]*1000.));
    text_fitParam_sigma->SetTextSize(0.035);
    text_fitParam_sigma->SetNDC(kTRUE);
    text_fitParam_sigma->Draw();


    // float xmin_InvMass = H1D_DetectedV0s_O2_SmallpTinterval[iPt]->GetXaxis()->GetXmin();
    // float xmax_InvMass = H1D_DetectedV0s_O2_SmallpTinterval[iPt]->GetXaxis()->GetXmax();
    cout << "xmin_InvMass= " << liminf[ipart] << ", xmax_InvMass= " << limsup[ipart] << endl;

    int ibin_liminf = H1D_DetectedV0s_O2_SmallpTinterval[iPt]->GetXaxis()->FindBin(liminf[ipart]);
    int ibin_limsup = H1D_DetectedV0s_O2_SmallpTinterval[iPt]->GetXaxis()->FindBin(limsup[ipart]);
    int nbinInvMass = ibin_limsup - ibin_liminf;
    cout << "nbinInvMass = " << nbinInvMass << endl;

    int totalCount = H1D_DetectedV0s_O2_SmallpTinterval[iPt]->Integral(ibin_liminf+1,ibin_limsup-1); //integral of a hist counts the number of bins, does not care for bin width: not true integral
    double backgroundCount = bkgparab[iPt]->Integral(liminf[ipart],limsup[ipart])*1./((limsup[ipart]-liminf[ipart])/nbinInvMass); //bkgparab is not an histogram but a function; here we need to convert it to a count
    //issue I actually need the integral of the positive part of the background function; for now kept simple integral as background is close to 0 compared to signal
    // float totalCount = H1D_DetectedV0s_O2_SmallpTinterval[iPt]->Integral(1,-1); //integral of a hist counts the number of bins, does not care for bin width: not true integral
    // float backgroundCount = bkgparab[iPt]->Integral(liminf[ipart],limsup[ipart])*1./((limsup[ipart]-liminf[ipart])/nbinInvMass); //bkgparab is not an histogram but a function; here we need to convert it to a count
    // //issue I actually need the integral of the positive part of the background function; for now kept simple integral as background is close to 0 compared to signal



    double SignalCountFit = totalCount - backgroundCount;
    cout << "totalCount= " << totalCount << ", backgroundCount= " << backgroundCount << endl;
    cout << "SignalCountFit= " << SignalCountFit << endl;

    double dpT = hPtDiff_RawSpectrum_O2->GetXaxis()->GetBinWidth(ibinPt);
    float drapidity = 1.5;//1.5
    double dN_dpT = SignalCountFit *1./dpT*1./drapidity;
    hPtDiff_RawSpectrum_O2->SetBinContent(ibinPt,dN_dpT);
    hPtDiff_RawSpectrum_O2->SetBinError(ibinPt,sqrt(SignalCountFit*1./dpT*1./drapidity)); //Issue with errors; this error of the O2 data spectrum is false, not right formula with fit
    cout << "(ibinPt,dN_dpT) = (" << ibinPt << "," << dN_dpT*1./SelectedEventCount << ")" << endl;

    hPtDiff_RawSpectrum_HEP[ipart]->SetBinError(ibinPt,hPtDiff_RawSpectrum_HEP_statErrors[ipart]->GetBinContent(ibinPt));
  }

  canvasMassvsPt->SaveAs("Fits_invMass_ptDiff_"+NamePart[ipart]+".pdf");

  hPtDiff_RawSpectrum_O2->Scale(1./SelectedEventCount);
  hstat = (TH1D*)hPtDiff_RawSpectrum_O2->Clone("hstat");
  hstat->Reset("M");
  hstat->Divide(hPtDiff_RawSpectrum_O2,hPtDiff_RawSpectrum_HEP[ipart]);

  //////Error Bars///////
  hsyst     = new TH1D("hsyst", "hsyst", 100, 0, 10);
  hsystCorr = new TH1D("hsystCorr", "hsystCorr", 100, 0, 10);
  hstat     ->Sumw2(); 
  hsyst     ->Sumw2();
  hsystCorr ->Sumw2(); 
  Int_t nbinx = hstat->GetNbinsX();
  
  for(Int_t ibinx = 1; ibinx <= nbinx; ibinx++){
    hsyst->SetBinContent (ibinx, hstat->GetBinContent(ibinx));
    hsyst->SetBinError   (ibinx, hstat->GetBinContent(ibinx)*0.08);
    hsystCorr->SetBinContent (ibinx, hstat->GetBinContent(ibinx));
    hsystCorr->SetBinError   (ibinx, hstat->GetBinContent(ibinx)*0.15);
  }

  *SaveAs_Title += "PseudoEfficiency_HEPcomparison_withFit";
  texXtitle = texPtX;
  texYtitle = texptDifferentialYield_HEPratio_pseudoEff;

}


void InvMassDistributions_withFit(TH1D* &hstat, TH1D* &hsyst, TH1D* &hsystCorr, Int_t ipart, TFile* file_O2Analysis, TString* SaveAs_Title, TString* &texXtitle, TString* &texYtitle, Double_t* pTbins, Int_t nbinpT) {

  //HEP file
  TFile* file_HEPAnalysis = new TFile("HEPData-ins881474-v1-root.root");
  TH1D* hPtDiff_RawSpectrum_HEP[numPart]; //= (TH1D*)file_HEPAnalysis->Get("Table 1/Hist1D_y1");
  TH1D* hPtDiff_RawSpectrum_HEP_statErrors[numPart]; // = (TH1D*)file_HEPAnalysis->Get("Table 1/Hist1D_y1_e1");
  hPtDiff_RawSpectrum_HEP[0] = (TH1D*)file_HEPAnalysis->Get("Table 1/Hist1D_y1"); //K0S
  hPtDiff_RawSpectrum_HEP_statErrors[0] = (TH1D*)file_HEPAnalysis->Get("Table 1/Hist1D_y1_e1");
  hPtDiff_RawSpectrum_HEP[1] = (TH1D*)file_HEPAnalysis->Get("Table 2/Hist1D_y1"); //Lambda
  hPtDiff_RawSpectrum_HEP_statErrors[1] = (TH1D*)file_HEPAnalysis->Get("Table 2/Hist1D_y1_e1");
  hPtDiff_RawSpectrum_HEP[2] = (TH1D*)file_HEPAnalysis->Get("Table 3/Hist1D_y1"); //AntiLambda
  hPtDiff_RawSpectrum_HEP_statErrors[2] = (TH1D*)file_HEPAnalysis->Get("Table 3/Hist1D_y1_e1");

  //O2 file
  TH1I* H1I_SelectedEventCount = (TH1I*)file_O2Analysis->Get("lambdakzero-analysis/hSelectedEventCounter");
  int SelectedEventCount = H1I_SelectedEventCount->GetEntries();

  TH3D* H3D_DetectedV0s_O2 = (TH3D*)file_O2Analysis->Get("lambdakzero-analysis/h3dMass"+NamePart[ipart]);

  TH1D* H1D_DetectedV0s_PtAxis_O2 = (TH1D*)H3D_DetectedV0s_O2->ProjectionY(NamehistoInvMass[ipart],0,-1,0,-1);

  //Rebinning InvMass axis if necessary
  TH3D* H3D_DetectedV0s_O2_rebinnedZ = (TH3D*)H3D_DetectedV0s_O2->RebinZ(1,"H3D_DetectedV0s_O2_rebinnedZ");

  //Fit tools initialisation
  TF1 *gauss[nbinpT];
  TF1 *GaussPlusPolynom[nbinpT];
  TF1 *bkgparab[nbinpT];
  
  Float_t fFitResult[nbinpT];
  Double_t parGaussParab[nbinpT][6];  
  Double_t parGaussParabError_mu[nbinpT];  
  Double_t parGaussParabError_sigma[nbinpT];  

  TH1D *H1D_DetectedV0s_O2_SmallpTinterval[nbinpT];
  TH1D *hPtDiff_RawSpectrum_O2 = new TH1D("hPtDiff_RawSpectrum_O2","hPtDiff_RawSpectrum_O2",nbinpT,pTbins);

  TCanvas *canvasMassvsPt = new TCanvas ("canvasMassPt", NamehistoInvMass[ipart], 800, 600);
  // TCanvas *canvasMassvsPt = new TCanvas ("canvasMassPt", NamehistoInvMass[ipart], 1600, 1200);
  canvasMassvsPt->Divide(std::ceil(std::sqrt(nbinpT)),std::ceil(std::sqrt(nbinpT)));
  TH1 *hFrame[nbinpT];

  for(int ibinPt = 1; ibinPt <= nbinpT; ibinPt++){
    int iPt = ibinPt-1; //different convention for bins that start at 1 and arrays that start at 0

    canvasMassvsPt->cd(ibinPt); //ATTENTION: +1 because canvas divide starts counting from 1

    int ibinPt_originalAxis_low = H3D_DetectedV0s_O2->GetYaxis()->FindBin(pTbins[iPt]);
    int ibinPt_originalAxis_up = H3D_DetectedV0s_O2->GetYaxis()->FindBin(pTbins[iPt+1]) - 1;
    cout << "ibinPt: " << ibinPt << "; originalAxis [low;up] = [" << ibinPt_originalAxis_low << ";" << ibinPt_originalAxis_up << "]" << endl;

    H1D_DetectedV0s_O2_SmallpTinterval[iPt] = (TH1D*)H3D_DetectedV0s_O2_rebinnedZ->ProjectionZ("InvMass_"+NamePart[ipart]+Form("_Pt%i", ibinPt),0,-1,ibinPt_originalAxis_low,ibinPt_originalAxis_up);
    hFrame[ibinPt] = canvasMassvsPt->DrawFrame(LowerLimitDisplay_Mass[ipart],0,UpperLimitDisplay_Mass[ipart],1.45*H1D_DetectedV0s_O2_SmallpTinterval[iPt]->GetMaximum());
    hFrame[ibinPt]->SetXTitle(texInvMassDecays_titles[ipart]->Data());
    hFrame[ibinPt]->SetYTitle("Counts");
    H1D_DetectedV0s_O2_SmallpTinterval[iPt]->Draw("E,same");

    gauss[iPt] = new TF1("gauss_"+NamePart[ipart]+Form("_Pt%i", ibinPt), "gaus", min_range_signal[ipart], max_range_signal[ipart]);
    gauss[iPt]->SetLineColor(kBlue);   
    gauss[iPt]->SetParName(0, "norm");
    gauss[iPt]->SetParName(1, "mean");
    gauss[iPt]->SetParName(2, "sigma");
    gauss[iPt]->SetParameters(1., MassPart[ipart], 0.0002);
    gauss[iPt]->SetParLimits(0, 0., 1.1*H1D_DetectedV0s_O2_SmallpTinterval[iPt]->GetBinContent(H1D_DetectedV0s_O2_SmallpTinterval[iPt]->GetMaximumBin()));
    gauss[iPt]->SetParLimits(1, min_range_signal[ipart], max_range_signal[ipart]);
    gauss[iPt]->SetParLimits(2, 0.0001, 0.05);

    bkgparab[iPt] = new TF1("parab_"+NamePart[ipart]+Form("_Pt%i", ibinPt), fparab, liminf[ipart], limsup[ipart], 4);
    bkgparab[iPt]->SetLineColor(kGreen);
    bkgparab[iPt]->FixParameter(3, ipart);

    GaussPlusPolynom[iPt] = new TF1("GaussPlusPolynom_"+NamePart[ipart]+Form("_Pt%i", ibinPt), "gaus(0)+pol2(3)", liminf[ipart], limsup[ipart]);
    GaussPlusPolynom[iPt]->SetLineColor(kRed); 
    GaussPlusPolynom[iPt]->SetParName(0, "norm");
    GaussPlusPolynom[iPt]->SetParName(1, "mean");
    GaussPlusPolynom[iPt]->SetParName(2, "sigma");
    GaussPlusPolynom[iPt]->SetParName(3, "p0");
    GaussPlusPolynom[iPt]->SetParName(4, "p1");
    GaussPlusPolynom[iPt]->SetParName(5, "p2");

    GaussPlusPolynom[iPt]->SetParameters(1., MassPart[ipart], 0.0002);
    GaussPlusPolynom[iPt]->SetParLimits(0, 0., 1.1*H1D_DetectedV0s_O2_SmallpTinterval[iPt]->GetBinContent(H1D_DetectedV0s_O2_SmallpTinterval[iPt]->GetMaximumBin()));
    GaussPlusPolynom[iPt]->SetParLimits(1, min_range_signal[ipart], max_range_signal[ipart]);
    GaussPlusPolynom[iPt]->SetParLimits(2, 0.0001, 0.05);

    H1D_DetectedV0s_O2_SmallpTinterval[iPt]->Fit(gauss[iPt], "R0QL");
    reject = true; // for background fit: ignore signal range when fitting
    H1D_DetectedV0s_O2_SmallpTinterval[iPt]->Fit(bkgparab[iPt], "R0QL");
    reject = false; // to allow plotting it without the gap

    gauss[iPt]->GetParameters(&parGaussParab[iPt][0]);
    bkgparab[iPt]->GetParameters(&parGaussParab[iPt][3]);
    GaussPlusPolynom[iPt]->SetParameters(parGaussParab[iPt]);

    fFitResult[iPt] = H1D_DetectedV0s_O2_SmallpTinterval[iPt]->Fit(GaussPlusPolynom[iPt],"SRB+QL"); //L option is log likelyhood: would supposedly be better in case of counts

    //give background the final parameters from that last signal+background fit
    GaussPlusPolynom[iPt]->GetParameters(&parGaussParab[iPt][0]);
    bkgparab[iPt]->SetParameters(&parGaussParab[iPt][3]);

    parGaussParabError_mu[iPt] = GaussPlusPolynom[iPt]->GetParError(1);
    parGaussParabError_sigma[iPt] = GaussPlusPolynom[iPt]->GetParError(2);

    //plot Inv Mass histogram and the background and signal fits for each iPt bin
    // gPad->DrawFrame(LowerLimitDisplay_Mass[ipart],0,UpperLimitDisplay_Mass[ipart],1.45*H1D_DetectedV0s_O2_SmallpTinterval[iPt]->GetMaximum());
    // gPad->SetXTitle(texInvMassDecays_titles[ipart]->Data());
    // gPad->SetYTitle("Counts");
    // H1D_DetectedV0s_O2_SmallpTinterval[iPt]->Draw("E,same");
    H1D_DetectedV0s_O2_SmallpTinterval[iPt]->GetXaxis()->SetTitle(texInvMassDecays_titles[ipart]->Data());
    H1D_DetectedV0s_O2_SmallpTinterval[iPt]->GetYaxis()->SetTitle("Counts");
    gauss[iPt]->GetHistogram()->GetYaxis()->SetTitle("efafesf");
    GaussPlusPolynom[iPt]->GetHistogram()->GetYaxis()->SetTitle("del_Phi");
    Int_t icolor=0;
    H1D_DetectedV0s_O2_SmallpTinterval[iPt]->SetMarkerStyle(markers[0]);
    H1D_DetectedV0s_O2_SmallpTinterval[iPt]->SetMarkerColor(colors [icolor]);
    H1D_DetectedV0s_O2_SmallpTinterval[iPt]->SetLineColor  (colors [icolor]);
    // bkgparab[iPt]->Draw("same");
    // gauss[iPt]->Draw("same");


    // TLatex * text_pTbin = new TLatex (0.62,0.75,Form("pT [%.1f;%.1f] GeV", pTbins[iPt], pTbins[iPt+1]));
    // text_pTbin->SetNDC(kTRUE);
    // text_pTbin->Draw();
    TLatex * textContext = new TLatex (0.18,0.82,"#bf{ALICE Performance}"); //BOLD
    textContext->SetTextSize(0.05);
    textContext->SetNDC(kTRUE); //remove if I want x,y in TLatex to be in the coordinate system of the histogram
    textContext->Draw();
    TLatex * textColl = new TLatex (0.18,0.75,"pp #sqrt{#it{s}} = 0.9 TeV, pilot beam 2021");
    textColl->SetTextSize(0.04);
    textColl->SetNDC(kTRUE); //remove if I want x,y in TLatex to be in the coordinate system of the histogram
    textColl->Draw();
    TLatex * text_part = new TLatex (0.18,0.70,NamePart_Latex[ipart]+",   0 < #it{p}_{T} < 10 GeV/#it{c}");
    text_part->SetTextSize(0.04);
    text_part->SetNDC(kTRUE); //remove if I want x,y in TLatex to be in the coordinate system of the histogram
    text_part->Draw();
    TLatex * text_fit = new TLatex (0.18,0.63,Form("Gaussian Fit:"));
    text_fit->SetTextSize(0.035);
    text_fit->SetNDC(kTRUE);
    text_fit->Draw();
    // TLatex * text_fitParam_mu = new TLatex (0.18,0.59,Form("#it{#mu} = %.0f #pm %.2f MeV/#it{c}^{2}", parGaussParab[iPt][1]*1000., parGaussParabError_mu[iPt]*1000.));
    TLatex * text_fitParam_mu = new TLatex (0.18,0.59,Form("#it{#mu} = %.0f MeV/#it{c}^{2}", parGaussParab[iPt][1]*1000.));
    text_fitParam_mu->SetTextSize(0.035);
    text_fitParam_mu->SetNDC(kTRUE);
    text_fitParam_mu->Draw();
    // TLatex * text_fitParam_sigma = new TLatex (0.18,0.55,Form("#it{#sigma} = %.1f #pm %.2f MeV/#it{c}^{2}", parGaussParab[iPt][2]*1000., parGaussParabError_sigma[iPt]*1000.));
    TLatex * text_fitParam_sigma = new TLatex (0.18,0.55,Form("#it{#sigma} = %.1f MeV/#it{c}^{2}", parGaussParab[iPt][2]*1000.));
    text_fitParam_sigma->SetTextSize(0.035);
    text_fitParam_sigma->SetNDC(kTRUE);
    text_fitParam_sigma->Draw();

    // float xmin_InvMass = H1D_DetectedV0s_O2_SmallpTinterval[iPt]->GetXaxis()->GetXmin();
    // float xmax_InvMass = H1D_DetectedV0s_O2_SmallpTinterval[iPt]->GetXaxis()->GetXmax();
    cout << "xmin_InvMass= " << liminf[ipart] << ", xmax_InvMass= " << limsup[ipart] << endl;

    int ibin_liminf = H1D_DetectedV0s_O2_SmallpTinterval[iPt]->GetXaxis()->FindBin(liminf[ipart]);
    int ibin_limsup = H1D_DetectedV0s_O2_SmallpTinterval[iPt]->GetXaxis()->FindBin(limsup[ipart]);
    int nbinInvMass = ibin_limsup - ibin_liminf;
    cout << "nbinInvMass = " << nbinInvMass << endl;

    int totalCount = H1D_DetectedV0s_O2_SmallpTinterval[iPt]->Integral(ibin_liminf+1,ibin_limsup-1); //integral of a hist counts the number of bins, does not care for bin width: not true integral
    double backgroundCount = bkgparab[iPt]->Integral(liminf[ipart],limsup[ipart])*1./((limsup[ipart]-liminf[ipart])/nbinInvMass); //bkgparab is not an histogram but a function; here we need to convert it to a count
    //issue I actually need the integral of the positive part of the background function; for now kept simple integral as background is close to 0 compared to signal
    // float totalCount = H1D_DetectedV0s_O2_SmallpTinterval[iPt]->Integral(1,-1); //integral of a hist counts the number of bins, does not care for bin width: not true integral
    // float backgroundCount = bkgparab[iPt]->Integral(liminf[ipart],limsup[ipart])*1./((limsup[ipart]-liminf[ipart])/nbinInvMass); //bkgparab is not an histogram but a function; here we need to convert it to a count
    // //issue I actually need the integral of the positive part of the background function; for now kept simple integral as background is close to 0 compared to signal



    double SignalCountFit = totalCount - backgroundCount;
    cout << "totalCount= " << totalCount << ", backgroundCount= " << backgroundCount << endl;
    cout << "SignalCountFit= " << SignalCountFit << endl;

    double dpT = hPtDiff_RawSpectrum_O2->GetXaxis()->GetBinWidth(ibinPt);
    float drapidity = 1.5;//1.5
    double dN_dpT = SignalCountFit *1./dpT*1./drapidity;
    hPtDiff_RawSpectrum_O2->SetBinContent(ibinPt,dN_dpT);
    hPtDiff_RawSpectrum_O2->SetBinError(ibinPt,sqrt(SignalCountFit*1./dpT*1./drapidity)); //Issue with errors; this error of the O2 data spectrum is false, not right formula with fit
    cout << "(ibinPt,dN_dpT) = (" << ibinPt << "," << dN_dpT*1./SelectedEventCount << ")" << endl;

    hPtDiff_RawSpectrum_HEP[ipart]->SetBinError(ibinPt,hPtDiff_RawSpectrum_HEP_statErrors[ipart]->GetBinContent(ibinPt));
  }

  canvasMassvsPt->SaveAs("InvMass_ptDiff_"+NamePart[ipart]+".pdf");

  hPtDiff_RawSpectrum_O2->Scale(1./SelectedEventCount);
  hstat = (TH1D*)hPtDiff_RawSpectrum_O2->Clone("hstat");
  hstat->Reset("M");
  // hstat->Divide(hPtDiff_RawSpectrum_O2,hPtDiff_RawSpectrum_HEP[ipart]);

  //////Error Bars///////
  hsyst     = new TH1D("hsyst", "hsyst", 100, 0, 10);
  hsystCorr = new TH1D("hsystCorr", "hsystCorr", 100, 0, 10);
  hstat     ->Sumw2(); 
  hsyst     ->Sumw2();
  hsystCorr ->Sumw2(); 
  Int_t nbinx = hstat->GetNbinsX();
  
  for(Int_t ibinx = 1; ibinx <= nbinx; ibinx++){
    hsyst->SetBinContent (ibinx, hstat->GetBinContent(ibinx));
    hsyst->SetBinError   (ibinx, hstat->GetBinContent(ibinx)*0.08);
    hsystCorr->SetBinContent (ibinx, hstat->GetBinContent(ibinx));
    hsystCorr->SetBinError   (ibinx, hstat->GetBinContent(ibinx)*0.15);
  }

  *SaveAs_Title += "ToBeDeleted";
  texXtitle = texPtX;
  texYtitle = texptDifferentialYield_HEPratio_pseudoEff;

}

void InvMassDistributions_withFit_MC_datalike(TH1D* &hstat, TH1D* &hsyst, TH1D* &hsystCorr, Int_t ipart, TFile* file_O2Analysis, TString* SaveAs_Title, TString* &texXtitle, TString* &texYtitle, Double_t* pTbins, Int_t nbinpT) {

  // // HEP file
  // TFile* file_HEPAnalysis = new TFile("HEPData-ins881474-v1-root.root");
  // TH1D* hPtDiff_RawSpectrum_HEP[numPart]; //= (TH1D*)file_HEPAnalysis->Get("Table 1/Hist1D_y1");
  // TH1D* hPtDiff_RawSpectrum_HEP_statErrors[numPart]; // = (TH1D*)file_HEPAnalysis->Get("Table 1/Hist1D_y1_e1");
  // hPtDiff_RawSpectrum_HEP[0] = (TH1D*)file_HEPAnalysis->Get("Table 1/Hist1D_y1"); //K0S
  // hPtDiff_RawSpectrum_HEP_statErrors[0] = (TH1D*)file_HEPAnalysis->Get("Table 1/Hist1D_y1_e1");
  // hPtDiff_RawSpectrum_HEP[1] = (TH1D*)file_HEPAnalysis->Get("Table 2/Hist1D_y1"); //Lambda
  // hPtDiff_RawSpectrum_HEP_statErrors[1] = (TH1D*)file_HEPAnalysis->Get("Table 2/Hist1D_y1_e1");
  // hPtDiff_RawSpectrum_HEP[2] = (TH1D*)file_HEPAnalysis->Get("Table 3/Hist1D_y1"); //AntiLambda
  // hPtDiff_RawSpectrum_HEP_statErrors[2] = (TH1D*)file_HEPAnalysis->Get("Table 3/Hist1D_y1_e1");

  //O2 file
  TH1I* H1I_SelectedEventCount = (TH1I*)file_O2Analysis->Get("lambdakzero-analysis-mc/hSelectedEventCounter");
  int SelectedEventCount = H1I_SelectedEventCount->GetEntries();

  TH3D* H3D_DetectedV0s_O2 = (TH3D*)file_O2Analysis->Get("lambdakzero-analysis-mc/h3dMass"+NamePart[ipart]);

  TH1D* H1D_DetectedV0s_PtAxis_O2 = (TH1D*)H3D_DetectedV0s_O2->ProjectionY(NamehistoInvMass[ipart],0,-1,0,-1);

  //Rebinning InvMass axis if necessary
  TH3D* H3D_DetectedV0s_O2_rebinnedZ = (TH3D*)H3D_DetectedV0s_O2->RebinZ(1,"H3D_DetectedV0s_O2_rebinnedZ");

  //Fit tools initialisation
  TF1 *gauss[nbinpT];
  TF1 *GaussPlusPolynom[nbinpT];
  TF1 *bkgparab[nbinpT];
  
  Float_t fFitResult[nbinpT];
  // Double_t parGaussParab[nbinpT][6]; //parabola background
  Double_t parGaussParab[nbinpT][5];   //line background
  Double_t parGaussParabError_mu[nbinpT];  
  Double_t parGaussParabError_sigma[nbinpT];  

  TH1D *H1D_DetectedV0s_O2_SmallpTinterval[nbinpT];
  TH1D *hPtDiff_RawSpectrum_O2 = new TH1D("hPtDiff_RawSpectrum_O2","hPtDiff_RawSpectrum_O2",nbinpT,pTbins);

  TCanvas *canvasMassvsPt = new TCanvas ("canvasMassPt", NamehistoInvMass[ipart], 800, 600);
  // TCanvas *canvasMassvsPt = new TCanvas ("canvasMassPt", NamehistoInvMass[ipart], 1600, 1200);
  canvasMassvsPt->Divide(std::ceil(std::sqrt(nbinpT)),std::ceil(std::sqrt(nbinpT)));
  TH1 *hFrame[nbinpT];

  for(int ibinPt = 1; ibinPt <= nbinpT; ibinPt++){
    int iPt = ibinPt-1; //different convention for bins that start at 1 and arrays that start at 0

    canvasMassvsPt->cd(ibinPt); //ATTENTION: +1 because canvas divide starts counting from 1

    int ibinPt_originalAxis_low = H3D_DetectedV0s_O2->GetYaxis()->FindBin(pTbins[iPt]);
    int ibinPt_originalAxis_up = H3D_DetectedV0s_O2->GetYaxis()->FindBin(pTbins[iPt+1]) - 1;
    cout << "ibinPt: " << ibinPt << "; originalAxis [low;up] = [" << ibinPt_originalAxis_low << ";" << ibinPt_originalAxis_up << "]" << endl;

    H1D_DetectedV0s_O2_SmallpTinterval[iPt] = (TH1D*)H3D_DetectedV0s_O2_rebinnedZ->ProjectionZ("InvMass_"+NamePart[ipart]+Form("_Pt%i", ibinPt),0,-1,ibinPt_originalAxis_low,ibinPt_originalAxis_up);
    hFrame[ibinPt] = canvasMassvsPt->DrawFrame(LowerLimitDisplay_Mass[ipart],0,UpperLimitDisplay_Mass[ipart],1.45*H1D_DetectedV0s_O2_SmallpTinterval[iPt]->GetMaximum());
    hFrame[ibinPt]->SetXTitle(texInvMassDecays_titles[ipart]->Data());
    hFrame[ibinPt]->SetYTitle("Counts");
    H1D_DetectedV0s_O2_SmallpTinterval[iPt]->Draw("E,same");

    gauss[iPt] = new TF1("gauss_"+NamePart[ipart]+Form("_Pt%i", ibinPt), "gaus", min_range_signal[ipart], max_range_signal[ipart]);
    gauss[iPt]->SetLineColor(kBlue);   
    gauss[iPt]->SetParName(0, "norm");
    gauss[iPt]->SetParName(1, "mean");
    gauss[iPt]->SetParName(2, "sigma");
    gauss[iPt]->SetParameters(1., MassPart[ipart], 0.0002);
    gauss[iPt]->SetParLimits(0, 0., 1.1*H1D_DetectedV0s_O2_SmallpTinterval[iPt]->GetBinContent(H1D_DetectedV0s_O2_SmallpTinterval[iPt]->GetMaximumBin()));
    gauss[iPt]->SetParLimits(1, min_range_signal[ipart], max_range_signal[ipart]);
    gauss[iPt]->SetParLimits(2, 0.0001, 0.05);

    // //parabola background
    // bkgparab[iPt] = new TF1("parab_"+NamePart[ipart]+Form("_Pt%i", ibinPt), fparab, liminf[ipart], limsup[ipart], 4);
    // bkgparab[iPt]->SetLineColor(kGreen);
    // bkgparab[iPt]->FixParameter(3, ipart);

    // line background
    bkgparab[iPt] = new TF1("line_"+NamePart[ipart]+Form("_Pt%i", ibinPt), fline, liminf[ipart], limsup[ipart], 3);
    bkgparab[iPt]->SetLineColor(kGreen);
    bkgparab[iPt]->FixParameter(2, ipart);

    // GaussPlusPolynom[iPt] = new TF1("GaussPlusPolynom_"+NamePart[ipart]+Form("_Pt%i", ibinPt), "gaus(0)+pol2(3)", liminf[ipart], limsup[ipart]);
    GaussPlusPolynom[iPt] = new TF1("GaussPlusPolynom_"+NamePart[ipart]+Form("_Pt%i", ibinPt), "gaus(0)+pol1(3)", liminf[ipart], limsup[ipart]);
    GaussPlusPolynom[iPt]->SetLineColor(kRed); 
    GaussPlusPolynom[iPt]->SetParName(0, "norm");
    GaussPlusPolynom[iPt]->SetParName(1, "mean");
    GaussPlusPolynom[iPt]->SetParName(2, "sigma");
    GaussPlusPolynom[iPt]->SetParName(3, "p0");
    GaussPlusPolynom[iPt]->SetParName(4, "p1");
    // GaussPlusPolynom[iPt]->SetParName(5, "p2"); //parabola background

    GaussPlusPolynom[iPt]->SetParameters(1., MassPart[ipart], 0.0002);
    GaussPlusPolynom[iPt]->SetParLimits(0, 0., 1.1*H1D_DetectedV0s_O2_SmallpTinterval[iPt]->GetBinContent(H1D_DetectedV0s_O2_SmallpTinterval[iPt]->GetMaximumBin()));
    GaussPlusPolynom[iPt]->SetParLimits(1, min_range_signal[ipart], max_range_signal[ipart]);
    GaussPlusPolynom[iPt]->SetParLimits(2, 0.0001, 0.05);

    H1D_DetectedV0s_O2_SmallpTinterval[iPt]->Fit(gauss[iPt], "R0QL");
    reject = true; // for background fit: ignore signal range when fitting
    H1D_DetectedV0s_O2_SmallpTinterval[iPt]->Fit(bkgparab[iPt], "R0QL");
    reject = false; // to allow plotting it without the gap

    gauss[iPt]->GetParameters(&parGaussParab[iPt][0]);
    bkgparab[iPt]->GetParameters(&parGaussParab[iPt][3]);
    GaussPlusPolynom[iPt]->SetParameters(parGaussParab[iPt]);

    fFitResult[iPt] = H1D_DetectedV0s_O2_SmallpTinterval[iPt]->Fit(GaussPlusPolynom[iPt],"SRB+QL"); //L option is log likelyhood: would supposedly be better in case of counts

    //give background the final parameters from that last signal+background fit
    GaussPlusPolynom[iPt]->GetParameters(&parGaussParab[iPt][0]);
    bkgparab[iPt]->SetParameters(&parGaussParab[iPt][3]);

    parGaussParabError_mu[iPt] = GaussPlusPolynom[iPt]->GetParError(1);
    parGaussParabError_sigma[iPt] = GaussPlusPolynom[iPt]->GetParError(2);

    //plot Inv Mass histogram and the background and signal fits for each iPt bin
    // gPad->DrawFrame(LowerLimitDisplay_Mass[ipart],0,UpperLimitDisplay_Mass[ipart],1.45*H1D_DetectedV0s_O2_SmallpTinterval[iPt]->GetMaximum());
    // gPad->SetXTitle(texInvMassDecays_titles[ipart]->Data());
    // gPad->SetYTitle("Counts");
    // H1D_DetectedV0s_O2_SmallpTinterval[iPt]->Draw("E,same");
    H1D_DetectedV0s_O2_SmallpTinterval[iPt]->GetXaxis()->SetTitle(texInvMassDecays_titles[ipart]->Data());
    H1D_DetectedV0s_O2_SmallpTinterval[iPt]->GetYaxis()->SetTitle("Counts");
    gauss[iPt]->GetHistogram()->GetYaxis()->SetTitle("efafesf");
    GaussPlusPolynom[iPt]->GetHistogram()->GetYaxis()->SetTitle("del_Phi");
    Int_t icolor=0;
    H1D_DetectedV0s_O2_SmallpTinterval[iPt]->SetMarkerStyle(markers[0]);
    H1D_DetectedV0s_O2_SmallpTinterval[iPt]->SetMarkerColor(colors [icolor]);
    H1D_DetectedV0s_O2_SmallpTinterval[iPt]->SetLineColor  (colors [icolor]);
    // bkgparab[iPt]->Draw("same");
    // gauss[iPt]->Draw("same");


    // TLatex * text_pTbin = new TLatex (0.62,0.75,Form("pT [%.1f;%.1f] GeV", pTbins[iPt], pTbins[iPt+1]));
    // text_pTbin->SetNDC(kTRUE);
    // text_pTbin->Draw();
    TLatex * textContext = new TLatex (0.18,0.82,"#bf{ALICE Performance}"); //BOLD
    textContext->SetTextSize(0.05);
    textContext->SetNDC(kTRUE); //remove if I want x,y in TLatex to be in the coordinate system of the histogram
    textContext->Draw();
    TLatex * textColl = new TLatex (0.18,0.75,"pp #sqrt{#it{s}} = 0.9 TeV, pilot beam 2021");
    textColl->SetTextSize(0.04);
    textColl->SetNDC(kTRUE); //remove if I want x,y in TLatex to be in the coordinate system of the histogram
    textColl->Draw();
    TLatex * text_part = new TLatex (0.18,0.70,NamePart_Latex[ipart]+",   0 < #it{p}_{T} < 10 GeV/#it{c}");
    text_part->SetTextSize(0.04);
    text_part->SetNDC(kTRUE); //remove if I want x,y in TLatex to be in the coordinate system of the histogram
    text_part->Draw();
    TLatex * text_fit = new TLatex (0.18,0.63,Form("Gaussian Fit:"));
    text_fit->SetTextSize(0.035);
    text_fit->SetNDC(kTRUE);
    text_fit->Draw();
    // TLatex * text_fitParam_mu = new TLatex (0.18,0.59,Form("#it{#mu} = %.0f #pm %.2f MeV/#it{c}^{2}", parGaussParab[iPt][1]*1000., parGaussParabError_mu[iPt]*1000.));
    TLatex * text_fitParam_mu = new TLatex (0.18,0.59,Form("#it{#mu} = %.0f MeV/#it{c}^{2}", parGaussParab[iPt][1]*1000.));
    text_fitParam_mu->SetTextSize(0.035);
    text_fitParam_mu->SetNDC(kTRUE);
    text_fitParam_mu->Draw();
    // TLatex * text_fitParam_sigma = new TLatex (0.18,0.55,Form("#it{#sigma} = %.1f #pm %.2f MeV/#it{c}^{2}", parGaussParab[iPt][2]*1000., parGaussParabError_sigma[iPt]*1000.));
    TLatex * text_fitParam_sigma = new TLatex (0.18,0.55,Form("#it{#sigma} = %.1f MeV/#it{c}^{2}", parGaussParab[iPt][2]*1000.));
    text_fitParam_sigma->SetTextSize(0.035);
    text_fitParam_sigma->SetNDC(kTRUE);
    text_fitParam_sigma->Draw();

    // float xmin_InvMass = H1D_DetectedV0s_O2_SmallpTinterval[iPt]->GetXaxis()->GetXmin();
    // float xmax_InvMass = H1D_DetectedV0s_O2_SmallpTinterval[iPt]->GetXaxis()->GetXmax();
    cout << "xmin_InvMass= " << liminf[ipart] << ", xmax_InvMass= " << limsup[ipart] << endl;

    int ibin_liminf = H1D_DetectedV0s_O2_SmallpTinterval[iPt]->GetXaxis()->FindBin(liminf[ipart]);
    int ibin_limsup = H1D_DetectedV0s_O2_SmallpTinterval[iPt]->GetXaxis()->FindBin(limsup[ipart]);
    int nbinInvMass = ibin_limsup - ibin_liminf;
    cout << "nbinInvMass = " << nbinInvMass << endl;

    int totalCount = H1D_DetectedV0s_O2_SmallpTinterval[iPt]->Integral(ibin_liminf+1,ibin_limsup-1); //integral of a hist counts the number of bins, does not care for bin width: not true integral
    double backgroundCount = bkgparab[iPt]->Integral(liminf[ipart],limsup[ipart])*1./((limsup[ipart]-liminf[ipart])/nbinInvMass); //bkgparab is not an histogram but a function; here we need to convert it to a count
    //issue I actually need the integral of the positive part of the background function; for now kept simple integral as background is close to 0 compared to signal
    // float totalCount = H1D_DetectedV0s_O2_SmallpTinterval[iPt]->Integral(1,-1); //integral of a hist counts the number of bins, does not care for bin width: not true integral
    // float backgroundCount = bkgparab[iPt]->Integral(liminf[ipart],limsup[ipart])*1./((limsup[ipart]-liminf[ipart])/nbinInvMass); //bkgparab is not an histogram but a function; here we need to convert it to a count
    // //issue I actually need the integral of the positive part of the background function; for now kept simple integral as background is close to 0 compared to signal



    double SignalCountFit = totalCount - backgroundCount;
    cout << "totalCount= " << totalCount << ", backgroundCount= " << backgroundCount << endl;
    cout << "SignalCountFit= " << SignalCountFit << endl;

    double dpT = hPtDiff_RawSpectrum_O2->GetXaxis()->GetBinWidth(ibinPt);
    float drapidity = 1.5;//1.5
    double dN_dpT = SignalCountFit *1./dpT*1./drapidity;
    hPtDiff_RawSpectrum_O2->SetBinContent(ibinPt,dN_dpT);
    hPtDiff_RawSpectrum_O2->SetBinError(ibinPt,sqrt(SignalCountFit*1./dpT*1./drapidity)); //Issue with errors; this error of the O2 data spectrum is false, not right formula with fit
    cout << "(ibinPt,dN_dpT) = (" << ibinPt << "," << dN_dpT*1./SelectedEventCount << ")" << endl;

    // hPtDiff_RawSpectrum_HEP[ipart]->SetBinError(ibinPt,hPtDiff_RawSpectrum_HEP_statErrors[ipart]->GetBinContent(ibinPt));
  }

  canvasMassvsPt->SaveAs("InvMass_ptDiff_"+NamePart[ipart]+".pdf");

  hPtDiff_RawSpectrum_O2->Scale(1./SelectedEventCount);
  hstat = (TH1D*)hPtDiff_RawSpectrum_O2->Clone("hstat");
  hstat->Reset("M");
  // hstat->Divide(hPtDiff_RawSpectrum_O2,hPtDiff_RawSpectrum_HEP[ipart]);

  //////Error Bars///////
  hsyst     = new TH1D("hsyst", "hsyst", 100, 0, 10);
  hsystCorr = new TH1D("hsystCorr", "hsystCorr", 100, 0, 10);
  hstat     ->Sumw2(); 
  hsyst     ->Sumw2();
  hsystCorr ->Sumw2(); 
  Int_t nbinx = hstat->GetNbinsX();
  
  for(Int_t ibinx = 1; ibinx <= nbinx; ibinx++){
    hsyst->SetBinContent (ibinx, hstat->GetBinContent(ibinx));
    hsyst->SetBinError   (ibinx, hstat->GetBinContent(ibinx)*0.08);
    hsystCorr->SetBinContent (ibinx, hstat->GetBinContent(ibinx));
    hsystCorr->SetBinError   (ibinx, hstat->GetBinContent(ibinx)*0.15);
  }

  *SaveAs_Title += "ToBeDeleted";
  texXtitle = texPtX;
  texYtitle = texptDifferentialYield_HEPratio_pseudoEff;

}

void PtDifferential_Fit(TH1D* &hPtDiff_RawSpectrum_O2, TH3D* H3D_DetectedV0s_O2, TH3D* H3D_DetectedV0s_O2_rebinnedZ, Int_t ipart, Double_t* pTbins, Int_t nbinpT, Double_t* parGaussianParab_Sigma, Double_t* parGaussianParab_Sigma_error) {

  //Fit tools initialisation
  TF1 *gauss[nbinpT];
  TF1 *GaussPlusPolynom[nbinpT];
  TF1 *bkgparab[nbinpT];
  
  Float_t fFitResult[nbinpT];
  Double_t parGaussParab[nbinpT][6];  

  TH1D *H1D_DetectedV0s_O2_SmallpTinterval[nbinpT];
  // TH1D *hPtDiff_RawSpectrum_O2 = new TH1D("hPtDiff_RawSpectrum_O2","hPtDiff_RawSpectrum_O2",nbinpT,pTbins);

  // TCanvas *canvasMassvsPt = new TCanvas ("canvasMassPt", NamehistoInvMass[ipart], 800, 600);
  TCanvas *canvasMassvsPt = new TCanvas ("canvasMassPt", NamehistoInvMass[ipart], 1600, 1200);
  canvasMassvsPt->Divide(4,4);

  for(int ibinPt = 1; ibinPt <= nbinpT; ibinPt++){
    int iPt = ibinPt-1; //different convention for bins that start at 1 and arrays that start at 0

    canvasMassvsPt->cd(ibinPt); //ATTENTION: +1 because canvas divide starts counting from 1

    int ibinPt_originalAxis_low = H3D_DetectedV0s_O2->GetYaxis()->FindBin(pTbins[iPt]);
    int ibinPt_originalAxis_up = H3D_DetectedV0s_O2->GetYaxis()->FindBin(pTbins[iPt+1]) - 1;
    cout << "ibinPt: " << ibinPt << "; originalAxis [low;up] = [" << ibinPt_originalAxis_low << ";" << ibinPt_originalAxis_up << "]" << endl;

    H1D_DetectedV0s_O2_SmallpTinterval[iPt] = (TH1D*)H3D_DetectedV0s_O2_rebinnedZ->ProjectionZ("InvMass_"+NamePart[ipart]+Form("_Pt%i", ibinPt),0,-1,ibinPt_originalAxis_low,ibinPt_originalAxis_up);

    gauss[iPt] = new TF1("gauss_"+NamePart[ipart]+Form("_Pt%i", ibinPt), "gaus", min_range_signal[ipart], max_range_signal[ipart]);
    gauss[iPt]->SetLineColor(kBlue);   
    gauss[iPt]->SetParName(0, "norm");
    gauss[iPt]->SetParName(1, "mean");
    gauss[iPt]->SetParName(2, "sigma");
    gauss[iPt]->SetParameters(1., MassPart[ipart], 0.0002);
    gauss[iPt]->SetParLimits(0, 0., 1.1*H1D_DetectedV0s_O2_SmallpTinterval[iPt]->GetBinContent(H1D_DetectedV0s_O2_SmallpTinterval[iPt]->GetMaximumBin()));
    gauss[iPt]->SetParLimits(1, min_range_signal[ipart], max_range_signal[ipart]);
    gauss[iPt]->SetParLimits(2, 0.0001, 0.05);

    bkgparab[iPt] = new TF1("parab_"+NamePart[ipart]+Form("_Pt%i", ibinPt), fparab, liminf[ipart], limsup[ipart], 4);
    bkgparab[iPt]->SetLineColor(kGreen);
    bkgparab[iPt]->FixParameter(3, ipart);

    GaussPlusPolynom[iPt] = new TF1("GaussPlusPolynom_"+NamePart[ipart]+Form("_Pt%i", ibinPt), "gaus(0)+pol2(3)", liminf[ipart], limsup[ipart]);
    GaussPlusPolynom[iPt]->SetLineColor(kRed); 
    GaussPlusPolynom[iPt]->SetParName(0, "norm");
    GaussPlusPolynom[iPt]->SetParName(1, "mean");
    GaussPlusPolynom[iPt]->SetParName(2, "sigma");
    GaussPlusPolynom[iPt]->SetParName(3, "p0");
    GaussPlusPolynom[iPt]->SetParName(4, "p1");
    GaussPlusPolynom[iPt]->SetParName(5, "p2");

    GaussPlusPolynom[iPt]->SetParameters(1., MassPart[ipart], 0.0002);
    GaussPlusPolynom[iPt]->SetParLimits(0, 0., 1.1*H1D_DetectedV0s_O2_SmallpTinterval[iPt]->GetBinContent(H1D_DetectedV0s_O2_SmallpTinterval[iPt]->GetMaximumBin()));
    GaussPlusPolynom[iPt]->SetParLimits(1, min_range_signal[ipart], max_range_signal[ipart]);
    GaussPlusPolynom[iPt]->SetParLimits(2, 0.0001, 0.05);

    H1D_DetectedV0s_O2_SmallpTinterval[iPt]->Fit(gauss[iPt], "R0QL");
    reject = true; // for background fit: ignore signal range when fitting
    H1D_DetectedV0s_O2_SmallpTinterval[iPt]->Fit(bkgparab[iPt], "R0QL");
    reject = false; // to allow plotting it without the gap

    gauss[iPt]->GetParameters(&parGaussParab[iPt][0]);
    bkgparab[iPt]->GetParameters(&parGaussParab[iPt][3]);
    GaussPlusPolynom[iPt]->SetParameters(parGaussParab[iPt]);

    fFitResult[iPt] = H1D_DetectedV0s_O2_SmallpTinterval[iPt]->Fit(GaussPlusPolynom[iPt],"SRB+QL"); //L option is log likelyhood: would supposedly be better in case of counts

    parGaussianParab_Sigma[iPt] = parGaussParab[iPt][2];
    parGaussianParab_Sigma_error[iPt] = GaussPlusPolynom[iPt]->GetParError(2);

    //give background the final parameters from that last signal+background fit
    GaussPlusPolynom[iPt]->GetParameters(&parGaussParab[iPt][0]);
    bkgparab[iPt]->SetParameters(&parGaussParab[iPt][3]);

    //plot Inv Mass histogram and the background and signal fits for each iPt bin
    H1D_DetectedV0s_O2_SmallpTinterval[iPt]->GetXaxis()->SetTitle("inv. mass (GeV/c^{2})");
    H1D_DetectedV0s_O2_SmallpTinterval[iPt]->GetYaxis()->SetTitle("Counts");
    gPad->DrawFrame(LowerLimitDisplay_Mass[ipart],0,UpperLimitDisplay_Mass[ipart],1.2*H1D_DetectedV0s_O2_SmallpTinterval[iPt]->GetMaximum());
    H1D_DetectedV0s_O2_SmallpTinterval[iPt]->Draw("E,same");
    Int_t icolor=0;
    // H1D_DetectedV0s_O2_SmallpTinterval[iPt]->SetMarkerStyle(markers[0]);
    // // use the same color for markers and lines
    // H1D_DetectedV0s_O2_SmallpTinterval[iPt]->SetMarkerColor(colors [icolor]);
    H1D_DetectedV0s_O2_SmallpTinterval[iPt]->SetLineColor  (colors [icolor]);
    // bkgparab[iPt]->Draw("same");
    // gauss[iPt]->Draw("same");

    // TLegend *legt = new TLegend( 0.62, 0.6, 0.87, 0.8, Form("pT [%.1f;%.1f]", pTbins[iPt], pTbins[iPt+1]));
    // legt->Draw();
    TLatex * text_pTbin = new TLatex (0.62,0.75,Form("pT [%.1f;%.1f] GeV", pTbins[iPt], pTbins[iPt+1]));
    text_pTbin->SetNDC(kTRUE); //remove if I want x,y in TLatex to be in the coordinate system of the histogram
    text_pTbin->Draw();
    TLatex * text_fitParam_mu = new TLatex (0.62,0.65,Form("#mu = %.5f GeV/c^{2}", parGaussParab[iPt][1]));
    text_fitParam_mu->SetNDC(kTRUE); //remove if I want x,y in TLatex to be in the coordinate system of the histogram
    text_fitParam_mu->Draw();
    TLatex * text_fitParam_sigma = new TLatex (0.62,0.55,Form("#sigma = %.5f GeV/c^{2}", parGaussParab[iPt][2]));
    text_fitParam_sigma->SetNDC(kTRUE); //remove if I want x,y in TLatex to be in the coordinate system of the histogram
    text_fitParam_sigma->Draw();

    // float xmin_InvMass = H1D_DetectedV0s_O2_SmallpTinterval[iPt]->GetXaxis()->GetXmin();
    // float xmax_InvMass = H1D_DetectedV0s_O2_SmallpTinterval[iPt]->GetXaxis()->GetXmax();
    cout << "xmin_InvMass= " << liminf[ipart] << ", xmax_InvMass= " << limsup[ipart] << endl;

    int ibin_liminf = H1D_DetectedV0s_O2_SmallpTinterval[iPt]->GetXaxis()->FindBin(liminf[ipart]);
    int ibin_limsup = H1D_DetectedV0s_O2_SmallpTinterval[iPt]->GetXaxis()->FindBin(limsup[ipart]);
    int nbinInvMass = ibin_limsup - ibin_liminf;
    cout << "nbinInvMass = " << nbinInvMass << endl;

    int totalCount = H1D_DetectedV0s_O2_SmallpTinterval[iPt]->Integral(ibin_liminf+1,ibin_limsup-1); //integral of a hist counts the number of bins, does not care for bin width: not true integral
    double backgroundCount = bkgparab[iPt]->Integral(liminf[ipart],limsup[ipart])*1./((limsup[ipart]-liminf[ipart])/nbinInvMass); //bkgparab is not an histogram but a function; here we need to convert it to a count
    //issue I actually need the integral of the positive part of the background function; for now kept simple integral as background is close to 0 compared to signal

    // float totalCount = H1D_DetectedV0s_O2_SmallpTinterval[iPt]->Integral(1,-1); //integral of a hist counts the number of bins, does not care for bin width: not true integral
    // float backgroundCount = bkgparab[iPt]->Integral(liminf[ipart],limsup[ipart])*1./((limsup[ipart]-liminf[ipart])/nbinInvMass); //bkgparab is not an histogram but a function; here we need to convert it to a count
    // //issue I actually need the integral of the positive part of the background function; for now kept simple integral as background is close to 0 compared to signal



    double SignalCountFit = totalCount - backgroundCount;
    cout << "totalCount= " << totalCount << ", backgroundCount= " << backgroundCount << endl;
    cout << "SignalCountFit= " << SignalCountFit << endl;

    double dpT = hPtDiff_RawSpectrum_O2->GetXaxis()->GetBinWidth(ibinPt);
    float drapidity = 1.5; //1.5
    double dN_dpT = SignalCountFit *1./dpT*1./drapidity;
    hPtDiff_RawSpectrum_O2->SetBinContent(ibinPt,dN_dpT);
    hPtDiff_RawSpectrum_O2->SetBinError(ibinPt,sqrt(SignalCountFit*1./dpT*1./drapidity)); //Issue with errors; this error of the O2 data spectrum is false, not right formula with fit
  }

  canvasMassvsPt->SaveAs("Fits_invMass_ptDiff_"+NamePart[ipart]+".pdf");
}

void PtDifferential_SigmaOfFit(TH1D* &hstat, TH1D* &hsyst, TH1D* &hsystCorr, Int_t ipart, TFile* file_O2Analysis, TString* SaveAs_Title, TString* &texXtitle, TString* &texYtitle, Double_t* pTbins, Int_t nbinpT, bool isMC) {

  // //HEP file
  // TFile* file_HEPAnalysis = new TFile("HEPData-ins881474-v1-root.root");
  // TH1D* hPtDiff_RawSpectrum_HEP[numPart]; //= (TH1D*)file_HEPAnalysis->Get("Table 1/Hist1D_y1");
  // TH1D* hPtDiff_RawSpectrum_HEP_statErrors[numPart]; // = (TH1D*)file_HEPAnalysis->Get("Table 1/Hist1D_y1_e1");
  // hPtDiff_RawSpectrum_HEP[0] = (TH1D*)file_HEPAnalysis->Get("Table 1/Hist1D_y1"); //K0S
  // hPtDiff_RawSpectrum_HEP_statErrors[0] = (TH1D*)file_HEPAnalysis->Get("Table 1/Hist1D_y1_e1");
  // hPtDiff_RawSpectrum_HEP[1] = (TH1D*)file_HEPAnalysis->Get("Table 2/Hist1D_y1"); //Lambda
  // hPtDiff_RawSpectrum_HEP_statErrors[1] = (TH1D*)file_HEPAnalysis->Get("Table 2/Hist1D_y1_e1");
  // hPtDiff_RawSpectrum_HEP[2] = (TH1D*)file_HEPAnalysis->Get("Table 3/Hist1D_y1"); //AntiLambda
  // hPtDiff_RawSpectrum_HEP_statErrors[2] = (TH1D*)file_HEPAnalysis->Get("Table 3/Hist1D_y1_e1");

  //O2 file
  TH1I* H1I_SelectedEventCount;
  TH3D* H3D_DetectedV0s_O2;

  if (isMC){
    H1I_SelectedEventCount = (TH1I*)file_O2Analysis->Get("lambdakzero-analysis-mc/hSelectedEventCounter");
    H3D_DetectedV0s_O2 = (TH3D*)file_O2Analysis->Get("lambdakzero-analysis-mc/h3dMass"+NamePart[ipart]);
  }
  else {
    H1I_SelectedEventCount = (TH1I*)file_O2Analysis->Get("lambdakzero-analysis/hSelectedEventCounter");
    H3D_DetectedV0s_O2 = (TH3D*)file_O2Analysis->Get("lambdakzero-analysis/h3dMass"+NamePart[ipart]);
  }
  int SelectedEventCount = H1I_SelectedEventCount->GetEntries();

  TH1D* H1D_DetectedV0s_PtAxis_O2 = (TH1D*)H3D_DetectedV0s_O2->ProjectionY(NamehistoInvMass[ipart],0,-1,0,-1);

  //Rebinning InvMass axis if necessary
  TH3D* H3D_DetectedV0s_O2_rebinnedZ = (TH3D*)H3D_DetectedV0s_O2->RebinZ(1,"H3D_DetectedV0s_O2_rebinnedZ");

  TH1D *H1D_DetectedV0s_O2_SmallpTinterval[nbinpT];
  TH1D *hPtDiff_RawSpectrum_O2 = new TH1D("hPtDiff_RawSpectrum_O2","hPtDiff_RawSpectrum_O2",nbinpT,pTbins);

  Double_t parGaussianParab_Sigma[nbinpT];  
  Double_t parGaussianParab_Sigma_error[nbinpT];  
  PtDifferential_Fit(hPtDiff_RawSpectrum_O2, H3D_DetectedV0s_O2, H3D_DetectedV0s_O2_rebinnedZ, ipart, pTbins, nbinpT, parGaussianParab_Sigma, parGaussianParab_Sigma_error);

  hstat = (TH1D*)hPtDiff_RawSpectrum_O2->Clone("hstat");
  hstat->Reset("M");
  int iPt = 0;
  for(int ibinPt = 1; ibinPt <= nbinpT; ibinPt++){
    iPt = ibinPt - 1;
    hstat->SetBinContent(ibinPt,parGaussianParab_Sigma[iPt]);
    hstat->SetBinError(ibinPt,parGaussianParab_Sigma_error[iPt]);
    
    cout << "hstat (ibinPt; Sigma) = (" << ibinPt << "; "<< parGaussianParab_Sigma[ibinPt] << ")" << endl;
  }

  //////Error Bars///////
  hsyst     = new TH1D("hsyst", "hsyst", 100, 0, 10);
  hsystCorr = new TH1D("hsystCorr", "hsystCorr", 100, 0, 10);
  hstat     ->Sumw2(); 
  hsyst     ->Sumw2();
  hsystCorr ->Sumw2(); 
  Int_t nbinx = hstat->GetNbinsX();
  
  for(Int_t ibinx = 1; ibinx <= nbinx; ibinx++){
    hsyst->SetBinContent (ibinx, hstat->GetBinContent(ibinx));
    hsyst->SetBinError   (ibinx, hstat->GetBinContent(ibinx)*0.08);
    hsystCorr->SetBinContent (ibinx, hstat->GetBinContent(ibinx));
    hsystCorr->SetBinError   (ibinx, hstat->GetBinContent(ibinx)*0.15);
  }

  *SaveAs_Title += "SigmaFit_pTdiff";
  texXtitle = texPtX;
  texYtitle = texSgimaGaussFit;
}




void Efficiency_O2data_allCountSignalRegion(TH1D* &hstat, TH1D* &hsyst, TH1D*&hsystCorr, Int_t ipart, TFile* file_O2Analysis, TString* SaveAs_Title, TString* &texXtitle, TString* &texYtitle, Double_t* pTbins, Int_t nbinpT) {

  //O2 file
  // TFile* file_O2Analysis = new TFile("AnalysisResults_WithSel8Cut.root");

  TH1I* H1I_SelectedEventCount = (TH1I*)file_O2Analysis->Get("lambdakzero-analysis-mc/hSelectedEventCounter");
  float SelectedEventCount = H1I_SelectedEventCount->GetEntries();

  TH3D* H3D_DetectedV0s_O2 = (TH3D*)file_O2Analysis->Get("lambdakzero-analysis-mc/h3dMass"+NamePart[ipart]);
  TH1D* H1D_DetectedV0s_PtAxis_O2 = (TH1D*)H3D_DetectedV0s_O2->ProjectionY("InvMass_K0s",0,-1,0,-1);
  // TH1D* hPtDiff_RawSpectrum_O2 = (TH1D*)H1D_DetectedV0s_PtAxis_O2->Clone("hPtDiff_RawSpectrum");
  // hPtDiff_RawSpectrum_O2->Reset("M"); //should be empty, only used that to get same Pt binning as used in analysis
 
  TH1D* TrueV0PtSpectrum = (TH1D*)file_O2Analysis->Get("lambdakzero-particle-count-mc/h"+NamePart[ipart]+"Count_PtDiff");
 
  //Mass window in which we count K0S
  float MassSignalInterval_lowEdge = 0.48;
  float MassSignalInterval_upEdge = 0.5;
  int MassSignalInterval_lowEdgeBin = H3D_DetectedV0s_O2->GetZaxis()->FindBin(MassSignalInterval_lowEdge);
  int MassSignalInterval_upEdgeBin = H3D_DetectedV0s_O2->GetZaxis()->FindBin(MassSignalInterval_upEdge);

  // Double_t xbins[9] = {0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.6, 2.0, 3.0};
  // Int_t nbinpT = 8;
  TH1D *H1D_DetectedV0s_O2_SmallpTinterval[nbinpT];
  TH1D *hPtDiff_RawSpectrum_O2 = new TH1D("hPtDiff_RawSpectrum_O2","hPtDiff_RawSpectrum_O2",nbinpT,pTbins);
  TH1D *TrueV0PtSpectrum_rebinned = new TH1D("TrueV0PtSpectrum_rebinned","TrueV0PtSpectrum_rebinned",nbinpT,pTbins);

  for(int ibinPt = 1; ibinPt <= nbinpT; ibinPt++){
    int iPt = ibinPt-1; //different convention for bins that start at 1 and arrays that start at 0

    int ibinPt_originalAxis_low = H3D_DetectedV0s_O2->GetYaxis()->FindBin(pTbins[iPt]);
    int ibinPt_originalAxis_up = H3D_DetectedV0s_O2->GetYaxis()->FindBin(pTbins[iPt+1]) - 1;
    cout << "ibinPt: " << ibinPt << "; originalAxis [low;up] = [" << ibinPt_originalAxis_low << ";" << ibinPt_originalAxis_up << "]" << endl;

    H1D_DetectedV0s_O2_SmallpTinterval[iPt] = (TH1D*)H3D_DetectedV0s_O2->ProjectionZ("InvMass_"+NamePart[ipart]+Form("_Pt%i", ibinPt),0,-1,ibinPt_originalAxis_low,ibinPt_originalAxis_up);

    float SignalCountRough = H1D_DetectedV0s_O2_SmallpTinterval[iPt]->Integral(MassSignalInterval_lowEdgeBin,MassSignalInterval_upEdgeBin);

    // float dpT = hPtDiff_RawSpectrum_O2->GetXaxis()->GetBinWidth(ibinPt);
    // float dN_dpT = SignalCountRough *1./dpT;
    hPtDiff_RawSpectrum_O2->SetBinContent(ibinPt,SignalCountRough);
    hPtDiff_RawSpectrum_O2->SetBinError(ibinPt,sqrt(SignalCountRough));

    int trueV0count_ibinpT = TrueV0PtSpectrum->Integral(ibinPt_originalAxis_low,ibinPt_originalAxis_up);
    TrueV0PtSpectrum_rebinned->SetBinContent(ibinPt,trueV0count_ibinpT);
    TrueV0PtSpectrum_rebinned->SetBinError(ibinPt,sqrt(trueV0count_ibinpT));

    cout << "(ibinPt, SignalCountRough, trueV0count_ibinpT) = (" << ibinPt << ", " << SignalCountRough << ", " << trueV0count_ibinpT << ")" << endl;
    cout << "efficiency = " << SignalCountRough*1./trueV0count_ibinpT << endl;
  }

  // H3D_DetectedV0s_O2->RebinY(10);
  // hPtDiff_RawSpectrum_O2->Rebin(10);

  // Int_t nbinPt = hPtDiff_RawSpectrum_O2->GetNbinsX();
  // for(Int_t ibinPt = 1; ibinPt <= nbinPt; ibinPt++){
  //   float SignalCountRough = H3D_DetectedV0s_O2->Integral(0,-1,ibinPt,ibinPt+1,MassSignalInterval_lowEdgeBin,MassSignalInterval_upEdgeBin);

  //   float dpT = hPtDiff_RawSpectrum_O2->GetXaxis()->GetBinWidth(ibinPt);
  //   float dN_dpT = SignalCountRough *1./dpT;
  //   hPtDiff_RawSpectrum_O2->SetBinContent(ibinPt,dN_dpT);
  //   hPtDiff_RawSpectrum_O2->SetBinError(ibinPt,sqrt(SignalCountRough*1./dpT));
  // }

  // hPtDiff_RawSpectrum_O2->Scale(1./SelectedEventCount);
  hstat = (TH1D*)hPtDiff_RawSpectrum_O2->Clone("hstat");
  hstat->Divide(hPtDiff_RawSpectrum_O2,TrueV0PtSpectrum_rebinned);


  //////Error Bars///////
  hsyst     = new TH1D("hsyst", "hsyst", 100, 0, 10);
  hsystCorr = new TH1D("hsystCorr", "hsystCorr", 100, 0, 10);
  hstat     ->Sumw2(); 
  hsyst     ->Sumw2();
  hsystCorr ->Sumw2(); 
  Int_t nbinx = hstat->GetNbinsX();
  
  for(Int_t ibinx = 1; ibinx <= nbinx; ibinx++){
    hsyst->SetBinContent (ibinx, hstat->GetBinContent(ibinx));
    hsyst->SetBinError   (ibinx, hstat->GetBinContent(ibinx)*0.08);
    hsystCorr->SetBinContent (ibinx, hstat->GetBinContent(ibinx));
    hsystCorr->SetBinError   (ibinx, hstat->GetBinContent(ibinx)*0.15);
  } 

  *SaveAs_Title += "Efficiency_O2data_allCountSignalRegion";
  texXtitle = texPtX;
  texYtitle = texEfficiency;
}



void Efficiency_O2MCdata_withFit(TH1D* &hstat, TH1D* &hsyst, TH1D*&hsystCorr, Int_t ipart, TFile* file_O2Analysis, TString* SaveAs_Title, TString* &texXtitle, TString* &texYtitle, Double_t* pTbins, Int_t nbinpT) {

  //O2 file
  TH1I* H1I_SelectedEventCount = (TH1I*)file_O2Analysis->Get("lambdakzero-analysis-mc/hSelectedEventCounter");
  float SelectedEventCount = H1I_SelectedEventCount->GetEntries();

  TH3D* H3D_DetectedV0s_O2 = (TH3D*)file_O2Analysis->Get("lambdakzero-analysis-mc/h3dMass"+NamePart[ipart]);
  TH1D* H1D_DetectedV0s_PtAxis_O2 = (TH1D*)H3D_DetectedV0s_O2->ProjectionY(NamehistoInvMass[ipart],0,-1,0,-1);

  TH1D* TrueV0PtSpectrum = (TH1D*)file_O2Analysis->Get("lambdakzero-particle-count-mc/h"+NamePart[ipart]+"Count_PtDiff");

  //Rebinning InvMass axis if necessary
  TH3D* H3D_DetectedV0s_O2_rebinnedZ = (TH3D*)H3D_DetectedV0s_O2->RebinZ(3,"H3D_DetectedV0s_O2_rebinnedZ");

  //Fit tools initialisation
  TF1 *gauss[nbinpT];
  TF1 *GaussPlusPolynom[nbinpT];
  TF1 *bkgparab[nbinpT];
  Float_t fFitResult[nbinpT];
  Double_t parGaussParab[nbinpT][6];  

  TH1D *H1D_DetectedV0s_O2_SmallpTinterval[nbinpT];
  TH1D *hPtDiff_RawSpectrum_O2 = new TH1D("hPtDiff_RawSpectrum_O2","hPtDiff_RawSpectrum_O2",nbinpT,pTbins);
  TH1D *TrueV0PtSpectrum_rebinned = new TH1D("TrueV0PtSpectrum_rebinned","TrueV0PtSpectrum_rebinned",nbinpT,pTbins);

  TCanvas *canvasMassvsPt = new TCanvas ("canvasMassPt", NamehistoInvMass[ipart], 1200, 800);
  // canvasMassvsPt->Divide(4,4);
  canvasMassvsPt->Divide(std::ceil(std::sqrt(nbinpT)),std::ceil(std::sqrt(nbinpT)));

  for(int ibinPt = 1; ibinPt <= nbinpT; ibinPt++){
    int iPt = ibinPt-1; //different convention for bins that start at 1 and arrays that start at 0

    canvasMassvsPt->cd(ibinPt); //ATTENTION: +1 because canvas divide starts counting from 1

    int ibinPt_originalAxis_low = H3D_DetectedV0s_O2->GetYaxis()->FindBin(pTbins[iPt]);
    int ibinPt_originalAxis_up = H3D_DetectedV0s_O2->GetYaxis()->FindBin(pTbins[iPt+1]) - 1;
    cout << "ibinPt: " << ibinPt << "; originalAxis [low;up] = [" << ibinPt_originalAxis_low << ";" << ibinPt_originalAxis_up << "]" << endl;

    H1D_DetectedV0s_O2_SmallpTinterval[iPt] = (TH1D*)H3D_DetectedV0s_O2_rebinnedZ->ProjectionZ("InvMass_"+NamePart[ipart]+Form("_Pt%i", ibinPt),0,-1,ibinPt_originalAxis_low,ibinPt_originalAxis_up);

    gauss[iPt] = new TF1("gauss_"+NamePart[ipart]+Form("_Pt%i", ibinPt), "gaus", min_range_signal[ipart], max_range_signal[ipart]);
    gauss[iPt]->SetLineColor(kRed);   
    gauss[iPt]->SetParName(0, "norm");
    gauss[iPt]->SetParName(1, "mean");
    gauss[iPt]->SetParName(2, "sigma");
    gauss[iPt]->SetParameters(1., MassPart[ipart], 0.0002);
    gauss[iPt]->SetParLimits(0, 0., 1.1*H1D_DetectedV0s_O2_SmallpTinterval[iPt]->GetBinContent(H1D_DetectedV0s_O2_SmallpTinterval[iPt]->GetMaximumBin()));
    gauss[iPt]->SetParLimits(1, min_range_signal[ipart], max_range_signal[ipart]);
    gauss[iPt]->SetParLimits(2, 0.0001, 0.05);

    bkgparab[iPt] = new TF1("parab_"+NamePart[ipart]+Form("_Pt%i", ibinPt), fparab, liminf[ipart], limsup[ipart], 4);
    bkgparab[iPt]->SetLineColor(kGreen);
    bkgparab[iPt]->FixParameter(3, ipart);

    GaussPlusPolynom[iPt] = new TF1("GaussPlusPolynom_"+NamePart[ipart]+Form("_Pt%i", ibinPt), "gaus(0)+pol2(3)", liminf[ipart], limsup[ipart]);
    GaussPlusPolynom[iPt]->SetLineColor(kAzure); 
    GaussPlusPolynom[iPt]->SetParName(0, "norm");
    GaussPlusPolynom[iPt]->SetParName(1, "mean");
    GaussPlusPolynom[iPt]->SetParName(2, "sigma");
    GaussPlusPolynom[iPt]->SetParName(3, "p0");
    GaussPlusPolynom[iPt]->SetParName(4, "p1");
    GaussPlusPolynom[iPt]->SetParName(5, "p2");

    GaussPlusPolynom[iPt]->SetParameters(1., MassPart[ipart], 0.0002);
    GaussPlusPolynom[iPt]->SetParLimits(0, 0., 1.1*H1D_DetectedV0s_O2_SmallpTinterval[iPt]->GetBinContent(H1D_DetectedV0s_O2_SmallpTinterval[iPt]->GetMaximumBin()));
    GaussPlusPolynom[iPt]->SetParLimits(1, min_range_signal[ipart], max_range_signal[ipart]);
    GaussPlusPolynom[iPt]->SetParLimits(2, 0.0001, 0.05);

    H1D_DetectedV0s_O2_SmallpTinterval[iPt]->Fit(gauss[iPt], "R0QL");
    reject = true; // for background fit: ignore signal range when fitting
    H1D_DetectedV0s_O2_SmallpTinterval[iPt]->Fit(bkgparab[iPt], "R0QL");
    reject = false;

    gauss[iPt]->GetParameters(&parGaussParab[iPt][0]);
    bkgparab[iPt]->GetParameters(&parGaussParab[iPt][3]);
    GaussPlusPolynom[iPt]->SetParameters(parGaussParab[iPt]);

    fFitResult[iPt] = H1D_DetectedV0s_O2_SmallpTinterval[iPt]->Fit(GaussPlusPolynom[iPt],"SRB+QL"); //L option is log likelyhood: would supposedly be better in case of counts

    //give background the final parameters from that last signal+background fit
    GaussPlusPolynom[iPt]->GetParameters(&parGaussParab[iPt][0]);
    bkgparab[iPt]->SetParameters(&parGaussParab[iPt][3]);

    //plot Inv Mass histogram and the background and signal fits for each iPt bin
    H1D_DetectedV0s_O2_SmallpTinterval[iPt]->GetXaxis()->SetTitle("inv. mass (GeV/c^{2})");
    H1D_DetectedV0s_O2_SmallpTinterval[iPt]->GetYaxis()->SetTitle("Counts");
    gPad->DrawFrame(LowerLimitDisplay_Mass[ipart],0,UpperLimitDisplay_Mass[ipart],1.2*H1D_DetectedV0s_O2_SmallpTinterval[iPt]->GetMaximum());
    H1D_DetectedV0s_O2_SmallpTinterval[iPt]->Draw("E,same");
    Int_t icolor=0;
    // H1D_DetectedV0s_O2_SmallpTinterval[iPt]->SetMarkerStyle(markers[0]);
    // // use the same color for markers and lines
    // H1D_DetectedV0s_O2_SmallpTinterval[iPt]->SetMarkerColor(colors [icolor]);
    H1D_DetectedV0s_O2_SmallpTinterval[iPt]->SetLineColor  (colors [icolor]);
    bkgparab[iPt]->Draw("same");
    // gauss[iPt]->Draw("same");

    TLatex * text_pTbin = new TLatex (0.62,0.75,Form("pT [%.1f;%.1f] GeV", pTbins[iPt], pTbins[iPt+1]));
    text_pTbin->SetNDC(kTRUE);
    text_pTbin->Draw();
    TLatex * text_fitParam_mu = new TLatex (0.62,0.65,Form("#mu = %.5f GeV/c^{2}", parGaussParab[iPt][1]));
    text_fitParam_mu->SetNDC(kTRUE);
    text_fitParam_mu->Draw();
    TLatex * text_fitParam_sigma = new TLatex (0.62,0.55,Form("#sigma = %.5f GeV/c^{2}", parGaussParab[iPt][2]));
    text_fitParam_sigma->SetNDC(kTRUE);
    text_fitParam_sigma->Draw();

    cout << "xmin_InvMass= " << liminf[ipart] << ", xmax_InvMass= " << limsup[ipart] << endl;
    int ibin_liminf = H1D_DetectedV0s_O2_SmallpTinterval[iPt]->GetXaxis()->FindBin(liminf[ipart]);
    int ibin_limsup = H1D_DetectedV0s_O2_SmallpTinterval[iPt]->GetXaxis()->FindBin(limsup[ipart]);
    int nbinInvMass = ibin_limsup - ibin_liminf;
    cout << "nbinInvMass = " << nbinInvMass << endl;

    int totalCount = H1D_DetectedV0s_O2_SmallpTinterval[iPt]->Integral(ibin_liminf+1,ibin_limsup-1); //integral of a hist counts the number of bins, does not care for bin width: not true integral
    double backgroundCount = bkgparab[iPt]->Integral(liminf[ipart],limsup[ipart])*1./((limsup[ipart]-liminf[ipart])/nbinInvMass); //bkgparab is not an histogram but a function; here we need to convert it to a count


    float SignalCountFit = totalCount - backgroundCount;
    cout << "totalCount= " << totalCount << ", backgroundCount= " << backgroundCount << endl;
    cout << "SignalCountFit= " << SignalCountFit << endl;

    hPtDiff_RawSpectrum_O2->SetBinContent(ibinPt,SignalCountFit);
    hPtDiff_RawSpectrum_O2->SetBinError(ibinPt,sqrt(SignalCountFit));

    int trueV0count_ibinpT = TrueV0PtSpectrum->Integral(ibinPt_originalAxis_low,ibinPt_originalAxis_up);
    TrueV0PtSpectrum_rebinned->SetBinContent(ibinPt,trueV0count_ibinpT);
    TrueV0PtSpectrum_rebinned->SetBinError(ibinPt,sqrt(trueV0count_ibinpT));

    cout << "(ibinPt, SignalCountFit, trueV0count_ibinpT) = (" << ibinPt << ", " << SignalCountFit << ", " << trueV0count_ibinpT << ")" << endl;
    cout << "efficiency = " << SignalCountFit*1./trueV0count_ibinpT << endl;
  }
  canvasMassvsPt->SaveAs("Fits_invMass_ptDiff_"+NamePart[ipart]+".pdf");

  // hPtDiff_RawSpectrum_O2->Scale(1./SelectedEventCount);
  hstat = (TH1D*)hPtDiff_RawSpectrum_O2->Clone("hstat");
  hstat->Divide(hPtDiff_RawSpectrum_O2,TrueV0PtSpectrum_rebinned);

  //////Error Bars///////
  hsyst     = new TH1D("hsyst", "hsyst", 100, 0, 10);
  hsystCorr = new TH1D("hsystCorr", "hsystCorr", 100, 0, 10);
  hstat     ->Sumw2(); 
  hsyst     ->Sumw2();
  hsystCorr ->Sumw2(); 
  Int_t nbinx = hstat->GetNbinsX();
  
  for(Int_t ibinx = 1; ibinx <= nbinx; ibinx++){
    hsyst->SetBinContent (ibinx, hstat->GetBinContent(ibinx));
    hsyst->SetBinError   (ibinx, hstat->GetBinContent(ibinx)*0.08);
    hsystCorr->SetBinContent (ibinx, hstat->GetBinContent(ibinx));
    hsystCorr->SetBinError   (ibinx, hstat->GetBinContent(ibinx)*0.15);
  } 

  *SaveAs_Title += "Efficiency_O2data_withFit";
  texXtitle = texPtX;
  texYtitle = texEfficiency;
}


void Efficiency_O2data_truePt_trueV0s(TH1D* &hstat, TH1D* &hsyst, TH1D*&hsystCorr, Int_t ipart, TFile* file_O2Analysis, TString* SaveAs_Title, TString* &texXtitle, TString* &texYtitle, Double_t* pTbins, Int_t nbinpT) {

  //O2 file
  // TFile* file_O2Analysis = new TFile("AnalysisResults_WithSel8Cut.root");

  TH1I* H1I_SelectedEventCount = (TH1I*)file_O2Analysis->Get("lambdakzero-analysis-mc/hSelectedEventCounter");
  float SelectedEventCount = H1I_SelectedEventCount->GetEntries();

  TH3D* H3D_DetectedV0s_O2 = (TH3D*)file_O2Analysis->Get("lambdakzero-analysis-mc/h3dMass"+NamePart[ipart]+"_MC_truePt");
  TH1D* H1D_DetectedV0s_PtAxis_O2 = (TH1D*)H3D_DetectedV0s_O2->ProjectionY("InvMass_K0s",0,-1,0,-1);
  // TH1D* hPtDiff_RawSpectrum_O2 = (TH1D*)H1D_DetectedV0s_PtAxis_O2->Clone("hPtDiff_RawSpectrum");
  // hPtDiff_RawSpectrum_O2->Reset("M"); //should be empty, only used that to get same Pt binning as used in analysis
 
  TH1D* TrueV0PtSpectrum = (TH1D*)file_O2Analysis->Get("lambdakzero-particle-count-mc/h"+NamePart[ipart]+"Count_PtDiff");

  // Double_t xbins[9] = {0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.6, 2.0, 3.0};
  // Int_t nbinpT = 8;
  TH1D *H1D_DetectedV0s_O2_SmallpTinterval[nbinpT];
  TH1D *hPtDiff_RawSpectrum_O2 = new TH1D("hPtDiff_RawSpectrum_O2","hPtDiff_RawSpectrum_O2",nbinpT,pTbins);
  TH1D *TrueV0PtSpectrum_rebinned = new TH1D("TrueV0PtSpectrum_rebinned","TrueV0PtSpectrum_rebinned",nbinpT,pTbins);

  for(int ibinPt = 1; ibinPt <= nbinpT; ibinPt++){
    int iPt = ibinPt-1; //different convention for bins that start at 1 and arrays that start at 0

    int ibinPt_originalAxis_low = H3D_DetectedV0s_O2->GetYaxis()->FindBin(pTbins[iPt]);
    int ibinPt_originalAxis_up = H3D_DetectedV0s_O2->GetYaxis()->FindBin(pTbins[iPt+1]) - 1;
    cout << "ibinPt: " << ibinPt << "; originalAxis [low;up] = [" << ibinPt_originalAxis_low << ";" << ibinPt_originalAxis_up << "]" << endl;

    H1D_DetectedV0s_O2_SmallpTinterval[iPt] = (TH1D*)H3D_DetectedV0s_O2->ProjectionZ("InvMass_"+NamePart[ipart]+Form("_Pt%i", ibinPt),0,-1,ibinPt_originalAxis_low,ibinPt_originalAxis_up);

    float SignalCountRough = H1D_DetectedV0s_O2_SmallpTinterval[iPt]->Integral(1,H1D_DetectedV0s_O2_SmallpTinterval[iPt]->GetNbinsX()); // avoids overflow bins 0 and -1

    // float dpT = hPtDiff_RawSpectrum_O2->GetXaxis()->GetBinWidth(ibinPt);
    // float dN_dpT = SignalCountRough *1./dpT;
    hPtDiff_RawSpectrum_O2->SetBinContent(ibinPt,SignalCountRough);
    hPtDiff_RawSpectrum_O2->SetBinError(ibinPt,sqrt(SignalCountRough));

    int trueV0count_ibinpT = TrueV0PtSpectrum->Integral(ibinPt_originalAxis_low,ibinPt_originalAxis_up);
    TrueV0PtSpectrum_rebinned->SetBinContent(ibinPt,trueV0count_ibinpT);
    TrueV0PtSpectrum_rebinned->SetBinError(ibinPt,sqrt(trueV0count_ibinpT));

    cout << "(ibinPt, SignalCountRough, trueV0count_ibinpT) = (" << ibinPt << ", " << SignalCountRough << ", " << trueV0count_ibinpT << ")" << endl;
    cout << "efficiency = " << SignalCountRough*1./trueV0count_ibinpT << endl;
  }

  // H3D_DetectedV0s_O2->RebinY(10);
  // hPtDiff_RawSpectrum_O2->Rebin(10);

  // Int_t nbinPt = hPtDiff_RawSpectrum_O2->GetNbinsX();
  // for(Int_t ibinPt = 1; ibinPt <= nbinPt; ibinPt++){
  //   float SignalCountRough = H3D_DetectedV0s_O2->Integral(0,-1,ibinPt,ibinPt+1,MassSignalInterval_lowEdgeBin,MassSignalInterval_upEdgeBin);

  //   float dpT = hPtDiff_RawSpectrum_O2->GetXaxis()->GetBinWidth(ibinPt);
  //   float dN_dpT = SignalCountRough *1./dpT;
  //   hPtDiff_RawSpectrum_O2->SetBinContent(ibinPt,dN_dpT);
  //   hPtDiff_RawSpectrum_O2->SetBinError(ibinPt,sqrt(SignalCountRough*1./dpT));
  // }

  // hPtDiff_RawSpectrum_O2->Scale(1./SelectedEventCount);
  hstat = (TH1D*)hPtDiff_RawSpectrum_O2->Clone("hstat");
  hstat->Divide(hPtDiff_RawSpectrum_O2,TrueV0PtSpectrum_rebinned);


  //////Error Bars///////
  hsyst     = new TH1D("hsyst", "hsyst", 100, 0, 10);
  hsystCorr = new TH1D("hsystCorr", "hsystCorr", 100, 0, 10);
  hstat     ->Sumw2(); 
  hsyst     ->Sumw2();
  hsystCorr ->Sumw2(); 
  Int_t nbinx = hstat->GetNbinsX();
  
  for(Int_t ibinx = 1; ibinx <= nbinx; ibinx++){
    hsyst->SetBinContent (ibinx, hstat->GetBinContent(ibinx));
    hsyst->SetBinError   (ibinx, hstat->GetBinContent(ibinx)*0.08);
    hsystCorr->SetBinContent (ibinx, hstat->GetBinContent(ibinx));
    hsystCorr->SetBinError   (ibinx, hstat->GetBinContent(ibinx)*0.15);
  } 

  *SaveAs_Title += "Efficiency_O2data_truePt_trueV0s";
  texXtitle = texPtX;
  texYtitle = texEfficiency;
}

void RawSpectrum_O2data_truePt_trueV0s(TH1D* &hstat, TH1D* &hsyst, TH1D*&hsystCorr, Int_t ipart, TFile* file_O2Analysis, TString* SaveAs_Title, TString* &texXtitle, TString* &texYtitle, Double_t* pTbins, Int_t nbinpT) {

  //O2 file
  // TFile* file_O2Analysis = new TFile("AnalysisResults_WithSel8Cut.root");

  TH1I* H1I_SelectedEventCount = (TH1I*)file_O2Analysis->Get("lambdakzero-analysis-mc/hSelectedEventCounter");
  float SelectedEventCount = H1I_SelectedEventCount->GetEntries();

  TH3D* H3D_DetectedV0s_O2 = (TH3D*)file_O2Analysis->Get("lambdakzero-analysis-mc/h3dMass"+NamePart[ipart]+"_MC_truePt");
  TH1D* H1D_DetectedV0s_PtAxis_O2 = (TH1D*)H3D_DetectedV0s_O2->ProjectionY("InvMass_K0s",0,-1,0,-1);
  // TH1D* hPtDiff_RawSpectrum_O2 = (TH1D*)H1D_DetectedV0s_PtAxis_O2->Clone("hPtDiff_RawSpectrum");
  // hPtDiff_RawSpectrum_O2->Reset("M"); //should be empty, only used that to get same Pt binning as used in analysis
 
  TH1D* TrueV0PtSpectrum = (TH1D*)file_O2Analysis->Get("lambdakzero-particle-count-mc/h"+NamePart[ipart]+"Count_PtDiff");

  // Double_t xbins[9] = {0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.6, 2.0, 3.0};
  // Int_t nbinpT = 8;
  TH1D *H1D_DetectedV0s_O2_SmallpTinterval[nbinpT];
  TH1D *hPtDiff_RawSpectrum_O2 = new TH1D("hPtDiff_RawSpectrum_O2","hPtDiff_RawSpectrum_O2",nbinpT,pTbins);
  TH1D *TrueV0PtSpectrum_rebinned = new TH1D("TrueV0PtSpectrum_rebinned","TrueV0PtSpectrum_rebinned",nbinpT,pTbins);

  for(int ibinPt = 1; ibinPt <= nbinpT; ibinPt++){
    int iPt = ibinPt-1; //different convention for bins that start at 1 and arrays that start at 0

    int ibinPt_originalAxis_low = H3D_DetectedV0s_O2->GetYaxis()->FindBin(pTbins[iPt]);
    int ibinPt_originalAxis_up = H3D_DetectedV0s_O2->GetYaxis()->FindBin(pTbins[iPt+1]) - 1;
    cout << "ibinPt: " << ibinPt << "; originalAxis [low;up] = [" << ibinPt_originalAxis_low << ";" << ibinPt_originalAxis_up << "]" << endl;

    H1D_DetectedV0s_O2_SmallpTinterval[iPt] = (TH1D*)H3D_DetectedV0s_O2->ProjectionZ("InvMass_"+NamePart[ipart]+Form("_Pt%i", ibinPt),0,-1,ibinPt_originalAxis_low,ibinPt_originalAxis_up);

    float SignalCountRough = H1D_DetectedV0s_O2_SmallpTinterval[iPt]->Integral(1,H1D_DetectedV0s_O2_SmallpTinterval[iPt]->GetNbinsX());

    float dpT = hPtDiff_RawSpectrum_O2->GetXaxis()->GetBinWidth(ibinPt);
    float dN_dpT = SignalCountRough *1./dpT;
    hPtDiff_RawSpectrum_O2->SetBinContent(ibinPt,dN_dpT);
    hPtDiff_RawSpectrum_O2->SetBinError(ibinPt,sqrt(SignalCountRough*1./dpT));

    int trueV0count_ibinpT = TrueV0PtSpectrum->Integral(ibinPt_originalAxis_low,ibinPt_originalAxis_up);
    TrueV0PtSpectrum_rebinned->SetBinContent(ibinPt,trueV0count_ibinpT);
    TrueV0PtSpectrum_rebinned->SetBinError(ibinPt,sqrt(trueV0count_ibinpT));

    cout << "(ibinPt, SignalCountRough, trueV0count_ibinpT) = (" << ibinPt << ", " << SignalCountRough << ", " << trueV0count_ibinpT << ")" << endl;
    cout << "yield = " << SignalCountRough*1./dpT*1./SelectedEventCount << endl;
  }

  // H3D_DetectedV0s_O2->RebinY(10);
  // hPtDiff_RawSpectrum_O2->Rebin(10);

  // Int_t nbinPt = hPtDiff_RawSpectrum_O2->GetNbinsX();
  // for(Int_t ibinPt = 1; ibinPt <= nbinPt; ibinPt++){
  //   float SignalCountRough = H3D_DetectedV0s_O2->Integral(0,-1,ibinPt,ibinPt+1,MassSignalInterval_lowEdgeBin,MassSignalInterval_upEdgeBin);

  //   float dpT = hPtDiff_RawSpectrum_O2->GetXaxis()->GetBinWidth(ibinPt);
  //   float dN_dpT = SignalCountRough *1./dpT;
  //   hPtDiff_RawSpectrum_O2->SetBinContent(ibinPt,dN_dpT);
  //   hPtDiff_RawSpectrum_O2->SetBinError(ibinPt,sqrt(SignalCountRough*1./dpT));
  // }

  hPtDiff_RawSpectrum_O2->Scale(1./SelectedEventCount);
  hstat = (TH1D*)hPtDiff_RawSpectrum_O2->Clone("hstat");


  //////Error Bars///////
  hsyst     = new TH1D("hsyst", "hsyst", 100, 0, 10);
  hsystCorr = new TH1D("hsystCorr", "hsystCorr", 100, 0, 10);
  hstat     ->Sumw2(); 
  hsyst     ->Sumw2();
  hsystCorr ->Sumw2(); 
  Int_t nbinx = hstat->GetNbinsX();
  
  for(Int_t ibinx = 1; ibinx <= nbinx; ibinx++){
    hsyst->SetBinContent (ibinx, hstat->GetBinContent(ibinx));
    hsyst->SetBinError   (ibinx, hstat->GetBinContent(ibinx)*0.08);
    hsystCorr->SetBinContent (ibinx, hstat->GetBinContent(ibinx));
    hsystCorr->SetBinError   (ibinx, hstat->GetBinContent(ibinx)*0.15);
  } 

  *SaveAs_Title += "RawSpectrum_O2data_truePt_trueV0s";
  texXtitle = texPtX;
  texYtitle = texEfficiency;
}




void V0Daughters_TrackingEfficiencies(TH1D* &hstat, TH1D* &hsyst, TH1D*&hsystCorr, Int_t ipart, TFile* file_O2Analysis, TString* SaveAs_Title, Double_t* pTbins, Int_t nbinpT) {

  //O2 file
  // TFile* file_O2Analysis = new TFile("AnalysisResults_WithSel8Cut.root");

  TH1D* H1D_K0S_gen_pionPlus = (TH1D*)file_O2Analysis->Get("v0daughters-tracking-efficiency/hPtGenerated_PionPlusFromV0K0S");
  TH1D* H1D_K0S_gen_pionNeg = (TH1D*)file_O2Analysis->Get("v0daughters-tracking-efficiency/hPtGenerated_PionNegFromV0K0S");
  TH1D* H1D_K0S_rec_pionPlus = (TH1D*)file_O2Analysis->Get("v0daughters-tracking-efficiency/hPtReconstructed_PionPlusFromV0K0S");
  TH1D* H1D_K0S_rec_pionNeg = (TH1D*)file_O2Analysis->Get("v0daughters-tracking-efficiency/hPtReconstructed_PionNegFromV0K0S");

  TH1D* H1D_Lambda_gen_pionNeg = (TH1D*)file_O2Analysis->Get("v0daughters-tracking-efficiency/hPtGenerated_PionNegFromV0Lambda");
  TH1D* H1D_Lambda_gen_protonPlus = (TH1D*)file_O2Analysis->Get("v0daughters-tracking-efficiency/hPtGenerated_ProtonPlusFromV0Lambda");
  TH1D* H1D_Lambda_rec_pionNeg = (TH1D*)file_O2Analysis->Get("v0daughters-tracking-efficiency/hPtReconstructed_PionNegFromV0Lambda");
  TH1D* H1D_Lambda_rec_protonPlus = (TH1D*)file_O2Analysis->Get("v0daughters-tracking-efficiency/hPtReconstructed_ProtonPlusFromV0Lambda");

  TH1D* H1D_AntiLambda_gen_pionPlus = (TH1D*)file_O2Analysis->Get("v0daughters-tracking-efficiency/hPtGenerated_PionPlusFromV0AntiLambda");
  TH1D* H1D_AntiLambda_gen_protonNeg = (TH1D*)file_O2Analysis->Get("v0daughters-tracking-efficiency/hPtGenerated_ProtonNegFromV0AntiLambda");
  TH1D* H1D_AntiLambda_rec_pionPlus = (TH1D*)file_O2Analysis->Get("v0daughters-tracking-efficiency/hPtReconstructed_PionPlusFromV0AntiLambda");
  TH1D* H1D_AntiLambda_rec_protonNeg = (TH1D*)file_O2Analysis->Get("v0daughters-tracking-efficiency/hPtReconstructed_ProtonNegFromV0AntiLambda");

  TH1D* H1D_K0S_gen_pionPlus_rebin = (TH1D*)H1D_K0S_gen_pionPlus->Rebin(nbinpT,"H1D_K0S_gen_pionPlus_rebin",pTbins);
  TH1D* H1D_K0S_gen_pionNeg_rebin = (TH1D*)H1D_K0S_gen_pionNeg->Rebin(nbinpT,"H1D_K0S_gen_pionNeg_rebin",pTbins);
  TH1D* H1D_K0S_rec_pionPlus_rebin = (TH1D*)H1D_K0S_rec_pionPlus->Rebin(nbinpT,"H1D_K0S_rec_pionPlus_rebin",pTbins);
  TH1D* H1D_K0S_rec_pionNeg_rebin = (TH1D*)H1D_K0S_rec_pionNeg->Rebin(nbinpT,"H1D_K0S_rec_pionNeg_rebin",pTbins);

  TH1D* H1D_Lambda_gen_pionNeg_rebin = (TH1D*)H1D_Lambda_gen_pionNeg->Rebin(nbinpT,"H1D_Lambda_gen_pionNeg_rebin",pTbins);
  TH1D* H1D_Lambda_gen_protonPlus_rebin = (TH1D*)H1D_Lambda_gen_protonPlus->Rebin(nbinpT,"H1D_Lambda_gen_protonPlus_rebin",pTbins);
  TH1D* H1D_Lambda_rec_pionNeg_rebin = (TH1D*)H1D_Lambda_rec_pionNeg->Rebin(nbinpT,"H1D_Lambda_rec_pionNeg_rebin",pTbins);
  TH1D* H1D_Lambda_rec_protonPlus_rebin = (TH1D*)H1D_Lambda_rec_protonPlus->Rebin(nbinpT,"H1D_Lambda_rec_protonPlus_rebin",pTbins);

  TH1D* H1D_AntiLambda_gen_pionPlus_rebin = (TH1D*)H1D_AntiLambda_gen_pionPlus->Rebin(nbinpT,"H1D_AntiLambda_gen_pionPlus_rebin",pTbins);
  TH1D* H1D_AntiLambda_gen_protonNeg_rebin = (TH1D*)H1D_AntiLambda_gen_protonNeg->Rebin(nbinpT,"H1D_AntiLambda_gen_protonNeg_rebin",pTbins);
  TH1D* H1D_AntiLambda_rec_pionPlus_rebin = (TH1D*)H1D_AntiLambda_rec_pionPlus->Rebin(nbinpT,"H1D_AntiLambda_rec_pionPlus_rebin",pTbins);
  TH1D* H1D_AntiLambda_rec_protonNeg_rebin = (TH1D*)H1D_AntiLambda_rec_protonNeg->Rebin(nbinpT,"H1D_AntiLambda_rec_protonNeg_rebin",pTbins);

  TCanvas *canvasV0DaughterTrackEff = new TCanvas ("canvasV0DaughterTrackEff", "canvasV0DaughterTrackEff", 1200, 800);
  canvasV0DaughterTrackEff->Divide(2,1);
  Int_t icolor=0;

  if (ipart == 0) {
    TH1D* H1D_Eff_K0S_pionPlus = (TH1D*)H1D_K0S_gen_pionPlus_rebin->Clone("H1D_Eff_K0S_pionPlus");
    H1D_Eff_K0S_pionPlus->Divide(H1D_K0S_rec_pionPlus,H1D_K0S_gen_pionPlus);
    canvasV0DaughterTrackEff->cd(1);
    H1D_Eff_K0S_pionPlus->Draw("E");
    TH1D* H1D_Eff_K0S_pionNeg = (TH1D*)H1D_K0S_gen_pionNeg_rebin->Clone("H1D_Eff_K0S_pionNeg");
    H1D_Eff_K0S_pionNeg->Divide(H1D_K0S_rec_pionNeg,H1D_K0S_gen_pionNeg);
    canvasV0DaughterTrackEff->cd(2);
    H1D_Eff_K0S_pionNeg->Draw("E");

    H1D_Eff_K0S_pionPlus->SetMarkerStyle(markers[0]);
    H1D_Eff_K0S_pionPlus->SetMarkerColor(colors [icolor]);
    H1D_Eff_K0S_pionPlus->SetLineColor  (colors [icolor]);
    H1D_Eff_K0S_pionNeg->SetMarkerStyle(markers[0]);
    H1D_Eff_K0S_pionNeg->SetMarkerColor(colors [icolor]);
    H1D_Eff_K0S_pionNeg->SetLineColor  (colors [icolor]);
  }

  if (ipart == 0) {
    TH1D* H1D_Eff_Lambda_pionNeg = (TH1D*)H1D_Lambda_gen_pionNeg_rebin->Clone("H1D_Eff_Lambda_pionNeg");
    H1D_Eff_Lambda_pionNeg->Divide(H1D_Lambda_rec_pionNeg,H1D_Lambda_gen_pionNeg);
    canvasV0DaughterTrackEff->cd(1);
    H1D_Eff_Lambda_pionNeg->Draw("E");
    TH1D* H1D_Eff_Lambda_protonPlus = (TH1D*)H1D_Lambda_gen_protonPlus_rebin->Clone("H1D_Eff_Lambda_protonPlus");
    H1D_Eff_Lambda_protonPlus->Divide(H1D_Lambda_rec_protonPlus,H1D_Lambda_gen_protonPlus);
    canvasV0DaughterTrackEff->cd(2);
    H1D_Eff_Lambda_protonPlus->Draw("E");

    H1D_Eff_Lambda_pionNeg->SetMarkerStyle(markers[0]);
    H1D_Eff_Lambda_pionNeg->SetMarkerColor(colors [icolor]);
    H1D_Eff_Lambda_pionNeg->SetLineColor  (colors [icolor]);
    H1D_Eff_Lambda_protonPlus->SetMarkerStyle(markers[0]);
    H1D_Eff_Lambda_protonPlus->SetMarkerColor(colors [icolor]);
    H1D_Eff_Lambda_protonPlus->SetLineColor  (colors [icolor]);
  }

  if (ipart == 0) {
    TH1D* H1D_Eff_AntiLambda_pionPlus = (TH1D*)H1D_AntiLambda_gen_pionPlus_rebin->Clone("H1D_Eff_AntiLambda_pionPlus");
    H1D_Eff_AntiLambda_pionPlus->Divide(H1D_AntiLambda_rec_pionPlus,H1D_AntiLambda_gen_pionPlus);
    canvasV0DaughterTrackEff->cd(1);
    H1D_Eff_AntiLambda_pionPlus->Draw("E");
    TH1D* H1D_Eff_AntiLambda_protonNeg = (TH1D*)H1D_AntiLambda_gen_protonNeg_rebin->Clone("H1D_Eff_AntiLambda_protonNeg");
    H1D_Eff_AntiLambda_protonNeg->Divide(H1D_AntiLambda_rec_protonNeg,H1D_AntiLambda_gen_protonNeg);
    canvasV0DaughterTrackEff->cd(2);
    H1D_Eff_AntiLambda_protonNeg->Draw("E");

    H1D_Eff_AntiLambda_pionPlus->SetMarkerStyle(markers[0]);
    H1D_Eff_AntiLambda_pionPlus->SetMarkerColor(colors [icolor]);
    H1D_Eff_AntiLambda_pionPlus->SetLineColor  (colors [icolor]);
    H1D_Eff_AntiLambda_protonNeg->SetMarkerStyle(markers[0]);
    H1D_Eff_AntiLambda_protonNeg->SetMarkerColor(colors [icolor]);
    H1D_Eff_AntiLambda_protonNeg->SetLineColor  (colors [icolor]);
  }


  // // hPtDiff_RawSpectrum_O2->Scale(1./SelectedEventCount);
  // hstat = (TH1D*)H1D_Eff_AntiLambda_protonNeg->Clone("hstat");
  // hstat->Divide(H1D_Eff_AntiLambda_protonNeg,H1D_Eff_AntiLambda_protonNeg);
  hstat     = new TH1D("hsyst", "hsyst", 100, 0, 10); 

  //////Error Bars///////
  hsyst     = new TH1D("hsyst", "hsyst", 100, 0, 10);
  hsystCorr = new TH1D("hsystCorr", "hsystCorr", 100, 0, 10);
  hstat     ->Sumw2(); 
  hsyst     ->Sumw2();
  hsystCorr ->Sumw2(); 
  Int_t nbinx = hstat->GetNbinsX();
  
  for(Int_t ibinx = 1; ibinx <= nbinx; ibinx++){
    hsyst->SetBinContent (ibinx, hstat->GetBinContent(ibinx));
    hsyst->SetBinError   (ibinx, hstat->GetBinContent(ibinx)*0.08);
    hsystCorr->SetBinContent (ibinx, hstat->GetBinContent(ibinx));
    hsystCorr->SetBinError   (ibinx, hstat->GetBinContent(ibinx)*0.15);
  } 

  *SaveAs_Title += "ToBeDeleted";
}

void V0Daughters_TrackingEfficiencies_HistogramProcessing(TH1D* &daughter1, TH1D* &daughter2, Int_t ipart, TFile* file_O2Analysis, Double_t* pTbinsProtons, Int_t nbinpTProtons, Double_t* pTbinsPions, Int_t nbinpTPions) {

  TH1D* H1D_K0S_gen_pionPlus = (TH1D*)file_O2Analysis->Get("v0daughters-tracking-efficiency/hPtGenerated_PionPlusFromV0K0S");
  TH1D* H1D_K0S_gen_pionNeg = (TH1D*)file_O2Analysis->Get("v0daughters-tracking-efficiency/hPtGenerated_PionNegFromV0K0S");
  TH1D* H1D_K0S_rec_pionPlus = (TH1D*)file_O2Analysis->Get("v0daughters-tracking-efficiency/hPtReconstructed_PionPlusFromV0K0S");
  TH1D* H1D_K0S_rec_pionNeg = (TH1D*)file_O2Analysis->Get("v0daughters-tracking-efficiency/hPtReconstructed_PionNegFromV0K0S");
  H1D_K0S_gen_pionPlus->Sumw2();
  H1D_K0S_gen_pionNeg->Sumw2();
  H1D_K0S_rec_pionPlus->Sumw2();
  H1D_K0S_rec_pionNeg->Sumw2();

  TH1D* H1D_Lambda_gen_pionNeg = (TH1D*)file_O2Analysis->Get("v0daughters-tracking-efficiency/hPtGenerated_PionNegFromV0Lambda");
  TH1D* H1D_Lambda_gen_protonPlus = (TH1D*)file_O2Analysis->Get("v0daughters-tracking-efficiency/hPtGenerated_ProtonPlusFromV0Lambda");
  TH1D* H1D_Lambda_rec_pionNeg = (TH1D*)file_O2Analysis->Get("v0daughters-tracking-efficiency/hPtReconstructed_PionNegFromV0Lambda");
  TH1D* H1D_Lambda_rec_protonPlus = (TH1D*)file_O2Analysis->Get("v0daughters-tracking-efficiency/hPtReconstructed_ProtonPlusFromV0Lambda");
  H1D_Lambda_gen_pionNeg->Sumw2();
  H1D_Lambda_gen_protonPlus->Sumw2();
  H1D_Lambda_rec_pionNeg->Sumw2();
  H1D_Lambda_rec_protonPlus->Sumw2();

  TH1D* H1D_AntiLambda_gen_pionPlus = (TH1D*)file_O2Analysis->Get("v0daughters-tracking-efficiency/hPtGenerated_PionPlusFromV0AntiLambda");
  TH1D* H1D_AntiLambda_gen_protonNeg = (TH1D*)file_O2Analysis->Get("v0daughters-tracking-efficiency/hPtGenerated_ProtonNegFromV0AntiLambda");
  TH1D* H1D_AntiLambda_rec_pionPlus = (TH1D*)file_O2Analysis->Get("v0daughters-tracking-efficiency/hPtReconstructed_PionPlusFromV0AntiLambda");
  TH1D* H1D_AntiLambda_rec_protonNeg = (TH1D*)file_O2Analysis->Get("v0daughters-tracking-efficiency/hPtReconstructed_ProtonNegFromV0AntiLambda");
  H1D_AntiLambda_gen_pionPlus->Sumw2();
  H1D_AntiLambda_gen_protonNeg->Sumw2();
  H1D_AntiLambda_rec_pionPlus->Sumw2();
  H1D_AntiLambda_rec_protonNeg->Sumw2();

  TH1D* H1D_K0S_gen_pionPlus_rebin = (TH1D*)H1D_K0S_gen_pionPlus->Rebin(nbinpTPions,"H1D_K0S_gen_pionPlus_rebin",pTbinsPions);
  TH1D* H1D_K0S_gen_pionNeg_rebin = (TH1D*)H1D_K0S_gen_pionNeg->Rebin(nbinpTPions,"H1D_K0S_gen_pionNeg_rebin",pTbinsPions);
  TH1D* H1D_K0S_rec_pionPlus_rebin = (TH1D*)H1D_K0S_rec_pionPlus->Rebin(nbinpTPions,"H1D_K0S_rec_pionPlus_rebin",pTbinsPions);
  TH1D* H1D_K0S_rec_pionNeg_rebin = (TH1D*)H1D_K0S_rec_pionNeg->Rebin(nbinpTPions,"H1D_K0S_rec_pionNeg_rebin",pTbinsPions);

  TH1D* H1D_Lambda_gen_pionNeg_rebin = (TH1D*)H1D_Lambda_gen_pionNeg->Rebin(nbinpTPions,"H1D_Lambda_gen_pionNeg_rebin",pTbinsPions);
  TH1D* H1D_Lambda_gen_protonPlus_rebin = (TH1D*)H1D_Lambda_gen_protonPlus->Rebin(nbinpTProtons,"H1D_Lambda_gen_protonPlus_rebin",pTbinsProtons);
  TH1D* H1D_Lambda_rec_pionNeg_rebin = (TH1D*)H1D_Lambda_rec_pionNeg->Rebin(nbinpTPions,"H1D_Lambda_rec_pionNeg_rebin",pTbinsPions);
  TH1D* H1D_Lambda_rec_protonPlus_rebin = (TH1D*)H1D_Lambda_rec_protonPlus->Rebin(nbinpTProtons,"H1D_Lambda_rec_protonPlus_rebin",pTbinsProtons);

  TH1D* H1D_AntiLambda_gen_pionPlus_rebin = (TH1D*)H1D_AntiLambda_gen_pionPlus->Rebin(nbinpTPions,"H1D_AntiLambda_gen_pionPlus_rebin",pTbinsPions);
  TH1D* H1D_AntiLambda_gen_protonNeg_rebin = (TH1D*)H1D_AntiLambda_gen_protonNeg->Rebin(nbinpTProtons,"H1D_AntiLambda_gen_protonNeg_rebin",pTbinsProtons);
  TH1D* H1D_AntiLambda_rec_pionPlus_rebin = (TH1D*)H1D_AntiLambda_rec_pionPlus->Rebin(nbinpTPions,"H1D_AntiLambda_rec_pionPlus_rebin",pTbinsPions);
  TH1D* H1D_AntiLambda_rec_protonNeg_rebin = (TH1D*)H1D_AntiLambda_rec_protonNeg->Rebin(nbinpTProtons,"H1D_AntiLambda_rec_protonNeg_rebin",pTbinsProtons);

  if (ipart == 0) {
    daughter1 = (TH1D*)H1D_K0S_gen_pionPlus_rebin->Clone("H1D_Eff_K0S_pionPlus");
    daughter1->Divide(H1D_K0S_rec_pionPlus_rebin,H1D_K0S_gen_pionPlus_rebin);
    daughter2 = (TH1D*)H1D_K0S_gen_pionNeg_rebin->Clone("H1D_Eff_K0S_pionNeg");
    daughter2->Divide(H1D_K0S_rec_pionNeg_rebin,H1D_K0S_gen_pionNeg_rebin); 
  }

  if (ipart == 1) {
    daughter1 = (TH1D*)H1D_Lambda_gen_pionNeg_rebin->Clone("H1D_Eff_Lambda_pionNeg");
    daughter1->Divide(H1D_Lambda_rec_pionNeg_rebin,H1D_Lambda_gen_pionNeg_rebin);
    daughter2 = (TH1D*)H1D_Lambda_gen_protonPlus_rebin->Clone("H1D_Eff_Lambda_protonPlus");
    daughter2->Divide(H1D_Lambda_rec_protonPlus_rebin,H1D_Lambda_gen_protonPlus_rebin);
  }

  if (ipart == 2) {
    daughter1 = (TH1D*)H1D_AntiLambda_gen_pionPlus_rebin->Clone("H1D_Eff_AntiLambda_pionPlus");
    daughter1->Divide(H1D_AntiLambda_rec_pionPlus_rebin,H1D_AntiLambda_gen_pionPlus_rebin);
    daughter2 = (TH1D*)H1D_AntiLambda_gen_protonNeg_rebin->Clone("H1D_Eff_AntiLambda_protonNeg");
    daughter2->Divide(H1D_AntiLambda_rec_protonNeg_rebin,H1D_AntiLambda_gen_protonNeg_rebin);
  }
}


void V0Daughters_TrackingEfficiencies_Daughter1(TH1D* &hstat, TH1D* &hsyst, TH1D*&hsystCorr, Int_t ipart, TFile* file_O2Analysis, TString* SaveAs_Title, TString* &texXtitle, TString* &texYtitle, TString* &Extra, Double_t* pTbinsProtons, Int_t nbinpTProtons, Double_t* pTbinsPions, Int_t nbinpTPions) {

  TH1D* daughter1 = new TH1D("daughter1", "daughter1", 100, 0, 10);
  TH1D* daughter2 = new TH1D("daughter2", "daughter2", 100, 0, 10);
  V0Daughters_TrackingEfficiencies_HistogramProcessing(daughter1, daughter2, ipart, file_O2Analysis, pTbinsProtons, nbinpTProtons, pTbinsPions, nbinpTPions);
  // // hPtDiff_RawSpectrum_O2->Scale(1./SelectedEventCount);
  hstat = (TH1D*)daughter1->Clone("hstat");

  //////Error Bars///////
  hsyst     = new TH1D("hsyst", "hsyst", 100, 0, 10);
  hsystCorr = new TH1D("hsystCorr", "hsystCorr", 100, 0, 10);
  hstat     ->Sumw2(); 
  hsyst     ->Sumw2();
  hsystCorr ->Sumw2(); 
  Int_t nbinx = hstat->GetNbinsX();
  
  for(Int_t ibinx = 1; ibinx <= nbinx; ibinx++){
    hsyst->SetBinContent (ibinx, hstat->GetBinContent(ibinx));
    hsyst->SetBinError   (ibinx, hstat->GetBinContent(ibinx)*0.08);
    hsystCorr->SetBinContent (ibinx, hstat->GetBinContent(ibinx));
    hsystCorr->SetBinError   (ibinx, hstat->GetBinContent(ibinx)*0.15);
  } 

  if (ipart == 0) {
    *SaveAs_Title += "TrackEff_K0S_pionPlus";
    Extra = texK0S_daughter_PiPlus;
  }
  if (ipart == 1) {
    *SaveAs_Title += "TrackEff_Lambda_pionNeg";
    Extra = texLambda_daughter_PiMinus;
  }
  if (ipart == 2) {
    *SaveAs_Title += "TrackEff_AntiLambda_pionPlus";
    Extra = texAntiLambda_daughter_PiPlus;
  }

  texXtitle = texPtX;
  texYtitle = texDaughterRecoEfficiency;
}

void V0Daughters_TrackingEfficiencies_Daughter2(TH1D* &hstat, TH1D* &hsyst, TH1D*&hsystCorr, Int_t ipart, TFile* file_O2Analysis, TString* SaveAs_Title, TString* &texXtitle, TString* &texYtitle, TString* &Extra, Double_t* pTbinsProtons, Int_t nbinpTProtons, Double_t* pTbinsPions, Int_t nbinpTPions) {

  TH1D* daughter1 = new TH1D("daughter1", "daughter1", 100, 0, 10);
  TH1D* daughter2 = new TH1D("daughter2", "daughter2", 100, 0, 10);
  V0Daughters_TrackingEfficiencies_HistogramProcessing(daughter1, daughter2, ipart, file_O2Analysis, pTbinsProtons, nbinpTProtons, pTbinsPions, nbinpTPions);
  // // hPtDiff_RawSpectrum_O2->Scale(1./SelectedEventCount);
  hstat = (TH1D*)daughter2->Clone("hstat");;

  //////Error Bars///////
  hsyst     = new TH1D("hsyst", "hsyst", 100, 0, 10);
  hsystCorr = new TH1D("hsystCorr", "hsystCorr", 100, 0, 10);
  hstat     ->Sumw2(); 
  hsyst     ->Sumw2();
  hsystCorr ->Sumw2(); 
  Int_t nbinx = hstat->GetNbinsX();
  
  for(Int_t ibinx = 1; ibinx <= nbinx; ibinx++){
    hsyst->SetBinContent (ibinx, hstat->GetBinContent(ibinx));
    hsyst->SetBinError   (ibinx, hstat->GetBinContent(ibinx)*0.08);
    hsystCorr->SetBinContent (ibinx, hstat->GetBinContent(ibinx));
    hsystCorr->SetBinError   (ibinx, hstat->GetBinContent(ibinx)*0.15);
  } 

  if (ipart == 0) {
    *SaveAs_Title += "TrackEff_K0S_pionNeg";
    Extra = texK0S_daughter_PiMinus;
  }
  if (ipart == 1) {
    *SaveAs_Title += "TrackEff_Lambda_protonPlus";
    Extra = texLambda_daughter_protonPlus;
  }
  if (ipart == 2) {
    *SaveAs_Title += "TrackEff_AntiLambda_protonNeg";
    Extra = texAntiLambda_daughter_protonMinus;
  }

  texXtitle = texPtX;
  texYtitle = texDaughterRecoEfficiency;
}




void pT_Spectrum_postAnalyserCuts(TH1D* &hstat, TH1D* &hsyst, TH1D*&hsystCorr, Int_t ipart, TFile* file_O2Analysis, TString* SaveAs_Title, TString* &texXtitle, TString* &texYtitle, Double_t* pTbins, Int_t nbinpT, bool isMC) {

  //O2 file
  // TFile* file_O2Analysis = new TFile("AnalysisResults_WithSel8Cut.root");
  TH1I* H1I_SelectedEventCount;
  TH3D* H3D_DetectedV0s_O2;

  if (isMC) {
    H1I_SelectedEventCount = (TH1I*)file_O2Analysis->Get("lambdakzero-analysis-mc/hSelectedEventCounter");
    H3D_DetectedV0s_O2 = (TH3D*)file_O2Analysis->Get("lambdakzero-analysis-mc/h3dMass"+NamePart[ipart]);
  }
  else {
    H1I_SelectedEventCount = (TH1I*)file_O2Analysis->Get("lambdakzero-analysis/hSelectedEventCounter");
    H3D_DetectedV0s_O2 = (TH3D*)file_O2Analysis->Get("lambdakzero-analysis/h3dMass"+NamePart[ipart]);
  }

  float SelectedEventCount = H1I_SelectedEventCount->GetEntries();

  TH1D* H1D_DetectedV0s_PtAxis_O2 = (TH1D*)H3D_DetectedV0s_O2->ProjectionY("InvMass_K0s",0,-1,0,-1);
  // TH1D* hPtDiff_RawSpectrum_O2 = (TH1D*)H1D_DetectedV0s_PtAxis_O2->Clone("hPtDiff_RawSpectrum");
  // hPtDiff_RawSpectrum_O2->Reset("M"); //should be empty, only used that to get same Pt binning as used in analysis
 
  TH1D* H1D_DetectedV0s_PtAxis_O2_rebin = (TH1D*)H1D_DetectedV0s_PtAxis_O2->Rebin(nbinpT,"H1D_DetectedV0s_PtAxis_O2_rebin",pTbins); 

  for(int ibinPt = 1; ibinPt <= nbinpT; ibinPt++){

    double dpT = H1D_DetectedV0s_PtAxis_O2_rebin->GetXaxis()->GetBinWidth(ibinPt);
    double Count = H1D_DetectedV0s_PtAxis_O2_rebin->GetBinContent(ibinPt);

    H1D_DetectedV0s_PtAxis_O2_rebin->SetBinContent(ibinPt,Count*1./dpT);
    H1D_DetectedV0s_PtAxis_O2_rebin->SetBinError(ibinPt,sqrt(Count*1./dpT));
  }

  H1D_DetectedV0s_PtAxis_O2_rebin->Scale(1./SelectedEventCount);
  hstat = (TH1D*)H1D_DetectedV0s_PtAxis_O2_rebin->Clone("hstat");


  //////Error Bars///////
  hsyst     = new TH1D("hsyst", "hsyst", 100, 0, 10);
  hsystCorr = new TH1D("hsystCorr", "hsystCorr", 100, 0, 10);
  hstat     ->Sumw2(); 
  hsyst     ->Sumw2();
  hsystCorr ->Sumw2(); 
  Int_t nbinx = hstat->GetNbinsX();
  
  for(Int_t ibinx = 1; ibinx <= nbinx; ibinx++){
    hsyst->SetBinContent (ibinx, hstat->GetBinContent(ibinx));
    hsyst->SetBinError   (ibinx, hstat->GetBinContent(ibinx)*0.08);
    hsystCorr->SetBinContent (ibinx, hstat->GetBinContent(ibinx));
    hsystCorr->SetBinError   (ibinx, hstat->GetBinContent(ibinx)*0.15);
  } 

  *SaveAs_Title += "pT_Spectrum_DirtySignal_postAnalyserCuts";
  texXtitle = texPtX;
  texYtitle = texptDifferentialYield;
}

void pT_Spectrum_preAnalyserCutsK0S(TH1D* &hstat, TH1D* &hsyst, TH1D*&hsystCorr, Int_t ipart, TFile* file_O2Analysis, TString* SaveAs_Title, TString* &texXtitle, TString* &texYtitle, Double_t* pTbins, Int_t nbinpT, bool isMC) {

  int useless = ipart;
  //O2 file
  TFile* file_O2AnalysisQA = new TFile("QAResults.root");

  TH1I* H1I_SelectedEventCount;

  if (isMC) {
    H1I_SelectedEventCount = (TH1I*)file_O2Analysis->Get("lambdakzero-analysis-mc/hSelectedEventCounter");
  }
  else {
    H1I_SelectedEventCount = (TH1I*)file_O2Analysis->Get("lambdakzero-analysis/hSelectedEventCounter");
  }
  float SelectedEventCount = H1I_SelectedEventCount->GetEntries();

  TH2D* H3D_DetectedV0s_O2 = (TH2D*)file_O2AnalysisQA->Get("v0cascades-q-a/histos-V0/InvMassK0S");
  TH1D* H1D_DetectedV0s_PtAxis_O2 = (TH1D*)H3D_DetectedV0s_O2->ProjectionX("InvMass_K0s",0,-1);
  // TH1D* hPtDiff_RawSpectrum_O2 = (TH1D*)H1D_DetectedV0s_PtAxis_O2->Clone("hPtDiff_RawSpectrum");
  // hPtDiff_RawSpectrum_O2->Reset("M"); //should be empty, only used that to get same Pt binning as used in analysis
 
  TH1D* H1D_DetectedV0s_PtAxis_O2_rebin = (TH1D*)H1D_DetectedV0s_PtAxis_O2->Rebin(nbinpT,"H1D_DetectedV0s_PtAxis_O2_rebin",pTbins); 

  for(int ibinPt = 1; ibinPt <= nbinpT; ibinPt++){

    double dpT = H1D_DetectedV0s_PtAxis_O2_rebin->GetXaxis()->GetBinWidth(ibinPt);
    double Count = H1D_DetectedV0s_PtAxis_O2_rebin->GetBinContent(ibinPt);

    H1D_DetectedV0s_PtAxis_O2_rebin->SetBinContent(ibinPt,Count*1./dpT);
    H1D_DetectedV0s_PtAxis_O2_rebin->SetBinError(ibinPt,sqrt(Count*1./dpT));
  }

  H1D_DetectedV0s_PtAxis_O2_rebin->Scale(1./SelectedEventCount);
  hstat = (TH1D*)H1D_DetectedV0s_PtAxis_O2_rebin->Clone("hstat");


  //////Error Bars///////
  hsyst     = new TH1D("hsyst", "hsyst", 100, 0, 10);
  hsystCorr = new TH1D("hsystCorr", "hsystCorr", 100, 0, 10);
  hstat     ->Sumw2(); 
  hsyst     ->Sumw2();
  hsystCorr ->Sumw2(); 
  Int_t nbinx = hstat->GetNbinsX();
  
  for(Int_t ibinx = 1; ibinx <= nbinx; ibinx++){
    hsyst->SetBinContent (ibinx, hstat->GetBinContent(ibinx));
    hsyst->SetBinError   (ibinx, hstat->GetBinContent(ibinx)*0.08);
    hsystCorr->SetBinContent (ibinx, hstat->GetBinContent(ibinx));
    hsystCorr->SetBinError   (ibinx, hstat->GetBinContent(ibinx)*0.15);
  } 

  *SaveAs_Title += "pT_Spectrum_DirtySignal_preAnalyserCuts";
  texXtitle = texPtX;
  texYtitle = texptDifferentialYield;
}

void Get_McEfficiency_vsPt(TH1D* &hMcEfficiency_vsPt, TH3D* H3D_detectedV0s_TruePt, TH1D* TrueV0PtSpectrum, Int_t ipart, Double_t* pTbins, Int_t nbinpT, Int_t ibinXaxisCut) {

  // what is fed as input: TH3D* H3D_detectedV0s_TruePt = (TH3D*)file_O2Analysis->Get("lambdakzero-analysis-mc/h3dMass"+NamePart[ipart]+"_MC_truePt");
  // what is fed as input: TH1D* TrueV0PtSpectrum = (TH1D*)file_O2Analysis->Get("lambdakzero-particle-count-mc/h"+NamePart[ipart]+"Count_PtDiff");


  TH1D *H1D_detectedV0s_TruePt_SmallpTinterval[nbinpT];
  TH1D *hMcSignalCount_vsPt = new TH1D("hMcSignalCount_vsPt","hMcSignalCount_vsPt",nbinpT,pTbins);
  TH1D *TrueV0PtSpectrum_AnalysisBins = new TH1D("TrueV0PtSpectrum_AnalysisBins","TrueV0PtSpectrum_AnalysisBins",nbinpT,pTbins);

  for(int ibinPt = 1; ibinPt <= nbinpT; ibinPt++){
    int iPt = ibinPt-1; //different convention for bins that start at 1 and arrays that start at 0

    int ibinPt_originalAxis_low = H3D_detectedV0s_TruePt->GetYaxis()->FindBin(pTbins[iPt]);
    int ibinPt_originalAxis_up = H3D_detectedV0s_TruePt->GetYaxis()->FindBin(pTbins[iPt+1]) - 1;
    cout << "ibinPt: " << ibinPt << "; originalAxis [low;up] = [" << ibinPt_originalAxis_low << ";" << ibinPt_originalAxis_up << "]" << endl;

    H1D_detectedV0s_TruePt_SmallpTinterval[iPt] = (TH1D*)H3D_detectedV0s_TruePt->ProjectionZ("InvMass_"+NamePart[ipart]+Form("_Pt%i", ibinPt),0,ibinXaxisCut,ibinPt_originalAxis_low,ibinPt_originalAxis_up);

    float SignalCountRough = H1D_detectedV0s_TruePt_SmallpTinterval[iPt]->Integral(1,H1D_detectedV0s_TruePt_SmallpTinterval[iPt]->GetNbinsX()); // avoids overflow bins 0 and -1

    hMcSignalCount_vsPt->SetBinContent(ibinPt,SignalCountRough);
    hMcSignalCount_vsPt->SetBinError(ibinPt,sqrt(SignalCountRough));

    int trueV0count_ibinpT = TrueV0PtSpectrum->Integral(ibinPt_originalAxis_low,ibinPt_originalAxis_up);
    TrueV0PtSpectrum_AnalysisBins->SetBinContent(ibinPt,trueV0count_ibinpT);
    TrueV0PtSpectrum_AnalysisBins->SetBinError(ibinPt,sqrt(trueV0count_ibinpT));

    cout << "(ibinPt, SignalCountRough, trueV0count_ibinpT) = (" << ibinPt << ", " << SignalCountRough << ", " << trueV0count_ibinpT << ")" << endl;
    cout << "efficiency = " << SignalCountRough*1./trueV0count_ibinpT << endl;
  }

  hMcEfficiency_vsPt = (TH1D*)hMcSignalCount_vsPt->Clone("hMcEfficiency_vsPt");
  hMcEfficiency_vsPt->Divide(hMcSignalCount_vsPt,TrueV0PtSpectrum_AnalysisBins);
  delete hMcSignalCount_vsPt; // avoids memory leaks when Get_McEfficiency_vsPt is called in a loop
  delete TrueV0PtSpectrum_AnalysisBins; // avoids memory leaks when Get_McEfficiency_vsPt is called in a loop
}

void Get_McSignalBackgroundRatio_vsPt(TH1D* &hSignalBackgroundRatio_vsPt, TH3D* H3D_detectedV0s_DataLike, Int_t ipart, Double_t* pTbins, Int_t nbinpT, Int_t ibinXaxisCut) {
  //Looks within the signal area of the fit, 6Sigma

  // what is fed as input: TH3D* H3D_detectedV0s_DataLike = (TH3D*)file_O2Analysis->Get("lambdakzero-analysis-mc/h3dMass"+NamePart[ipart]);

  ////////////////////////////////// Fit Start //////////////////////////////////
  //Fit tools initialisation
  TF1 *gauss[nbinpT];
  TF1 *GaussPlusPolynom[nbinpT];
  TF1 *bkgparab[nbinpT];
  Float_t fFitResult[nbinpT];
  Double_t parGaussParab[nbinpT][6];  

  TH1D *H1D_detectedV0s_DataLike_SmallpTinterval[nbinpT];
  TH1D *hSignalCountFit_vsPt = new TH1D("hSignalCountFit_vsPt","hSignalCountFit_vsPt",nbinpT,pTbins);
  TH1D *hBackgroundCountFit_vsPt = new TH1D("hBackgroundCountFit_vsPt","hBackgroundCountFit_vsPt",nbinpT,pTbins);

  for(int ibinPt = 1; ibinPt <= nbinpT; ibinPt++){
    int iPt = ibinPt-1; //different convention for bins that start at 1 and arrays that start at 0

    int ibinPt_originalAxis_low = H3D_detectedV0s_DataLike->GetYaxis()->FindBin(pTbins[iPt]);
    int ibinPt_originalAxis_up = H3D_detectedV0s_DataLike->GetYaxis()->FindBin(pTbins[iPt+1]) - 1;
    cout << "ibinPt: " << ibinPt << "; originalAxis [low;up] = [" << ibinPt_originalAxis_low << ";" << ibinPt_originalAxis_up << "]" << endl;

    H1D_detectedV0s_DataLike_SmallpTinterval[iPt] = (TH1D*)H3D_detectedV0s_DataLike->ProjectionZ("InvMass_"+NamePart[ipart]+Form("_Pt%i", ibinPt),0,ibinXaxisCut,ibinPt_originalAxis_low,ibinPt_originalAxis_up);

    gauss[iPt] = new TF1("gauss_"+NamePart[ipart]+Form("_Pt%i", ibinPt), "gaus", min_range_signal[ipart], max_range_signal[ipart]);
    gauss[iPt]->SetLineColor(kRed);   
    gauss[iPt]->SetParName(0, "norm");
    gauss[iPt]->SetParName(1, "mean");
    gauss[iPt]->SetParName(2, "sigma");
    gauss[iPt]->SetParameters(1., MassPart[ipart], 0.0002);
    gauss[iPt]->SetParLimits(0, 0., 1.1*H1D_detectedV0s_DataLike_SmallpTinterval[iPt]->GetBinContent(H1D_detectedV0s_DataLike_SmallpTinterval[iPt]->GetMaximumBin()));
    gauss[iPt]->SetParLimits(1, min_range_signal[ipart], max_range_signal[ipart]);
    gauss[iPt]->SetParLimits(2, 0.0001, 0.05);

    bkgparab[iPt] = new TF1("parab_"+NamePart[ipart]+Form("_Pt%i", ibinPt), fparab, liminf[ipart], limsup[ipart], 4);
    bkgparab[iPt]->SetLineColor(kGreen);
    bkgparab[iPt]->FixParameter(3, ipart);

    GaussPlusPolynom[iPt] = new TF1("GaussPlusPolynom_"+NamePart[ipart]+Form("_Pt%i", ibinPt), "gaus(0)+pol2(3)", liminf[ipart], limsup[ipart]);
    GaussPlusPolynom[iPt]->SetLineColor(kAzure); 
    GaussPlusPolynom[iPt]->SetParName(0, "norm");
    GaussPlusPolynom[iPt]->SetParName(1, "mean");
    GaussPlusPolynom[iPt]->SetParName(2, "sigma");
    GaussPlusPolynom[iPt]->SetParName(3, "p0");
    GaussPlusPolynom[iPt]->SetParName(4, "p1");
    GaussPlusPolynom[iPt]->SetParName(5, "p2");

    GaussPlusPolynom[iPt]->SetParameters(1., MassPart[ipart], 0.0002);
    GaussPlusPolynom[iPt]->SetParLimits(0, 0., 1.1*H1D_detectedV0s_DataLike_SmallpTinterval[iPt]->GetBinContent(H1D_detectedV0s_DataLike_SmallpTinterval[iPt]->GetMaximumBin()));
    GaussPlusPolynom[iPt]->SetParLimits(1, min_range_signal[ipart], max_range_signal[ipart]);
    GaussPlusPolynom[iPt]->SetParLimits(2, 0.0001, 0.05);

    H1D_detectedV0s_DataLike_SmallpTinterval[iPt]->Fit(gauss[iPt], "R0QL");
    reject = true; // for background fit: ignore signal range when fitting
    H1D_detectedV0s_DataLike_SmallpTinterval[iPt]->Fit(bkgparab[iPt], "R0QL");
    reject = false;

    gauss[iPt]->GetParameters(&parGaussParab[iPt][0]);
    bkgparab[iPt]->GetParameters(&parGaussParab[iPt][3]);
    GaussPlusPolynom[iPt]->SetParameters(parGaussParab[iPt]);

    fFitResult[iPt] = H1D_detectedV0s_DataLike_SmallpTinterval[iPt]->Fit(GaussPlusPolynom[iPt],"SRB+QL"); //L option is log likelyhood: would supposedly be better in case of counts

    //give background the final parameters from that last signal+background fit
    GaussPlusPolynom[iPt]->GetParameters(&parGaussParab[iPt][0]);
    bkgparab[iPt]->SetParameters(&parGaussParab[iPt][3]);

    cout << "xmin_InvMass= " << liminf[ipart] << ", xmax_InvMass= " << limsup[ipart] << endl;
    int ibin_liminf = H1D_detectedV0s_DataLike_SmallpTinterval[iPt]->GetXaxis()->FindBin(liminf[ipart]);
    int ibin_limsup = H1D_detectedV0s_DataLike_SmallpTinterval[iPt]->GetXaxis()->FindBin(limsup[ipart]);
    int nbinInvMass = ibin_limsup - ibin_liminf;
    cout << "nbinInvMass = " << nbinInvMass << endl;

    int totalCount = H1D_detectedV0s_DataLike_SmallpTinterval[iPt]->Integral(ibin_liminf+1,ibin_limsup-1); //integral of a hist counts the number of bins, does not care for bin width: not true integral

    double backgroundCount = bkgparab[iPt]->Integral(liminf[ipart],limsup[ipart])*1./((limsup[ipart]-liminf[ipart])/nbinInvMass); //bkgparab is not an histogram but a function; here we need to convert it to a count
    hBackgroundCountFit_vsPt->SetBinContent(ibinPt,backgroundCount);
    hBackgroundCountFit_vsPt->SetBinError(ibinPt,sqrt(backgroundCount));

    float SignalCountFit = totalCount - backgroundCount;
    hSignalCountFit_vsPt->SetBinContent(ibinPt,SignalCountFit);
    hSignalCountFit_vsPt->SetBinError(ibinPt,sqrt(SignalCountFit));
    // cout << "totalCount= " << totalCount << ", backgroundCount= " << backgroundCount << endl;
    // cout << "SignalCountFit= " << SignalCountFit << endl;

    // cout << "(ibinPt, SignalCountFit, trueV0count_ibinpT) = (" << ibinPt << ", " << SignalCountFit << ", " << trueV0count_ibinpT << ")" << endl;
    // cout << "efficiency = " << SignalCountFit*1./trueV0count_ibinpT << endl;
  }
  ////////////////////////////////// Fit End //////////////////////////////////
  hSignalBackgroundRatio_vsPt = (TH1D*)hSignalCountFit_vsPt->Clone("hSignalBackgroundRatio_vsPt");
  hSignalBackgroundRatio_vsPt->Divide(hSignalCountFit_vsPt,hBackgroundCountFit_vsPt);

  delete hSignalCountFit_vsPt; // avoids memory leaks when Get_McEfficiency_vsPt is called in a loop
  delete hBackgroundCountFit_vsPt; // avoids memory leaks when Get_McEfficiency_vsPt is called in a loop
}


void Efficiency_DcaScan(TH1D* &hstat, TH1D* &hsyst, TH1D*&hsystCorr, Int_t ipart, TFile* file_O2Analysis, TString* SaveAs_Title, TString* &texXtitle, TString* &texYtitle, Double_t* pTbins, Int_t nbinpT, Double_t maxDCAcut) {

  TH3D* H3D_detectedV0s_TruePt_DcaHist = (TH3D*)file_O2Analysis->Get("lambdakzero-analysis-mc/h3dMass"+NamePart[ipart]+"Dca_MC_truePt");
  TH3D* H3D_detectedV0s_DataLike_DcaHist = (TH3D*)file_O2Analysis->Get("lambdakzero-analysis-mc/h3dMass"+NamePart[ipart]+"Dca");
  TH1D* TrueV0PtSpectrum = (TH1D*)file_O2Analysis->Get("lambdakzero-particle-count-mc/h"+NamePart[ipart]+"Count_PtDiff");

  Int_t maxDCAcutBin = H3D_detectedV0s_TruePt_DcaHist->GetXaxis()->FindBin(maxDCAcut);
  TH1D *hMcEfficiency_vsPt[maxDCAcutBin];
  TH1D *hSignalBackgroundRatio_vsPt[maxDCAcutBin];
  TH1D *hMcEfficiency_vsDca[nbinpT];
  TH1D *hSignalBackgroundRatio_vsDca[nbinpT];

  //Initialise hMcEfficiency_vsDca histograms for each pT bin
  for(int ibinPt = 1; ibinPt <= nbinpT; ibinPt++){
    int itablePt = ibinPt-1;
    // hMcEfficiency_vsDca[itablePt] = new TH1D(Form("hMcEfficiency_vsDca_Pt%i", ibinPt),Form("hMcEfficiency_vsDca_Pt%i", ibinPt),maxDCAcutBin,DCAbins_cut);
    hMcEfficiency_vsDca[itablePt] = (TH1D*)H3D_detectedV0s_TruePt_DcaHist->ProjectionX(Form("hMcEfficiency_vsDca_Pt%i", ibinPt),0,maxDCAcutBin,0,-1);
    hMcEfficiency_vsDca[itablePt]->Reset("M");

    hSignalBackgroundRatio_vsDca[itablePt] = (TH1D*)H3D_detectedV0s_DataLike_DcaHist->ProjectionX(Form("hSignalBackgroundRatio_vsDca_Pt%i", ibinPt),0,maxDCAcutBin,0,-1);
    hSignalBackgroundRatio_vsDca[itablePt]->Reset("M");
  }

  //loops over DCA cuts to get the "efficiency vs pT" histograms and use them to fill "efficiency vs DCAcut" 
  for(Int_t ibinDcaCut = 1; ibinDcaCut <= maxDCAcutBin; ibinDcaCut++){
    int itableDcaCut = ibinDcaCut-1;
    Get_McEfficiency_vsPt(hMcEfficiency_vsPt[itableDcaCut], H3D_detectedV0s_TruePt_DcaHist, TrueV0PtSpectrum, ipart, pTbins, nbinpT, ibinDcaCut);
    Get_McSignalBackgroundRatio_vsPt(hSignalBackgroundRatio_vsPt[itableDcaCut], H3D_detectedV0s_DataLike_DcaHist, ipart, pTbins, nbinpT, ibinDcaCut);
    for(int ibinPt = 1; ibinPt <= nbinpT; ibinPt++){
      int itablePt = ibinPt-1;

      hMcEfficiency_vsDca[itablePt]->SetBinContent(ibinDcaCut,hMcEfficiency_vsPt[itableDcaCut]->GetBinContent(ibinPt));
      hMcEfficiency_vsDca[itablePt]->SetBinError(ibinDcaCut,hMcEfficiency_vsPt[itableDcaCut]->GetBinError(ibinPt));

      hSignalBackgroundRatio_vsDca[itablePt]->SetBinContent(ibinDcaCut,hSignalBackgroundRatio_vsPt[itableDcaCut]->GetBinContent(ibinPt));
      hSignalBackgroundRatio_vsDca[itablePt]->SetBinError(ibinDcaCut,hSignalBackgroundRatio_vsPt[itableDcaCut]->GetBinError(ibinPt));
    }
  }

  //Plots the "efficiency vs dca" for each pT bin
  TCanvas *canvasDcaScanvsPt_eff = new TCanvas ("canvasDcaScanvsPt_eff", NamehistoInvMass[ipart], 800, 600);
  // TCanvas *canvasDcaScanvsPt_eff = new TCanvas ("canvasDcaScanvsPt_eff", NamehistoInvMass[ipart], 1600, 1200);
  canvasDcaScanvsPt_eff->Divide(std::ceil(std::sqrt(nbinpT)),std::ceil(std::sqrt(nbinpT)));
  TH1 *hFrame_eff[nbinpT];
  for(int ibinPt = 1; ibinPt <= nbinpT; ibinPt++){
    int itablePt = ibinPt-1; //different convention for bins that start at 1 and arrays that start at 0

    canvasDcaScanvsPt_eff->cd(ibinPt); //ATTENTION: +1 because canvas divide starts counting from 1

    hFrame_eff[ibinPt] = canvasDcaScanvsPt_eff->DrawFrame(0,0,maxDCAcut,1.45*hMcEfficiency_vsDca[itablePt]->GetMaximum());
    hFrame_eff[ibinPt]->SetXTitle("DCA cut (cm)");
    hFrame_eff[ibinPt]->SetYTitle("Efficiency");
    hMcEfficiency_vsDca[itablePt]->Draw("E,same");
  }
  canvasDcaScanvsPt_eff->SaveAs("Efficiency_vsPt_"+NamePart[ipart]+".pdf");

  //Plots the "SignalToBackgroundRatio vs dca" for each pT bin
  TCanvas *canvasDcaScanvsPt_SBratio = new TCanvas ("canvasDcaScanvsPt_SBratio", NamehistoInvMass[ipart], 800, 600);
  // TCanvas *canvasDcaScanvsPt_SBratio = new TCanvas ("canvasDcaScanvsPt_SBratio", NamehistoInvMass[ipart], 1600, 1200);
  canvasDcaScanvsPt_SBratio->Divide(std::ceil(std::sqrt(nbinpT)),std::ceil(std::sqrt(nbinpT)));
  TH1 *hFrame_SBratio[nbinpT];
  for(int ibinPt = 1; ibinPt <= nbinpT; ibinPt++){
    int itablePt = ibinPt-1; //different convention for bins that start at 1 and arrays that start at 0

    canvasDcaScanvsPt_SBratio->cd(ibinPt); //ATTENTION: +1 because canvas divide starts counting from 1

    hFrame_SBratio[ibinPt] = canvasDcaScanvsPt_SBratio->DrawFrame(0,0,maxDCAcut,1.45*hSignalBackgroundRatio_vsDca[itablePt]->GetMaximum());
    hFrame_SBratio[ibinPt]->SetXTitle("DCA cut (cm)");
    hFrame_SBratio[ibinPt]->SetYTitle("Signal to Background Ratio");
    hSignalBackgroundRatio_vsDca[itablePt]->Draw("E,same");
  }
  canvasDcaScanvsPt_SBratio->SaveAs("SignalBackgroundRatio_vsPt_"+NamePart[ipart]+".pdf");

  ///////// standard stuff just to be able to use same global function at beginning of this file:

  hstat = (TH1D*)hMcEfficiency_vsDca[0]->Clone("hstat");
  hstat->Reset("M");
  // hstat->Divide(hPtDiff_RawSpectrum_O2,hPtDiff_RawSpectrum_HEP[ipart]);

  //////Error Bars///////
  hsyst     = new TH1D("hsyst", "hsyst", 100, 0, 10);
  hsystCorr = new TH1D("hsystCorr", "hsystCorr", 100, 0, 10);
  hstat     ->Sumw2(); 
  hsyst     ->Sumw2();
  hsystCorr ->Sumw2(); 
  Int_t nbinx = hstat->GetNbinsX();
  
  for(Int_t ibinx = 1; ibinx <= nbinx; ibinx++){
    hsyst->SetBinContent (ibinx, hstat->GetBinContent(ibinx));
    hsyst->SetBinError   (ibinx, hstat->GetBinContent(ibinx)*0.08);
    hsystCorr->SetBinContent (ibinx, hstat->GetBinContent(ibinx));
    hsystCorr->SetBinError   (ibinx, hstat->GetBinContent(ibinx)*0.15);
  }

  *SaveAs_Title += "ToBeDeleted";
  texXtitle = texPtX;
  texYtitle = texptDifferentialYield_HEPratio_pseudoEff;

}


// void Systematics_Graphs_cutVariation(TH1D* &hstat, TH1D* &hsyst, TH1D*&hsystCorr, Int_t ipart, TFile* file_O2Analysis_ref, TFile* file_O2Analysis_CutVariation_array, TString* SaveAs_Title, TString* &texXtitle, TString* &texYtitle, Double_t* pTbins, Int_t nbinpT, bool isMC) {
//   // saves the systematics for particle ipart from file file_O2Analysis into hsyst, and produces a graph that will be ploted (hstat, x-axis title texXtitle and y-axis title texYtitle with pTbins)
//   // saves with SaveAs_Title name

//   // barlow check

//   TH1I* H1I_SelectedEventCount;
//   TH3D* H3D_DetectedV0s_O2;

//   H1I_SelectedEventCount = (TH1I*)file_O2Analysis_ref->Get("lambdakzero-analysis-mc/hSelectedEventCounter");
//   H3D_DetectedV0s_ref = (TH3D*)file_O2Analysis_ref->Get("lambdakzero-analysis-mc/h3dMass"+NamePart[ipart]);

//   float SelectedEventCount = H1I_SelectedEventCount->GetEntries();

//   for(int cutIterator = 0; ibinPt <= 3; ibinPt++){
//     H3D_DetectedV0s_cutVariation = (TH3D*)file_O2Analysis_CutVariation_array[cutIterator]->Get("lambdakzero-analysis-mc/h3dMass"+NamePart[ipart]);

//     TH1D* H1D_Efficiency = GetEfficiency(H3D_DetectedV0s_cutVariation);
//   }


//   //////////////////// ONGOING WORK
//   TH1D* H1D_DetectedV0s_PtAxis_O2 = (TH1D*)H3D_DetectedV0s_O2->ProjectionY("InvMass_K0s",0,-1,0,-1);
//   // TH1D* hPtDiff_RawSpectrum_O2 = (TH1D*)H1D_DetectedV0s_PtAxis_O2->Clone("hPtDiff_RawSpectrum");
//   // hPtDiff_RawSpectrum_O2->Reset("M"); //should be empty, only used that to get same Pt binning as used in analysis
 
//   TH1D* H1D_DetectedV0s_PtAxis_O2_rebin = (TH1D*)H1D_DetectedV0s_PtAxis_O2->Rebin(nbinpT,"H1D_DetectedV0s_PtAxis_O2_rebin",pTbins); 

//   for(int ibinPt = 1; ibinPt <= nbinpT; ibinPt++){

//     double dpT = H1D_DetectedV0s_PtAxis_O2_rebin->GetXaxis()->GetBinWidth(ibinPt);
//     double Count = H1D_DetectedV0s_PtAxis_O2_rebin->GetBinContent(ibinPt);

//     H1D_DetectedV0s_PtAxis_O2_rebin->SetBinContent(ibinPt,Count*1./dpT);
//     H1D_DetectedV0s_PtAxis_O2_rebin->SetBinError(ibinPt,sqrt(Count*1./dpT));
//   }

//   H1D_DetectedV0s_PtAxis_O2_rebin->Scale(1./SelectedEventCount);
//   hstat = (TH1D*)H1D_DetectedV0s_PtAxis_O2_rebin->Clone("hstat");


//   //////Error Bars///////
//   hsyst     = new TH1D("hsyst", "hsyst", 100, 0, 10);
//   hsystCorr = new TH1D("hsystCorr", "hsystCorr", 100, 0, 10);
//   hstat     ->Sumw2(); 
//   hsyst     ->Sumw2();
//   hsystCorr ->Sumw2(); 
//   Int_t nbinx = hstat->GetNbinsX();
  
//   for(Int_t ibinx = 1; ibinx <= nbinx; ibinx++){
//     hsyst->SetBinContent (ibinx, hstat->GetBinContent(ibinx));
//     hsyst->SetBinError   (ibinx, hstat->GetBinContent(ibinx)*0.08);
//     hsystCorr->SetBinContent (ibinx, hstat->GetBinContent(ibinx));
//     hsystCorr->SetBinError   (ibinx, hstat->GetBinContent(ibinx)*0.15);
//   } 

//   *SaveAs_Title += "pT_Spectrum_DirtySignal_postAnalyserCuts";
//   texXtitle = texPtX;
//   texYtitle = texptDifferentialYield;

// }