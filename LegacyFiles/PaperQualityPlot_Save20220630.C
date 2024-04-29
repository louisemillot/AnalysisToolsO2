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

#include<array>
using namespace std;

void SetStyle(Bool_t graypalette=kFALSE);
void myPadSetUp(TPad *currentPad, float currentLeft=0.11, float currentTop=0.04, float currentRight=0.04, float currentBottom=0.15);
void DrawLogo (Int_t logo=0, Double_t xmin =  0.28, Double_t ymin= 0.68) ;
void FakeHistosOnlyForExample(TH1*&hstat, TH1*&hsyst, TH1*&hsystCorr);
void LoadLibs();

void test_SetDefaultSumw2();

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
void Efficiency_TpcCrossedRowsScan(TH1D* &hstat, TH1D* &hsyst, TH1D*&hsystCorr, Int_t ipart, TFile* file_O2Analysis, TString* SaveAs_Title, TString* &texXtitle, TString* &texYtitle, Double_t* pTbins, Int_t nbinpT, Double_t minTpcCrossedRowscut, Double_t maxTpcCrossedRowscut);

void Systematics_CutVariations_Graphs(TH1D* &hstat, TH1D* &hsyst, TH1D* &hsystCorr, Int_t ipart, TFile** file_O2Analysis_CutVariation_array, Int_t id_referenceAnalysis, Int_t nCutsIterations, TString* SaveAs_Title, TString* &texXtitle, TString* &texYtitle, Double_t* pTbins, Int_t nbinpT);
void CorrectedYield_withSystematics(TH1D* &hstat, TH1D* &hsyst, TH1D* &hsystCorr, Int_t ipart, TFile** file_O2Analysis_CutVariation_array, Int_t id_referenceAnalysis, Int_t nCutsIterations, TString* SaveAs_Title, TString* &texXtitle, TString* &texYtitle, Double_t* pTbins, Int_t nbinpT);
void Draw_FeedDownMatrix(TH1D* &hstat, TH1D* &hsyst, TH1D* &hsystCorr, Int_t ipart, TFile *file_O2Analysis, TString* SaveAs_Title, TString* &texXtitle, TString* &texYtitle, Double_t* pTbins, Int_t nbinpT);
void Draw_XsiPlus_MassPlot(TH1D* &hstat, TH1D* &hsyst, TH1D* &hsystCorr, Int_t ipart, TFile *file_O2Analysis, TString* SaveAs_Title, TString* &texXtitle, TString* &texYtitle, Double_t* pTbins, Int_t nbinpT);
void Draw_FeedDownMatrix(TH1D* &hstat, TH1D* &hsyst, TH1D* &hsystCorr, Int_t ipart, TFile *file_O2Analysis, TString* SaveAs_Title, TString* &texXtitle, TString* &texYtitle, Double_t* pTbins, Int_t nbinpT);
void Draw_FeedDown(TH1D* &hstat, TH1D* &hsyst, TH1D* &hsystCorr, Int_t ipart, TFile *file_O2Analysis, TString* SaveAs_Title, TString* &texXtitle, TString* &texYtitle, Double_t* pTbins, Int_t nbinpT);
void Draw_TH1_Histograms(TH1** histograms_collection, TString** pdfLegend_collection, Int_t collectionSize, TString* pdfName, TString* &texXtitle, TString* &texYtitle, Double_t* bins, Int_t nbins, Int_t ratioOption);
void Draw_TH1_Histograms_withOneHistSystematics(TH1** histograms_collection, TString** pdfLegend_collection, Int_t collectionSize, TString* pdfName, TString* &texXtitle, TString* &texYtitle, Double_t* bins, Int_t nbins, Int_t ratioOption);
void Draw_TH1_Histograms_withOneHistSystematics_ratio(TH1** histograms_collection, TString** pdfLegend_collection, Int_t collectionSize, TString* pdfName, TString* &texXtitle, TString* &texYtitle, Double_t* bins, Int_t nbins, Int_t ratioOption);
void Test_LambdaTrueMCYieldvsHEP(TH1D* &hstat, TH1D* &hsyst, TH1D* &hsystCorr, Int_t ipart, TFile* file_O2Analysis, Int_t nCutsIterations, TString* SaveAs_Title, TString* &texXtitle, TString* &texYtitle, Double_t* pTbins, Int_t nbinpT);
void Test_V0DaugtherPairsTrackingEfficiency(TH1D* &hstat, TH1D* &hsyst, TH1D* &hsystCorr, Int_t ipart, TFile* file_O2Analysis, TString* SaveAs_Title, TString* &texXtitle, TString* &texYtitle, Double_t* pTbins, Int_t nbinpT);
void Test_V0DaugtherPairsTrackingEfficiency_PostDcaFitter(TH1D* &hstat, TH1D* &hsyst, TH1D* &hsystCorr, Int_t ipart, TFile* file_O2Analysis, TString* SaveAs_Title, TString* &texXtitle, TString* &texYtitle, Double_t* pTbins, Int_t nbinpT);

bool reject = false ;
float fparab(float *x, float *par);
float fline(float *x, float *par);
Double_t myLevyPtXsi(Double_t *pt, Double_t *par);
// Histogram processing functions; 

// returns pt diff raw spectrum in hPtDiff_RawSpectrum_O2
// void PtDifferential_Fit(TH1D* &hPtDiff_RawSpectrum_O2, TH3D* H3D_DetectedV0s_O2, TH3D* H3D_DetectedV0s_O2_rebinnedZ, Int_t ipart, Double_t* pTbins, Int_t nbinpT, Double_t* parGaussianParab_Sigma);
void PtDifferential_Fit(TH1D* &hPtDiff_RawSpectrum_O2, TH3D* H3D_DetectedV0s_O2, TH3D* H3D_DetectedV0s_O2_rebinnedZ, Int_t ipart, Double_t* pTbins, Int_t nbinpT, Double_t* parGaussianParab_Sigma, Double_t* parGaussianParab_Sigma_error);

//Custom Get functions
void Get_McEfficiency_vsPt(TH1D* &hMcEfficiency_vsPt, TH3D* H3D_detectedV0s_TruePt, TH1D* TrueV0PtSpectrum, Int_t ipart, Double_t* pTbins, Int_t nbinpT, Int_t ibinXaxisCut_low, Int_t ibinXaxisCut_high);
void Get_McSignalBackgroundRatio_vsPt(TH1D* &hSignalBackgroundRatio_vsPt, TH3D* H3D_detectedV0s_DataLike, Int_t ipart, Double_t* pTbins, Int_t nbinpT, Int_t ibinXaxisCut_low, Int_t ibinXaxisCut_high);
void Get_Systematics_OneCut(TH1D** hSystematicUncertainty, TH1D** hSystematicUncertainty_PreBarlow, Int_t ipart, TFile **file_O2Analysis_CutVariation_array, Int_t id_referenceAnalysis, Int_t nCutsIterations, Double_t* pTbins, Int_t nbinpT, Int_t InvMassRebinFactor);
void Get_Systematics_OneCut_McTrueAndMcDatalikeDifferentCuts(TH1D* hSystematicUncertainty, TH1D* hSystematicUncertainty_PreBarlow, Int_t ipart, TFile **file_O2Analysis_CutVariation_Datalike_array, TFile **file_O2Analysis_CutVariation_True_array, Int_t id_referenceAnalysis, Int_t nCutsIterations, Double_t* pTbins, Int_t nbinpT, Int_t InvMassRebinFactor);
void Get_Systematics_Fit_SignalExtractionType(TH1D* hSystematicUncertainty, TH1D* hSystematicUncertainty_PreBarlow, Int_t ipart, TH3D* H3D_detectedV0s_DataLike, TH1I* H1I_SelectedEventCount, TH3D* H3D_detectedV0s_TruePt, TH1D* TrueV0PtSpectrum, Double_t* pTbins, Int_t nbinpT, Int_t InvMassRebinFactor);
void Get_Systematics_Fit_SideBandVariation(TH1D* hSystematicUncertainty, TH1D* hSystematicUncertainty_PreBarlow, Int_t ipart, TH3D* H3D_detectedV0s_DataLike, TH1I* H1I_SelectedEventCount, TH3D* H3D_detectedV0s_TruePt, TH1D* TrueV0PtSpectrum, Double_t* pTbins, Int_t nbinpT, Int_t InvMassRebinFactor);
void Get_Systematics_Fit_Binning(TH1D* hSystematicUncertainty, TH1D* hSystematicUncertainty_PreBarlow, Int_t ipart, TH3D* H3D_detectedV0s_DataLike, TH1I* H1I_SelectedEventCount, TH3D* H3D_detectedV0s_TruePt, TH1D* TrueV0PtSpectrum, Double_t* pTbins, Int_t nbinpT);
void Get_RawYield_vsPt(TH1D* &hRawYield_vsPt, TH3D* H3D_detectedV0s_DataLike, TH1I* H1I_SelectedEventCount, Int_t ipart, Double_t* pTbins, Int_t nbinpT, Int_t ibinXaxisCut_low,  Int_t ibinXaxisCut_high, Int_t warning_cutArry_ID, Double_t SideBandSizeMultiplierModifier, Int_t SignalExtractionType, Int_t InvMassRebinFactor);
void Get_FeedDown(TH1D * H1D_FeedDownCorrection, Int_t ipart, TFile *file_O2Analysis, Double_t* pTbins, Int_t nbinpT);


std::array<Double_t, 3> GetTotalAndBackgroundCount_SidebandBackgroundExtrapolation(TH1D* H1D_detectedV0s_DataLike_SmallpTinterval, Double_t SignalMean, Double_t SignalStandardDev, Double_t SideBandSizeMultiplierModifier);
std::array<Double_t, 3> GetTotalAndBackgroundCount_BackgroundIntegration(TH1D* H1D_detectedV0s_DataLike_SmallpTinterval, Double_t SignalMean, Double_t SignalStandardDev, Double_t SideBandSizeMultiplierModifier, TFitResultPtr FitResult, TF1 *bkgparab);

// Preferred colors and markers
const Int_t fillColors[] = {kGray+1,  kRed-10, kBlue-9, kGreen-8, kMagenta-9, kOrange-9,kCyan-8,kYellow-7}; // for syst bands
const Int_t colors[]     = {kBlack, kRed+1 , kBlue+1, kGreen+3, kMagenta+1, kOrange-1,kCyan+2,kYellow+2};
const Int_t markers[]    = {kFullCircle, kFullSquare,kOpenCircle,kOpenSquare,kOpenDiamond,kOpenCross,kFullCross,kFullDiamond,kFullStar,kOpenStar};

//Options:
// TFile* file_O2Analysis = new TFile("AnalysisResults_run505673.root");
// TFile* file_O2Analysis = new TFile("AnalysisResults_run505600.root");
TFile* file_O2Analysis = new TFile("AnalysisResult_TpcXrows_Ref.root");
// TFile* file_O2Analysis = new TFile("AnalysisResult_testTrackEfficiency.root");
// TFile* file_O2Analysis = new TFile("AnalysisResults_largeK0Speak.root");
// TFile* file_O2Analysis = new TFile("AnalysisResults_thinK0Speak.root");

// TFile* file_O2Analysis = new TFile("AnalysisResults_noMC.root");
// TFile* file_O2Analysis = new TFile("AnalysisResultsMC.root");
const Int_t numPart = 3;
const Int_t ipart = 0;

const  Int_t nCutsIterations = 3;
// TFile* file_O2Analysis1 = new TFile("AnalysisResults_run505673.root");
// TFile* file_O2Analysis2 = new TFile("AnalysisResults_run505673.root");
// TFile* file_O2Analysis1 = new TFile("AnalysisResult_testTrackEfficiency.root");
// TFile* file_O2Analysis2 = new TFile("AnalysisResult_testTrackEfficiency.root");
TFile* file_O2Analysis1 = new TFile("AnalysisResult_TpcXrows_Ref.root");
TFile* file_O2Analysis2 = new TFile("AnalysisResult_TpcXrows_Tight.root");
TFile* file_O2Analysis3 = new TFile("AnalysisResult_TpcXrows_veryTight.root");
// TFile* file_O2Analysis_CutVariation_array[nCutsIterations] = {file_O2Analysis1, file_O2Analysis2};
TFile* file_O2Analysis_CutVariation_array[nCutsIterations] = {file_O2Analysis1, file_O2Analysis2, file_O2Analysis3};
const  Int_t id_referenceAnalysis = 0;
const Int_t N_SigmaBarlow = 2; //often 2 in ALICE LF analyses
const Int_t N_SideBandSizeVariation = 3;
const Int_t id_referenceSideband = 0.;
const Float_t SideBandSizeMultiplierModifier_array[N_SideBandSizeVariation] = {0, -0.5, -1}; // can't go above 1: 3+1 is 4 and that's the higher we can go with lambda having 4MeV sigma in rare cases and 2*4*4MeV going below the lower cutoff of pi+p masses; can't go below 0 because then we start losing signal that's within the guassian area; +/-3Sigma is 99.9% of gaussian area 
const Int_t SignalExtractionType_default = 1; //default is 1 for integration as opposed to 0 for bin counting; if using bin counting (0) then make sure H3D inv mass window is large enough to accomodate the 12sigma on both sides!
const Int_t InvMassRebinFactor_standard[numPart] = {10,5,5}; // I like 10 best for K0s at InvMass 0.4 to 0.6 with 400 bins, do variations around that
const Int_t N_BinningVariation = 5;
const Int_t InvMassRebinFactor_variation_array[numPart][N_BinningVariation] = {{5, 10, 2, 1, 20}, {5, 2, 1, 10, 4}, {5, 2, 1, 10, 4}};

const TString NamehistoInvMass[numPart] = {"InvMassK0S", "InvMassLambda", "InvMassAntiLambda"};//,"InvMassXiPlus", "InvMassXiMinus", "InvMassOmegaPlus", "InvMassOmegaMinus"};
const TString NamePart[numPart] = {"K0Short", "Lambda", "AntiLambda"};//, "XiPlu", "XiMin", "OmPlu", "OmMin"};
const TString NamePart_multiStrange[2] = {"XsiMinus", "XsiPlus"};
const TString NamePart_Latex[numPart] = {"K^{0}_{S}", "#Lambda", "#bar{#Lambda}"};//, "XiPlu", "XiMin", "OmPlu", "OmMin"};
// const Float_t LowerLimitDisplay_Mass[numPart] = {0.455, 1.085, 1.085};
// const Float_t UpperLimitDisplay_Mass[numPart] = {0.525, 1.145, 1.145};
const Float_t LowerLimitDisplay_Mass[numPart] = {0.4, 1.082, 1.082};
const Float_t UpperLimitDisplay_Mass[numPart] = {0.6, 1.146, 1.146};

const Float_t MassPart[numPart] = {0.497611, 1.115683, 1.115683};

const Float_t min_range_signal[numPart] = {0.47, 1.105, 1.105};
const Float_t max_range_signal[numPart] = {0.51, 1.1205, 1.1205};
const Float_t liminf[numPart] = {0.45, 1.08, 1.08};
const Float_t limsup[numPart] = {0.55, 1.18, 1.18};

// pT binning for K0S, Lambda and Antilambda
Double_t pTbins[numPart][20] = {{0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.6, 2.0, 2.4, 3.0},{0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 2.0, 2.4, 3.0},{0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 2.0, 2.4, 3.0}};
Int_t nbinpT[numPart] = {16,8,8};
// Double_t pTbins[numPart][20] = {{0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.6, 2.0, 2.4, 3.0},{0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 2.0, 2.4, 3.0},{0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 2.0, 2.4, 3.0}};
// Int_t nbinpT[numPart] = {16,8,8};
// Double_t pTbins[numPart][20] = {{0.0, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.6, 2.0, 2.4, 3.0},{0.0, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 2.0, 2.4, 3.0, 3.5},{0.0, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 2.0, 2.4, 3.0, 3.5}};
// Int_t nbinpT[numPart] = {17,10,10};
// Double_t pTbins[numPart][20] = {{0.0, 0.2, 0.4, 0.6, 0.8, 1., 1.2, 1.6, 2., 2.4, 3.0},{0.0, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 2.0, 2.4, 3.5},{0.0, 0.6, 0.8, 1.2, 1.6, 2.0, 2.4, 3.5}};
// Int_t nbinpT[numPart] = {10,9,7};
// Double_t pTbins[numPart][20] = {{0.0, 0.4, 0.8, 1.2, 1.6, 3.0},{0.0, 0.6, 1.0, 1.4, 2.0, 3.5},{0.0, 0.6, 1.0, 1.4, 2.0, 3.5}};
// Int_t nbinpT[numPart] = {5,5,5};
// Double_t pTbins[numPart][20] = {{0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.6, 2.0, 2.4, 3.0},{0.0, 0.6, 1.0, 1.4, 2.0, 3.0},{0.0, 0.6, 1.0, 1.4, 2.0, 3.0}};
// Int_t nbinpT[numPart] = {16,5,5};
// Double_t pTbins[numPart][4] = {{0., 10.},{0., 3.},{0., 3.}};
// Int_t nbinpT[numPart] = {1,1,1};


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
TString* texEfficiency = new TString("Efficiency #it{V0}_{detected}/#it{V0}_{MC} (ratio)");
TString* texDaughterRecoEfficiency = new TString("Track Reco Efficiency #it{N}_{reco}/#it{N}_{genMC} (ratio)");
TString* texDaughterPairsRecoEfficiency = new TString("Track Reco Efficiency V0Pairs (ratio)");

TString* texK0S_daughter_PiPlus = new TString("#pi^{+}");
TString* texLambda_daughter_PiMinus = new TString("#pi^{-}");
TString* texAntiLambda_daughter_PiPlus = new TString("#pi^{+}");

TString* texK0S_daughter_PiMinus = new TString("#pi^{-}");
TString* texLambda_daughter_protonPlus = new TString("p");
TString* texAntiLambda_daughter_protonMinus = new TString("#bar{p}");

TString* texRatio = new TString("ratio");

const TString* texInvMass_K0ShortDecay = new TString("M_{#pi^{+}#pi^{-}} (GeV/#it{c}^{2})");
const TString* texInvMass_LambdaDecay = new TString("M_{p#pi^{-}} (GeV/#it{c}^{2})");
const TString* texInvMass_AntiLambdaDecay = new TString("M_{#bar{p}#pi^{+}} (GeV/#it{c}^{2})");

const TString* texInvMassDecays_titles[numPart] = {texInvMass_K0ShortDecay,texInvMass_LambdaDecay,texInvMass_AntiLambdaDecay};

void PaperQualityPlot_Current() {
  // Load necessary libraries
  LoadLibs();
  // Set the default style
  SetStyle();

  // Automatically creates error structure when creating new histograms
  // TH1::SetDefaultSumw2(kTRUE); // unfortunately seems to remove error bars on my graphs; not sure why
  // test_SetDefaultSumw2();

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
  // InvMassDistributions_withFit_MC_datalike(hstat, hsyst, hsystCorr, ipart, file_O2Analysis, SaveAs_Title, texXtitle, texYtitle, pTbins[ipart], nbinpT[ipart]);


  // PseudoEfficiency_HEPcomparison_allCountSignalRegion(hstat, hsyst, hsystCorr, ipart, file_O2Analysis, SaveAs_Title, texXtitle, texYtitle, pTbins[ipart], nbinpT[ipart]);
  // PseudoEfficiency_HEPcomparison_allCountSignalRegionOLD_withBugs(hstat, hsyst, hsystCorr, ipart, file_O2Analysis, SaveAs_Title, pTbins[ipart], nbinpT[ipart]);
  // RawSpectrum_O2data_withFit(hstat, hsyst, hsystCorr, ipart, file_O2Analysis, SaveAs_Title, texXtitle, texYtitle, pTbins[ipart], nbinpT[ipart]);
  // Efficiency_O2MCdata_withFit(hstat, hsyst, hsystCorr, ipart, file_O2Analysis, SaveAs_Title, texXtitle, texYtitle, pTbins[ipart], nbinpT[ipart]);
  // Efficiency_O2data_allCountSignalRegion(hstat, hsyst, hsystCorr, ipart, file_O2Analysis, SaveAs_Title, texXtitle, texYtitle, pTbins[ipart], nbinpT[ipart]);
  // Efficiency_O2data_truePt_trueV0s(hstat, hsyst, hsystCorr, ipart, file_O2Analysis, SaveAs_Title, texXtitle, texYtitle, pTbins[ipart], nbinpT[ipart]);
  // RawSpectrum_O2data_truePt_trueV0s(hstat, hsyst, hsystCorr, ipart, file_O2Analysis, SaveAs_Title, texXtitle, texYtitle, pTbins[ipart], nbinpT[ipart]);
  // V0Daughters_TrackingEfficiencies_Daughter2(hstat, hsyst, hsystCorr, ipart, file_O2Analysis, SaveAs_Title, texXtitle, texYtitle, Extra, pTbinsProtons, nbinpTProtons, pTbinsPions, nbinpTPions);

  bool isMC = true;

  Double_t maxDCAcut = 1;
  Double_t minTpcCrossedRowscut = 0;
  Double_t maxTpcCrossedRowscut = 200;
  TH3D* H3D_detectedV0s_TruePt_DcaHist = (TH3D*)file_O2Analysis->Get("lambdakzero-analysis-mc/h3dMass"+NamePart[ipart]+"Dca");
  Double_t* DCAbins_fullFileRange = (Double_t*)H3D_detectedV0s_TruePt_DcaHist->GetXaxis();
  // Efficiency_DcaScan(hstat, hsyst, hsystCorr, ipart, file_O2Analysis, SaveAs_Title, texXtitle, texYtitle, pTbins[ipart], nbinpT[ipart], maxDCAcut);
  // Efficiency_TpcCrossedRowsScan(hstat, hsyst, hsystCorr, ipart, file_O2Analysis, SaveAs_Title, texXtitle, texYtitle, pTbins[ipart], nbinpT[ipart], minTpcCrossedRowscut, maxTpcCrossedRowscut);

  // Systematics_CutVariations_Graphs(hstat, hsyst, hsystCorr, ipart, file_O2Analysis_CutVariation_array, id_referenceAnalysis, nCutsIterations, SaveAs_Title, texXtitle, texYtitle, pTbins[ipart], nbinpT[ipart]);
  CorrectedYield_withSystematics(hstat, hsyst, hsystCorr, ipart, file_O2Analysis_CutVariation_array, id_referenceAnalysis, nCutsIterations, SaveAs_Title, texXtitle, texYtitle, pTbins[ipart], nbinpT[ipart]);
  // Draw_FeedDownMatrix(hstat, hsyst, hsystCorr, ipart, file_O2Analysis, SaveAs_Title, texXtitle, texYtitle, pTbins[ipart], nbinpT[ipart]);
  // Draw_XsiPlus_MassPlot(hstat, hsyst, hsystCorr, ipart, file_O2Analysis, SaveAs_Title, texXtitle, texYtitle, pTbins[ipart], nbinpT[ipart]);
  // Draw_FeedDown(hstat, hsyst, hsystCorr, ipart, file_O2Analysis, SaveAs_Title, texXtitle, texYtitle, pTbins[ipart], nbinpT[ipart]);
  // Test_LambdaTrueMCYieldvsHEP(hstat, hsyst, hsystCorr, ipart, file_O2Analysis, nCutsIterations, SaveAs_Title, texXtitle, texYtitle, pTbins[ipart], nbinpT[ipart]);

  // Test_V0DaugtherPairsTrackingEfficiency(hstat, hsyst, hsystCorr, ipart, file_O2Analysis, SaveAs_Title, texXtitle, texYtitle, pTbins[ipart], nbinpT[ipart]);
  // Test_V0DaugtherPairsTrackingEfficiency_PostDcaFitter(hstat, hsyst, hsystCorr, ipart, file_O2Analysis, SaveAs_Title, texXtitle, texYtitle, pTbins[ipart], nbinpT[ipart]);

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
  // TH1 * h = cfig->DrawFrame(0,0,3,0.017); // PtDifferential_SigmaOfFit
  // TH1 * h = cfig->DrawFrame(0,0,3,0.1); // PseudoEfficiency_HEPcomparison
  // TH1 * h = cfig->DrawFrame(0,0,3,0.06); // Efficiency_K0S
  // TH1 * h = cfig->DrawFrame(0,0,3,0.7); // Efficiency_K0S new MC
  // TH1 * h = cfig->DrawFrame(0,0,3,1); // Efficiency Daughter Pairs reco
  TH1 * h = cfig->DrawFrame(0,0,3,0.4); // CorrectedYield K0S
  // TH1 * h = cfig->DrawFrame(0,0,3,0.05); // CorrectedYield Lambda
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
  hsyst->Draw("E2,same"); //SHOULD BE PLOTTED EVENTUALLY
  hstat->Draw("E,same");
  hstat->SetMarkerStyle(markers[0]);
  // use the same color for markers and lines
  hstat->SetMarkerColor(colors [icolor]);
  hstat->SetLineColor  (colors [icolor]);

  // Draw the logo   
  //  0: Just "ALICE" (for final data), to be added only if ALICE does not appear otherwise (e.g. in the legend)
  //  >0: ALICE Preliminary
  // DrawLogo(1, 0.59, 0.81);

  // // You should always specify the colliding system
  // // NOTATION: pp, p-Pb, Pb-Pb. 
  // // Don't forget to use #sqrt{s_{NN}} for p-Pb and Pb-Pb
  // // You can change the position of this with
  // TLatex * textContext = new TLatex (0.18,0.82,"#bf{ALICE Run 3 Performance - unofficial}"); //BOLD
  // textContext->SetTextSize(0.05);
  // textContext->SetNDC(kTRUE); //remove if I want x,y in TLatex to be in the coordinate system of the histogram
  // textContext->Draw();
  // // TLatex * textColl = new TLatex (0.18,0.75,"pp #sqrt{#it{s}} = 900 GeV, #it{B_{z}} = 0.2 T");
  // TLatex * textColl = new TLatex (0.18,0.75,"pp #sqrt{#it{s}} = 900 GeV, MC anchored pilot beam 2021");
  // textColl->SetTextSize(0.04);
  // textColl->SetNDC(kTRUE); //remove if I want x,y in TLatex to be in the coordinate system of the histogram
  // textColl->Draw();
  // TLatex * text_part = new TLatex (0.18,0.70,NamePart_Latex[ipart]);
  // text_part->SetTextSize(0.04);
  // text_part->SetNDC(kTRUE); //remove if I want x,y in TLatex to be in the coordinate system of the histogram
  // text_part->Draw();
  // TLatex * text_extra = new TLatex (0.65,0.55,*Extra);
  // text_extra->SetNDC(kTRUE); //remove if I want x,y in TLatex to be in the coordinate system of the histogram
  // text_extra->Draw();
  // // TLatex * text2 = new TLatex (0.55,55,"V0A Multiplicity Classes (Pb-Side)");
  // // text2->SetTextSizePixels(24);
  // // text2->Draw();

  TLatex * textColl = new TLatex (0.18,0.82,"pp #sqrt{#it{s}} = 900 GeV, MC anchored pilot beam 2021");
  textColl->SetTextSize(0.04);
  textColl->SetNDC(kTRUE); //remove if I want x,y in TLatex to be in the coordinate system of the histogram
  textColl->Draw();
  TLatex * text_part = new TLatex (0.18,0.75,NamePart_Latex[ipart]);
  text_part->SetTextSize(0.04);
  text_part->SetNDC(kTRUE); //remove if I want x,y in TLatex to be in the coordinate system of the histogram
  text_part->Draw();

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

void test_SetDefaultSumw2() {
  if (TH1::GetDefaultSumw2()) {
    cout << "SetDefaultSumw2 outside the function has GetDefaultSumw2() true inside test_SetDefaultSumw2() function)"<< endl;
  }
  else {
    cout << "SetDefaultSumw2 outside the function has GetDefaultSumw2() false inside test_SetDefaultSumw2() function)"<< endl;
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
  hstat->Divide(hPtDiff_RawSpectrum_O2_rebinned,hPtDiff_RawSpectrum_HEP[ipart], 1., 1., "b");

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

  if (ipart>2) {
    cout << "ipart>2, HEP data downloaded and coded only goes to 2" << endl;
  } 
  for(int ibinPt = 1; ibinPt <= nbinpT; ibinPt++){
    hPtDiff_RawSpectrum_HEP[ipart]->SetBinError(ibinPt, hPtDiff_RawSpectrum_HEP_statErrors[ipart]->GetBinContent(ibinPt));
  }

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
  
  TFitResultPtr fFitResult[nbinpT];
  // Double_t parGaussParab[nbinpT][6]; //parab backround
  Double_t parGaussParab[nbinpT][5]; //linear background
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
    // hFrame[ibinPt] = canvasMassvsPt->DrawFrame(0.45,0,0.55,1.45*H1D_DetectedV0s_O2_SmallpTinterval[iPt]->GetMaximum());
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

    // // Parabola background; with low statistics produces bad stuff; better for PbPb
    // bkgparab[iPt] = new TF1("parab_"+NamePart[ipart]+Form("_Pt%i", ibinPt), fparab, liminf[ipart], limsup[ipart], 4);
    // bkgparab[iPt]->SetLineColor(kGreen);
    // bkgparab[iPt]->FixParameter(3, ipart);

    // Linear background; less issues with brackground fit going negative; better for pp
    bkgparab[iPt] = new TF1("line_"+NamePart[ipart]+Form("_Pt%i", ibinPt), fline, liminf[ipart], limsup[ipart], 3);
    bkgparab[iPt]->SetLineColor(kGreen);
    bkgparab[iPt]->FixParameter(2, ipart);

    GaussPlusPolynom[iPt] = new TF1("GaussPlusPolynom_"+NamePart[ipart]+Form("_Pt%i", ibinPt), "gaus(0)+pol1(3)", liminf[ipart], limsup[ipart]);
    GaussPlusPolynom[iPt]->SetLineColor(kRed); 
    GaussPlusPolynom[iPt]->SetParName(0, "norm");
    GaussPlusPolynom[iPt]->SetParName(1, "mean");
    GaussPlusPolynom[iPt]->SetParName(2, "sigma");
    GaussPlusPolynom[iPt]->SetParName(3, "p0");
    GaussPlusPolynom[iPt]->SetParName(4, "p1");
    // GaussPlusPolynom[iPt]->SetParName(5, "p2"); //parab background only

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


  //   // float xmin_InvMass = H1D_DetectedV0s_O2_SmallpTinterval[iPt]->GetXaxis()->GetXmin();
  //   // float xmax_InvMass = H1D_DetectedV0s_O2_SmallpTinterval[iPt]->GetXaxis()->GetXmax();
  //   cout << "xmin_InvMass= " << liminf[ipart] << ", xmax_InvMass= " << limsup[ipart] << endl;

  //   int ibin_liminf = H1D_DetectedV0s_O2_SmallpTinterval[iPt]->GetXaxis()->FindBin(liminf[ipart]);
  //   int ibin_limsup = H1D_DetectedV0s_O2_SmallpTinterval[iPt]->GetXaxis()->FindBin(limsup[ipart]);
  //   int nbinInvMass = ibin_limsup - ibin_liminf;
  //   cout << "nbinInvMass = " << nbinInvMass << endl;

  //   int totalCount = H1D_DetectedV0s_O2_SmallpTinterval[iPt]->Integral(ibin_liminf+1,ibin_limsup-1); //integral of a hist counts the number of bins, does not care for bin width: not true integral
  //   double backgroundCount = bkgparab[iPt]->Integral(liminf[ipart],limsup[ipart])*1./((limsup[ipart]-liminf[ipart])/nbinInvMass); //bkgparab is not an histogram but a function; here we need to convert it to a count
  //   //issue I actually need the integral of the positive part of the background function; for now kept simple integral as background is close to 0 compared to signal
  //   // float totalCount = H1D_DetectedV0s_O2_SmallpTinterval[iPt]->Integral(1,-1); //integral of a hist counts the number of bins, does not care for bin width: not true integral
  //   // float backgroundCount = bkgparab[iPt]->Integral(liminf[ipart],limsup[ipart])*1./((limsup[ipart]-liminf[ipart])/nbinInvMass); //bkgparab is not an histogram but a function; here we need to convert it to a count
  //   // //issue I actually need the integral of the positive part of the background function; for now kept simple integral as background is close to 0 compared to signal



  //   double SignalCountFit = totalCount - backgroundCount;
  //   cout << "totalCount= " << totalCount << ", backgroundCount= " << backgroundCount << endl;
  //   cout << "SignalCountFit= " << SignalCountFit << endl;

  //   double dpT = hPtDiff_RawSpectrum_O2->GetXaxis()->GetBinWidth(ibinPt);
  //   float drapidity = 1.5;//1.5
  //   double d2N_dpTdy = SignalCountFit *1./dpT*1./drapidity;
  //   hPtDiff_RawSpectrum_O2->SetBinContent(ibinPt,d2N_dpTdy);
  //   hPtDiff_RawSpectrum_O2->SetBinError(ibinPt,sqrt(SignalCountFit*1./dpT*1./drapidity)); //Issue with errors; this error of the O2 data spectrum is false, not right formula with fit
  //   cout << "(ibinPt,d2N_dpTdy) = (" << ibinPt << "," << d2N_dpTdy*1./SelectedEventCount << ")" << endl;

  //   hPtDiff_RawSpectrum_HEP[ipart]->SetBinError(ibinPt,hPtDiff_RawSpectrum_HEP_statErrors[ipart]->GetBinContent(ibinPt));
  // }

    //////////////////////////////////////////////////////////////////////////////////
    //////////////////////////// Signal Extraction ///////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////

    // //OLD signal extraction, my naive one
    // cout << "xmin_InvMass= " << liminf[ipart] << ", xmax_InvMass= " << limsup[ipart] << endl;
    // int ibin_liminf = H1D_detectedV0s_DataLike_SmallpTinterval[iPt]->GetXaxis()->FindBin(liminf[ipart]);
    // int ibin_limsup = H1D_detectedV0s_DataLike_SmallpTinterval[iPt]->GetXaxis()->FindBin(limsup[ipart]);
    // int nbinInvMass = ibin_limsup - ibin_liminf;
    // cout << "nbinInvMass = " << nbinInvMass << endl;
    // int totalCount = H1D_detectedV0s_DataLike_SmallpTinterval[iPt]->Integral(ibin_liminf+1,ibin_limsup-1); //integral of a hist counts the number of bins, does not care for bin width: not true integral
    // double backgroundCount = bkgparab[iPt]->Integral(liminf[ipart],limsup[ipart])*1./((limsup[ipart]-liminf[ipart])/nbinInvMass); //bkgparab is not an histogram but a function; here we need to convert it to a count
    // float SignalCountFit = totalCount - backgroundCount;
    // if (backgroundCount < 0) {
    //   cout << "WARNING - negative background count from fit - (cutArry_ID,ibinPt) = (" << ibinPt << "," << warning_cutArry_ID << ")" << endl;
    // }
    //NEW signal extraction, based on pp paper from David:

    float SignalMean = parGaussParab[iPt][1];
    float SignalStandardDev = parGaussParab[iPt][2];
    cout << "(SignalMean, SignalStandardDev) = (" << SignalMean << ", " << SignalStandardDev << ")" << endl;
    std::array<Double_t, 3> Signal_vector;
    Float_t SideBandSizeMultiplierModifier = SideBandSizeMultiplierModifier_array[0];

    Int_t SignalExtractionType = SignalExtractionType_default;
    if (SignalExtractionType == 0) {
      Signal_vector = GetTotalAndBackgroundCount_SidebandBackgroundExtrapolation(H1D_DetectedV0s_O2_SmallpTinterval[iPt], SignalMean, SignalStandardDev, SideBandSizeMultiplierModifier);
    }
    if (SignalExtractionType == 1) {
      Signal_vector = GetTotalAndBackgroundCount_BackgroundIntegration(H1D_DetectedV0s_O2_SmallpTinterval[iPt], SignalMean, SignalStandardDev, SideBandSizeMultiplierModifier, fFitResult[iPt], bkgparab[iPt]);
    }

    Double_t SignalCount = Signal_vector[0];
    Double_t SignalCount_StatError = Signal_vector[1];
    if (SignalCount < 0) {
      cout << "WARNING - negative SignalCount count from fit - ibinPt = " << ibinPt << endl;
    }
    // cout << "totalCount= " << totalCount << ", backgroundCount= " << backgroundCount << endl;
    // cout << "SignalCount= " << SignalCount << endl;

    Double_t dpT = hPtDiff_RawSpectrum_O2->GetXaxis()->GetBinWidth(ibinPt);
    Double_t drapidity = 1.5; //1.5
    Double_t d2N_dpTdy = SignalCount *1./dpT*1./drapidity;
    hPtDiff_RawSpectrum_O2->SetBinContent(ibinPt,d2N_dpTdy);
    hPtDiff_RawSpectrum_O2->SetBinError(ibinPt,sqrt(SignalCount_StatError*1./dpT*1./drapidity));  // error on d2N_dpTdy
    
    cout << "(ibinPt, SignalCount, TotalCount) = (" << ibinPt << ", " << SignalCount << ", " << Signal_vector[2] << ")" << endl;
  }

  canvasMassvsPt->SaveAs("Fits_invMass_ptDiff_"+NamePart[ipart]+".pdf");

  hPtDiff_RawSpectrum_O2->Scale(1./SelectedEventCount);
  hstat = (TH1D*)hPtDiff_RawSpectrum_O2->Clone("hstat");
  hstat->Reset("M");
  hstat->Divide(hPtDiff_RawSpectrum_O2,hPtDiff_RawSpectrum_HEP[ipart], 1., 1., "b");

  //////Error Bars///////
  hsyst     = new TH1D("hsyst", "hsyst", nbinpT, pTbins);
  hsystCorr = new TH1D("hsystCorr", "hsystCorr", nbinpT, pTbins);
  hstat     ->Sumw2(); 
  hsyst     ->Sumw2();
  hsystCorr ->Sumw2(); 
  Int_t nbinx = hstat->GetNbinsX();
  
  for(Int_t ibinx = 1; ibinx <= nbinx; ibinx++){
    hsyst->SetBinContent (ibinx, hstat->GetBinContent(ibinx));
    hsyst->SetBinError   (ibinx, hstat->GetBinContent(ibinx)*0.000001);
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
    // TLatex * textColl = new TLatex (0.18,0.75,"pp #sqrt{#it{s}} = 0.9 TeV, pilot beam 2021");
    TLatex * textColl = new TLatex (0.18,0.75,"pp #sqrt{#it{s}} = 0.9 TeV, Monte Carlo");
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
  TH3D* H3D_DetectedV0s_O2_rebinnedZ = (TH3D*)H3D_DetectedV0s_O2->RebinZ(5,"H3D_DetectedV0s_O2_rebinnedZ");

  //Fit tools initialisation
  TF1 *gauss[nbinpT];
  TF1 *GaussPlusPolynom[nbinpT];
  TF1 *bkgparab[nbinpT];
  
  Float_t fFitResult[nbinpT];
  // Double_t parGaussParab[nbinpT][6]; //parabola background
  Double_t parGaussParab[nbinpT][5]; //line background
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

    // GaussPlusPolynom[iPt] = new TF1("GaussPlusPolynom_"+NamePart[ipart]+Form("_Pt%i", ibinPt), "gaus(0)+pol2(3)", liminf[ipart], limsup[ipart]); //parabola background
    GaussPlusPolynom[iPt] = new TF1("GaussPlusPolynom_"+NamePart[ipart]+Form("_Pt%i", ibinPt), "gaus(0)+pol1(3)", liminf[ipart], limsup[ipart]); //line background
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

    fFitResult[iPt] = H1D_DetectedV0s_O2_SmallpTinterval[iPt]->Fit(GaussPlusPolynom[iPt],"S0RB+QL"); //L option is log likelyhood: would supposedly be better in case of counts

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
    bkgparab[iPt]->Draw("same");
    // gauss[iPt]->Draw("same");
    GaussPlusPolynom[iPt]->Draw("same");

    // TLatex * text_pTbin = new TLatex (0.62,0.75,Form("pT [%.1f;%.1f] GeV", pTbins[iPt], pTbins[iPt+1]));
    // text_pTbin->SetNDC(kTRUE);
    // text_pTbin->Draw();
    TLatex * textContext = new TLatex (0.18,0.82,"#bf{ALICE Performance}"); //BOLD
    textContext->SetTextSize(0.05);
    textContext->SetNDC(kTRUE); //remove if I want x,y in TLatex to be in the coordinate system of the histogram
    textContext->Draw();
    // TLatex * textColl = new TLatex (0.18,0.75,"pp #sqrt{#it{s}} = 0.9 TeV, pilot beam 2021");
    TLatex * textColl = new TLatex (0.18,0.75,"pp #sqrt{#it{s}} = 0.9 TeV, Monte Carlo");
    textColl->SetTextSize(0.04);
    textColl->SetNDC(kTRUE); //remove if I want x,y in TLatex to be in the coordinate system of the histogram
    textColl->Draw();
    TLatex * text_part = new TLatex (0.18,0.70,NamePart_Latex[ipart]+Form(",   %.1f < #it{p}_{T} < %.1f GeV/#it{c}", pTbins[iPt], pTbins[iPt+1])); //temporary
    // TLatex * text_part = new TLatex (0.18,0.70,NamePart_Latex[ipart]+",   0 < #it{p}_{T} < 10 GeV/#it{c}"); for performance figure
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
  hstat->Divide(hPtDiff_RawSpectrum_O2,TrueV0PtSpectrum_rebinned, 1., 1., "b");


  //////Error Bars///////
  hsyst     = (TH1D*)hstat->Clone("hsyst");
  hsyst->Reset("M");
  hsystCorr     = (TH1D*)hstat->Clone("hsystCorr");
  hsystCorr->Reset("M");
  hstat     ->Sumw2(); 
  hsyst     ->Sumw2();
  hsystCorr ->Sumw2(); 
  Int_t nbinx = hstat->GetNbinsX();
  
  for(Int_t ibinx = 1; ibinx <= nbinx; ibinx++){
    hsyst->SetBinContent (ibinx, hstat->GetBinContent(ibinx));
    hsyst->SetBinError   (ibinx, hstat->GetBinContent(ibinx)*0.0000001);
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
  hstat->Divide(hPtDiff_RawSpectrum_O2,TrueV0PtSpectrum_rebinned, 1., 1., "b");

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
  hstat->Divide(hPtDiff_RawSpectrum_O2,TrueV0PtSpectrum_rebinned, 1., 1., "b");

  hsyst = (TH1D*)hstat->Clone("hsyst");
  hsyst->Reset("M");

  //////Error Bars///////
  hsystCorr = new TH1D("hsystCorr", "hsystCorr", 100, 0, 10);
  hstat     ->Sumw2(); 
  hsyst     ->Sumw2();
  hsystCorr ->Sumw2(); 
  Int_t nbinx = hstat->GetNbinsX();
  
  for(Int_t ibinx = 1; ibinx <= nbinx; ibinx++){
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
    H1D_Eff_K0S_pionPlus->Divide(H1D_K0S_rec_pionPlus,H1D_K0S_gen_pionPlus, 1., 1., "b");
    canvasV0DaughterTrackEff->cd(1);
    H1D_Eff_K0S_pionPlus->Draw("E");
    TH1D* H1D_Eff_K0S_pionNeg = (TH1D*)H1D_K0S_gen_pionNeg_rebin->Clone("H1D_Eff_K0S_pionNeg");
    H1D_Eff_K0S_pionNeg->Divide(H1D_K0S_rec_pionNeg,H1D_K0S_gen_pionNeg, 1., 1., "b");
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
    H1D_Eff_Lambda_pionNeg->Divide(H1D_Lambda_rec_pionNeg,H1D_Lambda_gen_pionNeg, 1., 1., "b");
    canvasV0DaughterTrackEff->cd(1);
    H1D_Eff_Lambda_pionNeg->Draw("E");
    TH1D* H1D_Eff_Lambda_protonPlus = (TH1D*)H1D_Lambda_gen_protonPlus_rebin->Clone("H1D_Eff_Lambda_protonPlus");
    H1D_Eff_Lambda_protonPlus->Divide(H1D_Lambda_rec_protonPlus,H1D_Lambda_gen_protonPlus, 1., 1., "b");
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
    H1D_Eff_AntiLambda_pionPlus->Divide(H1D_AntiLambda_rec_pionPlus,H1D_AntiLambda_gen_pionPlus, 1., 1., "b");
    canvasV0DaughterTrackEff->cd(1);
    H1D_Eff_AntiLambda_pionPlus->Draw("E");
    TH1D* H1D_Eff_AntiLambda_protonNeg = (TH1D*)H1D_AntiLambda_gen_protonNeg_rebin->Clone("H1D_Eff_AntiLambda_protonNeg");
    H1D_Eff_AntiLambda_protonNeg->Divide(H1D_AntiLambda_rec_protonNeg,H1D_AntiLambda_gen_protonNeg, 1., 1., "b");
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
    daughter1->Divide(H1D_K0S_rec_pionPlus_rebin,H1D_K0S_gen_pionPlus_rebin, 1., 1., "b");
    daughter2 = (TH1D*)H1D_K0S_gen_pionNeg_rebin->Clone("H1D_Eff_K0S_pionNeg");
    daughter2->Divide(H1D_K0S_rec_pionNeg_rebin,H1D_K0S_gen_pionNeg_rebin, 1., 1., "b"); 
  }

  if (ipart == 1) {
    daughter1 = (TH1D*)H1D_Lambda_gen_pionNeg_rebin->Clone("H1D_Eff_Lambda_pionNeg");
    daughter1->Divide(H1D_Lambda_rec_pionNeg_rebin,H1D_Lambda_gen_pionNeg_rebin, 1., 1., "b");
    daughter2 = (TH1D*)H1D_Lambda_gen_protonPlus_rebin->Clone("H1D_Eff_Lambda_protonPlus");
    daughter2->Divide(H1D_Lambda_rec_protonPlus_rebin,H1D_Lambda_gen_protonPlus_rebin, 1., 1., "b");
  }

  if (ipart == 2) {
    daughter1 = (TH1D*)H1D_AntiLambda_gen_pionPlus_rebin->Clone("H1D_Eff_AntiLambda_pionPlus");
    daughter1->Divide(H1D_AntiLambda_rec_pionPlus_rebin,H1D_AntiLambda_gen_pionPlus_rebin, 1., 1., "b");
    daughter2 = (TH1D*)H1D_AntiLambda_gen_protonNeg_rebin->Clone("H1D_Eff_AntiLambda_protonNeg");
    daughter2->Divide(H1D_AntiLambda_rec_protonNeg_rebin,H1D_AntiLambda_gen_protonNeg_rebin, 1., 1., "b");
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










//------------------------------- Newer functions more rigorous error wise -------------------------------//
void Get_McEfficiency_vsPt(TH1D* &hMcEfficiency_vsPt, TH3D* H3D_detectedV0s_TruePt, TH1D* TrueV0PtSpectrum, Int_t ipart, Double_t* pTbins, Int_t nbinpT, Int_t ibinXaxisCut_low, Int_t ibinXaxisCut_high) {

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

    H1D_detectedV0s_TruePt_SmallpTinterval[iPt] = (TH1D*)H3D_detectedV0s_TruePt->ProjectionZ("InvMass_"+NamePart[ipart]+Form("_Pt%i", ibinPt),ibinXaxisCut_low,ibinXaxisCut_high,ibinPt_originalAxis_low,ibinPt_originalAxis_up);

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
  hMcEfficiency_vsPt->Divide(hMcSignalCount_vsPt,TrueV0PtSpectrum_AnalysisBins, 1., 1., "b"); // option b for binomial because efficiency: https://twiki.cern.ch/twiki/bin/view/ALICE/PWGLFPAGSTRANGENESSEfficiency
  delete hMcSignalCount_vsPt; // avoids memory leaks when Get_McEfficiency_vsPt is called in a loop
  delete TrueV0PtSpectrum_AnalysisBins; // avoids memory leaks when Get_McEfficiency_vsPt is called in a loop
}

void Get_McSignalBackgroundRatio_vsPt(TH1D* &hSignalBackgroundRatio_vsPt, TH3D* H3D_detectedV0s_DataLike, Int_t ipart, Double_t* pTbins, Int_t nbinpT, Int_t ibinXaxisCut_low,  Int_t ibinXaxisCut_high) {
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

    H1D_detectedV0s_DataLike_SmallpTinterval[iPt] = (TH1D*)H3D_detectedV0s_DataLike->ProjectionZ("InvMass_"+NamePart[ipart]+Form("_Pt%i", ibinPt),ibinXaxisCut_low,ibinXaxisCut_high,ibinPt_originalAxis_low,ibinPt_originalAxis_up);

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

std::array<Double_t, 3> GetTotalAndBackgroundCount_SidebandBackgroundExtrapolation(TH1D* H1D_detectedV0s_DataLike_SmallpTinterval, Double_t SignalMean, Double_t SignalStandardDev, Double_t SideBandSizeMultiplierModifier) {
  // iwth our larger peaks due to lower pagnetic field, we can't go higher than 8sigma for the lower background bound for lambda (due to cutoff from pi+p masses), or 10sigma for k0s (with current window limit)
  // 8 sigma limit is 2*4, so we sample the signal in [-4sigma;+4sigma] interval
  // for a start, the sidebandModifier should be floats, ints is too large a change
  // could have 2*4, 2*3.5, 2*3? 3 is probably the lowest we can go lest we start missing lot of signal (99.9% of guassian area is within -3;+3 sigma)
  Double_t LeftBackgroundLowerBound = SignalMean - 2*(4 + SideBandSizeMultiplierModifier)*SignalStandardDev;
  Double_t LeftBackgroundUpperBound = SignalMean - (4 + SideBandSizeMultiplierModifier)*SignalStandardDev;
  Double_t RightBackgroundLowerBound = SignalMean + (4 + SideBandSizeMultiplierModifier)*SignalStandardDev;
  Double_t RightBackgroundUpperBound = SignalMean + 2*(4 + SideBandSizeMultiplierModifier)*SignalStandardDev;
  cout <<"Sidebands: (" << LeftBackgroundLowerBound << ", " << LeftBackgroundUpperBound << ", " << RightBackgroundLowerBound << ", " << RightBackgroundUpperBound << ")" << endl;

  Int_t LeftBackgroundLowerBound_bin = H1D_detectedV0s_DataLike_SmallpTinterval->GetXaxis()->FindBin(LeftBackgroundLowerBound);
  Int_t LeftBackgroundUpperBound_bin = H1D_detectedV0s_DataLike_SmallpTinterval->GetXaxis()->FindBin(LeftBackgroundUpperBound);
  Int_t RightBackgroundLowerBound_bin = H1D_detectedV0s_DataLike_SmallpTinterval->GetXaxis()->FindBin(RightBackgroundLowerBound);
  Int_t RightBackgroundUpperBound_bin = H1D_detectedV0s_DataLike_SmallpTinterval->GetXaxis()->FindBin(RightBackgroundUpperBound);

  Int_t TotalCountInSignalRegion = H1D_detectedV0s_DataLike_SmallpTinterval->Integral(LeftBackgroundUpperBound_bin, RightBackgroundLowerBound_bin);

  // samples background in the two sidebands chosen: here mean-12sigma to mean-6sigma and mean+6sigma to mean+12sigma for sideBandmodifier=0
  // only do a count of the histogram, not an actual intgral of the background fit function to avoid issues with background line function going below 0 values due to very low background
  Int_t BackgroundCountInSidebands = H1D_detectedV0s_DataLike_SmallpTinterval->Integral(LeftBackgroundLowerBound_bin, LeftBackgroundUpperBound_bin) + H1D_detectedV0s_DataLike_SmallpTinterval->Integral(RightBackgroundLowerBound_bin, RightBackgroundUpperBound_bin); // symmetric sampling windows: mean-12sigma to mean-6sigma and mean+6sigma to mean+12Sigma
  
  // linear background function and symmetric sampling windows and same total size for background sampling window (12sigma) and for signal window (12sigma) means it's equal
  // can integrate a linear function to check if one wants to
  Int_t BackgroundExtrapolationInSignalRegion = BackgroundCountInSidebands; 

  Double_t SignalCount = TotalCountInSignalRegion - BackgroundExtrapolationInSignalRegion;
  cout << "leftBackgroundCount = " << H1D_detectedV0s_DataLike_SmallpTinterval->Integral(LeftBackgroundLowerBound_bin, LeftBackgroundUpperBound_bin) << ", rightBackgroundCount = "<< H1D_detectedV0s_DataLike_SmallpTinterval->Integral(RightBackgroundLowerBound_bin, RightBackgroundUpperBound_bin) << endl;
  cout << "TotalCountInSignalRegion = " << TotalCountInSignalRegion << endl;
  // Stat uncertainty on signal count:
  // Sigma_signalCount**2 = Sigma_backgroundCount**2 + Sigma_totalCount**2 = sqrt(Sigma_backgroundCount)**2 + sqrt(Sigma_totalCount)**2
  // see https://twiki.cern.ch/twiki/bin/view/ALICE/PWGLFPAGSTRANGENESSSignalExtraction
  Double_t SignalCount_StatError = sqrt(TotalCountInSignalRegion + BackgroundCountInSidebands);

  std::array<Double_t, 3> Signal_vector = {SignalCount, SignalCount_StatError, (Double_t)TotalCountInSignalRegion};
  return Signal_vector;
}

std::array<Double_t, 3> GetTotalAndBackgroundCount_BackgroundIntegration(TH1D* H1D_detectedV0s_DataLike_SmallpTinterval, Double_t SignalMean, Double_t SignalStandardDev, Double_t SideBandSizeMultiplierModifier, TFitResultPtr FitResult, TF1 *bkgparab) {
  // iwth our larger peaks due to lower pagnetic field, we can't go higher than 8sigma for the lower background bound for lambda (due to cutoff from pi+p masses), or 10sigma for k0s (with current window limit)
  // 8 sigma limit is 2*4, so we sample the signal in [-4sigma;+4sigma] interval
  // for a start, the sidebandModifier should be floats, ints is too large a change
  // could have 2*4, 2*3.5, 2*3? 3 is probably the lowest we can go lest we start missing lot of signal (99.9% of guassian area is within -3;+3 sigma)
  Double_t LeftBackgroundLowerBound = SignalMean - 2*(4 + SideBandSizeMultiplierModifier)*SignalStandardDev;
  Double_t LeftBackgroundUpperBound = SignalMean - (4 + SideBandSizeMultiplierModifier)*SignalStandardDev;
  Double_t RightBackgroundLowerBound = SignalMean + (4 + SideBandSizeMultiplierModifier)*SignalStandardDev;
  Double_t RightBackgroundUpperBound = SignalMean + 2*(4 + SideBandSizeMultiplierModifier)*SignalStandardDev;

  Int_t LeftBackgroundLowerBound_bin = H1D_detectedV0s_DataLike_SmallpTinterval->GetXaxis()->FindBin(LeftBackgroundLowerBound);
  Int_t LeftBackgroundUpperBound_bin = H1D_detectedV0s_DataLike_SmallpTinterval->GetXaxis()->FindBin(LeftBackgroundUpperBound);
  Int_t RightBackgroundLowerBound_bin = H1D_detectedV0s_DataLike_SmallpTinterval->GetXaxis()->FindBin(RightBackgroundLowerBound);
  Int_t RightBackgroundUpperBound_bin = H1D_detectedV0s_DataLike_SmallpTinterval->GetXaxis()->FindBin(RightBackgroundUpperBound);

  //OLD signal extraction, my naive one
  Int_t nbinInvMass = RightBackgroundLowerBound_bin - LeftBackgroundUpperBound_bin;
  Int_t TotalCountInSignalRegion = H1D_detectedV0s_DataLike_SmallpTinterval->Integral(LeftBackgroundUpperBound_bin,RightBackgroundLowerBound_bin); //integral of a hist counts the number of bins, does not care for bin width: not true integral
  Double_t BackgroundCountInSignalRegion = bkgparab->Integral(LeftBackgroundUpperBound,RightBackgroundLowerBound)*1./((RightBackgroundLowerBound-LeftBackgroundUpperBound)/nbinInvMass); //bkgparab is not an histogram but a function; here we need to convert it to a count
  if (BackgroundCountInSignalRegion < 0) {
    // cout << "WARNING - negative background count from fit - (cutArry_ID,ibinPt) = (" << ibinPt << "," << warning_cutArry_ID << ")" << endl;
    cout << "WARNING - negative background count from fit" << endl;
  }
  Double_t TotalCountInSignalRegion_SystError = sqrt(TotalCountInSignalRegion*TotalCountInSignalRegion);

  ///////////// Error calculation for Background integration /////////////
  // (https://root.cern/doc/master/ErrorIntegral_8C.html)
  // estimated integral  and error analytically

  Double_t * param = bkgparab->GetParameters();
  TMatrixDSym covMatrix = FitResult->GetCovarianceMatrix();

  Double_t Sigma_integral_fromRootFunction = bkgparab->IntegralError(LeftBackgroundUpperBound,RightBackgroundLowerBound, param, covMatrix.GetMatrixArray());// have to include the scaling factor used to get a pTdifferential curve from simple background count
  Double_t BackgroundCountInSignalRegion_SystError = Sigma_integral_fromRootFunction * sqrt(1./((RightBackgroundLowerBound-LeftBackgroundUpperBound)/nbinInvMass));

  ///////////// Outputs /////////////
  Double_t SignalCount = TotalCountInSignalRegion - BackgroundCountInSignalRegion;
  Double_t SignalCount_StatError = sqrt(TotalCountInSignalRegion_SystError*TotalCountInSignalRegion_SystError + BackgroundCountInSignalRegion_SystError*BackgroundCountInSignalRegion_SystError); //assuming uncorrelated; worst case anyway; one is count one is integration of fit

  std::array<Double_t, 3> Signal_vector = {SignalCount, SignalCount_StatError, (Double_t)TotalCountInSignalRegion};

  ////// Analytical error by hand: gives about 1/2 of the error that IntegralError() gives
  // // double Integral  = a/2*B**2 + b*B - a/2*A**2 - b*A; for our line ax+b between A and B
  // Double_t dI_da = 1./2*(RightBackgroundLowerBound*RightBackgroundLowerBound - LeftBackgroundUpperBound*LeftBackgroundUpperBound); // partial derivative = B**2 - A**2 for our line
  // Double_t dI_db = RightBackgroundLowerBound - LeftBackgroundUpperBound; // partial derivative = B - A for our line
  // // estimated error with correlations
  // Double_t Sigma_Integral = std::sqrt(dI_da*dI_da * covMatrix(0,0) + dI_db*dI_db * covMatrix(1,1) + 2.* dI_da*dI_db * covMatrix(0,1));// have to include the scaling factor used to get a pTdifferential curve from simple background count
  
  // if ( std::fabs(Sigma_integral_fromRootFunction - Sigma_Integral) > 1.E-6*Sigma_Integral ) {
  //   std::cout << " ERROR: test failed : different analytical  integral - (RootAutocalc, PersoDerivation) = (" << Sigma_integral_fromRootFunction << ", " << Sigma_Integral << ") - CHECK demo file again" << std::endl;
  // }

  return Signal_vector;
}



void Get_RawYield_vsPt(TH1D* &hRawYield_vsPt, TH3D* H3D_detectedV0s_DataLike, TH1I* H1I_SelectedEventCount, Int_t ipart, Double_t* pTbins, Int_t nbinpT, Int_t ibinXaxisCut_low,  Int_t ibinXaxisCut_high, Int_t warning_cutArry_ID, Double_t SideBandSizeMultiplierModifier, Int_t SignalExtractionType, Int_t InvMassRebinFactor) {
  //Looks within the signal area of the fit, 6Sigma

  // what is fed as input: TH3D* H3D_detectedV0s_DataLike = (TH3D*)file_O2Analysis->Get("lambdakzero-analysis-mc/h3dMass"+NamePart[ipart]);
  TH3D* H3D_detectedV0s_DataLike_rebinnedZ = (TH3D*)H3D_detectedV0s_DataLike->RebinZ(InvMassRebinFactor,"H3D_detectedV0s_DataLike_rebinnedZ");
  cout << "---------------------------------------- InvMassRebinFactor = " << InvMassRebinFactor << endl;
  ////////////////////////////////// Fit initialisation //////////////////////////////////
  //Fit tools initialisation
  TF1 *gauss[nbinpT];
  TF1 *GaussPlusPolynom[nbinpT];
  TF1 *bkgparab[nbinpT];
  TFitResultPtr fFitResult[nbinpT];
  // Double_t parGaussParab[nbinpT][6]; //parab backround
  Double_t parGaussParab[nbinpT][5]; //linear background

  TH1D *H1D_detectedV0s_DataLike_SmallpTinterval[nbinpT];
  TH1D *hPtDiff_RawYield_temp = new TH1D("hPtDiff_RawYield_temp","hPtDiff_RawYield_temp",nbinpT,pTbins);

  float SelectedEventCount = H1I_SelectedEventCount->GetEntries();
  cout << "*--------------------- Get_RawYield_vsPt ---------------------*" << endl;

  TCanvas *canvasMassvsPtFits = new TCanvas ("canvasMassPt_fits", NamehistoInvMass[ipart], 800, 600);
  // TCanvas *canvasMassvsPt = new TCanvas ("canvasMassPt", NamehistoInvMass[ipart], 1600, 1200);
  canvasMassvsPtFits->Divide(std::ceil(std::sqrt(nbinpT)),std::ceil(std::sqrt(nbinpT)));
  TH1 *hFrameFits[nbinpT];
  cout << "testAAA" << endl;
  for(int ibinPt = 1; ibinPt <= nbinpT; ibinPt++){
    int iPt = ibinPt-1; //different convention for bins that start at 1 and arrays that start at 0
    
    canvasMassvsPtFits->cd(ibinPt); //ATTENTION: +1 because canvas divide starts counting from 1

    int ibinPt_originalAxis_low = H3D_detectedV0s_DataLike->GetYaxis()->FindBin(pTbins[iPt]);
    int ibinPt_originalAxis_up = H3D_detectedV0s_DataLike->GetYaxis()->FindBin(pTbins[iPt+1]) - 1;
    // cout << "ibinPt: " << ibinPt << "; originalAxis [low;up] = [" << ibinPt_originalAxis_low << ";" << ibinPt_originalAxis_up << "]" << endl;

    H1D_detectedV0s_DataLike_SmallpTinterval[iPt] = (TH1D*)H3D_detectedV0s_DataLike_rebinnedZ->ProjectionZ("InvMass_"+NamePart[ipart]+Form("_Pt%i", ibinPt),ibinXaxisCut_low,ibinXaxisCut_high,ibinPt_originalAxis_low,ibinPt_originalAxis_up);
    hFrameFits[iPt] = canvasMassvsPtFits->DrawFrame(LowerLimitDisplay_Mass[ipart],0,UpperLimitDisplay_Mass[ipart],1.45*H1D_detectedV0s_DataLike_SmallpTinterval[iPt]->GetMaximum());
    hFrameFits[iPt]->SetXTitle(texInvMassDecays_titles[ipart]->Data());
    hFrameFits[iPt]->SetYTitle("Counts");
    H1D_detectedV0s_DataLike_SmallpTinterval[iPt]->Draw("E,same");
    ////////////////////////////////////////////////////////////////////
    //////////////////////////// Fit start /////////////////////////////
    ////////////////////////////////////////////////////////////////////

    gauss[iPt] = new TF1("gauss_"+NamePart[ipart]+Form("_Pt%i", ibinPt), "gaus", min_range_signal[ipart], max_range_signal[ipart]);
    gauss[iPt]->SetLineColor(kAzure);   
    gauss[iPt]->SetParName(0, "norm");
    gauss[iPt]->SetParName(1, "mean");
    gauss[iPt]->SetParName(2, "sigma");
    gauss[iPt]->SetParameters(1., MassPart[ipart], 0.0002);
    gauss[iPt]->SetParLimits(0, 0., 1.1*H1D_detectedV0s_DataLike_SmallpTinterval[iPt]->GetBinContent(H1D_detectedV0s_DataLike_SmallpTinterval[iPt]->GetMaximumBin()));
    gauss[iPt]->SetParLimits(1, min_range_signal[ipart], max_range_signal[ipart]);
    gauss[iPt]->SetParLimits(2, 0.0001, 0.05);

    // cout << "par0 upper limit" <<  1.1*H1D_detectedV0s_DataLike_SmallpTinterval[iPt]->GetBinContent(H1D_detectedV0s_DataLike_SmallpTinterval[iPt]->GetMaximumBin()) << endl;
    // // Parabola background; with low statistics produces bad stuff; better for PbPb
    // bkgparab[iPt] = new TF1("parab_"+NamePart[ipart]+Form("_Pt%i", ibinPt), fparab, liminf[ipart], limsup[ipart], 4);
    // bkgparab[iPt]->SetLineColor(kGreen);
    // bkgparab[iPt]->FixParameter(3, ipart);

    // Linear background; less issues with brackground fit going negative; better for pp
    bkgparab[iPt] = new TF1("line_"+NamePart[ipart]+Form("_Pt%i", ibinPt), fline, liminf[ipart], limsup[ipart], 3);
    bkgparab[iPt]->SetLineColor(kGreen);
    bkgparab[iPt]->FixParameter(2, ipart);

    GaussPlusPolynom[iPt] = new TF1("GaussPlusPolynom_"+NamePart[ipart]+Form("_Pt%i", ibinPt), "gaus(0)+pol1(3)", liminf[ipart], limsup[ipart]);
    GaussPlusPolynom[iPt]->SetLineColor(kRed); 
    GaussPlusPolynom[iPt]->SetParName(0, "norm");
    GaussPlusPolynom[iPt]->SetParName(1, "mean");
    GaussPlusPolynom[iPt]->SetParName(2, "sigma");
    GaussPlusPolynom[iPt]->SetParName(3, "p0");
    GaussPlusPolynom[iPt]->SetParName(4, "p1");
    // GaussPlusPolynom[iPt]->SetParName(5, "p2"); //parab background only

    GaussPlusPolynom[iPt]->SetParameters(1., MassPart[ipart], 0.0002);
    GaussPlusPolynom[iPt]->SetParLimits(0, 0., 1.1*H1D_detectedV0s_DataLike_SmallpTinterval[iPt]->GetBinContent(H1D_detectedV0s_DataLike_SmallpTinterval[iPt]->GetMaximumBin()));
    GaussPlusPolynom[iPt]->SetParLimits(1, min_range_signal[ipart], max_range_signal[ipart]);
    GaussPlusPolynom[iPt]->SetParLimits(2, 0.0001, 0.05);
    // cout << "par0 upper limit" <<  1.1*H1D_detectedV0s_DataLike_SmallpTinterval[iPt]->GetBinContent(H1D_detectedV0s_DataLike_SmallpTinterval[iPt]->GetMaximumBin()) << endl;

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


    //////////////////////////////////////////////////////////////////////////////////
    //////////////////////////// Plot fits for checks /////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////

    //plot Inv Mass histogram and the background and signal fits for each iPt bin
    // gPad->DrawFrame(LowerLimitDisplay_Mass[ipart],0,UpperLimitDisplay_Mass[ipart],1.45*H1D_detectedV0s_DataLike_SmallpTinterval[iPt]->GetMaximum());
    // gPad->SetXTitle(texInvMassDecays_titles[ipart]->Data());
    // gPad->SetYTitle("Counts");
    // H1D_detectedV0s_DataLike_SmallpTinterval[iPt]->Draw("E,same");
    H1D_detectedV0s_DataLike_SmallpTinterval[iPt]->GetXaxis()->SetTitle(texInvMassDecays_titles[ipart]->Data());
    H1D_detectedV0s_DataLike_SmallpTinterval[iPt]->GetYaxis()->SetTitle("Counts");
    gauss[iPt]->GetHistogram()->GetYaxis()->SetTitle("efafesf");
    GaussPlusPolynom[iPt]->GetHistogram()->GetYaxis()->SetTitle("del_Phi");
    Int_t icolor=0;
    H1D_detectedV0s_DataLike_SmallpTinterval[iPt]->SetMarkerStyle(markers[0]);
    H1D_detectedV0s_DataLike_SmallpTinterval[iPt]->SetMarkerColor(colors [icolor]);
    H1D_detectedV0s_DataLike_SmallpTinterval[iPt]->SetLineColor  (colors [icolor]);
    // bkgparab[iPt]->Draw("same");
    // gauss[iPt]->Draw("same");


    // TLatex * text_pTbin = new TLatex (0.62,0.75,Form("pT [%.1f;%.1f] GeV", pTbins[iPt], pTbins[iPt+1]));
    // text_pTbin->SetNDC(kTRUE);
    // text_pTbin->Draw();
    TLatex * textContext = new TLatex (0.18,0.82,"#bf{ALICE Performance}"); //BOLD
    textContext->SetTextSize(0.05);
    textContext->SetNDC(kTRUE); //remove if I want x,y in TLatex to be in the coordinate system of the histogram
    textContext->Draw();
    // TLatex * textColl = new TLatex (0.18,0.75,"pp #sqrt{#it{s}} = 0.9 TeV, pilot beam 2021");
    TLatex * textColl = new TLatex (0.18,0.75,"pp #sqrt{#it{s}} = 0.9 TeV, Monte Carlo");
    textColl->SetTextSize(0.04);
    textColl->SetNDC(kTRUE); //remove if I want x,y in TLatex to be in the coordinate system of the histogram
    textColl->Draw();
    // TLatex * text_part = new TLatex (0.18,0.70,NamePart_Latex[ipart]+",   0 < #it{p}_{T} < 10 GeV/#it{c}");
    TLatex * text_part = new TLatex (0.18,0.70,NamePart_Latex[ipart]+Form(",   %.1f < #it{p}_{T} < %.1f GeV/#it{c}", pTbins[iPt], pTbins[iPt+1]));    
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

    //////////////////////////////////////////////////////////////////////////////////
    //////////////////////////// Signal Extraction ///////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////

    // //OLD signal extraction, my naive one
    // cout << "xmin_InvMass= " << liminf[ipart] << ", xmax_InvMass= " << limsup[ipart] << endl;
    // int ibin_liminf = H1D_detectedV0s_DataLike_SmallpTinterval[iPt]->GetXaxis()->FindBin(liminf[ipart]);
    // int ibin_limsup = H1D_detectedV0s_DataLike_SmallpTinterval[iPt]->GetXaxis()->FindBin(limsup[ipart]);
    // int nbinInvMass = ibin_limsup - ibin_liminf;
    // cout << "nbinInvMass = " << nbinInvMass << endl;
    // int totalCount = H1D_detectedV0s_DataLike_SmallpTinterval[iPt]->Integral(ibin_liminf+1,ibin_limsup-1); //integral of a hist counts the number of bins, does not care for bin width: not true integral
    // double backgroundCount = bkgparab[iPt]->Integral(liminf[ipart],limsup[ipart])*1./((limsup[ipart]-liminf[ipart])/nbinInvMass); //bkgparab is not an histogram but a function; here we need to convert it to a count
    // float SignalCountFit = totalCount - backgroundCount;
    // if (backgroundCount < 0) {
    //   cout << "WARNING - negative background count from fit - (cutArry_ID,ibinPt) = (" << ibinPt << "," << warning_cutArry_ID << ")" << endl;
    // }
    //NEW signal extraction, based on pp paper from David:

    float SignalMean = parGaussParab[iPt][1];
    float SignalStandardDev = parGaussParab[iPt][2];

    std::array<Double_t, 3> Signal_vector;

    if (SignalExtractionType == 0) {
      Signal_vector = GetTotalAndBackgroundCount_SidebandBackgroundExtrapolation(H1D_detectedV0s_DataLike_SmallpTinterval[iPt], SignalMean, SignalStandardDev, SideBandSizeMultiplierModifier);
    }
    if (SignalExtractionType == 1) {
      Signal_vector = GetTotalAndBackgroundCount_BackgroundIntegration(H1D_detectedV0s_DataLike_SmallpTinterval[iPt], SignalMean, SignalStandardDev, SideBandSizeMultiplierModifier, fFitResult[iPt], bkgparab[iPt]);
    }

    Double_t SignalCount = Signal_vector[0];
    Double_t SignalCount_StatError = Signal_vector[1];
    if (SignalCount < 0) {
      cout << "WARNING - negative SignalCount count from fit - (cutArry_ID,ibinPt) = (" << ibinPt << "," << warning_cutArry_ID << ")" << endl;
    }
    // cout << "totalCount= " << totalCount << ", backgroundCount= " << backgroundCount << endl;
    // cout << "SignalCount= " << SignalCount << endl;

    Double_t dpT = hPtDiff_RawYield_temp->GetXaxis()->GetBinWidth(ibinPt);
    Double_t drapidity = 1.5; //1.5
    Double_t d2N_dpTdy = SignalCount *1./dpT*1./drapidity;
    hPtDiff_RawYield_temp->SetBinContent(ibinPt,d2N_dpTdy);
    hPtDiff_RawYield_temp->SetBinError(ibinPt,sqrt(SignalCount_StatError*1./dpT*1./drapidity));  // error on d2N_dpTdy
    
    cout << "(ibinPt, SignalCount, TotalCount) = (" << ibinPt << ", " << SignalCount << ", " << Signal_vector[2] << ")" << endl;
  }

  canvasMassvsPtFits->SaveAs("Fits_invMass_ptDiff_"+NamePart[ipart]+".pdf");

  ////////////////////////////////// Fit End //////////////////////////////////
  // hSignalBackgroundRatio_vsPt = (TH1D*)hSignalCountFit_vsPt->Clone("hSignalBackgroundRatio_vsPt");
  // hSignalBackgroundRatio_vsPt->Divide(hSignalCountFit_vsPt,hBackgroundCountFit_vsPt);
  hRawYield_vsPt = (TH1D*)hPtDiff_RawYield_temp->Clone("hRawYield_vsPt");
  hRawYield_vsPt->Scale(1./SelectedEventCount);

  delete hPtDiff_RawYield_temp;
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
    Get_McEfficiency_vsPt(hMcEfficiency_vsPt[itableDcaCut], H3D_detectedV0s_TruePt_DcaHist, TrueV0PtSpectrum, ipart, pTbins, nbinpT, 0, ibinDcaCut);
    Get_McSignalBackgroundRatio_vsPt(hSignalBackgroundRatio_vsPt[itableDcaCut], H3D_detectedV0s_DataLike_DcaHist, ipart, pTbins, nbinpT, 0, ibinDcaCut);
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


void Efficiency_TpcCrossedRowsScan(TH1D* &hstat, TH1D* &hsyst, TH1D*&hsystCorr, Int_t ipart, TFile* file_O2Analysis, TString* SaveAs_Title, TString* &texXtitle, TString* &texYtitle, Double_t* pTbins, Int_t nbinpT, Double_t minTpcCrossedRowscut, Double_t maxTpcCrossedRowscut) {

  TH3D* H3D_detectedV0s_TruePt_TpcCrossedRowsHist = (TH3D*)file_O2Analysis->Get("lambdakzero-analysis-mc/h3dMass"+NamePart[ipart]+"TpcRows_MC_truePt");
  TH3D* H3D_detectedV0s_DataLike_TpcCrossedRowsHist = (TH3D*)file_O2Analysis->Get("lambdakzero-analysis-mc/h3dMass"+NamePart[ipart]+"TpcRows");
  TH1D* TrueV0PtSpectrum = (TH1D*)file_O2Analysis->Get("lambdakzero-particle-count-mc/h"+NamePart[ipart]+"Count_PtDiff");

  Int_t minTpcCrossedRowscutBin = H3D_detectedV0s_TruePt_TpcCrossedRowsHist->GetXaxis()->FindBin(minTpcCrossedRowscut);
  Int_t maxTpcCrossedRowscutBin = H3D_detectedV0s_TruePt_TpcCrossedRowsHist->GetXaxis()->FindBin(maxTpcCrossedRowscut);
  TH1D *hMcEfficiency_vsPt[maxTpcCrossedRowscutBin]; // exact size should be maxTpcCrossedRowscutBin - minTpcCrossedRowscutBin
  TH1D *hSignalBackgroundRatio_vsPt[maxTpcCrossedRowscutBin]; // but it shouldn't matter
  TH1D *hMcEfficiency_vsTpcCrossedRows[nbinpT];
  TH1D *hSignalBackgroundRatio_vsTpcCrossedRows[nbinpT];

  //Initialise hMcEfficiency_vsDca histograms for each pT bin
  for(int ibinPt = 1; ibinPt <= nbinpT; ibinPt++){
    int itablePt = ibinPt-1;
    // hMcEfficiency_vsDca[itablePt] = new TH1D(Form("hMcEfficiency_vsDca_Pt%i", ibinPt),Form("hMcEfficiency_vsDca_Pt%i", ibinPt),maxDCAcutBin,DCAbins_cut);
    hMcEfficiency_vsTpcCrossedRows[itablePt] = (TH1D*)H3D_detectedV0s_TruePt_TpcCrossedRowsHist->ProjectionX(Form("hMcEfficiency_vsTpcCrossedRows_Pt%i", ibinPt),minTpcCrossedRowscutBin,maxTpcCrossedRowscutBin,0,-1);
    hMcEfficiency_vsTpcCrossedRows[itablePt]->Reset("M");

    hSignalBackgroundRatio_vsTpcCrossedRows[itablePt] = (TH1D*)H3D_detectedV0s_DataLike_TpcCrossedRowsHist->ProjectionX(Form("hSignalBackgroundRatio_vsTpcCrossedRows_Pt%i", ibinPt),minTpcCrossedRowscutBin,maxTpcCrossedRowscutBin,0,-1);
    hSignalBackgroundRatio_vsTpcCrossedRows[itablePt]->Reset("M");
  }

  //loops over DCA cuts to get the "efficiency vs pT" histograms and use them to fill "efficiency vs DCAcut" 
  for(Int_t ibinTpcCrossedRowsCut = minTpcCrossedRowscutBin; ibinTpcCrossedRowsCut <= maxTpcCrossedRowscutBin; ibinTpcCrossedRowsCut++){
    int itableTpcCrossedRowsCut = ibinTpcCrossedRowsCut-1;
    Get_McEfficiency_vsPt(hMcEfficiency_vsPt[itableTpcCrossedRowsCut], H3D_detectedV0s_TruePt_TpcCrossedRowsHist, TrueV0PtSpectrum, ipart, pTbins, nbinpT, ibinTpcCrossedRowsCut, maxTpcCrossedRowscutBin);
    Get_McSignalBackgroundRatio_vsPt(hSignalBackgroundRatio_vsPt[itableTpcCrossedRowsCut], H3D_detectedV0s_DataLike_TpcCrossedRowsHist, ipart, pTbins, nbinpT, ibinTpcCrossedRowsCut, maxTpcCrossedRowscutBin);
    for(int ibinPt = 1; ibinPt <= nbinpT; ibinPt++){
      int itablePt = ibinPt-1;

      hMcEfficiency_vsTpcCrossedRows[itablePt]->SetBinContent(ibinTpcCrossedRowsCut,hMcEfficiency_vsPt[itableTpcCrossedRowsCut]->GetBinContent(ibinPt));
      hMcEfficiency_vsTpcCrossedRows[itablePt]->SetBinError(ibinTpcCrossedRowsCut,hMcEfficiency_vsPt[itableTpcCrossedRowsCut]->GetBinError(ibinPt));

      hSignalBackgroundRatio_vsTpcCrossedRows[itablePt]->SetBinContent(ibinTpcCrossedRowsCut,hSignalBackgroundRatio_vsPt[itableTpcCrossedRowsCut]->GetBinContent(ibinPt));
      hSignalBackgroundRatio_vsTpcCrossedRows[itablePt]->SetBinError(ibinTpcCrossedRowsCut,hSignalBackgroundRatio_vsPt[itableTpcCrossedRowsCut]->GetBinError(ibinPt));
    }
  }

  //Plots the "efficiency vs dca" for each pT bin
  TCanvas *canvasTpcCrossedRowsScanvsPt_eff = new TCanvas ("canvasTpcCrossedRowsScanvsPt_eff", NamehistoInvMass[ipart], 800, 600);
  // TCanvas *canvasDcaScanvsPt_eff = new TCanvas ("canvasDcaScanvsPt_eff", NamehistoInvMass[ipart], 1600, 1200);
  canvasTpcCrossedRowsScanvsPt_eff->Divide(std::ceil(std::sqrt(nbinpT)),std::ceil(std::sqrt(nbinpT)));
  TH1 *hFrame_eff[nbinpT];
  for(int ibinPt = 1; ibinPt <= nbinpT; ibinPt++){
    int itablePt = ibinPt-1; //different convention for bins that start at 1 and arrays that start at 0

    canvasTpcCrossedRowsScanvsPt_eff->cd(ibinPt); //ATTENTION: +1 because canvas divide starts counting from 1

    hFrame_eff[ibinPt] = canvasTpcCrossedRowsScanvsPt_eff->DrawFrame(minTpcCrossedRowscut,0,maxTpcCrossedRowscut,1.45*hMcEfficiency_vsTpcCrossedRows[itablePt]->GetMaximum());
    hFrame_eff[ibinPt]->SetXTitle("TpcCrossedRows cut (cm)");
    hFrame_eff[ibinPt]->SetYTitle("Efficiency");
    hMcEfficiency_vsTpcCrossedRows[itablePt]->Draw("E,same");
  }
  canvasTpcCrossedRowsScanvsPt_eff->SaveAs("Efficiency_vsPt_"+NamePart[ipart]+".pdf");

  //Plots the "SignalToBackgroundRatio vs dca" for each pT bin
  TCanvas *canvasTpcCrossedRowsScanvsPt_SBratio = new TCanvas ("canvasTpcCrossedRowsScanvsPt_SBratio", NamehistoInvMass[ipart], 800, 600);
  // TCanvas *canvasDcaScanvsPt_SBratio = new TCanvas ("canvasDcaScanvsPt_SBratio", NamehistoInvMass[ipart], 1600, 1200);
  canvasTpcCrossedRowsScanvsPt_SBratio->Divide(std::ceil(std::sqrt(nbinpT)),std::ceil(std::sqrt(nbinpT)));
  TH1 *hFrame_SBratio[nbinpT];
  for(int ibinPt = 1; ibinPt <= nbinpT; ibinPt++){
    int itablePt = ibinPt-1; //different convention for bins that start at 1 and arrays that start at 0

    canvasTpcCrossedRowsScanvsPt_SBratio->cd(ibinPt); //ATTENTION: +1 because canvas divide starts counting from 1

    hFrame_SBratio[ibinPt] = canvasTpcCrossedRowsScanvsPt_SBratio->DrawFrame(minTpcCrossedRowscutBin,0,maxTpcCrossedRowscut,1.45*hSignalBackgroundRatio_vsTpcCrossedRows[itablePt]->GetMaximum());
    hFrame_SBratio[ibinPt]->SetXTitle("TpcCrossedRows cut (cm)");
    hFrame_SBratio[ibinPt]->SetYTitle("Signal to Background Ratio");
    hSignalBackgroundRatio_vsTpcCrossedRows[itablePt]->Draw("E,same");
  }
  canvasTpcCrossedRowsScanvsPt_SBratio->SaveAs("SignalBackgroundRatio_vsPt_"+NamePart[ipart]+".pdf");

  ///////// standard stuff just to be able to use same global function at beginning of this file:

  hstat = (TH1D*)hMcEfficiency_vsTpcCrossedRows[0]->Clone("hstat");
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


void Get_Systematics_OneCut(TH1D* hSystematicUncertainty, TH1D* hSystematicUncertainty_PreBarlow, Int_t ipart, TFile **file_O2Analysis_CutVariation_array, Int_t id_referenceAnalysis, Int_t nCutsIterations, Double_t* pTbins, Int_t nbinpT, Int_t InvMassRebinFactor) {
  // get systematics from cut variation for the array of AnalysisResults.root files file_O2Analysis_CutVariation_array, all one variation of a single cut (cospa for example, or dcav0dau)
  // nCutsIterations is the number of variations on the cut being looked at; ie the size of the array file_O2Analysis_CutVariation_array

  TH3D* H3D_DetectedV0s_O2; // ? doesn't seem useful
  
  TH3D *H3D_DetectedV0s_MCtrue_cutVariation[nCutsIterations];
  TH3D *H3D_DetectedV0s_MCdatalike_cutVariation[nCutsIterations];
  TH1I *H1I_SelectedEventCount[nCutsIterations];
  TH1D *H1D_TrueV0PtSpectrum_cutVariation[nCutsIterations];

  TH1D *hMcEfficiency_vsPt[nCutsIterations];
  TH1D *hRawYield_vsPt[nCutsIterations];

  TH1D *hCorrectedYield_vsPt[nCutsIterations];

  for(int cutsArrayIterator = 0; cutsArrayIterator < nCutsIterations; cutsArrayIterator++){
    H3D_DetectedV0s_MCtrue_cutVariation[cutsArrayIterator] = (TH3D*)file_O2Analysis_CutVariation_array[cutsArrayIterator]->Get("lambdakzero-analysis-mc/h3dMass"+NamePart[ipart]+"_MC_truePt");
    H1I_SelectedEventCount[cutsArrayIterator] = (TH1I*)file_O2Analysis_CutVariation_array[cutsArrayIterator]->Get("lambdakzero-analysis-mc/hSelectedEventCounter");
    H1D_TrueV0PtSpectrum_cutVariation[cutsArrayIterator] = (TH1D*)file_O2Analysis_CutVariation_array[cutsArrayIterator]->Get("lambdakzero-particle-count-mc/h"+NamePart[ipart]+"Count_PtDiff");

    H3D_DetectedV0s_MCdatalike_cutVariation[cutsArrayIterator] = (TH3D*)file_O2Analysis_CutVariation_array[cutsArrayIterator]->Get("lambdakzero-analysis-mc/h3dMass"+NamePart[ipart]);

    Get_McEfficiency_vsPt(hMcEfficiency_vsPt[cutsArrayIterator], H3D_DetectedV0s_MCtrue_cutVariation[cutsArrayIterator], H1D_TrueV0PtSpectrum_cutVariation[cutsArrayIterator], ipart, pTbins, nbinpT, 0, -1);
    Get_RawYield_vsPt(hRawYield_vsPt[cutsArrayIterator], H3D_DetectedV0s_MCdatalike_cutVariation[cutsArrayIterator], H1I_SelectedEventCount[cutsArrayIterator], ipart, pTbins, nbinpT, 0, -1, cutsArrayIterator, SideBandSizeMultiplierModifier_array[0], SignalExtractionType_default, InvMassRebinFactor); ///real line that need to take the place of the one above once it's been written
    // Get_McEfficiency_vsPt(hRawYield_vsPt[cutsArrayIterator], H3D_DetectedV0s_MCdatalike_cutVariation[cutsArrayIterator], H1D_TrueV0PtSpectrum_cutVariation[cutsArrayIterator], ipart, pTbins, nbinpT, 0, -1); ///test purposes

    // //debug loop
    // for(int ibinPt = 1; ibinPt <= nbinpT; ibinPt++){
    //   cout << "cutID: " << cutsArrayIterator << ", pT bin: " << ibinPt << endl;
    //   cout << "hRawYield = " << hRawYield_vsPt[cutsArrayIterator]->GetBinContent(ibinPt) << ", hRawYield_error = " << hRawYield_vsPt[cutsArrayIterator]->GetBinError(ibinPt) << endl;
    //   cout << "hEfficiency = " << hMcEfficiency_vsPt[cutsArrayIterator]->GetBinContent(ibinPt) << ", hEfficiency_error = " << hMcEfficiency_vsPt[cutsArrayIterator]->GetBinError(ibinPt) << endl;
    // }
    hCorrectedYield_vsPt[cutsArrayIterator] = (TH1D*)hRawYield_vsPt[cutsArrayIterator]->Clone(Form("hCorrectedYield_vsPt_cutsArrayIterator%i", cutsArrayIterator));
    hCorrectedYield_vsPt[cutsArrayIterator]->Reset("M");
    // hCorrectedYield_vsPt[cutsArrayIterator]->Sumw2();
    // hRawYield_vsPt[cutsArrayIterator]->Sumw2();
    hCorrectedYield_vsPt[cutsArrayIterator]->Divide(hRawYield_vsPt[cutsArrayIterator],hMcEfficiency_vsPt[cutsArrayIterator], 1., 1., "b");
  // Efficiency errors are binomial-ish;  TGraphAsymmErrors::BayesDivide (https://root.cern.ch/doc/master/classTH1.html#ac782a09c31b4f7de40f8fd4f77efa090)
  // and https://inspirehep.net/files/57287ac8e45a976ab423f3dd456af694
  // PWGLF strangeness PAG says to just use option "b" for divide, and just see it as a binomial https://twiki.cern.ch/twiki/bin/view/ALICE/PWGLFPAGSTRANGENESSEfficiency
  }

  TH1D* hCorrectedYield_vsPt_REF = hCorrectedYield_vsPt[id_referenceAnalysis]; //this is assuming that's where I indeed put it
  // TH1D* hSystematicUncertainty[cutsIterator];
  Double_t hSigmaBarlow[nbinpT];
  Double_t CorrectedYieldDifference;
  int id_cutsArray_maxDeviation = 2;
  for(int ibinPt = 1; ibinPt <= nbinpT; ibinPt++){
    int PtArrayIterator = ibinPt - 1;
    Double_t MaxYieldDifference_cutsIteratorArray = 0;
    for(int cutsArrayIterator = 0; cutsArrayIterator < nCutsIterations; cutsArrayIterator++){
      CorrectedYieldDifference = abs(hCorrectedYield_vsPt[cutsArrayIterator]->GetBinContent(ibinPt) - hCorrectedYield_vsPt_REF->GetBinContent(ibinPt)); 
      if (CorrectedYieldDifference > MaxYieldDifference_cutsIteratorArray) {
        // MaxDifference_cutsIteratorArray = CorrectedYieldDifference;
        MaxYieldDifference_cutsIteratorArray = CorrectedYieldDifference;
        id_cutsArray_maxDeviation = cutsArrayIterator;
      }
    }

    Double_t SystUncertainty = MaxYieldDifference_cutsIteratorArray;

    // Barlow condition for cut variation systematics (Systematic Errors: facts and fictions, by Roger Barlow, https://arxiv.org/abs/hep-ex/0207026)
    Double_t StatUncertainty_REF = hCorrectedYield_vsPt_REF->GetBinError(ibinPt);
    Double_t StatUncertainty_MaxDeviationCase = hCorrectedYield_vsPt[id_cutsArray_maxDeviation]->GetBinError(ibinPt);
    
    hSigmaBarlow[PtArrayIterator] = sqrt(abs(StatUncertainty_MaxDeviationCase*StatUncertainty_MaxDeviationCase - StatUncertainty_REF*StatUncertainty_REF)); //stat error of the difference in the case of subsample
    
    hSystematicUncertainty_PreBarlow->SetBinContent(ibinPt,SystUncertainty);
    hSystematicUncertainty_PreBarlow->SetBinError(ibinPt,hSigmaBarlow[PtArrayIterator]);

    if (SystUncertainty > N_SigmaBarlow*hSigmaBarlow[PtArrayIterator]) { //Could ask for 1Sigma, 4Sigma or whatever depending on how conservative we want to be; one suggested in PWGLF note is 2Sigma
      hSystematicUncertainty->SetBinContent(ibinPt,SystUncertainty);
      hSystematicUncertainty->SetBinError(ibinPt,hSigmaBarlow[PtArrayIterator]);
    }
    else {
      hSystematicUncertainty->SetBinContent(ibinPt,0.);
    }
    // cout << "pT bin: " << ibinPt  << ", maxDeviation cut ID: " << id_cutsArray_maxDeviation << endl;
    // cout << "StatUncertainty_MaxDeviationCase = " << StatUncertainty_MaxDeviationCase  << ", StatUncertainty_REF = " << StatUncertainty_REF << endl;
    // here I ve got a nan (Not a number) problem with StatUncertainty_REF at pT bin 10;
    // cutID: 0, pT bin: 10
    // hRawYield = -0.00573374, hRawYield_error = nan    ; negative yield should probably be safeguarded; in fact negative background should be safeguarded

    cout << "Systematics CutVariation" << endl;
    cout << "(SystUncertainty, SigmaBarlow) = (" << SystUncertainty << "," << hSigmaBarlow[PtArrayIterator] << ") and check is " << (SystUncertainty > N_SigmaBarlow*hSigmaBarlow[PtArrayIterator]) << endl;
    cout << "CorrectedYield = " << hCorrectedYield_vsPt_REF->GetBinContent(ibinPt) << endl;
    cout << "Systematics/CorrectedYield = " << SystUncertainty/hCorrectedYield_vsPt_REF->GetBinContent(ibinPt) << endl;


  }
}


void Get_Systematics_OneCut_McTrueAndMcDatalikeDifferentCuts(TH1D* hSystematicUncertainty, TH1D* hSystematicUncertainty_PreBarlow, Int_t ipart, TFile **file_O2Analysis_CutVariation_Datalike_array, TFile **file_O2Analysis_CutVariation_True_array, Int_t id_referenceAnalysis, Int_t nCutsIterations, Double_t* pTbins, Int_t nbinpT, Int_t InvMassRebinFactor) {
  // get systematics from cut variation for the array of AnalysisResults.root files file_O2Analysis_CutVariation_array, all one variation of a single cut (cospa for example, or dcav0dau)
  // nCutsIterations is the number of variations on the cut being looked at; ie the size of the array file_O2Analysis_CutVariation_array

  TH3D* H3D_DetectedV0s_O2; // ? doesn't seem useful
  
  TH3D *H3D_DetectedV0s_MCtrue_cutVariation[nCutsIterations];
  TH3D *H3D_DetectedV0s_MCdatalike_cutVariation[nCutsIterations];
  TH1I *H1I_SelectedEventCount[nCutsIterations];
  TH1D *H1D_TrueV0PtSpectrum_cutVariation[nCutsIterations];

  TH1D *hMcEfficiency_vsPt[nCutsIterations];
  TH1D *hRawYield_vsPt[nCutsIterations];

  TH1D *hCorrectedYield_vsPt[nCutsIterations];

  for(int cutsArrayIterator = 0; cutsArrayIterator < nCutsIterations; cutsArrayIterator++){
    H3D_DetectedV0s_MCtrue_cutVariation[cutsArrayIterator] = (TH3D*)file_O2Analysis_CutVariation_True_array[cutsArrayIterator]->Get("lambdakzero-analysis-mc/h3dMass"+NamePart[ipart]+"_MC_truePt");
    H1I_SelectedEventCount[cutsArrayIterator] = (TH1I*)file_O2Analysis_CutVariation_Datalike_array[cutsArrayIterator]->Get("lambdakzero-analysis-mc/hSelectedEventCounter");
    H1D_TrueV0PtSpectrum_cutVariation[cutsArrayIterator] = (TH1D*)file_O2Analysis_CutVariation_True_array[cutsArrayIterator]->Get("lambdakzero-particle-count-mc/h"+NamePart[ipart]+"Count_PtDiff");

    H3D_DetectedV0s_MCdatalike_cutVariation[cutsArrayIterator] = (TH3D*)file_O2Analysis_CutVariation_Datalike_array[cutsArrayIterator]->Get("lambdakzero-analysis-mc/h3dMass"+NamePart[ipart]);

    Get_McEfficiency_vsPt(hMcEfficiency_vsPt[cutsArrayIterator], H3D_DetectedV0s_MCtrue_cutVariation[cutsArrayIterator], H1D_TrueV0PtSpectrum_cutVariation[cutsArrayIterator], ipart, pTbins, nbinpT, 0, -1);
    Get_RawYield_vsPt(hRawYield_vsPt[cutsArrayIterator], H3D_DetectedV0s_MCdatalike_cutVariation[cutsArrayIterator], H1I_SelectedEventCount[cutsArrayIterator], ipart, pTbins, nbinpT, 0, -1, cutsArrayIterator, SideBandSizeMultiplierModifier_array[0], SignalExtractionType_default, InvMassRebinFactor); ///real line that need to take the place of the one above once it's been written
    // Get_McEfficiency_vsPt(hRawYield_vsPt[cutsArrayIterator], H3D_DetectedV0s_MCdatalike_cutVariation[cutsArrayIterator], H1D_TrueV0PtSpectrum_cutVariation[cutsArrayIterator], ipart, pTbins, nbinpT, 0, -1); ///test purposes

    // //debug loop
    // for(int ibinPt = 1; ibinPt <= nbinpT; ibinPt++){
    //   cout << "cutID: " << cutsArrayIterator << ", pT bin: " << ibinPt << endl;
    //   cout << "hRawYield = " << hRawYield_vsPt[cutsArrayIterator]->GetBinContent(ibinPt) << ", hRawYield_error = " << hRawYield_vsPt[cutsArrayIterator]->GetBinError(ibinPt) << endl;
    //   cout << "hEfficiency = " << hMcEfficiency_vsPt[cutsArrayIterator]->GetBinContent(ibinPt) << ", hEfficiency_error = " << hMcEfficiency_vsPt[cutsArrayIterator]->GetBinError(ibinPt) << endl;
    // }
    hCorrectedYield_vsPt[cutsArrayIterator] = (TH1D*)hRawYield_vsPt[cutsArrayIterator]->Clone(Form("hCorrectedYield_vsPt_cutsArrayIterator%i", cutsArrayIterator));
    hCorrectedYield_vsPt[cutsArrayIterator]->Reset("M");
    // hCorrectedYield_vsPt[cutsArrayIterator]->Sumw2();
    // hRawYield_vsPt[cutsArrayIterator]->Sumw2();
    hCorrectedYield_vsPt[cutsArrayIterator]->Divide(hRawYield_vsPt[cutsArrayIterator],hMcEfficiency_vsPt[cutsArrayIterator], 1., 1., "b");
  // Efficiency errors are binomial-ish;  TGraphAsymmErrors::BayesDivide (https://root.cern.ch/doc/master/classTH1.html#ac782a09c31b4f7de40f8fd4f77efa090)
  // and https://inspirehep.net/files/57287ac8e45a976ab423f3dd456af694
  // PWGLF strangeness PAG says to just use option "b" for divide, and just see it as a binomial https://twiki.cern.ch/twiki/bin/view/ALICE/PWGLFPAGSTRANGENESSEfficiency
  }

  TH1D* hCorrectedYield_vsPt_REF = hCorrectedYield_vsPt[id_referenceAnalysis]; //this is assuming that's where I indeed put it
  // TH1D* hSystematicUncertainty[cutsIterator];
  Double_t hSigmaBarlow[nbinpT];
  Double_t CorrectedYieldDifference;
  int id_cutsArray_maxDeviation = 2;
  for(int ibinPt = 1; ibinPt <= nbinpT; ibinPt++){
    int PtArrayIterator = ibinPt - 1;
    Double_t MaxYieldDifference_cutsIteratorArray = 0;
    for(int cutsArrayIterator = 0; cutsArrayIterator < nCutsIterations; cutsArrayIterator++){
      CorrectedYieldDifference = abs(hCorrectedYield_vsPt[cutsArrayIterator]->GetBinContent(ibinPt) - hCorrectedYield_vsPt_REF->GetBinContent(ibinPt)); 
      if (CorrectedYieldDifference > MaxYieldDifference_cutsIteratorArray) {
        // MaxDifference_cutsIteratorArray = CorrectedYieldDifference;
        MaxYieldDifference_cutsIteratorArray = CorrectedYieldDifference;
        id_cutsArray_maxDeviation = cutsArrayIterator;
      }
    }

    Double_t SystUncertainty = MaxYieldDifference_cutsIteratorArray;

    // Barlow condition for cut variation systematics (Systematic Errors: facts and fictions, by Roger Barlow, https://arxiv.org/abs/hep-ex/0207026)
    Double_t StatUncertainty_REF = hCorrectedYield_vsPt_REF->GetBinError(ibinPt);
    Double_t StatUncertainty_MaxDeviationCase = hCorrectedYield_vsPt[id_cutsArray_maxDeviation]->GetBinError(ibinPt);
    
    hSigmaBarlow[PtArrayIterator] = sqrt(abs(StatUncertainty_MaxDeviationCase*StatUncertainty_MaxDeviationCase - StatUncertainty_REF*StatUncertainty_REF)); //stat error of the difference in the case of subsample
    
    hSystematicUncertainty_PreBarlow->SetBinContent(ibinPt,SystUncertainty);
    hSystematicUncertainty_PreBarlow->SetBinError(ibinPt,hSigmaBarlow[PtArrayIterator]);

    if (SystUncertainty > N_SigmaBarlow*hSigmaBarlow[PtArrayIterator]) { //Could ask for 1Sigma, 4Sigma or whatever depending on how conservative we want to be; one suggested in PWGLF note is 2Sigma
      hSystematicUncertainty->SetBinContent(ibinPt,SystUncertainty);
      hSystematicUncertainty->SetBinError(ibinPt,hSigmaBarlow[PtArrayIterator]);
    }
    else {
      hSystematicUncertainty->SetBinContent(ibinPt,0.);
    }
    // cout << "pT bin: " << ibinPt  << ", maxDeviation cut ID: " << id_cutsArray_maxDeviation << endl;
    // cout << "StatUncertainty_MaxDeviationCase = " << StatUncertainty_MaxDeviationCase  << ", StatUncertainty_REF = " << StatUncertainty_REF << endl;
    // here I ve got a nan (Not a number) problem with StatUncertainty_REF at pT bin 10;
    // cutID: 0, pT bin: 10
    // hRawYield = -0.00573374, hRawYield_error = nan    ; negative yield should probably be safeguarded; in fact negative background should be safeguarded

    cout << "Systematics CutVariation" << endl;
    cout << "(SystUncertainty, SigmaBarlow) = (" << SystUncertainty << "," << hSigmaBarlow[PtArrayIterator] << ") and check is " << (SystUncertainty > N_SigmaBarlow*hSigmaBarlow[PtArrayIterator]) << endl;
    cout << "CorrectedYield = " << hCorrectedYield_vsPt_REF->GetBinContent(ibinPt) << endl;
    cout << "Systematics/CorrectedYield = " << SystUncertainty/hCorrectedYield_vsPt_REF->GetBinContent(ibinPt) << endl;


  }
}

void Get_Systematics_Fit_SideBandVariation(TH1D* hSystematicUncertainty, TH1D* hSystematicUncertainty_PreBarlow, Int_t ipart, TH3D* H3D_detectedV0s_DataLike, TH1I* H1I_SelectedEventCount, TH3D* H3D_detectedV0s_TruePt, TH1D* TrueV0PtSpectrum, Double_t* pTbins, Int_t nbinpT, Int_t InvMassRebinFactor) {
// ongoing work: can I really use barlow there??
  // get systematics from cut variation for the array of AnalysisResults.root files file_O2Analysis_CutVariation_array, all one variation of a single cut (cospa for example, or dcav0dau)
  // nCutsIterations is the number of variations on the cut being looked at; ie the size of the array file_O2Analysis_CutVariation_array

  // Int_t SignalExtractionType = 1; default already defined beginning of this file

  TH1D *hMcEfficiency_vsPt[N_SideBandSizeVariation];
  TH1D *hRawYield_vsPt[N_SideBandSizeVariation];

  TH1D *hCorrectedYield_vsPt[N_SideBandSizeVariation];


  for(int id_SideBandSizeVariation = 0; id_SideBandSizeVariation < N_SideBandSizeVariation; id_SideBandSizeVariation++){
    cout << "SideBandSizeMultiplierModifier_array[id_SideBandSizeVariation]" << SideBandSizeMultiplierModifier_array[id_SideBandSizeVariation] << endl;
    Get_McEfficiency_vsPt(hMcEfficiency_vsPt[id_SideBandSizeVariation], H3D_detectedV0s_TruePt, TrueV0PtSpectrum, ipart, pTbins, nbinpT, 0, -1);
    Get_RawYield_vsPt(hRawYield_vsPt[id_SideBandSizeVariation], H3D_detectedV0s_DataLike, H1I_SelectedEventCount, ipart, pTbins, nbinpT, 0, -1, 0, SideBandSizeMultiplierModifier_array[id_SideBandSizeVariation], SignalExtractionType_default, InvMassRebinFactor); ///real line that need to take the place of the one above once it's been written

    // //debug loop
    // for(int ibinPt = 1; ibinPt <= nbinpT; ibinPt++){
    //   cout << "cutID: " << id_SideBandSizeVariation << ", pT bin: " << ibinPt << endl;
    //   cout << "hRawYield = " << hRawYield_vsPt[id_SideBandSizeVariation]->GetBinContent(ibinPt) << ", hRawYield_error = " << hRawYield_vsPt[id_SideBandSizeVariation]->GetBinError(ibinPt) << endl;
    //   cout << "hEfficiency = " << hMcEfficiency_vsPt[id_SideBandSizeVariation]->GetBinContent(ibinPt) << ", hEfficiency_error = " << hMcEfficiency_vsPt[id_SideBandSizeVariation]->GetBinError(ibinPt) << endl;
    // }
    hCorrectedYield_vsPt[id_SideBandSizeVariation] = (TH1D*)hRawYield_vsPt[id_SideBandSizeVariation]->Clone(Form("hCorrectedYield_vsPt_id_SideBandSizeVariation%i", id_SideBandSizeVariation));
    hCorrectedYield_vsPt[id_SideBandSizeVariation]->Reset("M");
    

    hCorrectedYield_vsPt[id_SideBandSizeVariation]->Divide(hRawYield_vsPt[id_SideBandSizeVariation],hMcEfficiency_vsPt[id_SideBandSizeVariation], 1., 1., "b");
  // Efficiency errors are binomial-ish;  TGraphAsymmErrors::BayesDivide (https://root.cern.ch/doc/master/classTH1.html#ac782a09c31b4f7de40f8fd4f77efa090)
  // and https://inspirehep.net/files/57287ac8e45a976ab423f3dd456af694
  // PWGLF strangeness PAG says to just use option "b" for divide, and just see it as a binomial https://twiki.cern.ch/twiki/bin/view/ALICE/PWGLFPAGSTRANGENESSEfficiency
  }

  TH1D* hCorrectedYield_vsPt_REF = hCorrectedYield_vsPt[id_referenceSideband]; //this is assuming that's where I indeed put it
  // TH1D* hSystematicUncertainty[cutsIterator];
  Double_t hSigmaBarlow[nbinpT];
  Double_t CorrectedYieldDifference;
  int id_SideBandSizeVariation_maxDeviation = 0;
  for(int ibinPt = 1; ibinPt <= nbinpT; ibinPt++){
    int PtArrayIterator = ibinPt - 1;
    Double_t MaxYieldDifference_SideBandVariation = 0;
    for(int id_SideBandSizeVariation = 0; id_SideBandSizeVariation < N_SideBandSizeVariation; id_SideBandSizeVariation++){
      CorrectedYieldDifference = abs(hCorrectedYield_vsPt[id_SideBandSizeVariation]->GetBinContent(ibinPt) - hCorrectedYield_vsPt_REF->GetBinContent(ibinPt)); 
      if (CorrectedYieldDifference > MaxYieldDifference_SideBandVariation) {
        // MaxDifference_SideBand = CorrectedYieldDifference;
        MaxYieldDifference_SideBandVariation = CorrectedYieldDifference;
        id_SideBandSizeVariation_maxDeviation = id_SideBandSizeVariation;
      }
    }
    cout << "MaxYieldDifference_SideBandVariation = " << MaxYieldDifference_SideBandVariation << endl;
    Double_t SystUncertainty = MaxYieldDifference_SideBandVariation;

    // Barlow condition for cut variation systematics (Systematic Errors: facts and fictions, by Roger Barlow, https://arxiv.org/abs/hep-ex/0207026)
    Double_t StatUncertainty_REF = hCorrectedYield_vsPt_REF->GetBinError(ibinPt);
    Double_t StatUncertainty_MaxDeviationCase = hCorrectedYield_vsPt[id_SideBandSizeVariation_maxDeviation]->GetBinError(ibinPt);
    
    hSigmaBarlow[PtArrayIterator] = sqrt(abs(StatUncertainty_MaxDeviationCase*StatUncertainty_MaxDeviationCase - StatUncertainty_REF*StatUncertainty_REF)); //stat error of the difference in the case of subsample
    
    hSystematicUncertainty_PreBarlow->SetBinContent(ibinPt,SystUncertainty);
    hSystematicUncertainty_PreBarlow->SetBinError(ibinPt,hSigmaBarlow[PtArrayIterator]);

    if (SystUncertainty > N_SigmaBarlow*hSigmaBarlow[PtArrayIterator]) { //Could ask for 1Sigma, 4Sigma or whatever depending on how conservative we want to be; one suggested in PWGLF note is 2Sigma
      hSystematicUncertainty->SetBinContent(ibinPt,SystUncertainty);
      hSystematicUncertainty->SetBinError(ibinPt,hSigmaBarlow[PtArrayIterator]);
    }
    else {
      hSystematicUncertainty->SetBinContent(ibinPt,0.);
    }
    // cout << "pT bin: " << ibinPt  << ", maxDeviation cut ID: " << id_cutsArray_maxDeviation << endl;
    // cout << "StatUncertainty_MaxDeviationCase = " << StatUncertainty_MaxDeviationCase  << ", StatUncertainty_REF = " << StatUncertainty_REF << endl;
    // here I ve got a nan (Not a number) problem with StatUncertainty_REF at pT bin 10;
    // cutID: 0, pT bin: 10
    // hRawYield = -0.00573374, hRawYield_error = nan    ; negative yield should probably be safeguarded; in fact negative background should be safeguarded

    cout << "Systematics SideBandVariation" << endl;
    cout << "(SystUncertainty, SigmaBarlow) = (" << SystUncertainty << "," << hSigmaBarlow[PtArrayIterator] << ") and check is " << (SystUncertainty > N_SigmaBarlow*hSigmaBarlow[PtArrayIterator]) << endl;
    cout << "CorrectedYield = " << hCorrectedYield_vsPt_REF->GetBinContent(ibinPt) << endl;
    cout << "Systematics/CorrectedYield = " << SystUncertainty/hCorrectedYield_vsPt_REF->GetBinContent(ibinPt) << endl;


  }
}

void Get_Systematics_Fit_SignalExtractionType(TH1D* hSystematicUncertainty, TH1D* hSystematicUncertainty_PreBarlow, Int_t ipart, TH3D* H3D_detectedV0s_DataLike, TH1I* H1I_SelectedEventCount, TH3D* H3D_detectedV0s_TruePt, TH1D* TrueV0PtSpectrum, Double_t* pTbins, Int_t nbinpT, Int_t InvMassRebinFactor) {
  // ongoing work: can I really use barlow there??
  // get systematics from cut variation for the array of AnalysisResults.root files file_O2Analysis_CutVariation_array, all one variation of a single cut (cospa for example, or dcav0dau)
  // nCutsIterations is the number of variations on the cut being looked at; ie the size of the array file_O2Analysis_CutVariation_array
  Int_t N_SignalExtractionType = 2;

  TH1D *hMcEfficiency_vsPt[N_SignalExtractionType];
  TH1D *hRawYield_vsPt[N_SignalExtractionType];

  TH1D *hCorrectedYield_vsPt[N_SignalExtractionType];


  for(int SignalExtractionType_loopID = 0; SignalExtractionType_loopID < N_SignalExtractionType; SignalExtractionType_loopID++){

    Get_McEfficiency_vsPt(hMcEfficiency_vsPt[SignalExtractionType_loopID], H3D_detectedV0s_TruePt, TrueV0PtSpectrum, ipart, pTbins, nbinpT, 0, -1);
    Get_RawYield_vsPt(hRawYield_vsPt[SignalExtractionType_loopID], H3D_detectedV0s_DataLike, H1I_SelectedEventCount, ipart, pTbins, nbinpT, 0, -1, 0, SideBandSizeMultiplierModifier_array[0], SignalExtractionType_loopID, InvMassRebinFactor); ///real line that need to take the place of the one above once it's been written

    // //debug loop
    // for(int ibinPt = 1; ibinPt <= nbinpT; ibinPt++){
    //   cout << "cutID: " << SignalExtractionType << ", pT bin: " << ibinPt << endl;
    //   cout << "hRawYield = " << hRawYield_vsPt[SignalExtractionType_loopID]->GetBinContent(ibinPt) << ", hRawYield_error = " << hRawYield_vsPt[SignalExtractionType_loopID]->GetBinError(ibinPt) << endl;
    //   cout << "hEfficiency = " << hMcEfficiency_vsPt[SignalExtractionType_loopID]->GetBinContent(ibinPt) << ", hEfficiency_error = " << hMcEfficiency_vsPt[SignalExtractionType_loopID]->GetBinError(ibinPt) << endl;
    // }
    hCorrectedYield_vsPt[SignalExtractionType_loopID] = (TH1D*)hRawYield_vsPt[SignalExtractionType_loopID]->Clone(Form("hCorrectedYield_vsPt_SignalExtractionType_loopID%i", SignalExtractionType_loopID));
    hCorrectedYield_vsPt[SignalExtractionType_loopID]->Reset("M");
    

    hCorrectedYield_vsPt[SignalExtractionType_loopID]->Divide(hRawYield_vsPt[SignalExtractionType_loopID],hMcEfficiency_vsPt[SignalExtractionType_loopID], 1., 1., "b");
  // Efficiency errors are binomial-ish;  TGraphAsymmErrors::BayesDivide (https://root.cern.ch/doc/master/classTH1.html#ac782a09c31b4f7de40f8fd4f77efa090)
  // and https://inspirehep.net/files/57287ac8e45a976ab423f3dd456af694
  // PWGLF strangeness PAG says to just use option "b" for divide, and just see it as a binomial https://twiki.cern.ch/twiki/bin/view/ALICE/PWGLFPAGSTRANGENESSEfficiency
  }

  TH1D* hCorrectedYield_vsPt_REF = hCorrectedYield_vsPt[id_referenceAnalysis]; //this is assuming that's where I indeed put it
  // TH1D* hSystematicUncertainty[cutsIterator];
  Double_t hSigmaBarlow[nbinpT];
  Double_t CorrectedYieldDifference;
  int id_SignalExtractionType_maxDeviation = 0;
  for(int ibinPt = 1; ibinPt <= nbinpT; ibinPt++){
    int PtArrayIterator = ibinPt - 1;
    Double_t MaxYieldDifference_SignalExtractionVariation = 0;
    for(int SignalExtractionType_loopID = 0; SignalExtractionType_loopID < N_SignalExtractionType; SignalExtractionType_loopID++){
      CorrectedYieldDifference = abs(hCorrectedYield_vsPt[SignalExtractionType_loopID]->GetBinContent(ibinPt) - hCorrectedYield_vsPt_REF->GetBinContent(ibinPt)); 
      if (CorrectedYieldDifference > MaxYieldDifference_SignalExtractionVariation) {
        // MaxYieldDifference_SignalExtractionVariation = CorrectedYieldDifference;
        MaxYieldDifference_SignalExtractionVariation = CorrectedYieldDifference;
        id_SignalExtractionType_maxDeviation = SignalExtractionType_loopID;
      }
    }

    Double_t SystUncertainty = MaxYieldDifference_SignalExtractionVariation;

    // Barlow condition for cut variation systematics (Systematic Errors: facts and fictions, by Roger Barlow, https://arxiv.org/abs/hep-ex/0207026)
    Double_t StatUncertainty_REF = hCorrectedYield_vsPt_REF->GetBinError(ibinPt);
    Double_t StatUncertainty_MaxDeviationCase = hCorrectedYield_vsPt[id_SignalExtractionType_maxDeviation]->GetBinError(ibinPt);
    
    hSigmaBarlow[PtArrayIterator] = sqrt(abs(StatUncertainty_MaxDeviationCase*StatUncertainty_MaxDeviationCase - StatUncertainty_REF*StatUncertainty_REF)); //stat error of the difference in the case of subsample
    
    hSystematicUncertainty_PreBarlow->SetBinContent(ibinPt,SystUncertainty);
    hSystematicUncertainty_PreBarlow->SetBinError(ibinPt,hSigmaBarlow[PtArrayIterator]);

    if (SystUncertainty > N_SigmaBarlow*hSigmaBarlow[PtArrayIterator]) { //Could ask for 1Sigma, 4Sigma or whatever depending on how conservative we want to be; one suggested in PWGLF note is 2Sigma
      hSystematicUncertainty->SetBinContent(ibinPt,SystUncertainty);
      hSystematicUncertainty->SetBinError(ibinPt,hSigmaBarlow[PtArrayIterator]);
    }
    else {
      hSystematicUncertainty->SetBinContent(ibinPt,0.);
    }
    // cout << "pT bin: " << ibinPt  << ", maxDeviation cut ID: " << id_cutsArray_maxDeviation << endl;
    // cout << "StatUncertainty_MaxDeviationCase = " << StatUncertainty_MaxDeviationCase  << ", StatUncertainty_REF = " << StatUncertainty_REF << endl;
    // here I ve got a nan (Not a number) problem with StatUncertainty_REF at pT bin 10;
    // cutID: 0, pT bin: 10
    // hRawYield = -0.00573374, hRawYield_error = nan    ; negative yield should probably be safeguarded; in fact negative background should be safeguarded

    cout << "Systematics SignalExtractionType" << endl;
    cout << "(SystUncertainty, SigmaBarlow) = (" << SystUncertainty << "," << hSigmaBarlow[PtArrayIterator] << ") and check is " << (SystUncertainty > N_SigmaBarlow*hSigmaBarlow[PtArrayIterator]) << endl;
    cout << "CorrectedYield = " << hCorrectedYield_vsPt_REF->GetBinContent(ibinPt) << endl;
    cout << "Systematics/CorrectedYield = " << SystUncertainty/hCorrectedYield_vsPt_REF->GetBinContent(ibinPt) << endl;


  }
}


void Get_Systematics_Fit_Binning(TH1D* hSystematicUncertainty, TH1D* hSystematicUncertainty_PreBarlow, Int_t ipart, TH3D* H3D_detectedV0s_DataLike, TH1I* H1I_SelectedEventCount, TH3D* H3D_detectedV0s_TruePt, TH1D* TrueV0PtSpectrum, Double_t* pTbins, Int_t nbinpT) {
// ongoing work: can I really use barlow there??
  // get systematics from cut variation for the array of AnalysisResults.root files file_O2Analysis_CutVariation_array, all one variation of a single cut (cospa for example, or dcav0dau)
  // nCutsIterations is the number of variations on the cut being looked at; ie the size of the array file_O2Analysis_CutVariation_array

  // Int_t SignalExtractionType = 1; default already defined beginning of this file

  TH1D *hMcEfficiency_vsPt[N_BinningVariation];
  TH1D *hRawYield_vsPt[N_BinningVariation];

  TH1D *hCorrectedYield_vsPt[N_BinningVariation];


  for(int id_BinningVariation = 0; id_BinningVariation < N_BinningVariation; id_BinningVariation++){
    cout << "InvMassRebinFactor_variation_array[id_BinningVariation]" << InvMassRebinFactor_variation_array[id_BinningVariation] << endl;
    Get_McEfficiency_vsPt(hMcEfficiency_vsPt[id_BinningVariation], H3D_detectedV0s_TruePt, TrueV0PtSpectrum, ipart, pTbins, nbinpT, 0, -1);
    Get_RawYield_vsPt(hRawYield_vsPt[id_BinningVariation], H3D_detectedV0s_DataLike, H1I_SelectedEventCount, ipart, pTbins, nbinpT, 0, -1, 0, 0, SignalExtractionType_default, InvMassRebinFactor_variation_array[ipart][id_BinningVariation]); ///real line that need to take the place of the one above once it's been written

    // //debug loop
    // for(int ibinPt = 1; ibinPt <= nbinpT; ibinPt++){
    //   cout << "cutID: " << id_BinningVariation << ", pT bin: " << ibinPt << endl;
    //   cout << "hRawYield = " << hRawYield_vsPt[id_BinningVariation]->GetBinContent(ibinPt) << ", hRawYield_error = " << hRawYield_vsPt[id_BinningVariation]->GetBinError(ibinPt) << endl;
    //   cout << "hEfficiency = " << hMcEfficiency_vsPt[id_BinningVariation]->GetBinContent(ibinPt) << ", hEfficiency_error = " << hMcEfficiency_vsPt[id_BinningVariation]->GetBinError(ibinPt) << endl;
    // }
    hCorrectedYield_vsPt[id_BinningVariation] = (TH1D*)hRawYield_vsPt[id_BinningVariation]->Clone(Form("hCorrectedYield_vsPt_id_BinningVariation%i", id_BinningVariation));
    hCorrectedYield_vsPt[id_BinningVariation]->Reset("M");
    

    hCorrectedYield_vsPt[id_BinningVariation]->Divide(hRawYield_vsPt[id_BinningVariation],hMcEfficiency_vsPt[id_BinningVariation], 1., 1., "b");
  // Efficiency errors are binomial-ish;  TGraphAsymmErrors::BayesDivide (https://root.cern.ch/doc/master/classTH1.html#ac782a09c31b4f7de40f8fd4f77efa090)
  // and https://inspirehep.net/files/57287ac8e45a976ab423f3dd456af694
  // PWGLF strangeness PAG says to just use option "b" for divide, and just see it as a binomial https://twiki.cern.ch/twiki/bin/view/ALICE/PWGLFPAGSTRANGENESSEfficiency
  }

  TH1D* hCorrectedYield_vsPt_REF = hCorrectedYield_vsPt[id_referenceSideband]; //this is assuming that's where I indeed put it
  // TH1D* hSystematicUncertainty[cutsIterator];
  Double_t hSigmaBarlow[nbinpT];
  Double_t CorrectedYieldDifference;
  int id_BinningVariation_maxDeviation = 0;
  for(int ibinPt = 1; ibinPt <= nbinpT; ibinPt++){
    int PtArrayIterator = ibinPt - 1;
    Double_t MaxYieldDifference_BinningVariation = 0;
    for(int id_BinningVariation = 0; id_BinningVariation < N_BinningVariation; id_BinningVariation++){
      CorrectedYieldDifference = abs(hCorrectedYield_vsPt[id_BinningVariation]->GetBinContent(ibinPt) - hCorrectedYield_vsPt_REF->GetBinContent(ibinPt)); 
      if (CorrectedYieldDifference > MaxYieldDifference_BinningVariation) {
        // MaxDifference_SideBand = CorrectedYieldDifference;
        MaxYieldDifference_BinningVariation = CorrectedYieldDifference;
        id_BinningVariation_maxDeviation = id_BinningVariation;
      }
    }
    cout << "MaxYieldDifference_BinningVariation = " << MaxYieldDifference_BinningVariation << endl;
    Double_t SystUncertainty = MaxYieldDifference_BinningVariation;

    // Barlow condition for cut variation systematics (Systematic Errors: facts and fictions, by Roger Barlow, https://arxiv.org/abs/hep-ex/0207026)
    Double_t StatUncertainty_REF = hCorrectedYield_vsPt_REF->GetBinError(ibinPt);
    Double_t StatUncertainty_MaxDeviationCase = hCorrectedYield_vsPt[id_BinningVariation_maxDeviation]->GetBinError(ibinPt);
    
    hSigmaBarlow[PtArrayIterator] = sqrt(abs(StatUncertainty_MaxDeviationCase*StatUncertainty_MaxDeviationCase - StatUncertainty_REF*StatUncertainty_REF)); //stat error of the difference in the case of subsample
    
    hSystematicUncertainty_PreBarlow->SetBinContent(ibinPt,SystUncertainty);
    hSystematicUncertainty_PreBarlow->SetBinError(ibinPt,hSigmaBarlow[PtArrayIterator]);

    if (SystUncertainty > N_SigmaBarlow*hSigmaBarlow[PtArrayIterator]) { //Could ask for 1Sigma, 4Sigma or whatever depending on how conservative we want to be; one suggested in PWGLF note is 2Sigma
      hSystematicUncertainty->SetBinContent(ibinPt,SystUncertainty);
      hSystematicUncertainty->SetBinError(ibinPt,hSigmaBarlow[PtArrayIterator]);
    }
    else {
      hSystematicUncertainty->SetBinContent(ibinPt,0.);
    }
    // cout << "pT bin: " << ibinPt  << ", maxDeviation cut ID: " << id_cutsArray_maxDeviation << endl;
    // cout << "StatUncertainty_MaxDeviationCase = " << StatUncertainty_MaxDeviationCase  << ", StatUncertainty_REF = " << StatUncertainty_REF << endl;
    // here I ve got a nan (Not a number) problem with StatUncertainty_REF at pT bin 10;
    // cutID: 0, pT bin: 10
    // hRawYield = -0.00573374, hRawYield_error = nan    ; negative yield should probably be safeguarded; in fact negative background should be safeguarded

    cout << "Systematics BinningVariation" << endl;
    cout << "(SystUncertainty, SigmaBarlow) = (" << SystUncertainty << "," << hSigmaBarlow[PtArrayIterator] << ") and check is " << (SystUncertainty > N_SigmaBarlow*hSigmaBarlow[PtArrayIterator]) << endl;
    cout << "CorrectedYield = " << hCorrectedYield_vsPt_REF->GetBinContent(ibinPt) << endl;
    cout << "Systematics/CorrectedYield = " << SystUncertainty/hCorrectedYield_vsPt_REF->GetBinContent(ibinPt) << endl;


  }
}


void Systematics_CutVariations_Graphs(TH1D* &hstat, TH1D* &hsyst, TH1D* &hsystCorr, Int_t ipart, TFile** file_O2Analysis_CutVariation_array, Int_t id_referenceAnalysis, Int_t nCutsIterations, TString* SaveAs_Title, TString* &texXtitle, TString* &texYtitle, Double_t* pTbins, Int_t nbinpT) {
  Int_t aaaa = nCutsIterations; //to remove warning

  TH3D* H3D_detectedV0s_TruePt = (TH3D*)file_O2Analysis_CutVariation_array[id_referenceAnalysis]->Get("lambdakzero-analysis-mc/h3dMass"+NamePart[ipart]+"_MC_truePt");
  TH3D* H3D_detectedV0s_DataLike = (TH3D*)file_O2Analysis_CutVariation_array[id_referenceAnalysis]->Get("lambdakzero-analysis-mc/h3dMass"+NamePart[ipart]);
  TH1D* TrueV0PtSpectrum = (TH1D*)file_O2Analysis_CutVariation_array[id_referenceAnalysis]->Get("lambdakzero-particle-count-mc/h"+NamePart[ipart]+"Count_PtDiff");
  TH1I* H1I_SelectedEventCount = (TH1I*)file_O2Analysis_CutVariation_array[id_referenceAnalysis]->Get("lambdakzero-analysis-mc/hSelectedEventCounter");

  TH1D* hSystematicUncertainty_CutVariation = new TH1D("hSystematicUncertainty_CutVariation", "hSystematicUncertainty_CutVariation", nbinpT, pTbins);
  TH1D* hSystematicUncertainty_CutVariation_PreBarlow = new TH1D("hSystematicUncertainty_CutVariation_PreBarlow", "hSystematicUncertainty_CutVariation_PreBarlow", nbinpT, pTbins);

  TH1D* hSystematicUncertainty_Fit_SidebandVariation = new TH1D("hSystematicUncertainty_Fit_SidebandVariation", "hSystematicUncertainty_Fit_SidebandVariation", nbinpT, pTbins);
  TH1D* hSystematicUncertainty_Fit_SidebandVariation_PreBarlow = new TH1D("hSystematicUncertainty_Fit_SidebandVariation_PreBarlow", "hSystematicUncertainty_Fit_SidebandVariation_PreBarlow", nbinpT, pTbins);
  TH1D* hSystematicUncertainty_Fit_SignalExtractionType = new TH1D("hSystematicUncertainty_Fit_SignalExtractionType", "hSystematicUncertainty_Fit_SignalExtractionType", nbinpT, pTbins);
  TH1D* hSystematicUncertainty_Fit_SignalExtractionType_PreBarlow = new TH1D("hSystematicUncertainty_Fit_SignalExtractionType_PreBarlow", "hSystematicUncertainty_Fit_SignalExtractionType_PreBarlow", nbinpT, pTbins);
  TH1D* hSystematicUncertainty_Fit_BinningVariation = new TH1D("hSystematicUncertainty_Fit_BinningVariation", "hSystematicUncertainty_Fit_BinningVariation", nbinpT, pTbins);
  TH1D* hSystematicUncertainty_Fit_BinningVariation_PreBarlow = new TH1D("hSystematicUncertainty_Fit_BinningVariation_PreBarlow", "hSystematicUncertainty_Fit_BinningVariation_PreBarlow", nbinpT, pTbins);

  TH1D* hCorrectedYield_vsPt = new TH1D("hCorrectedYield_vsPt", "hCorrectedYield_vsPt", nbinpT, pTbins);
  TH1D* hSystematicUncertainty_Total = new TH1D("hSystematicUncertainty_Total", "hSystematicUncertainty_Total", nbinpT, pTbins);
  TH1D* hSystematicUncertainty_Total_PreBarlow = new TH1D("hSystematicUncertainty_Total_PreBarlow", "hSystematicUncertainty_Total_PreBarlow", nbinpT, pTbins);

  hSystematicUncertainty_CutVariation->Sumw2();
  hSystematicUncertainty_CutVariation_PreBarlow->Sumw2();
  hSystematicUncertainty_Fit_SidebandVariation->Sumw2();
  hSystematicUncertainty_Fit_SidebandVariation_PreBarlow->Sumw2();
  hSystematicUncertainty_Fit_SignalExtractionType->Sumw2();
  hSystematicUncertainty_Fit_SignalExtractionType_PreBarlow->Sumw2();
  hSystematicUncertainty_Fit_BinningVariation->Sumw2();
  hSystematicUncertainty_Fit_BinningVariation_PreBarlow->Sumw2();
  hSystematicUncertainty_Total->Sumw2();
  hSystematicUncertainty_Total_PreBarlow->Sumw2();

  TH1D *hMcEfficiency_vsPt;
  TH1D *hRawYield_vsPt;
  
  Int_t InvMassRebinFactor = InvMassRebinFactor_standard[ipart];

  Get_McEfficiency_vsPt(hMcEfficiency_vsPt, H3D_detectedV0s_TruePt, TrueV0PtSpectrum, ipart, pTbins, nbinpT, 0, -1);
  Get_RawYield_vsPt(hRawYield_vsPt, H3D_detectedV0s_DataLike, H1I_SelectedEventCount, ipart, pTbins, nbinpT, 0, -1, 0, SideBandSizeMultiplierModifier_array[0], SignalExtractionType_default, InvMassRebinFactor); ///real line that need to take the place of the one above once it's been written

  Get_Systematics_OneCut(hSystematicUncertainty_CutVariation, hSystematicUncertainty_CutVariation_PreBarlow, ipart, file_O2Analysis_CutVariation_array, id_referenceAnalysis, nCutsIterations, pTbins, nbinpT, InvMassRebinFactor);
  Get_Systematics_Fit_SideBandVariation(hSystematicUncertainty_Fit_SidebandVariation, hSystematicUncertainty_Fit_SidebandVariation_PreBarlow, ipart, H3D_detectedV0s_DataLike, H1I_SelectedEventCount, H3D_detectedV0s_TruePt, TrueV0PtSpectrum, pTbins, nbinpT, InvMassRebinFactor);
  Get_Systematics_Fit_SignalExtractionType(hSystematicUncertainty_Fit_SignalExtractionType, hSystematicUncertainty_Fit_SignalExtractionType_PreBarlow, ipart, H3D_detectedV0s_DataLike, H1I_SelectedEventCount, H3D_detectedV0s_TruePt, TrueV0PtSpectrum, pTbins, nbinpT, InvMassRebinFactor); //very large systematics as is
  Get_Systematics_Fit_Binning(hSystematicUncertainty_Fit_BinningVariation, hSystematicUncertainty_Fit_BinningVariation_PreBarlow, ipart, H3D_detectedV0s_DataLike, H1I_SelectedEventCount, H3D_detectedV0s_TruePt, TrueV0PtSpectrum, pTbins, nbinpT);


  hCorrectedYield_vsPt->Sumw2();
  hCorrectedYield_vsPt->Divide(hRawYield_vsPt,hMcEfficiency_vsPt, 1., 1., "b");

  hSystematicUncertainty_CutVariation->Divide(hCorrectedYield_vsPt); //get it as a ratio of ref corrected yield
  hSystematicUncertainty_CutVariation_PreBarlow->Divide(hCorrectedYield_vsPt); //get it as a ratio of ref corrected yield
  hSystematicUncertainty_Fit_SidebandVariation->Divide(hCorrectedYield_vsPt); //get it as a ratio of ref corrected yield
  hSystematicUncertainty_Fit_SidebandVariation_PreBarlow->Divide(hCorrectedYield_vsPt); //get it as a ratio of ref corrected yield
  hSystematicUncertainty_Fit_SignalExtractionType->Divide(hCorrectedYield_vsPt); //get it as a ratio of ref corrected yield
  hSystematicUncertainty_Fit_SignalExtractionType_PreBarlow->Divide(hCorrectedYield_vsPt); //get it as a ratio of ref corrected yield
  hSystematicUncertainty_Fit_BinningVariation->Divide(hCorrectedYield_vsPt); //get it as a ratio of ref corrected yield
  hSystematicUncertainty_Fit_BinningVariation_PreBarlow->Divide(hCorrectedYield_vsPt); //get it as a ratio of ref corrected yield

  hSystematicUncertainty_Total->Add(hSystematicUncertainty_CutVariation);
  hSystematicUncertainty_Total->Add(hSystematicUncertainty_Fit_SidebandVariation);
  hSystematicUncertainty_Total->Add(hSystematicUncertainty_Fit_SignalExtractionType);
  hSystematicUncertainty_Total->Add(hSystematicUncertainty_Fit_BinningVariation);

  hSystematicUncertainty_Total_PreBarlow->Add(hSystematicUncertainty_CutVariation_PreBarlow);
  hSystematicUncertainty_Total_PreBarlow->Add(hSystematicUncertainty_Fit_SidebandVariation_PreBarlow);
  hSystematicUncertainty_Total_PreBarlow->Add(hSystematicUncertainty_Fit_SignalExtractionType_PreBarlow);
  hSystematicUncertainty_Total_PreBarlow->Add(hSystematicUncertainty_Fit_BinningVariation_PreBarlow);

  // //to check statistical uncertainties
  // for(int ibinPt = 1; ibinPt <= nbinpT; ibinPt++){
  //   hSystematicUncertainty_CutVariation_PreBarlow->SetBinContent(ibinPt,hCorrectedYield_vsPt->GetBinContent(ibinPt));
  //   hSystematicUncertainty_CutVariation_PreBarlow->SetBinError(ibinPt,hCorrectedYield_vsPt->GetBinError(ibinPt));
  //   hSystematicUncertainty_CutVariation->SetBinContent(ibinPt,hCorrectedYield_vsPt->GetBinContent(ibinPt));
  //   hSystematicUncertainty_CutVariation->SetBinError(ibinPt,hCorrectedYield_vsPt->GetBinError(ibinPt));
  // }

  //////////////////////////////////
  //////////// Plots ///////////////
  //////////////////////////////////
  //Post Barlow
  TCanvas *canvasSyst = new TCanvas ("canvasSyst", NamehistoInvMass[ipart], 1200, 800);
  canvasSyst->cd(0); //ATTENTION: +1 because canvas divide starts counting from 1
  TH1 *hFrame = canvasSyst->DrawFrame(pTbins[0],0,pTbins[nbinpT],1.2*hSystematicUncertainty_Total->GetMaximum());
  hFrame->SetXTitle("#it{p}_{T} (GeV/#it{c})");
  hFrame->SetYTitle("uncertainty (ratio)");

  hSystematicUncertainty_CutVariation->Draw("E,same");
  hSystematicUncertainty_Fit_SidebandVariation->Draw("E,same");
  hSystematicUncertainty_Fit_SignalExtractionType->Draw("E,same");
  hSystematicUncertainty_Fit_BinningVariation->Draw("E,same");
  hSystematicUncertainty_Total->Draw("E,same");

  Int_t icolor_CutVariation = 1;
  Int_t icolor_Fit_SidebandVariation = 2;
  Int_t icolor_Fit_SignalExtractionType = 3;
  Int_t icolor_Fit_BinningVariation = 4;
  Int_t icolor_Total = 0;

  hSystematicUncertainty_CutVariation->SetMarkerStyle(markers[0]);
  hSystematicUncertainty_CutVariation->SetMarkerColor(colors [icolor_CutVariation]);
  hSystematicUncertainty_CutVariation->SetLineColor  (colors [icolor_CutVariation]);

  hSystematicUncertainty_Fit_SidebandVariation->SetMarkerStyle(markers[0]);
  hSystematicUncertainty_Fit_SidebandVariation->SetMarkerColor(colors [icolor_Fit_SidebandVariation]);
  hSystematicUncertainty_Fit_SidebandVariation->SetLineColor  (colors [icolor_Fit_SidebandVariation]);

  hSystematicUncertainty_Fit_SignalExtractionType->SetMarkerStyle(markers[0]);
  hSystematicUncertainty_Fit_SignalExtractionType->SetMarkerColor(colors [icolor_Fit_SignalExtractionType]);
  hSystematicUncertainty_Fit_SignalExtractionType->SetLineColor  (colors [icolor_Fit_SignalExtractionType]);

  hSystematicUncertainty_Fit_BinningVariation->SetMarkerStyle(markers[0]);
  hSystematicUncertainty_Fit_BinningVariation->SetMarkerColor(colors [icolor_Fit_BinningVariation]);
  hSystematicUncertainty_Fit_BinningVariation->SetLineColor  (colors [icolor_Fit_BinningVariation]);

  hSystematicUncertainty_Total->SetMarkerStyle(markers[0]);
  hSystematicUncertainty_Total->SetMarkerColor(colors [icolor_Total]);
  hSystematicUncertainty_Total->SetLineColor  (colors [icolor_Total]);

  //Legend, if needed
  TLegend * leg = new TLegend(0.62, 0.6, 0.87, 0.8);
  leg->AddEntry(hSystematicUncertainty_CutVariation_PreBarlow,     "CutVariation",   "LP");
  leg->AddEntry(hSystematicUncertainty_Fit_SidebandVariation,     "SidebandVariation",  "LP");
  leg->AddEntry(hSystematicUncertainty_Fit_SignalExtractionType, "SignalExtractionType",    "LP" );
  leg->AddEntry(hSystematicUncertainty_Fit_BinningVariation, "BinningVariation",    "LP" );
  leg->AddEntry(hSystematicUncertainty_Total, "Total",    "LP" );
  // leg->SetFillColor(0);
  leg->SetTextSize(gStyle->GetTextSize()*0.8);
  leg->Draw("same");

  canvasSyst->SaveAs("Systematics_CutVariations_vsPt_"+NamePart[ipart]+".pdf");

  //Pre Barlow
  TCanvas *canvasSyst_preBarlow = new TCanvas ("canvasSyst_preBarlow", NamehistoInvMass[ipart], 1200, 800);
  canvasSyst_preBarlow->cd(0); //ATTENTION: +1 because canvas divide starts counting from 1
  TH1 *hFrame_preBarlow = canvasSyst_preBarlow->DrawFrame(pTbins[0],0,pTbins[nbinpT],1.2*hSystematicUncertainty_Total_PreBarlow->GetMaximum());
  hFrame_preBarlow->SetXTitle("#it{p}_{T} (GeV/#it{c})");
  hFrame_preBarlow->SetYTitle("uncertainty (ratio)");

  hSystematicUncertainty_CutVariation_PreBarlow->Draw("E,same");
  hSystematicUncertainty_Fit_SignalExtractionType_PreBarlow->Draw("E,same");
  hSystematicUncertainty_Fit_SidebandVariation_PreBarlow->Draw("E,same");
  hSystematicUncertainty_Fit_BinningVariation_PreBarlow->Draw("E,same");
  hSystematicUncertainty_Total_PreBarlow->Draw("E,same");

  hSystematicUncertainty_CutVariation_PreBarlow->SetMarkerStyle(markers[0]);
  hSystematicUncertainty_CutVariation_PreBarlow->SetMarkerColor(colors [icolor_CutVariation]);
  hSystematicUncertainty_CutVariation_PreBarlow->SetLineColor  (colors [icolor_CutVariation]);

  hSystematicUncertainty_Fit_SignalExtractionType_PreBarlow->SetMarkerStyle(markers[1]);
  hSystematicUncertainty_Fit_SignalExtractionType_PreBarlow->SetMarkerColor(colors [icolor_Fit_SignalExtractionType]);
  hSystematicUncertainty_Fit_SignalExtractionType_PreBarlow->SetLineColor  (colors [icolor_Fit_SignalExtractionType]);

  hSystematicUncertainty_Fit_SidebandVariation_PreBarlow->SetMarkerStyle(markers[1]);
  hSystematicUncertainty_Fit_SidebandVariation_PreBarlow->SetMarkerColor(colors [icolor_Fit_SidebandVariation]);
  hSystematicUncertainty_Fit_SidebandVariation_PreBarlow->SetLineColor  (colors [icolor_Fit_SidebandVariation]);

  hSystematicUncertainty_Fit_BinningVariation_PreBarlow->SetMarkerStyle(markers[1]);
  hSystematicUncertainty_Fit_BinningVariation_PreBarlow->SetMarkerColor(colors [icolor_Fit_BinningVariation]);
  hSystematicUncertainty_Fit_BinningVariation_PreBarlow->SetLineColor  (colors [icolor_Fit_BinningVariation]);

  hSystematicUncertainty_Total_PreBarlow->SetMarkerStyle(markers[0]);
  hSystematicUncertainty_Total_PreBarlow->SetMarkerColor(colors [icolor_Total]);
  hSystematicUncertainty_Total_PreBarlow->SetLineColor  (colors [icolor_Total]);

  //Legend, if needed
  TLegend * leg_preBarlow = new TLegend(0.62, 0.6, 0.87, 0.8);
  leg_preBarlow->AddEntry(hSystematicUncertainty_CutVariation_PreBarlow,     "CutVariation PreBarlow",   "LP");
  leg_preBarlow->AddEntry(hSystematicUncertainty_Fit_SidebandVariation_PreBarlow,     "SidebandVariation PreBarlow",  "LP");
  leg_preBarlow->AddEntry(hSystematicUncertainty_Fit_SignalExtractionType_PreBarlow, "SignalExtractionType PreBarlow",    "LP" );
  leg_preBarlow->AddEntry(hSystematicUncertainty_Fit_BinningVariation_PreBarlow, "BinningVariation PreBarlow",    "LP" );
  leg_preBarlow->AddEntry(hSystematicUncertainty_Total_PreBarlow, "Total PreBarlow",    "LP" );
  // leg->SetFillColor(0);
  leg_preBarlow->SetTextSize(gStyle->GetTextSize()*0.8);
  leg_preBarlow->Draw("same");

  canvasSyst_preBarlow->SaveAs("Systematics_CutVariations_PreBarlow_vsPt_"+NamePart[ipart]+".pdf");

  /////////////////////////////////Does Nothing - only there to respect template/////////////////////////////////////////////
  hstat = (TH1D*)hSystematicUncertainty_Total->Clone("hstat");
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

  // delete hCorrectedYield_vsPt;
  // delete hSystematicUncertainty;
  // delete hSystematicUncertainty_PreBarlow;
}

Double_t myLevyPtXsi(Double_t *pt, Double_t *par) { //https://twiki.cern.ch/twiki/bin/view/ALICE/PWG2SpectraTopical7TeVLambdas#Levy_fit

  Double_t lMass = 1.32171; //from PDG
  Double_t ldNdy = par[0];
  // Double_t l2pi = 2*TMath::Pi();
  Double_t lTemp = par[1];
  Double_t lPower = par[2]; // the n in the pp900GeV paper
  Double_t lEventCount= par[3]; // the n in the pp900GeV paper
  Double_t lBigCoef = ((lPower-1)*(lPower-2)) / (lPower*lTemp*(lPower*lTemp+lMass*(lPower-2)));
  Double_t lInPower = 1 + (TMath::Sqrt(pt[0]*pt[0]+lMass*lMass)-lMass) / (lPower*lTemp);

  return lEventCount * ldNdy * pt[0] * lBigCoef * TMath::Power(lInPower,(-1)*lPower);
}

void Get_FeedDown(TH1D * H1D_FeedDownCorrection, Int_t ipart, TFile *file_O2Analysis, Double_t* pTbins, Int_t nbinpT) {

  TH2D* H2D_Feeddown_numerator = (TH2D*)file_O2Analysis->Get("lambdakzero-analysis-mc/h"+NamePart[ipart]+"FeedDownMatrix");
  TH1D* H1D_Feeddown_denominator = (TH1D*)file_O2Analysis->Get("lambdakzero-particle-count-mc/hXsiMinusCount_PtDiff");
  // Rebinning of those two histograms
  TH1D* H1D_Feeddown_denominator_rebinned = (TH1D*)H1D_Feeddown_denominator->Rebin(nbinpT,"H1D_Feeddown_denominator_rebinned",pTbins);

  TH2D* H2D_Feeddown_numerator_rebinned = new TH2D("H2D_Feeddown_numerator_rebinned", "H2D_Feeddown_numerator_rebinned", nbinpT, pTbins, nbinpT, pTbins);

  TH1I* H1I_SelectedEventCount = (TH1I*)file_O2Analysis->Get("lambdakzero-analysis-mc/hSelectedEventCounter");

  for(Int_t i = 1; i <= nbinpT; i++){
    Int_t arrayID_i = i - 1;
    Int_t LowerBound_bin_i = H2D_Feeddown_numerator->GetXaxis()->FindBin(pTbins[arrayID_i]);
    Int_t UpperBound_bin_i = H2D_Feeddown_numerator->GetXaxis()->FindBin(pTbins[arrayID_i+1]);
    for(Int_t j = 1; j <= nbinpT; j++){
      Int_t arrayID_j = j - 1;
      Int_t LowerBound_bin_j = H2D_Feeddown_numerator->GetYaxis()->FindBin(pTbins[arrayID_j]);
      Int_t UpperBound_bin_j = H2D_Feeddown_numerator->GetYaxis()->FindBin(pTbins[arrayID_j+1]);
      H2D_Feeddown_numerator_rebinned->SetBinContent(i, j, H2D_Feeddown_numerator->Integral(LowerBound_bin_i,UpperBound_bin_i,LowerBound_bin_j,UpperBound_bin_j));
    }
  }

  H2D_Feeddown_numerator_rebinned->Sumw2();

  TH2D* hFij = new TH2D("hFij", "hFij", nbinpT, pTbins, nbinpT, pTbins);
  for(Int_t i = 1; i <= nbinpT; i++){
    for(Int_t j = 1; j <= nbinpT; j++){
      hFij->SetBinContent(i,j,H2D_Feeddown_numerator_rebinned->GetBinContent(i,j)*1./H1D_Feeddown_denominator_rebinned->GetBinContent(j));
      hFij->SetBinError(i,j,sqrt(hFij->GetBinContent(i,j)*(1-hFij->GetBinContent(i,j))*H1D_Feeddown_denominator_rebinned->GetBinContent(j))*1./H1D_Feeddown_denominator_rebinned->GetBinContent(j)); // this is effectively an efficiency: xsi- decays to (pi- + Lambda) 99.9% of the time (see PDG): so error is binomial: 1/n*sqrt(np(1-p))
    }
  }

  //Xsi yield from pp 900GeV paper fit function (Xsi- + AntiXsi+)
  //dNdy = 0.0101 ± 0.0020
  //T = 175 ± 50
  //n = 5.2 ± 2.3
  Float_t dNdy[2] = {3.05470e-04, 0.0020};
  Float_t lTemp[2] = {1.00055e-01,50};
  Float_t lPower[2] = {3.6,2.3};
  TF1* LevisTsalisFunction_Xsipp900GeV = new TF1("LevisTsalisFunction_Xsipp900GeV", myLevyPtXsi, 0, pTbins[nbinpT], 4);

  // Get generated Xsi distribution
  // to get more statistics, get Xsi- + Xsi+ given the ratio Xsi-/Xsi+ is very close to 1, within a few percents; DAMN only true for Xsi-/Xsi0
  TH1D* H1D_Xsi_PtDistrib;
  if (ipart > 2) {
    cout << "ipart used is not covered by Get_FeedDown()" << endl;
  }
  H1D_Xsi_PtDistrib = (TH1D*)file_O2Analysis->Get("lambdakzero-particle-count-mc/h"+NamePart_multiStrange[ipart-1]+"Count_PtDiff");


  H1D_Xsi_PtDistrib->Sumw2();
  H1D_Xsi_PtDistrib->Rebin(5.);
  // H1D_XsiMinus_PtDistrib->Scale(1./H1I_SelectedEventCount->GetEntries()); //scale to per event yield; done outside function because fit likes it better if not done before, for errors.

  LevisTsalisFunction_Xsipp900GeV->SetParameter(0, dNdy[0]);
  LevisTsalisFunction_Xsipp900GeV->SetParameter(1, lTemp[0]);
  LevisTsalisFunction_Xsipp900GeV->SetParameter(2, lPower[0]);
  LevisTsalisFunction_Xsipp900GeV->FixParameter(3,H1I_SelectedEventCount->GetEntries());

  // LevisTsalisFunction_Xsipp900GeV->SetParError(0, dNdy[1]);
  // LevisTsalisFunction_Xsipp900GeV->SetParError(1, lTemp[1]);
  // LevisTsalisFunction_Xsipp900GeV->SetParError(2, lPower[1]);

  // LevisTsalisFunction_Xsipp900GeV->SetParLimits(0, 0.0001, 0.1);
  // LevisTsalisFunction_Xsipp900GeV->SetParLimits(1, 0, 900);
  // LevisTsalisFunction_Xsipp900GeV->SetParLimits(2, 0, 20);



  // TFitResultPtr fFitResult = H1D_XsiMinus_PtDistrib->Fit(LevisTsalisFunction_Xsipp900GeV,"S+WLR"); //
  TFitResultPtr fFitResult = H1D_Xsi_PtDistrib->Fit(LevisTsalisFunction_Xsipp900GeV,"S+L0R"); //

  TCanvas *canvasFeedDown = new TCanvas ("canvasFeedDown", NamehistoInvMass[ipart], 1200, 800);
  canvasFeedDown->cd(0); //ATTENTION: +1 because canvas divide starts counting from 1
  TH1 *hFrameFeeddown = canvasFeedDown->DrawFrame(0.,0,pTbins[nbinpT],1.2*H1D_Xsi_PtDistrib->GetMaximum());
  hFrameFeeddown->SetXTitle("#it{p}_{T} (GeV/#it{c})");
  hFrameFeeddown->SetYTitle("XsiMinusDistribution");
  H1D_Xsi_PtDistrib->Draw("E, same");
  LevisTsalisFunction_Xsipp900GeV->Draw("E, same");
  canvasFeedDown->SaveAs("FeedDownLevyFit.pdf");

  Double_t * param = LevisTsalisFunction_Xsipp900GeV->GetParameters();
  TMatrixDSym covMatrix = fFitResult->GetCovarianceMatrix();

  for(Int_t i = 1; i <= nbinpT; i++){
    Double_t FeedDownCorrection_i = 0;
    Double_t FeedDownCorrection_i_StatUncertainty = 0;
    cout << "*---------------- i: " << i << " ----------------*" << endl;
    for(Int_t j = 1; j <= nbinpT; j++){
      Double_t CorrectedYield_Xsi_j = LevisTsalisFunction_Xsipp900GeV->Integral(pTbins[j-1],pTbins[j])*1./((pTbins[j]-pTbins[j-1]));
      // Double_t CorrectedYield_Xsi_j = H1D_Xsi_PtDistrib->Integral(pTbins[j-1],pTbins[j]); // temporary
      cout << "j: " << j << ", CorrectedYield_Xsi_j = " << CorrectedYield_Xsi_j << endl; //surprisingly is most high at j=9 with 400  -> because pT bins increase in length

      Double_t CorrectedYield_Xsi_j_uncertainty = LevisTsalisFunction_Xsipp900GeV->IntegralError(pTbins[j-1],pTbins[j], param, covMatrix.GetMatrixArray());
      // Double_t CorrectedYield_Xsi_j_uncertainty = 0; test
      FeedDownCorrection_i = FeedDownCorrection_i + hFij->GetBinContent(i,j)*CorrectedYield_Xsi_j;
      FeedDownCorrection_i_StatUncertainty = sqrt(FeedDownCorrection_i_StatUncertainty*FeedDownCorrection_i_StatUncertainty + hFij->GetBinError(i,j)*hFij->GetBinError(i,j) + CorrectedYield_Xsi_j_uncertainty*CorrectedYield_Xsi_j_uncertainty); //product and sums
      cout << "      numerator = " << H2D_Feeddown_numerator_rebinned->GetBinContent(i,j) << endl;
      cout << "      denominator = " << H1D_Feeddown_denominator_rebinned->GetBinContent(j) << endl;
      cout << "      hFij = " << hFij->GetBinContent(i,j) << endl;
      cout << "      hFij Error = " << hFij->GetBinError(i,j) << endl;
    }
    H1D_FeedDownCorrection->SetBinContent(i, FeedDownCorrection_i);
    H1D_FeedDownCorrection->SetBinError(i, FeedDownCorrection_i_StatUncertainty);
    cout << "           FeedDownCorrection_i = " << FeedDownCorrection_i << endl;
    cout << "                          error = " << FeedDownCorrection_i_StatUncertainty << endl;
  }

  H1D_FeedDownCorrection->Scale(2.); // need to multiply contribution by two given Xsi0 is also contributing to Lambda feeddown and has not been put in the Fij matrix
  //                                    ideally we could do the same with Xsi0 as done here for Xsi-/+ but Xsi0 can't be measured given it decays into pi0 which is neutral (and Lambda obviously that we measure but isn't enough on its own to reconstruct the decay and that polutes our primary lambda result)
  H1D_FeedDownCorrection->Scale(1./H1I_SelectedEventCount->GetEntries()); //scale to per event yield; 
  delete hFij;
}

void Draw_FeedDownMatrix(TH1D* &hstat, TH1D* &hsyst, TH1D* &hsystCorr, Int_t ipart, TFile *file_O2Analysis, TString* SaveAs_Title, TString* &texXtitle, TString* &texYtitle, Double_t* pTbins, Int_t nbinpT) {

  TH2D* H2D_Feeddown_numerator = (TH2D*)file_O2Analysis->Get("lambdakzero-analysis-mc/h"+NamePart[ipart]+"FeedDownMatrix");
  
  TH2D* H2D_Feeddown_numerator_rebinned = new TH2D("H2D_Feeddown_numerator_rebinned", "H2D_Feeddown_numerator_rebinned", nbinpT, pTbins, nbinpT, pTbins);

  for(Int_t i = 1; i <= nbinpT; i++){
    Int_t arrayID_i = i - 1;
    Int_t LowerBound_bin_i = H2D_Feeddown_numerator->GetXaxis()->FindBin(pTbins[arrayID_i]);
    Int_t UpperBound_bin_i = H2D_Feeddown_numerator->GetXaxis()->FindBin(pTbins[arrayID_i+1]);
    for(Int_t j = 1; j <= nbinpT; j++){
      Int_t arrayID_j = j - 1;
      Int_t LowerBound_bin_j = H2D_Feeddown_numerator->GetYaxis()->FindBin(pTbins[arrayID_j]);
      Int_t UpperBound_bin_j = H2D_Feeddown_numerator->GetYaxis()->FindBin(pTbins[arrayID_j+1]);
      H2D_Feeddown_numerator_rebinned->SetBinContent(i, j, H2D_Feeddown_numerator->Integral(LowerBound_bin_i,UpperBound_bin_i,LowerBound_bin_j,UpperBound_bin_j));
    }
  }

  TCanvas *canvasFeedDown = new TCanvas ("canvasFeedDown", NamehistoInvMass[ipart], 1200, 800);
  canvasFeedDown->cd(0); //ATTENTION: +1 because canvas divide starts counting from 1
  // TH1 *hFrame = canvasFeedDown->DrawFrame(pTbins[0],0,pTbins[nbinpT],1.2*hSystematicUncertainty_Total->GetMaximum());
  // hFrame->SetXTitle("#it{p}_{T} (GeV/#it{c})");
  // hFrame->SetYTitle("uncertainty (ratio)");

  H2D_Feeddown_numerator_rebinned->Draw("COLZ");
  canvasFeedDown->SaveAs("FeedDown_Matrix_"+NamePart[ipart]+".pdf");

  /////////////////////////////////Does Nothing - only there to respect template/////////////////////////////////////////////
  hstat = (TH1D*)H2D_Feeddown_numerator_rebinned->Clone("hstat");
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

  // delete hCorrectedYield_vsPt;
  // delete hSystematicUncertainty;
  // delete hSystematicUncertainty_PreBarlow;
}


void Draw_XsiPlus_MassPlot(TH1D* &hstat, TH1D* &hsyst, TH1D* &hsystCorr, Int_t ipart, TFile *file_O2Analysis, TString* SaveAs_Title, TString* &texXtitle, TString* &texYtitle, Double_t* pTbins, Int_t nbinpT) {

  TH1D* H1D_XsiPlus_PtDistrib = (TH1D*)file_O2Analysis->Get("lambdakzero-particle-count-mc/hXsiPlusCount_PtDiff");
  TH1D* H1D_XsiPlus_PtDistrib_rebinned = (TH1D*)H1D_XsiPlus_PtDistrib->Rebin(nbinpT,"H1D_XsiPlus_PtDistrib_rebinned",pTbins);



  TCanvas *canvasFeedDown = new TCanvas ("canvasFeedDown", NamehistoInvMass[ipart], 1200, 800);
  canvasFeedDown->cd(0); //ATTENTION: +1 because canvas divide starts counting from 1
  // TH1 *hFrame = canvasFeedDown->DrawFrame(pTbins[0],0,pTbins[nbinpT],1.2*hSystematicUncertainty_Total->GetMaximum());
  // hFrame->SetXTitle("#it{p}_{T} (GeV/#it{c})");
  // hFrame->SetYTitle("uncertainty (ratio)");

  H1D_XsiPlus_PtDistrib_rebinned->Draw("");
  canvasFeedDown->SaveAs("XsiPlus_InvMassPlot.pdf");

  /////////////////////////////////Does Nothing - only there to respect template/////////////////////////////////////////////
  hstat = (TH1D*)H1D_XsiPlus_PtDistrib_rebinned->Clone("hstat");
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

  // delete hCorrectedYield_vsPt;
  // delete hSystematicUncertainty;
  // delete hSystematicUncertainty_PreBarlow;
}

void Draw_FeedDown(TH1D* &hstat, TH1D* &hsyst, TH1D* &hsystCorr, Int_t ipart, TFile *file_O2Analysis, TString* SaveAs_Title, TString* &texXtitle, TString* &texYtitle, Double_t* pTbins, Int_t nbinpT) {

  TH1D* H1D_XsiPlus_InvMassPlot = (TH1D*)file_O2Analysis->Get("lambdakzero-particle-count-mc/hXsiPlusCount_PtDiff");
  TH1D* H1D_XsiPlus_InvMassPlot_rebinned = (TH1D*)H1D_XsiPlus_InvMassPlot->Rebin(nbinpT,"H1D_XsiPlus_InvMassPlot_rebinned",pTbins);

  //Event count
  TH1I* H1I_SelectedEventCount = (TH1I*)file_O2Analysis->Get("lambdakzero-analysis-mc/hSelectedEventCounter");

  //Feeddown
  TH1D* H1D_Feeddown_Correction = new TH1D("H1D_Feeddown_Correction", "H1D_Feeddown_Correction", nbinpT, pTbins);
  Get_FeedDown(H1D_Feeddown_Correction, ipart, file_O2Analysis, pTbins, nbinpT);

  //Raw yield
  TH3D* H3D_detectedV0s_DataLike = (TH3D*)file_O2Analysis->Get("lambdakzero-analysis-mc/h3dMass"+NamePart[ipart]);
  TH1D *hRawYield_vsPt;
  Int_t InvMassRebinFactor = InvMassRebinFactor_standard[ipart]; //used if signalextractiontype = 1
  Get_RawYield_vsPt(hRawYield_vsPt, H3D_detectedV0s_DataLike, H1I_SelectedEventCount, ipart, pTbins, nbinpT, 0, -1, 0, SideBandSizeMultiplierModifier_array[0], SignalExtractionType_default, InvMassRebinFactor); ///real line that need to take the place of the one above once it's been written

  TH1D* H1D_FractionRawYieldRemoved = new TH1D("H1D_FractionRawYieldRemoved", "H1D_FractionRawYieldRemoved", nbinpT, pTbins);
  H1D_FractionRawYieldRemoved->Divide(H1D_Feeddown_Correction,hRawYield_vsPt);

  TCanvas *canvasFeedDown = new TCanvas ("canvasFeedDown", NamehistoInvMass[ipart], 1200, 800);
  canvasFeedDown->cd(0); //ATTENTION: +1 because canvas divide starts counting from 1
  TH1 *hFrame = canvasFeedDown->DrawFrame(pTbins[0],0,pTbins[nbinpT],1.2*H1D_FractionRawYieldRemoved->GetMaximum());
  hFrame->SetXTitle("#it{p}_{T} (GeV/#it{c})");
  hFrame->SetYTitle("FractionRawYieldRemoved from feeddown");
  H1D_FractionRawYieldRemoved->Draw("E,same");
  canvasFeedDown->SaveAs("FractionRawYieldRemoved.pdf");

  /////////////////////////////////Does Nothing - only there to respect template/////////////////////////////////////////////
  hstat = (TH1D*)H1D_XsiPlus_InvMassPlot_rebinned->Clone("hstat");
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

  // delete hCorrectedYield_vsPt;
  // delete hSystematicUncertainty;
  // delete hSystematicUncertainty_PreBarlow;
}



void CorrectedYield_withSystematics(TH1D* &hstat, TH1D* &hsyst, TH1D* &hsystCorr, Int_t ipart, TFile** file_O2Analysis_CutVariation_array, Int_t id_referenceAnalysis, Int_t nCutsIterations, TString* SaveAs_Title, TString* &texXtitle, TString* &texYtitle, Double_t* pTbins, Int_t nbinpT) {
  Int_t aaaa = nCutsIterations; //to remove warning

  //HEP references for comparison
  TFile* file_HEPAnalysis = new TFile("HEPData-ins881474-v1-root.root");
  TH1D* hPtDiff_RawSpectrum_HEP[numPart]; //= (TH1D*)file_HEPAnalysis->Get("Table 1/Hist1D_y1");
  TH1D* hPtDiff_RawSpectrum_HEP_statErrors[numPart]; // = (TH1D*)file_HEPAnalysis->Get("Table 1/Hist1D_y1_e1");
  hPtDiff_RawSpectrum_HEP[0] = (TH1D*)file_HEPAnalysis->Get("Table 1/Hist1D_y1"); //K0S
  hPtDiff_RawSpectrum_HEP_statErrors[0] = (TH1D*)file_HEPAnalysis->Get("Table 1/Hist1D_y1_e1");
  hPtDiff_RawSpectrum_HEP[1] = (TH1D*)file_HEPAnalysis->Get("Table 2/Hist1D_y1"); //Lambda
  hPtDiff_RawSpectrum_HEP_statErrors[1] = (TH1D*)file_HEPAnalysis->Get("Table 2/Hist1D_y1_e1");
  hPtDiff_RawSpectrum_HEP[2] = (TH1D*)file_HEPAnalysis->Get("Table 3/Hist1D_y1"); //AntiLambda
  hPtDiff_RawSpectrum_HEP_statErrors[2] = (TH1D*)file_HEPAnalysis->Get("Table 3/Hist1D_y1_e1");

  if (ipart>2) {
    cout << "ipart>2, HEP data downloaded and coded only goes to 2" << endl;
  } 
  TH1D* hPtDiff_RawSpectrum_HEP_rebinned = new TH1D("hPtDiff_RawSpectrum_HEP_rebinned", "hPtDiff_RawSpectrum_HEP_rebinned", nbinpT, pTbins); //in case we want to discard one or more bins like the last one for lambda between pT 3.0 and 3.5 GeV
  for(int ibinPt = 1; ibinPt <= nbinpT; ibinPt++){
    hPtDiff_RawSpectrum_HEP_rebinned->SetBinContent(ibinPt, hPtDiff_RawSpectrum_HEP[ipart]->GetBinContent(ibinPt));
    hPtDiff_RawSpectrum_HEP_rebinned->SetBinError(ibinPt, hPtDiff_RawSpectrum_HEP_statErrors[ipart]->GetBinContent(ibinPt));
  }

  TH3D* H3D_detectedV0s_TruePt = (TH3D*)file_O2Analysis_CutVariation_array[id_referenceAnalysis]->Get("lambdakzero-analysis-mc/h3dMass"+NamePart[ipart]+"_MC_truePt");
  TH3D* H3D_detectedV0s_DataLike = (TH3D*)file_O2Analysis_CutVariation_array[id_referenceAnalysis]->Get("lambdakzero-analysis-mc/h3dMass"+NamePart[ipart]);
  TH1D* TrueV0PtSpectrum = (TH1D*)file_O2Analysis_CutVariation_array[id_referenceAnalysis]->Get("lambdakzero-particle-count-mc/h"+NamePart[ipart]+"Count_PtDiff");
  TH1I* H1I_SelectedEventCount = (TH1I*)file_O2Analysis_CutVariation_array[id_referenceAnalysis]->Get("lambdakzero-analysis-mc/hSelectedEventCounter");

  TH1D* hSystematicUncertainty_Total = new TH1D("hSystematicUncertainty", "hSystematicUncertainty", nbinpT, pTbins);
  TH1D* hSystematicUncertainty_Total_PreBarlow = new TH1D("hSystematicUncertainty_PreBarlow", "hSystematicUncertainty_PreBarlow", nbinpT, pTbins);
  hSystematicUncertainty_Total->Sumw2();
  hSystematicUncertainty_Total_PreBarlow->Sumw2();

  TH1D* hSystematicUncertainty_CutVariation = new TH1D("hSystematicUncertainty_CutVariation", "hSystematicUncertainty_CutVariation", nbinpT, pTbins);
  TH1D* hSystematicUncertainty_CutVariation_PreBarlow = new TH1D("hSystematicUncertainty_CutVariation_PreBarlow", "hSystematicUncertainty_CutVariation_PreBarlow", nbinpT, pTbins);
  TH1D* hSystematicUncertainty_Fit_SidebandVariation = new TH1D("hSystematicUncertainty_Fit_SidebandVariation", "hSystematicUncertainty_Fit_SidebandVariation", nbinpT, pTbins);
  TH1D* hSystematicUncertainty_Fit_SidebandVariation_PreBarlow = new TH1D("hSystematicUncertainty_Fit_SidebandVariation_PreBarlow", "hSystematicUncertainty_Fit_SidebandVariation_PreBarlow", nbinpT, pTbins);
  TH1D* hSystematicUncertainty_Fit_SignalExtractionType = new TH1D("hSystematicUncertainty_Fit_SignalExtractionType", "hSystematicUncertainty_Fit_SignalExtractionType", nbinpT, pTbins);
  TH1D* hSystematicUncertainty_Fit_SignalExtractionType_PreBarlow = new TH1D("hSystematicUncertainty_Fit_SignalExtractionType_PreBarlow", "hSystematicUncertainty_Fit_SignalExtractionType_PreBarlow", nbinpT, pTbins);
  TH1D* hSystematicUncertainty_Fit_BinningVariation = new TH1D("hSystematicUncertainty_Fit_BinningVariation", "hSystematicUncertainty_Fit_BinningVariation", nbinpT, pTbins);
  TH1D* hSystematicUncertainty_Fit_BinningVariation_PreBarlow = new TH1D("hSystematicUncertainty_Fit_BinningVariation_PreBarlow", "hSystematicUncertainty_Fit_BinningVariation_PreBarlow", nbinpT, pTbins);

  hSystematicUncertainty_CutVariation->Sumw2();
  hSystematicUncertainty_CutVariation_PreBarlow->Sumw2();
  hSystematicUncertainty_Fit_SidebandVariation->Sumw2();
  hSystematicUncertainty_Fit_SidebandVariation_PreBarlow->Sumw2();
  hSystematicUncertainty_Fit_SignalExtractionType->Sumw2();
  hSystematicUncertainty_Fit_SignalExtractionType_PreBarlow->Sumw2();
  hSystematicUncertainty_Fit_BinningVariation->Sumw2();
  hSystematicUncertainty_Fit_BinningVariation_PreBarlow->Sumw2();

  TH1D *hMcEfficiency_vsPt;
  TH1D *hRawYield_vsPt;
  TH1D *H1D_Feeddown_Correction = new TH1D("H1D_Feeddown_Correction", "H1D_Feeddown_Correction", nbinpT, pTbins);
  Int_t InvMassRebinFactor = InvMassRebinFactor_standard[ipart];
  

  Get_McEfficiency_vsPt(hMcEfficiency_vsPt, H3D_detectedV0s_TruePt, TrueV0PtSpectrum, ipart, pTbins, nbinpT, 0, -1);
  Get_RawYield_vsPt(hRawYield_vsPt, H3D_detectedV0s_DataLike, H1I_SelectedEventCount, ipart, pTbins, nbinpT, 0, -1, 0, SideBandSizeMultiplierModifier_array[0], SignalExtractionType_default, InvMassRebinFactor); ///real line that need to take the place of the one above once it's been written
  Get_Systematics_OneCut(hSystematicUncertainty_CutVariation, hSystematicUncertainty_CutVariation_PreBarlow, ipart, file_O2Analysis_CutVariation_array, id_referenceAnalysis, nCutsIterations, pTbins, nbinpT, InvMassRebinFactor);
  Get_Systematics_Fit_SideBandVariation(hSystematicUncertainty_Fit_SidebandVariation, hSystematicUncertainty_Fit_SidebandVariation_PreBarlow, ipart, H3D_detectedV0s_DataLike, H1I_SelectedEventCount, H3D_detectedV0s_TruePt, TrueV0PtSpectrum, pTbins, nbinpT, InvMassRebinFactor);
  Get_Systematics_Fit_SignalExtractionType(hSystematicUncertainty_Fit_SignalExtractionType, hSystematicUncertainty_Fit_SignalExtractionType_PreBarlow, ipart, H3D_detectedV0s_DataLike, H1I_SelectedEventCount, H3D_detectedV0s_TruePt, TrueV0PtSpectrum, pTbins, nbinpT, InvMassRebinFactor); //very large systematics as is
  Get_Systematics_Fit_Binning(hSystematicUncertainty_Fit_BinningVariation, hSystematicUncertainty_Fit_BinningVariation_PreBarlow, ipart, H3D_detectedV0s_DataLike, H1I_SelectedEventCount, H3D_detectedV0s_TruePt, TrueV0PtSpectrum, pTbins, nbinpT);
  // Get_Systematics_OneCut_McTrueAndMcDatalikeDifferentCuts(hSystematicUncertainty, hSystematicUncertainty_PreBarlow, ipart, file_O2Analysis_CutVariation_Datalike_array, file_O2Analysis_CutVariation_True_array, id_referenceAnalysis, nCutsIterations, pTbins, nbinpT, InvMassRebinFactor);

  hSystematicUncertainty_Total->Add(hSystematicUncertainty_CutVariation);
  hSystematicUncertainty_Total->Add(hSystematicUncertainty_Fit_SidebandVariation);
  hSystematicUncertainty_Total->Add(hSystematicUncertainty_Fit_SignalExtractionType);
  hSystematicUncertainty_Total->Add(hSystematicUncertainty_Fit_BinningVariation);

  // hSystematicUncertainty_Total->Add(hSystematicUncertainty_CutVariation_PreBarlow);
  // hSystematicUncertainty_Total->Add(hSystematicUncertainty_Fit_SidebandVariation_PreBarlow);
  // hSystematicUncertainty_Total->Add(hSystematicUncertainty_Fit_SignalExtractionType_PreBarlow);
  // hSystematicUncertainty_Total->Add(hSystematicUncertainty_Fit_BinningVariation_PreBarlow);

  // Get_Systematics_Full(); // gather all systematics; eventually

  // TH1D *hCorrectedYield_vsPt;
  // hCorrectedYield_vsPt = (TH1D*)hRawYield_vsPt->Clone("hCorrectedYield_vsPt");
  // hCorrectedYield_vsPt->Reset("M");
  TH1D* hCorrectedYield_vsPt = new TH1D("hCorrectedYield_vsPt", "hCorrectedYield_vsPt",nbinpT,pTbins);
  hCorrectedYield_vsPt->Sumw2();
  
  // Lambda and AntiLambda Feeddown
  if (ipart == 1 || ipart == 2) {
    Get_FeedDown(H1D_Feeddown_Correction, ipart, file_O2Analysis_CutVariation_array[id_referenceAnalysis], pTbins, nbinpT);
    H1D_Feeddown_Correction->Scale(-1.);
    hRawYield_vsPt->Add(H1D_Feeddown_Correction);
    // cout << "test00000" << endl;
  }

  hCorrectedYield_vsPt->Divide(hRawYield_vsPt,hMcEfficiency_vsPt, 1., 1., "b");

  hstat = (TH1D*)hCorrectedYield_vsPt->Clone("hstat");
  hsyst = (TH1D*)hCorrectedYield_vsPt->Clone("hsyst");

  TH1* hCorrectedYield_vsPt_syst = (TH1D*)hCorrectedYield_vsPt->Clone("hCorrectedYield_vsPt_syst");
  for(int ibinPt = 1; ibinPt <= nbinpT; ibinPt++){
    hsyst->SetBinError(ibinPt,hSystematicUncertainty_Total->GetBinContent(ibinPt));
    // hsyst->SetBinError(ibinPt,hSystematicUncertainty_PreBarlow->GetBinContent(ibinPt));
    // hsyst->SetBinError(ibinPt,0.00001);
    cout << "ibinPt = "<< ibinPt << ", rawYield = " << hRawYield_vsPt->GetBinContent(ibinPt) << ", Efficiency = " << hMcEfficiency_vsPt->GetBinContent(ibinPt) << endl;
    hCorrectedYield_vsPt_syst->SetBinError(ibinPt, hSystematicUncertainty_Total->GetBinContent(ibinPt));
  }

  /////////// Draw Efficiency /////////////
  TString* pdfName0 = new TString("Efficiency_"+NamePart[ipart]);;
  // *pdfName0 += "Efficiency_"+NamePart[ipart];

  TH1* histograms_collection_0[1] = {hMcEfficiency_vsPt};
  TString* pdfLegend_collection_0[1] = {new TString("Efficiency "+NamePart_Latex[ipart])};
  Int_t ratioOption0 = 0;
  texXtitle = texPtX;
  texYtitle = texEfficiency;
  Draw_TH1_Histograms(histograms_collection_0, pdfLegend_collection_0, 1, pdfName0, texXtitle, texYtitle, pTbins, nbinpT, ratioOption0);


  /////////// Draw Corrected yield comparison O2 vs HEPdata /////////////
  TString* pdfName1 = new TString("CorrectewdYields_O2_vs_HEPdata_vs_MCgenYield_"+NamePart[ipart]);;
  // *pdfName1 += ;
  // O2 MC true yield
  // TH1I* H1I_TotalEventCount = (TH1I*)file_O2Analysis->Get("lambdakzero-builder/hEventCounter");
  TH1I* H1I_TotalEventCount = (TH1I*)file_O2Analysis_CutVariation_array[id_referenceAnalysis]->Get("lambdakzero-particle-count-mc/hSelAndRecoMcCollCounter");

  TrueV0PtSpectrum->Sumw2();
  TH1D* TrueV0PtSpectrum_rebinned = (TH1D*)TrueV0PtSpectrum->Rebin(nbinpT,"TrueV0PtSpectrum_rebinned",pTbins);
  TH1D* TrueV0PtSpectrum_rebinned_syst = (TH1D*)TrueV0PtSpectrum->Rebin(nbinpT,"TrueV0PtSpectrum_rebinned_syst",pTbins);
  for (Int_t ibinPt = 1; ibinPt <= nbinpT; ibinPt++)
  {
    Double_t dpT = TrueV0PtSpectrum_rebinned->GetXaxis()->GetBinWidth(ibinPt);
    Double_t drapidity = 1.5; //1.5
    Double_t d2N_dpTdy = TrueV0PtSpectrum_rebinned->GetBinContent(ibinPt) *1./dpT*1./drapidity;
    Double_t d2N_dpTdy_error = sqrt(TrueV0PtSpectrum_rebinned->GetBinContent(ibinPt)) * sqrt(1./dpT*1./drapidity);
    TrueV0PtSpectrum_rebinned->SetBinContent(ibinPt, d2N_dpTdy);
    TrueV0PtSpectrum_rebinned->SetBinError(ibinPt, d2N_dpTdy_error);
    TrueV0PtSpectrum_rebinned_syst->SetBinContent(ibinPt, d2N_dpTdy);
    TrueV0PtSpectrum_rebinned_syst->SetBinError(ibinPt, 0.);
  }
  TrueV0PtSpectrum_rebinned->Scale(1./H1I_TotalEventCount->GetEntries());
  TrueV0PtSpectrum_rebinned_syst->Scale(1./H1I_TotalEventCount->GetEntries());

  // TH1* histograms_collection_1[4] = {hCorrectedYield_vsPt_syst, hCorrectedYield_vsPt, hPtDiff_RawSpectrum_HEP_rebinned, TrueV0PtSpectrum_rebinned};
  // TString* pdfLegend_collection_1[4] = {new TString("hCorrectedYield_vsPt_syst"), new TString("Corrected Yield"), new TString("published pp #sqrt{#it{s}} = 900 GeV"),new TString("True MC")};
  // Int_t ratioOption1 = 0;
  // texXtitle = texPtX;
  // texYtitle = texptDifferentialYield;
  // Draw_TH1_Histograms_withOneHistSystematics(histograms_collection_1, pdfLegend_collection_1, 4, pdfName1, texXtitle, texYtitle, pTbins, nbinpT, ratioOption1);
  TH1* histograms_collection_1[3] = {hCorrectedYield_vsPt_syst, hCorrectedYield_vsPt, TrueV0PtSpectrum_rebinned};
  TString* pdfLegend_collection_1[3] = {new TString("hCorrectedYield_vsPt_syst"), new TString("Corrected Yield"),new TString("True MC")};
  Int_t ratioOption1 = 0;
  texXtitle = texPtX;
  texYtitle = texptDifferentialYield;
  Draw_TH1_Histograms_withOneHistSystematics(histograms_collection_1, pdfLegend_collection_1, 3, pdfName1, texXtitle, texYtitle, pTbins, nbinpT, ratioOption1);

  ////////////////// Ratio_CorrectewdYield_to_MCgenYield and HEP //////////////////
  TString* pdfName2 = new TString("Ratio_CorrectewdYield_to_MCgenYield_"+NamePart[ipart]);;
  // *pdfName2 += "Ratio_CorrectewdYield_to_MCgenYield";

  TH1* Ratio_CorrectewdYield_to_MCgenYield;
  TH1* Ratio_CorrectewdYield_to_MCgenYield_syst;
  TH1* Ratio_HEP_to_CorrectewdYield;

  Ratio_CorrectewdYield_to_MCgenYield = (TH1D*)TrueV0PtSpectrum_rebinned->Clone("Ratio_CorrectewdYield_to_MCgenYield");
  Ratio_CorrectewdYield_to_MCgenYield->Reset("M");
  Ratio_CorrectewdYield_to_MCgenYield->Sumw2();
  Ratio_CorrectewdYield_to_MCgenYield->Divide(hCorrectedYield_vsPt,TrueV0PtSpectrum_rebinned);
  Ratio_CorrectewdYield_to_MCgenYield_syst = (TH1D*)TrueV0PtSpectrum_rebinned->Clone("Ratio_CorrectewdYield_to_MCgenYield_syst");
  Ratio_CorrectewdYield_to_MCgenYield_syst->Reset("M");
  Ratio_CorrectewdYield_to_MCgenYield_syst->Sumw2();
  Ratio_CorrectewdYield_to_MCgenYield_syst->Divide(hCorrectedYield_vsPt_syst,TrueV0PtSpectrum_rebinned_syst);

  Ratio_HEP_to_CorrectewdYield = (TH1D*)TrueV0PtSpectrum_rebinned->Clone("Ratio_HEP_to_CorrectewdYield");
  Ratio_HEP_to_CorrectewdYield->Reset("M");
  Ratio_HEP_to_CorrectewdYield->Sumw2();
  Ratio_HEP_to_CorrectewdYield->Divide(hPtDiff_RawSpectrum_HEP_rebinned,hCorrectedYield_vsPt);
  // Ratio_HEP_to_CorrectewdYield->Divide(hPtDiff_RawSpectrum_HEP[ipart],TrueV0PtSpectrum_rebinned);

  // TH1* histograms_collection_2[3] = {Ratio_CorrectewdYield_to_MCgenYield_syst, Ratio_CorrectewdYield_to_MCgenYield, Ratio_HEP_to_CorrectewdYield};
  // // TString* pdfLegend_collection_2[3] = {new TString("Ratio_CorrectewdYield_to_MCgenYield_syst"), new TString("#frac{corrected yield}{true MC yield}"), new TString("#frac{yield published pp #sqrt{#it{s}} = 900 GeV}{corrected yield}")};
  // TString* pdfLegend_collection_2[3] = {new TString("Ratio_CorrectewdYield_to_MCgenYield_syst"), new TString("CorrectedYield/TrueMC"), new TString("PublishedYield/CorrectedYield")};
  TH1* histograms_collection_2[2] = {Ratio_CorrectewdYield_to_MCgenYield_syst, Ratio_CorrectewdYield_to_MCgenYield};
  TString* pdfLegend_collection_2[2] = {new TString("Ratio_CorrectewdYield_to_MCgenYield_syst"), new TString("CorrectedYield/TrueMC")};
  // TH1* histograms_collection_2[1] = {Ratio_HEP_to_CorrectewdYield};
  // TString* pdfLegend_collection_2[1] = {new TString("Ratio_HEP_to_CorrectewdYield")};
  Int_t ratioOption2 = 1; //to draw line at ratio=1
  texXtitle = texPtX;
  texYtitle = texRatio;
  // Draw_TH1_Histograms_withOneHistSystematics(histograms_collection_2, pdfLegend_collection_2, 3, pdfName2, texXtitle, texYtitle, pTbins, nbinpT, ratioOption2);
  Draw_TH1_Histograms_withOneHistSystematics_ratio(histograms_collection_2, pdfLegend_collection_2, 2, pdfName2, texXtitle, texYtitle, pTbins, nbinpT, ratioOption2);

  /////////////////////////////////Does Nothing - only there to respect template/////////////////////////////////////////////
  //////Error Bars///////
  hsystCorr = new TH1D("hsystCorr", "hsystCorr", 100, 0, 10);
  // hstat     ->Sumw2(); 
  // hsyst     ->Sumw2();
  hsystCorr ->Sumw2(); 
  Int_t nbinx = hstat->GetNbinsX();
  
  for(Int_t ibinx = 1; ibinx <= nbinx; ibinx++){
    hsystCorr->SetBinContent (ibinx, hstat->GetBinContent(ibinx));
    hsystCorr->SetBinError   (ibinx, hstat->GetBinContent(ibinx)*0.0);
  }
  *SaveAs_Title += "CorrectedYield_withSystematics";
  texXtitle = texPtX;
  texYtitle = texptDifferentialYield;
}

void Draw_TH1_Histograms(TH1** histograms_collection, TString** pdfLegend_collection, Int_t collectionSize, TString* pdfName, TString* &texXtitle, TString* &texYtitle, Double_t* bins, Int_t nbins, Int_t ratioOption) {
  
  Float_t MaximumY = 0;
  for (Int_t i = 0; i < collectionSize; i++) {
    Float_t maxTemporary = histograms_collection[i]->GetMaximum();
    if (maxTemporary > MaximumY) {
      MaximumY = maxTemporary;
    }
  }

  ////plots
  TCanvas *canvas = new TCanvas ("canvas"+*pdfName, NamehistoInvMass[ipart], 1200, 800);
  canvas->cd(0); //ATTENTION: +1 because canvas divide starts counting from 1
  TH1 *hFrame = canvas->DrawFrame(bins[0],0,bins[nbins],1.4*MaximumY);
  hFrame->SetXTitle(texXtitle->Data());
  hFrame->SetYTitle(texYtitle->Data());
  TLegend * leg = new TLegend(0.62, 0.77, 0.87, 0.87);

  //draw line at ratio=1
  if (ratioOption == 1) {
    // TLine *l=new TLine(histograms_collection[0]->GetXaxis()->GetBinLowEdge(0), 1.0, histograms_collection[0]->GetXaxis()->GetBinUpEdge(nbins), 1.0);
    TLine *l=new TLine(bins[0], 1.0, bins[nbins], 1.0);
    l->SetLineColor(kGray);
    l->Draw();
  }

  //draw histograms from collection
  for (Int_t i = 0; i < collectionSize; i++) {
    histograms_collection[i]->Draw("same");
    histograms_collection[i]->SetMarkerStyle(markers[0]);
    histograms_collection[i]->SetMarkerColor(colors[i]);
    histograms_collection[i]->SetLineColor(colors[i]);

    leg->AddEntry(histograms_collection[i], *pdfLegend_collection[i],   "LP");
  }

  leg->SetTextSize(gStyle->GetTextSize()*0.6);
  if (collectionSize >= 2){
    leg->Draw("same");
  }


  TLatex * textColl = new TLatex (0.18,0.82,"pp #sqrt{#it{s}} = 900 GeV, pilot beam 2021 MC");
  textColl->SetTextSize(0.04);
  textColl->SetNDC(kTRUE); //remove if I want x,y in TLatex to be in the coordinate system of the histogram
  textColl->Draw();
  TLatex * text_part = new TLatex (0.18,0.75,NamePart_Latex[ipart]);
  text_part->SetTextSize(0.04);
  text_part->SetNDC(kTRUE); //remove if I want x,y in TLatex to be in the coordinate system of the histogram
  text_part->Draw();

  canvas->SaveAs(*pdfName+".pdf");
}


void Draw_TH1_Histograms_withOneHistSystematics(TH1** histograms_collection, TString** pdfLegend_collection, Int_t collectionSize, TString* pdfName, TString* &texXtitle, TString* &texYtitle, Double_t* bins, Int_t nbins, Int_t ratioOption) {
  //draws only one hist with systematics, has to be the second one in the collection, and the associated systematics is taken from first hist in collection
  Float_t MaximumY = 0;
  for (Int_t i = 1; i < collectionSize; i++) {
    Float_t maxTemporary = histograms_collection[i]->GetMaximum();
    if (maxTemporary > MaximumY) {
      MaximumY = maxTemporary;
    }
  }

  ////plots
  TCanvas *canvas = new TCanvas ("canvas"+*pdfName, NamehistoInvMass[ipart], 1200, 800);
  canvas->cd(0); //ATTENTION: +1 because canvas divide starts counting from 1
  TH1 *hFrame = canvas->DrawFrame(bins[0],0,bins[nbins],1.4*MaximumY);
  hFrame->SetXTitle(texXtitle->Data());
  hFrame->SetYTitle(texYtitle->Data());
  TLegend * leg = new TLegend(0.59, 0.74, 0.87, 0.87);

  //draw line at ratio=1
  if (ratioOption == 1) {
    // TLine *l=new TLine(histograms_collection[0]->GetXaxis()->GetBinLowEdge(0), 1.0, histograms_collection[0]->GetXaxis()->GetBinUpEdge(nbins), 1.0);
    TLine *l=new TLine(bins[0], 1.0, bins[nbins], 1.0);
    l->SetLineColor(kGray);
    l->Draw();
  }

  //draw histograms from collection, ignoring first one that is the systematics for second one
  for (Int_t i = 1; i < collectionSize; i++) {
    histograms_collection[i]->Draw("same");
    histograms_collection[i]->SetMarkerStyle(markers[0]);
    histograms_collection[i]->SetMarkerColor(colors[i-1]);
    histograms_collection[i]->SetLineColor(colors[i-1]);

    leg->AddEntry(histograms_collection[i], *pdfLegend_collection[i],   "LP");
  }
  //draw systematics for 2d hist using 1st hist in collection
  histograms_collection[0]->Draw("E2,same");
  histograms_collection[0]->SetMarkerStyle(markers[0]);
  histograms_collection[0]->SetMarkerColor(colors[0]);
  histograms_collection[0]->SetLineColor(colors[0]);
  histograms_collection[0]->SetFillStyle(0); // To draw empty boxes

  leg->SetTextSize(gStyle->GetTextSize()*0.5);
  if (collectionSize >= 3) {
    leg->Draw("same");
  }

  TLatex * textColl = new TLatex (0.18,0.82,"pp #sqrt{#it{s}} = 900 GeV, pilot beam 2021 MC");
  textColl->SetTextSize(0.04);
  textColl->SetNDC(kTRUE); //remove if I want x,y in TLatex to be in the coordinate system of the histogram
  textColl->Draw();
  TLatex * text_part = new TLatex (0.18,0.75,NamePart_Latex[ipart]);
  text_part->SetTextSize(0.04);
  text_part->SetNDC(kTRUE); //remove if I want x,y in TLatex to be in the coordinate system of the histogram
  text_part->Draw();

  canvas->SaveAs(*pdfName+".pdf");
}

void Draw_TH1_Histograms_withOneHistSystematics_ratio(TH1** histograms_collection, TString** pdfLegend_collection, Int_t collectionSize, TString* pdfName, TString* &texXtitle, TString* &texYtitle, Double_t* bins, Int_t nbins, Int_t ratioOption) {
  //draws only one hist with systematics, has to be the second one in the collection, and the associated systematics is taken from first hist in collection

  ////plots
  TCanvas *canvas = new TCanvas ("canvas"+*pdfName, NamehistoInvMass[ipart], 1200, 800);
  canvas->cd(0); //ATTENTION: +1 because canvas divide starts counting from 1
  TH1 *hFrame = canvas->DrawFrame(bins[0],0.5,bins[nbins],1.5);
  hFrame->SetXTitle(texXtitle->Data());
  hFrame->SetYTitle(texYtitle->Data());
  TLegend * leg = new TLegend(0.59, 0.74, 0.87, 0.87);

  //draw line at ratio=1
  if (ratioOption == 1) {
    // TLine *l=new TLine(histograms_collection[0]->GetXaxis()->GetBinLowEdge(0), 1.0, histograms_collection[0]->GetXaxis()->GetBinUpEdge(nbins), 1.0);
    TLine *l=new TLine(bins[0], 1.0, bins[nbins], 1.0);
    l->SetLineColor(kGray);
    l->Draw();
  }

  //draw histograms from collection, ignoring first one that is the systematics for second one
  for (Int_t i = 1; i < collectionSize; i++) {
    histograms_collection[i]->Draw("same");
    histograms_collection[i]->SetMarkerStyle(markers[0]);
    histograms_collection[i]->SetMarkerColor(colors[i-1]);
    histograms_collection[i]->SetLineColor(colors[i-1]);

    leg->AddEntry(histograms_collection[i], *pdfLegend_collection[i],   "LP");
  }
  //draw systematics for 2d hist using 1st hist in collection
  histograms_collection[0]->Draw("E2,same");
  histograms_collection[0]->SetMarkerStyle(markers[0]);
  histograms_collection[0]->SetMarkerColor(colors[0]);
  histograms_collection[0]->SetLineColor(colors[0]);
  histograms_collection[0]->SetFillStyle(0); // To draw empty boxes

  leg->SetTextSize(gStyle->GetTextSize()*0.5);
  if (collectionSize >= 3) {
    leg->Draw("same");
  }

  TLatex * textColl = new TLatex (0.18,0.82,"pp #sqrt{#it{s}} = 900 GeV, pilot beam 2021 MC");
  textColl->SetTextSize(0.04);
  textColl->SetNDC(kTRUE); //remove if I want x,y in TLatex to be in the coordinate system of the histogram
  textColl->Draw();
  TLatex * text_part = new TLatex (0.18,0.75,NamePart_Latex[ipart]);
  text_part->SetTextSize(0.04);
  text_part->SetNDC(kTRUE); //remove if I want x,y in TLatex to be in the coordinate system of the histogram
  text_part->Draw();

  canvas->SaveAs(*pdfName+".pdf");
}

void Test_LambdaTrueMCYieldvsHEP(TH1D* &hstat, TH1D* &hsyst, TH1D* &hsystCorr, Int_t ipart, TFile* file_O2Analysis, Int_t nCutsIterations, TString* SaveAs_Title, TString* &texXtitle, TString* &texYtitle, Double_t* pTbins, Int_t nbinpT) {

  //HEP references for comparison
  TFile* file_HEPAnalysis = new TFile("HEPData-ins881474-v1-root.root");
  TH1D* hPtDiff_RawSpectrum_HEP[numPart]; //= (TH1D*)file_HEPAnalysis->Get("Table 1/Hist1D_y1");
  TH1D* hPtDiff_RawSpectrum_HEP_statErrors[numPart]; // = (TH1D*)file_HEPAnalysis->Get("Table 1/Hist1D_y1_e1");
  hPtDiff_RawSpectrum_HEP[0] = (TH1D*)file_HEPAnalysis->Get("Table 1/Hist1D_y1"); //K0S
  hPtDiff_RawSpectrum_HEP_statErrors[0] = (TH1D*)file_HEPAnalysis->Get("Table 1/Hist1D_y1_e1");
  hPtDiff_RawSpectrum_HEP[1] = (TH1D*)file_HEPAnalysis->Get("Table 2/Hist1D_y1"); //Lambda
  hPtDiff_RawSpectrum_HEP_statErrors[1] = (TH1D*)file_HEPAnalysis->Get("Table 2/Hist1D_y1_e1");
  hPtDiff_RawSpectrum_HEP[2] = (TH1D*)file_HEPAnalysis->Get("Table 3/Hist1D_y1"); //AntiLambda
  hPtDiff_RawSpectrum_HEP_statErrors[2] = (TH1D*)file_HEPAnalysis->Get("Table 3/Hist1D_y1_e1");

  if (ipart>2) {
    cout << "ipart>2, HEP data downloaded and coded only goes to 2" << endl;
  } 
  for(int ibinPt = 1; ibinPt <= nbinpT; ibinPt++){
    hPtDiff_RawSpectrum_HEP[ipart]->SetBinError(ibinPt, hPtDiff_RawSpectrum_HEP_statErrors[ipart]->GetBinContent(ibinPt));
  }

  TH1D* TrueV0PtSpectrum = (TH1D*)file_O2Analysis->Get("lambdakzero-particle-count-mc/h"+NamePart[ipart]+"Count_PtDiff");
  // TH1I* H1I_TotalEventCount = (TH1I*)file_O2Analysis->Get("lambdakzero-builder/hEventCounter");
  TH1I* H1I_TotalEventCount = (TH1I*)file_O2Analysis->Get("lambdakzero-particle-count-mc/hSelAndRecoMcCollCounter");

  TrueV0PtSpectrum->Scale(1./H1I_TotalEventCount->GetEntries());
  TH1D* TrueV0PtSpectrum_rebinned = (TH1D*)TrueV0PtSpectrum->Rebin(nbinpT,"TrueV0PtSpectrum_rebinned",pTbins);
  for (Int_t ibinPt = 0; ibinPt < nbinpT; ibinPt++)
  {
    Double_t dpT = TrueV0PtSpectrum_rebinned->GetXaxis()->GetBinWidth(ibinPt);
    Double_t drapidity = 1.5; //1.5
    Double_t d2N_dpTdy = TrueV0PtSpectrum_rebinned->GetBinContent(ibinPt) *1./dpT*1./drapidity;
    TrueV0PtSpectrum_rebinned->SetBinContent(ibinPt, d2N_dpTdy);
  }
  
  TString* pdfName = new TString("");;
  *pdfName += "CorrectewdYields_O2_vs_HEPdata";

  TH1* histograms_collection[2] = {TrueV0PtSpectrum_rebinned, hPtDiff_RawSpectrum_HEP[ipart]};
  TString* pdfLegend_collection[2] = {new TString("TrueV0PtSpectrum_rebinned"),new TString("hPtDiff_RawSpectrum_HEP")};
  Int_t ratioOption = 0;
  Draw_TH1_Histograms(histograms_collection, pdfLegend_collection, 2, pdfName, texXtitle, texYtitle, pTbins, nbinpT, ratioOption);

  /////////////////////////////////Does Nothing - only there to respect template/////////////////////////////////////////////
  hstat = (TH1D*)TrueV0PtSpectrum->Clone("hstat");
  hstat->Reset("M");
  // hstat->Divide(hPtDiff_RawSpectrum_O2,hPtDiff_RawSpectrum_HEP[ipart]);

  //////Error Bars///////
  hsyst     = (TH1D*)hstat->Clone("hsyst");
  hsyst->Reset("M");
  hsystCorr = (TH1D*)hstat->Clone("hsystCorr");
  hsystCorr->Reset("M");
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

void Test_V0DaugtherPairsTrackingEfficiency(TH1D* &hstat, TH1D* &hsyst, TH1D* &hsystCorr, Int_t ipart, TFile* file_O2Analysis, TString* SaveAs_Title, TString* &texXtitle, TString* &texYtitle, Double_t* pTbins, Int_t nbinpT) {


  TH1D* ReconstructedDaughterPairs = (TH1D*)file_O2Analysis->Get("aimeric-count-reco-pions-from-k0-s/hPtDiffK0ShortFullyReco");
  TH1D* TrueV0PtSpectrum = (TH1D*)file_O2Analysis->Get("aimeric-count-reco-pions-from-k0-s/hPtDiffK0ShortGen");
  ReconstructedDaughterPairs->Sumw2();
  TrueV0PtSpectrum->Sumw2();
  /////////////////////////////////Does Nothing - only there to respect template/////////////////////////////////////////////
  hstat = (TH1D*)ReconstructedDaughterPairs->Clone("hstat");
  hstat->Sumw2();
  hstat->Reset("M");
  hstat->Divide(ReconstructedDaughterPairs,TrueV0PtSpectrum);

  //////Error Bars///////
  hsyst     = (TH1D*)hstat->Clone("hsyst");
  hsyst->Reset("M");
  hsystCorr = (TH1D*)hstat->Clone("hsystCorr");
  hsystCorr->Reset("M");
  hstat     ->Sumw2(); 
  hsyst     ->Sumw2();
  hsystCorr ->Sumw2(); 
  Int_t nbinx = hstat->GetNbinsX();
  
  for(Int_t ibinx = 1; ibinx <= nbinx; ibinx++){
    hsyst->SetBinContent (ibinx, hstat->GetBinContent(ibinx));
    hsyst->SetBinError   (ibinx, hstat->GetBinContent(ibinx)*0.000001);
    hsystCorr->SetBinContent (ibinx, hstat->GetBinContent(ibinx));
    hsystCorr->SetBinError   (ibinx, hstat->GetBinContent(ibinx)*0.15);
  }
  *SaveAs_Title += "DaughterPairsRecoEfficiency";
  texXtitle = texPtX;
  texYtitle = texDaughterPairsRecoEfficiency;
}


void Test_V0DaugtherPairsTrackingEfficiency_PostDcaFitter(TH1D* &hstat, TH1D* &hsyst, TH1D* &hsystCorr, Int_t ipart, TFile* file_O2Analysis, TString* SaveAs_Title, TString* &texXtitle, TString* &texYtitle, Double_t* pTbins, Int_t nbinpT) {


  TH1D* ReconstructedDaughterPairs = (TH1D*)file_O2Analysis->Get("k0s-reco-study/hPtDiffK0ShortFullyReco_PassCuts_debug9");
  TH1D* TrueV0PtSpectrum = (TH1D*)file_O2Analysis->Get("k0s-reco-study/hPtDiffK0ShortGen");
  ReconstructedDaughterPairs->Sumw2();
  TrueV0PtSpectrum->Sumw2();
  /////////////////////////////////Does Nothing - only there to respect template/////////////////////////////////////////////
  hstat = (TH1D*)ReconstructedDaughterPairs->Clone("hstat");
  hstat->Sumw2();
  hstat->Reset("M");
  hstat->Divide(ReconstructedDaughterPairs,TrueV0PtSpectrum);

  //////Error Bars///////
  hsyst     = (TH1D*)hstat->Clone("hsyst");
  hsyst->Reset("M");
  hsystCorr = (TH1D*)hstat->Clone("hsystCorr");
  hsystCorr->Reset("M");
  hstat     ->Sumw2(); 
  hsyst     ->Sumw2();
  hsystCorr ->Sumw2(); 
  Int_t nbinx = hstat->GetNbinsX();
  
  for(Int_t ibinx = 1; ibinx <= nbinx; ibinx++){
    hsyst->SetBinContent (ibinx, hstat->GetBinContent(ibinx));
    hsyst->SetBinError   (ibinx, hstat->GetBinContent(ibinx)*0.000001);
    hsystCorr->SetBinContent (ibinx, hstat->GetBinContent(ibinx));
    hsystCorr->SetBinError   (ibinx, hstat->GetBinContent(ibinx)*0.15);
  }
  *SaveAs_Title += "DaughterPairsRecoEfficiency_PostDcaFitter";
  texXtitle = texPtX;
  texYtitle = texDaughterPairsRecoEfficiency;
}