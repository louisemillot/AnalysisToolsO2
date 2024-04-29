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

#include<array>
using namespace std;

void SetStyle(Bool_t graypalette=kFALSE);
void myPadSetUp(TPad *currentPad, float currentLeft=0.11, float currentTop=0.04, float currentRight=0.04, float currentBottom=0.15);
void DrawLogo (int logo=0, double xmin =  0.28, double ymin= 0.68) ;
void FakeHistosOnlyForExample(TH1*&hstat, TH1*&hsyst, TH1*&hsystCorr);
void LoadLibs();

void test_SetDefaultSumw2();

void PseudoEfficiency_HEPcomparison_allCountSignalRegion(TH1D* &hstat, TH1D* &hsyst, TH1D*&hsystCorr, int ipart, TFile* file_O2Analysis, TString* SaveAs_Title, TString* &texXtitle, TString* &texYtitle, double* pTbins, int nbinpT);
void PseudoEfficiency_HEPcomparison_allCountSignalRegionOLD_withBugs(TH1D* &hstat, TH1D* &hsyst, TH1D*&hsystCorr, int ipart, TFile* file_O2Analysis, TString* SaveAs_Title, double* pTbins, int nbinpT);
void RawSpectrum_O2data_allCountSignalRegion(TH1D* &hstat, TH1D* &hsyst, TH1D*&hsystCorr, int ipart, TFile* file_O2Analysis, TString* SaveAs_Title, TString* &texXtitle, TString* &texYtitle, double* pTbins, int nbinpT);
void InvMass_Plot(TH1D* &hstat, TH1D* &hsyst, TH1D*&hsystCorr, int ipart, TFile* file_O2Analysis, TString* SaveAs_Title, TString* &texXtitle, TString* &texYtitle, bool isMC);
void PtDistribution_O2data(TH1D* &hstat, TH1D* &hsyst, TH1D*&hsystCorr, int ipart, TFile* file_O2Analysis, TString* SaveAs_Title, TString* &texXtitle, TString* &texYtitle);
void DcaHistoProcessing_O2data(TH1D* &hstat, TH1D* &hsyst, TH1D*&hsystCorr, int ipart, TFile* file_O2Analysis);
void Efficiency_O2data_allCountSignalRegion(TH1D* &hstat, TH1D* &hsyst, TH1D*&hsystCorr, int ipart, TFile* file_O2Analysis, TString* SaveAs_Title, TString* &texXtitle, TString* &texYtitle, double* pTbins, int nbinpT);
void InvMassDistributions_withFit(TH1D* &hstat, TH1D* &hsyst, TH1D* &hsystCorr, int ipart, TFile* file_O2Analysis, TString* SaveAs_Title, TString* &texXtitle, TString* &texYtitle, double* pTbins, int nbinpT);
void InvMassDistributions_withFit_MC_datalike(TH1D* &hstat, TH1D* &hsyst, TH1D* &hsystCorr, int ipart, TFile* file_O2Analysis, TString* SaveAs_Title, TString* &texXtitle, TString* &texYtitle, double* pTbins, int nbinpT);


void PtDifferential_SigmaOfFit(TH1D* &hstat, TH1D* &hsyst, TH1D*&hsystCorr, int ipart, TFile* file_O2Analysis, TString* SaveAs_Title, TString* &texXtitle, TString* &texYtitle, double* pTbins, int nbinpT, bool isMC);
void RawSpectrum_O2data_withFit(TH1D* &hstat, TH1D* &hsyst, TH1D*&hsystCorr, int ipart, TFile* file_O2Analysis, TString* SaveAs_Title, TString* &texXtitle, TString* &texYtitle, double* pTbins, int nbinpT);
void Efficiency_O2MCdata_withFit(TH1D* &hstat, TH1D* &hsyst, TH1D*&hsystCorr, int ipart, TFile* file_O2Analysis, TString* SaveAs_Title, TString* &texXtitle, TString* &texYtitle, double* pTbins, int nbinpT);
void Efficiency_O2data_truePt_trueV0s(TH1D* &hstat, TH1D* &hsyst, TH1D*&hsystCorr, int ipart, TFile* file_O2Analysis, TString* SaveAs_Title, TString* &texXtitle, TString* &texYtitle, double* pTbins, int nbinpT);
void RawSpectrum_O2data_truePt_trueV0s(TH1D* &hstat, TH1D* &hsyst, TH1D*&hsystCorr, int ipart, TFile* file_O2Analysis, TString* SaveAs_Title, TString* &texXtitle, TString* &texYtitle, double* pTbins, int nbinpT);

void V0Daughters_TrackingEfficiencies(TH1D* &hstat, TH1D* &hsyst, TH1D*&hsystCorr, int ipart, TFile* file_O2Analysis, TString* SaveAs_Title, double* pTbins, int nbinpT);
void V0Daughters_TrackingEfficiencies_HistogramProcessing(TH1D* &daughter1, TH1D* &daughter2, int ipart, TFile* file_O2Analysis, double* pTbinsProtons, int nbinpTProtons, double* pTbinsPions, int nbinpTPions);
void V0Daughters_TrackingEfficiencies_Daughter1(TH1D* &hstat, TH1D* &hsyst, TH1D*&hsystCorr, int ipart, TFile* file_O2Analysis, TString* SaveAs_Title, TString* &texXtitle, TString* &texYtitle, TString* &Extra, double* pTbinsProtons, int nbinpTProtons, double* pTbinsPions, int nbinpTPions);
void V0Daughters_TrackingEfficiencies_Daughter2(TH1D* &hstat, TH1D* &hsyst, TH1D*&hsystCorr, int ipart, TFile* file_O2Analysis, TString* SaveAs_Title, TString* &texXtitle, TString* &texYtitle, TString* &Extra, double* pTbinsProtons, int nbinpTProtons, double* pTbinsPions, int nbinpTPions);

void pT_Spectrum_postAnalyserCuts(TH1D* &hstat, TH1D* &hsyst, TH1D*&hsystCorr, int ipart, TFile* file_O2Analysis, TString* SaveAs_Title, TString* &texXtitle, TString* &texYtitle, double* pTbins, int nbinpT, bool isMC);
void pT_Spectrum_preAnalyserCutsK0S(TH1D* &hstat, TH1D* &hsyst, TH1D*&hsystCorr, int ipart, TFile* file_O2Analysis, TString* SaveAs_Title, TString* &texXtitle, TString* &texYtitle, double* pTbins, int nbinpT, bool isMC);

void Efficiency_DcaScan(TH1D* &hstat, TH1D* &hsyst, TH1D*&hsystCorr, int ipart, TFile* file_O2Analysis, TString* SaveAs_Title, TString* &texXtitle, TString* &texYtitle, double* pTbins, int nbinpT, double maxDCAcut);
void Efficiency_TpcCrossedRowsScan(TH1D* &hstat, TH1D* &hsyst, TH1D*&hsystCorr, int ipart, TFile* file_O2Analysis, TString* SaveAs_Title, TString* &texXtitle, TString* &texYtitle, double* pTbins, int nbinpT, double minTpcCrossedRowscut, double maxTpcCrossedRowscut);

void Systematics_CutVariations_Graphs(TH1D* &hstat, TH1D* &hsyst, TH1D* &hsystCorr, int ipart, TString* SaveAs_Title, TString* &texXtitle, TString* &texYtitle, double* pTbins, int nbinpT);
void CorrectedYield_withSystematics(TH1D* &hstat, TH1D* &hsyst, TH1D* &hsystCorr, int ipart, TString* SaveAs_Title, TString* &texXtitle, TString* &texYtitle, double* pTbins, int nbinpT);
void Draw_XsiPlus_MassPlot(TH1D* &hstat, TH1D* &hsyst, TH1D* &hsystCorr, int ipart, TFile *file_O2Analysis, TString* SaveAs_Title, TString* &texXtitle, TString* &texYtitle, double* pTbins, int nbinpT);
void Draw_FeedDownMatrix(TH1D* &hstat, TH1D* &hsyst, TH1D* &hsystCorr, int ipart, TFile *file_O2Analysis, TString* SaveAs_Title, TString* &texXtitle, TString* &texYtitle, double* pTbins, int nbinpT);
void Draw_FeedDown_ipart(TH1D* &hstat, TH1D* &hsyst, TH1D* &hsystCorr, int ipart, TFile *file_O2Analysis, TString* SaveAs_Title, TString* &texXtitle, TString* &texYtitle, double* pTbins, int nbinpT);
void Draw_FeedDown_LambdaAntiLambda(TH1D* &hstat, TH1D* &hsyst, TH1D* &hsystCorr, int ipart, TFile *file_O2Analysis, TString* SaveAs_Title, TString* &texXtitle, TString* &texYtitle, double* pTbins, int nbinpT);
void Draw_TH1_Histograms(TH1** histograms_collection, TString** pdfLegend_collection, int collectionSize, TString* pdfName, TString* &texXtitle, TString* &texYtitle, double* bins, int nbins, int ratioOption);
void Draw_TH1_Histograms_Efficiency(TH1** histograms_collection, TString** pdfLegend_collection, int collectionSize, TString* pdfName, TString* &texXtitle, TString* &texYtitle, double* bins, int nbins, int ratioOption);
void Draw_TH1_Histograms_withOneHistSystematics(TH1** histograms_collection, TString** pdfLegend_collection, int collectionSize, TString* pdfName, TString* &texXtitle, TString* &texYtitle, double* bins, int nbins, int ratioOption);
void Draw_TH1_Histograms_withOneHistSystematics_ratio(TH1** histograms_collection, TString** pdfLegend_collection, int collectionSize, TString* pdfName, TString* &texXtitle, TString* &texYtitle, double* bins, int nbins, int ratioOption);
void Test_LambdaTrueMCYieldvsHEP(TH1D* &hstat, TH1D* &hsyst, TH1D* &hsystCorr, int ipart, TFile* file_O2Analysis, TString* SaveAs_Title, TString* &texXtitle, TString* &texYtitle, double* pTbins, int nbinpT);
void Test_V0DaugtherPairsTrackingEfficiency(TH1D* &hstat, TH1D* &hsyst, TH1D* &hsystCorr, int ipart, TFile* file_O2Analysis, TString* SaveAs_Title, TString* &texXtitle, TString* &texYtitle, double* pTbins, int nbinpT);
void Test_V0DaugtherPairsTrackingEfficiency_PostDcaFitter(TH1D* &hstat, TH1D* &hsyst, TH1D* &hsystCorr, int ipart, TFile* file_O2Analysis, TString* SaveAs_Title, TString* &texXtitle, TString* &texYtitle, double* pTbins, int nbinpT);
void Draw_Fit_Mass_Histograms(TH1D* &hstat, TH1D* &hsyst, TH1D* &hsystCorr, TFile *file_O2Analysis, TString* SaveAs_Title, TString* &texXtitle, TString* &texYtitle, int ipart, double* pTbins, int nbinpT);
void Draw_Mass_Histograms(TH1D* &hstat, TH1D* &hsyst, TH1D* &hsystCorr, TFile *file_O2Analysis, TString* SaveAs_Title, TString* &texXtitle, TString* &texYtitle, int ipart, double* pTbins, int nbinpT);
void Draw_Rapidity_MC(TH1D* &hstat, TH1D* &hsyst, TH1D* &hsystCorr, TFile *file_O2Analysis, TString* SaveAs_Title, TString* &texXtitle, TString* &texYtitle, int ipart);
void PseudoEfficiency_O2vsHEP(TH1D* &hstat, TH1D* &hsyst, TH1D* &hsystCorr, int ipart, TString* SaveAs_Title, TString* &texXtitle, TString* &texYtitle, double* pTbins, int nbinpT);
void Draw_pT_vs_pTmc(TH1D* &hstat, TH1D* &hsyst, TH1D* &hsystCorr, TFile *file_O2Analysis, TString* SaveAs_Title, TString* &texXtitle, TString* &texYtitle, int ipart);
void Draw_FinderEfficiencies(TH1D* &hstat, TH1D* &hsyst, TH1D* &hsystCorr, TFile *file_O2Analysis_Run1Cuts, TFile *file_O2Analysis_SVCuts, TString* SaveAs_Title, TString* &texXtitle, TString* &texYtitle, int ipart, double* pTbins, int nbinpT);
void Draw_XiChargedToNeutralRatio(TH1D* &hstat, TH1D* &hsyst, TH1D* &hsystCorr, TFile *file_O2Analysis, TString* SaveAs_Title, TString* &texXtitle, TString* &texYtitle, int ipart, double* pTbins, int nbinpT);

//PhD_Amendments
void Draw_Mass_Histograms_multipleRuns(TH1D* &hstat, TH1D* &hsyst, TH1D* &hsystCorr, TString* SaveAs_Title, TString* &texXtitle, TString* &texYtitle, int ipart);
void Draw_TH1_Histograms_Runs(TH1D** histograms_collection, int collectionSize, TString* pdfName, TString* &texXtitle, TString* &texYtitle);
void Draw_EfficiencyWithSigmoidFit(TH1D* &hstat, TH1D* &hsyst, TH1D* &hsystCorr, TFile *file_O2Analysis, TString* SaveAs_Title, TString* &texXtitle, TString* &texYtitle, int ipart, double* pTbins, int nbinpT);
void Get_McEfficiency_vsPt_alwaysFitSmooth(TH1D* &hMcEfficiency_vsPt, TH3D* H3D_detectedV0s_TruePt, TH1D* TrueV0PtSpectrum, int ipart, double* pTbins, int nbinpT, int ibinXaxisCut_low, int ibinXaxisCut_high, double SideBandSizeMultiplierModifier, float SignalMean, float SignalStandardDev, int InvMassRebinFactor);
void Draw_PtResponseMatrix(TH1D* &hstat, TH1D* &hsyst, TH1D* &hsystCorr, TFile *file_O2Analysis, TString* SaveAs_Title, TString* &texXtitle, TString* &texYtitle, double* pTbins, int nbinpT);
void Get_PtResponseMatrix(TH2D * H2D_PtResponseMatrix, TFile *file_O2Analysis, double* pTbins, int nbinpT);
void Get_PtResponseMatrix_Density(TH2D * H2D_PtResponseMatrix, TFile *file_O2Analysis, double* pTbins, int nbinpT);

bool reject = false ;
float fparab(float *x, float *par);
float fline(float *x, float *par);
double myLevyPtXsi(double *pt, double *par);
float fsigmoid(float *x, float *par);
double fgenerallogistic(double *x, double *par);

// Histogram processing functions; 

// returns pt diff raw spectrum in hPtDiff_RawSpectrum_O2
// void PtDifferential_Fit(TH1D* &hPtDiff_RawSpectrum_O2, TH3D* H3D_DetectedV0s_O2, TH3D* H3D_DetectedV0s_O2_rebinnedZ, int ipart, double* pTbins, int nbinpT, double* parGaussianParab_Sigma);
void PtDifferential_Fit(TH1D* &hPtDiff_RawSpectrum_O2, TH3D* H3D_DetectedV0s_O2, TH3D* H3D_DetectedV0s_O2_rebinnedZ, int ipart, double* pTbins, int nbinpT, double* parGaussianParab_Sigma, double* parGaussianParab_Sigma_error);

//Custom Get functions
void Get_McEfficiency_vsPt(TH1D* &hMcEfficiency_vsPt, TH3D* H3D_detectedV0s_TruePt, TH1D* TrueV0PtSpectrum, int ipart, double* pTbins, int nbinpT, int ibinXaxisCut_low, int ibinXaxisCut_high, double SideBandSizeMultiplierModifier, float SignalMean, float SignalStandardDev, int InvMassRebinFactor);
void Get_McSignalBackgroundRatio_vsPt(TH1D* &hSignalBackgroundRatio_vsPt, TH3D* H3D_detectedV0s_DataLike, int ipart, double* pTbins, int nbinpT, int ibinXaxisCut_low, int ibinXaxisCut_high);
void Get_Systematics_OneCut(TH1D** hSystematicUncertainty, TH1D** hSystematicUncertainty_PreBarlow, int ipart, TFile **file_O2Analysis_CutVariation_array, int nCutsIterations, double* pTbins, int nbinpT, int InvMassRebinFactor);
void Get_Systematics_OneCut_McTrueAndMcDatalikeDifferentCuts(TH1D* hSystematicUncertainty, TH1D* hSystematicUncertainty_PreBarlow, int ipart, TFile **file_O2Analysis_CutVariation_Datalike_array, TFile **file_O2Analysis_CutVariation_True_array, int nCutsIterations, double* pTbins, int nbinpT, int InvMassRebinFactor);
void Get_Systematics_Fit_SignalExtractionType(TH1D* hSystematicUncertainty, TH1D* hSystematicUncertainty_PreBarlow, int ipart, TH3D* H3D_detectedV0s_DataLike, TH1I* H1I_SelectedEventCount, TH3D* H3D_detectedV0s_TruePt, TH1D* TrueV0PtSpectrum, TFile* file_O2Analysis, double* pTbins, int nbinpT, int InvMassRebinFactor);
void Get_Systematics_Fit_SideBandVariation(TH1D* hSystematicUncertainty, TH1D* hSystematicUncertainty_PreBarlow, int ipart, TH3D* H3D_detectedV0s_DataLike, TH1I* H1I_SelectedEventCount, TH3D* H3D_detectedV0s_TruePt, TH1D* TrueV0PtSpectrum, TFile* file_O2Analysis, double* pTbins, int nbinpT, int InvMassRebinFactor);
void Get_Systematics_Fit_Binning(TH1D* hSystematicUncertainty, TH1D* hSystematicUncertainty_PreBarlow, int ipart, TH3D* H3D_detectedV0s_DataLike, TH1I* H1I_SelectedEventCount, TH3D* H3D_detectedV0s_TruePt, TH1D* TrueV0PtSpectrum, TFile* file_O2Analysis, double* pTbins, int nbinpT);
void Get_RawYield_vsPt(TH1D* &hRawYield_vsPt, TH3D* H3D_detectedV0s_DataLike, TH1I* H1I_SelectedEventCount, TFile *file_O2Analysis, double* SignalMean, double* SignalStandardDev, int ipart, double* pTbins, int nbinpT, int ibinXaxisCut_low,  int ibinXaxisCut_high, int warning_cutArry_ID, double SideBandSizeMultiplierModifier, int SignalExtractionType, int InvMassRebinFactor);
void Get_RawYield_vsPt_FeeddownCorrected(TH1D* &hRawYield_vsPt, TH3D* H3D_detectedV0s_DataLike, TH1I* H1I_SelectedEventCount, TFile *file_O2Analysis, double* SignalMean, double* SignalStandardDev, int ipart, double* pTbins, int nbinpT, int ibinXaxisCut_low,  int ibinXaxisCut_high, int warning_cutArry_ID, double SideBandSizeMultiplierModifier, int SignalExtractionType, int InvMassRebinFactor);
void Get_FeedDown(TH1D * H1D_FeedDownCorrection, int ipart_option, TFile *file_O2Analysis, double* pTbins, int nbinpT);
void Get_FeedDown_RawCount(TH1D * H1D_FeedDownCorrection, int ipart_option, TFile *file_O2Analysis, double* pTbins, int nbinpT);
void Get_RawYield_vsPt_FeeddownCorrected_RawCount(TH1D* &hRawYield_vsPt, TH3D* H3D_detectedV0s_DataLike, TH1I* H1I_SelectedEventCount, TFile *file_O2Analysis, double* SignalMean, double* SignalStandardDev, int ipart, double* pTbins, int nbinpT, int ibinXaxisCut_low,  int ibinXaxisCut_high, int warning_cutArry_ID, double SideBandSizeMultiplierModifier, int SignalExtractionType, int InvMassRebinFactor);
void Get_RawYield_vsPt_FeeddownCorrected_withUnfolding(TH1D* &hRawYield_vsPt, TH3D* H3D_detectedV0s_DataLike, TH1I* H1I_SelectedEventCount, TFile *file_O2Analysis, double* SignalMean, double* SignalStandardDev, int ipart, double* pTbins, int nbinpT, int ibinXaxisCut_low,  int ibinXaxisCut_high, int warning_cutArry_ID, double SideBandSizeMultiplierModifier, int SignalExtractionType, int InvMassRebinFactor);


std::array<double, 3> GetTotalAndBackgroundCount_SidebandBackgroundExtrapolation(TH1D* H1D_detectedV0s_DataLike_SmallpTinterval, double SignalMean, double SignalStandardDev, double SideBandSizeMultiplierModifier);
std::array<double, 3> GetTotalAndBackgroundCount_BackgroundIntegration(TH1D* H1D_detectedV0s_DataLike_SmallpTinterval, double SignalMean, double SignalStandardDev, double SideBandSizeMultiplierModifier, TFitResultPtr FitResult, TF1 *bkgparab);

// Preferred colors and markers
const int fillColors[] = {kGray+1,  kRed-10, kBlue-9, kGreen-8, kMagenta-9, kOrange-9,kCyan-8,kYellow-7}; // for syst bands
const int colors[]     = {kBlack, kRed+1 , kBlue+1, kGreen+3, kMagenta+1, kOrange-1,kCyan+2,kYellow+2};
const int markers[]    = {kFullSquare,kFullCircle,kOpenCircle,kOpenSquare,kOpenDiamond,kOpenCross,kFullCross,kFullDiamond,kFullStar,kOpenStar};

//Options:
// TFile* file_O2Analysis = new TFile("AnalysisResult_TpcXrows_Ref.root");
// TFile* file_O2Analysis = new TFile("AnalysisResult_TpcXrows_Tight.root");
// TFile* file_O2Analysis = new TFile("AnalysisResult_TpcXrows_veryTight.root");
// TFile* file_O2Analysis = new TFile("AnalysisResult_DCA000.root");

// TFile* file_O2Analysis = new TFile("AnalysisResults_noMC.root");
// TFile* file_O2Analysis = new TFile("AnalysisResultsMC.root");
const int numPart = 3;
const int ipart = 1;

// //****************************  Apass4  ***********************************//
// TFile* file_O2Analysis = new TFile("AnalysisResult_Ref.root");
// const bool UseFinder = false;
// const int N_variableCut = 1;
// const TString VariableCut[N_variableCut] = {"TpcXrows"};
// const  int nCutsIterations_max = 2;
// const  int nCutsIterations_array[N_variableCut] = {2};
// bool doFitEfficiencySmooth = false;

// TFile* file_O2Analysis_collection_TpcXrows[nCutsIterations_array[0]] = {file_O2Analysis, file_O2Analysis};
// TFile** file_O2Analysis_CutVariation_Collection_array[N_variableCut] = {file_O2Analysis_collection_TpcXrows};

// //****************************  Finder  ***********************************//
// TFile* file_O2Analysis = new TFile("AnalysisResults_Finder_Run1Cuts.root");
// TFile* file_O2Analysis1 = new TFile("AnalysisResults_Finder_Run1Cuts.root");
// TFile* file_O2Analysis2 = new TFile("AnalysisResults_Finder_SVcuts.root");
// const bool UseFinder = true;
// const int N_variableCut = 1;
// const TString VariableCut[N_variableCut] = {"TpcXrows"};
// const  int nCutsIterations_max = 2;
// const  int nCutsIterations_array[N_variableCut] = {2};
// // TFile* file_O2Analysis01 = new TFile("AnalysisResult_TpcXrows_Ref.root");
// TFile* file_O2Analysis02 = new TFile("AnalysisResult_TpcXrows_Tight.root");
// TFile* file_O2Analysis03 = new TFile("AnalysisResult_TpcXrows_veryTight.root");
// TFile* file_O2Analysis_collection_TpcXrows[nCutsIterations_array[0]] = {file_O2Analysis02, file_O2Analysis03};
// TFile** file_O2Analysis_CutVariation_Collection_array[N_variableCut] = {file_O2Analysis_collection_TpcXrows};

// //****************************  TpcXrows  ***********************************//
// TFile* file_O2Analysis = new TFile("AnalysisResult_TpcXrows_Ref.root");
// const bool UseFinder = false;
// const int N_variableCut = 1;
// const TString VariableCut[N_variableCut] = {"TpcXrows"};
// const  int nCutsIterations_max = 2;
// const  int nCutsIterations_array[N_variableCut] = {2};
// // TFile* file_O2Analysis01 = new TFile("AnalysisResult_TpcXrows_Ref.root");
// TFile* file_O2Analysis02 = new TFile("AnalysisResult_TpcXrows_Tight.root");
// TFile* file_O2Analysis03 = new TFile("AnalysisResult_TpcXrows_veryTight.root");
// TFile* file_O2Analysis_collection_TpcXrows[nCutsIterations_array[0]] = {file_O2Analysis02, file_O2Analysis03};
// TFile** file_O2Analysis_CutVariation_Collection_array[N_variableCut] = {file_O2Analysis_collection_TpcXrows};

// ****************************  FullAnalysis  ***********************************//
TFile* file_O2Analysis = new TFile("AnalysisResult_Ref.root");
TFile* file_O2Analysis_unfoldingInfo = new TFile("AnalysisResult_Ref_unfoldingInfo.root");

// TFile* file_O2Analysis01 = new TFile("AnalysisResult_Ref.root");

TFile* file_O2Analysis02 = new TFile("AnalysisResult_TpcXrows_Tight.root");
TFile* file_O2Analysis03 = new TFile("AnalysisResult_TpcXrows_veryTight.root");

TFile* file_O2Analysis04 = new TFile("AnalysisResult_dcav0dau_Tight.root");
TFile* file_O2Analysis05 = new TFile("AnalysisResult_dcav0dau_veryTight.root");

TFile* file_O2Analysis06 = new TFile("AnalysisResult_dcatopv_Tight.root");
TFile* file_O2Analysis07 = new TFile("AnalysisResult_dcatopv_veryTight.root");


TFile* file_O2Analysis08 = new TFile("AnalysisResult_v0radius_Tight.root");
TFile* file_O2Analysis09 = new TFile("AnalysisResult_v0radius_veryTight.root");

TFile* file_O2Analysis10 = new TFile("AnalysisResult_cosPA_veryLoose.root");
TFile* file_O2Analysis11 = new TFile("AnalysisResult_cosPA_Loose.root");
TFile* file_O2Analysis12 = new TFile("AnalysisResult_cosPA_Tight.root");
TFile* file_O2Analysis13 = new TFile("AnalysisResult_cosPA_veryTight.root");

TFile* file_O2Analysis14 = new TFile("AnalysisResult_TpcPID_Tight.root");
TFile* file_O2Analysis15 = new TFile("AnalysisResult_TpcPID_veryTight.root");
TFile* file_O2Analysis16 = new TFile("AnalysisResult_TpcPID_Loose.root");
TFile* file_O2Analysis17 = new TFile("AnalysisResult_TpcPID_veryLoose.root");

bool doFitEfficiencySmooth = false;
bool doUnfolding = true; // true doesn't work well for now
const bool UseFinder = false;
const int N_variableCut = 6;
const TString VariableCut[N_variableCut] = {"TpcXrows", "dcav0dau", "dcatopv", "v0radius", "cosPA", "TpcPID"};
const int nCutsIterations_max = 4;
// const  int nCutsIterations_array[N_variableCut] = {2, 2, 2, 2, 4, 4};
const int nCutsIterations_TpcXrows = 2;
const int nCutsIterations_dcav0dau = 2;
const int nCutsIterations_dcatopv = 2;
const int nCutsIterations_v0radius = 2;
const int nCutsIterations_cosPA = 4;
const int nCutsIterations_TpcPID = 2;
const int nCutsIterations_array[N_variableCut] = {2, 2, 2, 2, 4, 2};
// TFile* file_O2Analysis_collection_TpcXrows[nCutsIterations_array[0]] = {file_O2Analysis02, file_O2Analysis03};
// TFile* file_O2Analysis_collection_dcav0dau[nCutsIterations_array[1]] = {file_O2Analysis04, file_O2Analysis05};
// TFile* file_O2Analysis_collection_dcatopv[nCutsIterations_array[2]] = {file_O2Analysis06, file_O2Analysis07};
// TFile* file_O2Analysis_collection_v0radius[nCutsIterations_array[3]] = {file_O2Analysis08, file_O2Analysis09};
// TFile* file_O2Analysis_collection_cosPA[nCutsIterations_array[4]] = {file_O2Analysis10, file_O2Analysis11, file_O2Analysis12, file_O2Analysis13};
// // TFile* file_O2Analysis_collection_TpcPID[nCutsIterations_array[5]] = {file_O2Analysis14, file_O2Analysis15, file_O2Analysis16, file_O2Analysis17};
// TFile* file_O2Analysis_collection_TpcPID[nCutsIterations_array[5]] = {file_O2Analysis15, file_O2Analysis16};
TFile* file_O2Analysis_collection_TpcXrows[nCutsIterations_TpcXrows] = {file_O2Analysis02, file_O2Analysis03};
TFile* file_O2Analysis_collection_dcav0dau[nCutsIterations_dcav0dau] = {file_O2Analysis04, file_O2Analysis05};
TFile* file_O2Analysis_collection_dcatopv[nCutsIterations_dcatopv] = {file_O2Analysis06, file_O2Analysis07};
TFile* file_O2Analysis_collection_v0radius[nCutsIterations_v0radius] = {file_O2Analysis08, file_O2Analysis09};
TFile* file_O2Analysis_collection_cosPA[nCutsIterations_cosPA] = {file_O2Analysis10, file_O2Analysis11, file_O2Analysis12, file_O2Analysis13};
// TFile* file_O2Analysis_collection_TpcPID[nCutsIterations_array[5]] = {file_O2Analysis14, file_O2Analysis15, file_O2Analysis16, file_O2Analysis17};
TFile* file_O2Analysis_collection_TpcPID[nCutsIterations_TpcPID] = {file_O2Analysis15, file_O2Analysis16};
TFile** file_O2Analysis_CutVariation_Collection_array[N_variableCut] = {file_O2Analysis_collection_TpcXrows, file_O2Analysis_collection_dcav0dau, file_O2Analysis_collection_dcatopv, file_O2Analysis_collection_v0radius, file_O2Analysis_collection_cosPA, file_O2Analysis_collection_TpcPID};

const float N_SigmaBarlow = 1.9; //often 2 in ALICE LF analyses
const int N_SideBandSizeVariation = 3;
const int id_referenceSideband = 0.;
const float DefaultSideBandSize_array[numPart] = {4., 4., 4.};
const float SideBandSizeMultiplierModifier_array[N_SideBandSizeVariation] = {0, -0.5, -1}; // can't go above 1: 3+1 is 4 and that's the higher we can go with lambda having 4MeV sigma in rare cases and 2*4*4MeV going below the lower cutoff of pi+p masses; can't go below 0 because then we start losing signal that's within the guassian area; +/-3Sigma is 99.9% of gaussian area 
const int SignalExtractionType_default = 1; //default is 1 for integration as opposed to 0 for bin counting; if using bin counting (0) then make sure H3D inv mass window is large enough to accomodate the 12sigma on both sides!
const int InvMassRebinFactor_standard[numPart] = {10,2,2}; // I like 10 best for K0s at InvMass 0.4 to 0.6 with 400 bins, do variations around that
const int N_BinningVariation = 4;
const int InvMassRebinFactor_variation_array[numPart][N_BinningVariation] = {{5, 10, 2, 20}, {1, 2, 4, 5}, {1, 2, 4, 5}};

const TString NamehistoInvMass[numPart] = {"InvMassK0S", "InvMassLambda", "InvMassAntiLambda"};//,"InvMassXiPlus", "InvMassXiMinus", "InvMassOmegaPlus", "InvMassOmegaMinus"};
const TString NamePart[numPart] = {"K0Short", "Lambda", "AntiLambda"};//, "XiPlu", "XiMin", "OmPlu", "OmMin"};
const TString NamePart_multiStrange[2] = {"XsiMinus", "XsiPlus"};
const TString NamePart_multiStrange_Latex[2] = {"#Xi^{-}", "#Xi^{+}"};
const TString NamePart_Latex[numPart] = {"K^{0}_{S}", "#Lambda", "#bar{#Lambda}"};//, "XiPlu", "XiMin", "OmPlu", "OmMin"};
// const float LowerLimitDisplay_Mass[numPart] = {0.455, 1.085, 1.085};
// const float UpperLimitDisplay_Mass[numPart] = {0.525, 1.145, 1.145};
const float LowerLimitDisplay_Mass[numPart] = {0.4, 1.082, 1.082};
const float UpperLimitDisplay_Mass[numPart] = {0.6, 1.146, 1.146};

const float MassPart[numPart] = {0.497611, 1.115683, 1.115683};

const float min_range_signal[numPart] = {0.47, 1.105, 1.105};
const float max_range_signal[numPart] = {0.51, 1.1205, 1.1205};
const float liminf[numPart] = {0.45, 1.08, 1.08};
const float limsup[numPart] = {0.55, 1.18, 1.18};

// pT binning for K0S, Lambda and Antilambda
double pTbins[numPart][20] = {{0.0, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.6, 2.0, 2.4, 3.0},{0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 2.0, 2.4, 3.0},{0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 2.0, 2.4, 3.0}};
int nbinpT[numPart] = {17,8,8};
// double pTbins[numPart][20] = {{0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.6, 2.0, 2.4, 3.0},{0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 2.0, 2.4, 3.0},{0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 2.0, 2.4, 3.0}};
// int nbinpT[numPart] = {16,8,8};
// double pTbins[numPart][20] = {{0.0, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.6, 2.0, 2.4, 3.0},{0.6, 0.8, 1.0, 1.2, 1.6, 2.0, 2.4, 3.0},{0.6, 0.8, 1.0, 1.2, 1.6, 2.0, 2.4, 3.0}};
// int nbinpT[numPart] = {17,7,7};
// double pTbins[numPart][4] = {{0.6, 0.7},{1., 1.2},{0.6, 1.2}};
// int nbinpT[numPart] = {1,1,1};
// double pTbins[numPart][4] = {{0., 3.},{0., 3.},{0., 3.}};
// int nbinpT[numPart] = {1,1,1};

double pTbinsXi[20] = {0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 2.0, 2.4, 3.0};
int nbinpT_Xi = 8;
// double pTbinsXi[20] = {0.6, 1.4, 2.0, 3.0};
// int nbinpT_Xi = 3;

//For V0 daughter tracking efficiencies
int daughterSpecies = 0; //0 pion, 1 proton
double pTbinsPions[20] = {0.0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.4, 0.5, 0.6, 0.8, 1.0, 1.4, 1.8, 2.2, 2.6};
int nbinpTPions = 15;
double pTbinsProtons[20] = {0.0, 0.1, 0.2, 0.3, 0.5, 0.7, 0.9, 1.2, 1.4, 1.6, 2.0, 3.0, 4.0};
int nbinpTProtons = 12;
double daughterWindow[2] = {pTbinsPions[nbinpTPions], pTbinsProtons[nbinpTProtons]};

// Run List and run-specific files
// const int nRuns = 8;
// const int runList[nRuns] = {505548, 505582, 505600, 505629, 505637, 505645, 505658, 505669};
// const TString runList_string[nRuns] = {"run 505548", "run 505582", "run 505600", "run 505629", "run 505637", "run 505645", "run 505658", "run 505669"};
// TFile* file_O2Analysis1 = new TFile("AnalysisResult_Ref_run505548.root");
// TFile* file_O2Analysis2 = new TFile("AnalysisResult_Ref_run505582.root");
// TFile* file_O2Analysis3 = new TFile("AnalysisResult_Ref_run505600.root");
// TFile* file_O2Analysis4 = new TFile("AnalysisResult_Ref_run505629.root");
// TFile* file_O2Analysis5 = new TFile("AnalysisResult_Ref_run505637.root");
// TFile* file_O2Analysis6 = new TFile("AnalysisResult_Ref_run505645.root");
// TFile* file_O2Analysis7 = new TFile("AnalysisResult_Ref_run505658.root");
// TFile* file_O2Analysis8 = new TFile("AnalysisResult_Ref_run505669.root");

// TFile* file_O2Analysis_RunComparison_array[nRuns] = {file_O2Analysis1, file_O2Analysis2, file_O2Analysis3, file_O2Analysis4, file_O2Analysis5, file_O2Analysis6, file_O2Analysis7, file_O2Analysis8};
// const int colors_Runs[] = {kCyan, kAzure+7, kBlue, kMagenta+3, kPink+5, kRed, kOrange+1, kOrange-2}; // for syst bands
// const int markers_Runs[] = {kFullSquare, kFullSquare, kFullSquare, kFullCross,kFullCross,kFullCross,kFullCross,kFullCross,kFullCross,kFullCross};


const int nRuns = 2; // 8;
const int runList[nRuns] = {505548, 505582};
const TString runList_string[nRuns] = {"bad runs", "good runs"};
// TFile* file_O2Analysis1 = new TFile("AnalysisResult_BadRuns.root");if running the bad/goodrun comparison this is the correct file
// TFile* file_O2Analysis2 = new TFile("AnalysisResult_GoodRuns.root");
TFile* file_O2Analysis1 = new TFile("AnalysisResult_Ref.root"); // placeholders
TFile* file_O2Analysis2 = new TFile("AnalysisResult_Ref.root");

TFile* file_O2Analysis_RunComparison_array[nRuns] = {file_O2Analysis1, file_O2Analysis2}; //, file_O2Analysis3, file_O2Analysis4, file_O2Analysis5, file_O2Analysis6, file_O2Analysis7, file_O2Analysis8};
const int colors_Runs[] = {kBlue, kRed}; // for syst bands
const int markers_Runs[] = {kFullSquare, kFullCircle};





// const float MassPart[numPart] = {0.497611, 1.115683, 1.115683};// 1.32171, 1.32171, 1.67245, 1.67245};
// const float WidthPartLimitsFit[2][numPart] = {{0.003, 0.001, 0.001}, {0.01, 0.005, 0.005}};//{{0.003, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001}, {0.01, 0.005, 0.005, 0.005, 0.005, 0.008, 0.008}};

// === Commonly used x/ titles: ===
// rapidity
TString* texRapidity = new TString("y");
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
TString* texFeeddownRemoved = new TString("Feeddown removed - ratio");

TString* texPseudoEfficiency = new TString("Pseudo efficiency #it{V0}_{detected}/#it{V0}_{expected} (ratio)");

TString* texCountRelativeToMax = new TString("Count normalised to peak maximum");
TString* texInvMass_K0ShortDecay_runs = new TString("M_{#pi^{+}#pi^{-}} (GeV/#it{c}^{2})");


TString* texPtMeasured = new TString("#it{p}_{T}^{measured} (GeV/#it{c})");
TString* texPtMC = new TString("#it{p}_{T}^{MC} (GeV/#it{c})");

const TString* texInvMass_K0ShortDecay = new TString("M_{#pi^{+}#pi^{-}} (GeV/#it{c}^{2})");
const TString* texInvMass_LambdaDecay = new TString("M_{p#pi^{-}} (GeV/#it{c}^{2})");
const TString* texInvMass_AntiLambdaDecay = new TString("M_{#bar{p}#pi^{+}} (GeV/#it{c}^{2})");

const TString* texInvMassDecays_titles[numPart] = {texInvMass_K0ShortDecay,texInvMass_LambdaDecay,texInvMass_AntiLambdaDecay};

const TString* texPt_K0S = new TString("#it{p}_{T}(K^{0}_{S}) (GeV/#it{c})");
const TString* texPt_Lambda = new TString("#it{p}_{T}(#Lambda) (GeV/#it{c})");
const TString* texPt_AntiLambda = new TString("#it{p}_{T}(#bar{#Lambda}) (GeV/#it{c})");
const TString* texPt_particle[numPart] = {texPt_K0S,texPt_Lambda,texPt_AntiLambda};
const TString* texPt_Omega = new TString("#it{p}_{T}(#Omega) (GeV/#it{c})");
const TString* texPt_Xi = new TString("#it{p}_{T}(#Xi) (GeV/#it{c})");

void PaperQualityPlot_PilotBeam_withUnfolding() {
  // Load necessary libraries
  LoadLibs();
  // Set the default style
  SetStyle();

  // Automatically creates error structure when creating new histograms
  // TH1::SetDefaultSumw2(kTRUE); // unfortunately seems to remove error bars on my graphs; not sure why
  // test_SetDefaultSumw2();

  //histograms definition and processing
  TH1D* hstat, *hsyst, *hsystCorr;
  int icolor=0;

  TString* SaveAs_Title = new TString("");
  TString* texXtitle = new TString("");
  TString* texYtitle = new TString("");
  TString* Extra = new TString("");

  // PtDistribution_O2data(hstat, hsyst, hsystCorr, ipart, file_O2Analysis, SaveAs_Title);
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

  // double maxDCAcut = 1;
  // double minTpcCrossedRowscut = 0;
  // double maxTpcCrossedRowscut = 200;
  // TH3D* H3D_detectedV0s_TruePt_DcaHist = (TH3D*)file_O2Analysis->Get("lambdakzero-analysis-mc/h3dMass"+NamePart[ipart]+"Dca");
  // double* DCAbins_fullFileRange = (double*)H3D_detectedV0s_TruePt_DcaHist->GetXaxis();
  // Efficiency_TpcCrossedRowsScan(hstat, hsyst, hsystCorr, ipart, file_O2Analysis, SaveAs_Title, texXtitle, texYtitle, pTbins[ipart], nbinpT[ipart], minTpcCrossedRowscut, maxTpcCrossedRowscut);

  // Efficiency_DcaScan(hstat, hsyst, hsystCorr, ipart, file_O2Analysis, SaveAs_Title, texXtitle, texYtitle, pTbins[ipart], nbinpT[ipart], maxDCAcut);

  // Systematics_CutVariations_Graphs(hstat, hsyst, hsystCorr, ipart, SaveAs_Title, texXtitle, texYtitle, pTbins[ipart], nbinpT[ipart]);
  CorrectedYield_withSystematics(hstat, hsyst, hsystCorr, ipart, SaveAs_Title, texXtitle, texYtitle, pTbins[ipart], nbinpT[ipart]);
  // Draw_FeedDownMatrix(hstat, hsyst, hsystCorr, ipart, file_O2Analysis, SaveAs_Title, texXtitle, texYtitle, pTbins[ipart], nbinpT[ipart]);
  // Draw_XsiPlus_MassPlot(hstat, hsyst, hsystCorr, ipart, file_O2Analysis, SaveAs_Title, texXtitle, texYtitle, pTbins[ipart], nbinpT[ipart]);
  // Draw_FeedDown_LambdaAntiLambda(hstat, hsyst, hsystCorr, ipart, file_O2Analysis, SaveAs_Title, texXtitle, texYtitle, pTbins[ipart], nbinpT[ipart]);
  // Draw_FeedDown_ipart(hstat, hsyst, hsystCorr, ipart, file_O2Analysis, SaveAs_Title, texXtitle, texYtitle, pTbins[ipart], nbinpT[ipart]);
  // Test_LambdaTrueMCYieldvsHEP(hstat, hsyst, hsystCorr, ipart, file_O2Analysis, nCutsIterations, SaveAs_Title, texXtitle, texYtitle, pTbins[ipart], nbinpT[ipart]);
  // Draw_Fit_Mass_Histograms(hstat, hsyst, hsystCorr, file_O2Analysis, SaveAs_Title, texXtitle, texYtitle, ipart, pTbins[ipart], nbinpT[ipart]);
  // Draw_Mass_Histograms(hstat, hsyst, hsystCorr, file_O2Analysis, SaveAs_Title, texXtitle, texYtitle, ipart, pTbins[ipart], nbinpT[ipart]);
  // Draw_Rapidity_MC(hstat, hsyst, hsystCorr, file_O2Analysis, SaveAs_Title, texXtitle, texYtitle, ipart);
  // PseudoEfficiency_O2vsHEP(hstat, hsyst, hsystCorr, ipart, SaveAs_Title, texXtitle, texYtitle, pTbins[ipart], nbinpT[ipart]);
  // Draw_pT_vs_pTmc(hstat, hsyst, hsystCorr, file_O2Analysis, SaveAs_Title, texXtitle, texYtitle, ipart);
  // Draw_FinderEfficiencies(hstat, hsyst, hsystCorr, file_O2Analysis1, file_O2Analysis2, SaveAs_Title, texXtitle, texYtitle, ipart, pTbins[ipart], nbinpT[ipart]);
  // Draw_XiChargedToNeutralRatio(hstat, hsyst, hsystCorr, file_O2Analysis, SaveAs_Title, texXtitle, texYtitle, ipart, pTbins[ipart], nbinpT[ipart]);

  // Test_V0DaugtherPairsTrackingEfficiency(hstat, hsyst, hsystCorr, ipart, file_O2Analysis, SaveAs_Title, texXtitle, texYtitle, pTbins[ipart], nbinpT[ipart]);
  // Test_V0DaugtherPairsTrackingEfficiency_PostDcaFitter(hstat, hsyst, hsystCorr, ipart, file_O2Analysis, SaveAs_Title, texXtitle, texYtitle, pTbins[ipart], nbinpT[ipart]);

  // PtDifferential_SigmaOfFit(hstat, hsyst, hsystCorr, ipart, file_O2Analysis, SaveAs_Title, texXtitle, texYtitle, pTbins[ipart] ,nbinpT[ipart], isMC);
  // InvMass_Plot(hstat, hsyst, hsystCorr, ipart, file_O2Analysis, SaveAs_Title, texXtitle, texYtitle, isMC);
  // pT_Spectrum_postAnalyserCuts(hstat, hsyst, hsystCorr, ipart, file_O2Analysis, SaveAs_Title, texXtitle, texYtitle, pTbins[ipart], nbinpT[ipart], isMC);
  // pT_Spectrum_preAnalyserCutsK0S(hstat, hsyst, hsystCorr, ipart, file_O2Analysis, SaveAs_Title, texXtitle, texYtitle, pTbins[ipart], nbinpT[ipart], isMC);

  //////// PhD Amendments
  // Draw_Mass_Histograms_multipleRuns(hstat, hsyst, hsystCorr, SaveAs_Title, texXtitle, texYtitle, ipart);
  // Draw_EfficiencyWithSigmoidFit(hstat, hsyst, hsystCorr, file_O2Analysis, SaveAs_Title, texXtitle, texYtitle, ipart, pTbins[ipart], nbinpT[ipart]);
  // Draw_PtResponseMatrix(hstat, hsyst, hsystCorr, file_O2Analysis, SaveAs_Title, texXtitle, texYtitle, pTbins[ipart], nbinpT[ipart]);

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
  TH1 * h = cfig->DrawFrame(0,0,3,0.1); // PseudoEfficiency_HEPcomparison
  // TH1 * h = cfig->DrawFrame(0,0,3,0.06); // Efficiency_K0S
  // TH1 * h = cfig->DrawFrame(0,0,3,0.7); // Efficiency_K0S new MC
  // TH1 * h = cfig->DrawFrame(0,0,3,1); // Efficiency Daughter Pairs reco
  // TH1 * h = cfig->DrawFrame(0,0,3,0.25); // CorrectedYield K0S
  // TH1 * h = cfig->DrawFrame(0,0,3,0.2); // Feeddown ratio removed
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

  TLatex * textColl = new TLatex (0.18,0.82,"pp #sqrt{#it{s}} = 900 GeV, pilot beam 2021");
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

void DrawLogo (int logo, double xmin, double ymin) {

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
  // double AliLogo_LowX =xmin;
  // double AliLogo_LowY = ymin;
  // double AliLogo_Height = size;
  // //ALICE logo is a  file that is 821x798 pixels->should be wider than a square
  // double AliLogo_Width  = (821./798.) * AliLogo_Height * gPad->GetWh() / gPad->GetWw();
  
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
  int nbinx = hstat->GetNbinsX();
  
  for(int ibinx = 1; ibinx <= nbinx; ibinx++){
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

void RawSpectrum_O2data_allCountSignalRegion(TH1D* &hstat, TH1D* &hsyst, TH1D*&hsystCorr, int ipart, TFile* file_O2Analysis, TString* SaveAs_Title, TString* &texXtitle, TString* &texYtitle, double* pTbins, int nbinpT) {

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

  // double xbins[9] = {0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.6, 2.0, 3.0};
  // int nbinpT = 8;
  TH1D *H1D_DetectedV0s_O2_SmallpTinterval[nbinpT];
  TH1D *hPtDiff_RawSpectrum_O2 = new TH1D("hPtDiff_RawSpectrum_O2","hPtDiff_RawSpectrum_O2",nbinpT,pTbins);

  for(int ibinPt = 1; ibinPt <= nbinpT; ibinPt++){
    int iPt = ibinPt-1; //different convention for bins that start at 1 and arrays that start at 0

    int ibinPt_originalAxis_low = H3D_DetectedV0s_O2->GetYaxis()->FindBin(pTbins[iPt]);
    int ibinPt_originalAxis_up = H3D_DetectedV0s_O2->GetYaxis()->FindBin(pTbins[iPt+1]) - 1;
    cout << "ibinPt: " << ibinPt << "; originalAxis [low;up] = [" << ibinPt_originalAxis_low << ";" << ibinPt_originalAxis_up << "]" << endl;

    H1D_DetectedV0s_O2_SmallpTinterval[iPt] = (TH1D*)H3D_DetectedV0s_O2->ProjectionZ("RawSpectrum_O2data_allCountSignalRegion_InvMass_"+NamePart[ipart]+Form("_Pt%i", ibinPt),0,-1,ibinPt_originalAxis_low,ibinPt_originalAxis_up);

    float SignalCountRough = H1D_DetectedV0s_O2_SmallpTinterval[iPt]->Integral(MassSignalInterval_lowEdgeBin,MassSignalInterval_upEdgeBin);

    float dpT = hPtDiff_RawSpectrum_O2->GetXaxis()->GetBinWidth(ibinPt);
    float dN_dpT = SignalCountRough *1./dpT;
    hPtDiff_RawSpectrum_O2->SetBinContent(ibinPt,dN_dpT);
    hPtDiff_RawSpectrum_O2->SetBinError(ibinPt,sqrt(SignalCountRough)*1./dpT);
  }

  // H3D_DetectedV0s_O2->RebinY(10);
  // hPtDiff_RawSpectrum_O2->Rebin(10);

  // int nbinPt = hPtDiff_RawSpectrum_O2->GetNbinsX();
  // for(int ibinPt = 1; ibinPt <= nbinPt; ibinPt++){
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
  int nbinx = hstat->GetNbinsX();
  
  for(int ibinx = 1; ibinx <= nbinx; ibinx++){
    hsyst->SetBinContent (ibinx, hstat->GetBinContent(ibinx));
    hsyst->SetBinError   (ibinx, hstat->GetBinContent(ibinx)*0.08);
    hsystCorr->SetBinContent (ibinx, hstat->GetBinContent(ibinx));
    hsystCorr->SetBinError   (ibinx, hstat->GetBinContent(ibinx)*0.15);
  } 

  *SaveAs_Title += "RawSpectrum_O2data_allCountSignalRegion";
  texXtitle = texPtX;
  texYtitle = texFeeddownRemoved;
}

void PtDistribution_O2data(TH1D* &hstat, TH1D* &hsyst, TH1D*&hsystCorr, int ipart, TFile* file_O2Analysis, TString* SaveAs_Title, TString* &texXtitle, TString* &texYtitle) {

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
  int nbinx = hstat->GetNbinsX();
  
  for(int ibinx = 1; ibinx <= nbinx; ibinx++){
    hsyst->SetBinContent (ibinx, hstat->GetBinContent(ibinx));
    hsyst->SetBinError   (ibinx, hstat->GetBinContent(ibinx)*0.08);
    hsystCorr->SetBinContent (ibinx, hstat->GetBinContent(ibinx));
    hsystCorr->SetBinError   (ibinx, hstat->GetBinContent(ibinx)*0.15);
  } 

  *SaveAs_Title += "Pt_Distrib_O2data";
  texXtitle = texPtX;
  texYtitle = texRawYield;
}

void InvMass_Plot(TH1D* &hstat, TH1D* &hsyst, TH1D*&hsystCorr, int ipart, TFile* file_O2Analysis, TString* SaveAs_Title, TString* &texXtitle, TString* &texYtitle, bool isMC) {

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
  // float ScalingFactor = SelectedEventCount;
  hstat->Rebin(10);   

  // float ScalingFactor = hstat->GetMaximum();  
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
  int nbinx = hstat->GetNbinsX();
  
  for(int ibinx = 1; ibinx <= nbinx; ibinx++){
    hsyst->SetBinContent (ibinx, hstat->GetBinContent(ibinx));
    hsyst->SetBinError   (ibinx, hstat->GetBinContent(ibinx)*0.08);
    hsystCorr->SetBinContent (ibinx, hstat->GetBinContent(ibinx));
    hsystCorr->SetBinError   (ibinx, hstat->GetBinContent(ibinx)*0.15);
  } 

  *SaveAs_Title += "InvMass";
  texXtitle = texPtX;
  texYtitle = texCount;
}

void DcaHistoProcessing_O2data(TH1D* &hstat, TH1D* &hsyst, TH1D*&hsystCorr, int ipart, TFile* file_O2Analysis) {

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
  int nbinx = hstat->GetNbinsX();
  
  for(int ibinx = 1; ibinx <= nbinx; ibinx++){
    hsyst->SetBinContent (ibinx, hstat->GetBinContent(ibinx));
    hsyst->SetBinError   (ibinx, hstat->GetBinContent(ibinx)*0.08);
    hsystCorr->SetBinContent (ibinx, hstat->GetBinContent(ibinx));
    hsystCorr->SetBinError   (ibinx, hstat->GetBinContent(ibinx)*0.15);
  } 

}

void PseudoEfficiency_HEPcomparison_allCountSignalRegionOLD_withBugs(TH1D* &hstat, TH1D* &hsyst, TH1D*&hsystCorr, int ipart, TFile* file_O2Analysis, TString* SaveAs_Title, double* pTbins, int nbinpT) {
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

  int nbinPt = hPtDiff_RawSpectrum_O2->GetNbinsX();
  for(int ibinPt = 1; ibinPt <= nbinPt; ibinPt++){
    float SignalCountRough = H3D_DetectedV0s_O2->Integral(0,-1,ibinPt,ibinPt+1,MassSignalInterval_lowEdgeBin,MassSignalInterval_upEdgeBin);

    hPtDiff_RawSpectrum_O2->SetBinContent(ibinPt,SignalCountRough); //hPtDiff_RawSpectrum_O2 is actually a count, not the 1/Nev*dN/dpT
  }

  double xbins[17] = {0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.6, 2.0, 2.4, 3.0};
  TH1D* hPtDiff_RawSpectrum_O2_rebinned = (TH1D*)hPtDiff_RawSpectrum_O2->Rebin(16,"hPtDiff_RawSpectrum_O2_rebinned",xbins);
  hPtDiff_RawSpectrum_O2_rebinned->Sumw2();

  double placeholder = pTbins[0];
  int placeholder2 = nbinpT;

  // ->Divide(hPtDiff_RawSpectrum_HEP);
  int nbinpT_2 = hPtDiff_RawSpectrum_O2_rebinned->GetNbinsX();
  for(int ibinPt = 1; ibinPt <= nbinpT_2; ibinPt++){
    float SignalCount_ibinPt = hPtDiff_RawSpectrum_O2_rebinned->GetBinContent(ibinPt);
    float dpT = hPtDiff_RawSpectrum_O2_rebinned->GetXaxis()->GetBinWidth(ibinPt);
    // float ratio_HEP_O2 = (SignalCount_ibinPt*1./dpT)*1./hPtDiff_RawSpectrum_HEP->GetBinContent(ibinPt);// = dN/dpT_O2 / dN/dpT_HEP)
    hPtDiff_RawSpectrum_O2_rebinned->SetBinContent(ibinPt,SignalCount_ibinPt*1./dpT);

    hPtDiff_RawSpectrum_O2_rebinned->SetBinError(ibinPt,sqrt(SignalCount_ibinPt)*1./dpT);//*1./hPtDiff_RawSpectrum_HEP->GetBinContent(ibinPt)
    hPtDiff_RawSpectrum_HEP[ipart]->SetBinError(ibinPt,hPtDiff_RawSpectrum_HEP_statErrors[ipart]->GetBinContent(ibinPt));
  }

// 'MIGHT BE MISSING a division by rapidity interval for HEP DATA because it is dN**2/dpT/dy and rapidity interval is [-0.75,0.75]'

  hstat = (TH1D*)hPtDiff_RawSpectrum_O2->Rebin(16,"hstat",xbins);
  hPtDiff_RawSpectrum_O2_rebinned->Scale(1./SelectedEventCount);
  hstat->Reset("M");
  hstat->Divide(hPtDiff_RawSpectrum_O2_rebinned,hPtDiff_RawSpectrum_HEP[ipart], 1., 1., "b");

  //Error Bars Processing
  int nbinx = hstat->GetNbinsX();

  for(int ibinx = 1; ibinx <= nbinx; ibinx++){
    hsyst->SetBinContent (ibinx, hstat->GetBinContent(ibinx));
    hsyst->SetBinError   (ibinx, hstat->GetBinContent(ibinx)*0.08);
    hsystCorr->SetBinContent (ibinx, hstat->GetBinContent(ibinx));
    hsystCorr->SetBinError   (ibinx, hstat->GetBinContent(ibinx)*0.15);
  } 

  *SaveAs_Title += "PseudoEfficiency_HEPcomparison_allCountSignalRegion";
}

void PseudoEfficiency_HEPcomparison_allCountSignalRegion(TH1D* &hstat, TH1D* &hsyst, TH1D*&hsystCorr, int ipart, TFile* file_O2Analysis, TString* SaveAs_Title, TString* &texXtitle, TString* &texYtitle, double* pTbins, int nbinpT) {
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
  // int nbinPt = hPtDiff_RawSpectrum_O2->GetNbinsX();

  for(int ibinPt = 1; ibinPt <= nbinpT; ibinPt++){
    int iPt = ibinPt-1; //different convention for bins that start at 1 and arrays that start at 0

    int ibinPt_originalAxis_low = H3D_DetectedV0s_O2->GetYaxis()->FindBin(pTbins[iPt]);
    int ibinPt_originalAxis_up = H3D_DetectedV0s_O2->GetYaxis()->FindBin(pTbins[iPt+1]) - 1;
    cout << "ibinPt: " << ibinPt << "; originalAxis [low;up] = [" << ibinPt_originalAxis_low << ";" << ibinPt_originalAxis_up << "]" << endl;

    H1D_DetectedV0s_O2_SmallpTinterval[iPt] = (TH1D*)H3D_DetectedV0s_O2->ProjectionZ("PseudoEfficiency_HEPcomparison_allCountSignalRegion_InvMass_"+NamePart[ipart]+Form("_Pt%i", ibinPt),0,-1,ibinPt_originalAxis_low,ibinPt_originalAxis_up);

    // for(int ibinPt = 1; ibinPt <= nbinPt; ibinPt++){

    float SignalCountRough = H1D_DetectedV0s_O2_SmallpTinterval[iPt]->Integral(MassSignalInterval_lowEdgeBin,MassSignalInterval_upEdgeBin);

    hPtDiff_RawSpectrum_O2->SetBinContent(ibinPt,SignalCountRough); //hPtDiff_RawSpectrum_O2 is actually a count, not the 1/Nev*dN/dpT

    double dpT = hPtDiff_RawSpectrum_O2->GetXaxis()->GetBinWidth(ibinPt);
    float drapidity = 1.5;//1.5
    double dN_dpT = SignalCountRough *1./dpT*1./drapidity;
    hPtDiff_RawSpectrum_O2->SetBinContent(ibinPt,dN_dpT);
    hPtDiff_RawSpectrum_O2->SetBinError(ibinPt,sqrt(SignalCountRough)*1./dpT*1./drapidity); //Issue with errors; this error of the O2 data spectrum is false, not right formula with fit
    cout << "(ibinPt,dN_dpT) = (" << ibinPt << "," << dN_dpT*1./SelectedEventCount << ")" << endl;

    hPtDiff_RawSpectrum_HEP[ipart]->SetBinError(ibinPt,hPtDiff_RawSpectrum_HEP_statErrors[ipart]->GetBinContent(ibinPt));
  }

  // double xbins[17] = {0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.6, 2.0, 2.4, 3.0};
  // TH1D* hPtDiff_RawSpectrum_O2_rebinned = (TH1D*)hPtDiff_RawSpectrum_O2->Rebin(16,"hPtDiff_RawSpectrum_O2_rebinned",xbins);
  // hPtDiff_RawSpectrum_O2_rebinned->Sumw2();

  // double placeholder = pTbins[0];
  // int placeholder2 = nbinpT;

  // // ->Divide(hPtDiff_RawSpectrum_HEP);
  // int nbinpT_2 = hPtDiff_RawSpectrum_O2_rebinned->GetNbinsX();
  // for(int ibinPt = 1; ibinPt <= nbinpT_2; ibinPt++){
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
  int nbinx = hstat->GetNbinsX();

  for(int ibinx = 1; ibinx <= nbinx; ibinx++){
    hsyst->SetBinContent (ibinx, hstat->GetBinContent(ibinx));
    hsyst->SetBinError   (ibinx, hstat->GetBinContent(ibinx)*0.08);
    hsystCorr->SetBinContent (ibinx, hstat->GetBinContent(ibinx));
    hsystCorr->SetBinError   (ibinx, hstat->GetBinContent(ibinx)*0.15);
  } 

  *SaveAs_Title += "PseudoEfficiency_HEPcomparison_allCountSignalRegion";
  texXtitle = texPtX;
  texYtitle = texptDifferentialYield_HEPratio_pseudoEff;
}



double fparab(double *x, double *par) {
  const int numPart=7;
  // float liminf[numPart]={0.48, 1.11, 1.11, 1.31,  1.31,  1.665, 1.665};
  // float limsup[numPart]={0.54, 1.13, 1.13, 1.335, 1.335, 1.685, 1.685}; Original
  // float limsup[numPart]={0.515, 1.13, 1.13, 1.335, 1.335, 1.685, 1.685}; //Tweaked
  int part=par[3];
  if (reject && x[0] > min_range_signal[part] && x[0] < max_range_signal[part]) { //Original was using this duplicate but different definition of limsup and liminf
    TF1::RejectPoint();
    return 0;
  }
  return par[0] + par[1]*x[0] + par[2]*x[0]*x[0];
}
double fline(double *x, double *par) {
  const int numPart=7;
  // float liminf[numPart]={0.48, 1.11, 1.11, 1.31,  1.31,  1.665, 1.665};
  // float limsup[numPart]={0.54, 1.13, 1.13, 1.335, 1.335, 1.685, 1.685}; Original
  // float limsup[numPart]={0.515, 1.13, 1.13, 1.335, 1.335, 1.685, 1.685}; //Tweaked
  int part=par[2];
  if (reject && x[0] > min_range_signal[part] && x[0] < max_range_signal[part]) { //Original was using this duplicate but different definition of limsup and liminf
    TF1::RejectPoint();
    return 0;
  }
  return par[0] + par[1]*x[0];
}

void RawSpectrum_O2data_withFit(TH1D* &hstat, TH1D* &hsyst, TH1D*&hsystCorr, int ipart, TFile* file_O2Analysis, TString* SaveAs_Title, TString* &texXtitle, TString* &texYtitle, double* pTbins, int nbinpT) {

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
  float fFitResult[nbinpT];
  double parGaussParab[nbinpT][6];  

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

    H1D_DetectedV0s_O2_SmallpTinterval[iPt] = (TH1D*)H3D_DetectedV0s_O2_rebinnedZ->ProjectionZ("RawSpectrum_O2data_withFit_InvMass_"+NamePart[ipart]+Form("_Pt%i", ibinPt),0,-1,ibinPt_originalAxis_low,ibinPt_originalAxis_up);

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
    hPtDiff_RawSpectrum_O2->SetBinError(ibinPt,sqrt(SignalCountFit)*1./dpT);
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
  int nbinx = hstat->GetNbinsX();
  
  for(int ibinx = 1; ibinx <= nbinx; ibinx++){
    hsyst->SetBinContent (ibinx, hstat->GetBinContent(ibinx));
    hsyst->SetBinError   (ibinx, hstat->GetBinContent(ibinx)*0.08);
    hsystCorr->SetBinContent (ibinx, hstat->GetBinContent(ibinx));
    hsystCorr->SetBinError   (ibinx, hstat->GetBinContent(ibinx)*0.15);
  } 

  *SaveAs_Title += "RawSpectrum_O2data_withFit";
  texXtitle = texPtX;
  texYtitle = texRawYield;
}

void InvMassDistributions_withFit(TH1D* &hstat, TH1D* &hsyst, TH1D* &hsystCorr, int ipart, TFile* file_O2Analysis, TString* SaveAs_Title, TString* &texXtitle, TString* &texYtitle, double* pTbins, int nbinpT) {

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
  
  float fFitResult[nbinpT];
  double parGaussParab[nbinpT][6];  
  double parGaussParabError_mu[nbinpT];  
  double parGaussParabError_sigma[nbinpT];  

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

    H1D_DetectedV0s_O2_SmallpTinterval[iPt] = (TH1D*)H3D_DetectedV0s_O2_rebinnedZ->ProjectionZ("InvMassDistributions_withFit_InvMass_"+NamePart[ipart]+Form("_Pt%i", ibinPt),0,-1,ibinPt_originalAxis_low,ibinPt_originalAxis_up);
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
    int icolor=0;
    H1D_DetectedV0s_O2_SmallpTinterval[iPt]->SetMarkerStyle(markers[0]);
    H1D_DetectedV0s_O2_SmallpTinterval[iPt]->SetMarkerColor(colors [icolor]);
    H1D_DetectedV0s_O2_SmallpTinterval[iPt]->SetLineColor  (colors [icolor]);
    // bkgparab[iPt]->Draw("same");
    // gauss[iPt]->Draw("same");


    // TLatex * text_pTbin = new TLatex (0.62,0.75,Form("pT [%.1f;%.1f] GeV", pTbins[iPt], pTbins[iPt+1]));
    // text_pTbin->SetNDC(kTRUE);
    // text_pTbin->Draw();
    // TLatex * textContext = new TLatex (0.18,0.82,"#bf{ALICE Performance}"); //BOLD
    // textContext->SetTextSize(0.05);
    // textContext->SetNDC(kTRUE); //remove if I want x,y in TLatex to be in the coordinate system of the histogram
    // textContext->Draw();
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
    hPtDiff_RawSpectrum_O2->SetBinError(ibinPt,sqrt(SignalCountFit)*1./dpT*1./drapidity); //Issue with errors; this error of the O2 data spectrum is false, not right formula with fit
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
  int nbinx = hstat->GetNbinsX();
  
  for(int ibinx = 1; ibinx <= nbinx; ibinx++){
    hsyst->SetBinContent (ibinx, hstat->GetBinContent(ibinx));
    hsyst->SetBinError   (ibinx, hstat->GetBinContent(ibinx)*0.08);
    hsystCorr->SetBinContent (ibinx, hstat->GetBinContent(ibinx));
    hsystCorr->SetBinError   (ibinx, hstat->GetBinContent(ibinx)*0.15);
  }

  *SaveAs_Title += "ToBeDeleted";
  texXtitle = texPtX;
  texYtitle = texptDifferentialYield_HEPratio_pseudoEff;

}

void InvMassDistributions_withFit_MC_datalike(TH1D* &hstat, TH1D* &hsyst, TH1D* &hsystCorr, int ipart, TFile* file_O2Analysis, TString* SaveAs_Title, TString* &texXtitle, TString* &texYtitle, double* pTbins, int nbinpT) {

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
  
  float fFitResult[nbinpT];
  // double parGaussParab[nbinpT][6]; //parabola background
  double parGaussParab[nbinpT][5]; //line background
  double parGaussParabError_mu[nbinpT];  
  double parGaussParabError_sigma[nbinpT];  

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

    H1D_DetectedV0s_O2_SmallpTinterval[iPt] = (TH1D*)H3D_DetectedV0s_O2_rebinnedZ->ProjectionZ("InvMassDistributions_withFit_MC_datalike_InvMass_"+NamePart[ipart]+Form("_Pt%i", ibinPt),0,-1,ibinPt_originalAxis_low,ibinPt_originalAxis_up);
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
    int icolor=0;
    H1D_DetectedV0s_O2_SmallpTinterval[iPt]->SetMarkerStyle(markers[0]);
    H1D_DetectedV0s_O2_SmallpTinterval[iPt]->SetMarkerColor(colors [icolor]);
    H1D_DetectedV0s_O2_SmallpTinterval[iPt]->SetLineColor  (colors [icolor]);
    bkgparab[iPt]->Draw("same");
    // gauss[iPt]->Draw("same");
    GaussPlusPolynom[iPt]->Draw("same");

    // TLatex * text_pTbin = new TLatex (0.62,0.75,Form("pT [%.1f;%.1f] GeV", pTbins[iPt], pTbins[iPt+1]));
    // text_pTbin->SetNDC(kTRUE);
    // text_pTbin->Draw();
    // TLatex * textContext = new TLatex (0.18,0.82,"#bf{ALICE Performance}"); //BOLD
    // textContext->SetTextSize(0.05);
    // textContext->SetNDC(kTRUE); //remove if I want x,y in TLatex to be in the coordinate system of the histogram
    // textContext->Draw();
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
    hPtDiff_RawSpectrum_O2->SetBinError(ibinPt,sqrt(SignalCountFit)*1./dpT*1./drapidity); //Issue with errors; this error of the O2 data spectrum is false, not right formula with fit
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
  int nbinx = hstat->GetNbinsX();
  
  for(int ibinx = 1; ibinx <= nbinx; ibinx++){
    hsyst->SetBinContent (ibinx, hstat->GetBinContent(ibinx));
    hsyst->SetBinError   (ibinx, hstat->GetBinContent(ibinx)*0.08);
    hsystCorr->SetBinContent (ibinx, hstat->GetBinContent(ibinx));
    hsystCorr->SetBinError   (ibinx, hstat->GetBinContent(ibinx)*0.15);
  }

  *SaveAs_Title += "ToBeDeleted";
  texXtitle = texPtX;
  texYtitle = texptDifferentialYield_HEPratio_pseudoEff;

}

void PtDifferential_Fit(TH1D* &hPtDiff_RawSpectrum_O2, TH3D* H3D_DetectedV0s_O2, TH3D* H3D_DetectedV0s_O2_rebinnedZ, int ipart, double* pTbins, int nbinpT, double* parGaussianParab_Sigma, double* parGaussianParab_Sigma_error) {

  //Fit tools initialisation
  TF1 *gauss[nbinpT];
  TF1 *GaussPlusPolynom[nbinpT];
  TF1 *bkgparab[nbinpT];
  
  float fFitResult[nbinpT];
  double parGaussParab[nbinpT][6];  

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

    H1D_DetectedV0s_O2_SmallpTinterval[iPt] = (TH1D*)H3D_DetectedV0s_O2_rebinnedZ->ProjectionZ("PtDifferential_Fit_InvMass_"+NamePart[ipart]+Form("_Pt%i", ibinPt),0,-1,ibinPt_originalAxis_low,ibinPt_originalAxis_up);

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
    int icolor=0;
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
    hPtDiff_RawSpectrum_O2->SetBinError(ibinPt,sqrt(SignalCountFit)*1./dpT*1./drapidity); //Issue with errors; this error of the O2 data spectrum is false, not right formula with fit
  }

  canvasMassvsPt->SaveAs("Fits_invMass_ptDiff_"+NamePart[ipart]+".pdf");
}

void PtDifferential_SigmaOfFit(TH1D* &hstat, TH1D* &hsyst, TH1D* &hsystCorr, int ipart, TFile* file_O2Analysis, TString* SaveAs_Title, TString* &texXtitle, TString* &texYtitle, double* pTbins, int nbinpT, bool isMC) {

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

  double parGaussianParab_Sigma[nbinpT];  
  double parGaussianParab_Sigma_error[nbinpT];  
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
  int nbinx = hstat->GetNbinsX();
  
  for(int ibinx = 1; ibinx <= nbinx; ibinx++){
    hsyst->SetBinContent (ibinx, hstat->GetBinContent(ibinx));
    hsyst->SetBinError   (ibinx, hstat->GetBinContent(ibinx)*0.08);
    hsystCorr->SetBinContent (ibinx, hstat->GetBinContent(ibinx));
    hsystCorr->SetBinError   (ibinx, hstat->GetBinContent(ibinx)*0.15);
  }

  *SaveAs_Title += "SigmaFit_pTdiff";
  texXtitle = texPtX;
  texYtitle = texSgimaGaussFit;
}




void Efficiency_O2data_allCountSignalRegion(TH1D* &hstat, TH1D* &hsyst, TH1D*&hsystCorr, int ipart, TFile* file_O2Analysis, TString* SaveAs_Title, TString* &texXtitle, TString* &texYtitle, double* pTbins, int nbinpT) {

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

  // double xbins[9] = {0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.6, 2.0, 3.0};
  // int nbinpT = 8;
  TH1D *H1D_DetectedV0s_O2_SmallpTinterval[nbinpT];
  TH1D *hPtDiff_RawSpectrum_O2 = new TH1D("hPtDiff_RawSpectrum_O2","hPtDiff_RawSpectrum_O2",nbinpT,pTbins);
  TH1D *TrueV0PtSpectrum_rebinned = new TH1D("TrueV0PtSpectrum_rebinned","TrueV0PtSpectrum_rebinned",nbinpT,pTbins);

  for(int ibinPt = 1; ibinPt <= nbinpT; ibinPt++){
    int iPt = ibinPt-1; //different convention for bins that start at 1 and arrays that start at 0

    int ibinPt_originalAxis_low = H3D_DetectedV0s_O2->GetYaxis()->FindBin(pTbins[iPt]);
    int ibinPt_originalAxis_up = H3D_DetectedV0s_O2->GetYaxis()->FindBin(pTbins[iPt+1]) - 1;
    cout << "ibinPt: " << ibinPt << "; originalAxis [low;up] = [" << ibinPt_originalAxis_low << ";" << ibinPt_originalAxis_up << "]" << endl;

    H1D_DetectedV0s_O2_SmallpTinterval[iPt] = (TH1D*)H3D_DetectedV0s_O2->ProjectionZ("Efficiency_O2data_allCountSignalRegion_InvMass_"+NamePart[ipart]+Form("_Pt%i", ibinPt),0,-1,ibinPt_originalAxis_low,ibinPt_originalAxis_up);

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

  // int nbinPt = hPtDiff_RawSpectrum_O2->GetNbinsX();
  // for(int ibinPt = 1; ibinPt <= nbinPt; ibinPt++){
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
  int nbinx = hstat->GetNbinsX();
  
  for(int ibinx = 1; ibinx <= nbinx; ibinx++){
    hsyst->SetBinContent (ibinx, hstat->GetBinContent(ibinx));
    hsyst->SetBinError   (ibinx, hstat->GetBinContent(ibinx)*0.0000001);
    hsystCorr->SetBinContent (ibinx, hstat->GetBinContent(ibinx));
    hsystCorr->SetBinError   (ibinx, hstat->GetBinContent(ibinx)*0.15);
  } 

  *SaveAs_Title += "Efficiency_O2data_allCountSignalRegion";
  texXtitle = texPtX;
  texYtitle = texEfficiency;
}



void Efficiency_O2MCdata_withFit(TH1D* &hstat, TH1D* &hsyst, TH1D*&hsystCorr, int ipart, TFile* file_O2Analysis, TString* SaveAs_Title, TString* &texXtitle, TString* &texYtitle, double* pTbins, int nbinpT) {

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
  float fFitResult[nbinpT];
  double parGaussParab[nbinpT][6];  

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

    H1D_DetectedV0s_O2_SmallpTinterval[iPt] = (TH1D*)H3D_DetectedV0s_O2_rebinnedZ->ProjectionZ("Efficiency_O2MCdata_withFit_InvMass_"+NamePart[ipart]+Form("_Pt%i", ibinPt),0,-1,ibinPt_originalAxis_low,ibinPt_originalAxis_up);

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
    int icolor=0;
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
  int nbinx = hstat->GetNbinsX();
  
  for(int ibinx = 1; ibinx <= nbinx; ibinx++){
    hsyst->SetBinContent (ibinx, hstat->GetBinContent(ibinx));
    hsyst->SetBinError   (ibinx, hstat->GetBinContent(ibinx)*0.08);
    hsystCorr->SetBinContent (ibinx, hstat->GetBinContent(ibinx));
    hsystCorr->SetBinError   (ibinx, hstat->GetBinContent(ibinx)*0.15);
  } 

  *SaveAs_Title += "Efficiency_O2data_withFit";
  texXtitle = texPtX;
  texYtitle = texEfficiency;
}


void Efficiency_O2data_truePt_trueV0s(TH1D* &hstat, TH1D* &hsyst, TH1D*&hsystCorr, int ipart, TFile* file_O2Analysis, TString* SaveAs_Title, TString* &texXtitle, TString* &texYtitle, double* pTbins, int nbinpT) {

  //O2 file
  // TFile* file_O2Analysis = new TFile("AnalysisResults_WithSel8Cut.root");

  TH1I* H1I_SelectedEventCount = (TH1I*)file_O2Analysis->Get("lambdakzero-analysis-mc/hSelectedEventCounter");
  float SelectedEventCount = H1I_SelectedEventCount->GetEntries();

  TH3D* H3D_DetectedV0s_O2 = (TH3D*)file_O2Analysis->Get("lambdakzero-analysis-mc/h3dMass"+NamePart[ipart]+"_MC_truePt");
  TH1D* H1D_DetectedV0s_PtAxis_O2 = (TH1D*)H3D_DetectedV0s_O2->ProjectionY("InvMass_K0s",0,-1,0,-1);
  // TH1D* hPtDiff_RawSpectrum_O2 = (TH1D*)H1D_DetectedV0s_PtAxis_O2->Clone("hPtDiff_RawSpectrum");
  // hPtDiff_RawSpectrum_O2->Reset("M"); //should be empty, only used that to get same Pt binning as used in analysis
 
  TH1D* TrueV0PtSpectrum = (TH1D*)file_O2Analysis->Get("lambdakzero-particle-count-mc/h"+NamePart[ipart]+"Count_PtDiff");

  // double xbins[9] = {0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.6, 2.0, 3.0};
  // int nbinpT = 8;
  TH1D *H1D_DetectedV0s_O2_SmallpTinterval[nbinpT];
  TH1D *hPtDiff_RawSpectrum_O2 = new TH1D("hPtDiff_RawSpectrum_O2","hPtDiff_RawSpectrum_O2",nbinpT,pTbins);
  TH1D *TrueV0PtSpectrum_rebinned = new TH1D("TrueV0PtSpectrum_rebinned","TrueV0PtSpectrum_rebinned",nbinpT,pTbins);

  for(int ibinPt = 1; ibinPt <= nbinpT; ibinPt++){
    int iPt = ibinPt-1; //different convention for bins that start at 1 and arrays that start at 0

    int ibinPt_originalAxis_low = H3D_DetectedV0s_O2->GetYaxis()->FindBin(pTbins[iPt]);
    int ibinPt_originalAxis_up = H3D_DetectedV0s_O2->GetYaxis()->FindBin(pTbins[iPt+1]) - 1;
    cout << "ibinPt: " << ibinPt << "; originalAxis [low;up] = [" << ibinPt_originalAxis_low << ";" << ibinPt_originalAxis_up << "]" << endl;

    H1D_DetectedV0s_O2_SmallpTinterval[iPt] = (TH1D*)H3D_DetectedV0s_O2->ProjectionZ("Efficiency_O2data_truePt_trueV0s_InvMass_"+NamePart[ipart]+Form("_Pt%i", ibinPt),0,-1,ibinPt_originalAxis_low,ibinPt_originalAxis_up);

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

  // int nbinPt = hPtDiff_RawSpectrum_O2->GetNbinsX();
  // for(int ibinPt = 1; ibinPt <= nbinPt; ibinPt++){
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
  int nbinx = hstat->GetNbinsX();
  
  for(int ibinx = 1; ibinx <= nbinx; ibinx++){
    hsystCorr->SetBinContent (ibinx, hstat->GetBinContent(ibinx));
    hsystCorr->SetBinError   (ibinx, hstat->GetBinContent(ibinx)*0.15);
  } 

  *SaveAs_Title += "Efficiency_O2data_truePt_trueV0s";
  texXtitle = texPtX;
  texYtitle = texEfficiency;
}

void RawSpectrum_O2data_truePt_trueV0s(TH1D* &hstat, TH1D* &hsyst, TH1D*&hsystCorr, int ipart, TFile* file_O2Analysis, TString* SaveAs_Title, TString* &texXtitle, TString* &texYtitle, double* pTbins, int nbinpT) {

  //O2 file
  // TFile* file_O2Analysis = new TFile("AnalysisResults_WithSel8Cut.root");

  TH1I* H1I_SelectedEventCount = (TH1I*)file_O2Analysis->Get("lambdakzero-analysis-mc/hSelectedEventCounter");
  float SelectedEventCount = H1I_SelectedEventCount->GetEntries();

  TH3D* H3D_DetectedV0s_O2 = (TH3D*)file_O2Analysis->Get("lambdakzero-analysis-mc/h3dMass"+NamePart[ipart]+"_MC_truePt");
  TH1D* H1D_DetectedV0s_PtAxis_O2 = (TH1D*)H3D_DetectedV0s_O2->ProjectionY("InvMass_K0s",0,-1,0,-1);
  // TH1D* hPtDiff_RawSpectrum_O2 = (TH1D*)H1D_DetectedV0s_PtAxis_O2->Clone("hPtDiff_RawSpectrum");
  // hPtDiff_RawSpectrum_O2->Reset("M"); //should be empty, only used that to get same Pt binning as used in analysis
 
  TH1D* TrueV0PtSpectrum = (TH1D*)file_O2Analysis->Get("lambdakzero-particle-count-mc/h"+NamePart[ipart]+"Count_PtDiff");

  // double xbins[9] = {0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.6, 2.0, 3.0};
  // int nbinpT = 8;
  TH1D *H1D_DetectedV0s_O2_SmallpTinterval[nbinpT];
  TH1D *hPtDiff_RawSpectrum_O2 = new TH1D("hPtDiff_RawSpectrum_O2","hPtDiff_RawSpectrum_O2",nbinpT,pTbins);
  TH1D *TrueV0PtSpectrum_rebinned = new TH1D("TrueV0PtSpectrum_rebinned","TrueV0PtSpectrum_rebinned",nbinpT,pTbins);

  for(int ibinPt = 1; ibinPt <= nbinpT; ibinPt++){
    int iPt = ibinPt-1; //different convention for bins that start at 1 and arrays that start at 0

    int ibinPt_originalAxis_low = H3D_DetectedV0s_O2->GetYaxis()->FindBin(pTbins[iPt]);
    int ibinPt_originalAxis_up = H3D_DetectedV0s_O2->GetYaxis()->FindBin(pTbins[iPt+1]) - 1;
    cout << "ibinPt: " << ibinPt << "; originalAxis [low;up] = [" << ibinPt_originalAxis_low << ";" << ibinPt_originalAxis_up << "]" << endl;

    H1D_DetectedV0s_O2_SmallpTinterval[iPt] = (TH1D*)H3D_DetectedV0s_O2->ProjectionZ("RawSpectrum_O2data_truePt_trueV0s_InvMass_"+NamePart[ipart]+Form("_Pt%i", ibinPt),0,-1,ibinPt_originalAxis_low,ibinPt_originalAxis_up);

    float SignalCountRough = H1D_DetectedV0s_O2_SmallpTinterval[iPt]->Integral(1,H1D_DetectedV0s_O2_SmallpTinterval[iPt]->GetNbinsX());

    float dpT = hPtDiff_RawSpectrum_O2->GetXaxis()->GetBinWidth(ibinPt);
    float dN_dpT = SignalCountRough *1./dpT;
    hPtDiff_RawSpectrum_O2->SetBinContent(ibinPt,dN_dpT);
    hPtDiff_RawSpectrum_O2->SetBinError(ibinPt,sqrt(SignalCountRough)*1./dpT);

    int trueV0count_ibinpT = TrueV0PtSpectrum->Integral(ibinPt_originalAxis_low,ibinPt_originalAxis_up);
    TrueV0PtSpectrum_rebinned->SetBinContent(ibinPt,trueV0count_ibinpT);
    TrueV0PtSpectrum_rebinned->SetBinError(ibinPt,sqrt(trueV0count_ibinpT));

    cout << "(ibinPt, SignalCountRough, trueV0count_ibinpT) = (" << ibinPt << ", " << SignalCountRough << ", " << trueV0count_ibinpT << ")" << endl;
    cout << "yield = " << SignalCountRough*1./dpT*1./SelectedEventCount << endl;
  }

  // H3D_DetectedV0s_O2->RebinY(10);
  // hPtDiff_RawSpectrum_O2->Rebin(10);

  // int nbinPt = hPtDiff_RawSpectrum_O2->GetNbinsX();
  // for(int ibinPt = 1; ibinPt <= nbinPt; ibinPt++){
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
  int nbinx = hstat->GetNbinsX();
  
  for(int ibinx = 1; ibinx <= nbinx; ibinx++){
    hsyst->SetBinContent (ibinx, hstat->GetBinContent(ibinx));
    hsyst->SetBinError   (ibinx, hstat->GetBinContent(ibinx)*0.08);
    hsystCorr->SetBinContent (ibinx, hstat->GetBinContent(ibinx));
    hsystCorr->SetBinError   (ibinx, hstat->GetBinContent(ibinx)*0.15);
  } 

  *SaveAs_Title += "RawSpectrum_O2data_truePt_trueV0s";
  texXtitle = texPtX;
  texYtitle = texEfficiency;
}




void V0Daughters_TrackingEfficiencies(TH1D* &hstat, TH1D* &hsyst, TH1D*&hsystCorr, int ipart, TFile* file_O2Analysis, TString* SaveAs_Title, double* pTbins, int nbinpT) {

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
  int icolor=0;

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
  int nbinx = hstat->GetNbinsX();
  
  for(int ibinx = 1; ibinx <= nbinx; ibinx++){
    hsyst->SetBinContent (ibinx, hstat->GetBinContent(ibinx));
    hsyst->SetBinError   (ibinx, hstat->GetBinContent(ibinx)*0.08);
    hsystCorr->SetBinContent (ibinx, hstat->GetBinContent(ibinx));
    hsystCorr->SetBinError   (ibinx, hstat->GetBinContent(ibinx)*0.15);
  } 

  *SaveAs_Title += "ToBeDeleted";
}

void V0Daughters_TrackingEfficiencies_HistogramProcessing(TH1D* &daughter1, TH1D* &daughter2, int ipart, TFile* file_O2Analysis, double* pTbinsProtons, int nbinpTProtons, double* pTbinsPions, int nbinpTPions) {

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


void V0Daughters_TrackingEfficiencies_Daughter1(TH1D* &hstat, TH1D* &hsyst, TH1D*&hsystCorr, int ipart, TFile* file_O2Analysis, TString* SaveAs_Title, TString* &texXtitle, TString* &texYtitle, TString* &Extra, double* pTbinsProtons, int nbinpTProtons, double* pTbinsPions, int nbinpTPions) {

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
  int nbinx = hstat->GetNbinsX();
  
  for(int ibinx = 1; ibinx <= nbinx; ibinx++){
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

void V0Daughters_TrackingEfficiencies_Daughter2(TH1D* &hstat, TH1D* &hsyst, TH1D*&hsystCorr, int ipart, TFile* file_O2Analysis, TString* SaveAs_Title, TString* &texXtitle, TString* &texYtitle, TString* &Extra, double* pTbinsProtons, int nbinpTProtons, double* pTbinsPions, int nbinpTPions) {

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
  int nbinx = hstat->GetNbinsX();
  
  for(int ibinx = 1; ibinx <= nbinx; ibinx++){
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




void pT_Spectrum_postAnalyserCuts(TH1D* &hstat, TH1D* &hsyst, TH1D*&hsystCorr, int ipart, TFile* file_O2Analysis, TString* SaveAs_Title, TString* &texXtitle, TString* &texYtitle, double* pTbins, int nbinpT, bool isMC) {

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
    H1D_DetectedV0s_PtAxis_O2_rebin->SetBinError(ibinPt,sqrt(Count)*1./dpT);
  }

  H1D_DetectedV0s_PtAxis_O2_rebin->Scale(1./SelectedEventCount);
  hstat = (TH1D*)H1D_DetectedV0s_PtAxis_O2_rebin->Clone("hstat");


  //////Error Bars///////
  hsyst     = new TH1D("hsyst", "hsyst", 100, 0, 10);
  hsystCorr = new TH1D("hsystCorr", "hsystCorr", 100, 0, 10);
  hstat     ->Sumw2(); 
  hsyst     ->Sumw2();
  hsystCorr ->Sumw2(); 
  int nbinx = hstat->GetNbinsX();
  
  for(int ibinx = 1; ibinx <= nbinx; ibinx++){
    hsyst->SetBinContent (ibinx, hstat->GetBinContent(ibinx));
    hsyst->SetBinError   (ibinx, hstat->GetBinContent(ibinx)*0.08);
    hsystCorr->SetBinContent (ibinx, hstat->GetBinContent(ibinx));
    hsystCorr->SetBinError   (ibinx, hstat->GetBinContent(ibinx)*0.15);
  } 

  *SaveAs_Title += "pT_Spectrum_DirtySignal_postAnalyserCuts";
  texXtitle = texPtX;
  texYtitle = texptDifferentialYield;
}

void pT_Spectrum_preAnalyserCutsK0S(TH1D* &hstat, TH1D* &hsyst, TH1D*&hsystCorr, int ipart, TFile* file_O2Analysis, TString* SaveAs_Title, TString* &texXtitle, TString* &texYtitle, double* pTbins, int nbinpT, bool isMC) {

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
    H1D_DetectedV0s_PtAxis_O2_rebin->SetBinError(ibinPt,sqrt(Count)*1./dpT);
  }

  H1D_DetectedV0s_PtAxis_O2_rebin->Scale(1./SelectedEventCount);
  hstat = (TH1D*)H1D_DetectedV0s_PtAxis_O2_rebin->Clone("hstat");


  //////Error Bars///////
  hsyst     = new TH1D("hsyst", "hsyst", 100, 0, 10);
  hsystCorr = new TH1D("hsystCorr", "hsystCorr", 100, 0, 10);
  hstat     ->Sumw2(); 
  hsyst     ->Sumw2();
  hsystCorr ->Sumw2(); 
  int nbinx = hstat->GetNbinsX();
  
  for(int ibinx = 1; ibinx <= nbinx; ibinx++){
    hsyst->SetBinContent (ibinx, hstat->GetBinContent(ibinx));
    hsyst->SetBinError   (ibinx, hstat->GetBinContent(ibinx)*0.08);
    hsystCorr->SetBinContent (ibinx, hstat->GetBinContent(ibinx));
    hsystCorr->SetBinError   (ibinx, hstat->GetBinContent(ibinx)*0.15);
  } 

  *SaveAs_Title += "pT_Spectrum_DirtySignal_preAnalyserCuts";
  texXtitle = texPtX;
  texYtitle = texptDifferentialYield;
}










// void Get_McEfficiency_vsPt(TH1D* &hMcEfficiency_vsPt, TH3D* H3D_detectedV0s_TruePt, TH1D* TrueV0PtSpectrum, int ipart, double* pTbins, int nbinpT, int ibinXaxisCut_low, int ibinXaxisCut_high, double SideBandSizeMultiplierModifier, float SignalMean, float SignalStandardDev) {

//   // what is fed as input: TH3D* H3D_detectedV0s_TruePt = (TH3D*)file_O2Analysis->Get("lambdakzero-analysis-mc/h3dMass"+NamePart[ipart]+"_MC_truePt");
//   // what is fed as input: TH1D* TrueV0PtSpectrum = (TH1D*)file_O2Analysis->Get("lambdakzero-particle-count-mc/h"+NamePart[ipart]+"Count_PtDiff");

//   TH1D *H1D_detectedV0s_TruePt_SmallpTinterval[nbinpT];
//   TH1D *hMcSignalCount_vsPt = new TH1D("hMcSignalCount_vsPt","hMcSignalCount_vsPt",nbinpT,pTbins);
//   TH1D *TrueV0PtSpectrum_AnalysisBins = new TH1D("TrueV0PtSpectrum_AnalysisBins","TrueV0PtSpectrum_AnalysisBins",nbinpT,pTbins);

//   for(int ibinPt = 1; ibinPt <= nbinpT; ibinPt++){
//     int iPt = ibinPt-1; //different convention for bins that start at 1 and arrays that start at 0

//     int ibinPt_originalAxis_low = H3D_detectedV0s_TruePt->GetYaxis()->FindBin(pTbins[iPt]);
//     int ibinPt_originalAxis_up = H3D_detectedV0s_TruePt->GetYaxis()->FindBin(pTbins[iPt+1]) - 1;
//     cout << "ibinPt: " << ibinPt << "; originalAxis [low;up] = [" << ibinPt_originalAxis_low << ";" << ibinPt_originalAxis_up << "]" << endl;

//     H1D_detectedV0s_TruePt_SmallpTinterval[iPt] = (TH1D*)H3D_detectedV0s_TruePt->ProjectionZ("EfficiencyFunction_InvMass_"+NamePart[ipart]+Form("_Pt%i", ibinPt),ibinXaxisCut_low,ibinXaxisCut_high,ibinPt_originalAxis_low,ibinPt_originalAxis_up);

//     float SignalCountFull = H1D_detectedV0s_TruePt_SmallpTinterval[iPt]->Integral(1,H1D_detectedV0s_TruePt_SmallpTinterval[iPt]->GetNbinsX()); // avoids overflow bins 0 and -1 ; without fit

//     int ibinMass_signalRange_low = H1D_detectedV0s_TruePt_SmallpTinterval[iPt]->GetXaxis()->FindBin(SignalMean - (DefaultSideBandSize_array[ipart] + SideBandSizeMultiplierModifier)*SignalStandardDev);
//     int ibinMass_signalRange_up = H1D_detectedV0s_TruePt_SmallpTinterval[iPt]->GetXaxis()->FindBin(SignalMean + (DefaultSideBandSize_array[ipart] + SideBandSizeMultiplierModifier)*SignalStandardDev);
//     cout << "(ibinMass_signalRange_low, ibinMass_signalRange_up) = (" << ibinMass_signalRange_low << ", " << ibinMass_signalRange_up << ")" << endl;

//     float SignalCount_SignalRange = H1D_detectedV0s_TruePt_SmallpTinterval[iPt]->Integral(ibinMass_signalRange_low, ibinMass_signalRange_up); // avoids overflow bins 0 and -1 ; without fit


//     hMcSignalCount_vsPt->SetBinContent(ibinPt,SignalCount_SignalRange);
//     hMcSignalCount_vsPt->SetBinError(ibinPt,sqrt(SignalCount_SignalRange));
//     // hMcSignalCount_vsPt->SetBinContent(ibinPt,SignalCountFull); // for test purposes
//     // hMcSignalCount_vsPt->SetBinError(ibinPt,sqrt(SignalCountFull));

//     int trueV0count_ibinpT = TrueV0PtSpectrum->Integral(ibinPt_originalAxis_low,ibinPt_originalAxis_up);
//     TrueV0PtSpectrum_AnalysisBins->SetBinContent(ibinPt,trueV0count_ibinpT);
//     TrueV0PtSpectrum_AnalysisBins->SetBinError(ibinPt,sqrt(trueV0count_ibinpT));

//     cout << "(ibinPt, SignalCountFull, SignalCount_SignalRange, trueV0count_ibinpT) = (" << ibinPt << ", " << SignalCountFull << ", " << SignalCount_SignalRange << ", " << trueV0count_ibinpT << ")" << endl;
//     cout << "efficiency = " << SignalCount_SignalRange*1./trueV0count_ibinpT << endl;
//   }

//   hMcEfficiency_vsPt = (TH1D*)hMcSignalCount_vsPt->Clone("hMcEfficiency_vsPt");
//   hMcEfficiency_vsPt->Divide(hMcSignalCount_vsPt,TrueV0PtSpectrum_AnalysisBins, 1., 1., "b"); // option b for binomial because efficiency: https://twiki.cern.ch/twiki/bin/view/ALICE/PWGLFPAGSTRANGENESSEfficiency
//   delete hMcSignalCount_vsPt; // avoids memory leaks when Get_McEfficiency_vsPt is called in a loop
//   delete TrueV0PtSpectrum_AnalysisBins; // avoids memory leaks when Get_McEfficiency_vsPt is called in a loop

// }

void Get_McEfficiency_vsPt(TH1D* &hMcEfficiency_vsPt, TH3D* H3D_detectedV0s_TruePt, TH1D* TrueV0PtSpectrum, int ipart, double* pTbins, int nbinpT, int ibinXaxisCut_low, int ibinXaxisCut_high, double SideBandSizeMultiplierModifier, float SignalMean, float SignalStandardDev, int InvMassRebinFactor) {

  // what is fed as input: TH3D* H3D_detectedV0s_TruePt = (TH3D*)file_O2Analysis->Get("lambdakzero-analysis-mc/h3dMass"+NamePart[ipart]+"_MC_truePt");
  // what is fed as input: TH1D* TrueV0PtSpectrum = (TH1D*)file_O2Analysis->Get("lambdakzero-particle-count-mc/h"+NamePart[ipart]+"Count_PtDiff");
  TH1D *H3D_detectedV0s_TruePt_rebinnedZ = (TH1D*)H3D_detectedV0s_TruePt->RebinZ(InvMassRebinFactor,"H1D_detectedV0s_TruePt_SmallpTinterval_rebinned"); 

  TH1D *H1D_detectedV0s_TruePt_SmallpTinterval[nbinpT];
  TH1D *hMcSignalCount_vsPt = new TH1D("hMcSignalCount_vsPt","hMcSignalCount_vsPt",nbinpT,pTbins);
  TH1D *TrueV0PtSpectrum_AnalysisBins = new TH1D("TrueV0PtSpectrum_AnalysisBins","TrueV0PtSpectrum_AnalysisBins",nbinpT,pTbins);

  for(int ibinPt = 1; ibinPt <= nbinpT; ibinPt++){
    int iPt = ibinPt-1; //different convention for bins that start at 1 and arrays that start at 0

    int ibinPt_originalAxis_low = H3D_detectedV0s_TruePt->GetYaxis()->FindBin(pTbins[iPt]);
    int ibinPt_originalAxis_up = H3D_detectedV0s_TruePt->GetYaxis()->FindBin(pTbins[iPt+1]) - 1;
    cout << "ibinPt: " << ibinPt << "; originalAxis [low;up] = [" << ibinPt_originalAxis_low << ";" << ibinPt_originalAxis_up << "]" << endl;

    H1D_detectedV0s_TruePt_SmallpTinterval[iPt] = (TH1D*)H3D_detectedV0s_TruePt->ProjectionZ("EfficiencyFunction_InvMass_"+NamePart[ipart]+Form("_Pt%i", ibinPt),ibinXaxisCut_low,ibinXaxisCut_high,ibinPt_originalAxis_low,ibinPt_originalAxis_up);

    float SignalCountFull = H1D_detectedV0s_TruePt_SmallpTinterval[iPt]->Integral(1,H1D_detectedV0s_TruePt_SmallpTinterval[iPt]->GetNbinsX()); // avoids overflow bins 0 and -1 ; without fit

    int ibinMass_signalRange_low = H1D_detectedV0s_TruePt_SmallpTinterval[iPt]->GetXaxis()->FindBin(SignalMean - (DefaultSideBandSize_array[ipart] + SideBandSizeMultiplierModifier)*SignalStandardDev);
    int ibinMass_signalRange_up = H1D_detectedV0s_TruePt_SmallpTinterval[iPt]->GetXaxis()->FindBin(SignalMean + (DefaultSideBandSize_array[ipart] + SideBandSizeMultiplierModifier)*SignalStandardDev);
    cout << "(ibinMass_signalRange_low, ibinMass_signalRange_up) = (" << ibinMass_signalRange_low << ", " << ibinMass_signalRange_up << ")" << endl;

    float SignalCount_SignalRange = H1D_detectedV0s_TruePt_SmallpTinterval[iPt]->Integral(ibinMass_signalRange_low, ibinMass_signalRange_up); // avoids overflow bins 0 and -1 ; without fit


    hMcSignalCount_vsPt->SetBinContent(ibinPt,SignalCount_SignalRange);
    hMcSignalCount_vsPt->SetBinError(ibinPt,sqrt(SignalCount_SignalRange));
    // hMcSignalCount_vsPt->SetBinContent(ibinPt,SignalCountFull); // for test purposes
    // hMcSignalCount_vsPt->SetBinError(ibinPt,sqrt(SignalCountFull));

    int trueV0count_ibinpT = TrueV0PtSpectrum->Integral(ibinPt_originalAxis_low,ibinPt_originalAxis_up);
    TrueV0PtSpectrum_AnalysisBins->SetBinContent(ibinPt,trueV0count_ibinpT);
    TrueV0PtSpectrum_AnalysisBins->SetBinError(ibinPt,sqrt(trueV0count_ibinpT));

    cout << "(ibinPt, SignalCountFull, SignalCount_SignalRange, trueV0count_ibinpT) = (" << ibinPt << ", " << SignalCountFull << ", " << SignalCount_SignalRange << ", " << trueV0count_ibinpT << ")" << endl;
    cout << "efficiency = " << SignalCount_SignalRange*1./trueV0count_ibinpT << endl;
  }

  hMcEfficiency_vsPt = (TH1D*)hMcSignalCount_vsPt->Clone("hMcEfficiency_vsPt");
  hMcEfficiency_vsPt->Divide(hMcSignalCount_vsPt,TrueV0PtSpectrum_AnalysisBins, 1., 1., "b"); // option b for binomial because efficiency: https://twiki.cern.ch/twiki/bin/view/ALICE/PWGLFPAGSTRANGENESSEfficiency
  delete hMcSignalCount_vsPt; // avoids memory leaks when Get_McEfficiency_vsPt is called in a loop
  delete TrueV0PtSpectrum_AnalysisBins; // avoids memory leaks when Get_McEfficiency_vsPt is called in a loop

  if (doFitEfficiencySmooth) {
    ////// fit with generalised logistics
    // float a = 1.7;
    float K = 1.2;
    float Q = 0.46;
    float B = 2.7;
    float Nu = 20;
    cout << "DEBUG 7" <<endl;
    TF1* sigmoidFit_Efficiency = new TF1("sigmoidFit_Efficiency", fgenerallogistic, 0, 3, 2);
    sigmoidFit_Efficiency->SetParameter(0, hMcEfficiency_vsPt->GetMaximum()); //*H1I_SelectedEventCount_SVCuts->GetEntries()
    sigmoidFit_Efficiency->SetParameter(1, Q);
    // sigmoidFit_Efficiency->SetParameter(2, c);
    // sigmoidFit_Efficiency->SetParameter(3, d);
    // sigmoidFit_Efficiency->SetParameter(4, e);
    sigmoidFit_Efficiency->SetParLimits(0, 0., 1.5*hMcEfficiency_vsPt->GetMaximum());
    // sigmoidFit_Efficiency->SetParLimits(1, 0, 10);
    // sigmoidFit_Efficiency->SetParLimits(2, 0, 10);
    // sigmoidFit_Efficiency->SetParLimits(3, 0, 10);

    cout << "DEBUG 8" <<endl;
    TFitResultPtr fFitResult_neutral = hMcEfficiency_vsPt->Fit(sigmoidFit_Efficiency,"SME+R0"); // the option G helps avoid the "Warning in <Fit>: Abnormal termination of minimization" warning and combined with ME gives "ERROR MATRIX ACCURATE " instead of "ERROR MATRIX UNCERTAINTY  24.6 per cent" or "ERROR NOT POS. DEFINED"

    // TF1* sigmoidFit_Efficiency_errorPlus = new TF1("sigmoidFit_Efficiency_errorPlus", fgenerallogistic, 0, 3, 2);
    // TF1* sigmoidFit_Efficiency_errorMinus = new TF1("sigmoidFit_Efficiency_errorMinus", fgenerallogistic, 0, 3, 2);

    // sigmoidFit_Efficiency_errorPlus->SetParameter(0, sigmoidFit_Efficiency->GetParameter(0) + sigmoidFit_Efficiency->GetParError(0)); //*H1I_SelectedEventCount_SVCuts->GetEntries()
    // sigmoidFit_Efficiency_errorPlus->SetParameter(1, sigmoidFit_Efficiency->GetParameter(1) + sigmoidFit_Efficiency->GetParError(1));
    // sigmoidFit_Efficiency_errorMinus->SetParameter(0, sigmoidFit_Efficiency->GetParameter(0) - sigmoidFit_Efficiency->GetParError(0)); //*H1I_SelectedEventCount_SVCuts->GetEntries()
    // sigmoidFit_Efficiency_errorMinus->SetParameter(1, sigmoidFit_Efficiency->GetParameter(1) - sigmoidFit_Efficiency->GetParError(1));
    double x, dist_from_histEff;
    for(int ibinPt = 1; ibinPt <= nbinpT; ibinPt++){
      x = hMcEfficiency_vsPt->GetBinCenter(ibinPt);
      // hMcEfficiency_vsPt->SetBinContent(ibinPt, sigmoidFit_Efficiency->Eval(hMcEfficiency_vsPt->GetBinCenter(ibinPt)));
      // double fitted_efficiency_errorPlus = abs(abs(sigmoidFit_Efficiency->Eval(hMcEfficiency_vsPt->GetBinCenter(ibinPt))) - abs(sigmoidFit_Efficiency_errorPlus->Eval(hMcEfficiency_vsPt->GetBinCenter(ibinPt))));
      // double fitted_efficiency_errorMinus = abs(abs(sigmoidFit_Efficiency->Eval(hMcEfficiency_vsPt->GetBinCenter(ibinPt))) - abs(sigmoidFit_Efficiency_errorMinus->Eval(hMcEfficiency_vsPt->GetBinCenter(ibinPt))));
      // hMcEfficiency_vsPt->SetBinError(ibinPt, max(fitted_efficiency_errorPlus, fitted_efficiency_errorMinus));
      // cout << "ibinPt = " << ibinPt << ", hMcEfficiency_vsPt error+ = " << fitted_efficiency_errorPlus << ", Eff(params) = " << sigmoidFit_Efficiency->Eval(hMcEfficiency_vsPt->GetBinCenter(ibinPt)) << ", Eff(params-dparam) = " << sigmoidFit_Efficiency_errorPlus->Eval(hMcEfficiency_vsPt->GetBinCenter(ibinPt)) << endl;
      // cout << "         " << ibinPt << ", hMcEfficiency_vsPt error- = " << fitted_efficiency_errorMinus << ", Eff(params) = " << sigmoidFit_Efficiency->Eval(hMcEfficiency_vsPt->GetBinCenter(ibinPt)) << ", Eff(params-dparam) = " << sigmoidFit_Efficiency_errorMinus->Eval(hMcEfficiency_vsPt->GetBinCenter(ibinPt)) << endl;
      x = hMcEfficiency_vsPt->GetBinCenter(ibinPt);

      // TMatrixDSym covMatrix = fFitResult->GetCovarianceMatrix();
      // K = sigmoidFit_Efficiency->GetParameter(0);
      // Q = sigmoidFit_Efficiency->GetParameter(1);
      // double dsigmoid_dK = 1./pow(1 + Q*exp(-3*x),10);
      // double dsigmoid_dQ = -10 * K * exp(-3*x) /pow(1 + Q*exp(-3*x),11);
      // double fitted_efficiency_error = sqrt(dsigmoid_dK*dsigmoid_dK * covMatrix(0,0) + dsigmoid_dQ*dsigmoid_dQ * covMatrix(1,1) + 2*dsigmoid_dK*dsigmoid_dQ * covMatrix(0,1));
      dist_from_histEff = abs(abs(sigmoidFit_Efficiency->Eval(x)) - abs(hMcEfficiency_vsPt->GetBinContent(ibinPt)));
      // cout << "ibinPt = " << ibinPt << ", fitted_efficiency_error = " << fitted_efficiency_error << endl;
      // cout << "           dist_from_histEff = " << dist_from_histEff << endl;

      hMcEfficiency_vsPt->SetBinContent(ibinPt, sigmoidFit_Efficiency->Eval(x));
      hMcEfficiency_vsPt->SetBinError(ibinPt, dist_from_histEff);
    }
  }


  // and now I should try and fit it with the sigmoid
}
void Get_McEfficiency_vsPt_alwaysFitSmooth(TH1D* &hMcEfficiency_vsPt, TH3D* H3D_detectedV0s_TruePt, TH1D* TrueV0PtSpectrum, int ipart, double* pTbins, int nbinpT, int ibinXaxisCut_low, int ibinXaxisCut_high, double SideBandSizeMultiplierModifier, float SignalMean, float SignalStandardDev, int InvMassRebinFactor) {

  // what is fed as input: TH3D* H3D_detectedV0s_TruePt = (TH3D*)file_O2Analysis->Get("lambdakzero-analysis-mc/h3dMass"+NamePart[ipart]+"_MC_truePt");
  // what is fed as input: TH1D* TrueV0PtSpectrum = (TH1D*)file_O2Analysis->Get("lambdakzero-particle-count-mc/h"+NamePart[ipart]+"Count_PtDiff");
  TH1D *H3D_detectedV0s_TruePt_rebinnedZ = (TH1D*)H3D_detectedV0s_TruePt->RebinZ(InvMassRebinFactor,"H1D_detectedV0s_TruePt_SmallpTinterval_rebinned"); 

  TH1D *H1D_detectedV0s_TruePt_SmallpTinterval[nbinpT];
  TH1D *hMcSignalCount_vsPt = new TH1D("hMcSignalCount_vsPt","hMcSignalCount_vsPt",nbinpT,pTbins);
  TH1D *TrueV0PtSpectrum_AnalysisBins = new TH1D("TrueV0PtSpectrum_AnalysisBins","TrueV0PtSpectrum_AnalysisBins",nbinpT,pTbins);

  for(int ibinPt = 1; ibinPt <= nbinpT; ibinPt++){
    int iPt = ibinPt-1; //different convention for bins that start at 1 and arrays that start at 0

    int ibinPt_originalAxis_low = H3D_detectedV0s_TruePt->GetYaxis()->FindBin(pTbins[iPt]);
    int ibinPt_originalAxis_up = H3D_detectedV0s_TruePt->GetYaxis()->FindBin(pTbins[iPt+1]) - 1;
    cout << "ibinPt: " << ibinPt << "; originalAxis [low;up] = [" << ibinPt_originalAxis_low << ";" << ibinPt_originalAxis_up << "]" << endl;

    H1D_detectedV0s_TruePt_SmallpTinterval[iPt] = (TH1D*)H3D_detectedV0s_TruePt->ProjectionZ("EfficiencyFunction_InvMass_"+NamePart[ipart]+Form("_Pt%i", ibinPt),ibinXaxisCut_low,ibinXaxisCut_high,ibinPt_originalAxis_low,ibinPt_originalAxis_up);

    float SignalCountFull = H1D_detectedV0s_TruePt_SmallpTinterval[iPt]->Integral(1,H1D_detectedV0s_TruePt_SmallpTinterval[iPt]->GetNbinsX()); // avoids overflow bins 0 and -1 ; without fit

    int ibinMass_signalRange_low = H1D_detectedV0s_TruePt_SmallpTinterval[iPt]->GetXaxis()->FindBin(SignalMean - (DefaultSideBandSize_array[ipart] + SideBandSizeMultiplierModifier)*SignalStandardDev);
    int ibinMass_signalRange_up = H1D_detectedV0s_TruePt_SmallpTinterval[iPt]->GetXaxis()->FindBin(SignalMean + (DefaultSideBandSize_array[ipart] + SideBandSizeMultiplierModifier)*SignalStandardDev);
    cout << "(ibinMass_signalRange_low, ibinMass_signalRange_up) = (" << ibinMass_signalRange_low << ", " << ibinMass_signalRange_up << ")" << endl;

    float SignalCount_SignalRange = H1D_detectedV0s_TruePt_SmallpTinterval[iPt]->Integral(ibinMass_signalRange_low, ibinMass_signalRange_up); // avoids overflow bins 0 and -1 ; without fit


    hMcSignalCount_vsPt->SetBinContent(ibinPt,SignalCount_SignalRange);
    hMcSignalCount_vsPt->SetBinError(ibinPt,sqrt(SignalCount_SignalRange));
    // hMcSignalCount_vsPt->SetBinContent(ibinPt,SignalCountFull); // for test purposes
    // hMcSignalCount_vsPt->SetBinError(ibinPt,sqrt(SignalCountFull));

    int trueV0count_ibinpT = TrueV0PtSpectrum->Integral(ibinPt_originalAxis_low,ibinPt_originalAxis_up);
    TrueV0PtSpectrum_AnalysisBins->SetBinContent(ibinPt,trueV0count_ibinpT);
    TrueV0PtSpectrum_AnalysisBins->SetBinError(ibinPt,sqrt(trueV0count_ibinpT));

    cout << "(ibinPt, SignalCountFull, SignalCount_SignalRange, trueV0count_ibinpT) = (" << ibinPt << ", " << SignalCountFull << ", " << SignalCount_SignalRange << ", " << trueV0count_ibinpT << ")" << endl;
    cout << "efficiency = " << SignalCount_SignalRange*1./trueV0count_ibinpT << endl;
  }

  hMcEfficiency_vsPt = (TH1D*)hMcSignalCount_vsPt->Clone("hMcEfficiency_vsPt");
  hMcEfficiency_vsPt->Divide(hMcSignalCount_vsPt,TrueV0PtSpectrum_AnalysisBins, 1., 1., "b"); // option b for binomial because efficiency: https://twiki.cern.ch/twiki/bin/view/ALICE/PWGLFPAGSTRANGENESSEfficiency
  delete hMcSignalCount_vsPt; // avoids memory leaks when Get_McEfficiency_vsPt is called in a loop
  delete TrueV0PtSpectrum_AnalysisBins; // avoids memory leaks when Get_McEfficiency_vsPt is called in a loop

  ////// fit with generalised logistics
  // float a = 1.7;
  float K = 1.2;
  float Q = 0.46;
  float B = 2.7;
  float Nu = 20;
  cout << "DEBUG 7" <<endl;
  TF1* sigmoidFit_Efficiency = new TF1("sigmoidFit_Efficiency", fgenerallogistic, 0, 3, 2);
  sigmoidFit_Efficiency->SetParameter(0, hMcEfficiency_vsPt->GetMaximum()); //*H1I_SelectedEventCount_SVCuts->GetEntries()
  sigmoidFit_Efficiency->SetParameter(1, Q);
  // sigmoidFit_Efficiency->SetParameter(2, c);
  // sigmoidFit_Efficiency->SetParameter(3, d);
  // sigmoidFit_Efficiency->SetParameter(4, e);
  sigmoidFit_Efficiency->SetParLimits(0, 0., 1.5*hMcEfficiency_vsPt->GetMaximum());
  // sigmoidFit_Efficiency->SetParLimits(1, 0, 10);
  // sigmoidFit_Efficiency->SetParLimits(2, 0, 10);
  // sigmoidFit_Efficiency->SetParLimits(3, 0, 10);

  cout << "DEBUG 8" <<endl;
  TFitResultPtr fFitResult = hMcEfficiency_vsPt->Fit(sigmoidFit_Efficiency,"MES+R0"); // the option G helps avoid the "Warning in <Fit>: Abnormal termination of minimization" warning and combined with ME gives "ERROR MATRIX ACCURATE " instead of "ERROR MATRIX UNCERTAINTY  24.6 per cent" or "ERROR NOT POS. DEFINED"

  TF1* sigmoidFit_Efficiency_errorPlus = new TF1("sigmoidFit_Efficiency_errorPlus", fgenerallogistic, 0, 3, 2);
  TF1* sigmoidFit_Efficiency_errorMinus = new TF1("sigmoidFit_Efficiency_errorMinus", fgenerallogistic, 0, 3, 2);

  sigmoidFit_Efficiency_errorPlus->SetParameter(0, sigmoidFit_Efficiency->GetParameter(0) + sigmoidFit_Efficiency->GetParError(0)); //*H1I_SelectedEventCount_SVCuts->GetEntries()
  sigmoidFit_Efficiency_errorPlus->SetParameter(1, sigmoidFit_Efficiency->GetParameter(1) + sigmoidFit_Efficiency->GetParError(1));
  // sigmoidFit_Efficiency_errorPlus->SetParameter(2, sigmoidFit_Efficiency->GetParameter(2) + sigmoidFit_Efficiency->GetParError(3));
  // sigmoidFit_Efficiency_errorPlus->SetParameter(3, sigmoidFit_Efficiency->GetParameter(3) + sigmoidFit_Efficiency->GetParError(4));
  sigmoidFit_Efficiency_errorMinus->SetParameter(0, sigmoidFit_Efficiency->GetParameter(0) - sigmoidFit_Efficiency->GetParError(0)); //*H1I_SelectedEventCount_SVCuts->GetEntries()
  sigmoidFit_Efficiency_errorMinus->SetParameter(1, sigmoidFit_Efficiency->GetParameter(1) - sigmoidFit_Efficiency->GetParError(1));
  // sigmoidFit_Efficiency_errorMinus->SetParameter(2, sigmoidFit_Efficiency->GetParameter(2) - sigmoidFit_Efficiency->GetParError(3));
  // sigmoidFit_Efficiency_errorMinus->SetParameter(3, sigmoidFit_Efficiency->GetParameter(3) - sigmoidFit_Efficiency->GetParError(4));  // sigmoidFit_Efficiency_error->SetParameter(4, sigmoidFit_Efficiency->GetParameter(4));
  double x;
  // TMatrixDSym covMatrix;
  for(int ibinPt = 1; ibinPt <= nbinpT; ibinPt++){
    x = hMcEfficiency_vsPt->GetBinCenter(ibinPt);

    // TMatrixDSym covMatrix = fFitResult->GetCovarianceMatrix();
    // K = sigmoidFit_Efficiency->GetParameter(0);
    // Q = sigmoidFit_Efficiency->GetParameter(1);
    // double dsigmoid_dK = 1./pow(1 + Q*exp(-3*x),10);
    // double dsigmoid_dQ = -10 * K * exp(-3*x) /pow(1 + Q*exp(-3*x),11);
    // double fitted_efficiency_error = sqrt(dsigmoid_dK*dsigmoid_dK * covMatrix(0,0) + dsigmoid_dQ*dsigmoid_dQ * covMatrix(1,1) + 2*dsigmoid_dK*dsigmoid_dQ * covMatrix(0,1));
    double dist_from_histEff = abs(abs(sigmoidFit_Efficiency->Eval(x)) - abs(hMcEfficiency_vsPt->GetBinContent(ibinPt)));
    // cout << "ibinPt = " << ibinPt << ", fitted_efficiency_error = " << fitted_efficiency_error << endl;
    // cout << "           dist_from_histEff = " << dist_from_histEff << endl;

    hMcEfficiency_vsPt->SetBinContent(ibinPt, sigmoidFit_Efficiency->Eval(x));
    hMcEfficiency_vsPt->SetBinError(ibinPt, dist_from_histEff);

  }




  // and now I should try and fit it with the sigmoid
}

void Get_McSignalBackgroundRatio_vsPt(TH1D* &hSignalBackgroundRatio_vsPt, TH3D* H3D_detectedV0s_DataLike, int ipart, double* pTbins, int nbinpT, int ibinXaxisCut_low,  int ibinXaxisCut_high) {
  //Looks within the signal area of the fit, 6Sigma

  // what is fed as input: TH3D* H3D_detectedV0s_DataLike = (TH3D*)file_O2Analysis->Get("lambdakzero-analysis-mc/h3dMass"+NamePart[ipart]);

  ////////////////////////////////// Fit Start //////////////////////////////////
  //Fit tools initialisation
  TF1 *gauss[nbinpT];
  TF1 *GaussPlusPolynom[nbinpT];
  TF1 *bkgparab[nbinpT];
  float fFitResult[nbinpT];
  double parGaussParab[nbinpT][6];  

  TH1D *H1D_detectedV0s_DataLike_SmallpTinterval[nbinpT];
  TH1D *hSignalCountFit_vsPt = new TH1D("hSignalCountFit_vsPt","hSignalCountFit_vsPt",nbinpT,pTbins);
  TH1D *hBackgroundCountFit_vsPt = new TH1D("hBackgroundCountFit_vsPt","hBackgroundCountFit_vsPt",nbinpT,pTbins);

  for(int ibinPt = 1; ibinPt <= nbinpT; ibinPt++){
    int iPt = ibinPt-1; //different convention for bins that start at 1 and arrays that start at 0

    int ibinPt_originalAxis_low = H3D_detectedV0s_DataLike->GetYaxis()->FindBin(pTbins[iPt]);
    int ibinPt_originalAxis_up = H3D_detectedV0s_DataLike->GetYaxis()->FindBin(pTbins[iPt+1]) - 1;
    cout << "ibinPt: " << ibinPt << "; originalAxis [low;up] = [" << ibinPt_originalAxis_low << ";" << ibinPt_originalAxis_up << "]" << endl;

    H1D_detectedV0s_DataLike_SmallpTinterval[iPt] = (TH1D*)H3D_detectedV0s_DataLike->ProjectionZ("Get_McSignalBackgroundRatio_vsPt_function_InvMass_"+NamePart[ipart]+Form("_Pt%i", ibinPt),ibinXaxisCut_low,ibinXaxisCut_high,ibinPt_originalAxis_low,ibinPt_originalAxis_up);

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

std::array<double, 3> GetTotalAndBackgroundCount_SidebandBackgroundExtrapolation(TH1D* H1D_detectedV0s_DataLike_SmallpTinterval, double SignalMean, double SignalStandardDev, double SideBandSizeMultiplierModifier) {
  // iwth our larger peaks due to lower pagnetic field, we can't go higher than 8sigma for the lower background bound for lambda (due to cutoff from pi+p masses), or 10sigma for k0s (with current window limit)
  // 8 sigma limit is 2*4, so we sample the signal in [-4sigma;+4sigma] interval
  // for a start, the sidebandModifier should be floats, ints is too large a change
  // could have 2*4, 2*3.5, 2*3? 3 is probably the lowest we can go lest we start missing lot of signal (99.9% of guassian area is within -3;+3 sigma)
  double LeftBackgroundLowerBound = SignalMean - 2*(DefaultSideBandSize_array[ipart] + SideBandSizeMultiplierModifier)*SignalStandardDev;
  double LeftBackgroundUpperBound = SignalMean - (DefaultSideBandSize_array[ipart] + SideBandSizeMultiplierModifier)*SignalStandardDev;
  double RightBackgroundLowerBound = SignalMean + (DefaultSideBandSize_array[ipart] + SideBandSizeMultiplierModifier)*SignalStandardDev;
  double RightBackgroundUpperBound = SignalMean + 2*(DefaultSideBandSize_array[ipart] + SideBandSizeMultiplierModifier)*SignalStandardDev;
  cout <<"Sidebands: (" << LeftBackgroundLowerBound << ", " << LeftBackgroundUpperBound << ", " << RightBackgroundLowerBound << ", " << RightBackgroundUpperBound << ")" << endl;

  int LeftBackgroundLowerBound_bin = H1D_detectedV0s_DataLike_SmallpTinterval->GetXaxis()->FindBin(LeftBackgroundLowerBound);
  int LeftBackgroundUpperBound_bin = H1D_detectedV0s_DataLike_SmallpTinterval->GetXaxis()->FindBin(LeftBackgroundUpperBound);
  int RightBackgroundLowerBound_bin = H1D_detectedV0s_DataLike_SmallpTinterval->GetXaxis()->FindBin(RightBackgroundLowerBound);
  int RightBackgroundUpperBound_bin = H1D_detectedV0s_DataLike_SmallpTinterval->GetXaxis()->FindBin(RightBackgroundUpperBound);

  Long_t TotalCountInSignalRegion = H1D_detectedV0s_DataLike_SmallpTinterval->Integral(LeftBackgroundUpperBound_bin, RightBackgroundLowerBound_bin);

  // samples background in the two sidebands chosen: here mean-12sigma to mean-6sigma and mean+6sigma to mean+12sigma for sideBandmodifier=0
  // only do a count of the histogram, not an actual intgral of the background fit function to avoid issues with background line function going below 0 values due to very low background
  Long_t BackgroundCountInSidebands = H1D_detectedV0s_DataLike_SmallpTinterval->Integral(LeftBackgroundLowerBound_bin, LeftBackgroundUpperBound_bin) + H1D_detectedV0s_DataLike_SmallpTinterval->Integral(RightBackgroundLowerBound_bin, RightBackgroundUpperBound_bin); // symmetric sampling windows: mean-12sigma to mean-6sigma and mean+6sigma to mean+12Sigma
  
  // linear background function and symmetric sampling windows and same total size for background sampling window (12sigma) and for signal window (12sigma) means it's equal
  // can integrate a linear function to check if one wants to
  Long_t BackgroundExtrapolationInSignalRegion = BackgroundCountInSidebands; 

  double SignalCount = TotalCountInSignalRegion - BackgroundExtrapolationInSignalRegion;
  cout << "leftBackgroundCount = " << H1D_detectedV0s_DataLike_SmallpTinterval->Integral(LeftBackgroundLowerBound_bin, LeftBackgroundUpperBound_bin) << ", rightBackgroundCount = "<< H1D_detectedV0s_DataLike_SmallpTinterval->Integral(RightBackgroundLowerBound_bin, RightBackgroundUpperBound_bin) << endl;
  cout << "TotalCountInSignalRegion = " << TotalCountInSignalRegion << endl;
  // Stat uncertainty on signal count:
  // Sigma_signalCount**2 = Sigma_backgroundCount**2 + Sigma_totalCount**2 = sqrt(Sigma_backgroundCount)**2 + sqrt(Sigma_totalCount)**2
  // see https://twiki.cern.ch/twiki/bin/view/ALICE/PWGLFPAGSTRANGENESSSignalExtraction
  double SignalCount_StatError = sqrt(TotalCountInSignalRegion + BackgroundCountInSidebands);

  std::array<double, 3> Signal_vector = {SignalCount, SignalCount_StatError, (double)TotalCountInSignalRegion};
  return Signal_vector;
}

std::array<double, 3> GetTotalAndBackgroundCount_BackgroundIntegration(TH1D* H1D_detectedV0s_DataLike_SmallpTinterval, double SignalMean, double SignalStandardDev, double SideBandSizeMultiplierModifier, TFitResultPtr FitResult, TF1 *bkgparab) {
  // iwth our larger peaks due to lower pagnetic field, we can't go higher than 8sigma for the lower background bound for lambda (due to cutoff from pi+p masses), or 10sigma for k0s (with current window limit)
  // 8 sigma limit is 2*4, so we sample the signal in [-4sigma;+4sigma] interval
  // for a start, the sidebandModifier should be floats, ints is too large a change
  // could have 2*4, 2*3.5, 2*3? 3 is probably the lowest we can go lest we start missing lot of signal (99.9% of guassian area is within -3;+3 sigma)
  double LeftBackgroundLowerBound = SignalMean - 2*(DefaultSideBandSize_array[ipart] + SideBandSizeMultiplierModifier)*SignalStandardDev;
  double LeftBackgroundUpperBound = SignalMean - (DefaultSideBandSize_array[ipart] + SideBandSizeMultiplierModifier)*SignalStandardDev;
  double RightBackgroundLowerBound = SignalMean + (DefaultSideBandSize_array[ipart] + SideBandSizeMultiplierModifier)*SignalStandardDev;
  double RightBackgroundUpperBound = SignalMean + 2*(DefaultSideBandSize_array[ipart] + SideBandSizeMultiplierModifier)*SignalStandardDev;

  int LeftBackgroundLowerBound_bin = H1D_detectedV0s_DataLike_SmallpTinterval->GetXaxis()->FindBin(LeftBackgroundLowerBound);
  int LeftBackgroundUpperBound_bin = H1D_detectedV0s_DataLike_SmallpTinterval->GetXaxis()->FindBin(LeftBackgroundUpperBound);
  int RightBackgroundLowerBound_bin = H1D_detectedV0s_DataLike_SmallpTinterval->GetXaxis()->FindBin(RightBackgroundLowerBound);
  int RightBackgroundUpperBound_bin = H1D_detectedV0s_DataLike_SmallpTinterval->GetXaxis()->FindBin(RightBackgroundUpperBound);

  Long_t TotalCountInSignalRegion = H1D_detectedV0s_DataLike_SmallpTinterval->Integral(LeftBackgroundUpperBound_bin,RightBackgroundLowerBound_bin); //integral of a hist counts the number of bins, does not care for bin width: not true integral

  // Two cases: 0 of bkgparab inside signal region or outside
  //bkgparab = ax+b
  double a = bkgparab->GetParameter(0);
  double b = bkgparab->GetParameter(1);
  double BkgZero = -b/a;
  int BkgZero_bin = H1D_detectedV0s_DataLike_SmallpTinterval->GetXaxis()->FindBin(BkgZero);

  Long_t nbinInvMass = RightBackgroundLowerBound_bin - LeftBackgroundUpperBound_bin;

  double BackgroundCountInSignalRegion = 0;
  double Sigma_integral_fromRootFunction = 0;
  double BackgroundCountInSignalRegion_SystError = 0;
  double * param = bkgparab->GetParameters();
  TMatrixDSym covMatrix = FitResult->GetCovarianceMatrix();

  if (BkgZero < LeftBackgroundUpperBound || RightBackgroundLowerBound < BkgZero){
    //the 0 of the linear function is outside the Signal region: bkgparab doesn't change sign in that region: no issue, integrate bkgparab over the whole region
    BackgroundCountInSignalRegion = bkgparab->Integral(LeftBackgroundUpperBound,RightBackgroundLowerBound)*1./((RightBackgroundLowerBound-LeftBackgroundUpperBound)/nbinInvMass); //bkgparab is not an histogram but a function; here we need to convert it to a count

    Sigma_integral_fromRootFunction = bkgparab->IntegralError(LeftBackgroundUpperBound,RightBackgroundLowerBound, param, covMatrix.GetMatrixArray());// have to include the scaling factor used to get a pTdifferential curve from simple background count
    BackgroundCountInSignalRegion_SystError = Sigma_integral_fromRootFunction * sqrt(1./((RightBackgroundLowerBound-LeftBackgroundUpperBound)/nbinInvMass));
  }
  else {
    //the 0 of the linear function is inside the Signal region: bkgparab changes sign in that region
    //only add the integral part of bkgparab that is positive
    Long_t nbinInvMass_LeftOfBkgZero = BkgZero_bin - LeftBackgroundUpperBound_bin;
    Long_t nbinInvMass_RightOfBkgZero = RightBackgroundLowerBound_bin - BkgZero_bin;
    double BackgroundCountInSignalRegion_LeftOfBkgZero = bkgparab->Integral(LeftBackgroundUpperBound,BkgZero_bin)*1./((BkgZero-LeftBackgroundUpperBound)/nbinInvMass_LeftOfBkgZero); //bkgparab is not an histogram but a function; here we need to convert it to a count
    double BackgroundCountInSignalRegion_RightOfBkgZero = bkgparab->Integral(BkgZero_bin,RightBackgroundLowerBound)*1./((RightBackgroundLowerBound-BkgZero)/nbinInvMass_RightOfBkgZero); //bkgparab is not an histogram but a function; here we need to convert it to a count

    //errors
    double Sigma_integral_fromRootFunction_LeftOfBkgZero = bkgparab->IntegralError(LeftBackgroundUpperBound,BkgZero_bin, param, covMatrix.GetMatrixArray());// have to include the scaling factor used to get a pTdifferential curve from simple background count
    double BackgroundCountInSignalRegion_SystError_LeftOfBkgZero = Sigma_integral_fromRootFunction_LeftOfBkgZero * sqrt(1./((BkgZero_bin-LeftBackgroundUpperBound)/nbinInvMass_LeftOfBkgZero));
    double Sigma_integral_fromRootFunction_RightOfBkgZero = bkgparab->IntegralError(BkgZero_bin,RightBackgroundLowerBound, param, covMatrix.GetMatrixArray());// have to include the scaling factor used to get a pTdifferential curve from simple background count
    double BackgroundCountInSignalRegion_SystError_RightOfBkgZero = Sigma_integral_fromRootFunction_RightOfBkgZero * sqrt(1./((RightBackgroundLowerBound-BkgZero_bin)/nbinInvMass_RightOfBkgZero));

    if (BackgroundCountInSignalRegion_LeftOfBkgZero > 0) {
      BackgroundCountInSignalRegion = BackgroundCountInSignalRegion_LeftOfBkgZero;
      BackgroundCountInSignalRegion_SystError = BackgroundCountInSignalRegion_SystError_LeftOfBkgZero;
    }
    else{
      if (BackgroundCountInSignalRegion_RightOfBkgZero > 0) {
        BackgroundCountInSignalRegion = BackgroundCountInSignalRegion_RightOfBkgZero;
        BackgroundCountInSignalRegion_SystError = BackgroundCountInSignalRegion_SystError_RightOfBkgZero;
      }
      else{
        cout << "BackgroundCountInSignalRegion_LeftOfBkgZero =0 AND BackgroundCountInSignalRegion_RightOfBkgZero =0; shouldn't be happening; they thus all stay at their initialisation value of 0 " << endl;
      }
    }
  }
  if (BackgroundCountInSignalRegion < 0) {
    // cout << "WARNING - negative background count from fit - (cutArry_ID,ibinPt) = (" << ibinPt << "," << warning_cutArry_ID << ")" << endl;
    cout << "WARNING - negative background count from fit" << endl;
  }
//   double TotalCountInSignalRegion_SystError = sqrt(TotalCountInSignalRegion*TotalCountInSignalRegion); // AIMERIC du futur: pretty sure this should be sqrt(TotalCountInSignalRegion) simply
  double TotalCountInSignalRegion_SystError = sqrt(TotalCountInSignalRegion); // AIMERIC du futur: pretty sure this should be sqrt(TotalCountInSignalRegion) simply

  ///////////// Error calculation for Background integration /////////////
  // (https://root.cern/doc/master/ErrorIntegral_8C.html)
  // estimated integral  and error analytically

  ///////////// Outputs /////////////
  double SignalCount = TotalCountInSignalRegion - BackgroundCountInSignalRegion;
  double SignalCount_StatError = sqrt(TotalCountInSignalRegion_SystError*TotalCountInSignalRegion_SystError + BackgroundCountInSignalRegion_SystError*BackgroundCountInSignalRegion_SystError); //assuming uncorrelated; worst case anyway; one is count one is integration of fit

  std::array<double, 3> Signal_vector = {SignalCount, SignalCount_StatError, (double)TotalCountInSignalRegion};

  ////// Analytical error by hand: gives about 1/2 of the error that IntegralError() gives
  // // double Integral  = a/2*B**2 + b*B - a/2*A**2 - b*A; for our line ax+b between A and B
  // double dI_da = 1./2*(RightBackgroundLowerBound*RightBackgroundLowerBound - LeftBackgroundUpperBound*LeftBackgroundUpperBound); // partial derivative = B**2 - A**2 for our line
  // double dI_db = RightBackgroundLowerBound - LeftBackgroundUpperBound; // partial derivative = B - A for our line
  // // estimated error with correlations
  // double Sigma_Integral = std::sqrt(dI_da*dI_da * covMatrix(0,0) + dI_db*dI_db * covMatrix(1,1) + 2.* dI_da*dI_db * covMatrix(0,1));// have to include the scaling factor used to get a pTdifferential curve from simple background count
  
  // if ( std::fabs(Sigma_integral_fromRootFunction - Sigma_Integral) > 1.E-6*Sigma_Integral ) {
  //   std::cout << " ERROR: test failed : different analytical  integral - (RootAutocalc, PersoDerivation) = (" << Sigma_integral_fromRootFunction << ", " << Sigma_Integral << ") - CHECK demo file again" << std::endl;
  // }

  return Signal_vector;
}



void Get_RawYield_vsPt(TH1D* &hRawYield_vsPt, TH3D* H3D_detectedV0s_DataLike, TH1I* H1I_SelectedEventCount, TFile *file_O2Analysis, double* SignalMean, double* SignalStandardDev, int ipart, double* pTbins, int nbinpT, int ibinXaxisCut_low,  int ibinXaxisCut_high, int warning_cutArry_ID, double SideBandSizeMultiplierModifier, int SignalExtractionType, int InvMassRebinFactor) {
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
  // double parGaussParab[nbinpT][6]; //parab backround
  double parGaussParab[nbinpT][5]; //linear background

  TH1D *H1D_detectedV0s_DataLike_SmallpTinterval[nbinpT];
  TH1D *hPtDiff_RawYield_temp = new TH1D("hPtDiff_RawYield_temp","hPtDiff_RawYield_temp",nbinpT,pTbins);

  float SelectedEventCount = H1I_SelectedEventCount->GetEntries();
  // float SelectedEventCount = H1I_SelectedEventCount->GetBinContent(3);
  cout << "SelectedEventCount = " << SelectedEventCount << endl;
  cout << "*--------------------- Get_RawYield_vsPt ---------------------*" << endl;

  TCanvas *canvasMassvsPtFits = new TCanvas ("canvasMassPt_fits", NamehistoInvMass[ipart], 800, 600);
  // TCanvas *canvasMassvsPt = new TCanvas ("canvasMassPt", NamehistoInvMass[ipart], 1600, 1200);
  if (ipart == 0) {
  canvasMassvsPtFits->Divide(std::ceil(std::sqrt(nbinpT)),std::ceil(std::sqrt(nbinpT))-1);
  }
  else {
  canvasMassvsPtFits->Divide(std::ceil(std::sqrt(nbinpT)),std::ceil(std::sqrt(nbinpT)));
  }
  TH1 *hFrameFits[nbinpT];
  cout << "testAAA" << endl;
  for(int ibinPt = 1; ibinPt <= nbinpT; ibinPt++){
    int iPt = ibinPt-1; //different convention for bins that start at 1 and arrays that start at 0
    
    canvasMassvsPtFits->cd(ibinPt); //ATTENTION: +1 because canvas divide starts counting from 1

    int ibinPt_originalAxis_low = H3D_detectedV0s_DataLike->GetYaxis()->FindBin(pTbins[iPt]);
    int ibinPt_originalAxis_up = H3D_detectedV0s_DataLike->GetYaxis()->FindBin(pTbins[iPt+1]) - 1;
    // cout << "ibinPt: " << ibinPt << "; originalAxis [low;up] = [" << ibinPt_originalAxis_low << ";" << ibinPt_originalAxis_up << "]" << endl;

    H1D_detectedV0s_DataLike_SmallpTinterval[iPt] = (TH1D*)H3D_detectedV0s_DataLike_rebinnedZ->ProjectionZ("Get_RawYield_vsPt_InvMass_"+NamePart[ipart]+Form("_Pt%i", ibinPt),ibinXaxisCut_low,ibinXaxisCut_high,ibinPt_originalAxis_low,ibinPt_originalAxis_up);
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

    H1D_detectedV0s_DataLike_SmallpTinterval[iPt]->Fit(gauss[iPt], "R0QL");

    float lowEdge_fit = gauss[iPt]->GetParameter(1) - 2* DefaultSideBandSize_array[ipart] * gauss[iPt]->GetParameter(2);
    float upEdge_fit = gauss[iPt]->GetParameter(1) + 2* DefaultSideBandSize_array[ipart] * gauss[iPt]->GetParameter(2);

    // Linear background; less issues with brackground fit going negative; better for pp
    // bkgparab[iPt] = new TF1("line_"+NamePart[ipart]+Form("_Pt%i", ibinPt), fline, liminf[ipart], limsup[ipart], 3);
    bkgparab[iPt] = new TF1("line_"+NamePart[ipart]+Form("_Pt%i", ibinPt), fline, lowEdge_fit, upEdge_fit, 3);
    bkgparab[iPt]->SetLineColor(kGreen);
    bkgparab[iPt]->FixParameter(2, ipart);

    // GaussPlusPolynom[iPt] = new TF1("GaussPlusPolynom_"+NamePart[ipart]+Form("_Pt%i", ibinPt), "gaus(0)+pol1(3)", liminf[ipart], limsup[ipart]);
    GaussPlusPolynom[iPt] = new TF1("GaussPlusPolynom_"+NamePart[ipart]+Form("_Pt%i", ibinPt), "gaus(0)+pol1(3)", lowEdge_fit, upEdge_fit);
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
    int icolor=0;
    H1D_detectedV0s_DataLike_SmallpTinterval[iPt]->SetMarkerStyle(markers[0]);
    H1D_detectedV0s_DataLike_SmallpTinterval[iPt]->SetMarkerColor(colors [icolor]);
    H1D_detectedV0s_DataLike_SmallpTinterval[iPt]->SetLineColor  (colors [icolor]);
    H1D_detectedV0s_DataLike_SmallpTinterval[iPt]->SetMarkerSize (0.5);
    // bkgparab[iPt]->Draw("same");
    // gauss[iPt]->Draw("same");


    // TLatex * text_pTbin = new TLatex (0.62,0.75,Form("pT [%.1f;%.1f] GeV", pTbins[iPt], pTbins[iPt+1]));
    // text_pTbin->SetNDC(kTRUE);
    // text_pTbin->Draw();
    // TLatex * textContext = new TLatex (0.18,0.82,"#bf{ALICE Performance}"); //BOLD
    // textContext->SetTextSize(0.05);
    // textContext->SetNDC(kTRUE); //remove if I want x,y in TLatex to be in the coordinate system of the histogram
    // textContext->Draw();
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

    *SignalMean = parGaussParab[iPt][1];
    *SignalStandardDev = parGaussParab[iPt][2];

    std::array<double, 3> Signal_vector;

    if (SignalExtractionType == 0) {
      Signal_vector = GetTotalAndBackgroundCount_SidebandBackgroundExtrapolation(H1D_detectedV0s_DataLike_SmallpTinterval[iPt], *SignalMean, *SignalStandardDev, SideBandSizeMultiplierModifier);
    }
    if (SignalExtractionType == 1) {
      Signal_vector = GetTotalAndBackgroundCount_BackgroundIntegration(H1D_detectedV0s_DataLike_SmallpTinterval[iPt], *SignalMean, *SignalStandardDev, SideBandSizeMultiplierModifier, fFitResult[iPt], bkgparab[iPt]);
    }

    double SignalCount = Signal_vector[0];
    double SignalCount_StatError = Signal_vector[1];
    if (SignalCount < 0) {
      cout << "WARNING - negative SignalCount count from fit - (cutArry_ID,ibinPt) = (" << ibinPt << "," << warning_cutArry_ID << ")" << endl;
    }
    // cout << "totalCount= " << totalCount << ", backgroundCount= " << backgroundCount << endl;
    // cout << "SignalCount= " << SignalCount << endl;

    double dpT = hPtDiff_RawYield_temp->GetXaxis()->GetBinWidth(ibinPt);
    double drapidity = 1.5; //1.5
    double d2N_dpTdy = SignalCount *1./dpT*1./drapidity;
    hPtDiff_RawYield_temp->SetBinContent(ibinPt,d2N_dpTdy);
    hPtDiff_RawYield_temp->SetBinError(ibinPt,SignalCount_StatError*1./dpT*1./drapidity);  // error on d2N_dpTdy
    
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

void Get_RawYield_vsPt_FeeddownCorrected(TH1D* &hRawYield_vsPt, TH3D* H3D_detectedV0s_DataLike, TH1I* H1I_SelectedEventCount, TFile *file_O2Analysis, double* SignalMean, double* SignalStandardDev, int ipart, double* pTbins, int nbinpT, int ibinXaxisCut_low,  int ibinXaxisCut_high, int warning_cutArry_ID, double SideBandSizeMultiplierModifier, int SignalExtractionType, int InvMassRebinFactor) {
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
  // double parGaussParab[nbinpT][6]; //parab backround
  double parGaussParab[nbinpT][5]; //linear background

  TH1D *H1D_detectedV0s_DataLike_SmallpTinterval[nbinpT];
  TH1D *hPtDiff_RawYield_temp = new TH1D("hPtDiff_RawYield_temp","hPtDiff_RawYield_temp",nbinpT,pTbins);

  float SelectedEventCount = H1I_SelectedEventCount->GetEntries();
  cout << "*--------------------- Get_RawYield_vsPt ---------------------*" << endl;

  TCanvas *canvasMassvsPtFits = new TCanvas ("canvasMassPt_fits", NamehistoInvMass[ipart], 800, 600);
  // TCanvas *canvasMassvsPt = new TCanvas ("canvasMassPt", NamehistoInvMass[ipart], 1600, 1200);
  if (ipart == 0) {
  canvasMassvsPtFits->Divide(std::ceil(std::sqrt(nbinpT)),std::ceil(std::sqrt(nbinpT))-1);
  }
  else {
  canvasMassvsPtFits->Divide(std::ceil(std::sqrt(nbinpT)),std::ceil(std::sqrt(nbinpT)));
  }
  TH1 *hFrameFits[nbinpT];
  for(int ibinPt = 1; ibinPt <= nbinpT; ibinPt++){
    int iPt = ibinPt-1; //different convention for bins that start at 1 and arrays that start at 0
    
    canvasMassvsPtFits->cd(ibinPt); //ATTENTION: +1 because canvas divide starts counting from 1

    int ibinPt_originalAxis_low = H3D_detectedV0s_DataLike->GetYaxis()->FindBin(pTbins[iPt]);
    int ibinPt_originalAxis_up = H3D_detectedV0s_DataLike->GetYaxis()->FindBin(pTbins[iPt+1]) - 1;
    // cout << "ibinPt: " << ibinPt << "; originalAxis [low;up] = [" << ibinPt_originalAxis_low << ";" << ibinPt_originalAxis_up << "]" << endl;

    H1D_detectedV0s_DataLike_SmallpTinterval[iPt] = (TH1D*)H3D_detectedV0s_DataLike_rebinnedZ->ProjectionZ("Get_RawYield_vsPt_InvMass_"+NamePart[ipart]+Form("_Pt%i", ibinPt),ibinXaxisCut_low,ibinXaxisCut_high,ibinPt_originalAxis_low,ibinPt_originalAxis_up);
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
    H1D_detectedV0s_DataLike_SmallpTinterval[iPt]->Fit(gauss[iPt], "R0QL");

    float lowEdge_fit = gauss[iPt]->GetParameter(1) - 2* DefaultSideBandSize_array[ipart] * gauss[iPt]->GetParameter(2);
    float upEdge_fit = gauss[iPt]->GetParameter(1) + 2* DefaultSideBandSize_array[ipart] * gauss[iPt]->GetParameter(2);
    // cout << "lowEdge_fit = " << lowEdge_fit << ", upEdge_fit = " << upEdge_fit << endl;

    // Linear background; less issues with brackground fit going negative; better for pp
    // bkgparab[iPt] = new TF1("line_"+NamePart[ipart]+Form("_Pt%i", ibinPt), fline, liminf[ipart], limsup[ipart], 3);
    bkgparab[iPt] = new TF1("line_"+NamePart[ipart]+Form("_Pt%i", ibinPt), fline, lowEdge_fit, upEdge_fit, 3);
    bkgparab[iPt]->SetLineColor(kGreen);
    bkgparab[iPt]->FixParameter(2, ipart);

    // GaussPlusPolynom[iPt] = new TF1("GaussPlusPolynom_"+NamePart[ipart]+Form("_Pt%i", ibinPt), "gaus(0)+pol1(3)", liminf[ipart], limsup[ipart]);
    GaussPlusPolynom[iPt] = new TF1("GaussPlusPolynom_"+NamePart[ipart]+Form("_Pt%i", ibinPt), "gaus(0)+pol1(3)", lowEdge_fit, upEdge_fit);
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
    GaussPlusPolynom[iPt]->SetParLimits(2, 0.001, 0.05);
    // cout << "par0 upper limit" <<  1.1*H1D_detectedV0s_DataLike_SmallpTinterval[iPt]->GetBinContent(H1D_detectedV0s_DataLike_SmallpTinterval[iPt]->GetMaximumBin()) << endl;

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
    int icolor=0;
    H1D_detectedV0s_DataLike_SmallpTinterval[iPt]->SetMarkerStyle(markers[0]);
    H1D_detectedV0s_DataLike_SmallpTinterval[iPt]->SetMarkerColor(colors [icolor]);
    H1D_detectedV0s_DataLike_SmallpTinterval[iPt]->SetLineColor  (colors [icolor]);
    H1D_detectedV0s_DataLike_SmallpTinterval[iPt]->SetMarkerSize (0.5);
    // bkgparab[iPt]->Draw("Bsame");
    // gauss[iPt]->Draw("same");


    // TLatex * text_pTbin = new TLatex (0.62,0.75,Form("pT [%.1f;%.1f] GeV", pTbins[iPt], pTbins[iPt+1]));
    // text_pTbin->SetNDC(kTRUE);
    // text_pTbin->Draw();
    // TLatex * textContext = new TLatex (0.18,0.82,"#bf{ALICE Performance}"); //BOLD
    // textContext->SetTextSize(0.05);
    // textContext->SetNDC(kTRUE); //remove if I want x,y in TLatex to be in the coordinate system of the histogram
    // textContext->Draw();
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

    *SignalMean = parGaussParab[iPt][1];
    *SignalStandardDev = parGaussParab[iPt][2];

    std::array<double, 3> Signal_vector;

    if (SignalExtractionType == 0) {
      Signal_vector = GetTotalAndBackgroundCount_SidebandBackgroundExtrapolation(H1D_detectedV0s_DataLike_SmallpTinterval[iPt], *SignalMean, *SignalStandardDev, SideBandSizeMultiplierModifier);
    }
    if (SignalExtractionType == 1) {
      Signal_vector = GetTotalAndBackgroundCount_BackgroundIntegration(H1D_detectedV0s_DataLike_SmallpTinterval[iPt], *SignalMean, *SignalStandardDev, SideBandSizeMultiplierModifier, fFitResult[iPt], bkgparab[iPt]);
    }

    double SignalCount = Signal_vector[0];
    double SignalCount_StatError = Signal_vector[1];
    if (SignalCount < 0) {
      cout << "WARNING - negative SignalCount count from fit - (cutArry_ID,ibinPt) = (" << ibinPt << "," << warning_cutArry_ID << ")" << endl;
    }
    // cout << "totalCount= " << totalCount << ", backgroundCount= " << backgroundCount << endl;
    // cout << "SignalCount= " << SignalCount << endl;

    double dpT = hPtDiff_RawYield_temp->GetXaxis()->GetBinWidth(ibinPt);
    double drapidity = 1.5; //1.5
    double d2N_dpTdy = SignalCount *1./dpT*1./drapidity;
    hPtDiff_RawYield_temp->SetBinContent(ibinPt,d2N_dpTdy);
    hPtDiff_RawYield_temp->SetBinError(ibinPt,SignalCount_StatError*1./dpT*1./drapidity);  // error on d2N_dpTdy 
    
    cout << "(ibinPt, SignalCount, TotalCount) = (" << ibinPt << ", " << SignalCount << ", " << Signal_vector[2] << ")" << endl;
    cout << "         SignalCount_StatError = " << SignalCount_StatError << endl;
  }

  canvasMassvsPtFits->SaveAs("Fits_invMass_ptDiff_"+NamePart[ipart]+".pdf");

  ////////////////////////////////// Fit End //////////////////////////////////
  // hSignalBackgroundRatio_vsPt = (TH1D*)hSignalCountFit_vsPt->Clone("hSignalBackgroundRatio_vsPt");
  // hSignalBackgroundRatio_vsPt->Divide(hSignalCountFit_vsPt,hBackgroundCountFit_vsPt);
  hRawYield_vsPt = (TH1D*)hPtDiff_RawYield_temp->Clone("hRawYield_vsPt");
  hRawYield_vsPt->Scale(1./SelectedEventCount);

  // Lambda and AntiLambda Feeddown
  TH1D *H1D_Feeddown_Correction = new TH1D("H1D_Feeddown_Correction", "H1D_Feeddown_Correction", nbinpT, pTbins);
  if (ipart == 1 || ipart == 2) {
    Get_FeedDown(H1D_Feeddown_Correction, ipart, file_O2Analysis, pTbins, nbinpT);
    cout << "DEBUG Get_RawYield_vsPt_FeeddownCorrected 1" <<endl;
    H1D_Feeddown_Correction->Scale(-1.);
    cout << "DEBUG Get_RawYield_vsPt_FeeddownCorrected 2" <<endl;
    hRawYield_vsPt->Add(H1D_Feeddown_Correction);
    cout << "DEBUG Get_RawYield_vsPt_FeeddownCorrected 3" <<endl;
    // cout << "test00000" << endl;
    // for(int ibinPt = 1; ibinPt <= nbinpT; ibinPt++){
    //   cout << "H1D_Feeddown_Correction(" << ibinPt << ") = " << H1D_Feeddown_Correction->GetBinContent(ibinPt) << endl;
    // }
  }

  delete H1D_Feeddown_Correction;
    cout << "DEBUG Get_RawYield_vsPt_FeeddownCorrected 4" <<endl;

  delete hPtDiff_RawYield_temp;
    cout << "DEBUG Get_RawYield_vsPt_FeeddownCorrected 5" <<endl;
}

void Get_RawYield_vsPt_FeeddownCorrected_withUnfolding(TH1D* &hRawYield_vsPt, TH3D* H3D_detectedV0s_DataLike, TH1I* H1I_SelectedEventCount, TFile *file_O2Analysis, double* SignalMean, double* SignalStandardDev, int ipart, double* pTbins, int nbinpT, int ibinXaxisCut_low,  int ibinXaxisCut_high, int warning_cutArry_ID, double SideBandSizeMultiplierModifier, int SignalExtractionType, int InvMassRebinFactor) {
  float SelectedEventCount = H1I_SelectedEventCount->GetEntries();

  TH1D* truth = (TH1D*)file_O2Analysis->Get("lambdakzero-particle-count-mc/h"+NamePart[ipart]+"Count_PtDiff");
  TH1D* truth_rebinned = (TH1D*)truth->Rebin(nbinpT,"truth_rebinned",pTbins);
  
  Get_RawYield_vsPt_FeeddownCorrected_RawCount(hRawYield_vsPt, H3D_detectedV0s_DataLike, H1I_SelectedEventCount, file_O2Analysis, SignalMean, SignalStandardDev, ipart, pTbins, nbinpT, ibinXaxisCut_low, ibinXaxisCut_high, warning_cutArry_ID, SideBandSizeMultiplierModifier, SignalExtractionType, InvMassRebinFactor);
  TH2D *H2D_PtResponseMatrix = new TH2D("H2D_PtResponseMatrix", "H2D_PtResponseMatrix", nbinpT, pTbins, nbinpT, pTbins);
  Get_PtResponseMatrix(H2D_PtResponseMatrix, file_O2Analysis, pTbins, nbinpT);

  // Get_RawYield_vsPt_FeeddownCorrected(hRawYield_vsPt, H3D_detectedV0s_DataLike, H1I_SelectedEventCount, file_O2Analysis, SignalMean, SignalStandardDev, ipart, pTbins, nbinpT, ibinXaxisCut_low, ibinXaxisCut_high, warning_cutArry_ID, SideBandSizeMultiplierModifier, SignalExtractionType, InvMassRebinFactor);
  // TH2D *H2D_PtResponseMatrix_Density = new TH2D("H2D_PtResponseMatrix", "H2D_PtResponseMatrix", nbinpT, pTbins, nbinpT, pTbins);
  // Get_PtResponseMatrix_Density(H2D_PtResponseMatrix, file_O2Analysis, pTbins, nbinpT);
  // for(int ibinPt = 1; ibinPt <= nbinpT; ibinPt++){
  //   double dpT = truth_rebinned->GetXaxis()->GetBinWidth(ibinPt);
  //   double drapidity = 1.5; //1.5
  //   double d2N_dpTdy = truth_rebinned->GetBinContent(ibinPt) *1./dpT*1./drapidity;
  //   truth_rebinned->SetBinContent(ibinPt,d2N_dpTdy);
  //   truth_rebinned->SetBinError(ibinPt,truth_rebinned->GetBinError(ibinPt) *1./dpT*1./drapidity);  // error on d2N_dpTdy 
  // }    
  // truth_rebinned->Scale(1./SelectedEventCount);

  TH1D* measured = (TH1D*)hRawYield_vsPt->Clone("measured");

  // TH1D* projX = (TH1D*)H2D_PtResponseMatrix->ProjectionX("projX",0,-1);
  // TH1D* projY = (TH1D*)H2D_PtResponseMatrix->ProjectionY("projY",0,-1);
  cout << "--------------------------------------------------------------------- PRE UNFOLDING TEST" << endl;
  for(int ibinPt = 1; ibinPt <= nbinpT; ibinPt++){
    cout << "measured(" << ibinPt << ") = " << measured->GetBinContent(ibinPt) << endl;
    cout << "truth_rebinned(" << ibinPt << ") = " << truth_rebinned->GetBinContent(ibinPt) << endl;
    for(int jbinPt = 1; jbinPt <= nbinpT; jbinPt++){
      cout << "H2D_PtResponseMatrix(" << ibinPt << "," << jbinPt << ") = " << H2D_PtResponseMatrix->GetBinContent(ibinPt,jbinPt) << endl;
    }
  }

  RooUnfoldResponse* response = new RooUnfoldResponse(measured, truth_rebinned, H2D_PtResponseMatrix);
  // RooUnfoldResponse* response = new RooUnfoldResponse(projX, projY, H2D_PtResponseMatrix, "", "");
  // RooUnfold* unfold = new RooUnfoldBayes(response, measured, 4);
  RooUnfold* unfold = new RooUnfoldBinByBin(response, measured);

  // TH1D* hist_unfold = (TH1D*)unfold.Hunfold();
  // TH1D* hist_unfold = (TH1D*)unfold->Hreco()->Clone("hist_unfold");
  TH1D* hist_unfold = static_cast<TH1D*>(unfold->Hreco());


  hRawYield_vsPt = (TH1D*)hist_unfold->Clone("hRawYield_vsPt");

  for(int ibinPt = 1; ibinPt <= nbinpT; ibinPt++){
    double dpT = hRawYield_vsPt->GetXaxis()->GetBinWidth(ibinPt);
    double drapidity = 1.5; //1.5
    double d2N_dpTdy = hRawYield_vsPt->GetBinContent(ibinPt) *1./dpT*1./drapidity;
    hRawYield_vsPt->SetBinContent(ibinPt,d2N_dpTdy);
    hRawYield_vsPt->SetBinError(ibinPt,hRawYield_vsPt->GetBinError(ibinPt) *1./dpT*1./drapidity);  // error on d2N_dpTdy 
  }    

  hRawYield_vsPt->Scale(1./SelectedEventCount);

  // Lambda and AntiLambda Feeddown
  TH1D *H1D_Feeddown_Correction = new TH1D("H1D_Feeddown_Correction", "H1D_Feeddown_Correction", nbinpT, pTbins);
  if (ipart == 1 || ipart == 2) {
    Get_FeedDown(H1D_Feeddown_Correction, ipart, file_O2Analysis, pTbins, nbinpT);
    H1D_Feeddown_Correction->Scale(-1.);
    hRawYield_vsPt->Add(H1D_Feeddown_Correction);
    // cout << "test00000" << endl;
    for(int ibinPt = 1; ibinPt <= nbinpT; ibinPt++){
    cout << "H1D_Feeddown_Correction(" << ibinPt << ") = " << H1D_Feeddown_Correction->GetBinContent(ibinPt) << endl;
    }
  }
  delete H1D_Feeddown_Correction;


  cout << "--------------------------------------------------------------------- POST UNFOLDING TEST" << hist_unfold->GetEntries() << endl;
  for(int ibinPt = 1; ibinPt <= nbinpT; ibinPt++){
    cout << "truth_rebinned(" << ibinPt << ") = " << truth_rebinned->GetBinContent(ibinPt) << endl;
    cout << "H2D_PtResponseMatrix(" << ibinPt << ") = " << H2D_PtResponseMatrix->GetBinContent(ibinPt,ibinPt) << endl;
    cout << "hist_unfold(" << ibinPt << ") = " << hist_unfold->GetBinContent(ibinPt) << endl;
    cout << "measured(" << ibinPt << ") = " << measured->GetBinContent(ibinPt) << endl;
    cout << "measuredPostUnfoldPostCorrections(" << ibinPt << ") = " << hRawYield_vsPt->GetBinContent(ibinPt) << endl;
  }

  // TUnfold unfold(H2D_PtResponseMatrix, TUnfold::kHistMapOutputHoriz, TUnfold::kRegModeNone); // kRegModeNone regularisation method choice, here none, will see what's best later
  // // set input distribution and bias scale (=0)
  // cout << "AMIERIC DEBUG unfolding 4" << endl;
  // if(unfold.SetInput(hRawYield_vsPt, 0.0)>=10000) {
  //   std::cout<<"Unfolding result may be wrong\n";
  // }
  // cout << "AMIERIC DEBUG unfolding 5" << endl;

  // // do the unfolding here
  // double tauMin=0.0;
  // double tauMax=0.0;
  // int nScan=30;
  // int iBest;
  // TSpline *logTauX,*logTauY;
  // TGraph *lCurve;
  // // this method scans the parameter tau
  // // finally, the unfolding is done for the "best" choice of tau
  // iBest=unfold.ScanLcurve(nScan,tauMin,tauMax,&lCurve,&logTauX,&logTauY);
  // std::cout<<"tau="<<unfold.GetTau()<<"\n";
  // std::cout<<"chi**2="<<unfold.GetChi2A()<<"+"<<unfold.GetChi2L()
  //          <<" / "<<unfold.GetNdf()<<"\n";
  // // save point corresponding to the kink in the L curve as TGraph
  // double t[1],x[1],y[1];
  // logTauX->GetKnot(iBest,t[0],x[0]);
  // logTauY->GetKnot(iBest,t[0],y[0]);
  // TGraph *bestLcurve=new TGraph(1,x,y);
  // TGraph *bestLogTauX=new TGraph(1,t,x);

  // //============================================================
  // // extract unfolding results into histograms
  // // set up a bin map, excluding underflow and overflow bins
  // // the binMap relates the the output of the unfolding to the final
  // // histogram bins
  // int *binMap=new int[nbinpT+2];
  // for(int i=1;i<=nbinpT;i++) binMap[i]=i;
  // binMap[0]=-1;
  // binMap[nbinpT+1]=-1;
  // TH1D *histPtUnfold=new TH1D("Unfolded",";pT(gen)",nbinpT,pTbins);
  // unfold.GetOutput(histPtUnfold,binMap);
  // TH1D *histPtDetFold=new TH1D("FoldedBack","pT(det)",nbinpT,pTbins);
  // unfold.GetFoldedOutput(histPtDetFold); // Aimeric: probably for comparison? some kind of closure test
  // // store global correlation coefficients
  // TH1D *histRhoi=new TH1D("rho_I","mass",nbinpT,pTbins);
  // unfold.GetRhoI(histRhoi,binMap);
  // delete[] binMap;
  // binMap=0;
}

void Get_RawYield_vsPt_FeeddownCorrected_RawCount(TH1D* &hRawYield_vsPt, TH3D* H3D_detectedV0s_DataLike, TH1I* H1I_SelectedEventCount, TFile *file_O2Analysis, double* SignalMean, double* SignalStandardDev, int ipart, double* pTbins, int nbinpT, int ibinXaxisCut_low,  int ibinXaxisCut_high, int warning_cutArry_ID, double SideBandSizeMultiplierModifier, int SignalExtractionType, int InvMassRebinFactor) {
  // same as Get_RawYield_vsPt_FeeddownCorrected but does not divide by Nevent or delta_pT or delta_y; better for unfolding

  // what is fed as input: TH3D* H3D_detectedV0s_DataLike = (TH3D*)file_O2Analysis->Get("lambdakzero-analysis-mc/h3dMass"+NamePart[ipart]);
  TH3D* H3D_detectedV0s_DataLike_rebinnedZ = (TH3D*)H3D_detectedV0s_DataLike->RebinZ(InvMassRebinFactor,"H3D_detectedV0s_DataLike_rebinnedZ");
  cout << "---------------------------------------- InvMassRebinFactor = " << InvMassRebinFactor << endl;
  ////////////////////////////////// Fit initialisation //////////////////////////////////
  //Fit tools initialisation
  TF1 *gauss[nbinpT];
  TF1 *GaussPlusPolynom[nbinpT];
  TF1 *bkgparab[nbinpT];
  TFitResultPtr fFitResult[nbinpT];
  // double parGaussParab[nbinpT][6]; //parab backround
  double parGaussParab[nbinpT][5]; //linear background

  TH1D *H1D_detectedV0s_DataLike_SmallpTinterval[nbinpT];
  TH1D *hPtDiff_RawYield_temp = new TH1D("hPtDiff_RawYield_temp","hPtDiff_RawYield_temp",nbinpT,pTbins);

  float SelectedEventCount = H1I_SelectedEventCount->GetEntries();
  cout << "*--------------------- Get_RawYield_vsPt ---------------------*" << endl;

  TCanvas *canvasMassvsPtFits = new TCanvas ("canvasMassPt_fits", NamehistoInvMass[ipart], 800, 600);
  // TCanvas *canvasMassvsPt = new TCanvas ("canvasMassPt", NamehistoInvMass[ipart], 1600, 1200);
  if (ipart == 0) {
  canvasMassvsPtFits->Divide(std::ceil(std::sqrt(nbinpT)),std::ceil(std::sqrt(nbinpT))-1);
  }
  else {
  canvasMassvsPtFits->Divide(std::ceil(std::sqrt(nbinpT)),std::ceil(std::sqrt(nbinpT)));
  }
  TH1 *hFrameFits[nbinpT];
  for(int ibinPt = 1; ibinPt <= nbinpT; ibinPt++){
    int iPt = ibinPt-1; //different convention for bins that start at 1 and arrays that start at 0
    
    canvasMassvsPtFits->cd(ibinPt); //ATTENTION: +1 because canvas divide starts counting from 1

    int ibinPt_originalAxis_low = H3D_detectedV0s_DataLike->GetYaxis()->FindBin(pTbins[iPt]);
    int ibinPt_originalAxis_up = H3D_detectedV0s_DataLike->GetYaxis()->FindBin(pTbins[iPt+1]) - 1;
    // cout << "ibinPt: " << ibinPt << "; originalAxis [low;up] = [" << ibinPt_originalAxis_low << ";" << ibinPt_originalAxis_up << "]" << endl;

    H1D_detectedV0s_DataLike_SmallpTinterval[iPt] = (TH1D*)H3D_detectedV0s_DataLike_rebinnedZ->ProjectionZ("Get_RawYield_vsPt_InvMass_"+NamePart[ipart]+Form("_Pt%i", ibinPt),ibinXaxisCut_low,ibinXaxisCut_high,ibinPt_originalAxis_low,ibinPt_originalAxis_up);
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
    H1D_detectedV0s_DataLike_SmallpTinterval[iPt]->Fit(gauss[iPt], "R0QL");

    float lowEdge_fit = gauss[iPt]->GetParameter(1) - 2* DefaultSideBandSize_array[ipart] * gauss[iPt]->GetParameter(2);
    float upEdge_fit = gauss[iPt]->GetParameter(1) + 2* DefaultSideBandSize_array[ipart] * gauss[iPt]->GetParameter(2);
    // cout << "lowEdge_fit = " << lowEdge_fit << ", upEdge_fit = " << upEdge_fit << endl;

    // Linear background; less issues with brackground fit going negative; better for pp
    // bkgparab[iPt] = new TF1("line_"+NamePart[ipart]+Form("_Pt%i", ibinPt), fline, liminf[ipart], limsup[ipart], 3);
    bkgparab[iPt] = new TF1("line_"+NamePart[ipart]+Form("_Pt%i", ibinPt), fline, lowEdge_fit, upEdge_fit, 3);
    bkgparab[iPt]->SetLineColor(kGreen);
    bkgparab[iPt]->FixParameter(2, ipart);

    // GaussPlusPolynom[iPt] = new TF1("GaussPlusPolynom_"+NamePart[ipart]+Form("_Pt%i", ibinPt), "gaus(0)+pol1(3)", liminf[ipart], limsup[ipart]);
    GaussPlusPolynom[iPt] = new TF1("GaussPlusPolynom_"+NamePart[ipart]+Form("_Pt%i", ibinPt), "gaus(0)+pol1(3)", lowEdge_fit, upEdge_fit);
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
    GaussPlusPolynom[iPt]->SetParLimits(2, 0.001, 0.05);
    // cout << "par0 upper limit" <<  1.1*H1D_detectedV0s_DataLike_SmallpTinterval[iPt]->GetBinContent(H1D_detectedV0s_DataLike_SmallpTinterval[iPt]->GetMaximumBin()) << endl;

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
    int icolor=0;
    H1D_detectedV0s_DataLike_SmallpTinterval[iPt]->SetMarkerStyle(markers[0]);
    H1D_detectedV0s_DataLike_SmallpTinterval[iPt]->SetMarkerColor(colors [icolor]);
    H1D_detectedV0s_DataLike_SmallpTinterval[iPt]->SetLineColor  (colors [icolor]);
    H1D_detectedV0s_DataLike_SmallpTinterval[iPt]->SetMarkerSize (0.5);
    // bkgparab[iPt]->Draw("Bsame");
    // gauss[iPt]->Draw("same");


    // TLatex * text_pTbin = new TLatex (0.62,0.75,Form("pT [%.1f;%.1f] GeV", pTbins[iPt], pTbins[iPt+1]));
    // text_pTbin->SetNDC(kTRUE);
    // text_pTbin->Draw();
    // TLatex * textContext = new TLatex (0.18,0.82,"#bf{ALICE Performance}"); //BOLD
    // textContext->SetTextSize(0.05);
    // textContext->SetNDC(kTRUE); //remove if I want x,y in TLatex to be in the coordinate system of the histogram
    // textContext->Draw();
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

    *SignalMean = parGaussParab[iPt][1];
    *SignalStandardDev = parGaussParab[iPt][2];

    std::array<double, 3> Signal_vector;

    if (SignalExtractionType == 0) {
      Signal_vector = GetTotalAndBackgroundCount_SidebandBackgroundExtrapolation(H1D_detectedV0s_DataLike_SmallpTinterval[iPt], *SignalMean, *SignalStandardDev, SideBandSizeMultiplierModifier);
    }
    if (SignalExtractionType == 1) {
      Signal_vector = GetTotalAndBackgroundCount_BackgroundIntegration(H1D_detectedV0s_DataLike_SmallpTinterval[iPt], *SignalMean, *SignalStandardDev, SideBandSizeMultiplierModifier, fFitResult[iPt], bkgparab[iPt]);
    }

    double SignalCount = Signal_vector[0];
    double SignalCount_StatError = Signal_vector[1];
    if (SignalCount < 0) {
      cout << "WARNING - negative SignalCount count from fit - (cutArry_ID,ibinPt) = (" << ibinPt << "," << warning_cutArry_ID << ")" << endl;
    }
    // cout << "totalCount= " << totalCount << ", backgroundCount= " << backgroundCount << endl;
    // cout << "SignalCount= " << SignalCount << endl;

    // double dpT = hPtDiff_RawYield_temp->GetXaxis()->GetBinWidth(ibinPt);
    // double drapidity = 1.5; //1.5
    // double d2N_dpTdy = SignalCount *1./dpT*1./drapidity;
    // hPtDiff_RawYield_temp->SetBinContent(ibinPt,d2N_dpTdy);
    hPtDiff_RawYield_temp->SetBinContent(ibinPt,SignalCount);
    // hPtDiff_RawYield_temp->SetBinError(ibinPt,SignalCount_StatError*1./dpT*1./drapidity);  // error on d2N_dpTdy
    hPtDiff_RawYield_temp->SetBinError(ibinPt,SignalCount_StatError);  // error on signalCount
    
    cout << "(ibinPt, SignalCount, TotalCount) = (" << ibinPt << ", " << SignalCount << ", " << Signal_vector[2] << ")" << endl;
  }

  canvasMassvsPtFits->SaveAs("Fits_invMass_ptDiff_"+NamePart[ipart]+".pdf");

  ////////////////////////////////// Fit End //////////////////////////////////
  // hSignalBackgroundRatio_vsPt = (TH1D*)hSignalCountFit_vsPt->Clone("hSignalBackgroundRatio_vsPt");
  // hSignalBackgroundRatio_vsPt->Divide(hSignalCountFit_vsPt,hBackgroundCountFit_vsPt);
  hRawYield_vsPt = (TH1D*)hPtDiff_RawYield_temp->Clone("hRawYield_vsPt");
  // hRawYield_vsPt->Scale(1./SelectedEventCount);

  // Lambda and AntiLambda Feeddown
  TH1D *H1D_Feeddown_Correction = new TH1D("H1D_Feeddown_Correction", "H1D_Feeddown_Correction", nbinpT, pTbins);
  if (ipart == 1 || ipart == 2) {
    Get_FeedDown_RawCount(H1D_Feeddown_Correction, ipart, file_O2Analysis, pTbins, nbinpT);
    H1D_Feeddown_Correction->Scale(-1.);
    hRawYield_vsPt->Add(H1D_Feeddown_Correction);
    // cout << "test00000" << endl;
  }
  delete H1D_Feeddown_Correction;

  delete hPtDiff_RawYield_temp;
}


void Efficiency_DcaScan(TH1D* &hstat, TH1D* &hsyst, TH1D*&hsystCorr, int ipart, TFile* file_O2Analysis, TString* SaveAs_Title, TString* &texXtitle, TString* &texYtitle, double* pTbins, int nbinpT, double maxDCAcut) {

  TH3D* H3D_detectedV0s_TruePt_DcaHist = (TH3D*)file_O2Analysis->Get("lambdakzero-analysis-mc/h3dMass"+NamePart[ipart]+"Dca_MC_truePt");
  TH3D* H3D_detectedV0s_DataLike_DcaHist = (TH3D*)file_O2Analysis->Get("lambdakzero-analysis-mc/h3dMass"+NamePart[ipart]+"Dca");
  TH1D* TrueV0PtSpectrum = (TH1D*)file_O2Analysis->Get("lambdakzero-particle-count-mc/h"+NamePart[ipart]+"Count_PtDiff");

  int maxDCAcutBin = H3D_detectedV0s_TruePt_DcaHist->GetXaxis()->FindBin(maxDCAcut);
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
  for(int ibinDcaCut = 1; ibinDcaCut <= maxDCAcutBin; ibinDcaCut++){
    int itableDcaCut = ibinDcaCut-1;
    Get_McEfficiency_vsPt(hMcEfficiency_vsPt[itableDcaCut], H3D_detectedV0s_TruePt_DcaHist, TrueV0PtSpectrum, ipart, pTbins, nbinpT, 0, ibinDcaCut, SideBandSizeMultiplierModifier_array[0], min_range_signal[ipart], max_range_signal[ipart], InvMassRebinFactor_standard[ipart]);
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
  int nbinx = hstat->GetNbinsX();
  
  for(int ibinx = 1; ibinx <= nbinx; ibinx++){
    hsyst->SetBinContent (ibinx, hstat->GetBinContent(ibinx));
    hsyst->SetBinError   (ibinx, hstat->GetBinContent(ibinx)*0.08);
    hsystCorr->SetBinContent (ibinx, hstat->GetBinContent(ibinx));
    hsystCorr->SetBinError   (ibinx, hstat->GetBinContent(ibinx)*0.15);
  }

  *SaveAs_Title += "ToBeDeleted";
  texXtitle = texPtX;
  texYtitle = texptDifferentialYield_HEPratio_pseudoEff;

}


void Efficiency_TpcCrossedRowsScan(TH1D* &hstat, TH1D* &hsyst, TH1D*&hsystCorr, int ipart, TFile* file_O2Analysis, TString* SaveAs_Title, TString* &texXtitle, TString* &texYtitle, double* pTbins, int nbinpT, double minTpcCrossedRowscut, double maxTpcCrossedRowscut) {

  TH3D* H3D_detectedV0s_TruePt_TpcCrossedRowsHist = (TH3D*)file_O2Analysis->Get("lambdakzero-analysis-mc/h3dMass"+NamePart[ipart]+"TpcRows_MC_truePt");
  TH3D* H3D_detectedV0s_DataLike_TpcCrossedRowsHist = (TH3D*)file_O2Analysis->Get("lambdakzero-analysis-mc/h3dMass"+NamePart[ipart]+"TpcRows");
  TH1D* TrueV0PtSpectrum = (TH1D*)file_O2Analysis->Get("lambdakzero-particle-count-mc/h"+NamePart[ipart]+"Count_PtDiff");

  int minTpcCrossedRowscutBin = H3D_detectedV0s_TruePt_TpcCrossedRowsHist->GetXaxis()->FindBin(minTpcCrossedRowscut);
  int maxTpcCrossedRowscutBin = H3D_detectedV0s_TruePt_TpcCrossedRowsHist->GetXaxis()->FindBin(maxTpcCrossedRowscut);
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
  for(int ibinTpcCrossedRowsCut = minTpcCrossedRowscutBin; ibinTpcCrossedRowsCut <= maxTpcCrossedRowscutBin; ibinTpcCrossedRowsCut++){
    int itableTpcCrossedRowsCut = ibinTpcCrossedRowsCut-1;
    Get_McEfficiency_vsPt(hMcEfficiency_vsPt[itableTpcCrossedRowsCut], H3D_detectedV0s_TruePt_TpcCrossedRowsHist, TrueV0PtSpectrum, ipart, pTbins, nbinpT, ibinTpcCrossedRowsCut, maxTpcCrossedRowscutBin, SideBandSizeMultiplierModifier_array[0], min_range_signal[ipart], max_range_signal[ipart], InvMassRebinFactor_standard[ipart]);
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
  int nbinx = hstat->GetNbinsX();
  
  for(int ibinx = 1; ibinx <= nbinx; ibinx++){
    hsyst->SetBinContent (ibinx, hstat->GetBinContent(ibinx));
    hsyst->SetBinError   (ibinx, hstat->GetBinContent(ibinx)*0.08);
    hsystCorr->SetBinContent (ibinx, hstat->GetBinContent(ibinx));
    hsystCorr->SetBinError   (ibinx, hstat->GetBinContent(ibinx)*0.15);
  }

  *SaveAs_Title += "ToBeDeleted";
  texXtitle = texPtX;
  texYtitle = texptDifferentialYield_HEPratio_pseudoEff;

}


void Get_Systematics_OneCut(TH1D* hSystematicUncertainty, TH1D* hSystematicUncertainty_PreBarlow, int ipart, TFile **file_O2Analysis_CutVariation_array, int nCutsIterations, double* pTbins, int nbinpT, int InvMassRebinFactor) {
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

    double SignalMean, SignalStandardDev;
    Get_RawYield_vsPt_FeeddownCorrected(hRawYield_vsPt[cutsArrayIterator], H3D_DetectedV0s_MCdatalike_cutVariation[cutsArrayIterator], H1I_SelectedEventCount[cutsArrayIterator], file_O2Analysis_CutVariation_array[cutsArrayIterator], &SignalMean, &SignalStandardDev, ipart, pTbins, nbinpT, 0, -1, cutsArrayIterator, SideBandSizeMultiplierModifier_array[0], SignalExtractionType_default, InvMassRebinFactor); ///real line that need to take the place of the one above once it's been written
    Get_McEfficiency_vsPt(hMcEfficiency_vsPt[cutsArrayIterator], H3D_DetectedV0s_MCtrue_cutVariation[cutsArrayIterator], H1D_TrueV0PtSpectrum_cutVariation[cutsArrayIterator], ipart, pTbins, nbinpT, 0, -1, SideBandSizeMultiplierModifier_array[0], SignalMean, SignalStandardDev, InvMassRebinFactor);

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
    hCorrectedYield_vsPt[cutsArrayIterator]->Divide(hRawYield_vsPt[cutsArrayIterator],hMcEfficiency_vsPt[cutsArrayIterator], 1., 1., "");
  // Efficiency errors are binomial-ish;  TGraphAsymmErrors::BayesDivide (https://root.cern.ch/doc/master/classTH1.html#ac782a09c31b4f7de40f8fd4f77efa090)
  // and https://inspirehep.net/files/57287ac8e45a976ab423f3dd456af694
  // PWGLF strangeness PAG says to just use option "b" for divide, and just see it as a binomial https://twiki.cern.ch/twiki/bin/view/ALICE/PWGLFPAGSTRANGENESSEfficiency
  }

  /////// REFERENCE YIELD ///////
  TH1D* hCorrectedYield_vsPt_REF;
  TH3D* H3D_DetectedV0s_MCtrue_cutVariation_REF = (TH3D*)file_O2Analysis->Get("lambdakzero-analysis-mc/h3dMass"+NamePart[ipart]+"_MC_truePt");
  TH1I* H1I_SelectedEventCount_REF = (TH1I*)file_O2Analysis->Get("lambdakzero-analysis-mc/hSelectedEventCounter");
  TH1D* H1D_TrueV0PtSpectrum_cutVariation_REF = (TH1D*)file_O2Analysis->Get("lambdakzero-particle-count-mc/h"+NamePart[ipart]+"Count_PtDiff");
  TH3D* H3D_DetectedV0s_MCdatalike_cutVariation_REF = (TH3D*)file_O2Analysis->Get("lambdakzero-analysis-mc/h3dMass"+NamePart[ipart]);
  TH1D *hMcEfficiency_vsPt_REF;
  TH1D *hRawYield_vsPt_REF;
  double SignalMean, SignalStandardDev;
  Get_RawYield_vsPt_FeeddownCorrected(hRawYield_vsPt_REF, H3D_DetectedV0s_MCdatalike_cutVariation_REF, H1I_SelectedEventCount_REF, file_O2Analysis, &SignalMean, &SignalStandardDev, ipart, pTbins, nbinpT, 0, -1, 0, SideBandSizeMultiplierModifier_array[0], SignalExtractionType_default, InvMassRebinFactor); ///real line that need to take the place of the one above once it's been written
  Get_McEfficiency_vsPt(hMcEfficiency_vsPt_REF, H3D_DetectedV0s_MCtrue_cutVariation_REF, H1D_TrueV0PtSpectrum_cutVariation_REF, ipart, pTbins, nbinpT, 0, -1, SideBandSizeMultiplierModifier_array[0], SignalMean, SignalStandardDev, InvMassRebinFactor);

  // //debug loop
  // for(int ibinPt = 1; ibinPt <= nbinpT; ibinPt++){
  //   cout << "cutID: " << cutsArrayIterator << ", pT bin: " << ibinPt << endl;
  //   cout << "hRawYield = " << hRawYield_vsPt[cutsArrayIterator]->GetBinContent(ibinPt) << ", hRawYield_error = " << hRawYield_vsPt[cutsArrayIterator]->GetBinError(ibinPt) << endl;
  //   cout << "hEfficiency = " << hMcEfficiency_vsPt[cutsArrayIterator]->GetBinContent(ibinPt) << ", hEfficiency_error = " << hMcEfficiency_vsPt[cutsArrayIterator]->GetBinError(ibinPt) << endl;
  // }
  hCorrectedYield_vsPt_REF = (TH1D*)hRawYield_vsPt_REF->Clone("hCorrectedYield_vsPt_REF");
  hCorrectedYield_vsPt_REF->Reset("M");
  // hCorrectedYield_vsPt[cutsArrayIterator]->Sumw2();
  // hRawYield_vsPt[cutsArrayIterator]->Sumw2();
  hCorrectedYield_vsPt_REF->Divide(hRawYield_vsPt_REF,hMcEfficiency_vsPt_REF, 1., 1., "");
  ///////REF YIELD END /////

  // TH1D* hSystematicUncertainty[cutsIterator];
  double hSigmaBarlow[nbinpT];
  double CorrectedYieldDifference;
  int id_cutsArray_maxDeviation = 2;
  for(int ibinPt = 1; ibinPt <= nbinpT; ibinPt++){
    int PtArrayIterator = ibinPt - 1;
    double MaxYieldDifference_cutsIteratorArray = 0;
    for(int cutsArrayIterator = 0; cutsArrayIterator < nCutsIterations; cutsArrayIterator++){
      CorrectedYieldDifference = abs(hCorrectedYield_vsPt[cutsArrayIterator]->GetBinContent(ibinPt) - hCorrectedYield_vsPt_REF->GetBinContent(ibinPt));
      if (CorrectedYieldDifference > MaxYieldDifference_cutsIteratorArray) {
        // MaxDifference_cutsIteratorArray = CorrectedYieldDifference;
        MaxYieldDifference_cutsIteratorArray = CorrectedYieldDifference;
        id_cutsArray_maxDeviation = cutsArrayIterator;
      }
    }

    double SystUncertainty = MaxYieldDifference_cutsIteratorArray;

    // Barlow condition for cut variation systematics (Systematic Errors: facts and fictions, by Roger Barlow, https://arxiv.org/abs/hep-ex/0207026)
    double StatUncertainty_REF = hCorrectedYield_vsPt_REF->GetBinError(ibinPt);
    double StatUncertainty_MaxDeviationCase = hCorrectedYield_vsPt[id_cutsArray_maxDeviation]->GetBinError(ibinPt);
    
    hSigmaBarlow[PtArrayIterator] = sqrt(abs(StatUncertainty_MaxDeviationCase*StatUncertainty_MaxDeviationCase - StatUncertainty_REF*StatUncertainty_REF)); //stat error of the difference in the case of subsample
    cout << "         StatUncertainty_MaxDeviationCase = " << StatUncertainty_MaxDeviationCase << ", StatUncertainty_REF = " << StatUncertainty_REF << endl;
    cout << "                      hRawYield_vsPt_REF = " << hRawYield_vsPt_REF->GetBinError(ibinPt) << ", hMcEfficiency_vsPt_REF = " << hMcEfficiency_vsPt_REF->GetBinError(ibinPt) << endl;

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


// void Get_Systematics_OneCut_McTrueAndMcDatalikeDifferentCuts(TH1D* hSystematicUncertainty, TH1D* hSystematicUncertainty_PreBarlow, int ipart, TFile **file_O2Analysis_CutVariation_Datalike_array, TFile **file_O2Analysis_CutVariation_True_array, int nCutsIterations, double* pTbins, int nbinpT, int InvMassRebinFactor) {
//   // get systematics from cut variation for the array of AnalysisResults.root files file_O2Analysis_CutVariation_array, all one variation of a single cut (cospa for example, or dcav0dau)
//   // nCutsIterations is the number of variations on the cut being looked at; ie the size of the array file_O2Analysis_CutVariation_array

//   TH3D* H3D_DetectedV0s_O2; // ? doesn't seem useful
  
//   TH3D *H3D_DetectedV0s_MCtrue_cutVariation[nCutsIterations];
//   TH3D *H3D_DetectedV0s_MCdatalike_cutVariation[nCutsIterations];
//   TH1I *H1I_SelectedEventCount[nCutsIterations];
//   TH1D *H1D_TrueV0PtSpectrum_cutVariation[nCutsIterations];

//   TH1D *hMcEfficiency_vsPt[nCutsIterations];
//   TH1D *hRawYield_vsPt[nCutsIterations];

//   TH1D *hCorrectedYield_vsPt[nCutsIterations];

//   for(int cutsArrayIterator = 0; cutsArrayIterator < nCutsIterations; cutsArrayIterator++){
//     H3D_DetectedV0s_MCtrue_cutVariation[cutsArrayIterator] = (TH3D*)file_O2Analysis_CutVariation_True_array[cutsArrayIterator]->Get("lambdakzero-analysis-mc/h3dMass"+NamePart[ipart]+"_MC_truePt");
//     H1I_SelectedEventCount[cutsArrayIterator] = (TH1I*)file_O2Analysis_CutVariation_Datalike_array[cutsArrayIterator]->Get("lambdakzero-analysis-mc/hSelectedEventCounter");
//     H1D_TrueV0PtSpectrum_cutVariation[cutsArrayIterator] = (TH1D*)file_O2Analysis_CutVariation_True_array[cutsArrayIterator]->Get("lambdakzero-particle-count-mc/h"+NamePart[ipart]+"Count_PtDiff");

//     H3D_DetectedV0s_MCdatalike_cutVariation[cutsArrayIterator] = (TH3D*)file_O2Analysis_CutVariation_Datalike_array[cutsArrayIterator]->Get("lambdakzero-analysis-mc/h3dMass"+NamePart[ipart]);

//     double SignalMean, SignalStandardDev;
//     Get_RawYield_vsPt_FeeddownCorrected(hRawYield_vsPt[cutsArrayIterator], H3D_DetectedV0s_MCdatalike_cutVariation[cutsArrayIterator], H1I_SelectedEventCount[cutsArrayIterator], file_O2Analysis_CutVariation_Datalike_array[cutsArrayIterator], &SignalMean, &SignalStandardDev, ipart, pTbins, nbinpT, 0, -1, cutsArrayIterator, SideBandSizeMultiplierModifier_array[0], SignalExtractionType_default, InvMassRebinFactor); ///real line that need to take the place of the one above once it's been written
//     Get_McEfficiency_vsPt(hMcEfficiency_vsPt[cutsArrayIterator], H3D_DetectedV0s_MCtrue_cutVariation[cutsArrayIterator], H1D_TrueV0PtSpectrum_cutVariation[cutsArrayIterator], ipart, pTbins, nbinpT, 0, -1, SideBandSizeMultiplierModifier_array[0], SignalMean, SignalStandardDev);

//     // //debug loop
//     // for(int ibinPt = 1; ibinPt <= nbinpT; ibinPt++){
//     //   cout << "cutID: " << cutsArrayIterator << ", pT bin: " << ibinPt << endl;
//     //   cout << "hRawYield = " << hRawYield_vsPt[cutsArrayIterator]->GetBinContent(ibinPt) << ", hRawYield_error = " << hRawYield_vsPt[cutsArrayIterator]->GetBinError(ibinPt) << endl;
//     //   cout << "hEfficiency = " << hMcEfficiency_vsPt[cutsArrayIterator]->GetBinContent(ibinPt) << ", hEfficiency_error = " << hMcEfficiency_vsPt[cutsArrayIterator]->GetBinError(ibinPt) << endl;
//     // }
//     hCorrectedYield_vsPt[cutsArrayIterator] = (TH1D*)hRawYield_vsPt[cutsArrayIterator]->Clone(Form("hCorrectedYield_vsPt_cutsArrayIterator%i", cutsArrayIterator));
//     hCorrectedYield_vsPt[cutsArrayIterator]->Reset("M");
//     // hCorrectedYield_vsPt[cutsArrayIterator]->Sumw2();
//     // hRawYield_vsPt[cutsArrayIterator]->Sumw2();
//     hCorrectedYield_vsPt[cutsArrayIterator]->Divide(hRawYield_vsPt[cutsArrayIterator],hMcEfficiency_vsPt[cutsArrayIterator], 1., 1., "");
//   // Efficiency errors are binomial-ish;  TGraphAsymmErrors::BayesDivide (https://root.cern.ch/doc/master/classTH1.html#ac782a09c31b4f7de40f8fd4f77efa090)
//   // and https://inspirehep.net/files/57287ac8e45a976ab423f3dd456af694
//   // PWGLF strangeness PAG says to just use option "b" for divide, and just see it as a binomial https://twiki.cern.ch/twiki/bin/view/ALICE/PWGLFPAGSTRANGENESSEfficiency
//   }

//   TH1D* hCorrectedYield_vsPt_REF = hCorrectedYield_vsPt[id_referenceAnalysis]; //this is assuming that's where I indeed put it
//   // TH1D* hSystematicUncertainty[cutsIterator];
//   double hSigmaBarlow[nbinpT];
//   double CorrectedYieldDifference;
//   int id_cutsArray_maxDeviation = 2;
//   for(int ibinPt = 1; ibinPt <= nbinpT; ibinPt++){
//     int PtArrayIterator = ibinPt - 1;
//     double MaxYieldDifference_cutsIteratorArray = 0;
//     for(int cutsArrayIterator = 0; cutsArrayIterator < nCutsIterations; cutsArrayIterator++){
//       CorrectedYieldDifference = abs(hCorrectedYield_vsPt[cutsArrayIterator]->GetBinContent(ibinPt) - hCorrectedYield_vsPt_REF->GetBinContent(ibinPt)); 
//       if (CorrectedYieldDifference > MaxYieldDifference_cutsIteratorArray) {
//         // MaxDifference_cutsIteratorArray = CorrectedYieldDifference;
//         MaxYieldDifference_cutsIteratorArray = CorrectedYieldDifference;
//         id_cutsArray_maxDeviation = cutsArrayIterator;
//       }
//     }

//     double SystUncertainty = MaxYieldDifference_cutsIteratorArray;

//     // Barlow condition for cut variation systematics (Systematic Errors: facts and fictions, by Roger Barlow, https://arxiv.org/abs/hep-ex/0207026)
//     double StatUncertainty_REF = hCorrectedYield_vsPt_REF->GetBinError(ibinPt);
//     double StatUncertainty_MaxDeviationCase = hCorrectedYield_vsPt[id_cutsArray_maxDeviation]->GetBinError(ibinPt);
    
//     hSigmaBarlow[PtArrayIterator] = sqrt(abs(StatUncertainty_MaxDeviationCase*StatUncertainty_MaxDeviationCase - StatUncertainty_REF*StatUncertainty_REF)); //stat error of the difference in the case of subsample
    
//     hSystematicUncertainty_PreBarlow->SetBinContent(ibinPt,SystUncertainty);
//     hSystematicUncertainty_PreBarlow->SetBinError(ibinPt,hSigmaBarlow[PtArrayIterator]);

//     if (SystUncertainty > N_SigmaBarlow*hSigmaBarlow[PtArrayIterator]) { //Could ask for 1Sigma, 4Sigma or whatever depending on how conservative we want to be; one suggested in PWGLF note is 2Sigma
//       hSystematicUncertainty->SetBinContent(ibinPt,SystUncertainty);
//       hSystematicUncertainty->SetBinError(ibinPt,hSigmaBarlow[PtArrayIterator]);
//     }
//     else {
//       hSystematicUncertainty->SetBinContent(ibinPt,0.);
//     }
//     // cout << "pT bin: " << ibinPt  << ", maxDeviation cut ID: " << id_cutsArray_maxDeviation << endl;
//     // cout << "StatUncertainty_MaxDeviationCase = " << StatUncertainty_MaxDeviationCase  << ", StatUncertainty_REF = " << StatUncertainty_REF << endl;
//     // here I ve got a nan (Not a number) problem with StatUncertainty_REF at pT bin 10;
//     // cutID: 0, pT bin: 10
//     // hRawYield = -0.00573374, hRawYield_error = nan    ; negative yield should probably be safeguarded; in fact negative background should be safeguarded

//     cout << "Systematics CutVariation" << endl;
//     cout << "(SystUncertainty, SigmaBarlow) = (" << SystUncertainty << "," << hSigmaBarlow[PtArrayIterator] << ") and check is " << (SystUncertainty > N_SigmaBarlow*hSigmaBarlow[PtArrayIterator]) << endl;
//     cout << "CorrectedYield = " << hCorrectedYield_vsPt_REF->GetBinContent(ibinPt) << endl;
//     cout << "Systematics/CorrectedYield = " << SystUncertainty/hCorrectedYield_vsPt_REF->GetBinContent(ibinPt) << endl;


//   }
// }




void Get_Systematics_Fit_SideBandVariation(TH1D* hSystematicUncertainty, TH1D* hSystematicUncertainty_PreBarlow, int ipart, TH3D* H3D_detectedV0s_DataLike, TH1I* H1I_SelectedEventCount, TH3D* H3D_detectedV0s_TruePt, TH1D* TrueV0PtSpectrum, TFile* file_O2Analysis, double* pTbins, int nbinpT, int InvMassRebinFactor) {



  cout <<"Aaaaaa good this point AAAAAAAAA 0"<<endl;
  // int SignalExtractionType = 1; default already defined beginning of this file

  TH1D *hMcEfficiency_vsPt[N_SideBandSizeVariation];
  TH1D *hRawYield_vsPt[N_SideBandSizeVariation];

  TH1D *hCorrectedYield_vsPt[N_SideBandSizeVariation];


  for(int id_SideBandSizeVariation = 0; id_SideBandSizeVariation < N_SideBandSizeVariation; id_SideBandSizeVariation++){
    cout << "SideBandSizeMultiplierModifier_array[id_SideBandSizeVariation]" << SideBandSizeMultiplierModifier_array[id_SideBandSizeVariation] << endl;

    double SignalMean, SignalStandardDev;
    Get_RawYield_vsPt_FeeddownCorrected(hRawYield_vsPt[id_SideBandSizeVariation], H3D_detectedV0s_DataLike, H1I_SelectedEventCount, file_O2Analysis, &SignalMean, &SignalStandardDev, ipart, pTbins, nbinpT, 0, -1, 0, SideBandSizeMultiplierModifier_array[id_SideBandSizeVariation], SignalExtractionType_default, InvMassRebinFactor); ///real line that need to take the place of the one above once it's been written
    cout << "debug Get_Systematics_Fit_SideBandVariation 11" << endl;
    Get_McEfficiency_vsPt(hMcEfficiency_vsPt[id_SideBandSizeVariation], H3D_detectedV0s_TruePt, TrueV0PtSpectrum, ipart, pTbins, nbinpT, 0, -1, SideBandSizeMultiplierModifier_array[id_SideBandSizeVariation], SignalMean, SignalStandardDev, InvMassRebinFactor);

    // //debug loop
    // for(int ibinPt = 1; ibinPt <= nbinpT; ibinPt++){
    //   cout << "cutID: " << id_SideBandSizeVariation << ", pT bin: " << ibinPt << endl;
    //   cout << "hRawYield = " << hRawYield_vsPt[id_SideBandSizeVariation]->GetBinContent(ibinPt) << ", hRawYield_error = " << hRawYield_vsPt[id_SideBandSizeVariation]->GetBinError(ibinPt) << endl;
    //   cout << "hEfficiency = " << hMcEfficiency_vsPt[id_SideBandSizeVariation]->GetBinContent(ibinPt) << ", hEfficiency_error = " << hMcEfficiency_vsPt[id_SideBandSizeVariation]->GetBinError(ibinPt) << endl;
    // }
    hCorrectedYield_vsPt[id_SideBandSizeVariation] = (TH1D*)hRawYield_vsPt[id_SideBandSizeVariation]->Clone(Form("hCorrectedYield_vsPt_id_SideBandSizeVariation%i", id_SideBandSizeVariation));
    hCorrectedYield_vsPt[id_SideBandSizeVariation]->Reset("M");
    

    hCorrectedYield_vsPt[id_SideBandSizeVariation]->Divide(hRawYield_vsPt[id_SideBandSizeVariation],hMcEfficiency_vsPt[id_SideBandSizeVariation], 1., 1., "");
  // Efficiency errors are binomial-ish;  TGraphAsymmErrors::BayesDivide (https://root.cern.ch/doc/master/classTH1.html#ac782a09c31b4f7de40f8fd4f77efa090)
  // and https://inspirehep.net/files/57287ac8e45a976ab423f3dd456af694
  // PWGLF strangeness PAG says to just use option "b" for divide, and just see it as a binomial https://twiki.cern.ch/twiki/bin/view/ALICE/PWGLFPAGSTRANGENESSEfficiency
  }
  cout <<"Aaaaaa good this point AAAAAAAAA 1"<<endl;

  TH1D* hCorrectedYield_vsPt_REF = hCorrectedYield_vsPt[id_referenceSideband]; //this is assuming that's where I indeed put it
  // TH1D* hSystematicUncertainty[cutsIterator];
  double hSigmaBarlow[nbinpT];
  double CorrectedYieldDifference;
  int id_SideBandSizeVariation_maxDeviation = 0;
  for(int ibinPt = 1; ibinPt <= nbinpT; ibinPt++){
    int PtArrayIterator = ibinPt - 1;
    double MaxYieldDifference_SideBandVariation = 0;
    for(int id_SideBandSizeVariation = 0; id_SideBandSizeVariation < N_SideBandSizeVariation; id_SideBandSizeVariation++){
      CorrectedYieldDifference = abs(hCorrectedYield_vsPt[id_SideBandSizeVariation]->GetBinContent(ibinPt) - hCorrectedYield_vsPt_REF->GetBinContent(ibinPt)); 
      if (CorrectedYieldDifference > MaxYieldDifference_SideBandVariation) {
        // MaxDifference_SideBand = CorrectedYieldDifference;
        MaxYieldDifference_SideBandVariation = CorrectedYieldDifference;
        id_SideBandSizeVariation_maxDeviation = id_SideBandSizeVariation;
      }
    }
    cout << "MaxYieldDifference_SideBandVariation = " << MaxYieldDifference_SideBandVariation << endl;
    double SystUncertainty = MaxYieldDifference_SideBandVariation;
    cout <<"Aaaaaa good this point AAAAAAAAA 2:" << ibinPt <<endl;

    // Barlow condition for cut variation systematics (Systematic Errors: facts and fictions, by Roger Barlow, https://arxiv.org/abs/hep-ex/0207026)
    double StatUncertainty_REF = hCorrectedYield_vsPt_REF->GetBinError(ibinPt);
    double StatUncertainty_MaxDeviationCase = hCorrectedYield_vsPt[id_SideBandSizeVariation_maxDeviation]->GetBinError(ibinPt);
    
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
    cout <<"Aaaaaa good this point AAAAAAAAA 3:" << ibinPt <<endl;

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
  cout <<"Aaaaaa good this point AAAAAAAAA 4"<<endl;

}

void Get_Systematics_Fit_SignalExtractionType(TH1D* hSystematicUncertainty, TH1D* hSystematicUncertainty_PreBarlow, int ipart, TH3D* H3D_detectedV0s_DataLike, TH1I* H1I_SelectedEventCount, TH3D* H3D_detectedV0s_TruePt, TH1D* TrueV0PtSpectrum, TFile* file_O2Analysis, double* pTbins, int nbinpT, int InvMassRebinFactor) {
  int N_SignalExtractionType = 2;

  TH1D *hMcEfficiency_vsPt[N_SignalExtractionType];
  TH1D *hRawYield_vsPt[N_SignalExtractionType];

  TH1D *hCorrectedYield_vsPt[N_SignalExtractionType];


  for(int SignalExtractionType_loopID = 0; SignalExtractionType_loopID < N_SignalExtractionType; SignalExtractionType_loopID++){

    double SignalMean, SignalStandardDev;
    Get_RawYield_vsPt_FeeddownCorrected(hRawYield_vsPt[SignalExtractionType_loopID], H3D_detectedV0s_DataLike, H1I_SelectedEventCount, file_O2Analysis, &SignalMean, &SignalStandardDev, ipart, pTbins, nbinpT, 0, -1, 0, SideBandSizeMultiplierModifier_array[0], SignalExtractionType_loopID, InvMassRebinFactor); ///real line that need to take the place of the one above once it's been written
    Get_McEfficiency_vsPt(hMcEfficiency_vsPt[SignalExtractionType_loopID], H3D_detectedV0s_TruePt, TrueV0PtSpectrum, ipart, pTbins, nbinpT, 0, -1, SideBandSizeMultiplierModifier_array[0], SignalMean, SignalStandardDev, InvMassRebinFactor);

    // //debug loop
    // for(int ibinPt = 1; ibinPt <= nbinpT; ibinPt++){
    //   cout << "cutID: " << SignalExtractionType << ", pT bin: " << ibinPt << endl;
    //   cout << "hRawYield = " << hRawYield_vsPt[SignalExtractionType_loopID]->GetBinContent(ibinPt) << ", hRawYield_error = " << hRawYield_vsPt[SignalExtractionType_loopID]->GetBinError(ibinPt) << endl;
    //   cout << "hEfficiency = " << hMcEfficiency_vsPt[SignalExtractionType_loopID]->GetBinContent(ibinPt) << ", hEfficiency_error = " << hMcEfficiency_vsPt[SignalExtractionType_loopID]->GetBinError(ibinPt) << endl;
    // }
    hCorrectedYield_vsPt[SignalExtractionType_loopID] = (TH1D*)hRawYield_vsPt[SignalExtractionType_loopID]->Clone(Form("hCorrectedYield_vsPt_SignalExtractionType_loopID%i", SignalExtractionType_loopID));
    hCorrectedYield_vsPt[SignalExtractionType_loopID]->Reset("M");
    

    hCorrectedYield_vsPt[SignalExtractionType_loopID]->Divide(hRawYield_vsPt[SignalExtractionType_loopID],hMcEfficiency_vsPt[SignalExtractionType_loopID], 1., 1., "");
  // Efficiency errors are binomial-ish;  TGraphAsymmErrors::BayesDivide (https://root.cern.ch/doc/master/classTH1.html#ac782a09c31b4f7de40f8fd4f77efa090)
  // and https://inspirehep.net/files/57287ac8e45a976ab423f3dd456af694
  // PWGLF strangeness PAG says to just use option "b" for divide, and just see it as a binomial https://twiki.cern.ch/twiki/bin/view/ALICE/PWGLFPAGSTRANGENESSEfficiency
  }

  TH1D* hCorrectedYield_vsPt_REF = hCorrectedYield_vsPt[SignalExtractionType_default]; //this is assuming that's where I indeed put it
  // TH1D* hSystematicUncertainty[cutsIterator];
  double hSigmaBarlow[nbinpT];
  double CorrectedYieldDifference;
  int id_SignalExtractionType_maxDeviation = 0;
  for(int ibinPt = 1; ibinPt <= nbinpT; ibinPt++){
    int PtArrayIterator = ibinPt - 1;
    double MaxYieldDifference_SignalExtractionVariation = 0;
    for(int SignalExtractionType_loopID = 0; SignalExtractionType_loopID < N_SignalExtractionType; SignalExtractionType_loopID++){
      CorrectedYieldDifference = abs(hCorrectedYield_vsPt[SignalExtractionType_loopID]->GetBinContent(ibinPt) - hCorrectedYield_vsPt_REF->GetBinContent(ibinPt)); 
      if (CorrectedYieldDifference > MaxYieldDifference_SignalExtractionVariation) {
        // MaxYieldDifference_SignalExtractionVariation = CorrectedYieldDifference;
        MaxYieldDifference_SignalExtractionVariation = CorrectedYieldDifference;
        id_SignalExtractionType_maxDeviation = SignalExtractionType_loopID;
      }
    }

    double SystUncertainty = MaxYieldDifference_SignalExtractionVariation;

    // Barlow condition for cut variation systematics (Systematic Errors: facts and fictions, by Roger Barlow, https://arxiv.org/abs/hep-ex/0207026)
    double StatUncertainty_REF = hCorrectedYield_vsPt_REF->GetBinError(ibinPt);
    double StatUncertainty_MaxDeviationCase = hCorrectedYield_vsPt[id_SignalExtractionType_maxDeviation]->GetBinError(ibinPt);
    
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


void Get_Systematics_Fit_Binning(TH1D* hSystematicUncertainty, TH1D* hSystematicUncertainty_PreBarlow, int ipart, TH3D* H3D_detectedV0s_DataLike, TH1I* H1I_SelectedEventCount, TH3D* H3D_detectedV0s_TruePt, TH1D* TrueV0PtSpectrum, TFile* file_O2Analysis, double* pTbins, int nbinpT) {

  // int SignalExtractionType = 1; default already defined beginning of this file

  TH1D *hMcEfficiency_vsPt[N_BinningVariation];
  TH1D *hRawYield_vsPt[N_BinningVariation];

  TH1D *hCorrectedYield_vsPt[N_BinningVariation];


  for(int id_BinningVariation = 0; id_BinningVariation < N_BinningVariation; id_BinningVariation++){
    cout << "InvMassRebinFactor_variation_array[id_BinningVariation]" << InvMassRebinFactor_variation_array[ipart][id_BinningVariation] << endl;
    double SignalMean, SignalStandardDev;
    Get_RawYield_vsPt_FeeddownCorrected(hRawYield_vsPt[id_BinningVariation], H3D_detectedV0s_DataLike, H1I_SelectedEventCount, file_O2Analysis, &SignalMean, &SignalStandardDev, ipart, pTbins, nbinpT, 0, -1, 0, 0, SignalExtractionType_default, InvMassRebinFactor_variation_array[ipart][id_BinningVariation]); ///real line that need to take the place of the one above once it's been written
    Get_McEfficiency_vsPt(hMcEfficiency_vsPt[id_BinningVariation], H3D_detectedV0s_TruePt, TrueV0PtSpectrum, ipart, pTbins, nbinpT, 0, -1, SideBandSizeMultiplierModifier_array[0], SignalMean, SignalStandardDev, InvMassRebinFactor_variation_array[ipart][id_BinningVariation]);

    // //debug loop
    // for(int ibinPt = 1; ibinPt <= nbinpT; ibinPt++){
    //   cout << "cutID: " << id_BinningVariation << ", pT bin: " << ibinPt << endl;
    //   cout << "hRawYield = " << hRawYield_vsPt[id_BinningVariation]->GetBinContent(ibinPt) << ", hRawYield_error = " << hRawYield_vsPt[id_BinningVariation]->GetBinError(ibinPt) << endl;
    //   cout << "hEfficiency = " << hMcEfficiency_vsPt[id_BinningVariation]->GetBinContent(ibinPt) << ", hEfficiency_error = " << hMcEfficiency_vsPt[id_BinningVariation]->GetBinError(ibinPt) << endl;
    // }
    hCorrectedYield_vsPt[id_BinningVariation] = (TH1D*)hRawYield_vsPt[id_BinningVariation]->Clone(Form("hCorrectedYield_vsPt_id_BinningVariation%i", id_BinningVariation));
    hCorrectedYield_vsPt[id_BinningVariation]->Reset("M");
    hCorrectedYield_vsPt[id_BinningVariation]->Divide(hRawYield_vsPt[id_BinningVariation],hMcEfficiency_vsPt[id_BinningVariation], 1., 1., "");
  // Efficiency errors are binomial-ish;  TGraphAsymmErrors::BayesDivide (https://root.cern.ch/doc/master/classTH1.html#ac782a09c31b4f7de40f8fd4f77efa090)
  // and https://inspirehep.net/files/57287ac8e45a976ab423f3dd456af694
  // PWGLF strangeness PAG says to just use option "b" for divide, and just see it as a binomial https://twiki.cern.ch/twiki/bin/view/ALICE/PWGLFPAGSTRANGENESSEfficiency
  }

  // TH1D* hCorrectedYield_vsPt_REF = hCorrectedYield_vsPt[id_referenceSideband]; //this is assuming that's where I indeed put it
  TH1D *hRawYield_vsPt_REF;
  TH1D *hMcEfficiency_vsPt_REF;
  double SignalMean, SignalStandardDev;
  Get_RawYield_vsPt_FeeddownCorrected(hRawYield_vsPt_REF, H3D_detectedV0s_DataLike, H1I_SelectedEventCount, file_O2Analysis, &SignalMean, &SignalStandardDev, ipart, pTbins, nbinpT, 0, -1, 0, 0, SignalExtractionType_default, InvMassRebinFactor_standard[ipart]); ///real line that need to take the place of the one above once it's been written
  Get_McEfficiency_vsPt(hMcEfficiency_vsPt_REF, H3D_detectedV0s_TruePt, TrueV0PtSpectrum, ipart, pTbins, nbinpT, 0, -1, SideBandSizeMultiplierModifier_array[0], SignalMean, SignalStandardDev, InvMassRebinFactor_standard[ipart]);
  TH1D* hCorrectedYield_vsPt_REF = (TH1D*)hRawYield_vsPt_REF->Clone("hCorrectedYield_vsPt_REF");

  hCorrectedYield_vsPt_REF->Reset("M");
  hCorrectedYield_vsPt_REF->Divide(hRawYield_vsPt_REF,hMcEfficiency_vsPt_REF, 1., 1., "");

  // TH1D* hSystematicUncertainty[cutsIterator];
  double hSigmaBarlow[nbinpT];
  double CorrectedYieldDifference;
  int id_BinningVariation_maxDeviation = 0;
  for(int ibinPt = 1; ibinPt <= nbinpT; ibinPt++){
    int PtArrayIterator = ibinPt - 1;
    double MaxYieldDifference_BinningVariation = 0;
    for(int id_BinningVariation = 0; id_BinningVariation < N_BinningVariation; id_BinningVariation++){
      CorrectedYieldDifference = abs(hCorrectedYield_vsPt[id_BinningVariation]->GetBinContent(ibinPt) - hCorrectedYield_vsPt_REF->GetBinContent(ibinPt)); 
      if (CorrectedYieldDifference > MaxYieldDifference_BinningVariation) {
        // MaxDifference_SideBand = CorrectedYieldDifference;
        MaxYieldDifference_BinningVariation = CorrectedYieldDifference;
        id_BinningVariation_maxDeviation = id_BinningVariation;
      }
    }
    cout << "MaxYieldDifference_BinningVariation = " << MaxYieldDifference_BinningVariation << endl;
    double SystUncertainty = MaxYieldDifference_BinningVariation;

    // Barlow condition for cut variation systematics (Systematic Errors: facts and fictions, by Roger Barlow, https://arxiv.org/abs/hep-ex/0207026)
    double StatUncertainty_REF = hCorrectedYield_vsPt_REF->GetBinError(ibinPt);
    double StatUncertainty_MaxDeviationCase = hCorrectedYield_vsPt[id_BinningVariation_maxDeviation]->GetBinError(ibinPt);
    
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


void Systematics_CutVariations_Graphs(TH1D* &hstat, TH1D* &hsyst, TH1D* &hsystCorr, int ipart, TString* SaveAs_Title, TString* &texXtitle, TString* &texYtitle, double* pTbins, int nbinpT) {

  TH3D* H3D_detectedV0s_TruePt = (TH3D*)file_O2Analysis->Get("lambdakzero-analysis-mc/h3dMass"+NamePart[ipart]+"_MC_truePt");
  TH3D* H3D_detectedV0s_DataLike = (TH3D*)file_O2Analysis->Get("lambdakzero-analysis-mc/h3dMass"+NamePart[ipart]);
  TH1D* TrueV0PtSpectrum = (TH1D*)file_O2Analysis->Get("lambdakzero-particle-count-mc/h"+NamePart[ipart]+"Count_PtDiff");
  TH1I* H1I_SelectedEventCount = (TH1I*)file_O2Analysis->Get("lambdakzero-analysis-mc/hSelectedEventCounter");

  TH1D* hSystematicUncertainty_Fit_SidebandVariation = new TH1D("hSystematicUncertainty_Fit_SidebandVariation", "hSystematicUncertainty_Fit_SidebandVariation", nbinpT, pTbins);
  TH1D* hSystematicUncertainty_Fit_SidebandVariation_PreBarlow = new TH1D("hSystematicUncertainty_Fit_SidebandVariation_PreBarlow", "hSystematicUncertainty_Fit_SidebandVariation_PreBarlow", nbinpT, pTbins);
  TH1D* hSystematicUncertainty_Fit_SignalExtractionType = new TH1D("hSystematicUncertainty_Fit_SignalExtractionType", "hSystematicUncertainty_Fit_SignalExtractionType", nbinpT, pTbins);
  TH1D* hSystematicUncertainty_Fit_SignalExtractionType_PreBarlow = new TH1D("hSystematicUncertainty_Fit_SignalExtractionType_PreBarlow", "hSystematicUncertainty_Fit_SignalExtractionType_PreBarlow", nbinpT, pTbins);
  TH1D* hSystematicUncertainty_Fit_BinningVariation = new TH1D("hSystematicUncertainty_Fit_BinningVariation", "hSystematicUncertainty_Fit_BinningVariation", nbinpT, pTbins);
  TH1D* hSystematicUncertainty_Fit_BinningVariation_PreBarlow = new TH1D("hSystematicUncertainty_Fit_BinningVariation_PreBarlow", "hSystematicUncertainty_Fit_BinningVariation_PreBarlow", nbinpT, pTbins);

  TH1D* hCorrectedYield_vsPt = new TH1D("hCorrectedYield_vsPt", "hCorrectedYield_vsPt", nbinpT, pTbins);
  TH1D* hSystematicUncertainty_Total = new TH1D("hSystematicUncertainty_Total", "hSystematicUncertainty_Total", nbinpT, pTbins);
  TH1D* hSystematicUncertainty_Total_PreBarlow = new TH1D("hSystematicUncertainty_Total_PreBarlow", "hSystematicUncertainty_Total_PreBarlow", nbinpT, pTbins);


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
  
  int InvMassRebinFactor = InvMassRebinFactor_standard[ipart];

  double SignalMean, SignalStandardDev;

  cout << "Calculating Raw Yield for Ref" << endl;
  Get_RawYield_vsPt_FeeddownCorrected(hRawYield_vsPt, H3D_detectedV0s_DataLike, H1I_SelectedEventCount, file_O2Analysis, &SignalMean, &SignalStandardDev, ipart, pTbins, nbinpT, 0, -1, 0, SideBandSizeMultiplierModifier_array[0], SignalExtractionType_default, InvMassRebinFactor); ///real line that need to take the place of the one above once it's been written
  cout << "Calculating Efficiency" << endl;
  Get_McEfficiency_vsPt(hMcEfficiency_vsPt, H3D_detectedV0s_TruePt, TrueV0PtSpectrum, ipart, pTbins, nbinpT, 0, -1, SideBandSizeMultiplierModifier_array[0], SignalMean, SignalStandardDev, InvMassRebinFactor);

  cout << "Calculating Systematics - Fit sideband variation" << endl;
  Get_Systematics_Fit_SideBandVariation(hSystematicUncertainty_Fit_SidebandVariation, hSystematicUncertainty_Fit_SidebandVariation_PreBarlow, ipart, H3D_detectedV0s_DataLike, H1I_SelectedEventCount, H3D_detectedV0s_TruePt, TrueV0PtSpectrum, file_O2Analysis, pTbins, nbinpT, InvMassRebinFactor);
  cout << "Calculating Systematics - Fit signal extraction type" << endl;
  Get_Systematics_Fit_SignalExtractionType(hSystematicUncertainty_Fit_SignalExtractionType, hSystematicUncertainty_Fit_SignalExtractionType_PreBarlow, ipart, H3D_detectedV0s_DataLike, H1I_SelectedEventCount, H3D_detectedV0s_TruePt, TrueV0PtSpectrum, file_O2Analysis, pTbins, nbinpT, InvMassRebinFactor); //very large systematics as is
  cout << "Calculating Systematics - Fit binning" << endl;
  Get_Systematics_Fit_Binning(hSystematicUncertainty_Fit_BinningVariation, hSystematicUncertainty_Fit_BinningVariation_PreBarlow, ipart, H3D_detectedV0s_DataLike, H1I_SelectedEventCount, H3D_detectedV0s_TruePt, TrueV0PtSpectrum, file_O2Analysis, pTbins, nbinpT);

  //systematics on cut variations over the N_variableCut variables that are being tested
  TH1D* hSystematicUncertainty_CutVariation[N_variableCut];
  TH1D* hSystematicUncertainty_CutVariation_PreBarlow[N_variableCut];
  for (int i_variableCut = 0; i_variableCut < N_variableCut; i_variableCut++) {
    hSystematicUncertainty_CutVariation[i_variableCut] = new TH1D("hSystematicUncertainty_CutVariation_"+VariableCut[i_variableCut], "hSystematicUncertainty_CutVariation_"+VariableCut[i_variableCut], nbinpT, pTbins);
    hSystematicUncertainty_CutVariation_PreBarlow[i_variableCut] = new TH1D("hSystematicUncertainty_CutVariation_PreBarlow_"+VariableCut[i_variableCut], "hSystematicUncertainty_CutVariation_PreBarlow_"+VariableCut[i_variableCut], nbinpT, pTbins);
    hSystematicUncertainty_CutVariation[i_variableCut]->Sumw2();
    hSystematicUncertainty_CutVariation_PreBarlow[i_variableCut]->Sumw2();
    cout << "Calculating Systematics - Cut variation number " << i_variableCut << endl;
    Get_Systematics_OneCut(hSystematicUncertainty_CutVariation[i_variableCut], hSystematicUncertainty_CutVariation_PreBarlow[i_variableCut], ipart, file_O2Analysis_CutVariation_Collection_array[i_variableCut], nCutsIterations_array[i_variableCut], pTbins, nbinpT, InvMassRebinFactor);
  }

  // Add uncertainties in quadrature: multiply hists by themselves then add
  hSystematicUncertainty_Fit_SidebandVariation->Multiply(hSystematicUncertainty_Fit_SidebandVariation);
  hSystematicUncertainty_Fit_SignalExtractionType->Multiply(hSystematicUncertainty_Fit_SignalExtractionType);
  hSystematicUncertainty_Fit_BinningVariation->Multiply(hSystematicUncertainty_Fit_BinningVariation);


  hSystematicUncertainty_Fit_SidebandVariation_PreBarlow->Multiply(hSystematicUncertainty_Fit_SidebandVariation_PreBarlow);
  hSystematicUncertainty_Fit_SignalExtractionType_PreBarlow->Multiply(hSystematicUncertainty_Fit_SignalExtractionType_PreBarlow);
  hSystematicUncertainty_Fit_BinningVariation_PreBarlow->Multiply(hSystematicUncertainty_Fit_BinningVariation_PreBarlow);

  TH1D* hSystematicUncertainty_CutVariation_Total = new TH1D("hSystematicUncertainty", "hSystematicUncertainty", nbinpT, pTbins);
  TH1D* hSystematicUncertainty_CutVariation_Total_PreBarlow = new TH1D("hSystematicUncertainty_PreBarlow", "hSystematicUncertainty_PreBarlow", nbinpT, pTbins);
  hSystematicUncertainty_CutVariation_Total->Sumw2();
  hSystematicUncertainty_CutVariation_Total_PreBarlow->Sumw2();
  for (int i_variableCut = 0; i_variableCut < N_variableCut; i_variableCut++) {
    hSystematicUncertainty_CutVariation[i_variableCut]->Multiply(hSystematicUncertainty_CutVariation[i_variableCut]);
    hSystematicUncertainty_CutVariation_Total->Add(hSystematicUncertainty_CutVariation[i_variableCut]);
    hSystematicUncertainty_CutVariation_PreBarlow[i_variableCut]->Multiply(hSystematicUncertainty_CutVariation_PreBarlow[i_variableCut]);
    hSystematicUncertainty_CutVariation_Total_PreBarlow->Add(hSystematicUncertainty_CutVariation_PreBarlow[i_variableCut]);
  }
  hSystematicUncertainty_Total->Add(hSystematicUncertainty_Fit_SidebandVariation);
  hSystematicUncertainty_Total->Add(hSystematicUncertainty_Fit_SignalExtractionType);
  hSystematicUncertainty_Total->Add(hSystematicUncertainty_Fit_BinningVariation);
  hSystematicUncertainty_Total->Add(hSystematicUncertainty_CutVariation_Total);

  // hSystematicUncertainty_Total_PreBarlow->Add(hSystematicUncertainty_CutVariation_PreBarlow);
  hSystematicUncertainty_Total_PreBarlow->Add(hSystematicUncertainty_Fit_SidebandVariation_PreBarlow);
  hSystematicUncertainty_Total_PreBarlow->Add(hSystematicUncertainty_Fit_SignalExtractionType_PreBarlow);
  hSystematicUncertainty_Total_PreBarlow->Add(hSystematicUncertainty_Fit_BinningVariation_PreBarlow);
  // hSystematicUncertainty_Total_PreBarlow->Add(hSystematicUncertainty_CutVariation_Total_PreBarlow);

  // back to Sigma from Sigma^2
  for(int ibinPt = 1; ibinPt <= nbinpT; ibinPt++){
      // wrong errors; can probably get away with not looking at the errors of those systematic errors
    if (abs(hSystematicUncertainty_Total->GetBinContent(ibinPt)) > 1E-9 ) { //can't divide by 0
      // hSystematicUncertainty_Total->SetBinError(ibinPt, 1./2 * 1./sqrt(hSystematicUncertainty_Total->GetBinContent(ibinPt)) * hSystematicUncertainty_Total->GetBinError(ibinPt)); // SigmaTot^2 = Sum(Sigma^2) -> error(SigmaTot) = 1/(2*SigmaTot)*Error(SigmaTot^2)  
      hSystematicUncertainty_Total->SetBinContent(ibinPt,sqrt(hSystematicUncertainty_Total->GetBinContent(ibinPt)));
    }
    if (abs(hSystematicUncertainty_Total_PreBarlow->GetBinContent(ibinPt)) > 1E-9 ) {
      // hSystematicUncertainty_Total_PreBarlow->SetBinError(ibinPt, 1./2 * 1./sqrt(hSystematicUncertainty_Total_PreBarlow->GetBinContent(ibinPt)) * hSystematicUncertainty_Total_PreBarlow->GetBinError(ibinPt));    
      hSystematicUncertainty_Total_PreBarlow->SetBinContent(ibinPt,sqrt(hSystematicUncertainty_Total_PreBarlow->GetBinContent(ibinPt)));
    }
    if (abs(hSystematicUncertainty_CutVariation_Total->GetBinContent(ibinPt)) > 1E-9 ) {
      // hSystematicUncertainty_CutVariation_Total->SetBinError(ibinPt, 1./2 * 1./sqrt(hSystematicUncertainty_CutVariation_Total->GetBinContent(ibinPt)) * hSystematicUncertainty_CutVariation_Total->GetBinError(ibinPt));    
      hSystematicUncertainty_CutVariation_Total->SetBinContent(ibinPt,sqrt(hSystematicUncertainty_CutVariation_Total->GetBinContent(ibinPt)));
    }
    if (abs(hSystematicUncertainty_CutVariation_Total_PreBarlow->GetBinContent(ibinPt)) > 1E-9 ) {
      hSystematicUncertainty_CutVariation_Total_PreBarlow->SetBinError(ibinPt, 0.);  
      // hSystematicUncertainty_CutVariation_Total_PreBarlow->SetBinError(ibinPt, 1./2 * 1./sqrt(hSystematicUncertainty_CutVariation_Total_PreBarlow->GetBinContent(ibinPt)) * hSystematicUncertainty_CutVariation_Total_PreBarlow->GetBinError(ibinPt));    
      hSystematicUncertainty_CutVariation_Total_PreBarlow->SetBinContent(ibinPt,sqrt(hSystematicUncertainty_CutVariation_Total_PreBarlow->GetBinContent(ibinPt)));
    }
    if (abs(hSystematicUncertainty_Fit_SidebandVariation->GetBinContent(ibinPt)) > 1E-9 ) {
      // hSystematicUncertainty_Fit_SidebandVariation->SetBinError(ibinPt, 1./2 * 1./sqrt(hSystematicUncertainty_Fit_SidebandVariation->GetBinContent(ibinPt)) * hSystematicUncertainty_Fit_SidebandVariation->GetBinError(ibinPt));    
      hSystematicUncertainty_Fit_SidebandVariation->SetBinContent(ibinPt,sqrt(hSystematicUncertainty_Fit_SidebandVariation->GetBinContent(ibinPt)));
    }
    if (abs(hSystematicUncertainty_Fit_SidebandVariation_PreBarlow->GetBinContent(ibinPt)) > 1E-9 ) {
      // hSystematicUncertainty_Fit_SidebandVariation_PreBarlow->SetBinError(ibinPt, 1./2 * 1./sqrt(hSystematicUncertainty_Fit_SidebandVariation_PreBarlow->GetBinContent(ibinPt)) * hSystematicUncertainty_Fit_SidebandVariation_PreBarlow->GetBinError(ibinPt));    
      hSystematicUncertainty_Fit_SidebandVariation_PreBarlow->SetBinContent(ibinPt,sqrt(hSystematicUncertainty_Fit_SidebandVariation_PreBarlow->GetBinContent(ibinPt)));
    }
    if (abs(hSystematicUncertainty_Fit_SignalExtractionType->GetBinContent(ibinPt)) > 1E-9 ) {
      // hSystematicUncertainty_Fit_SignalExtractionType->SetBinError(ibinPt, 1./2 * 1./sqrt(hSystematicUncertainty_Fit_SignalExtractionType->GetBinContent(ibinPt)) * hSystematicUncertainty_Fit_SignalExtractionType->GetBinError(ibinPt));    
      hSystematicUncertainty_Fit_SignalExtractionType->SetBinContent(ibinPt,sqrt(hSystematicUncertainty_Fit_SignalExtractionType->GetBinContent(ibinPt)));
    }
    if (abs(hSystematicUncertainty_Fit_SignalExtractionType_PreBarlow->GetBinContent(ibinPt)) > 1E-9 ) {
      // hSystematicUncertainty_Fit_SignalExtractionType_PreBarlow->SetBinError(ibinPt, 1./2 * 1./sqrt(hSystematicUncertainty_Fit_SignalExtractionType_PreBarlow->GetBinContent(ibinPt)) * hSystematicUncertainty_Fit_SignalExtractionType_PreBarlow->GetBinError(ibinPt));    
      hSystematicUncertainty_Fit_SignalExtractionType_PreBarlow->SetBinContent(ibinPt,sqrt(hSystematicUncertainty_Fit_SignalExtractionType_PreBarlow->GetBinContent(ibinPt)));
    }
    if (abs(hSystematicUncertainty_Fit_BinningVariation->GetBinContent(ibinPt)) > 1E-9 ) {
      // hSystematicUncertainty_Fit_BinningVariation->SetBinError(ibinPt, 1./2 * 1./sqrt(hSystematicUncertainty_Fit_BinningVariation->GetBinContent(ibinPt)) * hSystematicUncertainty_Fit_BinningVariation->GetBinError(ibinPt));    
      hSystematicUncertainty_Fit_BinningVariation->SetBinContent(ibinPt,sqrt(hSystematicUncertainty_Fit_BinningVariation->GetBinContent(ibinPt)));
    }
    if (abs(hSystematicUncertainty_Fit_BinningVariation_PreBarlow->GetBinContent(ibinPt)) > 1E-9 ) {
      // hSystematicUncertainty_Fit_BinningVariation_PreBarlow->SetBinError(ibinPt, 1./2 * 1./sqrt(hSystematicUncertainty_Fit_BinningVariation_PreBarlow->GetBinContent(ibinPt)) * hSystematicUncertainty_Fit_BinningVariation_PreBarlow->GetBinError(ibinPt));    
      hSystematicUncertainty_Fit_BinningVariation_PreBarlow->SetBinContent(ibinPt,sqrt(hSystematicUncertainty_Fit_BinningVariation_PreBarlow->GetBinContent(ibinPt)));
    }
  }
  //divide by corrected yield for ratio

  hCorrectedYield_vsPt->Sumw2();
  hCorrectedYield_vsPt->Divide(hRawYield_vsPt,hMcEfficiency_vsPt, 1., 1., "");

  // for (int i_variableCut = 0; i_variableCut < N_variableCut; i_variableCut++) {
  //   hSystematicUncertainty_CutVariation_Total->Divide(hCorrectedYield_vsPt);//get it as a ratio of ref corrected yield  
  //   hSystematicUncertainty_CutVariation_Total_PreBarlow->Divide(hCorrectedYield_vsPt);//get it as a ratio of ref corrected yield  
  // }
  hSystematicUncertainty_CutVariation_Total->Divide(hCorrectedYield_vsPt);//get it as a ratio of ref corrected yield  
  hSystematicUncertainty_CutVariation_Total_PreBarlow->Divide(hCorrectedYield_vsPt);//get it as a ratio of ref corrected yield  
  hSystematicUncertainty_Fit_SidebandVariation->Divide(hCorrectedYield_vsPt); //get it as a ratio of ref corrected yield
  hSystematicUncertainty_Fit_SidebandVariation_PreBarlow->Divide(hCorrectedYield_vsPt); //get it as a ratio of ref corrected yield
  hSystematicUncertainty_Fit_SignalExtractionType->Divide(hCorrectedYield_vsPt); //get it as a ratio of ref corrected yield
  hSystematicUncertainty_Fit_SignalExtractionType_PreBarlow->Divide(hCorrectedYield_vsPt); //get it as a ratio of ref corrected yield
  hSystematicUncertainty_Fit_BinningVariation->Divide(hCorrectedYield_vsPt); //get it as a ratio of ref corrected yield
  hSystematicUncertainty_Fit_BinningVariation_PreBarlow->Divide(hCorrectedYield_vsPt); //get it as a ratio of ref corrected yield


  hSystematicUncertainty_Total->Divide(hCorrectedYield_vsPt); //get it as a ratio of ref corrected yield
  hSystematicUncertainty_Total_PreBarlow->Divide(hCorrectedYield_vsPt); //get it as a ratio of ref corrected yield
  
  cout << "POSTbarlow ---------- TOTAL --------- Syst = " << hSystematicUncertainty_Total->GetBinContent(1) << ", SystError = " << hSystematicUncertainty_Total->GetBinError(1) << endl;
  cout << "POSTbarlow ---------- Binning Var --------- Syst = " << hSystematicUncertainty_Fit_BinningVariation->GetBinContent(1) << ", SystError = " << hSystematicUncertainty_Fit_BinningVariation->GetBinError(1) << endl;
  cout << "PREBARLOW ---------- TOTAL --------- Syst = " << hSystematicUncertainty_Total_PreBarlow->GetBinContent(1) << ", SystError = " << hSystematicUncertainty_Total_PreBarlow->GetBinError(1) << endl;
  cout << "PREBARLOW ---------- Binning Var --------- Syst = " << hSystematicUncertainty_Fit_BinningVariation_PreBarlow->GetBinContent(1) << ", SystError = " << hSystematicUncertainty_Fit_BinningVariation_PreBarlow->GetBinError(1) << endl;

  //////////////////////////////////
  //////////// Plots ///////////////
  //////////////////////////////////
  //Post Barlow
  TCanvas *canvasSyst = new TCanvas ("canvasSyst", NamehistoInvMass[ipart], 1200, 800);
  canvasSyst->cd(0); //ATTENTION: +1 because canvas divide starts counting from 1
  TH1 *hFrame = canvasSyst->DrawFrame(pTbins[0],0,pTbins[nbinpT],1.2*hSystematicUncertainty_Total->GetMaximum());
  hFrame->SetXTitle("#it{p}_{T} (GeV/#it{c})");
  hFrame->SetYTitle("uncertainty (fraction of spectrum)");

  // hSystematicUncertainty_Total->Draw("SAME");
  // hSystematicUncertainty_CutVariation_Total->Draw("SAME");
  hSystematicUncertainty_Fit_SidebandVariation->Draw("SAME");
  hSystematicUncertainty_Fit_SignalExtractionType->Draw("SAME");
  hSystematicUncertainty_Fit_BinningVariation->Draw("SAME");

  // hSystematicUncertainty_Total->Draw("HIST SAME P");
  // hSystematicUncertainty_CutVariation_Total->Draw("HIST SAME P");
  // hSystematicUncertainty_Fit_SidebandVariation->Draw("HIST SAME P");
  // hSystematicUncertainty_Fit_SignalExtractionType->Draw("HIST SAME P");
  // hSystematicUncertainty_Fit_BinningVariation->Draw("HIST SAME P");

  int icolor_Fit_BinningVariation = 0;
  int icolor_Fit_SidebandVariation = 1;
  int icolor_Fit_SignalExtractionType = 2;
  int icolor_CutVariation = 3;
  int icolor_Total = 4;


  hSystematicUncertainty_Fit_BinningVariation->SetMarkerStyle(markers[0]);
  hSystematicUncertainty_Fit_BinningVariation->SetMarkerColor(colors [icolor_Fit_BinningVariation]);
  hSystematicUncertainty_Fit_BinningVariation->SetLineColor  (colors [icolor_Fit_BinningVariation]);

  hSystematicUncertainty_Fit_SidebandVariation->SetMarkerStyle(markers[1]);
  hSystematicUncertainty_Fit_SidebandVariation->SetMarkerColor(colors [icolor_Fit_SidebandVariation]);
  hSystematicUncertainty_Fit_SidebandVariation->SetLineColor  (colors [icolor_Fit_SidebandVariation]);

  hSystematicUncertainty_Fit_SignalExtractionType->SetMarkerStyle(markers[2]);
  hSystematicUncertainty_Fit_SignalExtractionType->SetMarkerColor(colors [icolor_Fit_SignalExtractionType]);
  hSystematicUncertainty_Fit_SignalExtractionType->SetLineColor  (colors [icolor_Fit_SignalExtractionType]);

  // hSystematicUncertainty_CutVariation_Total->SetMarkerStyle(markers[3]);
  // hSystematicUncertainty_CutVariation_Total->SetMarkerColor(colors [icolor_CutVariation]);
  // hSystematicUncertainty_CutVariation_Total->SetLineColor  (colors [icolor_CutVariation]);

  hSystematicUncertainty_Total->SetMarkerStyle(markers[4]);
  hSystematicUncertainty_Total->SetMarkerColor(colors [icolor_Total]);
  hSystematicUncertainty_Total->SetLineColor  (colors [icolor_Total]);

  //Legend, if needed
  TLegend * leg = new TLegend(0.62, 0.6, 0.87, 0.8);
  leg->AddEntry(hSystematicUncertainty_CutVariation_Total,     "CutVariation_total",   "LP");
  leg->AddEntry(hSystematicUncertainty_Fit_SidebandVariation,     "SidebandVariation",  "LP");
  leg->AddEntry(hSystematicUncertainty_Fit_SignalExtractionType, "SignalExtractionType",    "LP" );
  leg->AddEntry(hSystematicUncertainty_Fit_BinningVariation, "BinningVariation",    "LP" );
  // leg->AddEntry(hSystematicUncertainty_Total, "Total",    "LP" );
  // leg->SetFillColor(0);
  leg->SetTextSize(gStyle->GetTextSize()*0.8);
  // leg->Draw("same");

  canvasSyst->SaveAs("Systematics_CutVariations_vsPt_"+NamePart[ipart]+".pdf");

  /////////////////////////// Pre Barlow ///////////////////////////////////
  TCanvas *canvasSyst_preBarlow = new TCanvas ("canvasSyst_preBarlow", NamehistoInvMass[ipart], 1200, 800);
  canvasSyst_preBarlow->cd(0); //ATTENTION: +1 because canvas divide starts counting from 1
  TH1 *hFrame_preBarlow = canvasSyst_preBarlow->DrawFrame(pTbins[0],0,pTbins[nbinpT],1.2*hSystematicUncertainty_Total_PreBarlow->GetMaximum());
  hFrame_preBarlow->SetXTitle("#it{p}_{T} (GeV/#it{c})");
  hFrame_preBarlow->SetYTitle("uncertainty (fraction of spectrum)");

  hSystematicUncertainty_CutVariation_Total_PreBarlow->Draw("HISTSAME");
  hSystematicUncertainty_Fit_SignalExtractionType_PreBarlow->Draw("HISTSAME");
  hSystematicUncertainty_Fit_SidebandVariation_PreBarlow->Draw("HISTSAME");
  hSystematicUncertainty_Fit_BinningVariation_PreBarlow->Draw("HISTSAME");
  hSystematicUncertainty_Total_PreBarlow->Draw("HISTSAME");

  hSystematicUncertainty_CutVariation_Total_PreBarlow->SetMarkerStyle(markers[0]);
  hSystematicUncertainty_CutVariation_Total_PreBarlow->SetMarkerColor(colors [icolor_CutVariation]);
  hSystematicUncertainty_CutVariation_Total_PreBarlow->SetLineColor  (colors [icolor_CutVariation]);

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
  leg_preBarlow->AddEntry(hSystematicUncertainty_CutVariation_Total_PreBarlow,     "CutVariation_total PreBarlow",   "L");
  leg_preBarlow->AddEntry(hSystematicUncertainty_Fit_SidebandVariation_PreBarlow,     "SidebandVariation PreBarlow",  "L");
  leg_preBarlow->AddEntry(hSystematicUncertainty_Fit_SignalExtractionType_PreBarlow, "SignalExtractionType PreBarlow",    "L" );
  leg_preBarlow->AddEntry(hSystematicUncertainty_Fit_BinningVariation_PreBarlow, "BinningVariation PreBarlow",    "L" );
  leg_preBarlow->AddEntry(hSystematicUncertainty_Total_PreBarlow, "Total PreBarlow",    "L" );
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
  int nbinx = hstat->GetNbinsX();
  
  for(int ibinx = 1; ibinx <= nbinx; ibinx++){
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

double myLevyPtXsi(double *pt, double *par) { //https://twiki.cern.ch/twiki/bin/view/ALICE/PWG2SpectraTopical7TeVLambdas#Levy_fit

  double lMass = 1.32171; //from PDG
  double ldNdy = par[0];
  // double l2pi = 2*TMath::Pi();
  double lTemp = par[1];
  double lPower = par[2]; // the n in the pp900GeV paper
  // double lEventCount= par[3]; // the n in the pp900GeV paper
  double lBigCoef = ((lPower-1)*(lPower-2)) / (lPower*lTemp*(lPower*lTemp+lMass*(lPower-2)));
  double lInPower = 1 + (TMath::Sqrt(pt[0]*pt[0]+lMass*lMass)-lMass) / (lPower*lTemp);

  // return lEventCount * ldNdy * pt[0] * lBigCoef * TMath::Power(lInPower,(-1)*lPower);
  return ldNdy * pt[0] * lBigCoef * TMath::Power(lInPower,(-1)*lPower);
}

void Get_FeedDown(TH1D * H1D_FeedDownCorrection, int ipart_option, TFile *file_O2Analysis, double* pTbins, int nbinpT) {
  TH2D* H2D_Feeddown_numerator = (TH2D*)file_O2Analysis->Get("lambdakzero-analysis-mc/h"+NamePart[ipart_option]+"FeedDownMatrix");
  TH1D* H1D_Feeddown_denominator = (TH1D*)file_O2Analysis->Get("lambdakzero-particle-count-mc/h"+NamePart_multiStrange[ipart_option-1]+"Count_PtDiff");
  // Rebinning of those two histograms
  // TH1D* H1D_Feeddown_denominator_rebinned = (TH1D*)H1D_Feeddown_denominator->Rebin(nbinpT,"H1D_Feeddown_denominator_rebinned_"+NamePart_multiStrange[ipart_option],pTbins); could do it that way for 1 dim hist but we decide to do it the same way as for the 2d hist here

  // int nbinpT_j = nbinpT;
  // double* pTbins_j = pTbins;
  int nbinpT_j = nbinpT_Xi;
  double* pTbins_j = pTbinsXi;

  TH2D* H2D_Feeddown_numerator_rebinned = new TH2D("H2D_Feeddown_numerator_rebinned_"+NamePart[ipart_option], "H2D_Feeddown_numerator_rebinned_"+NamePart[ipart_option], nbinpT, pTbins, nbinpT_j, pTbinsXi);
  TH1D* H1D_Feeddown_denominator_rebinned = new TH1D("H1D_Feeddown_denominator_rebinned_"+NamePart[ipart_option], "H1D_Feeddown_denominator_rebinned_"+NamePart[ipart_option], nbinpT_j, pTbinsXi);

  TH1I* H1I_SelectedEventCount = (TH1I*)file_O2Analysis->Get("lambdakzero-analysis-mc/hSelectedEventCounter");

  cout << "______Rebinning:" << endl;
  for(int i = 1; i <= nbinpT; i++){
    int arrayID_i = i - 1;
    int LowerBound_bin_i = H2D_Feeddown_numerator->GetXaxis()->FindBin(pTbins[arrayID_i]);
    int UpperBound_bin_i = H2D_Feeddown_numerator->GetXaxis()->FindBin(pTbins[arrayID_i+1]) - 1;
    cout << "LowerBound bin_i " << i << " = " << LowerBound_bin_i << ", UpperBound bin_i " << i << " = " << UpperBound_bin_i << endl;
    // for(int j = 1; j <= nbinpT; j++){
    for(int j = 1; j <= nbinpT_j; j++){
      int arrayID_j = j - 1;
      int LowerBound_bin_j = H2D_Feeddown_numerator->GetYaxis()->FindBin(pTbins_j[arrayID_j]);
      int UpperBound_bin_j = H2D_Feeddown_numerator->GetYaxis()->FindBin(pTbins_j[arrayID_j+1]) - 1;
      cout << "LowerBound bin_j " << j << " = " << LowerBound_bin_j << ", UpperBound bin_j " << j << " = " << UpperBound_bin_j << endl;
      H2D_Feeddown_numerator_rebinned->SetBinContent(i, j, H2D_Feeddown_numerator->Integral(LowerBound_bin_i,UpperBound_bin_i,LowerBound_bin_j,UpperBound_bin_j));
      H1D_Feeddown_denominator_rebinned->SetBinContent(j, H1D_Feeddown_denominator->Integral(LowerBound_bin_j,UpperBound_bin_j));
    }
  }

  
  H2D_Feeddown_numerator_rebinned->Sumw2();

  TH2D* hFij = new TH2D("hFij_"+NamePart[ipart_option], "hFij_"+NamePart[ipart_option], nbinpT, pTbins, nbinpT_j, pTbins_j);
  for(int i = 1; i <= nbinpT; i++){
    // for(int j = 1; j <= nbinpT; j++){
    for(int j = 1; j <= nbinpT_j; j++){
      hFij->SetBinContent(i,j,H2D_Feeddown_numerator_rebinned->GetBinContent(i,j)*1./H1D_Feeddown_denominator_rebinned->GetBinContent(j));
      hFij->SetBinError(i,j,sqrt(hFij->GetBinContent(i,j)*(1-hFij->GetBinContent(i,j))/H1D_Feeddown_denominator_rebinned->GetBinContent(j))); // this is effectively an efficiency: xsi- decays to (pi- + Lambda) 99.9% of the time (see PDG): so error is binomial: sqrt(np(1-p))
      // hFij->SetBinError(i,j,sqrt(hFij->GetBinContent(i,j)*(1-hFij->GetBinContent(i,j))*H1D_Feeddown_denominator_rebinned->GetBinContent(j))*1./H1D_Feeddown_denominator_rebinned->GetBinContent(j)); // OLD error calculation: 1/n*sqrt(np(1-p)); not sure why we multiply by 1/n
    }
  }

  //Xsi yield from pp 900GeV paper fit function (Xsi- + AntiXsi+)
  //dNdy = 0.0101 ± 0.0020
  //T = 175 ± 50
  //n = 5.2 ± 2.3
  float dNdy[2] = {3.05470e-04, 0.0020};
  float lTemp[2] = {1.00055e-01,50};
  float lPower[2] = {3.6,2.3};
  TF1* LevisTsalisFunction_Xsipp900GeV = new TF1("LevisTsalisFunction_Xsipp900GeV_"+NamePart[ipart_option], myLevyPtXsi, 0, pTbins_j[nbinpT_j], 3);

  // Get generated Xsi distribution
  // to get more statistics, get Xsi- + Xsi+ given the ratio Xsi-/Xsi+ is very close to 1, within a few percents; DAMN only true for Xsi-/Xsi0
  // TH1D* H1D_Xsi_PtDistrib;
  if (ipart_option > 2) {
    cout << "ipart_option used is not covered by Get_FeedDown()" << endl;
  }
  TH1D* H1D_Xsi_PtDistrib = (TH1D*)file_O2Analysis->Get("lambdakzero-particle-count-mc/h"+NamePart_multiStrange[ipart_option-1]+"Count_PtDiff");

  int XiDistrib_RebinFactor = 5.0;
  TH1F *H1D_Xsi_PtDistrib_rebinned = (TH1F*)H1D_Xsi_PtDistrib->Rebin(XiDistrib_RebinFactor,"H1D_Xsi_PtDistrib_rebinned");
  double dpT_Xi = H1D_Xsi_PtDistrib_rebinned->GetXaxis()->GetBinWidth(1);

  H1D_Xsi_PtDistrib_rebinned->Sumw2();
  // H1D_XsiMinus_PtDistrib->Scale(1./H1I_SelectedEventCount->GetEntries()); //scale to per event yield; done outside function because fit likes it better if not done before, for errors.

  LevisTsalisFunction_Xsipp900GeV->SetParameter(0, dNdy[0]*H1I_SelectedEventCount->GetEntries());
  LevisTsalisFunction_Xsipp900GeV->SetParameter(1, lTemp[0]);
  LevisTsalisFunction_Xsipp900GeV->SetParameter(2, lPower[0]);
  // LevisTsalisFunction_Xsipp900GeV->FixParameter(3,H1I_SelectedEventCount->GetEntries());

  // LevisTsalisFunction_Xsipp900GeV->SetParError(0, dNdy[1]);
  // LevisTsalisFunction_Xsipp900GeV->SetParError(1, lTemp[1]);
  // LevisTsalisFunction_Xsipp900GeV->SetParError(2, lPower[1]);

  // LevisTsalisFunction_Xsipp900GeV->SetParLimits(0, 0.0001, 0.1);
  // LevisTsalisFunction_Xsipp900GeV->SetParLimits(1, 0, 900);
  // LevisTsalisFunction_Xsipp900GeV->SetParLimits(2, 0, 20);



  // TFitResultPtr fFitResult = H1D_XsiMinus_PtDistrib->Fit(LevisTsalisFunction_Xsipp900GeV,"S+WLR"); //
  TFitResultPtr fFitResult = H1D_Xsi_PtDistrib_rebinned->Fit(LevisTsalisFunction_Xsipp900GeV,"S+L0R"); //

  TCanvas *canvasFeedDown = new TCanvas ("canvasFeedDown", NamehistoInvMass[ipart_option], 1200, 800);
  canvasFeedDown->cd(0); //ATTENTION: +1 because canvas divide starts counting from 1
  TH1 *hFrameFeeddown = canvasFeedDown->DrawFrame(0.,0,pTbins_j[nbinpT_j],1.2*H1D_Xsi_PtDistrib_rebinned->GetMaximum());
  hFrameFeeddown->SetXTitle("#it{p}_{T} (GeV/#it{c})");
  hFrameFeeddown->SetYTitle("d#it{N}/d#it{p}_{T}("+NamePart_multiStrange_Latex[ipart_option-1]+")");
  H1D_Xsi_PtDistrib_rebinned->Draw("E, same");
  LevisTsalisFunction_Xsipp900GeV->Draw("E, same");
  canvasFeedDown->SaveAs("FeedDownLevyFit"+NamePart_multiStrange_Latex[ipart_option-1]+".pdf");

  double * param = LevisTsalisFunction_Xsipp900GeV->GetParameters();
  TMatrixDSym covMatrix = fFitResult->GetCovarianceMatrix();

  double FeedDownCorrection_i_Count = 0;
  double FeedDownCorrection_i_Count_StatUncertainty = 0;
  for(int i = 1; i <= nbinpT; i++){
    FeedDownCorrection_i_Count = 0;
    FeedDownCorrection_i_Count_StatUncertainty = 0;
    cout << "*---------------- i: " << i << " ----------------*" << endl;
    for(int j = 1; j <= nbinpT_j; j++){
      // int pT_Xi_lowEdge = H1D_Xsi_PtDistrib_rebinned->GetXaxis()->GetBinLowEdge(j);
      // int pT_Xi_upEdge = H1D_Xsi_PtDistrib_rebinned->GetXaxis()->GetBinLowEdge(j)+dpT_Xi;
      double FittedYield_Xi_pTbin_j = LevisTsalisFunction_Xsipp900GeV->Integral(pTbins_j[j-1],pTbins_j[j])*1./dpT_Xi; 
      // double FittedYield_Xi_pTbin_j = LevisTsalisFunction_Xsipp900GeV->Integral(pTbins[j-1],pTbins[j])*1./((pTbins[j]-pTbins[j-1])); 

      // double FittedYield_Xi_pTbin_j = H1D_Xsi_PtDistrib_rebinned->Integral(pTbins[j-1],pTbins[j]); // temporary
      cout << "j: " << j << ", FittedYield_Xi_pTbin_j = " << FittedYield_Xi_pTbin_j << endl; //surprisingly is most high at j=9 with 400  -> because pT bins increase in length

      double FittedYield_Xi_pTbin_j_uncertainty = 1./dpT_Xi * LevisTsalisFunction_Xsipp900GeV->IntegralError(pTbins_j[j-1],pTbins_j[j], param, covMatrix.GetMatrixArray());
      // double FittedYield_Xi_pTbin_j_uncertainty = 0; test

      // double MatrixTimesYield_value = hFij->GetBinContent(i,j)*FittedYield_Xi_pTbin_j;
      // double MatrixTimesYield_error_Square = 0;
      // if (abs(MatrixTimesYield_value) < 1E-9){ //avoid division by 0; sets to 0 if MatrixTimesYield_value = 0
      // }
      // else {
      //   MatrixTimesYield_error_Square = MatrixTimesYield_value * ( pow(hFij->GetBinError(i,j)/hFij->GetBinContent(i,j),2) + pow((FittedYield_Xi_pTbin_j_uncertainty)/(FittedYield_Xi_pTbin_j),2) );
      // }
      double MatrixTimesYield_error_Square = pow(FittedYield_Xi_pTbin_j*hFij->GetBinError(i,j) + pow(hFij->GetBinContent(i,j),2)*FittedYield_Xi_pTbin_j_uncertainty,2); //  Error propagation for product: (Sigma_AB)^2 = B^2 * Sigma_A + A^2 * Sigma_B

      FeedDownCorrection_i_Count = FeedDownCorrection_i_Count + hFij->GetBinContent(i,j)*FittedYield_Xi_pTbin_j;
      FeedDownCorrection_i_Count_StatUncertainty = sqrt(pow(FeedDownCorrection_i_Count_StatUncertainty,2) + MatrixTimesYield_error_Square); //product and sums
      cout << "      numerator = " << H2D_Feeddown_numerator_rebinned->GetBinContent(i,j) << endl;
      cout << "      denominator = " << H1D_Feeddown_denominator_rebinned->GetBinContent(j) << endl;
      cout << "      hFij = " << hFij->GetBinContent(i,j) << endl;
      cout << "      hFij Error = " << hFij->GetBinError(i,j) << endl;
      cout << "      FeedDownCorrection_i_Count = " << FeedDownCorrection_i_Count << endl;
      cout << "      FeedDownCorrection_i_Count_StatUncertainty = " << FeedDownCorrection_i_Count_StatUncertainty << endl;
    }

    double dpT = H1D_FeedDownCorrection->GetXaxis()->GetBinWidth(i);
    double drapidity = 1.5;
    double Feeddown_d2NdpTdy =  FeedDownCorrection_i_Count *1./dpT*1./drapidity; //from count to yield d2N/dpTdy
    double Feeddown_d2NdpTdy_StatUncertainty =  FeedDownCorrection_i_Count_StatUncertainty * 1./dpT*1./drapidity; //from count to yield d2N/dpTdy
    H1D_FeedDownCorrection->SetBinContent(i,Feeddown_d2NdpTdy);
    H1D_FeedDownCorrection->SetBinError(i, Feeddown_d2NdpTdy_StatUncertainty);
    cout << "           Feeddown_d2NdpTdy_i = " << Feeddown_d2NdpTdy << endl;
    cout << "                          error = " << Feeddown_d2NdpTdy_StatUncertainty << endl;
  }

  cout << "DEBUG GetFeeddown end" << endl;
  // H1D_FeedDownCorrection->Scale(2.); // need to multiply contribution by two given Xsi0 is also contributing to Lambda feeddown and has not been put in the Fij matrix; this is done based on the assumption, verified to 2% according to analysis note, that yield of Xi0 is about the same as the one of Xi- and the same as the one of Xi+
  //                                    ideally we could do the same with Xsi0 as done here for Xsi-/+ but Xsi0 can't be measured given it decays into pi0 which is neutral (and Lambda obviously that we measure but isn't enough on its own to reconstruct the decay and that polutes our primary lambda result)
  // ACTUALLY we do put xsi0 in the fij (see lambdaanalysisMC code)

  H1D_FeedDownCorrection->Scale(1./H1I_SelectedEventCount->GetEntries()); //scale to per event yield; 
  delete hFij;
  delete H2D_Feeddown_numerator_rebinned;
  canvasFeedDown->Close();
  delete LevisTsalisFunction_Xsipp900GeV;
}


void Get_FeedDown_RawCount(TH1D * H1D_FeedDownCorrection, int ipart_option, TFile *file_O2Analysis, double* pTbins, int nbinpT) {
  TH2D* H2D_Feeddown_numerator = (TH2D*)file_O2Analysis->Get("lambdakzero-analysis-mc/h"+NamePart[ipart_option]+"FeedDownMatrix");
  TH1D* H1D_Feeddown_denominator = (TH1D*)file_O2Analysis->Get("lambdakzero-particle-count-mc/h"+NamePart_multiStrange[ipart_option-1]+"Count_PtDiff");
  // Rebinning of those two histograms
  // TH1D* H1D_Feeddown_denominator_rebinned = (TH1D*)H1D_Feeddown_denominator->Rebin(nbinpT,"H1D_Feeddown_denominator_rebinned_"+NamePart_multiStrange[ipart_option],pTbins); could do it that way for 1 dim hist but we decide to do it the same way as for the 2d hist here

  // int nbinpT_j = nbinpT;
  // double* pTbins_j = pTbins;
  int nbinpT_j = nbinpT_Xi;
  double* pTbins_j = pTbinsXi;

  TH2D* H2D_Feeddown_numerator_rebinned = new TH2D("H2D_Feeddown_numerator_rebinned_"+NamePart[ipart_option], "H2D_Feeddown_numerator_rebinned_"+NamePart[ipart_option], nbinpT, pTbins, nbinpT_j, pTbinsXi);
  TH1D* H1D_Feeddown_denominator_rebinned = new TH1D("H1D_Feeddown_denominator_rebinned_"+NamePart[ipart_option], "H1D_Feeddown_denominator_rebinned_"+NamePart[ipart_option], nbinpT_j, pTbinsXi);

  TH1I* H1I_SelectedEventCount = (TH1I*)file_O2Analysis->Get("lambdakzero-analysis-mc/hSelectedEventCounter");

  cout << "______Rebinning:" << endl;
  for(int i = 1; i <= nbinpT; i++){
    int arrayID_i = i - 1;
    int LowerBound_bin_i = H2D_Feeddown_numerator->GetXaxis()->FindBin(pTbins[arrayID_i]);
    int UpperBound_bin_i = H2D_Feeddown_numerator->GetXaxis()->FindBin(pTbins[arrayID_i+1]) - 1;
    cout << "LowerBound bin_i " << i << " = " << LowerBound_bin_i << ", UpperBound bin_i " << i << " = " << UpperBound_bin_i << endl;
    // for(int j = 1; j <= nbinpT; j++){
    for(int j = 1; j <= nbinpT_j; j++){
      int arrayID_j = j - 1;
      int LowerBound_bin_j = H2D_Feeddown_numerator->GetYaxis()->FindBin(pTbins_j[arrayID_j]);
      int UpperBound_bin_j = H2D_Feeddown_numerator->GetYaxis()->FindBin(pTbins_j[arrayID_j+1]) - 1;
      cout << "LowerBound bin_j " << j << " = " << LowerBound_bin_j << ", UpperBound bin_j " << j << " = " << UpperBound_bin_j << endl;
      H2D_Feeddown_numerator_rebinned->SetBinContent(i, j, H2D_Feeddown_numerator->Integral(LowerBound_bin_i,UpperBound_bin_i,LowerBound_bin_j,UpperBound_bin_j));
      H1D_Feeddown_denominator_rebinned->SetBinContent(j, H1D_Feeddown_denominator->Integral(LowerBound_bin_j,UpperBound_bin_j));
    }
  }

  
  H2D_Feeddown_numerator_rebinned->Sumw2();

  TH2D* hFij = new TH2D("hFij_"+NamePart[ipart_option], "hFij_"+NamePart[ipart_option], nbinpT, pTbins, nbinpT_j, pTbins_j);
  for(int i = 1; i <= nbinpT; i++){
    // for(int j = 1; j <= nbinpT; j++){
    for(int j = 1; j <= nbinpT_j; j++){
      hFij->SetBinContent(i,j,H2D_Feeddown_numerator_rebinned->GetBinContent(i,j)*1./H1D_Feeddown_denominator_rebinned->GetBinContent(j));
      hFij->SetBinError(i,j,sqrt(hFij->GetBinContent(i,j)*(1-hFij->GetBinContent(i,j))/H1D_Feeddown_denominator_rebinned->GetBinContent(j))); // this is effectively an efficiency: xsi- decays to (pi- + Lambda) 99.9% of the time (see PDG): so error is binomial: sqrt(np(1-p))
      // hFij->SetBinError(i,j,sqrt(hFij->GetBinContent(i,j)*(1-hFij->GetBinContent(i,j))*H1D_Feeddown_denominator_rebinned->GetBinContent(j))*1./H1D_Feeddown_denominator_rebinned->GetBinContent(j)); // OLD error calculation: 1/n*sqrt(np(1-p)); not sure why we multiply by 1/n
    }
  }

  //Xsi yield from pp 900GeV paper fit function (Xsi- + AntiXsi+)
  //dNdy = 0.0101 ± 0.0020
  //T = 175 ± 50
  //n = 5.2 ± 2.3
  float dNdy[2] = {3.05470e-04, 0.0020};
  float lTemp[2] = {1.00055e-01,50};
  float lPower[2] = {3.6,2.3};
  TF1* LevisTsalisFunction_Xsipp900GeV = new TF1("LevisTsalisFunction_Xsipp900GeV_"+NamePart[ipart_option], myLevyPtXsi, 0, pTbins_j[nbinpT_j], 3);

  // Get generated Xsi distribution
  // to get more statistics, get Xsi- + Xsi+ given the ratio Xsi-/Xsi+ is very close to 1, within a few percents; DAMN only true for Xsi-/Xsi0
  // TH1D* H1D_Xsi_PtDistrib;
  if (ipart_option > 2) {
    cout << "ipart_option used is not covered by Get_FeedDown()" << endl;
  }
  TH1D* H1D_Xsi_PtDistrib = (TH1D*)file_O2Analysis->Get("lambdakzero-particle-count-mc/h"+NamePart_multiStrange[ipart_option-1]+"Count_PtDiff");

  int XiDistrib_RebinFactor = 5.0;
  TH1F *H1D_Xsi_PtDistrib_rebinned = (TH1F*)H1D_Xsi_PtDistrib->Rebin(XiDistrib_RebinFactor,"H1D_Xsi_PtDistrib_rebinned");
  double dpT_Xi = H1D_Xsi_PtDistrib_rebinned->GetXaxis()->GetBinWidth(1);

  H1D_Xsi_PtDistrib_rebinned->Sumw2();
  // H1D_XsiMinus_PtDistrib->Scale(1./H1I_SelectedEventCount->GetEntries()); //scale to per event yield; done outside function because fit likes it better if not done before, for errors.

  LevisTsalisFunction_Xsipp900GeV->SetParameter(0, dNdy[0]*H1I_SelectedEventCount->GetEntries());
  LevisTsalisFunction_Xsipp900GeV->SetParameter(1, lTemp[0]);
  LevisTsalisFunction_Xsipp900GeV->SetParameter(2, lPower[0]);
  // LevisTsalisFunction_Xsipp900GeV->FixParameter(3,H1I_SelectedEventCount->GetEntries());

  // LevisTsalisFunction_Xsipp900GeV->SetParError(0, dNdy[1]);
  // LevisTsalisFunction_Xsipp900GeV->SetParError(1, lTemp[1]);
  // LevisTsalisFunction_Xsipp900GeV->SetParError(2, lPower[1]);

  // LevisTsalisFunction_Xsipp900GeV->SetParLimits(0, 0.0001, 0.1);
  // LevisTsalisFunction_Xsipp900GeV->SetParLimits(1, 0, 900);
  // LevisTsalisFunction_Xsipp900GeV->SetParLimits(2, 0, 20);



  // TFitResultPtr fFitResult = H1D_XsiMinus_PtDistrib->Fit(LevisTsalisFunction_Xsipp900GeV,"S+WLR"); //
  TFitResultPtr fFitResult = H1D_Xsi_PtDistrib_rebinned->Fit(LevisTsalisFunction_Xsipp900GeV,"S+L0R"); //

  TCanvas *canvasFeedDown = new TCanvas ("canvasFeedDown", NamehistoInvMass[ipart_option], 1200, 800);
  canvasFeedDown->cd(0); //ATTENTION: +1 because canvas divide starts counting from 1
  TH1 *hFrameFeeddown = canvasFeedDown->DrawFrame(0.,0,pTbins_j[nbinpT_j],1.2*H1D_Xsi_PtDistrib_rebinned->GetMaximum());
  hFrameFeeddown->SetXTitle("#it{p}_{T} (GeV/#it{c})");
  hFrameFeeddown->SetYTitle("d#it{N}/d#it{p}_{T}("+NamePart_multiStrange_Latex[ipart_option-1]+")");
  H1D_Xsi_PtDistrib_rebinned->Draw("E, same");
  LevisTsalisFunction_Xsipp900GeV->Draw("E, same");
  canvasFeedDown->SaveAs("FeedDownLevyFit"+NamePart_multiStrange_Latex[ipart_option-1]+".pdf");

  double * param = LevisTsalisFunction_Xsipp900GeV->GetParameters();
  TMatrixDSym covMatrix = fFitResult->GetCovarianceMatrix();

  double FeedDownCorrection_i_Count = 0;
  double FeedDownCorrection_i_Count_StatUncertainty = 0;
  for(int i = 1; i <= nbinpT; i++){
    FeedDownCorrection_i_Count = 0;
    FeedDownCorrection_i_Count_StatUncertainty = 0;
    cout << "*---------------- i: " << i << " ----------------*" << endl;
    for(int j = 1; j <= nbinpT_j; j++){
      // int pT_Xi_lowEdge = H1D_Xsi_PtDistrib_rebinned->GetXaxis()->GetBinLowEdge(j);
      // int pT_Xi_upEdge = H1D_Xsi_PtDistrib_rebinned->GetXaxis()->GetBinLowEdge(j)+dpT_Xi;
      double FittedYield_Xi_pTbin_j = LevisTsalisFunction_Xsipp900GeV->Integral(pTbins_j[j-1],pTbins_j[j])*1./dpT_Xi; 
      // double FittedYield_Xi_pTbin_j = LevisTsalisFunction_Xsipp900GeV->Integral(pTbins[j-1],pTbins[j])*1./((pTbins[j]-pTbins[j-1])); 

      // double FittedYield_Xi_pTbin_j = H1D_Xsi_PtDistrib_rebinned->Integral(pTbins[j-1],pTbins[j]); // temporary
      cout << "j: " << j << ", FittedYield_Xi_pTbin_j = " << FittedYield_Xi_pTbin_j << endl; //surprisingly is most high at j=9 with 400  -> because pT bins increase in length

      double FittedYield_Xi_pTbin_j_uncertainty = 1./dpT_Xi * LevisTsalisFunction_Xsipp900GeV->IntegralError(pTbins_j[j-1],pTbins_j[j], param, covMatrix.GetMatrixArray());
      // double FittedYield_Xi_pTbin_j_uncertainty = 0; test

      // double MatrixTimesYield_value = hFij->GetBinContent(i,j)*FittedYield_Xi_pTbin_j;
      // double MatrixTimesYield_error_Square = 0;
      // if (abs(MatrixTimesYield_value) < 1E-9){ //avoid division by 0; sets to 0 if MatrixTimesYield_value = 0
      // }
      // else {
      //   MatrixTimesYield_error_Square = MatrixTimesYield_value * ( pow(hFij->GetBinError(i,j)/hFij->GetBinContent(i,j),2) + pow((FittedYield_Xi_pTbin_j_uncertainty)/(FittedYield_Xi_pTbin_j),2) );
      // }
      double MatrixTimesYield_error_Square = pow(FittedYield_Xi_pTbin_j*hFij->GetBinError(i,j) + pow(hFij->GetBinContent(i,j),2)*FittedYield_Xi_pTbin_j_uncertainty,2); //  Error propagation for product: (Sigma_AB)^2 = B^2 * Sigma_A + A^2 * Sigma_B

      FeedDownCorrection_i_Count = FeedDownCorrection_i_Count + hFij->GetBinContent(i,j)*FittedYield_Xi_pTbin_j;
      FeedDownCorrection_i_Count_StatUncertainty = sqrt(pow(FeedDownCorrection_i_Count_StatUncertainty,2) + MatrixTimesYield_error_Square); //product and sums
      cout << "      numerator = " << H2D_Feeddown_numerator_rebinned->GetBinContent(i,j) << endl;
      cout << "      denominator = " << H1D_Feeddown_denominator_rebinned->GetBinContent(j) << endl;
      cout << "      hFij = " << hFij->GetBinContent(i,j) << endl;
      cout << "      hFij Error = " << hFij->GetBinError(i,j) << endl;
      cout << "      FeedDownCorrection_i_Count = " << FeedDownCorrection_i_Count << endl;
      cout << "      FeedDownCorrection_i_Count_StatUncertainty = " << FeedDownCorrection_i_Count_StatUncertainty << endl;
    }

    double Feeddown_Count =  FeedDownCorrection_i_Count; //
    double Feeddown_Count_StatUncertainty =  FeedDownCorrection_i_Count_StatUncertainty; //
    H1D_FeedDownCorrection->SetBinContent(i,Feeddown_Count);
    H1D_FeedDownCorrection->SetBinError(i, Feeddown_Count_StatUncertainty);
    cout << "           Feeddown_Count_i = " << Feeddown_Count << endl;
    cout << "                          error = " << Feeddown_Count_StatUncertainty << endl;
  }

  // H1D_FeedDownCorrection->Scale(2.); // need to multiply contribution by two given Xsi0 is also contributing to Lambda feeddown and has not been put in the Fij matrix; this is done based on the assumption, verified to 2% according to analysis note, that yield of Xi0 is about the same as the one of Xi- and the same as the one of Xi+
  //                                    ideally we could do the same with Xsi0 as done here for Xsi-/+ but Xsi0 can't be measured given it decays into pi0 which is neutral (and Lambda obviously that we measure but isn't enough on its own to reconstruct the decay and that polutes our primary lambda result)
  // ACTUALLY we do put xsi0 in the fij (see lambdaanalysisMC code)

  // H1D_FeedDownCorrection->Scale(1./H1I_SelectedEventCount->GetEntries()); //scale to per event yield; 
  delete hFij;
  delete H2D_Feeddown_numerator_rebinned;
  canvasFeedDown->Close();
  delete LevisTsalisFunction_Xsipp900GeV;
}

void Draw_FeedDownMatrix(TH1D* &hstat, TH1D* &hsyst, TH1D* &hsystCorr, int ipart_option, TFile *file_O2Analysis, TString* SaveAs_Title, TString* &texXtitle, TString* &texYtitle, double* pTbins, int nbinpT) {

  TH2D* H2D_Feeddown_numerator = (TH2D*)file_O2Analysis->Get("lambdakzero-analysis-mc/h"+NamePart[ipart_option]+"FeedDownMatrix");
  TH1D* H1D_Feeddown_denominator = (TH1D*)file_O2Analysis->Get("lambdakzero-particle-count-mc/h"+NamePart_multiStrange[ipart_option-1]+"Count_PtDiff");

  double pTbins_start0[nbinpT+2];
  pTbins_start0[0] = 0.;
  for(int i = 1; i <= nbinpT+1; i++){
    pTbins_start0[i] = pTbins[i-1];
    cout << pTbins_start0[i] << endl;
  }

  TH2D* H2Frame = new TH2D("H2Frame", "H2Frame", nbinpT+1, pTbins_start0, nbinpT+1, pTbins_start0);
  TH2D* H2D_Feeddown_numerator_rebinned = new TH2D("H2D_Feeddown_numerator_rebinned", "H2D_Feeddown_numerator_rebinned", nbinpT+1, pTbins_start0, nbinpT+1, pTbins_start0);
  TH1D* H1D_Feeddown_denominator_rebinned = new TH1D("H1D_Feeddown_denominator_rebinned_"+NamePart[ipart_option], "H1D_Feeddown_denominator_rebinned_"+NamePart[ipart_option], nbinpT, pTbins);

  for(int i = 1; i <= nbinpT; i++){
    int arrayID_i = i - 1;
    int LowerBound_bin_i = H2D_Feeddown_numerator->GetXaxis()->FindBin(pTbins[arrayID_i]);
    int UpperBound_bin_i = H2D_Feeddown_numerator->GetXaxis()->FindBin(pTbins[arrayID_i+1]) - 1;
    for(int j = 1; j <= nbinpT; j++){
      int arrayID_j = j - 1;
      int LowerBound_bin_j = H2D_Feeddown_numerator->GetYaxis()->FindBin(pTbins[arrayID_j]);
      int UpperBound_bin_j = H2D_Feeddown_numerator->GetYaxis()->FindBin(pTbins[arrayID_j+1]) - 1;
      H2D_Feeddown_numerator_rebinned->SetBinContent(i+1, j+1, H2D_Feeddown_numerator->Integral(LowerBound_bin_i,UpperBound_bin_i,LowerBound_bin_j,UpperBound_bin_j));
      H1D_Feeddown_denominator_rebinned->SetBinContent(j, H1D_Feeddown_denominator->Integral(LowerBound_bin_j,UpperBound_bin_j));
    }
  }

  TH2D* hFij = new TH2D("hFij_"+NamePart[ipart_option], "hFij_"+NamePart[ipart_option], nbinpT, pTbins, nbinpT, pTbins);
  for(int i = 1; i <= nbinpT; i++){
    for(int j = 1; j <= nbinpT; j++){
      hFij->SetBinContent(i,j,H2D_Feeddown_numerator_rebinned->GetBinContent(i,j)*1./H1D_Feeddown_denominator_rebinned->GetBinContent(j));
      hFij->SetBinError(i,j,sqrt(hFij->GetBinContent(i,j)*(1-hFij->GetBinContent(i,j))/H1D_Feeddown_denominator_rebinned->GetBinContent(j))); // this is effectively an efficiency: xsi- decays to (pi- + Lambda) 99.9% of the time (see PDG): so error is binomial: sqrt(np(1-p))
      // hFij->SetBinError(i,j,sqrt(hFij->GetBinContent(i,j)*(1-hFij->GetBinContent(i,j))*H1D_Feeddown_denominator_rebinned->GetBinContent(j))*1./H1D_Feeddown_denominator_rebinned->GetBinContent(j)); // OLD error calculation: 1/n*sqrt(np(1-p)); not sure why we multiply by 1/n
    }
  }

  TCanvas *canvasFeedDown = new TCanvas ("canvasFeedDown", NamehistoInvMass[ipart_option], 800, 800);
  canvasFeedDown->cd(0); //ATTENTION: +1 because canvas divide starts counting from 1
  // TH1 *hFrame = canvasFeedDown->DrawFrame(pTbins[0],0,pTbins[nbinpT],1.2*hSystematicUncertainty_Total->GetMaximum());
  // hFrame->SetXTitle("#it{p}_{T} (GeV/#it{c})");
  // hFrame->SetYTitle("uncertainty (ratio)");
  hFij->SetXTitle(texPt_particle[ipart_option]->Data());
  hFij->SetYTitle(texPt_Xi->Data());
  hFij->Draw("COLZ");
  canvasFeedDown->SetRightMargin(0.15);
  canvasFeedDown->SaveAs("FeedDown_Matrix_"+NamePart[ipart_option]+".pdf");

  /////////////////////////////////Does Nothing - only there to respect template/////////////////////////////////////////////
  hstat = (TH1D*)H2D_Feeddown_numerator_rebinned->Clone("hstat");
  hstat->Reset("M");
  // hstat->Divide(hPtDiff_RawSpectrum_O2,hPtDiff_RawSpectrum_HEP[ipart_option]);

  //////Error Bars///////
  hsyst     = new TH1D("hsyst", "hsyst", 100, 0, 10);
  hsystCorr = new TH1D("hsystCorr", "hsystCorr", 100, 0, 10);
  hstat     ->Sumw2(); 
  hsyst     ->Sumw2();
  hsystCorr ->Sumw2(); 
  int nbinx = hstat->GetNbinsX();
  
  for(int ibinx = 1; ibinx <= nbinx; ibinx++){
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


void Draw_XsiPlus_MassPlot(TH1D* &hstat, TH1D* &hsyst, TH1D* &hsystCorr, int ipart, TFile *file_O2Analysis, TString* SaveAs_Title, TString* &texXtitle, TString* &texYtitle, double* pTbins, int nbinpT) {

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
  int nbinx = hstat->GetNbinsX();
  
  for(int ibinx = 1; ibinx <= nbinx; ibinx++){
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

void Draw_FeedDown_ipart(TH1D* &hstat, TH1D* &hsyst, TH1D* &hsystCorr, int ipart, TFile *file_O2Analysis, TString* SaveAs_Title, TString* &texXtitle, TString* &texYtitle, double* pTbins, int nbinpT) {

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
  int InvMassRebinFactor = InvMassRebinFactor_standard[ipart]; //used if signalextractiontype = 1
  double SignalMean, SignalStandardDev;
  Get_RawYield_vsPt(hRawYield_vsPt, H3D_detectedV0s_DataLike, H1I_SelectedEventCount, file_O2Analysis, &SignalMean, &SignalStandardDev, ipart, pTbins, nbinpT, 0, -1, 0, SideBandSizeMultiplierModifier_array[0], SignalExtractionType_default, InvMassRebinFactor); ///real line that need to take the place of the one above once it's been written

  // ratio of feeddown removed to raw yield
  TH1D* H1D_FractionRawYieldRemoved = new TH1D("H1D_FractionRawYieldRemoved", "H1D_FractionRawYieldRemoved", nbinpT, pTbins);
  H1D_FractionRawYieldRemoved->Divide(H1D_Feeddown_Correction,hRawYield_vsPt);

  TCanvas *canvasFeedDown = new TCanvas ("canvasFeedDown", NamehistoInvMass[ipart], 1200, 800);
  canvasFeedDown->cd(0); //ATTENTION: +1 because canvas divide starts counting from 1
  TH1 *hFrame = canvasFeedDown->DrawFrame(0.,0,pTbins[nbinpT],1.2*H1D_FractionRawYieldRemoved->GetMaximum());
  hFrame->SetXTitle("#it{p}_{T} (GeV/#it{c})");
  hFrame->SetYTitle("FractionRawYieldRemoved from feeddown");
  H1D_FractionRawYieldRemoved->Draw("E,same");
  H1D_FractionRawYieldRemoved->SetMarkerStyle(markers[0]);
  H1D_FractionRawYieldRemoved->SetMarkerColor(colors[0]);
  H1D_FractionRawYieldRemoved->SetLineColor(colors[0]);
  H1D_FractionRawYieldRemoved->SetFillStyle(0); // To draw empty boxes
  canvasFeedDown->SaveAs("FractionRawYieldRemoved.pdf");

  /////////////////////////////////Does Nothing - only there to respect template/////////////////////////////////////////////
  hstat = (TH1D*)H1D_FractionRawYieldRemoved->Clone("hstat");
  hstat->Reset("M");
  // hstat->Divide(hPtDiff_RawSpectrum_O2,hPtDiff_RawSpectrum_HEP[ipart]);

  //////Error Bars///////
  hsyst     = new TH1D("hsyst", "hsyst", 100, 0, 10);
  hsystCorr = new TH1D("hsystCorr", "hsystCorr", 100, 0, 10);
  hstat     ->Sumw2(); 
  hsyst     ->Sumw2();
  hsystCorr ->Sumw2(); 
  int nbinx = hstat->GetNbinsX();
  
  for(int ibinx = 1; ibinx <= nbinx; ibinx++){
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


void Draw_FeedDown_LambdaAntiLambda(TH1D* &hstat, TH1D* &hsyst, TH1D* &hsystCorr, int ipart, TFile *file_O2Analysis, TString* SaveAs_Title, TString* &texXtitle, TString* &texYtitle, double* pTbins, int nbinpT) {

  int ipart_Lambda = 1;
  int ipart_AntiLambda = 2;

  TH1D* H1D_XsiPlus_InvMassPlot = (TH1D*)file_O2Analysis->Get("lambdakzero-particle-count-mc/hXsiPlusCount_PtDiff");
  TH1D* H1D_XsiPlus_InvMassPlot_rebinned = (TH1D*)H1D_XsiPlus_InvMassPlot->Rebin(nbinpT,"H1D_XsiPlus_InvMassPlot_rebinned",pTbins);

  //Event count
  TH1I* H1I_SelectedEventCount = (TH1I*)file_O2Analysis->Get("lambdakzero-analysis-mc/hSelectedEventCounter");

  //Feeddown
  TH1D* H1D_Feeddown_Correction_Lambda = new TH1D("H1D_Feeddown_Correction_Lambda", "H1D_Feeddown_Correction_Lambda", nbinpT, pTbins);
  TH1D* H1D_Feeddown_Correction_AntiLambda = new TH1D("H1D_Feeddown_Correction_AntiLambda", "H1D_Feeddown_Correction_AntiLambda", nbinpT, pTbins);

  cout << "--------------------------- LAMBDA Get_Feeddown ---------------------------" << endl;
  Get_FeedDown(H1D_Feeddown_Correction_Lambda, ipart_Lambda, file_O2Analysis, pTbins, nbinpT);
  cout << "--------------------------- antiLAMBDA Get_Feeddown ---------------------------" << endl;
  Get_FeedDown(H1D_Feeddown_Correction_AntiLambda, ipart_AntiLambda, file_O2Analysis, pTbins, nbinpT);
  // Get_FeedDown(H1D_Feeddown_Correction_Lambda, ipart_Lambda, file_O2Analysis, pTbins, nbinpT);
  // Get_FeedDown(H1D_Feeddown_Correction_AntiLambda, ipart_Lambda, file_O2Analysis, pTbins, nbinpT);

  //Raw yield
  TH3D* H3D_detectedV0s_DataLike_Lambda = (TH3D*)file_O2Analysis->Get("lambdakzero-analysis-mc/h3dMass"+NamePart[ipart_Lambda]);
  TH3D* H3D_detectedV0s_DataLike_AntiLambda = (TH3D*)file_O2Analysis->Get("lambdakzero-analysis-mc/h3dMass"+NamePart[ipart_AntiLambda]);
  TH1D *hRawYield_vsPt_Lambda;
  TH1D *hRawYield_vsPt_AntiLambda;
  int InvMassRebinFactor_Lambda = InvMassRebinFactor_standard[ipart_Lambda]; //used if signalextractiontype = 1
  int InvMassRebinFactor_AntiLambda = InvMassRebinFactor_standard[ipart_AntiLambda]; //used if signalextractiontype = 1
  double SignalMean_Lambda, SignalStandardDev_Lambda;
  double SignalMean_AntiLambda, SignalStandardDev_AntiLambda;
  Get_RawYield_vsPt(hRawYield_vsPt_Lambda, H3D_detectedV0s_DataLike_Lambda, H1I_SelectedEventCount, file_O2Analysis, &SignalMean_Lambda, &SignalStandardDev_Lambda, ipart_Lambda, pTbins, nbinpT, 0, -1, 0, SideBandSizeMultiplierModifier_array[0], SignalExtractionType_default, InvMassRebinFactor_Lambda); ///real line that need to take the place of the one above once it's been written
  Get_RawYield_vsPt(hRawYield_vsPt_AntiLambda, H3D_detectedV0s_DataLike_AntiLambda, H1I_SelectedEventCount, file_O2Analysis, &SignalMean_AntiLambda, &SignalStandardDev_AntiLambda, ipart_AntiLambda, pTbins, nbinpT, 0, -1, 0, SideBandSizeMultiplierModifier_array[0], SignalExtractionType_default, InvMassRebinFactor_AntiLambda); ///real line that need to take the place of the one above once it's been written

  cout << "Lambda RawYield = "<< hRawYield_vsPt_Lambda->GetBinContent(1) << ", FeeddownRemoved = " << H1D_Feeddown_Correction_Lambda->GetBinContent(1) << endl;
  cout << "antiLambda RawYield = "<< hRawYield_vsPt_AntiLambda->GetBinContent(1) << ", FeeddownRemoved = " << H1D_Feeddown_Correction_AntiLambda->GetBinContent(1) << endl;
  // ratio of feeddown removed to raw yield
  TH1D* H1D_FractionRawYieldRemoved_Lambda = new TH1D("H1D_FractionRawYieldRemoved_Lambda", "H1D_FractionRawYieldRemoved_Lambda", nbinpT, pTbins);
  H1D_FractionRawYieldRemoved_Lambda->Divide(H1D_Feeddown_Correction_Lambda,hRawYield_vsPt_Lambda);
  TH1D* H1D_FractionRawYieldRemoved_AntiLambda = new TH1D("H1D_FractionRawYieldRemoved_AntiLambda", "H1D_FractionRawYieldRemoved_AntiLambda", nbinpT, pTbins);
  H1D_FractionRawYieldRemoved_AntiLambda->Divide(H1D_Feeddown_Correction_AntiLambda,hRawYield_vsPt_AntiLambda);
  
  TCanvas *canvasFeedDown = new TCanvas ("canvasFeedDown", NamehistoInvMass[ipart], 1200, 800);
  canvasFeedDown->cd(0); //ATTENTION: +1 because canvas divide starts counting from 1

  float MaximumFraction = (H1D_FractionRawYieldRemoved_Lambda->GetMaximum() > H1D_FractionRawYieldRemoved_AntiLambda->GetMaximum()) ? MaximumFraction = H1D_FractionRawYieldRemoved_Lambda->GetMaximum() : MaximumFraction = H1D_FractionRawYieldRemoved_AntiLambda->GetMaximum();
  TH1 *hFrame = canvasFeedDown->DrawFrame(0.,0,pTbins[nbinpT],1.4*MaximumFraction);

  TLegend * leg = new TLegend(0.75, 0.74, 0.88, 0.87);

  hFrame->SetXTitle("#it{p}_{T} (GeV/#it{c})");
  hFrame->SetYTitle("ratio removed");

  H1D_FractionRawYieldRemoved_Lambda->Draw("E,same");
  H1D_FractionRawYieldRemoved_Lambda->SetMarkerStyle(markers[1]);
  H1D_FractionRawYieldRemoved_Lambda->SetMarkerColor(colors[0]);
  H1D_FractionRawYieldRemoved_Lambda->SetLineColor(colors[0]);
  H1D_FractionRawYieldRemoved_Lambda->SetFillStyle(1); // To draw empty boxes
  leg->AddEntry(H1D_FractionRawYieldRemoved_Lambda, "#Lambda",   "LP");

  H1D_FractionRawYieldRemoved_AntiLambda->Draw("E,same");
  H1D_FractionRawYieldRemoved_AntiLambda->SetMarkerStyle(markers[2]);
  H1D_FractionRawYieldRemoved_AntiLambda->SetMarkerColor(colors[0]);
  H1D_FractionRawYieldRemoved_AntiLambda->SetLineColor(colors[0]);
  H1D_FractionRawYieldRemoved_AntiLambda->SetFillStyle(0); // To draw empty boxes
  leg->AddEntry(H1D_FractionRawYieldRemoved_AntiLambda, "#bar{#Lambda}",   "LP");

  leg->SetTextSize(gStyle->GetTextSize()*1.);    
  leg->Draw("same");

  canvasFeedDown->SaveAs("FractionRawYieldRemoved_LambdaVsAntiLambda.pdf");

  /////////////////////////////////Does Nothing - only there to respect template/////////////////////////////////////////////
  hstat = (TH1D*)H1D_FractionRawYieldRemoved_AntiLambda->Clone("hstat");
  hstat->Reset("M");
  // hstat->Divide(hPtDiff_RawSpectrum_O2,hPtDiff_RawSpectrum_HEP[ipart]);

  //////Error Bars///////
  hsyst     = new TH1D("hsyst", "hsyst", 100, 0, 10);
  hsystCorr = new TH1D("hsystCorr", "hsystCorr", 100, 0, 10);
  hstat     ->Sumw2(); 
  hsyst     ->Sumw2();
  hsystCorr ->Sumw2(); 
  int nbinx = hstat->GetNbinsX();
  
  for(int ibinx = 1; ibinx <= nbinx; ibinx++){
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

void CorrectedYield_withSystematics(TH1D* &hstat, TH1D* &hsyst, TH1D* &hsystCorr, int ipart, TString* SaveAs_Title, TString* &texXtitle, TString* &texYtitle, double* pTbins, int nbinpT) {
  cout << "DEBUG CorrectedYield_withSystematics 1" << endl;

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

  TH3D* H3D_detectedV0s_TruePt = (TH3D*)file_O2Analysis->Get("lambdakzero-analysis-mc/h3dMass"+NamePart[ipart]+"_MC_truePt");
  TH3D* H3D_detectedV0s_DataLike = (TH3D*)file_O2Analysis->Get("lambdakzero-analysis-mc/h3dMass"+NamePart[ipart]);
  TH1D* TrueV0PtSpectrum = (TH1D*)file_O2Analysis->Get("lambdakzero-particle-count-mc/h"+NamePart[ipart]+"Count_PtDiff");
  TH1I* H1I_SelectedEventCount = (TH1I*)file_O2Analysis->Get("lambdakzero-analysis-mc/hSelectedEventCounter");

  TH1D* hSystematicUncertainty_Total = new TH1D("hSystematicUncertainty", "hSystematicUncertainty", nbinpT, pTbins);
  TH1D* hSystematicUncertainty_Total_PreBarlow = new TH1D("hSystematicUncertainty_PreBarlow", "hSystematicUncertainty_PreBarlow", nbinpT, pTbins);
  hSystematicUncertainty_Total->Sumw2();
  hSystematicUncertainty_Total_PreBarlow->Sumw2();

  TH1D* hSystematicUncertainty_Fit_SidebandVariation = new TH1D("hSystematicUncertainty_Fit_SidebandVariation", "hSystematicUncertainty_Fit_SidebandVariation", nbinpT, pTbins);
  TH1D* hSystematicUncertainty_Fit_SidebandVariation_PreBarlow = new TH1D("hSystematicUncertainty_Fit_SidebandVariation_PreBarlow", "hSystematicUncertainty_Fit_SidebandVariation_PreBarlow", nbinpT, pTbins);
  TH1D* hSystematicUncertainty_Fit_SignalExtractionType = new TH1D("hSystematicUncertainty_Fit_SignalExtractionType", "hSystematicUncertainty_Fit_SignalExtractionType", nbinpT, pTbins);
  TH1D* hSystematicUncertainty_Fit_SignalExtractionType_PreBarlow = new TH1D("hSystematicUncertainty_Fit_SignalExtractionType_PreBarlow", "hSystematicUncertainty_Fit_SignalExtractionType_PreBarlow", nbinpT, pTbins);
  TH1D* hSystematicUncertainty_Fit_BinningVariation = new TH1D("hSystematicUncertainty_Fit_BinningVariation", "hSystematicUncertainty_Fit_BinningVariation", nbinpT, pTbins);
  TH1D* hSystematicUncertainty_Fit_BinningVariation_PreBarlow = new TH1D("hSystematicUncertainty_Fit_BinningVariation_PreBarlow", "hSystematicUncertainty_Fit_BinningVariation_PreBarlow", nbinpT, pTbins);

  hSystematicUncertainty_Fit_SidebandVariation->Sumw2();
  hSystematicUncertainty_Fit_SidebandVariation_PreBarlow->Sumw2();
  hSystematicUncertainty_Fit_SignalExtractionType->Sumw2();
  hSystematicUncertainty_Fit_SignalExtractionType_PreBarlow->Sumw2();
  hSystematicUncertainty_Fit_BinningVariation->Sumw2();
  hSystematicUncertainty_Fit_BinningVariation_PreBarlow->Sumw2();

  TH1D *hMcEfficiency_vsPt;
  TH1D *hRawYield_vsPt;
  int InvMassRebinFactor = InvMassRebinFactor_standard[ipart];

  cout << "DEBUG CorrectedYield_withSystematics 2" << endl;

  //Systematics on corrected yield

  cout << "Calculating Systematics - Fit sideband variation" << endl;
  Get_Systematics_Fit_SideBandVariation(hSystematicUncertainty_Fit_SidebandVariation, hSystematicUncertainty_Fit_SidebandVariation_PreBarlow, ipart, H3D_detectedV0s_DataLike, H1I_SelectedEventCount, H3D_detectedV0s_TruePt, TrueV0PtSpectrum, file_O2Analysis, pTbins, nbinpT, InvMassRebinFactor);
  cout << "Calculating Systematics - Fit signal extraction type" << endl;
  Get_Systematics_Fit_SignalExtractionType(hSystematicUncertainty_Fit_SignalExtractionType, hSystematicUncertainty_Fit_SignalExtractionType_PreBarlow, ipart, H3D_detectedV0s_DataLike, H1I_SelectedEventCount, H3D_detectedV0s_TruePt, TrueV0PtSpectrum, file_O2Analysis, pTbins, nbinpT, InvMassRebinFactor); //very large systematics as is
  cout << "Calculating Systematics - Fit binning" << endl;
  Get_Systematics_Fit_Binning(hSystematicUncertainty_Fit_BinningVariation, hSystematicUncertainty_Fit_BinningVariation_PreBarlow, ipart, H3D_detectedV0s_DataLike, H1I_SelectedEventCount, H3D_detectedV0s_TruePt, TrueV0PtSpectrum, file_O2Analysis, pTbins, nbinpT);
  
  //systematics on cut variations over the N_variableCut variables that are being tested
  TH1D* hSystematicUncertainty_CutVariation[N_variableCut];
  TH1D* hSystematicUncertainty_CutVariation_PreBarlow[N_variableCut];
  for (int i_variableCut = 0; i_variableCut < N_variableCut; i_variableCut++) {
    hSystematicUncertainty_CutVariation[i_variableCut] = new TH1D("hSystematicUncertainty_CutVariation_"+VariableCut[i_variableCut], "hSystematicUncertainty_CutVariation_"+VariableCut[i_variableCut], nbinpT, pTbins);
    hSystematicUncertainty_CutVariation_PreBarlow[i_variableCut] = new TH1D("hSystematicUncertainty_CutVariation_PreBarlow_"+VariableCut[i_variableCut], "hSystematicUncertainty_CutVariation_PreBarlow_"+VariableCut[i_variableCut], nbinpT, pTbins);
    hSystematicUncertainty_CutVariation[i_variableCut]->Sumw2();
    hSystematicUncertainty_CutVariation_PreBarlow[i_variableCut]->Sumw2();
    cout << "*********************************** Calculating Systematics - Cut variation number " << i_variableCut << endl;
    Get_Systematics_OneCut(hSystematicUncertainty_CutVariation[i_variableCut], hSystematicUncertainty_CutVariation_PreBarlow[i_variableCut], ipart, file_O2Analysis_CutVariation_Collection_array[i_variableCut], nCutsIterations_array[i_variableCut], pTbins, nbinpT, InvMassRebinFactor);
  }

  cout << "DEBUG CorrectedYield_withSystematics 3" << endl;


  // Get_Systematics_OneCut_McTrueAndMcDatalikeDifferentCuts(hSystematicUncertainty, hSystematicUncertainty_PreBarlow, ipart, file_O2Analysis_CutVariation_Datalike_array, file_O2Analysis_CutVariation_True_array, id_referenceAnalysis, nCutsIterations, pTbins, nbinpT, InvMassRebinFactor);

  //Reference Corrected Yield
  double SignalMean, SignalStandardDev;
  cout << "Calculating Raw Yield for Ref" << endl;
  if (doUnfolding == false){
    Get_RawYield_vsPt_FeeddownCorrected(hRawYield_vsPt, H3D_detectedV0s_DataLike, H1I_SelectedEventCount, file_O2Analysis, &SignalMean, &SignalStandardDev, ipart, pTbins, nbinpT, 0, -1, 0, SideBandSizeMultiplierModifier_array[0], SignalExtractionType_default, InvMassRebinFactor); ///real line that need to take the place of the one above once it's been written
  }
  if (doUnfolding == true){
    Get_RawYield_vsPt_FeeddownCorrected_withUnfolding(hRawYield_vsPt, H3D_detectedV0s_DataLike, H1I_SelectedEventCount, file_O2Analysis, &SignalMean, &SignalStandardDev, ipart, pTbins, nbinpT, 0, -1, 0, SideBandSizeMultiplierModifier_array[0], SignalExtractionType_default, InvMassRebinFactor); ///real line that need to take the place of the one above once it's been written
  }
  cout << "Calculating Efficiency" << endl;
  Get_McEfficiency_vsPt(hMcEfficiency_vsPt, H3D_detectedV0s_TruePt, TrueV0PtSpectrum, ipart, pTbins, nbinpT, 0, -1, SideBandSizeMultiplierModifier_array[0], SignalMean, SignalStandardDev, InvMassRebinFactor);

  cout << "DEBUG CorrectedYield_withSystematics 4" << endl;

  // Add uncertainties in quadrature: multiply hists by themselves then add
  hSystematicUncertainty_Fit_SidebandVariation->Multiply(hSystematicUncertainty_Fit_SidebandVariation);
  hSystematicUncertainty_Fit_SignalExtractionType->Multiply(hSystematicUncertainty_Fit_SignalExtractionType);
  hSystematicUncertainty_Fit_BinningVariation->Multiply(hSystematicUncertainty_Fit_BinningVariation);
  hSystematicUncertainty_Total->Add(hSystematicUncertainty_Fit_SidebandVariation);
  hSystematicUncertainty_Total->Add(hSystematicUncertainty_Fit_SignalExtractionType);
  hSystematicUncertainty_Total->Add(hSystematicUncertainty_Fit_BinningVariation);
  //cut variation is not added for now
  // for (int i_variableCut = 0; i_variableCut < N_variableCut; i_variableCut++) {
  //   hSystematicUncertainty_CutVariation[i_variableCut]->Multiply(hSystematicUncertainty_CutVariation[i_variableCut]);
  //   hSystematicUncertainty_Total->Add(hSystematicUncertainty_CutVariation[i_variableCut]);
  // }
  cout << "DEBUG CorrectedYield_withSystematics 5" << endl;

  for(int ibinPt = 1; ibinPt <= nbinpT; ibinPt++){
    if (abs(hSystematicUncertainty_Total->GetBinContent(ibinPt)) > 1E-9) { //can't divide by 0
      // wrong errors; can probably get away with not looking at the errors of those systematic errors
      // hSystematicUncertainty_Total->SetBinError(ibinPt, 1./2 * 1./sqrt(hSystematicUncertainty_Total->GetBinContent(ibinPt)) * hSystematicUncertainty_Total->GetBinError(ibinPt)); // SigmaTot^2 = Sum(Sigma^2) -> error(SigmaTot) = 1/(2*SigmaTot)*Error(SigmaTot^2)  
      hSystematicUncertainty_Total->SetBinError(ibinPt, 0.); 
      hSystematicUncertainty_Total->SetBinContent(ibinPt,sqrt(hSystematicUncertainty_Total->GetBinContent(ibinPt)));   
    }
  }
  cout << "DEBUG CorrectedYield_withSystematics 6" << endl;


  // Get_Systematics_Full(); // gather all systematics; eventually

  // TH1D *hCorrectedYield_vsPt;
  // hCorrectedYield_vsPt = (TH1D*)hRawYield_vsPt->Clone("hCorrectedYield_vsPt");
  // hCorrectedYield_vsPt->Reset("M");
  TH1D* hCorrectedYield_vsPt = new TH1D("hCorrectedYield_vsPt", "hCorrectedYield_vsPt",nbinpT,pTbins);
  hCorrectedYield_vsPt->Sumw2();
  
  // // Lambda and AntiLambda Feeddown
  // TH1D *H1D_Feeddown_Correction = new TH1D("H1D_Feeddown_Correction", "H1D_Feeddown_Correction", nbinpT, pTbins);
  // if (ipart == 1 || ipart == 2) {
  //   Get_FeedDown(H1D_Feeddown_Correction, ipart, file_O2Analysis_CutVariation_array[id_referenceAnalysis], pTbins, nbinpT);
  //   H1D_Feeddown_Correction->Scale(-1.);
  //   hRawYield_vsPt->Add(H1D_Feeddown_Correction);
  //   // cout << "test00000" << endl;
  // }

  hstat = (TH1D*)hCorrectedYield_vsPt->Clone("hstat");
  hsyst = (TH1D*)hCorrectedYield_vsPt->Clone("hsyst");

  TH1* hCorrectedYield_vsPt_syst;

  if (doUnfolding == false){
    hCorrectedYield_vsPt->Divide(hRawYield_vsPt,hMcEfficiency_vsPt, 1., 1., ""); // not sure why I had b as an option (binomial) for the corrected yield
    hCorrectedYield_vsPt_syst = (TH1D*)hCorrectedYield_vsPt->Clone("hCorrectedYield_vsPt_syst");
  }
  if (doUnfolding == true){
    for(int ibinPt = 1; ibinPt <= nbinpT; ibinPt++){
      hCorrectedYield_vsPt->SetBinContent(ibinPt, hRawYield_vsPt->GetBinContent(ibinPt));
      hCorrectedYield_vsPt->SetBinError(ibinPt, hRawYield_vsPt->GetBinError(ibinPt));
    }
    hCorrectedYield_vsPt_syst = (TH1D*)hRawYield_vsPt->Clone("hCorrectedYield_vsPt_syst");
  }
  cout << "DEBUG CorrectedYield_withSystematics 7" << endl;

  cout << "final corrected distribution" << endl;
  for(int ibinPt = 1; ibinPt <= nbinpT; ibinPt++){
    hsyst->SetBinError(ibinPt,hSystematicUncertainty_Total->GetBinContent(ibinPt));
    // hsyst->SetBinError(ibinPt,hSystematicUncertainty_PreBarlow->GetBinContent(ibinPt));
    // hsyst->SetBinError(ibinPt,0.00001);
    cout << "ibinPt = "<< ibinPt << ", rawYield = " << hRawYield_vsPt->GetBinContent(ibinPt) << ", Efficiency = " << hMcEfficiency_vsPt->GetBinContent(ibinPt) << endl;
    hCorrectedYield_vsPt_syst->SetBinError(ibinPt, hSystematicUncertainty_Total->GetBinContent(ibinPt));
    hCorrectedYield_vsPt->SetBinError(ibinPt, sqrt(pow(hCorrectedYield_vsPt->GetBinError(ibinPt),2) + pow(hSystematicUncertainty_Total->GetBinContent(ibinPt),2)));
  }
  cout << "DEBUG CorrectedYield_withSystematics 8" << endl;

  /////////// Draw Efficiency /////////////
  TString* pdfName0 = new TString("Efficiency_"+NamePart[ipart]);;
  // *pdfName0 += "Efficiency_"+NamePart[ipart];

  TH1* histograms_collection_0[1] = {hMcEfficiency_vsPt};
  TString* pdfLegend_collection_0[1] = {new TString("Efficiency "+NamePart_Latex[ipart])};
  int ratioOption0 = 0;
  texXtitle = texPtX;
  texYtitle = texEfficiency;
  Draw_TH1_Histograms_Efficiency(histograms_collection_0, pdfLegend_collection_0, 1, pdfName0, texXtitle, texYtitle, pTbins, nbinpT, ratioOption0);


  /////////// Draw Corrected yield comparison O2 vs HEPdata /////////////
  TString* pdfName1 = new TString("CorrectewdYields_O2_vs_MCgenYield_"+NamePart[ipart]);;
  // *pdfName1 += ;
  // O2 MC true yield
  // TH1I* H1I_TotalEventCount = (TH1I*)file_O2Analysis->Get("lambdakzero-builder/hEventCounter");
  TH1I* H1I_TotalEventCount = (TH1I*)file_O2Analysis->Get("lambdakzero-particle-count-mc/hSelAndRecoMcCollCounter");

  TrueV0PtSpectrum->Sumw2();
  TH1D* TrueV0PtSpectrum_rebinned = (TH1D*)TrueV0PtSpectrum->Rebin(nbinpT,"TrueV0PtSpectrum_rebinned",pTbins);
  TH1D* TrueV0PtSpectrum_rebinned_syst = (TH1D*)TrueV0PtSpectrum->Rebin(nbinpT,"TrueV0PtSpectrum_rebinned_syst",pTbins);
  for (int ibinPt = 1; ibinPt <= nbinpT; ibinPt++)
  {
    double dpT = TrueV0PtSpectrum_rebinned->GetXaxis()->GetBinWidth(ibinPt);
    double drapidity = 1.5; //1.5
    double d2N_dpTdy = TrueV0PtSpectrum_rebinned->GetBinContent(ibinPt) *1./dpT*1./drapidity;
    double d2N_dpTdy_error = sqrt(TrueV0PtSpectrum_rebinned->GetBinContent(ibinPt)) * 1./dpT*1./drapidity;
    TrueV0PtSpectrum_rebinned->SetBinContent(ibinPt, d2N_dpTdy);
    TrueV0PtSpectrum_rebinned->SetBinError(ibinPt, d2N_dpTdy_error);
    TrueV0PtSpectrum_rebinned_syst->SetBinContent(ibinPt, d2N_dpTdy);
    TrueV0PtSpectrum_rebinned_syst->SetBinError(ibinPt, 0.);
  }
  TrueV0PtSpectrum_rebinned->Scale(1./H1I_TotalEventCount->GetEntries());
  TrueV0PtSpectrum_rebinned_syst->Scale(1./H1I_TotalEventCount->GetEntries());

  // TH1* histograms_collection_1[4] = {hCorrectedYield_vsPt_syst, hCorrectedYield_vsPt, hPtDiff_RawSpectrum_HEP_rebinned, TrueV0PtSpectrum_rebinned};
  // TString* pdfLegend_collection_1[4] = {new TString("hCorrectedYield_vsPt_syst"), new TString("Corrected Yield"), new TString("published pp #sqrt{#it{s}} = 900 GeV"),new TString("True MC")};
  // int ratioOption1 = 0;
  // texXtitle = texPtX;
  // texYtitle = texptDifferentialYield;
  // Draw_TH1_Histograms_withOneHistSystematics(histograms_collection_1, pdfLegend_collection_1, 4, pdfName1, texXtitle, texYtitle, pTbins, nbinpT, ratioOption1);
  TH1* histograms_collection_1[3] = {hCorrectedYield_vsPt_syst, hCorrectedYield_vsPt, TrueV0PtSpectrum_rebinned};
  TString* pdfLegend_collection_1[3] = {new TString("hCorrectedYield_vsPt_syst"), new TString("Corrected Spectrum"),new TString("True MC")};
  int ratioOption1 = 0;
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
  int ratioOption2 = 1; //to draw line at ratio=1
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
  int nbinx = hstat->GetNbinsX();
  
  for(int ibinx = 1; ibinx <= nbinx; ibinx++){
    hsystCorr->SetBinContent (ibinx, hstat->GetBinContent(ibinx));
    hsystCorr->SetBinError   (ibinx, hstat->GetBinContent(ibinx)*0.0);
  }
  *SaveAs_Title += "CorrectedYield_withSystematics";
  texXtitle = texPtX;
  texYtitle = texptDifferentialYield;
}


void Draw_TH1_Histograms(TH1** histograms_collection, TString** pdfLegend_collection, int collectionSize, TString* pdfName, TString* &texXtitle, TString* &texYtitle, double* bins, int nbins, int ratioOption) {
  
  float MaximumY = 0;
  for (int i = 0; i < collectionSize; i++) {
    float maxTemporary = histograms_collection[i]->GetMaximum();
    if (maxTemporary > MaximumY) {
      MaximumY = maxTemporary;
    }
  }

  ////plots
  TCanvas *canvas = new TCanvas ("canvas"+*pdfName, NamehistoInvMass[ipart], 1200, 800);
  canvas->cd(0); //ATTENTION: +1 because canvas divide starts counting from 1
  TH1 *hFrame = canvas->DrawFrame(bins[0],0,bins[nbins],1.4*MaximumY);
  // TH1 *hFrame = canvas->DrawFrame(bins[0],0,bins[nbins],0.5);
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
  for (int i = 0; i < collectionSize; i++) {
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


void Draw_TH1_Histograms_Efficiency(TH1** histograms_collection, TString** pdfLegend_collection, int collectionSize, TString* pdfName, TString* &texXtitle, TString* &texYtitle, double* bins, int nbins, int ratioOption) {
  
  float MaximumY = 0;
  for (int i = 0; i < collectionSize; i++) {
    float maxTemporary = histograms_collection[i]->GetMaximum();
    if (maxTemporary > MaximumY) {
      MaximumY = maxTemporary;
    }
  }

  ////plots
  TCanvas *canvas = new TCanvas ("canvas"+*pdfName, NamehistoInvMass[ipart], 1200, 800);
  canvas->cd(0); //ATTENTION: +1 because canvas divide starts counting from 1
  // TH1 *hFrame = canvas->DrawFrame(bins[0],0,bins[nbins],1.4*MaximumY);
  TH1 *hFrame;
  if (UseFinder) {
    hFrame = canvas->DrawFrame(0,0,bins[nbins],0.73);
  }
  else {
    hFrame = canvas->DrawFrame(0,0,bins[nbins],0.25);
  }

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
  for (int i = 0; i < collectionSize; i++) {
    histograms_collection[i]->Draw("same");
    histograms_collection[i]->SetMarkerStyle(markers[i]);
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


void Draw_TH1_Histograms_withOneHistSystematics(TH1** histograms_collection, TString** pdfLegend_collection, int collectionSize, TString* pdfName, TString* &texXtitle, TString* &texYtitle, double* bins, int nbins, int ratioOption) {
  //draws only one hist with systematics, has to be the second one in the collection, and the associated systematics is taken from first hist in collection
  float MaximumY = 0;
  for (int i = 1; i < collectionSize; i++) {
    float maxTemporary = histograms_collection[i]->GetMaximum();
    if (maxTemporary > MaximumY) {
      MaximumY = maxTemporary;
    }
  }

  ////plots
  TCanvas *canvas = new TCanvas ("canvas"+*pdfName, NamehistoInvMass[ipart], 1200, 800);
  canvas->cd(0); //ATTENTION: +1 because canvas divide starts counting from 1
  TH1 *hFrame = canvas->DrawFrame(0.,0,bins[nbins],1.4*MaximumY);
  hFrame->SetXTitle(texXtitle->Data());
  hFrame->SetYTitle(texYtitle->Data());
  TLegend * leg = new TLegend(0.59, 0.74, 0.87, 0.87);

  //draw line at ratio=1
  if (ratioOption == 1) {
    // TLine *l=new TLine(histograms_collection[0]->GetXaxis()->GetBinLowEdge(0), 1.0, histograms_collection[0]->GetXaxis()->GetBinUpEdge(nbins), 1.0);
    TLine *l=new TLine(0, 1.0, bins[nbins], 1.0);
    l->SetLineColor(kGray);
    l->Draw();
  }

  //draw histograms from collection, ignoring first one that is the systematics for second one
  for (int i = 1; i < collectionSize; i++) {
    histograms_collection[i]->Draw("same");
    histograms_collection[i]->SetMarkerStyle(markers[i-1]);
    histograms_collection[i]->SetMarkerColor(colors[i-1]);
    histograms_collection[i]->SetLineColor(colors[i-1]);

    leg->AddEntry(histograms_collection[i], *pdfLegend_collection[i],   "LP");
  }

  // if (ipart == 0) {
  //   //draw systematics for 2d hist using 1st hist in collection
  //   histograms_collection[0]->Draw("E2,same");
  //   histograms_collection[0]->SetMarkerStyle(markers[0]);
  //   histograms_collection[0]->SetMarkerColor(colors[0]);
  //   histograms_collection[0]->SetLineColor(colors[0]);
  //   histograms_collection[0]->SetFillStyle(0); // To draw empty boxes
  // }

  // leg->SetTextSize(gStyle->GetTextSize()*0.5);
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

void Draw_TH1_Histograms_withOneHistSystematics_ratio(TH1** histograms_collection, TString** pdfLegend_collection, int collectionSize, TString* pdfName, TString* &texXtitle, TString* &texYtitle, double* bins, int nbins, int ratioOption) {
  //draws only one hist with systematics, has to be the second one in the collection, and the associated systematics is taken from first hist in collection

  ////plots
  TCanvas *canvas = new TCanvas ("canvas"+*pdfName, NamehistoInvMass[ipart], 1200, 800);
  canvas->cd(0); //ATTENTION: +1 because canvas divide starts counting from 1
  TH1 *hFrame = canvas->DrawFrame(0,0.,bins[nbins],2);
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
  for (int i = 1; i < collectionSize; i++) {
    histograms_collection[i]->Draw("same");
    histograms_collection[i]->SetMarkerStyle(markers[i-1]);
    histograms_collection[i]->SetMarkerColor(colors[i-1]);
    histograms_collection[i]->SetLineColor(colors[i-1]);

    leg->AddEntry(histograms_collection[i], *pdfLegend_collection[i],   "LP");
  }

  // if (ipart == 0) {
  //   //draw systematics for 2d hist using 1st hist in collection
  //   histograms_collection[0]->Draw("E2,same");
  //   histograms_collection[0]->SetMarkerStyle(markers[0]);
  //   histograms_collection[0]->SetMarkerColor(colors[0]);
  //   histograms_collection[0]->SetLineColor(colors[0]);
  //   histograms_collection[0]->SetFillStyle(0); // To draw empty boxes
  // }

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

void Test_LambdaTrueMCYieldvsHEP(TH1D* &hstat, TH1D* &hsyst, TH1D* &hsystCorr, int ipart, TFile* file_O2Analysis, TString* SaveAs_Title, TString* &texXtitle, TString* &texYtitle, double* pTbins, int nbinpT) {

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
  for (int ibinPt = 0; ibinPt < nbinpT; ibinPt++)
  {
    double dpT = TrueV0PtSpectrum_rebinned->GetXaxis()->GetBinWidth(ibinPt);
    double drapidity = 1.5; //1.5
    double d2N_dpTdy = TrueV0PtSpectrum_rebinned->GetBinContent(ibinPt) *1./dpT*1./drapidity; //note from Aimeric du futur: should be sqrt(getbincontent)
    TrueV0PtSpectrum_rebinned->SetBinContent(ibinPt, d2N_dpTdy);
  }
  
  TString* pdfName = new TString("");;
  *pdfName += "CorrectewdYields_O2_vs_HEPdata";

  TH1* histograms_collection[2] = {TrueV0PtSpectrum_rebinned, hPtDiff_RawSpectrum_HEP[ipart]};
  TString* pdfLegend_collection[2] = {new TString("TrueV0PtSpectrum_rebinned"),new TString("hPtDiff_RawSpectrum_HEP")};
  int ratioOption = 0;
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
  int nbinx = hstat->GetNbinsX();
  
  for(int ibinx = 1; ibinx <= nbinx; ibinx++){
    hsyst->SetBinContent (ibinx, hstat->GetBinContent(ibinx));
    hsyst->SetBinError   (ibinx, hstat->GetBinContent(ibinx)*0.08);
    hsystCorr->SetBinContent (ibinx, hstat->GetBinContent(ibinx));
    hsystCorr->SetBinError   (ibinx, hstat->GetBinContent(ibinx)*0.15);
  }
  *SaveAs_Title += "ToBeDeleted";
  texXtitle = texPtX;
  texYtitle = texptDifferentialYield_HEPratio_pseudoEff;
}

void Test_V0DaugtherPairsTrackingEfficiency(TH1D* &hstat, TH1D* &hsyst, TH1D* &hsystCorr, int ipart, TFile* file_O2Analysis, TString* SaveAs_Title, TString* &texXtitle, TString* &texYtitle, double* pTbins, int nbinpT) {


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
  int nbinx = hstat->GetNbinsX();
  
  for(int ibinx = 1; ibinx <= nbinx; ibinx++){
    hsyst->SetBinContent (ibinx, hstat->GetBinContent(ibinx));
    hsyst->SetBinError   (ibinx, hstat->GetBinContent(ibinx)*0.000001);
    hsystCorr->SetBinContent (ibinx, hstat->GetBinContent(ibinx));
    hsystCorr->SetBinError   (ibinx, hstat->GetBinContent(ibinx)*0.15);
  }
  *SaveAs_Title += "DaughterPairsRecoEfficiency";
  texXtitle = texPtX;
  texYtitle = texDaughterPairsRecoEfficiency;
}


void Test_V0DaugtherPairsTrackingEfficiency_PostDcaFitter(TH1D* &hstat, TH1D* &hsyst, TH1D* &hsystCorr, int ipart, TFile* file_O2Analysis, TString* SaveAs_Title, TString* &texXtitle, TString* &texYtitle, double* pTbins, int nbinpT) {


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
  int nbinx = hstat->GetNbinsX();
  
  for(int ibinx = 1; ibinx <= nbinx; ibinx++){
    hsyst->SetBinContent (ibinx, hstat->GetBinContent(ibinx));
    hsyst->SetBinError   (ibinx, hstat->GetBinContent(ibinx)*0.000001);
    hsystCorr->SetBinContent (ibinx, hstat->GetBinContent(ibinx));
    hsystCorr->SetBinError   (ibinx, hstat->GetBinContent(ibinx)*0.15);
  }
  *SaveAs_Title += "DaughterPairsRecoEfficiency_PostDcaFitter";
  texXtitle = texPtX;
  texYtitle = texDaughterPairsRecoEfficiency;
}


void Draw_Fit_Mass_Histograms(TH1D* &hstat, TH1D* &hsyst, TH1D* &hsystCorr, TFile *file_O2Analysis, TString* SaveAs_Title, TString* &texXtitle, TString* &texYtitle, int ipart, double* pTbins, int nbinpT) {
  //Looks within the signal area of the fit, 6Sigma

  // what is fed as input: TH3D* H3D_detectedV0s_DataLike = (TH3D*)file_O2Analysis->Get("lambdakzero-analysis-mc/h3dMass"+NamePart[ipart]);
  TH3D* H3D_detectedV0s_DataLike = (TH3D*)file_O2Analysis->Get("lambdakzero-analysis-mc/h3dMass"+NamePart[ipart]);
  // TH3D* H3D_detectedV0s_DataLike = (TH3D*)file_O2Analysis->Get("lambdakzero-analysis-mc/h3dMass"+NamePart[ipart]+"_MC_truePt");

  TH3D* H3D_detectedV0s_DataLike_rebinnedZ = (TH3D*)H3D_detectedV0s_DataLike->RebinZ(InvMassRebinFactor_standard[ipart],"H3D_detectedV0s_DataLike_rebinnedZ");
  cout << "---------------------------------------- InvMassRebinFactor_standard = " << InvMassRebinFactor_standard[ipart] << endl;
  ////////////////////////////////// Fit initialisation //////////////////////////////////
  //Fit tools initialisation
  TF1 *gauss[nbinpT];
  TF1 *GaussPlusPolynom[nbinpT];
  TF1 *bkgparab[nbinpT];
  TF1 *bkgparab_left[nbinpT];
  TF1 *bkgparab_right[nbinpT];
  TFitResultPtr fFitResult[nbinpT];
  // double parGaussParab[nbinpT][6]; //parab backround
  double parGaussParab[nbinpT][5]; //linear background

  TH1D *H1D_detectedV0s_DataLike_SmallpTinterval[nbinpT];
  TH1D *hPtDiff_RawYield_temp = new TH1D("hPtDiff_RawYield_temp","hPtDiff_RawYield_temp",nbinpT,pTbins);

  TCanvas *canvasMassvsPtFits = new TCanvas ("canvasMassPt_fits", NamehistoInvMass[ipart], 800, 600);
  // TCanvas *canvasMassvsPt = new TCanvas ("canvasMassPt", NamehistoInvMass[ipart], 1600, 1200);
  if (ipart == 0) {
  canvasMassvsPtFits->Divide(std::ceil(std::sqrt(nbinpT)),std::ceil(std::sqrt(nbinpT))-1);
  }
  else {
  canvasMassvsPtFits->Divide(std::ceil(std::sqrt(nbinpT)),std::ceil(std::sqrt(nbinpT)));
  }
  TH1 *hFrameFits[nbinpT];
  cout << "testAAA" << endl;
  for(int ibinPt = 1; ibinPt <= nbinpT; ibinPt++){
    int iPt = ibinPt-1; //different convention for bins that start at 1 and arrays that start at 0
    
    canvasMassvsPtFits->cd(ibinPt); //ATTENTION: +1 because canvas divide starts counting from 1

    int ibinPt_originalAxis_low = H3D_detectedV0s_DataLike->GetYaxis()->FindBin(pTbins[iPt]);
    int ibinPt_originalAxis_up = H3D_detectedV0s_DataLike->GetYaxis()->FindBin(pTbins[iPt+1]) - 1;
    // cout << "ibinPt: " << ibinPt << "; originalAxis [low;up] = [" << ibinPt_originalAxis_low << ";" << ibinPt_originalAxis_up << "]" << endl;

    H1D_detectedV0s_DataLike_SmallpTinterval[iPt] = (TH1D*)H3D_detectedV0s_DataLike_rebinnedZ->ProjectionZ("Get_RawYield_vsPt_InvMass_"+NamePart[ipart]+Form("_Pt%i", ibinPt),0,-1,ibinPt_originalAxis_low,ibinPt_originalAxis_up);
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

    H1D_detectedV0s_DataLike_SmallpTinterval[iPt]->Fit(gauss[iPt], "R0QL");

    float lowEdge_fit = gauss[iPt]->GetParameter(1) - 2* DefaultSideBandSize_array[ipart] * gauss[iPt]->GetParameter(2);
    float upEdge_fit = gauss[iPt]->GetParameter(1) + 2* DefaultSideBandSize_array[ipart] * gauss[iPt]->GetParameter(2);
    // cout << "lowEdge_fit = " << lowEdge_fit << ", upEdge_fit = " << upEdge_fit << endl;

    // Linear background; less issues with brackground fit going negative; better for pp
    // bkgparab[iPt] = new TF1("line_"+NamePart[ipart]+Form("_Pt%i", ibinPt), fline, liminf[ipart], limsup[ipart], 3);
    bkgparab[iPt] = new TF1("line_"+NamePart[ipart]+Form("_Pt%i", ibinPt), fline, lowEdge_fit, upEdge_fit, 3);
    bkgparab[iPt]->SetLineColor(kGreen);
    bkgparab[iPt]->FixParameter(2, ipart);

    // GaussPlusPolynom[iPt] = new TF1("GaussPlusPolynom_"+NamePart[ipart]+Form("_Pt%i", ibinPt), "gaus(0)+pol1(3)", liminf[ipart], limsup[ipart]);
    GaussPlusPolynom[iPt] = new TF1("GaussPlusPolynom_"+NamePart[ipart]+Form("_Pt%i", ibinPt), "gaus(0)+pol1(3)", lowEdge_fit, upEdge_fit);
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
    GaussPlusPolynom[iPt]->SetParLimits(2, 0.001, 0.05);
    // cout << "par0 upper limit" <<  1.1*H1D_detectedV0s_DataLike_SmallpTinterval[iPt]->GetBinContent(H1D_detectedV0s_DataLike_SmallpTinterval[iPt]->GetMaximumBin()) << endl;

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
    int icolor=0;
    H1D_detectedV0s_DataLike_SmallpTinterval[iPt]->SetMarkerStyle(markers[0]);
    H1D_detectedV0s_DataLike_SmallpTinterval[iPt]->SetMarkerColor(colors [icolor]);
    H1D_detectedV0s_DataLike_SmallpTinterval[iPt]->SetLineColor  (colors [icolor]);
    H1D_detectedV0s_DataLike_SmallpTinterval[iPt]->SetMarkerSize (0.5);

    bkgparab[iPt]->SetLineColor(kGreen);
    bkgparab[iPt]->SetFillColor(kGreen);
    bkgparab[iPt]->Draw("BSAME");
    H1D_detectedV0s_DataLike_SmallpTinterval[iPt]->Draw("E,same");


    // for two background regions:
    float lowerleftupEdge_fit = gauss[iPt]->GetParameter(1) - 2* DefaultSideBandSize_array[ipart] * gauss[iPt]->GetParameter(2);
    float upperleftupEdge_fit = gauss[iPt]->GetParameter(1) - 1* DefaultSideBandSize_array[ipart] * gauss[iPt]->GetParameter(2);
    float lowerrightlowEdge_fit = gauss[iPt]->GetParameter(1) + 1* DefaultSideBandSize_array[ipart] * gauss[iPt]->GetParameter(2);
    float upperrightlowEdge_fit = gauss[iPt]->GetParameter(1) + 2* DefaultSideBandSize_array[ipart] * gauss[iPt]->GetParameter(2);
        
    TLine *lowerleft_separation=new TLine(lowerleftupEdge_fit, 0, lowerleftupEdge_fit, 0.3*H1D_detectedV0s_DataLike_SmallpTinterval[iPt]->GetMaximum());//(x1,y1) (x2,y2) are the two points connected by the line
    TLine *upperleft_separation=new TLine(upperleftupEdge_fit, 0, upperleftupEdge_fit, H1D_detectedV0s_DataLike_SmallpTinterval[iPt]->GetMaximum());//(x1,y1) (x2,y2) are the two points connected by the line
    TLine *lowerright_separation=new TLine(lowerrightlowEdge_fit, 0, lowerrightlowEdge_fit, H1D_detectedV0s_DataLike_SmallpTinterval[iPt]->GetMaximum());//(x1,y1) (x2,y2) are the two points connected by the line
    TLine *upperright_separation=new TLine(upperrightlowEdge_fit, 0, upperrightlowEdge_fit, 0.3*H1D_detectedV0s_DataLike_SmallpTinterval[iPt]->GetMaximum());//(x1,y1) (x2,y2) are the two points connected by the line
    lowerleft_separation->SetLineColor(kGray);
    upperleft_separation->SetLineColor(kGray);
    lowerright_separation->SetLineColor(kGray);
    upperright_separation->SetLineColor(kGray);
    lowerleft_separation->Draw("same");
    upperleft_separation->Draw("same");
    lowerright_separation->Draw("same");
    upperright_separation->Draw("same");


    // bkgparab_left[iPt] = new TF1("line_"+NamePart[ipart]+Form("_Pt%i", ibinPt), fline, lowEdge_fit, leftupEdge_fit, 3);
    // bkgparab_right[iPt] = new TF1("line_"+NamePart[ipart]+Form("_Pt%i", ibinPt), fline, rightlowEdge_fit, upEdge_fit, 3);
    // bkgparab_left[iPt]->SetParameters(&parGaussParab[iPt][3]);
    // bkgparab_right[iPt]->SetParameters(&parGaussParab[iPt][3]);

    // bkgparab_left[iPt]->SetLineColor(kGreen);
    // bkgparab_left[iPt]->SetFillColor(kGreen);
    // bkgparab_left[iPt]->Draw("BSAME");
    // bkgparab_right[iPt]->SetLineColor(kGreen);
    // bkgparab_right[iPt]->SetFillColor(kGreen);
    // bkgparab_right[iPt]->Draw("BSAME");


    // TLatex * text_pTbin = new TLatex (0.62,0.75,Form("pT [%.1f;%.1f] GeV", pTbins[iPt], pTbins[iPt+1]));
    // text_pTbin->SetNDC(kTRUE);
    // text_pTbin->Draw();
    // TLatex * textContext = new TLatex (0.18,0.82,"#bf{ALICE Performance}"); //BOLD
    // textContext->SetTextSize(0.05);
    // textContext->SetNDC(kTRUE); //remove if I want x,y in TLatex to be in the coordinate system of the histogram
    // textContext->Draw();
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


    //test count

    int ibin_liminf = H1D_detectedV0s_DataLike_SmallpTinterval[iPt]->GetXaxis()->FindBin(upperleftupEdge_fit);
    int ibin_limsup = H1D_detectedV0s_DataLike_SmallpTinterval[iPt]->GetXaxis()->FindBin(upperrightlowEdge_fit);
    long int SignalCountInSignalRegion = H1D_detectedV0s_DataLike_SmallpTinterval[iPt]->Integral(ibin_liminf,ibin_limsup);
    long int SignalCountFull = H1D_detectedV0s_DataLike_SmallpTinterval[iPt]->Integral(0,-1);
    cout << "iPt = " << iPt<< endl;
    cout << "SignalCountInSignalRegion = " << SignalCountInSignalRegion << ", SignalCountFull = " << SignalCountFull<< endl;
    cout << "RatioInside = " << SignalCountInSignalRegion*1./SignalCountFull << endl;
  }

  canvasMassvsPtFits->SaveAs("SignalExtraction_"+NamePart[ipart]+".pdf");

  delete hPtDiff_RawYield_temp;

  /////////////////////////////////Does Nothing - only there to respect template/////////////////////////////////////////////
  hstat = (TH1D*)H1D_detectedV0s_DataLike_SmallpTinterval[0]->Clone("hstat");
  hstat->Sumw2();
  hstat->Reset("M");

  //////Error Bars///////
  hsyst     = (TH1D*)hstat->Clone("hsyst");
  hsyst->Reset("M");
  hsystCorr = (TH1D*)hstat->Clone("hsystCorr");
  hsystCorr->Reset("M");
  hstat     ->Sumw2(); 
  hsyst     ->Sumw2();
  hsystCorr ->Sumw2(); 
  int nbinx = hstat->GetNbinsX();
  
  for(int ibinx = 1; ibinx <= nbinx; ibinx++){
    hsyst->SetBinContent (ibinx, hstat->GetBinContent(ibinx));
    hsyst->SetBinError   (ibinx, hstat->GetBinContent(ibinx)*0.000001);
    hsystCorr->SetBinContent (ibinx, hstat->GetBinContent(ibinx));
    hsystCorr->SetBinError   (ibinx, hstat->GetBinContent(ibinx)*0.15);
  }
  *SaveAs_Title += "delete";
  texXtitle = texPtX;
  texYtitle = texDaughterPairsRecoEfficiency;
}



void Draw_Mass_Histograms(TH1D* &hstat, TH1D* &hsyst, TH1D* &hsystCorr, TFile *file_O2Analysis, TString* SaveAs_Title, TString* &texXtitle, TString* &texYtitle, int ipart, double* pTbins, int nbinpT) {
  //Looks within the signal area of the fit, 6Sigma

  // what is fed as input: TH3D* H3D_detectedV0s_DataLike = (TH3D*)file_O2Analysis->Get("lambdakzero-analysis-mc/h3dMass"+NamePart[ipart]);
  TH3D* H3D_detectedV0s_DataLike = (TH3D*)file_O2Analysis->Get("lambdakzero-analysis-mc/h3dMass"+NamePart[ipart]);

  TH3D* H3D_detectedV0s_DataLike_rebinnedZ = (TH3D*)H3D_detectedV0s_DataLike->RebinZ(InvMassRebinFactor_standard[ipart],"H3D_detectedV0s_DataLike_rebinnedZ");
  cout << "---------------------------------------- InvMassRebinFactor_standard = " << InvMassRebinFactor_standard[ipart] << endl;


  TH1D *H1D_detectedV0s_DataLike_SmallpTinterval[nbinpT];

  TCanvas *canvasMassvsPtFits = new TCanvas ("canvasMassPt_fits", NamehistoInvMass[ipart], 800, 600);
  // TCanvas *canvasMassvsPt = new TCanvas ("canvasMassPt", NamehistoInvMass[ipart], 1600, 1200);
  if (ipart == 0) {
  canvasMassvsPtFits->Divide(std::ceil(std::sqrt(nbinpT)),std::ceil(std::sqrt(nbinpT))-1);
  }
  else {
  canvasMassvsPtFits->Divide(std::ceil(std::sqrt(nbinpT)),std::ceil(std::sqrt(nbinpT)));
  }
  TH1 *hFrameFits[nbinpT];
  for(int ibinPt = 1; ibinPt <= nbinpT; ibinPt++){
    int iPt = ibinPt-1; //different convention for bins that start at 1 and arrays that start at 0
    
    canvasMassvsPtFits->cd(ibinPt); //ATTENTION: +1 because canvas divide starts counting from 1

    int ibinPt_originalAxis_low = H3D_detectedV0s_DataLike->GetYaxis()->FindBin(pTbins[iPt]);
    int ibinPt_originalAxis_up = H3D_detectedV0s_DataLike->GetYaxis()->FindBin(pTbins[iPt+1]) - 1;
    // cout << "ibinPt: " << ibinPt << "; originalAxis [low;up] = [" << ibinPt_originalAxis_low << ";" << ibinPt_originalAxis_up << "]" << endl;

    H1D_detectedV0s_DataLike_SmallpTinterval[iPt] = (TH1D*)H3D_detectedV0s_DataLike_rebinnedZ->ProjectionZ("Get_RawYield_vsPt_InvMass_"+NamePart[ipart]+Form("_Pt%i", ibinPt),0,-1,ibinPt_originalAxis_low,ibinPt_originalAxis_up);
    hFrameFits[iPt] = canvasMassvsPtFits->DrawFrame(LowerLimitDisplay_Mass[ipart],0,UpperLimitDisplay_Mass[ipart],1.45*H1D_detectedV0s_DataLike_SmallpTinterval[iPt]->GetMaximum());
    hFrameFits[iPt]->SetXTitle(texInvMassDecays_titles[ipart]->Data());
    hFrameFits[iPt]->SetYTitle("Counts");
    H1D_detectedV0s_DataLike_SmallpTinterval[iPt]->Draw("E,same");

    H1D_detectedV0s_DataLike_SmallpTinterval[iPt]->GetXaxis()->SetTitle(texInvMassDecays_titles[ipart]->Data());
    H1D_detectedV0s_DataLike_SmallpTinterval[iPt]->GetYaxis()->SetTitle("Counts");
    int icolor=0;
    H1D_detectedV0s_DataLike_SmallpTinterval[iPt]->SetMarkerStyle(markers[0]);
    H1D_detectedV0s_DataLike_SmallpTinterval[iPt]->SetMarkerColor(colors [icolor]);
    H1D_detectedV0s_DataLike_SmallpTinterval[iPt]->SetLineColor  (colors [icolor]);
    H1D_detectedV0s_DataLike_SmallpTinterval[iPt]->SetMarkerSize (0.5);
    // TLatex * text_pTbin = new TLatex (0.62,0.75,Form("pT [%.1f;%.1f] GeV", pTbins[iPt], pTbins[iPt+1]));
    // text_pTbin->SetNDC(kTRUE);
    // text_pTbin->Draw();
    // TLatex * textContext = new TLatex (0.18,0.82,"#bf{ALICE Performance}"); //BOLD
    // textContext->SetTextSize(0.05);
    // textContext->SetNDC(kTRUE); //remove if I want x,y in TLatex to be in the coordinate system of the histogram
    // textContext->Draw();
    // TLatex * textColl = new TLatex (0.18,0.75,"pp #sqrt{#it{s}} = 0.9 TeV, pilot beam 2021");
    TLatex * textColl = new TLatex (0.18,0.75,"pp #sqrt{#it{s}} = 0.9 TeV, Monte Carlo");
    textColl->SetTextSize(0.04);
    textColl->SetNDC(kTRUE); //remove if I want x,y in TLatex to be in the coordinate system of the histogram
    textColl->Draw();
    // TLatex * text_part = new TLatex (0.18,0.70,NamePart_Latex[ipart]+",   0 < #it{p}_{T} < 10 GeV/#it{c}");
    TLatex * text_part = new TLatex (0.18,0.70,NamePart_Latex[ipart]);    
    text_part->SetTextSize(0.04);
    text_part->SetNDC(kTRUE); //remove if I want x,y in TLatex to be in the coordinate system of the histogram
    text_part->Draw();
    // TLatex * text_fit = new TLatex (0.18,0.63,Form("Gaussian Fit:"));
    // text_fit->SetTextSize(0.035);
    // text_fit->SetNDC(kTRUE);
    // text_fit->Draw();
    // // TLatex * text_fitParam_mu = new TLatex (0.18,0.59,Form("#it{#mu} = %.0f #pm %.2f MeV/#it{c}^{2}", parGaussParab[iPt][1]*1000., parGaussParabError_mu[iPt]*1000.));
    // TLatex * text_fitParam_mu = new TLatex (0.18,0.59,Form("#it{#mu} = %.0f MeV/#it{c}^{2}", parGaussParab[iPt][1]*1000.));
    // text_fitParam_mu->SetTextSize(0.035);
    // text_fitParam_mu->SetNDC(kTRUE);
    // text_fitParam_mu->Draw();
    // // TLatex * text_fitParam_sigma = new TLatex (0.18,0.55,Form("#it{#sigma} = %.1f #pm %.2f MeV/#it{c}^{2}", parGaussParab[iPt][2]*1000., parGaussParabError_sigma[iPt]*1000.));
    // TLatex * text_fitParam_sigma = new TLatex (0.18,0.55,Form("#it{#sigma} = %.1f MeV/#it{c}^{2}", parGaussParab[iPt][2]*1000.));
    // text_fitParam_sigma->SetTextSize(0.035);
    // text_fitParam_sigma->SetNDC(kTRUE);
    // text_fitParam_sigma->Draw();

  }

  canvasMassvsPtFits->SaveAs("InvMass_"+NamePart[ipart]+".pdf");

  /////////////////////////////////Does Nothing - only there to respect template/////////////////////////////////////////////
  hstat = (TH1D*)H1D_detectedV0s_DataLike_SmallpTinterval[0]->Clone("hstat");
  hstat->Sumw2();
  hstat->Reset("M");

  //////Error Bars///////
  hsyst     = (TH1D*)hstat->Clone("hsyst");
  hsyst->Reset("M");
  hsystCorr = (TH1D*)hstat->Clone("hsystCorr");
  hsystCorr->Reset("M");
  hstat     ->Sumw2(); 
  hsyst     ->Sumw2();
  hsystCorr ->Sumw2(); 
  int nbinx = hstat->GetNbinsX();
  
  for(int ibinx = 1; ibinx <= nbinx; ibinx++){
    hsyst->SetBinContent (ibinx, hstat->GetBinContent(ibinx));
    hsyst->SetBinError   (ibinx, hstat->GetBinContent(ibinx)*0.000001);
    hsystCorr->SetBinContent (ibinx, hstat->GetBinContent(ibinx));
    hsystCorr->SetBinError   (ibinx, hstat->GetBinContent(ibinx)*0.15);
  }
  *SaveAs_Title += "delete";
  texXtitle = texPtX;
  texYtitle = texDaughterPairsRecoEfficiency;
}



void Draw_Rapidity_MC(TH1D* &hstat, TH1D* &hsyst, TH1D* &hsystCorr, TFile *file_O2Analysis, TString* SaveAs_Title, TString* &texXtitle, TString* &texYtitle, int ipart) {
  //Looks within the signal area of the fit, 6Sigma

  // what is fed as input: TH3D* H3D_detectedV0s_DataLike = (TH3D*)file_O2Analysis->Get("lambdakzero-analysis-mc/h3dMass"+NamePart[ipart]);
  TH1D* H1D_rapidity = (TH1D*)file_O2Analysis->Get("lambdakzero-particle-count-mc/h"+NamePart[ipart]+"Count_Ydistrib");
  H1D_rapidity->Sumw2();
  TH1D* H1D_rapidity_rebinned = (TH1D*)H1D_rapidity->Rebin(5.,"H1D_rapidity_rebinned");
  H1D_rapidity_rebinned->Sumw2();

  TH1D* H1D_Ncoll_MC = (TH1D*)file_O2Analysis->Get("lambdakzero-particle-count-mc/hSelAndRecoMcCollCounter");
  float Ncoll_MC = H1D_Ncoll_MC->GetEntries();
  H1D_rapidity_rebinned->Scale(1./Ncoll_MC); //per collision

  TCanvas *canvas = new TCanvas ("canvas", NamehistoInvMass[ipart], 800, 600);

  TH1 *hFrameFits;

  canvas->cd(); //ATTENTION: +1 because canvas divide starts counting from 1
  hFrameFits = canvas->DrawFrame(-1.,0,1.,2.*H1D_rapidity_rebinned->GetMaximum());
  hFrameFits->SetXTitle("y");
  hFrameFits->SetYTitle("Count per collision");
  H1D_rapidity_rebinned->Draw("E,same");

  H1D_rapidity_rebinned->GetXaxis()->SetTitle(texInvMassDecays_titles[ipart]->Data());
  H1D_rapidity_rebinned->GetYaxis()->SetTitle("Counts");
  int icolor=0;
  H1D_rapidity_rebinned->SetMarkerStyle(markers[0]);
  H1D_rapidity_rebinned->SetMarkerColor(colors [icolor]);
  H1D_rapidity_rebinned->SetLineColor  (colors [icolor]);
  H1D_rapidity_rebinned->SetMarkerSize (0.5);

  // TLatex * text_pTbin = new TLatex (0.62,0.75,Form("pT [%.1f;%.1f] GeV", pTbins[iPt], pTbins[iPt+1]));
  // text_pTbin->SetNDC(kTRUE);
  // text_pTbin->Draw();
  // TLatex * textContext = new TLatex (0.18,0.82,"#bf{ALICE Performance}"); //BOLD
  // textContext->SetTextSize(0.05);
  // textContext->SetNDC(kTRUE); //remove if I want x,y in TLatex to be in the coordinate system of the histogram
  // textContext->Draw();
  // TLatex * textColl = new TLatex (0.18,0.75,"pp #sqrt{#it{s}} = 0.9 TeV, pilot beam 2021");
  TLatex * textColl = new TLatex (0.18,0.75,"pp #sqrt{#it{s}} = 0.9 TeV, Monte Carlo");
  textColl->SetTextSize(0.04);
  textColl->SetNDC(kTRUE); //remove if I want x,y in TLatex to be in the coordinate system of the histogram
  textColl->Draw();
  // TLatex * text_part = new TLatex (0.18,0.70,NamePart_Latex[ipart]+",   0 < #it{p}_{T} < 10 GeV/#it{c}");
  TLatex * text_part = new TLatex (0.18,0.70,NamePart_Latex[ipart]);    
  text_part->SetTextSize(0.04);
  text_part->SetNDC(kTRUE); //remove if I want x,y in TLatex to be in the coordinate system of the histogram
  text_part->Draw();
  // TLatex * text_fit = new TLatex (0.18,0.63,Form("Gaussian Fit:"));
  // text_fit->SetTextSize(0.035);
  // text_fit->SetNDC(kTRUE);
  // text_fit->Draw();
  // // TLatex * text_fitParam_mu = new TLatex (0.18,0.59,Form("#it{#mu} = %.0f #pm %.2f MeV/#it{c}^{2}", parGaussParab[iPt][1]*1000., parGaussParabError_mu[iPt]*1000.));
  // TLatex * text_fitParam_mu = new TLatex (0.18,0.59,Form("#it{#mu} = %.0f MeV/#it{c}^{2}", parGaussParab[iPt][1]*1000.));
  // text_fitParam_mu->SetTextSize(0.035);
  // text_fitParam_mu->SetNDC(kTRUE);
  // text_fitParam_mu->Draw();
  // // TLatex * text_fitParam_sigma = new TLatex (0.18,0.55,Form("#it{#sigma} = %.1f #pm %.2f MeV/#it{c}^{2}", parGaussParab[iPt][2]*1000., parGaussParabError_sigma[iPt]*1000.));
  // TLatex * text_fitParam_sigma = new TLatex (0.18,0.55,Form("#it{#sigma} = %.1f MeV/#it{c}^{2}", parGaussParab[iPt][2]*1000.));
  // text_fitParam_sigma->SetTextSize(0.035);
  // text_fitParam_sigma->SetNDC(kTRUE);
  // text_fitParam_sigma->Draw();


  canvas->SaveAs("Rapidity_MC_"+NamePart[ipart]+".pdf");

  /////////////////////////////////Does Nothing - only there to respect template/////////////////////////////////////////////
  hstat = (TH1D*)H1D_rapidity_rebinned->Clone("hstat");
  hstat->Sumw2();
  hstat->Reset("M");

  //////Error Bars///////
  hsyst     = (TH1D*)hstat->Clone("hsyst");
  hsyst->Reset("M");
  hsystCorr = (TH1D*)hstat->Clone("hsystCorr");
  hsystCorr->Reset("M");
  hstat     ->Sumw2(); 
  hsyst     ->Sumw2();
  hsystCorr ->Sumw2(); 
  int nbinx = hstat->GetNbinsX();
  
  for(int ibinx = 1; ibinx <= nbinx; ibinx++){
    hsyst->SetBinContent (ibinx, hstat->GetBinContent(ibinx));
    hsyst->SetBinError   (ibinx, hstat->GetBinContent(ibinx)*0.000001);
    hsystCorr->SetBinContent (ibinx, hstat->GetBinContent(ibinx));
    hsystCorr->SetBinError   (ibinx, hstat->GetBinContent(ibinx)*0.15);
  }
  *SaveAs_Title += "delete";
  texXtitle = texPtX;
  texYtitle = texDaughterPairsRecoEfficiency;
}




void PseudoEfficiency_O2vsHEP(TH1D* &hstat, TH1D* &hsyst, TH1D* &hsystCorr, int ipart, TString* SaveAs_Title, TString* &texXtitle, TString* &texYtitle, double* pTbins, int nbinpT) {

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

  TH3D* H3D_detectedV0s_DataLike = (TH3D*)file_O2Analysis->Get("lambdakzero-analysis/h3dMass"+NamePart[ipart]);
  // TH1I* H1I_SelectedEventCount = (TH1I*)file_O2Analysis->Get("lambdakzero-analysis/hSelectedEventCounter");
  TH1I* H1I_SelectedEventCount = (TH1I*)file_O2Analysis->Get("lambdakzero-analysis/hEventSelection");

  TH1D *hRawYield_vsPt;
  int InvMassRebinFactor = InvMassRebinFactor_standard[ipart];

  //Reference Corrected Yield
  double SignalMean, SignalStandardDev;
  cout << "Calculating Raw Yield for Ref" << endl;
  Get_RawYield_vsPt(hRawYield_vsPt, H3D_detectedV0s_DataLike, H1I_SelectedEventCount, file_O2Analysis, &SignalMean, &SignalStandardDev, ipart, pTbins, nbinpT, 0, -1, 0, SideBandSizeMultiplierModifier_array[0], SignalExtractionType_default, InvMassRebinFactor); ///real line that need to take the place of the one above once it's been written

  hstat = (TH1D*)hRawYield_vsPt->Clone("hstat");
  hstat->Reset("M");
  hstat->Divide(hRawYield_vsPt,hPtDiff_RawSpectrum_HEP_rebinned, 1., 1., "");

  /////////////////////////////////Does Nothing - only there to respect template/////////////////////////////////////////////

  //////Error Bars///////
  hsyst     = (TH1D*)hstat->Clone("hsyst");
  hsyst->Reset("M");
  hsystCorr = (TH1D*)hstat->Clone("hsystCorr");
  hsystCorr->Reset("M");
  hstat     ->Sumw2(); 
  hsyst     ->Sumw2();
  hsystCorr ->Sumw2(); 
  int nbinx = hstat->GetNbinsX();
  
  for(int ibinx = 1; ibinx <= nbinx; ibinx++){
    hsyst->SetBinContent (ibinx, hstat->GetBinContent(ibinx));
    hsyst->SetBinError   (ibinx, hstat->GetBinContent(ibinx)*0.000001);
    hsystCorr->SetBinContent (ibinx, hstat->GetBinContent(ibinx));
    hsystCorr->SetBinError   (ibinx, hstat->GetBinContent(ibinx)*0.15);
  }
  *SaveAs_Title += "PseudoEfficiency_HEPcomparison_withFit";
  texXtitle = texPtX;
  texYtitle = texPseudoEfficiency;
}


void Draw_pT_vs_pTmc(TH1D* &hstat, TH1D* &hsyst, TH1D* &hsystCorr, TFile *file_O2Analysis, TString* SaveAs_Title, TString* &texXtitle, TString* &texYtitle, int ipart) {
  //Looks within the signal area of the fit, 6Sigma

  // what is fed as input: TH3D* H3D_detectedV0s_DataLike = (TH3D*)file_O2Analysis->Get("lambdakzero-analysis-mc/h3dMass"+NamePart[ipart]);
  TH1D* H2D_pT = (TH1D*)file_O2Analysis->Get("lambdakzero-analysis-mc/hMC_vs_Data_pT_K0S");
  // TH1D* H2D_pT = (TH1D*)file_O2Analysis->Get("lambdakzero-qa/hMC_vs_Data_pT_allTracks");
  // H2D_pT->Sumw2();
  TH1D* H1D_pT_MCminusData = new TH1D("H1D_pT_MCminusData", "H1D_pT_MCminusData", 1000, -.5, .5);
  H1D_pT_MCminusData->Reset();

  double Delta_pT_mean = 0;
  long Total_Count = 0;
  for (int i = 0; i < H2D_pT->GetNbinsX(); i++) {
    for (int j = 0; j < H2D_pT->GetNbinsY(); j++) {
      // cout << "pT_MC = " << H2D_pT->GetYaxis()->GetBinCenter(j) << ", pT_data = " << H2D_pT->GetYaxis()->GetBinCenter(j) << endl;
      double pT_Delta = H2D_pT->GetYaxis()->GetBinCenter(j) - H2D_pT->GetXaxis()->GetBinCenter(i); // MC - data
      H1D_pT_MCminusData->Fill(pT_Delta, H2D_pT->GetBinContent(i,j));
      Delta_pT_mean = Delta_pT_mean + pT_Delta*H2D_pT->GetBinContent(i,j);
      Total_Count = Total_Count + H2D_pT->GetBinContent(i,j);
    }
  }
  cout << "Delta_pT_mean = " << Delta_pT_mean <<endl;
  cout << "Total_Count = " << Total_Count <<endl;
  Delta_pT_mean = Delta_pT_mean * 1./Total_Count;

  TCanvas *canvas = new TCanvas ("canvas", NamehistoInvMass[ipart], 800, 600);

  TH1 *hFrameFits;

  canvas->cd(); //ATTENTION: +1 because canvas divide starts counting from 1
  hFrameFits = canvas->DrawFrame(-0.3,0,0.3,1.4*H1D_pT_MCminusData->GetMaximum());
  hFrameFits->SetXTitle("pT_MC - pT_data");
  hFrameFits->SetYTitle("Count per collision");
  H1D_pT_MCminusData->Draw("same");

  int icolor=0;
  H1D_pT_MCminusData->SetMarkerStyle(markers[0]);
  H1D_pT_MCminusData->SetMarkerColor(colors [icolor]);
  H1D_pT_MCminusData->SetLineColor  (colors [icolor]);
  H1D_pT_MCminusData->SetMarkerSize (0.5);

  // TLatex * text_pTbin = new TLatex (0.62,0.75,Form("pT [%.1f;%.1f] GeV", pTbins[iPt], pTbins[iPt+1]));
  // text_pTbin->SetNDC(kTRUE);
  // text_pTbin->Draw();
  // TLatex * textContext = new TLatex (0.18,0.82,"#bf{ALICE Performance}"); //BOLD
  // textContext->SetTextSize(0.05);
  // textContext->SetNDC(kTRUE); //remove if I want x,y in TLatex to be in the coordinate system of the histogram
  // textContext->Draw();
  // TLatex * textColl = new TLatex (0.18,0.75,"pp #sqrt{#it{s}} = 0.9 TeV, pilot beam 2021");
  TLatex * textColl = new TLatex (0.18,0.75,"pp #sqrt{#it{s}} = 0.9 TeV, Monte Carlo");
  textColl->SetTextSize(0.04);
  textColl->SetNDC(kTRUE); //remove if I want x,y in TLatex to be in the coordinate system of the histogram
  textColl->Draw();
  // TLatex * text_part = new TLatex (0.18,0.70,NamePart_Latex[ipart]+",   0 < #it{p}_{T} < 10 GeV/#it{c}");
  TLatex * text_part = new TLatex (0.18,0.70,NamePart_Latex[ipart]);    
  text_part->SetTextSize(0.04);
  text_part->SetNDC(kTRUE); //remove if I want x,y in TLatex to be in the coordinate system of the histogram
  text_part->Draw();
  TLatex * text_fit = new TLatex (0.18,0.63,Form("#Delta#it{p}_{T} mean = %.2f MeV/#it{c}", Delta_pT_mean*1000.));
  text_fit->SetTextSize(0.035);
  text_fit->SetNDC(kTRUE);
  text_fit->Draw();
  // // TLatex * text_fitParam_mu = new TLatex (0.18,0.59,Form("#it{#mu} = %.0f #pm %.2f MeV/#it{c}^{2}", parGaussParab[iPt][1]*1000., parGaussParabError_mu[iPt]*1000.));
  // TLatex * text_fitParam_mu = new TLatex (0.18,0.59,Form("#it{#mu} = %.0f MeV/#it{c}^{2}", parGaussParab[iPt][1]*1000.));
  // text_fitParam_mu->SetTextSize(0.035);
  // text_fitParam_mu->SetNDC(kTRUE);
  // text_fitParam_mu->Draw();
  // // TLatex * text_fitParam_sigma = new TLatex (0.18,0.55,Form("#it{#sigma} = %.1f #pm %.2f MeV/#it{c}^{2}", parGaussParab[iPt][2]*1000., parGaussParabError_sigma[iPt]*1000.));
  // TLatex * text_fitParam_sigma = new TLatex (0.18,0.55,Form("#it{#sigma} = %.1f MeV/#it{c}^{2}", parGaussParab[iPt][2]*1000.));
  // text_fitParam_sigma->SetTextSize(0.035);
  // text_fitParam_sigma->SetNDC(kTRUE);
  // text_fitParam_sigma->Draw();


  canvas->SaveAs("pT_MC_minus_pT_data_"+NamePart[ipart]+".pdf");

  /////////////////////////////////Does Nothing - only there to respect template/////////////////////////////////////////////
  hstat = (TH1D*)H1D_pT_MCminusData->Clone("hstat");
  hstat->Sumw2();
  hstat->Reset("M");

  //////Error Bars///////
  hsyst     = (TH1D*)hstat->Clone("hsyst");
  hsyst->Reset("M");
  hsystCorr = (TH1D*)hstat->Clone("hsystCorr");
  hsystCorr->Reset("M");
  hstat     ->Sumw2(); 
  hsyst     ->Sumw2();
  hsystCorr ->Sumw2(); 
  int nbinx = hstat->GetNbinsX();
  
  for(int ibinx = 1; ibinx <= nbinx; ibinx++){
    hsyst->SetBinContent (ibinx, hstat->GetBinContent(ibinx));
    hsyst->SetBinError   (ibinx, hstat->GetBinContent(ibinx)*0.000001);
    hsystCorr->SetBinContent (ibinx, hstat->GetBinContent(ibinx));
    hsystCorr->SetBinError   (ibinx, hstat->GetBinContent(ibinx)*0.15);
  }
  *SaveAs_Title += "delete";
  texXtitle = texPtX;
  texYtitle = texDaughterPairsRecoEfficiency;
}


void Draw_FinderEfficiencies(TH1D* &hstat, TH1D* &hsyst, TH1D* &hsystCorr, TFile *file_O2Analysis_Run1Cuts, TFile *file_O2Analysis_SVCuts, TString* SaveAs_Title, TString* &texXtitle, TString* &texYtitle, int ipart, double* pTbins, int nbinpT){


  int InvMassRebinFactor = InvMassRebinFactor_standard[ipart];

  double SignalMean_Run1Cuts, SignalStandardDev_Run1Cuts;
  TH3D* H3D_detectedV0s_DataLike_Run1Cuts = (TH3D*)file_O2Analysis_Run1Cuts->Get("lambdakzero-analysis-mc/h3dMass"+NamePart[ipart]);
  TH3D* H3D_detectedV0s_TruePt_Run1Cuts = (TH3D*)file_O2Analysis_Run1Cuts->Get("lambdakzero-analysis-mc/h3dMass"+NamePart[ipart]+"_MC_truePt");
  TH1D* TrueV0PtSpectrum_Run1Cuts = (TH1D*)file_O2Analysis_Run1Cuts->Get("lambdakzero-particle-count-mc/h"+NamePart[ipart]+"Count_PtDiff");
  TH1I* H1I_SelectedEventCount_Run1Cuts = (TH1I*)file_O2Analysis_Run1Cuts->Get("lambdakzero-analysis-mc/hSelectedEventCounter");
  TH1D *hRawYield_vsPt_Run1Cuts;  
  TH1D *hMcEfficiency_vsPt_Run1Cuts;  
  Get_RawYield_vsPt_FeeddownCorrected(hRawYield_vsPt_Run1Cuts, H3D_detectedV0s_DataLike_Run1Cuts, H1I_SelectedEventCount_Run1Cuts, file_O2Analysis_Run1Cuts, &SignalMean_Run1Cuts, &SignalStandardDev_Run1Cuts, ipart, pTbins, nbinpT, 0, -1, 0, SideBandSizeMultiplierModifier_array[0], SignalExtractionType_default, InvMassRebinFactor);
  Get_McEfficiency_vsPt(hMcEfficiency_vsPt_Run1Cuts, H3D_detectedV0s_TruePt_Run1Cuts, TrueV0PtSpectrum_Run1Cuts, ipart, pTbins, nbinpT, 0, -1, SideBandSizeMultiplierModifier_array[0], SignalMean_Run1Cuts, SignalStandardDev_Run1Cuts, InvMassRebinFactor);

  double SignalMean_SVCuts, SignalStandardDev_SVCuts;
  TH3D* H3D_detectedV0s_DataLike_SVCuts = (TH3D*)file_O2Analysis_SVCuts->Get("lambdakzero-analysis-mc/h3dMass"+NamePart[ipart]);
  TH3D* H3D_detectedV0s_TruePt_SVCuts = (TH3D*)file_O2Analysis_SVCuts->Get("lambdakzero-analysis-mc/h3dMass"+NamePart[ipart]+"_MC_truePt");
  TH1D* TrueV0PtSpectrum_SVCuts = (TH1D*)file_O2Analysis_SVCuts->Get("lambdakzero-particle-count-mc/h"+NamePart[ipart]+"Count_PtDiff");
  TH1I* H1I_SelectedEventCount_SVCuts = (TH1I*)file_O2Analysis_SVCuts->Get("lambdakzero-analysis-mc/hSelectedEventCounter");
  TH1D *hRawYield_vsPt_SVCuts;  
  TH1D *hMcEfficiency_vsPt_SVCuts;  
  Get_RawYield_vsPt_FeeddownCorrected(hRawYield_vsPt_SVCuts, H3D_detectedV0s_DataLike_SVCuts, H1I_SelectedEventCount_SVCuts, file_O2Analysis_SVCuts, &SignalMean_SVCuts, &SignalStandardDev_SVCuts, ipart, pTbins, nbinpT, 0, -1, 0, SideBandSizeMultiplierModifier_array[0], SignalExtractionType_default, InvMassRebinFactor);
  Get_McEfficiency_vsPt(hMcEfficiency_vsPt_SVCuts, H3D_detectedV0s_TruePt_SVCuts, TrueV0PtSpectrum_SVCuts, ipart, pTbins, nbinpT, 0, -1, SideBandSizeMultiplierModifier_array[0], SignalMean_SVCuts, SignalStandardDev_SVCuts, InvMassRebinFactor);

  /////////// Draw Efficiency /////////////
  TString* pdfName0 = new TString("Finder_Efficiency_"+NamePart[ipart]);;
  // *pdfName0 += "Efficiency_"+NamePart[ipart];

  TH1* histograms_collection_0[2] = {hMcEfficiency_vsPt_Run1Cuts, hMcEfficiency_vsPt_SVCuts};
  TString* pdfLegend_collection_0[2] = {new TString("Run 1 cuts "+NamePart_Latex[ipart]), new TString("Run 3 cuts "+NamePart_Latex[ipart])};
  int ratioOption0 = 0;
  texXtitle = texPtX;
  texYtitle = texEfficiency;
  Draw_TH1_Histograms_Efficiency(histograms_collection_0, pdfLegend_collection_0, 2, pdfName0, texXtitle, texYtitle, pTbins, nbinpT, ratioOption0);





  /////////////////////////////////Does Nothing - only there to respect template/////////////////////////////////////////////
  hstat = (TH1D*)TrueV0PtSpectrum_SVCuts->Clone("hstat");
  hstat->Sumw2();
  hstat->Reset("M");

  //////Error Bars///////
  hsyst     = (TH1D*)hstat->Clone("hsyst");
  hsyst->Reset("M");
  hsystCorr = (TH1D*)hstat->Clone("hsystCorr");
  hsystCorr->Reset("M");
  hstat     ->Sumw2(); 
  hsyst     ->Sumw2();
  hsystCorr ->Sumw2(); 
  int nbinx = hstat->GetNbinsX();
  
  for(int ibinx = 1; ibinx <= nbinx; ibinx++){
    hsyst->SetBinContent (ibinx, hstat->GetBinContent(ibinx));
    hsyst->SetBinError   (ibinx, hstat->GetBinContent(ibinx)*0.000001);
    hsystCorr->SetBinContent (ibinx, hstat->GetBinContent(ibinx));
    hsystCorr->SetBinError   (ibinx, hstat->GetBinContent(ibinx)*0.15);
  }
  *SaveAs_Title += "delete";
  texXtitle = texPtX;
  texYtitle = texDaughterPairsRecoEfficiency;
}


void Draw_XiChargedToNeutralRatio(TH1D* &hstat, TH1D* &hsyst, TH1D* &hsystCorr, TFile *file_O2Analysis, TString* SaveAs_Title, TString* &texXtitle, TString* &texYtitle, int ipart, double* pTbins, int nbinpT){


  TH1D* H1D_XiCharged = (TH1D*)file_O2Analysis->Get("lambdakzero-particle-count-mc/hXsiMinusCount_PtDiff");
  TH1D* H1D_Xi0 = (TH1D*)file_O2Analysis->Get("lambdakzero-particle-count-mc/hXsi0Count_PtDiff");
  TH1I* H1I_SelectedEventCount = (TH1I*)file_O2Analysis->Get("lambdakzero-analysis-mc/hSelectedEventCounter");

  /////////// Draw Efficiency /////////////
  // TString* pdfName = new TString("Xi_Ratio_chargedToNeutral");;
  // *pdfName0 += "Efficiency_"+NamePart[ipart];

  // TH1* histograms_collection[2] = {H1D_XiCharged, H1D_Xi0};
  // TString* pdfLegend_collection[2] = {new TString("XiCharged"), new TString("Xi0")};
  // int ratioOption = 1; //to draw line at ratio=1
  // texXtitle = texPtX;
  // texYtitle = texRatio;
  // Draw_TH1_Histograms_withOneHistSystematics_ratio(histograms_collection, pdfLegend_collection, 2, pdfName, texXtitle, texYtitle, pTbins, nbinpT, ratioOption);

  //Xsi yield from pp 900GeV paper fit function (Xsi- + AntiXsi+)
  //dNdy = 0.0101 ± 0.0020
  //T = 175 ± 50
  //n = 5.2 ± 2.3
  cout << "coucou" << endl;
  float dNdy[2] = {3.05470e-04, 0.0020};
  float lTemp[2] = {1.00055e-01,50};
  float lPower[2] = {3.6,2.3};
  cout << "coucou" << endl;
  TF1* LevisTsalisFunction_XsiCharged_pp900GeV = new TF1("LevisTsalisFunction_XsiCharged_pp900GeV", myLevyPtXsi, 0, 3, 3);
  TF1* LevisTsalisFunction_XsiNeutral_pp900GeV = new TF1("LevisTsalisFunction_XsiNeutral_pp900GeV", myLevyPtXsi, 0, 3, 3);
  cout << "coucou" << endl;

  // Xi charged  
  int XiDistrib_RebinFactor = 5.0;
  TH1F *H1D_XiCharged_rebinned = (TH1F*)H1D_XiCharged->Rebin(XiDistrib_RebinFactor,"H1D_XiCharged_rebinned");
  cout << "coucou" << endl;
  H1D_XiCharged_rebinned->Sumw2();
  cout << "coucou" << endl;

  LevisTsalisFunction_XsiCharged_pp900GeV->SetParameter(0, dNdy[0]*H1I_SelectedEventCount->GetEntries());
  LevisTsalisFunction_XsiCharged_pp900GeV->SetParameter(1, lTemp[0]);
  LevisTsalisFunction_XsiCharged_pp900GeV->SetParameter(2, lPower[0]);

  TFitResultPtr fFitResult_charged = H1D_XiCharged_rebinned->Fit(LevisTsalisFunction_XsiCharged_pp900GeV,"S+L0R"); //


  // Xi neutral
  TH1F *H1D_Xi0_rebinned = (TH1F*)H1D_Xi0->Rebin(XiDistrib_RebinFactor,"H1D_Xi0_rebinned");
  H1D_Xi0_rebinned->Sumw2();

  LevisTsalisFunction_XsiNeutral_pp900GeV->SetParameter(0, dNdy[0]*H1I_SelectedEventCount->GetEntries());
  LevisTsalisFunction_XsiNeutral_pp900GeV->SetParameter(1, lTemp[0]);
  LevisTsalisFunction_XsiNeutral_pp900GeV->SetParameter(2, lPower[0]);

  TFitResultPtr fFitResult_neutral = H1D_Xi0_rebinned->Fit(LevisTsalisFunction_XsiNeutral_pp900GeV,"S+L0R"); //
  
  TH1D* ratio_XiChargedToNeutral = (TH1D*)H1D_XiCharged->Clone("ratio_XiChargedToNeutral");
  ratio_XiChargedToNeutral->Reset("M");
  int nbinpTxi = ratio_XiChargedToNeutral->GetNbinsX();
  for(int ibinPt = 1; ibinPt <= nbinpTxi; ibinPt++){
    float pT = ratio_XiChargedToNeutral->GetXaxis()->GetBinCenter(ibinPt);
    ratio_XiChargedToNeutral->SetBinContent(ibinPt,LevisTsalisFunction_XsiCharged_pp900GeV->Eval(pT,0,0)/LevisTsalisFunction_XsiNeutral_pp900GeV->Eval(pT,0,0));
  }
  // TF1* ratio_XiChargedToNeutral = new TF1("ratio_XiChargedToNeutral", "LevisTsalisFunction_XsiCharged_pp900GeV/LevisTsalisFunction_XsiNeutral_pp900GeV", 0, 3);

  TCanvas *canvasRatioXiChargedNeutral = new TCanvas ("canvasRatioXiChargedNeutral", NamehistoInvMass[ipart], 1200, 800);
  canvasRatioXiChargedNeutral->cd(0); //ATTENTION: +1 because canvas divide starts counting from 1
  ratio_XiChargedToNeutral->Draw();
  canvasRatioXiChargedNeutral->SaveAs("Ratio_Xi_chargedToNeutral_GenDistribs.pdf");



  /////////////////////////////////Does Nothing - only there to respect template/////////////////////////////////////////////
  hstat = (TH1D*)H1D_XiCharged->Clone("hstat");
  hstat->Sumw2();
  hstat->Reset("M");

  //////Error Bars///////
  hsyst     = (TH1D*)hstat->Clone("hsyst");
  hsyst->Reset("M");
  hsystCorr = (TH1D*)hstat->Clone("hsystCorr");
  hsystCorr->Reset("M");
  hstat     ->Sumw2(); 
  hsyst     ->Sumw2();
  hsystCorr ->Sumw2(); 
  int nbinx = hstat->GetNbinsX();
  
  for(int ibinx = 1; ibinx <= nbinx; ibinx++){
    hsyst->SetBinContent (ibinx, hstat->GetBinContent(ibinx));
    hsyst->SetBinError   (ibinx, hstat->GetBinContent(ibinx)*0.000001);
    hsystCorr->SetBinContent (ibinx, hstat->GetBinContent(ibinx));
    hsystCorr->SetBinError   (ibinx, hstat->GetBinContent(ibinx)*0.15);
  }
  *SaveAs_Title += "delete";
  texXtitle = texPtX;
  texYtitle = texDaughterPairsRecoEfficiency;
}

void Draw_Mass_Histograms_multipleRuns(TH1D* &hstat, TH1D* &hsyst, TH1D* &hsystCorr, TString* SaveAs_Title, TString* &texXtitle, TString* &texYtitle, int ipart) {
  //Looks within the signal area of the fit, 6Sigma; NOT SURE ZHAT I MEANT THERE 6 motnhs ago

  TH3D *H3D_DetectedV0s_DataLike[nRuns];
  TH3D *H3D_detectedV0s_DataLike_rebinnedZ[nRuns];
  TH1D* H1D_detectedV0s_DataLike_rebinnedZ_projectedZ[nRuns];

  // TCanvas *canvasMassvsPtFits = new TCanvas ("canvasMassPt_Runs", NamehistoInvMass[ipart], 800, 600);

  // TH1 *hFrameFits;

  for(int nRunsIterator = 0; nRunsIterator < nRuns; nRunsIterator++){
    cout << "------------- run = " << runList[nRunsIterator] << endl;
    H3D_DetectedV0s_DataLike[nRunsIterator] = (TH3D*)file_O2Analysis_RunComparison_array[nRunsIterator]->Get("lambdakzero-analysis-mc/h3dMass"+NamePart[ipart]);
    // H3D_detectedV0s_DataLike_rebinnedZ[nRunsIterator] = (TH3D*)H3D_DetectedV0s_DataLike[nRunsIterator]->RebinZ(InvMassRebinFactor_standard[ipart],"H3D_detectedV0s_DataLike_rebinnedZ");
    H3D_detectedV0s_DataLike_rebinnedZ[nRunsIterator] = (TH3D*)H3D_DetectedV0s_DataLike[nRunsIterator]->RebinZ(2.,"H3D_detectedV0s_DataLike_rebinnedZ");
    cout << "---------------------------------------- InvMassRebinFactor_standard = " << InvMassRebinFactor_standard[ipart] << endl;

    H1D_detectedV0s_DataLike_rebinnedZ_projectedZ[nRunsIterator] = (TH1D*)H3D_detectedV0s_DataLike_rebinnedZ[nRunsIterator]->ProjectionZ("Get_RawYield_vsPt_InvMass_"+NamePart[ipart]+Form("_run%i", runList[nRunsIterator]),0,-1,0,-1);
    H1D_detectedV0s_DataLike_rebinnedZ_projectedZ[nRunsIterator]->Scale(1./H1D_detectedV0s_DataLike_rebinnedZ_projectedZ[nRunsIterator]->GetMaximum()); //1./H1D_detectedV0s_DataLike_rebinnedZ_projectedZ[nRunsIterator]->GetEntries()*
    // hFrameFits = canvasMassvsPtFits->DrawFrame(LowerLimitDisplay_Mass[ipart],0,UpperLimitDisplay_Mass[ipart],1.45*H1D_detectedV0s_DataLike_rebinnedZ_projectedZ[nRunsIterator]->GetMaximum());
    // hFrameFits->SetXTitle(texInvMassDecays_titles[ipart]->Data());
    // hFrameFits->SetYTitle("Counts");
    // H1D_detectedV0s_DataLike_rebinnedZ_projectedZ[nRunsIterator]->Draw("E,same");

    // int icolor=0;
    // H1D_detectedV0s_DataLike_rebinnedZ_projectedZ[nRunsIterator]->SetMarkerStyle(markers[0]);
    // H1D_detectedV0s_DataLike_rebinnedZ_projectedZ[nRunsIterator]->SetMarkerColor(colors [nRunsIterator]);
    // H1D_detectedV0s_DataLike_rebinnedZ_projectedZ[nRunsIterator]->SetLineColor  (colors [nRunsIterator]);
    // H1D_detectedV0s_DataLike_rebinnedZ_projectedZ[nRunsIterator]->SetMarkerSize (0.5);
  }

  // TLatex * textColl = new TLatex (0.18,0.75,"pp #sqrt{#it{s}} = 0.9 TeV, Monte Carlo");
  // textColl->SetTextSize(0.04);
  // textColl->SetNDC(kTRUE); //remove if I want x,y in TLatex to be in the coordinate system of the histogram
  // textColl->Draw();
  // // TLatex * text_part = new TLatex (0.18,0.70,NamePart_Latex[ipart]+",   0 < #it{p}_{T} < 10 GeV/#it{c}");
  // TLatex * text_part = new TLatex (0.18,0.70,NamePart_Latex[ipart]);    
  // text_part->SetTextSize(0.04);
  // text_part->SetNDC(kTRUE); //remove if I want x,y in TLatex to be in the coordinate system of the histogram
  // text_part->Draw();
  // canvasMassvsPtFits->SaveAs("InvMass_"+NamePart[ipart]+"_RunComparison.pdf");

  texXtitle = texInvMass_K0ShortDecay_runs;
  texYtitle = texCountRelativeToMax;
  TString* pdfName = new TString("Inv_Mass_RunComparison_"+NamePart[ipart]+"_NormalisedToMax");;
  Draw_TH1_Histograms_Runs(H1D_detectedV0s_DataLike_rebinnedZ_projectedZ, nRuns, pdfName, texXtitle, texYtitle);


  // TH1* histograms_collection_2[2] = {Ratio_CorrectewdYield_to_MCgenYield_syst, Ratio_CorrectewdYield_to_MCgenYield};
  // Draw_TH1_Histograms_withOneHistSystematics_ratio(histograms_collection_2, pdfLegend_collection_2, 2, pdfName2, texXtitle, texYtitle, pTbins, nbinpT, ratioOption2);

  /////////////////////////////////Does Nothing - only there to respect template/////////////////////////////////////////////
  hstat = (TH1D*)H1D_detectedV0s_DataLike_rebinnedZ_projectedZ[0]->Clone("hstat");
  hstat->Sumw2();
  hstat->Reset("M");

  //////Error Bars///////
  hsyst     = (TH1D*)hstat->Clone("hsyst");
  hsyst->Reset("M");
  hsystCorr = (TH1D*)hstat->Clone("hsystCorr");
  hsystCorr->Reset("M");
  hstat     ->Sumw2(); 
  hsyst     ->Sumw2();
  hsystCorr ->Sumw2(); 
  int nbinx = hstat->GetNbinsX();
  
  for(int ibinx = 1; ibinx <= nbinx; ibinx++){
    hsyst->SetBinContent (ibinx, hstat->GetBinContent(ibinx));
    hsyst->SetBinError   (ibinx, hstat->GetBinContent(ibinx)*0.000001);
    hsystCorr->SetBinContent (ibinx, hstat->GetBinContent(ibinx));
    hsystCorr->SetBinError   (ibinx, hstat->GetBinContent(ibinx)*0.15);
  }
  *SaveAs_Title += "delete";
  texXtitle = texPtX;
  texYtitle = texDaughterPairsRecoEfficiency;
}

void Draw_TH1_Histograms_Runs(TH1D** histograms_collection, int collectionSize, TString* pdfName, TString* &texXtitle, TString* &texYtitle) {
  //draws only one hist with systematics, has to be the second one in the collection, and the associated systematics is taken from first hist in collection

  ////plots
  TCanvas *canvas = new TCanvas ("canvas"+*pdfName, NamehistoInvMass[ipart], 1200, 800);
  canvas->cd(0); //ATTENTION: +1 because canvas divide starts counting from 1
  TH1 *hFrame = canvas->DrawFrame(0.455,0,0.525,1.2*histograms_collection[0]->GetMaximum());
  hFrame->SetXTitle(texXtitle->Data());
  hFrame->SetYTitle(texYtitle->Data());
  TLegend * leg = new TLegend(0.7, 0.65, 0.87, 0.87);

  //draw histograms from collection, ignoring first one that is the systematics for second one
  for (int i = 0; i < collectionSize; i++) {
    // histograms_collection[i]->Draw("hist same p");
    histograms_collection[i]->Draw("same");
    histograms_collection[i]->SetMarkerStyle(markers_Runs[i]);
    histograms_collection[i]->SetMarkerColor(colors_Runs[i]);
    histograms_collection[i]->SetLineColor(colors_Runs[i]);

    leg->AddEntry(histograms_collection[i], runList_string[i],   "LP");
  }

  // if (ipart == 0) {
  //   //draw systematics for 2d hist using 1st hist in collection
  //   histograms_collection[0]->Draw("E2,same");
  //   histograms_collection[0]->SetMarkerStyle(markers[0]);
  //   histograms_collection[0]->SetMarkerColor(colors[0]);
  //   histograms_collection[0]->SetLineColor(colors[0]);
  //   histograms_collection[0]->SetFillStyle(0); // To draw empty boxes
  // }

  leg->SetTextSize(gStyle->GetTextSize()*0.5);
  if (collectionSize >= 2) {
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



double fpolysixthorder(double *x, double *par) {
  return par[0] + par[1]*x[0] + par[2]*x[0]*x[0] + par[3]*x[0]*x[0]*x[0] + par[4]*x[0]*x[0]*x[0]*x[0] + par[5]*x[0]*x[0]*x[0]*x[0]*x[0] + par[6]*x[0]*x[0]*x[0]*x[0]*x[0]*x[0];
}
double fsigmoid(double *x, double *par) {
  return par[0] / (1 + exp(-par[1] * (x[0]-par[2]))); //    y = L / (1 + exp(-k*(x-x0))) + b
}
double fgompertz(double *x, double *par) {
  return par[0] * exp(-par[1] * exp(-par[2]*x[0])); //    
}
double fgenerallogistic(double *x, double *par) {
  // return par[0] / std::pow((1 + par[1] * exp(-par[2]*x[0])),par[3]); //   works 
  return par[0] / std::pow((1 + par[1] * exp(-3*x[0])),10); //   test
}
  // parameter seeds that work for generalised logistic function
  // float a = 0.15; 
  // float b = 0;
  // float c = 1;
  // float d = 0.7;

void Draw_EfficiencyWithSigmoidFit(TH1D* &hstat, TH1D* &hsyst, TH1D* &hsystCorr, TFile *file_O2Analysis, TString* SaveAs_Title, TString* &texXtitle, TString* &texYtitle, int ipart, double* pTbins, int nbinpT){

  int InvMassRebinFactor = InvMassRebinFactor_standard[ipart];

  double SignalMean_SVCuts, SignalStandardDev_SVCuts;
  TH3D* H3D_detectedV0s_DataLike_SVCuts = (TH3D*)file_O2Analysis->Get("lambdakzero-analysis-mc/h3dMass"+NamePart[ipart]);
  TH3D* H3D_detectedV0s_TruePt_SVCuts = (TH3D*)file_O2Analysis->Get("lambdakzero-analysis-mc/h3dMass"+NamePart[ipart]+"_MC_truePt");
  TH1D* TrueV0PtSpectrum_SVCuts = (TH1D*)file_O2Analysis->Get("lambdakzero-particle-count-mc/h"+NamePart[ipart]+"Count_PtDiff");
  TH1I* H1I_SelectedEventCount_SVCuts = (TH1I*)file_O2Analysis->Get("lambdakzero-analysis-mc/hSelectedEventCounter"); // old 2022 files
  // TH1I* H1I_SelectedEventCount_SVCuts = (TH1I*)file_O2Analysis->Get("lambdakzero-analysis-mc/hEventSelection"); // new 2023 files
  TH1D *hRawYield_vsPt_SVCuts;  
  TH1D *hMcEfficiency_vsPt_SVCuts;  
  TH1D *hMcEfficiency_vsPt_SVCuts_withFit;  
  cout << "DEBUG 1" <<endl;

  Get_RawYield_vsPt_FeeddownCorrected(hRawYield_vsPt_SVCuts, H3D_detectedV0s_DataLike_SVCuts, H1I_SelectedEventCount_SVCuts, file_O2Analysis, &SignalMean_SVCuts, &SignalStandardDev_SVCuts, ipart, pTbins, nbinpT, 0, -1, 0, SideBandSizeMultiplierModifier_array[0], SignalExtractionType_default, InvMassRebinFactor);
  cout << "DEBUG 2" <<endl;
  Get_McEfficiency_vsPt(hMcEfficiency_vsPt_SVCuts, H3D_detectedV0s_TruePt_SVCuts, TrueV0PtSpectrum_SVCuts, ipart, pTbins, nbinpT, 0, -1, SideBandSizeMultiplierModifier_array[0], SignalMean_SVCuts, SignalStandardDev_SVCuts, InvMassRebinFactor);
  Get_McEfficiency_vsPt_alwaysFitSmooth(hMcEfficiency_vsPt_SVCuts_withFit, H3D_detectedV0s_TruePt_SVCuts, TrueV0PtSpectrum_SVCuts, ipart, pTbins, nbinpT, 0, -1, SideBandSizeMultiplierModifier_array[0], SignalMean_SVCuts, SignalStandardDev_SVCuts, InvMassRebinFactor);

  cout << "DEBUG 3" <<endl;

  /////////// Draw Efficiency /////////////
  TString* pdfName = new TString("Finder_Efficiency_"+NamePart[ipart]);;
  // *pdfName0 += "Efficiency_"+NamePart[ipart];

  TString* pdfLegend = new TString("Run 3 cuts "+NamePart_Latex[ipart]);
  int ratioOption0 = 0;
  texXtitle = texPtX;
  texYtitle = texEfficiency;

  // Draw_TH1_Histograms_Efficiency(histograms_collection_0, pdfLegend_collection_0, 2, pdfName0, texXtitle, texYtitle, pTbins, nbinpT, ratioOption0);
  // void Draw_TH1_Histograms_Efficiency(TH1** histograms_collection, TString** pdfLegend_collection, int collectionSize, TString* pdfName, TString* &texXtitle, TString* &texYtitle, double* bins, int nbins, int ratioOption) {
  float MaximumY = hMcEfficiency_vsPt_SVCuts->GetMaximum();

  cout << "DEBUG 4" <<endl;

  ////plots
  TCanvas *canvas = new TCanvas ("canvas"+*pdfName, NamehistoInvMass[ipart], 1200, 800);
  canvas->cd(0); //ATTENTION: +1 because canvas divide starts counting from 1
  // TH1 *hFrame = canvas->DrawFrame(bins[0],0,bins[nbins],1.4*MaximumY);
  TH1 *hFrame;

  hFrame = canvas->DrawFrame(0,0,pTbins[nbinpT],0.25);


  hFrame->SetXTitle(texXtitle->Data());
  hFrame->SetYTitle(texYtitle->Data());
  TLegend * leg = new TLegend(0.62, 0.77, 0.87, 0.87);
  cout << "DEBUG 5" <<endl;


  // ////fit with sigmoid
  // float b = 1;
  // float c = 1;
  // float d = 2.7;
  // cout << "DEBUG 7" <<endl;
  // TF1* sigmoidFit_Efficiency = new TF1("sigmoidFit_Efficiency", fsigmoid, 0, 3, 3);
  // sigmoidFit_Efficiency->SetParameter(0, hMcEfficiency_vsPt_SVCuts->GetMaximum()); //*H1I_SelectedEventCount_SVCuts->GetEntries()
  // sigmoidFit_Efficiency->SetParameter(1, b);
  // sigmoidFit_Efficiency->SetParameter(2, c);
  // // sigmoidFit_Efficiency->SetParameter(3, d);
  // // sigmoidFit_Efficiency->SetParameter(4, e);
  // sigmoidFit_Efficiency->SetParLimits(0, 0., 1.5*hMcEfficiency_vsPt_SVCuts->GetMaximum());
  // // sigmoidFit_Efficiency->SetParLimits(1, 0, 10);
  // // sigmoidFit_Efficiency->SetParLimits(2, 0, 10);
  // // sigmoidFit_Efficiency->SetParLimits(3, 0, 10);

  ////// fit with generalised logistics
  // float a = 1.7;
  float b = 1.2;
  float c = 0.46;
  float d = 2.7;
  float e = 20;
  cout << "DEBUG 7" <<endl;
  TF1* sigmoidFit_Efficiency = new TF1("sigmoidFit_Efficiency", fgenerallogistic, 0, 3, 2);
  sigmoidFit_Efficiency->SetParameter(0, hMcEfficiency_vsPt_SVCuts->GetMaximum()); //*H1I_SelectedEventCount_SVCuts->GetEntries()
  sigmoidFit_Efficiency->SetParameter(1, b);
  // sigmoidFit_Efficiency->SetParameter(2, c);
  // sigmoidFit_Efficiency->SetParameter(3, d);
  // sigmoidFit_Efficiency->SetParameter(4, e);
  sigmoidFit_Efficiency->SetParLimits(0, 0., 1.5*hMcEfficiency_vsPt_SVCuts->GetMaximum());
  // sigmoidFit_Efficiency->SetParLimits(1, 0, 10);
  // sigmoidFit_Efficiency->SetParLimits(2, 0, 10);
  // sigmoidFit_Efficiency->SetParLimits(3, 0, 10);

  cout << "DEBUG 8" <<endl;
  TFitResultPtr fFitResult_neutral = hMcEfficiency_vsPt_SVCuts->Fit(sigmoidFit_Efficiency,"ME+R0"); // the option G helps avoid the "Warning in <Fit>: Abnormal termination of minimization" warning and combined with ME gives "ERROR MATRIX ACCURATE " instead of "ERROR MATRIX UNCERTAINTY  24.6 per cent" or "ERROR NOT POS. DEFINED"
  // fFitResult_neutral = hMcEfficiency_vsPt_SVCuts->Fit(sigmoidFit_Efficiency,"S+R0"); // the option G helps avoid the "Warning in <Fit>: Abnormal termination of minimization" warning and combined with ME gives "ERROR MATRIX ACCURATE " instead of "ERROR MATRIX UNCERTAINTY  24.6 per cent" or "ERROR NOT POS. DEFINED"
  // gotta maybe calculate errors manually, looking at a plot of the fit function with the errors, it doesnt seem as bit as what the drawn function shows  

  sigmoidFit_Efficiency->SetLineColor(kRed);
  // sigmoidFit_Efficiency->SetFillColor(kRed);
  sigmoidFit_Efficiency->Draw("same");
  
  //draw histograms from collection
  int i = 0;
  hMcEfficiency_vsPt_SVCuts->Draw("E,same");
  hMcEfficiency_vsPt_SVCuts->SetMarkerStyle(markers[i]);
  hMcEfficiency_vsPt_SVCuts->SetMarkerColor(colors[i]);
  hMcEfficiency_vsPt_SVCuts->SetLineColor(colors[i]);
  hMcEfficiency_vsPt_SVCuts_withFit->Draw("E,same");
  hMcEfficiency_vsPt_SVCuts_withFit->SetMarkerStyle(markers[i+1]);
  hMcEfficiency_vsPt_SVCuts_withFit->SetMarkerColor(colors[i+1]);
  hMcEfficiency_vsPt_SVCuts_withFit->SetLineColor(colors[i+1]);

  leg->AddEntry(hMcEfficiency_vsPt_SVCuts, *pdfLegend,   "LP");

  leg->SetTextSize(gStyle->GetTextSize()*0.6);
  leg->Draw("same");

  cout << "DEBUG 6" <<endl;

  TLatex * textColl = new TLatex (0.18,0.82,"pp #sqrt{#it{s}} = 900 GeV, pilot beam 2021 MC");
  textColl->SetTextSize(0.04);
  textColl->SetNDC(kTRUE); //remove if I want x,y in TLatex to be in the coordinate system of the histogram
  textColl->Draw();
  TLatex * text_part = new TLatex (0.18,0.75,NamePart_Latex[ipart]);
  text_part->SetTextSize(0.04);
  text_part->SetNDC(kTRUE); //remove if I want x,y in TLatex to be in the coordinate system of the histogram
  text_part->Draw();


  cout << "DEBUG 9" <<endl;

  canvas->SaveAs(*pdfName+".pdf");


  cout << "DEBUG 10" <<endl;

  /////////////////////////////////Does Nothing - only there to respect template/////////////////////////////////////////////
  hstat = (TH1D*)TrueV0PtSpectrum_SVCuts->Clone("hstat");
  hstat->Sumw2();
  hstat->Reset("M");

  //////Error Bars///////
  hsyst     = (TH1D*)hstat->Clone("hsyst");
  hsyst->Reset("M");
  hsystCorr = (TH1D*)hstat->Clone("hsystCorr");
  hsystCorr->Reset("M");
  hstat     ->Sumw2(); 
  hsyst     ->Sumw2();
  hsystCorr ->Sumw2(); 
  int nbinx = hstat->GetNbinsX();
  
  for(int ibinx = 1; ibinx <= nbinx; ibinx++){
    hsyst->SetBinContent (ibinx, hstat->GetBinContent(ibinx));
    hsyst->SetBinError   (ibinx, hstat->GetBinContent(ibinx)*0.000001);
    hsystCorr->SetBinContent (ibinx, hstat->GetBinContent(ibinx));
    hsystCorr->SetBinError   (ibinx, hstat->GetBinContent(ibinx)*0.15);
  }
  *SaveAs_Title += "delete";
  texXtitle = texPtX;
  texYtitle = texDaughterPairsRecoEfficiency;
}


void Draw_PtResponseMatrix(TH1D* &hstat, TH1D* &hsyst, TH1D* &hsystCorr, TFile *file_O2Analysis, TString* SaveAs_Title, TString* &texXtitle, TString* &texYtitle, double* pTbins, int nbinpT) {
  TH2D* H2D_PtDataVsMC = (TH2D*)file_O2Analysis_unfoldingInfo->Get("thesis-plots-cut-distribs-signal-bkg/hPtDataVsMC_"+NamePart[ipart]);
  // TH2D* H2D_PtDataVsMC = (TH2D*)file_O2Analysis->Get("thesis-plots-cut-distribs-signal-bkg/hPtDataVsMC");

  TH2D* H2D_PtDataVsMC_rebinned = new TH2D("H2D_PtDataVsMC_rebinned", "H2D_PtDataVsMC_rebinned", nbinpT, pTbins, nbinpT, pTbins);

  TH2D* H2D_PtDataVsMC_rebinned_density = new TH2D("H2D_PtDataVsMC_rebinned", "H2D_PtDataVsMC_rebinned", nbinpT, pTbins, nbinpT, pTbins);


  cout << "______Rebinning_H2D_PtDataVsMC:" << endl;
  for(int i = 1; i <= nbinpT; i++){
    int arrayID_i = i - 1;
    int LowerBound_bin_i = H2D_PtDataVsMC->GetXaxis()->FindBin(pTbins[arrayID_i]);
    int UpperBound_bin_i = H2D_PtDataVsMC->GetXaxis()->FindBin(pTbins[arrayID_i+1]) - 1;
    cout << "LowerBound bin_i " << i << " = " << LowerBound_bin_i << ", UpperBound bin_i " << i << " = " << UpperBound_bin_i << endl;
    for(int j = 1; j <= nbinpT; j++){
      int arrayID_j = j - 1;
      int LowerBound_bin_j = H2D_PtDataVsMC->GetYaxis()->FindBin(pTbins[arrayID_j]);
      int UpperBound_bin_j = H2D_PtDataVsMC->GetYaxis()->FindBin(pTbins[arrayID_j+1]) - 1;
      cout << "LowerBound bin_j " << j << " = " << LowerBound_bin_j << ", UpperBound bin_j " << j << " = " << UpperBound_bin_j << endl;
      H2D_PtDataVsMC_rebinned->SetBinContent(i, j, H2D_PtDataVsMC->Integral(LowerBound_bin_i,UpperBound_bin_i,LowerBound_bin_j,UpperBound_bin_j));
    }
  }  
  H2D_PtDataVsMC_rebinned->Sumw2();
  for(int i = 1; i <= nbinpT; i++){
    for(int j = 1; j <= nbinpT; j++){
      double A_pT = H2D_PtDataVsMC_rebinned->GetXaxis()->GetBinWidth(i) * H2D_PtDataVsMC_rebinned->GetYaxis()->GetBinWidth(j);
      cout << "H2D_PtDataVsMC_rebinned(" << i << "," << j << ") = " << H2D_PtDataVsMC_rebinned->GetBinContent(i, j) << ", A_pT = " << A_pT << endl;
      H2D_PtDataVsMC_rebinned_density->SetBinContent(i, j, 1./A_pT * 1./H2D_PtDataVsMC_rebinned->GetEntries() * H2D_PtDataVsMC_rebinned->GetBinContent(i, j));
      H2D_PtDataVsMC_rebinned_density->SetBinError(i, j, sqrt(1./A_pT * 1./H2D_PtDataVsMC_rebinned->GetEntries()) * H2D_PtDataVsMC_rebinned->GetBinError(i, j));
    }
  }
  TCanvas *canvas = new TCanvas ("canvas", NamehistoInvMass[ipart], 800, 600);
  canvas->cd();

  H2D_PtDataVsMC_rebinned_density->SetXTitle(texPtMC->Data());
  H2D_PtDataVsMC_rebinned_density->SetYTitle(texPtMeasured->Data());
  
  H2D_PtDataVsMC_rebinned_density->Draw("COLZ");
  canvas->SetRightMargin(0.15);
  canvas->SaveAs("PtResponseMatrix_true"+NamePart[ipart]+".pdf");

  /////////////////////////////////Does Nothing - only there to respect template/////////////////////////////////////////////
  hstat = (TH1D*)H2D_PtDataVsMC->ProjectionX("placeholder", 0, -1)->Clone("hstat");
  hstat->Sumw2();
  hstat->Reset("M");

  //////Error Bars///////
  hsyst     = (TH1D*)hstat->Clone("hsyst");
  hsyst->Reset("M");
  hsystCorr = (TH1D*)hstat->Clone("hsystCorr");
  hsystCorr->Reset("M");
  hstat     ->Sumw2(); 
  hsyst     ->Sumw2();
  hsystCorr ->Sumw2(); 
  int nbinx = hstat->GetNbinsX();
  
  for(int ibinx = 1; ibinx <= nbinx; ibinx++){
    hsyst->SetBinContent (ibinx, hstat->GetBinContent(ibinx));
    hsyst->SetBinError   (ibinx, hstat->GetBinContent(ibinx)*0.000001);
    hsystCorr->SetBinContent (ibinx, hstat->GetBinContent(ibinx));
    hsystCorr->SetBinError   (ibinx, hstat->GetBinContent(ibinx)*0.15);
  }
  *SaveAs_Title += "delete";
  texXtitle = texPtX;
  texYtitle = texDaughterPairsRecoEfficiency;
}



void Get_PtResponseMatrix(TH2D * H2D_PtResponseMatrix, TFile *file_O2Analysis, double* pTbins, int nbinpT){
  TH2D* H2D_PtDataVsMC = (TH2D*)file_O2Analysis_unfoldingInfo->Get("thesis-plots-cut-distribs-signal-bkg/hPtDataVsMC_"+NamePart[ipart]);




  cout << "______Rebinning_H2D_PtDataVsMC:" << endl;
  for(int i = 1; i <= nbinpT; i++){
    int arrayID_i = i - 1;
    int LowerBound_bin_i = H2D_PtDataVsMC->GetXaxis()->FindBin(pTbins[arrayID_i]);
    int UpperBound_bin_i = H2D_PtDataVsMC->GetXaxis()->FindBin(pTbins[arrayID_i+1]) - 1;
    cout << "LowerBound bin_i " << i << " = " << LowerBound_bin_i << ", UpperBound bin_i " << i << " = " << UpperBound_bin_i << endl;
    for(int j = 1; j <= nbinpT; j++){
      int arrayID_j = j - 1;
      int LowerBound_bin_j = H2D_PtDataVsMC->GetYaxis()->FindBin(pTbins[arrayID_j]);
      int UpperBound_bin_j = H2D_PtDataVsMC->GetYaxis()->FindBin(pTbins[arrayID_j+1]) - 1;
      cout << "LowerBound bin_j " << j << " = " << LowerBound_bin_j << ", UpperBound bin_j " << j << " = " << UpperBound_bin_j << endl;
      H2D_PtResponseMatrix->SetBinContent(i, j, H2D_PtDataVsMC->Integral(LowerBound_bin_i,UpperBound_bin_i,LowerBound_bin_j,UpperBound_bin_j));
      cout << "PtResponseMatrix(" << i << ", " << j << ") = " << j << " = " << H2D_PtResponseMatrix->GetBinContent(i, j) << endl;
    }
  }  
  for(int i = 1; i <= nbinpT; i++){
    for(int j = 1; j <= nbinpT; j++){
      cout << "PtResponseMatrix(" << i << ", " << j << ") = " << j << " = " << H2D_PtResponseMatrix->GetBinContent(i, j) << endl;
    }
  }
  H2D_PtResponseMatrix->Sumw2();
  // // Density matrix
  // for(int i = 1; i <= nbinpT; i++){
  //   for(int j = 1; j <= nbinpT; j++){
  //     double A_pT = H2D_PtResponseMatrix->GetXaxis()->GetBinWidth(i) * H2D_PtResponseMatrix->GetYaxis()->GetBinWidth(j);
  //     cout << "H2D_PtResponseMatrix(" << i << "," << j << ") = " << H2D_PtResponseMatrix->GetBinContent(i, j) << ", A_pT = " << A_pT << endl;
  //     H2D_PtResponseMatrix_density->SetBinContent(i, j, 1./A_pT * 1./H2D_PtResponseMatrix->GetEntries() * H2D_PtResponseMatrix->GetBinContent(i, j));
  //     H2D_PtResponseMatrix_density->SetBinError(i, j, sqrt(1./A_pT * 1./H2D_PtResponseMatrix->GetEntries()) * H2D_PtResponseMatrix->GetBinError(i, j));
  //   }
  // }
}




void Get_PtResponseMatrix_Density(TH2D * H2D_PtResponseMatrix, TFile *file_O2Analysis, double* pTbins, int nbinpT){
  TH2D* H2D_PtDataVsMC = (TH2D*)file_O2Analysis_unfoldingInfo->Get("thesis-plots-cut-distribs-signal-bkg/hPtDataVsMC_"+NamePart[ipart]);

  TH1I* H1I_SelectedEventCount = (TH1I*)file_O2Analysis->Get("lambdakzero-analysis/hSelectedEventCounter");
  float SelectedEventCount = H1I_SelectedEventCount->GetEntries();



  cout << "______Rebinning_H2D_PtDataVsMC:" << endl;
  for(int i = 1; i <= nbinpT; i++){
    int arrayID_i = i - 1;
    int LowerBound_bin_i = H2D_PtDataVsMC->GetXaxis()->FindBin(pTbins[arrayID_i]);
    int UpperBound_bin_i = H2D_PtDataVsMC->GetXaxis()->FindBin(pTbins[arrayID_i+1]) - 1;
    cout << "LowerBound bin_i " << i << " = " << LowerBound_bin_i << ", UpperBound bin_i " << i << " = " << UpperBound_bin_i << endl;
    for(int j = 1; j <= nbinpT; j++){
      int arrayID_j = j - 1;
      int LowerBound_bin_j = H2D_PtDataVsMC->GetYaxis()->FindBin(pTbins[arrayID_j]);
      int UpperBound_bin_j = H2D_PtDataVsMC->GetYaxis()->FindBin(pTbins[arrayID_j+1]) - 1;
      cout << "LowerBound bin_j " << j << " = " << LowerBound_bin_j << ", UpperBound bin_j " << j << " = " << UpperBound_bin_j << endl;
      H2D_PtResponseMatrix->SetBinContent(i, j, H2D_PtDataVsMC->Integral(LowerBound_bin_i,UpperBound_bin_i,LowerBound_bin_j,UpperBound_bin_j));
      cout << "PtResponseMatrix(" << i << ", " << j << ") = " << j << " = " << H2D_PtResponseMatrix->GetBinContent(i, j) << endl;
    }
  }  
  for(int i = 1; i <= nbinpT; i++){
    for(int j = 1; j <= nbinpT; j++){
      cout << "PtResponseMatrix(" << i << ", " << j << ") = " << j << " = " << H2D_PtResponseMatrix->GetBinContent(i, j) << endl;
    }
  }
  H2D_PtResponseMatrix->Sumw2();

  // Transform into density matrix
  for(int i = 1; i <= nbinpT; i++){
    for(int j = 1; j <= nbinpT; j++){
      double A_pT = H2D_PtResponseMatrix->GetXaxis()->GetBinWidth(i) * H2D_PtResponseMatrix->GetYaxis()->GetBinWidth(j);
      cout << "H2D_PtResponseMatrix(" << i << "," << j << ") = " << H2D_PtResponseMatrix->GetBinContent(i, j) << ", A_pT = " << A_pT << endl;
      H2D_PtResponseMatrix->SetBinContent(i, j, 1./A_pT * 1./H2D_PtResponseMatrix->GetEntries() * H2D_PtResponseMatrix->GetBinContent(i, j));
      H2D_PtResponseMatrix->SetBinError(i, j, sqrt(1./A_pT * 1./H2D_PtResponseMatrix->GetEntries()) * H2D_PtResponseMatrix->GetBinError(i, j));
    }
  }
  H2D_PtResponseMatrix->Scale(1./SelectedEventCount);

}

