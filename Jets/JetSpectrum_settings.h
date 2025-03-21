#ifndef JETSPECTRUM_SETTINGS_H
#define JETSPECTRUM_SETTINGS_H

// Analysis settings
const int nJetType = 3;
const TString jetType[nJetType] = {"charged", "neutral", "full"};
const int nJetLevel = 3;
const TString jetLevel[nJetLevel] = {"data", "mcd", "mcp"};
const int nRadius = 3;
const TString RadiusLegend[nRadius] = {"R = 0.2", "R = 0.4", "R = 0.6"};
const float arrayRadius[nRadius] = {0.2, 0.4, 0.6};
const int nRandomConeTypes = 5;
const TString randomConeTypeList[nRandomConeTypes] = {"", "withoutleadingjet", "randomtrackdirection", "randomtrackdirectionwithoutoneleadingjets", "randomtrackdirectionwithouttwoleadingjets"};

const double etaAnalysisRange[2] = {-0.5, 0.5};
// Choice of jet type (charged, neutral, full) and level (data, detector level, particle level)
const int iJetType = 0;
const int iJetLevel = 0;
const int randomConeType = 1;

TFile* file_O2Analysis_run2ComparisonFile = new TFile("Datasets/Run2_Unfolding_AreaBased_HannahMethod_R020_Nominal_ExtendedPtRange/Unfolding_AreaBased_HannahMethod_R020_Nominal_ExtendedPtRange.root");


const float centralityRange[2] = {0, 10};

////////////////////////////////////////////////
////////// Unfolding settings - start //////////
////////////////////////////////////////////////



char mergingPrior[] = "noPriorMerging";     // prior options: mcpPriorMerging, mcdPriorMerging, measuredPriorMerging, noPriorMerging, testAliPhysics
char unfoldingPrior[] = "mcpPriorUnfolding";     // prior options: mcpPriorUnfolding, mcdPriorUnfolding, measuredPriorUnfolding, noPriorUnfolding, testAliPhysics /////// if using mcp as prior, should have the leading track cut like data
const bool doYSliceNormToOneDetResp = true; // should be true (done by marta)
const bool doYSliceNormToOneCombinedResp = false; // should be false (not done by marta); breaks unfolding with svd if true
const bool doUnfoldingPriorDivision = false; // unfolding doesn't work anymore if this is done, gives almost flat pT distribution, though refolding test is good; I am not sure why, might be because roounfold already does that; one good reason to avoid it anyway, is that roounfold already seems to deal with errors; my error propagation doesn't take into account off-diagonal covariance elements, and so can only be worse 
const bool scaleRespByWidth = false; 
const bool matrixTransformationOrder = 0; //0: reweight with unfoldingPrior, then rebin with merging prior, then do YSliceNorm and scaleRespByWidth if set to true (0 seems to work well); 1: rebin, then YSliceNorm and scaleRespByWidth, then reweight; 2: rebin, then reweight, then YSliceNorm and scaleRespByWidth

char unfoldingMethod[] = "Bayes"; // unfolding method options: Bayes, Svd
char optionsAnalysis[100] = "";

const bool isDataPbPb = true; // if false -> pp
const bool doBkgSubtractionInData = true;
const bool doBkgSubtractionInMC = false;
const bool useFactorisedMatrix = true; // use factorised response matrix for unfolding, or not
const bool mcIsWeighted = true; // use if the MC has been weighted to have more high pt jets?
int applyEfficiencies = 3; // for test purposes: 0: no efficiency correction, 1: kine only, 2: jet finding efficiency only, 3: both active; only applied if useManualRespMatrixSettingMethod is true
bool applyFakes = true; // only applied if useManualRespMatrixSettingMethod is true; 18/03: if false?
const bool useFineBinningTest = false; //looks like this gives the same flat distrib as when using coarse binning: so rebinning isnt the issue; need to change finBinning back to start at 0 when I dont use this

const bool doWidthScalingEarly = false;                         // to avoid pT bin width having an influence on spectrum; which one should be done? early or end? for now will be done at end
const bool doWidthScalingAtEnd = true;                          //


// all three below should probably be true;
// but then it breaks svd convergence! find out why;
const bool normDetRespByNEvts = false; //that's what breaks svd; https://arxiv.org/pdf/hep-ph/9509307 seems to say one should use the number of events matrix (see last paragraph of conclusion) instead of a probability matrix, to further reduce errors
const bool normGenAndMeasByNEvts = false;

const bool normaliseDistribsAfterUnfolding = true;   //both normaliseDistribsAfterUnfolding and normaliseDistribsBeforeUnfolding should be the same, else refolding test fails; without the counts are 1Ei, with they are 1E-j, so should be set to true
const bool normaliseDistribsBeforeUnfolding = true;   //both normaliseDistribsAfterUnfolding and normaliseDistribsBeforeUnfolding should be the same, else refolding test fails; without the counts are 1Ei, with they are 1E-j, so should be set to true

const bool useManualRespMatrixSettingMethod = true; // 18/03 remark: if set to true, both refold methods are constistent; if set to false, the roounfold one is identical as if true, but the manual one becomes different and bad
const bool normaliseRespYSliceForRefold = true; // ??????? THAT IS APPARENTLY REQUIRED TO REFOLD MANUALLY! even though the initial resp matrix used for the unfolding isn't normalised like this

bool controlMC = false; // use file_O2Analysis_ppSimDetectorEffect_unfoldingControl MC file as input to unfolding (with h_jet_pt_rhoareasubtracted distrib on file), rather than real data, and as comparison to gen (with h_jet_pt_part distrib on file); weighted control MC, and control for PbPb are not yet implemented
bool comparePbPbWithRun2 = false; // if isDataPbPb == true, then do the comparison with file_O2Analysis_run2ComparisonFile (Nevents for this is hardcoded to what Laura told me: see mattermost discussion)

bool automaticBestSvdParameter = false;

const bool drawIntermediateResponseMatrices = false;

// 18/03 remarks:
// - applyFakes and applyEfficiencies:
//      checked if useManualRespMatrixSettingMethod is true, but the non-initial response method does take fake into account; 
//      maybe rename applyFakes/EfficienciesIfInitialResponseMethod; 
//      also, is also used if useManualRespMatrixSettingMethod is false in Get_Pt_spectrum_dataUnfoldedThenRefolded_RooUnfoldMethod_preWidthScalingAndEvtNorm --> UNDERSTAND WHY
// - useManualRespMatrixSettingMethod if false, manual and roounfold refolding not consistent; what seems to be failing is the manual (manual multiplication of matrices) refolding in this case 
//   (though it works fine if useManualRespMatrixSettingMethod is true); if useManualRespMatrixSettingMethod is false, then the response matrix is filled weirdly (joonsuk Unfolding) with roounfold methods: RooUnfoldResponse->Fill(), Fake() and Miss()
// -----> for now only use useManualRespMatrixSettingMethod = true; if time, find out what's wrong in refolding in other method
// - manual refold has error issues; if I use "Joonsuk binning for pp" it's fine, but if I use "PbPb Aimeric old" then errors explode in the refolding plot; 
//   and it looks like bins 8 and 9 have weirdly immense errors (not the case with joonsuk binning)
//   differences between two binnings: 1) joonsuk has same bin count for gen and rec, mine has less for gen; 2) joonsuk rec starts at same value as gen, I start rec above gen
//   "PbPb Aimeric default" works much better;


float ptWindowDisplay[2] = {10, 140};
std::array<std::array<float, 2>, 2> drawnWindowUnfoldedMeasurement = {{{ptWindowDisplay[0], ptWindowDisplay[1]}, {-999, -999}}}; // {{xmin, xmax}, {ymin, ymax}}


// char mergingPrior[] = "noPriorMerging";     // prior options: mcpPriorMerging, mcdPriorMerging, measuredPriorMerging, noPriorMerging, testAliPhysics
// char unfoldingPrior[] = "noPriorUnfolding";     // prior options: mcpPriorUnfolding, mcdPriorUnfolding, measuredPriorUnfolding, noPriorUnfolding, testAliPhysics /////// if using mcp as prior, should have the leading track cut like data
// const bool doYSliceNormToOneDetResp = false; //SHOULD BE TRUE IF USING PRIOR
// char unfoldingMethod[] = "Bayes"; // unfolding method options: Bayes, Svd
// char optionsAnalysis[100] = "";

// const bool isDataPbPb = true; // if false -> pp
// const bool doBkgSubtractionInData = true;
// const bool doBkgSubtractionInMC = false;
// const bool useFactorisedMatrix = true; // use factorised response matrix for unfolding, or not
// const bool mcIsWeighted = true; // use if the MC has been weighted to have more high pt jets?
// int applyEfficiencies = 3; // for test purposes: 0: no efficiency correction, 1: kine only, 2: jet finding efficiency only, 3: both active; only applied if useManualRespMatrixSettingMethod is true
// bool applyFakes = true; // only applied if useManualRespMatrixSettingMethod is true
// const bool useFineBinningTest = true; //looks like this gives the same flat distrib as when using coarse binning: so rebinning isnt the issue; need to change finBinning back to start at 0 when I dont use this

// const bool scaleRespByWidth = false;
// const bool doWidthScalingEarly = false;                          //  doesn't seem to have any effect, so I can probably use it: doesn't change the ratios (at least measured/unfolded and mcp/unfolded, haven't checked folded/unfolded)
// const bool doWidthScalingAtEnd = true;                          //  doesn't seem to have any effect, so I can probably use it: doesn't change the ratios (at least measured/unfolded and mcp/unfolded, haven't checked folded/unfolded)


// // all three below should probably be true;
// // but then it breaks svd convergence! find out why;
// const bool normDetRespByNEvts = false; //that's what breaks svd; https://arxiv.org/pdf/hep-ph/9509307 seems to say one should use the number of events matrix (see last paragraph of conclusion) instead of a probability matrix
// const bool normGenAndMeasByNEvts = false;

// const bool normaliseDistribsAfterUnfolding = true;   //both normaliseDistribsAfterUnfolding and normaliseDistribsBeforeUnfolding should be the same, else refolding test fails; without the counts are 1Ei, with they are 1E-j, so should be set to true
// const bool normaliseDistribsBeforeUnfolding = true;   //both normaliseDistribsAfterUnfolding and normaliseDistribsBeforeUnfolding should be the same, else refolding test fails; without the counts are 1Ei, with they are 1E-j, so should be set to true

// const bool useManualRespMatrixSettingMethod = false; // discrepancy false vs true here seems to be that I do not model fakes in my initial method
// const bool normaliseRespYSliceForRefold = true; // ??????? THAT IS APPARENTLY REQUIRED TO REFOLD MANUALLY! even though the initial resp matrix used for the unfolding isn't normalised like this

// bool controlMC = false; // use file_O2Analysis_ppSimDetectorEffect_unfoldingControl MC file as input to unfolding (with h_jet_pt_rhoareasubtracted distrib on file), rather than real data, and as comparison to gen (with h_jet_pt_part distrib on file); weighted control MC, and control for PbPb are not yet implemented
// bool comparePbPbWithRun2 = false; // if isDataPbPb == true, then do the comparison with file_O2Analysis_run2ComparisonFile (Nevents for this is hardcoded to what Laura told me: see mattermost discussion)

// const bool drawIntermediateResponseMatrices = false;


// float ptWindowDisplay[2] = {5, 200};
// std::array<std::array<float, 2>, 2> drawnWindow = {{{ptWindowDisplay[0], ptWindowDisplay[1]}, {-999, -999}}}; // {{xmin, xmax}, {ymin, ymax}}


////////////////////////////////////////////////
////////// Unfolding settings - end ////////////
////////////////////////////////////////////////

// Lessons:
// - if the rec distrib is cut at some 5GeV or something else, using the prior after normYslice won't work well if I don't set the binning to start AFTER this cut, not before! it works perfectly if I start at the cut

////////////////////////////////////////////////
//////////////// pt binning options ////////////
////////////////////////////////////////////////


// // pT binning for jets - gen = rec
// double ptBinsJetsRec[nRadius][30] = {{0.0, 5., 10., 15., 20., 25., 30., 40., 50., 60., 70., 80., 100., 120., 140., 200.},{0.0, 5., 10., 15., 20., 25., 30., 40., 50., 60., 70., 80., 100., 120., 140., 200.},{0.0, 5., 10., 15., 20., 25., 30., 40., 50., 60., 70., 80., 100., 120., 140., 200.}};
// int nBinPtJetsRec[nRadius] = {15,15,15};
// double ptBinsJetsGen[nRadius][30] = {{0.0, 5., 10., 15., 20., 25., 30., 40., 50., 60., 70., 80., 100., 120., 140., 200.},{0.0, 5., 10., 15., 20., 25., 30., 40., 50., 60., 70., 80., 100., 120., 140., 200.},{0.0, 5., 10., 15., 20., 25., 30., 40., 50., 60., 70., 80., 100., 120., 140., 200.}};
// int nBinPtJetsGen[nRadius] = {15,15,15};



// // tests
// // WORKS FINE WITH SVD PbPb
// // WORKS FINE WITH SVD
// // pT binning for jets - gen = rec - start at 10 // does work for svd even though 
// double ptBinsJetsRec[nRadius][30] = {{20., 30., 40., 50., 60., 70., 80., 90., 100., 110., 120., 130., 140.},{20., 30., 40., 50., 60., 70., 80., 90., 100., 110., 120., 130., 140.},{20., 30., 40., 50., 60., 70., 80., 90., 100., 110., 120., 130., 140.}};
// int nBinPtJetsRec[nRadius] = {12,12,12};
// double ptBinsJetsGen[nRadius][30] = {{5., 10., 20., 30., 40., 50., 60., 70., 80., 90., 100., 110., 120., 130., 140., 150., 160., 170., 180., 190., 200.},{5., 10., 20., 30., 40., 50., 60., 70., 80., 90., 100., 110., 120., 130., 140., 150., 160., 170., 180., 190., 200.},{5., 10., 20., 30., 40., 50., 60., 70., 80., 90., 100., 110., 120., 130., 140., 150., 160., 170., 180., 190., 200.}};
// int nBinPtJetsGen[nRadius] = {20,20,20};


// // PbPb
// // Hannah bossi identical fo run 2 comparison
// double ptBinsJetsRec[nRadius][30] = {{20., 25., 30., 35., 40., 45., 50., 55., 60., 65., 70., 75., 80., 85., 90., 95., 100., 110., 120., 140.},{20., 25., 30., 35., 40., 45., 50., 55., 60., 65., 70., 75., 80., 85., 90., 95., 100., 110., 120., 140.},{20., 25., 30., 35., 40., 45., 50., 55., 60., 65., 70., 75., 80., 85., 90., 95., 100., 110., 120., 140.}};
// int nBinPtJetsRec[nRadius] = {19,19,19};
// double ptBinsJetsGen[nRadius][30] = {{10., 20., 40., 60., 70., 85., 100., 120., 140., 200.},{10., 20., 40., 60., 70., 85., 100., 120., 140., 200.},{10., 20., 40., 60., 70., 85., 100., 120., 140., 200.}};
// int nBinPtJetsGen[nRadius] = {9,9,9};


// PbPb Aimeric default
double ptBinsJetsRec[nRadius][30] = {{10., 15., 20., 25., 30., 35., 40., 45., 50., 55., 60., 65., 70., 75., 80., 85., 90., 95., 100., 110., 120., 140.},{5., 10, 20., 25., 30., 35., 40., 45., 50., 55., 60., 65., 70., 75., 80., 85., 90., 95., 100., 110., 120., 140., 200.},{5., 10, 20., 25., 30., 35., 40., 45., 50., 55., 60., 65., 70., 75., 80., 85., 90., 95., 100., 110., 120., 140., 200.}};
// double ptBinsJetsRec[nRadius][30] = {{40., 45., 50., 55., 60., 65., 70., 75., 80., 85., 90., 95., 100., 110., 120., 140.},{5., 10, 20., 25., 30., 35., 40., 45., 50., 55., 60., 65., 70., 75., 80., 85., 90., 95., 100., 110., 120., 140., 200.},{5., 10, 20., 25., 30., 35., 40., 45., 50., 55., 60., 65., 70., 75., 80., 85., 90., 95., 100., 110., 120., 140., 200.}};
int nBinPtJetsRec[nRadius] = {21,22,22};
// int nBinPtJetsRec[nRadius] = {15,22,22};
double ptBinsJetsGen[nRadius][30] = {{0., 10., 20., 30., 40., 50., 60., 70., 80., 90., 100., 120., 140., 200.},{0., 10., 20., 30., 40., 50., 60., 70., 80., 90., 100., 120., 140., 200.},{0., 10., 20., 30., 40., 50., 60., 70., 80., 90., 100., 120., 140., 200.}};
int nBinPtJetsGen[nRadius] = {13,13,13};

// // PbPb Aimeric old
// double ptBinsJetsRec[nRadius][30] = {{20., 25., 30., 35., 40., 45., 50., 55., 60., 65., 70., 75., 80., 85., 90., 95., 100., 110., 120., 140., 200.},{20., 25., 30., 35., 40., 45., 50., 55., 60., 65., 70., 75., 80., 85., 90., 95., 100., 110., 120., 140., 200.},{20., 25., 30., 35., 40., 45., 50., 55., 60., 65., 70., 75., 80., 85., 90., 95., 100., 110., 120., 140., 200.}};
// int nBinPtJetsRec[nRadius] = {20,20,20};
// double ptBinsJetsGen[nRadius][30] = {{10., 20., 40., 60., 70., 85., 100., 120., 140., 200.},{10., 20., 40., 60., 70., 85., 100., 120., 140., 200.},{10., 20., 40., 60., 70., 85., 100., 120., 140., 200.}};
// int nBinPtJetsGen[nRadius] = {9,9,9};

// ///////////////Wenhui binning for Data unfolding GP PbPb MC/////////////////////////////
// double ptBinsJetsRec[nRadius][30] = {{-5, 0, 5, 10, 15, 20, 25, 30, 40, 50, 60, 70, 80, 90, 100, 120, 140, 200},
//                                      {-5, 0, 5, 10, 15, 20, 25, 30, 40, 50, 60, 70, 80, 90, 100, 120, 140, 200},
//                                      {-5, 0, 5, 10, 15, 20, 25, 30, 40, 50, 60, 70, 80, 90, 100, 120, 140, 200}};
// int nBinPtJetsRec[nRadius] = {17,17,17};
// double ptBinsJetsGen[nRadius][30] = {{-5, 0, 5, 10, 15, 20, 25, 30, 40, 50, 60, 70, 80, 90, 100, 120, 140, 200},
//                                     {-5, 0, 5, 10, 15, 20, 25, 30, 40, 50, 60, 70, 80, 90, 100, 120, 140, 200} ,
//                                      {-5, 0, 5, 10, 15, 20, 25, 30, 40, 50, 60, 70, 80, 90, 100, 120, 140, 200}};
// int nBinPtJetsGen[nRadius] = {17,17,17};


// // Joonsuk binning for pp
// double ptBinsJetsRec[nRadius][30] = {{5,  6,  7,  8,  9,  10, 12, 14,  16,  18, 20, 25, 30, 40, 50, 60, 70, 85, 100, 140, 200},{5,  6,  7,  8,  9,  10, 12, 14,  16,  18, 20, 25, 30, 40, 50, 60, 70, 85, 100, 140, 200},{5,  6,  7,  8,  9,  10, 12, 14,  16,  18, 20, 25, 30, 40, 50, 60, 70, 85, 100, 140, 200}};
// int nBinPtJetsRec[nRadius] = {20,20,20};
// double ptBinsJetsGen[nRadius][30] = {{5,  6,  7,  8,  9,  10, 12, 14,  16,  18, 20, 25, 30, 40, 50, 60, 70, 85, 100, 140, 200},{5,  6,  7,  8,  9,  10, 12, 14,  16,  18, 20, 25, 30, 40, 50, 60, 70, 85, 100, 140, 200},{5,  6,  7,  8,  9,  10, 12, 14,  16,  18, 20, 25, 30, 40, 50, 60, 70, 85, 100, 140, 200}};
// int nBinPtJetsGen[nRadius] = {20,20,20};


// Double_t ptbin[21] = {5,  6,  7,  8,  9,  10, 12, 14,  16,  18, 20, 25, 30, 40, 50, 60, 70, 85, 100, 140, 200};
// Double_t ptbinGen[26] = {0, 1, 2, 3, 4, 5,  6,  7,  8,  9,  10, 12, 14,  16,  18, 20, 25, 30, 40, 50, 60, 70, 85, 100, 140, 200};
// Int_t nptBins = sizeof(ptbin) / sizeof(ptbin[0]) - 1;
// Int_t nptBinsGen = sizeof(ptbinGen) / sizeof(ptbinGen[0]) - 1;


// // pT binning for jets - gen = rec - start at 5 but rec has a smaller window ; good to check stuff without worrying about a badly setup normalisation by pt bin width
// double ptBinsJetsRec[nRadius][30] = {{30., 40., 50., 60., 70., 80., 100., 120.},{30., 40., 50., 60., 70., 80., 100., 120.},{30., 40., 50., 60., 70., 80., 100., 120.}};
// int nBinPtJetsRec[nRadius] = {7,7,7};
// double ptBinsJetsGen[nRadius][30] = {{0., 10., 15., 20., 25., 30., 40., 50., 60., 70., 80., 100., 120., 140., 200.},{0., 10., 15., 20., 25., 30., 40., 50., 60., 70., 80., 100., 120., 140., 200.},{0., 10., 15., 20., 25., 30., 40., 50., 60., 70., 80., 100., 120., 140., 200.}};
// int nBinPtJetsGen[nRadius] = {14,14,14};


// // CURRENT VERSION TO USE OUTSIDE OF TESTS
// // Bayes version
// // pT binning for jets - hiroki tweak gen start at 10GeV , bin 140-200 subdivided; MARTA's version as well
// double ptBinsJetsRec[nRadius][30] = {{30., 35., 40., 45., 50., 55., 60., 65., 70., 75., 80., 85., 90., 95., 100., 110., 120.},{30., 35., 40., 45., 50., 55., 60., 65., 70., 75., 80., 85., 90., 95., 100., 110., 120.},{30., 35., 40., 45., 50., 55., 60., 65., 70., 75., 80., 85., 90., 95., 100., 110., 120.}};
// int nBinPtJetsRec[nRadius] = {16,16,16};
// double ptBinsJetsGen[nRadius][30] = {{0., 10., 20., 30., 40., 50., 60., 70., 80., 90., 100., 110., 120., 130., 140., 200.},{0., 10., 20., 30., 40., 50., 60., 70., 80., 90., 100., 110., 120., 130., 140., 200.},{0., 10., 20., 30., 40., 50., 60., 70., 80., 90., 100., 110., 120., 130., 140., 200.}};
// int nBinPtJetsGen[nRadius] = {15,15,15};
// // double ptBinsJetsGen[nRadius][30] = {{10., 20., 30., 40., 50., 60., 70., 90., 120., 200.},{10., 20., 30., 40., 50., 60., 70., 90., 120., 200.},{10., 20., 30., 40., 50., 60., 70., 90., 120., 200.}};
// // int nBinPtJetsGen[nRadius] = {9,9,9};


// // OLD
// // // WORKS FINE WITH SVD
// // // WORKS FINE WITH SVD?
// double ptBinsJetsRec[nRadius][30] = {{30., 35., 40., 45., 50., 55., 60., 65., 70., 75., 80., 85., 90., 95., 100., 110., 120.},{30., 35., 40., 45., 50., 55., 60., 65., 70., 75., 80., 85., 90., 95., 100., 110., 120.},{30., 35., 40., 45., 50., 55., 60., 65., 70., 75., 80., 85., 90., 95., 100., 110., 120.}};
// int nBinPtJetsRec[nRadius] = {16,16,16};
// double ptBinsJetsGen[nRadius][30] = {{0., 10., 20., 30., 40., 50., 60., 70., 90., 120., 200.},{0., 10., 20., 30., 40., 50., 60., 70., 90., 120., 200.},{0., 10., 20., 30., 40., 50., 60., 70., 90., 120., 200.}};
// int nBinPtJetsGen[nRadius] = {10,10,10};
// // double ptBinsJetsGen[nRadius][30] = {{10., 20., 30., 40., 50., 60., 70., 90., 120., 200.},{10., 20., 30., 40., 50., 60., 70., 90., 120., 200.},{10., 20., 30., 40., 50., 60., 70., 90., 120., 200.}};
// // int nBinPtJetsGen[nRadius] = {9,9,9};


// // test1 constant bin size
// double ptBinsJetsRec[nRadius][30] = {{30., 35., 40., 45., 50., 55., 60., 65., 70., 75., 80., 85., 90., 95., 100., 105., 110., 115., 120., 150.},{30., 35., 40., 45., 50., 55., 60., 65., 70., 75., 80., 85., 90., 95., 100., 105., 110., 115., 120., 150.},{30., 35., 40., 45., 50., 55., 60., 65., 70., 75., 80., 85., 90., 95., 100., 105., 110., 115., 120., 150.}};
// int nBinPtJetsRec[nRadius] = {19,19,19};
// double ptBinsJetsGen[nRadius][30] = {{0., 10., 20., 30., 40., 50., 60., 70., 80., 90., 100., 110., 120., 130., 140., 150., 160., 170., 180., 190., 200.},{0., 10., 20., 30., 40., 50., 60., 70., 80., 90., 100., 110., 120., 130., 140., 150., 160., 170., 180., 190., 200.},{0., 10., 20., 30., 40., 50., 60., 70., 80., 90., 100., 110., 120., 130., 140., 150., 160., 170., 180., 190., 200.}};
// int nBinPtJetsGen[nRadius] = {20,20,20};
// // double ptBinsJetsGen[nRadius][30] = {{10., 20., 30., 40., 50., 60., 70., 90., 120., 200.},{10., 20., 30., 40., 50., 60., 70., 90., 120., 200.},{10., 20., 30., 40., 50., 60., 70., 90., 120., 200.}};
// // int nBinPtJetsGen[nRadius] = {9,9,9};


// // pT binninb for jets tests
// int nBinPtJetsRec[nRadius] = {14,14,14};
// double ptBinsJetsRec[nRadius][20] = {{5., 10., 15., 20., 25., 30., 40., 50., 60., 70., 80., 100., 120., 140., 200.},{5., 10., 15., 20., 25., 30., 40., 50., 60., 70., 80., 100., 120., 140., 200.},{5., 10., 15., 20., 25., 30., 40., 50., 60., 70., 80., 100., 120., 140., 200.}};
// int nBinPtJetsGen[nRadius] = {26,26,26};
// double ptBinsJetsGen[nRadius][201] = {{  0.,  5.,
//                                         10., 15.,
//                                         20., 25.,
//                                         30., 35.,
//                                         40., 45.,
//                                         50., 55.,
//                                         60., 65.,
//                                         70., 75.,
//                                         80., 85.,
//                                         90., 95.,
//                                        100.,
//                                        110.,
//                                        120.,
//                                        130.,
//                                        140.,
//                                        150.,
//                                        200},               
//                                      {   0.,  5.,
//                                         10., 15.,
//                                         20., 25.,
//                                         30., 35.,
//                                         40., 45.,
//                                         50., 55.,
//                                         60., 65.,
//                                         70., 75.,
//                                         80., 85.,
//                                         90., 95.,
//                                        100.,
//                                        110.,
//                                        120.,
//                                        130.,
//                                        140.,
//                                        150.,
//                                        200},               
//                                      {   0.,  5.,
//                                         10., 15.,
//                                         20., 25.,
//                                         30., 35.,
//                                         40., 45.,
//                                         50., 55.,
//                                         60., 65.,
//                                         70., 75.,
//                                         80., 85.,
//                                         90., 95.,
//                                        100.,
//                                        110.,
//                                        120.,
//                                        130.,
//                                        140.,
//                                        150.,
//                                        200}};


// // pT binninb for jets tests
// int nBinPtJetsRec[nRadius] = {14,14,14};
// double ptBinsJetsRec[nRadius][20] = {{5., 10., 15., 20., 25., 30., 40., 50., 60., 70., 80., 100., 120., 140., 200.},{5., 10., 15., 20., 25., 30., 40., 50., 60., 70., 80., 100., 120., 140., 200.},{5., 10., 15., 20., 25., 30., 40., 50., 60., 70., 80., 100., 120., 140., 200.}};
// int nBinPtJetsGen[nRadius] = {50,50,50};
// double ptBinsJetsGen[nRadius][201] = {{  0.,  2.5,  5.,  7.5,
//                                         10., 12.5, 15., 17.5,
//                                         20., 22.5, 25., 27.5,
//                                         30., 32.5, 35., 37.5,
//                                         40., 42.5, 45., 47.5,
//                                         50., 52.5, 55., 57.5,
//                                         60., 62.5, 65., 67.5,
//                                         70., 72.5, 75., 77.5,
//                                         80., 82.5, 85., 87.5,
//                                         90., 92.5, 95., 97.5,
//                                        100.,
//                                        110.,
//                                        120.,
//                                        130.,
//                                        140.,
//                                        150.,
//                                        160.,
//                                        170.,
//                                        180.,
//                                        190.,
//                                        200},               
//                                      {   0.,  2.5,  5.,  7.5,
//                                         10., 12.5, 15., 17.5,
//                                         20., 22.5, 25., 27.5,
//                                         30., 32.5, 35., 37.5,
//                                         40., 42.5, 45., 47.5,
//                                         50., 52.5, 55., 57.5,
//                                         60., 62.5, 65., 67.5,
//                                         70., 72.5, 75., 77.5,
//                                         80., 82.5, 85., 87.5,
//                                         90., 92.5, 95., 97.5,
//                                        100.,
//                                        110.,
//                                        120.,
//                                        130.,
//                                        140.,
//                                        150.,
//                                        160.,
//                                        170.,
//                                        180.,
//                                        190.,
//                                        200},               
//                                      {   0.,  2.5,  5.,  7.5,
//                                         10., 12.5, 15., 17.5,
//                                         20., 22.5, 25., 27.5,
//                                         30., 32.5, 35., 37.5,
//                                         40., 42.5, 45., 47.5,
//                                         50., 52.5, 55., 57.5,
//                                         60., 62.5, 65., 67.5,
//                                         70., 72.5, 75., 77.5,
//                                         80., 82.5, 85., 87.5,
//                                         90., 92.5, 95., 97.5,
//                                        100.,
//                                        110.,
//                                        120.,
//                                        130.,
//                                        140.,
//                                        150.,
//                                        160.,
//                                        170.,
//                                        180.,
//                                        190.,
//                                        200}};
                            

// int nBinPtJetsFine[nRadius] = {40,40,40};
// // int nBinPtJetsFine[nRadius] = {195,195,195};
// // double ptBinsJetsFine[nRadius][201] = {{05., 06., 07., 08., 09.,
// double ptBinsJetsFine[nRadius][201] = {{   0., 05.,
//                                         10., 15.,
//                                         20., 25.,
//                                         30., 35.,
//                                         40., 45.,
//                                         50., 55.,
//                                         60., 65.,
//                                         70., 75.,
//                                         80., 85.,
//                                         90., 95.,
//                                        100., 105.,
//                                        110., 115.,
//                                        120., 125.,
//                                        130., 135.,
//                                        140., 145.,
//                                        150., 155.,
//                                        160., 165.,
//                                        170., 175.,
//                                        180., 185.,
//                                        190., 195.,
//                                        200.},
//                                      {   0., 05.,
//                                         10., 15.,
//                                         20., 25.,
//                                         30., 35.,
//                                         40., 45.,
//                                         50., 55.,
//                                         60., 65.,
//                                         70., 75.,
//                                         80., 85.,
//                                         90., 95.,
//                                        100., 105.,
//                                        110., 115.,
//                                        120., 125.,
//                                        130., 135.,
//                                        140., 145.,
//                                        150., 155.,
//                                        160., 165.,
//                                        170., 175.,
//                                        180., 185.,
//                                        190., 195.,
//                                        200.},
//                                      {   0., 05.,
//                                         10., 15.,
//                                         20., 25.,
//                                         30., 35.,
//                                         40., 45.,
//                                         50., 55.,
//                                         60., 65.,
//                                         70., 75.,
//                                         80., 85.,
//                                         90., 95.,
//                                        100., 105.,
//                                        110., 115.,
//                                        120., 125.,
//                                        130., 135.,
//                                        140., 145.,
//                                        150., 155.,
//                                        160., 165.,
//                                        170., 175.,
//                                        180., 185.,
//                                        190., 195.,
//                                        200.}}; // shift+option+left click hold lets one edit columns in vs code

// // fine binning standard, for pp and Pb-Pb factorised
// // int nBinPtJetsFine[nRadius] = {120,120,120};
// int nBinPtJetsFine[nRadius] = {115,115,115};
// // double ptBinsJetsFine[nRadius][201] = {{05., 06., 07., 08., 09.,
// double ptBinsJetsFine[nRadius][201] = {{ 0., 01., 02., 03., 04., 05., 06., 07., 08., 09.,
//                                         10., 11., 12., 13., 14., 15., 16., 17., 18., 19.,
//                                         20., 21., 22., 23., 24., 25., 26., 27., 28., 29.,
//                                         30., 31., 32., 33., 34., 35., 36., 37., 38., 39.,
//                                         40., 41., 42., 43., 44., 45., 46., 47., 48., 49.,
//                                         50., 51., 52., 53., 54., 55., 56., 57., 58., 59.,
//                                         60., 61., 62., 63., 64., 65., 66., 67., 68., 69.,
//                                         70., 71., 72., 73., 74., 75., 76., 77., 78., 79.,
//                                         80., 81., 82., 83., 84., 85., 86., 87., 88., 89.,
//                                         90., 91., 92., 93., 94., 95., 96., 97., 98., 99.,
//                                        100., 105.,
//                                        110., 115.,
//                                        120., 125.,
//                                        130., 135.,
//                                        140., 145.,
//                                        150., 155.,
//                                        160., 165.,
//                                        170., 175.,
//                                        180., 185.,
//                                        190., 195.,
//                                        200.},
//                                     //  {05., 06., 07., 08., 09.,
//                                      {   0., 01., 02., 03., 04., 05., 06., 07., 08., 09.,
//                                         10., 11., 12., 13., 14., 15., 16., 17., 18., 19.,
//                                         20., 21., 22., 23., 24., 25., 26., 27., 28., 29.,
//                                         30., 31., 32., 33., 34., 35., 36., 37., 38., 39.,
//                                         40., 41., 42., 43., 44., 45., 46., 47., 48., 49.,
//                                         50., 51., 52., 53., 54., 55., 56., 57., 58., 59.,
//                                         60., 61., 62., 63., 64., 65., 66., 67., 68., 69.,
//                                         70., 71., 72., 73., 74., 75., 76., 77., 78., 79.,
//                                         80., 81., 82., 83., 84., 85., 86., 87., 88., 89.,
//                                         90., 91., 92., 93., 94., 95., 96., 97., 98., 99.,
//                                        100., 105.,
//                                        110., 115.,
//                                        120., 125.,
//                                        130., 135.,
//                                        140., 145.,
//                                        150., 155.,
//                                        160., 165.,
//                                        170., 175.,
//                                        180., 185.,
//                                        190., 195.,
//                                        200.},
//                                     //  {05., 06., 07., 08., 09.,
//                                      {   0., 01., 02., 03., 04., 05., 06., 07., 08., 09.,
//                                         10., 11., 12., 13., 14., 15., 16., 17., 18., 19.,
//                                         20., 21., 22., 23., 24., 25., 26., 27., 28., 29.,
//                                         30., 31., 32., 33., 34., 35., 36., 37., 38., 39.,
//                                         40., 41., 42., 43., 44., 45., 46., 47., 48., 49.,
//                                         50., 51., 52., 53., 54., 55., 56., 57., 58., 59.,
//                                         60., 61., 62., 63., 64., 65., 66., 67., 68., 69.,
//                                         70., 71., 72., 73., 74., 75., 76., 77., 78., 79.,
//                                         80., 81., 82., 83., 84., 85., 86., 87., 88., 89.,
//                                         90., 91., 92., 93., 94., 95., 96., 97., 98., 99.,
//                                        100., 105.,
//                                        110., 115.,
//                                        120., 125.,
//                                        130., 135.,
//                                        140., 145.,
//                                        150., 155.,
//                                        160., 165.,
//                                        170., 175.,
//                                        180., 185.,
//                                        190., 195.,
//                                        200.},}; // shift+option+left click hold lets one edit columns in vs code



// // fine binning standard, for Pb-Pb non-factorised
// // int nBinPtJetsFine[nRadius] = {120,120,120};
// int nBinPtJetsFine[nRadius] = {240,240,240};
// // double ptBinsJetsFine[nRadius][201] = {{05., 06., 07., 08., 09.,
// double ptBinsJetsFine[nRadius][402] = {
//                                     { -200., -195.,
//                                       -190., -185.,
//                                       -180., -175.,
//                                       -170., -165.,
//                                       -160., -155.,
//                                       -150., -145.,
//                                       -140., -135.,
//                                       -130., -125.,
//                                       -120., -115.,
//                                       -110., -105.,
//                                       -100.,-99.,-98.,-97.,-96.,-95.,-94.,-93.,-92.,-91.,
//                                        -90.,-89.,-88.,-87.,-86.,-85.,-84.,-83.,-82.,-81.,
//                                        -80.,-79.,-78.,-77.,-76.,-75.,-74.,-73.,-72.,-71.,
//                                        -70.,-69.,-68.,-67.,-66.,-65.,-64.,-63.,-62.,-61.,
//                                        -60.,-59.,-58.,-57.,-56.,-55.,-54.,-53.,-52.,-51.,
//                                        -50.,-49.,-48.,-47.,-46.,-45.,-44.,-43.,-42.,-41.,
//                                        -40.,-39.,-38.,-37.,-36.,-35.,-34.,-33.,-32.,-31.,
//                                        -30.,-29.,-28.,-27.,-26.,-25.,-24.,-23.,-22.,-21.,
//                                        -20.,-19.,-18.,-17.,-16.,-15.,-14.,-13.,-12.,-11.,
//                                        -10.,-09.,-08.,-07.,-06.,-05.,-04.,-03.,-02.,-01.,
//                                          0., 01., 02., 03., 04., 05., 06., 07., 08., 09.,
//                                         10., 11., 12., 13., 14., 15., 16., 17., 18., 19.,
//                                         20., 21., 22., 23., 24., 25., 26., 27., 28., 29.,
//                                         30., 31., 32., 33., 34., 35., 36., 37., 38., 39.,
//                                         40., 41., 42., 43., 44., 45., 46., 47., 48., 49.,
//                                         50., 51., 52., 53., 54., 55., 56., 57., 58., 59.,
//                                         60., 61., 62., 63., 64., 65., 66., 67., 68., 69.,
//                                         70., 71., 72., 73., 74., 75., 76., 77., 78., 79.,
//                                         80., 81., 82., 83., 84., 85., 86., 87., 88., 89.,
//                                         90., 91., 92., 93., 94., 95., 96., 97., 98., 99.,
//                                        100., 105.,
//                                        110., 115.,
//                                        120., 125.,
//                                        130., 135.,
//                                        140., 145.,
//                                        150., 155.,
//                                        160., 165.,
//                                        170., 175.,
//                                        180., 185.,
//                                        190., 195.,
//                                        200.},
//                                     { -200., -195.,
//                                       -190., -185.,
//                                       -180., -175.,
//                                       -170., -165.,
//                                       -160., -155.,
//                                       -150., -145.,
//                                       -140., -135.,
//                                       -130., -125.,
//                                       -120., -115.,
//                                       -110., -105.,
//                                       -100.,-99.,-98.,-97.,-96.,-95.,-94.,-93.,-92.,-91.,
//                                        -90.,-89.,-88.,-87.,-86.,-85.,-84.,-83.,-82.,-81.,
//                                        -80.,-79.,-78.,-77.,-76.,-75.,-74.,-73.,-72.,-71.,
//                                        -70.,-69.,-68.,-67.,-66.,-65.,-64.,-63.,-62.,-61.,
//                                        -60.,-59.,-58.,-57.,-56.,-55.,-54.,-53.,-52.,-51.,
//                                        -50.,-49.,-48.,-47.,-46.,-45.,-44.,-43.,-42.,-41.,
//                                        -40.,-39.,-38.,-37.,-36.,-35.,-34.,-33.,-32.,-31.,
//                                        -30.,-29.,-28.,-27.,-26.,-25.,-24.,-23.,-22.,-21.,
//                                        -20.,-19.,-18.,-17.,-16.,-15.,-14.,-13.,-12.,-11.,
//                                        -10.,-09.,-08.,-07.,-06.,-05.,-04.,-03.,-02.,-01.,
//                                          0., 01., 02., 03., 04., 05., 06., 07., 08., 09.,
//                                         10., 11., 12., 13., 14., 15., 16., 17., 18., 19.,
//                                         20., 21., 22., 23., 24., 25., 26., 27., 28., 29.,
//                                         30., 31., 32., 33., 34., 35., 36., 37., 38., 39.,
//                                         40., 41., 42., 43., 44., 45., 46., 47., 48., 49.,
//                                         50., 51., 52., 53., 54., 55., 56., 57., 58., 59.,
//                                         60., 61., 62., 63., 64., 65., 66., 67., 68., 69.,
//                                         70., 71., 72., 73., 74., 75., 76., 77., 78., 79.,
//                                         80., 81., 82., 83., 84., 85., 86., 87., 88., 89.,
//                                         90., 91., 92., 93., 94., 95., 96., 97., 98., 99.,
//                                        100., 105.,
//                                        110., 115.,
//                                        120., 125.,
//                                        130., 135.,
//                                        140., 145.,
//                                        150., 155.,
//                                        160., 165.,
//                                        170., 175.,
//                                        180., 185.,
//                                        190., 195.,
//                                        200.},
//                                     { -200., -195.,
//                                       -190., -185.,
//                                       -180., -175.,
//                                       -170., -165.,
//                                       -160., -155.,
//                                       -150., -145.,
//                                       -140., -135.,
//                                       -130., -125.,
//                                       -120., -115.,
//                                       -110., -105.,
//                                       -100.,-99.,-98.,-97.,-96.,-95.,-94.,-93.,-92.,-91.,
//                                        -90.,-89.,-88.,-87.,-86.,-85.,-84.,-83.,-82.,-81.,
//                                        -80.,-79.,-78.,-77.,-76.,-75.,-74.,-73.,-72.,-71.,
//                                        -70.,-69.,-68.,-67.,-66.,-65.,-64.,-63.,-62.,-61.,
//                                        -60.,-59.,-58.,-57.,-56.,-55.,-54.,-53.,-52.,-51.,
//                                        -50.,-49.,-48.,-47.,-46.,-45.,-44.,-43.,-42.,-41.,
//                                        -40.,-39.,-38.,-37.,-36.,-35.,-34.,-33.,-32.,-31.,
//                                        -30.,-29.,-28.,-27.,-26.,-25.,-24.,-23.,-22.,-21.,
//                                        -20.,-19.,-18.,-17.,-16.,-15.,-14.,-13.,-12.,-11.,
//                                        -10.,-09.,-08.,-07.,-06.,-05.,-04.,-03.,-02.,-01.,
//                                          0., 01., 02., 03., 04., 05., 06., 07., 08., 09.,
//                                         10., 11., 12., 13., 14., 15., 16., 17., 18., 19.,
//                                         20., 21., 22., 23., 24., 25., 26., 27., 28., 29.,
//                                         30., 31., 32., 33., 34., 35., 36., 37., 38., 39.,
//                                         40., 41., 42., 43., 44., 45., 46., 47., 48., 49.,
//                                         50., 51., 52., 53., 54., 55., 56., 57., 58., 59.,
//                                         60., 61., 62., 63., 64., 65., 66., 67., 68., 69.,
//                                         70., 71., 72., 73., 74., 75., 76., 77., 78., 79.,
//                                         80., 81., 82., 83., 84., 85., 86., 87., 88., 89.,
//                                         90., 91., 92., 93., 94., 95., 96., 97., 98., 99.,
//                                        100., 105.,
//                                        110., 115.,
//                                        120., 125.,
//                                        130., 135.,
//                                        140., 145.,
//                                        150., 155.,
//                                        160., 165.,
//                                        170., 175.,
//                                        180., 185.,
//                                        190., 195.,
//                                        200.}}; // shift+option+left click hold lets one edit columns in vs code



// int nBinPtJetsFine[nRadius] = {200,200,200};
int nBinPtJetsFine[nRadius] = {195,195,195};
double ptBinsJetsFine[nRadius][201] = {{05., 06., 07., 08., 09.,
// double ptBinsJetsFine[nRadius][201] = {{ 0., 01., 02., 03., 04., 05., 06., 07., 08., 09.,
                                        10., 11., 12., 13., 14., 15., 16., 17., 18., 19.,
                                        20., 21., 22., 23., 24., 25., 26., 27., 28., 29.,
                                        30., 31., 32., 33., 34., 35., 36., 37., 38., 39.,
                                        40., 41., 42., 43., 44., 45., 46., 47., 48., 49.,
                                        50., 51., 52., 53., 54., 55., 56., 57., 58., 59.,
                                        60., 61., 62., 63., 64., 65., 66., 67., 68., 69.,
                                        70., 71., 72., 73., 74., 75., 76., 77., 78., 79.,
                                        80., 81., 82., 83., 84., 85., 86., 87., 88., 89.,
                                        90., 91., 92., 93., 94., 95., 96., 97., 98., 99.,
                                       100.,101.,102.,103.,104.,105.,106.,107.,108.,109.,
                                       110.,111.,112.,113.,114.,115.,116.,117.,118.,119.,
                                       120.,121.,122.,123.,124.,125.,126.,127.,128.,129.,
                                       130.,131.,132.,133.,134.,135.,136.,137.,138.,139.,
                                       140.,141.,142.,143.,144.,145.,146.,147.,148.,149.,
                                       150.,151.,152.,153.,154.,155.,156.,157.,158.,159.,
                                       160.,161.,162.,163.,164.,165.,166.,167.,168.,169.,
                                       170.,171.,172.,173.,174.,175.,176.,177.,178.,179.,
                                       180.,181.,182.,183.,184.,185.,186.,187.,188.,189.,
                                       190.,191.,192.,193.,194.,195.,196.,197.,198.,199.,
                                       200},
                                     {05., 06., 07., 08., 09.,
                                    //  {   0., 01., 02., 03., 04., 05., 06., 07., 08., 09.,
                                        10., 11., 12., 13., 14., 15., 16., 17., 18., 19.,
                                        20., 21., 22., 23., 24., 25., 26., 27., 28., 29.,
                                        30., 31., 32., 33., 34., 35., 36., 37., 38., 39.,
                                        40., 41., 42., 43., 44., 45., 46., 47., 48., 49.,
                                        50., 51., 52., 53., 54., 55., 56., 57., 58., 59.,
                                        60., 61., 62., 63., 64., 65., 66., 67., 68., 69.,
                                        70., 71., 72., 73., 74., 75., 76., 77., 78., 79.,
                                        80., 81., 82., 83., 84., 85., 86., 87., 88., 89.,
                                        90., 91., 92., 93., 94., 95., 96., 97., 98., 99.,
                                       100.,101.,102.,103.,104.,105.,106.,107.,108.,109.,
                                       110.,111.,112.,113.,114.,115.,116.,117.,118.,119.,
                                       120.,121.,122.,123.,124.,125.,126.,127.,128.,129.,
                                       130.,131.,132.,133.,134.,135.,136.,137.,138.,139.,
                                       140.,141.,142.,143.,144.,145.,146.,147.,148.,149.,
                                       150.,151.,152.,153.,154.,155.,156.,157.,158.,159.,
                                       160.,161.,162.,163.,164.,165.,166.,167.,168.,169.,
                                       170.,171.,172.,173.,174.,175.,176.,177.,178.,179.,
                                       180.,181.,182.,183.,184.,185.,186.,187.,188.,189.,
                                       190.,191.,192.,193.,194.,195.,196.,197.,198.,199.,
                                       200},
                                     {05., 06., 07., 08., 09.,
                                    //  {   0., 01., 02., 03., 04., 05., 06., 07., 08., 09.,
                                        10., 11., 12., 13., 14., 15., 16., 17., 18., 19.,
                                        20., 21., 22., 23., 24., 25., 26., 27., 28., 29.,
                                        30., 31., 32., 33., 34., 35., 36., 37., 38., 39.,
                                        40., 41., 42., 43., 44., 45., 46., 47., 48., 49.,
                                        50., 51., 52., 53., 54., 55., 56., 57., 58., 59.,
                                        60., 61., 62., 63., 64., 65., 66., 67., 68., 69.,
                                        70., 71., 72., 73., 74., 75., 76., 77., 78., 79.,
                                        80., 81., 82., 83., 84., 85., 86., 87., 88., 89.,
                                        90., 91., 92., 93., 94., 95., 96., 97., 98., 99.,
                                       100.,101.,102.,103.,104.,105.,106.,107.,108.,109.,
                                       110.,111.,112.,113.,114.,115.,116.,117.,118.,119.,
                                       120.,121.,122.,123.,124.,125.,126.,127.,128.,129.,
                                       130.,131.,132.,133.,134.,135.,136.,137.,138.,139.,
                                       140.,141.,142.,143.,144.,145.,146.,147.,148.,149.,
                                       150.,151.,152.,153.,154.,155.,156.,157.,158.,159.,
                                       160.,161.,162.,163.,164.,165.,166.,167.,168.,169.,
                                       170.,171.,172.,173.,174.,175.,176.,177.,178.,179.,
                                       180.,181.,182.,183.,184.,185.,186.,187.,188.,189.,
                                       190.,191.,192.,193.,194.,195.,196.,197.,198.,199.,
                                       200}}; // shift+option+left click hold lets one edit columns in vs code


#endif