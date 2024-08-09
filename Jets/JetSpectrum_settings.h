// Analysis settings
const int nJetType = 3;
const TString jetType[nJetType] = {"charged", "neutral", "full"};
const int nJetLevel = 3;
const TString jetLevel[nJetLevel] = {"data", "mcd", "mcp"};
const int nRadius = 3;
const TString RadiusLegend[nRadius] = {"R = 0.2", "R = 0.4", "R = 0.6"};
float arrayRadius[nRadius] = {0.2, 0.4, 0.6};
const double etaAnalysisRange[2] = {-0.5, 0.5};
// Choice of jet type (charged, neutral, full) and level (data, detector level, particle level)
const int iJetType = 0;
const int iJetLevel = 0;
// centrality binning for collisions
// const int nCentralityBins = 2;
// const float arrayCentralityIntervals[nCentralityBins][2] = {{0, 10}, {50, 90}};
const int nCentralityBins = 1;
const float arrayCentralityIntervals[nCentralityBins][2] = {{50, 90}};






////////////////////////////////////////////////
////////// Unfolding settings - start //////////
////////////////////////////////////////////////


const bool doEvtNorm = 1;                               //  as of now, only changes the scale in the result, nothing more, so I should be fine using it
const bool doWidthScaling = 1;                          //  doesn't seem to have any effect, so I can probably use it: doesn't change the ratios (at least measured/unfolded and mcp/unfolded, haven't checked folded/unfolded)

char mergingPrior[] = "noMergingPrior";     // prior options: mcpPriorMerging, mcdPriorMerging, measuredPriorMerging, noMergingPrior, testAliPhysics
char unfoldingPrior[] = "noPrior";     // prior options: mcpPriorUnfolding, mcdPriorUnfolding, measuredPriorUnfolding, noPrior, testAliPhysics
char unfoldingMethod[] = "Svd"; // unfolding method options: Bayes, Svd
char normMethod[] = "evtNorm"; // evtNorm, noNorm
char optionsAnalysis[100] = "";

const bool isPbPb = false; // if false -> pp
const bool ppMcIsWeighted = false; // use if the MC has been weighted to have more high pt jets?
int applyEfficiencies = 3; // for test purposes: 0: no efficiency correction, 1: kine only, 2: jet finding efficiency only, 3: both active; only applied if useInitialResponseMethod is true
bool applyFakes = true; // only applied if useInitialResponseMethod is true
const bool useFineBinningTest = false; //looks like this gives the same flat distrib as when using coarse binning: so rebinning isnt the issue; need to change finBinning back to start at 0 when I dont use this
const bool scaleRespByWidth = false;

const bool useYSliceNorm = false; // doesn't change result for svd, doesn't work well for bayes though; SHOULD BE TRUE IF USING PRIOR for sure

// all three below should probably be true;
// but then it breaks svd convergence! find out why;
const bool normDetRespByNEvts = false; //that's what breaks svd; https://arxiv.org/pdf/hep-ph/9509307 seems to say one should use the number of events matrix (see last paragraph of conclusion) instead of a probability matrix
const bool normGenAndMeasByNEvts = false;
const bool normunfoldedByNEvts = false;

const bool useInitialResponseMethod = false; // discrepancy false true here seems to be that I do not model fakes in my initial method
const bool normaliseRespYSliceForRefold = true; // ??????? THAT IS APPARENTLY REQUIRED TO REFOLD MANUALLY! even though the initial resp matrix used for the unfolding isn't normalised like this

bool controlMC = true; // not yet implemented: weighted control MC, and control for PbPb

const bool drawIntermediateResponseMatrices = false;

// // save for PbPb that works, ish?
// const bool doEvtNorm = 1;                               //  as of now, only changes the scale in the result, nothing more, so I should be fine using it
// const bool doWidthScaling = 1;                          //  doesn't seem to have any effect, so I can probably use it: doesn't change the ratios (at least measured/unfolded and mcp/unfolded, haven't checked folded/unfolded)

// char mergingPrior[] = "mcpPriorMerging";     // prior options: mcpPriorMerging, mcdPriorMerging, measuredPriorMerging, noMergingPrior, testAliPhysics
// char unfoldingPrior[] = "noPrior";     // prior options: mcpPriorUnfolding, mcdPriorUnfolding, measuredPriorUnfolding, noPrior, testAliPhysics
// char unfoldingMethod[] = "Bayes"; // unfolding method options: Bayes, Svd, BinByBin
// char normMethod[] = "evtNorm"; // evtNorm, noNorm
// char optionsAnalysis[100] = "";

// const bool isPbPb = true; // if false -> pp
// const bool ppMcIsWeighted = true; // use if the MC has been weighted to have more high pt jets?
// int applyEfficiencies = 3; // for test purposes: 0: no efficiency correction, 1: kine only, 2: jet finding efficiency only, 3: both active
// const bool useFineBinningTest = false; //looks like this gives the same flat distrib as when using coarse binning: so rebinning isnt the issue; need to change finBinning back to start at 0 when I dont use this
// const bool useYSliceNorm = false; // helps measuredPriorUnfolding with useInitialResponseMethod = true; but mcpPriorUnfolding isn't great
// const bool useInitialResponseMethod = true; // if false, applyEfficiencies must be 0, if true it must be 3
// bool normaliseRespYSliceForRefold = true; // ??????? THAT IS APPARENTLY REQUIRED TO REFOLD MANUALLY! even though the initial resp matrix used for the unfolding isn't normalised like this

// // less used, test purposes
// const bool scaleRespByWidth = false; // test purposes
// const bool normDetRespByNEvts = false; // test purposes
// const bool normGenAndMeasByNEvts = false; // test purposes
// // end of save

// // Save of settings working for pp with // Joonsuk binning
// const bool doEvtNorm = 1;                               //  as of now, only changes the scale in the result, nothing more, so I should be fine using it
// const bool doWidthScaling = 1;                          //  doesn't seem to have any effect, so I can probably use it: doesn't change the ratios (at least measured/unfolded and mcp/unfolded, haven't checked folded/unfolded)

// char mergingPrior[] = "noMergingPrior";     // prior options: mcpPriorMerging, mcdPriorMerging, measuredPriorMerging, noMergingPrior, testAliPhysics
// char unfoldingPrior[] = "noPrior";     // prior options: mcpPriorUnfolding, mcdPriorUnfolding, measuredPriorUnfolding, noPrior, testAliPhysics
// char unfoldingMethod[] = "Bayes"; // unfolding method options: Bayes, Svd, BinByBin
// char normMethod[] = "evtNorm"; // evtNorm, noNorm
// char optionsAnalysis[100] = "";

// const bool isPbPb = false; // if false -> pp
// const bool ppMcIsWeighted = false; // use if the MC has been weighted to have more high pt jets?
// int applyEfficiencies = 3; // for test purposes: 0: no efficiency correction, 1: kine only, 2: jet finding efficiency only, 3: both active
// const bool useFineBinningTest = false; //looks like this gives the same flat distrib as when using coarse binning: so rebinning isnt the issue; need to change finBinning back to start at 0 when I dont use this
// const bool useYSliceNorm = false;
// const bool scaleRespByWidth = false;
// const bool normDetRespByNEvts = false;

// const bool useInitialResponseMethod = true;
// const bool normaliseRespYSliceForRefold = true; // ??????? THAT IS APPARENTLY REQUIRED TO REFOLD MANUALLY! even though the initial resp matrix used for the unfolding isn't normalised like this
// // end of save

////////////////////////////////////////////////
////////// Unfolding settings - end ////////////
////////////////////////////////////////////////



////////////////////////////////////////////////
//////////////// pt binning options ////////////
////////////////////////////////////////////////


// // pT binning for jets - gen = rec
// double ptBinsJetsRec[nRadius][30] = {{0.0, 5., 10., 15., 20., 25., 30., 40., 50., 60., 70., 80., 100., 120., 140., 200.},{0.0, 5., 10., 15., 20., 25., 30., 40., 50., 60., 70., 80., 100., 120., 140., 200.},{0.0, 5., 10., 15., 20., 25., 30., 40., 50., 60., 70., 80., 100., 120., 140., 200.}};
// int nBinPtJetsRec[nRadius] = {15,15,15};
// double ptBinsJetsGen[nRadius][30] = {{0.0, 5., 10., 15., 20., 25., 30., 40., 50., 60., 70., 80., 100., 120., 140., 200.},{0.0, 5., 10., 15., 20., 25., 30., 40., 50., 60., 70., 80., 100., 120., 140., 200.},{0.0, 5., 10., 15., 20., 25., 30., 40., 50., 60., 70., 80., 100., 120., 140., 200.}};
// int nBinPtJetsGen[nRadius] = {15,15,15};



// // tests
// // WORKS FINE WITH SVD
// // WORKS FINE WITH SVD
// // pT binning for jets - gen = rec - start at 10 // does work for svd even though 
// double ptBinsJetsRec[nRadius][30] = {{20., 30., 40., 50., 60., 70., 80., 90., 100., 110., 120., 130., 140.},{20., 30., 40., 50., 60., 70., 80., 90., 100., 110., 120., 130., 140.},{20., 30., 40., 50., 60., 70., 80., 90., 100., 110., 120., 130., 140.}};
// int nBinPtJetsRec[nRadius] = {12,12,12};
// double ptBinsJetsGen[nRadius][30] = {{5., 10., 20., 30., 40., 50., 60., 70., 80., 90., 100., 110., 120., 130., 140., 150., 160., 170., 180., 190., 200.},{5., 10., 20., 30., 40., 50., 60., 70., 80., 90., 100., 110., 120., 130., 140., 150., 160., 170., 180., 190., 200.},{5., 10., 20., 30., 40., 50., 60., 70., 80., 90., 100., 110., 120., 130., 140., 150., 160., 170., 180., 190., 200.}};
// int nBinPtJetsGen[nRadius] = {20,20,20};

// Joonsuk binning
double ptBinsJetsRec[nRadius][30] = {{5,  6,  7,  8,  9,  10, 12, 14,  16,  18, 20, 25, 30, 40, 50, 60, 70, 85, 100, 140, 200},{5,  6,  7,  8,  9,  10, 12, 14,  16,  18, 20, 25, 30, 40, 50, 60, 70, 85, 100, 140, 200},{5,  6,  7,  8,  9,  10, 12, 14,  16,  18, 20, 25, 30, 40, 50, 60, 70, 85, 100, 140, 200}};
int nBinPtJetsRec[nRadius] = {20,20,20};
double ptBinsJetsGen[nRadius][30] = {{0, 1, 2, 3, 4, 5,  6,  7,  8,  9,  10, 12, 14,  16,  18, 20, 25, 30, 40, 50, 60, 70, 85, 100, 140, 200},{0, 1, 2, 3, 4, 5,  6,  7,  8,  9,  10, 12, 14,  16,  18, 20, 25, 30, 40, 50, 60, 70, 85, 100, 140, 200},{0, 1, 2, 3, 4, 5,  6,  7,  8,  9,  10, 12, 14,  16,  18, 20, 25, 30, 40, 50, 60, 70, 85, 100, 140, 200}};
int nBinPtJetsGen[nRadius] = {25,25,25};


// Double_t ptbin[21] = {5,  6,  7,  8,  9,  10, 12, 14,  16,  18, 20, 25, 30, 40, 50, 60, 70, 85, 100, 140, 200};
// Double_t ptbinGen[26] = {0, 1, 2, 3, 4, 5,  6,  7,  8,  9,  10, 12, 14,  16,  18, 20, 25, 30, 40, 50, 60, 70, 85, 100, 140, 200};
// Int_t nptBins = sizeof(ptbin) / sizeof(ptbin[0]) - 1;
// Int_t nptBinsGen = sizeof(ptbinGen) / sizeof(ptbinGen[0]) - 1;


// // pT binning for jets - gen = rec - start at 5 but rec has a smaller window ; good to check stuff without worrying about a badly setup normalisation by pt bin width
// double ptBinsJetsRec[nRadius][30] = {{30., 40., 50., 60., 70., 80., 100., 120.},{30., 40., 50., 60., 70., 80., 100., 120.},{30., 40., 50., 60., 70., 80., 100., 120.}};
// int nBinPtJetsRec[nRadius] = {7,7,7};
// double ptBinsJetsGen[nRadius][30] = {{0., 10., 15., 20., 25., 30., 40., 50., 60., 70., 80., 100., 120., 140., 200.},{0., 10., 15., 20., 25., 30., 40., 50., 60., 70., 80., 100., 120., 140., 200.},{0., 10., 15., 20., 25., 30., 40., 50., 60., 70., 80., 100., 120., 140., 200.}};
// int nBinPtJetsGen[nRadius] = {14,14,14};


// // WORKS FINE WITH SVD
// // WORKS FINE WITH SVD
// // CURRENT VERSION TO USE OUTSIDE OF TESTS
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
                            

int nBinPtJetsFine[nRadius] = {200,200,200};
// int nBinPtJetsFine[nRadius] = {195,195,195};
// double ptBinsJetsFine[nRadius][201] = {{05., 06., 07., 08., 09.,
double ptBinsJetsFine[nRadius][201] = {{ 0., 01., 02., 03., 04., 05., 06., 07., 08., 09.,
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
                                    //  {05., 06., 07., 08., 09.,
                                     {   0., 01., 02., 03., 04., 05., 06., 07., 08., 09.,
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
                                    //  {05., 06., 07., 08., 09.,
                                     {   0., 01., 02., 03., 04., 05., 06., 07., 08., 09.,
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









////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////       file access choice       ////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////

// //////// -------- Pb-Pb -------- ////////
// TString* texCollisionDataInfo = new TString("Pb-Pb #sqrt{#it{s}} = 5.36 TeV"); 
// const TString* texDatasetsComparisonType = new TString("");
// const TString* texDatasetsComparisonCommonDenominator = new TString("LHC23zzh");
// const int nDatasets = 1;
// const TString Datasets[nDatasets] = {"LHC23zzh"};
// const TString DatasetsNames[nDatasets] = {""};
// TFile* file_O2Analysis_list[nDatasets] = {new TFile("Datasets/"+Datasets[0]+"/AnalysisResults.root")
//                                       };
// // TFile* file_O2Analysis_ppSimDetectorEffect[nDatasets] = new TFile("Datasets/ppSim_LHC23d4/AnalysisResults.root");
// TFile* file_O2Analysis_ppSimDetectorEffect = new TFile("Datasets/ppSim_LHC23d4_weighted/AnalysisResults.root");

// const TString trainId = "_id12832";
// // const TString trainId = "_id12832";
// // const TString analysisWorkflowData = "jet-finder-charged-qa_central_5090_lead5"+trainId;
// const TString analysisWorkflowData = "jet-finder-charged-qa_central_5090_lead5"+trainId;

// const TString analysisWorkflowMC = "jet-finder-charged-qa";





//////// -------- pp -------- ////////
TString* texCollisionDataInfo = new TString("pp #sqrt{#it{s}} = 13.6 TeV"); 
const TString* texDatasetsComparisonType = new TString("");
const TString* texDatasetsComparisonCommonDenominator = new TString("LHC22o pass6");
const int nDatasets = 1;
const TString Datasets[nDatasets] = {"LHC22o_pass6_train238827"};
const TString DatasetsNames[nDatasets] = {"LHC22o_pass6"};
TFile* file_O2Analysis_list[nDatasets] = {new TFile("Datasets/"+Datasets[0]+"/AnalysisResults.root")
                                      };
TFile* file_O2Analysis_ppSimDetectorEffect = new TFile("Datasets/LHC24b1b_sel8MC_train239181/AnalysisResults.root");

// TFile* file_O2Analysis_ppSimDetectorEffect_unfoldingControl = {new TFile("Datasets/LHC24f3_sel8MC_train240962/AnalysisResults.root")};
TFile* file_O2Analysis_ppSimDetectorEffect_unfoldingControl = {new TFile("Datasets/LHC24b1b_sel8MC_train239181/OneRun/AnalysisResults.root")};

// TFile* file_O2Analysis_ppSimDetectorEffect = new TFile("Datasets/LHC24b1b_sel8Full_train239409/AnalysisResults.root");
// TFile* file_O2Analysis_ppSimDetectorEffect = new TFile("Datasets/ppSim_LHC23d4_weighted/AnalysisResults.root");

const TString trainId = "";
// const TString trainId = "_id12832";
// const TString analysisWorkflowData = "jet-finder-charged-qa_central_5090_lead5"+trainId;
const TString analysisWorkflowData = "jet-finder-charged-qa"+trainId;

const TString analysisWorkflowMC = "jet-finder-charged-qa";
// const TString analysisWorkflowMC = "jet-finder-charged-qa_global_CollMatch";






















// Sum(Mik * Gen_k, k) = Rec_i 

// I want Sum(Rec_i, i) = Sum(Gen_i, i)

// Sum(Sum(Mik * Gen_k, k), i) = Sum(Rec_i, i)
// Sum(Sum(Mik, i) * Gen_k, k) = Sum(Rec_i, i)
// so I want Sum(Mik, i) = 1

