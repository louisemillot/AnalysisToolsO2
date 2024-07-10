

// Analysis settings
const int nJetType = 3;
const TString jetType[nJetType] = {"charged", "neutral", "full"};
const int nJetLevel = 3;
const TString jetLevel[nJetLevel] = {"data", "mcd", "mcp"};
const int nRadius = 3;
const TString RadiusLegend[nRadius] = {"R = 0.2", "R = 0.4", "R = 0.6"};
float arrayRadius[nRadius] = {0.2, 0.4, 0.6};
// const int nRadius = 1;
// const TString RadiusLegend[nRadius] = {"R = 0.4"};
// float arrayRadius[nRadius] = {0.4};
// const int nRadius = 9;
// const TString RadiusLegend[nRadius] = {"R = 0.2", "R = 0.25", "R = 0.3", "R = 0.35", "R = 0.4", "R = 0.45", "R = 0.5", "R = 0.55", "R = 0.6"};
// float arrayRadius[nRadius] = {0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6};

// Choice of jet type (charged, neutral, full) and level (data, detector level, particle level)
const int iJetType = 0;
const int iJetLevel = 0;


const int ptRebinValue = 1;


// centrality binning for collisions
// const int nCentralityBins = 2;
// const float arrayCentralityIntervals[nCentralityBins][2] = {{0, 10}, {50, 90}};
const int nCentralityBins = 1;
const float arrayCentralityIntervals[nCentralityBins][2] = {{50, 90}}; // for now uses this because pp sim I got has all its h3_jet_r_jet_pt_centrality entries below 0
const bool doNormalisation = 0;

// // pT binning for jets - gen = rec
// double ptBinsJetsRec[nRadius][30] = {{0.0, 5., 10., 15., 20., 25., 30., 40., 50., 60., 70., 80., 100., 120., 140., 200.},{0.0, 5., 10., 15., 20., 25., 30., 40., 50., 60., 70., 80., 100., 120., 140., 200.},{0.0, 5., 10., 15., 20., 25., 30., 40., 50., 60., 70., 80., 100., 120., 140., 200.}};
// int nBinPtJetsRec[nRadius] = {15,15,15};
// double ptBinsJetsGen[nRadius][30] = {{0.0, 5., 10., 15., 20., 25., 30., 40., 50., 60., 70., 80., 100., 120., 140., 200.},{0.0, 5., 10., 15., 20., 25., 30., 40., 50., 60., 70., 80., 100., 120., 140., 200.},{0.0, 5., 10., 15., 20., 25., 30., 40., 50., 60., 70., 80., 100., 120., 140., 200.}};
// int nBinPtJetsGen[nRadius] = {15,15,15};


// // pT binning for jets - hiroki layout
// double ptBinsJetsRec[nRadius][30] = {{30., 35., 40., 45., 50., 55., 60., 65., 70., 75., 80., 85., 90., 100., 120.},{30., 35., 40., 45., 50., 55., 60., 65., 70., 75., 80., 85., 90., 100., 120.},{30., 35., 40., 45., 50., 55., 60., 65., 70., 75., 80., 85., 90., 100., 120.}};
// int nBinPtJetsRec[nRadius] = {14,14,14};
// double ptBinsJetsGen[nRadius][30] = {{0.0, 10., 20., 30., 40., 50., 60., 70., 90., 120., 200.},{0.0, 10., 20., 30., 40., 50., 60., 70., 90., 120., 200.},{0.0, 10., 20., 30., 40., 50., 60., 70., 90., 120., 200.}};
// int nBinPtJetsGen[nRadius] = {10,10,10};


// pT binning for jets - hiroki for gen but one more rec bin at 110
double ptBinsJetsRec[nRadius][30] = {{30., 35., 40., 45., 50., 55., 60., 65., 70., 75., 80., 85., 90., 100., 110., 120.},{30., 35., 40., 45., 50., 55., 60., 65., 70., 75., 80., 85., 90., 100., 110., 120.},{30., 35., 40., 45., 50., 55., 60., 65., 70., 75., 80., 85., 90., 100., 110., 120.}};
int nBinPtJetsRec[nRadius] = {15,15,15};
double ptBinsJetsGen[nRadius][30] = {{5.0, 10., 20., 30., 40., 50., 60., 70., 90., 120., 200.},{5.0, 10., 20., 30., 40., 50., 60., 70., 90., 120., 200.},{5.0, 10., 20., 30., 40., 50., 60., 70., 90., 120., 200.}};
int nBinPtJetsGen[nRadius] = {10,10,10};
// double ptBinsJetsGen[nRadius][30] = {{10., 20., 30., 40., 50., 60., 70., 90., 120., 200.},{10., 20., 30., 40., 50., 60., 70., 90., 120., 200.},{10., 20., 30., 40., 50., 60., 70., 90., 120., 200.}};
// int nBinPtJetsGen[nRadius] = {9,9,9};

// // pT binning for jets - hiroki for rec but finer gen binning (shouldn't have finer gen binning for unfolding, as that gen binning is what the unfolded spectrum will have)
// double ptBinsJetsRec[nRadius][30] = {{30., 35., 40., 45., 50., 55., 60., 65., 70., 75., 80., 85., 90., 100., 120.},{30., 35., 40., 45., 50., 55., 60., 65., 70., 75., 80., 85., 90., 100., 120.},{30., 35., 40., 45., 50., 55., 60., 65., 70., 75., 80., 85., 90., 100., 120.}};
// int nBinPtJetsRec[nRadius] = {14,14,14};
// double ptBinsJetsGen[nRadius][22] = {{0.0, 5., 10., 15., 20., 25., 30., 35., 40., 45., 50., 55., 60., 65., 70., 75., 80., 85., 90., 100., 120., 200.},{0.0, 5., 10., 15., 20., 25., 30., 35., 40., 45., 50., 55., 60., 65., 70., 75., 80., 85., 90., 100., 120., 200.},{0.0, 5., 10., 15., 20., 25., 30., 35., 40., 45., 50., 55., 60., 65., 70., 75., 80., 85., 90., 100., 120., 200.}};
// int nBinPtJetsGen[nRadius] = {21,21,21};

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
                            

int nBinPtJetsFine[nRadius] = {195,195,195};
double ptBinsJetsFine[nRadius][201] = {{ 5.,  6.,  7.,  8.,  9.,
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
                                     {   5.,  6.,  7.,  8.,  9.,
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
                                     {   5.,  6.,  7.,  8.,  9.,
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

// Options to be set:


// //////// -------- Full cpass0 Analysis -------- ////////
// TString* texCollisionDataInfo = new TString("Pb-Pb #sqrt{#it{s}} = 5.36 TeV"); 
// const TString* texDatasetsComparisonType = new TString("Period"); 
// const TString* texDatasetsComparisonCommonDenominator = new TString("cpass0 periods");
// const int nDatasets = 10;
// const TString Datasets[nDatasets] = {"LHC23zzh_cpass0","LHC23zzi_cpass0", "LHC23zzf_cpass0", "LHC23zzg_cpass0", "LHC23zx_cpass0", "LHC23zy_cpass0", "LHC23zz_cpass0", "LHC23zza_cpass0", "LHC23zzb_cpass0", "LHC23zze_cpass0"};
// const TString DatasetsNames[nDatasets] = {"LHC23zzh","LHC23zzi", "LHC23zzf", "LHC23zzg", "LHC23zx", "LHC23zy", "LHC23zz", "LHC23zza", "LHC23zzb", "LHC23zze"};
// TFile* file_O2Analysis_list[nDatasets] = {new TFile("Datasets/"+Datasets[0]+"/AnalysisResults.root"),
//                                       new TFile("Datasets/"+Datasets[1]+"/AnalysisResults.root"),
//                                       new TFile("Datasets/"+Datasets[2]+"/AnalysisResults.root"),
//                                       new TFile("Datasets/"+Datasets[3]+"/AnalysisResults.root"),
//                                       new TFile("Datasets/"+Datasets[4]+"/AnalysisResults.root"),
//                                       new TFile("Datasets/"+Datasets[5]+"/AnalysisResults.root"),
//                                       new TFile("Datasets/"+Datasets[6]+"/AnalysisResults.root"),
//                                       new TFile("Datasets/"+Datasets[7]+"/AnalysisResults.root"),
//                                       new TFile("Datasets/"+Datasets[8]+"/AnalysisResults.root"),
//                                       new TFile("Datasets/"+Datasets[9]+"/AnalysisResults.root")
//                                       };
// TString analysisWorkflow = "jet-finder-"+jetType[iJetType]+"-qa";


//////// -------- TestingFile -------- ////////
TString* texCollisionDataInfo = new TString("Pb-Pb #sqrt{#it{s}} = 5.36 TeV"); 
const TString* texDatasetsComparisonType = new TString("");
const TString* texDatasetsComparisonCommonDenominator = new TString("LHC23zzh");
const int nDatasets = 1;
const TString Datasets[nDatasets] = {"LHC23zzh"};
const TString DatasetsNames[nDatasets] = {""};
TFile* file_O2Analysis_list[nDatasets] = {new TFile("Datasets/"+Datasets[0]+"/AnalysisResults.root")
                                      };
// TFile* file_O2Analysis_ppSimDetectorEffect[nDatasets] = new TFile("Datasets/ppSim_LHC23d4/AnalysisResults.root");
TFile* file_O2Analysis_ppSimDetectorEffect = new TFile("Datasets/ppSim_LHC23d4_weighted/AnalysisResults.root");

const TString trainId = "_id12436";
const TString analysisWorkflowData = "jet-finder-charged-qa"+trainId;

const TString analysisWorkflowMC = "jet-finder-charged-qa";