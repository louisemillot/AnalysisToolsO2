
// // Options to be set:
// // TFile* file_O2Analysis = new TFile("AnalysisResults_LHC23d4_fullProdHyperloop.root");
// TFile* file_O2Analysis = new TFile("AnalysisResults.root");
// // TString* texCollisionDataInfo = new TString("#splitline{pp #sqrt{#it{s}} = 13.6 TeV}{Jet-Jet simulation LHC23d4}");
// TString* texCollisionDataInfo = new TString("pp #sqrt{#it{s}} = 13.6 TeV - Jet-Jet sim LHC23d4");

// const int iJetType = 1;
// const int iJetLevel = 2;
// const float PtCutLow  = 20;
// const float PtCutHigh = 40;

// // Permanent Options
// const int nJetType = 4;
// const TString jetType[nJetType] = {"charged", "d0", "neutral", "full"};
// const int nJetLevel = 3;
// const TString jetLevel[nJetLevel] = {"data", "mcd", "mcp"};
// const int nAxisTypes = 3;
// const TString AxisType[nAxisTypes] = {"Kt", "CA", "SD"};
// const TString AxisTypeLegend[nAxisTypes] = {"#it{k}_{T} axes", "C/A axes", "Soft Drop axes"};

// int nBinSubRatio[nJetLevel] = {100, 100, 100};
// double SubRatioBinEdges[nJetLevel][200]; // = {{0},{0},{0}}; // 200 just gives some margin
// const float nSubRatioMax = 1.2;

// // const float drapidity = 1.; // for now I'm not looking at rapidity differential

// // double pTbins[nJetType][20] = {{0.0, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.6, 2.0, 2.4, 3.0},{0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 2.0, 2.4, 3.0},{0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 2.0, 2.4, 3.0}};
// // int nbinpT[nJetType] = {17,8,8};


// Analysis settings
const int nJetType = 3;
const TString jetType[nJetType] = {"charged", "neutral", "full"};
const int nJetLevel = 3;
const TString jetLevel[nJetLevel] = {"data", "mcd", "mcp"};
const int nRadius = 3;
const int nAxisTypes = 3;
const TString AxisType[nAxisTypes] = {"Kt", "CA", "CASD"};
const TString AxisTypeLegend[nAxisTypes] = {"#it{k}_{T} axes", "C/A axes", "C/A SD axes"};

// const TString RadiusLegend[nRadius] = {"R = 0.2", "R = 0.4", "R = 0.6"};
// float arrayRadius[nRadius] = {0.2, 0.4, 0.6};
// const int nRadius = 1;
// const TString RadiusLegend[nRadius] = {"R = 0.4"};
// float arrayRadius[nRadius] = {0.4};
// const int nRadius = 9;
// const TString RadiusLegend[nRadius] = {"R = 0.2", "R = 0.25", "R = 0.3", "R = 0.35", "R = 0.4", "R = 0.45", "R = 0.5", "R = 0.55", "R = 0.6"};
// float arrayRadius[nRadius] = {0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6};

// Choice of jet type (charged, neutral, full) and level (data, detector level, particle level)
const int iJetType = 0;
const int iJetLevel = 0;
const float PtCutLow  = 0;
const float PtCutHigh = 300;


const int nCentralityBins = 7;
const float arrayCentralityBinning[nCentralityBins+1] = {0, 10, 20, 30, 40, 60, 80, 100};
// const int nCentralityBins = 3;
// const float arrayCentralityBinning[nCentralityBins+1] = {0, 10, 50, 90};


int nBinSubRatio[nJetLevel] = {100, 100, 100};
double SubRatioBinEdges[nJetLevel][200]; // = {{0},{0},{0}}; // 200 just gives some margin
const float nSubRatioMax = 1.2;

// const float drapidity = 1.; // for now I'm not looking at rapidity differential

// double pTbins[nJetType][20] = {{0.0, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.6, 2.0, 2.4, 3.0},{0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 2.0, 2.4, 3.0},{0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 2.0, 2.4, 3.0}};
// int nbinpT[nJetType] = {17,8,8};

// TFile* file_AliAnalysis = new TFile("../AnalysisResults_Run2_merged_Jaime.root");


// //////// -------- Inclusive vs D0 -------- ////////
// TString* texCollisionDataInfo = new TString("pp #sqrt{#it{s}} = 13.6 TeV"); 
// const TString* texDatasetsComparisonType = new TString("Inclusive/HF");
// const TString* texDatasetsComparisonCommonDenominator = new TString("LHC22o apass6");
// const int nDatasets = 2;
// const TString Datasets[nDatasets] = {"Inclusive/JE_LHC22o_pass6_small_C_R6_18_F_R6_8", "HF_D0/JE_HF_LHC22_pass6_medium_2P3PDstar_D0C_R4_4"};
// const TString DatasetsNames[nDatasets] = {"inclusive", "D0"};
// TFile* file_O2Analysis_list[nDatasets] = {new TFile("Datasets/"+Datasets[0]+"/AnalysisResults.root"),
//                                         new TFile("Datasets/"+Datasets[1]+"/AnalysisResults.root")
//                                         };
// const TString analysisWorkflow = "jet-nsubjettiness";


//////// -------- Inclusive - before/after cleaning of workflow -------- ////////
TString* texCollisionDataInfo = new TString("pp #sqrt{#it{s}} = 13.6 TeV"); 
const TString* texDatasetsComparisonType = new TString("Inclusive/HF");
const TString* texDatasetsComparisonCommonDenominator = new TString("LHC22o apass6");
const int nDatasets = 2;
const TString Datasets[nDatasets] = {"Inclusive/JE_LHC22o_pass6_small_C_R6_18_F_R6_8", "Inclusive/cleaningOfWorkflow"};
const TString DatasetsNames[nDatasets] = {"beforeCleaning", "afterCleaning"};
TFile* file_O2Analysis_list[nDatasets] = {new TFile("Datasets/"+Datasets[0]+"/AnalysisResults.root"),
                                        new TFile("Datasets/"+Datasets[1]+"/AnalysisResults.root")
                                        };
const TString analysisWorkflow = "jet-nsubjettiness";
