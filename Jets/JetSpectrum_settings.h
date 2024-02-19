

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
const int nCentralityBins = 3;
const float arrayCentralityBinning[nCentralityBins+1] = {0, 10, 50, 90};

// pT binning for jets
double ptBinsJetsRec[nRadius][20] = {{1.0, 5., 10., 15., 20., 25., 30., 40., 50., 60., 70., 80., 100., 120., 140., 200.},{1.0, 5., 10., 15., 20., 25., 30., 40., 50., 60., 70., 80., 100., 120., 140., 200.},{1.0, 5., 10., 15., 20., 25., 30., 40., 50., 60., 70., 80., 100., 120., 140., 200.}};
int nBinPtJetsRec[nRadius] = {15,15,15};
double ptBinsJetsGen[nRadius][20] = {{1.0, 5., 10., 15., 20., 25., 30., 40., 50., 60., 70., 80., 100., 120., 140., 200.},{1.0, 5., 10., 15., 20., 25., 30., 40., 50., 60., 70., 80., 100., 120., 140., 200.},{1.0, 5., 10., 15., 20., 25., 30., 40., 50., 60., 70., 80., 100., 120., 140., 200.}};
int nBinPtJetsGen[nRadius] = {15,15,15};

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
const TString analysisWorkflow = "jet-finder-charged-qa";