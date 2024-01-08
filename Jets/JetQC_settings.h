
// Bin edge control
float e_binEdge = 0.00001;

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
// const Int_t nRadius = 9;
// const TString RadiusLegend[nRadius] = {"R = 0.2", "R = 0.25", "R = 0.3", "R = 0.35", "R = 0.4", "R = 0.45", "R = 0.5", "R = 0.55", "R = 0.6"};
// Float_t arrayRadius[nRadius] = {0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6};

// Choice of jet type (charged, neutral, full) and level (data, detector level, particle level)
const Int_t iJetType = 0;
const Int_t iJetLevel = 0;


// Options to be set:
// //////// -------- Full cpass0 Analysis -------- ////////
// TString* texCollisionDataInfo = new TString("Pb-Pb #sqrt{#it{s}} = 5.36 TeV"); 
// const TString* texDatasetsComparisonType = new TString("Period"); 
// const TString* texDatasetsComparisonCommonDenominator = new TString("cpass0 periods");
// const Int_t nDatasets = 10;
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

// //////// -------- Flat Phi Periods -------- ////////
// TString* texCollisionDataInfo = new TString("Pb-Pb #sqrt{#it{s}} = 5.36 TeV"); 
// const TString* texDatasetsComparisonType = new TString("Period"); 
// const TString* texDatasetsComparisonCommonDenominator = new TString("good Phi cpass0 periods");
// const Int_t nDatasets = 4;
// const TString Datasets[nDatasets] = {"LHC23zzh_cpass0", "LHC23zzf_cpass0", "LHC23zzg_cpass0", "LHC23zzi_cpass0"};
// const TString DatasetsNames[nDatasets] = {"LHC23zzh", "LHC23zzf", "LHC23zzg", "LHC23zzi"};
// TFile* file_O2Analysis_list[nDatasets] = {new TFile("Datasets/"+Datasets[0]+"/AnalysisResults.root"),
//                                       new TFile("Datasets/"+Datasets[1]+"/AnalysisResults.root"),
//                                       new TFile("Datasets/"+Datasets[2]+"/AnalysisResults.root"),
//                                       new TFile("Datasets/"+Datasets[3]+"/AnalysisResults.root")
//                                       };
// const TString analysisWorkflow = "jet-finder-charged-qa";

// //////// -------- Hole Phi Periods -------- ////////
// TString* texCollisionDataInfo = new TString("Pb-Pb #sqrt{#it{s}} = 5.36 TeV"); 
// const TString* texDatasetsComparisonType = new TString("Period"); 
// const TString* texDatasetsComparisonCommonDenominator = new TString("bad Phi cpass0 periods");
// const Int_t nDatasets = 6;
// const TString Datasets[nDatasets] = {"LHC23zx_cpass0", "LHC23zy_cpass0", "LHC23zz_cpass0", "LHC23zza_cpass0", "LHC23zzb_cpass0", "LHC23zze_cpass0"};
// const TString DatasetsNames[nDatasets] = {"LHC23zx", "LHC23zy", "LHC23zz", "LHC23zza", "LHC23zzb", "LHC23zze"};
// TFile* file_O2Analysis_list[nDatasets] = {new TFile("Datasets/"+Datasets[0]+"/AnalysisResults.root"),
//                                       new TFile("Datasets/"+Datasets[1]+"/AnalysisResults.root"),
//                                       new TFile("Datasets/"+Datasets[2]+"/AnalysisResults.root"),
//                                       new TFile("Datasets/"+Datasets[3]+"/AnalysisResults.root"),
//                                       new TFile("Datasets/"+Datasets[4]+"/AnalysisResults.root"),
//                                       new TFile("Datasets/"+Datasets[5]+"/AnalysisResults.root")
//                                       };
// const TString analysisWorkflow = "jet-finder-charged-qa";

// //////// -------- Cpass comparison -------- ////////
// TString* texCollisionDataInfo = new TString("Pb-Pb #sqrt{#it{s}} = 5.36 TeV"); 
// const TString* texDatasetsComparisonType = new TString("Cpass");
// const TString* texDatasetsComparisonCommonDenominator = new TString("LHC23zzh run 544122");
// const Int_t nDatasets = 5;
// const TString Datasets[nDatasets] = {"LHC23zzh_cpass11", "LHC23zzh_cpass10", "LHC23zzh_cpass8", "LHC23zzh_cpass1", "LHC23zzh_cpass0"};
// const TString DatasetsNames[nDatasets] = {"cpass11", "cpass10", "cpass8", "cpass1", "cpass0"};
// TFile* file_O2Analysis_list[nDatasets] = {new TFile("Datasets/"+Datasets[0]+"/AnalysisResults.root"),
//                                       new TFile("Datasets/"+Datasets[1]+"/AnalysisResults.root"),
//                                       new TFile("Datasets/"+Datasets[2]+"/AnalysisResults.root"),
//                                       new TFile("Datasets/"+Datasets[3]+"/AnalysisResults.root"),
//                                       new TFile("Datasets/"+Datasets[4]+"/AnalysisResults.root")
//                                       };
// const TString analysisWorkflow = "jet-finder-charged-qa";

// //////// -------- Run Comparison - LHC23zzh -------- ////////
// TString* texCollisionDataInfo = new TString("Pb-Pb #sqrt{#it{s}} = 5.36 TeV"); 
// const TString* texDatasetsComparisonType = new TString("Run"); 
// const TString* texDatasetsComparisonCommonDenominator = new TString("LHC23zzh_cpass0");
// const Int_t nDatasets = 8;
// const TString Datasets[nDatasets] = {"LHC23zzh_run544122","LHC23zzh_run544126","LHC23zzh_run544124","LHC23zzh_run544121","LHC23zzh_run544116","LHC23zzh_run544098","LHC23zzh_run544095","LHC23zzh_run544091"};
// const TString DatasetsNames[nDatasets] = {"run 544122","run 544126","run 544124","run 544121","run 544116","run 544098","run 544095","run 544091"};
// TFile* file_O2Analysis_list[nDatasets] = {new TFile("Datasets/"+Datasets[0]+"/AnalysisResults.root"),
//                                       new TFile("Datasets/"+Datasets[1]+"/AnalysisResults.root"),
//                                       new TFile("Datasets/"+Datasets[2]+"/AnalysisResults.root"),
//                                       new TFile("Datasets/"+Datasets[3]+"/AnalysisResults.root"),
//                                       new TFile("Datasets/"+Datasets[4]+"/AnalysisResults.root"),
//                                       new TFile("Datasets/"+Datasets[5]+"/AnalysisResults.root"),
//                                       new TFile("Datasets/"+Datasets[6]+"/AnalysisResults.root"),
//                                       new TFile("Datasets/"+Datasets[7]+"/AnalysisResults.root")
//                                       };
// TString analysisWorkflow = "jet-finder-"+jetType[iJetType]+"-qa";

// //////// -------- Track Eta cut comparison -------- ////////
// TString* texCollisionDataInfo = new TString("Pb-Pb #sqrt{#it{s}} = 5.36 TeV"); 
// const TString* texDatasetsComparisonType = new TString("#eta_{track} cut");
// const TString* texDatasetsComparisonCommonDenominator = new TString("LHC23zzh cpass 11");
// const Int_t nDatasets = 5;
// const TString Datasets[nDatasets] = {"LHC23zzh_eta0.9", "LHC23zzh_eta0.8", "LHC23zzh_eta0.7", "LHC23zzh_eta0.6", "LHC23zzh_eta0.5"};
// const TString DatasetsNames[nDatasets] = {"-0.9 < #eta_{Tracks} < 0.9", "-0.8 < #eta_{Tracks} < 0.8", "-0.7 < #eta_{Tracks} < 0.7", "-0.6 < #eta_{Tracks} < 0.6", "-0.5 < #eta_{Tracks} < 0.5"};
// TFile* file_O2Analysis_list[nDatasets] = {new TFile("Datasets/"+Datasets[0]+"/AnalysisResults.root"),
//                                       new TFile("Datasets/"+Datasets[1]+"/AnalysisResults.root"),
//                                       new TFile("Datasets/"+Datasets[2]+"/AnalysisResults.root"),
//                                       new TFile("Datasets/"+Datasets[3]+"/AnalysisResults.root"),
//                                       new TFile("Datasets/"+Datasets[4]+"/AnalysisResults.root")
//                                       };
// const TString analysisWorkflow = "jet-finder-charged-qa";

// //////// -------- Jet Radius comparison -------- ////////
// TString* texCollisionDataInfo = new TString("Pb-Pb #sqrt{#it{s}} = 5.36 TeV"); 
// const TString* texDatasetsComparisonType = new TString("");
// const TString* texDatasetsComparisonCommonDenominator = new TString("LHC23zzh cpass 11");
// const Int_t nDatasets = 1;
// const TString Datasets[nDatasets] = {"LHC23zzh_cpass11"};
// const TString DatasetsNames[nDatasets] = {""};
// TFile* file_O2Analysis_list[nDatasets] = {new TFile("Datasets/"+Datasets[0]+"/AnalysisResults.root")
//                                       };
// const TString analysisWorkflow = "jet-finder-charged-qa";

// //////// -------- TRD tracks check -------- ////////
// TString* texCollisionDataInfo = new TString("Pb-Pb #sqrt{#it{s}} = 5.36 TeV"); 
// const TString* texDatasetsComparisonType = new TString("");
// const TString* texDatasetsComparisonCommonDenominator = new TString("LHC23zzh");
// const Int_t nDatasets = 1;
// const TString Datasets[nDatasets] = {"LHC23zzh"};
// const TString DatasetsNames[nDatasets] = {""};
// TFile* file_O2Analysis_list[nDatasets] = {new TFile("Datasets/"+Datasets[0]+"/AnalysisResults.root")
//                                       };
// const TString analysisWorkflow = "jet-finder-charged-qa";

//////// -------- Track Se check -------- ////////
TString* texCollisionDataInfo = new TString("Pb-Pb #sqrt{#it{s}} = 5.36 TeV"); 
const TString* texDatasetsComparisonType = new TString("TrackSel");
const TString* texDatasetsComparisonCommonDenominator = new TString("LHC23zzh");
const Int_t nDatasets = 3;
const TString Datasets[nDatasets] = {"LHC23zzh_AntonTrackSel", "LHC23zzh_GlobalWithOptionalTpc", "LHC23zzh_GlobalDefault"};
const TString DatasetsNames[nDatasets] = {"AntonTrackSel", "GlobalWithOptionalTpc", "GlobalDefault"};
TFile* file_O2Analysis_list[nDatasets] = {new TFile("Datasets/"+Datasets[0]+"/AnalysisResults.root"),
                                        new TFile("Datasets/"+Datasets[1]+"/AnalysisResults.root"),
                                        new TFile("Datasets/"+Datasets[2]+"/AnalysisResults.root")
                                        };
const TString analysisWorkflow = "jet-finder-charged-qa";


