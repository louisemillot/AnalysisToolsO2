
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

//////// -------- Track Sel check -------- ////////
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
