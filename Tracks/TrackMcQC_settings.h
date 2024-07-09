
// Bin edge control
// float GLOBAL_epsilon = 0.00001;

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

const double deltaEtaMcVsTrackEfficiency = 0;



// //////// -------- Track Sel check rough -------- ////////
// TString* texCollisionDataInfo = new TString("Pb-Pb #sqrt{#it{s}} = 5.36 TeV"); 
// const TString* texDatasetsComparisonType = new TString("TrackSel");
// const TString* texDatasetsComparisonCommonDenominator = new TString("LHC23zzh");
// const int nDatasets = 3;
// const TString Datasets[nDatasets] = {"LHC23zzh_AntonTrackSel", "LHC23zzh_GlobalWithOptionalTpc", "LHC23zzh_GlobalDefault"};
// const TString DatasetsNames[nDatasets] = {"MMTrackSel", "GlobalWithOptionalTpc", "GlobalTrackSel"};
// TFile* file_O2Analysis_list[nDatasets] = {new TFile("Datasets/"+Datasets[0]+"/AnalysisResults.root"),
//                                         new TFile("Datasets/"+Datasets[1]+"/AnalysisResults.root"),
//                                         new TFile("Datasets/"+Datasets[2]+"/AnalysisResults.root")
//                                         };
// const TString analysisWorkflow = "jet-finder-charged-qa";




// //////// -------- Track Sel check Detailed-------- ////////
// TString* texCollisionDataInfo = new TString("Pb-Pb #sqrt{#it{s}} = 5.36 TeV"); 
// const TString* texDatasetsComparisonType = new TString("TrackSel");
// const TString* texDatasetsComparisonCommonDenominator = new TString("LHC23zzh");
// const int nDatasets = 6;
// const TString Datasets[nDatasets] = {"LHC23zzh_AntonTrackSel", "LHC23zzh_intermTrackSel2", "LHC23zzh_intermTrackSel1bis", "LHC23zzh_intermTrackSel1", "LHC23zzh_GlobalTpcOptional", "LHC23zzh_GlobalDefault"};
// const TString DatasetsNames[nDatasets] = {"MMTrackSel", "intermTrackSel2", "intermTrackSel1bis", "intermTrackSel1", "GlobalTpcOptional", "GlobalTrackSel"};
// TFile* file_O2Analysis_list[nDatasets] = {new TFile("Datasets/"+Datasets[0]+"/AnalysisResults.root"),
//                                         new TFile("Datasets/"+Datasets[1]+"/AnalysisResults.root"),
//                                         new TFile("Datasets/"+Datasets[2]+"/AnalysisResults.root"),
//                                         new TFile("Datasets/"+Datasets[3]+"/AnalysisResults.root"),
//                                         new TFile("Datasets/"+Datasets[4]+"/AnalysisResults.root"),
//                                         new TFile("Datasets/"+Datasets[5]+"/AnalysisResults.root")
//                                         };
// const TString analysisWorkflow = "jet-finder-charged-qa";

// //////// -------- Track Sel check small -------- ////////
// TString* texCollisionDataInfo = new TString("Pb-Pb #sqrt{#it{s}} = 5.36 TeV"); 
// const TString* texDatasetsComparisonType = new TString("TrackSel");
// const TString* texDatasetsComparisonCommonDenominator = new TString("LHC23zzh");
// const int nDatasets = 2;
// const TString Datasets[nDatasets] = {"LHC23zzh_MMTracks", "LHC23zzh_GlobalTracks"};
// const TString DatasetsNames[nDatasets] = {"MMTrackSel", "GlobalTrackSel"};
// TFile* file_O2Analysis_list[nDatasets] = {new TFile("Datasets/"+Datasets[0]+"/AnalysisResults.root"),
//                                         new TFile("Datasets/"+Datasets[1]+"/AnalysisResults.root")
//                                         };
// const TString analysisWorkflow = "track-efficiency";


// //////// -------- Track Sel check three -------- ////////
// TString* texCollisionDataInfo = new TString("Pb-Pb #sqrt{#it{s}} = 5.36 TeV"); 
// const TString* texDatasetsComparisonType = new TString("TrackSel");
// const TString* texDatasetsComparisonCommonDenominator = new TString("LHC23zzh");
// const int nDatasets = 3;
// const TString Datasets[nDatasets] = {"LHC23zzh_MMTracks", "LHC23zzh_GlobalTpcOptional", "LHC23zzh_GlobalTracks"};
// const TString DatasetsNames[nDatasets] = {"MMTrackSel", "GlobalTpcOptional", "GlobalTrackSel"};
// TFile* file_O2Analysis_list[nDatasets] = {new TFile("Datasets/"+Datasets[0]+"/AnalysisResults.root"),
//                                         new TFile("Datasets/"+Datasets[1]+"/AnalysisResults.root"),
//                                         new TFile("Datasets/"+Datasets[2]+"/AnalysisResults.root")
//                                         };
// const TString analysisWorkflow = "track-efficiency";



// //////// -------- Apass vs Cpasses comparison -------- ////////
// TString* texCollisionDataInfo = new TString("Pb-Pb #sqrt{#it{s}} = 5.36 TeV"); 
// const TString* texDatasetsComparisonType = new TString("Cpass");
// const TString* texDatasetsComparisonCommonDenominator = new TString("LHC23zzh run 544122");
// const int nDatasets = 6;
// const TString Datasets[nDatasets] = {"LHC23zzh_apass1", "LHC23zzh_cpass11", "LHC23zzh_cpass10", "LHC23zzh_cpass8", "LHC23zzh_cpass1", "LHC23zzh_cpass0"};
// const TString DatasetsNames[nDatasets] = {"apass1", "cpass11", "cpass10", "cpass8", "cpass1", "cpass0"};
// TFile* file_O2Analysis_list[nDatasets] = {new TFile("Datasets/"+Datasets[0]+"/AnalysisResults.root"),
//                                       new TFile("Datasets/"+Datasets[1]+"/AnalysisResults.root"),
//                                       new TFile("Datasets/"+Datasets[2]+"/AnalysisResults.root"),
//                                       new TFile("Datasets/"+Datasets[3]+"/AnalysisResults.root"),
//                                       new TFile("Datasets/"+Datasets[4]+"/AnalysisResults.root"),
//                                       new TFile("Datasets/"+Datasets[5]+"/AnalysisResults.root")
//                                       };
// const TString analysisWorkflow = "jet-finder-charged-qa";


// //////// -------- Apass - single -------- ////////
// TString* texCollisionDataInfo = new TString("Pb-Pb #sqrt{#it{s}} = 5.36 TeV"); 
// const TString* texDatasetsComparisonType = new TString("Apass2");
// const TString* texDatasetsComparisonCommonDenominator = new TString("LHC23zzh run 544122");
// const int nDatasets = 1;
// const TString Datasets[nDatasets] = {"LHC23zzh_run544122"};
// const TString DatasetsNames[nDatasets] = {"apass2"};
// TFile* file_O2Analysis_list[nDatasets] = {new TFile("Datasets/"+Datasets[0]+"/AnalysisResults.root")
//                                         };
// const TString analysisWorkflow = "jet-finder-charged-qa";


// //////// -------- Track efficeincy calculation comparison -------- ////////
// TString* texCollisionDataInfo = new TString("Pb-Pb #sqrt{#it{s}} = 5.36 TeV"); 
// const TString* texDatasetsComparisonType = new TString("EfficiencyCalc");
// const TString* texDatasetsComparisonCommonDenominator = new TString("LHC23k6d");
// const int nDatasets = 2;
// const TString Datasets[nDatasets] = {"LHC23k6d_AllChecksAndSel8Global", "LHC23k6d_AllChecksAndSel8"};
// const TString DatasetsNames[nDatasets] = {"globalTracks", "MMsel"};
// TFile* file_O2Analysis_list[nDatasets] = {new TFile("Datasets/"+Datasets[0]+"/AnalysisResults.root"),
//                                         new TFile("Datasets/"+Datasets[1]+"/AnalysisResults.root")
//                                         };
// const TString analysisWorkflow = "track-efficiency";

// //////// -------- Track efficeincy calculation comparison -------- ////////
// TString* texCollisionDataInfo = new TString("Pb-Pb #sqrt{#it{s}} = 5.36 TeV"); 
// const TString* texDatasetsComparisonType = new TString("EfficiencyCalc");
// const TString* texDatasetsComparisonCommonDenominator = new TString("LHC23k6d");
// const int nDatasets = 4;
// const TString Datasets[nDatasets] = {"LHC23zzh_MMsel_Cent0010", "LHC23zzh_MMsel_Cent5090", "LHC23zzh_Global_Cent0010", "LHC23zzh_Global_Cent5090"};
// const TString DatasetsNames[nDatasets] = {"MMsel 0-10%", "MMsel 50-90%", "globalTracks 0-10%", "globalTracks 50-90%"};
// TFile* file_O2Analysis_list[nDatasets] = {new TFile("Datasets/"+Datasets[0]+"/AnalysisResults.root"),
//                                         new TFile("Datasets/"+Datasets[1]+"/AnalysisResults.root"),
//                                         new TFile("Datasets/"+Datasets[2]+"/AnalysisResults.root"),
//                                         new TFile("Datasets/"+Datasets[3]+"/AnalysisResults.root")
//                                         };
// const TString analysisWorkflow = "track-efficiency";



// //////// -------- Track efficeincy calculation comparison -------- ////////
// TString* texCollisionDataInfo = new TString("Pb-Pb #sqrt{#it{s}} = 5.36 TeV"); 
// const TString* texDatasetsComparisonType = new TString("EfficiencyCalc");
// const TString* texDatasetsComparisonCommonDenominator = new TString("LHC24d2");
// const int nDatasets = 2;
// const TString Datasets[nDatasets] = {"LHC24d2", "LHC24d2"};
// const TString DatasetsNames[nDatasets] = {"LHC24d2 - globalTracks", "LHC24d2 - uniformTracks"};
// TFile* file_O2Analysis_list[nDatasets] = {new TFile("Datasets/"+Datasets[0]+"/AnalysisResults.root"),
//                                         new TFile("Datasets/"+Datasets[1]+"/AnalysisResults.root")
//                                         };
// const TString analysisWorkflow[nDatasets] = {"track-efficiency",
//                                           "track-efficiency_uniformTracks"
//                                           };


// //////// -------- Track efficiency calculation centrality comparison -------- ////////
// TString* texCollisionDataInfo = new TString("Pb-Pb #sqrt{#it{s}} = 5.36 TeV"); 
// const TString* texDatasetsComparisonType = new TString("EfficiencyCalc");
// const TString* texDatasetsComparisonCommonDenominator = new TString("LHC23k6d");
// const int nDatasets = 4;
// const TString Datasets[nDatasets] = {"LHC23k6d_Global_cent0010", "LHC23k6d_Global_cent5090", "LHC23k6d_Uniform_cent0010", "LHC23k6d_Uniform_cent5090"};
// const TString DatasetsNames[nDatasets] = {"Global 0-10%", "Global 50-90%", "Uniform 0-10%", "Uniform 50-90%"};
// TFile* file_O2Analysis_list[nDatasets] = {new TFile("Datasets/"+Datasets[0]+"/AnalysisResults.root"),
//                                         new TFile("Datasets/"+Datasets[1]+"/AnalysisResults.root"),
//                                         new TFile("Datasets/"+Datasets[2]+"/AnalysisResults.root"),
//                                         new TFile("Datasets/"+Datasets[3]+"/AnalysisResults.root")
//                                         };
// const TString analysisWorkflow[nDatasets] = {"track-efficiency",
//                                             "track-efficiency",
//                                             "track-efficiency",
//                                             "track-efficiency"
//                                           };


// //////// -------- LHC24b1 vs ppJetAnchored -------- ////////
// TString* texCollisionDataInfo = new TString("Pb-Pb #sqrt{#it{s}} = 5.36 TeV"); 
// const TString* texDatasetsComparisonType = new TString("TrackSel");
// const TString* texDatasetsComparisonCommonDenominator = new TString("LHC23zzh_apass3");
// const int nDatasets = 3;
// const TString Datasets[nDatasets] = {"AnalysisResults_LHC24b1_noTrackSel", "AnalysisResults_ppJets_likeLHC24b1_noTrackSel", "AnalysisResults_unanchored_ppJets"};
// const TString DatasetsNames[nDatasets] = {"LHC24b1", "ppJetsAnchored", "ppJetsUnanchored"};
// TFile* file_O2Analysis_list[nDatasets] = {new TFile(Datasets[0]+".root"),
//                                         new TFile(Datasets[1]+".root"),
//                                         new TFile(Datasets[2]+".root")
//                                         };
// const TString analysisWorkflow[nDatasets] = {"qa-event-track/Tracks",
//                                           "qa-event-track/Tracks",
//                                           "qa-event-track/Tracks"
//                                           };

// //////// -------- Track efficiency calculation dcaz cut comparison - avec DCAz cut 2cm-------- ////////
// TString* texCollisionDataInfo = new TString("Pb-Pb #sqrt{#it{s}} = 5.36 TeV"); 
// const TString* texDatasetsComparisonType = new TString("EfficiencyCalc");
// const TString* texDatasetsComparisonCommonDenominator = new TString("LHC24d2");
// const int nDatasets = 15;
// const TString Datasets[nDatasets] = {"LHC24d2_3", "LHC24d2_3", "LHC24d2_3", "LHC24d2_2", "LHC24d2_1", "LHC24d2_2", "LHC24d2_1", "LHC24d2_2", "LHC24d2_1", "LHC24d2_3", "LHC24d2_3", "LHC24d2_3", "LHC24d2_2", "LHC24d2_2", "LHC24d2_1"};
// const TString DatasetsNames[nDatasets] = {"uniform DCAz 0.005", "uniform DCAz 0.01", "uniform DCAz 0.03", "uniform DCAz 0.05", "uniform DCAz 0.1", "uniform DCAz 0.2", "uniform DCAz 0.5", "uniform DCAz 1.0", "uniform DCAz 2.0", "global DCAz 0.005", "global DCAz 0.01", "global DCAz 0.03", "global DCAz 0.05", "global DCAz 1.0", "global DCAz 2.0"};
// TFile* file_O2Analysis_list[nDatasets] = {new TFile("Datasets/"+Datasets[0]+"/AnalysisResults.root"),
//                                         new TFile("Datasets/"+Datasets[1]+"/AnalysisResults.root"),
//                                         new TFile("Datasets/"+Datasets[2]+"/AnalysisResults.root"),
//                                         new TFile("Datasets/"+Datasets[3]+"/AnalysisResults.root"),
//                                         new TFile("Datasets/"+Datasets[4]+"/AnalysisResults.root"),
//                                         new TFile("Datasets/"+Datasets[5]+"/AnalysisResults.root"),
//                                         new TFile("Datasets/"+Datasets[6]+"/AnalysisResults.root"),
//                                         new TFile("Datasets/"+Datasets[7]+"/AnalysisResults.root"),
//                                         new TFile("Datasets/"+Datasets[8]+"/AnalysisResults.root"),
//                                         new TFile("Datasets/"+Datasets[9]+"/AnalysisResults.root"),
//                                         new TFile("Datasets/"+Datasets[10]+"/AnalysisResults.root"),
//                                         new TFile("Datasets/"+Datasets[11]+"/AnalysisResults.root"),
//                                         new TFile("Datasets/"+Datasets[12]+"/AnalysisResults.root"),
//                                         new TFile("Datasets/"+Datasets[13]+"/AnalysisResults.root"),
//                                         new TFile("Datasets/"+Datasets[14]+"/AnalysisResults.root")
//                                         };
// const TString analysisWorkflow[nDatasets] = {"track-efficiency_uniformTracks_DCAz0005",
//                                             "track-efficiency_uniformTracks_DCAxyz001",
//                                             "track-efficiency_uniformTracks_DCAz003",
//                                             "track-efficiency_uniformTracks_DCAxy005_id13159",
//                                             "track-efficiency_uniformTracks_DCAxy01_id13159",
//                                             "track-efficiency_uniformTracks_DCAxy02_id13159",
//                                             "track-efficiency_uniformTracks_DCAxy05_id13159",
//                                             "track-efficiency_uniformTracks_DCAxy1_id13159",
//                                             "track-efficiency_id13159",
//                                             "track-efficiency_globalTracks_DCAz0005",
//                                             "track-efficiency_globalTracks_DCAz001",
//                                             "track-efficiency_globalTracks_DCAz003",
//                                             "track-efficiency_globalTracks_DCAxy005_id13171",
//                                             "track-efficiency_globalTracks_DCAxy1_id13171",
//                                             "track-efficiency_id13171",
//                                           };


//////// -------- Track efficiency calculation dcaz cut comparison - avec DCAz cut sans 2cm-------- ////////
TString* texCollisionDataInfo = new TString("Pb-Pb #sqrt{#it{s}} = 5.36 TeV"); 
const TString* texDatasetsComparisonType = new TString("EfficiencyCalc");
const TString* texDatasetsComparisonCommonDenominator = new TString("LHC24d2");
const int nDatasets = 13;
const TString Datasets[nDatasets] = {"LHC24d2_3", "LHC24d2_3", "LHC24d2_3", "LHC24d2_2", "LHC24d2_1", "LHC24d2_2", "LHC24d2_1", "LHC24d2_2", "LHC24d2_3", "LHC24d2_3", "LHC24d2_3", "LHC24d2_2", "LHC24d2_2"};
const TString DatasetsNames[nDatasets] = {"uniform DCAz 0.005", "uniform DCAz 0.01", "uniform DCAz 0.03", "uniform DCAz 0.05", "uniform DCAz 0.1", "uniform DCAz 0.2", "uniform DCAz 0.5", "uniform DCAz 1.0", "global DCAz 0.005", "global DCAz 0.01", "global DCAz 0.03", "global DCAz 0.05", "global DCAz 1.0"};
TFile* file_O2Analysis_list[nDatasets] = {new TFile("Datasets/"+Datasets[0]+"/AnalysisResults.root"),
                                        new TFile("Datasets/"+Datasets[1]+"/AnalysisResults.root"),
                                        new TFile("Datasets/"+Datasets[2]+"/AnalysisResults.root"),
                                        new TFile("Datasets/"+Datasets[3]+"/AnalysisResults.root"),
                                        new TFile("Datasets/"+Datasets[4]+"/AnalysisResults.root"),
                                        new TFile("Datasets/"+Datasets[5]+"/AnalysisResults.root"),
                                        new TFile("Datasets/"+Datasets[6]+"/AnalysisResults.root"),
                                        new TFile("Datasets/"+Datasets[7]+"/AnalysisResults.root"),
                                        new TFile("Datasets/"+Datasets[8]+"/AnalysisResults.root"),
                                        new TFile("Datasets/"+Datasets[9]+"/AnalysisResults.root"),
                                        new TFile("Datasets/"+Datasets[10]+"/AnalysisResults.root"),
                                        new TFile("Datasets/"+Datasets[11]+"/AnalysisResults.root"),
                                        new TFile("Datasets/"+Datasets[12]+"/AnalysisResults.root")
                                        };
const TString analysisWorkflow[nDatasets] = {"track-efficiency_uniformTracks_DCAz0005",
                                            "track-efficiency_uniformTracks_DCAxyz001",
                                            "track-efficiency_uniformTracks_DCAz003",
                                            "track-efficiency_uniformTracks_DCAxy005_id13159",
                                            "track-efficiency_uniformTracks_DCAxy01_id13159",
                                            "track-efficiency_uniformTracks_DCAxy02_id13159",
                                            "track-efficiency_uniformTracks_DCAxy05_id13159",
                                            "track-efficiency_uniformTracks_DCAxy1_id13159",
                                            "track-efficiency_globalTracks_DCAz0005",
                                            "track-efficiency_globalTracks_DCAz001",
                                            "track-efficiency_globalTracks_DCAz003",
                                            "track-efficiency_globalTracks_DCAxy005_id13171",
                                            "track-efficiency_globalTracks_DCAxy1_id13171"
                                          };