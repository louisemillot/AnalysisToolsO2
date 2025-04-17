// To avoid your local _inputs.h files to be replaced by the template in the git repository, you should ask git to ignore them with git update-index --assume-unchanged (for more details see https://stackoverflow.com/questions/3319479/can-i-git-commit-a-file-and-ignore-its-content-changes)

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


// //////// -------- Track efficiency calculation dcaz cut comparison - avec DCAz cut sans 2cm-------- ////////
// TString* texCollisionDataInfo = new TString("Pb-Pb #sqrt{#it{s}} = 5.36 TeV"); 
// const TString* texDatasetsComparisonType = new TString("EfficiencyCalc");
// const TString* texDatasetsComparisonCommonDenominator = new TString("LHC24d2");
// const int nDatasets = 13;
// const TString Datasets[nDatasets] = {"LHC24d2_3", "LHC24d2_3", "LHC24d2_3", "LHC24d2_2", "LHC24d2_1", "LHC24d2_2", "LHC24d2_1", "LHC24d2_2", "LHC24d2_3", "LHC24d2_3", "LHC24d2_3", "LHC24d2_2", "LHC24d2_2"};
// const TString DatasetsNames[nDatasets] = {"uniform DCAz 0.005", "uniform DCAz 0.01", "uniform DCAz 0.03", "uniform DCAz 0.05", "uniform DCAz 0.1", "uniform DCAz 0.2", "uniform DCAz 0.5", "uniform DCAz 1.0", "global DCAz 0.005", "global DCAz 0.01", "global DCAz 0.03", "global DCAz 0.05", "global DCAz 1.0"};
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
//                                         new TFile("Datasets/"+Datasets[12]+"/AnalysisResults.root")
//                                         };
// const TString analysisWorkflow[nDatasets] = {"track-efficiency_uniformTracks_DCAz0005",
//                                             "track-efficiency_uniformTracks_DCAxyz001",
//                                             "track-efficiency_uniformTracks_DCAz003",
//                                             "track-efficiency_uniformTracks_DCAxy005_id13159",
//                                             "track-efficiency_uniformTracks_DCAxy01_id13159",
//                                             "track-efficiency_uniformTracks_DCAxy02_id13159",
//                                             "track-efficiency_uniformTracks_DCAxy05_id13159",
//                                             "track-efficiency_uniformTracks_DCAxy1_id13159",
//                                             "track-efficiency_globalTracks_DCAz0005",
//                                             "track-efficiency_globalTracks_DCAz001",
//                                             "track-efficiency_globalTracks_DCAz003",
//                                             "track-efficiency_globalTracks_DCAxy005_id13171",
//                                             "track-efficiency_globalTracks_DCAxy1_id13171"
//                                           };

                                        

// //////// -------- Track efficiency calculation dcaz cut comparison - avec DCAz cut sans 2cm-------- ////////
// TString* texCollisionDataInfo = new TString("Pb-Pb #sqrt{#it{s}} = 5.36 TeV"); 
// const TString* texDatasetsComparisonType = new TString("EfficiencyCalc");
// const TString* texDatasetsComparisonCommonDenominator = new TString("LHC24d2");
// const int nDatasets = 2;
// const TString Datasets[nDatasets] = {"LHC24d2_2", "LHC24d2_2"};
// const TString DatasetsNames[nDatasets] = {"uniform", "global"};
// TFile* file_O2Analysis_list[nDatasets] = {new TFile("Datasets/"+Datasets[0]+"/AnalysisResults.root"),
//                                         new TFile("Datasets/"+Datasets[1]+"/AnalysisResults.root")
//                                         };
// const TString analysisWorkflow[nDatasets] = {"track-efficiency_uniformTracks_DCAxy1_id13159",
//                                             "track-efficiency_globalTracks_DCAxy1_id13171"
//                                           };


// //////// -------- Run 3 MC local -------- ////////
// TString* texCollisionDataInfo = new TString("Pb-Pb #sqrt{#it{s}} = 5.36 TeV"); 
// const TString* texDatasetsComparisonType = new TString("global tracks");
// const TString* texDatasetsComparisonCommonDenominator = new TString("LHC24d2 - global tracks");
// const int nDatasets = 1;
// const TString Datasets[nDatasets] = {"LHC24d2_global"};
// const TString DatasetsNames[nDatasets] = {"global"};
// TFile* file_O2Analysis_list[nDatasets] = {new TFile("Datasets/"+Datasets[0]+"/AnalysisResults.root")
//                                         };
// const TString analysisWorkflow[nDatasets] = {"track-efficiency_globalTracks_DCAz20"
//                                           };



// //////// -------- MC pp anchored PbPB vs pp -------- ////////
// TString* texCollisionDataInfo = new TString("pp #sqrt{#it{s}} = 13.6 TeV");
// const TString* texDatasetsComparisonType = new TString("anchor");
// const TString* texDatasetsComparisonCommonDenominator = new TString("jet-jet sim");
// const int nDatasets = 3;
// // const TString Datasets[nDatasets] = {"LHC23zzh_apass3_global", "LHC23zzh_apass3_uniform", "LHC23zzh_apass3_itsOnly", "LHC23zzh_apass3_globalWithoutTpcCrossedRows"};
// // const TString Datasets[nDatasets] = {"anchoredJetJetPbPb", "anchoredJetJetpp", "LHC24f3_train240962", "LHC24f3_local", "unanchoredJetJet_train230486"};
// const TString Datasets[nDatasets] = {"anchoredJetJetPbPb", "anchoredJetJetpp", "LHC24f3_local"};
// const TString DatasetsNames[nDatasets] = {"PbPbAnchored", "ppAnchored", "LHC24f3 pp"};
// TFile* file_O2Analysis_list[nDatasets] = {new TFile("Datasets/"+Datasets[0]+"/AnalysisResults.root"),
//                                         new TFile("Datasets/"+Datasets[1]+"/AnalysisResults.root"),
//                                         new TFile("Datasets/"+Datasets[2]+"/AnalysisResults.root")
//                                         };
// const TString analysisWorkflow[nDatasets] = {"track-efficiency",
//                                             "track-efficiency",
//                                             "track-efficiency"
//                                           };

// const TString trainId = "";
// const bool isDatasetWeighted[nDatasets] = {true, true, false};

// // jaime testJetAnchoredGapGen-20240801-155530

// //////// -------- MC pp anchored PbPB vs pp -------- ////////
// TString* texCollisionDataInfo = new TString("pp #sqrt{#it{s}} = 13.6 TeV");
// const TString* texDatasetsComparisonType = new TString("includePartonEvent");
// const TString* texDatasetsComparisonCommonDenominator = new TString("jet-jet sim");
// const int nDatasets = 2;
// // const TString Datasets[nDatasets] = {"LHC23zzh_apass3_global", "LHC23zzh_apass3_uniform", "LHC23zzh_apass3_itsOnly", "LHC23zzh_apass3_globalWithoutTpcCrossedRows"};
// // const TString Datasets[nDatasets] = {"anchoredJetJetPbPb", "anchoredJetJetpp", "LHC24f3_train240962", "LHC24f3_local", "unanchoredJetJet_train230486"};
// const TString Datasets[nDatasets] = {"anchoredJetJetPbPb", "anchoredJetJetPbPb_includePartonEvent"};
// const TString DatasetsNames[nDatasets] = {"includePartonEvent off", "includePartonEvent on"};
// TFile* file_O2Analysis_list[nDatasets] = {new TFile("Datasets/"+Datasets[0]+"/AnalysisResults.root")
//                                         new TFile("Datasets/"+Datasets[1]+"/AnalysisResults.root")
//                                         };
// const TString analysisWorkflow[nDatasets] = {"jet-finder-charged-qa",
//                                             "jet-finder-charged-qa"
//                                           };

// const TString trainId = "";
// const bool isDatasetWeighted[nDatasets] = {true, true};

// // jaime testJetAnchoredGapGen-20240801-155530


// //////// -------- PbPb MC anchored to LHC23zzh apass4 -------- ////////
// TString* texCollisionDataInfo = new TString("PbPb MC #sqrt{#it{s}} = 5.36 TeV");
// const TString* texDatasetsComparisonType = new TString("centrality");
// const TString* texDatasetsComparisonCommonDenominator = new TString("LHC24d2 - anchored to apass3");
// const int nDatasets = 4;
// // const TString Datasets[nDatasets] = {"LHC23zzh_apass3_global", "LHC23zzh_apass3_uniform", "LHC23zzh_apass3_itsOnly", "LHC23zzh_apass3_globalWithoutTpcCrossedRows"};
// // const TString Datasets[nDatasets] = {"anchoredJetJetPbPb", "anchoredJetJetpp", "LHC24f3_train240962", "LHC24f3_local", "unanchoredJetJet_train230486"};
// const TString Datasets[nDatasets] = {"LHC24d2_PbPbMCapass4Anchored_train258203", "LHC24d2_PbPbMCapass4Anchored_train258203", "LHC24d2_PbPbMCapass4Anchored_train258203", "LHC24d2_PbPbMCapass4Anchored_train258203"};
// const TString DatasetsNames[nDatasets] = {"0-10%", "10-30%", "30-50%", "50-70%"};
// TFile* file_O2Analysis_list[nDatasets] = {new TFile("Datasets/"+Datasets[0]+"/AnalysisResults.root"),
//                                         new TFile("Datasets/"+Datasets[1]+"/AnalysisResults.root"),
//                                         new TFile("Datasets/"+Datasets[2]+"/AnalysisResults.root"),
//                                         new TFile("Datasets/"+Datasets[3]+"/AnalysisResults.root")
//                                         };
// const TString analysisWorkflow[nDatasets] = {"track-efficiency_central0010",
//                                             "track-efficiency_central1030",
//                                             "track-efficiency_central3050",
//                                             "track-efficiency_central5070"
//                                           };

// const TString trainId = "";
// const bool isDatasetWeighted[nDatasets] = {false, false, false, false};

// // jaime testJetAnchoredGapGen-20240801-155530



// //////// -------- accept split vs not accepting split collisions - PbPb MC anchored to LHC23zzh apass4 -------- ////////
// TString* texCollisionDataInfo = new TString("PbPb MC #sqrt{#it{s}} = 5.36 TeV");
// const TString* texDatasetsComparisonType = new TString("splitCollisions");
// const TString* texDatasetsComparisonCommonDenominator = new TString("noSplit");
// const int nDatasets = 8;
// // const TString Datasets[nDatasets] = {"LHC23zzh_apass3_global", "LHC23zzh_apass3_uniform", "LHC23zzh_apass3_itsOnly", "LHC23zzh_apass3_globalWithoutTpcCrossedRows"};
// // const TString Datasets[nDatasets] = {"anchoredJetJetPbPb", "anchoredJetJetpp", "LHC24f3_train240962", "LHC24f3_local", "unanchoredJetJet_train230486"};
// const TString Datasets[nDatasets] = {"LHC24g3", "LHC24g3", "LHC24g3", "LHC24g3", "LHC24g3", "LHC24g3", "LHC24g3", "LHC24g3"};
// const TString DatasetsNames[nDatasets] = {"noSplit 0-10%", "noSplit 10-30%", "noSplit 30-50%", "noSplit 50-70%", "okSplit 0-10%", "okSplit 10-30%", "okSplit 30-50%", "okSplit 50-70%"};
// TFile* file_O2Analysis_list[nDatasets] = {new TFile("Datasets/"+Datasets[0]+"/AnalysisResults.root"),
//                                         new TFile("Datasets/"+Datasets[1]+"/AnalysisResults.root"),
//                                         new TFile("Datasets/"+Datasets[2]+"/AnalysisResults.root"),
//                                         new TFile("Datasets/"+Datasets[3]+"/AnalysisResults.root"),
//                                         new TFile("Datasets/"+Datasets[4]+"/AnalysisResults.root"),
//                                         new TFile("Datasets/"+Datasets[5]+"/AnalysisResults.root"),
//                                         new TFile("Datasets/"+Datasets[6]+"/AnalysisResults.root"),
//                                         new TFile("Datasets/"+Datasets[7]+"/AnalysisResults.root")
//                                         };
// const TString analysisWorkflow[nDatasets] = {"track-efficiency_central0010_id16386",
//                                             "track-efficiency_central1030_id16386",
//                                             "track-efficiency_central3050_id16840",
//                                             "track-efficiency_central5070_id16840",
//                                             "track-efficiency_central0010_id16837",
//                                             "track-efficiency_central1030_id16837",
//                                             "track-efficiency_central3050_id16841",
//                                             "track-efficiency_central5070_id16841",
//                                           };

// const TString trainId = "";
// const bool isDatasetWeighted[nDatasets] = {false, false, false, false, false, false, false, false};

// // jaime testJetAnchoredGapGen-20240801-155530


//////// -------- accept split vs not accepting split collisions - PbPb MC anchored to LHC23zzh apass4 -------- ////////
// TString* texCollisionDataInfo = new TString("PbPb MC #sqrt{#it{s}} = 5.36 TeV");
// const TString* texDatasetsComparisonType = new TString("splitCollisions");
// const TString* texDatasetsComparisonCommonDenominator = new TString("noSplit");
// const int nDatasets = 1;
// // const TString Datasets[nDatasets] = {"LHC23zzh_apass3_global", "LHC23zzh_apass3_uniform", "LHC23zzh_apass3_itsOnly", "LHC23zzh_apass3_globalWithoutTpcCrossedRows"};
// // const TString Datasets[nDatasets] = {"anchoredJetJetPbPb", "anchoredJetJetpp", "LHC24f3_train240962", "LHC24f3_local", "unanchoredJetJet_train230486"};
// const TString Datasets[nDatasets] = {"LHC24g3"}; //DatasetFiles
// const TString DatasetsNames[nDatasets] = {"noSplit 0-10%"};
// TFile* file_O2Analysis_list[nDatasets] = {new TFile(Datasets[0]+"/AnalysisResults.root")};
// const TString analysisWorkflow[nDatasets] = {"track-efficiency_central0010_id16386"};

// const TString trainId = "";
// const bool isDatasetWeighted[nDatasets] = {false};
// const std::string histDatasetComparisonStructure = "twoByTwoDatasetPairs";

// //////// -------- local pp test sim anchored to PbPb comparison to pp anchor -------- ////////
// TString* texCollisionDataInfo = new TString("PbPb MC #sqrt{#it{s}} = 5.36 TeV");
// const TString* texDatasetsComparisonType = new TString("ppAnchorPbPb");
// const TString* texDatasetsComparisonCommonDenominator = new TString("sim");
// const int nDatasets = 1;
// // const TString Datasets[nDatasets] = {"LHC23zzh_apass3_global", "LHC23zzh_apass3_uniform", "LHC23zzh_apass3_itsOnly", "LHC23zzh_apass3_globalWithoutTpcCrossedRows"};
// // const TString Datasets[nDatasets] = {"anchoredJetJetPbPb", "anchoredJetJetpp", "LHC24f3_train240962", "LHC24f3_local", "unanchoredJetJet_train230486"};
// const TString Datasets[nDatasets] = {"pp_sim_anchored_to_PbPb_5360GeV"}; //DatasetFiles
// const TString DatasetsNames[nDatasets] = {"PbPb anchor (5.36TeV)"};
// TFile* file_O2Analysis_list[nDatasets] = {new TFile("Datasets/"+Datasets[0]+"/AnalysisResults.root")
//                                           };
// const TString analysisWorkflow[nDatasets] = {"track-efficiency"
//                                             };

// //////// -------- accept split vs not accepting split collisions - PbPb MC anchored to LHC23zzh apass4 -------- ////////
// TString* texCollisionDataInfo = new TString("Pb-Pb #sqrt{#it{s}} = 5.36 TeV");
// const TString* texDatasetsComparisonType = new TString("Occupancy analysis");
// const TString* texDatasetsComparisonCommonDenominator = new TString("LHC24g3 - sel8FullPbPb");
// const int nDatasets = 8;
// // const TString Datasets[nDatasets] = {"LHC24g3","LHC24d2"};
// // const TString Datasets[nDatasets] = {"LHC24g3"};
//  const TString Datasets[nDatasets] = {"LHC24g3","LHC24d2","LHC24g3","LHC24d2","LHC24g3","LHC24d2","LHC24g3","LHC24d2"};
// const TString DatasetsNames[nDatasets] = {"0-10% pass4","0-10% pass3","10-30% pass4","10-30% pass3","30-50% pass4","30-50% pass3","50-70% pass4","50-70% pass3"};
// TFile* file_O2Analysis_list[nDatasets] = {new TFile(Datasets[0]+"/AnalysisResults_LHC24g3_centrality.root"),
//                                         new TFile(Datasets[1]+"/AnalysisResults_LHC24d2_centrality.root"),
//                                         new TFile(Datasets[2]+"/AnalysisResults_LHC24g3_centrality.root"),
//                                         new TFile(Datasets[3]+"/AnalysisResults_LHC24d2_centrality.root"),
//                                         new TFile(Datasets[4]+"/AnalysisResults_LHC24g3_centrality.root"),
//                                         new TFile(Datasets[5]+"/AnalysisResults_LHC24d2_centrality.root"),
//                                         new TFile(Datasets[4]+"/AnalysisResults_LHC24g3_centrality.root"),
//                                         new TFile(Datasets[5]+"/AnalysisResults_LHC24d2_centrality.root")
//                                         };
// // const TString analysisWorkflow[nDatasets] = {"track-efficiency","track-efficiency"};
// const TString analysisWorkflow[nDatasets] = {"track-efficiency_central0010","track-efficiency_central0010","track-efficiency_central1030","track-efficiency_central1030","track-efficiency_central3050","track-efficiency_central3050","track-efficiency_central5070","track-efficiency_central5070"};

// const TString trainId = "";
// const bool isDatasetWeighted[nDatasets] = {false,false,false,false,false,false,false,false};
// const std::string histDatasetComparisonStructure = "twoByTwoDatasetPairs";



// //////// -------- pp sim test Anchor Pb-Pb gap comparison-------- ////////
// TString* texCollisionDataInfo = new TString("pp #sqrt{#it{s}} = 5.36 TeV");
// const TString* texDatasetsComparisonType = new TString("gap value");
// const TString* texDatasetsComparisonCommonDenominator = new TString("MC anchored to Pb-Pb");
// const int nDatasets = 4;
// const TString Datasets[nDatasets] = {"pp_sim_anchored_to_PbPb_5360GeV_gap2", "pp_sim_anchored_to_PbPb_5360GeV_gap3", "pp_sim_anchored_to_PbPb_5360GeV_gap4", "pp_sim_anchored_to_PbPb_5360GeV_gap5"};
// // const TString DatasetsNames[nDatasets] = {"0-10%", "50-90%"};
// const TString DatasetsNames[nDatasets] = {"gap2", "gap3", "gap4", "gap5"};
// TFile* file_O2Analysis_list[nDatasets] = {new TFile("Datasets/"+Datasets[0]+"/AnalysisResults.root"),
//                                         new TFile("Datasets/"+Datasets[1]+"/AnalysisResults.root"),
//                                         new TFile("Datasets/"+Datasets[2]+"/AnalysisResults.root"),
//                                         new TFile("Datasets/"+Datasets[3]+"/AnalysisResults.root")
//                                         };
// const TString analysisWorkflow[nDatasets] = {"track-efficiency",
//                                           "track-efficiency",
//                                           "track-efficiency",
//                                           "track-efficiency"
//                                           };

// const TString trainId = "";
// const bool isDatasetWeighted[nDatasets] = {false, false, false, false};
// const std::string histDatasetComparisonStructure = "";
// const bool datasetsAreSubsetsofId0 = false;


//////// -------- LHC25b6 - pp sim anchored to PbPb 10% ////////
// TString* texCollisionDataInfo = new TString("pp #sqrt{#it{s}} = 5.36 TeV");
// const TString* texDatasetsComparisonType = new TString("");
// const TString* texDatasetsComparisonCommonDenominator = new TString("");
// const int nDatasets = 2;
// const TString Datasets[nDatasets] = {"LHC25b6_pp_sim_PbPbAnchor_10percent_train370115", "LHC25b4a_pp_ref_Tracks_train371362"};
// // const TString DatasetsNames[nDatasets] = {"0-10%", "50-90%"};
// const TString DatasetsNames[nDatasets] = {"jet-jet MC Pb-Pb anchor", "ppRef gen.purp. MC"};
// TFile* file_O2Analysis_list[nDatasets] = {new TFile("Datasets/"+Datasets[0]+"/AnalysisResults.root"),
//                                           new TFile("Datasets/"+Datasets[1]+"/AnalysisResults.root")
//                                         };

// const TString analysisWorkflow[nDatasets] = {"track-efficiency",
//                                              "track-efficiency"
//                                           };

// const TString wagonId[nDatasets] = {"",
//                                     ""
//                                     };
// const bool isDatasetWeighted[nDatasets] = {true, false};
// const std::string histDatasetComparisonStructure = "";
// const bool datasetsAreSubsetsofId0 = false;

/////
// TString* texCollisionDataInfo = new TString("Pb-Pb #sqrt{#it{s}} = 5.36 TeV");
// const TString* texDatasetsComparisonType = new TString("MC");
// const TString* texDatasetsComparisonCommonDenominator = new TString("LHC24g3");
// const int nDatasets = 8;
// const TString Datasets[nDatasets] = {"LHC24g3","LHC24g3","LHC24g3","LHC24g3","LHC24g3","LHC24g3","LHC24g3","LHC24g3"};
// const TString DatasetsNames[nDatasets] = {"00-10% 1k cut ","00-10% no cut","10-30% 1k cut ","10-30% no cut ","30-50% 1k cut","30-50% no cut ","50-70% 1k cut ","50-70% no cut"};
// TFile* file_O2Analysis_list[nDatasets] = {new TFile(Datasets[0]+"/AnalysisResults_LHC24g3_medium_trackefficiency_occupancy1k.root"),
//                                         new TFile(Datasets[1]+"/AnalysisResults_LHC24g3_medium.root"),
//                                         new TFile(Datasets[2]+"/AnalysisResults_LHC24g3_medium_trackefficiency_occupancy1k.root"),
//                                         new TFile(Datasets[3]+"/AnalysisResults_LHC24g3_medium.root"),
//                                         new TFile(Datasets[4]+"/AnalysisResults_LHC24g3_medium_trackefficiency_occupancy1k.root"),
//                                         new TFile(Datasets[5]+"/AnalysisResults_LHC24g3_medium.root"),
//                                         new TFile(Datasets[6]+"/AnalysisResults_LHC24g3_medium_trackefficiency_occupancy1k.root"),
//                                         new TFile(Datasets[7]+"/AnalysisResults_LHC24g3_medium.root"),
//                                         };
// const TString analysisWorkflow[nDatasets] = {"track-efficiency_central0010","track-efficiency_central0010","track-efficiency_central1030","track-efficiency_central1030","track-efficiency_central3050","track-efficiency_central3050","track-efficiency_central5070","track-efficiency_central5070"
                                         
//                                           };

// const TString trainId = "";
// const bool isDatasetWeighted[nDatasets] = {false,false,false,false,false,false,false,false};
// const bool trackHistsObsoleteVersion[nDatasets] = {false,false,false,false,false,false,false,false};
// const std::string histDatasetComparisonStructure = "twoByTwoDatasetPairs";
// const bool datasetsAreSubsetsofId0 = false;




// TString* texCollisionDataInfo = new TString("Pb-Pb #sqrt{#it{s}} = 5.36 TeV");
// const TString* texDatasetsComparisonType = new TString("MC");
// const TString* texDatasetsComparisonCommonDenominator = new TString("LHC24g3_full 1k occupancy cut");
// const int nDatasets = 4;
// const TString Datasets[nDatasets] = {"LHC24g3","LHC24g3","LHC24g3","LHC24g3"};
// const TString DatasetsNames[nDatasets] = {"00-10%","10-30%","30-50%","50-70%"};
// TFile* file_O2Analysis_list[nDatasets] = {new TFile(Datasets[0]+"/AnalysisResults_LHC24g3_fulltrain_1kcut.root"),
//                                         new TFile(Datasets[1]+"/AnalysisResults_LHC24g3_fulltrain_1kcut.root"),
//                                         new TFile(Datasets[2]+"/AnalysisResults_LHC24g3_fulltrain_1kcut.root"),
//                                         new TFile(Datasets[3]+"/AnalysisResults_LHC24g3_fulltrain_1kcut.root")
//                                         };
// const TString analysisWorkflow[nDatasets] = {"track-efficiency_central0010","track-efficiency_central1030","track-efficiency_central3050","track-efficiency_central5070"
                                         
//                                           };

// const TString trainId = "";
// const bool isDatasetWeighted[nDatasets] = {false,false,false,false};
// const bool trackHistsObsoleteVersion[nDatasets] = {false,false,false,false};
// const std::string histDatasetComparisonStructure = "twoByTwoDatasetPairs";
// const bool datasetsAreSubsetsofId0 = false;



TString* texCollisionDataInfo = new TString("Pb-Pb #sqrt{#it{s}} = 5.36 TeV");
const TString* texDatasetsComparisonType = new TString("MC");
const TString* texDatasetsComparisonCommonDenominator = new TString("LHC24g3_full 1k occupancy cut 0-10%");
// const TString* texDatasetsComparisonCommonDenominator = new TString("LHC24g3_full 1k occupancy cut 50-70%");
const int nDatasets = 26;
const TString Datasets[nDatasets] = {"LHC24g3_1kcut","LHC24g3_1kcut","LHC24g3_1kcut","LHC24g3_1kcut","LHC24g3_1kcut","LHC24g3_1kcut","LHC24g3_1kcut","LHC24g3_1kcut","LHC24g3_1kcut","LHC24g3_1kcut","LHC24g3_1kcut","LHC24g3_1kcut","LHC24g3_1kcut","LHC24g3_1kcut","LHC24g3_1kcut","LHC24g3_1kcut","LHC24g3_1kcut","LHC24g3_1kcut","LHC24g3_1kcut","LHC24g3_1kcut","LHC24g3_1kcut","LHC24g3_1kcut","LHC24g3_1kcut","LHC24g3_1kcut","LHC24g3_1kcut","LHC24g3_1kcut"};
const TString DatasetsNames[nDatasets] = {"6,6 kHz - run544124","6,3 kHz - run544013","11,6 kHz - run544123","12,6 kHz - run544392","13,2 kHz - run544477","14,7 kHz - run544391","15,1 kHz - run544492","16,5 kHz - run544122","17,9 kHz - run544390","18 kHz - run544098","19,4 kHz - run544475","19,4 kHz - run544454","22,5 kHz - run544121","23,7 kHz - run544032","24,2 kHz - run544491","25,1 kHz - run544095","26,8 kHz - run544389","27,9 kHz - run544451","28,7 kHz - run544185","29,3 kHz - run544474","29,5 kHz - run544510","30,3 kHz - run544028","32,9 kHz -run544184","38,3 kHz -run544116","39,6 kHz - run544508","43,1 kHz - run544490"};
TFile* file_O2Analysis_list[nDatasets] = {new TFile(Datasets[0]+"/hy_1575327/AnalysisResults.root"),
                                          new TFile(Datasets[1]+"/hy_1575317/AnalysisResults.root"),
                                          new TFile(Datasets[2]+"/hy_1575326/AnalysisResults.root"),
                                          new TFile(Datasets[3]+"/hy_1575335/AnalysisResults.root"),
                                          new TFile(Datasets[4]+"/hy_1575341/AnalysisResults.root"),
                                          new TFile(Datasets[5]+"/hy_1575334/AnalysisResults.root"),
                                          new TFile(Datasets[6]+"/hy_1575344/AnalysisResults.root"),
                                          new TFile(Datasets[7]+"/hy_1575325/AnalysisResults.root"),
                                          new TFile(Datasets[8]+"/hy_1575333/AnalysisResults.root"),
                                          new TFile(Datasets[9]+"/hy_1575322/AnalysisResults.root"),
                                          new TFile(Datasets[10]+"/hy_1575339/AnalysisResults.root"),
                                          new TFile(Datasets[11]+"/hy_1575337/AnalysisResults.root"),
                                          new TFile(Datasets[12]+"/hy_1575324/AnalysisResults.root"),
                                          new TFile(Datasets[13]+"/hy_1575319/AnalysisResults.root"),
                                          new TFile(Datasets[14]+"/hy_1575343/AnalysisResults.root"),
                                          new TFile(Datasets[15]+"/hy_1575321/AnalysisResults.root"),
                                          new TFile(Datasets[16]+"/hy_1575332/AnalysisResults.root"),
                                          new TFile(Datasets[17]+"/hy_1575336/AnalysisResults.root"),
                                          new TFile(Datasets[18]+"/hy_1575331/AnalysisResults.root"),
                                          new TFile(Datasets[19]+"/hy_1575338/AnalysisResults.root"),
                                          new TFile(Datasets[20]+"/hy_1575346/AnalysisResults.root"),
                                          new TFile(Datasets[21]+"/hy_1575318/AnalysisResults.root"),
                                          new TFile(Datasets[22]+"/hy_1575330/AnalysisResults.root"),
                                          new TFile(Datasets[23]+"/hy_1575323/AnalysisResults.root"),
                                          new TFile(Datasets[24]+"/hy_1575345/AnalysisResults.root"),
                                          new TFile(Datasets[25]+"/hy_1575342/AnalysisResults.root")
                                        };
const TString analysisWorkflow[nDatasets] = {"track-efficiency_central0010","track-efficiency_central0010","track-efficiency_central0010","track-efficiency_central0010","track-efficiency_central0010","track-efficiency_central0010","track-efficiency_central0010","track-efficiency_central0010","track-efficiency_central0010","track-efficiency_central0010","track-efficiency_central0010","track-efficiency_central0010","track-efficiency_central0010","track-efficiency_central0010","track-efficiency_central0010","track-efficiency_central0010","track-efficiency_central0010","track-efficiency_central0010","track-efficiency_central0010","track-efficiency_central0010","track-efficiency_central0010","track-efficiency_central0010","track-efficiency_central0010","track-efficiency_central0010","track-efficiency_central0010","track-efficiency_central0010"};
// const TString analysisWorkflow[nDatasets] = {"track-efficiency_central5070","track-efficiency_central5070","track-efficiency_central5070","track-efficiency_central5070","track-efficiency_central5070","track-efficiency_central5070","track-efficiency_central5070","track-efficiency_central5070","track-efficiency_central5070","track-efficiency_central5070","track-efficiency_central5070","track-efficiency_central5070","track-efficiency_central5070","track-efficiency_central5070","track-efficiency_central5070","track-efficiency_central5070","track-efficiency_central5070","track-efficiency_central5070","track-efficiency_central5070","track-efficiency_central5070","track-efficiency_central5070","track-efficiency_central5070","track-efficiency_central5070","track-efficiency_central5070","track-efficiency_central5070","track-efficiency_central5070"};

const TString trainId = "";
const bool isDatasetWeighted[nDatasets] = {false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false};
const bool trackHistsObsoleteVersion[nDatasets] = {false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false};
const std::string histDatasetComparisonStructure = "";
const bool datasetsAreSubsetsofId0 = false;


// TString* texCollisionDataInfo = new TString("pp anchored to PbPb 2023 apass4");
// const TString* texDatasetsComparisonType = new TString("MC");
// const TString* texDatasetsComparisonCommonDenominator = new TString("LHC25b6 no occupancy cut");
// const int nDatasets = 15;
// const TString Datasets[nDatasets] = {"LHC25b6_nocut","LHC25b6_nocut","LHC25b6_nocut","LHC25b6_nocut","LHC25b6_nocut","LHC25b6_nocut","LHC25b6_nocut","LHC25b6_nocut","LHC25b6_nocut","LHC25b6_nocut","LHC25b6_nocut","LHC25b6_nocut","LHC25b6_nocut","LHC25b6_nocut","LHC25b6_nocut"};
// const TString DatasetsNames[nDatasets] = {"6,3 kHz - run544013","6,6 kHz - run544124","13,2 kHz - run544477","15,1 kHz - run544492","16,1 kHz - run544476","19,4 kHz - run544475","19,4 kHz - run544454","23,7 kHz - run544032","24,2 kHz - run544491","25,1 kHz - run544095","28,7 kHz - run544185","29,32 kHz - run544091","30,3 kHz - run544028","38,3 kHz -run544116","43,1 kHz - run544490"};
// TFile* file_O2Analysis_list[nDatasets] = {new TFile(Datasets[0]+"/hy_1614812/AnalysisResults.root"),
//                                           new TFile(Datasets[1]+"/hy_1614822/AnalysisResults.root"),
//                                           new TFile(Datasets[2]+"/hy_1614834/AnalysisResults.root"),           
//                                           new TFile(Datasets[3]+"/hy_1614837/AnalysisResults.root"),
//                                           new TFile(Datasets[4]+"/hy_1614833/AnalysisResults.root"),
//                                           new TFile(Datasets[5]+"/hy_1614832/AnalysisResults.root"),
//                                           new TFile(Datasets[6]+"/hy_1614830/AnalysisResults.root"),                                   
//                                           new TFile(Datasets[7]+"/hy_1614814/AnalysisResults.root"),
//                                           new TFile(Datasets[8]+"/hy_1614836/AnalysisResults.root"),
//                                           new TFile(Datasets[9]+"/hy_1614816/AnalysisResults.root"),                                                                                
//                                           new TFile(Datasets[10]+"/hy_1614824/AnalysisResults.root"),
//                                           new TFile(Datasets[11]+"/hy_1614815/AnalysisResults.root"),                                                                           
//                                           new TFile(Datasets[12]+"/hy_1614813/AnalysisResults.root"),                                  
//                                           new TFile(Datasets[13]+"/hy_1614818/AnalysisResults.root"),                                   
//                                           new TFile(Datasets[14]+"/hy_1614835/AnalysisResults.root")
//                                         };
// const TString analysisWorkflow[nDatasets] = {"track-efficiency","track-efficiency","track-efficiency","track-efficiency","track-efficiency","track-efficiency","track-efficiency","track-efficiency","track-efficiency","track-efficiency","track-efficiency","track-efficiency","track-efficiency","track-efficiency","track-efficiency"
                                         
//                                           };

// const TString trainId = "";
// const bool isDatasetWeighted[nDatasets] = {false,false,false,false,false,false,false,false,false,false,false,false,false,false,false};
// const bool trackHistsObsoleteVersion[nDatasets] = {false,false,false,false,false,false,false,false,false,false,false,false,false,false,false};
// const std::string histDatasetComparisonStructure = "";
// const bool datasetsAreSubsetsofId0 = false;



// TString* texCollisionDataInfo = new TString("Pb-Pb #sqrt{#it{s}} = 5.36 TeV");
// const TString* texDatasetsComparisonType = new TString("MC");
// const TString* texDatasetsComparisonCommonDenominator = new TString("LHC24g3_full 1k occupancy cut 0-10%");
// // const TString* texDatasetsComparisonCommonDenominator = new TString("LHC24g3_full 1k occupancy cut 50-70%");
// const int nDatasets = 2;
// const TString Datasets[nDatasets] = {"LHC24g3_1kcut","LHC24g3_1kcut"};
// const TString DatasetsNames[nDatasets] = {"6,3 kHz - run544013_louise","6,3 kHz - run544013_aimeric"};
// TFile* file_O2Analysis_list[nDatasets] = {new TFile(Datasets[0]+"/hy_1575317/AnalysisResults.root"),                                   
//                                           new TFile(Datasets[1]+"/AnalysisResults_run544013.root")
//                                         };
// const TString analysisWorkflow[nDatasets] = {"track-efficiency_central0010","track-efficiency_central0010"};
// // const TString analysisWorkflow[nDatasets] = {"track-efficiency_central5070","track-efficiency_central5070","track-efficiency_central5070","track-efficiency_central5070","track-efficiency_central5070","track-efficiency_central5070","track-efficiency_central5070","track-efficiency_central5070","track-efficiency_central5070","track-efficiency_central5070","track-efficiency_central5070","track-efficiency_central5070","track-efficiency_central5070","track-efficiency_central5070","track-efficiency_central5070","track-efficiency_central5070","track-efficiency_central5070","track-efficiency_central5070","track-efficiency_central5070","track-efficiency_central5070","track-efficiency_central5070","track-efficiency_central5070","track-efficiency_central5070","track-efficiency_central5070","track-efficiency_central5070","track-efficiency_central5070"};


// const TString trainId = "";
// const bool isDatasetWeighted[nDatasets] = {false,false};
// const bool trackHistsObsoleteVersion[nDatasets] = {false,false};
// const std::string histDatasetComparisonStructure = "";
// const bool datasetsAreSubsetsofId0 = false;





