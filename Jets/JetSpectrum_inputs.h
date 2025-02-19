
// To avoid your local _inputs.h files to be replaced by the template in the git repository, you should ask git to ignore them with git update-index --assume-unchanged (for more details see https://stackoverflow.com/questions/3319479/can-i-git-commit-a-file-and-ignore-its-content-changes)

#ifndef JETSPECTRUM_INPUTS_H
#define JETSPECTRUM_INPUTS_H

//////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////       file access choice       ////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////

// //////// -------- Pb-Pb -------- ////////
// TString* texCollisionDataInfo = new TString("Pb-Pb #sqrt{#it{s}} = 5.36 TeV"); 
// const TString* texDatasetsComparisonType = new TString("");
// const TString* texDatasetsComparisonCommonDenominator = new TString("LHC23zzh");
// const int nDatasets = 1;
// const TString Datasets[nDatasets] = {"LHC23zzh_apass4_occupancy01000_train297793"};
// const TString DatasetsNames[nDatasets] = {""};
// TFile* file_O2Analysis_list[nDatasets] = {new TFile("Datasets/"+Datasets[0]+"/AnalysisResults.root")
//                                       };
// // TFile* file_O2Analysis_MCfileForMatrix[nDatasets] = new TFile("Datasets/ppSim_LHC23d4/AnalysisResults.root");
// // TFile* file_O2Analysis_MCfileForMatrix = new TFile("Datasets/ppSim_LHC23d4_weighted_withLeadingTrackCut/AnalysisResults.root");
// TFile* file_O2Analysis_MCfileForMatrix = new TFile("Datasets/ppSim_LHC23d4_weighted_train_256548/AnalysisResults.root");
// TFile* file_O2Analysis_ppSimDetectorEffect_unfoldingControl = {new TFile("Datasets/LHC24b1b_sel8MC_train239181/OneRun/AnalysisResults.root")};

// // const TString trainId = "_id12832";
// // const TString analysisWorkflowData = "jet-finder-charged-qa_central_5090_lead5"+trainId;
// // const TString trainId = "_id12436";
// // const TString analysisWorkflowData = "jet-finder-charged-qa_central_0010_lead5"+trainId;
// const TString trainId = "";
// const TString analysisWorkflowData = "jet-finder-charged-qa"+trainId;

// const TString analysisWorkflowMC = "jet-finder-charged-qa";





// //////// -------- pp -------- ////////
// TString* texCollisionDataInfo = new TString("pp #sqrt{#it{s}} = 13.6 TeV"); 
// const TString* texDatasetsComparisonType = new TString("");
// const TString* texDatasetsComparisonCommonDenominator = new TString("LHC22o pass7");
// const int nDatasets = 1;
// const TString Datasets[nDatasets] = {"LHC22o_pass7_train261733"};
// const TString DatasetsNames[nDatasets] = {"LHC22o_pass7"};
// TFile* file_O2Analysis_list[nDatasets] = {new TFile("Datasets/"+Datasets[0]+"/AnalysisResults.root")
//                                       };
// TFile* file_O2Analysis_MCfileForMatrix = new TFile("Datasets/LHC24f3b_train261768/half1/AnalysisResults.root");

// // TFile* file_O2Analysis_ppSimDetectorEffect_unfoldingControl = {new TFile("Datasets/LHC24f3_sel8MC_train240962/AnalysisResults.root")};
// TFile* file_O2Analysis_ppSimDetectorEffect_unfoldingControl = {new TFile("Datasets/LHC24f3b_train261768/half2/AnalysisResults.root")};

// // TFile* file_O2Analysis_MCfileForMatrix = new TFile("Datasets/LHC24b1b_sel8Full_train239409/AnalysisResults.root");
// // TFile* file_O2Analysis_MCfileForMatrix = new TFile("Datasets/ppSim_LHC23d4_weighted/AnalysisResults.root");

// const TString trainId = "";
// // const TString trainId = "_id12832";
// // const TString analysisWorkflowData = "jet-finder-charged-qa_central_5090_lead5"+trainId;
// const TString analysisWorkflowData = "jet-finder-charged-qa"+trainId;

// const TString analysisWorkflowMC = "jet-finder-charged-qa";
// // const TString analysisWorkflowMC = "jet-finder-charged-qa_global_CollMatch";





// //////// -------- pp spectrum with Joonsuk files tests -------- ////////
// TString* texCollisionDataInfo = new TString("pp #sqrt{#it{s}} = 13.6 TeV"); 
// const TString* texDatasetsComparisonType = new TString("");
// const TString* texDatasetsComparisonCommonDenominator = new TString("");
// const int nDatasets = 1;
// const TString Datasets[nDatasets] = {"Data_LHC22o_minBias_small_Joonsuk"};
// const TString DatasetsNames[nDatasets] = {""};
// TFile* file_O2Analysis_list[nDatasets] = {new TFile("Datasets/"+Datasets[0]+"/AnalysisResults.root")
//                                       };
// // TFile* file_O2Analysis_MCfileForMatrix[nDatasets] = new TFile("Datasets/ppSim_LHC23d4/AnalysisResults.root");
// // TFile* file_O2Analysis_MCfileForMatrix = new TFile("Datasets/ppSim_LHC23d4_weighted_withLeadingTrackCut/AnalysisResults.root");
// TFile* file_O2Analysis_MCfileForMatrix = new TFile("Datasets/ppSim_LHC24f3b_withTrackTuner_Joonsuk/AnalysisResults.root");
// TFile* file_O2Analysis_ppSimDetectorEffect_unfoldingControl = {new TFile("Datasets/ppSim_LHC24f3b_withTrackTuner_Joonsuk/AnalysisResults.root")};

// // const TString trainId = "_id12832";
// // const TString analysisWorkflowData = "jet-finder-charged-qa_central_5090_lead5"+trainId;
// // const TString trainId = "_id12436";
// // const TString analysisWorkflowData = "jet-finder-charged-qa_central_0010_lead5"+trainId;
// const TString trainId = "";
// const TString analysisWorkflowData = "jet-finder-charged-qa"+trainId;
// const TString analysisWorkflow_unfoldingControl = "jet-finder-charged-qa"+trainId;

// const TString analysisWorkflowMC = "jet-finder-charged-qa";







//////// -------- Pb-Pb spectrum with Wenhui Angantyr files tests -------- ////////
TString* texCollisionDataInfo = new TString("Pb-Pb #sqrt{#it{s}} = 5.36 TeV"); 
const TString* texDatasetsComparisonType = new TString("");
const TString* texDatasetsComparisonCommonDenominator = new TString("");
const int nDatasets = 1;
const TString Datasets[nDatasets] = {"Data_halfMCAngantyr"};
const TString DatasetsNames[nDatasets] = {""};
TFile* file_O2Analysis_list[nDatasets] = {new TFile("Datasets/"+Datasets[0]+"/AnalysisResults.root")
                                      };
// TFile* file_O2Analysis_MCfileForMatrix[nDatasets] = new TFile("Datasets/ppSim_LHC23d4/AnalysisResults.root");
// TFile* file_O2Analysis_MCfileForMatrix = new TFile("Datasets/ppSim_LHC23d4_weighted_withLeadingTrackCut/AnalysisResults.root");
// TFile* file_O2Analysis_MCfileForMatrix = new TFile("Datasets/MC_halfMCAngantyr/AnalysisResults.root");
TFile* file_O2Analysis_MCfileForMatrix = new TFile("Datasets/MC_halfMCAngantyr/AnalysisResults.root");
TFile* file_O2Analysis_ppSimDetectorEffect_unfoldingControl = {new TFile("Datasets/MC_halfMCAngantyr/AnalysisResults.root")};

// const TString trainId = "_id12832";
// const TString analysisWorkflowData = "jet-finder-charged-qa_central_5090_lead5"+trainId;
// const TString trainId = "_id12436";
// const TString analysisWorkflowData = "jet-finder-charged-qa_central_0010_lead5"+trainId;
const TString trainId = "";
const TString analysisWorkflowData = "jet-spectra-charged"+trainId;
const TString analysisWorkflowBkg = "jet-background-analysis"+trainId;
const TString analysisWorkflow_unfoldingControl = "jet-spectra-charged"+trainId;

const TString analysisWorkflowMC = "jet-spectra-charged";


changelist: 
- adding analysisWorkflowBkg variable for fluctuation response
- divided isPbPb bool into isPbPb and useFactorisedMatrix
- replace file_O2Analysis_ppSimDetectorEffect name by file_O2Analysis_MCfileForMatrix as it can be used for Pb-Pb with gen purpose mc
- trying to allow both jet-finder-charged-qa and jet-spectra-charged as inputs despite different histograms ...


TODO: 
- fcontrolMC only works for pp for now
- keep changing Get_Pt_spectrum_mcpMatched_preWidthScalingAndEvtNorm and other matching functions with correct hist names, try adding compatibility with wenhui (here I can understand the change of wording as base tag aint great compared to reco/gen) following structure:
    "if (!isPbPb || useFactorisedMatrix == true){
      if (!fcontrolMC){"
#endif
