
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
const TString Datasets[nDatasets] = {"Data_halfMCAngantyr_train356485"};
const TString DatasetsNames[nDatasets] = {""};
TFile* file_O2Analysis_list[nDatasets] = {new TFile("Datasets/"+Datasets[0]+"/AnalysisResults.root")
                                      };
// TFile* file_O2Analysis_MCfileForMatrix[nDatasets] = new TFile("Datasets/ppSim_LHC23d4/AnalysisResults.root");
// TFile* file_O2Analysis_MCfileForMatrix = new TFile("Datasets/ppSim_LHC23d4_weighted_withLeadingTrackCut/AnalysisResults.root");
// TFile* file_O2Analysis_MCfileForMatrix = new TFile("Datasets/MC_halfMCAngantyr/AnalysisResults.root");
TFile* file_O2Analysis_MCfileForMatrix = new TFile("Datasets/MC_halfMCAngantyr_train356485/AnalysisResults.root");
TFile* file_O2Analysis_ppSimDetectorEffect_unfoldingControl = {new TFile("Datasets/MC_halfMCAngantyr_train356485/AnalysisResults.root")}; // use this MC file as input to unfolding (with h_jet_pt_rhoareasubtracted distrib on file) and as comparison to gen (with h_jet_pt_part distrib on file)

// const TString trainId = "_id12832";
// const TString analysisWorkflowData = "jet-finder-charged-qa_central_5090_lead5"+trainId;
// const TString trainId = "_id12436";
// const TString analysisWorkflowData = "jet-finder-charged-qa_central_0010_lead5"+trainId;
const TString trainIdData = "_id24998";
const TString analysisWorkflowData = "jet-spectra-charged"+trainIdData;
const TString trainIdBkg = "_id24998";
const TString analysisWorkflowBkg = "jet-background-analysis"+trainIdBkg;
const TString trainIdUnfoldingControl = "_id24998";
const TString analysisWorkflow_unfoldingControl = "jet-spectra-charged"+trainIdUnfoldingControl;

const TString trainIdMC = "_id24998";
const TString analysisWorkflowMC = "jet-spectra-charged"+trainIdMC;


// changelist: 
// - adding analysisWorkflowBkg variable for fluctuation response
// - changed isDataPbPb bool to useFactorisedMatrix, and reworked the isDataPbPb, ppMcIsWeighted (now mcIsWeighted), useFactorisedMatrix to be used only when necessary
// - replace file_O2Analysis_ppSimDetectorEffect name by file_O2Analysis_MCfileForMatrix as it can be used for Pb-Pb with gen purpose mc
// - adding compatibility with jet-spectra-charged workflow
// - removed all centralty mention; now we assume the centrality has been selected before making the AO2D
// - normunfoldedByNEvts was set to !normunfoldedByNEvts in tests, instead of normunfoldedByNEvts
// - renamed doEvtNorm to normaliseDistribsBeforeUnfolding, and normUnfoldedByNEvts to normaliseDistribsAfterUnfolding
// - removed option "evtNorm" check as it was doing the same as normaliseDistribsBeforeUnfolding

// TODO: 
// - normUnfoldedByNEvts vs doEvtNorm vs options.find("evtNorm") ?????
//     - normUnfoldedByNEvts: 
//         - appears three times, in Get_Pt_spectrum_unfolded_preWidthScaling, Get_Pt_spectrum_unfoldedThenRefolded_RooUnfoldMethod_preWidthScaling and Get_Pt_spectrum_unfoldedThenRefolded_preWidthScaling
//         - if false, no normalisation is done at all (not even entries norm)
//         -> probably should be renamed normaliseDistribsAfterUnfolding
//     - doEvtNorm: 
//         - appears where normUnfoldedByNEvts does (get_pt_spectrum_unfolded.....)
//         - appears in: Draw_Pt_spectrum_raw, Draw_Pt_spectrum_mcp, Draw_Pt_spectrum_mcdMatched and Draw_Pt_spectrum_unfolded, only for LABELS
//         - appears in Get_Pt_spectrum_bkgCorrected_rec/gen/fineBinning_preWidthScaling, Get_Pt_spectrum_mcp_gen/fineBinning_preWidthScaling, Get_Pt_spectrum_mcd_rec/fine/genBinning_preWidthScaling, Get_Pt_spectrum_mcdMatched_gen/rec/fineBinning_preWidthScaling, Get_Pt_spectrum_mcpMatched_gen/fine/recBinning_preWidthScaling
//         -> probably should be renamed normaliseDistribsBeforeUnfolding, and only affect the Get_Pt_spectrum_ nonunfolded functions
//       and draw raw/mcp/mcdMatched and unfolded labels should be set with something like normaliseDistribsBeforeUnfolding || normaliseDistribsAfterUnfolding (although what happens if both are normalised?)

//       "evtNorm" option is always used with doEvtNorm -> probably should just be removed; also entriesnorm never used


#endif
