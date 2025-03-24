#ifndef JETSPECTRUM_UNFOLDING_C
#define JETSPECTRUM_UNFOLDING_C

#include "JetSpectrum_Unfolding.h"

#include "./JetSpectrum_ResponseMatrixFunctions.h"
#include "./JetSpectrum_ResponseMatrixFunctions.C"
#include "./JetSpectrum_SpectraGetters.h"
#include "./JetSpectrum_SpectraGetters.C"
#include "./JetSpectrum_EfficiencyPurityGetters.h"
#include "./JetSpectrum_EfficiencyPurityGetters.C"

#include "../Settings/AxisTitles.h"
#include "../Settings/GlobalSettings.h"
#include "../Utilities/AnalysisUtilities.h"
#include "../Utilities/HistogramUtilities.h"
#include "../Utilities/HistogramPlotting.h"
#include "../Utilities/AnalysisUtilities.C" 
#include "../Utilities/HistogramUtilities.C"
#include "../Utilities/HistogramPlotting.C" 

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////// RooUnfold Custom Utilities ///////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



int GetSvdBestRegularisationParameter_notYetSatisfying(TSVDUnfold* unfoldTSvd){
  //alternative ideas: average over 3-4 k values, find when ratio of those as k increases becomes 1ish
  TH1D* H1D_D = unfoldTSvd->GetD();
  double minimum = 999999;
  double tempContent;
  int k = -99;
  for (int iBinX = 1; iBinX < H1D_D->GetNbinsX(); iBinX++) {
    tempContent = abs(abs(H1D_D->GetBinContent(iBinX)) - 1); // one wants |k| closest to 1 as possible
    if (minimum > tempContent) {
      minimum = tempContent;
      k = iBinX - 1; // or is it iBinX?
    }
  }
  cout << "k = " << k << ", |H1D-1| minimum = " << minimum << endl;
  return k;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////// Spectrum Unfolding functions ////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


std::pair<int, RooUnfold*> Get_Pt_spectrum_unfolded_preWidthScalingAndEvtNorm(TH1D* &H1D_jetPt_unfolded, int iDataset, int iRadius, int unfoldParameterInput, std::string options) {
  //for now makes the assumption gen and rec have the same pT binning

  TString partialUniqueSpecifier = Datasets[iDataset]+"_R="+Form("%.1f",arrayRadius[iRadius]);

  TH1D* measured;
  TH1D* mcp;
  TH1D* mcd;

  if (!normGenAndMeasByNEvts) {
    Get_Pt_spectrum_bkgCorrected_recBinning_preWidthScalingAndEvtNorm(measured, iDataset, iRadius, options);
    Get_Pt_spectrum_mcp_genBinning_preWidthScalingAndEvtNorm(mcp, iDataset, iRadius, false, options);
    Get_Pt_spectrum_mcd_recBinning_preWidthScalingAndEvtNorm(mcd, iDataset, iRadius, options);

    if (useFineBinningTest) {
      Get_Pt_spectrum_bkgCorrected_fineBinning_preWidthScalingAndEvtNorm(measured, iDataset, iRadius, options);
      Get_Pt_spectrum_mcp_fineBinning_preWidthScalingAndEvtNorm(mcp, iDataset, iRadius, false, options);
      Get_Pt_spectrum_mcd_fineBinning_preWidthScalingAndEvtNorm(mcd, iDataset, iRadius, options);
    }
  } else{
    Get_Pt_spectrum_bkgCorrected_recBinning_preWidthScaling(measured, iDataset, iRadius, options);
    Get_Pt_spectrum_mcp_genBinning_preWidthScaling(mcp, iDataset, iRadius, false, options);
    Get_Pt_spectrum_mcd_recBinning_preWidthScaling(mcd, iDataset, iRadius, options);

    if (useFineBinningTest) {
      Get_Pt_spectrum_bkgCorrected_fineBinning_preWidthScaling(measured, iDataset, iRadius, options);
      Get_Pt_spectrum_mcp_fineBinning_preWidthScaling(mcp, iDataset, iRadius, false, options);
      Get_Pt_spectrum_mcd_fineBinning_preWidthScaling(mcd, iDataset, iRadius, options);
    }
  } 

  bool divideSuccessFakes;
  TH1D* H1D_jetFakes;
  if (useManualRespMatrixSettingMethod) {
    if (applyFakes) {
      if (!useFineBinningTest){ 
        divideSuccessFakes = Get_Pt_JetFakes(H1D_jetFakes, iDataset, iRadius, options);
      } else {
        divideSuccessFakes = Get_Pt_JetFakes_fineBinning(H1D_jetFakes, iDataset, iRadius, options);
      }
      if (divideSuccessFakes){
        measured->Multiply(H1D_jetFakes);
      } else {
        cout << "################## measured->Multiply(H1D_jetFakes) failed because Get_Pt_JetFakes() FAILED!!!!! ##################" << endl;
      }
    }
  }


  // TH1D* measured;
  // Get_Pt_spectrum_bkgCorrected_recBinning(measured, iDataset, iRadius, options);
  // TH1D* mcp;
  // Get_Pt_spectrum_mcp_genBinning(mcp, iDataset, iRadius, options);
  // TH1D* mcdMatched;
  // Get_Pt_spectrum_mcdMatched_genBinning(mcdMatched, iDataset, iRadius, options);

  TH1D* H1D_kinematicEfficiency;
  TH2D* H2D_jetPtResponseMatrix_fluctuations;
  TH2D* H2D_jetPtResponseMatrix_detectorResponse;
  TH2D* H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined;
  TH2D* H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined_fineBinning;
  
  Get_PtResponseMatrix_Fluctuations(H2D_jetPtResponseMatrix_fluctuations, iDataset, iRadius);
  Get_PtResponseMatrix_detectorResponse(H2D_jetPtResponseMatrix_detectorResponse, iDataset, iRadius);
  // compute matrixFluctuations times matrixDetector

  Get_PtResponseMatrix_DetectorAndFluctuationsCombined_fineBinning(H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined_fineBinning, H2D_jetPtResponseMatrix_detectorResponse, H2D_jetPtResponseMatrix_fluctuations, iDataset, iRadius, options);
  Get_ResponseMatrix_Pt_KinematicEffiency(H1D_kinematicEfficiency, H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined_fineBinning, partialUniqueSpecifier, iRadius); // I want the efficiency before the reweighting and normalisation
  if (useFineBinningTest) {
    ReweightResponseMatrixWithPrior_fineBinning(H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined_fineBinning, iDataset, iRadius, options);
  }
  
  Get_PtResponseMatrix_DetectorAndFluctuationsCombined(H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined, H2D_jetPtResponseMatrix_detectorResponse, H2D_jetPtResponseMatrix_fluctuations, iDataset, iRadius, options);
  // FinaliseResponseMatrix(H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined, iDataset, iRadius, options);

  if (drawIntermediateResponseMatrices) {
    TH2D* H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined_postWeighting = (TH2D*)H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined->Clone("H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined_postWeighting"+partialUniqueSpecifier);

    TString* pdfNamePost = new TString("responseMatrix_combined_postWeighting"+(TString)"_R="+Form("%.1f",arrayRadius[iRadius])+"_"+Datasets[iDataset]+DatasetsNames[iDataset]);
    TString* pdfNamePost_logz = new TString("responseMatrix_combined_postWeighting"+(TString)"_R="+Form("%.1f",arrayRadius[iRadius])+"_"+Datasets[iDataset]+DatasetsNames[iDataset]+"_logz");

    TString textContextPost(contextCustomTwoFields(*texDatasetsComparisonCommonDenominator, contextJetRadius(arrayRadius[iRadius]), ""));

    Draw_TH2_Histogram(H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined_postWeighting, textContextPost, pdfNamePost, texPtJetRecX, texPtJetGenX, texCollisionDataInfo, drawnWindowAuto, th2ContoursNone, contourNumberNone, "");
    Draw_TH2_Histogram(H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined_postWeighting, textContextPost, pdfNamePost_logz, texPtJetRecX, texPtJetGenX, texCollisionDataInfo, drawnWindowAuto, th2ContoursNone, contourNumberNone, "logz");
  }

  cout << "RooUnfoldResponse setting - start" << endl;
  RooUnfoldResponse* response;
  if (useManualRespMatrixSettingMethod) {

    // // based on Marta's work: https://twiki.cern.ch/twiki/bin/viewauth/ALICE/JEJetSpectrumUnfolding

    if (useFineBinningTest) {
      response = new RooUnfoldResponse(0, 0, H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined_fineBinning); // measured and mcp_rebinned are here to take inneficiencies and fakes into account; or is it really what's happening? 'Alternatively, the response matrix can be constructed from a pre-existing TH2D 2-dimensional histogram (with truth and measured distribution TH1D histograms for normalisation).' from https://hepunx.rl.ac.uk/~adye/software/unfold/RooUnfold.html, and I'm already normalising so maybe tehre's no need for more normalisation
    } else {
      response = new RooUnfoldResponse(0, 0, H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined); // measured and mcp_rebinned are here to take inneficiencies and fakes into account; or is it really what's happening? 'Alternatively, the response matrix can be constructed from a pre-existing TH2D 2-dimensional histogram (with truth and measured distribution TH1D histograms for normalisation).' from https://hepunx.rl.ac.uk/~adye/software/unfold/RooUnfold.html, and I'm already normalising so maybe tehre's no need for more normalisation
    }

    cout << "Get_Pt_spectrum_unfolded(): should I use response->UseOverflow() ? using it gives a ratio unfolded/mcp much higher than without using it" << endl;
    // response->UseOverflow();

  } else {
    TH2D *Respt;
    if (useFineBinningTest){
      Respt = H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined_fineBinning;
    } else {
      Respt = H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined;
    }
    TH1D *mcdMatched = (TH1D*)Respt->ProjectionX("Respt_projX", 1, Respt->GetNbinsY());
    TH1D *mcpMatched = (TH1D*)Respt->ProjectionY("Respt_projY", 1, Respt->GetNbinsX());
    TH1D *fake = (TH1D*)mcd->Clone();
    if (normDetRespByNEvts) {
      if (mcIsWeighted) {
        fake->Scale(1./GetNEventsSelected_JetFramework_weighted(file_O2Analysis_MCfileForMatrix, analysisWorkflowMC));
      } else {
        fake->Scale(1./GetNEventsSelected_JetFramework(file_O2Analysis_MCfileForMatrix, analysisWorkflowMC));
      }
    }
    for (auto i = 1; i <= fake->GetNbinsX(); i++) { 
      cout << "Fake i = " << i << " content: "<< fake->GetBinContent(i) << " - " << mcdMatched->GetBinContent(i) << " = " << fake->GetBinContent(i) - mcdMatched->GetBinContent(i) << " ------ error of mcd = " << fake->GetBinError(i) << ", error of mcdMatched = " << mcdMatched->GetBinError(i) << endl;
    }
    fake->Add(mcdMatched, -1);
    TH1D *miss = (TH1D*)mcp->Clone(); // for each bin of gen-level jet distribution, how many are not matched to a rec-level jet
    if (normDetRespByNEvts) {
      if (mcIsWeighted) {
        miss->Scale(1./GetNEventsSelected_JetFramework_weighted(file_O2Analysis_MCfileForMatrix, analysisWorkflowMC));
      } else {
        miss->Scale(1./GetNEventsSelected_JetFramework(file_O2Analysis_MCfileForMatrix, analysisWorkflowMC));
      }
    }
    for (auto i = 1; i <= fake->GetNbinsX(); i++) { 
      cout << "Miss i = " << i << " content: "<< miss->GetBinContent(i) << " - " << mcpMatched->GetBinContent(i) << " = " << miss->GetBinContent(i) - mcpMatched->GetBinContent(i) << " ------ error of mcp = " << miss->GetBinError(i) << ", error of mcpMatched = " << mcpMatched->GetBinError(i) << endl;
    }
    miss->Add(mcpMatched, -1); // for each bin of rec-level jet distribution, how many are not matched to a gen-level jet

    response = new RooUnfoldResponse(mcd, mcp);
    for (auto i = 1; i <= Respt->GetNbinsX(); i++) {
      for (auto j = 1; j <= Respt->GetNbinsY(); j++) { // ptpair
        Double_t bincenx = Respt->GetXaxis()->GetBinCenter(i);
        Double_t binceny = Respt->GetYaxis()->GetBinCenter(j);
        Double_t bincont = Respt->GetBinContent(i, j);
        response->Fill(bincenx, binceny, bincont);
        // cout << "response i,j = " << i << "," << j << " content: "<< bincont << endl;
      }
    }
    for (auto i = 1; i <= miss->GetNbinsX(); i++) { 
      Double_t bincenx = miss->GetXaxis()->GetBinCenter(i);
      Double_t bincont = miss->GetBinContent(i);
      response->Miss(bincenx, bincont);
      // cout << "Miss i = " << i << " content: "<< bincont << endl;
    }
    for (auto i = 1; i <= fake->GetNbinsX(); i++) {
      Double_t bincenx = fake->GetXaxis()->GetBinCenter(i);
      Double_t bincont = fake->GetBinContent(i);
      response->Fake(bincenx, bincont);
      // cout << "Fake i = " << i << " content: "<< bincont << endl;
    }


    // Get_ResponseMatrix_Pt_KinematicEffiency(H1D_kinematicEfficiency, (TH2D*)response->Hresponse(), partialUniqueSpecifier, iRadius); // I want the efficiency before the reweighting and normalisation

  }
  cout << "RooUnfoldResponse setting - end" << endl;



  // note about unfolding: if SUM(Mkj over k)=1 then it's also true for the inverse of M (can easily demosntrate it, just write it);
  RooUnfoldBayes* unfoldBayes = new RooUnfoldBayes(response, measured, unfoldParameterInput);
  int unfoldParameterSvdInitial = 1;
  RooUnfoldSvd* unfoldSvdInitialiser = new RooUnfoldSvd(response, measured, unfoldParameterSvdInitial); // instance of RooUnfoldSvd only used to find the best regularisation parameter; the unfolded spectrum returned by it is not retrived
  RooUnfold* unfold = unfoldBayes; // default Bayes
  TH1D* hist_unfold;

  int unfoldParameter;
  if (options.find("Svd") != std::string::npos) {
    unfoldSvdInitialiser->Hreco(); // necessary to have GetD() give a meaningful output
    TSVDUnfold *tsvdUnfold = (TSVDUnfold*)unfoldSvdInitialiser->Impl();
    if (automaticBestSvdParameter) {
      unfoldParameter = GetSvdBestRegularisationParameter_notYetSatisfying(tsvdUnfold);
    } else {
      unfoldParameter = unfoldParameterInput;
    }
    RooUnfoldSvd* unfoldSvd = new RooUnfoldSvd(response, measured, unfoldParameter); // the RooUnfoldSvd instance that is actually used to unfold, with the best regularisation parameter

    // unfoldSvd->SetRegParm(-200);  // doesnt seem to change anything, even if I put 0 here ... annoying; but accoring to mailing list of roounfold it doesn't work well
    hist_unfold = (TH1D*)(unfoldSvd->Hreco());
    unfold = unfoldSvd;

    // plot svd d distribution
    TH1D* H1D_D = tsvdUnfold->GetD();
    TString* pdfName_regparam = new TString("Svd_regularisationd_distribution_"+(TString)"_R="+Form("%.1f",arrayRadius[iRadius])+"_"+Datasets[iDataset]+DatasetsNames[iDataset]+"_"+(TString)mergingPrior+"_"+(TString)unfoldingPrior);
    TString textContext(contextCustomTwoFields(*texDatasetsComparisonCommonDenominator, contextJetRadius(arrayRadius[iRadius]), ""));
    std::array<std::array<float, 2>, 2> drawnWindowSvdParam = {{{0, 30}, {0.01, 10000}}}; // {{xmin, xmax}, {ymin, ymax}}
    Draw_TH1_Histogram(H1D_D, textContext, pdfName_regparam, texSvdDvector, texSvdK, texCollisionDataInfo, drawnWindowSvdParam, legendPlacementAuto, contextPlacementAuto, "logy");

  } else if (options.find("Bayes") != std::string::npos) {
    unfoldParameter = unfoldParameterInput;
    hist_unfold = (TH1D*)(unfoldBayes->Hreco());
    unfold = unfoldBayes;
  }  


  H1D_jetPt_unfolded = (TH1D*)hist_unfold->Clone("H1D_jetPt_unfolded"+partialUniqueSpecifier);

  bool divideSuccessEff;
  TH1D* H1D_jetEfficiency;
  if (!useFineBinningTest){ 
    divideSuccessEff = Get_Pt_JetEfficiency(H1D_jetEfficiency, iDataset, iRadius, options);
  } else {
    divideSuccessEff = Get_Pt_JetEfficiency_fineBinning(H1D_jetEfficiency, iDataset, iRadius, options);
  }

  if (useManualRespMatrixSettingMethod) {
    if (applyEfficiencies > 0) {
      if (applyEfficiencies == 2 || applyEfficiencies == 3) {
        if (divideSuccessEff){
          H1D_jetPt_unfolded->Divide(H1D_jetEfficiency);
        } else {
          cout << "################## Get_Pt_JetEfficiency FAILED!!!!! in Get_Pt_spectrum_unfolded_preWidthScalingAndEvtNorm ##################" << endl;
        }
      }
      if ((applyEfficiencies == 1 || applyEfficiencies == 3) && !useFineBinningTest) {
        H1D_jetPt_unfolded->Divide(H1D_kinematicEfficiency);
      }
    }
  }

  if (drawIntermediateResponseMatrices) {
    TH2D* H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined_postUnfolding = (TH2D*)unfold->response()->Hresponse()->Clone("H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined_postUnfolding"+partialUniqueSpecifier);

    TString* pdfName = new TString("responseMatrix_combined_postUnfolding"+(TString)"_R="+Form("%.1f",arrayRadius[iRadius])+"_"+Datasets[iDataset]+DatasetsNames[iDataset]);
    TString* pdfName_logz = new TString("responseMatrix_combined_postUnfolding"+(TString)"_R="+Form("%.1f",arrayRadius[iRadius])+"_"+Datasets[iDataset]+DatasetsNames[iDataset]+"_logz");

    TString textContext(contextCustomTwoFields(*texDatasetsComparisonCommonDenominator, contextJetRadius(arrayRadius[iRadius]), ""));

    Draw_TH2_Histogram(H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined_postUnfolding, textContext, pdfName, texPtJetRecX, texPtJetGenX, texCollisionDataInfo, drawnWindowAuto, th2ContoursNone, contourNumberNone, "");
    Draw_TH2_Histogram(H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined_postUnfolding, textContext, pdfName_logz, texPtJetRecX, texPtJetGenX, texCollisionDataInfo, drawnWindowAuto, th2ContoursNone, contourNumberNone, "logz");
  }


  // if (doWidthScalingEarly) {
  //   TransformRawHistToYield(H1D_jetPt_unfolded);
  // }

  std::pair<int, RooUnfold*> unfoldInfo(unfoldParameter, unfold);
  return unfoldInfo;
}
std::pair<int, RooUnfold*> Get_Pt_spectrum_unfolded_preWidthScaling(TH1D* &H1D_jetPt_unfolded, int iDataset, int iRadius, int unfoldParameterInput, std::string options) {
  std::pair<int, RooUnfold*> unfoldInfo = Get_Pt_spectrum_unfolded_preWidthScalingAndEvtNorm(H1D_jetPt_unfolded, iDataset, iRadius, unfoldParameterInput, options);
  cout << "Get_Pt_spectrum_unfolded_preWidthScaling test 1" << endl;
  if (normaliseDistribsAfterUnfolding){
    if (!controlMC) {
      NormaliseRawHistToNEvents(H1D_jetPt_unfolded, GetNEventsSelected_JetFramework(file_O2Analysis_list[iDataset], analysisWorkflowData));
    } else {
      if (!mcIsWeighted) {
        NormaliseRawHistToNEvents(H1D_jetPt_unfolded, GetNEventsSelected_JetFramework_weighted(file_O2Analysis_ppSimDetectorEffect_unfoldingControl, analysisWorkflow_unfoldingControl));
      } else {
        NormaliseRawHistToNEvents(H1D_jetPt_unfolded, GetNEventsSelected_JetFramework(file_O2Analysis_ppSimDetectorEffect_unfoldingControl, analysisWorkflow_unfoldingControl));
      }
    }
  }
   cout << "Get_Pt_spectrum_unfolded_preWidthScaling test 2" << endl;

  return unfoldInfo;
}
std::pair<int, RooUnfold*> Get_Pt_spectrum_unfolded(TH1D* &H1D_jetPt_unfolded, int iDataset, int iRadius, int unfoldParameterInput, std::string options) {
  std::pair<int, RooUnfold*> unfoldInfo = Get_Pt_spectrum_unfolded_preWidthScaling(H1D_jetPt_unfolded, iDataset, iRadius, unfoldParameterInput, options);

  if (doWidthScalingAtEnd) {
    TransformRawHistToYield(H1D_jetPt_unfolded);
  }
  return unfoldInfo;
}


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////// Spectrum Refolding tests //////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void Get_Pt_spectrum_dataUnfoldedThenRefolded_preWidthScalingAndEvtNorm(TH1D* &H1D_jetPt_unfoldedThenRefolded, int iDataset, int iRadius, int unfoldParameterInput, std::string options) {
  TH1D* H1D_jetPt_unfolded;
  // TH1D* H1D_jetPt_raw[nRadius];
  TString partialUniqueSpecifier = Datasets[iDataset]+"_R="+Form("%.1f",arrayRadius[iRadius]);


  RooUnfold* unfold = Get_Pt_spectrum_unfolded_preWidthScalingAndEvtNorm(H1D_jetPt_unfolded, iDataset, iRadius, unfoldParameterInput, options).second;
  // Get_Pt_spectrum_unfolded(H1D_jetPt_unfolded, iDataset, iRadius, unfoldParameterInput, options);

  // cout << "((((((((((((((((()))))))))))))))))" << endl;
  // cout << "REFOLDING TEST: " << endl;
  // for(int iBinX = 0; iBinX <= H1D_jetPt_unfolded->GetNbinsX()+1; iBinX++){
  //   cout << "H1D_jetPt_unfolded(" << iBinX << ") = " << H1D_jetPt_unfolded->GetBinContent(iBinX) << ", error = "<< H1D_jetPt_unfolded->GetBinError(iBinX) << endl;
  // }

  TH2D* H2D_jetPtResponseMatrix_fluctuations;
  TH2D* H2D_jetPtResponseMatrix_detectorResponse;
  TH2D* H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined;
  TH2D* H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined_fineBinning;
  TH1D* H1D_kinematicEfficiency;

  Get_PtResponseMatrix_Fluctuations(H2D_jetPtResponseMatrix_fluctuations, iDataset, iRadius);
  Get_PtResponseMatrix_detectorResponse(H2D_jetPtResponseMatrix_detectorResponse, iDataset, iRadius);

  Get_PtResponseMatrix_DetectorAndFluctuationsCombined_fineBinning(H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined_fineBinning, H2D_jetPtResponseMatrix_detectorResponse, H2D_jetPtResponseMatrix_fluctuations, iDataset, iRadius, options);
  Get_ResponseMatrix_Pt_KinematicEffiency(H1D_kinematicEfficiency, H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined_fineBinning, partialUniqueSpecifier, iRadius); // I want the efficiency before the reweighting and normalisation
  if (useFineBinningTest) {
    ReweightResponseMatrixWithPrior_fineBinning(H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined_fineBinning, iDataset, iRadius, options);
  }

  Get_PtResponseMatrix_DetectorAndFluctuationsCombined(H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined, H2D_jetPtResponseMatrix_detectorResponse, H2D_jetPtResponseMatrix_fluctuations, iDataset, iRadius, options);
  // FinaliseResponseMatrix(H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined, iDataset, iRadius, options);
  

  if (drawIntermediateResponseMatrices) {
    TH2D* H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined_duringRefolding = (TH2D*)H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined->Clone("H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined_duringRefolding"+partialUniqueSpecifier);

    TString* pdfName = new TString("responseMatrix_combined_duringRefolding"+(TString)"_R="+Form("%.1f",arrayRadius[iRadius])+"_"+Datasets[iDataset]+DatasetsNames[iDataset]);
    TString* pdfName_logz = new TString("responseMatrix_combined_duringRefolding"+(TString)"_R="+Form("%.1f",arrayRadius[iRadius])+"_"+Datasets[iDataset]+DatasetsNames[iDataset]+"_logz");

    TString textContext(contextCustomTwoFields(*texDatasetsComparisonCommonDenominator, contextJetRadius(arrayRadius[iRadius]), ""));

    Draw_TH2_Histogram(H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined_duringRefolding, textContext, pdfName, texPtJetRecX, texPtJetGenX, texCollisionDataInfo, drawnWindowAuto, th2ContoursNone, contourNumberNone, "");
    Draw_TH2_Histogram(H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined_duringRefolding, textContext, pdfName_logz, texPtJetRecX, texPtJetGenX, texCollisionDataInfo, drawnWindowAuto, th2ContoursNone, contourNumberNone, "logz");
  }

  cout << "((((((((((((((((()))))))))))))))))" << endl;
  cout << "REFOLDING TEST: pre efficiency" << endl;
  for(int iBinX = 0; iBinX <= H1D_jetPt_unfolded->GetNbinsX()+1; iBinX++){
    cout << "  H1D_jetPt_unfolded(" << iBinX << ") = " << H1D_jetPt_unfolded->GetBinContent(iBinX) << endl;
    cout << "e_H1D_jetPt_unfolded(" << iBinX << ") = " << H1D_jetPt_unfolded->GetBinError(iBinX) << endl;
  }
  

  bool divideSuccessEff;
  TH1D* H1D_jetEfficiency;
  if (!useFineBinningTest){ 
    divideSuccessEff = Get_Pt_JetEfficiency(H1D_jetEfficiency, iDataset, iRadius, options);
  } else {
    divideSuccessEff = Get_Pt_JetEfficiency_fineBinning(H1D_jetEfficiency, iDataset, iRadius, options);
  }

  if (applyEfficiencies > 0) {
    if (applyEfficiencies == 2 || applyEfficiencies == 3) {
      if (divideSuccessEff){
        H1D_jetPt_unfolded->Multiply(H1D_jetEfficiency);
      } else {
        cout << "################## H1D_jetPt_unfolded->Multiply(H1D_jetEfficiency) failed because Get_Pt_JetEfficiency() FAILED!!!!! ##################" << endl;
      }
    }
    if ((applyEfficiencies == 1 || applyEfficiencies == 3) && !useFineBinningTest) {
      H1D_jetPt_unfolded->Multiply(H1D_kinematicEfficiency);
    }
  }

  TH2D* refoldingResponseMatrix;
  if (useManualRespMatrixSettingMethod) {
    if (useFineBinningTest) {
      refoldingResponseMatrix = (TH2D*)H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined_fineBinning->Clone("refoldingResponseMatrix"+partialUniqueSpecifier);
    } else {
      refoldingResponseMatrix = (TH2D*)H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined->Clone("refoldingResponseMatrix"+partialUniqueSpecifier);
    }
  } else {
    refoldingResponseMatrix = (TH2D*)unfold->response()->HresponseNoOverflow();
    // refoldingResponseMatrix = (TH2D*)unfold->response()->Hresponse();
    cout << "refold why use HresponseNoOverflow instead of Hresponse?" << endl;
  }

  if (normaliseRespYSliceForRefold){
    NormaliseYSlicesToOne(refoldingResponseMatrix);
  }
  // H1D_jetPt_unfoldedThenRefolded = (TH1D*)GetMatrixVectorProductTH2xTH1((TH2D*)unfold->response()->Hresponse(), H1D_jetPt_unfolded).Clone("Get_Pt_spectrum_dataUnfoldedThenRefolded_preWidthScalingAndEvtNorm"+partialUniqueSpecifier); //gives the same
  // H1D_jetPt_unfoldedThenRefolded = (TH1D*)GetMatrixVectorProductTH2xTH1(H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined, H1D_jetPt_unfolded).Clone("Get_Pt_spectrum_dataUnfoldedThenRefolded_preWidthScalingAndEvtNorm"+partialUniqueSpecifier); //gives the same
  cout << "--UnfoldedJetPt in RooUnfold Refold Method--" << endl;
  H1D_jetPt_unfoldedThenRefolded = (TH1D*)GetMatrixVectorProductTH2xTH1(refoldingResponseMatrix, H1D_jetPt_unfolded).Clone("Get_Pt_spectrum_dataUnfoldedThenRefolded_preWidthScalingAndEvtNorm"+partialUniqueSpecifier); //gives the same


  bool divideSuccessFakes;
  TH1D* H1D_jetFakes;
  if (applyFakes) { // for useManualRespMatrixSettingMethod set to false, it will give a slightly different result as Get_Pt_JetFakes doesn't get fakes from the same histogram as the one used to encode fakes in the response matrix
    if (!useFineBinningTest){ 
      divideSuccessFakes = Get_Pt_JetFakes(H1D_jetFakes, iDataset, iRadius, options);
    } else {
      divideSuccessFakes = Get_Pt_JetFakes_fineBinning(H1D_jetFakes, iDataset, iRadius, options);
    }
    if (divideSuccessFakes){
      H1D_jetPt_unfoldedThenRefolded->Divide(H1D_jetFakes);
    } else {
      cout << "################## H1D_jetPt_unfoldedThenRefolded->Divide(H1D_jetFakes) failed because Get_Pt_JetFakes() FAILED!!!!! ##################" << endl;
    }
  }

  // cout << "((((((((((((((((()))))))))))))))))" << endl;
  // cout << "REFOLDING TEST: post efficiency" << endl;
  // for(int iBinX = 0; iBinX <= H1D_jetPt_unfoldedThenRefolded->GetNbinsX()+1; iBinX++){
  //   cout << "H1D_jetPt_unfoldedThenRefolded(" << iBinX << ") = " << H1D_jetPt_unfoldedThenRefolded->GetBinContent(iBinX) << ", error = "<< H1D_jetPt_unfoldedThenRefolded->GetBinError(iBinX) << endl;
  // }
  

  // cout << "Refolding check of Combined matrix = " << endl; 
  // for(int iBinY = 1; iBinY <= H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined->GetNbinsY(); iBinY++){
  //   cout << "combinedMatrix integral of slice iBinY = " << iBinY << " is: " << H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined->Integral(0, -1, iBinY, iBinY) << endl;
  // }


  // if (useFineBinningTest) {
  //   H1D_jetPt_unfoldedThenRefolded = (TH1D*)GetMatrixVectorProductTH2xTH1(H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined_fineBinning, H1D_jetPt_unfolded).Clone("Get_Pt_spectrum_dataUnfoldedThenRefolded_preWidthScalingAndEvtNorm"+partialUniqueSpecifier);
  // } else {
  //   H1D_jetPt_unfoldedThenRefolded = (TH1D*)GetMatrixVectorProductTH2xTH1(H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined, H1D_jetPt_unfolded).Clone("Get_Pt_spectrum_dataUnfoldedThenRefolded_preWidthScalingAndEvtNorm"+partialUniqueSpecifier);
  // }
  

  // if (doWidthScalingEarly) {
  //   TransformRawHistToYield(H1D_jetPt_unfoldedThenRefolded);
  // }

  // cout << "still not giving back the measured used as input to the unfolding; got an issue somewhere" << endl;
}
void Get_Pt_spectrum_dataUnfoldedThenRefolded_preWidthScaling(TH1D* &H1D_jetPt_unfoldedThenRefolded, int iDataset, int iRadius, int unfoldParameterInput, std::string options) {
  Get_Pt_spectrum_dataUnfoldedThenRefolded_preWidthScalingAndEvtNorm(H1D_jetPt_unfoldedThenRefolded, iDataset, iRadius, unfoldParameterInput, options);

  if (normaliseDistribsAfterUnfolding){
    if (!controlMC) {
      NormaliseRawHistToNEvents(H1D_jetPt_unfoldedThenRefolded, GetNEventsSelected_JetFramework(file_O2Analysis_list[iDataset], analysisWorkflowData));
    } else {
      if (!mcIsWeighted) {
        NormaliseRawHistToNEvents(H1D_jetPt_unfoldedThenRefolded, GetNEventsSelected_JetFramework_weighted(file_O2Analysis_ppSimDetectorEffect_unfoldingControl, analysisWorkflow_unfoldingControl));
      } else {
        NormaliseRawHistToNEvents(H1D_jetPt_unfoldedThenRefolded, GetNEventsSelected_JetFramework(file_O2Analysis_ppSimDetectorEffect_unfoldingControl, analysisWorkflow_unfoldingControl));
      }
    }
  }
}
void Get_Pt_spectrum_dataUnfoldedThenRefolded(TH1D* &H1D_jetPt_unfoldedThenRefolded, int iDataset, int iRadius, int unfoldParameterInput, std::string options) {
  Get_Pt_spectrum_dataUnfoldedThenRefolded_preWidthScaling(H1D_jetPt_unfoldedThenRefolded, iDataset, iRadius, unfoldParameterInput, options);

  if (doWidthScalingAtEnd) {
    TransformRawHistToYield(H1D_jetPt_unfoldedThenRefolded);
  }
}


void Get_Pt_spectrum_dataUnfoldedThenRefolded_RooUnfoldMethod_preWidthScalingAndEvtNorm(TH1D* &H1D_jetPt_unfoldedThenRefolded, int iDataset, int iRadius, int unfoldParameterInput, std::string options) {
  // ApplyToTruth function doesn't apply errors on folding

  // Matches exactly with the manual method IF NO PRIOR
  // if I have a non flat prior, then the roounfold method gives me a good closure test, but not the manual method!

  TH1D* H1D_jetPt_unfolded;
  // TH1D* H1D_jetPt_raw[nRadius];
  TString partialUniqueSpecifier = Datasets[iDataset]+"_R="+Form("%.1f",arrayRadius[iRadius]);


  RooUnfold* unfold = Get_Pt_spectrum_unfolded_preWidthScalingAndEvtNorm(H1D_jetPt_unfolded, iDataset, iRadius, unfoldParameterInput, options).second;

  TH2D* H2D_jetPtResponseMatrix_fluctuations;
  TH2D* H2D_jetPtResponseMatrix_detectorResponse;
  TH2D* H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined;
  TH2D* H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined_fineBinning;
  TH1D* H1D_kinematicEfficiency;

  Get_PtResponseMatrix_Fluctuations(H2D_jetPtResponseMatrix_fluctuations, iDataset, iRadius);
  Get_PtResponseMatrix_detectorResponse(H2D_jetPtResponseMatrix_detectorResponse, iDataset, iRadius);

  Get_PtResponseMatrix_DetectorAndFluctuationsCombined_fineBinning(H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined_fineBinning, H2D_jetPtResponseMatrix_detectorResponse, H2D_jetPtResponseMatrix_fluctuations, iDataset, iRadius, options);
  Get_ResponseMatrix_Pt_KinematicEffiency(H1D_kinematicEfficiency, H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined_fineBinning, partialUniqueSpecifier, iRadius); // I want the efficiency before the reweighting and normalisation

  Get_PtResponseMatrix_DetectorAndFluctuationsCombined(H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined, H2D_jetPtResponseMatrix_detectorResponse, H2D_jetPtResponseMatrix_fluctuations, iDataset, iRadius, options);
  // FinaliseResponseMatrix(H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined, iDataset, iRadius, options);

  bool divideSuccessEff;
  TH1D* H1D_jetEfficiency;
  if (!useFineBinningTest){ 
    divideSuccessEff = Get_Pt_JetEfficiency(H1D_jetEfficiency, iDataset, iRadius, options);
  } else {
    divideSuccessEff = Get_Pt_JetEfficiency_fineBinning(H1D_jetEfficiency, iDataset, iRadius, options);
  }

  if (useManualRespMatrixSettingMethod) { // this condition isn't present in the manual refolding method as in the manual case the response matrix object is just a th2 and not a roounfoldresponse object and doesn't have the efficiencies/fakes encoded in it
    if (applyEfficiencies > 0) {
      if (applyEfficiencies == 2 || applyEfficiencies == 3) {
        if (divideSuccessEff){
          H1D_jetPt_unfolded->Multiply(H1D_jetEfficiency);
        } else {
          cout << "################## H1D_jetPt_unfolded->Multiply(H1D_jetEfficiency) failed because Get_Pt_JetEfficiency() FAILED!!!!! ##################" << endl;
        }
      }
      if ((applyEfficiencies == 1 || applyEfficiencies == 3) && !useFineBinningTest) {
        H1D_jetPt_unfolded->Multiply(H1D_kinematicEfficiency);
      }
    }
  }

  // cout << "Refolding check of Combined matrix = " << endl; 
  // for(int iBinY = 1; iBinY <= H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined->GetNbinsY(); iBinY++){
  //   cout << "combinedMatrix integral of slice iBinY = " << iBinY << " is: " << H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined->Integral(0, -1, iBinY, iBinY) << endl;
  // }


  // RooUnfoldResponse* rooresponse_postUnfold = unfold->response();
  // H1D_jetPt_unfoldedThenRefolded = (TH1D*)unfold->response()->ApplyToTruth(H1D_jetPt_unfolded, "H1D_jetPt_unfoldedThenRefolded_withApplyToTruth")->Clone("Get_Pt_spectrum_dataUnfoldedThenRefolded_preWidthScalingAndEvtNorm"+partialUniqueSpecifier);

  // RooUnfoldResponse* response = new RooUnfoldResponse(0, 0, H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined); // measured and mcp_rebinned are here to take inneficiencies and fakes into account; or is it really what's happening? 'Alternatively, the response matrix can be constructed from a pre-existing TH2D 2-dimensional histogram (with truth and measured distribution TH1D histograms for normalisation).' from https://hepunx.rl.ac.uk/~adye/software/unfold/RooUnfold.html, and I'm already normalising so maybe tehre's no need for more normalisation
  H1D_jetPt_unfoldedThenRefolded = (TH1D*)unfold->response()->ApplyToTruth(H1D_jetPt_unfolded, "H1D_jetPt_unfoldedThenRefolded_withApplyToTruth")->Clone("Get_Pt_spectrum_dataUnfoldedThenRefolded_preWidthScalingAndEvtNorm"+partialUniqueSpecifier);


  cout << "--UnfoldedJetPt in RooUnfold Refold Method--" << endl;
  for(int iBinK = 1; iBinK <= H1D_jetPt_unfolded->GetNbinsX(); iBinK++){ // 0 and n+1 take underflow and overflow into account
    cout << "-----" << endl;
    cout << " U(" << iBinK << ") = " << H1D_jetPt_unfolded->GetBinContent(iBinK) << endl;
    cout << "eU(" << iBinK << ") = " << H1D_jetPt_unfolded->GetBinError(iBinK) << endl;
  }

  bool divideSuccessFakes;
  TH1D* H1D_jetFakes;
  if (useManualRespMatrixSettingMethod) { // this condition isn't present in the manual refolding method as in the manual case the response matrix object is just a th2 and not a roounfoldresponse object and doesn't have the efficiencies/fakes encoded in it
    if (applyFakes) {
      if (!useFineBinningTest){ 
        divideSuccessFakes = Get_Pt_JetFakes(H1D_jetFakes, iDataset, iRadius, options);
      } else {
        divideSuccessFakes = Get_Pt_JetFakes_fineBinning(H1D_jetFakes, iDataset, iRadius, options);
      }
      if (divideSuccessFakes){
        H1D_jetPt_unfoldedThenRefolded->Divide(H1D_jetFakes);
      } else {
        cout << "################## H1D_jetPt_unfoldedThenRefolded->Divide(H1D_jetFakes) failed because Get_Pt_JetFakes() FAILED!!!!! ##################" << endl;
      }
    }
  }

  //not sure I should normalise --> probably not as H1D_jetPt_unfolded is already normalised in Get_Pt_spectrum_unfolded_preWidthScalingAndEvtNorm

  // cout << "still not giving back the measured used as input to the unfolding; got an issue somewhere" << endl;
}
void Get_Pt_spectrum_dataUnfoldedThenRefolded_RooUnfoldMethod_preWidthScaling(TH1D* &H1D_jetPt_unfoldedThenRefolded, int iDataset, int iRadius, int unfoldParameterInput, std::string options) {
  Get_Pt_spectrum_dataUnfoldedThenRefolded_RooUnfoldMethod_preWidthScalingAndEvtNorm(H1D_jetPt_unfoldedThenRefolded, iDataset, iRadius, unfoldParameterInput, options);

  if (normaliseDistribsAfterUnfolding){
    if (!controlMC) {
      NormaliseRawHistToNEvents(H1D_jetPt_unfoldedThenRefolded, GetNEventsSelected_JetFramework(file_O2Analysis_list[iDataset], analysisWorkflowData));
    } else {
      if (!mcIsWeighted) {
        NormaliseRawHistToNEvents(H1D_jetPt_unfoldedThenRefolded, GetNEventsSelected_JetFramework_weighted(file_O2Analysis_ppSimDetectorEffect_unfoldingControl, analysisWorkflow_unfoldingControl));
      } else {
        NormaliseRawHistToNEvents(H1D_jetPt_unfoldedThenRefolded, GetNEventsSelected_JetFramework(file_O2Analysis_ppSimDetectorEffect_unfoldingControl, analysisWorkflow_unfoldingControl));
      }
    }
  }
}
void Get_Pt_spectrum_dataUnfoldedThenRefolded_RooUnfoldMethod(TH1D* &H1D_jetPt_unfoldedThenRefolded, int iDataset, int iRadius, int unfoldParameterInput, std::string options) {
  Get_Pt_spectrum_dataUnfoldedThenRefolded_RooUnfoldMethod_preWidthScaling(H1D_jetPt_unfoldedThenRefolded, iDataset, iRadius, unfoldParameterInput, options);

  if (doWidthScalingAtEnd) {
    TransformRawHistToYield(H1D_jetPt_unfoldedThenRefolded);
  }
}



// void Get_Pt_spectrum_mcpFoldedWithFluctuationsThenUnfolded_preWidthScalingAndEvtNorm(TH1D* &H1D_jetPt_mcpFoldedThenUnfolded, int iDataset, int iRadius, int unfoldParameterInput, std::string options) {

//   TString partialUniqueSpecifier = Datasets[iDataset]+"_R="+Form("%.1f",arrayRadius[iRadius]);

//   TH1D* H1D_jetPt_mcp_control;
//   Get_Pt_spectrum_mcp_genBinning(H1D_jetPt_mcp_control, iDataset, iRadius, true, options)
//   TH2D* refoldingResponseMatrix;
//   Get_PtResponseMatrix_Fluctuations(refoldingResponseMatrix, iDataset, iRadius);

//   TH1D* H1D_jetPt_mcp_control_folded = (TH1D*)GetMatrixVectorProductTH2xTH1(refoldingResponseMatrix, H1D_jetPt_mcp_control).Clone("Get_Pt_spectrum_mcpFoldedWithFluctuationsThenUnfolded_preWidthScalingAndEvtNorm"+partialUniqueSpecifier); //gives the same








//   RooUnfold* unfold = Get_Pt_spectrum_unfolded_preWidthScalingAndEvtNorm(H1D_jetPt_unfolded, iDataset, iRadius, unfoldParameterInput, options).second;
//   // Get_Pt_spectrum_unfolded(H1D_jetPt_unfolded, iDataset, iRadius, unfoldParameterInput, options);

//   // cout << "((((((((((((((((()))))))))))))))))" << endl;
//   // cout << "REFOLDING TEST: " << endl;
//   // for(int iBinX = 0; iBinX <= H1D_jetPt_unfolded->GetNbinsX()+1; iBinX++){
//   //   cout << "H1D_jetPt_unfolded(" << iBinX << ") = " << H1D_jetPt_unfolded->GetBinContent(iBinX) << ", error = "<< H1D_jetPt_unfolded->GetBinError(iBinX) << endl;
//   // }

//   TH2D* H2D_jetPtResponseMatrix_fluctuations;
//   TH2D* H2D_jetPtResponseMatrix_detectorResponse;
//   TH2D* H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined;
//   TH2D* H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined_fineBinning;
//   TH1D* H1D_kinematicEfficiency;

//   Get_PtResponseMatrix_Fluctuations(H2D_jetPtResponseMatrix_fluctuations, iDataset, iRadius);
//   Get_PtResponseMatrix_detectorResponse(H2D_jetPtResponseMatrix_detectorResponse, iDataset, iRadius);

//   Get_PtResponseMatrix_DetectorAndFluctuationsCombined_fineBinning(H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined_fineBinning, H2D_jetPtResponseMatrix_detectorResponse, H2D_jetPtResponseMatrix_fluctuations, iDataset, iRadius, options);
//   Get_ResponseMatrix_Pt_KinematicEffiency(H1D_kinematicEfficiency, H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined_fineBinning, partialUniqueSpecifier, iRadius); // I want the efficiency before the reweighting and normalisation
//   if (useFineBinningTest) {
//     ReweightResponseMatrixWithPrior_fineBinning(H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined_fineBinning, iDataset, iRadius, options);
//   }

//   Get_PtResponseMatrix_DetectorAndFluctuationsCombined(H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined, H2D_jetPtResponseMatrix_detectorResponse, H2D_jetPtResponseMatrix_fluctuations, iDataset, iRadius, options);
// //  FinaliseResponseMatrix(H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined, iDataset, iRadius, options);
  

//   if (drawIntermediateResponseMatrices) {
//     TH2D* H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined_duringRefolding = (TH2D*)H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined->Clone("H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined_duringRefolding"+partialUniqueSpecifier);

//     TString* pdfName = new TString("responseMatrix_combined_duringRefolding"+jetType[iJetType]+"_"+Datasets[iDataset]+DatasetsNames[iDataset]);
//     TString* pdfName_logz = new TString("responseMatrix_combined_duringRefolding"+(TString)"_R="+Form("%.1f",arrayRadius[iRadius])+"_"+Datasets[iDataset]+DatasetsNames[iDataset]+"_logz");

//     TString textContext(contextCustomOneField(*texDatasetsComparisonCommonDenominator, ""));

//     Draw_TH2_Histogram(H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined_duringRefolding, textContext, pdfName, texPtJetRecX, texPtJetGenX, texCollisionDataInfo, drawnWindowAuto, th2ContoursNone, contourNumberNone, "");
//     Draw_TH2_Histogram(H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined_duringRefolding, textContext, pdfName_logz, texPtJetRecX, texPtJetGenX, texCollisionDataInfo, drawnWindowAuto, th2ContoursNone, contourNumberNone, "logz");
//   }

//   cout << "((((((((((((((((()))))))))))))))))" << endl;
//   cout << "REFOLDING TEST: pre efficiency" << endl;
//   for(int iBinX = 0; iBinX <= H1D_jetPt_unfolded->GetNbinsX()+1; iBinX++){
//     cout << "  H1D_jetPt_unfolded(" << iBinX << ") = " << H1D_jetPt_unfolded->GetBinContent(iBinX) << endl;
//     cout << "e_H1D_jetPt_unfolded(" << iBinX << ") = " << H1D_jetPt_unfolded->GetBinError(iBinX) << endl;
//   }
  

//   bool divideSuccessEff;
//   TH1D* H1D_jetEfficiency;
//   if (!useFineBinningTest){ 
//     divideSuccessEff = Get_Pt_JetEfficiency(H1D_jetEfficiency, iDataset, iRadius, options);
//   } else {
//     divideSuccessEff = Get_Pt_JetEfficiency_fineBinning(H1D_jetEfficiency, iDataset, iRadius, options);
//   }

//   if (applyEfficiencies > 0) {
//     if (applyEfficiencies == 2 || applyEfficiencies == 3) {
//       if (divideSuccessEff){
//         H1D_jetPt_unfolded->Multiply(H1D_jetEfficiency);
//       } else {
//         cout << "################## H1D_jetPt_unfolded->Multiply(H1D_jetEfficiency) failed because Get_Pt_JetEfficiency() FAILED!!!!! ##################" << endl;
//       }
//     }
//     if ((applyEfficiencies == 1 || applyEfficiencies == 3) && !useFineBinningTest) {
//       H1D_jetPt_unfolded->Multiply(H1D_kinematicEfficiency);
//     }
//   }

//   TH2D* refoldingResponseMatrix;
//   if (useManualRespMatrixSettingMethod) {
//     if (useFineBinningTest) {
//       refoldingResponseMatrix = (TH2D*)H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined_fineBinning->Clone("refoldingResponseMatrix"+partialUniqueSpecifier);
//     } else {
//       refoldingResponseMatrix = (TH2D*)H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined->Clone("refoldingResponseMatrix"+partialUniqueSpecifier);
//     }
//   } else {
//     refoldingResponseMatrix = (TH2D*)unfold->response()->HresponseNoOverflow();
//     // refoldingResponseMatrix = (TH2D*)unfold->response()->Hresponse();
//     cout << "refold why use HresponseNoOverflow instead of Hresponse?" << endl;
//   }

//   if (normaliseRespYSliceForRefold){
//     NormaliseYSlicesToOne(refoldingResponseMatrix);
//   }
//   // H1D_jetPt_unfoldedThenRefolded = (TH1D*)GetMatrixVectorProductTH2xTH1((TH2D*)unfold->response()->Hresponse(), H1D_jetPt_unfolded).Clone("Get_Pt_spectrum_dataUnfoldedThenRefolded_preWidthScalingAndEvtNorm"+partialUniqueSpecifier); //gives the same
//   // H1D_jetPt_unfoldedThenRefolded = (TH1D*)GetMatrixVectorProductTH2xTH1(H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined, H1D_jetPt_unfolded).Clone("Get_Pt_spectrum_dataUnfoldedThenRefolded_preWidthScalingAndEvtNorm"+partialUniqueSpecifier); //gives the same
//   cout << "--UnfoldedJetPt in RooUnfold Refold Method--" << endl;
//   H1D_jetPt_unfoldedThenRefolded = (TH1D*)GetMatrixVectorProductTH2xTH1(refoldingResponseMatrix, H1D_jetPt_unfolded).Clone("Get_Pt_spectrum_dataUnfoldedThenRefolded_preWidthScalingAndEvtNorm"+partialUniqueSpecifier); //gives the same


//   bool divideSuccessFakes;
//   TH1D* H1D_jetFakes;
//   if (applyFakes) { // for useManualRespMatrixSettingMethod set to false, it will give a slightly different result as Get_Pt_JetFakes doesn't get fakes from the same histogram as the one used to encode fakes in the response matrix
//     if (!useFineBinningTest){ 
//       divideSuccessFakes = Get_Pt_JetFakes(H1D_jetFakes, iDataset, iRadius, options);
//     } else {
//       divideSuccessFakes = Get_Pt_JetFakes_fineBinning(H1D_jetFakes, iDataset, iRadius, options);
//     }
//     if (divideSuccessFakes){
//       H1D_jetPt_unfoldedThenRefolded->Divide(H1D_jetFakes);
//     } else {
//       cout << "################## H1D_jetPt_unfoldedThenRefolded->Divide(H1D_jetFakes) failed because Get_Pt_JetFakes() FAILED!!!!! ##################" << endl;
//     }
//   }
// }
// void Get_Pt_spectrum_mcpFoldedWithFluctuationsThenUnfolded_preWidthScaling(TH1D* &H1D_jetPt_mcpFoldedThenUnfolded, int iDataset, int iRadius, int unfoldParameterInput, std::string options) {
//   Get_Pt_spectrum_dataUnfoldedThenRefolded_RooUnfoldMethod_preWidthScalingAndEvtNorm(H1D_jetPt_mcpFoldedThenUnfolded, iDataset, iRadius, unfoldParameterInput, options);

//   if (normaliseDistribsAfterUnfolding){
//     // if (!controlMC) {
//       NormaliseRawHistToNEvents(H1D_jetPt_mcpFoldedThenUnfolded, GetNEventsSelected_JetFramework(file_O2Analysis_list[iDataset], analysisWorkflowData));
//     // } else {
//     //   if (!mcIsWeighted) {
//     //     NormaliseRawHistToNEvents(H1D_jetPt_mcpFoldedThenUnfolded, GetNEventsSelected_JetFramework_weighted(file_O2Analysis_ppSimDetectorEffect_unfoldingControl, analysisWorkflow_unfoldingControl));
//     //   } else {
//     //     NormaliseRawHistToNEvents(H1D_jetPt_mcpFoldedThenUnfolded, GetNEventsSelected_JetFramework(file_O2Analysis_ppSimDetectorEffect_unfoldingControl, analysisWorkflow_unfoldingControl));
//     //   }
//     // }
//   }
// }
// void Get_Pt_spectrum_mcpFoldedWithFluctuationsThenUnfolded(TH1D* &H1D_jetPt_mcpFoldedThenUnfolded, int iDataset, int iRadius, int unfoldParameterInput, std::string options) {
//   Get_Pt_spectrum_dataUnfoldedThenRefolded_RooUnfoldMethod_preWidthScaling(H1D_jetPt_mcpFoldedThenUnfolded, iDataset, iRadius, unfoldParameterInput, options);

//   if (doWidthScalingAtEnd) {
//     TransformRawHistToYield(H1D_jetPt_mcpFoldedThenUnfolded);
//   }
// }



#endif