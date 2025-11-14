#ifndef JETSPECTRU_RESPONSEMATRIXFUNCTIONS_C
#define JETSPECTRU_RESPONSEMATRIXFUNCTIONS_C

#include "JetSpectrum_ResponseMatrixFunctions.h"

#include "./JetSpectrum_SpectraGetters.h"
#include "./JetSpectrum_SpectraGetters.C"
#include "./JetSpectrum_Unfolding.h"
#include "./JetSpectrum_Unfolding.C"
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
//////////////////////////////////////////////////////////////////////////// Response matrix functions ///////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void Get_PtResponseMatrix_DetectorAndFluctuationsCombined_fineBinning(TH2D* &H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined_fineBinning, TH2D* H2D_jetPtResponseMatrix_detectorResponse, TH2D* H2D_jetPtResponseMatrix_fluctuations, int iDataset, int iRadius, __attribute__ ((unused)) std::string options) {
  // https://github.com/alisw/AliPhysics/blob/master/PWGJE/PWGJE/AliAnaChargedJetResponseMaker.cxx for ann example that works, by marta verveij

  TString partialUniqueSpecifier = Datasets[iDataset]+"_R="+Form("%.1f",arrayRadius[iRadius]);

  // matrix product of fluct response times det response; assumes the two are of the same size binning wise
  // Careful: xy of hist and ij of Resp(i,j) are inverted ! hist(j,i) = matrix(i,j) and so if matrix(i,j)=SUM(matrixA(i,k)matrixB(k,j)) then hist(j,i)=SUM(histA(k,i)histB(j,k)), and if we replace j,i by gen,rec we get hist(gen,rec)=SUM(histA(k,rec)histB(gen,k))
  H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined_fineBinning = (TH2D*)GetMatrixProductTH2xTH2(H2D_jetPtResponseMatrix_fluctuations, H2D_jetPtResponseMatrix_detectorResponse).Clone("Get_PtResponseMatrix_DetectorAndFluctuationsCombined_fineBinning"+partialUniqueSpecifier);

  if (drawIntermediateResponseMatrices) {
    struct stat st1{};
    if (stat("pdfFolder/ResponseMatrices", &st1) == -1) {
        mkdir("pdfFolder/ResponseMatrices", 0700);
    }
    struct stat st2{};
    if (stat("pngFolder/ResponseMatrices", &st2) == -1) {
        mkdir("pngFolder/ResponseMatrices", 0700);
    }
    
    TH2D* H2D_jetPtResponseMatrix_fineBinningPreTransforms = (TH2D*)H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined_fineBinning->Clone("H2D_jetPtResponseMatrix_fineBinningPreTransforms"+partialUniqueSpecifier);

    TString priorInfo = (TString)(TString)mergingPrior+"-"+(TString)unfoldingPrior;
    TString* pdfName_preRebin = new TString("ResponseMatrices/responseMatrix_combined_fineBinningPreTransforms"+(TString)"_R="+Form("%.1f",arrayRadius[iRadius])+"_"+Datasets[iDataset]+DatasetsNames[iDataset]+"_"+priorInfo);
    TString* pdfName_preRebin_logz = new TString("ResponseMatrices/responseMatrix_combined_fineBinningPreTransforms"+(TString)"_R="+Form("%.1f",arrayRadius[iRadius])+"_"+Datasets[iDataset]+DatasetsNames[iDataset]+"_"+priorInfo+"_logz");

    TString texCombinedMatrix = contextCustomOneField((TString)"Combined matrix - "+(TString)*texEnergy, "");
    TString textContextMatrixDetails = contextCustomFourFields((TString)"Detector response: "+(TString)*texCollisionMCType, "", (TString)"Fluctuations response: "+*texCollisionDataType, contextJetRadius(arrayRadius[iRadius]), "");

    // the matrix natural visualisation is actually the NON transposed histograms, rotated by 90° anti trigonometrically
    TH2D* MatrixResponse;
    TString* xLabel;
    TString* yLabel;
    if (transposeResponseHistogramsInDrawing) {
      MatrixResponse = (TH2D*)GetTransposeHistogram(H2D_jetPtResponseMatrix_fineBinningPreTransforms).Clone("responseMatrix_combined_fineBinningPreTransforms"+(TString)"_R="+Form("%.1f",arrayRadius[iRadius])+"_"+Datasets[iDataset]+DatasetsNames[iDataset]+"_"+priorInfo);
      xLabel = texPtJetGen;
      yLabel = texPtJetRec;
    } else {
      MatrixResponse = (TH2D*)H2D_jetPtResponseMatrix_fineBinningPreTransforms->Clone("responseMatrix_combined_fineBinningPreTransforms"+(TString)"_R="+Form("%.1f",arrayRadius[iRadius])+"_"+Datasets[iDataset]+DatasetsNames[iDataset]+"_"+priorInfo);
      xLabel = texPtJetRec;
      yLabel = texPtJetGen;
    }

    Draw_TH2_Histogram(MatrixResponse, textContextMatrixDetails, pdfName_preRebin, xLabel, yLabel, &texCombinedMatrix, drawnWindow2DAuto, th2ContoursNone, contourNumberNone, "");
    Draw_TH2_Histogram(MatrixResponse, textContextMatrixDetails, pdfName_preRebin_logz, xLabel, yLabel, &texCombinedMatrix, drawnWindow2DAuto, th2ContoursNone, contourNumberNone, "logz");
  }
  // cout << "bin(topleft 1,N) = " << H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined_fineBinning->GetBinContent(1,H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined_fineBinning->GetNbinsY()) << endl;
  // cout << "bin(bottom left 1,1) = " << H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined_fineBinning->GetBinContent(1,1) << endl;
}


void Get_PtResponseMatrix_DetectorAndFluctuationsCombined(TH2D* &H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined, TH2D* H2D_jetPtResponseMatrix_detectorResponse, TH2D* H2D_jetPtResponseMatrix_fluctuations, int iDataset, int iRadius, std::string options) {
  // https://github.com/alisw/AliPhysics/blob/master/PWGJE/PWGJE/AliAnaChargedJetResponseMaker.cxx for ann example that works, by marta verveij

  // function should be removed when time allows as it doesn't do anything more than that Get_PtResponseMatrix_DetectorAndFluctuationsCombined_fineBinning()

  TString partialUniqueSpecifier = Datasets[iDataset]+"_R="+Form("%.1f",arrayRadius[iRadius]);

  Get_PtResponseMatrix_DetectorAndFluctuationsCombined_fineBinning(H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined, H2D_jetPtResponseMatrix_detectorResponse, H2D_jetPtResponseMatrix_fluctuations, iDataset, iRadius, options);
  FinaliseResponseMatrix(H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined, iDataset, iRadius, options);
}

void ReweightResponseMatrixWithPrior_modular(TH2D* &H2D_jetPtResponseMatrix, int iDataset, int iRadius, std::string options) {
   TString partialUniqueSpecifier = Datasets[iDataset]+"_R="+Form("%.1f",arrayRadius[iRadius]);
 
  // before this, all y-slices (ie pt gen slices) have been normalised to 1;means each pt gen slice has a proba of 1
  // withthis function, we give each slice a weight so that they have different normalisation values, corresponding to the prior 

  // TH2D* H2D_jetPtResponseMatrix_preReweightWithPrior = H2D_jetPtResponseMatrix->Clone();

  // prior choice; none by default (flat)
  TH1D* priorSpectrumWeighting;
  if (options.find("mcpPriorUnfolding") != std::string::npos) {
    if (!normGenAndMeasByNEvtsBeforeUnfolding) {
      if (matrixTransformationOrder == 0 || matrixTransformationOrder == 3) {
        Get_Pt_spectrum_mcp_fineBinning_preWidthScalingAtEndAndEvtNorm(priorSpectrumWeighting, iDataset, iRadius, false, options); 
      } else {
        Get_Pt_spectrum_mcp_genBinning_preWidthScalingAtEndAndEvtNorm(priorSpectrumWeighting, iDataset, iRadius, false, options); 
      }
    } else {
       if (matrixTransformationOrder == 0 || matrixTransformationOrder == 3) {
        Get_Pt_spectrum_mcp_fineBinning_preWidthScalingAtEnd(priorSpectrumWeighting, iDataset, iRadius, false, options); 
      } else {
        Get_Pt_spectrum_mcp_genBinning_preWidthScalingAtEnd(priorSpectrumWeighting, iDataset, iRadius, false, options); 
      }
    }
    WeightMatrixWithPrior(H2D_jetPtResponseMatrix, priorSpectrumWeighting, doUnfoldingPriorDivision);
    // for (int i = 1; i < priorSpectrumWeighting->GetNbinsX(); i++)
    // {
    //   cout << "prior(" << i << ")" << priorSpectrumWeighting->GetBinContent(i)<< endl;
    //   cout << "responseIntegralLine(" << i << ")" << H2D_jetPtResponseMatrix->Integral(1,H2D_jetPtResponseMatrix->GetNbinsX(), i, i)<< endl;
    // }
    
  }
  if (options.find("mcdPriorUnfolding") != std::string::npos) {
    if (!normGenAndMeasByNEvtsBeforeUnfolding) {
      if (matrixTransformationOrder == 0 || matrixTransformationOrder == 3) {
        Get_Pt_spectrum_mcd_fineBinning_preWidthScalingAtEndAndEvtNorm(priorSpectrumWeighting, iDataset, iRadius, options);
      } else {
        Get_Pt_spectrum_mcd_genBinning_preWidthScalingAtEndAndEvtNorm(priorSpectrumWeighting, iDataset, iRadius, options);
      } 
    } else {
      if (matrixTransformationOrder == 0 || matrixTransformationOrder == 3) {
        Get_Pt_spectrum_mcd_fineBinning_preWidthScalingAtEnd(priorSpectrumWeighting, iDataset, iRadius, options);
      } else {
        Get_Pt_spectrum_mcd_genBinning_preWidthScalingAtEnd(priorSpectrumWeighting, iDataset, iRadius, options); 
      }
    }
    WeightMatrixWithPrior(H2D_jetPtResponseMatrix, priorSpectrumWeighting, doUnfoldingPriorDivision);
  }
  if (options.find("measuredPriorUnfolding") != std::string::npos) {
    if (!normGenAndMeasByNEvtsBeforeUnfolding) {
      if (matrixTransformationOrder == 0 || matrixTransformationOrder == 3) {
        Get_Pt_spectrum_bkgCorrected_fineBinning_preWidthScalingAtEndAndEvtNorm(priorSpectrumWeighting, iDataset, iRadius, options);
      } else {
        Get_Pt_spectrum_bkgCorrected_genBinning_preWidthScalingAtEndAndEvtNorm(priorSpectrumWeighting, iDataset, iRadius, options);
      }
    } else {
      if (matrixTransformationOrder == 0 || matrixTransformationOrder == 3) {
        Get_Pt_spectrum_bkgCorrected_fineBinning_preWidthScalingAtEnd(priorSpectrumWeighting, iDataset, iRadius, options);
      } else {
        Get_Pt_spectrum_bkgCorrected_genBinning_preWidthScalingAtEnd(priorSpectrumWeighting, iDataset, iRadius, options); 
      }
    }
    WeightMatrixWithPrior(H2D_jetPtResponseMatrix, priorSpectrumWeighting, doUnfoldingPriorDivision);
  }
  if (options.find("testAliPhysics") != std::string::npos) {
    if (!normGenAndMeasByNEvtsBeforeUnfolding) {
      if (matrixTransformationOrder == 0 || matrixTransformationOrder == 3) {
        Get_Pt_spectrum_bkgCorrected_fineBinning_preWidthScalingAtEndAndEvtNorm(priorSpectrumWeighting, iDataset, iRadius, options);
      } else {
        Get_Pt_spectrum_bkgCorrected_genBinning_preWidthScalingAtEndAndEvtNorm(priorSpectrumWeighting, iDataset, iRadius, options); 
      }
    } else {
      if (matrixTransformationOrder == 0 || matrixTransformationOrder == 3) {
        Get_Pt_spectrum_bkgCorrected_fineBinning_preWidthScalingAtEnd(priorSpectrumWeighting, iDataset, iRadius, options);
      } else {
        Get_Pt_spectrum_bkgCorrected_genBinning_preWidthScalingAtEnd(priorSpectrumWeighting, iDataset, iRadius, options); 
      }
    }
    H2D_jetPtResponseMatrix = (TH2D*)NormalizeResponsMatrixYaxisWithPrior(H2D_jetPtResponseMatrix, priorSpectrumWeighting)->Clone(H2D_jetPtResponseMatrix->GetName()+(TString)"_testAliPhysics");
  }
  // cout << "((((((((((((((((((((((((()))))))))))))))))))))))))" << endl;
  // cout << "pre norm that shouldn't be" << endl;
  // cout << "H2D_jetPtResponseMatrix->Integral(1, N, 1, 1)" << H2D_jetPtResponseMatrix->Integral(1, H2D_jetPtResponseMatrix->GetNbinsX(), 1, 1) << endl;
  // cout << "H2D_jetPtResponseMatrix->Integral(0, -1, 1, 1)" << H2D_jetPtResponseMatrix->Integral(0, -1, 1, 1) << endl;
  // // NormaliseXSlicesToOneNoUnderOverFlows(H2D_jetPtResponseMatrix);
  // // NormaliseYSlicesToOneNoUnderOverFlows(H2D_jetPtResponseMatrix);
  // cout << "post norm that shouldn't be" << endl;
  // cout << "H2D_jetPtResponseMatrix->Integral(1, N, 1, 1)" << H2D_jetPtResponseMatrix->Integral(1, H2D_jetPtResponseMatrix->GetNbinsX(), 1, 1) << endl;
  // cout << "H2D_jetPtResponseMatrix->Integral(0, -1, 1, 1)" << H2D_jetPtResponseMatrix->Integral(0, -1, 1, 1) << endl;
  // cout << "((((((((((((((((((((((((()))))))))))))))))))))))))" << endl;


  if (drawIntermediateResponseMatrices) {
    TH2D* H2D_jetPtResponseMatrix_postReweightWithPrior = (TH2D*)H2D_jetPtResponseMatrix->Clone("H2D_jetPtResponseMatrix_postReweightWithPrior"+partialUniqueSpecifier);

    struct stat st1{};
    if (stat("pdfFolder/ResponseMatrices", &st1) == -1) {
        mkdir("pdfFolder/ResponseMatrices", 0700);
    }
    struct stat st2{};
    if (stat("pngFolder/ResponseMatrices", &st2) == -1) {
        mkdir("pngFolder/ResponseMatrices", 0700);
    }
    
    TString priorInfo = (TString)(TString)mergingPrior+"-"+(TString)unfoldingPrior;
    TString* pdfName_preRebin = new TString("ResponseMatrices/responseMatrix_combined_postReweightWithPrior"+(TString)"_R="+Form("%.1f",arrayRadius[iRadius])+"_"+Datasets[iDataset]+DatasetsNames[iDataset]+"_"+priorInfo);
    TString* pdfName_preRebin_logz = new TString("ResponseMatrices/responseMatrix_combined_postReweightWithPrior"+(TString)"_R="+Form("%.1f",arrayRadius[iRadius])+"_"+Datasets[iDataset]+DatasetsNames[iDataset]+"_"+priorInfo+"_logz");


    TString texCombinedMatrix = contextCustomOneField((TString)"Combined matrix - "+(TString)*texEnergy, "");
    TString textContextMatrixDetails = contextCustomFourFields((TString)"Detector response: "+(TString)*texCollisionMCType, "", (TString)"Fluctuations response: "+*texCollisionDataType, contextJetRadius(arrayRadius[iRadius]), "");

    // the matrix natural visualisation is actually the NON transposed histograms, rotated by 90° anti trigonometrically
    TH2D* MatrixResponse;
    TString* xLabel;
    TString* yLabel;
    if (transposeResponseHistogramsInDrawing) {
      MatrixResponse = (TH2D*)GetTransposeHistogram(H2D_jetPtResponseMatrix_postReweightWithPrior).Clone("responseMatrix_combined_postReweightWithPrior"+(TString)"_R="+Form("%.1f",arrayRadius[iRadius])+"_"+Datasets[iDataset]+DatasetsNames[iDataset]+"_"+priorInfo);
      xLabel = texPtJetGen;
      yLabel = texPtJetRec;
    } else {
      MatrixResponse = (TH2D*)H2D_jetPtResponseMatrix_postReweightWithPrior->Clone("responseMatrix_combined_postReweightWithPrior"+(TString)"_R="+Form("%.1f",arrayRadius[iRadius])+"_"+Datasets[iDataset]+DatasetsNames[iDataset]+"_"+priorInfo);
      xLabel = texPtJetRec;
      yLabel = texPtJetGen;
    }

    Draw_TH2_Histogram(MatrixResponse, textContextMatrixDetails, pdfName_preRebin, xLabel, yLabel, &texCombinedMatrix, drawnWindow2DAuto, th2ContoursNone, contourNumberNone, "");
    Draw_TH2_Histogram(MatrixResponse, textContextMatrixDetails, pdfName_preRebin_logz, xLabel, yLabel, &texCombinedMatrix, drawnWindow2DAuto, th2ContoursNone, contourNumberNone, "logz");
  }
}



void MergeResponseMatrixBins(TH2D* &H2D_jetPtResponseMatrix, int iDataset, int iRadius, std::string options) {
  // Merge bins of the combined response matrix with fine binning to get the coarse one
  TH1D* priorSpectrumMerging;
  bool debugBool = false;

  TH2D* H2D_jetPtResponseMatrix_postBinMerge;

  TString partialUniqueSpecifier = Datasets[iDataset]+"_R="+Form("%.1f",arrayRadius[iRadius]);

  Get_Pt_spectrum_mcp_fineBinning(priorSpectrumMerging, iDataset, iRadius, false, options); //take mcp as prior by default
  if (options.find("mcpPriorMerging") != std::string::npos) {
    priorSpectrumMerging->Reset("M");
    Get_Pt_spectrum_mcp_fineBinning(priorSpectrumMerging, iDataset, iRadius, false, options); 
    H2D_jetPtResponseMatrix_postBinMerge = (TH2D*)RebinVariableBins2D_PriorWeightedBinMerging(H2D_jetPtResponseMatrix, nBinPtJetsRec[iRadius], nBinPtJetsGen[iRadius], ptBinsJetsRec[iRadius], ptBinsJetsGen[iRadius], priorSpectrumMerging, debugBool).Clone("MergeResponseMatrixBins"+partialUniqueSpecifier);
  }
  if (options.find("mcdPriorMerging") != std::string::npos) {
    priorSpectrumMerging->Reset("M");
    Get_Pt_spectrum_mcd_fineBinning(priorSpectrumMerging, iDataset, iRadius, options);
    H2D_jetPtResponseMatrix_postBinMerge = (TH2D*)RebinVariableBins2D_PriorWeightedBinMerging(H2D_jetPtResponseMatrix, nBinPtJetsRec[iRadius], nBinPtJetsGen[iRadius], ptBinsJetsRec[iRadius], ptBinsJetsGen[iRadius], priorSpectrumMerging, debugBool).Clone("MergeResponseMatrixBins"+partialUniqueSpecifier);
  }
  if (options.find("measuredPriorMerging") != std::string::npos) {
    priorSpectrumMerging->Reset("M");
    Get_Pt_spectrum_bkgCorrected_fineBinning(priorSpectrumMerging, iDataset, iRadius, options);
    H2D_jetPtResponseMatrix_postBinMerge = (TH2D*)RebinVariableBins2D_PriorWeightedBinMerging(H2D_jetPtResponseMatrix, nBinPtJetsRec[iRadius], nBinPtJetsGen[iRadius], ptBinsJetsRec[iRadius], ptBinsJetsGen[iRadius], priorSpectrumMerging, debugBool).Clone("MergeResponseMatrixBins"+partialUniqueSpecifier);
  }
  if (options.find("noPriorMerging") != std::string::npos) {
    H2D_jetPtResponseMatrix_postBinMerge = (TH2D*)RebinVariableBins2D(H2D_jetPtResponseMatrix, nBinPtJetsRec[iRadius], nBinPtJetsGen[iRadius], ptBinsJetsRec[iRadius], ptBinsJetsGen[iRadius], debugBool).Clone("MergeResponseMatrixBins"+partialUniqueSpecifier);
  }
  if (options.find("testAliPhysics") != std::string::npos) {
    H2D_jetPtResponseMatrix_postBinMerge = (TH2D*)RebinVariableBins2D_aliPhysics(H2D_jetPtResponseMatrix, nBinPtJetsRec[iRadius], nBinPtJetsGen[iRadius], ptBinsJetsRec[iRadius], ptBinsJetsGen[iRadius], debugBool)->Clone("MergeResponseMatrixBins"+partialUniqueSpecifier);
  }
  // normalising priorSpectrum with evtNorm doesn't change anything as the weighting does prior_content(i)/prior_integral()
  // dividing priorSpectrum by binwidth doesn't change anything for the same reason

  H2D_jetPtResponseMatrix = (TH2D*)H2D_jetPtResponseMatrix_postBinMerge->Clone("RespMatrix_MergeResponseMatrixBins_output_"+partialUniqueSpecifier);
  // H2D_jetPtResponseMatrix_postBinMerge = (TH2D*)RebinVariableBins2D(H2D_jetPtResponseMatrix_preRebin, nBinPtJetsRec[iRadius], nBinPtJetsGen[iRadius], ptBinsJetsRec[iRadius], ptBinsJetsGen[iRadius], true).Clone("Get_PtResponseMatrix_DetectorAndFluctuationsCombined"+partialUniqueSpecifier);


  // When looking at combined response matrix before normalisation, large bins in y will look strange, and out of place compared to other bin slices of same size; this is because it potentially merges A LOT of bins together; it'll look a lot better after normalisation:

  if (!useMatrixOverflows){
    H2D_jetPtResponseMatrix->SetBinContent(0, 0);
    H2D_jetPtResponseMatrix->SetBinError(0, 0);
    H2D_jetPtResponseMatrix->SetBinContent(nBinPtJetsRec[iRadius], 0);
    H2D_jetPtResponseMatrix->SetBinError(nBinPtJetsRec[iRadius], 0);
  }

  // debug
  // for(int iBinY = 0; iBinY <= H2D_jetPtResponseMatrix_detectorResponse->GetNbinsY()+1; iBinY++){ // 0 and n+1 take underflow and overflow into account
  //   for(int iBinX = 0; iBinX <= H2D_jetPtResponseMatrix_detectorResponse->GetNbinsX()+1; iBinX++){ // 0 and n+1 take underflow and overflow into account
  //     if (iBinX == 0) {
  //       cout << "iBinX = 0 --> H2D_jetPtResponseMatrix_detectorResponse->GetBinContent(0,"<< iBinY << ") = " << H2D_jetPtResponseMatrix_detectorResponse->GetBinContent(0,iBinY) << endl;
  //       cout << "iBinX = 0 --> H2D_jetPtResponseMatrix_fluctuations->GetBinContent(0,"<< iBinY << ") = " << H2D_jetPtResponseMatrix_fluctuations->GetBinContent(0,iBinY) << endl;
  //     }
  //   }
  // }
  if (drawIntermediateResponseMatrices) {
    // TH2D* H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined_postBinMerge = (TH2D*)H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined->Clone("H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined_postBinMerge"+partialUniqueSpecifier);

    struct stat st1{};
    if (stat("pdfFolder/ResponseMatrices", &st1) == -1) {
        mkdir("pdfFolder/ResponseMatrices", 0700);
    }
    struct stat st2{};
    if (stat("pngFolder/ResponseMatrices", &st2) == -1) {
        mkdir("pngFolder/ResponseMatrices", 0700);
    }

    TString priorInfo = (TString)(TString)mergingPrior+"-"+(TString)unfoldingPrior;
    TString* pdfName = new TString("ResponseMatrices/responseMatrix_combined_postBinMerge"+(TString)"_R="+Form("%.1f",arrayRadius[iRadius])+"_"+Datasets[iDataset]+DatasetsNames[iDataset]+"_"+priorInfo);
    TString* pdfName_logz = new TString("ResponseMatrices/responseMatrix_combined_postBinMerge"+(TString)"_R="+Form("%.1f",arrayRadius[iRadius])+"_"+Datasets[iDataset]+DatasetsNames[iDataset]+"_"+priorInfo+"_logz");


    TString texCombinedMatrix = contextCustomOneField((TString)"Combined matrix - "+(TString)*texEnergy, "");
    TString textContextMatrixDetails = contextCustomFourFields((TString)"Detector response: "+(TString)*texCollisionMCType, "", (TString)"Fluctuations response: "+*texCollisionDataType, contextJetRadius(arrayRadius[iRadius]), "");

    // the matrix natural visualisation is actually the NON transposed histograms, rotated by 90° anti trigonometrically
    TH2D* MatrixResponse;
    TString* xLabel;
    TString* yLabel;
    if (transposeResponseHistogramsInDrawing) {
      MatrixResponse = (TH2D*)GetTransposeHistogram(H2D_jetPtResponseMatrix_postBinMerge).Clone("responseMatrix_combined_postReweightWithPrior"+(TString)"_R="+Form("%.1f",arrayRadius[iRadius])+"_"+Datasets[iDataset]+DatasetsNames[iDataset]+"_"+priorInfo);
      xLabel = texPtJetGen;
      yLabel = texPtJetRec;
    } else {
      MatrixResponse = (TH2D*)H2D_jetPtResponseMatrix_postBinMerge->Clone("responseMatrix_combined_postReweightWithPrior"+(TString)"_R="+Form("%.1f",arrayRadius[iRadius])+"_"+Datasets[iDataset]+DatasetsNames[iDataset]+"_"+priorInfo);
      xLabel = texPtJetRec;
      yLabel = texPtJetGen;
    }

    Draw_TH2_Histogram(MatrixResponse, textContextMatrixDetails, pdfName, xLabel, yLabel, &texCombinedMatrix, drawnWindow2DAuto, th2ContoursNone, contourNumberNone, "");
    Draw_TH2_Histogram(MatrixResponse, textContextMatrixDetails, pdfName_logz, xLabel, yLabel, &texCombinedMatrix, drawnWindow2DAuto, th2ContoursNone, contourNumberNone, "logz");
  }

  // TransformRawResponseToYieldResponse(H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined);
  // cout << "Should I normalise the combined matrix?" << endl;
  // cout << "     Marta doesn't do it, but it looks like I need it due to merging (some bins are very large)" << endl;
  // cout << "     it's actually done in AliAnaChargedJetResponseMaker::MakeResponseMatrixRebin" << endl;
}



void NormYSlicesAndScaleRespByWidth(TH2D* &H2D_jetPtResponseMatrix, int iDataset, int iRadius, std::string options) {
  TString partialUniqueSpecifier = Datasets[iDataset]+"_R="+Form("%.1f",arrayRadius[iRadius]);

  if (doYSliceNormToOneCombinedResp) {
    NormaliseYSlicesToOne(H2D_jetPtResponseMatrix); // Marta doesn't do it, but it looks like I need it due to merging (some bins are very large) ; actually marta probably uses it: AliAnaChargedJetResponseMaker::MakeResponseMatrixRebin does it by default inside the rebinning function, and it takes into account the whole Fine range (1fine, Nfine)
  } 
  if (scaleRespByXYWidth) {
    H2D_jetPtResponseMatrix->Scale(1., "width");
  }
  if (scaleRespByYWidth) {
    double binWidthY;
    for(int iBinGen = 1; iBinGen <= H2D_jetPtResponseMatrix->GetNbinsY(); iBinGen++){
      binWidthY = H2D_jetPtResponseMatrix->GetYaxis()->GetBinWidth(iBinGen);
      for(int iBinRec = 1; iBinRec <= H2D_jetPtResponseMatrix->GetNbinsX(); iBinRec++){
        H2D_jetPtResponseMatrix->SetBinContent(iBinRec, iBinGen, 1./binWidthY * H2D_jetPtResponseMatrix->GetBinContent(iBinRec, iBinGen));
        H2D_jetPtResponseMatrix->SetBinError(iBinRec, iBinGen, 1./binWidthY * H2D_jetPtResponseMatrix->GetBinError(iBinRec, iBinGen));
      }
    }
  }

  if (drawIntermediateResponseMatrices) {
    TH2D* H2D_jetPtResponseMatrix_postYSliceNormAndWidthNorm = (TH2D*)H2D_jetPtResponseMatrix->Clone("NormYSlicesAndScaleRespByWidth"+partialUniqueSpecifier);

    struct stat st1{};
    if (stat("pdfFolder/ResponseMatrices", &st1) == -1) {
        mkdir("pdfFolder/ResponseMatrices", 0700);
    }
    struct stat st2{};
    if (stat("pngFolder/ResponseMatrices", &st2) == -1) {
        mkdir("pngFolder/ResponseMatrices", 0700);
    }

    TString priorInfo = (TString)(TString)mergingPrior+"-"+(TString)unfoldingPrior;
    TString* pdfName = new TString("ResponseMatrices/responseMatrix_combined_postYSliceNormAndWidthNorm"+(TString)"_R="+Form("%.1f",arrayRadius[iRadius])+"_"+Datasets[iDataset]+DatasetsNames[iDataset]+"_"+priorInfo);
    TString* pdfName_logz = new TString("ResponseMatrices/responseMatrix_combined_postYSliceNormAndWidthNorm"+(TString)"_R="+Form("%.1f",arrayRadius[iRadius])+"_"+Datasets[iDataset]+DatasetsNames[iDataset]+"_"+priorInfo+"_logz");


    TString texCombinedMatrix = contextCustomOneField((TString)"Combined matrix - "+(TString)*texEnergy, "");
    TString textContextMatrixDetails = contextCustomFourFields((TString)"Detector response: "+(TString)*texCollisionMCType, "", (TString)"Fluctuations response: "+*texCollisionDataType, contextJetRadius(arrayRadius[iRadius]), "");

    // the matrix natural visualisation is actually the NON transposed histograms, rotated by 90° anti trigonometrically
    TH2D* MatrixResponse;
    TString* xLabel;
    TString* yLabel;
    if (transposeResponseHistogramsInDrawing) {
      MatrixResponse = (TH2D*)GetTransposeHistogram(H2D_jetPtResponseMatrix_postYSliceNormAndWidthNorm).Clone("responseMatrix_combined_postReweightWithPrior"+(TString)"_R="+Form("%.1f",arrayRadius[iRadius])+"_"+Datasets[iDataset]+DatasetsNames[iDataset]+"_"+priorInfo);
      xLabel = texPtJetGen;
      yLabel = texPtJetRec;
    } else {
      MatrixResponse = (TH2D*)H2D_jetPtResponseMatrix_postYSliceNormAndWidthNorm->Clone("responseMatrix_combined_postReweightWithPrior"+(TString)"_R="+Form("%.1f",arrayRadius[iRadius])+"_"+Datasets[iDataset]+DatasetsNames[iDataset]+"_"+priorInfo);
      xLabel = texPtJetRec;
      yLabel = texPtJetGen;
    }

    Draw_TH2_Histogram(MatrixResponse, textContextMatrixDetails, pdfName, xLabel, yLabel, &texCombinedMatrix, drawnWindow2DAuto, th2ContoursNone, contourNumberNone, "");
    Draw_TH2_Histogram(MatrixResponse, textContextMatrixDetails, pdfName_logz, xLabel, yLabel, &texCombinedMatrix, drawnWindow2DAuto, th2ContoursNone, contourNumberNone, "logz");
  }
}

void FinaliseResponseMatrix(TH2D* &H2D_jetPtResponseMatrix, int iDataset, int iRadius, std::string options) {
  if (matrixTransformationOrder == 0) {
    ReweightResponseMatrixWithPrior_modular(H2D_jetPtResponseMatrix, iDataset, iRadius, options);
    MergeResponseMatrixBins(H2D_jetPtResponseMatrix, iDataset, iRadius, options);
    NormYSlicesAndScaleRespByWidth(H2D_jetPtResponseMatrix, iDataset, iRadius, options);
  } else if (matrixTransformationOrder == 1) {
    MergeResponseMatrixBins(H2D_jetPtResponseMatrix, iDataset, iRadius, options);
    NormYSlicesAndScaleRespByWidth(H2D_jetPtResponseMatrix, iDataset, iRadius, options);
    ReweightResponseMatrixWithPrior_modular(H2D_jetPtResponseMatrix, iDataset, iRadius, options);
  } else if (matrixTransformationOrder == 2) {
    MergeResponseMatrixBins(H2D_jetPtResponseMatrix, iDataset, iRadius, options);
    ReweightResponseMatrixWithPrior_modular(H2D_jetPtResponseMatrix, iDataset, iRadius, options);
    NormYSlicesAndScaleRespByWidth(H2D_jetPtResponseMatrix, iDataset, iRadius, options);
  } else if (matrixTransformationOrder == 3) {
    ReweightResponseMatrixWithPrior_modular(H2D_jetPtResponseMatrix, iDataset, iRadius, options);
    NormYSlicesAndScaleRespByWidth(H2D_jetPtResponseMatrix, iDataset, iRadius, options);
    MergeResponseMatrixBins(H2D_jetPtResponseMatrix, iDataset, iRadius, options);
  }
}

void ReweightResponseMatrixWithPrior_fineBinningOnly(TH2D* &H2D_jetPtResponseMatrix, int iDataset, int iRadius, std::string options) {
  //not 
  // before this, all y-slices (ie pt gen slices) have been normalised to 1;means each pt gen slice has a proba of 1
  // withthis function, we give each slice a weight so that they have different normalisation values, corresponding to the prior 

  // prior choice; none by default (flat)
  TH1D* priorSpectrumWeighting;
  if (options.find("mcpPriorUnfolding") != std::string::npos) {
    if (!normGenAndMeasByNEvtsBeforeUnfolding) {
      Get_Pt_spectrum_mcp_fineBinning_preWidthScalingAtEndAndEvtNorm(priorSpectrumWeighting, iDataset, iRadius, false, options); 
    } else {
      Get_Pt_spectrum_mcp_fineBinning_preWidthScalingAtEnd(priorSpectrumWeighting, iDataset, iRadius, false, options); 
    }
    WeightMatrixWithPrior(H2D_jetPtResponseMatrix, priorSpectrumWeighting, doUnfoldingPriorDivision);
  }
  if (options.find("mcdPriorUnfolding") != std::string::npos) {
    if (!normGenAndMeasByNEvtsBeforeUnfolding) {
      Get_Pt_spectrum_mcd_fineBinning_preWidthScalingAtEndAndEvtNorm(priorSpectrumWeighting, iDataset, iRadius, options); 
    } else {
      Get_Pt_spectrum_mcd_fineBinning_preWidthScalingAtEnd(priorSpectrumWeighting, iDataset, iRadius, options); 
    }
    WeightMatrixWithPrior(H2D_jetPtResponseMatrix, priorSpectrumWeighting, doUnfoldingPriorDivision);
  }
  if (options.find("measuredPriorUnfolding") != std::string::npos) {
    if (!normGenAndMeasByNEvtsBeforeUnfolding) {
      Get_Pt_spectrum_bkgCorrected_fineBinning_preWidthScalingAtEndAndEvtNorm(priorSpectrumWeighting, iDataset, iRadius, options); 
    } else {
      Get_Pt_spectrum_bkgCorrected_fineBinning_preWidthScalingAtEnd(priorSpectrumWeighting, iDataset, iRadius, options); 
    }
    WeightMatrixWithPrior(H2D_jetPtResponseMatrix, priorSpectrumWeighting, doUnfoldingPriorDivision);
  }
  if (options.find("testAliPhysics") != std::string::npos) {
    if (!normGenAndMeasByNEvtsBeforeUnfolding) {
      Get_Pt_spectrum_bkgCorrected_fineBinning_preWidthScalingAtEndAndEvtNorm(priorSpectrumWeighting, iDataset, iRadius, options); 
    } else {
      Get_Pt_spectrum_bkgCorrected_fineBinning_preWidthScalingAtEnd(priorSpectrumWeighting, iDataset, iRadius, options); 
    }
    H2D_jetPtResponseMatrix = (TH2D*)NormalizeResponsMatrixYaxisWithPrior(H2D_jetPtResponseMatrix, priorSpectrumWeighting)->Clone(H2D_jetPtResponseMatrix->GetName()+(TString)"_testAliPhysics");
  }
}

void Get_thetagMatrix_detectorResponse4D(TH2D* &H2D_thetagMatrix_detectorResponse,
                                       int iDataset,
                                       int iRadius,
                                       double ptmin,
                                       double ptmax)
{
  // =============================
  // === 0. Chargement du histo ===
  // =============================
  TString partialUniqueSpecifier = Datasets[iDataset] + "_R=" + Form("%.1f", arrayRadius[iRadius]);
  TString histoname = "h4_ptMCD_ptMCP_thetagMCD_thetagMCP_norange_eventwise";

  cout << "📂 Lecture du THnSparse 4D : " << histoname << endl;

  THnSparse* h4 = (THnSparse*) file_O2Analysis_MCfileForMatrix[iDataset]->Get(analysisWorkflowMC + "/" + histoname);
  if (!h4) {
    cerr << "❌ ERREUR : impossible de trouver " << histoname << endl;
    return;
  }

  cout << "✅ Histogramme chargé : " << histoname << endl;
  cout << "   - Entries : " << h4->GetEntries() << endl;
  cout << "   - Dimensions : " << h4->GetNdimensions() << endl;

  // =============================
  // === 1. Définition des axes ===
  // =============================
  // Vérifie que l’ordre correspond bien :
  // 0 : ptMCD  (detector level)
  // 1 : ptMCP  (particle level)
  // 2 : thetagMCD
  // 3 : thetagMCP
  int axisPtMCD = 0;
  int axisPtMCP = 1;
  int axisThetagMCD = 2;
  int axisThetagMCP = 3;

  // =============================
  // === 2. Sélection en pT (MCP) ===
  // =============================
  cout << Form("🔧 Application du cut en pT^MCP : %.1f < pT < %.1f GeV/c", ptmin, ptmax) << endl;
  h4->GetAxis(axisPtMCP)->SetRangeUser(ptmin, ptmax);

  // =============================
  // === 3. Projection en 2D ===
  // =============================
  // Projection : X = thetagMCD, Y = thetagMCP
  TH2D* h2 = (TH2D*) h4->Projection(axisThetagMCP, axisThetagMCD);
  h2->SetName(Form("H2D_thetagMatrix_detectorResponse_R%.1f_%s_pt%.0f_%.0f",
                   arrayRadius[iRadius],
                   Datasets[iDataset].Data(),
                   ptmin, ptmax));

  h2->SetTitle(Form("#theta_{g}^{MCD} vs #theta_{g}^{MCP}  (%.0f < p_{T}^{MCP} < %.0f GeV/c)",
                    ptmin, ptmax));
  h2->GetXaxis()->SetTitle("#theta_{g}^{MCD}");
  h2->GetYaxis()->SetTitle("#theta_{g}^{MCP}");

  // =============================
  // === 4. Normalisation éventuelle ===
  // =============================
  if (doYSliceNormToOneDetResp) {
    cout << "🧮 Normalisation slice par slice (Y)" << endl;
    NormaliseYSlicesToOne(h2);
  }

  if (normDetRespByNEvts) {
    double Nevents = GetNEventsSelected_JetFramework(file_O2Analysis_MCfileForMatrix[iDataset], analysisWorkflowMC);
    cout << "🧮 Normalisation par N_events = " << Nevents << endl;
    if (Nevents > 0) h2->Scale(1. / Nevents);
  }

  // =============================
  // === 5. Sortie finale ===
  // =============================
  H2D_thetagMatrix_detectorResponse = (TH2D*) h2->Clone("H2D_thetagMatrix_detectorResponse" + partialUniqueSpecifier);

  cout << "✅ Matrice thetag (2D) prête : " << H2D_thetagMatrix_detectorResponse->GetName() << endl;
}


void Get_PtResponseMatrix_Fluctuations(TH2D* &H2D_jetPtResponseMatrix_fluctuations, int iDataset, int iRadius) { 
  // see Hiroki Yokoyama thesis
  // iRadius is for chosing the pT binning
  TString partialUniqueSpecifier = Datasets[iDataset]+"_R="+Form("%.1f",arrayRadius[iRadius]);

  if (useFactorisedMatrix == false){
    TH2D H2D_identity = TH2D("H2D_response_"+partialUniqueSpecifier, "H2D_response_"+partialUniqueSpecifier, nBinPtJetsFine[iRadius], ptBinsJetsFine[iRadius], nBinPtJetsFine[iRadius], ptBinsJetsFine[iRadius]);
    for(int iBinRec = 0; iBinRec <= H2D_identity.GetNbinsX()+1; iBinRec++){
      H2D_identity.SetBinContent(iBinRec, iBinRec, 1);
      H2D_identity.SetBinError(iBinRec, iBinRec, 0);
    }
    H2D_jetPtResponseMatrix_fluctuations = (TH2D*)H2D_identity.Clone("H2D_jetPtResponseMatrix_fluctuations"+partialUniqueSpecifier);
  } else {
    cout << "I should check that the average of each ptGen slice is as displaced to the diagonal as the randomCone distrib is; ie should I use GetBinLowEdge or GetBinLowEdge+width for ptGen" << endl;


    TH2D* H2D_fluctuations_centrality;
    TH1D* H1D_fluctuations;

    H2D_fluctuations_centrality = (TH2D*)((TH2D*)file_O2Analysis_list[iDataset]->Get(analysisWorkflowBkg+"/h2_centrality_rhorandomcone"+randomConeTypeList[randomConeType]))->Clone("Get_PtResponseMatrix_Fluctuations"+Datasets[iDataset]+DatasetsNames[iDataset]);
    H2D_fluctuations_centrality->Sumw2();


    int ibinCent_low = H2D_fluctuations_centrality->GetXaxis()->FindBin(centralityRange[0]);
    int ibinCent_high = H2D_fluctuations_centrality->GetXaxis()->FindBin(centralityRange[1]);
    H1D_fluctuations = (TH1D*)H2D_fluctuations_centrality->ProjectionY("bkgFluctuationCentrality_highRes_"+partialUniqueSpecifier, ibinCent_low, ibinCent_high, "e");

    NormaliseRawHistToIntegral(H1D_fluctuations); // normalising fluctuations to 1
    // cout << "Integral H1D_fluctuations: " << H1D_fluctuations->Integral(1, H1D_fluctuations->GetNbinsX()) << endl;

    TH2D H2D_response = TH2D("H2D_response_"+partialUniqueSpecifier, "H2D_response_"+partialUniqueSpecifier, nBinPtJetsFine[iRadius], ptBinsJetsFine[iRadius], nBinPtJetsFine[iRadius], ptBinsJetsFine[iRadius]); 

    //==================== Build response matrix: shift deltaPt by pT gen along the pT rec axis ====================//
    int ibinZeroFluct= H1D_fluctuations->FindBin(0+GLOBAL_epsilon);
    double integralError;
    double ptGen, ptRec_low, ptRec_up;
    int iBin_fluct_low, iBin_fluct_high;
    for(int iBinRec = 0; iBinRec <= H2D_response.GetNbinsX()+1; iBinRec++){
      for(int iBinGen = 0; iBinGen <= H2D_response.GetNbinsY()+1; iBinGen++){
        ptGen = H2D_response.GetYaxis()->GetBinCenter(iBinGen); // was bincenter before but then it'd give .5 values of GeV, and 
        ptRec_low = H2D_response.GetXaxis()->GetBinLowEdge(iBinRec);
        ptRec_up = H2D_response.GetXaxis()->GetBinLowEdge(iBinRec+1);
        // int iBin_fluct_low = H1D_fluctuations->GetXaxis()->FindBin(ptRec_low - ptGen + GLOBAL_epsilon);
        // int iBin_fluct_high = H1D_fluctuations->GetXaxis()->FindBin(ptRec_up - ptGen - GLOBAL_epsilon);
        iBin_fluct_low = H1D_fluctuations->GetXaxis()->FindBin(ptRec_low - ptGen + GLOBAL_epsilon);
        if (iBinGen == 10 && (iBin_fluct_low >= iBin_fluct_high)) { // checks iBinRec =10 so that the message doesn't appear NbinsX*NBinsY times
          cout << "Get_PtResponseMatrix_Fluctuations: some bins are counted twice in the integral, binning needs to be looked at, right now the fluctuation matrix has too many entries" << endl;
        }
        iBin_fluct_high = H1D_fluctuations->GetXaxis()->FindBin(ptRec_up - ptGen - GLOBAL_epsilon)-1;
        H2D_response.SetBinContent(iBinRec, iBinGen, H1D_fluctuations->IntegralAndError(iBin_fluct_low, iBin_fluct_high, integralError)); 
        H2D_response.SetBinError(iBinRec, iBinGen, integralError); 
        // if (iBinRec == 10) {
        //   cout << "DeltaBin = " << iBin_fluct_high - iBin_fluct_low << ", iBin_low = " << iBin_fluct_low << ", iBin_high = " << iBin_fluct_high << endl;
        // }
      }
    }
    float integralCheck = H2D_response.Integral(1, H2D_response.GetNbinsX(), 10, 10);
    if (!(integralCheck > 0.999)) {
      cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl;
      cout << "!!!!!!Response matrix of fluctuations does not have line integral equal to 1; something is wrong with binning!!!!!!" << endl;
      cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl;
    }

    //========================================= Build response matrix end =========================================//

    H2D_jetPtResponseMatrix_fluctuations = (TH2D*)H2D_response.Clone("H2D_jetPtResponseMatrix_fluctuations"+partialUniqueSpecifier);
  }
}

#endif