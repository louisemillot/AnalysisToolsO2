#ifndef JETSPECTRU_RESPONSEMATRIXFUNCTIONS_H
#define JETSPECTRU_RESPONSEMATRIXFUNCTIONS_H

void Get_PtResponseMatrix_Fluctuations(TH2D* &H2D_jetPtResponseMatrix_fluctuations, int iDataset, int iRadius);
void Get_PtResponseMatrix_detectorResponse(TH2D* &H2D_jetPtResponseMatrix_detectorResponse, int iDataset, int iRadius);
void Get_PtResponseMatrix_DetectorAndFluctuationsCombined(TH2D* &H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined, TH2D* H2D_jetPtResponseMatrix_detectorResponse, TH2D* H2D_jetPtResponseMatrix_fluctuations, int iDataset, int iRadius, std::string options);
void Get_PtResponseMatrix_DetectorAndFluctuationsCombined_fineBinning(TH2D* &H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined_fineBinning, TH2D* H2D_jetPtResponseMatrix_detectorResponse, TH2D* H2D_jetPtResponseMatrix_fluctuations, int iDataset, int iRadius, std::string options);
void ReweightResponseMatrixWithPrior(TH2D* &H2D_jetPtResponseMatrix, int iDataset, int iRadius, std::string options);
void ReweightResponseMatrixWithPrior_fineBinning(TH2D* &H2D_jetPtResponseMatrix, int iDataset, int iRadius, std::string options);
void NormYSlicesAndScaleRespByWidth(TH2D* &H2D_jetPtResponseMatrix, int iDataset, int iRadius, std::string options);
void MergeResponseMatrixBins(TH2D* &H2D_jetPtResponseMatrix, int iDataset, int iRadius, std::string options);
void FinaliseResponseMatrix(TH2D* &H2D_jetPtResponseMatrix, int iDataset, int iRadius, std::string options);


#endif