#ifndef JETSPECTRUM_SPECTRAGETTERS_H
#define JETSPECTRUM_SPECTRAGETTERS_H

void Get_Pt_spectrum_bkgCorrected_recBinning_preWidthScalingAndEvtNorm(TH1D* &H1D_jetPt_raw, int iDataset, int iRadius, std::string options);
void Get_Pt_spectrum_bkgCorrected_recBinning_preWidthScaling(TH1D* &H1D_jetPt, int iDataset, int iRadius, std::string options);
void Get_Pt_spectrum_bkgCorrected_recBinning(TH1D* &H1D_jetPt_raw, int iDataset, int iRadius, std::string options);

void Get_Pt_spectrum_bkgCorrected_genBinning_preWidthScalingAndEvtNorm(TH1D* &H1D_jetPt_raw, int iDataset, int iRadius, std::string options);
void Get_Pt_spectrum_bkgCorrected_genBinning_preWidthScaling(TH1D* &H1D_jetPt, int iDataset, int iRadius, std::string options);
void Get_Pt_spectrum_bkgCorrected_genBinning(TH1D* &H1D_jetPt_raw, int iDataset, int iRadius, std::string options);

void Get_Pt_spectrum_bkgCorrected_fineBinning_preWidthScalingAndEvtNorm(TH1D* &H1D_jetPt_raw, int iDataset, int iRadius, std::string options);
void Get_Pt_spectrum_bkgCorrected_fineBinning_preWidthScaling(TH1D* &H1D_jetPt, int iDataset, int iRadius, std::string options);
void Get_Pt_spectrum_bkgCorrected_fineBinning(TH1D* &H1D_jetPt_raw, int iDataset, int iRadius, std::string options);

void Get_Pt_spectrum_mcp_genBinning_preWidthScalingAndEvtNorm(TH1D* &H1D_jetPt_mcp, int iDataset, int iRadius, bool fcontrolMC, std::string options);
void Get_Pt_spectrum_mcp_genBinning_preWidthScaling(TH1D* &H1D_jetPt_mcp, int iDataset, int iRadius, bool fcontrolMC, std::string options);
void Get_Pt_spectrum_mcp_genBinning(TH1D* &H1D_jetPt_mcp, int iDataset, int iRadius, bool fcontrolMC, std::string options);

void Get_Pt_spectrum_mcp_recBinning_preWidthScalingAndEvtNorm(TH1D* &H1D_jetPt_mcp, int iDataset, int iRadius, bool fcontrolMC, std::string options);
void Get_Pt_spectrum_mcp_recBinning_preWidthScaling(TH1D* &H1D_jetPt_mcp, int iDataset, int iRadius, bool fcontrolMC, std::string options);
void Get_Pt_spectrum_mcp_recBinning(TH1D* &H1D_jetPt_mcp, int iDataset, int iRadius, bool fcontrolMC, std::string options);

void Get_Pt_spectrum_mcdMatched_genBinning_preWidthScalingAndEvtNorm(TH1D* &H1D_jetPt_mcd, int iDataset, int iRadius, std::string options);
void Get_Pt_spectrum_mcdMatched_genBinning_preWidthScaling(TH1D* &H1D_jetPt_mcdMatched, int iDataset, int iRadius, std::string options);
void Get_Pt_spectrum_mcdMatched_genBinning(TH1D* &H1D_jetPt_mcd, int iDataset, int iRadius, std::string options);

void Get_Pt_spectrum_mcpMatched_genBinning_preWidthScalingAndEvtNorm(TH1D* &H1D_jetPt_mcpMatched, int iDataset, int iRadius, std::string options);
void Get_Pt_spectrum_mcpMatched_genBinning_preWidthScaling(TH1D* &H1D_jetPt_mcpMatched, int iDataset, int iRadius, std::string options);
void Get_Pt_spectrum_mcpMatched_genBinning(TH1D* &H1D_jetPt_mcp, int iDataset, int iRadius, std::string options);

void Get_Pt_spectrum_mcp_fineBinning_preWidthScalingAndEvtNorm(TH1D* &H1D_jetPt_mcp, int iDataset, int iRadius, bool fcontrolMC, std::string options);
void Get_Pt_spectrum_mcp_fineBinning_preWidthScaling(TH1D* &H1D_jetPt_mcp, int iDataset, int iRadius, bool fcontrolMC, std::string options);
void Get_Pt_spectrum_mcp_fineBinning(TH1D* &H1D_jetPt_mcp, int iDataset, int iRadius, bool fcontrolMC, std::string options);

void Get_Pt_spectrum_mcd_fineBinning_preWidthScalingAndEvtNorm(TH1D* &H1D_jetPt_mcp, int iDataset, int iRadius, std::string options);
void Get_Pt_spectrum_mcd_fineBinning_preWidthScaling(TH1D* &H1D_jetPt_mcd, int iDataset, int iRadius, std::string options);
void Get_Pt_spectrum_mcd_fineBinning(TH1D* &H1D_jetPt_mcp, int iDataset, int iRadius, std::string options);

void Get_Pt_spectrum_mcd_recBinning_preWidthScalingAndEvtNorm(TH1D* &H1D_jetPt_mcp, int iDataset, int iRadius, std::string options);
void Get_Pt_spectrum_mcd_recBinning_preWidthScaling(TH1D* &H1D_jetPt_mcd, int iDataset, int iRadius, std::string options);
void Get_Pt_spectrum_mcd_recBinning(TH1D* &H1D_jetPt_mcp, int iDataset, int iRadius, std::string options);

void Get_Pt_spectrum_mcpMatched_fineBinning(TH1D* &H1D_jetPt_mcpMatched, int iDataset, int iRadius, std::string options);
void Get_Pt_spectrum_mcdMatched_recBinning(TH1D* &H1D_jetPt_mcdMatched, int iDataset, int iRadius, std::string options);
void Get_Pt_spectrum_mcdMatched_fineBinning(TH1D* &H1D_jetPt_mcdMatched, int iDataset, int iRadius, std::string options);


#endif