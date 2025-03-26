#ifndef JETSPECTRUM_UNFOLDING_H
#define JETSPECTRUM_UNFOLDING_H


int GetSvdBestRegularisationParameter_notYetSatisfying(TSVDUnfold* unfoldTSvd);
std::pair<int, RooUnfold*> Get_Pt_spectrum_unfolded_preWidthScalingAtEndAndEvtNorm(TH1D* &H1D_jetPt_unfolded, TH1D* &measuredInput, int iDataset, int iRadius, int unfoldParameterInput, std::string options);
std::pair<int, RooUnfold*> Get_Pt_spectrum_unfolded_preWidthScalingAtEnd(TH1D* &H1D_jetPt_unfolded, TH1D* &measuredInput, int iDataset, int iRadius, int unfoldParameterInput, std::string options);
std::pair<int, RooUnfold*> Get_Pt_spectrum_unfolded(TH1D* &H1D_jetPt_unfolded, TH1D* &measuredInput, int iDataset, int iRadius, int unfoldParameterInput, std::string options);


void Get_Pt_spectrum_dataUnfoldedThenRefolded_preWidthScalingAtEndAndEvtNorm(TH1D* &H1D_jetPt_unfoldedThenRefolded, int iDataset, int iRadius, int unfoldParameterInput, std::string options);
void Get_Pt_spectrum_dataUnfoldedThenRefolded_preWidthScalingAtEnd(TH1D* &H1D_jetPt_unfoldedThenRefolded, int iDataset, int iRadius, int unfoldParameterInput, std::string options);
void Get_Pt_spectrum_dataUnfoldedThenRefolded(TH1D* &H1D_jetPt_unfoldedThenRefolded, int iDataset, int iRadius, int unfoldParameterInput, std::string options);
void Get_Pt_spectrum_dataUnfoldedThenRefolded_RooUnfoldMethod_preWidthScalingAtEndAndEvtNorm(TH1D* &H1D_jetPt_unfoldedThenRefolded, int iDataset, int iRadius, int unfoldParameterInput, std::string options);
void Get_Pt_spectrum_dataUnfoldedThenRefolded_RooUnfoldMethod_preWidthScalingAtEnd(TH1D* &H1D_jetPt_unfoldedThenRefolded, int iDataset, int iRadius, int unfoldParameterInput, std::string options);
void Get_Pt_spectrum_dataUnfoldedThenRefolded_RooUnfoldMethod(TH1D* &H1D_jetPt_unfoldedThenRefolded, int iDataset, int iRadius, int unfoldParameterInput, std::string options);


#endif