#include "TStyle.h"
#include "TFile.h"
#include "TFitResult.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TF1.h"
#include "TLatex.h"
#include "TLine.h"
#include "TRatioPlot.h"
#include "TLegend.h"
#include "TH1.h"
// #include <RooUnfold.h>
// #include "RooUnfoldResponse.h"
// #include "RooUnfoldBayes.h"
// #include "RooUnfoldBinByBin.h"

//My Libraries
#include "./JetSpectrum_settings.h"
#include "../Settings/AxisTitles.h"
#include "../Utilities/AnalysisUtilities.h"
#include "../Utilities/HistogramUtilities.h"
#include "../Utilities/AnalysisUtilities.C" // bizarre but only including the .h fils doesn't work for the standard 'root macro.C+' method, I need to include the .C as well
#include "../Utilities/HistogramUtilities.C" // bizarre but only including the .h fils doesn't work for the standard 'root macro.C+' method, I need to include the .C as well

#include<array>
#include <iomanip>
#include <sstream>
#include <string.h>
using namespace std;

// Misc utilities
void SetStyle(Bool_t graypalette=kFALSE);
void LoadLibs();

// Plot Utilities
TString contextDataset1D(int iDataset, float* variableRange, const char options[]);
TString contextDatasetComp(float jetRadius, float* variableRange, const char options[]);
TString contextPtRange(float* PtRange);
TString contextEtaRange(float* PtRange);
TString contextJetRadius(float jetRadius);

//////////// Pt Spectrum analysis functions
void Get_Pt_spectrum_raw(TH1D** H1D_jetPt_raw, int iDataset, float* etaRange, const char options[]);
void Draw_Pt_spectrum_raw(int iDataset, float* etaRange);

/////////////////////////////////////////////////////
///////////////////// Main Macro ////////////////////
/////////////////////////////////////////////////////

void JetSpectrum() {
  // Load necessary libraries
  LoadLibs();
  // Set the default style
  SetStyle();

  // TString* SaveAs_Title = new TString("");
  TString* texXtitle = new TString("");
  TString* texYtitle = new TString("");
  // TString* Extra = new TString("");

  int iDataset = 0;

  float etaRangeSym[2] = {-0.5, 0.5};
  float etaRangeNeg[2] = {-0.5, 0};
  float etaRangePos[2] = {0, 0.5};

  float jetRadiusForDataComp = 0.4;
  float jetR02 = 0.2;
  float jetR06 = 0.6;
  
  float ptRange[2] = {1, 200};

  Draw_Pt_spectrum_raw(iDataset, etaRangeSym);


}

/////////////////////////////////////////////////////
/////////////////// Misc utilities //////////////////
/////////////////////////////////////////////////////

void LoadLibs() {
  // gSystem->Load("libCore.so");  
  // gSystem->Load("libGeom.so");
  // gSystem->Load("libPhysics.so");
  // gSystem->Load("libVMC");
  // gSystem->Load("libTree");
  // gSystem->Load("libMinuit");
  // gSystem->Load("libSTEERBase");
  // gSystem->Load("libESD");
  // gSystem->Load("libAOD");
  // gSystem->Load("libANALYSIS");
  // gSystem->Load("libANALYSISalice");
  // gSystem->Load("libCORRFW");
  // gSystem->Load("libPWGTools");
}

void SetStyle(Bool_t graypalette) {
  cout << "Setting style!" << endl;
  
  gStyle->Reset("Plain");
  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);
  if(graypalette) gStyle->SetPalette(8,0);
  else gStyle->SetPalette(1);
  gStyle->SetCanvasColor(10);
  gStyle->SetCanvasBorderMode(0);
  gStyle->SetFrameLineWidth(1);
  gStyle->SetFrameFillColor(kWhite);
  gStyle->SetPadColor(10);
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);
  gStyle->SetPadBottomMargin(0.15);
  gStyle->SetPadLeftMargin(0.15);
  gStyle->SetHistLineWidth(1);
  gStyle->SetHistLineColor(kRed);
  gStyle->SetFuncWidth(2);
  gStyle->SetFuncColor(kGreen);
  gStyle->SetLineWidth(2);
  gStyle->SetLineScalePS(1);
  gStyle->SetLabelSize(0.045,"xyz");
  gStyle->SetLabelOffset(0.01,"y");
  gStyle->SetLabelOffset(0.01,"x");
  gStyle->SetLabelColor(kBlack,"xyz");
  gStyle->SetTitleSize(0.05,"xyz");
  gStyle->SetTitleOffset(1.4,"y");
  gStyle->SetTitleOffset(1.2,"x");
  gStyle->SetTitleFillColor(kWhite);
  gStyle->SetTextSizePixels(26);
  gStyle->SetTextFont(42);
  //  gStyle->SetTickLength(0.04,"X");  gStyle->SetTickLength(0.04,"Y"); 

  gStyle->SetLegendBorderSize(0);
  gStyle->SetLegendFillColor(kWhite);
  //  gStyle->SetFillColor(kWhite);
  gStyle->SetLegendFont(42);
}




//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////// Context Utilities /////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

TString contextDatasetRadiusComp(int iDataset, float* variableRange, const char options[]){
  TString texContextDatasetRadiusComp;
  if (strstr(options, "pt") != NULL) { //  || strstr(options, "ratio") != NULL not sure why I had this here
    texContextDatasetRadiusComp = "#splitline{"+*texDatasetsComparisonCommonDenominator+" "+DatasetsNames[iDataset]+"}{#splitline{2023 QC}{"+contextPtRange(variableRange)+"}}";
  }
  if (strstr(options, "eta") != NULL) { //  || strstr(options, "ratio") != NULL not sure why I had this here
    texContextDatasetRadiusComp = "#splitline{"+*texDatasetsComparisonCommonDenominator+" "+DatasetsNames[iDataset]+"}{#splitline{2023 QC}{"+contextEtaRange(variableRange)+"}}";
  }

  return texContextDatasetRadiusComp;
}

TString contextDatasetComp(float jetRadius, float* variableRange, const char options[]){
  TString texcontextDatasetComp;
  if (strstr(options, "pt") != NULL) { //  || strstr(options, "ratio") != NULL not sure why I had this here
    texcontextDatasetComp = "#splitline{"+*texDatasetsComparisonCommonDenominator+"}{#splitline{"+contextJetRadius(jetRadius)+"}{"+contextPtRange(variableRange)+"}}";
  }
  if (strstr(options, "eta") != NULL) { //  || strstr(options, "ratio") != NULL not sure why I had this here
    texcontextDatasetComp = "#splitline{"+*texDatasetsComparisonCommonDenominator+"}{#splitline{"+contextJetRadius(jetRadius)+"}{"+contextEtaRange(variableRange)+"}}";
  }

  return texcontextDatasetComp;
}

TString contextPtRange(float* PtRange){
  std::stringstream ss;
  ss << PtRange[0] << " < #it{p}_{T} < " << PtRange[1];
  TString textContext((TString)ss.str());
  // TString texDataset(Form("%.0f", PtRange[0])+" < #it{p}_{T} < "+Form("%.0f", PtRange[1]));
  return textContext;
}

TString contextEtaRange(float* EtaRange){
  std::stringstream ss;
  ss << EtaRange[0] << " < #eta < " << EtaRange[1];
  TString textContext((TString)ss.str());
  return textContext;
}

TString contextJetRadius(float jetRadius){
  std::stringstream ss;
  ss << " R = " << jetRadius;
  TString textContext((TString)ss.str());
  // TString texDataset(Form("%.0f", PtRange[0])+" < #it{p}_{T} < "+Form("%.0f", PtRange[1]));
  return textContext;
}


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////// Spectrum  plot functions ///////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



void Get_Pt_spectrum_raw(TH1D** H1D_jetPt_raw, int iDataset, float* etaRange, const char options[]) {
  TH3D* H3D_jetRjetPtjetEta;
  TH1D* H1D_jetPt[nRadius];
  // TH1D* H1D_jetPt_raw[nRadius];

  H3D_jetRjetPtjetEta = (TH3D*)((TH3D*)file_O2Analysis_list[iDataset]->Get(analysisWorkflow+"/h3_jet_r_jet_pt_jet_eta"))->Clone("Draw_Pt_RadiusComparison"+Datasets[iDataset]+Form("%.1f", etaRange[0])+"<eta<"+Form("%.1f", etaRange[1]));
  H3D_jetRjetPtjetEta->Sumw2();

  int binEtaEdges[2] = {H3D_jetRjetPtjetEta->GetZaxis()->FindBin(etaRange[0]), H3D_jetRjetPtjetEta->GetZaxis()->FindBin(etaRange[1])};
  if (binEtaEdges[0] == 0) 
    cout << "WARNING: Get_Pt_spectrum_raw is counting the underflow with the chosen etaRange" << endl;
  if (binEtaEdges[1] == H3D_jetRjetPtjetEta->GetZaxis()->GetNbins()+1) 
    cout << "WARNING: Get_Pt_spectrum_raw is counting the overflow with the chosen etaRange" << endl;
  int ibinJetRadius = 0;

  for(int iRadius = 0; iRadius < nRadius; iRadius++){
    ibinJetRadius = H3D_jetRjetPtjetEta->GetXaxis()->FindBin(arrayRadius[iRadius]+e_binEdge);

    H1D_jetPt[iRadius] = (TH1D*)H3D_jetRjetPtjetEta->ProjectionY("jetPt_"+RadiusLegend[iRadius]+Form("%.1f", etaRange[0])+"<eta<"+Form("%.1f", etaRange[1]), ibinJetRadius, ibinJetRadius, binEtaEdges[0], binEtaEdges[1]);
    H1D_jetPt_raw[iRadius] = (TH1D*)H1D_jetPt[iRadius]->Rebin(PtRebinValue,"jetPt_rebinned_"+RadiusLegend[iRadius]+Form("%.1f", etaRange[0])+"<eta<"+Form("%.1f", etaRange[1]));

    if (strstr(options, "evtNorm") != NULL) {
      NormaliseYieldToNEvents(H1D_jetPt_raw[iRadius], GetNEventsSel8(file_O2Analysis_list[iDataset]));
    }
    if (strstr(options, "entriesNorm") != NULL) {
      NormaliseYieldToNEntries(H1D_jetPt_raw[iRadius]);
    }

  }
}

void Draw_Pt_spectrum_raw(int iDataset, float* etaRange) {

  TH1D* H1D_jetPt_raw[nRadius];
  Get_Pt_spectrum_raw(H1D_jetPt_raw, iDataset, etaRange, "evtNorm");

  TString* pdfName = new TString("jet_"+jetType[iJetType]+"_"+jetLevel[iJetLevel]+"_"+Datasets[iDataset]+"_Pt_@eta["+Form("%.1f", etaRange[0])+","+Form("%.1f", etaRange[1])+"]");

  TString textContext(contextDatasetRadiusComp(iDataset, etaRange, "eta"));
  
  Draw_TH1_Histograms_in_one(H1D_jetPt_raw, RadiusLegend, nRadius, textContext, pdfName, texPtX, texJetNormPtYield, texCollisionDataInfo, "logy");
}

// // PhD fuction for V0s
// void Get_V0_RawYield_vsPt_FeeddownCorrected_withUnfolding(TH1D* &hRawYield_vsPt, TH3D* H3D_detectedV0s_DataLike, TH1I* H1I_SelectedEventCount, TFile *file_O2Analysis, double* SignalMean, double* SignalStandardDev, int ipart, double* pTbins, int nbinpT, int ibinXaxisCut_low,  int ibinXaxisCut_high, int warning_cutArry_ID, double SideBandSizeMultiplierModifier, int SignalExtractionType, int InvMassRebinFactor) {
//   float SelectedEventCount = H1I_SelectedEventCount->GetEntries();

//   TH1D* truth = (TH1D*)file_O2Analysis->Get("lambdakzero-particle-count-mc/h"+NamePart[ipart]+"Count_PtDiff");
//   TH1D* truth_rebinned = (TH1D*)truth->Rebin(nbinpT,"truth_rebinned",pTbins);
  
//   Get_RawYield_vsPt_FeeddownCorrected_RawCount(hRawYield_vsPt, H3D_detectedV0s_DataLike, H1I_SelectedEventCount, file_O2Analysis, SignalMean, SignalStandardDev, ipart, pTbins, nbinpT, ibinXaxisCut_low, ibinXaxisCut_high, warning_cutArry_ID, SideBandSizeMultiplierModifier, SignalExtractionType, InvMassRebinFactor);
//   TH2D *H2D_PtResponseMatrix = new TH2D("H2D_PtResponseMatrix", "H2D_PtResponseMatrix", nbinpT, pTbins, nbinpT, pTbins);
//   Get_PtResponseMatrix(H2D_PtResponseMatrix, file_O2Analysis, pTbins, nbinpT);

//   // Get_RawYield_vsPt_FeeddownCorrected(hRawYield_vsPt, H3D_detectedV0s_DataLike, H1I_SelectedEventCount, file_O2Analysis, SignalMean, SignalStandardDev, ipart, pTbins, nbinpT, ibinXaxisCut_low, ibinXaxisCut_high, warning_cutArry_ID, SideBandSizeMultiplierModifier, SignalExtractionType, InvMassRebinFactor);
//   // TH2D *H2D_PtResponseMatrix_Density = new TH2D("H2D_PtResponseMatrix", "H2D_PtResponseMatrix", nbinpT, pTbins, nbinpT, pTbins);
//   // Get_PtResponseMatrix_Density(H2D_PtResponseMatrix, file_O2Analysis, pTbins, nbinpT);
//   // for(int ibinPt = 1; ibinPt <= nbinpT; ibinPt++){
//   //   double dpT = truth_rebinned->GetXaxis()->GetBinWidth(ibinPt);
//   //   double drapidity = 1.5; //1.5
//   //   double d2N_dpTdy = truth_rebinned->GetBinContent(ibinPt) *1./dpT*1./drapidity;
//   //   truth_rebinned->SetBinContent(ibinPt,d2N_dpTdy);
//   //   truth_rebinned->SetBinError(ibinPt,truth_rebinned->GetBinError(ibinPt) *1./dpT*1./drapidity);  // error on d2N_dpTdy 
//   // }    
//   // truth_rebinned->Scale(1./SelectedEventCount);

//   TH1D* measured = (TH1D*)hRawYield_vsPt->Clone("measured");

//   // TH1D* projX = (TH1D*)H2D_PtResponseMatrix->ProjectionX("projX",0,-1);
//   // TH1D* projY = (TH1D*)H2D_PtResponseMatrix->ProjectionY("projY",0,-1);
//   cout << "--------------------------------------------------------------------- PRE UNFOLDING TEST" << endl;
//   for(int ibinPt = 1; ibinPt <= nbinpT; ibinPt++){
//     cout << "measured(" << ibinPt << ") = " << measured->GetBinContent(ibinPt) << endl;
//     cout << "truth_rebinned(" << ibinPt << ") = " << truth_rebinned->GetBinContent(ibinPt) << endl;
//     for(int jbinPt = 1; jbinPt <= nbinpT; jbinPt++){
//     cout << "H2D_PtResponseMatrix(" << ibinPt << "," << jbinPt << ") = " << H2D_PtResponseMatrix->GetBinContent(ibinPt,jbinPt) << endl;
//     }
//   }

//   RooUnfoldResponse* response = new RooUnfoldResponse(measured, truth_rebinned, H2D_PtResponseMatrix);
//   // RooUnfoldResponse* response = new RooUnfoldResponse(projX, projY, H2D_PtResponseMatrix, "", "");
//   // RooUnfold* unfold = new RooUnfoldBayes(response, measured, 4);
//   RooUnfold* unfold = new RooUnfoldBinByBin(response, measured);

//   // TH1D* hist_unfold = (TH1D*)unfold.Hunfold();
//   // TH1D* hist_unfold = (TH1D*)unfold->Hreco()->Clone("hist_unfold");
//   TH1D* hist_unfold = static_cast<TH1D*>(unfold->Hreco());


//   hRawYield_vsPt = (TH1D*)hist_unfold->Clone("hRawYield_vsPt");

//   for(int ibinPt = 1; ibinPt <= nbinpT; ibinPt++){
//     double dpT = hRawYield_vsPt->GetXaxis()->GetBinWidth(ibinPt);
//     double drapidity = 1.5; //1.5
//     double d2N_dpTdy = hRawYield_vsPt->GetBinContent(ibinPt) *1./dpT*1./drapidity;
//     hRawYield_vsPt->SetBinContent(ibinPt,d2N_dpTdy);
//     hRawYield_vsPt->SetBinError(ibinPt,hRawYield_vsPt->GetBinError(ibinPt) *1./dpT*1./drapidity);  // error on d2N_dpTdy 
//   }    

//   hRawYield_vsPt->Scale(1./SelectedEventCount);

//   // Lambda and AntiLambda Feeddown
//   TH1D *H1D_Feeddown_Correction = new TH1D("H1D_Feeddown_Correction", "H1D_Feeddown_Correction", nbinpT, pTbins);
//   if (ipart == 1 || ipart == 2) {
//     Get_FeedDown(H1D_Feeddown_Correction, ipart, file_O2Analysis, pTbins, nbinpT);
//     H1D_Feeddown_Correction->Scale(-1.);
//     hRawYield_vsPt->Add(H1D_Feeddown_Correction);
//     // cout << "test00000" << endl;
//     for(int ibinPt = 1; ibinPt <= nbinpT; ibinPt++){
//     cout << "H1D_Feeddown_Correction(" << ibinPt << ") = " << H1D_Feeddown_Correction->GetBinContent(ibinPt) << endl;
//     }
//   }
//   delete H1D_Feeddown_Correction;


//   cout << "--------------------------------------------------------------------- POST UNFOLDING TEST" << hist_unfold->GetEntries() << endl;
//   for(int ibinPt = 1; ibinPt <= nbinpT; ibinPt++){
//     cout << "truth_rebinned(" << ibinPt << ") = " << truth_rebinned->GetBinContent(ibinPt) << endl;
//     cout << "H2D_PtResponseMatrix(" << ibinPt << ") = " << H2D_PtResponseMatrix->GetBinContent(ibinPt,ibinPt) << endl;
//     cout << "hist_unfold(" << ibinPt << ") = " << hist_unfold->GetBinContent(ibinPt) << endl;
//     cout << "measured(" << ibinPt << ") = " << measured->GetBinContent(ibinPt) << endl;
//     cout << "measuredPostUnfoldPostCorrections(" << ibinPt << ") = " << hRawYield_vsPt->GetBinContent(ibinPt) << endl;
//   }

//   // TUnfold unfold(H2D_PtResponseMatrix, TUnfold::kHistMapOutputHoriz, TUnfold::kRegModeNone); // kRegModeNone regularisation method choice, here none, will see what's best later
//   // // set input distribution and bias scale (=0)
//   // cout << "AMIERIC DEBUG unfolding 4" << endl;
//   // if(unfold.SetInput(hRawYield_vsPt, 0.0)>=10000) {
//   //   std::cout<<"Unfolding result may be wrong\n";
//   // }
//   // cout << "AMIERIC DEBUG unfolding 5" << endl;

//   // // do the unfolding here
//   // double tauMin=0.0;
//   // double tauMax=0.0;
//   // int nScan=30;
//   // int iBest;
//   // TSpline *logTauX,*logTauY;
//   // TGraph *lCurve;
//   // // this method scans the parameter tau
//   // // finally, the unfolding is done for the "best" choice of tau
//   // iBest=unfold.ScanLcurve(nScan,tauMin,tauMax,&lCurve,&logTauX,&logTauY);
//   // std::cout<<"tau="<<unfold.GetTau()<<"\n";
//   // std::cout<<"chi**2="<<unfold.GetChi2A()<<"+"<<unfold.GetChi2L()
//   //          <<" / "<<unfold.GetNdf()<<"\n";
//   // // save point corresponding to the kink in the L curve as TGraph
//   // double t[1],x[1],y[1];
//   // logTauX->GetKnot(iBest,t[0],x[0]);
//   // logTauY->GetKnot(iBest,t[0],y[0]);
//   // TGraph *bestLcurve=new TGraph(1,x,y);
//   // TGraph *bestLogTauX=new TGraph(1,t,x);

//   // //============================================================
//   // // extract unfolding results into histograms
//   // // set up a bin map, excluding underflow and overflow bins
//   // // the binMap relates the the output of the unfolding to the final
//   // // histogram bins
//   // int *binMap=new int[nbinpT+2];
//   // for(int i=1;i<=nbinpT;i++) binMap[i]=i;
//   // binMap[0]=-1;
//   // binMap[nbinpT+1]=-1;
//   // TH1D *histPtUnfold=new TH1D("Unfolded",";pT(gen)",nbinpT,pTbins);
//   // unfold.GetOutput(histPtUnfold,binMap);
//   // TH1D *histPtDetFold=new TH1D("FoldedBack","pT(det)",nbinpT,pTbins);
//   // unfold.GetFoldedOutput(histPtDetFold); // Aimeric: probably for comparison? some kind of closure test
//   // // store global correlation coefficients
//   // TH1D *histRhoi=new TH1D("rho_I","mass",nbinpT,pTbins);
//   // unfold.GetRhoI(histRhoi,binMap);
//   // delete[] binMap;
//   // binMap=0;
// }
