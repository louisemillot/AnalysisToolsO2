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
#include "./nSubJettiness_settings.h"
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

void SetStyle(Bool_t graypalette=kFALSE);
void LoadLibs();

void Draw_TH1_Histograms_in_one(TH1D** histograms_collection, const TString* legendList_string, int collectionSize, const TString Context, TString* pdfName, TString* &texXtitle, TString* &texYtitle);
void Draw_DeltaR_nSubAxes(int iDataset, double* subRatioRange, TString* &texXtitle, TString* &texYtitle);
void Draw_NSubjettiness(int iDataset, double* deltaRRange, TString* &texXtitle, TString* &texYtitle);


/////////////////////////////////////////////////////
///////////////////// Main Macro ////////////////////
/////////////////////////////////////////////////////

void nSubJettiness() {
  // Load necessary libraries
  LoadLibs();
  // Set the default style
  SetStyle();

  TString* SaveAs_Title = new TString("");
  TString* texXtitle = new TString("");
  TString* texYtitle = new TString("");
  TString* Extra = new TString("");

  double subRatioRange [2] = {0, 0.98};
  // double subRatioRange [2] = {0, 1.5};
  double deltaRRange [2] = {0, 0.5};

  for(int iJetType = 0; iJetType < nJetType; iJetType++) {
    for(int j = 0; j < nBinSubRatio[iJetType] + 1; j++) {
      SubRatioBinEdges[iJetType][j] = j * nSubRatioMax/nBinSubRatio[iJetType];  
    }
  }

  for(int iDataset = 0; iDataset < nDatasets; iDataset++) {
    Draw_NSubjettiness(iDataset, subRatioRange, texXtitle, texYtitle);
    Draw_DeltaR_nSubAxes(iDataset, subRatioRange, texXtitle, texYtitle);
  }
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
  gStyle->SetLabelSize(0.045,"xyz");
  gStyle->SetLabelOffset(0.01,"y");
  gStyle->SetLabelOffset(0.01,"x");
  gStyle->SetLabelColor(kBlack,"xyz");
  gStyle->SetTitleSize(0.05,"xyz");
  gStyle->SetTitleOffset(1.25,"y");
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

//should eventually have textContext variables be written using functions in here. See JetQC.C for examples

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////// nSubjettiness  plot functions ///////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void Draw_NSubjettiness(int iDataset, double* subRatioRange, TString* &texXtitle, TString* &texYtitle) {

  TH2D* H2D_NSubjettiness[nAxisTypes];
  TH1D* H1D_NSubjettiness_projectedY[nAxisTypes];
  TH1D* H1D_NSubjettiness_projectedY_rebinnedY[nAxisTypes];

  for(int iAxisType = 0; iAxisType < nAxisTypes; iAxisType++){
    // H2D_NSubjettiness[iAxisType] = (TH2D*)file_O2Analysis_list[iDataset]->Get("jet-nsubjettiness-"+jetType[iJetType]+"-"+jetLevel[iJetLevel]+"/hNSubRatio21VsPt_"+AxisType[iAxisType]);
    H2D_NSubjettiness[iAxisType] = (TH2D*)file_O2Analysis_list[iDataset]->Get("jet-nsubjettiness/hNSubRatio21VsPt_"+AxisType[iAxisType]);
    // H2D_NSubjettiness[iAxisType]->Sumw2();

    int ibinPt_low = H2D_NSubjettiness[iAxisType]->GetXaxis()->FindBin(PtCutLow+GLOBAL_epsilon);
    int ibinPt_high = H2D_NSubjettiness[iAxisType]->GetXaxis()->FindBin(PtCutHigh+GLOBAL_epsilon);
    int ibinRatio_low = H2D_NSubjettiness[iAxisType]->GetYaxis()->FindBin(subRatioRange[0]+GLOBAL_epsilon);
    int ibinRatio_high = H2D_NSubjettiness[iAxisType]->GetYaxis()->FindBin(subRatioRange[1]+GLOBAL_epsilon);
    cout << "ibinRatio_low = " << ibinRatio_low << ", ibinRatio_high = " << ibinRatio_high << endl;

    H1D_NSubjettiness_projectedY[iAxisType] = (TH1D*)H2D_NSubjettiness[iAxisType]->ProjectionY("jet-nsubjettiness-"+jetType[iJetType]+"-"+jetLevel[iJetLevel]+DatasetsNames[iDataset]+"/hNSubRatio21VsPt_"+AxisType[iAxisType]+"_projectedY",ibinPt_low,ibinPt_high);

    // H1D_NSubjettiness_projectedY_rebinnedY[iAxisType] = (TH1D*)H1D_NSubjettiness_projectedY[iAxisType]->Rebin(nBinSubRatio[iJetType],"H1D_NSubjettiness_projectedY_rebinnedY",SubRatioBinEdges[iJetType]);
    H1D_NSubjettiness_projectedY_rebinnedY[iAxisType] = (TH1D*)H1D_NSubjettiness_projectedY[iAxisType]->Rebin(1,"H1D_NSubjettiness_projectedY_rebinnedY"+AxisType[iAxisType]+DatasetsNames[iDataset]);

    H1D_NSubjettiness_projectedY_rebinnedY[iAxisType]->GetXaxis()->SetRange(ibinRatio_low,ibinRatio_high);

    cout << "1 getEntries = " << H1D_NSubjettiness_projectedY_rebinnedY[iAxisType]->GetEntries() << endl;
    cout << "1   integral = " << H1D_NSubjettiness_projectedY_rebinnedY[iAxisType]->Integral(ibinRatio_low, ibinRatio_high) << endl;

    double dsubratio, dN_dsubratio;
    for(int iBinSubRatio = 0; iBinSubRatio < nBinSubRatio[iJetType]; iBinSubRatio++){
      dsubratio = H1D_NSubjettiness_projectedY_rebinnedY[iAxisType]->GetXaxis()->GetBinWidth(iBinSubRatio);
      dN_dsubratio = H1D_NSubjettiness_projectedY_rebinnedY[iAxisType]->GetBinContent(iBinSubRatio) *1./dsubratio;
      H1D_NSubjettiness_projectedY_rebinnedY[iAxisType]->SetBinContent(iBinSubRatio, dN_dsubratio);
      H1D_NSubjettiness_projectedY_rebinnedY[iAxisType]->SetBinError(iBinSubRatio, H1D_NSubjettiness_projectedY_rebinnedY[iAxisType]->GetBinError(iBinSubRatio) *1./dsubratio);
    }

    // H1D_NSubjettiness_projectedY_rebinnedY[iAxisType]->Scale(1./H1D_NSubjettiness_projectedY_rebinnedY[iAxisType]->GetEntries());
    // H1D_NSubjettiness_projectedY_rebinnedY[iAxisType]->Scale(1./H1D_NSubjettiness_projectedY_rebinnedY[iAxisType]->Integral(ibinRatio_low, ibinRatio_high));
  }

  texXtitle = texNSubJettinessRatio;
  texYtitle = texdN_dsubratio;
  TString* pdfName = new TString("nSubjettiness_"+jetType[iJetType]+"_"+jetLevel[iJetLevel]+"_AxisTypeComparison_"+Form("%.0f", PtCutLow)+"<Pt<"+Form("%.0f", PtCutHigh)+"GeV_"+DatasetsNames[iDataset]);
  const TString textContext("#splitline{"+jetType[iJetType]+" "+jetLevel[iJetLevel]+" #minus "+DatasetsNames[iDataset]+"}{"+Form("%.0f", PtCutLow)+" < Pt < "+Form("%.0f", PtCutHigh)+" GeV}");

  Draw_TH1_Histograms_in_one(H1D_NSubjettiness_projectedY_rebinnedY, AxisTypeLegend, nAxisTypes, textContext, pdfName, texXtitle, texYtitle, texCollisionDataInfo, "");
  cout << "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX" << endl;
  cout << "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX" << endl;
  cout << "NORMALISATION IS WRONG ; also naming ; im saying 1/Nevts but I dont get the number of events anywhere" << endl;
  cout << "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX" << endl;
  cout << "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX" << endl;
}


void Draw_DeltaR_nSubAxes(int iDataset, double* deltaRRange, TString* &texXtitle, TString* &texYtitle) {

  TH2D* H2D_DeltaR[nAxisTypes];
  TH1D* H1D_DeltaR_projectedY[nAxisTypes];
  TH1D* H1D_DeltaR_projectedY_rebinnedY[nAxisTypes];


  for(int iAxisType = 0; iAxisType < nAxisTypes; iAxisType++){
    // H2D_DeltaR[iAxisType] = (TH2D*)file_O2Analysis_list[iDataset]->Get("jet-nsubjettiness-"+jetType[iJetType]+"-"+jetLevel[iJetLevel]+"/hDeltaRVsPt_"+AxisType[iAxisType]);
    H2D_DeltaR[iAxisType] = (TH2D*)file_O2Analysis_list[iDataset]->Get("jet-nsubjettiness/hDeltaRVsPt_"+AxisType[iAxisType]);
    // H2D_DeltaR[iAxisType]->Sumw2();

    int ibinPt_low = H2D_DeltaR[iAxisType]->GetXaxis()->FindBin(PtCutLow+GLOBAL_epsilon);
    int ibinPt_high = H2D_DeltaR[iAxisType]->GetXaxis()->FindBin(PtCutHigh+GLOBAL_epsilon);
    int ibinDeltaR_low = H2D_DeltaR[iAxisType]->GetYaxis()->FindBin(deltaRRange[0]+GLOBAL_epsilon);
    int ibinDeltaR_high = H2D_DeltaR[iAxisType]->GetYaxis()->FindBin(deltaRRange[1]+GLOBAL_epsilon);
    // cout << "ibinDeltaR_low = " << ibinDeltaR_low << ", ibinDeltaR_high = " << ibinDeltaR_high << endl;

    H1D_DeltaR_projectedY[iAxisType] = (TH1D*)H2D_DeltaR[iAxisType]->ProjectionY("jet-nsubjettiness-"+jetType[iJetType]+"-"+jetLevel[iJetLevel]+DatasetsNames[iDataset]+"/hNSubRatio21VsPt_"+AxisType[iAxisType]+"_projectedY",ibinPt_low,ibinPt_high);

    // H1D_DeltaR_projectedY_rebinnedY[iAxisType] = (TH1D*)H1D_NSubjettiness_projectedY[iAxisType]->Rebin(nBinSubRatio[iJetType],"H1D_DeltaR_projectedY_rebinnedY",SubRatioBinEdges[iJetType]);
    H1D_DeltaR_projectedY_rebinnedY[iAxisType] = (TH1D*)H1D_DeltaR_projectedY[iAxisType]->Rebin(1,"H1D_DeltaR_projectedY_rebinnedY"+AxisType[iAxisType]+DatasetsNames[iDataset]);
    // cout << "ibinPt: " << ibinPt << "; originalAxis [low;up] = [" << ibinPt_originalAxis_low << ";" << ibinPt_originalAxis_up << "]" << endl;

    H1D_DeltaR_projectedY_rebinnedY[iAxisType]->GetXaxis()->SetRange(ibinDeltaR_low,ibinDeltaR_high);

    double dsubratio, dN_dsubratio;
    for(int iBinSubRatio = 0; iBinSubRatio < nBinSubRatio[iJetType]; iBinSubRatio++){
      dsubratio = H1D_DeltaR_projectedY_rebinnedY[iAxisType]->GetXaxis()->GetBinWidth(iBinSubRatio); // bin width
      dN_dsubratio = H1D_DeltaR_projectedY_rebinnedY[iAxisType]->GetBinContent(iBinSubRatio) *1./dsubratio;
      H1D_DeltaR_projectedY_rebinnedY[iAxisType]->SetBinContent(iBinSubRatio, dN_dsubratio);
      H1D_DeltaR_projectedY_rebinnedY[iAxisType]->SetBinError(iBinSubRatio, H1D_DeltaR_projectedY_rebinnedY[iAxisType]->GetBinError(iBinSubRatio) *1./dsubratio);
    }
    // H1D_DeltaR_projectedY_rebinnedY[iAxisType]->Scale(1./H1D_DeltaR_projectedY_rebinnedY[iAxisType]->GetEntries());
    // H1D_DeltaR_projectedY_rebinnedY[iAxisType]->Scale(1./H1D_DeltaR_projectedY_rebinnedY[iAxisType]->Integral(ibinDeltaR_low, ibinDeltaR_high));
  }

  texXtitle = texDeltaR;
  texYtitle = texdN_dDeltaR;
  TString* pdfName = new TString("DeltaR_"+jetType[iJetType]+"_"+jetLevel[iJetLevel]+"_AxisTypeComparison_"+Form("%.0f", PtCutLow)+"<Pt<"+Form("%.0f", PtCutHigh)+"GeV_"+DatasetsNames[iDataset]);
  const TString textContext("#splitline{"+jetType[iJetType]+" "+jetLevel[iJetLevel]+" #minus "+DatasetsNames[iDataset]+"}{"+Form("%.0f", PtCutLow)+" < Pt < "+Form("%.0f", PtCutHigh)+" GeV}");

  Draw_TH1_Histograms_in_one(H1D_DeltaR_projectedY_rebinnedY, AxisTypeLegend, nAxisTypes, textContext, pdfName, texXtitle, texYtitle, texCollisionDataInfo, "");
  cout << "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX" << endl;
  cout << "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX" << endl;
  cout << "NORMALISATION IS WRONG ; also naming ; im saying 1/Nevts but I dont get the number of events anywhere" << endl;
  cout << "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX" << endl;
  cout << "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX" << endl;

}

// obsolete, now in HistogramUtilities.C
// void Draw_TH1_Histograms_in_one(TH1D** histograms_collection, const TString* legendList_string, int collectionSize, const TString Context, TString* pdfName, TString* &texXtitle, TString* &texYtitle) {

//   ////plots
//   TCanvas *canvas = new TCanvas ("canvas"+*pdfName, "canvas"+*pdfName, 800, 800);
//   canvas->cd(0);

//   float maxY = 0;
//   float minX = 9999999;
//   float maxX = 0;
//   for (int i = 0; i < collectionSize; i++) {
//     if (maxY < histograms_collection[i]->GetMaximum()) maxY = histograms_collection[i]->GetMaximum();
//     if (minX > histograms_collection[i]->GetXaxis()->GetXmin()) minX = histograms_collection[i]->GetXaxis()->GetXmin();
//     if (maxX < histograms_collection[i]->GetXaxis()->GetXmax()) maxX = histograms_collection[i]->GetXaxis()->GetXmax();
//   }
//   TH1 *hFrame = canvas->DrawFrame(minX,0,maxX,1.8*maxY);
//   hFrame->SetXTitle(texXtitle->Data());
//   hFrame->SetYTitle(texYtitle->Data());
//   TLegend * leg = new TLegend(0.7, 0.75, 0.87, 0.87);

//   //draw histograms from collection, ignoring first one that is the systematics for second one
//   for (int i = 0; i < collectionSize; i++) {
//     // histograms_collection[i]->Draw("hist same p");
//     histograms_collection[i]->Draw("same");
//     histograms_collection[i]->SetMarkerStyle(markers[i]);
//     histograms_collection[i]->SetMarkerColor(colors[i]);
//     histograms_collection[i]->SetLineColor(colors[i]);

//     leg->AddEntry(histograms_collection[i], legendList_string[i], "LP");
//   }

//   leg->SetTextSize(gStyle->GetTextSize()*0.3);
//   cout << "AIMERIC - legend: I should try and change this 0.3 factor to have a nice legend" << endl;
//   if (collectionSize >= 2) {
//     leg->Draw("same");
//   }

//   TLatex * textColl = new TLatex (0.18,0.82,texCollisionDataInfo->Data());
//   textColl->SetTextSize(0.04);
//   textColl->SetNDC(kTRUE); //remove if I want x,y in TLatex to be in the coordinate system of the histogram
//   textColl->Draw();
//   TLatex * text_part = new TLatex (0.18,0.75,Context);
//   text_part->SetTextSize(0.04);
//   text_part->SetNDC(kTRUE); //remove if I want x,y in TLatex to be in the coordinate system of the histogram
//   text_part->Draw();

//   canvas->SaveAs(*pdfName+".pdf");
// }

