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

#include<array>
using namespace std;

void SetStyle(Bool_t graypalette=kFALSE);
void myPadSetUp(TPad *currentPad, float currentLeft=0.11, float currentTop=0.04, float currentRight=0.04, float currentBottom=0.15);
void DrawLogo (Int_t logo=0, Double_t xmin =  0.28, Double_t ymin= 0.68) ;
void FakeHistosOnlyForExample(TH1*&hstat, TH1*&hsyst, TH1*&hsystCorr);
void LoadLibs();

void Draw_TH1_Histograms_in_one(TH1D** histograms_collection, const TString* legendList_string, Int_t collectionSize, const TString Context, TString* pdfName, TString* &texXtitle, TString* &texYtitle);
void Draw_DeltaR_nSubAxes(TString* &texXtitle, TString* &texYtitle);
void Draw_NSubjettiness(TString* &texXtitle, TString* &texYtitle);

// Preferred colors and markers
const Int_t fillColors[] = {kGray+1,  kRed-10, kBlue-9, kGreen-8, kMagenta-9, kOrange-9,kCyan-8,kYellow-7}; // for syst bands
const Int_t colors[]     = {kRed+1, kBlack, kBlue+1, kGreen+3, kMagenta+1, kOrange-1,kCyan+2,kYellow+2};
const Int_t markers[]    = {kOpenCircle,kFullCircle,kFullSquare,kOpenSquare,kOpenDiamond,kOpenCross,kFullCross,kFullDiamond,kFullStar,kOpenStar};

// Options to be set:
// TFile* file_O2Analysis = new TFile("AnalysisResults_LHC23d4_fullProdHyperloop.root");
TFile* file_O2Analysis = new TFile("AnalysisResults.root");
TString* texCollisionDataInfo = new TString("#splitline{pp #sqrt{#it{s}} = 13.6 TeV}{Jet-Jet simulation LHC23d4}");
const Int_t iJetType = 1;
const Int_t iJetLevel = 2;
const Float_t PtCutLow  = 20;
const Float_t PtCutHigh = 40;

// Permanent Options
const Int_t nJetType = 4;
const TString jetType[nJetType] = {"charged", "d0", "neutral", "full"};
const Int_t nJetLevel = 3;
const TString jetLevel[nJetLevel] = {"data", "mcd", "mcp"};
const Int_t nAxisTypes = 3;
const TString AxisType[nAxisTypes] = {"Kt", "CA", "SD"};
const TString AxisTypeLegend[nAxisTypes] = {"#it{k}_{T} axes", "C/A axes", "Soft Drop axes"};

Int_t nBinSubRatio[nJetLevel] = {100, 100, 100};
Double_t SubRatioBinEdges[nJetLevel][200]; // = {{0},{0},{0}}; // 200 just gives some margin
const Float_t nSubRatioMax = 1.2;

// const Float_t drapidity = 1.; // for now I'm not looking at rapidity differential

// Double_t pTbins[nJetType][20] = {{0.0, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.6, 2.0, 2.4, 3.0},{0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 2.0, 2.4, 3.0},{0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 2.0, 2.4, 3.0}};
// Int_t nbinpT[nJetType] = {17,8,8};




// === Commonly used x/ titles: ===
// rapidity
TString* texRapidity = new TString("y");
// pt invariant yields
TString* texPtX = new TString("#it{p}_{T} (GeV/#it{c})");
TString* texPtY = new TString("1/#it{N}_{ev} 1/(2#pi#it{p}_{T}) d#it{N}/(d#it{p}_{T}d#it{y}) ((GeV/#it{c})^{-2})");
// mt invariant yields
TString* texMtX = new TString("#it{m}_{T} (GeV/#it{c}^{2})");
TString* texMtY = new TString("1/#it{N}_{ev} 1/(2#pi#it{m}_{T}) d#it{N}/(d#it{m}_{T}d#it{y}) ((GeV/#it{c}^{2})^{-2})"); 
// Invariant mass with decay products K and pi
TString* texMassX = new TString("#it{M}_{K#pi} (GeV/#it{c}^{2})");
TString* texMassY = new TString("d#it{N}/(d#it{M}_{K#pi})");
// Invariant mass with decay products pi+ and pi-
TString* texMassPiPiX = new TString("#it{M}_{#it{#pi}^{+}#it{#pi}^{-}} (GeV/#it{c}^{2})");
TString* texMassPiPiY = new TString("d#it{N}/(d#it{M}_{#it{#pi}^{+}#it{#pi}^{-}})");
// <pt>, npart
TString* texMeanPt = new TString("#LT#it{p}_{T}#GT (GeV/#it{c})");
TString* texMeanNpart = new TString("#LT#it{N}_{part}#GT");
//AIMERIC TILTES
TString* texptDifferentialYield = new TString("1/#it{N}_{ev} d#it{N}/d#it{p}_{T} (GeV/#it{c})^{-1}");
TString* texptDifferentialYield_HEPratio_pseudoEff = new TString("O2/HEP ratio of 1/#it{N}_{ev} d#it{N}/d#it{p}_{T} (GeV/#it{c})^{-1}");
// TString* texRawYield = new TString("raw yield 1/#it{N}_{ev} d#it{N}/d#it{p}_{T}");
TString* texRawYield = new TString("d#it{N}/d#it{p}_{T}");
TString* texCount = new TString("count");
TString* texSgimaGaussFit = new TString("#sigma of Gaussian fit (GeV/#it{c}^{2})");
TString* texEfficiency = new TString("Efficiency #it{V0}_{detected}/#it{V0}_{MC} (ratio)");
TString* texDaughterRecoEfficiency = new TString("Track Reco Efficiency #it{N}_{reco}/#it{N}_{genMC} (ratio)");
TString* texDaughterPairsRecoEfficiency = new TString("Track Reco Efficiency V0Pairs (ratio)");

TString* texPiMinus = new TString("#pi^{-}");
TString* texPiPlus = new TString("#pi^{+}");

TString* texProton = new TString("p");
TString* texAntiProton = new TString("#bar{p}");

TString* texRatio = new TString("ratio");
TString* texFeeddownRemoved = new TString("Feeddown removed - ratio");

TString* texPseudoEfficiency = new TString("Pseudo efficiency #it{V0}_{detected}/#it{V0}_{expected} (ratio)");

TString* texCountRelativeToMax = new TString("Count normalised to peak maximum");

TString* texPtMeasured = new TString("#it{p}_{T}^{measured} (GeV/#it{c})");
TString* texPtMC = new TString("#it{p}_{T}^{MC} (GeV/#it{c})");

TString* texNSubJettinessRatio = new TString("#it{#tau_{2}/#tau_{1}}");
TString* texDeltaR = new TString("#it{#Delta R}");
TString* texdN_dsubratio = new TString("1/N^{jet} d#it{N}/d(#it{#tau_{2}/#tau_{1}})");
TString* texdN_dDeltaR = new TString("1/N^{jet} d#it{N}/d(#it{#Delta R})");

void nSubJettiness() {
  // Load necessary libraries
  LoadLibs();
  // Set the default style
  SetStyle();

  TString* SaveAs_Title = new TString("");
  TString* texXtitle = new TString("");
  TString* texYtitle = new TString("");
  TString* Extra = new TString("");

  for(Int_t iJetType = 0; iJetType < nJetType; iJetType++) {
    for(Int_t j = 0; j < nBinSubRatio[iJetType] + 1; j++) {
      SubRatioBinEdges[iJetType][j] = j * nSubRatioMax/nBinSubRatio[iJetType];  
    }
  }

  Draw_NSubjettiness(texXtitle, texYtitle);
  Draw_DeltaR_nSubAxes(texXtitle, texYtitle);
}

//________________________________
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

void DrawLogo (Int_t logo, Double_t xmin, Double_t ymin) {

  // Logo is not needed anymore, now we only write alice preliminary
  // Logo:
  // 0: Justr writes "ALICE" (for final data)
  // Anything eles: writes "ALICE Preliminary"

  TLatex *   tex = new TLatex(xmin,ymin, logo ? "ALICE Preliminary" : "ALICE");
  tex->SetNDC();
  tex->SetTextFont(42);
  tex->Draw();

  // OLD logo
  //  TPad * currentPad = gPad;
  // Double_t AliLogo_LowX =xmin;
  // Double_t AliLogo_LowY = ymin;
  // Double_t AliLogo_Height = size;
  // //ALICE logo is a  file that is 821x798 pixels->should be wider than a square
  // Double_t AliLogo_Width  = (821./798.) * AliLogo_Height * gPad->GetWh() / gPad->GetWw();
  
  // TPad *myPadLogo = new TPad("myPadLogo", "Pad for ALICE Logo",AliLogo_LowX,AliLogo_LowY,AliLogo_LowX+AliLogo_Width,AliLogo_LowY+AliLogo_Height);
  // myPadSetUp(myPadLogo,0,0,0,0);
  // //myPadLogo->SetFixedAspectRatio(1);
  // myPadLogo->Draw();
  // myPadLogo->cd();
  // if (logo == 0) {
  //   myPadLogo->SetFillColor(2); // color to first figure out where is the pad then comment !
  // } else if (logo == 1){
  //   TASImage *myAliceLogo = new TASImage(performanceLogoPath);
  //   myAliceLogo->Draw();
  // } else if (logo == 2) {
  //   TASImage *myAliceLogo = new TASImage(preliminaryLogoPath);
  //   myAliceLogo->Draw();
  // }
  // // go back to the old pad
  // currentPad->cd();

}

void myPadSetUp(TPad *currentPad, float currentLeft, float currentTop, float currentRight, float currentBottom){
  currentPad->SetLeftMargin(currentLeft);
  currentPad->SetTopMargin(currentTop);
  currentPad->SetRightMargin(currentRight);
  currentPad->SetBottomMargin(currentBottom);
  return;
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











void Draw_NSubjettiness(TString* &texXtitle, TString* &texYtitle) {

  TH2D* H2D_NSubjettiness[nAxisTypes];
  TH1D* H1D_NSubjettiness_projectedY[nAxisTypes];
  TH1D* H1D_NSubjettiness_projectedY_rebinnedY[nAxisTypes];

  for(int iAxisType = 0; iAxisType < nAxisTypes; iAxisType++){
    H2D_NSubjettiness[iAxisType] = (TH2D*)file_O2Analysis->Get("jet-nsubjettiness-"+jetType[iJetType]+"-"+jetLevel[iJetLevel]+"/hNSubRatio21VsPt_"+AxisType[iAxisType]);
    // H2D_NSubjettiness[iAxisType]->Sumw2();

    int ibinPt_low = H2D_NSubjettiness[iAxisType]->GetXaxis()->FindBin(PtCutLow);
    int ibinPt_high = H2D_NSubjettiness[iAxisType]->GetXaxis()->FindBin(PtCutHigh) - 1;
    H1D_NSubjettiness_projectedY[iAxisType] = (TH1D*)H2D_NSubjettiness[iAxisType]->ProjectionY("jet-nsubjettiness-"+jetType[iJetType]+"-"+jetLevel[iJetLevel]+"/hNSubRatio21VsPt_"+AxisType[iAxisType]+"_projectedY",ibinPt_low,ibinPt_high);

    // H1D_NSubjettiness_projectedY_rebinnedY[iAxisType] = (TH1D*)H1D_NSubjettiness_projectedY[iAxisType]->Rebin(nBinSubRatio[iJetType],"H1D_NSubjettiness_projectedY_rebinnedY",SubRatioBinEdges[iJetType]);
    H1D_NSubjettiness_projectedY_rebinnedY[iAxisType] = (TH1D*)H1D_NSubjettiness_projectedY[iAxisType]->Rebin(1,"H1D_NSubjettiness_projectedY_rebinnedY");


    Double_t dsubratio, dN_dsubratio;
    for(int iBinSubRatio = 0; iBinSubRatio < nBinSubRatio[iJetType]; iBinSubRatio++){
      dsubratio = H1D_NSubjettiness_projectedY_rebinnedY[iAxisType]->GetXaxis()->GetBinWidth(iBinSubRatio);
      dN_dsubratio = H1D_NSubjettiness_projectedY_rebinnedY[iAxisType]->GetBinContent(iBinSubRatio) *1./dsubratio;
      H1D_NSubjettiness_projectedY_rebinnedY[iAxisType]->SetBinContent(iBinSubRatio, dN_dsubratio);
      H1D_NSubjettiness_projectedY_rebinnedY[iAxisType]->SetBinError(iBinSubRatio, H1D_NSubjettiness_projectedY_rebinnedY[iAxisType]->GetBinError(iBinSubRatio) *1./dsubratio);
    }

    H1D_NSubjettiness_projectedY_rebinnedY[iAxisType]->Scale(1./H1D_NSubjettiness_projectedY_rebinnedY[iAxisType]->GetEntries());
  }

  texXtitle = texNSubJettinessRatio;
  texYtitle = texdN_dsubratio;
  TString* pdfName = new TString("nSubjettiness_"+jetType[iJetType]+"_"+jetLevel[iJetLevel]+"_AxisTypeComparison_"+Form("%.0f", PtCutLow)+"<Pt<"+Form("%.0f", PtCutHigh)+"GeV");
  const TString textContext(jetType[iJetType]+" "+jetLevel[iJetLevel]+" #minus "+Form("%.0f", PtCutLow)+" < Pt < "+Form("%.0f", PtCutHigh)+" GeV");

  Draw_TH1_Histograms_in_one(H1D_NSubjettiness_projectedY_rebinnedY, AxisTypeLegend, nAxisTypes, textContext, pdfName, texXtitle, texYtitle);

}


void Draw_DeltaR_nSubAxes(TString* &texXtitle, TString* &texYtitle) {

  TH2D* H2D_DeltaR[nAxisTypes];
  TH1D* H1D_DeltaR_projectedY[nAxisTypes];
  TH1D* H1D_DeltaR_projectedY_rebinnedY[nAxisTypes];


  for(int iAxisType = 0; iAxisType < nAxisTypes; iAxisType++){
    H2D_DeltaR[iAxisType] = (TH2D*)file_O2Analysis->Get("jet-nsubjettiness-"+jetType[iJetType]+"-"+jetLevel[iJetLevel]+"/hDeltaRVsPt_"+AxisType[iAxisType]);
    // H2D_DeltaR[iAxisType]->Sumw2();

    int ibinPt_low = H2D_DeltaR[iAxisType]->GetXaxis()->FindBin(PtCutLow);
    int ibinPt_high = H2D_DeltaR[iAxisType]->GetXaxis()->FindBin(PtCutHigh) - 1;
    H1D_DeltaR_projectedY[iAxisType] = (TH1D*)H2D_DeltaR[iAxisType]->ProjectionY("jet-nsubjettiness-"+jetType[iJetType]+"-"+jetLevel[iJetLevel]+"/hNSubRatio21VsPt_"+AxisType[iAxisType]+"_projectedY",ibinPt_low,ibinPt_high);

    // H1D_DeltaR_projectedY_rebinnedY[iAxisType] = (TH1D*)H1D_NSubjettiness_projectedY[iAxisType]->Rebin(nBinSubRatio[iJetType],"H1D_DeltaR_projectedY_rebinnedY",SubRatioBinEdges[iJetType]);
    H1D_DeltaR_projectedY_rebinnedY[iAxisType] = (TH1D*)H1D_DeltaR_projectedY[iAxisType]->Rebin(1,"H1D_DeltaR_projectedY_rebinnedY");
    // cout << "ibinPt: " << ibinPt << "; originalAxis [low;up] = [" << ibinPt_originalAxis_low << ";" << ibinPt_originalAxis_up << "]" << endl;

    Double_t dsubratio, dN_dsubratio;
    for(int iBinSubRatio = 0; iBinSubRatio < nBinSubRatio[iJetType]; iBinSubRatio++){
      dsubratio = H1D_DeltaR_projectedY_rebinnedY[iAxisType]->GetXaxis()->GetBinWidth(iBinSubRatio); // bin width
      dN_dsubratio = H1D_DeltaR_projectedY_rebinnedY[iAxisType]->GetBinContent(iBinSubRatio) *1./dsubratio;
      H1D_DeltaR_projectedY_rebinnedY[iAxisType]->SetBinContent(iBinSubRatio, dN_dsubratio);
      H1D_DeltaR_projectedY_rebinnedY[iAxisType]->SetBinError(iBinSubRatio, H1D_DeltaR_projectedY_rebinnedY[iAxisType]->GetBinError(iBinSubRatio) *1./dsubratio);
    }
    H1D_DeltaR_projectedY_rebinnedY[iAxisType]->Scale(1./H1D_DeltaR_projectedY_rebinnedY[iAxisType]->GetEntries());
  }

  texXtitle = texDeltaR;
  texYtitle = texdN_dDeltaR;
  TString* pdfName = new TString("DeltaR_"+jetType[iJetType]+"_"+jetLevel[iJetLevel]+"_AxisTypeComparison_"+Form("%.0f", PtCutLow)+"<Pt<"+Form("%.0f", PtCutHigh)+"GeV");
  const TString textContext(jetType[iJetType]+" "+jetLevel[iJetLevel]+" #minus "+Form("%.0f", PtCutLow)+" < Pt < "+Form("%.0f", PtCutHigh)+" GeV");

  Draw_TH1_Histograms_in_one(H1D_DeltaR_projectedY_rebinnedY, AxisTypeLegend, nAxisTypes, textContext, pdfName, texXtitle, texYtitle);

}

void Draw_TH1_Histograms_in_one(TH1D** histograms_collection, const TString* legendList_string, Int_t collectionSize, const TString Context, TString* pdfName, TString* &texXtitle, TString* &texYtitle) {

  ////plots
  TCanvas *canvas = new TCanvas ("canvas"+*pdfName, "canvas"+*pdfName, 800, 800);
  canvas->cd(0);

  float maxY = 0;
  float minX = 9999999;
  float maxX = 0;
  for (Int_t i = 0; i < collectionSize; i++) {
    if (maxY < histograms_collection[i]->GetMaximum()) maxY = histograms_collection[i]->GetMaximum();
    if (minX > histograms_collection[i]->GetXaxis()->GetXmin()) minX = histograms_collection[i]->GetXaxis()->GetXmin();
    if (maxX < histograms_collection[i]->GetXaxis()->GetXmax()) maxX = histograms_collection[i]->GetXaxis()->GetXmax();
  }
  TH1 *hFrame = canvas->DrawFrame(minX,0,maxX,1.8*maxY);
  hFrame->SetXTitle(texXtitle->Data());
  hFrame->SetYTitle(texYtitle->Data());
  TLegend * leg = new TLegend(0.7, 0.75, 0.87, 0.87);

  //draw histograms from collection, ignoring first one that is the systematics for second one
  for (Int_t i = 0; i < collectionSize; i++) {
    // histograms_collection[i]->Draw("hist same p");
    histograms_collection[i]->Draw("same");
    histograms_collection[i]->SetMarkerStyle(markers[i]);
    histograms_collection[i]->SetMarkerColor(colors[i]);
    histograms_collection[i]->SetLineColor(colors[i]);

    leg->AddEntry(histograms_collection[i], legendList_string[i], "LP");
  }

  leg->SetTextSize(gStyle->GetTextSize()*0.3);
  cout << "AIMERIC - legend: I should try and change this 0.3 factor to have a nice legend" << endl;
  if (collectionSize >= 2) {
    leg->Draw("same");
  }

  TLatex * textColl = new TLatex (0.18,0.82,texCollisionDataInfo->Data());
  textColl->SetTextSize(0.04);
  textColl->SetNDC(kTRUE); //remove if I want x,y in TLatex to be in the coordinate system of the histogram
  textColl->Draw();
  TLatex * text_part = new TLatex (0.18,0.75,Context);
  text_part->SetTextSize(0.04);
  text_part->SetNDC(kTRUE); //remove if I want x,y in TLatex to be in the coordinate system of the histogram
  text_part->Draw();

  canvas->SaveAs(*pdfName+".pdf");
}

