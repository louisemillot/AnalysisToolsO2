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

void SetStyle(Bool_t graypalette=kFALSE);
void myPadSetUp(TPad *currentPad, float currentLeft=0.11, float currentTop=0.04, float currentRight=0.04, float currentBottom=0.15);
void DrawLogo (Int_t logo=0, Double_t xmin =  0.28, Double_t ymin= 0.68) ;
void FakeHistosOnlyForExample(TH1*&hstat, TH1*&hsyst, TH1*&hsystCorr);
void LoadLibs();
void ArmenterosPodolanski_PreAnalyserCuts(TH2D* &hstat, TH2D* &hsyst, TH2D*&hsystCorr, TFile* file_O2Analysis, TString* SaveAs_Title, TString* &texXtitle, TString* &texYtitle);
void ArmenterosPodolanski_PostAnalyserCuts(TH2D* &hstat, TH2D* &hsyst, TH2D*&hsystCorr, TFile* file_O2Analysis, TString* SaveAs_Title, TString* &texXtitle, TString* &texYtitle);
void ArmenterosPodolanski_PostAnalyserCuts_MC(TH2D* &hstat, TH2D* &hsyst, TH2D*&hsystCorr, TFile* file_O2Analysis, TString* SaveAs_Title, TString* &texXtitle, TString* &texYtitle);
void ArmenterosPodolanski_apass2Format(TH2D* &hstat, TH2D* &hsyst, TH2D*&hsystCorr, TFile* file_O2Analysis, TString* SaveAs_Title, TString* &texXtitle, TString* &texYtitle);

// Preferred colors and markers
const Int_t fillColors[] = {kGray+1,  kRed-10, kBlue-9, kGreen-8, kMagenta-9, kOrange-9,kCyan-8,kYellow-7}; // for syst bands
const Int_t colors[]     = {kBlack, kRed+1 , kBlue+1, kGreen+3, kMagenta+1, kOrange-1,kCyan+2,kYellow+2};
const Int_t markers[]    = {kFullCircle, kFullSquare,kOpenCircle,kOpenSquare,kOpenDiamond,kOpenCross,kFullCross,kFullDiamond,kFullStar,kOpenStar};

//Options:
// TFile* file_O2Analysis = new TFile("AnalysisResults_LHC21k6.root");
TFile* file_O2Analysis = new TFile("AnalysisResult_TpcXrows_Ref.root");

Int_t RebinFactor = 5;


// const Float_t MassPart[numPart] = {0.497611, 1.115683, 1.115683};// 1.32171, 1.32171, 1.67245, 1.67245};
// const Float_t WidthPartLimitsFit[2][numPart] = {{0.003, 0.001, 0.001}, {0.01, 0.005, 0.005}};//{{0.003, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001}, {0.01, 0.005, 0.005, 0.005, 0.005, 0.008, 0.008}};

// === Commonly used x/ titles: ===
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
TString* texRawYield = new TString("O2 raw yield 1/#it{N}_{ev} d#it{N}/d#it{p}_{T}");
TString* texCount = new TString("count");
TString* texSgimaGaussFit = new TString("#sigma of Gaussian fit (GeV/#it{c}^{2})");
TString* texEfficiency = new TString("Efficiency #it{V0}_{detected}/#it{V0}_{MC} (%)");

TString* texQt = new TString("#it{q}_{T} (GeV/#it{c})");
TString* texAlpha = new TString("#alpha");

void ArmenterosPodolanskiPlot_Original() {
  // Load necessary libraries
  LoadLibs();
  // Set the default style
  SetStyle();

  //histograms definition and processing
  TH2D* hstat, *hsyst, *hsystCorr;
  Int_t icolor=0;

  TString* SaveAs_Title = new TString("") ;
  TString* texXtitle = new TString("") ;
  TString* texYtitle = new TString("") ;
  // PtDifferential_SigmaOfFit(hstat, hsyst, hsystCorr, ipart, file_O2Analysis, SaveAs_Title, texXtitle, texYtitle, pTbins[ipart] ,nbinpT[ipart]);
  // PtDistribution_O2data(hstat, hsyst, hsystCorr, ipart, file_O2Analysis, SaveAs_Title);
  // PseudoEfficiency_HEPcomparison_withFit(hstat, hsyst, hsystCorr, ipart, file_O2Analysis, SaveAs_Title, texXtitle, texYtitle, pTbins[ipart], nbinpT[ipart]);
  // PseudoEfficiency_HEPcomparison_allCountSignalRegion(hstat, hsyst, hsystCorr, ipart, file_O2Analysis, SaveAs_Title, texXtitle, texYtitle, pTbins[ipart], nbinpT[ipart]);
  // PseudoEfficiency_HEPcomparison_allCountSignalRegionOLD_withBugs(hstat, hsyst, hsystCorr, ipart, file_O2Analysis, SaveAs_Title, pTbins[ipart], nbinpT[ipart]);
  // RawSpectrum_O2data_withFit(hstat, hsyst, hsystCorr, ipart, file_O2Analysis, SaveAs_Title, texXtitle, texYtitle, pTbins[ipart], nbinpT[ipart]);
  // Efficiency_O2MCdata_withFit(hstat, hsyst, hsystCorr, ipart, file_O2Analysis, SaveAs_Title, texXtitle, texYtitle, pTbins[ipart], nbinpT[ipart]);
  // Efficiency_O2data_allCountSignalRegion(hstat, hsyst, hsystCorr, ipart, file_O2Analysis, SaveAs_Title, texXtitle, texYtitle, pTbins[ipart], nbinpT[ipart]);
  // ArmenterosPodolanski_apass2Format(hstat, hsyst, hsystCorr, file_O2Analysis, SaveAs_Title, texXtitle, texYtitle);
  ArmenterosPodolanski_PostAnalyserCuts(hstat, hsyst, hsystCorr, file_O2Analysis, SaveAs_Title, texXtitle, texYtitle);
  // ArmenterosPodolanski_PreAnalyserCuts(hstat, hsyst, hsystCorr, file_O2Analysis, SaveAs_Title, texXtitle, texYtitle);
  
  // RawSpectrum_O2data_truePt_trueV0s(hstat, hsyst, hsystCorr, ipart, file_O2Analysis, SaveAs_Title, texXtitle, texYtitle, pTbins[ipart], nbinpT[ipart]);

 // Prepare Figure, please stick to the default canvas size(s) unless absolutely necessary in your case
  // Rectangular
  TCanvas *cfig = new TCanvas("cfig", "Alice Figure Template", 800, 600); 
  // Square
  //TCanvas *cfig = new TCanvas("cfig", "Alice Figure Template", 800, 800); 
  // cfig->SetLogy();

  cfig->cd();
  // Set Titles etc..
  // DrawFrame(xmin, ymin, xmax, ymax)
  // TH1 * h = cfig->DrawFrame(0,0.0001,5,0.02);// RawSpectrum_O2data_allCountSignalRegion LOGY scale
  // TH1 * h = cfig->DrawFrame(0,0,3,0.1); // RawSpectrum_O2data_K0S
  // TH1 * h = cfig->DrawFrame(0,0,3,0.001); // RawSpectrum_O2data_Lambda
  // TH1 * h = cfig->DrawFrame(0,0,3,0.06); // PseudoEfficiency_HEPcomparison_withFit
  // TH1 * h = cfig->DrawFrame(0,0,3,0.017); // PtDifferential_SigmaOfFit
  // TH1 * h = cfig->DrawFrame(0,0,3,0.15); // PseudoEfficiency_HEPcomparison_allCountSignalRegion
  // TH2 * h = cfig->DrawFrame(0,0,3,0.7); // PseudoEfficiency_HEPcomparison_allCountSignalRegion
  // TH1 * h = cfig->DrawFrame(0,0,10,0.05); // Momentum Distribution
  // TH1 * h = cfig->DrawFrame(0.45,0,0.55,60); // InvMass_O2data K0S
  // TH1 * h = cfig->DrawFrame(LowerLimitDisplay_Mass[ipart],0,UpperLimitDisplay_Mass[ipart],1.1*hstat->GetMaximum()); // InvMass_O2data Lambda
  // TH1 * h = cfig->DrawFrame(0,0,0.1,500); // DcaHistoProcessing_O2data

  // if (0)
  //   {
  //     TLegend *legt = new TLegend( 0.46, 0.6, 0.72, 0.8, "Test Labels Legend");
  //     legt->AddEntry((TObject*)0, texMtX, "");
  //     legt->AddEntry((TObject*)0, texMtY, "");
  //     legt->AddEntry((TObject*)0, texMassX, "");
  //     legt->AddEntry((TObject*)0, texMeanPt, "");
  //     legt->AddEntry((TObject*)0, texMeanNpart, "");
  //     legt->SetFillColor(0);
  //     legt->SetTextSize(gStyle->GetTextSize()*0.6);
  //     legt->Draw();
  //   }

  // // Set titles
  // // h->SetXTitle(texPtX);
  // cout << "title: " << texXtitle->Data() <<endl;
  // h->SetXTitle(texXtitle->Data());
  // // Please be consistent on the y label
  // // h->SetYTitle(texptDifferentialYield);//(texPtY);
  // h->SetYTitle(texYtitle->Data());//(texPtY);

//   h->GetXaxis()->SetTitle(texXtitle);
//   h->GetYaxis()->SetTitle(texYtitle);

  // Draw your histos here:
  // hsystCorr->SetFillColor(fillColors[icolor]);
  // hsyst    ->SetFillColor(fillColors[icolor]);
  // // hsystCorr->Draw("E3,same"); //SHOULD BE PLOTTED EVENTUALLY
  // hsyst->SetFillStyle(0); // To draw empty boxes
  // hsyst->SetLineColor(colors[icolor]); // To draw empty boxes
  // hsystCorr->SetLineColor(colors[icolor]); // To draw empty boxes
  // hsyst->Draw("E2,same"); //SHOULD BE PLOTTED EVENTUALLY
  hstat->Draw("COLZ");
  hstat->SetMarkerStyle(markers[0]);
  hstat->SetXTitle(texXtitle->Data());
  hstat->SetYTitle(texYtitle->Data());
  hstat->GetXaxis()->SetLabelSize(0.047);
  hstat->GetYaxis()->SetLabelSize(0.047);
  hstat->SetTitleSize(0.055,"xy");

  // use the same color for markers and lines
  hstat->SetMarkerColor(colors [icolor]);
  hstat->SetLineColor  (colors [icolor]);

  // Draw the logo   
  //  0: Just "ALICE" (for final data), to be added only if ALICE does not appear otherwise (e.g. in the legend)
  //  >0: ALICE Preliminary
  // DrawLogo(1, 0.59, 0.81);

  // You should always specify the colliding system
  // NOTATION: pp, p-Pb, Pb-Pb. 
  // Don't forget to use #sqrt{s_{NN}} for p-Pb and Pb-Pb
  // You can change the position of this with
  TPave *Pave = new TPave(0.17,0.87,0.61,0.73,1,"NDC,NB");
  Pave->Draw();
  Pave->SetFillColor(0);
  TLatex * textContext = new TLatex (0.18,0.82,"#bf{ALICE Performance}"); //BOLD
  textContext->SetTextSize(0.05);
  textContext->SetNDC(kTRUE); //remove if I want x,y in TLatex to be in the coordinate system of the histogram
  textContext->Draw();
  TLatex * textColl = new TLatex (0.18,0.75,"pp #sqrt{#it{s}} = 900 GeV, pilot beam 2021");
  textColl->SetTextSize(0.04);
  textColl->SetNDC(kTRUE); //remove if I want x,y in TLatex to be in the coordinate system of the histogram
  textColl->Draw();
  // TLatex * text2 = new TLatex (0.55,55,"V0A Multiplicity Classes (Pb-Side)");
  // text2->SetTextSizePixels(24);
  // text2->Draw();
 
  // //Legend, if needed
  // TLegend * leg = new TLegend(  0.19,  0.19,  0.57, 0.42);
  // leg->AddEntry(hstat,     "0-5\%, stat errors",   "LPE");
  // leg->AddEntry(hsyst,     "syst error (Uncorrelated)",  "F");
  // leg->AddEntry(hsystCorr, "syst error (Correlated)",    "F" );
  // leg->SetFillColor(0);
  // leg->SetTextSize(gStyle->GetTextSize()*0.8);
  // leg->Draw();

  // // Save to HEP data

  // AliHEPDataParser * hepParser = new AliHEPDataParser(hstat, hsyst);
  // hepParser->SetTitle("pt distribution of pi+-, arXiv:XXXX.YYYY");
  // hepParser->SetName("1/Nev 1/p_T 1/2pi d^2N/(dp_Tdy) (GeV/c)^{-1}"); 
  // hepParser->SetXaxisName("PT IN GEV/c");
  // hepParser->SetReaction("RE: P PB --> PI + X");
  // hepParser->SetEnergy("SQRT(SNN) : 5020.0 GeV");
  // hepParser->SetRapidityRange("YRAP : -0.5 - +0.5");
  // hepParser->SaveHEPDataFile("figTemplateHEPData.txt");    // it must be specified explicity if graphs are to be used



  // cfig->SaveAs("HEPcomp_ratio_RawSpectrum_PtDifferential_K0sCount_WithSel8.pdf","pdf"); //HEPcomp_Ratio_
  // cfig->SaveAs("RawSpectrumLog_PtDifferential_K0sCount_WithSel8.pdf","pdf"); //HEPcomp_Ratio_
  // cfig->SaveAs("DCA_K0S_WithSel8.pdf","pdf"); 
  // cfig->SaveAs(strcat(SaveAs_Title, ".pdf"),"pdf"); 
  cfig->SaveAs(*SaveAs_Title+".pdf","pdf"); 

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

void FakeHistosOnlyForExample(TH1* &hstat, TH1* &hsyst, TH1*&hsystCorr) {

  TF1 * fExpo = new TF1 ("fExpo", "expo");
  fExpo->SetParameters(10, -0.3);
  hstat     = new TH1F("hstat", "hstat", 100, 0, 10);
  hsyst     = new TH1F("hsyst", "hsyst", 100, 0, 10);
  hsystCorr = new TH1F("hsystCorr", "hsystCorr", 100, 0, 10);
  hstat     ->Sumw2(); 
  hsyst     ->Sumw2();
  hsystCorr ->Sumw2(); 
  hstat->FillRandom("fExpo",20000);
  Int_t nbinx = hstat->GetNbinsX();
  
  for(Int_t ibinx = 1; ibinx <= nbinx; ibinx++){
    hsyst->SetBinContent (ibinx, hstat->GetBinContent(ibinx));
    hsyst->SetBinError   (ibinx, hstat->GetBinContent(ibinx)*0.08);
    hsystCorr->SetBinContent (ibinx, hstat->GetBinContent(ibinx));
    hsystCorr->SetBinError   (ibinx, hstat->GetBinContent(ibinx)*0.15);
  } 
}

void ArmenterosPodolanski_apass2Format(TH2D* &hstat, TH2D* &hsyst, TH2D*&hsystCorr, TFile* file_O2Analysis, TString* SaveAs_Title, TString* &texXtitle, TString* &texYtitle) {

  //O2 file
  // TFile* file_O2Analysis = new TFile("AnalysisResults_WithSel8Cut.root");

  TH2D* H2_Armenteros = (TH2D*)file_O2Analysis->Get("lambdakzero-qa/hArmenteros");

  TH2D* H2_Armenteros_RebinX = (TH2D*)H2_Armenteros->RebinX(5,"H2_Armenteros_RebinX");
  TH2D* H2_Armenteros_RebinXY = (TH2D*)H2_Armenteros_RebinX->RebinY(5,"H2_Armenteros_RebinXY");

  hstat = (TH2D*)H2_Armenteros_RebinXY->Clone("hstat");


  //////Error Bars///////
  hsyst     = new TH2D("hsyst", "hsyst", 100, 0, 10,100,0,10);
  hsystCorr = new TH2D("hsystCorr", "hsystCorr", 100, 0, 10,100,0,10);
  hstat     ->Sumw2(); 
  hsyst     ->Sumw2();
  hsystCorr ->Sumw2(); 
  Int_t nbinx = hstat->GetNbinsX();
  
  // for(Int_t ibinx = 1; ibinx <= nbinx; ibinx++){
  //   hsyst->SetBinContent (ibinx, hstat->GetBinContent(ibinx));
  //   hsyst->SetBinError   (ibinx, hstat->GetBinContent(ibinx)*0.08);
  //   hsystCorr->SetBinContent (ibinx, hstat->GetBinContent(ibinx));
  //   hsystCorr->SetBinError   (ibinx, hstat->GetBinContent(ibinx)*0.15);
  // } 

  *SaveAs_Title += "ArmenterosPodolanski_PreAnalyserCuts";
  texXtitle = texAlpha;
  texYtitle = texQt;
}

void ArmenterosPodolanski_PreAnalyserCuts(TH2D* &hstat, TH2D* &hsyst, TH2D*&hsystCorr, TFile* file_O2Analysis, TString* SaveAs_Title, TString* &texXtitle, TString* &texYtitle) {

  //O2 file
  // TFile* file_O2Analysis = new TFile("AnalysisResults_WithSel8Cut.root");

  TH2D* H2_Armenteros_preCuts = (TH2D*)file_O2Analysis->Get("lambdakzero-qa/hArmenterosPreAnalyserCuts");

  TH2D* H2_Armenteros_preCuts_RebinX = (TH2D*)H2_Armenteros_preCuts->RebinX(RebinFactor,"H2_Armenteros_preCuts_RebinX");
  TH2D* H2_Armenteros_preCuts_RebinXY = (TH2D*)H2_Armenteros_preCuts_RebinX->RebinY(RebinFactor,"H2_Armenteros_preCuts_RebinXY");

  hstat = (TH2D*)H2_Armenteros_preCuts_RebinXY->Clone("hstat");


  //////Error Bars///////
  hsyst     = new TH2D("hsyst", "hsyst", 100, 0, 10,100,0,10);
  hsystCorr = new TH2D("hsystCorr", "hsystCorr", 100, 0, 10,100,0,10);
  hstat     ->Sumw2(); 
  hsyst     ->Sumw2();
  hsystCorr ->Sumw2(); 
  Int_t nbinx = hstat->GetNbinsX();
  
  // for(Int_t ibinx = 1; ibinx <= nbinx; ibinx++){
  //   hsyst->SetBinContent (ibinx, hstat->GetBinContent(ibinx));
  //   hsyst->SetBinError   (ibinx, hstat->GetBinContent(ibinx)*0.08);
  //   hsystCorr->SetBinContent (ibinx, hstat->GetBinContent(ibinx));
  //   hsystCorr->SetBinError   (ibinx, hstat->GetBinContent(ibinx)*0.15);
  // } 

  *SaveAs_Title += "ArmenterosPodolanski_PreAnalyserCuts";
  texXtitle = texAlpha;
  texYtitle = texQt;//texQt;
}

void ArmenterosPodolanski_PostAnalyserCuts(TH2D* &hstat, TH2D* &hsyst, TH2D*&hsystCorr, TFile* file_O2Analysis, TString* SaveAs_Title, TString* &texXtitle, TString* &texYtitle) {

  //O2 file
  // TFile* file_O2Analysis = new TFile("AnalysisResults_WithSel8Cut.root");

  TH2D* H2_Armenteros_postCuts = (TH2D*)file_O2Analysis->Get("lambdakzero-analysis-mc/hArmenterosPostAnalyserCuts");

  TH2D* H2_Armenteros_postCuts_RebinX = (TH2D*)H2_Armenteros_postCuts->RebinX(RebinFactor,"H2_Armenteros_postCuts_RebinX");
  TH2D* H2_Armenteros_postCuts_RebinXY = (TH2D*)H2_Armenteros_postCuts_RebinX->RebinY(RebinFactor,"H2_Armenteros_postCuts_RebinXY");

  hstat = (TH2D*)H2_Armenteros_postCuts_RebinXY->Clone("hstat");


  //////Error Bars///////
  hsyst     = new TH2D("hsyst", "hsyst", 100, 0, 10,100,0,10);
  hsystCorr = new TH2D("hsystCorr", "hsystCorr", 100, 0, 10,100,0,10);
  hstat     ->Sumw2(); 
  hsyst     ->Sumw2();
  hsystCorr ->Sumw2(); 
  Int_t nbinx = hstat->GetNbinsX();
  
  // for(Int_t ibinx = 1; ibinx <= nbinx; ibinx++){
  //   hsyst->SetBinContent (ibinx, hstat->GetBinContent(ibinx));
  //   hsyst->SetBinError   (ibinx, hstat->GetBinContent(ibinx)*0.08);
  //   hsystCorr->SetBinContent (ibinx, hstat->GetBinContent(ibinx));
  //   hsystCorr->SetBinError   (ibinx, hstat->GetBinContent(ibinx)*0.15);
  // } 

  *SaveAs_Title += "ArmenterosPodolanski_PostAnalyserCuts";
  texXtitle = texAlpha;
  texYtitle = texQt;
}

void ArmenterosPodolanski_PostAnalyserCuts_MC(TH2D* &hstat, TH2D* &hsyst, TH2D*&hsystCorr, TFile* file_O2Analysis, TString* SaveAs_Title, TString* &texXtitle, TString* &texYtitle) {

  //O2 file
  // TFile* file_O2Analysis = new TFile("AnalysisResults_WithSel8Cut.root");

  TH2D* H2_Armenteros_postCuts_MC = (TH2D*)file_O2Analysis->Get("lambdakzero-analysis-mc/hArmenterosPostAnalyserCuts_MC");

  TH2D* H2_Armenteros_postCuts_MC_RebinX = (TH2D*)H2_Armenteros_postCuts_MC->RebinX(RebinFactor,"H2_Armenteros_postCuts_MC_RebinX");
  TH2D* H2_Armenteros_postCuts_MC_RebinXY = (TH2D*)H2_Armenteros_postCuts_MC_RebinX->RebinY(RebinFactor,"H2_Armenteros_postCuts_MC_RebinXY");

  hstat = (TH2D*)H2_Armenteros_postCuts_MC_RebinXY->Clone("hstat");


  //////Error Bars///////
  hsyst     = new TH2D("hsyst", "hsyst", 100, 0, 10,100,0,10);
  hsystCorr = new TH2D("hsystCorr", "hsystCorr", 100, 0, 10,100,0,10);
  hstat     ->Sumw2(); 
  hsyst     ->Sumw2();
  hsystCorr ->Sumw2(); 
  Int_t nbinx = hstat->GetNbinsX();
  
  // for(Int_t ibinx = 1; ibinx <= nbinx; ibinx++){
  //   hsyst->SetBinContent (ibinx, hstat->GetBinContent(ibinx));
  //   hsyst->SetBinError   (ibinx, hstat->GetBinContent(ibinx)*0.08);
  //   hsystCorr->SetBinContent (ibinx, hstat->GetBinContent(ibinx));
  //   hsystCorr->SetBinError   (ibinx, hstat->GetBinContent(ibinx)*0.15);
  // } 

  *SaveAs_Title += "ArmenterosPodolanski_PostAnalyserCuts_MC";
  texXtitle = texAlpha;
  texYtitle = texQt;
}