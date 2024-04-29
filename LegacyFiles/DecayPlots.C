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
void DecayLength(TH1D* &hstat, TH1D* &hsyst, TH1D*&hsystCorr, TFile* file_O2Analysis, TString* SaveAs_Title, TString* &texXtitle, TString* &texYtitle);
void CTau(TH1D* &hstat, TH1D* &hsyst, TH1D*&hsystCorr, TFile* file_O2Analysis, TString* SaveAs_Title, TString* &texXtitle, TString* &texYtitle);

// Preferred colors and markers
const Int_t fillColors[] = {kGray+1,  kRed-10, kBlue-9, kGreen-8, kMagenta-9, kOrange-9,kCyan-8,kYellow-7}; // for syst bands
const Int_t colors[]     = {kBlack, kRed+1 , kBlue+1, kGreen+3, kMagenta+1, kOrange-1,kCyan+2,kYellow+2};
const Int_t markers[]    = {kFullCircle, kFullSquare,kOpenCircle,kOpenSquare,kOpenDiamond,kOpenCross,kFullCross,kFullDiamond,kFullStar,kOpenStar};


//Options:
TFile* file_O2Analysis = new TFile("QAResults.root");
// TFile* file_O2Analysis = new TFile("AnalysisResults_noMC.root");
// TFile* file_O2Analysis = new TFile("AnalysisResultsMC.root");

const Int_t numPart = 3;
const Int_t ipart = 0;
const TString NamePart[numPart] = {"K0s", "Lambda", "AntiLambda"};//, "XiPlu", "XiMin", "OmPlu", "OmMin"};


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

TString* texDecayLength = new TString("Decay Length (cm)");
TString* texCTau = new TString("CTau (cm)");

void DecayPlots() {
  // Load necessary libraries
  LoadLibs();
  // Set the default style
  SetStyle();

  //histograms definition and processing
  TH1D* hstat, *hsyst, *hsystCorr;
  Int_t icolor=0;

  TString* SaveAs_Title = new TString("") ;
  TString* texXtitle = new TString("") ;
  TString* texYtitle = new TString("") ;

  // DecayLength(hstat, hsyst, hsystCorr, file_O2Analysis, SaveAs_Title, texXtitle, texYtitle);
  CTau(hstat, hsyst, hsystCorr, file_O2Analysis, SaveAs_Title, texXtitle, texYtitle);
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
  TH1 * h = cfig->DrawFrame(0,0,40,3000); // DecayLength

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
  TLatex * text = new TLatex (0.55,0.75,"p-p #sqrt{#it{s}_{NN}} = 0.900 TeV");
  text->SetNDC(kTRUE); //remove if I want x,y in TLatex to be in the coordinate system of the histogram
  text->Draw();
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

void DecayLength(TH1D* &hstat, TH1D* &hsyst, TH1D*&hsystCorr, TFile* file_O2Analysis, TString* SaveAs_Title, TString* &texXtitle, TString* &texYtitle) {

  //O2 file
  // TFile* file_O2Analysis = new TFile("AnalysisResults_WithSel8Cut.root");

  TH1D* H1_DecayLength = (TH1D*)file_O2Analysis->Get("v0cascades-q-a/histos-V0/DecayLength"+NamePart[numPart]);

  // TH1D* H1_DecayLength_Rebin = (TH1D*)H2_Armenteros->Rebin(5,"H2_Armenteros_RebinX");

  hstat = (TH1D*)H1_DecayLength->Clone("hstat");


  //////Error Bars///////
  hsyst     = new TH1D("hsyst", "hsyst", 100, 0, 10);
  hsystCorr = new TH1D("hsystCorr", "hsystCorr", 100, 0, 10);
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

  *SaveAs_Title += "DecayLength";
  texXtitle = texDecayLength;
  texYtitle = texCount;
}

void CTau(TH1D* &hstat, TH1D* &hsyst, TH1D*&hsystCorr, TFile* file_O2Analysis, TString* SaveAs_Title, TString* &texXtitle, TString* &texYtitle) {

  //O2 file
  // TFile* file_O2Analysis = new TFile("AnalysisResults_WithSel8Cut.root");

  TH1D* H1_CTau = (TH1D*)file_O2Analysis->Get("v0cascades-q-a/histos-V0/CTau"+NamePart[numPart]);

  hstat = (TH1D*)H1_CTau->Clone("hstat");


  //////Error Bars///////
  hsyst     = new TH1D("hsyst", "hsyst", 100, 0, 10);
  hsystCorr = new TH1D("hsystCorr", "hsystCorr", 100, 0, 10);
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

  *SaveAs_Title += "CTau";
  texXtitle = texCTau;
  texYtitle = texCount;
}
