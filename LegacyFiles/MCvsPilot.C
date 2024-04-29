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
void V0DCAV0Daughters(TH1D* &hstat_Data, TH1D* &hsyst_Data, TH1D*&hsystCorr_Data, TH1D* &hstat_MC, TH1D* &hsyst_MC, TH1D*&hsystCorr_MC, TFile* file_O2Analysis_Data, TFile* file_O2Analysis_MC, TString* SaveAs_Title, TString* &texXtitle, TString* &texYtitle);
void V0CosPA(TH1D* &hstat_Data, TH1D* &hsyst_Data, TH1D*&hsystCorr_Data, TH1D* &hstat_MC, TH1D* &hsyst_MC, TH1D*&hsystCorr_MC, TFile* file_O2Analysis_Data, TFile* file_O2Analysis_MC, TString* SaveAs_Title, TString* &texXtitle, TString* &texYtitle);
void DCAPosToPV(TH1D* &hstat_Data, TH1D* &hsyst_Data, TH1D*&hsystCorr_Data, TH1D* &hstat_MC, TH1D* &hsyst_MC, TH1D*&hsystCorr_MC, TFile* file_O2Analysis_Data, TFile* file_O2Analysis_MC, TString* SaveAs_Title, TString* &texXtitle, TString* &texYtitle);
void V0Radius(TH1D* &hstat_Data, TH1D* &hsyst_Data, TH1D*&hsystCorr_Data, TH1D* &hstat_MC, TH1D* &hsyst_MC, TH1D*&hsystCorr_MC, TFile* file_O2Analysis_Data, TFile* file_O2Analysis_MC, TString* SaveAs_Title, TString* &texXtitle, TString* &texYtitle);
void TpcNClsFound_PosTrack(TH1D* &hstat_Data, TH1D* &hsyst_Data, TH1D*&hsystCorr_Data, TH1D* &hstat_MC, TH1D* &hsyst_MC, TH1D*&hsystCorr_MC, TFile* file_O2Analysis_Data, TFile* file_O2Analysis_MC, TFile* file_QAResults_Data, TFile* file_QAResults_MC, TString* SaveAs_Title, TString* &texXtitle, TString* &texYtitle);
void DecayLength(TH1D* &hstat_Data, TH1D* &hsyst_Data, TH1D*&hsystCorr_Data, TH1D* &hstat_MC, TH1D* &hsyst_MC, TH1D*&hsystCorr_MC, TFile* file_O2Analysis_Data, TFile* file_O2Analysis_MC, TFile* file_QAResults_Data, TFile* file_QAResults_MC, TString* SaveAs_Title, TString* &texXtitle, TString* &texYtitle);
void Ctau(TH1D* &hstat_Data, TH1D* &hsyst_Data, TH1D*&hsystCorr_Data, TH1D* &hstat_MC, TH1D* &hsyst_MC, TH1D*&hsystCorr_MC, TFile* file_O2Analysis_Data, TFile* file_O2Analysis_MC, TFile* file_QAResults_Data, TFile* file_QAResults_MC, TString* SaveAs_Title, TString* &texXtitle, TString* &texYtitle);
void pTspectrum(TH1D* &hstat_Data, TH1D* &hsyst_Data, TH1D*&hsystCorr_Data, TH1D* &hstat_MC, TH1D* &hsyst_MC, TH1D*&hsystCorr_MC, TFile* file_O2Analysis_Data, TFile* file_O2Analysis_MC, TString* SaveAs_Title, TString* &texXtitle, TString* &texYtitle);

// Preferred colors and markers
const Int_t fillColors[] = {kGray+1,  kRed-10, kBlue-9, kGreen-8, kMagenta-9, kOrange-9,kCyan-8,kYellow-7}; // for syst bands
const Int_t colors[]     = {kBlack, kRed+1 , kBlue+1, kGreen+3, kMagenta+1, kOrange-1,kCyan+2,kYellow+2};
const Int_t markers[]    = {kFullCircle, kFullSquare,kOpenCircle,kOpenSquare,kOpenDiamond,kOpenCross,kFullCross,kFullDiamond,kFullStar,kOpenStar};


//Options:
TFile* file_O2Analysis_Data = new TFile("AnalysisResults_Data.root");
// TFile* file_O2Analysis = new TFile("AnalysisResults_noMC.root");
TFile* file_O2Analysis_MC = new TFile("AnalysisResults_MC.root");
TFile* file_QAResults_Data = new TFile("QAResults_Data.root");
TFile* file_QAResults_MC = new TFile("QAResults_MC.root");

const Int_t normalisation = 1; //0 for per event, 1 for per V0
const Int_t numPart = 3;
const Int_t ipart = 0;
const TString NamePart[numPart] = {"K0s", "Lambda", "AntiLambda"};//, "XiPlu", "XiMin", "OmPlu", "OmMin"};
const TString NamePart2[numPart] = {"K0Short", "Lambda", "AntiLambda"};//, "XiPlu", "XiMin", "OmPlu", "OmMin"};
const TString NamePart_Latex[numPart] = {"K^{0}_{S}", "#Lambda", "#bar{#Lambda}"};//, "XiPlu", "XiMin", "OmPlu", "OmMin"};


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
TString* texCtau = new TString("Ctau (cm)");


TString* texCountPerEvent = new TString("count per event");
TString* texCountPerV0 = new TString("count per V0");
TString* texDCA = new TString("DCA (cm^{2})");
TString* texCosPA = new TString("Cosine Pointing Angle");
TString* texV0Radius = new TString("V0 Radius (cm)");
TString* texRapidity = new TString("Rapidity");
TString* texTpcNClsFound = new TString("TPC Number of Clusters Found");

void MCvsPilot() {
  // Load necessary libraries
  LoadLibs();
  // Set the default style
  SetStyle();

  //histograms definition and processing
  TH1D* hstat_Data, *hsyst_Data, *hsystCorr_Data;
  TH1D* hstat_MC, *hsyst_MC, *hsystCorr_MC;
  Int_t icolor=0;
  Int_t icolor_Data=2;
  Int_t icolor_MC=1;

  TString* SaveAs_Title = new TString("") ;
  TString* texXtitle = new TString("") ;
  TString* texYtitle = new TString("") ;

  // V0DCAV0Daughters(hstat_Data, hsyst_Data, hsystCorr_Data, hstat_MC, hsyst_MC, hsystCorr_MC, file_O2Analysis_Data, file_O2Analysis_MC, SaveAs_Title, texXtitle, texYtitle);
  // V0CosPA(hstat_Data, hsyst_Data, hsystCorr_Data, hstat_MC, hsyst_MC, hsystCorr_MC, file_O2Analysis_Data, file_O2Analysis_MC, SaveAs_Title, texXtitle, texYtitle);
  // DCAPosToPV(hstat_Data, hsyst_Data, hsystCorr_Data, hstat_MC, hsyst_MC, hsystCorr_MC, file_O2Analysis_Data, file_O2Analysis_MC, SaveAs_Title, texXtitle, texYtitle);
  // V0Radius(hstat_Data, hsyst_Data, hsystCorr_Data, hstat_MC, hsyst_MC, hsystCorr_MC, file_O2Analysis_Data, file_O2Analysis_MC, SaveAs_Title, texXtitle, texYtitle);
  // TpcNClsFound_PosTrack(hstat_Data, hsyst_Data, hsystCorr_Data, hstat_MC, hsyst_MC, hsystCorr_MC, file_O2Analysis_Data, file_O2Analysis_MC, file_QAResults_Data, file_QAResults_MC, SaveAs_Title, texXtitle, texYtitle);
  // DecayLength(hstat_Data, hsyst_Data, hsystCorr_Data, hstat_MC, hsyst_MC, hsystCorr_MC, file_O2Analysis_Data, file_O2Analysis_MC, file_QAResults_Data, file_QAResults_MC, SaveAs_Title, texXtitle, texYtitle);
  // Ctau(hstat_Data, hsyst_Data, hsystCorr_Data, hstat_MC, hsyst_MC, hsystCorr_MC, file_O2Analysis_Data, file_O2Analysis_MC, file_QAResults_Data, file_QAResults_MC, SaveAs_Title, texXtitle, texYtitle);
  pTspectrum(hstat_Data, hsyst_Data, hsystCorr_Data, hstat_MC, hsyst_MC, hsystCorr_MC, file_O2Analysis_Data, file_O2Analysis_MC, SaveAs_Title, texXtitle, texYtitle);

  // RawSpectrum_O2data_truePt_trueV0s(hstat, hsyst, hsystCorr, ipart, file_O2Analysis, SaveAs_Title, texXtitle, texYtitle, pTbins[ipart], nbinpT[ipart]);

 // Prepare Figure, please stick to the default canvas size(s) unless absolutely necessary in your case
  // Rectangular
  TCanvas *cfig = new TCanvas("cfig", "Alice Figure Template", 800, 600); 
  // Square
  //TCanvas *cfig = new TCanvas("cfig", "Alice Figure Template", 800, 800); 
  cfig->SetLogy();

  cfig->cd();
  // Set Titles etc..
  // DrawFrame(xmin, ymin, xmax, ymax)
  // TH1 * h = cfig->DrawFrame(0,0.000001,1,0.05); // V0DCAV0Daughters per Event
  // TH1 * h = cfig->DrawFrame(0.99,0.000001,1,0.05); // CosPA per Event
  // TH1 * h = cfig->DrawFrame(-5,0.000001,5,0.005); // DCAPosToPV per Event
  // TH1 * h = cfig->DrawFrame(0,0.000001,5,0.05); // V0Radius per Event
  // TH1 * h = cfig->DrawFrame(0,0.000001,160,0.1); // TpcNClsFound_PosTrack per Event
  // TH1 * h = cfig->DrawFrame(0,0.000001,30,0.1); // DecayLength K0Short per Event
  // TH1 * h = cfig->DrawFrame(0,0.000001,15,0.1); // Ctau K0Short per Event

  // TH1 * h = cfig->DrawFrame(0,0.0001,1,10); // V0DCAV0Daughters per V0
  // TH1 * h = cfig->DrawFrame(0.99,0.0001,1,1); // CosPA per V0
  // TH1 * h = cfig->DrawFrame(-3,0.0001,3,0.1); // DCAPosToPV per V0
  // TH1 * h = cfig->DrawFrame(0,0.00001,10,5); // V0Radius per V0
  // TH1 * h = cfig->DrawFrame(0,0.0001,160,1); // TpcNClsFound_PosTrack per V0
  // TH1 * h = cfig->DrawFrame(0,0.00001,30,10); // DecayLength K0Short per V0
  // TH1 * h = cfig->DrawFrame(0,0.0001,15,10); // Ctau K0Short per V0
  TH1 * h = cfig->DrawFrame(0,0.0001,4,1); // pTspectrum K0Short per V0

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
  cout << "title: " << texXtitle->Data() <<endl;
  h->SetXTitle(texXtitle->Data());
  // // Please be consistent on the y label
  // // h->SetYTitle(texptDifferentialYield);//(texPtY);
  h->SetYTitle(texYtitle->Data());//(texPtY);

//   h->GetXaxis()->SetTitle(texXtitle);
//   h->GetYaxis()->SetTitle(texYtitle);

  // Draw your histos here:
  // Data
  hsystCorr_Data->SetFillColor(fillColors[icolor_Data]);
  hsyst_Data    ->SetFillColor(fillColors[icolor_Data]);
  // hsystCorr->Draw("E3,same"); //SHOULD BE PLOTTED EVENTUALLY
  hsyst_Data->SetFillStyle(0); // To draw empty boxes
  hsyst_Data->SetLineColor(colors[icolor_Data]); // To draw empty boxes
  hsystCorr_Data->SetLineColor(colors[icolor_Data]); // To draw empty boxes
  // hsyst_Data->Draw("E2,same"); //SHOULD BE PLOTTED EVENTUALLY
  hstat_Data->Draw("E,same");
  hstat_Data->SetMarkerStyle(markers[0]);
  // use the same color for markers and lines
  hstat_Data->SetMarkerColor(colors [icolor_Data]);
  hstat_Data->SetLineColor  (colors [icolor_Data]);

  // MC
  hsystCorr_MC->SetFillColor(fillColors[icolor_MC]);
  hsyst_MC    ->SetFillColor(fillColors[icolor_MC]);
  // hsystCorr->Draw("E3,same"); //SHOULD BE PLOTTED EVENTUALLY
  hsyst_MC->SetFillStyle(0); // To draw empty boxes
  hsyst_MC->SetLineColor(colors[icolor_MC]); // To draw empty boxes
  hsystCorr_MC->SetLineColor(colors[icolor_MC]); // To draw empty boxes
  // hsyst_MC->Draw("E2,same"); //SHOULD BE PLOTTED EVENTUALLY
  hstat_MC->Draw("E,same");
  hstat_MC->SetMarkerStyle(markers[0]);
  // use the same color for markers and lines
  hstat_MC->SetMarkerColor(colors [icolor_MC]);
  hstat_MC->SetLineColor  (colors [icolor_MC]);

  // Draw the logo   
  //  0: Just "ALICE" (for final data), to be added only if ALICE does not appear otherwise (e.g. in the legend)
  //  >0: ALICE Preliminary
  // DrawLogo(1, 0.59, 0.81);

  // You should always specify the colliding system
  // NOTATION: pp, p-Pb, Pb-Pb. 
  // Don't forget to use #sqrt{s_{NN}} for p-Pb and Pb-Pb
  // You can change the position of this with
  TLatex * textContext = new TLatex (0.18,0.82,"#bf{ALICE Run 3 Performance}"); //BOLD
  textContext->SetTextSize(0.05);
  textContext->SetNDC(kTRUE); //remove if I want x,y in TLatex to be in the coordinate system of the histogram
  textContext->Draw();
  TLatex * textColl = new TLatex (0.18,0.75,"pp #sqrt{#it{s}} = 900 GeV, #it{B_{z}} = 0.2 T");
  textColl->SetTextSize(0.04);
  textColl->SetNDC(kTRUE); //remove if I want x,y in TLatex to be in the coordinate system of the histogram
  textColl->Draw();
  TLatex * text_part = new TLatex (0.18,0.70,NamePart_Latex[ipart]);
  text_part->SetTextSize(0.04);
  text_part->SetNDC(kTRUE); //remove if I want x,y in TLatex to be in the coordinate system of the histogram
  text_part->Draw();
  // TLatex * text_extra = new TLatex (0.65,0.55,*Extra);
  // text_extra->SetNDC(kTRUE); //remove if I want x,y in TLatex to be in the coordinate system of the histogram
  // text_extra->Draw();
  // TLatex * text2 = new TLatex (0.55,55,"V0A Multiplicity Classes (Pb-Side)");
  // text2->SetTextSizePixels(24);
  // text2->Draw();
//  copy new layout of PaperPlot for comments and axis titles

  //Legend, if needed
  TLegend * leg = new TLegend(  0.19,  0.19,  0.3, 0.26); //x1 y1 x2 y2
  leg->AddEntry(hstat_Data,     "Data",   "LPE");
  leg->AddEntry(hstat_MC,     "MC",   "LPE");
  // leg->AddEntry(hsyst,     "syst error (Uncorrelated)",  "F");
  // leg->AddEntry(hsystCorr, "syst error (Correlated)",    "F" );
  leg->SetFillColor(0);
  leg->SetTextSize(gStyle->GetTextSize()*0.8);
  leg->Draw();

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
  
   
  if (normalisation == 0) {
    cfig->SaveAs(*SaveAs_Title+"_perEvent.pdf","pdf");
  }
  else if (normalisation == 1) {
    cfig->SaveAs(*SaveAs_Title+"_perV0.pdf","pdf");
  }
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

void V0DCAV0Daughters(TH1D* &hstat_Data, TH1D* &hsyst_Data, TH1D*&hsystCorr_Data, TH1D* &hstat_MC, TH1D* &hsyst_MC, TH1D*&hsystCorr_MC, TFile* file_O2Analysis_Data, TFile* file_O2Analysis_MC, TString* SaveAs_Title, TString* &texXtitle, TString* &texYtitle) {

  TH1I* H1I_SelectedEventCount_Data = (TH1I*)file_O2Analysis_Data->Get("lambdakzero-builder/hEventCounter");
  float SelectedEventCount_Data = H1I_SelectedEventCount_Data->GetEntries();
  TH1I* H1I_SelectedEventCount_MC = (TH1I*)file_O2Analysis_MC->Get("lambdakzero-builder/hEventCounter");
  float SelectedEventCount_MC = H1I_SelectedEventCount_MC->GetEntries();

  TH1D* H1_V0DCAV0Daughters_Data = (TH1D*)file_O2Analysis_Data->Get("lambdakzero-qa/hDCAV0Dau");
  TH1D* H1_V0DCAV0Daughters_MC = (TH1D*)file_O2Analysis_MC->Get("lambdakzero-qa/hDCAV0Dau");
  float V0count_Data = H1_V0DCAV0Daughters_Data->GetEntries();
  float V0count_MC = H1_V0DCAV0Daughters_MC->GetEntries();
  
  if (normalisation == 0) {
    H1_V0DCAV0Daughters_Data->Scale(1./SelectedEventCount_Data);
    H1_V0DCAV0Daughters_MC->Scale(1./SelectedEventCount_MC);
  }
  else if (normalisation == 1) {
    H1_V0DCAV0Daughters_Data->Scale(1./V0count_Data);
    H1_V0DCAV0Daughters_MC->Scale(1./V0count_MC);
  }
  
  hstat_Data = (TH1D*)H1_V0DCAV0Daughters_Data->Clone("hstat_Data");
  hstat_MC = (TH1D*)H1_V0DCAV0Daughters_MC->Clone("hstat_MC");

  // cout << "V0DCAV0Daughters max Data = " << hstat_Data->GetMaximum() << endl;
  // cout << "V0DCAV0Daughters max MC = " << hstat_MC->GetMaximum() << endl;

  //////Error Bars///////
  //Data
  hsyst_Data     = new TH1D("hsyst_Data", "hsyst_Data", 100, 0, 10);
  hsystCorr_Data = new TH1D("hsystCorr_Data", "hsystCorr_Data", 100, 0, 10);
  hstat_Data     ->Sumw2(); 
  hsyst_Data     ->Sumw2();
  hsystCorr_Data ->Sumw2(); 
  Int_t nbinx_Data = hstat_Data->GetNbinsX();
  //MC
  hsyst_MC     = new TH1D("hsyst_MC", "hsyst_MC", 100, 0, 10);
  hsystCorr_MC = new TH1D("hsystCorr_MC", "hsystCorr_MC", 100, 0, 10);
  hstat_MC     ->Sumw2(); 
  hsyst_MC     ->Sumw2();
  hsystCorr_MC ->Sumw2(); 
  Int_t nbinx_MC = hstat_MC->GetNbinsX();
  // for(Int_t ibinx = 1; ibinx <= nbinx; ibinx++){
  //   hsyst->SetBinContent (ibinx, hstat->GetBinContent(ibinx));
  //   hsyst->SetBinError   (ibinx, hstat->GetBinContent(ibinx)*0.08);
  //   hsystCorr->SetBinContent (ibinx, hstat->GetBinContent(ibinx));
  //   hsystCorr->SetBinError   (ibinx, hstat->GetBinContent(ibinx)*0.15);
  // } 

  *SaveAs_Title += "DCAV0Dau";
  texXtitle = texDCA;
  if (normalisation == 0) {
    texYtitle = texCountPerEvent;
  }
  else if (normalisation == 1) {
    texYtitle = texCountPerV0;
  }
}

void V0CosPA(TH1D* &hstat_Data, TH1D* &hsyst_Data, TH1D*&hsystCorr_Data, TH1D* &hstat_MC, TH1D* &hsyst_MC, TH1D*&hsystCorr_MC, TFile* file_O2Analysis_Data, TFile* file_O2Analysis_MC, TString* SaveAs_Title, TString* &texXtitle, TString* &texYtitle) {

  TH1I* H1I_SelectedEventCount_Data = (TH1I*)file_O2Analysis_Data->Get("lambdakzero-builder/hEventCounter");
  float SelectedEventCount_Data = H1I_SelectedEventCount_Data->GetEntries();
  TH1I* H1I_SelectedEventCount_MC = (TH1I*)file_O2Analysis_MC->Get("lambdakzero-builder/hEventCounter");
  float SelectedEventCount_MC = H1I_SelectedEventCount_MC->GetEntries();

  TH1D* H1_V0CosPA_Data = (TH1D*)file_O2Analysis_Data->Get("lambdakzero-qa/hV0CosPA");
  TH1D* H1_V0CosPA_MC = (TH1D*)file_O2Analysis_MC->Get("lambdakzero-qa/hV0CosPA");
  float V0count_Data = H1_V0CosPA_Data->GetEntries();
  float V0count_MC = H1_V0CosPA_MC->GetEntries();
  
  if (normalisation == 0) {
    H1_V0CosPA_Data->Scale(1./SelectedEventCount_Data);
    H1_V0CosPA_MC->Scale(1./SelectedEventCount_MC);
  }
  else if (normalisation == 1) {
    H1_V0CosPA_Data->Scale(1./V0count_Data);
    H1_V0CosPA_MC->Scale(1./V0count_MC);
  }

  hstat_Data = (TH1D*)H1_V0CosPA_Data->Clone("hstat_Data");
  hstat_MC = (TH1D*)H1_V0CosPA_MC->Clone("hstat_MC");


  //////Error Bars///////
  //Data
  hsyst_Data     = new TH1D("hsyst_Data", "hsyst_Data", 100, 0, 10);
  hsystCorr_Data = new TH1D("hsystCorr_Data", "hsystCorr_Data", 100, 0, 10);
  hstat_Data     ->Sumw2(); 
  hsyst_Data     ->Sumw2();
  hsystCorr_Data ->Sumw2(); 
  Int_t nbinx_Data = hstat_Data->GetNbinsX();
  //MC
  hsyst_MC     = new TH1D("hsyst_MC", "hsyst_MC", 100, 0, 10);
  hsystCorr_MC = new TH1D("hsystCorr_MC", "hsystCorr_MC", 100, 0, 10);
  hstat_MC     ->Sumw2(); 
  hsyst_MC     ->Sumw2();
  hsystCorr_MC ->Sumw2(); 
  Int_t nbinx_MC = hstat_MC->GetNbinsX();
  // for(Int_t ibinx = 1; ibinx <= nbinx; ibinx++){
  //   hsyst->SetBinContent (ibinx, hstat->GetBinContent(ibinx));
  //   hsyst->SetBinError   (ibinx, hstat->GetBinContent(ibinx)*0.08);
  //   hsystCorr->SetBinContent (ibinx, hstat->GetBinContent(ibinx));
  //   hsystCorr->SetBinError   (ibinx, hstat->GetBinContent(ibinx)*0.15);
  // } 

  *SaveAs_Title += "V0CosPA";
  texXtitle = texCosPA;
  if (normalisation == 0) {
    texYtitle = texCountPerEvent;
  }
  else if (normalisation == 1) {
    texYtitle = texCountPerV0;
  }
}

void DCAPosToPV(TH1D* &hstat_Data, TH1D* &hsyst_Data, TH1D*&hsystCorr_Data, TH1D* &hstat_MC, TH1D* &hsyst_MC, TH1D*&hsystCorr_MC, TFile* file_O2Analysis_Data, TFile* file_O2Analysis_MC, TString* SaveAs_Title, TString* &texXtitle, TString* &texYtitle) {

  TH1I* H1I_SelectedEventCount_Data = (TH1I*)file_O2Analysis_Data->Get("lambdakzero-builder/hEventCounter");
  float SelectedEventCount_Data = H1I_SelectedEventCount_Data->GetEntries();
  TH1I* H1I_SelectedEventCount_MC = (TH1I*)file_O2Analysis_MC->Get("lambdakzero-builder/hEventCounter");
  float SelectedEventCount_MC = H1I_SelectedEventCount_MC->GetEntries();

  TH1D* H1_DCAPosToPV_Data = (TH1D*)file_O2Analysis_Data->Get("lambdakzero-qa/hDCAPosToPV");
  TH1D* H1_DCAPosToPV_MC = (TH1D*)file_O2Analysis_MC->Get("lambdakzero-qa/hDCAPosToPV");
  float V0count_Data = H1_DCAPosToPV_Data->GetEntries();
  float V0count_MC = H1_DCAPosToPV_MC->GetEntries();
  
  if (normalisation == 0) {
    H1_DCAPosToPV_Data->Scale(1./SelectedEventCount_Data);
    H1_DCAPosToPV_MC->Scale(1./SelectedEventCount_MC);
  }
  else if (normalisation == 1) {
    H1_DCAPosToPV_Data->Scale(1./V0count_Data);
    H1_DCAPosToPV_MC->Scale(1./V0count_MC);
  }

  hstat_Data = (TH1D*)H1_DCAPosToPV_Data->Clone("hstat_Data");
  hstat_MC = (TH1D*)H1_DCAPosToPV_MC->Clone("hstat_MC");


  //////Error Bars///////
  //Data
  hsyst_Data     = new TH1D("hsyst_Data", "hsyst_Data", 100, 0, 10);
  hsystCorr_Data = new TH1D("hsystCorr_Data", "hsystCorr_Data", 100, 0, 10);
  hstat_Data     ->Sumw2(); 
  hsyst_Data     ->Sumw2();
  hsystCorr_Data ->Sumw2(); 
  Int_t nbinx_Data = hstat_Data->GetNbinsX();
  //MC
  hsyst_MC     = new TH1D("hsyst_MC", "hsyst_MC", 100, 0, 10);
  hsystCorr_MC = new TH1D("hsystCorr_MC", "hsystCorr_MC", 100, 0, 10);
  hstat_MC     ->Sumw2(); 
  hsyst_MC     ->Sumw2();
  hsystCorr_MC ->Sumw2(); 
  Int_t nbinx_MC = hstat_MC->GetNbinsX();
  // for(Int_t ibinx = 1; ibinx <= nbinx; ibinx++){
  //   hsyst->SetBinContent (ibinx, hstat->GetBinContent(ibinx));
  //   hsyst->SetBinError   (ibinx, hstat->GetBinContent(ibinx)*0.08);
  //   hsystCorr->SetBinContent (ibinx, hstat->GetBinContent(ibinx));
  //   hsystCorr->SetBinError   (ibinx, hstat->GetBinContent(ibinx)*0.15);
  // } 

  *SaveAs_Title += "DCAPosToPV";
  texXtitle = texDCA;
  if (normalisation == 0) {
    texYtitle = texCountPerEvent;
  }
  else if (normalisation == 1) {
    texYtitle = texCountPerV0;
  }
}


void V0Radius(TH1D* &hstat_Data, TH1D* &hsyst_Data, TH1D*&hsystCorr_Data, TH1D* &hstat_MC, TH1D* &hsyst_MC, TH1D*&hsystCorr_MC, TFile* file_O2Analysis_Data, TFile* file_O2Analysis_MC, TString* SaveAs_Title, TString* &texXtitle, TString* &texYtitle) {

  TH1I* H1I_SelectedEventCount_Data = (TH1I*)file_O2Analysis_Data->Get("lambdakzero-builder/hEventCounter");
  float SelectedEventCount_Data = H1I_SelectedEventCount_Data->GetEntries();
  TH1I* H1I_SelectedEventCount_MC = (TH1I*)file_O2Analysis_MC->Get("lambdakzero-builder/hEventCounter");
  float SelectedEventCount_MC = H1I_SelectedEventCount_MC->GetEntries();

  TH1D* H1_V0Radius_Data = (TH1D*)file_O2Analysis_Data->Get("lambdakzero-qa/hV0Radius");
  TH1D* H1_V0Radius_MC = (TH1D*)file_O2Analysis_MC->Get("lambdakzero-qa/hV0Radius");
  float V0count_Data = H1_V0Radius_Data->GetEntries();
  float V0count_MC = H1_V0Radius_MC->GetEntries();
  
  if (normalisation == 0) {
    H1_V0Radius_Data->Scale(1./SelectedEventCount_Data);
    H1_V0Radius_MC->Scale(1./SelectedEventCount_MC);
  }
  else if (normalisation == 1) {
    H1_V0Radius_Data->Scale(1./V0count_Data);
    H1_V0Radius_MC->Scale(1./V0count_MC);
  }

  hstat_Data = (TH1D*)H1_V0Radius_Data->Clone("hstat_Data");
  hstat_MC = (TH1D*)H1_V0Radius_MC->Clone("hstat_MC");


  //////Error Bars///////
  //Data
  hsyst_Data     = new TH1D("hsyst_Data", "hsyst_Data", 100, 0, 10);
  hsystCorr_Data = new TH1D("hsystCorr_Data", "hsystCorr_Data", 100, 0, 10);
  hstat_Data     ->Sumw2(); 
  hsyst_Data     ->Sumw2();
  hsystCorr_Data ->Sumw2(); 
  Int_t nbinx_Data = hstat_Data->GetNbinsX();
  //MC
  hsyst_MC     = new TH1D("hsyst_MC", "hsyst_MC", 100, 0, 10);
  hsystCorr_MC = new TH1D("hsystCorr_MC", "hsystCorr_MC", 100, 0, 10);
  hstat_MC     ->Sumw2(); 
  hsyst_MC     ->Sumw2();
  hsystCorr_MC ->Sumw2(); 
  Int_t nbinx_MC = hstat_MC->GetNbinsX();
  // for(Int_t ibinx = 1; ibinx <= nbinx; ibinx++){
  //   hsyst->SetBinContent (ibinx, hstat->GetBinContent(ibinx));
  //   hsyst->SetBinError   (ibinx, hstat->GetBinContent(ibinx)*0.08);
  //   hsystCorr->SetBinContent (ibinx, hstat->GetBinContent(ibinx));
  //   hsystCorr->SetBinError   (ibinx, hstat->GetBinContent(ibinx)*0.15);
  // } 

  *SaveAs_Title += "V0Radius";
  texXtitle = texDCA;
  if (normalisation == 0) {
    texYtitle = texCountPerEvent;
  }
  else if (normalisation == 1) {
    texYtitle = texCountPerV0;
  }
}

void TpcNClsFound_PosTrack(TH1D* &hstat_Data, TH1D* &hsyst_Data, TH1D*&hsystCorr_Data, TH1D* &hstat_MC, TH1D* &hsyst_MC, TH1D*&hsystCorr_MC, TFile* file_O2Analysis_Data, TFile* file_O2Analysis_MC, TFile* file_QAResults_Data, TFile* file_QAResults_MC, TString* SaveAs_Title, TString* &texXtitle, TString* &texYtitle) {

  TH1I* H1I_SelectedEventCount_Data = (TH1I*)file_O2Analysis_Data->Get("lambdakzero-builder/hEventCounter");
  float SelectedEventCount_Data = H1I_SelectedEventCount_Data->GetEntries();
  TH1I* H1I_SelectedEventCount_MC = (TH1I*)file_O2Analysis_MC->Get("lambdakzero-builder/hEventCounter");
  float SelectedEventCount_MC = H1I_SelectedEventCount_MC->GetEntries();

  TH1D* H1_TpcNClsFound_PosTrack_Data = (TH1D*)file_QAResults_Data->Get("v0cascades-q-a/histos-V0/TpcNClsFound_Pos_preBuilder");
  TH1D* H1_TpcNClsFound_PosTrack_MC = (TH1D*)file_QAResults_MC->Get("v0cascades-q-a/histos-V0/TpcNClsFound_Pos_preBuilder");
  float V0count_Data = H1_TpcNClsFound_PosTrack_Data->GetEntries();
  float V0count_MC = H1_TpcNClsFound_PosTrack_MC->GetEntries();
  
  if (normalisation == 0) {
    H1_TpcNClsFound_PosTrack_Data->Scale(1./SelectedEventCount_Data);
    H1_TpcNClsFound_PosTrack_MC->Scale(1./SelectedEventCount_MC);
  }
  else if (normalisation == 1) {
    H1_TpcNClsFound_PosTrack_Data->Scale(1./V0count_Data);
    H1_TpcNClsFound_PosTrack_MC->Scale(1./V0count_MC);
  }

  hstat_Data = (TH1D*)H1_TpcNClsFound_PosTrack_Data->Clone("hstat_Data");
  hstat_MC = (TH1D*)H1_TpcNClsFound_PosTrack_MC->Clone("hstat_MC");


  //////Error Bars///////
  //Data
  hsyst_Data     = new TH1D("hsyst_Data", "hsyst_Data", 100, 0, 10);
  hsystCorr_Data = new TH1D("hsystCorr_Data", "hsystCorr_Data", 100, 0, 10);
  hstat_Data     ->Sumw2(); 
  hsyst_Data     ->Sumw2();
  hsystCorr_Data ->Sumw2(); 
  Int_t nbinx_Data = hstat_Data->GetNbinsX();
  //MC
  hsyst_MC     = new TH1D("hsyst_MC", "hsyst_MC", 100, 0, 10);
  hsystCorr_MC = new TH1D("hsystCorr_MC", "hsystCorr_MC", 100, 0, 10);
  hstat_MC     ->Sumw2(); 
  hsyst_MC     ->Sumw2();
  hsystCorr_MC ->Sumw2(); 
  Int_t nbinx_MC = hstat_MC->GetNbinsX();
  // for(Int_t ibinx = 1; ibinx <= nbinx; ibinx++){
  //   hsyst->SetBinContent (ibinx, hstat->GetBinContent(ibinx));
  //   hsyst->SetBinError   (ibinx, hstat->GetBinContent(ibinx)*0.08);
  //   hsystCorr->SetBinContent (ibinx, hstat->GetBinContent(ibinx));
  //   hsystCorr->SetBinError   (ibinx, hstat->GetBinContent(ibinx)*0.15);
  // } 

  *SaveAs_Title += "TpcNClsFound_PosTrack";
  texXtitle = texTpcNClsFound;
  if (normalisation == 0) {
    texYtitle = texCountPerEvent;
  }
  else if (normalisation == 1) {
    texYtitle = texCountPerV0;
  }
}


void DecayLength(TH1D* &hstat_Data, TH1D* &hsyst_Data, TH1D*&hsystCorr_Data, TH1D* &hstat_MC, TH1D* &hsyst_MC, TH1D*&hsystCorr_MC, TFile* file_O2Analysis_Data, TFile* file_O2Analysis_MC, TFile* file_QAResults_Data, TFile* file_QAResults_MC, TString* SaveAs_Title, TString* &texXtitle, TString* &texYtitle) {

  TH1I* H1I_SelectedEventCount_Data = (TH1I*)file_O2Analysis_Data->Get("lambdakzero-builder/hEventCounter");
  float SelectedEventCount_Data = H1I_SelectedEventCount_Data->GetEntries();
  TH1I* H1I_SelectedEventCount_MC = (TH1I*)file_O2Analysis_MC->Get("lambdakzero-builder/hEventCounter");
  float SelectedEventCount_MC = H1I_SelectedEventCount_MC->GetEntries();

  TH1D* H1_DecayLength_Data = (TH1D*)file_QAResults_Data->Get("v0cascades-q-a/histos-V0/DecayLength"+NamePart[ipart]);
  TH1D* H1_DecayLength_MC = (TH1D*)file_QAResults_MC->Get("v0cascades-q-a/histos-V0/DecayLength"+NamePart[ipart]);
  float V0count_Data = H1_DecayLength_Data->GetEntries();
  float V0count_MC = H1_DecayLength_MC->GetEntries();
  
  if (normalisation == 0) {
    H1_DecayLength_Data->Scale(1./SelectedEventCount_Data);
    H1_DecayLength_MC->Scale(1./SelectedEventCount_MC);
  }
  else if (normalisation == 1) {
    H1_DecayLength_Data->Scale(1./V0count_Data);
    H1_DecayLength_MC->Scale(1./V0count_MC);
  }

  hstat_Data = (TH1D*)H1_DecayLength_Data->Clone("hstat_Data");
  hstat_MC = (TH1D*)H1_DecayLength_MC->Clone("hstat_MC");


  //////Error Bars///////
  //Data
  hsyst_Data     = new TH1D("hsyst_Data", "hsyst_Data", 100, 0, 10);
  hsystCorr_Data = new TH1D("hsystCorr_Data", "hsystCorr_Data", 100, 0, 10);
  hstat_Data     ->Sumw2(); 
  hsyst_Data     ->Sumw2();
  hsystCorr_Data ->Sumw2(); 
  Int_t nbinx_Data = hstat_Data->GetNbinsX();
  //MC
  hsyst_MC     = new TH1D("hsyst_MC", "hsyst_MC", 100, 0, 10);
  hsystCorr_MC = new TH1D("hsystCorr_MC", "hsystCorr_MC", 100, 0, 10);
  hstat_MC     ->Sumw2(); 
  hsyst_MC     ->Sumw2();
  hsystCorr_MC ->Sumw2(); 
  Int_t nbinx_MC = hstat_MC->GetNbinsX();
  // for(Int_t ibinx = 1; ibinx <= nbinx; ibinx++){
  //   hsyst->SetBinContent (ibinx, hstat->GetBinContent(ibinx));
  //   hsyst->SetBinError   (ibinx, hstat->GetBinContent(ibinx)*0.08);
  //   hsystCorr->SetBinContent (ibinx, hstat->GetBinContent(ibinx));
  //   hsystCorr->SetBinError   (ibinx, hstat->GetBinContent(ibinx)*0.15);
  // } 

  *SaveAs_Title += "DecayLength";
  texXtitle = texDecayLength;
  if (normalisation == 0) {
    texYtitle = texCountPerEvent;
  }
  else if (normalisation == 1) {
    texYtitle = texCountPerV0;
  }
}


void Ctau(TH1D* &hstat_Data, TH1D* &hsyst_Data, TH1D*&hsystCorr_Data, TH1D* &hstat_MC, TH1D* &hsyst_MC, TH1D*&hsystCorr_MC, TFile* file_O2Analysis_Data, TFile* file_O2Analysis_MC, TFile* file_QAResults_Data, TFile* file_QAResults_MC, TString* SaveAs_Title, TString* &texXtitle, TString* &texYtitle) {

  TH1I* H1I_SelectedEventCount_Data = (TH1I*)file_O2Analysis_Data->Get("lambdakzero-builder/hEventCounter");
  float SelectedEventCount_Data = H1I_SelectedEventCount_Data->GetEntries();
  TH1I* H1I_SelectedEventCount_MC = (TH1I*)file_O2Analysis_MC->Get("lambdakzero-builder/hEventCounter");
  float SelectedEventCount_MC = H1I_SelectedEventCount_MC->GetEntries();

  TH1D* H1_Ctau_Data = (TH1D*)file_QAResults_Data->Get("v0cascades-q-a/histos-V0/Ctau"+NamePart[ipart]);
  TH1D* H1_Ctau_MC = (TH1D*)file_QAResults_MC->Get("v0cascades-q-a/histos-V0/Ctau"+NamePart[ipart]);
  float V0count_Data = H1_Ctau_Data->GetEntries();
  float V0count_MC = H1_Ctau_MC->GetEntries();
  
  if (normalisation == 0) {
    H1_Ctau_Data->Scale(1./SelectedEventCount_Data);
    H1_Ctau_MC->Scale(1./SelectedEventCount_MC);
  }
  else if (normalisation == 1) {
    H1_Ctau_Data->Scale(1./V0count_Data);
    H1_Ctau_MC->Scale(1./V0count_MC);
  }

  hstat_Data = (TH1D*)H1_Ctau_Data->Clone("hstat_Data");
  hstat_MC = (TH1D*)H1_Ctau_MC->Clone("hstat_MC");


  //////Error Bars///////
  //Data
  hsyst_Data     = new TH1D("hsyst_Data", "hsyst_Data", 100, 0, 10);
  hsystCorr_Data = new TH1D("hsystCorr_Data", "hsystCorr_Data", 100, 0, 10);
  hstat_Data     ->Sumw2(); 
  hsyst_Data     ->Sumw2();
  hsystCorr_Data ->Sumw2(); 
  Int_t nbinx_Data = hstat_Data->GetNbinsX();
  //MC
  hsyst_MC     = new TH1D("hsyst_MC", "hsyst_MC", 100, 0, 10);
  hsystCorr_MC = new TH1D("hsystCorr_MC", "hsystCorr_MC", 100, 0, 10);
  hstat_MC     ->Sumw2(); 
  hsyst_MC     ->Sumw2();
  hsystCorr_MC ->Sumw2(); 
  Int_t nbinx_MC = hstat_MC->GetNbinsX();
  // for(Int_t ibinx = 1; ibinx <= nbinx; ibinx++){
  //   hsyst->SetBinContent (ibinx, hstat->GetBinContent(ibinx));
  //   hsyst->SetBinError   (ibinx, hstat->GetBinContent(ibinx)*0.08);
  //   hsystCorr->SetBinContent (ibinx, hstat->GetBinContent(ibinx));
  //   hsystCorr->SetBinError   (ibinx, hstat->GetBinContent(ibinx)*0.15);
  // } 

  *SaveAs_Title += "Ctau";
  texXtitle = texCtau;
  if (normalisation == 0) {
    texYtitle = texCountPerEvent;
  }
  else if (normalisation == 1) {
    texYtitle = texCountPerV0;
  }
}

void pTspectrum(TH1D* &hstat_Data, TH1D* &hsyst_Data, TH1D*&hsystCorr_Data, TH1D* &hstat_MC, TH1D* &hsyst_MC, TH1D*&hsystCorr_MC, TFile* file_O2Analysis_Data, TFile* file_O2Analysis_MC, TString* SaveAs_Title, TString* &texXtitle, TString* &texYtitle) {

  TH1I* H1I_SelectedEventCount_Data = (TH1I*)file_O2Analysis_Data->Get("lambdakzero-builder/hEventCounter");
  float SelectedEventCount_Data = H1I_SelectedEventCount_Data->GetEntries();
  TH1I* H1I_SelectedEventCount_MC = (TH1I*)file_O2Analysis_MC->Get("lambdakzero-builder/hEventCounter");
  float SelectedEventCount_MC = H1I_SelectedEventCount_MC->GetEntries();

  TH3D* H3_pTspectrum_Data = (TH3D*)file_O2Analysis_Data->Get("lambdakzero-analysis/h3dMass"+NamePart2[ipart]);
  TH3D* H3_pTspectrum_MC = (TH3D*)file_O2Analysis_MC->Get("lambdakzero-analysis-mc/h3dMass"+NamePart2[ipart]);
  TH1D* H1_pTspectrum_Data = (TH1D*)H3_pTspectrum_Data->ProjectionY("pTspectrum_Data",0,-1,0,-1);
  TH1D* H1_pTspectrum_MC = (TH1D*)H3_pTspectrum_MC->ProjectionY("pTspectrum_MC",0,-1,0,-1);

  float V0count_Data = H1_pTspectrum_Data->GetEntries();
  float V0count_MC = H1_pTspectrum_MC->GetEntries();
  
  if (normalisation == 0) {
    H1_pTspectrum_Data->Scale(1./SelectedEventCount_Data);
    H1_pTspectrum_MC->Scale(1./SelectedEventCount_MC);
  }
  else if (normalisation == 1) {
    H1_pTspectrum_Data->Scale(1./V0count_Data);
    H1_pTspectrum_MC->Scale(1./V0count_MC);
  }

  hstat_Data = (TH1D*)H1_pTspectrum_Data->Clone("hstat_Data");
  hstat_MC = (TH1D*)H1_pTspectrum_MC->Clone("hstat_MC");


  //////Error Bars///////
  //Data
  hsyst_Data     = new TH1D("hsyst_Data", "hsyst_Data", 100, 0, 10);
  hsystCorr_Data = new TH1D("hsystCorr_Data", "hsystCorr_Data", 100, 0, 10);
  hstat_Data     ->Sumw2(); 
  hsyst_Data     ->Sumw2();
  hsystCorr_Data ->Sumw2(); 
  Int_t nbinx_Data = hstat_Data->GetNbinsX();
  //MC
  hsyst_MC     = new TH1D("hsyst_MC", "hsyst_MC", 100, 0, 10);
  hsystCorr_MC = new TH1D("hsystCorr_MC", "hsystCorr_MC", 100, 0, 10);
  hstat_MC     ->Sumw2(); 
  hsyst_MC     ->Sumw2();
  hsystCorr_MC ->Sumw2(); 
  Int_t nbinx_MC = hstat_MC->GetNbinsX();
  // for(Int_t ibinx = 1; ibinx <= nbinx; ibinx++){
  //   hsyst->SetBinContent (ibinx, hstat->GetBinContent(ibinx));
  //   hsyst->SetBinError   (ibinx, hstat->GetBinContent(ibinx)*0.08);
  //   hsystCorr->SetBinContent (ibinx, hstat->GetBinContent(ibinx));
  //   hsystCorr->SetBinError   (ibinx, hstat->GetBinContent(ibinx)*0.15);
  // } 

  *SaveAs_Title += "pTspectrum";
  texXtitle = texPtX;
  if (normalisation == 0) {
    texYtitle = texCountPerEvent;
  }
  else if (normalisation == 1) {
    texYtitle = texCountPerV0;
  }
}


