#include "HistogramUtilities.h"
#include "../Settings/GlobalSettings.h"
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <sstream>
#include<array>

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////// Various Utilities /////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

float findMinFloat(float* array, int length){
  float min = array[0];
  for (int i = 1; i < length; i++)
      if (array[i] < min)
          min = array[i];
  return min;
}
float findMaxFloat(float* array, int length){
  float max = array[0];
  for (int i = 1; i < length; i++)
      if (array[i] > max)
          max = array[i];
  return max;
}


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////// Histogram Operations //////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

std::vector<double> GetTH1Bins(TH1* H1_histo) {
  std::vector<double> bins;
  for(int iBin = 1; iBin <= H1_histo->GetNbinsX(); iBin++){
    bins.push_back(H1_histo->GetBinLowEdge(iBin));
    // cout << "ibin " << iBin << ": lowEdge = " << H1_histo->GetBinLowEdge(iBin) << endl;
  }
  bins.push_back(H1_histo->GetXaxis()->GetBinUpEdge(H1_histo->GetNbinsX()));
  // cout << "ibin " << H1_histo->GetNbinsX() << ": uppEdge = " << H1_histo->GetXaxis()->GetBinUpEdge(H1_histo->GetNbinsX()) << endl;
  return bins;
}

TH2D RebinVariableBins2D(TH2D* H2D_hist, int nBinsX, int nBinsY, double* binsX, double* binsY){

  TH2D H2D_hist_rebinned("H2D_hist_rebinned", "H2D_hist_rebinned", nBinsX, binsX, nBinsY, binsY);

  for(int iBinX = 1; iBinX <= H2D_hist->GetNbinsX(); iBinX++){
    for(int iBinY = 1; iBinY <= H2D_hist->GetNbinsY(); iBinY++){
      H2D_hist_rebinned.Fill(H2D_hist->GetXaxis()->GetBinCenter(iBinX), H2D_hist->GetYaxis()->GetBinCenter(iBinY), H2D_hist->GetBinContent(iBinX, iBinY));
    }
  }
  cout << "RebinVariableBins2D - what of the errors" << endl;
  cout << "H2D_hist_rebinned nbins X =" << H2D_hist_rebinned.GetNbinsX() << endl;

  std::stringstream ss;
  ss << H2D_hist->GetName() << "_RebinVariableBins2D";
  TString histName((TString)ss.str());

  // H2D_hist->Reset("M");
  // H2D_hist = (TH2D*)(&H2D_hist_rebinned)->Clone(histName);
  // cout << "H2D_hist nbins X =" << H2D_hist->GetNbinsX() << endl;
  return H2D_hist_rebinned;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////// Histogram Context /////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

TString contextCustomThreeFields(TString mainContext, TString secondaryContext, TString tertiaryContext, const char options[]){
  TString texContextFinal;
  texContextFinal = "#splitline{"+mainContext+" "+secondaryContext+"}{#splitline{2023 QC}{"+tertiaryContext+"}}";
  // texContextFinal = "#splitline{"+mainContext+" "+secondaryContext+"}{#splitline{2023 QC}{test"+tertiaryContext+"}}";
  // texContextFinal = "testtesttestest";
  return texContextFinal;
}

TString contextCustomTwoFields(TString mainContext, TString secondaryContext, const char options[]){
  return contextCustomThreeFields(mainContext, (TString)"" , secondaryContext, options);
}

TString contextCustomOneField(TString mainContext, const char options[]){
  return contextCustomTwoFields(mainContext, (TString)"", options);
}

TString contextPtRange(float* PtRange){
  std::stringstream ss;
  ss.setf(std::ios::fixed);
  ss.precision(2);
  ss << PtRange[0] << " < #it{p}_{T} < ";
  ss.precision(0);
  ss << PtRange[1];
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






TString contextDatasetRadiusCompAndVarRange(TString* mainContext, int iDataset, float* variableRange, const char options[]){
  TString texcontextDatasetRadiusCompAndVarRange;
  if (strstr(options, "pt") != NULL) { //  || strstr(options, "ratio") != NULL not sure why I had this here
    texcontextDatasetRadiusCompAndVarRange = "#splitline{"+*mainContext+" "+DatasetsNames[iDataset]+"}{#splitline{2023 QC}{"+contextPtRange(variableRange)+"}}";
  }
  if (strstr(options, "eta") != NULL) { //  || strstr(options, "ratio") != NULL not sure why I had this here
    texcontextDatasetRadiusCompAndVarRange = "#splitline{"+*mainContext+" "+DatasetsNames[iDataset]+"}{#splitline{2023 QC}{"+contextEtaRange(variableRange)+"}}";
  }

  return texcontextDatasetRadiusCompAndVarRange;
}

TString contextDatasetCompAndRadiusAndVarRange(TString* mainContext, float jetRadius, float* variableRange, const char options[]){
  TString texcontextDatasetCompAndRadiusAndVarRange;
  if (strstr(options, "pt") != NULL) { //  || strstr(options, "ratio") != NULL not sure why I had this here
    texcontextDatasetCompAndRadiusAndVarRange = "#splitline{"+*mainContext+"}{#splitline{"+contextJetRadius(jetRadius)+"}{"+contextPtRange(variableRange)+"}}";
  }
  if (strstr(options, "eta") != NULL) { //  || strstr(options, "ratio") != NULL not sure why I had this here
    texcontextDatasetCompAndRadiusAndVarRange = "#splitline{"+*mainContext+"}{#splitline{"+contextJetRadius(jetRadius)+"}{"+contextEtaRange(variableRange)+"}}";
  }

  return texcontextDatasetCompAndRadiusAndVarRange;
}

TString contextDatasetCompAndRadius(TString* mainContext, float jetRadius, const char options[]){
  TString texcontextDatasetCompAndRadius;
  texcontextDatasetCompAndRadius = "#splitline{"+*mainContext+"}{"+contextJetRadius(jetRadius)+"}";

  return texcontextDatasetCompAndRadius;
}

TString contextDatasetComp(TString* mainContext, const char options[]){
  TString texcontextDatasetComp;
  texcontextDatasetComp = *mainContext;

  return texcontextDatasetComp;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////// Histogram Drawing /////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


void Draw_TH1_Histograms_in_one(TH1D** histograms_collection, const TString* legendList_string, int collectionSize, TString Context, TString* pdfName, TString* &texXtitle, TString* &texYtitle, TString* texCollisionDataInfo, std::array<std::array<float, 2>, 2> drawnWindow, const char options[], TF1** optionalFitCollection) {
  // has options:
  // - "autoratio" : if in the options string, the Y range is chosen automatically based on the difference to 1
  // - "standardratio" : if in the options string, the Y range is [0,2.2]
  // - "logy" : if in the options string, then the y axis of the plot is set to a log scale (except for the ratio plot)
  // - "avoidFirst" : if in the options string, then the first histogram of the collection isn't plotted

  // canvas settings
  TCanvas *canvas = new TCanvas ("canvas"+*pdfName, "canvas"+*pdfName, 800, 800);
  canvas->cd(0);

  float minY_collection[collectionSize];
  float maxY_collection[collectionSize];
  float minX_collection[collectionSize];
  float maxX_collection[collectionSize];

  // cout << "test1" << endl;
  // char* ratioOptionCheck = "ratio";
  // char* OptionPtr = &options;

  for (int i = 0; i < collectionSize; i++) {
    if (strstr(options, "autoXrange") != NULL) {
      maxX_collection[i] = histograms_collection[i]->FindLastBinAbove(0, 1);
      minX_collection[i] = histograms_collection[i]->FindFirstBinAbove(0,1);
      // cout << "test1.1a" << endl;
    }
    else {
      maxX_collection[i] = histograms_collection[i]->GetXaxis()->GetXmax();
      minX_collection[i] = histograms_collection[i]->GetXaxis()->GetXmin();
      // cout << "test1.1b" << endl;
    }
    maxY_collection[i] = histograms_collection[i]->GetMaximum();
    // cout << "test1.2" << endl;
    
    if (strstr(options, "logy") != NULL) { //  || strstr(options, "ratio") != NULL not sure why I had this here
      minY_collection[i] = histograms_collection[i]->GetMinimum(1E-10);
      // cout << "test1.3a" << endl;
    }
    else if (strstr(options, "minYnotZero") != NULL) {
      minY_collection[i] = histograms_collection[i]->GetMinimum(GLOBAL_epsilon); // asks for the first/last bin on the y axis (axis number 2) to have strictly more than 1 entry)
      // cout << "test1.3b" << endl;
    }
    else {
      minY_collection[i] = histograms_collection[i]->GetMinimum();
      // cout << "test1.3c" << endl;
    }
  }
  // cout << "test2" << endl;

  float yUpMarginScaling, yDownMarginScaling;
  float maxX, minX, maxY, minY;
  if (std::equal(std::begin(drawnWindow), std::end(drawnWindow), std::begin(drawnWindowAuto), std::end(drawnWindowAuto))) {
    yUpMarginScaling = 1.4;
    yDownMarginScaling = 1;
    maxX = findMaxFloat(maxX_collection, collectionSize);
    minX = findMinFloat(minX_collection, collectionSize);
    maxY = findMaxFloat(maxY_collection, collectionSize);
    minY = findMinFloat(minY_collection, collectionSize);

    if (strstr(options, "logy") != NULL) {
      yUpMarginScaling = 100;
      if (minY < 0) {
      } 
      else if (minY < 1E-10) { // if minY is 0 logY won't like it
        minY = minY+1E-10 ;
      }
      else {
        minY = minY-1E-10 ; // to be sure to see the point?
      }
    } else {
      if (minY < 0 || strstr(options, "minYnotZero") != NULL) {
        minY > 0 ? yDownMarginScaling = 0.95 : yDownMarginScaling = 1.05;
        yUpMarginScaling = 1.1;
      } else {
        minY = 0.;
      }
    // cout << "test3" << endl;

      if (strstr(options, "autoratio") != NULL) {
        float deltaMax = max(1-minY, maxY-1);
        minY = max((float)0., 1-deltaMax);
        maxY = 1+deltaMax;
        yUpMarginScaling = 1.3;
      }
      if (strstr(options, "standardratio") != NULL) {
        minY = 0;
        maxY = 2;
        yUpMarginScaling = 1.1;
      }
      if (strstr(options, "efficiency") != NULL) {
        minY = 0;
        maxY = 1.5;
        yUpMarginScaling = 1.1;
      }
    }
  } else { // if not drawnWindowAuto, then set the window to match the one entered in drawnWindow
    yUpMarginScaling = 1;
    yDownMarginScaling = 1;
    minX = drawnWindow[0][0];
    maxX = drawnWindow[0][1];
    minY = drawnWindow[1][0];
    maxY = drawnWindow[1][1];
  }




  TH1 *hFrame = canvas->DrawFrame(minX, yDownMarginScaling*minY, maxX, yUpMarginScaling*maxY);
  if (strstr(options, "logy") != NULL) {
    canvas->SetLogy();
  }  
  if (strstr(options, "logx") != NULL) {
    canvas->SetLogx();
  }
    // cout << "test4" << endl;

  hFrame->SetXTitle(texXtitle->Data());
  hFrame->SetYTitle(texYtitle->Data());
  // hFrame->GetYaxis()->SetTitleOffset(2);

  // legend settings
  TLegend * leg = new TLegend(0.7, 0.75, 0.8, 0.87);
  leg->SetTextSize(gStyle->GetTextSize()*0.7);
  if (collectionSize >= 6) { // maybe fine tune that
    leg->SetTextSize(gStyle->GetTextSize()*0.3);
  }
  if (strstr(options, "fit") != NULL) {
    leg->SetTextSize(gStyle->GetTextSize()*0.3);
  }
  if (collectionSize >= 6) {
    gStyle->SetPalette(kRainbow); // for the choice of marker's colours; only use this if we have many histograms in the same plot
  }
  // cout << "test5" << endl;

  // draws histograms from collection
  for (int i = 0; i < collectionSize; i++) {
    if (i!=0 || strstr(options, "avoidFirst") == NULL) { // if i=0 requires that the option avoidFirst isn't there
      if (collectionSize >= 6) {
        histograms_collection[i]->Draw("same PMC PLC"); // PMC uses the palette chosen with gStyle->SetPalette() to chose the colours of the markers, PLC for the lines
      }
      else {
        histograms_collection[i]->Draw("same");
        histograms_collection[i]->SetMarkerColor(colors[i]);
        histograms_collection[i]->SetLineColor(colors[i]);
      }
      histograms_collection[i]->SetMarkerStyle(markers[i]);

      leg->AddEntry(histograms_collection[i], legendList_string[i], "LP");
    }
  }
  if (strstr(options, "fit") != NULL) {
    for (int i = 0; i < collectionSize; i++) {
      optionalFitCollection[i]->SetNpx(2000);
      optionalFitCollection[i]->Draw("same");
      optionalFitCollection[i]->SetLineColor(histograms_collection[i]->GetLineColor());
    }
  }

  if (collectionSize >= 2) {
    leg->Draw("same");
  }
  // cout << "test6" << endl;

  if (strstr(options, "ratioLine") != NULL) {
    TLine myline(minX,1,maxX,1);
    myline.SetLineColor(kBlack);
    myline.SetLineWidth(1);
    // myline.SetLineStyle(2);
    myline.DrawLine(minX,1,maxX,1);
		canvas->Modified();
		canvas->Update();
  }

  if (strstr(options, "150MevLine") != NULL) {
    float lineEdgesX[4] = {0.150, 0.150};
    float lineEdgesY[4] = {0, yUpMarginScaling*maxY};
    TPolyLine* Line150Mev = new TPolyLine(2, lineEdgesX, lineEdgesY);
    cout << "Line150Mev->GetN() = " << Line150Mev->GetN() << endl;
    if (Line150Mev->GetN() > 0) {
      Line150Mev->Draw("");
      Line150Mev->SetLineColor(kBlack);
    }
  }


  // adds some text on the plot
  TLatex* textInfo = new TLatex();
  textInfo->SetTextSize(0.04);
  textInfo->SetNDC(kTRUE); //remove if I want x,y in TLatex to be in the coordinate system of the histogram
  textInfo->DrawLatex(0.18,0.82,texCollisionDataInfo->Data());
  textInfo->DrawLatex(0.18,0.75,Context);
  // cout << "test7" << endl;


  struct stat st1{};
  if (stat("pdfFolder/", &st1) == -1) {
      mkdir("pdfFolder/", 0700);
  }
  canvas->SaveAs("pdfFolder/"+*pdfName+".pdf");

  struct stat st2{};
  if (stat("pngFolder/", &st2) == -1) {
      mkdir("pngFolder/", 0700);
  }
  canvas->SaveAs("pngFolder/"+*pdfName+".png");

  // for(int iCentralityBin = 0; iCentralityBin < nCentralityBins; iCentralityBin++){
  //   cout << "histograms_collection[0]->GetBinContent(iCentralityBin) = " << histograms_collection[0]->GetBinContent(iCentralityBin) << endl;
  // }
}

void Draw_TH1_Histograms_in_one(TH1D** histograms_collection, const TString* legendList_string, int collectionSize, TString Context, TString* pdfName, TString* &texXtitle, TString* &texYtitle, TString* texCollisionDataInfo, std::array<std::array<float, 2>, 2> drawnWindow, const char options[]) {
  // is here to make optionalFitCollection an actual optional parameter; Draw_TH1_Histograms_in_one can be called without, and in that case optionalFitCollection is created empty for use by the actual Draw_TH1_Histograms_in_one function; it will only be used if 'options' has fit in it
  TF1* optionalFitCollectionDummy[collectionSize];
  Draw_TH1_Histograms_in_one(histograms_collection, legendList_string, collectionSize, Context, pdfName, texXtitle, texYtitle, texCollisionDataInfo, drawnWindow, options, optionalFitCollectionDummy);
}

void Draw_TH1_Histogram(TH1D* H1D_Sigma_asFunctionOf_Centrality, TString Context, TString* pdfName, TString* &texXtitle, TString* &texYtitle, TString* texCollisionDataInfo, std::array<std::array<float, 2>, 2> drawnWindow, const char options[]) {
  TH1D* singleHistArray[1] = {H1D_Sigma_asFunctionOf_Centrality};
  TString dummyLegend[1] = {(TString)""};
  int dummyCollectionSize = 1;
  Draw_TH1_Histograms_in_one(singleHistArray, dummyLegend, dummyCollectionSize, Context, pdfName, texXtitle, texYtitle, texCollisionDataInfo, drawnWindow, options);


  // for(int iCentralityBin = 0; iCentralityBin < nCentralityBins; iCentralityBin++){
  //   cout << "H1D_Sigma_asFunctionOf_Centrality->GetBinContent(iCentralityBin) = " << H1D_Sigma_asFunctionOf_Centrality->GetBinContent(iCentralityBin) << endl;
  // }
}

void Draw_TH2_Histograms(TH2D** histograms_collection, const TString* legendList_string, int collectionSize, TString Context, TString* pdfName, TString* &texXtitle, TString* &texYtitle, TString* texCollisionDataInfo, std::array<std::array<float, 2>, 2> drawnWindow, const char options[], TPolyLine* optionalLine) {

  double width = collectionSize*900;
  double height = 800;
  auto canvas = new TCanvas("canvas"+*pdfName, "canvas"+*pdfName, width, height);
  canvas->SetWindowSize(width + (width - canvas->GetWw()), height + (height - canvas->GetWh()));

  canvas->Divide(collectionSize,1);

  // draws histograms from collection
  for (int i = 0; i < collectionSize; i++) {
    canvas->cd(i+1);
    histograms_collection[i]->Draw("colz");
    histograms_collection[i]->SetXTitle(texXtitle->Data());
    histograms_collection[i]->SetYTitle(texYtitle->Data());
    canvas->cd(i+1)->SetRightMargin(0.18); // if the z-axis ever gets hidden, one can play with this
    if (strstr(options, "logz") != NULL) {
      gPad->SetLogz(); // sets log scale for the current pad
    }
    if (strstr(options, "logy") != NULL) {
      gPad->SetLogy();
    }
    if (strstr(options, "logx") != NULL) {
      gPad->SetLogx();
    }
    // leg->AddEntry(histograms_collection[i], legendList_string[i], "LP");
  }

  if (std::equal(std::begin(drawnWindow), std::end(drawnWindow), std::begin(drawnWindowAuto), std::end(drawnWindowAuto))) {
    if (strstr(options, "autoRangeSame") != NULL) {
      int maxXbin = 0;
      int maxYbin = 0;
      int symBinLimitMin = 9999999;
      for (int i = 0; i < collectionSize; i++) {
        if (maxXbin < histograms_collection[i]->FindLastBinAbove(1, 1)) {
          maxXbin = histograms_collection[i]->FindLastBinAbove(1, 1);// (asks for the first/last bin on the x axis (axis number 1) to have strictly more than 1 entry)
        }
        if (maxYbin < histograms_collection[i]->FindLastBinAbove(1, 2)) {
          maxYbin = histograms_collection[i]->FindLastBinAbove(1, 2);// (asks for the first/last bin on the y axis (axis number 2) to have strictly more than 1 entry)
        }
        if (strstr(options, "autoRangeSameSym") != NULL) {
          int symBinLimit = min(histograms_collection[i]->FindFirstBinAbove(1, 2), abs(histograms_collection[i]->GetNbinsY() - histograms_collection[i]->FindLastBinAbove(1, 2)));
          if (symBinLimitMin > symBinLimit) {
            symBinLimitMin = symBinLimit;
          }
        }
      }
      for (int i = 0; i < collectionSize; i++) {
        if (strstr(options, "autoRangeSameSym") != NULL) {
          int symBinLimit = min(histograms_collection[i]->FindFirstBinAbove(1, 2), abs(histograms_collection[i]->GetNbinsY() - histograms_collection[i]->FindLastBinAbove(1, 2))); //(asks for the first/last bin on the y axis (axis number 2) to have strictly more than 1 entry)
          histograms_collection[i]->GetYaxis()->SetRange(symBinLimitMin, histograms_collection[i]->GetNbinsY() - symBinLimitMin); //getting symmetric window around 0 on Y axis
        }
        else {
          histograms_collection[i]->GetXaxis()->SetRange(1, maxXbin);
          histograms_collection[i]->GetYaxis()->SetRange(1, maxYbin);
        }
      }
    }
    // else draws with unchanged histogram window
  } else { // if not drawnWindowAuto, then set the window to match the one entered in drawnWindow
    int minX = drawnWindow[0][0];
    int maxX = drawnWindow[0][1];
    int minY = drawnWindow[1][0];
    int maxY = drawnWindow[1][1];
    for (int i = 0; i < collectionSize; i++) {
      histograms_collection[i]->GetXaxis()->SetRangeUser(minX, maxX);
      histograms_collection[i]->GetYaxis()->SetRangeUser(minY, maxY);
    }
  }

  if (strstr(options, "drawLines") != NULL) {
    cout << "optionalLine->GetN() = " << optionalLine->GetN() << endl;
    if (optionalLine->GetN() > 0) {
      optionalLine->Draw("");
      optionalLine->SetLineColor(kRed);
    }
  }



  gStyle->SetPalette(kBird); // a better palette than the kRainbow that was used by default; https://root.cern.ch/doc/master/classTColor.html lists it as one of the better palettes for Colour Vision Deficiencies 

  // // adds some text on the plot
  TLatex* textInfo = new TLatex();
  textInfo->SetTextSize(0.04);
  textInfo->SetNDC(kTRUE); //remove if I want x,y in TLatex to be in the coordinate system of the histogram
  for (int i = 0; i < collectionSize; i++) {
    canvas->cd(i+1);
    textInfo->DrawLatex(0.18,0.82,texCollisionDataInfo->Data());
    textInfo->DrawLatex(0.18,0.75,Context);
    textInfo->DrawLatex(0.18,0.65,legendList_string[i]);
  }

  struct stat st1{};
  if (stat("pdfFolder/", &st1) == -1) {
      mkdir("pdfFolder/", 0700);
  }
  canvas->SaveAs("pdfFolder/"+*pdfName+".pdf");

  struct stat st2{};
  if (stat("pngFolder/", &st2) == -1) {
      mkdir("pngFolder/", 0700);
  }
  canvas->SaveAs("pngFolder/"+*pdfName+".png");
}

void Draw_TH2_Histograms(TH2D** histograms_collection, const TString* legendList_string, int collectionSize, TString Context, TString* pdfName, TString* &texXtitle, TString* &texYtitle, TString* texCollisionDataInfo, std::array<std::array<float, 2>, 2> drawnWindow, const char options[]) {
  // is here to make optionalFitCollection an actual optional parameter; Draw_TH1_Histograms_in_one can be called without, and in that case optionalFitCollection is created empty for use by the actual Draw_TH1_Histograms_in_one function; it will only be used if 'options' has fit in it
  TPolyLine* optionalLine;
  Draw_TH2_Histograms(histograms_collection, legendList_string, collectionSize, Context, pdfName, texXtitle, texYtitle, texCollisionDataInfo, drawnWindow, options, optionalLine);
}