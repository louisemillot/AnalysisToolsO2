
#ifndef HISTOGRAM_PLOTTING_C
#define HISTOGRAM_PLOTTING_C


#include "HistogramPlotting.h"
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
////////////////////////////////////////////////////////////////////////////// Histogram Context /////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

TString contextCustomThreeFields(TString mainContext, TString secondaryContext, TString tertiaryContext, __attribute__ ((unused)) std::string options){
  TString texContextFinal;
  texContextFinal = "#splitline{"+mainContext+" "+secondaryContext+"}{"+tertiaryContext+"}";
  // texContextFinal = "#splitline{"+mainContext+" "+secondaryContext+"}{#splitline{2023 QC}{test"+tertiaryContext+"}}";
  // texContextFinal = "testtesttestest";
  return texContextFinal;
}

TString contextCustomTwoFields(TString mainContext, TString secondaryContext, std::string options){
  return contextCustomThreeFields(mainContext, (TString)"" , secondaryContext, options);
}

TString contextCustomOneField(TString mainContext, std::string options){
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

TString contextCentRange(float* CentRange){
  std::stringstream ss;
    ss << "" << CentRange[0] << "-" << CentRange[1] << "%";
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






TString contextDatasetRadiusCompAndVarRange(TString mainContext, int iDataset, float* variableRange, std::string options){
  TString texcontextDatasetRadiusCompAndVarRange;
  if (options.find("pt") != std::string::npos) { //  || options.find("ratio") != NULL not sure why I had this here
    texcontextDatasetRadiusCompAndVarRange = "#splitline{"+mainContext+" "+DatasetsNames[iDataset]+"}{#splitline{2023 QC}{"+contextPtRange(variableRange)+"}}";
  }
  if (options.find("eta") != std::string::npos) { //  || options.find("ratio") != NULL not sure why I had this here
    texcontextDatasetRadiusCompAndVarRange = "#splitline{"+mainContext+" "+DatasetsNames[iDataset]+"}{#splitline{2023 QC}{"+contextEtaRange(variableRange)+"}}";
  }

  return texcontextDatasetRadiusCompAndVarRange;
}

TString contextDatasetCompAndRadiusAndVarRange(TString mainContext, float jetRadius, float* variableRange, std::string options){
  TString texcontextDatasetCompAndRadiusAndVarRange;
  if (options.find("pt") != std::string::npos) { //  || options.find("ratio") != NULL not sure why I had this here
    texcontextDatasetCompAndRadiusAndVarRange = "#splitline{"+mainContext+"}{#splitline{"+contextJetRadius(jetRadius)+"}{"+contextPtRange(variableRange)+"}}";
  }
  if (options.find("eta") != std::string::npos) { //  || options.find("ratio") != NULL not sure why I had this here
    texcontextDatasetCompAndRadiusAndVarRange = "#splitline{"+mainContext+"}{#splitline{"+contextJetRadius(jetRadius)+"}{"+contextEtaRange(variableRange)+"}}";
  }
  if (options.find("centrality") != std::string::npos) { //  || options.find("ratio") != NULL not sure why I had this here
    texcontextDatasetCompAndRadiusAndVarRange = "#splitline{"+mainContext+"}{#splitline{"+contextJetRadius(jetRadius)+"}{"+contextCentRange(variableRange)+"}}";
  }

  return texcontextDatasetCompAndRadiusAndVarRange;
}

TString contextDatasetCompAndRadius(TString mainContext, float jetRadius, __attribute__ ((unused)) std::string options){
  TString texcontextDatasetCompAndRadius;
  texcontextDatasetCompAndRadius = "#splitline{"+mainContext+"}{"+contextJetRadius(jetRadius)+"}";

  return texcontextDatasetCompAndRadius;
}

TString contextDatasetComp(TString mainContext, __attribute__ ((unused)) std::string options){
  TString texcontextDatasetComp;
  texcontextDatasetComp = mainContext;

  return texcontextDatasetComp;
}

void CentralityLegend(TString* centralityLegend, const float arrayCentralityIntervals[][2], int nCentralityBins){
  std::stringstream ss;
  ss.precision(2);
  for(int iCentralityBin = 0; iCentralityBin < nCentralityBins; iCentralityBin++){
    ss << "" << arrayCentralityIntervals[iCentralityBin][0] << "-" << arrayCentralityIntervals[iCentralityBin][1] << "%";
    centralityLegend[iCentralityBin] = (TString)ss.str();
    ss.str("");
    ss.clear();
  }
}

void IterationLegend(TString* iterationLegend, int unfoldIterationMin, int unfoldIterationMax, int step){
  const int nUnfoldIteration = std::floor((unfoldIterationMax - unfoldIterationMin + 1)/step);
  std::stringstream ss;
  ss.precision(2);
  for(int iUnfoldIteration = 0; iUnfoldIteration < nUnfoldIteration; iUnfoldIteration++){
    ss << "k_{unfold} = " << iUnfoldIteration * step + unfoldIterationMin;
    iterationLegend[iUnfoldIteration] = (TString)ss.str();
    ss.str("");
    ss.clear();
  }
}





//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////// Histogram Drawing /////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


void Draw_TH1_Histograms_in_one(TH1D** histograms_collection, const TString* legendList_string, int collectionSize, TString Context, TString* pdfName, TString* &texXtitle, TString* &texYtitle, TString* texCollisionDataInfo, std::array<std::array<float, 2>, 2> drawnWindow, std::array<std::array<float, 2>, 2> legendPlacement, std::array<float, 2> contextPlacement, std::string options, TF1** optionalFitCollection) {
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
    if (options.find("autoXrange") != std::string::npos) {
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
    
    if (options.find("logy") != std::string::npos) { //  || options.find("ratio") != NULL not sure why I had this here
      minY_collection[i] = histograms_collection[i]->GetMinimum(1E-10);
      // cout << "test1.3a" << endl;
    }
    else if (options.find("minYnotZero") != std::string::npos) {
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

  yDownMarginScaling = 1;
  maxX = findMaxFloat(maxX_collection, collectionSize);
  minX = findMinFloat(minX_collection, collectionSize);
  maxY = findMaxFloat(maxY_collection, collectionSize);
  minY = findMinFloat(minY_collection, collectionSize);
  maxY > 0 ? yUpMarginScaling = 1.4 : yUpMarginScaling = 0.6;

  if (options.find("logy") != std::string::npos) {
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
    if (minY < 0 || options.find("minYnotZero") != std::string::npos) {
      minY > 0 ? yDownMarginScaling = 0.95 : yDownMarginScaling = 1.05;
      yUpMarginScaling = 1.1;
    } else {
      minY = 0.;
    }
  // cout << "test3" << endl;

    if (options.find("autoratio") != std::string::npos) {
      float deltaMax = max(1-minY, maxY-1);
      minY = max((float)0., 1-deltaMax);
      maxY = 1+deltaMax;
      yUpMarginScaling = 1.3;
    }
    if (options.find("standardratio") != std::string::npos) {
      minY = 0;
      maxY = 2;
      yUpMarginScaling = 1.1;
    }
    if (options.find("zoomedratio") != std::string::npos) {
      minY = 0.8;
      maxY = 1.2;
      yUpMarginScaling = 1.1;
    }
    if (options.find("zoomedextraratio") != std::string::npos) {
      minY = 0.9;
      maxY = 1.1;
      yUpMarginScaling = 1.1;
    }
    if (options.find("efficiency") != std::string::npos) {
      minY = 0;
      maxY = 1.5;
      yUpMarginScaling = 1.1;
    }
  }

  if (!std::equal(std::begin(drawnWindow[0]), std::end(drawnWindow[0]), std::begin(drawnWindowAuto[0]), std::end(drawnWindowAuto[0]))) {
    minX = drawnWindow[0][0];
    maxX = drawnWindow[0][1];
  }
  if (!std::equal(std::begin(drawnWindow[1]), std::end(drawnWindow[1]), std::begin(drawnWindowAuto[1]), std::end(drawnWindowAuto[1]))) {
    yUpMarginScaling = 1;
    yDownMarginScaling = 1;
    minY = drawnWindow[1][0];
    maxY = drawnWindow[1][1];
  }




  TH1 *hFrame = canvas->DrawFrame(minX, yDownMarginScaling*minY, maxX, yUpMarginScaling*maxY);
  if (options.find("logy") != std::string::npos) {
    canvas->SetLogy();
  }  
  if (options.find("logx") != std::string::npos) {
    canvas->SetLogx();
  }
    // cout << "test4" << endl;

  hFrame->SetXTitle(texXtitle->Data());
  hFrame->SetYTitle(texYtitle->Data());
  // hFrame->GetYaxis()->SetTitleOffset(2);

  // legend settings
  double xLeftLegend = 0.7;
  double xRightLegend = 0.75;
  double yLowLegend = 0.8;
  double yUpLegend = 0.87;
  if (!std::equal(std::begin(legendPlacement[0]), std::end(legendPlacement[0]), std::begin(legendPlacementAuto[0]), std::end(legendPlacementAuto[0]))) {
    xLeftLegend = legendPlacement[0][0];
    xRightLegend = legendPlacement[0][1];
  }
  if (!std::equal(std::begin(legendPlacement[1]), std::end(legendPlacement[1]), std::begin(legendPlacementAuto[1]), std::end(legendPlacementAuto[1]))) {
    // leg = TLegend(0.7, 0.75, legendPlacement[1][0], legendPlacement[1][1]);
    // leg->SetY1NDC(legendPlacement[1][0]);
    // leg->SetY2NDC(legendPlacement[1][1]);
    yLowLegend = legendPlacement[1][0];
    yUpLegend = legendPlacement[1][1];
  }
  TLegend * leg = new TLegend(xLeftLegend, xRightLegend, yLowLegend, yUpLegend);

  leg->SetTextSize(gStyle->GetTextSize()*0.7);
  if (collectionSize >= 6) { // maybe fine tune that
    leg->SetTextSize(gStyle->GetTextSize()*0.3);
  }
  if (options.find("fit") != std::string::npos) {
    leg->SetTextSize(gStyle->GetTextSize()*0.3);
  }
  if (collectionSize >= 6) {
    gStyle->SetPalette(kRainbow); // for the choice of marker's colours; only use this if we have many histograms in the same plot
  }
  // cout << "test5" << endl;

  // draws histograms from collection, and setting the colors
  for (int i = 0; i < collectionSize; i++) {
    if (i!=0 || options.find("avoidFirst") == std::string::npos) { // if i=0 requires that the option avoidFirst isn't there (== std::string::npos means it didn't find it in the elements 0 to npos-1, where npos is the size of the string options)

      if (options.find("colorPairs") == std::string::npos) {
        if (collectionSize >= 6) {
          histograms_collection[i]->Draw("same PMC PLC"); // PMC uses the palette chosen with gStyle->SetPalette() to chose the colours of the markers, PLC for the lines
          if (options.find("histWithLine") != std::string::npos) {
            histograms_collection[i]->Draw("][ Hist same PMC PLC"); // PMC uses the palette chosen with gStyle->SetPalette() to chose the colours of the markers, PLC for the lines
          }
        } else {
          histograms_collection[i]->Draw("same");
          if (options.find("histWithLine") != std::string::npos) {
            histograms_collection[i]->Draw("][ Hist same");
          }
          histograms_collection[i]->SetMarkerColor(colors[i]);
          histograms_collection[i]->SetLineColor(colors[i]);
        }
        histograms_collection[i]->SetMarkerStyle(markers[i]);

      } else {
        if (collectionSize >= 2*6) {
          int nColors = gStyle->GetNumberOfColors();
          int histoColor = (float)nColors / collectionSize * (int)i/2;
          histograms_collection[i]->SetLineColor(gStyle->GetColorPalette(histoColor));
          histograms_collection[i]->SetMarkerColor(gStyle->GetColorPalette(histoColor));
          histograms_collection[i]->Draw("same"); // PMC uses the palette chosen with gStyle->SetPalette() to chose the colours of the markers, PLC for the lines
          if (options.find("histWithLine") != std::string::npos) {
            histograms_collection[i]->Draw("][ Hist same"); // PMC uses the palette chosen with gStyle->SetPalette() to chose the colours of the markers, PLC for the lines
          }
        } else {
          histograms_collection[i]->Draw("same");
          if (options.find("histWithLine") != std::string::npos) {
            histograms_collection[i]->Draw("][ Hist same");
          }
          histograms_collection[i]->SetMarkerColor(colors[(int)i/2]);
          histograms_collection[i]->SetLineColor(colors[(int)i/2]);
        }
        histograms_collection[i]->SetMarkerStyle(markersColorPairs[i]);
      }
      leg->AddEntry(histograms_collection[i], legendList_string[i], "LP");
    }
  }
  if (options.find("fit") != std::string::npos) {
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

  if (options.find("ratioLine") != std::string::npos) {
    TLine myline(minX,1,maxX,1);
    myline.SetLineColor(kBlack);
    myline.SetLineWidth(1);
    // myline.SetLineStyle(2);
    myline.DrawLine(minX,1,maxX,1);
		canvas->Modified();
		canvas->Update();
    cout << "minX = " << minX << endl;
  }

  if (options.find("150MevLine") != std::string::npos) {
    float lineEdgesX[4] = {0.150, 0.150};
    float lineEdgesY[4] = {minY, yUpMarginScaling*maxY};
    TPolyLine* Line150Mev = new TPolyLine(2, lineEdgesX, lineEdgesY);
    // cout << "Line150Mev->GetN() = " << Line150Mev->GetN() << endl;
    if (Line150Mev->GetN() > 0) {
      Line150Mev->Draw("");
      Line150Mev->SetLineColor(kBlack);
    }
  }

  if (options.find("datasetXaxisBinLabels") != std::string::npos) {
    // ChangeLabel(labNum, labAngle, labSize, labAlign, labColor, labFont, labText)
    // [in]	labNum	Number of the label to be changed, negative numbers start from the end
    // [in]	labAngle	New angle value
    // [in]	labSize	New size (0 erase the label)
    // [in]	labAlign	New alignment value
    // [in]	labColor	New label color
    // [in]	labFont	New label font
    // [in]	labText	New label text
    int nBins = histograms_collection[0]->GetNbinsX();
    for (int i = 0; i < nBins; i++) {
      hFrame->GetXaxis()->SetBinLabel(hFrame->GetNbinsX()/(2*nBins) * (2*i + 1), DatasetsNames[i]); // DrawFrame creates a histo with 1000 bins; takes bin number as input
      hFrame->GetXaxis()->ChangeLabel(1, 30, -1, -1, -1, -1, -1); // didn't manage to get it to work, but also didnt spend much time; maybe because it takes label number as input?
    }
  }

  // Context drawing 
  TLatex* textInfo = new TLatex();
  textInfo->SetTextSize(0.04);
  textInfo->SetNDC(kTRUE); //remove if I want x,y in TLatex to be in the coordinate system of the histogram  

  double xTopLeftCornerContext = 0.18;
  double yTopLeftCornerContext = 0.82;
  double deltaYContextsPosition = 0.07;
  if (contextPlacement[0] != contextPlacementAuto[0]) {
    xTopLeftCornerContext = contextPlacement[0];
  }
  if (contextPlacement[1] != contextPlacementAuto[1]) {
    yTopLeftCornerContext = contextPlacement[1];
  }
  textInfo->DrawLatex(xTopLeftCornerContext,yTopLeftCornerContext,texCollisionDataInfo->Data());
  textInfo->DrawLatex(xTopLeftCornerContext,yTopLeftCornerContext - deltaYContextsPosition,Context);




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

void Draw_TH1_Histograms_in_one(TH1D** histograms_collection, const TString* legendList_string, int collectionSize, TString Context, TString* pdfName, TString* &texXtitle, TString* &texYtitle, TString* texCollisionDataInfo, std::array<std::array<float, 2>, 2> drawnWindow, std::array<std::array<float, 2>, 2> legendPlacement, std::array<float, 2> contextPlacement, std::string options) {
  // is here to make optionalFitCollection an actual optional parameter; Draw_TH1_Histograms_in_one can be called without, and in that case optionalFitCollection is created empty for use by the actual Draw_TH1_Histograms_in_one function; it will only be used if 'options' has fit in it
  TF1* optionalFitCollectionDummy[collectionSize];
  Draw_TH1_Histograms_in_one(histograms_collection, legendList_string, collectionSize, Context, pdfName, texXtitle, texYtitle, texCollisionDataInfo, drawnWindow, legendPlacement, contextPlacement, options, optionalFitCollectionDummy);
}

void Draw_TH1_Histogram(TH1D* histogram, TString Context, TString* pdfName, TString* &texXtitle, TString* &texYtitle, TString* texCollisionDataInfo, std::array<std::array<float, 2>, 2> drawnWindow, std::array<std::array<float, 2>, 2> legendPlacement, std::array<float, 2> contextPlacement, std::string options) {
  TH1D* singleHistArray[1] = {histogram};
  TString dummyLegend[1] = {(TString)""};
  int dummyCollectionSize = 1;
  Draw_TH1_Histograms_in_one(singleHistArray, dummyLegend, dummyCollectionSize, Context, pdfName, texXtitle, texYtitle, texCollisionDataInfo, drawnWindow, legendPlacement, contextPlacement, options);
}

void Draw_TH2_Histograms(TH2D** histograms_collection, const TString* legendList_string, int collectionSize, TString Context, TString* pdfName, TString* &texXtitle, TString* &texYtitle, TString* texCollisionDataInfo, std::array<std::array<float, 2>, 2> drawnWindow, std::string options, TPolyLine* optionalLine) {

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
    if (options.find("logz") != std::string::npos) {
      gPad->SetLogz(); // sets log scale for the current pad
    }
    if (options.find("logy") != std::string::npos) {
      gPad->SetLogy();
    }
    if (options.find("logx") != std::string::npos) {
      gPad->SetLogx();
    }
    // leg->AddEntry(histograms_collection[i], legendList_string[i], "LP");
  }

  if (std::equal(std::begin(drawnWindow), std::end(drawnWindow), std::begin(drawnWindowAuto), std::end(drawnWindowAuto))) {
    if (options.find("autoRangeSame") != std::string::npos) {
      int maxXbin = 0;
      int maxYbin = 0;
      int symBinLimitMin = 9999999;
      for (int i = 0; i < collectionSize; i++) {
        if (maxXbin < histograms_collection[i]->FindLastBinAbove(GLOBAL_epsilon, 1)) {
          maxXbin = histograms_collection[i]->FindLastBinAbove(GLOBAL_epsilon, 1);// (asks for the first/last bin on the x axis (axis number 1) to have strictly more than 1 entry)
        }
        if (maxYbin < histograms_collection[i]->FindLastBinAbove(GLOBAL_epsilon, 2)) {
          maxYbin = histograms_collection[i]->FindLastBinAbove(GLOBAL_epsilon, 2);// (asks for the first/last bin on the y axis (axis number 2) to have strictly more than 1 entry)
        }
        if (options.find("autoRangeSameSym") != std::string::npos) {
          int symBinLimit = min(histograms_collection[i]->FindFirstBinAbove(GLOBAL_epsilon, 2), abs(histograms_collection[i]->GetNbinsY() - histograms_collection[i]->FindLastBinAbove(GLOBAL_epsilon, 2)));
          if (symBinLimitMin > symBinLimit) {
            symBinLimitMin = symBinLimit;
          }
        }
      }
      for (int i = 0; i < collectionSize; i++) {
        if (options.find("autoRangeSameSym") != std::string::npos) {
          int symBinLimit = min(histograms_collection[i]->FindFirstBinAbove(GLOBAL_epsilon, 2), abs(histograms_collection[i]->GetNbinsY() - histograms_collection[i]->FindLastBinAbove(GLOBAL_epsilon, 2))); //(asks for the first/last bin on the y axis (axis number 2) to have strictly more than 1 entry)
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

  if (options.find("drawLines") != std::string::npos) {
    // cout << "optionalLine->GetN() = " << optionalLine->GetN() << endl;
    if (optionalLine->GetN() > 0) {
      optionalLine->Draw("");
      optionalLine->SetLineColor(kRed);
    }
  }

  // cout << "---------------- TH2 drawing; getting min and max test:" << endl;
  for (int i = 0; i < collectionSize; i++) {
    histograms_collection[i]->GetZaxis()->SetRangeUser(histograms_collection[i]->GetMinimum(GLOBAL_epsilon), histograms_collection[i]->GetMaximum());
    // cout << "min = " << histograms_collection[i]->GetMinimum(GLOBAL_epsilon) << ", max = " << histograms_collection[i]->GetMaximum() << endl;
  }

  gStyle->SetPalette(kBird); // a better palette than the kRainbow that was used by default; https://root.cern.ch/doc/master/classTColor.html lists it as one of the better palettes for Colour Vision Deficiencies 
  gStyle->SetNumberContours(100);
  
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

void Draw_TH2_Histograms(TH2D** histograms_collection, const TString* legendList_string, int collectionSize, TString Context, TString* pdfName, TString* &texXtitle, TString* &texYtitle, TString* texCollisionDataInfo, std::array<std::array<float, 2>, 2> drawnWindow, std::string options) {
  // is here to make optionalFitCollection an actual optional parameter; Draw_TH1_Histograms_in_one can be called without, and in that case optionalFitCollection is created empty for use by the actual Draw_TH1_Histograms_in_one function; it will only be used if 'options' has fit in it
  TPolyLine* optionalLine;
  Draw_TH2_Histograms(histograms_collection, legendList_string, collectionSize, Context, pdfName, texXtitle, texYtitle, texCollisionDataInfo, drawnWindow, options, optionalLine);
}


void Draw_TH2_Histogram(TH2D* histogram, TString Context, TString* pdfName, TString* &texXtitle, TString* &texYtitle, TString* texCollisionDataInfo, std::array<std::array<float, 2>, 2> drawnWindow, std::string options) {
  TH2D* singleHistArray[1] = {histogram};
  TString dummyLegend[1] = {(TString)""};
  int dummyCollectionSize = 1;
  Draw_TH2_Histograms(singleHistArray, dummyLegend, dummyCollectionSize, Context, pdfName, texXtitle, texYtitle, texCollisionDataInfo, drawnWindow, options);


  // for(int iCentralityBin = 0; iCentralityBin < nCentralityBins; iCentralityBin++){
  //   cout << "histogram->GetBinContent(iCentralityBin) = " << histogram->GetBinContent(iCentralityBin) << endl;
  // }
}









#endif