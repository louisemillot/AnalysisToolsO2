#include "HistogramUtilities.h"
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////// Histogram Utilities ///////////////////////////////////////////////////////////////////////////////////
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

void Draw_TH1_Histograms_in_one(TH1D** histograms_collection, const TString* legendList_string, int collectionSize, TString Context, TString* pdfName, TString* &texXtitle, TString* &texYtitle, TString* texCollisionDataInfo, const char options[]) {
  // has options:
  // - "standardratio" : if in the options string, then an additional plot is drawn, the ratio of the histograms in the collection to the first histogram of the collection; the Y range is chosen automatically based on the difference to 1
  // - "autoratio" : if in the options string, then an additional plot is drawn, the ratio of the histograms in the collection to the first histogram of the collection; the Y range is [0,2.2]
  // - "logy" : if in the options string, then the y axis of the plot is set to a log scale (except for the ratio plot)
  // - "avoidFirst" : if in the options string, then the first histogram of the collection isn't plotted

  // canvas settings
  TCanvas *canvas = new TCanvas ("canvas"+*pdfName, "canvas"+*pdfName, 800, 800);
  canvas->cd(0);

  float minY_collection[collectionSize];
  float maxY_collection[collectionSize];
  float minX_collection[collectionSize];
  float maxX_collection[collectionSize];

  // char* ratioOptionCheck = "ratio";
  // char* OptionPtr = &options;

  for (int i = 0; i < collectionSize; i++) {
    maxX_collection[i] = histograms_collection[i]->GetXaxis()->GetXmax();
    minX_collection[i] = histograms_collection[i]->GetXaxis()->GetXmin();
    maxY_collection[i] = histograms_collection[i]->GetMaximum();
    if (strstr(options, "logy") != NULL) { //  || strstr(options, "ratio") != NULL not sure why I had this here
      minY_collection[i] = histograms_collection[i]->GetMinimum(1E-10);
    }
    else {
      minY_collection[i] = histograms_collection[i]->GetMinimum();
    }
  }

  float maxX = findMaxFloat(maxX_collection, collectionSize);
  float minX = findMinFloat(minX_collection, collectionSize);
  float maxY = findMaxFloat(maxY_collection, collectionSize);
  float minY = findMinFloat(minY_collection, collectionSize);
  float yUpMarginScaling = 1.4;


  if (strstr(options, "customTest") != NULL) {
    minY = 0.6;
    yUpMarginScaling = 1.1;
  }

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
  }

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

  TH1 *hFrame = canvas->DrawFrame(minX,minY,maxX,yUpMarginScaling*maxY);
  if (strstr(options, "logy") != NULL) {
    canvas->SetLogy();
  }  
  if (strstr(options, "logx") != NULL) {
    canvas->SetLogx();
  }
  
  hFrame->SetXTitle(texXtitle->Data());
  hFrame->SetYTitle(texYtitle->Data());
  // hFrame->GetYaxis()->SetTitleOffset(2);

  // legend settings
  TLegend * leg = new TLegend(0.7, 0.75, 0.8, 0.87);
  leg->SetTextSize(gStyle->GetTextSize()*0.7);
  if (collectionSize > 5) { // maybe fine tune that
  leg->SetTextSize(gStyle->GetTextSize()*0.3);
  }

  if (collectionSize >= 6) {
    gStyle->SetPalette(kRainbow); // for the choice of marker's colours; only use this if we have many histograms in the same plot
  }

  // draws histograms from collection
  for (int i = 0; i < collectionSize; i++) {
    if (i!=0 || strstr(options, "avoidFirst") == NULL) { // if i=0 requires that the option avoidFirst isn't there
      if (collectionSize >= 6) {
        histograms_collection[i]->Draw("same PMC PLC"); // PMC uses the palette chosen with gStyle->SetPalette() to chose the colours of the markers, PLC for the lines
      }
      else {
        if (strstr(options, "text") != NULL) {
          histograms_collection[i]->Draw("same TEXT");
        }
        else {
          histograms_collection[i]->Draw("same");
        }
        histograms_collection[i]->SetMarkerColor(colors[i]);
        histograms_collection[i]->SetLineColor(colors[i]);
      }
      histograms_collection[i]->SetMarkerStyle(markers[i]);

      leg->AddEntry(histograms_collection[i], legendList_string[i], "LP");
    }
  }

  if (collectionSize >= 2) {
    leg->Draw("same");
  }

  // adds some text on the plot
  TLatex* textInfo = new TLatex();
  textInfo->SetTextSize(0.04);
  textInfo->SetNDC(kTRUE); //remove if I want x,y in TLatex to be in the coordinate system of the histogram
  textInfo->DrawLatex(0.18,0.82,texCollisionDataInfo->Data());
  textInfo->DrawLatex(0.18,0.75,Context);

  canvas->SaveAs(*pdfName+".pdf");

  struct stat st{};
  if (stat("pngFolder/", &st) == -1) {
      mkdir("pngFolder/", 0700);
  }
  canvas->SaveAs("pngFolder/"+*pdfName+".png");
}



void Draw_TH2_Histograms(TH2D** histograms_collection, const TString* legendList_string, int collectionSize, TString Context, TString* pdfName, TString* &texXtitle, TString* &texYtitle, TString* texCollisionDataInfo, const char options[]) {

  double width = 2300;
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
  canvas->SaveAs(*pdfName+".pdf");

  struct stat st{};
  if (stat("pngFolder/", &st) == -1) { // checks if pngFolder exists in the current directory 
      mkdir("pngFolder/", 0777); // the argument dictates the permissions; here should give full permissions 
  }
  canvas->SaveAs("pngFolder/"+*pdfName+".png");
}