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

std::vector<double> MakeVariableBinning_twoWidths(double xMin, int nLeft, double xMiddle, double xMax, int nRight) {
  if (xMin > xMiddle || xMiddle > xMax){
    cout << "ERROR in MakeVariableBinning_twoWidths(): it should be xMin < xMiddle < xMax but it is not the case" << endl;
  }
  std::vector<double> bins;
  double binWidth;
  for(int i = 0; i < nLeft; i++){
    binWidth = (xMiddle - xMin) / nLeft;
    bins.push_back(xMin + i * binWidth);
    // cout << "ibin " << iBin << ": lowEdge = " << H1_histo->GetBinLowEdge(iBin) << endl;
  }
  for(int i = 0; i < nRight+1; i++){ // using nRight +1 here so that the rightmost edge is added. But due to double prec, might not be exactly xMax; alternative could be to forego the +1 and manually pushback xMax
    binWidth = (xMax - xMiddle) / nRight;
    bins.push_back(xMiddle + i * binWidth);
    // cout << "ibin " << iBin << ": lowEdge = " << H1_histo->GetBinLowEdge(iBin) << endl;
  }
  return bins;
}

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

TH2D RebinVariableBins2D(TH2D* H2D_hist, int nBinsX, int nBinsY, double* binsX, double* binsY, bool debug){
  if (debug == true) {cout << "___________ begin RebinVariableBins2D, nBinsX = " << nBinsX << ", nBinsY = " << nBinsY << endl;}
  // int debugCount = 0;

  TH2D H2D_hist_rebinned("H2D_hist_rebinned", "H2D_hist_rebinned", nBinsX, binsX, nBinsY, binsY);


  // int ibinY_low_test1 = H2D_hist->GetYaxis()->FindBin(H2D_hist_rebinned.GetYaxis()->GetBinLowEdge(1) );
  // int ibinY_high_test1 = H2D_hist->GetYaxis()->FindBin(H2D_hist_rebinned.GetYaxis()->GetBinLowEdge(1+1) ) - 1;
  // int ibinY_low_test5 = H2D_hist->GetYaxis()->FindBin(H2D_hist_rebinned.GetYaxis()->GetBinLowEdge(5) );
  // int ibinY_high_test5 = H2D_hist->GetYaxis()->FindBin(H2D_hist_rebinned.GetYaxis()->GetBinLowEdge(5+1) ) - 1;
  // cout << "pre rebin: H2D_hist Integral(0,-1, 1, 1)" << H2D_hist->Integral(0,-1, ibinY_low_test1, ibinY_high_test1) << endl;
  // cout << "pre rebin: H2D_hist Integral(0,-1, 5, 5)" << H2D_hist->Integral(0,-1, ibinY_high_test5, ibinY_high_test5) << endl;
  // cout << "pre rebin: H2D_hist Integral(0,-1, 0, -1)" << H2D_hist->Integral(0,-1, 0, -1) << endl;

  int ibinX_low, ibinX_high, ibinY_low, ibinY_high;
  double xLowEdgeRebin, xHighEdgeRebin, yLowEdgeRebin, yHighEdgeRebin;
  double H2D_hist_content, H2D_hist_contentError;
  for(int iBinX = 0; iBinX <= nBinsX+1; iBinX++){ // 0 and n+1 take underflow and overflow into account
    if (debug == true) {cout << "iBinX = " << iBinX << endl;}
    xLowEdgeRebin = H2D_hist_rebinned.GetXaxis()->GetBinLowEdge(iBinX);
    xHighEdgeRebin = H2D_hist_rebinned.GetXaxis()->GetBinLowEdge(iBinX+1);
    xLowEdgeRebin < binsX[0] ? ibinX_low = 0 : ibinX_low = H2D_hist->GetXaxis()->FindBin(xLowEdgeRebin); // according to definition of getbinlowedge, this takes into account under/overflows thanks to FindBin() giving the low edge of bin 1 minus the bin width (calculated assuming bins of equal sizes in all the axis); in fact, it works regardless of the value of iBin even for Nbin+2 or more (https://root.cern.ch/doc/master/classTAxis.html#a6d2d17ef00382842f58e88f356516a0d)
    xHighEdgeRebin > binsX[nBinsX] ? ibinX_high = H2D_hist->GetNbinsX()+1 : ibinX_high = H2D_hist->GetXaxis()->FindBin(xHighEdgeRebin) - 1;
    for(int iBinY = 0; iBinY <= nBinsY+1; iBinY++){ // 0 and n+1 take underflow and overflow into account
      yLowEdgeRebin = H2D_hist_rebinned.GetYaxis()->GetBinLowEdge(iBinY);
      yHighEdgeRebin = H2D_hist_rebinned.GetYaxis()->GetBinLowEdge(iBinY+1);
      yLowEdgeRebin < binsY[0] ? ibinY_low = 0 : ibinY_low = H2D_hist->GetYaxis()->FindBin(yLowEdgeRebin); // according to definition of getbinlowedge, it does work for under/overflows by giving the low edge of bin 1 minus the bin width (calculated assuming bins of equal sizes in all the axis); in fact, it works regardless of the value of iBin even for Nbin+2 or more (https://root.cern.ch/doc/master/classTAxis.html#a6d2d17ef00382842f58e88f356516a0d)
      yHighEdgeRebin > binsY[nBinsY] ? ibinY_high =  H2D_hist->GetNbinsY()+1 : ibinY_high = H2D_hist->GetYaxis()->FindBin(yHighEdgeRebin) - 1;

      H2D_hist_content = H2D_hist->IntegralAndError(ibinX_low, ibinX_high, ibinY_low, ibinY_high, H2D_hist_contentError);
      H2D_hist_rebinned.SetBinContent(iBinX, iBinY, H2D_hist_content);
      H2D_hist_rebinned.SetBinError(iBinX, iBinY, H2D_hist_contentError);
      if (debug == true) {cout << "ibinX_low = " << ibinX_low << ", ibinX_high = " << ibinX_high << ", ibinY_low = " << ibinY_low << ", ibinY_high = " << ibinY_high << "         --------          H2D_hist_rebinned(" << iBinX <<", "<< iBinY <<") = " << H2D_hist_rebinned.GetBinContent(iBinX, iBinY) << endl;}
    }
  }

  // cout << "post rebin: H2D_hist Integral(0,-1, 1, 1)" << H2D_hist_rebinned.Integral(0,-1, 1, 1) << endl;
  // cout << "post rebin: H2D_hist Integral(0,-1, 5, 5)" << H2D_hist_rebinned.Integral(0,-1, 5, 5) << endl;
  // cout << "post rebin: H2D_hist Integral(0,-1, 0, -1)" << H2D_hist_rebinned.Integral(0,-1, 0, -1) << endl;

  // if poisson      : bins A+B -> e= sqrt(N_A +N_B)
  // if propag errors: bins A+B -> e= sqrt(e(A)2 + e(B)2) = sqrt(N_A + N_B)  is the same

  H2D_hist_rebinned.ResetStats(); // setbincontent interacts weirdly with getentries(); resetstats makes it so that getentries gives the sum of bin contents correctly
  if (debug == true) {cout << "RebinVariableBins2D - what of the errors" << endl;}
  if (debug == true) {cout << "H2D_hist nbins X =" << H2D_hist->GetNbinsX() << endl;}
  if (debug == true) {cout << "H2D_hist_rebinned nbins X =" << H2D_hist_rebinned.GetNbinsX() << endl;}

  std::stringstream ss;
  ss << H2D_hist->GetName() << "_RebinVariableBins2D";
  TString histName((TString)ss.str());

  if (debug == true) {cout << "___________ end RebinVariableBins2D" << endl;}

  return H2D_hist_rebinned;
}
TH2D RebinVariableBins2D(TH2D* H2D_hist, int nBinsX, int nBinsY, double* binsX, double* binsY){
  return RebinVariableBins2D(H2D_hist, nBinsX, nBinsY, binsX, binsY, false);
}

TH2D RebinVariableBins2D_PriorWeightedBinMerging(TH2D* H2D_hist, int nBinsX, int nBinsY, double* binsX, double* binsY, TH1D* H1D_priorSpectrum, bool debug){
  // H1D_priorSpectrum is assumed to have same binning as H2D_hist
  if (debug == true) {cout << "___________ begin RebinVariableBins2D, nBinsX = " << nBinsX << ", nBinsY = " << nBinsY << endl;}

  TH2D H2D_hist_rebinned("H2D_hist_rebinned", "H2D_hist_rebinned", nBinsX, binsX, nBinsY, binsY);
  H2D_hist_rebinned.Sumw2();

  TH2D* H2D_hist_temp = (TH2D*)H2D_hist->Clone(H2D_hist->GetName()+(TString)"_temp"); // needs this as its content will be changed for the weighting, we do not want to modify H2D_hist

  // declaration of variables used when calculating content and error of histograms
  double H2D_hist_content;
  double H2D_hist_contentError = 0;
  double H2D_hist_contentErrorA = 0;
  double H2D_hist_contentErrorB = 0;


  // prior weighting
  double priorWeightedSpectrumContent = 0;
  double priorWeightedSpectrumError = 0;
  for(int iBinXFine = 0; iBinXFine <= H2D_hist_temp->GetNbinsX()+1; iBinXFine++){ // 0 and n+1 take underflow and overflow into account
    // if (debug == true) {cout << "iBinXFine = " << iBinXFine << endl;}
    priorWeightedSpectrumContent = H1D_priorSpectrum->GetBinContent(iBinXFine);
    priorWeightedSpectrumError = H1D_priorSpectrum->GetBinError(iBinXFine);
    for(int iBinYFine = 0; iBinYFine <= H2D_hist_temp->GetNbinsY()+1; iBinYFine++){ // 0 and n+1 take underflow and overflow into account
      // if (debug == true) {cout << "     iBinYFine = " << iBinYFine << endl;}
      H2D_hist_content = H2D_hist_temp->GetBinContent(iBinXFine, iBinYFine);
      H2D_hist_contentError = H2D_hist_temp->GetBinError(iBinXFine, iBinYFine);
      H2D_hist_temp->SetBinContent(iBinXFine, iBinYFine, H2D_hist_content * priorWeightedSpectrumContent);
      H2D_hist_temp->SetBinError(iBinXFine, iBinYFine, sqrt( priorWeightedSpectrumContent*priorWeightedSpectrumContent * H2D_hist_contentError*H2D_hist_contentError + H2D_hist_content*H2D_hist_content * priorWeightedSpectrumError*priorWeightedSpectrumError));
      // cout << "iBinXFine = " << iBinXFine << ", iBinYFine = " << iBinYFine << "         --------          priorWeightedSpectrumContent = " << priorWeightedSpectrumContent << ", H2D_hist_contentError = " << H2D_hist_contentError << endl;
      // cout << "                                         --------          priorWeightedSpectrumError = " << priorWeightedSpectrumError << ", H2D_hist_content = " << H2D_hist_content << endl;
      // cout << "                                         --------          H2hist resulting error (B2*SigmaA2 + A2*SigmaB2) = " << H2D_hist->GetBinError(iBinXFine, iBinYFine) << endl;
    }
  }

  cout << "##############################################" << endl;
  cout << "######### reweighting with prior #############" << endl;
  cout << "##############################################" << endl;

  int ibinX_low, ibinX_high, ibinY_low, ibinY_high;
  double xLowEdgeRebin, xHighEdgeRebin, yLowEdgeRebin, yHighEdgeRebin;
  double priorBinNorm, priorBinNormError;
  double hist2D_rebinned_content;
  double hist2D_rebinned_error;
  cout << "to think about: do I really give the value 0 if denominator is 0 ? when normalising by 1./priorBinNorm" << endl;
  for(int iBinY = 0; iBinY <= nBinsY+1; iBinY++){ // overflow and underflow are taken care of later in the function (in the truncation part)
    if (debug == true) {cout << "iBinY = " << iBinY << endl;}
    yLowEdgeRebin = H2D_hist_rebinned.GetYaxis()->GetBinLowEdge(iBinY);
    yHighEdgeRebin = H2D_hist_rebinned.GetYaxis()->GetBinLowEdge(iBinY+1);
    yLowEdgeRebin < binsY[0] ? ibinY_low = 0 : ibinY_low = H2D_hist_temp->GetYaxis()->FindBin(yLowEdgeRebin); // according to definition of getbinlowedge, it does work for under/overflows by giving the low edge of bin 1 minus the bin width (calculated assuming bins of equal sizes in all the axis); in fact, it works regardless of the value of iBin even for Nbin+2 or more (https://root.cern.ch/doc/master/classTAxis.html#a6d2d17ef00382842f58e88f356516a0d)
    yHighEdgeRebin > binsY[nBinsY] ? ibinY_high =  H2D_hist_temp->GetNbinsY()+1 : ibinY_high = H2D_hist_temp->GetYaxis()->FindBin(yHighEdgeRebin) - 1;

    for(int iBinX = 0; iBinX <= nBinsX+1; iBinX++){ // overflow and underflow are taken care of later in the function (in the truncation part)
    // for(int iBinX = 0; iBinX <= nBinsX+1; iBinX++){ // 0 and n+1 take underflow and overflow into account
      // if (debug == true) {cout << "     iBinX = " << iBinX << endl;}
      xLowEdgeRebin = H2D_hist_rebinned.GetXaxis()->GetBinLowEdge(iBinX);
      xHighEdgeRebin = H2D_hist_rebinned.GetXaxis()->GetBinLowEdge(iBinX+1);
      xLowEdgeRebin < binsX[0] ? ibinX_low = 0 : ibinX_low = H2D_hist_temp->GetXaxis()->FindBin(xLowEdgeRebin); // according to definition of getbinlowedge, this takes into account under/overflows thanks to FindBin() giving the low edge of bin 1 minus the bin width (calculated assuming bins of equal sizes in all the axis); in fact, it works regardless of the value of iBin even for Nbin+2 or more (https://root.cern.ch/doc/master/classTAxis.html#a6d2d17ef00382842f58e88f356516a0d)
      xHighEdgeRebin > binsX[nBinsX] ? ibinX_high =  H2D_hist_temp->GetNbinsX()+1 : ibinX_high = H2D_hist_temp->GetXaxis()->FindBin(xHighEdgeRebin) - 1;

      // dividing by prior integral over the xbin to make the weighting done before a probability
      priorBinNorm = H1D_priorSpectrum->IntegralAndError(ibinX_low, ibinX_high, priorBinNormError);
      for(int iBinYFine = ibinY_low; iBinYFine <= ibinY_high; iBinYFine++){
        for(int iBinXFine = ibinX_low; iBinXFine <= ibinX_high; iBinXFine++){
          H2D_hist_content = H2D_hist_temp->GetBinContent(iBinXFine, iBinYFine);
          H2D_hist_contentError = H2D_hist_temp->GetBinError(iBinXFine, iBinYFine);
          H2D_hist_content == 0 ? H2D_hist_contentErrorA = 0 : H2D_hist_contentErrorA = H2D_hist_contentError*H2D_hist_contentError / (H2D_hist_content*H2D_hist_content);
          priorBinNorm == 0     ? H2D_hist_contentErrorB = 0 : H2D_hist_contentErrorB = priorBinNormError*priorBinNormError / (priorBinNorm*priorBinNorm);
          priorBinNorm == 0 ? H2D_hist_content = 0 : H2D_hist_content = H2D_hist_content * 1./priorBinNorm; // do I really give the value 0 if denominator is 0 ?
          H2D_hist_temp->SetBinContent(iBinXFine, iBinYFine, H2D_hist_content);
          H2D_hist_temp->SetBinError(iBinXFine, iBinYFine, sqrt(H2D_hist_content*H2D_hist_content * (H2D_hist_contentErrorA + H2D_hist_contentErrorB))); // sigma(A/B)2 / (A/B) = sigma(A)2 /A2 + sigma(B)2 /B2
          // cout << "iBinXFine = " << iBinXFine << ", iBinYFine = " << iBinYFine << "         --------          H2D_histPostReweight_contentError = " << H2D_hist_contentError << ", priorBinNormError = " << priorBinNormError << ", H2D_histPreReweight_contentError = " << debugTest11 << endl;
          // cout << "                                         --------          H2D_histReweighted_content = " << H2D_hist_content << ", priorBinNorm = " << priorBinNorm << ", H2D_histPreReweight_content = " << debugTest12 << endl;
        }
      }
      hist2D_rebinned_content = H2D_hist_temp->IntegralAndError(ibinX_low, ibinX_high, ibinY_low, ibinY_high, hist2D_rebinned_error);
      H2D_hist_rebinned.SetBinContent(iBinX, iBinY, hist2D_rebinned_content);
      H2D_hist_rebinned.SetBinError(iBinX, iBinY, hist2D_rebinned_error);
      if (debug == true) {cout << "ibinX_low = " << ibinX_low << ", ibinX_high = " << ibinX_high << ", ibinY_low = " << ibinY_low << ", ibinY_high = " << ibinY_high << "         --------          H2D_hist_rebinned(" << iBinX <<", "<< iBinY <<") = " << H2D_hist_rebinned.GetBinContent(iBinX, iBinY) << "                      --- priorBinNorm = " << priorBinNorm << "                                (xLowEdgeRebin, xHighEdgeRebin) = (" << xLowEdgeRebin << ", " << xHighEdgeRebin << "), (yLowEdgeRebin, yHighEdgeRebin) = (" << yLowEdgeRebin << ", " << yHighEdgeRebin << ")" << endl;}
    }
  }


  if (debug == true) {
  // il semble se passer quelque chose à 30GeV et 120GeV, limites des recbins;
  // en fait, il y a des démarcations pour chaque bin de recbins; elles sont juste plus visible pour les grand bins puisque le steeply falling spectrum as plus de temps pour diminuer lorsque l'on va de gauche à droite
  //  ça me semble ok et normal, pas pathologique

  TString* pdfName_preRebin = new TString("H2D_hist_temp");
  TString* pdfName_preRebin_logz = new TString("H2D_hist_temp_logz");

  TString textContext_preRebin((TString)"");

  // Draw_TH2_Histograms(H2D_jetPtResponseMatrix_detectorResponse, centralityLegend, nCentralityBins, textContext, pdfName, texPtJetRecX, texPtJetGenX, texCollisionDataInfo, "logz");


  Draw_TH2_Histogram(H2D_hist_temp, textContext_preRebin, pdfName_preRebin, texPtJetRecX, texPtJetGenX, texCollisionDataInfo, drawnWindowAuto, "");
  Draw_TH2_Histogram(H2D_hist_temp, textContext_preRebin, pdfName_preRebin_logz, texPtJetRecX, texPtJetGenX, texCollisionDataInfo, drawnWindowAuto, "logz");
  }

  // cout << "##############################################" << endl;
  // cout << "######## taking care of truncation ###########" << endl;
  // cout << "##############################################" << endl;
  // // Store gen truncations in x-underflow and overflow for each pT_gen (y-axis); needed for NormaliseYSlicesAsProbaDensity
  // // I don't like the fact I couldn't add it to the above, but whatever, for now
  // double hist2D_rebinned_underflowContent, hist2D_rebinned_overflowContent;
  // double hist2D_rebinned_underflowError, hist2D_rebinned_underflowErrorA, hist2D_rebinned_underflowErrorB;
  // double hist2D_rebinned_overflowError, hist2D_rebinned_overflowErrorA, hist2D_rebinned_overflowErrorB;
  // double testDebug;
  // // truncation in each y slice
  // // for(int iBinY = 1; iBinY <= nBinsY; iBinY++){
  // for(int iBinY = 0; iBinY <= nBinsY+1; iBinY++){
  //   ibinY_low = H2D_hist->GetYaxis()->FindBin(H2D_hist_rebinned.GetYaxis()->GetBinLowEdge(iBinY) );
  //   ibinY_high = H2D_hist->GetYaxis()->FindBin(H2D_hist_rebinned.GetYaxis()->GetBinLowEdge(iBinY+1) ) - 1;
  //   ibinX_low = H2D_hist->GetXaxis()->FindBin(H2D_hist_rebinned.GetXaxis()->GetBinLowEdge(1) ); // according to definition of getbinlowedge, it does work for under/overflows by giving the low edge of bin 1 minus the bin width (calculated assuming bins of equal sizes in all the axis); in fact, it works regardless of the value of iBin even for Nbin+2 or more (https://root.cern.ch/doc/master/classTAxis.html#a6d2d17ef00382842f58e88f356516a0d)
  //   ibinX_high = H2D_hist->GetXaxis()->FindBin(H2D_hist_rebinned.GetXaxis()->GetBinLowEdge(nBinsX+1) ); // according to definition of getbinlowedge, it does work for under/overflows by giving the low edge of bin 1 minus the bin width (calculated assuming bins of equal sizes in all the axis); in fact, it works regardless of the value of iBin even for Nbin+2 or more (https://root.cern.ch/doc/master/classTAxis.html#a6d2d17ef00382842f58e88f356516a0d)


  //   hist2D_rebinned_underflowContent = H2D_hist->IntegralAndError(0, ibinX_low, ibinY_low, ibinY_high, hist2D_rebinned_underflowError);
  //   priorBinNorm = H1D_priorSpectrum->IntegralAndError(0, ibinX_low, priorBinNormError);
  //   testDebug = hist2D_rebinned_underflowError;
  //   hist2D_rebinned_underflowContent == 0 ? hist2D_rebinned_underflowErrorA = 0 : hist2D_rebinned_underflowErrorA = hist2D_rebinned_underflowError*hist2D_rebinned_underflowError / (hist2D_rebinned_underflowContent*hist2D_rebinned_underflowContent);
  //   priorBinNorm == 0                     ? hist2D_rebinned_underflowErrorB = 0 : hist2D_rebinned_underflowErrorB = priorBinNormError*priorBinNormError / (priorBinNorm*priorBinNorm);
  //   priorBinNorm == 0 ? hist2D_rebinned_underflowContent = 0 : hist2D_rebinned_underflowContent = hist2D_rebinned_underflowContent * 1./priorBinNorm; // do I really give the value 0 if denominator is 0 ?
  //   H2D_hist_rebinned.SetBinContent(0, iBinY, hist2D_rebinned_underflowContent);
  //   H2D_hist_rebinned.SetBinError(0, iBinY, sqrt(hist2D_rebinned_underflowContent*hist2D_rebinned_underflowContent * (hist2D_rebinned_underflowErrorA + hist2D_rebinned_underflowErrorB))); // sigma(A/B)2 / (A/B) = sigma(A)2 /A2 + sigma(B)2 /B2

  //   // cout << "xunderflow at ibiny="<< iBinY<< ": error = " << H2D_hist_rebinned.GetBinError(0, iBinY) << endl;
  //   // cout << "                                   priorBinNormError = " << priorBinNormError << ", hist2DintegralError = " << testDebug << endl;


  //   hist2D_rebinned_overflowContent = H2D_hist->IntegralAndError(ibinX_high, H2D_hist->GetNbinsX()+1, ibinY_low, ibinY_high, hist2D_rebinned_overflowError);
  //   priorBinNorm = H1D_priorSpectrum->IntegralAndError(ibinX_high, H2D_hist->GetNbinsX()+1, priorBinNormError);
  //   hist2D_rebinned_overflowContent == 0 ? hist2D_rebinned_overflowErrorA = 0 : hist2D_rebinned_overflowErrorA = hist2D_rebinned_overflowError*hist2D_rebinned_overflowError / (hist2D_rebinned_overflowContent*hist2D_rebinned_overflowContent);
  //   priorBinNorm == 0                    ? hist2D_rebinned_overflowErrorB = 0 : hist2D_rebinned_overflowErrorB = priorBinNormError*priorBinNormError / (priorBinNorm*priorBinNorm);
  //   priorBinNorm == 0 ? hist2D_rebinned_overflowContent = 0 : hist2D_rebinned_overflowContent = hist2D_rebinned_overflowContent * 1./priorBinNorm; // do I really give the value 0 if denominator is 0 ?
  //   if (priorBinNorm == 0) { cout << "WARNING WARNING priorBinNorm = 0 in iBinY = " << iBinY << endl;}
  //   H2D_hist_rebinned.SetBinContent(nBinsX+1, iBinY, hist2D_rebinned_overflowContent);
  //   H2D_hist_rebinned.SetBinError(nBinsX+1, iBinY, sqrt(hist2D_rebinned_overflowContent*hist2D_rebinned_overflowContent * (hist2D_rebinned_overflowErrorA + hist2D_rebinned_overflowErrorB))); // sigma(A/B)2 / (A/B) = sigma(A)2 /A2 + sigma(B)2 /B2
    
  //   // cout << "xoverflow at ibiny="<< iBinY<< ": error = " << H2D_hist_rebinned.GetBinError(nBinsX+1, iBinY) << endl;
  // }
  // // truncation in each x slice
  // cout << "truncation in each x slice (for RebinVariableBins2D_PriorWeightedBinMerging) hasn't been really checked, need to do so" << endl;
  // // for(int iBinX = 1; iBinX <= nBinsX; iBinX++){
  // for(int iBinX = 0; iBinX <= nBinsX+1; iBinX++){
  //   ibinY_low = H2D_hist->GetYaxis()->FindBin(H2D_hist_rebinned.GetYaxis()->GetBinLowEdge(1) ); // according to definition of getbinlowedge, it does work for under/overflows by giving the low edge of bin 1 minus the bin width (calculated assuming bins of equal sizes in all the axis); in fact, it works regardless of the value of iBin even for Nbin+2 or more (https://root.cern.ch/doc/master/classTAxis.html#a6d2d17ef00382842f58e88f356516a0d)
  //   ibinY_high = H2D_hist->GetYaxis()->FindBin(H2D_hist_rebinned.GetYaxis()->GetBinLowEdge(nBinsY+1) ) - 1; // according to definition of getbinlowedge, it does work for under/overflows by giving the low edge of bin 1 minus the bin width (calculated assuming bins of equal sizes in all the axis); in fact, it works regardless of the value of iBin even for Nbin+2 or more (https://root.cern.ch/doc/master/classTAxis.html#a6d2d17ef00382842f58e88f356516a0d)
  //   ibinX_low = H2D_hist->GetXaxis()->FindBin(H2D_hist_rebinned.GetXaxis()->GetBinLowEdge(iBinX) );
  //   ibinX_high = H2D_hist->GetXaxis()->FindBin(H2D_hist_rebinned.GetXaxis()->GetBinLowEdge(iBinX+1) );


  //   hist2D_rebinned_underflowContent = H2D_hist->IntegralAndError(ibinX_low, ibinX_low, 0, ibinY_low, hist2D_rebinned_underflowError);
  //   priorBinNorm = H1D_priorSpectrum->IntegralAndError(ibinX_low, ibinX_low, priorBinNormError);
  //   testDebug = hist2D_rebinned_underflowError;
  //   hist2D_rebinned_underflowContent == 0 ? hist2D_rebinned_underflowErrorA = 0 : hist2D_rebinned_underflowErrorA = hist2D_rebinned_underflowError*hist2D_rebinned_underflowError / (hist2D_rebinned_underflowContent*hist2D_rebinned_underflowContent);
  //   priorBinNorm == 0                     ? hist2D_rebinned_underflowErrorB = 0 : hist2D_rebinned_underflowErrorB = priorBinNormError*priorBinNormError / (priorBinNorm*priorBinNorm);
  //   priorBinNorm == 0 ? hist2D_rebinned_underflowContent = 0 : hist2D_rebinned_underflowContent = hist2D_rebinned_underflowContent * 1./priorBinNorm; // do I really give the value 0 if denominator is 0 ?
  //   H2D_hist_rebinned.SetBinContent(ibinX_low, 0, hist2D_rebinned_underflowContent);
  //   H2D_hist_rebinned.SetBinError(ibinX_low, 0, sqrt(hist2D_rebinned_underflowContent*hist2D_rebinned_underflowContent * (hist2D_rebinned_underflowErrorA + hist2D_rebinned_underflowErrorB))); // sigma(A/B)2 / (A/B) = sigma(A)2 /A2 + sigma(B)2 /B2

  //   // cout << "xunderflow at ibiny="<< iBinY<< ": error = " << H2D_hist_rebinned.GetBinError(0, iBinY) << endl;
  //   // cout << "                                   priorBinNormError = " << priorBinNormError << ", hist2DintegralError = " << testDebug << endl;


  //   hist2D_rebinned_overflowContent = H2D_hist->IntegralAndError(ibinX_low, ibinX_low, ibinY_high, H2D_hist->GetNbinsY()+1, hist2D_rebinned_overflowError);
  //   priorBinNorm = H1D_priorSpectrum->IntegralAndError(ibinX_low, ibinX_low, priorBinNormError);
  //   hist2D_rebinned_overflowContent == 0 ? hist2D_rebinned_overflowErrorA = 0 : hist2D_rebinned_overflowErrorA = hist2D_rebinned_overflowError*hist2D_rebinned_overflowError / (hist2D_rebinned_overflowContent*hist2D_rebinned_overflowContent);
  //   priorBinNorm == 0                    ? hist2D_rebinned_overflowErrorB = 0 : hist2D_rebinned_overflowErrorB = priorBinNormError*priorBinNormError / (priorBinNorm*priorBinNorm);
  //   priorBinNorm == 0 ? hist2D_rebinned_overflowContent = 0 : hist2D_rebinned_overflowContent = hist2D_rebinned_overflowContent * 1./priorBinNorm; // do I really give the value 0 if denominator is 0 ?
  //   H2D_hist_rebinned.SetBinContent(ibinX_low, nBinsY+1, hist2D_rebinned_overflowContent);
  //   H2D_hist_rebinned.SetBinError(ibinX_low, nBinsY+1, sqrt(hist2D_rebinned_overflowContent*hist2D_rebinned_overflowContent * (hist2D_rebinned_overflowErrorA + hist2D_rebinned_overflowErrorB))); // sigma(A/B)2 / (A/B) = sigma(A)2 /A2 + sigma(B)2 /B2
    
  //   // cout << "xoverflow at ibiny="<< iBinY<< ": error = " << H2D_hist_rebinned.GetBinError(nBinsX+1, iBinY) << endl;
  // }



  // normalisation
  // H2D_hist_rebinned.ResetStats(); // setbincontent interacts weirdly with getentries(); resetstats makes it so that getentries gives the sum of bin contents correctly
  // if (debug == true) {cout << "RebinVariableBins2D - what of the errors" << endl;}
  // if (debug == true) {cout << "H2D_hist nbins X =" << H2D_hist->GetNbinsX() << endl;}
  // if (debug == true) {cout << "H2D_hist_rebinned nbins X =" << H2D_hist_rebinned.GetNbinsX() << endl;}

  std::stringstream ss;
  ss << H2D_hist->GetName() << "_RebinVariableBins2D";
  TString histName((TString)ss.str());

  if (debug == true) {cout << "___________ end RebinVariableBins2D" << endl;}

  return H2D_hist_rebinned;
}
TH2D RebinVariableBins2D_PriorWeightedBinMerging(TH2D* H2D_hist, int nBinsX, int nBinsY, double* binsX, double* binsY, TH1D* H1D_priorSpectrum){
  return RebinVariableBins2D_PriorWeightedBinMerging(H2D_hist, nBinsX, nBinsY, binsX, binsY, H1D_priorSpectrum, false);
}




void NormaliseYSlicesAsProbaDensity(TH2D* H2D_hist){
  cout << "Normalising y-slices" << endl;
  // Normalisation of the H2D_hist 2D histogram: each pt gen slice is normalised to unity, takes underflow into account if it's present
  double genSliceNorm = 1;
  double genSliceNormError;
  // double dpT;

  // for(int iBinY = 0; iBinY <= H2D_hist->GetNbinsY()+1; iBinY++){ // 0 and n+1 take underflow and overflow into account
  //   for(int iBinX = 0; iBinX <= H2D_hist->GetNbinsX()+1; iBinX++){ // 0 and n+1 take underflow and overflow into account
  //   if (iBinY<10 && iBinX<10) {
  //     cout << "H2D_hist(" << iBinX <<", "<< iBinY <<") = " << H2D_hist->GetBinContent(iBinX, iBinY) << endl;
  //   }
  //   }
  // }
  double binContent, binError, binErrorA, binErrorB;
  for(int iBinY = 0; iBinY <= H2D_hist->GetNbinsY()+1; iBinY++){ // 0 and n+1 take underflow and overflow into account
    genSliceNorm = H2D_hist->IntegralAndError(0, H2D_hist->GetNbinsX()+1, iBinY, iBinY, genSliceNormError); // what I want ideally
    // genSliceNorm = H2D_hist->IntegralAndError(1, H2D_hist->GetNbinsX(), iBinY, iBinY, genSliceNormError);
    // cout << "iBinY = " << iBinY << "         --------  genSliceNorm = " << genSliceNorm << endl;
    // genSliceNorm = H2D_hist->Integral(1, H2D_hist->GetNbinsX(), iBinY, iBinY);
    // dpT = H2D_hist->GetYaxis()->GetBinWidth(iBinY);
    // cout << "iBinY = " << iBinY << "         --------          genSliceNorm = " << genSliceNorm << " ----  1/dpT*genSliceNorm = " << genSliceNorm*1./dpT << ", dpT = " << dpT << endl;
    for(int iBinX = 0; iBinX <= H2D_hist->GetNbinsX()+1; iBinX++){ // 0 and n+1 take underflow and overflow into account
      H2D_hist->GetBinContent(iBinX, iBinY) == 0 ? binErrorB = 0 : binErrorB = H2D_hist->GetBinError(iBinX, iBinY)*H2D_hist->GetBinError(iBinX, iBinY) / (H2D_hist->GetBinContent(iBinX, iBinY)*H2D_hist->GetBinContent(iBinX, iBinY));
      genSliceNorm == 0                          ? binErrorA = 0 : binErrorA = genSliceNormError*genSliceNormError / (genSliceNorm*genSliceNorm);
      genSliceNorm == 0 ? binContent = 0 : binContent = H2D_hist->GetBinContent(iBinX, iBinY) *1./genSliceNorm; // do I really give the value 0 if denominator is 0 ? 
      H2D_hist->SetBinContent(iBinX, iBinY, binContent);
      H2D_hist->SetBinError(iBinX, iBinY, sqrt(binContent*binContent * (binErrorA + binErrorB))); // sigma(A/B)2 / (A/B) = sigma(A)2 /A2 + sigma(B)2 /B2
      // cout << "iBinX = " << iBinX << ", iBinY = " << iBinY << "         --------  H2D_hist->GetBinContent(iBinX, iBinY) = " << H2D_hist->GetBinContent(iBinX, iBinY) << endl;
    }
    cout << "iBinY = " << iBinY << " -------- checking kinematic efficiency: Integral(1, N) = " << H2D_hist->Integral(1, H2D_hist->GetNbinsX(), iBinY, iBinY) << endl;
    cout << "iBinY = " << iBinY << " -------- checking kinematic efficiency: Integral(0,-1) = " << H2D_hist->Integral(0, -1, iBinY, iBinY) << endl;
    cout << "genSliceNorm = " << genSliceNorm << ", underflow = " << H2D_hist->GetBinContent(0, iBinY) << ", overflow = " << H2D_hist->GetBinContent(H2D_hist->GetNbinsX()+1, iBinY) << endl;

  }
} 



void WeightMatrixWithPrior(TH2D* H2D_hist, TH1D* priorSpectrum){
  cout << "Weighting y-slices" << endl;
  // ACTUALLY I DONT THINK THIS SHOULD BE USED, THE PRIOR SEEMS TO BE SET UP DIFFERENTLY
  // looks like I shouldn't wieght the under and overflows, or it puts too much weight on the overflow
  // maybe I should make the kinematic efficiency matrix before any re-weighting

  // actually what's important here is that the y-projection is the prior

  // weights y-slices with priorSpectrum 
  double genSliceNormError;
  // double genSliceNorm = priorSpectrum->IntegralAndError(0, -1, genSliceNormError);
  double genSliceNorm = priorSpectrum->IntegralAndError(1, priorSpectrum->GetNbinsX(), genSliceNormError);

  double binContent, binError, binErrorA, binErrorB;

  // multiplication
  // for(int iBinY = 0; iBinY <= H2D_hist->GetNbinsY()+1; iBinY++){ // 0 and n+1 take underflow and overflow into account
  //   for(int iBinX = 0; iBinX <= H2D_hist->GetNbinsX()+1; iBinX++){ // 0 and n+1 take underflow and overflow into account
  for(int iBinY = 1; iBinY <= H2D_hist->GetNbinsY(); iBinY++){ // 0 and n+1 take underflow and overflow into account
    for(int iBinX = 1; iBinX <= H2D_hist->GetNbinsX(); iBinX++){ // 0 and n+1 take underflow and overflow into account
      binContent = H2D_hist->GetBinContent(iBinX, iBinY) * priorSpectrum->GetBinContent(iBinY); // do I really give the value 0 if denominator is 0 ? 
      H2D_hist->SetBinContent(iBinX, iBinY, binContent);
      H2D_hist->SetBinError(iBinX, iBinY, sqrt(priorSpectrum->GetBinContent(iBinY)*priorSpectrum->GetBinContent(iBinY) * H2D_hist->GetBinError(iBinX, iBinY)*H2D_hist->GetBinError(iBinX, iBinY) + binContent*binContent * priorSpectrum->GetBinError(iBinY)*priorSpectrum->GetBinError(iBinY))); // sigma(A/B)2 / (A/B) = sigma(A)2 /A2 + sigma(B)2 /B2
      // cout << "iBinX = " << iBinX << ", iBinY = " << iBinY << "         --------  H2D_hist->GetBinContent(iBinX, iBinY) = " << H2D_hist->GetBinContent(iBinX, iBinY) << endl;
    }
  }

  // note: the matrix is not normalised to 1 in y-slices yet!

  // division
  // for(int iBinY = 0; iBinY <= H2D_hist->GetNbinsY()+1; iBinY++){ // 0 and n+1 take underflow and overflow into account
  //   for(int iBinX = 0; iBinX <= H2D_hist->GetNbinsX()+1; iBinX++){ // 0 and n+1 take underflow and overflow into account
  for(int iBinY = 1; iBinY <= H2D_hist->GetNbinsY(); iBinY++){ // 0 and n+1 take underflow and overflow into account
    for(int iBinX = 1; iBinX <= H2D_hist->GetNbinsX(); iBinX++){ // 0 and n+1 take underflow and overflow into account
      H2D_hist->GetBinContent(iBinX, iBinY) == 0 ? binErrorB = 0 : binErrorB = H2D_hist->GetBinError(iBinX, iBinY)*H2D_hist->GetBinError(iBinX, iBinY) / (H2D_hist->GetBinContent(iBinX, iBinY)*H2D_hist->GetBinContent(iBinX, iBinY));
      genSliceNorm == 0                          ? binErrorA = 0 : binErrorA = genSliceNormError*genSliceNormError / (genSliceNorm*genSliceNorm);
      genSliceNorm == 0 ? binContent = 0 : binContent = H2D_hist->GetBinContent(iBinX, iBinY) *1./genSliceNorm; // do I really give the value 0 if denominator is 0 ? 
      H2D_hist->SetBinContent(iBinX, iBinY, binContent);
      H2D_hist->SetBinError(iBinX, iBinY, sqrt(binContent*binContent * (binErrorA + binErrorB))); // sigma(A/B)2 / (A/B) = sigma(A)2 /A2 + sigma(B)2 /B2
      // cout << "iBinX = " << iBinX << ", iBinY = " << iBinY << "         --------  H2D_hist->GetBinContent(iBinX, iBinY) = " << H2D_hist->GetBinContent(iBinX, iBinY) << endl;
    }
    cout << "iBinY = " << iBinY << " -------- checking kinematic efficiency: Integral(1, N) = " << H2D_hist->Integral(1, H2D_hist->GetNbinsX(), iBinY, iBinY) << endl;
    cout << "iBinY = " << iBinY << " -------- checking kinematic efficiency: Integral(0,-1) = " << H2D_hist->Integral(0, -1, iBinY, iBinY) << endl;
    cout << "genSliceNorm = " << genSliceNorm << ", underflow = " << H2D_hist->GetBinContent(0, iBinY) << ", overflow = " << H2D_hist->GetBinContent(H2D_hist->GetNbinsX()+1, iBinY) << endl;
  }
} 


// TH2* AliAnaChargedJetResponseMaker::NormalizeResponsMatrixYaxisWithPrior(TH2 *h2RM, TH1 *hPrior) {
//   //
//   // Normalize such that the Y projection is the prior
//   //

//   //  TH1D *hProjRespY = (TH1D*)h2RM->ProjectionY("hProjRespY");
//   double intPrior = hPrior->Integral();//"width");
//   for (Int_t jbin = 1; jbin <= h2RM->GetNbinsY(); jbin++) {
//     //    double corr = hPrior->GetBinContent(jbin)/hProjRespY->GetBinContent(jbin);
//       for (Int_t ibin = 1; ibin <= h2RM->GetNbinsX(); ibin++) {
// 	double content = h2RM->GetBinContent(ibin,jbin);
// 	//	h2RM->SetBinContent(ibin,jbin,content*corr);
// 	h2RM->SetBinContent(ibin,jbin,hPrior->GetBinContent(jbin)/hPrior->GetBinWidth(jbin)/intPrior*content);
//     }
//   }

//   return h2RM;
// }

TH2D GetTransposeHistogram(TH2D* inputHist){
  int nBinsX = inputHist->GetNbinsX();
  int nBinsY = inputHist->GetNbinsY();
  std::vector<double> vectorBinsX = GetTH1Bins(inputHist->ProjectionX(inputHist->GetName()+(TString)"projX", 1, nBinsX));
  std::vector<double> vectorBinsY = GetTH1Bins(inputHist->ProjectionY(inputHist->GetName()+(TString)"projY", 1, nBinsY));
  double* binsX = &vectorBinsX[0];
  double* binsY = &vectorBinsY[0];

  TH2D transposedHist("inputHist_transposed", "inputHist_transposed", nBinsY, binsY, nBinsX, binsX);
  transposedHist.Sumw2();
  for(int iBinX = 0; iBinX <= nBinsX+1; iBinX++){ // 0 and n+1 take underflow and overflow into account
    for(int iBinY = 0; iBinY <= nBinsY+1; iBinY++){ // 0 and n+1 take underflow and overflow into account
      transposedHist.SetBinContent(iBinY, iBinX, inputHist->GetBinContent(iBinX, iBinY));
      transposedHist.SetBinError(iBinY, iBinX, inputHist->GetBinError(iBinX, iBinY));
    }
  }
  return transposedHist;
}

TH2D GetMatrixProductTH2xTH2(TH2D* histA, TH2D* histB){
  int nBinsX_A = histA->GetNbinsX();
  int nBinsY_A = histA->GetNbinsY();
  int nBinsX_B = histB->GetNbinsX();
  int nBinsY_B = histB->GetNbinsY();

    cout << "doing product matrix A times matrix B, with A(" << nBinsY_A << " x " << nBinsX_A << ") and B(" << nBinsY_B << " x " << nBinsX_B << ")" << endl; 
  if (nBinsY_A != nBinsX_B) {
    cout << "#########################################################################################" << endl;
    cout << "###################### MATRIX PRODUCT IMPOSSIBLE DUE TO DIMENSIONS ######################" << endl;
    cout << "#########################################################################################" << endl;
    cout << "matrix A times B, with A(" << nBinsY_A << " x " << nBinsX_A << ") and B(" << nBinsY_B << " x " << nBinsX_B << ")" << endl; 
  }

  int nBinsK = nBinsY_A;
  
  std::vector<double> vectorBinsX_AB = GetTH1Bins(histB->ProjectionX(histB->GetName()+(TString)"projX", 1, nBinsX_B));
  std::vector<double> vectorBinsY_AB = GetTH1Bins(histA->ProjectionY(histA->GetName()+(TString)"projY", 1, nBinsY_A));
  double* binsX_AB = &vectorBinsX_AB[0];
  double* binsY_AB = &vectorBinsY_AB[0];
  int nBinsX_AB = nBinsX_B;
  int nBinsY_AB = nBinsY_A;

  TH2D histAB("histAB", "histAB", nBinsX_B, binsX_AB, nBinsY_A, binsY_AB);
  histAB.Sumw2();
  // Roounfold sees the element of a hist histogram(i,j) as the element R(i,j) of a matrix; thus no transformation has to be done for the indices or axes of the histograms when doing the multiplication

  double productContent_iBinX_iBinY;
  double productError2_iBinX_iBinY;
  // for(int iBinX = 1; iBinX <= nBinsX_AB; iBinX++){
  for(int iBinX = 0; iBinX <= nBinsX_AB+1; iBinX++){ // 0 and n+1 take underflow and overflow into account
    // for(int iBinY = 1; iBinY <= nBinsY_AB; iBinY++){
    for(int iBinY = 0; iBinY <= nBinsY_AB+1; iBinY++){ // 0 and n+1 take underflow and overflow into account
      productContent_iBinX_iBinY = 0;
      productError2_iBinX_iBinY = 0;
      if (0 < iBinX && iBinX < nBinsX_AB+1 && 0 < iBinY && iBinY < nBinsY_AB+1) { // reason we do this, is because we don't want the under/ofverflows spread inside the matrix when we multiply, which happens when we do ; by separating it, we can still have the under/overflows calculated correctly in the else case, so that we can use them for the kinematic efficiency
        // ideally I want to take the over/underflows into account; or at least have them save how much of the distrib is cut
        // maybe I have them be recalculated with the product, but I don't let them be part of the sum for the non-under/overflow parts of the final matrix
        for(int iBinK = 1; iBinK <= nBinsK; iBinK++){ // 0 and n+1 take underflow and overflow into account
        // for(int iBinK = 0; iBinK <= nBinsK+1; iBinK++){ // 0 and n+1 take underflow and overflow into account
          productContent_iBinX_iBinY += histA->GetBinContent(iBinX, iBinK) * histB->GetBinContent(iBinK, iBinY); 
          productError2_iBinX_iBinY += pow(histB->GetBinContent(iBinK, iBinY), 2)*pow(histA->GetBinError(iBinX, iBinK), 2) + pow(histA->GetBinContent(iBinX, iBinK), 2)*pow(histB->GetBinError(iBinK, iBinY), 2); // simple sigma_ab = a2sigma_a2 + b2sigma_b2 ; that assumes there are no correlations; here it s background fluct from PbPB sim, and detector effects from a pp sim, so we can maybe say theyre not correlated    
          // cout << "iBinX = " << iBinX << ", iBinY = " << iBinY << "         --------  iBinK = " << iBinK << "        B = " << histB->GetBinContent(iBinX, iBinK) << ", sigmaA = " << histA->GetBinError(iBinX, iBinK) << ", A = " << histA->GetBinContent(iBinK, iBinY) << ", sigmaB = " << histB->GetBinError(iBinK, iBinY) << endl;
        }
      } else {
        // for(int iBinK = 1; iBinK <= nBinsK; iBinK++){ // 0 and n+1 take underflow and overflow into account
        for(int iBinK = 0; iBinK <= nBinsK+1; iBinK++){ // 0 and n+1 take underflow and overflow into account
          productContent_iBinX_iBinY += histA->GetBinContent(iBinX, iBinK) * histB->GetBinContent(iBinK, iBinY); 
          productError2_iBinX_iBinY += pow(histB->GetBinContent(iBinK, iBinY), 2)*pow(histA->GetBinError(iBinX, iBinK), 2) + pow(histA->GetBinContent(iBinX, iBinK), 2)*pow(histB->GetBinError(iBinK, iBinY), 2); // simple sigma_ab = a2sigma_a2 + b2sigma_b2 ; that assumes there are no correlations; here it s background fluct from PbPB sim, and detector effects from a pp sim, so we can maybe say theyre not correlated    
          // cout << "iBinX = " << iBinX << ", iBinY = " << iBinY << "         --------  iBinK = " << iBinK << "        B = " << histB->GetBinContent(iBinX, iBinK) << ", sigmaA = " << histA->GetBinError(iBinX, iBinK) << ", A = " << histA->GetBinContent(iBinK, iBinY) << ", sigmaB = " << histB->GetBinError(iBinK, iBinY) << endl;
        }
      }

      histAB.SetBinContent(iBinX, iBinY, productContent_iBinX_iBinY);
      histAB.SetBinError(iBinX, iBinY, sqrt(productError2_iBinX_iBinY));  
      // cout << "iBinX = " << iBinX << ", iBinY = " << iBinY << "         --------                    productContent_iBinX_iBinY = " << productContent_iBinX_iBinY << ", productError2_iBinX_iBinY = " << productError2_iBinX_iBinY << endl;
    }
  }
  return histAB;


}

TH1D GetMatrixVectorProductTH2xTH1(TH2D* histA, TH1D* histU){
  // bit different from TH2 x TH2 as TH1 is a (1,N) TH2 and not a (N,1) one
  int nBinsX_A = histA->GetNbinsX();
  int nBinsY_A = histA->GetNbinsY();
  int nBinsX_U = histU->GetNbinsX();

  cout << "doing product matrix A times vector U, with A(" << nBinsX_A << " x " << nBinsY_A << ") and U(" << nBinsX_U << " x 1)" << endl; 
  if (nBinsY_A != nBinsX_U) {
    cout << "#########################################################################################" << endl;
    cout << "###################### MATRIX PRODUCT IMPOSSIBLE DUE TO DIMENSIONS ######################" << endl;
    cout << "#########################################################################################" << endl;
  }

  int nBinsK = nBinsY_A;
  // int nBinsX_AB = nBinsX_B;
  // int nBinsY_AB = nBinsY_A;
  
  std::vector<double> vectorBinsX_AU = GetTH1Bins(histA->ProjectionX(histA->GetName()+(TString)"projY", 1, nBinsY_A));
  double* binsX_AU = &vectorBinsX_AU[0];
  int nBinsX_AU = nBinsX_A;

  TH1D histAU("histAU", "histAU", nBinsX_AU, binsX_AU);
  histAU.Sumw2();

  // Roounfold sees the element of a hist histogram(i,j) as the element R(i,j) of a matrix; thus no transformation has to be done for the indices or axes of the histograms when doing the multiplication

  double productContent_iBinX;
  double productError2_iBinX;
  // for(int iBinX = 1; iBinX <= nBinsX_AU; iBinX++){
  for(int iBinX = 0; iBinX <= nBinsX_AU+1; iBinX++){ // 0 and n+1 take underflow and overflow into account
    productContent_iBinX = 0;
    productError2_iBinX = 0;
    // for(int iBinK = 1; iBinK <= nBinsK; iBinK++){
    for(int iBinK = 0; iBinK <= nBinsK+1; iBinK++){ // 0 and n+1 take underflow and overflow into account
      productContent_iBinX += histA->GetBinContent(iBinX, iBinK) * histU->GetBinContent(iBinK); 
      productError2_iBinX += pow(histU->GetBinContent(iBinK), 2)*pow(histA->GetBinError(iBinX, iBinK), 2) + pow(histA->GetBinContent(iBinX, iBinK), 2)*pow(histU->GetBinError(iBinK), 2); // simple sigma_ab = a2sigma_a2 + b2sigma_b2 ; that assumes there are no correlations; here it s background fluct from PbPB sim, and detector effects from a pp sim, so we can maybe say theyre not correlated          
    }
    histAU.SetBinContent(iBinX, productContent_iBinX);
    histAU.SetBinError(iBinX, sqrt(productError2_iBinX));  
  }
  return histAU;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////// Histogram Context /////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

TString contextCustomThreeFields(TString mainContext, TString secondaryContext, TString tertiaryContext, __attribute__ ((unused)) const char options[]){
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

TString contextDatasetCompAndRadius(TString* mainContext, float jetRadius, __attribute__ ((unused)) const char options[]){
  TString texcontextDatasetCompAndRadius;
  texcontextDatasetCompAndRadius = "#splitline{"+*mainContext+"}{"+contextJetRadius(jetRadius)+"}";

  return texcontextDatasetCompAndRadius;
}

TString contextDatasetComp(TString* mainContext, __attribute__ ((unused)) const char options[]){
  TString texcontextDatasetComp;
  texcontextDatasetComp = *mainContext;

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

void IterationLegend(TString* iterationLegend, int nIterationmax){
  std::stringstream ss;
  ss.precision(2);
  for(int iIteration = 0; iIteration < nIterationmax; iIteration++){
    ss << "k_{unfold} = " << iIteration;
    iterationLegend[iIteration] = (TString)ss.str();
    ss.str("");
    ss.clear();
  }
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

void Draw_TH1_Histogram(TH1D* histogram, TString Context, TString* pdfName, TString* &texXtitle, TString* &texYtitle, TString* texCollisionDataInfo, std::array<std::array<float, 2>, 2> drawnWindow, const char options[]) {
  TH1D* singleHistArray[1] = {histogram};
  TString dummyLegend[1] = {(TString)""};
  int dummyCollectionSize = 1;
  Draw_TH1_Histograms_in_one(singleHistArray, dummyLegend, dummyCollectionSize, Context, pdfName, texXtitle, texYtitle, texCollisionDataInfo, drawnWindow, options);


  // for(int iCentralityBin = 0; iCentralityBin < nCentralityBins; iCentralityBin++){
  //   cout << "histogram->GetBinContent(iCentralityBin) = " << histogram->GetBinContent(iCentralityBin) << endl;
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

  cout << "---------------- TH2 drawing; getting min and max test:" << endl;
  for (int i = 0; i < collectionSize; i++) {
    histograms_collection[i]->GetZaxis()->SetRangeUser(histograms_collection[i]->GetMinimum(GLOBAL_epsilon), histograms_collection[i]->GetMaximum());
    cout << "min = " << histograms_collection[i]->GetMinimum(GLOBAL_epsilon) << ", max = " << histograms_collection[i]->GetMaximum() << endl;
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

void Draw_TH2_Histograms(TH2D** histograms_collection, const TString* legendList_string, int collectionSize, TString Context, TString* pdfName, TString* &texXtitle, TString* &texYtitle, TString* texCollisionDataInfo, std::array<std::array<float, 2>, 2> drawnWindow, const char options[]) {
  // is here to make optionalFitCollection an actual optional parameter; Draw_TH1_Histograms_in_one can be called without, and in that case optionalFitCollection is created empty for use by the actual Draw_TH1_Histograms_in_one function; it will only be used if 'options' has fit in it
  TPolyLine* optionalLine;
  Draw_TH2_Histograms(histograms_collection, legendList_string, collectionSize, Context, pdfName, texXtitle, texYtitle, texCollisionDataInfo, drawnWindow, options, optionalLine);
}


void Draw_TH2_Histogram(TH2D* histogram, TString Context, TString* pdfName, TString* &texXtitle, TString* &texYtitle, TString* texCollisionDataInfo, std::array<std::array<float, 2>, 2> drawnWindow, const char options[]) {
  TH2D* singleHistArray[1] = {histogram};
  TString dummyLegend[1] = {(TString)""};
  int dummyCollectionSize = 1;
  Draw_TH2_Histograms(singleHistArray, dummyLegend, dummyCollectionSize, Context, pdfName, texXtitle, texYtitle, texCollisionDataInfo, drawnWindow, options);


  // for(int iCentralityBin = 0; iCentralityBin < nCentralityBins; iCentralityBin++){
  //   cout << "histogram->GetBinContent(iCentralityBin) = " << histogram->GetBinContent(iCentralityBin) << endl;
  // }
}