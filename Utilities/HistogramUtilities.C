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

TH2D RebinVariableBins2D(TH2D* H2D_hist, int nBinsX, int nBinsY, double* binsX, double* binsY, bool debug){
  if (debug == true) {cout << "___________ begin RebinVariableBins2D, nBinsX = " << nBinsX << ", nBinsY = " << nBinsY << endl;}
  // int debugCount = 0;

  TH2D H2D_hist_rebinned("H2D_hist_rebinned", "H2D_hist_rebinned", nBinsX, binsX, nBinsY, binsY);
  // cout << "H2D_hist->GetEntries() = " << H2D_hist->GetEntries() << endl;
  // for(int iBinX = 1; iBinX <= H2D_hist->GetNbinsX(); iBinX++){
  //   for(int iBinY = 1; iBinY <= H2D_hist->GetNbinsY(); iBinY++){
  //     for (int iEntries_bin = 1; iEntries_bin <= H2D_hist->GetBinContent(iBinX, iBinY); iEntries_bin++){ // doing this loop because using fill(x, y, weight = H2D_hist->GetBinContent(iBinX, iBinY)) instead doesn't increment the number of entries by N=weight but by N=1.
  //       H2D_hist_rebinned.Fill(H2D_hist_rebinned.FindBin(H2D_hist->GetXaxis()->GetBinCenter(iBinX)), H2D_hist_rebinned.FindBin(H2D_hist->GetYaxis()->GetBinCenter(iBinY))); // only works if histogram isn't weighted
  //       debugCount += 1;
  //     }
  //   }
  // }
  // cout << "debugCount = " << debugCount << endl;
  int ibinX_low, ibinX_high, ibinY_low, ibinY_high;
  // for(int iBinX = 1; iBinX <= nBinsX; iBinX++){
  for(int iBinX = 0; iBinX <= nBinsX+1; iBinX++){ // 0 and n+1 take underflow and overflow into account
    if (debug == true) {cout << "iBinX = " << iBinX << endl;}
    ibinX_low = H2D_hist->GetXaxis()->FindBin(H2D_hist_rebinned.GetXaxis()->GetBinLowEdge(iBinX) ); // according to definition of getbinlowedge, it does work for under/overflows by giving the low edge of bin 1 minus the bin width (calculated assuming bins of equal sizes in all the axis); in fact, it works regardless of the value of iBin even for Nbin+2 or more (https://root.cern.ch/doc/master/classTAxis.html#a6d2d17ef00382842f58e88f356516a0d)
    ibinX_high = H2D_hist->GetXaxis()->FindBin(H2D_hist_rebinned.GetXaxis()->GetBinLowEdge(iBinX+1)) - 1;
    // for(int iBinY = 1; iBinY <= nBinsY; iBinY++){
    for(int iBinY = 0; iBinY <= nBinsY+1; iBinY++){ // 0 and n+1 take underflow and overflow into account
      ibinY_low = H2D_hist->GetYaxis()->FindBin(H2D_hist_rebinned.GetYaxis()->GetBinLowEdge(iBinY) );
      ibinY_high = H2D_hist->GetYaxis()->FindBin(H2D_hist_rebinned.GetYaxis()->GetBinLowEdge(iBinY+1) ) - 1;
      H2D_hist_rebinned.SetBinContent(iBinX, iBinY, H2D_hist->Integral(ibinX_low, ibinX_high, ibinY_low, ibinY_high));
      if (debug == true) {cout << "ibinX_low = " << ibinX_low << ", ibinX_high = " << ibinX_high << ", ibinY_low = " << ibinY_low << ", ibinY_high = " << ibinY_high << "         --------          H2D_hist_rebinned(iBinX, iBinY) = " << H2D_hist_rebinned.GetBinContent(iBinX, iBinY) << endl;}
    }
  }

  // Store gen truncations in x-underflow for each pT_gen (y-axis); needed for NormaliseYSlicesAsProbaDensity
  // I don't like the fact I couldn't add it to the above, but whatever, for now
  double hist2D_rebinned_underflowContent, hist2D_rebinned_overflowContent;
  double hist2D_rebinned_underflowError, hist2D_rebinned_overflowError;
  for(int iBinY = 0; iBinY <= nBinsY+1; iBinY++){ // 0 and n+1 take underflow and overflow into account
    ibinX_low = H2D_hist->GetXaxis()->FindBin(H2D_hist_rebinned.GetXaxis()->GetBinLowEdge(1) ); // according to definition of getbinlowedge, it does work for under/overflows by giving the low edge of bin 1 minus the bin width (calculated assuming bins of equal sizes in all the axis); in fact, it works regardless of the value of iBin even for Nbin+2 or more (https://root.cern.ch/doc/master/classTAxis.html#a6d2d17ef00382842f58e88f356516a0d)
    hist2D_rebinned_underflowContent = H2D_hist->IntegralAndError(0, ibinX_low, iBinY, iBinY, hist2D_rebinned_underflowError);
    H2D_hist_rebinned.SetBinContent(0, iBinY, hist2D_rebinned_underflowContent);
    H2D_hist_rebinned.SetBinError(0, iBinY, hist2D_rebinned_underflowError);
    cout << "iBinY = " << iBinY << ", ibinX_low = " << ibinX_low <<  "         --------          H2D_hist_rebinned(iBinX, iBinY) = " << H2D_hist_rebinned.GetBinContent(0, iBinY) << endl;

    ibinX_high = H2D_hist->GetXaxis()->FindBin(H2D_hist_rebinned.GetXaxis()->GetBinLowEdge(nBinsX+1) ); // according to definition of getbinlowedge, it does work for under/overflows by giving the low edge of bin 1 minus the bin width (calculated assuming bins of equal sizes in all the axis); in fact, it works regardless of the value of iBin even for Nbin+2 or more (https://root.cern.ch/doc/master/classTAxis.html#a6d2d17ef00382842f58e88f356516a0d)
    hist2D_rebinned_overflowContent = H2D_hist->IntegralAndError(ibinX_high, H2D_hist->GetNbinsX()+1, iBinY, iBinY, hist2D_rebinned_overflowError);
    H2D_hist_rebinned.SetBinContent(nBinsX+1, iBinY, hist2D_rebinned_overflowContent);
    H2D_hist_rebinned.SetBinError(nBinsX+1, iBinY, hist2D_rebinned_overflowError);
    cout << "iBinY = " << iBinY << ", ibinX_high = " << ibinX_high <<  "         --------          H2D_hist_rebinned(iBinX, iBinY) = " << H2D_hist_rebinned.GetBinContent(nBinsX+1, iBinY) << endl;
    
  }

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

TH2D RebinVariableBins2D_ySlicePriorWeighted(TH2D* H2D_hist, int nBinsX, int nBinsY, double* binsX, double* binsY, TH1D* H1D_priorSpectrum, bool debug){
  if (debug == true) {cout << "___________ begin RebinVariableBins2D, nBinsX = " << nBinsX << ", nBinsY = " << nBinsY << endl;}

  TH2D H2D_hist_rebinned("H2D_hist_rebinned", "H2D_hist_rebinned", nBinsX, binsX, nBinsY, binsY);
  
  // prior weighting
  double H2D_hist_content, priorWeightedSpectrumContent;
  double H2D_hist_contentError = 0;
  double priorWeightedSpectrumError = 0;
  for(int iBinXFine = 0; iBinXFine <= H2D_hist->GetNbinsX()+1; iBinXFine++){ // 0 and n+1 take underflow and overflow into account
    // if (debug == true) {cout << "iBinXFine = " << iBinXFine << endl;}
    priorWeightedSpectrumContent = H1D_priorSpectrum->GetBinContent(iBinXFine);
    priorWeightedSpectrumError = H1D_priorSpectrum->GetBinError(iBinXFine);
    for(int iBinYFine = 0; iBinYFine <= H2D_hist->GetNbinsY()+1; iBinYFine++){ // 0 and n+1 take underflow and overflow into account
      // if (debug == true) {cout << "     iBinYFine = " << iBinYFine << endl;}
      H2D_hist_content = H2D_hist->GetBinContent(iBinXFine, iBinYFine);
      H2D_hist_contentError = H2D_hist->GetBinError(iBinXFine, iBinYFine);
      H2D_hist->SetBinContent(iBinXFine, iBinYFine, H2D_hist_content * priorWeightedSpectrumContent);
      H2D_hist->SetBinError(iBinXFine, iBinYFine, sqrt( priorWeightedSpectrumContent*priorWeightedSpectrumContent * H2D_hist_contentError*H2D_hist_contentError + H2D_hist_content*H2D_hist_content * priorWeightedSpectrumError*priorWeightedSpectrumError));
    }
  }
  cout << "##############################################" << endl;
  cout << "######### reweighting with prior #############" << endl;
  cout << "##############################################" << endl;

  int ibinX_low, ibinX_high, ibinY_low, ibinY_high;
  double priorBinNorm, priorBinNormError;
  double hist2D_rebinned_content, hist2D_rebinned_error;
  for(int iBinY = 1; iBinY <= nBinsY; iBinY++){
  // for(int iBinY = 0; iBinY <= nBinsY+1; iBinY++){ // 0 and n+1 take underflow and overflow into account
    if (debug == true) {cout << "iBinY = " << iBinY << endl;}
    ibinY_low = H2D_hist->GetYaxis()->FindBin(H2D_hist_rebinned.GetYaxis()->GetBinLowEdge(iBinY) );
    ibinY_high = H2D_hist->GetYaxis()->FindBin(H2D_hist_rebinned.GetYaxis()->GetBinLowEdge(iBinY+1) ) - 1;
    for(int iBinX = 1; iBinX <= nBinsX; iBinX++){
    // for(int iBinX = 0; iBinX <= nBinsX+1; iBinX++){ // 0 and n+1 take underflow and overflow into account
      if (debug == true) {cout << "     iBinX = " << iBinX << endl;}
      ibinX_low = H2D_hist->GetXaxis()->FindBin(H2D_hist_rebinned.GetXaxis()->GetBinLowEdge(iBinX) ); // according to definition of getbinlowedge, it does work for under/overflows by giving the low edge of bin 1 minus the bin width (calculated assuming bins of equal sizes in all the axis); in fact, it works regardless of the value of iBin even for Nbin+2 or more (https://root.cern.ch/doc/master/classTAxis.html#a6d2d17ef00382842f58e88f356516a0d)
      ibinX_high = H2D_hist->GetXaxis()->FindBin(H2D_hist_rebinned.GetXaxis()->GetBinLowEdge(iBinX+1)) - 1;

      // dividing by prior integral over the xbin to make the weighting done before a probability
      priorBinNorm = H1D_priorSpectrum->IntegralAndError(ibinX_low, ibinX_high, priorBinNormError);
      for(int iBinYFine = ibinY_low; iBinYFine <= ibinY_high; iBinYFine++){
        for(int iBinXFine = ibinX_low; iBinXFine <= ibinX_high; iBinXFine++){
          priorBinNorm == 0 ? H2D_hist_content = 0 : H2D_hist_content = H2D_hist->GetBinContent(iBinXFine, iBinYFine)* 1./priorBinNorm;
          priorBinNorm == 0 ? H2D_hist_contentError = 0 : H2D_hist_contentError = sqrt( H2D_hist_contentError*H2D_hist_contentError / (H2D_hist_content*H2D_hist_content) + priorBinNormError*priorBinNormError / (priorBinNorm*priorBinNorm) );
          H2D_hist->SetBinContent(iBinXFine, iBinYFine, H2D_hist_content);
          H2D_hist->SetBinError(iBinXFine, iBinYFine, H2D_hist_contentError);
        }
      }
      hist2D_rebinned_content = H2D_hist->IntegralAndError(ibinX_low, ibinX_high, ibinY_low, ibinY_high, hist2D_rebinned_error);
      H2D_hist_rebinned.SetBinContent(iBinX, iBinY, hist2D_rebinned_content);
      H2D_hist_rebinned.SetBinError(iBinX, iBinY, hist2D_rebinned_error);
      if (debug == true) {cout << "ibinX_low = " << ibinX_low << ", ibinX_high = " << ibinX_high << ", ibinY_low = " << ibinY_low << ", ibinY_high = " << ibinY_high << "         --------          H2D_hist_rebinned(iBinX, iBinY) = " << H2D_hist_rebinned.GetBinContent(iBinX, iBinY) << endl;}
    }
  }
  cout << "##############################################" << endl;
  cout << "######## taking care of truncation ###########" << endl;
  cout << "##############################################" << endl;
  // Store gen truncations in x-underflow and overflow for each pT_gen (y-axis); needed for NormaliseYSlicesAsProbaDensity
  // I don't like the fact I couldn't add it to the above, but whatever, for now
  double hist2D_rebinned_underflowContent, hist2D_rebinned_overflowContent;
  double hist2D_rebinned_underflowError, hist2D_rebinned_overflowError;
  for(int iBinY = 1; iBinY <= nBinsY; iBinY++){
    ibinY_low = H2D_hist->GetYaxis()->FindBin(H2D_hist_rebinned.GetYaxis()->GetBinLowEdge(iBinY) );
    ibinY_high = H2D_hist->GetYaxis()->FindBin(H2D_hist_rebinned.GetYaxis()->GetBinLowEdge(iBinY+1) ) - 1;
    ibinX_low = H2D_hist->GetXaxis()->FindBin(H2D_hist_rebinned.GetXaxis()->GetBinLowEdge(1) ); // according to definition of getbinlowedge, it does work for under/overflows by giving the low edge of bin 1 minus the bin width (calculated assuming bins of equal sizes in all the axis); in fact, it works regardless of the value of iBin even for Nbin+2 or more (https://root.cern.ch/doc/master/classTAxis.html#a6d2d17ef00382842f58e88f356516a0d)
    ibinX_high = H2D_hist->GetXaxis()->FindBin(H2D_hist_rebinned.GetXaxis()->GetBinLowEdge(nBinsX+1) ); // according to definition of getbinlowedge, it does work for under/overflows by giving the low edge of bin 1 minus the bin width (calculated assuming bins of equal sizes in all the axis); in fact, it works regardless of the value of iBin even for Nbin+2 or more (https://root.cern.ch/doc/master/classTAxis.html#a6d2d17ef00382842f58e88f356516a0d)


    hist2D_rebinned_underflowContent = H2D_hist->IntegralAndError(0, ibinX_low, ibinY_low, ibinY_high, hist2D_rebinned_underflowError);
    priorBinNorm = H1D_priorSpectrum->IntegralAndError(0, ibinX_low, priorBinNormError);
    priorBinNorm == 0 ? hist2D_rebinned_underflowContent = 0 : hist2D_rebinned_underflowContent = hist2D_rebinned_underflowContent * 1./priorBinNorm;
    priorBinNorm == 0 ? hist2D_rebinned_underflowError = 0 : hist2D_rebinned_underflowError = sqrt( hist2D_rebinned_underflowError*hist2D_rebinned_underflowError / (hist2D_rebinned_underflowContent*hist2D_rebinned_underflowContent) + priorBinNormError*priorBinNormError / (priorBinNorm*priorBinNorm) );
    H2D_hist_rebinned.SetBinContent(0, iBinY, hist2D_rebinned_underflowContent);
    H2D_hist_rebinned.SetBinError(0, iBinY, hist2D_rebinned_underflowError);



    hist2D_rebinned_overflowContent = H2D_hist->IntegralAndError(ibinX_high, H2D_hist->GetNbinsX()+1, ibinY_low, ibinY_high, hist2D_rebinned_overflowError);
    priorBinNorm = H1D_priorSpectrum->IntegralAndError(ibinX_high, H2D_hist->GetNbinsX()+1, priorBinNormError);
    priorBinNorm == 0 ? hist2D_rebinned_overflowContent = 0 : hist2D_rebinned_overflowContent = hist2D_rebinned_overflowContent * 1./priorBinNorm;
    priorBinNorm == 0 ? hist2D_rebinned_overflowError = 0 : hist2D_rebinned_overflowError = sqrt( hist2D_rebinned_overflowError*hist2D_rebinned_overflowError / (hist2D_rebinned_overflowContent*hist2D_rebinned_overflowContent) + priorBinNormError*priorBinNormError / (priorBinNorm*priorBinNorm) );
    H2D_hist_rebinned.SetBinContent(nBinsX+1, iBinY, hist2D_rebinned_overflowContent);
    H2D_hist_rebinned.SetBinError(nBinsX+1, iBinY, hist2D_rebinned_overflowError);
    
  }
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
TH2D RebinVariableBins2D_ySlicePriorWeighted(TH2D* H2D_hist, int nBinsX, int nBinsY, double* binsX, double* binsY, TH1D* H1D_priorSpectrum){
  return RebinVariableBins2D_ySlicePriorWeighted(H2D_hist, nBinsX, nBinsY, binsX, binsY, H1D_priorSpectrum, false);
}

void NormaliseYSlicesAsProbaDensity(TH2D* H2D_hist){
  // Normalisation of the H2D_hist 2D histogram: each pt gen slice is normalised to unity, takes underflow into account if it's present
  double genSliceNorm = 1;
  // double dpT;
  double binContent, binError;
  for(int iBinY = 0; iBinY <= H2D_hist->GetNbinsY()+1; iBinY++){ // 0 and n+1 take underflow and overflow into account
    genSliceNorm = H2D_hist->Integral(0, H2D_hist->GetNbinsX()+1, iBinY, iBinY);
      cout << "iBinY = " << iBinY << "         --------  genSliceNorm = " << genSliceNorm << endl;
    // genSliceNorm = H2D_hist->Integral(1, H2D_hist->GetNbinsX(), iBinY, iBinY);
    // dpT = H2D_hist->GetYaxis()->GetBinWidth(iBinY);
    // cout << "iBinY = " << iBinY << "         --------          genSliceNorm = " << genSliceNorm << " ----  1/dpT*genSliceNorm = " << genSliceNorm*1./dpT << ", dpT = " << dpT << endl;
    for(int iBinX = 0; iBinX <= H2D_hist->GetNbinsX()+1; iBinX++){ // 0 and n+1 take underflow and overflow into account
      genSliceNorm == 0 ? binContent = 0 : binContent = H2D_hist->GetBinContent(iBinX, iBinY) *1./genSliceNorm;
      genSliceNorm == 0 ? binError = 0 : binError = H2D_hist->GetBinContent(iBinX, iBinY) *1./genSliceNorm;
      H2D_hist->SetBinContent(iBinX, iBinY, binContent);
      H2D_hist->SetBinError(iBinX, iBinY, binError);
      // cout << "iBinX = " << iBinX << ", iBinY = " << iBinY << "         --------  H2D_hist->GetBinContent(iBinX, iBinY) = " << H2D_hist->GetBinContent(iBinX, iBinY) << endl;
    }
  }
  for(int iBinY = 0; iBinY <= H2D_hist->GetNbinsY()+1; iBinY++){ // 0 and n+1 take underflow and overflow into account
    cout << "iBinY = " << iBinY << " -------- checking kinematic efficiency: sum of bins in slice iBinY = " << H2D_hist->Integral(1, H2D_hist->GetNbinsX(), iBinY, iBinY) << endl;
  } // I don't seem to have much overflow ...
  // IL ME MANQUE L OVERFLOW POUR UNE RAISON QUI M ECHAPPE

}

  // for(int ibinPt = 1; ibinPt <= nbinpT; ibinPt++){
  //   double dpT = hRawYield_vsPt->GetXaxis()->GetBinWidth(ibinPt);
  //   double drapidity = 1.5; //1.5
  //   double d2N_dpTdy = hRawYield_vsPt->GetBinContent(ibinPt) *1./dpT*1./drapidity;
  //   hRawYield_vsPt->SetBinContent(ibinPt,d2N_dpTdy);
  //   hRawYield_vsPt->SetBinError(ibinPt,hRawYield_vsPt->GetBinError(ibinPt) *1./dpT*1./drapidity);  // error on d2N_dpTdy 
  // }    


TH2D GetTransposeHistogram(TH2D* inputHist){
  int nBinsX = inputHist->GetNbinsX();
  int nBinsY = inputHist->GetNbinsY();
  std::vector<double> vectorBinsX = GetTH1Bins(inputHist->ProjectionX(inputHist->GetName()+(TString)"projX", 1, nBinsX));
  std::vector<double> vectorBinsY = GetTH1Bins(inputHist->ProjectionY(inputHist->GetName()+(TString)"projY", 1, nBinsY));
  double* binsX = &vectorBinsX[0];
  double* binsY = &vectorBinsY[0];

  TH2D transposedHist("inputHist_transposed", "inputHist_transposed", nBinsY, binsY, nBinsX, binsX);
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
  double* binsX_AB = &vectorBinsX[0];
  double* binsY_AB = &vectorBinsY[0];
  int nBinsX_AB = nBinsX_B;
  int nBinsY_AB = nBinsY_A;

  TH2D histAB("histAB", "histAB", nBinsX_B, binsX_AB, nBinsY_A, binsY_AB);

  // Roounfold sees the element of a hist histogram(i,j) as the element R(i,j) of a matrix; thus no transformation has to be done for the indices or axes of the histograms when doing the multiplication

  double productContent_iBinX_iBinY;
  double productError2_iBinX_iBinY;
  for(int iBinX = 1; iBinX <= nBinsX_AB; iBinX++){ // 0 and n+1 would take underflow and overflow into account, don't want that
    for(int iBinY = 1; iBinY <= nBinsY_AB; iBinY++){ // 0 and n+1 would take underflow and overflow into account, don't want that
      productContent_iBinX_iBinY = 0;
      productError2_iBinX_iBinY = 0;
      for(int iBinK = 1; iBinK <= nBinsK; iBinK++){ // 0 and n+1 would take underflow and overflow into account, don't want that
        productContent_iBinX_iBinY += histA->GetBinContent(iBinX, iBinK) * histB->GetBinContent(iBinK, iBinY); 
        productError2_iBinX_iBinY += pow(histA->GetBinError(iBinX, iBinK), 2) + pow(histB->GetBinError(iBinK, iBinY), 2); // simple sigma_ab = a2sigma_a2 + b2sigma_b2 ; that assumes there are no correlations; here it s background fluct from PbPB sim, and detector effects from a pp sim, so we can maybe say theyre not correlated    
      }
      histAB.SetBinContent(iBinX, iBinY, productContent_iBinX_iBinY);
      histAB.SetBinError(iBinX, iBinY, sqrt(productError2_iBinX_iBinY));  
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

  // Roounfold sees the element of a hist histogram(i,j) as the element R(i,j) of a matrix; thus no transformation has to be done for the indices or axes of the histograms when doing the multiplication

  double productContent_iBinX;
  double productError2_iBinX;
  for(int iBinX = 1; iBinX <= nBinsX_AU; iBinX++){ // 0 and n+1 would take underflow and overflow into account, don't want that
    productContent_iBinX = 0;
    productError2_iBinX = 0;
    for(int iBinK = 1; iBinK <= nBinsK; iBinK++){ // 0 and n+1 would take underflow and overflow into account, don't want that
      productContent_iBinX += histA->GetBinContent(iBinX, iBinK) * histU->GetBinContent(iBinK); 
      productError2_iBinX += pow(histA->GetBinError(iBinX, iBinK), 2) + pow(histU->GetBinError(iBinK), 2); // simple sigma_ab = a2sigma_a2 + b2sigma_b2 ; that assumes there are no correlations; here it s background fluct from PbPB sim, and detector effects from a pp sim, so we can maybe say theyre not correlated    
    }
    histAU.SetBinContent(iBinX, productContent_iBinX);
    histAU.SetBinError(iBinX, sqrt(productError2_iBinX));  
  }
  return histAU;
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