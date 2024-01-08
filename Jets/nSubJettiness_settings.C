
// Options to be set:
// TFile* file_O2Analysis = new TFile("AnalysisResults_LHC23d4_fullProdHyperloop.root");
TFile* file_O2Analysis = new TFile("AnalysisResults.root");
// TString* texCollisionDataInfo = new TString("#splitline{pp #sqrt{#it{s}} = 13.6 TeV}{Jet-Jet simulation LHC23d4}");
TString* texCollisionDataInfo = new TString("pp #sqrt{#it{s}} = 13.6 TeV - Jet-Jet sim LHC23d4");

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
