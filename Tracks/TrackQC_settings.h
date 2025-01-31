
// Bin edge control
// float GLOBAL_epsilon = 0.00001;

// Analysis settings
const int nJetType = 3;
const TString jetType[nJetType] = {"charged", "neutral", "full"};
const int nJetLevel = 3;
const TString jetLevel[nJetLevel] = {"data", "mcd", "mcp"};
const int nRadius = 3;
const TString RadiusLegend[nRadius] = {"R = 0.2", "R = 0.4", "R = 0.6"};
float arrayRadius[nRadius] = {0.2, 0.4, 0.6};
// const int nRadius = 1;
// const TString RadiusLegend[nRadius] = {"R = 0.4"};
// float arrayRadius[nRadius] = {0.4};
// const int nRadius = 9;
// const TString RadiusLegend[nRadius] = {"R = 0.2", "R = 0.25", "R = 0.3", "R = 0.35", "R = 0.4", "R = 0.45", "R = 0.5", "R = 0.55", "R = 0.6"};
// float arrayRadius[nRadius] = {0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6};
// Choice of jet type (charged, neutral, full) and level (data, detector level, particle level)
const int iJetType = 0;
const int iJetLevel = 0;


const int nCentralityBins = 7;
const float arrayCentralityBinning[nCentralityBins+1] = {0, 10, 20, 30, 40, 60, 80, 100};
// const int nCentralityBins = 3;
// const float arrayCentralityBinning[nCentralityBins+1] = {0, 10, 50, 90};