// Commonly used titles

////////////////////////////////////////////
////////////////// x-axis //////////////////
////////////////////////////////////////////

TString* texRapidity = new TString("y");
TString* texPtX = new TString("#it{p}_{T} (GeV/#it{c})");
TString* texEtaX = new TString("#it{#eta}");
TString* texPhiX = new TString("#it{#phi} (rad)");
TString* texNTracksX = new TString("#it{N}_{tracks}");
TString* texLeadPt = new TString("#it{p}_{T,track}^{leading}");
TString* texJetArea = new TString("#it{A}_{jet}");
TString* texJetNTracks = new TString("#it{N}_{tracks}");
TString* texJetNTrdTracks = new TString("#it{N}_{TRD tracks}");
TString* texJetNTrdTracksRatio = new TString("#frac{#it{N}_{TRD tracks}}{#it{N}_{tracks}}");


////////////////////////////////////////////
////////////////// y-axis //////////////////
////////////////////////////////////////////

TString* texPtDifferentialYield = new TString("1/#it{N}_{ev} d#it{N}/d#it{p}_{T} (GeV/#it{c})^{-1}");
TString* texCount = new TString("count");

TString* texPtMeasured = new TString("#it{p}_{T}^{measured} (GeV/#it{c})");
TString* texPtMC = new TString("#it{p}_{T}^{MC} (GeV/#it{c})");

TString* texRatioDatasets = new TString(*texDatasetsComparisonType+"/"+*texDatasetsComparisonType+" "+DatasetsNames[0]); 

////////// Jets
TString* texNSubJettinessRatio = new TString("#it{#tau_{2}/#tau_{1}}");
TString* texDeltaR = new TString("#it{#Delta R}");
TString* texdN_dsubratio = new TString("1/N^{ev} d#it{N}/d(#it{#tau_{2}/#tau_{1}})");
TString* texdN_dDeltaR = new TString("1/N^{ev} d#it{N}/d(#it{#Delta R})");

TString* texJetNormPtYield = new TString("1/#it{N}_{ev} d#it{N}_{jet}/d#it{p}_{T} (GeV/#it{c})^{-1}");
TString* texJetNormEtaYield = new TString("1/#it{N}_{ev} d#it{N}_{jet}/d#it{#eta}");
TString* texJetNormPhiYield = new TString("1/#it{N}_{ev} d#it{N}_{jet}/d#it{#phi}");
TString* texJetNormNTracksYield = new TString("1/#it{N}_{ev} d#it{N}_{jet}/d#it{N}_{tracks}");
TString* texJetNormAreaYield = new TString("1/#it{N}_{ev} d#it{N}_{jet}/d#it{A}_{jet}");

TString* texJetRatioEtaComparison = new TString("#it{N}_{jet}^{ #eta > 0} / #it{N}_{jet}^{ #eta < 0}");
TString* texJetRatioTRDvsNoTRD = new TString("(#it{N}_{jet}^{ #eta > 0} / #it{N}_{jet}^{ #eta < 0})_{onlyTRDtracks} / (#it{N}_{jet}^{ #eta > 0} / #it{N}_{jet}^{ #eta < 0})_{noTRDtracks}");


////////// Tracks
TString* texTrackNormPtYield = new TString("1/#it{N}_{ev} d#it{N}_{track}/d#it{p}_{T} (GeV/#it{c})^{-1}");
TString* texTrackNormEtaYield = new TString("1/#it{N}_{ev} d#it{N}_{track}/d#it{#eta}");
TString* texTrackNormPhiYield = new TString("1/#it{N}_{ev} d#it{N}_{track}/d#it{#phi}");
