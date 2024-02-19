// Commonly used titles

////////////////////////////////////////////
////////////////// x-axis //////////////////
////////////////////////////////////////////

TString* texRapidity = new TString("y");
TString* texPtX = new TString("#it{p}_{T} (GeV/#it{c})");
TString* texEtaX = new TString("#it{#eta}");
TString* texPhiX = new TString("#it{#phi} (rad)");
TString* texNTracksX = new TString("#it{N}_{tracks in jet}");
TString* texLeadPt = new TString("#it{p}_{T,track}^{leading}");
TString* texJetArea = new TString("#it{A}_{jet}");
TString* texJetNTracks = new TString("#it{N}_{tracks}");
TString* texJetNTrdTracks = new TString("#it{N}_{TRD tracks}");
TString* texJetNTrdTracksRatio = new TString("#frac{#it{N}_{TRD tracks}}{#it{N}_{tracks}}");
TString* texCentrality = new TString("Centrality (%)");
TString* texSelectedMultiplicity = new TString("#it{N}_{tracks}^{selected}");
TString* texLeadJetPt = new TString("#it{p}_{T,lead jet} (GeV/#it{c})");

TString* texPtJetRecX = new TString("#it{p}_{T,jet rec} (GeV/#it{c})");
TString* texPtJetGenX = new TString("#it{p}_{T,jet gen} (GeV/#it{c})");



////////////////////////////////////////////
////////////////// y-axis //////////////////
////////////////////////////////////////////

TString* texPtDifferentialYield = new TString("1/#it{N}_{ev} d#it{N}/d#it{p}_{T} (GeV/#it{c})^{-1}");
TString* texCount = new TString("count");

TString* texPtMeasured = new TString("#it{p}_{T}^{measured} (GeV/#it{c})");
TString* texPtMC = new TString("#it{p}_{T}^{MC} (GeV/#it{c})");

TString* texRatioDatasets = new TString(*texDatasetsComparisonType+"/"+*texDatasetsComparisonType+" "+DatasetsNames[0]); 

TString* texRho = new TString("#it{#rho} (GeV/#it{c})");
TString* texBkgFluctuationRandomCone = new TString("#delta #it{p}_{T} = #it{p}_{T,random cone} - #it{A}_{random cone} * #it{#rho} (GeV/c)");
TString* texEntriesNorm_BkgFluctuationYield = new TString("1/#it{N}_{coll}^{in cent} d#it{N}_{coll}/d#it{#delta pt}");
TString* texCollNorm_BkgFluctuationYield = new TString("1/#it{N}_{coll} d#it{N}_{coll}/d#it{#delta pt}");
TString* texCollNorm_BkgFluctuationYield_CentWindow = new TString("1/#it{N}_{coll}^{in cent} d#it{N}_{coll}/d#it{#delta pt}");
TString* texEntriesNormRho = new TString("1/#it{N}_{coll}^{in cent} d#it{N}_{coll}/d#it{#rho}");
TString* texCollNorm_RhoYield = new TString("1/#it{N}_{coll} d#it{N}_{coll}/d#it{#rho}");
TString* texCollNorm_RhoYield_CentWindow = new TString("1/#it{N}_{coll}^{in cent} d#it{N}_{coll}/d#it{#rho}");
TString* texCollNorm_selMultYield = new TString("1/#it{N}_{coll} d#it{N}_{coll}/d#it{N}_{tracks}^{selected}");
TString* texCollNorm_selMultYield_CentWindow = new TString("1/#it{N}_{coll}^{in cent} d#it{N}_{coll}/d#it{N}_{tracks}^{selected}");


TString* texSigmaBkgFluctFit = new TString("#sigma_{#it{#delta pt}}^{fit}");
TString* texMeanBkgFluctFit = new TString("<#it{#delta pt}>^{fit}");

TString* texRhoMean = new TString("<#it{#rho}>");


//////////////////////////
////////// Jets //////////
TString* texNSubJettinessRatio = new TString("#it{#tau_{2}/#tau_{1}}");
TString* texDeltaR = new TString("#it{#Delta R}");
TString* texdN_dsubratio = new TString("1/N^{ev} d#it{N}/d(#it{#tau_{2}/#tau_{1}})");
TString* texdN_dDeltaR = new TString("1/N^{ev} d#it{N}/d(#it{#Delta R})");

TString* texJetPtYield_EventNorm = new TString("1/#it{N}_{ev} d#it{N}_{jet}/d#it{p}_{T} (GeV/#it{c})^{-1}");
TString* texJetEtaYield_EventNorm = new TString("1/#it{N}_{ev} d#it{N}_{jet}/d#it{#eta}");
TString* texJetPhiYield_EventNorm = new TString("1/#it{N}_{ev} d#it{N}_{jet}/d#it{#phi}");
TString* texJetPtYield_EntriesNorm = new TString("1/#it{N}_{jet} d#it{N}_{jet}/d#it{p}_{T} (GeV/#it{c})^{-1}");
TString* texJetEtaYield_EntriesNorm = new TString("1/#it{N}_{jet} d#it{N}_{jet}/d#it{#eta}");
TString* texJetPhiYield_EntriesNorm = new TString("1/#it{N}_{jet} d#it{N}_{jet}/d#it{#phi}");
TString* texJetPtYield_EventNorm_CentWindow = new TString("1/#it{N}_{ev}^{in cent} d#it{N}_{jet}/d#it{p}_{T} (GeV/#it{c})^{-1}");
TString* texJetEtaYield_EventNorm_CentWindow = new TString("1/#it{N}_{ev}^{in cent} d#it{N}_{jet}/d#it{#eta}");
TString* texJetPhiYield_EventNorm_CentWindow = new TString("1/#it{N}_{ev}^{in cent} d#it{N}_{jet}/d#it{#phi}");
TString* texJetNormNTracksYield = new TString("1/#it{N}_{ev} d#it{N}_{jet}/d#it{N}_{tracks}");
TString* texJetNormAreaYield = new TString("1/#it{N}_{ev} d#it{N}_{jet}/d#it{A}_{jet}");


TString* texJetRatioEtaComparison = new TString("#it{N}_{jet}^{ #eta > 0} / #it{N}_{jet}^{ #eta < 0}");
TString* texJetRatioTRDvsNoTRD = new TString("(#it{N}_{jet}^{ #eta > 0} / #it{N}_{jet}^{ #eta < 0})_{onlyTRDtracks} / (#it{N}_{jet}^{ #eta > 0} / #it{N}_{jet}^{ #eta < 0})_{noTRDtracks}");

TString* texNGhost = new TString("#it{N}_{ghost}");
TString* texBinSpacing = new TString("Area steps");
TString* texJetRatioAreaJetVsNGhost = new TString("#it{A}_{jet}^{O2}/ #it{A}_{jet}^{0.005*nGhost}");

////////////////////////////
////////// Tracks //////////
TString* texTrackPtYield_EventNorm = new TString("1/#it{N}_{ev} d#it{N}_{track}/d#it{p}_{T} (GeV/#it{c})^{-1}");
TString* texTrackEtaYield_EventNorm = new TString("1/#it{N}_{ev} d#it{N}_{track}/d#it{#eta}");
TString* texTrackPhiYield_EventNorm = new TString("1/#it{N}_{ev} d#it{N}_{track}/d#it{#phi}");
TString* texTrackPtYield_EntriesNorm = new TString("1/#it{N}_{track} d#it{N}_{track}/d#it{p}_{T} (GeV/#it{c})^{-1}");
TString* texTrackEtaYield_EntriesNorm = new TString("1/#it{N}_{track} d#it{N}_{track}/d#it{#eta}");
TString* texTrackPhiYield_EntriesNorm = new TString("1/#it{N}_{track} d#it{N}_{track}/d#it{#phi}");
