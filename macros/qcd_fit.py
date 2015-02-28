import ROOT

def makeFit(varname, varmin, varmax, signalHist, backgroundHist, dataHist, plotName):
	# RooFit variables
	sihihVar = ROOT.RooRealVar(varname, varname, varmin, varmax)
	sihihArgList = ROOT.RooArgList()
	sihihArgList.add(sihihVar)
	sihihArgSet = ROOT.RooArgSet()
	sihihArgSet.add(sihihVar)

	# create PDFs
	signalDataHist = ROOT.RooDataHist('signalDataHist','signal RooDataHist', sihihArgList, signalHist)
	signalPdf = ROOT.RooHistPdf('signalPdf',varname+' of signal', sihihArgSet, signalDataHist)

	backgroundDataHist = ROOT.RooDataHist('backgroundDataHist','background RooDataHist', sihihArgList, backgroundHist)
	backgroundPdf = ROOT.RooHistPdf('backgroundPdf',varname+' of background', sihihArgSet, backgroundDataHist)

	# data
	dataDataHist = ROOT.RooDataHist('data '+varname, varname+' in Data', sihihArgList, dataHist)

	# signal fraction parameter
	sfname = 'signal fraction'
	if 'MET' in varname:
		sfname = 'multijet fraction'
	if 'M3' in varname:
		sfname = 'ttbar fraction'
	signalFractionVar = ROOT.RooRealVar(sfname,sfname, 0.5, 0.0, 1.0)
	sumPdf = ROOT.RooAddPdf('totalPdf','signal and background', signalPdf, backgroundPdf, signalFractionVar)
	
	# fit
	sumPdf.fitTo( dataDataHist, ROOT.RooFit.SumW2Error(ROOT.kFALSE), ROOT.RooFit.PrintLevel(-1) )
	
	if plotName!='':
		# plot results
		c1 = ROOT.TCanvas('c1', 'c1', 800, 600)
		plotter = ROOT.RooPlot('myplot','',sihihVar,varmin,varmax,20) # nBins is dummy
		dataDataHist.plotOn(plotter, ROOT.RooFit.Name('data'))
		sumPdf.plotOn(plotter, ROOT.RooFit.Name('sum'), ROOT.RooFit.LineColor(ROOT.kRed))
		sumPdf.plotOn(plotter, ROOT.RooFit.Components('signalPdf'), ROOT.RooFit.Name('signal'), 
			ROOT.RooFit.LineColor(ROOT.kGreen))
		sumPdf.plotOn(plotter, ROOT.RooFit.Components('backgroundPdf'), ROOT.RooFit.Name('background'), 
			ROOT.RooFit.LineColor(ROOT.kBlue))
		sumPdf.paramOn(plotter) # fix

		plotter.Draw()
		plotter.GetYaxis().SetTitleOffset(1.4)
		c1.SaveAs(plotName)
	print 'fit returned value ',signalFractionVar.getVal(),' +- ',signalFractionVar.getError()
	return (signalFractionVar.getVal(),signalFractionVar.getError())

openfiles = {}

def get1DHist(filename, histname):
	if filename not in openfiles:
		openfiles[filename] = ROOT.TFile(filename,'READ')
	file = openfiles[filename]
		
	hist = file.Get(histname)
	hist.SetDirectory(0)
	hist.SetFillColor(0)
	#hist.Sumw2()
	return hist

qcdMETfile = 'templates_presel_nomet_qcd.root'
normMETfile = 'templates_presel_nomet.root'

M3file = 'templates_presel.root'

M3file_photon = 'templates_barrel_scaled.root'
M3file_presel_scaled = 'templates_presel_scaled.root'

def doM3fit():
	print 'now do M3 fit'

	varToFit = 'M3'

	DataHist = get1DHist(M3file, 'Data_'+varToFit)

	TopHist = get1DHist(M3file, 'TTJets_'+varToFit)
	TopHist.Add(get1DHist(M3file, 'TTGamma_'+varToFit))

	WJHist = get1DHist(M3file, 'WJets_'+varToFit)
	#WJHist.Add(get1DHist(M3file, 'ZJets_'+varToFit)) #     fix
	#WJHist.Add(get1DHist(M3file, 'SingleTop_'+varToFit)) #     fix
	#WJHist.Add(get1DHist(M3file, 'QCD_'+varToFit)) #     fix
	#WJHist.Add(get1DHist(M3file, 'Vgamma_'+varToFit)) #     fix
	
	# remove other suspects from data
	DataHist.Add(get1DHist(M3file, 'ZJets_'+varToFit), -1.0)
	DataHist.Add(get1DHist(M3file, 'SingleTop_'+varToFit), -1.0)
	DataHist.Add(get1DHist(M3file, 'Vgamma_'+varToFit), -1.0)

	(m3TopFrac, m3TopFracErr) = makeFit(varToFit+'(GeV)', 40.0, 660.0, TopHist, WJHist, DataHist, varToFit+'_fit.png')
	lowfitBin = DataHist.FindBin(40.01)
	highfitBin = DataHist.FindBin(659.99)
			
	dataInt = DataHist.Integral(lowfitBin,highfitBin)
	topInt = TopHist.Integral(lowfitBin,highfitBin)
	WJInt = WJHist.Integral(lowfitBin,highfitBin)
	
	TopSF = m3TopFrac * dataInt / topInt
	TopSFerror = m3TopFracErr * dataInt / topInt
	print '#'*80
	print 'Correction to the Top scale factor: ', TopSF, ' +-', TopSFerror, '(fit error only)'
	WJetsSF = (1.0-m3TopFrac) * dataInt / WJInt
	WJetsSFerror = m3TopFracErr * dataInt / WJInt
	print 'Correction to WJets scale factor: ', WJetsSF, ' +-',WJetsSFerror,'(fit error only)'
	print '#'*80
	return (TopSF, TopSFerror, WJetsSF, WJetsSFerror)

