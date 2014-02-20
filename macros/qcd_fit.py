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
	signalFractionVar = ROOT.RooRealVar('signal fraction','signal fraction', 0.5, 0.0, 1.0)
	sumPdf = ROOT.RooAddPdf('totalPdf','signal and background', signalPdf, backgroundPdf, signalFractionVar)
	
	# fit
	sumPdf.fitTo( dataDataHist, ROOT.RooFit.SumW2Error(ROOT.kFALSE), ROOT.RooFit.PrintLevel(-1) )
	
	if plotName!='':
		# plot results
		c1 = ROOT.TCanvas('c1', 'c1', 800, 600)
		plotter = ROOT.RooPlot(sihihVar,varmin,varmax,20) # nBins is dummy 
		dataDataHist.plotOn(plotter, ROOT.RooFit.Name('data'))
		sumPdf.plotOn(plotter, ROOT.RooFit.Name('sum'), ROOT.RooFit.LineColor(ROOT.kRed))
		sumPdf.plotOn(plotter, ROOT.RooFit.Components('signalPdf'), ROOT.RooFit.Name('signal'), 
			ROOT.RooFit.LineColor(ROOT.kGreen))
		sumPdf.plotOn(plotter, ROOT.RooFit.Components('backgroundPdf'), ROOT.RooFit.Name('background'), 
			ROOT.RooFit.LineColor(ROOT.kBlue))
		sumPdf.paramOn(plotter)

		plotter.Draw()
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

def doQCDfit():
	varToFit = 'MET'

	qcdDataHist = get1DHist(qcdMETfile, 'Data_'+varToFit)
	# remove MC contribution
	qcdDataHist.Add(get1DHist(qcdMETfile, 'TTJets_'+varToFit), -1)
	#qcdDataHist.Add(get1DHist(qcdMETfile, 'WJets_'+varToFit), -1)


	DataHist = get1DHist(normMETfile, 'Data_'+varToFit)

	MCHist = get1DHist(normMETfile, 'TTJets_'+varToFit)
	MCHist.Add(get1DHist(normMETfile, 'WHIZARD_'+varToFit))
	MCHist.Add(get1DHist(normMETfile, 'WJets_'+varToFit))
	MCHist.Add(get1DHist(normMETfile, 'ZJets_'+varToFit))
	MCHist.Add(get1DHist(normMETfile, 'SingleTop_'+varToFit))

	(metFrac, metFracErr) = makeFit(varToFit, 0.0, 300.0, qcdDataHist, MCHist, DataHist, varToFit+'_QCD_fit.png')

	# recalculate data-driven QCD normalization
	lowbin = 1   #DataHist.FindBin(20.01)
	highbin = DataHist.GetNbinsX()+1 # overflow bin included

	print 'Will calculate integral in the bin range:', lowbin, highbin
	dataInt = DataHist.Integral(lowbin, highbin)
	print 'Integral of Data in the desired range: ', dataInt
	qcdInt = qcdDataHist.Integral(lowbin, highbin)
	mcInt = MCHist.Integral(lowbin, highbin)
	print 'Integral of data-driven QCD in the desired range: ', qcdInt
	print '#'*80
	# take into account only fit error
	# stat errors on histograms are treated while calculating the final answer
	QCDSF = metFrac*dataInt/qcdInt
	QCDSFerror = metFracErr*dataInt/qcdInt
	print 'Scale factor for QCD in nominal MET range: ', QCDSF,' +-',QCDSFerror,'(fit error only)'
	print 'Correction to all MC scale factors: ', (1-metFrac)*dataInt/mcInt, ' +-',metFracErr*dataInt/mcInt,'(fit error only)'
	print '#'*80
	return (QCDSF, QCDSFerror)


def doM3fit():
	print 'now do M3 fit'

	varToFit = 'M3'

	DataHist = get1DHist(M3file, 'Data_'+varToFit)

	TopHist = get1DHist(M3file, 'TTJets_'+varToFit)
	TopHist.Add(get1DHist(M3file, 'WHIZARD_'+varToFit))

	WJHist = get1DHist(M3file, 'WJets_'+varToFit)

	# remove other suspects from data
	DataHist.Add(get1DHist(M3file, 'ZJets_'+varToFit), -1.0)
	DataHist.Add(get1DHist(M3file, 'SingleTop_'+varToFit), -1.0)
	DataHist.Add(get1DHist(M3file, 'QCD_'+varToFit), -1.0)
	DataHist.Add(get1DHist(M3file, 'Vgamma_'+varToFit), -1.0)

	(m3TopFrac, m3TopFracErr) = makeFit(varToFit, 70.0, 500.0, TopHist, WJHist, DataHist, varToFit+'_fit.png')
	dataInt = DataHist.Integral()
	topInt = TopHist.Integral()
	WJInt = WJHist.Integral()
	TopSF = m3TopFrac * dataInt / topInt
	TopSFerror = m3TopFracErr * dataInt / topInt
	print '#'*80
	print 'Correction to the Top scale factor: ', TopSF, ' +-', TopSFerror, '(fit error only)'
	WJetsSF = (1.0-m3TopFrac) * dataInt / WJInt
	WJetsSFerror = m3TopFracErr * dataInt / WJInt
	print 'Correction to WJets scale factor: ', WJetsSF, ' +-',WJetsSFerror,'(fit error only)'
	print '#'*80
	return (TopSF, TopSFerror, WJetsSF, WJetsSFerror)


