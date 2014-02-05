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


varToFit = 'MET'
qcdDataHist = get1DHist('templates_presel_nomet_qcd.root', 'Data_'+varToFit)
# remove MC contribution
qcdDataHist.Add(get1DHist('templates_presel_nomet_qcd.root', 'TTJets_'+varToFit), -1)
qcdDataHist.Add(get1DHist('templates_presel_nomet_qcd.root', 'WJets_'+varToFit), -1)


DataHist = get1DHist('templates_presel_nomet.root', 'Data_'+varToFit)

MCHist = get1DHist('templates_presel_nomet.root', 'TTJets_'+varToFit)
MCHist.Add(get1DHist('templates_presel_nomet.root', 'WHIZARD_'+varToFit))
MCHist.Add(get1DHist('templates_presel_nomet.root', 'WJets_'+varToFit))
MCHist.Add(get1DHist('templates_presel_nomet.root', 'ZJets_'+varToFit))
MCHist.Add(get1DHist('templates_presel_nomet.root', 'SingleTop_'+varToFit))

(metFrac, metFracErr) = makeFit(varToFit, 0.0, 300.0, qcdDataHist, MCHist, DataHist, varToFit+'_QCD_fit.png')

# recalculate data-driven QCD normalization
lowbin = 1   #DataHist.FindBin(20.01)
highbin = DataHist.GetNbinsX()+1 # overflow bin included

print 'Will calculate integral in the bin range:', lowbin, highbin
dataInt = DataHist.Integral(lowbin, highbin)
print 'Integral of Data in the desired range: ', dataInt
qcdInt = qcdDataHist.Integral(lowbin, highbin)
print 'Integral of data-driven QCD in the desired range: ', qcdInt
print '#'*80
# take into account only fit error
# stat errors on histograms are treated while calculating the final answer
print 'Scale factor for QCD in nominal MET range: ', metFrac*dataInt/qcdInt,' +-',metFracErr*dataInt/qcdInt,'(fit error only)'
print '#'*80

print 'now do M3 fit'

varToFit = 'M3'

DataHist = get1DHist('templates_presel.root', 'Data_'+varToFit)

TopHist = get1DHist('templates_presel.root', 'TTJets_'+varToFit)
TopHist.Add(get1DHist('templates_presel.root', 'WHIZARD_'+varToFit))

WJHist = get1DHist('templates_presel.root', 'WJets_'+varToFit)

# remove other suspects from data
DataHist.Add(get1DHist('templates_presel.root', 'ZJets_'+varToFit), -1.0)
DataHist.Add(get1DHist('templates_presel.root', 'SingleTop_'+varToFit), -1.0)
DataHist.Add(get1DHist('templates_presel.root', 'QCD_'+varToFit), -1.0)

(m3TopFrac, m3TopFracErr) = makeFit(varToFit, 70.0, 500.0, TopHist, WJHist, DataHist, varToFit+'_fit.png')
dataInt = DataHist.Integral()
topInt = TopHist.Integral()
WJInt = WJHist.Integral()

print 'Correction to the Top scale factor: ', m3TopFrac * dataInt / topInt, ' +-',m3TopFracErr * dataInt / topInt,'(fit error only)'
print 'Correction to WJets scale factor: ', (1.0-m3TopFrac) * dataInt / WJInt, ' +-',m3TopFracErr * dataInt / WJInt,'(fit error only)'

