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


M3file_photon = 'templates_barrel_scaled.root'
M3file_presel_scaled = 'templates_presel_scaled.root'


def integral_bins(hist,bin1,bin2):
	err = ROOT.Double(0.0)
	integr = hist.IntegralAndError(bin1, bin2, err)
	print integr,err
	return integr,err


def doM3fit_photon():
	print 'M3 fit after photon selection'
	varToFit = 'M3'
	DataHist = get1DHist(M3file_photon, 'Data_'+varToFit)

	lowfitBin = DataHist.FindBin(50.01)
	highfitBin = DataHist.FindBin(499.99)
	
	TopHist = get1DHist(M3file_photon, 'TTJets_'+varToFit)
	TopHist.Add(get1DHist(M3file_photon, 'TTGamma_'+varToFit))
	
	# remove extra contributions from data
	WJetsHist = get1DHist(M3file_presel_scaled, 'WJets_'+varToFit)
	WJetsHist.Scale(get1DHist(M3file_photon, 'WJets_'+varToFit).Integral() / WJetsHist.Integral())
	DataHist.Add(WJetsHist, -1.0)
	
	ZJetsHist = get1DHist(M3file_presel_scaled, 'ZJets_'+varToFit)
	ZJetsHist.Scale(get1DHist(M3file_photon, 'ZJets_'+varToFit).Integral() / ZJetsHist.Integral())
	DataHist.Add(ZJetsHist, -1.0)
	
	SingleTopHist = get1DHist(M3file_presel_scaled, 'SingleTop_'+varToFit)
	SingleTopHist.Scale(get1DHist(M3file_photon, 'SingleTop_'+varToFit).Integral() / SingleTopHist.Integral())
	DataHist.Add(SingleTopHist, -1.0)
	
	
	#DataHist.Add(get1DHist(M3file_photon, 'Diboson_'+varToFit), -1.0)
	QCDHist = get1DHist(M3file_presel_scaled, 'QCD_'+varToFit)
	QCDHist.Scale(get1DHist(M3file_photon, 'QCD_'+varToFit).Integral() / QCDHist.Integral())
	DataHist.Add(QCDHist, -1.0)
	
	
	#BGHist = get1DHist(M3file_photon, 'WJets_'+varToFit)
	#BGHist.Add(get1DHist(M3file_photon, 'ZJets_'+varToFit))
	#BGHist.Add(get1DHist(M3file_photon, 'SingleTop_'+varToFit))
	#BGHist.Add(get1DHist(M3file_photon, 'Diboson_'+varToFit))
	#BGHist.Add(get1DHist(M3file_photon, 'QCD_'+varToFit))
	
	BGHist = get1DHist(M3file_presel_scaled, 'Vgamma_'+varToFit)
	BGphotonSel = get1DHist(M3file_photon, 'Vgamma_'+varToFit)

	dataInt,dataIntErr = integral_bins(DataHist,lowfitBin,highfitBin)
	topInt,topIntErr = integral_bins(TopHist,lowfitBin,highfitBin)
	bgInt,bgIntErr = integral_bins(BGphotonSel,lowfitBin,highfitBin)
	
	#DataHist.Rebin(2)
	#TopHist.Rebin(2)
	#BGHist.Rebin(2)
	
	(m3TopFrac, m3TopFracErr) = makeFit(varToFit+'(GeV), photon selection', 50.0, 500.0, TopHist, BGHist, DataHist, varToFit+'_photon_fit.png')
				
	
	bgSF = (1.0 - m3TopFrac) * dataInt / bgInt
	bgSFerror = bgSF * ( (m3TopFracErr / (1.0 - m3TopFrac))**2 + (bgIntErr/bgInt)**2 + (dataIntErr/dataInt)**2 )**0.5
	
	topSF = m3TopFrac * dataInt / topInt
	topSFErr = topSF * ( (m3TopFracErr / m3TopFrac)**2 + (topIntErr/topInt)**2 + (dataIntErr/dataInt)**2 )**0.5
	
	print '#'*80
	print 'Correction to Vgamma samples after M3 fit: ', bgSF, ' +-', bgSFerror, '(fit + stat error)'
	print 'Correction to ttbar samples after M3 fit: ', topSF, ' +-', topSFErr, '(fit + stat error)'
	print '#'*80
	return (bgSF,bgSFerror)

#doM3fit_photon()

