import ROOT
import array

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
	

def get1DHist(filename, histname):
	file = ROOT.TFile(filename,'READ')
	hist = file.Get(histname)
	hist.SetDirectory(0)
	file.Close()
	hist.Sumw2()
	return hist

def getTemplFrom2Dhist(filename, histname, auxAxis, minval, maxval):
	file = ROOT.TFile(filename,'READ')
	hist2d = file.Get(histname)
	hist2d.SetDirectory(0)
	file.Close()
	
	if auxAxis == 'Y':
		firstBin = hist2d.GetYaxis().FindBin(minval)
		lastBin = hist2d.GetYaxis().FindBin(maxval) - 1
	elif auxAxis == 'X':
		firstBin = hist2d.GetXaxis().FindBin(minval)
		lastBin = hist2d.GetXaxis().FindBin(maxval) - 1
	else:
		print 'ERROR: valid values for auxAxis are X or Y'
		return None

	#print 'Getting histogram ',filename,' : ',histname
	#print 'Projection with auxAxis ', auxAxis, '  range ', minval, maxval, '  bins ', firstBin, lastBin
	
	if auxAxis == 'Y':
		templhist = hist2d.ProjectionX('projX'+filename+histname+'_'+str(minval)+'_'+str(maxval), firstBin,lastBin)
	else:
		templhist = hist2d.ProjectionY('projY'+filename+histname+'_'+str(minval)+'_'+str(maxval), firstBin,lastBin)
	templhist.SetDirectory(0)
	templhist.Sumw2()
	return templhist

def optimizeBinBoundaries(hist, minAcc, firstBinValue=-9999, lastBinValue=9999):
	nbins = hist.GetNbinsX()
	accu = 0
	binlist = []
	
	# ignore empty bins in the beginning
	for firstind in xrange(1,nbins+1):
		if hist.GetBinContent(firstind) > 0.0001 or hist.GetBinLowEdge(firstind)>=firstBinValue:
			binlist.append(hist.GetBinLowEdge(firstind))
			break
	
	for ind in xrange(firstind,nbins+1):
		if hist.GetBinLowEdge(ind)>=lastBinValue:
			binlist.append(lastBinValue)
			break
		if accu > minAcc:
			accu=0
			binlist.append(hist.GetBinLowEdge(ind))
		accu+=hist.GetBinContent(ind)
		#print 'bin ', ind, '  value ', hist.GetBinContent(ind), '  low and high edges ', hist.GetBinLowEdge(ind), hist.GetBinLowEdge(ind+1)
	binlist.append(hist.GetBinLowEdge(nbins+1))
	print binlist
	newarray = array.array('d')
	newarray.fromlist(binlist)
	return newarray
	
def makeUniformHist(hist):
	nbins = hist.GetNbinsX()
	#print 'nbins',nbins
	outhist = ROOT.TH1F(hist.GetName()+'_uni',hist.GetName()+'_uni',nbins,1,nbins+1)
	for ind in xrange(1,nbins+1):
		outhist.SetBinContent(ind,hist.GetBinContent(ind))
		outhist.SetBinError(ind,hist.GetBinError(ind))
	return outhist


random = ROOT.TRandom3()

def makeRandUniformHist(hist):
	nbins = hist.GetNbinsX()
	#print 'nbins',nbins
	outhist = ROOT.TH1F(hist.GetName()+'_uni',hist.GetName()+'_uni',nbins,1,nbins)
	for ind in xrange(1,nbins+1):
		outhist.SetBinContent(ind, random.Poisson(hist.GetBinContent(ind)))
	return outhist


# scale factors
WHIZARD_SF = 19.63*1000 * 0.9081*2/1074860
TTJets1l_SF = 19.63*1000 * 99.44 / 24849108
TTJets2l_SF = 19.63*1000 * 23.83 / 12086717
#TTJetsInc_SF = 19.63*1000 * 227/6474753

def getWeightedHist(selection, histName, leftCut=None, rightCut=None):
	suffix = None
	if selection=='': suffix = '_'
	elif selection=='rs': suffix = '_rs_'
	elif selection=='fjrb': suffix = '_fjrb_'
	if suffix==None: 
		print 'wrong parameter'
		return None
	WHIZARD_filename =  'newhist/hist_1pho'+suffix+'barrel_top_WHIZARD.root'
	TTJets1l_filename = 'newhist/hist_1pho'+suffix+'barrel_top_TTJets1l.root'
	TTJets2l_filename = 'newhist/hist_1pho'+suffix+'barrel_top_TTJets2l.root'
	
	if leftCut==None and rightCut==None:
		# 1d histogram 
		whiz_hist = get1DHist(WHIZARD_filename, histName)
		whiz_hist.Scale(WHIZARD_SF)
		TTjets1l_hist = get1DHist(TTJets1l_filename, histName)
		TTjets1l_hist.Scale(TTJets1l_SF)
		TTjets2l_hist = get1DHist(TTJets2l_filename, histName)
		TTjets2l_hist.Scale(TTJets2l_SF)

		sumHist = whiz_hist
		sumHist.Add(TTjets1l_hist)
		sumHist.Add(TTjets2l_hist)
		return sumHist
	else:
		# 2d histogram projection
		whiz_histp = getTemplFrom2Dhist(WHIZARD_filename, histName, 'X', leftCut, rightCut)
		whiz_histp.Scale(WHIZARD_SF)
		TTJets1l_histp = getTemplFrom2Dhist(TTJets1l_filename, histName, 'X', leftCut, rightCut)
		TTJets1l_histp.Scale(TTJets1l_SF)
		TTJets2l_histp = getTemplFrom2Dhist(TTJets2l_filename, histName, 'X', leftCut, rightCut)
		TTJets2l_histp.Scale(TTJets2l_SF)

		sumHistp = whiz_histp
		sumHistp.Add(TTJets1l_histp)
		sumHistp.Add(TTJets2l_histp)
		return sumHistp


# fit range
lowFitrange = -0.5
highFitrange = 19.99
#highFitrange = 9.99

# number of pseudo experiments. off by default
NpseudoExp = 0

# photon Et range
#phoEtrange = '_60_up_'
phoEtrange = '_'
FitVarname = 'SCRChHadIso'
hist2d_name = 'photon1'+phoEtrange+'Sigma_ChSCRIso'
if phoEtrange=='_': phoEtrange = ''
hist_sig_templ_name = 'photon1'+phoEtrange+'ChHadRandIso'
hist_sig_name = 'photon1ChHadSCRIso'

usePhoIso = False
#usePhoIso = True
# for Pho Iso
if usePhoIso:
	FitVarname = 'SCRPhoIso'
	hist2d_name = 'photon1_Sigma_PhoSCRIso'
	hist_sig_templ_name = 'photon1PhoRandIso'
	hist_sig_name = 'photon1PhoSCRIso'
	lowFitrange = -3.0
	highFitrange = 9.99


# sideband in sigmaIetaIeta
sbLeft = 0.012
sbRight = 0.016

# signal template, data driven, random cone isolation
sig_templ = getWeightedHist('', hist_sig_templ_name)
# background template, data driven, sideband in sihih
bckg_templ = getWeightedHist('', hist2d_name, sbLeft, sbRight)
# data, no cut on isolation, cut on sihih
pseudodata = getWeightedHist('', hist2d_name, 0.0, 0.012)


### MC truth for SCR isolation selection, for true signal fraction calculation
pseudosignal = getWeightedHist('rs', hist2d_name, 0.0, 0.012)

leftbin = pseudodata.FindBin(lowFitrange)
rightbin = pseudodata.FindBin(highFitrange)
print 'calculating integrals in bins',leftbin,rightbin
MCtrueSelSignalFraction = pseudosignal.Integral(leftbin, rightbin)/pseudodata.Integral(leftbin, rightbin)
########################################

if phoEtrange:
	print 'Photon Et range: ',phoEtrange
else:
	print 'Photon Et range: all'

print 'Signal Integral: ',sig_templ.Integral()
print 'Background Integral: ',bckg_templ.Integral()
print 'Data Integral: ',pseudodata.Integral()

# make equal binning in all histograms
bckg_templ.Rebin(2)
pseudodata.Rebin(2)

#boundaries = optimizeBinBoundaries(pseudodata,15,lowFitrange,highFitrange)
boundaries = optimizeBinBoundaries(bckg_templ, 10, lowFitrange, highFitrange)

print 'new number of boundaries ',len(boundaries)
pseudodataR = pseudodata.Rebin(len(boundaries)-1,pseudodata.GetName()+'rebin',boundaries)
bckg_templR = bckg_templ.Rebin(len(boundaries)-1,bckg_templ.GetName()+'rebin',boundaries)
sig_templR  = sig_templ.Rebin(len(boundaries)-1,sig_templ.GetName()+'rebin',boundaries)
# to avoid signal going below fit range in the plot
print 'sig_templR underflow ', sig_templR.GetBinContent(0)
print 'settin to 0'
sig_templR.SetBinContent(0,0)

pseudodataU = makeUniformHist(pseudodataR)
bckg_templU = makeUniformHist(bckg_templR)
sig_templU  = makeUniformHist(sig_templR)

print '#'*50

# does not work well
#makeFit(FitVarname,lowFitrange,highFitrange, sig_templR, bckg_templR, pseudodataR, 'fit_'+FitVarname+'.png')

# do fit for fixed bin size histogram
lowUFitRange = pseudodataR.FindBin(lowFitrange)
if lowUFitRange == 0:
	print 'lower bin is underfow, setting to 1'
	lowUFitRange = 1
highUFitRange = pseudodataR.FindBin(highFitrange) + 1
if pseudodataR.GetBinLowEdge(highUFitRange) == highFitrange:
	print 'upper bin in on the border of fit range, reducing'
	highUFitRange -= 1

print 'fitting in the bin range ',lowUFitRange, highUFitRange
(fitSigFrac,fitSigFracErr) = makeFit(FitVarname+' bin number',lowUFitRange,highUFitRange, sig_templU, bckg_templU, pseudodataU, 'fit_'+phoEtrange+FitVarname+'_uniform.png')

# pseudo-experiments:
pe_results = ROOT.TH1F('pe_results','Pseudo-experiments',100,0,1)
#pe_results = ROOT.TH1F('pe_results','Pseudo-experiments',50,0.4,0.9)
for npi in xrange(NpseudoExp):
	print '-'*80
	print 'pe # ',npi
	pseudodataUrand = makeRandUniformHist(pseudodataR)
	bckg_templUrand = makeRandUniformHist(bckg_templR)
	sig_templUrand  = makeRandUniformHist(sig_templR)
	(fitSigFracRand,fitSigFracErrRand) = makeFit(FitVarname+' bin number',lowUFitRange, highUFitRange, sig_templUrand, bckg_templUrand, pseudodataUrand, '')
	pe_results.Fill(fitSigFracRand)

# draw the result
ROOT.gStyle.SetOptFit(111)
c1 = ROOT.TCanvas('c1','c1',800,800)
pe_results.GetXaxis().SetTitle('signal fraction')
pe_results.Draw()
pe_results.Fit('gaus')

if NpseudoExp > 0:
	c1.SaveAs('pe_results.png')

ROOT.gStyle.SetOptStat(0)

leg = ROOT.TLegend(0.6,0.7,0.99,0.94)
leg.SetFillColor(0)


leftbin = lowUFitRange
rightbin = highUFitRange

print 'bins for normalization ',leftbin,rightbin
bckg_templR.Scale((1.0-fitSigFrac)*pseudodataR.Integral(leftbin,rightbin)/bckg_templR.Integral(leftbin,rightbin))
sig_templR.Scale(fitSigFrac*pseudodataR.Integral(leftbin,rightbin)/sig_templR.Integral(leftbin,rightbin))

sig_templR.SetFillColor(ROOT.kGreen)
bckg_templR.SetFillColor(ROOT.kBlue)
sig_templR.SetLineColor(ROOT.kGreen)
bckg_templR.SetLineColor(ROOT.kBlue)
sig_templR.GetXaxis().SetRangeUser(lowFitrange,highFitrange)
bckg_templR.GetXaxis().SetRangeUser(lowFitrange,highFitrange)

stack = ROOT.THStack(FitVarname+'_stack',FitVarname+'_stack')
stack.Add(bckg_templR)
stack.Add(sig_templR)
stack.SetTitle('')
stack.Draw('hist')
stack.GetXaxis().SetTitle('photon '+FitVarname+' (GeV)')
stack.SetMaximum(1.2*stack.GetMaximum())
stack.GetXaxis().SetRangeUser(lowFitrange,highFitrange)

leg.AddEntry(pseudodataR, 'Data', 'lfp')
leg.AddEntry(sig_templR, 'Signal', 'lf')
leg.AddEntry(bckg_templR, 'Background', 'lf')
leg.Draw()
pseudodataR.SetMarkerStyle(8)
pseudodataR.Draw('esame')
c1.SaveAs('plot_'+phoEtrange+FitVarname+'.png')

leg.Clear()
## compare template shapes here

# signal template: compare MC truth photon isolation with random cone isolation
trueSignalIso = getWeightedHist('rs', hist2d_name, 0.0, 0.012)
trueSignalIso.Rebin(2)
trueSignalIso.GetXaxis().SetRangeUser(lowFitrange,highFitrange)
trueSignalIsoInt = trueSignalIso.Integral()
trueSignalIso.Scale(1.0/trueSignalIso.Integral())

randCone_Iso = getWeightedHist('', hist_sig_templ_name)
randCone_Iso.Scale(1.0/randCone_Iso.Integral())
randCone_Iso.SetLineColor(2)

leg.AddEntry(trueSignalIso, 'MC truth signal','lf')
leg.AddEntry(randCone_Iso, 'Random cone', 'lf')
randCone_Iso.SetTitle('')
randCone_Iso.GetXaxis().SetTitle('photon '+FitVarname+' (GeV)')
randCone_Iso.Draw()
trueSignalIso.Draw('same')
leg.Draw()
c1.SaveAs(hist_sig_name+'_sig'+phoEtrange+'templ.png')
c1.SetLogy(1)
c1.SaveAs(hist_sig_name+'_sig'+phoEtrange+'templ_log.png')
c1.SetLogy(0)
leg.Clear()

# background template comparison:
# pseudodata sideband
sbIso = getWeightedHist('', hist2d_name, sbLeft, sbRight)
sbIso.Rebin(2)
# MC truth signal iso in side-band
sbIsoMCtrue = getWeightedHist('rs', hist2d_name, sbLeft, sbRight)
sbIsoMCtrue.Rebin(2)
# calculate true signal fraction in side-band
print 'Side-band used: ',sbLeft,sbRight
print 'Fit range considered: ',lowFitrange, highFitrange
leftbin = sbIso.FindBin(lowFitrange)
rightbin = sbIso.FindBin(highFitrange)
print 'Fraction of MC truth signal in the side-band with respect to pseudodata: ',\
	sbIsoMCtrue.Integral(leftbin,rightbin)/sbIso.Integral(leftbin,rightbin)

sbIso.Scale(1.0/sbIso.Integral())
sbIso.SetLineColor(2)

# MC truth background (fjrb) shape
true_bckg = getWeightedHist('fjrb',hist2d_name, 0.0, 0.012)
true_bckg.Rebin(2)
true_bckg.Scale(1.0/true_bckg.Integral())


#true_sel_bckg = getWeightedHist('fjrb', hist_sig_name)
#true_sel_bckg.Scale(1.0/true_sel_bckg.Integral())
#true_sel_bckg.SetLineColor(3)

leg.AddEntry(true_bckg,'MC truth backg','lf')
leg.AddEntry(sbIso,'sihih side band', 'lf')
#leg.AddEntry(true_sel_bckg,'MC truth backg, nom. selection', 'lf')

#true_sel_bckg.Draw()
true_bckg.SetTitle('')
true_bckg.GetXaxis().SetTitle('photon '+FitVarname+' (GeV)')
true_bckg.Draw()
sbIso.Draw('same')
leg.Draw()
c1.SaveAs(hist_sig_name+phoEtrange+'_bckg.png')
c1.SetLogy(1)
c1.SaveAs(hist_sig_name+phoEtrange+'_bckg_log.png')
c1.SetLogy(0)

true_bckg.GetXaxis().SetRangeUser(lowFitrange,highFitrange)
true_bckg.Scale(1.0/true_bckg.Integral())
sbIso.GetXaxis().SetRangeUser(lowFitrange,highFitrange)
sbIso.Scale(1.0/sbIso.Integral())
#true_sel_bckg.GetXaxis().SetRangeUser(lowFitrange,highFitrange)
#true_sel_bckg.Scale(1.0/true_sel_bckg.Integral())

#true_sel_bckg.Draw()
true_bckg.Draw()
sbIso.Draw('same')
leg.Draw()
c1.SaveAs(hist_sig_name+phoEtrange+'_bckg_fitrange.png')
leg.Clear()


# find MC truth signal fraction in nominal selection #############################################################################
# extract "rs" parts, divide by "photon1ChHadSCRIso" sum
total_sig = getWeightedHist('', hist_sig_name)

print 'MC truth signal fraction afetr full selection is ', trueSignalIsoInt/total_sig.Integral()
print 'MC truth signal fraction if SCR iso in fit range used as selection: ',MCtrueSelSignalFraction
