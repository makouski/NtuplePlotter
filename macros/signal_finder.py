import ROOT
from math import exp

barrelFileName = 'templates_barrel_scaled.root'
myfile = None

def integral(histname):
	global myfile
	if myfile == None:
		myfile = ROOT.TFile(barrelFileName,'READ')
	err = ROOT.Double(0.0)
	hist = myfile.Get(histname)
	integr = hist.IntegralAndError(0,hist.GetNbinsX() + 1, err)
	print integr,err
	return integr,err

eleFakeSF = 1.5 
eleFakeSFErr = 0.2

# parameters used in chi square calculation
photnPurity = 0.562
photnPurityErr = 0.065

M3_photon_topFrac = 0.619
M3_photon_topFracErr = 0.089

Ndata = 915.0
NdataErr = 30.2


def readSamples(suffix):
	var = 'MET'
	samplnames = ['TTGamma', 'TTJets', 'Vgamma', 'SingleTop', 'WJets', 'ZJets', 'Diboson', 'QCD']
	samples = {}
	errors = {}
	for n in samplnames:
		if n=='QCD' and (suffix=='signal' or suffix=='electron'):
			samples[n] = 0.0
			errors[n] = 0.0
			continue
		print 'getting counts for ',n,suffix
		if n=='QCD':
			int,err = integral(n+'_'+var)
		else:
			int,err = integral(n+'_'+suffix+'_'+var)
		samples[n] = int
		errors[n] = err
		
	# deal with errors later
	return samples

def Negamma(pho,ele,fake, ttgSF, VgSF, jgSF):
	sum = 0.0
	for sampln in pho:
		npho = pho[sampln]
		if sampln == 'TTGamma':
			npho *= ttgSF
		if sampln == 'Vgamma':
			npho *= VgSF
		sum += npho
	
	for sampln in ele:
		nele = ele[sampln]
		nele *= eleFakeSF
		if sampln == 'TTGamma':
			nele *= ttgSF
		if sampln == 'Vgamma':
			nele *= VgSF
		sum += nele
	
	return sum

def NallMC(pho,ele,fake, ttgSF, VgSF, jgSF):
	# all MC is egamma plus fakes
	sum = Negamma(pho,ele,fake, ttgSF, VgSF, jgSF)
	
	for sampln in fake:
		nfake = fake[sampln]
		nfake *= jgSF
		if sampln == 'TTGamma':
			nfake *= ttgSF
		if sampln == 'Vgamma':
			nfake *= VgSF
		sum += nfake
	
	return sum


def Ntop(pho,ele,fake, ttgSF, VgSF, jgSF):
	toppho = {}
	topele = {}
	topfake = {}
	for t in ['TTGamma', 'TTJets']:
		toppho[t] = pho[t]
		topele[t] = ele[t]
		topfake[t] = fake[t]
	
	return NallMC(toppho,topele,topfake, ttgSF, VgSF, jgSF)



def getChiSq(pho,ele,fake, ttgSF, VgSF, jgSF):
	allmc = NallMC(pho,ele,fake, ttgSF, VgSF, jgSF)
	top = Ntop(pho,ele,fake, ttgSF, VgSF, jgSF)
	egamma = Negamma(pho,ele,fake, ttgSF, VgSF, jgSF)
	cs = 0.0
	cs += (allmc - Ndata)**2 / NdataErr**2
	cs += (egamma/allmc - photnPurity)**2 / photnPurityErr**2
	cs += (top/allmc - M3_photon_topFrac)**2 / M3_photon_topFracErr**2
	return cs

def getLikelihood(pho,ele,fake, ttgSF, VgSF, jgSF):
	return exp(-0.5*getChiSq(pho,ele,fake, ttgSF, VgSF, jgSF))

def seq(start, stop, step=1):
    n = int(round((stop - start)/float(step)))
    if n > 1:
        return([start + step*i for i in range(n+1)])
    else:
        return([])

def calculateTTGamma():
	# read all the numbers from the file
	pho = readSamples('signal')
	ele = readSamples('electron')
	fake = readSamples('fake')
	
	ttghist = ROOT.TH1F('ttghist','Marginalized Likelihood',140, 0.2+0.005,1.6+0.005)
	vghist = ROOT.TH1F('vghist','Marginalized Likelihood',130, -2.0+0.025,4.5+0.025)
	jghist = ROOT.TH1F('jghist','Marginalized Likelihood',120, 0.6+0.005,1.8+0.005)
	maxlk = -1.0
	bestttgSF = -1.0
	bestVgSF = -1.0
	bestjgSF = -1.0
	
	for ttgSF in seq(0.2, 1.6, 0.01):
		for VgSF in seq(-2.0, 4.5, 0.05):
			for jgSF in seq(0.6, 1.8, 0.01):
				lk = getLikelihood(pho,ele,fake,ttgSF, VgSF, jgSF)
				ttghist.Fill(ttgSF,lk*0.001)
				vghist.Fill(VgSF,lk*0.0001)
				jghist.Fill(jgSF,lk*0.001)
				if lk > maxlk:
					maxlk = lk
					bestttgSF = ttgSF
					bestVgSF = VgSF
					bestjgSF = jgSF
	
	print 'max likelihood ',maxlk
	print 'best ttgSF',bestttgSF
	print 'best VgSF',bestVgSF
	print 'best jgSF',bestjgSF
	b_allmc = NallMC(pho,ele,fake, bestttgSF, bestVgSF, bestjgSF)
	b_top = Ntop(pho,ele,fake, bestttgSF, bestVgSF, bestjgSF)
	b_egamma = Negamma(pho,ele,fake, bestttgSF, bestVgSF, bestjgSF)
	print 'total number of MC events ',b_allmc, '  Data:',Ndata
	print 'photon purity in MC ',b_egamma/b_allmc, '  template fit:',photnPurity
	print 'top fraction in MC ',b_top/b_allmc, '  M3 fit:',M3_photon_topFrac
	
	ROOT.gStyle.SetOptFit(111)
	ccc = ROOT.TCanvas('ccc','ccc',800,800)
	ttghist.Draw()
	ttghist.GetXaxis().SetTitle('TTGamma Scale Factor')
	ttghist.Fit('gaus')
	ccc.SaveAs('TTGamma_SF_Lkhood.png')
	fit = ttghist.GetFunction('gaus')
	bestttgSFErr = fit.GetParameter(2)
	
	vghist.Draw()
	vghist.GetXaxis().SetTitle('Vgamma Scale Factor')
	vghist.Fit('gaus')
	ccc.SaveAs('Vgamma_SF_Lkhood.png')
	
	jghist.Draw()
	jghist.GetXaxis().SetTitle('Jet to Photon Scale Factor')
	jghist.Fit('gaus')
	ccc.SaveAs('jet_gamma_SF_Lkhood.png')
	
	ttgammaSig = bestttgSF*pho['TTGamma']
	ttgammaSigErr = bestttgSFErr*pho['TTGamma']
	print 'number of signal events', ttgammaSig, ' +/-',ttgammaSigErr
	return ttgammaSig,ttgammaSigErr

#calculateTTGamma()
	