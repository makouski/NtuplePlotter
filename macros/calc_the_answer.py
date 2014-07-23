import ROOT

def integral(file,histname):
	err = ROOT.Double(0.0)
	hist = file.Get(histname)
	integr = hist.IntegralAndError(0,hist.GetNbinsX() + 1, err)
	print integr,err
	return integr,err

## >=4j
#photnPurity = 0.5706
#photnPurityErr = 0.073

#QCDSF = 0.1048
#QCDSFErr = 0.0086

#M3TopSF = 0.933
#M3TopSFErr = 0.0066

#M3WJetsSF = 1.61
#M3WJetsSFErr = 0.064

## >=3j
photnPurity = 0.561
photnPurityErr = 0.053

QCDSF = 1.338 #0.1429
QCDSFErr = 0.032 #0.0033

M3TopSF = 0.971 #0.931
M3TopSFErr = 0.0095 #0.0089

M3WJetsSF = 1.66 #1.595
M3WJetsSFErr = 0.04

preselFileName = 'templates_presel_scaled.root'
barrelFileName = 'templates_barrel_scaled.root'

def doTheCalculation():
	print '#$'*50
	print 'input data'
	print 'photon Purity: ',photnPurity,' +-',photnPurityErr
	print 'QCDSF: ',QCDSF,' +-',QCDSFErr
	print 'TopSF: ',M3TopSF,' +-',M3TopSFErr
	print 'WJetsSF: ',M3WJetsSF,' +-',M3WJetsSFErr
	
	print
	print 'Opening file with pre-selection histograms'
	preselfile = ROOT.TFile(preselFileName,'READ')
	print 'File was last modified: ',preselfile.GetModificationDate().AsString()
	print
	topPreselInt,topPreselErr = integral(preselfile,'TTJets_MET')
	print 'TTJets presel events expected ',topPreselInt,' +-',topPreselErr,'(MC stat)'
	topPreselErr = (topPreselErr**2 + (topPreselInt*M3TopSFErr/M3TopSF)**2 )**0.5
	print 'taking into account fit error: ',topPreselInt,' +-',topPreselErr,'(MC stat+fit)'
	print
	whizPreselInt,whizPreselErr = integral(preselfile,'TTGamma_MET')
	print 'WHIZARD presel events expected ',whizPreselInt,' +-',whizPreselErr,'(MC stat)'
	print
	print 'Opening file with photon in barrel'
	barrelfile = ROOT.TFile(barrelFileName,'READ')
	print 'File was last modified: ',barrelfile.GetModificationDate().AsString()
	print
	whizBarrInt,whizBarrErr = integral(barrelfile,'TTGamma_MET')
	print 'WHIZARD barrel events expected ',whizBarrInt, ' +-',whizBarrErr,'(MC stat)'
	print

	phoAcc = whizBarrInt/whizPreselInt
	# These errors are not completely independent. Treat them if they are for now.
	phoAccErr = phoAcc * ((whizBarrErr/whizBarrInt)**2 + (whizPreselErr/whizPreselInt)**2)**0.5
	print 'signal photon acceptance',phoAcc,' +-',phoAccErr

	print

	DataBarrInt,DataBarrErr = integral(barrelfile,'Data_MET')
	print 'Data events',DataBarrInt,' +-',DataBarrErr

	print 'Now getting expected background events from MC'
	print 'Vgamma:'
	VgammaBarrInt,VgammaBarrErr = integral(barrelfile,'Vgamma_MET')
	print 'SingleTop:'
	SingleTopBarrInt,SingleTopBarrErr = integral(barrelfile,'SingleTop_MET')
	print 'WJets:'
	WJetsBarrInt,WJetsBarrErr = integral(barrelfile,'WJets_MET')
	print 'ZJets:'
	ZJetsBarrInt,ZJetsBarrErr = integral(barrelfile,'ZJets_MET')
	print 'Diboson:'
	OtherBarrInt,OtherBarrErr = integral(barrelfile,'Diboson_MET')
	print 'QCD:'
	QCDBarrInt,QCDBarrErr = integral(barrelfile,'QCD_MET')
	bckgExp = VgammaBarrInt+SingleTopBarrInt+WJetsBarrInt+ZJetsBarrInt+OtherBarrInt+QCDBarrInt
	SqBckgErr = 0.0
	for bc,bcerr in [(VgammaBarrInt,VgammaBarrErr),
						(SingleTopBarrInt,SingleTopBarrErr),
						(WJetsBarrInt, (WJetsBarrErr**2 + (WJetsBarrInt*M3WJetsSFErr/M3WJetsSF)**2)**0.5 ),
						(ZJetsBarrInt,ZJetsBarrErr),
						(OtherBarrInt,OtherBarrErr),
						(QCDBarrInt, (QCDBarrErr**2 + (QCDBarrInt*QCDSFErr/QCDSF)**2)**0.5 )
					]:
		SqBckgErr += (bcerr)**2
	bckgExpErr = (SqBckgErr)**0.5
	print 'total background expected',bckgExp,' +-',bckgExpErr
	print
	xsRatio = (DataBarrInt - bckgExp) * photnPurity / phoAcc / topPreselInt
	xsRatioRelErr = ( (DataBarrErr**2 + bckgExpErr**2)/(DataBarrInt - bckgExp)**2 + 
						(photnPurityErr/photnPurity)**2 + 
						(phoAccErr/phoAcc)**2 + 
						(topPreselErr/topPreselInt)**2 
					)**0.5
	print '*'*80
	print 'final answer: cross section ratio:'
	print xsRatio,' +-',xsRatio*xsRatioRelErr
	print '*'*80
	
	vis_xsRatio = (DataBarrInt - bckgExp) * photnPurity / topPreselInt
	vis_xsRatioErr = ( (DataBarrErr**2 + bckgExpErr**2)/(DataBarrInt - bckgExp)**2 + 
						(photnPurityErr/photnPurity)**2 + 
						(topPreselErr/topPreselInt)**2 
					)**0.5
	print 'visible cross section ratio:'
	print vis_xsRatio,' +-',vis_xsRatio*vis_xsRatioErr
	print '*'*80
	return (xsRatio, xsRatio*xsRatioRelErr)
	
	
