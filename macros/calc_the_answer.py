import ROOT

def integral(file,histname):
	err = ROOT.Double(0.0)
	hist = file.Get(histname)
	integr = hist.IntegralAndError(0,hist.GetNbinsX() + 1, err)
	print integr,err
	return integr,err

def integral_bins(file,histname,bin1,bin2):
	err = ROOT.Double(0.0)
	hist = file.Get(histname)
	integr = hist.IntegralAndError(bin1, bin2, err)
	print integr,err
	return integr,err

def integral_one_bin(file,histname):
	err = ROOT.Double(0.0)
	hist = file.Get(histname)
	integr = hist.IntegralAndError(1, 1, err)
	print integr,err
	return integr,err

# for ttbar acceptance
TTJets1l_num        = 24849110
TTJets2l_num        = 12086717
TTJetsHad_num       = 31178278

## >=3j
photnPurity = 0.562
photnPurityErr = 0.065

#QCDSF = 1.338 #0.1429
#QCDSFErr = 0.032 #0.0033

M3TopSF = 0.964
M3TopSFErr = 0.0082

M3WJetsSF = 2.36
M3WJetsSFErr = 0.05

eleFakeSF = 1.5 
eleFakeSFErr = 0.2


preselFileName = 'templates_presel_scaled.root'
barrelFileName = 'templates_barrel_scaled.root'

def doTheCalculation():
	print '#$'*50
	print 'input data'
	print 'photon Purity: ',photnPurity,' +-',photnPurityErr
	#print 'QCDSF: ',QCDSF,' +-',QCDSFErr
	print 'TopSF: ',M3TopSF,' +-',M3TopSFErr
	print 'WJetsSF: ',M3WJetsSF,' +-',M3WJetsSFErr
	print 'eleFakeSF: ',eleFakeSF,' +-',eleFakeSFErr
	
	print
	print 'Opening file with pre-selection histograms'
	preselfile = ROOT.TFile(preselFileName,'READ')
	print 'File was last modified: ',preselfile.GetModificationDate().AsString()
	print
	topPreselInt,topPreselErr = integral(preselfile,'TTJets_MET')
	print 'TTJets presel events expected ',topPreselInt,' +-',topPreselErr,'(MC stat)'
	print
	whizPreselInt,whizPreselErr = integral(preselfile,'TTGamma_MET')
	print 'WHIZARD presel events expected ',whizPreselInt,' +-',whizPreselErr,'(MC stat)'
	print
	print
	# total number of events with top quark is the sum of them
	topPreselInt += whizPreselInt
	topPreselErr = (topPreselErr**2 + whizPreselErr**2 + (topPreselInt*M3TopSFErr/M3TopSF)**2 )**0.5
	print 'taking into account fit error and ttgamma contribution: ',topPreselInt,' +-',topPreselErr,'(MC stat+fit)'
	print
	# number of signal events in pre-selection
	ttgammaSigPreselInt, ttgammaSigPreselErr = integral_bins(preselfile,'TTGamma_MCcategory',3,4)
	print 'number of signal events with one or two leptons in preselection ', ttgammaSigPreselInt, ' +-', ttgammaSigPreselErr,'(MC stat)'
	
	print 'Opening file with photon in barrel'
	barrelfile = ROOT.TFile(barrelFileName,'READ')
	print 'File was last modified: ',barrelfile.GetModificationDate().AsString()
	print
	ttgammaBarrInt,ttgammaBarrErr = integral_bins(barrelfile,'TTGamma_signal_MCcategory',3,4)
	print 'ttgamma signal events with one or two leptons and matched photon in barrel expected ',ttgammaBarrInt, ' +-',ttgammaBarrErr,'(MC stat)'
	print
	
	phoAcc = ttgammaBarrInt/ttgammaSigPreselInt
	# These errors are not independent
	phoAccErr = phoAcc * ((ttgammaSigPreselErr/ttgammaSigPreselInt) + (ttgammaBarrErr/ttgammaBarrInt))
	print 'signal photon acceptance',phoAcc,' +-',phoAccErr
	print
	
	print 'reconstruction efficiency: matched signal events in barrel'
	ttgammaSigBarrInt,ttgammaSigBarrErr = integral(barrelfile,'TTGamma_signal_MET')
	
	print 'fiducial photons in barrel after pre-selection'
	genPhoInAcc,genPhoInAccErr = integral_one_bin(preselfile,'TTGamma_genPhoRegionWeight')
	phoRecoEff = ttgammaSigBarrInt/genPhoInAcc
	phoRecoEffErr = phoRecoEff * ((genPhoInAccErr/genPhoInAcc)**2 + (ttgammaSigBarrErr/ttgammaSigBarrInt)**2)**0.5
	print 'signal photon reconstruction efficiency (with truth matching to photon)', phoRecoEff,' +-',phoRecoEffErr
	print

	print 'TTJets acceptance calculation'
	accfile = ROOT.TFile('ttbar_acceptance.root','READ')
	print 'File was last modified: ',accfile.GetModificationDate().AsString()
	ttjets1l_top, ttjets1l_topErr = integral_one_bin(accfile,'TTJets1l_presel_MCcategory')
	ttjets2l_top, ttjets2l_topErr = integral_one_bin(accfile,'TTJets2l_presel_MCcategory')
	ttjetsHad_top, ttjetsHad_topErr = integral_one_bin(accfile,'TTJetsHad_presel_MCcategory')
	Whad = 0.676
	TTJets_topEffAcc = Whad*(1.0-Whad)*2 * ttjets1l_top / TTJets1l_num + (1.0-Whad)*(1.0-Whad) * ttjets2l_top / TTJets2l_num + Whad*Whad * ttjetsHad_top / TTJetsHad_num
	TTJets_topEffAccErr = ( (Whad*(1.0-Whad)*2 * ttjets1l_topErr / TTJets1l_num)**2 + ((1.0-Whad)*(1.0-Whad) * ttjets2l_topErr / TTJets2l_num)**2 + (Whad*Whad * ttjetsHad_topErr / TTJetsHad_num)**2 )**0.5
	print 'Inclusive ttbar acceptance ',TTJets_topEffAcc,' +-',TTJets_topEffAccErr
	print 
	print 'TTgamma acceptance calculation'
	sigaccfile = ROOT.TFile('signalAcc.root','READ')
	print 'File was last modified: ',sigaccfile.GetModificationDate().AsString()
	ttgamma_1l_2l_num, ttgamma_1l_2l_numErr = integral_bins(sigaccfile,'allCategory',3,4)
	ttgamma_1l_2l_Visnum, ttgamma_1l_2l_VisnumErr = integral_bins(sigaccfile,'VisAllCategory',3,4)
	
	ttgamma_1l_2l_presel, ttgamma_1l_2l_preselErr = integral_bins(accfile, 'TTGamma_presel_MCcategory',3,4)
	ttgamma_1l_2l_sig, ttgamma_1l_2l_sigErr = integral_bins(accfile, 'TTGamma_signal_MCcategory',3,4)
	
	TTGamma_topEffAcc = ttgamma_1l_2l_presel / ttgamma_1l_2l_num
	TTGamma_topEffAccErr = TTGamma_topEffAcc * ((ttgamma_1l_2l_preselErr/ttgamma_1l_2l_presel)**2 + (ttgamma_1l_2l_numErr/ttgamma_1l_2l_num)**2)**0.5
	print 'Acceptance for ttgamma after top selection ',TTGamma_topEffAcc, ' +/-',TTGamma_topEffAccErr
	
	TTGammaVis_topAcc = ttgamma_1l_2l_presel / ttgamma_1l_2l_Visnum
	TTGammaVis_topAccErr = TTGammaVis_topAcc * ((ttgamma_1l_2l_preselErr/ttgamma_1l_2l_presel)**2 + (ttgamma_1l_2l_VisnumErr/ttgamma_1l_2l_Visnum)**2)**0.5
	print 'Acceptance for visible ttgamma after top selection ',TTGammaVis_topAcc, ' +/-',TTGammaVis_topAccErr
	print
	print 'Product of ttgamma acceptance for top and photon calculated at once',ttgamma_1l_2l_sig/ttgamma_1l_2l_num
	print 
	
	DataBarrInt,DataBarrErr = integral(barrelfile,'Data_MET')
	print 'Data events',DataBarrInt,' +-',DataBarrErr

	print 'Getting expected background events with genuine photons from MC'
	print 'Vgamma:'
	VgammaBarrInt,VgammaBarrErr = integral(barrelfile,'Vgamma_signal_MET')
	
	print 'SingleTop:'
	SingleTopBarrInt,SingleTopBarrErr = integral(barrelfile,'SingleTop_signal_MET')
	print 'WJets:'
	WJetsBarrInt,WJetsBarrErr = integral(barrelfile,'WJets_signal_MET')
	print 'ZJets:'
	ZJetsBarrInt,ZJetsBarrErr = integral(barrelfile,'ZJets_signal_MET')
	print 'Diboson:'
	OtherBarrInt,OtherBarrErr = integral(barrelfile,'Diboson_signal_MET')
	
	bckgPhoExp = VgammaBarrInt + SingleTopBarrInt + WJetsBarrInt + ZJetsBarrInt + OtherBarrInt
	SqBckgPhoErr = VgammaBarrErr**2 + SingleTopBarrErr**2 + WJetsBarrErr**2 + (WJetsBarrInt*M3WJetsSFErr/M3WJetsSF)**2 + ZJetsBarrErr**2 + OtherBarrErr**2
	print 'total background with genuine photons expected', bckgPhoExp,' +-',SqBckgPhoErr**0.5
	
	
	print 'Now getting expected background events with photon faked by electron from MC'
	print 'TTJets:'
	TTJetsEleInt,TTJetsEleErr = integral(barrelfile,'TTJets_electron_MET')
	print 'Vgamma:'
	VgammaEleInt,VgammaEleErr = integral(barrelfile,'Vgamma_electron_MET')
	
	print 'SingleTop:'
	SingleTopEleInt,SingleTopEleErr = integral(barrelfile,'SingleTop_electron_MET')
	print 'WJets:'
	WJetsEleInt,WJetsEleErr = integral(barrelfile,'WJets_electron_MET')
	print 'ZJets:'
	ZJetsEleInt,ZJetsEleErr = integral(barrelfile,'ZJets_electron_MET')
	print 'Diboson:'
	OtherEleInt,OtherEleErr = integral(barrelfile,'Diboson_electron_MET')
	bckgEleExp = TTJetsEleInt + VgammaEleInt + SingleTopEleInt + WJetsEleInt + ZJetsEleInt + OtherEleInt
	SqBckgEleErr = TTJetsEleErr**2 + VgammaEleErr**2 + SingleTopEleErr**2 + WJetsEleErr**2 + (WJetsEleInt*M3WJetsSFErr/M3WJetsSF)**2 + ZJetsEleErr**2 + OtherEleErr**2
	print 'background with photon fakes expected',bckgEleExp,' +-',SqBckgEleErr**0.5
	
	bckgExp = bckgPhoExp + bckgEleExp*eleFakeSF
	bckgExpErr = ( SqBckgPhoErr + SqBckgEleErr*eleFakeSF*eleFakeSF )**0.5
	print 'total background expected',bckgExp,' +-',bckgExpErr
	print

	xsRatio = (DataBarrInt * photnPurity - bckgExp) / phoAcc / TTGamma_topEffAcc / topPreselInt * TTJets_topEffAcc
	xsRatioRelErr = (  ( bckgExpErr**2 + (DataBarrErr*photnPurity)**2 ) / (DataBarrInt * photnPurity - bckgExp)**2 + 
						(phoAccErr/phoAcc)**2 + 
						(topPreselErr/topPreselInt)**2 
					)**0.5
	print '*'*80
	print 'final answer: cross section ratio:'
	print xsRatio,' +-',xsRatio*xsRatioRelErr
	print '*'*80
	
	
	vis_xsRatio = (DataBarrInt * photnPurity - bckgExp) / phoRecoEff / TTGammaVis_topAcc / topPreselInt * TTJets_topEffAcc
	vis_xsRatioErr = ( ( bckgExpErr**2 + (DataBarrErr*photnPurity)**2 ) / (DataBarrInt * photnPurity - bckgExp)**2 + 
						(topPreselErr/topPreselInt)**2 +
						(phoRecoEffErr/phoRecoEff)**2
					)**0.5
	print 'visible cross section ratio:'
	print vis_xsRatio,' +-',vis_xsRatio*vis_xsRatioErr
	print '*'*80
	return (xsRatio, xsRatio*xsRatioRelErr)
	
	
